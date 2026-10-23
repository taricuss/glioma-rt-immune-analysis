
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')
from statsmodels.regression.mixed_linear_model import MixedLM
import json, csv, sys
from math import sqrt
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
from lmm_interaction_pipeline import merge_inputs  # noqa: E402

PROJ = SCRIPT_DIR.parent
INTER = PROJ / "outputs" / "intermediate"

OUT_CSV = str(INTER / "post_data_power_table.csv")
OUT_JSON = str(INTER / "post_data_power_audit.json")

# Effective discovery design is DERIVED from the actual merged inputs (never
# hardcoded). After the QC fix the modelled set is the full clinical_canonical
# cohort, but the derivation is retained so the table can never drift from the
# data if the modelled set changes again.
_fr = pd.read_csv(INTER / "sn_proportions_per_sample.csv")
if "sample_id" not in _fr.columns:
    _fr = _fr.rename(columns={_fr.columns[0]: "sample_id"})
_clin = pd.read_csv(INTER / "clinical_canonical.csv")
_pur = INTER / "estimate_scores.csv"
_pur = pd.read_csv(_pur) if Path(_pur).exists() else None
_merged = merge_inputs(_fr, _clin, _pur)

N_PATIENTS = int(_merged["patient_id"].nunique())
_split = _merged.groupby("patient_id")["treatment_arm"].first().value_counts()
N_RT = int(_split.get("Any_RT", 0))
N_NORT = int(_split.get("No_RT", 0))

with open(INTER / "lmm_interaction_summary.csv") as f:
    rows = list(csv.DictReader(f))

with open(INTER / "lmm_summary.json") as f:
    _lmm_meta = json.load(f)
N_CONVERGED = int(_lmm_meta.get("n_converged", 0))

beta_map = {}
se_map = {}
for r in rows:
    cf = r.get("beta_interaction","").strip()
    se = r.get("se_interaction","").strip()
    beta_map[r["cell_family"]] = float(cf) if cf else np.nan
    se_map[r["cell_family"]] = float(se) if se else np.nan

PRIMARY_CELLS = [
    "CD8_T", "CD4_conv", "Treg_FOXP3", "NK", "B", "Plasma",
    "Microglia_P2RY12_TMEM119", "MoMF_CD163", "DC", "Neutrophils"
]

patients_vec = np.repeat(np.arange(N_PATIENTS), 2)
treatment_vec = np.repeat(np.concatenate([np.ones(N_RT), np.zeros(N_NORT)]), 2)
timepoint_vec = np.tile([0, 1], N_PATIENTS)
df_synth = pd.DataFrame({
    "patient_id": patients_vec,
    "treatment": treatment_vec,
    "timepoint": timepoint_vec,
})

X4 = np.column_stack([
    np.ones(len(df_synth)),
    treatment_vec,
    timepoint_vec,
    treatment_vec * timepoint_vec,
])
ss_int = np.sum((X4[:, 3] - X4[:, 3].mean()) ** 2)
N_OBS = X4.shape[0]
K_COV = 7

ICC_SET = [0.40, 0.60, 0.70]
F2_ANCHORS = [0.28, 0.40]

NSIM = 800
rng = np.random.default_rng(7331)

def sim_run_post(n_patients, n_rt, n_nort, f2, icc, nsim, seed_offset=0):
    loc_rng = np.random.default_rng(7331 + seed_offset)
    p_vec = np.repeat(np.arange(n_patients), 2)
    t_vec = np.repeat(np.concatenate([np.ones(n_rt), np.zeros(n_nort)]), 2)
    tp_vec = np.tile([0, 1], n_patients)
    X = np.column_stack([np.ones(len(p_vec)), t_vec, tp_vec, t_vec * tp_vec])
    ss = np.sum((X[:, 3] - X[:, 3].mean()) ** 2)
    n_obs = X.shape[0]
    beta_int = np.sqrt(f2 * (1.0 - icc) * (n_obs - 4) / ss)
    pvals = np.full(nsim, np.nan)
    for s in range(nsim):
        re = loc_rng.normal(0, sqrt(icc), n_patients)
        re_vec = re[p_vec]
        eps = loc_rng.normal(0, sqrt(1.0 - icc), n_obs)
        y = X[:, 3] * beta_int + re_vec + eps
        tmp = pd.DataFrame({"patient_id": p_vec, "treatment": t_vec,
                            "timepoint": tp_vec, "y": y})
        try:
            fit = MixedLM.from_formula(
                "y ~ treatment * timepoint",
                groups="patient_id", data=tmp
            ).fit(disp=False, reml=False)
            pv = fit.pvalues.get("treatment:timepoint", np.nan)
            pvals[s] = pv
        except Exception:
            pass
    p05 = np.nanmean(pvals < 0.05)
    se_p = sqrt(max(0.0, p05 * (1 - p05) / max(1, nsim)))
    return p05, se_p, beta_int

def compute_f2_obs(beta_obs, se_obs, icc, n_obs=N_OBS, k_cov=K_COV):
    if np.isnan(beta_obs) or np.isnan(se_obs):
        return np.nan
    z_beta = beta_obs / se_obs
    var_beta_asym = se_obs ** 2
    sigma2_resid_approx = var_beta_asym * ss_int / (1.0 - icc)
    if sigma2_resid_approx <= 0:
        return np.nan
    f2_obs = (beta_obs ** 2 / ss_int) * (n_obs - k_cov) / (2.0 * sigma2_resid_approx)
    return f2_obs

rows_out = []
cache = {}

def cached_sim(n_patients, n_rt, n_nort, f2, icc, nsim, tag):
    key = (n_patients, n_rt, n_nort, round(f2, 10), round(icc, 10), nsim, tag)
    if key not in cache:
        se = int(abs(hash(str(key)) % (2**20)))
        cache[key] = sim_run_post(n_patients, n_rt, n_nort, f2, icc, nsim, seed_offset=se)
    return cache[key]

dc_f2obs_by_icc = {}

for idx, cell in enumerate(PRIMARY_CELLS):
    beta = beta_map.get(cell, np.nan)
    se = se_map.get(cell, np.nan)
    has_fit = not (np.isnan(beta) or np.isnan(se))

    row = {"cell_type": cell, "has_observed_fit": has_fit}

    if has_fit:
        row["obs_beta_interaction"] = beta
        row["obs_se_interaction"] = se

    for icc in ICC_SET:
        if has_fit:
            f2_obs = compute_f2_obs(beta, se, icc)
            row[f"f2_obs_icc{int(icc*100):02d}"] = f2_obs
            dc_f2obs_by_icc[(cell, icc)] = f2_obs
            if not np.isnan(f2_obs) and f2_obs > 0:
                pC, seC, _ = cached_sim(N_PATIENTS, N_RT, N_NORT,
                                        f2_obs, icc, NSIM, tag=f"obs_{cell}_{icc}")
            else:
                pC, seC = np.nan, np.nan
            row[f"power_f2obs_icc{int(icc*100):02d}"] = pC
            row[f"power_f2obs_icc{int(icc*100):02d}_se"] = seC

        for f2a in F2_ANCHORS:
            pA, seA, _ = cached_sim(N_PATIENTS, N_RT, N_NORT,
                                    f2a, icc, NSIM, tag=f"anc_{f2a}_{icc}")
            row[f"power_f2{int(f2a*100):03d}_icc{int(icc*100):02d}"] = pA
            row[f"power_f2{int(f2a*100):03d}_icc{int(icc*100):02d}_se"] = seA

    def headline_from_row(row, icc_ref=0.60):
        # A cell type with no observed fit has no observed f^2 to headline.
        # Reporting the f^2=0.40 anchor fallback as "ADEQUATE" here would
        # misattribute a design-level anchor to a cell that was never modelled.
        if not row.get("has_observed_fit", False):
            return "N/A (not estimable)"
        key = f"power_f2obs_icc{int(icc_ref*100):02d}"
        if key in row and not np.isnan(row.get(key, np.nan)):
            val = row[key]
        else:
            key_f2 = f"power_f2040_icc{int(icc_ref*100):02d}"
            val = row.get(key_f2, np.nan)
        if np.isnan(val):
            return "N/A (no converged fit)"
        if val >= 0.80:
            return "ADEQUATE"
        elif val >= 0.60:
            return "BORDERLINE"
        else:
            return "UNDERPOWERED"

    row["headline_label_icc60"] = headline_from_row(row, 0.60)
    row["headline_label_icc40"] = headline_from_row(row, 0.40)
    row["headline_label_icc70"] = headline_from_row(row, 0.70)
    rows_out.append(row)

df_out = pd.DataFrame(rows_out)
df_out.to_csv(OUT_CSV, index=False)

def pct(x):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "N/A"
    return f"{x:.1%}"

print("=" * 90)
print(f" POST-DATA POWER CHECK (n_patients={N_PATIENTS}, Any_RT={N_RT}, No_RT={N_NORT})")
print(f" statsmodels MixedLM REML=False, nsim={NSIM}")
print("=" * 90)

hdr = (f"{'Cell':<22} | {'Fit?':>4} | "
       f"{'Pwr f2=0.28':>12} | {'Pwr f2=0.40':>12} | "
       f"{'Pwr f2-obs':>12} | {'Label ICC60':>16}")
print(hdr)
print("-" * len(hdr))

for r in rows_out:
    print(
        f"{r['cell_type']:<22} | "
        f"{'Y' if r['has_observed_fit'] else 'N':>4} | "
        f"{pct(r.get('power_f2028_icc60', np.nan)):>12} | "
        f"{pct(r.get('power_f2040_icc60', np.nan)):>12} | "
        f"{pct(r.get('power_f2_obs_icc60', np.nan)):>12} | "
        f"{r['headline_label_icc60']:>16}"
    )

print()
print(" PRE-DATA baseline (target n=33 patients, ICC=0.60) for comparison:")
print("   f2=0.28 -> 62.6% BORDERLINE")
print("   f2=0.40 -> 74.4% BORDERLINE")
print()
print(" POST-DATA effective design: n={} patients (full recovery = pre-data target)".format(N_PATIENTS))
print(" ICC sensitivity at ICC={0.40, 0.60, 0.70}; observed ICC not identifiable from a paired 2-timepoint design")
print(" Models converged in lmm_summary.json: n_converged={}".format(N_CONVERGED))
print()

audit = {
    "generated_at": pd.Timestamp.utcnow().isoformat(),
    "n_patients": N_PATIENTS,
    "n_rt": N_RT,
    "n_nort": N_NORT,
    "n_obs": N_OBS,
    "icc_set": ICC_SET,
    "f2_anchors": F2_ANCHORS,
    "nsim": NSIM,
    "seed": 7331,
    "cells_with_observed_fit": [r["cell_type"] for r in rows_out if r["has_observed_fit"]],
    "n_converged_lmm": N_CONVERGED,
    "design_split_note": "Any_RT/No_RT split ({}/{}) derived from the merged modelled set, not the full clinical cohort.".format(N_RT, N_NORT),
    "note_icc_unobserved": ("Observed ICC is not identifiable from a paired 2-timepoint design; "
                            "ICC range assumed from glioma immune literature midpoints 0.40/0.60/0.70."),
}
with open(OUT_JSON, "w") as f:
    json.dump(audit, f, indent=2)

print(f"Wrote: {OUT_CSV}")
print(f"Wrote: {OUT_JSON}")
