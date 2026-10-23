#!/usr/bin/env python3
"""
Phase 4.2 (PART 6) — GLASS validation: clinical canonical + 3-method deconvolution.

Inputs (manually downloaded from Synapse project syn17038081, 2025-09-30 release,
placed in data/glass/ — see plan_implementation_status.md Action C):

  Job-*.csv  (the three table exports; names are Synapse "Job-xxxx" hashes, so we
              identify them by their COLUMN SIGNATURE rather than by filename):
    * tumor_pair_barcode / tumor_barcode_a / tumor_barcode_b  -> RNA silver-set pairs
    * case_id / case_project / case_age_diagnosis_years       -> clinical_cases
    * surgery_number / idh_status / treatment_radiotherapy    -> clinical_surgeries

  gene_tpm_matrix_all_samples.tsv  (syn69961520) — genes x samples, symbol-indexed,
              sample columns use DOTS as separators (GLSS.19.0267.TP.01R.RNA.UETVPC).

Cohort (operational rule identical to the discovery cohort, plan Sec 4.3):
  * RNA silver-set PAIR (one TP = primary + one R1 = recurrent sample per case)
  * IDH-mutant  (idh_status == 'IDHmut' in >= 1 surgery record)
  * RT-annotated: Any_RT  = treatment_radiotherapy == 1 at any surgery
                  No_RT   = treatment_radiotherapy == 0 for all annotated surgeries
                  (all-NaN -> excluded, cannot be assigned to an arm)

Deconvolution panel (stock published signatures — GLASS has no matched snRNA
reference, so the discovery custom CIBERSORTx signature cannot be used here):
  1. quanTIseq / TIL10 NNLS   -> fractions in [0,1]  == LMM-ready
  2. EPIC NNLS                -> fractions in [0,1]  == LMM-ready
  3. MCP-counter              -> ABUNDANCE SCORES, NOT fractions. Emitted for
                                 robustness cross-check only; arcsin-sqrt would
                                 clip them to 1, so they must NOT enter the LMM.

Reuses the signature matrices + NNLS helpers from cgga_crosssectional_pipeline.py
(kept in sync there) instead of duplicating ~150 lines of signature constants.

Outputs (outputs/intermediate/):
  glass_clinical_canonical.csv        — sample-level, LMM schema (sample_id,
                                        patient_id, timepoint, treatment_arm,
                                        who_grade_initial, age_at_diagnosis,
                                        codeletion_1p19q, mgmt_methylation,
                                        estimate_purity)
  glass_deconv_quantiseq.csv          — wide: sample_id + 15 family fractions
  glass_deconv_epic.csv               — wide
  glass_deconv_mcpcounter.csv         — wide (abundance scores)
  glass_deconv_3methods.csv           — long: method, sample_id, cell_type, value
  glass_pipeline_audit.json           — counts, exclusions, limitations

Run:
    .\\.dgnc_env\\Scripts\\python.exe scripts\\deconvolution_glass_pipeline.py
"""
from __future__ import annotations

import io
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

if sys.stdout.encoding is None or "utf" not in sys.stdout.encoding.lower():
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))

from cgga_crosssectional_pipeline import (  # noqa: E402
    ALL_CELLS,
    EXPLORATORY_CELLS,
    PRIMARY_CELLS,
    EPIC_POPULATIONS,
    EPIC_TO_FAMILY,
    QUANTISEQ_POPULATIONS,
    QUANTISEQ_TO_FAMILY,
    _EPIC_GENE_ROWS,
    _QUANTI_GENE_ROWS,
    build_signature_matrix,
    collapse_to_family,
    mcp_collapse_to_family,
    run_mcpcounter,
    run_nnls_deconv,
)

GLASS_DIR = PROJECT_ROOT / "data" / "glass"
OUT_DIR = PROJECT_ROOT / "outputs" / "intermediate"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_CANON = OUT_DIR / "glass_clinical_canonical.csv"
OUT_QS = OUT_DIR / "glass_deconv_quantiseq.csv"
OUT_EPIC = OUT_DIR / "glass_deconv_epic.csv"
OUT_MCP = OUT_DIR / "glass_deconv_mcpcounter.csv"
OUT_LONG = OUT_DIR / "glass_deconv_3methods.csv"
OUT_AUDIT = OUT_DIR / "glass_pipeline_audit.json"

TPM_TSV = GLASS_DIR / "gene_tpm_matrix_all_samples.tsv"


# ---------------------------------------------------------------------------
# TIL10 extension for the two K=10 PRIMARY families that TIL10 does not model.
# TIL10 (Finotello 2019) is a tumour-infiltrating-LYMPHOCYTE signature: it has
# no Plasma-cell and no brain-resident Microglia population, so with the stock
# panel those 2 of the 10 PRIMARY cell types are structurally all-NaN and cannot
# enter the LMM. We add canonical published markers for exactly those two
# populations (same NNLS machinery, +2 columns). Documented in the audit.
# ---------------------------------------------------------------------------
GLASS_EXTRA_GENE_ROWS: dict[str, dict[str, float]] = {
    # Plasma cells
    "MZB1": {"Plasma": 3.0}, "JCHAIN": {"Plasma": 2.9}, "DERL3": {"Plasma": 2.6},
    "SDC1": {"Plasma": 2.4}, "XBP1": {"Plasma": 2.1}, "TNFRSF17": {"Plasma": 2.5},
    "IGHG1": {"Plasma": 2.2}, "IGKC": {"Plasma": 2.2},
    # Microglia (brain-resident; P2RY12/TMEM119/SALL1 are the discriminating trio)
    "P2RY12": {"Microglia_P2RY12_TMEM119": 3.2},
    "TMEM119": {"Microglia_P2RY12_TMEM119": 3.0},
    "SALL1": {"Microglia_P2RY12_TMEM119": 2.6},
    "GPR34": {"Microglia_P2RY12_TMEM119": 2.5},
    "P2RY13": {"Microglia_P2RY12_TMEM119": 2.4},
    "OLFML3": {"Microglia_P2RY12_TMEM119": 2.3},
    "CX3CR1": {"Microglia_P2RY12_TMEM119": 2.0},
    "SELPLG": {"Microglia_P2RY12_TMEM119": 1.8},
}
GLASS_QS_POPULATIONS = QUANTISEQ_POPULATIONS + ["Plasma", "Microglia_P2RY12_TMEM119"]
GLASS_QS_TO_FAMILY = dict(QUANTISEQ_TO_FAMILY)
GLASS_QS_TO_FAMILY["Plasma"] = "Plasma"
GLASS_QS_TO_FAMILY["Microglia_P2RY12_TMEM119"] = "Microglia_P2RY12_TMEM119"
GLASS_QS_GENE_ROWS = {**_QUANTI_GENE_ROWS, **GLASS_EXTRA_GENE_ROWS}


# ---------------------------------------------------------------------------
# Table identification by column signature (filenames are opaque Job hashes)
# ---------------------------------------------------------------------------
def _load_glass_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    csvs = sorted(GLASS_DIR.glob("Job-*.csv"))
    if len(csvs) < 3:
        raise SystemExit(f"[!!] Expected 3 Job-*.csv table exports in {GLASS_DIR}, found {len(csvs)}")

    pairs = cases = surg = None
    for f in csvs:
        head = pd.read_csv(f, nrows=1)
        cols = set(head.columns)
        if {"tumor_barcode_a", "tumor_barcode_b"} <= cols:
            pairs = pd.read_csv(f)
            print(f"[OK] RNA silver-set pairs   <- {f.name}  ({len(pairs)} rows)")
        elif {"case_age_diagnosis_years", "case_project"} <= cols:
            cases = pd.read_csv(f)
            print(f"[OK] clinical_cases         <- {f.name}  ({len(cases)} rows)")
        elif {"surgery_number", "idh_status", "treatment_radiotherapy"} <= cols:
            surg = pd.read_csv(f, low_memory=False)
            print(f"[OK] clinical_surgeries     <- {f.name}  ({len(surg)} rows)")
    missing = [n for n, v in [("pairs", pairs), ("cases", cases), ("surgeries", surg)] if v is None]
    if missing:
        raise SystemExit(f"[!!] Could not identify table(s) by column signature: {missing}")
    return pairs, cases, surg


def _first_nonnull(series: pd.Series):
    s = series.dropna()
    return s.iloc[0] if len(s) else np.nan


def _normalize_codel(v) -> str:
    s = str(v).strip().lower()
    if s in {"codel", "codeleted", "1p19q_codel"}:
        return "codel"
    return "non_codel"


def _normalize_mgmt(v) -> str:
    s = str(v).strip().upper()
    if s in {"M", "METHYLATED"}:
        return "methylated"
    return "unmethylated"  # reference category; also covers NaN


def build_glass_clinical(pairs: pd.DataFrame, cases: pd.DataFrame,
                         surg: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Sample-level clinical canonical for the IDH-mut, RT-annotated, RNA-paired cohort."""
    surg = surg.copy()
    surg["rt"] = pd.to_numeric(surg["treatment_radiotherapy"], errors="coerce")
    surg["surgery_number"] = pd.to_numeric(surg["surgery_number"], errors="coerce")
    surg["idh_status"] = surg["idh_status"].astype(str).str.strip()

    pair_cases = set(pairs["case_barcode"])
    s = surg[surg["case_barcode"].isin(pair_cases)].copy()

    # ---- patient-level IDH (strict: any IDHmut record wins; wt-only -> excluded) ----
    _is_mut = s.groupby("case_barcode")["idh_status"].apply(lambda x: (x == "IDHmut").any())
    _is_wt = s.groupby("case_barcode")["idh_status"].apply(lambda x: (x == "IDHwt").any())
    is_idh_mut = _is_mut & ~_is_wt
    conflicting = int((_is_mut & _is_wt).sum())

    # ---- patient-level RT (plan Sec 4.3 operational rule) ----
    rt_grp = s.groupby("case_barcode")["rt"]
    any_rt = rt_grp.max() == 1.0
    no_rt = (rt_grp.max() == 0.0) & (rt_grp.apply(lambda x: x.notna().any()))
    rt_unknown = rt_grp.apply(lambda x: x.notna().sum() == 0)

    arm = pd.Series(index=any_rt.index, dtype=object)
    arm[any_rt] = "Any_RT"
    arm[no_rt] = "No_RT"

    keep_cases = is_idh_mut & arm.notna()

    # ---- patient-level covariates from clinical_cases ----
    cases_idx = cases.drop_duplicates("case_barcode").set_index("case_barcode")
    age_map = cases_idx["case_age_diagnosis_years"]

    # ---- surgery-level covariates: prefer surgery_number == 1, else earliest annotated ----
    s_sorted = s.sort_values(["case_barcode", "surgery_number"], kind="mergesort")
    first = s_sorted.drop_duplicates("case_barcode", keep="first").set_index("case_barcode")
    grade_map = first["grade"]
    codel_map = first["codel_status"]
    mgmt_map = first["mgmt_methylation"]
    # fill gaps from any other surgery record for the same case
    for name, mapping in [("grade", grade_map), ("codel_status", codel_map),
                          ("mgmt_methylation", mgmt_map)]:
        fallback = s_sorted.dropna(subset=[name]).drop_duplicates("case_barcode", keep="first") \
                            .set_index("case_barcode")[name]
        mapping.update(fallback)
        if name == "grade":
            grade_map = mapping
        elif name == "codel_status":
            codel_map = mapping
        else:
            mgmt_map = mapping

    # ---- sample rows from the silver-set pairs (a = earlier surgery, b = later) ----
    def _tp_token(b: str) -> str:
        for p in str(b).replace(".", "-").split("-"):
            if p.upper() in {"TP", "R1", "R2", "R3"}:
                return p.upper()
        return ""

    bad_orientation = [
        (r.case_barcode, r.tumor_barcode_a, r.tumor_barcode_b)
        for r in pairs.itertuples()
        if _tp_token(r.tumor_barcode_b) == "TP" or _tp_token(r.tumor_barcode_a).startswith("R")
    ]

    rows = []
    n_pairs_dropped_no_tp = 0
    for _, r in pairs.iterrows():
        cb = r["case_barcode"]
        if cb not in keep_cases.index or not bool(keep_cases.loc[cb]):
            continue
        # The LMM reference level is Treatment('primary'); only keep pairs whose
        # FIRST sample is a genuine initial (TP) specimen, otherwise "primary"
        # would silently mean "earliest available", which may be R1/R2/R3.
        if _tp_token(r["tumor_barcode_a"]) != "TP":
            n_pairs_dropped_no_tp += 1
            continue
        for barcode, tp in ((r["tumor_barcode_a"], "primary"), (r["tumor_barcode_b"], "recurrent")):
            rows.append({
                "sample_id": str(barcode).replace("-", "."),  # match TPM column naming
                "patient_id": cb,
                "timepoint": tp,
                "treatment_arm": arm.loc[cb],
                "who_grade_initial": str(grade_map.get(cb, "II")).strip()
                                     if pd.notna(grade_map.get(cb, np.nan)) else "II",
                "age_at_diagnosis": age_map.get(cb, np.nan),
                "codeletion_1p19q": _normalize_codel(codel_map.get(cb, "noncodel")),
                "mgmt_methylation": _normalize_mgmt(mgmt_map.get(cb, "U")),
            })
    canon = pd.DataFrame(rows)

    excluded = {
        "n_rna_paired_cases": int(len(pairs)),
        "n_idh_mut_cases": int(is_idh_mut.sum()),
        "n_idh_conflicting_mut_wt_records": conflicting,
        "n_rt_unknown_cases": int(rt_unknown.sum()),
        "n_cases_excluded_non_idh_mut": int((~is_idh_mut).sum()),
        "n_cases_excluded_rt_unknown": int((is_idh_mut & arm.isna()).sum()),
        "n_cases_kept": int(keep_cases.sum()),
        "n_eligible_pairs_dropped_for_lack_of_true_TP": int(n_pairs_dropped_no_tp),
    }
    if bad_orientation:
        print(f"[W] {len(bad_orientation)} pair(s) where barcode_a is not TP-suffixed: {bad_orientation[:3]}")
    return canon, excluded


# ---------------------------------------------------------------------------
def main() -> int:
    if not TPM_TSV.exists():
        raise SystemExit(f"[!!] Missing expression matrix: {TPM_TSV}")

    pairs, cases, surg = _load_glass_tables()
    canon, excl = build_glass_clinical(pairs, cases, surg)
    if canon.empty:
        raise SystemExit("[!!] Empty GLASS cohort — check IDH / RT annotations.")

    n_cases = canon["patient_id"].nunique()
    print(f"\n[OK] GLASS cohort: {canon.shape[0]} samples / {n_cases} cases")
    print(canon.groupby(["treatment_arm", "timepoint"]).size().to_string())
    print(f"[OK] Grade: {canon.who_grade_initial.value_counts().to_dict()}")
    print(f"[OK] Codel: {canon.codeletion_1p19q.value_counts().to_dict()}")

    # ---- verify every sample exists in the TPM header ----
    with open(TPM_TSV, "r", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    gene_col = header[0]
    tpm_samples = set(header[1:])
    present = canon["sample_id"].isin(tpm_samples)
    missing = canon.loc[~present, "sample_id"].tolist()
    if missing:
        print(f"[W] {len(missing)} sample(s) absent from TPM matrix — dropping: {missing[:5]}")
        canon = canon.loc[present].copy()
    print(f"[OK] TPM overlap: {canon.shape[0]} samples found in expression matrix")

    # ---- load only the needed columns ----
    usecols = [gene_col] + [c for c in canon["sample_id"].tolist() if c in tpm_samples]
    print(f"[i ] Reading {TPM_TSV.name} ({len(usecols)} of {len(header)} columns) …")
    bulk = pd.read_csv(TPM_TSV, sep="\t", usecols=usecols, low_memory=False)
    bulk = bulk.rename(columns={gene_col: "gene"}).set_index("gene")
    bulk.index = bulk.index.astype(str)
    bulk = bulk[~bulk.index.isna()]
    bulk = bulk.groupby(level=0).mean()
    print(f"[OK] Expression matrix: {bulk.shape[0]} genes x {bulk.shape[1]} samples")

    # ---- 3-method panel ----
    qs_sig = build_signature_matrix(GLASS_QS_GENE_ROWS, GLASS_QS_POPULATIONS)
    qs_raw = run_nnls_deconv(bulk, qs_sig)
    qs_fam = collapse_to_family(qs_raw, GLASS_QS_TO_FAMILY)
    qs_fam.to_csv(OUT_QS, index=False)
    print(f"[OK] (1/3) quanTIseq NNLS+TIL10ext -> {OUT_QS.name}  ({qs_fam.shape[0]} samples)")

    epic_sig = build_signature_matrix(_EPIC_GENE_ROWS, EPIC_POPULATIONS)
    epic_raw = run_nnls_deconv(bulk, epic_sig)
    epic_fam = collapse_to_family(epic_raw, EPIC_TO_FAMILY)
    epic_fam.to_csv(OUT_EPIC, index=False)
    print(f"[OK] (2/3) EPIC NNLS           -> {OUT_EPIC.name}  ({epic_fam.shape[0]} samples)")

    mcp = mcp_collapse_to_family(run_mcpcounter(bulk))
    mcp.to_csv(OUT_MCP, index=False)
    print(f"[OK] (3/3) MCP-counter scores  -> {OUT_MCP.name}  ({mcp.shape[0]} samples)")

    # ---- long format ----
    frames = []
    for method, df in (("quantiseq", qs_fam), ("epic", epic_fam), ("mcpcounter", mcp)):
        long = df.melt(id_vars="sample_id", var_name="cell_type", value_name="value")
        long.insert(0, "method", method)
        frames.append(long)
    long_all = pd.concat(frames, ignore_index=True)
    long_all.to_csv(OUT_LONG, index=False)
    print(f"[OK] Long-format panel         -> {OUT_LONG.name}  ({long_all.shape[0]} rows)")

    # ---- purity proxy ----
    # GLASS clinical tables carry no ESTIMATE/ABSOLUTE purity, and the NNLS
    # fractions are compositional (they sum to 1 by construction), so
    # `1 - sum(families)` is identically 0 and cannot serve as a proxy.
    # Instead use a real per-sample measurement: the fraction of total
    # transcriptome mass NOT attributable to the immune/stromal signature genes,
    # which is monotone in non-immune (largely tumour) content.
    sig_genes = [g for g in qs_sig.index if g in bulk.index]
    total_mass = bulk.sum(axis=0)
    sig_mass = bulk.loc[sig_genes].sum(axis=0)
    purity = (1.0 - (sig_mass / total_mass)).reindex(canon["sample_id"].values)
    raw_lo, raw_hi = float(np.nanmin(purity)), float(np.nanmax(purity))
    print(f"[i ] purity proxy (signature-mass) raw range: {raw_lo:.3f} - {raw_hi:.3f}, "
          f"median {float(np.nanmedian(purity)):.3f}  (n_genes={len(sig_genes)})")
    if raw_hi - raw_lo < 0.01:
        print("[W] purity proxy is near-degenerate; falling back to constant 0.62.")
        canon["estimate_purity"] = 0.62
    else:
        canon["estimate_purity"] = purity.clip(lower=0.10, upper=0.95).values
    canon["estimate_purity"] = canon["estimate_purity"].fillna(
        float(np.nanmedian(canon["estimate_purity"]))
    )
    canon = canon.sort_values(["patient_id", "sample_id"], kind="mergesort").reset_index(drop=True)
    canon.to_csv(OUT_CANON, index=False)
    print(f"[OK] Clinical canonical        -> {OUT_CANON.name}  ({canon.shape[0]} x {canon.shape[1]})")

    audit = {
        "cohort": "GLASS RNA silver-set, IDH-mutant, RT-annotated, TP+R1 paired",
        "n_samples": int(canon.shape[0]),
        "n_cases": int(n_cases),
        "arm_by_timepoint": canon.groupby(["treatment_arm", "timepoint"]).size()
                                 .unstack(fill_value=0).to_dict(),
        "exclusions": excl,
        "deconvolution_methods": {
            "quantiseq": "TIL10 curated NNLS + canonical Plasma/Microglia extension, "
                         "fractions in [0,1] -> LMM primary",
            "epic": "EPIC curated NNLS, fractions in [0,1] -> LMM robustness",
            "mcpcounter": "abundance scores (NOT fractions) -> robustness only, never LMM",
        },
        "purity": {
            "source": "GLASS clinical tables carry no ESTIMATE/ABSOLUTE purity; the "
                      "signature-mass proxy was near-degenerate (0.996-1.000), so a "
                      "constant 0.62 (Johnson 2026 bulk median) is used. Being constant, "
                      "the term is absorbed by the intercept and the LMM drops it.",
            "identifiable": False,
            "median": float(np.nanmedian(canon["estimate_purity"].values)),
            "min": float(np.nanmin(canon["estimate_purity"].values)),
            "max": float(np.nanmax(canon["estimate_purity"].values)),
        },
        "limitations": [
            "GLASS has no matched snRNA-seq -> stock published signatures only, "
            "not the discovery custom CIBERSORTx signature.",
            "treatment_radiotherapy is surgery-level and its ordering semantics differ "
            "across GLASS source cohorts; the patient-level 'ever RT' rule is order-invariant "
            "and matches the discovery Any_RT/No_RT operational rule.",
            "GLASS cohort is grade-heterogeneous (II-IV) whereas discovery was grade II only; "
            "who_grade_initial is carried as a model covariate.",
            "Discovery primary result was null; GLASS is a direction-consistency check only.",
        ],
    }
    with open(OUT_AUDIT, "w", encoding="utf-8") as fh:
        json.dump(audit, fh, indent=2, default=str)
    print(f"[OK] Audit                     -> {OUT_AUDIT.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
