#!/usr/bin/env python3
"""
Phase 3.1 — Primary statistical test: interaction LMM.

Per cell type (K=10 primary + K=5 exploratory), fit the locked model from
plan §4.5(b) / pre-spec §4.5(b):

    arcsin_sqrt_fraction ~ timepoint * treatment_arm
                          + estimate_purity
                          + who_grade_initial
                          + age_at_diagnosis
                          + codeletion_1p19q
                          + mgmt_methylation
                          + (1 | patient_id)

Response transform = arcsin(sqrt(fraction)) — avoids logit(0)=-Inf without
pseudo-counting (critical for low-infiltration glioma with many exact-zero
deconvolution fractions).

Hypothesis test on `timepoint:treatment_arm` coefficient (the interaction).
95% CIs via statsmodels profile-likelihood (not Wald — Wald can be
anti-conservative at small n). BH-FDR restricted to K=10 PRIMARY family only;
exploratory 5 always reported with raw p + EXPLORATORY flag.

Outputs (to outputs/intermediate/):
  lmm_all_coefficients.csv   — all model terms × cell type
  lmm_interaction_summary.csv — interaction-focused with BH-FDR
  lmm_summary.json           — run metadata (n_samples, n_patients, converge flags)
"""
from __future__ import annotations

import io
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

if sys.stdout.encoding is None or "utf" not in sys.stdout.encoding.lower():
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PRE_SPEC_PATH = PROJECT_ROOT / "outputs" / "sec4_prespec.json"
DEFAULT_FRACTIONS = PROJECT_ROOT / "outputs" / "intermediate" / "deconv_custom_cibersortx.csv"
DEFAULT_CLINICAL = PROJECT_ROOT / "outputs" / "intermediate" / "clinical_canonical.csv"
DEFAULT_PURITY = PROJECT_ROOT / "outputs" / "intermediate" / "estimate_scores.csv"

OUT_ALL_COEF = PROJECT_ROOT / "outputs" / "intermediate" / "lmm_all_coefficients.csv"
OUT_INTERACTION = PROJECT_ROOT / "outputs" / "intermediate" / "lmm_interaction_summary.csv"
OUT_JSON = PROJECT_ROOT / "outputs" / "intermediate" / "lmm_summary.json"

PRIMARY_CELLS = [
    "CD8_T","CD4_conv","Treg_FOXP3","NK","B","Plasma",
    "Microglia_P2RY12_TMEM119","MoMF_CD163","DC","Neutrophils"
]
EXPLORATORY_CELLS = [
    "T_exhausted_PDCD1_HAVCR2","M1","M2","Endothelial","Fibroblast_pericyte"
]
ALL_CELLS = PRIMARY_CELLS + EXPLORATORY_CELLS

# Cell-name alias map: canonical snRNA-cluster friendly long names produced by
# acquire_snrnaseq_discovery.py annotation -> LMM PRIMARY_CELLS / exploratory short names.
# Applied only when --method-tag starts with "snrna" or flag explicitly set.
# Values not in the map pass through unchanged.
SNRNA_DIRECT_CELL_ALIASES: dict[str, str] = {
    "CD8_T_cells": "CD8_T",
    "CD4_T_cells_conventional": "CD4_conv",
    "T_regulatory_cells_FOXP3": "Treg_FOXP3",
    "NK_cells": "NK",
    "B_cells": "B",
    "Plasma_cells": "Plasma",
    "Microglia_P2RY12_TMEM119": "Microglia_P2RY12_TMEM119",
    "Monocyte_derived_macrophages_CD163": "MoMF_CD163",
    "Dendritic_cells": "DC",
    "Neutrophils": "Neutrophils",
    "Endothelial_cells": "Endothelial",
    "Pericytes_Fibroblasts": "Fibroblast_pericyte",
    "Oligodendrocytes_OPC": "Oligodendrocytes_OPC",
    "Astrocytes_non_malignant": "Astrocytes_non_malignant",
    "Malignant_astrocytes_oligo": "Malignant_astrocytes_oligo",
}
DEFAULT_SNRNA_FRACTIONS = PROJECT_ROOT / "outputs" / "intermediate" / "sn_proportions_per_sample.csv"

# Locked formula from pre-spec §4.5(b), MINUS the mgmt_methylation term.
# 2026-10-19 deviation v1.0: MGMT is absent from all 8 Johnson 2026 supplementary
# sheets and the discovery covariates are no longer placeholder-filled (see
# build_clinical_canonical.py), so the term is dropped rather than fitted to
# random draws. Every other term is unchanged.
FORMULA_FIXED = (
    "y ~ C(timepoint, Treatment('primary')) * C(treatment_arm, Treatment('No_RT'))"
    " + estimate_purity"
    " + C(who_grade_initial, Treatment('II'))"
    " + age_at_diagnosis"
    " + C(codeletion_1p19q, Treatment('non_codel'))"
)


def arcsin_sqrt(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    x = np.clip(x, 0.0, 1.0)
    return np.arcsin(np.sqrt(x))


def fit_with_fallback(md) -> tuple[object, str]:
    """Fit a MixedLM, retrying with gradient-free/quasi-Newton optimizers.

    statsmodels' lbfgs raises `LinAlgError: Singular matrix` for some cohorts
    even when the model is well specified (observed on the GLASS validation
    cohort). bfgs/powell/cg reach the identical optimum, so this is a purely
    numerical failure. lbfgs is tried first so discovery-cohort numerics are
    bit-identical to the original run.
    """
    chain = [("lbfgs", 500), ("bfgs", 2000), ("powell", 2000)]
    last_exc: Exception | None = None
    for method, maxiter in chain:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                res = md.fit(reml=False, method=method, disp=False, maxiter=maxiter)
            return res, method
        except Exception as exc:  # noqa: BLE001 - escalate only if all fail
            last_exc = exc
    raise last_exc if last_exc is not None else RuntimeError("mixedlm fit failed")


def effective_formula(df: pd.DataFrame) -> tuple[str, list[str]]:
    """Drop numeric covariates that are constant in THIS cohort.

    A zero-variance numeric regressor is exactly collinear with the intercept,
    which makes the fixed-effects design rank-deficient and makes MixedLM raise
    LinAlgError: Singular matrix. Omitting the term is statistically identical
    (the constant is absorbed by the intercept). On the discovery cohort both
    covariates vary, so the locked formula is returned unchanged.
    """
    formula = FORMULA_FIXED
    dropped: list[str] = []
    for col in ("estimate_purity", "age_at_diagnosis"):
        if col in df.columns and pd.to_numeric(df[col], errors="coerce").nunique(dropna=True) <= 1:
            formula = formula.replace(f" + {col}", "", 1)
            dropped.append(col)
    return formula, dropped


def bh_fdr_primary_only(summary_df: pd.DataFrame) -> pd.DataFrame:
    mask = summary_df["cell_family"].isin(PRIMARY_CELLS)
    p_primary = summary_df.loc[mask, "asym_p_interaction"].values.astype(float)
    if len(p_primary) > 0 and np.isfinite(p_primary).sum() > 0:
        finite_mask = np.isfinite(p_primary)
        p_finite = p_primary[finite_mask]
        _, q_finite, _, _ = multipletests(p_finite, method="fdr_bh")
        q_primary = np.full(len(p_primary), np.nan)
        q_primary[finite_mask] = q_finite
        summary_df.loc[mask, "asym_q_interaction_bh"] = q_primary
    summary_df["exploratory_flag"] = ~summary_df["cell_family"].isin(PRIMARY_CELLS)
    return summary_df


def safe_profile_ci(model_result, param_name: str, alpha: float = 0.05) -> tuple[float, float]:
    """statsmodels profile confidence interval; fallback to Wald if profile fails."""
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ci = model_result.conf_int(alpha=alpha, method="profile")
        row = ci.loc[param_name].values
        return float(row[0]), float(row[1])
    except Exception:
        try:
            ci = model_result.conf_int(alpha=alpha)
            row = ci.loc[param_name].values
            return float(row[0]), float(row[1])
        except Exception:
            return float("nan"), float("nan")


def interaction_param_name(levels_timepoint: list[str], levels_treat: list[str]) -> str:
    """Match the exact patsy param-name syntax statsmodels generates for the
    interaction term given the baseline levels."""
    tp_levels_sorted = sorted(levels_timepoint, key=lambda s: (0 if s == "primary" else 1, s))
    tr_levels_sorted = sorted(levels_treat, key=lambda s: (0 if s == "No_RT" else 1, s))
    tp_nonbaseline = [x for x in tp_levels_sorted if x != "primary"]
    tr_nonbaseline = [x for x in tr_levels_sorted if x != "No_RT"]
    if not tp_nonbaseline or not tr_nonbaseline:
        return ""
    tp = tp_nonbaseline[0]
    tr = tr_nonbaseline[0]
    name = (
        f"C(timepoint, Treatment('primary'))[T.{tp}]:"
        f"C(treatment_arm, Treatment('No_RT'))[T.{tr}]"
    )
    return name


def merge_inputs(fractions_df: pd.DataFrame, clinical_df: pd.DataFrame,
                 purity_df: pd.DataFrame | None) -> pd.DataFrame:
    df = clinical_df.merge(fractions_df, on="sample_id", how="inner")
    if purity_df is not None and not purity_df.empty and "tumor_purity" in purity_df.columns:
        if "estimate_purity" not in df.columns:
            df = df.merge(purity_df[["sample_id","tumor_purity"]].rename(
                columns={"tumor_purity": "estimate_purity"}
            ), on="sample_id", how="left")
            df["estimate_purity"] = df["estimate_purity"].fillna(0.5)
    # Back-fill required covariates with safe defaults if they are missing;
    # downstream runs get real values from clinical_canonical.csv (grade/subtype/
    # purity/age from Johnson 2026 Supp Tables 1-3). MGMT is not modelled at all
    # (pre-spec deviation v1.0). Missing -> global median / reference category.
    if "estimate_purity" not in df.columns:
        df["estimate_purity"] = 0.5
    if "who_grade_initial" not in df.columns:
        df["who_grade_initial"] = "II"
    df["who_grade_initial"] = df["who_grade_initial"].astype(str).fillna("II")
    if "age_at_diagnosis" not in df.columns:
        df["age_at_diagnosis"] = np.nanmedian(df.get("age_at_diagnosis", 45.0))
    df["age_at_diagnosis"] = pd.to_numeric(df["age_at_diagnosis"], errors="coerce")
    df["age_at_diagnosis"] = df["age_at_diagnosis"].fillna(float(np.nanmedian(df["age_at_diagnosis"])))
    if "codeletion_1p19q" not in df.columns:
        df["codeletion_1p19q"] = "non_codel"
    df["codeletion_1p19q"] = df["codeletion_1p19q"].astype(str).fillna("non_codel")
    # mgmt_methylation: intentionally not created (pre-spec deviation v1.0).
    if "timepoint" not in df.columns:
        raise ValueError("clinical_canonical.csv must have column `timepoint` with values 'primary' and 'recurrent'.")
    df["timepoint"] = df["timepoint"].astype(str)
    # Normalize treatment labels: plan §4.3 = Any_RT vs No_RT.
    if "treatment_arm" not in df.columns:
        raise ValueError("clinical_canonical.csv must have column `treatment_arm`.")
    ta = df["treatment_arm"].astype(str).str.strip()
    ta = ta.replace({"Any_RT":"Any_RT","Any RT":"Any_RT","ANY_RT":"Any_RT","RT":"Any_RT",
                     "No_RT":"No_RT","No RT":"No_RT","NO_RT":"No_RT"})
    df["treatment_arm"] = ta
    if "patient_id" not in df.columns:
        raise ValueError("clinical_canonical.csv must have column `patient_id`.")
    return df


def main() -> int:
    if not PRE_SPEC_PATH.exists():
        print("[!!] Rule: never run analysis before the test list is locked.")
        print(f"     Missing {PRE_SPEC_PATH}. Run scripts/write_sec4_prespec.py first.")
        return 2

    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--fractions", type=str, default=None,
                    help="Override fractions CSV (default: custom CIBERSORTx)")
    ap.add_argument("--clinical", type=str, default=None)
    ap.add_argument("--purity", type=str, default=None)
    ap.add_argument("--method-tag", type=str, default="custom_cibersortx",
                    help="Method label written into output JSON (e.g. quantiseq / epic / mcpcounter)")
    ap.add_argument("--out-suffix", type=str, default="",
                    help="Suffix appended to the 3 output filenames so a validation-cohort "
                         "run does not overwrite the discovery outputs (e.g. --out-suffix _glass)")
    args = ap.parse_args()

    # Redirect outputs when a suffix is given (discovery files must stay intact).
    global OUT_ALL_COEF, OUT_INTERACTION, OUT_JSON
    if args.out_suffix:
        OUT_ALL_COEF = OUT_ALL_COEF.with_name(OUT_ALL_COEF.stem + args.out_suffix + OUT_ALL_COEF.suffix)
        OUT_INTERACTION = OUT_INTERACTION.with_name(OUT_INTERACTION.stem + args.out_suffix + OUT_INTERACTION.suffix)
        OUT_JSON = OUT_JSON.with_name(OUT_JSON.stem + args.out_suffix + OUT_JSON.suffix)
        print(f"[i ] --out-suffix {args.out_suffix!r} -> {OUT_INTERACTION.name}")

    method_is_snrna = args.method_tag.lower().startswith("snrna")

    frac_path = Path(args.fractions) if args.fractions else (
        DEFAULT_SNRNA_FRACTIONS if (method_is_snrna and DEFAULT_SNRNA_FRACTIONS.exists()) else DEFAULT_FRACTIONS
    )
    clin_path = Path(args.clinical) if args.clinical else DEFAULT_CLINICAL
    pur_path = Path(args.purity) if args.purity else (
        DEFAULT_PURITY if DEFAULT_PURITY.exists() else None
    )

    if not args.fractions and method_is_snrna:
        print(f"[i ] snrna_direct mode detected -> fractions default override: {frac_path.name}")

    if not frac_path.exists() or not clin_path.exists():
        print("[--] Required inputs missing. Expected:")
        print(f"     Fractions: {frac_path}  (snrna_direct auto-defaulted={method_is_snrna and not args.fractions})")
        print(f"     Clinical:  {clin_path}")
        print(f"     Purity:    {DEFAULT_PURITY} (optional, back-filled to median if absent)")
        print("[i ] This script runs after: (a) acquire_clinical_metadata.py ran,")
        print("     (b) deconvolution produced fraction CSVs, (c) estimate_purity.py ran.")
        pd.DataFrame(columns=["cell_family","term","estimate","se","z","p","ci_low","ci_high","converged","method_tag"]
                     ).to_csv(OUT_ALL_COEF, index=False)
        pd.DataFrame(columns=[
            "cell_family","beta_interaction","se_interaction","ci_low_interaction",
            "ci_high_interaction","asym_p_interaction","asym_q_interaction_bh",
            "exploratory_flag","converged","method_tag"
        ]).to_csv(OUT_INTERACTION, index=False)
        with open(OUT_JSON, "w", encoding="utf-8") as fh:
            json.dump({"status":"placeholder_required_inputs_missing",
                       "method_tag":args.method_tag}, fh, indent=2)
        print("[OK] Empty placeholders written (3 files).")
        return 0

    fractions_df = pd.read_csv(frac_path)
    clinical_df = pd.read_csv(clin_path)
    purity_df = pd.read_csv(pur_path) if pur_path is not None and pur_path.exists() else None
    if "sample_id" not in fractions_df.columns:
        fractions_df = fractions_df.rename(columns={fractions_df.columns[0]: "sample_id"})

    # Apply snRNA long-name -> LMM short-name alias map if snrna mode or aliases hit.
    applied_aliases: list[tuple[str, str]] = []
    rename_map: dict[str, str] = {}
    for src_col, dst_col in SNRNA_DIRECT_CELL_ALIASES.items():
        if src_col in fractions_df.columns and dst_col not in fractions_df.columns and src_col != dst_col:
            rename_map[src_col] = dst_col
            applied_aliases.append((src_col, dst_col))
        elif src_col in fractions_df.columns and src_col == dst_col:
            pass  # identity
    if applied_aliases:
        fractions_df = fractions_df.rename(columns=rename_map)
        print(f"[OK] snRNA cell-name alias map applied ({len(applied_aliases)} columns): " +
              ", ".join(f"{s}→{d}" for s, d in applied_aliases[:6]) +
              ("..." if len(applied_aliases) > 6 else ""))

    merged = merge_inputs(fractions_df, clinical_df, purity_df)
    n_samples = len(merged)
    n_patients = merged["patient_id"].nunique()
    print(f"[OK] Merged data: n_samples={n_samples}, n_patients={n_patients}")
    print(f"[OK] Timepoints: {sorted(merged['timepoint'].unique().tolist())}")
    print(f"[OK] Treatment arms: {sorted(merged['treatment_arm'].unique().tolist())}")

    tp_levels = sorted(merged["timepoint"].unique().tolist())
    treat_levels = sorted(merged["treatment_arm"].unique().tolist())
    ip_name = interaction_param_name(tp_levels, treat_levels)
    if not ip_name:
        print(f"[!!] Cannot determine interaction param name from tp={tp_levels}, treat={treat_levels}")
        return 4

    fit_formula, dropped_covariates = effective_formula(merged)
    if dropped_covariates:
        print(f"[W] Zero-variance covariate(s) dropped from this cohort (collinear with "
              f"intercept): {dropped_covariates}. Locked formula otherwise unchanged.")
        print(f"[i ] Effective formula: {fit_formula}")
    if merged["estimate_purity"].nunique() <= 1:
        print("[W] NOTE: estimate_purity is constant here -> purity adjustment is NOT "
              "applied in this cohort (no purity estimate available).")

    all_coef_rows = []
    interaction_rows = []
    conv_flags = {}
    optimizers_used: dict[str, str] = {}

    for cell in ALL_CELLS:
        if cell not in merged.columns:
            interaction_rows.append({
                "cell_family": cell,
                "beta_interaction": np.nan,
                "se_interaction": np.nan,
                "ci_low_interaction": np.nan,
                "ci_high_interaction": np.nan,
                "asym_p_interaction": np.nan,
                "asym_q_interaction_bh": np.nan,
                "exploratory_flag": cell not in PRIMARY_CELLS,
                "converged": False,
                "method_tag": args.method_tag,
            })
            continue
        work = merged[["patient_id","timepoint","treatment_arm",
                       "estimate_purity","who_grade_initial","age_at_diagnosis",
                       "codeletion_1p19q", cell]].copy()
        work = work.dropna(subset=[cell])
        work["y"] = arcsin_sqrt(work[cell].values)
        if work["y"].var(ddof=0) < 1e-16 or len(work) < 6:
            all_coef_rows.append({
                "cell_family": cell, "term": "singular_dropped",
                "estimate": np.nan, "se": np.nan, "z": np.nan, "p": np.nan,
                "ci_low": np.nan, "ci_high": np.nan,
                "converged": False, "method_tag": args.method_tag,
            })
            interaction_rows.append({
                "cell_family": cell,
                "beta_interaction": np.nan, "se_interaction": np.nan,
                "ci_low_interaction": np.nan, "ci_high_interaction": np.nan,
                "asym_p_interaction": np.nan, "asym_q_interaction_bh": np.nan,
                "exploratory_flag": cell not in PRIMARY_CELLS,
                "converged": False, "method_tag": args.method_tag,
            })
            continue
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                md = smf.mixedlm(fit_formula, work, groups=work["patient_id"], re_formula="~1")
                mdf, opt_used = fit_with_fallback(md)
            if opt_used != "lbfgs":
                print(f"[i ] {cell:30s}  fit recovered with optimizer '{opt_used}' "
                      f"(lbfgs raised LinAlgError)")
                optimizers_used[cell] = opt_used
            converged = bool(getattr(mdf, "converged", True)) and not bool(getattr(mdf, "mle_retvals", {}).get("warnflag", 0) if isinstance(getattr(mdf, "mle_retvals", None), dict) else 0)
            conv_flags[cell] = converged
            # All terms coef table
            for term, coef in mdf.params.items():
                se = mdf.bse.get(term, np.nan)
                z = coef / se if se and np.isfinite(se) and se > 0 else np.nan
                pval = mdf.pvalues.get(term, np.nan)
                ci_l, ci_h = safe_profile_ci(mdf, term)
                all_coef_rows.append({
                    "cell_family": cell, "term": term,
                    "estimate": float(coef), "se": float(se), "z": float(z), "p": float(pval),
                    "ci_low": float(ci_l), "ci_high": float(ci_h),
                    "converged": converged, "method_tag": args.method_tag,
                })
            # Interaction term focus
            beta_int = float(mdf.params.get(ip_name, np.nan))
            se_int = float(mdf.bse.get(ip_name, np.nan))
            p_int = float(mdf.pvalues.get(ip_name, np.nan))
            ci_lo, ci_hi = safe_profile_ci(mdf, ip_name)
            interaction_rows.append({
                "cell_family": cell,
                "beta_interaction": beta_int,
                "se_interaction": se_int,
                "ci_low_interaction": ci_lo,
                "ci_high_interaction": ci_hi,
                "asym_p_interaction": p_int,
                "asym_q_interaction_bh": np.nan,  # filled below
                "exploratory_flag": cell not in PRIMARY_CELLS,
                "converged": converged,
                "method_tag": args.method_tag,
            })
            print(f"[OK] {cell:30s}  beta_int={beta_int: .4f}  p={p_int:.4f}  conv={converged}")
        except Exception as exc:
            print(f"[--] {cell:30s}  LMM fit failed: {type(exc).__name__}: {exc}")
            interaction_rows.append({
                "cell_family": cell,
                "beta_interaction": np.nan, "se_interaction": np.nan,
                "ci_low_interaction": np.nan, "ci_high_interaction": np.nan,
                "asym_p_interaction": np.nan, "asym_q_interaction_bh": np.nan,
                "exploratory_flag": cell not in PRIMARY_CELLS,
                "converged": False, "method_tag": args.method_tag,
            })
            conv_flags[cell] = False

    all_coef_df = pd.DataFrame(all_coef_rows)
    interaction_df = pd.DataFrame(interaction_rows)
    interaction_df = bh_fdr_primary_only(interaction_df)

    OUT_ALL_COEF.parent.mkdir(parents=True, exist_ok=True)
    all_coef_df.to_csv(OUT_ALL_COEF, index=False)
    interaction_df.to_csv(OUT_INTERACTION, index=False)
    with open(OUT_JSON, "w", encoding="utf-8") as fh:
        json.dump({
            "status": "completed",
            "method_tag": args.method_tag,
            "n_samples": n_samples,
            "n_patients": n_patients,
            "formula_fixed": FORMULA_FIXED,
            "formula_effective": fit_formula,
            "covariates_dropped_zero_variance": dropped_covariates,
            "response_transform": "arcsin_sqrt",
            "interaction_param_name": ip_name,
            "bh_fdr_family": "K=10 PRIMARY ONLY, exploratory excluded from correction",
            "n_converged": sum(1 for v in conv_flags.values() if v),
            "convergence_flags": conv_flags,
            "optimizers_used": optimizers_used,
        }, fh, indent=2, default=str)
    print(f"\n[OK] Wrote {OUT_ALL_COEF}")
    print(f"[OK] Wrote {OUT_INTERACTION}")
    print(f"[OK] Wrote {OUT_JSON}")
    print("\nInteraction summary (PRIMARY family BH-FDR q < 0.05 highlighted):")
    print(interaction_df[["cell_family","beta_interaction","ci_low_interaction","ci_high_interaction",
                          "asym_p_interaction","asym_q_interaction_bh","exploratory_flag","converged"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
