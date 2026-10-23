#!/usr/bin/env python3
"""
Phase 3.2 — Permutation null for the LMM interaction test.

Operationalization (locked per pre-spec §4.5):
  * Shuffle ONLY at the patient level: within each shuffled assignment, a
    patient keeps the same shuffled treatment_arm for BOTH timepoints so the
    paired (primary, recurrent) correlation structure is preserved exactly.
  * n_perm = 10000 iterations (locked)
  * np.random.seed(7331) global; each joblib worker uses child-seeded rng from
    the parent sequence to preserve exact reproducibility.
  * For each cell type present in lmm_interaction_summary.csv with a
    non-NaN observed interaction beta, re-fit the EXACT same MixedLM model
    (formula + arcsin_sqrt transform + REML=False + lbfgs + maxiter=500) on
    every permuted treatment column. Extract interaction beta_i^perm.
  * Two-tailed EMPIRICAL p-value:
        perm_p = (1 + sum_{i=1..N} I(|beta_i^perm| >= |beta^obs|))
                 / (1 + N)
    The +1 in numerator/denominator avoids p = 0 (Phipson & Smyth 2010).

Outputs:
  outputs/intermediate/permutation_null_all_betas.npz   (all permuted betas
     per cell type; large; retained for re-use if alpha changes)
  outputs/intermediate/snrna_direct_perm_pvalues.csv    (summary table w/
     obs_beta + perm_p + diagnostic fields)
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
from joblib import Parallel, delayed

from lmm_interaction_pipeline import (
    FORMULA_FIXED,
    PRIMARY_CELLS,
    SNRNA_DIRECT_CELL_ALIASES,
    arcsin_sqrt,
    effective_formula,
    fit_with_fallback,
    interaction_param_name,
    merge_inputs,
)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = PROJECT_ROOT / "outputs" / "intermediate"
OUT_DIR.mkdir(parents=True, exist_ok=True)

NPZ_PATH = OUT_DIR / "permutation_null_all_betas.npz"
CSV_PATH = OUT_DIR / "snrna_direct_perm_pvalues.csv"
CKPT_DIR = OUT_DIR / "perm_checkpoints"

DEFAULT_CLINICAL = OUT_DIR / "clinical_canonical.csv"
DEFAULT_SNRNA_FRACTIONS = OUT_DIR / "sn_proportions_per_sample.csv"
DEFAULT_PURITY = OUT_DIR / "estimate_scores.csv"
LMM_SUMMARY_CSV = OUT_DIR / "lmm_interaction_summary.csv"

N_PERM = 10_000
SEED = 7331
N_JOBS = 4  # 4 physical cores, ~4 GB free RAM (each worker's working set is a 71-row frame). Checkpoints make this resumable.


def load_inputs_snrna() -> pd.DataFrame:
    fractions_df = pd.read_csv(DEFAULT_SNRNA_FRACTIONS)
    if "sample_id" not in fractions_df.columns:
        fractions_df = fractions_df.rename(columns={fractions_df.columns[0]: "sample_id"})
    applied_aliases = []
    rename_map = {}
    for src, dst in SNRNA_DIRECT_CELL_ALIASES.items():
        if src in fractions_df.columns and dst not in fractions_df.columns and src != dst:
            rename_map[src] = dst
            applied_aliases.append((src, dst))
    if applied_aliases:
        fractions_df = fractions_df.rename(columns=rename_map)
        print(f"[OK] snRNA aliases: {len(applied_aliases)} applied")

    clinical_df = pd.read_csv(DEFAULT_CLINICAL)
    purity_df = pd.read_csv(DEFAULT_PURITY) if DEFAULT_PURITY.exists() else None
    merged = merge_inputs(fractions_df, clinical_df, purity_df)
    print(f"[OK] Loaded merged data: n_samples={len(merged)}, n_patients={merged['patient_id'].nunique()}")
    return merged


def get_target_cells() -> list[tuple[str, float]]:
    """Return list of (cell_name, obs_beta) for non-NaN cells in summary."""
    summ = pd.read_csv(LMM_SUMMARY_CSV)
    out: list[tuple[str, float]] = []
    for _, row in summ.iterrows():
        cell = row["cell_family"]
        b = row["beta_interaction"]
        if pd.notna(b) and np.isfinite(b):
            out.append((cell, float(b)))
    print(f"[OK] Target cells for permutation: {len(out)} "
          + ", ".join(f"{c} (obs_b={b:.4f})" for c, b in out))
    return out


def fit_one_perm(merged: pd.DataFrame, ip_name: str, cell: str, formula: str,
                 shuffled_treatment: dict[str, str]) -> float:
    work = merged[["patient_id","timepoint","treatment_arm","estimate_purity",
                    "who_grade_initial","age_at_diagnosis","codeletion_1p19q",
                    cell]].copy()
    work["treatment_arm"] = work["patient_id"].map(shuffled_treatment).fillna(work["treatment_arm"])
    work[cell] = work[cell].fillna(0.0)
    work["y"] = arcsin_sqrt(work[cell].values)
    if work["y"].var(ddof=0) < 1e-16 or len(work) < 6:
        return float("nan")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            md = smf.mixedlm(formula, work, groups=work["patient_id"], re_formula="~1")
            # Same lbfgs-first optimizer-escalation chain as the observed fit
            # (lmm_interaction_pipeline.fit_with_fallback). lbfgs alone raises
            # LinAlgError for several discovery cell types even though the model
            # is well specified; bfgs reaches the identical optimum.
            mdf, _ = fit_with_fallback(md)
        b = float(mdf.params.get(ip_name, np.nan))
        return b
    except Exception:
        return float("nan")


def _ckpt_path(cell: str) -> Path:
    return CKPT_DIR / f"perm_ckpt_{cell}.npz"


def load_checkpoint(cell: str) -> np.ndarray | None:
    """Return the cached beta array for a cell if it matches seed/n_perm, else None.

    Partitioned-execution hard constraint: a crash mid-run must not force a full
    restart. Each completed cell is persisted independently and validated against
    the locked SEED/N_PERM before being trusted.
    """
    path = _ckpt_path(cell)
    if not path.exists():
        return None
    try:
        with np.load(path) as data:
            arr = data["arr"]
            ckpt_seed = int(data["seed"][0])
            ckpt_n = int(data["n_perm"][0])
        if ckpt_seed != SEED or ckpt_n != N_PERM or arr.shape[0] != N_PERM:
            print(f"[W] Stale checkpoint for {cell} (seed={ckpt_seed}, n={ckpt_n}); recomputing.")
            return None
        return arr
    except Exception as exc:  # noqa: BLE001 - corrupt checkpoint -> recompute
        print(f"[W] Unreadable checkpoint for {cell} ({exc}); recomputing.")
        return None


def save_checkpoint(cell: str, arr: np.ndarray) -> None:
    CKPT_DIR.mkdir(parents=True, exist_ok=True)
    np.savez(_ckpt_path(cell), arr=arr,
             seed=np.array([SEED]), n_perm=np.array([N_PERM]))


def summarize_cell(cell: str, arr: np.ndarray, obs_b: float) -> dict:
    n_ok = int(np.isfinite(arr).sum())
    if n_ok == 0:
        perm_p = float("nan")
        cnt = 0
    else:
        finite = arr[np.isfinite(arr)]
        cnt = int(np.sum(np.abs(finite) >= abs(obs_b)))
        perm_p = (1.0 + cnt) / (1.0 + len(finite))
    return {
        "cell_family": cell,
        "is_primary": cell in PRIMARY_CELLS,
        "obs_beta": obs_b,
        "perm_p": perm_p,
        "n_perm_ok": n_ok,
        "n_perm_total": N_PERM,
        "seed": SEED,
    }


def main() -> int:
    if sys.stdout.encoding is None or "utf" not in sys.stdout.encoding.lower():
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

    if not LMM_SUMMARY_CSV.exists() or not DEFAULT_SNRNA_FRACTIONS.exists() or not DEFAULT_CLINICAL.exists():
        print("[!!] Missing inputs. Run lmm_interaction_pipeline.py --method-tag snrna_direct first.")
        return 2

    merged = load_inputs_snrna()
    targets = get_target_cells()
    if not targets:
        print("[--] No cell types with observed beta. Writing empty CSV + done.")
        pd.DataFrame(columns=["cell_family","is_primary","obs_beta","perm_p","n_perm_ok","n_perm_total"]
                     ).to_csv(CSV_PATH, index=False)
        np.savez(NPZ_PATH)
        return 0

    tp_levels = sorted(merged["timepoint"].unique().tolist())
    treat_levels = sorted(merged["treatment_arm"].unique().tolist())
    ip_name = interaction_param_name(tp_levels, treat_levels)
    print(f"[OK] Interaction param name = {ip_name}")

    # Same zero-variance-covariate handling as the observed fit.
    fit_formula, dropped = effective_formula(merged)
    if dropped:
        print(f"[W] Zero-variance covariate(s) dropped: {dropped} -> {fit_formula}")
    else:
        print("[OK] Effective formula == locked formula (no zero-variance covariate)")

    patients = sorted(merged["patient_id"].unique().tolist())
    n_patients = len(patients)
    orig_treat = merged.groupby("patient_id")["treatment_arm"].first().to_dict()
    orig_treat_arr = np.array([orig_treat[p] for p in patients])

    rng = np.random.default_rng(SEED)
    perm_seeds = rng.integers(0, 2**31 - 1, size=N_PERM)
    print(f"[OK] Seeded rng seed={SEED}; n_perm={N_PERM}; n_jobs={N_JOBS}")

    all_betas: dict[str, np.ndarray] = {}
    rows = []

    for cell, obs_b in targets:
        cached = load_checkpoint(cell)
        if cached is not None:
            arr = cached
            print(f"\n[==] {cell}: resume from checkpoint (n_ok={int(np.isfinite(arr).sum())}/{N_PERM})")
        else:
            print(f"\n[>>] Permutation for {cell} (obs_beta={obs_b:.4f}) ...")
            def task(i: int) -> float:
                prng = np.random.default_rng(perm_seeds[i])
                shuffled_idx = prng.permutation(n_patients)
                new_vals = orig_treat_arr[shuffled_idx]
                new_treat_dict = {patients[k]: new_vals[k] for k in range(n_patients)}
                return fit_one_perm(merged, ip_name, cell, fit_formula, new_treat_dict)

            results = Parallel(n_jobs=N_JOBS, verbose=1, prefer="processes")(
                delayed(task)(i) for i in range(N_PERM)
            )
            arr = np.array(results, dtype=np.float64)
            save_checkpoint(cell, arr)
            print(f"[OK] {cell}: checkpoint saved -> {_ckpt_path(cell).name}")

        all_betas[cell] = arr
        row = summarize_cell(cell, arr, obs_b)
        print(f"[OK] {cell}: n_ok={row['n_perm_ok']}/{N_PERM} | obs_b={obs_b:.4f} | "
              f"perm_p={row['perm_p']:.6f}")
        rows.append(row)

        # Incremental CSV write: survives a crash on a later cell.
        pd.DataFrame(rows).to_csv(CSV_PATH, index=False)

    out_df = pd.DataFrame(rows)
    out_df.to_csv(CSV_PATH, index=False)
    print(f"\n[OK] Wrote {CSV_PATH}")

    np.savez(NPZ_PATH, **{c: all_betas[c] for c in all_betas},
             _seed=np.array([SEED]),
             _n_perm=np.array([N_PERM]))
    print(f"[OK] Wrote {NPZ_PATH} ({len(all_betas)} cell beta arrays saved)")

    print("\n== Permutation summary ==")
    with pd.option_context("display.max_columns", None, "display.width", 120):
        print(out_df.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
