#!/usr/bin/env python3
"""
Phase 2.2 Step 1 — Custom-reference CIBERSORTx deconvolution (primary method).

Two operational modes:

MODE A (post-snRNA pipeline, preferred):
    Signature matrix built from discovery snRNA-seq cluster markers.
    Signature rows = genes (HGNC symbol), cols = cell populations.
    Values are log2(mean_cluster_expr / mean_rest_expr) or similar —
    NNLS operates on the relative expression scale; sign flipped to positive
    values via abs() on the signature (NNLS requires non-negative coeffs).

MODE B (interim before snRNA clusters finalized, or direct snRNA proportions):
    Read per-sample snRNA proportions from
    outputs/intermediate/snrnaseq_cluster_proportions.csv (produced by
    scripts/acquire_snrnaseq_discovery.py) and map cluster column names
    onto the locked cell-type family names.

Primary NNLS: scipy.optimize.nnls on z-scored expression.
Optional L2 ridge (CIBERSORTx nu parameter): augmented system
    [sig; sqrt(lambda)*I] * f = [expr; 0]
then NNLS on the augmented problem.

Output: outputs/intermediate/deconv_custom_cibersortx.csv
        cols = sample_id + K_10_primary + K_5_exploratory
        Values = fractions (non-negative, row-sum <= 1; not forced to sum 1 so
        missing/unmodeled cell types leave residual sum-headroom.)

Locked cell-type family (pre-spec §4.5(a)), K=10 PRIMARY (BH-FDR across these):
  1 CD8_T
  2 CD4_conv
  3 Treg_FOXP3
  4 NK
  5 B
  6 Plasma
  7 Microglia_P2RY12_TMEM119
  8 MoMF_CD163
  9 DC
 10 Neutrophils

K=5 EXPLORATORY (raw p only, label EXPLORATORY):
 11 T_exhausted_PDCD1_HAVCR2
 12 M1
 13 M2
 14 Endothelial
 15 Fibroblast_pericyte
"""
from __future__ import annotations

import io
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls as _nnls

if sys.stdout.encoding is None or "utf" not in sys.stdout.encoding.lower():
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PRE_SPEC_PATH = PROJECT_ROOT / "outputs" / "sec4_prespec.json"
DEFAULT_BULK_GZ = PROJECT_ROOT / "outputs" / "intermediate" / "bulk_log2tpm_matrix.csv.gz"
DEFAULT_BULK_PQ = PROJECT_ROOT / "outputs" / "intermediate" / "bulk_log2tpm_matrix.parquet"
DEFAULT_SNRNA_PROP = PROJECT_ROOT / "outputs" / "intermediate" / "snrnaseq_cluster_proportions.csv"
DEFAULT_SIGMAT = PROJECT_ROOT / "outputs" / "intermediate" / "custom_signature_matrix.csv"
DEFAULT_OUTPUT = PROJECT_ROOT / "outputs" / "intermediate" / "deconv_custom_cibersortx.csv"

PRIMARY_CELLS = [
    "CD8_T","CD4_conv","Treg_FOXP3","NK","B","Plasma",
    "Microglia_P2RY12_TMEM119","MoMF_CD163","DC","Neutrophils"
]
EXPLORATORY_CELLS = [
    "T_exhausted_PDCD1_HAVCR2","M1","M2","Endothelial","Fibroblast_pericyte"
]
ALL_CELLS = PRIMARY_CELLS + EXPLORATORY_CELLS

# Default mapping from generic cluster labels (Leiden cluster_N or immune cell
# names) onto the locked family. Used in Mode B. User can override by passing
# a JSON mapping via --cluster-map CLI.
DEFAULT_CLUSTER_FAMILY_MAP = {
    "CD8": "CD8_T", "CD8_T": "CD8_T", "CD8T": "CD8_T", "CD8_positive": "CD8_T",
    "CD4": "CD4_conv", "CD4_T": "CD4_conv", "CD4T": "CD4_conv",
    "Treg": "Treg_FOXP3", "Treg_FOXP3": "Treg_FOXP3", "T_regulatory": "Treg_FOXP3", "FOXP3": "Treg_FOXP3",
    "NK": "NK", "NK_cell": "NK", "NK_cells": "NK",
    "B": "B", "B_cell": "B", "B_cells": "B",
    "Plasma": "Plasma", "Plasma_cell": "Plasma", "Plasmablast": "Plasma",
    "Microglia": "Microglia_P2RY12_TMEM119", "Microglia_P2RY12": "Microglia_P2RY12_TMEM119",
    "P2RY12": "Microglia_P2RY12_TMEM119", "TMEM119": "Microglia_P2RY12_TMEM119",
    "Macrophage": "MoMF_CD163", "MoMF": "MoMF_CD163", "CD163": "MoMF_CD163",
    "Mono": "MoMF_CD163", "Monocyte": "MoMF_CD163",
    "DC": "DC", "Dendritic": "DC", "Dendritic_cell": "DC",
    "Neutrophil": "Neutrophils", "Neutrophils": "Neutrophils", "Neut": "Neutrophils",
    "Tex": "T_exhausted_PDCD1_HAVCR2", "T_exhausted": "T_exhausted_PDCD1_HAVCR2",
    "PDCD1": "T_exhausted_PDCD1_HAVCR2", "HAVCR2": "T_exhausted_PDCD1_HAVCR2",
    "M1": "M1", "M1_macrophage": "M1",
    "M2": "M2", "M2_macrophage": "M2",
    "Endothelial": "Endothelial", "Endo": "Endothelial", "Vascular": "Endothelial",
    "Fibroblast": "Fibroblast_pericyte", "Pericyte": "Fibroblast_pericyte",
    "Fibro": "Fibroblast_pericyte", "Stromal": "Fibroblast_pericyte",
}


def load_expr(path: Path) -> pd.DataFrame:
    suf = path.suffix.lower()
    if suf in {".parquet", ".pq"}:
        return pd.read_parquet(path)
    if suf in {".gz", ".csv"}:
        kws = {"index_col": 0, "sep": ","}
        if suf == ".gz":
            kws["compression"] = "gzip"
        return pd.read_csv(path, **kws)
    raise ValueError(f"Unknown suffix {suf}")


def zscore_rows(mat: np.ndarray) -> np.ndarray:
    mu = mat.mean(axis=1, keepdims=True)
    sd = mat.std(axis=1, ddof=1, keepdims=True)
    sd[sd < 1e-12] = 1.0
    return (mat - mu) / sd


def nnls_with_l2(A: np.ndarray, b: np.ndarray, l2_lambda: float = 0.0):
    """scipy.optimize.nnls with optional L2 ridge via augmented system."""
    if l2_lambda <= 0:
        x, rnorm = _nnls(A, b)
        return x, rnorm
    K = A.shape[1]
    A_aug = np.vstack([A, np.sqrt(l2_lambda) * np.eye(K)])
    b_aug = np.concatenate([b, np.zeros(K)])
    x, rnorm = _nnls(A_aug, b_aug)
    return x, rnorm


def run_deconv(bulk: pd.DataFrame, sigmat: pd.DataFrame, l2_lambda: float = 0.0) -> pd.DataFrame:
    bulk_genes = set(bulk.index.astype(str))
    sig_genes = set(sigmat.index.astype(str))
    common = sorted(bulk_genes & sig_genes)
    if len(common) < 50:
        print(f"[!!] Only {len(common)} signature genes overlap with bulk (need >=50).")
        print("     Check gene-symbol space (HGNC vs Ensembl).")
    B = bulk.loc[common].values.astype(np.float64)
    S = sigmat.loc[common].values.astype(np.float64)
    Bz = zscore_rows(B)
    Sz = zscore_rows(S)

    samples = list(bulk.columns)
    cells = list(sigmat.columns)
    out = np.zeros((len(samples), len(cells)), dtype=np.float64)
    for j in range(len(samples)):
        x, _ = nnls_with_l2(Sz, Bz[:, j], l2_lambda=l2_lambda)
        out[j, :] = x
    df = pd.DataFrame(out, columns=cells, index=samples)
    df.index.name = "sample_id"
    return df.reset_index()


def collapse_to_family(raw: pd.DataFrame, family_map: dict[str, str]) -> pd.DataFrame:
    """Sum raw deconvolution columns onto the locked family names."""
    samples = raw["sample_id"].tolist()
    result = pd.DataFrame({"sample_id": samples})
    for fam in ALL_CELLS:
        cols = [c for c in raw.columns if c != "sample_id" and family_map.get(c, c) == fam]
        # also accept direct exact column match even if not in explicit map
        if fam in raw.columns and fam not in cols:
            cols.append(fam)
        if cols:
            result[fam] = raw[cols].sum(axis=1).values
        else:
            result[fam] = 0.0
    return result


def load_cluster_map(path: str | None) -> dict[str, str]:
    if path is None:
        return dict(DEFAULT_CLUSTER_FAMILY_MAP)
    p = Path(path)
    if not p.exists():
        print(f"[--] Cluster map {p} not found; using built-in default.")
        return dict(DEFAULT_CLUSTER_FAMILY_MAP)
    with open(p, "r", encoding="utf-8") as fh:
        j = json.load(fh)
    merged = dict(DEFAULT_CLUSTER_FAMILY_MAP)
    merged.update(j)
    return merged


def main() -> int:
    if not PRE_SPEC_PATH.exists():
        print("[!!] Rule: never run analysis before the test list is locked.")
        print(f"     Missing {PRE_SPEC_PATH}. Run scripts/write_sec4_prespec.py first.")
        return 2

    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--bulk", type=str, default=None)
    ap.add_argument("--signature", type=str, default=None,
                    help="Path to custom signature matrix CSV (rows=genes, cols=populations). "
                         "If omitted, falls back to snRNA direct proportions mode.")
    ap.add_argument("--snrna-props", type=str, default=None)
    ap.add_argument("--cluster-map", type=str, default=None,
                    help="Optional JSON mapping cluster-column -> family name.")
    ap.add_argument("--l2", type=float, default=0.0,
                    help="L2 ridge lambda for augmented-NNLS CIBERSORTx nu equivalent.")
    ap.add_argument("-o", "--output", type=str, default=None)
    args = ap.parse_args()

    out_path = Path(args.output) if args.output else DEFAULT_OUTPUT
    family_map = load_cluster_map(args.cluster_map)

    # --- Resolve bulk ---
    bulk_path = None
    if args.bulk:
        bulk_path = Path(args.bulk)
    elif DEFAULT_BULK_PQ.exists():
        bulk_path = DEFAULT_BULK_PQ
    elif DEFAULT_BULK_GZ.exists():
        bulk_path = DEFAULT_BULK_GZ

    sigmat_path = Path(args.signature) if args.signature else (
        DEFAULT_SIGMAT if DEFAULT_SIGMAT.exists() else None
    )

    # --- Mode B: snRNA direct proportions (before signature mat built or as cross-check) ---
    if sigmat_path is None:
        snrna_path = Path(args.snrna_props) if args.snrna_props else DEFAULT_SNRNA_PROP
        if not snrna_path.exists():
            print("[--] No signature matrix and no snRNA proportions found. Expected:")
            print(f"     Signature: {DEFAULT_SIGMAT}")
            print(f"     SnRNA props: {DEFAULT_SNRNA_PROP}")
            print(f"     Bulk: {DEFAULT_BULK_PQ}")
            print("[i ] This script becomes runnable AFTER:")
            print("     (a) acquire_bulk_discovery.py ran (bulk matrix exists), AND")
            print("     (b1) acquire_snrnaseq_discovery.py ran (snRNA proportions exist for Mode B), OR")
            print("     (b2) build_custom_signature_matrix.py ran (Mode A signature matrix exists).")
            empty = pd.DataFrame(columns=["sample_id"] + ALL_CELLS)
            out_path.parent.mkdir(parents=True, exist_ok=True)
            empty.to_csv(out_path, index=False)
            print(f"[OK] Placeholder written -> {out_path}")
            return 0
        raw = pd.read_csv(snrna_path)
        if "sample_id" not in raw.columns:
            raw = raw.rename(columns={raw.columns[0]: "sample_id"})
        fam = collapse_to_family(raw, family_map)
        fam.to_csv(out_path, index=False)
        print(f"[OK] Mode-B (direct snRNA proportions -> family). {len(fam)} samples.")
        print(f"[OK] Output -> {out_path}")
        return 0

    # --- Mode A: signature-matrix based NNLS ---
    if bulk_path is None or not bulk_path.exists():
        print(f"[!!] Signature matrix exists ({sigmat_path}) but bulk expression missing.")
        print(f"     Expected bulk path: {DEFAULT_BULK_PQ}")
        return 3

    bulk_df = load_expr(bulk_path)
    sigmat_df = pd.read_csv(sigmat_path, index_col=0)
    sigmat_df = sigmat_df.abs()  # NNLS requires non-negative; sign info in logFC direction handled upstream
    sigmat_df = sigmat_df.replace([np.inf, -np.inf], np.nan).fillna(0.0)

    raw = run_deconv(bulk_df, sigmat_df, l2_lambda=args.l2)
    fam = collapse_to_family(raw, family_map)
    fam.to_csv(out_path, index=False)
    print(f"[OK] Mode-A (custom CIBERSORTx NNLS{' L2='+str(args.l2) if args.l2>0 else ''}).")
    print(f"[OK] {fam.shape[0]} samples x {len(ALL_CELLS)} families -> {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
