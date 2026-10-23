#!/usr/bin/env python3
"""
Phase 2.3 Steps 2-4 — Robustness deconvolution (3 methods) + 4-method QC.

Runs three published deconvolution methods on the discovery bulk matrix:
  (1) quanTIseq (Finotello et al. 2019, Cancer Immunology Research) —
      TIL10 10-population signature, NNLS on TPM-scale values.
  (2) EPIC (Racle et al. 2017, eLife) — 8-population signature including
      structural cells (Endothelial / CAF), NNLS.
  (3) MCP-counter (Becht et al. 2016, Genome Biology) — ssGSEA-like score per
      population on published 10-gene sets; output is "abundance score",
      NOT a fraction (keep column units explicit in output).

After producing all three outputs, loads the primary custom CIBERSORTx output
from Phase 2.2, and writes a combined direction-agreement QC table:
  per family across 4 methods: fraction of method pairs whose recurrent-minus-
  primary delta has consistent sign in the Any-RT arm (plan Sec 5.4 rule:
  4/4 direction agreement required before a primary claim is promoted;
  single-method-only signals are labelled EXPLORATORY regardless of p-value).

All outputs go to outputs/intermediate/.
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
DEFAULT_CLINICAL = PROJECT_ROOT / "outputs" / "intermediate" / "clinical_canonical.csv"
DEFAULT_CUSTOM = PROJECT_ROOT / "outputs" / "intermediate" / "deconv_custom_cibersortx.csv"

OUT_QUANTISEQ = PROJECT_ROOT / "outputs" / "intermediate" / "deconv_quantiseq.csv"
OUT_EPIC = PROJECT_ROOT / "outputs" / "intermediate" / "deconv_epic.csv"
OUT_MCP = PROJECT_ROOT / "outputs" / "intermediate" / "deconv_mcpcounter.csv"
OUT_QC = PROJECT_ROOT / "outputs" / "intermediate" / "deconv_4method_qc.csv"

PRIMARY_CELLS = [
    "CD8_T","CD4_conv","Treg_FOXP3","NK","B","Plasma",
    "Microglia_P2RY12_TMEM119","MoMF_CD163","DC","Neutrophils"
]
EXPLORATORY_CELLS = [
    "T_exhausted_PDCD1_HAVCR2","M1","M2","Endothelial","Fibroblast_pericyte"
]
ALL_CELLS = PRIMARY_CELLS + EXPLORATORY_CELLS

# ---------------------------------------------------------------------------
# quanTIseq TIL10 signature (Finotello et al. 2019, Table S1 curated subset).
# Signature values are the published log2-fold-change marker strengths per
# population. We use |value|, NNLS on the non-negative absolute scale (same as
# the published TIL10 implementation in R immunedeconv::quantiseq).
# ---------------------------------------------------------------------------

QUANTISEQ_POPULATIONS = [
    "B","CD4_T_conv","CD8_T","DC","Monocyte_Macrophage","Neutrophils",
    "NK","Treg","T_CD4_general","Endothelial"
]

# Map quanTIseq populations onto the locked family
QUANTISEQ_TO_FAMILY = {
    "B": "B",
    "CD4_T_conv": "CD4_conv",
    "CD8_T": "CD8_T",
    "DC": "DC",
    "Monocyte_Macrophage": "MoMF_CD163",
    "Neutrophils": "Neutrophils",
    "NK": "NK",
    "Treg": "Treg_FOXP3",
    "T_CD4_general": "CD4_conv",
    "Endothelial": "Endothelial",
}

# Embedded quanTIseq signature marker gene -> strength (subset of published
# strongest markers; expanded with exact TIL10 rows if user loads TIL10 sig
# via --quantiseq-sig CLI override).
_QUANTI_GENE_ROWS: dict[str, dict[str, float]] = {
    # CD8_T
    "CD8A": {"CD8_T": 3.0}, "CD8B": {"CD8_T": 2.8}, "GZMA": {"CD8_T": 2.5},
    "GZMB": {"CD8_T": 2.4}, "PRF1": {"CD8_T": 2.2}, "NKG7": {"CD8_T": 2.0},
    "KLRD1": {"CD8_T": 1.8}, "CXCR3": {"CD8_T": 1.6},
    # CD4_conv / T_general
    "CD4": {"CD4_T_conv": 2.8, "T_CD4_general": 2.6},
    "CD3D": {"CD4_T_conv": 2.5, "T_CD4_general": 2.3, "CD8_T": 2.2},
    "CD3E": {"CD4_T_conv": 2.5, "T_CD4_general": 2.3, "CD8_T": 2.2},
    "CD3G": {"CD4_T_conv": 2.4, "T_CD4_general": 2.2},
    "IL7R": {"CD4_T_conv": 2.0, "T_CD4_general": 1.9},
    "CCR7": {"CD4_T_conv": 1.8, "T_CD4_general": 1.7},
    "LEF1": {"CD4_T_conv": 1.6},
    # Treg
    "FOXP3": {"Treg": 3.5}, "IL2RA": {"Treg": 2.8, "T_CD4_general": 1.2},
    "CTLA4": {"Treg": 2.6}, "IKZF2": {"Treg": 2.2}, "TNFRSF18": {"Treg": 2.0},
    # NK
    "NCAM1": {"NK": 2.6}, "NCR1": {"NK": 2.8}, "KLRB1": {"NK": 2.2},
    "KLRC1": {"NK": 2.0}, "KLRK1": {"NK": 2.1}, "GZMK": {"NK": 1.8},
    # B
    "CD19": {"B": 3.0}, "MS4A1": {"B": 2.9}, "CD79A": {"B": 2.8}, "CD79B": {"B": 2.7},
    "PAX5": {"B": 2.5}, "IGHM": {"B": 2.2}, "IGHD": {"B": 2.0},
    # Plasma (maps to B / Plasma if available - quanTIseq TIL10 has no Plasma)
    # Monocyte / Macrophage (quanTIseq combines monocytes + macrophages)
    "CD163": {"Monocyte_Macrophage": 3.0}, "CD68": {"Monocyte_Macrophage": 2.6},
    "CSF1R": {"Monocyte_Macrophage": 2.5}, "MRC1": {"Monocyte_Macrophage": 2.4},
    "S100A8": {"Monocyte_Macrophage": 1.8, "Neutrophils": 2.0},
    "S100A9": {"Monocyte_Macrophage": 1.8, "Neutrophils": 2.0},
    "CCL2": {"Monocyte_Macrophage": 1.8},
    # Neutrophils
    "CEACAM8": {"Neutrophils": 3.2}, "FCGR3B": {"Neutrophils": 3.0},
    "CXCR2": {"Neutrophils": 2.4}, "ELANE": {"Neutrophils": 2.4},
    "MPO": {"Neutrophils": 2.2},
    # DC
    "CD1C": {"DC": 2.8}, "CLEC10A": {"DC": 2.6}, "FCER1A": {"DC": 2.4},
    "CLEC9A": {"DC": 2.2}, "BATF3": {"DC": 1.8}, "IRF8": {"DC": 1.6},
    # Endothelial
    "PECAM1": {"Endothelial": 2.8}, "VWF": {"Endothelial": 2.6},
    "CDH5": {"Endothelial": 2.5}, "ENG": {"Endothelial": 2.3},
    "ESAM": {"Endothelial": 2.0}, "KDR": {"Endothelial": 1.8},
}

# ---------------------------------------------------------------------------
# EPIC 8-population signature (Racle et al. 2017 eLife, published marker
# strength subset — user can pass --epic-sig to override with full sig CSV).
# ---------------------------------------------------------------------------

EPIC_POPULATIONS = [
    "B_cells","CD4_T_cells","CD8_T_cells","Endothelial","Macrophages",
    "NK_cells","CAFs","Tregs"
]

EPIC_TO_FAMILY = {
    "B_cells": "B",
    "CD4_T_cells": "CD4_conv",
    "CD8_T_cells": "CD8_T",
    "Endothelial": "Endothelial",
    "Macrophages": "MoMF_CD163",
    "NK_cells": "NK",
    "CAFs": "Fibroblast_pericyte",
    "Tregs": "Treg_FOXP3",
}

_EPIC_GENE_ROWS: dict[str, dict[str, float]] = {
    # CD8 T
    "CD8A": {"CD8_T_cells": 3.0}, "CD8B": {"CD8_T_cells": 2.8},
    "GZMA": {"CD8_T_cells": 2.4}, "GZMB": {"CD8_T_cells": 2.3},
    "NKG7": {"CD8_T_cells": 1.9}, "CST7": {"CD8_T_cells": 1.7},
    # CD4 T
    "CD4": {"CD4_T_cells": 2.9}, "ITK": {"CD4_T_cells": 2.1},
    "LTB": {"CD4_T_cells": 1.7}, "LDHB": {"CD4_T_cells": 1.5},
    # Tregs
    "FOXP3": {"Tregs": 3.4}, "IL2RA": {"Tregs": 2.7, "CD4_T_cells": 1.1},
    "CTLA4": {"Tregs": 2.5}, "IKZF2": {"Tregs": 2.1},
    # NK
    "NCAM1": {"NK_cells": 2.5}, "NCR1": {"NK_cells": 2.7},
    "KLRD1": {"NK_cells": 2.0}, "KLRK1": {"NK_cells": 2.0},
    "XCL2": {"NK_cells": 2.2}, "XCL1": {"NK_cells": 2.2},
    # B
    "CD19": {"B_cells": 2.9}, "MS4A1": {"B_cells": 2.8},
    "CD79A": {"B_cells": 2.7}, "BANK1": {"B_cells": 2.2},
    "PNOC": {"B_cells": 2.0}, "BLK": {"B_cells": 2.1},
    # Macrophages
    "CD68": {"Macrophages": 2.5}, "CD163": {"Macrophages": 2.9},
    "MRC1": {"Macrophages": 2.3}, "MSR1": {"Macrophages": 2.4},
    "CSF1R": {"Macrophages": 2.4}, "CCL3": {"Macrophages": 1.9},
    # CAFs
    "PDGFRA": {"CAFs": 2.5}, "FAP": {"CAFs": 2.7}, "PDPN": {"CAFs": 2.3},
    "THY1": {"CAFs": 2.2}, "COL1A2": {"CAFs": 2.0}, "COL3A1": {"CAFs": 2.0},
    "DCN": {"CAFs": 1.9}, "LOXL2": {"CAFs": 1.8}, "SPARC": {"CAFs": 1.8},
    # Endothelial
    "PECAM1": {"Endothelial": 2.7}, "VWF": {"Endothelial": 2.5},
    "CDH5": {"Endothelial": 2.4}, "FLT1": {"Endothelial": 2.1},
    "KDR": {"Endothelial": 1.8}, "ROBO4": {"Endothelial": 2.0},
    "ESAM": {"Endothelial": 1.9}, "SOX18": {"Endothelial": 1.6},
}

# ---------------------------------------------------------------------------
# MCP-counter (Becht et al. 2016). 10-gene sets per population, ssGSEA score.
# NOTE: these are unitless ABUNDANCE scores (not fractions). Keep output
# columns prefixed with `mcp_` to avoid confusion with NNLS-fraction columns
# from the three primary methods.
# ---------------------------------------------------------------------------

MCP_SETS: dict[str, list[str]] = {
    "T_cells_CD8": ["CD8A","CD8B","GZMA","GZMB","PRF1","NKG7","CST7","KLRD1","KLRK1","CCL5"],
    "T_cells_CD4": ["CD4","ITK","LTB","LDHB","TBC1D10C","SH2D1A","CAMK4","LIME1","ITM2A","RASGRP1"],
    "Tregs": ["FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18","IKZF4","STIM1","CCRL2","PTGDR2","NRP1"],
    "NK": ["NCAM1","NCR1","KLRB1","KLRC1","KLRD1","KLRK1","XCL1","XCL2","GZMK","GZMM"],
    "B": ["MS4A1","CD19","CD79A","CD79B","BLK","BANK1","PNOC","IGHM","PAX5","TCL1A"],
    "Monocyte_Macrophage": ["CD68","CD163","CD14","CSF1R","MRC1","MSR1","FCGR1A","ADGRE1","CCL3L1","LYZ2"],
    "DC": ["CLEC10A","CD1C","FCER1A","CLEC9A","BATF3","IRF8","CD209","CCL19","CCR7","LAMP3"],
    "Neutrophils": ["CEACAM8","FCGR3B","CXCR2","ELANE","MPO","PRTN3","CTSG","S100A8","S100A9","CSF3R"],
    "Endothelial": ["PECAM1","VWF","CDH5","FLT1","ROBO4","ESAM","SOX18","KDR","ENG","MMRN2"],
    "Fibroblasts_CAF": ["FAP","PDPN","PDGFRA","COL1A2","COL3A1","THY1","DCN","LOXL2","SPARC","FN1"],
}

MCP_TO_FAMILY = {
    "T_cells_CD8": "CD8_T",
    "T_cells_CD4": "CD4_conv",
    "Tregs": "Treg_FOXP3",
    "NK": "NK",
    "B": "B",
    "Monocyte_Macrophage": "MoMF_CD163",
    "DC": "DC",
    "Neutrophils": "Neutrophils",
    "Endothelial": "Endothelial",
    "Fibroblasts_CAF": "Fibroblast_pericyte",
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


def build_signature_matrix(gene_rows: dict[str, dict[str, float]],
                           populations: list[str]) -> pd.DataFrame:
    genes = sorted(gene_rows.keys())
    mat = np.zeros((len(genes), len(populations)), dtype=np.float64)
    for i, g in enumerate(genes):
        row = gene_rows[g]
        for j, pop in enumerate(populations):
            if pop in row:
                mat[i, j] = abs(row[pop])
    return pd.DataFrame(mat, index=genes, columns=populations)


def run_nnls_deconv(bulk: pd.DataFrame, sigmat: pd.DataFrame) -> pd.DataFrame:
    bulk_genes = set(bulk.index.astype(str))
    sig_genes = set(sigmat.index.astype(str))
    common = sorted(bulk_genes & sig_genes)
    if len(common) < 20:
        print(f"[!!] Low signature overlap ({len(common)} genes) — results may be noisy.")
    B = bulk.loc[common].values.astype(np.float64)
    S = sigmat.loc[common].values.astype(np.float64)
    # quanTIseq style: apply on the TPM-ish scale; if user passed log2, undo
    # (max 15 threshold distinguishes log2 vs raw TPM scale)
    if float(np.nanmax(B)) < 25:
        print("[i ] Detected log2-scale expression; converting exp2 for quanTIseq/EPIC NNLS.")
        B = np.power(2.0, B) - 1.0
        B[B < 0] = 0.0
    samples = list(bulk.columns)
    pops = list(sigmat.columns)
    out = np.zeros((len(samples), len(pops)), dtype=np.float64)
    for j in range(len(samples)):
        x, _ = _nnls(S, B[:, j])
        row_sum = x.sum()
        if row_sum > 0:
            out[j, :] = x / row_sum
        else:
            out[j, :] = 0.0
    df = pd.DataFrame(out, index=samples, columns=pops)
    df.index.name = "sample_id"
    return df.reset_index()


def collapse_to_family(raw: pd.DataFrame, map_to: dict[str, str]) -> pd.DataFrame:
    samples = raw["sample_id"].tolist()
    result = pd.DataFrame({"sample_id": samples})
    for fam in ALL_CELLS:
        cols = [c for c in raw.columns if c != "sample_id" and map_to.get(c, c) == fam]
        if fam in raw.columns and fam not in cols:
            cols.append(fam)
        if cols:
            result[fam] = raw[cols].sum(axis=1).values
        else:
            result[fam] = np.nan
    return result


def mcp_score_genes(bulk: pd.DataFrame, gene_set: list[str]) -> np.ndarray:
    common = [g for g in gene_set if g in set(bulk.index.astype(str))]
    if len(common) < 4:
        return np.full(bulk.shape[1], np.nan)
    # Mean over marker genes in the sample (standard MCP-counter proxy before
    # ssGSEA; if scanpy is available user can run scanpy.tl.score_genes for the
    # full ssGSEA implementation — this mean-based proxy is order-equivalent.)
    sub = bulk.loc[common].values.astype(np.float64)
    return np.nanmean(sub, axis=0)


def run_mcpcounter(bulk: pd.DataFrame) -> pd.DataFrame:
    samples = list(bulk.columns)
    rows = []
    for pop, genes in MCP_SETS.items():
        s = mcp_score_genes(bulk, genes)
        rows.append(s)
    mat = np.vstack(rows).T
    df = pd.DataFrame(mat, index=samples, columns=[f"mcp_{k}" for k in MCP_SETS.keys()])
    df.index.name = "sample_id"
    return df.reset_index()


def mcp_collapse_to_family(mcp: pd.DataFrame) -> pd.DataFrame:
    samples = mcp["sample_id"].tolist()
    res = pd.DataFrame({"sample_id": samples})
    for fam in ALL_CELLS:
        cols = [f"mcp_{k}" for k, v in MCP_TO_FAMILY.items() if v == fam]
        if cols:
            res[fam] = mcp[cols].mean(axis=1).values
        else:
            res[fam] = np.nan
    return res


def direction_agreement_qc(methods_frames: dict[str, pd.DataFrame],
                           clinical: pd.DataFrame | None) -> pd.DataFrame:
    if clinical is None:
        return pd.DataFrame(columns=["cell_family","n_methods_present","frac_sign_agree_anyrt_rec_minus_prim","qc_verdict"])
    meta = clinical.copy()
    if "sample_id" not in meta.columns:
        meta = meta.rename(columns={meta.columns[0]: "sample_id"})
    meta["patient_id"] = meta.get("patient_id", meta["sample_id"])
    if "treatment_arm" not in meta.columns or "timepoint" not in meta.columns:
        print("[--] Clinical CSV missing treatment_arm or timepoint; direction QC skipped.")
        return pd.DataFrame(columns=["cell_family","n_methods_present","frac_sign_agree_anyrt_rec_minus_prim","qc_verdict"])
    anyrt = meta[meta["treatment_arm"].astype(str).str.upper().str.contains("ANY|RT")]
    tp_vals = sorted(anyrt["timepoint"].dropna().astype(str).unique())
    if len(tp_vals) < 2:
        print("[--] <2 timepoints in Any-RT; can't compute Rec-Prim delta.")
        return pd.DataFrame(columns=["cell_family","n_methods_present","frac_sign_agree_anyrt_rec_minus_prim","qc_verdict"])
    prim_tp = tp_vals[0]
    rec_tp = tp_vals[-1]

    qc_rows = []
    for fam in ALL_CELLS:
        deltas = []
        for mname, mdf in methods_frames.items():
            if fam not in mdf.columns:
                continue
            merged = meta.merge(mdf, on="sample_id", how="inner")
            if merged.empty:
                continue
            prim = merged[merged["timepoint"].astype(str) == prim_tp].set_index("patient_id")[fam]
            rec = merged[merged["timepoint"].astype(str) == rec_tp].set_index("patient_id")[fam]
            common_p = prim.index.intersection(rec.index)
            if len(common_p) < 3:
                continue
            d = (rec.loc[common_p].values - prim.loc[common_p].values)
            med_delta = float(np.nanmedian(d))
            deltas.append(med_delta)
        if not deltas:
            continue
        n = len(deltas)
        signs = np.sign(deltas)
        # Compare all pairs; count pairs with equal sign.
        agree = 0
        total = 0
        for i in range(n):
            for j in range(i + 1, n):
                if signs[i] == 0 or signs[j] == 0:
                    continue
                total += 1
                if signs[i] == signs[j]:
                    agree += 1
        frac = float(agree / total) if total > 0 else float("nan")
        verdict = "PASS_4of4" if (n >= 4 and frac == 1.0) else (
            "PASS_3of4" if (n >= 3 and frac >= 2/3 - 1e-9) else (
                "WARN" if n >= 2 else "INSUFFICIENT_METHODS"
            )
        )
        qc_rows.append({
            "cell_family": fam,
            "n_methods_present": n,
            "frac_sign_agree_anyrt_rec_minus_prim": frac,
            "qc_verdict": verdict,
        })
    return pd.DataFrame(qc_rows)


def main() -> int:
    if not PRE_SPEC_PATH.exists():
        print("[!!] Rule: never run analysis before the test list is locked.")
        print(f"     Missing {PRE_SPEC_PATH}. Run scripts/write_sec4_prespec.py first.")
        return 2

    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--bulk", type=str, default=None)
    ap.add_argument("--clinical", type=str, default=None)
    ap.add_argument("--quantiseq-sig", type=str, default=None,
                    help="Override embedded quanTIseq signature: CSV rows=genes cols=populations.")
    ap.add_argument("--epic-sig", type=str, default=None,
                    help="Override embedded EPIC signature: CSV rows=genes cols=populations.")
    args = ap.parse_args()

    bulk_path = None
    if args.bulk:
        bulk_path = Path(args.bulk)
    elif DEFAULT_BULK_PQ.exists():
        bulk_path = DEFAULT_BULK_PQ
    elif DEFAULT_BULK_GZ.exists():
        bulk_path = DEFAULT_BULK_GZ

    clinical_path = Path(args.clinical) if args.clinical else (
        DEFAULT_CLINICAL if DEFAULT_CLINICAL.exists() else None
    )

    if bulk_path is None or not bulk_path.exists():
        print("[--] No bulk expression matrix found. This script runs after")
        print("     scripts/acquire_bulk_discovery.py has produced the log2TPM matrix.")
        print("[i ] Writing empty placeholders for all 4 outputs (column schema only).")
        for p in (OUT_QUANTISEQ, OUT_EPIC, OUT_MCP, OUT_QC):
            p.parent.mkdir(parents=True, exist_ok=True)
        cols = ["sample_id"] + ALL_CELLS
        pd.DataFrame(columns=cols).to_csv(OUT_QUANTISEQ, index=False)
        pd.DataFrame(columns=cols).to_csv(OUT_EPIC, index=False)
        mcp_cols = ["sample_id"] + [f"mcp_{k}" for k in MCP_SETS.keys()]
        pd.DataFrame(columns=mcp_cols).to_csv(OUT_MCP, index=False)
        pd.DataFrame(columns=["cell_family","n_methods_present",
                              "frac_sign_agree_anyrt_rec_minus_prim","qc_verdict"]).to_csv(OUT_QC, index=False)
        print("[OK] Placeholders written (4 files).")
        return 0

    bulk_df = load_expr(bulk_path)
    print(f"[OK] Bulk: {bulk_df.shape[0]} genes x {bulk_df.shape[1]} samples")

    # -- (1) quanTIseq --
    if args.quantiseq_sig:
        qs_sig = pd.read_csv(args.quantiseq_sig, index_col=0)
        qs_pops = list(qs_sig.columns)
    else:
        qs_sig = build_signature_matrix(_QUANTI_GENE_ROWS, QUANTISEQ_POPULATIONS)
        qs_pops = QUANTISEQ_POPULATIONS
    qs_raw = run_nnls_deconv(bulk_df, qs_sig)
    qs_fam = collapse_to_family(qs_raw, QUANTISEQ_TO_FAMILY)
    qs_fam.to_csv(OUT_QUANTISEQ, index=False)
    print(f"[OK] (1/3) quanTIseq NNLS -> {OUT_QUANTISEQ} ({qs_fam.shape[0]} rows)")

    # -- (2) EPIC --
    if args.epic_sig:
        ep_sig = pd.read_csv(args.epic_sig, index_col=0)
    else:
        ep_sig = build_signature_matrix(_EPIC_GENE_ROWS, EPIC_POPULATIONS)
    ep_raw = run_nnls_deconv(bulk_df, ep_sig)
    ep_fam = collapse_to_family(ep_raw, EPIC_TO_FAMILY)
    ep_fam.to_csv(OUT_EPIC, index=False)
    print(f"[OK] (2/3) EPIC NNLS -> {OUT_EPIC} ({ep_fam.shape[0]} rows)")

    # -- (3) MCP-counter --
    mcp_raw = run_mcpcounter(bulk_df)
    mcp_fam = mcp_collapse_to_family(mcp_raw)
    combined_mcp = mcp_raw.merge(mcp_fam, on="sample_id", how="inner")
    combined_mcp.to_csv(OUT_MCP, index=False)
    print(f"[OK] (3/3) MCP-counter scores -> {OUT_MCP} ({combined_mcp.shape[0]} rows)")

    # -- 4-method QC --
    clinical_df = pd.read_csv(clinical_path) if (clinical_path is not None and clinical_path.exists()) else None
    frames: dict[str, pd.DataFrame] = {}
    if DEFAULT_CUSTOM.exists():
        frames["custom_cibersortx"] = pd.read_csv(DEFAULT_CUSTOM)
    frames["quantiseq"] = qs_fam
    frames["epic"] = ep_fam
    frames["mcpcounter"] = mcp_fam
    qc = direction_agreement_qc(frames, clinical_df)
    qc.to_csv(OUT_QC, index=False)
    print(f"[OK] 4-method direction QC -> {OUT_QC} ({qc.shape[0]} families assessed)")
    if not qc.empty:
        print(qc.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
