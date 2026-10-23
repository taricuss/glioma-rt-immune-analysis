#!/usr/bin/env python3
"""
Phase 4.3 — CGGA cross-sectional support (IDH-mutant recurrent only).

Plan Sec 5.8 Step 4.3:
  If CGGA 'therapy' field is usable after CSV download, run cross-sectional
  treatment-vs-untreated fraction contrast at the RECURRENT TIMEPOINT ONLY.
  Label output CROSS-SECTIONAL SUPPORT (not longitudinal validation).

Pipeline:
  1. Load CGGA clinical → filter IDH-mutant AND Recurrent → classify Any-RT/No-RT
  2. Load CGGA RSEM-genes expression (TPM linear scale → gene rows, sample cols)
  3. Run 3-method deconvolution (quanTIseq / EPIC / MCP-counter) using the same
     embedded signature matrices as deconvolution_robustness.py (NNLS + family map)
  4. Cross-sectional stats per cell type (Recurrent samples only):
       - arcsin(sqrt(fraction)) response transform (same as discovery LMM)
       - Any_RT  n=106  vs  No_RT  n=23
       - Mann-Whitney U (non-parametric cross-sectional) + Welch t (parametric)
       - BH-FDR restricted to K=10 PRIMARY family only
       - Direction consistency with discovery Part-3 LMM interaction beta sign
  5. Write CSV artifacts + summary JSON to outputs/intermediate/.
"""
from __future__ import annotations

import io
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls as _nnls
from scipy.stats import mannwhitneyu, ttest_ind

if sys.stdout.encoding is None or "utf" not in sys.stdout.encoding.lower():
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PRE_SPEC_PATH = PROJECT_ROOT / "outputs" / "sec4_prespec.json"

DEFAULT_CLINICAL = PROJECT_ROOT / "data" / "cgga" / "CGGA.mRNAseq_693_clinical.20200506.txt"
DEFAULT_EXPR = PROJECT_ROOT / "data" / "cgga" / "CGGA.mRNAseq_693.RSEM-genes.20200506.txt"
DEFAULT_DISCOVERY_LMM = PROJECT_ROOT / "outputs" / "intermediate" / "lmm_interaction_summary.csv"

OUT_DIR = PROJECT_ROOT / "outputs" / "intermediate"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_CLINICAL_CANONICAL = OUT_DIR / "cgga_clinical_canonical.csv"
OUT_DECONV_QUANTISEQ = OUT_DIR / "cgga_deconv_quantiseq.csv"
OUT_DECONV_EPIC = OUT_DIR / "cgga_deconv_epic.csv"
OUT_DECONV_MCP = OUT_DIR / "cgga_deconv_mcpcounter.csv"
OUT_STATS = OUT_DIR / "cgga_crosssectional_stats.csv"
OUT_AUDIT = OUT_DIR / "cgga_pipeline_audit.json"

PRIMARY_CELLS = [
    "CD8_T","CD4_conv","Treg_FOXP3","NK","B","Plasma",
    "Microglia_P2RY12_TMEM119","MoMF_CD163","DC","Neutrophils"
]
EXPLORATORY_CELLS = [
    "T_exhausted_PDCD1_HAVCR2","M1","M2","Endothelial","Fibroblast_pericyte"
]
ALL_CELLS = PRIMARY_CELLS + EXPLORATORY_CELLS

# =============================================================================
# Signature matrices (copied from deconvolution_robustness.py — keep in sync)
# =============================================================================

QUANTISEQ_POPULATIONS = [
    "B","CD4_T_conv","CD8_T","DC","Monocyte_Macrophage","Neutrophils",
    "NK","Treg","T_CD4_general","Endothelial"
]
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

_QUANTI_GENE_ROWS: dict[str, dict[str, float]] = {
    "CD8A": {"CD8_T": 3.0}, "CD8B": {"CD8_T": 2.8}, "GZMA": {"CD8_T": 2.5},
    "GZMB": {"CD8_T": 2.4}, "PRF1": {"CD8_T": 2.2}, "NKG7": {"CD8_T": 2.0},
    "KLRD1": {"CD8_T": 1.8}, "CXCR3": {"CD8_T": 1.6},
    "CD4": {"CD4_T_conv": 2.8, "T_CD4_general": 2.6},
    "CD3D": {"CD4_T_conv": 2.5, "T_CD4_general": 2.3, "CD8_T": 2.2},
    "CD3E": {"CD4_T_conv": 2.5, "T_CD4_general": 2.3, "CD8_T": 2.2},
    "CD3G": {"CD4_T_conv": 2.4, "T_CD4_general": 2.2},
    "IL7R": {"CD4_T_conv": 2.0, "T_CD4_general": 1.9},
    "CCR7": {"CD4_T_conv": 1.8, "T_CD4_general": 1.7},
    "LEF1": {"CD4_T_conv": 1.6},
    "FOXP3": {"Treg": 3.5}, "IL2RA": {"Treg": 2.8, "T_CD4_general": 1.2},
    "CTLA4": {"Treg": 2.6}, "IKZF2": {"Treg": 2.2}, "TNFRSF18": {"Treg": 2.0},
    "NCAM1": {"NK": 2.6}, "NCR1": {"NK": 2.8}, "KLRB1": {"NK": 2.2},
    "KLRC1": {"NK": 2.0}, "KLRK1": {"NK": 2.1}, "GZMK": {"NK": 1.8},
    "CD19": {"B": 3.0}, "MS4A1": {"B": 2.9}, "CD79A": {"B": 2.8}, "CD79B": {"B": 2.7},
    "PAX5": {"B": 2.5}, "IGHM": {"B": 2.2}, "IGHD": {"B": 2.0},
    "CD163": {"Monocyte_Macrophage": 3.0}, "CD68": {"Monocyte_Macrophage": 2.6},
    "CSF1R": {"Monocyte_Macrophage": 2.5}, "MRC1": {"Monocyte_Macrophage": 2.4},
    "S100A8": {"Monocyte_Macrophage": 1.8, "Neutrophils": 2.0},
    "S100A9": {"Monocyte_Macrophage": 1.8, "Neutrophils": 2.0},
    "CCL2": {"Monocyte_Macrophage": 1.8},
    "CEACAM8": {"Neutrophils": 3.2}, "FCGR3B": {"Neutrophils": 3.0},
    "CXCR2": {"Neutrophils": 2.4}, "ELANE": {"Neutrophils": 2.4},
    "MPO": {"Neutrophils": 2.2},
    "CD1C": {"DC": 2.8}, "CLEC10A": {"DC": 2.6}, "FCER1A": {"DC": 2.4},
    "CLEC9A": {"DC": 2.2}, "BATF3": {"DC": 1.8}, "IRF8": {"DC": 1.6},
    "PECAM1": {"Endothelial": 2.8}, "VWF": {"Endothelial": 2.6},
    "CDH5": {"Endothelial": 2.5}, "ENG": {"Endothelial": 2.3},
    "ESAM": {"Endothelial": 2.0}, "KDR": {"Endothelial": 1.8},
}

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
    "CD8A": {"CD8_T_cells": 3.0}, "CD8B": {"CD8_T_cells": 2.8},
    "GZMA": {"CD8_T_cells": 2.4}, "GZMB": {"CD8_T_cells": 2.3},
    "NKG7": {"CD8_T_cells": 1.9}, "CST7": {"CD8_T_cells": 1.7},
    "CD4": {"CD4_T_cells": 2.9}, "ITK": {"CD4_T_cells": 2.1},
    "LTB": {"CD4_T_cells": 1.7}, "LDHB": {"CD4_T_cells": 1.5},
    "FOXP3": {"Tregs": 3.4}, "IL2RA": {"Tregs": 2.7, "CD4_T_cells": 1.1},
    "CTLA4": {"Tregs": 2.5}, "IKZF2": {"Tregs": 2.1},
    "NCAM1": {"NK_cells": 2.5}, "NCR1": {"NK_cells": 2.7},
    "KLRD1": {"NK_cells": 2.0}, "KLRK1": {"NK_cells": 2.0},
    "XCL2": {"NK_cells": 2.2}, "XCL1": {"NK_cells": 2.2},
    "CD19": {"B_cells": 2.9}, "MS4A1": {"B_cells": 2.8},
    "CD79A": {"B_cells": 2.7}, "BANK1": {"B_cells": 2.2},
    "PNOC": {"B_cells": 2.0}, "BLK": {"B_cells": 2.1},
    "CD68": {"Macrophages": 2.5}, "CD163": {"Macrophages": 2.9},
    "MRC1": {"Macrophages": 2.3}, "MSR1": {"Macrophages": 2.4},
    "CSF1R": {"Macrophages": 2.4}, "CCL3": {"Macrophages": 1.9},
    "PDGFRA": {"CAFs": 2.5}, "FAP": {"CAFs": 2.7}, "PDPN": {"CAFs": 2.3},
    "THY1": {"CAFs": 2.2}, "COL1A2": {"CAFs": 2.0}, "COL3A1": {"CAFs": 2.0},
    "DCN": {"CAFs": 1.9}, "LOXL2": {"CAFs": 1.8}, "SPARC": {"CAFs": 1.8},
    "PECAM1": {"Endothelial": 2.7}, "VWF": {"Endothelial": 2.5},
    "CDH5": {"Endothelial": 2.4}, "FLT1": {"Endothelial": 2.1},
    "KDR": {"Endothelial": 1.8}, "ROBO4": {"Endothelial": 2.0},
    "ESAM": {"Endothelial": 1.9}, "SOX18": {"Endothelial": 1.6},
}

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


# =============================================================================
# Helpers (copied & adapted from deconvolution_robustness.py)
# =============================================================================

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
    """NNLS deconvolution on TPM-scale values. CGGA input is RSEM TPM
    (linear, not log) — pass through directly; auto-detect log via max<25."""
    bulk_genes = set(bulk.index.astype(str))
    sig_genes = set(sigmat.index.astype(str))
    common = sorted(bulk_genes & sig_genes)
    if len(common) < 20:
        print(f"[!!] Low signature overlap ({len(common)} genes) — results may be noisy.")
    B = bulk.loc[common].values.astype(np.float64)
    S = sigmat.loc[common].values.astype(np.float64)
    if float(np.nanmax(B)) < 25:
        print("[i ] Detected log2-scale expression; converting exp2 for NNLS.")
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


def run_mcpcounter(bulk: pd.DataFrame) -> pd.DataFrame:
    samples = list(bulk.columns)
    rows = []
    for pop, genes in MCP_SETS.items():
        common = [g for g in genes if g in set(bulk.index.astype(str))]
        if len(common) < 4:
            rows.append(np.full(len(samples), np.nan))
            continue
        sub = bulk.loc[common].values.astype(np.float64)
        rows.append(np.nanmean(sub, axis=0))
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


# =============================================================================
# CGGA-specific helpers
# =============================================================================

def arcsin_sqrt(x: np.ndarray) -> np.ndarray:
    xc = np.clip(np.asarray(x, dtype=np.float64), 0.0, 1.0)
    return np.arcsin(np.sqrt(xc))


def load_cgga_clinical(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.strip() for c in df.columns]
    df = df.rename(columns={
        "Radio_status (treated=1;un-treated=0)": "radio_status",
        "Chemo_status (TMZ treated=1;un-treated=0)": "chemo_status",
        "IDH_mutation_status": "idh_status",
        "1p19q_codeletion_status": "codeletion_1p19q",
        "MGMTp_methylation_status": "mgmt_status",
        "PRS_type": "sample_type",
        "Histology": "histology",
        "Grade": "grade",
        "Age": "age",
        "Gender": "gender",
        "OS": "os_months",
        "Censor (alive=0; dead=1)": "censor",
    })
    def classify_treatment(row):
        r = row["radio_status"]
        c = row["chemo_status"]
        if pd.isna(r) or pd.isna(c):
            return np.nan
        if r == 1.0:
            return "Any_RT"
        elif r == 0.0:
            return "No_RT"
        return np.nan
    df["treatment_arm"] = df.apply(classify_treatment, axis=1)
    def remap_grade(g):
        if pd.isna(g):
            return g
        s = str(g).strip().lower()
        if "ii" in s and "iii" not in s and "iv" not in s:
            return "WHO II"
        if "iii" in s and "iv" not in s:
            return "WHO III"
        if "iv" in s:
            return "WHO IV"
        return g
    df["who_grade_initial"] = df["grade"].apply(remap_grade)
    return df


def filter_cgga_idh_recurrent(clin: pd.DataFrame) -> pd.DataFrame:
    m1 = clin["idh_status"] == "Mutant"
    m2 = clin["sample_type"] == "Recurrent"
    sub = clin[m1 & m2].copy()
    return sub


def load_cgga_expr(path: Path, keep_samples: list[str]) -> pd.DataFrame:
    """CGGA RSEM-genes TPM file: first column = Gene_Name, then 693 sample cols.
    Returns DataFrame: index = HGNC gene symbol, columns = CGGA sample IDs (filtered)."""
    df = pd.read_csv(path, sep="\t")
    first_col = df.columns[0]
    df = df.rename(columns={first_col: "Gene_Name"})
    df["Gene_Name"] = df["Gene_Name"].astype(str).str.strip()
    # Deduplicate gene symbols: keep row with max row sum (highest-expressed isoform)
    if df["Gene_Name"].duplicated().any():
        df["_row_sum"] = df.drop(columns=["Gene_Name"]).sum(axis=1)
        df = df.sort_values("_row_sum", ascending=False).drop_duplicates(subset=["Gene_Name"], keep="first")
        df = df.drop(columns=["_row_sum"])
    df = df.set_index("Gene_Name")
    if keep_samples:
        keep = [s for s in keep_samples if s in df.columns]
        dropped = len(keep_samples) - len(keep)
        if dropped:
            print(f"[--] {dropped}/{len(keep_samples)} requested sample IDs not found in expression columns.")
        df = df[keep]
    return df


def bh_fdr_primary_only(pvals: pd.Series, is_primary_mask: pd.Series) -> pd.Series:
    """Benjamini-Hochberg restricted to PRIMARY cell types only.
    Non-PRIMARY rows → NaN q-values (per pre-spec §4.5(c))."""
    out = pd.Series(np.nan, index=pvals.index, dtype=np.float64)
    mask = is_primary_mask & pvals.notna()
    if mask.sum() == 0:
        return out
    sub_p = pvals[mask].values.astype(np.float64)
    n = len(sub_p)
    order = np.argsort(sub_p)
    ranks = np.empty(n, dtype=np.int64)
    ranks[order] = np.arange(1, n + 1)
    q_sub = sub_p * n / ranks
    # Ensure monotonicity: cummin from largest-rank to smallest
    q_sub_ordered = q_sub[order]
    q_sub_ordered = np.minimum.accumulate(q_sub_ordered[::-1])[::-1]
    q_sub[order] = q_sub_ordered
    q_sub = np.clip(q_sub, 0.0, 1.0)
    out.loc[mask] = q_sub
    return out


# =============================================================================
# Main
# =============================================================================

def main() -> int:
    if not PRE_SPEC_PATH.exists():
        print("[!!] Rule: never run analysis before the test list is locked.")
        print(f"     Missing {PRE_SPEC_PATH}. Run scripts/write_sec4_prespec.py first.")
        return 2

    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--clinical", type=str, default=str(DEFAULT_CLINICAL))
    ap.add_argument("--expr", type=str, default=str(DEFAULT_EXPR))
    ap.add_argument("--discovery-lmm", type=str, default=str(DEFAULT_DISCOVERY_LMM),
                    help="Path to discovery lmm_interaction_summary.csv for sign consistency check.")
    args = ap.parse_args()

    clinical_path = Path(args.clinical)
    expr_path = Path(args.expr)
    discovery_lmm_path = Path(args.discovery_lmm)

    if not clinical_path.exists():
        print(f"[!!] Clinical file not found: {clinical_path}")
        return 1
    if not expr_path.exists():
        print(f"[!!] Expression file not found: {expr_path}")
        return 1

    # ---- (1) Clinical harmonization & filter ----
    clin_raw = load_cgga_clinical(clinical_path)
    n_all = len(clin_raw)
    print(f"[OK] CGGA clinical raw: {n_all} samples")

    clin_idh_rec = filter_cgga_idh_recurrent(clin_raw)
    n_idh_rec_total = len(clin_idh_rec)
    print(f"[OK] IDH-mutant × Recurrent: {n_idh_rec_total}")

    clin = clin_idh_rec.dropna(subset=["treatment_arm"]).copy()
    n_arm = len(clin)
    arm_counts = clin["treatment_arm"].value_counts().to_dict()
    print(f"[OK] IDH-mutant × Recurrent × treatment non-NA: n={n_arm}  {arm_counts}")

    clin = clin.rename(columns={"CGGA_ID": "sample_id"})
    clin.to_csv(OUT_CLINICAL_CANONICAL, index=False)
    print(f"[OK] Written canonical clinical -> {OUT_CLINICAL_CANONICAL}")

    # ---- (2) Expression ----
    sample_ids = clin["sample_id"].tolist()
    expr = load_cgga_expr(expr_path, keep_samples=sample_ids)
    print(f"[OK] Expression matrix: {expr.shape[0]} genes × {expr.shape[1]} samples")

    # ---- (3a) quanTIseq ----
    qs_sig = build_signature_matrix(_QUANTI_GENE_ROWS, QUANTISEQ_POPULATIONS)
    qs_raw = run_nnls_deconv(expr, qs_sig)
    qs_fam = collapse_to_family(qs_raw, QUANTISEQ_TO_FAMILY)
    qs_fam.to_csv(OUT_DECONV_QUANTISEQ, index=False)
    print(f"[OK] (1/3) quanTIseq NNLS -> {OUT_DECONV_QUANTISEQ}  (n={qs_fam.shape[0]})")

    # ---- (3b) EPIC ----
    ep_sig = build_signature_matrix(_EPIC_GENE_ROWS, EPIC_POPULATIONS)
    ep_raw = run_nnls_deconv(expr, ep_sig)
    ep_fam = collapse_to_family(ep_raw, EPIC_TO_FAMILY)
    ep_fam.to_csv(OUT_DECONV_EPIC, index=False)
    print(f"[OK] (2/3) EPIC NNLS -> {OUT_DECONV_EPIC}  (n={ep_fam.shape[0]})")

    # ---- (3c) MCP-counter ----
    mcp_raw = run_mcpcounter(expr)
    mcp_fam = mcp_collapse_to_family(mcp_raw)
    combined_mcp = mcp_raw.merge(mcp_fam, on="sample_id", how="inner")
    combined_mcp.to_csv(OUT_DECONV_MCP, index=False)
    print(f"[OK] (3/3) MCP-counter scores -> {OUT_DECONV_MCP}  (n={combined_mcp.shape[0]})")

    # ---- (4) Cross-sectional stats per cell type × per method ----
    discovery_signs: dict[str, float] = {}
    if discovery_lmm_path.exists():
        dlmm = pd.read_csv(discovery_lmm_path)
        if "beta_interaction" in dlmm.columns and "cell_family" in dlmm.columns:
            for _, r in dlmm.iterrows():
                cf = r["cell_family"]
                b = r["beta_interaction"]
                if pd.notna(b):
                    discovery_signs[cf] = np.sign(b)
        print(f"[OK] Loaded discovery sign map for sign-consistency: {len(discovery_signs)} cells")
    else:
        print(f"[i ] Discovery LMM summary not found at {discovery_lmm_path}; sign-consistency skipped.")

    method_frame_map: dict[str, pd.DataFrame] = {
        "quantiseq": qs_fam,
        "epic": ep_fam,
        "mcpcounter": mcp_fam,
    }

    anyrt_mask = (clin.set_index("sample_id")["treatment_arm"] == "Any_RT")
    nort_mask = (clin.set_index("sample_id")["treatment_arm"] == "No_RT")
    n_anyrt = int(anyrt_mask.sum())
    n_nort = int(nort_mask.sum())

    stats_rows = []
    for fam in ALL_CELLS:
        is_primary = fam in PRIMARY_CELLS
        for method_name, mdf in method_frame_map.items():
            mdf_idx = mdf.set_index("sample_id")
            if fam not in mdf_idx.columns:
                continue
            merged = mdf_idx[[fam]].join(anyrt_mask.rename("is_anyrt"), how="inner")
            merged = merged.dropna(subset=[fam, "is_anyrt"])
            if len(merged) < 10:
                continue
            anyrt_vals = merged.loc[merged["is_anyrt"], fam].values.astype(np.float64)
            nort_vals = merged.loc[~merged["is_anyrt"], fam].values.astype(np.float64)
            if len(anyrt_vals) < 5 or len(nort_vals) < 3:
                continue

            # For non-MCP methods: arcsin-sqrt before testing; for MCP (scores not fractions): use raw
            if method_name in ("quantiseq", "epic"):
                anyrt_x = arcsin_sqrt(anyrt_vals)
                nort_x = arcsin_sqrt(nort_vals)
            else:
                anyrt_x = anyrt_vals
                nort_x = nort_vals

            # Effect size: delta (Any_RT median - No_RT median), arcsin-scale when applicable
            med_anyrt = float(np.nanmedian(anyrt_x))
            med_nort = float(np.nanmedian(nort_x))
            delta_median = med_anyrt - med_nort

            # Mann-Whitney U (two-sided, default alternative)
            try:
                if len(set(np.concatenate([anyrt_x, nort_x]))) <= 1:
                    mw_p = np.nan
                    mw_u = np.nan
                else:
                    mw_u, mw_p = mannwhitneyu(anyrt_x, nort_x, alternative="two-sided", method="asymptotic")
            except Exception:
                mw_p = np.nan
                mw_u = np.nan

            # Welch t-test (two-sided)
            try:
                _, tt_p = ttest_ind(anyrt_x, nort_x, equal_var=False, nan_policy="omit")
            except Exception:
                tt_p = np.nan

            # Sign consistency with discovery
            cgga_sign = np.sign(delta_median) if delta_median != 0 else 0.0
            disc_sign = discovery_signs.get(fam, np.nan)
            sign_agree = bool(cgga_sign != 0 and not np.isnan(disc_sign) and cgga_sign == disc_sign)

            stats_rows.append({
                "cell_family": fam,
                "is_primary": is_primary,
                "method": method_name,
                "n_anyrt": n_anyrt,
                "n_nort": n_nort,
                "median_anyrt_transformed": med_anyrt,
                "median_nort_transformed": med_nort,
                "delta_median_anyrt_minus_nort": delta_median,
                "cgga_delta_sign": cgga_sign,
                "discovery_interaction_sign": disc_sign,
                "sign_consistent_with_discovery": sign_agree,
                "mannwhitney_u_stat": mw_u,
                "mannwhitney_p_2tail": mw_p,
                "welch_t_p_2tail": tt_p,
            })

    stats_df = pd.DataFrame(stats_rows)
    if not stats_df.empty:
        # BH-FDR restricted to PRIMARY family × per-method
        for method in stats_df["method"].unique():
            mask_m = stats_df["method"] == method
            pvals = stats_df.loc[mask_m, "mannwhitney_p_2tail"].reset_index(drop=True)
            prim_mask = stats_df.loc[mask_m, "is_primary"].reset_index(drop=True).astype(bool)
            qs = bh_fdr_primary_only(pvals, prim_mask)
            stats_df.loc[mask_m, "mannwhitney_q_primary_bh"] = qs.values

            pvals_t = stats_df.loc[mask_m, "welch_t_p_2tail"].reset_index(drop=True)
            qs_t = bh_fdr_primary_only(pvals_t, prim_mask)
            stats_df.loc[mask_m, "welch_t_q_primary_bh"] = qs_t.values

        stats_df.to_csv(OUT_STATS, index=False)
        print(f"\n[OK] Cross-sectional stats -> {OUT_STATS}  (n_rows={stats_df.shape[0]})")

        # Compact summary per cell family × method (p-value table)
        summary_pivot = stats_df.pivot_table(
            index="cell_family", columns="method",
            values=["mannwhitney_p_2tail", "sign_consistent_with_discovery"],
            aggfunc="first"
        )
        print("\n--- Mann-Whitney p (Any_RT vs No_RT, IDH-mut Recurrent only) ---")
        print(summary_pivot["mannwhitney_p_2tail"].round(4).reindex(ALL_CELLS).to_string())

        if "sign_consistent_with_discovery" in summary_pivot.columns.get_level_values(0):
            print("\n--- Sign consistency with discovery interaction beta (True=same sign) ---")
            print(summary_pivot["sign_consistent_with_discovery"].reindex(ALL_CELLS).to_string())
    else:
        pd.DataFrame(columns=[
            "cell_family","is_primary","method","n_anyrt","n_nort",
            "median_anyrt_transformed","median_nort_transformed",
            "delta_median_anyrt_minus_nort","cgga_delta_sign",
            "discovery_interaction_sign","sign_consistent_with_discovery",
            "mannwhitney_u_stat","mannwhitney_p_2tail","welch_t_p_2tail",
            "mannwhitney_q_primary_bh","welch_t_q_primary_bh"
        ]).to_csv(OUT_STATS, index=False)
        print("[--] No valid stats rows produced.")

    # ---- (5) Audit ----
    audit = {
        "cgga_version": "mRNAseq_693 (20200506)",
        "samples_total_clinical": n_all,
        "samples_idh_mut_recurrent_total": n_idh_rec_total,
        "samples_idh_mut_recurrent_treatment_nonna": n_arm,
        "treatment_arm_counts": arm_counts,
        "expression_shape_genes_samples": list(expr.shape),
        "n_methods": 3,
        "stats_rows": len(stats_df),
        "sign_consistency_overall": {},
        "output_paths": {
            "clinical_canonical": str(OUT_CLINICAL_CANONICAL),
            "deconv_quantiseq": str(OUT_DECONV_QUANTISEQ),
            "deconv_epic": str(OUT_DECONV_EPIC),
            "deconv_mcpcounter": str(OUT_DECONV_MCP),
            "stats": str(OUT_STATS),
        },
    }
    if not stats_df.empty:
        # Aggregate sign consistency: majority across methods per cell, then overall PRIMARY
        sc = stats_df.groupby("cell_family")["sign_consistent_with_discovery"].mean()
        audit["sign_consistency_per_cell"] = {k: float(v) for k, v in sc.items()}
        primary_cells = [c for c in PRIMARY_CELLS if c in sc.index]
        if primary_cells:
            audit["sign_consistency_overall"]["frac_primary_n_ge7_of_10_same_sign"] = bool(
                sum(float(sc.get(c, 0.0)) >= 0.5 for c in primary_cells) >= 7
            )
            audit["sign_consistency_overall"]["n_primary_10_with_majority_same_sign"] = int(
                sum(float(sc.get(c, 0.0)) >= 0.5 for c in primary_cells)
            )

    with open(OUT_AUDIT, "w", encoding="utf-8") as fh:
        json.dump(audit, fh, indent=2, default=str)
    print(f"\n[OK] Audit JSON -> {OUT_AUDIT}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
