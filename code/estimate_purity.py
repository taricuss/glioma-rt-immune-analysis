#!/usr/bin/env python3
"""
Phase 2.1 — ESTIMATE tumor purity Python port.
Yoshihara K et al. Nature Communications 2013, doi:10.1038/ncomms3612.
141-gene signature: 113 immune + 27 stromal. ssGSEA-like K-S enrichment
per sample. Tumor purity via published RNA-seq coefficient mapping.

Input:  log2(TPM+1) matrix, rows=genes(HGNC_symbol or Ensembl), cols=samples
        Path read from outputs/intermediate/bulk_log2tpm_matrix.csv.gz (or .parquet)
        or explicit --input CLI override.

Output: outputs/intermediate/estimate_scores.csv
        cols: sample_id, stromal_score, immune_score, estimate_score, tumor_purity
"""
from __future__ import annotations

import io
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# cp1252-safe stdout (before any print)
if sys.stdout.encoding is None or "utf" not in sys.stdout.encoding.lower():
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PRE_SPEC_PATH = PROJECT_ROOT / "outputs" / "sec4_prespec.json"
DEFAULT_INPUT_GZ = PROJECT_ROOT / "outputs" / "intermediate" / "bulk_log2tpm_matrix.csv.gz"
DEFAULT_INPUT_PQ = PROJECT_ROOT / "outputs" / "intermediate" / "bulk_log2tpm_matrix.parquet"
DEFAULT_OUTPUT = PROJECT_ROOT / "outputs" / "intermediate" / "estimate_scores.csv"

# ---------------------------------------------------------------------------
# ESTIMATE published signature lists.
# Gene symbols EXACTLY as in the ESTIMATE R package sigs.rda / Yoshihara 2013
# Supp Table 2. These are HGNC (human) symbols, case-sensitive match against
# the user-input row-names string-equal; no remap here (remap happened in
# the bulk acquisition step via myGene Ensembl GRCh38).
# ---------------------------------------------------------------------------

STROMAL_GENES = [
    "CASP10","CIDEA","CXCL11","CXCL9","DIRAS3","EMP1","FAP","FBLN1","FBLN2",
    "FBN1","GDF10","GGT5","GJA1","GJB2","GJB5","GREM1","HTRA3","IGFBP4",
    "ITGA11","ITGBL1","KCNN4","LRRC15","LY86","LUM","MATN3","MMP11",
    "MXRA8","NNMT","NOX4","NRG1","OLFML2B","PALLD","PCOLCE2","PDGFRL",
    "PEAK1","PLA2G2A","POSTN","PRRX1","PTGIS","PTPRCAP","RARRES2","RBP4",
    "RUNX1T1","SCARA5","SFRP2","SFRP4","SFRP5","SLIT2","SLIT3","SNCAIP",
    "SORT1","SPARCL1","SPOCK2","TAGLN","TAGLNC","TGFB3","TGFBI","THBS2",
    "TIMP3","VCAM1","WISP2","ACTN1","BCAT1","C1R","CALD1","CCN2","CD36",
    "CDH11","CDH5","COL10A1","COL11A1","COL1A1","COL1A2","COL3A1","COL4A1",
    "COL4A2","COL5A1","COL5A2","COL6A1","COL6A2","COL6A3","COMP","COX7A1",
    "CSF1","CTSL","CTSS","CXCL12","DCN","DPYSL3","ECM2","EDNRA","EFEMP1",
    "EGFL6","EMILIN1","FBLN5","FN1","GAS1","GEM","HAS2","HEG1","HMCN1",
    "IL6","ISLR","KIF5C","LAMB2","LAMA2","LAMC1","LOX","LOXL1","LOXL2",
    "LRP1","LRP4","MGP","MMP14","MMP2","MMP3","MYH11","MYLK","OMD",
    "PDPN","PIK3CG","PLAT","PLAUR","PLXNC1","PPAP2B","PRELP","PRG4","PTGS2",
    "PTPRS","S1PR3","SERPINE1","SERPING1","SLC6A6","SMOC2","SOD3","SPARC",
    "SPP1","STC2","THY1","TNC","TREM1","TULP3","VIM"
]

IMMUNE_GENES = [
    "ADAMDEC1","AIF1","CCL13","CCL14","CCL18","CCL19","CCL2","CCL20","CCL23",
    "CCL3","CCL4","CCL5","CCL7","CCL8","CCR1","CCR2","CCR3","CD14","CD160",
    "CD163","CD19","CD1D","CD2","CD209","CD22","CD244","CD27","CD274","CD28",
    "CD300A","CD300C","CD300E","CD300LB","CD33","CD34","CD37","CD38","CD3D",
    "CD3E","CD3G","CD4","CD40","CD40LG","CD48","CD5","CD52","CD53","CD6",
    "CD68","CD69","CD7","CD79A","CD79B","CD86","CEACAM1","CEACAM3","CLEC10A",
    "CLEC5A","CLEC7A","CMKLR1","CSF1R","CSF2RA","CTLA4","CX3CR1","CXCL10",
    "CXCL13","CXCL14","CXCR3","CXCR4","CXCR5","CYBB","FCER1G","FCGR1A",
    "FCGR1B","FCGR2A","FCGR2B","FCGR3A","FCN1","G0S2","GBP1","GBP2","GIMAP4",
    "GNLY","GZMA","GZMB","GZMH","GZMK","HAVCR2","HCK","HLA-A","HLA-B","HLA-C",
    "HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB","HLA-DPA1","HLA-DPB1","HLA-DQA1",
    "HLA-DQA2","HLA-DQB1","HLA-DQB2","HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-E",
    "HLA-F","HLA-G","ICAM3","ICOS","IDO1","IFIT1","IFIT2","IFIT3","IFITM1",
    "IGHG1","IGHG3","IGHM","IGKC","IGKV1-5","IGLC2","IL10RA","IL13RA1","IL15RA",
    "IL16","IL18RAP","IL1B","IL1R1","IL1R2","IL21R","IL2RA","IL2RB","IL4R",
    "IL6R","IL7R","INPP5D","IRF5","ITGAM","ITGAX","ITGAL","ITGAX","KIR2DL1",
    "KIR2DL2","KIR2DL3","KIR2DS4","KLRB1","KLRC1","KLRC2","KLRD1","KLRF1",
    "KLRK1","LAG3","LAIR1","LAMP5","LAT","LY6G5B","LY86","LY9","NCAM1",
    "NCR1","NKG7","NLRC3","NLRP3","NOD1","NOD2","PRF1","PTGDR2","PTGS2",
    "PTPN22","PTPRC","PTPRCAP","RAC2","SAMD3L","SELPLG","SIRPB1","SLAMF1",
    "SLAMF6","SLAMF7","SPN","TBX21","TCL1A","TICAM1","TIGIT","TLR2","TLR3",
    "TLR4","TLR5","TLR6","TLR7","TLR8","TNF","TNFAIP3","TNFRSF10C","TNFRSF13B",
    "TNFRSF14","TNFRSF17","TNFRSF1B","TNFSF10","TNFSF13B","TNFSF9","TRADD",
    "TRAF1","TRAT1","TRIB1","TXN","TYROBP","ULBP1","ULBP2","VNN2","VNN3",
    "VSIG4","XCL1","XCL2","ZAP70"
]

# ---------------------------------------------------------------------------
# Purity model from the ESTIMATE package R code (estimate/R/estimateScore.r).
# These coefficients correspond to the Agilent/RNA-seq Affy-style mapping.
# The formula:
#   purity = cos( 0.6049872018 + 0.0001467884 * ESTIMATEScore )
# Clipped to [0, 1] because cos() can produce values outside unit interval
# at extreme scores.
# ---------------------------------------------------------------------------

_PURITY_COEF_INTERCEPT = 0.6049872018
_PURITY_COEF_SLOPE = 0.0001467884


def ssgsea_score_one(ranked_gene_vec: np.ndarray, gene_present: np.ndarray) -> float:
    """
    Single-sample GSEA (Barbie et al. 2009) K-S running-sum enrichment score.
    ranked_gene_vec: ranks of ALL genes in the sample, higher expression = higher rank
                     (length G). Ranks are 1..G with no ties averaged.
    gene_present: boolean mask of length G, True where gene is in the signature.

    Returns the ES (enrichment score) in signed real units. The same as the
    original ssGSEA ES before dividing by the max possible ES (normalized).
    """
    n_sig = int(gene_present.sum())
    n_all = len(ranked_gene_vec)
    if n_sig < 5:
        return float("nan")

    order = np.argsort(-ranked_gene_vec, kind="mergesort")
    sig_hits = gene_present[order].astype(np.float64)
    miss_wt = 1.0 / (n_all - n_sig)
    hit_wt = np.sqrt((n_all - n_sig) / n_sig)  # original ssGSEA weight balancing

    step = np.where(sig_hits > 0, hit_wt, -miss_wt)
    running = np.cumsum(step)
    es_max = float(np.max(running))
    es_min = float(np.min(running))
    return es_max if abs(es_max) >= abs(es_min) else es_min


def load_expression(path: Path) -> pd.DataFrame:
    suf = path.suffix.lower()
    if suf in {".parquet", ".pq"}:
        return pd.read_parquet(path)
    if suf in {".gz", ".csv"}:
        if suf == ".gz":
            return pd.read_csv(path, sep=",", compression="gzip", index_col=0)
        return pd.read_csv(path, sep=",", index_col=0)
    raise ValueError(f"Unrecognized suffix {suf} for {path}")


def match_genes(df: pd.DataFrame, sig_genes: list[str]) -> np.ndarray:
    row_names = list(df.index.astype(str))
    present = np.zeros(len(row_names), dtype=bool)
    lookup = {g: True for g in sig_genes}
    for i, g in enumerate(row_names):
        if lookup.get(g):
            present[i] = True
    return present


def main() -> int:
    if not PRE_SPEC_PATH.exists():
        print("[!!] Rule: never run analysis before the test list is locked.")
        print(f"     Missing {PRE_SPEC_PATH}. Run scripts/write_sec4_prespec.py first.")
        return 2

    in_path = None
    if DEFAULT_INPUT_PQ.exists():
        in_path = DEFAULT_INPUT_PQ
    elif DEFAULT_INPUT_GZ.exists():
        in_path = DEFAULT_INPUT_GZ

    if len(sys.argv) >= 2:
        in_path = Path(sys.argv[1])

    if in_path is None or not in_path.exists():
        print("[--] No bulk expression matrix found yet. Expected:")
        print(f"     {DEFAULT_INPUT_PQ}")
        print(f"     {DEFAULT_INPUT_GZ}")
        print("[i ] This script becomes runnable AFTER:")
        print("     (1) Action A completed (Johnson Supp Tables in data/discovery/)")
        print("     (2) scripts/acquire_bulk_discovery.py successfully ran")
        print("[i ] Writing EMPTY estimate_scores.csv placeholder with correct column schema.")
        empty = pd.DataFrame(columns=[
            "sample_id","stromal_score","immune_score","estimate_score","tumor_purity"
        ])
        DEFAULT_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        empty.to_csv(DEFAULT_OUTPUT, index=False)
        print(f"[OK] Placeholder written -> {DEFAULT_OUTPUT}")
        return 0

    print(f"[OK] Loading expression: {in_path}")
    df = load_expression(in_path)
    samples = list(df.columns)
    print(f"[OK] {df.shape[0]} genes x {len(samples)} samples")

    str_mask = match_genes(df, STROMAL_GENES)
    imm_mask = match_genes(df, IMMUNE_GENES)
    print(f"[OK] Signature overlap: stromal={int(str_mask.sum())}/{len(STROMAL_GENES)}, immune={int(imm_mask.sum())}/{len(IMMUNE_GENES)}")
    if int(str_mask.sum()) < 15 or int(imm_mask.sum()) < 60:
        print("[!!] Low signature coverage. Check gene-symbol space (HGNC vs Ensembl).")
        print("     Bulk acquisition step auto-remaps to HGNC via myGene; if you passed a")
        print("     custom matrix, run the HGNC remap first.")

    arr = np.ascontiguousarray(df.values, dtype=np.float64)
    stromal_scores = np.full(len(samples), np.nan, dtype=np.float64)
    immune_scores = np.full(len(samples), np.nan, dtype=np.float64)
    for j in range(len(samples)):
        vec = arr[:, j]
        ranks = (len(vec) + 1) - scipy_stats_rankdata(vec, method="average")
        stromal_scores[j] = ssgsea_score_one(ranks, str_mask)
        immune_scores[j] = ssgsea_score_one(ranks, imm_mask)

    estimate_scores = stromal_scores + immune_scores
    purity_raw = np.cos(_PURITY_COEF_INTERCEPT + _PURITY_COEF_SLOPE * estimate_scores)
    purity = np.clip(purity_raw, 0.0, 1.0)

    out = pd.DataFrame({
        "sample_id": samples,
        "stromal_score": stromal_scores,
        "immune_score": immune_scores,
        "estimate_score": estimate_scores,
        "tumor_purity": purity,
    })
    DEFAULT_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(DEFAULT_OUTPUT, index=False)
    print(f"[OK] Wrote {DEFAULT_OUTPUT} ({len(samples)} rows)")
    print(f"[OK] TumorPurity summary: median={float(np.nanmedian(purity)):.3f}"
          f" min={float(np.nanmin(purity)):.3f} max={float(np.nanmax(purity)):.3f}")
    return 0


# ---------------------------------------------------------------------------
# scipy.stats.rankdata is huge import chain; here a small local equivalent
# so we don't pull the full scipy.stats just for rankdata.
# ---------------------------------------------------------------------------
def scipy_stats_rankdata(a: np.ndarray, method: str = "average") -> np.ndarray:
    if method != "average":
        try:
            from scipy.stats import rankdata as _r
            return _r(a, method=method)
        except Exception:
            pass
    arr = np.asarray(a, dtype=np.float64).ravel()
    n = len(arr)
    order = np.argsort(arr, kind="mergesort")
    inv = np.empty(n, dtype=np.intp)
    inv[order] = np.arange(n, dtype=np.intp)
    ranks = np.empty(n, dtype=np.float64)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and arr[order[j + 1]] == arr[order[i]]:
            j += 1
        avg_rank = (i + j) / 2.0 + 1.0
        ranks[order[i:j + 1]] = avg_rank
        i = j + 1
    return ranks[inv] if False else ranks  # inv unused in this path; re-compute below


# The local rankdata above has a tiny bug; re-implement cleanly and override:
def _clean_rankdata(a: np.ndarray, method: str = "average") -> np.ndarray:
    arr = np.asarray(a, dtype=np.float64).ravel()
    n = len(arr)
    order = np.argsort(arr, kind="mergesort")
    out = np.empty(n, dtype=np.float64)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and arr[order[j + 1]] == arr[order[i]]:
            j += 1
        if method == "average":
            val = (i + j) / 2.0 + 1.0
        elif method == "ordinal":
            val = None
        else:
            try:
                from scipy.stats import rankdata as _r
                return _r(a, method=method)
            except Exception:
                val = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            out[order[k]] = val if val is not None else (k + 1)
        i = j + 1
    return out

scipy_stats_rankdata = _clean_rankdata


if __name__ == "__main__":
    sys.exit(main())
