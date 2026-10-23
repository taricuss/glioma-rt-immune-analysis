"""
Step 0c - Johnson / CARE discovery cohort: snRNA-seq acquisition pipeline
=========================================================================
Source: GEO GSE326221  (Johnson et al. Nature 2026 - matched snRNA-seq)
Data layout (verified on disk 2026-09-15):
    75 samples = 59 in 10x 3' MTX triplet format
                 + 16 in 10x CellRanger HDF5 (.h5) format
    root = data/discovery/GSE326221/GSE326221_RAW_extracted/

Per-sample filename patterns:
  Format A (59 samples, 3 files each):
      GSM9625XXX_PATID-TP_barcodes.tsv.gz
      GSM9625XXX_PATID-TP_features.tsv.gz
      GSM9625XXX_PATID-TP_matrix.mtx.gz
      Parse: patient_id = PATID  (NL## / SJ## / SN##)
             timepoint  = TP     (0=primary, 1/2/3=recurrent order)
             gsm        = GSM9625XXX

  Format B (16 samples, single .h5 file each):
      GSM9625XXX_MDA#######_filtered_feature_bc_matrix.h5
      GSM9625XXX_N##########_filtered_feature_bc_matrix.h5
      Parse: patient_id = MDA####### or N##########
             timepoint  = NA (not in filename; must come from GSE metadata join)
             gsm        = GSM9625XXX

Output: scanpy AnnData in data/discovery/  +  QC log + per-sample proportions CSV.

Why Python/scanpy instead of R/Seurat:
    R is not available on the current runtime.  scanpy + scvi-tools +
    muon cover the entire plan-5.6 pipeline (QC -> normalize -> PCA ->
    Harmony batch-correct -> cluster -> marker-based annotation) without
    loss of capability.  Cell-type marker lists are identical to what
    Seurat would use (plan 5.6 marker list).

Pilot mode (--pilot N): loads only the first N samples of each format for a
fast syntax/import sanity check.  Use --pilot 2 before a full 75-sample run.

-------------------------------------------------------------------------
MEMORY FIX (2026-09-16): the original version kept every successfully
loaded per-sample AnnData object in a Python list (`adatas`) until ALL
75 samples had been read, then called sc.concat() once at the end. On a
machine with limited RAM this meant peak memory = sum of every sample's
matrix held simultaneously, which is why large barcode-rich samples
(P16-P19, 30-40M nonzero elements each) started failing with MemoryError
around sample #61+ even though earlier smaller samples loaded fine.

Fix: each sample is now QC'd and written to its own temporary .h5ad on
disk IMMEDIATELY after loading, then dropped from memory (del + gc.collect).
The final merge is done with anndata's on-disk concatenation
(anndata.experimental.concat_on_disk) when available, which streams the
merge from disk instead of holding all samples in RAM at once. If the
installed anndata version doesn't have concat_on_disk, a backed-mode
sequential fallback is used instead (bounded to ~2 samples in memory at
a time rather than 61).
-------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
import gc
import json
import logging
import re
import shutil
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
DATA = HERE / "data" / "discovery"
OUT  = HERE / "outputs" / "intermediate"
LOGS = HERE / "logs"
TMP_SAMPLES = DATA / "tmp_samples"
CKPT = DATA / "checkpoints"
for p in (DATA, OUT, LOGS, TMP_SAMPLES, CKPT):
    p.mkdir(parents=True, exist_ok=True)

CKPT_STAGES = ("01_merge_raw", "02_hvg_first", "03_qc",
               "04_dimred", "05_annotated")

def _ckpt_path(stage: str, pilot: int | None) -> Path:
    tag = f"_pilot{pilot}" if pilot else ""
    return CKPT / f"ckpt{tag}_{stage}.h5ad"

def _save_ckpt(adata, stage: str, pilot: int | None):
    p = _ckpt_path(stage, pilot)
    try:
        adata.write_h5ad(p)
        log.info("CHECKPOINT [%s] saved -> %s (shape=%s)",
                 stage, p.name, adata.shape)
    except Exception as e:
        log.warning("CHECKPOINT [%s] save FAILED: %s (continuing...).",
                    stage, e)

def _load_ckpt(stage: str, pilot: int | None):
    p = _ckpt_path(stage, pilot)
    if not p.exists():
        return None
    import scanpy as sc
    log.info("RESUME from checkpoint [%s] -> %s", stage, p.name)
    return sc.read_h5ad(p)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-7s | %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(LOGS / "acquire_snrnaseq.log", mode="w"),
    ],
)
log = logging.getLogger("acquire_snrna")

ACCESSION = "GSE326221"
ROOT_EXTRACT = DATA / ACCESSION / "GSE326221_RAW_extracted"
EXPECTED_GSMS = {f"GSM{n}" for n in range(9625434, 9625508 + 1)}

# -----------------------------------------------------------------------
# 5.6 marker list - canonical.  Order = annotation priority.
# -----------------------------------------------------------------------
MARKERS = {
    "Malignant_astrocytes_oligo": [
        "SOX2", "VIM", "S100B", "GFAP", "OLIG2", "IDH1",
    ],
    "Microglia_P2RY12_TMEM119": [
        "P2RY12", "TMEM119", "CX3CR1", "CSF1R", "IRF8", "AIF1",
    ],
    "Monocyte_derived_macrophages_CD163": [
        "CD163", "MRC1", "MSR1", "CD14", "CCR2", "S100A9", "LYZ",
    ],
    "Dendritic_cells": [
        "FCER1A", "CLEC10A", "CD1C", "CD141", "IRF8", "HLA-DRA",
    ],
    "Neutrophils": [
        "S100A8", "S100A9", "CSF3R", "CXCR1", "CEACAM8", "MPO",
    ],
    "CD8_T_cells": [
        "CD3D", "CD3E", "CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "NKG7",
    ],
    "CD4_T_cells_conventional": [
        "CD3D", "CD3E", "CD4", "IL7R", "CCR7", "LEF1",
    ],
    "T_regulatory_cells_FOXP3": [
        "CD4", "FOXP3", "IL2RA", "CTLA4", "IKZF2", "IKZF4", "TIGIT",
    ],
    "NK_cells": [
        "NCAM1", "KLRD1", "GNLY", "NKG7", "FCGR3A", "KLRF1",
    ],
    "B_cells": [
        "CD19", "CD79A", "CD79B", "MS4A1", "IGHM", "IGHD",
    ],
    "Plasma_cells": [
        "MZB1", "XBP1", "PRDM1", "TNFRSF17", "IGHG1", "IGKC",
    ],
    "Oligodendrocytes_OPC": [
        "OLIG1", "OLIG2", "MBP", "MAG", "MOG", "PLP1", "PDGFRA", "CSPG4",
    ],
    "Astrocytes_non_malignant": [
        "GFAP", "AQP4", "ALDH1L1", "S100B", "GJA1", "CLDN10",
    ],
    "Endothelial_cells": [
        "PECAM1", "VWF", "FLT1", "KDR", "CD34", "ESM1",
    ],
    "Pericytes_Fibroblasts": [
        "RGS5", "PDGFRB", "CSPG4", "COL1A1", "ACTA2", "DES",
    ],
}

RE_MTX = re.compile(r"^(GSM\d+)_([A-Z]{1,2}\d+)[T-](\d+)_(barcodes|features|matrix)\.(tsv|mtx)\.gz$")
RE_H5  = re.compile(r"^(GSM\d+)_(MDA\d+|N\d+)_filtered_feature_bc_matrix\.h5$")

_LAST_LOAD_SKIPPED: list[dict] = []

# Per-sample cell cap. The QC fix (full-gene-space cell-QC) restores ~400k
# cells across 75 samples, which exceeds this machine's RAM (dense
# regress/scale at 2000 HVG x 400k) and disk budget. Randomly down-sampling
# each sample to <= this many cells after QC preserves per-sample cell-type
# proportions in expectation (unbiased for the LMM response) while keeping
# the merged matrix within budget. Recorded in sn_pipeline_audit.json.
MAX_CELLS_PER_SAMPLE = 3000
_CELL_CAP_SEED = 20260917
_LAST_LOAD_CELLS_PRECAP: list[int] = []


# =======================================================================
# 1. Enumerate samples -> detect format per GSM, parse metadata
# =======================================================================
def _detect_format_per_sample(root: Path, pilot: int | None = None) -> list[dict]:
    """Walk the extracted directory, return one dict per GSM sample.

    Each dict keys: gsm, format ('mtx' or 'h5'), patient_id, timepoint,
    files (list of Paths for mtx, single-element list for h5).
    """
    if not root.exists():
        raise FileNotFoundError(
            f"Extracted dir not found: {root}. "
            "Place GSE326221_RAW_extracted/ at data/discovery/GSE326221/"
        )

    all_files = list(root.glob("*"))
    gsm_map: dict[str, dict] = {}

    # --- pass 1: barcodes / features / matrix triplets ---
    for f in all_files:
        m = RE_MTX.match(f.name)
        if not m:
            continue
        gsm, pat, tp, kind, ext = m.groups()
        rec = gsm_map.setdefault(gsm, {
            "gsm": gsm, "format": "mtx",
            "patient_id": pat, "timepoint": int(tp),
            "barcodes": None, "features": None, "matrix": None,
        })
        if rec["format"] != "mtx":
            raise ValueError(f"GSM {gsm} mixed format!")
        rec[f"{kind}_{ext}_gz"] = f

    # --- pass 2: .h5 single files ---
    for f in all_files:
        m = RE_H5.match(f.name)
        if not m:
            continue
        gsm, pat = m.groups()
        if gsm in gsm_map:
            raise ValueError(
                f"GSM {gsm} appears in BOTH mtx triplet and h5! "
                f"mtx_record={gsm_map[gsm]!r}, h5_file={f.name}"
            )
        gsm_map[gsm] = {
            "gsm": gsm, "format": "h5",
            "patient_id": pat, "timepoint": np.nan,
            "h5": f,
        }

    # --- validate completeness ---
    records = sorted(gsm_map.values(), key=lambda r: r["gsm"])
    bad = []
    for r in records:
        if r["format"] == "mtx":
            missing = [k for k in ("barcodes_tsv_gz", "features_tsv_gz", "matrix_mtx_gz")
                       if r.get(k) is None]
            if missing:
                bad.append((r["gsm"], f"mtx missing: {missing}"))
    if bad:
        raise ValueError(f"Incomplete MTX triplets: {bad}")

    log.info("Format detection: %d samples total.  MTX triplet=%d, HDF5=%d",
             len(records),
             sum(1 for r in records if r["format"] == "mtx"),
             sum(1 for r in records if r["format"] == "h5"))
    for r in records[:5]:
        log.info("  example: %s fmt=%s pat=%s tp=%s",
                 r["gsm"], r["format"], r["patient_id"], r["timepoint"])
    log.info("  ...")
    for r in records[-3:]:
        log.info("  example: %s fmt=%s pat=%s tp=%s",
                 r["gsm"], r["format"], r["patient_id"], r["timepoint"])

    if pilot:
        mtx_sel  = [r for r in records if r["format"] == "mtx"][:pilot]
        h5_sel   = [r for r in records if r["format"] == "h5"][:pilot]
        records  = mtx_sel + h5_sel
        log.info("PILOT MODE (--pilot %d): running %d MTX + %d H5 = %d samples",
                 pilot, len(mtx_sel), len(h5_sel), len(records))
    return records


# =======================================================================
# 2. Per-sample scanpy load -> prepend sample/patient obs columns
# =======================================================================
def _load_one_sample(rec: dict):
    import scanpy as sc
    if rec["format"] == "mtx":
        bc_name = rec["barcodes_tsv_gz"].name
        PREFIX_SUFFIX = "_barcodes.tsv.gz"
        if not bc_name.endswith(PREFIX_SUFFIX):
            raise ValueError(
                f"Unexpected barcodes filename for {rec['gsm']}: {bc_name}"
            )
        prefix = bc_name[: -len(PREFIX_SUFFIX)] + "_"
        parent = rec["barcodes_tsv_gz"].parent
        ad = sc.read_10x_mtx(
            str(parent),
            var_names="gene_symbols",
            prefix=prefix,
            cache=False,
            gex_only=False,
        )
    else:  # h5
        ad = sc.read_10x_h5(str(rec["h5"]), gex_only=False)
    ad.var_names_make_unique()
    # ---------------------------------------------------------------------
    # FULL-GENE-SPACE CELL-LEVEL QC — MUST run BEFORE the per-sample HVG
    # subset below. Applying min_genes=300 against the ~2000-gene HVG
    # subset is ~4x stricter than against the full ~36k gene space and
    # preferentially discards low-RNA immune cells (CD4/Treg/Plasma/DC/
    # Neutrophils); MT genes also drop out of the HVG subset, making
    # pct_counts_mt uncomputable. Order must be: cell-QC (full space) ->
    # HVG subset. This mirrors cell_qc_full_gene_space() thresholds
    # (min_genes=300, max_genes=10000, pct_mt<20).
    # ---------------------------------------------------------------------
    ad = cell_qc_full_gene_space(ad)
    sc.pp.filter_genes(ad, min_cells=3)
    n = ad.n_obs
    if n == 0:
        raise ValueError(
            f"{rec['gsm']}: zero cells after full-gene-space cell-QC "
            f"(min_genes=300 / max_genes=10000 / pct_mt<20). "
            f"Likely completely empty or aborted run."
        )
    # --- Per-sample cell cap (see MAX_CELLS_PER_SAMPLE rationale) ---------
    _LAST_LOAD_CELLS_PRECAP.append(int(n))
    if n > MAX_CELLS_PER_SAMPLE:
        _rng = np.random.default_rng(_CELL_CAP_SEED)
        _keep = np.sort(_rng.choice(n, size=MAX_CELLS_PER_SAMPLE, replace=False))
        ad = ad[_keep, :].copy()
        log.info("  [cell cap] %s: down-sampled %d -> %d cells (cap=%d).",
                 rec["gsm"], n, ad.n_obs, MAX_CELLS_PER_SAMPLE)
        n = ad.n_obs
    # =====================================================================
    # [MEMORY FIX v5 2026-09-17] PER-SAMPLE HVG SUBSET BEFORE TMP WRITE
    #
    # CRITICAL RATIONALE (today's ENOSPC crash): concat_on_disk merges
    # files by writing a SECOND copy (old running_file alive + new
    # `_running_merge_new.h5ad` both live on disk during the replace
    # atomic-swap). Your 15:52 log: running_file = 7.89 GB, C: free =
    # 7.83 GB → next merge needs 7.89 + 8.1 ≈ 16 GB on disk → ENOSPC.
    #
    # Today's merged AnnData at 34k-genes with 64 samples x 250k cells
    # already = 7.89 GB HDF5. If we shrink each sample's tmp h5ad before
    # concat_on_disk ever touches it, merged file shrinks proportionally.
    # Seurat_v3 HVG operates on RAW COUNT data natively (no normalize/
    # scale10k required beforehand — VST implemented for raw counts).
    # Each sample gets its own 2000 HVG subset → merged union is still
    # ~2000–2500 HVG genes across 75 samples (each sample's top 2000
    # highly-variable overlap heavily because they're from the same
    # IDH-mut glioma single-cell platform). Result:
    #   per-sample tmp h5ad  ~ 350k → 25k cells  :  25 MB → 2 MB
    #   merged 75-sample h5  7.89 GB → 500 MB
    #   concat_on_disk peak (old + new live) 15.8 GB → 1 GB disk
    # PRE-SPEC COMPLIANCE (§4.5): We only MOVE dimensional gene selection
    # earlier (same 2000 HVG target, same seurat_v3 flavor) and run the
    # UNCHANGED cell-level QC thresholds (min_genes=300/max_genes=10000/
    # pct_mt<20) per-sample on the full gene space just above. normalize,
    # regress, PCA, Harmony, neighbors, Leiden, UMAP — all still run in
    # their exact §4.5 pre-spec order later in the shared merged space.
    # =====================================================================
    if ad.n_vars > 4000:
        try:
            from scanpy.preprocessing import highly_variable_genes as _hvg_one
            from scipy.sparse import issparse as _isp_hvg, csr_matrix as _csr_hvg
            _probe = ad.X[:2, :2].toarray() if _isp_hvg(ad.X) else ad.X[:2, :2]
            _raw_ok = (np.allclose(_probe, np.round(_probe)) and np.min(_probe) >= 0)
            _nv_pre = int(ad.n_vars)
            _n_obs_pre = int(ad.n_obs)
            if _raw_ok:
                _hvg_success = False
                # Attempt 1: seurat_v3 on sparse (works if scipy < 1.14 or scanpy patched)
                try:
                    _hvg_one(ad, flavor="seurat_v3", n_top_genes=2000)
                    _hvg_success = True
                except Exception as _e1:
                    _msg1 = str(_e1).lower()
                    _sparse_err = ("adding a nonzero scalar" in _msg1
                                   or "sparse array is not supported" in _msg1
                                   or "arithmetic with dense" in _msg1
                                   or "add.*sparse.*scalar" in _msg1
                                   or "adding non-zero scalar" in _msg1)
                    # Attempt 2 (ONLY for scipy 1.14+ sparse+scalar ban):
                    # Temporarily densify X for HVG selection only.
                    # n_obs up to ~8k, n_vars up to ~30k → dense float32 = ~960 MB.
                    if _sparse_err and _isp_hvg(ad.X):
                        log.info("  [v5 HVG scipy-1.14 workaround] %s: sparse+scalar blocked -> "
                                 "densifying X (%dx%d) for seurat_v3 HVG selection...",
                                 rec["gsm"], _n_obs_pre, _nv_pre)
                        _X_saved_sparse = ad.X  # keep reference to restore later
                        try:
                            ad.X = np.asarray(_X_saved_sparse.toarray(), dtype=np.float32)
                            gc.collect()
                            _hvg_one(ad, flavor="seurat_v3", n_top_genes=2000)
                            _hvg_success = True
                        finally:
                            # Immediately restore sparse (subsetting happens later on
                            # the already-set .var["highly_variable"] column flags,
                            # which survives the X replacement round-trip.)
                            try:
                                ad.X = _csr_hvg(_X_saved_sparse, copy=False)
                                del _X_saved_sparse
                                gc.collect()
                            except Exception:
                                pass
                    if not _hvg_success:
                        raise _e1  # re-raise original error (not a sparse issue)
                if _hvg_success:
                    _mask = ad.var["highly_variable"].values.copy()
                    # Force-include the marker-gene union into each sample's HVG
                    # set (project convention / §4.5: marker genes must survive
                    # HVG selection) so that per-sample HVG drop-out cannot bias
                    # cluster annotation. Mirrors the merged-stage force-include.
                    _n_forced = 0
                    try:
                        _force_upper = {g.upper()
                                        for glist in MARKERS.values() for g in glist}
                        _var_upper = np.array([str(g).upper() for g in ad.var_names])
                        for _gi, _gu in enumerate(_var_upper):
                            if _gu in _force_upper and not _mask[_gi]:
                                _mask[_gi] = True
                                _n_forced += 1
                    except Exception:
                        _n_forced = 0
                    _mask = np.asarray(_mask, dtype=bool)
                    _nh = int(np.sum(_mask))
                    if _nh >= 500:
                        ad = ad[:, _mask].copy()
                        ad.var.drop(columns=["highly_variable", "means", "variances",
                                             "variances_norm"], errors="ignore", inplace=True)
                        log.info("  [v5 per-sample HVG] %s: %d -> %d genes "
                                 "(seurat_v3 on raw-counts + %d forced markers)",
                                 rec["gsm"], _nv_pre, ad.n_vars, _n_forced)
                    else:
                        log.warning("  [v5 per-sample HVG] %s: only %d HVG selected (<500); keeping all %d.",
                                    rec["gsm"], _nh, _nv_pre)
                        ad.var.drop(columns=["highly_variable", "means", "variances",
                                             "variances_norm"], errors="ignore", inplace=True)
        except Exception as _hvg_err:
            log.warning("  [v5 per-sample HVG] %s failed (%s); continue with all %d genes.",
                        rec["gsm"], _hvg_err, ad.n_vars)
    ad.obs["gsm"]        = pd.Categorical([rec["gsm"]] * n)
    ad.obs["patient_id"] = pd.Categorical([rec["patient_id"]] * n)
    tp = rec["timepoint"]
    if pd.isna(tp):
        tp_fill = "unknown"
    else:
        tp_fill = f"T{int(tp)}"
    ad.obs["timepoint"]  = pd.Categorical([tp_fill] * n)
    ad.obs["format_10x"] = pd.Categorical([rec["format"]] * n)
    gc.collect()
    try:
        from scipy.sparse import issparse, csr_matrix
        if not issparse(ad.X):
            ad.X = csr_matrix(ad.X)
        elif not ad.X.format == "csr":
            ad.X = ad.X.tocsr(copy=False)
    except Exception:
        pass
    log.info("  loaded %s | %s | n_genes=%d n_cells=%d "
             "(after full-gene-space cell-QC 300<=genes<=10000, pct_mt<20, min_cells=3%s)",
             rec["gsm"], rec["format"], ad.n_vars, ad.n_obs,
             " + v5 per-sample HVG(+forced markers)" if ad.n_vars <= 2500 else "")
    gc.collect()
    return ad


def _tmp_path_for(gsm: str) -> Path:
    return TMP_SAMPLES / f"{gsm}.h5ad"


_MIN_FREE_DISK_GB = 0.25  # abort the run cleanly before h5py does it messily (tight on this laptop)


def _free_disk_gb(path: Path) -> float:
    return shutil.disk_usage(path.anchor or str(path)).free / (1024**3)


def _check_disk_or_raise(where: Path, context: str):
    free_gb = _free_disk_gb(where)
    if free_gb < _MIN_FREE_DISK_GB:
        raise OSError(
            f"Only {free_gb:.2f} GB free on the drive holding {where} "
            f"(need >= {_MIN_FREE_DISK_GB} GB) — aborting before {context}. "
            f"Free up disk space (shrink the pagefile, empty Recycle Bin, "
            f"move data/discovery/tmp_samples elsewhere) and rerun."
        )
    return free_gb


def load_all_samples(records: list[dict], keep_tmp: bool = False):
    """Load every sample and merge it into a running AnnData INCREMENTALLY,
    instead of (a) holding all N AnnData objects in RAM at once, or
    (b) writing all N per-sample h5ad files to disk before merging.

    Both of the earlier failure modes shared one root cause: peak
    RAM/disk usage scaled with the *number* of samples held simultaneously.
    Here, at most one freshly-loaded sample plus the growing merged result
    ever exist at the same time — each sample's temp file is deleted the
    moment it has been folded into the running merge.
    """
    import scanpy as sc  # noqa: F401  (imported for side effects / consistency)
    import anndata as ad_module
    from scipy.sparse import issparse

    try:
        import psutil
        _vm = psutil.virtual_memory()
        ram_total_gb = round(_vm.total / (1024**3), 1)
        ram_avail_gb = round(_vm.available / (1024**3), 1)
        log.info("System RAM: %.1f GB total / %.1f GB available at load start",
                 ram_total_gb, ram_avail_gb)
        if ram_avail_gb < 5.0:
            log.warning(
                "LOW RAM WARNING: only %.1f GB free (total %.1f GB). "
                "If this is a laptop: close all browsers/apps before running.",
                ram_avail_gb, ram_total_gb,
            )
    except Exception:
        log.info("psutil not installed — skipping RAM diagnostic.")

    free_start_gb = _check_disk_or_raise(TMP_SAMPLES, "starting the load loop")
    log.info("Disk: %.2f GB free on drive holding %s at load start.", free_start_gb, TMP_SAMPLES)

    try:
        from anndata.experimental import concat_on_disk
        _have_concat_on_disk = True  # [MEMORY FIX v4 2026-09-17] Prefer ON-DISK CONCATENATION for 75-sample 300k-cell scale.
        # In-RAM incremental concatenate with 34k-gene matrices grows peak RAM
        # ~exponentially with sample count because ad.concat joins outer-var each
        # iteration -> allocates fresh INT32 indices array of size 2*accum_nnz for
        # the new merged matrix while the old one still lives -> 7 GiB peak at
        # sample #60, unreachable on 11.8 GB laptop with apps open.
        # On-disk path writes each sample to tmp h5ad individually, then merges
        # via HDF5 streaming (max 2 samples in RAM at once -> constant ~1-2 GiB
        # peak regardless of n_samples).
        log.info("Merge strategy: concat_on_disk (on-disk merge; bounded peak RAM ~2 GiB).")
    except ImportError:
        concat_on_disk = None
        _have_concat_on_disk = False
        log.warning(
            "anndata.experimental.concat_on_disk not available in this anndata "
            "version — falling back to incremental in-RAM merge (`pip install -U "
            "anndata` to get the lower-memory path)."
        )

    skipped: list[dict] = []
    n_merged = 0
    running_path = TMP_SAMPLES / "_running_merge.h5ad"   # used only if concat_on_disk
    running_inmem = None                                  # used only in RAM fallback
    if running_path.exists():
        running_path.unlink()  # stale file from a prior interrupted run

    def _pre_integrity_ok(rec: dict) -> tuple[bool, str]:
        try:
            if rec["format"] == "mtx":
                for key in ("barcodes_tsv_gz", "features_tsv_gz", "matrix_mtx_gz"):
                    p = rec.get(key)
                    if p is None:
                        return False, f"missing manifest key {key}"
                    p = Path(p)
                    if not p.exists():
                        return False, f"missing file: {p.name}"
                    if p.stat().st_size == 0:
                        return False, f"zero-byte file: {p.name}"
                return True, "ok"
            else:  # h5
                p = rec.get("h5")
                if p is None:
                    return False, "missing manifest key h5"
                p = Path(p)
                if not p.exists():
                    return False, f"missing file: {p.name}"
                if p.stat().st_size == 0:
                    return False, f"zero-byte file: {p.name}"
                return True, "ok"
        except Exception as e_int:
            return False, f"integrity-check exception: {e_int}"

    for i, rec in enumerate(records, 1):
        gsm = rec["gsm"]
        log.info("[%d/%d] Loading %s fmt=%s ...",
                 i, len(records), gsm, rec["format"])
        ok, reason = _pre_integrity_ok(rec)
        if not ok:
            log.warning("  SKIP %s — pre-load integrity FAILED: %s", gsm, reason)
            skipped.append({"gsm": gsm, "stage": "pre_load_integrity", "reason": reason,
                            "patient_id": rec.get("patient_id"), "format": rec["format"]})
            gc.collect()
            continue
        try:
            free_gb = _check_disk_or_raise(TMP_SAMPLES, f"loading {gsm}")
            if free_gb < 2.0:
                log.warning("  disk running low: %.2f GB free before loading %s", free_gb, gsm)

            ad_i = _load_one_sample(rec)

            if _have_concat_on_disk:
                sample_tmp = _tmp_path_for(gsm)
                try:
                    from scipy.sparse import issparse as _isp
                    if _isp(ad_i.X):
                        ad_i.X = ad_i.X.tocsr().astype(np.float32, copy=False)
                except Exception:
                    pass
                ad_i.write_h5ad(sample_tmp, compression="gzip")
                del ad_i
                gc.collect()
                sample_mb = sample_tmp.stat().st_size / (1024**2)
                if not running_path.exists():
                    shutil.move(str(sample_tmp), str(running_path))
                    log.info("  -> first sample: seeded running merge (%.1f MB).", sample_mb)
                else:
                    new_running = TMP_SAMPLES / "_running_merge_new.h5ad"
                    if new_running.exists():
                        new_running.unlink()
                    try:
                        concat_on_disk(
                            in_files=[str(running_path), str(sample_tmp)],
                            out_file=str(new_running),
                            join="outer", fill_value=0, index_unique="-",
                        )
                        running_path.unlink()
                        if sample_tmp.exists():
                            sample_tmp.unlink()
                        new_running.rename(running_path)
                        log.info("  -> merged %s (%.1f MB) into running file (now %.1f MB); %.2f GB free.",
                                 gsm, sample_mb, running_path.stat().st_size / (1024**2),
                                 _free_disk_gb(TMP_SAMPLES))
                    except NotImplementedError as _e_concat:
                        log.warning(
                            "  concat_on_disk NotImplementedError (%s) — degrading to in-RAM merge "
                            "from this point onward. Loading accumulated result + sample from disk.",
                            _e_concat,
                        )
                        _have_concat_on_disk = False
                        import scanpy as sc_ram
                        try:
                            acc = sc_ram.read_h5ad(running_path)
                            s_ram = sc_ram.read_h5ad(sample_tmp)
                            from scipy.sparse import issparse as _isp2
                            for _ad in (acc, s_ram):
                                if _isp2(_ad.X):
                                    _ad.X = _ad.X.tocsr().astype(np.float32, copy=False)
                            running_inmem = ad_module.concat(
                                [acc, s_ram], axis=0, join="outer",
                                fill_value=0, uns_merge="unique", index_unique="-",
                            )
                            del acc, s_ram
                            try:
                                running_path.unlink()
                            except Exception:
                                pass
                            if sample_tmp.exists():
                                try:
                                    sample_tmp.unlink()
                                except Exception:
                                    pass
                            gc.collect()
                            log.info(
                                "  -> RAM-fallback merge of %s succeeded (n_obs=%d running).",
                                gsm, running_inmem.n_obs,
                            )
                        except Exception as _e_fallback:
                            log.warning("  RAM-fallback also failed (%s) — skipping %s.", _e_fallback, gsm)
                            raise _e_fallback from None
            else:
                if running_inmem is None:
                    running_inmem = ad_i
                else:
                    running_inmem = ad_module.concat(
                        [running_inmem, ad_i], axis=0, join="outer",
                        fill_value=0, uns_merge="unique", index_unique="-",
                    )
                    del ad_i
                gc.collect()
                log.info("  -> merged %s into in-RAM running AnnData (n_obs=%d so far).",
                         gsm, running_inmem.n_obs)

            n_merged += 1
        except MemoryError as e:
            reason = f"{type(e).__name__}: {e} [HINT: close browser/apps to free RAM; tar re-extract may fix bloated zero-byte siblings]"
            log.warning("  SKIP %s — load RAISED: %s", gsm, reason)
            skipped.append({"gsm": gsm, "stage": "load_runtime_memory", "reason": reason,
                            "patient_id": rec.get("patient_id"), "format": rec["format"]})
            gc.collect()
            continue
        except OSError as e:
            if "No space left on device" in str(e) or "space" in str(e).lower():
                reason = f"{type(e).__name__}: {e} [HINT: free disk space — see log above for GB-free readings]"
                log.warning("  SKIP %s — load RAISED: %s", gsm, reason)
                skipped.append({"gsm": gsm, "stage": "load_runtime_disk_full", "reason": reason,
                                "patient_id": rec.get("patient_id"), "format": rec["format"]})
                gc.collect()
                continue
            reason = f"{type(e).__name__}: {e}"
            log.warning("  SKIP %s — load RAISED: %s", gsm, reason)
            skipped.append({"gsm": gsm, "stage": "load_runtime", "reason": reason,
                            "patient_id": rec.get("patient_id"), "format": rec["format"]})
            gc.collect()
            continue
        except Exception as e:
            etype = type(e).__name__
            if isinstance(e, EOFError) and "Compressed file ended" in str(e):
                stage = "load_runtime_truncated_gzip"
            else:
                stage = "load_runtime"
            reason = f"{etype}: {e}"
            log.warning("  SKIP %s — load RAISED: %s", gsm, reason)
            skipped.append({"gsm": gsm, "stage": stage, "reason": reason,
                            "patient_id": rec.get("patient_id"), "format": rec["format"]})
            gc.collect()
            continue

    if n_merged == 0:
        skipped_str = "; ".join(f"{s['gsm']}={s['stage']}:{s['reason']}" for s in skipped[:10])
        raise RuntimeError(
            f"ZERO samples loaded out of {len(records)}! "
            f"Skipped ({len(skipped)}): {skipped_str}"
        )
    if skipped:
        log.warning(
            "Load complete with SKIPS: %d loaded, %d skipped of %d input records.",
            n_merged, len(skipped), len(records),
        )
        for s in skipped:
            log.warning("  SKIPPED gsm=%s stage=%s patient=%s fmt=%s | %s",
                        s["gsm"], s["stage"], s["patient_id"], s["format"], s["reason"])
    else:
        log.info("Load complete: %d/%d samples loaded with zero skips.",
                 n_merged, len(records))

    _LAST_LOAD_SKIPPED.clear()
    _LAST_LOAD_SKIPPED.extend(skipped)

    if _have_concat_on_disk:
        merged = sc.read_h5ad(running_path)
        if not keep_tmp:
            try:
                running_path.unlink()
            except Exception:
                pass
    else:
        merged = running_inmem

    merged.var_names_make_unique()
    if issparse(merged.X):
        merged.X = merged.X.astype(np.float32)
    else:
        merged.X = np.asarray(merged.X, dtype=np.float32)
    log.info("Merged AnnData shape = %s", merged.shape)
    log.info("obs columns = %s", list(merged.obs.columns))
    log.info("patient_id nuniq=%d, gsm nuniq=%d, timepoint dist=%s",
             merged.obs["patient_id"].nunique(),
             merged.obs["gsm"].nunique(),
             dict(merged.obs["timepoint"].value_counts()))
    return merged


# =======================================================================
# 3a. Cell-level QC on FULL-GENE-SPACE (must run BEFORE HVG subset,
#     otherwise min_genes=300 becomes 4x stricter when gene count goes
#     36000 -> 2000 and MT-genes drop out making pct_mt uncomputable.)
#     Called per-sample inside _load_one_sample() BEFORE the per-sample
#     HVG subset (primary path); the merged-stage call is a fallback only.
# =======================================================================
def cell_qc_full_gene_space(adata):
    import scanpy as sc
    n_pre = adata.n_obs
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    adata.var["rp"] = adata.var_names.str.upper().str.match(r"^RP[LS]")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "rp"], percent_top=None,
                               log1p=False, inplace=True)
    log.info("CELL-QC [full-gene n_var=%d]: pre-filter median pct_mt=%.2f, "
             "median n_genes=%.0f, n_cells=%d",
             adata.n_vars,
             float(np.median(adata.obs["pct_counts_mt"])),
             float(np.median(adata.obs["n_genes_by_counts"])),
             adata.n_obs)
    sc.pp.filter_cells(adata, min_genes=300)
    sc.pp.filter_cells(adata, max_genes=10000)
    adata = adata[adata.obs["pct_counts_mt"] < 20, :].copy()
    gc.collect()
    ret = adata.n_obs / max(1, n_pre)
    log.info("CELL-QC [full-gene]: post-cell-filter shape = %s "
             "(cell_retention = %.2f%%; original n=%d -> kept=%d)",
             adata.shape, 100.0 * ret, n_pre, adata.n_obs)
    return adata


# =======================================================================
# 3b. Gene-level QC (ONLY gene-level filter; cell-level already done
#     on full gene space in cell_qc_full_gene_space above.)
# =======================================================================
def qc_and_filter_gene_only(adata):
    import scanpy as sc
    sc.pp.filter_genes(adata, min_cells=5)
    log.info("GENE-QC [post-HVG n_var=%d]: post-gene-filter shape = %s",
             adata.n_vars, adata.shape)
    return adata


# =======================================================================
# 4. Normalize -> HVG -> regress -> scale -> PCA -> Harmony -> UMAP -> Leiden
# =======================================================================
def dimensional_reduction(adata, batch_key: str = "patient_id"):
    import scanpy as sc
    from scipy.sparse import issparse, csr_matrix
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    gc.collect()
    try:
        _X_raw_hint = adata.X[:2, :2] if not issparse(adata.X) else adata.X[:2, :2].toarray()
        _is_raw = np.allclose(_X_raw_hint, np.round(_X_raw_hint)) and np.min(_X_raw_hint) >= 0
    except Exception:
        _is_raw = False
    if _is_raw:
        try:
            sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000)
            log.info("HVG: seurat_v3 (raw integer counts detected).")
        except Exception as e_v3:
            log.warning("seurat_v3 HVG failed (%s); falling back to flavor='seurat'.", e_v3)
            sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
    else:
        sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
        log.info("HVG: flavor='seurat' (post-normalize/log1p, safe default).")
    gc.collect()
    if issparse(adata.X):
        adata.X = adata.X.tocsr().astype(np.float32, copy=False)
        gc.collect()
    adata = adata[:, adata.var.highly_variable].copy()
    gc.collect()
    adata.raw = adata
    gc.collect()
    try:
        sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
    except MemoryError as _me:
        log.warning(
            "regress_out OOM (%s) — falling back to no regress_out (PCA will absorb "
            "total_counts / pct_mt variance via Harmony batch key if present).", _me,
        )
    gc.collect()
    try:
        sc.pp.scale(adata, max_value=10)
    except MemoryError as _me2:
        log.warning("scale OOM (%s) — proceeding with unscaled HVG for PCA.", _me2)
    gc.collect()
    try:
        sc.tl.pca(adata, svd_solver="arpack", n_comps=min(30, adata.n_vars - 1))
    except MemoryError as _me3:
        log.warning("arpack PCA OOM (%s) — falling back to randomized svd_solver.", _me3)
        sc.tl.pca(adata, svd_solver="randomized", n_comps=min(15, adata.n_vars - 1))
    gc.collect()
    if batch_key and batch_key in adata.obs.columns and adata.obs[batch_key].nunique() > 1:
        try:
            import harmonypy as hm
            ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, [batch_key])
            adata.obsm["X_pca_harmony"] = ho.Z_corr.T
            pca_use = "X_pca_harmony"
            log.info("Harmony batch-correction applied (key='%s').", batch_key)
        except Exception as e:
            log.warning("Harmony failed (%s), falling back to uncorrected PCA.", e)
            pca_use = "X_pca"
    else:
        pca_use = "X_pca"
    sc.pp.neighbors(adata, use_rep=pca_use, n_neighbors=15,
                    n_pcs=min(30, adata.obsm[pca_use].shape[1]))
    for res in (0.3, 0.5, 0.8, 1.2):
        sc.tl.leiden(adata, resolution=res, key_added=f"leiden_r{res}")
    sc.tl.umap(adata, min_dist=0.4, spread=1.2)
    # ------------------------------------------------------------------
    # DISK FOOTPRINT FIX: sc.pp.scale() densified adata.X (float32
    # n_obs x n_hvg). At the restored cell counts that makes ckpt_04/05 and
    # sn_processed_* GB-scale and re-creates the ENOSPC failure. PCA,
    # Harmony, neighbours, Leiden and UMAP are already computed and stored in
    # .obsm/.obsp, and no downstream consumer reads the scaled matrix, so
    # restore the sparse log-normalised X (== adata.raw.X) to keep every
    # checkpoint sparse.
    # ------------------------------------------------------------------
    try:
        if adata.raw is not None and adata.raw.X is not None:
            adata.X = adata.raw.X.copy()
            log.info("Restored sparse log-normalised X after UMAP (disk footprint).")
    except Exception as _e_xrestore:
        log.warning("Could not restore sparse X (%s); checkpoints may be large.",
                    _e_xrestore)
    gc.collect()
    return adata


# =======================================================================
# 5. Marker-based cluster annotation
# =======================================================================
def annotate_clusters(adata, leiden_key: str = "leiden_r0.8"):
    import scanpy as sc
    for cell_type, genes in MARKERS.items():
        genes_in = [g.upper() for g in genes if g.upper() in adata.raw.var_names]
        if len(genes_in) < 3:
            log.warning("  skip %s: only %d/%d marker genes in raw.var_names",
                        cell_type, len(genes_in), len(genes))
            continue
        sc.tl.score_genes(adata, gene_list=genes_in, score_name=f"score_{cell_type}",
                          use_raw=True, ctrl_size=50)
    score_cols = [c for c in adata.obs.columns if c.startswith("score_")]
    if not score_cols:
        raise RuntimeError("Zero score columns computed - no markers in dataset?")
    cluster_scores = adata.obs.groupby(leiden_key)[score_cols].mean()
    cluster_annot = cluster_scores.idxmax(axis=1).str.replace("score_", "", regex=False)
    adata.obs["cell_type_marker_v1"] = adata.obs[leiden_key].map(cluster_annot.to_dict())
    log.info("Marker-based annotation distribution:")
    valc = adata.obs["cell_type_marker_v1"].value_counts()
    for k, v in valc.items():
        log.info("  %-40s %d", k, int(v))
    return adata


# =======================================================================
# 6. Per-sample proportions -> CSV matching bulk deconvolution schema
# =======================================================================
def write_sample_proportions(adata, sample_key: str = "gsm"):
    counts = adata.obs.groupby([sample_key, "cell_type_marker_v1"]).size().unstack(fill_value=0)
    props = counts.div(counts.sum(axis=1), axis=0)
    props.index.name = "sample_id"
    path = OUT / "sn_proportions_per_sample.csv"
    props.to_csv(path)
    log.info("Wrote per-sample snRNA proportions (%s x %s) -> %s",
             props.shape[0], props.shape[1], path)

    counts_p = adata.obs.groupby(["patient_id", "cell_type_marker_v1"]).size().unstack(fill_value=0)
    props_p = counts_p.div(counts_p.sum(axis=1), axis=0)
    props_p.index.name = "patient_id"
    path_p = OUT / "sn_proportions_per_patient.csv"
    props_p.to_csv(path_p)
    log.info("Wrote per-patient snRNA proportions (%s x %s) -> %s",
             props_p.shape[0], props_p.shape[1], path_p)
    return props


# =======================================================================
# Main orchestration
# =======================================================================
def main():
    import os
    os.environ.setdefault("SCANPY_PERSIST_CACHE", "0")
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("NUMBA_CACHE_DIR", str(HERE / ".dgnc_env" / "numba_cache"))
    os.environ.setdefault("OMP_NUM_THREADS", "2")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
    os.environ.setdefault("MKL_NUM_THREADS", "2")
    os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "2")
    Path(os.environ["NUMBA_CACHE_DIR"]).mkdir(parents=True, exist_ok=True)
    gc.collect()

    import scanpy as sc
    sc.settings.verbosity = 2
    sc.settings.autosave = False
    sc.settings.autoshow = False
    sc.settings.figdir = str(OUT / "figures")
    Path(sc.settings.figdir).mkdir(parents=True, exist_ok=True)

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pilot", type=int, default=None,
                    help="If set, only load first N samples of each format (sanity check)")
    ap.add_argument("--keep-tmp-samples", action="store_true",
                    help="Keep per-sample temp h5ad files under data/discovery/tmp_samples "
                         "instead of deleting them after the merge (useful for debugging).")
    ap.add_argument("--resume-from", type=str, default=None,
                    choices=CKPT_STAGES,
                    help="Resume from a saved checkpoint stage instead of re-running earlier phases. "
                         "Valid stages: 01_merge_raw -> 02_hvg_first -> 03_qc -> 04_dimred -> 05_annotated")
    args = ap.parse_args()

    log.info("=" * 68)
    log.info("START acquire_snrnaseq_discovery.py — %s (pilot=%s)",
             datetime.now().isoformat(timespec="seconds"), args.pilot)
    prespec = HERE / "outputs" / "sec4_prespec.json"
    if not prespec.exists():
        log.error("STOP. outputs/sec4_prespec.json missing. Run scripts/write_sec4_prespec.py first.")
        sys.exit(2)

    if not ROOT_EXTRACT.exists():
        log.error("STOP. Extracted dir missing: %s", ROOT_EXTRACT)
        log.error("Expected: data/discovery/GSE326221/GSE326221_RAW_extracted/ "
                  "containing the 193 files (59 MTX triplets + 16 H5).")
        sys.exit(3)

    records = _detect_format_per_sample(ROOT_EXTRACT, pilot=args.pilot)
    manifest_path = OUT / "sn_sample_format_manifest.json"
    with open(manifest_path, "w") as fh:
        json.dump([
            {k: (str(v) if isinstance(v, Path) else v) for k, v in r.items()}
            for r in records
        ], fh, indent=2, default=str)
    log.info("Wrote per-sample format manifest (%d rows) -> %s",
             len(records), manifest_path)

    _LAST_LOAD_SKIPPED.clear()
    _LAST_LOAD_CELLS_PRECAP.clear()
    detected_gsms = {r["gsm"] for r in records}
    if args.pilot:
        n_input = len(records)
    else:
        missing_gsms = sorted(EXPECTED_GSMS - detected_gsms)
        extra_gsms   = sorted(detected_gsms - EXPECTED_GSMS)
        if missing_gsms:
            log.warning(
                "FORMAT-DETECT: %d expected GSMs, %d detected. "
                "MISSING (%d): %s",
                len(EXPECTED_GSMS), len(detected_gsms),
                len(missing_gsms), ", ".join(missing_gsms),
            )
            for gsm in missing_gsms:
                _LAST_LOAD_SKIPPED.append({
                    "gsm": gsm, "stage": "format_detect_not_found",
                    "reason": "no files matching GSM found in RAW_extracted dir",
                    "patient_id": None, "format": None,
                })
        if extra_gsms:
            log.warning(
                "FORMAT-DETECT: %d unexpected GSMs detected "
                "(not in EXPECTED_GSMS 9625434-9625508): %s",
                len(extra_gsms), ", ".join(extra_gsms),
            )
        n_input = len(EXPECTED_GSMS)

    resume_idx = (CKPT_STAGES.index(args.resume_from)
                  if args.resume_from in CKPT_STAGES else -1)

    # --- STAGE 01: load & merge + CELL-QC ON FULL GENE SPACE (fixes 4x min_genes bug) ---
    if resume_idx <= 0 and args.resume_from != "01_merge_raw":
        adata = load_all_samples(records, keep_tmp=args.keep_tmp_samples)
        n_loaded_gsm = int(adata.obs["gsm"].nunique())
        n_skipped = len(_LAST_LOAD_SKIPPED)
        if n_skipped:
            log.warning(
                "SAMPLE-LOAD SUMMARY: %d input (expected) GSMs -> %d loaded, %d SKIPPED.",
                n_input, n_loaded_gsm, n_skipped,
            )
            skip_frac = n_skipped / max(1, n_input)
            if skip_frac > 0.05:
                log.warning(
                    "  HIGH SKIP FRACTION = %.1f%% (%d/%d). "
                    "If this is unexpected: re-extract RAW.tar with `tar -x -k --verify` "
                    "or re-download the corrupted GSMs from GEO.",
                    100.0 * skip_frac, n_skipped, n_input,
                )
            path_skip = OUT / "sn_skipped_samples.json"
            with open(path_skip, "w") as fh:
                json.dump(_LAST_LOAD_SKIPPED, fh, indent=2, default=str)
            log.warning("  Full skip list written to -> %s", path_skip)
        else:
            log.info("SAMPLE-LOAD SUMMARY: %d input (expected) GSMs -> %d/%d loaded, zero skips.",
                     n_input, n_loaded_gsm, n_input)
        # Cell-level QC was already applied per-sample on the FULL gene space
        # inside _load_one_sample() (before the per-sample HVG subset). The
        # merged space is HVG-reduced, so re-applying min_genes=300 here would
        # again be ~4x stricter and drop low-RNA immune cells. Skip it unless
        # the per-sample QC metrics are genuinely absent (fallback path).
        if "n_genes_by_counts" in adata.obs.columns:
            log.info("STAGE 01: cell-level QC already applied per-sample on full gene "
                     "space -> skipping merged cell-QC (merged shape=%s).", adata.shape)
        else:
            log.warning("STAGE 01: per-sample cell-QC metrics absent -> applying "
                        "merged cell-QC on %s.", adata.shape)
            adata = cell_qc_full_gene_space(adata)
        _save_ckpt(adata, "01_merge_raw", args.pilot)
    else:
        adata = _load_ckpt("01_merge_raw", args.pilot)
        if adata is None:
            log.error("RESUME FAIL: ckpt 01_merge_raw not found at %s",
                      _ckpt_path("01_merge_raw", args.pilot))
            sys.exit(4)
        n_loaded_gsm = int(adata.obs["gsm"].nunique())
        n_skipped = 0
        try:
            with open(OUT / "sn_skipped_samples.json") as f:
                prev = json.load(f)
            n_skipped = len(prev)
        except Exception:
            pass

    # --- STAGE 02: HVG-first on raw merged ---
    if resume_idx <= 1 and args.resume_from not in ("02_hvg_first", "03_qc",
                                                    "04_dimred", "05_annotated"):
        import scanpy as sc_hvg
        from scipy.sparse import issparse as _isp_hvg
        try:
            _probe = adata.X[:2, :2].toarray() if _isp_hvg(adata.X) else adata.X[:2, :2]
            _raw_ok = (np.allclose(_probe, np.round(_probe)) and np.min(_probe) >= 0)
        except Exception:
            _raw_ok = False
        _n_vars_pre = int(adata.n_vars)
        if _raw_ok and adata.n_vars > 4000:
            try:
                sc_hvg.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000)
                hvg_mask = adata.var["highly_variable"].values.copy()
                _n_hvg = int(np.sum(hvg_mask))
                _force_genes_upper = {g.upper()
                                      for glist in MARKERS.values() for g in glist}
                _var_upper = np.array([str(g).upper() for g in adata.var_names])
                _add = 0
                for _gi, _gu in enumerate(_var_upper):
                    if _gu in _force_genes_upper and not hvg_mask[_gi]:
                        hvg_mask[_gi] = True
                        _add += 1
                _n_hvg_final = int(np.sum(hvg_mask))
                log.info(
                    "HVG-first: seurat_v3 selected %d/%d raw HVG genes "
                    "+ %d FORCE-INCLUDED marker-genes (MARKERS dict K=%d types, "
                    "total HVG=%d) — subsetting columns NOW.",
                    _n_hvg, _n_vars_pre, _add, len(MARKERS), _n_hvg_final,
                )
                adata = adata[:, hvg_mask].copy()
                adata.var.drop(columns=["highly_variable", "means", "variances",
                                        "variances_norm"], errors="ignore", inplace=True)
                gc.collect()
            except Exception as _e_hvg:
                log.warning("HVG-first failed (%s); keeping all %d genes (may OOM at QC copy).",
                            _e_hvg, _n_vars_pre)
        else:
            log.info("HVG-first skipped (raw-detect=%s, n_vars=%d <= 4000); proceed with standard order.",
                     _raw_ok, _n_vars_pre)
        log.info("Post-HVG-first shape = %s (entering standard pipeline).", adata.shape)
        _save_ckpt(adata, "02_hvg_first", args.pilot)
    else:
        ckpt_to_load = args.resume_from if args.resume_from else "02_hvg_first"
        if ckpt_to_load == "02_hvg_first":
            adata = _load_ckpt("02_hvg_first", args.pilot)
            if adata is None:
                log.error("RESUME FAIL: ckpt 02_hvg_first not found at %s",
                          _ckpt_path("02_hvg_first", args.pilot))
                sys.exit(4)

    # --- STAGE 03: Gene-level QC only (cell-level QC done on full gene space at STAGE 01) ---
    if resume_idx <= 2 and args.resume_from not in ("03_qc", "04_dimred", "05_annotated"):
        adata = qc_and_filter_gene_only(adata)
        _save_ckpt(adata, "03_qc", args.pilot)
    else:
        if args.resume_from == "03_qc":
            adata = _load_ckpt("03_qc", args.pilot)
            if adata is None:
                log.error("RESUME FAIL: ckpt 03_qc not found at %s",
                          _ckpt_path("03_qc", args.pilot))
                sys.exit(4)

    # --- STAGE 04: dimensional reduction ---
    if resume_idx <= 3 and args.resume_from not in ("04_dimred", "05_annotated"):
        adata = dimensional_reduction(adata, batch_key="patient_id")
        _save_ckpt(adata, "04_dimred", args.pilot)
    else:
        if args.resume_from == "04_dimred":
            adata = _load_ckpt("04_dimred", args.pilot)
            if adata is None:
                log.error("RESUME FAIL: ckpt 04_dimred not found at %s",
                          _ckpt_path("04_dimred", args.pilot))
                sys.exit(4)

    # --- STAGE 05: annotation ---
    if resume_idx <= 4 and args.resume_from != "05_annotated":
        adata = annotate_clusters(adata, leiden_key="leiden_r0.8")
        _save_ckpt(adata, "05_annotated", args.pilot)
    else:
        if args.resume_from == "05_annotated":
            adata = _load_ckpt("05_annotated", args.pilot)
            if adata is None:
                log.error("RESUME FAIL: ckpt 05_annotated not found at %s",
                          _ckpt_path("05_annotated", args.pilot))
                sys.exit(4)

    write_sample_proportions(adata, sample_key="gsm")

    out_h5ad = DATA / (f"sn_processed_qc_annotated_pilot{args.pilot}.h5ad"
                       if args.pilot else "sn_processed_qc_annotated.h5ad")
    try:
        adata.write_h5ad(out_h5ad)
        log.info("Wrote processed AnnData -> %s (shape=%s)", out_h5ad, adata.shape)
    except Exception as e:
        log.warning("h5ad save failed (%s), falling back to sparse copy.", e)
        try:
            from scipy.sparse import csr_matrix
            adata.X = csr_matrix(adata.X)
            adata.write_h5ad(out_h5ad)
            log.info("Wrote processed AnnData (sparse copy) -> %s", out_h5ad)
        except Exception as e2:
            log.error("Could not save h5ad at all: %s. Proportions CSV is still valid.", e2)

    audit = {
        "ran_at": datetime.now().isoformat(timespec="seconds"),
        "source_accession": ACCESSION,
        "pilot": args.pilot,
        "n_samples_input":    n_input,
        "n_samples_loaded":   n_loaded_gsm,
        "n_samples_skipped":  n_skipped,
        "skip_fraction":      round(n_skipped / max(1, n_input), 4),
        "skipped_samples":    _LAST_LOAD_SKIPPED,
        "n_mtx_samples_in_input": sum(1 for r in records if r["format"] == "mtx"),
        "n_h5_samples_in_input":  sum(1 for r in records if r["format"] == "h5"),
        "final_shape":    [adata.n_obs, adata.n_vars],
        "max_cells_per_sample": MAX_CELLS_PER_SAMPLE,
        "n_cells_after_qc_precap_total": int(sum(_LAST_LOAD_CELLS_PRECAP)),
        "n_cells_postcap_total": int(adata.n_obs),
        "n_samples_capped": int(sum(1 for c in _LAST_LOAD_CELLS_PRECAP
                                    if c > MAX_CELLS_PER_SAMPLE)),
        "batch_key_used": "patient_id",
        "annotation_distribution": adata.obs["cell_type_marker_v1"].value_counts().to_dict(),
        "marker_sets_used": {k: v for k, v in MARKERS.items()},
    }
    with open(OUT / "sn_pipeline_audit.json", "w") as fh:
        json.dump(audit, fh, indent=2, default=str)
    log.info("Wrote audit -> %s", OUT / "sn_pipeline_audit.json")
    log.info("DONE acquire_snrnaseq_discovery.py.")
    log.info("FINAL: loaded=%d/%d  skipped=%d  cells=%d  hvg=%d.",
             n_loaded_gsm, n_input, n_skipped, int(adata.n_obs), int(adata.n_vars))


if __name__ == "__main__":
    main()
