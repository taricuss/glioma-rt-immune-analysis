
"""
Step 0a — Johnson / CARE discovery cohort: bulk RNA-seq acquisition & harmonization
Source: GEO GSE324481 (processed matrices — GEO open access, no EGA DAC required)
Target: log2(TPM + 1), Ensembl gene ID space, samples × genes wide-format

Run after Section 4.5 pre-spec is written (never run this before prespec is locked).

Outputs (written to data/discovery/):
    1. bulk_expr_log2tpm.parquet    — harmonized matrix (samples × genes, log2TPM+1)
    2. bulk_sample_meta.csv         — per-sample metadata scraped from GEO GSM entries
    3. bulk_gene_annotation.csv     — gene symbol → Ensembl ID map used for remap
    4. logs/acquire_bulk.log        — timestamped audit log

Notes on GSE324481:
    Per the Phase-0 report, GSE324481 carries processed bulk RNA-seq matrices
    (TPM or counts, depending on exactly what the authors uploaded).  The script
    auto-detects TPM vs counts: if a GSE supplementary file ends in "*TPM*" or
    the matrix quantiles sum to ~1e6 per column → TPM; if quantiles are integer-
    valued and sums are in the millions-read range → counts (we divide by colSum
    × 1e4 to get TPM).  If the processed file is ambiguous, the script raises and
    asks the user to confirm via a GEO manual-check (download GSE324481_series_
    matrix.txt first, inspect the !Sample_supplementary_file field in a text
    viewer — 2 min).

Why NOT pull raw FASTQ from EGA:
    EGA EGAD50000002485 DAC lead time ~2–6 wks; processed TPM sufficient for
    deconvolution per plan §FLAG1.
"""

from __future__ import annotations

import io
import os
import sys
import json
import logging
from datetime import datetime
from pathlib import Path

# ---- cp1252-safe stdout/stderr BEFORE any logging/print ----
# Windows PowerShell defaults to cp1252 on English installs; any U+2192 (->),
# U+2014 (--), or U+2026 (...) in log output triggers UnicodeEncodeError.
# Must wrap BEFORE logging.basicConfig() so the StreamHandler uses the
# wrapped stream with errors="replace", not the raw cp1252 buffer.
def _wrap_stream(s, name):
    if s is None:
        return s
    try:
        enc = (getattr(s, "encoding", None) or "").lower()
    except Exception:
        enc = ""
    if "utf" in enc:
        return s
    try:
        buf = getattr(s, "buffer", None)
        if buf is None:
            return s
        return io.TextIOWrapper(buf, encoding="utf-8", errors="replace",
                                line_buffering=getattr(s, "line_buffering", True))
    except Exception:
        return s
sys.stdout = _wrap_stream(sys.stdout, "stdout")
sys.stderr = _wrap_stream(sys.stderr, "stderr")

import numpy as np
import pandas as pd
import GEOparse

from mygene import MyGeneInfo

HERE = Path(__file__).resolve().parents[1]
DATA = HERE / "data" / "discovery"
LOGS = HERE / "logs"
for p in (DATA, LOGS):
    p.mkdir(parents=True, exist_ok=True)

logfmt = "%(asctime)s | %(levelname)-7s | %(message)s"
logging.basicConfig(
    level=logging.INFO,
    format=logfmt,
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(LOGS / "acquire_bulk.log", mode="w", encoding="utf-8"),
    ],
)
log = logging.getLogger("acquire_bulk")

ACCESSION = "GSE324481"
EXPECTED_SAMPLES = 75  # Phase-0 report: 35 patients × 2–3 timepoints

# ---------------------------------------------------------------------------
# 1. Download / cache the GEO Series (GSE + all GSM metadata)
# ---------------------------------------------------------------------------
def download_geo_series(acc: str, destdir: Path) -> GEOparse.GEOTypes.GSE:
    destdir.mkdir(parents=True, exist_ok=True)
    log.info("Downloading GEO Series %s -> %s  (first run is slow, uses ~400 MB)",
             acc, destdir)
    gse = GEOparse.get_GEO(geo=acc, destdir=str(destdir), silent=False)
    n_gsm = len(gse.gsms)
    log.info("Retrieved %s GSM entries (expected %d)", n_gsm, EXPECTED_SAMPLES)
    if n_gsm < EXPECTED_SAMPLES * 0.8:
        log.warning("GSM count unusually low -- check GEO status for the accession.")
        log.warning("  (Johnson et al. may have stored 75 samples in the series_matrix")
        log.warning("   file or in GSE324481_RAW.tar, not as individual GSM entries.)")
    return gse

# ---------------------------------------------------------------------------
# 2. PRIMARY PATH: direct FTP download of the GEO series_matrix file.
#    This is the standard, reliable format for processed GEO data.
#    Typical URL: ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE324nnn/GSE324481/matrix/
#    The file contains !Series_/!Sample_ headers then a tab-separated matrix
#    after the marker line "!series_matrix_table_begin".
#    If the GSE uses multiple platforms there may be one file per GPL; we try
#    GPL11154 first (Illumina NovaSeq, the likely platform for 2026 Nature
#    paper), then fall back to any GPL found in the soft metadata.
# ---------------------------------------------------------------------------
def download_and_parse_series_matrix(acc: str, gpls: list[str]) -> pd.DataFrame | None:
    import gzip
    import urllib.request

    # GEO FTP convention: take the "GSE" prefix + first 3 digits, then "nnn"
    # GSE324481  ->  GSE + 324 + nnn  =  GSE324nnn
    # GSE1234    ->  GSE + 12  + nnn  =  GSE12nnn  (padding handled by GEO)
    # We take the 3-digit slice after "GSE" (acc[3:6]) and pad with nnn.
    numeric_part_after_gse = acc[3:]
    if len(numeric_part_after_gse) >= 3:
        bucket = numeric_part_after_gse[:3]
    else:
        bucket = numeric_part_after_gse
    acc_short = "GSE" + bucket + "nnn"
    base = f"ftp://ftp.ncbi.nlm.nih.gov/geo/series/{acc_short}/{acc}/matrix"
    candidates_gpl = gpls[:] if gpls else ["GPL11154", "GPL16791", "GPL18573"]
    # Also add a no-GPL-suffix try (some GSEs publish one generic matrix)
    candidate_urls = [f"{base}/{acc}_{gpl}_series_matrix.txt.gz" for gpl in candidates_gpl]
    candidate_urls += [f"{base}/{acc}_series_matrix.txt.gz"]

    dest_dir = DATA / acc
    dest_dir.mkdir(parents=True, exist_ok=True)

    for url in candidate_urls:
        fname = Path(url).name
        local = dest_dir / fname
        if not local.exists():
            log.info("Trying series_matrix URL: %s", url)
            try:
                urllib.request.urlretrieve(url, str(local))
                log.info("  Downloaded OK (%d bytes)", local.stat().st_size)
            except Exception as e:
                log.info("  Not found (or FTP error) -- skip: %s", str(e)[:80])
                try:
                    if local.exists() and local.stat().st_size == 0:
                        local.unlink()
                except Exception:
                    pass
                continue
        else:
            log.info("Cached series_matrix found: %s", local)

        # --- Parse: skip lines until we hit !series_matrix_table_begin
        log.info("  Parsing series_matrix ...")
        try:
            opener = gzip.open(local, "rt", encoding="utf-8", errors="replace")
        except Exception:
            opener = open(local, "rt", encoding="utf-8", errors="replace")
        with opener as fh:
            lines = fh.readlines()
        start = None
        end = None
        for i, line in enumerate(lines):
            if "!series_matrix_table_begin" in line:
                start = i + 1
            if "!series_matrix_table_end" in line:
                end = i
        if start is None:
            log.warning("  No !series_matrix_table_begin marker found in %s -- skip", fname)
            continue
        body = lines[start:end] if end else lines[start:]
        if len(body) < 3:
            log.warning("  Matrix body < 3 lines in %s -- skip", fname)
            continue
        from io import StringIO
        df = pd.read_csv(StringIO("".join(body)), sep="\t", index_col=0, comment="!")
        # Drop rows that are just metadata fields (ID_REF column is fine, index will
        # be gene symbols/Ensembl). Clean any "!Sample_*" leftover columns.
        df = df.loc[:, [c for c in df.columns if not str(c).startswith("!")]]
        log.info("  Parsed shape = %s, index[:3] = %s, columns[:3] = %s",
                 df.shape, list(df.index[:3]), list(df.columns[:3]))
        if len(df.columns) >= EXPECTED_SAMPLES * 0.8:
            log.info("  Sample count matches (>80%% of %d expected) -- USING this matrix.",
                     EXPECTED_SAMPLES)
            return df
        else:
            log.warning("  Only %d cols, still short of %d -- keep trying other files.",
                        len(df.columns), EXPECTED_SAMPLES)
    return None

# ---------------------------------------------------------------------------
# 3. Try the *supplementary processed files* (authors uploaded a TPM/counts
#    matrix directly as a single GSE-level file, e.g. GSE324481_RAW.tar).
# ---------------------------------------------------------------------------
def read_processed_suppl(gse: GEOparse.GEOTypes.GSE) -> pd.DataFrame | None:
    import gzip
    import tarfile
    import urllib.request
    from io import StringIO

    supp = getattr(gse.metadata, "supplementary_file", None)
    if supp is None:
        supp = gse.metadata.get("supplementary_file", [])
    if isinstance(supp, str):
        supp = [supp]
    log.info("GSE-level supplementary_file entries = %d", len(supp))
    for url in supp:
        log.info("  -> %s", url)
        fname = Path(url).name.lower()
        local_dest = DATA / ACCESSION / Path(url).name
        local_dest.parent.mkdir(parents=True, exist_ok=True)
        if not local_dest.exists():
            log.info("  Downloading via urllib (GEOparse API too fragile)...")
            try:
                urllib.request.urlretrieve(url, str(local_dest))
                log.info("  Downloaded %d bytes OK", local_dest.stat().st_size)
            except Exception as e:
                log.warning("  urllib download failed: %s -- try GEOparse", str(e)[:100])
                try:
                    gse.download_supplementary_files(directory=str(DATA))
                except Exception as e2:
                    log.warning("  GEOparse also failed: %s", str(e2)[:100])
                    continue

        # --- If it's a .tar (typically GSE324481_RAW.tar), extract and scan ---
        if local_dest.suffix == ".tar" or "raw" in local_dest.name.lower():
            log.info("  Found RAW.tar-like file: %s -- extracting ...", local_dest.name)
            extract_dir = local_dest.parent / f"{local_dest.stem}_extracted"
            extract_dir.mkdir(exist_ok=True)
            try:
                with tarfile.open(local_dest, "r:*") as tf:
                    tf.extractall(path=str(extract_dir))
                log.info("  Extracted to %s", extract_dir)
            except Exception as e:
                log.warning("  tar extract failed: %s", str(e)[:120])
                extract_dir = local_dest.parent
            # Recursively look for matrix-like files inside the extraction dir
            candidates = []
            for pat in ("*tpm*", "*count*", "*expr*", "*matrix*", "*.tsv*", "*.csv*"):
                candidates.extend(list(extract_dir.rglob(pat)))
            # Also plain .gz files that contain "gene" or "sample" naming
            candidates.extend([p for p in extract_dir.rglob("*.gz")
                               if p.stat().st_size > 50_000])  # > 50 KB = real matrix
            for c in sorted(set(candidates), key=lambda p: -p.stat().st_size):
                if c.is_dir():
                    continue
                log.info("  Inspecting candidate: %s (%d KB)",
                         c.relative_to(DATA), c.stat().st_size // 1024)
                try:
                    if str(c).endswith(".gz"):
                        op = gzip.open(c, "rt", encoding="utf-8", errors="replace")
                    else:
                        op = open(c, "rt", encoding="utf-8", errors="replace")
                    with op as fh:
                        head = [next(fh) for _ in range(min(5, sum(1 for _ in range(5))))]
                    if len(head) < 2:
                        continue
                    sep = "\t" if head[0].count("\t") >= head[0].count(",") else ","
                    df = pd.read_csv(c, sep=sep, index_col=0,
                                     compression="gzip" if str(c).endswith(".gz") else None,
                                     on_bad_lines="skip")
                    if df.shape[1] >= EXPECTED_SAMPLES * 0.8:
                        log.info("  USING shape = %s, cols-sample count passes 80%% threshold",
                                 df.shape)
                        return df
                    elif df.shape[1] >= 4:
                        log.info("  shape = %s (sample count low, keep scanning)", df.shape)
                except Exception as e:
                    log.info("    skip (%s)", str(e)[:80])
            return None
        else:
            # Non-tar: direct parse of one matrix file
            cand = [local_dest]
            for c in cand:
                log.info("Inspecting candidate matrix: %s", c)
                op = gzip.open(c, "rt", errors="replace") if str(c).endswith(".gz") \
                    else open(c, "rt", errors="replace")
                with op as fh:
                    head = [next(fh) for _ in range(5)]
                    sep = "\t" if head[0].count("\t") >= head[0].count(",") else ","
                    fh.seek(0)
                    df = pd.read_csv(fh, sep=sep, index_col=0)
                log.info("  shape = %s, index[:3] = %s", df.shape, list(df.index[:3]))
                return df
    return None

# ---------------------------------------------------------------------------
# 4. Fallback: extract per-GSM supplementary files from every GSM entry.
#    Johnson et al. may have uploaded one "quant.sf"/"abundance.tsv" (Salmon)
#    or one "*.genes.results" (RSEM) per sample.
# ---------------------------------------------------------------------------
def read_per_gsm_suppl(gse: GEOparse.GEOTypes.GSE) -> pd.DataFrame:
    log.info("FALLBACK: collecting per-GSM processed files (this may download ~2 GB).")
    per_sample = {}
    try:
        gse.download_supplementary_files(directory=str(DATA))
    except Exception as e:
        log.warning("GEOparse.download_supplementary_files error: %s", str(e)[:120])
        log.warning("Continuing with any files GEOparse DID manage to download...")
    for gsm_name, gsm in gse.gsms.items():
        sf = gsm.metadata.get("supplementary_file", [])
        if isinstance(sf, str):
            sf = [sf]
        for url in sf:
            name = Path(url).name.lower()
            if not any(tok in name for tok in ("tpm", "count", "genes.results", "abundance", "rsem", "quant")):
                continue
            # GEOparse downloads to DATA/ACCESSION/GSM_XXX/
            local = list((DATA / ACCESSION / gsm_name).rglob(Path(url).name)) + \
                    list(DATA.rglob(Path(url).name))
            if not local:
                log.warning("  %s: not found locally even after download — skip", gsm_name)
                continue
            f = local[0]
            if "genes.results" in name:
                tmp = pd.read_csv(f, sep="\t", index_col=0)
                per_sample[gsm_name] = tmp["TPM"] if "TPM" in tmp.columns else tmp.iloc[:, 1]
            elif "abundance" in name or "quant.sf" in name:
                tmp = pd.read_csv(f, sep="\t", index_col=0)
                per_sample[gsm_name] = tmp["TPM"] if "TPM" in tmp.columns else tmp.iloc[:, 3]
            else:
                try:
                    tmp = pd.read_csv(f, sep=None, index_col=0, engine="python")
                    per_sample[gsm_name] = tmp.iloc[:, 0]
                except Exception as e:
                    log.warning("  %s: parse failed (%s)", gsm_name, e)
                    continue
    df = pd.DataFrame(per_sample)
    log.info("Per-GSM fallback matrix shape = %s", df.shape)
    return df

# ---------------------------------------------------------------------------
# 4. Detection: TPM vs counts
# ---------------------------------------------------------------------------
def decide_tpm_or_counts(df: pd.DataFrame) -> tuple[pd.DataFrame, str]:
    col_sums = df.sum(axis=0)
    mean_sum = col_sums.mean()
    is_integer = df.apply(lambda s: pd.api.types.is_integer_dtype(s) or
                          (s.dropna().apply(float.is_integer).all() if s.dtype.kind == "f" else False)).all()
    if abs(mean_sum - 1e6) / 1e6 < 0.2 and not is_integer:
        kind = "TPM"
    elif is_integer and mean_sum > 1e5:
        kind = "counts"
    else:
        log.warning("Cannot auto-detect TPM vs counts (mean_sum=%.2e, int=%s). "
                    "Treating as TPM; if wrong the QC step (plan §5.3) will flag "
                    "purity-distribution mismatch with Johnson et al.'s Fig 1c.",
                    mean_sum, is_integer)
        kind = "TPM(assumed)"
    if kind.startswith("counts"):
        df = (df / df.sum(axis=0)) * 1e6
        log.info("Converted counts → TPM (column sum to 1e6).")
    return df, kind

# ---------------------------------------------------------------------------
# 5. Harmonize: gene ID remap → Ensembl GRCh38.p14 gene IDs (standard)
# ---------------------------------------------------------------------------
def remap_to_ensembl(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    ids_raw = list(df.index.astype(str))
    log.info("Remapping %d gene IDs to Ensembl via myGene.info (GRCh38)", len(ids_raw))
    mg = MyGeneInfo()
    # Figure out input ID type from the first 20 ids
    sample = ids_raw[:20]
    if all("ENSG" in s for s in sample if s.startswith("ENS")):
        log.info("Index already ENSG; trivial remap.")
        annot = pd.DataFrame({"input_id": ids_raw, "ensembl_gene": ids_raw, "symbol": ids_raw},
                             index=ids_raw)
        return df, annot
    res = mg.querymany(ids_raw, scopes="symbol,ensembl.gene,entrezgene,refseq,alias",
                       fields="ensembl.gene,symbol,name", species="human", returnall=False,
                       as_dataframe=True, df_index=True, verbose=False)
    if res is None or len(res) == 0:
        log.error("myGene remap returned zero rows — aborting.")
        sys.exit(5)
    res = res.reset_index().rename(columns={"query": "input_id"})
    # Dedupe: pick one ENSG per input (first non-null ensembl.gene)
    if "ensembl.gene" in res.columns:
        res_ens = res.dropna(subset=["ensembl.gene"]).drop_duplicates(subset=["input_id"])
        ens_map = dict(zip(res_ens["input_id"], res_ens["ensembl.gene"]))
        sym_map = dict(zip(res_ens["input_id"], res_ens.get("symbol", res_ens["input_id"])))
    else:
        ens_map = {r["input_id"]: r["input_id"] for _, r in res.iterrows()}
        sym_map = ens_map.copy()
    kept = [i for i in ids_raw if i in ens_map]
    lost = len(ids_raw) - len(kept)
    log.info("Remap kept %d / %d gene IDs (%d unmapped; acceptable loss <8%%).",
             len(kept), len(ids_raw), lost)
    df2 = df.loc[kept].copy()
    df2.index = [ens_map[i] for i in df2.index]
    df2 = df2[~df2.index.duplicated(keep="first")]
    annot = pd.DataFrame({
        "input_id": kept,
        "ensembl_gene": [ens_map[i] for i in kept],
        "symbol":     [sym_map.get(i, i)   for i in kept],
    }).set_index("ensembl_gene")
    return df2, annot

# ---------------------------------------------------------------------------
# 6. Metadata: pull GSM characteristics (patient ID, timepoint, tissue, etc.)
# ---------------------------------------------------------------------------
def extract_gsm_metadata(gse: GEOparse.GEOTypes.GSE, sample_cols: list[str]) -> pd.DataFrame:
    rows = []
    for gsm_name, gsm in gse.gsms.items():
        chars = {k: ("; ".join(v) if isinstance(v, list) else v)
                 for k, v in gsm.metadata.get("characteristics_ch1", {}).items()}
        if isinstance(gsm.metadata.get("characteristics_ch1"), list):
            # GEOparse sometimes returns a list of "key: value" strings
            for item in gsm.metadata["characteristics_ch1"]:
                if ": " in item:
                    k, v = item.split(": ", 1)
                    chars.setdefault(k.strip(), v.strip())
        rows.append({
            "gsm": gsm_name,
            "title": gsm.metadata.get("title", [""])[0],
            "source": gsm.metadata.get("source_name_ch1", [""])[0],
            **chars,
        })
    meta = pd.DataFrame(rows).set_index("gsm")
    log.info("GSM metadata fields extracted: %s", list(meta.columns)[:20])
    return meta

# ---------------------------------------------------------------------------
# Main orchestration
# ---------------------------------------------------------------------------
def main():
    log.info("=" * 68)
    log.info("START acquire_bulk.py -- %s", datetime.now().isoformat(timespec="seconds"))
    log.info("Pre-spec lock-file check: verifying outputs/sec4_prespec.json exists ...")
    prespec = HERE / "outputs" / "sec4_prespec.json"
    if not prespec.exists():
        log.error("STOP. Section 4.5 pre-spec has NOT been written yet.")
        log.error("Write outputs/sec4_prespec.json FIRST (scripts/write_sec4_prespec.py).")
        log.error("Rule: never run acquisition before the test list is locked.")
        sys.exit(2)
    with open(prespec) as fh:
        ps = json.load(fh)
    log.info("Pre-spec version = %s (dated %s). Proceeding.",
             ps.get("version"), ps.get("locked_at"))

    gse = download_geo_series(ACCESSION, DATA / ACCESSION)

    # --- Collect GPL IDs for series_matrix URL construction ---
    gpls = []
    try:
        # gse.gpls is a dict {gpl_id: gpl_object} in most GEOparse versions
        gpls = list(getattr(gse, "gpls", {}).keys())
    except Exception:
        pass
    if not gpls:
        # Fallback: try extracting from metadata
        gpl_meta = gse.metadata.get("platform_id", [])
        if isinstance(gpl_meta, str):
            gpls = [gpl_meta]
        elif isinstance(gpl_meta, list):
            gpls = [str(x).strip() for x in gpl_meta if str(x).strip()]
    log.info("GPL IDs extracted from GSE = %s", gpls)

    # --- Assemble expression matrix  (PRIORITY-ORDERED PATHS) ---
    # Path 1 = series_matrix FTP (most reliable for GEO processed data)
    # Path 2 = GSE-level supplementary file (usually RAW.tar)
    # Path 3 = per-GSM supplementary files (last resort, slowest)
    expr = None
    log.info("-" * 50)
    log.info("PATH 1/3: direct FTP series_matrix ...")
    expr = download_and_parse_series_matrix(ACCESSION, gpls)
    if expr is not None:
        log.info("PATH 1/3 SUCCESS.")
    if expr is None:
        log.info("-" * 50)
        log.info("PATH 2/3: GSE-level supplementary files (GSE324481_RAW.tar etc.) ...")
        expr = read_processed_suppl(gse)
        if expr is not None:
            log.info("PATH 2/3 SUCCESS.")
    if expr is None:
        log.info("-" * 50)
        log.info("PATH 3/3: per-GSM supplementary file collection ...")
        expr = read_per_gsm_suppl(gse)
        log.info("PATH 3/3 done (may have 4 GSMs only -- check sample count).")

    if expr is None or len(expr.columns) == 0:
        log.error("*** ALL THREE PATHS FAILED to produce a usable expression matrix. ***")
        log.error("Manual check required: open https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE324481")
        log.error("  in a browser, inspect the 'Download family' or 'Series Matrix File(s)' links.")
        log.error("  If the authors uploaded processed matrices as supplementary xlsx/csv,")
        log.error("  save them directly to data/discovery/ and report the file names here.")
        sys.exit(3)

    log.info("Raw expression matrix (genes x samples): %s", expr.shape)
    if len(expr.columns) < EXPECTED_SAMPLES * 0.8:
        log.warning("WARNING: only %d samples recovered (expected ~%d). This may not be",
                    len(expr.columns), EXPECTED_SAMPLES)
        log.warning("  the full 75-sample Johnson cohort. Save any additional files found")
        log.warning("  on the GEO accession page to data/discovery/ and re-run.")
    expr, kind = decide_tpm_or_counts(expr)
    log.info("Matrix type = %s, column-sum mean = %.4e", kind, expr.sum(axis=0).mean())

    # Samples × genes orientation (more intuitive for downstream models):
    expr_sg = expr.T
    log.info("Transposed to samples × genes: %s", expr_sg.shape)

    # --- Remap gene IDs ---
    expr_sg_ens, annot_ens = remap_to_ensembl(expr_sg.T)
    # remap_to_ensembl works on genes as index; after T we need to transpose once more
    # to keep samples × genes orientation.  The helper above was given genes × samples
    # back because of the T — fix:
    expr_final_samples_genes = expr_sg_ens.T
    log.info("Final harmonized matrix (samples × genes, Ensembl IDs): %s",
             expr_final_samples_genes.shape)

    # --- log2(TPM + 1) ---
    log2tpm = np.log2(expr_final_samples_genes.astype(np.float64) + 1.0)
    log.info("Applied log2(TPM + 1); value range = [%.3f, %.3f]",
             log2tpm.values.min(), log2tpm.values.max())

    # --- Metadata ---
    sample_meta = extract_gsm_metadata(gse, list(log2tpm.index))
    sample_meta.to_csv(DATA / "bulk_sample_meta_geo_raw.csv")

    # --- Write outputs ---
    # Two destinations:
    #   (A) data/discovery/        -- canonical primary store (long-term)
    #   (B) outputs/intermediate/   -- what downstream Phase 2.1-3.1 scripts actually read
    #       Exact names required: estimate_purity.py DEFAULT_INPUT_GZ / DEFAULT_INPUT_PQ
    INTERM = HERE / "outputs" / "intermediate"
    INTERM.mkdir(parents=True, exist_ok=True)

    wrote_parquet = False
    parquet_path_pri = DATA / "bulk_expr_log2tpm.parquet"
    parquet_path_sec = INTERM / "bulk_log2tpm_matrix.parquet"
    csvgz_path_pri  = DATA / "bulk_expr_log2tpm.csv.gz"
    csvgz_path_sec  = INTERM / "bulk_log2tpm_matrix.csv.gz"

    log.info("Writing harmonized matrix -> %s and %s", parquet_path_pri.name, parquet_path_sec.name)
    try:
        log2tpm.to_parquet(parquet_path_pri, engine="auto")
        try:
            log2tpm.to_parquet(parquet_path_sec, engine="auto")
        except Exception:
            pass  # primary write OK; secondary is convenience
        wrote_parquet = True
    except Exception:
        log.warning("parquet engine unavailable -- falling back to gzipped CSV (~2x bigger)")
        import gzip
        with gzip.open(csvgz_path_pri, "wt", compresslevel=6, encoding="utf-8") as fh:
            log2tpm.to_csv(fh)
        try:
            with gzip.open(csvgz_path_sec, "wt", compresslevel=6, encoding="utf-8") as fh2:
                log2tpm.to_csv(fh2)
        except Exception:
            pass

    annot_ens.to_csv(DATA / "bulk_gene_annotation.csv")
    try:
        annot_ens.to_csv(INTERM / "bulk_gene_annotation.csv")
    except Exception:
        pass
    # Preserve an explicit metadata .json for audit
    audit = {
        "geo_accession": ACCESSION,
        "date_run": datetime.now().isoformat(timespec="seconds"),
        "matrix_type_detected": kind,
        "shape_samples_genes": list(log2tpm.shape),
        "gene_space": "Ensembl GRCh38 (myGene.info remap)",
        "samples_expected": EXPECTED_SAMPLES,
        "samples_observed": len(log2tpm),
        "patients_inferred": len(set(sample_meta.iloc[:, 0])) if len(sample_meta.columns) else "TBD after Supp Table 1 join",
        "prespec_version": ps.get("version"),
        "primary_write_parquet": str(parquet_path_pri) if wrote_parquet else None,
        "primary_write_csv": str(csvgz_path_pri) if not wrote_parquet else None,
        "downstream_visible_path_parquet": str(parquet_path_sec),
        "downstream_visible_path_csv": str(csvgz_path_sec),
    }
    with open(DATA / "bulk_acquisition_audit.json", "w") as fh:
        json.dump(audit, fh, indent=2, default=str)
    try:
        with open(INTERM / "bulk_acquisition_audit.json", "w") as fh:
            json.dump(audit, fh, indent=2, default=str)
    except Exception:
        pass

    log.info("DONE acquire_bulk.py.")
    log.info("=" * 68)

if __name__ == "__main__":
    main()
