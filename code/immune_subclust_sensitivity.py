#!/usr/bin/env python
"""
SENSITIVITY: immune-compartment subclustering on the existing cap=3000 snRNA
checkpoints.  This is an ANNOTATION DEVIATION from pre-spec v1.0, not a primary
analysis.  Read-only w.r.t. every locked/primary output; each artefact written
by this script carries the "_immune_subclust" suffix.

The design below was pre-declared before any result existed and must not be
adjusted after seeing results:

  1. compartment (union):
     (a) all cells in leiden_r0.8 clusters whose primary label is
         Microglia / MoMF / NK / B / CD8_T; plus
     (b) any cell, in any cluster, with PTPRC log-normalised expression > 0
         (skipped, and reported, if PTPRC is absent from the 2000-gene space).
     STOP if the compartment exceeds 80,000 cells (RAM guard).
  2. re-embedding on the subset only, from the existing 2000 genes:
     scale(max_value=10) -> PCA 30 -> Harmony(patient_id) -> neighbours k=15 ->
     Leiden at resolutions 0.5 / 1.0 / 1.5.  PRIMARY resolution = 1.0;
     0.5 and 1.5 are stability checks only.  Seed = 7331 everywhere.
  3. annotation with the pipeline's own MARKERS dict (imported, not retyped);
     sc.tl.score_genes(ctrl_size=50).  A primary type counts as RESOLVED at a
     resolution iff >= 1 cluster satisfies ALL of:
       - highest cluster-mean score in that cluster AND margin >= 0.10 over the
         second-best type;
       - cluster has >= 50 cells;
       - >= 2 of the type's marker genes are DE-up vs the rest
         (rank_genes_groups, wilcoxon, logFC > 0.25, BH padj < 0.05);
       - not a doublet/mixture: second-best score < 50% of the winner's score.
     A type is EVALUABLE iff RESOLVED AND mapped cells >= 50 AND >= 10 samples
     each carry >= 10 of those cells.
  4. cluster labels are mapped back to cells; per-sample proportions use each
     sample's TOTAL cap=3000 cell count (221,016-cell denominator, not the
     subset), so they are comparable with sn_proportions_per_sample.csv.

PRE-DECLARED GATE (applied mechanically, printed at the end):
     EVALUABLE >= 6/10 primary at res 1.0 AND >= 5/10 at BOTH 0.5 and 1.5
     -> report and ASK USER first.
     otherwise -> STOP: subclustering does not recover the family.

No LMM and no permutation null are run here.
"""
from __future__ import annotations

import gc
import hashlib
import json
import logging
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

OUT = ROOT / "outputs" / "intermediate"
LOGS = ROOT / "logs"
CKPT = ROOT / "data" / "discovery" / "checkpoints" / "ckpt_05_annotated.h5ad"
PRESPEC = ROOT / "outputs" / "sec4_prespec.json"

for _p in (OUT, LOGS):
    _p.mkdir(parents=True, exist_ok=True)

# Resource guards must be set before numpy/scanpy grab their thread pools.
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("NUMBA_CACHE_DIR", str(ROOT / ".dgnc_env" / "numba_cache"))
Path(os.environ["NUMBA_CACHE_DIR"]).mkdir(parents=True, exist_ok=True)

LOG_PATH = LOGS / "immune_subclust.log"
# NOTE: logging is configured BEFORE importing the pipeline module.  The module
# calls logging.basicConfig(..., FileHandler(logs/acquire_snrnaseq.log, mode="w"))
# at import time; because the root logger already has handlers here, that call
# is a no-op and the pipeline log is NOT truncated.
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    handlers=[logging.StreamHandler(sys.stdout),
              logging.FileHandler(LOG_PATH, mode="w")],
)
log = logging.getLogger("immune_subclust")

MIN_FREE_DISK_GB = 3.0
COMPARTMENT_STOP = 80_000
SEED = 7331
RESOLUTIONS = (0.5, 1.0, 1.5)
PRIMARY_RES = 1.0
MIN_CLUSTER_CELLS = 50
MIN_MARGIN = 0.10
MAX_RATIO = 0.50
MIN_DE_MARKERS = 2
DE_LOGFC = 0.25
DE_PADJ = 0.05
MIN_MAPPED_CELLS = 50
MIN_SAMPLES_GE10 = 10

SUFFIX = "_immune_subclust"
T0 = time.time()


def _peak_rss_gb() -> float:
    try:
        import psutil
        return psutil.Process().memory_info().rss / 1e9
    except Exception:
        return float("nan")


_rss_peak = [_peak_rss_gb()]


def step(msg: str) -> None:
    """Log a progress marker and track peak RSS."""
    _rss_peak[0] = max(_rss_peak[0], _peak_rss_gb())
    log.info("[%6.1fs | rss~%.2f GB] %s", time.time() - T0, _rss_peak[0], msg)


def sha256_file(path: Path, chunk: int = 8 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def free_disk_gb() -> float:
    return shutil.disk_usage(str(ROOT.anchor or "<PROJECT_ROOT>")).free / 1e9


def git_head() -> str | None:
    try:
        r = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                           capture_output=True, text=True, timeout=20)
        return r.stdout.strip() or None if r.returncode == 0 else None
    except Exception:
        return None


# ---------------------------------------------------------------------------
# 0. Pre-flight: disk gate
# ---------------------------------------------------------------------------
free_start = free_disk_gb()
log.info("=" * 78)
log.info("START immune_subclust_sensitivity.py  (SENSITIVITY: annotation "
         "deviation from pre-spec v1.0)")
log.info("free disk at start: %.2f GB (gate %.2f GB)", free_start, MIN_FREE_DISK_GB)
if free_start < MIN_FREE_DISK_GB:
    log.error("ABORT: < %.2f GB free on C:.", MIN_FREE_DISK_GB)
    sys.exit(2)

for _p, _label in ((CKPT, "input checkpoint"), (PRESPEC, "pre-spec")):
    if not _p.exists():
        log.error("ABORT: %s missing: %s", _label, _p)
        sys.exit(2)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.sparse import issparse  # noqa: E402

# The pipeline module builds its basicConfig() argument list at import time,
# which evaluates logging.FileHandler(<acquire log>, mode="w") EAGERLY and thus
# opens+truncates that log regardless of whether basicConfig later installs the
# handler.  Pre-configuring logging (above) therefore CANNOT protect the file.
# Back the log up immediately before the import and restore it afterwards.
_ACQ_LOG = LOGS / "acquire_snrnaseq.log"
_acq_backup = None
if _ACQ_LOG.exists() and _ACQ_LOG.stat().st_size > 0:
    _acq_backup = LOGS / "_acquire_snrnaseq.log.presubclust.bak"
    shutil.copy2(_ACQ_LOG, _acq_backup)
try:
    import acquire_snrnaseq_discovery as mainpipe  # noqa: E402
finally:
    if _acq_backup is not None:
        shutil.copy2(_acq_backup, _ACQ_LOG)
        _acq_backup.unlink()
        log.info("restored logs/acquire_snrnaseq.log from pre-import backup "
                 "(%d bytes)", _ACQ_LOG.stat().st_size)
MARKERS: dict[str, list[str]] = mainpipe.MARKERS

with open(PRESPEC) as fh:
    prespec = json.load(fh)
PRIMARY: list[str] = list(prespec["a_cell_types_primary"])
EXPLORATORY: list[str] = list(prespec["a_cell_types_exploratory"])
PRE_SPEC_SEED = prespec["d_permutation"]["seed"]
assert len(PRIMARY) == 10, PRIMARY
missing_in_markers = [t for t in PRIMARY if t not in MARKERS]
assert not missing_in_markers, missing_in_markers
step(f"MARKERS imported from pipeline module ({len(MARKERS)} types); "
     f"K=10 primary family verified against sec4_prespec.json "
     f"(pre-spec permutation seed={PRE_SPEC_SEED}, using SEED={SEED})")

# ---------------------------------------------------------------------------
# 1. Load the annotated checkpoint (in memory; 221,016 x 2,000)
# ---------------------------------------------------------------------------
import scanpy as sc  # noqa: E402

sc.settings.verbosity = 1
sc.settings.autosave = False
sc.settings.autoshow = False

INPUT_SHA = sha256_file(CKPT)
step(f"sha256(ckpt_05_annotated.h5ad) = {INPUT_SHA[:16]}... ({CKPT.stat().st_size/1e6:.1f} MB)")

adata = sc.read_h5ad(CKPT)
step(f"loaded {CKPT.name}: shape={adata.shape}")

obs_needed = ["gsm", "patient_id", "leiden_r0.8", "cell_type_marker_v1"]
absent = [c for c in obs_needed if c not in adata.obs.columns]
if absent:
    log.error("ABORT: obs columns missing: %s (have: %s)", absent, list(adata.obs.columns))
    sys.exit(2)

if issparse(adata.X):
    adata.X = adata.X.tocsr().astype(np.float32, copy=False)

raw_same = None
if adata.raw is not None:
    a = adata.X[:200, :]
    b = adata.raw.X[:200, :]
    a = a.toarray() if issparse(a) else np.asarray(a)
    b = b.toarray() if issparse(b) else np.asarray(b)
    raw_same = bool(np.allclose(a, b))
    log.info("X == raw.X on first 200 cells: %s  (X sparse=%s, raw sparse=%s)",
             raw_same, issparse(adata.X), issparse(adata.raw.X))
    if raw_same:
        adata.raw = None
        gc.collect()
    else:
        log.warning("X != raw.X -> using adata.raw as the expression source.")

src = adata
if not raw_same and adata.raw is not None:
    src = adata.raw.to_adata()

var_names = list(src.var_names)
step(f"expression source: shape={src.shape}, n_genes={len(var_names)}")

log.info("cell_type_marker_v1 value_counts: %s",
         adata.obs["cell_type_marker_v1"].value_counts().to_dict())
log.info("leiden_r0.8 n_clusters=%d ; gsm n=%d ; patient_id n=%d",
         adata.obs["leiden_r0.8"].nunique(), adata.obs["gsm"].nunique(),
         adata.obs["patient_id"].nunique())

# Per-sample TOTAL cell count over the FULL checkpoint (all 221,016 cells).
# This is the required proportion denominator; it is captured here because
# `adata` is deleted once the immune subset has been built.  NOTE: the primary
# sn_proportions_per_sample.csv holds PROPORTIONS (its rows sum to ~1.0), so its
# row sums are NOT cell counts and must not be used as the denominator.
SAMPLE_TOTALS = adata.obs["gsm"].value_counts()

# ---------------------------------------------------------------------------
# 2. Compartment rule
# ---------------------------------------------------------------------------
SEED_TYPES = ["Microglia_P2RY12_TMEM119", "Monocyte_derived_macrophages_CD163",
              "NK_cells", "B_cells", "CD8_T_cells"]
cluster_label = (adata.obs.groupby("leiden_r0.8", observed=True)["cell_type_marker_v1"]
                 .agg(lambda s: s.iloc[0]))
seed_clusters = cluster_label.index[cluster_label.isin(SEED_TYPES)].tolist()
mask_a = adata.obs["leiden_r0.8"].isin(seed_clusters).to_numpy()
step(f"compartment (a): {len(seed_clusters)} leiden_r0.8 clusters labelled "
     f"{SEED_TYPES} -> {int(mask_a.sum())} cells")

ptprc_present = "PTPRC" in var_names
if ptprc_present:
    j = var_names.index("PTPRC")
    col = src.X[:, j]
    col = np.asarray(col.toarray() if issparse(col) else col).ravel()
    mask_b = col > 0
    step(f"compartment (b): PTPRC present in gene space -> {int(mask_b.sum())} cells "
         f"with log-norm PTPRC > 0")
else:
    mask_b = np.zeros(adata.n_obs, dtype=bool)
    log.warning("PTPRC NOT in var_names -> rule (b) SKIPPED (reported as such).")

comp_mask = mask_a | mask_b
n_comp = int(comp_mask.sum())
comp_fraction = n_comp / adata.n_obs
step(f"compartment = {n_comp} cells ({comp_fraction:.1%} of {adata.n_obs}); "
     f"rule(a)={int(mask_a.sum())}, rule(b)={int(mask_b.sum())}, "
     f"overlap(a&b)={int((mask_a & mask_b).sum())}")

if n_comp > COMPARTMENT_STOP:
    log.error("ABORT (RAM guard): compartment %d > %d cells.", n_comp, COMPARTMENT_STOP)
    with open(OUT / f"immune_subclust_audit{SUFFIX}.json", "w") as fh:
        json.dump({"aborted": "compartment_gt_80000", "compartment_cells": n_comp,
                   "rule_a_cells": int(mask_a.sum()), "rule_b_cells": int(mask_b.sum()),
                   "ptprc_present": ptprc_present}, fh, indent=2)
    sys.exit(3)

comp_labels = adata.obs.loc[comp_mask, "cell_type_marker_v1"].value_counts().to_dict()
log.info("compartment primary-label composition: %s", comp_labels)

# ---------------------------------------------------------------------------
# 3. Subset -> re-embed
# ---------------------------------------------------------------------------
idx = np.flatnonzero(comp_mask)
sub = src[idx, :].copy()
if src is not adata:
    del src
del adata
gc.collect()
step(f"subset built: {sub.shape}; sparse={issparse(sub.X)}")

# raw = independent copy of the log-normalised matrix so that scale() can never
# touch what score_genes / rank_genes_groups read.
sub.raw = sub.copy()
gc.collect()

if issparse(sub.X):
    sub.X = sub.X.tocsr().astype(np.float32, copy=False)

for cell_type, genes in MARKERS.items():
    genes_in = [g.upper() for g in genes if g.upper() in sub.var_names]
    if len(genes_in) < 3:
        log.warning("  skip %s: only %d/%d markers in var_names", cell_type,
                    len(genes_in), len(genes))
        continue
    sc.tl.score_genes(sub, gene_list=genes_in, score_name=f"score_{cell_type}",
                      use_raw=True, ctrl_size=50)
score_cols = [c for c in sub.obs.columns if c.startswith("score_")]
if not score_cols:
    log.error("ABORT: zero score columns.")
    sys.exit(2)
step(f"score_genes done: {len(score_cols)} score columns")

sc.pp.scale(sub, max_value=10)
gc.collect()
step("scale(max_value=10) done")

sc.tl.pca(sub, n_comps=30, svd_solver="arpack", random_state=SEED)
gc.collect()
step("PCA(30, arpack) done")

pca_use = "X_pca"
harmony_ok = False
try:
    import harmonypy as hm
    try:
        ho = hm.run_harmony(sub.obsm["X_pca"], sub.obs, ["patient_id"], random_state=SEED)
    except TypeError:
        ho = hm.run_harmony(sub.obsm["X_pca"], sub.obs, ["patient_id"])
    sub.obsm["X_pca_harmony"] = ho.Z_corr.T
    pca_use = "X_pca_harmony"
    harmony_ok = True
except Exception as exc:  # pragma: no cover
    log.warning("Harmony failed (%s) -> uncorrected PCA.", exc)
step(f"Harmony on patient_id done (ok={harmony_ok}); use_rep={pca_use}")

sc.pp.neighbors(sub, use_rep=pca_use, n_neighbors=15, n_pcs=30, random_state=SEED)
gc.collect()
step("neighbors(k=15) done")


def _leiden(key: str, res: float) -> None:
    try:
        sc.tl.leiden(sub, resolution=res, key_added=key, random_state=SEED)
    except TypeError:
        sc.tl.leiden(sub, resolution=res, key_added=key)


RES_KEYS = {r: f"leiden_r{r}" for r in RESOLUTIONS}
for r in RESOLUTIONS:
    _leiden(RES_KEYS[r], r)
    step(f"leiden res={r} -> {sub.obs[RES_KEYS[r]].nunique()} clusters")
gc.collect()

# ---------------------------------------------------------------------------
# 4. Annotation / resolution criteria
# ---------------------------------------------------------------------------
names_all = [c.replace("score_", "", 1) for c in score_cols]
markers_present = {t: [g for g in MARKERS[t] if g.upper() in sub.var_names]
                   for t in MARKERS}
log.info("marker genes present in gene space: %s",
         {t: len(v) for t, v in markers_present.items()})


def evaluate_resolution(res: float) -> dict:
    lk = RES_KEYS[res]
    obs = sub.obs
    g = obs.groupby(lk, observed=True)[score_cols].mean()
    sizes = obs[lk].value_counts()
    S = g.to_numpy(dtype=np.float64)
    order = np.argsort(-S, axis=1)

    sc.tl.rank_genes_groups(sub, groupby=lk, method="wilcoxon", use_raw=True,
                            n_genes=sub.n_vars, key_added=f"rgg_{lk}")
    clusters_out = []
    for ci, cl in enumerate(g.index):
        cl = str(cl)
        winner = names_all[order[ci, 0]]
        second = names_all[order[ci, 1]]
        w_s = float(S[ci, order[ci, 0]])
        s_s = float(S[ci, order[ci, 1]])
        margin = w_s - s_s
        ratio = (s_s / w_s) if w_s != 0 else float("inf")
        n_cells_cl = int(sizes.get(cl, 0))
        flags = []
        if n_cells_cl < MIN_CLUSTER_CELLS:
            flags.append("cluster_lt50")
        if margin < MIN_MARGIN:
            flags.append(f"margin_lt_{MIN_MARGIN}")
        if ratio >= MAX_RATIO:
            flags.append("possible_doublet")
        if w_s == 0:
            flags.append("winner_zero")
        if w_s <= 0:
            flags.append("winner_nonpositive")
        df = sc.get.rank_genes_groups_df(sub, group=cl, key=f"rgg_{lk}")
        wm = set(markers_present.get(winner, []))
        dff = df[df["names"].isin(wm)].dropna(subset=["logfoldchanges", "pvals_adj"])
        de_markers = dff[(dff["logfoldchanges"] > DE_LOGFC) &
                         (dff["pvals_adj"] < DE_PADJ)]["names"].tolist()
        if len(de_markers) < MIN_DE_MARKERS:
            flags.append(f"de_markers_lt{MIN_DE_MARKERS}")
        clusters_out.append({
            "cluster": cl, "n_cells": n_cells_cl, "winner": winner,
            "winner_score": round(w_s, 4), "second": second,
            "second_score": round(s_s, 4), "margin": round(margin, 4),
            "ratio": (round(ratio, 4) if np.isfinite(ratio) else None),
            "n_markers_in_var": len(markers_present.get(winner, [])),
            "n_de_up_markers": len(de_markers),
            "de_up_markers": de_markers,
            "passes": len(flags) == 0, "flags": flags,
        })

    per_type = {}
    for t in PRIMARY:
        win = [c for c in clusters_out if c["winner"] == t]
        res_clusters = [c for c in win if c["passes"]]
        cell_mask = obs[lk].isin([c["cluster"] for c in res_clusters]).to_numpy()
        n_cells = int(cell_mask.sum())
        per_sample = (obs.loc[cell_mask, "gsm"].value_counts() if n_cells else
                      pd.Series(dtype=int))
        n_s_ge10 = int((per_sample >= 10).sum()) if n_cells else 0
        de_top = []
        for c in res_clusters:
            de_top.extend(c["de_up_markers"])
        seen = set()
        de_top = [x for x in de_top if not (x in seen or seen.add(x))]
        per_type[t] = {
            "resolved": len(res_clusters) >= 1,
            "winner_clusters": [c["cluster"] for c in win],
            "resolved_clusters": [c["cluster"] for c in res_clusters],
            "n_resolved_clusters": len(res_clusters),
            "n_cells": n_cells,
            "n_samples_with_ge10": n_s_ge10,
            "n_marker_genes_in_var": len(markers_present.get(t, [])),
            "de_up_markers": de_top,
            "flags": sorted({f for c in win for f in c["flags"]}),
            "evaluable": bool(len(res_clusters) >= 1 and n_cells >= MIN_MAPPED_CELLS
                              and n_s_ge10 >= MIN_SAMPLES_GE10),
        }
    return {"clusters": clusters_out, "per_type": per_type}


results: dict[float, dict] = {}
for r in RESOLUTIONS:
    results[r] = evaluate_resolution(r)
    n_ev = sum(v["evaluable"] for v in results[r]["per_type"].values())
    step(f"res={r}: {sum(v['resolved'] for v in results[r]['per_type'].values())}/10 "
         f"resolved, {n_ev}/10 evaluable")

# Contamination check: compartment clusters whose argmax winner is a NON-immune
# label (a MARKERS key outside the K=10 primary immune family).
NON_IMMUNE = [t for t in MARKERS if t not in PRIMARY]
contamination: dict[str, dict] = {}
for r in RESOLUTIONS:
    cl = results[r]["clusters"]
    winners = [c["winner"] for c in cl]
    non = [w for w in winners if w in NON_IMMUNE]
    contamination[str(r)] = {
        "n_clusters": len(cl),
        "n_non_immune_argmax_clusters": len(non),
        "non_immune_winner_counts": (pd.Series(non).value_counts().to_dict()
                                     if non else {}),
        "non_immune_labels": NON_IMMUNE,
    }
    log.info("res=%s contamination: %d/%d compartment clusters argmax a NON-immune "
             "label %s", r, len(non), len(cl),
             contamination[str(r)]["non_immune_winner_counts"])

# ---------------------------------------------------------------------------
# 5. Map back to cells; per-sample proportions on 221,016 denominators
# ---------------------------------------------------------------------------
# Denominator = each sample's TOTAL cap=3000 cell count over the full 221,016
# cells (captured at load time), NOT the primary proportions CSV (whose rows sum
# to ~1.0 and are therefore not counts).  The CSV is read read-only, solely for
# its sample ordering.
props_index = pd.read_csv(OUT / "sn_proportions_per_sample.csv", index_col=0).index
TOTALS = SAMPLE_TOTALS.reindex(props_index).astype(int)
TOTALS.index.name = "sample_id"
log.info("denominator = per-sample TOTAL cap=3000 cells: n_samples=%d, sum=%d "
         "(expect 221016)", len(TOTALS), int(TOTALS.sum()))

lk = RES_KEYS[PRIMARY_RES]
argmax_map = (sub.obs.groupby(lk, observed=True)[score_cols].mean()
              .idxmax(axis=1).str.replace("score_", "", regex=False).to_dict())
assigned = sub.obs[lk].map(argmax_map)
counts = (pd.DataFrame({"sample_id": sub.obs["gsm"].to_numpy(),
                        "type": assigned.to_numpy()})
          .groupby(["sample_id", "type"], observed=True).size().unstack(fill_value=0))
counts = counts.reindex(columns=list(MARKERS.keys()), fill_value=0)
counts = counts.reindex(TOTALS.index, fill_value=0)
props = counts.div(TOTALS, axis=0)
props.index.name = "sample_id"
props_path = OUT / f"sn_proportions_per_sample{SUFFIX}.csv"
props.to_csv(props_path)
step(f"wrote per-sample proportions {props.shape} -> {props_path.name}")

# ---------------------------------------------------------------------------
# 6. Evaluability table + audit
# ---------------------------------------------------------------------------
rows = []
for r in RESOLUTIONS:
    for t in PRIMARY:
        d = results[r]["per_type"][t]
        rows.append({
            "cell_type": t, "resolution": r, "resolved": bool(d["resolved"]),
            "n_resolved_clusters": d["n_resolved_clusters"],
            "resolved_clusters": ";".join(d["resolved_clusters"]),
            "n_cells": d["n_cells"],
            "n_samples_with_ge10": d["n_samples_with_ge10"],
            "n_marker_genes_in_var": d["n_marker_genes_in_var"],
            "top_de_markers": ";".join(d["de_up_markers"]),
            "flags": ";".join(d["flags"]),
            "evaluable": bool(d["evaluable"]),
        })
eval_df = pd.DataFrame(rows)
eval_path = OUT / f"immune_subclust_evaluability{SUFFIX}.csv"
eval_df.to_csv(eval_path, index=False)
step(f"wrote evaluability table {eval_df.shape} -> {eval_path.name}")

n_eval = {r: int(sum(v["evaluable"] for v in results[r]["per_type"].values()))
          for r in RESOLUTIONS}
gate_pass = bool(n_eval[1.0] >= 6 and n_eval[0.5] >= 5 and n_eval[1.5] >= 5)
free_end = free_disk_gb()

audit = {
    "label": "SENSITIVITY: annotation deviation from pre-spec v1.0 "
             "(immune-compartment subclustering, cap=3000 input unchanged)",
    "ran_at": datetime.now().isoformat(timespec="seconds"),
    "script": "scripts/immune_subclust_sensitivity.py",
    "input_path": str(CKPT), "input_sha256": INPUT_SHA,
    "input_bytes": CKPT.stat().st_size,
    "repo_git_head": git_head(),
    "compartment": {
        "rule": "union((a) leiden_r0.8 clusters labelled Microglia/MoMF/NK/B/CD8_T, "
                "(b) PTPRC log-norm > 0 in any cluster)",
        "rule_a_cells": int(mask_a.sum()),
        "rule_b_cells": int(mask_b.sum()),
        "overlap_cells": int((mask_a & mask_b).sum()),
        "ptprc_present": bool(ptprc_present),
        "rule_b_skipped": (not bool(ptprc_present)),
        "n_cells": n_comp, "fraction_of_all_cells": round(comp_fraction, 4),
        "seed_clusters": [str(c) for c in seed_clusters],
        "primary_label_composition": comp_labels,
    },
    "predeclared": {
        "seed": SEED, "resolutions": list(RESOLUTIONS),
        "primary_resolution": PRIMARY_RES, "min_cluster_cells": MIN_CLUSTER_CELLS,
        "min_margin": MIN_MARGIN, "max_second_over_winner_ratio": MAX_RATIO,
        "min_de_markers": MIN_DE_MARKERS, "de_logfc": DE_LOGFC, "de_padj": DE_PADJ,
        "min_mapped_cells": MIN_MAPPED_CELLS, "min_samples_ge10": MIN_SAMPLES_GE10,
        "harmony_applied": harmony_ok, "use_rep": pca_use,
        "regress_out": "NOT applied (not part of the pre-declared design)",
    },
    "marker_genes_present": {t: len(v) for t, v in markers_present.items()},
    "n_cells_total_cap3000": int(TOTALS.sum()),
    "n_samples_total_cap3000": int(len(TOTALS)),
    "evaluable_counts": {str(r): n_eval[r] for r in RESOLUTIONS},
    "gate": {"rule": ">=6/10 at 1.0 AND >=5/10 at both 0.5 and 1.5",
             "pass": gate_pass,
             "decision": ("ASK_USER_BEFORE_ANYTHING_ELSE" if gate_pass
                          else "STOP_subclustering_does_not_recover_family")},
    "clusters_per_resolution": {str(r): results[r]["clusters"] for r in RESOLUTIONS},
    "contamination_check": contamination,
    "notes": [
        "SENSITIVITY: annotation deviation from pre-spec v1.0 - NOT a primary result.",
        "All columns of sn_proportions_per_sample_immune_subclust.csv are the full "
        "MARKERS key names (15 columns, no aliases).",
        "Per-sample proportions use each sample's TOTAL cap=3000 cell count over the "
        "full 221,016 cells; cells outside the immune compartment belong to no type, "
        "so rows do NOT sum to 1.",
        "regress_out deliberately NOT applied here (the main pipeline applied it).",
        f"Leiden random_state={SEED} here vs scanpy default 0 in the main pipeline.",
        "score_genes and rank_genes_groups use_raw=True on the unscaled log-norm matrix.",
    ],
    "resource": {
        "seed": SEED, "free_disk_gb_before": round(free_start, 2),
        "free_disk_gb_after": round(free_end, 2),
        "peak_rss_gb": round(_rss_peak[0], 2),
        "runtime_s": round(time.time() - T0, 1),
        "subset_shape": list(sub.shape),
    },
    "outputs": [props_path.name, eval_path.name],
}
audit_path = OUT / f"immune_subclust_audit{SUFFIX}.json"
with open(audit_path, "w") as fh:
    json.dump(audit, fh, indent=2, default=str)
step(f"wrote audit -> {audit_path.name}")

# ---------------------------------------------------------------------------
# 7. Report
# ---------------------------------------------------------------------------
log.info("=" * 78)
log.info("GATE TABLE (evaluable = resolved AND >=50 cells AND >=10 samples with >=10 cells)")
hdr = f"{'primary cell type':38s} " + " ".join(f"res{r:<4}" for r in RESOLUTIONS)
log.info(hdr)
for t in PRIMARY:
    cells = " ".join(
        f"{int(results[r]['per_type'][t]['n_cells']):>4d}"
        f"{'E' if results[r]['per_type'][t]['evaluable'] else ('R' if results[r]['per_type'][t]['resolved'] else '.')}"
        for r in RESOLUTIONS)
    log.info("%-38s %s", t, cells)
log.info("evaluable: " + ", ".join(f"res{r}={n_eval[r]}/10" for r in RESOLUTIONS))
log.info("GATE %s -> %s", "PASS" if gate_pass else "FAIL", audit["gate"]["decision"])
log.info("peak rss %.2f GB | runtime %.1f s | free disk %.2f -> %.2f GB",
         _rss_peak[0], time.time() - T0, free_start, free_end)
log.info("DONE immune_subclust_sensitivity.py (no h5ad written).")
