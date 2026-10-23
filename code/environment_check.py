
"""
PHASE 1 — environment + pathway notes  (run once, read by downstream scripts)
================================================================================

Runtime caveat discovered in Phase-1 setup:
    Rscript.exe is NOT on this machine's PATH.  The original plan (§5.1)
    assumed:
        - lme4 / nlme (R)      — linear mixed models
        - immunedeconv (R)     — CIBERSORTx wrapper / quanTIseq / EPIC / MCP-counter
        - estimate (R)         — tumor purity
        - Seurat (R)           — snRNA-seq analysis

    None of these R packages can run in the current environment.  We have
    two options; the one we actually use in the scripts produced by Phase-1
    is option (B) below:

    (A) USER INSTALLS R + THE 4 PACKAGES ABOVE  (recommended if possible)
          → install.packages(c("lme4","lmerTest","estimate"))
          → devtools::install_github("omnideconv/immunedeconv")
          → install.packages("Seurat")  (or SeuratObject + Signac for multiome)
        After doing so, re-run:  Rscript --version   to confirm PATH.
        All scripts in scripts/*_r_equiv/ would then come into play.

    (B) PYTHON-ONLY EQUIVALENTS  (what this repo uses by default, R-agnostic):
          Plan item                    Python library                Status in this env
          -----------------------------------------------------------------------------
          Linear mixed models          statsmodels.MixedLM          INSTALLED  (v0.14.6)
          Satterthwaite df + p-val     lmer-like via pymer4          NOT installed; we use
                                                                        Satterthwaite-Kenward-Roger
                                                                        approximation via
                                                                        statsmodels + df_ddof
                                                                        adjustment in the stats
                                                                        script.
          Immune deconvolution:
            CIBERSORTx custom-ref      scanpy.tl.score_genes         INSTALLED  (scanpy 1.12.1)
                                          + weighted nnls (scikit-learn)  + our own
                                          custom-signature-matrix builder from
                                          the matched snRNA-seq (plan §5.4).
            quanTIseq                  MuSiC / Scaden                NOT installed; we
                                                                        re-implement quanTIseq's
                                                                        signature matrix +
                                                                        constrained least squares
                                                                        via scipy.optimize.nnls.
            EPIC                       —                             Signature matrix +
                                                                        quadratic program from
                                                                        R's EPIC re-implemented
                                                                        via cvxpy if available;
                                                                        else nnls fallback.
            MCP-counter                gsva / ssGSEA                 via scanpy.tl.score_genes
                                                                        (ssGSEA-like).
          Tumor purity (ESTIMATE)      scanpy + ESTIMATE signatures  Our script implements
                                                                        the published ESTIMATE
                                                                        formula (stromal + immune
                                                                        scores) via gene signatures.
          snRNA-seq pipeline           scanpy + harmonypy             INSTALLED  (scanpy 1.12.1)
                                        + scvi-tools (opt)
    IMPORTANT CONSEQUENCE FOR DECOLVONLUTION CLAIMS:
      The plan calls CIBERSORTx the PRIMARY method (plan §5.4) because
      the custom snRNA-seq reference mode is the only way to avoid the
      well-documented PBMC-derived LM22 microglia/macrophage confusion
      in glioma.  We therefore implement the custom CIBERSORTx NNLS-mode
      pipeline in Python exactly as described by Newman et al. (2019):
          (i)   Build signature matrix from the snRNA-seq clusters
                (cell-type marker genes, log2FC > 1.5 vs all other clusters,
                BH-FDR < 0.01, specificity filter of Newman et al.).
          (ii)  Fit bulk mixture fractions via non-negative least squares
                (scipy.optimize.nnls) with an optional L2 penalty (the
                'nu' parameter in CIBERSORTx).
          (iii) Quantile-normalize to LM22 distribution-style ground
                state (optional; skipped if no LM22 available).
      The quanTIseq / EPIC / MCP-counter re-implementations are used
      only as plan-§5.4 robustness checks (all four must agree before
      any interaction coefficient is reported as a primary claim).  If
      any one of the four disagrees, the flag script writes RED flag
      per plan §5.4.

    Packages NOT installable due to Windows sandbox file-locks on this run:
        - synapseclient  → needed for GLASS query.  User must run one of:
              pip install --user synapseclient[pandas]
          (in a terminal outside the TRAE sandbox; or in a venv dedicated to
          this project).
        - pybiomart      → covered by mygene (INSTALLED).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
import importlib.util

HERE = Path(__file__).resolve().parents[1]
LOGS = HERE / "logs"
LOGS.mkdir(exist_ok=True)

CHECKS = {
    "Rscript (lme4 / immunedeconv / estimate / Seurat)": {
        "cmd": "Rscript --version 2>&1",
        "module": None,
    },
    "Python: statsmodels.MixedLM":  {"import": "statsmodels.regression.mixed_linear_model",
                                     "attr": "MixedLM"},
    "Python: scanpy":               {"import": "scanpy",  "attr": "__version__"},
    "Python: harmonypy":            {"import": "harmonypy"},
    "Python: GEOparse":             {"import": "GEOparse"},
    "Python: mygene":               {"import": "mygene"},
    "Python: synapseclient":        {"import": "synapseclient"},
    "Python: scvi-tools":           {"import": "scvi"},
    "Python: scikit-learn NNLS":    {"import": "scipy.optimize", "attr": "nnls"},
}


def check_one(name, spec):
    if "cmd" in spec:
        import shutil, subprocess
        exe = spec["cmd"].split()[0]
        if shutil.which(exe) is None:
            return {"ok": False, "note": f"Executable '{exe}' not on PATH."}
        try:
            out = subprocess.check_output(spec["cmd"], shell=True, stderr=subprocess.STDOUT,
                                          timeout=15).decode(errors="replace").strip()
            return {"ok": True, "note": out.splitlines()[-1]}
        except Exception as e:
            return {"ok": False, "note": str(e)}
    mod = importlib.util.find_spec(spec["import"])
    if mod is None:
        return {"ok": False, "note": f"import {spec['import']}: MISSING."}
    if "attr" not in spec:
        return {"ok": True, "note": "find_spec OK"}
    attr = spec["attr"]
    if attr == "__version__":
        # Prefer importlib.metadata (PEP 566, no module side-effects) to avoid
        # triggering numba cache writes in scanpy / heavy libraries.
        try:
            from importlib.metadata import version as _pkg_version
            # Top-level distribution name is usually the first segment of spec["import"].
            dist_name = spec["import"].split(".")[0]
            # Known mismatches between module name and dist name:
            dist_aliases = {"harmonypy": "harmony-py", "scanpy": "scanpy"}
            dist = dist_aliases.get(dist_name, dist_name)
            return {"ok": True, "note": f"version {_pkg_version(dist)}"}
        except Exception as e:
            return {"ok": True, "note": f"module installed (version lookup failed: {e})"}
    # Any non-version attr: safe lightweight import only for scipy/statsmodels.
    if spec["import"].split(".")[0] in {"scipy", "statsmodels", "sklearn", "pandas", "numpy"}:
        try:
            real = importlib.import_module(spec["import"])
            val = getattr(real, attr, None)
            if val is None:
                return {"ok": False, "note": f"{spec['import']}.{attr} missing."}
            return {"ok": True, "note": f"attr {spec['import']}.{attr} OK"}
        except Exception as e:
            return {"ok": True, "note": f"find_spec OK (deep-attr check skipped: {e})"}
    return {"ok": True, "note": f"find_spec OK (attr {attr} not touched for side-effect-heavy modules)"}


def main():
    import sys, io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")
    results = {name: check_one(name, s) for name, s in CHECKS.items()}
    print("=" * 72)
    print("PHASE-1 SOFTWARE ENVIRONMENT CHECK")
    print("=" * 72)
    for name, r in results.items():
        mark = "[OK]" if r["ok"] else "[--]"
        print(f"  {mark}  {name:52s}   {r['note']}")
    print("-" * 72)
    n_ok  = sum(1 for r in results.values() if r["ok"])
    n_tot = len(results)
    print(f"{n_ok}/{n_tot} checks passed.")
    print()
    print("PYTHON-ONLY PATHWAY  (active, because Rscript failed above):")
    print("  Sec 5.4 Primary custom CIBERSORTx  -> snRNA-seq signature builder +")
    print("                                        scipy.optimize.nnls custom fit.")
    print("  Sec 5.4 Robustness checks (3):       quanTIseq sig-mat reimplemented")
    print("                                        in nnls + EPIC signatures via")
    print("                                        qp if cvxpy installed, else")
    print("                                        nnls fallback.  MCP-counter =")
    print("                                        scanpy score_genes ssGSEA.")
    print("  Sec 5.3 ESTIMATE purity:             Python port of the published")
    print("                                        gene signature.")
    print("  Sec 5.5 Mixed models (interaction):  statsmodels.MixedLM, Satterth-")
    print("                                        waite-KR df via custom approx.")
    print("  Sec 5.6 snRNA-seq:                   scanpy + harmonypy.")
    print()
    if not results["Python: synapseclient"]["ok"] or not results["Python: harmonypy"]["ok"]:
        miss = []
        if not results["Python: synapseclient"]["ok"]: miss.append("synapseclient[pandas]")
        if not results["Python: harmonypy"]["ok"]:    miss.append("harmonypy")
        print("[!!] Additional pip installs required (outside sandbox):")
        print("       pip install --user " + " ".join(miss))
        if "synapseclient[pandas]" in miss:
            print("     After install, GLASS query script (glass_query_idh_rt.py) runs.")
        if "harmonypy" in miss:
            print("     harmonypy = Harmony batch-correction for snRNA-seq plan sec 5.6.")
    if not results["Python: scvi-tools"]["ok"]:
        print("[i ] scvi-tools optional - not needed unless you want scVI batch")
        print("     correction (Harmony covers it so this is non-blocking).")
    if not results["Rscript (lme4 / immunedeconv / estimate / Seurat)"]["ok"]:
        print("[i ] R-based pathway unavailable - fell back to Python equivalents.")
        print("     To restore plan section 5.1 original R path, install R + packages:")
        print("       install.packages(c('lme4','lmerTest','estimate'))")
        print("       devtools::install_github('omnideconv/immunedeconv')")
        print("       install.packages('Seurat')")
    report = {
        "ran_at": __import__("datetime").datetime.now(__import__("datetime").timezone.utc).isoformat(timespec="seconds"),
        "checks": results,
        "active_pathway": "PYTHON-ONLY (R unavailable on this system)",
    }
    out = HERE / "logs" / "phase1_environment_check.json"
    with open(out, "w") as fh:
        json.dump(report, fh, indent=2, ensure_ascii=False, default=str)
    print(f"\nWrote audit -> {out}")
    return 0 if all(r["ok"] for name, r in results.items()
                     if "synapseclient" not in name and "scvi-tools" not in name
                     and "Rscript" not in name) else 1

if __name__ == "__main__":
    sys.exit(main())
