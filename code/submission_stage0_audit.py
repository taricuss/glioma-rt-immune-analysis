# -*- coding: utf-8 -*-
"""Stage 0: verify project state for the submission package build.

Read-only with respect to the project:
  * records modification times of every file under outputs/, data/, logs/
    plus the locked pre-spec files (used again in Stage 6 to prove nothing changed)
  * copies the three existing manuscript artifacts to %TEMP% (backup, per task rule)
  * prints the timestamp table requested by Stage 0
No pipeline is imported or executed.
"""
import hashlib
import io
import json
import os
import shutil
import sys
import time

ROOT = r"<PROJECT_ROOT>"
TMP = os.environ.get("TEMP", r"<PROJECT_ROOT>")
os.chdir(ROOT)

WATCH_FILES = [
    r"outputs\intermediate\lmm_interaction_summary.csv",
    r"outputs\intermediate\snrna_direct_perm_pvalues.csv",
    r"outputs\intermediate\lmm_interaction_summary_glass.csv",
    r"outputs\intermediate\glass_direction_consistency.json",
    r"outputs\intermediate\sn_proportions_per_sample.csv",
    r"outputs\intermediate\immune_subclust_evaluability_immune_subclust.csv",
    r"outputs\intermediate\immune_subclust_audit_immune_subclust.json",
    r"outputs\intermediate\post_data_power_table.csv",
    r"outputs\intermediate\post_data_power_audit.json",
    r"outputs\intermediate\lmm_summary.json",
    r"outputs\intermediate\lmm_summary_glass.json",
    r"outputs\intermediate\cgga_crosssectional_stats.csv",
    r"outputs\intermediate\cgga_pipeline_audit.json",
    r"outputs\intermediate\clinical_canonical_adapter_audit.json",
    r"outputs\intermediate\glass_pipeline_audit.json",
    r"outputs\intermediate\sn_pipeline_audit.json",
    r"outputs\sec4_prespec.json",
    r"outputs\sec4_prespec.md",
    r"phase0_summary.md",
    r"plan.md",
    r"plan_implementation_status.md",
]

BACKUP_FILES = [
    r"outputs\manuscript_skeleton.html",
    r"outputs\manuscript_audit.json",
    r"scripts\build_manuscript_skeleton.py",
]

WATCH_DIRS = [r"outputs", r"data", r"logs"]


def sha256(path, limit=None):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(1024 * 1024)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def walk(root):
    out = []
    for dp, dn, fn in os.walk(root):
        dn[:] = [d for d in dn if d != "__pycache__"]
        for name in fn:
            out.append(os.path.join(dp, name))
    return out


def main():
    print("=" * 78)
    print("STAGE 0 - project state audit")
    print("=" * 78)
    print("python:", sys.version.split()[0])
    print("cwd   :", os.getcwd())
    print()

    # ---- 1. requested timestamp table -------------------------------------
    print("-" * 78)
    print("Timestamp table (requested files)")
    print("-" * 78)
    missing = []
    rows = []
    for rel in WATCH_FILES:
        p = os.path.join(ROOT, rel)
        if not os.path.exists(p):
            missing.append(rel)
            print("MISSING  %-62s" % rel)
            rows.append({"path": rel, "exists": False})
            continue
        st = os.stat(p)
        ts = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(st.st_mtime))
        print("%s  %10d bytes  %s" % (ts, st.st_size, rel))
        rows.append({
            "path": rel, "exists": True, "mtime": st.st_mtime,
            "mtime_str": ts, "size": st.st_size,
        })
    print()

    # ---- 2. full baseline of watched dirs + pre-spec files ----------------
    baseline = {"files": {}, "generated": time.strftime("%Y-%m-%d %H:%M:%S")}
    n = 0
    for d in WATCH_DIRS:
        for p in walk(os.path.join(ROOT, d)):
            rel = os.path.relpath(p, ROOT)
            st = os.stat(p)
            baseline["files"][rel] = {
                "mtime": st.st_mtime, "size": st.st_size, "sha256": sha256(p),
            }
            n += 1
    for rel in (r"outputs\sec4_prespec.json", r"outputs\sec4_prespec.md"):
        p = os.path.join(ROOT, rel)
        st = os.stat(p)
        baseline["files"][rel] = {
            "mtime": st.st_mtime, "size": st.st_size, "sha256": sha256(p),
        }
        n += 1

    base_path = os.path.join(TMP, "dgnc_stage0_baseline.json")
    with io.open(base_path, "w", encoding="utf-8") as f:
        json.dump(baseline, f, indent=1)
    print("baseline files recorded:", n, "->", base_path)
    print()

    # ---- 3. %TEMP% backups of the existing manuscript artifacts -----------
    print("-" * 78)
    print("%TEMP% backups")
    print("-" * 78)
    for rel in BACKUP_FILES:
        src = os.path.join(ROOT, rel)
        dst = os.path.join(TMP, "dgnc_bak_" + os.path.basename(rel))
        if not os.path.exists(src):
            print("MISSING  %-52s" % rel)
            continue
        shutil.copy2(src, dst)
        print("backed up %-52s -> %s (%d bytes)"
              % (rel, os.path.basename(dst), os.path.getsize(dst)))
    print()

    if missing:
        print("!! missing requested files:", missing)
    else:
        print("All requested files present.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
