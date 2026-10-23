# -*- coding: utf-8 -*-
"""Stage 5 - assemble the submission package.

Copies/redacts only aggregate artefacts, keeps the locked pre-specification
byte-identical, writes the README, the congress abstract and the TODO list, and
writes checks/manifest.sha256 over every file in the package.
No analysis is run and no existing project file is modified.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import os
import re
import shutil

from docx import Document
from docx.shared import Pt

ROOT = r"<PROJECT_ROOT>"
SUB = os.path.join(ROOT, "submission")
INT = os.path.join(ROOT, "outputs", "intermediate")
SCRIPTS = os.path.join(ROOT, "scripts")
SUP = os.path.join(SUB, "supplementary")
AUD = os.path.join(SUP, "audits")
CODE = os.path.join(SUB, "code")
CHK = os.path.join(SUB, "checks")
LST = os.path.join(SUB, "checklists")
for d in (SUP, AUD, CODE, CHK, LST):
    os.makedirs(d, exist_ok=True)

PACKAGE_DATE = "2026-10-20"


# ------------------------------------------------------------------ registry
def _nest(flat):
    out = {}
    for k, v in flat.items():
        cur = out
        parts = k.split(".")
        for p in parts[:-1]:
            cur = cur.setdefault(p, {})
        cur[parts[-1]] = v
    return out


NUM = _nest(json.load(io.open(os.path.join(SUP, "results_numbers.json"),
                              "r", encoding="utf-8"))["values"])


def g(dotted):
    parts = dotted.split(".")
    for i in range(len(parts), 0, -1):
        prefix = ".".join(parts[:i])
        if prefix in NUM:
            cur, rest = NUM[prefix], parts[i:]
            while rest:
                for j in range(len(rest), 0, -1):
                    key = ".".join(rest[:j])
                    if isinstance(cur, dict) and key in cur:
                        cur, rest = cur[key], rest[j:]
                        break
                else:
                    raise KeyError(dotted)
            return cur
    raise KeyError(dotted)


def n0(d):
    return "%d" % int(round(float(g(d))))


def n3(d):
    return "%.3f" % float(g(d))


def pct(d):
    return "%.1f%%" % (100.0 * float(g(d)))


def python_version():
    """The interpreter version, taken from Table S5 (never hand-typed)."""
    with io.open(os.path.join(SUB, "tables", "TableS5.csv"), "r", encoding="utf-8",
                 newline="") as f:
        for r in csv.reader(f):
            if len(r) >= 2 and r[0].lower() == "python":
                return r[1]
    raise RuntimeError("no python row in TableS5.csv")


PYVER = python_version()


# ------------------------------------------------------------------ redaction
PATH_RE = re.compile(r"[A-Za-z]:\\\\?[^\s\"']+|/[Uu]sers/[^\s\"']+")
ID_RE = re.compile(r"\bGSM\d{6,}\b|\b(?:SJ|MDA)\d{3,}\b|\bP(?:1[0-9]|[2-9][0-9])\b(?![A-Za-z0-9])")
DROP_KEYS = {"excluded_patients", "skipped_samples", "sample_ids", "patients",
             "per_sample", "samples", "gsm_ids"}


def redact(obj, dropped):
    """Recursively strip local paths and patient/sample identifiers."""
    if isinstance(obj, dict):
        out = {}
        for k, v in obj.items():
            if k in DROP_KEYS and isinstance(v, (list, dict)):
                dropped.add(k)
                n = len(v)
                out[k + "_n"] = n
                continue
            out[k] = redact(v, dropped)
        return out
    if isinstance(obj, list):
        return [redact(v, dropped) for v in obj]
    if isinstance(obj, str):
        s = PATH_RE.sub("[PATH REMOVED]", obj)
        s = ID_RE.sub("[ID REMOVED]", s)
        return s
    return obj


AUDIT_FILES = [
    "lmm_summary.json",
    "lmm_summary_glass.json",
    "sn_pipeline_audit.json",
    "glass_pipeline_audit.json",
    "cgga_pipeline_audit.json",
    "cgga_crosssectional_stats.csv",
    "post_data_power_audit.json",
    "clinical_canonical_adapter_audit.json",
    "immune_subclust_audit_immune_subclust.json",
    "immune_subclust_evaluability_immune_subclust.csv",
]


def copy_audits():
    kept, skipped = [], []
    for name in AUDIT_FILES:
        src = os.path.join(INT, name)
        if not os.path.exists(src):
            skipped.append(name + " (absent)")
            continue
        with io.open(src, "r", encoding="utf-8-sig", newline="") as f:
            raw = f.read()
        dropped = set()
        if name.endswith(".json"):
            obj = redact(json.loads(raw), dropped)
            txt = json.dumps(obj, indent=2, ensure_ascii=False)
        else:
            rows = list(csv.reader(io.StringIO(raw)))
            rows = [[redact(c, dropped) for c in r] for r in rows]
            buf = io.StringIO()
            csv.writer(buf, lineterminator="\n").writerows(rows)
            txt = buf.getvalue()
        if ID_RE.search(txt) or PATH_RE.search(txt):
            skipped.append(name + " (identifiers or paths could not be removed cleanly)")
            continue
        with io.open(os.path.join(AUD, name), "w", encoding="utf-8",
                     newline="\n") as f:
            f.write(txt if txt.endswith("\n") else txt + "\n")
        kept.append(name + (" [keys reduced to counts: %s]"
                            % ", ".join(sorted(dropped)) if dropped else ""))
    return kept, skipped


# ------------------------------------------------------------------ other files
def copy_prespec():
    out = []
    for name in ("sec4_prespec.json", "sec4_prespec.md"):
        src = os.path.join(ROOT, "outputs", name)
        dst = os.path.join(SUP, name)
        shutil.copyfile(src, dst)          # byte-identical, source untouched
        out.append(name)
    return out


def copy_code():
    n = 0
    for name in sorted(os.listdir(SCRIPTS)):
        if not name.endswith(".py") or name.startswith("_"):
            continue
        with io.open(os.path.join(SCRIPTS, name), "r", encoding="utf-8") as f:
            src = f.read()
        src = PATH_RE.sub("<PROJECT_ROOT>", src)
        src = re.sub(r"[USER]", "[USER]", src, flags=re.IGNORECASE)
        with io.open(os.path.join(CODE, name), "w", encoding="utf-8",
                     newline="\n") as f:
            f.write(src)
        n += 1
    pkgs = []
    with io.open(os.path.join(SUB, "tables", "TableS5.csv"), "r", encoding="utf-8",
                 newline="") as f:
        for r in list(csv.reader(f))[1:]:
            if len(r) >= 2 and r[0] not in ("note", "python"):
                pkgs.append("%s==%s" % (r[0], r[1]))
    with io.open(os.path.join(CODE, "requirements.txt"), "w", encoding="utf-8",
                 newline="\n") as f:
        f.write("\n".join(pkgs) + "\n")
    with io.open(os.path.join(CODE, "README.md"), "w", encoding="utf-8",
                 newline="\n") as f:
        f.write(
            "# Code\n\n"
            "Copies of the analysis, table and figure scripts, excluded files whose "
            "names begin with an underscore (internal diagnostics) and the virtual "
            "environment. `requirements.txt` is the frozen dependency list of the "
            "analysis environment (Python %s).\n\n"
            "The scripts are the analysis code as run, except that local absolute "
            "paths were replaced by `<PROJECT_ROOT>` and the machine account name by "
            "`[USER]`. They expect the project's local directory layout. They read "
            "aggregate intermediate artefacts from `outputs/`; those artefacts are "
            "not redistributed here. No script in this folder is executed by the "
            "package build, and no pipeline was re-run to produce any number in the "
            "manuscript.\n\n"
            "Reproducing the tables and figures from the aggregate files requires "
            "only `submission/tables/*.csv`, `submission/figures/figure_legends.txt` "
            "and `submission/supplementary/results_numbers.json`; see the package "
            "README.\n" % PYVER)
    return n, len(pkgs)


# ------------------------------------------------------------------ readme
README = """# DGNC submission package - {date}

Pre-specified negative-discovery analysis: does radiotherapy exposure moderate the
change in immune cell composition between primary and recurrent IDH-mutant glioma?

**Framing (fixed).** The discovery single-nucleus arm can evaluate {ev_n} of the
{ev_tot} pre-specified primary cell types (NK cells and microglia), and both are
null with a design that is underpowered at the pre-specified anchor effect sizes.
The GLASS bulk arm and the CGGA arm are EXPLORATORY. The pre-declared
immune-subclustering analysis is a SENSITIVITY analysis that failed its gate. This
package claims no replication, no agreement between cohorts and no causal effect of
radiotherapy.

## Contents

| Path | What it is |
| --- | --- |
| `manuscript/manuscript.docx` | Manuscript (Python-docx), figures and Tables 1-6 embedded |
| `manuscript/manuscript.html` | Same content, self-contained CSS, figures as relative images |
| `manuscript/manuscript_metrics.json` | Title options, title used, abstract and main-text word counts |
| `figures/Figure1..4.{{png,tiff,pdf,svg}}` | Main figures, 600 dpi raster, LZW TIFF, vector PDF/SVG |
| `figures/FigureS1..S3.{{png,tiff,pdf,svg}}` | Supplementary figures |
| `figures/figure_legends.txt` / `.docx` | Full legends and alt text for every figure |
| `tables/Table1..6.{{csv,xlsx}}` | Main tables (xlsx: bold header, set column widths, 3 decimals) |
| `tables/TableS1..S6.{{csv,xlsx}}` | Supplementary tables |
| `tables/tables_index.csv` | Table id, title and row count |
| `supplementary/results_numbers.json` | Single source of truth: every reported number |
| `supplementary/numbers_provenance.csv` | Value, source file and source key for each number |
| `supplementary/sec4_prespec.json` / `.md` | Locked pre-specification v{psver}, copied byte-identical |
| `supplementary/audits/*` | Analysis audit artefacts, aggregate content only, paths and identifiers removed |
| `code/*.py`, `code/requirements.txt` | Analysis/table/figure scripts and frozen dependencies |
| `checks/qa_report.md` | Quality-assurance report for this package |
| `checks/manifest.sha256` | sha256 of every file in the package |
| `checklists/DGNC_abstract.docx` | One-page congress abstract (same numbers) |
| `checklists/submission_todo.md` | Everything that must be filled in by hand before submission |

## How the numbers were produced

No analysis was re-run. `supplementary/results_numbers.json` was derived from the
existing aggregate files in `outputs/` and `outputs/intermediate/`; the file and key
of each value are recorded in `numbers_provenance.csv`. Figures and tables read that
JSON (tables also read each other's CSV), and the manuscript reads the tables' CSV
files, so the manuscript, the tables and the figures cannot disagree.

## Reproducing the tables and figures

1. `submission/tables/*.csv` are plain CSV and stand alone.
2. Figures are plotted from `submission/supplementary/results_numbers.json` plus the
   CSVs in `submission/tables/`; legends live in `figures/figure_legends.txt`.
3. `code/submission_build_numbers.py` rebuilds the registry, `code/submission_make_figures.py`
   rebuilds the figures, `code/submission_make_tables.py` rebuilds the tables and
   `code/submission_build_manuscript.py` rebuilds the manuscript. These scripts expect
   the project's local directory layout.

## Software and seeds

Python {pyver}; the full frozen dependency list is `code/requirements.txt` and
TableS5. Random seeds: {seed} for the patient-level permutation null and for the
post-data power simulation (800 simulation draws); the per-sample cell-cap
subsampling seed is recorded in the analysis scripts.

## Data

No raw or patient-level data are redistributed. Discovery single-nucleus RNA
sequencing: GEO GSE326221, with an additional controlled-access bulk arm (EGA
EGAD50000002485) that was not accessed. External arm: GLASS, Synapse project
syn17038081, subject to its own data-use terms. Cross-sectional arm: CGGA
mRNAseq_693 (20200506), subject to its own data-use terms.

## Local files deliberately not included

Raw and intermediate matrices (`.h5ad`, `.mtx`, tar archives), per-sample and
per-patient tables (`clinical_canonical.csv`, per-sample cell-type proportions,
GLASS and CGGA patient-level outputs), run logs, temporary diagnostics
(`_tmp_*`), the earlier draft manuscript skeleton and its audit trail, and the
project planning/status documents. Files that could not be redacted cleanly are
excluded and listed in `checks/qa_report.md`.

## Deviations from the pre-specification

Listed in Table 6 and summarised in the manuscript Methods. The locked
pre-specification itself is included unmodified in `supplementary/`.
"""


def write_readme(kept, skipped):
    txt = README.format(
        date=PACKAGE_DATE,
        ev_n=n0("evaluability.n_primary_evaluable"),
        ev_tot=n0("evaluability.n_primary_types"),
        psver=str(g("prespec.version")),
        pyver=PYVER,
        seed=n0("power.seed"))
    txt += ("\n## Audit artefacts included\n\n"
            + "\n".join("- " + k for k in kept)
            + "\n\n## Audit artefacts excluded\n\n"
            + ("\n".join("- " + s for s in skipped) if skipped else "- none")
            + "\n")
    with io.open(os.path.join(SUB, "README.md"), "w", encoding="utf-8",
                 newline="\n") as f:
        f.write(txt)
    return os.path.join(SUB, "README.md")


# ------------------------------------------------------------------ checklist
def write_abstract_docx():
    doc = Document()
    st = doc.styles["Normal"]
    st.font.name = "Calibri"
    st.font.size = Pt(11)
    p = doc.add_paragraph()
    r = p.add_run("DGNC 2027 - abstract draft (working copy)")
    r.bold = True
    r.font.size = Pt(14)
    doc.add_paragraph("Title: see manuscript/manuscript_metrics.json (title option 1 "
                      "is used throughout the package).")
    doc.add_paragraph("Authors: [AUTHOR LIST]")
    doc.add_paragraph("Affiliations: [AFFILIATIONS]")

    ev_n, ev_tot = n0("evaluability.n_primary_evaluable"), n0("evaluability.n_primary_types")
    body = [
        ("Background:", "Radiotherapy exposure is a candidate moderator of the change "
         "in immune cell composition between primary and recurrent IDH-mutant glioma, "
         "and a pre-specified test of that interaction has not been reported."),
        ("Methods:", "In a pre-specified analysis of paired single-nucleus RNA "
         "sequencing (%s samples, %s patients; %s with any radiotherapy exposure, %s "
         "with none) we fitted patient-random-intercept mixed models of "
         "arcsin-square-root-transformed cell-type fractions with a timepoint by "
         "treatment-arm interaction, under a pre-declared rule that a primary cell "
         "type must have a non-zero proportion in at least %s of %s samples. We added "
         "a 10,000-shuffle patient-level permutation null (seed %s), a post-data power "
         "calculation (%s simulation draws), two external cohorts designated "
         "exploratory (GLASS bulk deconvolution, %s samples from %s cases; CGGA "
         "cross-sectional, %s IDH-mutant recurrent tumours) and a pre-declared "
         "immune-subclustering sensitivity analysis."
         % (n0("cohorts.discovery.n_samples"), n0("cohorts.discovery.n_patients"),
            n0("cohorts.discovery.n_any_rt_patients"),
            n0("cohorts.discovery.n_no_rt_patients"), n0("evaluability.threshold_samples"),
            n0("evaluability.n_samples_with_proportions"), n0("power.seed"),
            n0("power.nsim"), n0("cohorts.glass.n_samples"), n0("cohorts.glass.n_cases"),
            n0("cohorts.cgga.n_analysed"))),
        ("Results:", "Only %s of the %s pre-specified primary cell types were "
         "evaluable: NK cells (interaction beta %s, 95%% confidence interval %s to %s, "
         "Benjamini-Hochberg q %s, permutation p %s) and microglia (beta %s, 95%% "
         "confidence interval %s to %s, q %s, permutation p %s); neither differed from "
         "its permutation null. The design was underpowered, with %s power at f2 = %s "
         "and %s at f2 = %s at an intraclass correlation of 0.60. Five primary types "
         "formed no distinct cluster, and three (B cells, monocyte-derived "
         "macrophages, CD8 T cells) fell below the evaluability threshold. GLASS "
         "estimates (exploratory) reached q < 0.05 for B cells and microglia, and both "
         "discovery-evaluable types change direction between cohorts. The "
         "immune-subclustering sensitivity analysis failed its gate (%s, %s and %s of "
         "%s primary types at resolutions 0.5, 1.0 and 1.5)."
         % (ev_n, ev_tot, n3("discovery.lmm_by_cell.NK.beta"),
            n3("discovery.lmm_by_cell.NK.ci_low"), n3("discovery.lmm_by_cell.NK.ci_high"),
            n3("discovery.lmm_by_cell.NK.bh_q"), n3("discovery.lmm_by_cell.NK.perm_p"),
            n3("discovery.lmm_by_cell.Microglia_P2RY12_TMEM119.beta"),
            n3("discovery.lmm_by_cell.Microglia_P2RY12_TMEM119.ci_low"),
            n3("discovery.lmm_by_cell.Microglia_P2RY12_TMEM119.ci_high"),
            n3("discovery.lmm_by_cell.Microglia_P2RY12_TMEM119.bh_q"),
            n3("discovery.lmm_by_cell.Microglia_P2RY12_TMEM119.perm_p"),
            pct("power.headline_f2028_icc60"), "%.2f" % float(g("power.f2_anchors")[0]),
            pct("power.headline_f2040_icc60"), "%.2f" % float(g("power.f2_anchors")[1]),
            n0("sensitivity.evaluable_counts.0.5"), n0("sensitivity.evaluable_counts.1.0"),
            n0("sensitivity.evaluable_counts.1.5"), ev_tot)),
        ("Conclusion:", "This cohort cannot support or exclude a radiotherapy-associated "
         "immune shift at recurrence; annotation depth, not sample number, is the "
         "first constraint to address."),
    ]
    words = 0
    for head, txt in body:
        p = doc.add_paragraph()
        r = p.add_run(head + " ")
        r.bold = True
        p.add_run(txt)
        words += len(head.split()) + len(txt.split())
    doc.add_paragraph("Abstract word count: %d. Word limit for the congress abstract: "
                      "[VERIFY on the congress site]." % words)
    out = os.path.join(LST, "DGNC_abstract.docx")
    doc.save(out)
    return out, words


TODO = """# Submission checklist - fill in by hand

Everything below is a placeholder in the package and must be completed before
submission. Nothing here is sourced from the project files.

## Metadata
- [ ] Author list (order, initials, surnames) - replaces `[AUTHOR LIST]`
- [ ] Affiliations, with the affiliation index for each author - `[AFFILIATIONS]`
- [ ] Corresponding author and e-mail - `[CORRESPONDING AUTHOR] [EMAIL]`
- [ ] Acknowledgements - `[ACKNOWLEDGEMENTS]`
- [ ] Author contributions (CRediT roles) - `[AUTHOR CONTRIBUTIONS]`
- [ ] Funding statement, including grant numbers - `[FUNDING]`
- [ ] Competing interests for every author - `[CONFLICTS]`
- [ ] Ethics approval and consent to participate, with the responsible committee
      and approval number - `[ETHICS]`
- [ ] Consent for publication - `[CONSENT]`
- [ ] Generative-AI use statement - `[AI-USE STATEMENT]`
- [ ] Corresponding-author confirmation that all authors approved the submitted
      version

## Decisions
- [ ] Journal choice (the package is written to be journal-agnostic): candidate
      targets are a precision-oncology or neuro-oncology journal - decide, then
      apply that journal's structure, reference style and word limits
- [ ] Title: three options are listed in
      `manuscript/manuscript_metrics.json`; the first is used in the files
- [ ] Abstract format for the chosen journal (the 5-sentence version may need
      restructuring; the numbers are fixed)
- [ ] Congress abstract word limit and format - `[VERIFY on the congress site]`;
      see `checklists/DGNC_abstract.docx`
- [ ] Figure count and file-format requirements of the chosen journal (the package
      ships PNG 600 dpi, LZW TIFF 600 dpi, vector PDF and SVG for each figure)
- [ ] Whether the journal requires the pre-specification as a supplement (it is
      included unmodified)

## Citations to complete (`[VERIFY CITATION]`)
- [ ] Tang et al., Cancer Cell - year, volume, pages, DOI
- [ ] GLASS-NL cohort publication - full citation
- [ ] CGGA mRNAseq_693 primary dataset publication
- [ ] quanTIseq method publication
- [ ] EPIC method publication
- [ ] MCP-counter method publication
- [ ] Software citations (scanpy, anndata, statsmodels, matplotlib, numpy, pandas),
      or a single "software versions" citation as the journal allows
- [ ] Verify the two transcriptions already taken from the project resource file:
      Johnson et al. Nature 2026 (doi:10.1038/s41586-026-10612-6, PMID 42236943) and
      Makarevic et al. Int J Mol Sci 2020 (doi:10.3390/ijms21207801, PMID 33096928)

## Data-use and approvals
- [ ] GLASS / Synapse project syn17038081 - confirm the data-use terms permit this
      publication and record the accession wording the journal requires
- [ ] CGGA - confirm the data-use terms permit this publication
- [ ] GEO GSE326221 - confirm the accession citation wording
- [ ] Confirm that no controlled-access arm (EGA EGAD50000002485) contributed data
      to this report
- [ ] Institutional approval / institutional sign-off for the reuse of published
      aggregate data, if the institution requires it
- [ ] Confirm that the submitted files contain no patient-level or raw data (the
      package was built to exclude them; re-check after any manual edit)

## Package mechanics
- [ ] Re-run `code/submission_build_numbers.py`, `code/submission_make_figures.py`,
      `code/submission_make_tables.py` and `code/submission_build_manuscript.py`
      after any manual edit, then regenerate `checks/manifest.sha256`
- [ ] Re-open `manuscript/manuscript.docx` in Word and confirm line numbering,
      table grids and figure placement survived
- [ ] Confirm the manuscript word counts in the document still match
      `manuscript/manuscript_metrics.json` after editing
"""


def write_todo():
    p = os.path.join(LST, "submission_todo.md")
    with io.open(p, "w", encoding="utf-8", newline="\n") as f:
        f.write(TODO)
    return p


# ------------------------------------------------------------------ manifest
def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def write_manifest():
    p = os.path.join(CHK, "manifest.sha256")
    lines, n = [], 0
    for base, dirs, files in os.walk(SUB):
        dirs[:] = [d for d in dirs if d != "__pycache__"]
        for name in sorted(files):
            full = os.path.join(base, name)
            rel = os.path.relpath(full, SUB).replace("\\", "/")
            if rel == "checks/manifest.sha256":
                continue
            lines.append("%s  %s" % (sha256(full), rel))
            n += 1
    with io.open(p, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(sorted(lines)) + "\n")
    return p, n


def main():
    kept, skipped = copy_audits()
    ps = copy_prespec()
    ncode, npkg = copy_code()
    rd = write_readme(kept, skipped)
    ab, words = write_abstract_docx()
    td = write_todo()
    mf, nfiles = write_manifest()
    print("audits kept:    %d (%s)" % (len(kept), ", ".join(kept)))
    print("audits skipped: %d (%s)" % (len(skipped), "; ".join(skipped) or "none"))
    print("pre-spec copied unmodified: %s" % ", ".join(ps))
    print("code: %d scripts + %d pinned packages" % (ncode, npkg))
    print("README:         %s" % rd)
    print("congress abstract: %s (%d words)" % (ab, words))
    print("todo:           %s" % td)
    print("manifest:       %s (%d files hashed)" % (mf, nfiles))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
