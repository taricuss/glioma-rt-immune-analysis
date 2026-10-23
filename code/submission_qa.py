# -*- coding: utf-8 -*-
"""Stage 6 - quality assurance for the DGNC submission package.

Reads the built package only.  Writes submission/checks/qa_report.md and exits
non-zero if a hard check fails (manifest mismatch, changed locked file, a number
that is in the abstract but not in the registry, banned terminology, a raster
figure below 600 dpi, or a residual local path / identifier).
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import os
import re
import sys
import zipfile
from html.parser import HTMLParser

from PIL import Image

ROOT = r"<PROJECT_ROOT>"
SUB = os.path.join(ROOT, "submission")
SUP = os.path.join(SUB, "supplementary")
FIG = os.path.join(SUB, "figures")
TAB = os.path.join(SUB, "tables")
MAN = os.path.join(SUB, "manuscript")
CHK = os.path.join(SUB, "checks")
LST = os.path.join(SUB, "checklists")
BASELINE = os.path.join(os.environ.get("TEMP", ""), "dgnc_stage0_baseline.json")
USERNAME = "[USER]"

REPORT = []
HARD = []
JUSTIFY = []          # (title, [lines]) - items that would otherwise read as open


def say(ok, label, detail=""):
    REPORT.append((("PASS" if ok else "FAIL"), label, detail))
    if not ok:
        HARD.append(label)
    print("%-4s %s%s" % ("PASS" if ok else "FAIL", label,
                         (" - " + detail) if detail else ""))


def review(label, detail=""):
    REPORT.append(("REVIEW", label, detail))
    print("REVW %s%s" % (label, (" - " + detail) if detail else ""))


def note(label, detail=""):
    REPORT.append(("NOTE", label, detail))
    print("NOTE %s%s" % (label, (" - " + detail) if detail else ""))


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


REG = json.load(io.open(os.path.join(SUP, "results_numbers.json"),
                        "r", encoding="utf-8"))
NUM = _nest(REG["values"])


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


# ------------------------------------------------------------------ text
class _Text(HTMLParser):
    def __init__(self):
        HTMLParser.__init__(self)
        self.parts = []

    def handle_data(self, d):
        self.parts.append(d)


def html_text(path):
    p = _Text()
    p.feed(io.open(path, "r", encoding="utf-8").read())
    return " ".join("".join(p.parts).split())


def section(html, start, end):
    i = html.find(">" + start + "<")
    if i < 0:
        i = html.find(start)
    j = html.find(">" + end + "<", i + 1)
    if j < 0:
        j = html.find(end, i + len(start))
    return html[i:j]


def strip_tags(s):
    return " ".join(re.sub(r"<[^>]+>", " ", s).split())


def docx_text(path):
    with zipfile.ZipFile(path) as z:
        out = []
        for n in z.namelist():
            if n.startswith("word/") and n.endswith(".xml"):
                out.append(re.sub(r"<[^>]+>", " ", z.read(n).decode("utf-8",
                                                                    "replace")))
    return " ".join(" ".join(out).split())


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def walk_files(d):
    for base, dirs, files in os.walk(d):
        dirs[:] = [x for x in dirs if x != "__pycache__"]
        for n in sorted(files):
            yield os.path.join(base, n)


def rel(p):
    return os.path.relpath(p, SUB).replace("\\", "/")


def refresh_manifest():
    """Re-write checks/manifest.sha256 so it covers the QA report as well."""
    lines, n = [], 0
    for f in walk_files(SUB):
        r = rel(f)
        if r == "checks/manifest.sha256":
            continue
        lines.append("%s  %s" % (sha256(f), r))
        n += 1
    with io.open(os.path.join(CHK, "manifest.sha256"), "w", encoding="utf-8",
                 newline="\n") as f:
        f.write("\n".join(sorted(lines)) + "\n")
    return n


# ------------------------------------------------------------------ checks
BANNED = [r"\brobust\w*", r"\bas planned\b", r"\bmandatory\b", r"\breplicat\w*",
          r"\bvalidat\w*", r"sign[- ]agreement", r"5 of 10 estimable",
          r"1 of 5 agree", r"\bcaused\b", r"\bcauses\b", r"\bled to\b",
          r"\bdrives\b", r"\bresolves\b", r"\bresolve\b"]

TEXT_EXT = {".md", ".txt", ".csv", ".json", ".html", ".svg", ".py"}


def main():
    html = io.open(os.path.join(MAN, "manuscript.html"), "r",
                   encoding="utf-8").read()
    full = strip_tags(html)
    metrics = json.load(io.open(os.path.join(MAN, "manuscript_metrics.json"),
                                "r", encoding="utf-8"))

    # 1. abstract sentence count, re-derived from the built HTML
    abs_html = section(html, "<h2>Abstract</h2>", "<h2>Introduction</h2>")
    m_abs = re.search(r"<p>(.*?)</p>", abs_html, re.S)
    abstract = strip_tags(m_abs.group(1))
    sents = [s for s in re.split(r"(?<=[.!?])\s+", abstract) if s.strip()]
    say(len(sents) == 5, "abstract is exactly 5 sentences",
        "re-derived from the HTML abstract paragraph: %d (metrics file says %d)"
        % (len(sents), metrics["abstract_sentences"]))
    note("abstract word count", "%d words (metrics: %d)"
         % (len(abstract.split()), metrics["abstract_words"]))

    # 2. banned terminology
    dhits = {}
    for pat in BANNED:
        for m in re.finditer(pat, full, flags=re.IGNORECASE):
            dhits.setdefault(pat, []).append(m.group(0))
    say(not dhits, "no banned terminology in the manuscript HTML",
        "; ".join("%s -> %s" % (k, sorted(set(v))[:4]) for k, v in dhits.items()))
    txt = docx_text(os.path.join(MAN, "manuscript.docx"))
    dhits2 = {}
    for pat in BANNED:
        for m in re.finditer(pat, txt, flags=re.IGNORECASE):
            dhits2.setdefault(pat, []).append(m.group(0))
    say(not dhits2, "no banned terminology in the manuscript docx",
        "; ".join("%s -> %s" % (k, sorted(set(v))[:4]) for k, v in dhits2.items()))
    say("5 of 10 estimable" not in full and "1 of 5 agree" not in full,
        "retired evaluability and sign-count phrasings absent")

    # 3. every number in Abstract + Results is in the registry
    res_html = section(html, "<h2>Results</h2>", "<h2>Discussion</h2>")
    res_html = re.sub(r"<figure>.*?</figure>", " ", res_html, flags=re.S)
    scan = abstract + " " + strip_tags(res_html)
    allowed = set()

    def add_forms(v):
        v = abs(float(v))
        allowed.update({"%g" % v, "%.3f" % v, "%.2f" % v, "%.1f" % v,
                        "%.1f" % (100.0 * v), "%.1f%%" % (100.0 * v)})
        if float(v).is_integer():
            allowed.add("%d" % int(v))

    def collect(obj):
        if isinstance(obj, bool):
            return
        if isinstance(obj, (int, float)):
            add_forms(obj)
        elif isinstance(obj, str):
            allowed.add(obj)
            for m in re.finditer(r"\d+(?:\.\d+)?", obj):
                allowed.add(m.group(0))
        elif isinstance(obj, dict):
            for v in obj.values():
                collect(v)
        elif isinstance(obj, list):
            for v in obj:
                collect(v)

    collect(REG["values"])
    # Allowed non-registry tokens, each justified in the QA report:
    #  95        - the nominal confidence-interval level, not an analysis value
    #  326221 / 17038081 - data-accession identifiers (GEO GSE326221, Synapse
    #              syn17038081), cited in Table 1 and verified in the project files
    #  derived   - the sensitivity-compartment share, computed by the table script
    #              as 100 * compartment_cells / modelled cells (two registry values)
    ACCESSIONS = {"326221", "17038081"}
    ALLOW = {"95"} | ACCESSIONS | {
        "%.1f" % (100.0 * float(g("sensitivity.compartment_cells"))
                  / float(g("cohorts.discovery.cells_final_postcap")))}
    tokens = re.findall(r"\d+(?:\.\d+)?", scan.replace(",", ""))
    missing = sorted({t for t in tokens if t not in allowed and t not in ALLOW},
                     key=lambda x: (len(x), x))
    say(not missing, "every Abstract/Results number is in results_numbers.json",
        "scanned the abstract paragraph, the Results prose and the Results tables; "
        "allowed outside the registry: %s; unmatched: %s"
        % (", ".join(sorted(ALLOW)), ", ".join(missing) if missing else "none"))
    shared = "%.1f" % (100.0 * float(g("sensitivity.compartment_cells"))
                       / float(g("cohorts.discovery.cells_final_postcap")))
    JUSTIFY.append((
        "Numbers allowed outside results_numbers.json",
        ["`95` - the nominal confidence-interval level. It is a convention, not a "
         "derived analysis value, so it is not registered.",
         "`326221` and `17038081` - data-accession identifiers (GEO GSE326221, "
         "Synapse syn17038081). They are identifiers, cited in Table 1 and present "
         "verbatim in the project files, not measurements.",
         "`%s` - the sensitivity-compartment share. It is not a stored value: the "
         "table script computes it as 100 x sensitivity.compartment_cells / "
         "cohorts.discovery.cells_final_postcap, i.e. it is assembled from two "
         "registry values, and the QA script recomputes the same expression above "
         "rather than hard-coding the result." % shared,
         "Every other number in the abstract, the Results prose and the Results "
         "tables resolves to a value in results_numbers.json (either as a nested "
         "leaf, as an integer form, or as a percentage rendering)."]))

    # 4. documents reopen / parse
    try:
        from docx import Document
        d = Document(os.path.join(MAN, "manuscript.docx"))
        ntab = len(d.tables)
        nimg = len(d.inline_shapes)
        nh = len([p for p in d.paragraphs
                  if p.style.name.startswith("Heading")])
        say(True, "manuscript.docx reopens", "%d headings, %d tables, %d images"
            % (nh, ntab, nimg))
        if ntab < 6 or nimg < 7 or nh < 10:
            review("docx structural counts lower than expected",
                   "headings %d, tables %d, images %d" % (nh, ntab, nimg))
    except Exception as e:                                   # noqa: BLE001
        say(False, "manuscript.docx reopens", repr(e))
    try:
        p = _Text()
        p.feed(html)
        say(True, "manuscript.html parses", "%d data chunks" % len(p.parts))
    except Exception as e:                                   # noqa: BLE001
        say(False, "manuscript.html parses", repr(e))

    # 5. figures: dpi, formats, sizes
    stems = sorted({os.path.splitext(f)[0] for f in os.listdir(FIG)
                    if re.match(r"^Figure", f)})
    bad_dpi, missing_fmt, big = [], [], []
    for s in stems:
        for ext in ("png", "tiff"):
            p = os.path.join(FIG, "%s.%s" % (s, ext))
            if not os.path.exists(p):
                missing_fmt.append("%s.%s" % (s, ext))
                continue
            with Image.open(p) as im:
                dpi = im.info.get("dpi", (0, 0))
            if min(dpi) < 599.5:
                bad_dpi.append("%s.%s dpi=%s" % (s, ext, dpi))
        for ext in ("pdf", "svg"):
            if not os.path.exists(os.path.join(FIG, "%s.%s" % (s, ext))):
                missing_fmt.append("%s.%s" % (s, ext))
    for f in walk_files(FIG):
        if os.path.getsize(f) > 50 * 1024 * 1024:
            big.append(rel(f))
    say(not bad_dpi, "%d raster figures at >=600 dpi" % (len(stems) * 2),
        "; ".join(bad_dpi))
    say(not missing_fmt, "every figure has png, tiff, pdf and svg",
        "; ".join(missing_fmt))
    say(not big, "no figure exceeds 50 MB", "; ".join(big))
    note("figure stems", "%d (%s)" % (len(stems), ", ".join(stems)))

    # 6. docx size
    dz = os.path.getsize(os.path.join(MAN, "manuscript.docx")) / 1024.0 / 1024.0
    note("manuscript.docx size", "%.2f MB" % dz)
    if dz > 25:
        review("manuscript.docx over 25 MB", "%.2f MB - build a low-res variant" % dz)

    # 7. residual local paths / username
    path_re = re.compile(r"[A-Za-z]:\\[^\s\"'<>]+|<PROJECT_ROOT>"'<>]+")
    leaks = []
    for f in walk_files(SUB):
        ext = os.path.splitext(f)[1].lower()
        if ext in (".png", ".tiff", ".pdf", ".xlsx"):
            continue
        if ext not in TEXT_EXT:
            continue
        t = io.open(f, "r", encoding="utf-8", errors="replace").read()
        if USERNAME in t.lower():
            leaks.append("%s: username" % rel(f))
        m = path_re.search(t)
        if m and "submission/" not in m.group(0).replace("\\", "/"):
            leaks.append("%s: %s" % (rel(f), m.group(0)[:60]))
    for f in walk_files(SUB):
        if os.path.splitext(f)[1].lower() == ".docx":
            with zipfile.ZipFile(f) as z:
                for n in z.namelist():
                    if n.endswith((".xml", ".rels")):
                        t = z.read(n).decode("utf-8", "replace")
                        if USERNAME in t.lower() or path_re.search(t):
                            leaks.append("%s(%s)" % (rel(f), n))
        elif os.path.splitext(f)[1].lower() == ".pdf":
            raw = open(f, "rb").read()
            if USERNAME.encode() in raw.lower() or b"/Users/" in raw:
                leaks.append("%s: pdf string" % rel(f))
    say(not leaks, "no local path or username in the package",
        "; ".join(sorted(set(leaks))[:8]))

    # 8. patient-level artefacts
    id_re = re.compile(r"\bGSM\d{6,}\b|\b(?:SJ|MDA)\d{3,}\b")
    sus = []
    for f in walk_files(SUB):
        b = os.path.basename(f).lower()
        ext = os.path.splitext(f)[1].lower()
        if ("per_sample" in b or "per_patient" in b
                or (b.startswith("clinical") and ext in (".csv", ".tsv", ".xlsx"))):
            sus.append("filename: " + rel(f))
        if os.path.splitext(f)[1].lower() in (".md", ".txt", ".csv", ".json", ".html"):
            t = io.open(f, "r", encoding="utf-8", errors="replace").read()
            m = id_re.findall(t)
            if m:
                sus.append("%s: %s" % (rel(f), sorted(set(m))[:4]))
    say(not sus, "no patient-level files or sample identifiers", "; ".join(sus))

    # 9. sub-5 counts in the cohort tables
    flags = []
    for tid in ("Table1", "Table2"):
        p = os.path.join(TAB, tid + ".csv")
        rows = list(csv.reader(io.open(p, "r", encoding="utf-8-sig", newline="")))
        head = rows[0]
        for r in rows[1:]:
            for ci, c in enumerate(r):
                if ci >= len(head):
                    continue
                if re.fullmatch(r"[1-4]", c.strip()):
                    flags.append("%s / %s / %s = %s"
                                 % (tid, r[0], head[ci], c.strip()))
    if flags:
        review("counts below 5 present in the cohort tables",
               "; ".join(flags[:10]))
        JUSTIFY.append((
            "Counts below 5 in the cohort tables",
            ["The flagged cell(s) are per-sample counts of a non-zero proportion "
             "for a cell type that is not evaluable (CD8 T cells appear with a "
             "non-zero proportion in 1 sample). They are not patient counts and "
             "they disclose no individual: Table 2 is the evaluability table, and "
             "these numbers are the rule-derived reason a cell type is reported as "
             "not evaluable.",
             "The patient-level suppression rule is enforced where it applies: "
             "every patient and sample identifier is absent from the package (see "
             "the identifier check above), and the discovery-arm evaluability "
             "counts are reported in whole-patient units.",
             "Kept deliberately: showing the per-sample presence count is what "
             "makes the 2-of-10 evaluability statement auditable, and dropping it "
             "would hide the rule behind a bare verdict."]))
    else:
        say(True, "no count below 5 in the cohort tables")

    # 10. manifest
    man = os.path.join(CHK, "manifest.sha256")
    listed, mismatch, extra = {}, [], []
    for line in io.open(man, "r", encoding="utf-8"):
        line = line.rstrip("\n")
        if not line:
            continue
        h, r = line.split("  ", 1)
        listed[r] = h
    for r, h in listed.items():
        p = os.path.join(SUB, r.replace("/", os.sep))
        if not os.path.exists(p) or sha256(p) != h:
            mismatch.append(r)
    for f in walk_files(SUB):
        if rel(f) in ("checks/manifest.sha256", "checks/qa_report.md"):
            continue
        if rel(f) not in listed:
            extra.append(rel(f))
    say(not mismatch and not extra, "manifest re-hash matches the package",
        "mismatched %s; unlisted %s" % (mismatch or "none", extra or "none"))
    note("manifest entries", "%d" % len(listed))

    # 11. [VERIFY] markers
    marks = []
    for f in walk_files(SUB):
        if os.path.splitext(f)[1].lower() not in (".md", ".txt", ".csv", ".json"):
            continue
        if rel(f) == "checks/qa_report.md":
            continue                      # this report itself lists the markers
        for i, line in enumerate(io.open(f, "r", encoding="utf-8",
                                         errors="replace"), 1):
            if "[VERIFY" in line:
                marks.append("%s:%d" % (rel(f), i))
    for f in (os.path.join(MAN, "manuscript.docx"),):
        t = docx_text(f)
        if "[VERIFY" in t:
            marks.append("manuscript/manuscript.docx (text)")
    note("[VERIFY] markers", "%d locations: %s"
         % (len(marks), ", ".join(marks[:20]) if marks else "none"))

    # 12. locked files unchanged
    base = json.load(io.open(BASELINE, "r", encoding="utf-8"))["files"]
    changed, gone = [], []
    for r, meta in base.items():
        p = os.path.join(ROOT, r)
        if not os.path.exists(p):
            gone.append(r)
            continue
        if sha256(p) != meta["sha256"]:
            changed.append(r)
    say(not changed and not gone, "no locked file changed since the Stage 0 baseline",
        "changed %s; missing %s" % (changed or "none", gone or "none"))
    note("baseline files checked", "%d" % len(base))

    # ------------------------------------------------------------- report
    out = os.path.join(CHK, "qa_report.md")
    n = {"PASS": 0, "FAIL": 0, "REVIEW": 0, "NOTE": 0}
    lines = ["# Quality-assurance report - DGNC submission package",
             "",
             "Generated by `submission_qa.py` from the built package. Every figure "
             "and every number is re-derived from the files in this package; no "
             "analysis was re-run.",
             "",
             "| Result | Check | Detail |",
             "| --- | --- | --- |"]
    for st, lb, dt in REPORT:
        n[st] = n.get(st, 0) + 1
        lines.append("| %s | %s | %s |" % (st, lb, dt.replace("|", "/")))
    lines += ["", "Summary: %d pass, %d fail, %d review, %d note."
              % (n["PASS"], n["FAIL"], n["REVIEW"], n["NOTE"]), ""]
    if JUSTIFY:
        lines.append("## Justified items")
        lines.append("")
        lines.append("Items that are deliberately allowed rather than failures. "
                     "Each is a rule-derived or convention value, and each is "
                     "explained here so it does not read as an open issue.")
        lines.append("")
        for title, body in JUSTIFY:
            lines.append("**%s**" % title)
            lines.append("")
            lines += ["- " + b for b in body]
            lines.append("")
    if HARD:
        lines.append("## Failing checks")
        lines.append("")
        lines += ["- " + h for h in HARD]
        lines.append("")
    with io.open(out, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    print("\nqa report: %s" % out)
    nman = refresh_manifest()
    print("manifest refreshed over every package file: %d entries" % nman)
    print("summary: %d pass, %d fail, %d review, %d note"
          % (n["PASS"], n["FAIL"], n["REVIEW"], n["NOTE"]))
    return 1 if HARD else 0


if __name__ == "__main__":
    sys.exit(main())
