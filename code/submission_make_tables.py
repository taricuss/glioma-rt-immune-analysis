# -*- coding: utf-8 -*-
"""Stage 3 - tables for the DGNC submission package.

Reads ONLY existing aggregate artefacts + the Stage 1 numbers registry.
Never imports a pipeline script; the marker dictionary is read with `ast`.

Outputs to submission/tables/: Table1..Table6, TableS1..TableS6 as .csv and .xlsx.
The manuscript builder reads the same .csv files, so text and document cannot diverge.
"""
from __future__ import annotations

import ast
import csv
import io
import json
import os
import subprocess
import sys

ROOT = r"<PROJECT_ROOT>"
INT = os.path.join(ROOT, "outputs", "intermediate")
OUT = os.path.join(ROOT, "submission")
TAB = os.path.join(OUT, "tables")
os.makedirs(TAB, exist_ok=True)

PRIMARY_SHORT = ["CD8_T", "CD4_conv", "Treg_FOXP3", "NK", "B", "Plasma",
                 "Microglia_P2RY12_TMEM119", "MoMF_CD163", "DC", "Neutrophils"]
PRIMARY_LONG = {
    "CD8_T": "CD8_T_cells", "CD4_conv": "CD4_T_cells_conventional",
    "Treg_FOXP3": "T_regulatory_cells_FOXP3", "NK": "NK_cells",
    "B": "B_cells", "Plasma": "Plasma_cells",
    "Microglia_P2RY12_TMEM119": "Microglia_P2RY12_TMEM119",
    "MoMF_CD163": "Monocyte_derived_macrophages_CD163",
    "DC": "Dendritic_cells", "Neutrophils": "Neutrophils",
}
DISP = {"CD8_T": "CD8 T cells", "CD4_conv": "CD4 conventional T cells",
        "Treg_FOXP3": "Regulatory T cells (FOXP3)", "NK": "NK cells",
        "B": "B cells", "Plasma": "Plasma cells",
        "Microglia_P2RY12_TMEM119": "Microglia (P2RY12/TMEM119)",
        "MoMF_CD163": "Monocyte-derived macrophages (CD163)",
        "DC": "Dendritic cells", "Neutrophils": "Neutrophils"}

def _nest(flat):
    """The registry stores dotted flat keys; rebuild the tree so table code can
    address values as nested dicts (keys may themselves contain dots)."""
    out = {}
    for k, v in flat.items():
        cur = out
        parts = k.split(".")
        for p in parts[:-1]:
            cur = cur.setdefault(p, {})
        cur[parts[-1]] = v
    return out


NUM = _nest(json.load(io.open(os.path.join(OUT, "supplementary",
                                           "results_numbers.json"),
                              "r", encoding="utf-8"))["values"])


def g(name):
    if name in NUM:
        return NUM[name]
    parts = name.split(".")
    for i in range(len(parts), 0, -1):
        pre = ".".join(parts[:i])
        if pre in NUM:
            cur = NUM[pre]
            rest = parts[i:]
            while rest:
                for j in range(len(rest), 0, -1):
                    k = ".".join(rest[:j])
                    if isinstance(cur, dict) and k in cur:
                        cur, rest = cur[k], rest[j:]
                        break
                else:
                    raise KeyError(name)
            return cur
    raise KeyError(name)


def rd_csv(rel):
    with io.open(os.path.join(ROOT, rel), "r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def f3(x):
    return "" if x is None or x == "" else "%.3f" % float(x)


TABLES = {}


def add(tid, title, header, rows, notes=""):
    TABLES[tid] = {"title": title, "header": header, "rows": rows, "notes": notes}


def write_csv(tid, obj):
    p = os.path.join(TAB, tid + ".csv")
    with io.open(p, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(obj["header"])
        for r in obj["rows"]:
            w.writerow(r)
        if obj["notes"]:
            w.writerow([])
            w.writerow(["note", obj["notes"]])
    return p


def write_xlsx(tid, obj):
    from openpyxl import Workbook
    from openpyxl.styles import Font, Alignment
    from openpyxl.utils import get_column_letter
    wb = Workbook()
    ws = wb.active
    ws.title = tid[:31]
    ws.append(obj["header"])
    for c in ws[1]:
        c.font = Font(bold=True)
        c.alignment = Alignment(vertical="top", wrap_text=True)
    for r in obj["rows"]:
        ws.append(list(r))
    for col in range(1, len(obj["header"]) + 1):
        letter = get_column_letter(col)
        width = max(len(str(obj["header"][col - 1])) + 2,
                    *(len(str(r[col - 1])) + 2 for r in obj["rows"])) if obj["rows"] else 12
        ws.column_dimensions[letter].width = min(max(width, 10), 52)
        for row in range(2, len(obj["rows"]) + 2):
            cell = ws.cell(row=row, column=col)
            if isinstance(cell.value, float):
                cell.number_format = "0.000"
    if obj["notes"]:
        ws.append([])
        ws.append(["note", obj["notes"]])
    p = os.path.join(TAB, tid + ".xlsx")
    wb.save(p)
    return p


# =============================================================== Table 1
def t1():
    d, gl, cg = g("cohorts.discovery"), g("cohorts.glass"), g("cohorts.cgga")
    rows = [
        ["Discovery snRNA-seq (Johnson et al. 2026; GSE326221 / syn17038081)",
         "PRIMARY - pre-specified analysis set",
         d["n_samples"], d["n_patients"], d["n_any_rt_patients"], d["n_no_rt_patients"],
         "%d snRNA-seq samples passed QC; %d joined to treatment annotation and were "
         "modelled. Grade II at initial diagnosis. Per-sample cap %d cells."
         % (d["samples_input_snrna"], d["n_samples"], d["cap_cells_per_sample"])],
        ["GLASS RNA silver-set (IDH-mutant, RT-annotated, TP+R1 paired)",
         "EXPLORATORY - external direction check",
         gl["n_samples"], gl["n_cases"], gl["anyrt_pairs"], gl["nort_pairs"],
         "Bulk RNA-seq deconvolved with stock published signatures; %d / %d paired "
         "samples at primary and recurrent. No matched single-nucleus data."
         % (gl["anyrt_pairs"], gl["nort_pairs"])],
        ["CGGA mRNAseq_693 (%s)" % cg["version"],
         "EXPLORATORY - cross-sectional arm",
         cg["n_analysed"], "not recorded separately", cg["anyrt_samples"], cg["nort_samples"],
         "IDH-mutant recurrent samples with non-missing treatment; cross-sectional "
         "(no matched primary tumour), so no interaction can be estimated. "
         "Of %d IDH-mutant recurrent samples, %d had non-missing treatment."
         % (cg["idh_mut_recurrent_total"], cg["n_analysed"])],
    ]
    add("Table1", "Cohorts contributing to the analysis",
        ["Cohort", "Role in this study", "Samples (n)", "Patients or cases (n)",
         "Any_RT (n)", "No_RT (n)", "Notes"], rows,
        "Counts only; no patient-level or per-sample values are reported. Where a "
        "stratum contained fewer than five patients it would be reported as "
        "'fewer than 5'; no reported stratum falls below this. The bulk RNA arm of the "
        "EGA Yale DAC was not used (raw reads exceeded local processing capacity and "
        "added little over the open single-nucleus arm).")


# =============================================================== Table 2
def t2():
    prim = g("evaluability.by_primary_cell_type")
    dist = g("cohorts.discovery.annotation_distribution")
    rows = []
    for s in PRIMARY_SHORT:
        v = prim[s]
        cells = dist.get(PRIMARY_LONG[s], 0)
        if v["evaluable"]:
            status = "Evaluable"
            reason = "Present in >=10 of 75 samples"
        elif not v["has_proportion_column"]:
            status = "Not evaluable"
            reason = ("No distinct Leiden cluster formed under cluster-level argmax "
                      "marker-score annotation, so no per-sample proportion exists")
        else:
            status = "Not evaluable"
            reason = "Present in fewer than 10 of 75 samples"
        rows.append([DISP[s], PRIMARY_LONG[s], cells, v["present_samples"],
                     v["evaluable"], status, reason])
    add("Table2", "Evaluability of the 10 pre-specified primary cell types (discovery)",
        ["Cell type", "Annotation label", "Cells assigned (n)",
         "Samples with a non-zero proportion (of 75)", "Evaluable (>=10 samples)",
         "Status", "Reason"], rows,
        "A limitation of marker-score annotation rather than absence of the cells: "
        "NK and CD8 T cells share markers and NK cells were absorbed into "
        "CD8-T-winning clusters, which themselves "
        "failed the margin and marker thresholds. The cell counts are cluster-level "
        "assignments, not patient counts.")


# =============================================================== Table 3
def t3():
    disc = g("discovery.lmm_by_cell")
    rows = []
    for s in PRIMARY_SHORT:
        v = disc[s]
        if v["beta"] is None:
            continue
        tag = "PRIMARY - evaluable" if v["evaluable"] else "NOT EVALUABLE; descriptive only"
        rows.append([DISP[s], f3(v["beta"]), f3(v["ci_low"]), f3(v["ci_high"]),
                     f3(v["asym_p"]), f3(v["bh_q"]), f3(v["perm_p"]), tag])
    add("Table3", "Discovery interaction estimates (arcsin-square-root fraction)",
        ["Cell type", "Beta (timepoint x Any_RT)", "95% CI low", "95% CI high",
         "Asymptotic p", "BH q (K=10)", "Permutation p", "Status"], rows,
        "Model: arcsin-sqrt(fraction) ~ timepoint * treatment_arm + purity + grade + "
        "age + 1p19q + (1|patient). Confidence intervals are profile-likelihood "
        "(Wald fallback). BH q is the Benjamini-Hochberg value over the 10 primary "
        "types; the permutation p is from 10,000 patient-level treatment-label "
        "shuffles (seed 7331). Only NK and Microglia are evaluable; rows for B, MoMF "
        "and CD8 T cells carry finite model estimates but are below the evaluability "
        "threshold and are descriptive only. The remaining five primary types have no "
        "discovery estimate and are omitted here (see Table 2).")


# =============================================================== Table 4
def t4():
    glass = g("glass.lmm_by_cell")
    rows = []
    for s in PRIMARY_SHORT:
        v = glass[s]
        if v["beta"] is None:
            rows.append([DISP[s], "", "", "", "", "", "No estimate"])
            continue
        sig = "EXPLORATORY - BH q < 0.05" if (v["bh_q"] is not None and v["bh_q"] < 0.05) \
            else "EXPLORATORY"
        rows.append([DISP[s], f3(v["beta"]), f3(v["ci_low"]), f3(v["ci_high"]),
                     f3(v["asym_p"]), f3(v["bh_q"]), sig])
    add("Table4", "GLASS external-cohort interaction estimates (EXPLORATORY)",
        ["Cell type", "Beta (timepoint x Any_RT)", "95% CI low", "95% CI high",
         "Asymptotic p", "BH q (K=10)", "Status"], rows,
        "EXPLORATORY. Bulk deconvolution with stock published signatures; no matched "
        "single-nucleus reference. Purity was not identifiable in this cohort (a "
        "constant proxy was absorbed by the intercept). These estimates are not "
        "confirmatory and do not stand as independent evidence for the discovery "
        "analysis. "
        "In discovery only NK and Microglia were evaluable; both flip sign relative to "
        "GLASS, which is why a signed count of agreeing cell types is not used.")


# =============================================================== Table 5
def t5():
    ea = g("sensitivity.evaluable_by_resolution")
    rows = []
    for s in PRIMARY_SHORT:
        longn = PRIMARY_LONG[s]
        cells = []
        for r in ("0.5", "1.0", "1.5"):
            hit = [e for e in ea[r] if e["cell_type"] == longn]
            cells.append("evaluable" if (hit and hit[0]["evaluable"]) else "not evaluable")
        rows.append([DISP[s]] + cells)
    rows.append(["Evaluable cell types (n of 10)"] +
                [str(g("sensitivity.evaluable_counts")[r]) for r in ("0.5", "1.0", "1.5")])
    rows.append(["Gate threshold (n of 10)", "5", "6", "5"])
    rows.append(["Gate outcome", "fail", "fail", "fail"])
    add("Table5", "Pre-declared immune-subclustering sensitivity analysis",
        ["Cell type", "r = 0.5", "r = 1.0", "r = 1.5"], rows,
        "SENSITIVITY. The compartment comprised %d cells (%.1f%% of all cells); PTPRC "
        "was absent from the 2,000-gene highly variable space, so cells outside the "
        "annotated immune clusters could not be recovered. The pre-declared gate "
        "required at least 6 of 10 cell types at r = 1.0 and at least 5 of 10 at both "
        "r = 0.5 and r = 1.5; it failed, and no mechanistic analysis was performed on "
        "this annotation."
        % (g("sensitivity.compartment_cells"),
           100.0 * g("sensitivity.compartment_cells") / g("cohorts.discovery.cells_final_postcap")))


# =============================================================== Table 6
def t6():
    rows = [
        ["MGMT methylation covariate dropped",
         "The mgmt_methylation fixed effect in the locked v1.0 discovery formula was "
         "removed.",
         "MGMT is absent from all eight of the Johnson supplementary covariate sheets "
         "used to build the clinical canonical table.",
         "All discovery estimates are fitted without MGMT adjustment; the effective "
         "formula is reported alongside the locked formula."],
        ["1p19q status derived as an IDH-subtype proxy",
         "codeletion_1p19q was mapped from the recorded molecular subtype "
         "(Oligodendroglioma -> codel; Astrocytoma -> non_codel) rather than a "
         "molecular 1p19q call.",
         "No molecular 1p19q call is available in any supplementary sheet.",
         "1p19q enters the model as a two-level proxy covariate (43 non-codel / 28 "
         "codel patients)."],
        ["Tumour purity derived from single-nucleus malignant-cell fraction",
         "estimate_purity was computed as rna_malignant_cell_number / rna_cell_number "
         "instead of a bulk ESTIMATE or ABSOLUTE purity.",
         "No bulk purity estimate is available for the single-nucleus discovery set.",
         "Purity enters the discovery model as a continuous covariate (median 0.729)."],
        ["Cell-level QC applied on the full gene space",
         "Cell QC (min_genes = 300) is evaluated on the full 36,601-gene expression "
         "matrix before the ~2,000 highly variable gene subset is taken.",
         "Applying the same thresholds after HVG subsetting removed roughly four times "
         "more cells than intended and lost 33 of 75 samples.",
         "All discovery proportions derive from the recovered 75-sample design "
         "(%d cells)." % g("cohorts.discovery.cells_final_postcap")],
        ["Per-sample cell cap of 3,000 cells",
         "Each sample was uniformly subsampled to at most %d cells after QC "
         "(%d of %d samples capped), reducing %d cells to %d."
         % (g("cohorts.discovery.cap_cells_per_sample"),
            g("cohorts.discovery.n_samples_capped"),
            g("cohorts.discovery.samples_input_snrna"),
            g("cohorts.discovery.cells_precap_total"),
            g("cohorts.discovery.cells_final_postcap")),
         "Bounded peak memory for the merge step on the available hardware.",
         "Rare cell types retain fewer cells, but the cap is applied uniformly and "
         "does not preferentially remove them."],
        ["Cluster-level mean argmax marker-score annotation",
         "Cell types were assigned from the mean marker score of each Leiden cluster "
         "(r = 0.8) rather than per-cell scoring.",
         "Preserves cluster coherence and keeps cell-type assignment auditable at the "
         "cluster level.",
         "Only 5 of 10 primary types receive any cells; NK and Microglia are the only "
         "evaluable types."],
        ["Optimizer fallback and zero-variance covariate dropping",
         "The mixed model is fitted with lbfgs, retrying with bfgs then powell when "
         "lbfgs raises a singular-matrix error; covariates with zero variance in a "
         "given cohort are dropped from that cohort's effective formula.",
         "Numerical optimisation failure and collinear terms, not a change in the "
         "hypothesis being tested.",
         "In GLASS, the constant purity proxy was dropped; discovery retained all "
         "terms. Optimizer identity is recorded per cell type in the model audit."],
        ["Post-data power simulated with 800 simulation draws",
         "Power was estimated by Monte-Carlo simulation with NSIM = 800 and a "
         "patient-level design (seed %d), at ICC 0.40/0.60/0.70."
         % g("power.seed"),
         "A post-data power check gives an honest ceiling on what the recovered "
         "design could have detected.",
         "At ICC 0.60 the design reaches %.1f%% power at f2 = 0.28 and %.1f%% at "
         "f2 = 0.40; the headline label is %s."
         % (g("power.headline_f2028_icc60") * 100,
            g("power.headline_f2040_icc60") * 100, g("power.headline_label"))],
        ["GLASS deconvolution uses stock signatures with an added Plasma and "
         "Microglia row set, and no purity adjustment",
         "GLASS fractions come from stock published signatures (quanTIseq/TIL10 for "
         "the primary fit) rather than the discovery custom signature; TIL10 was "
         "extended with canonical Plasma and Microglia marker rows; no purity term is "
         "fitted.",
         "GLASS has no matched single-nucleus reference for a custom signature, and "
         "TIL10 models neither plasma cells nor brain-resident microglia; purity is "
         "not identifiable because the signature-mass proxy is near-degenerate.",
         "GLASS is reported as EXPLORATORY only. Verified in "
         "scripts/deconvolution_glass_pipeline.py (TIL10 extension block) and "
         "glass_pipeline_audit.json."],
        ["Forest plot markers filled by BH q, no pooled diamond",
         "Markers are filled when BH q < 0.05 and hollow otherwise, and no pooled "
         "meta-analytic estimate is drawn.",
         "Pooling opposite-signed estimates across cohorts would be read as agreement "
         "between cohorts, which the data do not support; this deviates from "
         "pre-specification "
         "section 4.5(g), which allowed a diamond only when I2 heterogeneity was "
         "non-significant.",
         "Figure 3. Both cohorts' estimates are shown separately; no between-cohort "
         "heterogeneity statistic is reported because no pooled estimate is drawn."],
        ["Signed direction consistency restricted to the discovery-evaluable types",
         "The pre-specified >=7/10 direction-agreement test was retired; the GLASS "
         "direction comparison is descriptive and is limited to the intersection of "
         "cell types evaluable in discovery (NK, Microglia).",
         "The plan-literal test cannot be evaluated when only two types are evaluable, "
         "and a 2-cell agreement count cannot support a claim of agreement.",
         "Both evaluable types flip sign relative to GLASS. No agreement count is "
         "reported anywhere in this package."],
        ["Pre-declared immune-subclustering sensitivity analysis",
         "The immune compartment was re-clustered at resolutions 0.5, 1.0 and 1.5 and "
         "evaluability re-assessed under the same >=10-sample rule.",
         "Pre-declared to test whether the annotation step, rather than the data, "
         "limited which cell types could be resolved.",
         "The gate failed (%d/%d/%d of 10 evaluable); the primary annotation is "
         "unchanged and no mechanistic analysis followed."
         % (g("sensitivity.evaluable_counts")["0.5"],
            g("sensitivity.evaluable_counts")["1.0"],
            g("sensitivity.evaluable_counts")["1.5"])],
    ]
    add("Table6", "Deviations from pre-specification v1.0.0",
        ["Item", "What changed", "Reason", "Affected results"], rows,
        "Each deviation was verified against the named script or audit artefact. Where "
        "a claim could not be confirmed from project files it is marked [VERIFY].")


# =============================================================== Table S1
def ts1():
    pt = {r["cell_type"]: r for r in rd_csv(r"outputs\intermediate\post_data_power_table.csv")}
    rows = []
    for s in PRIMARY_SHORT:
        r = pt[s]
        if r["has_observed_fit"] == "True":
            f2o = f3(r["f2_obs_icc60"])
            pfo = "%.1f%%" % (100 * float(r["power_f2obs_icc60"]))
            tag = r["headline_label_icc60"]
        else:
            f2o, pfo, tag = "", "n/a (no observed fit)", "n/a (not estimable)"
        rows.append([DISP[s], f2o, pfo,
                     "%.1f%%" % (100 * float(r["power_f2028_icc60"])),
                     "%.1f%%" % (100 * float(r["power_f2040_icc60"])), tag])
    add("TableS1", "Post-data power (NSIM = %d, ICC = 0.60)" % g("power.nsim"),
        ["Cell type", "Observed f2", "Power at observed f2",
         "Power at f2 = 0.28", "Power at f2 = 0.40", "Headline label"], rows,
        "Simulated with seed %d at the recovered design (%d patients; %d Any_RT, "
        "%d No_RT). f2-anchored power is a property of the design and is therefore "
        "identical across cell types."
        % (g("power.seed"), g("power.n_patients"), g("power.design_anyrt"),
           g("power.design_nort")))


# =============================================================== Table S2
def ts2():
    st = rd_csv(r"outputs\intermediate\cgga_crosssectional_stats.csv")
    rows = []
    for r in st:
        rows.append([r["cell_family"],
                     "primary" if r["is_primary"] == "True" else "exploratory",
                     r["method"], r["n_anyrt"], r["n_nort"],
                     f3(r["delta_median_anyrt_minus_nort"]),
                     f3(r["mannwhitney_p_2tail"]),
                     f3(r["mannwhitney_q_primary_bh"])])
    add("TableS2", "CGGA cross-sectional statistics (EXPLORATORY)",
        ["Cell type", "Family", "Method", "Any_RT (n)", "No_RT (n)",
         "Delta median (Any_RT - No_RT)", "Mann-Whitney p (2-tailed)",
         "BH q over primary types"], rows,
        "EXPLORATORY cross-sectional comparison; no matched primary tumour, so no "
        "interaction can be estimated. MCP-counter returns abundance scores rather "
        "than fractions and is descriptive only. BH q applies to the primary family "
        "within each method.")


# =============================================================== Table S3
def ts3():
    ea = json.load(io.open(os.path.join(
        INT, "immune_subclust_audit_immune_subclust.json"), "r", encoding="utf-8"))
    rows = []
    for res in ("0.5", "1.0", "1.5"):
        for c in ea["clusters_per_resolution"][res]:
            rows.append([res, c["cluster"], c["n_cells"], c.get("winner", ""),
                         f3(c.get("winner_score")), c.get("second", ""),
                         f3(c.get("second_score")), f3(c.get("margin")),
                         f3(c.get("ratio")), c.get("n_de_up_markers", ""),
                         "yes" if c.get("passes") else "no",
                         ";".join(c.get("flags", []))])
    add("TableS3", "Per-cluster diagnostics of the immune-subclustering sensitivity analysis",
        ["Resolution", "Cluster", "Cells (n)", "Argmax winner", "Winner score",
         "Second label", "Second score", "Margin", "Second/winner ratio",
         "DE up-markers (n)", "Passes all criteria", "Flags"], rows,
        "SENSITIVITY. All clusters at all three pre-declared resolutions are listed "
        "so the gate outcome can be reproduced.")


# =============================================================== Table S4
def ts4():
    """Marker dictionary parsed from source with ast - the pipeline is never imported."""
    src = os.path.join(ROOT, "scripts", "acquire_snrnaseq_discovery.py")
    with io.open(src, "r", encoding="utf-8") as f:
        tree = ast.parse(f.read(), filename=src)
    markers = None
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            for t in node.targets:
                if isinstance(t, ast.Name) and t.id == "MARKERS":
                    markers = ast.literal_eval(node.value)
    if markers is None:
        raise SystemExit("[!!] MARKERS assignment not found by ast parse")
    rows = [[k, len(v), ", ".join(v)] for k, v in markers.items()]
    add("TableS4", "Marker gene sets used for cluster-level annotation",
        ["Cell type / population", "Genes (n)", "Genes"], rows,
        "Read from scripts/acquire_snrnaseq_discovery.py with the ast module (the "
        "module is never imported). Genes absent from the 2,000-gene highly variable "
        "space are force-included before annotation.")


# =============================================================== Table S5
def ts5():
    try:
        out = subprocess.run([sys.executable, "-m", "pip", "freeze"],
                             capture_output=True, text=True, timeout=180)
        lines = [l for l in out.stdout.splitlines() if l.strip()]
    except Exception as e:  # pragma: no cover
        lines = ["[VERIFY] pip freeze failed: %r" % e]
    rows = []
    for l in lines:
        if "==" in l:
            n, v = l.split("==", 1)
            rows.append([n, v])
        else:
            rows.append([l, ""])
    rows.append(["python", sys.version.split()[0]])
    add("TableS5", "Software versions",
        ["Package", "Version"], rows,
        "Read-only listing of the project virtual environment. Python %s."
        % sys.version.split()[0])


# =============================================================== Table S6
def ts6():
    rows = [
        ["Version", g("prespec.version")],
        ["Locked at (UTC)", g("prespec.locked_at")],
        ["Integrity sha256 (sec4_prespec.json)",
         g("prespec.sha256_json")[:16] + " (prefix; full hash in the copied file)"],
        ["Integrity sha256 (sec4_prespec.md)", g("prespec.sha256_md")[:16]],
        ["Primary cell types (BH-FDR family K)", "10"],
        ["Exploratory cell types (raw p only)", "5"],
        ["Response transform", "arcsin_sqrt of the cell-type fraction"],
        ["Hypothesis term", "timepoint x treatment_arm (binary Any_RT vs No_RT)"],
        ["Random effects", "(1 | patient_id)"],
        ["Permutation null", "10,000 patient-level treatment-label shuffles, seed %d"
         % g("power.seed")],
        ["Discovery arm sizes", "Any_RT %d patients; No_RT %d patients"
         % (g("cohorts.discovery.n_any_rt_patients"), g("cohorts.discovery.n_no_rt_patients"))],
        ["Forest plot (section 4.5(g))", "profile-likelihood 95% CI; pooled diamond "
         "only if I2 p > 0.10 - superseded, see Table 6"],
    ]
    add("TableS6", "Pre-specification v1.0.0 summary",
        ["Item", "Value"], rows,
        "The locked pre-specification is copied unmodified into "
        "submission/supplementary/. Deviations are listed in Table 6.")


def write_index():
    p = os.path.join(TAB, "tables_index.csv")
    with io.open(p, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["table_id", "title", "n_rows", "csv", "xlsx"])
        for tid, obj in TABLES.items():
            w.writerow([tid, obj["title"], len(obj["rows"]),
                        tid + ".csv", tid + ".xlsx"])
    return p


def main():
    for fn in (t1, t2, t3, t4, t5, t6, ts1, ts2, ts3, ts4, ts5, ts6):
        fn()
    for tid, obj in TABLES.items():
        write_csv(tid, obj)
        write_xlsx(tid, obj)
    print("index written to", write_index())
    print("tables written to", TAB)
    for tid in TABLES:
        print("  %-9s %-62s rows=%d" % (tid, TABLES[tid]["title"][:62],
                                        len(TABLES[tid]["rows"])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
