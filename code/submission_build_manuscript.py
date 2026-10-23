# -*- coding: utf-8 -*-
"""Stage 4 - manuscript for the DGNC submission package.

Builds submission/manuscript/manuscript.docx and manuscript.html from ONE
content structure, so the two renderings cannot diverge.

Every number is read from submission/supplementary/results_numbers.json
(Stage 1); table contents are read from submission/tables/*.csv (Stage 3);
figure legends are read from submission/figures/figure_legends.txt (Stage 2).
Nothing is re-analysed and no pipeline script is imported.
"""
from __future__ import annotations

import csv
import io
import json
import os
import re

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor

ROOT = r"<PROJECT_ROOT>"
SUB = os.path.join(ROOT, "submission")
TAB = os.path.join(SUB, "tables")
FIG = os.path.join(SUB, "figures")
MAN = os.path.join(SUB, "manuscript")
os.makedirs(MAN, exist_ok=True)

TITLE_OPTIONS = [
    "Limits of single-nucleus data for detecting radiotherapy-associated immune "
    "remodeling at recurrence in IDH-mutant glioma: a pre-specified analysis with "
    "exploratory external cohorts",
    "Two of ten: a pre-specified negative-discovery analysis of radiotherapy-associated "
    "immune change at IDH-mutant glioma recurrence",
    "What single-nucleus profiling of paired IDH-mutant glioma samples can and cannot "
    "test about radiotherapy-associated immune remodeling: a pre-specified negative report",
]
TITLE = TITLE_OPTIONS[0]
RUNNING_TITLE = ("Limits of single-nucleus data for a radiotherapy interaction in "
                 "IDH-mutant glioma")

BANNED = [r"\brobust\w*", r"\bas planned\b", r"\bmandatory\b", r"\breplicat\w*",
          r"\bvalidat\w*", r"sign[- ]agreement", r"5 of 10 estimable",
          r"1 of 5 agree", r"\bcaused\b", r"\bcauses\b", r"\bled to\b",
          r"\bdrives\b", r"\bresolves\b", r"\bresolve\b"]


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


NUM = _nest(json.load(io.open(os.path.join(SUB, "supplementary", "results_numbers.json"),
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


def n0(dotted):
    return "%d" % int(round(float(g(dotted))))


def n3(dotted):
    return "%.3f" % float(g(dotted))


def pct(dotted):
    return "%.1f%%" % (100.0 * float(g(dotted)))


# ------------------------------------------------------------------ inputs
def read_table(tid):
    with io.open(os.path.join(TAB, tid + ".csv"), "r", encoding="utf-8",
                 newline="") as f:
        rows = list(csv.reader(f))
    header, body, notes = rows[0], [], ""
    for r in rows[1:]:
        if not r:
            continue
        if len(r) >= 2 and r[0] == "note":
            notes = r[1]
            continue
        body.append(r)
    return {"header": header, "rows": body, "notes": notes}


TAB_TITLES = {}
with io.open(os.path.join(TAB, "tables_index.csv"), "r", encoding="utf-8",
             newline="") as f:
    for r in csv.DictReader(f):
        TAB_TITLES[r["table_id"]] = r["title"]


def read_legends():
    txt = io.open(os.path.join(FIG, "figure_legends.txt"), "r",
                  encoding="utf-8").read()
    out = {}
    for block in re.split(r"\n-{5,}\n", txt):
        lines = [l.strip() for l in block.strip().split("\n")]
        if not lines or not lines[0]:
            continue
        m = re.match(r"^(Figure S?\d+)\.\s+(.*)$", lines[0])
        if not m:
            continue
        key = m.group(1).replace(" ", "")
        paras = [p for p in re.split(r"\n\s*\n", "\n".join(lines[1:])) if p.strip()]
        alt = ""
        legend_parts = []
        for p in paras:
            p = " ".join(p.split())
            if p.lower().startswith("alt text:"):
                alt = p.split(":", 1)[1].strip()
            else:
                legend_parts.append(p)
        out[key] = {"title": m.group(2).strip(),
                    "legend": " ".join(legend_parts),
                    "alt": alt}
    return out


LEG = read_legends()


# ------------------------------------------------------------------ numbers
D_NS = n0("cohorts.discovery.n_samples")
D_NP = n0("cohorts.discovery.n_patients")
D_ANY = n0("cohorts.discovery.n_any_rt_patients")
D_NO = n0("cohorts.discovery.n_no_rt_patients")
D_IN = n0("cohorts.discovery.samples_input_snrna")
D_EXCL = n0("cohorts.discovery.n_excluded_patients")
D_GENES = "{:,}".format(int(g("cohorts.discovery.n_genes_full")))
D_HVG = "{:,}".format(int(g("cohorts.discovery.n_hvg")))
D_CAP = "{:,}".format(int(g("cohorts.discovery.cap_cells_per_sample")))
D_NCAP = n0("cohorts.discovery.n_samples_capped")
D_PRECAP = "{:,}".format(int(g("cohorts.discovery.cells_precap_total")))
D_POSTCAP = "{:,}".format(int(g("cohorts.discovery.cells_final_postcap")))

EV_TOT = n0("evaluability.n_primary_types")
EV_N = n0("evaluability.n_primary_evaluable")
EV_THR = n0("evaluability.threshold_samples")
EV_SAMP = n0("evaluability.n_samples_with_proportions")
NK_S = n0("evaluability.by_primary_cell_type.NK.present_samples")
MG_S = n0("evaluability.by_primary_cell_type.Microglia_P2RY12_TMEM119.present_samples")
B_S = n0("evaluability.by_primary_cell_type.B.present_samples")
MOMF_S = n0("evaluability.by_primary_cell_type.MoMF_CD163.present_samples")
CD8_S = n0("evaluability.by_primary_cell_type.CD8_T.present_samples")

NK_B, NK_L, NK_H = n3("discovery.lmm_by_cell.NK.beta"), \
    n3("discovery.lmm_by_cell.NK.ci_low"), n3("discovery.lmm_by_cell.NK.ci_high")
NK_P, NK_Q = n3("discovery.lmm_by_cell.NK.asym_p"), n3("discovery.lmm_by_cell.NK.bh_q")
NK_PERM = n3("discovery.lmm_by_cell.NK.perm_p")
MG_B, MG_L, MG_H = n3("discovery.lmm_by_cell.Microglia_P2RY12_TMEM119.beta"), \
    n3("discovery.lmm_by_cell.Microglia_P2RY12_TMEM119.ci_low"), \
    n3("discovery.lmm_by_cell.Microglia_P2RY12_TMEM119.ci_high")
MG_P = n3("discovery.lmm_by_cell.Microglia_P2RY12_TMEM119.asym_p")
MG_Q = n3("discovery.lmm_by_cell.Microglia_P2RY12_TMEM119.bh_q")
MG_PERM = n3("discovery.lmm_by_cell.Microglia_P2RY12_TMEM119.perm_p")

PW_NSIM = n0("power.nsim")
PW_F2A = "%.2f" % float(g("power.f2_anchors")[0])
PW_F2B = "%.2f" % float(g("power.f2_anchors")[1])
PW_PA = pct("power.headline_f2028_icc60")
PW_PB = pct("power.headline_f2040_icc60")
PW_LAB = str(g("power.headline_label"))

GL_NS = n0("cohorts.glass.n_samples")
GL_NC = n0("cohorts.glass.n_cases")
GL_ANY = n0("cohorts.glass.anyrt_pairs")
GL_NO = n0("cohorts.glass.nort_pairs")
GL_EST = n0("glass.n_primary_estimable")
GL_BE, GL_BQ = n3("glass.lmm_by_cell.B.beta"), n3("glass.lmm_by_cell.B.bh_q")
GL_ME = n3("glass.lmm_by_cell.Microglia_P2RY12_TMEM119.beta")
GL_MQ = n3("glass.lmm_by_cell.Microglia_P2RY12_TMEM119.bh_q")
GL_NK = n3("glass.lmm_by_cell.NK.beta")
GL_SIG = n0("glass.n_cells_bh_q_lt_0.05")

CG_CLN = n0("cohorts.cgga.samples_clinical_total")
CG_IDH = n0("cohorts.cgga.idh_mut_recurrent_total")
CG_N = n0("cohorts.cgga.n_analysed")
CG_ANY = n0("cohorts.cgga.anyrt_samples")
CG_NO = n0("cohorts.cgga.nort_samples")
CG_M = n0("cohorts.cgga.n_methods")
CG_Q = n3("cgga.best_primary_mw_q")
CG_P = n3("cgga.best_primary_mw_p")
CG_EQ = n3("cgga.endothelial_by_method.quantiseq.mw_p")
CG_EE = n3("cgga.endothelial_by_method.epic.mw_p")
CG_EM = n3("cgga.endothelial_by_method.mcpcounter.mw_p")

SE_C = "{:,}".format(int(g("sensitivity.compartment_cells")))
SE_05 = n0("sensitivity.evaluable_counts.0.5")
SE_10 = n0("sensitivity.evaluable_counts.1.0")
SE_15 = n0("sensitivity.evaluable_counts.1.5")
SE_C05 = "%d of %d" % (g("sensitivity.contamination_by_resolution.0.5.non_immune_clusters"),
                       g("sensitivity.contamination_by_resolution.0.5.total_clusters"))
SE_C10 = "%d of %d" % (g("sensitivity.contamination_by_resolution.1.0.non_immune_clusters"),
                       g("sensitivity.contamination_by_resolution.1.0.total_clusters"))
SE_C15 = "%d of %d" % (g("sensitivity.contamination_by_resolution.1.5.non_immune_clusters"),
                       g("sensitivity.contamination_by_resolution.1.5.total_clusters"))

PS_VER = str(g("prespec.version"))
PS_LOCK = str(g("prespec.locked_at"))[:10]
PS_SHA = str(g("prespec.sha256_json"))[:16]

SW = {}
with io.open(os.path.join(TAB, "TableS5.csv"), "r", encoding="utf-8",
             newline="") as f:
    for r in list(csv.reader(f))[1:]:
        if len(r) >= 2 and r[0] not in ("note",):
            SW[r[0]] = r[1]


def sw(name):
    return SW.get(name, "[VERIFY]")


# ------------------------------------------------------------------ abstract
ABSTRACT = (
    f"Radiotherapy exposure is a candidate moderator of the change in immune cell "
    f"composition between primary and recurrent IDH-mutant glioma, and a pre-specified "
    f"test of that interaction has not been reported. "
    f"We fitted patient-random-intercept mixed models of arcsin-square-root-transformed "
    f"cell-type fractions with a timepoint by treatment-arm interaction to paired "
    f"single-nucleus RNA sequencing from {D_NS} samples from {D_NP} patients "
    f"({D_ANY} with any radiotherapy exposure, {D_NO} with none), under a pre-declared "
    f"evaluability rule of a non-zero proportion in at least {EV_THR} of {EV_SAMP} "
    f"samples, with a 10,000-shuffle patient-level permutation null, a post-data power "
    f"calculation, two external cohorts designated exploratory, and a pre-declared "
    f"immune-subclustering sensitivity analysis. "
    f"Only {EV_N} of the {EV_TOT} pre-specified primary cell types met the evaluability "
    f"rule: NK cells (interaction beta {NK_B}, 95% confidence interval {NK_L} to {NK_H}, "
    f"Benjamini-Hochberg q {NK_Q}, permutation p {NK_PERM}) and microglia (beta {MG_B}, "
    f"95% confidence interval {MG_L} to {MG_H}, q {MG_Q}, permutation p {MG_PERM}), "
    f"while post-data simulation at an intraclass correlation of 0.60 gave {PW_PA} power "
    f"at f2 = {PW_F2A} and {PW_PB} power at f2 = {PW_F2B}, which the power audit labels "
    f"{PW_LAB}. "
    f"In the exploratory GLASS bulk deconvolution arm ({GL_NS} samples, {GL_NC} cases) "
    f"the interaction was estimable for all {GL_EST} primary types and reached q < 0.05 "
    f"for B cells and microglia, and both discovery-evaluable types change direction "
    f"between cohorts, so this arm is reported as exploratory and is not treated as "
    f"evidence about the discovery estimates. "
    f"The single-nucleus design therefore neither supports nor excludes a "
    f"radiotherapy-associated immune shift at recurrence, and a design able to answer "
    f"the question would need a larger radiotherapy-stratified paired cohort with "
    f"annotation that separates the T and NK lineages."
)

KEYWORDS = ("IDH-mutant glioma; glioma recurrence; radiotherapy; single-nucleus RNA "
            "sequencing; immune cell composition; pre-specified analysis; negative "
            "results; deconvolution")


# ------------------------------------------------------------------ content
def build_blocks():
    """One content structure rendered to docx and html."""
    B = {"title": [], "abstract": [], "main": [], "back": [], "supp": []}

    B["title"] += [
        ("title", TITLE),
        ("p", "Running title: " + RUNNING_TITLE),
        ("p", "Authors: [AUTHOR LIST]"),
        ("p", "Affiliations: [AFFILIATIONS]"),
        ("p", "Corresponding author: [CORRESPONDING AUTHOR] [EMAIL]"),
        ("p", "Funding: [FUNDING]"),
        ("p", "Competing interests: [CONFLICTS]"),
        ("p", "Ethics approval and consent to participate: [ETHICS]"),
        ("p", "This is a pre-specified negative-discovery report: it makes no causal "
              "claim about radiotherapy and no claim of agreement between cohorts."),
        ("pb", ""),
    ]

    B["abstract"] += [
        ("h1", "Abstract"),
        ("p", ABSTRACT),
        ("p", "Keywords: " + KEYWORDS),
    ]

    M = B["main"]
    M += [
        ("h1", "Introduction"),
        ("p", "IDH-mutant gliomas recur, and whether radiotherapy given before "
              "recurrence leaves a measurable imprint on the immune composition of the "
              "recurrent tumour is unresolved. The informative question is not whether "
              "irradiated and non-irradiated tumours differ in composition, because such "
              "a comparison is confounded by indication and by time, but whether "
              "radiotherapy exposure moderates the change in composition between the "
              "primary and the recurrent sample of the same patient."),
        ("p", "The evidence base for that interaction is thin. A longitudinal analysis of "
              "60 IDH-mutant whole-transcriptome pairs reported a pooled change towards "
              "an immune-activated phenotype at recurrence that was not significant and "
              "was not stratified by treatment (Tang et al., Cancer Cell; full citation "
              "details [VERIFY CITATION]). A matched-pair histological series of 27 "
              "radiotherapy-only recurrent pairs reported a change in lymphocyte density "
              "(Makarevic et al., Int J Mol Sci 2020, doi:10.3390/ijms21207801; PMID "
              "33096928), but histology does not separate lymphoid lineages and cannot "
              "adjust for tumour purity. Neither study tested a treatment-by-timepoint "
              "interaction."),
        ("p", "Longitudinal single-nucleus RNA sequencing with treatment annotation now "
              "makes that test possible in principle (Johnson et al., Nature 2026, "
              "doi:10.1038/s41586-026-10612-6; GSE326221). We therefore wrote and locked "
              f"a pre-specification (version {PS_VER}, locked {PS_LOCK}, integrity hash "
              f"prefix {PS_SHA}, summarised in Table S6) before inspecting the "
              "interaction estimates. It fixed the cell types to be tested, an "
              "evaluability rule, the model, a permutation null, a post-data power "
              "calculation and a sensitivity analysis, and it fixed the report as a "
              "negative-discovery analysis of what the design can and cannot test."),
    ]

    M += [
        ("h1", "Methods"),
        ("h2", "Cohorts and samples"),
        ("p", f"The discovery arm was the open single-nucleus RNA sequencing arm of the "
              f"Johnson/CARE longitudinal glioma cohort (Nature 2026, "
              f"doi:10.1038/s41586-026-10612-6). {D_IN} samples were loaded, and after "
              f"joining harmonised clinical covariates the modelled design contained "
              f"{D_NS} samples from {D_NP} patients, {D_ANY} with any radiotherapy "
              f"exposure and {D_NO} with none; {D_EXCL} patients were excluded because "
              f"treatment exposure was not recorded. The unit of analysis was the "
              f"per-sample proportion of cells assigned to each cell type."),
        ("p", f"The first external arm was GLASS, used as bulk RNA sequencing: {GL_NS} "
              f"paired samples from {GL_NC} cases, comprising {GL_ANY} pairs with any "
              f"radiotherapy exposure and {GL_NO} without (Synapse project "
              f"syn17038081). The second external arm was a cross-sectional sample of "
              f"the CGGA mRNAseq_693 release (20200506): {CG_CLN} clinical records, of "
              f"which {CG_IDH} IDH-mutant recurrent tumours were identified and {CG_N} "
              f"had expression data ({CG_ANY} any radiotherapy, {CG_NO} none)."),
        ("h2", "Single-nucleus processing and cell-type annotation"),
        ("p", f"Cell-level quality control (at least 300 genes per cell) was applied on "
              f"the full {D_GENES}-gene expression space before selecting the {D_HVG} "
              f"most highly variable genes, because applying the same thresholds after "
              f"gene selection removed approximately four times as many cells as "
              f"intended. Each sample was then uniformly subsampled to at most {D_CAP} "
              f"cells ({D_NCAP} of {D_IN} samples capped), reducing {D_PRECAP} cells to "
              f"{D_POSTCAP} (Table 6)."),
        ("p", "Cell types were assigned by the argmax of the mean marker score of each "
              "Leiden cluster at resolution 0.8, so that every cell in a cluster "
              "receives that cluster's label. The marker sets are listed in Table S4 "
              "and cover the 10 pre-specified primary types plus five exploratory "
              "types. Because NK and CD8 T cells share markers, NK cells can be absorbed "
              "into CD8-T-labelled clusters; the consequence for evaluability is "
              "reported under Results."),
        ("h2", "Evaluability rule"),
        ("p", f"A pre-specified primary cell type was evaluable only if it had a "
              f"non-zero proportion in at least {EV_THR} of the {EV_SAMP} samples with "
              f"proportions (Table 2). The rule was fixed before the interaction "
              f"estimates were seen and was applied unchanged to the sensitivity "
              f"analysis. Estimator output that falls below the rule is retained in the "
              f"tables but tagged NOT EVALUABLE; descriptive only, and no conclusion is "
              f"drawn from it."),
        ("h2", "Interaction mixed model"),
        ("p", "For each cell type the per-sample fraction was transformed with "
              "arcsin(sqrt(fraction)) and fitted as arcsin_sqrt_fraction ~ timepoint * "
              "treatment_arm + purity + WHO grade at diagnosis + age at diagnosis + "
              "1p19q status + (1 | patient), with treatment coded binary as any "
              "radiotherapy exposure versus none. The hypothesis term was "
              "timepoint:treatment_arm. Confidence intervals are profile-likelihood "
              "intervals with a Wald fallback if the profile search did not converge. "
              "Where an optimiser failed, the fit was retried with successive optimisers, "
              "and covariates with zero variance in a given arm were dropped from that "
              "arm's effective formula (Table 6 records the optimiser chain and the "
              "dropped terms)."),
        ("h2", "Permutation null"),
        ("p", f"To avoid relying on the asymptotic distribution of the interaction "
              f"coefficient, its null distribution was generated by 10,000 patient-level "
              f"shuffles of the treatment label with the two timepoints of a patient kept "
              f"together, refitting the same model each time (seed "
              f"{n0('power.seed')}). Permutation p values are reported alongside the "
              f"asymptotic p values for every cell type with a discovery estimate."),
        ("h2", "Post-data power"),
        ("p", f"Power was estimated after the data were seen, by Monte-Carlo simulation "
              f"with {PW_NSIM} simulation draws each at the recovered design ({D_NP} "
              f"patients; {D_ANY} any radiotherapy and {D_NO} none), for assumed "
              f"intraclass correlations of 0.40, 0.60 and 0.70 and effect sizes anchored "
              f"at f2 = {PW_F2A} and f2 = {PW_F2B}, with the patient-level correlation "
              f"included "
              f"(seed {n0('power.seed')}). This is a ceiling on what the recovered design "
              f"could have detected, not a pre-study sample-size calculation."),
        ("h2", "Exploratory external arms"),
        ("p", "In GLASS, cell-type fractions were estimated by constrained least-squares "
              "deconvolution with a stock published signature matrix, because the cohort "
              "has no matched single-nucleus reference from which a custom signature "
              "could be built. The matrix was extended with canonical plasma-cell and "
              "microglia rows, since the stock marker set models neither population. "
              "Purity was not identifiable in this cohort and was therefore not adjusted "
              "for. All GLASS estimates are labelled EXPLORATORY."),
        ("p", f"Cell-type abundance in the CGGA arm was estimated with three stock "
              f"methods (quanTIseq, EPIC and MCP-counter) and compared cross-sectionally "
              f"between tumours with and without radiotherapy exposure using two-tailed "
              f"Mann-Whitney tests; MCP-counter returns abundance scores rather than "
              f"fractions and is reported descriptively. This arm has no second "
              f"timepoint, so it cannot address the interaction; all CGGA results are "
              f"labelled EXPLORATORY."),
        ("h2", "Pre-declared immune-subclustering sensitivity analysis"),
        ("p", f"The sensitivity analysis re-clustered the annotated immune compartment "
              f"({SE_C} cells) at Leiden resolutions 0.5, 1.0 and 1.5 and re-applied the "
              f"same evaluability rule. Its pre-declared gate required at least 6 of 10 "
              f"primary cell types evaluable at resolution 1.0 and at least 5 of 10 at "
              f"both 0.5 and 1.5; the outcome is reported under Results (Table 5, "
              f"Table S3)."),
        ("h2", "Statistics, multiplicity and reporting labels"),
        ("p", f"Benjamini-Hochberg correction was applied across the {EV_TOT} primary "
              f"cell types (K = {EV_TOT}); exploratory types are reported with raw p "
              f"values only, and the exploratory arms were not corrected. Every result in "
              f"this report carries one of four labels: PRIMARY (discovery arm, "
              f"pre-specified, evaluable), EXPLORATORY (external cohort designated "
              f"exploratory in the pre-specification), DIRECTION-ONLY (descriptive "
              f"comparison, no test) or SENSITIVITY (pre-declared sensitivity analysis). "
              f"Figures were produced with matplotlib {sw('matplotlib')}, tables and "
              f"counts with pandas {sw('pandas')} and numpy {sw('numpy')}, "
              f"single-nucleus processing with scanpy {sw('scanpy')} and annData "
              f"{sw('anndata')}, and mixed models with statsmodels {sw('statsmodels')}; "
              f"the full frozen dependency list is Table S5."),
        ("h2", "Deviations from the pre-specification (v1.0.0)"),
        ("p", f"The analysis deviates from the locked pre-specification ({PS_VER}, locked "
              f"{PS_LOCK}, integrity hash prefix {PS_SHA}) in the ways listed in Table 6, "
              f"each with its reason and the results it affects: MGMT methylation dropped "
              f"as a covariate; 1p19q status entered as an IDH-subtype proxy; purity "
              f"derived from the single-nucleus malignant-cell fraction; cell-level "
              f"quality control applied on the full gene space; a per-sample cap of "
              f"{D_CAP} cells; cluster-level argmax annotation; optimiser fallback with "
              f"dropping of zero-variance covariates; post-data power simulated with "
              f"{PW_NSIM} simulation draws; GLASS deconvolution with stock signatures "
              f"and no purity adjustment; forest-plot markers filled by BH q with no "
              f"pooled estimate, deviating from section 4.5(g) of the pre-specification, "
              f"which permitted a pooled diamond when between-cohort heterogeneity was "
              f"not significant; the signed direction-consistency count retired because "
              f"two cell types are evaluable and a two-cell count cannot support a claim "
              f"of agreement; and the pre-declared immune-subclustering sensitivity "
              f"analysis."),
    ]

    M += [
        ("h1", "Results"),
        ("h2", "Design and evaluability [PRIMARY]"),
        ("p", f"The discovery arm loaded {D_IN} samples and supported a modelled design of "
              f"{D_NS} samples from {D_NP} patients ({D_ANY} with any radiotherapy "
              f"exposure, {D_NO} with none) and {D_POSTCAP} cells (Figure 1, Table 1). Of "
              f"the {EV_TOT} pre-specified primary cell types, five received at least one "
              f"cell in the annotation, and five (CD4 conventional T cells, regulatory T "
              f"cells, plasma cells, dendritic cells and neutrophils) formed no distinct "
              f"cluster and therefore have no proportion column at all (Table 2). Under "
              f"the pre-declared rule, only NK cells ({NK_S} of {EV_SAMP} samples) and "
              f"microglia ({MG_S} of {EV_SAMP} samples) were evaluable; B cells ({B_S}), "
              f"monocyte-derived macrophages ({MOMF_S}) and CD8 T cells ({CD8_S}) fell "
              f"below the threshold of {EV_THR} samples (Figure 2, Table 2). "
              f"{EV_N} of {EV_TOT} primary cell types are therefore evaluable. This is a "
              f"limit of marker-score annotation rather than evidence that the "
              f"populations are missing: NK and CD8 T cells share markers, and NK cells "
              f"were absorbed into CD8-T-labelled clusters that themselves failed the "
              f"cluster margin and marker-count criteria."),
        ("fig", "Figure1"),
        ("h2", "Discovery interaction estimates [PRIMARY]"),
        ("p", f"For NK cells the timepoint by treatment interaction was beta = {NK_B} "
              f"(profile-likelihood 95% confidence interval {NK_L} to {NK_H}; asymptotic "
              f"p = {NK_P}; BH q = {NK_Q}; permutation p = {NK_PERM}). For microglia it "
              f"was beta = {MG_B} (95% confidence interval {MG_L} to {MG_H}; p = {MG_P}; "
              f"q = {MG_Q}; permutation p = {MG_PERM}) (Table 3, Figure 3 panel A). "
              f"Neither coefficient is distinguished from its permutation null by the "
              f"asymptotic or the permutation test, and both confidence intervals include "
              f"zero with appreciable width in both directions. The interaction is "
              f"unresolved in both evaluable cell types, and no radiotherapy effect is "
              f"claimed. Of the remaining eight primary types, five have no discovery "
              f"estimate at all and three (B cells, monocyte-derived macrophages and CD8 "
              f"T cells) carry finite model estimates that fall below the evaluability "
              f"threshold; those three are tagged NOT EVALUABLE; descriptive only in "
              f"Table 3."),
        ("h2", "Power of the recovered design [PRIMARY]"),
        ("p", f"At an assumed intraclass correlation of 0.60, the recovered design "
              f"reaches {PW_PA} power at f2 = {PW_F2A} and {PW_PB} power at f2 = "
              f"{PW_F2B} with {PW_NSIM} simulation draws each. The power audit labels "
              f"the design {PW_LAB} at the anchor effect sizes used in the "
              f"pre-specification (Table S1, Figure S1). A null interaction in this "
              f"design therefore does not exclude an effect of the size the anchors "
              f"describe."),
        ("h2", "Exploratory external cohort (GLASS) [EXPLORATORY]"),
        ("p", f"In the GLASS arm ({GL_NS} samples from {GL_NC} cases; {GL_ANY} pairs with "
              f"any radiotherapy exposure and {GL_NO} without) the interaction was "
              f"estimable for all {GL_EST} primary cell types. BH q < 0.05 was reached "
              f"for B cells (beta = {GL_BE}, q = {GL_BQ}) and microglia (beta = {GL_ME}, "
              f"q = {GL_MQ}), which is {GL_SIG} of {GL_EST} primary types (Table 4, "
              f"Figure 3 panel B). These estimates are EXPLORATORY: the fractions come "
              f"from stock deconvolution signatures without a matched single-nucleus "
              f"reference and without purity adjustment, and the group without "
              f"radiotherapy exposure contains only {GL_NO} pairs (Figure S3)."),
        ("h2", "Direction comparison between cohorts [DIRECTION-ONLY]"),
        ("p", f"Descriptively, both cell types that are evaluable in the discovery arm "
              f"change direction between the two cohorts: NK cells beta = {NK_B} in "
              f"discovery and {GL_NK} in GLASS, and microglia beta = {MG_B} in discovery "
              f"and {GL_ME} in GLASS. No count of agreeing cell types is computed, no "
              f"pooled estimate across cohorts is reported, and no agreement between "
              f"cohorts is claimed; this comparison is DIRECTION-ONLY and the difference "
              f"is expected when both estimates are compatible with zero and the cohorts "
              f"differ in assay, deconvolution method and treatment coding."),
        ("h2", "Exploratory cross-sectional cohort (CGGA) [EXPLORATORY]"),
        ("p", f"In the CGGA arm ({CG_N} IDH-mutant recurrent tumours; {CG_ANY} with any "
              f"radiotherapy exposure and {CG_NO} without) no primary cell type reached "
              f"BH q < 0.05 with any of the {CG_M} deconvolution methods (best primary "
              f"q = {CG_Q}, corresponding to p = {CG_P} in the method where it was "
              f"smallest; Table S2). Endothelial cells, an exploratory type, differed "
              f"between exposure groups in the same direction with all {CG_M} methods "
              f"(quanTIseq p = {CG_EQ}, EPIC p = {CG_EE}, MCP-counter p = {CG_EM}), and "
              f"this is reported as EXPLORATORY only (Figure S2). The arm is "
              f"cross-sectional and cannot address the timepoint by treatment "
              f"interaction."),
        ("h2", "Pre-declared immune-subclustering sensitivity analysis [SENSITIVITY]"),
        ("p", f"The immune compartment contained {SE_C} cells. Under the same evaluability "
              f"rule, {SE_05}, {SE_10} and {SE_15} of the {EV_TOT} primary types were "
              f"evaluable at resolutions 0.5, 1.0 and 1.5, against a gate requiring at "
              f"least 6 of 10 at 1.0 and at least 5 of 10 at both 0.5 and 1.5, so the "
              f"gate failed (Table 5, Figure 4). Microglia were evaluable at all three "
              f"resolutions and monocyte-derived macrophages at 1.0 and 1.5; the other "
              f"eight types were not. Subclusters whose argmax label was a non-immune "
              f"cell type accounted for {SE_C05}, {SE_C10} and {SE_C15} of the "
              f"subclusters at the three resolutions, and PTPRC was absent from the "
              f"{D_HVG}-gene highly variable space used for this step, so immune cells "
              f"sitting in non-immune clusters could not be recovered. The sensitivity "
              f"analysis does not replace the primary annotation and no mechanistic "
              f"analysis was performed on it."),
        ("fig", "Figure2"),
        ("fig", "Figure3"),
        ("fig", "Figure4"),
    ]

    M += [
        ("h1", "Discussion"),
        ("p", f"The pre-specified question was whether radiotherapy exposure moderates "
              f"the change in immune cell composition between primary and recurrent "
              f"IDH-mutant glioma. In this single-nucleus cohort it could be tested for "
              f"{EV_N} of {EV_TOT} pre-specified primary cell types, and in neither of "
              f"them did the interaction differ from its permutation null, in a design "
              f"whose power at the anchor effect sizes is {PW_PA} at f2 = {PW_F2A} and "
              f"{PW_PB} at f2 = {PW_F2B}. The design therefore neither supports nor "
              f"excludes a radiotherapy-associated immune shift at recurrence, and the "
              f"report is framed accordingly."),
        ("p", "The binding constraint was annotation rather than recruitment. Cluster-level "
              "argmax marker scoring assigns a whole cluster one label; five primary "
              "types formed no distinct cluster, and NK cells share markers with CD8 T "
              "cells. Re-clustering the immune compartment at three resolutions did not "
              "repair this, because the marker space did not contain the gene needed to "
              "separate immune from non-immune cells in the first place. Reporting these "
              "populations as unresolved by the annotation is more accurate than "
              "reporting them as if the cells were not there, and it changes what the "
              "cohort can be used for: the same samples could in principle support the "
              "interaction test for the myeloid compartment alone."),
        ("p", f"The GLASS arm was included because it is the largest openly available "
              f"longitudinal IDH-mutant cohort with treatment annotation, but it can only "
              f"be exploratory here: the fractions come from stock signature matrices "
              f"that were not trained on this tumour type or platform, no matched "
              f"single-nucleus reference exists from which to build a custom signature, "
              f"purity was not identifiable, and only {GL_NO} pairs lack radiotherapy "
              f"exposure. Its two q < 0.05 results are not treated as evidence about the "
              f"discovery estimates, and the two cell types that can be compared "
              f"descriptively change direction between cohorts. That difference is not "
              f"interpretable as disagreement, because both estimates are compatible with "
              f"zero and the cohorts differ in assay, deconvolution method and treatment "
              f"coding."),
        ("p", "A design able to answer the question would need three things that this "
              "cohort does not have simultaneously: a radiotherapy-only stratum large "
              "enough to separate radiotherapy from concomitant chemotherapy; a sample "
              "size chosen for a target interaction effect with the patient-level "
              "correlation accounted for; and annotation that separates the T and NK "
              "lineages, for example single-cell libraries with transcript plus surface "
              "protein or clonotype evidence, or a marker scheme with explicit NK versus "
              "CD8 T separation carried through the cluster-level decision. Cross-cohort "
              "comparison would additionally need a deconvolution signature built from "
              "the same compartment in matched tissue."),
        ("p", "Reporting this result is useful despite the null outcome, because it bounds "
              "what the current generation of openly available single-nucleus glioma "
              "data can support for treatment-stratified immune questions and it "
              "identifies annotation, not sample number, as the first constraint to fix."),
    ]

    M += [
        ("h1", "Limitations"),
        ("bullets", [
            f"The design is underpowered at the pre-specified anchor effect sizes "
            f"({PW_PA} at f2 = {PW_F2A} and {PW_PB} at f2 = {PW_F2B} at an intraclass "
            f"correlation of 0.60), so a null interaction is compatible with a modest "
            f"effect.",
            f"Only {EV_N} of {EV_TOT} pre-specified primary cell types are evaluable; "
            f"five had no proportion column and three fell below the threshold of "
            f"{EV_THR} of {EV_SAMP} samples.",
            "Annotation assigns labels to clusters rather than to cells: NK cells were "
            "absorbed into "
            "CD8-T-labelled clusters, and the pre-declared immune-subclustering "
            "sensitivity analysis did not repair this.",
            f"The GLASS arm uses stock deconvolution signatures without a matched "
            f"single-nucleus reference and without purity adjustment, and contains only "
            f"{GL_NO} pairs without radiotherapy exposure.",
            "Radiotherapy exposure is coded as any exposure and is co-exposed with "
            "chemotherapy, so a radiotherapy-specific effect cannot be isolated, and no "
            "causal language about radiotherapy is used anywhere in this report.",
            f"The per-sample cap of {D_CAP} cells ({D_NCAP} of {D_IN} samples capped) "
            f"leaves fewer cells for rare populations, although it is applied uniformly "
            f"and does not preferentially remove them.",
            "Covariates are proxies: 1p19q status enters as an IDH-subtype proxy, purity "
            "is the single-nucleus malignant-cell fraction, and MGMT methylation is "
            "unavailable.",
            "The CGGA arm is cross-sectional and cannot address the interaction; it is "
            "reported as an exploratory consistency check only.",
            "No raw or patient-level data are reported or redistributed in this package; "
            "all reported values are aggregate.",
        ]),
    ]

    B["back"] += [
        ("h1", "Data availability"),
        ("p", "Discovery single-nucleus RNA sequencing: GEO GSE326221 (Johnson/CARE "
              "cohort). The external bulk arm: GLASS, Synapse project syn17038081. The "
              "cross-sectional arm: CGGA mRNAseq_693 release (20200506). An additional "
              "controlled-access bulk arm of the discovery project (EGA "
              "EGAD50000002485) was not accessed, so bulk data are not part of this "
              "report. All values reported here are aggregate; no patient-level or raw "
              "data are included in the package, and cohort tables give counts and "
              "summaries only."),
        ("h1", "Code availability"),
        ("p", "Analysis, table and figure code is included in the submission package "
              "(submission/code) together with a frozen dependency list (Table S5). The "
              "scripts contain local absolute paths from the analysis environment and are "
              "provided for transparency rather than as a portable workflow. Every "
              "reported value is recorded in submission/supplementary/"
              "results_numbers.json with its source file and key in "
              "numbers_provenance.csv. Random seeds: the permutation null and the power "
              "simulation use seed 7331, and the subsampling seed used for the per-sample "
              "cell cap is recorded in the analysis scripts."),
        ("h1", "Acknowledgements and declarations"),
        ("p", "Acknowledgements: [ACKNOWLEDGEMENTS]."),
        ("p", "Author contributions: [AUTHOR CONTRIBUTIONS]."),
        ("p", "Funding: [FUNDING]."),
        ("p", "Competing interests: [CONFLICTS]."),
        ("p", "Ethics approval and consent to participate: [ETHICS]."),
        ("p", "Consent for publication: [CONSENT]."),
        ("p", "Use of generative artificial intelligence: [AI-USE STATEMENT]."),
        ("h1", "References"),
        ("refs", [
            "Johnson et al. Longitudinal single-nucleus and bulk profiling of IDH-mutant "
            "glioma (CARE consortium). Nature 2026. doi:10.1038/s41586-026-10612-6; "
            "PMID 42236943. [Source of the discovery cohort; citation transcribed from "
            "the project resource table.]",
            "Makarevic et al. Matched-pair histological analysis of recurrent IDH-mutant "
            "glioma after radiotherapy. Int J Mol Sci 2020. doi:10.3390/ijms21207801; "
            "PMID 33096928; PMC7590222.",
            "Varn et al. Longitudinal glioma evolution (GLASS). Cell 2022. "
            "doi:10.1016/j.cell.2022.04.038.",
            "GLASS-NL cohort (Vallentgoed et al., conference abstract 2024) - full "
            "citation to be supplied from the source publication [VERIFY CITATION].",
            "Tang et al. Longitudinal whole-transcriptome analysis of IDH-mutant glioma "
            "recurrence. Cancer Cell - year, volume and pages to be supplied "
            "[VERIFY CITATION].",
            "CGGA (Chinese Glioma Genome Atlas) mRNAseq_693 dataset, release 20200506 - "
            "primary dataset publication to be supplied [VERIFY CITATION].",
            "quanTIseq deconvolution method - method publication to be supplied "
            "[VERIFY CITATION].",
            "EPIC deconvolution method - method publication to be supplied "
            "[VERIFY CITATION].",
            "MCP-counter deconvolution method - method publication to be supplied "
            "[VERIFY CITATION].",
            "Single-cell analysis software (scanpy " + sw("scanpy") + ", anndata " +
            sw("anndata") + ", statsmodels " + sw("statsmodels") + ", matplotlib " +
            sw("matplotlib") + ", numpy " + sw("numpy") + ", pandas " + sw("pandas") +
            ") - versions are fixed in Table S5; software publications to be supplied "
            "[VERIFY CITATION].",
        ]),
    ]

    S = B["supp"]
    S += [("pb", ""), ("h1", "Supplementary figures")]
    for stem in ("FigureS1", "FigureS2", "FigureS3"):
        S += [("fig", stem)]
    S += [("h1", "Supplementary tables")]
    for tid in ("TableS1", "TableS2", "TableS3", "TableS4", "TableS5", "TableS6"):
        S += [("table", tid)]

    # inline tables: Table 1/2 in Results, 3/4/5 in Results, 6 in Methods
    insert = {
        "Table1": "Design and evaluability [PRIMARY]",
        "Table2": "Design and evaluability [PRIMARY]",
        "Table3": "Discovery interaction estimates [PRIMARY]",
        "Table4": "Exploratory external cohort (GLASS) [EXPLORATORY]",
        "Table5": "Pre-declared immune-subclustering sensitivity analysis [SENSITIVITY]",
        "Table6": "Deviations from the pre-specification (v1.0.0)",
    }
    for tid, anchor in insert.items():
        for i, (k, v) in enumerate(M):
            if k == "h2" and v == anchor:
                M.insert(i + 1, ("table", tid))
                break
        else:
            raise RuntimeError("anchor not found for " + tid)
    return B


def blocks_of(doc_blocks):
    return doc_blocks


def word_count(blocks):
    n = 0
    for kind, val in blocks:
        if kind in ("p", "h1", "h2", "h3"):
            n += len(str(val).split())
        elif kind in ("bullets", "refs"):
            for it in val:
                n += len(str(it).split())
    return n


# ------------------------------------------------------------------ renderers
def _text_of(blocks):
    out = []
    for kind, val in blocks:
        if kind in ("p", "h1", "h2", "h3", "title"):
            out.append(str(val))
        elif kind in ("bullets", "refs"):
            out.extend(str(x) for x in val)
        elif kind == "table":
            t = read_table(val)
            out.append(TAB_TITLES.get(val, val) + " " + " ".join(t["header"]) + " "
                       + " ".join(" ".join(r) for r in t["rows"]) + " " + t["notes"])
        elif kind == "fig":
            out.append(LEG[val]["legend"])
    return "\n".join(out)


def _line_numbers(doc):
    sectPr = doc.sections[0]._sectPr
    ln = OxmlElement("w:lnNumType")
    ln.set(qn("w:countBy"), "1")
    ln.set(qn("w:restart"), "newPage")
    ln.set(qn("w:distance"), "360")
    cols = sectPr.find(qn("w:cols"))
    if cols is not None:
        cols.addprevious(ln)
    else:
        sectPr.append(ln)


def _docx_table(doc, tid):
    t = read_table(tid)
    cap = "Table %s. %s" % (tid.replace("Table", ""), TAB_TITLES.get(tid, tid))
    p = doc.add_paragraph()
    r = p.add_run(cap)
    r.bold = True
    r.font.size = Pt(9)
    tab = doc.add_table(rows=1, cols=len(t["header"]))
    tab.style = "Table Grid"
    for i, h in enumerate(t["header"]):
        cell = tab.rows[0].cells[i]
        cell.text = ""
        run = cell.paragraphs[0].add_run(h)
        run.bold = True
        run.font.size = Pt(8)
    for row in t["rows"]:
        cells = tab.add_row().cells
        for i, v in enumerate(row[:len(t["header"])]):
            cells[i].text = ""
            run = cells[i].paragraphs[0].add_run(v)
            run.font.size = Pt(8)
    if t["notes"]:
        p = doc.add_paragraph()
        r = p.add_run("Note. " + t["notes"])
        r.italic = True
        r.font.size = Pt(8)


def render_docx(b, path):
    doc = Document()
    st = doc.styles["Normal"]
    st.font.name = "Calibri"
    st.font.size = Pt(11)
    st.paragraph_format.line_spacing = 1.15
    st.paragraph_format.space_after = Pt(6)
    for lvl, size in ((1, 14), (2, 12), (3, 11)):
        h = doc.styles["Heading %d" % lvl]
        h.font.name = "Calibri"
        h.font.size = Pt(size)
        h.font.bold = True
        h.font.color.rgb = RGBColor(0, 0, 0)
    _line_numbers(doc)

    seq = (b["title"] + b["abstract"] + b["main"] + b["back"] + b["supp"])
    for kind, val in seq:
        if kind == "title":
            p = doc.add_paragraph()
            r = p.add_run(val)
            r.bold = True
            r.font.size = Pt(16)
        elif kind in ("h1", "h2", "h3"):
            doc.add_heading(val, level=int(kind[1]))
        elif kind == "p":
            doc.add_paragraph(val)
        elif kind == "bullets":
            for it in val:
                doc.add_paragraph(it, style="List Bullet")
        elif kind == "refs":
            for it in val:
                p = doc.add_paragraph(it)
                p.paragraph_format.left_indent = Inches(0.25)
                p.paragraph_format.space_after = Pt(4)
                for r in p.runs:
                    r.font.size = Pt(9)
        elif kind == "pb":
            doc.add_page_break()
        elif kind == "fig":
            png = os.path.join(FIG, val + ".png")
            doc.add_picture(png, width=Inches(6.5))
            doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
            lg = LEG[val]
            p = doc.add_paragraph()
            r = p.add_run("%s. %s %s" % (val.replace("Figure", "Figure ").strip(),
                                         lg["title"], lg["legend"]))
            r.font.size = Pt(8.5)
            p.paragraph_format.space_after = Pt(10)
        elif kind == "table":
            _docx_table(doc, val)
            doc.add_paragraph("")
    doc.save(path)
    return path


def _html_text(blocks):
    out = []
    for kind, val in blocks:
        if kind == "title":
            out.append("<h1 class='title'>%s</h1>" % _esc(val))
        elif kind == "h1":
            out.append("<h2>%s</h2>" % _esc(val))
        elif kind == "h2":
            out.append("<h3>%s</h3>" % _esc(val))
        elif kind == "h3":
            out.append("<h4>%s</h4>" % _esc(val))
        elif kind == "p":
            out.append("<p>%s</p>" % _esc(val))
        elif kind == "bullets":
            out.append("<ul>" + "".join("<li>%s</li>" % _esc(x) for x in val) + "</ul>")
        elif kind == "refs":
            out.append("<ol class='refs'>"
                       + "".join("<li>%s</li>" % _esc(x) for x in val) + "</ol>")
        elif kind == "pb":
            out.append("<div class='pagebreak'></div>")
        elif kind == "fig":
            lg = LEG[val]
            out.append(
                "<figure><img src='../figures/%s.png' alt=\"%s\">"
                "<figcaption><strong>%s. %s</strong> %s</figcaption></figure>"
                % (val, _esc(lg["alt"]), val.replace("Figure", "Figure ").strip(),
                   _esc(lg["title"]), _esc(lg["legend"])))
        elif kind == "table":
            t = read_table(val)
            n = val.replace("Table", "")
            head = "".join("<th>%s</th>" % _esc(h) for h in t["header"])
            body = "".join("<tr>" + "".join("<td>%s</td>" % _esc(c) for c in r) + "</tr>"
                           for r in t["rows"])
            note = ("<p class='tnote'>Note. %s</p>" % _esc(t["notes"])) if t["notes"] else ""
            out.append("<div class='tbl'><p class='tcap'><strong>Table %s. %s</strong></p>"
                       "<table><thead><tr>%s</tr></thead><tbody>%s</tbody></table>%s</div>"
                       % (n, _esc(TAB_TITLES.get(val, val)), head, body, note))
    return "\n".join(out)


def _esc(s):
    return (str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


CSS = """
body{font-family:Calibri,Carlito,"Segoe UI",Arial,sans-serif;font-size:11pt;
 line-height:1.15;color:#111;max-width:820px;margin:32px auto;padding:0 18px}
h1.title{font-size:17pt;line-height:1.2}
h2{font-size:14pt;margin-top:22px}
h3{font-size:12pt;margin-top:16px}
h4{font-size:11pt;margin-top:12px}
figure{margin:18px 0;page-break-inside:avoid}
img{width:100%;max-width:100%;height:auto;border:1px solid #ddd}
figcaption{font-size:9pt;color:#222;margin-top:6px}
table{border-collapse:collapse;font-size:8.5pt;width:100%}
th,td{border:1px solid #999;padding:3px 4px;text-align:left;vertical-align:top}
th{background:#f0f0f0}
.tcap{font-size:9pt;margin-bottom:2px}
.tnote{font-size:8.5pt;font-style:italic;color:#333}
.refs{font-size:9.5pt}
.pagebreak{page-break-after:always}
@media print{body{max-width:none}}
"""


def render_html(b, path):
    seq = (b["title"] + b["abstract"] + b["main"] + b["back"] + b["supp"])
    html = ("<!DOCTYPE html>\n<html lang='en'>\n<head>\n<meta charset='utf-8'>\n"
            "<meta name='viewport' content='width=device-width, initial-scale=1'>\n"
            "<title>%s</title>\n<style>%s</style>\n</head>\n<body>\n%s\n</body>\n</html>\n"
            % (_esc(TITLE), CSS, _html_text(seq)))
    with io.open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(html)
    return path


# ------------------------------------------------------------------ checks
def check_abstract():
    sents = [s for s in re.split(r"(?<=[.!?])\s+", ABSTRACT.strip()) if s]
    print("abstract sentences: %d" % len(sents))
    for i, s in enumerate(sents, 1):
        print("  S%d words=%d" % (i, len(s.split())))
    if len(sents) != 5:
        raise SystemExit("STOP: the abstract must be exactly 5 sentences")
    return len(sents)


def check_banned(text):
    hits = []
    for pat in BANNED:
        for m in re.finditer(pat, text, flags=re.IGNORECASE):
            hits.append((pat, text[max(0, m.start() - 45):m.start() + 45].replace("\n", " ")))
    print("banned-term hits: %d" % len(hits))
    for w, ctx in hits:
        print("  [%s] ...%s..." % (w, ctx))
    return hits


def main():
    b = build_blocks()

    wc_abs = len(ABSTRACT.split())
    wc_main = word_count(b["main"])
    metrics = ("Manuscript metrics. Abstract: %d words in 5 sentences. Main text "
               "(Introduction through Limitations, excluding the abstract, keywords, "
               "tables, figure captions, references and supplementary material): %d "
               "words. Tables 1-6 and Figures 1-4 are embedded; Tables S1-S6 and "
               "Figures S1-S3 are supplementary." % (wc_abs, wc_main))
    b["back"].insert(b["back"].index(("h1", "References")), ("p", metrics))

    docx_path = render_docx(b, os.path.join(MAN, "manuscript.docx"))
    html_path = render_html(b, os.path.join(MAN, "manuscript.html"))
    print("wrote", docx_path)
    print("wrote", html_path)
    print(metrics)

    n_sent = check_abstract()
    full = _text_of(b["title"] + b["abstract"] + b["main"] + b["back"])
    check_banned(full)
    with io.open(os.path.join(MAN, "manuscript_metrics.json"), "w",
                 encoding="utf-8") as f:
        json.dump({"title_options": TITLE_OPTIONS, "title_used": TITLE,
                   "abstract_sentences": n_sent, "abstract_words": wc_abs,
                   "main_text_words": wc_main}, f, indent=2)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
