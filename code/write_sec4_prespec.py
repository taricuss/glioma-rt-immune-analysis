
"""
Section 4.5 — Pre-specification (LOCKED before any analysis runs)
=================================================================
Plan §4.5 and §5.5 are explicit: exact cell-type list, model formula,
multiple-testing correction method, permutation scheme, and — where
applicable — mechanistic marker sets must be written & date-stamped
BEFORE running deconvolution.  This script writes two artifacts:

    1. outputs/sec4_prespec.json   (machine-readable, versioned, sha256 of itself)
    2. outputs/sec4_prespec.md     (human-readable, same content + rationale)

Rule for any deviation later:
  If a cell type, covariate, model, or correction is changed AFTER partial
  results are seen, every downstream figure that uses it is labelled
  "EXPLORATORY" in the manuscript draft.  Primary claims come only from
  this locked list.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
OUT = HERE / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

VERSION = "1.0.0"
LOCKED_AT = datetime.now(timezone.utc).isoformat(timespec="seconds")

# =========================================================================
# (a) EXACT CELL-TYPE LIST TESTED (primary)
#
# Rationale for this exact list:
#   Johnson et al. (Nature 2026) Fig 5 stratified myeloid populations on RT
#   status and reported significance for the myeloid panel.  The novelty
#   gap is that they did NOT repeat the same RT-stratified analysis for
#   the lymphocyte panel.  We therefore test every immune population
#   that can be resolved both by CIBERSORTx custom-reference (built from
#   the matched snRNA-seq in §5.4) and by quanTIseq/EPIC/MCP-counter for
#   the §5.4 robustness-check.
#
#   Ordering = lymphocyte-focused first (the gap), then myeloid comparators
#   (to reproduce Johnson's result as a positive-control sanity check,
#   labelled as such in any write-up — not primary claims).
# =========================================================================
CELL_TYPES_PRIMARY = [
    # ---- Lymphocyte panel (the gap) ----
    "CD8_T_cells",
    "CD4_T_cells_conventional",
    "T_regulatory_cells_FOXP3",
    "NK_cells",
    "B_cells",
    "Plasma_cells",
    # ---- Myeloid comparators (positive control, Johnson Fig 5 already
    #      showed signal here — report as replication, not primary claim) ----
    "Microglia_P2RY12_TMEM119",
    "Monocyte_derived_macrophages_CD163",
    "Dendritic_cells",
    "Neutrophils",
]
CELL_TYPES_EXPLORATORY = [
    "T_cells_exhausted_PDCD1_HAVCR2",   # requires Step-5 mechanistic flag trigger
    "M1_macrophages",                   # only if Step 3 signal hits myeloid
    "M2_macrophages",
    "Endothelial_cells",
    "Fibroblasts_pericytes",
]

# =========================================================================
# (b) EXACT MODEL FORMULA (per cohort)
#
# §5.5 rule: `fraction ~ timepoint * treatment_arm + purity + grade + age
#             + (1 | patient_id)` with fraction arcsine- or logit-
#             transformed.  We write it twice — discovery + validation —
#             because the treatment variable has different granularity
#             (plan §4.3: discovery = binary Any_RT / No_RT; GLASS =
#             four-level if n supports reporting per-cell).
# =========================================================================
MODEL_DISCOVERY = {
    "description": (
        "Per-cell-type linear mixed model.  Response Y is arcsin(sqrt(fraction))"
        " — arcsin-square-root is preferred over logit for fractions because"
        " the deconvolution output has many exact zeros in low-infiltration"
        " glioma samples; logit(0) = -Inf forces a pseudo-count that can"
        " distort the interaction coefficient."
    ),
    "formula_r_style": (
        "arcsin_sqrt_fraction ~ timepoint * treatment_arm"
        "   + estimate_purity"
        "   + who_grade_initial"
        "   + age_at_diagnosis"
        "   + codeletion_1p19q"
        "   + mgmt_methylation"
        "   + (1 | patient_id)"
    ),
    "response_transform": "arcsin_sqrt",  # alternative: "logit" — recorded for robustness check
    "hypothesis_term": "timepoint:treatment_arm",
    "treatment_definition": "BINARY — Any_RT (n=20) vs No_RT (n=13), exactly per §4.3 lock below.",
    "timepoint_levels": ["initial", "recurrent"],
    "treatment_levels": ["No_RT", "Any_RT"],
    "covariates_mandatory": [
        "estimate_purity",         # plan §5.3: non-negotiable
        "who_grade_initial",
        "age_at_diagnosis",
        "codeletion_1p19q",
        "mgmt_methylation",
    ],
    "random_effect_structure": "(1 | patient_id)",
    "notes": (
        "Do NOT add chemo as a within-arm covariate: RT-only n=3 cannot support"
        " that estimate (plan §4.3 hard rule)."
    ),
}

MODEL_GLASS = {
    "description": "Mirror of discovery formula; treatment can be four-level (untreated / TMZ_only / RT_only / RT_TMZ) if per-arm n > 3 after IDH filter.",
    "formula_r_style": MODEL_DISCOVERY["formula_r_style"],
    "response_transform": MODEL_DISCOVERY["response_transform"],
    "hypothesis_term": MODEL_DISCOVERY["hypothesis_term"],
    "treatment_definition_fallback": "If four-level per-cell n collapses below 3 in any arm, fall back to binary Any_RT/No_RT discovery definition.",
    "covariates_mandatory": MODEL_DISCOVERY["covariates_mandatory"],
}

# =========================================================================
# (c) EXACT MULTIPLE-TESTING CORRECTION
#
# BH-FDR across K = len(CELL_TYPES_PRIMARY) = 10 tests in the primary
# family.  Exploratory cell types are NOT included in the same correction
# pass (plan §4.5 pre-spec rule).  Any p-values from the exploratory list
# are reported as raw p with an explicit "EXPLORATORY — uncorrected"
# footnote in every table that carries them.
#
# Why not Bonferroni: 10 tests → Bonferroni threshold 0.005 is overly
#   conservative for a confirmatory study targeting IF ≥ 7 journals,
#   every one of which explicitly prefers BH-FDR over Bonferroni as the
#   default across-family correction for deconvolution studies.
# =========================================================================
MULTIPLE_TESTING = {
    "family": "CELL_TYPES_PRIMARY",
    "k_tests": len(CELL_TYPES_PRIMARY),
    "method": "Benjamini-Hochberg step-up (BH-FDR)",
    "alpha_family": 0.05,
    "rejection_threshold_convention": "q-value < 0.05 per cell type within the primary family.",
    "exploratory_policy": (
        "CELL_TYPES_EXPLORATORY are a separate family.  Never merge them into"
        " the primary BH-FDR pass; report raw p + 'EXPLORATORY' label."
    ),
    "robustness_checks": [
        "Bonferroni-adjusted p (0.05 / 10 = 0.005) reported alongside BH q as a sensitivity row.",
    ],
}

# =========================================================================
# (d) EXACT PERMUTATION SCHEME (plan §5.5 + §4.5)
#
# Goal: empirical null for the interaction coefficient, because at n=33
# patients the asymptotic F-test's tail accuracy is questionable.
# =========================================================================
PERMUTATION = {
    "units_shuffled": "patient-level treatment labels (timepoints stay paired within patient).",
    "n_permutations": 10000,
    "seed": 7331,
    "how": (
        "For each of 10000 iterations, sample an independent permutation of"
        " the length-N_patients treatment vector WITHOUT stratification (no"
        " block preservation by subtype — conservative; subtype IS a covariate"
        " in the model so residual subtype-balance-drift is already modelled)."
        " Fit the exact same LMM formula each time, extract the interaction"
        " coefficient, build the two-tailed empirical CDF, compute empirical"
        " p = (1 + count(|perm_coef| >= |observed_coef|)) / (1 + n_perm)."
    ),
    "reporting": (
        "For every primary cell type, four p-values appear side-by-side:"
        " (1) asymptotic LMM p, (2) LMM Satterthwaite-KE df adjusted p,"
        " (3) permutation-based empirical p, (4) BH-FDR q derived from (3)."
        " Primary claim uses (4); (1) and (2) shown for reviewers who ask."
    ),
    "parallelism": "Use joblib loky backend, n_jobs = min(cpu_count, 8).",
}

# =========================================================================
# (e) MECHANISTIC MARKER SETS — pre-specified, fishing-expedition ban
#
# §5.7 hard rule: only check marker sets IF the primary test (Step 3
# and/or Step 4) hits a specific population with q < 0.05.  If it does,
# these exact lists, nothing else.
# =========================================================================
MECHANISTIC_MARKERS = {
    "T_cell_exhaustion": {
        "trigger_cell_type": "CD8_T_cells",
        "genes": ["PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "TOX", "LAYN"],
        "test": "At RECURRENT timepoint, Any_RT vs No_RT DESeq2/limma contrast restricted to CD8+ deconvolution-pseudo-bulk (or snRNA CD8+ cluster if available).",
    },
    "T_regulatory_suppression": {
        "trigger_cell_type": "T_regulatory_cells_FOXP3",
        "genes": ["FOXP3", "IL2RA", "CTLA4", "TIGIT", "IL10", "TGFB1", "IKZF2", "IKZF4"],
        "test": "Same recurrent-only Any_RT vs No_RT design.",
    },
    "Myeloid_polarization_M1_M2": {
        "trigger_cell_type_any": ["Microglia_P2RY12_TMEM119", "Monocyte_derived_macrophages_CD163"],
        "M1_genes": ["CD80", "CD86", "IL1B", "IL6", "TNF", "NOS2", "IRF5"],
        "M2_genes": ["CD163", "MRC1", "MSR1", "ARG1", "IL10", "TGFB1", "CCL22", "CHI3L1"],
        "test": "Module score (AddModuleScore in Seurat) per sample, Any_RT vs No_RT at recurrent timepoint, Mann–Whitney U.",
    },
    "GSEA_restriction": {
        "allowed_databases": [
            "MSigDB HALLMARK",
            "ImmuneSigDB C7 (immunologic signatures) — restricted to sub-collections with keywords: glioma, lymphocyte, tcell, exhaustion, myeloid, radiotherapy.",
        ],
        "forbidden": (
            "Never run the full MSigDB C2/C5/C7 unrestricted.  Never use KEGG"
            " as a primary pathway layer (too coarse)."
        ),
    },
}

# =========================================================================
# (f) DISCOVERY-COHORT TREATMENT-LABEL LOCK (plan §4.3 + Phase-0 report)
#
# Two patients (of 35 total) have UNCLEAR / NA treatment labels.  They
# are EXCLUDED from all primary models, listed by ID below.  Any later
# sensitivity-analysis re-inclusion must carry an explicit
# "post-hoc sensitivity" label and NEVER appear in the primary forest
# plot or primary abstract numbers.
# =========================================================================
TREATMENT_LOCK = {
    "total_patients_in_johnson_cohort": 35,
    "patients_excluded_treatment_na": 2,
    "excluded_note": (
        "Patient IDs in this block are to be filled from Johnson et al."
        " Nature 2026 Supplementary Table 1 column 'treatment_between_timepoints'"
        " = 'Not Reported' / NA.  The acquire_clinical_metadata.py script"
        " reads Supp Table 1 directly and cross-checks these counts."
    ),
    "n_any_rt": 20,
    "any_rt_breakdown": {
        "rt_only": 3,
        "rt_plus_chemo_confirmed": 16,
        "rt_plus_chemo_unconfirmed_chemo": 1,
    },
    "n_no_rt": 13,
    "no_rt_breakdown": {
        "chemo_only": 6,
        "untreated_no_adjuvant": 7,
    },
    "rule_never_break": (
        "RT-only subgroup (n=3) must NEVER be reported as a standalone"
        " inferential test in the primary section.  Descriptive proportions"
        " only; any RT-only contrast goes to Appendix labelled"
        " 'n-too-small — descriptive only, no p-value reported'."
    ),
    "same_rule_applies_to_glass_val_if_n_too_small": True,
}

# =========================================================================
# (g) EXTERNAL VALIDATION — pre-specified forest-plot spec (plan §5.8)
# =========================================================================
FOREST_PLOT_SPEC = {
    "x_axis": "Interaction coefficient (timepoint × treatment), arcsin-sqrt-fraction units.",
    "rows_order": "Same order as CELL_TYPES_PRIMARY.",
    "two_columns": "Discovery (left) | GLASS validation (right).",
    "error_bars": "95% profile-likelihood CI for the interaction term; if profile-likelihood fails → 95% Wald CI with explicit footnote.",
    "diamond": "Fixed-effects meta-analytic pooled estimate across the two cohorts ONLY if I² heterogeneity p > 0.10 (i.e. homogeneous).  If heterogeneous, draw separate circles + no pooled diamond.",
    "annotation": "Every row carries BH-FDR q-value from the discovery-cohort permutation test.",
}

PRESPEC = {
    "version": VERSION,
    "locked_at": LOCKED_AT,
    "a_cell_types_primary": CELL_TYPES_PRIMARY,
    "a_cell_types_exploratory": CELL_TYPES_EXPLORATORY,
    "b_model_discovery": MODEL_DISCOVERY,
    "b_model_glass": MODEL_GLASS,
    "c_multiple_testing": MULTIPLE_TESTING,
    "d_permutation": PERMUTATION,
    "e_mechanistic_markers": MECHANISTIC_MARKERS,
    "f_treatment_lock": TREATMENT_LOCK,
    "g_forest_plot": FOREST_PLOT_SPEC,
    "integrity_placeholder": None,  # filled below after first hash
}

# Compute integrity hash over the whole doc (minus its own placeholder field)
def self_hash(d: dict) -> str:
    d2 = {k: v for k, v in d.items() if k != "integrity_placeholder"}
    blob = json.dumps(d2, sort_keys=True, ensure_ascii=False, default=str).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()

PRESPEC["integrity_placeholder"] = self_hash(PRESPEC)

# ------------------
# Write JSON artifact
# ------------------
json_path = OUT / "sec4_prespec.json"
with open(json_path, "w", encoding="utf-8") as fh:
    json.dump(PRESPEC, fh, indent=2, ensure_ascii=False)
print(f"WROTE {json_path}  (sha256 = {PRESPEC['integrity_placeholder'][:16]}…)")

# ------------------
# Write Markdown artifact (same content, human-readable)
# ------------------
md_lines = []
md_lines.append(f"# Section 4.5 — Pre-specification  (LOCKED v{VERSION})")
md_lines.append("")
md_lines.append(f"**Locked at (UTC):** {LOCKED_AT}  ")
md_lines.append(f"**Integrity sha256:** `{PRESPEC['integrity_placeholder']}`  ")
md_lines.append("")
md_lines.append("> **Rule for any deviation later:** If any field in this file is altered ")
md_lines.append("> *after* partial results from deconvolution/statistics have been ")
md_lines.append("> inspected, every downstream figure or table that uses the changed ")
md_lines.append("> setting must be labelled **EXPLORATORY** in the manuscript.  Primary ")
md_lines.append("> claims come only from this locked v1.0.0 document.")
md_lines.append("")
md_lines.append("---")
md_lines.append("")

def _kv_list(title: str, items: list[str]):
    md_lines.append(f"## {title}")
    md_lines.append("")
    for i, it in enumerate(items, 1):
        md_lines.append(f"{i}. `{it}`")
    md_lines.append("")

_kv_list("(a) Primary cell types tested (BH-FDR family K=10)", CELL_TYPES_PRIMARY)
_kv_list("(a) Exploratory cell types (raw p only, NEVER merged into primary BH-FDR)",
         CELL_TYPES_EXPLORATORY)

md_lines.append("## (b) Exact model formulas")
md_lines.append("")
md_lines.append("### Discovery cohort (Johnson CARE / Nature 2026)")
md_lines.append("")
md_lines.append(f"```r")
md_lines.append(MODEL_DISCOVERY["formula_r_style"])
md_lines.append("```")
md_lines.append("")
md_lines.append(f"- Response transform: **{MODEL_DISCOVERY['response_transform']}** (arcsin-square-root preferred over logit because many deconvolution fractions are exactly 0 in glioma).")
md_lines.append(f"- Hypothesis term: `{MODEL_DISCOVERY['hypothesis_term']}`")
md_lines.append(f"- Treatment: **{MODEL_DISCOVERY['treatment_definition']}** (binary levels: {MODEL_DISCOVERY['treatment_levels']}; timepoint levels: {MODEL_DISCOVERY['timepoint_levels']}).")
md_lines.append(f"- Mandatory covariates (modeled as fixed effects every run): **{', '.join(MODEL_DISCOVERY['covariates_mandatory'])}**.")
md_lines.append(f"- Random effects: `{MODEL_DISCOVERY['random_effect_structure']}`.")
md_lines.append(f"- Note: {MODEL_DISCOVERY['notes']}")
md_lines.append("")
md_lines.append("### GLASS validation cohort")
md_lines.append("")
md_lines.append(f"- Mirror formula; four-level treatment (untreated / TMZ_only / RT_only / RT_TMZ) used only if each non-empty cell has n ≥ 3.  Otherwise fall back to binary discovery definition (**{MODEL_GLASS['treatment_definition_fallback']}**).")
md_lines.append("")

md_lines.append("## (c) Multiple-testing correction")
md_lines.append("")
md_lines.append(f"- Family: CELL_TYPES_PRIMARY, K = **{MULTIPLE_TESTING['k_tests']}** tests.")
md_lines.append(f"- Method: **{MULTIPLE_TESTING['method']}**.")
md_lines.append(f"- α (family-wide): 0.05.  Rejection convention: **{MULTIPLE_TESTING['rejection_threshold_convention']}**.")
md_lines.append(f"- Exploratory policy: {MULTIPLE_TESTING['exploratory_policy']}")
md_lines.append(f"- Robustness: {MULTIPLE_TESTING['robustness_checks'][0]}")
md_lines.append("")

md_lines.append("## (d) Permutation-based empirical null for the interaction coefficient")
md_lines.append("")
md_lines.append(f"- Units shuffled: {PERMUTATION['units_shuffled']}")
md_lines.append(f"- Iterations: **{PERMUTATION['n_permutations']}**, seed = `{PERMUTATION['seed']}`.")
md_lines.append(f"- How: {PERMUTATION['how']}")
md_lines.append(f"- Reporting: {PERMUTATION['reporting']}")
md_lines.append(f"- Parallelism: {PERMUTATION['parallelism']}")
md_lines.append("")

md_lines.append("## (e) Mechanistic marker sets — CONDITIONAL (only if primary hit triggers them)")
md_lines.append("")
for k, m in MECHANISTIC_MARKERS.items():
    if k == "GSEA_restriction":
        md_lines.append(f"### GSEA restriction (§5.7 hard rule)")
        md_lines.append(f"- Allowed: {', '.join(m['allowed_databases'])}.")
        md_lines.append(f"- Forbidden: {m['forbidden']}")
        continue
    trig = m.get("trigger_cell_type", m.get("trigger_cell_type_any", "—"))
    md_lines.append(f"### {k}")
    md_lines.append(f"- Trigger: {trig} (q < 0.05 in primary test).")
    if "genes" in m:
        md_lines.append(f"- Gene list: `{', '.join(m['genes'])}`.")
    if "M1_genes" in m:
        md_lines.append(f"- M1 list: `{', '.join(m['M1_genes'])}`.")
        md_lines.append(f"- M2 list: `{', '.join(m['M2_genes'])}`.")
    if "test" in m:
        md_lines.append(f"- Test design: {m['test']}")
    md_lines.append("")

md_lines.append("## (f) Discovery-cohort treatment-label lock (plan §4.3)")
md_lines.append("")
md_lines.append(f"- Total cohort: 35 patients.")
md_lines.append(f"- Excluded (unclear / NA treatment): **2 patients** (IDs populated from Johnson et al. Supp Table 1 by `scripts/acquire_clinical_metadata.py`; listed in every output that carries the primary fit).")
md_lines.append(f"- **Any-RT arm: n = {TREATMENT_LOCK['n_any_rt']}** patients (RT-only={TREATMENT_LOCK['any_rt_breakdown']['rt_only']}, RT+chemo confirmed={TREATMENT_LOCK['any_rt_breakdown']['rt_plus_chemo_confirmed']}, RT+chemo chemo-unconfirmed={TREATMENT_LOCK['any_rt_breakdown']['rt_plus_chemo_unconfirmed_chemo']}).")
md_lines.append(f"- **No-RT arm: n = {TREATMENT_LOCK['n_no_rt']}** patients (chemo-only={TREATMENT_LOCK['no_rt_breakdown']['chemo_only']}, untreated={TREATMENT_LOCK['no_rt_breakdown']['untreated_no_adjuvant']}).")
md_lines.append(f"- Rule never break: {TREATMENT_LOCK['rule_never_break']}")
md_lines.append("")

md_lines.append("## (g) Forest plot — pre-specified visual")
md_lines.append("")
md_lines.append(f"- X-axis: {FOREST_PLOT_SPEC['x_axis']}")
md_lines.append(f"- Rows order: CELL_TYPES_PRIMARY order.")
md_lines.append(f"- Two panels: Discovery | GLASS.  Error bars: 95% profile-likelihood CI; Wald fallback with footnote if profile fails.")
md_lines.append(f"- Pooled diamond only if I² p > 0.10 (no significant cross-cohort heterogeneity).")
md_lines.append(f"- Annotations: BH-FDR q (from permutation-based p) on every discovery row.")
md_lines.append("")
md_lines.append("---")
md_lines.append(f"*End of pre-spec v{VERSION}.  Last edit MUST pre-date first `scripts/acquire_bulk_discovery.py` run; the latter script aborts if this file does not exist.*")

md_path = OUT / "sec4_prespec.md"
with open(md_path, "w", encoding="utf-8") as fh:
    fh.write("\n".join(md_lines))
print(f"WROTE {md_path}  ({len(md_lines)} lines)")
