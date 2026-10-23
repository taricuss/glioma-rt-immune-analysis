# Stage 2 plan — second-stage re-annotation and re-analysis of the RT-interaction question in IDH-mutant glioma

**Status:** locked before execution. Lock record (SHA-256 of this file + UTC timestamp): `stage2_plan.lock.json` in this directory.

**Plan locked (UTC):** 2026-10-21T08:57:05Z
**Companion documents:** `revisionplan.md` (operative revision spec), `jno_requirements.md` (journal requirements), `audit_report.md` (Phase 1 reproduction + treatment-blind diagnosis), `outputs/sec4_prespec.md` / `.json` (**Stage 1 pre-specification, locked 2026-09-13, never modified**), `submission/supplementary/results_numbers.json` + `numbers_provenance.csv` (numbers registry).

---

## 1. Why this stage exists, and what it is allowed to claim

The Stage 1 pre-specified analysis found that only 2 of 10 pre-specified cell types were evaluable (NK and Microglia; both non-significant). The Phase 1 audit showed that this failure was **caused by our own annotation step**, not by the biology: Stage 1 scored marker panels in a 2,000-gene highly variable space in which `PTPRC` and 16 other panel genes were absent, so no immune gate was possible and several panels lost enough markers that `annotate_clusters()` skipped them entirely (`< 3` genes present → no score column → the type could never win a cluster).

**This is a second-stage (post-hoc) analysis, not a pre-specified one.** It is designed after Stage 1 results were known. It earns credibility only from three things, all of which are fixed in this document before any Stage 2 interaction estimate is computed or viewed:

1. the annotation is **treatment-blind and timepoint-blind** by construction (§3, §4);
2. every QC gate and the pass/fail decision rule are **pre-declared** (§4, §5);
3. no parameter is adjusted after seeing a Stage 2 interaction estimate (§13).

Nothing in Stage 1 is changed. Stage 1 results are reported in the manuscript exactly as locked.

### 1.1 Decisions recorded at STOP 2 (user-confirmed 2026-10-21)

| Decision | Resolution |
|---|---|
| Stage 1 power figure is not reproducible (61.9% → 63.25%) | Report it as an **approximate Monte-Carlo estimate with the simulation size disclosed** ("≈62%, nsim = 800"). No Stage 1 number is re-locked; the locked pre-spec and its hash are untouched. See §7.4. |
| B and plasma cells | **Pre-declared as one pooled primary type `B_plasma`** (§5.2). Pooled positivity reaches the ≥10-sample rule in 14 of 75 samples; separately, neither would. |
| Stage 2 evaluability gate | **≥6 of the 9 Stage 2 primary types evaluable, and `CD8_T`, `NK` and `Microglia_P2RY12_TMEM119` each individually evaluable** (§5.3). |
| Stage 1 permutation null re-run | **Not re-run.** Both of its inputs reproduce bit-exactly (`audit_report.md` §1.1). Locked Stage 1 values remain the reported ones. |

---

## 2. Inputs, cell set and gene space

| Item | Value / rule |
|---|---|
| Expression matrix | `data/discovery/checkpoints/ckpt_01_merge_raw.h5ad` — **221,016 cells × 12,314 genes** |
| Cell set | **Identical to Stage 1.** Every checkpoint (`ckpt_01`…`ckpt_05`) holds the same 221,016 cells in the same order, so the Stage 2 analysis reuses the Stage 1 cells, the Stage 1 per-sample cap (3,000 cells/sample, seed 20260917) and the Stage 1 sample membership, changing **only the cell-type labels** |
| Gene space change | 2,000 → 12,314 genes. This is the single deliberate change and the reason the stage exists |
| Barcode check | Before scoring, the `obs/_index` of `ckpt_01` and `ckpt_05` must be compared element-wise. If they differ, **stop and report**; do not proceed |
| Fraction-table denominator | 75 snRNA samples (`sn_proportions_per_sample.csv` row count, verified) |
| Clinical table | `outputs/intermediate/clinical_canonical.csv` (71 samples / 33 patients / 8 columns) |
| Purity covariate | `outputs/intermediate/estimate_scores.csv` → `estimate_purity` |

No raw GEO matrix is re-merged. This was declined for resource reasons (14 GiB, `data/discovery/GSE326221/GSE326221_RAW_extracted/`) and is unnecessary: all 17 genes missing from the HVG space are present in `ckpt_01`.

---

## 3. Treatment-blind hierarchical annotation

### 3.1 Scoring function (fixed in advance)

For every cell and every panel `P`:

```
score_P = mean_{g in P, g present in ckpt_01} log1p( CP10K(count_g) )
          − mean_{c in C} log1p( CP10K(count_c) )
```

* `CP10K(count_g) = count_g / total_counts(cell) * 1e4`, computed from the raw counts in `ckpt_01`.
* **Control set `C`** (a fixed, treatment-blind, data-derived set of 50 genes): all genes present in `ckpt_01` that appear in **no** panel below, ranked by total detection count (number of cells with count > 0) across all 221,016 cells, descending, ties broken alphabetically; the top 50 are used. `C` is computed once, from the whole cohort, before any behaviour of any panel is inspected, and is shared across all panels so that it cancels in every comparison.
* No clustering step is used. Annotation is **per cell**, which removes the cluster-mean `argmax` mechanism that produced the Stage 1 failure.
* `P` genes absent from `ckpt_01` are recorded in `stage2_annotation_qc.md`, and a panel with fewer than 3 present genes fails gate **G1** (§4).

### 3.2 Panel definitions

Every Stage 1 `MARKERS` gene is **retained** in the corresponding Stage 2 panel — no marker is dropped, only restored. Restored genes (absent from the 2,000-HVG space) are marked `*`.

**Layer 1 — immune vs non-immune**

* `IMMUNE_PAN_LEUKO`: PTPRC `*`, CD53, LAPTM5, CORO1A, ARHGDIB, LSP1, DOCK2, IKZF1, SPN, CD48
* Non-immune panels (Stage 1 panels, unchanged): `Malignant_astrocytes_oligo` (SOX2, VIM, S100B, GFAP, OLIG2, IDH1), `Astrocytes_non_malignant` (GFAP, AQP4, ALDH1L1, S100B, GJA1, CLDN10), `Oligodendrocytes_OPC` (OLIG1, OLIG2, MBP, MAG, MOG, PLP1, PDGFRA, CSPG4), `Endothelial_cells` (PECAM1, VWF, FLT1, KDR, CD34, ESM1), `Pericytes_Fibroblasts` (RGS5, PDGFRB, CSPG4, COL1A1, ACTA2, DES)

Rule: the cell is **immune** if `score_IMMUNE_PAN_LEUKO` exceeds every non-immune panel score; otherwise it takes the winning non-immune panel label. (Non-immune labels do not enter the primary estimand, but they are retained for QC and for the composition table.)

**Layer 2 — within immune cells**

* `TNK`: CD3D, CD3E, CD3G, TRAC, TRBC2, CD2, CD247, LCK, IL7R
* `MYELOID`: P2RY12, TMEM119, CX3CR1, CSF1R, IRF8, AIF1, CD163, MRC1, MSR1, CD14, CCR2, S100A9, LYZ, ITGAM `*`, FCGR3A, ITGAX `*`, TYROBP, C1QA, C1QB, SPI1
* `BPLASMA`: CD19, CD79A, CD79B, MS4A1, IGHM, IGHD, MZB1, XBP1, PRDM1, TNFRSF17, IGHG1, IGKC, BANK1 `*`, JCHAIN, SDC1 `*`, CD37

**Layer 3 — within `TNK`**

* `CD8_T`: CD3D, CD3E, CD8A, CD8B, GZMA, GZMB, PRF1, NKG7, GZMK
* `CD4_conv`: CD3D, CD3E, CD4, IL7R, CCR7, LEF1, TRAT1, SELL, LTB
* `Treg_FOXP3`: CD4, FOXP3, IL2RA, CTLA4, IKZF2, IKZF4, TIGIT, TNFRSF4, TNFRSF18
* `NK`: NCAM1 `*`, KLRD1, GNLY, NKG7, FCGR3A, KLRF1, NCR1 `*`, KLRB1, SH2D1B, TYROBP, PRF1, GZMA, GZMB

**Hard CD3 rule (fixed in advance):** a cell may only be labelled `NK` if `count(CD3D) + count(CD3E) + count(CD3G) == 0`. If the `NK` panel wins but CD3 is detected, the cell is reassigned to the winning CD3-defining panel (`CD8_T` / `CD4_conv` / `Treg_FOXP3`). The number of cells moved by this rule is reported per sample, treatment-blind.

**Layer 4 — within `MYELOID`**

* `Microglia_P2RY12_TMEM119`: P2RY12, TMEM119, CX3CR1, CSF1R, IRF8, AIF1
* `MoMF_CD163`: CD163, MRC1, MSR1, CD14, CCR2, S100A9, LYZ, ITGAM `*`
* `DC`: FCER1A, CLEC10A, CD1C, CD141, IRF8, HLA-DRA, ITGAX `*`, CLEC9A `*`, LILRA4 `*`
* `Neutrophils`: S100A8, S100A9, CSF3R, CXCR1, CEACAM8, MPO, FCGR3B `*`

### 3.3 Ambiguity (doublet) rule

Within each layer, let `best` and `second` be the top two panel scores for the cell. A cell is flagged `ambiguous_<layer>` if `best − second < δ` with **δ = 0.05** (score units; pre-declared, not tuned).

* Ambiguous cells are **excluded from the numerator of every type** but remain in the denominator ("fraction of all cells"), so the denominators are unaffected by the rule.
* The ambiguous fraction is reported per layer and per type.
* Pre-declared sensitivity: the whole annotation is repeated at **δ = 0.025** and **δ = 0.10**. These are sensitivity checks on a fixed parameter, not tuning.

### 3.4 Independent second method (internal check)

No external reference atlas is available offline for this dataset under the current resource limits; that absence is stated as a limitation in the manuscript and is not worked around by importing a reference.

The independent check is an internal method with a **different scoring principle**:

* Method B: per cell, the **binary detection vector** over the same panels (`count > 0`), scored as the fraction of panel genes detected, control-corrected the same way; nearest-centroid assignment.
* Method A = §3.1. Agreement is quantified with Cohen's κ, overall and per layer (gate **G3**).

---

## 4. Annotation QC gates (pre-declared, treatment-blind)

All gates are evaluated on columns that cannot depend on `timepoint` or `treatment_arm`. None of the diagnostics reads either column. Every gate outcome is reported in `stage2_annotation_qc.md` whether it passes or fails.

| Gate | Rule | Fail action |
|---|---|---|
| **G1** panel completeness | each panel entering the primary set retains ≥3 genes present in `ckpt_01` | the affected type is **declared non-evaluable with a recorded reason**; it is never silently dropped (Stage 1's silent skip is the defect being corrected) |
| **G2** minimum cells & spread | ≥100 cells in the cohort **and** non-zero in ≥10 of 75 samples (unchanged pre-specified evaluability rule) | type is non-evaluable |
| **G3** inter-method agreement | κ(A, B) ≥ 0.60 for the Layer-1 immune/non-immune decision, ≥ 0.50 for the Layer-2 three-way split, and ≥ 0.40 for each primary type individually | gate false → annotated type set rejected; stop |
| **G4** marker specificity | (i) for ≥70% of a type's cells, that type's own panel is the argmax within its layer; (ii) the type's own-panel mean score in its own cells exceeds that panel's mean score in every other type's cells (one-vs-rest separation) | gate false for a type → that type is non-evaluable |
| **G5** ambiguity | ambiguous fraction ≤ 20% per type (§3.3) | type is non-evaluable |
| **G6** continuity diagnostic (reported, not a gate) | concordance of the Stage 2 label with the Stage 1 label for cells Stage 1 assigned to `Microglia_P2RY12_TMEM119`, `NK_cells` or `B_cells`; expected ≥75% for Microglia | reported; a large drop is flagged in the QC file and in the manuscript's limitations |

One additional pre-declared diagnostic: the number of cells moved off `NK` by the hard CD3 rule (§3.2), per sample.

---

## 5. Evaluability, primary type set, and the Stage 2 gate

### 5.1 Evaluability rule (unchanged from Stage 1)

A type is **evaluable** if its per-sample fraction is non-zero in **≥10 of the 75 snRNA samples** in the Stage 2 fraction table. This is the Stage 1 rule, applied unchanged.

### 5.2 Stage 2 primary type set — 9 types

`CD8_T`, `CD4_conv`, `Treg_FOXP3`, `NK`, `B_plasma` (B and plasma pooled, pre-declared), `Microglia_P2RY12_TMEM119`, `MoMF_CD163`, `DC`, `Neutrophils`.

**Explicit arithmetic, stated for transparency:** Stage 1's primary list had 10 entries with `B` and `Plasma` separate. Pooling them, as pre-declared above, yields **9**. Consequently Stage 2 multiplicity control is Benjamini–Hochberg over **K = 9**, not K = 10. Stage 1's K = 10 is unchanged and remains the quantity reported for Stage 1. This is a pre-declared consequence of the B/plasma decision, not an adjustment made after seeing a result.

### 5.3 Stage 2 pass/fail gate (pre-declared)

Stage 2 **PASSES** if and only if all of the following hold, evaluated before any interaction model is fitted:

1. gates G1, G3, G4, G5 hold for every type entering the primary set;
2. **≥6 of the 9** primary types are evaluable by the §5.1 rule;
3. `CD8_T`, `NK` and `Microglia_P2RY12_TMEM119` are each **individually evaluable**.

**On FAIL: stop. Report the QC table, report that the second-stage analysis could not be executed, and state a bounded null.** Do not adjust panels, thresholds, δ, resolutions or the model, and do not re-run. No Stage 2 interaction estimate is computed at all if the gate fails.

### 5.4 Residual risk, stated now

Fixing the annotation addresses **evaluability, not power**. The discovery design is 33 patients (20 Any_RT / 13 No_RT) and Stage 1 already reported it as underpowered. A correctly annotated T/NK compartment may still yield a null interaction. This plan is written so that a bounded null is a reportable outcome, not a failure of execution.

---

## 6. Estimands and denominators (pre-declared)

| Level | Definition | Transform |
|---|---|---|
| **Primary** (continuity with Stage 1) | per sample, `count(type) / count(all cells in sample)` | `arcsin_sqrt` |
| **Secondary** | per sample, `count(type) / count(immune cells in sample)` where "immune" = the Layer-1 label | `arcsin_sqrt` |

Both are reported. BH and the §5.3 gate apply to the primary estimand; the secondary is reported as a pre-declared companion with the same model, and is labelled as secondary wherever it appears.

---

## 7. Statistical model, multiplicity, permutation null, power

### 7.1 Model (identical to Stage 1 — no term added, removed or re-ordered)

```
arcsin_sqrt_fraction ~ C(timepoint, Treatment('primary')) * C(treatment_arm, Treatment('No_RT'))
                       + estimate_purity + C(who_grade_initial, Treatment('II'))
                       + age_at_diagnosis + C(codeletion_1p19q, Treatment('non_codel'))
                       + (1 | patient_id)
```

* `REML = False`; profile-likelihood confidence intervals with Wald fallback on failure.
* Optimiser chain: `lbfgs` (maxiter 500) → `bfgs` (2000) → `powell` (2000). The chain actually used is recorded per cell type.
* `mgmt_methylation` is **not** in the model (deviation v1.0 from the pre-spec; unchanged and re-declared here).
* The interaction term `timepoint:treatment_arm` is the only reported effect of interest.
* Implemented by reusing `scripts/lmm_interaction_pipeline.py` with `--method-tag snrna_direct`, `--fractions <stage2 table>` and `--out-suffix _stage2`, so no Stage 1 output file is written.

### 7.2 Multiplicity

Benjamini–Hochberg across the **K = 9** primary types (§5.2). Reporting threshold q < 0.05, as in Stage 1.

### 7.3 Permutation null

Identical operationalisation to Stage 1 (`scripts/lmm_permutation_null.py`): patient-level shuffles only (a patient keeps the same permuted arm at both timepoints), 10,000 permutations, seed 7331, two-tailed empirical p `(1 + Σ I(|β_perm| ≥ |β_obs|)) / (1 + N)`, applied to the Stage 2 primary tables with a new output suffix.

### 7.4 Power

Design unchanged from Stage 1 `scripts/post_data_power_check.py`: the same derived design (n = 33 patients, 20/13 split, 66 observations), the same ICC grid (0.40/0.60/0.70), the same f² anchors (0.28/0.40), nsim = 800.

Known defect, handled without editing Stage 1: the script derives its effective seed from Python's salted `hash(str(key))`, so the generator seed varies between processes while the audit JSON still records `"seed": 7331`. Stage 2 therefore sets **`PYTHONHASHSEED=0`** for the process, which makes the salted hash deterministic and the run reproducible, with no change to any script.

Reporting rule (STOP 2 decision): the power result is presented as an **approximate Monte-Carlo quantity with nsim disclosed** — e.g. "≈62% (nsim = 800)" — and never as a precise percentage. Macro-reproducible quantities (the design, the convergence count, the classification labels) may be quoted exactly. What the figure does and does not support is stated in the manuscript's limitations.

### 7.5 Smallest effect size of interest, and the equivalence-style statement

**SESOI (pre-declared):** a difference between arms in the primary-to-recurrence change in cell-type fraction of **5 percentage points (0.05 in fraction units)**. This is declared before any Stage 2 estimate is seen.

Conversion to the modelled scale, per type `t`, using the treatment-blind cohort-mean primary fraction `f̄_t` of the type (mean over all 75 samples, computed once, arm-blind):

```
SESOI_t = asin( sqrt( min(f̄_t + 0.05, 1) ) ) − asin( sqrt( f̄_t ) )
```

**Decision rule for interpreting each primary estimate** (95% CI `[lo, hi]`, interaction β on the `arcsin_sqrt` scale):

1. `lo > −SESOI_t` **and** `hi < +SESOI_t` → "the data are consistent with no clinically important interaction at the pre-declared threshold."
2. `[lo, hi]` lies entirely outside `[−SESOI_t, +SESOI_t]` and excludes 0 → "the estimate is at least as large as the pre-declared threshold of clinical relevance."
3. otherwise → "**inconclusive**: the confidence interval neither excludes the null nor lies within the pre-declared equivalence band."

Wording for outcomes 1 and 2 is observational only; no causal statement about radiotherapy is permitted anywhere (§11).

---

## 8. Sensitivity analyses

| Analysis | Rule |
|---|---|
| **Cell cap** | Primary analysis uses the Stage 1 cap (3,000 cells/sample, seed 20260917). Sensitivity: **1,000 cells/sample**, seed 20260917. A *no-cap* run is **not** performed: the 437,454 pre-cap QC-passed cells are not retained by any Stage 1 checkpoint, and re-merging raw matrices was declined for resource reasons. This is stated as a limitation, not silently omitted. |
| **δ (ambiguity band)** | annotation re-derived at δ = 0.025 and δ = 0.10 (§3.3); estimates reported as sensitivity only, never as primary. |
| **Estimand** | secondary immune-compartment fraction (§6) reported alongside the primary. |
| **Chemotherapy** | **Pre-declared as not estimable.** The source metadata records a 4-level exposure breakdown (RT-only n = 3; RT+chemo confirmed n = 16; RT+chemo unconfirmed n = 1; chemo-only n = 6; untreated n = 7), so a within-arm chemotherapy contrast cannot separate radiotherapy from concurrent alkylating chemotherapy at n = 3. Pre-declared rule: such a contrast is run only if the RT-only subgroup reaches n ≥ 5; it does not. The 4-level breakdown is reported descriptively and the RT/chemo co-exposure is stated as an unaddressed confounder. |
| **Purity adjustment** | `estimate_purity` remains in the model (Stage 1 specification). No purity-stratified re-analysis is added. |

---

## 9. Exploratory arms — relabelling only

GLASS (bulk deconvolution) and CGGA are **exploratory** and are **not re-run and not re-annotated**. The Stage 2 annotation panel is a discovery snRNA panel and does not transfer to bulk deconvolution, so no new GLASS/CGGA computation is claimed to follow from Stage 2.

* GLASS: one short paragraph plus one supplementary table/figure, reporting the direction of effect descriptively. Its published null/discordant behaviour against discovery is not evidence about the Stage 2 estimates.
* CGGA: one short paragraph, cross-sectional, no interaction claim.
* Sign-agreement counting stays **retired** (it was removed before this revision and is not reintroduced; see `audit_report.md` §1.1 for the stale auxiliary field that made this decision concrete).
* No "significant" language may be applied to either exploratory arm.

---

## 10. Outputs and provenance

New files (all additive; no Stage 1 file is overwritten):

| Artefact | Path |
|---|---|
| per-cell Stage 2 labels (method A + B, ambiguity flags) | `outputs/intermediate/stage2_cell_labels.parquet` |
| Stage 2 fraction tables (primary, secondary, cap-1,000, δ variants) | `outputs/intermediate/stage2_proportions_per_sample*.csv` |
| annotation QC report | `stage2_annotation_qc.md` |
| QC machine-readable record | `outputs/intermediate/stage2_annotation_qc.json` |
| interaction results | `outputs/intermediate/lmm_interaction_summary_stage2.csv`, `lmm_all_coefficients_stage2.csv`, `lmm_summary_stage2.json` |
| permutation null | `outputs/intermediate/snrna_direct_perm_pvalues_stage2.csv` |
| power | `outputs/intermediate/post_data_power_table_stage2.csv` (+ audit JSON) |
| figures/tables | regenerated under `submission_JNO/` |

Every new numeric value is written into `submission/supplementary/results_numbers.json` under a `stage2.*` key, and a corresponding row is added to `numbers_provenance.csv`. **Stage 1 values keep their existing keys and are not edited**; the two stages are separated by key prefix (`stage1.*` / `stage2.*`).

Recorded per model run: the seed, the optimiser actually used, the convergence flag, and any dropped term.

---

## 11. Reporting language (fixed wording rules)

* **Stage 1** is always described as **"pre-specified (locked 2026-09-13; SHA-256 prefix `b7ba4c7bae5c62bd`)"**.
* **Stage 2** is always described as a **"second-stage analysis; its plan was locked on 2026-10-21 (SHA-256 in `stage2_plan.lock.json`), after Stage 1 results were known"**. The phrase "pre-specified" must never be applied to Stage 2.
* No causal language about radiotherapy anywhere: "caused", "induced by radiotherapy", "RT-driven", "radiotherapy effect" are prohibited. The exposure is *any-RT exposure including concurrent chemotherapy* versus none.
* "Significant" may be used **only** for Stage 1 primary results with the test named. For exploratory and Stage 2 results, effects are reported as estimates with confidence intervals, plus the §7.5 decision-rule wording.
* "Validated" and "confirmed" are not used.
* Every reported result carries exactly one label: `PRIMARY` / `STAGE 2` / `EXPLORATORY` / `DIRECTION-ONLY` / `SENSITIVITY`.

---

## 12. Deviation policy

Any departure from this document requires a dated addendum in this directory with its own SHA-256, written **before** the affected computation. Adjusting panels, genes, thresholds, δ, the model, the gate or the multiplicity rule **after** a Stage 2 interaction estimate has been computed or viewed is prohibited: if a Stage 2 result is inconvenient, the correct action is to report it, not to retune the plan.

---

## 13. Sign-off

Execution of §3–§8 is authorised only after this file and its lock record have been reviewed.

* Reviewer: ______________________  Date: ______________
