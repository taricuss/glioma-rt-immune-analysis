# Section 4.5 — Pre-specification  (LOCKED v1.0.0)

**Locked at (UTC):** 2026-09-13T12:29:40+00:00  
**Integrity sha256:** `e69cfd7aec4e5adec695a08edc6b5306b3336949ab54e6d34af6e2f2c10757da`  

> **Rule for any deviation later:** If any field in this file is altered 
> *after* partial results from deconvolution/statistics have been 
> inspected, every downstream figure or table that uses the changed 
> setting must be labelled **EXPLORATORY** in the manuscript.  Primary 
> claims come only from this locked v1.0.0 document.

---

## (a) Primary cell types tested (BH-FDR family K=10)

1. `CD8_T_cells`
2. `CD4_T_cells_conventional`
3. `T_regulatory_cells_FOXP3`
4. `NK_cells`
5. `B_cells`
6. `Plasma_cells`
7. `Microglia_P2RY12_TMEM119`
8. `Monocyte_derived_macrophages_CD163`
9. `Dendritic_cells`
10. `Neutrophils`

## (a) Exploratory cell types (raw p only, NEVER merged into primary BH-FDR)

1. `T_cells_exhausted_PDCD1_HAVCR2`
2. `M1_macrophages`
3. `M2_macrophages`
4. `Endothelial_cells`
5. `Fibroblasts_pericytes`

## (b) Exact model formulas

### Discovery cohort (Johnson CARE / Nature 2026)

```r
arcsin_sqrt_fraction ~ timepoint * treatment_arm   + estimate_purity   + who_grade_initial   + age_at_diagnosis   + codeletion_1p19q   + mgmt_methylation   + (1 | patient_id)
```

- Response transform: **arcsin_sqrt** (arcsin-square-root preferred over logit because many deconvolution fractions are exactly 0 in glioma).
- Hypothesis term: `timepoint:treatment_arm`
- Treatment: **BINARY — Any_RT (n=20) vs No_RT (n=13), exactly per §4.3 lock below.** (binary levels: ['No_RT', 'Any_RT']; timepoint levels: ['initial', 'recurrent']).
- Mandatory covariates (modeled as fixed effects every run): **estimate_purity, who_grade_initial, age_at_diagnosis, codeletion_1p19q, mgmt_methylation**.
- Random effects: `(1 | patient_id)`.
- Note: Do NOT add chemo as a within-arm covariate: RT-only n=3 cannot support that estimate (plan §4.3 hard rule).

### GLASS validation cohort

- Mirror formula; four-level treatment (untreated / TMZ_only / RT_only / RT_TMZ) used only if each non-empty cell has n ≥ 3.  Otherwise fall back to binary discovery definition (**If four-level per-cell n collapses below 3 in any arm, fall back to binary Any_RT/No_RT discovery definition.**).

## (c) Multiple-testing correction

- Family: CELL_TYPES_PRIMARY, K = **10** tests.
- Method: **Benjamini-Hochberg step-up (BH-FDR)**.
- α (family-wide): 0.05.  Rejection convention: **q-value < 0.05 per cell type within the primary family.**.
- Exploratory policy: CELL_TYPES_EXPLORATORY are a separate family.  Never merge them into the primary BH-FDR pass; report raw p + 'EXPLORATORY' label.
- Robustness: Bonferroni-adjusted p (0.05 / 10 = 0.005) reported alongside BH q as a sensitivity row.

## (d) Permutation-based empirical null for the interaction coefficient

- Units shuffled: patient-level treatment labels (timepoints stay paired within patient).
- Iterations: **10000**, seed = `7331`.
- How: For each of 10000 iterations, sample an independent permutation of the length-N_patients treatment vector WITHOUT stratification (no block preservation by subtype — conservative; subtype IS a covariate in the model so residual subtype-balance-drift is already modelled). Fit the exact same LMM formula each time, extract the interaction coefficient, build the two-tailed empirical CDF, compute empirical p = (1 + count(|perm_coef| >= |observed_coef|)) / (1 + n_perm).
- Reporting: For every primary cell type, four p-values appear side-by-side: (1) asymptotic LMM p, (2) LMM Satterthwaite-KE df adjusted p, (3) permutation-based empirical p, (4) BH-FDR q derived from (3). Primary claim uses (4); (1) and (2) shown for reviewers who ask.
- Parallelism: Use joblib loky backend, n_jobs = min(cpu_count, 8).

## (e) Mechanistic marker sets — CONDITIONAL (only if primary hit triggers them)

### T_cell_exhaustion
- Trigger: CD8_T_cells (q < 0.05 in primary test).
- Gene list: `PDCD1, HAVCR2, LAG3, TIGIT, CTLA4, TOX, LAYN`.
- Test design: At RECURRENT timepoint, Any_RT vs No_RT DESeq2/limma contrast restricted to CD8+ deconvolution-pseudo-bulk (or snRNA CD8+ cluster if available).

### T_regulatory_suppression
- Trigger: T_regulatory_cells_FOXP3 (q < 0.05 in primary test).
- Gene list: `FOXP3, IL2RA, CTLA4, TIGIT, IL10, TGFB1, IKZF2, IKZF4`.
- Test design: Same recurrent-only Any_RT vs No_RT design.

### Myeloid_polarization_M1_M2
- Trigger: ['Microglia_P2RY12_TMEM119', 'Monocyte_derived_macrophages_CD163'] (q < 0.05 in primary test).
- M1 list: `CD80, CD86, IL1B, IL6, TNF, NOS2, IRF5`.
- M2 list: `CD163, MRC1, MSR1, ARG1, IL10, TGFB1, CCL22, CHI3L1`.
- Test design: Module score (AddModuleScore in Seurat) per sample, Any_RT vs No_RT at recurrent timepoint, Mann–Whitney U.

### GSEA restriction (§5.7 hard rule)
- Allowed: MSigDB HALLMARK, ImmuneSigDB C7 (immunologic signatures) — restricted to sub-collections with keywords: glioma, lymphocyte, tcell, exhaustion, myeloid, radiotherapy..
- Forbidden: Never run the full MSigDB C2/C5/C7 unrestricted.  Never use KEGG as a primary pathway layer (too coarse).
## (f) Discovery-cohort treatment-label lock (plan §4.3)

- Total cohort: 35 patients.
- Excluded (unclear / NA treatment): **2 patients** (IDs populated from Johnson et al. Supp Table 1 by `scripts/acquire_clinical_metadata.py`; listed in every output that carries the primary fit).
- **Any-RT arm: n = 20** patients (RT-only=3, RT+chemo confirmed=16, RT+chemo chemo-unconfirmed=1).
- **No-RT arm: n = 13** patients (chemo-only=6, untreated=7).
- Rule never break: RT-only subgroup (n=3) must NEVER be reported as a standalone inferential test in the primary section.  Descriptive proportions only; any RT-only contrast goes to Appendix labelled 'n-too-small — descriptive only, no p-value reported'.

## (g) Forest plot — pre-specified visual

- X-axis: Interaction coefficient (timepoint × treatment), arcsin-sqrt-fraction units.
- Rows order: CELL_TYPES_PRIMARY order.
- Two panels: Discovery | GLASS.  Error bars: 95% profile-likelihood CI; Wald fallback with footnote if profile fails.
- Pooled diamond only if I² p > 0.10 (no significant cross-cohort heterogeneity).
- Annotations: BH-FDR q (from permutation-based p) on every discovery row.

---
*End of pre-spec v1.0.0.  Last edit MUST pre-date first `scripts/acquire_bulk_discovery.py` run; the latter script aborts if this file does not exist.*