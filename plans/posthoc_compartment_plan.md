# Post-hoc compartment-level analysis of the Stage 2 annotation

**Status:** POST-HOC. Locked and hashed before any compartment-level interaction estimate,
permutation p-value or equivalence classification was computed or viewed.

**Authority:** external review, `aplan.md` Step 1 ("fix the measurement"), and the decision
recorded in the submission README: the terminal-type annotation is not carried by the data at
the level at which the reviewer needs it, so the estimand is moved to the annotation level that
the pre-declared quality-control gates actually certified.

**Supersedes nothing.** The Stage 1 pre-specification (v1.0.0, SHA-256 prefix `b7ba4c7bae5c62bd`),
the Stage 2 plan (`006d3e47…1ced5e48`), its addendum and Amendment v1.2 (`c5295238…3fde1e`) are
unchanged and remain the authority for every estimate they already produced. No existing Stage 2
estimate, K = 9 family, panel, gene, threshold, δ, resolution, SESOI or model term is touched.

---

## 1. Why this analysis exists

The Stage 2 annotation is hierarchical. Each level was gated separately, and the gate outcome is
recorded in `stage2_annotation_qc.json`:

| Level | What it resolves | Pre-declared gate | Outcome |
|---|---|---|---|
| L1 | broad compartment / lineage family | κ ≥ 0.60 | κ = 0.9525 — **passed** |
| L2 | immune sub-compartment (T/NK, myeloid, B/plasma) | κ ≥ 0.50 | κ = 0.9959 — **passed** |
| L3/L4 | terminal type (CD8, CD4, Treg, NK, microglia, MoMF, DC, neutrophil) | AUROC ≥ 0.70, own-panel separation, ambiguity ≤ 20% | **failed for 5 of 9 types** |

The five terminal types whose specificity gate failed are exactly the types separated by markers
that the T/NK/B lineage shares. Their fraction is therefore a noisy read-out, and a narrow
confidence interval computed on a noisy read-out does not establish equivalence: classical
measurement error attenuates the coefficient towards zero and does not, by itself, widen the
interval enough to be visible. Reporting "no interaction larger than 5 percentage points" from
those populations invites the objection that the label, not the biology, produced the null.

The compartment level is the only level of the hierarchy whose agreement gate passed for every
unit it resolves. Its members are the *sums* of the types that compete with each other, so a
misassignment between two members does not change the compartment total. The compartment
estimand is therefore the one the annotation supports.

## 2. Estimand

Three disjoint, exhaustive compartments over the 14 annotation labels of `ALL_LABELS`:

| Compartment | Members |
|---|---|
| `LYMPHOID` | CD8_T, CD4_conv, Treg_FOXP3, NK, B_plasma |
| `MYELOID` | Microglia_P2RY12_TMEM119, MoMF_CD163, DC, Neutrophils |
| `NEU_NONIMMUNE` | Malignant_astrocytes_oligo, Astrocytes_non_malignant, Oligodendrocytes_OPC, Endothelial_cells, Pericytes_Fibroblasts |

`NEU_NONIMMUNE` is carried as a **negative control** (an orthogonality test): the exposure is
applied to the patient, not to the tumour compartment, so a compartment that is not the object of
the hypothesis should not move. It is reported without multiplicity correction and no efficacy
claim is made from it. The Benjamini–Hochberg family is **K = 2** (`LYMPHOID`, `MYELOID`).

## 3. Aggregation procedure (unchanged from the Stage 2 run)

Identical to `fractions()` in `scripts/stage2_annotate.py`, applied to the compartments:

* numerator: cells with `methodA_ambiguous_d0050 == 0`, counted per label and summed within the
  compartment (δ = 0.05 is unchanged);
* denominator: all cells of the sample (`denom_kind = "all"`);
* ambiguous cells are excluded from every numerator and retained in the denominator;
* per-sample compartment fraction = compartment numerator / sample denominator.

Because each cell carries exactly one label and the denominator is shared, the compartment
fraction is exactly the sum of its members' fraction columns in
`stage2_proportions_per_sample.csv`. The runner asserts this identity to < 1e-12 before fitting,
so the compartment table is a re-aggregation of the frozen Stage 2 output and not a new
annotation.

## 4. Model, multiplicity, SESOI and decision rule (verbatim reuse)

* model: `FORMULA_FIXED` from `scripts/lmm_interaction_pipeline.py`, imported not rewritten;
  `arcsin_sqrt_fraction ~ timepoint * treatment_arm + purity + grade + age + 1p19q + (1|patient)`;
  REML = False; profile-likelihood CI with Wald fallback; optimiser chain lbfgs → bfgs → powell.
* evaluability: unchanged — non-zero in ≥ 10 of 75 samples.
* SESOI: the §7.5 formula applied to each compartment,
  `SESOI_c = asin(sqrt(min(f̄_c + 0.05, 1))) − asin(sqrt(f̄_c))`, with `f̄_c` the arm-blind
  mean over all 75 samples. The 5 percentage-point clinical threshold is unchanged.
* decision rule: §7.5 three-way rule on the 95% CI.
* multiplicity: Benjamini–Hochberg over K = 2, q < 0.05.
* label on every output: POST-HOC (STAGE 2 annotation, post-hoc aggregation level).

## 5. What this analysis is not

* It is not a rescue of the terminal-type result. It is a lower-resolution estimand that is
  answerable with the annotation that exists. Both are reported.
* It is not pre-specified. It was designed after the Stage 1 and Stage 2 estimates were known,
  and it is labelled POST-HOC wherever it appears.
* It does not change the paired design, the cohort, the exposure contrast (any RT-containing
  treatment vs none), the covariates or the equivalence margin.

## 6. Honest-reporting clause

The result is reported as computed, in whichever direction it falls. If a compartment shows an
interaction at or beyond the ±5-percentage-point band, that is reported as a positive finding
with the exposure misclassification (radiotherapy versus chemotherapy-containing treatment)
stated as its principal limitation. No outcome-driven change to this plan, to K, to the SESOI or
to the model is permitted after the estimate has been seen; any further change requires a dated,
hashed amendment written before the affected computation.

## 7. Outputs

```
outputs/intermediate/stage2_compartment_proportions_per_sample.csv
outputs/intermediate/lmm_interaction_summary_compartment_posthoc.csv
outputs/intermediate/jno/compartment_posthoc_report.json
```

Registry keys are written under `compartment.posthoc.*` in
`submission_JNO/supplementary/results_numbers.json` with per-value provenance rows.
