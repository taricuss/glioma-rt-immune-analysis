# Amendment 2 to `stage2_plan.md` — authorising the second-stage interaction model

**Applies to:** `submission/supplementary/stage2_plan.md`, SHA-256 `006d3e47198c8d9bd51170b8ac1ae3e9d9f4c1a2cbf9f2d8155560691ced5e48`, locked 2026-10-21T08:57:05Z.
**Written (UTC):** 2026-10-21. The exact SHA-256 of this file and its UTC lock time are in `stage2_plan.amendment_2.lock.json`.
**Written before:** any Stage 2 interaction estimate, permutation p-value, power figure or equivalence classification had been computed or viewed. No such quantity exists on disk as of this writing.
**Companions:** `stage2_plan.addendum_1.md` (operational only), `submission/supplementary/stage2_annotation_qc.md` (the gate outcome of record).

---

## 1. Why this is an amendment and not another addendum

`stage2_plan.addendum_1.md` recorded three purely operational departures and changed no scientific parameter. **This document changes one**: the stop clause of §5.3. Per the plan's §12 deviation policy, a change of that kind requires a dated, hashed instrument written before the affected computation. That is this file. It is disclosed as an amendment, not folded into the plan silently.

---

## 2. What was known at the time of writing (full disclosure)

| Item | Seen? | Comment |
|---|---|---|
| Stage 1 (v1.0) interaction estimates | **Yes** | Locked, reported and unchanged by this amendment. Plan §1 already states that this stage is post-hoc and was designed with Stage 1 results known |
| Stage 2 annotation QC gate outcome (G1–G5, §5.3) | **Yes** | It is an *annotation* diagnostic, already on disk in `stage2_annotation_qc.md`, and it is the reason this amendment exists. Concealing it would be neither possible nor useful |
| Any Stage 2 interaction estimate | **No** | None exists. No model has been fitted on the Stage 2 fraction tables |

This is precisely the situation `JNOplan.md` Phase 0 anticipates in its own words: *"the honest framing is a measurement correction after a failed evaluability gate. It is not a fresh pre-registration."* The gate outcome is known; the estimand is not.

---

## 3. The single change: §5.3's stop clause

§5.3 as locked reads: *"On FAIL: stop. Report the QC table, report that the second-stage analysis could not be executed, and state a bounded null. … No Stage 2 interaction estimate is computed at all if the gate fails."*

**That clause is superseded as follows. Conditions 2 and 3 of §5.3 remain binding gates. Condition 1 is reclassified from a study-level stop to a per-type annotation caveat, reported next to every estimate for the affected type.**

Conditions 2 and 3, checked against the annotation QC of record, both hold:

| §5.3 condition | Rule | Outcome |
|---|---|---|
| 2 | ≥6 of the 9 primary types evaluable by the §5.1 rule | **9 of 9** — passes under all five pre-declared variant definitions (primary, δ = 0.025, δ = 0.10, 1,000-cell cap, secondary denominator); see `outputs/intermediate/stage2_evaluability_variants.json` |
| 3 | `CD8_T`, `NK`, `Microglia_P2RY12_TMEM119` each individually evaluable | **all three evaluable** |

Condition 1 is where the gate failed:

| Type | Failed gate(s) | Consequence under this amendment |
|---|---|---|
| `CD8_T` | G4(i), G4(ii), G5 | estimated and reported **with a marker-specificity and ambiguity caveat** |
| `Treg_FOXP3` | G4(i), G4(ii) | estimated and reported **with a marker-specificity caveat** |
| `NK` | G4(i), G4(ii) | estimated and reported **with a marker-specificity caveat** |
| `B_plasma` | G4(i), G4(ii) | estimated and reported **with a marker-specificity caveat** |
| `Neutrophils` | G5 | estimated and reported **with an ambiguity caveat** |
| `CD4_conv`, `Microglia_P2RY12_TMEM119`, `MoMF_CD163`, `DC` | none | reported without caveat |

### 3.1 Grounds, fixed in writing before the estimate exists

1. **The plan's own evaluability rule passes with margin** (9 of 9). `JNOplan.md`'s no-go point is defined on evaluability — *"If fewer than 6 of 10 types are evaluable, stop investing"* — and is therefore **not triggered**.
2. **The failures are properties of marker panels, not of the data at the estimand.** G4(i)/G4(ii) test whether a panel separates its own cells from neighbours; G5 tests the ambiguity fraction. Four non-specific panels do not make the cohort's composition data inadmissible, and suppressing every estimand for the whole study because of them discards admissible information about the five types that pass and about the lymphoid and myeloid compartments as wholes.
3. **The two estimates the manuscript already reports are untouched.** `NK` and `Microglia_P2RY12_TMEM119` — the two Stage 1 could estimate — pass G5, and Microglia passes G4 as well. This amendment does not revise a single number already in the manuscript; it adds estimates that did not previously exist.
4. **The alternative is indefensible in the other direction.** Holding §5.3 leaves the study's central question permanently untestable on data that can answer part of it, while a locked, treatment-blind, independently cross-checked annotation sits unused.

### 3.2 The cost, stated plainly

This is a **relaxation of a pre-declared gate, decided after that gate's outcome was seen** and after the Stage 1 estimates were seen. It is disclosed in the manuscript (Methods and Limitations), in the cover letter, and in this file. **No quantity produced under this amendment may be described as pre-specified.** Every Stage 2 quantity carries the `STAGE 2` label, and every type listed in §3 above carries its named caveat wherever its estimate appears.

---

## 4. What is explicitly NOT changed

Nothing in this amendment retunes anything. The following are frozen exactly as locked, and may not be adjusted after the estimate is seen:

| Frozen item | Where locked |
|---|---|
| Panels and their gene sets | §3.2 |
| Gene space (12,314 genes), cell set (221,016 cells), cap (3,000/sample, seed 20260917), sample membership | §2, §8 |
| Ambiguity rule and δ = 0.05 (sensitivities 0.025 / 0.10) | §3.3, §8 |
| §5.1 evaluability rule | §5.1 |
| Conditions 2 and 3 of the §5.3 gate | §5.3, restated in §3 above |
| Estimands and denominators | §6 |
| The confirmatory model specification (§7.1), term for term | §7.1 |
| BH multiplicity over K = 9 at q < 0.05 | §7.2 |
| Permutation null: patient-level shuffles, 10,000 permutations, seed 7331, two-tailed empirical p | §7.3 |
| SESOI = 5 percentage points and its conversion SESOI_t | §7.5 |
| The three-way interpretation rule | §7.5 |
| Reporting-language rules | §11 |

---

## 5. Confirmatory specification, fixed here

**The confirmatory model is §7.1, verbatim and unchanged** — `arcsin_sqrt_fraction ~ timepoint * treatment_arm + estimate_purity + grade + age + 1p19q + (1 | patient_id)`, no term added, removed or re-ordered, fitted to the Stage 2 primary and secondary fraction tables via `scripts/lmm_interaction_pipeline.py` with `--method-tag snrna_direct --fractions <stage2 table> --out-suffix _stage2`.

*Why the model is not also changed.* `JNOplan.md` Phase 2 argues for a count-based model, and its argument (arcsin-square-root transforms ignore the closure of compositional data and the information in absolute cell numbers) is sound. Adopting it **as the primary** at the same time as relaxing the gate would mean changing two pre-declared things at once, and would also break the like-for-like comparison with Stage 1 that gives the two-stage design its point. It is therefore pre-declared here as an **alternative specification** (§6), reported alongside, never instead.

---

## 6. Pre-declared sensitivity battery

Fixed here, before any estimate is seen, so that none of it is an outcome-driven fork.

| # | Analysis | Definition | Label |
|---|---|---|---|
| S1 | Non-malignant denominator | per sample, `count(type) / count(all cells − Malignant_astrocytes_oligo)`; same model as §7.1 | `STAGE 2 SENSITIVITY` |
| S2 | Lineage aggregates | `T_NK = CD8_T + CD4_conv + Treg_FOXP3 + NK`; `Myeloid = Microglia + MoMF + DC`; each as a fraction of all cells and of immune cells; same model | `STAGE 2 SENSITIVITY` |
| S3 | Per-patient change | per patient, `Δ = logit(f_recurrent) − logit(f_primary)` with pseudo-count 1e-6, regressed on arm; Welch t test and a patient-level permutation (10,000, seed 7331) | `STAGE 2 SENSITIVITY` |
| S4 | Interval-adjusted | §7.1 plus `surgical_interval_mo_t1t2` as one additional fixed effect | `STAGE 2 SENSITIVITY` |
| S5 | Interval-stratified permutation | permutation performed within interval tertiles | `STAGE 2 SENSITIVITY` |
| S6 | Count-based GLMM | binomial / beta-binomial GLMM on per-sample counts with a patient random intercept and the same fixed-effect structure as §7.1 | `STAGE 2 SENSITIVITY` |
| S7 | Negative-control types | `Oligodendrocytes_OPC`, `Endothelial_cells`, `Astrocytes_non_malignant` — types that biology does not predict to respond to RT | `STAGE 2 SENSITIVITY` |
| S8 | Pseudo-count sensitivity | transforms recomputed at pseudo-counts 1e-4 and 1e-8 | `STAGE 2 SENSITIVITY` |
| S9 | Leave-one-patient-out influence | §7.1 refitted with each of the 33 patients removed; report max abs Δβ and the patient responsible | `STAGE 2 SENSITIVITY` |
| S10 | Batch, site and region | §7.1 plus `platform_10x` and `surgery_location` from the author supplement | `STAGE 2 SENSITIVITY` |

**Multiplicity.** The primary family is unchanged: BH over **K = 9** at q < 0.05 (§7.2). The sensitivity battery is reported **unadjusted** and labelled as such; it does not form a second family and may not be used to promote any type to "detected".

**No additional sensitivity analysis may be added later** except through a further dated, hashed amendment written before it is run.

---

## 7. Interpretation, unchanged, plus one pre-declared addition

The §7.5 three-way decision rule applies to every type, using the same per-type `SESOI_t`:

1. 95% CI `[lo, hi]` entirely inside `(−SESOI_t, +SESOI_t)` → consistent with no clinically important interaction at the pre-declared threshold;
2. 95% CI entirely outside that band and excluding 0 → at least as large as the pre-declared threshold;
3. otherwise → inconclusive.

**Pre-declared addition (S11, this amendment).** For completeness, a classical **TOST** is also reported per type: two one-sided tests at α = 0.05 against `±SESOI_t`, equivalently a 90% CI lying entirely inside the band, with both one-sided p-values tabulated. Where the conventional TOST (90% CI) and rule 1 above (95% CI) disagree, **both are reported and the disagreement is stated**; the stricter rule is never silently replaced by the more favourable one.

No causal language about radiotherapy is permitted anywhere (§11).

---

## 8. Disclosure obligations created by this amendment

1. **Manuscript** — Methods must state that the pre-declared Stage 2 gate failed on marker specificity, and Limitations must state the amendment and its per-type caveats. Every affected estimate carries its caveat in text and in the corresponding table.
2. **Cover letter** — one line naming the two-stage design and the amendment.
3. **Tracking files** — logged in `JNOplan_implementation_status.md` and referenced from `revisionplan_implementation_status.md`.
4. **Submission package** — this file and `stage2_plan.amendment_2.lock.json` are added to `submission_JNO/supplementary/`, so the reviewer sees the gate failure and the amendment together.
5. **Numbers registry** — every new Stage 2 value is registered in `results_numbers.json` with its provenance row before it appears in any manuscript text.

---

## 9. Residual risk, restated

Power is unchanged by this amendment. The registered design is 33 patients (20 any-RT / 13 no-RT), which Stage 1 already reported as underpowered. A correctly annotated T/NK compartment **may still yield an interval that neither excludes zero nor lies inside the equivalence band**. That is a reportable outcome under §7 rule 3, and this amendment does not make the study powered. It makes the estimand testable and its answer explicit.

---

## 10. Sign-off

| Item | Value |
|---|---|
| Authorising decision | User, 2026-10-21, on the J1 fork in `JNOplan_implementation_status.md` §9 — option **J1-A, "Amend to v1.2 and estimate"** |
| Stage 2 QC gate of record | **FAIL** — G4(i)/G4(ii) for `CD8_T`, `Treg_FOXP3`, `NK`, `B_plasma`; G5 for `CD8_T`, `Neutrophils` |
| Evaluability of record | **9 of 9** primary types (all five pre-declared variant definitions) |
| Stage 2 estimates computed before this file was hashed | **none** |
| Stage 1 estimates, pre-specification and hashes | unchanged by this amendment |
