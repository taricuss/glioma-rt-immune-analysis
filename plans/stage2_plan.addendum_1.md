# Addendum 1 to `stage2_plan.md`

**Applies to:** `submission/supplementary/stage2_plan.md`, SHA-256 `006d3e47198c8d9bd51170b8ac1ae3e9d9f4c1a2cbf9f2d8155560691ced5e48`, locked 2026-10-21T08:57:05Z.

**Written (UTC):** 2026-10-21T09:12:00Z — **before** any Stage 2 annotation, QC gate or interaction estimate was computed or viewed. No Stage 2 result existed at the time of writing.

Per §12 of the plan, this addendum records three departures. None of them changes a scientific parameter: no panel, gene, threshold, δ, gate, estimand, model term, multiplicity rule, SESOI or decision rule is altered.

## A1.1 Output format for per-cell labels

`stage2_cell_labels.parquet` → **`stage2_cell_labels.csv.gz`**.

Reason: `pyarrow` is not installed in the frozen Stage 1 environment (`.dgnc_env`; see `submission/code/requirements.txt`), so parquet cannot be written without installing a new dependency, which the project's dependency rule disallows. gzip-compressed CSV carries identical content. No other output path changes.

## A1.2 Operationalisation of the ambiguity flag used for numerator exclusion

The plan (§3.3) says ambiguity is evaluated "within each layer" and that ambiguous cells are excluded from every type's numerator. It does not say which layer's margin governs the exclusion when a cell passes through several layers.

Operationalisation, fixed here: **the margin used for numerator exclusion is the margin of the layer that decides the cell's reported label** — Layer 3 for `CD8_T` / `CD4_conv` / `Treg_FOXP3` / `NK`, Layer 4 for `Microglia_P2RY12_TMEM119` / `MoMF_CD163` / `DC` / `Neutrophils`, and Layer 2 for `B_plasma`. Layer 1 and Layer 2 margins are still computed and their full distributions are reported in `stage2_annotation_qc.md` as diagnostics, but a Layer 1 immune/non-immune tie does not by itself exclude a cell from a terminal-type numerator.

Reason: the terminal layer is the decision that actually assigns the reported type; a Layer 1 tie is already resolved by the fact that the cell's winning terminal layer is a T/NK or myeloid panel. The alternative reading (exclude on any layer's tie) is also reported in the QC file as a diagnostic so the consequence of the choice is visible.

## A1.3 Operationalisation of gate G4(i)

The plan (§4, G4) defines G4(i) as "for ≥70% of a type's cells, that type's own panel is the argmax within its layer".

Operationalisation, fixed here: G4(i) is **a rank-based one-vs-rest AUROC** — for each type, the type's own panel score used to separate its cells from all other cells in the smallest layer containing at least one competing type (L3 for the T/NK subtypes, L4 for the myeloid subtypes, L2 for `B_plasma`) — required **≥ 0.70**. The layer in which each type is tested is stated in the QC file.

Reason: the written form is identically 1.0 for method A by construction (a cell is assigned to the argmax of that layer), so it is not a test and can never fail. The AUROC form tests the same property — that a type's panel genuinely separates its cells — while being falsifiable. The argmax-consistency statistic is also reported, for completeness.

## A1.4 Note on the Stage 1 LMM script

`scripts/lmm_interaction_pipeline.py` is **not edited**. The Stage 2 driver imports it and rebinds only module-level constants (`PRIMARY_CELLS`, `EXPLORATORY_CELLS`, `ALL_CELLS`), passing `--out-suffix _stage2`, so the Stage 1 numerics that produced the locked outputs are provably untouched. Recorded in `stage2_annotation_qc.md`.
