# Immune-cell composition in paired primary and recurrent IDH-mutant glioma by prior radiotherapy

Analysis code, the **locked analysis plans**, and the small derived result files behind the
manuscript submitted to the *Journal of Neuro-Oncology*.

This repository is the code-and-provenance companion to the paper. It contains the scripts that
produced the reported numbers, tables and figures; the two locked analysis plans in the **exact
bytes** that the manuscript's Table S2 and Table S3 cite by SHA-256; and the derived result files
that back every value in the text. It deliberately does not contain the manuscript text or any raw
or intermediate data.

---

## 1. What the analysis does

The question is whether prior radiotherapy (RT) is associated with a change in immune-cell
composition between paired primary and recurrent IDH-mutant glioma, and if so in which compartment.

The design is two-stage and the stages are kept separate on purpose:

| Stage | Plan document | Locked | Role |
|---|---|---|---|
| Stage 1 | [`plans/sec4_prespec.md`](plans/sec4_prespec.md) (+ `.json`) | 2026-09-13T12:29:40Z | The pre-specified analysis: cell-type list, model, multiplicity correction, permutation null. |
| Stage 2 | [`plans/stage2_plan.md`](plans/stage2_plan.md) (+ addendum 1, amendment 2) | 2026-10-21T08:57:05Z | A **post-hoc, second-stage** re-annotation and re-analysis. It is *not* pre-specified and does not claim to be. |

Both plans were locked before the estimates they govern were computed or viewed. The lock records
(`plans/*.lock.json`) carry the SHA-256 of the plan file they lock and the UTC timestamp.

The primary model, identical in both stages, is:

```text
arcsin_sqrt(fraction) ~ timepoint * treatment_arm
                        + purity + grade + age + 1p19q
                        + (1 | patient_id)
```

Multiplicity is controlled by Benjamini-Hochberg over the pre-declared family (K = 9 in Stage 2,
K = 10 in Stage 1). The permutation null is 10,000 patient-level treatment-label shuffles with
seed 7331. The smallest effect size of interest is 5 percentage points on the fraction scale.

## 2. Repository contents

| Path | What it is |
|---|---|
| `code/` | 25 Python scripts: acquisition, harmonisation, the LMM pipelines, the external-cohort checks and the submission builders. This is the code **as run**. |
| `plans/` | The two locked plans and the locked post-hoc compartment plan, with their lock records. Bytes unchanged — see §5. |
| `results/` | `results_numbers.json` (the numeric registry) and `numbers_provenance.csv` (per-value provenance), plus `results/audits/` with the 10 raw audit files the registry was built from. |
| `tables/` | `Table1.csv` and `TableS1.csv`–`TableS15.csv` — the derived tables as machine-readable CSV. |
| `requirements.txt` | The frozen dependency list of the analysis environment (Python 3.11). |

**No script here is required to read any table in this repository.** The tables, the registry and
the provenance file are produced by the scripts and are shipped as outputs so that each reported
value can be traced back to the artefact it came from.

### Path convention in `code/`

The scripts are the code as run, with two substitutions applied for redistribution: every local
absolute path was replaced by `<PROJECT_ROOT>`, and the machine account name by `[USER]`. They
therefore expect the project's own directory layout, with inputs under `data/` and intermediates
under `outputs/`. Those directories are **not** part of this repository.

## 3. Data you must obtain

Raw and intermediate data are not redistributed here. The scripts expect:

| Expected location | Source | Fetched by |
|---|---|---|
| `data/discovery/GSE326221/GSE326221_RAW_extracted/` | **GEO accession [GSE326221](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE326221)** — matched single-nucleus RNA-seq of paired primary/recurrent IDH-mutant glioma, 75 samples (59 in 10x MTX triplet format, 16 in CellRanger HDF5 format) | `code/acquire_snrnaseq_discovery.py` |
| `data/discovery/GSE326221_series_matrix.txt.gz` | GEO GSE326221 series matrix (sample-level GSM order and timepoint titles) | `code/acquire_snrnaseq_discovery.py` |
| `data/discovery/` bulk matrices | **GEO accession [GSE324481](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE324481)** — processed bulk RNA-seq matrices, GEO open access | `code/acquire_bulk_discovery.py` |
| `data/discovery/Johnson_2026_Nature_SupplTable*.xlsx` | Supplementary tables of Johnson KC, Spitzer A, Varn FS, et al, *Nature* 2026, doi:10.1038/s41586-026-10612-6 — patient-level treatment and clinical annotation | `code/acquire_clinical_metadata.py` |
| `data/cgga/CGGA.mRNAseq_693_*` | CGGA mRNAseq_693 cohort (external, cross-sectional check) | `code/cgga_crosssectional_pipeline.py` |
| `data/glass/` | GLASS-NL cohort files (external, longitudinal direction check) | `code/deconvolution_glass_pipeline.py` |
| `data/discovery/checkpoints/ckpt_05_annotated.h5ad` | Derived by `code/acquire_snrnaseq_discovery.py` | — |

The GEO accessions are open access; no controlled-access application is required for the
processed matrices used here. The cell-level objects (`*.h5ad`) are large and are not distributed.

## 4. Environment

```bash
python -m venv .venv
.venv\Scripts\activate          # Windows
pip install -r requirements.txt
```

Python 3.11. `code/environment_check.py` prints the resolved versions of the packages the
pipelines depend on. The versions quoted in the manuscript's software-inventory table
(Table S13) are the ones in `requirements.txt`.

## 5. Run order

The stage labels below are the ones written in each script's own header where a header carries
one (Steps 0a–0c, Phase 2.0/2.1, Phase 3.1/3.2, the submission Stages 1–2).

**Step 0 — acquisition and harmonisation**

1. `code/acquire_bulk_discovery.py` — GSE324481 processed matrices and per-sample metadata.
2. `code/acquire_clinical_metadata.py` — harmonises treatment status into `treatment_arm`
   (`Any_RT` vs `No_RT`) from the Johnson et al. supplementary tables.
3. `code/acquire_snrnaseq_discovery.py` — GSE326221 → per-sample AnnData, QC, per-sample
   cell-type proportions. Supports `--pilot N` for a fast import/syntax check.

**Pre-specification — must precede any analysis**

4. `code/write_sec4_prespec.py` — writes `outputs/sec4_prespec.json` / `.md` (v1.0.0) and the
   file's own SHA-256. The `acquire_*` scripts abort if this file is absent, so the plan cannot
   be written after the data are seen. The resulting file is the one copied into `plans/`.

**Phase 2 — cohort assembly and covariates**

5. `code/build_clinical_canonical.py` — `clinical_canonical.csv` and `estimate_scores.csv`.
6. `code/estimate_purity.py` — ESTIMATE (Yoshihara et al. 2013) purity, used as a covariate.

**Phase 3 — primary analysis**

7. `code/lmm_interaction_pipeline.py` — the interaction LMM, one fit per cell type.
8. `code/lmm_permutation_null.py` — the 10,000-shuffle patient-level permutation null.
9. `code/post_data_power_check.py` — observed-design power / minimum detectable effect.
10. `code/immune_subclust_sensitivity.py` — immune-only subclustering sensitivity.
11. `code/plot_forest_interaction.py` — the forest plot of the interaction estimates.

**External cohorts and method checks**

12. `code/cgga_crosssectional_pipeline.py` — cross-sectional CGGA check.
13. `code/deconvolution_custom_cibersortx.py`, `code/deconvolution_robustness.py` — custom
    signature deconvolution and its robustness variants.
14. `code/deconvolution_glass_pipeline.py`, `code/glass_direction_consistency.py` — GLASS-NL
    deconvolution and the direction-consistency check.

**Submission artefacts** (these read aggregate files only; no pipeline is re-run)

15. `code/submission_stage0_audit.py` — read-only audit of the project state.
16. `code/submission_build_numbers.py` — Stage 1: derives every reported value from the existing
    artefacts into `results/results_numbers.json` and `results/numbers_provenance.csv`. It exits
    non-zero if any of its pre-declared sanity sentinels disagree.
17. `code/submission_make_figures.py` — Stage 2: the figure files and legends.
18. `code/submission_make_tables.py` — the `Table1` / `TableS1`–`TableS15` files in `tables/`.
19. `code/submission_build_manuscript.py`, `code/build_manuscript_skeleton.py` — manuscript text.
20. `code/submission_build_package.py` — assembles the submission package.
21. `code/submission_qa.py` — the consistency checker: every numeral in the text is matched
    against the registry in `results/results_numbers.json`.

## 6. The locked plans, and how to verify them

The plan files in `plans/` are byte-identical to the ones cited in the manuscript; only their
directory differs (the manuscript refers to them as `submission/supplementary/…`).
Verify with:

```powershell
Get-ChildItem plans\sec4_prespec.json, plans\sec4_prespec.md, plans\stage2_plan.md |
    Get-FileHash -Algorithm SHA256 | Format-List
```

```bash
sha256sum plans/sec4_prespec.json plans/sec4_prespec.md plans/stage2_plan.md
```

| File | SHA-256 | Bytes | Cited in |
|---|---|---|---|
| `plans/sec4_prespec.json` | `B7BA4C7BAE5C62BD8C9B66ECAAB0E4ED87EBD0EF0AA3880E584E4B5FEB69F399` | 7,721 | Table S2 |
| `plans/sec4_prespec.md` | `7D44034E8FBA1DCA76CEB2DAB8B0C42F17E347CBF71B8C94793DA8ECB0B5BB58` | 6,548 | Table S2 |
| `plans/stage2_plan.md` | `006D3E47198C8D9BD51170B8AC1AE3E9D9F4C1A2CBF9F2D8155560691CED5E48` | 22,355 | Table S3 |
| `plans/stage2_plan.addendum_1.md` | `882A8C6F98C962B2C281FC4E7A97F764B5C218964CD7D6E53AC8D8E121D0E344` | 3,525 | — |
| `plans/stage2_plan.amendment_2.md` | `C52952385895D260F916469ABA6FF2E3C617F79B38F91B74D161CD31A5BFDE1E` | 12,707 | — |
| `plans/posthoc_compartment_plan.md` | `FAC49FCF41BEF4BE9CD352F64A29BF72313AA26D4C2716652AC2ED7379474A44` | 6,409 | — |

The `plans/*.lock.json` files are the lock records: each states the UTC timestamp at which the
plan was frozen and reproduces the SHA-256 of the file it locks. Decomposing the hash against the
registry in `results/numbers_provenance.csv` therefore gives a complete chain from the frozen
plan to each reported estimate.

## 7. What is deliberately not here

* The manuscript text, the cover letter and the supplementary document.
* The raw and intermediate data (`data/`, `outputs/`) and every `*.h5ad` object.
* Internal working notes, decision logs and session records.
* No credentials or tokens of any kind appear in this repository, and no local absolute path does
  either (see the `code/` path convention in §2).

## 8. Licence

MIT — see [LICENSE](LICENSE).

## 9. Citation

If you use this code, please cite the manuscript (reference to be added on publication) and the
discovery dataset:

> Johnson KC, Spitzer A, Varn FS, et al (2026) Acquired genetic and cell-state changes in
> IDH-mutant glioma progression. *Nature* 655:1048–1059. https://doi.org/10.1038/s41586-026-10612-6
