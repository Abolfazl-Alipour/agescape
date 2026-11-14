# Step 4 muscle analyses: Perez complete; Kedlian acquisition incomplete

## Verdict and completion status

**Step 4 is incomplete.** The required 4.03-GB Kedlian raw H5AD was not completely downloaded, so no Kedlian expression values, count validation, eligibility outcomes, or biological models were produced. The smaller curated CELLxGENE object was not substituted, parallel downloading was not used, and the resumable partial file was retained. Perez was then completed sequentially as the lightweight second arm.

Perez provides whole-nucleus donor summaries and composition-confounded donor-pseudobulk gene models only. No released barcode-level annotations were found. No study pooling, cross-cohort numerical comparison, or independent-replication label was made. Non-significance is not equivalence.

## Frozen plan and resource guard

The design was frozen before expression outcomes in `data/manifest/muscle_analysis_plan.json` at HEAD `6763430`: donor inference, separate assays, six donors/two per group/20 cells, full-rank models, assay-wide multiplicity families, primary exclusion of Kedlian 470BR, a declared all-donor sensitivity, established P1 with expression-matched null, donor CLR, donor-pseudobulk gene models, and cross-fitted released-young Hellinger only under the nonnumeric-age sensitivity label.

Kedlian acquisition used two guarded sequential attempts with a 30-GB reserve and 5-GB cap. Only 154,009,600 of 4,025,846,616 bytes (3.827%) were acquired in approximately 2315 seconds. A bounded 4-MiB range diagnostic took 44.31 seconds. At the observed rate, more than 11 hours remained. The partial file was not opened or analyzed.

All 17 Perez H5 downloads completed sequentially through the same guard in 19.907 seconds (118,346,194 bytes). Perez analysis used a 6-GB memory cap with swap disabled; observed peak RSS was 0.603 GB and elapsed time was 7.138 seconds. No FASTQ or bulk count table was downloaded.

## Kedlian: not run

- Required raw-count value validation: **not run**, because the complete raw object is unavailable.
- P1, CLR, gene, Hellinger, primary adult, and all-donor models: **not run**, not null results.
- Numeric-age analysis: blocked because exact donor ages are not public.
- Partial local input: `../agescape-revision-data/muscle/kedlian/SKM_human_raw_cells2nuclei_2023-06-22.h5ad`.
- Audit: `results/followup-muscle-kedlian-download-resource.json` and `data/manifest/E-MTAB-13874.json`.

## Perez verified inputs

All 17 filtered 10x H5 files contained only the expected `matrix` group and passed streamed finite/nonnegative/integer validation. There are 58,309,907 stored int32 values, 111,145,052 total raw UMI, and exactly 134,486 released filtered nuclei from 6 young and 11 old donors. Feature order is identical across files. The 72-man bulk arm is excluded.

| GSM | Donor | Group | Exact age | Nuclei | Mean UMI |
|---|---|---|---:|---:|---:|
| GSM5098737 | HM1 | Old | 75 | 6,587 | 147.051 |
| GSM5098738 | HM2 | Old | 76 | 7,626 | 141.231 |
| GSM5098739 | HM3 | Old | 75 | 4,206 | 229.463 |
| GSM5098740 | HM4 | Old | 79 | 11,512 | 70.6357 |
| GSM5098741 | HM5 | Young | 20 | 5,877 | 168.487 |
| GSM5098742 | HM6 | Old | 82 | 7,311 | 142.456 |
| GSM5098743 | HM7 | Old | 72 | 4,478 | 164.623 |
| GSM5098744 | HM8 | Young | 21 | 5,510 | 144.785 |
| GSM5098745 | HM9 | Old | 70 | 9,552 | 157.539 |
| GSM5098746 | HM10 | Young | 21 | 4,739 | 234.734 |
| GSM5098747 | HM11 | Old | 75 | 3,991 | 208.418 |
| GSM5098748 | HM12 | Old | 72 | 5,342 | 233.826 |
| GSM5098749 | HM13 | Young | 19 | 3,623 | 169.29 |
| GSM5098750 | HM15 | Old | 79 | 3,537 | 135.296 |
| GSM5098751 | HM21 | Young | 22 | 14,685 | 2266.17 |
| GSM5098752 | HM22 | Young | 21 | 13,829 | 1697.95 |
| GSM5098753 | HM23 | Old | 90 | 22,081 | 1866.08 |

## Perez whole-nucleus donor-pseudobulk models

Counts were summed by donor before log1p(CPM). The primary model is old+mean log raw UMI; the predeclared sensitivity is exact age per 10 years+mean log raw UMI. Effects include 95% donor-level OLS intervals. BH is applied across all 20,343 genes separately within each model. These associations cannot distinguish cell-intrinsic changes from changing cell-type composition.

| Model | Tests | BH q<0.05 | Minimum BH q |
|---|---:|---:|---:|
| Old-young group, depth adjusted | 20,343 | 14 | 0.000236699 |
| Exact age per 10 years, depth adjusted | 20,343 | 8 | 0.000229834 |

The five smallest unadjusted p-values per model are shown for orientation; complete estimates are in `results/followup-muscle-perez-gene-models.csv`.

### exact_age_10y_sensitivity

| Model | Gene | Effect [95% CI] | p | BH q |
|---|---|---:|---:|---:|
| exact_age_10y_sensitivity | HCN1 | 0.395949 [0.32413, 0.467768] | 1.12979e-08 | 0.000229834 |
| exact_age_10y_sensitivity | CHODL | -0.130748 [-0.165265, -0.0962305] | 1.14406e-06 | 0.0116368 |
| exact_age_10y_sensitivity | UNC13C | 0.406707 [0.292793, 0.520621] | 2.27031e-06 | 0.0153949 |
| exact_age_10y_sensitivity | CARMN | -0.198095 [-0.257537, -0.138652] | 4.95403e-06 | 0.025195 |
| exact_age_10y_sensitivity | DNAH11 | 0.33149 [0.227502, 0.435478] | 8.0998e-06 | 0.0264154 |

### primary_group_depth_adjusted

| Model | Gene | Effect [95% CI] | p | BH q |
|---|---|---:|---:|---:|
| primary_group_depth_adjusted | HCN1 | 2.30709 [1.88765, 2.72652] | 1.16354e-08 | 0.000236699 |
| primary_group_depth_adjusted | CHODL | -0.78561 [-0.954303, -0.616917] | 9.47716e-08 | 0.000963969 |
| primary_group_depth_adjusted | UNC13C | 2.34602 [1.6544, 3.03763] | 4.0614e-06 | 0.0202401 |
| primary_group_depth_adjusted | DMRT2 | -1.24533 [-1.61461, -0.876059] | 4.33652e-06 | 0.0202401 |
| primary_group_depth_adjusted | CARMN | -1.15439 [-1.50093, -0.807862] | 4.9747e-06 | 0.0202401 |

## Perez blocked endpoints

| endpoint | status | outcomes_examined | reason |
|---|---|---|---|
| whole_nucleus_descriptive | completed | True | 17 donor summaries from released filtered nuclei. |
| whole_nucleus_donor_pseudobulk_gene_models | completed_composition_confounded | True | Group and exact-age sensitivity models; cell-composition changes cannot be separated from cell-intrinsic expression. |
| P1_transport | blocked_composition_confounded | False | No released per-nucleus annotation; no whole-nucleus transport claim. |
| cell_type_models | blocked_no_released_annotations | False | No annotations were inferred or created. |
| donor_CLR_composition | blocked_no_released_annotations | False | No cell-type counts are available. |
| Hellinger_displacement | blocked_no_released_annotations | False | No defensible cell-type-specific young reference. |
| bulk_72_man_arm | excluded | False | Bulk RNA-seq is outside single-nucleus inference and was not downloaded. |

The GEO series/per-donor inventories, all 17 H5 root schemas, primary paper and supplements, and public processed-object/repository searches were checked. No released barcode-level annotation mapping was found. Secondary studies' reannotations were not substituted. P1, cell-type models, CLR, and Hellinger were not computed.

## Cross-cohort status

No numerical cross-cohort comparison was possible: Kedlian has no estimates, and Perez whole-nucleus estimates would not be compatible with cell-type-specific intercostal-muscle estimands in any case. Perez reuse by Kedlian and unresolved individual disjointness prohibit an independent-replication label.

## Delivered files and reproduction

- `data/manifest/muscle_analysis_plan.json`: frozen design.
- `scripts/run_followup_kedlian.py`: complete backed Kedlian pipeline, ready but not run on the incomplete file; synthetic end-to-end tested.
- `scripts/download_followup_perez.py` and `scripts/run_followup_perez.py`: guarded sequential acquisition, validation, models, and manifest hashes.
- `src/followup_muscle.py` and `tests/test_followup_muscle.py`: model/safeguard helpers and synthetic tests.
- `results/followup-muscle-perez-*`: donor summaries, 40,686 gene-model rows, endpoint status, download/analysis resource records, and verification.
- `results/followup-muscle-kedlian-download-resource.json`: incomplete acquisition audit; no Kedlian outcome CSV exists.
- `results/followup-muscle-family-summary.csv`, `-resource-summary.csv`, and `-cross-cohort-comparison.csv`: joint status tables.
- `data/manifest/E-MTAB-13874.json` and `GSE167186.json`: updated acquisition, validation, hash, analysis, and blocker records.

Resume Kedlian without deleting the partial file:

```bash
python scripts/resource_guard.py download --url https://cellgeni.cog.sanger.ac.uk/muscleageingcellatlas/raw/SKM_human_raw_cells2nuclei_2023-06-22.h5ad --output ../agescape-revision-data/muscle/kedlian/SKM_human_raw_cells2nuclei_2023-06-22.h5ad --min-free-gb 30 --max-file-gb 5
python scripts/resource_guard.py run --path ../agescape-revision-data/muscle --min-free-gb 30 --memory-gb 6 -- python scripts/run_followup_kedlian.py build-cache
python scripts/resource_guard.py run --path ../agescape-revision-data/muscle --min-free-gb 30 --memory-gb 6 -- python scripts/run_followup_kedlian.py analyze
python scripts/run_followup_kedlian.py finalize
python scripts/report_followup_muscle.py
```

Perez reproduction:

```bash
python scripts/download_followup_perez.py
python scripts/resource_guard.py run --path ../agescape-revision-data/muscle --min-free-gb 30 --memory-gb 6 -- python scripts/run_followup_perez.py run
python scripts/run_followup_perez.py finalize
```

## Validation record

- Full suite: **127 passed, 1 skipped**.
- Synthetic Kedlian and Perez workflows passed; no real Kedlian outcome is claimed.
- `git diff --check`: passed. No commit or push was made. Parent must decide whether to resume the long Kedlian transfer before accepting Step 4.

## Not attempted

No raw FASTQ, reclustering, numeric-age imputation, Perez bulk inference, alternate-object substitution, parallel download, study pooling, independent-replication claim, or PBMC file modification was performed.
