# Follow-up Step 2: Tabula Sapiens 2.0 metadata gate

## Status

The supplementary-table metadata gate is complete. It used only the two released workbooks below. The abandoned 45 GB all-cells file was not downloaded or opened, and no expression analysis was run. All resulting tables are marked `supplementary_tables_only_no_H5AD_verification`; supplementary cell counts are not raw expression counts.

## Sources and scope

- `/home/reza/projects/agescape-revision-data/ts2/media-18.xlsx` (178,962 bytes; SHA-256 `40379d7c8dded32d4710b39c29ce2eeca65456f9d410daaa1d5815a3f952d07b`). The script reads `all`, `tissue_donor_cell-type`, `tissue`, and `donor`. These sheets supply aggregate tissue, compartment, cell-ontology label, method, donor, and number-of-cells fields.
- `/home/reza/projects/agescape-revision-data/ts2/media-24.xlsx` (38,437 bytes; SHA-256 `7a4606ebf5f77dd489083a4806862424e843053eecfd78c6dd5929ba36449ce3`). The script reads `Manually_Curated_metadata` and `ChatTS_metadata` for donor age and sex reconciliation. The manually curated sheet is the primary source.

The workbook field `method` distinguishes `10X` and `smartseq`. It is not treated as verification of the H5AD `assay` field or of chemistry.

## What is verified from the supplementary tables

### Counts and donor-organ links

- The `all` sheet contains 3,223 unique donor × tissue × compartment × cell-type × method rows, covering 24 donors and 28 tissues. Its 1,136,218 cells reconcile exactly to the donor-aggregated `tissue_donor_cell-type` sheet and to both the `tissue` and `donor` total sheets.
- Method-specific totals are 1,093,048 cells labelled `10X` and 43,170 labelled `smartseq`. These are released summary counts, not a check of matrix dimensions or raw-count contents.
- The donor-organ crosswalk contains 123 observed pairings. Method-specific cell totals sum exactly to the corresponding all-method totals.
- The cell-type coverage table retains method and compartment rather than pooling 10x and Smart-seq. Its `meets_20_cell_metadata_floor` flag is only a supplementary-count screen, not full analysis eligibility.

### Donor ages and sex

- Curated ages and sex are available for 22 of the 24 donors represented in `media-18.xlsx`: 12 female and 10 male donors, aged 22–74 years, with mean 49.91 and sample SD 14.17 years.
- TSP9 (Pancreas) and TSP20 (Stomach) occur in the count workbook but are absent from both donor-metadata sheets. Their ages and sex remain unresolved, so a complete 24-donor age distribution is not verified.
- Twenty-one curated ages agree exactly with `ChatTS_metadata`. TSP30 is 26 in `Manually_Curated_metadata` but 321 in `ChatTS_metadata`; the curated value 26 is retained and the conflict is explicit in the donor CSV. Sex agrees after normalization for all 22 donors with curated metadata.
- Among the 22 age-verified donors, four are older than 60 and six are aged 60 or older. This does not verify the planned statement about all 24 donors, and “over 60” must not be used interchangeably with “60 or older.”

### Tissue, cell-type, and organ-pair screens

- Age summaries use one curated age per donor and exclude TSP9 and TSP20. Organ summaries aggregate methods for description; analysis eligibility must still be evaluated separately by method and cell type.
- Applying the six-age-verified-donor, sample age-SD ≥12, and per-donor 20-cell supplementary floors separately within each method leaves one organ-cell-type candidate: 10x blood natural killer cells, with seven donors and age SD 12.06. This is only a metadata candidate for a later analysis, not a result.
- The organ-pair table contains 378 10x pairs and 300 Smart-seq pairs. None passes the organ-level screen of at least six shared age-verified donors and shared-donor age SD ≥12. Five 10x pairs have six shared age-verified donors, but all five fail the age-SD floor; Smart-seq has at most three shared donors for any pair.
- These pair results are organ-level screens. They do not establish that a shared cell type clears the donor and cell-count floors in both organs.

The supplementary workbook has a Skin summary with four donors (TSP2, TSP10, TSP14, TSP21; ages 61, 33, 59, and 45) and 17,786 cells. Four donors do not meet the six-donor floor. The local `skin.h5ad` is known four-donor Tabula Sapiens 1.0 data: it was not read or used, is not represented as TS2, and cannot authenticate a TS2 skin object. No TS1-to-TS2 independent replication is claimed.

## Generated tables

- `results/followup-ts2-donor-summary.csv`: 24 rows; curated age/sex reconciliation and supplementary donor totals.
- `results/followup-ts2-organ-summary.csv`: 28 rows; descriptive tissue totals, age-verified donor counts, and age dispersion.
- `results/followup-ts2-donor-organ-crosswalk.csv`: 123 rows; observed donor-organ links with all-method and method-specific supplementary cell totals.
- `results/followup-ts2-celltype-coverage.csv`: 3,223 rows; method- and compartment-resolved supplementary cell-type counts and the 20-cell metadata flag.
- `results/followup-ts2-organ-pair-eligibility.csv`: 678 rows; method-specific shared-donor and age-SD screens.

Columns derived from `media-18.xlsx` are labelled `supplementary` or `method_supplementary`. No assay, primary-data fraction, raw-count, or expression fields have been imputed.

## Unresolved fields and blocked analyses

The workbooks do not verify cell barcodes, duplicate cells across release objects, exact assay/chemistry, H5AD disease harmonization, anatomical position, enrichment or suspension type, raw-count layers, `is_primary_data`, per-cell QC/depth, or concordance between workbook and H5AD labels. Curated clinical columns in `media-24.xlsx` are not substitutes for a per-cell H5AD disease field. TSP9/TSP20 ages and sex also remain unresolved.

Consequently, no P1/off-program expression, Hellinger displacement, direct Shannon entropy, CLR composition, covariate-adjusted aging model, or exploratory RNA organ-age model was run. Any otherwise eligible within-organ analysis remains blocked pending small per-organ TS2 release objects that permit cell-level metadata reconciliation, deduplication, raw-count checks, QC, and donor-level inference. Smart-seq must remain separate from the primary 10x arm. Confirmatory paired-organ interaction analyses are already blocked by the supplementary shared-donor/age screen. Skin-specific aging analyses are blocked by both the four-donor floor and the absence of an authenticated small TS2 skin object. No plasma-proteomic organ-age claim is possible from these transcriptomic summaries.

## Execution and validation

The script is `scripts/prepare_ts2_metadata.py`. The final run used the existing `agescape` environment, a 512 MiB virtual-memory limit, and single-threaded numerical libraries:

```bash
ulimit -v 524288
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /usr/bin/time -v /home/reza/miniconda3/envs/agescape/bin/python \
  scripts/prepare_ts2_metadata.py
```

It exited 0 in 0.56 seconds with maximum resident set size 127,160 kB (about 124.2 MiB). An independent CSV audit rechecked row counts, key uniqueness, all five cell-total reconciliations, method-specific totals, the two missing donor identities, the TSP30 conflict, sex agreement, and the absence of passing organ-pair screens; all assertions passed. No large file was accessed, and no download was started.

## Blood 10x natural-killer-cell expression gate

### Public asset and access decision

The public [CELLxGENE collection API](https://api.cellxgene.cziscience.com/curation/v1/collections/e5f58829-1a66-40b5-a624-9046778e74f5) was queried on 2026-11-13. Collection `e5f58829-1a66-40b5-a624-9046778e74f5`, collection version `61d34250-a255-4c7c-b2f0-0e012fbe2e78`, contains one organ-specific blood dataset:

- Dataset: `983d5ec9-40e8-4512-9e65-a572a9c486cb`, titled `Tabula Sapiens - Blood`.
- Dataset version: `1bf55482-ad49-446a-87c6-f3442bebcef1`, published 2023-06-16 and revised 2026-06-11.
- Sole H5AD asset: `https://datasets.cellxgene.cziscience.com/1bf55482-ad49-446a-87c6-f3442bebcef1.h5ad`.
- Exact API/HEAD size: 2,568,875,081 bytes (2.568875081 decimal GB; 2.392451354 GiB), with 85,233 cells in the catalog.

This is an authenticated public per-organ object, but it is not small enough for the required maximum 2 GiB file. The guarded command

```bash
/home/reza/miniconda3/envs/agescape/bin/python scripts/resource_guard.py download \
  --url https://datasets.cellxgene.cziscience.com/1bf55482-ad49-446a-87c6-f3442bebcef1.h5ad \
  --output /home/reza/projects/agescape-revision-data/ts2/blood.h5ad \
  --min-free-gb 30 --max-file-gb 2
```

refused the transfer before creating an output file: `download size 2.39 GB exceeds maximum file size 2.00 GB` (exit 1). No blood H5AD was downloaded or opened. The 45 GB all-cells object and old `skin.h5ad` were not used.

The catalog lists donors TSP1, TSP2, TSP7, TSP8, TSP10, TSP14, TSP21, TSP25, and TSP27; assays 10x 3' v3, 10x 5' v2, and Smart-seq2; disease `normal`; and cell suspensions. These are catalog leads only. Without the H5AD, raw integer counts, the raw-count layer, `is_primary_data`, cross-object deduplication, per-cell donor/age/sex/assay/disease concordance, chemistry, and actual UMI depth remain unverified.

### Frozen donor family and endpoint disposition

Before any expression outcome was accessed, the 10x blood natural-killer-cell donor family was frozen from `results/followup-ts2-celltype-coverage.csv`:

| Donor | Age (years) | Sex | Supplementary cells | 20-cell screen |
|---|---:|---|---:|---|
| TSP1 | 59 | F | 71 | pass |
| TSP2 | 61 | F | 243 | pass |
| TSP7 | 69 | F | 1,641 | pass |
| TSP8 | 56 | M | 14 | fail |
| TSP10 | 33 | M | 285 | pass |
| TSP14 | 59 | M | 64 | pass |
| TSP21 | 45 | F | 250 | pass |
| TSP25 | 60 | M | 406 | pass |
| TSP27 | 56 | F | 11 | fail |

The seven supplementary candidates are TSP1, TSP2, TSP7, TSP10, TSP14, TSP21, and TSP25 (4 female, 3 male; ages 33–69, mean 55.14, sample SD 12.0614 years; 2,960 supplementary cells). Any analysis with these donors was prespecified as exploratory, with donors—not cells—as inferential units. These counts and covariates remain supplementary-table values, not H5AD-verified analysis eligibility.

`results/followup-ts2-blood-nk-endpoint-status.csv` freezes the intended designs and records every endpoint as **untestable**, not null:

- Donor-pseudobulk gene models: mean cell-level log1p(CPM) per donor and gene, continuous-age OLS with effects per 10 years and 95% confidence intervals; sex, mean log actual UMI, and nonconstant assay terms only with a full-rank design and positive residual degrees of freedom; BH across all successfully fitted variable genes with finite donor values. No cell-level inference was permitted.
- P1 transport: the established P1_Clean signature and expression-matched random-null procedure only after detected up/down genes and valid null pools were verified.
- Off-program expression: the established off-program definition and technical exclusions only after its gene/signature floors were verified.
- Hellinger displacement: the established leave-one-donor-out young-reference endpoint only with at least 50 remaining reference cells, at least 20 cells per donor, and identifiable donor-level covariates.

The asset-size block prevented all raw-count and metadata checks, so P1 detection, off-program, and Hellinger reference floors were **not evaluated and do not pass by assumption**. No pseudobulk matrix, gene model, effect size, confidence interval, p-value, BH q-value, P1/off-program score, or Hellinger distance was computed. There is no expression result from this arm.

### Reproducible blocked-arm deliverables

- `scripts/prepare_ts2_blood_nk_gate.py`: metadata-only, fail-closed script for the exact collection/dataset/version, size gate, donor freeze, endpoint disposition, and manifest. It never downloads or opens an H5AD.
- `data/manifest/tabula_sapiens_2_0_blood.json`: exact asset identity/size, catalog metadata, coverage-file SHA-256, access decision, and unverified gates.
- `results/followup-ts2-blood-nk-asset.csv`: one-row asset inventory.
- `results/followup-ts2-blood-nk-donor-eligibility.csv`: nine-donor supplementary screen with all expression-eligibility and H5AD-verification flags false.
- `results/followup-ts2-blood-nk-endpoint-status.csv`: four frozen endpoints, all untestable with `outcomes_examined=False`.
- `tests/test_ts2_blood_nk_gate.py`: six offline tests covering exact asset resolution, the inclusive size boundary, fail-closed version/duplicate handling, the seven-donor freeze, and absence of outcome fields.

The script completed with exit 0 and printed the 2.392 GiB block. The focused test command `/home/reza/miniconda3/envs/agescape/bin/python -m pytest tests/test_ts2_blood_nk_gate.py -q` completed with `6 passed`. Base HEAD was `983b37e`; no commit or push was made. Other TS2 organs, Smart-seq, skin, and paired-organ expression analyses were not attempted. Unlocking this arm requires a different authenticated blood object that satisfies the 2 GiB ceiling or an explicit policy decision from Reza; no access request or external contact was made.
