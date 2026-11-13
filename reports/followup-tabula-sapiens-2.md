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
