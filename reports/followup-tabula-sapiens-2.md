# Follow-up Step 2: Tabula Sapiens 2.0 metadata gate

## Status

The supplementary-table metadata gate used only the two released workbooks below. The abandoned 45 GB all-cells file was not downloaded or opened. The supplementary tables are marked `supplementary_tables_only_no_H5AD_verification`; supplementary cell counts are not raw expression counts. The subsequently approved blood NK H5AD arm is reported separately below.

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

The workbooks alone did not verify cell barcodes, duplicate cells across release objects, exact assay/chemistry, H5AD disease harmonization, anatomical position, enrichment or suspension type, raw-count layers, `is_primary_data`, per-cell QC/depth, or concordance between workbook and H5AD labels. Curated clinical columns in `media-24.xlsx` are not substitutes for a per-cell H5AD disease field. The blood H5AD checks below resolve the relevant fields only for the blood NK arm; TSP9/TSP20 ages and sex and other-organ H5AD fields remain unresolved.

No expression analysis was run at the supplementary-only stage. Other within-organ analyses remain pending suitable per-organ objects and their own gates. Smart-seq must remain separate from the primary 10x arm. Confirmatory paired-organ interaction analyses are blocked by the supplementary shared-donor/age screen. Skin-specific aging analyses are blocked by both the four-donor floor and the absence of an authenticated TS2 skin analysis. No plasma-proteomic organ-age claim is possible from these transcriptomic summaries.

## Execution and validation

The script is `scripts/prepare_ts2_metadata.py`. The final run used the existing `agescape` environment, a 512 MiB virtual-memory limit, and single-threaded numerical libraries:

```bash
ulimit -v 524288
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /usr/bin/time -v /home/reza/miniconda3/envs/agescape/bin/python \
  scripts/prepare_ts2_metadata.py
```

It exited 0 in 0.56 seconds with maximum resident set size 127,160 kB (about 124.2 MiB). An independent CSV audit rechecked row counts, key uniqueness, all five cell-total reconciliations, method-specific totals, the two missing donor identities, the TSP30 conflict, sex agreement, and the absence of passing organ-pair screens; all assertions passed. No large file was accessed, and no download was started.

## Blood 10x natural-killer-cell expression analysis

### Approved asset and authentication

The public [CELLxGENE collection API](https://api.cellxgene.cziscience.com/curation/v1/collections/e5f58829-1a66-40b5-a624-9046778e74f5), collection version `61d34250-a255-4c7c-b2f0-0e012fbe2e78`, contains one organ-specific blood dataset: dataset `983d5ec9-40e8-4512-9e65-a572a9c486cb`, version `1bf55482-ad49-446a-87c6-f3442bebcef1`. Its sole H5AD is `https://datasets.cellxgene.cziscience.com/1bf55482-ad49-446a-87c6-f3442bebcef1.h5ad`.

After the 3 GiB maximum was explicitly approved, the file was downloaded with `scripts/resource_guard.py download --min-free-gb 30 --max-file-gb 3` to `/home/reza/projects/agescape-revision-data/ts2/blood.h5ad`. The download completed in about 23 seconds. Authentication found:

- Exact size: 2,568,875,081 bytes (2.392451354 GiB).
- SHA-256: `667680ecdb017d3b4f3393ad4b668ff82345a9297aa38d7ab2b2830df0d1aa52`.
- Dimensions: 85,233 cells × 60,606 features.
- `raw.X`: CSR with 201,330,996 stored values; every value was checked in chunks and was finite, nonnegative, and integer-valued. Normalized `X`, `decontXcounts`, and `scale_data` were not used for expression outcomes.

The 45 GB all-cells object and old `skin.h5ad` were not used.

### H5AD metadata and eligibility gates

The object contains 2,985 10x natural-killer cells. Donor counts exactly reproduce the supplementary table for all nine donors: TSP1 71, TSP2 243, TSP7 1,641, TSP8 14, TSP10 285, TSP14 64, TSP21 250, TSP25 406, and TSP27 11. Age and sex agree with the curated supplementary values for every donor. All NK cells are blood, immune compartment, normal disease, and cell suspension. Assays are 10x 3' v3 and 10x 5' v2; Smart-seq was not pooled into the analysis.

The frozen seven-donor family is unchanged: TSP1, TSP2, TSP7, TSP10, TSP14, TSP21, and TSP25, with 2,960 cells. TSP8 and TSP27 remain excluded by the 20-cell floor. Raw NK depths range from 2,601 to 38,867 UMIs, and depth is computed from `raw.X`, not assumed from the release's `total_counts` column.

All 85,233 observation indices and `observation_joinid` values are unique within this single blood object. However, `is_primary_data` is false for every row, including all NK cells. The analysis therefore makes no `is_primary_data=True` or primary-only claim: it analyzes each observation once in this per-organ object and never concatenates overlapping release objects. This scope is recorded in the gate CSV and manifest rather than silently filtering to zero cells or relabelling the flag.

The full donor design was checked before outcomes and is rank 5 of 5 with two residual degrees of freedom: intercept, continuous age per 10 years, sex, any 5′ assay contribution from that donor, and donor mean log raw UMI. All models are exploratory donor-level OLS; no cell-level inference was used.

### Gene-wise donor-pseudobulk results

For each donor and gene, the outcome is the mean of cell-level log1p(CPM), exactly as frozen. All 32,830 variable genes with finite donor values and positive standard error were tested. Effects are changes in mean cell-level log1p(CPM) per 10 years of age with ordinary OLS 95% confidence intervals; Benjamini-Hochberg correction was applied once across all 32,830 genes.

**No gene had BH q < 0.05.** The smallest q-value was 0.8876988605. Complete effect sizes, standard errors, intervals, p-values, q-values, and model diagnostics are in `results/followup-ts2-blood-nk-gene-models.csv`. With only seven donors and two residual degrees of freedom, these are imprecise exploratory estimates, not evidence of transcriptome-wide stability.

### P1 transport and matched-null result

All 76 P1_Clean up genes and 42 down genes were present in `raw.var`, and every expression-decile random-signature pool exceeded its required size. The established 200-draw expression-matched null used seed 20261105 and the same donor design.

- Covariate-adjusted P1 age effect: **+0.0317244765 per 10 years**.
- Standard error: 0.0384871379; residual df = 2.
- Ordinary 95% CI: **−0.1338723124 to +0.1973212655**; ordinary OLS p = **0.4964344685**.
- Matched-null mean effect: −0.0026493494; null SD: 0.0127062044.
- Centered two-sided empirical p: **0.009950248756** (one of 200 null draws met or exceeded the centered observed magnitude; `(1 + 1)/201`). The single-endpoint BH value is numerically the same.

These tests answer different questions. The ordinary regression does not distinguish the P1 age slope from zero. The matched null asks whether the P1 slope is unusual relative to expression-matched random gene sets and is small. The latter is a **signature-specificity result, not an FDR-significant chronological-aging association**, and it does not override the wide direct-age confidence interval. Donor scores, the fitted model, gene coverage, pool floors, and all 200 null effects are committed in the aggregate CSVs.

### Off-program and Hellinger dispositions

The established human young-reference rule requires at least two eligible donors aged 20–40 years. Only TSP10 (age 33) qualifies. Although TSP10 contributes 285 NK cells, leave-one-donor-out reference construction leaves zero young reference cells for TSP10 itself. The two-young-donor and minimum-50-remaining-reference-cell requirements therefore fail.

- **Hellinger displacement: untestable.** No distances or age model were computed.
- **Off-program expression: untestable.** No cross-fitted reference markers or off-program outcomes were computed.

These are eligibility failures, not null results, and criteria were not relaxed after seeing expression outcomes.

### Outputs, execution, and limitations

- `scripts/run_followup_ts2_blood_nk.py`: streamed raw-count validation, donor pseudobulk construction, gene models, P1 null, gates, and manifest updates.
- `results/followup-ts2-blood-nk-h5ad-gates.csv`: ten explicit H5AD and endpoint gates.
- `results/followup-ts2-blood-nk-donor-summary.csv`: seven donor rows with covariates, depth, and P1 scores.
- `results/followup-ts2-blood-nk-gene-models.csv`: 32,830 tested genes.
- `results/followup-ts2-blood-nk-p1-model.csv`, `-p1-null.csv`, `-p1-gene-coverage.csv`, and `-p1-null-pools.csv`: P1 estimates and reproducibility inputs/outputs.
- Updated asset, donor-eligibility, endpoint-status, and manifest files record the 3 GiB approval, exact hash, H5AD checks, and endpoint dispositions.
- External caches are under `/home/reza/projects/agescape-revision-data/ts2/`; only committed-size aggregate CSVs are in `results/`.

The final guarded run used `resource_guard.py run --min-free-gb 30 --memory-gb 4`, backed HDF5 group reads, 128-cell CSR batches, and single-threaded numerical libraries. The final cache-reuse verification run exited 0 in 5.04 seconds with maximum RSS 361,640 kB (about 353 MiB). Independent checks reproduced the P1 OLS estimate/CI/null p-value, recomputed all gene BH values, refitted four genes from the safe `allow_pickle=False` cache, and verified 32,830 unique tested genes.

The analysis remains cross-sectional and exploratory: seven donors, two residual degrees of freedom, unequal donor cell counts, mixed 3′/5′ assays, and multiple enrichment/sample-preparation histories limit inference. No equivalence, longitudinal, causal, or independent TS1-to-TS2 replication claim is made. Other organs, Smart-seq, skin, and paired-organ expression analyses were not attempted. No Lake file was modified, and no commit or push was made.
