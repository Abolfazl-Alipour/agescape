# Dataset 08: Human Retinal Pigment Epithelium Bulk RNA-seq (Butler et al., 2021, GSE159435)

## Verdict
**In bulk human RPE (`n = 13` donors, ages 31 to 93 yr, `8M / 5F`), batch- and sex-adjusted donor OLS detected `2` significant age-associated genes out of `16212` tested at `q_BH < 0.05` (`RUBCNL` `+0.2560/dec`, `q_BH = 0.0128`; `ARL6IP5` `-0.0891/dec`, `q_BH = 0.0370`), while `P1_Clean` mouse signature transport failed in the opposite (negative) direction (`composite_effect = -0.0969`, `95% CI [-0.1880, -0.0058]`, analytic `p = 0.0398`, 200-draw empirical `p = 0.0050`).** Neither mouse UP genes (`up_effect = -0.0624`, `p = 0.1932`) nor DOWN genes (`down_effect = +0.0345`, `p = 0.2419`) recapitulate mouse direction.

## Assay Scope and Excluded Families
- **Bulk RNA-seq design**: Pooled macular plus peripheral retinal pigment epithelium (RPE) from `13` human donors with no reported retinal pathology (Butler et al., 2021, PMID 33934486, PMCID PMC8184696; GEO `GSE159435`).
- **No Family A fractions and no cross-cell-type concordance**: Because these are bulk tissue counts across heterogeneous macular/peripheral RPE (plus potential trace choroid/photoreceptor carryover) rather than dissociated single cells, Family A compositional shifts and cross-cell-type directional concordance cannot be computed without confounding within-cell expression and tissue composition. Only Family B (gene-level age OLS) and Family C (`P1_Clean` ortholog signature) are run.
- **Count provenance and mapping**: Raw integer HTSeq counts (`noninteger_entries = 0`, `max = 655797`) were extracted from `GSE159435_RAW.tar` across all 13 GSMs (`GSM4829762` to `GSM4829774`), removing the 5 `__*` HTSeq summary lines and mapping Ensembl IDs to HGNC symbols via the GRCh38/ortholog table (`117 / 118` `P1_Clean` genes mapped; only `CD24` absent).

## Cohort Structure (`n = 13` donors, 2 Batches)
- Notice the batch-age imbalance in GEO metadata ([donors CSV](../results/followup-organ-08-eye-rpe-gse159435-donors.csv)): `batch s1` (`n = 7`) spans the full adult range `31 to 93 yr` (`31M, 40F, 51M, 67F, 77M, 86M, 93F`), whereas `batch s2` (`n = 6`) contains exclusively older donors aged `61 to 92 yr` (`61F, 62M, 71M, 73M, 73F, 92M`). Adjusting for `batch_s2` is essential and `Sens_BatchS1_Only` tests the homogeneous single-batch `31 to 93 yr` ladder.

| Donor (GEO GSM) | Age (yr) | Sex | Batch | Library Size (mapped counts) |
|---|---:|:---:|:---:|---:|
| donor1 (GSM4829762) | 31 | M | s1 | 9530400 |
| donor2 (GSM4829763) | 40 | F | s1 | 12830995 |
| donor3 (GSM4829764) | 51 | M | s1 | 12705594 |
| donor4 (GSM4829765) | 67 | F | s1 | 15024150 |
| donor5 (GSM4829766) | 77 | M | s1 | 11644779 |
| donor6 (GSM4829767) | 86 | M | s1 | 12696924 |
| donor7 (GSM4829768) | 93 | F | s1 | 13057453 |
| donor8 (GSM4829769) | 61 | F | s2 | 16750222 |
| donor9 (GSM4829770) | 62 | M | s2 | 9323324 |
| donor10 (GSM4829771) | 71 | M | s2 | 14303930 |
| donor11 (GSM4829772) | 73 | M | s2 | 14193664 |
| donor12 (GSM4829773) | 73 | F | s2 | 7414353 |
| donor13 (GSM4829774) | 92 | M | s2 | 14871766 |

## Family B: Gene-Level Donor Age OLS
- **Primary model (`Primary_BatchSexAdj`, `~ age_per_decade + sex_male + batch_s2 + log_lib`, `M_B = 16212`)**: **`2` genes at `q_BH < 0.05`** ([gene summary CSV](../results/followup-organ-08-eye-rpe-gse159435-gene-summary.csv), [top genes CSV](../results/followup-organ-08-eye-rpe-gse159435-top-genes.csv)).
- **Ensembl ID annotation (`top-genes.csv`)**: All `44` `ENSG*` rows (`28` unique `ENSG*` IDs across the 300-row top-100-per-suite file) are annotated via Ensembl GRCh38 REST API (`/lookup/id` + `/archive/id`) and Ensembl GRCh37 REST API (`grch37.rest.ensembl.org/lookup/id`): `9` unique active GRCh38 symbols (`11` rows), `3` retired GRCh38 IDs with active GRCh38 replacements (`6` rows: `ENSG00000255028` -> `LINC02715`, `ENSG00000237461` -> `NAMA`, `ENSG00000259126` -> `LINC02277`), `4` retired GRCh38 IDs resolved via GRCh37 (`7` rows including `ENSG00000196963` = `PCDHB16` and `ENSG00000132142` = `ACACA`), and `12` active GRCh38 novel lncRNAs/pseudogenes lacking HGNC symbols (`20` rows, headed by `ENSG00000235912`).
- **`ENSG00000235912` annotation (rank 10)**: Active GRCh38 `transcribed_processed_pseudogene` ("novel pseudogene", Ensembl GRCh38 REST `/lookup/id`), with no HGNC symbol (GRCh37 clone locus `RP1-159A19.3`, Ensembl GRCh37 REST `/lookup/id`); `q_BH = 0.1048`.

| Suite | Group | n_donors | Genes Tested | Global q<0.05 | Min p-value | Min q_BH |
|---|---|---:|---:|---:|---:|---:|
| Primary_BatchSexAdj | RPE_Bulk | 13 | 16212 | 2 | 7.88e-07 | 0.0128 |
| Sens_SexAdj_NoBatch | RPE_Bulk | 13 | 16212 | 3 | 8.54e-07 | 0.0138 |
| Sens_BatchS1_Only | RPE_Bulk_s1 | 7 | 16416 | 0 | 1.48e-05 | 0.1822 |

### Top 10 Ranked Genes (`Primary_BatchSexAdj`)
| Gene | Ensembl / Mapping Annotation (Source) | Mean log2CPM | Effect (log2CPM/dec) | SE | p-value | q_BH |
|---|---|---:|---:|---:|---:|---:|
| RUBCNL | input_symbol (GEO raw symbol) | 2.00 | +0.2560 | 0.0187 | 7.88e-07 | 0.0128 |
| ARL6IP5 | input_symbol (GEO raw symbol) | 7.42 | -0.0891 | 0.0082 | 4.56e-06 | 0.0370 |
| CACNA2D2 | input_symbol (GEO raw symbol) | 4.32 | -0.3002 | 0.0320 | 1.38e-05 | 0.0737 |
| MASTL | input_symbol (GEO raw symbol) | 4.31 | -0.1234 | 0.0138 | 1.95e-05 | 0.0737 |
| LINC03040 | input_symbol (GEO raw symbol) | 2.67 | +0.5035 | 0.0599 | 3.05e-05 | 0.0737 |
| F10 | input_symbol (GEO raw symbol) | 2.50 | -0.3515 | 0.0418 | 3.07e-05 | 0.0737 |
| PAX8 | input_symbol (GEO raw symbol) | 1.37 | -0.2609 | 0.0313 | 3.21e-05 | 0.0737 |
| LARGE2 | input_symbol (GEO raw symbol) | 3.41 | +0.4275 | 0.0521 | 3.64e-05 | 0.0737 |
| ANKRD54 | input_symbol (GEO raw symbol) | 3.15 | +0.1472 | 0.0185 | 4.48e-05 | 0.0808 |
| ENSG00000235912 | active_grch38_no_hgnc_symbol (`transcribed_processed_pseudogene`; GRCh37 clone `RP1-159A19.3`; Ensembl GRCh38/37 REST) | 1.97 | +0.2396 | 0.0318 | 6.78e-05 | 0.1048 |

## Family C: P1_Clean Mouse Signature Transport
- Across all three models (`Primary_BatchSexAdj`, `Sens_SexAdj_NoBatch`, `Sens_BatchS1_Only`; `66` up and `36` down `P1_Clean` genes expressed above `1 CPM` in `>=50%` of donors), **there is zero positive transport of `P1_Clean`**: all three composite point estimates are **negative** (`-0.0704` to `-0.0969` per decade), driven by inverted average slopes relative to mouse aging ([P1 signature CSV](../results/followup-organ-08-eye-rpe-gse159435-p1-signature.csv)). This reinforces the conclusion that `P1_Clean` does not positively generalize into aging human RPE.

| Suite | n_donors | P1 Up/Down | Composite Effect | Analytic p | Empirical p (200 draws) | q_BH (Comp) | q_BH (Emp) |
|---|---:|:---:|---:|---:|---:|---:|---:|
| Primary_BatchSexAdj | 13 | 66/36 | -0.0969 | 0.0398 | 0.0050 | 0.0398 | 0.0050 |
| Sens_SexAdj_NoBatch | 13 | 66/36 | -0.0930 | 0.0313 | 0.0050 | 0.0313 | 0.0050 |
| Sens_BatchS1_Only | 7 | 66/36 | -0.0704 | 0.1020 | 0.0050 | 0.1020 | 0.0050 |

## Key Limitations
1. **Bulk tissue confounding**: Macular and peripheral RPE were pooled together, and bulk RNA cannot distinguish intrinsic RPE aging from region ratio or trace adjacent tissue variation.
2. **Batch imbalance**: Batch `s2` contains only donors aged `61 to 92 yr` (zero young donors), reducing effective cross-lifespan information primarily to the `7` `batch s1` donors (`31 to 93 yr`) once `batch_s2` is controlled (`Sens_SexAdj_NoBatch` hits largely reflect batch `s2` vs `s1` effects correlated with age).
3. **Sample size**: With `n = 13` (`residual_df = 8`), power for moderate gene-level slopes is limited.

## Output Files
- [Donors](../results/followup-organ-08-eye-rpe-gse159435-donors.csv)
- [Gene Summary](../results/followup-organ-08-eye-rpe-gse159435-gene-summary.csv)
- [Top Genes](../results/followup-organ-08-eye-rpe-gse159435-top-genes.csv)
- [P1 Signature](../results/followup-organ-08-eye-rpe-gse159435-p1-signature.csv)
