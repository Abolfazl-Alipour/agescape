# Dataset 05: Human Bone Marrow Sorted HSCe Bulk RNA-seq (Adelman et al., GSE104406)

## 1. Study Verification and Data Provenance
- **Accession**: `GSE104406` (`GSM2797287` to `GSM2797306`), Primary source: Adelman et al., *Cancer Discovery* 2019 (`PMID: 31085557`).
- **Assay and Cell Population**: Bulk RNA-seq on FACS-sorted human bone marrow hematopoietic stem-cell-enriched populations (`Lin- CD34+ CD38-`), `n = 20` healthy donors without known hematological malignancy.
- **Age and Sex Breakdown**:
  - Young (`n = 10`, ages `18 to 28` plus `Young_25_RNA` labeled `<30`): `7M / 3F`.
  - Aged (`n = 10`, ages `65 to 75`): `1M / 9F`.
  - **Confounding note**: Strong age-by-sex imbalance (`Young` 70% male vs `Aged` 10% male), requiring explicit `sex_male` adjustment alongside `log_lib`.
- **Feature ID normalization**: Raw tar count ID style `ensembl`; `57,255` rows -> `57,247` unique symbols -> `57,100` non-technical genes (`147` removed).

- **Total cohort (`n = 20`)**: See [donors CSV](../results/followup-organ-05-bm-adelman-gse104406-donors.csv).

| GSM | Title | Age Raw | Group | Sex | Non-Tech Lib Size |
|---|---|---|---|---|---:|
| GSM2797287 | Young_1_RNA | 18 | young | M | 18,028,144 |
| GSM2797288 | Young_12_RNA | 25 | young | M | 18,467,611 |
| GSM2797289 | Young_16_RNA | 27 | young | F | 15,117,098 |
| GSM2797290 | Young_17_RNA | 28 | young | M | 35,279,564 |
| GSM2797291 | Young_25_RNA | <30 | young | M | 17,132,179 |
| GSM2797292 | Young_4_RNA | 21 | young | F | 25,710,609 |
| GSM2797293 | Young_5_RNA | 22 | young | M | 27,742,026 |
| GSM2797294 | Young_6_RNA | 22 | young | M | 17,724,927 |
| GSM2797295 | Young_7_RNA | 22 | young | F | 35,799,527 |
| GSM2797296 | Young_9_RNA | 22 | young | M | 49,799,166 |
| GSM2797297 | Old_11_RNA | 67 | old | F | 21,403,248 |
| GSM2797298 | Old_18_RNA | 70 | old | F | 36,814,349 |
| GSM2797299 | Old_19_RNA | 70 | old | F | 23,701,865 |
| GSM2797300 | Old_21_RNA | 71 | old | M | 14,183,782 |
| GSM2797301 | Old_22_RNA | 71 | old | F | 19,767,246 |
| GSM2797302 | Old_23_RNA | 72 | old | F | 22,669,286 |
| GSM2797303 | Old_25_RNA | 75 | old | F | 29,487,076 |
| GSM2797304 | Old_28_RNA | 74 | old | F | 20,495,188 |
| GSM2797305 | Old_20_RNA | 75 | old | F | 16,939,085 |
| GSM2797306 | Old_4_RNA | 65 | old | F | 30,356,603 |

## 2. Scope Rules: FACS-Sorted Bulk HSCe Only
- **Family A (CLR cell-type fractions)**: **NOT RUN**. Cells were FACS-sorted (`Lin- CD34+ CD38-`) and sequenced in bulk, so proportions of bone marrow cell types cannot be inferred.
- **Cross-cell-type concordance**: **NOT RUN**. Single sorted population in bulk RNA-seq.
- **Family B (Donor-level gene age OLS)**: Adjusting for `sex_male + log_lib` (`CPM >= 1.0` in `>= 50%` of donors):
  - **Primary (`n = 19` exact ages)**: `M_B = 14,115` genes (`residual_df = 15`).
  - **Sensitivity 1 (`n = 20`, `Young_25_RNA` coded `25`)**: `M_B = 14,296` genes.
  - **Sensitivity 2 (`n = 20`, binary `old_vs_young`)**: `M_B = 14,296` genes.
- **Family C (`P1_Clean` signature)**: `K_C = 1`, 200 expression-decile-matched random draws + analytic OLS.

## 3. Results

### Family B: Bulk Sorted HSCe Gene Age OLS
- **Primary (`n = 19` exact ages)**: **`72`** genes passed `q_BH < 0.05` (minimum `p = 9.55e-07`, minimum `q_BH = 0.0083`; `#1` hit is `ENSG00000119917` = **`IFIT3`**, `+0.5212 log2CPM/decade`, `p = 9.55e-07`, `q_BH = 0.0083`) ([gene summary CSV](../results/followup-organ-05-bm-adelman-gse104406-gene-summary.csv), [top genes CSV](../results/followup-organ-05-bm-adelman-gse104406-top-genes.csv)).
- **Ensembl ID annotation (`top-genes.csv`)**: All `41` `ENSG*` rows in the top-200 file (`14` among the `72` `q_BH < 0.05` hits) are annotated via Ensembl GRCh38 REST API (`/lookup/id` + `/archive/id`, `34` active HGNC symbols, e.g., `IFIT3`, `ZNF230`, `ZNF616`, `ZNF200`, `NUDT16L2P`, `GIMAP2`, `LINC02175`, `KRT8P12`, `RPL23AP7`, `SNHG1`, `RPL34P20`) and Ensembl GRCh37 REST API (`grch37.rest.ensembl.org/lookup/id`, `6` retired GRCh38 IDs including `ENSG00000203819` = `HIST2H2BC` and `ENSG00000233895` = `RP1-122P22.2`), plus `1` active GRCh38 novel lncRNA (`ENSG00000266965`, no HGNC symbol; GRCh37 clone `RP11-712P20.2`).

| Gene (Raw) | Ensembl Symbol | Mapping Status (Source) | Mean log2CPM | Effect (/10 yr) | SE | p-value | q_BH |
|---|---|---|---:|---:|---:|---:|---:|
| ENSG00000119917 | IFIT3 | active_grch38 (Ensembl GRCh38 REST `/lookup/id`) | 4.38 | +0.5212 | 0.0657 | 9.55e-07 | 0.0083 |
| HSPA2 | HSPA2 | input_symbol (GEO raw symbol) | 1.77 | +0.4668 | 0.0599 | 1.18e-06 | 0.0083 |
| TRIM68 | TRIM68 | input_symbol (GEO raw symbol) | 4.96 | +0.2553 | 0.0366 | 4.46e-06 | 0.0149 |
| LMNA | LMNA | input_symbol (GEO raw symbol) | 6.32 | -0.7818 | 0.1135 | 5.19e-06 | 0.0149 |
| AKAP8 | AKAP8 | input_symbol (GEO raw symbol) | 6.46 | -0.1450 | 0.0215 | 6.64e-06 | 0.0149 |
| TEKT2 | TEKT2 | input_symbol (GEO raw symbol) | 1.44 | -0.5129 | 0.0771 | 7.68e-06 | 0.0149 |
| ENSG00000159882 | ZNF230 | active_grch38 (Ensembl GRCh38 REST `/lookup/id`) | 3.25 | +0.6458 | 0.0978 | 8.36e-06 | 0.0149 |
| ESAM | ESAM | input_symbol (GEO raw symbol) | 3.49 | +0.2733 | 0.0414 | 8.47e-06 | 0.0149 |
| TNPO3 | TNPO3 | input_symbol (GEO raw symbol) | 6.15 | +0.1540 | 0.0244 | 1.39e-05 | 0.0219 |
| PPP1R10 | PPP1R10 | input_symbol (GEO raw symbol) | 7.75 | +0.1986 | 0.0327 | 2.11e-05 | 0.0263 |

- **Sensitivity (`n = 20` continuous, `Young_25_RNA = 25`)**: `134` genes passed `q_BH < 0.05` (min `q_BH = 0.0079`).
- **Sensitivity (`n = 20` binary `old_vs_young`)**: `218` genes passed `q_BH < 0.05` (min `q_BH = 0.0033`).

### Family C: Mouse `P1_Clean` Signature Check
([P1 signature CSV](../results/followup-organ-05-bm-adelman-gse104406-p1-signature.csv))
- **Primary (`n = 19` exact ages)**: `48` UP / `15` DOWN detected; UP slope `-0.0070` (`p = 0.7999`), DOWN slope `-0.0217` (`p = 0.5131`), Composite `+0.0147` (`analytic p = 0.7407`, 200-draw random null mean `-0.0073` +/- `0.0370`, **`empirical_null_pvalue = 0.5373`**, `q_BH = 0.5373`).
- **Sensitivity (`n = 20` continuous)**: Composite `+0.0078` (`analytic p = 0.8554`, **`empirical_null_pvalue = 0.6517`**).
- **Sensitivity (`n = 20` binary)**: Composite `+0.0094` (`analytic p = 0.9635`, **`empirical_null_pvalue = 0.7562`**).

## 4. Limitations and Verdict
- **Limits**: FACS-sorted (`Lin- CD34+ CD38-`) bulk RNA-seq cannot deconvolve HSC subpopulations or evaluate tissue composition, and sex is strongly imbalanced across age (`7M/3F` young vs `1M/9F` aged).
- **Verdict**: Null transport for the mouse `P1_Clean` aging signature in sorted human bone marrow HSCe (`empirical p = 0.5373` primary exact-age `n = 19`; `0.7562` binary `n = 20`).

## 5. Output Artifacts
- [Donors](../results/followup-organ-05-bm-adelman-gse104406-donors.csv)
- [Gene Summary](../results/followup-organ-05-bm-adelman-gse104406-gene-summary.csv)
- [Top Genes](../results/followup-organ-05-bm-adelman-gse104406-top-genes.csv)
- [P1 Signature](../results/followup-organ-05-bm-adelman-gse104406-p1-signature.csv)
