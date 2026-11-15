# Dataset 04: Human Pancreas Aging Smart-seq2 (Enge et al., GSE81547)

## 1. Study Verification and Data Provenance
- **Accession**: `GSE81547` (`2,544` single cells across `8` donors), Primary source: Enge et al., *Cell* 2017 (`PMID: 28965763`, `DOI: 10.1016/j.cell.2017.09.004`).
- **Critical Span and Protocol Limits**:
  1. **Span is development/pediatric to midlife ONLY**: Ages span `1 to 54` years (`GEO donor_age: 1, 5, 6, 21, 22, 38, 44, 54`; `1_M` is coded `1` in GEO `series_matrix` and reported as 1 month in Enge et al., `6M / 2F`). There are **zero old donors (`>= 60` yr)**.
  2. **FACS-gated Smart-seq2 capture**: Single cells were FACS-sorted into 96-well plates for Smart-seq2 (`2,544` cells parsed; `5` HTSeq summary rows and `92` ERCC spike-in rows excluded before CPM normalization).

- **Total cohort (`n = 8`)**: See [donors CSV](../results/followup-organ-04-pancreas-enge-gse81547-donors.csv).

| Donor ID | Age Raw (Years) | Sex | QC Cells |
|---|---:|---|---:|
| 1_M | 1 (1.000) | M | 221 |
| 5_M | 5 (5.000) | M | 331 |
| 6_M | 6 (6.000) | M | 178 |
| 21_M | 21 (21.000) | M | 505 |
| 22_M | 22 (22.000) | M | 286 |
| 38_F | 38 (38.000) | F | 474 |
| 44_F | 44 (44.000) | F | 277 |
| 54_M | 54 (54.000) | M | 272 |

## 2. Explicit Multiple-Testing Families
- **Inferential unit**: Donor (`n = 8` all donors; `n = 5` adults `>= 21` yr sensitivity). Continuous `age_per_decade = age_years / 10`.
- **Family A (CLR cell-type proportions)**: `7` FACS-captured pancreas cell types present in `>= 4` donors (`total >= 20` cells), adjusting for `sex_male`, BH across `K_A = 7`. Because cells were FACS-sorted rather than unbiased whole-tissue dissociation, proportions are library-capture fractions rather than true tissue abundances.
- **Family B (Donor-level gene age test + concordance)**: Donor pseudobulks (`>= 20` cells in `>= 5` donors; `CPM >= 1.0` in `>= 50%` donors) adjusting for `sex_male + log_lib` (`all_8`) or `log_lib` (`adult_5`). Global BH across `M_B = 46,655` `(cell_type, gene)` tests. Concordance measured across shared expressed genes with `1,000` donor-age permutations.
- **Family C (`P1_Clean` mouse aging signature)**: `76` UP / `42` DOWN mapped human genes, `K_C = 4` cell types (`all_8`) and `K_C = 1` (`adult_5`), 200 expression-decile-matched random draws + BH.

## 3. Results

### Family A: FACS Library CLR Proportions (All 8 Donors and Adult-Only Sensitivity)
In the full `n = 8` cohort (`1 to 54` yr, `K_A = 7`, `covariates = sex_male`, `residual_df = 5`), `0` cell types passed `q_BH < 0.05` ([fractions CSV](../results/followup-organ-04-pancreas-enge-gse81547-fractions.csv)).

| Cell Type | Total Cells | Median Frac | CLR Slope (/10 yr) | p-value | q_BH (K=7) |
|---|---:|---:|---:|---:|---:|
| alpha | 998 | 0.4320 | +0.267 | 0.0707 | 0.4950 |
| delta | 83 | 0.0244 | -0.197 | 0.2578 | 0.9022 |
| acinar | 411 | 0.1534 | -0.120 | 0.5285 | 0.9790 |
| unsure | 262 | 0.0510 | +0.131 | 0.7108 | 0.9790 |
| beta | 348 | 0.1390 | -0.039 | 0.8919 | 0.9790 |
| ductal | 389 | 0.1788 | -0.033 | 0.9314 | 0.9790 |
| mesenchymal | 53 | 0.0091 | -0.010 | 0.9790 | 0.9790 |

- **Adult-only sensitivity (`adult_5`, `n = 5`, ages `21 to 54` yr, `K_A = 6`, `covariates = sex_male`, `residual_df = 2`)** ([fractions CSV](../results/followup-organ-04-pancreas-enge-gse81547-fractions.csv)): `alpha` (`828` cells, median fraction `0.4768`) reached nominal/FDR significance (`CLR slope +0.150` per decade, `SE = 0.0082`, `p = 0.0030`, `q_BH = 0.0178`, `R^2 = 0.995`), while the remaining `5` adult cell types (`delta` `-0.411`, `p = 0.0919`, `q = 0.2758`; `unsure` `-0.301`, `p = 0.6574`, `q = 0.8201`; `beta` `+0.277`, `p = 0.6850`, `q = 0.8201`; `acinar` `+0.128`, `p = 0.7722`, `q = 0.8201`; `ductal` `-0.129`, `p = 0.8201`, `q = 0.8201`) were null.
- **Explicit limits on `adult_5` alpha fraction result**: Only `n = 5` donors (`residual_df = 2`), fragile, FACS-sorted Smart-seq2 library-capture proportions rather than unbiased whole-islet tissue composition, model adjusts `sex_male` only (all `8` non-diabetic donors in `GSE81547`, with no T2D covariate in the model), and there are zero elderly donors (`>= 60` yr).

| Cohort | Cell Type | Donors (n) | Total Cells | Median Frac | Covariates | Residual df | CLR Slope (/10 yr) | SE | p-value | q_BH (K=6) | R^2 |
|---|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|
| adult_5 | alpha | 5 | 828 | 0.4768 | sex_male (no T2D covariate) | 2 | +0.150 | 0.0082 | 0.0030 | 0.0178 | 0.995 |
| adult_5 | delta | 5 | 49 | 0.0277 | sex_male | 2 | -0.411 | 0.1340 | 0.0919 | 0.2758 | 0.831 |
| adult_5 | unsure | 5 | 219 | 0.0478 | sex_male | 2 | -0.301 | 0.5829 | 0.6574 | 0.8201 | 0.589 |
| adult_5 | beta | 5 | 222 | 0.1324 | sex_male | 2 | +0.277 | 0.5905 | 0.6850 | 0.8201 | 0.150 |
| adult_5 | acinar | 5 | 214 | 0.1324 | sex_male | 2 | +0.128 | 0.3860 | 0.7722 | 0.8201 | 0.074 |
| adult_5 | ductal | 5 | 265 | 0.1762 | sex_male | 2 | -0.129 | 0.4984 | 0.8201 | 0.8201 | 0.454 |

### Family B: Donor-Level Gene Age OLS and Cross-Cell-Type Concordance
- **All 8 donors (`1 to 54` yr, `M_B = 46,655`)** ([gene summary CSV](../results/followup-organ-04-pancreas-enge-gse81547-gene-summary.csv), [top genes CSV](../results/followup-organ-04-pancreas-enge-gse81547-top-genes.csv), [concordance CSV](../results/followup-organ-04-pancreas-enge-gse81547-concordance.csv)): `0` genes passed `q_BH < 0.05` (minimum `q_BH = 0.7551`). Concordance across `4` cell types (`9178` common genes): mean pairwise Pearson `r = +0.0166` (`perm_p_one_sided = 0.9880`).
- **Adult-only sensitivity (`n = 5`, `21 to 54` yr)**: Only 1 cell type (`alpha`, `11,831` genes, `0` genes `q_BH < 0.05`) met `>= 20` cells across all 5 adult donors, so cross-cell-type concordance is not applicable.

| Cell Type | Donors | Genes Tested | Global q<0.05 | Min p-value | Min Global q_BH |
|---|---:|---:|---:|---:|---:|
| acinar | 6 | 11050 | 0 | 2.61e-04 | 0.9673 |
| alpha | 8 | 12509 | 0 | 3.56e-05 | 0.7551 |
| beta | 7 | 11322 | 0 | 2.98e-05 | 0.7551 |
| ductal | 5 | 11774 | 0 | 1.25e-04 | 0.9673 |

### Family C: Mouse `P1_Clean` Signature Check (All 8 Donors)
Minimum empirical null `q_BH = 0.5174` (minimum analytic `q_BH = 0.7431`) across `K_C = 4` cell types; adult-only sensitivity (`n = 5`, `alpha` `K_C = 1`) minimum empirical `q_BH = 0.0846` ([P1 signature CSV](../results/followup-organ-04-pancreas-enge-gse81547-p1-signature.csv)).

| Cell Type | Donors | UP/DOWN Detected | Composite Slope | Analytic p | Empirical Null p (200) | Empirical q_BH |
|---|---:|---|---:|---:|---:|---:|
| beta | 7 | 48/19 | +0.1162 | 0.5457 | 0.1940 | 0.5174 |
| ductal | 5 | 59/25 | -0.3336 | 0.1873 | 0.2736 | 0.5174 |
| alpha | 8 | 54/21 | +0.0539 | 0.5574 | 0.3881 | 0.5174 |
| acinar | 6 | 60/22 | +0.0069 | 0.9008 | 0.8358 | 0.8358 |

## 4. Verdict
- **No mouse aging transport**: `P1_Clean` shows null aging association across pancreas cell types in both the `1 to 54` yr cohort (`min q_BH = 0.5174`) and adult-only `21 to 54` yr `alpha` subset (`q_BH = 0.0846`).
- **Adult-only composition caveated**: In `adult_5` (`n = 5`, `residual_df = 2`), FACS `alpha` proportion rises nominally/FDR (`+0.150/dec`, `p = 0.0030`, `q_BH = 0.0178`, `R^2 = 0.995`), but this rests on `5` non-diabetic donors (`df = 2`, `sex_male` adjusted only, zero T2D covariate, zero elderly donors `>= 60` yr) under FACS-gated Smart-seq2 capture and is fragile ([fractions CSV](../results/followup-organ-04-pancreas-enge-gse81547-fractions.csv)).
- **Explicit limits**: Cohort spans pediatric development to midlife (`1 to 54` yr) with no elderly donors (`>= 60` yr) and uses FACS-sorted Smart-seq2 capture.

## 5. Output Artifacts
- [Donors](../results/followup-organ-04-pancreas-enge-gse81547-donors.csv)
- [Fractions](../results/followup-organ-04-pancreas-enge-gse81547-fractions.csv)
- [Gene Summary](../results/followup-organ-04-pancreas-enge-gse81547-gene-summary.csv)
- [Top Genes](../results/followup-organ-04-pancreas-enge-gse81547-top-genes.csv)
- [Concordance](../results/followup-organ-04-pancreas-enge-gse81547-concordance.csv)
- [P1 Signature](../results/followup-organ-04-pancreas-enge-gse81547-p1-signature.csv)
