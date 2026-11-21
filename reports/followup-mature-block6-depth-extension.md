# Block 6 mature-to-old (30-plus): Depth versus age extension

## Status and scope

Complete. This is Phase 2 Item 13 of the Agescape mature-to-old rebuild. The analysis specification, MDE80 definition, and verdict rules are reused verbatim from the frozen original Block 6 report (`reports/followup-11pt-block6-depth-extension.md`), changing only the donor age restriction to human mature donors aged 30 years and older (`age_years >= 30`). Because Block 6 models library sizes rather than gene scores, the mature P1m signature does not enter this block.

## Inputs and provenance

- Gao donor table: `results/followup-between-subject-gao-donors.csv` (non-pooled cell-type rows restricted to `age_years >= 30`).
- Yazar donor table: `results/followup-between-subject-yazar-donors.csv` (non-pooled cell-type rows restricted to `age_years >= 30`).
- Original Block 6 tables (read-only comparison and Part B verification source): `results/followup-11pt-block6-depth-per-cell-type.csv` and `results/followup-11pt-block6-mssm-depth.csv`.
- MSSM h5ad verification source: `/home/reza/projects/agescape-revision-data/mssm.h5ad` (column-selective `obs` read via `h5py`; sparse expression matrix `X` not loaded).
- Pooled mature cross-check only (excluded from Part A): `results/followup-mature-nov17-depth.csv`. At 30-plus, pooled adjusted `mean_log_umi_ge1500` slopes are null in both Gao (`+0.0056` per decade, `p = 0.6041`, N = 109, lost 8) and Yazar (`+0.0004` per decade, `p = 0.7423`, N = 932, lost 49). Pooled adjusted `mean_log_umi_all` is borderline negative in Gao (`-0.0406`, `p = 0.0541`) and null in Yazar (`+0.0026`, `p = 0.1722`), while Yazar captured cell counts still decline (`-0.0292` log cells per decade, `p = 2.08e-06`).
- Block 9 and Item 7 interneuron reconciliation: `results/followup-11pt-block9-mature-interneuron.csv`, `reports/followup-11pt-block9-mature-interneuron.md`, `results/followup-mature-interneuron*.csv`, and `reports/followup-mature-interneuron.md` were checked. They contain per-gene differential expression and cross-cohort replication outputs (using donor `mean_log_umi` as a regression covariate) and contain no newer donor-level depth-versus-age regressions. Their reported MSSM mature donor counts (VIP 193, pvalb 185, sst 174, GABA 153; ages 31 to 89) match Part B identically.

## Frozen analysis specification applied at 30-plus

### Part A: per-cell-type depth slopes in Gao and Yazar (age 30-plus)

For each non-pooled cell type with at least 6 eligible donors aged 30 years and older, donor-level OLS regresses `mean_log_umi_ge1500` on `age_per_decade + sex_male + batch`. Gao batch is the observed study fixed effect (up to seven studies, represented by an intercept and up to six study indicators). Yazar batch is pool fixed effects (`study_or_pool`). A donor is eligible for the primary outcome with at least 20 cells in `n_cells_ge1500`, `age_years >= 30`, and nonmissing outcome, age, sex, and batch. Slopes are in natural-log UMI per decade.

Sensitivity 1 (`sensitivity_adjusted`) uses the same adjusted model with `mean_log_umi_all` and at least 20 cells in `n_cells_all` per donor at `age_years >= 30`. Sensitivity 2 (`unadjusted` and `sensitivity_unadjusted`) fits an unadjusted age-only model for each outcome on its outcome-specific eligible 30-plus donors. No pooled rows, cell-level fits, downsampling, or rarefaction are used.

BH multiple-testing families are defined separately within each dataset over all non-pooled cell types meeting the 6-donor eligibility floor at 30-plus: 13 Gao and 17 Yazar candidates for `mean_log_umi_ge1500` (down from 18 Yazar originally because `CD16-negative, CD56-bright natural killer cell, human` dropped from 6 to 5 donors), and 13 Gao and 19 Yazar candidates for `mean_log_umi_all`. Within each family, BH is calculated over successfully fitted candidates (`status == FIT`). Unadjusted age-only fits receive no q-values. Collinear covariates are dropped without dropping age; any model with zero residual degrees of freedom is labeled `UNTESTABLE` (`insufficient_residual_df`), never null.

### Part B: formal MSSM mature neuron test (age 30-plus)

The four pre-stated MSSM neuronal units are `GABAergic neuron`, `VIP GABAergic cortical interneuron`, `pvalb GABAergic cortical interneuron`, and `sst GABAergic cortical interneuron`. Per-donor mean natural-log UMI uses full-depth per-cell `n_counts` passing the `fetch_unit_h5ad` QC route: primary data, disease `normal`, `n_counts >= 1500`, parsed age 30 to 120 years, and at least 20 eligible cells per donor per unit (minimum 6 donors per unit). Primary donor-level OLS regresses `mean_log_umi_full_depth` on `age_per_decade + sex_male`. The declared mature primary BH family consists of exactly these four MSSM neuronal units. Sensitivity 1 (`unadjusted`) fits age-only OLS on the same 30-plus eligible donors and receives no q-values.

Direct refitting from `/home/reza/projects/agescape-revision-data/mssm.h5ad` verified all four existing `mature_age30_plus` rows in `results/followup-11pt-block6-mssm-depth.csv` to machine precision (exact integer matches on `n_donors`, `n_cells`, and `residual_df`; maximum absolute float difference `9.714e-17` across slope, SE, 95 percent CI, p-value, MDE80, outcome SD, and age summaries). Because all four rows verified, the primary Part B table reuses those verified rows verbatim and applies the declared four-unit BH family.

### MDE80 and frozen verdict rules

MDE80 is `(t_0.975 + t_0.80) * SE` on the residual degrees of freedom, reported in natural-log UMI per decade and in outcome SD per decade. Under the frozen power threshold, MDE80 at most 0.50 SD per decade is required for a null to count as adequately powered (at most 0.25 SD per decade is high power).

- Part A per dataset: `PROVES` if at least one batch-adjusted primary cell-type slope is negative with BH q below 0.05. Otherwise `DOES NOT PROVE` requires all fitted primary slopes to be null or nonnegative with MDE80 at most 0.50 SD per decade and no `UNTESTABLE` primary candidate. Otherwise `INCONCLUSIVE`.
- Part B (MSSM mature four-unit family): `PROVES` if at least one of the four mature primary neuronal slopes is negative with BH q below 0.05. Otherwise `DOES NOT PROVE` requires all four mature primary fits to be null or nonnegative with MDE80 at most 0.50 SD per decade. Otherwise `INCONCLUSIVE`.
- Fixed seed: `20261120`. No random sampling, permutation, or bootstrap is used; all fits are deterministic.

## Results

All slopes and MDE80 values are in natural-log UMI per decade. `FIT` and `UNTESTABLE` describe model evaluability, not scientific verdicts.

### Part A primary fits (age 30-plus, `mean_log_umi_ge1500`)

| Dataset | Cell type | N orig | N | N lost | Slope | 95% CI | p | BH q | MDE80 | MDE80 SD | Status |
|---|---|---|---|---|---|---|---|---|---|---|---|
| gao | Adipocyte | 37 | 35 | 2 | -0.0657 | -0.1473 to 0.0159 | 0.1103 | 0.4507 | 0.1156 | 0.2971 | FIT |
| gao | Cardiomyocyte | 109 | 102 | 7 | -0.0243 | -0.0776 to 0.0290 | 0.3678 | 0.5410 | 0.0760 | 0.1550 | FIT |
| gao | Endocardial | 72 | 66 | 6 | -0.0356 | -0.0802 to 0.0090 | 0.1151 | 0.4507 | 0.0635 | 0.2268 | FIT |
| gao | Endothelial | 87 | 80 | 7 | 0.0124 | -0.0080 to 0.0328 | 0.2286 | 0.4507 | 0.0291 | 0.1369 | FIT |
| gao | Epicardial | 20 | 18 | 2 | 0.0380 | -0.0512 to 0.1272 | 0.3745 | 0.5410 | 0.1251 | 0.5014 | FIT |
| gao | Fibroblast | 106 | 99 | 7 | -0.0216 | -0.0541 to 0.0108 | 0.1889 | 0.4507 | 0.0463 | 0.1684 | FIT |
| gao | LEC | 36 | 34 | 2 | -0.0219 | -0.0845 to 0.0406 | 0.4785 | 0.5997 | 0.0886 | 0.4388 | FIT |
| gao | Lymphoid | 59 | 56 | 3 | 0.0046 | -0.0129 to 0.0220 | 0.6011 | 0.6511 | 0.0248 | 0.1582 | FIT |
| gao | Mast | 18 | 17 | 1 | 0.0491 | 0.0048 to 0.0935 | 0.0330 | 0.4295 | 0.0618 | 0.6461 | FIT |
| gao | Myeloid | 106 | 99 | 7 | -0.0162 | -0.0435 to 0.0112 | 0.2427 | 0.4507 | 0.0390 | 0.1614 | FIT |
| gao | Neuronal | 60 | 56 | 4 | -0.0019 | -0.0264 to 0.0227 | 0.8794 | 0.8794 | 0.0349 | 0.1992 | FIT |
| gao | Pericyte | 95 | 88 | 7 | 0.0140 | -0.0067 to 0.0347 | 0.1830 | 0.4507 | 0.0295 | 0.1304 | FIT |
| gao | vSMC | 86 | 79 | 7 | 0.0097 | -0.0193 to 0.0386 | 0.5075 | 0.5997 | 0.0412 | 0.1947 | FIT |
| yazar | CD14-low, CD16-positive monocyte | 263 | 257 | 6 | -5.46e-04 | -0.0108 to 0.0097 | 0.9166 | 0.9630 | 0.0147 | 0.0808 | FIT |
| yazar | CD14-positive monocyte | 429 | 408 | 21 | -0.0047 | -0.0117 to 0.0024 | 0.1928 | 0.3855 | 0.0101 | 0.0648 | FIT |
| yazar | CD4-positive, alpha-beta cytotoxic T cell | 93 | 92 | 1 | -0.0195 | -0.0397 to 6.93e-04 | 0.0581 | 0.3855 | 0.0288 | 0.2498 | FIT |
| yazar | central memory CD4-positive, alpha-beta T cell | 955 | 909 | 46 | -0.0018 | -0.0044 to 7.18e-04 | 0.1595 | 0.3855 | 0.0036 | 0.0313 | FIT |
| yazar | central memory CD8-positive, alpha-beta T cell | 111 | 106 | 5 | 0.0012 | -0.0189 to 0.0213 | 0.9069 | 0.9630 | 0.0286 | 0.2388 | FIT |
| yazar | effector memory CD4-positive, alpha-beta T cell | 269 | 258 | 11 | 0.0012 | -0.0077 to 0.0102 | 0.7844 | 0.9630 | 0.0128 | 0.1161 | FIT |
| yazar | effector memory CD8-positive, alpha-beta T cell | 747 | 722 | 25 | -6.14e-04 | -0.0039 to 0.0026 | 0.7103 | 0.9630 | 0.0046 | 0.0454 | FIT |
| yazar | gamma-delta T cell | 129 | 120 | 9 | -2.57e-04 | -0.0113 to 0.0107 | 0.9630 | 0.9630 | 0.0157 | 0.1400 | FIT |
| yazar | memory B cell | 254 | 238 | 16 | 0.0086 | -0.0031 to 0.0203 | 0.1501 | 0.3855 | 0.0167 | 0.1134 | FIT |
| yazar | mucosal-associated invariant T cell | 46 | 38 | 8 | 0.0262 | -0.0321 to 0.0845 | 0.3508 | 0.6237 | 0.0819 | 0.6140 | FIT |
| yazar | naive B cell | 322 | 312 | 10 | -3.54e-04 | -0.0072 to 0.0065 | 0.9194 | 0.9630 | 0.0098 | 0.1373 | FIT |
| yazar | naive thymus-derived CD4-positive, alpha-beta T cell | 821 | 778 | 43 | -0.0022 | -0.0050 to 5.42e-04 | 0.1146 | 0.3855 | 0.0040 | 0.0416 | FIT |
| yazar | naive thymus-derived CD8-positive, alpha-beta T cell | 393 | 352 | 41 | 0.0052 | -8.18e-04 to 0.0113 | 0.0899 | 0.3855 | 0.0086 | 0.0698 | FIT |
| yazar | natural killer cell | 816 | 784 | 32 | -0.0034 | -0.0068 to 4.74e-05 | 0.0532 | 0.3855 | 0.0049 | 0.0505 | FIT |
| yazar | plasmablast | 14 | 14 | 0 | NA | NA to NA | NA | NA | NA | NA | UNTESTABLE |
| yazar | regulatory T cell | 69 | 69 | 0 | 0.0147 | -0.0076 to 0.0370 | 0.1889 | 0.3855 | 0.0316 | 0.3854 | FIT |
| yazar | transitional stage B cell | 213 | 195 | 18 | 0.0020 | -0.0103 to 0.0144 | 0.7447 | 0.9630 | 0.0176 | 0.1412 | FIT |

#### UNTESTABLE primary cell types

| Dataset | Cell type | N orig | N | N lost | In candidate family | Reason |
|---|---|---|---|---|---|---|
| yazar | plasmablast | 14 | 14 | 0 | Yes | insufficient_residual_df |
| yazar | CD16-negative, CD56-bright natural killer cell, human | 6 | 5 | 1 | No | fewer_than_6_donors |
| yazar | CD4-positive, alpha-beta T cell | 0 | 0 | 0 | No | fewer_than_6_donors |
| yazar | CD8-positive, alpha-beta T cell | 0 | 0 | 0 | No | fewer_than_6_donors |
| yazar | conventional dendritic cell | 5 | 5 | 0 | No | fewer_than_6_donors |
| yazar | dendritic cell | 0 | 0 | 0 | No | fewer_than_6_donors |
| yazar | double negative thymocyte | 0 | 0 | 0 | No | fewer_than_6_donors |
| yazar | erythrocyte | 0 | 0 | 0 | No | fewer_than_6_donors |
| yazar | hematopoietic precursor cell | 0 | 0 | 0 | No | fewer_than_6_donors |
| yazar | innate lymphoid cell | 0 | 0 | 0 | No | fewer_than_6_donors |
| yazar | plasmacytoid dendritic cell | 3 | 3 | 0 | No | fewer_than_6_donors |
| yazar | platelet | 0 | 0 | 0 | No | fewer_than_6_donors |

In Yazar primary fits, `plasmablast` (N = 14, N lost = 0, spanning 12 pools plus intercept, age, and sex) remains `UNTESTABLE` (`insufficient_residual_df`) inside the 17-type mature candidate family. `CD16-negative, CD56-bright natural killer cell, human` lost 1 donor under 30 (N = 6 to N = 5) and dropped below the 6-donor floor (`fewer_than_6_donors`), exiting the mature candidate family alongside the 10 types already below the floor in the original analysis.

### Part A sensitivity summary (age 30-plus)

| Dataset | Model | Fitted | UNTESTABLE | Negative hits | Positive hits | Threshold |
|---|---|---|---|---|---|---|
| gao | sensitivity_adjusted | 13 | 0 | 0 | 0 | BH q < 0.05 |
| gao | sensitivity_unadjusted | 13 | 0 | 10 | 0 | nominal p < 0.05 |
| gao | unadjusted | 13 | 0 | 3 | 0 | nominal p < 0.05 |
| yazar | sensitivity_adjusted | 17 | 11 | 0 | 0 | BH q < 0.05 |
| yazar | sensitivity_unadjusted | 19 | 9 | 0 | 14 | nominal p < 0.05 |
| yazar | unadjusted | 17 | 11 | 0 | 9 | nominal p < 0.05 |

#### Sensitivity 1 adjusted fits (age 30-plus, `mean_log_umi_all`)

| Dataset | Cell type | N orig | N | N lost | Slope | 95% CI | p | BH q | MDE80 | MDE80 SD | Status |
|---|---|---|---|---|---|---|---|---|---|---|---|
| gao | Adipocyte | 46 | 44 | 2 | -0.0931 | -0.1682 to -0.0180 | 0.0165 | 0.0715 | 0.1067 | 0.1292 | FIT |
| gao | Cardiomyocyte | 115 | 107 | 8 | -0.0384 | -0.0885 to 0.0117 | 0.1311 | 0.2434 | 0.0714 | 0.1019 | FIT |
| gao | Endocardial | 87 | 80 | 7 | -0.0529 | -0.1042 to -0.0015 | 0.0437 | 0.1136 | 0.0731 | 0.1175 | FIT |
| gao | Endothelial | 111 | 103 | 8 | -0.0382 | -0.0848 to 0.0084 | 0.1070 | 0.2318 | 0.0664 | 0.1104 | FIT |
| gao | Epicardial | 24 | 22 | 2 | 0.0570 | -0.0597 to 0.1737 | 0.3174 | 0.3751 | 0.1644 | 0.3089 | FIT |
| gao | Fibroblast | 117 | 109 | 8 | -0.0654 | -0.1117 to -0.0191 | 0.0061 | 0.0712 | 0.0660 | 0.1110 | FIT |
| gao | LEC | 62 | 59 | 3 | -0.0430 | -0.1110 to 0.0250 | 0.2097 | 0.2937 | 0.0967 | 0.1758 | FIT |
| gao | Lymphoid | 113 | 105 | 8 | -0.0169 | -0.0611 to 0.0272 | 0.4485 | 0.4858 | 0.0629 | 0.1138 | FIT |
| gao | Mast | 73 | 69 | 4 | -0.0052 | -0.0694 to 0.0591 | 0.8731 | 0.8731 | 0.0915 | 0.2217 | FIT |
| gao | Myeloid | 117 | 109 | 8 | -0.0560 | -0.0988 to -0.0131 | 0.0110 | 0.0712 | 0.0611 | 0.1135 | FIT |
| gao | Neuronal | 102 | 95 | 7 | -0.0440 | -0.0864 to -0.0016 | 0.0420 | 0.1136 | 0.0604 | 0.1314 | FIT |
| gao | Pericyte | 117 | 109 | 8 | -0.0297 | -0.0755 to 0.0161 | 0.2015 | 0.2937 | 0.0654 | 0.1146 | FIT |
| gao | vSMC | 115 | 107 | 8 | -0.0285 | -0.0750 to 0.0179 | 0.2259 | 0.2937 | 0.0663 | 0.1304 | FIT |
| yazar | CD14-low, CD16-positive monocyte | 293 | 287 | 6 | -0.0011 | -0.0133 to 0.0112 | 0.8637 | 0.9111 | 0.0175 | 0.0750 | FIT |
| yazar | CD14-positive monocyte | 493 | 468 | 25 | -0.0036 | -0.0128 to 0.0056 | 0.4431 | 0.9111 | 0.0132 | 0.0565 | FIT |
| yazar | CD16-negative, CD56-bright natural killer cell, human | 33 | 30 | 3 | -0.0157 | -0.5515 to 0.5200 | 0.9111 | 0.9111 | 0.6678 | 2.5130 | FIT |
| yazar | CD4-positive, alpha-beta cytotoxic T cell | 280 | 275 | 5 | -0.0101 | -0.0225 to 0.0024 | 0.1114 | 0.4486 | 0.0178 | 0.0589 | FIT |
| yazar | central memory CD4-positive, alpha-beta T cell | 980 | 931 | 49 | 9.11e-04 | -0.0035 to 0.0053 | 0.6828 | 0.9111 | 0.0062 | 0.0227 | FIT |
| yazar | central memory CD8-positive, alpha-beta T cell | 293 | 270 | 23 | 0.0023 | -0.0117 to 0.0162 | 0.7484 | 0.9111 | 0.0199 | 0.0676 | FIT |
| yazar | conventional dendritic cell | 6 | 6 | 0 | NA | NA to NA | NA | NA | NA | NA | UNTESTABLE |
| yazar | effector memory CD4-positive, alpha-beta T cell | 698 | 658 | 40 | 0.0081 | 0.0014 to 0.0149 | 0.0177 | 0.1506 | 0.0096 | 0.0332 | FIT |
| yazar | effector memory CD8-positive, alpha-beta T cell | 966 | 917 | 49 | 0.0017 | -0.0020 to 0.0054 | 0.3663 | 0.8896 | 0.0053 | 0.0205 | FIT |
| yazar | gamma-delta T cell | 329 | 298 | 31 | 0.0018 | -0.0078 to 0.0115 | 0.7127 | 0.9111 | 0.0138 | 0.0512 | FIT |
| yazar | memory B cell | 614 | 573 | 41 | -0.0071 | -0.0156 to 0.0015 | 0.1048 | 0.4486 | 0.0122 | 0.0399 | FIT |
| yazar | mucosal-associated invariant T cell | 112 | 93 | 19 | 0.0037 | -0.0187 to 0.0260 | 0.7417 | 0.9111 | 0.0318 | 0.1268 | FIT |
| yazar | naive B cell | 844 | 796 | 48 | -0.0017 | -0.0067 to 0.0034 | 0.5116 | 0.9111 | 0.0072 | 0.0258 | FIT |
| yazar | naive thymus-derived CD4-positive, alpha-beta T cell | 974 | 925 | 49 | 2.34e-04 | -0.0032 to 0.0037 | 0.8930 | 0.9111 | 0.0049 | 0.0176 | FIT |
| yazar | naive thymus-derived CD8-positive, alpha-beta T cell | 671 | 623 | 48 | 0.0045 | -0.0014 to 0.0104 | 0.1319 | 0.4486 | 0.0084 | 0.0315 | FIT |
| yazar | natural killer cell | 977 | 928 | 49 | 0.0025 | -0.0015 to 0.0066 | 0.2236 | 0.6337 | 0.0058 | 0.0236 | FIT |
| yazar | plasmablast | 14 | 14 | 0 | NA | NA to NA | NA | NA | NA | NA | UNTESTABLE |
| yazar | regulatory T cell | 606 | 576 | 30 | 0.0102 | 0.0019 to 0.0184 | 0.0158 | 0.1506 | 0.0118 | 0.0397 | FIT |
| yazar | transitional stage B cell | 493 | 456 | 37 | 8.72e-04 | -0.0083 to 0.0101 | 0.8521 | 0.9111 | 0.0131 | 0.0443 | FIT |

Within the Yazar `mean_log_umi_all` 19-type candidate family, 2 types are `UNTESTABLE`: yazar: conventional dendritic cell (N=6, lost=0, insufficient_residual_df), yazar: plasmablast (N=14, lost=0, insufficient_residual_df). Below the 6-donor `mean_log_umi_all` floor (`fewer_than_6_donors`), 9 Yazar types are `UNTESTABLE`: yazar: CD4-positive, alpha-beta T cell (N=0, lost=0), yazar: CD8-positive, alpha-beta T cell (N=0, lost=0), yazar: dendritic cell (N=0, lost=0), yazar: double negative thymocyte (N=1, lost=0), yazar: erythrocyte (N=0, lost=0), yazar: hematopoietic precursor cell (N=0, lost=0), yazar: innate lymphoid cell (N=0, lost=0), yazar: plasmacytoid dendritic cell (N=3, lost=0), yazar: platelet (N=2, lost=0).

In Gao `mean_log_umi_all`, 5 adjusted slopes are nominally negative (`Fibroblast` p = 0.0061, q = 0.0712; `Myeloid` p = 0.0110, q = 0.0712; `Adipocyte` p = 0.0165, q = 0.0715; `Neuronal` p = 0.0420, q = 0.1136; `Endocardial` p = 0.0437, q = 0.1136), matching the borderline pooled Gao `mean_log_umi_all` slope (-0.0406, p = 0.0541) in `results/followup-mature-nov17-depth.csv`, but none passes BH q < 0.05. In Yazar `mean_log_umi_all`, both original BH-significant positive adjusted hits (`effector memory CD4-positive, alpha-beta T cell`, originally q = 0.0273, and `naive thymus-derived CD8-positive, alpha-beta T cell`, originally q = 0.0273) lose BH significance at 30-plus (q = 0.1506 and q = 0.4486). Unadjusted age-only slopes remain strongly confounded by study in Gao (10 nominally negative hits in `mean_log_umi_all`, 3 in `mean_log_umi_ge1500`) and by pool in Yazar (14 nominally positive hits in `mean_log_umi_all`, 9 in `mean_log_umi_ge1500`). Full per-type unadjusted estimates are in `results/followup-mature-block6-depth-per-cell-type.csv`.

### Part B: MSSM formal mature four-unit test (age 30-plus)

| Cell type | Model | N orig | N | N lost | Slope | 95% CI | p | BH q | MDE80 | MDE80 SD | Status | Provenance |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GABAergic neuron | primary | 155 | 153 | 2 | -0.1039 | -0.1479 to -0.0598 | 6.90e-06 | 2.76e-05 | 0.0628 | 0.1573 | FIT | reused_verified |
| GABAergic neuron | unadjusted | 155 | 153 | 2 | -0.0954 | -0.1376 to -0.0532 | 1.55e-05 | NA | 0.0602 | 0.1507 | FIT | refit_sensitivity |
| VIP GABAergic cortical interneuron | primary | 197 | 193 | 4 | -0.0877 | -0.1284 to -0.0470 | 3.37e-05 | 6.74e-05 | 0.0581 | 0.1442 | FIT | reused_verified |
| VIP GABAergic cortical interneuron | unadjusted | 197 | 193 | 4 | -0.0859 | -0.1244 to -0.0474 | 1.82e-05 | NA | 0.0550 | 0.1364 | FIT | refit_sensitivity |
| pvalb GABAergic cortical interneuron | primary | 189 | 185 | 4 | -0.0839 | -0.1284 to -0.0393 | 2.74e-04 | 3.65e-04 | 0.0636 | 0.1475 | FIT | reused_verified |
| pvalb GABAergic cortical interneuron | unadjusted | 189 | 185 | 4 | -0.0815 | -0.1236 to -0.0395 | 1.79e-04 | NA | 0.0600 | 0.1391 | FIT | refit_sensitivity |
| sst GABAergic cortical interneuron | primary | 177 | 174 | 3 | -0.0524 | -0.0899 to -0.0150 | 0.0064 | 0.0064 | 0.0535 | 0.1545 | FIT | reused_verified |
| sst GABAergic cortical interneuron | unadjusted | 177 | 174 | 3 | -0.0485 | -0.0836 to -0.0134 | 0.0070 | NA | 0.0501 | 0.1446 | FIT | refit_sensitivity |

All four mature `primary` rows (`age_per_decade + sex_male`, age 30-plus; observed ages 31 to 89) were verified directly against `/home/reza/projects/agescape-revision-data/mssm.h5ad` (all reused: True; maximum absolute difference across all numeric fields = `9.714e-17`) and corrected with BH across the declared four-unit neuronal family. All four units show statistically significant negative depth slopes after four-unit BH correction (`q` from `2.76e-05` to `0.0064`), with high power (`MDE80 SD` from `0.1442` to `0.1573` SD per decade). The unadjusted age-only 30-plus sensitivity (`model = unadjusted`, directly fit from `mssm.h5ad`, no BH q) agrees in sign and magnitude across all four units.

## Independent BH verification

Every BH family was computed with `statsmodels.stats.multitest.multipletests(..., method='fdr_bh')` and checked against an independent step-up BH implementation (`bh_independent`). Maximum absolute q-value differences were zero to machine precision in all five families.

| BH family | Candidate size | Fitted (m) | Max abs q diff |
|---|---|---|---|
| gao_mean_log_umi_ge1500_adjusted_mature30 | 13 | 13 | 1.11e-16 |
| yazar_mean_log_umi_ge1500_adjusted_mature30 | 17 | 16 | 0.00e+00 |
| gao_mean_log_umi_all_adjusted_mature30 | 13 | 13 | 5.55e-17 |
| yazar_mean_log_umi_all_adjusted_mature30 | 19 | 17 | 1.11e-16 |
| mssm_neuronal_mature_primary | 4 | 4 | 0.00e+00 |

Manual spot-checks (at least 3 ranked p-values per family, showing `p * m / k` and step-up monotonic minimum `q_manual`):

| Family | Cell type | Rank k/m | p | p * m / k | q_manual | q_bh | Abs diff |
|---|---|---|---|---|---|---|---|
| gao_mean_log_umi_ge1500_adjusted_mature30 | Mast | 1/13 | 0.033040 | 0.429515 | 0.429515 | 0.429515 | 5.55e-17 |
| gao_mean_log_umi_ge1500_adjusted_mature30 | Adipocyte | 2/13 | 0.110289 | 0.716879 | 0.450733 | 0.450733 | 0.00e+00 |
| gao_mean_log_umi_ge1500_adjusted_mature30 | Myeloid | 7/13 | 0.242702 | 0.450733 | 0.450733 | 0.450733 | 0.00e+00 |
| gao_mean_log_umi_ge1500_adjusted_mature30 | Neuronal | 13/13 | 0.879401 | 0.879401 | 0.879401 | 0.879401 | 1.11e-16 |
| yazar_mean_log_umi_ge1500_adjusted_mature30 | natural killer cell | 1/16 | 0.053236 | 0.851783 | 0.385515 | 0.385515 | 0.00e+00 |
| yazar_mean_log_umi_ge1500_adjusted_mature30 | CD4-positive, alpha-beta cytotoxic T cell | 2/16 | 0.058058 | 0.464462 | 0.385515 | 0.385515 | 0.00e+00 |
| yazar_mean_log_umi_ge1500_adjusted_mature30 | mucosal-associated invariant T cell | 9/16 | 0.350813 | 0.623667 | 0.623667 | 0.623667 | 0.00e+00 |
| yazar_mean_log_umi_ge1500_adjusted_mature30 | gamma-delta T cell | 16/16 | 0.963030 | 0.963030 | 0.963030 | 0.963030 | 0.00e+00 |
| gao_mean_log_umi_all_adjusted_mature30 | Fibroblast | 1/13 | 0.006105 | 0.079369 | 0.071184 | 0.071184 | 0.00e+00 |
| gao_mean_log_umi_all_adjusted_mature30 | Myeloid | 2/13 | 0.010951 | 0.071184 | 0.071184 | 0.071184 | 0.00e+00 |
| gao_mean_log_umi_all_adjusted_mature30 | Cardiomyocyte | 7/13 | 0.131069 | 0.243413 | 0.243413 | 0.243413 | 2.78e-17 |
| gao_mean_log_umi_all_adjusted_mature30 | Mast | 13/13 | 0.873106 | 0.873106 | 0.873106 | 0.873106 | 0.00e+00 |
| yazar_mean_log_umi_all_adjusted_mature30 | regulatory T cell | 1/17 | 0.015823 | 0.268993 | 0.150650 | 0.150650 | 0.00e+00 |
| yazar_mean_log_umi_all_adjusted_mature30 | effector memory CD4-positive, alpha-beta T cell | 2/17 | 0.017724 | 0.150650 | 0.150650 | 0.150650 | 0.00e+00 |
| yazar_mean_log_umi_all_adjusted_mature30 | naive B cell | 9/17 | 0.511581 | 0.966320 | 0.911069 | 0.911069 | 0.00e+00 |
| yazar_mean_log_umi_all_adjusted_mature30 | CD16-negative, CD56-bright natural killer cell, human | 17/17 | 0.911069 | 0.911069 | 0.911069 | 0.911069 | 0.00e+00 |
| mssm_neuronal_mature_primary | GABAergic neuron | 1/4 | 6.90e-06 | 2.76e-05 | 2.76e-05 | 2.76e-05 | 0.00e+00 |
| mssm_neuronal_mature_primary | VIP GABAergic cortical interneuron | 2/4 | 3.37e-05 | 6.74e-05 | 6.74e-05 | 6.74e-05 | 0.00e+00 |
| mssm_neuronal_mature_primary | pvalb GABAergic cortical interneuron | 3/4 | 2.74e-04 | 3.65e-04 | 3.65e-04 | 3.65e-04 | 0.00e+00 |
| mssm_neuronal_mature_primary | sst GABAergic cortical interneuron | 4/4 | 0.006384 | 0.006384 | 0.006384 | 0.006384 | 0.00e+00 |

## Power after 30-plus restriction

Restricting to age 30-plus reduces donor counts and compresses the age span, inflating standardized MDE80 across every fitted primary model that lost donors: mature-to-original `MDE80 SD` ratios range from `1.167` to `1.351` across the 13 Gao primary types, `1.189` to `1.397` across the 15 fitted Yazar primary types that lost donors (`regulatory T cell` had minimum age 42, lost 0 donors, and stayed at ratio `1.000`), and `1.052` to `1.089` across the 4 MSSM units.

Crucially, restriction breaks the frozen `0.50 SD per decade` adequate-power threshold in three fitted primary cell types (no null claim is made for any of them):

- Gao `Epicardial`: N dropped from 20 to 18 (lost 2; `residual_df` 15 to 13); `MDE80` rose from `0.1026` to `0.1251` log UMI per decade, and `MDE80 SD` crossed the threshold from `0.4136` to `0.5014` SD per decade.
- Gao `Mast`: N dropped from 18 to 17 (lost 1; `residual_df` 11 to 10); `MDE80` rose from `0.0451` to `0.0618` log UMI per decade, and `MDE80 SD` crossed the threshold from `0.4780` to `0.6461` SD per decade.
- Yazar `mucosal-associated invariant T cell`: N dropped from 46 to 38 (lost 8; `residual_df` 17 to 14); `MDE80` rose from `0.0608` to `0.0819` log UMI per decade, and `MDE80 SD` crossed the threshold from `0.4448` to `0.6140` SD per decade.

In addition, Yazar `CD16-negative, CD56-bright natural killer cell, human` loses its primary candidate eligibility (`N = 6` to `N = 5`, `fewer_than_6_donors`), and in the `mean_log_umi_all` sensitivity its `MDE80 SD` inflates from `0.7928` (`residual_df = 3`) to `2.5130` (`residual_df = 2`). By contrast, the remaining 11 Gao primary types (`MDE80 SD` `0.1304` to `0.4388`), the remaining 15 fitted Yazar primary types (`MDE80 SD` `0.0313` to `0.3854`), and all 4 MSSM units (`MDE80 SD` `0.1442` to `0.1573`) stay below `0.50 SD per decade`.

## Verdicts under the frozen rules and comparison vs original Block 6

- **Part A, Gao**: **INCONCLUSIVE**. No adjusted primary slope passes BH for a decline, but DOES NOT PROVE is blocked by primary candidate(s) with MDE80 above 0.50 SD per decade: Epicardial (MDE80 SD=0.5014), Mast (MDE80 SD=0.6461).
- **Part A, Yazar**: **INCONCLUSIVE**. No adjusted primary slope passes BH for a decline, but DOES NOT PROVE is blocked by UNTESTABLE primary candidate(s): plasmablast (N=14, insufficient_residual_df); primary candidate(s) with MDE80 above 0.50 SD per decade: mucosal-associated invariant T cell (MDE80 SD=0.6140).
- **Part B, MSSM**: **PROVES**. The mature 30-plus formal four-unit test supports declining depth in GABAergic neuron, VIP GABAergic cortical interneuron, pvalb GABAergic cortical interneuron, sst GABAergic cortical interneuron.

| Arm | Original Block 6 | Mature 30-plus | What survived and what died |
|---|---|---|---|
| Part A, Gao | DOES NOT PROVE | INCONCLUSIVE | Verdict changed due to power loss: still 0 negative BH hits, but Epicardial (MDE80 SD 0.5014) and Mast (0.6461) cross above 0.50 SD per decade. Endocardial nominal negative slope died (p 0.0321 to 0.1151, q 0.4507); Mast remains nominally positive (p 0.0330, q 0.4295). |
| Part A, Yazar | INCONCLUSIVE | INCONCLUSIVE | Verdict survived as INCONCLUSIVE (plasmablast N=14 still UNTESTABLE, plus MAIT MDE80 SD crossed to 0.6140). Still 0 negative BH hits; nominal negative central memory CD4 (p 0.0087 to 0.1595) and NK (p 0.0298 to 0.0532) died. Original BH-significant positive primary hit in naive CD8 T died (q 0.0259 to 0.3855), and both positive mean_log_umi_all BH hits (naive CD8, effector memory CD4) died. CD16-negative CD56-bright NK dropped below 6 donors (6 to 5). |
| Part B, MSSM | PROVES | PROVES | Verdict survived (agrees with original PROVES). All 4 units pass four-unit BH at 30-plus: GABA (-0.1039, q=2.76e-05), VIP (-0.0877, q=6.74e-05), pvalb (-0.0839, q=3.65e-04), sst (-0.0524, q=0.0064). GABA, VIP, and pvalb slopes steepen slightly; sst attenuates slightly (-0.0560 to -0.0524). |

## Caveats

- **MSSM nested-unit caveat**: The four MSSM neuronal units (`GABAergic neuron`, `VIP GABAergic cortical interneuron`, `pvalb GABAergic cortical interneuron`, and `sst GABAergic cortical interneuron`) are nested and overlapping neuronal populations sharing donors rather than four statistically independent cell classes. The four-unit BH family reflects the pre-stated test structure, not independence across units.
- **MSSM missing-covariate caveat**: The MSSM model adjusts only for sex (`age_per_decade + sex_male`). No batch, pool, postmortem interval, RIN, or processing-cohort covariates enter the model. Direct inspection of `/home/reza/projects/agescape-revision-data/mssm.h5ad` `obs` confirmed that `Source` (`M`), `assay` (`10x 3' v3`), and `suspension_type` (`nucleus`) are constant across all 4,140,453 cells, and no technical batch or pool column is available to adjust for potential unmeasured confounding.
- **Batch and pool confounding in Gao and Yazar**: Unadjusted slopes remain heavily confounded by study in Gao (negative bias) and pool in Yazar (positive bias). Study or pool fixed effects remove these apparent trends in every primary cell type.
- **Power and non-equivalence**: MDE80 values are minimum detectable effects for 80 percent power at two-sided alpha 0.05, not equivalence proofs. Where MDE80 exceeds 0.50 SD per decade (`Gao Epicardial`, `Gao Mast`, `Yazar MAIT`), null results are explicitly treated as underpowered rather than evidence of absence.

## Output files

- `reports/followup-mature-block6-depth-extension.md`
- `results/followup-mature-block6-depth-per-cell-type.csv`
- `results/followup-mature-block6-mssm-depth.csv`
- `scripts/run_followup_mature_block6_depth_extension.py`
