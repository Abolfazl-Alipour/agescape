# Block 6: Depth versus age extension

## Status

Complete. The specification and correction families below were frozen before model fitting.

## Inputs

- Gao donor table: `results/followup-between-subject-gao-donors.csv`.
- Yazar donor table: `results/followup-between-subject-yazar-donors.csv`.
- MSSM h5ad: `/home/reza/workspace/rev2-mssm-null-scratch/mssm.h5ad`, a symlink to the data root.
- The MSSM neuronal unit list was enumerated before fitting from `results/reviewer-rev2-human-de-summary.csv` and matches the four MSSM rows in the drift-analysis audit `results/h3-human-vector-design-audit.csv`. It was also cross-checked against `INTERNEURON_TYPES` in `scripts/run_rev2_de_interneuron.py`: `GABAergic neuron`, `VIP GABAergic cortical interneuron`, `pvalb GABAergic cortical interneuron`, and `sst GABAergic cortical interneuron`.
- The human mature cutoff of 30 years and older comes from the header of `scripts/run_mature_only_composition.py`.

The November 17 headline is cited, not recomputed: Gao pooled adjusted slope -0.0021 per decade, p 0.81; Yazar pooled adjusted slope -0.0003 per decade, p 0.81. Yazar captured cell counts change by -2.7 percent per decade, p 6.9e-08. Those pooled rows are excluded from Part A.

## Exact analysis specification

### Part A: per-cell-type depth slopes in Gao and Yazar

For each non-pooled cell type with at least 6 eligible donors, fit donor-level OLS of `mean_log_umi_ge1500` on `age_per_decade + sex_male + batch`. Gao batch is the seven-study fixed effect, represented by an intercept and six study indicators. Yazar batch is pool fixed effects. A donor is eligible for the primary outcome with at least 20 cells in `n_cells_ge1500` and nonmissing outcome, age, sex, and batch. The primary slope is the change in natural-log UMI per decade. Report slope, 95 percent CI, p-value, N donors, and MDE80 for every cell type.

Sensitivity 1 uses the same adjusted model with `mean_log_umi_all` as outcome and at least 20 cells in `n_cells_all` per donor. Sensitivity 2 is an unadjusted age-only model for each outcome on its outcome-specific eligible donors. No pooled rows, cell-level fits, downsampling, or rarefaction are used.

The primary multiple-testing family is fixed before fitting as all eligible non-pooled cell types within each dataset, with BH correction separately in Gao and Yazar. Candidate family sizes are 13 Gao and 18 Yazar cell types. BH is calculated over successfully fitted members of these fixed candidate families. The adjusted `mean_log_umi_all` sensitivity is a separate family within each dataset, with candidate sizes 13 Gao and 19 Yazar. Unadjusted age-only fits receive no q-values and are not mixed into either family.

### Part B: formal MSSM neuron test

For each of the four neuronal units listed above, calculate per-donor mean natural-log UMI from full-depth per-cell `n_counts`, using all cells passing the existing QC route in `fetch_unit_h5ad`: primary data, disease `normal`, `n_counts >= 1500`, parsed age 20 to 120 years, and at least 20 eligible cells per donor per unit. There is no depth downsampling. The h5ad `n_counts` field is the per-cell raw library size used by the existing route, so the sparse expression matrix is not needed. A unit requires at least 6 eligible donors.

The primary donor-level OLS is mean log UMI on `age_per_decade + sex_male`. Report slope per decade, 95 percent CI, p-value, N donors, and MDE80. Sensitivity 1 is an unadjusted age-only model on all eligible donors. Sensitivity 2 repeats the primary age-plus-sex model restricted to donors age 30 years and older. The primary family is fixed as the four MSSM neuronal units, with BH correction across those units. Sensitivity fits are not added to this family.

For both parts, sex is omitted only when it is constant or rank-deficient. Batch indicators are included at all observed levels through an intercept and k minus 1 indicators. Collinear indicators are omitted without dropping age. A model without residual degrees of freedom is UNTESTABLE, not null. All observations are donors, and each donor contributes one value per unit.

## MDE80 and verdict rules fixed before fitting

MDE80 is the repository convention `(t_0.975 + t_0.80) * SE` on the residual degrees of freedom, reported in natural-log UMI per decade and in outcome SD per decade. Following the existing follow-up convention, MDE80 at most 0.50 SD per decade is small enough for a null to count as adequately powered; at most 0.25 SD per decade is high power. Both raw and standardized MDE80 values are retained.

For each dataset in Part A, PROVES means at least one batch-adjusted primary cell-type slope is negative with BH q below 0.05. If no cell type proves a decline, DOES NOT PROVE requires all fitted primary slopes to be null or nonnegative with MDE80 at most 0.50 SD per decade, with no UNTESTABLE primary candidate. Otherwise the dataset is INCONCLUSIVE. For Part B, PROVES means at least one of the four primary neuronal slopes is negative with BH q below 0.05. If no unit proves a decline, DOES NOT PROVE requires all four primary fits to be null or nonnegative with MDE80 at most 0.50 SD per decade. Otherwise Part B is INCONCLUSIVE. Sensitivities are reported alongside but do not replace these primary verdict rules.

Fixed seed: 20261120. No random sampling, permutation, or bootstrap is used, so every fit is deterministic. The human technical gene mask is not applicable because this block uses library sizes, not gene-level expression.

## Results

All slopes and MDE80 values are natural-log UMI per decade. FIT and UNTESTABLE describe model evaluability, not scientific verdicts.

### Part A primary fits

| Dataset | Cell type | N | Slope | 95% CI | p | BH q | MDE80 | MDE80 SD | Status |
|---|---|---|---|---|---|---|---|---|---|
| gao | Adipocyte | 37 | -0.0622 | -0.1279 to 0.0036 | 0.0631 | 0.2733 | 0.0933 | 0.2313 | FIT |
| gao | Cardiomyocyte | 109 | -0.0217 | -0.0635 to 0.0200 | 0.3042 | 0.4394 | 0.0595 | 0.1239 | FIT |
| gao | Endocardial | 72 | -0.0374 | -0.0716 to -0.0033 | 0.0321 | 0.2474 | 0.0486 | 0.1787 | FIT |
| gao | Endothelial | 87 | 0.0116 | -0.0039 to 0.0272 | 0.1391 | 0.4394 | 0.0221 | 0.1080 | FIT |
| gao | Epicardial | 20 | -0.0088 | -0.0818 to 0.0641 | 0.7995 | 0.7995 | 0.1026 | 0.4136 | FIT |
| gao | Fibroblast | 106 | -0.0156 | -0.0412 to 0.0101 | 0.2312 | 0.4394 | 0.0366 | 0.1369 | FIT |
| gao | LEC | 36 | -0.0150 | -0.0674 to 0.0374 | 0.5630 | 0.7217 | 0.0743 | 0.3751 | FIT |
| gao | Lymphoid | 59 | 0.0084 | -0.0064 to 0.0232 | 0.2616 | 0.4394 | 0.0210 | 0.1356 | FIT |
| gao | Mast | 18 | 0.0345 | 0.0023 to 0.0667 | 0.0381 | 0.2474 | 0.0451 | 0.4780 | FIT |
| gao | Myeloid | 106 | -0.0114 | -0.0329 to 0.0102 | 0.2969 | 0.4394 | 0.0307 | 0.1308 | FIT |
| gao | Neuronal | 60 | 0.0048 | -0.0159 to 0.0254 | 0.6463 | 0.7217 | 0.0294 | 0.1675 | FIT |
| gao | Pericyte | 95 | 0.0095 | -0.0064 to 0.0255 | 0.2381 | 0.4394 | 0.0227 | 0.1039 | FIT |
| gao | vSMC | 86 | 0.0047 | -0.0170 to 0.0265 | 0.6662 | 0.7217 | 0.0310 | 0.1525 | FIT |
| yazar | CD14-low, CD16-positive monocyte | 263 | -0.0031 | -0.0116 to 0.0055 | 0.4763 | 0.7908 | 0.0122 | 0.0668 | FIT |
| yazar | CD14-positive monocyte | 429 | -0.0023 | -0.0081 to 0.0034 | 0.4227 | 0.7908 | 0.0082 | 0.0529 | FIT |
| yazar | CD16-negative, CD56-bright natural killer cell, human | 6 | NA | NA to NA | NA | NA | NA | NA | UNTESTABLE |
| yazar | CD4-positive, alpha-beta cytotoxic T cell | 93 | -0.0117 | -0.0279 to 0.0045 | 0.1520 | 0.4863 | 0.0230 | 0.2004 | FIT |
| yazar | central memory CD4-positive, alpha-beta T cell | 955 | -0.0028 | -0.0049 to -7.09e-04 | 0.0087 | 0.0697 | 0.0030 | 0.0260 | FIT |
| yazar | central memory CD8-positive, alpha-beta T cell | 111 | 0.0027 | -0.0139 to 0.0193 | 0.7470 | 0.8614 | 0.0237 | 0.1917 | FIT |
| yazar | effector memory CD4-positive, alpha-beta T cell | 269 | -0.0011 | -0.0085 to 0.0062 | 0.7635 | 0.8614 | 0.0105 | 0.0957 | FIT |
| yazar | effector memory CD8-positive, alpha-beta T cell | 747 | -8.09e-04 | -0.0035 to 0.0019 | 0.5562 | 0.8089 | 0.0039 | 0.0380 | FIT |
| yazar | gamma-delta T cell | 129 | 0.0028 | -0.0053 to 0.0109 | 0.4943 | 0.7908 | 0.0115 | 0.1032 | FIT |
| yazar | memory B cell | 254 | 0.0074 | -0.0015 to 0.0163 | 0.1021 | 0.4086 | 0.0127 | 0.0880 | FIT |
| yazar | mucosal-associated invariant T cell | 46 | 0.0055 | -0.0376 to 0.0487 | 0.7897 | 0.8614 | 0.0608 | 0.4448 | FIT |
| yazar | naive B cell | 322 | 3.54e-04 | -0.0054 to 0.0061 | 0.9031 | 0.9031 | 0.0082 | 0.1155 | FIT |
| yazar | naive thymus-derived CD4-positive, alpha-beta T cell | 821 | -0.0013 | -0.0036 to 9.47e-04 | 0.2552 | 0.5834 | 0.0032 | 0.0342 | FIT |
| yazar | naive thymus-derived CD8-positive, alpha-beta T cell | 393 | 0.0073 | 0.0028 to 0.0118 | 0.0016 | 0.0259 | 0.0064 | 0.0520 | FIT |
| yazar | natural killer cell | 816 | -0.0032 | -0.0060 to -3.12e-04 | 0.0298 | 0.1588 | 0.0041 | 0.0422 | FIT |
| yazar | plasmablast | 14 | NA | NA to NA | NA | NA | NA | NA | UNTESTABLE |
| yazar | regulatory T cell | 69 | 0.0147 | -0.0076 to 0.0370 | 0.1889 | 0.5038 | 0.0316 | 0.3854 | FIT |
| yazar | transitional stage B cell | 213 | 0.0011 | -0.0078 to 0.0100 | 0.8076 | 0.8614 | 0.0127 | 0.1011 | FIT |

Cell types below the six-donor primary floor are UNTESTABLE: yazar: CD4-positive, alpha-beta T cell (N=0), yazar: CD8-positive, alpha-beta T cell (N=0), yazar: conventional dendritic cell (N=5), yazar: dendritic cell (N=0), yazar: double negative thymocyte (N=0), yazar: erythrocyte (N=0), yazar: hematopoietic precursor cell (N=0), yazar: innate lymphoid cell (N=0), yazar: plasmacytoid dendritic cell (N=3), yazar: platelet (N=0).

### Part A sensitivity summary

| Dataset | Model | Fitted | UNTESTABLE | Negative hits | Positive hits | Threshold |
|---|---|---|---|---|---|---|
| gao | sensitivity_adjusted | 13 | 0 | 0 | 0 | BH q < 0.05 |
| gao | sensitivity_unadjusted | 13 | 0 | 5 | 0 | nominal p < 0.05 |
| gao | unadjusted | 13 | 0 | 1 | 0 | nominal p < 0.05 |
| yazar | sensitivity_adjusted | 17 | 11 | 0 | 2 | BH q < 0.05 |
| yazar | sensitivity_unadjusted | 19 | 9 | 0 | 15 | nominal p < 0.05 |
| yazar | unadjusted | 18 | 10 | 0 | 9 | nominal p < 0.05 |

The adjusted `mean_log_umi_all` fits use their separate pre-stated BH families. Unadjusted fits are age-only and are labeled with nominal p-values. Full per-unit sensitivity estimates, CIs, and MDE80 values are in the Part A CSV.

### Part B MSSM formal test

| Cell type | Model | N | Slope | 95% CI | p | BH q | MDE80 | MDE80 SD | Status |
|---|---|---|---|---|---|---|---|---|---|
| GABAergic neuron | mature_age30_plus | 153 | -0.1039 | -0.1479 to -0.0598 | 6.90e-06 | NA | 0.0628 | 0.1573 | FIT |
| GABAergic neuron | primary | 155 | -0.0985 | -0.1402 to -0.0569 | 6.50e-06 | 2.60e-05 | 0.0595 | 0.1495 | FIT |
| GABAergic neuron | unadjusted | 155 | -0.0915 | -0.1317 to -0.0513 | 1.34e-05 | NA | 0.0574 | 0.1443 | FIT |
| VIP GABAergic cortical interneuron | mature_age30_plus | 193 | -0.0877 | -0.1284 to -0.0470 | 3.37e-05 | NA | 0.0581 | 0.1442 | FIT |
| VIP GABAergic cortical interneuron | primary | 197 | -0.0663 | -0.1044 to -0.0281 | 7.49e-04 | 0.0015 | 0.0545 | 0.1329 | FIT |
| VIP GABAergic cortical interneuron | unadjusted | 197 | -0.0684 | -0.1050 to -0.0318 | 2.99e-04 | NA | 0.0523 | 0.1275 | FIT |
| pvalb GABAergic cortical interneuron | mature_age30_plus | 185 | -0.0839 | -0.1284 to -0.0393 | 2.74e-04 | NA | 0.0636 | 0.1475 | FIT |
| pvalb GABAergic cortical interneuron | primary | 189 | -0.0599 | -0.1015 to -0.0183 | 0.0050 | 0.0050 | 0.0594 | 0.1354 | FIT |
| pvalb GABAergic cortical interneuron | unadjusted | 189 | -0.0620 | -0.1018 to -0.0222 | 0.0025 | NA | 0.0569 | 0.1297 | FIT |
| sst GABAergic cortical interneuron | mature_age30_plus | 174 | -0.0524 | -0.0899 to -0.0150 | 0.0064 | NA | 0.0535 | 0.1545 | FIT |
| sst GABAergic cortical interneuron | primary | 177 | -0.0560 | -0.0910 to -0.0210 | 0.0019 | 0.0025 | 0.0499 | 0.1435 | FIT |
| sst GABAergic cortical interneuron | unadjusted | 177 | -0.0528 | -0.0858 to -0.0198 | 0.0019 | NA | 0.0471 | 0.1354 | FIT |

The primary q-values are corrected across exactly four neuronal units. The age-only and mature 30 years and older fits are sensitivities and have no q-values.

## Verdicts

- Part A, gao: DOES NOT PROVE. No adjusted primary slope passes BH for a decline, and every candidate has MDE80 at most 0.50 SD per decade.
- Part A, yazar: INCONCLUSIVE. No adjusted primary slope passes BH for a decline, but at least one candidate is UNTESTABLE or has MDE80 above 0.50 SD per decade.
- Part B, MSSM: PROVES. The direct formal test supports declining depth in GABAergic neuron, VIP GABAergic cortical interneuron, pvalb GABAergic cortical interneuron, sst GABAergic cortical interneuron.

Does the MSSM neuron depth decline survive a direct formal test? PROVES. The direct formal test supports declining depth in GABAergic neuron, VIP GABAergic cortical interneuron, pvalb GABAergic cortical interneuron, sst GABAergic cortical interneuron.

The Part A verdicts concern per-cell-type depth slopes and do not repeat the pooled November 17 tests. The Part B verdict concerns the four pre-stated neuronal units only.

## Output files

- `results/followup-11pt-block6-depth-per-cell-type.csv`
- `results/followup-11pt-block6-mssm-depth.csv`
- `scripts/run_followup_11pt_block6_depth_extension.py`
