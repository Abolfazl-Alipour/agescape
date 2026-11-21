# November 17 CSV tests rerun at 30-plus

Date: 2026-11-21. Pipeline model: kimi-code/k3.1-pro.

This is Phase 2 Item 8. The original November 17 scripts and outputs were not modified. Human mature means age at least 30 years. The new categorical baseline is exactly 30 to 39 years, MID is 40 to 59, and OLD is at least 60. No P1 score was used.

## Depth versus age

### Methods and restriction

The inputs are the pooled donor rows in `results/followup-between-subject-gao-donors.csv` and `results/followup-between-subject-yazar-donors.csv`. Gao retains 109 of 117 donors, losing 8 under 30. Yazar retains 932 of 981, losing 49 under 30. Mature age ranges are 30 to 91 and 30 to 97 years, respectively.

All four original outcomes were refit under Model A, age only, and Model B, age plus sex plus batch. Gao Model B has 6 study dummies and Yazar Model B has 74 pool dummies. Models use OLS with classical standard errors. Slopes and 95 percent CIs are per decade. MDE80 is the absolute slope with 80 percent power for a two-sided 0.05 test, calculated as the sum of the 0.975 and 0.80 residual-df t quantiles times the age standard error. No BH family was used, matching the original uncorrected outcome-level analysis.

### Gao heart, N = 109

| Outcome | Model | Slope | CI low | CI high | p | R2 | MDE80 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| mean_log_umi_ge1500 | A | -0.062820 | -0.110601 | -0.015039 | 0.010455 | 0.059697 | 0.068147 |
| mean_log_umi_ge1500 | B | 0.005631 | -0.015846 | 0.027107 | 0.604105 | 0.843225 | 0.030626 |
| mean_log_umi_all | A | -0.145178 | -0.237192 | -0.053164 | 0.002269 | 0.083770 | 0.131236 |
| mean_log_umi_all | B | -0.040585 | -0.081901 | 0.000732 | 0.054119 | 0.847550 | 0.058919 |
| log_n_cells_all | A | 0.046479 | -0.025779 | 0.118738 | 0.205018 | 0.014969 | 0.103058 |
| log_n_cells_all | B | 0.046500 | -0.005468 | 0.098468 | 0.078901 | 0.579534 | 0.074107 |
| frac_cells_ge1500 | A | -0.047063 | -0.082940 | -0.011186 | 0.010625 | 0.059444 | 0.051169 |
| frac_cells_ge1500 | B | -0.017785 | -0.038526 | 0.002955 | 0.091989 | 0.740597 | 0.029576 |

The original adjusted primary-depth null SURVIVED: slope 0.005631, p 0.604105, compared with the original -0.002107, p 0.807104. Study confounding also SURVIVED: the primary unadjusted slope is -0.062820, p 0.010455, and study plus sex adjustment reverses its sign while R2 rises from 0.059697 to 0.843225. The original adjusted exploratory fraction signal DIED: -0.017785, p 0.091989, compared with -0.017215, p 0.039578. Adjusted captured-cell counts remain non-significant. The all-cell adjusted depth slope is borderline but not significant at 0.05, p 0.054119.

### Yazar blood, N = 932

| Outcome | Model | Slope | CI low | CI high | p | R2 | MDE80 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| mean_log_umi_ge1500 | A | 0.005861 | 0.002156 | 0.009567 | 0.001966 | 0.010255 | 0.005296 |
| mean_log_umi_ge1500 | B | 0.000433 | -0.002150 | 0.003015 | 0.742293 | 0.670420 | 0.003690 |
| mean_log_umi_all | A | 0.037259 | 0.025253 | 0.049266 | 1.64620e-09 | 0.038355 | 0.017157 |
| mean_log_umi_all | B | 0.002645 | -0.001155 | 0.006444 | 0.172210 | 0.933959 | 0.005429 |
| log_n_cells_all | A | -0.013824 | -0.025244 | -0.002405 | 0.017707 | 0.006033 | 0.016318 |
| log_n_cells_all | B | -0.029188 | -0.041178 | -0.017199 | 2.07969e-06 | 0.248593 | 0.017133 |
| frac_cells_ge1500 | A | 0.037999 | 0.026710 | 0.049290 | 6.64552e-11 | 0.044818 | 0.016133 |
| frac_cells_ge1500 | B | 0.001613 | -0.001706 | 0.004932 | 0.340547 | 0.943385 | 0.004743 |

The original adjusted primary-depth null SURVIVED: slope 0.000433, p 0.742293, compared with -0.000266, p 0.806791. Pool confounding SURVIVED: the primary unadjusted slope is positive, p 0.001966, and pool plus sex adjustment removes it while R2 rises from 0.010255 to 0.670420. The captured-cell decline SURVIVED: -0.029188 log cells per decade, p 2.07969e-06, equivalent to 2.877 percent fewer captured cells per decade. The original estimate was -0.027194, p 6.94845e-08. This remains a capture and loading association, not a per-cell depth effect. The adjusted all-cell and exploratory-fraction outcomes remain non-significant.

### Power after restriction

Restriction did not destroy power, but MDE80 increased in every refit. Across all 16 fits, the mature-to-original MDE80 ratio ranges from 1.161 to 1.266 in Gao and from 1.207 to 1.225 in Yazar. For adjusted primary depth, MDE80 rose from 0.024336 to 0.030626 log UMI per decade in Gao and from 0.003052 to 0.003690 in Yazar. For adjusted captured counts, it rose from 0.063101 to 0.074107 log cells in Gao and from 0.014026 to 0.017133 in Yazar. These are minimum detectable effects, not equivalence bounds or evidence for a null biological effect.

Full-precision refits are in `results/followup-mature-nov17-depth.csv`.

## Donor distance variance

### Distance rebuild and reference audit

Neither committed donor table contained a 30-to-39-reference distance. The columns `a1_dist_full_1500` use a 20-to-39 reference and `a1_dist_mature_ref_mid_loo` uses a 40-to-59 reference. New distances were therefore recomputed from the raw h5ads in `scripts/run_mature_nov17_donor_spread.py` using the original `reference_without_donor` and `hellinger_distance` functions.

The rebuild uses the original selected donors and cell types, count gate, gene projection and exclusions, post-exclusion 1,500-UMI cell filter, full-depth counts, and 20-cell donor floor. Each reference unit requires at least 2 eligible 30-to-39 donors and at least 50 cells after every leave-one-out subtraction. Donors aged 30 to 39 use leave-one-out references; older donors use the fixed 30-to-39 reference. Recomputed mature cell counts matched every corresponding committed donor-by-group count. Gao has 12 reference donors and Yazar has 59. All scored donors are at least 30. No donor under 30 enters a reference or scored arm.

Levene tests are median-centered. Variance-ratio CIs use 2,000 donor bootstraps within bins, seed 7. Yazar residualization subtracts the pool mean and adds the grand mean within each mature distance version. No mature pooled singleton pools required exclusion. The original OLD-versus-MID distance using the 40-to-59 reference is retained as a separate arm. No BH family was used, matching the original exploratory supplement.

### Main comparisons

| Dataset | Reference and adjustment | Contrast | N A | N B | Ratio | Bootstrap CI | Levene p |
| --- | --- | --- | ---: | ---: | ---: | --- | ---: |
| Gao | 30-39, raw | OLD vs 30-39 | 41 | 12 | 1.048484 | 0.387525 to 6.456563 | 0.622569 |
| Gao | 30-39, raw | MID vs 30-39 | 52 | 12 | 0.610971 | 0.250148 to 3.546018 | 0.705604 |
| Gao | 30-39, raw | OLD vs MID | 41 | 52 | 1.716095 | 0.772375 to 3.365498 | 0.149203 |
| Gao | retained 40-59, raw | OLD vs MID | 41 | 52 | 1.679106 | 0.832915 to 3.087699 | 0.128068 |
| Yazar | 30-39, raw | OLD vs 30-39 | 680 | 59 | 1.542756 | 1.086570 to 2.686931 | 3.46120e-05 |
| Yazar | 30-39, raw | MID vs 30-39 | 193 | 59 | 1.386762 | 0.940371 to 2.481605 | 0.126695 |
| Yazar | 30-39, raw | OLD vs MID | 680 | 193 | 1.112488 | 0.927739 to 1.367805 | 0.000823 |
| Yazar | retained 40-59, raw | OLD vs MID | 680 | 193 | 1.127973 | 0.941228 to 1.406568 | 0.000472 |
| Yazar | 30-39, pool-residualized | OLD vs 30-39 | 680 | 59 | 1.578469 | 1.074873 to 2.562159 | 0.136603 |
| Yazar | 30-39, pool-residualized | MID vs 30-39 | 193 | 59 | 1.217836 | 0.751805 to 2.073760 | 0.706118 |
| Yazar | 30-39, pool-residualized | OLD vs MID | 680 | 193 | 1.296126 | 0.981035 to 1.761226 | 0.059663 |
| Yazar | retained 40-59, pool-residualized | OLD vs MID | 680 | 193 | 1.287118 | 0.981989 to 1.738962 | 0.065199 |

Gao's original OLD-versus-YOUNG null SURVIVED as an OLD-versus-30-to-39 null, but precision is worse with only 12 baseline donors. The original ratio was 0.785256, p 0.847291, CI 0.332694 to 2.866719. The mature ratio is 1.048484, p 0.622569, CI 0.387525 to 6.456563. This is inconclusive about large variance changes, not proof of no change.

Yazar's raw OLD-versus-YOUNG widening SURVIVED as raw OLD-versus-30-to-39 widening: ratio 1.542756, p 3.46120e-05, compared with the original 1.620755, p 1.87336e-07. The original absence of pool-residualized Levene support also persisted: p 0.136603, compared with 0.430585. However, the new residualized ratio is 1.578469 and its bootstrap CI is 1.074873 to 2.562159, excluding 1, unlike the original CI 0.943458 to 1.756557. Thus the adjusted arm is method-discrepant and does not establish a robust Levene-supported increase. The bootstrap result is reported rather than hidden.

The retained OLD-versus-MID arm remains distinct. Gao retained-reference p 0.128068 and Yazar raw retained-reference p 0.000471846 reproduce their original point estimates and Levene results. The new 30-to-39-reference OLD-versus-MID arm has the same mixed pattern in Yazar: significant raw Levene p 0.000823 with a bootstrap CI including 1, and pool-residualized p 0.059663 with a CI including 1.

### Per-bin descriptives

| Dataset | Reference and adjustment | Bin | N | Mean | SD | IQR |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| Gao | 30-39, raw | 30-39 | 12 | 0.640185 | 0.053648 | 0.053864 |
| Gao | 30-39, raw | MID | 52 | 0.613396 | 0.041934 | 0.044987 |
| Gao | 30-39, raw | OLD | 41 | 0.632678 | 0.054933 | 0.077577 |
| Gao | retained 40-59, raw | MID | 52 | 0.607254 | 0.046971 | 0.056355 |
| Gao | retained 40-59, raw | OLD | 41 | 0.628931 | 0.060866 | 0.087062 |
| Yazar | 30-39, raw | 30-39 | 59 | 0.583393 | 0.009321 | 0.008018 |
| Yazar | 30-39, raw | MID | 193 | 0.582269 | 0.010976 | 0.015278 |
| Yazar | 30-39, raw | OLD | 680 | 0.579917 | 0.011577 | 0.019878 |
| Yazar | retained 40-59, raw | MID | 193 | 0.581020 | 0.010844 | 0.015059 |
| Yazar | retained 40-59, raw | OLD | 680 | 0.578485 | 0.011517 | 0.019566 |
| Yazar | 30-39, pool-residualized | 30-39 | 59 | 0.579514 | 0.003331 | 0.004413 |
| Yazar | 30-39, pool-residualized | MID | 193 | 0.579826 | 0.003676 | 0.003955 |
| Yazar | 30-39, pool-residualized | OLD | 680 | 0.580946 | 0.004185 | 0.004842 |
| Yazar | retained 40-59, pool-residualized | MID | 193 | 0.578255 | 0.003717 | 0.004354 |
| Yazar | retained 40-59, pool-residualized | OLD | 680 | 0.579270 | 0.004216 | 0.004864 |

Full quartiles and medians are in `results/followup-mature-nov17-donor-spread-descriptives.csv`; full-precision comparisons are in `results/followup-mature-nov17-donor-spread.csv`.

### Yazar per-cell-type supplement

Six cell types met the 20-donor-per-bin floor on the new raw distance. All six also have pool-residualized results, unlike the original raw-only supplement. Memory B, naive B, and transitional stage B no longer met the raw baseline floor. No Gao cell type met the floor. These tests remain exploratory and uncorrected.

| Cell type | Adjustment | N OLD | N 30-39 | Ratio | Bootstrap CI | Levene p |
| --- | --- | ---: | ---: | ---: | --- | ---: |
| CD14-positive monocyte | raw | 299 | 21 | 1.772442 | 1.074240 to 3.898490 | 0.125807 |
| CD14-positive monocyte | pool-residualized | 295 | 20 | 1.347622 | 0.970150 to 2.215343 | 0.817656 |
| central memory CD4 T | raw | 660 | 58 | 1.348278 | 0.992807 to 2.167575 | 2.72612e-05 |
| central memory CD4 T | pool-residualized | 660 | 58 | 1.496122 | 1.130346 to 2.116129 | 0.250524 |
| effector memory CD8 T | raw | 552 | 34 | 0.852559 | 0.627278 to 1.418115 | 0.706404 |
| effector memory CD8 T | pool-residualized | 552 | 33 | 0.781564 | 0.487876 to 1.461710 | 0.294171 |
| naive CD4 T | raw | 564 | 51 | 1.294862 | 0.953827 to 2.007567 | 0.012317 |
| naive CD4 T | pool-residualized | 564 | 50 | 0.937978 | 0.629331 to 1.598016 | 0.816469 |
| naive CD8 T | raw | 210 | 40 | 0.805379 | 0.571684 to 1.360348 | 0.562818 |
| naive CD8 T | pool-residualized | 206 | 39 | 1.442608 | 0.950757 to 2.560633 | 0.275995 |
| natural killer cell | raw | 598 | 38 | 1.039735 | 0.713882 to 2.003945 | 0.093842 |
| natural killer cell | pool-residualized | 598 | 38 | 1.435060 | 0.980301 to 2.424461 | 0.566448 |

The original raw central-memory CD4 and naive-CD4 nominal Levene signals persist at 2.72612e-05 and 0.012317. The original natural-killer signal DIED at p 0.093842. All six pool-residualized Levene tests are non-significant. Central-memory CD4 has a residualized bootstrap CI above 1 despite Levene p 0.250524, another method discrepancy. Singleton-pool exclusions in the supplement are listed per row in the CSV.

## Survivor-bias correlation

### Mature join and methods

Per-unit drift slopes were regenerated from the Item 4 full-depth mature donor endpoints with the original `fit_unit` implementation in `scripts/run_cross_species_expression_vector.py`. Mouse rows were restricted to full depth and mice aged at least 18 months from the 18-month reference table. All human scored donors were at least 30 years old and all human reference windows were exactly 30 to 39 years. The primary model regresses profile Hellinger distance on age, sex, and mean log UMI; the sensitivity model omits mean log UMI. Slopes are per year.

Composition inputs are the mature-only, age-plus-sex CLR slopes in `results/followup-mature-composition-mouse.csv` and `results/followup-mature-composition-human.csv`. The exact inner-join keys are `species`, `organ` or `tissue`, `cohort`, and `cell_type`. The join retained 51 of 114 mouse composition units and 87 of 285 human composition units, 138 units total. Every joined unit has both depth-adjusted and unadjusted drift slopes. The legacy joined CSV was neither read nor filtered.

Pearson, Spearman, absolute-slope, and unadjusted-sensitivity tests were rerun separately by species with two-sided p-values and 10,000 label shuffles, seed 20261117, using the original test order and permutation implementation. Share loss is minus the CLR fraction slope, so survivor bias predicts a negative correlation between share loss and drift. The complete 20-row result table is `results/followup-mature-nov17-survivor.csv`; the regenerated units are `results/followup-mature-nov17-survivor-joined.csv`; the pipeline is `scripts/run_mature_nov17_survivor.py`.

### Results and verdict

| Scope | Test | Pearson r | p | Permutation p | Spearman rho | p | Permutation p |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Mouse, N = 51 | Share loss vs drift | -0.057408 | 0.689053 | 0.668133 | +0.138733 | 0.331600 | 0.328867 |
| Mouse, N = 51 | Share loss vs absolute drift | +0.227395 | 0.108545 | 0.114689 | +0.072127 | 0.614983 | 0.616038 |
| Mouse, N = 51 | Fraction slope vs drift, unadjusted | +0.128605 | 0.368440 | 0.357164 | -0.036380 | 0.799921 | 0.803520 |
| Human, N = 87 | Share loss vs drift | -0.195875 | 0.069029 | 0.075092 | -0.100587 | 0.353929 | 0.353865 |
| Human, N = 87 | Share loss vs absolute drift | +0.287452 | 0.006942 | 0.007999 | +0.142925 | 0.186627 | 0.195580 |
| Human, N = 87 | Fraction slope vs drift, unadjusted | +0.051811 | 0.633652 | 0.617338 | -0.096741 | 0.372723 | 0.376662 |

The original verdict of no support for survivor bias SURVIVED. Neither species has a significant primary Pearson or Spearman correlation. Mouse Pearson is near zero and mouse Spearman is in the opposite direction from the prediction. Human Pearson is in the predicted direction but is not significant by either its analytic or permutation p-value, and human Spearman is null. The significant human Pearson correlation between share loss and absolute drift is positive, opposite the survivor-bias prediction; its Spearman counterpart is not significant.

The original nominal mouse share-loss Pearson signal DIED: r changed from +0.250201, p 0.049843, permutation p 0.048795, N 62, to -0.057408, p 0.689053, permutation p 0.668133, N 51. The original human share-loss Pearson was +0.078340, p 0.366446, permutation p 0.211879, N 135; it changed sign but remains non-significant. The original human share-loss Spearman was +0.151210, p 0.080008, permutation p 0.073193; the mature result is -0.100587, p 0.353929, permutation p 0.353865. The original absolute-drift tests were all null; the new human Pearson absolute-drift result is therefore new, not a survival of an original signal.

Descriptively, mouse has 30 losing and 21 gaining units. Mean signed drift is -0.002502345 per year in losing units versus -0.004861982 in gaining units; mean absolute drift is 0.017079600 versus 0.012472379. Human has 45 losing and 42 gaining units. Mean signed drift is -0.000106741 per year in losing units versus +0.000047542 in gaining units; mean absolute drift is 0.000431215 versus 0.000284267. Thus losing types have greater mean absolute drift in both species, despite the human signed-mean difference.

These remain ecological correlations across selected units that share donors, not causal evidence or a within-type cell-survival test. The smaller mature joined set and non-significant primary tests do not establish absence of survivor bias.

## Runtime and verification

| Pipeline | Start UTC | End UTC | Duration seconds | Status |
| --- | --- | --- | ---: | --- |
| depth | 2026-11-21T20:33:07.616654+00:00 | 2026-11-21T20:33:07.914705+00:00 | 0.298067 | completed |
| Gao distance rebuild | 2026-11-21T20:38:19.053456+00:00 | 2026-11-21T20:39:42.893177+00:00 | 83.839706 | completed |
| Yazar distance rebuild | 2026-11-21T20:39:43.701828+00:00 | 2026-11-21T20:41:40.825476+00:00 | 117.123645 | completed |
| complete donor-spread pipeline | 2026-11-21T20:38:15.609731+00:00 | 2026-11-21T20:41:42.779147+00:00 | 207.169438 | completed |
| survivor bias | 2026-11-21T20:58:20.202929+00:00 | 2026-11-21T20:59:02.815794+00:00 | 42.612863 | completed |

The first donor-spread launch failed before analysis because the repository import path was missing. The corrected run above completed. Maximum reported resident memory was 2,183.6 MB. The raw-data distance audit table is `results/followup-mature-nov17-donor-spread-distances.csv`. The survivor-bias pipeline used only committed composition and donor CSVs, loaded no count matrices, and did not read the legacy joined CSV. Its model was kimi-code/k3.1-pro.

One self-review pass checked all script age filters, reference construction, output ages, reference ranges, donor and reference floors, CI ordering, variance-ratio arithmetic, descriptives, supplement floors, and forbidden dash characters. New distances differ from the young-reference distances for every donor-group pair where both exist, so no full-range distance was relabeled as mature. All committed mature cell counts matched the raw-data rebuild. The survivor-bias review recomputed every joined slope from the mature composition and donor tables, verified the exact keys and one-to-one joins, confirmed all 51 mouse and 87 human fitted drift units matched composition, recomputed all analytic correlations and p-values, and checked mature ages, reference windows, output schemas, and forbidden dash characters. No age leakage, wrong join keys, duplicated full-range estimates, or other real errors remained after this pass. No completed arm is unverified.
