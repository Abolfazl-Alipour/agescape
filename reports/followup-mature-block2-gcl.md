# Block 2 (mature-to-old): Gene coordination, mature vs old, donor level

Status: complete. This specification was written before running the analysis. Mature-to-old rerun of reports/followup-11pt-block2-gcl.md: donors aged under 30 are excluded, detected-gene sets and fixed gene pairs are rebuilt from mature (30-plus) donors only, and the group comparison baseline YOUNG (20 to 39) is replaced by MATURE (30 to 39). Every other frozen rule is unchanged.

## Exact analysis specification

- Gene set, fixed before running: the detected-gene set per cell type, defined as genes detected in at least 10 percent of cells in that cell type pooled across mature (30-plus) donors only. Record N genes per cell type. Do not change the gene set after seeing results.
- Per donor per cell type, using the shared eligibility floors: compute the gene coordination level (GCL) as the mean pairwise Spearman correlation between genes across that donor's cells, using full-depth within-cell proportions. For cell types with more than 2000 genes in the set, subsample gene PAIRS, not genes, with a fixed seed so the estimator is identical across donors. Record the pair count and seed.
- Primary model, per dataset and cell type: OLS of per-donor GCL on age_per_decade + sex_male + batch, with Gao study indicators and Yazar pool fixed effects. Report slope per decade, 95 percent CI, p-value, N donors, and MDE80.
- Group comparison secondary: mean GCL in OLD vs MATURE donors per cell type, donor-permutation null with 2000 permutations and a fixed seed, using batch-residualized GCL after residualizing on batch first and then permuting donor age-bin labels.
- Split-half technical null: split each donor's cells into two halves with 5 fixed seeds, compute GCL per half, and report the half-to-half noise floor per cell type. A claimed age effect must exceed this floor. Report the floor next to the slope.
- Depth sensitivity: refit the primary model adding per-donor mean log UMI as a covariate.
- A labeled 1500 rarefaction sensitivity refits the primary model only.
- Family, stated before correction: BH q within each dataset across cell types, separately for the primary model and the group comparison.
- The shared rules additionally require an unadjusted sensitivity without batch covariates. It retains age_per_decade and sex_male.

## Fixed implementation details

- Inputs are the Gao and Yazar h5ad files and donor tables specified for Block 2. Yazar excludes rows labeled peripheral blood mononuclear cell. Donor is the unit for every model, permutation, and split-half summary.
- The human technical gene mask excludes uppercase MT-, RPL, RPS, MRPL, MRPS, XIST, TSIX, DDX3Y, EIF2S3Y, UTY, and KDM5D before depth calculation and normalization. ERCC controls are also excluded, consistent with the existing donor tables.
- Full-depth and rarefied analyses use cells with at least 1500 post-mask raw counts. The 1500 threshold is an eligibility filter for the full-depth primary, not downsampling. A donor-cell-type unit requires at least 20 eligible cells. A slope requires at least 6 eligible donors, and a group comparison requires at least 6 eligible donors in each of MATURE and OLD. Units below these floors are UNTESTABLE.
- Age bins: MATURE 30 to 39, MID 40 to 59, OLD 60 and older. Donors under 30 are excluded from all Block 2 calculations, including detected-gene-set derivation, fixed-pair sampling, GCL scoring, and every model. LT20 donors do not enter continuous-age models in this rerun.
- The original Block 2 count cache is reused for the projected per-cell counts: its values are donor-local (per-cell projection and UMI filtering only). Detected-gene sets are derived from the cached cells of mature donors only, fixed pairs are sampled from those mature gene sets, and checkpoints are kept in a separate mature cache directory.
- For a cell type with at most 2000 detected genes, all unique unordered gene pairs are used. Above 2000 genes, exactly 20000 unique unordered pairs are sampled without replacement using base seed 20261117 and stable dataset and cell-type offsets. The same pairs are used for every donor, half, and rarefied fit within that cell type.
- Spearman correlations use average ranks for ties. Constant gene vectors within a donor or half yield undefined pairwise correlations and those pairs are omitted for that donor or half. GCL is the arithmetic mean of the remaining fixed-pair correlations. A GCL with no defined pairs is UNTESTABLE.
- The split-half noise floor is the mean over donors and 5 seeds of the absolute GCL difference between the two within-donor halves. Split-half seeds are 20261116, 20261117, 20261118, 20261119, and 20261120. Each half uses the cell type's fixed gene pairs.
- Group comparison uses OLD minus MATURE. Batch residualization uses intercept plus batch fixed effects within eligible MATURE and OLD donors. Its CI and MDE80 use the OLS group coefficient and standard error; its p-value is the two-sided donor-permutation p-value. Permutations use base seed 20261118 and 2000 permutations.
- Rarefaction uses multinomial downsampling to exactly 1500 post-mask counts per eligible cell with base seed 20261119. It retains the primary gene set, fixed pairs, and eligible donors.
- MDE80 is (t_0.975 + t_0.80) times the effect standard error using model residual degrees of freedom. For a null primary, MDE80 at or below the split-half noise floor is small enough to exclude a change detectable above technical noise; MDE80 above that floor is large.
- Depth adjustment adds the donor table's mean_log_umi_ge1500 for the same donor-cell-type unit.
- Primary and group BH families include only evaluated cell-type rows within each dataset. Depth-adjusted, unadjusted, and rarefied sensitivities are reported with BH q within their own dataset-by-model slices for transparency; the two required confirmatory families remain primary and group comparison.
- This block does not score P1, so no P1/P1m change applies here.

## Pre-stated verdict rules

For each cell type, PROVES that coordination changes with age requires batch-adjusted primary q < 0.05, absolute primary slope greater than the split-half noise floor, and the same sign in the group comparison. DOES NOT PROVE applies to a null primary whose MDE80 is at or below the split-half noise floor. INCONCLUSIVE applies to a null primary with MDE80 above the noise floor or to disagreement in direction between primary and group comparison. The direction, coordination rises or falls, is reported with every verdict. Levy 2020 predicts a fall, so a rise is a different finding.

For each dataset summary, PROVES requires at least one PROVES cell type and no INCONCLUSIVE cell type with a significant primary effect failing a verdict condition. DOES NOT PROVE requires every evaluated cell type to meet DOES NOT PROVE. All other summaries are INCONCLUSIVE. Gao and Yazar are reported separately and are never pooled. The single commit-label verdict is PROVES only if both datasets are PROVES, DOES NOT PROVE only if both are DOES NOT PROVE, and otherwise INCONCLUSIVE.

## Comparison against the original Block 2

Original verdicts (full range, OLD vs YOUNG): Gao INCONCLUSIVE (0 PROVES, 9 DOES NOT PROVE, 4 INCONCLUSIVE, 0 UNTESTABLE; inconclusive Adipocyte, Cardiomyocyte, Epicardial, Mast); Yazar INCONCLUSIVE (0 PROVES, 11 DOES NOT PROVE, 5 INCONCLUSIVE, 12 UNTESTABLE; inconclusive CD4-positive, alpha-beta cytotoxic T cell, memory B cell, mucosal-associated invariant T cell, natural killer cell, regulatory T cell). Yazar natural killer cell was the closest call (primary q 0.03493, group q 0.01399, MDE80 0.0002367 vs noise floor 0.004101, INCONCLUSIVE because the slope is below the floor). The mature comparison table in the Results section records the mature NK numbers explicitly and, per cell type, whether each verdict survives or dies.

## Execution note

Delegated-model note: Kimi hit its 5-hour usage quota at about 13:13 PST Nov 21 and Antigravity Opus 5.5 High hit its individual quota (reset in 139h59m), then Antigravity Gemini 4 Neon XHigh hit its quota (reset in about 4h31m from 21:41 UTC) after spending its run on a multi-round plan review with no implementation. No delegated model produced any Block 1 to 3 output. This block was therefore executed directly by the runner from the frozen original script by copy-and-remap (mature-only detected-gene sets and pairs, mature bin relabeling, mature output paths, separate mature checkpoint cache), with the coordinator to verify before reporting. Model for this block: direct runner execution, no delegated model.
## Results

All effects are changes in mean pairwise Spearman correlation. Primary slopes are per decade. Group effects are OLD minus MATURE after batch residualization. NA denotes an UNTESTABLE or not applicable field.

| Dataset | Cell type | N genes | N pairs | N donors | Primary slope | Primary q | Noise floor | MDE80 | OLD-MATURE | Group q | Verdict | Direction |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gao | Adipocyte | 5974 | 20000 | 35 | -0.00101 | 0.7304 | 0.006894 | 0.003383 | NA | NA | INCONCLUSIVE | coordination falls |
| gao | Cardiomyocyte | 6396 | 20000 | 102 | -0.0001323 | 0.9135 | 0.001743 | 0.00344 | 0.001936 | 0.7236 | INCONCLUSIVE | coordination falls |
| gao | Endocardial | 5079 | 20000 | 66 | -0.002262 | 0.1012 | 0.004782 | 0.002586 | -0.00514 | 0.08096 | DOES NOT PROVE | coordination falls |
| gao | Endothelial | 4266 | 20000 | 80 | 0.0003921 | 0.7788 | 0.004953 | 0.002109 | -0.0009832 | 0.7236 | INCONCLUSIVE | coordination rises |
| gao | Epicardial | 5438 | 20000 | 18 | -0.0002482 | 0.8302 | 0.005369 | 0.002286 | NA | NA | INCONCLUSIVE | coordination falls |
| gao | Fibroblast | 5151 | 20000 | 99 | -0.001914 | 0.1012 | 0.001925 | 0.002078 | -0.006633 | 0.08096 | INCONCLUSIVE | coordination falls |
| gao | LEC | 4901 | 20000 | 34 | -0.0008688 | 0.7304 | 0.006314 | 0.002613 | NA | NA | INCONCLUSIVE | coordination falls |
| gao | Lymphoid | 3998 | 20000 | 56 | 0.0001378 | 0.8302 | 0.00419 | 0.001318 | 0.001005 | 0.7236 | DOES NOT PROVE | coordination rises |
| gao | Mast | 3801 | 20000 | 17 | 0.0006041 | 0.7304 | 0.003372 | 0.002027 | NA | NA | INCONCLUSIVE | coordination rises |
| gao | Myeloid | 4787 | 20000 | 99 | -0.000979 | 0.1868 | 0.002155 | 0.001352 | -0.004511 | 0.08096 | DOES NOT PROVE | coordination falls |
| gao | Neuronal | 4535 | 20000 | 56 | -0.001409 | 0.301 | 0.007334 | 0.002349 | -0.005393 | 0.1552 | DOES NOT PROVE | coordination falls |
| gao | Pericyte | 3946 | 20000 | 88 | -0.0001705 | 0.7788 | 0.002326 | 0.0009018 | -0.001041 | 0.5757 | DOES NOT PROVE | coordination falls |
| gao | vSMC | 4066 | 20000 | 79 | -0.0003695 | 0.7788 | 0.006079 | 0.001821 | -0.001167 | 0.5757 | DOES NOT PROVE | coordination falls |
| yazar | CD14-low, CD16-positive monocyte | 3911 | 20000 | 257 | 0.0003671 | 0.5911 | 0.005132 | 0.0006122 | 0.0002692 | 0.8096 | DOES NOT PROVE | coordination rises |
| yazar | CD14-positive monocyte | 3345 | 20000 | 408 | 0.000119 | 0.638 | 0.005376 | 0.0004236 | 0.0004524 | 0.7756 | DOES NOT PROVE | coordination rises |
| yazar | CD16-negative, CD56-bright natural killer cell, human | 3210 | 20000 | 5 | NA | NA | NA | NA | NA | NA | UNTESTABLE | UNTESTABLE |
| yazar | CD4-positive, alpha-beta T cell | 7045 | 20000 | 0 | NA | NA | NA | NA | NA | NA | UNTESTABLE | UNTESTABLE |
| yazar | CD4-positive, alpha-beta cytotoxic T cell | 2927 | 20000 | 92 | -0.0004581 | 0.638 | 0.005397 | 0.001678 | NA | NA | INCONCLUSIVE | coordination falls |
| yazar | CD8-positive, alpha-beta T cell | 5540 | 20000 | 0 | NA | NA | NA | NA | NA | NA | UNTESTABLE | UNTESTABLE |
| yazar | central memory CD4-positive, alpha-beta T cell | 3066 | 20000 | 909 | -1.839e-05 | 0.851 | 0.00331 | 0.0002013 | -0.0007614 | 0.1469 | DOES NOT PROVE | coordination falls |
| yazar | central memory CD8-positive, alpha-beta T cell | 3147 | 20000 | 106 | -0.0006003 | 0.5911 | 0.005495 | 0.001274 | NA | NA | INCONCLUSIVE | coordination falls |
| yazar | conventional dendritic cell | 5348 | 20000 | 5 | NA | NA | NA | NA | NA | NA | UNTESTABLE | UNTESTABLE |
| yazar | dendritic cell | 4970 | 20000 | 0 | NA | NA | NA | NA | NA | NA | UNTESTABLE | UNTESTABLE |
| yazar | double negative thymocyte | 4381 | 20000 | 0 | NA | NA | NA | NA | NA | NA | UNTESTABLE | UNTESTABLE |
| yazar | effector memory CD4-positive, alpha-beta T cell | 3017 | 20000 | 258 | 0.0002743 | 0.638 | 0.005361 | 0.0007159 | -0.001148 | 0.7756 | INCONCLUSIVE | coordination rises |
| yazar | effector memory CD8-positive, alpha-beta T cell | 2849 | 20000 | 722 | 2.918e-05 | 0.851 | 0.003738 | 0.0002662 | -0.0002787 | 0.7756 | INCONCLUSIVE | coordination rises |
| yazar | erythrocyte | 2209 | 20000 | 0 | NA | NA | NA | NA | NA | NA | UNTESTABLE | UNTESTABLE |
| yazar | gamma-delta T cell | 2906 | 20000 | 120 | 0.0003522 | 0.5911 | 0.004465 | 0.0006936 | -0.0002385 | 0.8096 | INCONCLUSIVE | coordination rises |
| yazar | hematopoietic precursor cell | 5575 | 20000 | 0 | NA | NA | NA | NA | NA | NA | UNTESTABLE | UNTESTABLE |
| yazar | innate lymphoid cell | 3225 | 20000 | 0 | NA | NA | NA | NA | NA | NA | UNTESTABLE | UNTESTABLE |
| yazar | memory B cell | 3080 | 20000 | 238 | 0.0003369 | 0.7318 | 0.01053 | 0.001601 | -0.001658 | 0.7756 | INCONCLUSIVE | coordination rises |
| yazar | mucosal-associated invariant T cell | 3015 | 20000 | 38 | 0.0006387 | 0.638 | 0.003854 | 0.002108 | NA | NA | INCONCLUSIVE | coordination rises |
| yazar | naive B cell | 2637 | 20000 | 312 | 4.249e-06 | 0.9825 | 0.004394 | 0.0005443 | 0.0002289 | 0.8096 | DOES NOT PROVE | coordination rises |
| yazar | naive thymus-derived CD4-positive, alpha-beta T cell | 2726 | 20000 | 778 | 3.204e-05 | 0.7318 | 0.002726 | 0.0001688 | -0.0002024 | 0.7756 | INCONCLUSIVE | coordination rises |
| yazar | naive thymus-derived CD8-positive, alpha-beta T cell | 2932 | 20000 | 352 | 0.0002268 | 0.5911 | 0.004559 | 0.0004319 | 0.0004368 | 0.7756 | DOES NOT PROVE | coordination rises |
| yazar | natural killer cell | 3035 | 20000 | 784 | -0.0003086 | 0.04182 | 0.004052 | 0.0002866 | -0.002035 | 0.01199 | INCONCLUSIVE | coordination falls |
| yazar | plasmablast | 5837 | 20000 | 14 | NA | NA | 0.02484 | NA | NA | NA | UNTESTABLE | UNTESTABLE |
| yazar | plasmacytoid dendritic cell | 4247 | 20000 | 3 | NA | NA | NA | NA | NA | NA | UNTESTABLE | UNTESTABLE |
| yazar | platelet | 3067 | 20000 | 0 | NA | NA | NA | NA | NA | NA | UNTESTABLE | UNTESTABLE |
| yazar | regulatory T cell | 2847 | 20000 | 69 | 0.0007791 | 0.638 | 0.006734 | 0.002449 | NA | NA | INCONCLUSIVE | coordination rises |
| yazar | transitional stage B cell | 3047 | 20000 | 195 | 0.0003026 | 0.638 | 0.006811 | 0.001021 | 0.001495 | 0.7756 | DOES NOT PROVE | coordination rises |

### Sensitivity models

| Dataset | Cell type | Depth adjusted slope | Depth q | Unadjusted slope | Unadjusted q | 1500 rarefied slope | 1500 q |
| --- | --- | --- | --- | --- | --- | --- | --- |
| gao | Adipocyte | 0.0002226 | 0.9987 | -0.001347 | 0.544 | -0.0002824 | 0.619 |
| gao | Cardiomyocyte | 0.0005119 | 0.878 | -0.0005819 | 0.8888 | 1.265e-05 | 0.7667 |
| gao | Endocardial | -0.001467 | 0.2862 | -0.002578 | 0.0527 | -7.839e-05 | 0.7667 |
| gao | Endothelial | -1.144e-06 | 0.9987 | -2.673e-05 | 0.9744 | 2.639e-05 | 0.8398 |
| gao | Epicardial | -0.0005623 | 0.8258 | -0.0008115 | 0.5453 | -0.0005227 | 0.02615 |
| gao | Fibroblast | -0.001341 | 0.2862 | -0.001962 | 0.1078 | -1.993e-05 | 0.7667 |
| gao | LEC | -0.0005262 | 0.8258 | -0.0006677 | 0.6811 | 0.0002957 | 0.619 |
| gao | Lymphoid | 5.473e-05 | 0.9987 | -0.0001362 | 0.8925 | -8.768e-05 | 0.7667 |
| gao | Mast | -3.639e-05 | 0.9987 | -0.0001127 | 0.8925 | 0.0006061 | 0.7667 |
| gao | Myeloid | -0.0005706 | 0.2885 | -0.001426 | 0.0527 | 5.024e-05 | 0.7667 |
| gao | Neuronal | -0.001341 | 0.2862 | -0.001757 | 0.1248 | 3.469e-05 | 0.8398 |
| gao | Pericyte | -0.0003847 | 0.4504 | -0.0003777 | 0.544 | 2.88e-05 | 0.7667 |
| gao | vSMC | -0.0005768 | 0.6759 | -0.0007138 | 0.544 | 0.0001047 | 0.7667 |
| yazar | CD14-low, CD16-positive monocyte | 0.0003667 | 0.508 | 0.0003274 | 0.4845 | 3.749e-05 | 0.9896 |
| yazar | CD14-positive monocyte | 0.0001646 | 0.6778 | 0.0003195 | 0.1553 | -9.06e-05 | 0.511 |
| yazar | CD16-negative, CD56-bright natural killer cell, human | NA | NA | NA | NA | NA | NA |
| yazar | CD4-positive, alpha-beta T cell | NA | NA | NA | NA | NA | NA |
| yazar | CD4-positive, alpha-beta cytotoxic T cell | -0.0001177 | 0.8982 | 0.0002644 | 0.7634 | -0.0001826 | 0.7781 |
| yazar | CD8-positive, alpha-beta T cell | NA | NA | NA | NA | NA | NA |
| yazar | central memory CD4-positive, alpha-beta T cell | 2.126e-05 | 0.8548 | 0.0001235 | 0.3164 | -7.102e-07 | 0.9896 |
| yazar | central memory CD8-positive, alpha-beta T cell | -0.0006243 | 0.508 | -0.0002262 | 0.7488 | -5.33e-05 | 0.9896 |
| yazar | conventional dendritic cell | NA | NA | NA | NA | NA | NA |
| yazar | dendritic cell | NA | NA | NA | NA | NA | NA |
| yazar | double negative thymocyte | NA | NA | NA | NA | NA | NA |
| yazar | effector memory CD4-positive, alpha-beta T cell | 0.0002496 | 0.6778 | 7.766e-05 | 0.8564 | 0.0003285 | 0.04175 |
| yazar | effector memory CD8-positive, alpha-beta T cell | 4.496e-05 | 0.8548 | 3.046e-05 | 0.8564 | 4.747e-05 | 0.511 |
| yazar | erythrocyte | NA | NA | NA | NA | NA | NA |
| yazar | gamma-delta T cell | 0.0003562 | 0.508 | 0.0002401 | 0.4845 | -1.645e-06 | 0.9896 |
| yazar | hematopoietic precursor cell | NA | NA | NA | NA | NA | NA |
| yazar | innate lymphoid cell | NA | NA | NA | NA | NA | NA |
| yazar | memory B cell | -0.0001568 | 0.8548 | 0.0002737 | 0.7634 | 2.967e-05 | 0.9896 |
| yazar | mucosal-associated invariant T cell | 0.0002643 | 0.8548 | 0.0001881 | 0.7634 | -0.0005946 | 0.5195 |
| yazar | naive B cell | 1.275e-05 | 0.9419 | 1.02e-05 | 0.9488 | 5.428e-05 | 0.7781 |
| yazar | naive thymus-derived CD4-positive, alpha-beta T cell | 6.648e-05 | 0.6778 | 1.442e-05 | 0.8569 | 5.139e-05 | 0.511 |
| yazar | naive thymus-derived CD8-positive, alpha-beta T cell | 0.0001387 | 0.6836 | 0.0001587 | 0.4845 | 0.0001037 | 0.511 |
| yazar | natural killer cell | -0.0002253 | 0.2501 | -0.0003045 | 0.01093 | -2.06e-05 | 0.8291 |
| yazar | plasmablast | NA | NA | -0.001006 | 0.8569 | NA | NA |
| yazar | plasmacytoid dendritic cell | NA | NA | NA | NA | NA | NA |
| yazar | platelet | NA | NA | NA | NA | NA | NA |
| yazar | regulatory T cell | 0.000674 | 0.7181 | 0.000449 | 0.7488 | 0.0003857 | 0.511 |
| yazar | transitional stage B cell | 0.0002477 | 0.7181 | 0.000215 | 0.7488 | 1.863e-05 | 0.9896 |

### Verdicts

- Gao: INCONCLUSIVE; no proven direction. PROVES 0, DOES NOT PROVE 6, INCONCLUSIVE 7, UNTESTABLE 0 cell types.
  Inconclusive units: Adipocyte, Cardiomyocyte, Endothelial, Epicardial, Fibroblast, LEC, Mast.
- Yazar: INCONCLUSIVE; no proven direction. PROVES 0, DOES NOT PROVE 6, INCONCLUSIVE 10, UNTESTABLE 12 cell types.
  Inconclusive units: CD4-positive, alpha-beta cytotoxic T cell, central memory CD8-positive, alpha-beta T cell, effector memory CD4-positive, alpha-beta T cell, effector memory CD8-positive, alpha-beta T cell, gamma-delta T cell, memory B cell, mucosal-associated invariant T cell, naive thymus-derived CD4-positive, alpha-beta T cell, natural killer cell, regulatory T cell.
  Untestable units: CD16-negative, CD56-bright natural killer cell, human, CD4-positive, alpha-beta T cell, CD8-positive, alpha-beta T cell, conventional dendritic cell, dendritic cell, double negative thymocyte, erythrocyte, hematopoietic precursor cell, innate lymphoid cell, plasmablast, plasmacytoid dendritic cell, platelet.
- Commit-label summary: INCONCLUSIVE. This does not pool donor data or replace the separate Gao and Yazar verdicts.
- A rise would differ from the Levy 2020 prediction of falling coordination. No significant rise is described as confirmation of that prediction.

### Mature natural killer cell (Yazar)

Yazar natural killer cell was the closest call in the original Block 2 (primary q 0.03493, group q 0.01399, MDE80 0.0002367 vs noise floor 0.004101, INCONCLUSIVE because the slope was below the floor). The mature numbers are: 3035 detected genes from mature donors only, 20000 fixed pairs, 784 donors in the primary model.
Primary slope -0.0003086 per decade, 95 percent CI -0.0005091 to -0.000108, p 0.002614, BH q 0.04182, split-half noise floor 0.004052, MDE80 0.0002866.
Group comparison OLD minus MATURE -0.002035, permutation p 0.0009995, BH q 0.01199, 38 MATURE and 598 OLD donors, batch-residualized means 0.008106 MATURE and 0.006071 OLD.
Sensitivities: depth-adjusted slope -0.0002253 q 0.2501; unadjusted slope -0.0003045 q 0.01093; 1500 rarefied slope -2.06e-05 q 0.8291.
Mature verdict: INCONCLUSIVE; coordination falls. Original verdict: INCONCLUSIVE with primary slope -0.0002595 and group difference -0.001437.

### Comparison against the original Block 2

Per cell type, the original verdict (full age range, OLD vs YOUNG) against the mature verdict (30-plus, OLD vs MATURE). SURVIVED means the verdict is identical in the mature rerun. DIED means the verdict changed. UNTESTABLE means the cell type is untestable in the original, in the mature rerun, or in both. Counts: SURVIVED 21, DIED 8, UNTESTABLE 12.

| Dataset | Cell type | Original verdict | Mature verdict | Comparison | Note |
| --- | --- | --- | --- | --- | --- | --- |
| gao | Adipocyte | INCONCLUSIVE | INCONCLUSIVE | SURVIVED | same verdict |
| gao | Cardiomyocyte | INCONCLUSIVE | INCONCLUSIVE | SURVIVED | same verdict |
| gao | Endocardial | DOES NOT PROVE | DOES NOT PROVE | SURVIVED | same verdict |
| gao | Endothelial | DOES NOT PROVE | INCONCLUSIVE | DIED | verdict changed from DOES NOT PROVE to INCONCLUSIVE |
| gao | Epicardial | INCONCLUSIVE | INCONCLUSIVE | SURVIVED | same verdict |
| gao | Fibroblast | DOES NOT PROVE | INCONCLUSIVE | DIED | verdict changed from DOES NOT PROVE to INCONCLUSIVE |
| gao | LEC | DOES NOT PROVE | INCONCLUSIVE | DIED | verdict changed from DOES NOT PROVE to INCONCLUSIVE |
| gao | Lymphoid | DOES NOT PROVE | DOES NOT PROVE | SURVIVED | same verdict |
| gao | Mast | INCONCLUSIVE | INCONCLUSIVE | SURVIVED | same verdict |
| gao | Myeloid | DOES NOT PROVE | DOES NOT PROVE | SURVIVED | same verdict |
| gao | Neuronal | DOES NOT PROVE | DOES NOT PROVE | SURVIVED | same verdict |
| gao | Pericyte | DOES NOT PROVE | DOES NOT PROVE | SURVIVED | same verdict |
| gao | vSMC | DOES NOT PROVE | DOES NOT PROVE | SURVIVED | same verdict |
| yazar | CD14-low, CD16-positive monocyte | DOES NOT PROVE | DOES NOT PROVE | SURVIVED | same verdict |
| yazar | CD14-positive monocyte | DOES NOT PROVE | DOES NOT PROVE | SURVIVED | same verdict |
| yazar | CD16-negative, CD56-bright natural killer cell, human | UNTESTABLE | UNTESTABLE | UNTESTABLE | untestable in both |
| yazar | CD4-positive, alpha-beta T cell | UNTESTABLE | UNTESTABLE | UNTESTABLE | untestable in both |
| yazar | CD4-positive, alpha-beta cytotoxic T cell | INCONCLUSIVE | INCONCLUSIVE | SURVIVED | same verdict |
| yazar | CD8-positive, alpha-beta T cell | UNTESTABLE | UNTESTABLE | UNTESTABLE | untestable in both |
| yazar | central memory CD4-positive, alpha-beta T cell | DOES NOT PROVE | DOES NOT PROVE | SURVIVED | same verdict |
| yazar | central memory CD8-positive, alpha-beta T cell | DOES NOT PROVE | INCONCLUSIVE | DIED | verdict changed from DOES NOT PROVE to INCONCLUSIVE |
| yazar | conventional dendritic cell | UNTESTABLE | UNTESTABLE | UNTESTABLE | untestable in both |
| yazar | dendritic cell | UNTESTABLE | UNTESTABLE | UNTESTABLE | untestable in both |
| yazar | double negative thymocyte | UNTESTABLE | UNTESTABLE | UNTESTABLE | untestable in both |
| yazar | effector memory CD4-positive, alpha-beta T cell | DOES NOT PROVE | INCONCLUSIVE | DIED | verdict changed from DOES NOT PROVE to INCONCLUSIVE |
| yazar | effector memory CD8-positive, alpha-beta T cell | DOES NOT PROVE | INCONCLUSIVE | DIED | verdict changed from DOES NOT PROVE to INCONCLUSIVE |
| yazar | erythrocyte | UNTESTABLE | UNTESTABLE | UNTESTABLE | untestable in both |
| yazar | gamma-delta T cell | DOES NOT PROVE | INCONCLUSIVE | DIED | verdict changed from DOES NOT PROVE to INCONCLUSIVE |
| yazar | hematopoietic precursor cell | UNTESTABLE | UNTESTABLE | UNTESTABLE | untestable in both |
| yazar | innate lymphoid cell | UNTESTABLE | UNTESTABLE | UNTESTABLE | untestable in both |
| yazar | memory B cell | INCONCLUSIVE | INCONCLUSIVE | SURVIVED | same verdict |
| yazar | mucosal-associated invariant T cell | INCONCLUSIVE | INCONCLUSIVE | SURVIVED | same verdict |
| yazar | naive B cell | DOES NOT PROVE | DOES NOT PROVE | SURVIVED | same verdict |
| yazar | naive thymus-derived CD4-positive, alpha-beta T cell | DOES NOT PROVE | INCONCLUSIVE | DIED | verdict changed from DOES NOT PROVE to INCONCLUSIVE |
| yazar | naive thymus-derived CD8-positive, alpha-beta T cell | DOES NOT PROVE | DOES NOT PROVE | SURVIVED | same verdict |
| yazar | natural killer cell | INCONCLUSIVE | INCONCLUSIVE | SURVIVED | same verdict |
| yazar | plasmablast | UNTESTABLE | UNTESTABLE | UNTESTABLE | untestable in both |
| yazar | plasmacytoid dendritic cell | UNTESTABLE | UNTESTABLE | UNTESTABLE | untestable in both |
| yazar | platelet | UNTESTABLE | UNTESTABLE | UNTESTABLE | untestable in both |
| yazar | regulatory T cell | INCONCLUSIVE | INCONCLUSIVE | SURVIVED | same verdict |
| yazar | transitional stage B cell | DOES NOT PROVE | DOES NOT PROVE | SURVIVED | same verdict |
### Reproducibility

Pair base seed 20261117, permutation base seed 20261118, rarefaction base seed 20261119, and split-half seeds 20261116, 20261117, 20261118, 20261119, 20261120. Exact derived per-cell-type pair and permutation seeds are in the result CSV and analysis checkpoints.
Primary and group comparison family sizes: gao primary 13, group 9, yazar primary 16, group 12.
No primary cells were downsampled. The only downsampling is the explicitly labeled 1500 rarefaction sensitivity.
