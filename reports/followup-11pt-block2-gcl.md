# Block 2: Gene coordination, young vs old, donor level

Status: complete. This specification was written before running the analysis.

## Exact analysis specification

- Gene set, fixed before running: the detected-gene set per cell type, defined as genes detected in at least 10 percent of cells in that cell type pooled across donors. Record N genes per cell type. Do not change the gene set after seeing results.
- Per donor per cell type, using the shared eligibility floors: compute the gene coordination level (GCL) as the mean pairwise Spearman correlation between genes across that donor's cells, using full-depth within-cell proportions. For cell types with more than 2000 genes in the set, subsample gene PAIRS, not genes, with a fixed seed so the estimator is identical across donors. Record the pair count and seed.
- Primary model, per dataset and cell type: OLS of per-donor GCL on age_per_decade + sex_male + batch, with Gao study indicators and Yazar pool fixed effects. Report slope per decade, 95 percent CI, p-value, N donors, and MDE80.
- Group comparison secondary: mean GCL in OLD vs YOUNG donors per cell type, donor-permutation null with 2000 permutations and a fixed seed, using batch-residualized GCL after residualizing on batch first and then permuting donor age-bin labels.
- Split-half technical null: split each donor's cells into two halves with 5 fixed seeds, compute GCL per half, and report the half-to-half noise floor per cell type. A claimed age effect must exceed this floor. Report the floor next to the slope.
- Depth sensitivity: refit the primary model adding per-donor mean log UMI as a covariate.
- A labeled 1500 rarefaction sensitivity refits the primary model only.
- Family, stated before correction: BH q within each dataset across cell types, separately for the primary model and the group comparison.
- The shared rules additionally require an unadjusted sensitivity without batch covariates. It retains age_per_decade and sex_male.

## Fixed implementation details

- Inputs are the Gao and Yazar h5ad files and donor tables specified for Block 2. Yazar excludes rows labeled peripheral blood mononuclear cell. Donor is the unit for every model, permutation, and split-half summary.
- The human technical gene mask excludes uppercase MT-, RPL, RPS, MRPL, MRPS, XIST, TSIX, DDX3Y, EIF2S3Y, UTY, and KDM5D before depth calculation and normalization. ERCC controls are also excluded, consistent with the existing donor tables.
- Full-depth and rarefied analyses use cells with at least 1500 post-mask raw counts. The 1500 threshold is an eligibility filter for the full-depth primary, not downsampling. A donor-cell-type unit requires at least 20 eligible cells. A slope requires at least 6 eligible donors, and a group comparison requires at least 6 eligible donors in each of YOUNG and OLD. Units below these floors are UNTESTABLE.
- Age bins reuse assign_age_bin from scripts/run_followup_between_subject_gao_yazar.py: LT20 under 20, YOUNG 20 to 39, MID 40 to 59, OLD 60 and older. LT20 donors remain in continuous-age models but not the OLD vs YOUNG comparison.
- For a cell type with at most 2000 detected genes, all unique unordered gene pairs are used. Above 2000 genes, exactly 20000 unique unordered pairs are sampled without replacement using base seed 20261117 and stable dataset and cell-type offsets. The same pairs are used for every donor, half, and rarefied fit within that cell type.
- Spearman correlations use average ranks for ties. Constant gene vectors within a donor or half yield undefined pairwise correlations and those pairs are omitted for that donor or half. GCL is the arithmetic mean of the remaining fixed-pair correlations. A GCL with no defined pairs is UNTESTABLE.
- The split-half noise floor is the mean over donors and 5 seeds of the absolute GCL difference between the two within-donor halves. Split-half seeds are 20261116, 20261117, 20261118, 20261119, and 20261120. Each half uses the cell type's fixed gene pairs.
- Group comparison uses OLD minus YOUNG. Batch residualization uses intercept plus batch fixed effects within eligible YOUNG and OLD donors. Its CI and MDE80 use the OLS group coefficient and standard error; its p-value is the two-sided donor-permutation p-value. Permutations use base seed 20261118 and 2000 permutations.
- Rarefaction uses multinomial downsampling to exactly 1500 post-mask counts per eligible cell with base seed 20261119. It retains the primary gene set, fixed pairs, and eligible donors.
- MDE80 is (t_0.975 + t_0.80) times the effect standard error using model residual degrees of freedom. For a null primary, MDE80 at or below the split-half noise floor is small enough to exclude a change detectable above technical noise; MDE80 above that floor is large.
- Depth adjustment adds the donor table's mean_log_umi_ge1500 for the same donor-cell-type unit.
- Primary and group BH families include only evaluated cell-type rows within each dataset. Depth-adjusted, unadjusted, and rarefied sensitivities are reported with BH q within their own dataset-by-model slices for transparency; the two required confirmatory families remain primary and group comparison.

## Pre-stated verdict rules

For each cell type, PROVES that coordination changes with age requires batch-adjusted primary q < 0.05, absolute primary slope greater than the split-half noise floor, and the same sign in the group comparison. DOES NOT PROVE applies to a null primary whose MDE80 is at or below the split-half noise floor. INCONCLUSIVE applies to a null primary with MDE80 above the noise floor or to disagreement in direction between primary and group comparison. The direction, coordination rises or falls, is reported with every verdict. Levy 2020 predicts a fall, so a rise is a different finding.

For each dataset summary, PROVES requires at least one PROVES cell type and no INCONCLUSIVE cell type with a significant primary effect failing a verdict condition. DOES NOT PROVE requires every evaluated cell type to meet DOES NOT PROVE. All other summaries are INCONCLUSIVE. Gao and Yazar are reported separately and are never pooled. The single commit-label verdict is PROVES only if both datasets are PROVES, DOES NOT PROVE only if both are DOES NOT PROVE, and otherwise INCONCLUSIVE.

## Results

All effects are changes in mean pairwise Spearman correlation. Primary slopes are per decade. Group effects are OLD minus YOUNG after batch residualization. NA denotes an UNTESTABLE or not applicable field.

| Dataset | Cell type | N genes | N pairs | N donors | Primary slope | Primary q | Noise floor | MDE80 | OLD-YOUNG | Group q | Verdict | Direction |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gao | Adipocyte | 5974 | 20000 | 37 | -0.0002986 | 0.7891 | 0.006752 | 0.002801 | NA | NA | INCONCLUSIVE | coordination falls |
| gao | Cardiomyocyte | 6406 | 20000 | 109 | -0.0003919 | 0.7891 | 0.001707 | 0.002693 | 0.0001989 | 0.949 | INCONCLUSIVE | coordination falls |
| gao | Endocardial | 5082 | 20000 | 72 | -0.001333 | 0.3347 | 0.004604 | 0.001965 | -0.003785 | 0.2274 | DOES NOT PROVE | coordination falls |
| gao | Endothelial | 4259 | 20000 | 87 | 0.000369 | 0.7444 | 0.004995 | 0.001602 | 0.0003491 | 0.949 | DOES NOT PROVE | coordination rises |
| gao | Epicardial | 5448 | 20000 | 20 | -0.00138 | 0.3347 | 0.005579 | 0.002165 | NA | NA | INCONCLUSIVE | coordination falls |
| gao | Fibroblast | 5139 | 20000 | 106 | -0.0009443 | 0.3347 | 0.001888 | 0.001676 | -0.003883 | 0.2274 | DOES NOT PROVE | coordination falls |
| gao | LEC | 4946 | 20000 | 36 | -0.0009088 | 0.5097 | 0.006399 | 0.002172 | -0.005667 | 0.2274 | DOES NOT PROVE | coordination falls |
| gao | Lymphoid | 3977 | 20000 | 59 | 0.0003368 | 0.624 | 0.004347 | 0.001095 | 0.001403 | 0.5335 | DOES NOT PROVE | coordination rises |
| gao | Mast | 3781 | 20000 | 18 | 0.0001423 | 0.7891 | 0.002677 | 0.001005 | NA | NA | INCONCLUSIVE | coordination rises |
| gao | Myeloid | 4776 | 20000 | 106 | -0.0004075 | 0.5331 | 0.002132 | 0.001077 | -0.0024 | 0.2432 | DOES NOT PROVE | coordination falls |
| gao | Neuronal | 4486 | 20000 | 60 | -0.001118 | 0.3347 | 0.006919 | 0.001989 | -0.004176 | 0.2409 | DOES NOT PROVE | coordination falls |
| gao | Pericyte | 3922 | 20000 | 95 | -6.334e-05 | 0.7891 | 0.002292 | 0.0006689 | -0.0007255 | 0.5335 | DOES NOT PROVE | coordination falls |
| gao | vSMC | 4056 | 20000 | 86 | -0.0008447 | 0.3347 | 0.006616 | 0.001561 | -0.003114 | 0.2274 | DOES NOT PROVE | coordination falls |
| yazar | CD14-low, CD16-positive monocyte | 3903 | 20000 | 263 | 0.000241 | 0.5992 | 0.005223 | 0.000505 | 0.0003491 | 0.9501 | DOES NOT PROVE | coordination rises |
| yazar | CD14-positive monocyte | 3344 | 20000 | 429 | 0.0001455 | 0.5992 | 0.005371 | 0.0003423 | 0.0007794 | 0.6973 | DOES NOT PROVE | coordination rises |
| yazar | CD16-negative, CD56-bright natural killer cell, human | 3206 | 20000 | 6 | NA | NA | 0.006201 | NA | NA | NA | UNTESTABLE | UNTESTABLE |
| yazar | CD4-positive, alpha-beta T cell | 7090 | 20000 | 0 | NA | NA | NA | NA | NA | NA | UNTESTABLE | UNTESTABLE |
| yazar | CD4-positive, alpha-beta cytotoxic T cell | 2924 | 20000 | 93 | -0.0006551 | 0.5992 | 0.005333 | 0.001387 | NA | NA | INCONCLUSIVE | coordination falls |
| yazar | CD8-positive, alpha-beta T cell | 5681 | 20000 | 0 | NA | NA | NA | NA | NA | NA | UNTESTABLE | UNTESTABLE |
| yazar | central memory CD4-positive, alpha-beta T cell | 3064 | 20000 | 955 | -5.85e-05 | 0.5992 | 0.003333 | 0.0001717 | -0.0005976 | 0.2204 | DOES NOT PROVE | coordination falls |
| yazar | central memory CD8-positive, alpha-beta T cell | 3138 | 20000 | 111 | -0.0003465 | 0.5992 | 0.005287 | 0.000935 | -0.0001817 | 0.967 | DOES NOT PROVE | coordination falls |
| yazar | conventional dendritic cell | 5338 | 20000 | 5 | NA | NA | NA | NA | NA | NA | UNTESTABLE | UNTESTABLE |
| yazar | dendritic cell | 4874 | 20000 | 0 | NA | NA | NA | NA | NA | NA | UNTESTABLE | UNTESTABLE |
| yazar | double negative thymocyte | 4366 | 20000 | 0 | NA | NA | NA | NA | NA | NA | UNTESTABLE | UNTESTABLE |
| yazar | effector memory CD4-positive, alpha-beta T cell | 3023 | 20000 | 269 | -1.5e-05 | 0.9816 | 0.005707 | 0.0006012 | -0.0008002 | 0.8279 | DOES NOT PROVE | coordination falls |
| yazar | effector memory CD8-positive, alpha-beta T cell | 2848 | 20000 | 747 | -3.369e-05 | 0.8189 | 0.003766 | 0.0002185 | -0.0003593 | 0.8279 | DOES NOT PROVE | coordination falls |
| yazar | erythrocyte | 2214 | 20000 | 0 | NA | NA | NA | NA | NA | NA | UNTESTABLE | UNTESTABLE |
| yazar | gamma-delta T cell | 2897 | 20000 | 129 | 0.0005176 | 0.06698 | 0.004509 | 0.0005432 | 0.0008494 | 0.7434 | DOES NOT PROVE | coordination rises |
| yazar | hematopoietic precursor cell | 5543 | 20000 | 0 | NA | NA | NA | NA | NA | NA | UNTESTABLE | UNTESTABLE |
| yazar | innate lymphoid cell | 3308 | 20000 | 0 | NA | NA | NA | NA | NA | NA | UNTESTABLE | UNTESTABLE |
| yazar | memory B cell | 3076 | 20000 | 254 | 0.0003965 | 0.5992 | 0.0106 | 0.001235 | -8.108e-05 | 0.967 | INCONCLUSIVE | coordination rises |
| yazar | mucosal-associated invariant T cell | 3010 | 20000 | 46 | 0.000385 | 0.5992 | 0.003904 | 0.001461 | -4.467e-05 | 0.967 | INCONCLUSIVE | coordination rises |
| yazar | naive B cell | 2639 | 20000 | 322 | -0.0001231 | 0.5992 | 0.004443 | 0.0004398 | -0.000348 | 0.9476 | DOES NOT PROVE | coordination falls |
| yazar | naive thymus-derived CD4-positive, alpha-beta T cell | 2725 | 20000 | 821 | -1.66e-05 | 0.8386 | 0.002743 | 0.0001369 | -0.0001395 | 0.8964 | DOES NOT PROVE | coordination falls |
| yazar | naive thymus-derived CD8-positive, alpha-beta T cell | 2914 | 20000 | 393 | 8.632e-05 | 0.5992 | 0.004405 | 0.0003203 | 0.0002531 | 0.8964 | DOES NOT PROVE | coordination rises |
| yazar | natural killer cell | 3028 | 20000 | 816 | -0.0002595 | 0.03493 | 0.004101 | 0.0002367 | -0.001437 | 0.01399 | INCONCLUSIVE | coordination falls |
| yazar | plasmablast | 5839 | 20000 | 14 | NA | NA | 0.02536 | NA | NA | NA | UNTESTABLE | UNTESTABLE |
| yazar | plasmacytoid dendritic cell | 4230 | 20000 | 3 | NA | NA | NA | NA | NA | NA | UNTESTABLE | UNTESTABLE |
| yazar | platelet | 3076 | 20000 | 0 | NA | NA | NA | NA | NA | NA | UNTESTABLE | UNTESTABLE |
| yazar | regulatory T cell | 2849 | 20000 | 69 | 0.0006919 | 0.5992 | 0.006972 | 0.002361 | NA | NA | INCONCLUSIVE | coordination rises |
| yazar | transitional stage B cell | 3041 | 20000 | 213 | 6.053e-06 | 0.9816 | 0.006987 | 0.0007397 | 0.0001707 | 0.967 | DOES NOT PROVE | coordination rises |

### Sensitivity models

| Dataset | Cell type | Depth adjusted slope | Depth q | Unadjusted slope | Unadjusted q | 1500 rarefied slope | 1500 q |
| --- | --- | --- | --- | --- | --- | --- | --- |
| gao | Adipocyte | 0.0007914 | 0.6355 | -1.077e-05 | 0.9922 | 4.705e-06 | 0.9759 |
| gao | Cardiomyocyte | 0.0001809 | 0.8855 | -0.0005802 | 0.8823 | 3.554e-05 | 0.3482 |
| gao | Endocardial | -0.0005277 | 0.6355 | -0.001643 | 0.2385 | 1.133e-05 | 0.9759 |
| gao | Endothelial | 1.222e-05 | 0.9814 | 7.802e-05 | 0.9922 | 6.806e-05 | 0.5962 |
| gao | Epicardial | -0.001262 | 0.2068 | -0.001713 | 0.2385 | 0.000149 | 0.8307 |
| gao | Fibroblast | -0.0005152 | 0.6355 | -0.001168 | 0.2385 | 3.205e-06 | 0.9759 |
| gao | LEC | -0.0006545 | 0.6355 | -0.0009484 | 0.4454 | 0.0003535 | 0.3482 |
| gao | Lymphoid | 0.0001508 | 0.8741 | 0.0001827 | 0.8823 | 0.000185 | 0.3482 |
| gao | Mast | -0.0001192 | 0.8855 | 2.56e-05 | 0.9922 | -0.0002858 | 0.6972 |
| gao | Myeloid | -0.0001207 | 0.8741 | -0.0006931 | 0.2385 | 7.82e-05 | 0.3482 |
| gao | Neuronal | -0.001297 | 0.2068 | -0.001185 | 0.2385 | -4.205e-05 | 0.9759 |
| gao | Pericyte | -0.0002094 | 0.6355 | -0.0001864 | 0.8043 | 6.213e-05 | 0.5791 |
| gao | vSMC | -0.0009643 | 0.2068 | -0.001102 | 0.2385 | -1.965e-05 | 0.9759 |
| yazar | CD14-low, CD16-positive monocyte | 0.0002338 | 0.7204 | 0.0002978 | 0.4124 | 5.44e-05 | 0.8347 |
| yazar | CD14-positive monocyte | 0.0001675 | 0.7204 | 0.0002937 | 0.05776 | 6.801e-05 | 0.7118 |
| yazar | CD16-negative, CD56-bright natural killer cell, human | NA | NA | 0.001125 | 0.5854 | NA | NA |
| yazar | CD4-positive, alpha-beta T cell | NA | NA | NA | NA | NA | NA |
| yazar | CD4-positive, alpha-beta cytotoxic T cell | -0.0004114 | 0.8464 | 9.15e-05 | 0.989 | -4.841e-05 | 0.861 |
| yazar | CD8-positive, alpha-beta T cell | NA | NA | NA | NA | NA | NA |
| yazar | central memory CD4-positive, alpha-beta T cell | 2.189e-06 | 0.9719 | 6.841e-05 | 0.5854 | -2.819e-05 | 0.7118 |
| yazar | central memory CD8-positive, alpha-beta T cell | -0.0003824 | 0.7204 | 2.101e-05 | 0.989 | -0.0001666 | 0.7263 |
| yazar | conventional dendritic cell | NA | NA | NA | NA | NA | NA |
| yazar | dendritic cell | NA | NA | NA | NA | NA | NA |
| yazar | double negative thymocyte | NA | NA | NA | NA | NA | NA |
| yazar | effector memory CD4-positive, alpha-beta T cell | 1.588e-05 | 0.9719 | -0.0002627 | 0.4124 | 7.517e-05 | 0.7263 |
| yazar | effector memory CD8-positive, alpha-beta T cell | -1.428e-05 | 0.9719 | -2.221e-05 | 0.989 | -1.322e-05 | 0.861 |
| yazar | erythrocyte | NA | NA | NA | NA | NA | NA |
| yazar | gamma-delta T cell | 0.0004708 | 0.1518 | 0.0003628 | 0.05776 | 0.000116 | 0.7263 |
| yazar | hematopoietic precursor cell | NA | NA | NA | NA | NA | NA |
| yazar | innate lymphoid cell | NA | NA | NA | NA | NA | NA |
| yazar | memory B cell | -4.875e-05 | 0.9719 | -4.817e-06 | 0.989 | -2.464e-05 | 0.861 |
| yazar | mucosal-associated invariant T cell | 0.0003066 | 0.8464 | 8.665e-05 | 0.989 | 0.0003014 | 0.7263 |
| yazar | naive B cell | -0.0001307 | 0.8464 | -0.0001089 | 0.8093 | 2.396e-05 | 0.861 |
| yazar | naive thymus-derived CD4-positive, alpha-beta T cell | 1.629e-06 | 0.9719 | -7.581e-06 | 0.989 | -3.918e-05 | 0.7118 |
| yazar | naive thymus-derived CD8-positive, alpha-beta T cell | -3.665e-05 | 0.9719 | 6.197e-05 | 0.8404 | 4.738e-05 | 0.7263 |
| yazar | natural killer cell | -0.0001791 | 0.1518 | -0.0002497 | 0.0159 | 8.078e-06 | 0.861 |
| yazar | plasmablast | NA | NA | -0.001103 | 0.989 | NA | NA |
| yazar | plasmacytoid dendritic cell | NA | NA | NA | NA | NA | NA |
| yazar | platelet | NA | NA | NA | NA | NA | NA |
| yazar | regulatory T cell | 0.0006122 | 0.8464 | 0.0005376 | 0.6104 | 0.000262 | 0.8007 |
| yazar | transitional stage B cell | -2.579e-05 | 0.9719 | 5.401e-06 | 0.989 | -1.423e-05 | 0.861 |

### Verdicts

- Gao: INCONCLUSIVE; no proven direction. PROVES 0, DOES NOT PROVE 9, INCONCLUSIVE 4, UNTESTABLE 0 cell types.
  Inconclusive units: Adipocyte, Cardiomyocyte, Epicardial, Mast.
- Yazar: INCONCLUSIVE; no proven direction. PROVES 0, DOES NOT PROVE 11, INCONCLUSIVE 5, UNTESTABLE 12 cell types.
  Inconclusive units: CD4-positive, alpha-beta cytotoxic T cell, memory B cell, mucosal-associated invariant T cell, natural killer cell, regulatory T cell.
  Untestable units: CD16-negative, CD56-bright natural killer cell, human, CD4-positive, alpha-beta T cell, CD8-positive, alpha-beta T cell, conventional dendritic cell, dendritic cell, double negative thymocyte, erythrocyte, hematopoietic precursor cell, innate lymphoid cell, plasmablast, plasmacytoid dendritic cell, platelet.
- Commit-label summary: INCONCLUSIVE. This does not pool donor data or replace the separate Gao and Yazar verdicts.
- A rise would differ from the Levy 2020 prediction of falling coordination. No significant rise is described as confirmation of that prediction.

### Reproducibility

Pair base seed 20261117, permutation base seed 20261118, rarefaction base seed 20261119, and split-half seeds 20261116, 20261117, 20261118, 20261119, 20261120. Exact derived per-cell-type pair and permutation seeds are in the result CSV and analysis checkpoints.
Primary and group comparison family sizes: gao primary 13, group 10, yazar primary 16, group 14.
No primary cells were downsampled. The only downsampling is the explicitly labeled 1500 rarefaction sensitivity.
