# Block 3 (mature-to-old): Full profile-cloud donor dispersion in Gao heart and Yazar blood

Status: complete. This specification was written before running the analysis. Mature-to-old rerun of reports/followup-11pt-block3-profile-cloud.md: donors aged under 30 are excluded, and the primary baseline YOUNG (20 to 39) is replaced by MATURE (30 to 39). Every other frozen rule is unchanged.

## Exact analysis specification

- Per dataset, cell type, and donor, using the shared eligibility floors: donor centroid = mean within-cell proportions across that donor's cells at full depth, over the technical-masked gene set.
- Primary outcome, per dataset, cell type, and age bin, with LT20 and ages 20 to 29 excluded: mean pairwise Hellinger distance among donor centroids WITHIN the bin. Secondary outcome: mean distance of each centroid to its own bin centroid.
- Primary comparison: OLD vs MATURE mean pairwise distance, with a donor-permutation null that permutes donor age-bin labels within dataset and cell type, uses 2000 permutations and a fixed seed, and recomputes the within-bin means each time. Report observed difference, permutation p-value, N donors per bin, and MDE80.
- Batch handling has two arms: (a) raw centroids and (b) batch-residualized centroids, residualizing each gene's centroid proportions on batch across donors, using Gao study indicators and Yazar pool means, then re-normalizing to proportions before computing distances. The batch-residualized arm is primary and the raw arm is the sensitivity.
- MID comparisons are secondary: OLD vs MID and MID vs MATURE use the same permutation procedure and the same arms. OLD vs MID is retained as a distinct existing mature contrast and is not overwritten by the new primary.
- Split-half technical null: split each donor's cells into two halves with 5 fixed seeds, build half-centroids, compute mean pairwise distances among half-centroids within bins, and report the noise floor. A claimed age difference must exceed it.
- Family, stated before correction: BH q within each (dataset, arm) slice across cell types for the OLD vs MATURE comparison. MID comparisons are a separate stated family.

## Fixed implementation details

- Inputs are the two Block 1 h5ad files and the donor tables `results/followup-between-subject-gao-donors.csv` and `results/followup-between-subject-yazar-donors.csv`. Yazar excludes rows labeled peripheral blood mononuclear cell. No new downloads are used.
- The human technical gene mask excludes uppercase MT-, RPL, RPS, MRPL, MRPS, XIST, TSIX, DDX3Y, EIF2S3Y, UTY, and KDM5D before depth calculation and normalization. ERCC controls are also excluded, consistent with the existing donor tables. All remaining genes are retained without a detection filter.
- Full-depth raw counts are converted to within-cell proportions. Cells with at least 1500 post-mask counts are eligible, with no downsampling. A donor-cell-type unit requires at least 20 eligible cells. Each comparison requires at least 6 eligible donors in each bin. Units below these floors are UNTESTABLE, never nulls.
- Age bins: MATURE 30 to 39, MID 40 to 59, OLD 60 and older. Donors under 30 (LT20 and 20 to 29) are excluded from all Block 3 calculations. The mature bin relabeling is applied to the reconciled donor metadata before any donor selection, so batch residualization, permutations, and split-half floors all use mature donors only.
- The original Block 2 count cache is reused for the projected per-cell counts: its values are donor-local (per-cell projection and UMI filtering only), and every downstream step in this script (donor selection, centroid computation, batch residualization, permutations, split-half) uses mature donors only. Checkpoints are kept in a separate mature cache directory and never overwrite the original checkpoints.
- Each cell contributes equally to its donor centroid. This is the arithmetic mean of within-cell proportion vectors, not a sum-of-counts pseudobulk. Donor is the unit for every comparison, permutation, and noise-floor summary.
- Hellinger distance uses the formula in `hellinger_distance` from `scripts/run_h3_expression_vector.py`, `sqrt(1 - sum_g sqrt(p_g q_g))`, on proportion vectors. The implementation uses an equivalent vectorized square-root affinity matrix for all pairs and validates it against the imported function. Distances from donors to their own bin centroid call the imported function directly. The bin centroid is the arithmetic mean of donor centroids in that bin.
- Batch residualization is fitted once per cell type and arm input across all eligible MATURE, MID, and OLD donors. The design is intercept plus study indicators for Gao and intercept plus pool indicators for Yazar, equivalent to removing batch means and restoring the global gene mean. Residualized values are clipped at zero and each donor profile is divided by its sum to obtain nonnegative proportions. The same procedure is applied independently to each split half.
- Comparisons use the first bin minus the second bin: OLD minus MATURE, OLD minus MID, and MID minus MATURE. For each comparison, only donors in those two bins enter the permutation. Each of 2000 permutations holds the two bin sizes fixed and recomputes both within-bin mean pairwise distances from the fixed donor distance matrix. P-values are two-sided with the plus-one correction. Permutation base seed is 20261118, with deterministic dataset, cell-type, arm, and comparison offsets recorded in the CSV.
- The secondary own-bin-centroid outcome is reported as observed bin means and observed comparison differences. The permutation p-value and BH q-value apply to the primary mean-pairwise outcome only.
- Split-half seeds are 20261116, 20261117, 20261118, 20261119, and 20261120. For each donor and seed, cells are randomly permuted and divided into halves of floor(N/2) and ceil(N/2). For each arm and comparison, the noise floor is the mean over 5 seeds of the absolute half-to-half change in that comparison's mean-pairwise difference. Bin-specific absolute half-to-half changes and the corresponding secondary-outcome noise floor are also recorded. This directly estimates technical instability of the age contrast.
- MDE80 is `(z_0.975 + z_0.80) * SD(permutation null)`, equal to 2.801585 times the permutation-null standard deviation. This is the smallest absolute difference with approximately 80 percent power for a two-sided 0.05 test under the donor-permutation scale. For a null or negative batch-residualized result, MDE80 at or below its split-half noise floor is pre-stated as small enough to exclude an increase detectable above technical noise; a larger MDE80 is underpowered.
- The OLD vs MATURE BH families contain evaluated cell types within each dataset and arm. The separate MID family contains both OLD vs MID and MID vs MATURE evaluated rows within each dataset and arm. UNTESTABLE rows are excluded and have missing q-values. Family sizes are reported after evaluation.
- No 1500 rarefaction is performed in Block 3 because it is not part of this block's specification. The 1500 threshold is only a cell eligibility filter.
- This block does not score P1, so no P1/P1m change applies here.

## Pre-stated verdict rules

For each cell type, PROVES that the donor cloud spreads with age requires a positive batch-residualized OLD vs MATURE difference, BH q < 0.05, and a difference greater than the split-half noise floor. DOES NOT PROVE applies when the batch-residualized result is null or negative and MDE80 is at or below the noise floor. INCONCLUSIVE applies to a null with MDE80 above the noise floor, a positive significant result that does not exceed the noise floor, or opposite signs in the raw and batch-residualized arms. UNTESTABLE units remain UNTESTABLE.

For each dataset summary, PROVES requires at least one PROVES cell type. DOES NOT PROVE requires at least one evaluated cell type and every evaluated cell type to meet DOES NOT PROVE. All other dataset summaries are INCONCLUSIVE. Gao and Yazar are reported separately and are never pooled.

## Comparison against the original Block 3

Original verdicts (full range, OLD vs YOUNG): Gao INCONCLUSIVE (0 PROVES, 0 DOES NOT PROVE, 10 INCONCLUSIVE, 3 UNTESTABLE); Yazar PROVES with 2 proven units (memory B cell difference 0.0089008, q 0.015548; naive thymus-derived CD8-positive, alpha-beta T cell difference 0.021979, q 0.0009995), 0 DOES NOT PROVE, 12 INCONCLUSIVE, 14 UNTESTABLE. The mature comparison table in the Results section records, per cell type, whether each PROVES verdict survives or dies and whether the contraction pattern (6 other Yazar types with significant negative differences) survives.

## Execution note

Delegated-model note: Kimi hit its 5-hour usage quota at about 13:13 PST Nov 21 and Antigravity Opus 5.5 High hit its individual quota (reset in 139h59m), then Antigravity Gemini 4 Neon XHigh hit its quota (reset in about 4h31m from 21:41 UTC) after spending its run on a multi-round plan review with no implementation. No delegated model produced any Block 1 to 3 output. This block was therefore executed directly by the runner from the frozen original script by copy-and-remap (mature bin relabeling, mature output paths, separate mature checkpoint cache), with the coordinator to verify before reporting. Model for this block: direct runner execution, no delegated model.
## Results

All pairwise effects are first-bin minus second-bin Hellinger distance differences. The batch-residualized arm is primary; raw is the sensitivity. NA denotes an UNTESTABLE or not applicable field.

### OLD vs MATURE primary results and raw sensitivity

| Dataset | Cell type | N MATURE | N OLD | Residualized MATURE | Residualized OLD | Residualized difference | p | BH q | Noise floor | MDE80 | Raw difference | Raw q | Raw noise | Raw MDE80 | Own-bin secondary difference | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gao | Adipocyte | 3 | 12 | 0.23591 | 0.2586 | 0.022685 | NA | NA | 0.005401 | NA | -0.049149 | NA | 0.0028521 | NA | 0.037884 | UNTESTABLE |
| gao | Cardiomyocyte | 11 | 39 | 0.17025 | 0.15216 | -0.018097 | 0.048476 | 0.21814 | 0.00058976 | 0.027057 | -0.0084993 | 0.58696 | 0.00043808 | 0.026511 | -0.009423 | INCONCLUSIVE |
| gao | Endocardial | 8 | 29 | 0.23539 | 0.2257 | -0.0096921 | 0.52974 | 0.98601 | 0.005604 | 0.043153 | 0.024366 | 0.31394 | 0.0050566 | 0.04979 | -0.00056036 | INCONCLUSIVE |
| gao | Endothelial | 10 | 30 | 0.20669 | 0.20585 | -0.00083192 | 0.93253 | 0.98601 | 0.0045133 | 0.031567 | 0.0034819 | 0.87356 | 0.0049522 | 0.034866 | 0.0046367 | INCONCLUSIVE |
| gao | Epicardial | 2 | 8 | 0.30604 | 0.24309 | -0.062949 | NA | NA | 0.01116 | NA | 0.043576 | NA | 0.011692 | NA | 0.00091216 | UNTESTABLE |
| gao | Fibroblast | 12 | 38 | 0.20311 | 0.20001 | -0.0030994 | 0.7916 | 0.98601 | 0.0015195 | 0.032033 | -0.0037783 | 0.87356 | 0.00096468 | 0.036679 | 0.0031159 | INCONCLUSIVE |
| gao | LEC | 4 | 13 | 0.28607 | 0.25481 | -0.031263 | NA | NA | 0.0053577 | NA | -0.010069 | NA | 0.0096037 | NA | -0.0032903 | UNTESTABLE |
| gao | Lymphoid | 9 | 22 | 0.31192 | 0.26887 | -0.043048 | 0.0074963 | 0.067466 | 0.0068243 | 0.052629 | -0.030798 | 0.16492 | 0.0052443 | 0.045779 | -0.023357 | INCONCLUSIVE |
| gao | Mast | 1 | 8 | NA | 0.26138 | NA | NA | NA | NA | NA | NA | NA | NA | NA | 0.17685 | UNTESTABLE |
| gao | Myeloid | 11 | 39 | 0.22286 | 0.22263 | -0.00022955 | 0.98601 | 0.98601 | 0.0029594 | 0.041033 | 0.041341 | 0.033733 | 0.0011468 | 0.046446 | 0.0022727 | INCONCLUSIVE |
| gao | Neuronal | 9 | 19 | 0.29677 | 0.28258 | -0.014183 | 0.33083 | 0.98601 | 0.0064753 | 0.039655 | 0.016111 | 0.31394 | 0.0046732 | 0.031963 | -0.0037966 | INCONCLUSIVE |
| gao | Pericyte | 11 | 33 | 0.2337 | 0.23616 | 0.0024578 | 0.88006 | 0.98601 | 0.0029359 | 0.041586 | 0.035066 | 0.033733 | 0.001932 | 0.036949 | 0.0057198 | INCONCLUSIVE |
| gao | vSMC | 9 | 30 | 0.26305 | 0.2673 | 0.0042479 | 0.76362 | 0.98601 | 0.0069279 | 0.038593 | 0.0027579 | 0.87356 | 0.0070198 | 0.044194 | 0.010134 | INCONCLUSIVE |
| yazar | CD14-low, CD16-positive monocyte | 9 | 204 | 0.16381 | 0.16449 | 0.00068097 | 0.89105 | 0.91204 | 0.0026158 | 0.014605 | 0.0041295 | 0.50555 | 0.0042485 | 0.013606 | 0.0064926 | INCONCLUSIVE |
| yazar | CD14-positive monocyte | 21 | 299 | 0.17758 | 0.17078 | -0.0067936 | 0.057971 | 0.099379 | 0.0025254 | 0.01002 | -0.015623 | 0.001999 | 0.0036283 | 0.012334 | -0.0014651 | INCONCLUSIVE |
| yazar | CD16-negative, CD56-bright natural killer cell, human | 0 | 4 | NA | 5.2684e-09 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE |
| yazar | CD4-positive, alpha-beta T cell | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE |
| yazar | CD4-positive, alpha-beta cytotoxic T cell | 1 | 78 | NA | 0.16258 | NA | NA | NA | NA | NA | NA | NA | NA | NA | 0.10943 | UNTESTABLE |
| yazar | CD8-positive, alpha-beta T cell | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE |
| yazar | central memory CD4-positive, alpha-beta T cell | 58 | 660 | 0.16398 | 0.1381 | -0.025878 | 0.00049975 | 0.0011994 | 0.00067815 | 0.0085716 | -0.026352 | 0.0011994 | 0.00097377 | 0.0094141 | -0.01877 | INCONCLUSIVE |
| yazar | central memory CD8-positive, alpha-beta T cell | 4 | 80 | 0.1174 | 0.17064 | 0.053235 | NA | NA | 0.0028044 | NA | -0.021743 | NA | 0.0071379 | NA | 0.04113 | UNTESTABLE |
| yazar | conventional dendritic cell | 0 | 5 | NA | 0.079515 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE |
| yazar | dendritic cell | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE |
| yazar | double negative thymocyte | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE |
| yazar | effector memory CD4-positive, alpha-beta T cell | 7 | 201 | 0.19382 | 0.20263 | 0.0088056 | 0.13443 | 0.20165 | 0.0035841 | 0.01735 | -0.016633 | 0.020561 | 0.009572 | 0.017828 | 0.01692 | INCONCLUSIVE |
| yazar | effector memory CD8-positive, alpha-beta T cell | 34 | 552 | 0.18741 | 0.171 | -0.016417 | 0.00049975 | 0.0011994 | 0.00078957 | 0.011395 | -0.025337 | 0.0011994 | 0.0014543 | 0.01436 | -0.010055 | INCONCLUSIVE |
| yazar | erythrocyte | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE |
| yazar | gamma-delta T cell | 7 | 88 | 0.18817 | 0.18064 | -0.0075328 | 0.31084 | 0.41446 | 0.0041533 | 0.022028 | -0.012044 | 0.19865 | 0.005295 | 0.022059 | -0.00033746 | INCONCLUSIVE |
| yazar | hematopoietic precursor cell | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE |
| yazar | innate lymphoid cell | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE |
| yazar | memory B cell | 10 | 176 | 0.18453 | 0.19881 | 0.014282 | 0.012494 | 0.024988 | 0.0025937 | 0.0152 | 0.0029383 | 0.67657 | 0.0034412 | 0.016424 | 0.017006 | PROVES |
| yazar | mucosal-associated invariant T cell | 3 | 23 | 0.17558 | 0.14197 | -0.033612 | NA | NA | 0.0033384 | NA | -0.013753 | NA | 0.005545 | NA | -0.0098918 | UNTESTABLE |
| yazar | naive B cell | 11 | 251 | 0.19975 | 0.19906 | -0.00068529 | 0.91204 | 0.91204 | 0.0014664 | 0.016023 | 0.0020351 | 0.7941 | 0.0020282 | 0.020959 | 0.0039198 | INCONCLUSIVE |
| yazar | naive thymus-derived CD4-positive, alpha-beta T cell | 51 | 564 | 0.17284 | 0.15443 | -0.018411 | 0.00049975 | 0.0011994 | 0.00063132 | 0.010692 | -0.023741 | 0.0011994 | 0.00087049 | 0.013087 | -0.012346 | INCONCLUSIVE |
| yazar | naive thymus-derived CD8-positive, alpha-beta T cell | 40 | 210 | 0.17539 | 0.18787 | 0.012476 | 0.00049975 | 0.0011994 | 0.00092446 | 0.0087659 | 0.014582 | 0.0011994 | 0.0016239 | 0.010685 | 0.010717 | PROVES |
| yazar | natural killer cell | 38 | 598 | 0.19671 | 0.16495 | -0.03176 | 0.00049975 | 0.0011994 | 0.0023343 | 0.010824 | -0.036202 | 0.0011994 | 0.0025394 | 0.013035 | -0.022651 | INCONCLUSIVE |
| yazar | plasmablast | 1 | 8 | NA | 0.10698 | NA | NA | NA | NA | NA | NA | NA | NA | NA | 0.063994 | UNTESTABLE |
| yazar | plasmacytoid dendritic cell | 0 | 2 | NA | 1.0537e-08 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE |
| yazar | platelet | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE |
| yazar | regulatory T cell | 0 | 59 | NA | 0.16894 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE |
| yazar | transitional stage B cell | 9 | 147 | 0.19696 | 0.19503 | -0.0019303 | 0.74913 | 0.89895 | 0.0015267 | 0.017441 | -0.005546 | 0.50555 | 0.0027846 | 0.01909 | 0.0036034 | INCONCLUSIVE |

### MID secondary comparisons

The MID BH family combines OLD vs MID and MID vs MATURE within each dataset and arm. Permutation p-values and q-values apply to mean pairwise distance, not the descriptive own-bin-centroid outcome.

| Dataset | Cell type | Arm | Comparison | N bin 1 | N bin 2 | Bin 1 pairwise | Bin 2 pairwise | Difference | p | BH q | Noise floor | MDE80 | Own-bin secondary difference | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gao | Adipocyte | batch_residualized | MID_vs_MATURE | 20 | 3 | 0.24833 | 0.23591 | 0.012419 | NA | NA | 0.0041934 | NA | 0.034267 | UNTESTABLE_BIN_FLOOR |
| gao | Adipocyte | batch_residualized | OLD_vs_MID | 12 | 20 | 0.2586 | 0.24833 | 0.010266 | 0.35932 | 0.55247 | 0.0038682 | 0.031246 | 0.003617 | EVALUATED |
| gao | Adipocyte | raw | MID_vs_MATURE | 20 | 3 | 0.29818 | 0.38912 | -0.090945 | NA | NA | 0.0059677 | NA | -0.019119 | UNTESTABLE_BIN_FLOOR |
| gao | Adipocyte | raw | OLD_vs_MID | 12 | 20 | 0.33998 | 0.29818 | 0.041796 | 0.0024988 | 0.032984 | 0.0056442 | 0.041425 | 0.024567 | EVALUATED |
| gao | Cardiomyocyte | batch_residualized | MID_vs_MATURE | 52 | 11 | 0.1545 | 0.17025 | -0.015756 | 0.063968 | 0.21361 | 0.00041883 | 0.023902 | -0.0074545 | EVALUATED |
| gao | Cardiomyocyte | batch_residualized | OLD_vs_MID | 39 | 52 | 0.15216 | 0.1545 | -0.0023409 | 0.65517 | 0.796 | 0.00043224 | 0.014293 | -0.0019685 | EVALUATED |
| gao | Cardiomyocyte | raw | MID_vs_MATURE | 52 | 11 | 0.19619 | 0.19734 | -0.0011553 | 0.90955 | 0.97251 | 0.00067236 | 0.032917 | 0.0040222 | EVALUATED |
| gao | Cardiomyocyte | raw | OLD_vs_MID | 39 | 52 | 0.18884 | 0.19619 | -0.007344 | 0.31634 | 0.54465 | 0.00024055 | 0.019296 | -0.0050791 | EVALUATED |
| gao | Endocardial | batch_residualized | MID_vs_MATURE | 29 | 8 | 0.25907 | 0.23539 | 0.02368 | 0.10945 | 0.28464 | 0.0076402 | 0.040735 | 0.023829 | EVALUATED |
| gao | Endocardial | batch_residualized | OLD_vs_MID | 29 | 29 | 0.2257 | 0.25907 | -0.033373 | 0.00049975 | 0.010995 | 0.0045295 | 0.027842 | -0.024389 | EVALUATED |
| gao | Endocardial | raw | MID_vs_MATURE | 29 | 8 | 0.32559 | 0.27653 | 0.049055 | 0.0029985 | 0.032984 | 0.0065125 | 0.047944 | 0.046086 | EVALUATED |
| gao | Endocardial | raw | OLD_vs_MID | 29 | 29 | 0.3009 | 0.32559 | -0.024689 | 0.04048 | 0.14843 | 0.0034826 | 0.03381 | -0.014673 | EVALUATED |
| gao | Endothelial | batch_residualized | MID_vs_MATURE | 40 | 10 | 0.21042 | 0.20669 | 0.0037353 | 0.71264 | 0.796 | 0.0024998 | 0.029476 | 0.0074239 | EVALUATED |
| gao | Endothelial | batch_residualized | OLD_vs_MID | 30 | 40 | 0.20585 | 0.21042 | -0.0045672 | 0.53573 | 0.6933 | 0.0035224 | 0.020314 | -0.0027872 | EVALUATED |
| gao | Endothelial | raw | MID_vs_MATURE | 40 | 10 | 0.27041 | 0.25252 | 0.017887 | 0.096452 | 0.26524 | 0.0027241 | 0.030954 | 0.020871 | EVALUATED |
| gao | Endothelial | raw | OLD_vs_MID | 30 | 40 | 0.256 | 0.27041 | -0.014405 | 0.066967 | 0.21047 | 0.0043212 | 0.022238 | -0.011365 | EVALUATED |
| gao | Epicardial | batch_residualized | MID_vs_MATURE | 8 | 2 | 0.29707 | 0.30604 | -0.0089746 | NA | NA | 0.013346 | NA | 0.041582 | UNTESTABLE_BIN_FLOOR |
| gao | Epicardial | batch_residualized | OLD_vs_MID | 8 | 8 | 0.24309 | 0.29707 | -0.053974 | 0.009995 | 0.073297 | 0.0089056 | 0.06043 | -0.04067 | EVALUATED |
| gao | Epicardial | raw | MID_vs_MATURE | 8 | 2 | 0.3316 | 0.31199 | 0.019608 | NA | NA | 0.013389 | NA | 0.066942 | UNTESTABLE_BIN_FLOOR |
| gao | Epicardial | raw | OLD_vs_MID | 8 | 8 | 0.35557 | 0.3316 | 0.023968 | 0.4088 | 0.59957 | 0.0097653 | 0.072717 | 0.016286 | EVALUATED |
| gao | Fibroblast | batch_residualized | MID_vs_MATURE | 49 | 12 | 0.21137 | 0.20311 | 0.0082611 | 0.35782 | 0.55247 | 0.0011436 | 0.024786 | 0.012142 | EVALUATED |
| gao | Fibroblast | batch_residualized | OLD_vs_MID | 38 | 49 | 0.20001 | 0.21137 | -0.01136 | 0.067966 | 0.21361 | 0.0011223 | 0.017339 | -0.0090258 | EVALUATED |
| gao | Fibroblast | raw | MID_vs_MATURE | 49 | 12 | 0.27645 | 0.27701 | -0.00056206 | 0.96402 | 0.97251 | 0.0016824 | 0.031473 | 0.0092684 | EVALUATED |
| gao | Fibroblast | raw | OLD_vs_MID | 38 | 49 | 0.27323 | 0.27645 | -0.0032162 | 0.65767 | 0.80382 | 0.0009638 | 0.019462 | -0.00095857 | EVALUATED |
| gao | LEC | batch_residualized | MID_vs_MATURE | 17 | 4 | 0.26653 | 0.28607 | -0.019537 | NA | NA | 0.0079715 | NA | 0.0061988 | UNTESTABLE_BIN_FLOOR |
| gao | LEC | batch_residualized | OLD_vs_MID | 13 | 17 | 0.25481 | 0.26653 | -0.011726 | 0.11644 | 0.28464 | 0.0038446 | 0.020703 | -0.0094891 | EVALUATED |
| gao | LEC | raw | MID_vs_MATURE | 17 | 4 | 0.3297 | 0.33432 | -0.0046194 | NA | NA | 0.0075551 | NA | 0.026637 | UNTESTABLE_BIN_FLOOR |
| gao | LEC | raw | OLD_vs_MID | 13 | 17 | 0.32425 | 0.3297 | -0.0054494 | 0.6082 | 0.78708 | 0.005093 | 0.029989 | -0.0059171 | EVALUATED |
| gao | Lymphoid | batch_residualized | MID_vs_MATURE | 25 | 9 | 0.28298 | 0.31192 | -0.028937 | 0.063968 | 0.21361 | 0.0077297 | 0.043358 | -0.012058 | EVALUATED |
| gao | Lymphoid | batch_residualized | OLD_vs_MID | 22 | 25 | 0.26887 | 0.28298 | -0.014111 | 0.18641 | 0.4068 | 0.0026491 | 0.029841 | -0.011299 | EVALUATED |
| gao | Lymphoid | raw | MID_vs_MATURE | 25 | 9 | 0.38595 | 0.38646 | -0.00051861 | 0.97251 | 0.97251 | 0.0046432 | 0.043451 | 0.010847 | EVALUATED |
| gao | Lymphoid | raw | OLD_vs_MID | 22 | 25 | 0.35567 | 0.38595 | -0.030279 | 0.007996 | 0.058637 | 0.003655 | 0.03313 | -0.021287 | EVALUATED |
| gao | Mast | batch_residualized | MID_vs_MATURE | 8 | 1 | 0.27883 | NA | NA | NA | NA | NA | NA | 0.18298 | UNTESTABLE_BIN_FLOOR |
| gao | Mast | batch_residualized | OLD_vs_MID | 8 | 8 | 0.26138 | 0.27883 | -0.017448 | 0.38631 | 0.55247 | 0.0024572 | 0.052868 | -0.006121 | EVALUATED |
| gao | Mast | raw | MID_vs_MATURE | 8 | 1 | 0.43906 | NA | NA | NA | NA | NA | NA | 0.31143 | UNTESTABLE_BIN_FLOOR |
| gao | Mast | raw | OLD_vs_MID | 8 | 8 | 0.40665 | 0.43906 | -0.032409 | 0.32184 | 0.54465 | 0.0074658 | 0.090054 | -0.013794 | EVALUATED |
| gao | Myeloid | batch_residualized | MID_vs_MATURE | 49 | 11 | 0.23148 | 0.22286 | 0.0086194 | 0.4018 | 0.55247 | 0.0027789 | 0.028969 | 0.01078 | EVALUATED |
| gao | Myeloid | batch_residualized | OLD_vs_MID | 39 | 49 | 0.22263 | 0.23148 | -0.008849 | 0.27086 | 0.49659 | 0.0031897 | 0.022541 | -0.0085077 | EVALUATED |
| gao | Myeloid | raw | MID_vs_MATURE | 49 | 11 | 0.29647 | 0.2659 | 0.030573 | 0.022489 | 0.12369 | 0.0032308 | 0.037559 | 0.032362 | EVALUATED |
| gao | Myeloid | raw | OLD_vs_MID | 39 | 49 | 0.30724 | 0.29647 | 0.010768 | 0.27086 | 0.54465 | 0.0030659 | 0.026876 | 0.01266 | EVALUATED |
| gao | Neuronal | batch_residualized | MID_vs_MATURE | 28 | 9 | 0.29372 | 0.29677 | -0.003042 | 0.8011 | 0.8066 | 0.0062825 | 0.033362 | 0.0066341 | EVALUATED |
| gao | Neuronal | batch_residualized | OLD_vs_MID | 19 | 28 | 0.28258 | 0.29372 | -0.01114 | 0.2034 | 0.4068 | 0.0061547 | 0.023936 | -0.010431 | EVALUATED |
| gao | Neuronal | raw | MID_vs_MATURE | 28 | 9 | 0.36224 | 0.34316 | 0.019079 | 0.12144 | 0.29685 | 0.0055261 | 0.034686 | 0.025245 | EVALUATED |
| gao | Neuronal | raw | OLD_vs_MID | 19 | 28 | 0.35927 | 0.36224 | -0.0029686 | 0.74113 | 0.85815 | 0.0051712 | 0.024438 | -0.0021333 | EVALUATED |
| gao | Pericyte | batch_residualized | MID_vs_MATURE | 44 | 11 | 0.23819 | 0.2337 | 0.0044846 | 0.72364 | 0.796 | 0.0023459 | 0.033245 | 0.008875 | EVALUATED |
| gao | Pericyte | batch_residualized | OLD_vs_MID | 33 | 44 | 0.23616 | 0.23819 | -0.0020268 | 0.8066 | 0.8066 | 0.0029124 | 0.023425 | -0.0031552 | EVALUATED |
| gao | Pericyte | raw | MID_vs_MATURE | 44 | 11 | 0.31105 | 0.2809 | 0.030146 | 0.028986 | 0.12754 | 0.0019595 | 0.039923 | 0.029253 | EVALUATED |
| gao | Pericyte | raw | OLD_vs_MID | 33 | 44 | 0.31597 | 0.31105 | 0.00492 | 0.5907 | 0.78708 | 0.0028762 | 0.024814 | 0.0085876 | EVALUATED |
| gao | vSMC | batch_residualized | MID_vs_MATURE | 40 | 9 | 0.2863 | 0.26305 | 0.023252 | 0.062969 | 0.21361 | 0.0069439 | 0.034475 | 0.025488 | EVALUATED |
| gao | vSMC | batch_residualized | OLD_vs_MID | 30 | 40 | 0.2673 | 0.2863 | -0.019004 | 0.009995 | 0.073297 | 0.0026261 | 0.021566 | -0.015354 | EVALUATED |
| gao | vSMC | raw | MID_vs_MATURE | 40 | 9 | 0.34087 | 0.32863 | 0.012235 | 0.35682 | 0.56072 | 0.0082883 | 0.036467 | 0.019129 | EVALUATED |
| gao | vSMC | raw | OLD_vs_MID | 30 | 40 | 0.33139 | 0.34087 | -0.009477 | 0.29535 | 0.54465 | 0.002486 | 0.025557 | -0.0056309 | EVALUATED |
| yazar | CD14-low, CD16-positive monocyte | batch_residualized | MID_vs_MATURE | 44 | 9 | 0.16831 | 0.16381 | 0.0045021 | 0.52824 | 0.64307 | 0.0023419 | 0.019624 | 0.0079201 | EVALUATED |
| yazar | CD14-low, CD16-positive monocyte | batch_residualized | OLD_vs_MID | 204 | 44 | 0.16449 | 0.16831 | -0.0038211 | 0.15092 | 0.26412 | 0.0024635 | 0.0074008 | -0.0014275 | EVALUATED |
| yazar | CD14-low, CD16-positive monocyte | raw | MID_vs_MATURE | 44 | 9 | 0.21571 | 0.20707 | 0.0086386 | 0.09995 | 0.1554 | 0.0053353 | 0.014752 | 0.014151 | EVALUATED |
| yazar | CD14-low, CD16-positive monocyte | raw | OLD_vs_MID | 204 | 44 | 0.2112 | 0.21571 | -0.0045092 | 0.074463 | 0.139 | 0.0031936 | 0.0070881 | -0.0022417 | EVALUATED |
| yazar | CD14-positive monocyte | batch_residualized | MID_vs_MATURE | 88 | 21 | 0.17672 | 0.17758 | -0.0008576 | 0.84808 | 0.87949 | 0.0026039 | 0.011758 | 0.0024832 | EVALUATED |
| yazar | CD14-positive monocyte | batch_residualized | OLD_vs_MID | 299 | 88 | 0.17078 | 0.17672 | -0.005936 | 0.001999 | 0.007996 | 0.00079426 | 0.0052097 | -0.0039484 | EVALUATED |
| yazar | CD14-positive monocyte | raw | MID_vs_MATURE | 88 | 21 | 0.21401 | 0.22561 | -0.011602 | 0.021489 | 0.046285 | 0.0029113 | 0.014065 | -0.0051174 | EVALUATED |
| yazar | CD14-positive monocyte | raw | OLD_vs_MID | 299 | 88 | 0.20999 | 0.21401 | -0.004021 | 0.087456 | 0.15305 | 0.00090239 | 0.0066012 | -0.0021587 | EVALUATED |
| yazar | CD16-negative, CD56-bright natural killer cell, human | batch_residualized | MID_vs_MATURE | 1 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | CD16-negative, CD56-bright natural killer cell, human | batch_residualized | OLD_vs_MID | 4 | 1 | 5.2684e-09 | NA | NA | NA | NA | NA | NA | 4.5195e-10 | UNTESTABLE_BIN_FLOOR |
| yazar | CD16-negative, CD56-bright natural killer cell, human | raw | MID_vs_MATURE | 1 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | CD16-negative, CD56-bright natural killer cell, human | raw | OLD_vs_MID | 4 | 1 | 0.27903 | NA | NA | NA | NA | NA | NA | 0.18207 | UNTESTABLE_BIN_FLOOR |
| yazar | CD4-positive, alpha-beta T cell | batch_residualized | MID_vs_MATURE | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | CD4-positive, alpha-beta T cell | batch_residualized | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | CD4-positive, alpha-beta T cell | raw | MID_vs_MATURE | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | CD4-positive, alpha-beta T cell | raw | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | CD4-positive, alpha-beta cytotoxic T cell | batch_residualized | MID_vs_MATURE | 13 | 1 | 0.14605 | NA | NA | NA | NA | NA | NA | 0.097061 | UNTESTABLE_BIN_FLOOR |
| yazar | CD4-positive, alpha-beta cytotoxic T cell | batch_residualized | OLD_vs_MID | 78 | 13 | 0.16258 | 0.14605 | 0.016527 | 0.061469 | 0.15647 | 0.00075149 | 0.025141 | 0.012364 | EVALUATED |
| yazar | CD4-positive, alpha-beta cytotoxic T cell | raw | MID_vs_MATURE | 13 | 1 | 0.23504 | NA | NA | NA | NA | NA | NA | 0.1678 | UNTESTABLE_BIN_FLOOR |
| yazar | CD4-positive, alpha-beta cytotoxic T cell | raw | OLD_vs_MID | 78 | 13 | 0.25544 | 0.23504 | 0.020398 | 0.00049975 | 0.001999 | 0.0023947 | 0.015858 | 0.022614 | EVALUATED |
| yazar | CD8-positive, alpha-beta T cell | batch_residualized | MID_vs_MATURE | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | CD8-positive, alpha-beta T cell | batch_residualized | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | CD8-positive, alpha-beta T cell | raw | MID_vs_MATURE | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | CD8-positive, alpha-beta T cell | raw | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | central memory CD4-positive, alpha-beta T cell | batch_residualized | MID_vs_MATURE | 191 | 58 | 0.14927 | 0.16398 | -0.01471 | 0.00049975 | 0.0023322 | 0.00058887 | 0.010271 | -0.010926 | EVALUATED |
| yazar | central memory CD4-positive, alpha-beta T cell | batch_residualized | OLD_vs_MID | 660 | 191 | 0.1381 | 0.14927 | -0.011168 | 0.00049975 | 0.0023322 | 0.00024843 | 0.0052409 | -0.0078439 | EVALUATED |
| yazar | central memory CD4-positive, alpha-beta T cell | raw | MID_vs_MATURE | 191 | 58 | 0.17691 | 0.19264 | -0.015735 | 0.00049975 | 0.001999 | 0.00083105 | 0.011059 | -0.012165 | EVALUATED |
| yazar | central memory CD4-positive, alpha-beta T cell | raw | OLD_vs_MID | 660 | 191 | 0.16629 | 0.17691 | -0.010617 | 0.00049975 | 0.001999 | 0.00045711 | 0.0056882 | -0.0072683 | EVALUATED |
| yazar | central memory CD8-positive, alpha-beta T cell | batch_residualized | MID_vs_MATURE | 22 | 4 | 0.15255 | 0.1174 | 0.035149 | NA | NA | 0.0035693 | NA | 0.026533 | UNTESTABLE_BIN_FLOOR |
| yazar | central memory CD8-positive, alpha-beta T cell | batch_residualized | OLD_vs_MID | 80 | 22 | 0.17064 | 0.15255 | 0.018086 | 0.012494 | 0.034983 | 0.0011725 | 0.020385 | 0.014597 | EVALUATED |
| yazar | central memory CD8-positive, alpha-beta T cell | raw | MID_vs_MATURE | 22 | 4 | 0.25002 | 0.27771 | -0.027682 | NA | NA | 0.0069254 | NA | 0.0017542 | UNTESTABLE_BIN_FLOOR |
| yazar | central memory CD8-positive, alpha-beta T cell | raw | OLD_vs_MID | 80 | 22 | 0.25596 | 0.25002 | 0.0059398 | 0.14443 | 0.18382 | 0.0013751 | 0.011306 | 0.0077182 | EVALUATED |
| yazar | conventional dendritic cell | batch_residualized | MID_vs_MATURE | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | conventional dendritic cell | batch_residualized | OLD_vs_MID | 5 | 0 | 0.079515 | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | conventional dendritic cell | raw | MID_vs_MATURE | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | conventional dendritic cell | raw | OLD_vs_MID | 5 | 0 | 0.21098 | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | dendritic cell | batch_residualized | MID_vs_MATURE | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | dendritic cell | batch_residualized | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | dendritic cell | raw | MID_vs_MATURE | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | dendritic cell | raw | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | double negative thymocyte | batch_residualized | MID_vs_MATURE | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | double negative thymocyte | batch_residualized | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | double negative thymocyte | raw | MID_vs_MATURE | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | double negative thymocyte | raw | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | effector memory CD4-positive, alpha-beta T cell | batch_residualized | MID_vs_MATURE | 50 | 7 | 0.20026 | 0.19382 | 0.0064404 | 0.34583 | 0.4611 | 0.0033097 | 0.019414 | 0.014628 | EVALUATED |
| yazar | effector memory CD4-positive, alpha-beta T cell | batch_residualized | OLD_vs_MID | 201 | 50 | 0.20263 | 0.20026 | 0.0023652 | 0.32134 | 0.4611 | 0.0010512 | 0.0067883 | 0.0022923 | EVALUATED |
| yazar | effector memory CD4-positive, alpha-beta T cell | raw | MID_vs_MATURE | 50 | 7 | 0.25655 | 0.27726 | -0.020701 | 0.0009995 | 0.0034983 | 0.0076087 | 0.018634 | -0.003149 | EVALUATED |
| yazar | effector memory CD4-positive, alpha-beta T cell | raw | OLD_vs_MID | 201 | 50 | 0.26062 | 0.25655 | 0.004068 | 0.12994 | 0.17325 | 0.0022549 | 0.0074667 | 0.0049896 | EVALUATED |
| yazar | effector memory CD8-positive, alpha-beta T cell | batch_residualized | MID_vs_MATURE | 136 | 34 | 0.18337 | 0.18741 | -0.0040462 | 0.33583 | 0.4611 | 0.0012859 | 0.011899 | -0.0009416 | EVALUATED |
| yazar | effector memory CD8-positive, alpha-beta T cell | batch_residualized | OLD_vs_MID | 552 | 136 | 0.171 | 0.18337 | -0.012371 | 0.00049975 | 0.0023322 | 0.00097453 | 0.0060846 | -0.0091138 | EVALUATED |
| yazar | effector memory CD8-positive, alpha-beta T cell | raw | MID_vs_MATURE | 136 | 34 | 0.21586 | 0.22876 | -0.012898 | 0.012994 | 0.030318 | 0.0023743 | 0.014462 | -0.0083107 | EVALUATED |
| yazar | effector memory CD8-positive, alpha-beta T cell | raw | OLD_vs_MID | 552 | 136 | 0.20342 | 0.21586 | -0.012439 | 0.00049975 | 0.001999 | 0.0013821 | 0.0076769 | -0.0097371 | EVALUATED |
| yazar | erythrocyte | batch_residualized | MID_vs_MATURE | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | erythrocyte | batch_residualized | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | erythrocyte | raw | MID_vs_MATURE | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | erythrocyte | raw | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | gamma-delta T cell | batch_residualized | MID_vs_MATURE | 25 | 7 | 0.18818 | 0.18817 | 5.763e-06 | 0.9995 | 0.9995 | 0.0035402 | 0.024816 | 0.0042734 | EVALUATED |
| yazar | gamma-delta T cell | batch_residualized | OLD_vs_MID | 88 | 25 | 0.18064 | 0.18818 | -0.0075385 | 0.095452 | 0.22272 | 0.0012727 | 0.012932 | -0.0046109 | EVALUATED |
| yazar | gamma-delta T cell | raw | MID_vs_MATURE | 25 | 7 | 0.27049 | 0.26516 | 0.0053339 | 0.51574 | 0.55541 | 0.0044473 | 0.022331 | 0.015474 | EVALUATED |
| yazar | gamma-delta T cell | raw | OLD_vs_MID | 88 | 25 | 0.25311 | 0.27049 | -0.017378 | 0.00049975 | 0.001999 | 0.0021939 | 0.012347 | -0.011176 | EVALUATED |
| yazar | hematopoietic precursor cell | batch_residualized | MID_vs_MATURE | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | hematopoietic precursor cell | batch_residualized | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | hematopoietic precursor cell | raw | MID_vs_MATURE | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | hematopoietic precursor cell | raw | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | innate lymphoid cell | batch_residualized | MID_vs_MATURE | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | innate lymphoid cell | batch_residualized | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | innate lymphoid cell | raw | MID_vs_MATURE | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | innate lymphoid cell | raw | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | memory B cell | batch_residualized | MID_vs_MATURE | 52 | 10 | 0.19512 | 0.18453 | 0.010598 | 0.14193 | 0.26412 | 0.0027917 | 0.020204 | 0.012034 | EVALUATED |
| yazar | memory B cell | batch_residualized | OLD_vs_MID | 176 | 52 | 0.19881 | 0.19512 | 0.0036842 | 0.18641 | 0.30702 | 0.001179 | 0.0078426 | 0.0049721 | EVALUATED |
| yazar | memory B cell | raw | MID_vs_MATURE | 52 | 10 | 0.2605 | 0.25286 | 0.0076428 | 0.23188 | 0.28229 | 0.0044061 | 0.017836 | 0.01359 | EVALUATED |
| yazar | memory B cell | raw | OLD_vs_MID | 176 | 52 | 0.2558 | 0.2605 | -0.0047046 | 0.11194 | 0.15672 | 0.0018212 | 0.0081747 | -0.0022957 | EVALUATED |
| yazar | mucosal-associated invariant T cell | batch_residualized | MID_vs_MATURE | 12 | 3 | 0.15089 | 0.17558 | -0.02469 | NA | NA | 0.0030019 | NA | -0.0030172 | UNTESTABLE_BIN_FLOOR |
| yazar | mucosal-associated invariant T cell | batch_residualized | OLD_vs_MID | 23 | 12 | 0.14197 | 0.15089 | -0.0089224 | 0.57071 | 0.66583 | 0.0014008 | 0.040871 | -0.0068746 | EVALUATED |
| yazar | mucosal-associated invariant T cell | raw | MID_vs_MATURE | 12 | 3 | 0.26961 | 0.27755 | -0.0079331 | NA | NA | 0.0053815 | NA | 0.024072 | UNTESTABLE_BIN_FLOOR |
| yazar | mucosal-associated invariant T cell | raw | OLD_vs_MID | 23 | 12 | 0.26379 | 0.26961 | -0.0058203 | 0.24938 | 0.29094 | 0.0029481 | 0.014177 | 8.5734e-05 | EVALUATED |
| yazar | naive B cell | batch_residualized | MID_vs_MATURE | 50 | 11 | 0.19755 | 0.19975 | -0.0021971 | 0.68716 | 0.74001 | 0.0015033 | 0.015716 | 0.0022404 | EVALUATED |
| yazar | naive B cell | batch_residualized | OLD_vs_MID | 251 | 50 | 0.19906 | 0.19755 | 0.0015118 | 0.6027 | 0.67502 | 0.00030343 | 0.0079275 | 0.0016794 | EVALUATED |
| yazar | naive B cell | raw | MID_vs_MATURE | 50 | 11 | 0.23854 | 0.24252 | -0.0039783 | 0.57821 | 0.59963 | 0.0019962 | 0.019684 | 0.0039198 | EVALUATED |
| yazar | naive B cell | raw | OLD_vs_MID | 251 | 50 | 0.24455 | 0.23854 | 0.0060134 | 0.10545 | 0.1554 | 0.00058209 | 0.010324 | 0.0053411 | EVALUATED |
| yazar | naive thymus-derived CD4-positive, alpha-beta T cell | batch_residualized | MID_vs_MATURE | 163 | 51 | 0.16645 | 0.17284 | -0.0063907 | 0.11844 | 0.23788 | 0.00075692 | 0.011482 | -0.0033603 | EVALUATED |
| yazar | naive thymus-derived CD4-positive, alpha-beta T cell | batch_residualized | OLD_vs_MID | 564 | 163 | 0.15443 | 0.16645 | -0.01202 | 0.00049975 | 0.0023322 | 0.00048152 | 0.0067052 | -0.0089857 | EVALUATED |
| yazar | naive thymus-derived CD4-positive, alpha-beta T cell | raw | MID_vs_MATURE | 163 | 51 | 0.19729 | 0.20962 | -0.012329 | 0.0089955 | 0.022898 | 0.0007797 | 0.01285 | -0.0083287 | EVALUATED |
| yazar | naive thymus-derived CD4-positive, alpha-beta T cell | raw | OLD_vs_MID | 564 | 163 | 0.18588 | 0.19729 | -0.011412 | 0.00049975 | 0.001999 | 0.00074117 | 0.0079 | -0.0090323 | EVALUATED |
| yazar | naive thymus-derived CD8-positive, alpha-beta T cell | batch_residualized | MID_vs_MATURE | 102 | 40 | 0.18129 | 0.17539 | 0.0059019 | 0.11894 | 0.23788 | 0.00098061 | 0.010565 | 0.0050545 | EVALUATED |
| yazar | naive thymus-derived CD8-positive, alpha-beta T cell | batch_residualized | OLD_vs_MID | 210 | 102 | 0.18787 | 0.18129 | 0.0065744 | 0.003998 | 0.013993 | 0.0003633 | 0.0062084 | 0.0056624 | EVALUATED |
| yazar | naive thymus-derived CD8-positive, alpha-beta T cell | raw | MID_vs_MATURE | 102 | 40 | 0.2293 | 0.22159 | 0.0077018 | 0.09995 | 0.1554 | 0.0019524 | 0.013138 | 0.0075608 | EVALUATED |
| yazar | naive thymus-derived CD8-positive, alpha-beta T cell | raw | OLD_vs_MID | 210 | 102 | 0.23618 | 0.2293 | 0.0068801 | 0.0089955 | 0.022898 | 0.00057516 | 0.0075701 | 0.0065099 | EVALUATED |
| yazar | natural killer cell | batch_residualized | MID_vs_MATURE | 148 | 38 | 0.18156 | 0.19671 | -0.015149 | 0.00049975 | 0.0023322 | 0.002419 | 0.012099 | -0.010797 | EVALUATED |
| yazar | natural killer cell | batch_residualized | OLD_vs_MID | 598 | 148 | 0.16495 | 0.18156 | -0.016611 | 0.00049975 | 0.0023322 | 0.00035989 | 0.0059422 | -0.011854 | EVALUATED |
| yazar | natural killer cell | raw | MID_vs_MATURE | 148 | 38 | 0.21495 | 0.23055 | -0.015606 | 0.001999 | 0.0062191 | 0.0026996 | 0.014015 | -0.01127 | EVALUATED |
| yazar | natural killer cell | raw | OLD_vs_MID | 598 | 148 | 0.19435 | 0.21495 | -0.020596 | 0.00049975 | 0.001999 | 0.00058316 | 0.0070215 | -0.01541 | EVALUATED |
| yazar | plasmablast | batch_residualized | MID_vs_MATURE | 5 | 1 | 0.052775 | NA | NA | NA | NA | NA | NA | 0.042935 | UNTESTABLE_BIN_FLOOR |
| yazar | plasmablast | batch_residualized | OLD_vs_MID | 8 | 5 | 0.10698 | 0.052775 | 0.054201 | NA | NA | 0.020478 | NA | 0.021058 | UNTESTABLE_BIN_FLOOR |
| yazar | plasmablast | raw | MID_vs_MATURE | 5 | 1 | 0.26659 | NA | NA | NA | NA | NA | NA | 0.17595 | UNTESTABLE_BIN_FLOOR |
| yazar | plasmablast | raw | OLD_vs_MID | 8 | 5 | 0.27261 | 0.26659 | 0.0060155 | NA | NA | 0.024977 | NA | 0.010295 | UNTESTABLE_BIN_FLOOR |
| yazar | plasmacytoid dendritic cell | batch_residualized | MID_vs_MATURE | 1 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | plasmacytoid dendritic cell | batch_residualized | OLD_vs_MID | 2 | 1 | 1.0537e-08 | NA | NA | NA | NA | NA | NA | 0 | UNTESTABLE_BIN_FLOOR |
| yazar | plasmacytoid dendritic cell | raw | MID_vs_MATURE | 1 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | plasmacytoid dendritic cell | raw | OLD_vs_MID | 2 | 1 | 0.25458 | NA | NA | NA | NA | NA | NA | 0.13387 | UNTESTABLE_BIN_FLOOR |
| yazar | platelet | batch_residualized | MID_vs_MATURE | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | platelet | batch_residualized | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | platelet | raw | MID_vs_MATURE | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | platelet | raw | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | regulatory T cell | batch_residualized | MID_vs_MATURE | 10 | 0 | 0.19965 | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | regulatory T cell | batch_residualized | OLD_vs_MID | 59 | 10 | 0.16894 | 0.19965 | -0.030718 | 0.011994 | 0.034983 | 0.0023384 | 0.031629 | -0.022551 | EVALUATED |
| yazar | regulatory T cell | raw | MID_vs_MATURE | 10 | 0 | 0.28456 | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | regulatory T cell | raw | OLD_vs_MID | 59 | 10 | 0.27515 | 0.28456 | -0.009413 | 0.029985 | 0.05997 | 0.0051824 | 0.01205 | 0.0024373 | EVALUATED |
| yazar | transitional stage B cell | batch_residualized | MID_vs_MATURE | 39 | 9 | 0.19092 | 0.19696 | -0.0060461 | 0.46127 | 0.58707 | 0.0020953 | 0.022064 | -0.0019561 | EVALUATED |
| yazar | transitional stage B cell | batch_residualized | OLD_vs_MID | 147 | 39 | 0.19503 | 0.19092 | 0.0041158 | 0.22889 | 0.35604 | 0.0022335 | 0.0095023 | 0.0055596 | EVALUATED |
| yazar | transitional stage B cell | raw | MID_vs_MATURE | 39 | 9 | 0.25873 | 0.26439 | -0.0056594 | 0.36282 | 0.40636 | 0.0032928 | 0.016968 | 0.003841 | EVALUATED |
| yazar | transitional stage B cell | raw | OLD_vs_MID | 147 | 39 | 0.25884 | 0.25873 | 0.00011347 | 0.98051 | 0.98051 | 0.0024061 | 0.0097255 | 0.0010453 | EVALUATED |

### Verdicts

- Gao: INCONCLUSIVE. PROVES 0, DOES NOT PROVE 0, INCONCLUSIVE 9, UNTESTABLE 4 cell types.
  Inconclusive units: Cardiomyocyte, Endocardial, Endothelial, Fibroblast, Lymphoid, Myeloid, Neuronal, Pericyte, vSMC.
  Untestable units: Adipocyte, Epicardial, LEC, Mast.
- Yazar: PROVES. PROVES 2, DOES NOT PROVE 0, INCONCLUSIVE 10, UNTESTABLE 16 cell types.
  Proven units: memory B cell, naive thymus-derived CD8-positive, alpha-beta T cell.
  Inconclusive units: CD14-low, CD16-positive monocyte, CD14-positive monocyte, central memory CD4-positive, alpha-beta T cell, effector memory CD4-positive, alpha-beta T cell, effector memory CD8-positive, alpha-beta T cell, gamma-delta T cell, naive B cell, naive thymus-derived CD4-positive, alpha-beta T cell, natural killer cell, transitional stage B cell.
  Untestable units: CD16-negative, CD56-bright natural killer cell, human, CD4-positive, alpha-beta T cell, CD4-positive, alpha-beta cytotoxic T cell, CD8-positive, alpha-beta T cell, central memory CD8-positive, alpha-beta T cell, conventional dendritic cell, dendritic cell, double negative thymocyte, erythrocyte, hematopoietic precursor cell, innate lymphoid cell, mucosal-associated invariant T cell, plasmablast, plasmacytoid dendritic cell, platelet, regulatory T cell.

### Distinction from the Nov 17 distance-variance test

The Nov 17 test reduced each donor to one distance to a young reference and then compared the variance of those donor distances. Block 3 instead measures dispersion within each age bin directly from pairwise distances among donor expression centroids. It does not use a young reference and does not analyze the variance of one distance per donor. These are different estimands.
- Gao: this profile-cloud test does not clearly agree or disagree with the Nov 17 batch-adjusted null because this test is inconclusive. The Nov 17 headline was Gao OLD to YOUNG distance-variance ratio 0.79, p 0.85; Yazar raw ratio 1.62, p 1.9e-07, changing to pool-residualized ratio 1.25, p 0.43.
- Yazar: this profile-cloud test does not agree with the Nov 17 batch-adjusted null, although the estimands differ. The Nov 17 headline was Gao OLD to YOUNG distance-variance ratio 0.79, p 0.85; Yazar raw ratio 1.62, p 1.9e-07, changing to pool-residualized ratio 1.25, p 0.43.

### Reproducibility and family sizes

Permutation base seed 20261118; split-half seeds 20261116, 20261117, 20261118, 20261119, 20261120. Exact comparison seeds are in the CSV. All 2000 permutations resample donor labels, never cells.
Evaluated BH family sizes: gao batch_residualized mid_secondary m=22; gao batch_residualized old_mature m=9; gao raw mid_secondary m=22; gao raw old_mature m=9; yazar batch_residualized mid_secondary m=28; yazar batch_residualized old_mature m=12; yazar raw mid_secondary m=28; yazar raw old_mature m=12.
No primary cells were downsampled. Full-depth within-cell proportions use only cells passing the 1500-count eligibility filter.

### Comparison against the original Block 3

Original values are the batch-residualized OLD vs YOUNG results from reports/followup-11pt-block3-profile-cloud.md. Mature values are the batch-residualized OLD vs MATURE results from the CSV in this run. N YOUNG is the original baseline size; N MATURE is the mature baseline size. Verdicts use the frozen rules in each run.

| Dataset | Cell type | Original diff (p / q) | Original N YOUNG | Original verdict | Mature diff (p / q) | Mature N MATURE | Mature verdict | Change |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gao | Adipocyte | NA (untestable) | 5 | UNTESTABLE | 0.022685 (NA / NA) | 3 | UNTESTABLE | none |
| gao | Cardiomyocyte | -0.01441 (0.05997 / 0.29985) | 18 | INCONCLUSIVE | -0.018097 (0.04848 / 0.2181) | 11 | INCONCLUSIVE | none |
| gao | Endocardial | -0.0048909 (0.68916 / 0.80293) | 14 | INCONCLUSIVE | -0.0096921 (0.5297 / 0.986) | 8 | INCONCLUSIVE | none |
| gao | Endothelial | 0.012335 (0.1909 / 0.47726) | 17 | INCONCLUSIVE | -0.00083192 (0.9325 / 0.986) | 10 | INCONCLUSIVE | sign flips, still null |
| gao | Epicardial | NA (untestable) | 4 | UNTESTABLE | -0.062949 (NA / NA) | 2 | UNTESTABLE | none |
| gao | Fibroblast | -0.0045053 (0.61869 / 0.80293) | 19 | INCONCLUSIVE | -0.0030994 (0.7916 / 0.986) | 12 | INCONCLUSIVE | none |
| gao | LEC | -0.00088755 (0.94553 / 0.94553) | 6 | INCONCLUSIVE | -0.031263 (NA / NA) | 4 | UNTESTABLE | lost to bin floor (MATURE n 4) |
| gao | Lymphoid | -0.041819 (0.0054973 / 0.054973) | 12 | INCONCLUSIVE | -0.043048 (0.007496 / 0.06747) | 9 | INCONCLUSIVE | none |
| gao | Mast | NA (untestable) | 2 | UNTESTABLE | NA (NA / NA) | 1 | UNTESTABLE | none |
| gao | Myeloid | 0.0041776 (0.72264 / 0.80293) | 18 | INCONCLUSIVE | -0.00022955 (0.986 / 0.986) | 11 | INCONCLUSIVE | sign flips, still null |
| gao | Neuronal | -0.015889 (0.18541 / 0.47726) | 13 | INCONCLUSIVE | -0.014183 (0.3308 / 0.986) | 9 | INCONCLUSIVE | none |
| gao | Pericyte | 0.01364 (0.24188 / 0.48376) | 18 | INCONCLUSIVE | 0.0024578 (0.8801 / 0.986) | 11 | INCONCLUSIVE | none |
| gao | vSMC | 0.0092678 (0.42479 / 0.70798) | 16 | INCONCLUSIVE | 0.0042479 (0.7636 / 0.986) | 9 | INCONCLUSIVE | none |
| yazar | CD14-low, CD16-positive monocyte | -0.0051431 (0.2094 / 0.2255) | 14 | INCONCLUSIVE | 0.00068097 (0.8911 / 0.912) | 9 | INCONCLUSIVE | sign flips, still null |
| yazar | CD14-positive monocyte | -0.016277 (0.00049975 / 0.0009995) | 40 | INCONCLUSIVE | -0.0067936 (0.05797 / 0.09938) | 21 | INCONCLUSIVE | contraction loses BH significance |
| yazar | central memory CD4-positive, alpha-beta T cell | -0.024145 (0.00049975 / 0.0009995) | 99 | INCONCLUSIVE | -0.025878 (0.0004998 / 0.001199) | 58 | INCONCLUSIVE | none, contraction survives |
| yazar | central memory CD8-positive, alpha-beta T cell | 0.015276 (0.13393 / 0.17046) | 7 | INCONCLUSIVE | 0.053235 (NA / NA) | 4 | UNTESTABLE | lost to bin floor (MATURE n 4) |
| yazar | effector memory CD4-positive, alpha-beta T cell | -0.0081061 (0.043478 / 0.06087) | 17 | INCONCLUSIVE | 0.0088056 (0.1344 / 0.2016) | 7 | INCONCLUSIVE | sign flips, still null |
| yazar | effector memory CD8-positive, alpha-beta T cell | -0.019741 (0.00049975 / 0.0009995) | 56 | INCONCLUSIVE | -0.016417 (0.0004998 / 0.001199) | 34 | INCONCLUSIVE | none, contraction survives |
| yazar | gamma-delta T cell | -0.016373 (0.0064968 / 0.011369) | 15 | INCONCLUSIVE | -0.0075328 (0.3108 / 0.4145) | 7 | INCONCLUSIVE | contraction loses BH significance |
| yazar | memory B cell | 0.0089008 (0.009995 / 0.015548) | 25 | PROVES | 0.014282 (0.01249 / 0.02499) | 10 | PROVES | survives, larger diff, smaller baseline n |
| yazar | mucosal-associated invariant T cell | 0.0047791 (0.7971 / 0.7971) | 10 | INCONCLUSIVE | -0.033612 (NA / NA) | 3 | UNTESTABLE | lost to bin floor (MATURE n 3) |
| yazar | naive B cell | 0.0052111 (0.2089 / 0.2255) | 21 | INCONCLUSIVE | -0.00068529 (0.912 / 0.912) | 11 | INCONCLUSIVE | sign flips, still null |
| yazar | naive thymus-derived CD4-positive, alpha-beta T cell | -0.018435 (0.00049975 / 0.0009995) | 89 | INCONCLUSIVE | -0.018411 (0.0004998 / 0.001199) | 51 | INCONCLUSIVE | none, contraction survives |
| yazar | naive thymus-derived CD8-positive, alpha-beta T cell | 0.021979 (0.00049975 / 0.0009995) | 77 | PROVES | 0.012476 (0.0004998 / 0.001199) | 40 | PROVES | survives, smaller diff |
| yazar | natural killer cell | -0.028636 (0.00049975 / 0.0009995) | 66 | INCONCLUSIVE | -0.03176 (0.0004998 / 0.001199) | 38 | INCONCLUSIVE | none, contraction survives |
| yazar | transitional stage B cell | 0.016709 (0.00049975 / 0.0009995) | 26 | INCONCLUSIVE | -0.0019303 (0.7491 / 0.899) | 9 | INCONCLUSIVE | positive signal disappears |

Headline: both original Yazar PROVES survive at 30-plus (memory B cell diff 0.014282, q 0.02499, floor 0.0025937, raw diff 0.0029383 same sign; naive CD8 diff 0.012476, q 0.001199, floor 0.00092446, raw diff 0.014582 same sign). Gao stays INCONCLUSIVE with no PROVES and no DOES NOT PROVE. The Yazar contraction pattern partly survives: 4 of the 6 originally BH-significant contracting types remain BH-significant (central memory CD4, effector memory CD8, naive CD4, natural killer cell); CD14-positive monocyte (q 0.09938) and gamma-delta T cell (q 0.4145) lose BH significance. Power loss from restriction: Gao LEC moves to UNTESTABLE (MATURE n 4) and Yazar central memory CD8 (n 4) and mucosal-associated invariant T cell (n 3) move to UNTESTABLE; mature MDE80 values are 0.027 to 0.053 in Gao and 0.0086 to 0.022 in Yazar for evaluated types, all above their split-half noise floors, so evaluated nulls remain INCONCLUSIVE rather than DOES NOT PROVE under the frozen rules. OLD vs MID is retained as a distinct mature contrast in the MID secondary table and is not overwritten by the new primary.
