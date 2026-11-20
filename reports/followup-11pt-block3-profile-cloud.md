# Block 3: Full profile-cloud donor dispersion in Gao heart and Yazar blood

Status: complete. This specification was written before running the analysis.

## Exact analysis specification

- Per dataset, cell type, and donor, using the shared eligibility floors: donor centroid = mean within-cell proportions across that donor's cells at full depth, over the technical-masked gene set.
- Primary outcome, per dataset, cell type, and age bin, with LT20 excluded: mean pairwise Hellinger distance among donor centroids WITHIN the bin. Secondary outcome: mean distance of each centroid to its own bin centroid.
- Primary comparison: OLD vs YOUNG mean pairwise distance, with a donor-permutation null that permutes donor age-bin labels within dataset and cell type, uses 2000 permutations and a fixed seed, and recomputes the within-bin means each time. Report observed difference, permutation p-value, N donors per bin, and MDE80.
- Batch handling has two arms: (a) raw centroids and (b) batch-residualized centroids, residualizing each gene's centroid proportions on batch across donors, using Gao study indicators and Yazar pool means, then re-normalizing to proportions before computing distances. The batch-residualized arm is primary and the raw arm is the sensitivity.
- MID comparisons are secondary: OLD vs MID and MID vs YOUNG use the same permutation procedure and the same arms.
- Split-half technical null: split each donor's cells into two halves with 5 fixed seeds, build half-centroids, compute mean pairwise distances among half-centroids within bins, and report the noise floor. A claimed age difference must exceed it.
- Family, stated before correction: BH q within each (dataset, arm) slice across cell types for the OLD vs YOUNG comparison. MID comparisons are a separate stated family.

## Fixed implementation details

- Inputs are the two Block 1 h5ad files and the donor tables `results/followup-between-subject-gao-donors.csv` and `results/followup-between-subject-yazar-donors.csv`. Yazar excludes rows labeled peripheral blood mononuclear cell. No new downloads are used.
- The human technical gene mask excludes uppercase MT-, RPL, RPS, MRPL, MRPS, XIST, TSIX, DDX3Y, EIF2S3Y, UTY, and KDM5D before depth calculation and normalization. ERCC controls are also excluded, consistent with the existing donor tables. All remaining genes are retained without a detection filter.
- Full-depth raw counts are converted to within-cell proportions. Cells with at least 1500 post-mask counts are eligible, with no downsampling. A donor-cell-type unit requires at least 20 eligible cells. Each comparison requires at least 6 eligible donors in each bin. Units below these floors are UNTESTABLE, never nulls.
- Age bins reuse `assign_age_bin` from `scripts/run_followup_between_subject_gao_yazar.py`: LT20 under 20, YOUNG 20 to 39, MID 40 to 59, OLD 60 and older. LT20 is excluded from all Block 3 calculations.
- Each cell contributes equally to its donor centroid. This is the arithmetic mean of within-cell proportion vectors, not a sum-of-counts pseudobulk. Donor is the unit for every comparison, permutation, and noise-floor summary.
- Hellinger distance uses the formula in `hellinger_distance` from `scripts/run_h3_expression_vector.py`, `sqrt(1 - sum_g sqrt(p_g q_g))`, on proportion vectors. The implementation uses an equivalent vectorized square-root affinity matrix for all pairs and validates it against the imported function. Distances from donors to their own bin centroid call the imported function directly. The bin centroid is the arithmetic mean of donor centroids in that bin.
- Batch residualization is fitted once per cell type and arm input across all eligible YOUNG, MID, and OLD donors. The design is intercept plus study indicators for Gao and intercept plus pool indicators for Yazar, equivalent to removing batch means and restoring the global gene mean. Residualized values are clipped at zero and each donor profile is divided by its sum to obtain nonnegative proportions. The same procedure is applied independently to each split half.
- Comparisons use the first bin minus the second bin: OLD minus YOUNG, OLD minus MID, and MID minus YOUNG. For each comparison, only donors in those two bins enter the permutation. Each of 2000 permutations holds the two bin sizes fixed and recomputes both within-bin mean pairwise distances from the fixed donor distance matrix. P-values are two-sided with the plus-one correction. Permutation base seed is 20261118, with deterministic dataset, cell-type, arm, and comparison offsets recorded in the CSV.
- The secondary own-bin-centroid outcome is reported as observed bin means and observed comparison differences. The permutation p-value and BH q-value apply to the primary mean-pairwise outcome only.
- Split-half seeds are 20261116, 20261117, 20261118, 20261119, and 20261120. For each donor and seed, cells are randomly permuted and divided into halves of floor(N/2) and ceil(N/2). For each arm and comparison, the noise floor is the mean over 5 seeds of the absolute half-to-half change in that comparison's mean-pairwise difference. Bin-specific absolute half-to-half changes and the corresponding secondary-outcome noise floor are also recorded. This directly estimates technical instability of the age contrast.
- MDE80 is `(z_0.975 + z_0.80) * SD(permutation null)`, equal to 2.801585 times the permutation-null standard deviation. This is the smallest absolute difference with approximately 80 percent power for a two-sided 0.05 test under the donor-permutation scale. For a null or negative batch-residualized result, MDE80 at or below its split-half noise floor is pre-stated as small enough to exclude an increase detectable above technical noise; a larger MDE80 is underpowered.
- The OLD vs YOUNG BH families contain evaluated cell types within each dataset and arm. The separate MID family contains both OLD vs MID and MID vs YOUNG evaluated rows within each dataset and arm. UNTESTABLE rows are excluded and have missing q-values. Family sizes are reported after evaluation.
- No 1500 rarefaction is performed in Block 3 because it is not part of this block's specification. The 1500 threshold is only a cell eligibility filter.

## Pre-stated verdict rules

For each cell type, PROVES that the donor cloud spreads with age requires a positive batch-residualized OLD vs YOUNG difference, BH q < 0.05, and a difference greater than the split-half noise floor. DOES NOT PROVE applies when the batch-residualized result is null or negative and MDE80 is at or below the noise floor. INCONCLUSIVE applies to a null with MDE80 above the noise floor, a positive significant result that does not exceed the noise floor, or opposite signs in the raw and batch-residualized arms. UNTESTABLE units remain UNTESTABLE.

For each dataset summary, PROVES requires at least one PROVES cell type. DOES NOT PROVE requires at least one evaluated cell type and every evaluated cell type to meet DOES NOT PROVE. All other dataset summaries are INCONCLUSIVE. Gao and Yazar are reported separately and are never pooled.

## Results

All pairwise effects are first-bin minus second-bin Hellinger distance differences. The batch-residualized arm is primary; raw is the sensitivity. NA denotes an UNTESTABLE or not applicable field.

### OLD vs YOUNG primary results and raw sensitivity

| Dataset | Cell type | N YOUNG | N OLD | Residualized YOUNG | Residualized OLD | Residualized difference | p | BH q | Noise floor | MDE80 | Raw difference | Raw q | Raw noise | Raw MDE80 | Own-bin secondary difference | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gao | Adipocyte | 5 | 12 | 0.25282 | 0.26096 | 0.0081446 | NA | NA | 0.005445 | NA | -0.049363 | NA | 0.0042088 | NA | 0.015901 | UNTESTABLE |
| gao | Cardiomyocyte | 18 | 39 | 0.16651 | 0.1521 | -0.01441 | 0.05997 | 0.29985 | 0.00023679 | 0.021956 | -0.00095368 | 0.90305 | 0.00039172 | 0.021214 | -0.0085518 | INCONCLUSIVE |
| gao | Endocardial | 14 | 29 | 0.23 | 0.22511 | -0.0048909 | 0.68916 | 0.80293 | 0.0033939 | 0.035148 | 0.041355 | 0.0049975 | 0.0039304 | 0.037117 | -0.0017681 | INCONCLUSIVE |
| gao | Endothelial | 17 | 30 | 0.19293 | 0.20527 | 0.012335 | 0.1909 | 0.47726 | 0.0030114 | 0.026023 | 0.029671 | 0.0087456 | 0.0030552 | 0.029506 | 0.010603 | INCONCLUSIVE |
| gao | Epicardial | 4 | 8 | 0.30609 | 0.24501 | -0.06108 | NA | NA | 0.0069177 | NA | 0.024925 | NA | 0.0094391 | NA | -0.032914 | UNTESTABLE |
| gao | Fibroblast | 19 | 38 | 0.20418 | 0.19967 | -0.0045053 | 0.61869 | 0.80293 | 0.0012759 | 0.024917 | 0.021917 | 0.044264 | 0.00051446 | 0.030039 | -0.001902 | INCONCLUSIVE |
| gao | LEC | 6 | 13 | 0.25587 | 0.25498 | -0.00088755 | 0.94553 | 0.94553 | 0.0037923 | 0.036247 | 0.028419 | 0.13681 | 0.010315 | 0.050709 | 0.0089901 | INCONCLUSIVE |
| gao | Lymphoid | 12 | 22 | 0.3103 | 0.26848 | -0.041819 | 0.0054973 | 0.054973 | 0.0048245 | 0.046018 | -0.0091708 | 0.58082 | 0.0046234 | 0.039214 | -0.026416 | INCONCLUSIVE |
| gao | Mast | 2 | 8 | 0.22659 | 0.26205 | 0.035468 | NA | NA | 0.0033288 | NA | 0.074183 | NA | 0.014593 | NA | 0.059939 | UNTESTABLE |
| gao | Myeloid | 18 | 39 | 0.21799 | 0.22217 | 0.0041776 | 0.72264 | 0.80293 | 0.0018116 | 0.033223 | 0.064204 | 0.0024988 | 0.001253 | 0.039142 | 0.0022148 | INCONCLUSIVE |
| gao | Neuronal | 13 | 19 | 0.2985 | 0.28261 | -0.015889 | 0.18541 | 0.47726 | 0.0054052 | 0.03287 | 0.02208 | 0.021989 | 0.0029296 | 0.025774 | -0.0093782 | INCONCLUSIVE |
| gao | Pericyte | 18 | 33 | 0.22229 | 0.23593 | 0.01364 | 0.24188 | 0.48376 | 0.0030433 | 0.032736 | 0.055869 | 0.0024988 | 0.0023902 | 0.030871 | 0.010916 | INCONCLUSIVE |
| gao | vSMC | 16 | 30 | 0.25766 | 0.26692 | 0.0092678 | 0.42479 | 0.70798 | 0.0050227 | 0.030917 | 0.026601 | 0.043312 | 0.0053941 | 0.034718 | 0.0097026 | INCONCLUSIVE |
| yazar | CD14-low, CD16-positive monocyte | 14 | 204 | 0.17067 | 0.16553 | -0.0051431 | 0.2094 | 0.2255 | 0.0021399 | 0.011669 | -0.0096247 | 0.033074 | 0.0035982 | 0.011804 | 3.0941e-05 | INCONCLUSIVE |
| yazar | CD14-positive monocyte | 40 | 299 | 0.18776 | 0.17149 | -0.016277 | 0.00049975 | 0.0009995 | 0.0018016 | 0.007096 | -0.02265 | 0.00087456 | 0.0019153 | 0.009262 | -0.010311 | INCONCLUSIVE |
| yazar | CD16-negative, CD56-bright natural killer cell, human | 1 | 4 | NA | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | 6.3951e-10 | UNTESTABLE |
| yazar | CD4-positive, alpha-beta T cell | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE |
| yazar | CD4-positive, alpha-beta cytotoxic T cell | 2 | 78 | 0.17986 | 0.16289 | -0.016974 | NA | NA | 0.005967 | NA | -0.013994 | NA | 0.010342 | NA | 0.016799 | UNTESTABLE |
| yazar | CD8-positive, alpha-beta T cell | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE |
| yazar | central memory CD4-positive, alpha-beta T cell | 99 | 660 | 0.16254 | 0.1384 | -0.024145 | 0.00049975 | 0.0009995 | 0.00030711 | 0.0068961 | -0.024848 | 0.00087456 | 0.00049973 | 0.007713 | -0.017493 | INCONCLUSIVE |
| yazar | central memory CD8-positive, alpha-beta T cell | 7 | 80 | 0.15725 | 0.17252 | 0.015276 | 0.13393 | 0.17046 | 0.0037066 | 0.03015 | -0.028242 | 0.00087456 | 0.0060473 | 0.018204 | 0.015676 | INCONCLUSIVE |
| yazar | conventional dendritic cell | 0 | 5 | NA | 0.079515 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE |
| yazar | dendritic cell | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE |
| yazar | double negative thymocyte | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE |
| yazar | effector memory CD4-positive, alpha-beta T cell | 17 | 201 | 0.21128 | 0.20317 | -0.0081061 | 0.043478 | 0.06087 | 0.0019504 | 0.011309 | -0.015761 | 0.00087456 | 0.0026959 | 0.011496 | -0.0021322 | INCONCLUSIVE |
| yazar | effector memory CD8-positive, alpha-beta T cell | 56 | 552 | 0.19114 | 0.1714 | -0.019741 | 0.00049975 | 0.0009995 | 0.00096094 | 0.0087064 | -0.028007 | 0.00087456 | 0.0016089 | 0.011391 | -0.013695 | INCONCLUSIVE |
| yazar | erythrocyte | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE |
| yazar | gamma-delta T cell | 15 | 88 | 0.19862 | 0.18225 | -0.016373 | 0.0064968 | 0.011369 | 0.002356 | 0.015379 | -0.01403 | 0.017491 | 0.003553 | 0.015946 | -0.011513 | INCONCLUSIVE |
| yazar | hematopoietic precursor cell | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE |
| yazar | innate lymphoid cell | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE |
| yazar | memory B cell | 25 | 176 | 0.19082 | 0.19972 | 0.0089008 | 0.009995 | 0.015548 | 0.001778 | 0.0098053 | 0.01182 | 0.0054417 | 0.0029861 | 0.010948 | 0.0084528 | PROVES |
| yazar | mucosal-associated invariant T cell | 10 | 23 | 0.13739 | 0.14217 | 0.0047791 | 0.7971 | 0.7971 | 0.0019576 | 0.044673 | -0.0084854 | 0.10495 | 0.003094 | 0.013989 | 0.0032067 | INCONCLUSIVE |
| yazar | naive B cell | 21 | 251 | 0.19435 | 0.19956 | 0.0052111 | 0.2089 | 0.2255 | 0.0014887 | 0.01139 | 0.0063761 | 0.25087 | 0.0018795 | 0.015479 | 0.0054435 | INCONCLUSIVE |
| yazar | naive thymus-derived CD4-positive, alpha-beta T cell | 89 | 564 | 0.17329 | 0.15486 | -0.018435 | 0.00049975 | 0.0009995 | 0.00057953 | 0.0085048 | -0.021264 | 0.00087456 | 0.0005807 | 0.010161 | -0.013344 | INCONCLUSIVE |
| yazar | naive thymus-derived CD8-positive, alpha-beta T cell | 77 | 210 | 0.16748 | 0.18946 | 0.021979 | 0.00049975 | 0.0009995 | 0.00084771 | 0.0070891 | 0.019777 | 0.00087456 | 0.0014574 | 0.0086311 | 0.018381 | PROVES |
| yazar | natural killer cell | 66 | 598 | 0.1939 | 0.16526 | -0.028636 | 0.00049975 | 0.0009995 | 0.0013637 | 0.0082159 | -0.031939 | 0.00087456 | 0.0017223 | 0.010063 | -0.020948 | INCONCLUSIVE |
| yazar | plasmablast | 1 | 8 | NA | 0.10698 | NA | NA | NA | NA | NA | NA | NA | NA | NA | 0.063994 | UNTESTABLE |
| yazar | plasmacytoid dendritic cell | 0 | 2 | NA | 1.0537e-08 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE |
| yazar | platelet | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE |
| yazar | regulatory T cell | 0 | 59 | NA | 0.16894 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE |
| yazar | transitional stage B cell | 26 | 147 | 0.18016 | 0.19687 | 0.016709 | 0.00049975 | 0.0009995 | 0.0011755 | 0.010807 | -0.005272 | 0.23627 | 0.0016438 | 0.011871 | 0.016393 | INCONCLUSIVE |

### MID secondary comparisons

The MID BH family combines OLD vs MID and MID vs YOUNG within each dataset and arm. Permutation p-values and q-values apply to mean pairwise distance, not the descriptive own-bin-centroid outcome.

| Dataset | Cell type | Arm | Comparison | N bin 1 | N bin 2 | Bin 1 pairwise | Bin 2 pairwise | Difference | p | BH q | Noise floor | MDE80 | Own-bin secondary difference | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gao | Adipocyte | batch_residualized | MID_vs_YOUNG | 20 | 5 | 0.24888 | 0.25282 | -0.0039389 | NA | NA | 0.0034983 | NA | 0.011059 | UNTESTABLE_BIN_FLOOR |
| gao | Adipocyte | batch_residualized | OLD_vs_MID | 12 | 20 | 0.26096 | 0.24888 | 0.012083 | 0.28686 | 0.43985 | 0.0046943 | 0.031715 | 0.0048415 | EVALUATED |
| gao | Adipocyte | raw | MID_vs_YOUNG | 20 | 5 | 0.29818 | 0.38934 | -0.091159 | NA | NA | 0.006244 | NA | -0.043769 | UNTESTABLE_BIN_FLOOR |
| gao | Adipocyte | raw | OLD_vs_MID | 12 | 20 | 0.33998 | 0.29818 | 0.041796 | 0.0024988 | 0.0095785 | 0.0056442 | 0.041425 | 0.024567 | EVALUATED |
| gao | Cardiomyocyte | batch_residualized | MID_vs_YOUNG | 52 | 18 | 0.15395 | 0.16651 | -0.012558 | 0.069965 | 0.20115 | 0.00036479 | 0.019399 | -0.0069245 | EVALUATED |
| gao | Cardiomyocyte | batch_residualized | OLD_vs_MID | 39 | 52 | 0.1521 | 0.15395 | -0.0018516 | 0.72764 | 0.76071 | 0.00043055 | 0.014181 | -0.0016273 | EVALUATED |
| gao | Cardiomyocyte | raw | MID_vs_YOUNG | 52 | 18 | 0.19619 | 0.1898 | 0.0063903 | 0.50375 | 0.6098 | 0.00046479 | 0.024874 | 0.0064949 | EVALUATED |
| gao | Cardiomyocyte | raw | OLD_vs_MID | 39 | 52 | 0.18884 | 0.19619 | -0.007344 | 0.31634 | 0.43543 | 0.00024055 | 0.019296 | -0.0050791 | EVALUATED |
| gao | Endocardial | batch_residualized | MID_vs_YOUNG | 29 | 14 | 0.25808 | 0.23 | 0.028084 | 0.010995 | 0.074713 | 0.0052993 | 0.032759 | 0.022276 | EVALUATED |
| gao | Endocardial | batch_residualized | OLD_vs_MID | 29 | 29 | 0.22511 | 0.25808 | -0.032975 | 0.00049975 | 0.011494 | 0.0048021 | 0.027815 | -0.024044 | EVALUATED |
| gao | Endocardial | raw | MID_vs_YOUNG | 29 | 14 | 0.32559 | 0.25954 | 0.066043 | 0.00049975 | 0.0022989 | 0.0048113 | 0.038925 | 0.052172 | EVALUATED |
| gao | Endocardial | raw | OLD_vs_MID | 29 | 29 | 0.3009 | 0.32559 | -0.024689 | 0.04048 | 0.084639 | 0.0034826 | 0.03381 | -0.014673 | EVALUATED |
| gao | Endothelial | batch_residualized | MID_vs_YOUNG | 40 | 17 | 0.20987 | 0.19293 | 0.016944 | 0.045477 | 0.14943 | 0.0016144 | 0.023751 | 0.013469 | EVALUATED |
| gao | Endothelial | batch_residualized | OLD_vs_MID | 30 | 40 | 0.20527 | 0.20987 | -0.0046092 | 0.53373 | 0.61379 | 0.0033993 | 0.020487 | -0.0028661 | EVALUATED |
| gao | Endothelial | raw | MID_vs_YOUNG | 40 | 17 | 0.27041 | 0.22633 | 0.044076 | 0.00049975 | 0.0022989 | 0.0016266 | 0.025244 | 0.035651 | EVALUATED |
| gao | Endothelial | raw | OLD_vs_MID | 30 | 40 | 0.256 | 0.27041 | -0.014405 | 0.066967 | 0.12835 | 0.0043212 | 0.022238 | -0.011365 | EVALUATED |
| gao | Epicardial | batch_residualized | MID_vs_YOUNG | 8 | 4 | 0.29486 | 0.30609 | -0.011227 | NA | NA | 0.0073804 | NA | 0.004698 | UNTESTABLE_BIN_FLOOR |
| gao | Epicardial | batch_residualized | OLD_vs_MID | 8 | 8 | 0.24501 | 0.29486 | -0.049853 | 0.012994 | 0.074713 | 0.0087187 | 0.060133 | -0.037612 | EVALUATED |
| gao | Epicardial | raw | MID_vs_YOUNG | 8 | 4 | 0.3316 | 0.33064 | 0.00095688 | NA | NA | 0.0074809 | NA | 0.018884 | UNTESTABLE_BIN_FLOOR |
| gao | Epicardial | raw | OLD_vs_MID | 8 | 8 | 0.35557 | 0.3316 | 0.023968 | 0.4088 | 0.52235 | 0.0097653 | 0.072717 | 0.016286 | EVALUATED |
| gao | Fibroblast | batch_residualized | MID_vs_YOUNG | 49 | 19 | 0.20991 | 0.20418 | 0.0057329 | 0.45377 | 0.5493 | 0.00098593 | 0.020251 | 0.0062763 | EVALUATED |
| gao | Fibroblast | batch_residualized | OLD_vs_MID | 38 | 49 | 0.19967 | 0.20991 | -0.010238 | 0.096952 | 0.24777 | 0.0011304 | 0.017135 | -0.0081783 | EVALUATED |
| gao | Fibroblast | raw | MID_vs_YOUNG | 49 | 19 | 0.27645 | 0.25132 | 0.025133 | 0.0069965 | 0.022989 | 0.00095227 | 0.025611 | 0.024417 | EVALUATED |
| gao | Fibroblast | raw | OLD_vs_MID | 38 | 49 | 0.27323 | 0.27645 | -0.0032162 | 0.65767 | 0.68757 | 0.0009638 | 0.019462 | -0.00095857 | EVALUATED |
| gao | LEC | batch_residualized | MID_vs_YOUNG | 17 | 6 | 0.26619 | 0.25587 | 0.010319 | 0.37931 | 0.48467 | 0.0058333 | 0.032106 | 0.018221 | EVALUATED |
| gao | LEC | batch_residualized | OLD_vs_MID | 13 | 17 | 0.25498 | 0.26619 | -0.011207 | 0.13293 | 0.28004 | 0.003683 | 0.020698 | -0.0092311 | EVALUATED |
| gao | LEC | raw | MID_vs_YOUNG | 17 | 6 | 0.3297 | 0.29583 | 0.033869 | 0.02049 | 0.047126 | 0.0067658 | 0.042286 | 0.041402 | EVALUATED |
| gao | LEC | raw | OLD_vs_MID | 13 | 17 | 0.32425 | 0.3297 | -0.0054494 | 0.6082 | 0.66612 | 0.005093 | 0.029989 | -0.0059171 | EVALUATED |
| gao | Lymphoid | batch_residualized | MID_vs_YOUNG | 25 | 12 | 0.28221 | 0.3103 | -0.028091 | 0.036982 | 0.14176 | 0.0062705 | 0.038897 | -0.015399 | EVALUATED |
| gao | Lymphoid | batch_residualized | OLD_vs_MID | 22 | 25 | 0.26848 | 0.28221 | -0.013728 | 0.18841 | 0.33333 | 0.0024097 | 0.029455 | -0.011017 | EVALUATED |
| gao | Lymphoid | raw | MID_vs_YOUNG | 25 | 12 | 0.38595 | 0.36484 | 0.021109 | 0.13043 | 0.23077 | 0.0035568 | 0.038963 | 0.022338 | EVALUATED |
| gao | Lymphoid | raw | OLD_vs_MID | 22 | 25 | 0.35567 | 0.38595 | -0.030279 | 0.007996 | 0.022989 | 0.003655 | 0.03313 | -0.021287 | EVALUATED |
| gao | Mast | batch_residualized | MID_vs_YOUNG | 8 | 2 | 0.27955 | 0.22659 | 0.052965 | NA | NA | 0.0039474 | NA | 0.066208 | UNTESTABLE_BIN_FLOOR |
| gao | Mast | batch_residualized | OLD_vs_MID | 8 | 8 | 0.26205 | 0.27955 | -0.017497 | 0.37631 | 0.48467 | 0.0024022 | 0.051939 | -0.006269 | EVALUATED |
| gao | Mast | raw | MID_vs_YOUNG | 8 | 2 | 0.43906 | 0.33247 | 0.10659 | NA | NA | 0.022007 | NA | 0.1378 | UNTESTABLE_BIN_FLOOR |
| gao | Mast | raw | OLD_vs_MID | 8 | 8 | 0.40665 | 0.43906 | -0.032409 | 0.32184 | 0.43543 | 0.0074658 | 0.090054 | -0.013794 | EVALUATED |
| gao | Myeloid | batch_residualized | MID_vs_YOUNG | 49 | 18 | 0.22983 | 0.21799 | 0.01184 | 0.16392 | 0.31418 | 0.0023193 | 0.023736 | 0.0096478 | EVALUATED |
| gao | Myeloid | batch_residualized | OLD_vs_MID | 39 | 49 | 0.22217 | 0.22983 | -0.0076622 | 0.34483 | 0.48467 | 0.0031699 | 0.022575 | -0.007433 | EVALUATED |
| gao | Myeloid | raw | MID_vs_YOUNG | 49 | 18 | 0.29647 | 0.24304 | 0.053435 | 0.00049975 | 0.0022989 | 0.0029209 | 0.030727 | 0.044656 | EVALUATED |
| gao | Myeloid | raw | OLD_vs_MID | 39 | 49 | 0.30724 | 0.29647 | 0.010768 | 0.27086 | 0.43543 | 0.0030659 | 0.026876 | 0.01266 | EVALUATED |
| gao | Neuronal | batch_residualized | MID_vs_YOUNG | 28 | 13 | 0.29332 | 0.2985 | -0.0051824 | 0.61219 | 0.6705 | 0.0050454 | 0.02798 | 0.00065282 | EVALUATED |
| gao | Neuronal | batch_residualized | OLD_vs_MID | 19 | 28 | 0.28261 | 0.29332 | -0.010707 | 0.21389 | 0.3514 | 0.0059931 | 0.023588 | -0.010031 | EVALUATED |
| gao | Neuronal | raw | MID_vs_YOUNG | 28 | 13 | 0.36224 | 0.33719 | 0.025049 | 0.0089955 | 0.022989 | 0.0027144 | 0.027519 | 0.023597 | EVALUATED |
| gao | Neuronal | raw | OLD_vs_MID | 19 | 28 | 0.35927 | 0.36224 | -0.0029686 | 0.74113 | 0.74113 | 0.0051712 | 0.024438 | -0.0021333 | EVALUATED |
| gao | Pericyte | batch_residualized | MID_vs_YOUNG | 44 | 18 | 0.2369 | 0.22229 | 0.014603 | 0.13393 | 0.28004 | 0.002171 | 0.026897 | 0.013175 | EVALUATED |
| gao | Pericyte | batch_residualized | OLD_vs_MID | 33 | 44 | 0.23593 | 0.2369 | -0.00096271 | 0.90955 | 0.90955 | 0.0028761 | 0.023327 | -0.0022587 | EVALUATED |
| gao | Pericyte | raw | MID_vs_YOUNG | 44 | 18 | 0.31105 | 0.2601 | 0.050949 | 0.00049975 | 0.0022989 | 0.002207 | 0.032728 | 0.039692 | EVALUATED |
| gao | Pericyte | raw | OLD_vs_MID | 33 | 44 | 0.31597 | 0.31105 | 0.00492 | 0.5907 | 0.66612 | 0.0028762 | 0.024814 | 0.0085876 | EVALUATED |
| gao | vSMC | batch_residualized | MID_vs_YOUNG | 40 | 16 | 0.28509 | 0.25766 | 0.027436 | 0.0094953 | 0.074713 | 0.0052339 | 0.028692 | 0.024329 | EVALUATED |
| gao | vSMC | batch_residualized | OLD_vs_MID | 30 | 40 | 0.26692 | 0.28509 | -0.018168 | 0.016992 | 0.078161 | 0.0025635 | 0.021673 | -0.014627 | EVALUATED |
| gao | vSMC | raw | MID_vs_YOUNG | 40 | 16 | 0.34087 | 0.30479 | 0.036078 | 0.00049975 | 0.0022989 | 0.0066626 | 0.02883 | 0.030069 | EVALUATED |
| gao | vSMC | raw | OLD_vs_MID | 30 | 40 | 0.33139 | 0.34087 | -0.009477 | 0.29535 | 0.43543 | 0.002486 | 0.025557 | -0.0056309 | EVALUATED |
| yazar | CD14-low, CD16-positive monocyte | batch_residualized | MID_vs_YOUNG | 44 | 14 | 0.16833 | 0.17067 | -0.0023432 | 0.67716 | 0.78134 | 0.002192 | 0.015663 | 0.00051054 | EVALUATED |
| yazar | CD14-low, CD16-positive monocyte | batch_residualized | OLD_vs_MID | 204 | 44 | 0.16553 | 0.16833 | -0.0028 | 0.27536 | 0.35917 | 0.0024205 | 0.0072176 | -0.0004796 | EVALUATED |
| yazar | CD14-low, CD16-positive monocyte | raw | MID_vs_YOUNG | 44 | 14 | 0.21571 | 0.22083 | -0.0051155 | 0.25787 | 0.29754 | 0.0026551 | 0.012532 | 0.00070082 | EVALUATED |
| yazar | CD14-low, CD16-positive monocyte | raw | OLD_vs_MID | 204 | 44 | 0.2112 | 0.21571 | -0.0045092 | 0.074463 | 0.1241 | 0.0031936 | 0.0070881 | -0.0022417 | EVALUATED |
| yazar | CD14-positive monocyte | batch_residualized | MID_vs_YOUNG | 88 | 40 | 0.17796 | 0.18776 | -0.0098018 | 0.001999 | 0.0066633 | 0.001805 | 0.0091909 | -0.0059628 | EVALUATED |
| yazar | CD14-positive monocyte | batch_residualized | OLD_vs_MID | 299 | 88 | 0.17149 | 0.17796 | -0.0064755 | 0.0014993 | 0.0056222 | 0.00081301 | 0.0052481 | -0.0043487 | EVALUATED |
| yazar | CD14-positive monocyte | raw | MID_vs_YOUNG | 88 | 40 | 0.21401 | 0.23264 | -0.018629 | 0.00049975 | 0.001363 | 0.0014301 | 0.011062 | -0.013236 | EVALUATED |
| yazar | CD14-positive monocyte | raw | OLD_vs_MID | 299 | 88 | 0.20999 | 0.21401 | -0.004021 | 0.087456 | 0.13809 | 0.00090239 | 0.0066012 | -0.0021587 | EVALUATED |
| yazar | CD16-negative, CD56-bright natural killer cell, human | batch_residualized | MID_vs_YOUNG | 1 | 1 | NA | NA | NA | NA | NA | NA | NA | 8.6034e-10 | UNTESTABLE_BIN_FLOOR |
| yazar | CD16-negative, CD56-bright natural killer cell, human | batch_residualized | OLD_vs_MID | 4 | 1 | 0 | NA | NA | NA | NA | NA | NA | -2.2082e-10 | UNTESTABLE_BIN_FLOOR |
| yazar | CD16-negative, CD56-bright natural killer cell, human | raw | MID_vs_YOUNG | 1 | 1 | NA | NA | NA | NA | NA | NA | NA | 1.0106e-07 | UNTESTABLE_BIN_FLOOR |
| yazar | CD16-negative, CD56-bright natural killer cell, human | raw | OLD_vs_MID | 4 | 1 | 0.27903 | NA | NA | NA | NA | NA | NA | 0.18207 | UNTESTABLE_BIN_FLOOR |
| yazar | CD4-positive, alpha-beta T cell | batch_residualized | MID_vs_YOUNG | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | CD4-positive, alpha-beta T cell | batch_residualized | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | CD4-positive, alpha-beta T cell | raw | MID_vs_YOUNG | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | CD4-positive, alpha-beta T cell | raw | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | CD4-positive, alpha-beta cytotoxic T cell | batch_residualized | MID_vs_YOUNG | 13 | 2 | 0.14604 | 0.17986 | -0.033821 | NA | NA | 0.0054303 | NA | 0.0041902 | UNTESTABLE_BIN_FLOOR |
| yazar | CD4-positive, alpha-beta cytotoxic T cell | batch_residualized | OLD_vs_MID | 78 | 13 | 0.16289 | 0.14604 | 0.016848 | 0.055972 | 0.099656 | 0.00075348 | 0.025225 | 0.012609 | EVALUATED |
| yazar | CD4-positive, alpha-beta cytotoxic T cell | raw | MID_vs_YOUNG | 13 | 2 | 0.23504 | 0.26943 | -0.034392 | NA | NA | 0.0092708 | NA | 0.025458 | UNTESTABLE_BIN_FLOOR |
| yazar | CD4-positive, alpha-beta cytotoxic T cell | raw | OLD_vs_MID | 78 | 13 | 0.25544 | 0.23504 | 0.020398 | 0.00049975 | 0.001363 | 0.0023947 | 0.015858 | 0.022614 | EVALUATED |
| yazar | CD8-positive, alpha-beta T cell | batch_residualized | MID_vs_YOUNG | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | CD8-positive, alpha-beta T cell | batch_residualized | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | CD8-positive, alpha-beta T cell | raw | MID_vs_YOUNG | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | CD8-positive, alpha-beta T cell | raw | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | central memory CD4-positive, alpha-beta T cell | batch_residualized | MID_vs_YOUNG | 191 | 99 | 0.14963 | 0.16254 | -0.012912 | 0.00049975 | 0.0021418 | 0.00035635 | 0.0084983 | -0.0096009 | EVALUATED |
| yazar | central memory CD4-positive, alpha-beta T cell | batch_residualized | OLD_vs_MID | 660 | 191 | 0.1384 | 0.14963 | -0.011233 | 0.00049975 | 0.0021418 | 0.00026163 | 0.0052926 | -0.0078919 | EVALUATED |
| yazar | central memory CD4-positive, alpha-beta T cell | raw | MID_vs_YOUNG | 191 | 99 | 0.17691 | 0.19114 | -0.014231 | 0.00049975 | 0.001363 | 0.00065305 | 0.0092779 | -0.011011 | EVALUATED |
| yazar | central memory CD4-positive, alpha-beta T cell | raw | OLD_vs_MID | 660 | 191 | 0.16629 | 0.17691 | -0.010617 | 0.00049975 | 0.001363 | 0.00045711 | 0.0056882 | -0.0072683 | EVALUATED |
| yazar | central memory CD8-positive, alpha-beta T cell | batch_residualized | MID_vs_YOUNG | 22 | 7 | 0.15892 | 0.15725 | 0.0016737 | 0.88406 | 0.9472 | 0.0041408 | 0.039309 | 0.0043795 | EVALUATED |
| yazar | central memory CD8-positive, alpha-beta T cell | batch_residualized | OLD_vs_MID | 80 | 22 | 0.17252 | 0.15892 | 0.013602 | 0.044478 | 0.088956 | 0.0012797 | 0.019324 | 0.011296 | EVALUATED |
| yazar | central memory CD8-positive, alpha-beta T cell | raw | MID_vs_YOUNG | 22 | 7 | 0.25002 | 0.28421 | -0.034181 | 0.00049975 | 0.001363 | 0.0071814 | 0.02051 | -0.015932 | EVALUATED |
| yazar | central memory CD8-positive, alpha-beta T cell | raw | OLD_vs_MID | 80 | 22 | 0.25596 | 0.25002 | 0.0059398 | 0.14443 | 0.18838 | 0.0013751 | 0.011306 | 0.0077182 | EVALUATED |
| yazar | conventional dendritic cell | batch_residualized | MID_vs_YOUNG | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | conventional dendritic cell | batch_residualized | OLD_vs_MID | 5 | 0 | 0.079515 | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | conventional dendritic cell | raw | MID_vs_YOUNG | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | conventional dendritic cell | raw | OLD_vs_MID | 5 | 0 | 0.21098 | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | dendritic cell | batch_residualized | MID_vs_YOUNG | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | dendritic cell | batch_residualized | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | dendritic cell | raw | MID_vs_YOUNG | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | dendritic cell | raw | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | double negative thymocyte | batch_residualized | MID_vs_YOUNG | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | double negative thymocyte | batch_residualized | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | double negative thymocyte | raw | MID_vs_YOUNG | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | double negative thymocyte | raw | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | effector memory CD4-positive, alpha-beta T cell | batch_residualized | MID_vs_YOUNG | 50 | 17 | 0.20326 | 0.21128 | -0.0080198 | 0.082459 | 0.1302 | 0.0021004 | 0.012973 | -0.0025729 | EVALUATED |
| yazar | effector memory CD4-positive, alpha-beta T cell | batch_residualized | OLD_vs_MID | 201 | 50 | 0.20317 | 0.20326 | -8.6361e-05 | 0.96902 | 0.96902 | 0.00099955 | 0.0068053 | 0.00044066 | EVALUATED |
| yazar | effector memory CD4-positive, alpha-beta T cell | raw | MID_vs_YOUNG | 50 | 17 | 0.25655 | 0.27638 | -0.019829 | 0.00049975 | 0.001363 | 0.0038437 | 0.012806 | -0.01203 | EVALUATED |
| yazar | effector memory CD4-positive, alpha-beta T cell | raw | OLD_vs_MID | 201 | 50 | 0.26062 | 0.25655 | 0.004068 | 0.12994 | 0.17718 | 0.0022549 | 0.0074667 | 0.0049896 | EVALUATED |
| yazar | effector memory CD8-positive, alpha-beta T cell | batch_residualized | MID_vs_YOUNG | 136 | 56 | 0.18372 | 0.19114 | -0.0074215 | 0.035482 | 0.076033 | 0.0012171 | 0.0097491 | -0.0045907 | EVALUATED |
| yazar | effector memory CD8-positive, alpha-beta T cell | batch_residualized | OLD_vs_MID | 552 | 136 | 0.1714 | 0.18372 | -0.012319 | 0.00049975 | 0.0021418 | 0.00095123 | 0.0061098 | -0.0091043 | EVALUATED |
| yazar | effector memory CD8-positive, alpha-beta T cell | raw | MID_vs_YOUNG | 136 | 56 | 0.21586 | 0.23143 | -0.015568 | 0.00049975 | 0.001363 | 0.0018519 | 0.011827 | -0.011543 | EVALUATED |
| yazar | effector memory CD8-positive, alpha-beta T cell | raw | OLD_vs_MID | 552 | 136 | 0.20342 | 0.21586 | -0.012439 | 0.00049975 | 0.001363 | 0.0013821 | 0.0076769 | -0.0097371 | EVALUATED |
| yazar | erythrocyte | batch_residualized | MID_vs_YOUNG | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | erythrocyte | batch_residualized | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | erythrocyte | raw | MID_vs_YOUNG | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | erythrocyte | raw | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | gamma-delta T cell | batch_residualized | MID_vs_YOUNG | 25 | 15 | 0.19563 | 0.19862 | -0.002991 | 0.6012 | 0.7515 | 0.0018039 | 0.015719 | -0.0014406 | EVALUATED |
| yazar | gamma-delta T cell | batch_residualized | OLD_vs_MID | 88 | 25 | 0.18225 | 0.19563 | -0.013382 | 0.0054973 | 0.014993 | 0.0012535 | 0.012575 | -0.010072 | EVALUATED |
| yazar | gamma-delta T cell | raw | MID_vs_YOUNG | 25 | 15 | 0.27049 | 0.26714 | 0.0033471 | 0.6052 | 0.67244 | 0.0032082 | 0.017721 | 0.0059398 | EVALUATED |
| yazar | gamma-delta T cell | raw | OLD_vs_MID | 88 | 25 | 0.25311 | 0.27049 | -0.017378 | 0.00049975 | 0.001363 | 0.0021939 | 0.012347 | -0.011176 | EVALUATED |
| yazar | hematopoietic precursor cell | batch_residualized | MID_vs_YOUNG | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | hematopoietic precursor cell | batch_residualized | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | hematopoietic precursor cell | raw | MID_vs_YOUNG | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | hematopoietic precursor cell | raw | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | innate lymphoid cell | batch_residualized | MID_vs_YOUNG | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | innate lymphoid cell | batch_residualized | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | innate lymphoid cell | raw | MID_vs_YOUNG | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | innate lymphoid cell | raw | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | memory B cell | batch_residualized | MID_vs_YOUNG | 52 | 25 | 0.19945 | 0.19082 | 0.0086338 | 0.062969 | 0.10495 | 0.0028177 | 0.01286 | 0.0067354 | EVALUATED |
| yazar | memory B cell | batch_residualized | OLD_vs_MID | 176 | 52 | 0.19972 | 0.19945 | 0.00026695 | 0.92604 | 0.95797 | 0.0012062 | 0.0076899 | 0.0017174 | EVALUATED |
| yazar | memory B cell | raw | MID_vs_YOUNG | 52 | 25 | 0.2605 | 0.24398 | 0.016525 | 0.0014993 | 0.0034598 | 0.0043413 | 0.013465 | 0.014754 | EVALUATED |
| yazar | memory B cell | raw | OLD_vs_MID | 176 | 52 | 0.2558 | 0.2605 | -0.0047046 | 0.11194 | 0.15992 | 0.0018212 | 0.0081747 | -0.0022957 | EVALUATED |
| yazar | mucosal-associated invariant T cell | batch_residualized | MID_vs_YOUNG | 12 | 10 | 0.16576 | 0.13739 | 0.028364 | 0.11694 | 0.16706 | 0.0010663 | 0.048816 | 0.020792 | EVALUATED |
| yazar | mucosal-associated invariant T cell | batch_residualized | OLD_vs_MID | 23 | 12 | 0.14217 | 0.16576 | -0.023585 | 0.090455 | 0.13568 | 0.0016144 | 0.039071 | -0.017585 | EVALUATED |
| yazar | mucosal-associated invariant T cell | raw | MID_vs_YOUNG | 12 | 10 | 0.26961 | 0.27228 | -0.0026651 | 0.67516 | 0.72339 | 0.0010754 | 0.01699 | -0.00079642 | EVALUATED |
| yazar | mucosal-associated invariant T cell | raw | OLD_vs_MID | 23 | 12 | 0.26379 | 0.26961 | -0.0058203 | 0.24938 | 0.29754 | 0.0029481 | 0.014177 | 8.5734e-05 | EVALUATED |
| yazar | naive B cell | batch_residualized | MID_vs_YOUNG | 50 | 21 | 0.2001 | 0.19435 | 0.0057509 | 0.13043 | 0.17787 | 0.0014054 | 0.010679 | 0.0059012 | EVALUATED |
| yazar | naive B cell | batch_residualized | OLD_vs_MID | 251 | 50 | 0.19956 | 0.2001 | -0.00053978 | 0.84308 | 0.93675 | 0.00044902 | 0.0076121 | -0.00045767 | EVALUATED |
| yazar | naive B cell | raw | MID_vs_YOUNG | 50 | 21 | 0.23854 | 0.23818 | 0.0003627 | 0.94353 | 0.97606 | 0.0019061 | 0.015469 | 0.0029831 | EVALUATED |
| yazar | naive B cell | raw | OLD_vs_MID | 251 | 50 | 0.24455 | 0.23854 | 0.0060134 | 0.10545 | 0.15817 | 0.00058209 | 0.010324 | 0.0053411 | EVALUATED |
| yazar | naive thymus-derived CD4-positive, alpha-beta T cell | batch_residualized | MID_vs_YOUNG | 163 | 89 | 0.16704 | 0.17329 | -0.006253 | 0.056472 | 0.099656 | 0.00091205 | 0.0093302 | -0.0042269 | EVALUATED |
| yazar | naive thymus-derived CD4-positive, alpha-beta T cell | batch_residualized | OLD_vs_MID | 564 | 163 | 0.15486 | 0.16704 | -0.012182 | 0.00049975 | 0.0021418 | 0.00049253 | 0.0067691 | -0.0091168 | EVALUATED |
| yazar | naive thymus-derived CD4-positive, alpha-beta T cell | raw | MID_vs_YOUNG | 163 | 89 | 0.19729 | 0.20714 | -0.0098514 | 0.0089955 | 0.016867 | 0.001226 | 0.010679 | -0.0070485 | EVALUATED |
| yazar | naive thymus-derived CD4-positive, alpha-beta T cell | raw | OLD_vs_MID | 564 | 163 | 0.18588 | 0.19729 | -0.011412 | 0.00049975 | 0.001363 | 0.00074117 | 0.0079 | -0.0090323 | EVALUATED |
| yazar | naive thymus-derived CD8-positive, alpha-beta T cell | batch_residualized | MID_vs_YOUNG | 102 | 77 | 0.18276 | 0.16748 | 0.015277 | 0.00049975 | 0.0021418 | 0.00069541 | 0.0091947 | 0.012214 | EVALUATED |
| yazar | naive thymus-derived CD8-positive, alpha-beta T cell | batch_residualized | OLD_vs_MID | 210 | 102 | 0.18946 | 0.18276 | 0.0067025 | 0.0029985 | 0.0089955 | 0.0003376 | 0.0061021 | 0.0061667 | EVALUATED |
| yazar | naive thymus-derived CD8-positive, alpha-beta T cell | raw | MID_vs_YOUNG | 102 | 77 | 0.2293 | 0.2164 | 0.012897 | 0.0014993 | 0.0034598 | 0.001786 | 0.01069 | 0.010663 | EVALUATED |
| yazar | naive thymus-derived CD8-positive, alpha-beta T cell | raw | OLD_vs_MID | 210 | 102 | 0.23618 | 0.2293 | 0.0068801 | 0.0089955 | 0.016867 | 0.00057516 | 0.0075701 | 0.0065099 | EVALUATED |
| yazar | natural killer cell | batch_residualized | MID_vs_YOUNG | 148 | 66 | 0.18212 | 0.1939 | -0.011779 | 0.00049975 | 0.0021418 | 0.0014377 | 0.0093255 | -0.0089053 | EVALUATED |
| yazar | natural killer cell | batch_residualized | OLD_vs_MID | 598 | 148 | 0.16526 | 0.18212 | -0.016857 | 0.00049975 | 0.0021418 | 0.00032907 | 0.0059684 | -0.012043 | EVALUATED |
| yazar | natural killer cell | raw | MID_vs_YOUNG | 148 | 66 | 0.21495 | 0.22629 | -0.011342 | 0.003998 | 0.0085671 | 0.0018825 | 0.011027 | -0.0088938 | EVALUATED |
| yazar | natural killer cell | raw | OLD_vs_MID | 598 | 148 | 0.19435 | 0.21495 | -0.020596 | 0.00049975 | 0.001363 | 0.00058316 | 0.0070215 | -0.01541 | EVALUATED |
| yazar | plasmablast | batch_residualized | MID_vs_YOUNG | 5 | 1 | 0.052775 | NA | NA | NA | NA | NA | NA | 0.042935 | UNTESTABLE_BIN_FLOOR |
| yazar | plasmablast | batch_residualized | OLD_vs_MID | 8 | 5 | 0.10698 | 0.052775 | 0.054201 | NA | NA | 0.020478 | NA | 0.021058 | UNTESTABLE_BIN_FLOOR |
| yazar | plasmablast | raw | MID_vs_YOUNG | 5 | 1 | 0.26659 | NA | NA | NA | NA | NA | NA | 0.17595 | UNTESTABLE_BIN_FLOOR |
| yazar | plasmablast | raw | OLD_vs_MID | 8 | 5 | 0.27261 | 0.26659 | 0.0060155 | NA | NA | 0.024977 | NA | 0.010295 | UNTESTABLE_BIN_FLOOR |
| yazar | plasmacytoid dendritic cell | batch_residualized | MID_vs_YOUNG | 1 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | plasmacytoid dendritic cell | batch_residualized | OLD_vs_MID | 2 | 1 | 1.0537e-08 | NA | NA | NA | NA | NA | NA | 0 | UNTESTABLE_BIN_FLOOR |
| yazar | plasmacytoid dendritic cell | raw | MID_vs_YOUNG | 1 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | plasmacytoid dendritic cell | raw | OLD_vs_MID | 2 | 1 | 0.25458 | NA | NA | NA | NA | NA | NA | 0.13387 | UNTESTABLE_BIN_FLOOR |
| yazar | platelet | batch_residualized | MID_vs_YOUNG | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | platelet | batch_residualized | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | platelet | raw | MID_vs_YOUNG | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | platelet | raw | OLD_vs_MID | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | regulatory T cell | batch_residualized | MID_vs_YOUNG | 10 | 0 | 0.19965 | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | regulatory T cell | batch_residualized | OLD_vs_MID | 59 | 10 | 0.16894 | 0.19965 | -0.030718 | 0.011994 | 0.027678 | 0.0023384 | 0.031629 | -0.022551 | EVALUATED |
| yazar | regulatory T cell | raw | MID_vs_YOUNG | 10 | 0 | 0.28456 | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_BIN_FLOOR |
| yazar | regulatory T cell | raw | OLD_vs_MID | 59 | 10 | 0.27515 | 0.28456 | -0.009413 | 0.029985 | 0.052915 | 0.0051824 | 0.01205 | 0.0024373 | EVALUATED |
| yazar | transitional stage B cell | batch_residualized | MID_vs_YOUNG | 39 | 26 | 0.19551 | 0.18016 | 0.015353 | 0.009995 | 0.024988 | 0.0018257 | 0.016783 | 0.013358 | EVALUATED |
| yazar | transitional stage B cell | batch_residualized | OLD_vs_MID | 147 | 39 | 0.19687 | 0.19551 | 0.0013558 | 0.67316 | 0.78134 | 0.0022126 | 0.0089001 | 0.0030346 | EVALUATED |
| yazar | transitional stage B cell | raw | MID_vs_YOUNG | 39 | 26 | 0.25873 | 0.26412 | -0.0053855 | 0.1909 | 0.23863 | 0.0018563 | 0.011685 | -0.0037219 | EVALUATED |
| yazar | transitional stage B cell | raw | OLD_vs_MID | 147 | 39 | 0.25884 | 0.25873 | 0.00011347 | 0.98051 | 0.98051 | 0.0024061 | 0.0097255 | 0.0010453 | EVALUATED |

### Verdicts

- Gao: INCONCLUSIVE. PROVES 0, DOES NOT PROVE 0, INCONCLUSIVE 10, UNTESTABLE 3 cell types.
  Inconclusive units: Cardiomyocyte, Endocardial, Endothelial, Fibroblast, LEC, Lymphoid, Myeloid, Neuronal, Pericyte, vSMC.
  Untestable units: Adipocyte, Epicardial, Mast.
- Yazar: PROVES. PROVES 2, DOES NOT PROVE 0, INCONCLUSIVE 12, UNTESTABLE 14 cell types.
  Proven units: memory B cell, naive thymus-derived CD8-positive, alpha-beta T cell.
  Inconclusive units: CD14-low, CD16-positive monocyte, CD14-positive monocyte, central memory CD4-positive, alpha-beta T cell, central memory CD8-positive, alpha-beta T cell, effector memory CD4-positive, alpha-beta T cell, effector memory CD8-positive, alpha-beta T cell, gamma-delta T cell, mucosal-associated invariant T cell, naive B cell, naive thymus-derived CD4-positive, alpha-beta T cell, natural killer cell, transitional stage B cell.
  Untestable units: CD16-negative, CD56-bright natural killer cell, human, CD4-positive, alpha-beta T cell, CD4-positive, alpha-beta cytotoxic T cell, CD8-positive, alpha-beta T cell, conventional dendritic cell, dendritic cell, double negative thymocyte, erythrocyte, hematopoietic precursor cell, innate lymphoid cell, plasmablast, plasmacytoid dendritic cell, platelet, regulatory T cell.

### Distinction from the Nov 17 distance-variance test

The Nov 17 test reduced each donor to one distance to a young reference and then compared the variance of those donor distances. Block 3 instead measures dispersion within each age bin directly from pairwise distances among donor expression centroids. It does not use a young reference and does not analyze the variance of one distance per donor. These are different estimands.
- Gao: this profile-cloud test does not clearly agree or disagree with the Nov 17 batch-adjusted null because this test is inconclusive. The Nov 17 headline was Gao OLD to YOUNG distance-variance ratio 0.79, p 0.85; Yazar raw ratio 1.62, p 1.9e-07, changing to pool-residualized ratio 1.25, p 0.43.
- Yazar: this profile-cloud test does not agree with the Nov 17 batch-adjusted null, although the estimands differ. The Nov 17 headline was Gao OLD to YOUNG distance-variance ratio 0.79, p 0.85; Yazar raw ratio 1.62, p 1.9e-07, changing to pool-residualized ratio 1.25, p 0.43.

### Reproducibility and family sizes

Permutation base seed 20261118; split-half seeds 20261116, 20261117, 20261118, 20261119, 20261120. Exact comparison seeds are in the CSV. All 2000 permutations resample donor labels, never cells.
Evaluated BH family sizes: gao batch_residualized mid_secondary m=23; gao batch_residualized old_young m=10; gao raw mid_secondary m=23; gao raw old_young m=10; yazar batch_residualized mid_secondary m=30; yazar batch_residualized old_young m=14; yazar raw mid_secondary m=30; yazar raw old_young m=14.
No primary cells were downsampled. Full-depth within-cell proportions use only cells passing the 1500-count eligibility filter.
