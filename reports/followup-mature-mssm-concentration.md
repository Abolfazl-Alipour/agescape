# Exploratory MSSM mature-to-old concentration test

This owner-approved exploratory run tests controls aged 30 years and older in MSSM DLPFC. Run model: kimi-code/k3.1-pro. It does not enter the manuscript automatically; Reza decides whether and how to use it after seeing the results. The manuscript TeX was not edited.

## Input, probe, and gene audit

Input: `/home/reza/projects/agescape-revision-data/mssm.h5ad`; 36,092,670,231 bytes; shape 4,140,453 cells x 34,176 genes. All cells are primary DLPFC data, assay ["10x 3' v3"], suspension ['nucleus']. The run used `obs['cell_type']`, not `subclass`.

Controls: 1,049,809 cells from 261 donors. Exactly 214 donors have a parseable age, range 19 to 89; 47 lack one. At age 30+, there are 209 donors: 129 male and 80 female. These exactly match the supplied probe counts. The 30-39 reference contains 8 donors.

Gene audit: 34,176 raw names; 167 duplicate occurrences collapsed; 33 MT, 659 RP/RPL/RPS, 95 MRP, 0 ERCC, and 5 sex genes excluded; 33,217 retained. Sex genes present: XIST, TSIX, DDX3Y, UTY, KDM5D; absent: EIF2S3Y.

## Validation before MSSM

Gao Lymphoid was independently read from raw count HDF5 slices with the adapted gene projection and per-cell functions. Full-range and mature 30+ donor scatter, Shannon, depth, and cell counts matched the committed donor rows, and the fitted slopes, standard errors, p-values, donor counts, and cell counts matched their respective committed model rows within rounding. Both validations PASS.

The requested `results/followup-fulldepth-gao-spread.csv` contains Shannon/KL rather than Hellinger scatter. Accordingly, scatter was checked against `results/followup-fulldepth-gao-cell-drift.csv`, and Shannon was checked against the requested spread CSV. The mature equivalents were also checked. No Yazar data were sought or used.

## Methods

Counts were converted to within-cell proportions after the exact original gene audit. Cells below 1,500 post-audit molecules were excluded as an eligibility filter only; no primary outcome was downsampled. Donor is the inferential unit. A donor-cell-type group needs at least 20 eligible cells, and a slope needs at least six donors.

Scatter is the mean per-cell Hellinger distance from the leave-one-out mean proportions of the donor's other cells of that type. Shannon is -sum(p log p). E1 is the original mean detected off-program gene count: young target detection <=5%, another-type detection >=20%, other-type CP10K >=20 and >=5 times target CP10K+1, with original reference and leave-one-out gates. E2 uses genes detected in <1% of young-reference cells and is evaluated in the original bounded strata [1500,2500), [2500,4500), [4500,9000), [9000,18000), with >=20 cells per donor-cell-type-stratum and the original >=200 leave-one-out young-reference cell gate. Counts are detected genes per cell, not UMIs.

Per-type OLS is outcome ~ age_per_decade + sex_male + mean(log UMI), exactly the original primary depth covariate, with the original minimum residual df=5. Pooled puts all of a donor's cells in one group; pooled E1/E2 combine the cell-type-specific gene-set scores. Stratified-global standardizes within cell type (or cell type x E2 stratum), centers sex and depth within group, includes cell-type/stratum dummies, and clusters SEs by donor. The original implementation uses linear mean-log-UMI in this regression, not nonlinear spline terms; the ambiguity was resolved in favor of the original code. Split-half uses the five original seeds [20261116, 20261117, 20261118, 20261119, 20261120]; scale-invariant outcomes use the original 0.5 sampling-SE factor and counts use 1.0. BH is across tested cell types within each measure (and within each E2 stratum).

## Eligibility and QC

The depth-gated slope rule yields 21 eligible cell types. The pre-depth donor/cell probe also yields 21 types, exactly matching the supplied value 21. Across all types, 823,904 of 825,793 cells pass the depth gate and 1,889 are excluded. The original residual-df>=5 model gate leaves 19 fitted scatter/Shannon types: T cells and smooth muscle cells each have seven donors and are not fitted. Donor-age handling exactly reproduces 209 mature donors, so it does not explain the Nov 23 count of 14.

| Cell type | Donors >=20 pre-depth | Donors >=20 post-depth | Cells >=1500 | Slope eligible | E1 markers |
| --- | --- | --- | --- | --- | --- |
| B cell | 2 | 2 | 93 | no | 0 |
| GABAergic neuron | 153 | 153 | 13219 | yes | 2 |
| L2/3 intratelencephalic projecting glutamatergic neuron | 195 | 195 | 88163 | yes | 1 |
| L2/3-6 intratelencephalic projecting glutamatergic neuron | 193 | 193 | 58536 | yes | 1 |
| L5 extratelencephalic projecting glutamatergic cortical neuron | 0 | 0 | 184 | no | 0 |
| L5/6 near-projecting glutamatergic neuron | 87 | 87 | 4175 | yes | 0 |
| L6 corticothalamic-projecting glutamatergic cortical neuron | 125 | 125 | 7526 | yes | 3 |
| L6 intratelencephalic projecting glutamatergic neuron | 123 | 123 | 10299 | yes | 1 |
| L6b glutamatergic neuron of the primary motor cortex | 74 | 74 | 4319 | yes | 0 |
| T cell | 7 | 7 | 778 | yes | 0 |
| VIP GABAergic cortical interneuron | 193 | 193 | 39232 | yes | 2 |
| astrocyte | 203 | 203 | 104411 | yes | 5 |
| endothelial cell | 89 | 87 | 4543 | yes | 2 |
| lamp5 GABAergic cortical interneuron | 175 | 175 | 14850 | yes | 2 |
| microglial cell | 196 | 196 | 36837 | yes | 8 |
| natural killer cell | 2 | 2 | 528 | no | 0 |
| oligodendrocyte | 206 | 206 | 309065 | yes | 9 |
| oligodendrocyte precursor cell | 205 | 205 | 55383 | yes | 5 |
| pericyte | 84 | 82 | 3910 | yes | 3 |
| perivascular macrophage | 24 | 24 | 1775 | yes | 0 |
| plasma cell | 0 | 0 | 15 | no | 0 |
| pvalb GABAergic cortical interneuron | 185 | 185 | 30761 | yes | 0 |
| smooth muscle cell | 7 | 7 | 779 | yes | 0 |
| sst GABAergic cortical interneuron | 174 | 174 | 30742 | yes | 4 |
| vascular leptomeningeal cell | 41 | 41 | 3781 | yes | 2 |

14-vs-21 reconciliation: 21 is the general `cell_type`-level six-donor/20-cell count both before and after depth gating. Endpoint-specific reference gates are narrower: 15 types have usable nonzero E1 markers, and the numbers of tested E2 types are s1_1500_2500: 2, s2_2500_4500: 4, s3_4500_9000: 6, s4_9000_18000: 8. The exact Nov 23 feasibility definition is not present in the supplied instruction, so its 14 cannot honestly be reproduced by relabeling the general count. The nearest endpoint-specific count is 15 for E1, indicating an additional endpoint/reference/type filter (or a different earlier type level), not a donor-age discrepancy. No count was forced to 14.

## Per-cell-type slopes

Entries are raw-outcome slopes per decade with p and BH q. Stratified-global effects are SD units, shown separately below.

### scatter

| Cell type | Slope | SE | p | BH q | Donors | Cells | MDE80 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GABAergic neuron | 0.00027984 | 0.00058263 | 0.6317 | 0.8573 | 153 | 12692 | 0.001643 |
| L2/3 intratelencephalic projecting glutamatergic neuron | -1.2244e-05 | 0.00065342 | 0.9851 | 0.9851 | 195 | 88093 | 0.00184 |
| L2/3-6 intratelencephalic projecting glutamatergic neuron | -0.0004223 | 0.00060396 | 0.4853 | 0.7926 | 193 | 58441 | 0.0017008 |
| L5/6 near-projecting glutamatergic neuron | -0.00064862 | 0.0008761 | 0.4612 | 0.7926 | 87 | 3394 | 0.0024837 |
| L6 corticothalamic-projecting glutamatergic cortical neuron | 0.0011176 | 0.00066567 | 0.09575 | 0.4548 | 125 | 6880 | 0.0018801 |
| L6 intratelencephalic projecting glutamatergic neuron | -0.0010507 | 0.0010463 | 0.3173 | 0.7569 | 123 | 9623 | 0.0029554 |
| L6b glutamatergic neuron of the primary motor cortex | -0.00012763 | 0.0010877 | 0.9069 | 0.9676 | 74 | 3373 | 0.0030903 |
| VIP GABAergic cortical interneuron | -5.3481e-05 | 0.00051032 | 0.9166 | 0.9676 | 193 | 39126 | 0.0014371 |
| astrocyte | 0.00072126 | 0.00057802 | 0.2136 | 0.7031 | 203 | 104382 | 0.0016273 |
| endothelial cell | 0.00053889 | 0.00086424 | 0.5346 | 0.7926 | 87 | 3456 | 0.0024501 |
| lamp5 GABAergic cortical interneuron | -0.00045395 | 0.00062894 | 0.4714 | 0.7926 | 175 | 14587 | 0.0017721 |
| microglial cell | 0.001707 | 0.00041754 | 6.386e-05 | 0.001213 | 196 | 36721 | 0.0011757 |
| oligodendrocyte | 0.0010884 | 0.00042466 | 0.0111 | 0.1055 | 206 | 309041 | 0.0011955 |
| oligodendrocyte precursor cell | 0.00013071 | 0.00057761 | 0.8212 | 0.9676 | 205 | 55369 | 0.0016261 |
| pericyte | 0.0015175 | 0.0006941 | 0.0318 | 0.2014 | 82 | 2648 | 0.0019692 |
| perivascular macrophage | -0.00030332 | 0.0016062 | 0.8521 | 0.9676 | 24 | 849 | 0.0047316 |
| pvalb GABAergic cortical interneuron | -0.00072364 | 0.00059055 | 0.222 | 0.7031 | 185 | 30605 | 0.0016634 |
| sst GABAergic cortical interneuron | -0.00029645 | 0.00048556 | 0.5423 | 0.7926 | 174 | 30475 | 0.0013682 |
| vascular leptomeningeal cell | 0.001033 | 0.001022 | 0.3187 | 0.7569 | 41 | 2633 | 0.0029409 |

### shannon

| Cell type | Slope | SE | p | BH q | Donors | Cells | MDE80 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GABAergic neuron | -0.0045071 | 0.0055764 | 0.4202 | 0.6654 | 153 | 12692 | 0.015726 |
| L2/3 intratelencephalic projecting glutamatergic neuron | -0.0087676 | 0.0052599 | 0.09718 | 0.3693 | 195 | 88093 | 0.014812 |
| L2/3-6 intratelencephalic projecting glutamatergic neuron | -0.0047374 | 0.0053792 | 0.3796 | 0.6654 | 193 | 58441 | 0.015149 |
| L5/6 near-projecting glutamatergic neuron | -0.003425 | 0.008652 | 0.6932 | 0.8745 | 87 | 3394 | 0.024528 |
| L6 corticothalamic-projecting glutamatergic cortical neuron | -0.0025033 | 0.006205 | 0.6873 | 0.8745 | 125 | 6880 | 0.017525 |
| L6 intratelencephalic projecting glutamatergic neuron | -0.0081061 | 0.0070732 | 0.2541 | 0.6559 | 123 | 9623 | 0.01998 |
| L6b glutamatergic neuron of the primary motor cortex | -0.0022585 | 0.0081489 | 0.7825 | 0.8745 | 74 | 3373 | 0.023153 |
| VIP GABAergic cortical interneuron | -0.012024 | 0.0053601 | 0.02604 | 0.3693 | 193 | 39126 | 0.015095 |
| astrocyte | 8.3243e-05 | 0.0040723 | 0.9837 | 0.9928 | 203 | 104382 | 0.011465 |
| endothelial cell | -5.7926e-05 | 0.006378 | 0.9928 | 0.9928 | 87 | 3456 | 0.018081 |
| lamp5 GABAergic cortical interneuron | -0.011099 | 0.0061459 | 0.0727 | 0.3693 | 175 | 14587 | 0.017317 |
| microglial cell | 0.0080233 | 0.0042748 | 0.06205 | 0.3693 | 196 | 36721 | 0.012037 |
| oligodendrocyte | 0.0008925 | 0.0030044 | 0.7667 | 0.8745 | 206 | 309041 | 0.0084579 |
| oligodendrocyte precursor cell | 0.0030237 | 0.0035629 | 0.3971 | 0.6654 | 205 | 55369 | 0.010031 |
| pericyte | 0.016529 | 0.0109 | 0.1334 | 0.4226 | 82 | 2648 | 0.030924 |
| perivascular macrophage | -0.021597 | 0.01929 | 0.2762 | 0.6559 | 24 | 849 | 0.056828 |
| pvalb GABAergic cortical interneuron | -0.010604 | 0.0061667 | 0.08721 | 0.3693 | 185 | 30605 | 0.01737 |
| sst GABAergic cortical interneuron | -0.0051076 | 0.0053068 | 0.3372 | 0.6654 | 174 | 30475 | 0.014953 |
| vascular leptomeningeal cell | 0.0046469 | 0.012179 | 0.705 | 0.8745 | 41 | 2633 | 0.035047 |

### E1

| Cell type | Slope | SE | p | BH q | Donors | Cells | MDE80 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GABAergic neuron | -0.0013037 | 0.0028229 | 0.6449 | 0.7909 | 153 | 12692 | 0.0079608 |
| L2/3 intratelencephalic projecting glutamatergic neuron | -0.0021177 | 0.0014018 | 0.1326 | 0.4971 | 191 | 86984 | 0.003948 |
| L2/3-6 intratelencephalic projecting glutamatergic neuron | -0.0025344 | 0.0012148 | 0.03832 | 0.3934 | 190 | 57297 | 0.0034212 |
| L6 corticothalamic-projecting glutamatergic cortical neuron | -0.00055001 | 0.0044312 | 0.9014 | 0.9014 | 125 | 6880 | 0.012515 |
| L6 intratelencephalic projecting glutamatergic neuron | -0.0055251 | 0.0028202 | 0.05245 | 0.3934 | 122 | 9421 | 0.0079668 |
| VIP GABAergic cortical interneuron | -0.0012554 | 0.0023938 | 0.6006 | 0.7909 | 193 | 39126 | 0.0067413 |
| astrocyte | -0.00092574 | 0.0027651 | 0.7381 | 0.7909 | 203 | 104382 | 0.0077849 |
| endothelial cell | -0.0033762 | 0.0041563 | 0.4189 | 0.7292 | 87 | 3456 | 0.011783 |
| lamp5 GABAergic cortical interneuron | 0.0007633 | 0.0019012 | 0.6886 | 0.7909 | 175 | 14587 | 0.0053568 |
| microglial cell | -0.0020449 | 0.0044206 | 0.6442 | 0.7909 | 196 | 36721 | 0.012448 |
| oligodendrocyte | -0.0049784 | 0.005505 | 0.3669 | 0.7292 | 206 | 309041 | 0.015497 |
| oligodendrocyte precursor cell | 0.002152 | 0.0023623 | 0.3634 | 0.7292 | 205 | 55369 | 0.0066503 |
| pericyte | 0.0096027 | 0.0054205 | 0.08037 | 0.4019 | 82 | 2648 | 0.015379 |
| sst GABAergic cortical interneuron | 0.0020837 | 0.0026775 | 0.4375 | 0.7292 | 174 | 30475 | 0.0075445 |
| vascular leptomeningeal cell | 0.0064925 | 0.006492 | 0.3238 | 0.7292 | 41 | 2633 | 0.018681 |

### E2 depth-stratum counts

| Cell type | Stratum | Slope | SE | p | BH q | Donors | Cells | MDE80 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| microglial cell | s1_1500_2500 | 1.2637 | 0.48591 | 0.01055 | 0.0211 | 117 | 6787 | 1.3732 |
| oligodendrocyte | s1_1500_2500 | 0.24351 | 0.46094 | 0.5983 | 0.5983 | 119 | 18549 | 1.3024 |
| astrocyte | s2_2500_4500 | 0.22459 | 0.47149 | 0.6345 | 0.6345 | 162 | 16428 | 1.3291 |
| microglial cell | s2_2500_4500 | 1.181 | 0.45466 | 0.01023 | 0.0409 | 171 | 15698 | 1.2813 |
| oligodendrocyte | s2_2500_4500 | 0.57381 | 0.4177 | 0.1711 | 0.2282 | 197 | 112034 | 1.1762 |
| oligodendrocyte precursor cell | s2_2500_4500 | 1.0655 | 0.70758 | 0.1346 | 0.2282 | 129 | 8789 | 1.998 |
| L2/3 intratelencephalic projecting glutamatergic neuron | s3_4500_9000 | -1.0858 | 0.76667 | 0.1594 | 0.239 | 121 | 8997 | 2.166 |
| astrocyte | s3_4500_9000 | 0.8205 | 0.46551 | 0.07954 | 0.239 | 199 | 40297 | 1.3107 |
| microglial cell | s3_4500_9000 | 0.97858 | 0.57348 | 0.09031 | 0.239 | 135 | 11354 | 1.6187 |
| oligodendrocyte | s3_4500_9000 | 0.44045 | 0.46702 | 0.3468 | 0.4162 | 199 | 138545 | 1.315 |
| oligodendrocyte precursor cell | s3_4500_9000 | 0.18709 | 0.61559 | 0.7615 | 0.7615 | 194 | 25376 | 1.7335 |
| sst GABAergic cortical interneuron | s3_4500_9000 | -1.0023 | 0.65707 | 0.1297 | 0.239 | 126 | 7858 | 1.8557 |
| GABAergic neuron | s4_9000_18000 | -1.3552 | 0.75631 | 0.07631 | 0.4187 | 100 | 4247 | 2.1406 |
| L2/3 intratelencephalic projecting glutamatergic neuron | s4_9000_18000 | -0.7051 | 0.496 | 0.157 | 0.4187 | 170 | 24142 | 1.3978 |
| VIP GABAergic cortical interneuron | s4_9000_18000 | -0.50738 | 0.57724 | 0.3807 | 0.4786 | 166 | 14854 | 1.627 |
| astrocyte | s4_9000_18000 | 0.45472 | 0.52753 | 0.3898 | 0.4786 | 186 | 34847 | 1.4859 |
| oligodendrocyte | s4_9000_18000 | 0.67241 | 0.5671 | 0.2377 | 0.4753 | 151 | 38058 | 1.5994 |
| oligodendrocyte precursor cell | s4_9000_18000 | 0.56647 | 0.69882 | 0.4188 | 0.4786 | 169 | 16286 | 1.9694 |
| pvalb GABAergic cortical interneuron | s4_9000_18000 | -0.99739 | 0.63951 | 0.1213 | 0.4187 | 132 | 9058 | 1.8054 |
| sst GABAergic cortical interneuron | s4_9000_18000 | -0.030571 | 0.69635 | 0.965 | 0.965 | 137 | 10656 | 1.9653 |

## Pooled and stratified-global

| Measure | Stratum | Type | Slope | SE | p | Donors | Cells | MDE80 | SH RMS | SH sampling SE | SH rho | SH SNR |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| scatter |  | pooled | 0.00041416 | 0.00065484 | 0.5278 | 209 | 823904 | 0.0018434 | 0.0001438 | 3.7577e-06 | 3.3012e-05 | 110.22 |
| scatter |  | stratified-global | -0.00085109 | 0.012446 | 0.9455 | 208 | 812388 | 0.035032 | 0.020646 | 0.0068901 | 0.0010165 | 0.12352 |
| shannon |  | pooled | -0.0024758 | 0.0033728 | 0.4638 | 209 | 823904 | 0.0094945 | 0.00099778 | 2.6205e-05 | 5.9909e-05 | 94.477 |
| shannon |  | stratified-global | -0.0095183 | 0.021774 | 0.6625 | 208 | 812388 | 0.06129 | 0.027843 | 0.0070839 | 0.00048479 | 1.3437 |
| E1 |  | pooled | 0.0013018 | 0.0031706 | 0.6818 | 209 | 778042 | 0.0089251 | 0.012021 | 0.00062199 | 0.039024 | 2.093 |
| E1 |  | stratified-global | -0.018329 | 0.026546 | 0.4907 | 208 | 771712 | 0.074724 | 0.58732 | 0.01441 | 0.42707 | 1.2719 |
| E2_stratum | s1_1500_2500 | pooled | 0.020323 | 0.41689 | 0.9612 | 155 | 26221 | 1.1755 | 0.76505 | 0.046999 | 0.012642 | 0.43241 |
| E2_stratum | s1_1500_2500 | stratified-global | 0.099993 | 0.061605 | 0.1068 | 143 | 25336 | 0.17379 | 0.12791 | 0.0061354 | 0.016719 | 16.298 |
| E2_stratum | s2_2500_4500 | pooled | 0.37487 | 0.41168 | 0.3636 | 204 | 154574 | 1.159 | 0.39361 | 0.020408 | 0.0024555 | 18.369 |
| E2_stratum | s2_2500_4500 | stratified-global | 0.084553 | 0.046776 | 0.07218 | 200 | 152949 | 0.13169 | 0.089061 | 0.0024919 | 0.0082206 | 33.932 |
| E2_stratum | s3_4500_9000 | pooled | 0.37847 | 0.47164 | 0.4232 | 207 | 234597 | 1.3277 | 0.27923 | 0.014701 | 0.00096664 | 25.745 |
| E2_stratum | s3_4500_9000 | stratified-global | 0.022069 | 0.046379 | 0.6347 | 205 | 232427 | 0.13056 | 0.073771 | 0.0016361 | 0.0056655 | 13.489 |
| E2_stratum | s4_9000_18000 | pooled | -0.37926 | 0.53304 | 0.4776 | 204 | 155583 | 1.5007 | 0.29928 | 0.015776 | 0.00089358 | 24.04 |
| E2_stratum | s4_9000_18000 | stratified-global | -0.018257 | 0.046298 | 0.6938 | 203 | 152148 | 0.13034 | 0.075201 | 0.0016173 | 0.0062177 | 11.289 |

Pooled slopes are raw outcome units; stratified-global slopes are within-cell-type SD units. Pooled q=p as a singleton; original stratified-global q is NA. SH RMS is the root-mean-square donor-level half difference, SH sampling SE uses the original scale factor, SH rho is the technical variance ratio, and SH SNR is |slope|/SH sampling SE.

## Split-half technical null and passing types

The pooled/global table reports the actual five-seed null summaries. Per-type RMS half differences, sampling SEs, technical variance ratios, and SNR are in the slopes CSV. The fixed-seed nulls quantify molecule-sampling repeatability; they are not biological permutation tests.

- BH PASS: microglial cell scatter : slope 0.001707, p 6.386e-05, q 0.001213.
- BH PASS: microglial cell E2_stratum s1_1500_2500: slope 1.2637, p 0.01055, q 0.0211.
- BH PASS: microglial cell E2_stratum s2_2500_4500: slope 1.181, p 0.01023, q 0.0409.

Nominal per-type p < 0.05: 7 rows; BH-passing: 3 rows.
- Nominal: microglial cell scatter : slope 0.001707, p 6.386e-05, q 0.001213.
- Nominal: oligodendrocyte scatter : slope 0.0010884, p 0.0111, q 0.1055.
- Nominal: pericyte scatter : slope 0.0015175, p 0.0318, q 0.2014.
- Nominal: VIP GABAergic cortical interneuron shannon : slope -0.012024, p 0.02604, q 0.3693.
- Nominal: L2/3-6 intratelencephalic projecting glutamatergic neuron E1 : slope -0.0025344, p 0.03832, q 0.3934.
- Nominal: microglial cell E2_stratum s1_1500_2500: slope 1.2637, p 0.01055, q 0.0211.
- Nominal: microglial cell E2_stratum s2_2500_4500: slope 1.181, p 0.01023, q 0.0409.

Nominally significant pooled/global rows (p < 0.05): none.

All three BH passes are microglial-cell results with positive slopes: more within-cell Hellinger scatter and more detected E2 off-program genes with age. This is the opposite direction from expression concentration. Their split-half SNRs are 76.61 (scatter), 18.27 (E2 s1), and 24.86 (E2 s2), so the point estimates are large relative to molecular sampling noise. However, no pooled or stratified-global row passes p < 0.05, and the E2 global summaries for s1 and s2 are p 0.1068 and 0.0722. These are localized exploratory microglial associations, not a general MSSM concentration replication.


## Power and interpretation

The pre-stated power context is MDE80 = 0.029 per decade. The test is powered for Gao-sized effects (-0.0723 per decade) and larger, not Yazar-sized effects (-0.0197 per decade, approximately 453 donors needed). A null here does not rule out a Yazar-sized fall. Row-specific model-based MDE80 values are also reported; they do not replace that pre-stated context.

## Limitations

- Exploratory, single-cohort, post-mortem, cross-sectional DLPFC analysis; donor age is not longitudinal change.
- Only normal controls aged 30+ enter. Eight 30-39 reference donors make E1/E2 marker estimation a potential source of instability; types without reference sets are untestable for those outcomes, not null. Usable E1 sets contain only 1 to 9 markers.
- No study/pool fixed effects exist in this single cohort. Available metadata have constant assay and suspension and no batch, postmortem interval, RIN, or processing-cohort adjustment; residual technical confounding remains possible.
- Depth adjustment is donor-level mean log UMI, exactly as in the original primary. Cell composition can affect pooled results; stratified-global is the within-type summary.
- E1/E2 detected-gene counts are not molecule counts. E2 strata have different donor sets and BH families and are not independent replications; the s1 and s2 BH families contain only 2 and 4 cell types.
- Split-half only estimates molecular sampling noise and cannot rule out ambient RNA, annotation error, or biological donor confounding.
- No result enters the paper automatically.
