# Block 8: Nonlinear / categorical age-group signature models

Status: specification frozen before running any Block 8 score or model.

## Goal

Test whether the signature effect is linear in age or concentrated in one age transition. All gene-level and signature tests so far use a linear age slope, which assumes linearity. The drift result itself is a nonlinearity finding (effect concentrated in young-to-mature, dead against the mature reference). This block asks the same question of the signature: is the effect a maturation contrast (MATURE to MID), an aging contrast (MID to OLD), or a steady trend?

## Inputs

- Gao and Yazar h5ad files (same paths as Block 1), same donor tables, same technical gene mask.
- P1m gene set (1,914 mouse genes, 1,119 up / 795 down; 1,724 human-mapped, 1,013 up / 711 down) from `results/followup-mature-p1m-genes.csv` (commit d3dc14c); verify counts before running. Do NOT re-derive it and never score the original young-derived P1 in this block.
- Scoring code: reuse the composite-score computation from the transport pipeline (`scripts/run_m9_pipeline.py` / follow-up transport scripts). Per-donor per-cell-type composite scores must be computed from full-depth within-cell proportions, the same way the transport models score donors.
- Age bins: donors under 30 are excluded from all models in this block. MATURE is 30 to 39 and replaces YOUNG (20 to 39) as the reference bin; MID stays 40 to 59 and OLD stays 60 and older.

## Exact analysis specification

- Per donor per cell type (eligibility floors per shared rule 6): P1m composite score = mean up-gene score minus mean down-gene score, donor pseudobulk, full depth, technical mask applied.
- Primary model, per dataset and cell type: OLS of composite score on categorical age bin (MATURE reference) + sex_male + batch (Gao study indicators; Yazar pool fixed effects). Report MID vs MATURE and OLD vs MATURE coefficients with 95 percent CIs and p-values, plus the OLD vs MID contrast (mature-to-old change, computed from the same fit, not a separate model).
- Trend test secondary: ordinal model with bin coded 1, 2, 3 (MATURE, MID, OLD), same covariates; report slope, CI, p-value. A significant categorical pattern with a null trend test means the effect is a step, not a slope; report that distinction explicitly.
- Unadjusted sensitivity: categorical model without batch covariates.
- Spline sensitivity (labeled, secondary): cubic B-spline on continuous age with 3 degrees of freedom, F test vs the linear model, mirroring the Phase 2 nonlinearity analysis in `reports/reviewer_response_phase2.md` section 2 and `scripts/run_reviewer_phase2_human_age.py`. Read that script and reuse its spline construction.
- Family (state before correction): BH q within each dataset across cell types, separately for the MID vs MATURE contrast, the OLD vs MID contrast, and the trend test. Three stated families per dataset; do not merge them.

## Fixed implementation details

- The P1m count was verified before model fitting from `results/followup-mature-p1m-genes.csv`: 1,914 unique mouse genes, 1,119 UP and 795 DOWN, with 1,724 human-mapped genes (1,013 UP and 711 DOWN) from `scripts/run_mature_p1m_transport.py` `load_programs()`. Mapped genes absent from a dataset are omitted from the corresponding directional denominator; detected mapped genes are not filtered by expression level.
- The existing Block 1 donor-sharded cache is reused. It contains raw integer counts after duplicate-symbol collapse and the human technical mask, restricted to post-mask UMI >= 1500 cells. The large h5ad matrices are not reread. The human technical gene mask excludes uppercase `MT-`, `RP[L/S]`, `MRP[L/S]`, `XIST`, `TSIX`, `DDX3Y`, `EIF2S3Y`, `UTY`, and `KDM5D`. ERCC controls are also excluded by the cache construction.
- For each cell and direction, the gene score is the mean within-cell proportion over the mapped genes present in that dataset. For each donor and cell type, up and down scores are arithmetic means of these cell scores across eligible cells, and composite = up minus down. This is an equal-cell-weight donor pseudobulk of full-depth within-cell proportions, with no log transform and no downsampling. Cells below 1500 post-mask UMIs were excluded as an eligibility filter only.
- A donor-cell-type pair requires at least 20 eligible cells. Categorical group models require at least 6 eligible donors in each of MATURE, MID, and OLD, both before and after batch pruning. Donors under 30 are excluded before any floor is checked. At least 6 donors total are required for every fit. The repository covariate pruning requires at least 2 donors per retained batch, full rank, and residual df floors of 5 for Gao and 10 for Yazar in adjusted models and 5 in unadjusted models. Units below a floor are UNTESTABLE, never nulls, and no unadjusted fallback replaces a failed primary fit.
- Age bins are mature bins derived from `age_years`: UNDER30 under 30, MATURE 30 to 39, MID 40 to 59, OLD 60 and older. UNDER30 scores may be present in the donor-score CSV but are excluded from every model. The Block 1 cache bins are still asserted unchanged as a cache-integrity check. Gao uses intercept plus 6 indicators for its 7 studies with ENCODE v4 (Snyder) as reference when present. Yazar uses intercept plus 74 indicators for its 75 pools. Both include sex_male.
- The OLD vs MID contrast is beta_OLD minus beta_MID from the categorical fit, with its standard error from that fit's covariance matrix. The ordinal trend estimate is per one-bin increment. OLD vs MATURE is reported descriptively and is not one of the three correction families. BH q-values are assigned only to evaluable primary rows in the three stated families. Unadjusted and spline sensitivities and OLD vs MATURE have no q-value because the specification defines exactly three families per dataset.
- The spline uses `dmatrix("bs(age, df=3, degree=3, include_intercept=False)", part)` exactly as in `scripts/run_reviewer_phase2_human_age.py`, on the same pruned donors and with the same sex and batch covariates as the primary. The nested linear model uses continuous age in years plus those covariates. The F test compares spline to linear and has 2 numerator degrees of freedom when both designs are full rank.
- For categorical contrasts and trend, MDE80 is the exact two-sided noncentral-t 80 percent power effect, computed as the fitted-contrast standard error times the required noncentrality. It is reported in composite-score units and divided by the donor-score SD in the fitted cell type. A small MDE80 is predefined as at most 0.50 donor-score SD per contrast or ordinal-bin step, matching the repository convention. For the spline F test, MDE80 is the exact noncentral-F noncentrality giving 80 percent power at alpha 0.05; its partial-R2 equivalent is also reported.
- No split-half technical null is required for these donor-score regressions. No random procedure is used, so no random seed is needed. No new downloads are used. Available memory is checked from `/proc/meminfo` MemAvailable before each cache-shard stage; below 20 Gi the process waits and retries every 60 seconds.

## Outputs

- `results/followup-mature-block8-categorical-signature.csv`: one row per dataset, cell type, contrast or test, model (primary, unadjusted, spline), with estimate, CI, p, q, N donors per bin, MDE80, and an age_window column recording 30-plus.
- `results/followup-mature-block8-donor-scores.csv`: per-donor per-cell-type P1m composite scores with age bin, sex, and batch, so the categorical fits are reproducible from the CSV alone. UNDER30 rows are present in this CSV but excluded from every model.
- `reports/followup-mature-block8-categorical.md`: spec, per-cell-type classification (maturation-concentrated, aging-concentrated, steady trend, null), and verdicts.
- `scripts/run_mature_block8_categorical_signature.py`.

## Pre-stated classification and verdict rules

A cell type is classified maturation-concentrated only if MID vs MATURE passes BH correction and OLD vs MID does not; aging-concentrated only if OLD vs MID passes correction and MID vs MATURE does not; steady trend only if the trend test passes correction and both steps agree in direction. Transition-concentrated classifications take precedence if more than one rule is met. Both steps agree in direction only when both estimated steps are nonzero and have the same sign. All other evaluable cell types are null. Unevaluable cell types are UNTESTABLE.

Block-level verdict: PROVES (signature effect is nonlinear / transition-concentrated) only if at least one cell type is classified maturation-concentrated or aging-concentrated after correction in the batch-adjusted primary model. DOES NOT PROVE if all passing cell types are steady-trend or if nothing passes with small MDE80. INCONCLUSIVE if nothing passes with large MDE80.

For a no-transition dataset or overall summary, any steady-trend cell type gives DOES NOT PROVE. If no cell type passes any of the three primary families, DOES NOT PROVE requires every evaluated primary MDE80 in those families to be at most 0.50 donor-score SD; otherwise the summary is INCONCLUSIVE. A significant pattern that meets no classification rule is INCONCLUSIVE. The overall summary is PROVES if either dataset is PROVES, DOES NOT PROVE only if both datasets are DOES NOT PROVE, and otherwise INCONCLUSIVE.

The classification will be compared against mature linear-slope results from `results/followup-mature-block5-gao-yazar-transport.csv` (primary arm, gene set P1m, composite BH q). For each passing cell type, the report states whether the mature linear model had already flagged it at composite BH q < 0.05. These mature transport fits use their log1p(CPM) depth-adjusted implementation while the categorical fits use full-depth proportion scores, so the comparison concerns prior flagging, not numerical identity with the new full-depth proportion score.

Exactly one deep review round follows completion of the outputs, followed by fixes, verification, commit, and push.

## Results

All estimates are on the full-depth within-cell-proportion composite scale. MID vs MATURE and OLD vs MATURE are coefficients with MATURE (30 to 39) as reference. OLD vs MID is OLD minus MID from the same categorical fit. Trend is per ordinal-bin increment. NA denotes UNTESTABLE or not applicable. BH q-values exist only for the three pre-stated primary families. Every fitted donor is at least 30 years old; donors under 30 are excluded from every model.

### Gao data and family sizes

The cache supplied 1,296 scored donor-cell-type rows, including 87 UNDER30 rows excluded from models. Scoring used 977 of 1013 mapped UP genes and 685 of 711 mapped DOWN genes. Absent UP genes: ABCB6, CASP8AP2, CCND2, CD177, CD24, CDR2, CEBPE, CIRSR, COXFA4L3, CTSE, CYBB, DDX3Y, EIF4G1, FAM53B, FEN1, FXN, GP9, H1-0, HEMK2, LCN2, LIME1, LUC7L2, LY75, MAGEB16, NME1, PDZK1IP1, PELO, PRDX2, PRP4K, PYCR2, RBM27, SH2D1A, SOWAHA, TRIM10, TSPO2, ZBED6. Absent DOWN genes: BSCL2, C3, CAMP, CCNO, CRLF2, CTSD, CXCR1, EIF4EBP3, GPR141, HIGD1A, KDELR2, LIMD2, LRG1, MT1B, PPP3R1, RAB4B, RBFA, SLC35H1, SLC9D1, SPN, SPNS2, TAP2, TGFB1, TM9SF1, TMEM120A, TOMM6. The source cache has 522,761 eligible full-depth cells after the technical mask and 1500-count eligibility filter.

| Primary family | Evaluated cell types | BH q < 0.05 |
| --- | --- | --- |
| MID_vs_MATURE | 9 | 0 |
| OLD_vs_MID | 9 | 0 |
| trend | 9 | 0 |

### Yazar data and family sizes

The cache supplied 20,858 scored donor-cell-type rows, including 1,070 UNDER30 rows excluded from models. Scoring used 1009 of 1013 mapped UP genes and 709 of 711 mapped DOWN genes. Absent UP genes: CIRSR, COXFA4L3, DDX3Y, HEMK2. Absent DOWN genes: SLC35H1, SLC9D1. The source cache has 617,760 eligible full-depth cells after the technical mask and 1500-count eligibility filter.

| Primary family | Evaluated cell types | BH q < 0.05 |
| --- | --- | --- |
| MID_vs_MATURE | 12 | 0 |
| OLD_vs_MID | 12 | 0 |
| trend | 12 | 1 |

### Primary categorical and trend results

#### Gao

| Cell type | Status | N MAT/MID/OLD | MID vs MATURE estimate (95 percent CI); p; q | OLD vs MATURE estimate (95 percent CI); p | OLD vs MID estimate (95 percent CI); p; q | Trend estimate (95 percent CI); p; q | MDE80 SD MID / OLD-MID / trend |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Adipocyte | UNTESTABLE_DONOR_OR_BIN_FLOOR | 3/20/12 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| Cardiomyocyte | EVALUATED | 11/51/39 | 6.560e-07 (-2.304e-06 to 3.616e-06); p=0.6608; q=0.96099 | 1.587e-06 (-1.455e-06 to 4.628e-06); p=0.30284 | 9.308e-07 (-9.831e-07 to 2.845e-06); p=0.33662; q=0.53481 | 8.383e-07 (-5.313e-07 to 2.208e-06); p=0.22726; q=0.40908 | 0.76799 / 0.49663 / 0.35541 |
| Endocardial | EVALUATED | 8/29/29 | 2.409e-07 (-5.016e-06 to 5.498e-06); p=0.92721; q=0.96099 | 1.560e-06 (-3.697e-06 to 6.816e-06); p=0.55471 | 1.319e-06 (-2.130e-06 to 4.768e-06); p=0.44701; q=0.57472 | 9.449e-07 (-1.424e-06 to 3.313e-06); p=0.42784; q=0.55008 | 1.1368 / 0.74593 / 0.51226 |
| Endothelial | EVALUATED | 10/39/30 | -1.612e-07 (-4.250e-06 to 3.927e-06); p=0.93754; q=0.96099 | 2.154e-06 (-2.094e-06 to 6.402e-06); p=0.3155 | 2.315e-06 (-5.238e-07 to 5.154e-06); p=0.10837; q=0.40505 | 1.417e-06 (-5.584e-07 to 3.393e-06); p=0.15705; q=0.35336 | 0.88279 / 0.61296 / 0.42658 |
| Epicardial | UNTESTABLE_DONOR_OR_BIN_FLOOR | 2/8/8 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| Fibroblast | EVALUATED | 11/49/38 | 1.265e-07 (-4.996e-06 to 5.249e-06); p=0.96099; q=0.96099 | 3.965e-06 (-1.298e-06 to 9.227e-06); p=0.13794 | 3.838e-06 (4.750e-07 to 7.202e-06); p=0.025777; q=0.23199 | 2.567e-06 (1.676e-07 to 4.966e-06); p=0.03629; q=0.26595 | 0.89587 / 0.58821 / 0.4196 |
| LEC | UNTESTABLE_DONOR_OR_BIN_FLOOR | 4/17/13 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| Lymphoid | EVALUATED | 9/25/22 | 2.870e-06 (-5.542e-06 to 1.128e-05); p=0.49602; q=0.96099 | 3.025e-06 (-5.497e-06 to 1.155e-05); p=0.47886 | 1.548e-07 (-6.128e-06 to 6.438e-06); p=0.96069; q=0.96069 | 1.225e-06 (-2.801e-06 to 5.252e-06); p=0.54373; q=0.61169 | 0.93042 / 0.69493 / 0.44541 |
| Mast | UNTESTABLE_DONOR_OR_BIN_FLOOR | 1/8/8 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| Myeloid | EVALUATED | 11/49/39 | 3.391e-06 (-1.475e-06 to 8.257e-06); p=0.16961; q=0.96099 | 4.939e-06 (-5.727e-08 to 9.936e-06); p=0.052627 | 1.548e-06 (-1.608e-06 to 4.705e-06); p=0.33243; q=0.53481 | 2.171e-06 (-8.516e-08 to 4.427e-06); p=0.0591; q=0.26595 | 0.89287 / 0.57919 / 0.41402 |
| Neuronal | EVALUATED | 9/28/19 | 2.430e-06 (-3.752e-06 to 8.612e-06); p=0.43323; q=0.96099 | 4.555e-06 (-2.019e-06 to 1.113e-05); p=0.16997 | 2.125e-06 (-2.465e-06 to 6.715e-06); p=0.35654; q=0.53481 | 2.243e-06 (-8.584e-07 to 5.344e-06); p=0.1525; q=0.35336 | 1.072 / 0.79589 / 0.53782 |
| Pericyte | EVALUATED | 11/44/33 | -1.165e-06 (-6.635e-06 to 4.304e-06); p=0.67267; q=0.96099 | 1.698e-06 (-3.947e-06 to 7.344e-06); p=0.55111 | 2.864e-06 (-9.107e-07 to 6.638e-06); p=0.13502; q=0.40505 | 1.408e-06 (-1.217e-06 to 4.033e-06); p=0.28891; q=0.43337 | 0.81524 / 0.56259 / 0.39126 |
| vSMC | EVALUATED | 9/40/30 | 2.433e-07 (-5.277e-06 to 5.763e-06); p=0.93021; q=0.96099 | 5.498e-07 (-5.120e-06 to 6.220e-06); p=0.84723 | 3.065e-07 (-3.322e-06 to 3.935e-06); p=0.86672; q=0.96069 | 2.848e-07 (-2.281e-06 to 2.851e-06); p=0.82548; q=0.82548 | 0.92122 / 0.60557 / 0.42821 |

#### Yazar

| Cell type | Status | N MAT/MID/OLD | MID vs MATURE estimate (95 percent CI); p; q | OLD vs MATURE estimate (95 percent CI); p | OLD vs MID estimate (95 percent CI); p; q | Trend estimate (95 percent CI); p; q | MDE80 SD MID / OLD-MID / trend |
| --- | --- | --- | --- | --- | --- | --- | --- |
| CD14-low, CD16-positive monocyte | EVALUATED | 8/41/196 | 6.546e-07 (-4.063e-06 to 5.372e-06); p=0.78462; q=0.85595 | 9.298e-07 (-3.521e-06 to 5.381e-06); p=0.68076 | 2.752e-07 (-1.950e-06 to 2.501e-06); p=0.80755; q=0.88096 | 3.701e-07 (-1.296e-06 to 2.036e-06); p=0.66177; q=0.80746 | 0.77516 / 0.36571 / 0.27375 |
| CD14-positive monocyte | EVALUATED | 20/87/295 | -4.624e-06 (-9.221e-06 to -2.711e-08); p=0.048672; q=0.29203 | -5.179e-06 (-9.510e-06 to -8.482e-07); p=0.019234 | -5.553e-07 (-2.921e-06 to 1.811e-06); p=0.64462; q=0.86754 | -1.679e-06 (-3.399e-06 to 4.103e-08); p=0.055688; q=0.22275 | 0.61748 / 0.31784 / 0.23106 |
| CD16-negative, CD56-bright natural killer cell, human | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/1/4 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| CD4-positive, alpha-beta T cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| CD4-positive, alpha-beta cytotoxic T cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 1/13/78 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| CD8-positive, alpha-beta T cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| central memory CD4-positive, alpha-beta T cell | EVALUATED | 58/191/660 | -1.075e-06 (-2.506e-06 to 3.572e-07); p=0.14108; q=0.42323 | -7.015e-07 (-2.029e-06 to 6.262e-07); p=0.30001 | 3.731e-07 (-4.483e-07 to 1.195e-06); p=0.37283; q=0.80222 | -8.070e-08 (-6.400e-07 to 4.787e-07); p=0.77712; q=0.80746 | 0.19717 / 0.11311 / 0.077025 |
| central memory CD8-positive, alpha-beta T cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 4/22/80 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| conventional dendritic cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/5 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| dendritic cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| double negative thymocyte | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| effector memory CD4-positive, alpha-beta T cell | EVALUATED | 6/49/195 | 7.989e-06 (1.751e-06 to 1.423e-05); p=0.012341; q=0.1481 | 1.122e-05 (5.295e-06 to 1.714e-05); p=0.00024536 | 3.227e-06 (8.261e-07 to 5.629e-06); p=0.0086887; q=0.10426 | 4.155e-06 (2.226e-06 to 6.083e-06); p=3.303e-05; q=0.00039635 | 0.77488 / 0.29826 / 0.23952 |
| effector memory CD8-positive, alpha-beta T cell | EVALUATED | 33/136/552 | -8.260e-07 (-2.830e-06 to 1.178e-06); p=0.41853; q=0.77574 | -1.194e-06 (-3.065e-06 to 6.778e-07); p=0.21087 | -3.675e-07 (-1.397e-06 to 6.619e-07); p=0.48355; q=0.82893 | -4.948e-07 (-1.234e-06 to 2.442e-07); p=0.18905; q=0.3781 | 0.29878 / 0.15348 / 0.11019 |
| erythrocyte | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| gamma-delta T cell | EVALUATED | 7/23/81 | 2.329e-06 (-5.130e-06 to 9.789e-06); p=0.53574; q=0.77574 | -1.142e-06 (-8.471e-06 to 6.187e-06); p=0.75708 | -3.471e-06 (-7.477e-06 to 5.344e-07); p=0.088392; q=0.35357 | -1.852e-06 (-4.871e-06 to 1.167e-06); p=0.22549; q=0.38656 | 1.1212 / 0.60205 / 0.45373 |
| hematopoietic precursor cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| innate lymphoid cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| memory B cell | EVALUATED | 9/48/170 | -1.141e-06 (-6.039e-06 to 3.758e-06); p=0.64645; q=0.77574 | 1.181e-06 (-3.385e-06 to 5.746e-06); p=0.61046 | 2.321e-06 (4.492e-08 to 4.598e-06); p=0.045694; q=0.27416 | 1.459e-06 (-2.313e-07 to 3.148e-06); p=0.090267; q=0.2708 | 0.71156 / 0.33064 / 0.24546 |
| mucosal-associated invariant T cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 3/12/23 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| naive B cell | EVALUATED | 11/49/247 | -1.090e-06 (-5.167e-06 to 2.986e-06); p=0.59895; q=0.77574 | -1.385e-06 (-5.249e-06 to 2.479e-06); p=0.48089 | -2.949e-07 (-2.300e-06 to 1.710e-06); p=0.77237; q=0.88096 | -5.024e-07 (-1.992e-06 to 9.870e-07); p=0.50715; q=0.76072 | 0.65874 / 0.324 / 0.24067 |
| naive thymus-derived CD4-positive, alpha-beta T cell | EVALUATED | 50/163/564 | -1.124e-06 (-2.506e-06 to 2.568e-07); p=0.11043; q=0.42323 | -1.470e-06 (-2.755e-06 to -1.851e-07); p=0.025 | -3.453e-07 (-1.152e-06 to 4.617e-07); p=0.40111; q=0.80222 | -5.937e-07 (-1.140e-06 to -4.727e-08); p=0.033257; q=0.19954 | 0.16476 / 0.096266 / 0.065185 |
| naive thymus-derived CD8-positive, alpha-beta T cell | EVALUATED | 39/99/206 | 1.169e-07 (-1.826e-06 to 2.060e-06); p=0.90586; q=0.90586 | -8.950e-07 (-2.790e-06 to 9.999e-07); p=0.35332 | -1.012e-06 (-2.343e-06 to 3.195e-07); p=0.13578; q=0.40733 | -5.981e-07 (-1.472e-06 to 2.754e-07); p=0.17879; q=0.3781 | 0.22856 / 0.15659 / 0.10274 |
| natural killer cell | EVALUATED | 38/147/598 | 1.051e-06 (-7.854e-07 to 2.888e-06); p=0.26148; q=0.62756 | 8.279e-07 (-8.650e-07 to 2.521e-06); p=0.33732 | -2.234e-07 (-1.192e-06 to 7.447e-07); p=0.65065; q=0.86754 | 1.442e-07 (-5.332e-07 to 8.215e-07); p=0.67616; q=0.80746 | 0.28659 / 0.15106 / 0.10569 |
| plasmablast | UNTESTABLE_DONOR_OR_BIN_FLOOR | 1/5/8 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| plasmacytoid dendritic cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/1/2 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| platelet | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| regulatory T cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/10/59 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| transitional stage B cell | EVALUATED | 9/35/141 | 1.550e-06 (-4.477e-06 to 7.576e-06); p=0.61196; q=0.77574 | 1.338e-06 (-4.499e-06 to 7.176e-06); p=0.65103 | -2.114e-07 (-3.507e-06 to 3.085e-06); p=0.89926; q=0.89926 | 2.973e-07 (-2.110e-06 to 2.705e-06); p=0.80746; q=0.80746 | 0.88114 / 0.48189 / 0.35203 |

### Unadjusted categorical sensitivity

#### Gao

| Cell type | Status | N MAT/MID/OLD | MID vs MATURE | OLD vs MATURE | OLD vs MID |
| --- | --- | --- | --- | --- | --- |
| Adipocyte | UNTESTABLE_DONOR_OR_BIN_FLOOR | 3/20/12 | NA; p=NA | NA; p=NA | NA; p=NA |
| Cardiomyocyte | EVALUATED | 11/52/39 | 2.527e-06; p=0.23172 | 1.528e-06; p=0.48214 | -9.987e-07; p=0.45597 |
| Endocardial | EVALUATED | 8/29/29 | 2.266e-06; p=0.39935 | 2.760e-06; p=0.30609 | 4.941e-07; p=0.77898 |
| Endothelial | EVALUATED | 10/40/30 | 1.180e-06; p=0.62292 | 1.854e-06; p=0.45683 | 6.738e-07; p=0.67961 |
| Epicardial | UNTESTABLE_DONOR_OR_BIN_FLOOR | 2/8/8 | NA; p=NA | NA; p=NA | NA; p=NA |
| Fibroblast | EVALUATED | 12/49/38 | -5.094e-07; p=0.85351 | 1.313e-06; p=0.64414 | 1.822e-06; p=0.32254 |
| LEC | UNTESTABLE_DONOR_OR_BIN_FLOOR | 4/17/13 | NA; p=NA | NA; p=NA | NA; p=NA |
| Lymphoid | EVALUATED | 9/25/22 | 3.609e-06; p=0.47115 | 3.075e-06; p=0.54478 | -5.347e-07; p=0.88501 |
| Mast | UNTESTABLE_DONOR_OR_BIN_FLOOR | 1/8/8 | NA; p=NA | NA; p=NA | NA; p=NA |
| Myeloid | EVALUATED | 11/49/39 | 3.837e-06; p=0.1443 | 4.430e-06; p=0.1003 | 5.933e-07; p=0.72229 |
| Neuronal | EVALUATED | 9/28/19 | 2.931e-06; p=0.34939 | 5.171e-06; p=0.12009 | 2.239e-06; p=0.34192 |
| Pericyte | EVALUATED | 11/44/33 | 2.789e-07; p=0.93198 | 1.471e-06; p=0.66177 | 1.192e-06; p=0.59089 |
| vSMC | EVALUATED | 9/40/30 | 8.536e-07; p=0.7898 | 3.419e-07; p=0.91712 | -5.117e-07; p=0.80613 |

#### Yazar

| Cell type | Status | N MAT/MID/OLD | MID vs MATURE | OLD vs MATURE | OLD vs MID |
| --- | --- | --- | --- | --- | --- |
| CD14-low, CD16-positive monocyte | EVALUATED | 9/44/204 | -7.753e-07; p=0.81013 | -9.923e-07; p=0.74121 | -2.170e-07; p=0.88249 |
| CD14-positive monocyte | EVALUATED | 21/88/299 | -4.505e-06; p=0.082956 | -5.856e-06; p=0.015486 | -1.350e-06; p=0.29736 |
| CD16-negative, CD56-bright natural killer cell, human | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/1/4 | NA; p=NA | NA; p=NA | NA; p=NA |
| CD4-positive, alpha-beta T cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA; p=NA | NA; p=NA | NA; p=NA |
| CD4-positive, alpha-beta cytotoxic T cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 1/13/78 | NA; p=NA | NA; p=NA | NA; p=NA |
| CD8-positive, alpha-beta T cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA; p=NA | NA; p=NA | NA; p=NA |
| central memory CD4-positive, alpha-beta T cell | EVALUATED | 58/191/660 | -9.412e-07; p=0.54293 | -3.620e-06; p=0.010496 | -2.679e-06; p=0.0016131 |
| central memory CD8-positive, alpha-beta T cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 4/22/80 | NA; p=NA | NA; p=NA | NA; p=NA |
| conventional dendritic cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/5 | NA; p=NA | NA; p=NA | NA; p=NA |
| dendritic cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA; p=NA | NA; p=NA | NA; p=NA |
| double negative thymocyte | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA; p=NA | NA; p=NA | NA; p=NA |
| effector memory CD4-positive, alpha-beta T cell | EVALUATED | 7/50/201 | 5.007e-06; p=0.28715 | 4.896e-06; p=0.27466 | -1.106e-07; p=0.95203 |
| effector memory CD8-positive, alpha-beta T cell | EVALUATED | 34/136/552 | 6.472e-07; p=0.72453 | -1.539e-06; p=0.36324 | -2.186e-06; p=0.017139 |
| erythrocyte | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA; p=NA | NA; p=NA | NA; p=NA |
| gamma-delta T cell | EVALUATED | 7/25/88 | 8.274e-07; p=0.83873 | 8.914e-07; p=0.81343 | 6.405e-08; p=0.97615 |
| hematopoietic precursor cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA; p=NA | NA; p=NA | NA; p=NA |
| innate lymphoid cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA; p=NA | NA; p=NA | NA; p=NA |
| memory B cell | EVALUATED | 10/52/176 | -1.190e-06; p=0.73275 | -8.277e-07; p=0.80085 | 3.624e-07; p=0.82002 |
| mucosal-associated invariant T cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 3/12/23 | NA; p=NA | NA; p=NA | NA; p=NA |
| naive B cell | EVALUATED | 11/50/251 | 1.767e-07; p=0.95308 | 1.644e-06; p=0.55421 | 1.468e-06; p=0.29326 |
| naive thymus-derived CD4-positive, alpha-beta T cell | EVALUATED | 51/163/564 | -2.312e-06; p=0.21955 | -7.407e-06; p=1.767e-05 | -5.095e-06; p=1.225e-06 |
| naive thymus-derived CD8-positive, alpha-beta T cell | EVALUATED | 40/102/210 | -8.139e-07; p=0.70548 | -8.505e-06; p=2.617e-05 | -7.691e-06; p=6.317e-08 |
| natural killer cell | EVALUATED | 38/148/598 | 8.758e-07; p=0.59982 | 1.111e-07; p=0.94232 | -7.647e-07; p=0.36362 |
| plasmablast | UNTESTABLE_DONOR_OR_BIN_FLOOR | 1/5/8 | NA; p=NA | NA; p=NA | NA; p=NA |
| plasmacytoid dendritic cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/1/2 | NA; p=NA | NA; p=NA | NA; p=NA |
| platelet | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA; p=NA | NA; p=NA | NA; p=NA |
| regulatory T cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/10/59 | NA; p=NA | NA; p=NA | NA; p=NA |
| transitional stage B cell | EVALUATED | 9/39/147 | -1.118e-06; p=0.75828 | -1.962e-06; p=0.56069 | -8.443e-07; p=0.63368 |

### Cubic B-spline sensitivity

| Dataset | Cell type | Status | F | p | df num | df den | MDE80 noncentrality | MDE80 partial R2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gao | Adipocyte | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| gao | Cardiomyocyte | EVALUATED | 0.62714 | 0.53641 | 2 | 91 | 9.959 | 0.098644 |
| gao | Endocardial | EVALUATED | 0.12693 | 0.88105 | 2 | 56 | 10.169 | 0.15369 |
| gao | Endothelial | EVALUATED | 1.5063 | 0.22883 | 2 | 70 | 10.059 | 0.12565 |
| gao | Epicardial | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| gao | Fibroblast | EVALUATED | 1.4591 | 0.23801 | 2 | 88 | 9.9704 | 0.10177 |
| gao | LEC | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| gao | Lymphoid | EVALUATED | 2.5732 | 0.087004 | 2 | 47 | 10.276 | 0.17941 |
| gao | Mast | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| gao | Myeloid | EVALUATED | 3.0398 | 0.052838 | 2 | 89 | 9.9665 | 0.10071 |
| gao | Neuronal | EVALUATED | 0.90909 | 0.40986 | 2 | 47 | 10.276 | 0.17941 |
| gao | Pericyte | EVALUATED | 0.86251 | 0.42604 | 2 | 79 | 10.01 | 0.11246 |
| gao | vSMC | EVALUATED | 0.44182 | 0.64465 | 2 | 70 | 10.059 | 0.12565 |
| yazar | CD14-low, CD16-positive monocyte | EVALUATED | 0.1847 | 0.8315 | 2 | 190 | 9.7882 | 0.048993 |
| yazar | CD14-positive monocyte | EVALUATED | 0.82035 | 0.44115 | 2 | 340 | 9.7201 | 0.027794 |
| yazar | CD16-negative, CD56-bright natural killer cell, human | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | CD4-positive, alpha-beta T cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | CD4-positive, alpha-beta cytotoxic T cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | CD8-positive, alpha-beta T cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | central memory CD4-positive, alpha-beta T cell | EVALUATED | 2.9974 | 0.050456 | 2 | 830 | 9.6695 | 0.011516 |
| yazar | central memory CD8-positive, alpha-beta T cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | conventional dendritic cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | dendritic cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | double negative thymocyte | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | effector memory CD4-positive, alpha-beta T cell | EVALUATED | 1.8767 | 0.15582 | 2 | 197 | 9.7827 | 0.047309 |
| yazar | effector memory CD8-positive, alpha-beta T cell | EVALUATED | 0.13318 | 0.87533 | 2 | 643 | 9.6797 | 0.014831 |
| yazar | erythrocyte | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | gamma-delta T cell | EVALUATED | 0.37309 | 0.68991 | 2 | 73 | 10.041 | 0.12092 |
| yazar | hematopoietic precursor cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | innate lymphoid cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | memory B cell | EVALUATED | 1.9618 | 0.14367 | 2 | 176 | 9.8006 | 0.052748 |
| yazar | mucosal-associated invariant T cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | naive B cell | EVALUATED | 0.87299 | 0.41891 | 2 | 263 | 9.7453 | 0.03573 |
| yazar | naive thymus-derived CD4-positive, alpha-beta T cell | EVALUATED | 2.0611 | 0.12809 | 2 | 699 | 9.6761 | 0.013654 |
| yazar | naive thymus-derived CD8-positive, alpha-beta T cell | EVALUATED | 0.14352 | 0.86637 | 2 | 282 | 9.7378 | 0.033379 |
| yazar | natural killer cell | EVALUATED | 2.4752 | 0.084875 | 2 | 705 | 9.6757 | 0.013539 |
| yazar | plasmablast | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | plasmacytoid dendritic cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | platelet | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | regulatory T cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | transitional stage B cell | EVALUATED | 2.2433 | 0.10999 | 2 | 137 | 9.8485 | 0.067066 |

### Classification and comparison with prior linear models

Transition classifications require the stated BH pattern in the batch-adjusted primary model. A significant categorical step with a null trend is a step pattern, not a steady slope. Linear flagged means the mature Block 5 linear model had composite BH q < 0.05 in its primary P1m arm.

#### Gao

| Cell type | Classification | Prior linear slope per decade | Prior linear q | Linear flagged |
| --- | --- | --- | --- | --- |
| Adipocyte | UNTESTABLE | 0.032853 | 0.25774 | no |
| Cardiomyocyte | null | 0.012263 | 0.25774 | no |
| Endocardial | null | -0.009473 | 0.45264 | no |
| Endothelial | null | 0.008606 | 0.3807 | no |
| Epicardial | UNTESTABLE | 0.0091301 | 0.61273 | no |
| Fibroblast | null | 0.012932 | 0.30955 | no |
| LEC | UNTESTABLE | 0.032052 | 0.30955 | no |
| Lymphoid | null | 0.014019 | 0.3807 | no |
| Mast | UNTESTABLE | 0.018537 | 0.3807 | no |
| Myeloid | null | 0.006274 | 0.45658 | no |
| Neuronal | null | 0.034161 | 0.21686 | no |
| Pericyte | null | 0.010728 | 0.3807 | no |
| vSMC | null | 0.013074 | 0.3807 | no |

#### Yazar

| Cell type | Classification | Prior linear slope per decade | Prior linear q | Linear flagged |
| --- | --- | --- | --- | --- |
| CD14-low, CD16-positive monocyte | null | -0.015428 | 0.014185 | yes |
| CD14-positive monocyte | null | -0.00016521 | 0.98512 | no |
| CD16-negative, CD56-bright natural killer cell, human | UNTESTABLE | NA | NA | not tested |
| CD4-positive, alpha-beta T cell | UNTESTABLE | NA | NA | not tested |
| CD4-positive, alpha-beta cytotoxic T cell | UNTESTABLE | -0.00015283 | 0.98512 | no |
| CD8-positive, alpha-beta T cell | UNTESTABLE | NA | NA | not tested |
| central memory CD4-positive, alpha-beta T cell | null | -0.00045647 | 0.98512 | no |
| central memory CD8-positive, alpha-beta T cell | UNTESTABLE | 0.00020668 | 0.98512 | no |
| conventional dendritic cell | UNTESTABLE | NA | NA | not tested |
| dendritic cell | UNTESTABLE | NA | NA | not tested |
| double negative thymocyte | UNTESTABLE | NA | NA | not tested |
| effector memory CD4-positive, alpha-beta T cell | steady trend | -0.00047321 | 0.98512 | no |
| effector memory CD8-positive, alpha-beta T cell | null | -0.0040325 | 0.18243 | no |
| erythrocyte | UNTESTABLE | NA | NA | not tested |
| gamma-delta T cell | null | -0.0015029 | 0.98512 | no |
| hematopoietic precursor cell | UNTESTABLE | NA | NA | not tested |
| innate lymphoid cell | UNTESTABLE | NA | NA | not tested |
| memory B cell | null | 0.0022362 | 0.98512 | no |
| mucosal-associated invariant T cell | UNTESTABLE | -0.012212 | 0.83441 | no |
| naive B cell | null | 0.0054655 | 0.27942 | no |
| naive thymus-derived CD4-positive, alpha-beta T cell | null | -0.00033234 | 0.98512 | no |
| naive thymus-derived CD8-positive, alpha-beta T cell | null | -0.0030132 | 0.83441 | no |
| natural killer cell | null | -0.0015554 | 0.83441 | no |
| plasmablast | UNTESTABLE | NA | NA | not tested |
| plasmacytoid dendritic cell | UNTESTABLE | NA | NA | not tested |
| platelet | UNTESTABLE | NA | NA | not tested |
| regulatory T cell | UNTESTABLE | -0.0057336 | 0.70148 | no |
| transitional stage B cell | null | -0.0052608 | 0.70148 | no |

- Gao verdict: INCONCLUSIVE. Classifications: 0 maturation-concentrated, 0 aging-concentrated, 0 steady trend, 9 null, 4 UNTESTABLE.
- Yazar verdict: DOES NOT PROVE. Classifications: 0 maturation-concentrated, 0 aging-concentrated, 1 steady trend, 11 null, 16 UNTESTABLE.

Block-level verdict: INCONCLUSIVE. This verdict applies to the P1m score in Gao heart and Yazar blood under the pre-stated categorical and trend rules.

The unadjusted and spline models are sensitivities and do not determine the primary classification. The spline F test addresses departure from continuous linearity but is not substituted for the three categorical BH families.

### Survived or died versus the original Block 8

The original Block 8 verdict was PROVES (Gao 1 maturation-concentrated, Yazar 2 maturation-concentrated and 1 steady trend). Under the mature P1m score with the MATURE 30 to 39 baseline, no maturation-concentrated or aging-concentrated label survives in either dataset. All three original transition labels died, and one new steady-trend label appeared.

| Cell type | Original label | Original key q | Mature label | Mature key q | Outcome |
| --- | --- | --- | --- | --- | --- |
| Gao Myeloid | maturation-concentrated | MID vs YOUNG q=0.028581 | null | MID vs MATURE q=0.96099; OLD vs MID q=0.53481; trend q=0.26595 | DIED |
| Yazar central memory CD4-positive, alpha-beta T cell | maturation-concentrated | MID vs YOUNG q=0.0057793 | null | MID vs MATURE q=0.42323; OLD vs MID q=0.80222; trend q=0.80746 | DIED |
| Yazar naive thymus-derived CD4-positive, alpha-beta T cell | maturation-concentrated | MID vs YOUNG q=0.026064 | null | MID vs MATURE q=0.42323; OLD vs MID q=0.80222; trend q=0.19954 | DIED |
| Yazar naive thymus-derived CD8-positive, alpha-beta T cell | steady trend | trend q=0.0087412 | null | MID vs MATURE q=0.90586; OLD vs MID q=0.40733; trend q=0.3781 | DIED |
| Yazar effector memory CD4-positive, alpha-beta T cell | null | MID p=0.094132; trend p=0.062154 | steady trend | trend q=0.00039635; MID vs MATURE q=0.1481; OLD vs MID q=0.10426 | NEW |

Notes on the four deaths: the Gao Myeloid mature MID step is 3.391e-06 with p=0.16961, against an original MID step of 2.704e-05 with p=0.0028581; the P1m composite scale differs from P1, so the comparison is on the corrected verdict, not on the raw estimate size. Both Yazar naive and central memory CD4 mature MID steps are negative (minus 1.124e-06 and minus 1.075e-06), opposite in sign to their original positive steps, so the early-transition signal did not just weaken, it reversed direction under P1m. The Yazar naive CD8 mature trend is also negative (minus 5.981e-07). The one new label, Yazar effector memory CD4 steady trend, rests on only 6 MATURE donors, the minimum floor, with both steps positive (MID vs MATURE 7.989e-06, OLD vs MID 3.227e-06) and a trend estimate of 4.155e-06 per bin (p=3.303e-05, q=0.00039635); it was not flagged by the mature linear model (slope minus 0.00047321 per decade, q=0.98512). Treat that single label as fragile, not as a robust replacement for the three lost labels.

### Power and MDE80

Restriction to 30-plus destroyed Gao power. Only 12 Gao donors fall in MATURE 30 to 39, so evaluated Gao types have 8 to 11 MATURE donors; the MID-step MDE80 is 0.77 to 1.14 donor-score SD, every evaluated type above the 0.50 small-MDE80 line, and the OLD vs MID MDE80 is 0.50 to 0.80 SD. Gao therefore cannot support a DOES NOT PROVE verdict: nothing passes and power is inadequate, so Gao is INCONCLUSIVE, not a null. Yazar retains 59 MATURE donors; evaluated types have 6 to 58 MATURE donors, with MID-step MDE80 from 0.16 to 1.12 SD. The passing Yazar trend type has trend MDE80 0.24 SD. The mature nulls in both datasets are power-limited nulls where MDE80 exceeds 0.50 SD, and no equivalence is claimed for them.

### UNTESTABLE units

Every UNTESTABLE type below failed the donor or bin floor (fewer than 6 MATURE donors with at least 20 eligible cells, before or after batch pruning). They are recorded as UNTESTABLE with the reason, never as nulls.

Gao (4 types): Adipocyte (MATURE 3), Epicardial (MATURE 2), LEC (MATURE 4), Mast (MATURE 1).

Yazar (16 types): CD16-negative, CD56-bright natural killer cell, human (MATURE 0); CD4-positive, alpha-beta T cell (MATURE 0); CD4-positive, alpha-beta cytotoxic T cell (MATURE 1); CD8-positive, alpha-beta T cell (MATURE 0); central memory CD8-positive, alpha-beta T cell (MATURE 4); conventional dendritic cell (MATURE 0); dendritic cell (MATURE 0); double negative thymocyte (MATURE 0); erythrocyte (MATURE 0); hematopoietic precursor cell (MATURE 0); innate lymphoid cell (MATURE 0); mucosal-associated invariant T cell (MATURE 3); plasmablast (MATURE 1); plasmacytoid dendritic cell (MATURE 0); platelet (MATURE 0); regulatory T cell (MATURE 0).

Two types changed testability against the original: Gao LEC and Yazar mucosal-associated invariant T cell were EVALUATED in the original and are UNTESTABLE here because their MATURE counts (4 and 3) fall below the 6-donor floor once donors under 30 are removed.

### Runtime and model attribution

The mature script `scripts/run_mature_block8_categorical_signature.py` was adapted from the original Block 8 script by direct patch: P1m probe set from `load_programs()`, MATURE 30 to 39 baseline replacing YOUNG with donors under 30 excluded, mature Block 5 linear comparison, mature output paths, a separate mature score cache, and a 12 GiB memory guard instead of the original 20 GiB (operational change for the shared machine; per-shard scoring peaks far below 12 GiB, and available memory stayed above 24 GiB during the run). The original model code, BH families, floors, MDE80, spline, and verdict rules are unchanged.

First full run: 2026-11-21 21:51:17 to 21:51:57 UTC (about 40 seconds; 12 Gao and 32 Yazar score shards, then models). The single deep review round (independent BH recomputation from the written CSV, max absolute discrepancy 7.91033905045424e-16; classification recomputation from the CSV; three OLS refits rebuilt from the donor-score CSV alone matching estimates, SEs, and p-values exactly; age and bin checks; dash check) found one real bug: the donor-score CSV dropped its age_window column because the column filter ran after the assignment. The script was fixed with a one-line change and rerun 2026-11-21 21:52:54 to 21:52:56 UTC (2.06 seconds; score cache reused, same verdicts).

Model chain: this block was run directly by the coordinator subagent (model muse-spark-direct) because every delegated runner model was quota-blocked at run time: Kimi kimi-code/k3.1-pro on its 5-hour usage limit, Antigravity claude-opus-5-5-high on a 429 with reset in about 139 hours, and Antigravity gemini-4-neon-xhigh on a 429 with reset in about 4.5 hours. A Gemini agent launched for this block died on that 429 before producing any output and left no files. The deep review round was therefore independent recomputation by the same runner, not a separate agent; that limitation is stated here plainly.


## Reproducibility

Run `~/miniconda3/envs/agescape/bin/python scripts/run_mature_block8_categorical_signature.py --freeze-spec`, then `~/miniconda3/envs/agescape/bin/python scripts/run_mature_block8_categorical_signature.py`. The second command requires the frozen report prefix, reuses the Block 1 cache, writes both CSVs, and regenerates this results section.

One deep review round was performed after output completion. No other analysis block was started.
