# Block 8: Nonlinear / categorical age-group signature models

Status: specification frozen before running any Block 8 score or model.

## Goal

Test whether the signature effect is linear in age or concentrated in one age transition. All gene-level and signature tests so far use a linear age slope, which assumes linearity. The drift result itself is a nonlinearity finding (effect concentrated in young-to-mature, dead against the mature reference). This block asks the same question of the signature: is the effect a maturation contrast (YOUNG to MID), an aging contrast (MID to OLD), or a steady trend?

## Inputs

- Gao and Yazar h5ad files (same paths as Block 1), same donor tables, same technical gene mask.
- P1_Clean gene set (139 genes, 89 up, 50 down) from `scripts/run_m5_heterogeneity.py` logic; verify counts before running.
- Scoring code: reuse the composite-score computation from the transport pipeline (`scripts/run_m9_pipeline.py` / follow-up transport scripts). Per-donor per-cell-type composite scores must be computed from full-depth within-cell proportions, the same way the transport models score donors.
- Age bins per shared rule 11 (LT20 excluded from all models in this block).

## Exact analysis specification

- Per donor per cell type (eligibility floors per shared rule 6): P1 composite score = mean up-gene score minus mean down-gene score, donor pseudobulk, full depth, technical mask applied.
- Primary model, per dataset and cell type: OLS of composite score on categorical age bin (YOUNG reference) + sex_male + batch (Gao study indicators; Yazar pool fixed effects). Report MID vs YOUNG and OLD vs YOUNG coefficients with 95 percent CIs and p-values, plus the OLD vs MID contrast (mature-to-old change, computed from the same fit, not a separate model).
- Trend test secondary: ordinal model with bin coded 1, 2, 3 (YOUNG, MID, OLD), same covariates; report slope, CI, p-value. A significant categorical pattern with a null trend test means the effect is a step, not a slope; report that distinction explicitly.
- Unadjusted sensitivity: categorical model without batch covariates.
- Spline sensitivity (labeled, secondary): cubic B-spline on continuous age with 3 degrees of freedom, F test vs the linear model, mirroring the Phase 2 nonlinearity analysis in `reports/reviewer_response_phase2.md` section 2 and `scripts/run_reviewer_phase2_human_age.py`. Read that script and reuse its spline construction.
- Family (state before correction): BH q within each dataset across cell types, separately for the MID vs YOUNG contrast, the OLD vs MID contrast, and the trend test. Three stated families per dataset; do not merge them.

## Fixed implementation details

- The P1_Clean count was verified before model fitting from `results/m5-program-hierarchy.csv`: 139 unique mouse genes, 89 UP and 50 DOWN. `results/m5-gene-level-programs.csv` verifies the 257-gene P0 source set and `results/m5-mouse-program-catalog.csv` was read as the unit-level source summary. One-to-one mouse-to-human mapping uses `load_program_tiers`, giving 76 UP and 42 DOWN human genes. Mapped genes absent from a dataset are omitted from the corresponding directional denominator; detected mapped genes are not filtered by expression level.
- The existing Block 1 donor-sharded cache is reused. It contains raw integer counts after duplicate-symbol collapse and the human technical mask, restricted to post-mask UMI >= 1500 cells. The large h5ad matrices are not reread. The human technical gene mask excludes uppercase `MT-`, `RP[L/S]`, `MRP[L/S]`, `XIST`, `TSIX`, `DDX3Y`, `EIF2S3Y`, `UTY`, and `KDM5D`. ERCC controls are also excluded by the cache construction.
- For each cell and direction, the gene score is the mean within-cell proportion over the mapped genes present in that dataset. For each donor and cell type, up and down scores are arithmetic means of these cell scores across eligible cells, and composite = up minus down. This is an equal-cell-weight donor pseudobulk of full-depth within-cell proportions, with no log transform and no downsampling. Cells below 1500 post-mask UMIs were excluded as an eligibility filter only.
- A donor-cell-type pair requires at least 20 eligible cells. Categorical group models require at least 6 eligible donors in each of YOUNG, MID, and OLD, both before and after batch pruning. At least 6 donors total are required for every fit. The repository covariate pruning requires at least 2 donors per retained batch, full rank, and residual df floors of 5 for Gao and 10 for Yazar in adjusted models and 5 in unadjusted models. Units below a floor are UNTESTABLE, never nulls, and no unadjusted fallback replaces a failed primary fit.
- Age bins reuse `assign_age_bin` unchanged: LT20 under 20, YOUNG 20 to 39, MID 40 to 59, OLD 60 and older. LT20 scores may be present in the donor-score CSV but are excluded from every model. Gao uses intercept plus 6 indicators for its 7 studies with ENCODE v4 (Snyder) as reference when present. Yazar uses intercept plus 74 indicators for its 75 pools. Both include sex_male.
- The OLD vs MID contrast is beta_OLD minus beta_MID from the categorical fit, with its standard error from that fit's covariance matrix. The ordinal trend estimate is per one-bin increment. OLD vs YOUNG is reported descriptively and is not one of the three correction families. BH q-values are assigned only to evaluable primary rows in the three stated families. Unadjusted and spline sensitivities and OLD vs YOUNG have no q-value because the specification defines exactly three families per dataset.
- The spline uses `dmatrix("bs(age, df=3, degree=3, include_intercept=False)", part)` exactly as in `scripts/run_reviewer_phase2_human_age.py`, on the same pruned donors and with the same sex and batch covariates as the primary. The nested linear model uses continuous age in years plus those covariates. The F test compares spline to linear and has 2 numerator degrees of freedom when both designs are full rank.
- For categorical contrasts and trend, MDE80 is the exact two-sided noncentral-t 80 percent power effect, computed as the fitted-contrast standard error times the required noncentrality. It is reported in composite-score units and divided by the donor-score SD in the fitted cell type. A small MDE80 is predefined as at most 0.50 donor-score SD per contrast or ordinal-bin step, matching the repository convention. For the spline F test, MDE80 is the exact noncentral-F noncentrality giving 80 percent power at alpha 0.05; its partial-R2 equivalent is also reported.
- No split-half technical null is required for these donor-score regressions. No random procedure is used, so no random seed is needed. No new downloads are used. Available memory is checked from `/proc/meminfo` MemAvailable before each cache-shard stage; below 20 Gi the process waits and retries every 60 seconds.

## Outputs

- `results/followup-11pt-block8-categorical-signature.csv`: one row per dataset, cell type, contrast or test, model (primary, unadjusted, spline), with estimate, CI, p, q, N donors per bin, MDE80.
- `results/followup-11pt-block8-donor-scores.csv`: per-donor per-cell-type composite scores with age bin, sex, and batch, so the categorical fits are reproducible from the CSV alone.
- `reports/followup-11pt-block8-categorical-signature.md`: spec, per-cell-type classification (maturation-concentrated, aging-concentrated, steady trend, null), and verdicts.
- `scripts/run_followup_11pt_block8_categorical_signature.py`.

## Pre-stated classification and verdict rules

A cell type is classified maturation-concentrated only if MID vs YOUNG passes BH correction and OLD vs MID does not; aging-concentrated only if OLD vs MID passes correction and MID vs YOUNG does not; steady trend only if the trend test passes correction and both steps agree in direction. Transition-concentrated classifications take precedence if more than one rule is met. Both steps agree in direction only when both estimated steps are nonzero and have the same sign. All other evaluable cell types are null. Unevaluable cell types are UNTESTABLE.

Block-level verdict: PROVES (signature effect is nonlinear / transition-concentrated) only if at least one cell type is classified maturation-concentrated or aging-concentrated after correction in the batch-adjusted primary model. DOES NOT PROVE if all passing cell types are steady-trend or if nothing passes with small MDE80. INCONCLUSIVE if nothing passes with large MDE80.

For a no-transition dataset or overall summary, any steady-trend cell type gives DOES NOT PROVE. If no cell type passes any of the three primary families, DOES NOT PROVE requires every evaluated primary MDE80 in those families to be at most 0.50 donor-score SD; otherwise the summary is INCONCLUSIVE. A significant pattern that meets no classification rule is INCONCLUSIVE. The overall summary is PROVES if either dataset is PROVES, DOES NOT PROVE only if both datasets are DOES NOT PROVE, and otherwise INCONCLUSIVE.

The classification will be compared against linear-slope results from `results/followup-organ-17-heart-gao-p1-signature.csv` (`pooled_7study_primary`, `q_bh_composite`) and `results/followup-hcata-pool-sensitivity-p1-signature.csv` (`study=yazar`, `arm=pool_fe`, non-pooled, `q_bh_all_composite`). For each passing cell type, the report states whether the linear model had already flagged it at its source CSV's BH q < 0.05. These existing transport fits use their original log2CPM and depth-adjusted implementation, so the comparison concerns prior flagging, not numerical identity with the new full-depth proportion score.

Exactly one deep review round follows completion of the outputs, followed by fixes, verification, commit, and push.

## Results

All estimates are on the full-depth within-cell-proportion composite scale. MID vs YOUNG and OLD vs YOUNG are coefficients with YOUNG as reference. OLD vs MID is OLD minus MID from the same categorical fit. Trend is per ordinal-bin increment. NA denotes UNTESTABLE or not applicable. BH q-values exist only for the three pre-stated primary families.

### Gao data and family sizes

The cache supplied 1,296 scored donor-cell-type rows, including 0 LT20 rows excluded from models. Scoring used 70 of 76 mapped UP genes and 39 of 42 mapped DOWN genes. Absent UP genes: CD24, FGG, FOXQ1, MT1B, PIGR, PLA2G2D. Absent DOWN genes: EPHB3, MMP13, MMP3. The source cache has 522,761 eligible full-depth cells after the technical mask and 1500-count eligibility filter.

| Primary family | Evaluated cell types | BH q < 0.05 |
| --- | --- | --- |
| MID_vs_YOUNG | 10 | 1 |
| OLD_vs_MID | 10 | 0 |
| trend | 10 | 0 |

### Yazar data and family sizes

The cache supplied 20,858 scored donor-cell-type rows, including 116 LT20 rows excluded from models. Scoring used 76 of 76 mapped UP genes and 42 of 42 mapped DOWN genes. Absent UP genes: none. Absent DOWN genes: none. The source cache has 617,760 eligible full-depth cells after the technical mask and 1500-count eligibility filter.

| Primary family | Evaluated cell types | BH q < 0.05 |
| --- | --- | --- |
| MID_vs_YOUNG | 13 | 2 |
| OLD_vs_MID | 13 | 0 |
| trend | 13 | 3 |

### Primary categorical and trend results

#### Gao

| Cell type | Status | N Y/M/O | MID vs YOUNG estimate (95 percent CI); p; q | OLD vs YOUNG estimate (95 percent CI); p | OLD vs MID estimate (95 percent CI); p; q | Trend estimate (95 percent CI); p; q | MDE80 SD MID / OLD-MID / trend |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Adipocyte | UNTESTABLE_DONOR_OR_BIN_FLOOR | 5/20/12 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| Cardiomyocyte | EVALUATED | 18/51/39 | 3.893e-06 (-2.514e-06 to 1.030e-05); p=0.23083; q=0.28854 | -2.201e-07 (-6.763e-06 to 6.323e-06); p=0.94693 | -4.113e-06 (-9.049e-06 to 8.226e-07); p=0.10139; q=0.9537 | -8.691e-07 (-4.046e-06 to 2.308e-06); p=0.58851; q=0.73563 | 0.80786 / 0.62233 / 0.40059 |
| Endocardial | EVALUATED | 14/29/29 | 1.695e-05 (-6.402e-06 to 4.030e-05); p=0.1519; q=0.28854 | 1.904e-05 (-3.891e-06 to 4.196e-05); p=0.10205 | 2.087e-06 (-1.652e-05 to 2.070e-05); p=0.82337; q=0.99111 | 8.342e-06 (-2.743e-06 to 1.943e-05); p=0.13766; q=0.37782 | 0.80051 / 0.63798 / 0.38006 |
| Endothelial | EVALUATED | 17/39/30 | 7.798e-06 (-4.327e-06 to 1.992e-05); p=0.20422; q=0.28854 | 3.165e-06 (-9.368e-06 to 1.570e-05); p=0.61657 | -4.633e-06 (-1.470e-05 to 5.437e-06); p=0.36251; q=0.9537 | 7.479e-07 (-5.425e-06 to 6.920e-06); p=0.81005; q=0.81005 | 0.74411 / 0.61798 / 0.37884 |
| Epicardial | UNTESTABLE_DONOR_OR_BIN_FLOOR | 4/8/8 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| Fibroblast | EVALUATED | 18/49/38 | 3.752e-05 (-1.287e-05 to 8.790e-05); p=0.14265; q=0.28854 | 4.588e-05 (-5.512e-06 to 9.727e-05); p=0.079553 | 8.365e-06 (-3.099e-05 to 4.772e-05); p=0.67403; q=0.9629 | 2.033e-05 (-4.433e-06 to 4.510e-05); p=0.10645; q=0.37782 | 0.68424 / 0.5345 / 0.33638 |
| LEC | EVALUATED | 6/17/13 | 1.282e-05 (-1.490e-05 to 4.054e-05); p=0.3526; q=0.39178 | 7.103e-06 (-2.033e-05 to 3.454e-05); p=0.60087 | -5.714e-06 (-2.615e-05 to 1.472e-05); p=0.57222; q=0.9537 | 1.544e-06 (-1.142e-05 to 1.451e-05); p=0.8097; q=0.81005 | 1.3162 / 0.97025 / 0.61569 |
| Lymphoid | EVALUATED | 12/25/22 | 1.993e-05 (-9.213e-06 to 4.907e-05); p=0.17579; q=0.28854 | 1.955e-05 (-9.846e-06 to 4.895e-05); p=0.18775 | -3.790e-07 (-2.440e-05 to 2.364e-05); p=0.97485; q=0.99111 | 8.384e-06 (-5.964e-06 to 2.273e-05); p=0.24634; q=0.37782 | 0.85436 / 0.70409 / 0.42067 |
| Mast | UNTESTABLE_DONOR_OR_BIN_FLOOR | 2/8/8 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| Myeloid | EVALUATED | 18/49/39 | 2.704e-05 (9.503e-06 to 4.457e-05); p=0.0028581; q=0.028581 | 2.143e-05 (3.547e-06 to 3.931e-05); p=0.019349 | -5.607e-06 (-1.915e-05 to 7.935e-06); p=0.41321; q=0.9537 | 7.656e-06 (-1.185e-06 to 1.650e-05); p=0.088854; q=0.37782 | 0.57683 / 0.44553 / 0.29086 |
| Neuronal | EVALUATED | 13/28/19 | 3.432e-05 (1.927e-06 to 6.671e-05); p=0.038267; q=0.19134 | 2.534e-05 (-8.881e-06 to 5.955e-05); p=0.14336 | -8.985e-06 (-3.575e-05 to 1.778e-05); p=0.50355; q=0.9537 | 9.614e-06 (-7.484e-06 to 2.671e-05); p=0.26447; q=0.37782 | 1.016 / 0.8395 / 0.53633 |
| Pericyte | EVALUATED | 18/44/33 | 5.299e-06 (-2.691e-05 to 3.751e-05); p=0.74447; q=0.74447 | 1.825e-05 (-1.481e-05 to 5.132e-05); p=0.27551 | 1.296e-05 (-1.339e-05 to 3.931e-05); p=0.33113; q=0.9537 | 9.682e-06 (-6.385e-06 to 2.575e-05); p=0.23431; q=0.37782 | 0.59879 / 0.48985 / 0.29869 |
| vSMC | EVALUATED | 16/40/30 | 1.800e-05 (-7.585e-06 to 4.358e-05); p=0.16531; q=0.28854 | 1.811e-05 (-8.033e-06 to 4.426e-05); p=0.1718 | 1.156e-07 (-2.049e-05 to 2.072e-05); p=0.99111; q=0.99111 | 7.663e-06 (-5.049e-06 to 2.037e-05); p=0.23377; q=0.37782 | 0.84881 / 0.68378 / 0.42182 |

#### Yazar

| Cell type | Status | N Y/M/O | MID vs YOUNG estimate (95 percent CI); p; q | OLD vs YOUNG estimate (95 percent CI); p | OLD vs MID estimate (95 percent CI); p; q | Trend estimate (95 percent CI); p; q | MDE80 SD MID / OLD-MID / trend |
| --- | --- | --- | --- | --- | --- | --- | --- |
| CD14-low, CD16-positive monocyte | EVALUATED | 13/41/197 | -9.851e-06 (-4.158e-05 to 2.188e-05); p=0.54106; q=0.76288 | -1.532e-07 (-2.855e-05 to 2.825e-05); p=0.99152 | 9.697e-06 (-7.761e-06 to 2.716e-05); p=0.27467; q=0.52468 | 3.700e-06 (-7.863e-06 to 1.526e-05); p=0.52872; q=0.57348 | 0.78882 / 0.43405 / 0.28747 |
| CD14-positive monocyte | EVALUATED | 39/87/295 | -3.840e-06 (-1.773e-05 to 1.004e-05); p=0.58683; q=0.76288 | 1.144e-06 (-1.140e-05 to 1.369e-05); p=0.85785 | 4.984e-06 (-4.123e-06 to 1.409e-05); p=0.28252; q=0.52468 | 1.802e-06 (-3.827e-06 to 7.431e-06); p=0.52936; q=0.57348 | 0.43007 / 0.28207 / 0.17435 |
| CD16-negative, CD56-bright natural killer cell, human | UNTESTABLE_DONOR_OR_BIN_FLOOR | 1/1/4 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| CD4-positive, alpha-beta T cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| CD4-positive, alpha-beta cytotoxic T cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 2/13/78 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| CD8-positive, alpha-beta T cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| central memory CD4-positive, alpha-beta T cell | EVALUATED | 99/191/660 | 1.249e-05 (5.535e-06 to 1.944e-05); p=0.00044456; q=0.0057793 | 1.827e-05 (1.207e-05 to 2.447e-05); p=1.027e-08 | 5.786e-06 (9.536e-07 to 1.062e-05); p=0.018995; q=0.12347 | 8.338e-06 (5.486e-06 to 1.119e-05); p=1.326e-08; q=1.724e-07 | 0.18815 / 0.13079 / 0.077207 |
| central memory CD8-positive, alpha-beta T cell | UNTESTABLE_POST_PRUNE_BIN_FLOOR | 5/17/68 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| conventional dendritic cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/5 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| dendritic cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| double negative thymocyte | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| effector memory CD4-positive, alpha-beta T cell | EVALUATED | 16/49/195 | 2.192e-05 (-3.777e-06 to 4.761e-05); p=0.094132; q=0.24746 | 2.590e-05 (1.930e-06 to 4.988e-05); p=0.034337 | 3.986e-06 (-1.137e-05 to 1.935e-05); p=0.60943; q=0.88587 | 9.827e-06 (-5.038e-07 to 2.016e-05); p=0.062154; q=0.19225 | 0.56764 / 0.33931 / 0.22824 |
| effector memory CD8-positive, alpha-beta T cell | EVALUATED | 55/136/552 | 1.769e-06 (-1.286e-05 to 1.640e-05); p=0.81243; q=0.87409 | 3.679e-06 (-9.443e-06 to 1.680e-05); p=0.58215 | 1.910e-06 (-7.223e-06 to 1.104e-05); p=0.68144; q=0.88587 | 1.861e-06 (-3.875e-06 to 7.598e-06); p=0.52427; q=0.57348 | 0.33335 / 0.20812 / 0.13072 |
| erythrocyte | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| gamma-delta T cell | EVALUATED | 15/24/81 | -3.464e-05 (-7.282e-05 to 3.546e-06); p=0.074812; q=0.24746 | -3.582e-05 (-7.209e-05 to 4.418e-07); p=0.05279 | -1.186e-06 (-3.084e-05 to 2.847e-05); p=0.9368; q=0.9851 | -1.503e-05 (-3.249e-05 to 2.432e-06); p=0.090637; q=0.19638 | 0.78305 / 0.60809 / 0.35811 |
| hematopoietic precursor cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| innate lymphoid cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| memory B cell | EVALUATED | 24/49/170 | 6.731e-06 (-1.391e-05 to 2.737e-05); p=0.5208; q=0.76288 | 7.554e-06 (-1.119e-05 to 2.630e-05); p=0.42775 | 8.225e-07 (-1.337e-05 to 1.502e-05); p=0.90913; q=0.9851 | 3.043e-06 (-5.536e-06 to 1.162e-05); p=0.48502; q=0.57348 | 0.54757 / 0.37661 / 0.2276 |
| mucosal-associated invariant T cell | EVALUATED | 6/9/15 | 8.803e-06 (-8.966e-05 to 0.00010727); p=0.85206; q=0.87409 | 2.709e-05 (-8.840e-05 to 0.00014259); p=0.62577 | 1.829e-05 (-6.008e-05 to 9.666e-05); p=0.62753; q=0.88587 | 1.447e-05 (-3.987e-05 to 6.880e-05); p=0.58167; q=0.58167 | 2.1398 / 1.703 / 1.1818 |
| naive B cell | EVALUATED | 21/50/248 | 1.057e-05 (-1.109e-05 to 3.222e-05); p=0.33756; q=0.73139 | 1.044e-05 (-9.466e-06 to 3.034e-05); p=0.30277 | -1.279e-07 (-1.359e-05 to 1.334e-05); p=0.9851; q=0.9851 | 3.519e-06 (-5.201e-06 to 1.224e-05); p=0.42764; q=0.57348 | 0.46438 / 0.28878 / 0.18701 |
| naive thymus-derived CD4-positive, alpha-beta T cell | EVALUATED | 89/163/564 | 8.906e-06 (2.849e-06 to 1.496e-05); p=0.0040098; q=0.026064 | 1.453e-05 (9.108e-06 to 1.995e-05); p=1.887e-07 | 5.625e-06 (1.315e-06 to 9.935e-06); p=0.010605; q=0.12347 | 6.899e-06 (4.382e-06 to 9.417e-06); p=9.994e-08; q=6.496e-07 | 0.22876 / 0.16278 / 0.095073 |
| naive thymus-derived CD8-positive, alpha-beta T cell | EVALUATED | 73/99/208 | 7.709e-06 (-1.353e-06 to 1.677e-05); p=0.095176; q=0.24746 | 1.355e-05 (4.933e-06 to 2.216e-05); p=0.002149 | 5.836e-06 (-1.679e-06 to 1.335e-05); p=0.12752; q=0.52468 | 6.655e-06 (2.450e-06 to 1.086e-05); p=0.0020172; q=0.0087412 | 0.28474 / 0.23612 / 0.13212 |
| natural killer cell | EVALUATED | 66/147/598 | 4.341e-06 (-7.843e-06 to 1.652e-05); p=0.48448; q=0.76288 | 8.763e-06 (-2.037e-06 to 1.956e-05); p=0.11162 | 4.422e-06 (-3.534e-06 to 1.238e-05); p=0.27555; q=0.52468 | 4.393e-06 (-4.263e-07 to 9.212e-06); p=0.073943; q=0.19225 | 0.29969 / 0.19571 / 0.11854 |
| plasmablast | UNTESTABLE_DONOR_OR_BIN_FLOOR | 1/5/8 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| plasmacytoid dendritic cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/1/2 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| platelet | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| regulatory T cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/10/59 | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA | NA (NA to NA); p=NA; q=NA | NA (NA to NA); p=NA; q=NA | NA / NA / NA |
| transitional stage B cell | EVALUATED | 22/36/142 | 1.990e-06 (-2.278e-05 to 2.676e-05); p=0.87409; q=0.87409 | 1.391e-05 (-8.890e-06 to 3.670e-05); p=0.22998 | 1.192e-05 (-6.635e-06 to 3.047e-05); p=0.20634; q=0.52468 | 7.911e-06 (-2.876e-06 to 1.870e-05); p=0.14941; q=0.27748 | 0.57723 / 0.43235 / 0.2514 |

### Unadjusted categorical sensitivity

#### Gao

| Cell type | Status | N Y/M/O | MID vs YOUNG | OLD vs YOUNG | OLD vs MID |
| --- | --- | --- | --- | --- | --- |
| Adipocyte | UNTESTABLE_DONOR_OR_BIN_FLOOR | 5/20/12 | NA; p=NA | NA; p=NA | NA; p=NA |
| Cardiomyocyte | EVALUATED | 18/52/39 | 3.856e-06; p=0.20917 | 1.654e-07; p=0.95874 | -3.690e-06; p=0.11925 |
| Endocardial | EVALUATED | 14/29/29 | 3.163e-06; p=0.8157 | 1.070e-05; p=0.43251 | 7.540e-06; p=0.49027 |
| Endothelial | EVALUATED | 17/40/30 | 3.418e-06; p=0.62101 | 1.373e-06; p=0.85045 | -2.045e-06; p=0.72073 |
| Epicardial | UNTESTABLE_DONOR_OR_BIN_FLOOR | 4/8/8 | NA; p=NA | NA; p=NA | NA; p=NA |
| Fibroblast | EVALUATED | 19/49/38 | 5.346e-05; p=0.063067 | 4.167e-05; p=0.16262 | -1.179e-05; p=0.60218 |
| LEC | EVALUATED | 6/17/13 | 8.263e-06; p=0.60041 | 6.125e-06; p=0.70586 | -2.139e-06; p=0.84915 |
| Lymphoid | EVALUATED | 12/25/22 | 1.099e-05; p=0.5353 | 1.613e-05; p=0.37253 | 5.136e-06; p=0.72358 |
| Mast | UNTESTABLE_DONOR_OR_BIN_FLOOR | 2/8/8 | NA; p=NA | NA; p=NA | NA; p=NA |
| Myeloid | EVALUATED | 18/49/39 | 1.726e-05; p=0.14239 | 1.414e-05; p=0.24454 | -3.116e-06; p=0.73005 |
| Neuronal | EVALUATED | 13/28/19 | 2.888e-05; p=0.077293 | 2.278e-05; p=0.18862 | -6.094e-06; p=0.65205 |
| Pericyte | EVALUATED | 18/44/33 | 1.939e-05; p=0.3782 | 2.059e-05; p=0.36891 | 1.201e-06; p=0.94648 |
| vSMC | EVALUATED | 16/40/30 | 2.179e-05; p=0.093486 | 1.898e-05; p=0.15888 | -2.804e-06; p=0.78753 |

#### Yazar

| Cell type | Status | N Y/M/O | MID vs YOUNG | OLD vs YOUNG | OLD vs MID |
| --- | --- | --- | --- | --- | --- |
| CD14-low, CD16-positive monocyte | EVALUATED | 14/44/204 | 8.011e-06; p=0.65314 | 7.652e-07; p=0.96194 | -7.246e-06; p=0.45344 |
| CD14-positive monocyte | EVALUATED | 40/88/299 | 3.008e-07; p=0.97327 | 6.713e-06; p=0.39708 | 6.412e-06; p=0.26169 |
| CD16-negative, CD56-bright natural killer cell, human | UNTESTABLE_DONOR_OR_BIN_FLOOR | 1/1/4 | NA; p=NA | NA; p=NA | NA; p=NA |
| CD4-positive, alpha-beta T cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA; p=NA | NA; p=NA | NA; p=NA |
| CD4-positive, alpha-beta cytotoxic T cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 2/13/78 | NA; p=NA | NA; p=NA | NA; p=NA |
| CD8-positive, alpha-beta T cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA; p=NA | NA; p=NA | NA; p=NA |
| central memory CD4-positive, alpha-beta T cell | EVALUATED | 99/191/660 | 1.864e-05; p=0.0035847 | 3.507e-05; p=4.157e-10 | 1.643e-05; p=0.00011241 |
| central memory CD8-positive, alpha-beta T cell | EVALUATED | 7/22/80 | 1.897e-05; p=0.37697 | 4.449e-05; p=0.02379 | 2.552e-05; p=0.034328 |
| conventional dendritic cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/5 | NA; p=NA | NA; p=NA | NA; p=NA |
| dendritic cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA; p=NA | NA; p=NA | NA; p=NA |
| double negative thymocyte | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA; p=NA | NA; p=NA | NA; p=NA |
| effector memory CD4-positive, alpha-beta T cell | EVALUATED | 17/50/201 | 3.242e-05; p=0.068798 | 5.197e-05; p=0.0012816 | 1.955e-05; p=0.051214 |
| effector memory CD8-positive, alpha-beta T cell | EVALUATED | 56/136/552 | 2.863e-06; p=0.77197 | 1.928e-05; p=0.027363 | 1.642e-05; p=0.0059422 |
| erythrocyte | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA; p=NA | NA; p=NA | NA; p=NA |
| gamma-delta T cell | EVALUATED | 15/25/88 | -1.302e-05; p=0.56654 | -1.012e-06; p=0.95881 | 1.200e-05; p=0.44726 |
| hematopoietic precursor cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA; p=NA | NA; p=NA | NA; p=NA |
| innate lymphoid cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA; p=NA | NA; p=NA | NA; p=NA |
| memory B cell | EVALUATED | 25/52/176 | 1.224e-05; p=0.35766 | 2.743e-05; p=0.019466 | 1.520e-05; p=0.07873 |
| mucosal-associated invariant T cell | EVALUATED | 10/12/23 | 2.065e-05; p=0.47101 | 4.667e-05; p=0.069846 | 2.602e-05; p=0.2763 |
| naive B cell | EVALUATED | 21/50/251 | 1.230e-05; p=0.48743 | 7.793e-06; p=0.61476 | -4.509e-06; p=0.66847 |
| naive thymus-derived CD4-positive, alpha-beta T cell | EVALUATED | 89/163/564 | 1.266e-05; p=0.0094752 | 2.530e-05; p=2.837e-09 | 1.264e-05; p=0.00012675 |
| naive thymus-derived CD8-positive, alpha-beta T cell | EVALUATED | 77/102/210 | 1.821e-05; p=0.007078 | 2.645e-05; p=1.157e-05 | 8.237e-06; p=0.12646 |
| natural killer cell | EVALUATED | 66/148/598 | 5.160e-06; p=0.54138 | 2.274e-05; p=0.0021943 | 1.758e-05; p=0.00082798 |
| plasmablast | UNTESTABLE_DONOR_OR_BIN_FLOOR | 1/5/8 | NA; p=NA | NA; p=NA | NA; p=NA |
| plasmacytoid dendritic cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/1/2 | NA; p=NA | NA; p=NA | NA; p=NA |
| platelet | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/0/0 | NA; p=NA | NA; p=NA | NA; p=NA |
| regulatory T cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | 0/10/59 | NA; p=NA | NA; p=NA | NA; p=NA |
| transitional stage B cell | EVALUATED | 26/39/147 | 2.244e-05; p=0.16211 | 2.780e-05; p=0.039692 | 5.362e-06; p=0.63825 |

### Cubic B-spline sensitivity

| Dataset | Cell type | Status | F | p | df num | df den | MDE80 noncentrality | MDE80 partial R2 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gao | Adipocyte | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| gao | Cardiomyocyte | EVALUATED | 1.3612 | 0.26115 | 2 | 98 | 9.9354 | 0.092049 |
| gao | Endocardial | EVALUATED | 0.7682 | 0.46821 | 2 | 62 | 10.116 | 0.14027 |
| gao | Endothelial | EVALUATED | 0.68582 | 0.50673 | 2 | 77 | 10.02 | 0.11514 |
| gao | Epicardial | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| gao | Fibroblast | EVALUATED | 0.14681 | 0.86365 | 2 | 95 | 9.9451 | 0.094765 |
| gao | LEC | EVALUATED | 0.54491 | 0.58572 | 2 | 29 | 10.703 | 0.26958 |
| gao | Lymphoid | EVALUATED | 0.76762 | 0.46951 | 2 | 50 | 10.236 | 0.16993 |
| gao | Mast | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| gao | Myeloid | EVALUATED | 0.72475 | 0.48708 | 2 | 96 | 9.9418 | 0.093842 |
| gao | Neuronal | EVALUATED | 2.9402 | 0.061873 | 2 | 51 | 10.224 | 0.16699 |
| gao | Pericyte | EVALUATED | 0.12265 | 0.88473 | 2 | 86 | 9.9783 | 0.10396 |
| gao | vSMC | EVALUATED | 0.27169 | 0.76282 | 2 | 77 | 10.02 | 0.11514 |
| yazar | CD14-low, CD16-positive monocyte | EVALUATED | 5.1556 | 0.006579 | 2 | 195 | 9.7843 | 0.047778 |
| yazar | CD14-positive monocyte | EVALUATED | 0.33316 | 0.71688 | 2 | 359 | 9.7155 | 0.02635 |
| yazar | CD16-negative, CD56-bright natural killer cell, human | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | CD4-positive, alpha-beta T cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | CD4-positive, alpha-beta cytotoxic T cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | CD8-positive, alpha-beta T cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | central memory CD4-positive, alpha-beta T cell | EVALUATED | 1.1576 | 0.31471 | 2 | 871 | 9.6679 | 0.010978 |
| yazar | central memory CD8-positive, alpha-beta T cell | UNTESTABLE_POST_PRUNE_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | conventional dendritic cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | dendritic cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | double negative thymocyte | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | effector memory CD4-positive, alpha-beta T cell | EVALUATED | 2.3438 | 0.098501 | 2 | 207 | 9.7755 | 0.045095 |
| yazar | effector memory CD8-positive, alpha-beta T cell | EVALUATED | 0.75758 | 0.4692 | 2 | 665 | 9.6782 | 0.014345 |
| yazar | erythrocyte | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | gamma-delta T cell | EVALUATED | 1.0663 | 0.34905 | 2 | 81 | 10 | 0.10989 |
| yazar | hematopoietic precursor cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | innate lymphoid cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | memory B cell | EVALUATED | 1.683 | 0.18856 | 2 | 191 | 9.7874 | 0.048745 |
| yazar | mucosal-associated invariant T cell | EVALUATED | 1.1971 | 0.32935 | 2 | 15 | 11.847 | 0.44128 |
| yazar | naive B cell | EVALUATED | 0.080776 | 0.92242 | 2 | 273 | 9.7412 | 0.034453 |
| yazar | naive thymus-derived CD4-positive, alpha-beta T cell | EVALUATED | 0.84861 | 0.42843 | 2 | 737 | 9.674 | 0.012956 |
| yazar | naive thymus-derived CD8-positive, alpha-beta T cell | EVALUATED | 2.1109 | 0.12285 | 2 | 314 | 9.7272 | 0.030048 |
| yazar | natural killer cell | EVALUATED | 2.1255 | 0.12011 | 2 | 733 | 9.6742 | 0.013026 |
| yazar | plasmablast | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | plasmacytoid dendritic cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | platelet | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | regulatory T cell | UNTESTABLE_DONOR_OR_BIN_FLOOR | NA | NA | NA | NA | NA | NA |
| yazar | transitional stage B cell | EVALUATED | 0.68055 | 0.50789 | 2 | 150 | 9.8297 | 0.061501 |

### Classification and comparison with prior linear models

Transition classifications require the stated BH pattern in the batch-adjusted primary model. A significant categorical step with a null trend is a step pattern, not a steady slope. Linear flagged means the source transport CSV had composite BH q < 0.05 in its original linear model.

#### Gao

| Cell type | Classification | Prior linear slope per decade | Prior linear q | Linear flagged |
| --- | --- | --- | --- | --- |
| Adipocyte | UNTESTABLE | -0.094066 | 0.20698 | no |
| Cardiomyocyte | null | 0.027703 | 0.23164 | no |
| Endocardial | null | 0.087917 | 0.015008 | yes |
| Endothelial | null | 0.030257 | 0.20698 | no |
| Epicardial | UNTESTABLE | 0.012795 | 0.78847 | no |
| Fibroblast | null | 0.034834 | 0.20698 | no |
| LEC | null | 0.050581 | 0.39667 | no |
| Lymphoid | null | 0.032978 | 0.23959 | no |
| Mast | UNTESTABLE | 0.07374 | 0.2551 | no |
| Myeloid | maturation-concentrated | 0.037938 | 0.1545 | no |
| Neuronal | null | 0.01091 | 0.7804 | no |
| Pericyte | null | 0.040029 | 0.13842 | no |
| vSMC | null | 0.017588 | 0.55938 | no |

#### Yazar

| Cell type | Classification | Prior linear slope per decade | Prior linear q | Linear flagged |
| --- | --- | --- | --- | --- |
| CD14-low, CD16-positive monocyte | null | 0.034808 | 0.10028 | no |
| CD14-positive monocyte | null | 0.022517 | 0.031382 | yes |
| CD16-negative, CD56-bright natural killer cell, human | UNTESTABLE | NA | NA | not tested |
| CD4-positive, alpha-beta T cell | UNTESTABLE | NA | NA | not tested |
| CD4-positive, alpha-beta cytotoxic T cell | UNTESTABLE | 0.036131 | 0.12856 | no |
| CD8-positive, alpha-beta T cell | UNTESTABLE | NA | NA | not tested |
| central memory CD4-positive, alpha-beta T cell | maturation-concentrated | 0.023014 | 3.382e-06 | yes |
| central memory CD8-positive, alpha-beta T cell | UNTESTABLE | 0.0092981 | 0.72555 | no |
| conventional dendritic cell | UNTESTABLE | NA | NA | not tested |
| dendritic cell | UNTESTABLE | NA | NA | not tested |
| double negative thymocyte | UNTESTABLE | NA | NA | not tested |
| effector memory CD4-positive, alpha-beta T cell | null | 0.020298 | 0.10077 | no |
| effector memory CD8-positive, alpha-beta T cell | null | 0.016432 | 0.019475 | yes |
| erythrocyte | UNTESTABLE | NA | NA | not tested |
| gamma-delta T cell | null | 0.033498 | 0.062151 | no |
| hematopoietic precursor cell | UNTESTABLE | NA | NA | not tested |
| innate lymphoid cell | UNTESTABLE | NA | NA | not tested |
| memory B cell | null | 0.037917 | 0.01011 | yes |
| mucosal-associated invariant T cell | null | 0.025891 | 0.45501 | no |
| naive B cell | null | 0.028607 | 0.011149 | yes |
| naive thymus-derived CD4-positive, alpha-beta T cell | maturation-concentrated | 0.026924 | 1.413e-05 | yes |
| naive thymus-derived CD8-positive, alpha-beta T cell | steady trend | 0.036442 | 0.019929 | yes |
| natural killer cell | null | 0.016819 | 0.0026521 | yes |
| plasmablast | UNTESTABLE | NA | NA | not tested |
| plasmacytoid dendritic cell | UNTESTABLE | NA | NA | not tested |
| platelet | UNTESTABLE | NA | NA | not tested |
| regulatory T cell | UNTESTABLE | 0.020819 | 0.26926 | no |
| transitional stage B cell | null | -0.0044494 | 0.79244 | no |

- Gao verdict: PROVES. Classifications: 1 maturation-concentrated, 0 aging-concentrated, 0 steady trend, 9 null, 3 UNTESTABLE.
- Yazar verdict: PROVES. Classifications: 2 maturation-concentrated, 0 aging-concentrated, 1 steady trend, 10 null, 15 UNTESTABLE.

Block-level verdict: PROVES. This verdict applies to the P1_Clean score in Gao heart and Yazar blood under the pre-stated categorical and trend rules.

The unadjusted and spline models are sensitivities and do not determine the primary classification. The spline F test addresses departure from continuous linearity but is not substituted for the three categorical BH families.

## Reproducibility

Run `~/miniconda3/envs/agescape/bin/python scripts/run_followup_11pt_block8_categorical_signature.py --freeze-spec`, then `~/miniconda3/envs/agescape/bin/python scripts/run_followup_11pt_block8_categorical_signature.py`. The second command requires the frozen report prefix, reuses the Block 1 cache, writes both CSVs, and regenerates this results section.

One deep review round was performed after output completion. No other analysis block was started.
