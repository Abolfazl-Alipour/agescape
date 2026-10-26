# Project Overview: Conserved Aging Explorer

**Subtitle:** within mouse; cross-species conservation not established

## What this project is

Agescape is a reproducible single-cell transcriptomics analysis and public interactive explorer. It asks a narrow question: where can aging-associated expression changes be estimated credibly, and what happens when a mouse-derived aging signature is evaluated in other cohorts and species?

The project is both a scientific analysis and a data-communication artifact. The public app is https://rezaalipour.com/agescape/. The scientific reports and committed aggregate tables are the source of its claims.

## What it is not

Agescape is not a longitudinal experiment, a perturbation study, a universal map of mammalian aging, or proof that naked mole-rats do not age. It does not provide matched human replication for each mouse lung/kidney cell-type unit. It does not claim that non-significant results demonstrate biological flatness.

## Data and analysis design

| Arm | Data | Replication unit | Purpose |
|---|---|---|---|
| Mouse | Tabula Muris Senis droplet atlas: 245,389 cells, 23 mice, 16 tissues | Mouse | QC, pseudobulk DE, cell-type-unit ranking, published-coefficient agreement |
| Initial human skin | Tabula Sapiens: 4 donors | Donor | Illustrate why the initial design was inadequate |
| Naked mole-rat | Savina et al. 2022 skin data: 5 animals | Animal | Descriptive score trajectory, TOST equivalence, and numeric-cutoff comparison |
| Public human cohorts | 11 CELLxGENE cohorts, 621 scored control donors | Donor | Whole-cohort and exploratory cell-type score regressions |
| Human lifespan atlas | DLPFC atlas: 284 control donors aged 0.2-97 years, 27 cell types | Donor | Lifespan score regressions with normalization and covariate robustness |
| Longevity context | AnAge Build 14, 4,141 species | Species record | Lifespan and body-mass context, not single-cell evidence for every species |

Mouse DE uses donor-level pseudobulk and the 3-month versus 18-month contrast with sex adjustment. Only 17 of 163 tissue-by-cell-type units meet the analysis floor. M9 and M10 use age as a continuous variable. M9 scores per-cell log1p counts; M10 applies log1p to donor-level summed counts. Their slope magnitudes are therefore not directly comparable.

## Main results

1. **Mouse signal is concentrated by atlas eligibility, not shown to be concentrated by biology.** Seventeen units are eligible, all in lung or kidney, and 14 have aging-associated DE genes. Lung adventitial cells have 55 and kidney proximal convoluted tubule cells have 53.
2. **Canonical positive controls do not reach nominal significance in the top unit.** In lung adventitial cells: Apoe (+0.099, p=0.630), Il6 (+0.778, p=0.267), Lgals3 (+0.467, p=0.344), and Igfbp3 (-0.905, p=0.312). Calibration rests on the permutation null (exact Fisher-Pitman p=0.111 marginal, p=0.133 linked animals; null mean 4.6 against 55 observed) and the sex floor.
3. **Effect estimates broadly agree with the published mouse atlas.** Thirteen units agree, four partially agree, and none disagree. Kidney proximal tubule agreement has Spearman rho 0.569 across 12,713 shared genes with 85.1% sign concordance.
4. **Signature hierarchy and external overlap.** P0 catalog contains 257 unique genes. P1_Clean contains 139 genes (82 up, 57 down); P2_Clean_SexTrimmed contains 114 genes; P3_Clean_SexIEGTrimmed contains 110 genes. Covariance-calibrated $Q^*$ identifies 54 genes with significant between-organ heterogeneity. P1_Clean significantly overlaps established aging signatures: Zhang et al. (2021) FACS 330 GAGs (21 shared, 5.83-fold enrichment, p=6.42e-11), Zhang et al. Droplet 59 GAGs (7 shared, 10.87-fold, p=3.37e-06), and Saul et al. (2022) SenMayo 118 (10 shared, 11.31-fold, p=1.68e-08).
5. **Mouse signature rises, robust across specifications.** Spec (a) whole-body 23 mice slope is 0.134 (p=0.015); Spec (c) REML + CR2 cluster-robust across 19 clusters is 0.205 (p=0.00064, df=8.48); Spec (d1) out-of-sample ages slope is 0.107 (p=0.012); Spec (d3) 14 non-Lung/Kidney holdout tissues slope is 0.144 (p=0.0167, df=5.51). Spec (d2) confirms a late-adult plateau (slope -0.0056, p=0.956).
6. **The NMR result is inconclusive.** No rise is detected (slope +0.0058 per normalized lifespan, p=0.9332). Five animals, two group-mean ages, skin-only sampling, and 31% signature coverage cannot establish equivalence or a species difference. TOST equivalence at mouse delta_min=0.2023 yields 90% CI [-0.1450, 0.1567] (p_TOST=0.0274; permutation p=1.0); the direct contrast with mouse whole-body is non-significant (p=0.1153). Whole-epidermis Igfbp3 pseudobulk expression is unchanged (221.8 CPM young vs 223.8 CPM aged).
7. **No supported whole-cohort rise appears in public human data.** MSSM has a nominal decline (p=0.028) with BH q=0.311 across 11 cohorts.
8. **Human lifespan signals depend on cell type and analytical robustness.** In the DLPFC Lifespan Atlas (284 donors), whole-cortex CPM unadjusted slope is +0.00032/yr (p=0.0118); after adjusting for sex, brain bank, PMI, and technical PCs 1-3, the slope attenuates to null (p=0.8626). Glial associations (Astro CPM q=0.0108, Oligo CPM q=0.0019, Adaptive CPM q=0.0412) are single-cohort and exploratory; none retain significance in the full confounder-adjusted adult stratum.

## How the app is kept accountable

`scripts/export_app_data.py` reshapes committed result tables into static JSON. `scripts/check_app_data.py` checks numeric anchors, semantic bindings, source checksums, figure hashes, and data-size limits. The app does not perform statistical analysis in the browser.

Graphs distinguish eligibility, point estimates, confidence intervals, nominal evidence, and unsupported directions. Source drawers provide public aggregate downloads. Controlled-access DLPFC donor-level values are excluded from the repository and app.

## Principal limitations

- Partial sex-by-age confounding and limited tissue eligibility in TMS.
- Cross-sectional rather than longitudinal sampling throughout.
- Small NMR sample, group-mean ages, incomplete signature coverage, and cell-level rather than animal-replicated DE in the Savina comparison.
- Different tissues, compositions, expression summaries, and normalization across human datasets.
- Composition-weighted whole-donor scores and missing-cell-type variation in whole-cortex summaries.
- Exploratory cell-type multiplicity in M9 and attenuation across confounder-adjusted M10 sensitivity arms.
- No direct cross-cell-type heterogeneity test and no common biological cross-species effect scale.

## Reviewer guide

Start with `README.md` for commands, this document for design, and the reports for scientific detail. `HANDOFF.md` describes current operation and deployment. `PLAN.md` and `KICKOFF.md` preserve the original pre-registration and setup instructions; later reports and deviations logs supersede them where the project evolved.
