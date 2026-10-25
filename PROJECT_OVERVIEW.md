# Project Overview

## What this project is

Agescape is a reproducible single-cell transcriptomics analysis and public interactive explorer. It asks a narrow question: where can aging-associated expression changes be estimated credibly, and what happens when a mouse-derived gene score is evaluated in other cohorts and species?

The project is both a scientific analysis and a data-communication artifact. The public app is https://rezaalipour.com/agescape/. The scientific reports and committed aggregate tables are the source of its claims.

## What it is not

Agescape is not a longitudinal experiment, a perturbation study, a universal map of mammalian aging, or proof that naked mole-rats do not age. It does not provide matched human replication for each mouse lung/kidney cell-type unit. It does not claim that non-significant results demonstrate biological flatness.

## Data and analysis design

| Arm | Data | Replication unit | Purpose |
|---|---|---|---|
| Mouse | Tabula Muris Senis droplet atlas: 245,389 cells, 23 mice, 16 tissues | Mouse | QC, pseudobulk DE, cell-type-unit ranking, published-coefficient agreement |
| Initial human skin | Tabula Sapiens: 4 donors | Donor | Illustrate why the initial design was inadequate |
| Naked mole-rat | Savina et al. 2022 skin data: 5 animals | Animal | Descriptive score trajectory and numeric-cutoff comparison |
| Public human cohorts | 11 CELLxGENE cohorts, 621 scored control donors | Donor | Whole-cohort and exploratory cell-type score regressions |
| Human lifespan atlas | DLPFC atlas: 284 control donors aged 0.2-97 years, 27 cell types | Donor | Lifespan score regressions with normalization and covariate robustness |
| Longevity context | AnAge Build 14, 4,141 species | Species record | Lifespan and body-mass context, not single-cell evidence for every species |

Mouse DE uses donor-level pseudobulk and the 3-month versus 18-month contrast with sex adjustment. Only 17 of 163 tissue-by-cell-type units meet the analysis floor. M9 and M10 use age as a continuous variable. M9 scores per-cell log1p counts; M10 applies log1p to donor-level summed counts. Their slope magnitudes are therefore not directly comparable.

## Main results

1. **Mouse signal is concentrated by atlas eligibility, not shown to be concentrated by biology.** Seventeen units are eligible, all in lung or kidney, and 14 have aging-associated DE genes. Lung adventitial cells have 55 and kidney proximal convoluted tubule cells have 53.
2. **Effect estimates broadly agree with the published mouse atlas.** Thirteen units agree, four partially agree, and none disagree. Kidney proximal tubule agreement has Spearman rho 0.866 across 13,019 shared genes.
3. **The catalog is mostly non-recurrent among significant lists.** It contains 261 direction-specific entries across 257 unique genes. Unequal power prevents interpreting non-overlap as proof of private biology.
4. **The mouse score rises, but validation is partly circular.** The slope is 0.235 per normalized-lifespan unit, p=0.008, n=23, in a score derived from the same atlas.
5. **The NMR result is inconclusive.** No rise is detected, but five animals, two group-mean ages, skin-only sampling, and 31% program coverage cannot establish equivalence or a species difference.
6. **No supported whole-cohort rise appears in public human data.** MSSM has a nominal decline (p=0.028) with BH q=0.311 across 11 cohorts.
7. **Human lifespan signals depend on cell type and analytical robustness.** Seven DLPFC cell types have same-direction nominal signals in raw-sum and CPM scoring. Six also have BH q<0.05 in both, but adult-only robustness is not universal. The corrected available-cell-type average detects no rise (p=0.3619).

## How the app is kept accountable

`scripts/export_app_data.py` reshapes committed result tables into static JSON. `scripts/check_app_data.py` checks numeric anchors, semantic bindings, source checksums, figure hashes, and data-size limits. The app does not perform statistical analysis in the browser.

Graphs distinguish eligibility, point estimates, confidence intervals, nominal evidence, and unsupported directions. Source drawers provide public aggregate downloads. Controlled-access DLPFC donor-level values are excluded from the repository and app.

## Principal limitations

- Partial sex-by-age confounding and limited tissue eligibility in TMS.
- Cross-sectional rather than longitudinal sampling throughout.
- Small NMR sample, group-mean ages, incomplete program coverage, and cell-level rather than animal-replicated DE in the Savina comparison.
- Different tissues, compositions, expression summaries, and normalization across human datasets.
- Composition-weighted whole-donor scores and missing-cell-type variation in whole-cortex summaries.
- Exploratory cell-type multiplicity in M9 and incomplete robustness across all M10 sensitivity arms.
- No direct cross-cell-type heterogeneity test and no common biological cross-species effect scale.

## Reviewer guide

Start with `README.md` for commands, this document for design, and the reports for scientific detail. `HANDOFF.md` describes current operation and deployment. `PLAN.md` and `KICKOFF.md` preserve the original pre-registration and setup instructions; later reports and deviations logs supersede them where the project evolved.
