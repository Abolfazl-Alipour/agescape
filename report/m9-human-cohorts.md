# M9 Addendum: Public Human Cohort Analysis

Date: 2026-10-21; revised 2026-10-24. Pipeline: `scripts/run_m9_pipeline.py`. Data: `results/m9-cohort-summary.csv`, `results/m9-regressions.csv`, and aggregate export `app/data/humancohorts.json`.

## Question and data

Does the mouse-derived aging score track chronological age in human donors? M9 analyzes 11 public CELLxGENE cohorts from brain, kidney, and lung: MSSM, HBCC, RADC, SEA-AD, four kidney cohorts, three lung cohorts, and immune-aging BAL. Only control donors with numeric ages are included: 624 selected and 621 scored after coverage filtering. Each cohort has 98.4% program-gene coverage.

The mouse catalog has 261 direction-specific entries across 257 unique genes. Human ortholog mapping yields 121 up and 70 down slots; CKB, ID2, and STMN1 occur in both directional sets and do not simply cancel because the up and down means contain different numbers of genes.

## Methods

Per cell, score = mean log1p(raw count) across detected up-program genes minus the corresponding down-program mean. Donor-by-cell-type scores average cells of that type, requiring at least 20 cells. Whole-donor scores average cell types and therefore depend on sampled composition. OLS regresses score on age with sex adjustment. M9 and M10 use different expression summaries: M9 averages per-cell log1p counts, while M10 applies log1p to donor-level summed counts. Their slope magnitudes are not directly comparable.

The M4 floor of at least 6 donors and age SD at least 12 years is an eligibility rule, not proof of adequate power. BH correction across the 11 whole-cohort tests is reported in `q_bh`. Cell-type p-values remain exploratory and uncorrected across hundreds of tests.

## Whole-cohort results

No cohort has a nominally significant positive whole-cohort score slope. Eight point estimates are positive and three are negative; non-significant intervals do not establish flatness.

MSSM DLPFC (n=213) has slope -0.000193 per year, 95% CI [-0.000366, -0.000021], nominal p=0.028. Across all 11 cohort tests its BH q=0.311, so this decline is not supported after cohort-level multiplicity correction. The separate M10 DLPFC analysis does not reproduce that direction in its available-cell-type average, but different scoring and cohort composition prevent treating this as an exact replication failure.

## Cell-type results

The exploratory ranking contains positive and negative associations. Examples include SEA-AD oligodendrocyte precursor cells and oligodendrocytes in the positive direction, and MSSM L2/3 IT excitatory neurons in the negative direction. These are hypothesis-generating results, not a confirmed hit list. The MSSM neuronal decline does not replicate as a same-direction robust association in M10; M10's glial findings are reported separately with raw, CPM, sex-adjusted, adult-only, and multiple-testing results.

Differences between significant and non-significant lists do not prove opposite underlying effects or cell-type-private biology. Direct heterogeneity tests and comparable sampling would be required.

## Limitations

1. The program was derived in mouse lung and kidney; human brain, kidney, and lung scoring is cross-species and often cross-tissue transport.
2. Control donors can differ across cohorts in recruitment, comorbidities, tissue sampling, sequencing, and cell-type composition.
3. Whole-donor scores are composition-weighted and can move without a shared within-cell-type effect.
4. Donor-level replication avoids treating cells as independent animals, but observational age associations remain vulnerable to cohort confounding.
5. Non-significant slopes are not evidence of equivalence or absence of aging.
6. M10 is a complementary analysis, not an identical replication protocol.

## Bottom line

Public human cohorts provide no supported evidence of a whole-cohort rise in this mouse-derived score. Cell-type estimates vary, but exploratory multiplicity and cohort differences preclude a universal or cell-type-private biological conclusion. The strongest subsequent evidence and its own normalization limitations are documented in `report/m10-dlpfc-lifespan-atlas.md`.
