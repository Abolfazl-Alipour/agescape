# Survivor-bias cross-type correlation test

Date: 2026-11-17. Committed result tables only. No count matrices were loaded. Read-only on the repo.

## Question

Across cell types, do the types losing share with age show less within-type change? Survivor bias says yes: the most-changed cells die, so the surviving cells look unchanged, and the type loses share at the same time.

## Sign convention (read this first)

Share loss is defined as minus the fraction slope, so a positive loss means the type is losing share. Survivor bias predicts that larger share loss goes with SMALLER drift, which is a NEGATIVE correlation between share loss and drift. Written with the raw fraction slope instead, the same prediction is a POSITIVE correlation between fraction slope and drift, because a falling type has a negative fraction slope. Both framings are reported below so the direction cannot be misread.

## Inputs

Fraction slopes: `~/workspace/agescape-repo/results/reviewer-phase1-composition-global-bh.csv`, 431 rows, column `effect` with `effect_type` CLR_slope_per_year (centered log ratio slope per year, the change in a cell type's share within its organ). The 624-row file `rev3-composition-global-bh-624.csv` was checked: it contains all 431 rows with identical effects, and its 193 extra rows are fly, worm, and zebrafish, which have no rows in the drift table used here. The 431-row file is the primary input.

Drift slopes: `~/workspace/agescape-repo/results/h3-cross-species-vector-tests.csv`, filtered to `metric` profile_hellinger and `model` depth_adjusted, column `slope_per_year` (Hellinger distance from the young reference per year; Hellinger distance is a 0 to 1 measure of how different two expression proportion vectors are). The unadjusted model was carried as a sensitivity. The `standardized_effect` column is empty in this file, so raw slopes are the only option.

Join key, exact match, no fuzzy matching: species, organ (called `tissue` in the drift table), cohort, cell_type.

## Join result

Mouse: 62 of 136 composition units joined; all 62 drift units matched. 74 composition units have no drift row, including all of Pancreas, Thymus, and Fat, which are absent from the drift table.

Human: 135 of 291 composition units joined; all 135 drift units matched. 156 composition units have no drift row. Four cohorts in the composition table (Kidney_NormalAdult, RADC_DLPFC, Kidney_SN_v15, SEAAD_DLPFC) are absent from the drift table.

Total joined: 197 units (62 mouse, 135 human). No join problem on keys. No missing slopes in the joined set.

Scopes were kept separate. Mouse and human slopes are both per year, but a year is a very different fraction of lifespan in each species, so no pooled cross-species correlation is reported.

## Results

Primary endpoint: Hellinger drift slope, depth-adjusted model. P-values are two-sided. Permutation p-values use 10,000 label shuffles, seed 20261117.

Mouse, N = 62 units:

| Test | Pearson r | p | Permutation p | Spearman rho | p | Permutation p |
|---|---|---|---|---|---|---|
| Fraction slope vs drift slope | -0.250 | 0.0498 | 0.0559 | -0.153 | 0.2340 | 0.2298 |
| Fraction slope vs absolute drift slope | -0.200 | 0.1192 | 0.1132 | 0.018 | 0.8895 | 0.8893 |
| Share loss vs drift slope | 0.250 | 0.0498 | 0.0488 | 0.153 | 0.2340 | 0.2346 |
| Share loss vs absolute drift slope | 0.200 | 0.1192 | 0.1137 | -0.018 | 0.8895 | 0.8880 |
| Fraction slope vs drift slope, unadjusted model | -0.114 | 0.3770 | 0.3670 | 0.013 | 0.9178 | 0.9112 |

Human, N = 135 units:

| Test | Pearson r | p | Permutation p | Spearman rho | p | Permutation p |
|---|---|---|---|---|---|---|
| Fraction slope vs drift slope | -0.078 | 0.3664 | 0.2122 | -0.151 | 0.0800 | 0.0774 |
| Fraction slope vs absolute drift slope | -0.071 | 0.4124 | 0.2335 | 0.007 | 0.9334 | 0.9301 |
| Share loss vs drift slope | 0.078 | 0.3664 | 0.2119 | 0.151 | 0.0800 | 0.0732 |
| Share loss vs absolute drift slope | 0.071 | 0.4124 | 0.2383 | -0.007 | 0.9334 | 0.9313 |
| Fraction slope vs drift slope, unadjusted model | -0.142 | 0.1013 | 0.0921 | -0.092 | 0.2872 | 0.2880 |

Descriptive split, mean drift slope in units losing share vs units gaining share: mouse losing (41 units) 0.00363 vs gaining (21 units) 0.00297; human losing (67 units) 0.00028 vs gaining (68 units) 0.00001. Mean absolute drift: mouse losing 0.00489 vs gaining 0.00488; human losing 0.00038 vs gaining 0.00018. Losing units drift the same or slightly more, not less.

## Verdict: no support for survivor bias in this test

The share-loss correlations are POSITIVE in point estimate in both species (mouse Pearson 0.25, Spearman 0.15; human Pearson 0.08, Spearman 0.15), which is the OPPOSITE direction from the survivor-bias prediction, and the rank-based tests are null in both species. The single nominally significant result, mouse Pearson p 0.0498, does not survive its permutation p value (0.0559), does not appear in the Spearman test, does not appear in the unadjusted-model sensitivity, and does not appear in the absolute-drift version. It is also likely pulled by one unit: the single mouse liver unit carries the most extreme fraction slope in the joined set (-0.162 per year) with high drift (0.0159 per year), and the rank-based test that downweights it is null. The absolute-drift versions, which test magnitude of change regardless of direction, are null everywhere (all p values above 0.11).

In plain words: types losing share do not show less within-type change. If anything they show slightly more. This cross-type test does not support survivor bias.

## Per-organ description

Per-organ means are in `survivor-bias-per-organ-descriptive.csv`. Organ Ns are small (mouse 1 to 12 units per organ; human 15 to 80). No formal test per organ was run. Nothing in the per-organ means changes the verdict.

## Caveats

This is an ecological correlation across units that share donors, not causal evidence. Units from the same tissue share animals or donors, so the 62 and 135 units are not independent observations, and the p-values are optimistic. The drift table covers only the units that passed its eligibility filters, so the joined set is a selected subset (62 of 136 mouse units, 135 of 291 human units). Composition and drift slopes were estimated from the same donors. A within-type survivor-bias test (do the most-changed cells within one type disappear) needs cell-level data and is not tested here.

## Files

- `survivor-bias-joined-scatter-data.csv`: 197 joined units with fraction slope, drift slope, unadjusted drift slope, share loss, and absolute drift slope.
- `survivor-bias-correlation-results.csv`: 20 rows, all scope, test, method, correlation, p-value, permutation p-value, and N combinations.
- `survivor-bias-per-organ-descriptive.csv`: per-organ descriptive means, no tests.
- `analysis.py`: the exact script used, including the permutation seed.
