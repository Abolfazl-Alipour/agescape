# Donor-spread variance test: distance-variance version of donor-cloud spread

Date: 2026-11-17. Analysis of committed donor tables only. No count matrix
was loaded. Read-only on the repo.

## What this test is and is not

Question: are old donors MORE variable than young donors in how far their
expression profiles sit from the young reference?

This is the DISTANCE-VARIANCE version of donor-cloud spread. It compares the
spread (variance) of one number per donor, the Hellinger distance from that
donor's profile to the young reference, across age bins. It is NOT the full
profile-cloud test, which would compare how far old donors sit from each
other in full profile space. That test needs count matrices and was not run
here.

A pooled distance also carries cell-type composition differences between
donors, so spread here includes composition variation, not only within-type
profile spread. The per-cell-type supplement below checks the same question
inside single cell types.

## Inputs and data handling

- Gao heart: `results/followup-between-subject-gao-donors.csv`, pooled rows
  (`All_LV_Nuclei`), one row per donor, 117 donors.
- Yazar blood: `results/followup-between-subject-yazar-donors.csv`, pooled
  rows (`All_Cells`), one row per donor, 981 donors.
- The files are donor-by-cell-type rows, not one row per donor. Only the
  pooled-group rows were used, so the donor is the unit. Verified: one pooled
  row per donor in both files.
- Age bins come from the `bin` column: YOUNG ages 20 to 39, MID ages 40 to
  59, OLD ages 60 and up. Yazar has 5 extra donors aged 19 in bin LT20; they
  appear in the descriptives only and in no comparison.
- Primary distance: `a1_dist_full_1500`, Hellinger distance to the young
  reference, full depth, cells with at least 1500 counts.
- Gao: 4 donors lack the primary distance (1 MID, 3 OLD), matching the
  eligibility flag exactly. Analysis N: YOUNG 20, MID 52, OLD 41.
- Yazar: no missing primary distances. Analysis N: YOUNG 103, MID 193,
  OLD 680.
- Mature-ref distance (`a1_dist_mature_ref_mid_loo`) is missing for ALL young
  donors in both datasets by construction, because the reference is built
  from MID donors. The mature-ref sensitivity therefore compares OLD vs MID,
  never OLD vs YOUNG.

## Methods

Per age bin: N, mean, standard deviation (SD), interquartile range (IQR) of
the distance. Variance comparison: Levene test, median-centered
(Brown-Forsythe), from scipy, plus the variance ratio with a donor-bootstrap
95 percent confidence interval (resample donors within bins, 2000 resamples,
fixed seed 7). Variance ratio above 1 means the first group is more spread
out. Yazar pool-residualized version: distance minus its pool mean plus the
grand mean, using all 75 pools (sizes 9 to 14, no singleton pools, no donors
excluded).

## Results: Gao heart

Per-bin SD of distance to the young reference:

- YOUNG: N 20, mean 0.6419, SD 0.0619, IQR 0.0634
- MID: N 52, mean 0.6128, SD 0.0422, IQR 0.0453
- OLD: N 41, mean 0.6325, SD 0.0549, IQR 0.0798

Comparisons:

- OLD vs YOUNG: variance ratio 0.79, 95 percent CI 0.33 to 2.87, Levene
  p 0.85. No evidence of a difference.
- OLD vs MID: variance ratio 1.70, 95 percent CI 0.79 to 3.26, Levene
  p 0.16. Not significant.
- Mature-ref OLD vs MID: variance ratio 1.68, 95 percent CI 0.83 to 3.11,
  Levene p 0.13. Not significant.

Verdict, Gao: FLAT. No evidence that between-donor spread in
distance-to-young-reference differs with age. Precision is limited by only
20 young donors; the OLD vs YOUNG interval is wide and cannot rule out a
large ratio in either direction.

## Results: Yazar blood

Per-bin SD of distance to the young reference:

- YOUNG: N 103, mean 0.5835, SD 0.0092, IQR 0.0101
- MID: N 193, mean 0.5826, SD 0.0111, IQR 0.0158
- OLD: N 680, mean 0.5803, SD 0.0117, IQR 0.0201
- LT20: N 5, mean 0.5865, SD 0.0046 (descriptive only)

Comparisons, raw distance:

- OLD vs YOUNG: variance ratio 1.62, 95 percent CI 1.23 to 2.29, Levene
  p 1.9e-07. Significant: old donors are more spread out.
- OLD vs MID: variance ratio 1.11, 95 percent CI 0.94 to 1.36, Levene
  p 0.0008. Mixed: the Levene test is significant but the bootstrap interval
  includes 1, and the point estimate is small. The two methods test related
  but different quantities; do not read this as a strong effect.
- Mature-ref OLD vs MID: variance ratio 1.13, 95 percent CI 0.95 to 1.40,
  Levene p 0.0005. Same mixed pattern as the young-ref OLD vs MID result.

Comparisons, pool-residualized distance (pool means removed first):

- Per-bin SD after residualizing: YOUNG 0.0038, MID 0.0037, OLD 0.0043.
  Pool adjustment shrinks every SD by roughly 60 percent, showing how much
  of the raw spread is pool structure.
- OLD vs YOUNG: variance ratio 1.25, 95 percent CI 0.94 to 1.76, Levene
  p 0.43. Not significant.
- OLD vs MID: variance ratio 1.29, 95 percent CI 0.97 to 1.72, Levene
  p 0.058. Not significant.

Verdict, Yazar: the raw increase in spread with age DOES NOT SURVIVE pool
adjustment. On the raw distance, old donors look more variable. Once pool
means are removed, the difference is no longer significant. The apparent
widening is mostly a pooling artifact. The residualized point estimates
still lean above 1, but the intervals include 1, so no robust increase is
established.

## Per-cell-type supplement (Yazar only, exploratory, uncorrected)

Same OLD vs YOUNG variance comparison inside single cell types, for the 9
types with at least 20 donors with a distance in each bin. Raw distances,
not pool-residualized, so the pool caveat applies.

- 8 of 9 types have a variance ratio above 1 (range 0.94 to 1.66).
- 3 of 9 pass the Levene test at p below 0.05: central memory CD4 T
  (ratio 1.41, p 7.8e-07), naive CD4 T (ratio 1.53, p 2.4e-05), natural
  killer cells (ratio 1.19, p 0.0034).
- Naive CD8 T, the exception in the earlier spread analyses, is flat here
  (ratio 0.94, p 0.98).
- No Gao cell type met the 20-donor-per-bin floor, so Gao has no supplement.

This supplement is descriptive only. It was not pool-adjusted and carries no
multiple-testing correction. Full table:
`donor-spread-variance-by-cell-type.csv`.

## Caveats

1. This measures spread of distances, not spread of full profiles. A donor
   cloud can widen in profile space while distances to one reference stay
   equally spread, and vice versa.
2. YOUNG donors' distances use a leave-one-out young reference (each young
   donor is excluded from its own reference), while MID and OLD donors use a
   fixed young reference. The YOUNG arm's spread is therefore not perfectly
   comparable. The OLD vs MID comparisons avoid this problem, and they agree
   with the verdicts above.
3. Pooled distances mix cell types, so composition differences between donors
   contribute to spread. The per-cell-type supplement partly addresses this
   for Yazar.
4. Gao YOUNG has only 20 donors, so the Gao OLD vs YOUNG comparison is
   imprecise.
5. The Yazar OLD vs MID raw comparison shows a Levene/bootstrap disagreement
   (significant Levene p, interval including 1). Both numbers are reported;
   neither alone settles it, and the effect size is small either way.

## Overall verdict

- Gao heart: between-donor spread in distance-to-young-reference is FLAT
  with age. No evidence of an increase.
- Yazar blood: between-donor spread looks increased on the raw distance but
  the increase DOES NOT SURVIVE pool adjustment. No robust increase is
  established.

Neither dataset proves that old donors are more variable than young donors
in distance-to-young-reference once batch structure is accounted for.

## Output files

- `donor-spread-variance-results.csv`: all comparisons with variances,
  ratios, intervals, test statistics, p-values, N per group, method.
- `donor-spread-variance-descriptives.csv`: per-bin N, mean, SD, quartiles,
  IQR for each distance version.
- `donor-spread-variance-by-cell-type.csv`: per-cell-type supplement.
- `run_donor_spread_variance.py`: the analysis script, for reproduction.
