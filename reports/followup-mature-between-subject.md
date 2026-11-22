# Item 9b: mature between-subject donor spread at 30-plus

## Design

Between-subject spread asks whether donors become more different from each
other with age. It is the donor-level version of the randomness question:
if aging is random drift, old donors should sit farther apart in expression
space than mature donors do, and gene-level variance across donors should
rise. This rerun restricts everything to donors aged 30 and older, with
the reference group mature donors aged 30 to 39, the comparison group old
donors aged 60 and older, and mid donors aged 40 to 59 as a secondary
contrast. Donor is the unit in every test. No young-versus-mature
comparison is used anywhere.

Two parts. Part A is donor-distance spread: each donor gets a Hellinger
distance from the leave-one-out mature reference profile, then the spread
of those distances is compared between age groups with a Levene test
(median-centered, also called Brown-Forsythe) on the variance ratio, with
a 2000-resample bootstrap percentile interval for the ratio. Part B is
gene-level variance: for each cell type, the mean log variance ratio and
the total log variance ratio of old versus mature donors across genes,
plus a permutation eta2 for the age-group mean effect, each with a
Benjamini-Hochberg correction over the tested cell types in that contrast
and arm. Two residualization arms run in each part: raw, and batch-mean
residualized (in Yazar the batches are the pools, so this arm is the same
computation Item 8 called pool residualized).

## Cohorts

- Gao heart: 109 donors, age 30 to 91. Mature 30 to 39: 12 donors. Mid
  40 to 59: 53 donors. Old 60-plus: 44 donors. A1 distance fits use 104
  donors, A3 dispersion fits use 105 donors; the rest have no eligible
  cells at the 1500 molecule floor.
- Yazar blood: 932 donors, age 30 to 97. Mature 30 to 39: 59 donors. Mid
  40 to 59: 193 donors. Old 60-plus: 680 donors. A1 and A3 fits use 925
  donors.

Source files: `results/followup-mature-between-subject-gao-donors.csv`,
`results/followup-mature-between-subject-yazar-donors.csv`,
`results/followup-mature-between-subject-gao-models.csv` (508 rows: A1
Hellinger distance fits and A3 within-donor dispersion fits),
`results/followup-mature-between-subject-yazar-models.csv` (350 rows),
`results/followup-mature-between-subject-distance-variance-gao.csv` (112
rows), `results/followup-mature-between-subject-distance-variance-yazar.csv`
(232 rows),
`results/followup-mature-between-subject-variance-gao.csv` (156 rows),
`results/followup-mature-between-subject-variance-yazar.csv` (220 rows),
`results/followup-mature-between-subject-variance-gao-ineligible.csv`,
`results/followup-mature-between-subject-variance-yazar-ineligible.csv`,
`results/followup-mature-between-subject-gao-reference-audit.csv`,
`results/followup-mature-between-subject-yazar-reference-audit.csv`,
`results/followup-mature-between-subject-runtime.csv`.

## Part A: donor-distance spread

### A1 and A3 slopes first (stratified global, full depth, 1500 floor)

The A1 slope is the change in donor distance from the mature reference
per decade of age. The A3 slope is the change in within-donor dispersion
per decade. Both are negative in Yazar and null in Gao: donors move
closer to the mature reference with age, and cells become less dispersed
within a donor. That is the opposite of the widening prediction, and it
matches the Item 9a full-depth finding that spread falls with age.

- Gao A1 full_1500: slope -0.0075 SD per decade, p 0.767, n 104 donors.
- Gao A3 full_1500: slope -0.0619 SD per decade, p 0.0846, n 105 donors.
- Yazar A1 full_1500: slope -0.0261 SD per decade, p 3.55e-06, n 925
  donors. At 1500 molecules the slope is -0.0819, p 2.07e-10; at 500 it
  is -0.119, p 2.49e-12.
- Yazar A3 full_1500: slope -0.0776 SD per decade, p 1.97e-13, n 925
  donors. At 1500 molecules the slope is -0.1359, p 1.64e-14.

### Levene variance ratios, pooled donors, mature 30 to 39 reference

Ratio is old-group variance divided by mature-group variance. Ratio
above 1 means old donors are more spread out.

Gao, OLD 60-plus vs MATURE 30 to 39, 12 vs 41 donors:

- Raw: ratio 1.048, Levene p 0.623, bootstrap 95 percent interval 0.388
  to 6.457.
- Batch residualized: ratio 0.894, Levene p 0.824, interval 0.482 to
  2.426.

Gao, MID 40 to 59 vs MATURE 30 to 39, 12 vs 52 donors: raw ratio 0.611,
p 0.706; batch ratio 1.111, p 0.903. Gao, OLD vs MID, 52 vs 41 donors:
raw ratio 1.716, p 0.149; batch ratio 0.805, p 0.923. Every Gao pooled
row is null, and the OLD vs MATURE interval is very wide because the
mature reference has only 12 donors.

Yazar, OLD 60-plus vs MATURE 30 to 39, 680 vs 59 donors:

- Raw: ratio 1.543, Levene p 3.46e-05, bootstrap interval 1.075 to
  2.624. The raw widening survives the mature restriction.
- Batch (pool) residualized: ratio 1.578, Levene p 0.137, interval
  1.062 to 2.582. The Levene test is null after residualization while
  the bootstrap interval still excludes 1. The two methods disagree, so
  this arm is NOT claimed as a widening signal.

Yazar, MID 40 to 59 vs MATURE 30 to 39, 193 vs 59 donors: raw ratio
1.387, p 0.127; batch ratio 1.218, p 0.706. Yazar, OLD vs MID, 680 vs
193 donors: raw ratio 1.112, p 0.000823; batch ratio 1.296, p 0.0597.

### Per-cell-type distance variance, Yazar OLD vs MATURE

Only 6 of 19 Yazar cell types are testable for OLD vs MATURE at the
mature reference (the rest have fewer than 20 donors per group or no
eligible reference cells). Raw arm: central memory CD4 T cells ratio
1.348, p 2.7e-05, q 0.000164; naive CD4 T cells ratio 1.295, p 0.0123,
q 0.03695. Both die in the batch arm (central memory CD4 q 0.588,
naive CD4 q 0.818). No cell type passes in the batch arm. For OLD vs
MID, only central memory CD4 passes in the raw arm (ratio 1.063,
p 0.000838, q 0.0109) and nothing passes in the batch arm. Gao has no
testable per-cell-type row for OLD vs MATURE or MID vs MATURE (only the
pooled row is tested); its OLD vs MID per-cell-type rows are all null,
smallest q 0.411 in the batch arm.

### Retained 40 to 59 reference sensitivity

Keeping the original 40 to 59 reference instead of the mature 30 to 39
reference, OLD vs MID pooled: Yazar raw ratio 1.128, p 0.000472; batch
ratio 1.287, p 0.0652. Gao raw ratio 1.679, p 0.128; batch ratio 0.879,
p 0.945. The Yazar raw signal is present under either reference and the
batch arm is null under either reference.

## Part B: gene-level variance, OLD vs MATURE 30 to 39

Mean log var ratio is the average across genes of log old variance over
mature variance. Total log var ratio is the log ratio of summed
variances. Both are tested one-sided for widening (old greater than
mature) with BH correction over the cell types in each contrast and arm.
Each cell type appears twice, at pseudocount eps 1e-08 and 0.001; the
two rows share p and q values, so counts below are per cell type.

Gao: 10 cell types tested per arm. Zero significant in any of the four
q columns (mean log var ratio one-sided and two-sided, total log var
ratio one-sided, eta2) in any contrast. Global rows are null as well:
stratified global mean log var ratio +0.042 batch arm, p 0.583
one-sided, and +0.055 depth-only arm, p 0.568; total log var ratio
-0.058 batch arm, p 0.744, and +0.049 depth-only arm, p 0.319.

Yazar: 16 cell types tested per arm. Only two widen, and only in the
gene-level test:

- Naive CD4 T cells: mean log var ratio +0.278 batch arm, p 0.00350,
  q 0.02799; +0.318 depth-only arm, p 0.0009995, q 0.007996. Total log
  var ratio +0.269 in both arms, q 0.01599 batch and 0.007996
  depth-only. Reference 59 donors, comparison 675 donors.
- Naive CD8 T cells: mean log var ratio +0.456 batch arm, p 0.0004998,
  q 0.007996; +0.622 depth-only arm, p 0.0004998, q 0.007996. Total log
  var ratio +0.455 batch arm and +0.587 depth-only arm, q 0.007996 in
  both. Reference 59 donors, comparison 391 to 392 donors.

No other Yazar cell type widens. The stratified global goes the other
way: mean log var ratio -0.103 batch arm and -0.068 depth-only arm,
one-sided p 1.0 and 0.999. The pooled All_Cells total log var ratio is
+0.141 batch arm, p 0.0270, and +0.099 depth-only arm, p 0.0670, while
its mean log var ratio is null (+0.014, p 0.532 batch; +0.035,
p 0.405 depth-only), so the pooled total is carried by the two naive T
types, not by a global widening.

The eta2 permutation test is significant in 23 Yazar rows covering 12
cell types, but every one of them has a NEGATIVE mean log var ratio
(for example natural killer cells -0.369 batch arm, q eta2 0.001142;
regulatory T cells -0.158, q eta2 0.001142; transitional B cells -0.180,
q eta2 0.001142). Eta2 here measures age-group mean shifts, not spread:
old donors differ in mean expression and have LOWER variance. It does
not support widening.

For the other contrasts in Yazar, naive CD8 T cells also widen in OLD
vs MID (both arms, mean log var ratio q below 0.05) and in MID vs
MATURE (both arms), and naive CD4 T cells widen in OLD vs MID on the
total log var ratio. Gao has zero significant rows in OLD vs MID and
MID vs MATURE as well.

## Verdicts

- Gao: DOES NOT PROVE. Null on every arm and every contrast, with only
  12 mature reference donors, so the OLD vs MATURE interval is wide
  (0.388 to 6.457 raw) and per-cell-type tests are untestable.
- Yazar: DOES NOT PROVE for the headline claim that donors spread apart
  with mature age. The raw pooled widening (ratio 1.543, p 3.46e-05)
  does not survive batch and pool residualization (Levene p 0.137), and
  the two raw per-cell-type passes die in the batch arm. The A1 and A3
  slopes point the other way: donor distances shrink and within-donor
  dispersion falls with age.
- Secondary finding, Yazar only: naive CD4 and naive CD8 T cells show
  gene-level variance widening in both arms at 30-plus. This is a
  per-cell-type signal, not a global one: the stratified global is
  negative and no other cell type widens.

## Comparison with the November 17 original

Original, all ages with a 40 to 59 reference: Yazar raw widening ratio
1.62, Levene p 1.9e-07, which vanished after pool residualization to
ratio 1.25, p 0.43; Gao null at ratio 0.79, p 0.85. Mature restriction
with a 30 to 39 reference: Yazar raw ratio 1.543, p 3.46e-05, batch
residualized ratio 1.578, Levene p 0.137; Gao raw ratio 1.048,
p 0.623. The shape is unchanged: a raw Yazar signal that does not
survive residualization, and a Gao null. The mature version is stricter
in one way (the batch arm ratio stays above 1 at 1.578 instead of
falling to 1.25) and weaker in another (the batch Levene p of 0.137 is
not significant, and the bootstrap interval is the only positive
indicator in that arm).

Two bookkeeping notes on the numbers. First, the bootstrap intervals in
the 9b distance files differ slightly from the committed Item 8 file
`results/followup-mature-nov17-donor-spread.csv` for the same pooled
rows: Yazar OLD vs MATURE raw interval 1.075 to 2.624 here versus 1.087
to 2.687 there, and batch interval 1.062 to 2.582 here versus 1.075 to
2.562 there. Point estimates, Levene statistics, and p values are
identical. Both use 2000 resamples at seed 7; the small difference is
bootstrap resampling detail between the two implementations, and no
conclusion changes. Second, the 9b batch arm and the Item 8 pool arm
are the same computation in Yazar (pools are the batches): residualized
standard deviations 0.004185 and 0.003331, ratio 1.578469, Levene
p 0.136603 in both files.

## Limits

- Gao mature reference is 12 donors. Per-cell-type OLD vs MATURE and
  MID vs MATURE distance tests are untestable, and the pooled interval
  is wide. Any Gao widening below a ratio of about 2 would be missed.
- Yazar OLD vs MATURE per-cell-type distance tests cover only 6 of 19
  types. The rest have fewer than 20 donors per group or no eligible
  reference cells; see the ineligible lists.
- The batch residualized arm and the bootstrap interval disagree in
  Yazar OLD vs MATURE (Levene p 0.137, interval 1.062 to 2.582). The
  interval is not a significance test on its own, and the Levene result
  governs the verdict.
- The naive T cell gene-level widening is a variance result, not a
  transport or drift result. It does not show that the signature
  transports or that donors drift apart.
- The eta2 permutation p values bottom out at the 2000-draw floor
  (0.0004998). Small eta2 values (mean 0.0026 for naive CD4) are
  statistically significant against the permutation null but small in
  absolute terms.

## Runtime and provenance

- Gao main pipeline failed with KeyError `residualized_distance` in the
  distance step (cell types with zero eligible distances produced an
  empty frame) and was fixed by a rerun that initializes the column
  before residualization. Distance rerun: 112 rows, 40 tested, 72
  untestable, 1.9 seconds.
- Yazar expression A2 ran in parallel on 6 workers: 87 of 87 units,
  3616 seconds, 220 variance rows, 35 ineligible rows. Yazar distance
  A2: 232 rows, 82 tested, 150 untestable, 4.0 seconds.
- The `model_used` value in the runtime file for the Gao main run reads
  Gemini 4 Neon XHigh. That value is hardcoded in the driver script;
  the reruns were executed directly with no delegated model, and their
  runtime rows say so.

## Verification

- Independent BH recomputation: every q value in the two distance files
  recomputed from the p values with the correct families (pooled rows
  are families of 1; cell-type rows group by version, adjustment, and
  contrast). Zero mismatches. Every q value in the two gene-level
  variance files recomputed with the correct families (one entry per cell
  type per contrast and arm; the two eps rows share p and q). Zero
  mismatches across all 48 families in each file.
- Headline numbers in this report were recomputed from the CSVs, not
  copied from the run logs. A1 and A3 stratified-global slopes were
  checked against `results/followup-mature-between-subject-gao-models.csv`
  and `results/followup-mature-between-subject-yazar-models.csv`.
- `tests/test_mature_item9b.py` passes on beast in the agescape env. It
  checks the mature donor selection, the contrast definitions, the Part A
  output contract, and the PBMC P1m contract.
- The two `.pkl` intermediate files
  (`results/followup-mature-between-subject-yazar_a2_parallel.pkl`,
  `results/followup-mature-between-subject-yazar_a2_test_nk.pkl`) are
  not committed: they are parallel-run intermediates, reproducible from
  the scripts, and the assembled CSVs carry the results.
