# Held-out validation of the cleaned mature-to-old signature (P1m)

## Question

Does P1m (1,914 genes, 18 vs 30 months, kidney/lung/marrow) capture real age
signal, or is it noise fitted to 8 mice? The cross-species transport test
mostly fails, and without a positive control that failure is uninterpretable:
it could mean no conserved aging program, or it could mean the probe is bad.

## Method (frozen rules throughout)

- Derivation: pseudobulk counts per mouse per cell type, design sex plus age
  group, PyDESeq2, contrast 30m vs 18m, eligibility at least 3 mice and at
  least 50 cells on both sides, FDR below 0.05 and absolute log2 fold change
  at least 0.25, same clean-gene filtering and floor-10 rule as
  `scripts/run_mature_p1m_derivation.py` (whose functions were reused
  directly). 21m and 24m animals never enter any derivation.
- Score: per donor per cell type, mean log1p counts-per-million of up genes
  minus mean log1p counts-per-million of down genes, missing genes omitted.
  This is the same score the transport test uses.
- Composite test: ordinary least squares of score on age group plus sex
  (sex included only when it varies); two-sided p; Benjamini-Hochberg
  within the tested family.
- Empirical test: 2,000-draw matched-random null. Random sets match the
  probe on up/down counts and expression-decile distribution, probe genes
  excluded. Empirical p uses the centered rule from
  `src/followup_muscle.py` (two-sided on centered effects); Benjamini-Hochberg
  within the family. Seed 20261122.
- The random draws were generated with a vectorized permutation (argsort of
  uniforms), the same sampling distribution as choice-without-replacement.
  Validated against the verbatim loop from
  `scripts/run_rev2_transport_null.py` on one real kidney unit: per-draw-mean
  Kolmogorov-Smirnov D 0.070, p 0.713; identical empirical p 0.0050 from all
  three nulls. Two independent vectorized runs differ more from each other
  (D 0.155) than vectorized-vs-verbatim, so the residual difference is Monte
  Carlo noise.
- Implementation: `scripts/run_mature_heldout_validation.py` (variants 1-3)
  and `scripts/run_mature_heldout_validation_3b.py` (variant 3b supplement).
  Main run 529.9 seconds on beast. Model: Muse Spark (subagent).

## Variant 1: cross-organ 3-fold (main variant)

Derive a fresh signature from two organs, test it in the third organ it
never saw. Pass bar per unit: both tests q below 0.05 in the predicted
(positive) direction.

| Fold | Derive from | Held out | Signature genes | Units tested | Both-pass (predicted) | Empirical pass |
|---|---|---|---|---|---|---|
| 1 | Lung + Marrow | Kidney | 1,858 (1,109 up, 749 down) | 7 of 11 | 0 | 7 |
| 2 | Kidney + Marrow | Lung | 1,906 (1,115 up, 791 down) | 8 of 12 | 0 | 8 |
| 3 | Kidney + Lung | Marrow | 111 (36 up, 75 down) | 11 of 13 | 0 | 7 |
| Total | | | | 26 of 36 | 0 | 22 |

- All 26 tested slopes are positive (range 0.059 to 0.495, median 0.160).
  Every held-out organ moves in the predicted direction.
- 22 of 26 pass the empirical test, all 22 with positive slopes. The 4
  empirical failures are all in fold 3: granulocyte (empirical p 0.138),
  granulocytopoietic cell (0.132), proerythroblast (0.320), precursor B cell
  (0.102).
- 0 of 26 pass the composite test. Smallest composite q is 0.0595 (lung
  bronchial smooth muscle, composite p 0.00744). Six units have composite p
  below 0.05 uncorrected; median composite p is 0.0996. The composite test
  has almost no power here: 6 to 8 mice, 3 to 5 residual degrees of freedom.
- Scale check: median observed slope 0.160 against median random-slope
  standard deviation 0.016, a 10-fold gap. That is why the empirical p hits
  its floor (0.0005, i.e. 1 of 2,001) for 22 units while the composite test
  fails: the empirical test asks whether the gene set carries more age signal
  than a random set of the same size, and it does by an order of magnitude.
- 10 units untestable, all for the same reason: fewer than 3 mice per side
  after the 20-cell-per-donor floor (kidney T cell, lymphocyte, cortex artery
  cell, B cell; lung mature NK T cell, CD4 T cell, T cell, vein endothelial
  cell; marrow late pro-B cell, naive B cell).
- BH recomputation from the raw p-values matches the stored q-values to a
  maximum absolute difference of 4.4e-16.

## The signature is 94 percent marrow (found while running the folds)

- 1,848 of 1,914 P1m genes have at least one marrow supporting unit; 1,803
  are supported ONLY by marrow. Only 88 genes have any kidney support and 26
  have any lung support.
- Per-unit DE counts in the original 4v4 derivation: marrow
  granulocytopoietic cell 1,577; granulocyte 497; megakaryocyte-erythroid
  progenitor 344; naive T cell 324; proerythroblast 300; promonocyte 267;
  hematopoietic precursor 252; erythroblast 240; macrophage 77; monocyte 36.
  Best kidney unit 36, best lung unit 9.
- Consequence: P1m is a marrow-dominated signature with small kidney and
  lung contributions, not a three-organ signature. Fold 3's 111-gene
  signature is what kidney plus lung alone produce. Folds 1 and 2 test a
  marrow-heavy signature in kidney and lung; fold 3 tests a kidney/lung
  signature in marrow.

## Variant 2: intermediate ages (existing P1m, no new derivation)

Per-mouse mean P1m scores across eligible units (14 mice, 334 mouse-unit
rows):

| Age | Mice | Mean score |
|---|---|---|
| 18 months | 4 | -0.725 |
| 21 months | 2 | -0.718 |
| 24 months | 4 | -0.593 |
| 30 months | 4 | -0.467 |

Spearman rho of per-mouse mean score vs age in months is 0.733, p 0.00285.
The 18-to-21 step is flat (-0.725 to -0.718); the rise appears at 24 and 30
months. Per-mouse unit counts range from 9 to 36 (mice lack cells in some
units; e.g. 18-M-53 has no kidney cells). Full per-mouse table:
`results/followup-mature-heldout-intermediate-ages-per-mouse.csv`.

## Variant 3: donor split (descriptive only, as specified)

Hold out 18-M-53 and 30-M-5, derive on 3 mice per side.

- Signature: 84 genes (22 up, 62 down). Only 20 units eligible, all
  kidney/lung.
- Direction check: 23 of 23 held-out unit checks match the predicted
  direction.
- Limitation that defines this variant: all 13 marrow units drop out,
  because every eligible marrow unit has only 3 mice at 30m (30-M-2, 30-M-3,
  30-M-5), and holding out 30-M-5 leaves 2, below the frozen 3-mouse floor.
  The specified split therefore validates only the kidney/lung subset (84
  genes) and says nothing about the marrow 94 percent. It does not show the
  signature is donor-driven; it shows the specified holdout cannot test
  marrow.

## Variant 3b: donor split that keeps marrow eligible (supplement)

Hold out 18-M-53 and 30-M-4 instead. 30-M-4 has no marrow cells, so every
tissue keeps exactly 3 mice per side and the frozen >=3 rule still holds.
Same frozen derivation rules, same direction check.

- Signature: 687 genes (565 up, 122 down); 609 of them are in P1m.
  31 units eligible.
- Marrow DE survives at 3v3: granulocytopoietic 379 (vs 1,577 at 4v4),
  hematopoietic precursor 287 (vs 252), promonocyte 122 (vs 267), naive T
  114 (vs 324), granulocyte 92 (vs 497). The marrow signal drops 3- to
  4-fold but does not collapse, so it is not driven by a single mouse.
- Direction check: 32 of 32 match (kidney 10 of 10, lung 15 of 15, marrow
  7 of 7). Marrow checks are 7 because only 18-M-53 has marrow cells among
  the held-out mice.
- Caveat: PyDESeq2 warns at n=6 that residual degrees of freedom are below
  3 and the dispersion trend fit falls back to a mean-based trend in some
  units. Dispersion estimates are noisy at 3v3; the 687-gene count is a
  noisy estimate, but the direction check does not depend on DE precision.

## Variant 4: FACS arm — NOT RUN

- No Tabula Muris Senis FACS files exist locally or on beast. The only FACS
  artifact is `data/manifest/tms_facs_brain_nonmyeloid.json` (brain only).
- Download size checked first via the CZ CELLxGENE API (Tabula Muris Senis
  collection 0b9d8a04-bb9d-44da-aa27-705bb65b54eb): kidney Smart-seq2 33,234,927
  bytes (1,833 cells), lung Smart-seq2 108,076,230 bytes (5,218 cells), bone
  marrow Smart-seq2 339,648,372 bytes (14,517 cells). Total 480,959,529 bytes
  (about 481 MB). Beast has 38 GB free, so disk is not the blocker.
- Decisive: the FACS datasets have ages 3m, 18m, and 24m only. Donors at 18m
  are 18_45_M, 18_47_F, 18_46_F, 18_53_M; at 24m are 24_58_M, 24_59_M,
  24_60_M, 24_61_M. There are NO 30-month mice in FACS, so the specified
  18-vs-30 FACS test is impossible even after download. The only FACS
  contrast available would be 18m vs 24m, a different endpoint. No download
  was made.

## Verdict

Against the pre-specified strict bar (both tests q below 0.05 in the
predicted direction): 0 of 26. The bar is NOT met.

That binary outcome is a power failure of the composite test, not a signal
failure, and the numbers show why:

1. All 26 held-out slopes are positive, median 0.160, ten times the median
   random-null slope SD (0.016).
2. 22 of 26 pass the calibrated empirical test, all in the predicted
   direction. The paper itself treats the empirical q as the calibrated
   evidence number and the composite q as an uncalibrated screening
   statistic.
3. The dose-response is monotonic (Spearman 0.733, p 0.00285) across the
   14 mice including the held-out 21m and 24m animals.
4. Donor direction checks pass 23 of 23 (specified split, kidney/lung
   subset) and 32 of 32 (3b supplement, including 7 of 7 marrow).

Conclusion: the probe captures real mouse age signal. The cross-species
transport failure is the failure of a real, marrow-dominated mouse signal to
transport, not the failure of noise. The top-tier reviewer's positive-control
requirement is met in substance on the calibrated test, but NOT in the
letter of the both-tests bar, because the composite test cannot pass with 6
to 8 mice (smallest q 0.0595).

Two consequences for the manuscript (not done here; the tex was not edited):

1. "Three-organ signature" is wrong. P1m is 94 percent marrow-only
   (1,803 of 1,914 genes). The paper should call it marrow-dominated.
2. The disclosure "No leave-one-donor-out check was run on the cleaned
   signature" is now outdated. A donor check has been run (variants 3 and
   3b); the 310-of-765 leave-one-donor-out figure still belongs to P1, not
   P1m, and must not be attached to P1m.

## Files

- `scripts/run_mature_heldout_validation.py` (variants 1-3)
- `scripts/run_mature_heldout_validation_3b.py` (variant 3b)
- `results/followup-mature-heldout-cross-organ-transport.csv` (26 tested,
  10 untestable, all three folds)
- `results/followup-mature-heldout-fold1-leave-Kidney-out-genes.csv` (1,858),
  `...-transport.csv`, `...-scores.csv`, `...-unit-ranking.csv`
- `results/followup-mature-heldout-fold2-leave-Lung-out-genes.csv` (1,906),
  `...-transport.csv`, `...-scores.csv`, `...-unit-ranking.csv`
- `results/followup-mature-heldout-fold3-leave-Marrow-out-genes.csv` (111),
  `...-transport.csv`, `...-scores.csv`, `...-unit-ranking.csv`
- `results/followup-mature-heldout-intermediate-ages-per-unit.csv` (334 rows),
  `results/followup-mature-heldout-intermediate-ages-per-mouse.csv` (14 rows)
- `results/followup-mature-heldout-donor-split-genes.csv` (84),
  `results/followup-mature-heldout-donor-split.csv` (23 rows),
  `results/followup-mature-heldout-donor-split-unit-ranking.csv`
- `results/followup-mature-heldout-donor-split-3b-genes.csv` (687),
  `results/followup-mature-heldout-donor-split-3b.csv` (32 rows),
  `results/followup-mature-heldout-donor-split-3b-unit-ranking.csv`
- `results/followup-mature-heldout-summary.json`
- Expression-matrix caches (not committed): `scratch/heldout/expr-{Kidney,Lung,Marrow}.npz`
