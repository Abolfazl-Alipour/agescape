# Mature M10 PsychAD pipeline at 30-plus with P1m/P0m

## Scope: what was already done, what this adds

Phase 2 Item 2 (commit 3f5b57a) already reran the core PsychAD transport with
P1m/P0m at 30-plus: the fully adjusted per-cell-type slopes, the 2000-draw
expression-matched null, and the strict P0m rescue comparison. Those results
are reused here and never recomputed as new claims. The fully adjusted slopes
in this run reproduce the committed Item 2 transport slopes exactly
(maximum absolute difference 0, verified from the CSVs).

What Item 2 did not produce, and what this run supplies, is the rest of the
M10 pipeline itself, mirroring the original `scripts/run_m10_pipeline.py`
outputs with the mature probe sets and donors aged at least 30:

1. Per-cell-type unadjusted slopes and raw-count sensitivity slopes
   (`m10-celltype-regressions` equivalent).
2. Stratified fits at 30-plus: age at least 40, HBCC-only, MSSM-only
   (`m10-confounder-audit` equivalent). The adult-at-least-20 and
   developmental-under-20 strata are not applicable: every fitted donor is
   aged at least 30.
3. Whole-cohort pooled summary at 30-plus (`m10-cohort-summary` equivalent).
4. A 200-draw matched null as HISTORICAL SENSITIVITY ONLY, same seed and
   procedure as Item 2's 2000-draw null. The 2000-draw null stays primary.

Block 5 PsychAD rows were already run (commit 9eb7652): the 283-row mature
Block 5 comparison contains all 27 DLPFC_Lifespan_Atlas units. No PsychAD
arm from the Phase 1 audit remains unrun after this commit.

- Input: `data/external/dlpfc-atlas/lifespan_pseudobulk.rds`, 539,406,079
  bytes, MD5 `aef628e7bc2e5134b790d98c6f57ec48`, verified before the run.
  The controlled DUA file stays on beast; it is never copied to git.
- Donors: 201 at 30-plus (of 284 controls aged 0.2 to 97), age 30 to 97.
- Probe sets: human-mapped P1m (1,724 genes: 1,013 up, 711 down) and P0m
  (2,049 genes: 1,157 up, 892 down) from Item 1.
- Models: OLS per cell type on log1p(CPM) pseudobulk scores, up-gene mean
  minus down-gene mean. Fully adjusted model: intercept, age, male sex,
  HBCC brain-bank indicator, PMI, TechPC1, TechPC2, TechPC3, exactly as in
  M10. BH correction is applied separately per tier over the 27 PsychAD
  units, unless stated otherwise.
- Runtime: 5.7 minutes total (P0m 2.2, P1m 3.0). Model used:
  direct-runner, no delegated model.

## Headline verdicts

1. **Pooled PsychAD score rises with age at 30-plus.** Mean score across
   all 27 cell types, 201 donors: P1m unadjusted slope +0.000508 per year
   (p 0.000547), covariate-adjusted +0.000371 per year (p 0.0241). P0m
   unadjusted +0.000486 per year (p 0.00182), adjusted +0.000385 per year
   (p 0.0288). The original full-range pooled P1 slope was null (+0.000186
   per year, p 0.349, 284 donors); the mature pooled rise is new.
2. **11 of 27 units pass the fully adjusted test per tier, the same 11 in
   both tiers.** Ten have positive slopes (seven inhibitory interneuron
   types, three excitatory types); one, Micro, has a negative slope.
3. **Micro flips direction between full range and mature.** Original
   full-range P1 Micro slope +0.002381 per year (p 4.45e-06, q 0.000106);
   mature P1m slope -0.001786 per year unadjusted (p 1.61e-08,
   q 2.17e-07) and -0.001433 adjusted (p 8.06e-05, q 0.000725). The
   original positive Micro transport was young-driven. The mature negative
   Micro slope is significant on both the log1p-CPM and raw-count scales.
4. **The 10 positive mature slopes are log1p-CPM-only.** Nine of the ten
   have negative raw-count slopes (all null except IN_SST); only EN_L6_CT
   keeps a positive raw slope (null, p 0.63). IN_SST is a direct sign
   conflict: CPM +0.001201 adjusted (p 1.03e-05, q 0.000139) against raw
   -0.001704 (p 0.00139, q 0.00940). The raw scale is confounded by
   library size, so this does not refute the CPM results, but only Micro's
   negative slope is strengthened by agreement across both scales.
5. **The 200-draw null kills every empirical pass.** At 200 draws, zero
   PsychAD units pass the empirical test in either tier. At 2000 draws,
   four pass (P0m 3, P1m 1). All four 200-draw empirical p-values sit at
   or near the 1/201 floor (0.004975), with BH q 0.0597 and 0.0672, just
   above 0.05. The 2000-draw null is essential; the 200-draw run cannot
   resolve these passes.

## Pooled cohort summary (mature vs original full-range)

| Probe set | Donors | Age range | Unadjusted slope/yr | p | Adjusted slope/yr | p (adjusted) |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
| P1m mature | 201 | 30 to 97 | +0.000508 | 0.000547 | +0.000371 | 0.0241 |
| P1 original | 284 | 0.2 to 97 | +0.000186 | 0.349 | -0.000209 | 0.400 |
| P0m mature | 201 | 30 to 97 | +0.000486 | 0.00182 | +0.000385 | 0.0288 |
| P0 original | 284 | 0.2 to 97 | +0.000318 | 0.0118 | +0.000027 | 0.863 |

Slopes are score units per year. The pooled score is the mean donor score
across all 27 cell types. The mature pooled rise is significant on both
scales of the test (unadjusted and adjusted) for both probe sets; the
original full-range pooled P1 rise was null on both.

## Per-cell-type regressions: the 11 fully adjusted passes

BH q over the 27 PsychAD units per tier. The 11 passing units are identical
in both tiers. Slopes are score units per year.

| Unit | P1m adjusted slope | P1m adjusted p | P1m adjusted q | P0m adjusted slope | P0m adjusted q |
| --- | ---: | ---: | ---: | ---: | ---: |
| IN_ADARB2 | +0.001612 | 1.85e-06 | 4.99e-05 | +0.001574 | 3.71e-05 |
| IN_SST | +0.001201 | 1.03e-05 | 0.000139 | +0.001226 | 0.000166 |
| Micro | -0.001433 | 8.06e-05 | 0.000725 | -0.001523 | 0.000379 |
| EN_L6_CT | +0.001166 | 0.000503 | 0.00340 | +0.001161 | 0.00388 |
| IN_PVALB | +0.000738 | 0.00148 | 0.00798 | +0.000861 | 0.00335 |
| IN_LAMP5_LHX6 | +0.000960 | 0.00259 | 0.0116 | +0.000999 | 0.00850 |
| IN_PVALB_CHC | +0.000889 | 0.00365 | 0.0136 | +0.000969 | 0.0108 |
| IN_VIP | +0.000685 | 0.00403 | 0.0136 | +0.000764 | 0.00850 |
| IN_LAMP5_RELN | +0.000957 | 0.00538 | 0.0162 | +0.000980 | 0.0108 |
| EN_L3_5_IT_2 | +0.000716 | 0.0201 | 0.0498 | +0.000709 | 0.0457 |
| EN_L5_6_NP | +0.000876 | 0.0203 | 0.0498 | +0.001070 | 0.0116 |

Counts across all 27 units per tier (BH over 27, identical sets in both
tiers): unadjusted q below 0.05 in 18 units (original full-range: P1 10 of
27, P0 14 of 27); fully adjusted q below 0.05 in the 11 units above;
raw-count q below 0.05 in 5 units (Endo, IN_SST, Micro, PC, SMC). Fully
adjusted slope signs: 21 positive, 6 negative (Micro, SMC, PVM, VLMC,
Endo, Adaptive). The six negative units other than Micro are all null on
the adjusted test.

## Micro: the direction flip, in full

| Version | Scale | Slope/yr | p | q (27-unit family) |
| --- | --- | ---: | ---: | ---: |
| Original full-range P1 | log1p-CPM unadjusted | +0.002381 | 4.45e-06 | 0.000106 |
| Original full-range P1 | log1p-CPM adjusted | +0.002069 | 0.00150 | not computed |
| Original full-range P1 | raw counts | -0.000841 | 0.257 | not computed |
| Mature P1m | log1p-CPM unadjusted | -0.001786 | 1.61e-08 | 2.17e-07 |
| Mature P1m | log1p-CPM adjusted | -0.001433 | 8.06e-05 | 0.000725 |
| Mature P1m | raw counts | -0.001687 | 8.47e-09 | 2.29e-07 |
| Mature P0m | log1p-CPM adjusted | -0.001523 | 4.21e-05 | 0.000379 |
| Mature P0m | raw counts | -0.001730 | 6.99e-09 | 1.89e-07 |

Item 2's both-pass verdict for PsychAD Micro is a two-sided test: it
records that the score changes with age, not the direction. The direction
at 30-plus is negative, opposite to the original full-range positive
slope. Micro strata at 30-plus (P1m): age at least 40 slope -0.001765
(p 2.78e-05); HBCC-only -0.001912 (p 0.00153); MSSM-only -0.000711
(p 0.151, null). The negative slope holds in the older half and in HBCC
donors and is weaker in MSSM donors.

## Confounder audit: stratified agreement

For the 11 fully adjusted passes, all 11 agree in sign with the full
30-plus fit in every stratum (age at least 40, HBCC-only, MSSM-only), in
both tiers. How many of the 11 pass BH q below 0.05 within each stratum
(BH over 27 per tier per stratum):

| Tier | Age at least 40 (median n 174) | HBCC-only (median n 89) | MSSM-only (median n 107) |
| --- | ---: | ---: | ---: |
| P1m | 4 of 11 | 1 of 11 | 9 of 11 |
| P0m | 9 of 11 | 0 of 11 | 9 of 11 |

The HBCC stratum (median 89 donors) is underpowered: almost nothing
passes q below 0.05 there, but all 11 slopes agree in sign with the full
fit. MSSM-only carries most of the signal (9 of 11 pass in both tiers).
The age-at-least-40 stratum differs by probe set (P1m 4 of 11, P0m 9 of
11) while the full-fit pass sets are identical; the P1m older-stratum
slopes are the same sign but noisier.

## 200-draw null sensitivity (historical only)

Same seed (20261105) and same procedure as Item 2's 2000-draw null:
expression-decile-matched random gene sets, up and down matched
separately, pool excludes detected probe genes, fully adjusted donor
model, plus-one two-sided empirical p, BH within the 27-unit PsychAD
family per tier. One precision note: the per-tier RNG stream is consumed
sequentially across units in sorted order, so only the first unit
(Adaptive) shares its first 200 draws with the 2000-draw run. This is a
same-seed, same-procedure sensitivity, not a strict prefix of the
2000-draw draws.

| Tier | Empirical passes at 200 draws | Empirical passes at 2000 draws | Verdict changes |
| --- | ---: | ---: | --- |
| P1m | 0 of 27 | 1 of 27 | 1 PASS_2000_ONLY, 26 SAME_FAIL |
| P0m | 0 of 27 | 3 of 27 | 3 PASS_2000_ONLY, 24 SAME_FAIL |

The four units that pass at 2000 draws but not at 200:

| Unit | Tier | 200-draw p | 200-draw q | 2000-draw p | 2000-draw q |
| --- | --- | ---: | ---: | ---: | ---: |
| EN_L6_IT_2 | P0m | 0.004975 | 0.0597 | 0.002499 | 0.0337 |
| IN_SST | P0m | 0.009950 | 0.0597 | 0.003998 | 0.0360 |
| Micro | P0m | 0.004975 | 0.0597 | 0.000500 | 0.0135 |
| Micro | P1m | 0.004975 | 0.0672 | 0.000500 | 0.0135 |

The 200-draw empirical p floor is 1/201 = 0.004975. Three of the four
200-draw p-values sit exactly on that floor and the fourth is one step
above it, yet their BH q values (0.0597, 0.0672) land just above 0.05.
The 200-draw run cannot resolve any PsychAD empirical pass; the
2000-draw run can. This is the same Monte Carlo floor effect the
original Block 5 found for MSSM oligodendrocytes (200-draw q 0.0522 vs
2000-draw q 0.0105), reproduced at 30-plus.

## What this changes for the paper

- The mature PsychAD story is not a null. The pooled score rises with age
  at 30-plus, 11 units pass the fully adjusted test, and Micro passes
  both tests under P1m and P0m (Item 2). The rewrite must carry the
  directions: Micro falls with age at 30-plus; the other ten passing
  units rise.
- The Micro direction flip is the sharpest mature-to-old lesson in this
  cohort: a positive full-range transport that reverses inside the mature
  window. Any manuscript sentence that calls the PsychAD Micro result a
  positive age association is wrong at 30-plus.
- The ten positive passing slopes rest on log1p-CPM normalization only.
  The rewrite must not present them as robust to the raw-count
  sensitivity, and IN_SST must be flagged as a direct sign conflict
  between the two scales.
- The 200-draw sensitivity justifies keeping the 2000-draw null primary:
  at 200 draws the empirical family has zero passes and the floor hides
  all four real ones.

## Verification

- Every BH q in the four result files was recomputed by hand from the
  reported p-values (independent implementation, not the run script's
  `multipletests` call). Maximum absolute difference across all families:
  1.61e-15. Zero mismatches.
- Fully adjusted slopes match the committed Item 2 transport CSV exactly:
  maximum absolute difference 0 across all 54 unit-tier rows (the run
  script also asserted this in-run at rtol 1e-10, atol 1e-12).
- The `q_composite_mature170` join against the Item 2 transport CSV is
  exact: maximum absolute difference 0.
- The 200-draw empirical p floor is exactly 1/201 = 0.0049751243781094.
- Row counts: regressions 54 (27 units times 2 tiers), confounder audit
  54, cohort summary 2, null200 54.
- A test file (`tests/test_followup_mature_psychad_m10.py`) re-checks row
  counts, the hand-recomputed BH values, the transport cross-check, the
  null floor, and the verdict-change counts.

## Limits

- The raw-count sensitivity is confounded by library size (total counts
  per donor change with age), so a raw-sign flip does not refute a CPM
  result; it only withholds strengthening. IN_SST is the exception that
  must be flagged, because it is significant on both scales with opposite
  signs.
- The HBCC-only stratum (median 89 donors) is underpowered; its
  near-zero pass count is not evidence of a bank difference, since all 11
  slopes agree in sign.
- The pooled summary is computed on log1p-CPM scores only. No pooled
  raw-count summary was computed.
- The 200-draw run is a sensitivity only. No verdict in this report rests
  on it; the 2000-draw null from Item 2 stays primary.

## Outputs

- `results/followup-mature-psychad-m10-regressions.csv`: 54 rows, per-cell-type
  unadjusted, adjusted, and raw-count slopes with BH q over 27 units per tier.
- `results/followup-mature-psychad-m10-confounder-audit.csv`: 54 rows,
  stratified fits (age at least 40, HBCC-only, MSSM-only) with BH q per stratum.
- `results/followup-mature-psychad-m10-cohort-summary.csv`: 2 rows, pooled
  cohort slopes per tier.
- `results/followup-mature-psychad-m10-null200.csv`: 54 rows, 200-draw
  empirical p/q joined to the Item 2 2000-draw values with verdict-change flags.
- `results/followup-mature-psychad-m10-runtime.csv`: stage timings, MD5, seed.
- `figures/followup-mature-psychad-m10-forest-p1m.png` and
  `figures/followup-mature-psychad-m10-forest-p0m.png`: forest plots of the
  fully adjusted slopes.
- `scripts/run_mature_psychad_m10.py`: the run script.
- `tests/test_followup_mature_psychad_m10.py`: verification tests.
