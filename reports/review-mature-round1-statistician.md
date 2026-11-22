# Hostile statistical review, round 1: statistician

Manuscript: `report/plos-one/final-plos-one.tex`, commit `683a74d`
(tag `mature-to-old-rewrite-2026-11-21`), 721 lines.
Reviewer role: hostile statistician. Nothing is submitted anywhere.
Every number below was checked against the CSV in `results/` named with it,
or is marked as not rechecked. The tex was not edited.

## Scope and method

I read the tex in full and verified the load-bearing numbers against their
source CSVs by direct lookup and independent recomputation of the
Benjamini-Hochberg (BH) q values where the framing depends on them.
BH means Benjamini-Hochberg correction for multiple testing: q is the
adjusted p value, and q below 0.05 is the paper's pass threshold.
MDE80 means minimum detectable effect at 80 percent power: the smallest
true effect the test had an 80 percent chance of detecting.

The good news first: the arithmetic is mostly honest. I verified the P1m
derivation counts (1,914 genes; 1,119 up, 795 down; 1,724 human-mapped),
the P1 overlap (52 shared, 27 flips, 87 lost, 1,862 gained), the core
transport counts (170 tested, 82 untestable; P1m 7 composite passes; both
microglia passes with the exact slopes and q values), the locked 624
composition family (624 rows, 45 passes, 17 core plus 28 added; BH
recomputation matches stored q to 2e-16), both promoted tests with their
sensitivity failures, the five mature composition survivors with their
slopes and q values, the nine mouse composition deaths with SE inflation
ratios 1.30 to 5.53 (median 3.28), the mouse and human drift estimates with
their intervals, the interaction estimates, the Bayesian drift posteriors,
the Block 1 naive CD4 primary numbers, the Block 2 natural killer numbers,
the Block 3 cloud numbers, the Block 6 MSSM depth slopes and q values,
the Block 7 Part A slopes and q values and the 96 predicted pass-to-fail
count, the Block 8 zero passes, the Block 4 zero passes (0 of 472),
the Block 10 Harries (5) and Peters (16) Yazar counts and the Sood loss,
the PsychAD pooled and per-unit numbers including the IN_SST sign conflict
and the 200-draw versus 2,000-draw null counts (0 versus 4), the Thomson
PBMC CLR numbers, the GTEx Track 1 five passes and Track 2 whole blood
numbers, the RPE and Lake kidney numbers, the between-subject A1 and A3
slopes, the donor-spread ratios, the survivor-bias correlations, the depth
test slopes, the full-range drift and interaction numbers in Appendix A,
and the full-range Bayesian posteriors. The numbers check out. The
problems below are not arithmetic slips. They are contradictions between
what the text says and what its own reported q values mean, undeclared
correction families, a sensitivity variant presented as a primary result,
and significant results that are never reported.

## BLOCKER issues (must fix before submission)

### B1. Concordance is called "nominal only" in the same sentence that reports two BH q values below 0.05

Tex line 309 (Results): "Mouse within-tissue concordance passes at
nominal significance only: observed mean Pearson 0.2105 against a null
mean of 0.1465 (p = 0.0220, q = 0.0440 under the mature two-test family
and q = 0.0460 under the mature six-row mouse family), across 59 units
and 14 mice."
Tex line 333 (Discussion) and line 411 (Conclusions) repeat "nominal only".

Source files:
`results/followup-mature-concordance-mouse.csv`, row species Mouse,
scope within_tissue_or_cohort: observed 0.21047, joint null mean 0.14646,
p 0.021978, 59 units, 14 donors.
`results/followup-mature-concordance-bh-check.csv`: mature_two_test_family
mouse within q 0.043956; mature_mouse_six_rows mouse within q 0.045954.
`results/followup-mature-concordance-human.csv`: human p 0.8152.

Independent recomputation: two-test family, p values 0.0220 and 0.8152.
BH q for the smaller p is 0.0220 times 2 divided by 1, which is 0.0440.
Six-row family, p values 0.0220, 0.0230, 0.0230, 0.1029, 0.103, 0.105.
BH q for the smallest p is the cumulative minimum over ranks, which is
0.0230 times 6 divided by 3, which is 0.0460. Both reported q values are
correct. Both are below 0.05. Under either family the paper itself
declares, the mouse within-tissue result PASSES BH correction. Calling it
"nominal significance only" contradicts the paper's own numbers in the
same sentence.

The family choice decides the headline. A three-test family (mouse
within, mouse cross-tissue p 0.1029, human p 0.8152) would give
q 0.066, which fails. But no such family is computed or reported anywhere
in the repo. The two families the paper does report both pass. Two
further problems with the reported families: the six-row family counts
near-duplicate sensitivity gene sets (lenient and unfiltered; p 0.0230,
0.0230, 0.0220) as separate tests, which dilutes the correction; and the
two-test family excludes the mouse cross-tissue test (p 0.1029) run in
the same analysis without justification.

Fix: pick one primary family, report its q, and state the verdict that
follows from it. If the declared families stand, the result is
BH-significant and the "no conserved direction" framing for mouse in
Results, Discussion, and Conclusions must change. If the authors want
the nominal verdict, they must declare the three-test family as primary
and report q 0.066.

### B2. Methods contradicts itself on whether families were defined before correction

Tex line 226 (Methods): "Each analysis family was defined before
correction, and Benjamini-Hochberg correction was applied within the
stated family."
Tex line 143 (Methods): "Both the 431-test and 624-test outputs were
produced without selecting either one, and the choice of primary family
was deferred until after the mammalian p values and both correction
outcomes were available. The 624-test family is therefore the declared
primary policy for this report, not a preregistered or outcome-blind
scope."
Tex line 245 (Results) and line 397 (Discussion) repeat the post-hoc
selection admission: "The 624-test family was selected after both
correction outcomes were available."

These cannot both be true for the primary composition family. The blanket
pre-specification claim in line 226 is false for the family the whole
composition section rests on. A reviewer will find this in one read.
The honest statement is line 143. Fix: qualify line 226 to exclude the
composition family, or delete the blanket claim.

### B3. The ectopic headline numbers are a sensitivity variant, not the primary endpoint, and the primary endpoint is null

Tex line 299 (Results): "Ectopic gene counts fall with age in both
cohorts (Gao minus 0.127, p = 0.031; Yazar minus 0.0327, p = 0.0031),
again opposite the prediction."

Source: `results/followup-mature-fulldepth-gao-ectopic.csv` and
`results/followup-mature-fulldepth-yazar-ectopic.csv`,
gene_set E2_stratum, depth_mode bounded_depth_strata,
model stratified_global_primary: Gao slope -0.12687, p 0.0311; Yazar
slope -0.03269, p 0.00309. Those are the numbers in the text.

Primary endpoint, same files, gene_set E1, depth_mode full_depth,
model stratified_global_primary, analysis_role primary_full_depth:
Gao slope -0.06331, p 0.4683; Yazar slope -0.01868, p 0.0730. Both null.
The text never states that the quoted numbers come from the E2
sensitivity gene set and a bounded-depth stratum. It presents them as
the ectopic result, which inverts the primary endpoint's conclusion
(null in both cohorts) into significant falls in both cohorts.
Fix: report the primary E1 result as the headline and move the E2
numbers to a labeled sensitivity sentence.

## MAJOR issues (weaken a claim)

### M1. The mature composition q family is never declared, and the printed q values cannot be reproduced from the mature CSVs

Tex lines 247 to 251 (Results): the five surviving fraction shifts carry
q 0.0118 (four types) and q 0.0279 (multiciliated).
Source: `results/mature-only-composition-verdicts.csv`, column
mature_q_global_bh_orig_family: club cell 0.01178, MSSM L6b 0.01178,
oligodendrocyte 0.01178, OPC 0.01178, multiciliated 0.02785.
But the same file's refit-family q is 0.00858 for all five, and the
mature composition CSVs' own q_bh column gives 0.00091 for the three
MSSM types, 0.00424 for club, and 0.01022 for multiciliated
(`results/followup-mature-composition-mouse.csv`, 114 rows;
`results/followup-mature-composition-human.csv`, 285 rows).
A reader checking the mature CSVs cannot reproduce 0.0118. Reporting the
larger original-431-family q is conservative, but the family must be
stated in the text. Fix: one sentence naming the family for every
mature composition q in the paper.

### M2. Marrow is presented as a carrier of the concordance signal but fails BH correction across the seven tissues

Tex line 309: "The signal is carried by limb muscle (observed 0.428
against a null mean of 0.190, p = 0.001) and marrow (0.292 against
0.189, p = 0.019); no other tissue passes."
Source: `results/followup-mature-concordance-mouse-strata.csv`:
limb muscle observed 0.42838, null 0.19019, p 0.000999; marrow observed
0.29225, null 0.18871, p 0.01898. BH across the seven tissues: limb
muscle q 0.007 (passes); marrow q 0.0665 (fails, since 0.019 times 7
divided by 2 is 0.0665). No per-tissue family is declared anywhere in
the tex, so "passes" is judged at nominal p below 0.05 across seven
tests. Fix: declare the per-tissue family or drop marrow as a carrier.

### M3. The concordance signal is fragile: the excess over the null is smaller than the MDE80, and most of the observed correlation is shared-donor structure

Tex line 309. Source: `results/followup-mature-concordance-mouse.csv`
(observed 0.21047, null 0.14646, excess 0.0640) and
`results/followup-mature-concordance-mde80.csv` (mouse within mature:
observed minus null 0.0640, MDE80 0.0677, detectable threshold 0.2142).
The observed excess (0.064) is smaller than the minimum effect the test
was powered to detect at 80 percent (0.068). The joint null mean (0.146)
is 70 percent of the observed 0.210, so most of the correlation is
donor sharing across cell types, not age concordance. n is 14 mice with
1,000 permutations. This is a just-significant excess from an
underpowered test, and the text reports it without that context.

### M4. The naive CD4 headline slope is uncalibrated: finite-N calibration explains 56.4 percent of it

Tex line 277 (Results) and line 37 (Abstract): primary slope plus
0.00071 residual dispersion units per decade (p 5.174e-08,
q 8.279e-07, 777 donors).
Source for the primary numbers:
`results/followup-mature-block1-residual-dispersion.csv`, Yazar naive
thymus-derived CD4-positive alpha-beta T cell, primary arm: slope
0.00071, p 5.174e-08, q 8.279e-07, n 777. Those numbers are correct.
Source for the calibration: `reports/followup-mature-block1.md`,
section 5.7. Non-parametric subtraction of the real-cell finite-N
calibration curve accounts for 56.4 percent of the primary slope
(induced slope 0.000401 per decade, p 8.40e-07) and leaves a debiased
slope of 0.000313 (p 0.000638), and 0.000295 after depth adjustment
(p 0.00106). Adding log(N) leaves 0.000286 (p 0.00660); plus N plus depth
leaves 0.00024 (p 0.01893). The unadjusted arm (0.00003, p 0.8142) and
the rarefied 1500 arm (0.00008, p 0.7971) are null. The report's own
interpretation is that "the verdict rests on the frozen rule, not on a
resolved mechanism." The tex reports only the uncalibrated primary slope
and presents the signal as a clean surviving finding. Fix: report the
debiased slope (0.000313) beside the primary slope and state that 56
percent of the primary slope is finite-N calibration bias.

### M5. Two BH-significant Block 1 gene-fraction signals are never reported

Tex line 291: "Residual dispersion gives no dataset-level positive: Gao
does not prove and Yazar is inconclusive."
Source: `results/followup-mature-block1-gene-fractions.csv`: Gao Myeloid,
fraction of genes higher in old 0.6658, two-sided permutation p 0.002999,
q across cell types 0.02699; pooled Yazar All_Cells, fraction 0.5487,
p 0.002999, q 0.002999. Both pass BH. Neither appears anywhere in the
tex. The Yazar All_Cells q comes from a single-test family (q equals p),
which should be stated, but it is a declared family in the file and the
result passes it. The Gao Myeloid signal is plausibly ambient-RNA
artifact, but the tex never says so either; it is simply absent.
A prespecified secondary analysis with two BH-significant results cannot
be omitted while the text claims no dataset-level positive.

### M6. Composite q values are not calibrated and must not be quoted as evidence strength

The composite test is a parametric linear-model p value on a composite
score. The empirical test is the matched random-gene-set null. They
disagree by orders of magnitude on the same data, which means the
composite null model is wrong for these scores. Examples:
`results/followup-mature-block10-benchmark-results.csv`, Yazar All_Cells,
batch_adjusted: Harries composite q 7.51e-47 against empirical q 0.00425;
Peters composite q 2.60e-33 against empirical q 0.000654.
`results/followup-mature-organ-transports-c-gtex-track2-p1.csv`, Ovary,
P1m: composite q 1.27e-16 with empirical p 0.52. The composite test
alone would call the ovary effect overwhelmingly significant while the
matched null says it is chance.
Tex lines 265, 267, and 273 quote composite q values (1.29e-22,
7.51e-47, 4.03e-10) as if they measured evidence. The both-pass rule
protects the verdicts, but the quoted numbers mislead a reader about
strength. Fix: quote the empirical q as the calibrated number and give
the composite q only as the screening statistic.

### M7. The 200-draw empirical null cannot pass by construction, so it is not a sensitivity result

Tex line 267: "The 200-draw null kills every empirical pass (0 of 27,
against 4 at 2,000 draws)."
Source: `results/followup-mature-psychad-m10-null200.csv`: 0 pass at 200
draws, 4 pass at 2,000 draws (3 P0m, 1 P1m); 4 verdict changes,
all PASS_2000_ONLY. The counts are correct. The inference is not. With
200 draws the minimum achievable empirical p is 1 divided by 201, which
is 0.00498, and BH over 27 units gives a minimum q of 0.134, above
0.05. No unit can pass the 200-draw empirical test at q below 0.05,
whatever the data. The same floor applies to the PBMC 200-draw null,
where 26 of 35 rows sit at the 1/201 floor. Presenting the 200-draw
result as a sensitivity check that "kills" passes implies the data
failed a test; in fact the test was designed so it could not pass.
Fix: drop the 200-draw comparison or state explicitly that it cannot
pass by construction.

### M8. The naive CD8 scatter signal has undisclosed batch dominance and jackknife bias

Tex line 273: SD slope 0.0498 per decade, p 5.92e-05, q 4.74e-04,
344 donors, 17,878 cells.
Source: `results/followup-mature-fulldepth-yazar-cell-drift.csv`, part 1,
metric scatter, group naive thymus-derived CD8-positive alpha-beta
T cell, depth full_depth, model primary: slope_sd 0.0498023,
p 5.9205e-05, q 0.0004736, n_donors 344, n_cells 17878, MDE80 SD 0.03434.
All quoted numbers are correct. What the tex omits from the same row:
batch_partial_r2 0.8172 (batch explains 82 percent of the outcome
variance) and quenouille_bias_frac 0.2807 (the jackknife bias is 28
percent of the slope). The signal survives batch adjustment, but it
rides on an outcome that is overwhelmingly batch-driven, and the tex
discloses batch sensitivity for other signals while saying nothing here.
Fix: disclose both numbers in the same sentence.

### M9. Spleen is omitted from the per-organ drift list although it has the largest slope

Tex line 255: "per-organ decompositions disagree in sign: lung minus
0.194 (p = 0.021), limb muscle plus 0.286 (p = 0.043), kidney minus
0.231 (p = 0.057), none corrected as a family."
Source: `results/followup-mature-drift-mouse.csv`, full depth: spleen
slope plus 0.4528, p 0.0566; at 1,500 UMIs p 0.0076; at 3,000 UMIs
p 0.0054. Spleen has the largest slope of any tissue and is the only
tissue significant at fixed depths, yet it is not listed. The global
fixed-depth arms are null as stated, but the per-organ list is
selective. Fix: list spleen with the others, or state the selection rule
for the list.

### M10. The fourteen full-depth primary endpoints are never enumerated, and the stratified globals carry no q values

Tex line 299: "All fourteen primary endpoints do not prove in Gao heart
and Yazar blood."
The 14 endpoints are not enumerated anywhere in the tex, so the claim
cannot be checked from the manuscript. Source:
`results/followup-mature-fulldepth-gao-cell-drift.csv`,
`results/followup-mature-fulldepth-yazar-cell-drift.csv`,
`-spread.csv`, `-ectopic.csv`: every Stratified_CellType_Global row has
an empty q_bh_within_slice. The "do not prove" verdict therefore rests
on uncorrected p values plus an unstated directional rule; no
full-depth verdict rule appears in Table 3. Meanwhile the per-type rows
in the same files contain BH-significant passes the tex never mentions
except naive CD8, all in the falling direction: Yazar natural killer
scatter q 1.33e-07, CD14-positive monocyte scatter q 0.00360,
effector memory CD8 center shift q 0.00215, natural killer center shift
q 0.00215. The paper reports the one positive per-type exception and
silently drops the equally strong negative per-type signals. Fix:
enumerate the 14 endpoints, state the verdict rule, and report the
per-type BH-significant falls.

### M11. The Discussion 83 percent and 96 percent numbers are full-range, use changing test universes, and have no table citation

Tex line 347: "About 83 percent of gene hits did not persist after pool
adjustment and concordance fell by about half in the 981-donor blood
cohort, while over 96 percent did not persist after cohort adjustment
in the 99-donor cohort."
Sources: `reports/followup-hcata-pmc12511636.md` lines 61 and 89;
`results/followup-hcata-pool-sensitivity-genes-summary.csv`;
`results/followup-hcata-metadata-summary.json`. Yazar: 10,123
significant of 160,239 tests unadjusted, 1,741 of 134,036 after pool
adjustment; 1 minus 1741 divided by 10123 is 82.8 percent. Concordance
0.2430 to 0.1024. Perez: 99 donors, 10,021 to 391; 1 minus 391 divided
by 10021 is 96.1 percent. The percentages are arithmetically correct.
Three problems: both analyses use the full age range (Yazar 19 to 97,
Perez 21 to 75), not the 30-plus restriction the paper imposes everywhere
else; the test universe shrinks between arms (160,239 to 134,036 tests),
so part of the attenuation is donor and gene loss rather than
confounding removal; and the tex gives no table or file pointer for
either number. Fix: scope the numbers to the full range, give the
denominators, and cite the source table.

### M12. The RPE survival rests on 13 donors and a single-unit family

Tex line 265: "Retinal pigment epithelium is the one earlier both-pass
transport that survives under both tiers, in the negative direction
(effect minus 0.0211, q = 0.0443, empirical q = 0.0005)."
Source: `results/followup-mature-organ-transports-c-rpe.csv`,
transport_P1m row: n_donors 13, residual_df 8, mature_effect -0.02113,
mature_q 0.04432, mature_empirical_q 0.00050. The composite q equals
the p value because the family has one unit, and 0.0443 is marginal.
The tex should state the n and the single-unit family beside the claim.

### M13. The Gao Block 7 hits run opposite the prediction, and 14 of the 96 predicted gene changes did not attenuate

Tex line 297: "Gao gives two hits with zero predicted genes each."
Source: `results/followup-mature-block7-death-stress.csv`,
test1_primary_full, model primary: Gao Lymphoid slope -0.0292
(q 0.0318); Gao Pericyte slope -0.0242 (q 0.0378). Both are negative:
the high-stress-cell fraction FALLS with age in Gao, opposite the
survivor-bias prediction. The tex presents them neutrally as hits.
Source for the attenuation gap:
`results/followup-mature-block7-removal-sensitivity.csv`, summary rows:
Yazar central memory CD4 37 predicted with 4 not attenuated, naive CD4
56 predicted with 9 not attenuated, naive CD8 3 predicted with 1 not
attenuated; total 96 predicted, 14 not attenuated. The tex reports 96
without the 14. Fix: state the Gao direction and the 14 non-attenuated
genes.

### M14. The secondary-five q values printed as 0.050 hide a marginal pass

Tex line 251: two refit-family q values printed as "0.050".
Source: `results/followup-mature-composition-secondary-checks.csv`:
both are exactly 0.04989 (lung multiciliated columnar and SEA-AD L5).
Printing 0.050 reads as exactly at the threshold and hides that both
pass by 0.0001. Fix: print 0.0499.

## MINOR issues

### m1. Interaction bootstrap p values are unmentioned

Appendix A, line 449: mouse minus human interaction minus 0.0663
(p 0.0515) and minus 0.0447 after depth adjustment (p 0.0051).
Source: `results/reviewer-rev2-mouse-human-interaction.csv`: bootstrap
p values from 1,000 reps are 0.072 and 0.03. The depth-adjusted
bootstrap p (0.03) is six times the analytic p (0.0051). The tex
reports only the analytic p.

### m2. Table 1 says 479 vector donors; Appendix A says 480 donors for human concordance

Table 1, line 80: "479 vector donors". Appendix A, line 453: human
concordance "across 131 units and 480 donors". Source:
`results/h3-cross-species-vector-global.csv`, Human all cohorts,
n_donors 479; `results/reviewer-rev2-joint-concordance-human.csv`,
n_donors 480. One donor differs between the vector and concordance
eligibility. This is not an error, but the two numbers sit one page
apart with no note. Fix: one clarifying phrase.

## Items not rechecked in this round

I did not recheck: the HRCA 0 of 144 count and the amacrine fraction
claim (line 265); the Perez 6 of 10, 2 of 10, and 9 of 10 counts and the
Guilliams 0 of 8 (line 265); the full-range transport three passes and
the heart endocardium plus 0.0879 and five blood types (Appendix A);
the P1 118 human-mapped count; the method sensitivity counts 151, 312,
and 1 and the leave-one-donor-out 310 of 765 (Appendix A); the
off-program depth sweep table values and MDE80 0.6191 (Appendix A);
the nonlinear, annotation, and marker-panel numbers (Appendix A);
the exploratory vector numbers for fly, worm, zebrafish, and naked mole
rat (line 321); and the figure data files behind Figures 1 to 4.
The checker script passed (abstract 270 words, four figures, 68 cited
references, no forbidden dashes, no process residue), and the four
figure files exist in `report/plos-one/figures-v2/`.

## Verdict: MAJOR REVISION

The analysis itself is largely sound: I verified on the order of one
hundred numbers against their CSVs and nearly all match, including the
hardest ones (the 624 family, the transport counts, the drift nulls, the
Block 7 counts, the PsychAD numbers, the GTEx numbers). The paper does
not need new analysis. It needs its text made consistent with its own
statistics. Three issues matter most:

1. The concordance contradiction (B1). The paper calls a BH-passing
result "nominal only" in Results, Discussion, and Conclusions. The
family choice decides the headline "no conserved direction" claim for
mouse, and both families the paper declares pass. Either the claim or
the family must change.
2. The Methods self-contradiction on family pre-specification (B2).
Line 226 claims every family was defined before correction; line 143
admits the primary 624 family was selected after both correction
outcomes were available. One of the two sentences is false.
3. The ectopic headline numbers come from a sensitivity variant while
the primary endpoint is null in both cohorts (B3). The text reports
significant falls of minus 0.127 and minus 0.0327 without saying they
are the E2_stratum bounded-depth sensitivity; the primary E1 endpoint
is minus 0.0633 (p 0.468) in Gao and minus 0.0187 (p 0.073) in Yazar.

Fix B1, B2, and B3, declare the mature composition q family (M1),
report the naive CD4 debiased slope (M4), and either report or explain
away the unreported BH-significant Block 1 gene-fraction signals (M5).
The remaining MAJOR items are one-sentence fixes each. After those
changes the manuscript is submittable on statistical grounds.

Counts: 3 BLOCKER, 14 MAJOR, 2 MINOR.
