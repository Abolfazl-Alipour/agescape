# Bayesian refits for the mature Agescape analyses

Fitting script: `scripts/run_mature_bayesian.py` (closed-form conjugate fits,
items 1-4 and 6, item 5 mouse carry-over) and
`scripts/run_mature_bayesian_human_drift.py` (PyMC hierarchical fit, item 5
human extension, run on beast). Per-result CSVs:
`results/followup-mature-bayesian-composition.csv`,
`results/followup-mature-bayesian-narrowing.csv`,
`results/followup-mature-bayesian-naive-cd4.csv`,
`results/followup-mature-bayesian-transport.csv`,
`results/followup-mature-bayesian-drift.csv`,
`results/followup-mature-bayesian-concordance.csv`,
`results/followup-mature-bayesian-human-drift.csv` (posterior rows),
`results/followup-mature-bayesian-human-drift-diagnostics.csv`,
`results/followup-mature-bayesian-human-drift-frequentist.csv`.
Branch: `bayesian-additions`. Nothing in the tex, figures, or Q&A was edited
and no number in the paper was changed. Donor is the unit in every model;
the cell is never the unit.

## 1. Priors (fixed before fitting)

| Parameter | Primary prior | Sensitivity priors | Reason |
| --- | --- | --- | --- |
| Standardized age slope | Normal(0, 0.5) | Normal(0, 1.0), Normal(0, 0.25) | Same family as the tex "Bayesian refits" subsection; weakly informative, rules out absurd slopes, lets the data dominate anywhere plausible |
| Nuisance coefficients (intercept, sex, depth, cell-type dummies, batch fixed effects) | Normal(0, 1e4 x sigma2), near-flat | Same | Nuisance parameters carry no prior information |
| Residual variance sigma2 | Jeffreys p(sigma2) ~ 1/sigma2 | Same | Standard conjugate choice |
| Human drift, hierarchical (item 5) | Normal(0, 0.5) on global age slope; Normal(0, 0.5) on sex/depth; Gamma(2, 2) on group/unit/donor scales and unit sigma; Gamma(2, 0.1) on nu | Normal(0, 0.25), Normal(0, 1.0) on the depth-adjusted global slope | Same model as `scripts/bayesian_diagnostics_helper.py::build_hierarchical_student_t` used for the mature mouse refit |

Standardization. Design-matrix models divide the outcome by its own SD and
the age predictor by its SD, so the Normal(0, 0.5) prior sits on the fully
standardized slope (outcome SDs per predictor SD); posteriors are converted
back to the raw slope (per year or per decade). Summary-statistic models have
no design matrix, so the prior sits on the slope in outcome SDs per predictor
unit (per decade); the source CSVs supply the outcome SD (Shannon rows carry
outcome_sd = 1.0; the naive CD4 row carries outcome_sd_ref = 0.005688). The
one exception is PsychAD Micro (item 4), which has no committed per-donor
scores and therefore no outcome SD: its scale is borrowed from the MSSM
microglia P1m score SD (same score family, same tissue, same cell type),
with a half/double sensitivity reported in section 5.

ROPE (region of practical equivalence), defined before fitting. The ROPE
half-width is the slope that moves the outcome by less than 0.05 outcome SDs
over the paper's reporting span for that slope: one decade for items 1-4
(item 1 slopes are per year, so the ROPE is 0.05 x sd_y / 10 per year), 0.4
years (mouse) and 12.25 years (human) for item 5, and 0.05 correlation points
for item 6. ROPE sensitivity uses half width (0.025) and double width (0.10).
HDI-plus-ROPE rule: accept the null only when the full 95% interval falls
inside the ROPE, reject it only when the full interval falls outside it,
otherwise undecided. All intervals below are 95% equal-tailed credible
intervals from exact posterior draws, except the carried-over mouse drift
row, which reports the committed 94% highest-density interval (HDI).

Diagnostics. Closed-form fits have no sampler: the diagnostic is the
second-seed rerun (seed 20261123 vs 20261122). Largest move in any posterior
median across all closed-form fits: 0.00014 (item 4, Lung alveolar
macrophage, on a slope whose 95% interval is 0.09 wide); all other moves are
below 0.00006 in raw units. The OLS slope from the same design is reported
beside every Bayesian fit; items 1 and 2 reproduce the committed frequentist
slopes exactly (item 2 Gao OLS -0.0723, Yazar OLS -0.0197, matching the
frequentist stratified slopes to 4 decimals).

## 2. Headline numbers

### 2.1 Item 1: five headline fraction shifts (CLR on age + sex)

| Unit | n | Frequentist slope/yr (q) | Bayesian median/yr | 95% CI | P(claimed direction) | P(ROPE) |
| --- | --- | --- | --- | --- | --- | --- |
| Human lung club cell | 91 | -0.0804 (0.0118) | -0.0770 | -0.1144 to -0.0396 | 1.0000 | 0.001 |
| Human lung multiciliated epithelial cell | 91 | -0.0470 (0.0279) | -0.0450 | -0.0697 to -0.0203 | 0.9998 | 0.002 |
| MSSM L6b glutamatergic neuron | 209 | +0.0198 (0.0118) | +0.0193 | +0.0099 to +0.0288 | 1.0000 | 0.002 |
| MSSM oligodendrocyte | 209 | +0.0173 (0.0118) | +0.0169 | +0.0084 to +0.0254 | 0.9999 | 0.002 |
| MSSM oligodendrocyte precursor cell | 209 | -0.0119 (0.0118) | -0.0117 | -0.0174 to -0.0059 | 1.0000 | 0.002 |

Bayesian and frequentist agree on all five: same sign, same size (Bayesian
medians are 2-4% closer to zero than the OLS slopes, the expected shrinkage
from the Normal(0, 0.5) prior), every 95% interval excludes zero, every
direction probability is at least 0.9998, and every ROPE probability is at
most 0.002. Direction probabilities are unchanged under both sensitivity
priors (all stay above 0.9994). No disagreement.

### 2.2 Item 2: narrowing (scatter and Shannon entropy)

Stratified scatter model (per-decade slopes, outcome standardized within
cell type; sex and depth centered within cell type; cell-type dummies;
study/pool fixed effects; eligible cell types only, i.e. groups whose
per-type primary model kept at least one batch dummy, replicating the
original eligibility rule):

| Unit | n rows (donors) | Frequentist slope (p) | Bayesian median | 95% CI | P(<0) | P(ROPE) |
| --- | --- | --- | --- | --- | --- | --- |
| Gao scatter | 829 (105) | -0.0723 (0.041) | -0.0719 | -0.1053 to -0.0383 | 1.0000 | 0.101 |
| Yazar scatter | 5638 (925) | -0.0197 (6.69e-05) | -0.0196 | -0.0268 to -0.0125 | 1.0000 | 1.000 |

Shannon entropy (summary-statistic model on the Stratified_CellType_Global
rows; note the CSV model_spec is `stratified_global_primary`, not `primary`):

| Unit | n donors | Frequentist slope (p) | Bayesian median | 95% CI | P(<0) | P(ROPE) |
| --- | --- | --- | --- | --- | --- | --- |
| Gao Shannon, primary full-depth | 105 | -0.1001 (0.0539) | -0.0988 | -0.1998 to +0.0012 | 0.9736 | 0.168 |
| Gao Shannon, depth-by-cell-type spec | 105 | -0.1013 (0.0344) | -0.1002 | -0.1932 to -0.0080 | 0.9835 | 0.142 |
| Yazar Shannon, primary full-depth | 925 | -0.0298 (0.00365) | -0.0297 | -0.0499 to -0.0097 | 0.9982 | 0.976 |
| Yazar Shannon, depth-by-cell-type spec | 925 | -0.0274 (0.00701) | -0.0273 | -0.0474 to -0.0075 | 0.9965 | 0.987 |

Disagreements with the frequentist answer (see section 3):

1. Gao Shannon, primary spec: frequentist p = 0.0539 (not significant at
   0.05); Bayesian P(<0) = 0.9736, but the 95% interval still crosses zero
   (-0.1998 to +0.0012). Both answers say the same thing in different words
   (weak evidence, interval touches zero); the Bayesian number quantifies it
   as 97% probability negative rather than a p value just above 0.05.
2. Yazar scatter and Yazar Shannon: the direction is certain (P(<0) =
   1.0000 and 0.9982) but the slopes are small in absolute terms (-0.0197
   and -0.0297 SDs per decade), so under the pre-stated ROPE of 0.05 SDs per
   decade both sit inside the ROPE (P(ROPE) = 1.000 and 0.976). The
   frequentist p values (6.69e-05 and 0.00365) call these significant; the
   Bayesian ROPE verdict is that the per-decade effect is practically null,
   while the cumulative move over 50 years is about 0.10 SDs (scatter) and
   0.15 SDs (Shannon). The ROPE verdict is width-sensitive in the expected
   direction (at half width 0.025, Yazar scatter P(ROPE) = 0.930, Yazar
   Shannon P(ROPE) = 0.323), so the robust statement is the direction
   probability, not the ROPE label. Caveat: Yazar's MDE80 is 0.0138 SDs per
   decade for scatter, so the 0.05 ROPE is about 3.6x the smallest effect
   Yazar was powered to detect; the ROPE verdict for Yazar is driven by the
   wide band, not by uninformative data. Gao's scatter slope (-0.0719) is
   3.7x Yazar's and mostly outside the ROPE (P(ROPE) = 0.101).

Limitation. The conjugate stratified fit treats donor-by-cell-type rows as
independent given the design, while the frequentist fit clusters standard
errors by donor (the same donor appears in many cell types). The Bayesian
intervals are therefore narrower than the clustered frequentist intervals:
Gao posterior SD implies a 95% width of 0.067 vs the clustered frequentist
width of about 0.137 (cluster SE 0.0349); Yazar posterior width 0.014 vs
clustered width about 0.019 (cluster SE 0.0049). The direction
probabilities are not sensitive to this (both are 1.0000), but the Gao
interval width should be read with the clustering caveat in mind.

### 2.3 Item 3: naive CD4 dispersion

| Arm | Frequentist slope/decade (p) | Bayesian median | 95% CI | P(>0) | P(ROPE) |
| --- | --- | --- | --- | --- | --- |
| Primary | 0.000714 (5.17e-08, q 8.28e-07) | 0.000713 | 0.000457 to 0.000967 | 1.0000 | 0.001 |
| Debiased | 0.000313 (0.000638) | 0.000313 | 0.000132 to 0.000492 | 0.9997 | 0.377 |

The primary slope agrees exactly with the frequentist answer. The debiased
slope exists only in `reports/followup-mature-block1.md` section 5.7, in no
CSV, and no SE is committed anywhere for it. Its SE was derived here as
SE = slope / z with z = Phi^-1(1 - p/2) = 3.4149 for the two-sided
p = 0.000638, giving SE = 9.18e-05 per decade. That derivation is the only
input the debiased Bayesian fit adds beyond the committed slope and p value;
the posterior is then a standard Normal-likelihood-with-known-SE fit. The
debiased posterior stays positive (P(>0) = 0.9997) but its ROPE probability
is 0.377: the debiased slope (0.000313) sits almost exactly on the ROPE edge
(0.05 x outcome SD = 0.000284 per decade), so the data cannot separate "real
but small" from "practically null" for the debiased arm. At double ROPE
width P(ROPE) = 0.997; at half width it is 0.032. No prior sensitivity flip:
P(>0) stays above 0.9996 under both sensitivity priors for both arms.

### 2.4 Item 4: transport null

ROPE, defined before fitting: a score slope per decade that moves the
composite score by less than 0.05 within-unit score SDs (per-unit ROPE =
0.05 x that unit's score SD). ROPE sensitivity uses 0.025x and 0.10x SD.

| Unit | n | Frequentist slope/decade (p) | Bayesian median | 95% CI | P(claimed dir) | P(ROPE) |
| --- | --- | --- | --- | --- | --- | --- |
| MSSM microglial cell | 196 | -0.00830 (0.00205, q 0.0497) | -0.00813 | -0.01326 to -0.00297 | 0.9989 | 0.017 |
| PsychAD Micro | 201 | -0.01433 (8.06e-05, q 0.00457) | -0.01403 | -0.02114 to -0.00699 | 1.0000 | 0.001 |
| 10 null units (see CSV) | 46-925 | p 0.26-0.93 | -0.01516 to +0.02578 | all cross zero | 0.53-0.88 | 0.23-0.94 |

The two microglia passes agree with the frequentist answer: both slopes are
negative with 95% intervals excluding zero (MSSM P(<0) = 0.9989, PsychAD
P(<0) = 1.0000) and both sit outside the ROPE. No prior sensitivity flip.

The ten null units (3 MSSM, 3 HBCC, 2 Lung_Integrated from the P1m score
family, plus 1 gao and 1 yazar from the block8 file) all agree with the
frequentist null: every 95% interval crosses zero and every direction
probability is between 0.53 and 0.88. But under the HDI-plus-ROPE rule every
one of the ten is UNDECIDED, not accepted as null: no 95% interval falls
entirely inside the ROPE (closest is the yazar block8 CD14-low monocyte at
P(ROPE) = 0.935, whose interval still extends outside the ROPE), and none
falls entirely outside it either. The Bayesian answer therefore does not
upgrade these frequentist nulls to confirmed nulls; it says the data are
compatible with small effects in either direction. ROPE sensitivity does not
change any verdict: at double width the highest P(ROPE) is 1.0000 (yazar
block8) but its interval still crosses the ROPE boundary, so it stays
undecided under the strict rule.

Input discrepancy, recorded here because it changes where the numbers come
from. The task names `results/followup-mature-block8-donor-scores.csv` for
this item, but that file holds only gao and yazar rows: no PsychAD
(DLPFC_Lifespan_Atlas) and no MSSM rows exist in it, so the two microglia
passes cannot be fit from it. The per-donor P1m scores for the core cohorts
live in `results/followup-mature-p1m-m9-donor-scores.csv` (MSSM, HBCC,
kidney, lung, BAL), and the MSSM microglia fit uses that file. Per-donor
PsychAD scores are committed nowhere, so PsychAD Micro is fit as a
summary-statistic model from its committed slope (-0.001433 per year,
fully_adjusted model) and p value (8.06e-05) in
`results/followup-mature-p1m-transport.csv`, with SE = |slope| / z,
z = Phi^-1(1 - p/2) = 3.9426, giving SE = 0.003635 per decade, and the score
SD borrowed from the MSSM microglia P1m scores (0.0515). Borrowed-scale
sensitivity: P(<0) = 0.9999 at half scale and 1.0000 at double scale, so the
PsychAD verdict does not depend on the borrowed scale. The block8 file is
still used for two of the ten null units (one gao, one yazar). Only MSSM,
HBCC and Lung_Integrated have units with at least 40 donors in the per-donor
score file (kidney cohorts max 19 donors, Lung_CellRef max 10), so the
eight P1m nulls spread across those three cohorts.

### 2.5 Item 5: drift nulls

Mouse (carried over, verified from the committed file). The depth-adjusted
mature mouse refit in
`results/rev3-bayesian-vector-mature-only-posterior.csv`, spec
`mature_ge_18m-profile_hellinger_centroid_mean-depth_adjusted-prior_0.50`,
species scope: posterior mean +0.0389 per 0.4 years, 94% HDI -0.0998 to
+0.1773, P(positive) = 0.71 (file values 0.038944, -0.099779, 0.177334,
0.7068). Prior-0.25 sensitivity row: mean +0.0338, P(positive) = 0.6950. No
prior-1.00 mature mouse refit is committed, so the wide-prior sensitivity
for the mouse row is not available. Under the ROPE of 0.05 standardized
Hellinger SDs per 0.4 years, P(ROPE) = 0.449 (0.485 under prior 0.25):
undecided. The frequentist depth-adjusted slope is -0.0137 (p = 0.816), so
the two point estimates have opposite signs while both intervals are wide
and include zero; the tex already notes the two need not have the same sign.

Human (PyMC extension): NOT DONE - all five fitted specs failed
convergence, so no human Bayesian numbers can be read from them. The env
built on beast in about 2 minutes (well inside the 30-minute timebox) and
the run completed all five specs (unadjusted prior 0.50; depth-adjusted
priors 0.50, 0.25, 1.00; plus a second-seed rerun of depth-adjusted 0.50),
4 chains x 2000 kept draws each, on 4657 rows, 407 donors, 87 units. Every
spec fails the skill's convergence bar (split R-hat below 1.01, bulk ESS
above 400, zero divergences):

| Spec | Max R-hat | Min bulk ESS | Divergences | convergence_ok |
| --- | --- | --- | --- | --- |
| Unadjusted, prior 0.50 | 1.0711 | 44 | 50 | False |
| Depth-adjusted, prior 0.50 | 2.4782 | 5 | 26 | False |
| Depth-adjusted, prior 0.25 | 1.1848 | 16 | 52 | False |
| Depth-adjusted, prior 1.00 | 2.3317 | 5 | 163 | False |
| Depth-adjusted, prior 0.50, seed 2 | 1.2527 | 14 | 32 | False |

The sampler output is committed as an artifact
(`results/followup-mature-bayesian-human-drift.csv` and the diagnostics and
frequentist CSVs) so the failure is auditable, but its posterior medians
(e.g. depth-adjusted prior 0.50 median +0.0068 with mean +0.0659 and SD
0.241, vs prior 1.00 median -0.0226) disagree across chains and priors in
the way non-converged chains do, and they must not be cited as results. The
likely next step, NOT run here because the wrap-up order forbids starting a
new long run, is a rerun with the mouse script's stronger sampler settings
(target_accept 0.99, 2500 tuning draws, 2500 kept draws) and possibly a
reparameterization; the tex already notes human global models carried
convergence warnings. The frequentist check on the same prepared data
reproduces the drift report exactly (unadjusted +0.0626, p = 0.0816;
linear-depth -0.0071, p = 0.6106), which validates the data preparation even
though the Bayesian fit did not converge.


### 2.6 Item 6: mouse within-tissue concordance

Construction, exactly. Observed mean Pearson r = 0.21047 over n_pairs = 260
cell-type pairs; joint null mean 0.14646, SD 0.02724, from 1000 joint
permutations (one age assignment per mouse, shared across its cell types).
Fisher z of the observed correlation carries sampling error 1/sqrt(260-3)
under a flat prior on z. The null correlation is Normal(0.14646, 0.02724);
on the z scale its SE is 0.02724 / (1 - 0.14646^2) = 0.02784 by the delta
method. The estimand is delta_z = z_obs - z_null with a Normal(0, tau^2)
prior (tau = 0.5 primary, 1.0 and 0.25 sensitivity). delta_r is recovered by
drawing z_null from its Normal and back-transforming:
delta_r = tanh(z_null + delta_z) - tanh(z_null). Frequentist reference:
joint p_greater = 0.02198, q = 0.0440 in the two-test family.

| Estimand | Bayesian median | 95% CI | P(>0) | P(ROPE) |
| --- | --- | --- | --- | --- |
| delta_r (observed minus null correlation) | +0.0632 | -0.0682 to +0.1804 | 0.8319 | 0.374 |
| delta_z | +0.0653 | -0.0685 to +0.1978 | 0.8319 | - |

Disagreement with the frequentist answer, in strength but not direction.
The permutation test rejects the joint null at p = 0.022 (q = 0.044). The
Bayesian posterior puts only 83% probability on a positive difference, and
the 95% interval runs from -0.068 to +0.180. The gap comes from what each
method conditions on: the permutation p value treats the observed
correlation as fixed and asks how often the joint null produces something
that extreme, while the Bayesian construction also propagates the sampling
error of the observed correlation itself (SE about 0.062 on the z scale from
260 pairs) and the spread of the null distribution. With only 260 pairs the
observed correlation is noisy, and that noise dominates the posterior width.
Under both sensitivity priors P(>0) stays between 0.825 and 0.834, so the
conclusion is prior-robust: the data lean positive but do not establish the
difference. Under the HDI-plus-ROPE rule the result is undecided (the
interval straddles both the ROPE boundary and zero).

## 3. Every place the Bayesian answer disagrees with the frequentist answer

1. Item 2, Yazar scatter: frequentist p = 6.69e-05 (significant); Bayesian
   direction is certain (P(<0) = 1.0000) but the slope (-0.0196 SDs per
   decade) sits inside the pre-stated ROPE, so the per-decade effect is
   practically null. Same for Yazar Shannon (p = 0.00365; P(<0) = 0.9982;
   P(ROPE) = 0.976). The disagreement is about magnitude, not direction,
   and the ROPE verdict is band-dependent (see section 2.2).
2. Item 2, Gao Shannon primary spec: frequentist p = 0.0539 (just above
   0.05); Bayesian P(<0) = 0.9736 with a 95% interval that still crosses
   zero. Both say weak evidence; the Bayesian number quantifies it.
3. Item 3, debiased naive CD4: frequentist p = 0.000638 (significant);
   Bayesian P(>0) = 0.9997 but P(ROPE) = 0.377, so the debiased slope sits
   on the ROPE edge and the data cannot separate "real but small" from
   "practically null".
4. Item 4, ten null units: frequentist p values 0.26-0.93 (null); Bayesian
   intervals all cross zero, agreeing, but under the HDI-plus-ROPE rule all
   ten are undecided rather than accepted as null. The Bayesian answer does
   not upgrade these nulls to confirmed nulls.
5. Item 5, mouse drift: frequentist depth-adjusted slope -0.0137
   (p = 0.816); Bayesian posterior mean +0.0389, P(positive) = 0.71.
   Opposite signs on the point estimates; both intervals are wide and
   include zero.
6. Item 6, concordance: frequentist joint p = 0.022 (q = 0.044, rejects);
   Bayesian P(>0) = 0.83 with a 95% interval from -0.068 to +0.180. The
   Bayesian answer is more conservative because it propagates the sampling
   error of the observed correlation.

No disagreement in items 1 (all five headline shifts agree) or in the two
microglia transport passes (both agree).

## 4. What could not run, and why

1. Item 4, PsychAD Micro per-donor fit: not runnable as a design-matrix
   model. No per-donor PsychAD (DLPFC_Lifespan_Atlas) composite scores are
   committed anywhere in the repo (checked: the only per-donor score files
   are `results/followup-mature-block8-donor-scores.csv` with gao/yazar only
   and `results/followup-mature-p1m-m9-donor-scores.csv` with MSSM/HBCC/
   kidney/lung/BAL only). It was run instead as a summary-statistic model
   from its committed slope and p value, with the SE derived as |slope|/z
   and the score SD borrowed from MSSM microglia (section 2.4).
2. Item 4, block8 file as the sole input: not runnable as specified. The
   named file contains no PsychAD or MSSM rows, so the two microglia passes
   cannot be fit from it. The MSSM pass uses the P1m score file; the block8
   file contributes two of the ten null units.
3. Item 2, Shannon per-donor fit: not run. The task states there is no
   per-donor Shannon column in the donor CSVs; in fact the part2 donor files
   carry a `full_shannon` column per donor-by-cell-type row. The
   summary-statistic model on the Stratified_CellType_Global rows was run as
   specified, and no per-donor Shannon refit was added. Note: the CSV
   model_spec value is `stratified_global_primary`, not `primary` as written
   in the task.
4. Item 3, debiased naive CD4 SE: no SE is committed anywhere. The SE was
   derived from the slope and p value (section 2.3); the derivation is
   documented, not hidden.
5. Item 5, mouse prior-1.00 sensitivity: no prior-1.00 mature mouse refit
   is committed in `results/rev3-bayesian-vector-mature-only-posterior.csv`
   (only 0.25 and 0.50), so the wide-prior sensitivity for the mouse row is
   not available. The narrow-prior (0.25) row is reported instead.
6. Item 5, human PyMC extension: ran to completion but NOT DONE. All five
   specs failed convergence (max R-hat 1.07 to 2.48, min bulk ESS 5 to 44,
   26 to 163 divergences; skill bar is R-hat below 1.01, ESS above 400,
   zero divergences), so no human Bayesian numbers are reportable. The
   sampler output is committed as an auditable artifact. A rerun with
   stronger sampler settings was not started, per the wrap-up order.

## 5. Reproduction

```sh
python3 scripts/run_mature_bayesian.py
# human PyMC extension (beast, agescape-bayesian env):
#   copy scripts/run_mature_bayesian_human_drift.py,
#   scripts/bayesian_diagnostics_helper.py and
#   results/followup-mature-drift-human-donors.csv to beast, then run
#   ~/miniconda3/envs/agescape-bayesian/bin/python run_mature_bayesian_human_drift.py
```

Seeds: 20261122 (primary), 20261123 (second-seed check, closed-form fits);
human PyMC seed is `sum(map(ord, "agescape-mature-human-drift-bayes-v1"))`
with a second-seed rerun of the primary depth-adjusted spec.
