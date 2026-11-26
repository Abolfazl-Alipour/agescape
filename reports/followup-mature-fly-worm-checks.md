# Sensitivity checks on the mature-to-old fly and worm fraction rows

Branch: `followup-mature-fly-worm-checks` (from master `db40edc`).
Date: 2026-11-25. The scan, the tex, and master are untouched. No tex change, no deck change.

## What was checked

The pooled mature fraction scan (`results/mature-all-species-composition.csv`,
577 tests, 39 pass) has 26 non-mammalian passing rows: 25 fruit fly
(AFCA, body and head strata) plus 1 worm (ARSC, Gao GSE229022). The mammalian
rows already went through propeller, total-cell refit, leave-one-donor-out,
and study-aware checks. This report runs the same checks for fly and worm
wherever possible. The inferential unit for fly and worm is the POOL, not the
animal: fly 72 pools (36 head, 36 body; ages 30, 50, 70 days, 6 male and
6 female pools per age per stratum); worm 9 pools (days 6, 12, 14; 3 pools
per age). Mature windows: fly 30 days and older, worm 6 days and older.
Slopes are CLR slopes per day.

## 1. Starting point verified

The pooled BH q values were recomputed independently from the scan p-value
column (statsmodels `multipletests`, BH). Max absolute difference from the
committed `q_pooled_mature_bh` column: 3.1e-15. Exactly 39 rows pass at q
below 0.05, split 2 mouse, 11 human, 25 fly, 1 worm, and the pass flags agree
exactly with the committed `significant_pooled_05` column.

## 2. Study-aware check: IMPOSSIBLE for both species

The lung study-aware check worked only because Lung_Integrated pools donors
from multiple studies, so the effect could be re-tested inside each study.
No equivalent structure exists here.

Fly evidence, recomputed from the loaded AFCA obs
(`scripts/run_mature_fly_worm_analogs.py`, check `study_evidence`):

- Head file: 289,981 cells; dataset AFCA 189,670 cells, FCA 100,311 cells.
  Every FCA cell is age 5. After the mature filter (dataset AFCA, age 30 and
  older), 36 pools remain and every cell is dataset AFCA.
- Body file: 276,273 cells; dataset AFCA 179,679 cells, FCA 96,594 cells.
  Every FCA cell is age 5. After the mature filter, 36 pools remain and
  every cell is dataset AFCA.
- No study or batch column exists in the AFCA obs. The only other dataset in
  the files is FCA, which is entirely age 5 and excluded from the mature
  scan because study and age are inseparable there.

Worm evidence, recomputed from the loaded worm obs:

- One integrated file from one GEO series (GSE229022). The obs has no study,
  dataset, or batch column. Mature pools are wild-type N2 only. The
  `orig.ident` values are sample IDs inside the one study, not studies.

No substitute was faked and called study-aware. The two closest analogs are
in section 5, clearly labeled as analogs.

## 3. Propeller

Method: the repo's propeller run is not the classic two-group propeller test.
`scripts/run_reviewer_phase1_propeller.R` fits limma on logit-transformed
proportions with age (continuous) plus log total cells plus sex when it
varies, per cohort. It was reused exactly, with no change. Inputs were built
by `scripts/prepare_mature_fly_worm_propeller.py` (same layout as
`data/interim/mature-propeller`: one subdir per cohort with `counts.csv` and
`metadata.csv`, plus `manifest.csv`), reusing `load_strata()` from
`scripts/run_mature_all_species_fractions.py`. Three cohorts: fly head
(36 pools, 90 cell types), fly body (36 pools, 72 cell types), worm whole
animal (9 pools, 16 cell types). R, limma, and speckle were already installed
on Beast; nothing was installed. Output:
`results/followup-mature-fly-worm-propeller.csv`, 178 rows, same columns as
`results/followup-mature-propeller.csv` plus `propeller_q_pooled_26` (BH
across the 26 passing rows only, NA elsewhere).

Design notes: fly propeller fits include age, log total cells, and sex (sex
varies). Worm propeller fits include age, log total cells, and
`platform_nextseq` (it varies); sex is dropped because it is constant (all
worm pools are female/hermaphrodite).

Results for the 26 passing rows:

- Within-cohort BH q below 0.05: 23 of 26.
- Pooled-across-the-26 BH q below 0.05: 24 of 26.
- Propeller effect sign agrees with the scan sign: 24 of 26.

The three propeller failures (within-cohort q, the bar matching the
mammalian precedent):

1. `hemocyte_body` (fly body): propeller effect +0.0001 logit, p 0.960,
   within-cohort q 0.974, pooled q 0.960. The independent method finds
   nothing. The scan slope is -0.0110 per day.
2. `ejaculatory bulb` (fly body): propeller effect -0.0043 logit (same sign),
   p 0.301, within-cohort q 0.417, pooled q 0.313. The independent method
   finds nothing.
3. `adult ventral nervous system` (fly body): propeller effect +0.0042 logit,
   OPPOSITE sign to the scan (-0.0081 per day), p 0.0289, within-cohort q
   0.0605 (fail), pooled q 0.0313 (pass). This is the only row where the
   independent method points the other way.

The remaining 23 rows pass the within-cohort propeller bar. The worm row
ARSC passes clearly: propeller effect +0.320 logit, p 0.00078,
within-cohort q 0.0091, pooled q 0.0012.

Note: the pooled-26 q is a second correction over the already-selected 26
rows, not a fresh discovery scan. It is reported alongside the within-cohort
q, not instead of it.

## 4. Leave-one-pool-out (LODO)

Method: `scripts/run_mature_fly_worm_lodo.py`. For each of the 26 rows, the
exact scan model (CLR with 0.5 pseudocount, OLS on age plus sex, sex dropped
when rank-deficient, slopes per day) was refit with each pool removed in
turn. This mirrors the mammalian LODO: the CLR table is computed once on the
full stratum and restricted to the kept pools; the eligible cell-type set is
not recomputed. Fly: 25 rows times 36 pools. Worm: 1 row times 9 pools.
Output: `results/followup-mature-fly-worm-lodo.csv`.

Verification: every full fit reproduces the committed scan effect to a max
difference of 9.5e-17. One row (fly body, 16-cell germline cyst) was
hand-recomputed from scratch with numpy least squares on a CLR table computed
independently from the raw counts: full fit and one leave-one-out fit agree
with the pipeline to all printed digits.

Results: the sign holds in every refit for 26 of 26 rows. Every leave-one-out
refit stays nominally significant: the largest LODO p-value across all 26
rows times all pools is 0.0132 (ARSC); the largest fly LODO p-value is
0.0067 (neuron 1). No single pool drives any row. The largest relative move
is the ocellus retinula cell, max slope change 0.0044 on a scan slope of
-0.0214 (about 21 percent), with sign and significance intact. ARSC moves at
most 0.0212 on a scan slope of 0.2275 (about 9 percent); the driving pool is
N2D6, the day-6 pool with the most cells (13,166).

LODO ranges (scan slope, min, max, max absolute change, driving pool, max p):

| Organ | Cell type | Scan | Min | Max | Max change | Driver pool | Max p |
|---|---|---|---|---|---|---|---|
| body | 16-cell germline cyst in germarium region 2a and 2b | +0.0243 | +0.0200 | +0.0269 | 0.0042 | AFCA_female_body_30_S3 | 0.00063 |
| body | adult Malpighian tubule principal cell of lower ureter | -0.0255 | -0.0282 | -0.0229 | 0.0027 | AFCA_female_body_30_S6 | 0.00019 |
| body | adult differentiating enterocyte | +0.0249 | +0.0219 | +0.0284 | 0.0035 | AFCA_female_body_30_S3 | 0.0014 |
| body | adult fat body_body | -0.0185 | -0.0194 | -0.0178 | 0.0009 | AFCA_male_body_30_S1 | 6.2e-09 |
| body | adult ventral nervous system | -0.0081 | -0.0091 | -0.0074 | 0.0010 | AFCA_female_body_30_S2 | 0.0045 |
| body | ejaculatory bulb | -0.0161 | -0.0171 | -0.0144 | 0.0017 | AFCA_male_body_70_S5 | 0.00023 |
| body | enteroendocrine cell | +0.0378 | +0.0323 | +0.0433 | 0.0055 | AFCA_female_body_30_S6 | 0.0024 |
| body | escort cell | +0.0121 | +0.0098 | +0.0137 | 0.0023 | AFCA_female_body_30_S2 | 0.0036 |
| body | follicle stem cell and prefollicle cell | +0.0142 | +0.0130 | +0.0151 | 0.0012 | AFCA_female_body_30_S2 | 0.000048 |
| body | hemocyte_body | -0.0110 | -0.0124 | -0.0099 | 0.0014 | AFCA_female_body_30_S6 | 0.00068 |
| body | male accessory gland main cell | -0.0134 | -0.0143 | -0.0113 | 0.0021 | AFCA_male_body_70_S6 | 0.0050 |
| body | posterior midgut* | -0.0340 | -0.0390 | -0.0308 | 0.0049 | AFCA_male_body_30_S2 | 0.0020 |
| body | prefollicle cell/stalk follicle cell | +0.0179 | +0.0166 | +0.0191 | 0.0014 | AFCA_female_body_30_S1 | 0.0000098 |
| body | stalk follicle cell | +0.0219 | +0.0190 | +0.0234 | 0.0029 | AFCA_female_body_30_S2 | 0.000026 |
| head | T neuron T3 | +0.0113 | +0.0097 | +0.0123 | 0.0017 | AFCA_male_head_30_S4 | 0.0011 |
| head | centrifugal neuron C2 | +0.0102 | +0.0073 | +0.0110 | 0.0028 | AFCA_male_head_30_S5 | 0.0062 |
| head | distal medullary amacrine neuron Dm8 | +0.0104 | +0.0089 | +0.0112 | 0.0015 | AFCA_male_head_70_S5 | 0.0010 |
| head | epithelial cell_head | +0.0137 | +0.0130 | +0.0147 | 0.0010 | AFCA_female_head_30_S1 | 1.6e-08 |
| head | neuron 1 | +0.0088 | +0.0076 | +0.0100 | 0.0012 | AFCA_female_head_30_S6 | 0.0067 |
| head | ocellus retinula cell | -0.0214 | -0.0249 | -0.0170 | 0.0044 | AFCA_female_head_70_S1 | 0.0041 |
| head | olfactory receptor neuron | -0.0063 | -0.0069 | -0.0052 | 0.0011 | AFCA_male_head_70_S6 | 0.0015 |
| head | photoreceptor cell R7 | -0.0124 | -0.0133 | -0.0114 | 0.0010 | AFCA_male_head_30_S3 | 0.00013 |
| head | pigment cell | +0.0120 | +0.0108 | +0.0131 | 0.0012 | AFCA_male_head_30_S6 | 0.00016 |
| head | skeletal muscle of head | +0.0113 | +0.0103 | +0.0120 | 0.0009 | AFCA_female_head_70_S4 | 4.1e-08 |
| head | transmedullary neuron Tm29 | -0.0192 | -0.0205 | -0.0173 | 0.0018 | AFCA_female_head_70_S3 | 0.0016 |
| Whole animal | ARSC | +0.2275 | +0.2063 | +0.2471 | 0.0212 | N2D6 | 0.0132 |

## 5. Analogs (not the study-aware check)

### 5a. Worm platform-aware (ARSC only)

The day 6 and older worm data span two sequencing platforms, and platform is
balanced within each age (1 pool on platform 0 and 2 pools on platform 1 at
each of days 6, 12, 14), so platform is not confounded with age. All refits
use the exact scan model plus platform where stated.

- Platform-adjusted (9 pools, age plus platform): effect +0.2275 per day,
  p 0.0051, identical to the scan effect to all printed digits. Platform
  adds nothing; the age coefficient is unchanged because platform is balanced
  within age.
- Platform 1 only (6 pools): effect +0.2033, p 0.0023.
- Platform 0 only (3 pools): effect +0.2758, p 0.399, standard error 0.20,
  1 residual degree of freedom. Null only from low power; same sign and
  similar size. FRAGILE: 3 pools, 1 residual degree of freedom.

Verdict: platform does not explain ARSC. Both platforms agree in sign and
size.

### 5b. Fly sex-stratified (25 rows)

Fly pools are sexed (6 male and 6 female pools per age per stratum), so each
row was refit within male pools only (18 pools) and within female pools only
(18 pools), using the exact scan model (sex drops out inside one sex).

- Sign matches the scan in 48 of 50 refits.
- 23 of 25 rows hold the scan sign in both sexes.
- The 2 splits: `adult ventral nervous system` (male-only -0.0168, p 6e-06;
  female-only +0.0006, p 0.835, null) and `posterior midgut*` (male-only
  +0.0001, p 0.993, null; female-only -0.0681, p 3e-06). Both mismatches are
  nulls in one sex, not significant reversals. The ventral nervous system row
  is male-driven; the posterior midgut row is female-driven.

## 6. Per-row verdicts

Propeller verdict uses the within-cohort q (the bar matching the mammalian
precedent); the pooled-26 q is shown alongside. LODO verdict is sign-hold
across all refits. Analog verdict is both-sexes sign-hold (fly) or
platform-adjusted agreement (worm).

| Organ | Cell type | Scan q | Prop q within | Prop q pooled-26 | LODO | Analog | Overall |
|---|---|---|---|---|---|---|---|
| body | 16-cell germline cyst in germarium region 2a and 2b | 0.0086 | 0.00061 | 0.00024 | holds | both sexes | PASS |
| body | adult Malpighian tubule principal cell of lower ureter | 0.0042 | 0.0106 | 0.0044 | holds | both sexes | PASS |
| body | adult differentiating enterocyte | 0.0142 | 0.00096 | 0.00044 | holds | both sexes | PASS |
| body | adult fat body_body | 0.0000010 | 0.0000010 | 0.00000030 | holds | both sexes | PASS |
| body | adult ventral nervous system | 0.0323 | 0.0605 FAIL | 0.0313 | holds | male-only | FAIL |
| body | ejaculatory bulb | 0.0042 | 0.417 FAIL | 0.313 | holds | both sexes | FAIL |
| body | enteroendocrine cell | 0.0204 | 0.0010 | 0.00053 | holds | both sexes | PASS |
| body | escort cell | 0.0293 | 0.0017 | 0.00077 | holds | both sexes | PASS |
| body | follicle stem cell and prefollicle cell | 0.0016 | 0.00000081 | 0.000000098 | holds | both sexes | PASS |
| body | hemocyte_body | 0.0086 | 0.974 FAIL | 0.960 | holds | both sexes | FAIL |
| body | male accessory gland main cell | 0.0336 | 0.0239 | 0.0105 | holds | both sexes | PASS |
| body | posterior midgut* | 0.0175 | 0.00061 | 0.00024 | holds | female-only | FAIL |
| body | prefollicle cell/stalk follicle cell | 0.00048 | 0.0000010 | 0.00000030 | holds | both sexes | PASS |
| body | stalk follicle cell | 0.0012 | 0.0000036 | 0.0000011 | holds | both sexes | PASS |
| head | T neuron T3 | 0.0124 | 0.0055 | 0.00089 | holds | both sexes | PASS |
| head | centrifugal neuron C2 | 0.0418 | 0.0168 | 0.0036 | holds | both sexes | PASS |
| head | distal medullary amacrine neuron Dm8 | 0.0124 | 0.0083 | 0.0016 | holds | both sexes | PASS |
| head | epithelial cell_head | 0.0000013 | 0.00000035 | 0.000000098 | holds | both sexes | PASS |
| head | neuron 1 | 0.0417 | 0.0232 | 0.0049 | holds | both sexes | PASS |
| head | ocellus retinula cell | 0.0320 | 0.0184 | 0.0040 | holds | both sexes | PASS |
| head | olfactory receptor neuron | 0.0145 | 0.0055 | 0.00089 | holds | both sexes | PASS |
| head | photoreceptor cell R7 | 0.0040 | 0.00017 | 0.000025 | holds | both sexes | PASS |
| head | pigment cell | 0.0042 | 0.00017 | 0.000025 | holds | both sexes | PASS |
| head | skeletal muscle of head | 0.0000023 | 0.000000040 | 0.000000012 | holds | both sexes | PASS |
| head | transmedullary neuron Tm29 | 0.0145 | 0.0055 | 0.00075 | holds | both sexes | PASS |
| Whole animal | ARSC | 0.0417 | 0.0091 | 0.0012 | holds | platform agrees | PASS |

Overall: 22 of 26 rows pass every check run (21 fly plus ARSC). Four fly
rows carry at least one failed check: `hemocyte_body` (propeller finds
nothing), `ejaculatory bulb` (propeller null), `adult ventral nervous
system` (propeller points the other way, within-cohort q fails, male-only),
`posterior midgut*` (female-only). No row fails LODO.

## 7. What this does and does not show

- The fly and worm passes are not single-pool artifacts: LODO holds for all
  26 rows, and every leave-one-out refit stays nominally significant.
- The independent method (propeller) confirms 23 of 26 rows on the
  within-cohort bar and 24 of 26 on the pooled-26 bar. It finds nothing for
  hemocyte_body and ejaculatory bulb, and points the other way for adult
  ventral nervous system.
- The sex splits mean two fly rows are single-sex driven in the pooled-sex
  fit. They still pass the scan and LODO; the split is disclosed, not hidden.
- The study-aware check could not be run for either species, so a
  study-driven artifact cannot be excluded the way it was for lung. The
  analogs (platform for worm, sex for fly) are the closest available
  substitutes and both support the scan rows they cover.
- No standalone total-cell refit was run (it was not in the ordered steps).
  The propeller design already adjusts for log total cells per pool, so the
  total-count confound is covered inside the propeller check.
- The same honest limit the paper already carries applies here: these are
  shares of captured cells within an organ, not shares of true cells in the
  tissue. If old tissue releases some types worse than others when the tissue
  is dissolved, the share moves without any true change. The checks above test
  the measurable versions of that problem (total cells, platform, sex, single
  pool, independent method). They cannot rule out a capture bias that changes
  with age and hits one cell type only.

## 8. Files and reproduction

Branch `followup-mature-fly-worm-checks`, from master `db40edc`. Nothing on
master was touched. No data was downloaded. All intermediates are small
CSVs; no large files were written to the old disk.

New scripts:

- `scripts/prepare_mature_fly_worm_propeller.py`: builds
  `data/interim/mature-propeller-fly-worm/` (3 cohorts).
- `scripts/run_mature_fly_worm_lodo.py`: leave-one-pool-out for the 26 rows,
  with the hand recomputation.
- `scripts/run_mature_fly_worm_analogs.py`: study-aware evidence, worm
  platform analog, fly sex analog.
- `scripts/build_mature_fly_worm_summary.py`: joins all checks into the
  summary CSV and writes the pooled-26 propeller q back into the propeller
  CSV.

New results:

- `results/followup-mature-fly-worm-propeller.csv` (178 rows; 90 fly head,
  72 fly body, 16 worm).
- `results/followup-mature-fly-worm-lodo.csv` (26 rows).
- `results/followup-mature-fly-worm-analogs.csv` (56 rows: 3 study-evidence,
  50 fly sex refits, 3 worm platform refits).
- `results/followup-mature-fly-worm-checks-summary.csv` (26 rows, one per
  passing row, all checks joined).

Reproduction (on Beast, agescape conda env, R with limma and speckle):

1. `python scripts/prepare_mature_fly_worm_propeller.py`
2. `Rscript scripts/run_reviewer_phase1_propeller.R data/interim/mature-propeller-fly-worm results/followup-mature-fly-worm-propeller.csv`
3. `python scripts/run_mature_fly_worm_lodo.py`
4. `python scripts/run_mature_fly_worm_analogs.py`
5. `python scripts/build_mature_fly_worm_summary.py`
