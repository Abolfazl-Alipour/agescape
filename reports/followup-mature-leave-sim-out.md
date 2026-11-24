# Leave-Sim-out Gao scatter refit

## Method

The fits use `results/followup-mature-fulldepth-gao-donor-part1.csv` and the stratified-global scatter specification in `scripts/run_followup_mature_fulldepth.py:798-969`, with eligibility reconstructed as in the Part 1 pipeline. The model standardizes scatter within cell type, centers sex and depth within cell type, and includes cell-type strata and study/pool fixed effects. Inference uses donor-clustered standard errors with t-based p-values, with donor as the unit of analysis.

The leave-Sim-out fit repeats the identical specification after dropping these three Sim 2021 donors:

- `Sim 2021:GSM4742854_Adult1_processed` (age 35, female)
- `Sim 2021:GSM4742855_Adult2_processed` (age 42, male)
- `Sim 2021:GSM4742856_Adult3_processed` (age 42, male)

## Results

| Fit | Slope per decade | P-value | Donors | Rows |
| --- | ---: | ---: | ---: | ---: |
| Primary | -0.07228251880607431 | 0.041028425872899024 | 105 | 829 |
| Leave-Sim-out | -0.09351698226875774 | 0.02998223250331266 | 102 | 803 |

The leave-Sim-out result rounds to slope -0.0935 per decade and p = 0.0300, so the Sim donors do not drive the adjusted slope.
