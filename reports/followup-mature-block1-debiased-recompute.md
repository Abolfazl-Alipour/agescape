# Independent recompute: mature naive CD4 debiased slope

Method: fresh numpy QR OLS plus scipy two-sided t tests, applied to the
committed cache tables only (donor_group_stats.csv, finite_n_calib.csv,
donors_meta.csv). No report text was used as input except the raw-Y +N+D
pair, which lives in no diagnostics key. The design replicates the
singleton-pool prune and cross-checks the 73 pool dummies against the
committed primary row (777 donors, residual df 701).

## Calibration bias (subsampled minus full-N, high-N donors)

| n_sub | donor-seed rows | recomputed bias | diagnostics bias | abs diff |
| --- | --- | --- | --- | --- |
| 20 | 635 | 0.008139 | 0.008139 | 1.73e-18 |
| 50 | 635 | 0.011184 | 0.011184 | 1.73e-18 |
| 100 | 635 | 0.007137 | 0.007137 | 1.73e-18 |
| 200 | 635 | 0.001637 | 0.001637 | 0.00e+00 |

Max calibration bias difference: 1.73e-18.

## Slopes

| arm | recomputed slope | recomputed p | expected slope | expected p | abs slope diff |
| --- | --- | --- | --- | --- | --- |
| primary | 0.000714 | 5.174e-08 | 0.000714 | 5.174e-08 | 1.10e-17 |
| induced_by_calibration | 0.000401 | 8.396e-07 | 0.000401 | 8.396e-07 | 2.71e-19 |
| debiased | 0.000313 | 0.0006379 | 0.000313 | 0.0006379 | 1.08e-19 |
| debiased_plus_depth | 0.000295 | 0.001063 | 0.000295 | 0.001063 | 5.42e-20 |
| raw_plus_logN | 0.000286 | 0.006605 | 0.000286 | 0.006605 | 0.00e+00 |
| raw_plus_N_plus_depth | 0.000245 | 0.01893 | 0.000240 | 0.01893 | 4.64e-06 |
| debiased_plus_logN | 0.000238 | 0.01013 | new | new | new |
| debiased_plus_N_plus_depth | 0.000193 | 0.03253 | new | new | new |

Calibration share of the primary slope: 0.5613 (56.13 percent).
Max slope difference on checked arms: 1.10e-17.
Max p difference on checked arms: 3.144e-18.

## Verdict

VERIFIED: every checked recomputed slope and p matches its expected value (diagnostics arms to 1e-6 in slope; the report-sourced raw-Y +N+D arm to its 5-decimal rounding).
The debiased slope is 0.000313 per decade
(p = 0.0006379), still significant after calibration.
The exact calibration share is 56.13 percent, not 56.4 percent:
the report divided the induced slope by the rounded primary slope 0.00071
instead of the exact 0.000714. The tex wording (56 percent) is unaffected.
The true debiased +logN slope is 0.000238
(p = 0.01013) and the true debiased +N+D slope
is 0.000193
(p = 0.03253); the transcription CSV
labels the raw-Y fits 0.000286 and 0.00024 with those names instead.
