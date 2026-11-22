# Mature 30-plus full-depth drift, entropy, and identity

This is Phase 2 Item 9a. All derivation and tested contrasts use donors aged at least 30 years. The inherited YOUNG label means only the 30-to-39 cross-fitted baseline, with leave-one-out references for those donors. No young-derived P1 or original E1/E2 set enters this analysis.

## Frozen design

- Primary outcomes use full-depth raw counts as within-cell proportions; at least 1500 UMI is an eligibility filter only.
- Rarefied 1500 and 3000 UMI outcomes are labeled sensitivities. Every downsample output has explicit zeros eliminated before stored-entry counting.
- Donor is the unit for every p-value. Slopes are per decade. Full-depth primaries carry five-seed split-half SNR and MDE80.
- Per-cell-type BH families use the frozen dataset, part, metric, gene-set, depth, and model slices. Pooled rows are singleton families; stratified-global primary q is NaN by design.
- E2_stratum is primary E2 evidence. E2_global is only a depth-vulnerable sensitivity.

## Donor and cell restriction

Original selected donors include 117 Gao and 981 Yazar donors. Original global models contain 113 Gao donors and 522,761 eligible nuclei, and 973 Yazar donors and 617,760 eligible cells.

| Dataset | Selected donors original -> mature | Lost | Baseline donors 30-39 | Global donors original -> mature | Eligible cells original -> mature | Cells lost |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Gao heart | 117 -> 109 | 8 | 12 | 113 -> 104 (center_shift) / 105 (scatter/Part 2) | 522,761 -> 482,902 | 39,859 |
| Yazar blood | 981 -> 932 | 49 | 59 | 973 -> 925 | 617,760 -> 592,250 | 25,510 |

- Gao heart: 482,902 total eligible cells at >= 1500 UMI across 109 mature donors; 482,863 cells in pooled primary model; 477,460 cells across 104 donors in stratified-global center_shift (480,355 cells across 105 donors in stratified-global scatter).
- Yazar blood: 592,250 total eligible cells at >= 1500 UMI across 932 mature donors; 592,250 cells in pooled primary model; 503,423 cells across 925 donors in stratified-global center_shift (509,717 cells across 925 donors in stratified-global scatter).

## Primary verdicts and exact effects

Effects below are stratified-global SD per decade. Original values are full-range context only.

| Dataset | Endpoint | Original slope, p | Mature slope, p, SNR, MDE80 | Mature verdict | Reason |
| --- | --- | --- | --- | --- | --- |
| Gao heart | Part 1 center shift | 0.01977, p 2.778e-01 | -0.01298, p 5.945e-01, SNR 0.869, MDE80 0.0687 | DOES NOT PROVE | Global center_shift fails PROVES (slope=-0.0130, p=0.594, MDE80=0.069 [high_power<=0.25]) |
| Gao heart | Part 1 scatter | -0.0632, p 2.697e-02 | -0.07228, p 4.103e-02, SNR 4.85, MDE80 0.0988 | DOES NOT PROVE | Global scatter fails PROVES (slope=-0.0723, p=0.041, MDE80=0.099 [high_power<=0.25]) |
| Gao heart | Part 2 Shannon | -0.09984, p 1.281e-02 | -0.1001, p 5.390e-02, SNR 7.05, MDE80 0.145 | DOES NOT PROVE | shannon fails PROVES (slope=-0.1001, p=0.0539, MDE80=0.145 [high_power<=0.25]) |
| Gao heart | Part 2 KL methods | 0.09984, p 1.281e-02 | 0.1001, p 5.390e-02, SNR 7.05, MDE80 0.145 | DOES NOT PROVE | deman_kl_methods fails PROVES (slope=0.1001, p=0.0539, MDE80=0.145 [high_power<=0.25]) |
| Gao heart | Part 2 KL sqrt | 0.03824, p 1.828e-02 | 0.04171, p 5.436e-02, SNR 2.77, MDE80 0.0606 | DOES NOT PROVE | deman_kl_sqrt fails PROVES (slope=0.0417, p=0.0544, MDE80=0.061 [high_power<=0.25]) |
| Gao heart | Part 3 E1 count | 0.03071, p 5.643e-01 | -0.06331, p 4.683e-01, SNR 2.16, MDE80 0.246 | DOES NOT PROVE | E1 fails randomness PROVES (cnt p=0.468, SNR=2.16, MDE80=0.246<=0.5) |
| Gao heart | Part 3 E2 stratum count | -0.0179, p 6.581e-01 | -0.1269, p 3.112e-02, SNR 36.3, MDE80 0.164 | DOES NOT PROVE | E2 fails randomness PROVES (cnt p=0.0311, SNR=36.33, MDE80=0.164<=0.5) |
| Yazar blood | Part 1 center shift | -0.002807, p 3.839e-01 | -0.004388, p 2.659e-01, SNR 0.731, MDE80 0.0111 | DOES NOT PROVE | Global center_shift fails PROVES (slope=-0.0044, p=0.266, MDE80=0.011 [high_power<=0.25]) |
| Yazar blood | Part 1 scatter | -0.01689, p 1.481e-04 | -0.01967, p 6.685e-05, SNR 3.31, MDE80 0.0138 | DOES NOT PROVE | Global scatter fails PROVES (slope=-0.0197, p=6.69e-05, MDE80=0.014 [high_power<=0.25]) |
| Yazar blood | Part 2 Shannon | -0.02219, p 7.882e-03 | -0.02977, p 3.649e-03, SNR 5.04, MDE80 0.0287 | DOES NOT PROVE | shannon fails PROVES (slope=-0.0298, p=0.00365, MDE80=0.029 [high_power<=0.25]) |
| Yazar blood | Part 2 KL methods | 0.02219, p 7.882e-03 | 0.02977, p 3.649e-03, SNR 5.04, MDE80 0.0287 | DOES NOT PROVE | deman_kl_methods fails PROVES (slope=0.0298, p=0.00365, MDE80=0.029 [high_power<=0.25]) |
| Yazar blood | Part 2 KL sqrt | 0.01991, p 1.064e-04 | 0.02333, p 1.399e-04, SNR 3.81, MDE80 0.0171 | DOES NOT PROVE | deman_kl_sqrt fails PROVES (slope=0.0233, p=0.00014, MDE80=0.017 [high_power<=0.25]) |
| Yazar blood | Part 3 E1 count | -0.01207, p 8.246e-02 | -0.01868, p 7.301e-02, SNR 1.54, MDE80 0.0292 | DOES NOT PROVE | E1 fails randomness PROVES (cnt p=0.073, SNR=1.54, MDE80=0.029<=0.5) |
| Yazar blood | Part 3 E2 stratum count | -0.006566, p 4.564e-01 | -0.03269, p 3.086e-03, SNR 7.87, MDE80 0.0309 | DOES NOT PROVE | E2 fails randomness PROVES (cnt p=0.00309, SNR=7.87, MDE80=0.031<=0.5) |

## What survived or died

### Gao heart

- Part 1 center shift: The original non-significant result remained non-significant; original slope 0.01977 (p 2.778e-01), mature slope -0.01298 (p 5.945e-01, MDE80 0.0687).
- Part 1 scatter: The nominally significant direction survived; original slope -0.0632 (p 2.697e-02), mature slope -0.07228 (p 4.103e-02, MDE80 0.0988).
- Part 2 Shannon: The original nominal significance died; original slope -0.09984 (p 1.281e-02), mature slope -0.1001 (p 5.390e-02, MDE80 0.145).
- Part 2 KL methods: The original nominal significance died; original slope 0.09984 (p 1.281e-02), mature slope 0.1001 (p 5.390e-02, MDE80 0.145).
- Part 2 KL sqrt: The original nominal significance died; original slope 0.03824 (p 1.828e-02), mature slope 0.04171 (p 5.436e-02, MDE80 0.0606).
- Part 3 E1 count: The original non-significant result remained non-significant; original slope 0.03071 (p 5.643e-01), mature slope -0.06331 (p 4.683e-01, MDE80 0.246).
- Part 3 E2 stratum count: A new nominally significant association emerged; original slope -0.0179 (p 6.581e-01), mature slope -0.1269 (p 3.112e-02, MDE80 0.164).

### Yazar blood

- Part 1 center shift: The original non-significant result remained non-significant; original slope -0.002807 (p 3.839e-01), mature slope -0.004388 (p 2.659e-01, MDE80 0.0111).
- Part 1 scatter: The nominally significant direction survived; original slope -0.01689 (p 1.481e-04), mature slope -0.01967 (p 6.685e-05, MDE80 0.0138).
- Part 2 Shannon: The nominally significant direction survived; original slope -0.02219 (p 7.882e-03), mature slope -0.02977 (p 3.649e-03, MDE80 0.0287).
- Part 2 KL methods: The nominally significant direction survived; original slope 0.02219 (p 7.882e-03), mature slope 0.02977 (p 3.649e-03, MDE80 0.0287).
- Part 2 KL sqrt: The nominally significant direction survived; original slope 0.01991 (p 1.064e-04), mature slope 0.02333 (p 1.399e-04, MDE80 0.0171).
- Part 3 E1 count: The original non-significant result remained non-significant; original slope -0.01207 (p 8.246e-02), mature slope -0.01868 (p 7.301e-02, MDE80 0.0292).
- Part 3 E2 stratum count: A new nominally significant association emerged; original slope -0.006566 (p 4.564e-01), mature slope -0.03269 (p 3.086e-03, MDE80 0.0309).

Direction note on Part 3 E2: the two nominally significant mature E2 stratum count slopes are negative in both cohorts (Gao -0.1269 SD per decade, Yazar -0.0327 SD per decade), so ectopic counts fall with age rather than rise. That is the opposite direction from the randomness prediction, which is why both E2 verdicts stay DOES NOT PROVE despite the small p-values. The same holds for the E1 count slopes, which are negative in both cohorts.

## Pooled rows, reported separately

Pooled slopes are not within-cell-type evidence. A pooled positive with a discordant stratified slope is read as a composition effect.

| Dataset | Part and metric | Pooled raw slope | Pooled SD slope | p | Stratified SD slope | Stratified p |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Gao heart | Part 1 center_shift | -0.00205683 | -0.0416095 | 2.692e-01 | -0.012977 | 5.945e-01 |
| Gao heart | Part 1 scatter | -0.00504373 | -0.101133 | 3.466e-02 | -0.0722825 | 4.103e-02 |
| Gao heart | Part 2 shannon | -0.0216372 | -0.102704 | 1.206e-01 | -0.100142 | 5.390e-02 |
| Gao heart | Part 2 deman_kl_methods | 0.0216372 | 0.102704 | 1.206e-01 | 0.100142 | 5.390e-02 |
| Gao heart | Part 2 deman_kl_sqrt | 0.00928506 | 0.054861 | 1.962e-02 | 0.0417144 | 5.436e-02 |
| Yazar blood | Part 1 center_shift | 0.000455084 | 0.0400032 | 2.992e-05 | -0.00438789 | 2.659e-01 |
| Yazar blood | Part 1 scatter | 0.000171941 | 0.00940395 | 2.925e-01 | -0.0196652 | 6.685e-05 |
| Yazar blood | Part 2 shannon | 0.00292141 | 0.0261701 | 6.436e-02 | -0.0297734 | 3.649e-03 |
| Yazar blood | Part 2 deman_kl_methods | -0.00292141 | -0.0261701 | 6.436e-02 | 0.0297734 | 3.649e-03 |
| Yazar blood | Part 2 deman_kl_sqrt | 0.00066344 | 0.0119144 | 1.906e-01 | 0.0233329 | 1.399e-04 |

### Yazar naive CD8 exception retest

Full-range context only (never mature evidence): SD slope 0.0465764, p 6.948e-07, q 5.558e-06, donors 384, cells 21,073. Mature 30-plus retest: raw slope 0.00100582, SD slope 0.0498023, p 5.920e-05, q 4.736e-04, SNR 19, MDE80 0.0343, donors 344, cells 17,878.

## Sensitivity outcomes

| Dataset | Sensitivity | SD slope | p |
| --- | --- | ---: | ---: |
| Gao heart | Center equal reference size | 0.0377232 | 7.656e-02 |
| Gao heart | Center Pythagorean centroid | 0.0204142 | 6.315e-01 |
| Gao heart | Center spline residual | -0.0190183 | 7.286e-01 |
| Gao heart | Center rarefied 1500 | -0.0479176 | 2.553e-01 |
| Gao heart | Scatter rarefied 1500 | -0.140079 | 1.572e-02 |
| Gao heart | Shannon rarefied 1500 | -0.111861 | 6.479e-02 |
| Gao heart | Shannon rarefied 3000 | -0.101854 | 7.270e-02 |
| Gao heart | E1 count rarefied 1500 | -0.0598574 | 4.491e-01 |
| Gao heart | E1 count spline residual | -0.0778737 | 4.119e-01 |
| Gao heart | E1 stringent count | -0.063627 | 4.660e-01 |
| Gao heart | E2 global count, depth vulnerable | -0.0883833 | 2.752e-02 |
| Gao heart | E2 global count rarefied 1500 | -0.161039 | 6.152e-03 |
| Yazar blood | Center equal reference size | 0.00395865 | 3.461e-01 |
| Yazar blood | Center Pythagorean centroid | 0.013539 | 2.095e-01 |
| Yazar blood | Center spline residual | -0.026085 | 2.204e-02 |
| Yazar blood | Center rarefied 1500 | -0.0220972 | 2.078e-03 |
| Yazar blood | Scatter rarefied 1500 | -0.0320917 | 1.361e-05 |
| Yazar blood | Shannon rarefied 1500 | -0.0269001 | 7.975e-03 |
| Yazar blood | Shannon rarefied 3000 | -0.026752 | 2.811e-01 |
| Yazar blood | E1 count rarefied 1500 | 0.00518537 | 5.393e-01 |
| Yazar blood | E1 count spline residual | -0.0188946 | 8.750e-02 |
| Yazar blood | E1 stringent count | -0.0224149 | 2.784e-02 |
| Yazar blood | E2 global count, depth vulnerable | -0.0275011 | 3.267e-03 |
| Yazar blood | E2 global count rarefied 1500 | -0.0368215 | 7.228e-04 |

## Part 3 count gate and secondary metrics

| Dataset | Set and metric | Depth | SD slope | p | SNR | MDE80 |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| Gao heart | E1 singleton_fraction_20cells | subsampled_20cells | -0.00533377 | 9.174e-01 | NA | 0.145 |
| Gao heart | E1 jaccard_both_pos | full_depth | -0.0910174 | 6.064e-02 | 6.38 | 0.136 |
| Gao heart | E1 excess_jaccard | full_depth | -0.0592382 | 1.098e-01 | NA | 0.104 |
| Gao heart | E2_stratum singleton_fraction_20cells | bounded_depth_strata_20cells | 0.150908 | 1.160e-02 | NA | 0.166 |
| Gao heart | E2_stratum jaccard_both_pos | bounded_depth_strata | -0.0399652 | 5.279e-01 | 14 | 0.179 |
| Gao heart | E2_stratum excess_jaccard | bounded_depth_strata | 0.00629028 | 8.916e-01 | NA | 0.13 |
| Yazar blood | E1 singleton_fraction_20cells | subsampled_20cells | 0.0152983 | 1.615e-01 | NA | 0.0306 |
| Yazar blood | E1 jaccard_both_pos | full_depth | -0.00816447 | 5.135e-01 | 1.15 | 0.035 |
| Yazar blood | E1 excess_jaccard | full_depth | -0.00969114 | 4.364e-01 | NA | 0.0349 |
| Yazar blood | E2_stratum singleton_fraction_20cells | bounded_depth_strata_20cells | -0.00940226 | 4.575e-01 | NA | 0.0355 |
| Yazar blood | E2_stratum jaccard_both_pos | bounded_depth_strata | 0.00214723 | 8.803e-01 | 0.404 | 0.04 |
| Yazar blood | E2_stratum excess_jaccard | bounded_depth_strata | 0.00101252 | 9.302e-01 | NA | 0.0324 |

## Mature E1/E2 set sizes versus originals

| Dataset | E1 cell types original -> mature | E1 median markers original -> mature | E2 global cell types original -> mature | E2 global median genes original -> mature | Mature E2 stratum CT-strata, median genes | Genes kept |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Gao heart | 12 -> 11 | 13.5 -> 7 | NA -> 11 | 4677 -> 4523 | 20, 4843 | 15,930 |
| Yazar blood | 14 -> 13 | 7.5 -> 5 | NA -> 13 | 2.556e+04 -> 2.554e+04 | 21, 2.592e+04 | 35,307 |

Original E2 global cell-type counts and E2 stratum CT-stratum counts were not recorded in the original QC JSON, so those original fields are NA rather than inferred.

## Power after restriction

| Dataset | Endpoint | MDE80 original -> mature | Fold change |
| --- | --- | ---: | ---: |
| Gao heart | Part 1 center shift | 0.05123 -> 0.06873 | 1.34x |
| Gao heart | Part 1 scatter | 0.07969 -> 0.0988 | 1.24x |
| Gao heart | Part 2 Shannon | 0.1115 -> 0.1453 | 1.3x |
| Gao heart | Part 2 KL methods | 0.1115 -> 0.1453 | 1.3x |
| Gao heart | Part 2 KL sqrt | 0.04513 -> 0.06063 | 1.34x |
| Gao heart | Part 3 E1 count | 0.1501 -> 0.246 | 1.64x |
| Gao heart | Part 3 E2 stratum count | 0.114 -> 0.1642 | 1.44x |
| Yazar blood | Part 1 center shift | 0.009037 -> 0.01106 | 1.22x |
| Yazar blood | Part 1 scatter | 0.01244 -> 0.01377 | 1.11x |
| Yazar blood | Part 2 Shannon | 0.02337 -> 0.02865 | 1.23x |
| Yazar blood | Part 2 KL methods | 0.02337 -> 0.02865 | 1.23x |
| Yazar blood | Part 2 KL sqrt | 0.01434 -> 0.01711 | 1.19x |
| Yazar blood | Part 3 E1 count | 0.01947 -> 0.02919 | 1.5x |
| Yazar blood | Part 3 E2 stratum count | 0.02471 -> 0.0309 | 1.25x |

## Robustness, verification, and items not verified

- Rarefaction robustness: 150 of 159 rows pass, versus 152 of 159 originally.
- Independent BH recomputation matched all 7,603 eligible per-cell-type p-values across all six output CSVs. Pooled singleton (m = 1, q = p) and stratified-global NaN rules also passed.
- Regression suites: tests/test_followup_mature_fulldepth.py (10 passed, including explicit-zero rarefaction cleanup and under-30 leakage checks) and tests/test_followup_fulldepth.py (9 passed). Real metadata checks confirmed minimum mapped donor age 30.0 in both Gao and Yazar, and all reference donors strictly within 30 to 39 years.
- Explicit-zero rarefaction check: verified on both Pass 1 caches (y_1500_det strictly smaller than y_full_det; exact 1500 UMI per cell) and Part 3 rarefied_1500 outputs (distinct from full_depth outputs).
- Not verified / intentionally omitted: (1) Cross-check against original Part A donor CSVs was not performed because original Part A references and scored donors are full-range (20-plus / 19-plus). (2) Original E2 global cell-type count and original E2 stratum CT-stratum count were not stored in results/followup-fulldepth-qc.json and are reported as NA. (3) Gao Parts 1 and 2 and Yazar Pass 1 shards/pickles were produced by the prior kimi-code/k3.1-pro session and verified via shard byte-size, cell-count, donor-age, explicit-zero, and BH checks rather than recomputed from HDF5.

## Runtime and reproduction

- Current run model: Gemini 4 Neon (XHigh) under Antigravity.
- Prior cached Gao Parts 1-2 and Yazar Pass 1 model: kimi-code/k3.1-pro.
- Prior cached Pass 1/2 durations (verified from checkpoints): Gao Pass 1 565.92 s, Gao Pass 2 506.41 s, Yazar Pass 1 314.95 s.

| Run | Model | Start UTC | End UTC | Elapsed seconds |
| --- | --- | --- | --- | ---: |
| yazar parts12 | Gemini 4 Neon (XHigh) under Antigravity | 2026-11-21T21:22:46.156773+00:00 | 2026-11-21T21:32:59.757265+00:00 | 613.6 |
| all part3 | Gemini 4 Neon (XHigh) under Antigravity | 2026-11-21T21:33:33.561469+00:00 | 2026-11-21T22:07:17.745613+00:00 | 2024.184 |

Pass-level timings are in the pipeline logs and cache checkpoints. Commands:

```
PY=/home/reza/miniconda3/envs/agescape/bin/python
$PY scripts/resource_guard.py run --path /home/reza/projects/agescape-revision-data --min-free-gb 30 --memory-gb 12 -- $PY scripts/run_followup_mature_fulldepth.py --dataset yazar --stage parts12 --reuse-pass1
$PY scripts/run_followup_mature_fulldepth.py --dataset all --stage part3 --reuse-pass1
$PY -m pytest -p no:cacheprovider tests/test_followup_mature_fulldepth.py tests/test_followup_fulldepth.py -v
$PY scripts/report_followup_mature_fulldepth.py
```
