# Mature drift and mouse lung identity at 18-plus and 30-plus

## Scope

This is Phase 2 Item 4. Mouse mature models include only mice aged at least 18 months. Human mature models include only donors aged at least 30 years, with references built only from donors aged 30 through 39 years. No manuscript file was changed. The model used for all new pipeline runs was kimi-code/k3.1-pro.

## Models and reference construction

The Hellinger construction is unchanged: technical and sex-linked genes are removed, each cell count vector is converted to a probability vector, and distance to a reference probability vector is sqrt(1 - Bhattacharyya affinity). A donor by cell-type row requires at least 20 cells and at least 50 reference cells after removing that donor from the reference. The mouse reference contains only 18-month mice and is cross-fitted by donor.

Global drift models standardize the outcome within tissue, cohort, and cell-type unit; include unit intercepts, age, and sex; and cluster standard errors by donor. Mouse slopes are per 0.4 years and human slopes are per 12.25 years. Linear, inverse-square-root UMI, combined linear plus inverse-square-root, quadratic, and spline depth models retain the original definitions. Fixed-depth mouse references and scored cells are rarefied without replacement at 500, 1500, or 3000 UMIs with seed 42. Full-depth mouse global linear fits are reused rather than rerun.

Human references and scores were rebuilt from the seven local raw Census objects. Normal primary cells with at least 1500 raw UMIs were selected. Reference ages are 30 through 39 inclusive, not 20 through 40. The original eligible-unit universe is a mathematical superset of the newly eligible mature units, so restricting raw reads to that universe does not remove a newly eligible unit. Per-cohort models retain at least 6 donors per unit and at least 6 donor clusters globally.

The interaction uses the existing pooled OLS construction with unit fixed effects, species-specific sex and depth coefficients, and donor-cluster-robust inference. The human-minus-mouse interaction is reported on the published 0.4-year mouse and 12.25-year human scales. A 1000-replicate donor-cluster bootstrap resamples donors within species.

For mouse lung identity, off-program markers are derived only from 18-month cells with at least 1500 UMIs. Marker rules are unchanged: target detection at most 5 percent, another cell-type detection at least 20 percent, other-cell CP10K at least 20, and other-cell CP10K at least 5 times the quantity target CP10K plus 1. Marker derivation excludes the scored donor where applicable. New off-program fits cover 500, 1500, 3000, and full depth, with and without the actual-library covariate. Existing mature entropy results are reused without recalculation.

De Man controls are restricted to age 30-plus before rarefaction, nonzero-gene selection, cell balancing, and regression. The published square-root KL, full-depth Shannon, and 1500-UMI Shannon methods retain their original formulas and 20 seeds. Models adjust for age, sex, and log actual UMI.

## Restricted mouse drift

The correctly restricted full-depth frequentist result reuses 14 mice at ages 18, 21, 24, and 30 months: 4, 2, 4, and 4 mice respectively. The global fit contains 7 tissues, 51 eligible units, and 480 donor-unit rows.

| Endpoint | Model | Slope | SE | 95% interval | p | Status |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| Mean cell-to-reference Hellinger | Unadjusted | +0.109987 | 0.106069 | -0.119162 to +0.339135 | 0.318672 | REUSED |
| Mean cell-to-reference Hellinger | Linear log UMI | -0.013715 | 0.057853 | -0.138700 to +0.111270 | 0.816303 | REUSED |
| Mean cell-to-reference Hellinger | Spline log UMI | -0.012115 | 0.066194 | -0.155119 to +0.130890 | 0.857607 | NEW |
| Donor pseudobulk Hellinger | Unadjusted | +0.065764 | 0.082716 | -0.112933 to +0.244462 | 0.440854 | REUSED |
| Donor pseudobulk Hellinger | Linear log UMI | +0.035006 | 0.076443 | -0.130138 to +0.200150 | 0.654558 | REUSED |
| Donor pseudobulk Hellinger | Spline log UMI | +0.050985 | 0.076837 | -0.115011 to +0.216981 | 0.518571 | NEW |

The reused Bayesian depth-adjusted estimate is +0.038944 with 94 percent highest-density interval -0.099779 to +0.177334 and probability positive 0.7068. Its unadjusted prior-0.50 estimate is +0.177036 with 94 percent interval -0.094440 to +0.450339. Bayesian and frequentist point estimates need not have the same sign; their intervals remain compatible with continuation in either direction.

The historical 18-month-reference +0.007229 estimate, p = 0.665610, 95 percent interval -0.026995 to +0.041454, is unrestricted and mixed-age: it regresses all 23 mice and 754 donor-unit rows. It is not a mature estimate. The restricted depth-adjusted estimate is -0.013715 with a much wider interval and 14 mice. The historical unrestricted unadjusted slope was +0.042146, p = 0.037471; the restricted unadjusted slope is +0.109987, p = 0.318672. These are different estimands, not a correction of the historical number by relabeling.

### New mouse depth, per-organ, pseudobulk, and influence arms

| Depth | Linear depth-adjusted mean Hellinger slope | 95% interval | p | Donors | Units |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 500 | +0.065614 | -0.194870 to +0.326098 | 0.595526 | 14 | 51 |
| 1500 | +0.087384 | -0.193482 to +0.368250 | 0.513261 | 14 | 51 |
| 3000 | +0.073275 | -0.190450 to +0.337001 | 0.558657 | 14 | 46 |
| Full | -0.013715 | -0.138700 to +0.111270 | 0.816303 | 14 | 51 |

All fixed-depth mean Hellinger intervals cross zero under every fitted depth model. Pseudobulk slopes are negative at 500 UMIs, including -0.141990 with 95 percent interval -0.275187 to -0.008794 and p = 0.038445 under linear depth adjustment, but are imprecise at 1500, 3000, and full depth. This depth-sensitive negative pseudobulk result is not evidence for a general positive mature drift.

The full-depth lung per-organ slope is -0.193723, SE 0.069628, 95 percent interval -0.351234 to -0.036213, p = 0.021322, with 10 mice and 10 eligible units. Limb muscle is +0.285647, p = 0.042705, and kidney is -0.230895, p = 0.056719. Per-organ fits are decomposition and sensitivity arms, not a separately multiplicity-adjusted family.

Leave-one-donor-out full-depth linear mean-Hellinger slopes range from -0.057600 to +0.053734, with p values from 0.212079 to 0.863162. Full-depth unadjusted slopes range from -0.034974 to +0.176095. Full-depth pseudobulk linear slopes range from -0.019505 to +0.083515. Fixed-depth linear mean-Hellinger leave-one-out ranges are -0.052229 to +0.125543 at 500, -0.034206 to +0.169337 at 1500, and -0.015195 to +0.165773 at 3000. Reference endpoints are held fixed in these regression influence refits.

There are 36 explicit UNTESTABLE per-organ by depth rows. At full depth these are Bladder, Large Intestine, Liver, Mammary Gland, Pancreas, Skin, Thymus, Tongue, and Trachea, due to absent mature reference endpoints or failure of donor, unit, or rank floors. They are not null results. Complete per-organ and influence rows are in results/followup-mature-drift-mouse.csv.

## Mature human drift

The mature global fit contains 407 donors, 6 cohorts, 87 units, and 4657 donor-unit rows. ImmuneAging_BAL contributes no scored unit because only 1 donor aged 30-39 was available and no unit met the floor of 2 reference donors with at least 20 QC cells. This arm is UNTESTABLE, not null.

| Depth model | Slope | SE | 95% interval | p |
| --- | ---: | ---: | ---: | ---: |
| Unadjusted | +0.062550 | 0.035829 | -0.007883 to +0.132984 | 0.081601 |
| Linear log UMI | -0.007101 | 0.013933 | -0.034490 to +0.020289 | 0.610587 |
| Inverse sqrt UMI | +0.030186 | 0.015256 | +0.000195 to +0.060177 | 0.048538 |
| Linear plus inverse sqrt UMI | +0.009659 | 0.016782 | -0.023331 to +0.042649 | 0.565227 |
| Quadratic log UMI | -0.006376 | 0.013737 | -0.033381 to +0.020629 | 0.642785 |
| Spline log UMI | -0.006456 | 0.013664 | -0.033316 to +0.020405 | 0.636845 |

The historical mixed human estimate was +0.036483, p = 0.221842, with 95 percent interval -0.022121 to +0.095087, 479 donors, 135 units, and 5946 rows. The new unadjusted mature estimate is larger but its interval still crosses zero. The inverse-square-root depth model is nominally positive, whereas linear, quadratic, and spline estimates are slightly negative and imprecise. Depth-model disagreement remains explicit; no depth-adjusted positive continuation is claimed.

| Cohort | Scored donors | Eligible reference donors | Units | Unadjusted slope, p | Linear-depth slope, p |
| --- | ---: | ---: | ---: | --- | --- |
| MSSM_DLPFC | 208 | 8 | 19 | +0.126373, 0.009366 | +0.002346, 0.853756 |
| HBCC_DLPFC | 64 | 15 | 19 | +0.003632, 0.965249 | -0.002131, 0.907939 |
| Kidney_LivingDonor | 19 | 6 | 6 | -0.025317, 0.918606 | -0.019199, 0.878296 |
| Kidney_SCR_v15 | 18 | 6 | 9 | +0.367293, 0.250209 | +0.214614, 0.177030 |
| Lung_CellRef | 13 | 2 | 3 | +0.232355, 0.416077 | +0.010000, 0.908029 |
| Lung_Integrated | 85 | 22 | 31 | -0.067536, 0.434503 | -0.017977, 0.866462 |
| ImmuneAging_BAL | 0 | 0 eligible | 0 | UNTESTABLE | UNTESTABLE |

The MSSM unadjusted cohort slope loses significance under linear depth adjustment. Lung CellRef has only 2 reference donors and 3 eligible units, so its wide intervals are not stable evidence of a positive cohort effect. Per-tissue brain, kidney, and lung decompositions are retained in the CSV.

## Mature cross-species interaction

| Model | Human minus mouse interaction | SE | Analytic 95% interval | Analytic p | Bootstrap 95% interval | Bootstrap p |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Unadjusted | -0.047436 | 0.104372 | -0.252593 to +0.157721 | 0.649711 | -0.283902 to +0.265553 | 0.686 |
| Linear depth adjusted | +0.006614 | 0.055189 | -0.101867 to +0.115096 | 0.904660 | -0.135061 to +0.111142 | 0.978 |

Both fits use only 14 mature mice and 407 mature human donors, with 5137 donor-unit rows and 421 donor clusters. Both analytic and bootstrap interaction intervals cross zero. They are compatible with a larger or smaller mature human slope than the mature mouse slope; no cross-species equivalence or difference is claimed.

## Mouse lung identity and entropy

### New 18-month marker-set disorder

The primary all-cell global off-program results include 10 mature mice and 9 target cell types. Per-cell-type models require at least 6 donors; fibroblast of lung has only 5 and is explicitly UNTESTABLE in the per-cell-type arm. No per-cell-type off-program, stringent off-program, or detected-marker-count result has BH q below 0.05.

| Depth | Model | Effect per 15 months, SD | 95% interval | p | MDE80, SD |
| ---: | --- | ---: | ---: | ---: | ---: |
| 500 | Age plus sex | +0.081933 | -1.004424 to +1.168290 | 0.868303 | 1.510595 |
| 500 | Age, sex, actual library | -0.161994 | -1.115093 to +0.791106 | 0.709544 | 1.325298 |
| 1500 | Age plus sex | +0.015929 | -0.977929 to +1.009786 | 0.971870 | 1.381973 |
| 1500 | Age, sex, actual library | -0.146075 | -1.059490 to +0.767339 | 0.725873 | 1.270115 |
| 3000 | Age plus sex | +0.211903 | -0.627865 to +1.051671 | 0.582090 | 1.167709 |
| 3000 | Age, sex, actual library | +0.090756 | -0.671147 to +0.852659 | 0.793650 | 1.059437 |
| Full | Age plus sex | -0.149705 | -1.114297 to +0.814888 | 0.733601 | 1.341279 |
| Full | Age, sex, actual library | -0.255133 | -1.217017 to +0.706751 | 0.563293 | 1.337514 |

Balanced 20-cell sensitivities are also in the CSV. The intervals are wide and cross zero. The new marker-set analysis does not support a general increase in mature off-program expression, but it is underpowered and does not establish absence. The negative full-depth lung reference-distance slope is a separate endpoint and should not be conflated with off-program identity loss.

### Reused mature mouse entropy

All 144 existing global mature entropy rows are reused, covering 4 depths, 3 entropy definitions, 3 covariate specifications, and 4 eligibility families. No global entropy calculation was duplicated. In the manuscript-target 6-donor family, Shannon age-plus-sex effects per 15 months are +0.173321 at 500, +0.138029 at 1500, +0.138491 at 3000, and -0.171104 at full depth, with p = 0.712789, 0.778257, 0.776268, and 0.739454 respectively. Existing 80 percent-power MDEs are approximately 1.42 to 2.62 SD across the principal global age-plus-sex and age-only models. Isolated full-depth square-root KL signals remain formula, depth, and sex sensitive. Mature mouse entropy is inconclusive, not equivalent to zero.

### New 30-plus De Man entropy

All primary-seed intervals cross zero and no method has a BH q below 0.05. Across all 20 seeds, none of the 12 method by cell-type combinations has any p below 0.05.

| Method | Cell type | Donors | Slope per year | 95% interval | Primary p | Primary q |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Published KL | AT1 | 18 | -0.000809 | -0.003568 to +0.001950 | 0.539538 | 0.808493 |
| Published KL | AT2 | 17 | -0.000219 | -0.002134 to +0.001695 | 0.808493 | 0.808493 |
| Published KL | Aerocyte | 12 | +0.001034 | -0.003465 to +0.005533 | 0.610613 | 0.808493 |
| Published KL | gCap | 17 | -0.000873 | -0.003243 to +0.001496 | 0.440251 | 0.808493 |
| Shannon 1500 | AT1 | 19 | +0.001552 | -0.005201 to +0.008305 | 0.631348 | 0.808676 |
| Shannon 1500 | AT2 | 20 | +0.002668 | -0.002238 to +0.007574 | 0.265934 | 0.808676 |
| Shannon 1500 | Aerocyte | 13 | +0.004133 | -0.010337 to +0.018603 | 0.534307 | 0.808676 |
| Shannon 1500 | gCap | 17 | -0.000689 | -0.006709 to +0.005332 | 0.808676 | 0.808676 |
| Shannon full | AT1 | 20 | +0.003947 | -0.008837 to +0.016731 | 0.522074 | 0.696099 |
| Shannon full | AT2 | 20 | +0.002767 | -0.003697 to +0.009231 | 0.377654 | 0.696099 |
| Shannon full | Aerocyte | 13 | +0.006475 | -0.014391 to +0.027340 | 0.500459 | 0.696099 |
| Shannon full | gCap | 17 | -0.000005 | -0.007655 to +0.007645 | 0.998899 | 0.998899 |

The 30-plus De Man analysis therefore provides no statistically supported general entropy increase. It remains compatible with positive or negative age slopes within the reported intervals.

## Reused versus newly computed arms

| Arm | Disposition | Source or new output |
| --- | --- | --- |
| Mouse full-depth global mean Hellinger, unadjusted and linear depth | REUSED as-is | results/rev3-bayesian-vector-mature-only-frequentist.csv |
| Mouse full-depth pseudobulk, unadjusted and linear depth | REUSED as-is | results/rev3-bayesian-vector-mature-only-frequentist.csv |
| Ten Bayesian species-level mature specifications | REUSED as-is | results/rev3-bayesian-vector-mature-only-posterior.csv |
| All 144 mature global mouse lung entropy rows | REUSED as-is | results/rev3-entropy-mature-only-mouse-lung.csv |
| Full-depth mouse spline mean and pseudobulk models | NEW | results/followup-mature-drift-mouse.csv |
| Mouse 500, 1500, and 3000 cross-fitted 18-month references and all depth models | NEW | results/followup-mature-drift-mouse.csv |
| Mouse per-organ and leave-one-donor-out arms | NEW | results/followup-mature-drift-mouse.csv |
| All seven human reference and scored-donor rebuilds and depth models | NEW | results/followup-mature-drift-human.csv |
| Mature-only pooled interaction and 1000-replicate bootstrap | NEW | results/followup-mature-drift-interaction.csv |
| Mouse lung 18-month marker sets and off-program depth, balanced-cell, library, and power arms | NEW | results/followup-mature-h3-disorder.csv |
| De Man 30-plus rarefaction, entropy, and 20-seed models | NEW | results/followup-mature-h3-disorder.csv |

The reuse ledger contains 164 rows with file paths, values, intervals where available, donor counts, and specification identifiers: results/followup-mature-drift-reused.csv. The committed mixed-age 18-month donor table is used only as an input to filter mature full-depth endpoints and to compute new per-organ and influence models; its unrestricted global estimates are not reused as mature results.

## Runtime and run record

| Pipeline | Start UTC | End UTC | Duration seconds | Model | Status |
| --- | --- | --- | ---: | --- | --- |
| Mature H3 disorder | 2026-11-21T20:29:52.754153+00:00 | 2026-11-21T20:30:27.618228+00:00 | 34.864097 | kimi-code/k3.1-pro | completed |
| Reuse ledger | 2026-11-21T20:30:51.316482+00:00 | 2026-11-21T20:30:51.337442+00:00 | 0.020975 | kimi-code/k3.1-pro | completed |
| Mouse drift | 2026-11-21T20:30:51.338916+00:00 | 2026-11-21T20:33:41.178207+00:00 | 169.839303 | kimi-code/k3.1-pro | completed |
| Human drift | 2026-11-21T20:33:41.180019+00:00 | 2026-11-21T20:46:55.530683+00:00 | 794.350685 | kimi-code/k3.1-pro | completed |
| Interaction | 2026-11-21T20:46:55.532615+00:00 | 2026-11-21T20:46:56.955445+00:00 | 1.422848 | kimi-code/k3.1-pro | completed |

Combined active pipeline time was 1000.497908 seconds, or 16 minutes 40.498 seconds. Wall time from the H3 start through the interaction end was 17 minutes 4.201 seconds. Exact machine-readable runtimes are in results/followup-mature-drift-runtime.csv.

The post-run metadata-only rename wrote the required human donor CSV with reference_window_min and reference_window_max, then stopped while reading the expected zero-row ImmuneAging_BAL scratch checkpoint because that empty file has no CSV header. The required donor output had already been written and verified. No model was rerun and no numerical output changed. The scratch checkpoint is not a deliverable.

A small-kidney preflight completed reference and score calculations but its first direct function call failed when writing a checkpoint because scratch/mature-drift did not yet exist. No required result CSV existed at that point. Checkpoint directory creation was moved into the function, and the diagnosed preflight was rerun successfully before the full pipeline. This was a path initialization failure, not a model or memory failure.

The 38 targeted tests passed under the base analysis Python. An initial full-suite invocation with that interpreter stopped during collection because it lacks rdata; no tests ran in that invocation and it wrote no result files. The declared agescape project environment has rdata, and the full standard suite passed there: 195 passed, 1 skipped, with 1 unrelated implicit-modification warning.

## Verification limits and interpretation

The single self-review checked the scripts and every output table for age leakage, reference-window errors, duplicated entropy runs, untestable-row omissions, and historical-number wording. Programmatic checks passed: mouse donor ages are exactly 1.5, 1.75, 2.0, and 2.5 years; every human scored donor is at least 30; every human reference-window field is 30 to 39; every scored unit has at least 2 reference donors and at least 50 post-exclusion reference cells; mouse off-program donors are at least 18 months with an 18-month reference; and De Man model donors are at least 30. The historical +0.007229 appears only as an explicitly unrestricted mixed-age comparator. No Item 4 file contains an em dash or en dash, and no manuscript file changed.

The 36 mouse per-organ UNTESTABLE rows, 2 ImmuneAging_BAL UNTESTABLE rows, and 48 per-cell-type off-program UNTESTABLE rows remain explicit. The latter are the fibroblast 5-donor arm across depth, sampling, metric, and covariate combinations. None is interpreted as a null. Reference-window fields describe the design bounds, not the observed age extrema within each cohort.

Not independently verified: the local Census objects were not re-downloaded or compared with a remote second copy, and distances were not reimplemented in a second software package. Reused Bayesian chains and mature mouse entropy were not rerun by design. Existing Bayesian diagnostics all have convergence_ok = True, but their source table retains high Pareto-k warnings that were not independently resolved. Bootstrap interaction intervals use 1000 resamples and retain Monte Carlo resolution limits. Small donor counts, age-sex confounding, depth sensitivity, and wide MDEs preclude claims of equivalence or absence. Across the primary mouse, human, interaction, off-program, and De Man intervals, continuation in either direction remains compatible where reported.
