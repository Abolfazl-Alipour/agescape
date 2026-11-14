# Item 6: Bayesian Hierarchical Analysis of Whole-Vector Drift on Mature Mice Only (18 Months and Older)

## 1. Objective and Motivation

In the Item 3 audit (`results/rev3-young-vs-mature-audit-REPORT.md`), we identified an asymmetry in `results/reviewer-rev2-hellinger.csv` and `results/reviewer-rev2-hellinger-donors.csv` for `mouse_ref_18m`: although the Hellinger reference centroid was computed from the 18m mature mice, the donor-level regression in `mouse_ref_18m` still retained all 23 mice across 1m to 30m (mean age slope `+0.0072`, `p = 0.6656` after linear `log(UMI)` adjustment), while `results/h3-cross-species-bayesian-posterior.csv` modeled the full 1m to 30m lifespan against the young reference (`+0.1495`, `94% HDI [+0.0580, +0.2388]` unadjusted; `+0.0782`, `94% HDI [+0.0399, +0.1163]` depth-adjusted; `P(>0) = 0.9990`).

To isolate genuine post-maturity aging trajectories from early developmental changes (1m to 3m) without fitting an interior U-shaped centroid across young and old mice, this analysis restricts `results/reviewer-rev2-hellinger-donors.csv` (`table == "mouse_ref_18m"`) strictly to mature and old Tabula Muris Senis (TMS) mice and fits the exact same hierarchical Student-t regression framework (`scripts/bayesian_diagnostics_helper.py`, `scripts/run_rev3_bayesian_vector_mature_only.py`) alongside clustered frequentist OLS checks.

## 2. Dataset Composition and Design Structure

All mature mouse records in `mouse_ref_18m` (units with at least 6 mature donor observations) come from 14 TMS mice across 7 tissues (`Fat`, `Heart_and_Aorta`, `Kidney`, `Limb_Muscle`, `Lung`, `Marrow`, `Spleen`), 51 tissue-by-cell-type units, and 480 donor-cell-type observations:

| Subset | Ages included | Mice (Sex composition) | Tissues | Cell-type units | Observations |
| --- | --- | --- | --- | --- | --- |
| `mature_ge_18m` (primary) | 18m, 21m, 24m, 30m | 14 (4F, 10M): 18m (2F, 2M), 21m (2F), 24m (4M), 30m (4M) | 7 | 51 | 480 |
| `post_ref_gt_18m` (sensitivity) | 21m, 24m, 30m | 10 (2F, 8M): 21m (2F only), 24m (4M only), 30m (4M only) | 6 | 38 | 280 |
| `male_ge_24m` (unconfounded sensitivity) | 24m, 30m | 8 (8M only): 24m (4M), 30m (4M) | 4 | 20 | 128 |

Age is scaled identically to the paper's Bayesian vector model: units of 10% of the 4.0-year mouse maximum lifespan (`0.4` years = `4.8` months), with within-unit standardized outcomes (`(y - mean) / SD`), unit-centered `sex_male`, unit-centered `mean_log_umi` (in depth-adjusted models), donor random intercepts, tissue-level hierarchical age slopes (`group_age`), cell-type unit age slopes (`unit_age`), Student-t likelihood (`nu ~ Gamma(2, 0.1)`), and `global_age ~ Normal(0, prior_sd)` sensitivity across `prior_sd = 0.25` and `0.50` (`4` chains, `2,500` warmup + `2,500` sampling draws via `nutpie`, `target_accept = 0.99`).

## 3. Primary Posterior Results vs Frequentist Baselines

All numerical outputs are saved in `results/rev3-bayesian-vector-mature-only-posterior.csv`, `results/rev3-bayesian-vector-mature-only-frequentist.csv`, `results/rev3-bayesian-vector-mature-only-diagnostics.csv`, and `results/rev3-bayesian-vector-mature-only-prior-predictive.csv`.

### 3.1 Species-level Mature Drift Estimates (`SD` per 10% mouse lifespan = `4.8` months)

| Specification | Subset & Endpoint | Model & Prior SD | Bayes Mean (SD) | 94% HDI | `P(>0)` | Freq OLS (SE) | Freq `p` |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Primary | `mature_ge_18m` `profile_hellinger` | depth-adjusted, `0.25` | `+0.0338` (`0.0680`) | `[-0.0992, +0.1584]` | `0.6950` | `-0.0137` (`0.0579`) | `0.8163` |
| Primary | `mature_ge_18m` `profile_hellinger` | depth-adjusted, `0.50` | `+0.0389` (`0.0728`) | `[-0.0998, +0.1773]` | `0.7068` | `-0.0137` (`0.0579`) | `0.8163` |
| Sensitivity | `mature_ge_18m` `profile_hellinger` | unadjusted, `0.25` | `+0.1440` (`0.1299`) | `[-0.0898, +0.3997]` | `0.8699` | `+0.1100` (`0.1061`) | `0.3187` |
| Sensitivity | `mature_ge_18m` `profile_hellinger` | unadjusted, `0.50` | `+0.1770` (`0.1467`) | `[-0.0944, +0.4503]` | `0.8877` | `+0.1100` (`0.1061`) | `0.3187` |
| Sensitivity | `mature_ge_18m` `pseudobulk_hellinger` | depth-adjusted, `0.50` | `+0.0765` (`0.1597`) | `[-0.2251, +0.3813]` | `0.6959` | `+0.0350` (`0.0764`) | `0.6546` |
| Sensitivity | `mature_ge_18m` `pseudobulk_hellinger` | unadjusted, `0.50` | `+0.0959` (`0.1621`) | `[-0.2244, +0.3935]` | `0.7306` | `+0.0658` (`0.0827`) | `0.4409` |
| Sensitivity (sex-split caveat) | `post_ref_gt_18m` `profile_hellinger` | depth-adjusted, `0.50` | `+0.2300` (`0.0960`) | `[+0.0509, +0.4072]` | `0.9819` | `+0.1176` (`0.0629`) | `0.0945` |
| Sensitivity (sex-split caveat) | `post_ref_gt_18m` `profile_hellinger` | unadjusted, `0.50` | `+0.3070` (`0.1721`) | `[-0.0155, +0.6340]` | `0.9601` | `+0.1792` (`0.1222`) | `0.1765` |
| Sensitivity (unconfounded) | `male_ge_24m` `profile_hellinger` | depth-adjusted, `0.50` | `+0.2580` (`0.2931`) | `[-0.3078, +0.8059]` | `0.8259` | `+0.1634` (`0.1199`) | `0.2152` |
| Sensitivity (unconfounded) | `male_ge_24m` `profile_hellinger` | unadjusted, `0.50` | `+0.2758` (`0.3411`) | `[-0.3446, +0.9329]` | `0.8046` | `+0.3134` (`0.1660`) | `0.1010` |

For comparison, the full-range 1m to 30m mouse Bayesian estimates from `results/h3-cross-species-bayesian-posterior.csv` (`n = 23` donors) are `+0.1495` (`94% HDI [+0.0580, +0.2388]`, `P(>0) = 0.9986`) unadjusted and `+0.0782` (`94% HDI [+0.0399, +0.1163]`, `P(>0) = 0.9990`) depth-adjusted.

### 3.2 Key Biological and Statistical Conclusions

1. **No detectable global post-maturity whole-vector drift once sequencing depth is controlled**: Across all 14 mature mice (`18m to 30m`, 51 cell-type units, 480 observations), the primary depth-adjusted posterior mean is `+0.0389` (`94% HDI [-0.0998, +0.1773]`, `P(>0) = 0.7068`, `prior_sd = 0.50`) and `+0.0338` (`94% HDI [-0.0992, +0.1584]`, `P(>0) = 0.6950`, `prior_sd = 0.25`), matching the frequentist null (`-0.0137`, `SE = 0.0579`, `p = 0.8163`). The pseudobulk Hellinger endpoint agrees (`+0.0765`, `94% HDI [-0.2251, +0.3813]`, `P(>0) = 0.6959`; frequentist `p = 0.6546`). Thus, the Bayes framework confirms the core conclusion from Item 3: mouse whole-profile Hellinger drift is largely a young-to-mature transition rather than a universal mature-to-old trajectory.
2. **Library depth strongly dominates cell-level Hellinger distances**: Including `mean_log_umi` increases expected log pointwise predictive density (`elpd_loo`) on `mature_ge_18m` by `229.8` units (`-463.68` depth-adjusted vs `-693.48` unadjusted, `global_depth = -3.5191`, `94% HDI [-3.7907, -3.2568]`) and shrinks the apparent unadjusted slope from `+0.1770` (`P(>0) = 0.8877`) down to `+0.0389`.
3. **Why `post_ref_gt_18m` (`> 18m` only) diverges and why it is sex-confounded**: Excluding the four `18m` mice (2F, 2M) leaves only 2 females at `21m` and 8 males at `24m` and `30m`. Sex (`global_sex`) and age become strongly collinear across those 10 animals. In unpenalized OLS (`results/rev3-bayesian-vector-mature-only-frequentist.csv`), adjusting for sex and library depth yields a non-significant slope of `+0.1176` (`SE = 0.0629`, `95% CI [-0.0248, +0.2600]`, `p = 0.0945`). In the hierarchical Student-t model, the `Normal(0, 0.5)` prior regularizes `global_sex` (`+0.2653`, `94% HDI [-0.0107, +0.5386]`) and transfers variance onto the collinear `global_age` (`+0.2300`, `94% HDI [+0.0509, +0.4072]`). When sex confounding is removed completely in `male_ge_24m` (`24m vs 30m`, 8 male mice, 4 tissues, 20 units), the posterior interval spans zero broadly (`+0.2580`, `94% HDI [-0.3078, +0.8059]`, `P(>0) = 0.8259`; OLS `+0.1634`, `p = 0.2152`).

### 3.3 Tissue-level Heterogeneity (`mature_ge_18m`, depth-adjusted, `prior_sd = 0.50`)

Tissue-level posteriors are consistent with zero overall drift and show modest tissue-to-tissue spread centered near zero:
- `Limb_Muscle`: `+0.0974` (`94% HDI [-0.0492, +0.2602]`, `P(>0) = 0.8864`); top unit `Limb_Muscle | TMS | macrophage` (`+0.1073`, `P(>0) = 0.8894`). Note in unconfounded `male_ge_24m` (`24m vs 30m` males), `Limb_Muscle` is positive (`+0.8061`, `94% HDI [+0.3462, +1.2719]`, `P(>0) = 0.9982`) whereas `Kidney` is negative (`-0.0870`, `94% HDI [-0.6014, +0.4426]`).
- `Kidney`: `+0.0519` (`[-0.0767, +0.1828]`, `P(>0) = 0.7816`)
- `Spleen`: `+0.0428` (`[-0.0984, +0.1777]`, `P(>0) = 0.7284`)
- `Fat`: `+0.0404` (`[-0.1397, +0.2200]`, `P(>0) = 0.6697`)
- `Heart_and_Aorta`: `+0.0360` (`[-0.1301, +0.2090]`, `P(>0) = 0.6605`)
- `Marrow`: `+0.0357` (`[-0.0841, +0.1530]`, `P(>0) = 0.7112`)
- `Lung`: `-0.0284` (`[-0.1650, +0.1030]`, `P(>0) = 0.3328`)

## 4. Convergence, Prior Predictive, and Posterior Predictive Diagnostics

All 10 hierarchical models passed every MCMC convergence criterion (`results/rev3-bayesian-vector-mature-only-diagnostics.csv`):
- Maximum R-hat across all parameters and models: `<= 1.0070` (threshold `<= 1.01`; primary `mature_ge_18m` depth-adjusted `max_rhat = 1.0026`).
- Minimum bulk ESS: `579.9` (`1859.7` on primary `prior_0.50`, well above `400`).
- Minimum tail ESS: `800.2` (`2764.9` on primary `prior_0.50`, well above `400`).
- Divergent transitions: `0` across all 10 models (`min_bfmi` in `[0.8000, 0.8844]`).
- Prior predictive check: 100% of observed standardized values fall inside the central 99.8% prior predictive interval across all 10 models (`observed_within_prior_998 == True`).
- Posterior predictive 90% empirical coverage: `[0.9464, 0.9844]` across all 10 models (`0.9521` on primary `prior_0.50`). As expected with small `n = 6 to 9` observations per cell-type unit and Student-t likelihood with donor random intercepts, `8` to `19` pointwise LOO `pareto_k` values exceed `0.7` (`pareto_k_max` in `[0.9460, 1.1402]`).

## 5. Artifacts and Reproducibility

- Analysis script: `scripts/run_rev3_bayesian_vector_mature_only.py` (reusing `scripts/bayesian_diagnostics_helper.py`), executed under `scripts/resource_guard.py`.
- Result tables:
  - `results/rev3-bayesian-vector-mature-only-posterior.csv` (`507` rows)
  - `results/rev3-bayesian-vector-mature-only-diagnostics.csv` (`10` rows)
  - `results/rev3-bayesian-vector-mature-only-prior-predictive.csv` (`10` rows)
  - `results/rev3-bayesian-vector-mature-only-frequentist.csv` (`10` rows)
