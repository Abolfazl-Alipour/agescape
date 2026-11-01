# Bayesian analysis of whole-expression-vector displacement

## Bottom line

Bayesian analysis supports a positive mouse association and does not establish an overall human or NMR association. Human unadjusted estimates vary by model, while depth-adjusted human estimates are close to zero across prior choices. NMR has a negative posterior mean with a wide interval. All human global models and NMR remain exploratory because of residual convergence, influence, or calibration warnings. Frequentist results remain the reference analysis.

## Species-wide posterior estimates

Effects are standardized within cell-type units per 10 percent of species maximum lifespan. Intervals are 94 percent highest-density intervals.

| Species and model | Posterior mean | 94% HDI | Probability of positive effect |
|---|---:|---:|---:|
| Mouse, unadjusted | +0.1495 | +0.0580 to +0.2388 | 0.9986 |
| Mouse, depth-adjusted | +0.0782 | +0.0399 to +0.1163 | 0.9990 |
| Human global, unadjusted | +0.0893 | +0.0345 to +0.1447 | 0.9976 |
| Human global, depth-adjusted | +0.0133 | -0.0081 to +0.0344 | 0.8816 |
| NMR, unadjusted | -0.1581 | -0.7210 to +0.4798 | 0.2949 |

The positive human unadjusted global posterior is not robust to model choice. An alternative human hierarchical model with fixed cohort offsets and random cell-type age slopes gives a mean of +0.0278 with a 94% HDI from -0.1704 to +0.2152 and posterior probability of an increase of 0.6126. The frequentist unadjusted human effect is also non-significant. The human global result must therefore not be presented as established positive aging signal.

## Human depth and prior sensitivity

The depth-adjusted human global posterior is stable across global-age priors and every interval includes zero.

| Prior SD | Posterior mean | 94% HDI | Probability of positive effect |
|---:|---:|---:|---:|
| 0.25 | +0.0133 | -0.0081 to +0.0342 | 0.8828 |
| 0.50 | +0.0133 | -0.0081 to +0.0344 | 0.8816 |
| 1.00 | +0.0142 | -0.0070 to +0.0359 | 0.8996 |

This agrees with the frequentist and gene-level audits: depth accounts for much of the apparent overall human association. It does not negate the independently replicated interneuron gene-expression changes, which are a separate estimand.

## Diagnostics

All models have zero divergences and adequate effective sample sizes. Mouse hierarchical models pass the workflow's convergence assessment. The human global models have maximum R-hat values from 1.0115 to 1.0152, slightly above the 1.01 threshold. All models have influential-observation Pareto-k warnings. NMR also has poor predictive calibration, and its five animals make the hierarchical posterior especially exploratory.

The analysis used Student-t likelihoods, donor intercepts, and unit residual scales. Mouse models include tissue and cell-type age deviations. The human global model centers age within each cell-type unit and estimates a common within-unit slope. Human depth-adjusted random-slope models did not converge acceptably and are not used as evidence.

## Interpretation

Mouse shows the clearest whole-vector age association in both statistical frameworks. Human does not show a robust overall association after depth adjustment. NMR has a negative point estimate but no precise evidence of stability or a between-species difference. The strongest human result is not the global vector posterior; it is the depth-adjusted and independently replicated interneuron gene-expression analysis reported in `report/human-expression-vector-interpretation.md`.

Methods are in `notes/cross-species-bayesian-vector-methods.md`. Posterior estimates are in `results/h3-cross-species-bayesian-posterior.csv`, and diagnostics are in `results/h3-cross-species-bayesian-diagnostics.csv`.
