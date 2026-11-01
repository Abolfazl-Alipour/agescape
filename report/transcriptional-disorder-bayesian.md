# Frequentist and Bayesian analysis of transcriptional disorder

## Executive summary

The two frameworks agree on the direction and uncertainty of the main effects. Neither supports an overall increase in Shannon entropy with age. Both estimate a positive average off-program trend across cell types, but its uncertainty includes zero. The Bayesian model gives a posterior probability of 0.841 for a positive off-program effect, not proof that the biological hypothesis is true. Influential-animal and calibration warnings make this an exploratory signal.

For this dataset, the frequentist analyses should remain the confirmatory reference. The Bayesian hierarchical model is useful for exploratory effect estimation and partial pooling, but its predictive and influence diagnostics are not strong enough to make it the decisive analysis.

## Data and question

The analysis used animal-level measurements from Tabula Muris Senis lung cells. The entropy model contained 105 cell-type-by-animal observations from 16 animals and 9 cell types. The off-program model contained 96 observations from 16 animals and 8 cell types. The same animal could contribute to multiple cell types.

The question is whether transcriptional disorder increases over a 15-month age contrast. Outcomes were standardized within cell type, so effects are in cell-type standard-deviation units over 15 months. All models accounted for sex. Cells were not treated as independent biological replicates.

## Model specification

The Bayesian model included a global age effect, cell-type intercepts and age slopes, animal random intercepts, and cell-type residual scales. Non-centered parameterization was used for the hierarchical effects. A Student-t likelihood reduced the influence of outlying observations.

| Parameter | Prior | Reason |
|---|---|---|
| Global age effect | Normal(0, 0.5) | Weakly informative on the standardized 15-month scale |
| Global intercept and sex effect | Normal(0, 0.5) | Centered outcomes and centered sex make large values less plausible |
| Cell-type and animal scales | Gamma(2, 2) | Positive scales with shrinkage away from extreme values |
| Student-t degrees of freedom | Gamma(2, 0.1) | Allows heavy tails while retaining finite-variance mass |

The matching frequentist global model used cell-type fixed intercepts, age, and centered sex with animal-clustered standard errors. The separate cell-type frequentist tests retained Benjamini-Hochberg correction.

## Prior predictive check

Prior predictive samples covered the observed standardized outcomes. The entropy observations ranged from -2.158 to +1.921, while the central 99.8 percent prior predictive interval was -14.544 to +13.794. Off-program observations ranged from -1.626 to +2.802, with a prior interval of -13.073 to +12.552.

The priors therefore did not exclude the observed data, but they produced much more extreme values than observed. About 7 percent of prior predictions had magnitude above 5 standard deviations. This is consistent with heavy-tailed robust modeling, but it also warns against treating the posterior as prior-free evidence.

## Sampling and convergence

All 18 models in the prior and sampling sensitivity grid passed convergence checks with zero divergences. Longer tuning and sampling were required after warnings in an initial grid.

Convergence does not establish model adequacy. Every model had leave-one-out Pareto-k warnings for influential observations. In the primary off-program model, 8 observations had Pareto k above 0.7, with a maximum of 1.055. Influential observations occurred in both young and old animals and in pooled and unpooled identifiers. The primary off-program calibration check was under-confident, with mean coverage deviation 0.0979. These are substantive limitations on interpretation.

## Primary global effects

| Endpoint | Frequentist effect | Frequentist 95% CI | p | Bayesian posterior mean | Bayesian 94% HDI | P(effect > 0) |
|---|---:|---:|---:|---:|---:|---:|
| Entropy | -0.079 | -0.479 to +0.321 | 0.699 | -0.060 | -0.508 to +0.356 | 0.395 |
| Off-program UMI fraction | +0.244 | -0.196 to +0.684 | 0.278 | +0.243 | -0.224 to +0.702 | 0.841 |

The entropy estimate is centered slightly below zero in both frameworks. This does not support the predicted increase, but it is not proof of no change.

The off-program estimate is positive in both frameworks. The frequentist interval and Bayesian HDI include zero. The posterior probability of direction is approximately 84 percent, which is weak evidence for a positive trend under the fitted model and priors. It is not an 84 percent probability that epigenetic dysregulation causes aging.

## Prior and sampling sensitivity

Across prior standard deviations of 0.25, 0.5, and 1.0, the off-program posterior probability of an increase ranged from 0.793 to 0.857 in the primary sampling design. It ranged from 0.731 to 0.815 after excluding pooled animal identifiers and from 0.767 to 0.835 with 20 cells per animal. Every off-program HDI included zero.

The corresponding entropy probabilities of an increase ranged from 0.280 to 0.442 across all sensitivity models. Every entropy HDI included zero.

The separate frequentist cell-type tests gave no FDR-significant off-program result. Bayesian cell-type estimates were partially pooled toward the global mean, and every cell-type HDI also included zero. The positive Bayesian global trend is therefore a cross-cell-type average signal, not a confirmed effect in any particular cell type.

## Which framework makes more sense here?

For confirmatory claims, the frequentist analysis is the better reference for three reasons. The animal-clustered global model is comparatively simple, the cell-type analyses have explicit false-discovery control, and the Bayesian models retain influential-observation and calibration warnings despite clean convergence.

For exploratory estimation, the Bayesian hierarchical model adds useful structure. It estimates an overall effect and cell-type deviations simultaneously, accounts for repeated animals, and shrinks unstable cell-type estimates instead of interpreting each small sample independently. Its posterior probabilities describe effect direction under explicit assumptions, which is useful as long as those assumptions and diagnostic limitations are reported.

The appropriate choice is therefore frequentist confirmation with Bayesian exploratory estimation, not selection of whichever framework gives the stronger result.

## Conclusion and next steps

The third hypothesis is not established. The entropy prediction is unsupported in both frameworks. The off-program prediction has a weak positive average trend, but it is sensitive to animal influence, model assumptions, and sampling choices. The data do not establish DNA damage, failed chromatin restoration, or inflammation as a cause.

A stronger test requires more animals, preregistered effect definitions, technical covariates or empty-droplet ambient-RNA estimates, and chromatin measurements linked to the same cells. Numerical outputs are in `results/h3-bayesian-posterior.csv`, `results/h3-frequentist-global.csv`, and `results/h3-bayesian-diagnostics.csv`. Full methods are in `notes/transcriptional-disorder-methods.md`.
