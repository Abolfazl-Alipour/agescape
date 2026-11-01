# Bayesian whole-expression-vector models

## Estimands

The Bayesian analysis uses the committed animal- or donor-level Hellinger distances. Outcomes are standardized within each tissue-by-cohort-by-cell-type unit. Age is scaled per 10 percent of species maximum lifespan: 0.4 years for mouse, 12.25 years for human, and 3.1 years for NMR. Positive effects mean increasing distance from the cross-fitted young expression reference.

The primary Bayesian models use a Student-t observation model. Mouse and human are fitted with and without mean-log-UMI adjustment. Human depth-adjusted models are repeated with global-age prior standard deviations of 0.25, 0.5, and 1.0. NMR is fitted separately because it has five animals and two skin compartments.

## Hierarchy

Mouse and human models contain a species-wide age effect, tissue- or cohort-level age deviations, cell-type-unit age deviations, unit intercepts, donor intercepts, a global sex effect, and unit-specific residual scales. Mouse tissue slopes are partially pooled toward the species effect. Human cohort slopes use fixed offsets with Normal(0, 0.25) priors after the estimated cohort-scale model failed convergence checks. Cell-type-unit slopes remain partially pooled toward their tissue or cohort slope. Donor intercepts account for repeated measurements across cell types.

The NMR model omits the tissue hierarchy and sex adjustment. It retains unit intercepts, unit age deviations, animal intercepts, and unit residual scales. Its posterior should be interpreted together with the exact ten-assignment frequentist test and the small sample size.

Human depth-adjusted random-slope models did not converge acceptably and are not used as evidence. The human global analysis therefore centers age within each cell-type unit and estimates a common within-unit age slope with donor intercepts, global sex and depth effects, and unit residual scales. Unit intercepts are unnecessary after within-unit centering of the standardized outcome and age. The same global model is fitted without depth adjustment for comparison. Prior-sensitivity analyses use global-age standard deviations of 0.25, 0.5, and 1.0.

## Priors and computation

The primary global-age prior is Normal(0, 0.5) in standardized-outcome units per 10 percent of maximum lifespan. Other priors are Normal(0, 0.5) for global intercept, sex, and depth effects; Gamma(2, 2) for intercept, age-deviation, donor, and residual scales; and Gamma(2, 0.1) for Student-t degrees of freedom. Unit and donor effects use non-centered parameterizations.

Models are sampled with PyMC and nutpie using target acceptance 0.99 and descriptive reproducible seeds. Mouse and NMR use 2,000 tuning steps and 2,000 posterior draws per chain. Human age is centered at 40 years to reduce intercept-age correlation, and human models use 4,000 tuning steps and 4,000 draws after shorter initial fits had convergence flags. Prior-predictive sampling checks whether the priors produce plausible standardized outcomes. Full inference objects and diagnostic artifacts remain local under `results/h3-cross-species-bayesian/`.

## Diagnostics and interpretation

Diagnostics use the pinned `bayesian-workflow` revision already used by the project. Convergence, leave-one-out Pareto-k influence, posterior prediction, and probability-integral-transform calibration are recorded for every model. Any global effect with a convergence flag must be refit before interpretation. Influential-observation warnings limit confirmatory interpretation even when sampling converges.

Bayesian probabilities are posterior probabilities under the stated model, not frequentist p values. Partial pooling can reduce noisy cell-type estimates, but it cannot repair confounding by sequencing depth, assay, annotation, or cohort design. Unadjusted and depth-adjusted results must therefore be reported together.

## Outputs

- `results/h3-cross-species-bayesian-posterior.csv`: species, tissue or cohort, and cell-type-unit posterior effects.
- `results/h3-cross-species-bayesian-prior-predictive.csv`: prior-predictive range checks.
- `results/h3-cross-species-bayesian-diagnostics.csv`: aggregate convergence and model-diagnostic results.
