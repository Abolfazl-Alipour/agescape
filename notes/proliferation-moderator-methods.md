# Proliferation-class moderator analysis

## Question and status

This is a post hoc exploratory analysis suggested after the cross-species results were known. It asks whether annotated post-mitotic cells have different aging effects from proliferative-capable cells. It is not a preregistered hypothesis and does not measure proliferation directly.

## Biological classification

Post-mitotic means that the annotated mature cell has exited the cell cycle and normally does not divide. Quiescence is not equivalent: a quiescent cell can re-enter the cell cycle. Proliferative-capable includes stem cells, progenitors, and mature cells that can divide during homeostasis, activation, or injury.

The classification applies to the annotated cell, not to the turnover of its tissue. For example, an intestinal stem cell is proliferative-capable, while a mature enterocyte is post-mitotic even though the epithelium is continually renewed. Mixed or generic labels, including broad germ-cell and neural-cell labels, are excluded. Meiotic spermatocytes and mixed germline compartments form a separate excluded category because meiosis is not mitotic proliferation.

The deterministic label rules are in `src/proliferation.py`. Every resulting label and its rationale are exported in `results/proliferation-cell-type-classification.csv`. Classification uses only species and cell-type identity, never aging effect sizes or significance.

## Included analyses

Expression-vector analysis uses depth-adjusted Hellinger displacement in mouse, human, AFCA fly, Gao worm with the balanced day-6 reference, and fine zebrafish testis cell types. The unadjusted analysis is a sensitivity arm. Cell-fraction analysis uses CLR values in the same systems. The human observed-record-only composition analysis is a sensitivity arm. The fine zebrafish composition input is reconstructed from the same annotated atlas and filtering rules used by the original analysis.

NMR is not included in the moderator inference. Five animals and two evaluable expression compartments cannot support a general class comparison. The Roux worm arm is excluded because its libraries are technical replicates.

## Unit effects

The inferential unit is one cell type within one tissue and cohort. At least six donors, animals, or pools are required. Each unit receives an age-slope model with sex when it varies, sequencing platform for Gao worm, and mean log UMI in the depth-adjusted vector arm.

Effects are multiplied by the existing species age scale and divided by the unit's donor-level outcome standard deviation. The resulting unit is outcome SD per 10% of maximum lifespan. The age scales are 0.4 years for mouse, 12.25 years for human, 40 days for fly, 14 days for worm, and 22 months for zebrafish. This scaling aligns effect magnitudes but does not make assays or age coverage equivalent.

## Frequentist comparison

The primary contrast is post-mitotic minus proliferative-capable. Ordinary least squares gives each cell-type unit equal weight and includes tissue-by-cohort fixed intercepts. Only strata containing both classes contribute. HC3 heteroskedasticity-consistent standard errors and 95% confidence intervals are reported.

Signed effects test whether one class has a more positive aging slope. Absolute effects test whether one class has larger changes in either direction. Absolute effects are necessary for cell fractions because increases and decreases are both biologically relevant. They are exploratory because taking an absolute value introduces estimation bias for noisy null effects.

A label permutation test shuffles the biological class only within tissue-by-cohort strata. It uses 2,000 permutations and a fixed seed. This tests whether the observed contrast is unusual given the within-stratum class counts. Benjamini-Hochberg correction is applied across the five species separately for each analysis, sensitivity arm, and response. Combined-species models are secondary summaries and are not included in that correction.

## Bayesian comparison

The Bayesian analysis uses the signed unit effects and their standard errors. A normal measurement model represents the estimated unit slope. Tissue-by-cohort intercepts absorb average stratum effects. A unit-level random effect with scale tau allows biological heterogeneity beyond the reported slope uncertainty. The post-mitotic minus proliferative-capable contrast has a Normal prior centered at zero. Prior SDs 0.25, 0.5, and 1.0 assess prior sensitivity. Stratum intercepts have Normal(0, 1) priors and tau has a HalfNormal(0, 0.5) prior.

Posterior means, 94% highest-density intervals, posterior probabilities of a positive contrast, R-hat, effective sample size, and divergences are reported. The Bayesian model is precision-aware, while the frequentist model gives cell types equal weight. Agreement between them is more informative than either alone. Absolute effects do not receive the normal measurement-error model because their error distribution is not normal after taking an absolute value.

## Reproduction

```bash
python scripts/run_proliferation_moderator.py --permutations 2000
python scripts/run_proliferation_bayesian.py --draws 2000
python scripts/plot_proliferation_moderator.py
```

The Bayesian command uses the environment in `environment-bayesian.yml`. The frequentist command uses committed public donor-level vector results and local atlas inputs already required by the composition pipelines.

## Outputs

- `results/proliferation-cell-type-classification.csv`
- `results/proliferation-unit-effects.csv`
- `results/proliferation-group-summary.csv`
- `results/proliferation-moderator-tests.csv`
- `results/proliferation-bayesian-posterior.csv`
- `figures/proliferation-moderator.png`
