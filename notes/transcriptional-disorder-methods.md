# Transcriptional-disorder analysis

## Hypothesis and estimand

The biological hypothesis is that age-related epigenetic dysregulation reduces transcriptional fidelity, allowing broad or cell-type-inappropriate transcription. Single-cell RNA counts can test two downstream predictions: increased within-cell Shannon entropy and increased expression of genes characteristic of other cell types. They cannot identify DNA damage, chromatin-mark restoration, or any other upstream cause.

This differs from cell-to-cell variability. Cell-to-cell variability compares cells with other cells. Shannon entropy measures how broadly one cell distributes its UMIs across genes. Off-program expression measures the fraction of a cell's UMIs assigned to genes that are uncommon in young cells of its own type but characteristic of another cell type.

## Primary analysis

The analysis used mouse lung singlet cells with at least 1,500 UMIs. Counts were downsampled without replacement to exactly 1,500 UMIs per cell. Mitochondrial, ribosomal, and sex-linked marker genes were excluded from entropy calculations. Shannon entropy was calculated in natural-log units after renormalizing each cell over the retained genes. A fixed gene universe was used across ages.

For each target cell type, off-program genes were detected in at most 5 percent of young target cells and at least 20 percent of cells from another lung cell type. The other cell type also had to express the gene at least 20 counts per 10,000 and at least fivefold above the target level after adding one count per 10,000 to the target denominator. References used animals aged 1 or 3 months. Each young animal was left out of its own reference, including references for other cell types. Old animals used the full young reference.

Animals with at least 20 eligible cells of the target type were analyzed. Animal means were regressed on age and sex. Animals, not cells, were the inferential units. Benjamini-Hochberg correction was applied separately across cell types for entropy and off-program UMI fraction. The nine prespecified cell types matched the earlier cell-to-cell variability analysis.

## Sensitivity analyses

The analysis was repeated with exactly 20 randomly selected cells per animal, after excluding pooled animal identifiers, after additional exclusion of immunoglobulin, hemoglobin, and immediate-early genes, and using the number of detected off-program genes rather than their UMI fraction. These were sensitivity analyses, not replacements selected because of their results.

## Interpretation limits

Entropy can change because of shifts in metabolic state, cell size, RNA composition, or genuine regulatory changes. Off-program expression can arise from ambient RNA, imperfect doublet removal, annotation errors, or biological state transitions. Equal-depth sampling reduces UMI-depth artifacts but does not remove these alternatives. A null result does not demonstrate stable epigenetic regulation.

## Bayesian hierarchical analysis

Outcomes were standardized within each cell type. Age was expressed as a 15-month contrast from 3 months, and sex was centered within cell type. The model contained a global age effect, non-centered cell-type intercepts and age slopes, a non-centered animal intercept, and cell-type residual scales. A Student-t likelihood reduced the influence of outlying animal-cell-type measurements. The animal intercept accounted for measurements from the same animal across cell types.

The primary global age prior was Normal(0, 0.5) in standard-deviation units per 15 months. Prior sensitivity analyses used Normal(0, 0.25) and Normal(0, 1). Other priors were Normal(0, 0.5) for the global intercept and sex effect, Gamma(2, 2) for cell-type and animal scales, and Gamma(2, 0.1) for the Student-t degrees of freedom. These priors were chosen for standardized outcomes and were checked with prior predictive sampling before fitting.

Models were sampled with PyMC and nutpie using descriptive reproducible seeds. Posterior summaries report means, 94 percent highest-density intervals, and posterior probabilities of a positive age effect. Diagnostics used the pinned `bayesian-workflow` revision `18ab0de9b63b5bb03caa96146235e8f1534d45d4`, including convergence, leave-one-out influence, posterior predictive, and probability-integral-transform checks. Full NetCDF and diagnostic artifacts remain local under `results/h3-bayesian/`; aggregate outputs are committed separately.

For a matching frequentist global analysis, cell-type fixed intercepts, age, and centered sex were fitted by ordinary least squares with animal-clustered standard errors. This was compared with the Bayesian global effect rather than confusing the global effect with the separate cell-type tests.

## Outputs

- `results/h3-transcriptional-disorder-donors.csv`: animal-level measurements and sampling labels.
- `results/h3-transcriptional-disorder.csv`: frequentist cell-type effect estimates, confidence intervals, p values, and q values.
- `results/h3-offprogram-marker-provenance.csv`: off-program marker counts and source cell types for each cross-fitted reference.
- `results/h3-bayesian-posterior.csv`: global and cell-type posterior effects for all prior and sampling sensitivity analyses.
- `results/h3-bayesian-prior-predictive.csv`: prior predictive range checks.
- `results/h3-bayesian-diagnostics.csv`: convergence and leave-one-out diagnostics for every fitted model.
- `results/h3-frequentist-global.csv`: animal-clustered frequentist global effects matching the Bayesian estimand.
- `figures/h3-transcriptional-disorder.png`: primary frequentist cell-type estimates and 95 percent confidence intervals.
