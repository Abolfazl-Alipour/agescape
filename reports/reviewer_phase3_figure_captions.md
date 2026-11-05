# Phase 3 figure captions

## Mouse vector drift by cell type

Files: `figures/reviewer-phase3-mouse-vector-drift.png` and `.pdf`. Data: `results/reviewer-phase3-mouse-drift-figure-data.csv`.

Each point is a donor-level age slope in mean Hellinger distance from a cross-fitted young expression reference for one mouse tissue-by-cell-type unit. Bars are 95% confidence intervals. Positive values indicate increasing distance from the young reference. The figure does not show increased random noise, DNA damage, or a universal aging program, and it does not adjust the displayed unit-level intervals for multiplicity.

## Interneuron gene changes

Files: `figures/reviewer-phase3-interneuron-gene-changes.png` and `.pdf`. Data: `results/reviewer-phase3-interneuron-figure-data.csv`.

Each point is a gene-by-interneuron-subtype MSSM association from depth-adjusted donor-pseudobulk models. Color is the age slope in log counts per 10,000 per year, and point size is minus log10 of the within-cell-type BH q-value. Genes include the top ten FDR-significant genes per subtype and genes shared across VIP, pvalb, and sst sets. The figure displays MSSM effect estimates only. It does not show HBCC effect sizes, prove causation, or establish that every displayed gene replicates at FDR in HBCC.

## Composition changes with confidence intervals

Files: `figures/reviewer-phase3-composition-changes.png` and `.pdf`. Data: `results/reviewer-phase3-composition-figure-data.csv`.

Points are FDR-significant mouse or human cell-fraction age slopes after adjustment for log total recovered cells, sex when identifiable, and the fixed 0.5 pseudocount. Bars are 95% confidence intervals. The outcome is centered log-ratio abundance among recovered annotated cells. The figure does not measure absolute tissue cell numbers, prove fibrosis, or eliminate age-related recovery and annotation biases.

## Gene-effect concordance

Files: `figures/reviewer-phase3-universality-concordance.png` and `.pdf`. Data: `results/reviewer-phase3-universality-figure-data.csv`.

Blue circles show observed mean pairwise Pearson correlations between cell-type gene-level age-effect vectors within tissues or cohorts. Orange squares and bars show the mean and 95% interval of the age-permutation null. Ages were permuted within cell type while preserving sex or platform strata. NMR used all ten age assignments. The figure does not show that all genes change, that whole-vector distance increases universally, or that the same genes respond across species.
