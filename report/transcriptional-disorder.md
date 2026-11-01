# Does aging increase transcriptional disorder within cells?

## Bottom line

The mouse lung data do not establish a general age-related increase in transcriptional disorder. No cell type had an FDR-significant increase in within-cell Shannon entropy or off-program expression in the primary frequentist analysis. Entropy estimates changed with cell-number balancing, including a positive fibroblast result and negative bronchial smooth-muscle and vein endothelial results. A donor-adjusted Bayesian hierarchical analysis also found no entropy increase. It estimated a weak positive off-program trend across cell types, but all posterior intervals included zero and the trend weakened under some sensitivity analyses. The third hypothesis is not established, and transcriptomics does not test its proposed DNA-repair mechanism.

## What the hypothesis means

The proposed mechanism is that DNA damage and repair disturb chromatin regulation, causing cells to lose some transcriptional specificity. The transcriptomic prediction is broader expression across genes, increased expression of inappropriate cell-type programs, or both. Inflammation is a possible downstream expression phenotype, but it is not specific to this mechanism.

Shannon entropy measures how evenly a cell's transcripts are distributed across genes. Off-program expression measures transcripts from genes characteristic of other cell types but rarely detected in young cells of the target type. These are more direct readouts of the hypothesis than cell-to-cell variability, although they remain indirect and vulnerable to technical confounding.

## Primary results

All models used animal-level means, age and sex covariates, and cells downsampled without replacement to 1,500 UMIs. Slopes are per month. Off-program slopes are changes in UMI fraction per month.

| Cell type | Animals | Entropy slope | Entropy q | Off-program slope | Off-program q |
|---|---:|---:|---:|---:|---:|
| Lung fibroblast | 10 | +0.003010 | 0.5313 | +0.00002168 | 0.7954 |
| Bronchial smooth-muscle cell | 13 | -0.002166 | 0.1895 | +0.00000165 | 0.7954 |
| Alveolar macrophage | 14 | -0.001378 | 0.7664 | -0.00000009 | 0.9750 |
| Classical monocyte | 14 | +0.000720 | 0.7895 | +0.00000533 | 0.7954 |
| Intermediate monocyte | 8 | -0.000355 | 0.9464 | +0.00000603 | 0.7954 |
| B cell | 13 | +0.000115 | 0.9464 | +0.00000183 | 0.7979 |
| NK cell | 15 | -0.000297 | 0.9464 | -0.00000152 | 0.8213 |
| CD8-positive alpha-beta T cell | 10 | +0.003379 | 0.6392 | +0.00000123 | 0.7979 |
| Vein endothelial cell | 8 | -0.005760 | 0.1895 | not evaluable | not evaluable |

The vein endothelial reference produced no genes meeting the off-program definition, so only entropy was tested for that cell type. One CD8-positive T-cell animal also lacked an evaluable off-program reference, leaving nine animals for that endpoint.

## Sensitivity and interpretation

With exactly 20 cells per animal, fibroblast entropy increased with age (slope +0.004692 per month, q = 0.0487). The same analysis gave negative entropy slopes for bronchial smooth-muscle cells (-0.003136, q = 0.0464) and vein endothelial cells (-0.005370, q = 0.0464). The primary all-cell fibroblast slope was smaller and not significant. This dependence on cell-number balancing prevents a simple claim that fibroblast entropy increases with age.

No off-program result survived FDR correction in the primary, balanced-cell, pooled-identifier-exclusion, or stringent gene-exclusion analyses. Excluding pooled animal identifiers also did not produce a significant primary entropy increase.

The earlier distance-to-median analysis addressed a different outcome: variability between cells. Its nominal fibroblast increase also failed to reproduce under without-replacement UMI sampling. Taken together, the analyses do not provide consistent transcriptomic support for age-related loss of transcriptional control in these mouse lung cell types. A direct test of the proposed mechanism would require chromatin measurements and experimental manipulation of DNA damage or repair, not transcriptomics alone.

Full methods are in `notes/transcriptional-disorder-methods.md`. Frequentist cell-type estimates are in `results/h3-transcriptional-disorder.csv`. The combined Bayesian and frequentist global comparison is in `report/transcriptional-disorder-bayesian.md`.
