# Organ cell-fraction analysis

## Question

The test asks whether the fractions of cell types within an organ change with age. It does not test fibrosis, differential expression, or a universal aging program.

## Mouse

Each Tabula Muris Senis organ was analyzed separately. For each mouse and cell type, the input was the number of cells divided by the total analyzed cells in that mouse and organ. Cell types observed in at least six mice were tested. A centered log-ratio model with a 0.5 pseudocount was fitted against age and sex. Sex was dropped only when the design was rank-deficient. Benjamini-Hochberg correction was applied within each organ. A sensitivity analysis excluded pooled mouse identifiers.

## Human

Each public cohort was analyzed separately because tissue processing and cell-type annotations differ across cohorts. Donors younger than 20 years were excluded. Cell types observed in at least six adult donors were tested with the same centered log-ratio model used in mouse. Complete pre-threshold counts were used. A second model restricted each cell type to donors with at least one observed cell of that type.

## Naked mole-rat

The public GSE152886 matrix has no cell-type labels. Cells were re-annotated using the paper's processing outline and marker families. Sample 2V was excluded as in the paper. Cells required at least 2,000 UMIs and 700 detected genes. Expression was log-normalized, 2,000 highly variable genes were retained, and Leiden clustering used resolution 0.4. Immune markers were Cd207, Ctss, and Mfge8. Basal markers were Krt5, Trp63, and Krt15. Corneum markers were Krt10 and Krt1. Remaining keratinocyte clusters were labeled intermediate. Annotation was repeated at resolution 0.6 as a sensitivity check.

The test compared mean aged and young fractions with all ten possible assignments of two aged animals among five. This is an exact two-sided permutation test. The minimum possible p value is 0.1, so no NMR result can reach FDR 0.05 in this dataset. Non-significance must not be described as evidence that NMR composition is stable.

## Outputs

- `results/organ-composition-age-tests.csv`: all estimates, models, sensitivity labels, and q values.
- `results/organ-composition-summary.csv`: one primary summary per mouse organ, human cohort, and NMR skin.
- `results/nmr-skin-cluster-annotation.csv`: cluster-to-compartment assignments and marker means.
- `results/nmr-skin-cell-fractions-resolution-0.4.csv` and `results/nmr-skin-cell-fractions-resolution-0.6.csv`: animal-level NMR fractions.
