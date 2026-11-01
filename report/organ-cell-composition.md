# Do cell-type fractions change with age?

## Bottom line

Cell-type fractions change with age in many, but not all, tested organs. Eight of 13 adequately covered mouse organs had at least one cell-type fraction associated with age at FDR 0.05. Three of 11 human cohort datasets had such a change. Naked mole-rat skin had no significant change, but five animals provide only ten possible age-label assignments and cannot produce a two-sided permutation p value below 0.1. The NMR result is inconclusive, not evidence of stability.

## Results by organ or cohort

| Species | Organ or cohort | Donors | Tested cell types | FDR-significant types | Interpretation |
|---|---|---:|---:|---:|---|
| Mouse | Bladder | 8 | 4 | 1 | fraction change detected |
| Mouse | Fat | 6 | 5 | 0 | no change detected |
| Mouse | Heart and aorta | 11 | 8 | 0 | no change detected |
| Mouse | Kidney | 16 | 19 | 4 | fraction change detected |
| Mouse | Large intestine | 3 | 0 | 0 | insufficient donor coverage |
| Mouse | Limb muscle | 16 | 9 | 1 | fraction change detected |
| Mouse | Liver | 12 | 8 | 5 | fraction change detected |
| Mouse | Lung | 16 | 29 | 3 | fraction change detected |
| Mouse | Mammary gland | 6 | 7 | 0 | no change detected |
| Mouse | Marrow | 17 | 18 | 2 | fraction change detected |
| Mouse | Pancreas | 8 | 8 | 0 | no change detected |
| Mouse | Skin | 3 | 0 | 0 | insufficient donor coverage |
| Mouse | Spleen | 13 | 12 | 1 | fraction change detected |
| Mouse | Thymus | 6 | 6 | 0 | no change detected |
| Mouse | Tongue | 9 | 3 | 2 | fraction change detected |
| Mouse | Trachea | 2 | 0 | 0 | insufficient donor coverage |
| Human | HBCC DLPFC | 88 | 24 | 0 | no change detected |
| Human | MSSM DLPFC | 213 | 25 | 6 | fraction change detected |
| Human | RADC DLPFC | 18 | 23 | 0 | no change detected |
| Human | SEA-AD DLPFC | 21 | 18 | 0 | no change detected |
| Human | Kidney living donor | 19 | 17 | 1 | fraction change detected |
| Human | Kidney normal adult | 18 | 26 | 0 | no change detected |
| Human | Kidney SCR v15 | 18 | 23 | 0 | no change detected |
| Human | Kidney SN v15 | 9 | 22 | 0 | no change detected |
| Human | Immune-aging BAL | 12 | 31 | 0 | no change detected |
| Human | Lung CellRef | 17 | 32 | 0 | no change detected |
| Human | Lung integrated | 119 | 50 | 2 | fraction change detected |
| Naked mole-rat | Skin | 5 | 4 | 0 | inconclusive due to permutation floor |

## Main changes

Mouse kidney showed decreased capillary endothelial and mesangial fractions and an increased plasma-cell fraction. Mouse liver showed increased B-cell and myeloid-leukocyte fractions and a decreased hepatocyte fraction. Mouse lung showed an increased T-cell fraction and a decreased pulmonary interstitial fibroblast fraction. These are compositional associations, not evidence of fibrosis or organ destruction.

Human MSSM cortex showed six changing fractions, including increased oligodendrocytes and decreased oligodendrocyte precursor cells, GABAergic neurons, and SST interneurons. Human integrated lung showed decreased club-cell and multiciliated epithelial-cell fractions. The living-donor kidney cohort showed an increased kidney cortex collecting-duct epithelial-cell fraction.

## Sensitivity

Excluding pooled mouse identifiers retained most mouse findings, but liver NK cells and lung vein endothelial cells were no longer significant. Human MSSM findings remained significant with observed records only. Human integrated lung and living-donor kidney findings did not survive that sensitivity check. These human results should be read as cohort-specific and, outside MSSM, exploratory.

NMR basal, corneum, and immune compartments all had small age-group differences with permutation p values above 0.1. Basal direction changed between annotation resolutions 0.4 and 0.6. Only 46 cells were assigned to the immune compartment at the primary resolution. The NMR data do not establish either stability or change.

## Methods

Donor or animal was the replicate. Mouse and human tests used centered log-ratio fractions with age and, when identifiable, sex. Cell types required observations in at least six donors. FDR correction was applied within each organ or cohort. NMR cells were re-annotated from published marker families because GEO supplies no labels. Exact permutations compared aged and young fractions. Full methods are in `notes/organ-composition-methods.md`. All estimates are in `results/organ-composition-age-tests.csv`.
