# Whole-expression-vector displacement with age in mouse, human, and naked mole-rat

## Bottom line

The broad three-species hypothesis is not supported in its simplest form. Mouse shows a clear overall increase in distance from young expression profiles. Human does not show a significant overall increase, although several neuronal and epithelial cell types have significant positive or negative associations. NMR has a negative overall point estimate in two skin compartments, but five animals cannot establish expression stability or a significant between-species difference.

This analysis measures age-related displacement of complete expression vectors. It does not distinguish coordinated biological changes from random transcriptional noise, and it does not test DNA damage or chromatin restoration.

## Species-wide primary results

The primary endpoint is Hellinger distance between each cell's complete within-cell expression-proportion vector and a cross-fitted young reference. No UMI downsampling is used. Effects are in cell-type standard deviations per 10 percent of species maximum lifespan.

| Species | Tested units | Donors or animals | Effect | 95% CI | p value | Interpretation |
|---|---:|---:|---:|---:|---:|---|
| Mouse | 62 tissue-by-cell-type units | 23 | +0.1028 | +0.0678 to +0.1378 | 0.00000398 | Overall increase detected |
| Human | 135 cohort-by-cell-type units | 479 cohort-specific donors | +0.0365 | -0.0221 to +0.0951 | 0.2218 | No overall increase detected |
| Naked mole-rat | 2 skin compartments | 5 | -0.4744 | not estimated | exact p = 0.2 | Negative point estimate, unresolved |

Mouse and human models use unit fixed intercepts, sex, and donor-clustered t-based inference. The NMR value is the mean standardized aged-minus-young contrast rescaled from the observed 10.2-year gap to 3.1 years. Its exact test uses all ten assignments of two aged animals among five. The NMR effect is not directly comparable in precision to the regression-based mouse and human estimates.

## Human cohorts are individually unresolved

No human cohort-level primary summary has a significant age association.

| Cohort | Units | Donors | Effect per 10% lifespan | p value |
|---|---:|---:|---:|---:|
| MSSM DLPFC | 20 | 212 | +0.0894 | 0.0920 |
| HBCC DLPFC | 20 | 88 | -0.0007 | 0.9908 |
| Kidney living donor | 6 | 19 | +0.0096 | 0.9666 |
| Kidney SCR | 9 | 18 | +0.2076 | 0.5336 |
| Lung CellRef | 17 | 17 | -0.0182 | 0.8783 |
| Lung Integrated | 34 | 113 | +0.0367 | 0.5562 |
| Immune Aging BAL | 29 | 12 | -0.0404 | 0.7328 |

RADC DLPFC, SEA-AD DLPFC, Kidney Normal Adult, and Kidney SN lacked enough adult donors aged 20 to 40 for the prespecified young reference and were not tested.

## Cell-type results

Five mouse tissue-by-cell-type units pass Benjamini-Hochberg correction in the primary analysis. All five are positive.

| Tissue | Cell type | Animals | Slope per year | q value |
|---|---|---:|---:|---:|
| Bladder | bladder cell | 8 | +0.037168 | 0.0126 |
| Kidney | fenestrated cell | 14 | +0.008615 | 0.0314 |
| Kidney | podocyte | 12 | +0.021591 | 0.0126 |
| Limb muscle | mesenchymal stem cell | 16 | +0.022389 | 0.0059 |
| Lung | fibroblast of lung | 10 | +0.010948 | 0.0059 |

Six human cohort-by-cell-type units pass correction within their cohorts. Four are positive neuronal associations in MSSM and two are negative epithelial associations in integrated lung.

| Cohort | Cell type | Donors | Slope per year | q value |
|---|---|---:|---:|---:|
| MSSM DLPFC | GABAergic neuron | 155 | +0.000990 | 0.000071 |
| MSSM DLPFC | VIP GABAergic cortical interneuron | 197 | +0.000694 | 0.0037 |
| MSSM DLPFC | pvalb GABAergic cortical interneuron | 189 | +0.000566 | 0.0317 |
| MSSM DLPFC | sst GABAergic cortical interneuron | 177 | +0.000558 | 0.0108 |
| Lung Integrated | acinar cell | 16 | -0.005213 | 0.0035 |
| Lung Integrated | nasal mucosa goblet cell | 37 | -0.002443 | 0.0443 |

These cell-type results do not make the overall human summary significant. Follow-up analysis showed that all six lose significance after UMI-depth adjustment. Age-associated depth has the sign needed to create positive MSSM and negative lung associations. The unadjusted hits are therefore not established biological findings.

Normalized donor-pseudobulk analysis gives a different result. No gene survives depth-adjusted FDR in the two lung epithelial units, but MSSM interneurons have extensive depth-adjusted gene-expression changes. Twenty-four VIP, 20 pvalb, and 62 sst genes pass FDR in both MSSM and HBCC. These replicated interneuron changes support targeted neuronal biology rather than universal whole-vector disruption. See `report/human-expression-vector-interpretation.md` and `report/cross-species-bayesian-vector.md`.

## NMR compartment results

Both evaluable NMR compartments have negative aged-minus-young Hellinger distances.

| Compartment | Animals | Aged minus young distance | Standardized aged-minus-young effect | Exact p value |
|---|---:|---:|---:|---:|
| Basal | 5 | -0.009793 | -1.2731 | 0.2 |
| Corneum | 5 | -0.007285 | -1.8489 | 0.2 |

The direction is opposite to the primary mouse result. The exact permutation distribution has only ten assignments, so the minimum possible two-sided p value is 0.1. These data do not establish that NMR profiles are stable or that the NMR effect differs significantly from mouse or human effects.

## Sensitivity analyses

The mouse overall result remains positive after mean-log-UMI adjustment (effect +0.0679, p = 0.00000430) and after excluding pooled animal identifiers (unadjusted effect +0.1051, p = 0.00000867). The depth-adjusted human estimate is +0.0231 with p = 0.0503. This is a sensitivity result, not grounds to replace the non-significant primary human analysis.

The literal raw-count Euclidean analysis gives an overall mouse effect of +0.0563 (p = 0.0266), a human effect of -0.0121 (p = 0.6665), and an NMR effect of -0.2395 (exact p = 0.4). Raw Euclidean distance is strongly influenced by UMI depth and is not a pure biological endpoint.

## Interpretation

The analysis supports a general mouse atlas-wide association, not a universal mammalian association. Human aging shows cell-type-specific positive and negative vector displacement without an overall increase. NMR shows a reverse-direction point estimate in the available skin data, but the experiment is too small for a stability claim. The proposed mouse-and-human-versus-NMR pattern is therefore not observed as stated.

Full methods are in `notes/cross-species-expression-vector-methods.md`. Unit-level results are in `results/h3-cross-species-vector-tests.csv`, and species summaries are in `results/h3-cross-species-vector-global.csv`.
