# PBMC expression analysis

Analysis base revision at plan freeze: `5e8bc6fff4fbb372c80c89338fe2c1ab6c929cb2`. The statistical plan was frozen in `data/manifest/pbmc_analysis_plan.json` before expression download.

## Scope and main verdict

Only **GSE213516 (Zhu)** and **GSE214546 (Thomson)** were analyzed, sequentially and independently. GSE157007 remained descriptive/excluded. No matrices were concatenated across cohorts, no p-values were pooled across studies, and no cell-level inferential test was used. Strict young `<=40` versus old `>60` is primary; atlas age `>=60` and exact-age models are sensitivities.

- **Robustness warning—Thomson strict P1:** primary 11/13, depth 8/13, CMV 11/13, pool 9/13, but batch 0/13 cell types have q<0.05 (minimum batch q=0.109444). **Broad P1 is not robust to batch and must not be presented as confirmed aging biology.**
- **Composition replicated in direction:** CD8-naive T-cell CLR decreased in both cohorts, with FDR below 0.05 in both. Non-classical-monocyte CLR increased in both, but only Thomson passed FDR. Targeted Thomson batch estimates also remain directionally consistent, but both have exploratory all-label batch-family q=0.139044; no favorable targeted batch q is claimed.
- **CD8-naive gene replication was incomplete:** all six targeted effects were positive in Zhu; CLIC1, S100A4 and CCL5 passed its six-test BH family. Four were positive in Thomson, but none passed its BH family. PRF1 and GZMA were negative in Thomson.
- **P1 was cohort-dependent:** no Zhu strict-primary OLS test passed BH, whereas 11 Thomson labels passed. Zhu matched-null P1 survived for four labels; Thomson matched-null survived for 12. These are not pooled or interpreted as a universal effect.
- **Hellinger did not show a universal positive displacement:** three Zhu cell types passed strict depth-adjusted BH, all with negative old-minus-young effects; no Thomson cell type passed. Exploratory Thomson batch/CMV families do have survivors, but they do not reverse the primary conclusion or establish robust aging biology. Non-significance is not equivalence.
- **MAIT was untestable** in both released annotation universes. No integrated labels were imported to create it.

All numerical effects, confidence intervals and q values are in the linked aggregate CSVs. Tables below include every strict-primary P1 and Hellinger cell type, every targeted effect, all targeted atlas/exact-age sensitivities, Thomson targeted covariate sensitivities, and all ineligible endpoints. The full exploratory pseudobulk tables are gzip-compressed CSV rather than truncated hit lists.

## Input validation and label joins

| Cohort | Donors downloaded | Raw cells | Annotation rows | Matched rows | Unmatched annotation rows excluded | Counts | Joins |
|---|---:|---:|---:|---:|---:|---|---|
| GSE213516 | 17 | 139,761 | 131,800 | 131,601 | 199 | finite nonnegative integers | one-to-one |
| GSE214546 | 32 | 650,166 | 650,166 | 650,166 | 0 | finite nonnegative integers | one-to-one |

All 17 Zhu and 32 Thomson per-donor files passed full stored-value validation and byte-size checks. SHA-256 values and per-file dimensions are in `results/followup-pbmc-GSE213516-validation.json` and `results/followup-pbmc-GSE214546-validation.json`. Duplicate gene-symbol rows were excluded from P1 and targeted gene tests rather than artificially renamed.

Zhu annotation aliases were resolved by exact full barcode matches: HR→F40, HYD→F65, XMZ→F60 and ZY→F46. F55 uses `baseline_<barcode>-F55`; the other direct donors use `<donor>_<barcode>`. Across Zhu, 199 annotation-only rows were explicitly excluded; every retained annotation maps uniquely to a raw barcode.

**Zhu genome-build adjudication:** GEO says hg19 and the primary Methods say GRCh38. The actual released feature files contain Ensembl IDs, symbols and feature types, but no reference-genome identity; MatrixMarket metadata identifies Cell Ranger 3.0.0 without resolving the assembly. The conflict is therefore unresolved from public records. Zhu gene results are **symbol-level only**, with no coordinate, locus or assembly-specific claim.

Released-label inventories and the explicit common-label map are in `results/followup-pbmc-label-harmonization.csv`. Only CD8-naive and CD16/non-classical monocytes support the targeted cross-cohort comparisons. MAIT is absent/non-harmonizable in both.

## Strict targeted donor-pseudobulk replication

Effects are old-minus-young differences in donor log2CPM. BH is over six genes separately within each cohort and model variant.

| cohort | cell_type | gene | Effect [95% CI] | p | q | N | Status/reason |
|---|---|---|---|---|---|---|---|
| GSE213516 | CD8-Naive | CCL5 | 1.823 [0.567, 3.08] | 0.01096 | 0.02191 | 10 | fit |
| GSE214546 | CD8 Naive | CCL5 | 1.059 [0.1511, 1.967] | 0.02447 | 0.08544 | 23 | fit |
| GSE213516 | CD8-Naive | CLIC1 | 0.9544 [0.4969, 1.412] | 0.001689 | 0.005068 | 10 | fit |
| GSE214546 | CD8 Naive | CLIC1 | 0.231 [-0.08374, 0.5457] | 0.1415 | 0.2122 | 23 | fit |
| GSE213516 | CD8-Naive | GZMA | 0.8503 [-0.8697, 2.57] | 0.2807 | 0.2807 | 10 | fit |
| GSE214546 | CD8 Naive | GZMA | -0.165 [-0.84, 0.51] | 0.6156 | 0.7388 | 23 | fit |
| GSE213516 | CD8-Naive | PRF1 | 0.6979 [0.03135, 1.364] | 0.04247 | 0.05096 | 10 | fit |
| GSE214546 | CD8 Naive | PRF1 | -0.02647 [-0.4304, 0.3774] | 0.8926 | 0.8926 | 23 | fit |
| GSE213516 | CD8-Naive | PTGER2 | 1.029 [0.07664, 1.981] | 0.03783 | 0.05096 | 10 | fit |
| GSE214546 | CD8 Naive | PTGER2 | 0.5893 [0.06867, 1.11] | 0.02848 | 0.08544 | 23 | fit |
| GSE213516 | CD8-Naive | S100A4 | 0.9891 [0.5834, 1.395] | 0.0006879 | 0.004127 | 10 | fit |
| GSE214546 | CD8 Naive | S100A4 | 0.2576 [-0.07497, 0.5903] | 0.1218 | 0.2122 | 23 | fit |

### Thomson CMV/batch/pool/depth targeted sensitivities

All 23 strict donors supported full-rank sex, CMV, B-batch, P-pool and depth models at the donor level (9 residual df for the 12-level B-batch model). Each sensitivity is a separate six-test BH family.

| cell_type | gene | model | Effect [95% CI] | p | q | N | Status/reason |
|---|---|---|---|---|---|---|---|
| CD8 Naive | CCL5 | batch | 1.228 [-0.5515, 3.007] | 0.153 | 0.2295 | 23 | fit |
| CD8 Naive | CCL5 | cmv | 1.075 [0.05468, 2.095] | 0.03996 | 0.1101 | 23 | fit |
| CD8 Naive | CCL5 | depth | 1.033 [0.1288, 1.937] | 0.0273 | 0.1015 | 23 | fit |
| CD8 Naive | CCL5 | pool | 1.123 [0.1881, 2.058] | 0.02111 | 0.1266 | 23 | fit |
| CD8 Naive | CLIC1 | batch | 0.4397 [0.01209, 0.8674] | 0.04503 | 0.09007 | 23 | fit |
| CD8 Naive | CLIC1 | cmv | 0.3089 [-0.03227, 0.6501] | 0.0734 | 0.1101 | 23 | fit |
| CD8 Naive | CLIC1 | depth | 0.2332 [-0.09058, 0.557] | 0.1481 | 0.2222 | 23 | fit |
| CD8 Naive | CLIC1 | pool | 0.2435 [-0.08423, 0.5712] | 0.1364 | 0.2046 | 23 | fit |
| CD8 Naive | GZMA | batch | 0.2266 [-1.104, 1.557] | 0.709 | 0.709 | 23 | fit |
| CD8 Naive | GZMA | cmv | -0.1107 [-0.8666, 0.6451] | 0.7624 | 0.8335 | 23 | fit |
| CD8 Naive | GZMA | depth | -0.1857 [-0.8549, 0.4835] | 0.5682 | 0.6819 | 23 | fit |
| CD8 Naive | GZMA | pool | -0.07752 [-0.7447, 0.5897] | 0.8105 | 0.8105 | 23 | fit |
| CD8 Naive | PRF1 | batch | -0.1213 [-0.5885, 0.3459] | 0.5713 | 0.6856 | 23 | fit |
| CD8 Naive | PRF1 | cmv | 0.04539 [-0.4003, 0.4911] | 0.8335 | 0.8335 | 23 | fit |
| CD8 Naive | PRF1 | depth | -0.01727 [-0.4249, 0.3903] | 0.9303 | 0.9303 | 23 | fit |
| CD8 Naive | PRF1 | pool | 0.06316 [-0.2861, 0.4125] | 0.7093 | 0.8105 | 23 | fit |
| CD8 Naive | PTGER2 | batch | 0.9849 [0.1899, 1.78] | 0.02062 | 0.06187 | 23 | fit |
| CD8 Naive | PTGER2 | cmv | 0.7939 [0.2625, 1.325] | 0.005549 | 0.0333 | 23 | fit |
| CD8 Naive | PTGER2 | depth | 0.5844 [0.04962, 1.119] | 0.03383 | 0.1015 | 23 | fit |
| CD8 Naive | PTGER2 | pool | 0.5535 [0.01697, 1.09] | 0.04383 | 0.1315 | 23 | fit |
| CD8 Naive | S100A4 | batch | 0.4297 [0.1921, 0.6673] | 0.002715 | 0.01629 | 23 | fit |
| CD8 Naive | S100A4 | cmv | 0.3372 [-0.02431, 0.6987] | 0.0658 | 0.1101 | 23 | fit |
| CD8 Naive | S100A4 | depth | 0.2473 [-0.08211, 0.5768] | 0.1326 | 0.2222 | 23 | fit |
| CD8 Naive | S100A4 | pool | 0.2858 [-0.0542, 0.6258] | 0.0946 | 0.1892 | 23 | fit |

## Strict targeted CLR composition replication

Effects are old-minus-young donor CLR changes with a 0.5-cell pseudocount. Only eligible tests enter the three-target BH family; MAIT remains an explicit untestable row.

| cohort | cell_type | common_cell_type | Effect [95% CI] | p | q | N | Status/reason |
|---|---|---|---|---|---|---|---|
| GSE213516 | CD8-Naive | CD8 naive T | -1.606 [-2.628, -0.5833] | 0.007517 | 0.01503 | 10 | fit |
| GSE214546 | CD8 Naive | CD8 naive T | -1.424 [-1.981, -0.8665] | 3.236e-05 | 6.472e-05 | 23 | fit |
| GSE213516 | MAIT | MAIT | — | — | — | — | untestable: released_label_absent_or_not_harmonizable |
| GSE214546 | MAIT | MAIT | — | — | — | — | untestable: released_label_absent_or_not_harmonizable |
| GSE213516 | CD16-MC | non-classical monocytes | 0.5402 [-0.1071, 1.188] | 0.08905 | 0.08905 | 10 | fit |
| GSE214546 | CD16+ Monocytes | non-classical monocytes | 0.9213 [0.4274, 1.415] | 0.0009082 | 0.0009082 | 23 | fit |

### Thomson targeted CLR covariate robustness

These estimates come directly from `results/followup-pbmc-GSE214546-clr-models.csv`. Depth means log total annotated cells. The q values below are the existing exploratory all-label families, not a recomputed two-target or otherwise favorable targeted family. In particular, batch gives CD8 Naive -1.396839, 95% CI [-2.466337, -0.327341], and CD16+ Monocytes +1.040266, 95% CI [0.111162, 1.969369]; both have exploratory all-label batch-family q=0.139044. CMV, pool and depth results are shown without substituting their values into the primary targeted table.

| cell_type | model | Effect [95% CI] | p | q | N | Status/reason |
|---|---|---|---|---|---|---|
| CD16+ Monocytes | batch | 1.04 [0.1112, 1.969] | 0.03209 | 0.139 | 23 | fit |
| CD16+ Monocytes | cmv | 1.085 [0.566, 1.605] | 0.0003264 | 0.002122 | 23 | fit |
| CD16+ Monocytes | log_total_cells | 0.9174 [0.4125, 1.422] | 0.001202 | 0.007813 | 23 | fit |
| CD16+ Monocytes | pool | 0.8603 [0.3694, 1.351] | 0.001636 | 0.01063 | 23 | fit |
| CD8 Naive | batch | -1.397 [-2.466, -0.3273] | 0.0161 | 0.139 | 23 | fit |
| CD8 Naive | cmv | -1.423 [-2.05, -0.7972] | 0.0001371 | 0.001782 | 23 | fit |
| CD8 Naive | log_total_cells | -1.428 [-1.998, -0.8581] | 4.62e-05 | 0.0006005 | 23 | fit |
| CD8 Naive | pool | -1.45 [-2.029, -0.8704] | 4.683e-05 | 0.0006088 | 23 | fit |

### Thomson strict P1 robustness counts

Primary has 11/13, depth 8/13, CMV 11/13 and pool 9/13 cell types with q<0.05, but batch has 0/13 (minimum batch q=0.109444). Thus the broad P1 pattern is not robust to batch and is not confirmed aging biology.

| Model | Cell types with q<0.05 (of 13) | Minimum q |
|---|---:|---:|
| primary | 11 | 0.022422 |
| depth | 8 | 0.003517 |
| CMV | 11 | 0.014672 |
| pool | 9 | 0.032773 |
| batch | 0 | 0.109444 |

## Complete strict-primary P1 results

P1 uses per-cell log1p(CPM), arithmetic donor means and old-minus-young OLS with sex. `q` is the OLS BH value. The matched-null q values are in the aggregate P1 table and use 200 expression-decile-matched draws, seed 20261114.

| cohort | cell_type | Effect [95% CI] | p | q | N | Status/reason |
|---|---|---|---|---|---|---|
| GSE214546 | B cell progenitor | 0.0763 [0.0136, 0.139] | 0.01955 | 0.02938 | 23 | fit |
| GSE214546 | CD14+ Monocytes | 0.1029 [0.02919, 0.1766] | 0.008624 | 0.02242 | 23 | fit |
| GSE213516 | CD14-MC | -0.0168 [-0.1752, 0.1416] | 0.8092 | 0.9588 | 10 | fit |
| GSE213516 | CD14-MC-PPBP | 0.04307 [-0.05599, 0.1421] | 0.3381 | 0.9588 | 10 | fit |
| GSE214546 | CD16+ Monocytes | 0.08603 [0.02117, 0.1509] | 0.0119 | 0.02578 | 23 | fit |
| GSE213516 | CD16-MC | -0.03615 [-0.1593, 0.08697] | 0.5099 | 0.9588 | 10 | fit |
| GSE214546 | CD4 Memory | 0.0868 [0.03005, 0.1435] | 0.004593 | 0.02242 | 23 | fit |
| GSE214546 | CD4 Naive | 0.08114 [0.0309, 0.1314] | 0.003051 | 0.02242 | 23 | fit |
| GSE213516 | CD4-Naive | 0.004824 [-0.2081, 0.2178] | 0.9588 | 0.9588 | 10 | fit |
| GSE213516 | CD4-Tem | -0.006892 [-0.1621, 0.1483] | 0.9193 | 0.9588 | 10 | fit |
| GSE213516 | CD4-Tm | 0.005724 [-0.1823, 0.1937] | 0.9446 | 0.9588 | 10 | fit |
| GSE213516 | CD4-Treg | 0.04366 [-0.0918, 0.1791] | 0.4451 | 0.9588 | 8 | fit |
| GSE214546 | CD8 Naive | 0.08878 [0.02897, 0.1486] | 0.005694 | 0.02242 | 23 | fit |
| GSE214546 | CD8 effector | 0.08487 [0.01462, 0.1551] | 0.02034 | 0.02938 | 23 | fit |
| GSE213516 | CD8-CTL | 0.006569 [-0.1678, 0.1809] | 0.9315 | 0.9588 | 10 | fit |
| GSE213516 | CD8-Naive | 0.006419 [-0.1788, 0.1917] | 0.937 | 0.9588 | 10 | fit |
| GSE213516 | CD8-Tem | -0.01347 [-0.1719, 0.1449] | 0.8463 | 0.9588 | 10 | fit |
| GSE214546 | Dendritic cell | 0.08293 [0.02397, 0.1419] | 0.008206 | 0.02242 | 23 | fit |
| GSE214546 | Double negative T cell | 0.0899 [0.0175, 0.1623] | 0.01749 | 0.02938 | 23 | fit |
| GSE213516 | Intermed-MC | 0.04714 [-0.1637, 0.258] | 0.6134 | 0.9588 | 10 | fit |
| GSE213516 | Megakaryocytes | -0.02411 [-0.1648, 0.1166] | 0.6973 | 0.9588 | 10 | fit |
| GSE213516 | Memory-B-CRIP1 | -0.01975 [-0.2311, 0.1916] | 0.8314 | 0.9588 | 10 | fit |
| GSE213516 | Memory-B-IL32 | -0.09351 [-0.2521, 0.0651] | 0.1992 | 0.9588 | 9 | fit |
| GSE214546 | NK cell | 0.07999 [0.006213, 0.1538] | 0.03502 | 0.04138 | 23 | fit |
| GSE213516 | NK-FCER1G | -0.05327 [-0.2467, 0.1402] | 0.5358 | 0.9588 | 10 | fit |
| GSE213516 | NK-GZMH | -0.01278 [-0.1914, 0.1658] | 0.8704 | 0.9588 | 10 | fit |
| GSE213516 | NK-GZMK | -0.01695 [-0.2136, 0.1797] | 0.8443 | 0.9588 | 10 | fit |
| GSE213516 | NK-IL7R | -0.04177 [-0.2344, 0.1509] | 0.6148 | 0.9588 | 9 | fit |
| GSE213516 | NK-S100A2 | -0.04193 [-0.2023, 0.1184] | 0.5312 | 0.9588 | 8 | fit |
| GSE213516 | NK-SELL | -0.109 [-0.2688, 0.05085] | 0.1401 | 0.9588 | 8 | fit |
| GSE213516 | Naive-B | -0.02527 [-0.2087, 0.1582] | 0.7541 | 0.9588 | 10 | fit |
| GSE214546 | Platelets | -0.03407 [-0.07974, 0.0116] | 0.1349 | 0.1349 | 22 | fit |
| GSE213516 | cDC | -0.1086 [-0.2357, 0.01843] | 0.08137 | 0.9588 | 9 | fit |
| GSE213516 | pDC | -0.1933 [-0.4686, 0.08188] | 0.1114 | 0.9588 | 6 | fit |
| GSE214546 | pDC | 0.07742 [0.007284, 0.1476] | 0.03217 | 0.04138 | 23 | fit |
| GSE214546 | pre-B cell | 0.05714 [-0.006158, 0.1204] | 0.07433 | 0.08052 | 23 | fit |

### Strict P1 matched-null results

| cohort | cell_type | Effect [95% CI] | p | q | N | Status/reason |
|---|---|---|---|---|---|---|
| GSE214546 | B cell progenitor | 0.0763 [0.0136, 0.139] | 0.004975 | 0.00588 | 23 | fit |
| GSE214546 | CD14+ Monocytes | 0.1029 [0.02919, 0.1766] | 0.004975 | 0.00588 | 23 | fit |
| GSE213516 | CD14-MC | -0.0168 [-0.1752, 0.1416] | 0.1443 | 0.237 | 10 | fit |
| GSE213516 | CD14-MC-PPBP | 0.04307 [-0.05599, 0.1421] | 0.5871 | 0.6751 | 10 | fit |
| GSE214546 | CD16+ Monocytes | 0.08603 [0.02117, 0.1509] | 0.004975 | 0.00588 | 23 | fit |
| GSE213516 | CD16-MC | -0.03615 [-0.1593, 0.08697] | 0.05473 | 0.1399 | 10 | fit |
| GSE214546 | CD4 Memory | 0.0868 [0.03005, 0.1435] | 0.004975 | 0.00588 | 23 | fit |
| GSE214546 | CD4 Naive | 0.08114 [0.0309, 0.1314] | 0.004975 | 0.00588 | 23 | fit |
| GSE213516 | CD4-Naive | 0.004824 [-0.2081, 0.2178] | 0.5124 | 0.6203 | 10 | fit |
| GSE213516 | CD4-Tem | -0.006892 [-0.1621, 0.1483] | 0.1294 | 0.2289 | 10 | fit |
| GSE213516 | CD4-Tm | 0.005724 [-0.1823, 0.1937] | 0.2139 | 0.3029 | 10 | fit |
| GSE213516 | CD4-Treg | 0.04366 [-0.0918, 0.1791] | 0.6269 | 0.6866 | 8 | fit |
| GSE214546 | CD8 Naive | 0.08878 [0.02897, 0.1486] | 0.004975 | 0.00588 | 23 | fit |
| GSE214546 | CD8 effector | 0.08487 [0.01462, 0.1551] | 0.004975 | 0.00588 | 23 | fit |
| GSE213516 | CD8-CTL | 0.006569 [-0.1678, 0.1809] | 0.2736 | 0.3496 | 10 | fit |
| GSE213516 | CD8-Naive | 0.006419 [-0.1788, 0.1917] | 0.7861 | 0.7861 | 10 | fit |
| GSE213516 | CD8-Tem | -0.01347 [-0.1719, 0.1449] | 0.06468 | 0.1488 | 10 | fit |
| GSE214546 | Dendritic cell | 0.08293 [0.02397, 0.1419] | 0.004975 | 0.00588 | 23 | fit |
| GSE214546 | Double negative T cell | 0.0899 [0.0175, 0.1623] | 0.004975 | 0.00588 | 23 | fit |
| GSE213516 | Intermed-MC | 0.04714 [-0.1637, 0.258] | 0.7313 | 0.7646 | 10 | fit |
| GSE213516 | Megakaryocytes | -0.02411 [-0.1648, 0.1166] | 0.2239 | 0.3029 | 10 | fit |
| GSE213516 | Memory-B-CRIP1 | -0.01975 [-0.2311, 0.1916] | 0.03483 | 0.1144 | 10 | fit |
| GSE213516 | Memory-B-IL32 | -0.09351 [-0.2521, 0.0651] | 0.004975 | 0.02861 | 9 | fit |
| GSE214546 | NK cell | 0.07999 [0.006213, 0.1538] | 0.004975 | 0.00588 | 23 | fit |
| GSE213516 | NK-FCER1G | -0.05327 [-0.2467, 0.1402] | 0.01493 | 0.06866 | 10 | fit |
| GSE213516 | NK-GZMH | -0.01278 [-0.1914, 0.1658] | 0.1194 | 0.2289 | 10 | fit |
| GSE213516 | NK-GZMK | -0.01695 [-0.2136, 0.1797] | 0.07463 | 0.156 | 10 | fit |
| GSE213516 | NK-IL7R | -0.04177 [-0.2344, 0.1509] | 0.0398 | 0.1144 | 9 | fit |
| GSE213516 | NK-S100A2 | -0.04193 [-0.2023, 0.1184] | 0.2139 | 0.3029 | 8 | fit |
| GSE213516 | NK-SELL | -0.109 [-0.2688, 0.05085] | 0.004975 | 0.02861 | 8 | fit |
| GSE213516 | Naive-B | -0.02527 [-0.2087, 0.1582] | 0.0199 | 0.07629 | 10 | fit |
| GSE214546 | Platelets | -0.03407 [-0.07974, 0.0116] | 0.3433 | 0.3433 | 22 | fit |
| GSE213516 | cDC | -0.1086 [-0.2357, 0.01843] | 0.004975 | 0.02861 | 9 | fit |
| GSE213516 | pDC | -0.1933 [-0.4686, 0.08188] | 0.004975 | 0.02861 | 6 | fit |
| GSE214546 | pDC | 0.07742 [0.007284, 0.1476] | 0.01493 | 0.01617 | 23 | fit |
| GSE214546 | pre-B cell | 0.05714 [-0.006158, 0.1204] | 0.004975 | 0.00588 | 23 | fit |

## Complete strict depth-adjusted Hellinger results

Distances use actual raw UMI >=1500, technical-gene exclusions, strict-young cross-fitted references and donor mean per-cell distances. Negative effects mean old cells were closer to the young reference; they do not indicate a universal increase in displacement.

| cohort | cell_type | Effect [95% CI] | p | q | N | Status/reason |
|---|---|---|---|---|---|---|
| GSE214546 | B cell progenitor | 0.000569 [-0.003913, 0.005051] | 0.7933 | 0.7933 | 23 | fit |
| GSE214546 | CD14+ Monocytes | -0.002794 [-0.006114, 0.0005269] | 0.09434 | 0.2264 | 23 | fit |
| GSE213516 | CD14-MC | -0.00243 [-0.008524, 0.003664] | 0.3669 | 0.3844 | 10 | fit |
| GSE213516 | CD14-MC-PPBP | -0.003831 [-0.009344, 0.001681] | 0.1399 | 0.1951 | 10 | fit |
| GSE214546 | CD16+ Monocytes | -0.00367 [-0.008413, 0.001072] | 0.1218 | 0.2435 | 23 | fit |
| GSE213516 | CD16-MC | -0.004818 [-0.01197, 0.002337] | 0.1505 | 0.1951 | 10 | fit |
| GSE214546 | CD4 Memory | -0.003853 [-0.007101, -0.0006061] | 0.0225 | 0.103 | 23 | fit |
| GSE214546 | CD4 Naive | -0.002261 [-0.005766, 0.001245] | 0.1929 | 0.2647 | 23 | fit |
| GSE213516 | CD4-Naive | -0.004057 [-0.01215, 0.00404] | 0.2661 | 0.3081 | 10 | fit |
| GSE213516 | CD4-Tem | -0.004772 [-0.01186, 0.002319] | 0.1507 | 0.1951 | 10 | fit |
| GSE213516 | CD4-Tm | -0.006022 [-0.01358, 0.001536] | 0.09912 | 0.1677 | 10 | fit |
| GSE213516 | CD4-Treg | -0.01175 [-0.02284, -0.000655] | 0.04237 | 0.1036 | 8 | fit |
| GSE214546 | CD8 Naive | -0.002215 [-0.005592, 0.001163] | 0.1859 | 0.2647 | 23 | fit |
| GSE214546 | CD8 effector | -0.004873 [-0.008761, -0.0009843] | 0.01675 | 0.103 | 23 | fit |
| GSE213516 | CD8-CTL | -0.009703 [-0.01623, -0.003176] | 0.01087 | 0.05306 | 10 | fit |
| GSE213516 | CD8-Naive | -0.00318 [-0.01038, 0.00402] | 0.3213 | 0.3534 | 10 | fit |
| GSE213516 | CD8-Tem | -0.008345 [-0.01271, -0.003982] | 0.003395 | 0.0249 | 10 | fit |
| GSE214546 | Dendritic cell | -0.002154 [-0.007519, 0.003211] | 0.4111 | 0.4485 | 23 | fit |
| GSE214546 | Double negative T cell | -0.003687 [-0.006876, -0.0004972] | 0.02575 | 0.103 | 23 | fit |
| GSE213516 | Intermed-MC | -0.005115 [-0.02003, 0.009798] | 0.4335 | 0.4335 | 10 | fit |
| GSE213516 | Megakaryocytes | -0.008908 [-0.0261, 0.008289] | 0.2519 | 0.3079 | 10 | fit |
| GSE213516 | Memory-B-CRIP1 | -0.004661 [-0.01114, 0.001815] | 0.1287 | 0.1951 | 10 | fit |
| GSE213516 | Memory-B-IL32 | -0.01089 [-0.02201, 0.000241] | 0.05351 | 0.1177 | 9 | fit |
| GSE214546 | NK cell | -0.001428 [-0.00497, 0.002113] | 0.409 | 0.4485 | 23 | fit |
| GSE213516 | NK-FCER1G | -0.008723 [-0.01473, -0.002712] | 0.01206 | 0.05306 | 10 | fit |
| GSE213516 | NK-GZMH | -0.01039 [-0.01902, -0.001773] | 0.0256 | 0.07041 | 10 | fit |
| GSE213516 | NK-GZMK | -0.008626 [-0.01194, -0.005316] | 0.0006995 | 0.01539 | 10 | fit |
| GSE213516 | NK-IL7R | -0.01872 [-0.03251, -0.004932] | 0.01747 | 0.06407 | 9 | fit |
| GSE213516 | NK-S100A2 | -0.005094 [-0.01132, 0.001131] | 0.08552 | 0.1568 | 8 | fit |
| GSE213516 | NK-SELL | -0.009562 [-0.01681, -0.00231] | 0.02156 | 0.06777 | 8 | fit |
| GSE213516 | Naive-B | -0.004726 [-0.00974, 0.0002888] | 0.0606 | 0.1212 | 10 | fit |
| GSE214546 | Platelets | — | — | — | 13 | untestable: fewer_than_three_donors_per_group |
| GSE213516 | cDC | -0.005523 [-0.008053, -0.002994] | 0.002483 | 0.0249 | 9 | fit |
| GSE214546 | pDC | -0.002737 [-0.007036, 0.001563] | 0.1985 | 0.2647 | 23 | fit |
| GSE214546 | pre-B cell | -0.003168 [-0.006882, 0.0005459] | 0.09018 | 0.2264 | 23 | fit |

### Thomson Hellinger covariate sensitivity caveat

The strict depth-adjusted primary family has no survivors (minimum q=0.102988). The exploratory batch-depth family has 2 survivors and the exploratory CMV-depth family has 1; pool-depth has none. These sensitivity-family survivors have negative effects and do not reverse the primary result or establish robust aging biology.

| cell_type | model | Effect [95% CI] | p | q | N | Status/reason |
|---|---|---|---|---|---|---|
| CD4 Memory | batch_depth | -0.005363 [-0.008442, -0.002284] | 0.003859 | 0.04631 | 23 | fit |
| Double negative T cell | batch_depth | -0.005832 [-0.009684, -0.00198] | 0.008182 | 0.04909 | 23 | fit |
| Double negative T cell | cmv_depth | -0.004955 [-0.008105, -0.001805] | 0.003941 | 0.04729 | 23 | fit |

## Targeted atlas-age and exact-age sensitivities

Atlas groups add age-60 donors; exact-age models use all 17 Zhu and all 32 Thomson adults with continuous age. These are separate BH families, not replacements for the strict primary contrast.

| cohort | design | cell_type | gene | Effect [95% CI] | p | q | N | Status/reason |
|---|---|---|---|---|---|---|---|---|
| GSE213516 | atlas_group | CD8-Naive |  | -1.21 [-2.386, -0.03516] | 0.04476 | 0.08952 | 12 | fit |
| GSE214546 | atlas_group | CD8 Naive |  | -1.27 [-1.788, -0.7523] | 4.271e-05 | 8.543e-05 | 25 | fit |
| GSE213516 | atlas_group | CD8-Naive | CCL5 | 1.161 [-0.4394, 2.761] | 0.1352 | 0.2029 | 12 | fit |
| GSE214546 | atlas_group | CD8 Naive | CCL5 | 0.9835 [0.1682, 1.799] | 0.02029 | 0.06086 | 25 | fit |
| GSE213516 | atlas_group | CD8-Naive | CLIC1 | 0.717 [0.03432, 1.4] | 0.04151 | 0.1245 | 12 | fit |
| GSE214546 | atlas_group | CD8 Naive | CLIC1 | 0.3137 [0.03082, 0.5967] | 0.03132 | 0.06265 | 25 | fit |
| GSE213516 | atlas_group | CD8-Naive | GZMA | 0.4254 [-1.157, 2.008] | 0.5582 | 0.5582 | 12 | fit |
| GSE214546 | atlas_group | CD8 Naive | GZMA | -0.0933 [-0.6923, 0.5057] | 0.7497 | 0.8997 | 25 | fit |
| GSE213516 | atlas_group | CD8-Naive | PRF1 | 0.4395 [-0.3047, 1.184] | 0.2144 | 0.2573 | 12 | fit |
| GSE214546 | atlas_group | CD8 Naive | PRF1 | 0.01882 [-0.3306, 0.3682] | 0.9121 | 0.9121 | 25 | fit |
| GSE213516 | atlas_group | CD8-Naive | PTGER2 | 0.7841 [-0.2061, 1.774] | 0.1069 | 0.2029 | 12 | fit |
| GSE214546 | atlas_group | CD8 Naive | PTGER2 | 0.6221 [0.1522, 1.092] | 0.01181 | 0.06086 | 25 | fit |
| GSE213516 | atlas_group | CD8-Naive | S100A4 | 0.7553 [0.1435, 1.367] | 0.02096 | 0.1245 | 12 | fit |
| GSE214546 | atlas_group | CD8 Naive | S100A4 | 0.2632 [-0.03022, 0.5566] | 0.07627 | 0.1144 | 25 | fit |
| GSE213516 | atlas_group | MAIT |  | — | — | — | — | untestable: released_label_absent_or_not_harmonizable |
| GSE214546 | atlas_group | MAIT |  | — | — | — | — | untestable: released_label_absent_or_not_harmonizable |
| GSE213516 | atlas_group | CD16-MC |  | 0.4435 [-0.1291, 1.016] | 0.1136 | 0.1136 | 12 | fit |
| GSE214546 | atlas_group | CD16+ Monocytes |  | 0.7289 [0.2687, 1.189] | 0.003385 | 0.003385 | 25 | fit |
| GSE213516 | exact_age | CD8-Naive |  | -0.03713 [-0.06694, -0.007324] | 0.01823 | 0.03647 | 17 | fit |
| GSE214546 | exact_age | CD8 Naive |  | -0.03594 [-0.0505, -0.02137] | 2.227e-05 | 4.455e-05 | 32 | fit |
| GSE213516 | exact_age | CD8-Naive | CCL5 | 0.04411 [0.009647, 0.07856] | 0.01579 | 0.02369 | 17 | fit |
| GSE214546 | exact_age | CD8 Naive | CCL5 | 0.02432 [0.001954, 0.04668] | 0.0341 | 0.06198 | 32 | fit |
| GSE213516 | exact_age | CD8-Naive | CLIC1 | 0.02332 [0.008076, 0.03856] | 0.005463 | 0.01093 | 17 | fit |
| GSE214546 | exact_age | CD8 Naive | CLIC1 | 0.009962 [0.002279, 0.01764] | 0.01284 | 0.03852 | 32 | fit |
| GSE213516 | exact_age | CD8-Naive | GZMA | 0.0222 [-0.0231, 0.06751] | 0.311 | 0.311 | 17 | fit |
| GSE214546 | exact_age | CD8 Naive | GZMA | -0.0009608 [-0.0176, 0.01568] | 0.9068 | 0.9068 | 32 | fit |
| GSE213516 | exact_age | CD8-Naive | PRF1 | 0.01745 [-0.004594, 0.03949] | 0.1116 | 0.134 | 17 | fit |
| GSE214546 | exact_age | CD8 Naive | PRF1 | -0.003129 [-0.01333, 0.007076] | 0.5355 | 0.6426 | 32 | fit |
| GSE213516 | exact_age | CD8-Naive | PTGER2 | 0.03155 [0.0118, 0.0513] | 0.004095 | 0.01093 | 17 | fit |
| GSE214546 | exact_age | CD8 Naive | PTGER2 | 0.02202 [0.008548, 0.03548] | 0.002295 | 0.01377 | 32 | fit |
| GSE213516 | exact_age | CD8-Naive | S100A4 | 0.02595 [0.01419, 0.0377] | 0.0003196 | 0.001917 | 17 | fit |
| GSE214546 | exact_age | CD8 Naive | S100A4 | 0.007988 [0.0003366, 0.01564] | 0.04132 | 0.06198 | 32 | fit |
| GSE213516 | exact_age | MAIT |  | — | — | — | — | untestable: released_label_absent_or_not_harmonizable |
| GSE214546 | exact_age | MAIT |  | — | — | — | — | untestable: released_label_absent_or_not_harmonizable |
| GSE213516 | exact_age | CD16-MC |  | 0.01013 [-0.01259, 0.03285] | 0.3553 | 0.3553 | 17 | fit |
| GSE214546 | exact_age | CD16+ Monocytes |  | 0.0237 [0.01247, 0.03493] | 0.000169 | 0.000169 | 32 | fit |

The complete P1, exploratory CLR and Hellinger atlas/exact-age tables are in `results/followup-pbmc-replication-sensitivities.csv`. Thomson P1 and Hellinger CMV/batch/pool/depth results are in the cohort-specific model CSVs.

## Exploratory donor-pseudobulk families

All features detected in at least 20% of eligible donors were tested without effect-based prefiltering. BH was applied separately by cohort, design, model and cell type. Full effects/CIs/q values are in `results/followup-pbmc-GSE213516-pseudobulk-exploratory.csv.gz` and `results/followup-pbmc-GSE214546-pseudobulk-exploratory.csv.gz`; all q<0.05 rows are repeated in `results/followup-pbmc-pseudobulk-exploratory-survivors.csv.gz` for convenience.

| Cohort | Design | Model | Tested genes | q<0.05 |
|---|---|---|---:|---:|
| GSE213516 | atlas_group | primary | 337,321 | 2 |
| GSE213516 | exact_age | primary | 364,030 | 1 |
| GSE213516 | strict_group | depth | 344,557 | 0 |
| GSE213516 | strict_group | primary | 344,557 | 8 |
| GSE214546 | atlas_group | primary | 235,170 | 3,672 |
| GSE214546 | exact_age | primary | 231,278 | 3,912 |
| GSE214546 | strict_group | depth | 232,716 | 3,739 |
| GSE214546 | strict_group | primary | 232,716 | 3,632 |

Exploratory CLR results, including every eligible released label and raw-fraction sensitivities, are in `results/followup-pbmc-clr-models.csv`.

## Ineligible and untestable endpoints

Eligibility failures are not null effects. Generic failures below include cell types failing the 20-cell/three-donor/group/six-donor/reference floors. Model-specific untestable rows follow separately.

| cohort | endpoint | cell_type | design | reason | Effect [95% CI] | p | q | N | Status/reason |
|---|---|---|---|---|---|---|---|---|---|
| GSE213516 | eligibility | ABC | strict_group | donor_or_group_floor | — | — | — | 4 | untestable: donor_or_group_floor |
| GSE213516 | eligibility | ABC | atlas_group | donor_or_group_floor | — | — | — | 5 | untestable: donor_or_group_floor |
| GSE213516 | eligibility | CD4-Naive-NECTIN1 | strict_group | donor_or_group_floor | — | — | — | 0 | untestable: donor_or_group_floor |
| GSE213516 | eligibility | CD4-Naive-NECTIN1 | atlas_group | donor_or_group_floor | — | — | — | 0 | untestable: donor_or_group_floor |
| GSE213516 | eligibility | CD4-Naive-NECTIN1 | exact_age | donor_or_group_floor | — | — | — | 0 | untestable: donor_or_group_floor |
| GSE213516 | eligibility | Naive-B-NRGN | strict_group | donor_or_group_floor | — | — | — | 3 | untestable: donor_or_group_floor |
| GSE213516 | eligibility | Naive-B-NRGN | atlas_group | donor_or_group_floor | — | — | — | 4 | untestable: donor_or_group_floor |
| GSE213516 | eligibility | Naive-B-NRGN | exact_age | donor_or_group_floor | — | — | — | 5 | untestable: donor_or_group_floor |
| GSE213516 | eligibility | Plasma | strict_group | donor_or_group_floor | — | — | — | 4 | untestable: donor_or_group_floor |
| GSE213516 | eligibility | Plasma | atlas_group | donor_or_group_floor | — | — | — | 5 | untestable: donor_or_group_floor |
| GSE213516 | eligibility | RBC | strict_group | donor_or_group_floor | — | — | — | 1 | untestable: donor_or_group_floor |
| GSE213516 | eligibility | RBC | atlas_group | donor_or_group_floor | — | — | — | 1 | untestable: donor_or_group_floor |
| GSE213516 | eligibility | RBC | exact_age | donor_or_group_floor | — | — | — | 1 | untestable: donor_or_group_floor |
| GSE213516 | Hellinger_eligibility | ABC | strict_group | young_reference_or_scored_donor_floor | — | — | — | — | untestable: young_reference_or_scored_donor_floor |
| GSE213516 | Hellinger_eligibility | ABC | atlas_group | young_reference_or_scored_donor_floor | — | — | — | — | untestable: young_reference_or_scored_donor_floor |
| GSE213516 | Hellinger_eligibility | CD4-Naive-NECTIN1 | strict_group | young_reference_or_scored_donor_floor | — | — | — | — | untestable: young_reference_or_scored_donor_floor |
| GSE213516 | Hellinger_eligibility | CD4-Naive-NECTIN1 | atlas_group | young_reference_or_scored_donor_floor | — | — | — | — | untestable: young_reference_or_scored_donor_floor |
| GSE213516 | Hellinger_eligibility | CD4-Naive-NECTIN1 | exact_age | young_reference_or_scored_donor_floor | — | — | — | — | untestable: young_reference_or_scored_donor_floor |
| GSE213516 | Hellinger_eligibility | Naive-B-NRGN | strict_group | young_reference_or_scored_donor_floor | — | — | — | — | untestable: young_reference_or_scored_donor_floor |
| GSE213516 | Hellinger_eligibility | Naive-B-NRGN | atlas_group | young_reference_or_scored_donor_floor | — | — | — | — | untestable: young_reference_or_scored_donor_floor |
| GSE213516 | Hellinger_eligibility | Naive-B-NRGN | exact_age | young_reference_or_scored_donor_floor | — | — | — | — | untestable: young_reference_or_scored_donor_floor |
| GSE213516 | Hellinger_eligibility | Plasma | strict_group | young_reference_or_scored_donor_floor | — | — | — | — | untestable: young_reference_or_scored_donor_floor |
| GSE213516 | Hellinger_eligibility | Plasma | atlas_group | young_reference_or_scored_donor_floor | — | — | — | — | untestable: young_reference_or_scored_donor_floor |
| GSE213516 | Hellinger_eligibility | Plasma | exact_age | young_reference_or_scored_donor_floor | — | — | — | — | untestable: young_reference_or_scored_donor_floor |
| GSE213516 | Hellinger_eligibility | RBC | strict_group | young_reference_or_scored_donor_floor | — | — | — | — | untestable: young_reference_or_scored_donor_floor |
| GSE213516 | Hellinger_eligibility | RBC | atlas_group | young_reference_or_scored_donor_floor | — | — | — | — | untestable: young_reference_or_scored_donor_floor |
| GSE213516 | Hellinger_eligibility | RBC | exact_age | young_reference_or_scored_donor_floor | — | — | — | — | untestable: young_reference_or_scored_donor_floor |
| GSE213516 | Hellinger_eligibility | pDC | strict_group | young_reference_or_scored_donor_floor | — | — | — | — | untestable: young_reference_or_scored_donor_floor |

### Model-specific untestable rows

| cohort | endpoint | cell_type | design | model | reason | Effect [95% CI] | p | q | N | Status/reason |
|---|---|---|---|---|---|---|---|---|---|---|
| GSE213516 | CLR_targeted | MAIT | strict_group | primary | released_label_absent_or_not_harmonizable | — | — | — | — | untestable: released_label_absent_or_not_harmonizable |
| GSE213516 | CLR_targeted | MAIT | atlas_group | primary | released_label_absent_or_not_harmonizable | — | — | — | — | untestable: released_label_absent_or_not_harmonizable |
| GSE213516 | CLR_targeted | MAIT | exact_age | primary | released_label_absent_or_not_harmonizable | — | — | — | — | untestable: released_label_absent_or_not_harmonizable |
| GSE214546 | CLR_targeted | MAIT | strict_group | primary | released_label_absent_or_not_harmonizable | — | — | — | — | untestable: released_label_absent_or_not_harmonizable |
| GSE214546 | CLR_targeted | MAIT | atlas_group | primary | released_label_absent_or_not_harmonizable | — | — | — | — | untestable: released_label_absent_or_not_harmonizable |
| GSE214546 | CLR_targeted | MAIT | exact_age | primary | released_label_absent_or_not_harmonizable | — | — | — | — | untestable: released_label_absent_or_not_harmonizable |
| GSE213516 | Hellinger | pDC | atlas_group | primary | fewer_than_three_donors_per_group | — | — | — | 7 | untestable: fewer_than_three_donors_per_group |
| GSE214546 | Hellinger | Platelets | strict_group | primary | fewer_than_three_donors_per_group | — | — | — | 13 | untestable: fewer_than_three_donors_per_group |
| GSE214546 | Hellinger | Platelets | strict_group | unadjusted | fewer_than_three_donors_per_group | — | — | — | 13 | untestable: fewer_than_three_donors_per_group |
| GSE214546 | Hellinger | Platelets | strict_group | cmv_depth | fewer_than_three_donors_per_group | — | — | — | 13 | untestable: fewer_than_three_donors_per_group |
| GSE214546 | Hellinger | Platelets | strict_group | batch_depth | rank_or_residual_degrees_of_freedom | — | — | — | 13 | untestable: rank_or_residual_degrees_of_freedom |
| GSE214546 | Hellinger | Platelets | strict_group | pool_depth | fewer_than_three_donors_per_group | — | — | — | 13 | untestable: fewer_than_three_donors_per_group |

## Covariate identifiability

The complete rank table is `results/followup-pbmc-covariate-identifiability.csv`. At the donor level, every strict/atlas/exact-age Thomson primary, depth, CMV, B-batch and P-pool candidate was full rank. Cell-type-specific rank failures remain untestable, notably the strict Thomson platelet Hellinger model. CMV and pool adjustments are reported separately rather than hidden in a pooled model.

## Resource use

Every expression download and analysis/validation/summarization job used `scripts/resource_guard.py` with a 30 GB disk reserve, 3 GB per-file cap, 8 GB memory cap, swap disabled in the job cgroup and one heavy-job lock. No RAW tar was downloaded. Zhu completed before any Thomson download.

| Cohort | Download files | Download bytes | Analysis seconds | Peak RSS GiB | Donor-cache bytes |
|---|---:|---:|---:|---:|---:|
| GSE213516 | 52 | 895,410,109 | 47.2 | 1.740 | 192,593,066 |
| GSE214546 | 32 | 2,238,830,714 | 55.9 | 1.764 | 171,681,370 |

## Interpretation and limitations

The strict composition result supports a reproducible CD8-naive decrease and supports, less definitively across studies, increased non-classical-monocyte abundance. Their Thomson batch estimates remain directional but have exploratory all-label batch q=0.139044, so batch robustness is not established. The six-gene targeted program did not replicate uniformly at FDR. Broad Thomson P1 is not robust to batch (0/13 survivors) and is not confirmed aging biology. Hellinger effects were predominantly negative and did not establish greater old-cell displacement; exploratory batch/CMV survivors do not change that primary conclusion. Atlas-age and exact-age sensitivities must be read alongside the strict results rather than used to select a preferred design.

Zhu has only five donors per strict group and an unresolved assembly record; its symbol-level results should not be used for locus-specific claims. Thomson CMV/batch/pool models have limited residual degrees of freedom and observational covariate confounding. Released label definitions differ, and MAIT cannot be tested. No result establishes longitudinal change, causality or equivalence from non-significance.

## Reproduction

```bash
python scripts/fetch_followup_pbmc.py --cohort GSE213516
python scripts/resource_guard.py run --path ../agescape-revision-data --min-free-gb 30 --memory-gb 8 -- python scripts/run_followup_pbmc.py --cohort GSE213516
python scripts/fetch_followup_pbmc.py --cohort GSE214546
python scripts/resource_guard.py run --path ../agescape-revision-data --min-free-gb 30 --memory-gb 8 -- python scripts/run_followup_pbmc.py --cohort GSE214546
python scripts/resource_guard.py run --path ../agescape-revision-data --min-free-gb 30 --memory-gb 8 -- python scripts/summarize_followup_pbmc.py
python scripts/resource_guard.py run --path ../agescape-revision-data --min-free-gb 30 --memory-gb 8 -- python scripts/report_followup_pbmc.py
```

The fetch script invokes `resource_guard.py download` separately for every per-donor file. Raw data and donor caches remain outside Git under `../agescape-revision-data/pbmc/`.
