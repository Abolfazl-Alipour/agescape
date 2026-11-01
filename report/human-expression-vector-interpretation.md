# Why human expression-vector associations vary by cohort and cell type

## Bottom line

The current human results do not establish tissue-specific biological aging. No human cohort-level age association is significant. All six cell-type associations that pass FDR in the unadjusted analysis disappear after UMI-depth adjustment. Age-associated UMI depth has the sign needed to create positive associations in MSSM neurons and negative associations in integrated-lung epithelial cells. Depth, statistical power, and assay composition are therefore central explanations. Cell-type-specific biology remains possible, but it is not demonstrated by the unadjusted distance result.

## Depth explains the opposite directions

Hellinger distance is negatively associated with mean log UMI in all six unadjusted hits. Lower-depth cells tend to be farther from the young reference.

| Cohort and cell type | Depth change with age, per year | Depth-age p | Primary age slope | Depth-adjusted age slope | Depth-adjusted q |
|---|---:|---:|---:|---:|---:|
| MSSM GABAergic neuron | -0.009826 | 0.000007 | +0.000990 | +0.000036 | 0.4693 |
| MSSM VIP GABAergic interneuron | -0.006600 | 0.000772 | +0.000694 | +0.000044 | 0.2355 |
| MSSM pvalb GABAergic interneuron | -0.005959 | 0.005160 | +0.000566 | -0.000006 | 0.9062 |
| MSSM sst GABAergic interneuron | -0.005574 | 0.001922 | +0.000558 | +0.000017 | 0.7688 |
| Lung acinar cell | +0.048526 | 0.000275 | -0.005213 | -0.001015 | 0.8061 |
| Lung nasal mucosa goblet cell | +0.032337 | 0.000258 | -0.002443 | +0.000225 | 0.9778 |

In MSSM, older donors have lower UMI depth, which increases the distance measurement and produces a positive apparent aging association. In the lung epithelial units, older donors have higher UMI depth, which decreases the distance measurement and produces a negative apparent association. These opposite directions are exactly what sampling-depth confounding predicts.

Depth adjustment cannot determine whether the depth trend is caused by sequencing and sample processing or by a real change in cellular RNA content. It does show that the unadjusted distance associations are not independent evidence of biological expression disruption.

## Brain-cohort comparison

HBCC has the same positive unadjusted direction for all four GABAergic populations, but its estimates are smaller and non-significant. HBCC also has fewer evaluable donors per unit: 71 to 82 compared with 155 to 197 in MSSM. After depth adjustment, HBCC slopes are also approximately zero.

Direct cohort-by-age interaction tests do not establish different MSSM and HBCC slopes. Interaction p values are 0.0976 for GABAergic neurons, 0.2137 for VIP, 0.2492 for pvalb, and 0.8642 for sst. MSSM-only significance is therefore not evidence of a true MSSM-specific biological effect. Donor count and depth remain plausible explanations.

MSSM and HBCC both use 10x 3' v3 nuclei. Mixing different assay chemistries does not explain their depth trends. Sequencing depth, nuclear RNA content, and sample processing remain unresolved alternatives.

## Lung assay comparison

Integrated lung mixes seven assays and includes both cell and nucleus suspensions. The significant acinar analysis includes two assays and the goblet analysis includes five. Adding assay and suspension indicators to the depth-adjusted models does not recover a significant age association.

| Cell type | Age slope after depth and assay adjustment | p value |
|---|---:|---:|
| Acinar cell | -0.001064 | 0.1017 |
| Nasal mucosa goblet cell | +0.000633 | 0.1082 |

Both models are full rank. The results do not support an assay-independent age association in these unadjusted hits.

## Power and annotation

MSSM has the largest brain donor set, which increases precision. Kidney cohorts have 18 or 19 donors, Immune Aging BAL has 12, and Lung CellRef has 17. Several cohort and cell-type tests are therefore imprecise. A significant result in a larger cohort and a non-significant result in a smaller cohort do not by themselves demonstrate biological heterogeneity.

The significant lung labels are available in the integrated atlas but lack equivalent age-spanning tests in the other lung cohorts. Their replication cannot be evaluated from the current data. Broad cell-type changes and changes within annotated subtypes are also different biological questions.

## Normalized gene-level results

Donor-pseudobulk log-CP10K analysis separates gene-expression change from the distance artifact. No gene survives depth-adjusted FDR in acinar or nasal mucosa goblet cells. Acinar cells have 115 significant genes before depth adjustment, including FOS, DUSP1, and SOCS3, but none afterward. The lung results therefore lack normalized gene-level support in these models.

MSSM neurons show a different pattern. Depth-adjusted analysis retains 5 FDR-significant genes in broad GABAergic neurons, 294 in VIP, 333 in pvalb, and 123 in sst interneurons. Seven genes pass FDR in all three interneuron subtypes with consistent directions: GABRB1, RNLS, STAG1, and TTC28 increase, while CAMK2N1, OLFM1, and RNF152 decrease. Other repeated changes include CACNB2, SCN2B, MKRN3, and genes involved in neuronal signaling. This is evidence of coordinated transcriptional change within MSSM, not merely an increased distance summary.

The interneuron changes replicate in HBCC. Among MSSM-significant genes, 272 of 294 VIP genes, 295 of 333 pvalb genes, and 118 of 123 sst genes have the same direction in HBCC. Twenty-four VIP, 20 pvalb, and 62 sst genes pass FDR in both cohorts. All seven genes shared across the three MSSM subtypes have the same direction in HBCC. RNLS increases and passes FDR in all three HBCC subtypes. GABRB1 passes FDR in HBCC VIP and sst cells. These are reproducible cell-type-specific transcriptional associations, although they do not establish a causal mechanism or show that neurons are the only aging human cell type.

Gene-level summaries are in `results/h3-human-vector-gene-driver-summary.csv`, shared neuronal genes are in `results/h3-human-mssm-shared-neuronal-genes.csv`, and independent-cohort results are in `results/h3-human-neuronal-gene-replication-summary.csv`.

## Interpretation

The defensible conclusion has two parts. Human expression-vector distance is strongly depth-dependent, and apparent tissue or cell-type specificity can arise when depth also varies with age. Separately, interneurons show normalized gene-expression changes that replicate across MSSM and HBCC. The first result invalidates a simple biological reading of the unadjusted distance hits. The second supports targeted investigation of neuronal signaling in human brain aging, not a universal tissue-wide expression-disruption claim.

Design results are in `results/h3-human-vector-design-audit.csv`, brain-cohort comparisons are in `results/h3-human-brain-cohort-replication.csv`, and assay models are in `results/h3-human-lung-assay-audit.csv`.
