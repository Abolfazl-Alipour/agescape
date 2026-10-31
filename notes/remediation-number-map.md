# Remediation number map

This note is a verification record, not manuscript prose. Values were checked against committed result tables on 2026-10-31.

## Abstract number mapping

| Abstract claim or value | Committed source | Row | Column |
|---|---|---|---|
| 11 public human cohorts | `results/analysis-cohort-summary.csv` | `m9_cohorts` | `value` |
| 621 scored control donors | `results/analysis-cohort-summary.csv` | `m9_scored_control_donors` | `value` |
| 284 PsychAD control donors | `results/analysis-cohort-summary.csv` | `psychad_control_donors` | `value` |
| 1,000 methylation profiles | `results/analysis-cohort-summary.csv` | `methylation_array_profiles` | `value` |
| 139 primary genes, 257 sensitivity genes | `results/signature-tier-summary.csv` | `P1_Clean`, `P0_Raw257` | `n_mouse_genes` |
| Primary whole-cortex slope +0.000186, p = 0.3494 | `results/m10-cohort-summary.csv` | `DLPFC_Lifespan_Atlas` | `slope`, `p` |
| Adjusted whole-cortex slope -0.000209, p = 0.4004 | `results/m10-cohort-summary.csv` | `DLPFC_Lifespan_Atlas` | `slope_covar_adj`, `p_covar_adj` |
| 252 primary tests and 3 survivors | `results/analysis-cohort-summary.csv` | `primary_human_transport_tests`, `primary_human_transport_survivors` | `value` |
| MSSM L2/3 neuron slope -0.001303, q = 0.0208 | `results/human-transport-global-multiplicity.csv` | `transport_022` | `slope_per_year`, `q_global_bh` |
| MSSM oligodendrocyte slope +0.001152, q = 0.0346 | `results/human-transport-global-multiplicity.csv` | `transport_023` | `slope_per_year`, `q_global_bh` |
| PsychAD EN_L3_5_IT_3 slope -0.001383, q = 0.0208 | `results/human-transport-global-multiplicity.csv` | `transport_230` | `slope_per_year`, `q_global_bh` |
| Fibroblast Cdkn1a slope +0.1080, p = 0.0029, q = 0.0322 | `results/fibroblast-markers-audit.csv` | `Cdkn1a` | `slope_per_month`, `p_value`, `q_value_bh` |
| Human club-cell slope -0.0602, q = 0.0044 | `results/h2-composition-human-lung.csv` | `club cell` | `slope_per_year_or_month`, `q_bh` |
| Human multiciliated-cell slope -0.0342, q = 0.0225 | `results/h2-composition-human-lung.csv` | `multiciliated epithelial cell` | `slope_per_year_or_month`, `q_bh` |
| 44 human lung cell types | `results/analysis-cohort-summary.csv` | `human_lung_tested_cell_types` | `value` |
| Multiciliated detection caveat | `results/h2-composition-human-lung-presence.csv` | `multiciliated epithelial cell` | `presence_age_r`, `presence_age_p`, `q_bh` |
| Club-cell study-adjusted p = 0.1033 | `results/d3-human-composition-study-sensitivity.csv` | `club cell` | `study_adj_p` |
| 3 of 7 concordant, binomial p = 1.000, rho = 0.000 | `results/h2-composition-cross-species-test.csv` | first row | `n_sign_concordant`, `n_pairs`, `binomial_two_sided_p`, `spearman_rho` |
| Mouse clock MAE 0.36 years, predictive R-squared 0.3169 | `results/h1-epigenetic-clock-mouse.csv` | `Mouse`, `Liver` | `mae_years`, `predictive_r2` |
| Human clock MAE 5.06 years, predictive R-squared 0.6383, 661 profiles | `results/h1-epigenetic-clock-human.csv` | `Human`, `All Tissues` | `mae_years`, `predictive_r2`, `n_samples` |
| 724 of 816 mouse clock probes | `results/analysis-cohort-summary.csv` | `mouse_clock_probes_used`, `universal_clock_probes` | `value` |
| Five NMR animals and no detected slope | `results/m6-statistical-summary.csv` | `Naked mole-rat skin aging slope (n=5 animals)` | `n`, `slope`, `p` |

No abstract number remains without a committed CSV source.

## Main before-and-after values

| Analysis | Starting value at 3d02af8 | Remediated value |
|---|---|---|
| Human adult lung composition input | 117 donors, 25 cell types, thresholded rows zero-filled | 119 donors, 44 cell types, complete pre-threshold counts |
| Human lung FDR-significant types | 4 | 2 |
| Club-cell association | slope -0.067799, q = 0.0123 | slope -0.060227, q = 0.004423 |
| Multiciliated epithelial association | not in the previous 25-type table | slope -0.034154, q = 0.022513 |
| Club-cell study-adjusted result | slope -0.015436, p = 0.2415 | slope -0.015400, p = 0.103282 |
| Cross-species composition correlation | rho = -0.286, p = 0.535 | rho = 0.000, p = 1.000 |
| MSSM whole-cohort score | slope -0.000193, p = 0.02825, q = 0.3107 | slope -0.000404, p = 0.240302, q = 0.940170 |
| PsychAD primary whole-cortex model | P0 slope +0.000318, p = 0.01182 | P1 slope +0.000186, p = 0.349415 |
| PsychAD primary adjusted whole-cortex model | P0 slope +0.000027, p = 0.862622 | P1 slope -0.000209, p = 0.400449 |
| Primary human multiplicity | M9-only 225-test family and selected glial emphasis | one 252-test M9/M10 family, all survivors exported |
| Global survivors | SEA-AD OPC plus two neuronal results under the former scoring | MSSM L2/3 neuron, MSSM oligodendrocyte, PsychAD EN_L3_5_IT_3 |
| PsychAD microglia pooled adjusted | P0 slope +0.001417, p = 0.002022 | P1 slope +0.002069, p = 0.001496, global q = 0.062818 |
| PsychAD microglia adult stratum | P0 slope +0.001997, p = 0.001919, q = 0.051813 | P1 slope +0.003115, p = 0.000617, adult q = 0.016655 |
| Mouse clock R-squared | 0.851, squared Pearson correlation | 0.3169, predictive R-squared from per-sample predictions |
| Human clock R-squared | 0.642, squared Pearson correlation | 0.6383, predictive R-squared from per-sample predictions |
| Without-replacement noise family | 8 cell types | 9 cell types, intermediate monocyte restored |
| Intermediate monocyte noise | absent from sensitivity CSV | slope +0.000092, p = 0.710142, q = 0.798910 |

## Overlap audit

The audit compares only age, sex, and bank because direct identifiers are not shared. It found candidate matches for 160 of 214 public MSSM donors and 133 of 134 public HBCC donors. Fifty-two MSSM and 35 HBCC donors had unique candidates on these fields. These are not confirmed identity matches. Source: `results/m9-m10-donor-overlap-audit.csv`.

## Literature and data-use evidence

- Tyshkovskiy et al. full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC13233323/. Inspected abstract, study design, conserved modules, and data availability.
- Karakaslar et al. full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC10086525/. Inspected abstract, AP-1 mouse immune findings, human chromatin evidence, TMS analysis, and discussion.
- Lu et al. full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC10501909/. Inspected consortium training description, cross-validation distinction, and data availability.
- AD Knowledge Portal current data-use terms: https://adknowledgeportal.synapse.org/DataAccess/. Rendered and read in a browser. See `notes/adkp-data-use-verification.md`.

## Source discrepancy corrected during verification

The remediation brief described `P1_Clean` as excluding sex-dimorphic and immediate-early genes. The signature-construction code shows that P1 applies predicted-gene, pseudogene, biotype, ribosomal, mitochondrial, and ambient-RNA filters. Sex trimming occurs in P2 and immediate-early trimming occurs in P3. P1 remains the approved primary tier, but the manuscript now states its actual filters and does not claim independence from sex-dimorphic or AP-1-related immediate-early genes.
