# Cell-Type-Specific Transcriptional Signatures, Lung Compositional Change, and Epigenetic Clock Benchmarking Across Mammalian Aging

**Authors:** Abolfazl (Reza) Alipour  
**Correspondence:** `reza@agescape.org` | `https://rezaalipour.com/agescape/`  
**Target Repository:** `https://github.com/Abolfazl-Alipour/agescape`

---

## Abstract

We tested whether a mouse-derived transcriptional aging signature shows consistent age associations across human cohorts, cell types, and a naked mole-rat comparison. The analysis used Tabula Muris Senis, 621 control donors across 11 public human cohorts, 284 control donors in the PsychAD dorsolateral prefrontal cortex lifespan atlas, and 1,000 mammalian methylation array profiles. The primary signature contained 139 genes after technical, gene-biotype, ribosomal, mitochondrial, and ambient-RNA filters. The unfiltered 257-gene signature was retained as a sensitivity analysis.

In PsychAD cortex, the primary whole-cortex score was not associated with age before adjustment (slope +0.000186 per year, p = 0.3494) or after adjustment for sex, brain bank, postmortem interval, and technical covariates (slope -0.000209 per year, p = 0.4004). One Benjamini-Hochberg correction across 252 primary human cell-type tests left three survivors: MSSM L2/3 intratelencephalic projecting glutamatergic neurons (slope -0.001303, q = 0.0208), MSSM oligodendrocytes (slope +0.001152, q = 0.0346), and PsychAD EN_L3_5_IT_3 neurons (slope -0.001383, q = 0.0208). No oligodendrocyte precursor cell or microglial result survived this primary correction. Glial transport is therefore a candidate finding, not an established mechanism.

In mouse lung fibroblasts, Cdkn1a increased with age (slope +0.1080 per month, p = 0.0029, q = 0.0322). Fibroblast transcriptional noise did not show a reproducible increase after accounting for the UMI sampling model. Across 13 mouse organs with adequate donor coverage, 8 had at least one age-associated cell-type fraction change at FDR 0.05. Three of 11 human cohort datasets had such a change. In integrated human lung, club cells and multiciliated epithelial cells decreased (q = 0.0037 and q = 0.0173), but these associations were not retained in observed-record-only sensitivity analysis. Naked mole-rat skin showed no significant fraction change, but five animals permit only ten age-label assignments and cannot establish stability.

The published Universal Clock 2 coefficients yielded MAE 0.36 years in mouse liver and 5.06 years across 661 human tissue profiles. Predictive R-squared was 0.3169 in mouse liver and 0.6383 in the human multi-tissue evaluation. Mouse scoring used 724 of 816 clock probes. The human evaluation data were part of the consortium resource used to develop the published clocks, so this was a benchmark reproduction rather than independent validation. The five-animal naked mole-rat comparison did not detect a slope and could not establish senescence specificity.

Together, the results do not support a single whole-tissue transcriptional aging program that transports unchanged across species. They support targeted evaluation of cell-type-specific candidates while showing how detection, cohort design, normalization, and multiplicity alter apparent cross-species signals.

---

## 1. Introduction

A central question in geroscience is whether aging involves a conserved set of orthologous genes whose expression changes systematically across tissues and species. Bulk-tissue measurements can mix changes in cell-type abundance with changes in expression within cells. Single-cell and single-nucleus atlases allow these sources of variation to be examined separately, but they also introduce detection, integration, and unequal-sampling problems.

We asked four questions. First, does a transcriptional signature derived from mouse lung and kidney show whole-tissue or cell-type-specific age associations in human cohorts? Second, do cell-type abundances change concordantly in mouse and human lung? Third, does transcriptional noise increase with age after controlling for library sampling? Fourth, how does a published mammalian DNA methylation clock perform when its inverse link, probe availability, and predictive metrics are evaluated explicitly?

The primary analysis used a filtered 139-gene signature throughout. The unfiltered 257-gene discovery signature was analyzed separately as a sensitivity score. We treated cross-sectional associations as observational and did not interpret them as longitudinal or causal effects.

---

## 2. Results

### 2.1 Human transport is cell-type-specific and sensitive to the analysis model

The primary 139-gene score contained 89 mouse up-regulated genes and 50 mouse down-regulated genes. After one-to-one ortholog mapping, 76 up-regulated and 42 down-regulated human genes were available for scoring. The sensitivity score contained 166 up-regulated and 91 down-regulated mouse genes, with 120 up-regulated and 68 down-regulated human orthologs.

In the PsychAD lifespan atlas, the primary whole-cortex score had an unadjusted slope of +0.000186 per year (p = 0.3494). Adjustment for sex, brain bank, postmortem interval, and three technical principal components gave a slope of -0.000209 per year (p = 0.4004). Under the unfiltered sensitivity signature, the corresponding slopes were +0.000318 (p = 0.0118) before adjustment and +0.000027 (p = 0.8626) after adjustment. The primary analysis therefore did not detect a whole-cortex age association, and the nominal unadjusted association depended on the signature definition.

The primary multiplicity family contained 225 public-cohort cell-type tests and 27 fully adjusted PsychAD cell-type tests. Three of 252 tests survived Benjamini-Hochberg correction:

| Cohort and cell type | Age slope per year | Nominal p | Global BH q |
|---|---:|---:|---:|
| MSSM L2/3 intratelencephalic projecting glutamatergic neuron | -0.001303 | 0.000165 | 0.0208 |
| MSSM oligodendrocyte | +0.001152 | 0.000411 | 0.0346 |
| PsychAD EN_L3_5_IT_3 neuron | -0.001383 | 0.000082 | 0.0208 |

The PsychAD neuronal association attenuated in the adult stratum (slope -0.000445 per year, p = 0.3414, q = 0.6512), which indicates that the pooled lifespan result was influenced by developmental donors.

No oligodendrocyte precursor cell result survived the primary correction. SEA-AD OPCs had a nominal positive association (slope +0.003266 per year, p = 0.004583, global q = 0.1155). PsychAD OPCs had a pooled adjusted slope of +0.001028 (p = 0.034284, global q = 0.2700). Their adult-stratum slope was +0.001884 (p = 0.006017), but it did not survive the secondary 27-test adult correction (q = 0.0812).

No microglial result survived the primary correction. PsychAD microglia had a pooled adjusted slope of +0.002069 per year (p = 0.001496, global q = 0.0628). The adult-stratum slope was +0.003115 (p = 0.000617) and survived the secondary adult-only correction (q = 0.0167). MSSM microglia also had a nominal association (slope +0.002018, p = 0.003708, global q = 0.1155). These observations justify targeted follow-up, but they do not establish conserved glial aging.

The primary tier retains sex-dimorphic and immediate-early genes, so it does not test independence from those classes. The unfiltered score remained positive after removal of a generic inflammatory and stress subset (slope +0.001081 per year, p = 0.003876). This decomposition is descriptive and does not identify a causal regulator.

As a positive control in mouse lung adventitial cells, four canonical aging or senescence markers had directionally variable but non-significant effects: Apoe +0.099 log2 fold change (p = 0.630), Il6 +0.778 (p = 0.267), Lgals3 +0.467 (p = 0.344), and Igfbp3 -0.905 (p = 0.312). These null marker results do not validate the differential-expression ranking.

### 2.2 The naked mole-rat comparison remains inconclusive

We evaluated the primary signature in mice, human skin, and naked mole-rat skin using normalized maximum-lifespan units. The mouse whole-body ordinary least-squares slope was +0.1344 (95% CI +0.0293 to +0.2395, p = 0.0147, 23 mice). A model excluding pooled identifiers gave +0.0460 (95% CI -0.0185 to +0.1106, p = 0.1511, 20 mice). A tissue-adjusted cluster-robust model gave +0.2053 (95% CI +0.1159 to +0.2947, p = 0.00064, 152 observations in 19 physical-animal clusters). The late-adult subset did not show a slope (-0.0056, p = 0.9560, 10 mice). Its confidence interval was -0.2338 to +0.2225, so a late-adult plateau is not established.

The naked mole-rat skin slope was +0.0058 (95% CI -0.1981 to +0.2098, p = 0.9332, five animals). The direct mouse-versus-NMR slope difference was +0.1285 (95% CI -0.0314 to +0.2885, p = 0.1153). The exact aged-versus-young permutation test gave p = 1.0000.

The post hoc TOST margin was 0.2023 and exceeded the mouse whole-body ordinary least-squares estimate of 0.1344. The NMR 90% CI was -0.1450 to +0.1567 (TOST p = 0.0275). Equivalence at this margin does not exclude an effect of the observed mouse magnitude. Five animals, group-mean ages, skin-only sampling, and incomplete ortholog coverage prevent a definitive claim of senescence specificity or biological equivalence.

### 2.3 Transcriptional noise depends on the library-sampling model

We compared multinomial resampling and hypergeometric without-replacement subsampling to 1,500 UMI per cell in nine mouse lung cell types. The table reports the committed without-replacement cell counts and both sampling-model results.

| Cell type | Donors | Cells | Multinomial DM slope per month | p | Without-replacement DM slope per month | p | BH q |
|---|---:|---:|---:|---:|---:|---:|---:|
| fibroblast of lung | 10 | 1,475 | +0.00109 | 0.0443 | -0.000288 | 0.4748 | 0.6612 |
| bronchial smooth muscle cell | 13 | 2,320 | +0.000276 | 0.1038 | -0.000219 | 0.0690 | 0.3103 |
| alveolar macrophage | 14 | 1,324 | +0.000560 | 0.4518 | +0.000248 | 0.1265 | 0.3796 |
| classical monocyte | 14 | 7,899 | +0.000078 | 0.8923 | +0.000069 | 0.5142 | 0.6612 |
| intermediate monocyte | 8 | 1,678 | -0.000262 | 0.8439 | +0.000092 | 0.7101 | 0.7989 |
| B cell | 13 | 1,476 | +0.000499 | 0.3359 | -0.000458 | 0.0520 | 0.3103 |
| NK cell | 15 | 1,189 | -0.000830 | 0.0603 | +0.000254 | 0.2653 | 0.5761 |
| CD8-positive, alpha-beta T cell | 10 | 828 | +0.000174 | 0.5338 | -0.000360 | 0.3201 | 0.5761 |
| vein endothelial cell | 8 | 273 | +0.000602 | 0.4787 | +0.000017 | 0.9664 | 0.9664 |

The nominal fibroblast increase under multinomial resampling did not survive FDR correction and reversed sign under without-replacement subsampling. Classical monocytes showed no detected noise slope under either model despite nine differentially expressed genes in the 3-month versus 18-month comparison. Mean-expression change and cell-to-cell variability are therefore distinct measurements.

### 2.4 Lung fibroblast Cdkn1a increases with age

Donor-level pseudobulk analysis of 11 targeted markers showed increased fibroblast Cdkn1a expression with age (slope +0.1080 per month, SE 0.0242, p = 0.0029, BH q = 0.0322, model R-squared = 0.777, ten mice). Sirt6 decreased nominally (slope -0.0136 per month, p = 0.0245), but did not survive correction across the marker panel (q = 0.1350). Ccl2 also did not survive correction (slope +0.0481, p = 0.0475, q = 0.1741). These targeted associations do not establish cellular senescence by themselves.

### 2.5 Cell-type fractions change in some organs, but not universally

Eight of 13 adequately covered mouse organs had at least one age-associated cell-type fraction at FDR 0.05: bladder, kidney, limb muscle, liver, lung, marrow, spleen, and tongue. Fat, heart and aorta, mammary gland, pancreas, and thymus had no detected change. Large intestine, skin, and trachea had insufficient donor coverage for a test.

Mouse kidney showed decreased capillary endothelial and mesangial fractions and an increased plasma-cell fraction. Mouse liver showed increased B-cell and myeloid-leukocyte fractions and a decreased hepatocyte fraction. Mouse lung showed an increased T-cell fraction and a decreased pulmonary interstitial fibroblast fraction. These are changes in cell-type fractions, not evidence of fibrosis or organ destruction.

Three of 11 human cohort datasets had an age-associated fraction change. MSSM cortex had six, including increased oligodendrocytes and decreased oligodendrocyte precursor cells, GABAergic neurons, and SST interneurons. Integrated lung had decreased club-cell and multiciliated epithelial-cell fractions (q = 0.0037 and q = 0.0173). Living-donor kidney had an increased kidney cortex collecting-duct epithelial-cell fraction. The other eight cohort datasets had no detected change at FDR 0.05.

Excluding pooled mouse identifiers retained most mouse findings, but liver NK cells and lung vein endothelial cells were no longer significant. All six MSSM findings remained significant with observed records only. Integrated lung and living-donor kidney findings did not survive that sensitivity analysis, making those results more exploratory.

Naked mole-rat skin showed no significant fraction change. Basal, corneum, and immune compartments had small age-group differences with exact permutation p values above 0.1. Five animals permit only ten assignments of two aged animals, so no compartment can reach FDR 0.05 in this dataset. Basal direction also changed between annotation resolutions 0.4 and 0.6, and only 46 cells were assigned to the immune compartment at the primary resolution. NMR composition is therefore inconclusive rather than demonstrably stable.

### 2.6 Universal Clock 2 benchmark with predictive metrics

We applied the published Universal Clock 2 coefficients and double-exponential inverse link to 339 mouse liver profiles and 661 human tissue profiles. Predictive R-squared was calculated as one minus residual sum of squares divided by total sum of squares from per-sample chronological-age predictions. It was not calculated by squaring the Pearson correlation.

| Species and tissue | Samples | Probes used | Pearson r | Predictive R-squared | MAE in years |
|---|---:|---:|---:|---:|---:|
| Mouse liver | 339 | 724 of 816 | 0.923 | 0.3169 | 0.36 |
| Human all tissues | 661 | 816 of 816 | 0.801 | 0.6383 | 5.06 |
| Human lung | 113 | 816 of 816 | 0.839 | 0.6655 | 4.89 |
| Human kidney | 97 | 816 of 816 | 0.876 | 0.7560 | 4.10 |
| Human liver | 97 | 816 of 816 | 0.817 | 0.6474 | 4.96 |
| Human heart | 97 | 816 of 816 | 0.808 | 0.6348 | 5.18 |
| Human muscle | 57 | 816 of 816 | 0.728 | 0.5106 | 5.73 |
| Human blood | 32 | 816 of 816 | 0.799 | 0.4777 | 6.14 |
| Human adipose | 57 | 816 of 816 | 0.823 | 0.6733 | 4.50 |
| Human spleen | 62 | 816 of 816 | 0.838 | 0.6816 | 5.27 |
| Human lymph node | 27 | 816 of 816 | 0.539 | 0.1560 | 7.82 |
| Human bone marrow | 20 | 816 of 816 | 0.866 | 0.7348 | 4.39 |
| Human pituitary | 2 | 816 of 816 | 1.000 | 0.7685 | 3.88 |

The mouse implementation used only the 724 clock probes present in the array file, without imputation. The original publication identifies GSE184211 as Mammalian Methylation Consortium data and describes development of the clocks from consortium profiles. The human evaluation therefore overlaps the development resource and is not independent validation. Sample-level training membership for GSE199979 was not established by the inspected sources, so the mouse evaluation is also not presented as independent validation.

---

## 3. Discussion

The primary analysis did not support a whole-cortex transcriptional aging association that transported unchanged from mouse to human. Three cell-type tests survived the shared human multiplicity correction, but the PsychAD neuronal result attenuated in adults. Microglia showed an adult-stratum association in PsychAD, yet neither microglia nor OPCs survived the primary 252-test correction. Glial transport remains a candidate for directed validation rather than an established conserved mechanism.

The human lung analysis illustrates why atlas composition cannot be interpreted from thresholded score tables alone. Restoring counts below the expression-score threshold changed the evaluable cohort and cell-type set. Detection presence was itself associated with age, and study-center adjustment attenuated the leading club-cell result. These are constraints on inference, not evidence that the compositional changes are absent.

The naked mole-rat comparison did not detect a score slope. Its small sample, group-mean ages, tissue mismatch, and post hoc equivalence margin prevent a claim of negligible senescence or species-specific aging. The transcriptional-noise result also depended on the UMI sampling model. In contrast, fibroblast Cdkn1a showed an FDR-significant mean-expression increase, although this single marker is not sufficient to define a senescent cell state.

The methylation analysis was a reproduction of a published clock on consortium data. The inverse link yielded low chronological-age errors, but predictive R-squared was substantially below squared Pearson correlation in mouse liver. This distinction matters because correlation measures ordering, whereas predictive R-squared measures error relative to the observed-age variance.

### Limitations

1. All molecular cohorts were cross-sectional. Age associations may reflect cohort composition, technical variation, development, disease, or postmortem effects.
2. Mouse differential expression was concentrated in lung and kidney units that met the analysis floor. Partial sex-by-age confounding and pooled animal identifiers limit causal and lifespan-wide interpretation.
3. Public human scoring used per-cell library-size normalization but lacked the full technical covariates available in PsychAD. The surviving public-cohort associations require independent confirmation.
4. PsychAD and the public MSSM and HBCC cohorts draw from overlapping source populations. Shared identifiers were unavailable, so exact donor-level independence could not be established.
5. Human lung detection and study origin were associated with some compositional estimates. Integrated-atlas composition should not be interpreted as direct tissue cell counting.
6. The NMR comparison included five animals, two group-mean age categories, skin only, and incomplete signature coverage.
7. The methylation evaluations were not demonstrated to be independent of clock development data. The mouse clock implementation omitted 92 unavailable probes.

---

## 4. Methods

### 4.1 Data and inferential units

Mouse analyses used Tabula Muris Senis droplet data from 23 mice, with 16 lung-profiled mice for lung composition. Public human analyses used 11 CELLxGENE cohorts and 621 scored control donors. The PsychAD DLPFC lifespan analysis used 284 control donors aged 0.2 to 97 years across 27 cell types. Methylation benchmarking used GSE199979 and GSE184211 on GPL28271. The NMR comparison used five animals from GSE152886. Donors or animals were the inferential units unless a tissue-clustered model is explicitly named.

### 4.2 Signature definitions and scoring

The discovery signature was derived from mouse lung and kidney differential expression. The primary filtered signature contained 139 genes: 89 up-regulated and 50 down-regulated. Its filters removed predicted genes, pseudogenes, noncoding and non-protein-coding entries, ribosomal and mitochondrial genes, and specified ambient-RNA markers. It retained sex-dimorphic and immediate-early genes. The sensitivity signature contained 257 genes: 166 up-regulated and 91 down-regulated. Mouse-human mapping used one-to-one orthologs.

For donor or donor-cell-type pseudobulk, expression was transformed as log(1 + counts per million). The score was the mean transformed expression of up-regulated genes minus the mean transformed expression of down-regulated genes. Public cohort cells were library-size normalized before donor-cell-type averaging. PsychAD matrices were normalized within donor-cell-type pseudobulk. Missing genes were omitted rather than imputed. Every primary result used the filtered signature; the unfiltered score was always labeled as sensitivity.

### 4.3 Human transport models and multiplicity

Public-cohort cell-type regressions included age and, when variable, sex. PsychAD models included age, sex, brain bank, postmortem interval, and technical principal components 1 to 3. The primary PsychAD model used all donors. Adult models restricted donors to at least 20 years and were secondary sensitivity analyses.

One Benjamini-Hochberg family contained all 225 public-cohort cell-type tests and all 27 primary PsychAD cell-type tests. A result was called a primary survivor only at global q below 0.05. The 27 adult PsychAD tests formed a separate secondary sensitivity family. Nominal results that did not survive the relevant correction were not treated as validation.

### 4.4 Donor-overlap audit

PsychAD contained 112 MSSM and 172 HBCC donors. Public and PsychAD identifiers could not be matched directly. Age, sex, and brain-bank comparisons gave at least one candidate match for 160 of 214 public MSSM donors and 133 of 134 public HBCC donors. Fifty-two MSSM and 35 HBCC public donors had a unique candidate on those fields, but these matches were not unique identifiers and cannot confirm individual identity. PsychAD was therefore treated as a second analysis of overlapping source populations, not independent replication.

### 4.5 Composition and missing-data treatment

Each mouse organ and each human cohort was analyzed separately. Human donors younger than 20 years were excluded. Complete donor-by-cell-type counts were obtained before the 20-cell expression-score filter. Counts of 1 to 19 cells remained unchanged. A donor-cell-type combination absent from the complete atlas enumeration represented zero observed cells, not a row discarded by a score threshold. Cell types observed in at least six donors were tested.

Counts were transformed with a 0.5 pseudocount and centered log-ratio. Regression included age and sex, with sex dropped only when the design was rank-deficient. BH correction was applied within each organ or cohort. Sensitivity analyses excluded pooled mouse identifiers or restricted human tests to donors with at least one observed cell of the target type.

The NMR matrix lacked cell-type labels, so cells were re-annotated using the paper's processing outline and marker families. Sample 2V was excluded. Immune markers were Cd207, Ctss, and Mfge8; basal markers were Krt5, Trp63, and Krt15; corneum markers were Krt10 and Krt1. Remaining keratinocyte clusters were labeled intermediate. Exact permutations compared mean aged and young fractions across all ten assignments of two aged animals among five. Annotation was repeated at resolution 0.6 as a sensitivity check. The minimum possible two-sided NMR p value was 0.1.

### 4.6 Transcriptional noise and targeted markers

Singlet cells with at least 1,500 UMI were downsampled to 1,500 UMI using multinomial resampling or hypergeometric sampling without replacement. Distance-to-median and donor-centroid dispersion were regressed on age and sex. Eleven senescence, DNA-repair, and chromatin markers were evaluated in lung fibroblast donor pseudobulk with BH correction across the marker panel.

### 4.7 Positive controls

Apoe, Il6, Lgals3, and Igfbp3 were evaluated in the top-ranked mouse lung adventitial-cell unit using the same young-versus-old effect scale as the differential-expression analysis. Their outcomes were reported regardless of significance and were not used as evidence for the ranking.

### 4.8 Epigenetic clock evaluation

Published Universal Clock 2 coefficients were applied without refitting. Chronological age was scaled using species maximum lifespan and gestation time, and the double-exponential inverse link was applied to the linear predictor. Human metadata were aligned to methylation columns by array position before scoring. Mouse scoring used all clock probes present in the array file. Per-sample observed and predicted ages were saved. Predictive R-squared was calculated as 1 - SSE/SST on chronological-age predictions. Pearson correlation and MAE were reported separately.

### 4.9 Literature comparison

Tyshkovskiy and colleagues integrated bulk transcriptomes across mammalian species and tissues to derive conserved aging and mortality modules. Their work and ours both identify age-associated molecular structure across species, and CDKN1A is relevant to both analyses. The estimands differ: their conserved modules and mortality clocks are not equivalent to transport of our mouse lung-and-kidney single-cell signature. Our results add a direct test of how normalization, cell type, development, and multiplicity affect transport of one such signature. They do not refute their broader modules.

Karakaslar and colleagues identified increased Jun and Fos activity in aging mouse immune cells and chromatin-level evidence for conserved AP-1 regulation. Our mouse differential-expression analysis and their study both use Tabula Muris Senis, so shared atlas observations are not independent replication. Their functional immune experiments go beyond our observational analysis. The primary signature retains immediate-early genes, so our analysis does not test transport independently of that AP-1-related class. The PsychAD microglial candidate is not presented as validation of their immune AP-1 mechanism.

---

## References

1. Aitchison, J. (1986). *The Statistical Analysis of Compositional Data*. Chapman and Hall, London.
2. Gabitto, M. I., et al. (2023). Integrated single-cytometry and spatial transcriptomics reveals early cellular and molecular changes in Alzheimer's disease. *bioRxiv*, doi:10.1101/2023.05.08.539485.
3. Karakaslar, E. O., et al. (2023). Transcriptional activation of Jun and Fos members of the AP-1 complex is a conserved signature of immune aging that contributes to inflammaging. *Aging Cell*, 22, e13792. doi:10.1111/acel.13792.
4. Lu, A. T., et al. (2023). Universal DNA methylation age across mammalian tissues. *Nature Aging*, 3(9), 1144 to 1166. doi:10.1038/s43587-023-00462-6.
5. Martinez-Jimenez, C. P., et al. (2017). Aging increases cell-to-cell transcriptional variability upon immune stimulation. *Science*, 355(6332), 1433 to 1436.
6. Ruby, J. G., Smith, M., and Buffenstein, R. (2018). Naked mole-rat mortality rates defy Gompertzian laws by not increasing with age. *eLife*, 7, e31157.
7. Saul, D., et al. (2022). A new gene set identifies senescent cells and predicts senescence-associated pathologies across tissues. *Nature Aging*, 2(3), 218 to 227.
8. Sikkema, L., et al. (2023). An integrated cell atlas of the human lung in health and disease. *Nature Medicine*, 29(6), 1563 to 1577.
9. Tabula Muris Consortium. (2020). A single-cell transcriptomic atlas characterizes ageing tissues in the mouse. *Nature*, 583(7817), 590 to 595.
10. Tyshkovskiy, A., et al. (2026). Universal transcriptomic hallmarks of mammalian ageing and mortality. *Nature*, 654, 173 to 188. doi:10.1038/s41586-026-10542-3.
11. Zhang, X., et al. (2021). Characterization of the conserved mammalian aging signature through single-cell cross-species meta-analysis. *Cell Reports*, 36(12), 109763.

---

## Data and Code Availability

Analysis code and aggregate result tables are available at `https://github.com/Abolfazl-Alipour/agescape`. The interactive explorer is available at `https://rezaalipour.com/agescape/`.

The results published here are in whole or in part based on data obtained from the AD Knowledge Portal (NIH-NIA, U24AG061340). The PsychAD DLPFC Lifespan Atlas is available through the AD Knowledge Portal at `https://adknowledgeportal.synapse.org/` under Synapse `syn52396927` and pseudobulk resource `syn62064518`, subject to its controlled-access requirements. Only aggregate regression results are shared in this repository. Individual-level controlled data and derived donor-level values are not redistributed.

The remaining public datasets are identified by their atlas names and accession numbers in Methods. Public per-sample clock predictions and aggregate statistical outputs are included to permit numerical verification.
