# Draft Methods for author revision

This is a package input for the author. It is not an edit to the repository manuscript. All analyses were observational and cross-sectional unless noted otherwise.

## Datasets and inferential units

| Analysis arm | Dataset and accession | Inferential units and analysis scope |
|---|---|---|
| Mouse atlas | Tabula Muris Senis droplet data, distributed through CZ CELLxGENE | 23 mouse identifiers; 238,186 cleaned cells; tissue-cell-type analyses used eligible units with adequate donor coverage |
| Public human brain, kidney, and lung | CZ CELLxGENE Census release 2025-11-08, 11 dataset IDs listed below | 621 control donors for signature scoring; 479 cohort-specific donors with adequate young-reference coverage for expression vectors; composition cohorts used all eligible adult donors |
| Human lifespan cortex | PsychAD DLPFC Lifespan Atlas, Synapse syn52396927 and pseudobulk resource syn62064518 | 284 control donors, 27 cell types; donor-cell-type pseudobulk |
| Naked mole-rat skin | Savina and colleagues, GEO GSE152886 | Five animals after exclusion of sample 2V; young and aged groups contained three and two animals |
| Fruit fly | Aging Fly Cell Atlas, GEO GSE218661 | 36 age-by-sex pools measured in head and body; primary ages 30, 50, and 70 days |
| Gao worm atlas | GEO GSE229022 | 12 wild-type OP50 pools at days 1, 6, 12, and 14 |
| Roux worm atlas | GEO GSE208154 | Technical replicate series at days 1, 3, 5, 8, 11, and 15; descriptive only |
| Zebrafish testis | Sposato and colleagues, GEO GSE275361 | Ten males for broad compartments; six fish aged 5 to 22 months for fine-grained analyses |
| Mouse methylation | GEO GSE199979, GPL28271 | 339 liver array profiles; benchmark reproduction |
| Human methylation | GEO GSE184211, GPL28271 | 661 multi-tissue array profiles; benchmark reproduction |
| Longevity scaling | AnAge Build 14 | Species maximum-lifespan records used only for age scaling and context |

Public human dataset IDs were:

| Cohort | CELLxGENE dataset ID | Vector donors | Composition donors |
|---|---|---:|---:|
| MSSM DLPFC | 37a17b78-4864-4a42-b67b-31c00962795a | 212 | 213 |
| HBCC DLPFC | 5c97eeeb-7e52-44b3-b010-b832b1f5424c | 88 | 88 |
| RADC DLPFC | 5e57cd50-8e42-42d6-940d-5c1660d06864 | Not evaluable | 18 |
| SEA-AD DLPFC | 6f7fd0f1-a2ed-4ff1-80d3-33dde731cbc3 | Not evaluable | 21 |
| Kidney Normal Adult | 0b4a15a7-4e9e-4555-9733-2423e5c66469 | Not evaluable | 18 |
| Kidney Living Donor | 0bae7ebf-eb54-46a6-be9a-3461cecefa4c | 19 | 19 |
| Kidney SN | a12ccb9b-4fbe-457d-8590-ac78053259ef | Not evaluable | 9 |
| Kidney SCR | dea717d4-7bc0-4e46-950f-fd7e1cc8df7d | 18 | 18 |
| Lung CellRef | eb499fd8-7000-419f-8854-926b9b61b11f | 17 | 17 |
| Lung Integrated | 9f222629-9e39-47d0-b83f-e08d610c7479 | 113 | 119 |
| Immune Aging BAL | 1b350d0a-4535-4879-beb6-1142f3f94947 | 12 | 12 |

Only primary, normal-disease human cells were included. Public human analyses used donors aged at least 20 years. Donors, animals, fish, or pools were the inferential units, not individual cells. PsychAD and public MSSM/HBCC cohorts contain exact identifier-confirmed overlapping donors and were not treated as independent replications. The MSSM-to-HBCC interneuron comparison used different source banks.

## Mouse processing and differential expression

Mouse cells were filtered using the prespecified atlas quality-control rules: 500 to 5,000 detected genes, at least 1,000 UMIs, and doublet removal. Mitochondrial filtering was not applicable because the released droplet object lacked mitochondrial transcripts. Counts were normalized to 10,000 per cell and transformed with log1p. Highly variable genes, principal components, nearest-neighbor graphs, UMAP, and Leiden clustering were used for atlas reproduction. Dataset-provided cell-ontology labels were used for downstream cell-type analyses.

For differential expression, donor-cell-type pseudobulk counts were aggregated and log-transformed after counts-per-million scaling. Age effects were estimated within each cell type with sex included when identifiable. P-values used t-based regression inference and Benjamini-Hochberg correction within the stated family. Effect sizes were retained alongside p-values. Sensitivity analyses varied cell-count floors, pseudobulk specifications, and donor exclusions without selecting a preferred result.

## Signature derivation and transport

The primary P1_Clean signature contained 139 mouse genes, 89 up-regulated and 50 down-regulated. P0_Raw257 contained 257 genes and was retained as a sensitivity score. P1 used uniform-BH differential-expression candidates from adequately covered mouse units and removed predicted or pseudogene labels, non-protein-coding biotypes, ribosomal and mitochondrial genes, and specified ambient-RNA cases. Sex-dimorphic and immediate-early filters belonged to P2 and P3, not P1. Mouse signature derivation preceded the reported human transport evaluations; no human result was used to choose P1 membership.

One-to-one mouse-human orthologs were taken from the committed Ensembl-derived mapping. There were 118 mapped P1 genes, 76 up and 42 down. Missing genes were omitted rather than imputed. Public human cells were library-size normalized before donor-cell-type aggregation. PsychAD pseudobulk counts were CPM normalized before log1p transformation. The score was mean expression of up-regulated genes minus mean expression of down-regulated genes.

Public-cohort regressions included age and sex when variable. Primary PsychAD models included age, sex, brain bank, postmortem interval, and three technical principal components. Adult models restricted donors to at least 20 years and were secondary sensitivity analyses. One BH family covered 225 public-cohort cell-type tests and 27 primary PsychAD tests, for 252 tests. The 27 adult PsychAD tests formed a separate secondary family.

A matched-null analysis in PsychAD generated 200 random signatures per cell type. Each matched the detected P1 up/down counts and expression-decile distribution, excluded P1 genes from the candidate pool, and used the same up-minus-down score and fully adjusted regression. Empirical p-values and BH q-values were reported with gene coverage and score-depth correlations. Only aggregate results were exported from the controlled resource.

## Cell-type composition

For donor `i` and cell type `j`, counts were transformed with pseudocount `a` and centered log-ratio:

`CLR_ij = log((count_ij + a) / sum_j(count_ij + a)) - mean_j(log((count_ij + a) / sum_j(count_ij + a)))`.

The primary pseudocount was 0.5. Primary regressions included age and sex when identifiable. Existing sensitivity analyses excluded pooled mouse identifiers or used observed human records only. Reviewer-response models added log total recovered cells and evaluated pseudocounts 0.1, 0.5, and 1.0. Sequencing platform was included for Gao worm. NMR used exact age-label permutations because five animals provided only ten assignments of two aged animals.

Propeller provided a donor-level count-model check. Proportions were transformed as `logit((count + 0.5) / rowSum(count + 0.5))`, followed by limma empirical-Bayes regression with age, log total cells, sex when variable, and platform for Gao worm. Full effect sizes, confidence intervals, and FDR values were retained. Recovered composition was not interpreted as absolute tissue cell number because dissociation, nuclear isolation, quality control, survival, and annotation can vary by age.

Cross-species homolog comparisons used prespecified mouse-human cell-type mappings. Sign concordance was tested with an exact binomial test and slope concordance with Spearman correlation. The seven-pair Spearman result was a genuine exact zero with nonzero variance in both arms.

## Whole-expression-vector displacement

No UMI downsampling was used for the direct vector analysis. Mitochondrial, ribosomal, and sex-linked marker genes were excluded using species-appropriate symbols. Each cell's raw counts were converted to within-cell proportions `p`. A young reference probability vector `q` was constructed from pooled counts of the same annotated cell type. Hellinger distance was:

`H(p,q) = sqrt(1 - sum_g sqrt(p_g q_g))`.

Young donors or animals were excluded from their own reference. Each evaluated donor-cell-type unit required at least 20 cells and at least 50 remaining young-reference cells. Donor means were the inferential units. The literal raw-vector sensitivity endpoint used Euclidean distance to the young mean raw-count vector and was explicitly treated as depth-sensitive.

Mouse young references used animals aged 1 or 3 months. Human references used donors aged 20 to 40 years. NMR used animals aged 1.1 years as young and 11.3 years as aged. AFCA used day 30 as the young-adult reference and excluded day-5 cells from the primary analysis because age and study batch were inseparable. Gao worm used day 1 and day 6 in separate analyses. Day 6 included both sequencing platforms and was the balanced reviewer-response reference. Zebrafish fine analyses used six fish aged 5 to 22 months.

Continuous-age models included sex and, in sensitivity or primary reviewer-response arms as specified, mean log UMI. Platform was included for Gao worm, and incomplete spermatogenesis was included in zebrafish sensitivity analyses. Species summaries used unit fixed intercepts and donor-clustered standard errors. Effects were scaled to 10% of species maximum lifespan using AnAge values of 4 years for mouse, 122.5 years for human, and 31 years for NMR. These scales did not make different assays and tissues biologically equivalent.

The mouse leave-one-cell-type-out analysis removed each cell type across all tissues and refit the global model. Sign tests and unit-level slopes assessed whether the global association depended on one cell type.

## Gene-effect concordance across cell types

Donor-pseudobulk counts were aggregated per cell type, normalized to counts per 10,000, and log1p transformed. AFCA used `expm1` on the released log-normalized matrix before donor aggregation. Each cell type yielded a vector of gene-level age effects. Models included age, sex when variable, mean log UMI, and platform for Gao worm. At least six inferential units were required except for the five-animal NMR analysis.

Pairwise Pearson and Spearman correlations were computed between cell-type effect vectors within a tissue or cohort. The primary summary was mean pairwise Pearson correlation. Age labels were permuted within cell type while preserving sex and platform strata. Two hundred permutations were used. NMR used all ten assignments of two aged animals among five. Species summaries weighted strata by their number of cell-type pairs. This analysis tested shared gene-level age effects, not universal whole-transcriptome deterioration.

A random-effects meta-analysis combined the six global Hellinger estimates using inverse-variance weighting and a DerSimonian-Laird residual-heterogeneity estimate. Day 6 was used for the primary worm input and day 1 for sensitivity. NMR used a continuous-age animal-clustered model for its standard error, supplementing the exact group test.

## Entropy, off-program expression, and depth sensitivity

Within-cell Shannon entropy was `-sum_g p_g log(p_g)` after the stated gene exclusions. Off-program genes were uncommon in young target cells but characteristic of another lung cell type: detection in at most 5% of young target cells and at least 20% of another cell type, expression at least 20 counts per 10,000, and at least fivefold enrichment over the target. Young references excluded the evaluated young animal. Animal means were regressed on age and sex with BH correction across cell types.

The original analysis downsampled all cells without replacement to 1,500 UMIs. Reviewer-response analyses evaluated 500, 1,500, and 3,000 UMIs and full depth. Full-depth analysis used no downsampling and retained the original 1,500-UMI inclusion threshold. Young-reference marker definitions were fixed across depths. Off-program fractions used actual total UMIs as the denominator. The library-size arm added each animal's mean log actual pre-downsampling depth. Global models used cell-type intercepts and animal-clustered t-based inference. A positive full-depth off-program association was not detected. This did not establish equivalence or prove stable epigenetic regulation.

## Human nonlinear age models

Public human vector and composition inputs already contained only adults. Nonlinear age effects used a cubic B-spline with three degrees of freedom, fixed before examining results. F tests compared nonlinear and linear models on the same donor-cell-type units. Vector models included sex and depth, and composition models used CLR outcomes with sex when variable. BH correction was applied within cohort separately for linear and nonlinear comparisons. A nonlinear F test was not interpreted as evidence for a particular direction of aging.

## Annotation stability

No per-cell annotation-confidence probability was available. Excluded and unassigned cells were tabulated by age and tissue or cohort. Mouse exclusions compared raw and cleaned droplet objects. Human exclusions used primary normal-disease cells and the 1,500-UMI threshold. Spearman tests assessed age associations of excluded and unassigned fractions.

A fixed canonical marker panel was evaluated by donor-cell-type log-CP10K expression with age, sex, and depth covariates. At least six donors and 20 cells per donor-cell-type unit were required. BH correction was applied across targeted marker tests. Marker change was not equated with identity loss, and targeted markers could not validate all annotations.

## Bayesian hierarchical models

Bayesian disorder and vector outcomes were standardized within cell type. Models used global effects, non-centered cell-type intercepts and age slopes, animal or donor intercepts, cell-type residual scales, and Student-t residuals. The primary global age prior was Normal(0, 0.5), with prior SDs 0.25 and 1.0 as sensitivity analyses. Other priors were Normal(0, 0.5) for global intercepts and covariates, Gamma(2, 2) for cell-type and animal scales, and Gamma(2, 0.1) for Student-t degrees of freedom. Posterior means, 94% highest-density intervals, and probabilities of a positive effect were reported. Prior-predictive, convergence, posterior-predictive, calibration, and influential-observation diagnostics were retained. Bayesian estimates were exploratory when cohorts were small or diagnostics warned.

The proliferation-class Bayesian meta-analysis used signed unit effects and their standard errors, tissue-by-cohort intercepts, a unit-level heterogeneity scale, and a Normal(0, 0.5) prior for the post-mitotic minus proliferative-capable contrast. Prior SDs 0.25 and 1.0 were sensitivity analyses. Absolute effects did not receive a normal measurement-error model.

## Methylation benchmark

Published Universal Clock 2 coefficients were applied without refitting. Chronological age was scaled using species maximum lifespan and gestation time, and the published inverse-link correction was applied. Human metadata were aligned to methylation columns by array position. Mouse scoring used 724 of 816 clock probes present in the array file. Per-sample predictions were saved. Predictive R-squared was `1 - SSE/SST` on chronological-age predictions, distinct from squared Pearson correlation. MAE and correlations were reported separately. GSE184211 overlapped the clock-development resource, and GSE199979 sample-level training membership remained unresolved after checking the publisher's supplementary workbook. Neither evaluation was presented as independent validation.

## Literature comparison

Comparisons with Tyshkovskiy and Gladyshev distinguished bulk conserved aging and mortality modules from transport of this single-cell signature. Comparisons with Karakaslar distinguished shared TMS observations from independent replication and functional AP-1 experiments from observational associations. GO enrichment used static GO Biological Process libraries, with hypergeometric tests against measured annotated genes, term sizes 5 to 500, and BH correction. Pathway enrichment was treated as annotation-dependent and not causal evidence.

## Reproducibility and reporting

Public raw-data versions, download dates, and checksums were recorded in the data manifests. Aggregate result tables were committed for every reported number. Controlled PsychAD donor-level values and identifiers were not exported. Exploratory, secondary, and primary families were labeled explicitly. Non-significance was not treated as equivalence unless a stated equivalence test supported it.
