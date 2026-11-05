# Reviewer-response Phase 1 methods

## Scope

Phase 1 reanalyzes existing datasets and committed outputs. No new single-cell dataset was collected. The manuscript was not edited. All new outputs use `reviewer-phase1-` in their filenames. All analyses are reviewer-response analyses and are exploratory unless they merely reproduce an existing primary calculation.

## Mouse drift decomposition

The full mouse expression-vector estimate was recomputed with the existing donor-level model. Each cell type was removed in turn across all tissues. The output records the revised global estimate, standard error, p value, change from the full estimate, and number of tissue-cell-type units removed. The proportion of positive unit slopes was accompanied by an exact two-sided sign test.

## Gene-level concordance

Each cell type within each tissue or cohort contributed a vector of donor-pseudobulk gene-level age effects. Counts were aggregated by donor and cell type, normalized to counts per 10,000, and transformed with log1p. AFCA uses the released log-normalized matrix, with `expm1` applied before donor aggregation. Models included age, sex when variable, mean log UMI, and sequencing platform for Gao worm. At least six inferential units were required, except the five-animal NMR analysis.

Pairwise Pearson and Spearman correlations were calculated between cell-type gene-effect vectors within a tissue or cohort. The primary summary is the mean pairwise Pearson correlation, with median and Spearman summaries retained. Age labels were permuted independently within cell type while preserving sex and platform strata. Two hundred permutations were used. NMR used all ten assignments of two aged animals among five. Species summaries weight stratum summaries by their number of cell-type pairs. Permutation intervals describe the null, not confidence intervals for a biological effect.

A random-effects meta-analysis across species used existing global Hellinger effects and their standard errors. The balanced worm day-6 reference was primary and day 1 was retained as sensitivity. NMR used a continuous-age, animal-clustered model with cell-type intercepts because the exact two-group output has no standard error. Cross-species effect scales remain descriptive because tissues, assays, and age coverage differ.

## Nulls, power, and equivalence

Existing effect sizes and standard errors were used where available. Off-program and entropy standard errors were recovered from their t-based confidence intervals. NMR used the same clustered continuous-age model as the species meta-analysis. Zebrafish used the fine-grained depth-adjusted global model.

Minimum detectable effects at 80% power use the two-sided 0.05 t critical value plus the 80th-percentile t value, multiplied by the reported standard error. Equivalence uses two one-sided t tests with standardized margins of 0.1 and 0.25, both fixed before calculation and reported without selecting a preferred outcome. Exact NMR group-permutation inference remains in the original output; the clustered model is used only for interval and power calculations.

## Interneuron replication

MSSM and HBCC gene-level slopes were compared with Pearson and Spearman correlations, direction concordance, and exact binomial tests. Both-cohort FDR overlap was compared with a hypergeometric null using all 59,523 shared tested genes. VIP, pvalb, and sst subtype gene-set overlaps were tested similarly. Dot-plot data contain the top ten MSSM FDR genes per subtype plus all genes shared across the three subtypes.

GO Biological Process enrichment uses static annotation libraries from Enrichr and FlyEnrichr. The human library is GO Biological Process 2023; the fly library is GO Biological Process 2018 because the 2023 library is unavailable for fly. Hypergeometric tests use all measured genes with GO annotations as background, restrict terms to 5 to 500 background genes, and apply BH correction across all tested list-term pairs. Enrichment is annotation-dependent and does not establish mechanism.

## Fly hit coherence

The 38 AFCA cell types with depth-adjusted FDR-significant Hellinger effects were selected from the committed vector table. Donor-pseudobulk gene effects were computed independently for those units. Pairwise overlaps of their FDR-significant gene sets were compared with hypergeometric expectations. GO enrichment used the same rules as the interneuron analysis with the AFCA gene universe. A shared pathway annotation is not proof of a shared causal program.

## Multiplicity

The family inventory records hypothesis, family, correction, threshold, and inferential status. The core composition recheck applies one BH correction across all primary mouse, human, and NMR cell-fraction tests, then summarizes significant results by organ or cohort. This directly evaluates the original 8-of-13 and 3-of-11 summaries under a wider family.

The signature-transport recheck applies BH to the same 252 primary hypotheses. Overlapping donors affect independence claims, but they do not define which distinct cell-type hypotheses should be deleted. No tests were removed based on outcome. The recheck is labeled overlap-aware in interpretation, not a deduplicated joint-donor model.

## Composition sensitivity

CLR composition models were rerun with log total recovered cells, sex when variable, and platform for Gao worm. Pseudocounts 0.1, 0.5, and 1.0 were evaluated. FDR was calculated within organ or cohort and, for the pseudocount-0.5 arm, across all tests.

The donor-level count-model check used propeller's logit transformation `(count + 0.5) / rowSum(count + 0.5)` followed by limma empirical-Bayes regression. Models included age, log total recovered cells, sex when variable, and platform for Gao worm. Implementation used Bioconductor `speckle` 1.10.0 and limma in a temporary R environment. Empty or inadequately covered atlas strata were excluded before fitting.

## P1 signature and matched null

P1_Clean contains 139 mouse genes, 89 up and 50 down. It was derived from the mouse analysis before the reported human transport evaluations. P0_Raw257 is the 257-gene sensitivity tier. Mouse-human mapping uses the committed one-to-one Ensembl-derived table, with 76 human up genes and 42 human down genes. Missing genes are omitted, never imputed. PsychAD detected 113 of the 118 mapped genes in each reported cell type.

Two hundred random signatures per PsychAD cell type were generated with a fixed seed. Each matched P1's detected up and down counts and expression-decile distribution. P1 genes were excluded from the random candidate pool. Fully adjusted age models matched the primary PsychAD model. Outputs include empirical two-sided p values, BH q values, gene coverage, mean log library size, and score-depth correlations. Only aggregate results were exported from the controlled resource.

## Worm references

The day-1 reference uses the earliest young adult cells but is entangled with sequencing platform. Day 6 includes both platforms and is the balanced sensitivity reference. Both are reported. The day-6 result is used for the primary reviewer-response meta-analysis because it reduces reference-platform confounding; this is a design justification, not selection by effect size.

## Reproduction

```bash
python scripts/run_reviewer_phase1_core.py --equivalence-margin 0.1 0.25
python scripts/run_reviewer_phase1_gene_concordance.py --skip-human --permutations 200
python scripts/run_reviewer_phase1_gene_concordance.py --species Human --permutations 200 --output-suffix=-human
python scripts/run_reviewer_phase1_gene_concordance.py --species 'C. elegans' --permutations 200 --output-suffix=-worm
python scripts/run_reviewer_phase1_gene_concordance.py --species Zebrafish --permutations 200 --output-suffix=-zebrafish
python scripts/run_reviewer_phase1_gene_concordance.py --species 'Naked mole-rat' --output-suffix=-nmr
python scripts/combine_reviewer_phase1_concordance.py
python scripts/run_reviewer_phase1_meta.py
python scripts/run_reviewer_phase1_enrichment.py
python scripts/run_reviewer_phase1_composition.py
python scripts/prepare_reviewer_phase1_propeller.py
Rscript scripts/run_reviewer_phase1_propeller.R data/interim/reviewer-phase1-propeller results/reviewer-phase1-composition-propeller.csv
python scripts/run_reviewer_phase1_signature_null.py
python scripts/summarize_reviewer_phase1_composition.py
python scripts/summarize_reviewer_phase1.py
```
