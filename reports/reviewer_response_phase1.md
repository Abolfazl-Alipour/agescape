# Reviewer response Phase 1: reanalysis of existing outputs

## Scope and provenance

Phase 1 used existing datasets and committed outputs. The starting commit was `32227ce2075f4bbac9f0034c29f2e01a68ef4f6b`. No manuscript file was edited. No new single-cell dataset was collected. Exact analysis commands are listed in `notes/reviewer-phase1-methods.md`. All new result files begin with `reviewer-phase1-`.

## 1. Gene-effect concordance and species meta-analysis

The universality analysis compares gene-level age-effect vectors between cell types within a tissue or cohort. It does not test whether every gene changes, whether whole-vector distance increases, or whether the same genes respond across species. Age labels were permuted within cell type, preserving sex and platform strata. Two hundred permutations were used, except NMR, which used all ten age assignments.

Results across all six species:

| Species | Mean pairwise Pearson concordance | Null mean | Null 95% interval | Permutation p greater |
|---|---:|---:|---:|---:|
| Mouse | +0.1804 | +0.0011 | -0.0113 to +0.0180 | 0.004975 |
| Fruit fly | +0.1936 | +0.0001 | -0.0017 to +0.0029 | 0.004975 |
| Human | +0.1981 | -0.0003 | -0.0067 to +0.0092 | 0.004975 |
| C. elegans | +0.3537 | -0.0025 | -0.0219 to +0.0343 | 0.004975 |
| Zebrafish | +0.3393 | -0.0005 | -0.0321 to +0.0532 | 0.004975 |
| Naked mole-rat | +0.6146 | +0.6423 | +0.5346 to +0.7402 | 0.7000 |

Source rows for the mean Pearson statistic: mouse `results/reviewer-phase1-universality-summary.csv:15`, fly `:7`, human `:11`, worm `:3`, zebrafish `:23`, and NMR `:19`. The NMR result is unresolved because its null concordance is already high and only two compartments with five animals were evaluable. Positive concordance in the other five species indicates shared direction of gene-level effects across cell types. It does not establish a universal increase in transcriptional disorder.

The random-effects meta-analysis of global Hellinger effects used the balanced worm day-6 reference. The pooled effect was +0.0721 SD per species-specific age scale, with 95% CI -0.0232 to +0.1674 and p = 0.1383. Heterogeneity remained, with Q = 15.20 on 5 degrees of freedom, p = 0.00956, and tau-squared = 0.00545. With the worm day-1 reference, the pooled effect was +0.1454, 95% CI -0.0182 to +0.3090, p = 0.0815. Sources: `results/reviewer-phase1-species-random-effects.csv:8` and `results/reviewer-phase1-species-random-effects.csv:15`.

## 2. Mouse drift decomposition

Forty-five of 62 mouse tissue-cell-type units had positive unadjusted Hellinger age slopes, with an exact sign-test p = 0.000497. Leave-one-cell-type-out global effects ranged from +0.0970 to +0.1119 around the full estimate of +0.1028. The largest absolute change was 0.00912 SD. Sources: `results/reviewer-phase1-mouse-drift-summary.csv:2` and `results/reviewer-phase1-mouse-drift-leave-one-cell-type-out.csv:2-77`.

The largest reductions occurred after removing mesenchymal stem cells, podocytes, fenestrated cells, lung fibroblasts, or endothelial cells, giving effects of +0.0970, +0.0972, +0.0976, +0.0981, and +0.0982, respectively. No single cell type removes the global association. The signal is distributed across units, with some larger contributors.

## 3. Nulls, intervals, power, and equivalence

Minimum detectable effects use 80% power and two-sided alpha 0.05. Equivalence was evaluated with both 0.1-SD and 0.25-SD margins, fixed before calculation.

| Null claim | Effect | 95% CI | MDE at 80% power | Equivalent at 0.1 SD | Equivalent at 0.25 SD |
|---|---:|---:|---:|---|---|
| Human vector | +0.0365 | -0.0221 to +0.0951 | 0.0837 | Yes, p = 0.0169 | Yes, p = 1.54e-12 |
| Entropy | -0.0789 | -0.4786 to +0.3207 | 0.5621 | No | No |
| Off-program | +0.2438 | -0.1964 to +0.6840 | 0.6191 | No | No |
| Proliferation expression magnitude | -0.0051 | -0.1882 to +0.1781 | 0.2616 | No | Yes, p = 0.00450 |
| Proliferation fraction magnitude | -0.0033 | -0.0745 to +0.0679 | 0.1017 | Yes, p = 0.00393 | Yes, p = 1.43e-11 |
| NMR vector, clustered continuous-age model | -0.3582 | -0.8978 to +0.1814 | 0.7225 | No | No |
| Zebrafish fine vector | +0.2211 | -0.3588 to +0.8009 | 0.7873 | No | No |

Every value and source row is in `results/reviewer-phase1-nulls-intervals-power-equivalence.csv:2-15`. Only the human vector and proliferation fraction contrast establish equivalence at the narrower margin. Off-program, entropy, NMR, and zebrafish nulls remain imprecise rather than equivalent. The NMR interval here comes from a clustered continuous-age model and supplements, rather than replaces, the original exact group test.

## 4. Interneuron gene changes

MSSM-to-HBCC effect-size correlations among MSSM FDR genes were strong: VIP Pearson r = 0.8946 and Spearman rho = 0.9293; pvalb r = 0.8562 and rho = 0.8419; sst r = 0.9295 and rho = 0.9253. Same-direction concordance was 272 of 294 VIP, 295 of 333 pvalb, and 118 of 123 sst genes. Source rows: `results/reviewer-phase1-interneuron-replication-statistics.csv:5`, `:7`, and `:9`.

The both-cohort FDR overlaps of 24 VIP, 20 pvalb, and 62 sst genes greatly exceed hypergeometric expectations of 0.425, 0.442, and 1.455. The three-subtype overlap contains seven genes. Pairwise MSSM subtype overlaps were 34 VIP-pvalb, 14 VIP-sst, and 22 pvalb-sst genes, all above their chance expectations. Sources: `results/reviewer-phase1-interneuron-subtype-overlap.csv:2-5`.

GO enrichment produced three FDR-significant list-term results: chemical synaptic transmission in VIP and pvalb, and monocyte differentiation in sst. The synaptic terms are biologically coherent with neuronal identity. The small sst immune-annotation overlap should not be treated as evidence of immune conversion. Source: `results/reviewer-phase1-go-summary.csv:2` and `results/reviewer-phase1-interneuron-go-enrichment.csv`.

Dot-plot inputs are `results/reviewer-phase1-interneuron-dotplot-data.csv`. The shared-gene table is `results/h3-human-mssm-shared-neuronal-genes.csv`.

## 5. Fly bidirectional hits

Among 703 pairs of the 38 depth-adjusted vector-hit cell types, 537 pairs shared at least one FDR-significant gene and 480 pairs had hypergeometric q below 0.05. The median Jaccard index was 0.0130 and the largest overlap was 3,122 genes. Thirty-seven of 38 units had at least one gene-level FDR hit. Source: `results/reviewer-phase1-fly-hit-overlap-summary.csv:2`.

Thus the fly hits are not independent at the gene level. They contain shared components. GO enrichment was concentrated in specific units rather than one common pathway across all hits. Five terms passed FDR, led by mitochondrial ATP synthesis coupled electron transport in indirect flight muscle. Source: `results/reviewer-phase1-go-summary.csv:3`.

## 6. Multiplicity inventory and corrected summaries

The complete family-level inventory is `results/reviewer-phase1-multiplicity-inventory.csv`.

Applying one BH family across the core mouse, human, and NMR composition tests changed the original summaries as follows:

- Mouse: seven organs retained at least one FDR-significant fraction change, down from eight. Tongue no longer qualified.
- Human: two cohort datasets retained at least one change, down from three. Kidney Living Donor no longer qualified.
- MSSM retained three changes and integrated lung retained two under the global correction.

Sources: `results/reviewer-phase1-composition-global-bh-by-organ.csv:2-26` and the full test table `results/reviewer-phase1-composition-global-bh.csv`.

The 252 signature-transport p-values were re-corrected without outcome-based deletion. The three survivors and their q-values were unchanged. Donor overlap changes independence wording, not which distinct cell-type hypotheses exist. A deduplicated joint-donor model was not invented. Source: `results/reviewer-phase1-signature-global-bh-recheck.csv`.

## 7. Composition with totals, pseudocounts, and propeller

Total-cell-adjusted models included age, sex when variable, platform for Gao worm, and log total recovered cells. Within-stratum FDR summaries were 15 significant mouse tests across six organs, ten human tests across three cohort datasets, 52 fly tests across head and body, three worm tests, two zebrafish tests, and no NMR tests. Source: `results/reviewer-phase1-composition-total-summary.csv:2-7`. Full mouse and human effect sizes, confidence intervals, p-values, and q-values are in `results/reviewer-phase1-composition-total-adjusted.csv`.

Pseudocount 0.1, 0.5, and 1.0 effects were strongly correlated. Pearson correlations between the 0.1 and 1.0 arms were 0.9569 in mouse, 0.9434 in human, 0.9622 in fly, 0.9995 in worm, and effectively 1 in NMR and zebrafish. The fraction with the same sign in all three arms was 0.9412 in mouse, 0.9313 in human, 0.9444 in fly, and 1.0 in worm, NMR, and zebrafish. Source: `results/reviewer-phase1-composition-pseudocount-summary.csv:2-7`.

The propeller logit plus limma check showed high effect concordance with the CLR models: Spearman rho = 0.9738 in mouse, 0.9313 in human, 0.9711 in fly, 0.9971 in worm, and 1.0 in NMR and zebrafish. Ten mouse results and nine human results passed FDR in both the CLR and propeller models. Full results are in `results/reviewer-phase1-composition-propeller.csv`; the summary is `results/reviewer-phase1-composition-propeller-summary.csv:2-7`.

These models change some individual significance calls. They do not support treating recovered cell counts as direct measurements of absolute tissue composition.

## 8. Worm reference sensitivity

With the day-1 reference, the platform-adjusted global worm effect was +1.6755 SD over 14 days, with platform-stratified permutation p = 0.00250. With the balanced day-6 reference, it was +0.5319, with parametric p = 0.0189 and permutation p = 0.0535. Sources: `results/extended-worm-vector-global-platform-audit.csv:2` and `results/extended-worm-vector-global-platform-audit-reference6.csv:2`.

Day 6 is the primary reviewer-response reference because it includes both sequencing platforms and reduces reference-platform confounding. Both results are reported. The day-6 analysis supports a smaller and less secure age association, not a stable universal effect size.

## 9. P1 documentation and matched null

P1_Clean contains 139 mouse genes, 89 up and 50 down. P0_Raw257 contains 257 genes and remains a sensitivity tier. P1 was derived from the mouse differential-expression analysis before the reported human transport evaluations. Filters remove predicted or pseudogene labels, non-protein-coding biotypes, ribosomal and mitochondrial genes, and specified ambient-RNA cases. Sex-dimorphic and immediate-early genes are not removed from P1; those filters belong to P2 and P3. Exact criteria are implemented in `scripts/run_m5_heterogeneity.py:71-100`.

Mouse-human one-to-one mapping yields 76 up and 42 down genes. Missing genes are omitted rather than imputed. Source: `results/signature-tier-summary.csv:2-3`. PsychAD detected 113 of the 118 mapped genes in each reported cell type, giving coverage 0.9576.

The matched random-signature analysis used 200 expression-decile-matched random signatures per PsychAD cell type, with the correct up-minus-down score and the fully adjusted age model. No cell type survived BH correction against this null; the minimum q was 0.1007. Mean score-depth correlation was 0.1543 for actual P1 scores and 0.3021 for random scores. Source: `results/reviewer-phase1-signature-null-summary.csv:2`; cell-type results are in `results/reviewer-phase1-signature-random-null.csv`.

This null was run in PsychAD, not rerun across every public M9 cohort. Within PsychAD, the reported P1 cell-type associations should therefore be treated as not distinguishable from expression-matched random signatures under this test.

## Deliverables and stop boundary

The supplementary tables and figure-data files listed above are committed with this report. Local verification gave 61 passing tests and one Bayesian-only skip in the core environment; all 146 app-data anchors passed. Compilation and `git diff --check` passed. The manuscript remains unchanged, and nothing was submitted or externally reviewed.
