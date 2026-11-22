# Phase 3 rewrite brief: mature-to-old manuscript

Status: brief only. The manuscript file `report/plos-one/final-plos-one.tex` is still at the
optimized commit 0a022d8 state. This brief is the input to the Phase 3 rewrite job, which
starts after this brief lands. Do not edit the tex file from this brief.

Repo head this brief was built from: c44c2aa ("Phase 2: mature PBMC rerun at 30-plus
(Zhu and Thomson)"), on origin/master. Every number below comes from a committed Phase 2
report or result CSV at that head, except the Item 9b section, which is an explicit
placeholder. Nothing in this brief is submitted anywhere.

Terms used in this brief. Mature means age 30-plus in human donors and 18 months in mice.
Old means age 60-plus in human donors and 30 months in mice. P1 is the original 139-gene
mouse signature derived from differential expression across the full young-to-old range.
P1m is the mature-to-old signature: 18-month versus 30-month mouse differential
expression. P0m is the unfiltered mature-to-old signature. A both-pass unit passes two
tests at once: the composite test (parametric, q below 0.05) and the matched-random
empirical test (empirical q below 0.05). BH q is the Benjamini-Hochberg corrected
q-value. CLR is centered log-ratio, the composition transform used for cell-type
fractions. H2 is the cell-type fraction hypothesis. H3 is the off-program expression
hypothesis; it is never called entropy in the manuscript.

## 1. Frozen decisions the rewrite must obey

- Mature-to-old everywhere. Every analysis in the main text uses mature donors as the
  baseline and old donors as the comparison. Young-reference comparisons move to the
  appendix and never appear in the main text. This is item 3, decided 2026-11-14.
- One coherent story. Present the paper as one experiment. Never narrate the process
  by which conclusions were reached. Never frame any result as a follow-up, a rerun,
  a correction, or a second pass.
- The 624-test fraction family stays primary and locked. 45 passes (17 core plus 28
  added). Both promoted tests (mouse lung vein endothelial, human kidney cortex
  collecting duct) are disclosed as not established with core-only q, total-cell
  adjusted q (0.0627 and 0.4231), and propeller p (0.0528 and 0.0588). The locked
  full-range 624 numbers and disclosures are preserved, not overwritten. The mature
  mammalian composition table is separately declared. Fly, worm, and zebrafish rows
  are never claimed as mature-restricted.
- The matched-null primary analysis is the 2,000-draw Block 5 analysis. The 200-draw
  run is a sensitivity only.
- Item 9 stays in Limitations: zebrafish data is testis from 6 fish, 97.23 percent
  of cells germline.
- Quality bar is state of the art, cutting edge aging biology. Venue is irrelevant
  and stays Reza's call after the rewrite lands.
- Abstract carries p-values, q-values, and Bayesian results with specifics. Checker
  enforces only the 300-word PLOS limit.
- Checker constraints from scripts/check_plos_manuscript.py still apply: abstract at
  most 300 words; no em dashes or en dashes; "Figure 5" and "S1" forbidden; cited
  reference numbers must exactly equal listed numbers; exactly 4 figures with
  existing files; fixed section order; process-residue patterns forbidden; must
  contain the Kavox Labs affiliation, the AWS funding sentence, and the founder
  competing-interests sentence.
- Author fields are settled: affiliation Kavox Labs; funding "no external funding,
  Kavox Labs cloud computing infrastructure (Amazon Web Services)"; competing
  interests "The author is the founder of Kavox Labs".

## 2. Results inventory

One row per analysis. Verdict is the report's own verdict. Numbers are exact as
stated in the cited report; section 6 records which numbers were checked against
the result CSVs in this session.

| # | Analysis | Verdict | Headline numbers | Source files |
| --- | --- | --- | --- | --- |
| 1 | Mature-to-old audit (mandate) | n/a (audit) | P1 was derived from mouse differential expression across the full young-to-old range. Drift against young references is +0.1028 (p 0.000004) but against the 18-month mature reference it collapses to +0.0072 (p 0.67). No mature-to-old signature existed before Phase 2. Interneuron results, fractions, and Block 10 published signatures do not depend on P1. | reports/mature-to-old-audit.md |
| 2 | P1m derivation | n/a (derivation) | P1m 1,914 genes (1,119 up, 795 down), 1,724 human-mapped. P0m 2,598 genes, 2,049 human-mapped. Overlap with old P1: 52 of 139 genes, 27 flipped direction, 1,862 gained. 8 mice total: 4 at 18 months (2 female, 2 male), 4 at 30 months (all male). Only 36 of 151 units eligible. | reports/followup-mature-p1m-derivation.md |
| 3 | P1m transport, core cohorts | DOES NOT PROVE (original transport claim dead) | Core 48-unit comparison: P1m 2 both-pass, 13 empirical-only, 4 composite-only, 29 fail-both. P0m 3 both-pass. Only microglia pass both tests in the core cohorts (PsychAD Micro empirical q 0.0135, MSSM microglial empirical q 0.00262). All original both-pass units died, including Gao endocardial and MSSM oligodendrocyte. 170 of 252 units tested per probe set; 82 untestable (71 absent source files, 11 too few 30-plus donors). | reports/followup-mature-p1m-transport.md; results/followup-mature-block5-gao-yazar-transport.csv |
| 4 | Block 5 strict P0m rescue | DOES NOT PROVE | 79 null-evaluated units. P1m fails both tests in 48. P0m passes both tests in 0 of those 48. P1m both-pass 4 of 79; P0m both-pass 5 of 79. 4 of the 5 P0m both-pass units and all 4 P1m both-pass units have negative composite slopes; only DLPFC_Lifespan_Atlas IN_SST (P0m both-pass, P1m composite-only) has a positive slope. | reports/followup-mature-block5-p0m-comparison.md; results/followup-mature-block5-summary.csv |
| 5 | Composition, mature | PARTIAL | 399 rows: 114 mouse at 18-plus months, 285 human at 30-plus years. 5 of the original 15 survive, all human: Lung_Integrated club cell q 0.0118, multiciliated epithelial cell q 0.0279, MSSM L6b neuron, oligodendrocyte, OPC, all q 0.0118. 9 mouse died on power loss (median SE inflation 3.28x, range 1.30x to 5.53x). Bladder UNTESTABLE. Five new passes, all SECONDARY_NOT_PROMOTED. The two lung survivors do not hold within individual studies (meta-analysis p 0.58 and 0.39). Propeller NOT_VERIFIED. | reports/followup-mature-composition.md; results/followup-mature-composition-*.csv |
| 6 | Drift, mouse and human | DOES NOT PROVE (null) | Mouse depth-adjusted slope -0.0137, p 0.816, 14 donors, 51 units. Human -0.0071, p 0.611, 407 donors, 87 units. Bayesian depth-adjusted +0.0389, 94 percent interval -0.0998 to +0.1773, probability positive 0.7068. | reports/followup-mature-drift.md; results/followup-mature-drift-*.csv |
| 7 | Interneuron replication | PROVES on the unmoderated test only | Both-FDR counts VIP 4 vs 0.60 expected, pvalb 8 vs 1.00, sst 41 vs 11.66, block p 0.001. Mature three-subtype shared-gene intersection is TRMT11 only (original was 6 genes). Moderated intersection UNVERIFIED (recovery gate failed). SEA-AD null. | reports/followup-mature-interneuron.md; results/followup-mature-interneuron-*.csv |
| 8 | Block 1 residual dispersion | Gao DOES NOT PROVE; Yazar INCONCLUSIVE; naive CD4 PROVES | Yazar naive CD4 primary slope 0.00071 per decade, p 5.174e-08, q 8.279e-07, n 777. Debiased slope 0.000313, p 0.000638, still significant, so the verdict survives but the mechanism story is weaker. Pooled Yazar All_Cells flips INCONCLUSIVE to DOES NOT PROVE. Gao Myeloid gene fraction 0.6658, p 0.002999, q 0.02699, still explained as ENCODE ambient-RNA soup. | reports/followup-mature-block1.md; results/followup-mature-block1-*.csv |
| 9 | Block 2 gene coordination | INCONCLUSIVE in both datasets | Gao: 0 PROVES, 6 DOES NOT PROVE, 7 INCONCLUSIVE, 0 untestable. Yazar: 0 PROVES, 6 DOES NOT PROVE, 10 INCONCLUSIVE, 12 untestable. Closest call Yazar natural killer cell: primary slope -0.0003086 per decade, q 0.04182, group old minus mature -0.002035, group q 0.01199, stays INCONCLUSIVE because the slope is below the split-half noise floor 0.004052. Runtime Gao 39.6 min, Yazar 66.3 min. Vs original: 21 survived, 8 died (all does-not-prove to inconclusive on power loss), 12 untestable. | reports/followup-mature-block2-gcl.md; results/followup-mature-block2-*.csv |
| 10 | Block 3 profile cloud | Gao INCONCLUSIVE; Yazar PROVES (2 units) | Yazar memory B diff 0.014282, q 0.02499; naive CD8 diff 0.012476, q 0.001199. Both PROVES survive from the original. Gao 0 PROVES, 9 INCONCLUSIVE, 4 untestable. Contraction pattern partly survives (4 of 6). | reports/followup-mature-block3-profile-cloud.md; results/followup-mature-block3-*.csv |
| 11 | Block 4 reverse transport | DOES NOT PROVE (Yazar); INCONCLUSIVE (Gao) | Human signature tested in mouse. Yazar 0 of 38 units pass both on the 18-plus slope arm and 0 of 17 on the 18-vs-30 arm. Gao 0 testable units. | reports/followup-mature-block4-reverse-transport.md; results/followup-mature-block4-*.csv |
| 12 | Block 6 MSSM depth | Part A INCONCLUSIVE both datasets; Part B PROVES | MSSM neuron depth decline at 30-plus: GABA -0.1039 per decade, q 2.76e-05, N 153; VIP -0.0877, q 6.74e-05, N 193; pvalb -0.0839, q 3.65e-04, N 185; sst -0.0524, q 0.0064, N 174. | reports/followup-mature-block6-depth-extension.md; results/followup-mature-block6-*.csv |
| 13 | Block 7 survivor tests | Part A Gao INCONCLUSIVE, Yazar PROVES, overall PROVES; Part B DOES NOT PROVE | Yazar test 1: 3 BH hits (central memory CD4 slope +0.006713, q 8.269e-05, N 909; naive CD4 +0.008188, q 8.269e-05, N 777; naive CD8 +0.012181, q 2.504e-03, N 344). Predicted pass-to-fail genes total 96 (37, 56, 3). Gao test 1: 2 BH hits (Lymphoid -0.029237, q 3.179e-02, N 56; Pericyte -0.024178, q 3.777e-02, N 88) with 0 predicted genes each. Part B cross-type: Gao slope -0.00537, p 0.9241, N 13 types; Yazar -0.01384, p 0.1241, N 17 types; frozen prediction was positive, both null or opposite. Neither Part A nor Part B alone proves survivor bias as a mechanism. | reports/followup-mature-block7-survivor-tests.md; results/followup-mature-block7-*.csv |
| 14 | Block 8 categorical signature | INCONCLUSIVE (Gao INCONCLUSIVE; Yazar DOES NOT PROVE) | Only passing type: Yazar effector memory CD4 steady trend, trend q 0.0004, mature n 6. All original transition labels died (Gao Myeloid, Yazar central memory CD4, Yazar naive CD4 maturation-concentrated; Yazar naive CD8 steady trend). Gao MDE80 0.77 to 1.14 SD, so Gao is underpowered, not a null. | reports/followup-mature-block8-categorical.md; results/followup-mature-block8-*.csv |
| 15 | Block 10 published signatures | INCONCLUSIVE | Harries2011_6gene SURVIVES: Yazar adjusted both-pass 5 of 17 (original 4); All_Cells slope +0.1985 per decade, composite q 7.511e-47, empirical q 0.004248, N 932. Peters2015_1497gene SURVIVES: Yazar adjusted both-pass 16 of 17; All_Cells slope +0.09675, composite q 2.599e-33, empirical q 0.0006535, N 932. Sood2015_150probe LOST_WITH_BATCH: adjusted both-pass 0 of 14 Gao, 0 of 17 Yazar. deMagalhaes2009_73gene NO_MAJORITY_CATEGORY: 0 both-pass in every arm. The four signatures were not rederived from mature donors; this is a mature-tested sensitivity only. | reports/followup-mature-block10.md; results/followup-mature-block10-*.csv |
| 16 | Concordance | Mouse within-tissue passes at nominal significance; mouse cross-tissue null; human null | Mouse within-tissue: observed mean Pearson 0.2105, joint-null mean 0.1465, p 0.0220, q 0.0440 (mature two-test family) and q 0.0460 (mature six-row mouse family), 59 units, 14 mice. Mouse cross-tissue: 0.0604, p 0.1029, null. Human seven-cohort: 0.1793 vs joint-null mean 0.2126, p 0.8152, null, 127 units, 416 donors. Fly, worm, zebrafish, naked mole-rat NOT APPLICABLE as mature-to-old. | reports/followup-mature-concordance.md; results/followup-mature-concordance-*.csv |
| 17 | Nov 17 CSV tests at 30-plus | Depth nulls survive; survivor-bias null survives | Gao adjusted primary depth slope +0.00563, p 0.604, n 109 of 117. Yazar +0.000433, p 0.742, n 932 of 981. Yazar captured-cell decline survives: slope -0.02919 log cells per decade, p 2.08e-06, about 2.9 percent per decade. Gao exploratory fraction signal died (adjusted p 0.092 vs original 0.0396). Yazar raw donor-spread widening survives (ratio 1.543, Levene p 3.46e-05) but the pool-residualized arm is method-discrepant (Levene p 0.137, bootstrap interval 1.075 to 2.562 excludes 1) and is not claimed. Gao spread null (ratio 1.048, p 0.623). Survivor bias: original nominal mouse signal died (Pearson r +0.250, p 0.0498, N 62, to r -0.057, p 0.689, N 51). New human finding: share loss vs absolute drift Pearson +0.287, p 0.0069, permutation p 0.008, opposite the survivor-bias prediction, Spearman null, ecological correlation only. | reports/followup-mature-nov17-tests.md; results/followup-mature-nov17-*.csv |
| 18 | Organ transports Batch A | No positive-direction pass in any cohort | Gao 0 of 14 both-pass under either signature. Yazar P1m 2 of 17 and P0m 2 of 17, all negative (All_Cells P1m q 4.03e-10, emp q 0.00425, effect -0.0151; monocyte P1m q 0.0142, emp q 0.00425, effect -0.0154). Perez blood controls: multihot 6 of 10, level ge2 2 of 10, unadjusted 9 of 10, all negative. Guilliams liver: 0 of 8 both-pass in every arm under either signature, 9 mature donors. | reports/followup-mature-organ-transports-batch-a.md; results/followup-mature-organ-transports-a-*.csv |
| 19 | Organ transports Batch B (HRCA retina) | Clean null for transport | 0 of 144 both-pass under either signature in 9 suites. All 4 original both-pass units DIED, now EMPIRICAL_ONLY with flipped positive mature effects. Amacrine cell fraction signal survives in 3 snRNA suites. Genes: 1,153,712 tests, only 3 significant. Concordance null in all 9 suites. | reports/followup-mature-organ-transports-batch-b.md; results/followup-mature-organ-transports-b-hrca*.csv |
| 20 | Organ transports Batch C GTEx | Mixed, all Track 1 passes negative; Track 2 Whole Blood positive | Track 1 P1m 5 of 47 both-pass, all NEW, all negative: Artery-Aorta -0.0291, Artery-Tibial -0.0300, Brain-Hypothalamus -0.0432, Esophagus-Muscularis -0.0146, Pituitary -0.0218. P0m 6 of 47 (same 5 plus Brain-Nucleus accumbens -0.0256). All 6 original dual-pass tissues DIED under both tiers. Track 2 Whole Blood both-pass POSITIVE (+0.0119 P1m per decade bin, q 1.29e-22, emp q 0.0015; +0.0122 P0m, q 2.45e-20) and Testis both-pass negative. Muscle-Skeletal DIED both tiers. Bulk caveat stands: bulk shifts conflate within-cell aging with cell-type proportion shifts. | reports/followup-mature-organ-transports-batch-c-gtex.md; results/followup-mature-organ-transports-c-gtex*.csv |
| 21 | Organ transports Batch C small cohorts | RPE survives (negative); one new positive under P0m only; all Thomson PBMC both-pass rows died | 237 rows across 7 cohort CSVs. RPE bulk: the one original both-pass transport that SURVIVED under both tiers, negative direction (P1m effect -0.0211, q 0.0443, emp q 0.0005). Lake kidney scCv3 endothelial: one NEW both-pass under P0m only, POSITIVE (+0.0021 per year, q 0.0402, emp q 0.00064). All 11 original Thomson PBMC both-pass rows DIED. Everything else 0 both-pass. | reports/followup-mature-organ-transports-batch-c.md; results/followup-mature-organ-transports-c-*.csv |
| 22 | Full-depth drift, entropy, identity | DOES NOT PROVE on all 14 primary endpoints | All 14 primary endpoints DOES NOT PROVE in Gao heart and Yazar blood. One exception: Yazar naive CD8 T cells scatter rises with age (SD slope 0.0498, p 5.92e-05, q 4.74e-04, 344 donors, 17,878 cells). Rarefaction robustness 150 of 159 rows pass. Pooled Yazar center shift positive (p 2.99e-05) while within-cell-type result is null (p 0.266), read as a composition effect. | reports/followup-mature-fulldepth.md; results/followup-mature-fulldepth-*.csv |
| 23 | PsychAD M10 | Pooled score rises; per-unit picture mixed with a direction flip | Pooled score rises at 30-plus, 201 donors: P1m +0.000508 per year, p 0.000547, adjusted +0.000371, p 0.0241; P0m +0.000486, p 0.00182, adjusted +0.000385, p 0.0288. 11 of 27 units pass the fully adjusted test per tier, the same 11 in both tiers, 10 positive and Micro negative. Micro flips direction: original +0.002381 per year, mature -0.001786 unadjusted (p 1.61e-08, q 2.17e-07), adjusted -0.001433 (p 8.06e-05, q 0.000725). Caveat: the 10 positive slopes hold only on the log1p-CPM scale; 9 have negative raw-count slopes. IN_SST is a direct sign conflict (CPM +0.001201 adjusted vs raw -0.001704). 200-draw null kills every empirical pass (0 of 27 vs 4 at 2000). | reports/followup-mature-psychad.md; results/followup-mature-psychad-m10-*.csv |
| 24 | PBMC mature | MIXED | P1m transport DOES NOT PROVE, 0 of 35 pass both tests; smallest composite q 0.6299 (Zhu cDC). All 11 original Thomson both-pass rows died. 29 empirical-only passes not claimed (200-draw null miscalibrated, 26 of 35 rows at the 1/201 floor). CLR composition PARTIAL: Thomson CD8 naive falls (CLR -1.319, p 0.000444, q 0.00578; survives every sensitivity except batch, q 0.283); Zhu CD8-Naive shrinks to -0.954, p 0.196; Thomson CD4 Memory passes only the primary model (q 0.0391). Hellinger drift: 4 Zhu types pass, all negative (CD4-Tm -0.0126 q 0.0256, NK-GZMK -0.0137 q 0.0256, cDC -0.00653 q 0.0315, CD4-Naive -0.00993 q 0.0315); Thomson 0 of 13. Targeted genes in CD8 naive T: DOES NOT PROVE, 0 of 12 pass, best q 0.0844. Zhu 10 donors at 30-plus (3 mature, 7 old), Thomson 20 donors (11 mature, 9 old). P1m coverage 1,007 up plus 709 down of 1,724 human-mapped genes (99.5 percent). | reports/followup-mature-pbmc.md; results/followup-mature-pbmc-*.csv |
| 25 | Item 9b between-subject variance | PLACEHOLDER, see section 5 | Not yet committed at head c44c2aa. Yazar variance CSVs exist on beast uncommitted; report not yet written. | results/followup-mature-between-subject-*.csv (uncommitted at brief time) |

## 3. Coherent-story outline

The paper's claims in narrative order, as one experiment. Each claim names its
source files. Nothing here is framed as a follow-up or a reanalysis. The rewrite
must keep this order and this framing in the manuscript itself.

### Claim 1. A mature-to-old signature exists, and it is not the young-to-old signature

The signature is built from mouse differential expression between 18-month mature
animals and 30-month old animals: P1m, 1,914 genes (1,119 up, 795 down), 1,724
human-mapped; unfiltered P0m, 2,598 genes, 2,049 human-mapped. It barely overlaps
the full-range signature: 52 of 139 genes shared, 27 flipped direction, 1,862
gained. Sources: reports/followup-mature-p1m-derivation.md.

### Claim 2. Cell-type fractions change in a few human types, and in no mouse type that survives power

Five of fifteen fraction changes survive, all human: lung club cells and
multiciliated epithelial cells fall, and three MSSM neuron types (L6b, oligodendrocyte,
OPC) change share. Nine mouse fraction changes die on power loss alone, with median
standard-error inflation 3.28x. The two lung survivors do not hold within individual
studies. Sources: reports/followup-mature-composition.md; the locked 624-test family
stays as the fraction record with its disclosures (section 1).

### Claim 3. There is no uniform transcriptional drift from mature to old

Mouse drift is null (depth-adjusted slope -0.0137, p 0.816, 14 donors, 51 units) and
human drift is null (-0.0071, p 0.611, 407 donors, 87 units). The Bayesian model
agrees (+0.0389, 94 percent interval -0.0998 to +0.1773, probability positive 0.7068).
Sources: reports/followup-mature-drift.md.

### Claim 4. The mature signature does not transport across tissues, organs, or species in the predicted direction

In the core cohorts, only microglia pass both tests, and the strict unfiltered
rescue fails: P0m passes both in 0 of the 48 units where P1m fails both. Reverse
transport fails (0 of 38 Yazar units on the slope arm, 0 of 17 on the 18-vs-30 arm;
Gao untestable). Organ transports show no positive-direction pass in heart, blood,
or liver; Yazar and Perez both-pass units all score lower with age. Retina is a
clean null (0 of 144). GTEx Track 1 passes are all new and all negative; the six
original dual-pass tissues all die. RPE is the one original both-pass transport
that survives, in the negative direction. All eleven Thomson PBMC both-pass rows
die. Sources: reports/followup-mature-p1m-transport.md,
reports/followup-mature-block5-p0m-comparison.md,
reports/followup-mature-block4-reverse-transport.md,
reports/followup-mature-organ-transports-batch-a.md,
reports/followup-mature-organ-transports-batch-b.md,
reports/followup-mature-organ-transports-batch-c-gtex.md,
reports/followup-mature-organ-transports-batch-c.md,
reports/followup-mature-pbmc.md.

### Claim 5. Two signals survive scrutiny, and both are narrow

First, naive CD8 T cells. Their within-cell scatter rises with age in Yazar blood
(SD slope 0.0498, p 5.92e-05, q 4.74e-04, 344 donors, 17,878 cells), the one
exception among fourteen full-depth primary endpoints. Their donor cloud spreads
(Block 3 naive CD8 diff 0.012476, q 0.001199, and memory B 0.014282, q 0.02499).
Their share falls in Thomson PBMC (CLR -1.319, p 0.000444, q 0.00578). GTEx whole
blood is the one positive Track 2 transport (+0.0119 per decade bin, q 1.29e-22).
Second, MSSM neuron depth falls with age: GABA -0.1039, VIP -0.0877, pvalb
-0.0839, sst -0.0524 per decade, all q below 0.01. Sources:
reports/followup-mature-fulldepth.md, reports/followup-mature-block3-profile-cloud.md,
reports/followup-mature-pbmc.md,
reports/followup-mature-organ-transports-batch-c-gtex.md,
reports/followup-mature-block6-depth-extension.md.

### Claim 6. The interneuron pattern replicates on the unmoderated test, with one shared gene

VIP, pvalb, and sst interneuron subtypes replicate with both-FDR counts 4 vs 0.60,
8 vs 1.00, and 41 vs 11.66 expected, block p 0.001. The three-subtype shared-gene
intersection is TRMT11 only. Sources: reports/followup-mature-interneuron.md.

### Claim 7. Published aging signatures mostly survive batch adjustment, so the benchmark does not support a broad artifact claim

Harries 6-gene and Peters 1,497-gene signatures survive at 30-plus in Yazar blood.
Sood 150-probe loses with batch adjustment. deMagalhaes 73-gene has no majority
category. The scoped statement only is supported: some published signatures are
batch-sensitive; the broad statement that most aging signatures are batch or
pooling artifacts is not supported by this benchmark. Sources:
reports/followup-mature-block10.md.

### Claim 8. No general rise in disorder, coordination loss, or survivor-bias signal

Gene coordination is inconclusive in both datasets (Block 2: 0 proves in Gao and
Yazar). Categorical signature patterns are inconclusive (Block 8: one fragile
steady-trend label on 6 mature donors). Survivor tests do not prove survivor bias
as a mechanism (Block 7 Part B null in both datasets; Part A overall proves only
on the Yazar arm with 96 predicted pass-to-fail genes). The full-depth disorder
endpoints are all null except naive CD8 scatter. Depth per cell is null in both
datasets; only captured-cell counts fall in Yazar (about 2.9 percent per decade).
The original nominal mouse survivor-bias correlation dies at 30-plus. Sources:
reports/followup-mature-block2-gcl.md, reports/followup-mature-block8-categorical.md,
reports/followup-mature-block7-survivor-tests.md,
reports/followup-mature-fulldepth.md, reports/followup-mature-nov17-tests.md.

### Claim 9. Cross-tissue concordance is weak and species-specific

Mouse within-tissue concordance passes at nominal significance (mean Pearson
0.2105, p 0.0220, q 0.0440), carried by limb muscle and marrow. Mouse
cross-tissue is null (0.0604, p 0.1029). Human concordance is null (0.1793 vs
joint-null mean 0.2126, p 0.8152). Shared cell-type names across species are not
proven evolutionary matches; the paper says "no conserved direction", not proof.
Sources: reports/followup-mature-concordance.md.

### Limitations the story must carry

Power loss from the 30-plus restriction kills mouse fraction signals and leaves
Gao inconclusive in Blocks 2, 3, and 8. 82 of 252 transport units are untestable.
The P1m derivation has a sex confound (section 4). The PsychAD pooled rise holds
only on the log1p-CPM scale and microglia flip direction (section 4). The
moderated interneuron intersection is unverified (section 4). Zebrafish testis is
97.23 percent germline (item 9, Limitations). The De Man human lung entropy
replication stays blocked (IPF controls done, 28 subjects, no q below 0.05; the 32
Baylor lungs in GSE281219 have no public barcode-level labels) and stays a
limitation, not a blocker. The LACA NMR atlas verification is closed as
unverified (download link is a placeholder on both hosts).

## 4. Open questions the rewrite must handle honestly

1. P1m is 14x larger than P1 under the same derivation rules (1,914 genes vs 139).
   The size jump is unexplained in the derivation report. The rewrite must not
   invent an explanation. State the numbers, state that the cause is not
   established, and do not let the size difference pass without comment.
2. P1m derivation sex confound. All four 30-month mice are male; the 18-month
   group is 2 female and 2 male. Any mature-to-old mouse difference is tangled
   with sex. The rewrite must disclose this wherever P1m is introduced and must
   not present P1m as a clean age contrast.
3. Moderated interneuron intersection UNVERIFIED. The frozen recovery gate failed,
   so only the unmoderated TRMT11 claim is usable. The rewrite must claim TRMT11
   only, must not present a moderated intersection, and must say the moderated
   check did not verify.
4. Mouse concordance p 0.0220 is nominal until the rebuilt six-test family is
   applied. The report gives q 0.0440 under the mature two-test family and q
   0.0460 under the mature six-row mouse family. The rewrite must report the q
   values it relies on and must not present the p value alone as established.
5. Propeller checks NOT_VERIFIED. No Rscript, speckle, or limma exists on beast,
   so the composition propeller sensitivity was never run. The rewrite must not
   claim propeller support and must say the check is not verified.
6. Batch A runtime CSV model_used column wrongly reads "Gemini 4 Neon (XHigh)".
   The value is hardcoded in the script from the original run. The Yazar rerun
   and the Perez and Guilliams fixed reruns were executed directly on beast, not
   by Gemini. The rewrite must not cite the runtime label as the runner of the
   reruns, and any methods text must describe direct execution.
7. PsychAD IN_SST sign conflict. IN_SST is positive on the log1p-CPM scale
   (+0.001201 adjusted) and negative on the raw-count scale (-0.001704). Only
   Micro holds on both scales. The rewrite must disclose the scale dependence
   and must not present the PsychAD pooled rise as scale-free.
8. PBMC 200-draw null miscalibration. The 200-draw null puts 26 of 35 rows at
   the 1/201 floor, so the 29 empirical-only passes are not claimed. The rewrite
   must not cite any PBMC empirical-only pass, and must rest PBMC transport on
   the composite test only (0 of 35 both-pass).
9. Locked full-range 624 numbers and disclosures must be preserved, not
   overwritten. The mature composition table is a separate declaration. The
   rewrite must keep the 45-pass 624 record, the promoted-test sensitivity
   failures, and the Limitations disclosure exactly as locked, and must not mix
   mature counts into the 624 family.
10. TS2 prose slip, unverified in this session. A prior session noted the Batch C
    report says six empirical passes where the file shows four. The rewrite agent
    must check results/followup-mature-organ-transports-c-ts2.csv before citing
    any TS2 empirical-pass count.
11. Item 9b numbers are pending (section 5). The rewrite must not cite
    between-subject variance numbers until the Item 9b report lands and is
    verified.

## 5. PLACEHOLDER: Item 9b between-subject numbers (TO BE FILLED WHEN THAT COMMIT LANDS)

Status at brief time: NOT COMMITTED. Repo head is c44c2aa. The Item 9b report
file does not exist yet in the repo. Yazar variance CSVs exist on beast
uncommitted at results/followup-mature-between-subject-variance-yazar.csv
(220 rows) and results/followup-mature-between-subject-variance-yazar-ineligible.csv
(35 rows), plus Gao equivalents, donor files, and a Yazar A2 parallel pkl. The
Item 9b agent is writing the report and committing it separately.

FILL IN HERE WHEN THE ITEM 9B COMMIT LANDS (do not fill from the uncommitted
CSVs; fill only from the committed report after verifying its headline numbers
against the committed CSVs):

- Item 9b verdict (Gao / Yazar):
- Yazar variance headline numbers (ratio, Levene p, bootstrap interval, N):
- Gao variance headline numbers (ratio, Levene p, bootstrap interval, N):
- Ineligible counts and reasons:
- How Item 9b changes or leaves unchanged the Nov 17 donor-spread numbers in
  inventory row 17:
- Source files (committed report path and CSV paths):
- Commit hash:

## 6. Verification log

Checked by the brief author against the result CSVs in this session:

- Block 5 summary CSV: P1m both-pass 4, fail-both 48, strict rescue 0, verdict
  DOES NOT PROVE; P0m both-pass 5. Verified.
- Block 6 MSSM depth CSV: GABA slope -0.1039, q 2.76e-05, N 153. Verified.
- PsychAD cohort summary CSV: pooled P1m slope +0.0005083 per year, p 0.000547,
  adjusted +0.000371, p 0.0241, N 201. Verified.
- PBMC Thomson CLR CSV: CD8 Naive CLR -1.3189, p 0.000444. Zhu CLR CSV:
  CD8-Naive -0.9542, p 0.1961. Verified.
- Full-depth report line 89: Yazar naive CD8 scatter SD slope 0.0498023,
  p 5.920e-05, q 4.736e-04, donors 344, cells 17,878. Verified in the report;
  the value 0.0498 is present in results/followup-mature-fulldepth-yazar-spread.csv.
- Composition human CSV: Lung_Integrated club cell q 0.0118; multiciliated
  epithelial cell q 0.0279. Verified.
- Block 2 report read in full this session; Yazar natural killer cell numbers
  match the report table. Block 3 report read in full this session; Yazar
  memory B and naive CD8 numbers match the report table. Block 4 verdict table
  read in full. Block 6 verdicts read in full. Batch A, Batch B, Batch C GTEx,
  Batch C, full-depth, PsychAD, PBMC verdict sections read in full this session.

Checked by reader subagents against the result CSVs (their handoffs on file
with the parent agent; the brief author read the underlying reports directly
for the headline rows listed above):

- Block 7: Lymphoid slope and q, effector memory CD4 MID row verified in CSVs.
- Block 8: effector memory CD4 MID row verified in CSV.
- Block 10: Peters Yazar All_Cells row verified in CSV.
- Concordance: mouse within-tissue row verified in CSV.
- Nov 17 tests: depth rows verified in CSV.

Not independently re-verified cell by cell: every per-unit value quoted from a
report beyond the headline rows listed above. Those numbers are as stated in
their cited reports. The rewrite agent must re-check any number it promotes to
the abstract, the main-text claims, or a figure against the cited CSV before
using it.

Numbers the brief author could not verify: none of the headline numbers above
lacked a supporting CSV row or report table in the checks performed. The TS2
empirical-pass count (open question 10) is the one known unverified discrepancy
and is flagged, not resolved, here.

## 7. What the rewrite job must deliver

- Rewrite report/plos-one/final-plos-one.tex as one coherent mature-to-old story
  following section 3, with young-reference material appendix-only.
- Checker PASS, Tectonic PDF rebuild outside the repo, commit, push, tag.
- Issue-by-issue accounting in reports/mature-to-old-rewrite.md covering every
  open question in section 4 and the Item 9b fill-in from section 5.
- Refresh the submission package at ~/workspace/agescape-final/plos-one-rev2/
  from the new commit (it currently matches 0a022d8 and goes stale the moment
  the rewrite commits).
- Nothing is submitted anywhere. Submission needs Reza's explicit instruction
  AND his final preview check first.
