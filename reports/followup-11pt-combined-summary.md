# Agescape 11-point follow-up: combined summary

This report synthesizes the committed block reports and result CSVs through repository head `c8eb67f`. It adds no analysis and does not change the manuscript. Dataset verdicts are not pooled. All models use donors as the unit, and all q-values are BH-adjusted within the families defined by each block. MDE80 is the smallest detectable effect at 80 percent power. Effects are per decade unless another scale is stated.

## Verdicts

1. **Block 1, residual dispersion:** Gao **DOES NOT PROVE**; Yazar **INCONCLUSIVE**.
2. **Block 2, gene coordination:** Gao **INCONCLUSIVE**; Yazar **INCONCLUSIVE**. Neither dataset proves rising or falling coordination.
3. **Block 3, profile-cloud dispersion:** Gao **INCONCLUSIVE**; Yazar **PROVES** donor-cloud spreading in memory B cells and naive thymus-derived CD8-positive, alpha-beta T cells.
4. **Block 4, reverse transport:** Yazar **DOES NOT PROVE**; Gao **INCONCLUSIVE**.
5. **Block 5, unfiltered P0 transport:** **DOES NOT PROVE** that the P1 filter changed the strict transport verdict.
6. **Block 6, depth extension:** Part A Gao **DOES NOT PROVE**, Part A Yazar **INCONCLUSIVE**; Part B MSSM **PROVES** a direct neuronal depth decline.
7. **Block 7, survivor within-type tests:** Part A **PROVES** stress-cell dependence of the age signal in Yazar, with Gao **DOES NOT PROVE** and Yazar **PROVES**; Part B **DOES NOT PROVE** localization to fraction-losing types.
8. **Block 8, categorical signature models:** **PROVES** transition-concentrated P1 signature effects, with both Gao and Yazar **PROVES**.
9. **Block 9, mature interneuron refit:** **INCONCLUSIVE** because HBCC data are absent, so mature-only cross-cohort replication cannot be tested.
10. **Block 10, published-signature benchmark:** **DOES NOT PROVE** the broad claim that published signatures generally fail batch adjustment.

## New numbers that change or newly constrain a paper claim

| Block | Dataset and result | Verified numbers | Claim consequence | Source CSV |
| --- | --- | --- | --- | --- |
| 3 | Yazar memory B cell cloud | YOUNG 25 donors, mean distance 0.19082; OLD 176, 0.19972; OLD minus YOUNG 0.0089008, q 0.015548; noise floor 0.001778, MDE80 0.0098053 | Direct within-bin donor-cloud spreading qualifies the earlier distance-variance null for this different estimand | `followup-11pt-block3-profile-cloud.csv` |
| 3 | Yazar naive thymus-derived CD8-positive, alpha-beta T cell cloud | YOUNG 77, mean distance 0.16748; OLD 210, 0.18946; difference 0.021979, q 0.0009995; noise floor 0.00084771, MDE80 0.0070891 | A second cell type proves spreading above the split-half floor | `followup-11pt-block3-profile-cloud.csv` |
| 6 | MSSM GABAergic neuron depth | Primary N 155, slope -0.0985, 95 percent CI -0.1402 to -0.0569, q 2.60e-05, MDE80 0.1495 SD; mature N 153, slope -0.1039, p 6.90e-06 | The depth decline survives a direct formal test and the mature-donor sensitivity | `followup-11pt-block6-mssm-depth.csv` |
| 6 | MSSM VIP interneuron depth | Primary N 197, slope -0.0663, CI -0.1044 to -0.0281, q 0.0015, MDE80 0.1329 SD; mature N 193, slope -0.0877, p 3.37e-05 | Direct evidence is not limited to the broad GABAergic unit | `followup-11pt-block6-mssm-depth.csv` |
| 6 | MSSM pvalb interneuron depth | Primary N 189, slope -0.0599, CI -0.1015 to -0.0183, q 0.0050, MDE80 0.1354 SD; mature N 185, slope -0.0839, p 2.74e-04 | The decline persists after maturity in pvalb interneurons | `followup-11pt-block6-mssm-depth.csv` |
| 6 | MSSM sst interneuron depth | Primary N 177, slope -0.0560, CI -0.0910 to -0.0210, q 0.0025, MDE80 0.1435 SD; mature N 174, slope -0.0524, p 0.0064 | All four pre-stated MSSM neuronal units support declining depth | `followup-11pt-block6-mssm-depth.csv` |
| 7 | Yazar central memory CD4-positive, alpha-beta T cell stress dependence | High-stress fraction N 955, slope 0.006311, CI 0.003872 to 0.008751, q 7.484e-06, MDE80 0.042 SD; removal N 948, 10,777 common genes, FDR passes 310 to 232, including 77 predicted attenuated pass-to-fail changes | Removing high-stress cells changes gene-level age verdicts in the predicted direction | `followup-11pt-block7-death-stress.csv`; `followup-11pt-block7-removal-sensitivity.csv` |
| 7 | Yazar naive thymus-derived CD4-positive, alpha-beta T cell stress dependence | High-stress fraction N 821, slope 0.007257, CI 0.004318 to 0.010195, q 1.217e-05, MDE80 0.052 SD; removal N 791, 10,072 common genes, FDR passes 319 to 240, including 83 predicted changes | The age signal is partly dependent on captured high-stress cells | `followup-11pt-block7-death-stress.csv`; `followup-11pt-block7-removal-sensitivity.csv` |
| 7 | Yazar naive thymus-derived CD8-positive, alpha-beta T cell stress dependence | High-stress fraction N 384, slope 0.010752, CI 0.005743 to 0.015761, q 1.675e-04, MDE80 0.075 SD; removal N 349, 8,552 common genes, FDR passes 21 to 12, including 8 predicted changes | A third Yazar cell type meets both parts of the Part A rule | `followup-11pt-block7-death-stress.csv`; `followup-11pt-block7-removal-sensitivity.csv` |
| 8 | Gao Myeloid maturation contrast | N YOUNG/MID/OLD 18/49/39; MID minus YOUNG 2.704e-05, CI 9.503e-06 to 4.457e-05, q 0.028581; OLD minus MID -5.607e-06, q 0.9537; trend q 0.37782; prior linear q 0.1545 | A new maturation-concentrated effect was not flagged by the prior linear model | `followup-11pt-block8-categorical-signature.csv` |
| 8 | Yazar central memory CD4-positive, alpha-beta T cell maturation contrast | N 99/191/660; MID minus YOUNG 1.249e-05, CI 5.535e-06 to 1.944e-05, q 0.0057793; OLD minus MID 5.786e-06, q 0.12347; trend q 1.724e-07; prior linear q 3.382e-06 | The previously linear hit is concentrated in the young-to-mid transition under the pre-stated classification | `followup-11pt-block8-categorical-signature.csv` |
| 8 | Yazar naive thymus-derived CD4-positive, alpha-beta T cell maturation contrast | N 89/163/564; MID minus YOUNG 8.906e-06, CI 2.849e-06 to 1.496e-05, q 0.026064; OLD minus MID 5.625e-06, q 0.12347; trend q 6.496e-07; prior linear q 1.413e-05 | A second previously linear hit is transition-concentrated | `followup-11pt-block8-categorical-signature.csv` |
| 10 | Harries2011_6gene | Both-test passes, unadjusted to adjusted: Gao 0 to 0 of 13; Yazar 7 to 4 of 17; coverage 100 percent in both cohorts | This signature survives batch adjustment | `followup-11pt-block10-benchmark-results.csv`; `followup-11pt-block10-benchmark-signatures.csv` |
| 10 | Sood2015_150probe | Both-test passes: Gao 1 to 1 of 14; Yazar 0 to 0 of 17; coverage 78.38 percent Gao and 84.46 percent Yazar | A second signature survives, in Gao Adipocyte with adjusted slope 0.1807, composite q 0.003978, empirical q 0.006997 | `followup-11pt-block10-benchmark-results.csv`; `followup-11pt-block10-benchmark-signatures.csv` |
| 10 | Peters2015_1497gene | Both-test passes: Gao 0 to 0 of 14; Yazar 17 to 16 of 17; coverage 85.30 percent Gao and 89.38 percent Yazar | A third signature survives, so a majority of the frozen benchmark survives | `followup-11pt-block10-benchmark-results.csv`; `followup-11pt-block10-benchmark-signatures.csv` |
| 10 | deMagalhaes2009_73gene | Both-test passes: Gao 0 to 0 of 14; Yazar 0 to 0 of 17; coverage 90.00 percent Gao and 100 percent Yazar | NO_MAJORITY_CATEGORY because there is no unadjusted both-test pass whose loss could establish a batch artifact | `followup-11pt-block10-benchmark-results.csv`; `followup-11pt-block10-benchmark-signatures.csv` |

Block 7 Part A establishes dependence on captured high-stress cells, not survivor bias as a complete mechanism. Dying cells are under-captured, so removal is a lower bound. Block 8 gives transition classifications precedence when their BH pattern is met, even if an ordinal trend also passes. Block 10 covers four frozen signatures and two primary cohorts; it cannot prove a claim about all published signatures.

## Nulls and scoped checks that protect existing claims

| Block | Protected claim | Verified numbers and scope | Source CSV |
| --- | --- | --- | --- |
| 1 | No general Gao residual-dispersion rise | Zero of 13 cell types PROVES; 12 DOES NOT PROVE and one power-limited INCONCLUSIVE; minimum primary q 0.8885. Twelve nulls have MDE80 at most 0.50 SD, eight at most 0.25 SD. Pooled N 113, slope -0.0013, p 0.1282, MDE80 0.130 SD | `followup-11pt-block1-residual-dispersion.csv` |
| 1 | The Gao Myeloid secondary signal does not overturn the primary null | Gene fraction 0.6885 of 3,833 genes, permutation p 0.0010, q 0.0100, but primary slope 0.0004, p 0.4759, MDE80 0.190 SD. The block report localizes the secondary signal to an ENCODE ambient-RNA artifact | `followup-11pt-block1-gene-fractions.csv`; `followup-11pt-block1-residual-dispersion.csv` |
| 4 | Yazar-derived human-to-mouse transport is not supported | Discovery 14 genes, clean 13, mapped 10 for 76.9 percent coverage, 6 up and 4 down; zero of 51 tested TMS units passes both tests. Unadjusted sensitivity also passes zero of 54. Per-unit MDE80 is 0.3736 to 10.3 per 12 months in the adjusted arm | `followup-11pt-block4-reverse-signature-genes.csv`; `followup-11pt-block4-reverse-transport.csv` |
| 5 | The P1 filter does not change the strict transport verdict | P0 maps 188 of 257 genes, 73.15 percent; protein-coding coverage 188 of 227, 82.82 percent; noncoding or unknown coverage 0 of 30. P0 versus P1: 48 versus 47 nominal composite hits, 17 versus 11 composite passes, 15 versus 16 random passes, and 10 versus 7 both-test passes. P0 gains both-pass status in 4 units and loses it in 1, but passes both in zero of the 59 units where P1 fails both | `followup-11pt-block5-p0-coverage.csv`; `followup-11pt-block5-p0-transport.csv` |
| 6 | No Gao per-cell-type depth decline | All 13 candidate cell types fit; zero negative BH hits and minimum negative-slope q 0.2474. Every MDE80 is at most 0.4780 SD, below the 0.50 power threshold | `followup-11pt-block6-depth-per-cell-type.csv` |
| 7 | No Gao high-stress fraction age signal | Zero of 13 evaluated cell types passes BH; all primary MDE80 values are 0.168 to 0.356 SD, below 0.50 | `followup-11pt-block7-death-stress.csv` |
| 7 | No Gao narrowing localization to fraction-losing types | N 13 types; association -0.03739, CI -0.11536 to 0.04057, p 0.3138, MDE80 0.10898; mean spread slope -0.01430 in 7 losing types versus -0.01307 in 6 gaining types | `followup-11pt-block7-narrowing-localization.csv` |
| 7 | No Yazar narrowing localization to fraction-losing types | N 17 types; association -0.00901, CI -0.02408 to 0.00607, p 0.2224, MDE80 0.02121; mean spread slope -0.00687 in 10 losing types versus -0.00194 in 7 gaining types | `followup-11pt-block7-narrowing-localization.csv` |
| 8 | Not every passing signature effect is transition-concentrated | Yazar naive thymus-derived CD8-positive, alpha-beta T cell is steady trend: N 73/99/208, MID minus YOUNG 7.709e-06, OLD minus MID 5.836e-06, trend 6.655e-06, trend q 0.0087412; prior linear q 0.019929 | `followup-11pt-block8-categorical-signature.csv` |

The Block 4 reverse null is weaker than the mouse-to-human null because human batch noise enters signature derivation. The Block 5 verdict uses the strict rule: P1 must fail both corrected tests in a unit with an empirical result. Unit-level gains by P0 therefore do not reverse that verdict.

## Explicit INCONCLUSIVE register

### Block 1

- **Gao Mast, power:** N 17, primary slope -0.0003, q 0.8885, MDE80 0.590 SD. This exceeds the pre-stated 0.50 SD limit for a powered null.
- **Yazar overall, design and coverage:** One of 16 testable cell types PROVES, 14 DOES NOT PROVE, one is power-limited INCONCLUSIVE, and 12 of 28 are UNTESTABLE. The sole proving type is naive thymus-derived CD4-positive, alpha-beta T cell: primary N 821, slope 0.0005, CI 0.0002 to 0.0007, q 4.80e-04; reweighted slope 0.0007, p 2.78e-05. It has no evaluable internal within-state arm, and the reweighting uses global PBMC composition with residual imbalance. The block report's cell-count and depth diagnostics attenuate the signal; unadjusted, rarefied, and gene-fraction tests do not corroborate it. Pooled All_Cells also has design disagreement: primary slope 0.0001, p 0.4324; within-state 0.0002, p 0.0032; reweighted approximately zero, p 0.9286.
- **Yazar regulatory T cell, power:** N 54, primary slope -0.0009, q 0.6476, MDE80 0.505 SD, just above the powered-null threshold.

### Block 2

- **Gao, power:** Adipocyte, Cardiomyocyte, Epicardial, and Mast remain INCONCLUSIVE because primary MDE80 exceeds the split-half floor: respectively 0.002801 versus 0.006752, 0.002693 versus 0.001707, 0.002165 versus 0.005579, and 0.001005 versus 0.002677. Cardiomyocyte also has a design conflict: primary slope -0.0003919 but OLD minus YOUNG 0.0001989. The dataset summary is 0 PROVES, 9 DOES NOT PROVE, and 4 INCONCLUSIVE.
- **Yazar, power:** CD4-positive, alpha-beta cytotoxic T cell has MDE80 0.001387 versus noise 0.005333; regulatory T cell has 0.002361 versus 0.006972. Memory B cell and mucosal-associated invariant T cell add design conflicts: primary slopes 0.0003965 and 0.0003850, but group differences -0.00008108 and -0.00004467; their primary MDE80 values 0.001235 and 0.001461 are below their respective noise floors 0.0106 and 0.003904.
- **Yazar natural killer cell, technical-floor design:** Primary slope -0.0002595 passes BH at q 0.03493 and the group difference -0.001437 agrees at q 0.01399, but the absolute primary slope is far below the 0.004101 split-half floor. It therefore cannot prove falling coordination. The Yazar summary is 0 PROVES, 11 DOES NOT PROVE, 5 INCONCLUSIVE, and 12 UNTESTABLE.

### Block 3

- **Gao, power:** All 10 evaluated cell types remain INCONCLUSIVE because MDE80 exceeds the split-half floor: Cardiomyocyte, Endocardial, Endothelial, Fibroblast, LEC, Lymphoid, Myeloid, Neuronal, Pericyte, and vSMC. Across these rows MDE80 is 0.021956 to 0.046018, while noise floors are 0.00023679 to 0.0054052. Endocardial, Fibroblast, LEC, and Neuronal also have opposite raw and residualized directions. Three additional cell types are UNTESTABLE on bin floors, not null.
- **Yazar, power:** CD14-low, CD16-positive monocyte; CD14-positive monocyte; central memory CD4-positive, alpha-beta T cell; effector memory CD4-positive, alpha-beta T cell; effector memory CD8-positive, alpha-beta T cell; gamma-delta T cell; naive B cell; naive thymus-derived CD4-positive, alpha-beta T cell; and natural killer cell remain INCONCLUSIVE because their MDE80 exceeds the noise floor. Several have significant negative residualized differences, which do not meet the positive spreading rule and are not converted into powered nulls.
- **Yazar, design:** Central memory CD8-positive, alpha-beta T cell and mucosal-associated invariant T cell have opposite raw and residualized signs and MDE80 above their noise floors. Transitional stage B cell has residualized difference 0.016709, q 0.0009995, above its 0.0011755 noise floor, but raw difference -0.005272; the arm disagreement keeps it INCONCLUSIVE. These 12 inconclusive types coexist with the two proving types and do not cancel the dataset-level PROVES verdict.

### Block 4

- **Gao sensitivity, design:** Adjusted derivation yields 1 discovery gene, 1 clean gene, and 1 mapped gene for 100 percent coverage, but the signature has 0 up genes and 1 down gene. The up-minus-down composite is not estimable, so zero TMS units are testable. This is a missing-direction design limitation, not a mapping null.

### Block 6

- **Part A Yazar, coverage:** No fitted negative slope passes BH, and all 16 fitted candidates have MDE80 at most 0.4448 SD. However, two members of the fixed 18-type candidate family are UNTESTABLE for insufficient residual degrees of freedom: CD16-negative, CD56-bright natural killer cell, human with N 6, and plasmablast with N 14. Their presence prevents DOES NOT PROVE under the frozen rule. Ten additional rows below donor floors are UNTESTABLE, not null.

### Block 9

- **Mature MSSM-to-HBCC replication, coverage:** `/home/reza/projects/agescape-revision/scratch/hbcc.h5ad` is absent. VIP, pvalb, sst, and broad GABAergic replication are therefore all UNTESTABLE. MSSM retains 193, 185, 174, and 153 mature donors respectively, with mature WLS FDR counts 148, 133, 366, and 39, but these one-cohort results cannot estimate cross-cohort correlation, overlap, power, or MDE80.
- **Six shared genes, coverage:** None of FBN1, GABRB1, IGF1R, PHF3, TRMT11, and ZBTB16 can be confirmed as a surviving mature-only cross-cohort shared gene. In supplemental MSSM-only WLS fits, GABRB1, IGF1R, and TRMT11 pass FDR in VIP, pvalb, and sst; FBN1 and PHF3 pass in pvalb and sst but not VIP; ZBTB16 passes only in sst. Every HBCC entry is UNTESTABLE. This neither proves survival from maturity nor proves loss of the program.

All unresolved results above remain distinct from DOES NOT PROVE. Power-limited nulls, coverage failures, missing directional arms, technical-floor failures, and cross-arm design conflicts are not evidence of absence.
