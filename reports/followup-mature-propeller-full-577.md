# Propeller on the full pooled mature fraction family (577 tests)

Branch: followup-mature-propeller-full-577, off master feb9be2. Date: 2026-11-25.
Ordered by Reza on the Nov 25 21:16 PT voice call, after he caught that propeller
(the independent second method) had only ever been corrected over the 26 fly and
worm rows that passed the scan. This run corrects propeller once across all 577
tests of the pooled mature family: mouse 114, human 285, fly 162, worm 16.

## Headline numbers

1. Propeller passes 68 of 577 rows at q below 0.05 under the single correction
   q_propeller_pooled_577. The scan passes 39 of 577.
2. Per species, propeller passes: mouse 5 of 114, human 14 of 285, fly 44 of 162,
   worm 5 of 16.
3. Sign agreement between the scan CLR effect and the propeller logit effect:
   478 of 577 (82.8 percent). Mouse 91 of 114 (79.8 percent), human 239 of 285
   (83.9 percent), fly 134 of 162 (82.7 percent), worm 14 of 16 (87.5 percent).
4. Rows passing both the scan and propeller: 32 of 577. Rows that enter under
   propeller only: 36. Rows that drop: 7. Rows unchanged: 534.
5. The 4 headline mammalian rows all still pass propeller under the 577
   correction, so the abstract fraction sentence does not change. Mouse limb
   muscle macrophage is the closest to the line at q 0.0486.

## Method

One propeller p-value was assembled for every one of the 577 scan rows, from
two mature-only runs of the unchanged scripts/run_reviewer_phase1_propeller.R
(limma regression on logit-transformed proportions with a 0.5 pseudocount,
design age plus log total cells plus sex when variable, empirical Bayes on the
age coefficient, donor as the unit for mouse and human, pool as the unit for
fly and worm, mature windows only: mouse 18 months and older, human 30 years
and older, fly 30 days and older, worm 6 days and older):

1. results/followup-mature-propeller.csv: the 399 mammalian scan rows
   (mature inputs from scripts/prepare_mature_propeller.py plus the 9 HRCA
   retina suites from scripts/prepare_mature_propeller_hrca.py).
2. results/followup-mature-fly-worm-propeller.csv: the 178 fly and worm scan
   rows (mature inputs from scripts/prepare_mature_fly_worm_propeller.py).
   Worm rows carry the platform term the R script adds when it varies; fly
   rows have no platform column. This matches the committed fly and worm run
   exactly, so the p-values are reused as-is.

No fresh propeller fits were needed: every scan row matched exactly one row
in exactly one of the two files, with zero missing and zero duplicates, and
the two files do not overlap.

NOT reused: results/reviewer-phase1-composition-propeller.csv (612 rows). It
was built from scripts/prepare_reviewer_phase1_propeller.py, which loads
core_inputs and extended_inputs from scripts/run_reviewer_phase1_composition,
the full-range young-to-old inputs (it includes zebrafish, naked mole rat, and
mouse tissues the mature run skips for fewer than 5 mature donors).
Empirically its p-values differ from the mature runs on 398 of the 577 scan
rows, so it is a different design and is excluded.

The single correction is Benjamini-Hochberg across all 577 propeller
p-values, column q_propeller_pooled_577 in
results/followup-mature-propeller-full-577.csv. The build script is
scripts/build_mature_propeller_full_577.py. The per-cohort BH q from the R
output is kept alongside as propeller_q_within_cohort for traceability.

## Key rows

| Row | Scan effect | Scan q | Propeller effect | Propeller p | Propeller q577 | Propeller pass | Status vs scan |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Headline, MSSM L6b glutamatergic neuron | +0.0198 | 0.004211 | +0.0176 | 0.001093 | 0.01443 | PASS | unchanged |
| Headline, MSSM oligodendrocyte | +0.0173 | 0.004505 | +0.0215 | 5.376e-06 | 0.000282 | PASS | unchanged |
| Headline, MSSM oligodendrocyte precursor cell | -0.0119 | 0.004237 | -0.0157 | 2.888e-08 | 4.166e-06 | PASS | unchanged |
| Headline, mouse limb muscle macrophage | -0.1429 | 0.01243 | -0.0935 | 0.00539 | 0.0486 | PASS | unchanged |
| Disclosed failure, mouse limb muscle mesenchymal stem cell | -0.1005 | 0.008124 | -0.0694 | 0.02559 | 0.1554 | FAIL | dropped |
| Disclosed failure, SEA-AD L5 neuron | -0.0213 | 0.02488 | -0.0261 | 0.00133 | 0.01632 | PASS | unchanged |
| Flagged fly, hemocyte | -0.0110 | 0.00864 | +0.0001 | 0.9603 | 0.9845 | FAIL | dropped |
| Flagged fly, ejaculatory bulb | -0.0161 | 0.004237 | -0.0043 | 0.3009 | 0.6382 | FAIL | dropped |
| Flagged fly, adult ventral nervous system | -0.0081 | 0.03228 | +0.0042 | 0.02891 | 0.168 | FAIL | dropped |
| Flagged fly, posterior midgut | -0.0340 | 0.01752 | -0.0453 | 9.285e-05 | 0.002976 | PASS | unchanged |
| Worm, ARSC | +0.2275 | 0.04168 | +0.3205 | 0.000781 | 0.01288 | PASS | unchanged |

Notes on the key rows:

1. All 4 headline rows pass the 577 correction. Macrophage q 0.0486 is
   marginal: it passes only because the threshold is 0.05.
2. Mesenchymal stem cells drop: scan q 0.00812, propeller q577 0.1554. This
   row was already disclosed as failing propeller within cohort (q 0.0576);
   it now fails the pooled 577 correction as well.
3. SEA-AD L5 passes propeller q577 0.01632. Its disclosed failure is the
   total-cell refit (q 0.0647), not propeller, so its status in the paper
   does not change.
4. Of the 4 flagged fly rows, 3 drop under the 577 correction: hemocyte
   (q577 0.9845), ejaculatory bulb (0.6382), ventral nervous system (0.1680).
   Posterior midgut passes (q577 0.00298); its flag was the sex split
   (female-only), not propeller.
5. Worm ARSC passes (q577 0.01288), same sign, and passes every check that
   exists for it.
6. One new drop beyond the previously flagged rows: fly male accessory gland
   main cell (scan q 0.03364, propeller q577 0.07086).
7. Mammalian scan passes that survive propeller: 10 of 13. The 3 drops are
   mesenchymal stem cells, kidney collecting duct (q577 0.2544), and lung
   multiciliated columnar cells (q577 0.1792). Fly and worm scan passes that
   survive: 22 of 26.

## Every row that enters or drops

43 rows change status versus the scan's 39 passes (36 entered, 7 dropped).
Scan q is q_pooled_mature_bh; propeller q is q_propeller_pooled_577.

| Status | Species | Organ | Cohort | Cell type | Scan q | Propeller p | Propeller q577 |
| --- | --- | --- | --- | --- | ---: | ---: | ---: |
| dropped | Fruit fly | body | AFCA | adult ventral nervous system | 0.03228 | 0.02891 | 0.168 |
| dropped | Fruit fly | body | AFCA | ejaculatory bulb | 0.004237 | 0.3009 | 0.6382 |
| dropped | Fruit fly | body | AFCA | hemocyte_body | 0.00864 | 0.9603 | 0.9845 |
| dropped | Fruit fly | body | AFCA | male accessory gland main cell | 0.03364 | 0.009281 | 0.07086 |
| dropped | Human | kidney | Kidney_LivingDonor | kidney cortex collecting duct epithelial cell | 0.04085 | 0.05883 | 0.2544 |
| dropped | Human | lung | Lung_Integrated | multiciliated columnar cell of tracheobronchial tree | 0.02488 | 0.0344 | 0.1792 |
| dropped | Mouse | Limb_Muscle | TMS | mesenchymal stem cell | 0.008124 | 0.02559 | 0.1554 |
| entered | C. elegans | Whole animal | Gao_GSE229022 | Coelomocyte | 0.9482 | 0.004744 | 0.04415 |
| entered | C. elegans | Whole animal | Gao_GSE229022 | Embroynic_cells | 0.2179 | 0.002456 | 0.02674 |
| entered | C. elegans | Whole animal | Gao_GSE229022 | Germline | 0.1141 | 0.001134 | 0.01443 |
| entered | C. elegans | Whole animal | Gao_GSE229022 | Neuron | 0.8048 | 0.005609 | 0.04903 |
| entered | Fruit fly | body | AFCA | adult glial cell | 0.649 | 0.000176 | 0.004616 |
| entered | Fruit fly | body | AFCA | adult midgut-hindgut hybrid zone | 0.8766 | 0.001029 | 0.01443 |
| entered | Fruit fly | body | AFCA | adult peripheral nervous system | 0.7584 | 1.561e-06 | 0.0001001 |
| entered | Fruit fly | body | AFCA | adult reticular neuropil associated glial cell_body | 0.8695 | 0.0008115 | 0.01301 |
| entered | Fruit fly | body | AFCA | adult tracheal cell | 0.5905 | 0.00579 | 0.04987 |
| entered | Fruit fly | body | AFCA | copper cell | 0.1906 | 0.001992 | 0.02272 |
| entered | Fruit fly | body | AFCA | enterocyte of anterior adult midgut epithelium | 0.1718 | 4.724e-05 | 0.001817 |
| entered | Fruit fly | body | AFCA | enterocyte of posterior adult midgut epithelium | 0.2005 | 0.002008 | 0.02272 |
| entered | Fruit fly | body | AFCA | eo support cell | 0.586 | 0.0001738 | 0.004616 |
| entered | Fruit fly | body | AFCA | epithelial cell_body | 0.2179 | 0.0002443 | 0.005456 |
| entered | Fruit fly | body | AFCA | gustatory receptor neuron | 0.9969 | 0.0002429 | 0.005456 |
| entered | Fruit fly | body | AFCA | intestinal stem cell | 0.8812 | 3.67e-08 | 4.235e-06 |
| entered | Fruit fly | body | AFCA | midgut large flat cell | 0.2179 | 0.0006028 | 0.01054 |
| entered | Fruit fly | body | AFCA | oviduct | 0.1527 | 0.005226 | 0.04787 |
| entered | Fruit fly | body | AFCA | perineurial glial sheath | 0.2179 | 1.841e-05 | 0.0007587 |
| entered | Fruit fly | body | AFCA | subperineurial glial cell_body | 0.095 | 2.693e-06 | 0.0001554 |
| entered | Fruit fly | head | AFCA | cone cell | 0.2971 | 0.0006517 | 0.01106 |
| entered | Fruit fly | head | AFCA | distal medullary amacrine neuron Dm3 | 0.06125 | 0.0008898 | 0.01388 |
| entered | Fruit fly | head | AFCA | dopaminergic PAM neuron | 0.05167 | 0.002367 | 0.02627 |
| entered | Fruit fly | head | AFCA | gamma Kenyon cell | 0.06934 | 0.0004896 | 0.009417 |
| entered | Fruit fly | head | AFCA | neuron 6 | 0.1775 | 0.005877 | 0.04987 |
| entered | Fruit fly | head | AFCA | outer photoreceptor cell | 0.06359 | 0.0002671 | 0.005709 |
| entered | Fruit fly | head | AFCA | subperineurial glial cell_head | 0.2896 | 0.004493 | 0.0425 |
| entered | Human | brain | MSSM_DLPFC | pericyte | 0.1527 | 0.003411 | 0.03453 |
| entered | Human | brain | MSSM_DLPFC | smooth muscle cell | 0.08824 | 0.001018 | 0.01443 |
| entered | Human | lung | Lung_Integrated | pulmonary alveolar type 1 cell | 0.108 | 0.001543 | 0.01817 |
| entered | Human | lung | Lung_Integrated | pulmonary alveolar type 2 cell | 0.07498 | 0.001493 | 0.01795 |
| entered | Human | lung | Lung_Integrated | tracheobronchial goblet cell | 0.2179 | 0.003377 | 0.03453 |
| entered | Mouse | Limb_Muscle | TMS | B cell | 0.1517 | 0.0009667 | 0.01443 |
| entered | Mouse | Limb_Muscle | TMS | skeletal muscle cell | 0.5001 | 0.005546 | 0.04903 |
| entered | Mouse | Lung | TMS | neutrophil | 0.06934 | 0.001139 | 0.01443 |
| entered | Mouse | Marrow | TMS | plasma cell | 0.2901 | 0.0001475 | 0.004256 |

## Why propeller passes more rows than the scan

Propeller passes 68 rows while the scan passes 39, and 36 rows pass propeller
without passing the scan. This is not a data error. Propeller and the scan
test different quantities, and a decomposition refit on the same propeller
inputs (scripts/decompose_propeller_gap.py) shows which difference creates
the signal for the 6 most extreme entered rows (scan p above 0.5):

1. Transform. The scan fits CLR, which centers each type on the geometric
   mean across types, so a type whose share moves with the average gets a CLR
   slope near zero. Propeller fits the own-share logit, which keeps that
   movement. For gustatory receptor neuron the CLR slope is minus 0.00003
   (p 0.994) while the logit slope is plus 0.0115 per day (p 0.0038) with no
   other change. For intestinal stem cells the CLR slope is plus 0.0013
   (p 0.684) while the logit slope is plus 0.0130 (p 6.7e-07). The transform
   alone explains these two.
2. Log total cells. Propeller adds log total cells per donor or pool; the
   scan primary fit does not. For adult reticular neuropil glia the logit
   slope without log total is plus 0.0046 (p 0.71) and with log total is plus
   0.0360 (p 0.0010). For midgut-hindgut hybrid zone the same step moves p
   from 0.284 to 0.0012. For worm coelomocyte and neuron, log total plus the
   platform term moves p from 0.39 and 0.12 to 0.0074 and 0.0028. Total
   recovered cells change with age in these pools, so the adjustment matters.
3. Limma empirical Bayes shrinkage. The committed propeller p-values sit
   close to the plain OLS logit p-values with the same covariates (for
   example gustatory receptor neuron 0.000243 versus 0.000215), so shrinkage
   is a small extra step, not the source of the signal. The one partial
   exception is ARSC, where limma gives p 0.00078 against OLS 0.0048.

The plain CLR refit reproduces the scan effect exactly on these rows (for
example ARSC plus 0.2275 per day), which confirms both fits read the same
underlying counts. The honest statement is that propeller is a different test
of the same question, not a strict subset confirmation of the scan: it passes
rows the scan does not, because it measures the own-share slope and adjusts
for total cells, while the scan measures the slope relative to the
cross-type mean. The paper must present it that way and must not claim
propeller confirms only the 39.

## Headline implication

The abstract fraction sentence does not change. Its rule is that a row must
pass every check, and the 4 headline rows are the only mammalian rows that
do: they pass the scan, the pooled 577 propeller correction, the total-cell
refit, and leave-one-donor-out. The 36 entered rows failed the scan, so they
cannot pass every check. The 3 mammalian drops (mesenchymal stem cells,
kidney collecting duct, lung multiciliated columnar) were already out of the
headline with disclosed failures, and the pooled 577 correction confirms
those failures rather than changing them. No tex headline row changes.

## Verification

1. Independent BH recomputation from the delivered p-value column with
   scipy false_discovery_control matches q_propeller_pooled_577 to a maximum
   absolute difference of 1.1e-16, and gives the same pass count, 68. The
   build script also cross-checked against statsmodels multipletests to
   1.1e-16.
2. All 577 scan rows are present exactly once in the output; the key set is
   identical to the scan; 0 missing, 0 extra.
3. No row is missing a propeller p-value; p-values range from 4.4e-10 to
   0.9995.
4. The status column was recomputed from the scan and propeller pass flags:
   0 mismatches. 39 scan passes minus 7 dropped plus 36 entered equals 68
   propeller passes.
5. The decomposition refit reproduces the scan CLR effect exactly on the
   checked rows (ARSC plus 0.2275 per day), confirming both fits read the
   same counts.

## Files

1. results/followup-mature-propeller-full-577.csv: 577 rows, one per scan
   row, with scan effect and q_pooled_mature_bh, propeller effect, p-value,
   within-cohort q, q_propeller_pooled_577, pass flags, sign agreement, and
   status versus the scan.
2. scripts/build_mature_propeller_full_577.py: assembly and correction.
3. scripts/decompose_propeller_gap.py: the OLS decomposition refits above.
4. This report.
