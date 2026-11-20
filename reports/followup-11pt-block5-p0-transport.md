# Block 5: P0 257-gene unfiltered transport

Status: specification frozen before running. This analysis tests whether the P1 noncoding and related filters change transport verdicts by scoring the full P0_Raw257 set and the P1_Clean set with the same transport implementations.

## Exact analysis specification

- Gene sets are derived by the existing `scripts/run_m5_heterogeneity.py` logic. P0_Raw257 is the discovery contrast DE set with FDR below 0.05 and absolute log2FC at least 0.25. P1_Clean is the existing filtered set. Before transport, `results/m5-gene-level-programs.csv`, `results/m5-program-hierarchy.csv`, and `results/m5-mouse-program-catalog.csv` must verify 257 unique P0 genes and 139 P1 genes, including 89 P1 up and 50 P1 down genes. Any other counts stop the block.
- Ortholog mapping uses the committed one-to-one mouse-human table and rules in `scripts/build_reference_orthologs.py` and `src/aging_programs.py`. The mapped subset is scored. Coverage is reported as N mapped of 257, overall and split into protein-coding versus noncoding or unknown biotype classes.
- Scoring is identical to P1 in every respect except gene-set membership: the composite is the up-gene mean minus the down-gene mean, with the original normalization, donor eligibility, covariates, and per-unit OLS models. Gao uses sex, depth, and six study indicators. Yazar uses sex, depth, and 75-pool fixed effects. M9 uses age and sex, with sex included when variable. M10 uses age, sex, brain bank, PMI, and technical PCs 1 to 3. Donor is the unit of every model.
- The matched random-set comparison uses the existing expression-decile-matched implementation, 2000 iterations, fixed seed 20261105, matching separately on detected up and down set sizes and expression level. P0 matching uses the detected mapped P0 counts, not the P1 counts. Empirical p-values are two-sided around the random-null mean with the plus-one correction.
- The P1 random-set comparison exists only for MSSM, PsychAD, Gao, and Yazar. P0 and P1 empirical tests are run at 2000 iterations in those same units. The other M9 units retain their existing composite transport tests and have no empirical result; they are not assigned a new empirical family or treated as random-null failures.
- For every P1 unit, P0 and P1 are reported side by side with composite slope, composite q, empirical q where defined, pass or fail flags, N donors, and N genes mapped and detected. Global summaries report nominal composite hits, random-set passes among evaluated units, and units passing both tests.
- Shared rule 3 sensitivities refit Gao and Yazar without study or pool batch indicators while retaining sex and depth. These sensitivities are not used to change the primary verdict.

## Cohorts enumerated from existing P1 transport CSVs

- M9: MSSM_DLPFC, HBCC_DLPFC, RADC_DLPFC, SEAAD_DLPFC, Kidney_NormalAdult, Kidney_LivingDonor, Kidney_SN_v15, Kidney_SCR_v15, Lung_CellRef, Lung_Integrated, ImmuneAging_BAL.
- M10: DLPFC_Lifespan_Atlas, also called PsychAD in the P1 matched-null report.
- Follow-up: Gao_Heart and Yazar_Blood.

## P1 families copied before correction

- The primary M9 and M10 composite family is `primary_M9_M10_cell_type_transport`: one BH family across all 252 primary cell-type transport hypotheses, namely 225 M9 tests and 27 M10 tests, without outcome-based deletion. P0 composite q-values use this identical 252-test family.
- The PsychAD matched random-signature family contains its 27 cell types and applies BH across their empirical p-values.
- The MSSM matched random-signature family contains the 21 MSSM cell types with M9 primary tests and applies BH across their empirical p-values.
- Gao Family C contains 14 groups, 13 cell types plus All_LV_Nuclei, with separate BH correction for composite and empirical p-values.
- Yazar pool fixed-effects Family C contains 17 groups, 16 evaluable cell types plus All_Cells, with separate BH correction for composite and empirical p-values across all 17 groups.
- The unadjusted sensitivity families are separate Gao and Yazar slices and do not alter the primary families above.

## Pre-stated verdict rule

PROVES means the filter changed the verdict only if P0 passes both composite and empirical tests after correction in at least one primary unit where P1 fails both. The verdict uses the STRICT reading: fails both means fails the composite test AND fails the empirical test. Units without an empirical result cannot meet this condition. DOES NOT PROVE means no primary unit has P0 passing both tests while P1 fails both tests under this strict reading. INCONCLUSIVE applies if P0 mapping coverage is under 50 percent of 257 genes and P0 fails everywhere, in which case the result cannot distinguish filtering from mapping coverage.

## Results

P0 and P1 were evaluated in 283 primary units each, including 252 M9/M10 units, 14 Gao units, and 17 Yazar units. The primary output also contains 34 unadjusted Gao/Yazar sensitivity rows per gene set.

### Mapping coverage

| Biotype class | Mouse genes | Mapped | Unmapped | Coverage | Mapped up | Mapped down |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| all | 257 | 188 | 69 | 73.15% | 120 | 68 |
| noncoding_or_unknown | 30 | 0 | 30 | 0.00% | 0 | 0 |
| protein_coding | 227 | 188 | 39 | 82.82% | 120 | 68 |

All 30 P0 genes with noncoding or unknown biotypes lack a one-to-one human mapping. The mapped P0 set has 120 up and 68 down human genes, 188 unique genes total. P1 has 76 up and 42 down human genes, 118 total.

### Global comparison

Nominal hit means composite p < 0.05. Composite pass means BH q < 0.05 in the stated composite family. Random pass means empirical BH q < 0.05 among the 79 units per gene set with an existing P1 null family. Both pass requires both q-values below 0.05.

| Gene set | Primary units | Nominal composite hits | Composite passes | Null evaluated | Random passes | Both pass |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| P0_Raw257 | 283 | 48 | 17 | 79 | 15 | 10 |
| P1_Clean | 283 | 47 | 11 | 79 | 16 | 7 |

### Cohort comparison

| Cohort | Gene set | Units | Nominal composite | Composite q passes | Null evaluated | Random q passes | Both pass |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| DLPFC_Lifespan_Atlas | P0_Raw257 | 27 | 6 | 1 | 27 | 0 | 0 |
| DLPFC_Lifespan_Atlas | P1_Clean | 27 | 7 | 1 | 27 | 1 | 0 |
| Gao_Heart | P0_Raw257 | 14 | 6 | 3 | 14 | 8 | 3 |
| Gao_Heart | P1_Clean | 14 | 4 | 1 | 14 | 7 | 1 |
| HBCC_DLPFC | P0_Raw257 | 21 | 4 | 1 | 0 | 0 | 0 |
| HBCC_DLPFC | P1_Clean | 21 | 3 | 0 | 0 | 0 | 0 |
| ImmuneAging_BAL | P0_Raw257 | 29 | 0 | 0 | 0 | 0 | 0 |
| ImmuneAging_BAL | P1_Clean | 29 | 0 | 0 | 0 | 0 | 0 |
| Kidney_LivingDonor | P0_Raw257 | 9 | 1 | 0 | 0 | 0 | 0 |
| Kidney_LivingDonor | P1_Clean | 9 | 2 | 0 | 0 | 0 | 0 |
| Kidney_NormalAdult | P0_Raw257 | 19 | 0 | 0 | 0 | 0 | 0 |
| Kidney_NormalAdult | P1_Clean | 19 | 1 | 0 | 0 | 0 | 0 |
| Kidney_SCR_v15 | P0_Raw257 | 9 | 0 | 0 | 0 | 0 | 0 |
| Kidney_SCR_v15 | P1_Clean | 9 | 0 | 0 | 0 | 0 | 0 |
| Kidney_SN_v15 | P0_Raw257 | 16 | 0 | 0 | 0 | 0 | 0 |
| Kidney_SN_v15 | P1_Clean | 16 | 0 | 0 | 0 | 0 | 0 |
| Lung_CellRef | P0_Raw257 | 26 | 0 | 0 | 0 | 0 | 0 |
| Lung_CellRef | P1_Clean | 26 | 0 | 0 | 0 | 0 | 0 |
| Lung_Integrated | P0_Raw257 | 39 | 7 | 0 | 0 | 0 | 0 |
| Lung_Integrated | P1_Clean | 39 | 5 | 0 | 0 | 0 | 0 |
| MSSM_DLPFC | P0_Raw257 | 21 | 7 | 4 | 21 | 2 | 2 |
| MSSM_DLPFC | P1_Clean | 21 | 11 | 2 | 21 | 3 | 1 |
| RADC_DLPFC | P0_Raw257 | 18 | 3 | 0 | 0 | 0 | 0 |
| RADC_DLPFC | P1_Clean | 18 | 2 | 0 | 0 | 0 | 0 |
| SEAAD_DLPFC | P0_Raw257 | 18 | 7 | 1 | 0 | 0 | 0 |
| SEAAD_DLPFC | P1_Clean | 18 | 4 | 0 | 0 | 0 | 0 |
| Yazar_Blood | P0_Raw257 | 17 | 7 | 7 | 17 | 5 | 5 |
| Yazar_Blood | P1_Clean | 17 | 8 | 7 | 17 | 5 | 5 |

### Side-by-side units with a corrected composite pass or a verdict change

The required CSV contains the complete side-by-side table for every unit. Slopes below are per year for M9/M10 and per decade for Gao/Yazar, matching each source pipeline.

| Cohort | Unit | P0 slope | P0 composite q | P0 empirical q | P0 both | P1 slope | P1 composite q | P1 empirical q | P1 both |
| --- | --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | --- |
| DLPFC_Lifespan_Atlas | EN_L3_5_IT_3 | -0.001401 | 4.777e-05 | 0.1147 | False | -0.001383 | 0.02077 | 0.1957 | False |
| Gao_Heart | All_LV_Nuclei | 0.03071 | 0.04938 | 0.001749 | True | 0.03797 | 0.07721 | 0.006997 | False |
| Gao_Heart | Endocardial | 0.05918 | 0.02705 | 0.004664 | True | 0.08792 | 0.01501 | 0.009329 | True |
| Gao_Heart | Pericyte | 0.04291 | 0.0254 | 0.001749 | True | 0.04003 | 0.1384 | 0.01679 | False |
| HBCC_DLPFC | L6b glutamatergic neuron of the primary motor cortex | 0.001128 | 0.001038 | NA | False | 0.0008703 | 0.1904 | NA | False |
| MSSM_DLPFC | L2/3 intratelencephalic projecting glutamatergic neuron | -0.001236 | 4.231e-06 | 0.04723 | True | -0.001303 | 0.02079 | 0.1854 | False |
| MSSM_DLPFC | L2/3-6 intratelencephalic projecting glutamatergic neuron | -0.0006921 | 0.03484 | 0.5629 | False | -0.001162 | 0.06282 | 0.2399 | False |
| MSSM_DLPFC | oligodendrocyte | 0.0007702 | 0.03484 | 0.04723 | True | 0.001152 | 0.03456 | 0.01049 | True |
| MSSM_DLPFC | oligodendrocyte precursor cell | 0.0008438 | 0.04769 | 0.1484 | False | 0.0009621 | 0.2671 | 0.1154 | False |
| SEAAD_DLPFC | oligodendrocyte precursor cell | 0.002371 | 0.03484 | NA | False | 0.003266 | 0.1155 | NA | False |
| Yazar_Blood | All_Cells | 0.006526 | 0.02457 | 0.3914 | False | 0.005315 | 0.1247 | 0.5992 | False |
| Yazar_Blood | CD14-positive monocyte | 0.01892 | 0.04248 | 0.09204 | False | 0.02252 | 0.06669 | 0.06311 | False |
| Yazar_Blood | central memory CD4-positive, alpha-beta T cell | 0.0175 | 4.042e-05 | 0.02266 | True | 0.02301 | 5.75e-05 | 0.02549 | True |
| Yazar_Blood | effector memory CD8-positive, alpha-beta T cell | 0.01116 | 0.09968 | 0.1456 | False | 0.01643 | 0.0484 | 0.1136 | False |
| Yazar_Blood | memory B cell | 0.03147 | 0.01931 | 0.02266 | True | 0.03792 | 0.03791 | 0.04248 | True |
| Yazar_Blood | naive B cell | 0.02553 | 0.01235 | 0.03228 | True | 0.02861 | 0.03791 | 0.06311 | False |
| Yazar_Blood | naive thymus-derived CD4-positive, alpha-beta T cell | 0.02217 | 5.022e-05 | 0.02266 | True | 0.02692 | 0.0001201 | 0.03823 | True |
| Yazar_Blood | naive thymus-derived CD8-positive, alpha-beta T cell | 0.03205 | 0.02394 | 0.03228 | True | 0.03644 | 0.0484 | 0.04673 | True |
| Yazar_Blood | natural killer cell | 0.008755 | 0.09968 | 0.175 | False | 0.01682 | 0.01503 | 0.04928 | True |

### Both-pass disagreements

The filter changes unit-level both-pass status in 5 of 283 primary units: P0 gains both-pass status in 4 units and loses it in 1 unit. Composite-pass counts are 17 for P0 and 11 for P1. These unit-level changes do not alter the strict verdict below.

| Cohort | Unit | P0 composite q | P0 empirical q | P0 both pass | P1 composite q | P1 empirical q | P1 both pass |
| --- | --- | ---: | ---: | --- | ---: | ---: | --- |
| Gao_Heart | All_LV_Nuclei | 0.0494 | 0.0017 | True | 0.0772 | 0.0070 | False |
| Gao_Heart | Pericyte | 0.0254 | 0.0017 | True | 0.1384 | 0.0168 | False |
| MSSM_DLPFC | L2/3 intratelencephalic projecting glutamatergic neuron | 4.2e-06 | 0.0472 | True | 0.0208 | 0.1854 | False |
| Yazar_Blood | naive B cell | 0.0123 | 0.0323 | True | 0.0379 | 0.0631 | False |
| Yazar_Blood | natural killer cell | 0.0997 | 0.1750 | False | 0.0150 | 0.0493 | True |

### Unadjusted sensitivity

| Cohort | Gene set | Units | Composite q passes | Random q passes | Both pass |
| --- | --- | ---: | ---: | ---: | ---: |
| Gao_Heart | P0_Raw257 | 14 | 5 | 8 | 5 |
| Gao_Heart | P1_Clean | 14 | 0 | 4 | 0 |
| Yazar_Blood | P0_Raw257 | 20 | 14 | 10 | 10 |
| Yazar_Blood | P1_Clean | 20 | 12 | 9 | 9 |

## Verdict

**DOES NOT PROVE**. Under the strict reading, P0 passes both tests in no unit where P1 fails both the composite and the empirical test. P1 fails both tests under this definition in 59 units with an empirical result, and P0 passes both tests in none of them.

This verdict concerns the specified P0/P1 filter comparison, not universal cross-species transport. Empirical results are confined to the cohorts and units that already had P1 matched-null families; no new empirical family was defined for the other M9 units.

## Outputs and reproducibility

- `results/followup-11pt-block5-p0-transport.csv`: complete P0/P1 primary and unadjusted sensitivity comparison.
- `results/followup-11pt-block5-p0-coverage.csv`: mapping coverage overall and by biotype class.
- `scripts/run_followup_11pt_block5_p0_transport.py`: verification, transport, matched nulls, and report generation.
- Fixed seed: 20261105. Random iterations: 2000 per evaluated unit and gene set.
