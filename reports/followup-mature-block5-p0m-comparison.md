# Mature Block 5: P0m versus P1m filter comparison

## Spec

This is the full mature-to-old extension of frozen Block 5. Mouse derivation used exactly 18-versus-30-month animals and every fitted human donor in this report is at least 30 years old. P0m and P1m are the committed human-mapped probe sets, loaded by the same load_programs function used by the mature core implementation. P1m is not a subset of P0m, so their coverage and transport results are evaluated independently.

The committed 170-unit M9/M10 composite family, 48-unit MSSM/PsychAD matched-null comparison, and 82 core untestable dispositions from commit 3f5b57a are reused without refitting or recomputing core nulls. New work is restricted to Gao heart and Yazar blood. Both P0m and P1m primaries are computed fresh. Gao uses sex, depth, and the original six study indicators, at least 4 donors, residual df at least 2, and fresh 30-plus pseudobulk aggregation. Yazar uses the original cache, the same pool-pruning rule, residual df floor 10 for the primary and 2 for the unadjusted arm. Unadjusted arms retain sex and depth only.

Every tested Gao/Yazar unit and probe set has 2000 expression-decile-matched null draws at seed 20261105. Composite and empirical BH correction is separate within each cohort, arm, and probe set. Core composite and empirical corrections remain the committed 170-unit, 21-unit MSSM, and 27-unit PsychAD families. No empirical family was added for other M9 units.

The strict verdict is evaluated across all null-evaluated units: the 48 core units plus tested Gao and Yazar units. P0m must pass both corrected tests where P1m fails both. A unit without an empirical result cannot qualify and is not counted as a P1m failure. UNTESTABLE records an eligibility or input failure and is not a biological null.

## Mapping coverage

| Probe set | Biotype class | Mouse genes | Mapped | Unmapped | Coverage | Mapped up | Mapped down |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| P0m | all | 2598 | 2049 | 549 | 78.87% | 1157 | 892 |
| P0m | noncoding_or_unknown | 289 | 2 | 287 | 0.69% | 0 | 2 |
| P0m | protein_coding | 2309 | 2047 | 262 | 88.65% | 1157 | 890 |
| P1m | all | 1914 | 1724 | 190 | 90.07% | 1013 | 711 |
| P1m | protein_coding | 1914 | 1724 | 190 | 90.07% | 1013 | 711 |

## Global comparison

Primary units are tested units. Original design units also include all explicitly untestable units. Empirical-only, composite-only, and fail-both counts are among null-evaluated units only.

| Probe set | Original design | Primary tested | Untestable | Nominal composite | Composite pass | Null evaluated | Random pass | Both pass | Empirical only | Composite only | Fail both |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| P0m | 283 | 201 | 82 | 34 | 12 | 79 | 29 | 5 | 24 | 5 | 45 |
| P1m | 283 | 201 | 82 | 32 | 9 | 79 | 27 | 4 | 23 | 4 | 48 |

## Cohort comparison

| Cohort | Probe set | Tested | Untestable | Nominal composite | Composite pass | Null evaluated | Random pass | Both pass |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| DLPFC_Lifespan_Atlas | P0m | 27 | 0 | 15 | 7 | 27 | 3 | 2 |
| DLPFC_Lifespan_Atlas | P1m | 27 | 0 | 13 | 5 | 27 | 1 | 1 |
| Gao_Heart | P0m | 14 | 0 | 2 | 0 | 14 | 9 | 0 |
| Gao_Heart | P1m | 14 | 0 | 2 | 0 | 14 | 10 | 0 |
| HBCC_DLPFC | P0m | 20 | 1 | 0 | 0 | 0 | 0 | 0 |
| HBCC_DLPFC | P1m | 20 | 1 | 0 | 0 | 0 | 0 | 0 |
| ImmuneAging_BAL | P0m | 29 | 0 | 0 | 0 | 0 | 0 | 0 |
| ImmuneAging_BAL | P1m | 29 | 0 | 0 | 0 | 0 | 0 | 0 |
| Kidney_LivingDonor | P0m | 9 | 0 | 1 | 0 | 0 | 0 | 0 |
| Kidney_LivingDonor | P1m | 9 | 0 | 0 | 0 | 0 | 0 | 0 |
| Kidney_NormalAdult | P0m | 0 | 19 | 0 | 0 | 0 | 0 | 0 |
| Kidney_NormalAdult | P1m | 0 | 19 | 0 | 0 | 0 | 0 | 0 |
| Kidney_SCR_v15 | P0m | 9 | 0 | 3 | 0 | 0 | 0 | 0 |
| Kidney_SCR_v15 | P1m | 9 | 0 | 4 | 0 | 0 | 0 | 0 |
| Kidney_SN_v15 | P0m | 0 | 16 | 0 | 0 | 0 | 0 | 0 |
| Kidney_SN_v15 | P1m | 0 | 16 | 0 | 0 | 0 | 0 | 0 |
| Lung_CellRef | P0m | 18 | 8 | 0 | 0 | 0 | 0 | 0 |
| Lung_CellRef | P1m | 18 | 8 | 0 | 0 | 0 | 0 | 0 |
| Lung_Integrated | P0m | 37 | 2 | 8 | 2 | 0 | 0 | 0 |
| Lung_Integrated | P1m | 37 | 2 | 7 | 1 | 0 | 0 | 0 |
| MSSM_DLPFC | P0m | 21 | 0 | 3 | 1 | 21 | 15 | 1 |
| MSSM_DLPFC | P1m | 21 | 0 | 3 | 1 | 21 | 14 | 1 |
| RADC_DLPFC | P0m | 0 | 18 | 0 | 0 | 0 | 0 | 0 |
| RADC_DLPFC | P1m | 0 | 18 | 0 | 0 | 0 | 0 | 0 |
| SEAAD_DLPFC | P0m | 0 | 18 | 0 | 0 | 0 | 0 | 0 |
| SEAAD_DLPFC | P1m | 0 | 18 | 0 | 0 | 0 | 0 | 0 |
| Yazar_Blood | P0m | 17 | 0 | 2 | 2 | 17 | 2 | 2 |
| Yazar_Blood | P1m | 17 | 0 | 3 | 2 | 17 | 2 | 2 |

## Side-by-side units with a corrected composite pass or verdict change

Slopes are per year for M9/M10 and per decade for Gao/Yazar. Exact q values are shown.

| Cohort | Unit | P0m slope | P0m composite q | P0m empirical q | P0m both | P1m slope | P1m composite q | P1m empirical q | P1m both |
| --- | --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | --- |
| DLPFC_Lifespan_Atlas | EN_L6_CT | 0.00116086 | 0.024413208917675101 | 0.058470764617691101 | False | 0.00116586 | 0.021373766839891801 | 0.17541229385307339 | False |
| DLPFC_Lifespan_Atlas | IN_ADARB2 | 0.00157392 | 0.00023384516858100001 | 0.058470764617691101 | False | 0.00161241 | 0.00031422588647229998 | 0.1137288498607839 | False |
| DLPFC_Lifespan_Atlas | IN_LAMP5_LHX6 | 0.00099892 | 0.037878403118027697 | 0.078260869565217397 | False | 0.000959667 | 0.051739127888596498 | 0.22938530734632681 | False |
| DLPFC_Lifespan_Atlas | IN_PVALB | 0.000860604 | 0.021117097481106899 | 0.087187175642947701 | False | 0.00073828 | 0.049700865213445397 | 0.5468736220125231 | False |
| DLPFC_Lifespan_Atlas | IN_SST | 0.00122609 | 0.0010464026547363 | 0.035982008995502197 | True | 0.00120136 | 0.00087498957825019995 | 0.1137288498607839 | False |
| DLPFC_Lifespan_Atlas | IN_VIP | 0.000763628 | 0.037878403118027697 | 0.24812593703148419 | False | 0.000684604 | 0.052675447352253303 | 0.71064467766116934 | False |
| DLPFC_Lifespan_Atlas | Micro | -0.0015227 | 0.0023830267803199999 | 0.013493253373313301 | True | -0.00143296 | 0.0045676003596650996 | 0.013493253373313301 | True |
| Lung_Integrated | CD1c-positive myeloid dendritic cell | 0.00696531 | 0.037878403118027697 | NA | NA | 0.00621988 | 0.051739127888596498 | NA | NA |
| Lung_Integrated | multiciliated columnar cell of tracheobronchial tree | 0.00495269 | 0.037878403118027697 | NA | NA | 0.00426844 | 0.049700865213445397 | NA | NA |
| MSSM_DLPFC | microglial cell | -0.000847517 | 0.037878403118027697 | 0.0011660836248541999 | True | -0.000830419 | 0.049700865213445397 | 0.002623688155922 | True |
| Yazar_Blood | All_Cells | -0.0139159 | 8.917229529615622e-10 | 0.0042478760619690152 | True | -0.015066 | 4.0258054397214865e-10 | 0.0042478760619690152 | True |
| Yazar_Blood | CD14-low, CD16-positive monocyte | -0.0145362 | 0.012950914796221558 | 0.0042478760619690152 | True | -0.0154278 | 0.014184904914685619 | 0.0042478760619690152 | True |

## Both-pass disagreements

P0m and P1m disagree on both-pass status in 1 unit. P0m gains status in 1 and loses it in 0.

| Cohort | Unit | P0m composite q | P0m empirical q | P0m both | P1m composite q | P1m empirical q | P1m both |
| --- | --- | ---: | ---: | --- | ---: | ---: | --- |
| DLPFC_Lifespan_Atlas | IN_SST | 0.0010464026547363 | 0.035982008995502197 | True | 0.00087498957825019995 | 0.1137288498607839 | False |

## Unadjusted sensitivity

| Cohort | Probe set | Tested | Untestable | Nominal composite | Composite pass | Random pass | Both pass |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Gao_Heart | P0m | 14 | 0 | 1 | 0 | 12 | 0 |
| Gao_Heart | P1m | 14 | 0 | 1 | 0 | 12 | 0 |
| Yazar_Blood | P0m | 20 | 0 | 5 | 0 | 7 | 0 |
| Yazar_Blood | P1m | 20 | 0 | 5 | 2 | 6 | 2 |

## MDE80 and power

MDE80 per decade is (t_0.975 + t_0.80) times the composite standard error at the fitted residual df. The committed core transport CSV has no standard-error column, so no core MDE80 is fabricated.

| Cohort | Probe set | Arm | Tested | MDE80 median | MDE80 minimum | MDE80 maximum |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| Gao_Heart | P0m | primary | 14 | 0.0321148 | 0.0203577 | 0.0563217 |
| Gao_Heart | P0m | unadjusted | 14 | 0.0533185 | 0.033888 | 0.0828068 |
| Gao_Heart | P1m | primary | 14 | 0.0307497 | 0.0191907 | 0.0578915 |
| Gao_Heart | P1m | unadjusted | 14 | 0.0457309 | 0.0273646 | 0.0660861 |
| Yazar_Blood | P0m | primary | 17 | 0.0103461 | 0.00342667 | 0.040418 |
| Yazar_Blood | P0m | unadjusted | 20 | 0.0123411 | 0.00501832 | 0.173014 |
| Yazar_Blood | P1m | primary | 17 | 0.0109855 | 0.00364875 | 0.0410102 |
| Yazar_Blood | P1m | unadjusted | 20 | 0.0130581 | 0.00533439 | 0.22085 |

| Cohort | Original selected donors | Donors at 30-plus | Minimum fitted age | Original primary tested | Mature primary tested | Mature primary untestable |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Gao_Heart | 117 | 109 | 30 | 14 | 14 | 0 |
| Yazar_Blood | 981 | 932 | 30 | 17 | 17 | 0 |

Gao_Heart loses 8 of 117 selected donors at the 30-plus cutoff (109 retained; fitted unit donor counts range from 22 in Epicardial to 109 in All_LV_Nuclei, Fibroblast, Myeloid, and Pericyte). Primary MDE80 has median 0.0321148 (range 0.0203577 to 0.0563217) per decade for P0m and median 0.0307497 (range 0.0191907 to 0.0578915) per decade for P1m. All 14 primary and 14 unadjusted units remain testable, but the MDE80 range exceeds many observed mature composite slopes (for example 0.004 to 0.034 per decade), especially in smaller lineages such as Epicardial (22 donors, MDE80 0.0563 P0m / 0.0528 P1m), Adipocyte (44 donors, MDE80 0.0495 / 0.0492), and LEC (59 donors, MDE80 0.0559 / 0.0579). Non-significant Gao composite tests therefore reflect limited power as well as smaller mature slopes and are not equivalence claims.

Yazar_Blood loses 49 of 981 selected donors at the 30-plus cutoff (932 retained; primary pool-FE fitted donor counts range from 71 in mucosal-associated invariant T cell to 932 in All_Cells). Primary MDE80 has median 0.0103461 (range 0.00342667 to 0.040418) per decade for P0m and median 0.0109855 (range 0.00364875 to 0.0410102) per decade for P1m. Major blood lineages retain high precision (MDE80 down to 0.0034 to 0.0062 per decade), whereas sparse types remain low-powered: mucosal-associated invariant T cell has primary MDE80 0.0404 (P0m) and 0.0410 (P1m), and in the 20-unit unadjusted sensitivity arm conventional dendritic cell (6 mature donors, residual df 2) reaches MDE80 0.1730 (P0m) and 0.2209 (P1m) per decade.

For the reused core M9/M10 cohorts from commit 3f5b57a, the 30-plus cutoff reduced donor pools in HBCC_DLPFC (64 of 134 donors, 1 untestable), Lung_CellRef (13 of 21 donors, 8 untestable), Lung_Integrated (91 of 137 donors, 2 untestable), ImmuneAging_BAL (8 of 12 donors), and DLPFC_Lifespan_Atlas (201 of 284 donors), while 4 absent M9 objects account for 71 additional untestable units.

## Strict verdict

**DOES NOT PROVE**. Across 79 null-evaluated units, P1m fails both tests in 48. P0m passes both tests in 0 of those units. Core anchors are 29 P1m fail-both units and 0 rescues among 48 units. The additional Gao/Yazar units contribute 19 P1m fail-both units (4 in Gao_Heart and 15 in Yazar_Blood) and 0 rescues. Note that 4 of the 5 P0m both-pass units and all 4 P1m both-pass units (DLPFC_Lifespan_Atlas Micro, MSSM_DLPFC microglial cell, Yazar_Blood All_Cells, and Yazar_Blood CD14-low, CD16-positive monocyte) have negative composite slopes; only DLPFC_Lifespan_Atlas IN_SST (P0m both-pass, P1m composite-only) has a positive slope.

## Survived or died relative to original Block 5

SURVIVED means the same cohort and unit passes both tests with the corresponding mature probe set. DIED means it no longer passes both; an untestable mature unit is labeled DIED_UNTESTABLE and is not treated as a null.

The original full-range design had 10 P0 both-pass and 7 P1 both-pass units among 283, with 0 strict rescues among 59 P1 fail-both units.

| Probe set | Cohort | Unit | Original composite q | Original empirical q | Mature composite q | Mature empirical q | Disposition |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| P0 | Gao_Heart | All_LV_Nuclei | 0.049378552954749903 | 0.0017491254372813 | 0.26909933027856259 | 0.0023321672497084792 | DIED |
| P0 | Gao_Heart | Endocardial | 0.027050635893120099 | 0.0046643344994169002 | 0.701452974381364 | 0.92003998000999498 | DIED |
| P0 | Gao_Heart | Pericyte | 0.025399460578732401 | 0.0017491254372813 | 0.55578269264734736 | 0.014770392581487034 | DIED |
| P0 | MSSM_DLPFC | L2/3 intratelencephalic projecting glutamatergic neuron | 4.231011425369868e-06 | 0.047226386806596701 | 0.85714129149358476 | 0.039880059970014897 | DIED |
| P0 | MSSM_DLPFC | oligodendrocyte | 0.034843060039565402 | 0.047226386806596701 | 0.55875589715042107 | 0.41395968682325501 | DIED |
| P0 | Yazar_Blood | central memory CD4-positive, alpha-beta T cell | 4.0418165659300871e-05 | 0.022655338997168002 | 0.88455369802693473 | 0.97701149425287359 | DIED |
| P0 | Yazar_Blood | memory B cell | 0.0193075525247219 | 0.022655338997168002 | 0.88455369802693473 | 0.62125187406296845 | DIED |
| P0 | Yazar_Blood | naive B cell | 0.0123468831444218 | 0.032283858070964497 | 0.33930060437835513 | 0.1076128602365484 | DIED |
| P0 | Yazar_Blood | naive thymus-derived CD4-positive, alpha-beta T cell | 5.0215665547773583e-05 | 0.022655338997168002 | 0.88455369802693473 | 0.97314978874199265 | DIED |
| P0 | Yazar_Blood | naive thymus-derived CD8-positive, alpha-beta T cell | 0.0239377174854584 | 0.032283858070964497 | 0.88455369802693473 | 0.55950596130506181 | DIED |
| P1 | Gao_Heart | Endocardial | 0.0150077167664998 | 0.0093286689988338004 | 0.45264385921161515 | 0.62568715642178907 | DIED |
| P1 | MSSM_DLPFC | oligodendrocyte | 0.034557690720260402 | 0.010494752623688101 | 0.62514081597921367 | 0.88355822088955527 | DIED |
| P1 | Yazar_Blood | central memory CD4-positive, alpha-beta T cell | 5.7499682662783717e-05 | 0.025487256371812501 | 0.9851244129290041 | 0.93103448275862066 | DIED |
| P1 | Yazar_Blood | memory B cell | 0.037907733086793398 | 0.0424787606196897 | 0.9851244129290041 | 0.74677661169415288 | DIED |
| P1 | Yazar_Blood | naive thymus-derived CD4-positive, alpha-beta T cell | 0.0001201449616442 | 0.038230884557720397 | 0.9851244129290041 | 0.93103448275862066 | DIED |
| P1 | Yazar_Blood | naive thymus-derived CD8-positive, alpha-beta T cell | 0.048399571029992698 | 0.046726636681658801 | 0.83440637382946559 | 0.67753623188405787 | DIED |
| P1 | Yazar_Blood | natural killer cell | 0.015028571789439201 | 0.049275362318840499 | 0.83440637382946559 | 0.74677661169415288 | DIED |

## Untestable Gao and Yazar units (primary and unadjusted)

Zero units in the frozen Gao/Yazar primary rosters (14 Gao, 17 Yazar) and zero units in the Gao/Yazar unadjusted rosters (14 Gao, 20 Yazar) became untestable at 30-plus (130 of 130 rows in `results/followup-mature-block5-gao-yazar-transport.csv` are TESTED).

| Cohort | Arm | Unit | Probe set | Donors at 30-plus | Reason |
| --- | --- | --- | --- | ---: | --- |
| none | primary and unadjusted | none | none | 0 | none |

The 3 sparse Yazar cell types already excluded from the frozen 17-unit primary pool-FE design (evaluated in the 20-unit unadjusted sensitivity arm) also fail the pool-pruning rule at 30-plus: `CD16-negative, CD56-bright natural killer cell, human` (30 mature donors; 6 donors across 3 pools with at least 2 donors; estimated residual df 0 < 10), `conventional dendritic cell` (6 mature donors; 2 donors in 1 pool; estimated residual df -2 < 10), and `plasmablast` (14 mature donors; 4 donors across 2 pools; estimated residual df -1 < 10).

The comparison CSV separately retains all 82 committed core untestable units (71 from 4 absent M9 source objects and 11 from fewer than 4 donors at 30-plus) with their original reasons.

## Runtime and model attribution

- Script author: `kimi-code/k3.1-pro` (`scripts/run_mature_block5_gao_yazar.py`).
- Compute launcher: `claude-opus-5-5-high` (recorded in `results/followup-mature-block5-runtime.csv`).
- Monitoring and independent verification: `gemini-4-neon-xhigh`.
- Run window (UTC): `2026-11-21T21:16:13.935989+00:00` to `2026-11-21T21:33:29.844697+00:00` (total 1035.91 seconds, or 17 minutes 15.91 seconds: `load_inputs` 0.04 s, `gao_30plus_transport` 127.89 s using cached mature pseudobulk, `yazar_30plus_transport` 907.73 s, `comparison_and_verification` 0.21 s, `write_outputs` 0.04 s).

## Verification

In-memory `statsmodels` BH recomputation within every stated Gao/Yazar cohort, arm, and probe-set family gave maximum absolute q discrepancy `0`; independent recomputation from the written CSV (`results/followup-mature-block5-gao-yazar-transport.csv`) gives maximum absolute q discrepancy `1.27675647831893e-15` (float64 decimal round-trip precision). All 130 tested Gao/Yazar transport rows have 2000 null iterations and seed 20261105. Minimum fitted ages are 30 for Gao, 30 for Yazar, and 30 for the reused core M9/M10 units.

The rebuilt full comparison reproduces the committed core values exactly: P1m both-pass 2, P0m both-pass 3, P1m fail-both 29, and strict rescue 0 in 48 core units. All output CSVs, the report, and the script were verified to contain zero em dash or en dash characters.

## Files

- `scripts/run_mature_block5_gao_yazar.py`: new implementation; original and mature core scripts are unchanged.
- `results/followup-mature-block5-gao-yazar-transport.csv`: all Gao/Yazar primary and unadjusted rows, including untestable dispositions.
- `results/followup-mature-block5-p0m-comparison.csv`: complete 283-unit side-by-side comparison.
- `results/followup-mature-block5-coverage.csv`: independent P0m and P1m mapping coverage.
- `results/followup-mature-block5-summary.csv`: global counts and strict verdict.
- `results/followup-mature-block5-runtime.csv`: stage runtimes and null settings.
- `scratch/mature-blocks-458/mature_block5_gao_30plus_inputs.pkl`: large mature Gao pseudobulk cache.
