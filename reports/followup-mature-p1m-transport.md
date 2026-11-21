# Mature P1m and P0m core transport at 30-plus

## Scope and runtime

This is Phase 2 Item 2 only. It reruns the original M9 and PsychAD core transport with human-mapped P1m and P0m probe sets, restricts every human donor model to ages at least 30 years, and runs the original 2000-draw expression-decile-matched null. No follow-up cohort or manuscript file was included.

- Model used: `kimi-code/k3.1-pro`.
- Successful composite transport start: `2026-11-21T18:45:45.622485+00:00`.
- Successful composite transport end: `2026-11-21T18:55:54.255071+00:00`.
- Successful composite transport duration: 608.633 seconds, or 10 minutes 8.633 seconds.
- Successful matched-null start: `2026-11-21T19:09:15.521836+00:00`.
- Successful matched-null end: `2026-11-21T19:46:17.851248+00:00`.
- Successful matched-null duration: 2222.329 seconds, or 37 minutes 2.329 seconds.
- PsychAD input: `data/external/dlpfc-atlas/lifespan_pseudobulk.rds`, 539,406,079 bytes, MD5 `aef628e7bc2e5134b790d98c6f57ec48`, verified before both composite and null analyses.

An earlier transport attempt starting at `2026-11-21T18:40:15.549530+00:00` reached Lung Integrated but was terminated by a 300-second tool timeout. It wrote no final outputs. The successful rerun used the same scoring code without a task timeout. Earlier preflight attempts stopped before expression scoring on an invalid P1m-subset assumption, a metadata variable typo, and a duplicate-var assertion; none produced results used here.

The first null attempt started at `2026-11-21T18:57:29.279934+00:00` and stopped before random draws because the strict score check observed a 1.6187524654e-9 floating-point summation difference. The second started at `2026-11-21T19:03:26.160688+00:00` and stopped before random draws on a slope difference subsequently measured as 5.8004291474e-12. Neither attempt wrote a null or comparison result. The successful run used a 1e-7 score tolerance and a 1e-9 absolute slope tolerance. The strict observed-slope agreement check remained in place. The second attempt left an in-repo scratch pseudobulk cache, which the successful run loaded after pair-order verification.

## Probe sets and exact models

The probe sets are the human-mapped lists from Item 1, not the original P0 or P1 sets:

| Probe set | Mouse genes | Mouse up | Mouse down | Human-mapped genes | Human up | Human down |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| P0m | 2598 | 1474 | 1124 | 2049 | 1157 | 892 |
| P1m | 1914 | 1119 | 795 | 1724 | 1013 | 711 |

P1m is not a strict subset of P0m because their DE filtering rules differ. The implementation therefore uses the union of both mapped sets when streaming expression.

The original 252-test family contained 225 M9 cell types and 27 PsychAD cell types. The corresponding mature family retains every original unit that remains eligible, with no outcome-based deletion. It contains 170 tested units per probe set: 143 M9 and 27 PsychAD. Composite BH correction is applied separately to the same 170-unit P0m and P1m families. The other 82 original units are explicitly UNTESTABLE.

- M9 composite score: per-cell log1p(CPM), up-gene mean minus down-gene mean, then the unweighted mean across cells for each donor and cell type. A donor-cell-type pair requires at least 20 cells. The donor model is OLS with intercept and age, plus male sex when both sexes are present. At least 4 donors are required, exactly as in M9.
- PsychAD composite score: log1p(CPM) pseudobulk, up-gene mean minus down-gene mean. The fully adjusted OLS contains intercept, age, male sex, HBCC brain-bank indicator, PMI, and TechPC1, TechPC2, and TechPC3, exactly as in M10.
- All donors in all fitted models have age at least 30. The minimum fitted age in the transport table is 30.
- Local M9 donor-cell-type pairs and cell counts were asserted to match the original M9 score table exactly after applying the 30-plus filter. No original base script was edited.

## Original composite passes under P1m

SURVIVED means P1m composite BH q < 0.05 in the 170-unit mature family. All three original P1 composite passes DIED.

| Cohort | Unit | Original P1 composite q | P1m donors | P1m composite p | P1m composite q | P1m empirical p | P1m empirical q | Verdict |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| MSSM_DLPFC | L2/3 intratelencephalic projecting glutamatergic neuron | 0.0207906814973592 | 195 | 0.3319718865976601 | 0.6882343990439296 | 0.2073963018490754 | 0.2561954316959167 | DIED |
| MSSM_DLPFC | oligodendrocyte | 0.0345576907202604 | 206 | 0.2353471307215863 | 0.6251408159792137 | 0.8835582208895553 | 0.8835582208895553 | DIED |
| DLPFC_Lifespan_Atlas | EN_L3_5_IT_3 | 0.02077125268636786 | 196 | 0.0604830174193466 | 0.3213160300402791 | 0.1184407796101949 | 0.2332405225958449 | DIED |

The original MSSM oligodendrocyte both-pass result also DIED because neither its P1m composite q nor its P1m empirical q is below 0.05.

## Every new P1m composite pass

Seven units pass the P1m composite test, and none is one of the three original P1 composite passes. NA means that the original Block 5 design has no null family for that M9 unit; it does not mean a failed null.

| Cohort | Unit | Donors | Composite p | Composite q | Empirical p | Empirical q |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| DLPFC_Lifespan_Atlas | EN_L6_CT | 193 | 0.00050291216093860001 | 0.021373766839891801 | 0.051974012993503203 | 0.17541229385307339 |
| DLPFC_Lifespan_Atlas | IN_ADARB2 | 196 | 1.8483875674846107e-06 | 0.00031422588647229998 | 0.019990004997501198 | 0.1137288498607839 |
| DLPFC_Lifespan_Atlas | IN_PVALB | 197 | 0.0014777567426313999 | 0.049700865213445397 | 0.34432783608195899 | 0.5468736220125231 |
| DLPFC_Lifespan_Atlas | IN_SST | 197 | 1.029399503823854e-05 | 0.00087498957825019995 | 0.0154922538730634 | 0.1137288498607839 |
| DLPFC_Lifespan_Atlas | Micro | 201 | 8.0604712229384441e-05 | 0.0045676003596650996 | 0.00049975012493749999 | 0.013493253373313301 |
| Lung_Integrated | multiciliated columnar cell of tracheobronchial tree | 71 | 0.0018207876020091 | 0.049700865213445397 | NA | NA |
| MSSM_DLPFC | microglial cell | 196 | 0.0020465062146712001 | 0.049700865213445397 | 0.00049975012493749999 | 0.002623688155922 |

Only PsychAD Micro and MSSM microglial cell pass both tests. The Lung Integrated composite pass has no empirical verdict and is not promoted to both-pass status.

## Matched-null construction

The null implementation is the original procedure with only the required probe and age changes:

1. For each unit and probe set, donor pseudobulk expression is the mean over cells of log1p(CPM) per gene for MSSM and log1p(CPM) of the supplied pseudobulk for PsychAD.
2. Expression bins are deciles of mean logCPM across the included 30-plus donors, computed with `pd.qcut(..., 10, labels=False, duplicates="drop")`.
3. Up and down sets are matched separately. Each draw samples without replacement, within each decile, exactly the detected P1m or P0m number of genes for that direction and unit.
4. The pool excludes every detected gene in the current probe set. There is no additional detection filter, matching the original null.
5. Each tier uses a fresh `np.random.default_rng(20261105)` stream, sorted cell-type order, and 2000 draws. P0m precedes P1m as in the original Block 5 implementation.
6. Random scores use the same unit-specific donor model as the observed score. The empirical p-value is two-sided around the random-null mean with the plus-one correction.
7. BH correction is applied separately within the 21-unit MSSM and 27-unit PsychAD null families for each tier.

The null CSV has 96 rows: 48 units times two probe sets. All rows have 2000 iterations and seed 20261105. MSSM reconstructed scores differed from transport scores only by floating-point summation order; all maximum differences were below 1e-7, and observed slopes agreed within the stated tolerance.

## Block 5 core comparison

This comparison includes only the 48 core units with an original null family: 21 MSSM and 27 PsychAD. It deliberately excludes Gao and Yazar. Therefore its counts are not the original 283-unit or 79-null-unit counts. For context, the original core-only 48-unit comparison had one P1 both-pass, three P1 empirical-only, two P0 both-pass, and zero strict P0 rescues.

| Probe set | Units | Both pass | Empirical only | Composite only | Fail both |
| --- | ---: | ---: | ---: | ---: | ---: |
| P1m | 48 | 2 | 13 | 4 | 29 |
| P0m | 48 | 3 | 15 | 5 | 25 |

Every both-pass and empirical-only unit is listed below with exact p and q values.

| Probe set | Cohort | Unit | Status | Composite p | Composite q | Empirical p | Empirical q |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: |
| P1m | DLPFC_Lifespan_Atlas | Micro | BOTH_PASS | 8.0604712229384441e-05 | 0.0045676003596650996 | 0.00049975012493749999 | 0.013493253373313301 |
| P1m | MSSM_DLPFC | microglial cell | BOTH_PASS | 0.0020465062146712001 | 0.049700865213445397 | 0.00049975012493749999 | 0.002623688155922 |
| P1m | MSSM_DLPFC | GABAergic neuron | EMPIRICAL_ONLY | 0.0089607614370917996 | 0.094857126541885606 | 0.00099950024987499999 | 0.0034982508745627001 |
| P1m | MSSM_DLPFC | L2/3-6 intratelencephalic projecting glutamatergic neuron | EMPIRICAL_ONLY | 0.86182774602728596 | 0.93524885630032617 | 0.0014992503748125 | 0.0044977511244377001 |
| P1m | MSSM_DLPFC | L5/6 near-projecting glutamatergic neuron | EMPIRICAL_ONLY | 0.92877011559870204 | 0.96274951007182519 | 0.0054972513743128002 | 0.012826919873396601 |
| P1m | MSSM_DLPFC | L6 corticothalamic-projecting glutamatergic cortical neuron | EMPIRICAL_ONLY | 0.6468822583792172 | 0.86590538523202298 | 0.00049975012493749999 | 0.002623688155922 |
| P1m | MSSM_DLPFC | L6 intratelencephalic projecting glutamatergic neuron | EMPIRICAL_ONLY | 0.39752340608606751 | 0.72409382046876336 | 0.00099950024987499999 | 0.0034982508745627001 |
| P1m | MSSM_DLPFC | T cell | EMPIRICAL_ONLY | 0.053859744625465403 | 0.30961293291586728 | 0.026986506746626601 | 0.043593587821473802 |
| P1m | MSSM_DLPFC | VIP GABAergic cortical interneuron | EMPIRICAL_ONLY | 0.0036037187826214998 | 0.051739127888596498 | 0.0069965017491254002 | 0.0146926536731634 |
| P1m | MSSM_DLPFC | endothelial cell | EMPIRICAL_ONLY | 0.13617115970679741 | 0.53443315163818217 | 0.00049975012493749999 | 0.002623688155922 |
| P1m | MSSM_DLPFC | lamp5 GABAergic cortical interneuron | EMPIRICAL_ONLY | 0.25954322124039098 | 0.64070065845363333 | 0.0024987506246875999 | 0.0065592203898050001 |
| P1m | MSSM_DLPFC | perivascular macrophage | EMPIRICAL_ONLY | 0.10452138830969231 | 0.45560605160635143 | 0.00049975012493749999 | 0.002623688155922 |
| P1m | MSSM_DLPFC | pvalb GABAergic cortical interneuron | EMPIRICAL_ONLY | 0.056458828943481601 | 0.30961293291586728 | 0.015992003998000999 | 0.030530189450729099 |
| P1m | MSSM_DLPFC | smooth muscle cell | EMPIRICAL_ONLY | 0.6923700456697911 | 0.8794591177170189 | 0.026986506746626601 | 0.043593587821473802 |
| P1m | MSSM_DLPFC | sst GABAergic cortical interneuron | EMPIRICAL_ONLY | 0.153378229053379 | 0.56683258563205319 | 0.031984007996001998 | 0.047976011994002997 |
| P0m | DLPFC_Lifespan_Atlas | IN_SST | BOTH_PASS | 1.2310619467486497e-05 | 0.0010464026547363 | 0.0039980009995002003 | 0.035982008995502197 |
| P0m | DLPFC_Lifespan_Atlas | Micro | BOTH_PASS | 4.2053413770354472e-05 | 0.0023830267803199999 | 0.00049975012493749999 | 0.013493253373313301 |
| P0m | MSSM_DLPFC | microglial cell | BOTH_PASS | 0.0020641100096642001 | 0.037878403118027697 | 0.00049975012493749999 | 0.0011660836248541999 |
| P0m | DLPFC_Lifespan_Atlas | EN_L6_IT_2 | EMPIRICAL_ONLY | 0.031666526077831797 | 0.20705036281659259 | 0.0024987506246875999 | 0.033733133433283297 |
| P0m | MSSM_DLPFC | GABAergic neuron | EMPIRICAL_ONLY | 0.020953633500849299 | 0.16962465214973249 | 0.00049975012493749999 | 0.0011660836248541999 |
| P0m | MSSM_DLPFC | L2/3 intratelencephalic projecting glutamatergic neuron | EMPIRICAL_ONLY | 0.60012826543180198 | 0.85714129149358476 | 0.028485757121439199 | 0.039880059970014897 |
| P0m | MSSM_DLPFC | L2/3-6 intratelencephalic projecting glutamatergic neuron | EMPIRICAL_ONLY | 0.66758355129862235 | 0.88838079339896214 | 0.00049975012493749999 | 0.0011660836248541999 |
| P0m | MSSM_DLPFC | L5/6 near-projecting glutamatergic neuron | EMPIRICAL_ONLY | 0.84525257656669683 | 0.94581195043881683 | 0.0019990004997501002 | 0.0038162736813411001 |
| P0m | MSSM_DLPFC | L6 corticothalamic-projecting glutamatergic cortical neuron | EMPIRICAL_ONLY | 0.70419348139485294 | 0.89337978982929112 | 0.00049975012493749999 | 0.0011660836248541999 |
| P0m | MSSM_DLPFC | L6 intratelencephalic projecting glutamatergic neuron | EMPIRICAL_ONLY | 0.29249994717009919 | 0.7284547463477884 | 0.00049975012493749999 | 0.0011660836248541999 |
| P0m | MSSM_DLPFC | T cell | EMPIRICAL_ONLY | 0.084367125556834302 | 0.44820035452068219 | 0.015992003998000999 | 0.023988005997001498 |
| P0m | MSSM_DLPFC | VIP GABAergic cortical interneuron | EMPIRICAL_ONLY | 0.0102495604789479 | 0.1024956047894794 | 0.00099950024987499999 | 0.0020989505247376002 |
| P0m | MSSM_DLPFC | endothelial cell | EMPIRICAL_ONLY | 0.1297763883197236 | 0.51448191789753139 | 0.00049975012493749999 | 0.0011660836248541999 |
| P0m | MSSM_DLPFC | lamp5 GABAergic cortical interneuron | EMPIRICAL_ONLY | 0.4916451345102486 | 0.78205622109594131 | 0.00049975012493749999 | 0.0011660836248541999 |
| P0m | MSSM_DLPFC | perivascular macrophage | EMPIRICAL_ONLY | 0.1116962345019535 | 0.49969368066663422 | 0.00049975012493749999 | 0.0011660836248541999 |
| P0m | MSSM_DLPFC | pvalb GABAergic cortical interneuron | EMPIRICAL_ONLY | 0.13013366158584611 | 0.51448191789753139 | 0.00049975012493749999 | 0.0011660836248541999 |
| P0m | MSSM_DLPFC | smooth muscle cell | EMPIRICAL_ONLY | 0.69502732660036448 | 0.88838079339896214 | 0.010494752623688101 | 0.018365817091454201 |
| P0m | MSSM_DLPFC | sst GABAergic cortical interneuron | EMPIRICAL_ONLY | 0.17091356854012879 | 0.55875589715042107 | 0.014992503748125901 | 0.023988005997001498 |

### Strict P0m rescue verdict

**DOES NOT PROVE.** P1m fails both tests in 29 of the 48 null-evaluated units. P0m passes both tests in zero of those 29 units. Thus the strict filter-rescue criterion is unmet. Units without an empirical result are not eligible for rescue and are not treated as failures.

## UNTESTABLE units

All four named source objects were checked at their UUID paths under `/home/reza/projects/agescape-revision-data/human-h5ads/` and all four are absent. They account for 71 original cell-type units. The 30-plus donor counts below come from the committed original prethreshold metadata and do not represent new fits.

| Cohort | UUID | Confirmed local status | Original units now UNTESTABLE | Donors at 30-plus | Age range at 30-plus |
| --- | --- | --- | ---: | ---: | --- |
| Kidney_NormalAdult | 0b4a15a7-4e9e-4555-9733-2423e5c66469 | ABSENT | 19 | 18 | 43 to 74 |
| RADC_DLPFC | 5e57cd50-8e42-42d6-940d-5c1660d06864 | ABSENT | 18 | 18 | 64 to 89 |
| SEAAD_DLPFC | 6f7fd0f1-a2ed-4ff1-80d3-33dde731cbc3 | ABSENT | 18 | 20 | 42 to 89 |
| Kidney_SN_v15 | a12ccb9b-4fbe-457d-8590-ac78053259ef | ABSENT | 16 | 8 | 45 to 76 |

The base M9 script maps RADC to UUID `5e57cd50-8e42-42d6-940d-5c1660d06864` and SEAAD to `6f7fd0f1-a2ed-4ff1-80d3-33dde731cbc3`, while the task text pairs those two names with the reverse UUIDs. Both UUIDs were independently checked and are absent, so this label discrepancy does not change any test or family count.

Eleven additional units are UNTESTABLE because fewer than four donors remain at 30-plus with at least 20 cells:

| Cohort | Unit | Donors at 30-plus |
| --- | --- | ---: |
| HBCC_DLPFC | T cell | 0 |
| Lung_CellRef | CD8-positive, alpha-beta cytotoxic T cell | 3 |
| Lung_CellRef | airway submucosal gland collecting duct epithelial cell | 3 |
| Lung_CellRef | classical monocyte | 2 |
| Lung_CellRef | mucus secreting cell of tracheobronchial tree submucosal gland | 3 |
| Lung_CellRef | myoepithelial cell | 3 |
| Lung_CellRef | natural killer cell | 3 |
| Lung_CellRef | plasma cell | 3 |
| Lung_CellRef | serous secreting cell | 3 |
| Lung_Integrated | mucus secreting cell | 3 |
| Lung_Integrated | plasmacytoid dendritic cell | 3 |

No untestable unit has a composite or empirical p-value, and none is silently removed from the accounting.

## Donor counts and power at 30-plus

The table gives unique donors with scores in the original M9 output and donors at 30-plus in the new run. For absent objects, the mature count is metadata-only. PsychAD had 284 donors in the original summary and 201 at 30-plus; each PsychAD unit additionally requires nonzero pseudobulk library size.

| Cohort | Original scored donors | Donors at 30-plus | Tested units | Untestable units |
| --- | ---: | ---: | ---: | ---: |
| MSSM_DLPFC | 213 | 209 | 21 | 0 |
| HBCC_DLPFC | 134 | 64 | 20 | 1 |
| RADC_DLPFC | 18 | 18 | 0 | 18 |
| SEAAD_DLPFC | 21 | 20 | 0 | 18 |
| Kidney_NormalAdult | 18 | 18 | 0 | 19 |
| Kidney_LivingDonor | 19 | 19 | 9 | 0 |
| Kidney_SN_v15 | 10 | 8 | 0 | 16 |
| Kidney_SCR_v15 | 18 | 18 | 9 | 0 |
| Lung_CellRef | 21 | 13 | 18 | 8 |
| Lung_Integrated | 137 | 91 | 37 | 2 |
| ImmuneAging_BAL | 12 | 8 | 29 | 0 |
| DLPFC_Lifespan_Atlas | 284 | 201 | 27 | 0 |

The age restriction plainly reduced donor pools and power in HBCC, Lung CellRef, Lung Integrated, ImmuneAging_BAL, and PsychAD. It destroyed testability for the 11 listed present-source units. ImmuneAging_BAL retains only 8 donors cohort-wide and many units have only 4 donors, so those fits remain low-powered even though the original four-donor rule permits them. MSSM retained 209 of 213 scored donors. No new MDE80 or formal post-hoc power calculation was run, so a non-significant mature result is not an equivalence claim.

## Verification limits and remaining work

- The four absent M9 source objects were not restored, and their biological transport results remain unknown rather than null.
- Donor counts for absent objects are from committed prethreshold metadata because the expression objects cannot be re-read.
- Other M9 units outside MSSM retain composite tests only, exactly as in the original Block 5 design. No new empirical family was invented for them.
- The 30-plus mouse derivation remains age-sex confounded because all 30-month mice are male, as documented in Item 1.
- The numerical tolerances used for MSSM score and slope agreement are disclosed above. All composite p and q values in the deliverables come from the transport fits, and all empirical values come from the completed 2000-draw run.
- No Gao heart, Yazar blood, Perez, HRCA, GTEx, TS2, or other organ P1 transport was run. Those follow-up cohorts, their matched nulls, and the expanded 283-unit Block 5 comparison remain future work. Items 3 and 4 were not started.

## Outputs

- `results/followup-mature-p1m-transport.csv`: 340 tested rows, 170 per probe set, with exact composite p and q values and donor counts.
- `results/followup-mature-p1m-transport-null.csv`: 96 matched-null rows with exact empirical p and q values.
- `results/followup-mature-p1m-block5-p0m-comparison.csv`: 48 side-by-side core units and strict rescue indicators.
- `results/followup-mature-p1m-transport-untestable.csv`: all 82 untestable units with reasons and donor counts.
- `results/followup-mature-p1m-block5-summary.csv`: status counts and strict verdict.
- `results/followup-mature-p1m-donor-summary.csv`: cohort donor and family counts.
- `results/followup-mature-p1m-m9-donor-scores.csv`: thresholded mature M9 donor scores used by the models and null checks.
- `results/followup-mature-p1m-transport-runtime.csv` and `results/followup-mature-p1m-block5-runtime.csv`: successful run times and settings.
- `scripts/run_mature_p1m_transport.py` and `scripts/run_mature_p1m_block5.py`: new implementations; the three base scripts are unchanged.
