# Complete mature composition refit

## Scope and model

This is Phase 2 Item 3 only. Mouse fits use ages 18, 21, 24, and 30 months where present. Human fits use donors aged at least 30 years. Naked mole-rat rows are not eligible for a mature refit and retain their original p-values only for the 431-family correction.

The model is unchanged from run_organ_composition.py: donor by cell-type counts, pseudocount 0.5, centered log-ratio transformation, OLS on age plus sex, and age-only fallback when sex makes the design rank-deficient. A cell type is eligible when observed in at least 6 mature donors. Mouse effects are CLR slopes per month and human effects are CLR slopes per year, despite the inherited effect_type label CLR_slope_per_year.

The 624-family locked result was not used as an input or modified. Its SHA256 before and after this run is checked against 3f26242585017bd3a68a4a7c4dcb1d2684da0e9f662c37ed7e07e1bd6f054d96.

## Multiplicity and coverage

Mouse tested rows: 114. Human tested rows: 285. Total mature refit rows: 399.
q_431_family uses one BH correction over all 431 original primary tests after substituting mature p-values only at matching mouse or human keys. It is defined only for rows belonging to that original key set.
q_refit_family uses a separate BH correction over all 399 mature mouse and human tests actually run. The two q columns are never mixed. Within-stratum q_bh is also retained for continuity with the original outputs.

### Cohort coverage

| Species | Organ | Cohort | Mature donors | Tested cell types | Minimum age | Maximum age |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| Human | brain | HBCC_DLPFC | 64 | 24 | 30 | 85 |
| Human | brain | MSSM_DLPFC | 209 | 25 | 31 | 89 |
| Human | brain | RADC_DLPFC | 18 | 23 | 64 | 89 |
| Human | brain | SEAAD_DLPFC | 20 | 18 | 42 | 89 |
| Human | kidney | Kidney_LivingDonor | 19 | 17 | 33 | 65 |
| Human | kidney | Kidney_NormalAdult | 18 | 26 | 43 | 74 |
| Human | kidney | Kidney_SCR_v15 | 18 | 23 | 30 | 66 |
| Human | kidney | Kidney_SN_v15 | 8 | 19 | 45 | 76 |
| Human | lung | ImmuneAging_BAL | 8 | 31 | 33 | 75 |
| Human | lung | Lung_CellRef | 13 | 29 | 35 | 76 |
| Human | lung | Lung_Integrated | 91 | 50 | 30 | 81 |
| Mouse | Fat | TMS | 6 | 5 | 18 | 30 |
| Mouse | Heart_and_Aorta | TMS | 9 | 8 | 18 | 30 |
| Mouse | Kidney | TMS | 11 | 19 | 18 | 30 |
| Mouse | Limb_Muscle | TMS | 12 | 9 | 18 | 30 |
| Mouse | Liver | TMS | 7 | 6 | 18 | 30 |
| Mouse | Lung | TMS | 10 | 29 | 18 | 30 |
| Mouse | Marrow | TMS | 13 | 18 | 18 | 30 |
| Mouse | Pancreas | TMS | 8 | 8 | 18 | 30 |
| Mouse | Spleen | TMS | 10 | 12 | 18 | 30 |
| Mouse | Bladder | TMS | 4 | 0 | 18 | 24 |
| Mouse | Large_Intestine | TMS | 3 | 0 | 30 | 30 |
| Mouse | Mammary_Gland | TMS | 4 | 0 | 18 | 21 |
| Mouse | Skin | TMS | 3 | 0 | 18 | 21 |
| Mouse | Thymus | TMS | 5 | 0 | 18 | 24 |
| Mouse | Tongue | TMS | 4 | 0 | 18 | 24 |
| Mouse | Trachea | TMS | 0 | 0 | NA | NA |

## Original 15 global passes

SURVIVED means q_431_family < 0.05. DIED means the row was estimable but did not clear that bar. UNTESTABLE is not a null. Both mature q families are shown.

| Species | Organ | Cohort | Cell type | Original q | Mature effect | SE | SE inflation | p | q_431_family | q_refit_family | Donors original to mature | Verdict |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Mouse | Bladder | TMS | bladder cell | 0.03752322857 | NA | NA | NA | NA | NA | NA | 8 to NA | UNTESTABLE |
| Mouse | Kidney | TMS | kidney capillary endothelial cell | 0.00441905463 | -0.1567352404 | 0.0692955729 | 4.087235677 | 0.05356862952 | 0.3944180787 | 0.3989428543 | 16 to 11 | DIED |
| Mouse | Kidney | TMS | kidney mesangial cell | 0.00441905463 | -0.08452916779 | 0.06490343416 | 4.503457745 | 0.2290265586 | 0.7190107689 | 0.7148863862 | 16 to 11 | DIED |
| Mouse | Limb_Muscle | TMS | skeletal muscle cell | 0.01687035724 | 0.04753782612 | 0.03286980828 | 2.957900008 | 0.1820184906 | 0.6537497452 | 0.6542826823 | 16 to 12 | DIED |
| Mouse | Liver | TMS | B cell | 0.01845906587 | 0.2201534011 | 0.08756036112 | 5.534927369 | 0.06575475505 | 0.4292425545 | 0.4312862863 | 12 to 7 | DIED |
| Mouse | Lung | TMS | T cell | 0.00441905463 | 0.05523661379 | 0.01905975824 | 1.951243914 | 0.02304912249 | 0.3137855899 | 0.3186000534 | 16 to 10 | DIED |
| Mouse | Lung | TMS | pulmonary interstitial fibroblast | 0.01299657222 | -0.1045805728 | 0.03363253558 | 1.302131727 | 0.01709303549 | 0.3056530452 | 0.2947495288 | 16 to 10 | DIED |
| Mouse | Marrow | TMS | late pro-B cell | 0.02916220028 | 0.03281722854 | 0.04586454974 | 2.291385103 | 0.4906538608 | 0.8676436225 | 0.8645148375 | 17 to 13 | DIED |
| Mouse | Marrow | TMS | precursor B cell | 0.0289883515 | -0.0453676767 | 0.04047847473 | 3.278885302 | 0.2885731084 | 0.7579256251 | 0.7563269632 | 17 to 13 | DIED |
| Mouse | Spleen | TMS | plasma cell | 0.005348658021 | 0.07943311786 | 0.05525859864 | 4.860243714 | 0.1937394531 | 0.6700793466 | 0.6680581237 | 13 to 10 | DIED |
| Human | lung | Lung_Integrated | club cell | 0.005348658021 | -0.08040728205 | 0.01950646687 | 1.332118428 | 8.475795773e-05 | 0.01177686833 | 0.01090248367 | 119 to 91 | SURVIVED |
| Human | lung | Lung_Integrated | multiciliated epithelial cell | 0.0247881512 | -0.04700185429 | 0.01289610805 | 1.319137427 | 0.0004523849389 | 0.02785398695 | 0.02578594152 | 119 to 91 | SURVIVED |
| Human | brain | MSSM_DLPFC | L6b glutamatergic neuron of the primary motor cortex | 0.01845906587 | 0.01976532187 | 0.004882273583 | 1.081327246 | 7.297520376e-05 | 0.01177686833 | 0.01090248367 | 213 to 209 | SURVIVED |
| Human | brain | MSSM_DLPFC | oligodendrocyte | 0.006417349352 | 0.01726917662 | 0.004377200206 | 1.08453538 | 0.0001092980819 | 0.01177686833 | 0.01090248367 | 213 to 209 | SURVIVED |
| Human | brain | MSSM_DLPFC | oligodendrocyte precursor cell | 0.005348658021 | -0.01189556644 | 0.002988790224 | 1.098541 | 9.545292739e-05 | 0.01177686833 | 0.01090248367 | 213 to 209 | SURVIVED |

Totals: SURVIVED 5; DIED 9; UNTESTABLE 1.

For the 9 DIED rows, median SE inflation is 3.278885302x, with range 1.302131727x to 5.534927369x. These are losses of significance after donor and age-range restriction, not demonstrations of no mature effect. Exact per-row donor counts and SE ratios are in the table above.

## New mature passes

All rows below remain SECONDARY_NOT_PROMOTED. A new pass here means q_431_family or q_refit_family is below 0.05 and the row was not one of the original 15 passes.

| Species | Organ | Cohort | Cell type | Effect | SE | p | q_431_family | q_refit_family | Donors | Present | Cells |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Human | brain | SEAAD_DLPFC | L5 extratelencephalic projecting glutamatergic cortical neuron | -0.0212985384 | 0.005503460284 | 0.001229083311 | 0.05116803895 | 0.04988571378 | 20 | 20 | 1203 |
| Human | lung | Lung_Integrated | conventional dendritic cell | 0.02289666097 | 0.006442846999 | 0.000613187668 | 0.03303548562 | 0.03058273494 | 91 | 39 | 183 |
| Human | lung | Lung_Integrated | multiciliated columnar cell of tracheobronchial tree | -0.03807850571 | 0.01141795564 | 0.001250268516 | 0.05116803895 | 0.04988571378 | 91 | 88 | 21165 |
| Mouse | Limb_Muscle | TMS | macrophage | -0.1429451155 | 0.02649338419 | 0.0004355617472 | 0.02785398695 | 0.02578594152 | 12 | 12 | 1895 |
| Mouse | Limb_Muscle | TMS | mesenchymal stem cell | -0.1004571758 | 0.01683890821 | 0.0002111953918 | 0.01820504278 | 0.01685339227 | 12 | 12 | 9155 |

## Untestable rows and zero-test tissues

The untestable CSV has 96 rows, including cell-level eligibility failures, zero-test tissue summaries, and 4 naked mole-rat rows outside the mature definition. No untestable row enters either mature BH family.

| Species | Organ | Cohort | Mature donors | Candidate cell types | Reason |
| --- | --- | --- | ---: | ---: | --- |
| Mouse | Bladder | TMS | 4 | 4 | zero eligible tests: no cell type met the 6-donor presence rule and produced an estimable age model |
| Mouse | Large_Intestine | TMS | 3 | 5 | zero eligible tests: no cell type met the 6-donor presence rule and produced an estimable age model |
| Mouse | Mammary_Gland | TMS | 4 | 7 | zero eligible tests: no cell type met the 6-donor presence rule and produced an estimable age model |
| Mouse | Skin | TMS | 3 | 5 | zero eligible tests: no cell type met the 6-donor presence rule and produced an estimable age model |
| Mouse | Thymus | TMS | 5 | 6 | zero eligible tests: no cell type met the 6-donor presence rule and produced an estimable age model |
| Mouse | Tongue | TMS | 4 | 3 | zero eligible tests: no cell type met the 6-donor presence rule and produced an estimable age model |
| Mouse | Trachea | TMS | 0 | 12 | zero eligible tests: no cell type met the 6-donor presence rule and produced an estimable age model |

Bladder has 4 mature donors and zero eligible cell types, so the original bladder-cell pass is UNTESTABLE. Large_Intestine, Mammary_Gland, Skin, Thymus, Tongue, and Trachea also have zero eligible mouse tests. Cell-level reasons and mature presence and cell counts are in results/followup-mature-composition-untestable.csv.

## Study-aware human lung at 30-plus

The lung sensitivity uses the D3 donor-ID study proxy, the original adult-eligible Lung_Integrated cell universe, and only donors aged at least 30. Study groups with at least 5 mature donors enter fixed effects; smaller groups are pooled as other. Per-study fits require at least 5 donors and age SD at least 3 years. Meta-analysis uses inverse-variance fixed effects and DerSimonian-Laird random effects. LOSO retains the original rule of at least 20 remaining donors.

| Cell type | Arm | Model or group | Effect | SE | p | Donors | Studies |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| club cell | pooled_model | baseline_age_sex | -0.08045539892 | 0.0198913293 | 0.0001121061906 | 91 | NA |
| club cell | pooled_model | study_fixed_effects | -0.04256887571 | 0.01721855661 | 0.01554699013 | 91 | 9 |
| club cell | per_study_summary | per_study_count | NA | NA | NA | 69 | 7 |
| club cell | meta_combine | fixed_effect_ivw | -0.006115382077 | 0.01111635641 | 0.582233856 | 69 | 7 |
| club cell | meta_combine | random_effects_dl | -0.006115382077 | 0.01111635641 | 0.582233856 | 69 | 7 |
| multiciliated epithelial cell | pooled_model | baseline_age_sex | -0.04704997115 | 0.01335398374 | 0.0006784151757 | 91 | NA |
| multiciliated epithelial cell | pooled_model | study_fixed_effects | -0.0260657306 | 0.01320306292 | 0.05180808035 | 91 | 9 |
| multiciliated epithelial cell | per_study_summary | per_study_count | NA | NA | NA | 69 | 7 |
| multiciliated epithelial cell | meta_combine | fixed_effect_ivw | 0.009842855146 | 0.01148432946 | 0.391407119 | 69 | 7 |
| multiciliated epithelial cell | meta_combine | random_effects_dl | 0.01709793186 | 0.02080205215 | 0.4111139705 | 69 | 7 |

Per-study slopes for club cell:
- 2020-3173-NC: status UNTESTABLE; effect NA; SE NA; p NA; n 4; not estimated: fewer than 5 donors or within-study age SD below 3 years
- AN: status TESTED; effect -0.009693107033; SE 0.03893519276; p 0.82662814; n 5; mature donors only
- D: status UNTESTABLE; effect NA; SE NA; p NA; n 7; not estimated: fewer than 5 donors or within-study age SD below 3 years
- DD: status UNTESTABLE; effect NA; SE NA; p NA; n 1; not estimated: fewer than 5 donors or within-study age SD below 3 years
- Donor: status TESTED; effect 0.008768753707; SE 0.01956210402; p 0.697852364; n 5; mature donors only
- HC: status UNTESTABLE; effect NA; SE NA; p NA; n 1; not estimated: fewer than 5 donors or within-study age SD below 3 years
- LuT: status UNTESTABLE; effect NA; SE NA; p NA; n 3; not estimated: fewer than 5 donors or within-study age SD below 3 years
- NU: status UNTESTABLE; effect NA; SE NA; p NA; n 1; not estimated: fewer than 5 donors or within-study age SD below 3 years
- Patient: status TESTED; effect -0.001685477724; SE 0.09526414846; p 0.9869949265; n 6; mature donors only
- Shalek: status UNTESTABLE; effect NA; SE NA; p NA; n 1; not estimated: fewer than 5 donors or within-study age SD below 3 years
- T: status TESTED; effect -0.03547873096; SE 0.02378445572; p 0.1741250637; n 11; mature donors only
- THD: status TESTED; effect 0.03051885937; SE 0.02737173962; p 0.2909383778; n 13; mature donors only
- VUHD: status TESTED; effect -0.03981843831; SE 0.02790780943; p 0.1728649503; n 19; mature donors only
- lafy: status UNTESTABLE; effect NA; SE NA; p NA; n 4; not estimated: fewer than 5 donors or within-study age SD below 3 years
- muc: status TESTED; effect 0.02255725672; SE 0.05564742261; p 0.6973155912; n 10; mature donors only
LOSO for club cell:
- Drop VUHD: effect -0.0937210648; SE 0.02317142969; p 0.0001345714603; remaining donors 72; status TESTED.
- Drop THD: effect -0.09841032093; SE 0.02060485958; p 8.658385172e-06; remaining donors 78; status TESTED.
- Drop T: effect -0.09021505086; SE 0.02048566488; p 3.387503512e-05; remaining donors 80; status TESTED.
- Drop muc: effect -0.07856693616; SE 0.02264805141; p 0.0008536192238; remaining donors 81; status TESTED.
- Drop D: effect -0.05738691815; SE 0.02148935779; p 0.009151302915; remaining donors 84; status TESTED.
- Drop Patient: effect -0.06948245411; SE 0.0215134656; p 0.001783571119; remaining donors 85; status TESTED.
- Drop AN: effect -0.07984619708; SE 0.01912135105; p 7.291279645e-05; remaining donors 86; status TESTED.
- Drop Donor: effect -0.07949776675; SE 0.01990616508; p 0.0001400021866; remaining donors 86; status TESTED.

Per-study slopes for multiciliated epithelial cell:
- 2020-3173-NC: status UNTESTABLE; effect NA; SE NA; p NA; n 4; not estimated: fewer than 5 donors or within-study age SD below 3 years
- AN: status TESTED; effect 0.05746046971; SE 0.05997362684; p 0.4391197578; n 5; mature donors only
- D: status UNTESTABLE; effect NA; SE NA; p NA; n 7; not estimated: fewer than 5 donors or within-study age SD below 3 years
- DD: status UNTESTABLE; effect NA; SE NA; p NA; n 1; not estimated: fewer than 5 donors or within-study age SD below 3 years
- Donor: status TESTED; effect 0.0897584089; SE 0.02898054161; p 0.09034266653; n 5; mature donors only
- HC: status UNTESTABLE; effect NA; SE NA; p NA; n 1; not estimated: fewer than 5 donors or within-study age SD below 3 years
- LuT: status UNTESTABLE; effect NA; SE NA; p NA; n 3; not estimated: fewer than 5 donors or within-study age SD below 3 years
- NU: status UNTESTABLE; effect NA; SE NA; p NA; n 1; not estimated: fewer than 5 donors or within-study age SD below 3 years
- Patient: status TESTED; effect -0.001685477724; SE 0.06651202634; p 0.9813743961; n 6; mature donors only
- Shalek: status UNTESTABLE; effect NA; SE NA; p NA; n 1; not estimated: fewer than 5 donors or within-study age SD below 3 years
- T: status TESTED; effect -0.02078821039; SE 0.02006432558; p 0.3304655726; n 11; mature donors only
- THD: status TESTED; effect 0.01061523784; SE 0.0363098024; p 0.7759933676; n 13; mature donors only
- VUHD: status TESTED; effect -0.04446216385; SE 0.02648111959; p 0.1125685889; n 19; mature donors only
- lafy: status UNTESTABLE; effect NA; SE NA; p NA; n 4; not estimated: fewer than 5 donors or within-study age SD below 3 years
- muc: status TESTED; effect 0.04853286591; SE 0.02864350154; p 0.1340198378; n 10; mature donors only
LOSO for multiciliated epithelial cell:
- Drop VUHD: effect -0.04969443821; SE 0.01530893577; p 0.001806726932; remaining donors 72; status TESTED.
- Drop THD: effect -0.05642727386; SE 0.013949469; p 0.0001257320866; remaining donors 78; status TESTED.
- Drop T: effect -0.05394634813; SE 0.0130698338; p 9.206968193e-05; remaining donors 80; status TESTED.
- Drop muc: effect -0.04480644132; SE 0.01510575789; p 0.004000348158; remaining donors 81; status TESTED.
- Drop D: effect -0.04065806252; SE 0.01451056056; p 0.006351797256; remaining donors 84; status TESTED.
- Drop Patient: effect -0.03574840778; SE 0.01415786939; p 0.01349716923; remaining donors 85; status TESTED.
- Drop AN: effect -0.04936117873; SE 0.01299070871; p 0.0002751718766; remaining donors 86; status TESTED.
- Drop Donor: effect -0.0478516695; SE 0.01375207823; p 0.0008027567557; remaining donors 86; status TESTED.

## Secondary-hit checks

No secondary hit is promoted to a headline pass. Total-cell models add log total recovered eligible cells. Donor influence repeats the primary model leaving out one donor at a time. Study-aware status for Lung_Integrated hits is determined from the preceding lung output. Mature propeller ran on beast with R 4.3.3, speckle 1.2.0, and limma 3.58.1, using the same design as the original check (logit proportions; age plus log total cells plus sex when variable). The section below gives the results.

| Species | Organ | Cohort | Cell type | Total-cell p | Total-cell refit q | LODO fits | LODO maximum p | LODO same sign | Propeller |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: | --- | --- |
| Human | brain | SEAAD_DLPFC | L5 extratelencephalic projecting glutamatergic cortical neuron | 0.002107391674 | 0.06468071368 | 20 | 0.02430619094 | True | PASS |
| Human | lung | Lung_Integrated | conventional dendritic cell | 0.008066316875 | 0.1418663284 | 91 | 0.001500936388 | True | PASS |
| Human | lung | Lung_Integrated | multiciliated columnar cell of tracheobronchial tree | 0.01582111886 | 0.2054606442 | 91 | 0.003852883816 | True | FAIL |
| Mouse | Limb_Muscle | TMS | macrophage | 0.0003256323806 | 0.01913779075 | 12 | 0.003153531917 | True | PASS |
| Mouse | Limb_Muscle | TMS | mesenchymal stem cell | 0.0003357507149 | 0.01913779075 | 12 | 0.001800727759 | True | FAIL |

The three previously reported secondary findings are retained: mouse macrophage and mesenchymal stem cell declines, and human Lung_Integrated conventional dendritic-cell increase. All three are present in the new-pass table: True. They remain secondary and their check values are listed above.

Study-aware check for conventional dendritic cell: study fixed effect 0.01025419425, p 0.203111777; fixed-effect meta 0.007648050441, p 0.02540409979; random-effects meta 0.01182273511, p 0.1187782824; LOSO effects 0.01548513331 to 0.02605363413, p range 0.000410704574 to 0.0214228166 across 8 arms.
Study-aware check for multiciliated columnar cell of tracheobronchial tree: study fixed effect -0.02211660068, p 0.09511581062; fixed-effect meta -0.008133912336, p 0.5472720408; random-effects meta -0.008133912336, p 0.5472720408; LOSO effects -0.04564557616 to -0.02739201724, p range 0.0002316093175 to 0.03207667576 across 8 arms.

## Mature propeller check

R 4.3.3, speckle 1.2.0, and limma 3.58.1 were installed on beast, and the unchanged original runner scripts/run_reviewer_phase1_propeller.R was run on mature-only inputs (mice 18 months and older, human donors 30 years and older) prepared by scripts/prepare_mature_propeller.py. The design matches the original check: logit proportions with a 0.5 pseudocount, age plus log total cells plus sex when variable, limma empirical Bayes on the age coefficient, BH within cohort. 23 cohorts ran; 6 mouse tissues were skipped for fewer than 5 mature donors (Bladder, Large_Intestine, Mammary_Gland, Skin, Tongue, Trachea). 9 HRCA retina suites ran (88 rows). Independent BH recomputation matches the R output to a maximum absolute difference of 2.665e-15 across all 577 rows. Full output: results/followup-mature-propeller.csv. Key rows: results/followup-mature-propeller-key-rows.csv.

AGREE means propeller q below 0.05 with the same sign as the CLR effect. NOT_CONFIRMED means propeller q at or above 0.05 with the same sign. No row among the original 11 showed a sign conflict. Four HRCA rows have opposite signs between CLR and propeller, all on results that are null on both methods. Mouse effects are per month and human effects are per year, matching the CLR convention.

| Row | CLR effect | CLR q_431 | Propeller effect (logit) | Propeller p | Propeller q | Verdict |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| HEADLINE lung club | -0.0804 | 0.01178 | -0.0787 | 8.072e-05 | 0.002887 | AGREE |
| HEADLINE lung multiciliated | -0.0470 | 0.02785 | -0.0459 | 0.0001155 | 0.002887 | AGREE |
| HEADLINE MSSM L6b | 0.0198 | 0.01178 | 0.0176 | 0.001093 | 0.004792 | AGREE |
| HEADLINE MSSM oligodendrocyte | 0.0173 | 0.01178 | 0.0215 | 5.376e-06 | 6.72e-05 | AGREE |
| HEADLINE MSSM OPC | -0.0119 | 0.01178 | -0.0157 | 2.888e-08 | 7.219e-07 | AGREE |
| SECONDARY mouse limb muscle macrophage | -0.1429 | 0.02785 | -0.0935 | 0.00539 | 0.01664 | AGREE |
| SECONDARY mouse limb muscle mesenchymal stem cell | -0.1005 | 0.01821 | -0.0694 | 0.02559 | 0.05758 | NOT_CONFIRMED |
| SECONDARY lung conventional dendritic | 0.0229 | 0.03304 | 0.0240 | 0.003676 | 0.03063 | AGREE |
| SECONDARY SEA-AD L5 | -0.0213 | 0.05117 | -0.0261 | 0.00133 | 0.02393 | AGREE |
| SECONDARY lung multiciliated columnar | -0.0381 | 0.05117 | -0.0242 | 0.0344 | 0.1014 | NOT_CONFIRMED |

Thomson PBMC CD8 naive: CLR effect minus 1.318945, q 0.005778; propeller effect minus 0.0520, p 5.832e-05, q 0.0007582; AGREE. Propeller has no batch term, so this does not overturn the batch failure (batch q 0.2826).

HRCA retina amacrine is now run: 9 suites, 88 rows, all donors 30 years and older, inputs built by scripts/prepare_mature_propeller_hrca.py. Propeller confirms none of the 9 suites: every amacrine q is 0.43 or higher. The three CLR-significant suites (Sn_Adult_Unenriched, Sn_Adult_ExactAge, Sn_AllAges_Unenriched) are NOT_CONFIRMED. The other six suites are NULL on both methods.

For HRCA rows, the CLR effect is per decade (batch B used age per decade) and the CLR q is the batch-B family-A q within the suite, not a 431-family q. NULL means neither CLR nor propeller passes q below 0.05.

| Row | CLR effect (per decade) | CLR q (family A) | Propeller effect (logit, per year) | Propeller p | Propeller q | Verdict |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| HRCA Sn_Adult_MacularOnly amacrine | 0.1099 | 0.30928 | 0.0010 | 0.873 | 0.873 | NULL |
| HRCA Sn_Adult_PeripheralOnly amacrine | -0.0108 | 0.86493 | 0.0018 | 0.7525 | 0.9345 | NULL |
| HRCA Sn_Adult_Unenriched amacrine | 0.2479 | 0.0014338 | 0.0024 | 0.6794 | 0.8092 | NOT_CONFIRMED |
| HRCA Sn_Adult_ExactAge amacrine | 0.2930 | 0.00050533 | 0.0043 | 0.5115 | 0.5684 | NOT_CONFIRMED |
| HRCA Sn_ChenAncestry_Adult amacrine | 0.1038 | 0.16948 | -0.0016 | 0.702 | 0.78 | NULL |
| HRCA Sn_AllAges_Unenriched amacrine | 0.2479 | 0.0014338 | 0.0024 | 0.6794 | 0.8092 | NOT_CONFIRMED |
| HRCA Sc15_StudyAdj amacrine | 0.1876 | 0.91336 | -0.0044 | 0.8561 | 0.9622 | NULL |
| HRCA Sc9_External_StudyAdj amacrine | 0.4002 | 0.33447 | 0.0222 | 0.6224 | 0.9404 | NULL |
| HRCA Sc5_External_NoDisease amacrine | -0.0107 | 0.96409 | 0.0126 | 0.284 | 0.426 | NULL |

Why propeller disagrees: the propeller design has no study term, while the batch-B CLR model adjusts for study, sex, peripheral fraction, and log library size. A refit in Sn_Adult_Unenriched (scripts/diag_hrca_propeller_gap.py) shows the amacrine CLR signal depends on that adjustment: unadjusted slope plus 0.082 per decade (p 0.36); with study dummies plus sex, plus 0.218 (p 0.0017); full batch-B model, plus 0.248 (p 0.00019). The logit refit with propeller-like covariates reproduces the propeller estimate (plus 0.00245 per year), so the gap is the covariate set, not the transform or the inputs.

## Runtime and verification limits

Pipeline model: kimi-code/k3.1-pro.

| Pipeline | Start UTC | End UTC | Duration seconds | Status |
| --- | --- | --- | ---: | --- |
| complete_mature_composition | 2026-11-21T20:04:12.503379+00:00 | 2026-11-21T20:04:15.456583+00:00 | 2.953204 | completed |
| mature_lung_study_sensitivity | 2026-11-21T20:04:24.785610+00:00 | 2026-11-21T20:04:25.821547+00:00 | 1.035937 | completed |

One self-review independently recomputed both BH families, checked all output age minima, accounted for every original 431-family row as tested or untestable, checked all 44 lung cell types across 28 study-arm rows, and verified the locked 624-file hash. All 399 tested rows belong to the original family; 32 original rows are untestable. An in-memory call to the unchanged base mouse pipeline on the current h5ad files reproduces the new mouse effects to numerical precision. Seventeen rows in the older stored mature mouse CSV differ slightly, with a maximum absolute effect difference of 0.000697691 and at most 7 recovered cells, so that stored CSV reflects an earlier state of the per-tissue inputs; the original-survivor verdicts and three prior secondary hits are unchanged.

Mature propeller estimates are verified: independent BH recomputation matches the R output to a maximum absolute difference of 2.665e-15 across all 577 rows including the 88 HRCA rows. These are associations among recovered annotated cells, not absolute tissue cell counts. No manuscript file was read or changed for this item, and no non-mammalian extended-family row was refit.
