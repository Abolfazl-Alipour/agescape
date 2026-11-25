# All-species mature-to-old cell-fraction scan

Analysis ordered by Reza on 2026-11-24. Pipeline model: `kimi-code/k3.1-pro`.

## Scope, model, and eligibility

This is a mature-to-old fraction scan only. The inferential unit is the human donor, mouse or naked-mole-rat animal, zebrafish fish, or fly/worm pool; cells are never independent replicates. Fly head/body samples and Gao worm libraries are pools, as verified in `scripts/run_rev3_composition_fly_worm_zebrafish.py` and `scripts/run_gao_worm_aging.py`.

For each organ/cohort stratum, donor-by-cell-type counts receive a 0.5 pseudocount and a centered log-ratio (CLR) transform. Ordinary least squares fits CLR fraction on continuous age plus sex; sex is dropped when the design is rank-deficient. Slopes are per month for mouse and zebrafish, per year for human, and per day for fly and worm. The primary model does not adjust for total cells or sequencing platform.

The pipeline rule, not the tentative 3-unit/50-cell rule in the order, is used: a cell type must have count > 0 in at least 6 mature inferential units and the age design must be full rank with positive residual degrees of freedom. There is no additional total-cell >= 50 exclusion in the primary pipeline; that is a sensitivity arm in the small-species full-range runner. Untestable rows never enter a correction family.

## Mature windows and species decisions

Cutoffs were fixed from the repository record and source biology before fitting:

- **Mouse: 18 months and older.** Reuses the established mature composition window exactly.
- **Human: 30 years and older.** Reuses the established mature composition window exactly.
- **Fruit fly: 30 days and older.** AFCA's primary analyzed ages are 30, 50, and 70 days; day 30 is the repository's young-adult reference. Day-5 FCA cells remain excluded because study and age are inseparable. Thus the mature fly fit is the same age window as the primary AFCA fit.
- **C. elegans: 6 days and older.** Day 6 is the repository's balanced mature reference and includes both sequencing platforms. Day-1 pools are excluded, leaving days 6, 12, and 14.
- **Zebrafish: 12 months and older, but excluded.** The [source atlas](https://pmc.ncbi.nlm.nih.gov/articles/PMC11607696/) says males become fertile at 3 months, describes 5-month testes as approaching peak spermatogenesis, and identifies 12-month testes as closest to peak fertility. Twelve months is therefore the fully mature reference. Only four fish (two at 12 months, one at 20, one at 22) retain comparable fine labels, below the six-fish floor. The primary fine-label arm is not run; the separate broad-label reconciliation arm is not substituted to force a zebrafish result.
- **Naked mole-rat: excluded.** Direct inspection of the GSE152886 count-file header confirms the established samples. After excluding 2V, there are three animals at 1.1 years and two at 11.3 years, with no intermediate ages. There is no defensible mature-to-old window satisfying the six-animal rule.

## Multiple-testing families declared before correction

**Primary pooled mature family:** one Benjamini-Hochberg correction over all 577 tests actually run across all included species. `q_pooled_mature_bh < 0.05` is the only pass bar.

Transparency-only families: `q_species_bh` applies BH separately within each species, and `q_bh` applies BH separately within each organ/cohort stratum. Neither determines a pooled pass. Family sizes are recorded on every result row.

| Species | Mature tests | Pooled passes | Species-family passes |
| --- | ---: | ---: | ---: |
| Mouse | 114 | 2 | 2 |
| Human | 285 | 11 | 8 |
| Fruit fly | 162 | 25 | 36 |
| C. elegans | 16 | 1 | 1 |

## Donor, animal, and pool coverage

| Species | Organ | Cohort | Unit | Mature ages | Units | Tested cell types |
| --- | --- | --- | --- | --- | ---: | ---: |
| Mouse | Bladder | TMS | animal | 18;24 | 4 | 0 |
| Mouse | Fat | TMS | animal | 18;21;30 | 6 | 5 |
| Mouse | Heart_and_Aorta | TMS | animal | 18;21;24;30 | 9 | 8 |
| Mouse | Kidney | TMS | animal | 18;21;24;30 | 11 | 19 |
| Mouse | Large_Intestine | TMS | animal | 30 | 3 | 0 |
| Mouse | Limb_Muscle | TMS | animal | 18;21;24;30 | 12 | 9 |
| Mouse | Liver | TMS | animal | 18;21;24;30 | 7 | 6 |
| Mouse | Lung | TMS | animal | 18;21;30 | 10 | 29 |
| Mouse | Mammary_Gland | TMS | animal | 18;21 | 4 | 0 |
| Mouse | Marrow | TMS | animal | 18;21;24;30 | 13 | 18 |
| Mouse | Pancreas | TMS | animal | 18;21;24;30 | 8 | 8 |
| Mouse | Skin | TMS | animal | 18;21 | 3 | 0 |
| Mouse | Spleen | TMS | animal | 18;21;24;30 | 10 | 12 |
| Mouse | Thymus | TMS | animal | 18;21;24 | 5 | 0 |
| Mouse | Tongue | TMS | animal | 18;24 | 4 | 0 |
| Mouse | Trachea | TMS | animal |  | 0 | 0 |
| Human | brain | HBCC_DLPFC | donor | 30;31;32;35;36;37;38;40;41;42;43;44;45;46;48;50;51;52;53;54;55;56;57;58;59;60;61;63;65;67;71;78;83;85 | 64 | 24 |
| Human | lung | ImmuneAging_BAL | donor | 33;48;56;64;67;68;73;75 | 8 | 31 |
| Human | kidney | Kidney_LivingDonor | donor | 33;34;35;37;38;41;42;43;44;49;50;55;57;58;60;64;65 | 19 | 17 |
| Human | kidney | Kidney_NormalAdult | donor | 43;46;50;53;57;58;59;60;67;68;69;71;72;74 | 18 | 26 |
| Human | kidney | Kidney_SCR_v15 | donor | 30;31;35;39;40;41;42;47;48;49;52;53;56;57;66 | 18 | 23 |
| Human | kidney | Kidney_SN_v15 | donor | 45;54;57;59;63;66;69;76 | 8 | 19 |
| Human | lung | Lung_CellRef | donor | 35;39;42;45;50;51;52;57;59;68;69;70;76 | 13 | 29 |
| Human | lung | Lung_Integrated | donor | 30;31;33;34;35;36;37;38;40;41;42;44;45;47;49;52;53;54;55;56;57;58;59;60;61;62;63;64;65;66;67;68;69;71;72;73;75;76;81 | 91 | 50 |
| Human | brain | MSSM_DLPFC | donor | 31;35;36;38;39;40;41;42;43;45;46;47;48;49;50;51;52;53;54;55;56;57;58;59;60;61;62;63;64;65;66;67;68;69;70;71;72;73;74;75;76;77;78;79;80;81;82;83;84;85;86;87;88;89 | 209 | 25 |
| Human | brain | RADC_DLPFC | donor | 64;70;75;76;77;78;79;80;81;82;84;85;86;89 | 18 | 23 |
| Human | brain | SEAAD_DLPFC | donor | 42;50;72;75;78;80;81;82;83;84;85;86;87;88;89 | 20 | 18 |
| Fruit fly | head | AFCA | pool | 30;50;70 | 36 | 90 |
| Fruit fly | body | AFCA | pool | 30;50;70 | 36 | 72 |
| C. elegans | Whole animal | Gao_GSE229022 | pool | 6;12;14 | 9 | 16 |
| Zebrafish | Testis | GSE275361 | fish | 12;20;22 | 4 | 0 |
| Naked mole-rat | Skin | GSE152886 | animal | not defined; all available ages: 1.1;11.3 | 0 | 0 |

## Primary pooled results: 39 passes

Slopes are in CLR units per the stated age unit. Numeric values below are read back from `results/mature-all-species-composition.csv`.

| Species | Organ | Cohort | Cell type | Slope | Age unit | p | Species q | Pooled mature q |
| --- | --- | --- | --- | ---: | --- | ---: | ---: | ---: |
| Mouse | Limb_Muscle | TMS | macrophage | -0.14294511547 | month | 0.000435561747219 | 0.0248270195915 | 0.0124298147488 |
| Mouse | Limb_Muscle | TMS | mesenchymal stem cell | -0.100457175849 | month | 0.000211195391831 | 0.0240762746687 | 0.0081239827391 |
| Human | brain | MSSM_DLPFC | GABAergic neuron | -0.0117586248482 | year | 0.00158195848906 | 0.0500953521536 | 0.0294448402641 |
| Human | brain | MSSM_DLPFC | L6b glutamatergic neuron of the primary motor cortex | 0.0197653218673 | year | 7.29752037635e-05 | 0.00778748833623 | 0.00421066925715 |
| Human | brain | MSSM_DLPFC | oligodendrocyte | 0.0172691766153 | year | 0.000109298081912 | 0.00778748833623 | 0.00450464237594 |
| Human | brain | MSSM_DLPFC | oligodendrocyte precursor cell | -0.0118955664369 | year | 9.54529273902e-05 | 0.00778748833623 | 0.00423664146955 |
| Human | brain | MSSM_DLPFC | sst GABAergic cortical interneuron | -0.0126458758056 | year | 0.00316176335803 | 0.0819184142763 | 0.0467778835278 |
| Human | brain | SEAAD_DLPFC | L5 extratelencephalic projecting glutamatergic cortical neuron | -0.0212985383973 | year | 0.00122908331053 | 0.0445408158765 | 0.0248760321943 |
| Human | kidney | Kidney_LivingDonor | kidney cortex collecting duct epithelial cell | 0.0803463500705 | year | 0.00247762416665 | 0.0706122887495 | 0.0408454041188 |
| Human | lung | Lung_Integrated | club cell | -0.0804072820487 | year | 8.47579577302e-05 | 0.00778748833623 | 0.00423664146955 |
| Human | lung | Lung_Integrated | conventional dendritic cell | 0.0228966609694 | year | 0.000613187668043 | 0.0291264142321 | 0.0144885031039 |
| Human | lung | Lung_Integrated | multiciliated columnar cell of tracheobronchial tree | -0.0380785057062 | year | 0.00125026851583 | 0.0445408158765 | 0.0248760321943 |
| Human | lung | Lung_Integrated | multiciliated epithelial cell | -0.0470018542856 | year | 0.000452384938865 | 0.0257859415153 | 0.0124298147488 |
| Fruit fly | body | AFCA | 16-cell germline cyst in germarium region 2a and 2b | 0.0242663015547 | day | 0.00023822575474 | 0.00343657025065 | 0.00859101628032 |
| Fruit fly | body | AFCA | adult Malpighian tubule principal cell of lower ureter | -0.025515218896 | day | 7.08394326584e-05 | 0.00127510978785 | 0.00421066925715 |
| Fruit fly | body | AFCA | adult differentiating enterocyte | 0.0248601159417 | day | 0.000542252638366 | 0.00585632849435 | 0.0142218078335 |
| Fruit fly | body | AFCA | adult fat body_body | -0.0184790275693 | day | 1.70233992779e-09 | 2.75779068302e-07 | 9.82250138333e-07 |
| Fruit fly | body | AFCA | adult ventral nervous system | -0.00811540242914 | day | 0.00184621809006 | 0.0135948786632 | 0.0322808435747 |
| Fruit fly | body | AFCA | ejaculatory bulb | -0.0160598731814 | day | 9.01900715611e-05 | 0.00146107915929 | 0.00423664146955 |
| Fruit fly | body | AFCA | enteroendocrine cell | 0.0378402173986 | day | 0.000954467075408 | 0.00813808769558 | 0.0203973149078 |
| Fruit fly | body | AFCA | escort cell | 0.0121020924821 | day | 0.00152335115194 | 0.0123391443307 | 0.0292991204889 |
| Fruit fly | body | AFCA | follicle stem cell and prefollicle cell | 0.0141891774047 | day | 1.67877084784e-05 | 0.000453268128916 | 0.00161441796534 |
| Fruit fly | body | AFCA | hemocyte_body | -0.0109917080571 | day | 0.000254560759307 | 0.00343657025065 | 0.00864009165413 |
| Fruit fly | body | AFCA | male accessory gland main cell | -0.0133681854342 | day | 0.0019821449996 | 0.0139611952146 | 0.0336381666108 |
| Fruit fly | body | AFCA | posterior midgut* | -0.0340166526474 | day | 0.000789631799331 | 0.00710668619398 | 0.0175237518544 |
| Fruit fly | body | AFCA | prefollicle cell/stalk follicle cell | 0.0179125728938 | day | 3.3542152364e-06 | 0.000135845717074 | 0.000483845547851 |
| Fruit fly | body | AFCA | stalk follicle cell | 0.0219356651945 | day | 1.0573746542e-05 | 0.000342589387961 | 0.00122021035095 |
| Fruit fly | head | AFCA | T neuron T3 | 0.0113145718363 | day | 0.000420822273437 | 0.00486951487834 | 0.0124298147488 |
| Fruit fly | head | AFCA | centrifugal neuron C2 | 0.0101613907527 | day | 0.00275285653143 | 0.0178385103236 | 0.0417999531219 |
| Fruit fly | head | AFCA | distal medullary amacrine neuron Dm8 | 0.0104355393994 | day | 0.000403392716115 | 0.00486951487834 | 0.0124298147488 |
| Fruit fly | head | AFCA | epithelial cell_head | 0.0136624172515 | day | 4.4155378701e-09 | 3.57658567478e-07 | 1.27388267552e-06 |
| Fruit fly | head | AFCA | neuron 1 | 0.00879597663687 | day | 0.00263861506418 | 0.0178106516832 | 0.0416794598508 |
| Fruit fly | head | AFCA | ocellus retinula cell | -0.02144093805 | day | 0.00177264126818 | 0.0135948786632 | 0.0319629378669 |
| Fruit fly | head | AFCA | olfactory receptor neuron | -0.00628804359814 | day | 0.000608131858241 | 0.00598210190344 | 0.0144885031039 |
| Fruit fly | head | AFCA | photoreceptor cell R7 | -0.0123958760793 | day | 4.83495442102e-05 | 0.00111894659458 | 0.00398538385847 |
| Fruit fly | head | AFCA | pigment cell | 0.0119886830397 | day | 5.83724454103e-05 | 0.00118204201956 | 0.00421011262522 |
| Fruit fly | head | AFCA | skeletal muscle of head | 0.011268533147 | day | 1.20533015099e-08 | 6.50878281536e-07 | 2.31825165707e-06 |
| Fruit fly | head | AFCA | transmedullary neuron Tm29 | -0.0191562922341 | day | 0.000627751434312 | 0.00598210190344 | 0.0144885031039 |
| C. elegans | Whole animal | Gao_GSE229022 | ARSC | 0.227491777223 | day | 0.00267268633359 | 0.0427629813374 | 0.0416794598508 |

## Four prespecified headline rows

4 of 4 headline rows survive the pooled mature correction.

| Species | Cell type | Slope | p | Pooled mature q | Survives |
| --- | --- | ---: | ---: | ---: | --- |
| Human | L6b glutamatergic neuron of the primary motor cortex | 0.0197653218673 | 7.29752037635e-05 | 0.00421066925715 | True |
| Human | oligodendrocyte | 0.0172691766153 | 0.000109298081912 | 0.00450464237594 | True |
| Human | oligodendrocyte precursor cell | -0.0118955664369 | 9.54529273902e-05 | 0.00423664146955 | True |
| Mouse | macrophage | -0.14294511547 | 0.000435561747219 | 0.0124298147488 | True |

## Comparison with the full-range 624-test family

The locked full-range CSV has 45 passes. Of these, 32 pass the new pooled mature bar, 12 are tested but do not pass, and 1 are untestable in the declared mature window. The new pooled scan has 7 passes that were not full-range 624 passes.

| Species | Organ | Cell type | Full-range slope | Full-range p | Full-range q | Mature status | Mature pooled q |
| --- | --- | --- | ---: | ---: | ---: | --- | ---: |
| Mouse | Bladder | bladder cell | -0.0613267206959 | 0.00130591282714 | 0.0226358223371 | UNTESTABLE_MATURE | NA |
| Mouse | Kidney | kidney capillary endothelial cell | -0.110325132582 | 1.98099836692e-05 | 0.0015451787262 | NOT_RETAINED | 0.2532807925 |
| Mouse | Kidney | kidney mesangial cell | -0.0897633438043 | 3.07590809491e-05 | 0.00174487877384 | NOT_RETAINED | 0.570376900626 |
| Mouse | Limb_Muscle | skeletal muscle cell | 0.0532438089221 | 0.00035228124161 | 0.00955754325064 | NOT_RETAINED | 0.500117471668 |
| Mouse | Liver | B cell | 0.0847133157292 | 0.000459346977768 | 0.0108879460436 | NOT_RETAINED | 0.289595584675 |
| Mouse | Lung | T cell | 0.0618724647173 | 2.59978529012e-05 | 0.00174487877384 | NOT_RETAINED | 0.177469401013 |
| Mouse | Lung | pulmonary interstitial fibroblast | -0.129247913991 | 0.000241235679274 | 0.00716814589842 | NOT_RETAINED | 0.151733561188 |
| Mouse | Lung | vein endothelial cell | 0.0564900823989 | 0.00271978389096 | 0.0390405108093 | NOT_RETAINED | 0.915042094992 |
| Mouse | Marrow | late pro-B cell | 0.083432285188 | 0.000947264046254 | 0.0170167844301 | NOT_RETAINED | 0.804763819344 |
| Mouse | Marrow | precursor B cell | -0.0519663030205 | 0.000874358629967 | 0.0165333268212 | NOT_RETAINED | 0.623620537683 |
| Mouse | Spleen | plasma cell | 0.0758143269878 | 5.58106714048e-05 | 0.00260174328115 | NOT_RETAINED | 0.516743537949 |
| Human | kidney | kidney cortex collecting duct epithelial cell | 0.0803463500705 | 0.00247762416665 | 0.0377082312193 | RETAINED | 0.0408454041188 |
| Human | lung | club cell | -0.060165453717 | 7.4459276399e-05 | 0.00273309343959 | RETAINED | 0.00423664146955 |
| Human | lung | multiciliated epithelial cell | -0.0340921109919 | 0.000690157342 | 0.0138921994003 | RETAINED | 0.0124298147488 |
| Human | brain | L6b glutamatergic neuron of the primary motor cortex | 0.0160398644865 | 0.000471113049962 | 0.0108879460436 | RETAINED | 0.00421066925715 |
| Human | brain | oligodendrocyte | 0.015966479002 | 0.000104226091563 | 0.00342300427027 | RETAINED | 0.00450464237594 |
| Human | brain | oligodendrocyte precursor cell | -0.0110908992726 | 6.48554388903e-05 | 0.00269798625784 | RETAINED | 0.00423664146955 |
| Fruit fly | head | T neuron T3 | 0.0113145718363 | 0.000420822273437 | 0.010503723945 | RETAINED | 0.0124298147488 |
| Fruit fly | head | centrifugal neuron C2 | 0.0101613907527 | 0.00275285653143 | 0.0390405108093 | RETAINED | 0.0417999531219 |
| Fruit fly | head | distal medullary amacrine neuron Dm8 | 0.0104355393994 | 0.000403392716115 | 0.010488210619 | RETAINED | 0.0124298147488 |
| Fruit fly | head | dopaminergic PAM neuron | -0.0119093417944 | 0.00358223067016 | 0.0496735986262 | NOT_RETAINED | 0.051673677417 |
| Fruit fly | head | epithelial cell_head | 0.0136624172515 | 4.4155378701e-09 | 1.37764781547e-06 | RETAINED | 1.27388267552e-06 |
| Fruit fly | head | neuron 1 | 0.00879597663687 | 0.00263861506418 | 0.0390405108093 | RETAINED | 0.0416794598508 |
| Fruit fly | head | ocellus retinula cell | -0.02144093805 | 0.00177264126818 | 0.0291086355617 | RETAINED | 0.0319629378669 |
| Fruit fly | head | olfactory receptor neuron | -0.00628804359814 | 0.000608131858241 | 0.0130572298337 | RETAINED | 0.0144885031039 |
| Fruit fly | head | photoreceptor cell R7 | -0.0123958760793 | 4.83495442102e-05 | 0.00251417629893 | RETAINED | 0.00398538385847 |
| Fruit fly | head | pigment cell | 0.0119886830397 | 5.83724454103e-05 | 0.00260174328115 | RETAINED | 0.00421011262522 |
| Fruit fly | head | skeletal muscle of head | 0.011268533147 | 1.20533015099e-08 | 2.50708671406e-06 | RETAINED | 2.31825165707e-06 |
| Fruit fly | head | transmedullary neuron Tm29 | -0.0191562922341 | 0.000627751434312 | 0.0130572298337 | RETAINED | 0.0144885031039 |
| Fruit fly | body | 16-cell germline cyst in germarium region 2a and 2b | 0.0242663015547 | 0.00023822575474 | 0.00716814589842 | RETAINED | 0.00859101628032 |
| Fruit fly | body | adult Malpighian tubule principal cell of lower ureter | -0.025515218896 | 7.08394326584e-05 | 0.00273309343959 | RETAINED | 0.00421066925715 |
| Fruit fly | body | adult differentiating enterocyte | 0.0248601159417 | 0.000542252638366 | 0.0120844873693 | RETAINED | 0.0142218078335 |
| Fruit fly | body | adult fat body_body | -0.0184790275693 | 1.70233992779e-09 | 1.06226011494e-06 | RETAINED | 9.82250138333e-07 |
| Fruit fly | body | adult ventral nervous system | -0.00811540242914 | 0.00184621809006 | 0.029539489441 | RETAINED | 0.0322808435747 |
| Fruit fly | body | ejaculatory bulb | -0.0160598731814 | 9.01900715611e-05 | 0.00312658914745 | RETAINED | 0.00423664146955 |
| Fruit fly | body | enteroendocrine cell | 0.0378402173986 | 0.000954467075408 | 0.0170167844301 | RETAINED | 0.0203973149078 |
| Fruit fly | body | escort cell | 0.0121020924821 | 0.00152335115194 | 0.0256911113192 | RETAINED | 0.0292991204889 |
| Fruit fly | body | follicle stem cell and prefollicle cell | 0.0141891774047 | 1.67877084784e-05 | 0.00149650429864 | RETAINED | 0.00161441796534 |
| Fruit fly | body | hemocyte_body | -0.0109917080571 | 0.000254560759307 | 0.00722026880944 | RETAINED | 0.00864009165413 |
| Fruit fly | body | male accessory gland main cell | -0.0133681854342 | 0.0019821449996 | 0.0309214619937 | RETAINED | 0.0336381666108 |
| Fruit fly | body | posterior midgut* | -0.0340166526474 | 0.000789631799331 | 0.015397820087 | RETAINED | 0.0175237518544 |
| Fruit fly | body | prefollicle cell/stalk follicle cell | 0.0179125728938 | 3.3542152364e-06 | 0.000523257576878 | RETAINED | 0.000483845547851 |
| Fruit fly | body | stalk follicle cell | 0.0219356651945 | 1.0573746542e-05 | 0.00118877946187 | RETAINED | 0.00122021035095 |
| C. elegans | Whole animal | ARSC | 0.216003674045 | 1.14305717487e-05 | 0.00118877946187 | RETAINED | 0.0416794598508 |
| C. elegans | Whole animal | Spermatheca | -0.136623514036 | 2.94022267946e-05 | 0.00174487877384 | NOT_RETAINED | 0.171814134853 |

### Pooled mature passes not among the 45 full-range passes (7)

| Species | Organ | Cohort | Cell type | Mature slope | Mature p | Pooled mature q |
| --- | --- | --- | --- | ---: | ---: | ---: |
| Mouse | Limb_Muscle | TMS | macrophage | -0.14294511547 | 0.000435561747219 | 0.0124298147488 |
| Mouse | Limb_Muscle | TMS | mesenchymal stem cell | -0.100457175849 | 0.000211195391831 | 0.0081239827391 |
| Human | brain | MSSM_DLPFC | GABAergic neuron | -0.0117586248482 | 0.00158195848906 | 0.0294448402641 |
| Human | brain | MSSM_DLPFC | sst GABAergic cortical interneuron | -0.0126458758056 | 0.00316176335803 | 0.0467778835278 |
| Human | brain | SEAAD_DLPFC | L5 extratelencephalic projecting glutamatergic cortical neuron | -0.0212985383973 | 0.00122908331053 | 0.0248760321943 |
| Human | lung | Lung_Integrated | conventional dendritic cell | 0.0228966609694 | 0.000613187668043 | 0.0144885031039 |
| Human | lung | Lung_Integrated | multiciliated columnar cell of tracheobronchial tree | -0.0380785057062 | 0.00125026851583 | 0.0248760321943 |

## Comparison with the mammalian 399-row mature scan

The previous mammalian-only family had 10 passes. 10 of those 10 pass the new all-species pooled bar; 0 do not. All 399 mammalian tests are retained with the same fitted effects and p-values; only the correction family changes.

| Species | Organ | Cell type | Previous 399-row q | New pooled q | Retained |
| --- | --- | --- | ---: | ---: | --- |
| Mouse | Limb_Muscle | macrophage | 0.0257859415153 | 0.0124298147488 | True |
| Mouse | Limb_Muscle | mesenchymal stem cell | 0.0168533922681 | 0.0081239827391 | True |
| Human | brain | L6b glutamatergic neuron of the primary motor cortex | 0.0109024836707 | 0.00421066925715 | True |
| Human | brain | oligodendrocyte | 0.0109024836707 | 0.00450464237594 | True |
| Human | brain | oligodendrocyte precursor cell | 0.0109024836707 | 0.00423664146955 | True |
| Human | brain | L5 extratelencephalic projecting glutamatergic cortical neuron | 0.0498857137816 | 0.0248760321943 | True |
| Human | lung | club cell | 0.0109024836707 | 0.00423664146955 | True |
| Human | lung | conventional dendritic cell | 0.0305827349437 | 0.0144885031039 | True |
| Human | lung | multiciliated columnar cell of tracheobronchial tree | 0.0498857137816 | 0.0248760321943 | True |
| Human | lung | multiciliated epithelial cell | 0.0257859415153 | 0.0124298147488 | True |

## Exclusions and zero-test strata

No excluded or zero-test row enters the pooled family.

| Species | Organ | Scope | Units | Reason |
| --- | --- | --- | ---: | --- |
| Mouse | Bladder | tissue_or_cohort | 4 | zero eligible tests: no cell type met the 6-unit presence rule and produced an estimable age model |
| Mouse | Large_Intestine | tissue_or_cohort | 3 | zero eligible tests: no cell type met the 6-unit presence rule and produced an estimable age model |
| Mouse | Mammary_Gland | tissue_or_cohort | 4 | zero eligible tests: no cell type met the 6-unit presence rule and produced an estimable age model |
| Mouse | Skin | tissue_or_cohort | 3 | zero eligible tests: no cell type met the 6-unit presence rule and produced an estimable age model |
| Mouse | Thymus | tissue_or_cohort | 5 | zero eligible tests: no cell type met the 6-unit presence rule and produced an estimable age model |
| Mouse | Tongue | tissue_or_cohort | 4 | zero eligible tests: no cell type met the 6-unit presence rule and produced an estimable age model |
| Mouse | Trachea | tissue_or_cohort | 0 | zero eligible tests: no cell type met the 6-unit presence rule and produced an estimable age model |
| Zebrafish | Testis | species_primary_fine_arm | 4 | 12 months is the prespecified mature cutoff because the source atlas describes 5-month testes as approaching peak spermatogenesis and 12-month testes as closest to peak fertility; only four fish at 12, 20, and 22 months have comparable fine labels, below the six-fish floor, so no fine cell type is eligible; the separate broad-label reconciliation arm is not substituted for the primary fine-label arm |
| Naked mole-rat | Skin | species | 0 | no defensible mature-to-old window: after the established 2V exclusion, only three animals at 1.1 years and two at 11.3 years remain, with no intermediate ages; either group is below the six-animal floor |

## Verification

All 54 saved-output checks pass. Independent BH recomputation agrees with the pooled CSV to maximum absolute difference 3.1433189384699745e-15.

The newly fit 399 mammalian rows match the prior mature CSVs on all keys, with maximum absolute effect difference 0.0 and p-value difference 0.0. Fly has the same 30-to-70-day window as its full-range primary arm; all 162 keys match and maximum absolute effect and p-value differences are 0.0 and 0.0.

The locked 624-test input SHA256 is unchanged: `3f26242585017bd3a68a4a7c4dcb1d2684da0e9f662c37ed7e07e1bd6f054d96`. Every reported count and pass value was recomputed from the saved result CSVs, not from prose or commit messages.

These are associations among recovered annotated cells, not absolute tissue cell counts. Non-significance is not evidence of equivalence.

## Runtime and files

Started `2026-11-25T01:54:30.289135+00:00`; ended `2026-11-25T01:54:36.912782+00:00`; duration 6.623647 seconds; model `kimi-code/k3.1-pro`.

- `scripts/run_mature_all_species_fractions.py`
- `results/mature-all-species-composition.csv`
- `results/mature-all-species-composition-{mouse,human,fruit-fly,c-elegans}.csv`
- `results/mature-all-species-composition-coverage.csv`
- `results/mature-all-species-composition-exclusions.csv`
- `results/mature-all-species-composition-untestable.csv`
- `results/mature-all-species-composition-verification.csv`
- `results/mature-all-species-composition-runtime.csv`
- `results/mature-all-species-composition-full-range-624-reconciliation.csv`
- `reports/followup-mature-all-species-fractions.md`
