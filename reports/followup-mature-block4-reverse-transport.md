# Mature Block 4: Reverse signature transport

## Spec

This is a mature-to-old rerun of the frozen 11-point Block 4. Human derivation uses only donors aged at least 30 years before every model fit. Mouse scoring uses only donors aged at least 18 months. Both derivation and tested contrasts therefore satisfy the mature-to-old rule.

The original derivation is otherwise unchanged. For each human cell type and donor with at least 20 cells, full-depth raw counts received the human technical gene mask, were converted to per-cell log2(CPM + 1), and were averaged within donor and gene. The cached donor means were subset to 30-plus donors before OLS. Models use age per decade plus sex plus batch. Yazar uses the 75 pool fixed-effect candidates and Gao uses the 7 study indicators, with the original full-rank and residual-df handling. At least 6 donors are required per slope. Discovery requires BH FDR below 0.05 within each derivation cell type and absolute log2 fold change at least 0.25 per decade. A gene discovered in multiple cell types appears once with direction from the sign of its mean discovery slope.

The P1 clean filter retains protein-coding genes and excludes GM-prefixed, -PS, RIK, ribosomal, mitochondrial, and MRP genes. Mapping uses the Ensembl release 106 one-to-one table in the reverse human-to-mouse direction. No original cache was overwritten.

Two mouse scoring arms are prespecified for every derivation signature. The primary 18-plus arm fits continuous age in months plus sex plus log library depth over 18, 21, 24, and 30 months where present. The secondary 18-vs-30 arm uses only 18-month and 30-month donors and fits a binary covariate coded 0 at 18 months and 1 at 30 months, plus sex and log library depth. The coefficient is the adjusted 30-minus-18 contrast. Both arms require at least 20 cells per donor, at least 6 donors, exclusion of predicted doublets and unknown cell types, and the original expression and design floors. The pairwise arm additionally requires at least 3 donors on each side after the cell floor.

Composite scores are mean up-gene expression minus mean down-gene expression. Both arms reuse run_family_c_p1_signature and its matched random-set null unchanged: 2000 draws, direction-specific matching on gene count and donor-mean expression decile, seed 20261105, and a two-sided empirical p centered on the null mean. MDE80 is (t_0.975 + t_0.80) times composite SE using residual df.

BH correction is separate within each derivation cohort, derivation model, scoring arm, and test type, across tested TMS units only. A unit passes only when composite and empirical q are both below 0.05. The unadjusted sensitivity removes pool or study indicators during derivation, repeats both scoring arms, and remains in separate families.

The binary orientation was checked independently for every tested pairwise unit. An intercept-plus-age_group OLS with no covariates was fit to donor composite scores, and its coefficient was required to equal the raw mean score at 30 months minus the raw mean at 18 months to absolute tolerance 1e-10. The adjusted coefficient is also compared with that raw direction. Covariate adjustment can legitimately reverse a sign, so the exact equality assertion applies to the unadjusted orientation check.

All 30-month mice are male, as verified in the mature P1m derivation. Age and sex are not completely aliased because the 18-month group contains both sexes, but the 18-vs-30 contrast retains this age-sex confound and limited sex generalizability.

## Cache route and human donor coverage

The existing per-donor caches under tmp_organs/block4_reverse_transport_cache were reused. Both contained age_years and per-donor expression, so no h5ad re-streaming and no new cache directory were needed. Cached donor counts below are unique 30-plus donors with an eligible cached donor-cell-type mean. Fitted donors are the unique donors in cell types retaining at least 6 donors before model fitting.

| Cohort | 30-plus cached donors | Fitted donors | Cell types cached | Cell types fitted | Minimum fitted age | Maximum fitted age |
|---|---:|---:|---:|---:|---:|---:|
| Yazar | 932 | 932 | 19 | 19 | 30 | 97 |
| Gao | 109 | 109 | 13 | 13 | 30 | 91 |

## Derivation coverage and family sizes

| Cohort | Model | Scoring arm | Discovery | Clean | Mapped | Lost | Coverage | Up | Down | Composite family | Empirical family | Passing both |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Gao | batch_adjusted | 18-plus slope | 15 | 15 | 15 | 0 | 100.0% | 0 | 15 | 0 | 0 | 0 |
| Gao | batch_adjusted | 18-vs-30 | 15 | 15 | 15 | 0 | 100.0% | 0 | 15 | 0 | 0 | 0 |
| Gao | unadjusted | 18-plus slope | 1002 | 993 | 972 | 21 | 97.9% | 0 | 972 | 0 | 0 | 0 |
| Gao | unadjusted | 18-vs-30 | 1002 | 993 | 972 | 21 | 97.9% | 0 | 972 | 0 | 0 | 0 |
| Yazar | batch_adjusted | 18-plus slope | 7 | 7 | 5 | 2 | 71.4% | 2 | 3 | 38 | 38 | 0 |
| Yazar | batch_adjusted | 18-vs-30 | 7 | 7 | 5 | 2 | 71.4% | 2 | 3 | 17 | 17 | 0 |
| Yazar | unadjusted | 18-plus slope | 85 | 83 | 62 | 21 | 74.7% | 59 | 3 | 38 | 38 | 0 |
| Yazar | unadjusted | 18-vs-30 | 85 | 83 | 62 | 21 | 74.7% | 59 | 3 | 18 | 18 | 0 |

Discovery counts are unique genes after pooling cell-type discovery rows. Clean counts retain the P1 filter. Mapped counts and coverage are relative to the clean set. Per-gene mean slopes, minimum discovery p and q values, biotypes, directions, and mapping statuses are in the signature-genes CSV.

## Asymmetry

Human batch noise enters the DERIVATION in this block, so a reverse null is weaker evidence than the mouse-to-human null, where derivation came from one controlled experiment. A reverse null is therefore not equally strong evidence of symmetric transport failure.

## Verdict rules and sign disagreement

PROVES requires at least one TMS unit passing both corrected tests in that batch-adjusted cohort and scoring arm. DOES NOT PROVE requires no passing unit in that arm and per-unit MDE80 reporting. INCONCLUSIVE applies when under 50 percent of the clean set maps, when a mapped signature lacks an up or down direction, when cross-cohort adjusted slopes disagree in sign on every commonly tested unit, or when no unit is testable for another reason. Missing-direction units are UNTESTABLE, not nulls.

The sign-disagreement rule was evaluated before verdict calculation separately in each scoring arm.

| Scoring arm | Commonly tested units | Opposite signs on every shared unit |
|---|---:|---|
| 18-plus slope | 0 | False |
| 18-vs-30 | 0 | False |

| Derivation cohort | Scoring arm | Verdict | Tested units | Passing both |
|---|---|---|---:|---:|
| Yazar | 18-plus slope | DOES NOT PROVE | 38 | 0 |
| Yazar | 18-vs-30 | DOES NOT PROVE | 17 | 0 |
| Gao | 18-plus slope | INCONCLUSIVE | 0 | 0 |
| Gao | 18-vs-30 | INCONCLUSIVE | 0 | 0 |

### Yazar batch-adjusted 18-plus slope results

Passing both tests: 0 of 38 tested units. Composite and empirical family sizes are both 38.

| Tissue | Cell type | N donors | N genes | Slope per 12 months | SE | CI95 low | CI95 high | Composite p | Composite q | Empirical p | Empirical q | MDE80 per 12 months | Status |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| Fat | T cell | 6 | 3 | -0.528486 | 1.61177 | -7.46338 | 6.40641 | 0.774137 | 0.89059384368622829 | 0.532734 | 0.69806476072308665 | 8.64444 | TESTED |
| Fat | endothelial cell | 6 | 2 | 1.69716 | 1.10237 | -3.04597 | 6.44028 | 0.263551 | 0.60813891239076912 | 0.075962 | 0.26241424742174368 | 5.91236 | TESTED |
| Fat | mesenchymal stem cell of adipose | 6 | 2 | -0.802149 | 1.78861 | -8.49794 | 6.89364 | 0.697716 | 0.85526436203860756 | 0.428286 | 0.65251374312843569 | 9.5929 | TESTED |
| Fat | myeloid cell | 6 | 3 | -2.1576 | 0.42412 | -3.98245 | -0.332761 | 0.0365354 | 0.34708584060144282 | 0.0729635 | 0.26241424742174368 | 2.27469 | TESTED |
| Heart_and_Aorta | endocardial cell | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Heart_and_Aorta | endothelial cell of coronary artery | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Heart_and_Aorta | fibroblast of cardiac tissue | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Heart_and_Aorta | leukocyte | 8 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Kidney | B cell | 6 | 3 | 3.35812 | 1.05737 | -1.19139 | 7.90763 | 0.0864769 | 0.44469443890291149 | 0.0194903 | 0.1058042407367744 | 5.67103 | TESTED |
| Kidney | T cell | 9 | 3 | -0.0141107 | 1.55569 | -4.01315 | 3.98493 | 0.993114 | 0.9931137708515928 | 0.990505 | 0.99050474762618679 | 5.42957 | TESTED |
| Kidney | epithelial cell of proximal tubule | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Kidney | fenestrated cell | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Kidney | kidney collecting duct principal cell | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Kidney | kidney cortex artery cell | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Kidney | kidney distal convoluted tubule epithelial cell | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Kidney | kidney loop of Henle thick ascending limb epithelial cell | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Kidney | kidney proximal convoluted tubule epithelial cell | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Kidney | lymphocyte | 9 | 3 | 2.02081 | 1.95724 | -3.01042 | 7.05205 | 0.349177 | 0.62231474994571057 | 0.104948 | 0.28485757121439281 | 6.831 | TESTED |
| Kidney | macrophage | 11 | 2 | -1.00509 | 0.960291 | -3.27582 | 1.26564 | 0.33005 | 0.62231474994571057 | 0.455272 | 0.66539807019567132 | 3.13118 | TESTED |
| Kidney | podocyte | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Limb_Muscle | B cell | 8 | 2 | 1.27562 | 4.63653 | -11.5974 | 14.1487 | 0.796847 | 0.89059384368622829 | 0.474263 | 0.667481074277676 | 17.2359 | TESTED |
| Limb_Muscle | T cell | 8 | 3 | 0.190043 | 2.04039 | -5.47498 | 5.85507 | 0.930271 | 0.955413012027024 | 0.806097 | 0.85088011549780673 | 7.58496 | TESTED |
| Limb_Muscle | endothelial cell | 12 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Limb_Muscle | macrophage | 12 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Limb_Muscle | mesenchymal stem cell | 12 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Limb_Muscle | skeletal muscle satellite cell | 12 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Limb_Muscle | smooth muscle cell | 12 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Lung | B cell | 9 | 3 | 1.50693 | 0.711807 | -0.32283 | 3.33669 | 0.0878338 | 0.44469443890291149 | 0.0549725 | 0.26111944027986 | 2.4843 | TESTED |
| Lung | CD4-positive, alpha-beta T cell | 6 | 3 | 1.69227 | 1.6059 | -5.21735 | 8.6019 | 0.402498 | 0.62231474994571057 | 0.262869 | 0.55494474984729869 | 8.61293 | TESTED |
| Lung | CD8-positive, alpha-beta T cell | 8 | 3 | 0.496499 | 0.390109 | -0.586616 | 1.57961 | 0.272062 | 0.60813891239076912 | 0.526237 | 0.69806476072308665 | 1.45019 | TESTED |
| Lung | NK cell | 9 | 3 | 0.670675 | 0.587871 | -0.840494 | 2.18184 | 0.30561 | 0.61121981081077303 | 0.412794 | 0.65251374312843569 | 2.05174 | TESTED |
| Lung | alveolar macrophage | 9 | 2 | -0.954022 | 0.678873 | -2.69912 | 0.791076 | 0.218916 | 0.58151776329770455 | 0.35982 | 0.65251374312843569 | 2.36935 | TESTED |
| Lung | bronchial smooth muscle cell | 9 | 2 | -1.50676 | 0.959581 | -3.97344 | 0.959927 | 0.177161 | 0.56100875446716647 | 0.003998 | 0.037981009495252301 | 3.34906 | TESTED |
| Lung | classical monocyte | 10 | 2 | -1.93112 | 0.923733 | -4.19141 | 0.329178 | 0.0815356 | 0.44469443890291149 | 0.013993 | 0.088622355488922194 | 3.09692 | TESTED |
| Lung | intermediate monocyte | 6 | 3 | -1.35973 | 1.95851 | -9.7865 | 7.06704 | 0.559318 | 0.77473389853148533 | 0.277861 | 0.55572213893053468 | 10.5041 | TESTED |
| Lung | mature NK T cell | 6 | 3 | -2.21124 | 3.29094 | -16.371 | 11.9485 | 0.570857 | 0.77473389853148533 | 0.0724638 | 0.26241424742174368 | 17.6503 | TESTED |
| Lung | non-classical monocyte | 9 | 3 | -0.817423 | 0.585816 | -2.32331 | 0.688466 | 0.221714 | 0.58151776329770455 | 0.63918 | 0.73602592643072406 | 2.04457 | TESTED |
| Lung | vein endothelial cell | 7 | 2 | -1.50953 | 1.50449 | -6.29749 | 3.27844 | 0.389621 | 0.62231474994571057 | 0.177911 | 0.42253873063468272 | 6.26007 | TESTED |
| Marrow | erythroblast | 11 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Marrow | granulocyte | 13 | 2 | 0.652522 | 0.713699 | -0.961976 | 2.26702 | 0.384408 | 0.62231474994571057 | 0.707646 | 0.79089866831290234 | 2.24498 | TESTED |
| Marrow | granulocytopoietic cell | 13 | 2 | 1.71324 | 0.995769 | -0.539348 | 3.96582 | 0.119449 | 0.50433848606339471 | 0.37981 | 0.65251374312843569 | 3.13225 | TESTED |
| Marrow | hematopoietic precursor cell | 12 | 3 | -5.93278 | 1.7737 | -10.023 | -1.84261 | 0.0101576 | 0.19299483544261181 | 0.003998 | 0.037981009495252301 | 5.6668 | TESTED |
| Marrow | late pro-B cell | 8 | 3 | -1.28798 | 1.39927 | -5.17297 | 2.59701 | 0.409418 | 0.62231474994571057 | 0.382809 | 0.65251374312843569 | 5.20165 | TESTED |
| Marrow | macrophage | 10 | 2 | -0.336377 | 0.655552 | -1.94045 | 1.2677 | 0.626207 | 0.82054649912355593 | 0.429285 | 0.65251374312843569 | 2.19781 | TESTED |
| Marrow | megakaryocyte-erythroid progenitor cell | 13 | 3 | -0.908107 | 0.704513 | -2.50183 | 0.685613 | 0.229546 | 0.58151776329770455 | 0.171914 | 0.42253873063468272 | 2.21609 | TESTED |
| Marrow | monocyte | 12 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Marrow | naive B cell | 6 | 3 | -3.76669 | 7.62962 | -36.5943 | 29.0609 | 0.670412 | 0.849188579215157 | 0.604198 | 0.71748500749625188 | 40.9201 | TESTED |
| Marrow | naive T cell | 12 | 3 | -1.24071 | 0.247325 | -1.81105 | -0.670381 | 0.00103133 | 0.039190571568941797 | 0.0889555 | 0.28169248708978839 | 0.790178 | TESTED |
| Marrow | precursor B cell | 10 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Marrow | proerythroblast | 13 | 2 | 0.67993 | 0.603423 | -0.685107 | 2.04497 | 0.288977 | 0.61006358858376908 | 0.578711 | 0.71306282342699623 | 1.8981 | TESTED |
| Marrow | promonocyte | 13 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Spleen | B cell | 10 | 3 | 6.13607 | 3.0823 | -1.40606 | 13.6782 | 0.0936199 | 0.44469443890291149 | 0.00049975 | 0.018990504747626102 | 10.3338 | TESTED |
| Spleen | NK cell | 9 | 3 | -2.64663 | 3.82226 | -12.4721 | 7.1788 | 0.519506 | 0.75927767836810867 | 0.197901 | 0.44236705176823349 | 13.3402 | TESTED |
| Spleen | T cell | 9 | 3 | -2.37916 | 1.49213 | -6.2148 | 1.45649 | 0.171714 | 0.56100875446716647 | 0.0994503 | 0.28485757121439281 | 5.20772 | TESTED |
| Spleen | granulocyte | 6 | 3 | -0.189114 | 1.53695 | -6.80206 | 6.42383 | 0.913321 | 0.955413012027024 | 0.956022 | 0.98186042114078098 | 8.24312 | TESTED |
| Spleen | macrophage | 10 | 3 | 6.23454 | 3.75737 | -2.95941 | 15.4285 | 0.148129 | 0.56100875446716647 | 0.0009995 | 0.018990504747626102 | 12.597 | TESTED |
| Spleen | mature NK T cell | 9 | 3 | 0.306122 | 1.11216 | -2.55278 | 3.16503 | 0.794137 | 0.89059384368622829 | 0.792104 | 0.85088011549780673 | 3.88159 | TESTED |
| Spleen | plasma cell | 10 | 3 | -0.340541 | 1.83754 | -4.83683 | 4.15575 | 0.85908 | 0.93271559412261962 | 0.581709 | 0.71306282342699623 | 6.16055 | TESTED |
| Spleen | proerythroblast | 9 | 3 | 4.97157 | 1.65627 | 0.713992 | 9.22916 | 0.0300417 | 0.34708584060144282 | 0.0089955 | 0.068365817091454201 | 5.7806 | TESTED |

### Yazar batch-adjusted 18-vs-30 results

Passing both tests: 0 of 17 tested units. Composite and empirical family sizes are both 17.

| Tissue | Cell type | N donors | N genes | 30-minus-18 effect | SE | CI95 low | CI95 high | Composite p | Composite q | Empirical p | Empirical q | MDE80 30-minus-18 | Status |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| Fat | T cell | 4 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Fat | endothelial cell | 4 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Fat | mesenchymal stem cell of adipose | 4 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Fat | myeloid cell | 4 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Heart_and_Aorta | endocardial cell | 3 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Heart_and_Aorta | endothelial cell of coronary artery | 5 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Heart_and_Aorta | fibroblast of cardiac tissue | 5 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Heart_and_Aorta | leukocyte | 4 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Kidney | B cell | 3 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Kidney | T cell | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Kidney | epithelial cell of proximal tubule | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Kidney | fenestrated cell | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Kidney | kidney collecting duct principal cell | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Kidney | kidney cortex artery cell | 4 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Kidney | kidney distal convoluted tubule epithelial cell | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Kidney | kidney loop of Henle thick ascending limb epithelial cell | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Kidney | kidney proximal convoluted tubule epithelial cell | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Kidney | lymphocyte | 5 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Kidney | macrophage | 7 | 2 | -1.94044 | 0.986024 | -5.0784 | 1.19753 | 0.143735 | 0.38143583990132363 | 0.065967 | 0.2412793603198401 | 4.10276 | TESTED |
| Kidney | podocyte | 5 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Limb_Muscle | B cell | 2 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Limb_Muscle | T cell | 3 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Limb_Muscle | endothelial cell | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Limb_Muscle | macrophage | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Limb_Muscle | mesenchymal stem cell | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Limb_Muscle | skeletal muscle satellite cell | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Limb_Muscle | smooth muscle cell | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Lung | B cell | 7 | 3 | 1.44939 | 0.897091 | -1.40555 | 4.30434 | 0.204584 | 0.38643558888281532 | 0.0709645 | 0.2412793603198401 | 3.73272 | TESTED |
| Lung | CD4-positive, alpha-beta T cell | 4 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Lung | CD8-positive, alpha-beta T cell | 6 | 3 | 0.582384 | 0.559167 | -1.82352 | 2.98829 | 0.406997 | 0.53222734457260601 | 0.481759 | 0.58499321767687584 | 2.99899 | TESTED |
| Lung | NK cell | 7 | 3 | 0.731812 | 0.715675 | -1.54578 | 3.00941 | 0.381783 | 0.53222734457260601 | 0.405797 | 0.5748792270531401 | 2.97786 | TESTED |
| Lung | alveolar macrophage | 7 | 2 | -1.02386 | 0.829334 | -3.66318 | 1.61545 | 0.304878 | 0.4711748010845056 | 0.358321 | 0.5748792270531401 | 3.45079 | TESTED |
| Lung | bronchial smooth muscle cell | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Lung | classical monocyte | 8 | 2 | -2.026 | 1.14351 | -5.2009 | 1.1489 | 0.151134 | 0.38143583990132363 | 0.027986 | 0.23788105947026489 | 4.25091 | TESTED |
| Lung | intermediate monocyte | 6 | 3 | -1.35973 | 1.95851 | -9.7865 | 7.06704 | 0.559318 | 0.67917194449714924 | 0.277861 | 0.5748792270531401 | 10.5041 | TESTED |
| Lung | mature NK T cell | 5 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Lung | non-classical monocyte | 7 | 3 | -0.966283 | 0.763319 | -3.39551 | 1.46294 | 0.294921 | 0.4711748010845056 | 0.525237 | 0.59526903215059135 | 3.17611 | TESTED |
| Lung | vein endothelial cell | 5 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Marrow | erythroblast | 6 | 2 | -2.43988 | 1.28994 | -7.99006 | 3.11031 | 0.19911 | 0.38643558888281532 | 0.0994503 | 0.281775778777278 | 6.91838 | TESTED |
| Marrow | granulocyte | 7 | 2 | 0.396013 | 0.708206 | -1.85782 | 2.64984 | 0.61507 | 0.69707952648719707 | 0.98051 | 0.98050974512743638 | 2.94679 | TESTED |
| Marrow | granulocytopoietic cell | 7 | 2 | 1.61792 | 0.442167 | 0.210744 | 3.02509 | 0.0352682 | 0.19985303440784299 | 0.330335 | 0.5748792270531401 | 1.83982 | TESTED |
| Marrow | hematopoietic precursor cell | 6 | 3 | -5.5586 | 0.639969 | -8.31216 | -2.80504 | 0.0129973 | 0.19985303440784299 | 0.00549725 | 0.093453273363318295 | 3.43235 | TESTED |
| Marrow | late pro-B cell | 5 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Marrow | macrophage | 6 | 2 | -0.417753 | 0.981209 | -4.63956 | 3.80405 | 0.711727 | 0.75621046636452338 | 0.402299 | 0.5748792270531401 | 5.26253 | TESTED |
| Marrow | megakaryocyte-erythroid progenitor cell | 7 | 2 | 0.0464546 | 0.941346 | -2.94933 | 3.04224 | 0.963743 | 0.963742828043192 | 0.790605 | 0.84001749125437275 | 3.91687 | TESTED |
| Marrow | monocyte | 6 | 2 | -1.72806 | 0.700162 | -4.74061 | 1.28449 | 0.132345 | 0.38143583990132363 | 0.045977 | 0.2412793603198401 | 3.75519 | TESTED |
| Marrow | naive B cell | 4 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Marrow | naive T cell | 7 | 3 | -1.17179 | 0.282725 | -2.07155 | -0.272031 | 0.0255119 | 0.19985303440784299 | 0.172414 | 0.41871921182266009 | 1.1764 | TESTED |
| Marrow | precursor B cell | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Marrow | proerythroblast | 7 | 2 | 1.21853 | 0.648971 | -0.846784 | 3.28385 | 0.157062 | 0.38143583990132363 | 0.444778 | 0.5816322607926806 | 2.70032 | TESTED |
| Marrow | promonocyte | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Spleen | B cell | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Spleen | NK cell | 5 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Spleen | T cell | 5 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Spleen | granulocyte | 3 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Spleen | macrophage | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Spleen | mature NK T cell | 5 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Spleen | plasma cell | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |
| Spleen | proerythroblast | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_18_VS_30_DONOR_FLOOR |

### Gao batch-adjusted 18-plus slope results

Passing both tests: 0 of 0 tested units. Composite and empirical family sizes are both 0.

| Tissue | Cell type | N donors | N genes | Slope per 12 months | SE | CI95 low | CI95 high | Composite p | Composite q | Empirical p | Empirical q | MDE80 per 12 months | Status |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| Fat | T cell | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Fat | endothelial cell | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Fat | mesenchymal stem cell of adipose | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Fat | myeloid cell | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Heart_and_Aorta | endocardial cell | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Heart_and_Aorta | endothelial cell of coronary artery | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Heart_and_Aorta | fibroblast of cardiac tissue | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Heart_and_Aorta | leukocyte | 8 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | B cell | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | T cell | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | epithelial cell of proximal tubule | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | fenestrated cell | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | kidney collecting duct principal cell | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | kidney cortex artery cell | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | kidney distal convoluted tubule epithelial cell | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | kidney loop of Henle thick ascending limb epithelial cell | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | kidney proximal convoluted tubule epithelial cell | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | lymphocyte | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | macrophage | 11 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | podocyte | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Limb_Muscle | B cell | 8 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Limb_Muscle | T cell | 8 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Limb_Muscle | endothelial cell | 12 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Limb_Muscle | macrophage | 12 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Limb_Muscle | mesenchymal stem cell | 12 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Limb_Muscle | skeletal muscle satellite cell | 12 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Limb_Muscle | smooth muscle cell | 12 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | B cell | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | CD4-positive, alpha-beta T cell | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | CD8-positive, alpha-beta T cell | 8 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | NK cell | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | alveolar macrophage | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | bronchial smooth muscle cell | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | classical monocyte | 10 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | intermediate monocyte | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | mature NK T cell | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | non-classical monocyte | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | vein endothelial cell | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | erythroblast | 11 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | granulocyte | 13 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | granulocytopoietic cell | 13 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | hematopoietic precursor cell | 12 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | late pro-B cell | 8 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | macrophage | 10 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | megakaryocyte-erythroid progenitor cell | 13 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | monocyte | 12 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | naive B cell | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | naive T cell | 12 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | precursor B cell | 10 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | proerythroblast | 13 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | promonocyte | 13 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Spleen | B cell | 10 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Spleen | NK cell | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Spleen | T cell | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Spleen | granulocyte | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Spleen | macrophage | 10 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Spleen | mature NK T cell | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Spleen | plasma cell | 10 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Spleen | proerythroblast | 9 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |

### Gao batch-adjusted 18-vs-30 results

Passing both tests: 0 of 0 tested units. Composite and empirical family sizes are both 0.

| Tissue | Cell type | N donors | N genes | 30-minus-18 effect | SE | CI95 low | CI95 high | Composite p | Composite q | Empirical p | Empirical q | MDE80 30-minus-18 | Status |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| Fat | T cell | 4 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Fat | endothelial cell | 4 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Fat | mesenchymal stem cell of adipose | 4 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Fat | myeloid cell | 4 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Heart_and_Aorta | endocardial cell | 3 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Heart_and_Aorta | endothelial cell of coronary artery | 5 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Heart_and_Aorta | fibroblast of cardiac tissue | 5 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Heart_and_Aorta | leukocyte | 4 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | B cell | 3 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | T cell | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | epithelial cell of proximal tubule | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | fenestrated cell | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | kidney collecting duct principal cell | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | kidney cortex artery cell | 4 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | kidney distal convoluted tubule epithelial cell | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | kidney loop of Henle thick ascending limb epithelial cell | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | kidney proximal convoluted tubule epithelial cell | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | lymphocyte | 5 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | macrophage | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | podocyte | 5 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Limb_Muscle | B cell | 2 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Limb_Muscle | T cell | 3 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Limb_Muscle | endothelial cell | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Limb_Muscle | macrophage | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Limb_Muscle | mesenchymal stem cell | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Limb_Muscle | skeletal muscle satellite cell | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Limb_Muscle | smooth muscle cell | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | B cell | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | CD4-positive, alpha-beta T cell | 4 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | CD8-positive, alpha-beta T cell | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | NK cell | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | alveolar macrophage | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | bronchial smooth muscle cell | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | classical monocyte | 8 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | intermediate monocyte | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | mature NK T cell | 5 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | non-classical monocyte | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | vein endothelial cell | 5 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | erythroblast | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | granulocyte | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | granulocytopoietic cell | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | hematopoietic precursor cell | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | late pro-B cell | 5 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | macrophage | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | megakaryocyte-erythroid progenitor cell | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | monocyte | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | naive B cell | 4 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | naive T cell | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | precursor B cell | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | proerythroblast | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | promonocyte | 7 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Spleen | B cell | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Spleen | NK cell | 5 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Spleen | T cell | 5 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Spleen | granulocyte | 3 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Spleen | macrophage | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Spleen | mature NK T cell | 5 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Spleen | plasma cell | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Spleen | proerythroblast | 6 | 0 | NA | NA | NA | NA | NA | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |

## MDE80 by arm

| Cohort | Model | Scoring arm | Tested units | Median MDE80 | Minimum MDE80 | Maximum MDE80 | Units |
|---|---|---|---:|---:|---:|---:|---|
| Gao | batch_adjusted | 18-plus slope | 0 | NA | NA | NA | log2CPM per month |
| Gao | batch_adjusted | 18-vs-30 | 0 | NA | NA | NA | log2CPM for 30 minus 18 |
| Gao | unadjusted | 18-plus slope | 0 | NA | NA | NA | log2CPM per month |
| Gao | unadjusted | 18-vs-30 | 0 | NA | NA | NA | log2CPM for 30 minus 18 |
| Yazar | batch_adjusted | 18-plus slope | 38 | 0.462349 | 0.0658482 | 3.41 | log2CPM per month |
| Yazar | batch_adjusted | 18-vs-30 | 17 | 3.45079 | 1.1764 | 10.5041 | log2CPM for 30 minus 18 |
| Yazar | unadjusted | 18-plus slope | 38 | 0.301656 | 0.0568551 | 3.20691 | log2CPM per month |
| Yazar | unadjusted | 18-vs-30 | 18 | 3.41627 | 1.04457 | 9.30425 | log2CPM for 30 minus 18 |

These MDE80 values are the smallest composite effects with 80 percent power at the two-sided 0.05 level under the fitted t model. Large values show reduced resolution after age restriction. Non-passage is not equivalence or proof of no mature effect.

## Unadjusted sensitivity

| Cohort | Scoring arm | Tested units | Passing both | Minimum composite q | Minimum empirical q | Median MDE80 | MDE80 range |
|---|---|---:|---:|---:|---:|---:|---|
| Yazar | 18-plus slope | 38 | 0 | 0.1142620126837398 | 0.018990504747626102 | 0.301656 | 0.0568551 to 3.20691 |
| Yazar | 18-vs-30 | 18 | 0 | 0.36279958854853178 | 0.0089955022488755008 | 3.41627 | 1.04457 to 9.30425 |
| Gao | 18-plus slope | 0 | 0 | NA | NA | NA | NA |
| Gao | 18-vs-30 | 0 | 0 | NA | NA | NA | NA |

## Comparison with original Block 4

The original Yazar batch-adjusted derivation gave 14 discovered genes, 13 clean genes, 10 mapped genes, 6 up and 4 down, with 0 of 51 TMS units passing both corrected tests and a DOES NOT PROVE verdict. The original Gao batch-adjusted derivation gave only 1 down gene and no up gene, so that arm was UNTESTABLE_MISSING_UP_DIRECTION and INCONCLUSIVE, never a null.

| Original result | Mature scoring arm | Comparison | Evidence |
|---|---|---|---|
| Original Yazar 0 of 51 both-pass | 18-plus slope | SURVIVED | Mature result is 0 of 38 both-pass; No mature unit passes both corrected tests. |
| Original Yazar 0 of 51 both-pass | 18-vs-30 | SURVIVED | Mature result is 0 of 17 both-pass; No mature unit passes both corrected tests. |
| Original Gao 1 down, 0 up, inestimable | 18-plus slope | UNTESTABLE | Mature Gao still lacks a fully estimable mapped bidirectional signature: 0 up, 15 down, coverage 100.0%. |
| Original Gao 1 down, 0 up, inestimable | 18-vs-30 | UNTESTABLE | Mature Gao still lacks a fully estimable mapped bidirectional signature: 0 up, 15 down, coverage 100.0%. |

Mature family sizes differ from the original 51-unit family, so the comparison concerns the verdict and passage rule, not identical denominators or identical donor sets.

## Verification

BH q-values were recomputed from raw p-values in the serialized transport CSV with statsmodels multipletests. Maximum absolute discrepancy: composite 1.1934897514720433e-15; empirical 1.0998146837692957e-15. Passing counts and family sizes in this report were recomputed from that CSV.

Minimum fitted human ages: Yazar 30 years; Gao 30 years. Minimum fitted mouse age across TESTED rows: 18 months. Pairwise TESTED rows contain exactly ages 18 and 30 months and satisfy at least 6 total donors with at least 3 per side.

All 111 TESTED rows have n_random_sets = 2000 and random_seed = 20261105. The pairwise orientation check covered 35 units with maximum absolute unadjusted-coefficient error 1.7763568394002505e-15. Adjusted and raw mean-difference signs agreed in 23 of 35 pairwise units.

Forbidden Unicode em dash and en dash characters are rejected by the CSV and Markdown writers. A final source and output character scan is performed after the run.

## Runtime

Model: kimi-code/k3.1-pro. All matched nulls use 2000 iterations and seed 20261105. Exact stage start and end timestamps are in the runtime CSV.

| Stage | Seconds |
|---|---:|
| load_human_biotypes | 8.470 |
| load_and_restrict_Yazar_cache | 4.728 |
| derive_Yazar_batch_adjusted | 5.089 |
| derive_Yazar_unadjusted | 3.207 |
| load_and_restrict_Gao_cache | 0.145 |
| derive_Gao_batch_adjusted | 0.223 |
| derive_Gao_unadjusted | 0.225 |
| prepare_18plus_mouse_targets | 26.183 |
| score_Gao_batch_adjusted_mature_18plus_slope | 0.008 |
| score_Gao_unadjusted_mature_18plus_slope | 0.006 |
| score_Yazar_batch_adjusted_mature_18plus_slope | 15.014 |
| score_Yazar_unadjusted_mature_18plus_slope | 25.405 |
| score_Gao_batch_adjusted_mature_18_vs_30 | 0.063 |
| score_Gao_unadjusted_mature_18_vs_30 | 0.067 |
| score_Yazar_batch_adjusted_mature_18_vs_30 | 6.705 |
| score_Yazar_unadjusted_mature_18_vs_30 | 11.921 |
| write_gene_and_transport_csv | 0.014 |
| verify_serialized_statistics_and_ages | 0.009 |

## Files

- scripts/run_mature_block4_reverse_transport.py
- results/followup-mature-block4-reverse-signature-genes.csv
- results/followup-mature-block4-reverse-transport.csv
- results/followup-mature-block4-runtime.csv
- reports/followup-mature-block4-reverse-transport.md
