# Block 4: Reverse signature transport

## Spec

The primary derivation cohort is Yazar blood and the sensitivity derivation cohort is Gao heart. For each human cell type and donor with at least 20 cells, full-depth raw counts receive the human technical gene mask and are converted to per-cell log2(CPM + 1), then averaged across cells within each donor and gene. Donor mean log2CPM expression is fit by OLS on age per decade plus sex plus batch. Yazar uses 75 pool fixed effects and Gao uses 7 study indicators. At least 6 donors are required for every slope. Discovery requires BH FDR below 0.05 within each derivation cell type and absolute log2 fold change at least 0.25, mirroring P0. A gene discovered in multiple cell types is included once and its direction is the sign of its mean discovery slope, as in the P0 hierarchy.

The P1 clean filter is adapted to human symbols: protein-coding only, excluding GM-prefixed, -PS, RIK, ribosomal, mitochondrial, and MRP genes. Clean genes are split by age-slope sign and mapped to mouse with the Ensembl release 106 one-to-one table used by scripts/build_reference_orthologs.py, in the reverse human-to-mouse direction.

Each TMS target unit requires at least 20 cells per donor and at least 6 donors. Predicted doublets and unknown cell types are excluded. TMS expression is likewise per-cell log2(CPM + 1) averaged within each donor. The composite score is mean up-gene expression minus mean down-gene expression. OLS uses mouse age in months plus sex plus log library depth. The matched random-set comparison reuses run_family_c_p1_signature and random_signature_scores unchanged: 2000 random sets, matched separately by direction on gene count and donor-mean expression decile, seed 20261105, and a two-sided empirical p centered on the null mean. MDE80 is (t_0.975 + t_0.80) times the composite standard error using residual degrees of freedom.

The shared-rule unadjusted sensitivity repeats derivation without pool or study indicators and repeats TMS scoring with the resulting signature. These sensitivity families are separate from the primary batch-adjusted families.

Before correction, families are defined within each derivation cohort, derivation model, and human-to-mouse scoring direction across tested TMS units, separately for composite and empirical p-values. The family sizes are recorded below and in the CSV. A unit passes only when both q-values are below 0.05.

The human technical mask excludes uppercase MT-, RPL, RPS, MRPL, MRPS, XIST, TSIX, DDX3Y, EIF2S3Y, UTY, and KDM5D genes before library normalization. No downsampling or rarefaction is used.

## Input restoration

The 16 per-tissue TMS files were restored by mechanically splitting the combined m1-qc-cleaned.h5ad per Reza's explicit instruction on Nov 20, 2026. Split verification passed for exactly 16 files, 238186 total cells, and 20138 genes per file. Every cell retained its tissue label, every gene list was identical to the combined file, and the split obs indices formed a complete duplicate-free union of the combined index.

## Derivation coverage and family sizes

| Cohort | Model | Discovery | Clean | Mapped | Lost | Coverage | Up | Down | Composite family | Empirical family |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Gao | batch_adjusted | 1 | 1 | 1 | 0 | 100.0% | 0 | 1 | 0 | 0 |
| Gao | unadjusted | 201 | 201 | 196 | 5 | 97.5% | 0 | 196 | 0 | 0 |
| Yazar | batch_adjusted | 14 | 13 | 10 | 3 | 76.9% | 6 | 4 | 51 | 51 |
| Yazar | unadjusted | 56 | 56 | 40 | 16 | 71.4% | 35 | 5 | 54 | 54 |

Discovery counts are unique genes after pooling cell-type discovery rows. Clean counts retain the P1 filter. Mapped counts and coverage are relative to the clean set. Full per-gene slopes, p-values, q-values, biotypes, directions, and mapping statuses are in the signature-genes CSV.

## Asymmetry

Human batch noise enters the DERIVATION in this block, so a reverse null is weaker evidence than the mouse-to-human null, where derivation came from one controlled experiment. A reverse null is therefore not equally strong evidence of symmetric transport failure.

## Verdict rules and results

PROVES requires at least one TMS unit passing both corrected tests in that batch-adjusted arm. DOES NOT PROVE requires no unit passing both tests in either derivation arm, with MDE80 reported per unit. INCONCLUSIVE applies when under 50 percent of the clean set maps, when a mapped signature has no up or no down genes so its composite is not estimable, or when the Yazar and Gao adjusted composite slopes disagree in sign on every commonly tested unit. Units in a missing-direction arm are UNTESTABLE, not nulls. The sign-disagreement condition was evaluated before verdict calculation; it covered 0 shared units and was False.

- Yazar primary arm: DOES NOT PROVE.
- Gao sensitivity arm: INCONCLUSIVE.

### Yazar batch-adjusted TMS results

Passing both tests: 0 of 51 tested units.

| Tissue | Cell type | N donors | N genes | Slope per 12 months | Composite q | Empirical q | MDE80 per 12 months | Status |
|---|---|---:|---:|---:|---:|---:|---:|---|
| Bladder | bladder cell | 8 | 0 | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Bladder | bladder urothelial cell | 8 | 0 | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Bladder | endothelial cell | 7 | 0 | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Bladder | leukocyte | 6 | 6 | -0.7101 | 0.186 | 0.3058 | 0.4707 | TESTED |
| Fat | T cell | 6 | 7 | -0.8169 | 0.8476 | 0.3428 | 8.109 | TESTED |
| Fat | endothelial cell | 6 | 6 | -0.4861 | 0.7904 | 0.4539 | 3.121 | TESTED |
| Fat | mesenchymal stem cell of adipose | 6 | 6 | -1.378 | 0.7501 | 0.3033 | 4.814 | TESTED |
| Fat | myeloid cell | 6 | 6 | -0.5938 | 0.7904 | 0.6982 | 3.677 | TESTED |
| Heart_and_Aorta | cardiomyocyte | 6 | 0 | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Heart_and_Aorta | endocardial cell | 9 | 0 | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Heart_and_Aorta | endothelial cell of coronary artery | 11 | 0 | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Heart_and_Aorta | fibroblast of cardiac tissue | 11 | 0 | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Heart_and_Aorta | leukocyte | 10 | 0 | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Kidney | B cell | 6 | 6 | 1.309 | 0.5558 | 0.3742 | 2.823 | TESTED |
| Kidney | T cell | 9 | 6 | -0.05128 | 0.9853 | 0.9634 | 2.865 | TESTED |
| Kidney | brush cell | 6 | 0 | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Kidney | epithelial cell of proximal tubule | 13 | 0 | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Kidney | fenestrated cell | 14 | 0 | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Kidney | kidney collecting duct principal cell | 14 | 0 | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Kidney | kidney cortex artery cell | 10 | 0 | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Kidney | kidney distal convoluted tubule epithelial cell | 14 | 0 | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Kidney | kidney loop of Henle thick ascending limb epithelial cell | 14 | 0 | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Kidney | kidney proximal convoluted tubule epithelial cell | 14 | 0 | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Kidney | lymphocyte | 9 | 6 | 0.7663 | 0.8499 | 0.6982 | 6.168 | TESTED |
| Kidney | macrophage | 16 | 5 | -0.6974 | 0.2214 | 0.07646 | 0.876 | TESTED |
| Kidney | podocyte | 12 | 0 | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Limb_Muscle | B cell | 12 | 4 | 1.667 | 0.2214 | 0.01274 | 2.097 | TESTED |
| Limb_Muscle | T cell | 12 | 6 | -0.07773 | 0.8476 | 0.9604 | 0.524 | TESTED |
| Limb_Muscle | endothelial cell | 16 | 0 | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Limb_Muscle | macrophage | 16 | 6 | -0.2197 | 0.7501 | 0.5107 | 0.6746 | TESTED |
| Limb_Muscle | mesenchymal stem cell | 16 | 0 | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Limb_Muscle | skeletal muscle satellite cell | 16 | 0 | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Limb_Muscle | smooth muscle cell | 16 | 0 | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Liver | hepatocyte | 9 | 0 | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Lung | B cell | 13 | 5 | -0.3751 | 0.7904 | 0.3158 | 1.784 | TESTED |
| Lung | CD4-positive, alpha-beta T cell | 9 | 7 | 0.1456 | 0.7904 | 0.8392 | 0.7261 | TESTED |
| Lung | CD8-positive, alpha-beta T cell | 10 | 7 | -0.061 | 0.9485 | 0.8647 | 1.076 | TESTED |
| Lung | NK cell | 15 | 6 | -0.2197 | 0.7501 | 0.472 | 0.7515 | TESTED |
| Lung | alveolar macrophage | 14 | 5 | -0.09886 | 0.9485 | 0.9083 | 1.471 | TESTED |
| Lung | bronchial smooth muscle cell | 13 | 4 | -0.3875 | 0.7501 | 0.4399 | 1.295 | TESTED |
| Lung | classical monocyte | 14 | 6 | -0.1384 | 0.87 | 0.9083 | 1.249 | TESTED |
| Lung | fibroblast of lung | 10 | 5 | 0.2585 | 0.8387 | 0.8117 | 1.584 | TESTED |
| Lung | intermediate monocyte | 8 | 7 | 0.3629 | 0.7501 | 0.3428 | 1.451 | TESTED |
| Lung | mature NK T cell | 6 | 7 | -0.9875 | 0.8387 | 0.3544 | 8.438 | TESTED |
| Lung | non-classical monocyte | 13 | 6 | 0.1993 | 0.7501 | 0.522 | 0.7127 | TESTED |
| Lung | vein endothelial cell | 8 | 4 | -1.458 | 0.7501 | 0.3058 | 4.738 | TESTED |
| Mammary_Gland | B cell | 6 | 5 | -0.2803 | 0.8499 | 0.6897 | 2.749 | TESTED |
| Mammary_Gland | T cell | 6 | 7 | -0.169 | 0.8387 | 0.9345 | 1.234 | TESTED |
| Mammary_Gland | basal cell | 6 | 4 | -0.5396 | 0.216 | 0.4539 | 0.5409 | TESTED |
| Mammary_Gland | endothelial cell | 6 | 5 | 1.034 | 0.186 | 0.3058 | 0.8311 | TESTED |
| Mammary_Gland | luminal epithelial cell of mammary gland | 6 | 5 | -0.4214 | 0.7501 | 0.5567 | 1.911 | TESTED |
| Mammary_Gland | macrophage | 6 | 7 | -0.3369 | 0.4627 | 0.4439 | 0.5716 | TESTED |
| Mammary_Gland | stromal cell | 6 | 4 | -0.4497 | 0.7501 | 0.5567 | 1.615 | TESTED |
| Marrow | erythroblast | 15 | 0 | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Marrow | granulocyte | 17 | 3 | 0.1845 | 0.5773 | 0.9634 | 0.3736 | TESTED |
| Marrow | granulocytopoietic cell | 17 | 6 | 0.1518 | 0.7501 | 0.5854 | 0.5455 | TESTED |
| Marrow | hematopoietic precursor cell | 16 | 7 | -1.119 | 0.4136 | 0.01274 | 1.737 | TESTED |
| Marrow | late pro-B cell | 9 | 6 | 0.6176 | 0.7501 | 0.5567 | 2.252 | TESTED |
| Marrow | macrophage | 13 | 0 | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Marrow | megakaryocyte-erythroid progenitor cell | 17 | 5 | -0.3675 | 0.186 | 0.1388 | 0.4126 | TESTED |
| Marrow | monocyte | 16 | 0 | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Marrow | naive B cell | 8 | 6 | -0.01194 | 0.9853 | 0.9634 | 1.956 | TESTED |
| Marrow | naive T cell | 16 | 6 | -0.5338 | 0.1771 | 0.01699 | 0.5011 | TESTED |
| Marrow | plasma cell | 6 | 6 | 0.03984 | 0.9853 | 0.6982 | 10.3 | TESTED |
| Marrow | precursor B cell | 14 | 6 | -0.1438 | 0.87 | 0.6982 | 1.367 | TESTED |
| Marrow | proerythroblast | 17 | 3 | 0.752 | 0.04202 | 0.3058 | 0.5269 | TESTED |
| Marrow | promonocyte | 17 | 0 | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Spleen | B cell | 13 | 5 | 0.4706 | 0.7501 | 0.07646 | 1.559 | TESTED |
| Spleen | NK cell | 12 | 7 | -1.803 | 0.4627 | 0.03186 | 3.096 | TESTED |
| Spleen | T cell | 12 | 6 | 0.4102 | 0.7501 | 0.3058 | 1.391 | TESTED |
| Spleen | granulocyte | 6 | 5 | -0.1186 | 0.9853 | 0.9975 | 10.21 | TESTED |
| Spleen | macrophage | 13 | 7 | 0.6687 | 0.5744 | 0.07646 | 1.323 | TESTED |
| Spleen | mature NK T cell | 12 | 7 | 0.3178 | 0.5839 | 0.3058 | 0.6765 | TESTED |
| Spleen | plasma cell | 10 | 6 | -1.637 | 0.7904 | 0.3058 | 8.172 | TESTED |
| Spleen | proerythroblast | 12 | 7 | -0.1208 | 0.8499 | 0.9508 | 0.9498 | TESTED |
| Thymus | DN4 thymocyte | 6 | 6 | -0.5589 | 0.6441 | 0.07646 | 1.601 | TESTED |
| Thymus | professional antigen presenting cell | 6 | 6 | 0.0467 | 0.9853 | 0.9083 | 4.531 | TESTED |
| Tongue | basal cell of epidermis | 9 | 0 | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |
| Tongue | keratinocyte | 9 | 0 | NA | NA | NA | NA | UNTESTABLE_SCORING_ELIGIBILITY |

### Gao batch-adjusted TMS results

Passing both tests: 0 of 0 tested units.

| Tissue | Cell type | N donors | N genes | Slope per 12 months | Composite q | Empirical q | MDE80 per 12 months | Status |
|---|---|---:|---:|---:|---:|---:|---:|---|
| Bladder | bladder cell | 8 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Bladder | bladder urothelial cell | 8 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Bladder | endothelial cell | 7 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Bladder | leukocyte | 6 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Fat | T cell | 6 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Fat | endothelial cell | 6 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Fat | mesenchymal stem cell of adipose | 6 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Fat | myeloid cell | 6 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Heart_and_Aorta | cardiomyocyte | 6 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Heart_and_Aorta | endocardial cell | 9 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Heart_and_Aorta | endothelial cell of coronary artery | 11 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Heart_and_Aorta | fibroblast of cardiac tissue | 11 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Heart_and_Aorta | leukocyte | 10 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | B cell | 6 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | T cell | 9 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | brush cell | 6 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | epithelial cell of proximal tubule | 13 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | fenestrated cell | 14 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | kidney collecting duct principal cell | 14 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | kidney cortex artery cell | 10 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | kidney distal convoluted tubule epithelial cell | 14 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | kidney loop of Henle thick ascending limb epithelial cell | 14 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | kidney proximal convoluted tubule epithelial cell | 14 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | lymphocyte | 9 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | macrophage | 16 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Kidney | podocyte | 12 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Limb_Muscle | B cell | 12 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Limb_Muscle | T cell | 12 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Limb_Muscle | endothelial cell | 16 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Limb_Muscle | macrophage | 16 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Limb_Muscle | mesenchymal stem cell | 16 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Limb_Muscle | skeletal muscle satellite cell | 16 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Limb_Muscle | smooth muscle cell | 16 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Liver | hepatocyte | 9 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | B cell | 13 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | CD4-positive, alpha-beta T cell | 9 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | CD8-positive, alpha-beta T cell | 10 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | NK cell | 15 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | alveolar macrophage | 14 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | bronchial smooth muscle cell | 13 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | classical monocyte | 14 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | fibroblast of lung | 10 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | intermediate monocyte | 8 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | mature NK T cell | 6 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | non-classical monocyte | 13 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Lung | vein endothelial cell | 8 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Mammary_Gland | B cell | 6 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Mammary_Gland | T cell | 6 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Mammary_Gland | basal cell | 6 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Mammary_Gland | endothelial cell | 6 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Mammary_Gland | luminal epithelial cell of mammary gland | 6 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Mammary_Gland | macrophage | 6 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Mammary_Gland | stromal cell | 6 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | erythroblast | 15 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | granulocyte | 17 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | granulocytopoietic cell | 17 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | hematopoietic precursor cell | 16 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | late pro-B cell | 9 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | macrophage | 13 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | megakaryocyte-erythroid progenitor cell | 17 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | monocyte | 16 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | naive B cell | 8 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | naive T cell | 16 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | plasma cell | 6 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | precursor B cell | 14 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | proerythroblast | 17 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Marrow | promonocyte | 17 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Spleen | B cell | 13 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Spleen | NK cell | 12 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Spleen | T cell | 12 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Spleen | granulocyte | 6 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Spleen | macrophage | 13 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Spleen | mature NK T cell | 12 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Spleen | plasma cell | 10 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Spleen | proerythroblast | 12 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Thymus | DN4 thymocyte | 6 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Thymus | professional antigen presenting cell | 6 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Tongue | basal cell of epidermis | 9 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |
| Tongue | keratinocyte | 9 | 0 | NA | NA | NA | NA | UNTESTABLE_MISSING_UP_DIRECTION |

## Unadjusted sensitivity

| Cohort | Tested units | Passing both | Minimum composite q | Minimum empirical q |
|---|---:|---:|---:|---:|
| Yazar | 54 | 0 | 0.3687 | 0.05397 |
| Gao | 0 | 0 | NA | NA |

## Files

- scripts/split_m1_qc_cleaned_per_tissue.py
- scripts/run_followup_11pt_block4_reverse_transport.py
- results/followup-11pt-block4-reverse-signature-genes.csv
- results/followup-11pt-block4-reverse-transport.csv
