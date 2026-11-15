# Dataset 09: Human Retinal Aging Cohort (Ophthalmology Science 2026, PMC12892068, HRA014708 / Figshare 28757048)

## Verdict
Across the 12 unique Table 1 donors (18 retinal samples, ages 34 to 92 yr; Table 1 sex 8M / 4F, expression-verified sex 7M / 5F; 4 DR donors, 2 DM-without-DR donors, 6 non-diabetic donors), with paired eyes collapsed to the donor and diabetic status adjusted, the primary suite `All12_DiseaseAdj` gave 0 of 11 cell-class fraction shifts at q_BH < 0.05, 0 of 164089 global gene-level age hits at q_BH < 0.05, cross-cell-type concordance r = +0.3303, one-sided permutation p = 0.0110, and 0 of 48 `P1_Clean` signature tests (all suites) at q_BH < 0.05 on both the analytic and the 200-draw matched null (0 in the positive, transport-supporting direction; min q_BH analytic = 0.1715, min q_BH empirical = 0.0100). With n = 12 and heavy disease and surgical confounding, these are low-power results and do not establish absence of age effects.

## Notable Results and Caveats
- Cross-cell-type concordance reached one-sided permutation p < 0.05 in `All12_DiseaseAdj` (r = +0.3303, p = 0.0110), `All12_DiseaseAdj_ExprSex` (r = +0.3349, p = 0.0110), `All12_DiseaseLiveAdj` (r = +0.2846, p = 0.0020). These p values are not corrected across the 5 concordance tests in this report, and the two donor subsets that remove the disease or surgical strata (`NonDM6_AgeOnly` p = 0.1429, `Cadaver7_DiseaseAdj` p = 0.1608) did not reach p < 0.05. Shared donor-level factors that the covariates capture only coarsely (diabetic retinopathy, enucleation for ocular disease, cadaver handling) and ambient RNA from the dominant rod population can align age effects across marker-scored classes. This is a flag for review against the paper's human concordance result (p = 0.70 in the surveyed human datasets), not a replication of concordance: it comes from one 12-donor cohort with disease confounding.
- Single Family B hit: `Cadaver7_DiseaseAdj` Astrocyte MID1IP1 (effect -0.3611 log2CPM per decade, p = 4.92e-08, q_BH global = 0.0073, n = 7, residual df = 2). With residual df = 2 this is treated as unstable and is not replicated in the 12-donor suites.
- `P1_Clean` groups significant on the empirical null only (analytic q_BH >= 0.05): `All12_WholeRetina_DiseaseAdj` All_Retina_QC (composite -0.0603, q_BH empirical = 0.0100, q_BH analytic = 0.1715); `NonDM6_AgeOnly` Rod (composite -0.0835, q_BH empirical = 0.0199, q_BH analytic = 0.4259); `NonDM6_AgeOnly` HC (composite -0.2060, q_BH empirical = 0.0199, q_BH analytic = 0.4259); `NonDM6_AgeOnly` Cone (composite -0.1067, q_BH empirical = 0.0265, q_BH analytic = 0.4259). All of these are in the negative direction (opposite to mouse-signature transport) and none passes both criteria.

## Primary Source and Dual-Accession Audit
- Publication: Ophthalmology Science (2026), DOI `10.1016/j.xops.2025.101062`, PMID `41684508`, PMCID `PMC12892068`.
- GSA-Human `HRA014708`: controlled access. The browse page returned an anti-bot SecurityCheck interstitial; it was not bypassed and no access application was submitted.
- Figshare `doi:10.6084/m9.figshare.28757048`: openly hosts 18 Cell Ranger v8.0.1 count matrix triplets (`sample_1` to `sample_18`) plus `README.docx` and `Script.py`. The author-clustered `/RPE.h5ad` referenced by `Script.py` is not deposited. Counts are integer UMIs: 293,480 raw barcodes, 285,878 cells passing the paper's QC (at least 200 genes, more than 500 UMIs, under 15% mitochondrial UMIs).
- Cell classes were assigned by canonical marker scoring informed by Figure 1C (mean log1p(CP10k) of the marker set at least 0.5; otherwise Unassigned), because the author labels are not deposited. An unpartitioned whole-retina donor pseudobulk (`All12_WholeRetina_DiseaseAdj`) is reported as a label-free control.

## Metadata Discrepancies Found (reported, not resolved by assumption)
1. Donor sex: expression-inferred sex (mean log1p(CP10k) of XIST versus mean of RPS4Y1 and DDX3Y) matched Table 1 on 17 of 18 samples. Mismatch: donor D2 (sample(s) 2, 55 yr): Table 1 sex M, expression sex F (QC-cell pseudobulk XIST = 19198 UMIs, TSIX = 915 UMIs, RPS4Y1 + DDX3Y + UTY + KDM5D + EIF1AY = 0 UMIs). Primary models use Table 1 sex as printed; `All12_DiseaseAdj_ExprSex` repeats the primary model with expression-verified sex.
2. Cohort counts: the Figshare description states 18 retinal samples from 13 Chinese donors (8 male, 5 female, ages 34 to 92 years) with 7 DR, 3 diabetics without retinopathy, and 8 non-diabetic samples, whereas Table 1 of PMC12892068 lists 18 samples from 12 donors (D1 to D12; Table 1 sex 8M / 4F) with 8 DR, 3 DM, and 7 non-DM samples. The source of the difference cannot be determined from the deposited files. This analysis uses the Table 1 donor mapping (12 donor units); the cohort-level counts here therefore follow Table 1, and the sample-level sex check above is the only independent check applied.

## Donor Summary (n = 12 donors across 18 eyes)
| Donor | Age (yr) | Sex (Table 1) | Sex (expression) | XIST UMIs (QC cells) | Y-gene UMIs (5 genes) | Diabetic Status | Source | Eyes (Sample #) | QC Cells | Surgical Indication or Cause of Death |
|---|---:|:---:|:---:|---:|---:|---|---|---|---:|---|
| D6 | 34 | M | M | 45 | 60013 | Non-DM | Cadaver | 1 (6) | 16206 | Cerebral hemorrhage |
| D10 | 35 | M | M | 80 | 75848 | DR | Cadaver | 2 (13;14) | 36588 | Traumatic shock |
| D4 | 50 | F | F | 21581 | 4 | Non-DM | Live | 1 (4) | 6940 | Phthisis bulbi |
| D2 | 55 | M | F | 19198 | 0 | Non-DM | Live | 1 (2) | 11017 | Phthisis bulbi |
| D8 | 55 | F | F | 41091 | 13 | DM | Cadaver | 2 (9;10) | 33473 | Traumatic brain injury |
| D9 | 55 | M | M | 43 | 97348 | DR | Cadaver | 2 (11;12) | 36065 | Brainstem hemorrhage |
| D5 | 68 | F | F | 13664 | 0 | Non-DM | Live | 1 (5) | 6641 | Corneal perforation |
| D12 | 69 | M | M | 96 | 117434 | DR | Cadaver | 2 (17;18) | 32605 | Lung cancer |
| D3 | 71 | M | M | 20 | 18123 | Non-DM | Live | 1 (3) | 6813 | Corneal ulcer |
| D1 | 76 | M | M | 53 | 61848 | DM | Live | 1 (1) | 17107 | Corneal staphyloma |
| D11 | 84 | M | M | 72 | 107404 | DR | Cadaver | 2 (15;16) | 46333 | Prostate cancer |
| D7 | 92 | F | F | 46235 | 9 | Non-DM | Cadaver | 2 (7;8) | 36090 | Natural death |

## Multiple-Testing Families
- Family A (per suite): all cell classes passing presence thresholds (present in at least the suite's minimum donor count, at least 20 total cells), BH across classes within the suite. Model CLR fraction ~ age_per_decade + suite covariates (log_lib excluded).
- Family B (per suite): all gene by cell-class tests in the suite, BH applied globally across every gene in every class of that suite (q_BH global); within-class BH also reported. Genes kept if CPM at least 1 in at least half of donors. Concordance: mean pairwise Pearson r of age effects across classes, 1,000 donor-label permutations.
- Family C (per suite): one `P1_Clean` composite test per cell class (118 human orthologs, 76 up, 42 down), BH across classes within the suite for the analytic p and separately for the 200-draw expression-matched random-set empirical p.

## Suite Overview
| Suite | Family A q<0.05 | Family B global q<0.05 | Concordance |
|---|---|---|---|
| All12_DiseaseAdj | 0 of 11 | 0 of 164089 (min q_BH = 0.3621) | r = +0.3303, one-sided permutation p = 0.0110 |
| All12_DiseaseAdj_ExprSex | 0 of 11 | 0 of 164089 (min q_BH = 0.3419) | r = +0.3349, one-sided permutation p = 0.0110 |
| All12_DiseaseLiveAdj | 0 of 11 | 0 of 164089 (min q_BH = 0.1697) | r = +0.2846, one-sided permutation p = 0.0020 |
| All12_WholeRetina_DiseaseAdj | not run | 0 of 18317 (min q_BH = 0.1231) | not computed |
| NonDM6_AgeOnly | 0 of 11 | 0 of 133884 (min q_BH = 0.2209) | r = +0.1729, one-sided permutation p = 0.1429 |
| Cadaver7_DiseaseAdj | 0 of 11 | 1 of 149018 (min q_BH = 0.0073) | r = +0.1005, one-sided permutation p = 0.1608 |

Suite definitions: `All12_DiseaseAdj` covariates is_dr, is_dm_no_dr, sex_male (Table 1), log_lib; `All12_DiseaseAdj_ExprSex` same with expression-verified sex; `All12_DiseaseLiveAdj` adds is_live; `All12_WholeRetina_DiseaseAdj` unpartitioned whole-retina pseudobulk with is_dr, is_dm_no_dr, is_live, sex_male, log_lib; `NonDM6_AgeOnly` the 6 non-diabetic donors with log_lib only; `Cadaver7_DiseaseAdj` the cadaver donors with is_dr, is_dm_no_dr, log_lib. Cell classes require at least 10 cells per donor.

## Family A: Primary Suite Detail (`All12_DiseaseAdj`)
| Cell Class | Donors Present | QC Cells | Median Fraction | Effect (CLR/dec) | p-value | q_BH (Family A) |
|---|---:|---:|---:|---:|---:|---:|
| HC | 12 | 2241 | 0.0069 | -0.1225 | 0.1790 | 0.9696 |
| Cone | 12 | 4840 | 0.0129 | +0.1398 | 0.3874 | 0.9696 |
| Astrocyte | 12 | 951 | 0.0024 | +0.0754 | 0.4879 | 0.9696 |
| MGC | 12 | 25267 | 0.1030 | -0.0541 | 0.7046 | 0.9696 |
| BC | 12 | 16619 | 0.0547 | -0.0295 | 0.7344 | 0.9696 |
| RGC | 12 | 1143 | 0.0026 | +0.0481 | 0.7765 | 0.9696 |
| TC | 12 | 983 | 0.0010 | -0.0348 | 0.8897 | 0.9696 |
| Microglia | 12 | 8283 | 0.0151 | -0.0119 | 0.9286 | 0.9696 |
| AC | 12 | 1425 | 0.0046 | -0.0052 | 0.9607 | 0.9696 |
| Rod | 12 | 218978 | 0.7983 | -0.0081 | 0.9625 | 0.9696 |
| RPE | 12 | 132 | 0.0002 | +0.0027 | 0.9696 | 0.9696 |

## Family B: All Suites
| Suite | Group | n_donors | Genes Tested | Residual df | Within q<0.05 | Global q<0.05 | Min p | Min q_BH Global |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| All12_DiseaseAdj | AC | 11 | 17484 | 5 | 0 | 0 | 6.24e-05 | 0.4285 |
| All12_DiseaseAdj | Astrocyte | 12 | 14383 | 6 | 0 | 0 | 7.30e-05 | 0.4285 |
| All12_DiseaseAdj | BC | 12 | 17522 | 6 | 0 | 0 | 4.07e-06 | 0.3621 |
| All12_DiseaseAdj | Cone | 12 | 15923 | 6 | 0 | 0 | 3.59e-05 | 0.4285 |
| All12_DiseaseAdj | HC | 11 | 16032 | 5 | 0 | 0 | 1.30e-05 | 0.4285 |
| All12_DiseaseAdj | MGC | 12 | 18564 | 6 | 0 | 0 | 4.41e-06 | 0.3621 |
| All12_DiseaseAdj | Microglia | 12 | 16444 | 6 | 0 | 0 | 6.42e-05 | 0.4285 |
| All12_DiseaseAdj | RGC | 9 | 16724 | 3 | 0 | 0 | 5.61e-05 | 0.4285 |
| All12_DiseaseAdj | Rod | 12 | 17105 | 6 | 0 | 0 | 2.38e-05 | 0.4285 |
| All12_DiseaseAdj | TC | 6 | 13908 | 2 | 0 | 0 | 7.04e-06 | 0.3851 |
| All12_DiseaseAdj_ExprSex | AC | 11 | 17484 | 5 | 0 | 0 | 1.98e-05 | 0.4117 |
| All12_DiseaseAdj_ExprSex | Astrocyte | 12 | 14383 | 6 | 0 | 0 | 1.99e-04 | 0.4117 |
| All12_DiseaseAdj_ExprSex | BC | 12 | 17522 | 6 | 1 | 0 | 2.08e-06 | 0.3419 |
| All12_DiseaseAdj_ExprSex | Cone | 12 | 15923 | 6 | 0 | 0 | 4.81e-05 | 0.4117 |
| All12_DiseaseAdj_ExprSex | HC | 11 | 16032 | 5 | 0 | 0 | 2.08e-05 | 0.4117 |
| All12_DiseaseAdj_ExprSex | MGC | 12 | 18564 | 6 | 0 | 0 | 5.55e-05 | 0.4117 |
| All12_DiseaseAdj_ExprSex | Microglia | 12 | 16444 | 6 | 0 | 0 | 3.32e-05 | 0.4117 |
| All12_DiseaseAdj_ExprSex | RGC | 9 | 16724 | 3 | 0 | 0 | 5.61e-05 | 0.4117 |
| All12_DiseaseAdj_ExprSex | Rod | 12 | 17105 | 6 | 0 | 0 | 7.39e-05 | 0.4117 |
| All12_DiseaseAdj_ExprSex | TC | 6 | 13908 | 2 | 0 | 0 | 7.04e-06 | 0.4117 |
| All12_DiseaseLiveAdj | AC | 11 | 17484 | 4 | 0 | 0 | 1.35e-04 | 0.4480 |
| All12_DiseaseLiveAdj | Astrocyte | 12 | 14383 | 5 | 0 | 0 | 2.67e-04 | 0.4480 |
| All12_DiseaseLiveAdj | BC | 12 | 17522 | 5 | 2 | 0 | 3.39e-06 | 0.1697 |
| All12_DiseaseLiveAdj | Cone | 12 | 15923 | 5 | 0 | 0 | 3.97e-06 | 0.1697 |
| All12_DiseaseLiveAdj | HC | 11 | 16032 | 4 | 0 | 0 | 1.86e-05 | 0.4480 |
| All12_DiseaseLiveAdj | MGC | 12 | 18564 | 5 | 0 | 0 | 3.20e-05 | 0.4480 |
| All12_DiseaseLiveAdj | Microglia | 12 | 16444 | 5 | 0 | 0 | 6.02e-05 | 0.4480 |
| All12_DiseaseLiveAdj | RGC | 9 | 16724 | 2 | 0 | 0 | 6.66e-05 | 0.4480 |
| All12_DiseaseLiveAdj | Rod | 12 | 17105 | 5 | 0 | 0 | 4.14e-06 | 0.1697 |
| All12_DiseaseLiveAdj | TC | 6 | 13908 | 2 | 0 | 0 | 7.04e-06 | 0.2310 |
| All12_WholeRetina_DiseaseAdj | All_Retina_QC | 12 | 18317 | 5 | 0 | 0 | 8.16e-06 | 0.1231 |
| NonDM6_AgeOnly | AC | 5 | 16676 | 2 | 0 | 0 | 1.32e-04 | 0.8725 |
| NonDM6_AgeOnly | Astrocyte | 6 | 15094 | 3 | 0 | 0 | 2.27e-04 | 0.8725 |
| NonDM6_AgeOnly | BC | 6 | 17655 | 3 | 0 | 0 | 9.81e-05 | 0.8725 |
| NonDM6_AgeOnly | Cone | 6 | 16195 | 3 | 0 | 0 | 9.22e-05 | 0.8725 |
| NonDM6_AgeOnly | HC | 5 | 15862 | 2 | 0 | 0 | 2.90e-05 | 0.4411 |
| NonDM6_AgeOnly | MGC | 6 | 18534 | 3 | 0 | 0 | 4.95e-06 | 0.3312 |
| NonDM6_AgeOnly | Microglia | 6 | 16586 | 3 | 1 | 0 | 1.65e-06 | 0.2209 |
| NonDM6_AgeOnly | Rod | 6 | 17282 | 3 | 0 | 0 | 2.36e-05 | 0.4411 |
| Cadaver7_DiseaseAdj | AC | 7 | 17783 | 2 | 0 | 0 | 6.77e-06 | 0.2621 |
| Cadaver7_DiseaseAdj | Astrocyte | 7 | 13266 | 2 | 2 | 1 | 4.92e-08 | 0.0073 |
| Cadaver7_DiseaseAdj | BC | 7 | 17438 | 2 | 0 | 0 | 2.02e-05 | 0.3744 |
| Cadaver7_DiseaseAdj | Cone | 7 | 15629 | 2 | 0 | 0 | 9.06e-05 | 0.3929 |
| Cadaver7_DiseaseAdj | HC | 7 | 16141 | 2 | 0 | 0 | 3.53e-05 | 0.3744 |
| Cadaver7_DiseaseAdj | MGC | 7 | 18593 | 2 | 0 | 0 | 1.51e-05 | 0.3744 |
| Cadaver7_DiseaseAdj | Microglia | 7 | 16166 | 2 | 0 | 0 | 1.17e-05 | 0.3480 |
| Cadaver7_DiseaseAdj | RGC | 7 | 16989 | 2 | 0 | 0 | 4.87e-05 | 0.3816 |
| Cadaver7_DiseaseAdj | Rod | 7 | 17013 | 2 | 0 | 0 | 7.04e-06 | 0.2621 |

## Family C: P1_Clean Mouse Signature, All Suites
| Suite | Group | n_donors | P1 Up/Down | Composite Effect | Analytic p | Empirical p (200 draws) | q_BH (Analytic) | q_BH (Empirical) |
|---|---|---:|:---:|---:|---:|---:|---:|---:|
| All12_DiseaseAdj | Rod | 12 | 62/27 | -0.0417 | 0.4867 | 0.0498 | 0.8111 | 0.4975 |
| All12_DiseaseAdj | AC | 11 | 56/29 | +0.0350 | 0.2592 | 0.2239 | 0.6752 | 0.5188 |
| All12_DiseaseAdj | Microglia | 12 | 63/29 | -0.0367 | 0.1653 | 0.3234 | 0.6752 | 0.5188 |
| All12_DiseaseAdj | HC | 11 | 54/25 | -0.0427 | 0.6296 | 0.3284 | 0.8995 | 0.5188 |
| All12_DiseaseAdj | Astrocyte | 12 | 57/26 | +0.0470 | 0.2008 | 0.3532 | 0.6752 | 0.5188 |
| All12_DiseaseAdj | RGC | 9 | 63/28 | -0.0388 | 0.4405 | 0.3582 | 0.8111 | 0.5188 |
| All12_DiseaseAdj | TC | 6 | 55/25 | +0.1411 | 0.2701 | 0.3632 | 0.6752 | 0.5188 |
| All12_DiseaseAdj | BC | 12 | 62/29 | -0.0111 | 0.8464 | 0.5871 | 0.9354 | 0.6910 |
| All12_DiseaseAdj | Cone | 12 | 54/23 | -0.0137 | 0.8550 | 0.6219 | 0.9354 | 0.6910 |
| All12_DiseaseAdj | MGC | 12 | 65/32 | -0.0020 | 0.9354 | 0.9552 | 0.9354 | 0.9552 |
| All12_DiseaseAdj_ExprSex | Rod | 12 | 62/27 | -0.0431 | 0.4637 | 0.0498 | 0.9273 | 0.4975 |
| All12_DiseaseAdj_ExprSex | Microglia | 12 | 63/29 | -0.0358 | 0.1702 | 0.3333 | 0.9003 | 0.7711 |
| All12_DiseaseAdj_ExprSex | RGC | 9 | 63/28 | -0.0388 | 0.4405 | 0.3582 | 0.9273 | 0.7711 |
| All12_DiseaseAdj_ExprSex | TC | 6 | 55/25 | +0.1411 | 0.2701 | 0.3632 | 0.9003 | 0.7711 |
| All12_DiseaseAdj_ExprSex | Astrocyte | 12 | 57/26 | +0.0442 | 0.1815 | 0.3980 | 0.9003 | 0.7711 |
| All12_DiseaseAdj_ExprSex | AC | 11 | 56/29 | +0.0166 | 0.7540 | 0.4925 | 0.9816 | 0.7711 |
| All12_DiseaseAdj_ExprSex | HC | 11 | 54/25 | -0.0304 | 0.6936 | 0.5771 | 0.9816 | 0.7711 |
| All12_DiseaseAdj_ExprSex | BC | 12 | 62/29 | -0.0108 | 0.8441 | 0.6169 | 0.9816 | 0.7711 |
| All12_DiseaseAdj_ExprSex | Cone | 12 | 54/23 | +0.0046 | 0.9518 | 0.8060 | 0.9816 | 0.8955 |
| All12_DiseaseAdj_ExprSex | MGC | 12 | 65/32 | -0.0004 | 0.9816 | 0.9154 | 0.9816 | 0.9154 |
| All12_DiseaseLiveAdj | Rod | 12 | 62/27 | -0.0441 | 0.5139 | 0.0498 | 0.7661 | 0.4367 |
| All12_DiseaseLiveAdj | RGC | 9 | 63/28 | -0.0557 | 0.4362 | 0.1493 | 0.7661 | 0.4367 |
| All12_DiseaseLiveAdj | HC | 11 | 54/25 | -0.0510 | 0.5757 | 0.1990 | 0.7661 | 0.4367 |
| All12_DiseaseLiveAdj | Cone | 12 | 54/23 | -0.0288 | 0.6895 | 0.2587 | 0.7661 | 0.4367 |
| All12_DiseaseLiveAdj | BC | 12 | 62/29 | -0.0225 | 0.6779 | 0.2985 | 0.7661 | 0.4367 |
| All12_DiseaseLiveAdj | Astrocyte | 12 | 57/26 | +0.0466 | 0.2485 | 0.3632 | 0.7661 | 0.4367 |
| All12_DiseaseLiveAdj | TC | 6 | 55/25 | +0.1411 | 0.2701 | 0.3632 | 0.7661 | 0.4367 |
| All12_DiseaseLiveAdj | Microglia | 12 | 63/29 | -0.0335 | 0.2204 | 0.3881 | 0.7661 | 0.4367 |
| All12_DiseaseLiveAdj | AC | 11 | 56/29 | +0.0189 | 0.4735 | 0.3930 | 0.7661 | 0.4367 |
| All12_DiseaseLiveAdj | MGC | 12 | 65/32 | -0.0042 | 0.8754 | 0.9602 | 0.8754 | 0.9602 |
| All12_WholeRetina_DiseaseAdj | All_Retina_QC | 12 | 66/31 | -0.0603 | 0.1715 | 0.0100 | 0.1715 | 0.0100 |
| NonDM6_AgeOnly | Rod | 6 | 63/28 | -0.0835 | 0.3194 | 0.0050 | 0.4259 | 0.0199 |
| NonDM6_AgeOnly | HC | 5 | 53/25 | -0.2060 | 0.2660 | 0.0050 | 0.4259 | 0.0199 |
| NonDM6_AgeOnly | Cone | 6 | 55/23 | -0.1067 | 0.2799 | 0.0100 | 0.4259 | 0.0265 |
| NonDM6_AgeOnly | BC | 6 | 62/30 | -0.0412 | 0.5839 | 0.0995 | 0.6673 | 0.1910 |
| NonDM6_AgeOnly | AC | 5 | 56/27 | -0.1100 | 0.2138 | 0.1194 | 0.4259 | 0.1910 |
| NonDM6_AgeOnly | Microglia | 6 | 64/30 | -0.0430 | 0.1276 | 0.2189 | 0.4259 | 0.2919 |
| NonDM6_AgeOnly | Astrocyte | 6 | 57/27 | +0.0625 | 0.2977 | 0.3134 | 0.4259 | 0.3582 |
| NonDM6_AgeOnly | MGC | 6 | 65/32 | +0.0011 | 0.9831 | 0.9154 | 0.9831 | 0.9154 |
| Cadaver7_DiseaseAdj | RGC | 7 | 65/28 | -0.0460 | 0.4321 | 0.1294 | 0.7490 | 0.4366 |
| Cadaver7_DiseaseAdj | Cone | 7 | 52/20 | -0.0614 | 0.7094 | 0.1841 | 0.7490 | 0.4366 |
| Cadaver7_DiseaseAdj | Rod | 7 | 61/27 | -0.0460 | 0.0645 | 0.1841 | 0.5077 | 0.4366 |
| Cadaver7_DiseaseAdj | Astrocyte | 7 | 57/24 | +0.0727 | 0.3194 | 0.1940 | 0.7185 | 0.4366 |
| Cadaver7_DiseaseAdj | HC | 7 | 52/24 | -0.0471 | 0.1128 | 0.3134 | 0.5077 | 0.5642 |
| Cadaver7_DiseaseAdj | MGC | 7 | 65/32 | +0.0144 | 0.2937 | 0.4925 | 0.7185 | 0.7264 |
| Cadaver7_DiseaseAdj | BC | 7 | 61/28 | -0.0118 | 0.6452 | 0.6219 | 0.7490 | 0.7264 |
| Cadaver7_DiseaseAdj | Microglia | 7 | 61/26 | +0.0094 | 0.7490 | 0.6617 | 0.7490 | 0.7264 |
| Cadaver7_DiseaseAdj | AC | 7 | 57/29 | -0.0143 | 0.6974 | 0.7264 | 0.7490 | 0.7264 |

## Key Limitations
1. Missing author cell labels and controlled raw reads; marker scoring approximates cell classes.
2. Confounding: 4 of 12 donors have DR, 2 have DM without DR, and all 5 living donors had eyes enucleated for ocular disease (phthisis bulbi, corneal perforation, ulcer, or staphyloma). Disease is adjusted as covariates, but with n = 12 the adjustment is coarse.
3. Sample size: residual df is at most 6 in the primary suite and lower in sensitivity suites and in classes with fewer donors.
4. The Table 1 sex entry for one donor disagrees with expression, and the Figshare cohort counts disagree with Table 1 (see above).

## Output Files
- `results/followup-organ-09-eye-hra014708-donors.csv`
- `results/followup-organ-09-eye-hra014708-samples-qc.csv`
- `results/followup-organ-09-eye-hra014708-access-check.csv`
- `results/followup-organ-09-eye-hra014708-fractions.csv`
- `results/followup-organ-09-eye-hra014708-gene-summary.csv`
- `results/followup-organ-09-eye-hra014708-top-genes.csv`
- `results/followup-organ-09-eye-hra014708-concordance.csv`
- `results/followup-organ-09-eye-hra014708-p1-signature.csv`
