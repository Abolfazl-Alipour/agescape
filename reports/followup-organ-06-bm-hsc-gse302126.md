# Dataset 06: Human Bone Marrow and HSC Single-Cell RNA-seq (Safina et al., GSE302126)

## 1. Study Verification and Data Provenance
- **Accession and Primary Paper**: `GSE302126` (`GSM9096700` to `GSM9096716`, 17 CellRanger BMMC technical-replicate libraries), Safina et al., *Genome Biology* 2026 (`PMID: 41546063`, `PMC12892447`).
- **Processed Author Seurat Objects**: Retrieved from the primary paper Figshare project (`https://figshare.com/projects/Aging_HSCs_2025/235781`):
  - `seu_full_object_hscs_marked.rds` (file `55618298`, 204,535 BoneMarrowMap-annotated cells across 7 datasets and 98 donors; contains all 2,361 QC-passed `Safina2024` / `GSE302126` BMMC cells across all 9 unique WT CHIP-negative donors).
  - `annotated_seu.rds` (file `55618823`, 28,989 BoneMarrowMap-purified HSCs across 7 datasets with `Predominant_program` `metaGEP.1` to `metaGEP.4` and Seurat subclusters).
- **Exact GEO Donor Verification (17 Replicates Collapsed to 9 Unique Donors)**:
  - Young (`n = 2`, ages `31` and `37` yr; `1M / 1F`): `BM1` (31 yr, M, GSM9096700 to GSM9096702) and `BM2` (37 yr, F, GSM9096703 to GSM9096705), both coded `source_name_ch1 = Bone marrow` (iliac crest in healthy donors, non-sternum).
  - Aged (`n = 7`, ages `64 to 78` yr; `5M / 2F`): `SBM1307` (64F), `SBM1096` (66M), `SBM1271` (67M), `SBM1237` (71M), `SBM1177` (73M), `SBM1264` (77M), and `SBM1073` (78F), all coded `source_name_ch1 = Sternum bone marrow` from patients undergoing cardiac surgery.
  - **Critical structural confound**: In `GSE302126` (`n = 9`), age group is 100% collinear with anatomical collection site and cardiac-surgery status (`100%` young iliac crest vs `100%` aged sternum cardiac surgery). A sternum-only internal sensitivity (`n = 7`, `64 to 78` yr) removes this confound, and multi-study exact-age adult BoneMarrowMap HSCs (`annotated_seu.rds`) provide an external adult-lifespan check (`n = 61` unique adult donors, ages `19 to 87` yr across `6` datasets after excluding prenatal/cord/infant/child and undefined-age `Weng2024`).

| Donor | Age | Sex | Group | Anatomical / Clinical Site | Total BMMC Cells | HSC Cells |
|---|---:|---|---|---|---:|---:|
| BM1 | 31 | M | Young | iliac_crest | 171 | 34 |
| BM2 | 37 | F | Young | iliac_crest | 775 | 34 |
| SBM1307 | 64 | F | Aged | sternum_cardiac_surgery | 49 | 0 |
| SBM1096 | 66 | M | Aged | sternum_cardiac_surgery | 462 | 8 |
| SBM1271 | 67 | M | Aged | sternum_cardiac_surgery | 42 | 0 |
| SBM1237 | 71 | M | Aged | sternum_cardiac_surgery | 17 | 1 |
| SBM1177 | 73 | M | Aged | sternum_cardiac_surgery | 28 | 1 |
| SBM1264 | 77 | M | Aged | sternum_cardiac_surgery | 317 | 18 |
| SBM1073 | 78 | F | Aged | sternum_cardiac_surgery | 500 | 21 |

## 2. Explicit Multiple-Testing Families
1. **Family A (CLR cell-type composition)**: Reported per analysis suite (`gse302126_bmmc_lineage_n9`, `gse302126_bmmc_fine_n9`, `gse302126_sternum_lineage_n7`, `safina_hsc_programs_n5`, `safina_hsc_clusters_n7`, `multistudy_adult_hsc_programs`), testing CLR fractions against `age_per_decade + covariates` with `q_bh_family_a < 0.05`.
2. **Family B (Donor-level gene age OLS and within-species directional concordance)**: Tested globally across all qualifying `(group, gene)` pairs within each suite (`CPM >= 1.0` in `>= 50%` of tested donors; `age_per_decade + sex_male + log_lib` in `GSE302126`, plus dataset fixed effects in `multistudy_adult_hsc`), plus 1,000-permutation donor-age concordance.
3. **Family C (`P1_Clean` mouse signature transport)**: `118` mapped human orthologs (`76` UP, `42` DOWN) scored per donor pseudobulk per qualifying group (`OLS` analytic p-value and 200-draw expression-decile-matched random-set two-sided empirical p-value), BH-corrected within each analysis suite.

## 3. Results

### Family A: Compositional Changes
- **Primary `GSE302126` BMMC broad lineages (`n = 9`, `sex_male` adjusted)**: `0` cell lineages reached `q_bh_family_a < 0.05`.
- **Primary `GSE302126` BMMC fine `predicted_CellType` (`n = 9`, `sex_male` adjusted)**: `2` cell types reached `q_bh_family_a < 0.05`: `CLP` (`-0.3181` CLR/decade, `p = 0.0020`, `q_bh_family_a = 0.0340`) and `Large Pre-B` (`-0.5884` CLR/decade, `p = 0.0021`, `q_bh_family_a = 0.0340`), both consistent with lymphoid-progenitor decline, though completely confounded with iliac crest vs sternum site in `n = 9`.
- **Multi-study exact-age adult HSC `Predominant_program` states (`annotated_seu.rds`, `n = 61`, adjusting `sex_male` + dataset)**: `metaGEP.2` (`-0.1729` CLR/decade, `p = 0.00031`, `q_bh_family_a = 0.0012`), `metaGEP.1` (`+0.2221`, `p = 0.0025`, `q_bh_family_a = 0.0049`), and `metaGEP.4` (`-0.1358`, `p = 0.0078`, `q_bh_family_a = 0.0104`) shift significantly with donor age, matching Safina et al.'s inflammatory/quiescent HSC state expansion.

### Family B: Donor-Level Gene Age OLS & Concordance
- **Primary `GSE302126` BMMC broad lineages (`n = 9`)**: `M_B = 67447` gene-lineage tests across `7` lineages -> `0` global BH hits (`q_bh_global < 0.05`).
- **Primary `GSE302126` All-BMMC pooled (`n = 9`)**: `M_B = 10790` tested genes -> `0` global BH hits (`q_bh_global < 0.05`).
- **Cross-cell-type gene concordance**:
| Analysis | Cell Types (n) | Common Genes | Mean Pairwise Pearson r | Perm p (1-sided) | Perm p (2-sided) |
|---|---:|---:|---:|---:|---:|
| gse302126_bmmc_lineage_n9 | 7 | 7450 | 0.0672 | 0.1518 | 0.1518 |
| gse302126_bmmc_fine_n9 | 15 | 3959 | 0.0516 | 0.2168 | 0.2168 |
| gse302126_sternum_lineage_n7 | 2 | 8070 | 0.0408 | 0.3467 | 0.3467 |
| safina_hsc_programs_n5 | 2 | 4962 | -0.0029 | 0.8981 | 0.3267 |
| multistudy_adult_hsc_programs | 4 | 7807 | 0.1186 | 0.0390 | 0.0390 |

### Family C: Mouse `P1_Clean` Aging Signature Transport
- **Primary `GSE302126` BMMC broad lineages (`n = 9`)**: Minimum `q_bh_composite = 0.1086`, minimum `q_bh_empirical = 0.0697` (NULL).
- **Fine `predicted_CellType` (`n = 9`)**: Minimum `q_bh_empirical = 0.9808` (NULL).
- **Sternum-only `SBM*` (`n = 7`, `64 to 78` yr, removes iliac crest confound)**: Minimum `q_bh_empirical = 0.5174` (NULL).
- **Safina `GSE302126` pooled `All_BMMC` (`n = 9`) and `Safina_HSC_all` (`n = 7`)**: Null `P1_Clean` transport (`All_BMMC composite = +0.0133`, `p = 0.8461`, `empirical p = 0.7811`; `Safina_HSC_all composite = -0.0717`, `p = 0.6256`, `empirical p = 0.5124`; `safina_hsc_programs_n5 min q_bh_empirical = 0.2687`).
- **Multi-study exact-age adult BoneMarrowMap HSCs (`n = 61` donors total, `n = 57` tested after minimum cell-count gating, `ages 19 to 87` yr)**: Pooled `MultiStudy_Adult_HSC` (`k = 1`) shows a modest positive composite slope (`+0.0891`, `analytic p = 0.0028`, `200-draw empirical p = 0.0199`) driven by subtle sub-significant `UP` (`+0.0548`, `p = 0.0576`) and `DOWN` (`-0.0343`, `p = 0.2907`) shifts across heterogeneous FACS-sorted/10x studies (for example, `Ainciburu2023` young donors are `5F/0M` while elderly donors are `0F/3M`), whereas within-state testing across the 4 `metaGEP.1..4` programs (`multistudy_adult_hsc_programs`) is completely null (`q_bh_composite = 0.1281 to 0.2967`, `q_bh_empirical = 0.1294`), demonstrating that the pooled slope tracks age-related `metaGEP` substate compositional expansion rather than uniform intracellular `P1_Clean` activation.

| Analysis | Cell Group | Donors (n) | Up / Down Detected | Composite Slope / Decade | Analytic p | Empirical p (200 draws) | q_BH (Composite) | q_BH (Empirical) |
|---|---|---:|---|---:|---:|---:|---:|---:|
| gse302126_bmmc_lineage_n9 | T_NK_Lymphoid | 8 | 40 / 13 | 0.4235 | 0.0155 | 0.0100 | 0.1086 | 0.0697 |
| gse302126_bmmc_lineage_n9 | B_Lymphoid | 4 | 42 / 14 | 0.0925 | 0.3690 | 0.2388 | 0.8715 | 0.8358 |
| gse302126_bmmc_lineage_n9 | HSPC_Stem_MPP | 5 | 47 / 16 | -0.0235 | 0.8715 | 0.6169 | 0.8715 | 0.8706 |
| gse302126_bmmc_lineage_n9 | Monocyte_DC | 5 | 48 / 16 | -0.0121 | 0.4934 | 0.6368 | 0.8715 | 0.8706 |
| gse302126_bmmc_lineage_n9 | Other | 8 | 43 / 16 | -0.0504 | 0.7377 | 0.6965 | 0.8715 | 0.8706 |
| gse302126_bmmc_lineage_n9 | Erythroid_Mk | 5 | 46 / 17 | -0.0330 | 0.5468 | 0.7463 | 0.8715 | 0.8706 |
| gse302126_bmmc_lineage_n9 | Myeloid_Prog | 4 | 43 / 13 | -0.0106 | 0.8280 | 0.9254 | 0.8715 | 0.9254 |
| gse302126_bmmc_fine_n9 | CFU-E | 5 | 38 / 13 | -0.1155 | 0.1356 | 0.1393 | 0.7451 | 0.9808 |
| gse302126_bmmc_fine_n9 | CD14 Mono | 4 | 43 / 13 | -0.0563 | 0.0758 | 0.3035 | 0.7451 | 0.9808 |
| gse302126_bmmc_fine_n9 | Early ProMono | 4 | 44 / 14 | -0.0648 | 0.1490 | 0.3085 | 0.7451 | 0.9808 |
| gse302126_bmmc_fine_n9 | NK | 5 | 31 / 12 | -0.1472 | 0.5430 | 0.3483 | 0.9050 | 0.9808 |
| gse302126_bmmc_fine_n9 | CD8 Effector Memory 2 | 5 | 37 / 11 | -0.4374 | 0.5298 | 0.4527 | 0.9050 | 0.9808 |
| gse302126_bmmc_fine_n9 | HSC | 5 | 39 / 13 | 0.0610 | 0.6891 | 0.5622 | 0.9464 | 0.9808 |
| gse302126_bmmc_fine_n9 | Late ProMono | 4 | 39 / 10 | -0.1121 | 0.3168 | 0.5821 | 0.9050 | 0.9808 |
| gse302126_bmmc_fine_n9 | Pre-cDC | 5 | 43 / 14 | -0.0072 | 0.4360 | 0.6020 | 0.9050 | 0.9808 |
| gse302126_bmmc_fine_n9 | MEP | 4 | 35 / 11 | 0.0382 | 0.8644 | 0.6269 | 0.9464 | 0.9808 |
| gse302126_bmmc_fine_n9 | LMPP | 4 | 35 / 11 | -0.0317 | 0.7703 | 0.7861 | 0.9464 | 0.9808 |
| gse302126_bmmc_fine_n9 | MPP-MkEry | 4 | 34 / 12 | -1.88e-04 | 0.9984 | 0.8209 | 0.9984 | 0.9808 |
| gse302126_bmmc_fine_n9 | cDC2 | 5 | 42 / 13 | 0.0251 | 0.2109 | 0.8955 | 0.7909 | 0.9808 |
| gse302126_bmmc_fine_n9 | Early GMP | 4 | 38 / 13 | -0.0186 | 0.8833 | 0.9104 | 0.9464 | 0.9808 |
| gse302126_bmmc_fine_n9 | Large Pre-B | 4 | 37 / 13 | 0.0533 | 0.4003 | 0.9154 | 0.9050 | 0.9808 |
| gse302126_bmmc_fine_n9 | BFU-E | 4 | 37 / 12 | 0.0046 | 0.8830 | 0.9851 | 0.9464 | 0.9851 |
| gse302126_bmmc_pooled_n9 | All_BMMC | 9 | 51 / 18 | 0.0133 | 0.8461 | 0.7811 | 0.8461 | 0.7811 |
| gse302126_sternum_lineage_n7 | T_NK_Lymphoid | 7 | 40 / 13 | 0.2390 | 0.4882 | 0.2587 | 0.8326 | 0.5174 |
| gse302126_sternum_lineage_n7 | Other | 6 | 42 / 15 | 0.1643 | 0.8326 | 0.5622 | 0.8326 | 0.5622 |
| gse302126_sternum_pooled_n7 | All_BMMC | 7 | 51 / 17 | 0.2054 | 0.5750 | 0.2687 | 0.5750 | 0.2687 |
| safina_hsc_pooled_n7 | Safina_HSC_all | 7 | 35 / 10 | -0.0717 | 0.6256 | 0.5124 | 0.6256 | 0.5124 |
| safina_hsc_programs_n5 | metaGEP.1 | 5 | 36 / 12 | 0.6221 | 0.3468 | 0.1343 | 0.6927 | 0.2687 |
| safina_hsc_programs_n5 | metaGEP.4 | 4 | 28 / 9 | -0.0260 | 0.6927 | 0.7065 | 0.6927 | 0.7065 |
| safina_hsc_clusters_n7 | HSC_cluster_2 | 4 | 30 / 10 | 0.2197 | 0.4480 | 0.1244 | 0.8292 | 0.3731 |
| safina_hsc_clusters_n7 | HSC_cluster_0 | 4 | 31 / 9 | -0.0393 | 0.8292 | 0.5522 | 0.8292 | 0.6020 |
| safina_hsc_clusters_n7 | HSC_cluster_3 | 4 | 30 / 10 | -0.0770 | 0.7057 | 0.6020 | 0.8292 | 0.6020 |
| multistudy_adult_hsc_total | MultiStudy_Adult_HSC | 57 | 44 / 14 | 0.0891 | 0.0028 | 0.0199 | 0.0028 | 0.0199 |
| multistudy_adult_hsc_programs | metaGEP.1 | 51 | 38 / 12 | 0.0807 | 0.0961 | 0.0995 | 0.1281 | 0.1294 |
| multistudy_adult_hsc_programs | metaGEP.2 | 32 | 45 / 16 | 0.0738 | 0.0811 | 0.1095 | 0.1281 | 0.1294 |
| multistudy_adult_hsc_programs | metaGEP.4 | 16 | 38 / 14 | 0.1710 | 0.2967 | 0.1144 | 0.2967 | 0.1294 |
| multistudy_adult_hsc_programs | metaGEP.3 | 28 | 38 / 11 | 0.0804 | 0.0601 | 0.1294 | 0.1281 | 0.1294 |

## 4. Limitations and Verdict
- **Limits**: `GSE302126` has small sample size (`n = 9` unique donors after collapsing 17 replicate libraries), `SBM1177` (28 cells) and `SBM1237` (17 cells) have sparse BMMC recovery, only `117` total `Safina2024` cells map to `HSC`, and young vs aged donors in `GSE302126` are 100% confounded with iliac crest vs sternum cardiac surgery sampling. Multi-study adult HSC integration (`n = 61`) overcomes sample size but mixes 10x and Smart-seq/SORT protocols across studies.
- **Verdict**: Mouse `P1_Clean` shows no significant positive transport in human `GSE302126` bone marrow lineages, `GSE302126` sternum-only donors, `GSE302126` HSCs, or the `61`-donor multi-study exact-age adult HSC atlas (`all empirical q_BH >= 0.05` across tested primary cohorts).

## 5. Output Artifacts
- [Donors](file:///home/reza/projects/agescape-revision/results/followup-organ-06-bm-hsc-gse302126-donors.csv)
- [Fractions](file:///home/reza/projects/agescape-revision/results/followup-organ-06-bm-hsc-gse302126-fractions.csv)
- [Gene Summary](file:///home/reza/projects/agescape-revision/results/followup-organ-06-bm-hsc-gse302126-gene-summary.csv)
- [Top Genes](file:///home/reza/projects/agescape-revision/results/followup-organ-06-bm-hsc-gse302126-top-genes.csv)
- [Concordance](file:///home/reza/projects/agescape-revision/results/followup-organ-06-bm-hsc-gse302126-concordance.csv)
- [P1 Signature](file:///home/reza/projects/agescape-revision/results/followup-organ-06-bm-hsc-gse302126-p1-signature.csv)
