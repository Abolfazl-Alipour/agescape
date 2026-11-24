# Matched mouse-human P1m comparison

## Headline results

Pooled statistics weight every pair-gene row equally. P values are two-sided permutation values centered on the corresponding null mean, with the plus-one correction.

| arm | pairs | pairs with genes | pair-gene rows | sign agreement | gene-label p (sign) | donor-age p (sign) | Spearman | gene-label p (Spearman) | donor-age p (Spearman) | Pearson | gene-label p (Pearson) | donor-age p (Pearson) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| floored | 27 | 21 | 294 | 0.391156 | 0.025987 | 0.0549451 | 0.042512 | 0.492254 | 0.746254 | 0.0580543 | 0.321339 | 0.652348 |
| unfloored | 28 | 22 | 299 | 0.38796 | 0.0164918 | 0.0509491 | 0.0270615 | 0.616692 | 0.84016 | 0.0498349 | 0.365817 | 0.634366 |

Matching covered 16 of 36 mouse units (44.4%), yielding 28 matched pairs across cohorts. 20 mouse units were unmatched.

**Pooled verdict:** sign agreement is below 0.5 in both arms. Its gene-label permutation p values are below 0.05, but this does not mean agreement exceeds 0.5; the donor-age permutation p values are 0.0549451 (floored) and 0.0509491 (unfloored). Spearman and Pearson correlations are near zero, with no p value below 0.05 under either null in either arm. Thus the run does not show robust positive gene-effect concordance; the sign result is null-sensitive and is not evidence of conserved directionality.

## Design

The committed P1m derivation was read, not re-derived: 1,914 cleaned genes (1,119 up and 795 down), with 1,724 human-mapped genes. Each mouse unit contributes exactly the global-list genes whose `units_list` contains that unit; mouse effects are that unit's derivation log2FC for 30m versus 18m.

Human units were matched by tissue and exact cell-type identity after lowercasing and removing punctuation and whitespace. Marrow matches only to blood; no synonyms were introduced. Kidney Living Donor and SCR v15, lung Integrated, CellRef and ImmuneAging BAL, Yazar blood, PBMC Thomson/Zhu and Perez blood were audited. Normal Adult and SN v15 kidney sources are absent and were not fabricated.

Human expression is donor pseudobulk from donors aged 30 years or older: raw counts are summed within donor and matched cell type, normalized to counts per 10,000 within that donor-cell-type pseudobulk, transformed with log1p, and regressed on continuous age with OLS plus sex where sex is variable and identifiable. No depth, pool, batch, or detection covariate/filter was added. Duplicate feature symbols are summed before normalization. A gene is tested when its human ortholog is present and its slope is finite; zero slopes count as sign disagreement. One-donor and no-age-variation units have slope NA and contribute no pooled rows. Two-donor fits are age-only differences with no variance estimate, no slope SE/p, and are explicitly flagged.

The floored arm retains the 36 eligible mouse derivation units and, among their matched human pairs, requires at least 3 donors aged 30+ and at least 50 cells in total from those donors. This is the human analog of the mouse two-group floor: 'per side' applies to the mouse 18m/30m design, whereas the human fit is one continuous-age regression over donors aged 30+. The unfloored arm includes every matched pair regardless of donor or cell count. Removing this floor does not change the P1m gene list or add mouse units.

Gene-label nulls use 2,000 shuffles, within pair for pair-level p values and across all pair-gene rows for pooled p values. Donor-age nulls use 1,000 permutations within cohort and sex strata, sharing each donor's permuted age across all matched cell types from that donor, refitting every human slope, and recomputing pair-level and pooled statistics. Seeds are 20261124 for age permutations and deterministic pair-specific offsets from 20261125 for gene permutations.

## Match audit (all 36 mouse units)

| mouse unit | normalized name | decision | present matches | match decision / reason |
| --- | --- | --- | --- | --- |
| Lung::classical monocyte | classicalmonocyte | matched | 3 | Lung_Integrated::classical monocyte; Lung_CellRef::classical monocyte; ImmuneAging_BAL::classical monocyte |
| Marrow::granulocyte | granulocyte | unmatched_no_cell_type | 0 | no exactly normalized same-named cell type in same-tissue human roster |
| Marrow::granulocytopoietic cell | granulocytopoieticcell | unmatched_no_cell_type | 0 | no exactly normalized same-named cell type in same-tissue human roster |
| Kidney::epithelial cell of proximal tubule | epithelialcellofproximaltubule | matched | 2 | Kidney_LivingDonor::epithelial cell of proximal tubule; Kidney_SCR_v15::epithelial cell of proximal tubule; additional absent exact-name sources: Kidney_NormalAdult, Kidney_SN_v15 |
| Lung::intermediate monocyte | intermediatemonocyte | unmatched_no_cell_type | 0 | no exactly normalized same-named cell type in same-tissue human roster |
| Marrow::promonocyte | promonocyte | unmatched_no_cell_type | 0 | no exactly normalized same-named cell type in same-tissue human roster |
| Marrow::erythroblast | erythroblast | unmatched_no_cell_type | 0 | no exactly normalized same-named cell type in same-tissue human roster |
| Marrow::proerythroblast | proerythroblast | unmatched_no_cell_type | 0 | no exactly normalized same-named cell type in same-tissue human roster |
| Kidney::kidney proximal convoluted tubule epithelial cell | kidneyproximalconvolutedtubuleepithelialcell | unmatched_no_cell_type | 0 | no exactly normalized same-named cell type in same-tissue human roster |
| Marrow::megakaryocyte-erythroid progenitor cell | megakaryocyteerythroidprogenitorcell | unmatched_no_cell_type | 0 | no exactly normalized same-named cell type in same-tissue human roster |
| Kidney::macrophage | macrophage | unmatched_no_cell_type | 0 | no exactly normalized same-named cell type in same-tissue human roster |
| Marrow::hematopoietic precursor cell | hematopoieticprecursorcell | matched | 1 | Yazar_Blood::hematopoietic precursor cell |
| Lung::bronchial smooth muscle cell | bronchialsmoothmusclecell | unmatched_no_cell_type | 0 | no exactly normalized same-named cell type in same-tissue human roster |
| Marrow::monocyte | monocyte | unmatched_no_cell_type | 0 | no exactly normalized same-named cell type in same-tissue human roster |
| Marrow::naive T cell | naivetcell | unmatched_no_cell_type | 0 | no exactly normalized same-named cell type in same-tissue human roster |
| Kidney::kidney collecting duct principal cell | kidneycollectingductprincipalcell | matched | 1 | Kidney_SCR_v15::kidney collecting duct principal cell; additional absent exact-name sources: Kidney_NormalAdult, Kidney_SN_v15 |
| Kidney::kidney loop of Henle thick ascending limb epithelial cell | kidneyloopofhenlethickascendinglimbepithelialcell | matched | 1 | Kidney_SCR_v15::kidney loop of Henle thick ascending limb epithelial cell; additional absent exact-name sources: Kidney_NormalAdult, Kidney_SN_v15 |
| Lung::B cell | bcell | matched | 3 | Lung_Integrated::B cell; Lung_CellRef::B cell; ImmuneAging_BAL::B cell |
| Marrow::macrophage | macrophage | unmatched_no_cell_type | 0 | no exactly normalized same-named cell type in same-tissue human roster |
| Lung::alveolar macrophage | alveolarmacrophage | matched | 2 | Lung_Integrated::alveolar macrophage; Lung_CellRef::alveolar macrophage |
| Lung::non-classical monocyte | nonclassicalmonocyte | matched | 3 | Lung_Integrated::non-classical monocyte; Lung_CellRef::non-classical monocyte; ImmuneAging_BAL::non-classical monocyte |
| Kidney::fenestrated cell | fenestratedcell | unmatched_no_cell_type | 0 | no exactly normalized same-named cell type in same-tissue human roster |
| Kidney::T cell | tcell | matched | 2 | Kidney_LivingDonor::T cell; Kidney_SCR_v15::T cell; additional absent exact-name sources: Kidney_NormalAdult, Kidney_SN_v15 |
| Marrow::late pro-B cell | lateprobcell | unmatched_no_cell_type | 0 | no exactly normalized same-named cell type in same-tissue human roster |
| Lung::CD8-positive, alpha-beta T cell | cd8positivealphabetatcell | matched | 1 | Lung_Integrated::CD8-positive, alpha-beta T cell |
| Lung::mature NK T cell | maturenktcell | unmatched_no_cell_type | 0 | no exactly normalized same-named cell type in same-tissue human roster |
| Marrow::precursor B cell | precursorbcell | unmatched_no_cell_type | 0 | no exactly normalized same-named cell type in same-tissue human roster |
| Kidney::kidney distal convoluted tubule epithelial cell | kidneydistalconvolutedtubuleepithelialcell | matched | 2 | Kidney_LivingDonor::kidney distal convoluted tubule epithelial cell; Kidney_SCR_v15::kidney distal convoluted tubule epithelial cell; additional absent exact-name sources: Kidney_NormalAdult, Kidney_SN_v15 |
| Lung::NK cell | nkcell | unmatched_no_cell_type | 0 | no exactly normalized same-named cell type in same-tissue human roster |
| Kidney::lymphocyte | lymphocyte | unmatched_no_cell_type | 0 | no exactly normalized same-named cell type in same-tissue human roster |
| Marrow::naive B cell | naivebcell | matched | 1 | Yazar_Blood::naive B cell |
| Lung::CD4-positive, alpha-beta T cell | cd4positivealphabetatcell | matched | 1 | Lung_Integrated::CD4-positive, alpha-beta T cell |
| Kidney::kidney cortex artery cell | kidneycortexarterycell | unmatched_no_cell_type | 0 | no exactly normalized same-named cell type in same-tissue human roster |
| Kidney::B cell | bcell | matched | 2 | Kidney_LivingDonor::B cell; Kidney_SCR_v15::B cell; additional absent exact-name sources: Kidney_NormalAdult, Kidney_SN_v15 |
| Lung::T cell | tcell | matched | 1 | Lung_Integrated::T cell |
| Lung::vein endothelial cell | veinendothelialcell | matched | 2 | Lung_Integrated::vein endothelial cell; Lung_CellRef::vein endothelial cell |

Unmatched mouse units: Marrow::granulocyte; Marrow::granulocytopoietic cell; Lung::intermediate monocyte; Marrow::promonocyte; Marrow::erythroblast; Marrow::proerythroblast; Kidney::kidney proximal convoluted tubule epithelial cell; Marrow::megakaryocyte-erythroid progenitor cell; Kidney::macrophage; Lung::bronchial smooth muscle cell; Marrow::monocyte; Marrow::naive T cell; Marrow::macrophage; Kidney::fenestrated cell; Marrow::late pro-B cell; Lung::mature NK T cell; Marrow::precursor B cell; Lung::NK cell; Kidney::lymphocyte; Kidney::kidney cortex artery cell.

The match-rate assertion passed: matched plus unmatched mouse units = 16 + 20 = 36.

## Per-pair results

| mouse_unit | mouse_tissue | human_cohort | human_unit | arm | n_mouse_genes | n_human_genes_tested | sign_agreement_fraction | spearman | pearson | p_sign_gene_perm | p_spearman_gene_perm | p_sign_age_perm | p_spearman_age_perm | n_donors_human_30plus | n_cells_human_30plus | below_old_floor | notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Kidney::B cell | Kidney | Kidney_LivingDonor | B cell | floored | 4 | 4 | 0.5 | -0.8 | -0.2078 | 1 | 0.35832 | 1 | 0.20779 | 11 | 68 | false | ortholog_mapped=4;ortholog_collisions=0;human_genes_present=4;sex_adjusted=true;below_old_floor=false |
| Kidney::B cell | Kidney | Kidney_SCR_v15 | B cell | floored | 4 | 4 | 0.5 | 0.6 | 0.85298 | 1 | 0.38731 | 1 | 0.50849 | 13 | 131 | false | ortholog_mapped=4;ortholog_collisions=0;human_genes_present=4;sex_adjusted=true;below_old_floor=false |
| Kidney::T cell | Kidney | Kidney_LivingDonor | T cell | floored | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | 12 | 938 | false | ortholog_mapped=0;ortholog_collisions=0;human_genes_present=0;sex_adjusted=true;below_old_floor=false |
| Kidney::T cell | Kidney | Kidney_SCR_v15 | T cell | floored | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | 18 | 617 | false | ortholog_mapped=0;ortholog_collisions=0;human_genes_present=0;sex_adjusted=true;below_old_floor=false |
| Kidney::epithelial cell of proximal tubule | Kidney | Kidney_LivingDonor | epithelial cell of proximal tubule | floored | 16 | 15 | 0.46667 | -0.4 | -0.43453 | 0.61719 | 0.14193 | 0.78422 | 0.22378 | 19 | 20983 | false | ortholog_mapped=15;ortholog_collisions=0;human_genes_present=15;sex_adjusted=true;below_old_floor=false |
| Kidney::epithelial cell of proximal tubule | Kidney | Kidney_SCR_v15 | epithelial cell of proximal tubule | floored | 16 | 14 | 0.5 | 0.1033 | 0.068766 | 1 | 0.71564 | 1 | 0.73626 | 18 | 6429 | false | ortholog_mapped=15;ortholog_collisions=0;human_genes_present=14;sex_adjusted=true;below_old_floor=false |
| Kidney::kidney collecting duct principal cell | Kidney | Kidney_SCR_v15 | kidney collecting duct principal cell | floored | 22 | 17 | 0.47059 | 0.10049 | 0.12254 | 0.69615 | 0.73013 | 1 | 0.87712 | 18 | 2577 | false | ortholog_mapped=18;ortholog_collisions=0;human_genes_present=17;sex_adjusted=true;below_old_floor=false |
| Kidney::kidney distal convoluted tubule epithelial cell | Kidney | Kidney_LivingDonor | kidney distal convoluted tubule epithelial cell | floored | 9 | 7 | 0.14286 | 0.035714 | 0.045116 | 1 | 0.93453 | 0.11588 | 0.85115 | 18 | 109 | false | ortholog_mapped=7;ortholog_collisions=0;human_genes_present=7;sex_adjusted=true;below_old_floor=false |
| Kidney::kidney distal convoluted tubule epithelial cell | Kidney | Kidney_SCR_v15 | kidney distal convoluted tubule epithelial cell | floored | 9 | 6 | 0.83333 | -0.25714 | -0.14021 | 1 | 0.65917 | 0.2977 | 0.4016 | 17 | 2015 | false | ortholog_mapped=7;ortholog_collisions=0;human_genes_present=6;sex_adjusted=true;below_old_floor=false |
| Kidney::kidney loop of Henle thick ascending limb epithelial cell | Kidney | Kidney_SCR_v15 | kidney loop of Henle thick ascending limb epithelial cell | floored | 15 | 11 | 0.54545 | 0.5 | 0.62694 | 1 | 0.12744 | 0.8971 | 0.28272 | 18 | 4042 | false | ortholog_mapped=12;ortholog_collisions=0;human_genes_present=11;sex_adjusted=true;below_old_floor=false |
| Lung::B cell | Lung | ImmuneAging_BAL | B cell | floored | 2 | 1 | 1 | NA | NA | 1 | NA | 0.44855 | NA | 7 | 87 | false | ortholog_mapped=1;ortholog_collisions=0;human_genes_present=1;sex_adjusted=true;below_old_floor=false |
| Lung::B cell | Lung | Lung_CellRef | B cell | floored | 2 | 1 | 1 | NA | NA | 1 | NA | 1 | NA | 9 | 72 | false | ortholog_mapped=1;ortholog_collisions=0;human_genes_present=1;sex_adjusted=true;below_old_floor=false |
| Lung::B cell | Lung | Lung_Integrated | B cell | floored | 2 | 1 | 0 | NA | NA | 1 | NA | 1 | NA | 62 | 2341 | false | ortholog_mapped=1;ortholog_collisions=0;human_genes_present=1;sex_adjusted=true;below_old_floor=false |
| Lung::CD4-positive, alpha-beta T cell | Lung | Lung_Integrated | CD4-positive, alpha-beta T cell | floored | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | 73 | 10815 | false | ortholog_mapped=0;ortholog_collisions=0;human_genes_present=0;sex_adjusted=true;below_old_floor=false |
| Lung::CD8-positive, alpha-beta T cell | Lung | Lung_Integrated | CD8-positive, alpha-beta T cell | floored | 1 | 1 | 1 | NA | NA | 1 | NA | 1 | NA | 74 | 15164 | false | ortholog_mapped=1;ortholog_collisions=0;human_genes_present=1;sex_adjusted=true;below_old_floor=false |
| Lung::T cell | Lung | Lung_Integrated | T cell | floored | 2 | 1 | 1 | NA | NA | 1 | NA | 0.49151 | NA | 35 | 379 | false | ortholog_mapped=1;ortholog_collisions=0;human_genes_present=1;sex_adjusted=true;below_old_floor=false |
| Lung::alveolar macrophage | Lung | Lung_CellRef | alveolar macrophage | floored | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | 13 | 42464 | false | ortholog_mapped=0;ortholog_collisions=0;human_genes_present=0;sex_adjusted=true;below_old_floor=false |
| Lung::alveolar macrophage | Lung | Lung_Integrated | alveolar macrophage | floored | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | 71 | 50125 | false | ortholog_mapped=0;ortholog_collisions=0;human_genes_present=0;sex_adjusted=true;below_old_floor=false |
| Lung::classical monocyte | Lung | ImmuneAging_BAL | classical monocyte | floored | 7 | 4 | 1 | 0.8 | 0.94774 | 0.25437 | 0.35182 | 0.077922 | 0.22877 | 8 | 45354 | false | ortholog_mapped=4;ortholog_collisions=0;human_genes_present=4;sex_adjusted=true;below_old_floor=false |
| Lung::classical monocyte | Lung | Lung_CellRef | classical monocyte | floored | 7 | 4 | 0.25 | 0.4 | 0.55247 | 1 | 0.58821 | 1 | 0.65335 | 9 | 99 | false | ortholog_mapped=4;ortholog_collisions=0;human_genes_present=4;sex_adjusted=true;below_old_floor=false |
| Lung::classical monocyte | Lung | Lung_Integrated | classical monocyte | floored | 7 | 4 | 0.25 | -0.4 | -0.76276 | 1 | 0.57871 | 0.33766 | 0.82517 | 74 | 12991 | false | ortholog_mapped=4;ortholog_collisions=0;human_genes_present=4;sex_adjusted=true;below_old_floor=false |
| Lung::non-classical monocyte | Lung | ImmuneAging_BAL | non-classical monocyte | floored | 7 | 5 | 0 | -0.7 | -0.92229 | 0.093453 | 0.1984 | 0.13187 | 0.30669 | 8 | 7309 | false | ortholog_mapped=5;ortholog_collisions=0;human_genes_present=5;sex_adjusted=true;below_old_floor=false |
| Lung::non-classical monocyte | Lung | Lung_Integrated | non-classical monocyte | floored | 7 | 5 | 0.6 | -0.1 | 0.30437 | 1 | 0.95502 | 0.7023 | 0.9001 | 65 | 6931 | false | ortholog_mapped=5;ortholog_collisions=0;human_genes_present=5;sex_adjusted=true;below_old_floor=false |
| Lung::vein endothelial cell | Lung | Lung_CellRef | vein endothelial cell | floored | 1 | 1 | 0 | NA | NA | 1 | NA | 0.49251 | NA | 11 | 1635 | false | ortholog_mapped=1;ortholog_collisions=0;human_genes_present=1;sex_adjusted=true;below_old_floor=false |
| Lung::vein endothelial cell | Lung | Lung_Integrated | vein endothelial cell | floored | 1 | 1 | 0 | NA | NA | 1 | NA | 1 | NA | 63 | 5755 | false | ortholog_mapped=1;ortholog_collisions=0;human_genes_present=1;sex_adjusted=true;below_old_floor=false |
| Marrow::hematopoietic precursor cell | Marrow | Yazar_Blood | hematopoietic precursor cell | floored | 221 | 187 | 0.34225 | 0.031554 | -0.046851 | 0.009995 | 0.68716 | 0.017982 | 0.7013 | 670 | 1718 | false | ortholog_mapped=188;ortholog_collisions=0;human_genes_present=187;sex_adjusted=true;below_old_floor=false |
| Marrow::naive B cell | Marrow | Yazar_Blood | naive B cell | floored | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | 931 | 61914 | false | ortholog_mapped=0;ortholog_collisions=0;human_genes_present=0;sex_adjusted=true;below_old_floor=false |
| Kidney::B cell | Kidney | Kidney_LivingDonor | B cell | unfloored | 4 | 4 | 0.5 | -0.8 | -0.2078 | 1 | 0.35832 | 1 | 0.20779 | 11 | 68 | false | ortholog_mapped=4;ortholog_collisions=0;human_genes_present=4;sex_adjusted=true;below_old_floor=false |
| Kidney::B cell | Kidney | Kidney_SCR_v15 | B cell | unfloored | 4 | 4 | 0.5 | 0.6 | 0.85298 | 1 | 0.38731 | 1 | 0.50849 | 13 | 131 | false | ortholog_mapped=4;ortholog_collisions=0;human_genes_present=4;sex_adjusted=true;below_old_floor=false |
| Kidney::T cell | Kidney | Kidney_LivingDonor | T cell | unfloored | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | 12 | 938 | false | ortholog_mapped=0;ortholog_collisions=0;human_genes_present=0;sex_adjusted=true;below_old_floor=false |
| Kidney::T cell | Kidney | Kidney_SCR_v15 | T cell | unfloored | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | 18 | 617 | false | ortholog_mapped=0;ortholog_collisions=0;human_genes_present=0;sex_adjusted=true;below_old_floor=false |
| Kidney::epithelial cell of proximal tubule | Kidney | Kidney_LivingDonor | epithelial cell of proximal tubule | unfloored | 16 | 15 | 0.46667 | -0.4 | -0.43453 | 0.61719 | 0.14193 | 0.78422 | 0.22378 | 19 | 20983 | false | ortholog_mapped=15;ortholog_collisions=0;human_genes_present=15;sex_adjusted=true;below_old_floor=false |
| Kidney::epithelial cell of proximal tubule | Kidney | Kidney_SCR_v15 | epithelial cell of proximal tubule | unfloored | 16 | 14 | 0.5 | 0.1033 | 0.068766 | 1 | 0.71564 | 1 | 0.73626 | 18 | 6429 | false | ortholog_mapped=15;ortholog_collisions=0;human_genes_present=14;sex_adjusted=true;below_old_floor=false |
| Kidney::kidney collecting duct principal cell | Kidney | Kidney_SCR_v15 | kidney collecting duct principal cell | unfloored | 22 | 17 | 0.47059 | 0.10049 | 0.12254 | 0.69615 | 0.73013 | 1 | 0.87712 | 18 | 2577 | false | ortholog_mapped=18;ortholog_collisions=0;human_genes_present=17;sex_adjusted=true;below_old_floor=false |
| Kidney::kidney distal convoluted tubule epithelial cell | Kidney | Kidney_LivingDonor | kidney distal convoluted tubule epithelial cell | unfloored | 9 | 7 | 0.14286 | 0.035714 | 0.045116 | 1 | 0.93453 | 0.11588 | 0.85115 | 18 | 109 | false | ortholog_mapped=7;ortholog_collisions=0;human_genes_present=7;sex_adjusted=true;below_old_floor=false |
| Kidney::kidney distal convoluted tubule epithelial cell | Kidney | Kidney_SCR_v15 | kidney distal convoluted tubule epithelial cell | unfloored | 9 | 6 | 0.83333 | -0.25714 | -0.14021 | 1 | 0.65917 | 0.2977 | 0.4016 | 17 | 2015 | false | ortholog_mapped=7;ortholog_collisions=0;human_genes_present=6;sex_adjusted=true;below_old_floor=false |
| Kidney::kidney loop of Henle thick ascending limb epithelial cell | Kidney | Kidney_SCR_v15 | kidney loop of Henle thick ascending limb epithelial cell | unfloored | 15 | 11 | 0.54545 | 0.5 | 0.62694 | 1 | 0.12744 | 0.8971 | 0.28272 | 18 | 4042 | false | ortholog_mapped=12;ortholog_collisions=0;human_genes_present=11;sex_adjusted=true;below_old_floor=false |
| Lung::B cell | Lung | ImmuneAging_BAL | B cell | unfloored | 2 | 1 | 1 | NA | NA | 1 | NA | 0.44855 | NA | 7 | 87 | false | ortholog_mapped=1;ortholog_collisions=0;human_genes_present=1;sex_adjusted=true;below_old_floor=false |
| Lung::B cell | Lung | Lung_CellRef | B cell | unfloored | 2 | 1 | 1 | NA | NA | 1 | NA | 1 | NA | 9 | 72 | false | ortholog_mapped=1;ortholog_collisions=0;human_genes_present=1;sex_adjusted=true;below_old_floor=false |
| Lung::B cell | Lung | Lung_Integrated | B cell | unfloored | 2 | 1 | 0 | NA | NA | 1 | NA | 1 | NA | 62 | 2341 | false | ortholog_mapped=1;ortholog_collisions=0;human_genes_present=1;sex_adjusted=true;below_old_floor=false |
| Lung::CD4-positive, alpha-beta T cell | Lung | Lung_Integrated | CD4-positive, alpha-beta T cell | unfloored | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | 73 | 10815 | false | ortholog_mapped=0;ortholog_collisions=0;human_genes_present=0;sex_adjusted=true;below_old_floor=false |
| Lung::CD8-positive, alpha-beta T cell | Lung | Lung_Integrated | CD8-positive, alpha-beta T cell | unfloored | 1 | 1 | 1 | NA | NA | 1 | NA | 1 | NA | 74 | 15164 | false | ortholog_mapped=1;ortholog_collisions=0;human_genes_present=1;sex_adjusted=true;below_old_floor=false |
| Lung::T cell | Lung | Lung_Integrated | T cell | unfloored | 2 | 1 | 1 | NA | NA | 1 | NA | 0.49151 | NA | 35 | 379 | false | ortholog_mapped=1;ortholog_collisions=0;human_genes_present=1;sex_adjusted=true;below_old_floor=false |
| Lung::alveolar macrophage | Lung | Lung_CellRef | alveolar macrophage | unfloored | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | 13 | 42464 | false | ortholog_mapped=0;ortholog_collisions=0;human_genes_present=0;sex_adjusted=true;below_old_floor=false |
| Lung::alveolar macrophage | Lung | Lung_Integrated | alveolar macrophage | unfloored | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | 71 | 50125 | false | ortholog_mapped=0;ortholog_collisions=0;human_genes_present=0;sex_adjusted=true;below_old_floor=false |
| Lung::classical monocyte | Lung | ImmuneAging_BAL | classical monocyte | unfloored | 7 | 4 | 1 | 0.8 | 0.94774 | 0.25437 | 0.35182 | 0.077922 | 0.22877 | 8 | 45354 | false | ortholog_mapped=4;ortholog_collisions=0;human_genes_present=4;sex_adjusted=true;below_old_floor=false |
| Lung::classical monocyte | Lung | Lung_CellRef | classical monocyte | unfloored | 7 | 4 | 0.25 | 0.4 | 0.55247 | 1 | 0.58821 | 1 | 0.65335 | 9 | 99 | false | ortholog_mapped=4;ortholog_collisions=0;human_genes_present=4;sex_adjusted=true;below_old_floor=false |
| Lung::classical monocyte | Lung | Lung_Integrated | classical monocyte | unfloored | 7 | 4 | 0.25 | -0.4 | -0.76276 | 1 | 0.57871 | 0.33766 | 0.82517 | 74 | 12991 | false | ortholog_mapped=4;ortholog_collisions=0;human_genes_present=4;sex_adjusted=true;below_old_floor=false |
| Lung::non-classical monocyte | Lung | ImmuneAging_BAL | non-classical monocyte | unfloored | 7 | 5 | 0 | -0.7 | -0.92229 | 0.093453 | 0.1984 | 0.13187 | 0.30669 | 8 | 7309 | false | ortholog_mapped=5;ortholog_collisions=0;human_genes_present=5;sex_adjusted=true;below_old_floor=false |
| Lung::non-classical monocyte | Lung | Lung_CellRef | non-classical monocyte | unfloored | 7 | 5 | 0.2 | 0 | 0.34956 | 0.69765 | 1 | 0.12987 | 0.78422 | 4 | 4 | true | ortholog_mapped=5;ortholog_collisions=0;human_genes_present=5;sex_adjusted=false;below_old_floor=true |
| Lung::non-classical monocyte | Lung | Lung_Integrated | non-classical monocyte | unfloored | 7 | 5 | 0.6 | -0.1 | 0.30437 | 1 | 0.95502 | 0.7023 | 0.9001 | 65 | 6931 | false | ortholog_mapped=5;ortholog_collisions=0;human_genes_present=5;sex_adjusted=true;below_old_floor=false |
| Lung::vein endothelial cell | Lung | Lung_CellRef | vein endothelial cell | unfloored | 1 | 1 | 0 | NA | NA | 1 | NA | 0.49251 | NA | 11 | 1635 | false | ortholog_mapped=1;ortholog_collisions=0;human_genes_present=1;sex_adjusted=true;below_old_floor=false |
| Lung::vein endothelial cell | Lung | Lung_Integrated | vein endothelial cell | unfloored | 1 | 1 | 0 | NA | NA | 1 | NA | 1 | NA | 63 | 5755 | false | ortholog_mapped=1;ortholog_collisions=0;human_genes_present=1;sex_adjusted=true;below_old_floor=false |
| Marrow::hematopoietic precursor cell | Marrow | Yazar_Blood | hematopoietic precursor cell | unfloored | 221 | 187 | 0.34225 | 0.031554 | -0.046851 | 0.009995 | 0.68716 | 0.017982 | 0.7013 | 670 | 1718 | false | ortholog_mapped=188;ortholog_collisions=0;human_genes_present=187;sex_adjusted=true;below_old_floor=false |
| Marrow::naive B cell | Marrow | Yazar_Blood | naive B cell | unfloored | 0 | 0 | NA | NA | NA | NA | NA | NA | NA | 931 | 61914 | false | ortholog_mapped=0;ortholog_collisions=0;human_genes_present=0;sex_adjusted=true;below_old_floor=false |

## Effect of removing the floor

Removing the floor added 1 pairs, from 27 to 28 pairs, and changed the pooled pair-gene row count from 294 to 299. Sign agreement moved from 0.391156 to 0.38796; Spearman moved from 0.042512 to 0.0270615; Pearson moved from 0.0580543 to 0.0498349. The p values for both nulls are in the headline table.

Pairs entering only in the unfloored arm:

| mouse unit | human cohort | human unit | donors 30+ | cells 30+ | tested genes |
| --- | --- | --- | --- | --- | --- |
| Lung::non-classical monocyte | Lung_CellRef | non-classical monocyte | 4 | 4 | 5 |

## Independent spot-check and consistency checks

The spot-check pair is `Lung::non-classical monocyte` versus `ImmuneAging_BAL::non-classical monocyte`. Starting again from `results/followup-mature-p1m-genes.csv` and the unit DE CSV, the script checked every mouse log2FC and ortholog against the pair-gene audit rows, recomputed same-sign indicators as `sign(mouse_log2FC) == sign(human_slope)`, and obtained 0/5 = 0.000000000000. This exactly equals the unfloored CSV fraction, and the tested-gene count also equals 5. The non-committed audit table is `/home/reza/projects/agescape-revision-data/matched-p1m/pair-gene-rows.csv`.

The report's pooled values were computed from the same serialized pair-gene rows used for these checks, not copied by hand. Internal assertions confirmed 36 total match-audit rows, unique pair-arm keys, zero tested genes whenever donors <2, and exact spot-check agreement. After writing, the final CSV and pooled summaries were independently re-read and compared before completion.

## Limits

- Same-named cell types are not proven homologous or developmentally equivalent across mouse and human; exact-name matching is deliberately conservative and excludes many plausible synonyms.
- The P1m signature is marrow-dominated. Exact blood matching leaves two Yazar marrow pairs, and the hematopoietic-precursor pair alone supplies 187 of 299 unfloored pair-gene rows (187 of 294 floored). Thus pooled pair-gene statistics are dominated by that marrow pair despite 14 matched kidney/lung mouse units, and the 36 derivation units are not represented equally.
- Power is limited by small donor numbers, zero-gene pairs, and unequal gene counts per pair. Pair-gene rows are not independent biological replicates; permutation p values address the stated nulls, not all clustering or selection uncertainty.
- Kidney Normal Adult and SN v15 source files are absent. Their committed metadata can document exact-name matches, but no expression slopes were invented for them. Perez and PBMC contributed no exact normalized matches under this rule.
- Mouse log2FC is a selected two-group 18m-to-30m effect, while human slopes are continuous cross-sectional 30+ associations on a different expression scale. Correlation magnitudes are descriptive concordance measures, not evidence of conserved causal aging programs.
- The mouse 30m group is all male, as documented in the derivation, so mouse age and sex are confounded; human sex adjustment does not remove that derivation limitation.

## Reproducibility

Script: `scripts/run_mature_matched_p1m.py`. Started `2026-11-24T08:49:24.996674+00:00`; runtime 19.813 seconds. The full run used the `agescape` conda environment and sources under `/home/reza/projects/agescape-revision-data`. Human pseudobulks are cached outside the repository at `/home/reza/projects/agescape-revision-data/matched-p1m/human-units.pkl`; rerun with `--refresh-human-cache` to rebuild them. Outputs are `results/followup-mature-matched-p1m.csv` and this report.
