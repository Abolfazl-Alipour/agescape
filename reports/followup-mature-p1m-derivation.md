# Mature P1m and P0m derivation

## Scope and runtime

This is Phase 2 Item 1 only. The derivation uses exactly 18 months as the reference and 30 months as the target. The 21-month and 24-month animals are excluded before unit discovery, pseudobulk aggregation, and model fitting. No 21-month or 24-month animal enters any DE model.

- Model used: `kimi-code/k3.1-pro`.
- Successful run start: `2026-11-21T18:08:15.575552+00:00`.
- Successful run end: `2026-11-21T18:11:04.778714+00:00`.
- Successful run duration: 169.203 seconds, or 2 minutes 49.203 seconds.
- Inputs: all 16 `data/interim/m1-qc-cleaned-*.h5ad` files were present and read.
- Implementation: `scripts/run_mature_p1m_derivation.py`.

An initial run started at `2026-11-21T18:07:29.180040+00:00` and stopped on the first fitted unit with `KeyError: 'size_factors'`. PyDESeq2 0.5.4 stores size factors in `dds.obs`, not the `dds.obsm` location initially used. Only `results/followup-mature-p1m-eligibility.csv` existed after that failed run. The API was inspected, the orientation diagnostic was corrected, and the corrected full run was then performed. This was not a blind retry.

## Recovered construction rules

The new script preserves the original M2 and M5 construction rules with the requested age and eligibility changes:

1. Raw counts are summed per mouse by tissue and cell type to make pseudobulk samples.
2. Genes must have at least 10 total counts across the samples in a unit.
3. The requested design is `~ sex + age_group`. Sex is dropped only if it has no variation or the intercept, sex, and age design matrix is rank deficient.
4. The contrast is explicitly `['age_group', '30m', '18m']`. Thus positive log2FC means higher expression at 30m than at 18m.
5. DE eligibility requires at least 3 mice and at least 50 cells in both the 18m and 30m groups.
6. The significance rule is FDR < 0.05 and absolute log2FC >= 0.25.
7. P0m uses the original M2 PyDESeq2 independent-filtering result. P1m uses the original uniform-BH alternative: `independent_filter=False`, followed by independent BH recomputation over all non-missing raw p-values in each unit.
8. P1m additionally preserves the floor-10 rule: at least 3 mice per age must each contribute at least 10 cells. This rule does not restrict P0m beyond the initial eligibility rule.
9. P0m pools all significant genes across eligible units without clean-gene filtering. P1m pools significant genes only from floor-10 units and removes `Gm` predicted labels, `-ps` pseudogene labels, `Rik` labels matching the original regular expression, non-protein-coding or unknown biotypes, ribosomal and mitochondrial prefixes, `Rn45s`, lung immune `Mgp` and `Col3a1`, and kidney `Ins1` and `Ins2`.
10. P2/P3 sex-dimorphic and immediate-early filters are not applied.
11. For both tiers, gene direction is `UP` when the mean log2FC across significant units is positive and `DOWN` otherwise. There is no additional vote or sign-consistency requirement. The output also records the numbers of up and down supporting units.
12. Human mapping uses `data/reference/orthologs-mouse-human.csv`. All 17,063 rows are `ortholog_one2one`, with no duplicate mouse symbols and no missing human symbols. This table reproduces the original mapped counts of 188 P0 and 118 P1 genes.

The clean-gene implementation was also checked against all 243 original uniform-BH candidate rows. There were zero filter-decision mismatches.

## Sex verification and fitted designs

Direct inspection of mouse-level metadata across all 16 inputs verified that the 30m group is all male. The unique animals contributing to eligible units were:

| Age | Sex | Mouse IDs | Count |
| --- | --- | --- | ---: |
| 18m | female | 18-F-50; 18-F-51 | 2 |
| 18m | male | 18-M-52; 18-M-53 | 2 |
| 30m | female | none | 0 |
| 30m | male | 30-M-2; 30-M-3; 30-M-4; 30-M-5 | 4 |

The same counts were obtained for all source animals and for animals actually used in eligible units: 18m has 2 female and 2 male mice, while 30m has 0 female and 4 male mice.

All 36 eligible unit fits had sex variation and a full-rank design, so all 36 actually fit `~ sex + age_group`. Zero units dropped sex. PyDESeq2 did not drop or mishandle the sex term in these fits. Nevertheless, the all-male 30m group remains an age-sex confound and limits sex-generalizable interpretation.

## Gene counts and ortholog mapping

| Tier | Mouse total | Mouse up | Mouse down | Human-mapped total | Human-mapped up | Human-mapped down |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| P0m | 2598 | 1474 | 1124 | 2049 | 1157 | 892 |
| P1m | 1914 | 1119 | 795 | 1724 | 1013 | 711 |

The mapped counts are also the numbers of unique human symbols: 2,049 for P0m and 1,724 for P1m.

## Independent BH check

All 36 units passed. The script recomputed BH q-values from raw p-values with `statsmodels.stats.multitest.multipletests(..., method='fdr_bh')` and compared them with the PyDESeq2 `padj` values from the `independent_filter=False` fit. Missing-value patterns matched exactly. The maximum absolute discrepancy in the successful in-memory comparison was `2.220446049250313e-16`.

A second recomputation after CSV serialization and re-reading gave a maximum absolute discrepancy of `1.9673151996357774e-13`, consistent with floating-point serialization and still below the `1e-12` validation tolerance. No material mismatch was found. P0m's independent-filtered FDR is kept in a separate column and was not confused with uniform BH.

## Original P1 versus P1m

- Common genes: 52.
- Lost from original P1: 87.
- Gained in P1m: 1,862.
- Sign flips among common genes: 27.
- Common genes with unchanged direction: 25.

Gene-level status and both directions are recorded in `results/followup-mature-p1m-p1-overlap.csv`.

## Positive controls and orientation

The top unit by P0m DE count was `Marrow::granulocytopoietic cell`, with 1,730 P0m DE genes. Positive-control results in that unit were:

| Gene | log2FC, 30m vs 18m | Direction | Raw p | P0m FDR | Uniform FDR |
| --- | ---: | --- | ---: | ---: | ---: |
| Cdkn2a | 2.439756 | UP with age | 0.104126 | not available after independent filtering | 0.324636 |
| Spp1 | -1.015003 | DOWN with age | 0.025260 | 0.099148 | 0.128005 |
| Apoe | 0.570845 | UP with age | 0.289577 | 0.510963 | 0.582524 |
| Igfbp3 | UNTESTABLE | absent after the total-count filter | not available | not available | not available |
| Lgals3 | -0.706990 | DOWN with age | 0.061891 | 0.185502 | 0.235435 |

These are observed directions, not assertions that every control must be positive or significant.

The orientation assertion passed. The contrast was asserted to be `['age_group', '30m', '18m']`. As an independent numerical example, `Cdkn2a` in `Lung::classical monocyte` had model log2FC `+1.688323`; its mean size-factor-normalized pseudobulk expression increased from `0.221317` at 18m to `1.217836` at 30m. Thus positive log2FC is up with age.

## Eligibility and P1 floor-10 eligibility

Of 151 units observed at 18m or 30m, 36 were eligible for DE: 11 Kidney, 12 Lung, and 13 Marrow units. The other 115 were ineligible for the stated mouse or cell floors. Bladder and Tongue were directly verified to have no 30m mice. Mammary Gland, Skin, and Thymus also had no 30m mice, and individual Lung and Spleen units lacked one age side. Large Intestine had no 18m mice.

Thirty-four eligible units also passed the P1 floor-10 rule. The two exceptions were:

- `Marrow::naive B cell`: only 2 18m mice had at least 10 cells, versus 3 required.
- `Lung::CD4-positive, alpha-beta T cell`: only 2 30m mice had at least 10 cells, versus 3 required.

## Every ineligible unit and reason

| Unit | 18m mice | 18m cells | 30m mice | 30m cells | Reason |
| --- | ---: | ---: | ---: | ---: | --- |
| Spleen::B cell | 2 | 4563 | 4 | 4468 | 18m has 2 mice, requires at least 3 |
| Heart_and_Aorta::endothelial cell of coronary artery | 3 | 288 | 2 | 1961 | 30m has 2 mice, requires at least 3 |
| Liver::Kupffer cell | 1 | 221 | 3 | 1503 | 18m has 1 mouse, requires at least 3 |
| Limb_Muscle::endothelial cell | 4 | 1081 | 2 | 1182 | 30m has 2 mice, requires at least 3 |
| Large_Intestine::enterocyte of epithelium of large intestine | 0 | 0 | 3 | 1102 | no 18m mice or cells |
| Lung::lung macrophage | 2 | 16 | 4 | 1097 | 18m has 2 mice and 16 cells, requires at least 3 and 50 |
| Limb_Muscle::mesenchymal stem cell | 4 | 2943 | 2 | 1027 | 30m has 2 mice, requires at least 3 |
| Heart_and_Aorta::fibroblast of cardiac tissue | 3 | 416 | 2 | 993 | 30m has 2 mice, requires at least 3 |
| Spleen::mature NK T cell | 2 | 410 | 4 | 879 | 18m has 2 mice, requires at least 3 |
| Spleen::granulocyte | 2 | 15 | 4 | 758 | 18m has 2 mice and 15 cells, requires at least 3 and 50 |
| Spleen::T cell | 2 | 585 | 4 | 745 | 18m has 2 mice, requires at least 3 |
| Fat::mesenchymal stem cell of adipose | 3 | 829 | 1 | 655 | 30m has 1 mouse, requires at least 3 |
| Fat::endothelial cell | 3 | 630 | 1 | 632 | 30m has 1 mouse, requires at least 3 |
| Marrow::immature B cell | 4 | 37 | 3 | 627 | 18m has 37 cells, requires at least 50 |
| Heart_and_Aorta::leukocyte | 3 | 140 | 2 | 483 | 30m has 2 mice, requires at least 3 |
| Spleen::megakaryocyte-erythroid progenitor cell | 2 | 12 | 3 | 440 | 18m has 2 mice and 12 cells, requires at least 3 and 50 |
| Spleen::macrophage | 2 | 164 | 4 | 427 | 18m has 2 mice, requires at least 3 |
| Spleen::plasma cell | 2 | 120 | 4 | 367 | 18m has 2 mice, requires at least 3 |
| Large_Intestine::intestinal crypt stem cell | 0 | 0 | 3 | 313 | no 18m mice or cells |
| Fat::B cell | 3 | 330 | 1 | 301 | 30m has 1 mouse, requires at least 3 |
| Spleen::proerythroblast | 2 | 108 | 4 | 283 | 18m has 2 mice, requires at least 3 |
| Pancreas::pancreatic B cell | 2 | 1462 | 3 | 275 | 18m has 2 mice, requires at least 3 |
| Limb_Muscle::B cell | 3 | 13 | 2 | 272 | 18m has 13 cells and 30m has 2 mice, requires at least 50 and 3 |
| Pancreas::pancreatic ductal cel | 2 | 34 | 3 | 242 | 18m has 2 mice and 34 cells, requires at least 3 and 50 |
| Liver::myeloid leukocyte | 1 | 40 | 3 | 223 | 18m has 1 mouse and 40 cells, requires at least 3 and 50 |
| Heart_and_Aorta::smooth muscle cell | 3 | 34 | 2 | 213 | 18m has 34 cells and 30m has 2 mice, requires at least 50 and 3 |
| Marrow::plasma cell | 3 | 27 | 3 | 212 | 18m has 27 cells, requires at least 50 |
| Limb_Muscle::skeletal muscle satellite cell | 4 | 423 | 2 | 211 | 30m has 2 mice, requires at least 3 |
| Fat::myeloid cell | 3 | 1384 | 1 | 201 | 30m has 1 mouse, requires at least 3 |
| Heart_and_Aorta::endocardial cell | 3 | 128 | 2 | 190 | 30m has 2 mice, requires at least 3 |
| Large_Intestine::epithelial cell of large intestine | 0 | 0 | 3 | 188 | no 18m mice or cells |
| Pancreas::pancreatic acinar cell | 2 | 4 | 3 | 175 | 18m has 2 mice and 4 cells, requires at least 3 and 50 |
| Limb_Muscle::T cell | 3 | 39 | 2 | 167 | 18m has 39 cells and 30m has 2 mice, requires at least 50 and 3 |
| Limb_Muscle::macrophage | 4 | 733 | 2 | 166 | 30m has 2 mice, requires at least 3 |
| Liver::NK cell | 1 | 109 | 3 | 155 | 18m has 1 mouse, requires at least 3 |
| Kidney::podocyte | 3 | 48 | 4 | 145 | 18m has 48 cells, requires at least 50 |
| Limb_Muscle::smooth muscle cell | 4 | 242 | 2 | 143 | 30m has 2 mice, requires at least 3 |
| Fat::T cell | 3 | 338 | 1 | 137 | 30m has 1 mouse, requires at least 3 |
| Large_Intestine::large intestine goblet cell | 0 | 0 | 3 | 136 | no 18m mice or cells |
| Spleen::erythroblast | 0 | 0 | 3 | 128 | no 18m mice or cells |
| Heart_and_Aorta::erythrocyte | 2 | 13 | 2 | 122 | 18m has 2 mice and 13 cells; 30m has 2 mice; floors are 3 mice and 50 cells per age |
| Spleen::NK cell | 2 | 77 | 4 | 97 | 18m has 2 mice, requires at least 3 |
| Liver::B cell | 1 | 22 | 3 | 96 | 18m has 1 mouse and 22 cells, requires at least 3 and 50 |
| Liver::endothelial cell of hepatic sinusoid | 1 | 27 | 3 | 87 | 18m has 1 mouse and 27 cells, requires at least 3 and 50 |
| Marrow::basophil | 2 | 3 | 3 | 87 | 18m has 2 mice and 3 cells, requires at least 3 and 50 |
| Kidney::fibroblast | 3 | 35 | 4 | 84 | 18m has 35 cells, requires at least 50 |
| Heart_and_Aorta::cardiomyocyte | 2 | 52 | 2 | 82 | both ages have 2 mice, requires at least 3 |
| Lung::myeloid dendritic cell | 4 | 47 | 4 | 70 | 18m has 47 cells, requires at least 50 |
| Limb_Muscle::skeletal muscle cell | 4 | 29 | 2 | 58 | 18m has 29 cells and 30m has 2 mice, requires at least 50 and 3 |
| Kidney::plasma cell | 3 | 11 | 4 | 55 | 18m has 11 cells, requires at least 50 |
| Lung::regulatory T cell | 4 | 28 | 3 | 54 | 18m has 28 cells, requires at least 50 |
| Heart_and_Aorta::mast cell | 2 | 25 | 1 | 50 | 18m has 2 mice and 25 cells; 30m has 1 mouse; floors are 3 mice and 50 cells per age |
| Pancreas::pancreatic A cell | 2 | 155 | 3 | 50 | 18m has 2 mice, requires at least 3 |
| Limb_Muscle::Schwann cell | 4 | 46 | 2 | 46 | 18m has 46 cells; 30m has 2 mice and 46 cells; floors are 3 mice and 50 cells per age |
| Kidney::brush cell | 3 | 86 | 4 | 43 | 30m has 43 cells, requires at least 50 |
| Lung::basophil | 4 | 19 | 4 | 40 | 18m has 19 cells and 30m has 40 cells, requires at least 50 per age |
| Pancreas::leukocyte | 2 | 44 | 2 | 40 | 18m has 2 mice and 44 cells; 30m has 2 mice and 40 cells; floors are 3 mice and 50 cells per age |
| Fat::lymphocyte | 2 | 17 | 1 | 37 | 18m has 2 mice and 17 cells; 30m has 1 mouse and 37 cells; floors are 3 mice and 50 cells per age |
| Kidney::kidney loop of Henle ascending limb epithelial cell | 3 | 38 | 4 | 37 | 18m has 38 cells and 30m has 37 cells, requires at least 50 per age |
| Large_Intestine::hematopoietic stem cell | 0 | 0 | 3 | 36 | no 18m mice or cells; 30m has 36 cells, requires at least 50 |
| Pancreas::pancreatic D cell | 2 | 141 | 3 | 32 | 18m has 2 mice and 30m has 32 cells, requires at least 3 and 50 |
| Spleen::macrophage dendritic cell progenitor | 2 | 16 | 4 | 32 | 18m has 2 mice and 16 cells; 30m has 32 cells; floors are 3 mice and 50 cells per age |
| Marrow::NK cell | 3 | 24 | 3 | 31 | 18m has 24 cells and 30m has 31 cells, requires at least 50 per age |
| Liver::hepatocyte | 1 | 267 | 2 | 30 | 18m has 1 mouse; 30m has 2 mice and 30 cells; floors are 3 mice and 50 cells per age |
| Pancreas::endothelial cell | 2 | 50 | 3 | 28 | 18m has 2 mice and 30m has 28 cells, requires at least 3 and 50 |
| Lung::fibroblast of lung | 4 | 166 | 3 | 26 | 30m has 26 cells, requires at least 50 |
| Heart_and_Aorta::cardiac neuron | 3 | 19 | 2 | 23 | 18m has 19 cells; 30m has 2 mice and 23 cells; floors are 3 mice and 50 cells per age |
| Lung::adventitial cell | 4 | 56 | 4 | 23 | 30m has 23 cells, requires at least 50 |
| Lung::plasma cell | 2 | 4 | 2 | 22 | 18m has 2 mice and 4 cells; 30m has 2 mice and 22 cells; floors are 3 mice and 50 cells per age |
| Spleen::immature NKT cell | 2 | 14 | 4 | 22 | 18m has 2 mice and 14 cells; 30m has 22 cells; floors are 3 mice and 50 cells per age |
| Lung::neutrophil | 4 | 410 | 3 | 21 | 30m has 21 cells, requires at least 50 |
| Liver::plasmacytoid dendritic cell | 1 | 8 | 3 | 17 | 18m has 1 mouse and 8 cells; 30m has 17 cells; floors are 3 mice and 50 cells per age |
| Kidney::kidney mesangial cell | 3 | 27 | 3 | 15 | 18m has 27 cells and 30m has 15 cells, requires at least 50 per age |
| Fat::erythroblast | 2 | 31 | 1 | 14 | 18m has 2 mice and 31 cells; 30m has 1 mouse and 14 cells; floors are 3 mice and 50 cells per age |
| Lung::dendritic cell | 4 | 9 | 3 | 14 | 18m has 9 cells and 30m has 14 cells, requires at least 50 per age |
| Lung::plasmacytoid dendritic cell | 3 | 15 | 3 | 14 | 18m has 15 cells and 30m has 14 cells, requires at least 50 per age |
| Lung::smooth muscle cell of the pulmonary artery | 4 | 32 | 3 | 11 | 18m has 32 cells and 30m has 11 cells, requires at least 50 per age |
| Lung::type II pneumocyte | 3 | 35 | 4 | 11 | 18m has 35 cells and 30m has 11 cells, requires at least 50 per age |
| Lung::endothelial cell of lymphatic vessel | 3 | 11 | 4 | 10 | 18m has 11 cells and 30m has 10 cells, requires at least 50 per age |
| Pancreas::pancreatic stellate cell | 1 | 6 | 2 | 9 | 18m has 1 mouse and 6 cells; 30m has 2 mice and 9 cells; floors are 3 mice and 50 cells per age |
| Lung::pericyte cell | 4 | 7 | 3 | 8 | 18m has 7 cells and 30m has 8 cells, requires at least 50 per age |
| Kidney::kidney capillary endothelial cell | 3 | 19 | 2 | 7 | 18m has 19 cells; 30m has 2 mice and 7 cells; floors are 3 mice and 50 cells per age |
| Liver::hepatic stellate cell | 1 | 2 | 1 | 5 | 18m has 1 mouse and 2 cells; 30m has 1 mouse and 5 cells; floors are 3 mice and 50 cells per age |
| Fat::epithelial cell | 2 | 69 | 1 | 4 | 18m has 2 mice; 30m has 1 mouse and 4 cells; floors are 3 mice and 50 cells per age |
| Kidney::NK cell | 2 | 8 | 4 | 4 | 18m has 2 mice and 8 cells; 30m has 4 cells; floors are 3 mice and 50 cells per age |
| Lung::ciliated columnar cell of tracheobronchial tree | 3 | 9 | 3 | 3 | 18m has 9 cells and 30m has 3 cells, requires at least 50 per age |
| Lung::club cell of bronchiole | 3 | 4 | 2 | 3 | 18m has 4 cells; 30m has 2 mice and 3 cells; floors are 3 mice and 50 cells per age |
| Pancreas::pancreatic PP cell | 2 | 125 | 2 | 3 | both ages have 2 mice and 30m has 3 cells; floors are 3 mice and 50 cells per age |
| Marrow::erythroid progenitor | 1 | 1 | 2 | 2 | 18m has 1 mouse and 1 cell; 30m has 2 mice and 2 cells; floors are 3 mice and 50 cells per age |
| Bladder::bladder cell | 1 | 536 | 0 | 0 | 18m has 1 mouse; no 30m mice or cells |
| Bladder::bladder urothelial cell | 1 | 1381 | 0 | 0 | 18m has 1 mouse; no 30m mice or cells |
| Bladder::endothelial cell | 1 | 40 | 0 | 0 | 18m has 1 mouse and 40 cells; no 30m mice or cells |
| Bladder::leukocyte | 1 | 83 | 0 | 0 | 18m has 1 mouse; no 30m mice or cells |
| Lung::pulmonary interstitial fibroblast | 4 | 9 | 0 | 0 | 18m has 9 cells; no 30m mice or cells |
| Mammary_Gland::B cell | 2 | 233 | 0 | 0 | 18m has 2 mice; no 30m mice or cells |
| Mammary_Gland::T cell | 2 | 475 | 0 | 0 | 18m has 2 mice; no 30m mice or cells |
| Mammary_Gland::basal cell | 2 | 472 | 0 | 0 | 18m has 2 mice; no 30m mice or cells |
| Mammary_Gland::endothelial cell | 2 | 173 | 0 | 0 | 18m has 2 mice; no 30m mice or cells |
| Mammary_Gland::luminal epithelial cell of mammary gland | 2 | 919 | 0 | 0 | 18m has 2 mice; no 30m mice or cells |
| Mammary_Gland::macrophage | 2 | 134 | 0 | 0 | 18m has 2 mice; no 30m mice or cells |
| Mammary_Gland::stromal cell | 2 | 553 | 0 | 0 | 18m has 2 mice; no 30m mice or cells |
| Skin::T cell | 1 | 5 | 0 | 0 | 18m has 1 mouse and 5 cells; no 30m mice or cells |
| Skin::basal cell of epidermis | 1 | 40 | 0 | 0 | 18m has 1 mouse and 40 cells; no 30m mice or cells |
| Skin::epidermal cell | 1 | 23 | 0 | 0 | 18m has 1 mouse and 23 cells; no 30m mice or cells |
| Skin::keratinocyte stem cell | 1 | 20 | 0 | 0 | 18m has 1 mouse and 20 cells; no 30m mice or cells |
| Skin::stem cell of epidermis | 1 | 3 | 0 | 0 | 18m has 1 mouse and 3 cells; no 30m mice or cells |
| Thymus::DN3 thymocyte | 2 | 60 | 0 | 0 | 18m has 2 mice; no 30m mice or cells |
| Thymus::DN4 thymocyte | 2 | 343 | 0 | 0 | 18m has 2 mice; no 30m mice or cells |
| Thymus::double negative T cell | 2 | 339 | 0 | 0 | 18m has 2 mice; no 30m mice or cells |
| Thymus::immature T cell | 2 | 675 | 0 | 0 | 18m has 2 mice; no 30m mice or cells |
| Thymus::professional antigen presenting cell | 2 | 73 | 0 | 0 | 18m has 2 mice; no 30m mice or cells |
| Thymus::thymocyte | 2 | 531 | 0 | 0 | 18m has 2 mice; no 30m mice or cells |
| Tongue::Langerhans cell | 1 | 17 | 0 | 0 | 18m has 1 mouse and 17 cells; no 30m mice or cells |
| Tongue::basal cell of epidermis | 1 | 1830 | 0 | 0 | 18m has 1 mouse; no 30m mice or cells |
| Tongue::keratinocyte | 1 | 1816 | 0 | 0 | 18m has 1 mouse; no 30m mice or cells |

## Model warnings and verification limits

PyDESeq2 emitted three warnings that a dispersion trend curve did not converge and automatically switched to a mean-based dispersion trend. The run log did not retain the unit attribution for those three warnings. All fits completed, and no result was discarded because of these warnings. This fallback is a model limitation to carry forward.

No required input was missing and no required verification was left undone. Igfbp3 in the top unit is explicitly UNTESTABLE because it did not pass that unit's total-count filter. No biological null is claimed for that gene.

## Output files

- `results/followup-mature-p1m-eligibility.csv`
- `results/followup-mature-p1m-de_*.csv` for each of the 36 eligible units
- `results/followup-mature-p1m-genes.csv`
- `results/followup-mature-p1m-p0m-genes.csv`
- `results/followup-mature-p1m-summary.csv`
- `results/followup-mature-p1m-unit-ranking.csv`
- `results/followup-mature-p1m-bh-check.csv`
- `results/followup-mature-p1m-positive-controls.csv`
- `results/followup-mature-p1m-orientation-check.csv`
- `results/followup-mature-p1m-p1-overlap.csv`
