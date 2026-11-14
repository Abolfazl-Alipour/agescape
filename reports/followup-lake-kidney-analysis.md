# Restricted public Lake kidney analysis

## Verdict

The restricted public, assay-specific analysis is complete. It uses donors—not cells or libraries—as inferential units and does not pool scCv3 with snCv3. This is not the blocked full-45-reference-donor analysis and not controlled raw-sequencing work.

No result below is called replication merely because it is non-significant. The prior `Kidney_SCR_v15` and `Kidney_SN_v15` analyses use overlapping Lake releases and are release/method comparisons, not independent replication. The complete estimates, sensitivity models, per-30-year effects, deterministic `slope*30/122.5` lifespan-scaled effects, confidence intervals, and BH families are in the aggregate CSVs listed below.

## Frozen design and resources

The design was frozen before expression download or outcomes in `data/manifest/lake_kidney_analysis_plan.json` at HEAD `a335764`. Primary models include healthy exact-age adults aged 20–120 years. snCv3 donor 3490 (14 years) is excluded from every primary model and included only in the declared all-exact-age P1/CLR sensitivity. Decade-only donors 18-142, 18-162, and 18-312 are never assigned numeric ages.

Eligible cell types require at least six donors, 20 cells per donor, and identifiable covariates. P1 uses established `P1_Clean`, per-cell log1p(CPM), donor means, and age+sex OLS with a 200-draw expression-decile-matched random-signature null (seed 20261105). CLR uses all released non-unknown cell types per assay and a 0.5-cell pseudocount. Hellinger uses the established 20–40-year reference, cross-fitting by donor, 1,500-UMI cell filter, technical-gene exclusion, and reference/donor floors. No all-gene donor-pseudobulk model was predeclared or added.

Downloads and processing ran sequentially through `scripts/resource_guard.py` with a 30-GB disk reserve, 3-GB download cap, and 6-GB memory per processing job. scCv3 finished before snCv3 started. No commit or push was made.

## Verified inputs and exclusions

- **scCv3**: 107,344 rows, raw.X 107,344 x 26,837; 170,393,115 stored values all finite, nonnegative, and integer-valued. 18 healthy donors (18 exact-age; 0 decade-only), 18 primary adults. GEO reconciliation: 18 donors, 20 specimens, 20 libraries. Assay=["10x 3' v3"]; suspension=['cell']; release-wide is_primary_data=['False'].
- **snCv3**: 172,847 rows, raw.X 172,847 x 28,813; 326,679,849 stored values all finite, nonnegative, and integer-valued. 13 healthy donors (10 exact-age; 3 decade-only), 9 primary adults. GEO reconciliation: 13 donors, 18 specimens, 21 libraries. Assay=["10x 3' v3"]; suspension=['nucleus']; release-wide is_primary_data=['False'].

Sources are the public [CELLxGENE collection](https://cellxgene.cziscience.com/collections/bcb61471-2a44-4d00-a0af-ff085512674c) per-assay H5ADs and the public GEO metadata recorded in `data/manifest/lake_kidney.json`. Each object's own disease, development stage, donor ID, sex, assay, suspension, and released `cell_type` annotations were used. Feature IDs are unique; duplicate feature symbols are retained as distinct released features, and P1 symbol matching includes every matching feature as in the established scorer. GEO supplies specimen/library reconciliation but not ages.

The release-wide `is_primary_data=False` flag was recorded rather than silently discarding every cell. Each per-assay object was analyzed independently, with no cross-assay concatenation. Disease cells, decade-only ages, unknown annotations, and cell types failing the frozen donor/cell floors were excluded from the relevant models.

## P1 transport


**scCv3: P1 transport primary donor models**

| Cell type | Donors | Slope/year (95% CI) | p | BH q | Matched-null p | Null BH q |
|---|---:|---:|---:|---:|---:|---:|
| T cell | 8 | 0.0048668 [-0.0096509, 0.019384] | 0.42821 | 0.77078 | 0.094527 | 0.17015 |
| endothelial cell | 18 | 0.0077427 [-0.0020315, 0.017517] | 0.112 | 0.50401 | 0.014925 | 0.044776 |
| epithelial cell of proximal tubule | 18 | 0.00061036 [-0.013504, 0.014724] | 0.92778 | 0.95553 | 0.26368 | 0.35821 |
| kidney collecting duct intercalated cell | 14 | -0.012871 [-0.028311, 0.0025681] | 0.093686 | 0.50401 | 0.014925 | 0.044776 |
| kidney collecting duct principal cell | 17 | 0.0016139 [-0.012028, 0.015256] | 0.80339 | 0.95553 | 0.054726 | 0.12313 |
| kidney connecting tubule epithelial cell | 13 | -0.0067342 [-0.018299, 0.0048307] | 0.22361 | 0.58254 | 0.53234 | 0.53234 |
| kidney distal convoluted tubule epithelial cell | 14 | -0.011862 [-0.033791, 0.010068] | 0.25891 | 0.58254 | 0.0099502 | 0.044776 |
| kidney interstitial cell | 17 | 0.0019353 [-0.0057468, 0.0096175] | 0.59746 | 0.89619 | 0.27861 | 0.35821 |
| kidney loop of Henle thick ascending limb epithelial cell | 18 | -0.00046269 [-0.017857, 0.016931] | 0.95553 | 0.95553 | 0.48259 | 0.53234 |

**snCv3: P1 transport primary donor models**

| Cell type | Donors | Slope/year (95% CI) | p | BH q | Matched-null p | Null BH q |
|---|---:|---:|---:|---:|---:|---:|
| T cell | 6 | 0.0014176 [-0.0064462, 0.0092815] | 0.60635 | 0.90392 | 0.25373 | 0.47122 |
| endothelial cell | 9 | 0.00053477 [-0.0066077, 0.0076773] | 0.86067 | 0.90392 | 0.51244 | 0.66617 |
| epithelial cell of proximal tubule | 9 | 0.0038007 [-0.004925, 0.012526] | 0.32751 | 0.90392 | 0.0049751 | 0.032338 |
| kidney collecting duct intercalated cell | 9 | -0.0016746 [-0.0067847, 0.0034354] | 0.4532 | 0.90392 | 0.24378 | 0.47122 |
| kidney collecting duct principal cell | 9 | 0.0015641 [-0.015887, 0.019016] | 0.83368 | 0.90392 | 0.0049751 | 0.032338 |
| kidney connecting tubule epithelial cell | 6 | -0.0018764 [-0.023288, 0.019535] | 0.79845 | 0.90392 | 0.10945 | 0.28458 |
| kidney interstitial alternatively activated macrophage | 7 | 0.00051803 [-0.0068868, 0.0079229] | 0.85546 | 0.90392 | 0.59204 | 0.69968 |
| kidney interstitial fibroblast | 8 | 0.0022379 [-0.0092148, 0.013691] | 0.63679 | 0.90392 | 0.079602 | 0.25871 |
| kidney loop of Henle thick ascending limb epithelial cell | 9 | 0.0021493 [-0.0019554, 0.006254] | 0.2474 | 0.90392 | 0.0099502 | 0.043118 |
| kidney loop of Henle thin ascending limb epithelial cell | 7 | -0.0010684 [-0.024144, 0.022007] | 0.90392 | 0.90392 | 0.79104 | 0.79104 |
| kidney loop of Henle thin descending limb epithelial cell | 8 | -0.00091724 [-0.016827, 0.014993] | 0.88798 | 0.90392 | 0.73134 | 0.79104 |
| mononuclear phagocyte | 6 | 0.0017819 [-0.0073918, 0.010956] | 0.58022 | 0.90392 | 0.29353 | 0.47699 |
| renal interstitial pericyte | 8 | 0.00092217 [-0.010135, 0.011979] | 0.83871 | 0.90392 | 0.46766 | 0.66617 |

The conventional donor-model p-values and the matched-random-null empirical p-values are separate inferential summaries with separate assay-wise BH families. A small matched-null q does not make a non-significant conventional donor-model result a confirmed replication.

## Donor-level CLR composition


**scCv3: CLR composition primary donor models**

| Cell type | Donors | Slope/year (95% CI) | p | BH q |
|---|---:|---:|---:|---:|
| T cell | 18 | 0.055281 [0.00078947, 0.10977] | 0.047159 | 0.21221 |
| endothelial cell | 18 | 0.011062 [-0.025628, 0.047753] | 0.53016 | 0.68164 |
| epithelial cell of proximal tubule | 18 | -0.003836 [-0.061418, 0.053746] | 0.88897 | 0.88897 |
| kidney collecting duct intercalated cell | 18 | -0.04068 [-0.10429, 0.022932] | 0.19298 | 0.37325 |
| kidney collecting duct principal cell | 18 | 0.039794 [-0.024573, 0.10416] | 0.20736 | 0.37325 |
| kidney connecting tubule epithelial cell | 18 | -0.013621 [-0.10359, 0.076354] | 0.7514 | 0.84533 |
| kidney distal convoluted tubule epithelial cell | 18 | -0.054775 [-0.13977, 0.030217] | 0.18973 | 0.37325 |
| kidney interstitial cell | 18 | 0.034866 [0.0016007, 0.068132] | 0.041131 | 0.21221 |
| kidney loop of Henle thick ascending limb epithelial cell | 18 | 0.021986 [-0.029956, 0.073929] | 0.38121 | 0.57181 |

**snCv3: CLR composition primary donor models**

| Cell type | Donors | Slope/year (95% CI) | p | BH q |
|---|---:|---:|---:|---:|
| T cell | 9 | 0.013957 [-0.01712, 0.045033] | 0.31394 | 0.95374 |
| endothelial cell | 9 | 0.0011867 [-0.046831, 0.049204] | 0.95374 | 0.95374 |
| epithelial cell of proximal tubule | 9 | -0.026923 [-0.080194, 0.026348] | 0.26242 | 0.95374 |
| kidney collecting duct intercalated cell | 9 | -0.025633 [-0.065398, 0.014131] | 0.16579 | 0.95374 |
| kidney collecting duct principal cell | 9 | 0.0083579 [-0.044251, 0.060967] | 0.71089 | 0.95374 |
| kidney connecting tubule epithelial cell | 9 | -0.020145 [-0.10763, 0.067337] | 0.59353 | 0.95374 |
| kidney interstitial alternatively activated macrophage | 9 | 0.0048166 [-0.033263, 0.042897] | 0.76741 | 0.95374 |
| kidney interstitial fibroblast | 9 | 0.015793 [-0.074049, 0.10563] | 0.68212 | 0.95374 |
| kidney loop of Henle thick ascending limb epithelial cell | 9 | 0.0028676 [-0.049238, 0.054974] | 0.89728 | 0.95374 |
| kidney loop of Henle thin ascending limb epithelial cell | 9 | 0.0035061 [-0.12001, 0.12702] | 0.94688 | 0.95374 |
| kidney loop of Henle thin descending limb epithelial cell | 9 | 0.015609 [-0.095108, 0.12633] | 0.74189 | 0.95374 |
| mononuclear phagocyte | 9 | 0.013233 [-0.021382, 0.047848] | 0.38567 | 0.95374 |
| renal interstitial pericyte | 9 | 0.0067886 [-0.053295, 0.066872] | 0.79147 | 0.95374 |

These are associations among recovered annotated cell fractions, not absolute tissue cell counts. Zero cell-type counts remain in each donor's composition through the declared 0.5 pseudocount. Total-cell-adjusted and all-exact-age sensitivities are in `results/followup-lake-kidney-composition-models.csv`.

## Hellinger displacement


**scCv3: Hellinger primary donor models**

| Cell type | Donors | Slope/year (95% CI) | p | BH q |
|---|---:|---:|---:|---:|
| T cell | 8 | -0.00068055 [-0.0020274, 0.00066627] | 0.23329 | 0.34993 |
| endothelial cell | 18 | 0.00033432 [-0.00010916, 0.00077781] | 0.12821 | 0.30071 |
| epithelial cell of proximal tubule | 17 | -2.4102e-05 [-0.00056995, 0.00052174] | 0.92546 | 0.92546 |
| kidney collecting duct intercalated cell | 14 | 0.00068496 [-0.00033929, 0.0017092] | 0.16706 | 0.30071 |
| kidney collecting duct principal cell | 17 | 0.00035926 [-6.611e-05, 0.00078462] | 0.091119 | 0.30071 |
| kidney connecting tubule epithelial cell | 13 | 0.00080755 [-0.00010042, 0.0017155] | 0.075089 | 0.30071 |
| kidney distal convoluted tubule epithelial cell | 13 | 0.00019604 [-0.0003923, 0.00078439] | 0.47024 | 0.55032 |
| kidney interstitial cell | 16 | -0.00022302 [-0.00090403, 0.00045799] | 0.48917 | 0.55032 |
| kidney loop of Henle thick ascending limb epithelial cell | 18 | 0.00038353 [-0.00017284, 0.00093989] | 0.16142 | 0.30071 |

snCv3 Hellinger is **untestable** under the frozen reference rule: fewer_than_two_young_reference_donors_with_20_qc_cells. Only KRP446 (29 years) is an exact-age primary adult in the 20–40-year reference window, so the required two young reference donors are unavailable. This is an eligibility outcome, not a null result. scCv3 reference-floor details are in `results/followup-lake-kidney-hellinger-eligibility.csv`; donor-level distances are in `results/followup-lake-kidney-hellinger-donors.csv`.

## Sensitivity and multiplicity families

BH was applied separately by assay, endpoint, donor set, and labeled model. The table summarizes every fitted family; it does not pool assays or endpoints.

| Endpoint | Assay | Donor set | Model | Tests | Minimum BH q | BH q<0.05 |
|---|---|---|---|---:|---:|---:|
| CLR_composition | scCv3 | all_exact_age | primary_age_sex | 9 | 0.21221 | 0 |
| CLR_composition | scCv3 | primary_adult | primary_age_sex | 9 | 0.21221 | 0 |
| CLR_composition | scCv3 | primary_adult | total_cells_sensitivity | 9 | 0.24739 | 0 |
| CLR_composition | snCv3 | all_exact_age | primary_age_sex | 13 | 0.93014 | 0 |
| CLR_composition | snCv3 | primary_adult | primary_age_sex | 13 | 0.95374 | 0 |
| CLR_composition | snCv3 | primary_adult | total_cells_sensitivity | 13 | 0.99743 | 0 |
| Hellinger | scCv3 | primary_adult | primary_depth_adjusted | 9 | 0.30071 | 0 |
| Hellinger | scCv3 | primary_adult | unadjusted_sensitivity | 9 | 0.44808 | 0 |
| P1_transport | scCv3 | all_exact_age | primary_age_sex | 9 | 0.50401 | 0 |
| P1_transport | scCv3 | primary_adult | depth_sensitivity | 9 | 0.53786 | 0 |
| P1_transport | scCv3 | primary_adult | primary_age_sex | 9 | 0.50401 | 0 |
| P1_transport | snCv3 | all_exact_age | primary_age_sex | 13 | 0.73564 | 0 |
| P1_transport | snCv3 | primary_adult | depth_sensitivity | 13 | 0.99437 | 0 |
| P1_transport | snCv3 | primary_adult | primary_age_sex | 13 | 0.90392 | 0 |
| P1_transport matched null | scCv3 | primary_adult | expression-decile-matched random signatures | 9 | 0.044776 | 3 |
| P1_transport matched null | snCv3 | primary_adult | expression-decile-matched random signatures | 13 | 0.032338 | 3 |

The lifespan quantity is only the predeclared deterministic rescaling of each slope by 30/122.5, not an additional test. The all-exact-age snCv3 sensitivity is the only lifespan/age-range sensitivity beyond primary adult models and depth/total-cell adjustment.

## Comparison with existing kidney results

Only exact cell-type labels and compatible donor-level estimands were compared. For snCv3 P1, the all-exact-age model is used because the prior 10-donor `Kidney_SN_v15` model included the 14-year-old; the primary adult-only snCv3 model is not presented as the same donor set. CLR was not compared across releases because the released composition universes and QC differ.

| Endpoint | Assay/prior cohort | Cell type | Shared donor IDs | New slope | Prior slope | New q | Prior q |
|---|---|---|---:|---:|---:|---:|---:|
| P1_transport | scCv3/Kidney_SCR_v15 | T cell | 18 | 0.0048668 | 0.0023889 | 0.77078 | 0.91555 |
| P1_transport | scCv3/Kidney_SCR_v15 | endothelial cell | 18 | 0.0077427 | 0.0053901 | 0.50401 | 0.62827 |
| P1_transport | scCv3/Kidney_SCR_v15 | epithelial cell of proximal tubule | 18 | 0.00061036 | 0.0026012 | 0.95553 | 0.93928 |
| P1_transport | scCv3/Kidney_SCR_v15 | kidney collecting duct intercalated cell | 18 | -0.012871 | -0.0016307 | 0.50401 | 0.95278 |
| P1_transport | scCv3/Kidney_SCR_v15 | kidney collecting duct principal cell | 18 | 0.0016139 | -0.0082945 | 0.95553 | 0.60206 |
| P1_transport | scCv3/Kidney_SCR_v15 | kidney connecting tubule epithelial cell | 18 | -0.0067342 | -0.0024006 | 0.58254 | 0.91555 |
| P1_transport | scCv3/Kidney_SCR_v15 | kidney distal convoluted tubule epithelial cell | 18 | -0.011862 | -0.015348 | 0.58254 | 0.54731 |
| P1_transport | scCv3/Kidney_SCR_v15 | kidney interstitial cell | 18 | 0.0019353 | -0.00055868 | 0.89619 | 0.96419 |
| P1_transport | scCv3/Kidney_SCR_v15 | kidney loop of Henle thick ascending limb epithelial cell | 18 | -0.00046269 | -0.00014632 | 0.95553 | 0.98939 |
| P1_transport | snCv3/Kidney_SN_v15 | T cell | 10 | 0.0014176 | 0.0008266 | 0.90392 | 0.94313 |
| P1_transport | snCv3/Kidney_SN_v15 | endothelial cell | 10 | 0.00097755 | 0.00035204 | 0.90392 | 0.95278 |
| P1_transport | snCv3/Kidney_SN_v15 | epithelial cell of proximal tubule | 10 | 0.0038007 | 0.0036122 | 0.85152 | 0.79797 |
| P1_transport | snCv3/Kidney_SN_v15 | kidney collecting duct intercalated cell | 10 | -0.0021945 | -0.0030654 | 0.73564 | 0.47557 |
| P1_transport | snCv3/Kidney_SN_v15 | kidney collecting duct principal cell | 10 | 0.0010028 | 0.00046097 | 0.90392 | 0.98136 |
| P1_transport | snCv3/Kidney_SN_v15 | kidney connecting tubule epithelial cell | 10 | -0.0018764 | -0.0025269 | 0.90392 | 0.94313 |
| P1_transport | snCv3/Kidney_SN_v15 | kidney interstitial alternatively activated macrophage | 10 | 0.0023178 | 0.0014075 | 0.73564 | 0.80215 |
| P1_transport | snCv3/Kidney_SN_v15 | kidney interstitial fibroblast | 10 | 0.0045957 | 0.0038548 | 0.73564 | 0.62827 |
| P1_transport | snCv3/Kidney_SN_v15 | kidney loop of Henle thick ascending limb epithelial cell | 10 | 0.0016666 | -0.00015056 | 0.73564 | 0.98421 |
| P1_transport | snCv3/Kidney_SN_v15 | kidney loop of Henle thin ascending limb epithelial cell | 10 | -0.0010684 | -0.0012484 | 0.90392 | 0.96419 |
| P1_transport | snCv3/Kidney_SN_v15 | kidney loop of Henle thin descending limb epithelial cell | 10 | 0.0012344 | 0.00011007 | 0.90392 | 0.98939 |
| P1_transport | snCv3/Kidney_SN_v15 | mononuclear phagocyte | 10 | 0.0017819 | 0.0021344 | 0.90392 | 0.91555 |
| P1_transport | snCv3/Kidney_SN_v15 | renal interstitial pericyte | 10 | 0.00092217 | 0.0012198 | 0.90392 | 0.95278 |
| Hellinger | scCv3/Kidney_SCR_v15 | T cell | 18 | -0.00068055 | -0.00045778 | 0.34993 | 0.58756 |
| Hellinger | scCv3/Kidney_SCR_v15 | endothelial cell | 18 | 0.00033432 | 0.00026948 | 0.30071 | 0.58756 |
| Hellinger | scCv3/Kidney_SCR_v15 | epithelial cell of proximal tubule | 18 | -2.4102e-05 | 0.00057495 | 0.92546 | 0.45445 |
| Hellinger | scCv3/Kidney_SCR_v15 | kidney collecting duct intercalated cell | 18 | 0.00068496 | 5.132e-05 | 0.30071 | 0.97118 |
| Hellinger | scCv3/Kidney_SCR_v15 | kidney collecting duct principal cell | 18 | 0.00035926 | 0.00096261 | 0.30071 | 0.58756 |
| Hellinger | scCv3/Kidney_SCR_v15 | kidney connecting tubule epithelial cell | 18 | 0.00080755 | 1.3756e-05 | 0.30071 | 0.97118 |
| Hellinger | scCv3/Kidney_SCR_v15 | kidney distal convoluted tubule epithelial cell | 18 | 0.00019604 | -0.00061739 | 0.55032 | 0.10083 |
| Hellinger | scCv3/Kidney_SCR_v15 | kidney interstitial cell | 18 | -0.00022302 | -0.00037036 | 0.55032 | 0.58756 |
| Hellinger | scCv3/Kidney_SCR_v15 | kidney loop of Henle thick ascending limb epithelial cell | 18 | 0.00038353 | 4.7593e-05 | 0.30071 | 0.97118 |

Full compatible comparisons with the separate `Kidney_NormalAdult` and `Kidney_LivingDonor` studies are in `results/followup-lake-kidney-existing-comparison.csv`, labeled `separate_study`. They are not pooled with Lake and do not turn overlapping releases into independent evidence.

## Untestable, not attempted, and interpretation limits

- scCv3: 16 released cell types failed the six-donor/20-cell gate: B cell, conventional dendritic cell, cytotoxic T cell, kidney interstitial alternatively activated macrophage, kidney loop of Henle thin ascending limb epithelial cell, kidney loop of Henle thin descending limb epithelial cell, mature NK T cell, monocyte, mononuclear phagocyte, natural killer cell, non-classical monocyte, parietal epithelial cell, plasma cell, plasmacytoid dendritic cell, human, podocyte, regulatory T cell.
- snCv3: 13 released cell types failed the six-donor/20-cell gate: B cell, conventional dendritic cell, kidney distal convoluted tubule epithelial cell, mast cell, mature NK T cell, neural cell, neutrophil, non-classical monocyte, papillary tips cell, parietal epithelial cell, plasma cell, plasmacytoid dendritic cell, human, podocyte.
- snCv3 Hellinger: untestable because the young-reference floor fails, as detailed above.
- Full 45-reference-donor analysis: still blocked by incomplete public exact-age/common-mapping coverage; no ages were imputed.
- All-gene donor-pseudobulk models: not predeclared and not performed; there is therefore no gene-level multiplicity family or gene-level claim.
- Controlled KPMP/NIDDK/dbGaP raw sequencing: not requested, contacted, downloaded, or analyzed.
- SNARE2, Slide-seq2, Visium, muscle, and other follow-up cohorts: outside this restricted task and not touched.
- Non-significant P1, CLR, or Hellinger results are not equivalence, absence-of-effect, or replication claims. Small donor numbers and release differences limit inference.

## Reproducibility and delivered aggregates

- `scripts/run_followup_lake_kidney.py`: backed raw-count validation, donor aggregation, P1/null, CLR, Hellinger, comparison, and manifest update.
- `scripts/report_followup_lake_kidney.py`: renders this report from aggregate CSVs.
- `results/followup-lake-kidney-*-verification.json`: per-assay raw-count and metadata gates.
- `results/followup-lake-kidney-donor-crosswalk.csv`: exact ages, health, primary/sensitivity sets, GEO patient/specimen/library mapping.
- `results/followup-lake-kidney-celltype-eligibility.csv` and `-hellinger-eligibility.csv`: inclusion gates and untestable reasons.
- `results/followup-lake-kidney-p1-donor-scores.csv`, `-p1-models.csv`: donor scores, null summaries, effects/CIs, and both BH families.
- `results/followup-lake-kidney-composition-donor-counts.csv`, `-composition-models.csv`: zero-inclusive donor counts and CLR models.
- `results/followup-lake-kidney-hellinger-donors.csv`, `-hellinger-models.csv`: eligible donor distances and models.
- `results/followup-lake-kidney-existing-comparison.csv`: exact-label compatible prior/separate-study comparisons and shared donor counts.
- `tests/test_followup_lake_kidney.py`: helper and delivered-output validation.

Run with the Agescape environment:

```bash
python scripts/resource_guard.py run --path /home/reza/projects/agescape-revision-data/lake_kidney   --min-free-gb 30 --memory-gb 6 --   python scripts/run_followup_lake_kidney.py run --assay scCv3
# Repeat for snCv3 only after scCv3 completes, then:
python scripts/resource_guard.py run --path /home/reza/projects/agescape-revision-data/lake_kidney   --min-free-gb 30 --memory-gb 6 --   python scripts/run_followup_lake_kidney.py finalize
python scripts/report_followup_lake_kidney.py
```

## Validation record

- `python -m pytest tests/test_followup_lake_kidney.py -q`: **10 passed**.
- `python -m pytest -q`: **116 passed, 1 skipped** after the TS2 manifest portability fix in commit `1c527c0`. The Lake manifest passes the same path/security checks.
- `git diff --check`: passed. HEAD at the frozen design was `a335764`; no commit or push was performed and parent review is required.
