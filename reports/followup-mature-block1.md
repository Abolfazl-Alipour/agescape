# Mature Block 1 (30-plus): Long-Style Per-Donor Residual Dispersion and Composition Matching (Gao Heart + Yazar Blood)

Mature-to-old remap of the frozen 11pt Block 1. Identical rules except donors under 30 are excluded everywhere (detected sets, models, group tests, calibration, cell-count checks), the group baseline YOUNG (20 to 39) is replaced by MATURE (30 to 39), and continuous-age models use donors 30-plus only. Group tests compare MATURE (30 to 39) against OLD (60-plus). Donor is the unit in every model and permutation. Gao Heart and Yazar Blood are evaluated separately and never pooled.

Exclusions: Gao Heart excludes 8 donors under 30, leaving 109 primary donors (ages 30 to 91; 12 MATURE, 53 MID, 44 OLD). Yazar Blood excludes 49 donors under 30, leaving 932 primary donors (ages 30 to 97; 59 MATURE, 193 MID, 680 OLD). Eligible cells at post-mask UMI at least 1500 in the 30-plus sample: 482,902 Gao nuclei (was 522,761 full-range) and 592,250 Yazar cells (was 617,760 full-range). Verified after the run: zero under-30 donors appear in the mature donor-group stats (934 Gao rows, 105 distinct donors; 6,597 Yazar rows, 932 distinct donors), donor-state stats, or per-gene residual pieces for either dataset; minimum age in mature stats is 30.0 in both datasets.

Script: `scripts/run_mature_block1_residual_dispersion.py`, a copy-and-remap of `scripts/run_followup_11pt_block1_residual_dispersion.py` with the bin remap, the 30-plus exclusion, mature output paths, and the frozen narrative report builder replaced by a stub (the frozen Myeloid and naive CD4 narratives do not apply to mature data; this report is written by hand from verified CSVs). The original per-cell count cache was reused read-only through symlinks under `block1_cache_mature`; stats and diagnostics were written to the mature cache only, never to the frozen cache. Runtime: stats 78.8s Gao plus 75.1s Yazar; models and diagnostics 19.5s Gao plus 59.4s Yazar. A2 between-subject columns are omitted from the tables below because A2 is a full-range analysis and was not refit here.

## 1. Verdicts (Separate per Dataset)

| Dataset | Verdict | Non-Pooled Testable Breakdown (PROVES / DOES NOT PROVE / INCONCLUSIVE / UNTESTABLE) | Pooled Row Verdict | Summary Reason |
|---|---|---|---|---|
| Gao Heart | **DOES NOT PROVE** | 0 / 11 / 2 / 0 (of 13) | DOES NOT PROVE (powered_null(high_power<=0.25)) | 0/13_PROVES_and_11/13_DOES_NOT_PROVE(neg_sig=0) |
| Yazar Blood | **INCONCLUSIVE** (PROVES cell types: naive thymus-derived CD4-positive, alpha-beta T cell) | 1 / 13 / 2 / 12 (of 28) | DOES NOT PROVE (powered_null(high_power<=0.25)) | lineage_restricted_or_mixed(1/16_PROVES:[naive thymus-derived CD4-positive, alpha-beta T cell];13/16_DNP;2/16_INC) |

Headline changes against the original: Gao is unchanged (DOES NOT PROVE, 12 DNP to 11 DNP plus 2 INCONCLUSIVE). Yazar is unchanged (INCONCLUSIVE, 1 PROVES, 14 DNP to 13 DNP plus 2 INCONCLUSIVE). The pooled Yazar `All_Cells` row flips INCONCLUSIVE to DOES NOT PROVE (Section 6).

## 2. Pre-Specified Protocol (Frozen, Remapped to 30-plus)

- **Inputs (local, no new downloads):**
  - Gao Heart (`GSE290367_integrated_RNA_adata.h5ad`, `layers/counts`): 109 primary donors 30-plus (ages 30 to 91; 12 MATURE, 53 MID, 44 OLD), 7 studies (6 study dummy indicators plus intercept with `ENCODE v4 (Snyder)` reference), 15,930 retained genes (same technical-gene mask as the original), 13 `final_cell_type`s plus pooled `All_LV_Nuclei`, 62 `cell_type_plus_subcluster` states.
  - Yazar Blood (`yazar_1e44db10.h5ad`, `X`, reusing the frozen cache and dropping doublet rows labeled `peripheral blood mononuclear cell`): 932 donors 30-plus (ages 30 to 97; 59 MATURE, 193 MID, 680 OLD), 75 pools (`pool_number` fixed effects), 35,307 retained genes, 28 non-doublet `cell_type`s plus pooled `All_Cells`, 30 `predicted.celltype.l2` states.
- **Full-depth proportions and eligibility floors (unchanged):** primary and composition-matching arms use full-depth post-mask raw counts converted to within-cell proportions; post-mask `U_i >= 1500` is an eligibility filter only (never downsampled in primary). Floors: at least 20 eligible cells per donor per cell type (or state); at least 6 eligible donors per cell type for slopes plus `_prune_covariate_matrix` rank/df gates (`min_batch_donors=2`, `min_residual_df=5` Gao / `10` Yazar for batch-adjusted arms; `min_residual_df=5` for unadjusted); at least 6 eligible donors per age bin (enforced both pre- and post-singleton-batch pruning) for group comparisons and IPW reweighting. Age bins are MATURE 30 to 39, MID 40 to 59, OLD at least 60. Units below any floor are reported `UNTESTABLE`, never as nulls and never substituted with unadjusted fallbacks. Block 1 uses a donor-permutation null, not a split-half null.
- **Per-donor residual dispersion (unchanged):** exact two-pass float64 mean, unbiased variance (`ddof=1`), detection fraction per eligible `(dataset, group, donor)`; LAD (L1 median) linear fit of `log(v_g)` on `log(mu_g)` across `S_10%` (detection at least 0.10, non-constant, at least 10 genes); per-donor summary `Y_{d,c} = mean(r)` over `S_10%`.
- **Five required model arms (`results/followup-mature-block1-residual-dispersion.csv`):** 1. `primary`: donor-level OLS of `Y_{d,c}` on `age_per_decade + sex_male + batch`. 2. `within-state`: fixed MATURE-donor cell-proportion weights over fine-grained states, exact donor-clustered influence-function HC2 variance, `young_weight_coverage >= 0.50` (column name inherited from the frozen script; it is MATURE-donor coverage here). 3. `reweighted`: Hajek-stabilized truncated IPW on 3-bin (MATURE, MID, OLD) CLR composition with Pustejovsky-Tipton IPW-HC2 and Bell-McCaffrey df, plus Li and Li (2019) overlap weights and a 500-draw donor stratified pairs bootstrap. 4. `unadjusted`: OLS on `age_per_decade + sex_male`. 5. `rarefied`: 1500-UMI per-donor multinomial downsampling, primary model only. Seeds derive from `MASTER_SEED = 20261117`.
- **Gene-level secondary (`results/followup-mature-block1-gene-fractions.csv`):** batch-plus-sex-adjusted Freedman-Lane permutation (`B = 2000`, seed `20261117`) on the fraction of genes in age-blind complete-case `G_sec` whose adjusted residual dispersion is higher in OLD than MATURE, enforcing `n_mature >= 6` and `n_old >= 6` both pre- and post-pruning, alongside unadjusted raw donor-label permutation columns on the exact same row. BH correction across genes within each `(dataset, cell_type)` is NOT required and NOT applied; `q_bh_across_celltypes` summarizes multiplicity across testable non-pooled cell types per dataset.
- **BH families (mature evaluated non-pooled sizes):** Gao `primary` m=13, `within-state` m=9, `reweighted` m=4, `unadjusted` m=13, `rarefied` m=13, gene fractions m=9. Yazar `primary` m=16, `within-state` m=0, `reweighted` m=5, `unadjusted` m=17, `rarefied` m=16, gene fractions m=12. Pooled rows carry singleton m=1, q=p. `UNTESTABLE*` and `NOT_APPLICABLE*` rows are excluded from families with q=NaN.
- **MDE80 and verdict rule (unchanged):** exact two-sided noncentral-t MDE80 (`alpha = 0.05`, `power = 0.80`) standardized by the cell type's unweighted `primary` `outcome_sd`; powered-null threshold `mde80_sd_per_decade <= 0.50` (`<= 0.25` high power). Per-cell-type `PROVES` requires positive BH-significant `primary` (`q < 0.05`) and all evaluable composition-matching arms (at least 1 evaluable) positive; `DOES NOT PROVE` requires composition-artifact collapse, significantly negative `primary`, or null `primary` with `mde80_sd <= 0.50`, provided there is no nominally significant positive conflict under 2-arm composition disagreement; otherwise `INCONCLUSIVE`. Dataset-level `PROVES` requires at least 25 percent of testable non-pooled cell types to achieve `PROVES` and exceed significantly negative cell types; `DOES NOT PROVE` requires 0 `PROVES`, at least 75 percent `DOES NOT PROVE`, and under 25 percent composition disagreement; otherwise `INCONCLUSIVE`.
- **Column-name note:** the frozen column names `n_young`, `unadj_n_young`, `median_n_cells_young`, `young_weight_coverage`, and `ess_young` are inherited unchanged in the mature CSVs and refer to the MATURE (30 to 39) baseline bin throughout. The `ret_MATURE` keys in the diagnostics JSON likewise replace the original `ret_YOUNG`.
- **Independent BH check:** all `q_bh` values were recomputed from the reported p-values with independently written BH code (no statsmodels). `results/followup-mature-block1-bh-check.csv` covers 138 rows across 22 families (10 model families per dataset plus gene-fraction families, including pooled singletons): 0 mismatches, max absolute difference 9.99e-16.

## 3. Primary Results, Composition Matching, and Gene Fractions

### 3.1 Gao Heart

Evaluated non-pooled BH family sizes: `primary` m=13, `within-state` m=9, `reweighted` m=4, `unadjusted` m=13, `rarefied` m=13.

| Cell Type | Primary Slope (95% CI) | Pri p / q | N / MDE80 (SD) | Within-State Slope (p) | Reweighted Slope (p) | Unadj Slope (p) | Rarefied 1500 Slope (p) | Gene Frac OLD>M (p_2s) | Verdict |
|---|---|---|---|---|---|---|---|---|---|
| Adipocyte | -0.00216 (-0.00455 to 0.00022) | 0.07355 / 0.5263 | 35 / 0.307 | -0.00207 (0.2875) | UNTESTABLE_BIN_FLOOR | -0.00123 (0.4129) | -0.00449 (0.06989) | UNTESTABLE_BIN_FLOOR | DOES NOT PROVE (powered_null(moderate_power<=0.50)) |
| Cardiomyocyte | 0.00025 (-0.00144 to 0.00193) | 0.7726 / 0.8876 | 101 / 0.196 | -0.00035 (0.6594) | 0.00017 (0.8234) | -0.00126 (0.2524) | -0.00000 (0.9915) | 0.5657 (0.1419) | DOES NOT PROVE (powered_null(high_power<=0.25)) |
| Endocardial | 0.00046 (-0.00105 to 0.00197) | 0.5419 / 0.7828 | 66 / 0.269 | -0.00015 (0.8243) | UNTESTABLE_IMBALANCE_FAILURE | 0.00024 (0.7754) | -0.00089 (0.5336) | 0.5979 (0.02799) | DOES NOT PROVE (powered_null(moderate_power<=0.50)) |
| Endothelial | 0.00054 (-0.00084 to 0.00192) | 0.4378 / 0.7828 | 79 / 0.266 | 0.00091 (0.02009) | 0.00057 (0.4783) | -0.00002 (0.9794) | 0.00064 (0.6016) | 0.4997 (0.7756) | DOES NOT PROVE (powered_null(moderate_power<=0.50)) |
| Epicardial | 0.00082 (-0.00423 to 0.00587) | 0.7327 / 0.8876 | 17 / 0.606 | NOT_APPLICABLE_SINGLE_STATE | UNTESTABLE_BIN_FLOOR | -0.00030 (0.8784) | 0.00248 (0.6384) | UNTESTABLE_BIN_FLOOR | INCONCLUSIVE (underpowered_null(mde80_sd=0.606>0.50)) |
| Fibroblast | -0.00068 (-0.00230 to 0.00094) | 0.4073 / 0.7828 | 98 / 0.239 | -0.00025 (0.6186) | -0.00085 (0.3266) | -0.00089 (0.2897) | -0.00071 (0.2230) | 0.5662 (0.2789) | DOES NOT PROVE (powered_null(high_power<=0.25)) |
| LEC | 0.00016 (-0.00127 to 0.00160) | 0.8193 / 0.8876 | 34 / 0.432 | NOT_APPLICABLE_SINGLE_EVALUABLE_STATE | UNTESTABLE_BIN_FLOOR | 0.00018 (0.7882) | 0.00092 (0.7202) | UNTESTABLE_BIN_FLOOR | DOES NOT PROVE (powered_null(moderate_power<=0.50)) |
| Lymphoid | -0.00132 (-0.00281 to 0.00017) | 0.08097 / 0.5263 | 56 / 0.273 | -0.00058 (0.4811) | UNTESTABLE_IMBALANCE_FAILURE | -0.00122 (0.1419) | 0.00101 (0.3268) | 0.4695 (0.5187) | DOES NOT PROVE (powered_null(moderate_power<=0.50)) |
| Mast | -0.00120 (-0.00410 to 0.00171) | 0.3799 / 0.7828 | 16 / 0.795 | NOT_APPLICABLE_SINGLE_STATE | UNTESTABLE_BIN_FLOOR | -0.00101 (0.4122) | -0.00494 (0.3570) | UNTESTABLE_BIN_FLOOR | INCONCLUSIVE (underpowered_null(mde80_sd=0.795>0.50)) |
| Myeloid | 0.00068 (-0.00071 to 0.00207) | 0.3338 / 0.7828 | 99 / 0.237 | 0.00057 (0.1394) | UNTESTABLE_IMBALANCE_FAILURE | 0.00002 (0.9764) | 0.00030 (0.5821) | 0.6658 (0.002999) | DOES NOT PROVE (powered_null(high_power<=0.25)) |
| Neuronal | 0.00006 (-0.00105 to 0.00117) | 0.9101 / 0.9101 | 56 / 0.317 | NOT_APPLICABLE_SINGLE_EVALUABLE_STATE | UNTESTABLE_IMBALANCE_FAILURE | -0.00032 (0.5721) | -0.00368 (0.08768) | 0.5256 (0.4758) | DOES NOT PROVE (powered_null(moderate_power<=0.50)) |
| Pericyte | 0.00034 (-0.00063 to 0.00131) | 0.4876 / 0.7828 | 88 / 0.219 | 0.00050 (0.1183) | 0.00041 (0.3946) | -0.00016 (0.7861) | 0.00037 (0.4558) | 0.6126 (0.08496) | DOES NOT PROVE (powered_null(high_power<=0.25)) |
| vSMC | 0.00043 (-0.00070 to 0.00155) | 0.4526 / 0.7828 | 79 / 0.270 | 0.00105 (0.03600) | UNTESTABLE_IMBALANCE_FAILURE | 0.00053 (0.3608) | 0.00023 (0.8703) | 0.5253 (0.7836) | DOES NOT PROVE (powered_null(moderate_power<=0.50)) |
| **All_LV_Nuclei (pooled)** | -0.00101 (-0.00323 to 0.00120) | 0.3654 / 0.3654 | 105 / 0.167 | -0.00004 (0.9222) | UNTESTABLE_IMBALANCE_FAILURE | -0.00188 (0.2243) | -0.00122 (0.1712) | 0.6490 (0.1079) | DOES NOT PROVE (powered_null(high_power<=0.25)) |

Zero of 13 non-pooled Gao cell types pass BH (minimum primary q = 0.5263). Nominal single-arm blips do not survive multiplicity or cross-arm verification: Endothelial `within-state` p = 0.02009 (q = 0.1620), vSMC `within-state` p = 0.03600 (q = 0.1620), Adipocyte `rarefied` p = 0.06989. The `reweighted` family grows from m=2 to m=4 because Cardiomyocyte and Fibroblast now pass the median-SMD gate (Cardiomyocyte median SMD 0.27 to 0.24; Fibroblast 0.30 to 0.24); both new rows are null (Cardiomyocyte 0.00017, p = 0.8234; Fibroblast -0.00085, p = 0.3266).

### 3.2 Yazar Blood

Evaluated non-pooled BH family sizes: `primary` m=16, `within-state` m=0, `reweighted` m=5, `unadjusted` m=17, `rarefied` m=16.

| Cell Type | Primary Slope (95% CI) | Pri p / q | N / MDE80 (SD) | Within-State Slope (p) | Reweighted Slope (p) | Unadj Slope (p) | Rarefied 1500 Slope (p) | Gene Frac OLD>M (p_2s) | Verdict |
|---|---|---|---|---|---|---|---|---|---|
| CD14-low, CD16-positive monocyte | -0.00007 (-0.00104 to 0.00090) | 0.8918 / 0.9774 | 245 / 0.172 | NOT_APPLICABLE_SINGLE_STATE | UNTESTABLE_IMBALANCE_FAILURE | 0.00047 (0.2287) | 0.00058 (0.7083) | 0.5188 (0.3428) | DOES NOT PROVE (powered_null(high_power<=0.25)) |
| CD14-positive monocyte | 0.00022 (-0.00019 to 0.00062) | 0.2868 / 0.7383 | 402 / 0.110 | NOT_APPLICABLE_SINGLE_STATE | UNTESTABLE_IMBALANCE_FAILURE | 0.00022 (0.2297) | 0.00048 (0.5676) | 0.5053 (0.9225) | DOES NOT PROVE (powered_null(high_power<=0.25)) |
| CD4-positive, alpha-beta cytotoxic T cell | 0.00086 (-0.00038 to 0.00211) | 0.1703 / 0.6812 | 74 / 0.390 | NOT_APPLICABLE_SINGLE_STATE | UNTESTABLE_BIN_FLOOR | 0.00071 (0.0956) | 0.00144 (0.7124) | UNTESTABLE_BIN_FLOOR | DOES NOT PROVE (powered_null(moderate_power<=0.50)) |
| central memory CD4-positive, alpha-beta T cell | 0.00008 (-0.00010 to 0.00027) | 0.3692 / 0.7383 | 909 / 0.071 | NOT_APPLICABLE_SINGLE_STATE | 0.00007 (0.5304) | -0.00009 (0.3017) | -0.00002 (0.9359) | 0.4915 (0.4818) | DOES NOT PROVE (powered_null(high_power<=0.25)) |
| central memory CD8-positive, alpha-beta T cell | 0.00016 (-0.00112 to 0.00144) | 0.8032 / 0.9774 | 85 / 0.352 | NOT_APPLICABLE_SINGLE_STATE | UNTESTABLE_BIN_FLOOR | 0.00041 (0.3316) | -0.00720 (0.1083) | UNTESTABLE_BIN_FLOOR | DOES NOT PROVE (powered_null(moderate_power<=0.50)) |
| effector memory CD4-positive, alpha-beta T cell | 0.00001 (-0.00070 to 0.00072) | 0.9774 / 0.9774 | 250 / 0.182 | NOT_APPLICABLE_SINGLE_STATE | UNTESTABLE_ESS_COLLAPSE | 0.00011 (0.7187) | 0.00060 (0.7157) | 0.5034 (0.9955) | DOES NOT PROVE (powered_null(high_power<=0.25)) |
| effector memory CD8-positive, alpha-beta T cell | -0.00003 (-0.00032 to 0.00025) | 0.8225 / 0.9774 | 721 / 0.089 | NOT_APPLICABLE_SINGLE_STATE | -0.00006 (0.6844) | -0.00005 (0.6905) | 0.00035 (0.4084) | 0.5146 (0.4788) | DOES NOT PROVE (powered_null(high_power<=0.25)) |
| gamma-delta T cell | -0.00099 (-0.00179 to -0.00018) | 0.01665 / 0.1332 | 111 / 0.251 | NOT_APPLICABLE_SINGLE_STATE | UNTESTABLE_IMBALANCE_FAILURE | -0.00046 (0.1240) | 0.00219 (0.4059) | 0.5161 (0.3718) | DOES NOT PROVE (powered_null(moderate_power<=0.50)) |
| memory B cell | -0.00002 (-0.00090 to 0.00086) | 0.9654 / 0.9774 | 227 / 0.184 | NOT_APPLICABLE_SINGLE_STATE | UNTESTABLE_IMBALANCE_FAILURE | 0.00017 (0.6159) | 0.00066 (0.6689) | 0.4860 (0.4688) | DOES NOT PROVE (powered_null(high_power<=0.25)) |
| mucosal-associated invariant T cell | -0.00086 (-0.00264 to 0.00092) | 0.3165 / 0.7383 | 25 / 0.657 | NOT_APPLICABLE_SINGLE_STATE | UNTESTABLE_BIN_FLOOR | -0.00100 (0.03233) | -0.00403 (0.6378) | UNTESTABLE_BIN_FLOOR | INCONCLUSIVE (underpowered_null(mde80_sd=0.657>0.50)) |
| naive B cell | -0.00017 (-0.00066 to 0.00033) | 0.5054 / 0.8087 | 307 / 0.153 | NOT_APPLICABLE_SINGLE_STATE | UNTESTABLE_IMBALANCE_FAILURE | -0.00011 (0.6094) | 0.00005 (0.9415) | 0.4791 (0.2039) | DOES NOT PROVE (powered_null(high_power<=0.25)) |
| naive thymus-derived CD4-positive, alpha-beta T cell | 0.00071 (0.00046 to 0.00097) | 5.174e-08 / 8.279e-07 | 777 / 0.064 | NOT_APPLICABLE_SINGLE_STATE | 0.00080 (1.165e-06) | 0.00003 (0.8142) | 0.00008 (0.7971) | 0.5253 (0.1529) | PROVES (primary_pos_bh_sig_and_comp_agrees(reweighted)) |
| naive thymus-derived CD8-positive, alpha-beta T cell | 0.00040 (-0.00005 to 0.00085) | 0.08255 / 0.4403 | 344 / 0.134 | NOT_APPLICABLE_SINGLE_STATE | 0.00055 (0.05725) | 0.00011 (0.5486) | -0.00109 (0.2393) | 0.5211 (0.2619) | DOES NOT PROVE (powered_null(high_power<=0.25)) |
| natural killer cell | 0.00010 (-0.00016 to 0.00037) | 0.4407 / 0.7834 | 783 / 0.082 | NOT_APPLICABLE_SINGLE_EVALUABLE_STATE | 0.00007 (0.6739) | -0.00002 (0.8755) | 0.00016 (0.6399) | 0.5174 (0.4118) | DOES NOT PROVE (powered_null(high_power<=0.25)) |
| regulatory T cell | -0.00088 (-0.00283 to 0.00107) | 0.3643 / 0.7383 | 54 / 0.505 | NOT_APPLICABLE_SINGLE_STATE | UNTESTABLE_BIN_FLOOR | -0.00058 (0.3301) | 0.00515 (0.3860) | UNTESTABLE_BIN_FLOOR | INCONCLUSIVE (underpowered_null(mde80_sd=0.505>0.50)) |
| transitional stage B cell | -0.00008 (-0.00095 to 0.00078) | 0.8461 / 0.9774 | 185 / 0.219 | NOT_APPLICABLE_SINGLE_STATE | UNTESTABLE_IMBALANCE_FAILURE | -0.00016 (0.5851) | 0.00480 (0.02266) | 0.4968 (0.8376) | DOES NOT PROVE (powered_null(high_power<=0.25)) |
| **All_Cells (pooled)** | 0.00018 (-0.00012 to 0.00048) | 0.2333 / 0.2333 | 932 / 0.061 | 0.00021 (0.0006564) | 0.00015 (0.3336) | 0.00034 (0.0371) | -0.00025 (0.1688) | 0.5487 (0.002999) | DOES NOT PROVE (powered_null(high_power<=0.25)) |

Untestable non-pooled cell types in Yazar Blood (12, same set as the original): CD16-negative, CD56-bright natural killer cell, human (`UNTESTABLE_DONOR_FLOOR`, n=5); CD4-positive, alpha-beta T cell (`UNTESTABLE_DONOR_FLOOR`, n=0); CD8-positive, alpha-beta T cell (`UNTESTABLE_DONOR_FLOOR`, n=0); conventional dendritic cell (`UNTESTABLE_DONOR_FLOOR`, n=5); dendritic cell (`UNTESTABLE_DONOR_FLOOR`, n=0); double negative thymocyte (`UNTESTABLE_DONOR_FLOOR`, n=0); erythrocyte (`UNTESTABLE_DONOR_FLOOR`, n=0); hematopoietic precursor cell (`UNTESTABLE_DONOR_FLOOR`, n=0); innate lymphoid cell (`UNTESTABLE_DONOR_FLOOR`, n=0); plasmablast (`UNTESTABLE_NOT_EVALUABLE_WITHIN_BATCH_DF`, n=14); plasmacytoid dendritic cell (`UNTESTABLE_DONOR_FLOOR`, n=3); platelet (`UNTESTABLE_DONOR_FLOOR`, n=0).

The `reweighted` family grows from m=4 to m=5 because naive CD8 T now passes the median-SMD gate (0.26 to 0.20); its reweighted slope is 0.00055 (p = 0.05725, q = 0.1431), still null. Gamma-delta T has the only nominally significant negative primary in Yazar (-0.00099, p = 0.01665, q = 0.1332); it is DOES NOT PROVE by powered null, not by opposite direction, because the BH q is not significant.

## 4. Gene-Level Secondary: OLD vs MATURE Residual Dispersion Fraction (Batch-plus-Sex Freedman-Lane vs Unadjusted)

Source: `results/followup-mature-block1-gene-fractions.csv`. BH across genes is NOT required and NOT applied; `perm_p_two_sided` is the gene-level fraction test, with `q_bh_across_celltypes` summarizing multiplicity across testable non-pooled cell types per dataset. The `n_young` and `n_old` columns are inherited names and count MATURE and OLD donors. Counts are post-prune where the row passed the pre-prune floor and pre-prune where it failed early (for example Yazar central memory CD8 T shows 4/80 pre-prune at mature, while the original showed 4/65 post-prune).

| Dataset | Cell Type | Adj N_M / N_O | N Genes | Adj Frac OLD>M (null mean) | Adj p_2s / q_ct (orig p_2s) | Unadj Frac OLD>M (p_2s) |
|---|---|---|---|---|---|---|
| gao | Adipocyte | 3 / 12 | 0 | UNTESTABLE_BIN_FLOOR | NA (orig UNTESTABLE) | NA |
| gao | Cardiomyocyte | 10 / 38 | 6302 | 0.5657 (0.5072) | 0.1419 / 0.3193 (orig 0.1349) | 0.6038 (0.0560) |
| gao | Endocardial | 7 / 28 | 3895 | 0.5979 (0.5072) | 0.02799 / 0.1259 (orig 0.0390) | 0.6135 (0.0110) |
| gao | Endothelial | 9 / 29 | 3684 | 0.4997 (0.5072) | 0.7756 / 0.7836 (orig 0.8876) | 0.5337 (0.5567) |
| gao | Epicardial | 2 / 8 | 0 | UNTESTABLE_BIN_FLOOR | NA (orig UNTESTABLE) | NA |
| gao | Fibroblast | 10 / 36 | 4859 | 0.5662 (0.5072) | 0.2789 / 0.5019 (orig 0.3568) | 0.5477 (0.5747) |
| gao | LEC | 4 / 13 | 0 | UNTESTABLE_BIN_FLOOR | NA (orig 0.5416, p 0.2719) | NA |
| gao | Lymphoid | 8 / 21 | 3033 | 0.4695 (0.5072) | 0.5187 / 0.6670 (orig 0.5037) | 0.4856 (0.6877) |
| gao | Mast | 1 / 8 | 0 | UNTESTABLE_BIN_FLOOR | NA (orig UNTESTABLE) | NA |
| gao | Myeloid | 10 / 38 | 3842 | 0.6658 (0.5072) | 0.002999 / 0.02699 (orig 0.0010) | 0.6204 (0.0980) |
| gao | Neuronal | 8 / 18 | 3263 | 0.5256 (0.5072) | 0.4758 / 0.6670 (orig 0.5237) | 0.5451 (0.2399) |
| gao | Pericyte | 10 / 32 | 3500 | 0.6126 (0.5072) | 0.08496 / 0.2549 (orig 0.0510) | 0.6151 (0.0690) |
| gao | vSMC | 8 / 29 | 2766 | 0.5253 (0.5072) | 0.7836 / 0.7836 (orig 0.3888) | 0.5337 (0.5567) |
| gao | All_LV_Nuclei | 11 / 40 | 5327 | 0.6490 (0.5048) | 0.1079 / 0.1079 (orig 0.4238) | 0.7152 (0.0500) |
| yazar | CD14-low, CD16-positive monocyte | 8 / 194 | 1997 | 0.5188 (0.5019) | 0.3428 / 0.6423 (orig 0.7156) | 0.4787 (0.1909) |
| yazar | CD14-positive monocyte | 20 / 294 | 1615 | 0.5053 (0.5019) | 0.9225 / 0.9955 (orig 0.9105) | 0.4777 (0.2519) |
| yazar | CD16-negative, CD56-bright natural killer cell, human | 0 / 4 | 0 | UNTESTABLE_BIN_FLOOR | NA (orig UNTESTABLE) | NA |
| yazar | CD4-positive, alpha-beta T cell | 0 / 0 | 0 | UNTESTABLE_BIN_FLOOR | NA (orig UNTESTABLE) | NA |
| yazar | CD4-positive, alpha-beta cytotoxic T cell | 1 / 78 | 0 | UNTESTABLE_BIN_FLOOR | NA (orig UNTESTABLE) | NA |
| yazar | CD8-positive, alpha-beta T cell | 0 / 0 | 0 | UNTESTABLE_BIN_FLOOR | NA (orig UNTESTABLE) | NA |
| yazar | central memory CD4-positive, alpha-beta T cell | 58 / 660 | 1770 | 0.4915 (0.5019) | 0.4818 / 0.6423 (orig 0.7396) | 0.4565 (0.0050) |
| yazar | central memory CD8-positive, alpha-beta T cell | 4 / 80 | 0 | UNTESTABLE_BIN_FLOOR | NA (orig UNTESTABLE_POST_PRUNE_BIN_FLOOR) | 0.5366 (0.2299) |
| yazar | conventional dendritic cell | 0 / 5 | 0 | UNTESTABLE_BIN_FLOOR | NA (orig UNTESTABLE) | NA |
| yazar | dendritic cell | 0 / 0 | 0 | UNTESTABLE_BIN_FLOOR | NA (orig UNTESTABLE) | NA |
| yazar | double negative thymocyte | 0 / 0 | 0 | UNTESTABLE_BIN_FLOOR | NA (orig UNTESTABLE) | NA |
| yazar | effector memory CD4-positive, alpha-beta T cell | 6 / 190 | 1315 | 0.5034 (0.5019) | 0.9955 / 0.9955 (orig 0.6667) | 0.5190 (0.5037) |
| yazar | effector memory CD8-positive, alpha-beta T cell | 33 / 552 | 1335 | 0.5146 (0.5019) | 0.4788 / 0.6423 (orig 0.09295) | 0.5206 (0.2959) |
| yazar | erythrocyte | 0 / 0 | 0 | UNTESTABLE_BIN_FLOOR | NA (orig UNTESTABLE) | NA |
| yazar | gamma-delta T cell | 7 / 76 | 1736 | 0.5161 (0.5019) | 0.3718 / 0.6423 (orig 0.8956) | 0.5054 (0.9555) |
| yazar | hematopoietic precursor cell | 0 / 0 | 0 | UNTESTABLE_BIN_FLOOR | NA (orig UNTESTABLE) | NA |
| yazar | innate lymphoid cell | 0 / 0 | 0 | UNTESTABLE_BIN_FLOOR | NA (orig UNTESTABLE) | NA |
| yazar | memory B cell | 9 / 167 | 1428 | 0.4860 (0.5019) | 0.4688 / 0.6423 (orig 0.6797) | 0.4927 (0.6517) |
| yazar | mucosal-associated invariant T cell | 3 / 23 | 0 | UNTESTABLE_BIN_FLOOR | NA (orig UNTESTABLE_NOT_EVALUABLE_WITHIN_BATCH_DF) | 0.4820 (0.2499) |
| yazar | naive B cell | 11 / 246 | 1223 | 0.4791 (0.5019) | 0.2039 / 0.6423 (orig 0.1439) | 0.4946 (0.5277) |
| yazar | naive thymus-derived CD4-positive, alpha-beta T cell | 50 / 564 | 1283 | 0.5253 (0.5019) | 0.1529 / 0.6423 (orig 0.1749) | 0.5116 (0.6137) |
| yazar | naive thymus-derived CD8-positive, alpha-beta T cell | 35 / 205 | 1324 | 0.5211 (0.5019) | 0.2619 / 0.6423 (orig 0.5537) | 0.4778 (0.1799) |
| yazar | natural killer cell | 38 / 598 | 1492 | 0.5174 (0.5019) | 0.4118 / 0.6423 (orig 0.4108) | 0.5209 (0.2859) |
| yazar | plasmablast | 1 / 8 | 0 | UNTESTABLE_BIN_FLOOR | NA (orig UNTESTABLE) | NA |
| yazar | plasmacytoid dendritic cell | 0 / 2 | 0 | UNTESTABLE_BIN_FLOOR | NA (orig UNTESTABLE) | NA |
| yazar | platelet | 0 / 0 | 0 | UNTESTABLE_BIN_FLOOR | NA (orig UNTESTABLE) | NA |
| yazar | regulatory T cell | 0 / 59 | 0 | UNTESTABLE_BIN_FLOOR | NA (orig UNTESTABLE) | NA |
| yazar | transitional stage B cell | 8 / 139 | 1417 | 0.4968 (0.5019) | 0.8376 / 0.9955 (orig 0.8456) | 0.5088 (0.7946) |
| yazar | All_Cells | 59 / 680 | 3042 | 0.5487 (0.5019) | 0.002999 / 0.002999 (orig 0.0140) | 0.4722 (0.0970) |

Two gene-fraction signals survive the mature restriction with BH significance across cell types: Gao Myeloid (frac 0.6658, p = 0.002999, q_ct = 0.02699; original 0.6885, p = 0.0010, q_ct = 0.0100) and pooled Yazar All_Cells (frac 0.5487, p = 0.002999, q_ct = 0.002999; original 0.5432, p = 0.0140). Gao LEC drops from EVALUATED to UNTESTABLE_BIN_FLOOR (MATURE n = 4). The `perm_null_mean` sits slightly above 0.50 (0.5019 to 0.5072 across all evaluated rows) for the same finite-sample skew reason as the original; the Freedman-Lane permutation distribution calibrates for it.

## 5. Mandatory Non-Primary Diagnostics (Outside Primary CSVs and BH Families)

### 5.1 Methodological Cross-Checks (N-plus-Depth Adjusted, Shared MATURE Trend, Gamma-GLM, OLS-All)

| Cell Type | Pri LAD (p) | +N+D LAD (p) | WS +N+D (p) | Shared Trend (p) | Shared +N+D (p) | Gamma-GLM (p) | Gamma +N+D (p) | OLS-All (p) | OLS-All +N+D (p) |
|---|---|---|---|---|---|---|---|---|---|
| Gao Adipocyte | -0.00216 (0.07355) | -0.00042 (0.6396) | 0.00056 (0.7150) | 0.03603 (0.2378) | -0.00706 (0.2581) | 0.00334 (0.5161) | 0.00079 (0.7567) | -0.00001 (0.9962) | 0.00036 (0.7774) |
| Gao Cardiomyocyte | 0.00025 (0.7726) | 0.00073 (0.2323) | 0.00028 (0.5766) | 0.01857 (0.2951) | 0.00107 (0.8000) | -0.00094 (0.5853) | -0.00187 (0.1354) | -0.00030 (0.8100) | -0.00068 (0.1451) |
| Gao Endocardial | 0.00046 (0.5419) | 0.00071 (0.2811) | 0.00013 (0.8334) | 0.02197 (0.1454) | 0.00094 (0.8301) | 0.00156 (0.5080) | -0.00032 (0.7956) | 0.00151 (0.3846) | 0.00102 (0.4286) |
| Gao Endothelial | 0.00054 (0.4378) | 0.00063 (0.3149) | 0.00079 (0.03563) | -0.00999 (0.1646) | 0.00108 (0.7335) | -0.00215 (0.2056) | -0.00054 (0.5022) | -0.00147 (0.5028) | 0.00049 (0.6872) |
| Gao Epicardial | 0.00082 (0.7327) | -0.00273 (0.05445) | NA (NA) | -0.01128 (0.6588) | -0.01077 (0.3210) | 0.00875 (0.08632) | 0.00468 (0.0891) | 0.00380 (0.1872) | -0.00104 (0.6304) |
| Gao Fibroblast | -0.00068 (0.4073) | -0.00048 (0.5298) | -0.00002 (0.9645) | 0.00999 (0.3252) | -0.00158 (0.6466) | 0.00112 (0.3189) | 0.00159 (0.1341) | 0.00127 (0.2601) | 0.00032 (0.4924) |
| Gao LEC | 0.00016 (0.8193) | 0.00030 (0.6797) | NA (NA) | 0.00404 (0.8520) | -0.00696 (0.07556) | 0.00218 (0.4514) | 0.00218 (0.1300) | 0.00019 (0.9157) | -0.00013 (0.9343) |
| Gao Lymphoid | -0.00132 (0.08097) | -0.00138 (0.09358) | -0.00031 (0.7729) | -0.00015 (0.9863) | -0.00531 (0.3546) | 0.00548 (0.007079) | 0.00277 (0.05624) | 0.00333 (0.2043) | 0.00016 (0.9147) |
| Gao Mast | -0.00120 (0.3799) | 0.00034 (0.8584) | NA (NA) | NA (NA) | NA (NA) | 0.00535 (0.1359) | 0.00425 (0.2604) | -0.00105 (0.5543) | 0.00027 (0.9196) |
| Gao Myeloid | 0.00068 (0.3338) | 0.00112 (0.06974) | 0.00057 (0.09756) | 0.00885 (0.2744) | 0.00109 (0.7193) | -0.00161 (0.09624) | -0.00107 (0.2107) | -0.00058 (0.6567) | -0.00006 (0.9315) |
| Gao Neuronal | 0.00006 (0.9101) | 0.00036 (0.4370) | NA (NA) | -0.00168 (0.8389) | 0.00014 (0.9736) | -0.00128 (0.4301) | -0.00076 (0.4038) | 0.00383 (0.4918) | -0.00206 (0.5513) |
| Gao Pericyte | 0.00034 (0.4876) | -0.00026 (0.5861) | 0.00029 (0.2817) | -0.01482 (0.03215) | -0.00355 (0.3429) | -0.00093 (0.3911) | 0.00047 (0.5227) | -0.00007 (0.9692) | 0.00100 (0.1409) |
| Gao vSMC | 0.00043 (0.4526) | 0.00004 (0.9409) | 0.00050 (0.2502) | -0.00693 (0.3565) | -0.00232 (0.5054) | -0.00038 (0.7864) | 0.00001 (0.9910) | 0.00134 (0.4427) | -0.00039 (0.6950) |
| Gao All_LV_Nuclei | -0.00101 (0.3654) | -0.00115 (0.2916) | 0.00026 (0.3790) | -0.01243 (0.1517) | -0.00649 (0.1965) | 0.00108 (0.6053) | 0.00170 (0.3560) | -0.00185 (0.3286) | -0.00069 (0.5906) |
| Yazar CD14-low, CD16-positive monocyte | -0.00007 (0.8918) | 0.00005 (0.9164) | NA (NA) | 0.00126 (0.8041) | -0.00140 (0.02312) | 0.00027 (0.7429) | 0.00030 (0.4311) | -0.00110 (0.04654) | -0.00053 (0.2576) |
| Yazar CD14-positive monocyte | 0.00022 (0.2868) | 0.00014 (0.4707) | NA (NA) | 0.00757 (0.02124) | 0.00031 (0.3964) | 0.00113 (0.03798) | -0.00013 (0.7042) | -0.00053 (0.2257) | -0.00013 (0.7282) |
| Yazar CD4-positive, alpha-beta cytotoxic T cell | 0.00086 (0.1703) | 0.00104 (0.09176) | NA (NA) | NA (NA) | NA (NA) | 0.00077 (0.6161) | -0.00073 (0.2132) | -0.00037 (0.8549) | -0.00131 (0.2671) |
| Yazar central memory CD4-positive, alpha-beta T cell | 0.00008 (0.3692) | 0.00007 (0.4074) | NA (NA) | 0.00172 (0.1257) | 0.00032 (0.07765) | 0.00007 (0.8279) | 0.00011 (0.6463) | -0.00009 (0.7741) | -0.00008 (0.5286) |
| Yazar central memory CD8-positive, alpha-beta T cell | 0.00016 (0.8032) | -0.00025 (0.6876) | NA (NA) | -0.00289 (0.7440) | -0.00196 (0.1161) | -0.00011 (0.9256) | -0.00109 (0.0417) | 0.00071 (0.5849) | -0.00015 (0.9009) |
| Yazar effector memory CD4-positive, alpha-beta T cell | 0.00001 (0.9774) | 0.00033 (0.2222) | NA (NA) | -0.00264 (0.5139) | 0.00087 (0.1613) | -0.00139 (0.0138) | -0.00010 (0.6800) | 0.00001 (0.9878) | 0.00016 (0.8034) |
| Yazar effector memory CD8-positive, alpha-beta T cell | -0.00003 (0.8225) | 0.00002 (0.8986) | NA (NA) | 0.00076 (0.6263) | -0.00041 (0.09208) | 0.00046 (0.2731) | -0.00004 (0.8771) | 0.00001 (0.9828) | -0.00031 (0.1377) |
| Yazar gamma-delta T cell | -0.00099 (0.01665) | -0.00085 (0.0328) | NA (NA) | 0.00048 (0.9221) | -0.00136 (0.09909) | 0.00078 (0.4084) | -0.00025 (0.5227) | -0.00030 (0.7668) | -0.00004 (0.9572) |
| Yazar memory B cell | -0.00002 (0.9654) | 0.00019 (0.6086) | NA (NA) | -0.00554 (0.2709) | 0.00108 (0.1308) | -0.00141 (0.05956) | -0.00100 (0.04003) | -0.00012 (0.9182) | 0.00048 (0.5956) |
| Yazar mucosal-associated invariant T cell | -0.00086 (0.3165) | -0.00073 (0.4348) | NA (NA) | -0.02809 (0.3121) | -0.00374 (0.3191) | -0.00457 (0.06495) | -0.00029 (0.7937) | -0.00215 (0.4337) | -0.00183 (0.6437) |
| Yazar naive B cell | -0.00017 (0.5054) | -0.00020 (0.3987) | NA (NA) | 0.00186 (0.5577) | 0.00028 (0.5323) | 0.00009 (0.9068) | -0.00082 (0.2004) | 0.00134 (0.0819) | 0.00027 (0.5000) |
| Yazar naive thymus-derived CD4-positive, alpha-beta T cell | 0.00071 (5.174e-08) | 0.00024 (0.01893) | NA (NA) | 0.00046 (0.7370) | -0.00002 (0.9169) | -0.00081 (0.01073) | 0.00015 (0.5448) | -0.00085 (0.004412) | -0.00011 (0.4467) |
| Yazar naive thymus-derived CD8-positive, alpha-beta T cell | 0.00040 (0.08255) | 0.00006 (0.8179) | NA (NA) | -0.01028 (0.0003628) | 0.00007 (0.8551) | -0.00347 (5.843e-08) | -0.00033 (0.5673) | -0.00133 (0.02746) | -0.00016 (0.7232) |
| Yazar natural killer cell | 0.00010 (0.4407) | 0.00004 (0.7690) | NA (NA) | 0.00373 (0.0165) | -0.00002 (0.9456) | 0.00113 (0.04556) | 0.00024 (0.6585) | 0.00136 (0.001456) | -0.00022 (0.4196) |
| Yazar regulatory T cell | -0.00088 (0.3643) | -0.00043 (0.6359) | NA (NA) | NA (NA) | NA (NA) | 0.00173 (0.1186) | 0.00109 (0.09296) | 0.00266 (0.3970) | 0.00359 (0.1509) |
| Yazar transitional stage B cell | -0.00008 (0.8461) | 0.00004 (0.9210) | NA (NA) | 0.00432 (0.4215) | 0.00322 (0.001887) | 0.00107 (0.3400) | 0.00055 (0.5977) | -0.00015 (0.8882) | -0.00006 (0.9438) |
| Yazar All_Cells | 0.00018 (0.2333) | 0.00024 (0.08144) | 0.00009 (0.09397) | 0.00021 (0.8094) | 0.00134 (4.975e-05) | 0.00358 (0.01083) | 0.00239 (0.08194) | 0.00146 (0.01193) | 0.00075 (0.007494) |

Changes against the original cross-checks: Gao `Shared MATURE Trend + N+Depth` is null everywhere (p = 0.0756 to 0.9736; original 0.2093 to 0.9180). Yazar naive CD4 `+N+Depth` stays nominally significant (0.00024, p = 0.01893; original 0.0001, p = 0.3399). Gamma-GLM `+N+Depth` now has 2 of 16 non-pooled Yazar types at p under 0.05 (central memory CD8 T p = 0.0417, memory B p = 0.0400; original 0 of 16), both with negative slopes; OLS-All `+N+Depth` still has 0 of 16 (original 0 of 16). Central memory CD4 `Shared +N+Depth` is now p = 0.07765 (original 0.0140, q = 0.2107); at mature its BH q across the 14 evaluated non-pooled types is 0.2616. Pooled All_Cells `Shared +N+Depth` stays nominally positive (p = 4.975e-05; original 2.09e-05) but collapses after adjusting for `+N+Depth` plus top-15 state fractions (p = 0.7079; original 0.6351), and `+N+Depth` alone collapses the All_Cells gene fraction from 0.5487 to 0.5266.

### 5.2 Nuisance Covariates, Trend Parameters, and Batch-Partial Correlations with Y

Gao Heart: no `log(N)` slope reaches p under 0.02 except Lymphoid (0.2032, p = 0.02356) and Epicardial (0.4067, p = 0.04917, n = 17). Centered trend height `h_hat` is nominally negative in Mast (p = 0.02397) and Pericyte (p = 0.03056) and null elsewhere. Full Gao rows: Adipocyte logN -0.0491 (p = 0.69), log(E[1/U]) 0.0640 (p = 0.07558), Var(log U) -0.0053 (p = 0.6928), trend slope -0.00502 (p = 0.5561), height 0.04210 (p = 0.1888), pcorr(Y, logN) 0.493; Cardiomyocyte 0.0359 (0.4553), 0.0228 (0.3156), -0.0029 (0.8259), -0.00081 (0.8083), 0.01937 (0.3097), -0.233; Endocardial -0.0252 (0.7917), 0.0285 (0.1582), -0.0176 (0.03739), -0.00515 (0.3463), 0.02185 (0.1552), 0.471; Endothelial -0.0786 (0.4926), -0.0104 (0.2531), 0.0052 (0.3357), 0.00133 (0.6634), -0.01120 (0.1242), 0.315; Epicardial 0.4067 (0.04917), -0.0348 (0.3918), 0.0054 (0.5473), -0.01627 (0.09301), -0.01031 (0.7127), 0.556; Fibroblast -0.0964 (0.1008), 0.0167 (0.2553), -0.0134 (0.02927), -0.00091 (0.7355), 0.01203 (0.2648), 0.348; LEC -0.0434 (0.6402), 0.0164 (0.5589), -0.0119 (0.1546), -0.00762 (0.1785), 0.00488 (0.8255), 0.391; Lymphoid 0.2032 (0.02356), -0.0041 (0.5999), 0.0011 (0.7621), -0.01464 (0.02105), 0.00354 (0.6887), -0.068; Mast 0.1423 (0.1343), -0.0432 (0.03693), 0.0135 (0.1099), -0.00593 (0.3827), -0.03989 (0.02397), -0.257; Myeloid -0.1263 (0.05279), 0.0130 (0.2845), -0.0079 (0.09475), -0.00018 (0.9544), 0.00820 (0.3471), 0.368; Neuronal -0.0281 (0.6160), -0.0004 (0.9671), -0.0098 (0.2362), 0.00004 (0.9942), -0.00217 (0.7877), 0.586; Pericyte -0.0650 (0.3884), -0.0136 (0.1517), -0.0001 (0.9673), -0.00164 (0.6946), -0.01557 (0.03056), -0.144; vSMC 0.0212 (0.7579), -0.0106 (0.4127), -0.0031 (0.5136), -0.00422 (0.2726), -0.00667 (0.3838), 0.477; All_LV_Nuclei -0.0235 (0.4402), -0.0060 (0.5248), -0.0032 (0.7774), -0.00312 (0.3633), -0.01228 (0.1650), 0.167.

Yazar Blood: naive CD4 T `log(N)` falls steeply with age (-0.0897 per decade, p = 8.40e-08; original -0.0696, p = 2.90e-07) and Y correlates strongly negatively with `log(N)` (pcorr = -0.625; original -0.611). Naive CD8 T `log(N)` falls faster (-0.2097, p = 6.90e-16; original -0.1965, p = 1.25e-22). Natural killer cell `log(N)` rises (0.0777, p = 1.37e-06; original 0.0692, p = 1.97e-07). Full Yazar rows: CD14-low, CD16-positive monocyte logN 0.0213 (p = 0.3666), trend slope 0.00014 (p = 0.9266), height 0.00108 (p = 0.8385), pcorr -0.238; CD14-positive monocyte 0.0585 (0.002198), -0.00157 (0.0449), 0.00771 (0.0167), 0.362; CD4-positive, alpha-beta cytotoxic T cell 0.0346 (0.5913), 0.00077 (0.7749), 0.01760 (0.04563), 0.173; central memory CD4 T -0.0090 (0.5210), 0.00052 (0.1502), 0.00143 (0.2038), -0.331; central memory CD8 T 0.0344 (0.4647), 0.00042 (0.8275), -0.00358 (0.6873), 0.141; effector memory CD4 T -0.0607 (0.009074), 0.00226 (0.02237), -0.00323 (0.4077), 0.390; effector memory CD8 T 0.0352 (0.1156), -0.00036 (0.5162), 0.00065 (0.6791), -0.200; gamma-delta T 0.0558 (0.2001), -0.00079 (0.5554), 0.00152 (0.7569), -0.069; memory B -0.0036 (0.8850), 0.00294 (0.01458), -0.00583 (0.2272), 0.255; mucosal-associated invariant T -0.1408 (0.06103), -0.00076 (0.8663), -0.02727 (0.3410), 0.110; naive B 0.0623 (0.03225), -0.00092 (0.1950), 0.00224 (0.4661), -0.036; naive CD4 T -0.0897 (8.40e-08), 0.00144 (5.70e-05), -0.00068 (0.6285), -0.625; naive CD8 T -0.2097 (6.90e-16), 0.00312 (1.125e-06), -0.01121 (7.875e-05), -0.210; natural killer cell 0.0777 (1.37e-06), -0.00184 (3.21e-05), 0.00406 (0.008828), -0.016; regulatory T 0.0362 (0.3456), -0.00768 (0.006483), -0.00312 (0.7569), 0.020; transitional stage B cell 0.0835 (0.06857), 0.00050 (0.7343), 0.00417 (0.4357), 0.006; All_Cells -0.0229 (0.01379), 0.00172 (6.21e-07), -0.00043 (0.6370), 0.166.

### 5.3 Composition-Matching Audit: IPW Balance, Li and Li (2019) Overlap Weights

Gao Heart overlap weights are null across all 10 evaluable rows (p = 0.0897 to 0.9590): Cardiomyocyte 0.00016 (p = 0.8433), Endocardial 0.00041 (p = 0.5674), Endothelial 0.00050 (p = 0.5385), Fibroblast -0.00084 (p = 0.3292), Lymphoid -0.00103 (p = 0.0897), Myeloid 0.00079 (p = 0.2518), Neuronal 0.00003 (p = 0.9590), Pericyte 0.00036 (p = 0.4538), vSMC 0.00039 (p = 0.5092), All_LV_Nuclei -0.00088 (p = 0.4284). Adipocyte and LEC overlap p-values are suppressed by the same ESS/leverage floors as primary IPW (Adipocyte fails the MATURE bin floor; LEC suppressed at minESS under 6.0). IPW `reweighted` rows: Cardiomyocyte ESS 99.3/10.9, maxSMD 0.79 to 0.69, median SMD 0.24; Endothelial ESS 78.1/9.8, maxSMD 0.79 to 0.63, median SMD 0.19; Fibroblast ESS 96.3/10.9, maxSMD 0.75 to 0.62, median SMD 0.24; Pericyte ESS 86.9/10.8, maxSMD 0.66 to 0.52, median SMD 0.18. All four carry RESIDUAL_IMBALANCE_WARNING on max SMD.

Yazar Blood overlap weights: 12 evaluable rows, 1 NA (effector memory CD4 T, suppressed). Naive CD4 T overlap slope is now nominally significant (0.00054, p = 0.01236; original 0.0005, p = 0.0545), and gamma-delta T is nominally significant (-0.00079, p = 0.04323); all other evaluable rows are null (p = 0.0966 to 0.9668), including naive CD8 T (0.00053, p = 0.09656), natural killer cell (0.00015, p = 0.4569), and All_Cells (-0.00001, p = 0.9668). IPW `reweighted` rows: central memory CD4 T ESS 677.6/25.8, median SMD 0.23; effector memory CD8 T ESS 612.9/18.3, median SMD 0.19; naive CD4 T ESS 595.8/27.9, maxSMD 1.48 to 1.27, median SMD 0.21; naive CD8 T ESS 283.3/26.4, median SMD 0.20; natural killer cell ESS 761.9/30.7, maxSMD 0.72 to 0.17, median SMD 0.17 (the only ADEQUATELY_BALANCED row); All_Cells ESS 714.8/26.2, median SMD 0.24.

### 5.4 Post-Mask at-least-1500 UMI Retention Fraction by Age Bin (MATURE / MID / OLD)

Gao Heart: Adipocyte 0.71/0.86/0.53; Cardiomyocyte 0.98/0.92/0.93; Endocardial 0.87/0.50/0.66; Endothelial 0.51/0.28/0.41; Epicardial 0.90/0.87/0.83; Fibroblast 0.73/0.54/0.48; LEC 0.79/0.52/0.43; Lymphoid 0.20/0.19/0.22; Mast 0.16/0.13/0.32; Myeloid 0.49/0.46/0.41; Neuronal 0.46/0.36/0.36; Pericyte 0.40/0.24/0.35; vSMC 0.48/0.32/0.40; All_LV_Nuclei 0.66/0.58/0.55.

Yazar Blood: CD14-low, CD16-positive monocyte 0.93/0.91/0.92; CD14-positive monocyte 0.75/0.75/0.79; CD4-positive, alpha-beta cytotoxic T cell 0.29/0.39/0.43; central memory CD4-positive, alpha-beta T cell 0.37/0.46/0.53; central memory CD8-positive, alpha-beta T cell 0.45/0.52/0.62; effector memory CD4-positive, alpha-beta T cell 0.32/0.42/0.50; effector memory CD8-positive, alpha-beta T cell 0.35/0.42/0.51; gamma-delta T cell 0.39/0.44/0.55; memory B cell 0.40/0.47/0.55; mucosal-associated invariant T cell 0.39/0.54/0.61; naive B cell 0.23/0.27/0.37; naive thymus-derived CD4-positive, alpha-beta T cell 0.31/0.38/0.49; naive thymus-derived CD8-positive, alpha-beta T cell 0.41/0.49/0.58; natural killer cell 0.38/0.48/0.55; regulatory T cell 0.16/0.25/0.32; transitional stage B cell 0.43/0.49/0.56; All_Cells 0.37/0.45/0.53.

### 5.5 Real-Cell Finite-N Calibration Bias (Subsampled Minus Full-N Within High-N Donors, 5 Seeds)

| Cell Type | High-N Donors | Subsample N | LAD Mean Bias | LAD Slope beta Bias | OLS-All Mean Bias | Gamma-GLM Mean Bias |
|---|---|---|---|---|---|---|
| Gao Cardiomyocyte | 99 | 20 | 0.0022 | 0.1311 | 0.0069 | -0.1091 |
| Gao Cardiomyocyte | 99 | 50 | 0.0059 | 0.0849 | -0.0091 | -0.0664 |
| Gao Cardiomyocyte | 99 | 100 | 0.0039 | 0.0575 | -0.0116 | -0.0407 |
| Gao Cardiomyocyte | 99 | 200 | 0.0012 | 0.0355 | -0.0115 | -0.0224 |
| Yazar central memory CD4-positive, alpha-beta T cell | 148 | 20 | -0.0013 | 0.0403 | 0.0047 | -0.0393 |
| Yazar central memory CD4-positive, alpha-beta T cell | 148 | 50 | 0.0041 | 0.0316 | -0.0111 | -0.0188 |
| Yazar central memory CD4-positive, alpha-beta T cell | 148 | 100 | 0.0030 | 0.0203 | -0.0087 | -0.0085 |
| Yazar central memory CD4-positive, alpha-beta T cell | 148 | 200 | 0.0007 | 0.0069 | -0.0032 | -0.0021 |
| Yazar naive thymus-derived CD4-positive, alpha-beta T cell | 127 | 20 | 0.0081 | 0.0222 | 0.0003 | -0.0356 |
| Yazar naive thymus-derived CD4-positive, alpha-beta T cell | 127 | 50 | 0.0112 | 0.0284 | -0.0119 | -0.0176 |
| Yazar naive thymus-derived CD4-positive, alpha-beta T cell | 127 | 100 | 0.0071 | 0.0197 | -0.0083 | -0.0076 |
| Yazar naive thymus-derived CD4-positive, alpha-beta T cell | 127 | 200 | 0.0016 | 0.0063 | -0.0031 | -0.0015 |

The naive CD4 T calibration bias at N=50 (0.0112) is again nearly 3 times the bias in central memory CD4 T (0.0041), same as the original. High-N donor counts fall slightly at 30-plus (Gao Cardiomyocyte 106 to 99; Yazar central memory CD4 T 152 to 148; naive CD4 T 131 to 127).

### 5.6 Yazar Coarse-Lineage Diagnostic (Long et al. 2025 Reconciliation)

| Coarse Lineage | L2 States (eval / tot) | N Donors | Unmatched LAD Slope (p) | Within-L2 LAD Slope (p; MATURE cov) | Within-L2 +N+D Slope (p) | Unmatched Gamma-GLM Slope (p) |
|---|---|---|---|---|---|---|
| B | 3 / 4 | 708 | 0.00031 (0.4508) | -0.00009 (0.6653; 0.88) | 0.00001 (0.9789) | 0.00211 (0.01832) |
| CD4_T | 5 / 6 | 929 | 0.00040 (2.82e-05) | 0.00037 (1.50e-05; 1.00) | 0.00017 (0.01986) | -0.00012 (0.7566) |
| CD8_T | 3 / 4 | 856 | -0.00027 (0.02493) | 0.00018 (0.2180; 1.00) | 0.00002 (0.9137) | 0.00016 (0.6674) |
| Monocyte | 2 / 2 | 552 | 0.00036 (0.07798) | 0.00013 (0.5111; 1.00) | 0.00011 (0.5690) | 0.00104 (0.03842) |
| NK | 1 / 3 | 805 | 0.00004 (0.7269) | NA (NA; NA) | NA (NA) | 0.00088 (0.1064) |

Changes against the original: Monocyte unmatched is now null (0.00036, p = 0.0780; original 0.0004, p = 0.0196), so the pooling artifact it illustrated no longer even appears unmatched. CD4_T within-L2 plus N-plus-Depth stays nominally significant (0.00017, p = 0.01986; original 0.0001, p = 0.1604), so the original collapse under N-plus-Depth does not reproduce at 30-plus. CD8_T unmatched is now nominally significant negative (-0.00027, p = 0.02493; original -0.0002, p = 0.0527). The leave-one-state-out proof still holds: dropping only CD4 Naive collapses All_Cells within-state from 0.00021 (p = 0.0006564) to 0.00007 (p = 0.3219), and dropping CD4 Naive plus CD8 Naive yields 0.00001 (p = 0.8337); dropping CD4 Naive collapses coarse CD4_T within-L2 from 0.00037 (p = 1.50e-05) to 0.00006 (p = 0.5294). CD4 Naive carries a renormalized MATURE weight of 22.2 percent in All_Cells (raw MATURE share 21.2 percent; CD8 Naive adds 10.8 percent) and 46.6 percent of MATURE weight in coarse CD4_T.

### 5.7 Naive CD4 T Mechanistic and Causal Audit (Mature)

This is the one place where the mature rerun changes the original's mechanistic conclusion. The frozen PROVES rule still passes, and the signal is stronger: primary 0.00071 (p = 5.174e-08, q = 8.279e-07; original 0.00046, p = 3.00e-05, q = 4.80e-04) and reweighted 0.00080 (p = 1.165e-06; original 0.00068, p = 2.78e-05). But the finite-N calibration audit no longer kills it. Non-parametrically subtracting the interpolated real-cell finite-N calibration curve accounts for 56.4 percent of the primary slope (induced slope 0.000401 per decade, p = 8.40e-07; original 66.4 percent, 0.000302 per decade, p = 4.32e-06) and leaves a still-significant debiased age slope of 0.000313 (p = 0.000638; original 0.000153, p = 0.0469) and 0.000295 after depth adjustment (p = 0.00106; original 0.000131, p = 0.0818, null). Pure capture-yield sensitivity is still eliminated (`log_tot_elig` on `Y_debiased`: p = 0.5617). Adding `log(N)` to the primary leaves 0.000286 (p = 0.00660; original 0.0001, p = 0.1802), and `+N+Depth` leaves 0.00024 (p = 0.01893; original 0.0001, p = 0.3399). The `log(N)` coefficient itself is -0.00477 (p = 1.74e-72; original -0.00483, p = 4.76e-74).

The structural caveats are unchanged: decomposing `log(N) = log(N_tot_elig) + log(ct_frac)` still shows a genuine `log(ct_frac)`-specific residual coupling (beta = -0.00111, p = 4.48e-07 on `Y_debiased`; Wald test beta(log_tot_elig) not equal beta(log_ct_frac): t = 4.53, p = 6.84e-06; original t = 4.48, p = 8.84e-06). Naive CD4 T still has zero internal substate resolution (`within-state` N/A), and its confirming `reweighted` arm still conditions on global 30-state PBMC composition (`comp_scope = global_donor_composition`), not internal substates. 12 of 30 states remain imbalanced with SMD above 0.25 after IPW (CD8 Naive 1.27, CD4 TCM 0.73, NK 0.71; original 14 of 30 with CD8 Naive 1.35, CD4 TCM 0.74, NK 0.60), with median SMD 0.21 (original 0.23). IPW still steepens within-pool `log(N)` decline (-0.0897 to -0.1028 per decade, p = 1.96e-06; original -0.0696 to -0.0961, p = 4.80e-06), inflating the reweighted slope from 0.00071 to 0.00080. Controlling IPW for `+N+Depth` still leaves a nominally significant slope (0.000262, p = 0.02759; original 0.0002, p = 0.0823, null). Li and Li (2019) Overlap Weights are nominally significant before adjustment (0.00054, p = 0.01236; original 0.0005, p = 0.0545) but null after `+log(N)` (0.00024, p = 0.1722; original 0.0001, p = 0.3700) and after `+N+Depth` (0.00016, p = 0.3464; original 0.0001, p = 0.6152). Finally, `unadjusted` (0.00003, p = 0.8142), `rarefied 1500` (0.00008, p = 0.7971), `Gene Frac OLD>M` (0.5253, p = 0.1529), `Shared MATURE Trend` (p = 0.7370), and `Centered Height h_hat` (p = 0.6285) are all null, exactly as in the original.

Interpretation: under the frozen rule the naive CD4 T PROVES verdict survives the mature restriction and the signal is larger, but the original's central mechanistic claim (finite-N calibration plus depth eliminates the signal) does NOT reproduce at 30-plus. What survives is a smaller residual slope after calibration and depth adjustment, still significant, sitting on top of a still-present `log(ct_frac)` coupling and an IPW arm that amplifies rather than removes the `log(N)` decline. That is not proof of a biological dispersion rise; it is a weaker calibration story than the original had, and the verdict rests on the frozen rule, not on a resolved mechanism.

### 5.8 Gao Myeloid Gene-Level Secondary Audit (ENCODE Ambient-RNA Soup Artifact, Mature)

The Myeloid gene-fraction signal survives at 30-plus (frac 0.6658, p = 0.002999, q_ct = 0.02699, n_genes = 3842, evaluated MATURE plus OLD donors n = 48) while continuous `primary` (0.00068, p = 0.3338), `within-state` (p = 0.1394), and `Shared MATURE Trend` (p = 0.2744) are null, same pattern as the original. The per-study audit again localizes it to `ENCODE v4 (Snyder)`: ENCODE frac 0.6817 versus Koenig 2022 0.5088 and Chaffin 2022 0.5177 (original 0.7044 versus 0.4787 and 0.5181). Dropping `ENCODE v4 (Snyder)` (n = 23 of 48) collapses frac to 0.5075 and flips the `G_sec` mean difference negative (-0.00471, p = 0.1497; original -0.00446, p = 0.0561). The 20-gene sarcomeric cardiomyocyte UMI fraction inside Myeloid nuclei is again heavy specifically in ENCODE OLD donors (8.00 percent in ENCODE MATURE versus 12.80 percent, max 23.4 percent, in ENCODE OLD, versus 0.53 percent in non-ENCODE studies; original 6.91 percent versus 12.80 percent versus 0.58 percent), correlating r = 0.884 with Myeloid `mean(r)` on `G_sec` (original r = 0.886). Adding the single scalar cardiomyocyte ambient-RNA fraction collapses the `G_sec` OLD-minus-MATURE difference from 0.01702 (p = 0.01053) to 0.00637 (p = 0.07971) and frac to 0.5573 (original: from 0.01531, p = 0.0041, to 0.00372, p = 0.1940, and frac to 0.5372). The ENCODE ambient-RNA explanation of the Myeloid gene-fraction signal therefore reproduces at 30-plus.

## 6. Survived vs Died vs Untestable Against the Original (Per Cell Type)

Verdict comparison table. "Survived" means the same verdict under the frozen rule; "died" means the verdict changed or the unit became untestable; "untestable" means UNTESTABLE in both runs. MDE80 is in SD per decade.

| Dataset | Cell Type | Original Verdict | Mature Verdict | Status | Note |
|---|---|---|---|---|---|
| gao | Adipocyte | DOES NOT PROVE | DOES NOT PROVE | Survived | MDE80 0.255 to 0.307 |
| gao | Cardiomyocyte | DOES NOT PROVE | DOES NOT PROVE | Survived | MDE80 0.157 to 0.196; reweighted now EVALUATED, null |
| gao | Endocardial | DOES NOT PROVE | DOES NOT PROVE | Survived | MDE80 0.208 to 0.269; gene frac p 0.0390 to 0.02799 |
| gao | Endothelial | DOES NOT PROVE | DOES NOT PROVE | Survived | MDE80 0.210 to 0.266; within-state p 0.0230 to 0.02009, q 0.2070 to 0.1620, still a blip |
| gao | Epicardial | DOES NOT PROVE | INCONCLUSIVE | Died (power) | MDE80 0.461 to 0.606, n 19 to 17; no new signal, just lost power |
| gao | Fibroblast | DOES NOT PROVE | DOES NOT PROVE | Survived | MDE80 0.189 to 0.239; reweighted now EVALUATED, null |
| gao | LEC | DOES NOT PROVE | DOES NOT PROVE | Survived | MDE80 0.378 to 0.432; gene fraction died by floor (MATURE n = 4) |
| gao | Lymphoid | DOES NOT PROVE | DOES NOT PROVE | Survived | MDE80 0.239 to 0.273 |
| gao | Mast | INCONCLUSIVE | INCONCLUSIVE | Survived | MDE80 0.590 to 0.795, n 17 to 16 |
| gao | Myeloid | DOES NOT PROVE | DOES NOT PROVE | Survived | MDE80 0.190 to 0.237; gene frac significant in both, ENCODE artifact in both |
| gao | Neuronal | DOES NOT PROVE | DOES NOT PROVE | Survived | MDE80 0.272 to 0.317 |
| gao | Pericyte | DOES NOT PROVE | DOES NOT PROVE | Survived | MDE80 0.176 to 0.219 |
| gao | vSMC | DOES NOT PROVE | DOES NOT PROVE | Survived | MDE80 0.212 to 0.270; within-state p 0.0543 to 0.03600, q 0.2443 to 0.1620, still a blip |
| gao | All_LV_Nuclei (pooled) | DOES NOT PROVE | DOES NOT PROVE | Survived | MDE80 0.130 to 0.167 |
| yazar | CD14-low, CD16-positive monocyte | DOES NOT PROVE | DOES NOT PROVE | Survived | MDE80 0.143 to 0.172 |
| yazar | CD14-positive monocyte | DOES NOT PROVE | DOES NOT PROVE | Survived | MDE80 0.089 to 0.110 |
| yazar | CD4-positive, alpha-beta cytotoxic T cell | DOES NOT PROVE | DOES NOT PROVE | Survived | MDE80 0.311 to 0.390 |
| yazar | central memory CD4-positive, alpha-beta T cell | DOES NOT PROVE | DOES NOT PROVE | Survived | MDE80 0.059 to 0.071 |
| yazar | central memory CD8-positive, alpha-beta T cell | DOES NOT PROVE | DOES NOT PROVE | Survived | MDE80 0.290 to 0.352 |
| yazar | effector memory CD4-positive, alpha-beta T cell | DOES NOT PROVE | DOES NOT PROVE | Survived | MDE80 0.146 to 0.182 |
| yazar | effector memory CD8-positive, alpha-beta T cell | DOES NOT PROVE | DOES NOT PROVE | Survived | MDE80 0.074 to 0.089 |
| yazar | gamma-delta T cell | DOES NOT PROVE | DOES NOT PROVE | Survived | MDE80 0.189 to 0.251; nominally significant negative primary in both runs |
| yazar | memory B cell | DOES NOT PROVE | DOES NOT PROVE | Survived | MDE80 0.142 to 0.184 |
| yazar | mucosal-associated invariant T cell | DOES NOT PROVE | INCONCLUSIVE | Died (power) | MDE80 0.473 to 0.657, n 30 to 25; gene fraction also died by floor (MATURE n = 3) |
| yazar | naive B cell | DOES NOT PROVE | DOES NOT PROVE | Survived | MDE80 0.127 to 0.153 |
| yazar | naive thymus-derived CD4-positive, alpha-beta T cell | PROVES | PROVES | Survived (strengthened) | Primary p 3.00e-05 to 5.174e-08, q 4.80e-04 to 8.279e-07; slope 0.00046 to 0.00071; calibration kill does NOT reproduce (Section 5.7) |
| yazar | naive thymus-derived CD8-positive, alpha-beta T cell | DOES NOT PROVE | DOES NOT PROVE | Survived | MDE80 0.102 to 0.134; reweighted now EVALUATED (0.00055, p = 0.05725), null |
| yazar | natural killer cell | DOES NOT PROVE | DOES NOT PROVE | Survived | MDE80 0.068 to 0.082 |
| yazar | regulatory T cell | INCONCLUSIVE | INCONCLUSIVE | Survived | MDE80 0.505 in both runs, n = 54 in both runs, p = 0.3643 in both runs (no Treg donors under 30) |
| yazar | transitional stage B cell | DOES NOT PROVE | DOES NOT PROVE | Survived | MDE80 0.155 to 0.219 |
| yazar | All_Cells (pooled) | INCONCLUSIVE | DOES NOT PROVE | Died (flipped) | Reweighted slope flips sign (-0.00001 to 0.00015), removing the 2-arm disagreement that forced INCONCLUSIVE; primary stays a powered null |
| yazar | 12 untestable types (CD16-negative CD56-bright NK, CD4 T, CD8 T, conventional DC, dendritic cell, double negative thymocyte, erythrocyte, hematopoietic precursor cell, innate lymphoid cell, plasmablast, plasmacytoid DC, platelet) | UNTESTABLE | UNTESTABLE | Untestable in both | Same 12 types; plasmablast status code unchanged (NOT_EVALUABLE_WITHIN_BATCH_DF) |

Summary of the comparison: 27 of 29 non-pooled cell types with an original verdict keep the same verdict (12 of 13 Gao, 15 of 16 Yazar testable). Two die by power loss only (Gao Epicardial, Yazar MAIT), with no new positive signal in either. One pooled row flips (Yazar All_Cells) through the reweighted sign change, not through a new signal. One verdict strengthens (Yazar naive CD4 T PROVES). No cell type flips from DOES NOT PROVE to PROVES and none flips from PROVES to DOES NOT PROVE. Every MDE80 rises or stays flat at 30-plus (fewer donors), by 0.00 to 0.21 SD per decade; the only type with identical numbers in both runs is Yazar regulatory T cell, because it has no donors under 30.

## 7. Scientific Interpretation

- **Gao Heart (DOES NOT PROVE, unchanged):** Zero of 13 non-pooled cell types pass BH (minimum primary q = 0.5263; original 0.8885), 11 of 13 are powered nulls (4 high-power at most 0.25 SD per decade and 7 moderate-power 0.25 to 0.50 SD per decade, plus pooled `All_LV_Nuclei` at 0.167 SD per decade), and only sparse Epicardial (n = 17, MDE80 = 0.606) and Mast (n = 16, MDE80 = 0.795) are underpowered INCONCLUSIVE. `within-state` (m = 9) and `reweighted` (m = 4) produce zero BH hits, Li and Li (2019) Overlap Weights are null across all 10 evaluable rows (p = 0.0897 to 0.9590), and `Shared MATURE Trend + N+Depth` is null everywhere (p = 0.0756 to 0.9736). The Myeloid gene-fraction signal is again fully explained as ENCODE cardiomyocyte ambient-RNA soup (Section 5.8). This agrees with the original Gao verdict and with A2 between-subject variance (0 BH hits across all Gao cell types in the original).
- **Yazar Blood (INCONCLUSIVE, unchanged, 1 of 16 testable cell types PROVES, 13 DOES NOT PROVE, 2 INCONCLUSIVE, 12 of 28 UNTESTABLE):** 13 of 16 testable PBMC types are powered nulls (MDE80 = 0.064 to 0.390 SD per decade; 10 of 13 high-power at most 0.25 SD per decade). The sole cell type satisfying the pre-registered PROVES rule is again naive thymus-derived CD4-positive, alpha-beta T cell, with a larger slope and a smaller p-value than the original. The important change is mechanistic, not verdict-level: the original's finite-N calibration kill does not reproduce at 30-plus (Section 5.7). The calibration curve still explains about half the slope, but the debiased slope, the `+log(N)` slope, the `+N+Depth` slope, and the IPW `+N+Depth` slope all stay nominally significant at 30-plus, where the original killed each of them. Against that, the same structural caveats still apply (steep `log(N)` decline, strong negative `Y` to `log(N)` correlation, persistent `log(ct_frac)` coupling, global-composition IPW that steepens the `log(N)` decline, 12 of 30 states imbalanced after IPW, overlap weights null after N and depth adjustment, and null unadjusted, rarefied, gene-fraction, and shared-trend arms). The honest statement at 30-plus is: the naive CD4 T signal survives the mature restriction and survives finite-N calibration better than it did in the full-range run, but its mechanism is still unresolved, and the frozen PROVES verdict rests on the rule, not on a resolved mechanism.
- **Pooled All_Cells flip:** the pooled Yazar row moves INCONCLUSIVE to DOES NOT PROVE only because the reweighted slope changes sign from -0.00001 (p = 0.9286) to 0.00015 (p = 0.3336), so the two composition arms no longer disagree in sign and Branch 4b no longer fires. The within-state arm is still nominally significant (0.00021, p = 0.0006564) and is still driven by CD4 Naive (leave-one-out collapses it to p = 0.3219). Nothing about the pooled biology changed; the verdict rule is sensitive to the sign of a null reweighted slope.
- **Reconciliation with Long et al. (2025):** the coarse-lineage picture mostly reproduces, with two differences. Monocyte pooling no longer produces even an unmatched rise (p = 0.0780; original p = 0.0196). CD4_T within-L2 plus N-plus-Depth no longer collapses (p = 0.01986; original p = 0.1604), which is the same CD4 Naive finite-N story as Section 5.7, now at lineage level. B, CD8_T, and NK remain null on unmatched LAD except for a nominally significant negative CD8_T unmatched slope (-0.00027, p = 0.02493). Overall, full-depth per-donor residual dispersion at 30-plus still provides no replicated evidence of an intrinsic age-related transcriptional noise rise in either human heart or blood once composition, ambient-RNA, and cell-count effects are accounted for, with the single frozen-rule exception of Yazar naive CD4 T, whose mechanism is unresolved.

## 8. Reproduction

```bash
PY=<agescape-conda-env>/bin/python
$PY scripts/run_mature_block1_residual_dispersion.py --stage all --datasets gao,yazar
```

Cache is reused read-only from the frozen `block1_cache` through symlinks in `block1_cache_mature`; the script never rebuilds the cache (it exits the cache stage on the existing `CACHE_DONE.json`) and never writes to the frozen cache. Conda env `agescape` is activated from `/home/reza/miniconda3/etc/profile.d/conda.sh` (conda is not on PATH in non-interactive shells on beast).

Files:

- `scripts/run_mature_block1_residual_dispersion.py` (remapped script; frozen original untouched)
- `results/followup-mature-block1-residual-dispersion.csv` (215 rows: 14 Gao groups times 5 arms plus 29 Yazar groups times 5 arms)
- `results/followup-mature-block1-gene-fractions.csv` (43 rows: 14 Gao plus 29 Yazar)
- `results/followup-mature-block1-bh-check.csv` (138 rows: independent BH recomputation, 0 mismatches, max absolute difference 9.99e-16)
- `reports/followup-mature-block1.md` (this report)
- Diagnostics JSON: `block1_diagnostics_summary.json` under the mature cache root (`tmp_organs/block1_cache_mature/`), with retention, nuisance, crosschecks, calibration, coarse-lineage, and mechanism-audit values quoted in Section 5
- Mature stats: `tmp_organs/block1_cache_mature/gao/stats/` and `tmp_organs/block1_cache_mature/yazar/stats/` (donor-group stats, donor-state stats, finite-N calibration, per-gene residual pieces)

Nothing was submitted anywhere. The original report, script, results, and cache are untouched by this rerun.
