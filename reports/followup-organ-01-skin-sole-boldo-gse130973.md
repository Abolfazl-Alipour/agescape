# Dataset 01: Human Skin Single-Cell Aging (Sole-Boldo et al., GSE130973)

## 1. Study Verification and Data Provenance
- **Accession**: `GSE130973` (`GSM3758115` to `GSM3758119`), Primary source: Sole-Boldo et al., *Communications Biology* 2020 (`PMC7181753`, `DOI: 10.1038/s42003-020-0922-4`).
- **Tissue and Sampling**: Sun-protected inguinoiliac whole human skin (epidermis + dermis, unbiased single-cell droplet capture, 10x Chromium v2), excluding photoaged sites.
- **Verified Donor Cohort (`n = 5` males)**:
  - `S1` (`y1`, `GSM3758115`): age `25` yr, Male, Caucasian, `2,784` QC cells.
  - `S2` (`y2`, `GSM3758116`): age `27` yr, Male, Caucasian, `2,670` QC cells.
  - `S3` (`o1`, `GSM3758117`): age `53` yr, Male, Caucasian, `3,324` QC cells.
  - `S4` (`o2`, `GSM3758118`): age `70` yr, Male, Caucasian, `2,144` QC cells.
  - `S5` (`o3`, `GSM3758119`): age `69` yr, Male, Caucasian, `4,535` QC cells.
  - **Total verified cells**: `15,457` cells across `21,352` Seurat raw matrix genes (`199` technical `MT-/RP[LS]/MRP[LS]/sex-linked/ERCC` genes stripped prior to library-size normalization). Note all 5 donors are white Caucasian males (0 female donors).

## 2. Explicit Multiple-Testing Families and Methodological Design
- **Inferential unit**: Donor (`n = 5`). Age modeled continuously as `age_per_decade` (`age_years / 10`).
- **Family A (CLR cell-type fractions)**: OLS across the `9` major skin cell types (clusters mapped 1-to-1 to the 9 published identities via Supplementary Data 1 markers; pseudocount `0.5`). Unadjusted continuous age (`residual_df = 3`), BH correction across `K_A = 9`.
- **Family B (Donor-level gene age test + cross-cell-type concordance)**: Donor pseudobulk (`sum raw UMIs` -> `log2(CPM + 1)`, genes expressed at `CPM >= 1.0` in `>= 3` donors per cell type) across eligible major cell types (`>= 20` cells in `>= 4` donors). Primary model: `age_per_decade + log_lib` (`residual_df = 2` for `n=5`), sensitivity model: unadjusted `age_per_decade` (`residual_df = 3`). Global BH across all `M_B = 89,068` `(cell_type, gene)` tests. Cross-cell-type concordance tested on common expressed genes via mean pairwise Pearson `r` across cell-type `t`-statistics with `1,000` donor-age permutations.
- **Family C (`P1_Clean` 139-gene mouse aging signature)**: `76` mapped human UP and `42` DOWN symbols (`load_program_tiers`). Composite (`up_mean - down_mean`) OLS slope tested analytically and against `200` expression-decile-matched random gene sets (`matched_null_effect_pvalue`), with BH across `K_C = 7` cell types.

## 3. Results

### Family A: Cell-Type Composition (CLR Fractions, 9 Major Cell Types)
Minimum BH-corrected `q_BH = 0.0123` (`1` significant compositional shifts at `q < 0.05` in this `n = 5` sample; top nominal cell type `Pericytes` (slope `-0.312` CLR/decade, `p = 0.0014`, `q = 0.0123`)).

| Cell Type | Total Cells | Median Frac | CLR Slope (/10 yr) | p-value | q_BH (K=9) |
|---|---:|---:|---:|---:|---:|
| Pericytes | 1220 | 0.0665 | -0.312 | 0.0014 | 0.0123 |
| T_cells | 1281 | 0.0278 | -0.528 | 0.0392 | 0.1764 |
| Keratinocytes | 2586 | 0.1336 | +0.406 | 0.0738 | 0.2215 |
| Lymphatic_EC | 294 | 0.0176 | +0.180 | 0.1445 | 0.3252 |
| Vascular_EC | 1467 | 0.0924 | +0.083 | 0.2435 | 0.4382 |
| Melanocytes | 123 | 0.0067 | +0.083 | 0.3103 | 0.4655 |
| Fibroblasts | 5948 | 0.3330 | +0.080 | 0.4441 | 0.5709 |
| Macrophages_DC | 2228 | 0.1453 | +0.084 | 0.6311 | 0.7100 |
| Erythrocytes | 310 | 0.0042 | -0.077 | 0.8791 | 0.8791 |

### Family B: Donor-Level Gene Age Test and Cross-Cell-Type Concordance
Across `7` major cell types (`M_B = 89,068` tests, primary `age_per_decade + log_lib`), `n_sig_global_q05 = 0` genes passed `q_BH < 0.05` (minimum `q_BH = 0.4598`).
Cross-cell-type directional concordance across `9199` common expressed genes across `7` major cell types:
- Primary (`age_per_decade + log_lib`): mean pairwise Pearson `r = +0.0947` (`perm_p_one_sided = 0.0470`, `n_pairs = 21`).
- Sensitivity (`age_per_decade` unadjusted): mean pairwise Pearson `r = +0.2369` (`perm_p_one_sided = 0.1039`).

| Cell Type | Donors | Genes Tested | Global q<0.05 | Min p-value | Min Global q_BH |
|---|---:|---:|---:|---:|---:|
| Fibroblasts | 5 | 13766 | 0 | 1.26e-05 | 0.4598 |
| Keratinocytes | 5 | 13407 | 0 | 2.31e-04 | 0.4598 |
| Lymphatic_EC | 5 | 10973 | 0 | 2.57e-05 | 0.4598 |
| Macrophages_DC | 5 | 12870 | 0 | 5.92e-05 | 0.4598 |
| Pericytes | 5 | 13494 | 0 | 5.38e-04 | 0.4598 |
| T_cells | 5 | 11280 | 0 | 8.34e-06 | 0.4598 |
| Vascular_EC | 5 | 13278 | 0 | 9.78e-05 | 0.4598 |

### Family C: Mouse `P1_Clean` Signature Check
In `K_C = 7` major skin cell types (primary `age_per_decade + log_lib`), **no cell type showed statistically significant positive `P1_Clean` transport** (minimum `q_bh_composite = 0.7017`; `Keratinocytes` reached two-sided empirical `q_bh_empirical = 0.0348` in the **opposite/negative** direction `composite = -0.0940`, `95% CI [-0.330, +0.142]`, non-significant `analytic p = 0.2285`; `composite = -0.0600`, `analytic p = 0.2639`, `q_emp = 0.1625` unadjusted). Results in unadjusted sensitivity (`age_per_decade` only) and across the 17 fine clusters are saved in `results/followup-organ-01-skin-sole-boldo-gse130973-p1-signature.csv`.

| Cell Type | Donors | UP/DOWN Detected | Composite Slope | Analytic p | Empirical Null p (200) | Empirical q_BH |
|---|---:|---|---:|---:|---:|---:|
| Keratinocytes | 5 | 62/31 | -0.0940 | 0.2285 | 0.0050 | 0.0348 |
| Pericytes | 5 | 65/31 | +0.1477 | 0.6304 | 0.0547 | 0.1306 |
| Vascular_EC | 5 | 63/30 | -0.0561 | 0.5012 | 0.0647 | 0.1306 |
| T_cells | 5 | 64/29 | -0.2799 | 0.3180 | 0.0746 | 0.1306 |
| Lymphatic_EC | 5 | 57/25 | +0.0825 | 0.4737 | 0.1095 | 0.1509 |
| Macrophages_DC | 5 | 65/31 | -0.0606 | 0.3731 | 0.1294 | 0.1509 |
| Fibroblasts | 5 | 65/35 | +0.0147 | 0.8395 | 0.7114 | 0.7114 |

## 4. Limitations and Verdict
- **Sample-size and demographic limits**: Only `n = 5` male Caucasian donors (`25, 27, 53, 69, 70` yr), with 0 female donors and a gap between ages 27 and 53. Power is very limited (`df = 2` adjusted, `df = 3` unadjusted).
- **Consistent with Agescape core claims**: Zero global FDR-significant gene hits (`M_B = 89,068`, `min q_BH = 0.4598`) and no positive `P1_Clean` aging transport (`min q_bh_composite = 0.7017`; the two-sided empirical outlier in `Keratinocytes` is negative `-0.0940` with `analytic p = 0.2285` and vanishes unadjusted `q_emp = 0.1625`). One CLR compositional shift passed BH in Family A (`Pericytes` `-0.312` CLR/decade, `q_BH = 0.0123`, `n = 5` males).

## 5. Output Artifacts
- [Donors & Cell Counts](../results/followup-organ-01-skin-sole-boldo-gse130973-donors.csv)
- [Main Fractions](../results/followup-organ-01-skin-sole-boldo-gse130973-fractions-main.csv)
- [Fine Fractions](../results/followup-organ-01-skin-sole-boldo-gse130973-fractions-fine.csv)
- [Gene Summary](../results/followup-organ-01-skin-sole-boldo-gse130973-gene-summary.csv)
- [Top Genes](../results/followup-organ-01-skin-sole-boldo-gse130973-top-genes.csv)
- [Concordance](../results/followup-organ-01-skin-sole-boldo-gse130973-concordance.csv)
- [P1 Signature](../results/followup-organ-01-skin-sole-boldo-gse130973-p1-signature.csv)
