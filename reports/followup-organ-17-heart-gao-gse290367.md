# Follow-Up Organ 17 (Item #14): Heart (Gao et al. 2026 Integrated Atlas, GSE290367)

## Summary Verdict

In this unbiased snRNA-seq atlas of **117 deduplicated postnatal non-diseased human left ventricle (LV) donors** (**901,946 nuclei**, ages **20.0 to 91.0 years**, **69 male / 48 female**) across 7 exact-age component studies in `GSE290367`:

- **Family A (Cell-type composition, 13 CLR tests):** **0 significant shifts** after BH correction (`q_bh_family_a >= 0.3128`). The strongest nominal trend is `Endothelial` (`-0.1511` CLR/decade, `SE = 0.0712`, `p = 0.0362`, `q_bh = 0.3128`).
- **Family B (Donor-level gene expression, 145,873 tests across 13 cell types; 159,052 including `All_LV_Nuclei`):** **8 cell-type gene hits** at `q_bh_celltypes13 < 0.05` (`Fibroblast` 3: `LRRIQ1`, `EPHA6`, `PTCHD4`; `Cardiomyocyte` 3: `MPL`, `ADD1`, `FCHSD2`; `Neuronal` 1: `NDUFAF1`; `LEC` 1: `HELQ`) and **12 combined hits** at `q_bh_all14 < 0.05` (adding `Cardiomyocyte` `STRBP` plus `All_LV_Nuclei` `AACS`, `SLC10A6`, `TNFRSF19`).
- **Cross-cell-type concordance (`78` pairs across `13` cell types, `7,502` shared genes, `1,000` donor permutations):** **Not significant** (`mean Pearson r = 0.1309`, `median r = 0.1099`, `mean sign agreement = 0.5466`, `perm_p = 0.1998`). Consistent with the published paper's null human cross-cell-type concordance (`p = 0.70`).
- **Family C (`P1_Clean` 139-gene mouse aging signature, `14` groups, `2,000` matched-gene null iterations):** **Notable divergence from the paper's surveyed human brain/kidney/lung null** (`lowest q = 0.052` in 252 tests). In human LV, `Endocardial` (`N = 87`, `composite = +0.0879/decade`, `composite_p = 0.00107`, `q_bh_composite = 0.0150`, `empirical_p = 0.0020`, `q_bh_empirical = 0.0093`) is significant on both parametric and empirical tests, and `All_LV_Nuclei` (`+0.0380/decade`, `comp_p = 0.0110`, `emp_p = 0.00050`, `q_emp = 0.0070`), `Fibroblast` (`+0.0348`, `comp_p = 0.0751`, `emp_p = 0.0010`, `q_emp = 0.0070`), `Myeloid` (`+0.0379`, `emp_p = 0.0035`, `q_emp = 0.0122`), `Pericyte` (`+0.0400`, `emp_p = 0.0060`, `q_emp = 0.0168`), and `Endothelial` (`+0.0303`, `emp_p = 0.0125`, `q_emp = 0.0260`) exceed matched random-gene nulls positively (`Adipocyte` `N = 46` shifts oppositely at `-0.0941/decade`, `comp_p = 0.0994`, `q_emp = 0.0260`). Positive enrichment persists in 10x-only (`N = 108`), studies-with-`N >= 8` (`N = 110`), and standalone single-lab cohorts (`ENCODE v4` `N = 54`, `Koenig 2022` `N = 23`). This is **not a contradiction** of any existing heart claim in the paper (which did not analyze heart), but it is flagged for manuscript scope calibration because human LV carries modest `P1_Clean` signal driven primarily by `P1_DOWN` decline (`-0.0337/decade` in `All_LV_Nuclei`) and `Endocardial` `P1_UP` (`+0.0502/decade`).

---

## Primary-Source Verification and Donor Disposition

- **Accession & Primary Source:** GEO `GSE290367` (`GSE290367_integrated_RNA_adata.h5ad`, STARSolo v2.7.11a integer counts; `wgao688/Human-cardiac-multiome-analysis`), verified across all `901,946` primary nuclei (`1,259,162,328` nonzeros; `nonfinite = 0`, `negative = 0`, `noninteger = 0`; [`followup-organ-17-heart-gao-metadata-summary.json`](../results/followup-organ-17-heart-gao-metadata-summary.json)).
- **Complete 299-Donor Audit (`299 -> 117` primary donors; [`followup-organ-17-heart-gao-donors.csv`](../results/followup-organ-17-heart-gao-donors.csv)):**
  - Excluded **85** midpoint-bracket donors (`Reichart 2022` 67, `Litvinukova 2020` 14, `Kanemaru 2023` 4);
  - Excluded **71** disease donors (`DCM` 25, `HCM` 15, `ICM` 15, `pediatric HF` 9, `AMI` 7);
  - Excluded **13** fetal donors (`ENCODE v4` 8, `Sim 2021` 3, `Penn` 2);
  - Excluded **7** pediatric donors `< 18 years` (`Sim 2021` 3, `Hill 2022` 2, `Koenig 2022` 2);
  - Excluded **6** duplicate `Penn` Heart Bank Drop-seq donor copies whose matched donors appear as higher-depth 10x profiles in `Chaffin 2022` (`K1549/P1549`, `K1561/P1561`, `K1622/P1622`, `K1702/P1702`, `K1718/P1718`) and `Simonson 2023` (`K1739/P1739`). Dropping `Penn:K1561` also removes the GEO SOFT typo (`73` vs verified `65.0` years in both Penn and Chaffin metadata).
- **Verified Unbiased Sampling & Non-Overlap:** All 117 included donors have `region == LV`, `disease == ND`, and `age_status == postnatal`. Primary-source code audit confirmed `Penn` pan-nuclear DAPI+ FACS/nonFACS gating includes all 13 cardiac cell types (though Drop-seq captures fewer cardiomyocytes, `7 to 219` nuclei/donor vs 10x `128 to 8,267`), and `ENCODE v4 (Snyder)` biosamples (`54` unique `ENCSR*` accessions) have zero overlap with the other 6 studies or with Kanemaru/Litvinukova/Jia cohorts.
- **Included Component Studies (`117` donors, `901,946` nuclei):**

| Study | Assay | Donors | Nuclei | Age Range (yr) | Sex (M/F) | Source File |
|---|---|---:|---:|---|---|---|
| ENCODE v4 (Snyder) | Multiome-v1 | 54 | 446,726 | 20.0 to 69.0 | 35M / 19F | [`followup-organ-17-heart-gao-cellcounts.csv`](../results/followup-organ-17-heart-gao-cellcounts.csv) |
| Koenig 2022 | 5prime-v1 | 23 | 119,412 | 20.0 to 75.0 | 11M / 12F | [`followup-organ-17-heart-gao-cellcounts.csv`](../results/followup-organ-17-heart-gao-cellcounts.csv) |
| Chaffin 2022 | 3prime-v3 | 16 | 179,922 | 33.0 to 72.0 | 6M / 10F | [`followup-organ-17-heart-gao-cellcounts.csv`](../results/followup-organ-17-heart-gao-cellcounts.csv) |
| Penn (unique) | Dropseq | 9 | 50,160 | 22.0 to 91.0 | 8M / 1F | [`followup-organ-17-heart-gao-cellcounts.csv`](../results/followup-organ-17-heart-gao-cellcounts.csv) |
| Simonson 2023 | 3prime-v3 | 8 | 58,727 | 42.0 to 64.0 | 4M / 4F | [`followup-organ-17-heart-gao-cellcounts.csv`](../results/followup-organ-17-heart-gao-cellcounts.csv) |
| Kuppe 2022 | 3prime-v3 | 4 | 38,727 | 44.0 to 61.0 | 3M / 1F | [`followup-organ-17-heart-gao-cellcounts.csv`](../results/followup-organ-17-heart-gao-cellcounts.csv) |
| Sim 2021 | 3prime-v1 | 3 | 8,272 | 35.0 to 42.0 | 2M / 1F | [`followup-organ-17-heart-gao-cellcounts.csv`](../results/followup-organ-17-heart-gao-cellcounts.csv) |
| **Total Primary** | **7 cohorts** | **117** | **901,946** | **20.0 to 91.0** | **69M / 48F** | [`followup-organ-17-heart-gao-donors.csv`](../results/followup-organ-17-heart-gao-donors.csv) |

---

## Multiple-Testing Families and Statistical Design

Donor is the exclusive unit of replication across all families (`scripts/run_followup_organ_17_heart_gao.py`). All primary models regress on `age_per_decade` adjusting for `sex_male` and 6 study indicators (`reference = ENCODE v4 (Snyder)`), plus donor-by-group `log_lib` for expression and signature families.

1. **Family A (CLR Composition):** 13 cardiac cell types (`residual_df = 108`), BH corrected within Family A (`family_a_size = 13`).
2. **Family B (Donor-Level Gene Age Tests):** `145,873` OLS tests across the 13 cardiac cell types (`min_donors = 4`, `CPM >= 1.0` in `>= 50%` of donors within group; `15,930` non-technical genes evaluated from `16,115` raw genes), BH corrected across `145,873` (`q_bh_celltypes13`) and across all `159,052` tests including `All_LV_Nuclei` (`13,179` genes, `q_bh_all14`), plus within-group BH (`q_bh_within_group`). Cross-cell-type concordance (`13` cell types, `78` pairs, `7,502` common genes) uses `1,000` donor-age permutations.
3. **Family C (`P1_Clean` Signature):** 14 groups (13 cell types + `All_LV_Nuclei`), `family_c_size = 14`, using `2,000` expression-matched random gene-set nulls per group (`q_bh_composite` and `q_bh_empirical`).

---

## Results

### Family A: Cell-Type Fractions (`followup-organ-17-heart-gao-fractions.csv`)

No cell type reaches `q_bh_family_a < 0.05` (`residual_df = 108` across `N = 117` donors):

| Cell Type | Donors Present | Total Nuclei | Median Fraction | CLR Effect / Decade | SE | p-value | q_bh (13) | Source File |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| Endothelial | 117 | 91,384 | 0.1010 | -0.1511 | 0.0712 | 0.0362 | 0.3128 | [`followup-organ-17-heart-gao-fractions.csv`](../results/followup-organ-17-heart-gao-fractions.csv) |
| Epicardial | 57 | 2,964 | 0.0000 | +0.2026 | 0.1066 | 0.0601 | 0.3128 | [`followup-organ-17-heart-gao-fractions.csv`](../results/followup-organ-17-heart-gao-fractions.csv) |
| Lymphoid | 117 | 32,635 | 0.0201 | +0.0853 | 0.0489 | 0.0838 | 0.3128 | [`followup-organ-17-heart-gao-fractions.csv`](../results/followup-organ-17-heart-gao-fractions.csv) |
| Neuronal | 117 | 8,824 | 0.0098 | -0.0689 | 0.0428 | 0.1102 | 0.3128 | [`followup-organ-17-heart-gao-fractions.csv`](../results/followup-organ-17-heart-gao-fractions.csv) |
| Mast | 111 | 5,158 | 0.0038 | -0.1231 | 0.0786 | 0.1203 | 0.3128 | [`followup-organ-17-heart-gao-fractions.csv`](../results/followup-organ-17-heart-gao-fractions.csv) |
| vSMC | 117 | 19,223 | 0.0205 | +0.0361 | 0.0303 | 0.2352 | 0.5096 | [`followup-organ-17-heart-gao-fractions.csv`](../results/followup-organ-17-heart-gao-fractions.csv) |
| Pericyte | 117 | 108,979 | 0.1148 | -0.0334 | 0.0308 | 0.2807 | 0.5144 | [`followup-organ-17-heart-gao-fractions.csv`](../results/followup-organ-17-heart-gao-fractions.csv) |
| LEC | 114 | 4,899 | 0.0033 | +0.0692 | 0.0688 | 0.3166 | 0.5144 | [`followup-organ-17-heart-gao-fractions.csv`](../results/followup-organ-17-heart-gao-fractions.csv) |
| Cardiomyocyte | 117 | 267,584 | 0.2659 | -0.0221 | 0.0496 | 0.6559 | 0.8776 | [`followup-organ-17-heart-gao-fractions.csv`](../results/followup-organ-17-heart-gao-fractions.csv) |
| Endocardial | 115 | 20,471 | 0.0105 | -0.0270 | 0.0978 | 0.7828 | 0.8776 | [`followup-organ-17-heart-gao-fractions.csv`](../results/followup-organ-17-heart-gao-fractions.csv) |
| Adipocyte | 97 | 5,675 | 0.0014 | +0.0277 | 0.1064 | 0.7953 | 0.8776 | [`followup-organ-17-heart-gao-fractions.csv`](../results/followup-organ-17-heart-gao-fractions.csv) |
| Myeloid | 117 | 119,684 | 0.1358 | +0.0069 | 0.0285 | 0.8101 | 0.8776 | [`followup-organ-17-heart-gao-fractions.csv`](../results/followup-organ-17-heart-gao-fractions.csv) |
| Fibroblast | 117 | 214,466 | 0.2392 | -0.0021 | 0.0258 | 0.9358 | 0.9358 | [`followup-organ-17-heart-gao-fractions.csv`](../results/followup-organ-17-heart-gao-fractions.csv) |

### Family B: Gene-Level Aging and Concordance

- Summary counts and QC library depth metrics are in [`followup-organ-17-heart-gao-genes-summary.csv`](../results/followup-organ-17-heart-gao-genes-summary.csv) and [`followup-organ-17-heart-gao-qc.csv`](../results/followup-organ-17-heart-gao-qc.csv).
- Global significant gene hits (`q_bh_all14 < 0.05`, [`followup-organ-17-heart-gao-top-genes.csv`](../results/followup-organ-17-heart-gao-top-genes.csv)):

| Group | Gene | Donors | Effect / Decade | SE | p-value | q_bh_ct13 | q_bh_all14 | q_bh_within | Source File |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| Fibroblast | LRRIQ1 | 117 | +0.2479 | 0.0440 | 1.44e-07 | 0.0148 | 0.0150 | 0.0019 | [`followup-organ-17-heart-gao-top-genes.csv`](../results/followup-organ-17-heart-gao-top-genes.csv) |
| Neuronal | NDUFAF1 | 102 | +0.6726 | 0.1197 | 2.03e-07 | 0.0148 | 0.0150 | 0.0021 | [`followup-organ-17-heart-gao-top-genes.csv`](../results/followup-organ-17-heart-gao-top-genes.csv) |
| All_LV_Nuclei | AACS | 117 | -0.1160 | 0.0212 | 2.83e-07 | NA | 0.0150 | 0.0037 | [`followup-organ-17-heart-gao-top-genes.csv`](../results/followup-organ-17-heart-gao-top-genes.csv) |
| LEC | HELQ | 62 | +0.5194 | 0.0915 | 5.57e-07 | 0.0204 | 0.0178 | 0.0054 | [`followup-organ-17-heart-gao-top-genes.csv`](../results/followup-organ-17-heart-gao-top-genes.csv) |
| Cardiomyocyte | MPL | 115 | +0.2051 | 0.0385 | 5.58e-07 | 0.0204 | 0.0178 | 0.0057 | [`followup-organ-17-heart-gao-top-genes.csv`](../results/followup-organ-17-heart-gao-top-genes.csv) |
| All_LV_Nuclei | SLC10A6 | 117 | +0.1534 | 0.0295 | 9.62e-07 | NA | 0.0216 | 0.0063 | [`followup-organ-17-heart-gao-top-genes.csv`](../results/followup-organ-17-heart-gao-top-genes.csv) |
| Cardiomyocyte | ADD1 | 115 | -0.1035 | 0.0199 | 9.65e-07 | 0.0264 | 0.0216 | 0.0057 | [`followup-organ-17-heart-gao-top-genes.csv`](../results/followup-organ-17-heart-gao-top-genes.csv) |
| Fibroblast | EPHA6 | 117 | +0.3086 | 0.0597 | 1.09e-06 | 0.0264 | 0.0216 | 0.0070 | [`followup-organ-17-heart-gao-top-genes.csv`](../results/followup-organ-17-heart-gao-top-genes.csv) |
| Fibroblast | PTCHD4 | 117 | +0.3394 | 0.0681 | 2.41e-06 | 0.0481 | 0.0419 | 0.0104 | [`followup-organ-17-heart-gao-top-genes.csv`](../results/followup-organ-17-heart-gao-top-genes.csv) |
| Cardiomyocyte | FCHSD2 | 115 | -0.1145 | 0.0231 | 2.64e-06 | 0.0481 | 0.0419 | 0.0095 | [`followup-organ-17-heart-gao-top-genes.csv`](../results/followup-organ-17-heart-gao-top-genes.csv) |
| Cardiomyocyte | STRBP | 115 | -0.0915 | 0.0186 | 3.24e-06 | 0.0525 | 0.0468 | 0.0095 | [`followup-organ-17-heart-gao-top-genes.csv`](../results/followup-organ-17-heart-gao-top-genes.csv) |
| All_LV_Nuclei | TNFRSF19 | 117 | -0.1859 | 0.0381 | 3.69e-06 | NA | 0.0489 | 0.0162 | [`followup-organ-17-heart-gao-top-genes.csv`](../results/followup-organ-17-heart-gao-top-genes.csv) |

- **Concordance ([`followup-organ-17-heart-gao-concordance.csv`](../results/followup-organ-17-heart-gao-concordance.csv)):** Across `13` cardiac cell types (`78` pairs, `7,502` shared genes), `mean_pairwise_pearson_r = 0.1309` (`median = 0.1099`, `min = 0.0191`, `max = 0.4475`, `Spearman rho = 0.1264`, `mean sign agreement = 0.5466`, `perm_p_two_sided = 0.1998`).

### Family C: `P1_Clean` Mouse Signature (`followup-organ-17-heart-gao-p1-signature.csv`)

| Group | Donors | P1 Up / Down Detected | Up Effect (p) | Down Effect (p) | Composite Effect / Decade (p) | Empirical Null p (2,000 iters) | q_bh_composite | q_bh_empirical | Source File |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| All_LV_Nuclei | 117 | 66 / 33 | +0.0042 (0.8393) | -0.0337 (0.1003) | +0.0380 (0.0110) | 0.00050 | 0.0772 | 0.0070 | [`followup-organ-17-heart-gao-p1-signature.csv`](../results/followup-organ-17-heart-gao-p1-signature.csv) |
| Fibroblast | 117 | 65 / 32 | +0.0032 (0.8868) | -0.0316 (0.1202) | +0.0348 (0.0751) | 0.00100 | 0.2070 | 0.0070 | [`followup-organ-17-heart-gao-p1-signature.csv`](../results/followup-organ-17-heart-gao-p1-signature.csv) |
| Endocardial | 87 | 62 / 27 | +0.0502 (0.1515) | -0.0377 (0.2089) | +0.0879 (0.00107) | 0.00200 | 0.0150 | 0.0093 | [`followup-organ-17-heart-gao-p1-signature.csv`](../results/followup-organ-17-heart-gao-p1-signature.csv) |
| Myeloid | 117 | 65 / 31 | -0.0110 (0.5796) | -0.0489 (0.0467) | +0.0379 (0.0441) | 0.00350 | 0.1545 | 0.0122 | [`followup-organ-17-heart-gao-p1-signature.csv`](../results/followup-organ-17-heart-gao-p1-signature.csv) |
| Pericyte | 117 | 66 / 31 | +0.0030 (0.8951) | -0.0370 (0.0985) | +0.0400 (0.0297) | 0.00600 | 0.1384 | 0.0168 | [`followup-organ-17-heart-gao-p1-signature.csv`](../results/followup-organ-17-heart-gao-p1-signature.csv) |
| Endothelial | 111 | 67 / 31 | -0.0052 (0.8268) | -0.0354 (0.1176) | +0.0303 (0.1035) | 0.01249 | 0.2070 | 0.0260 | [`followup-organ-17-heart-gao-p1-signature.csv`](../results/followup-organ-17-heart-gao-p1-signature.csv) |
| Adipocyte | 46 | 60 / 25 | -0.0992 (0.1212) | -0.0051 (0.9288) | -0.0941 (0.0994) | 0.01299 | 0.2070 | 0.0260 | [`followup-organ-17-heart-gao-p1-signature.csv`](../results/followup-organ-17-heart-gao-p1-signature.csv) |
| Cardiomyocyte | 115 | 60 / 29 | +0.0384 (0.0893) | +0.0107 (0.6110) | +0.0277 (0.1324) | 0.12394 | 0.2316 | 0.1928 | [`followup-organ-17-heart-gao-p1-signature.csv`](../results/followup-organ-17-heart-gao-p1-signature.csv) |

---

## Sensitivity Analyses (`sensitivities-*.csv`)

- **10x-only (`N = 108`, dropping Drop-seq `Penn`), Studies `N >= 8` (`N = 110`, dropping `Kuppe 2022` and `Sim 2021`), and Five Standalone Cohorts (`ENCODE v4 N = 54`, `Koenig 2022 N = 23`, `Chaffin 2022 N = 16`, `Penn N = 9`, `Simonson 2023 N = 8`):**
  - **Family A fractions ([`followup-organ-17-heart-gao-sensitivities-fractions.csv`](../results/followup-organ-17-heart-gao-sensitivities-fractions.csv)):** 0 CLR shifts pass `q_bh_family_a < 0.05` in any sensitivity or single-study cohort (`min q_bh = 0.0812` in `Penn Pericyte`).
  - **Family B top genes & concordance ([`followup-organ-17-heart-gao-sensitivities-genes-summary.csv`](../results/followup-organ-17-heart-gao-sensitivities-genes-summary.csv), [`followup-organ-17-heart-gao-sensitivities-concordance.csv`](../results/followup-organ-17-heart-gao-sensitivities-concordance.csv)):** `Pericyte` and `All_LV_Nuclei` hits persist in `10x_only` (`N = 108`) and standalone `ENCODE v4` (`N = 54`). Cross-cell-type concordance stays null across all subsets (`10x_only` `r = 0.1900, p = 0.0919`; `ENCODE v4` `r = 0.2598, p = 0.0849`; `Koenig` `r = 0.1038, p = 0.5544`).
  - **Family C `P1_Clean` replication ([`followup-organ-17-heart-gao-sensitivities-p1-signature.csv`](../results/followup-organ-17-heart-gao-sensitivities-p1-signature.csv)):** Positive `P1_Clean` enrichment replicates in `10x_only` (`All_LV_Nuclei`, `Myeloid`, `Endocardial` `q_emp = 0.0070`; `Cardiomyocyte`, `Pericyte` `q_emp = 0.0084`), `studies_ge8` (`Endocardial` `q_comp = 0.0245, q_emp = 0.0035`; `All_LV_Nuclei`, `Myeloid`, `Fibroblast` `q_emp = 0.0035`), standalone `ENCODE v4` (`Myeloid` `+0.0577`, `Endocardial` `+0.0853`, both `q_emp = 0.0210`), and standalone `Koenig 2022` (`Endothelial` `+0.0543`, `All_LV_Nuclei` `+0.0309`, both `q_emp = 0.0035`; `Fibroblast` `q_emp = 0.0070`).

## Limitations

1. Only left ventricle (`LV`) adult non-diseased nuclei had exact ages across the 7 eligible studies (`Litvinukova`, `Kanemaru`, and `Reichart` encode bracket midpoints and were excluded).
2. `ENCODE v4 (Snyder)` contributes 54/117 donors (`46%`) on 10x Multiome-v1 RNA, while `Penn` uses Drop-seq with low cardiomyocyte capture and `Kuppe`/`Sim` contribute only 4 and 3 donors respectively; sensitivity analyses without `Penn` and without small studies show the findings do not hinge on either subset.
3. Because `P1_Clean` enrichment in LV comes predominantly from negative aging slopes in `P1_DOWN` genes plus `Endocardial` `P1_UP`, rather than uniform multi-lineage `P1_UP` induction, it contrasts modestly with the paper's brain/kidney/lung null without invalidating cell-type heterogeneity across organs.
