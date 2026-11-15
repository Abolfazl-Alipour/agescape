# GTEx v8 Open-Portal Follow-up Report (Bulk RNA-seq, Decade-Binned Age)

## 1. Access Scope, Data Verification, and Methodological Rules

- **Open-portal access only (no dbGaP application)**: Uses strictly open-access GTEx v8 files (`GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz`, `GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt`, and `GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt`). No protected-access phenotypes or exact chronological ages were requested or accessed.
- **File and integer-count verification (`results/followup-gtex-v8-summary.json`)**:
  - Archive size: `917,091,994` bytes; MD5 base64 `HW7Jcxikd9f9grpyFDkl7g==` (`md5_verified = true`).
  - Integer count gate across the full uncompressed GCT stream: `56,200` gene rows x `17,382` RNASEQ sample columns, `0` nonfinite values, `0` negative values, `0` noninteger lines, and `515,414,926` nonzero count entries.
- **Subject and sample structure (`results/followup-gtex-v8-tissue-audit.csv`)**:
  - `980` open-access subjects (`653` male, `327` female), of whom `948` unique donors appear across the `17,382` RNASEQ samples covering `54` annotated tissue sites (`SMTSD`). Inside every `SMTSD`, there is strictly `1` sample per donor (`n_donors == n_samples`).
  - `47` primary tissues meet `n_complete >= 50` (`16,632` complete tissue-donor observations; `n_donors` ranges from `85` in `Kidney - Cortex` to `802` in `Muscle - Skeletal`).
  - `2` isolated cell lines (`n_complete >= 50`) are separated into an independent stratum: `Cells - Cultured fibroblasts` (`n = 487 / 504`, `318M / 169F`) and `Cells - EBV-transformed lymphocytes` (`n = 174 / 174`, `112M / 62F`).
  - `5` low-sample sites (`n < 50`) were excluded prior to testing: `Bladder` (`n = 21`), `Cervix - Endocervix` (`n = 10`), `Cervix - Ectocervix` (`n = 9`), `Fallopian Tube` (`n = 9`), and `Kidney - Medulla` (`n = 4`).
- **Bulk RNA-seq analysis rules**:
  - GTEx v8 is bulk tissue (or bulk cultured cell line) RNA-seq. Per repository protocol, **gene-level age tests and `P1_Clean` signature tests ONLY are permitted; NO single-cell fraction/CLR models and NO cross-cell-type concordance analyses were run**, because bulk expression changes conflate within-cell aging with shifts in cell-type proportions.
  - **Coarse age encoding**: Open-access GTEx `AGE` reports `6` decade bins (`20-29`, `30-39`, `40-49`, `50-59`, `60-69`, `70-79`). All continuous linear and negative binomial slopes treat `age_bin_order` (`1.0` through `6.0`) as a per-decade ordinal predictor (`1` unit = `+1` age decade bin).

---

## 2. Explicit Multiple-Testing Families (Defined Before Testing)

All BH corrections were applied within pre-specified disjoint families in `scripts/run_followup_gtex_v8.py`:
1. **Track 1 Family B (OLS per-decade gene slopes across `47` primary tissues)**: `748,735` expressed-gene OLS models adjusted for `sex_male` (when both sexes exist), `log_lib`, `SMRIN`, `SMTSISCH_hr`, `SMEXNCRT`, and collection-center/isolation dummies (`SMCENTER`, `SMNABTCHT`).
2. **Track 1 Family B-cellline (OLS gene slopes across `2` isolated cell lines)**: `26,916` expressed-gene OLS models (`13,387` fibroblasts + `13,529` EBV-transformed lymphocytes).
3. **Track 1 Family C (`P1_Clean` signature across `47` primary tissues)**: `47` tissue-level `P1_Clean` composite t-tests and `47` two-sided empirical null tests (`2,000` 10-decile expression-matched random-gene draws per tissue).
4. **Track 1 Family C-cellline (`P1_Clean` across `2` isolated cell lines)**: `2` cell-line composite and `2` empirical null tests (`2,000` draws).
5. **Track 2 Family B (`PyDESeq2` NB-GLM across `6` focal primary tissues)**: `140,383` negative binomial likelihood-ratio/Wald tests across `Whole Blood`, `Muscle - Skeletal`, `Lung`, `Kidney - Cortex`, `Ovary`, and `Testis`.
6. **Track 2 Family C (`P1_Clean` NB-GLM signature across `6` focal primary tissues)**: `6` composite Wald z-tests (`mean(up_lfc) - mean(down_lfc)`) and `6` two-sided `2,000`-draw 10-decile expression-matched empirical null tests.

---

## 3. Track 1 Results: 49-Group OLS Scan (`log2(CPM + 1)`)

### 3.1 Gene-Level Summary (`results/followup-gtex-v8-track1-genes-summary.csv`, `results/followup-gtex-v8-track1-top-genes.csv`)

- **Primary tissues (`47` tissues, `748,735` tests)**: `86,495 / 748,735` (`11.55%`) pass Family B BH `q_bh_family < 0.05` (`93,950` within individual tissues).
- **Leading primary-tissue age genes (`q_bh_family`)**:
  - `Artery - Aorta` (`n = 431`): `ROBO2` (`effect = +0.7391/dec`, `t = +20.92`, `p = 3.85e-67`, `q_family = 2.88e-61`), `GLRB` (`+0.3518/dec`, `q_family = 5.70e-39`), `VANGL2` (`-0.3864/dec`, `q_family = 4.09e-38`), `TOX2` (`-0.4613/dec`, `q_family = 1.27e-35`), `PCLO` (`+0.3739/dec`, `q_family = 7.70e-34`), `PTCHD4` (`+0.2513/dec`, `q_family = 4.26e-31`), `PCDH10` (`+0.3762/dec`, `q_family = 1.14e-30`).
  - `Nerve - Tibial` (`n = 619`): `EDA2R` (`+0.2539/dec`, `t = +15.78`, `p = 2.86e-47`, `q_family = 1.07e-41`), `CDKN2A` (`+0.2806/dec`, `t = +13.99`, `p = 9.25e-39`, `q_family = 7.70e-34`).
  - `Artery - Tibial` (`n = 663`): `PCDH10` (`+0.4897/dec`, `p = 7.28e-47`, `q_family = 1.82e-41`), `MIR34AHG` (`+0.2439/dec`, `q_family = 5.84e-36`), `EDA2R` (`+0.2379/dec`, `q_family = 4.24e-30`).
  - `Thyroid` (`n = 652`): `PTCHD4` (`+0.2699/dec`, `p = 8.94e-36`, `q_family = 4.26e-31`).
  - `Adipose - Subcutaneous` (`n = 662`): `EYA4` (`+0.2411/dec`, `p = 5.33e-35`, `q_family = 2.22e-30`).
- **Isolated cell lines (`2` lines, `26,916` tests)**: `4,026 / 26,916` pass `q_bh_family < 0.05` (`Cells - Cultured fibroblasts n = 487`: `3,925` in family / `4,874` within-group, top `ZNF518B -0.1460/dec`, `p = 6.42e-13`; `Cells - EBV-transformed lymphocytes n = 174`: `101` in 2-line family / `0` within-group, top `OTUD7A -0.2756/dec`, `p = 2.11e-05`, `q_within = 0.2858`).

### 3.2 `P1_Clean` Across All 47 Primary Tissues and 2 Isolated Cell Lines (`results/followup-gtex-v8-track1-p1-signature.csv`)

- **Primary tissues (`47` tests)**: `13 / 47` tissues pass empirical BH `q_bh_empirical < 0.05` (`10` positive, `3` negative), while only **`6 / 47`** pass **BOTH** parametric `q_bh_composite < 0.05` and matched empirical `q_bh_empirical < 0.05` (`3` positive, `3` negative):

| Group | Class | n | Up/Down Genes | Composite Effect/dec (SE) | Composite p | q_comp | Emp Null p | q_emp |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Nerve - Tibial | primary_tissue | 619 | 67 / 34 | +0.0420 (0.0094) | 8.69e-06 | 0.00041 | 0.00050 | 0.00294 |
| Adrenal Gland | primary_tissue | 258 | 67 / 33 | +0.0523 (0.0156) | 0.00094 | 0.01108 | 0.00050 | 0.00294 |
| Thyroid | primary_tissue | 652 | 69 / 34 | +0.0269 (0.0077) | 0.00048 | 0.00748 | 0.00050 | 0.00294 |
| Adipose - Visceral (Omentum) | primary_tissue | 540 | 69 / 34 | -0.0429 (0.0117) | 0.00029 | 0.00682 | 0.00050 | 0.00294 |
| Colon - Sigmoid | primary_tissue | 373 | 69 / 34 | -0.0399 (0.0129) | 0.00215 | 0.01680 | 0.00050 | 0.00294 |
| Minor Salivary Gland | primary_tissue | 162 | 71 / 35 | -0.0617 (0.0205) | 0.00305 | 0.02048 | 0.00200 | 0.00854 |
| Brain - Caudate (basal ganglia) | primary_tissue | 246 | 69 / 36 | +0.0398 (0.0182) | 0.03003 | 0.15684 | 0.00050 | 0.00294 |
| Brain - Spinal cord (cervical c-1) | primary_tissue | 159 | 68 / 36 | +0.0471 (0.0238) | 0.04975 | 0.21258 | 0.00050 | 0.00294 |
| Brain - Anterior cingulate cortex (BA24) | primary_tissue | 176 | 67 / 35 | +0.0351 (0.0225) | 0.12129 | 0.31669 | 0.00050 | 0.00294 |
| Brain - Cortex | primary_tissue | 255 | 69 / 33 | +0.0328 (0.0194) | 0.09166 | 0.28720 | 0.00100 | 0.00470 |
| Brain - Hypothalamus | primary_tissue | 202 | 68 / 37 | +0.0347 (0.0239) | 0.14759 | 0.36509 | 0.00250 | 0.00979 |
| Brain - Amygdala | primary_tissue | 152 | 68 / 34 | +0.0280 (0.0242) | 0.24957 | 0.53317 | 0.00450 | 0.01626 |
| Heart - Left Ventricle | primary_tissue | 432 | 69 / 31 | +0.0244 (0.0148) | 0.09905 | 0.29096 | 0.00100 | 0.00470 |
| Whole Blood | primary_tissue | 746 | 56 / 24 | +0.0211 (0.0073) | 0.00414 | 0.02435 | 0.06347 | 0.12429 |
| Kidney - Cortex | primary_tissue | 85 | 72 / 37 | +0.0330 (0.0398) | 0.40969 | 0.64184 | 0.01499 | 0.05033 |
| Ovary | primary_tissue | 180 | 66 / 31 | +0.0332 (0.0173) | 0.05680 | 0.22245 | 0.02399 | 0.06632 |
| Muscle - Skeletal | primary_tissue | 802 | 62 / 30 | -0.0170 (0.0100) | 0.08799 | 0.28720 | 0.05547 | 0.11336 |
| Testis | primary_tissue | 361 | 69 / 38 | +0.0073 (0.0066) | 0.26731 | 0.54625 | 0.44178 | 0.53842 |
| Lung | primary_tissue | 578 | 72 / 35 | +0.0011 (0.0118) | 0.92773 | 0.94790 | 0.99150 | 0.99150 |
| Cells - Cultured fibroblasts | isolated_cell_line | 487 | 60 / 29 | -0.0326 (0.0095) | 0.00066 | 0.00132 | 0.00150 | 0.00150 |
| Cells - EBV-transformed lymphocytes | isolated_cell_line | 174 | 56 / 22 | +0.0357 (0.0175) | 0.04331 | 0.04331 | 0.00150 | 0.00150 |

- Key patterns in Track 1:
  - Six CNS regions (`Brain - Caudate`, `Spinal cord`, `BA24`, `Cortex`, `Hypothalamus`, `Amygdala`) and `Heart - Left Ventricle` pass `q_bh_empirical < 0.05` relative to matched background genes yet fail parametric `q_bh_composite < 0.05` (`0.1568` to `0.5332`) because donor-to-donor variance remains high relative to modest composite slopes (`+0.024` to `+0.047/dec`).
  - `Cells - Cultured fibroblasts` (`-0.0326/dec`) and `Cells - EBV-transformed lymphocytes` (`+0.0357/dec`) diverge in opposite directions, underscoring lineage/culture specificity rather than universal cellular drift.

---

## 4. Track 2 Results: 6-Tissue `PyDESeq2` Negative Binomial GLM + Sensitivities

### 4.1 Primary NB-GLM Gene Tests and `P1_Clean` (`results/followup-gtex-v8-track2-nbglm-summary.csv`, `results/followup-gtex-v8-track2-p1-nbglm.csv`, `results/followup-gtex-v8-track2-nbglm-top-genes.csv`)

Across `140,383` total NB-GLM tests in the `6` focal tissues matching established paper organs, `23,759` (`16.92%`) pass global `q_bh_track2_global < 0.05`:

| Tissue | n | df_approx | Tested | q_within < 0.05 (Up/Down) | q_global < 0.05 | NB P1 LFC/dec | Comp Z (Comp p, q) | Emp p | q_emp |
|---|---:|---:|---:|---|---:|---:|---|---:|---:|
| Whole Blood | 746 | 738 | 18,991 | 7,050 (3,345 / 3,705) | 6,007 | +0.0058 | +8.27 (1.30e-16, 7.81e-16) | 0.80960 | 0.93903 |
| Muscle - Skeletal | 802 | 793 | 19,795 | 6,682 (3,655 / 3,027) | 5,747 | -0.0560 | -7.04 (1.95e-12, 5.85e-12) | 0.00050 | 0.00300 |
| Lung | 578 | 569 | 23,680 | 3,141 (1,629 / 1,512) | 3,411 | +0.0029 | +1.00 (0.3187, 0.3187) | 0.93903 | 0.93903 |
| Kidney - Cortex | 85 | 77 | 23,249 | 7 (2 / 5) | 345 | +0.0131 | +2.26 (0.0236, 0.0291) | 0.54473 | 0.93903 |
| Ovary | 180 | 173 | 22,895 | 4,765 (2,323 / 2,442) | 4,509 | +0.0093 | +2.25 (0.0243, 0.0291) | 0.58621 | 0.93903 |
| Testis | 361 | 354 | 31,773 | 2,900 (2,025 / 875) | 3,740 | -0.0013 | +4.57 (4.79e-06, 9.57e-06) | 0.75512 | 0.93903 |

- **Zero positive empirical `P1_Clean` hits (`0 / 6`) in NB-GLM**:
  - `Muscle - Skeletal` shows a strong **negative** `P1_Clean` reversal (`nb_p1_effect_lfc = -0.0560/dec`, `composite_z = -7.038`, `q_bh_composite = 5.85e-12`, `empirical_null_p = 0.00050`, `q_bh_empirical = 0.00300`).
  - `Whole Blood`, `Testis`, `Kidney - Cortex`, and `Ovary` illustrate why matched random-gene permutations are essential: raw parametric composite Wald p-values look small (`1.30e-16` to `0.0243`) due to shared inter-gene correlation / bulk shifts across large `n`, yet `2,000` expression-matched random gene sets show zero empirical specificity (`empirical_null_p = 0.5447` to `0.8096`, `q_bh_empirical = 0.9390`).
  - `Lung` is flat null under both parametric (`p = 0.3187`) and empirical (`p = 0.9390`) tests, corroborating human lung single-cell findings.
- **Leading Track 2 NB-GLM genes (`results/followup-gtex-v8-track2-nbglm-top-genes.csv`)**:
  - `Testis`: `PTCHD4` (`lfc = +0.6523/dec`, `p = 1.36e-59`), `MFAP5` (`+0.3611`, `p = 1.21e-22`), `EDA2R` (`+0.3216`, `p = 1.75e-21`).
  - `Muscle - Skeletal`: `GALNTL6` (`+0.6683`, `p = 2.86e-44`), `UNC13C` (`+0.4845`, `p = 1.47e-39`), `CDKN2A` (`+0.5609`, `p = 6.75e-37`), `CLIC6` (`-0.4723`, `p = 1.09e-35`).
  - `Ovary`: `NLRP5` (`-2.4739`, `p = 1.30e-36`), `PLAG1` (`+0.5066`, `p = 2.13e-35`), `PADI6` (`-3.2132`, `p = 3.01e-31`), predominantly reflecting oocyte/follicle depletion across menopause in bulk tissue.
  - `Whole Blood`: `CD1A` (`-0.5149`, `p = 2.22e-32`), `FIGN` (`+0.5872`, `p = 7.13e-29`), `CCL3` (`+0.5737`, `p = 5.62e-28`), `RAG1` (`-0.6188`, `p = 8.35e-26`).
  - `Lung`: `ITGBL1` (`+0.3246`, `p = 2.17e-24`), `PTCHD4` (`+0.4191`, `p = 2.47e-24`), `EDA2R` (`+0.2923`, `p = 7.30e-23`), `CDKN2A` (`+0.3033`, `p = 2.12e-21`).
  - `Kidney - Cortex`: `PRRT4` (`-0.7562`, `p = 1.75e-06`), `KLK4` (`+0.9396`, `p = 2.08e-06`), `DGAT2` (`-0.4829`, `p = 4.85e-06`), `LPL` (`-0.5015`, `p = 8.83e-06`).

### 4.2 Sensitivity Analyses (`results/followup-gtex-v8-track2-sensitivities.csv`)

Across all `5` frozen Track 2 sensitivity arms:
1. **`oldest_bin_tail_excl_70_79`**: Dropping `70-79` donors (`23` blood, `28` muscle, `22` lung, `4` kidney, `4` ovary, `9` testis) leaves Wald rankings almost identical (`Spearman r = 0.9342` in `Kidney - Cortex` to `0.9903` in `Ovary`).
2. **`read_length_76bp`**: Restricting to `76 bp` reads is virtually identical (`r = 0.9993` to `1.0000`, dropping only `1` blood, `1` muscle, and `2` lung `101 bp` samples).
3. **`blood_signed_ischemic_time`**: Restricting `Whole Blood` (`n = 417`) to `SMTSISCH >= 0` cuts within-tissue DEGs from `7,050` to `3,996` (`r = 0.8474`).
4. **`extended_quality_and_death`**: Adjusting additionally for Hardy death scale (`DTHHRDY`) and autolysis score (`SMATSSCR`, omitted in blood where constant) sharply attenuates bulk DEG counts (`Whole Blood 7,050 -> 2,957`, `r = 0.8473`; `Muscle 6,682 -> 2,913`, `r = 0.8359`; `Lung 3,141 -> 1,734`, `r = 0.9264`; `Ovary 4,765 -> 2,260`, `r = 0.8796`; `Testis 2,900 -> 2,155`, `r = 0.9773`; `Kidney - Cortex 7 -> 1`, `r = 0.6854`), demonstrating strong post-mortem/manner-of-death confounding in unadjusted GTEx bulk age slopes.
5. **`categorical_age_lrt` (`df_diff = 5`)**: Nonlinear/menopause-step tissues such as `Ovary` (`2,394` passing LRT, `r = 0.7309`) and sparse `Kidney - Cortex` (`28` vs `7`, `r = 0.3910`) shift more strongly than mostly monotonic `Muscle` (`r = 0.7683`) or `Whole Blood` (`r = 0.7567`).

---

## 5. Limitations and Takeaway

- **Bulk cell-composition confounding**: Strongest GTEx single-gene age shifts include canonical parenchymal/stromal loss or compositional shifts (`Ovary NLRP5`/`PADI6`, `Whole Blood RAG1`/`CD1A`), so individual bulk hits cannot be attributed to cell-intrinsic aging without single-cell resolution.
- **Coarse decade bins and post-mortem heterogeneity**: `20-29` through `70-79` decade binning drops exact age variance, and `DTHHRDY`/`SMATSSCR` sensitivity shows manner-of-death explains a large fraction of apparent bulk DEGs.
- **Consistent with focal organ selectivity**: Neither `Lung`, `Kidney - Cortex`, `Whole Blood`, `Ovary`, nor `Testis` exhibits an empirically specific increase in the 139-gene murine aging signature (`0 / 6` positive in NB-GLM `P1_Clean`), while `Muscle - Skeletal` reverses significantly negatively (`-0.0560/dec`, `q_emp = 0.0030`), and only `3 / 47` primary tissues (`Nerve - Tibial`, `Adrenal Gland`, `Thyroid`) pass both parametric and empirical FDR positively versus `3 / 47` passing both negatively (`Adipose - Visceral`, `Colon - Sigmoid`, `Minor Salivary Gland`).
