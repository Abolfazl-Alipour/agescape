# HCATA (PMC12511636, Item #17) Reconciliation and Follow-up Execution

## Scope and Primary-Source Audit

This report resolves Item `#17` (Human Cell Age Transcriptomic Atlas, HCATA, `PMC12511636`, Zhang et al.) strictly as a discovery index against primary sources (original publications, GEO records, and CELLxGENE Census metadata) rather than treating Table S1/S2 summary labels (or unverified total cell counts) as ground truth.

Across all `76` rows in HCATA Table S1 (reconciled in `results/followup-hcata-pmc12511636-studies.csv`) and all `11` rows selected by HCATA for differential expression in Table S2 (audited in `results/followup-hcata-pmc12511636-table-s2-audit.csv`):
- **Total Table S1 rows**: `76` (`65` with `selected_for_deg_in_hcata = No`, `11` with `selected_for_deg_in_hcata = Yes`).
- **Ineligible by design / disease / sparsity / developmental age**: `60` rows (`15` disease, `11` sparse donors, `3` disease + sparse, `2` sparse + disease, `2` single donor, `2` single age, `2` cancer-adjacent + sparse, `2` pediatric + disease, `1` cancer-adjacent, `1` cancer-adjacent elderly, `1` pediatric, `1` fetal/pediatric, `1` pediatric/fetal, plus `16` other sparse/pediatric/controlled/pooled/binned categories).
- **Overlaps already analyzed cohorts (`10` rows)**: `4` already done in agescape (`already_done_agescape`), `2` overlapping HLCA (`overlaps_hlca`), `1` `overlaps_hlca_and_already_done`, `1` `already_covered_ts2`, `1` `overlaps_gao_heart`, and `1` `overlaps_hlca_mixed_platforms`.
- **Blocked for missing metadata or ages (`2` rows)**: `Oetjen 2018` (`30518681`, `GSE120221`, bone marrow, selected in Table S2; Matrix Market files only without public per-cell annotations or sex metadata) and `Williams 2021` (`34129837`, `GSE164241`, oral mucosa; no donor ages in GEO sample records).
- **Wrong species misattribution (`1` row)**: `Welch 2019` (`31178122`, Table S2 row 2: GEO `GSE126836` is `15` mouse BNST brain Drop-seq samples + LIGER method paper, not `38` human brain donors).
- **Newly executed open exact-age studies (`3` rows)**:
  1. `Guilliams 2022` (`35021063`, `GSE192742`, liver, Table S2 row 9): `eligible_sorted_execute` (FACS-sorted `CD45+`/`CD45-` enrichment, not unbiased tissue cell counting).
  2. `Yazar 2022` (`35389779`, `GSE196830`, OneK1K PBMC, not selected by HCATA despite top power): `eligible_open_execute`.
  3. `Perez 2022` (`35389781`, `GSE174188`, SLE/control PBMC, not selected by HCATA): `eligible_open_execute`.

### Specific Table S1 / Table S2 Forensic Findings

From `results/followup-hcata-pmc12511636-studies.csv` and `results/followup-hcata-pmc12511636-table-s2-audit.csv`:
- **Young 2018 (`30093597`, row 2)**: All `40,268` Census cells carry `is_primary_data = False` overlapping verbatim with Stewart 2019; macroscopically normal margins came predominantly from kidney tumor resections (`n = 10` tumor resections, `n = 2` transplant), leaving only `4` adult non-cancer transplant kidneys plus pediatric/tumor margins (`< 6` clean adult non-cancer donors).
- **Oetjen 2018 (`30518681`, row 6 / Table S2 row 1)**: `GSE120221` provides Matrix Market files only without public per-cell annotations or sex metadata.
- **Welch 2019 (`31178122`, row 11 / Table S2 row 2)**: Misattributed in HCATA Table S2 as `38` human brain donors (`18-66y`); `GSE126836` is `15` mouse BNST brain samples.
- **Drokhlyansky 2020 (`32888429`, row 24 / Table S2 row 6)**: Broad `SCP1038` contains `10` CRC cancer-adjacent normal colon donors (`35-90y`), not `59` healthy human donors.
- **Williams 2021 (`34129837`, row 41)**: Oral mucosa `GSE164241` lacks donor ages in GEO sample records.
- **Pelka 2021 (`34450029`, row 45 / Table S2 row 8)**: `GSE178341` `36` "normal" colon samples (`35-91y`) are tumor-adjacent normal mucosa from colorectal cancer surgical resections.
- **Madissoon 2023 (`36543915`, row 69 / Table S2 row 11)**: Overlaps HLCA (`290B`, `292B`, `296C`, `298C`), mixes `5` airway regions and `3` platforms across `n = 7` exact donors, and HCATA Table S1 erroneously pasted eye cell types onto this lung row.

---

## Study 1: Guilliams 2022 Liver (`GSE192742`, Sorted-Cell Rules)

### Cohort and Method Constraints

From `results/followup-hcata-guilliams-donors.csv`, the executable clean exact-age subset (`10x 3' v3` single-cell assay) comprises **`n = 10` donors**, aged `28.0` to `77.0` years (median `60.5` years; `8` female, `2` male). All cells in this dataset underwent FACS suspension enrichment (`CD45+`, `CD45-`, or `CD45+;CD14+;CD16+`); per sorted-cell rules, **no Family A CLR composition models and no cross-cell-type concordance inference were run** because sorted fractions distort whole-liver cellular proportions. Clinical comorbidities and resection pathology also mix with age (`H16` age `66` metastatic colon cancer/colorectal cancer; `H21` age `77` cystadenocarcinoma/hepatocellular carcinoma; `H13` age `73` thyroid tumor; `H11` `53F`, `H10` `55M`, `H23` `75F` with obesity/morbid obesity; `H25` `49F` obesity; `H18` `28F` hepatocellular adenoma), and donors span only `2` residual degrees of freedom in rare compartments (`endothelial cell n = 6, df = 3`). Note also that `4` donors (`H18`, `H11`, `H16`, `H13`) have `is_citeseq = 1.0` and `6` have `0.0`.

### Explicit Multiple-Testing Families and Results

- **Family B (Gene-level OLS on `age_per_decade + sex_male + log_lib`, per-cell-type)**: `81,881` gene-group tests across `7` cell types (`results/followup-hcata-guilliams-genes-summary.csv`). **Result: `0 / 81,881` genes pass global BH `q_bh_global < 0.05`** (minimum global `q = 0.6100`; leading nominal hits in `endothelial cell` `HERC6 p = 1.46e-05` with only `n = 6, df = 3` and `T cell COQ5 p = 3.09e-05`, `results/followup-hcata-guilliams-top-genes.csv`). Separate pooled groups (`All_Cells` and `All_CD45pos_Cells`, `n = 10`, `26,220` pooled tests combined): `0` global hits (`min q_bh_pool_only = 0.2058`, `GOT1 p = 7.85e-06`).
- **Family C (`P1_Clean` mouse aging signature, `2,000` matched random-gene draws)**: `7` cell types evaluated (`results/followup-hcata-guilliams-p1-signature.csv`). **Result: `0 / 7` pass BH `q_bh_composite < 0.05` (minimum `q = 0.1888`) and `0 / 7` pass BH `q_bh_empirical < 0.05` (minimum `q = 0.1259`)**. Every cell type and pooled group had a **negative** composite point estimate (`-0.0432` to `-0.1563` per decade; `All_CD45pos_Cells -0.0730, composite_p = 0.1512, empirical_null_p = 0.0085`; `All_Cells -0.0723, composite_p = 0.1614, empirical_null_p = 0.0165`).
- **Disease/Technical Sensitivity Arms (`4` arms, `results/followup-hcata-guilliams-sensitivities-genes.csv` and `results/followup-hcata-guilliams-sensitivities-p1.csv`)**: Excluding malignant resections (`no_malignant_resection_n8`, dropping `H16` and `H21`), excluding any tumor/adenoma (`no_any_tumor_n7`), restricting to females (`female_only_n8`), and adding `is_citeseq` (`adjust_citeseq_n10`) yielded **`0` cell-type gene hits at global `q < 0.05`** across every arm (`1` pooled `All_Cells` hit only in `no_malignant_resection_n8`, `q = 0.0439`, `df = 4`, fragile) and **`0` P1_Clean empirical hits (`min q_bh_empirical >= 0.0501` in `female_only_n8`, `>= 0.1201` in `adjust_citeseq_n10`, `>= 0.1249` in `no_any_tumor_n7`, `>= 0.2049` in `no_malignant_resection_n8`)**.

---

## Study 2: Yazar 2022 OneK1K PBMC (`GSE196830`, Full Framework + 75-Pool Fixed-Effect Sensitivity)

### Cohort and Multiplexing-Pool Confounding

From `results/followup-hcata-yazar-donors.csv`, the clean normal cohort contains **`n = 981` donors**, aged `19.0` to `97.0` years (median `67.0` years; `565` female, `416` male; `976` donors aged `>= 20` years), across `75` Chromium multiplexing pools (`pool_number`, `9-14` donors per pool). Critically, **donor age is strongly confounded with `pool_number` (one-way ANOVA `F = 3.7473`, `R^2 = 0.2343`, `p = 3.814e-21`)**. Azimuth doublets (`131` cells annotated under `"peripheral blood mononuclear cell"`) were dropped before CLR composition so only `28` bona fide blood cell types entered Family A.

Because `pool_number` explains `23.4%` of age variance, we report the primary age+sex+log_lib specification (`results/followup-hcata-yazar-*.csv`, plus `age_ge_20_n976`) **alongside a within-pool fixed-effect sensitivity (`pool_fe`, `results/followup-hcata-pool-sensitivity-*.csv`)** that restricts each cell type to donors in pools containing `>= 2` donors (`16` cell types plus `All_Cells` with `55-903` residual degrees of freedom; `3` ultra-rare cell types with `<= 1` within-pool df marked `not_evaluable_within_pool_df` in `results/followup-hcata-yazar-pool-fe-eligibility.csv`).

### Explicit Multiple-Testing Families and Results

1. **Family A (CLR composition, `28` non-doublet cell types across `981` donors)**:
   - **Primary (`results/followup-hcata-yazar-fractions.csv`)**: `21 / 28` cell types pass BH `q < 0.05` (`23 / 28` in `age_ge_20_n976`). Leading declines: `naive thymus-derived CD8-positive, alpha-beta T cell` (`effect_clr = -0.3787/dec`, `p = 8.90e-107`, `q = 2.49e-105`), `mucosal-associated invariant T cell` (`-0.1605/dec`, `q = 3.63e-16`), and `naive thymus-derived CD4-positive, alpha-beta T cell`; leading increases: `natural killer cell` (`+0.1296/dec`, `q = 7.70e-30`), `CD4-positive, alpha-beta cytotoxic T cell` (`+0.1913/dec`, `q = 3.75e-17`), and `effector memory CD8-positive, alpha-beta T cell` (`+0.1207/dec`, `q = 1.99e-15`).
   - **75-Pool Fixed Effects (`results/followup-hcata-pool-sensitivity-fractions.csv`)**: All **`21 / 28`** CLR shifts persist strongly (`df = 904`): `naive CD8+ T` (`-0.3827/dec`, `q = 3.11e-92`), `NK` (`+0.1028/dec`, `q = 4.46e-17`), `MAIT` (`-0.1622/dec`, `q = 1.05e-13`), `CD14-low, CD16-positive monocyte` (`+0.1596/dec`, `q = 2.61e-13`), `CD4 CTL` (`+0.1558/dec`, `q = 1.04e-10`), `regulatory T cell` (`+0.0673/dec`, `q = 1.04e-10`), `effector memory CD8+ T` (`+0.1099/dec`, `q = 1.82e-10`), `naive CD4+ T` (`-0.0784/dec`, `q = 8.80e-08`). Composition in PBMC robustly shifts away from naive/MAIT T cells toward NK, cytotoxic/effector T, Treg, and non-classical monocytes.

2. **Family B (Gene-level OLS and cross-cell-type concordance)**:
   - **Primary (`results/followup-hcata-yazar-genes-summary.csv`, `results/followup-hcata-yazar-top-genes.csv`, `results/followup-hcata-yazar-concordance.csv`)**: `10,123 / 160,239` cell-type gene tests pass global BH `q < 0.05` across `19` cell types (top: `naive CD4+ T SOX4`, `effect = -0.6341`, `p = 1.15e-67`, `q = 1.84e-62`; `LRRN3`, `q = 2.50e-52`), plus `7,925 / 13,103` in `All_Cells` (`LRRN3 q_pool_only = 5.37e-98`). Unadjusted concordance across `19` cell types (`171` pairs, `5,127` shared genes): `mean_r = +0.2430`, `median_r = +0.2259`, `perm_p = 0.000999` (`+0.2361` in `age_ge_20_n976`).
   - **75-Pool Fixed Effects (`results/followup-hcata-pool-sensitivity-genes-summary.csv`, `results/followup-hcata-pool-sensitivity-concordance.csv`)**: Adjusting for `pool_number` cuts cell-type DEG hits by **~83%**, from `10,123 / 160,239` to **`1,741 / 134,036`** across the `16` evaluable cell types (`All_Cells` drops from `7,925` to `6,241 / 13,103`), and cuts cross-cell-type concordance by **more than half**, to **`mean_r = +0.1024`** (`median_r = +0.0828`, `120` pairs across `16` cell types, `5,382` common genes, `perm_p = 0.000999`).

3. **Family C (`P1_Clean` mouse aging signature, `2,000` random draws)**:
   - **Primary (`results/followup-hcata-yazar-p1-signature.csv`)**: Whole-PBMC `All_Cells` (`n = 981`) is small (`+0.0086/dec`, `composite_p = 0.0032`) and **fails empirical gene-set permutation** (`empirical_null_p = 0.3323`; `+0.0084/dec, emp_p = 0.3413` in `age_ge_20_n976`). Within cell types, **`9 / 19`** pass `q_bh_empirical < 0.05` (and `11 / 19` pass `q_bh_composite < 0.05`, effects `+0.0253` to `+0.0833/dec`).
   - **75-Pool Fixed Effects (`results/followup-hcata-pool-sensitivity-p1-signature.csv`)**: Controlling `pool_number` drops `All_Cells` to **completely null** (`+0.0053/dec`, `composite_p = 0.0734`, `empirical_null_p = 0.5287`) and eliminates `4` unadjusted hits (leaving **`5 / 16`** cell types at `q_bh_empirical < 0.05`, attenuated slopes `+0.0168` to `+0.0379/dec`):
     - `central memory CD4-positive, alpha-beta T cell` (`n = 980, df = 902`): `+0.0230/dec`, `composite_p = 3.38e-06`, `empirical_null_p = 0.0015`, `q_bh_empirical = 0.0240`.
     - `naive thymus-derived CD4-positive, alpha-beta T cell` (`n = 974, df = 896`): `+0.0269/dec`, `composite_p = 1.41e-05`, `empirical_null_p = 0.0045`, `q_bh_empirical = 0.0360`.
     - `memory B cell` (`n = 614, df = 536`): `+0.0379/dec`, `composite_p = 0.0101`, `empirical_null_p = 0.0075`, `q_bh_empirical = 0.0400`.
     - `naive thymus-derived CD8-positive, alpha-beta T cell` (`n = 671, df = 593`): `+0.0364/dec` (halved from unadjusted `+0.0833`), `composite_p = 0.0199`, `empirical_null_p = 0.0110`, `q_bh_empirical = 0.0440`.
     - `natural killer cell` (`n = 977, df = 899`): `+0.0168/dec`, `composite_p = 0.0027`, `empirical_null_p = 0.0145`, `q_bh_empirical = 0.0464`.
     - *(By contrast, `naive B cell` drops from `+0.0509, q_emp = 0.0024` to `+0.0286, q_emp = 0.0594`, no longer passing.)*

---

## Study 3: Perez 2022 SLE/Control PBMC (`GSE174188`, Controls Only + Processing-Cohort Sensitivity)

### Cohort and Processing-Batch Confounding

Excluding `777,258` SLE patient cells (`486,418` normal cells retained) yields **`n = 99` healthy control donors**, aged `21.0` to `75.0` years (median `32.0` years; `97` female, `2` male). **Donor age is strongly confounded with `Processing_Cohort`** (`8` unique strings across cohorts `1.0-4.0`: one-way ANOVA `F = 4.0376`, `R^2 = 0.2370`, `p = 6.800e-04`). We therefore pair the unadjusted primary model (`primary` + `female_only_n97`) with two batch-adjusted sensitivity parameterizations (`results/followup-hcata-pool-sensitivity-*.csv`): **`cohort_multihot`** (indicators for presence in processing cohorts `1.0`, `2.0`, `3.0`, `4.0`, `n = 99`) and **`cohort_level_ge2`** (factor dummies for exact `Processing_Cohort` strata with `>= 2` donors, `n = 98`).

### Explicit Multiple-Testing Families and Results

1. **Family A (CLR composition, `11` cell types)**:
   - **Primary (`results/followup-hcata-perez-fractions.csv`)**: `2 / 11` pass BH `q < 0.05` (`plasmacytoid dendritic cell -0.1734/dec, q = 0.0026`; `CD8+ alpha-beta T -0.0738/dec, q = 0.0436`; same `2 / 11` in `female_only_n97`).
   - **Processing-Cohort Adjusted (`results/followup-hcata-pool-sensitivity-fractions.csv`)**: **Both collapse (`0 / 11` at `q < 0.05`)** under both `cohort_multihot` (`min q = 0.3710`, `non-classical mono p = 0.0514`, `pDC p = 0.1731`) and `cohort_level_ge2` (`min q = 0.5071`).

2. **Family B (Gene-level OLS and cross-cell-type concordance)**:
   - **Primary (`results/followup-hcata-perez-genes-summary.csv`, `results/followup-hcata-perez-concordance.csv`)**: `10,021 / 94,539` cell-type gene tests pass global BH `q < 0.05` (`CD4 T SOX4 effect = -0.6323, p = 3.42e-21, q = 3.23e-16`), and `4,785 / 12,421` in `All_Cells`, with inflated cross-cell-type concordance `mean_r = +0.4519` (`9` cell types, `36` pairs, `7,549` genes, `perm_p = 0.000999`).
   - **Processing-Cohort Adjusted (`results/followup-hcata-pool-sensitivity-genes-summary.csv`, `results/followup-hcata-pool-sensitivity-concordance.csv`)**: Adjusting for `Processing_Cohort` eliminates **>96% of cell-type DEGs**, collapsing `10,021` to **`391 / 94,539`** (`cohort_multihot`) and **`295 / 94,837`** (`cohort_level_ge2`), and collapsing `All_Cells` from `4,785` to **`46 / 12,421`** (`multihot`) and **`36 / 12,445`** (`level_ge2`). Cross-cell-type concordance plunges from `+0.4519` to **`mean_r = +0.2011` (`perm_p = 0.0280`)** under `cohort_multihot` and **`mean_r = +0.1775` (`perm_p = 0.0490`)** under `cohort_level_ge2`.

3. **Family C (`P1_Clean` mouse aging signature, `2,000` random draws)**:
   - **Primary (`results/followup-hcata-perez-p1-signature.csv`)**: **`0 / 9`** cell types pass `q_bh_empirical < 0.05` (`min q = 0.0787`; `8 / 9` negative composite effects, `classical monocyte -0.0355/dec, q_composite = 0.0299, emp_p = 0.0175`; `All_Cells -0.0410/dec, emp_p = 0.0045`).
   - **Processing-Cohort Adjusted (`results/followup-hcata-pool-sensitivity-p1-signature.csv`)**: **`0 / 9`** pass either `q_bh_composite < 0.05` or `q_bh_empirical < 0.05` (`cohort_multihot min q_emp = 0.2609`; `cohort_level_ge2 min q_emp = 0.2316`; `All_Cells` remains negative `-0.0263/dec, emp_p = 0.0075` and `-0.0234/dec, emp_p = 0.0110`).

---

## Comparison Against Established Paper Claims and Limits

- **Consistency with `report/manuscript-aging-mechanisms.md`**:
  1. **No whole-tissue PBMC P1_Clean increase**: Whole-PBMC `All_Cells` fails in all three blood analyses (`Yazar` pool-FE `+0.0053/dec, emp_p = 0.5287`; `Perez` cohort-adjusted negative `-0.0263/dec`; and prior Thomson PBMC `GSE214546`), consistent with the paper's conclusion against a single whole-tissue transported aging program.
  2. **Critical technical-batch dependence**: `report/manuscript-aging-mechanisms.md` explicitly notes that unadjusted public human cohorts lack full technical covariates. Here, controlling multiplexing pool (`Yazar`: `R^2 = 0.234`) and processing cohort (`Perez`: `R^2 = 0.237`) cuts DEGs by `83%` (`Yazar`) and `96-97%` (`Perez`), halves cross-cell-type concordance (`+0.243 -> +0.102` in `Yazar`; `+0.452 -> +0.178-0.201` in `Perez`), and removes `4 / 9` apparent `P1_Clean` empirical hits in `Yazar`.
  3. **Important cell-type caveat to flag for the morning checklist**: Unlike batch-adjusted Thomson PBMC (`0 / 13`) and smaller female-skewed Perez PBMC (`n = 99`, `0 / 9`), the much larger `Yazar` OneK1K PBMC cohort (`n = 981`) retains **`5 / 16` modest positive `P1_Clean` empirical hits** after 75-pool fixed effects (`+0.0168` to `+0.0379` per decade in central memory CD4, naive CD4, memory B, naive CD8, and NK cells). These within-pool PBMC cell-type effects are small in magnitude (and negative in Perez and Guilliams liver) and depend on large sample size (`n = 614-980`), so they sharpen rather than overturn the paper's conclusions, but should be cited explicitly if discussing human PBMC replication.
- **Limits**: Guilliams liver is small (`n = 10`, `2M/8F`), FACS-sorted, and heavily confounded by cancer resection and obesity; Perez PBMC is `97/99` female with median age `32` and strong batch-age confounding; Yazar OneK1K pools explain nearly a quarter of age variance (`residual pool effects` vs unmeasured donor-level cytometry/CMV covariates cannot be excluded).
