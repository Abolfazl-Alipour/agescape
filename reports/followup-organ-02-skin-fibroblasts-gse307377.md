# Dataset 02: Primary Human Dermal Fibroblast Bulk RNA-seq (GSE307377)

## 1. Study Verification and Data Provenance
- **Accession**: `GSE307377` (`GSM9222986` to `GSM9222994`).
- **Assay and Tissue**: Cultured primary human dermal fibroblasts (SD-NSC Human Cell Models of Aging Core; DMEM + 10% FBS for 48h prior to polyA bulk RNA-seq on NovaSeq 6000 SP, HOMER raw counts).
- **Verified Donor Cohort (`n = 9`, `5F / 4M`)**:
  - Young (`n = 4`, ages `23 to 33` yr): `23F, 26F, 23M, 33M`.
  - Old (`n = 5`, ages `64 to 72` yr): `72M, 70F, 69F, 64F, 64M`.
- **Technical Filtering**: Input genes `27,142` -> unique non-empty symbols `27,142` (`73187` fractional expected/HOMER entries, `0` negative) -> `26,920` non-technical genes (`222` `MT-/RP[LS]/MRP[LS]/sex-linked/ERCC` genes removed).

## 2. Scope Rules: Bulk Cultured Lineage Only
- **Family A (CLR cell-type fractions)**: **NOT RUN**. Bulk RNA-seq of a purified, in vitro cultured fibroblast monoculture measures one cell lineage and cannot estimate in vivo skin cell-type proportions.
- **Cross-cell-type concordance**: **NOT RUN**. Requires multiple distinct unbiased cell types from intact tissue.
- **Family B (Donor-level gene age test)**: Bulk OLS on `log2(CPM + 1)` adjusting for `sex_male + log_lib` (`residual_df = 5`), BH correction across the explicit family of `M_B = 12,695` expressed genes (`CPM >= 1.0` in `>= 5` of 9 donors). Both continuous `age_per_decade` (primary) and binary `old_vs_young` (sensitivity) were tested.
- **Family C (`P1_Clean` 139-gene mouse aging signature)**: `K_C = 1` group tested analytically and against `200` expression-decile-matched random gene sets (`matched_null_effect_pvalue`).

## 3. Results

| GSM | Age (yr) | Group | Sex | Non-Tech Lib Size |
|---|---:|---|---|---:|
| GSM9222986 | 23 | young | F | 26,970,200 |
| GSM9222987 | 26 | young | F | 28,201,761 |
| GSM9222988 | 72 | old | M | 21,738,520 |
| GSM9222989 | 70 | old | F | 31,009,748 |
| GSM9222990 | 69 | old | F | 20,965,942 |
| GSM9222991 | 23 | young | M | 24,780,320 |
| GSM9222992 | 64 | old | F | 23,638,859 |
| GSM9222993 | 33 | young | M | 23,012,347 |
| GSM9222994 | 64 | old | M | 26,744,536 |

### Family B: Bulk Dermal Fibroblast Gene Age OLS
- **Primary (`age_per_decade + sex_male + log_lib`, `M_B = 12,695`)**: `0` genes passed `q_BH < 0.05` (minimum `p = 1.33e-05`, minimum `q_BH = 0.1575`).
- **Sensitivity (`old_vs_young + sex_male + log_lib`, `M_B = 12,695`)**: `0` genes passed `q_BH < 0.05` (minimum `p = 1.07e-05`, minimum `q_BH = 0.1356`).

### Family C: Mouse `P1_Clean` Signature Check
- **Primary (`age_per_decade + sex_male + log_lib`)**: `59` UP / `28` DOWN detected; UP slope `+0.0121` (`p = 0.2388`), DOWN slope `+0.0169` (`p = 0.1458`), Composite slope `-0.0048` (`analytic p = 0.7729`, 200-draw random null mean `-0.0016` +/- `0.0167`, **`empirical_null_pvalue = 0.8408`**, `q_BH = 0.8408`).
- **Sensitivity (`old_vs_young + sex_male + log_lib`)**: Composite `-0.0022` (`analytic p = 0.9749`, **`empirical_null_pvalue = 0.9552`**, `q_BH = 0.9552`).

## 4. Limitations and Verdict
- **Culture and cohort limits**: Small donor count (`n = 9`), binary age sampling (`23 to 33` vs `64 to 72` with no middle-aged donors), and short-term in vitro passaging/48h DMEM culture can attenuate or alter in vivo aging signals.
- **Verdict**: Null transport for the mouse `P1_Clean` aging signature in human primary dermal fibroblasts (`empirical p = 0.8408` continuous; `0.9552` binary).

## 5. Output Artifacts
- [Donors](../results/followup-organ-02-skin-fibroblasts-gse307377-donors.csv)
- [Gene Summary](../results/followup-organ-02-skin-fibroblasts-gse307377-gene-summary.csv)
- [Top Genes](../results/followup-organ-02-skin-fibroblasts-gse307377-top-genes.csv)
- [P1 Signature](../results/followup-organ-02-skin-fibroblasts-gse307377-p1-signature.csv)
