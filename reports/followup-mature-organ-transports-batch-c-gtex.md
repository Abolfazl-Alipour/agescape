# Item 5 Batch C: GTEx v8 mature (30-plus bins) transport rerun

## 1. Scope, input verification, and frozen rules

GTEx v8 open-access gene-read counts only, no dbGaP application, no exact ages. File verified before any analysis: GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz downloaded from https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz ; size 917091994 bytes (expected 917091994); MD5 base64 HW7Jcxikd9f9grpyFDkl7g== (expected HW7Jcxikd9f9grpyFDkl7g==); md5_verified = true. SHA-256 148ab1c84609608a00772b0e1431a1f87d5bf40b63071997510181d0a684a110. Count gate on the full uncompressed stream: 56200 gene rows x 17382 RNASEQ sample columns, 0 nonfinite, 0 negative, 0 noninteger lines, 515414926 nonzero entries.

Age restriction: GTEx reports only ordinal decade bins 20-29, 30-39, 40-49, 50-59, 60-69, 70-79. Subject counts by bin verified from GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt: 20-29: 84, 30-39: 78, 40-49: 153, 50-59: 315, 60-69: 317, 70-79: 33, total 980 subjects. The rerun uses bins 30-39 through 70-79 only, with age_bin_order kept at its original values 2.0 to 6.0 so each slope stays per decade bin. No exact age was imputed. Donors in the 20-29 bin are excluded before count aggregation, detection filtering, and model fitting in both tracks.

Mature signatures from Item 1 (commit d3dc14c), loaded and asserted by load_verified_mature_programs: P1m 1,914 mouse genes (1,119 up, 795 down), 1,724 human-mapped (1,013 up, 711 down); P0m 2,598 mouse genes (1,474 up, 1,124 down), 2,049 human-mapped (1,157 up, 892 down). Seed 20261105. 2000-draw expression-decile-matched null as primary. Donor is the inferential unit, never the cell. Bulk RNA-seq only: gene-level age tests and signature tests only, no fraction/CLR models and no cross-tissue concordance, exactly as the original GTEx report states.

Multiple-testing families (declared before testing):
1. T1B-primary: all mature gene OLS p-values across the 47 primary tissues.
2. T1B-cellline: all mature gene OLS p-values across the 2 isolated cell lines.
3. T1C-primary-P1m and T1C-primary-P0m: 47 composite p-values and 47 empirical p-values each.
4. T1C-cellline-P1m and T1C-cellline-P0m: 2 composite and 2 empirical p-values each.
5. T2B-global: all mature NB-GLM gene p-values across the 6 focal tissues.
6. T2C-P1m and T2C-P0m: 6 composite Wald p-values and 6 empirical p-values each.
BH by bh_checked (explicit BH cross-checked against statsmodels, assert diff at most 1e-12); the independent recomputation is in results/followup-mature-organ-transports-c-gtex-bh-check.csv.

Verdict rule: verdict_standard from the Batch A script, applied to the original versus mature both-pass (transport) or any-family-hit (genes) status. SURVIVED means the original signal is still significant at mature q below 0.05; DIED means the original signal is gone; NEW means a mature-only signal; NULL means neither side is significant; UNTESTABLE means the row could not be tested at 30-plus. Note on the task wording: the Batch A partial report claimed 48 of 68 Gao rows have a blank verdict, but the Batch A file results/followup-mature-organ-transports-a-gao.csv as read here contains zero blank verdicts (48 NULL, 20 DIED), so no blank verdicts were invented here and the Batch A verdict_standard rule was used unchanged.

## 2. Track 1 transport (OLS, P1m and P0m)

Under P1m, 5 of 47 primary tissues pass both parametric and empirical q below 0.05, and 0 of 2 isolated cell lines pass both.
  Both-pass P1m primary: Artery - Aorta: effect -0.0291 per decade bin, q_comp 6.079e-05, q_emp 0.0007829, verdict NEW (original class COMPOSITE_ONLY).
  Both-pass P1m primary: Artery - Tibial: effect -0.0300 per decade bin, q_comp 1.295e-05, q_emp 0.0007829, verdict NEW (original class FAIL_BOTH).
  Both-pass P1m primary: Brain - Hypothalamus: effect -0.0432 per decade bin, q_comp 7.578e-05, q_emp 0.0007829, verdict NEW (original class EMPIRICAL_ONLY).
  Both-pass P1m primary: Esophagus - Muscularis: effect -0.0146 per decade bin, q_comp 0.04323, q_emp 0.0007829, verdict NEW (original class FAIL_BOTH).
  Both-pass P1m primary: Pituitary: effect -0.0218 per decade bin, q_comp 0.02152, q_emp 0.0007829, verdict NEW (original class FAIL_BOTH).
  Verdicts P1m Track 1 transport (49 rows): {'NULL': 36, 'DIED': 8, 'NEW': 5}.

Under P0m, 6 of 47 primary tissues pass both parametric and empirical q below 0.05, and 0 of 2 isolated cell lines pass both.
  Both-pass P0m primary: Artery - Aorta: effect -0.0292 per decade bin, q_comp 1.527e-05, q_emp 0.0007829, verdict NEW (original class COMPOSITE_ONLY).
  Both-pass P0m primary: Artery - Tibial: effect -0.0289 per decade bin, q_comp 1.325e-05, q_emp 0.0007829, verdict NEW (original class FAIL_BOTH).
  Both-pass P0m primary: Brain - Hypothalamus: effect -0.0424 per decade bin, q_comp 8.12e-05, q_emp 0.0007829, verdict NEW (original class EMPIRICAL_ONLY).
  Both-pass P0m primary: Brain - Nucleus accumbens (basal ganglia): effect -0.0256 per decade bin, q_comp 0.03762, q_emp 0.0007829, verdict NEW (original class FAIL_BOTH).
  Both-pass P0m primary: Esophagus - Muscularis: effect -0.0141 per decade bin, q_comp 0.04072, q_emp 0.0007829, verdict NEW (original class FAIL_BOTH).
  Both-pass P0m primary: Pituitary: effect -0.0215 per decade bin, q_comp 0.01989, q_emp 0.0007829, verdict NEW (original class FAIL_BOTH).
  Verdicts P0m Track 1 transport (49 rows): {'NULL': 35, 'DIED': 8, 'NEW': 6}.

Original Track 1 dual-pass tissues (3 positive, 3 negative) and their mature verdicts:
  Nerve - Tibial (positive originally), P1m: mature effect -0.0027, q_comp 0.5823, q_emp 0.2221, class FAIL_BOTH, verdict DIED.
  Nerve - Tibial (positive originally), P0m: mature effect -0.0019, q_comp 0.7257, q_emp 0.3326, class FAIL_BOTH, verdict DIED.
  Adrenal Gland (positive originally), P1m: mature effect -0.0107, q_comp 0.1372, q_emp 0.0007829, class EMPIRICAL_ONLY, verdict DIED.
  Adrenal Gland (positive originally), P0m: mature effect -0.0094, q_comp 0.1624, q_emp 0.0007829, class EMPIRICAL_ONLY, verdict DIED.
  Thyroid (positive originally), P1m: mature effect -0.0071, q_comp 0.132, q_emp 0.0007829, class EMPIRICAL_ONLY, verdict DIED.
  Thyroid (positive originally), P0m: mature effect -0.0080, q_comp 0.07376, q_emp 0.0007829, class EMPIRICAL_ONLY, verdict DIED.
  Adipose - Visceral (Omentum) (negative originally), P1m: mature effect -0.0006, q_comp 0.9493, q_emp 0.7762, class FAIL_BOTH, verdict DIED.
  Adipose - Visceral (Omentum) (negative originally), P0m: mature effect -0.0012, q_comp 0.8957, q_emp 0.4633, class FAIL_BOTH, verdict DIED.
  Colon - Sigmoid (negative originally), P1m: mature effect -0.0028, q_comp 0.7921, q_emp 0.07474, class FAIL_BOTH, verdict DIED.
  Colon - Sigmoid (negative originally), P0m: mature effect -0.0034, q_comp 0.7257, q_emp 0.02562, class EMPIRICAL_ONLY, verdict DIED.
  Minor Salivary Gland (negative originally), P1m: mature effect +0.0007, q_comp 0.9493, q_emp 0.9907, class FAIL_BOTH, verdict DIED.
  Minor Salivary Gland (negative originally), P0m: mature effect +0.0007, q_comp 0.9443, q_emp 0.96, class FAIL_BOTH, verdict DIED.

## 3. Track 1 gene-level OLS

Across 49 tested tissues, 776180 mature gene OLS tests, 36005 pass the mature family q below 0.05 (primary family and cell-line family separately). Original totals were 748,735 primary tests with 86,495 family-significant and 26,916 cell-line tests with 4,026 family-significant. Verdicts across the 49 genes_track1 rows: {'SURVIVED': 49}.
  Tissues whose gene-level family signal DIED at 30-plus: 0: none.
  Tissues with a NEW gene-level family signal at 30-plus: 0: none.

## 4. Track 2 NB-GLM (6 focal tissues)

Across 6 tested tissues, 140426 mature NB-GLM gene tests, 8920 pass the mature global q below 0.05. Original totals were 140,383 tests with 23,759 global-significant. Verdicts across the 6 genes_track2 rows: {'SURVIVED': 6}.
Under P1m, 2 of 6 focal tissues pass both composite and empirical q below 0.05. Verdicts P1m Track 2 transport (6 rows): {'NULL': 3, 'NEW': 2, 'DIED': 1}.
  Both-pass P1m Track 2: Whole Blood: effect +0.0119 LFC per decade bin, q_comp 1.289e-22, q_emp 0.001499, verdict NEW (original class COMPOSITE_ONLY).
  Both-pass P1m Track 2: Testis: effect -0.0144 LFC per decade bin, q_comp 1.887e-05, q_emp 0.001499, verdict NEW (original class COMPOSITE_ONLY).
  Muscle - Skeletal, P1m: mature effect -0.0019 LFC per decade bin, q_comp 0.1287, q_emp 0.6577, class FAIL_BOTH, verdict DIED (original P1 was a negative both-pass: effect -0.0560, q_comp 5.85e-12, q_emp 0.00300).
Under P0m, 2 of 6 focal tissues pass both composite and empirical q below 0.05. Verdicts P0m Track 2 transport (6 rows): {'NULL': 3, 'NEW': 2, 'DIED': 1}.
  Both-pass P0m Track 2: Whole Blood: effect +0.0122 LFC per decade bin, q_comp 2.453e-20, q_emp 0.001499, verdict NEW (original class COMPOSITE_ONLY).
  Both-pass P0m Track 2: Testis: effect -0.0132 LFC per decade bin, q_comp 1.115e-06, q_emp 0.001499, verdict NEW (original class COMPOSITE_ONLY).
  Muscle - Skeletal, P0m: mature effect -0.0015 LFC per decade bin, q_comp 0.201, q_emp 0.95, class FAIL_BOTH, verdict DIED (original P1 was a negative both-pass: effect -0.0560, q_comp 5.85e-12, q_emp 0.00300).

## 5. What was not done

No cross-tissue concordance: the original GTEx analysis ran none for bulk tissue, and none was run here. No Track 2 sensitivity arms (categorical-age LRT, oldest-bin exclusion, extended quality/death covariates, signed ischemic time, read-length restriction): the item asked for the primary per-gene and transport reruns at 30-plus, and the sensitivities were out of scope. No exact ages were imputed and no protected GTEx data were accessed. Tissues whose mature covariate-complete n falls below 50 would be UNTESTABLE; the per-tissue audit is in the Track 1 genes summary CSV.

Reading the empirical null: P1m (1,724 human-mapped genes) and P0m (2,049) are about 12 to 15 times larger than the original 139-gene P1_Clean, so the expression-matched random-gene null is much narrower and a small bulk shift can pass it where the same shift failed against the P1 null. The NEW both-pass verdicts above are statements about the mature signature under the frozen protocol, not proof of cell-intrinsic aging: the original report's central caveat still applies, because bulk expression changes conflate within-cell aging with shifts in cell-type proportions, and the Track 2 Whole Blood and Testis both-pass effects are small (+0.0119 to +0.0122 in blood and -0.0144 to -0.0132 in testis, LFC per decade bin).

## 6. Files

  Main table (Batch A schema, 165 rows): results/followup-mature-organ-transports-c-gtex.csv
  BH check (independent statsmodels recomputation): results/followup-mature-organ-transports-c-gtex-bh-check.csv
  Runtime (per-step seconds, peak RSS, donor counts): results/followup-mature-organ-transports-c-gtex-runtime.csv
  Track 1 gene-level OLS: results/followup-mature-organ-transports-c-gtex-track1-genes.csv
  Track 1 per-tissue summary: results/followup-mature-organ-transports-c-gtex-track1-genes-summary.csv
  Track 1 P1m/P0m transport: results/followup-mature-organ-transports-c-gtex-track1-p1.csv
  Track 2 NB-GLM gene-level: results/followup-mature-organ-transports-c-gtex-track2-nbglm-genes.csv
  Track 2 per-tissue summary: results/followup-mature-organ-transports-c-gtex-track2-nbglm-summary.csv
  Track 2 P1m/P0m transport: results/followup-mature-organ-transports-c-gtex-track2-p1.csv
  Script: scripts/run_followup_mature_organ_transports_c_gtex.py

