# Item 5 Batch A: COMPLETE (all four cohorts)

Status: COMPLETE. All four cohorts ran at 30-plus with mature P1m and P0m:
Gao heart, Yazar blood, Perez blood controls, Guilliams liver sorted cells.
Summary CSV: results/followup-mature-organ-transports-a-summary.csv (30 rows).
BH-check CSV: results/followup-mature-organ-transports-a-bh-check.csv (866 rows).
Every headline number below was read from a result CSV or the summary CSV and
independently recomputed; nothing is taken on trust from an agent log alone.

## What happened, in order

- Antigravity Opus 5.5 High: launch at 21:16 PST Nov 21 died immediately with
  429, quota resets in about 140 hours.
- Switched to Gemini 4 Neon XHigh per the task fallback rule. Batch A ran 21:18
  to about 21:49 PST. Gao completed and committed as 9a379c1. Yazar output was
  written but left uncommitted and unverified. Gemini died with 429 while
  starting Cohort 3 (Perez). Error text: "Individual quota reached. Please
  upgrade your subscription to increase your limits. Resets in 4h26m46s."
  Error id 2cc780e9-e6e7-4099-b989-23e96a443008-87. Kimi remains quota-dead
  since about 13:13 PST. No model calls were retried after the stop order.
- The takeover agent then found a q-value error in the Yazar file (see the fix
  section below) and stopped on an automated execution-risk flag attached to
  the Batch A script. Reza explicitly authorized running the script in chat on
  Nov 21 ("Yes, run the Batch A script"). The parent agent had read the script
  first: it downloads two public h5ad files via curl and computes statistics.
- The remaining work was then executed directly on beast, no coding-agent CLI:
  the script was fixed, Yazar was rerun, Perez and Guilliams were run and then
  rerun with the fixed script (downloads reused from the verified local files),
  and finalize was run. Repo HEAD at the fix was db1b836.

## The q-value fix (applies to Yazar, Perez, Guilliams; Gao unaffected)

- The bug: in the transport assembly the script computed BH q-values over the
  cell types only and assigned the pooled group (All_Cells) its raw p-value,
  uncorrected. For Yazar that meant BH over 16 types, with the All_Cells rows
  carrying q equal to raw p. The script already computed the correct all-groups
  BH and stored it in the notes column (q_all17_comp / q_all17_emp for Yazar,
  q_all_comp / q_all_emp for Perez and Guilliams); the mature_q columns just
  did not use it.
- The fix (scripts/run_followup_mature_organ_transports_a.py): the transport
  mature_q and mature_empirical_q columns now come from the all-groups BH for
  every tested unit, pooled group included. The summary both-pass counts now
  use the same all-groups family and count the pooled group in the denominator.
  The old cell-type-only BH entries remain in the BH-check CSV as a labeled
  sensitivity family (ct16 for Yazar). Gao needed no change: its 14-group
  transport family already includes every tested unit (All_LV_Nuclei is one
  of the 14).
- Corrected Yazar numbers (BH over all 17 tested units). The two NEW both-pass
  rows survive the correction, with different q-values:
  All_Cells P1m composite q 4.03e-10 (was 2.37e-11), empirical q 0.00425
  (was 0.0005); CD14-low, CD16-positive monocyte P1m composite q 0.0142
  (was 0.0267), empirical q 0.00425 (was 0.0080). P0m: All_Cells composite q
  8.92e-10, empirical q 0.00425; monocyte composite q 0.0130, empirical q
  0.00425. All four both-pass slopes are negative (wrong direction), see the
  Yazar section.
- Verification, independent of the script: a fresh BH implementation recomputed
  every q in all four cohort CSVs from the stored raw p-values over the
  declared families, re-derived every transport class and verdict, and checked
  every notes q_all value against its column. Result: 0 problems in every
  family (script scripts/verify_batcha_bh.py, run from /tmp on beast).
  Finalize also passed: summary 30 rows, BH-check 866 rows, max BH explicit
  vs statsmodels difference 1.110e-16, SHA-256 of all 36 original input files
  unchanged.
- Two wording corrections to the earlier partial report. First, NULL verdict
  rows are legitimate outcomes (neither original nor mature significant), not
  blank rows: 48 of 68 Gao rows and 26 of 90 Yazar rows are NULL. Second, the
  Gao Cardiomyocyte::STRBP genes_hit row has verdict NULL (original q 0.0525,
  mature q 0.5118), not a blank verdict.
- Genes q definition, confirmed from the cohort assembly code: each
  genes_celltype row's mature_q is the minimum per-gene global BH q within
  that cell type (the q of its smallest-p gene), not a BH across cell types.
  The BH-check CSV carries the top min-p gene per group for that check.

## Gao Heart (Cohort 1): COMPLETE, unchanged from 9a379c1

- Script run by Gemini 4 Neon XHigh, seed 20261105, 2000 draws. Runtime
  116.6 s, peak RSS 4226.93 MB. Source file:
  results/followup-mature-organ-transports-a-gao.csv (68 rows).
- Donors: 109 mature donors at 30-plus, 8 excluded under 30.
- Fractions (CLR, study-adjusted): 0 of 13 at mature q below 0.05. Minimum
  mature q 0.2539 (Lymphoid, p 0.01953). Original null held.
- Genes: 0 of 145,778 ct13 tests and 0 of 158,960 all14 tests significant at
  mature q below 0.05. Of the 12 listed original hit genes, 11 DIED and 1
  (Cardiomyocyte::STRBP) is NULL. Of the 5 cell-type rows with original q
  below 0.05, all 5 DIED.
- Concordance: null, mean r 0.0984, perm p 0.5095, 7428 common genes.
- Transport: 0 of 14 both-pass under P1m (10 empirical-only, 0
  composite-only), 0 of 14 under P0m (9 empirical-only, 0 composite-only).
  Original endocardial both-pass DIED: P1m slope -0.0095 per decade,
  q_comp 0.4526, q_emp 0.6257; P0m slope -0.0059, q_comp 0.7015,
  q_emp 0.9200.
- Verdicts: 48 NULL, 20 DIED, 0 SURVIVED, 0 NEW.
- Count reconciliation (closed): the per-row n_genes_tested sums are 145,778
  over the 13 cell types and 158,960 over all 14 groups, matching the runtime
  CSV exactly (ct13_tests 145,778). The agent log numbers (145,756 and
  158,938) were wrong by 22 in both cases. All sources agree 0 are
  significant.
- Independent BH recheck: all 6 Gao families (transport P1m/P0m composite and
  empirical, fractions, concordance) recompute exactly, 0 mismatches.

## Yazar Blood (Cohort 2): COMPLETE, fixed and rerun

- Script rerun directly on beast with the fixed script, seed 20261105, 2000
  draws. Runtime 686.7 s, peak RSS 5041.33 MB. Source file:
  results/followup-mature-organ-transports-a-yazar.csv (90 rows).
- Donors: 932 mature donors at 30-plus, 49 excluded under 30.
- Fractions (CLR, pool fixed effects): 20 of 28 at mature q below 0.05
  (original 21 of 28); unadjusted 18 of 28. Top hit naive thymus-derived
  CD8-positive, alpha-beta T cell, mature q 7.773e-72.
- Genes: 908 of 134,133 pool-FE cell-type tests significant at mature q below
  0.05 (original 1741 of 134,036); All_Cells 5445 of 13,112 (original 6241 of
  13,103). Unadjusted 7887 of 160,346.
- Concordance: pool-FE mean r 0.0794, perm p 0.0010 (original r 0.1024,
  p 0.0010), 16 types, 5407 genes; unadjusted r 0.2297.
- Transport: P1m 2 of 17 both-pass (0 positive), 0 empirical-only, 0
  composite-only; P0m 2 of 17 both-pass (0 positive), 0 empirical-only, 0
  composite-only. All four both-pass rows are NEW and all four slopes are
  negative: All_Cells P1m effect -0.0151 per decade, q_comp 4.03e-10,
  q_emp 0.00425; CD14-low, CD16-positive monocyte P1m effect -0.0154,
  q_comp 0.0142, q_emp 0.00425; All_Cells P0m effect -0.0139,
  q_comp 8.92e-10, q_emp 0.00425; monocyte P0m effect -0.0145,
  q_comp 0.0130, q_emp 0.00425.
- Verdicts: 38 SURVIVED, 5 NEW, 12 DIED, 26 NULL, 9 UNTESTABLE.
- Headline: no positive-direction mature transport pass in Yazar blood. The
  only mature both-pass units score lower with age, opposite to the
  signature's predicted direction.

## Perez Blood Controls (Cohort 3): COMPLETE

- Script run directly on beast with the fixed script, seed 20261105, 2000
  draws. The 12,218,105,530 byte h5ad was downloaded once from cellxgene
  (verified exact size) and reused for the fixed rerun. Fixed rerun runtime
  108.5 s, peak RSS 643.44 MB. Source file:
  results/followup-mature-organ-transports-a-perez.csv (126 rows).
- Donors: 61 mature donors at 30-plus, 38 excluded under 30.
- Fractions (CLR): 0 of 11 at mature q below 0.05 in all three arms
  (cohort multihot minimum q 0.2111, top non-classical monocyte p 0.02402;
  cohort level ge2 minimum q 0.2229; unadjusted minimum q 0.1322).
  Original unadjusted arm had 2 of 11.
- Genes: cohort multihot 3 of 93,345 cell-type tests (original 391 of 94,539),
  All_Cells 0 of 12,293 (original 46); cohort level ge2 3 of 93,648
  (original 295 of 94,837), All_Cells 0 of 12,293 (original 36); unadjusted
  2186 of 93,345 (original 10021 of 94,539), All_Cells 2554 of 12,293
  (original 4785).
- Concordance: cohort multihot mean r 0.1872, perm p 0.0120 (original r
  0.2011, p 0.0280), 9 types, 7381 genes; cohort level ge2 r 0.1695,
  p 0.0220 (original r 0.1775, p 0.0490); unadjusted r 0.4005, p 0.0010
  (original r 0.4519, p 0.0010).
- Transport: P1m cohort multihot 6 of 10 both-pass (0 positive), 2
  empirical-only; cohort level ge2 2 of 10 both-pass (0 positive), 6
  empirical-only; unadjusted 9 of 10 both-pass (0 positive), 1
  empirical-only. P0m cohort multihot 6 of 10 both-pass (0 positive), 3
  empirical-only; cohort level ge2 2 of 10 both-pass (0 positive), 7
  empirical-only; unadjusted 9 of 10 both-pass (0 positive). Every mature
  both-pass slope in every arm is negative. All_Cells P1m cohort multihot
  effect -0.0231, p_comp 0.04008, q_comp 0.1052, p_emp 0.0004998,
  class EMPIRICAL_ONLY.
- Verdicts: 57 NULL, 16 DIED, 20 SURVIVED, 33 NEW.
- Headline: same pattern as Yazar. Mature both-pass units exist in Perez
  blood controls, but all of them score lower with age, opposite to the
  signature's predicted direction. No positive-direction mature transport
  pass.

## Guilliams Liver Sorted Cells (Cohort 4): COMPLETE

- Script run directly on beast with the fixed script, seed 20261105, 2000
  draws. The 1,026,110,286 byte h5ad was downloaded once from cellxgene
  (verified exact size) and reused for the fixed rerun. Fixed rerun runtime
  48.4 s, peak RSS 384.05 MB. Source file:
  results/followup-mature-organ-transports-a-guilliams.csv (137 rows).
- Donors: 9 mature donors at 30-plus, 1 excluded under 30. Power is minimal;
  treat every null here as low power, not as evidence of absence.
- Fractions: UNTESTABLE. The single fractions row (All_Sorted_Cells) is
  untestable because FACS-sorted CD45-positive and CD45-negative enrichment
  distorts whole-liver proportions; not run in the original either.
- Concordance: UNTESTABLE (FACS-sorted exclusion retained verbatim; no
  cross-cell-type concordance in the original).
- Genes: 1 of 69,749 cell-type tests significant at mature q below 0.05
  (original 0 of 81,881), a monocyte hit, top gene ZNF570, q 0.0392;
  pooled 0 of 25,942 (original 0 of 26,220). 6 of 7 types tested, 1
  untestable.
- Transport: 0 of 8 both-pass under P1m in all five arms (primary 5
  empirical-only; adjust_citeseq 5; female_only 4 of 7 tested; no_any_tumor
  2 of 6 tested; no_malignant_resection 3 of 7 tested), 0 of 8 under P0m
  in all five arms (primary 4 empirical-only). Pooled groups: primary P1m
  All_CD45pos_Cells effect -0.0292, p_comp 0.4928, p_emp 0.0004998;
  All_Cells effect -0.0150, p_comp 0.7122, p_emp 0.02399.
- Verdicts: 105 NULL, 29 UNTESTABLE, 1 NEW (monocyte genes), 2 DIED.
- Headline: no mature transport pass of either sign in Guilliams liver at
  30-plus, on 9 donors.

## Batch A headline

No positive-direction mature transport pass exists in any Batch A cohort.
Gao 0 of 14 both-pass under either signature. Yazar 2 of 17 both-pass under
each signature, all four negative. Perez 6 of 10 (multihot), 2 of 10 (level
ge2), 9 of 10 (unadjusted) both-pass under each signature, all negative.
Guilliams 0 of 8 both-pass in every arm under either signature. The mature
signature does not transport in the predicted direction in heart, blood, or
liver at 30-plus. The only mature both-pass units anywhere in Batch A are
blood units (Yazar All_Cells and monocyte, Perez blood types) scoring lower
with age.

## Files

- scripts/run_followup_mature_organ_transports_a.py (fixed: all-groups BH
  for transport q columns in Yazar, Perez, Guilliams; summary counts use the
  same family)
- scripts/verify_batcha_bh.py (independent BH verification used for this
  report; stdlib only, nothing imported from the Batch A script)
- results/followup-mature-organ-transports-a-gao.csv (68 rows)
- results/followup-mature-organ-transports-a-yazar.csv (90 rows)
- results/followup-mature-organ-transports-a-perez.csv (126 rows)
- results/followup-mature-organ-transports-a-guilliams.csv (137 rows)
- results/followup-mature-organ-transports-a-runtime.csv
- results/followup-mature-organ-transports-a-summary.csv (30 rows)
- results/followup-mature-organ-transports-a-bh-check.csv (866 rows)
- reports/followup-mature-organ-transports-batch-a.md (this file)

## Deviations and open notes

- The runtime CSV model_used column still reads "Gemini 4 Neon (XHigh)" on
  the rerun rows. That value is hardcoded in the script from the original
  run; the Yazar rerun and the Perez and Guilliams fixed reruns were executed
  directly on beast by the takeover agent, not by Gemini. Do not cite the
  runtime label as the runner of the reruns.
- The first Perez and Guilliams runs used the unfixed script and were
  superseded by the fixed reruns; their outputs were overwritten, not kept.
- The runtime CSV holds one cohort_total row per cohort, from the latest run
  of each (Yazar rerun 686.7 s, Perez fixed rerun 108.5 s, Guilliams fixed
  rerun 48.4 s, Gao original 116.6 s).
- MDE80 values are present in the CSVs but are not quoted in this report;
  they are available per row in the mde80 column.
- Nothing under report/ was touched. No manuscript file was touched. No
  submission.
