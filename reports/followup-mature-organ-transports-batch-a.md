# Item 5 Batch A: PARTIAL result (Gao complete, Yazar present but unverified, Perez and Guilliams not started)

Status: PARTIAL. This is not the Batch A partial report the saved prompt describes: no BH-check CSV, no summary CSV, and no per-cohort report section were produced by the agent before it died. Every number below was read back from a CSV file by the task runner; nothing here is taken on trust from the agent log alone.

## What happened

- Antigravity Opus 5.5 High: launch at 21:16 PST died immediately with 429, quota resets in about 140 hours.
- Switched to Gemini 4 Neon XHigh per the task fallback rule. Batch A ran 21:18 to about 21:49 PST.
- Gemini died with 429 while starting Cohort 3 (Perez). Error text: "Individual quota reached. Please upgrade your subscription to increase your limits. Resets in 4h26m46s." Error id 2cc780e9-e6e7-4099-b989-23e96a443008-87. Kimi remains quota-dead since about 13:13 PST. No model calls were retried after the stop order.
- Repo HEAD was f91c455 at launch (the saved prompt said 4e43b6a, stale). HEAD was c693534 at commit time; other runners landed commits in parallel.

## Gao Heart (Cohort 1): COMPLETE

- Script: scripts/run_followup_mature_organ_transports_a.py (Gemini 4 Neon XHigh, seed 20261105, 2000 draws).
- Donors: 109 mature donors at 30-plus, 8 excluded under 30. Runtime 116.6 s, peak RSS 4226.93 MB. Source file: results/followup-mature-organ-transports-a-gao.csv (68 rows).
- Fractions (CLR, study-adjusted): 0 of 13 at mature q below 0.05. Original null held.
- Genes: 158,960 tests summed across 14 cell-type rows, 0 significant at mature q below 0.05. Of the 12 listed original hit genes, 11 DIED and 1 (Cardiomyocyte::STRBP, original q 0.0525) has a blank verdict. Of the 5 cell-type rows with original q below 0.05, all 5 DIED.
- Concordance: null, mean r 0.0984, p 0.5095.
- Transport: 0 of 14 both-pass under P1m, 0 of 14 under P0m. Original endocardial both-pass DIED: P1m slope -0.0095 per decade, q_comp 0.4526, q_emp 0.6257; P0m slope -0.0059, q_comp 0.7015, q_emp 0.9200.
- Count discrepancy to flag for the coordinator: the agent log says 145,756 ct13 tests and 158,938 all14 tests; the runtime CSV says ct13_tests 145,778; the per-row CSV sum is 158,960. All three agree that 0 are significant; reconcile the denominators during verification.
- 48 of 68 Gao rows have a blank verdict. Blank rows are those with no original q below 0.05 to compare against; the verdict rule needs a coordinator decision before the merge.

## Yazar Blood (Cohort 2): OUTPUT PRESENT, NOT VERIFIED

- Source file (uncommitted, on disk): results/followup-mature-organ-transports-a-yazar.csv (90 rows), written 21:49 PST. Runtime record says cohort_total ok, 712.7 s, 932 mature donors, 49 excluded under 30, peak RSS 5829.86 MB.
- File contents as read (NOT verified against a BH check, which was never produced): fractions_clr 28 rows, among rows with a verdict 19 SURVIVED, 2 DIED, 1 NEW; genes_celltype 20 rows, 17 SURVIVED, 3 UNTESTABLE; transport_P1m and transport_P0m 20 rows each, 5 DIED, 2 NEW, 3 UNTESTABLE each; concordance rows SURVIVED in both arms.
- Transport NEW both-pass rows in the file: All_Cells P1m q 2.37e-11, empirical q 0.0005; CD14-low, CD16-positive monocyte P1m q 0.0267, empirical q 0.0080. The same two rows are NEW under P0m (All_Cells q 5.25e-11, empirical q 0.0005; monocyte q 0.0244, empirical q 0.0080). These are agent-computed values only; the coordinator must verify them before any reporting.
- Quirk: 26 of 90 Yazar rows have a blank verdict (mostly FAIL_BOTH transport rows with no original pass). Same blank-verdict issue as Gao.
- Yazar is left UNCOMMITTED deliberately pending that verification.

## Perez Blood Controls (Cohort 3) and Guilliams Liver (Cohort 4): NOT STARTED

- No files exist. No numbers exist. Resume after the Gemini quota reset (about 18:30 PST Nov 21) or on another approved model.

## Deviations and unverified items

- No BH-check CSV and no summary CSV were produced; every q value above is quoted from the per-cohort CSVs, not from an independent BH recomputation.
- The saved prompt said "do not commit"; the parent task overrode that and ordered this commit.
- The saved prompt said HEAD was 4e43b6a; actual HEAD at launch was f91c455.
- MDE80 values are present in the CSVs but are not quoted here; the coordinator should check them during verification.
- Nothing under report/ was touched. No manuscript file was touched. No submission.
