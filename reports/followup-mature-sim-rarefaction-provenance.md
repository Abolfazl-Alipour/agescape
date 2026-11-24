# Sim leave-out refit and rarefaction split counts: provenance trace

Date: 2026-11-24. Base: a013ebe. Rule: every tex number must trace to a script plus a file.
Tex: report/plos-one/final-plos-one.tex lines 290, 294, 296. Tex not edited.

## 1. Leave-Sim-out refit: Gao minus 0.0735 per decade, p 0.038

1. What "Sim" is: the Sim 2021 study, one of the 7 Gao heart studies, 3 adult donors: GSM4742854_Adult1 (age 35, female), GSM4742855_Adult2 (age 42, male), GSM4742856_Adult3 (age 42, male). Source: results/followup-organ-17-heart-gao-donors.csv lines 108 to 110.
2. No source found in the repo: the exact value 0.0735 occurs only in the tex (lines 290, 296). No script, report, note, or result file contains it or any leave-Sim-out refit. Searched scripts/, results/, reports/, notes/, docs/, src/, tests/ for 0.0735, leave-Sim, leave-one-study-out, and sim_2021. Git history: the sentence entered the tex in commit da532ab, which touched only the tex file. The prior wording said explicitly "no leave-these-donors-out refit was run, so their influence is unknown."
3. Scratch search: no /home/reza/scratch or /home/reza/tmp exists. Exact 0.0735 grep over the data root and /tmp found only: GTEx sample-attributes metadata values (unrelated file, no Agescape analysis), my own trace prompt and recompute scripts in /tmp, and two /tmp copies of the tex itself. A "leave-Sim-out" grep over /tmp found only those same files. No refit script or output file exists there.
4. Recompute from committed data (new check, this trace): the primary Gao scatter fit was replicated exactly from results/followup-mature-fulldepth-gao-donor-part1.csv using the pipeline spec in scripts/run_followup_mature_fulldepth.py lines 798 to 969 (fit_stratified_global_with_batch: per-cell-type standardization, sex and depth centered within cell type, cell-type dummies, study dummies with ENCODE v4 (Snyder) as reference, donor-clustered standard errors). Replication: slope minus 0.0722825, standard error 0.0349365, p 0.0410284, 105 donors, 829 rows, 480,355 cells. That matches results/followup-mature-fulldepth-gao-cell-drift.csv row 599 to full precision.
5. Same spec with the 3 Sim donors dropped (26 donor-by-cell-type rows): slope minus 0.0935, p 0.0300, 102 donors, 803 rows, 472,556 cells. That does not match the tex (minus 0.0735, p 0.038).
6. Three other defensible specs also do not match: pooled per-donor OLS leave-Sim-out, standardized slope about minus 0.104, p 0.014; stratified fit with depth-by-cell-type interactions, minus 0.0866, p 0.0468; raw unadjusted OLS on the plotted standardized points, minus 0.0516, p 0.0873 (the all-donor raw fit is plus 0.0637, which reproduces the review's raw check, so the points construction is right).
7. Verdict: NO SOURCE FOUND, and the tex numbers do not reproduce under the primary spec or the three alternatives checked. The claim needs a filed refit or removal; that call is not made here.

## 2. Rarefaction split counts: 89 of 96 at 1,500 molecules, 61 of 63 at 3,000

1. Source found: results/followup-mature-fulldepth-robustness.csv. It has 159 rows. Depth rarefied_1500: 96 rows, 89 with passes_robustness true. Depth rarefied_3000: 63 rows, 61 true. Total 150 of 159.
2. Writer and definition: scripts/run_followup_mature_fulldepth.py lines 1300 to 1325. For each row the pipeline refits across seeds (Shannon 5 seeds, pass if at least 4 share the seed-0 sign; KL methods and KL sqrt 20 seeds, pass if at least 18 share the sign) and sets passes_robustness to n_same >= min_ok. The combined 150 of 159 is computed by scripts/report_followup_mature_fulldepth.py line 342 and stated in reports/followup-mature-fulldepth.md line 167.
3. Wording note: "pass" means sign-stable across rarefaction seeds, not literally unchanged values. All 159 rows are Part 2 rows (Shannon, KL methods, KL sqrt: 32 per metric at 1,500 molecules, 21 per metric at 3,000), while the tex sentence sits in the off-program paragraph. The counts are exact; the tex framing is loose.
4. Verdict: SOURCE FOUND. Counts match the file exactly.
