# Review round 4: statistician (hostile)

Manuscript reviewed: `report/plos-one/final-plos-one.tex` at commit `dd6f503`, tag `mature-to-old-donor-scatterplots-2026-11-22`. Four main figures in `report/plos-one/figures-v2/`, six SI figures in `report/plos-one/figures-si/`. Project Q and A checked against the same commit. This is the third hostile review in round 4. Sibling round 4 reviews: aging biology (MAJOR REVISION) and top tier (REJECT at top tier standard). This review is written from the statistical record only: committed CSVs, committed reports, and the figure builder scripts. No tex or figure was edited.

Method: every number in the Abstract, Results, Discussion, and Conclusions was checked against its source CSV or report. Per-donor refits were rerun independently for Figure 1 and Figure 2. All multiple-testing families were recomputed from the committed BH files. Tex line numbers are from `grep -n` on the `dd6f503` tex (656 lines).

## Verdict: MAJOR REVISION

Counts: 2 BLOCKER, 14 MAJOR, 9 MINOR.

## What verified clean

1. Figure 1: independent OLS refits of all five headline CLR slopes from `results/followup-mature-composition-donor-clr.csv` reproduce the committed values to machine precision (club minus 0.080407, SE 0.019506, p 8.476e-05, n 91; multiciliated minus 0.047002, SE 0.012896, p 0.000452, n 91; L6b plus 0.019765, SE 0.004882, p 7.298e-05, n 209; oligodendrocyte plus 0.017269, SE 0.004377, p 0.000109, n 209; OPC minus 0.011896, SE 0.002989, p 9.545e-05, n 209). The builder asserts refit equals committed to 1e-9 and plots partial residuals (sex removed).
2. Figure 3: counts match `results/followup-mature-p1m-block5-summary.csv` (P1m 2/13/4/29, P0m 3/15/5/25); builder asserts.
3. Figure 4: ring counts match the actual source `results/followup-mature-interneuron-genes.csv`, tier both_FDR (VIP 4, pvalb 8, sst 41, GABA 0). TRMT11 slopes verified (0.00205 to 0.00406 log counts per 10,000 per year; 2 to 4 percent per decade). The phase3 figure-data CSV is a stale file and must not be used.
4. The Abstract and Results numbers for drift, CD8 scatter, CD4 dispersion, block 10, interneuron, between-subject, and concordance p-values all verify against their source CSVs; they are cited as checked in the items below.
5. The 624 family counts verify: 624 rows (431 core plus 193 fly/worm/zebrafish), 45 passes (17 core plus 28 added), 15 under the 431 family, 30 pass only under 624 (`results/rev3-composition-global-bh-624.csv`).
6. The tex discloses the 624 selection timing, the P1m sex confound (all four 30-month mice male), the 14-fold P1m size difference, and the MSSM depth caveats honestly. The problems below are not a failure to disclose; they are wrong numbers, wrong families, and omitted sensitivities.

## BLOCKER items

### BLOCKER 1. The concordance "primary three-test family" does not exist in the analysis record, and the tex contradicts its own source report on a pass/fail verdict.

Tex lines 291 (Results), 213 (Methods), 337 (Discussion), 403 (Conclusions), 441 (S5 caption).

1. Tex line 291: mouse within-tissue is "nominal only" with "q = 0.066 under the primary three-test family of mouse within-tissue, mouse cross-tissue, and human tests". It dismisses the two-test family ("excludes the cross-tissue test without justification", q 0.044) and the six-row family ("counts near-duplicate sensitivity gene sets", q 0.046) as sensitivity families.
2. Source: `results/followup-mature-concordance-bh-check.csv` computes exactly four families: `original_six_species_primary_within`, `mature_mouse_six_rows`, `mature_two_test_family`, `original_mouse_six_rows`. NO three-test family exists. Grep for "three-test" and "0.066" in `reports/followup-mature-concordance.md` and all concordance CSVs returns nothing.
3. The source report's verdicts section (`reports/followup-mature-concordance.md` lines 14 to 19 and 212 to 215) says the mouse within-tissue result "passes 0.05 under both mature BH families" (two-test q 0.0440, six-row q 0.0460) and that "the nominal signal is now significant under both mature BH families, where the original was nominal-only".
4. The q 0.066 is the correct BH arithmetic for a three-test family over p 0.0220, 0.1029, 0.8152 (0.0220 times 3 divided by 1 equals 0.066). The arithmetic is right; the family was never computed, never prespecified, and appears nowhere in the record.
5. Adjudication: the correct primary family is the mature two-test family (mouse within-tissue primary, human), q 0.0440, PASSES 0.05. It is the faithful mature analogue of the prespecified original six-species primary family: one test per applicable species, within-tissue scope, primary gene set only. The original family had six species; the mature run has two applicable species because fly, worm, zebrafish, and NMR are "not applicable as mature-to-old tests under the mouse and human definitions" (tex line 291 itself). The cross-tissue test is a different scope and was a separate question in the original design; adding it post hoc creates a new family, not the primary one. The six-row family is invalid as primary because its three gene sets are near-duplicates of one hypothesis (p 0.0220, 0.0230, 0.0230; q 0.046 at BH rank 3 is triple-counting).
6. Consequence: "nominal only" is the wrong verdict in Results, Discussion, Conclusions, and the S5 caption. The correct verdict is that mouse within-tissue concordance passes the primary two-test family (q 0.0440), with the power caveat intact (excess over null 0.064 is smaller than MDE80 0.068; the signal is carried by limb muscle alone; the null mean is 70 percent of the observed value).
7. Fix: adopt the two-test family as primary and rewrite lines 291, 213, 337, 403, and 441, OR produce prespecification evidence for the three-test family (none exists in the record). The Q and A line 97 repeats the three-test claim and needs the same fix. Note: the aging-biology round 4 review accepted the three-test family as primary; it did not check whether the family exists in the analysis record.

### BLOCKER 2. PBMC reports one implementation's numbers while the paper's own batch C report says another; the two implementations use different unit sets and different BH families.

Tex lines 281 (Results), 507 (Appendix B table). Q and A line 55.

1. Tex line 281: "0 of 35 cell types pass both tests, with the smallest composite q 0.6299 in Zhu conventional dendritic cells"; "26 of 35 rows at the 1/201 floor"; "the 29 empirical-only passes are not claimed".
2. These three numbers verify against `results/followup-mature-pbmc-p1m-models.csv` (35 primary rows: 22 Zhu plus 13 Thomson; 26 at random_pvalue 1/201; 29 with random_q_bh below 0.05; smallest q_bh 0.6299, tied for Zhu CD16-MC, cDC, pDC).
3. But `results/followup-mature-organ-transports-c-pbmc.csv` has 36 units per tier (23 Zhu plus 13 Thomson), and `reports/followup-mature-organ-transports-batch-c.md` says "Both-pass 0 of 36 under P1m and 0 of 36 under P0m" with "Per-cohort BH families (23 Zhu, 13 Thomson) per tier". The transport CSV gives different numbers: P1m smallest composite q 0.6586 (CD16-MC, cDC, pDC tied), P0m smallest 0.6019; P1m empirical-only 32, P0m 30 (versus the tex's 29).
4. The missing unit in the p1m-models primary is Zhu CD4-Treg (present in the transport CSV, absent from the p1m-models primary). The two implementations also use different BH families (pooled 35 versus per-cohort 23/13), unreconciled anywhere in the tex.
5. Fix: reconcile the two implementations in one place, state which unit set and which BH family the manuscript's numbers come from, explain why Zhu CD4-Treg is absent from the primary, and correct lines 281 and 507 (and Q and A line 55) to the reconciled numbers. The headline "0 of 35" versus "0 of 36" is a factual inconsistency between the manuscript and its own analysis record.

## MAJOR items

### MAJOR 1. The 624-test composition family was selected after both correction outcomes were available, which invalidates FDR control; disclosure does not cure the design.

Tex lines 33 (Abstract), 225 (Results), 321 (Discussion), 403 (Conclusions).

1. Tex line 225: "The 624-test family was selected after both correction outcomes were available, and its experiment-wide rationale does not remove that selection timing."
2. Verified counts (`results/rev3-composition-global-bh-624.csv`): 624 rows (431 core plus 193 fly/worm/zebrafish); 45 passes under the 624 family (17 core plus 28 added); 15 passes under the 431 family; 30 rows pass only under 624 (2 core plus 28 fly/worm/zebrafish). The two promoted core tests are lung vein endothelial (p 0.00272, core q 0.0690, 624 q 0.0390, 320 cells, 14 donors) and kidney cortex collecting duct (p 0.00248, core q 0.0667, 624 q 0.0377, 45 cells, 11 donors). Both fail every sensitivity check (total-cell-adjusted q 0.0627 and 0.4231; propeller p 0.0528 and 0.0588).
3. Choosing the correction family after seeing both correction outcomes is outcome-dependent selection. FDR control is a property of a prespecified family; a family chosen because it promotes two tests that fail the core family does not control the false discovery rate at the stated level. The tex discloses the timing honestly, but disclosure does not make the 45-pass count a valid FDR-controlled result. The paper's headline composition number (45 passes, 17 core plus 28 added) rests on this family.
4. Fix: either demote the 624 family to a sensitivity analysis and report the 431-family 15 passes as the primary composition result, or keep 624 primary with a stated, prespecified rationale that was fixed before the correction outcomes were known (none exists in the record). The keep-624 decision was made by the author on 2026-11-20 and is not re-litigated here; the statistical invalidity stands as a review finding. Round 3 raised this as R3-M3 and it was rejected on the author's decision; this review records the statistical consequence, not a new decision.

### MAJOR 2. Figure 2 fitted lines contradict the plotted donor clouds; the two headline full-depth findings reverse sign without the covariate adjustment.

Tex lines 239 (Results), 241 (Fig 2 caption), 321 (Discussion), 403 (Conclusions).

1. Committed stratified slopes verify against the table CSVs: Gao scatter minus 0.072283, p 0.0410; Yazar scatter minus 0.019665, p 6.69e-05; pooled center shift plus 0.040003, p 2.99e-05; stratified center shift minus 0.004388, p 0.2659; Gao Shannon minus 0.100142, p 0.0539; Yazar Shannon minus 0.029773, p 0.00365.
2. But independent raw OLS on the plotted standardized donor points gives the OPPOSITE SIGN in two of five panels: Gao scatter raw plus 0.0637, p 0.030 (committed minus 0.0723); Yazar pooled center shift raw minus 0.0855, p 0.00023 (committed plus 0.0400). Yazar scatter raw minus 0.1434 (7 times the committed minus 0.0197); Gao Shannon raw minus 0.1952 (2 times the committed minus 0.1001); Yazar Shannon raw minus 0.1258 (4 times the committed minus 0.0298).
3. The builder (`scripts/build_optimization_figures.py`, fig2) plots RAW standardized points (not partial residuals) with the adjusted line. The Figure 1 builder plots partial residuals (sex removed), so Figure 1 does not have this defect.
4. The sign flips are driven by the adjustment: Gao stratified scatter covariates are sex_c, depth_c, cell-type stratum dummies, and batch fixed effects (6); Yazar pooled center shift covariates are sex_male, depth, and about 70 batch indicators.
5. The Gao low-scatter vertical stack at ages 35 to 42 is three Sim 2021 donors (GSM4742854_Adult1 age 35, mean scatter 0.532; GSM4742855_Adult2 age 42, mean 0.520; GSM4742856_Adult3 age 42, mean 0.530), the three lowest mean-scatter donors in the Gao cohort: a study effect, not an age effect.
6. The caption (line 241) discloses that the line is the stratified-model slope, so no number is wrong. But the figure is misleading: a reader who reads the cloud sees a rise where the paper claims a fall (Gao scatter) and a fall where the paper claims a rise (pooled center shift). The two headline full-depth findings are entirely adjustment-driven.
7. Fix: plot partial residuals (as Figure 1 does), or plot the unadjusted slope as well, or add a caption note stating the direction of the raw cloud and that the finding reverses without the batch and depth adjustment. Identify the Sim 2021 stack in the caption or remove the three donors in a sensitivity.

### MAJOR 3. The lung study-aware fixed-effects slopes are omitted; the tex reports only the meta-combination numbers and the headline "five hold" carries no qualification.

Tex lines 231 (Results), 33 (Abstract), 321 (Discussion), 403 (Conclusions).

1. `reports/followup-mature-composition.md` lines 110 to 118: club study fixed-effects slope minus 0.0426 (p 0.0155, n 91) versus meta-combination fixed effect minus 0.0061 (p 0.58, n 69); multiciliated study fixed effects minus 0.0261 (p 0.0518) versus meta plus 0.0098 (p 0.39).
2. Tex line 231 reports ONLY the meta-combination numbers ("The two lung survivors do not hold within individual studies") and omits the study fixed-effects slopes, which still show a halved but significant club effect.
3. Per-study club slopes (7 tested studies, none p below 0.05, mixed signs): AN minus 0.0097, p 0.83, n 5; Donor plus 0.0088, p 0.70, n 5; Patient minus 0.0017, p 0.99, n 6; T minus 0.0355, p 0.17, n 11; THD plus 0.0305, p 0.29, n 13; VUHD minus 0.0398, p 0.17, n 19; muc plus 0.0226, p 0.70, n 10.
4. The two lung rows are the only lung findings among the five headline fraction shifts, so their study-confounding bears directly on "Five human cell-type fraction shifts hold from maturity onward" in the Abstract (line 33) and Conclusions (line 403), which carry no qualification.
5. The fixed-effects result is itself fragile (the common-slope assumption is violated by the heterogeneous per-study slopes), so the honest reporting is both analyses with the heterogeneity noted, not one analysis alone.
6. Fix: report both the meta-combination and the study fixed-effects slopes for club and multiciliated cells, note the per-study heterogeneity, and qualify the "five hold" headline in the Abstract and Conclusions.

### MAJOR 4. The 399-row composition scan uses different multiple-testing families in different places; "clear the bar" has no single definition.

Tex lines 231 (Results), 227 (Fig 1 caption), 435 (S2 caption).

1. Pass counts depend on the family: q_refit below 0.05 gives 10 rows; q_431 below 0.05 gives 8 rows; q_bh (within-file) below 0.05 gives 14 rows (counted from `results/followup-mature-composition-mouse.csv` and `results/followup-mature-composition-human.csv`).
2. The 4 extra q_bh rows: human GABAergic neuron (q_bh 0.0099, q431 0.0568, qrefit 0.0574, p 0.00158); human smooth muscle cell (q_bh 0.0319, q431 0.1939, qrefit 0.1907, p 0.00765); human sst GABAergic cortical interneuron (q_bh 0.0158, q431 0.0973, qrefit 0.0970, p 0.00316); human kidney cortex collecting duct epithelial cell (q_bh 0.0421, q431 0.0821, qrefit 0.0824, p 0.00248).
3. S2 caption (line 435) says "The other 389 rows do not clear the bar" (10 pass), but under the within-file family 14 pass.
4. Tex line 231, "Five further cell types clear the bar in the same scan": the two q 0.051 rows (SEA-AD L5 extratelencephalic neuron and lung multiciliated columnar cell; q431 0.0512, qrefit 0.0499) clear only under the refit family, not the stated 431 standard.
5. Fig 1 caption (line 227) says q values use the 431 family, "the most conservative of the three families computed" (8 pass).
6. The five scan rows otherwise verify against `results/followup-mature-composition-secondary-checks.csv` (mouse limb muscle macrophage minus 0.143 per month, p 0.00044, 431 q 0.028, total-cell q 0.019, LODO holds, propeller AGREE q 0.0166; mesenchymal stem cell minus 0.100 per month, p 0.00021, q 0.018, total-cell q 0.019, LODO holds, propeller NOT_CONFIRMED q 0.0576; lung conventional dendritic cell plus 0.023 per year, p 0.00061, q 0.033, fails total-cell refit q 0.142, propeller AGREE q 0.0306).
7. Fix: pick one family as the definition of "the bar", state it once, and use it consistently in Results, the Figure 1 caption, and the S2 caption. Report the other families' counts as sensitivities, not as the headline.

### MAJOR 5. The two propeller p-values in the promoted-test sentence come from different analysis scopes (one full-range, one mature), undisclosed.

Tex line 225.

1. Lung vein endothelial propeller p 0.0528 comes from the FULL-RANGE propeller run (`results/reviewer-phase1-composition-propeller.csv` line 83, p 0.0527593732323824). The mature propeller run (`results/followup-mature-propeller.csv`) gives Mouse Lung TMS vein endothelial p 0.6413, q 0.9967.
2. Kidney cortex collecting duct propeller p 0.0588 IS from the mature run (p 0.0588313889343144).
3. So the two propeller p-values in the same sentence come from different analysis scopes, undisclosed. The composition report (`reports/followup-mature-composition.md`) never mentions "vein endothelial".
4. The total-cell-adjusted q values in the same sentence (0.0627, 0.4231) also come from the full-range file (`results/reviewer-phase1-composition-total-adjusted.csv`; lung vein endothelial q_bh 0.0645/0.0627/0.0571 across pseudocounts, kidney collecting duct 0.4167/0.4231/0.5301; the tex quotes the middle row in each case; all above 0.05, so the fail conclusion holds).
5. Fix: quote the mature propeller p (0.6413) for lung vein endothelial alongside the mature kidney value, or state explicitly that the lung propeller p is from the full-range run. The mature propeller check ran on mature-only inputs, 577 rows (`reports/followup-mature-composition.md` lines 191, 193). Note: the round 3 review also cited propeller p 0.0528 from the full-range file and did not catch the mature value.

### MAJOR 6. Gao Shannon entropy is presented as a fall in Discussion and Conclusions while the primary specification is non-significant, and the significant depth-by-cell-type specification is omitted.

Tex lines 239 (Results), 321 (Discussion), 403 (Conclusions).

1. Primary spec (`results/followup-mature-fulldepth-gao-spread.csv`, stratified_global_primary): Gao Shannon minus 0.100142, p 0.0539; Yazar minus 0.029773, p 0.00365.
2. Depth-by-cell-type spec (stratified_global_depth_x_ct, 105 donors): Gao minus 0.101305, p 0.0344; Yazar minus 0.027407, p 0.00701. The Gao depth-by-cell-type spec is significant and is nowhere in the tex.
3. Discussion line 321 and Conclusions line 403 say Shannon entropy "falls in heart in the same direction (minus 0.100, p 0.054)", using a non-significant primary result as corroboration while omitting the sensitivity that would make it significant.
4. Fix: report the range across defensible specifications (primary, depth-by-cell-type, spline-residual, rarefied) and stop using the Gao Shannon result as corroboration for the narrowing claim. The aging-biology round 4 review reached the same conclusion (its MAJOR 1); this review adds the verification that the depth-by-cell-type spec is significant and omitted. Q and A line 111 repeats the omission.

### MAJOR 7. Figure S2 plots mouse per-month and human per-year slopes on one shared axis, and the highlight numbering conflicts with the rank axis.

Tex line 435 (S2 caption); figure `report/plos-one/figures-si/s2-399-row-scan.pdf`.

1. Mouse slopes are per month and human slopes are per year (`scripts/run_mature_composition_complete.py` line 411), but both are plotted on one shared x-axis. Mouse limb muscle macrophage minus 0.143 per month equals minus 1.72 per year, 21 times the club slope, but appears only slightly more negative than the club cell.
2. The rank ordering is also computed across the mixed units, so the ranking is dimensionally incoherent.
3. The y-axis says "Rows sorted by slope (1 = most negative)" but the most negative point overall (grey dot near minus 0.158) is not the dot labeled "1" (macrophage at minus 0.143); the highlight numbering and the rank axis use different orderings.
4. Fix: convert all slopes to a common unit (per year) or split the figure by species, and make the highlight numbers match the rank axis. The aging-biology round 4 review covered the mixed-units axis (its MAJOR 3); the rank/numbering conflict is new here.

### MAJOR 8. No leave-one-donor-out check was ever run on P1m, the signature that underpins every transport claim.

Tex lines 433 (S1 caption), 269 (Results, transport), 273 (Fig 4 caption).

1. P1m is the cleaned mature-to-old signature: 1,914 genes (1,119 up, 795 down), 1,724 human-mapped, derived from kidney, lung, and marrow only, all four 30-month mice male (`reports/followup-mature-p1m-derivation.md`).
2. The only P1m robustness numbers in the record are the multi-unit support counts: 1,212 genes in exactly 1 unit, 702 in 2 or more, 361 in 3 or more, 132 in 5 or more, max 11, 57 with both directions (`results/followup-mature-p1m-genes.csv`; S1 bars sum to 1,914).
3. No leave-one-donor-out check was ever run on P1m. The 310-of-765 leave-one-donor-out sign-flip number in Appendix A (tex line 453) belongs to P1 (young-to-old), NOT P1m. A P1m leave-one-donor-out check would drop one of four old mice, so the derivation is donor-sensitive by construction and the check is the one that matters.
4. Every transport conclusion (170 of 252 testable, only microglia both-pass, unfiltered rescue 0 of 48, organ sweep, PBMC) rests on P1m. A signature whose donor-sensitivity was never tested cannot support a transport null with the strength the paper claims.
5. Fix: run a leave-one-donor-out check on the P1m derivation (drop each of the four old mice in turn, re-derive, and report the overlap with the full P1m), or state explicitly in Methods and Limitations that no leave-one-donor-out check was run on P1m and that the transport conclusions rest on a signature whose donor-robustness is untested.

### MAJOR 9. The Yazar pooled All_Cells gene-fraction signal "passes correction" only because it sits in a family of one.

Tex line 299.

1. Tex line 299: "in Yazar, pooled All_Cells show 0.5487 (q = 0.0030, single-test family)."
2. Source: `reports/followup-mature-block1.md` line 135: pooled Yazar All_Cells frac 0.5487, p 0.002999, q_ct 0.002999; family_size_m is 1. BH over 1 is no correction: q equals p.
3. If the pooled row were in the 16-cell-type family (17 tests), its q would be about 0.002999 times 17, about 0.051, just above 0.05, failing.
4. The tex labels it "single-test family" honestly, but the sentence frames it as one of "Two gene-fraction signals pass correction". A reader cannot tell from "single-test family" that no correction was applied. The Gao myeloid signal (frac 0.6658, p 0.002999, q_ct 0.02699, BH over 9) is a real corrected pass; the Yazar pooled signal is not.
5. Fix: state that the Yazar pooled gene-fraction signal is uncorrected (family of one) and would not survive inclusion in the cell-type family, and do not present it as a corrected pass.

### MAJOR 10. The drift claim "none corrected as a family" is true only under an unstated family definition; under the natural per-depth family, spleen at 3,000 molecules survives.

Tex line 285.

1. Tex line 285: spleen plus 0.453 (p 0.057; p 0.0076 at 1,500 molecules and p 0.0054 at 3,000 molecules), "none corrected as a family".
2. Verified from `results/followup-mature-drift-mouse.csv`: full depth, 7 tested tissues, BH smallest q 0.0993 (none survive); at 1,500 molecules, spleen p 0.00759, BH q 0.0532 (does not survive); at 3,000 molecules, spleen p 0.00536, BH q 0.0375 (SURVIVES over the 7-tissue per-depth family). Over all 28 tissue-by-depth tests, q 0.150 (does not survive).
3. The tex never defines the family. The claim is false under the most natural per-depth family (7 tissues at one depth), where spleen at 3,000 molecules passes at q 0.0375.
4. Fix: define the family explicitly and report the per-depth BH results alongside the pooled 28-test result. Mouse global drift minus 0.013715, p 0.8163, CI minus 0.1387 to plus 0.1113, and Bayesian plus 0.038944, 94 percent HDI minus 0.0998 to plus 0.1773, probability 0.7068, verified; human drift minus 0.007101, p 0.6106, CI minus 0.0345 to plus 0.0203, verified. The interaction numbers in Appendix A (minus 0.0663 unadjusted, p 0.0515; minus 0.0447 depth-adjusted, p 0.0051) are full-range, a different estimand, no contradiction.

### MAJOR 11. The naive CD8 gate contamination claim in Results has no citation.

Tex line 251.

1. Tex line 251: "The naive CD8 gate in old donors is also increasingly contaminated by memory-like, clonally expanded cells, a composition shift rather than per-cell randomness."
2. Grep for Venturi or Egorov in the tex returns nothing. The claim is literature-based, not from project data. The source is the round 2 review (`reports/review-mature-round2-aging-biology.md` line 87), citing Venturi et al., PNAS 2011, and Egorov et al., Frontiers in Immunology 2018.
3. Fix: add the Venturi and Egorov references to the tex and cite them at line 251, or delete the sentence. A substantive mechanistic claim in Results cannot stand uncited. Q and A line 139 repeats the claim uncited.

### MAJOR 12. The lung-versus-brain "average per-unit slope" claim is not reproducible from the committed data and contradicts the per-tissue slopes.

Tex line 341.

1. Tex line 341: "the average per-unit slope is slightly larger in lung than brain".
2. No committed file states this. Per-tissue linear-depth slopes from `results/followup-mature-drift-human.csv`: brain minus 0.00265 (p 0.799, 38 units), lung minus 0.01367 (p 0.893, 34 units), kidney plus 0.11393 (p 0.312, 15 units).
3. Independent per-unit computation from `results/followup-mature-drift-human-donors.csv` (profile Hellinger distance regressed on age plus sex plus mean log UMI per unit): brain mean plus 0.00001 per year, mean absolute 0.00007; lung mean minus 0.00026 per year, mean absolute 0.00066; kidney mean plus 0.00025.
4. Under the signed reading, lung is SMALLER (more negative) than brain, not larger. Under the magnitude reading, lung is 5 to 9 times larger, not "slightly" larger. Unadjusted per-tissue slopes: brain plus 0.0859, p 0.0387 (significant); lung minus 0.0556, p 0.507.
5. Fix: delete the sentence or replace it with the per-tissue slopes and their p-values, stating which estimand (signed mean or mean absolute) is being compared.

### MAJOR 13. Data Availability links point to a stale tag, not the manuscript of record.

Tex lines 411, 413.

1. The workflow manifest link and the tagged-source/archive links point to tag `mature-to-old-rewrite-2026-11-21`, not the current tag `mature-to-old-donor-scatterplots-2026-11-22` (confirmed to exist via `git tag`).
2. The links resolve, but to an older version of the repository, not the manuscript of record. A reader who follows them does not get the code and data that produced this manuscript.
3. Fix: update all links in lines 411 and 413 to `mature-to-old-donor-scatterplots-2026-11-22`. The aging-biology round 4 review flagged this as well.

### MAJOR 14. Figure S1 and its caption mislabel the P1m support units as "tissue-by-species" groups; they are tissue-by-cell-type groups, all mouse.

Tex line 433 (S1 caption); `scripts/build_si_figures.py` line 107.

1. S1 caption: "One test group is one tissue-by-species group."
2. The 36 units are tissue-by-cell-type (`results/followup-mature-p1m-genes.csv` units_list, e.g. "Marrow::hematopoietic precursor cell"), all mouse, one species. A tissue-by-species grouping would have at most 3 groups (kidney, lung, marrow), not 36.
3. The counts themselves verify (1,212 in exactly 1; 702 in 2 or more; 361 in 3 or more; 132 in 5 or more; max 11; 57 both directions; bars sum to 1,914).
4. Fix: change the S1 x-label and the caption to "one test group is one tissue-by-cell-type group". This is a one-word fix in the figure builder and the caption, but it is a factual error in a figure label, not a typo in prose.

## MINOR items

### MINOR 1. Figure 2 panel (b) caption misstates the point definition.

Tex line 241. Caption: "Each point is one donor by one cell type." Panel (b), pooled center shift, is one point per donor (All_Cells, 932 points), not per donor-by-cell-type. Fix: "Each point is one donor by one cell type, except panel (b), which is one point per donor."

### MINOR 2. Figure 4 panels use different axis ranges.

Tex line 273. The builder (`scripts/build_optimization_figures.py`, fig4_interneuron) sets each panel's limits from its own maximum slope (VIP/pvalb about plus/minus 0.008, SST about plus/minus 0.012, broad GABAergic about plus/minus 0.005). The broad GABAergic comparator's scatter is visually inflated relative to the interneuron panels. Fix: use one common scale.

### MINOR 3. Figure S5 error bars are unlabeled.

Tex line 441 (S5 caption). The error bars are the null 95 percent confidence intervals (`results/followup-mature-concordance-figure-data.csv`; mouse within null CI 0.1012 to 0.2055, observed 0.2105 just above the upper bound), but the figure does not label them. Fix: label the bars as the null 95 percent CI. The aging-biology round 4 review also flagged the S5 null label ("random-gene null" versus the joint donor-permutation null); this review adds the error-bar label.

### MINOR 4. Figure S4 shows only Yazar counts; Gao is absent and the Sood cohort split is unstated.

Tex line 439 (S4 caption). All counts are Yazar blood (Harries 5 after, 7 before; Peters 16 with and without; Sood 4 before, 0 after; de Magalhaes 0). Gao counts (all zero) are not shown. Sood "4 groups pass before adjustment" does not say 3 are Gao heart and 1 is Yazar blood (`reports/followup-mature-block10.md` lines 126 to 131). Fix: add the Gao zeros and the cohort split.

### MINOR 5. Between-subject log variance ratio point estimates are epsilon-sensitive and the tex does not state the epsilon.

Tex line 243. Duplicate rows in `results/followup-mature-between-subject-variance-yazar.csv` differ by eps_used (1e-08 versus 0.001): naive CD4 batch LVR 0.278 versus 0.213 (23 percent change); naive CD8 batch 0.456 versus 0.402. q-values are stable (CD4 batch q 0.02799, depth-only q 0.007996; CD8 batch q 0.007996, depth-only q 0.007996). The tex quotes the eps 1e-08 values without stating the epsilon. Fix: state the epsilon used. Also: Yazar A1 slope is depth-sensitive (full depth minus 0.0261, p 3.55e-06; at 1,500 molecules minus 0.0819, p 2.1e-10; at 500 minus 0.1192, p 2.5e-12; at 3,000 minus 0.0600, p 0.0176); the tex reports only full depth.

### MINOR 6. Methods line 53 says "the P1 analysis" where it means P1m.

Tex line 53: "so the P1 analysis neither confirms nor refutes that distinct claim". P1 is the young-to-old signature; the paper's analysis is P1m (mature-to-old). Fix: "so the P1m analysis" or "so this analysis".

### MINOR 7. Figure S3 title says "distances shrink" but the second row is within-donor spread, a different estimand.

Tex line 437 (S3 caption). Slopes verify (Yazar minus 0.026, p 3.6e-06; minus 0.078, p 2.0e-13; Gao minus 0.008, p 0.77 and minus 0.062, p 0.085). Fix: retitle to cover both rows, or split the title.

### MINOR 8. The debiased naive CD4 slope exists only in report prose, in no committed CSV.

Tex lines 33, 255, 403. The debiased slope (0.000313, p 0.000638) and the 56 percent calibration bias verify against `reports/followup-mature-block1.md` line 228 (primary 0.00071, p 5.174e-08, q 8.279e-07; reweighted 0.00080, p 1.165e-06; calibration accounts 56.4 percent, induced slope 0.000401, p 8.40e-07; debiased 0.000313, p 0.000638; debiased plus depth 0.000295, p 0.00106; plus log(N) 0.000286, p 0.00660; plus N plus depth 0.00024, p 0.01893; overlap weights null after plus log(N) 0.00024, p 0.1722 and plus N plus depth 0.00016, p 0.3464). But grep for 0.000313 in `results/` finds no committed CSV with the debiased slope. The report itself (line 232) says "That is not proof of a biological dispersion rise; it is a weaker calibration story than the original had, and the verdict rests on the frozen rule, not on a resolved mechanism." Fix: commit the debiased-slope computation as a CSV so the headline number is machine-reproducible.

### MINOR 9. Survivor-bias Spearman value is missing.

Tex line 309. "human Pearson plus 0.287, p 0.0069, permutation p 0.008; Spearman null" gives no Spearman r or p. Fix: give the Spearman value or delete the clause. Mouse Pearson minus 0.057, p 0.689, N 51 versus plus 0.250, p 0.0498, N 62 before restriction, verified as stated.

## What would change the verdict

1. Fix BLOCKER 1 and BLOCKER 2 (reconcile the concordance family and the PBMC implementations; correct the affected lines in Results, Discussion, Conclusions, Methods, the S5 caption, and Appendix B).
2. Fix MAJOR 2, 3, 4, 5, 6, 7, and 13 (Figure 2 partial residuals or caption, lung fixed-effects slopes, one scan family, mature propeller p, Gao Shannon specifications, S2 units, stale tag). These are all bounded edits, not new science, except MAJOR 8 (the P1m leave-one-donor-out check, which is new analysis; if it cannot be run before submission, state the gap explicitly in Methods and Limitations).
3. MAJOR 1 (624 family) is a design finding, not a fixable typo. The paper discloses the selection timing. If the author keeps 624 primary, the manuscript must say plainly that the 45-pass count is not a valid FDR-controlled result under a prespecified family, rather than presenting it as the primary composition result.
4. The remaining MAJOR and MINOR items are bounded wording, caption, and citation fixes.

## Notes

1. Overlap with the sibling round 4 reviews: the aging-biology review covered Gao Shannon (its MAJOR 1), the S2 mixed-units axis (its MAJOR 3), the S5 null label (its MAJOR 4), and the stale Data Availability tag; it accepted the three-test concordance family as primary without checking whether the family exists in the analysis record. The top-tier review noted the Figure 2 fitted lines are nearly flat and the Figure 2 caption ambiguity. The items in this review that are new are: BLOCKER 1 (the three-test family does not exist), BLOCKER 2 (PBMC 35 versus 36), MAJOR 2 (Figure 2 raw-OLS sign flips and the Sim 2021 donors), MAJOR 3 (lung fixed-effects slopes), MAJOR 4 (scan family incoherence), MAJOR 5 (propeller scope mix), MAJOR 7 (S2 rank/numbering conflict), MAJOR 8 (P1m leave-one-donor-out), MAJOR 9 (gene-fraction single-test family), MAJOR 10 (spleen 3,000 BH survival), MAJOR 11 (CD8 gate uncited), MAJOR 12 (lung-versus-brain slope), MAJOR 14 (S1 label), and MINOR 1, 5, 6, 8, 9.
2. The project Q and A inherits the tex's numbers and needs the same fixes: line 55 (PBMC 0 of 35), line 69 (Q3, propeller p 0.0528 and lung study-aware numbers), line 97 (Q5, three-test family q 0.066), line 111 (Q6, Gao Shannon depth-by-cell-type omission), line 139 (Q7, CD8 gate contamination uncited).
3. Items not verified in this review (stated as limits): the Thomson PBMC batch model parameter count (15 parameters on 20 donors, task item 4) was not confirmed from a committed file; the rarefaction "150 of 159 rows unchanged" (tex line 239) was not independently recomputed; the expanded 79-unit transport counts (P1m 4/79, P0m 5/79, rescue 0/48, tex line 273) were not independently recomputed from a CSV. None of these is load-bearing for the verdict above.
