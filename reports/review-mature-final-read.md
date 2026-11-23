# Final read review: mature Agescape manuscript

- Manuscript of record: `report/plos-one/final-plos-one.tex`, 676 lines, commit `28259b1` on `origin/master`.
- Scope of this review: the two new results incorporated after `abedea7` (held-out validation and Bayesian refits), checked number by number against `reports/followup-mature-heldout-validation.md`, `reports/followup-mature-bayesian.md`, and the committed `results/followup-mature-heldout-*.csv` / `results/followup-mature-bayesian-*.csv` files. Full tex read; all 4 main figures and all 6 SI figures read as rendered PNGs; project Q&A read at `~/workspace/goals/agescape-manuscript-submission/files/agescape-questions-and-answers.md`.
- Method: every new number recomputed from its source CSV with Python (no hand copying); the tex diff `abedea7..28259b1` inspected hunk by hunk; figure labels compared against caption and body values.

## Verdict

READY WITH MINOR FIXES

## Counts

1. BLOCKER: 0
2. MAJOR: 8
3. MINOR: 7

No number in the tex is wrong and no new number lacks a source. Every MAJOR item below is a wording, framing, or figure-label fix; none requires new analysis. All 8 MAJOR items should be fixed before submission, then re-checked in one pass.

## MAJOR items

### M1. Abstract drops the interneuron caveats the body carries

1. Tex line: 33 (Abstract).
2. What the tex says: "Interneuron correlations held in two cohorts."
3. What the body says: Results (interneuron section) states the moderated mature check did not verify (frozen recovery gate failed) and a SEA-AD test in 15 donors without dementia did not replicate the program; Discussion and Conclusions repeat both points.
4. Source: `reports/followup-mature-interneuron.md` (moderated gate failed; SEA-AD no replication); Q&A Q9 states the same.
5. History: the `abedea7` abstract said "Interneuron correlations held in two cohorts; moderated check failed, third failed." The final pass trimmed both caveats out of the abstract.
6. Exact fix: restore the `abedea7` wording in the abstract, within the 300-word cap (abstract is currently 275 words, so 25 words of headroom exist): "Interneuron correlations held in two cohorts; moderated check unverified, third cohort failed."

### M2. Abstract held-out sentence omits the strict-bar failure

1. Tex line: 33 (Abstract).
2. What the tex says: "Held-out validation found real mouse age signal: 26 of 26 slopes positive, 22 of 26 passed the calibrated test."
3. What the body says: Results line 275 states "None of the 26 passed the strict bar of both tests together. The smallest composite q value was 0.0595 and the median composite p value was 0.0996." The held-out report verdict is explicit: "Against the pre-specified strict bar (both tests q below 0.05 in the predicted direction): 0 of 26. The bar is NOT met."
4. Source: `results/followup-mature-heldout-cross-organ-transport.csv` (26 TESTED, 0 composite pass, 22 empirical pass, 0 both pass; min composite q 0.0595286; median composite p 0.0995991); `reports/followup-mature-heldout-validation.md`, Verdict section.
5. Judgment: the Results paragraph itself is balanced (it states the strict-bar failure in the same paragraph), and the "real signal" reading is defensible under the paper's own hierarchy (Methods: empirical q is the calibrated evidence number, composite q an uncalibrated screening statistic). But the abstract carries only the positive half, and the abstract is what most readers see. The pre-specified bar was not met, and the abstract must say so.
6. Exact fix: append to the abstract sentence: "though none of the 26 passed the pre-specified strict bar of both tests together (smallest composite q 0.0595)." Also qualify the Conclusions line 419 ("Held-out validation showed real mouse age signal") the same way, or leave Conclusions as is only if the abstract carries the qualifier, since Conclusions currently repeats the unqualified claim.

### M3. "Donor-held-out checks were run on the cleaned signature" misstates what the splits did

1. Tex lines: 115 (Methods, signature derivation), 267 (Results, signature section), 409 (Limitations). Same sentence in three places.
2. What the tex says: "Donor-held-out checks were run on the cleaned signature. ... A second split that kept marrow eligible matched in 32 of 32 checks, including 7 of 7 marrow checks."
3. What the sources say: both splits re-derived a fresh signature after holding out one mouse per side and checked direction in the held-out mice. Variant 3 derived an 84-gene signature (22 up, 62 down) and matched 23 of 23; variant 3b derived a 687-gene signature (565 up, 122 down; 609 of 687 in the cleaned signature) and matched 32 of 32, including 7 of 7 marrow. The cleaned 1,914-gene signature itself was never scored in held-out mice. The tex Methods held-out subsection (line 119) states the procedure correctly ("Donor-split variants derived a fresh signature after holding out one mouse per side"), so the tex contradicts itself.
4. Source: `reports/followup-mature-heldout-validation.md`, Variants 3 and 3b; `results/followup-mature-heldout-donor-split.csv` (23 rows, all `direction_match=True`); `results/followup-mature-heldout-donor-split-3b.csv` (32 rows, all True; 7 Marrow, 10 Kidney, 15 Lung).
5. Judgment: the numbers (23 of 23, 32 of 32, 7 of 7 marrow) are all correct. What is wrong is the characterization: the checks test donor-robustness of the derivation procedure, not of the cleaned signature's gene list. That distinction matters because the Limitations conclusion ("transport failure reads as failure of a real mouse signal") leans on donor robustness.
6. Exact fix (all three places): replace "Donor-held-out checks were run on the cleaned signature" with "Donor-held-out checks re-derived the signature after holding out one mouse per side". In Results line 267 and Limitations line 409, add "using a freshly derived 687-gene signature (609 genes shared with the cleaned signature)" after "32 of 32 checks".

### M4. Methods claims no human Bayesian drift numbers are reported; Appendix A reports one

1. Tex lines: 145 (Methods, Bayesian refits) and 467 (Appendix A).
2. What the tex says, line 145: "The human Bayesian drift model did not converge, so no human Bayesian drift numbers are reported."
3. What the tex says, line 467: "The full-range Bayesian refits gave a mouse depth-adjusted posterior mean of plus 0.0782 per 10 percent of maximum lifespan (94 percent interval plus 0.0399 to plus 0.1163, probability of a positive effect 0.999) and a human depth-adjusted posterior mean of plus 0.0133 (interval minus 0.0081 to plus 0.0344, probability 0.88)."
4. Source for the Appendix A numbers: `results/h3-human-bayesian-global-posterior.csv`, row species Human, model `global_depth_adjusted`, prior 0.5, scope species, unit ALL: posterior mean 0.0133191, HDI94 low -0.0080778, high 0.0343608, prob positive 0.881625. Mouse: `results/h3-cross-species-bayesian-posterior.csv`, species Mouse, model `depth_adjusted`, prior 0.5, scope species: mean 0.0782384, HDI 0.0398726 to 0.1163410, prob 0.999. Both numbers are sourced and correctly rounded; the problem is the Methods sentence, not the numbers.
5. Judgment: the non-convergence finding concerns the new mature human PyMC model (`reports/followup-mature-bayesian.md` section 2.5: all five specs failed convergence, max R-hat 1.07 to 2.48, min bulk ESS 5 to 44, 26 to 163 divergences; `results/followup-mature-bayesian-human-drift-diagnostics.csv`, all `convergence_ok=False`). The Appendix A number comes from the older full-range fit, which is a different model. The new Methods sentence is unscoped, so as written it is false.
6. Exact fix, line 145: replace "The human Bayesian drift model did not converge, so no human Bayesian drift numbers are reported." with "The mature human Bayesian drift model did not converge, so no mature human Bayesian drift numbers are reported; Appendix A reports the full-range human fit only."

### M5. Fig 1 donor panels show different q values than the caption and body

1. Tex line: 231 (Fig 1 caption). Figure: `report/plos-one/figures-v2/fig1-composition.png`, lower donor panels.
2. What the figure shows: club cell q 0.012, multiciliated q 0.028, L6b q 0.012, oligodendrocyte q 0.012, OPC q 0.012.
3. What the caption and body say: caption line 231 gives q 0.011, 0.026, 0.011, 0.011, 0.011; Results and Conclusions give refit-family q 0.0109 (club, L6b, oligodendrocyte, OPC) and 0.0258 (multiciliated).
4. Source of the discrepancy: `scripts/build_optimization_figures.py`, function `fig1_composition`, comment "q values are the original 431-test family (0.01178, 0.02785)" and code `"q": float(crow["q_431_family"])`. The figure panels use the 431-family q; the caption and body use the refit-family q on the 399-row mature scan. Both sets of numbers are correct under their own family; they disagree with each other for the same five rows.
5. Exact fix: rebuild Fig 1 donor panels with the refit-family q values (0.011, 0.026, 0.011, 0.011, 0.011), and add the family name to the panel labels or caption. Do not change the caption numbers; the refit family is the paper's stated bar ("The bar is the refit family on this 399-row mature scan").

### M6. Bayesian ROPE verdict for Yazar scatter/Shannon omits the source report's own qualifications

1. Tex line: 307 (Bayesian Results subsection).
2. What the tex says: "Yazar scatter and Yazar Shannon entropy sit inside the pre-stated region of practical equivalence, with region probabilities 1.000 and 0.976. Their directions are certain but their magnitudes are negligible by that region."
3. Numbers verified: `results/followup-mature-bayesian-narrowing.csv` gives Yazar scatter prob_rope 1.0 and Yazar Shannon (primary) prob_rope 0.9756; region half-width 0.05; 0.05 / 0.0138 = 3.62, so "about 3.6 times the smallest effect Yazar scatter was powered to detect at 80 percent power, 0.0138 standard deviations per decade" is correct.
4. What the source report adds and the tex drops: `reports/followup-mature-bayesian.md` section 2.2 states the cumulative move over 50 years is about 0.10 SDs (scatter) and 0.15 SDs (Shannon), and "the robust statement is the direction probability, not the ROPE label", with the caveat that "the ROPE verdict for Yazar is driven by the wide band, not by uninformative data" (ROPE 3.6x the MDE80). The tex reports the 3.6x ratio but not the cumulative move or the "robust statement is direction, not ROPE label" qualification.
5. Judgment: as written, "magnitudes are negligible" invites a reader to conclude the Bayesian analysis refutes the narrowing magnitude, while the source analysis says the per-decade increment is small but the cumulative move over the studied span is 0.10 to 0.15 SDs and the ROPE label is band-dependent.
6. Exact fix, line 307: after "Their directions are certain but their magnitudes are negligible by that region." add "Over 50 years the cumulative moves are about 0.10 standard deviations for scatter and 0.15 for Shannon entropy, so the robust Bayesian statement is the direction probability, not the ROPE label."

### M7. Naive CD4 debiased Bayesian result omits its ROPE qualification

1. Tex line: 305 (Bayesian Results subsection).
2. What the tex says: "The debiased slope was 0.000313, interval 0.000132 to 0.000492, probability 0.9997."
3. Numbers verified: `results/followup-mature-bayesian-naive-cd4.csv`, debiased row: posterior median 0.00031318, CI 0.00013208 to 0.00049249, prob_direction 0.999655. Correct.
4. What the tex omits: the same CSV row gives prob_rope 0.376945. `reports/followup-mature-bayesian.md` section 2.3 states the debiased slope "sits almost exactly on the ROPE edge (0.05 x outcome SD = 0.000284 per decade), so the data cannot separate 'real but small' from 'practically null' for the debiased arm" (at double ROPE width prob_rope 0.997; at half width 0.032).
5. Judgment: the debiased slope is the bias-corrected estimate of the paper's one surviving biological candidate, so presenting its direction probability 0.9997 without the ROPE-edge qualification overstates how firmly the Bayesian refit supports it. The primary arm (prob_rope 0.001) needs no such qualification.
6. Exact fix, line 305: after "probability 0.9997" add "(region probability 0.377: the debiased slope sits on the region edge, so the data cannot separate a real but small effect from a practically null one)."

### M8. Project Q&A is stale against the new manuscript

1. File: `~/workspace/goals/agescape-manuscript-submission/files/agescape-questions-and-answers.md` (outside the repo; no tex change fixes this).
2. What the Q&A says: header states "Source: the manuscript at commit abedea7". Q0 verdict states "P1m is a three-organ probe built on an all-male old group, and the paper states both limits." The Q&A contains no entry for held-out validation and no entry for the Bayesian refits.
3. What the manuscript now says: the signature is marrow-dominated (1,803 of 1,914 genes supported only by marrow), not a three-organ probe; donor-held-out checks and cross-organ held-out validation exist; Bayesian refits exist.
4. Exact fix: update the Q&A header source commit to 28259b1; correct Q0 to marrow-dominated framing; add a held-out-validation entry (26 tested, 22 empirical pass, 0 both-tests, strict bar not met) and a Bayesian-refits entry (headline directions agree; Yazar narrowing magnitudes inside ROPE; ten transport nulls undecided; mature human drift model did not converge).

## MINOR items

### m1. "Rose monotonically" overstates the intermediate-age pattern

1. Tex line: 275.
2. What the tex says: "Intermediate-age mice scored with the existing signature rose monotonically with age. Per-mouse mean scores were minus 0.724 at 18 months, minus 0.718 at 21 months, minus 0.593 at 24 months, and minus 0.467 at 30 months."
3. Numbers verified: `results/followup-mature-heldout-intermediate-ages-per-mouse.csv` gives group means -0.72441 (n=4), -0.71832 (n=2), -0.59293 (n=4), -0.46736 (n=4); Spearman rho 0.733, p 0.00285 (recomputed). The four means are monotonic, but the 18-to-21 step is +0.006, far smaller than within-group spread (18-month mice range -0.877 to -0.582; 21-month mice -0.765, -0.672), and the source report says "The 18-to-21 step is flat (-0.725 to -0.718); the rise appears at 24 and 30 months."
4. Exact fix: replace "rose monotonically with age" with "group means rose with age (flat from 18 to 21 months, rising at 24 and 30 months)".

### m2. Gao scatter direction probability 1.0000 comes from an unclustered fit

1. Tex line: 303.
2. What the tex says: "Gao scatter fell with median minus 0.0719 per decade, interval minus 0.1053 to minus 0.0383, probability 1.0000. ... The Gao scatter interval is narrower than the clustered frequentist interval because the conjugate fit does not cluster by donor."
3. Numbers verified: `results/followup-mature-bayesian-narrowing.csv`, Gao scatter row: median -0.0718651, CI -0.1053477 to -0.0383092, prob 0.99999. Correct.
4. What the tex discloses but does not quantify: `reports/followup-mature-bayesian.md` section 2.2 Limitation gives the widths: Bayesian 95% width 0.067 versus clustered frequentist width about 0.137 (cluster SE 0.0349); Yazar 0.014 versus about 0.019. The disclosure names the cause but not that the 1.0000 probability itself is optimistic for Gao (frequentist p is 0.041 on clustered errors).
5. Exact fix: append "so the Gao direction probability of 1.0000 is optimistic; the clustered frequentist interval is about twice as wide."

### m3. Held-out BH family is ambiguous

1. Tex line: 119 (Held-out Methods subsection).
2. What the tex says: "Both tests used Benjamini-Hochberg correction within the tested family."
3. What the source shows: `scripts/run_mature_heldout_validation.py` lines 420-421 applies BH per fold. Recomputed: smallest composite q 0.0595 reproduces only under per-fold BH (fold 2, 8 tested units: 0.00744 x 8 = 0.0595). Pooled BH across all 26 tested units would give 0.00744 x 26 = 0.193. The held-out report uses the same "within the tested family" wording, so the ambiguity is inherited, but a reader recomputing across 26 units gets a different smallest q.
4. Exact fix: replace "within the tested family" with "within each fold".

### m4. "Agreed on direction on every headline result" is overbroad, and the drift section says "agrees" despite opposite signs

1. Tex lines: 301 (Bayesian subsection lead) and 291 (drift section).
2. What the tex says, line 301: "Bayesian refits agreed with the frequentist direction on every headline result." Line 291: "A Bayesian hierarchical refit agrees: depth-adjusted posterior mean plus 0.0389 per 0.4 years (94 percent interval minus 0.0998 to plus 0.1773, probability of a positive effect 0.71)."
3. What the sources say: frequentist mouse drift slope is minus 0.0137 (p 0.816); Bayesian posterior mean is plus 0.0389. `reports/followup-mature-bayesian.md` section 2.5 states "the two point estimates have opposite signs while both intervals are wide and include zero". The Bayesian subsection covers composition, narrowing, CD4, transport, and concordance, where directions do agree, so line 301 is defensible only if read as scoped to that subsection; line 291's "agrees" is not defensible on sign.
4. Exact fix: scope line 301 to "Bayesian refits agreed with the frequentist direction on every headline result in this subsection." Change line 291 "A Bayesian hierarchical refit agrees:" to "A Bayesian hierarchical refit is likewise null but with the opposite sign on the point estimate:".

### m5. Fig S5 human observed bar does not match the legend color

1. Figure: `report/plos-one/figures-si/s5-concordance.png`.
2. What the figure shows: the legend defines dark blue as "Observed" and grey as "Joint donor-permutation null mean", but the Human within-cohort observed bar is rendered gold/brown while the two mouse observed bars are dark blue.
3. Numbers in the figure are correct (mouse within-tissue 0.2105 vs null 0.1465, q 0.044, p 0.022; mouse cross-tissue 0.0604 vs null 0.0297, p 0.103; human 0.1793 vs null 0.2126, q 0.815, p 0.815; footnote null 95% interval 0.101 to 0.206).
4. Exact fix: recolor the human observed bar dark blue to match the legend (or add the gold color to the legend).

### m6. Abstract drops the rescue denominator

1. Tex line: 33 (Abstract).
2. What the tex says: "unfiltered signature rescued none". The `abedea7` abstract said "unfiltered signature rescued none of 48".
3. What the body says: Results states the unfiltered signature passes both tests in zero of the 29 core units where the cleaned signature passes neither, and zero of the 48 units in the expanded 79-unit comparison.
4. Exact fix: restore "unfiltered signature rescued none of 48" (or "none of the 29 core units") in the abstract.

### m7. Abstract omits the Bayesian concordance disagreement

1. Tex line: 33 (Abstract).
2. What the tex says: "Bayesian refits agreed on headline directions and left transport nulls undecided."
3. What the body says, line 307: concordance Bayesian interval minus 0.0682 to plus 0.1804 includes zero, "against a frequentist q of 0.0440", and "The Bayesian answer is more conservative because it propagates the sampling error of the observed correlation."
4. Exact fix: extend the abstract sentence to "Bayesian refits agreed on headline directions, left transport nulls undecided, and did not establish the concordance difference (interval includes zero, against frequentist q 0.0440)." The abstract has 25 words of headroom under the 300-word cap.

## Verified clean: every new number checked against its source

1. Held-out cross-organ, `results/followup-mature-heldout-cross-organ-transport.csv` (26 TESTED, 10 UNTESTABLE, 36 rows): all 26 slopes positive; min 0.05917, max 0.49523, median 0.15967 (tex: 0.059 to 0.495, median 0.160); empirical pass 22, composite pass 0, both pass 0; min composite q 0.0595286 (tex 0.0595); median composite p 0.0995991 (tex 0.0996). Units per tested row have 6 to 8 mice total (3 to 5 residual degrees of freedom after age plus sex). All match tex lines 33 and 275.
2. Per-fold signature sizes: `results/followup-mature-heldout-fold1-leave-Kidney-out-genes.csv` 1,858 rows; fold2 1,906 rows; fold3 111 rows. Match tex line 275 and the held-out report table (fold 1: 1,109 up, 749 down; fold 2: 1,115 up, 791 down; fold 3: 36 up, 75 down).
3. Intermediate ages, `results/followup-mature-heldout-intermediate-ages-per-mouse.csv` (14 mice): means -0.72441, -0.71832, -0.59293, -0.46736 (tex: -0.724, -0.718, -0.593, -0.467); Spearman rho 0.73311, p 0.002853 (tex: 0.733, p 0.00285). Match tex line 275.
4. Donor splits: `results/followup-mature-heldout-donor-split.csv` 23 rows, 23 direction matches (9 Kidney, 14 Lung, 0 Marrow); `results/followup-mature-heldout-donor-split-3b.csv` 32 rows, 32 matches (10 Kidney, 15 Lung, 7 Marrow). Match tex lines 115, 267, 275, 409 (see M3 for the wording issue; the counts are right).
5. Marrow framing: `reports/followup-mature-heldout-validation.md` lines 79-80: 1,848 of 1,914 genes have at least one marrow supporting unit; 1,803 supported ONLY by marrow; 88 any kidney support; 26 any lung support; best kidney unit 36 genes, best lung unit 9. Match tex lines 115, 267, 351, 409.
6. Bayesian composition, `results/followup-mature-bayesian-composition.csv` (5 rows): club median -0.07702, CI -0.11436 to -0.03956, prob 0.99996 (tex: -0.0770, -0.1144 to -0.0396, 1.0000); multiciliated -0.04503, -0.06969 to -0.02029, 0.99975 (tex: -0.0450, -0.0697 to -0.0203, 0.9998); L6b 0.01933, 0.00988 to 0.02882, 0.99996 (tex: 0.0193, 0.0099 to 0.0288, 1.0000); oligodendrocyte 0.01689, 0.00842 to 0.02540, 0.99993 (tex: 0.0169, 0.0084 to 0.0254, 0.9999); OPC -0.01166, -0.01744 to -0.00585, prob 1.0 (tex: -0.0117, -0.0174 to -0.0059, 1.0000). All five intervals exclude zero. Match tex line 301.
7. Bayesian narrowing, `results/followup-mature-bayesian-narrowing.csv` (6 rows, primary rows used): Gao scatter -0.07187, -0.10535 to -0.03831, prob 1.0 (tex: -0.0719, -0.1053 to -0.0383, 1.0000); Yazar scatter -0.01965, -0.02680 to -0.01250, prob 1.0, prob_rope 1.0 (tex: -0.0196, -0.0268 to -0.0125, 1.0000); Gao Shannon primary -0.09882, -0.19982 to 0.00120, prob 0.9736 (tex: -0.0988, -0.1998 to 0.0012, 0.9736); Yazar Shannon primary -0.02970, -0.04990 to -0.00971, prob 0.9982, prob_rope 0.9756 (tex: -0.0297, -0.0499 to -0.0097, 0.9982). Match tex line 303. Note the Gao Shannon depth-by-cell-type row (-0.10015, CI -0.19316 to -0.00805, prob 0.9835) is not reported in the tex; the tex correctly uses the primary row.
8. Bayesian naive CD4, `results/followup-mature-bayesian-naive-cd4.csv`: primary median 0.00071318, CI 0.00045708 to 0.00096676, prob 1.0 (tex: 0.000713, 0.000457 to 0.000967, 1.0000); debiased median 0.00031318, CI 0.00013208 to 0.00049249, prob 0.999655 (tex: 0.000313, 0.000132 to 0.000492, 0.9997). Match tex line 305 (see M7 for the missing ROPE qualification; debiased prob_rope 0.376945, primary prob_rope 0.000505).
9. Bayesian transport, `results/followup-mature-bayesian-transport.csv` (12 rows): MSSM microglia median -0.00812933, CI -0.01326272 to -0.00297290, prob 0.998945 (tex: -0.00813, -0.01326 to -0.00297, 0.9989); PsychAD microglia median -0.01403036, CI -0.02114483 to -0.00698567, prob 0.999955 (tex: -0.01403, -0.02114 to -0.00699, 1.0000). The PsychAD borrowed-scale disclosure in tex line 305 matches `reports/followup-mature-bayesian.md` section 2.4 (score SD borrowed from MSSM microglia, 0.0515; half/double scale sensitivity leaves P(<0) at 0.9999/1.0000). The ten null units were verified programmatically: none has its full 95% interval inside its ROPE (closest is Yazar block8 CD14-low monocyte, prob_rope 0.9345, interval still extends outside the ROPE), so "ten transport null units are undecided" is correct. Match tex lines 305 and 307.
10. Bayesian concordance, `results/followup-mature-bayesian-concordance.csv`: delta_r median 0.0632008, CI -0.0682288 to 0.1804427, prob 0.831865, prob_rope 0.373915 (tex: 0.0632, -0.0682 to 0.1804, 0.8319, 0.374); frequentist q 0.0440. Match tex line 307.
11. Bayesian drift ROPE: `results/followup-mature-bayesian-drift.csv` gives mouse posterior mean 0.0389439, 94% HDI -0.0997789 to 0.1773341, prob positive 0.7068, prob_rope 0.4494122 (tex line 291: plus 0.0389, interval minus 0.0998 to plus 0.1773, probability 0.71, ROPE probability 0.449). Match.
12. ROPE half-width and MDE80 ratio: half-width 0.05 outcome SDs per the Bayesian report section 1 (ROPE definition) and tex Methods line 145; 0.05 / 0.0138 = 3.62 (tex line 307: "about 3.6 times"). The 0.0138 MDE80 value is stated in `reports/followup-mature-bayesian.md` section 2.2. Match.
13. No human Bayesian drift numbers from the mature model appear in the main text: the drift section (line 291) reports only frequentist human drift (-0.0071, p 0.611); the Bayesian subsection (lines 299-307) reports no human drift numbers. The only human Bayesian drift numbers anywhere in the tex are the full-range Appendix A numbers (line 467), which are the subject of M4.
14. Banned terms in prose: "P1m" (uppercase) 0 occurrences; "P0m" 0; "three-organ" 0; "MDE80" 0; hyphenated "both-pass" 0. Lowercase "p1m" survives only inside file paths and URLs (Data Availability links such as `reports/followup-mature-p1m-derivation.md`, line 427; S1 caption filename `figures-si/s1-p1m-multi-unit-support.pdf`, line 449). Those are paths, not prose; no change needed.
15. 310-of-765 figure: appears only in Appendix A, line 465, attached to P1 ("In leave-one-donor-out refits, 310 of 765 P1 gene-by-unit instances changed sign"), never attached to the mature signature. Correct; matches the held-out report's instruction that the figure "still belongs to P1, not P1m".
16. Abstract: 275 words by whitespace count, under the 300-word cap. No abstract number contradicts the body except the omissions in M1, M2, m6, m7.
17. 624-primary disclosure intact: stated in Results (composition section), Discussion, Limitations, and Conclusions: the family was chosen after both correction outcomes were available, FDR control is not guaranteed, and the two promoted tests fail sensitivity checks.
18. Diff `abedea7..28259b1` on the tex: changed hunks are the abstract, signature derivation, the new held-out Methods subsection, Bayesian Methods, signature Results, the new held-out Results paragraph, the drift ROPE sentence, the new Bayesian Results subsection, Discussion transport paragraph, Limitations, and Conclusions. Every changed number traces to a source listed above. No frequentist number changed except through the intended additions; the abstract trims in M1, M2, and m6 are the only silent deletions, and each is flagged above.

## Figure-by-figure consistency

1. Fig 1 (composition): forest panels match the caption (9 mouse rows open, 5 human rows filled, bladder untestable and not shown). Donor-panel slopes match the caption and body (-0.080, -0.047, +0.020, +0.017, -0.012). Donor-panel q labels do NOT match the caption/body (see M5).
2. Fig 2 (narrowing): all five panels match the caption: (a) Gao scatter -0.0723, p 0.041; (b) Yazar scatter -0.0197, p 6.69e-05; (c) pooled center shift +0.0400, p 2.99e-05, within-cell-type -0.0044, p 0.27 (caption p 0.266, rounding only); (d) Gao Shannon -0.1001, p 0.0539, not significant; (e) Yazar Shannon -0.0298, p 0.00365. The low cluster near -3 to -4.5 in panel (a) is visible as captioned.
3. Fig 3 (transport): rendered bars match the caption exactly: cleaned 2 both / 13 random-only / 4 composite-only / 29 neither; unfiltered 3 / 15 / 5 / 25. The in-figure title (0 of 29 rescue) matches the caption. Consistent with the new held-out text: the figure is about the core 48 units, not the held-out 26, and the tex does not confuse the two.
4. Fig 4 (interneuron): TRMT11 labeled in each of the three subtype panels; rings mark same-direction FDR genes in HBCC; broad GABAergic comparator grey with no rings. Consistent with the caption and body. Note the figure cannot show the moderated-check failure or the SEA-AD failure; those live only in text, which is why M1 matters.
5. Fig S1 (signature support): bars sum to 1,914 (1,212 + 341 + 163 + 66 + 48 + 27 + 24 + 18 + 9 + 5 + 1); 702 in 2 or more, 361 in 3 or more, 132 in 5 or more, 57 both directions. Matches the caption exactly and is consistent with the marrow-dominated framing (most genes rest on a single group).
6. Fig S2 (fraction scan): 10 highlighted rows (2 mouse: ranks 8, 15; 8 human: ranks 6, 21, 26, 51, 81, 219, 227, 237) with slopes matching the body (-0.143, -0.100, -0.080, -0.047, -0.038, -0.021, -0.012, +0.017, +0.020, +0.023). Matches the caption.
7. Fig S3 (between-subject shrinkage): slopes -0.026, p 3.6e-06 and -0.078, p 2.0e-13 match the caption and body (-0.0261, p 3.55e-06; -0.0776, p 1.97e-13); Gao nulls -0.008, p 0.77 and -0.062, p 0.085 match.
8. Fig S4 (published signatures): Harries Yazar 5 after / 7 before adjustment; Peters 16 / 16; Sood 4 before (3 Gao, 1 Yazar), 0 after; de Magalhaes 0. Matches the caption and body.
9. Fig S5 (concordance): observed vs null values and q/p labels match the caption and body; footnote null 95% interval 0.101 to 0.206 is consistent with observed 0.2105 sitting just above it. See m5 for the legend-color defect. Note the figure predates the Bayesian refits and therefore shows no Bayesian interval; the Bayesian concordance disagreement (line 307) is text-only, which is fine, but a reader of S5 alone sees only the frequentist pass.
10. Fig S6 (ectopic E1/E2): Gao primary -0.063, p 0.47; Gao stratum -0.127, p 0.031; Yazar primary -0.019, p 0.073; Yazar stratum -0.033, p 0.0031; donor counts 104 / 925 / 103 / 919. Matches the caption and body.

## Answers to the three judgment questions in the task

1. Does the held-out verdict paragraph overstate? The Results paragraph (line 275) does not overstate on its own: it reports 26 of 26 positive, 22 of 26 calibrated passes, 0 of 26 on the strict bar, smallest composite q 0.0595, median composite p 0.0996, and the power explanation, then draws the "real signal, transport failure not noise failure" inference, which follows the paper's own calibrated-vs-screening hierarchy and the held-out report's verdict. The overstatement is in the Abstract (M2) and Conclusions, which carry the positive half without the strict-bar qualifier, and in the repeated Discussion phrasing ("not a noise failure" stated twice in one paragraph, line 351) without noting the fold-3 weakness (only 7 of 11 marrow-held-out units pass even the calibrated test, including failures in granulocyte and granulocytopoietic cells, the top derivation units).
2. Does the abstract's trimmed wording drop any caveat the body carries? Yes: the interneuron moderated-check failure and SEA-AD failure (M1), the held-out strict-bar failure (M2), the rescue denominator (m6), and the Bayesian concordance disagreement (m7). The abstract is 275 words, so all four fit within the 300-word cap.
3. Is any Bayesian statement misleading? Three are, by omission rather than by wrong numbers: the Yazar ROPE "negligible" verdict without the cumulative-span qualification (M6); the debiased CD4 direction probability without its ROPE-edge qualification (M7); and the subsection lead "agreed on direction on every headline result" read against the drift sign disagreement (m4). The Gao clustering caveat and the PsychAD borrowed scale are disclosed in the tex (lines 303 and 305) and are not misleading as written, though the Gao probability 1.0000 remains optimistic (m2). The ten-null-units-undecided claim, the concordance Bayesian interval against frequentist q 0.0440, and the mouse drift ROPE probability 0.449 are all correct and fairly stated.
