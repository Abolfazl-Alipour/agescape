# Review round 3: top-tier journal reviewer (Nature Aging / Nature / Cell standard)

Manuscript: `report/plos-one/final-plos-one.tex` at commit 683a74d, tag `mature-to-old-rewrite-2026-11-21`, repo `github.com/Abolfazl-Alipour/agescape`.
Reviewer standard: Nature Aging / Nature / Cell. Date: 2026-11-21.
Checker result on this commit: PASS, abstract 270 words, four figures, 68 cited references.
I read the tex in full (721 lines), read Figures 1 and 4 as images, and verified every number cited below against a result CSV in the same repo. No number in this review is copied from another report without a CSV check.

## Verdict: REJECT at top-tier standard

This is a careful, honest, and rigorous negative result. It is not a Nature Aging, Nature, or Cell paper. The mature-to-old reframing is the scientifically right move, and it is also what kills the paper at this bar: it exposes that the original finding was young-driven, and it leaves nothing positive of top-tier scale behind. The paper's central negative (no transport, no drift, no uniform response) rests on a probe the paper itself says is sex-confounded and unstable, with no positive control anywhere in the study. The three narrow signals are real in the statistical sense and thin in the biological sense, and two of the three are arguably technical. No revision of the text alone fixes that. The work belongs in a specialty journal as a constraint paper, or as a resource, not at a top-tier journal.

Single most fatal weakness: the headline negative is uninterpretable because the probe is confounded and unvalidated. P1m is derived from four 30-month mice, all male, against two female and two male 18-month mice; it is fourteen times the size of the young-to-old signature under the same rules for a cause the paper says is not established; and the study contains no progeroid model, no accelerated-aging model, and no held-out-mouse validation. A null transport result from a probe like that cannot distinguish species divergence from a broken probe, and the paper says so itself (Discussion, lines 339 and 399; Limitations, line 387). A top-tier journal does not publish a null built on an admittedly confounded instrument.

Single most salvageable strength: the rigor and honesty of the negative audit. The batch and pool artifact analysis, the Block 10 published-signature benchmark, the disclosure that the 624-test family was selected after both correction outcomes were available, the reporting of negative-direction passes instead of hiding them, and the Limitations section are all exemplary. That material, plus the mature composition table, is a publishable constraint paper in a specialty venue. It is not a Nature Aging paper, because a constraint needs a discovery to constrain, and this paper's discovery is that its own earlier finding was young-driven.

## How I checked the numbers

I verified these against the CSVs, and they all check out: P1m 1,914 genes, 1,119 up and 795 down (`results/followup-mature-p1m-genes.csv`); overlap 52 shared, 27 reversed, 87 lost, 1,862 gained (`results/followup-mature-p1m-p1-overlap.csv`); 624 family 45 passes, 17 core plus 28 added, 26 of the added passes fly (`results/rev3-composition-global-bh-624.csv`); promoted tests vein endothelial p 0.00272, core q 0.0690, 624 q 0.0390, total-cell q 0.0627, propeller p 0.0528, and collecting duct p 0.00248, core q 0.0667, 624 q 0.0377, total-cell q 0.4231, propeller p 0.0588 (same file plus `results/reviewer-phase1-composition-total-adjusted.csv` and `results/reviewer-phase1-composition-propeller.csv`); mature composition five survivors with slopes plus 0.020, plus 0.017, minus 0.012, minus 0.080, minus 0.047 (`results/followup-mature-composition-human.csv`); mouse drift minus 0.0137, p 0.816, and human drift minus 0.0071, p 0.611 (`results/followup-mature-drift-mouse.csv`, `results/followup-mature-drift-human.csv`); microglia composite q 0.00457 and 0.0497, empirical q 0.0135 and 0.00262 (`results/followup-mature-p1m-transport.csv`, `results/followup-mature-p1m-transport-null.csv`); Block 5 counts P1m 2/13/4/29 and P0m 3/15/5/25, strict rescue 0 (`results/followup-mature-p1m-block5-summary.csv`); naive CD4 slope 0.00071, p 5.17e-08, q 8.28e-07 (`results/followup-mature-block1-residual-dispersion.csv`); naive CD8 scatter 0.0498, p 5.92e-05, q 4.74e-04, 344 donors, 17,878 cells (`results/followup-mature-fulldepth-yazar-cell-drift.csv`); MSSM depth slopes minus 0.1039, minus 0.0877, minus 0.0839, minus 0.0524 (`results/followup-mature-block6-mssm-depth.csv`); interneuron correlations 0.7523, 0.6299, 0.6996 and both-FDR counts 4, 8, 41 with TRMT11 the only shared gene (`results/followup-mature-interneuron-replication.csv`, `results/followup-mature-interneuron-genes.csv`); Harries All_Cells slope plus 0.1985, q 7.51e-47, and Peters plus 0.09675, q 2.60e-33 (`results/followup-mature-block10-benchmark-results.csv`); mouse concordance 0.2105 against null 0.1465, p 0.0220, and human 0.1793 against 0.2126, p 0.8152 (`results/followup-mature-concordance-mouse.csv`, `results/followup-mature-concordance-human.csv`); Yazar between-subject pooled ratio 1.543, p 3.46e-05, residualized 1.578, p 0.137 (`results/followup-mature-between-subject-distance-variance-yazar.csv`); PsychAD pooled P1m slope plus 0.000508, p 0.000547 (`results/followup-mature-psychad-m10-cohort-summary.csv`); GTEx whole blood plus 0.0119, q 1.29e-22 (`results/followup-mature-organ-transports-c-gtex-track2-p1.csv`). The arithmetic is sound. The problems below are about what the numbers mean and what is missing, not about arithmetic errors.

## Answers to the assigned questions

Is this paper important enough for Nature Aging, Nature, or Cell? No. A null result can clear the top-tier bar only when it overturns a major claim, closes a field, or carries a positive discovery of comparable weight. This paper overturns its own earlier finding, not a field consensus. Its target, the universal hallmarks claim of Tyshkovskiy and colleagues (ref 11), is not actually tested, because the paper never benchmarks their chronological or mortality clocks and says so (Introduction, line 54). The paper refutes a four-criterion construct it built itself, a construct stricter than the claim it argues against. That is a straw man at top-tier bar, however honestly disclosed.

Is the central claim novel, or a negative result dressed up? It is a negative result, presented honestly as one. The novelty is methodological: a mature-to-old re-audit that shows how much of the aging-transcriptome literature is young-driven, batch-driven, or pool-driven. That is genuinely useful. It is not a discovery about aging biology. The title says exactly what was found, and the abstract does not oversell the nulls. The dressing-up, such as it is, sits in the three narrow signals, which are promoted as the paper's positive content and do not carry that weight (see B3 and B4).

Does the mature-to-old reframing rescue the paper or expose that the original finding was young-driven? It exposes it, completely. The original three transport passes died (composite q 0.688, 0.625, 0.321). The mouse drift collapsed from plus 0.1028 (p 0.000004, young reference) to plus 0.0072 (p 0.67) against an 18-month reference, and the restricted mature estimate is minus 0.0137 (p 0.816). The signature is a different gene set: 52 of 139 shared, 27 of those reversed. The reframing rescued the paper's honesty and destroyed its positive result. That trade is correct science and fatal at top-tier bar.

Are the three narrow signals strong enough to carry a paper at this bar? No. Naive CD4 dispersion is a single cell type in one cohort whose own audit leaves the mechanism unresolved (see B4). Naive CD8 scatter is a single cell type in one cohort, and its share fall fails the batch model (q 0.283); the GTEx whole-blood positive is bulk tissue, which the paper itself says conflates composition with within-cell change. MSSM depth decline adjusts only for sex, has no postmortem interval, brain bank, RNA quality, batch, or technical principal components, uses four nested units rather than four independent confirmations, and the paper states it is "not evidence of altered biological transcriptional output". A top-tier paper cannot headline a signal it simultaneously disclaims as a biological finding. The three signals also have nothing to do with each other: no shared mechanism, no shared tissue, no shared direction. They are three leftovers, not a story.

Is the scope a strength or a sign of unfocused fishing? Both, and at this bar the second reading wins. The breadth (624 tests, dozens of cohorts, ten blocks, four non-mammalian species) is the design, and the donor-level discipline is real. But 28 of the 45 composition passes are fly and worm, 26 of them fly; the exploratory species contribute nothing to the mature-to-old claim and are explicitly never claimed as mature-restricted, yet they sit in the primary correction family and change core q values. A top-tier reviewer reads that as scope used to rescue two core tests that fail every sensitivity check. The mature table, 399 rows with its own families, is the honest primary result and should have been primary from the start.

Are the figures publication-quality, and do they show the story or hide it? Technically competent, editorially misaligned. Figure 1 is a clean forest plot, and Figure 4 is a clean four-panel replication scatter. But Figure 1 shows a different analysis set from the primary family and omits the two promoted tests, by its own caption (line 247). The mature composition table, the paper's actual primary composition result, has no figure. The drift null, the paper's central result, has no figure in the main text; the only drift figure (Figure 2) sits in Appendix A and shows the young-reference analysis the paper disowns. The transport failure, the P1m versus P1 overlap, and the death of all original passes have no figure at all. The figures show side results and hide the story.

Is the writing clear, is the abstract honest, and does the title match the results? The writing is unusually clear for a paper this large, and the Limitations section is exemplary. The title matches the results. The abstract is honest about the nulls but omits the single most informative fact about the one positive transport result: both microglia passes run in the negative direction. A reader of the abstract alone could conclude the signature transports in microglia. It does not; the score falls with age there.

## Issues

### BLOCKER 1. The central negative rests on a confounded, unstable probe with no positive control.

Sections: Methods, lines 120 to 131; Results, lines 237 to 241; Discussion, line 339; Limitations, lines 387 and 399.
Evidence: all four 30-month mice are male against two female and two male 18-month mice; P1m is fourteen times the size of P1 under the same rules for a cause the paper says is not established; only 36 of 151 mouse units were eligible for derivation, from kidney, lung, and marrow only; three PyDESeq2 dispersion warnings occurred with the affected units not recorded; classic aging-gene positive controls were not significant in the top derivation unit. The study has no progeroid model, no accelerated-aging model, and no held-out-mouse validation (Limitations, line 387; Predictions, line 361). The paper states the consequence itself: "Transport failure therefore remains a property of an unstable and sex-confounded probe as well as possible species divergence; derivation quality is not excluded" (line 339). At top-tier bar, a null that the authors themselves cannot assign to biology versus instrument is not publishable as a central claim. Fixing this needs new data or a new derivation, not new text: a mature signature derived with both sexes at old age, held-out validation across derivation methods, and a positive control showing the assay can detect transport when transport exists.

### BLOCKER 2. There is no discovery of top-tier scale.

Sections: Results, lines 271 to 277; Abstract, line 37; Conclusions, line 411.
Evidence: the three narrow signals are (a) one cell type in one cohort with an unresolved mechanism, (b) one cell type in one cohort whose share result fails batch adjustment, and (c) a depth decline the paper disclaims as biological evidence. The interneuron replication, the strongest positive in the paper, yields one shared gene, TRMT11, on the unmoderated test only, with the moderated check unverified and SEA-AD failing to replicate. The Block 10 benchmark is a useful audit, not a discovery. Nothing here changes how aging biology is studied, prioritizes a target, or opens a research program, except the negative audit itself, which belongs in a specialty venue.

### BLOCKER 3. The naive CD4 headline number is the un-debiased slope; the debiased number never appears in Results or abstract.

Sections: Abstract, line 37; Results, lines 271 to 277; Methods, line 212.
Evidence: the abstract and Results report slope 0.00071, p 5.17e-08, q 8.28e-07, which I verified in `results/followup-mature-block1-residual-dispersion.csv`. The paper's own Block 1 report (`reports/followup-mature-block1.md`, Section 5.7) shows finite-N calibration accounts for 56.4 percent of that slope, leaving a debiased slope of 0.000313 (p 0.000638), with a persistent log cell-fraction coupling (beta minus 0.00111, p 4.48e-07), an IPW arm that steepens rather than removes the log cell-count decline, and null unadjusted, rarefied, gene-fraction, and shared-trend arms. The report's own interpretation: "That is not proof of a biological dispersion rise; it is a weaker calibration story than the original had, and the verdict rests on the frozen rule, not on a resolved mechanism." The Methods state the calibration checks are binding on robustness (line 212). The Results section on naive CD4 mentions none of this. Reporting the largest estimate while the audit's debiased estimate sits in a repo report is selective reporting, and at top-tier bar it reads as such. The debiased slope, the calibration share, and the coupling must appear next to the primary slope in Results and abstract, or the signal must be demoted from the headline list.

### BLOCKER 4. The abstract omits the direction of the only transport passes.

Section: Abstract, line 37.
Evidence: the abstract says "only microglia passed both transport tests". The Results (lines 259 to 269) show both pass in the negative direction: PsychAD microglia slope minus 0.00143 per year, MSSM microglial minus 0.00083 per year. The Yazar both-pass units are likewise negative (All_Cells minus 0.0151, monocyte minus 0.0154 per decade), and the GTEx Track 1 passes are all negative. The direction is the finding: the signature does not transport, it anti-transports in the only units that pass. An abstract that omits this misleads a skimming reader, including an editor triaging the paper. This is fixable in one sentence, but at top-tier bar an abstract that needs the fix is a symptom of how the paper frames its one positive transport result.

### MAJOR 1. Figure 1 does not show the paper's result.

Section: Results, line 247, Figure 1 caption.
Evidence: the caption states the 25 displayed results "are distinct from the 45 passes in the 624-test primary family and do not display the two promoted primary tests", and that "the mature-to-old fraction table is reported in the text". I inspected the figure: a clean forest plot of full-range sensitivity-model slopes. The paper's lead figure therefore shows neither its primary family result nor its mature result. A top-tier paper's Figure 1 must show its headline claim. Replace it with the mature composition table as a forest plot, including the five survivors and the nine mouse losses with their inflated intervals, and move the current figure to the supplement.

### MAJOR 2. The two central results have no main-text figure.

Sections: Results, lines 253 to 269; Appendix A, Figure 2.
Evidence: the drift null (mouse minus 0.0137, human minus 0.0071) and the transport failure (all original passes dead, only negative microglia passes) are text-only. The only drift figure is Figure 2 in Appendix A, which plots young-reference slopes for 62 mouse units, the analysis the main text disowns. A reader cannot see the null. Add a main-text figure for mature drift by cohort and a transport figure showing original versus mature composite q for the core units, with the two negative microglia passes marked.

### MAJOR 3. The 624-test family is post-hoc scope-stretching, and the paper knows it.

Sections: Methods, lines 141 to 143; Results, lines 243 to 251.
Evidence: the tex states the family "was selected after both correction outcomes were available" and is "not a preregistered or outcome-blind scope". Of the 45 passes, 26 are fly, 2 are worm, 11 are mouse, 3 are MSSM brain, 2 are lung, 1 is kidney (verified in `results/rev3-composition-global-bh-624.csv`). The two promoted core tests pass only in the full family (core q 0.0690 and 0.0667) and fail total-cell adjustment (q 0.0627 and 0.4231) and propeller (p 0.0528 and 0.0588). The paper discloses all of this, which is to its credit, but disclosure does not cure the design: the family exists to rescue two tests that fail every sensitivity check. Make the mature 399-row table primary with its own correction families and its own figure, and report the 624 family as a sensitivity, or drop it.

### MAJOR 4. Two specific percentages in the Discussion have no traceable source.

Section: Discussion, line 347.
Evidence: "About 83 percent of gene hits did not persist after pool adjustment and concordance fell by about half in the 981-donor blood cohort, while over 96 percent did not persist after cohort adjustment in the 99-donor cohort." I searched `reports/followup-mature-*.md` and the Batch A result CSVs and could not reproduce either number. The gene-test counts I found point elsewhere: Gao 0 of 145,778 tests significant at mature against 295 of 94,837 originally (`reports/followup-mature-organ-transports-batch-a.md`, line 82), which is a different comparison. No table in the manuscript or repo supports 83 percent or 96 percent as stated. At top-tier bar, specific quantitative claims in the Discussion need a cited table or they come out.

### MAJOR 5. The paper does not test the literature claim it positions itself against.

Section: Introduction, line 54; ref 11.
Evidence: the Introduction invokes the universal hallmarks claim of Tyshkovskiy and colleagues (Nature, 2026) and then states: "We did not benchmark their chronological or mortality clocks, so the P1 analysis neither confirms nor refutes that distinct claim." The four-criterion construct is then defined as "deliberately more stringent than shared hallmarks, modular conserved signatures, or transportable clocks". A paper that defines a stricter test than its target's claim, fails it, and concedes it did not test the target's actual claim has not adjudicated anything in the literature. Either benchmark the Tyshkovskiy clocks directly in these cohorts, or remove the framing that this paper speaks to universal hallmarks and present it as a constraint on one specific operational construct.

### MAJOR 6. The interneuron replication rests on one gene, unmoderated, with a failed replication cohort.

Section: Results, lines 279 to 283; Figure 4.
Evidence: correlations are strong (VIP r 0.7523, pvalb 0.6299, sst 0.6996, verified in `results/followup-mature-interneuron-replication.csv`), but same-direction both-FDR counts are 4, 8, and 41, and the three-subtype intersection is TRMT11 only. The moderated mature intersection is unverified because the frozen recovery gate failed. SEA-AD (15 donors) did not replicate the program. PsychAD and the public MSSM and HBCC cohorts share overlapping donors, so the PsychAD arm is not independent replication of the public cohorts (Limitations). Figure 4 does not label TRMT11, the one gene the text asks the reader to care about. At top-tier bar, a one-gene intersection on the unmoderated test, with a failed replication cohort and an unverified moderated check, is a suggestive finding, not a headline result. Label TRMT11 in the figure, state the SEA-AD failure in the abstract-adjacent Results text rather than only in the section's last sentence, and do not present the moderated result as corroboration.

### MAJOR 7. A third of the transport family is untestable because source files are missing.

Section: Methods, line 129; Results, lines 259 to 269.
Evidence: 82 of 252 units are untestable at 30-plus, 71 because their source files are absent locally. I verified the count in `results/followup-mature-p1m-transport-untestable.csv` (82 rows). The missing units include Kidney Normal Adult, RADC, SEA-AD, and Kidney SN. A transport conclusion drawn with a third of the family missing for a recoverable reason (files not downloaded) is incomplete at top-tier bar. Restore the four missing cohorts or state in the abstract that the transport conclusion covers 170 of 252 units.

### MAJOR 8. P1m is a three-organ signature, not a mouse signature.

Sections: Methods, lines 120 to 131; Results, lines 237 to 241.
Evidence: only 36 of 151 mouse units were eligible for derivation (11 kidney, 12 lung, 13 marrow). Bladder, tongue, mammary gland, skin, and thymus had no 30-month mice at all. The signature the paper calls "the mature mouse signature" comes from three organs and eight mice. Its failure to transport to human brain, heart, or retina may reflect that it was never a whole-mouse signature. The text states the eligibility numbers but continues to call P1m a mouse signature throughout. Rename it or qualify every transport claim with its three-organ origin.

### MAJOR 9. MSSM depth decline should not be a headline signal.

Sections: Results, lines 271 to 277; Abstract, line 37; Limitations.
Evidence: the four MSSM neuron depth slopes are verified (minus 0.1039, minus 0.0877, minus 0.0839, minus 0.0524 per decade), but the models adjust only for sex; postmortem interval, brain bank, RNA quality, processing batch, and technical principal components were unavailable. The broad GABAergic unit contains the other three, so the four are nested, not independent. The paper states the finding is "an age-associated captured-library-depth decline, not evidence of altered biological transcriptional output". A signal the authors disclaim as biological should not be one of three headline findings in the abstract. Demote it to a technical observation, or remove it from the headline list and let the two blood signals carry the positive content.

### MINOR 1. The mature composition q values do not name their family.

Section: Results, lines 243 to 251.
Evidence: the five survivors are reported with q 0.0118, 0.0118, 0.0118, 0.0118, and 0.0279. In `results/followup-mature-composition-human.csv` these are the 431-family q values; the refit-family q values are 0.0109 and 0.0258, and the within-file q_bh values are 0.00091 and 0.0102. All three pass, so nothing is misreported, but the tex never states which family the printed q belongs to, while the Methods say the mature table has its own correction families. Name the family in the text.

### MINOR 2. Figures are raster, not vector.

Section: Results, Figures 1 to 4.
Evidence: the four figures are PNG files at 2100 by 1980 to 2700 by 1680 pixels. That resolution is adequate at print size, but Nature journals expect vector figures for line art and scatter plots. Supply vector versions.

### MINOR 3. Space devoted to a nominal mouse concordance result.

Section: Results, lines 307 to 311.
Evidence: mouse within-tissue concordance is p 0.0220 on 14 mice, with a joint null mean already at 0.1465, carried by limb muscle and marrow only; human concordance is null (p 0.8152). The paper frames it correctly as nominal only. It still occupies a full Results subsection and a main-text figure (Figure 3) for a result the paper does not claim. Condense to one paragraph and move Figure 3 to the supplement.

### MINOR 4. The 683a74d follow-up patch.

Note, not a manuscript defect: commit 683a74d, made after the rewrite, restored the naive CD4 signal to the abstract and corrected the Figure 4 caption. The rewrite had initially dropped a headline signal from the abstract ("Two narrow signals survived"). That a headline finding can disappear in a rewrite and need a patch is a process fragility the author should note, not a flaw in the current text, which now lists three signals.

## What would have to change for Nature Aging

Concretely, in order of importance:

1. Re-derive the mature signature with both sexes at old age, and validate transport in held-out mouse donors across derivation methods. Without this, the central negative stays uninterpretable.
2. Add a positive control: a progeroid or accelerated-aging model, or a signature already known to transport, showing the assay detects transport when it exists. Without this, every null in the paper is compatible with an insensitive assay.
3. Replicate at least one narrow signal in an independent cohort with a prespecified analysis: naive CD8 scatter with batch adjustment and cell-state resolution, or TRMT11 with same-cell multiome data. One replicated positive discovery is the minimum for this bar.
4. Benchmark the Tyshkovskiy chronological and mortality clocks directly in these cohorts, or drop the universal-hallmarks framing.
5. Make the mature 399-row composition table primary, give it a figure, and demote the 624 family to a sensitivity. Restore the four missing transport cohorts.
6. Fix the abstract (state the negative direction of the microglia passes; report the debiased naive CD4 slope alongside the primary slope), add main-text figures for the drift null and the transport failure, label TRMT11 in Figure 4, and cite or remove the 83 percent and 96 percent Discussion numbers.

If items 1 to 3 cannot be done, the right move is not a longer revision but a different venue: a specialty aging journal, or a negative-results and methods venue, where this paper's rigor and honesty are the point rather than a liability. As a constraint paper with a clear statement of what it rules out and what it does not, it is strong. At Nature Aging, Nature, or Cell, as written, it is not competitive.

## The three issues that matter most

1. The central negative rests on a sex-confounded, unstable probe with no positive control (BLOCKER 1). Everything else in the paper depends on this.
2. There is no discovery of top-tier scale: three thin, disconnected narrow signals, two arguably technical (BLOCKER 2, BLOCKER 3, MAJOR 9).
3. The paper does not adjudicate the literature claim it invokes, because it never tests the Tyshkovskiy clocks and defines a stricter construct than its target (MAJOR 5).
