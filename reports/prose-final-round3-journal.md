# Prose final round 3: journal style and structure review

Scope: prose and structure only, not science. No number, verdict, caveat, disclosure, or citation is challenged here. The 624 family, the marrow framing, the held-out and Bayesian results, and the literature framing are taken as settled, and no new analysis is proposed.

Manuscript: `report/plos-one/final-plos-one.tex`, base commit `d15befd` (tag `mature-to-old-final-base-2026-11-23`), 815 lines, read end to end. The tex was not edited.

The Abstract is fixed (owner's plain style, 293 words, no p-values, no jargon). Its content is not reviewed and no change to it is proposed. Two Abstract-body mismatches are flagged below (issues M12 and M13) with exact line numbers; both recommended fixes are body-side only.

## 1. Section balance (approximate word counts, body text only)

1. Abstract, lines 32-35: about 293 words, fixed.
2. Introduction, lines 36-63: about 1,794 words.
3. Materials and Methods, lines 64-235: about 6,851 words, 19 subsections.
4. Results, lines 236-435: about 8,943 words, 13 subsections.
5. Discussion, lines 436-469: about 2,990 words, no subheadings of its own.
6. Predictions, lines 470-485: about 375 words, nested as a Discussion subsection.
7. Limitations, lines 486-531: about 1,078 words, 46 lines, nested as a Discussion subsection.
8. Conclusions, lines 532-539: about 704 words, 8 lines.
9. Supporting Information, lines 564-650: 87 lines, including Appendices A-C.
10. References, lines 651-815: 165 lines.

The imbalance that matters is not Methods versus Results (both are long, as expected for this paper) but Discussion versus its own subsections: Discussion proper is 34 lines while Limitations alone is 46 lines, so the paper's caveat list is longer than its interpretation. Conclusions at 8 lines is the thinnest interpretive section in the paper. Within Results, the transport subsection is 30 lines while concordance and published signatures are 8 lines each, which is a second imbalance (issue M6).

## 2. BLOCKER

### B1. Figures are numbered out of order of first appearance: lines 300-322

1. Problem: the interneuron figure is Figure 4, cited in text at line 300 and placed at line 302, while the transport figure is Figure 3, cited at line 320 and placed at line 322, so a reader meets Figure 4 before Figure 3, which no published paper does.
2. Fix: renumber the interneuron figure to Figure 3 and the transport figure to Figure 4, and change the in-text citations to match (line 300 "(Figure 4)" becomes "(Figure 3)"; line 320 "(Figure 3)" becomes "(Figure 4)"); check that no other in-text reference to Figure 3 or Figure 4 exists before committing the swap.

## 3. MAJOR

### M1. Abstract and body tell the story in different orders: Abstract line 34; Results lines 240-394

1. Problem: the Abstract orders the findings as fractions, concentration, signature and microglia, GTEx whole blood, drift, concordance, interneuron, published signatures, then the 45-of-624 count, while the body orders them as fractions (240), concentration (256), limited signals (282), interneuron (296), signature (306), transport and microglia (316), drift (346), concordance (356), published signatures (364), so interneuron moves from near the end in the Abstract to fourth in the body, and the held-out check the Abstract mentions in its third sentence does not appear in the body until subsection 11 (372).
2. Fix: restructure the Results spine to match the Abstract order: fractions, concentration, limited signals, signature, held-out validation, transport, drift, concordance, interneuron, published signatures, Bayesian refits, disorder tests, exploratory arms, follow-up cohorts; then rewrite the roadmaps (see M3) to that order.

### M2. Held-out validation is placed after the transport failure it is needed to interpret: lines 316-379

1. Problem: the reader reads the whole transport failure (316-345) before learning at lines 372-379 that the probe carries real mouse signal, yet the Discussion (line 458) relies on that fact to call the failure a transport failure rather than a noise failure, so the transport section is uninterpretable on first reading in the way the paper itself needs it to be read.
2. Fix: move the Held-out validation subsection (372-379) to immediately before the transport subsection (316), directly after the signature subsection (306-315); keep its text unchanged.

### M3. Both roadmaps contradict the actual Results order: lines 58, 238

1. Problem: the Introduction roadmap (line 58) promises "Bayesian refits, published signatures, and disorder tests follow", but the actual order is published signatures (364), held-out validation (372), Bayesian refits (380), disorder tests (394), while the Results introduction (line 238) gives a second, different and incomplete order statement that omits published signatures, disorder tests, exploratory arms, and follow-up cohorts entirely.
2. Fix: delete the order sentence in the Results introduction (line 238) and keep only the Introduction roadmap (line 58), rewritten to the restructured order from M1; do not keep two roadmaps.

### M4. Discussion has no subheadings, restates Results, and is shorter than its own Limitations: lines 436-531

1. Problem: Discussion proper (436-469, 34 lines, one unbroken block) opens with three paragraphs that restate Results (438, 440, 442, 444) instead of interpreting them, while Limitations (486-531, 46 lines) is longer than the Discussion it belongs to, so the paper's interpretive core reads thinner than its caveat list.
2. Fix: cut paragraphs 440, 442, and 444 to a single short synthesis paragraph; add subheadings to the Discussion (findings, literature, theory, alternative explanations); move the batch-adjustment paragraph (466) into Limitations, where it belongs as a caveat rather than an interpretation.

### M5. Conclusions summarizes instead of concluding: lines 532-539

1. Problem: Conclusions paragraph 1 (line 534) reprints Results numbers (slopes and q values for all five fraction shifts, scatter, entropy, and off-program counts), paragraph 3 (line 538) is a third list of open items duplicating Predictions (470-485) and Limitations (486-531), and the take-home sentence is buried as the final sentence after that list, so the section restates rather than concludes.
2. Fix: rewrite Conclusions to two paragraphs with no slopes, p-values, or q values (name the five shifts and the concentration finding only); delete the "What stays open" list in line 538 and replace it with a pointer to Predictions and Limitations; end on the take-home claim plus one sentence stating what the field should do differently as a consequence.

### M6. The transport subsection is six findings in one subsection: lines 316-345

1. Problem: one 30-line subsection contains the core cohorts, the rescue test, reverse transport, the organ sweep, the PsychAD arm, and the PBMC cohorts, while concordance (356-363) and published signatures (364-371) get 8 lines each, so the most load-bearing negative result in the paper is the hardest to navigate and the two thin subsections look like afterthoughts next to it.
2. Fix: split the transport subsection into two subsections, "The mature signature does not transport in the core cohorts" (core cohorts, rescue test, reverse transport) and "The organ sweep finds no positive-direction pass" (organ sweep, PsychAD arm, PBMC cohorts), with the existing opening sentence (line 318) rewritten as the roadmap for the two.

### M7. The "Limited signals" subsection mixes three kinds of result and buries positives at its end: lines 282-295

1. Problem: one subsection holds a surviving biological candidate (naive CD4 dispersion), a technical observation (MSSM captured-depth decline), and an inconclusive scatter signal (naive CD8), and its final paragraph (line 292) buries three positive findings (CD8 donor-cloud spreading, the Thomson CD8 share fall, and the GTEx whole-blood positive) that belong with the analyses that produced them.
2. Fix: split the subsection into "One surviving candidate and one technical observation" (CD4 dispersion and MSSM depth) and "Naive CD8 scatter is inconclusive" (scatter verdict only); move the donor-cloud sentences from line 292 to the disorder subsection (394), the Thomson share-fall sentence to the fractions subsection (240), and the GTEx sentence to the transport subsection (316).

### M8. Methods does not follow the Results spine it promises: lines 66, 120-137, 160-171

1. Problem: the Methods introduction (line 66) promises "the tests in spine order", but Methods puts signature derivation and held-out validation (120-137) before composition (138-149) while Results puts composition first and the signature in the middle, and published-signature benchmarking has no Methods subsection of its own (it is one paragraph inside the targeted-dispersion subsection, 160-171) while exploratory arms have none beyond the dataset list.
2. Fix: reorder the Methods subsections into Results spine order (datasets, composition, concentration, limited-signal tests, signature, held-out, transport, drift, concordance, published signatures, Bayesian, disorder, exploratory, verdict rules, families); add one short dedicated subsection each for published-signature benchmarking (moved out of lines 160-171) and for the exploratory arms.

### M9. The Introduction does Methods' job and repeats itself: lines 40-62

1. Problem: the Introduction states the operational definition of aging (lines 50-51), the four uniformity criteria (line 52), the transcriptomic scope (line 54), and the design (line 56), then gives a roadmap (line 58) and two further closing meta-paragraphs (lines 60 and 62) that restate what the Introduction just did, and the cross-species divergence evidence (mouse methylation-clock sites distinct from human; muscle sharing under 5 percent of pathways) is stated twice, at line 40 and again at line 44.
2. Fix: delete the repeated clock and muscle sentences from either line 40 or line 44 (keep one copy); move the design paragraph (line 56) to the Methods introduction; delete the two closing meta-paragraphs (lines 60 and 62) and keep line 58 as the single roadmap.

### M10. The spine is restated five times across the manuscript: lines 58, 234, 238, 434-435, 468

1. Problem: the same spine summary appears in the Introduction roadmap (58), the Methods closing (234, "Taken together, this design establishes..."), the Results introduction (238), the Results closing (434-435, "Taken together, the Results establish the spine in order..."), and the Discussion closing (468, "What the Discussion establishes..."), so a reader is told the structure five times and the closers add no information the section has not already given.
2. Fix: keep one roadmap (Introduction line 58, corrected per M3); cut the Results closing paragraph (434-435) to one sentence; cut the Discussion closing paragraph (468) entirely; keep the Methods closing (234) only if it is shortened to a transition sentence without restating the findings.

### M11. The disorder subsection packs five analyses under a title that covers two of them: lines 394-411

1. Problem: one 18-line subsection covers residual dispersion, gene coordination, categorical signature patterns, depth per cell, and stress and survivor-bias tests, but its title ("No general rise in disorder, coordination loss, or survivor bias") does not cover the depth analysis, and the stress and survivor-bias material is a different question from disorder.
2. Fix: split into two subsections, "No general rise in disorder or coordination loss" (dispersion, coordination, categorical patterns) and "No survivor-bias signal" (stress analysis and cross-type localization, including the Part A and Part B material); move the depth-per-cell sentences to the limited-signals area (282) or give them their own short subsection.

### M12. Abstract-body mismatch 1 (Abstract fixed; fix is body-side only): Abstract line 34; Results lines 320-344

1. Problem: the Abstract states without qualification that "Only brain immune cells called microglia passed both tests", but the body counts three further both-test passes in the organ sweep (GTEx whole blood, line 332; retinal pigment epithelium in the negative direction and Lake kidney endothelium under the unfiltered signature only, line 334), and the body's own closing sentence (line 344) repeats the unqualified "Only microglia pass both tests" in the same breath as "whole blood stands as a bulk-composition candidate", so the Abstract's claim and the body's sweep count contradict each other as written.
2. Fix: in the body only, scope the claim wherever it appears: change line 318 to state the order and scope ("in the core cohorts") and change line 344 to "Only microglia pass both tests in the core cohorts, in the negative direction, and whole blood stands as a bulk-composition candidate"; do not change the Abstract.

### M13. Abstract-body mismatch 2 (Abstract fixed; fix is body-side only): Abstract line 34; Results lines 350, 354

1. Problem: the Abstract states without qualification that "Whole profiles did not drift from the mature reference", but the body qualifies the mouse result as "null but undecided under the region rule" with an interval that "includes meaningful continuation in either direction and does not establish absence" (line 350, repeated at line 354), so the Abstract's clean null overstates the mouse result as the body states it.
2. Fix: in the body only, add one sentence to the drift subsection closing (line 354) scoping the no-drift claim explicitly ("the no-drift claim is established in human; the mouse result is undecided, not a null"), so a reader comparing Abstract and body sees the scope; do not change the Abstract.

## 4. MINOR

### m1. Formulaic bookends in every Results and Methods subsection: lines 66-232 (Methods), 240-435 (Results)

1. Problem: all 13 Results subsections open "The question is..." and close "This subsection establishes...", and 19 Methods subsections open "This subsection states...", which reads as a report template rather than a journal paper.
2. Fix: keep the question openers, delete all 13 "This subsection establishes..." closing sentences (each repeats its subsection's first paragraph), and vary or delete the "This subsection states..." openers in Methods.

### m2. The Datasets subsection is dominated by tables and inline identifiers: lines 68-115

1. Problem: the 48-line Datasets subsection is mostly Table 1 and Table 2 plus 11 inline Census dataset identifiers, so the cohort logic a reader needs to follow Results is buried under accession strings.
2. Fix: move the 11 Census identifiers to Supporting Information and keep only cohort names, donor counts, and age windows in the text.

### m3. The verdict-rule table sits in Methods before the analyses it governs are described: lines 196-217

1. Problem: Table 3 (22 lines of unit-level and dataset-level decision rules for the ten follow-up analyses) appears in Methods before the Bayesian, nonlinear, and families subsections, so the reader must carry ten verdict rules through the rest of Methods before reaching the Results that use them.
2. Fix: move Table 3 to Supporting Information and keep a one-paragraph summary of the verdict-rule approach in Methods.

### m4. The follow-up cohorts subsection duplicates Appendix B: lines 426-435, 614-648

1. Problem: the Results follow-up subsection (426-435) is a cohort-by-cohort catalogue of nulls and blocked cohorts that repeats Appendix B (614-648), which lists the same cohorts with the same verdicts.
2. Fix: cut lines 426-435 to a three-sentence summary (what was tested, what was blocked, that no headline finding changed) and point to Appendix B for the cohort-by-cohort verdicts.

### m5. Exploratory and follow-up arms sit after the Results' logical close: lines 412-435

1. Problem: the exploratory arms (412-425) and follow-up cohorts (426-435) come after the disorder subsection (394-411), so the Results spine ends with catalogues rather than with the constraint tests, and the "Taken together" summary (434-435) has to re-cover material that belongs earlier.
2. Fix: combine lines 412-435 into a single "Exploratory and follow-up arms" subsection placed before the Bayesian refits subsection (380), after the published-signatures subsection (364).

### m6. Data Availability cites an older tag than the manuscript base: lines 542, 546, 548

1. Problem: Data Availability points the workflow manifest, the tagged source, the source archive, and every linked report and script at tag `mature-to-old-round4-optimization-2026-11-22`, while this manuscript's base is commit `d15befd` (tag `mature-to-old-final-base-2026-11-23`), so a reader following the links reaches an older version.
2. Fix: update every link in lines 542, 546, and 548 to the final base tag.

### m7. Appendices A-C are nested inside Supporting Information instead of standing alone: lines 564-650

1. Problem: Appendix A (young-reference comparisons, line 580), Appendix B (follow-up cohorts, line 614), and Appendix C (primary composition family, line 649) sit inside the Supporting Information section, although Appendix A holds comparisons the main text cites (line 186) and Appendix B is cited as the verdict record.
2. Fix: promote Appendices A-C to top-level sections placed after Supporting Information and before References.

### m8. Theory and the organ-schedules argument arrive late in the Discussion: lines 462-464

1. Problem: the theory paragraph (line 462) comes after the transport and concordance paragraphs (458, 460) although the theories were set up in the Introduction (line 42), and the organ-schedules paragraph (line 464) introduces a new explanatory hypothesis with five new literatures (organ proteomic clocks, ageotypes, organ-specific hallmark timing) in the second half of the Discussion.
2. Fix: move the theory paragraph (462) to directly after the findings synthesis (after line 444) and shorten the organ-schedules paragraph (464) to its three paper-internal points plus a single pointer sentence to the cited literatures.

### m9. The slogan "aging is local and plural" opens Conclusions and also closes it: lines 534, 538-539

1. Problem: the Conclusions opens with "From maturity to old age, aging is local and plural" (line 534) and ends with the same sentence, so the section is circular rather than progressive.
2. Fix: rewrite the Conclusions opening to state the conclusion with its scope (mature-to-old window, transcriptomic level, gene level) instead of repeating the slogan, and keep the slogan only as the final sentence.

### m10. The Results introduction forward-references an undefined probe: line 238

1. Problem: line 238 says the 1,914-gene signature "is defined in the signature subsection below", so a reader must hold an undefined probe through the fractions, concentration, and limited-signals subsections before reaching its definition at line 306.
2. Fix: add one sentence at line 238 defining the probe in place (1,914 genes, 18- versus 30-month mice, marrow-dominated) and keep the signature subsection where the M1 restructure puts it.

### m11. "Reproducibility and reporting" is too thin to be a subsection: lines 230-235

1. Problem: the final Methods subsection is three sentences (what was retained, what was not exported, how families were labeled) followed by a "Taken together" closer, which is a paragraph, not a subsection.
2. Fix: fold lines 230-235 into the "Statistical families and multiplicity" subsection (226-229) as a closing paragraph and delete the separate subheading.

### m12. The 624-family selection-timing caveat is repeated four times: lines 252, 440, 486-531, 534

1. Problem: the statement that the 624-test family was chosen after both correction outcomes were available, so false discovery rate control is not guaranteed, appears in the fractions subsection (line 252), the Discussion (line 440), Limitations, and Conclusions (line 534).
2. Fix: state it at first mention (line 252) and in Limitations only; delete it from the Discussion paragraph (440) and from Conclusions (534).

## 5. What holds up structurally

1. The narrative arc is sound in outline: problem (uniform versus plural), approach (mature-to-old window, donor-level inference, four criteria), findings (local fractions, concentration, limited candidates), constraints (transport failure, no drift, mouse-only concordance), conclusion (local and plural).
2. The Introduction sets up exactly the question the Results answer: the four criteria (line 52) map onto the constraint tests in Results, and the hypothesis (line 48) is stated in the same terms as the findings.
3. Methods gives a reader enough to follow Results: every Results subsection has a matching Methods subsection except published signatures and exploratory arms (issue M8), and the verdict rules are stated before use.
4. The Discussion's literature handling is well ordered in its middle: the lung-entropy distinction (446), the endpoint-separation literature (448), and the dominant-literature-against section (450-452) build the concentration claim's context in a logical sequence.

## 6. Counts

1. BLOCKER: 1 issue (B1).
2. MAJOR: 13 issues (M1-M13).
3. MINOR: 12 issues (m1-m12).
4. Total: 26 issues.
