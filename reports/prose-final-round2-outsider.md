# Prose final round 2: outsider-reader review

## Reviewer, base, and scope

1. Reviewer: a scientist outside aging and single-cell work, reading cold as a structural biologist / systems neuroscientist would.
2. Base: `~/workspace/agescape-repo`, `report/plos-one/final-plos-one.tex`, base commit `d15befd` (tag `mature-to-old-final-base-2026-11-23`), 815 lines read end to end.
3. Scope: prose only. No science was re-litigated. No number, verdict, caveat, disclosure, or citation was challenged. No new analysis was proposed. The 624-test family, the marrow framing, the held-out and Bayesian results, and the literature framing were treated as settled.
4. The abstract was treated as fixed. Its content was not reviewed and no change to it is proposed, except for the one direct abstract-body tension flagged in the Abstract check below.
5. Severity definitions used here: BLOCKER means a reader outside the field cannot follow the argument at that point; MAJOR means a passage is opaque without background the paper does not supply in place; MINOR means a term deserves one clarifying phrase.
6. Every issue below gives the exact tex line number, the problem in one sentence, and the exact fix recommended.
7. Tex lines are long paragraphs, so a line number points to a whole paragraph; the quoted phrase inside it locates the exact spot.

## Abstract check (fixed abstract; flag only)

1. Line 34 says: "Only brain immune cells called microglia passed both tests." Lines 332-334 say the organ sweep gives three positives to count: GTEx whole blood as a bulk-composition candidate, retinal pigment epithelium in the negative direction, and Lake kidney endothelium under the unfiltered signature only. Fix: reconcile the abstract sentence with lines 332-334 in the body, for example by scoping the "only microglia" claim to the core cohorts named at lines 320-322, and do not change the abstract wording in this review beyond flagging the tension.
2. Line 34 also says spread fell in heart and blood. Line 266 says Shannon entropy falls in blood (minus 0.0298, p = 0.0036) and is null in heart (minus 0.100, p = 0.054). Fix: in the body at line 266, state in one sentence whether "spread" in the abstract means within-cell scatter, which falls in both cohorts at lines 264-266, or Shannon entropy, which does not fall in heart, so an outsider does not read lines 34 and 266 as contradicting each other.

## BLOCKER issues

### B1. Lines 58, 66, and 238: three roadmaps disagree with the actual order of the paper
1. Problem: line 58, line 66 ("tests in spine order"), and line 238 each promise a different order, and none matches the actual Methods order (signature derivation at line 120 before composition at line 138, held-out validation at line 130 before the tests it validates) or the actual Results order (interneuron replication at line 296 and signature derivation at line 306 appear before transport at line 316, held-out validation at line 372 appears after transport, drift, and concordance).
2. Fix: rewrite lines 58, 66, and 238 to the same single list in the actual order used, namely fractions, concentration, limited signals, interneuron replication, signature derivation, transport, drift, concordance, published signatures, held-out validation, Bayesian refits, disorder tests, exploratory arms, follow-up cohorts, and delete the phrase "spine order" at line 66.

### B2. Lines 320-322 and 136: "pass both tests" is used as the central verdict before it is defined in plain language
1. Problem: lines 320-322 report that only two units "pass both tests" without saying in that paragraph what the two tests are, and the distinction between the composite test and the 2,000-draw matched-random test is buried in line 136 and only clarified as "uncalibrated screening" versus "calibrated evidence" at line 330, after many uses.
2. Fix: insert two sentences at the start of the transport Results at line 318: "A unit passes only if it clears two separate checks. First, the signature score must rise or fall with age in a standard regression (the composite test). Second, it must beat 2,000 random gene sets matched for size and expression level (the matched-random test); only the second is calibrated evidence."

### B3. Lines 244-252 and 246 (Figure 1 caption): the fraction counts cannot be reconstructed
1. Problem: lines 244-252 juggle 399 rows, 15 full-range associations, 17 core passes, 10 rows clearing the bar, 8 rows under the 431-test family, and 14 under within-file correction, so an outsider cannot say in one sentence how many fraction associations exist or which five are the headline.
2. Fix: rewrite lines 244-252 as a counting chain in four short sentences: "The mature scan has 399 rows. Fifteen rows passed in the full-range scan. Five of those fifteen still pass in the mature scan and are the headline. Five further rows pass in the same mature scan, for ten in total; the 8 and 14 are sensitivity counts under other corrections," and repeat that chain verbatim in the Figure 1 caption at line 246.

### B4. Lines 262-268: the 14 concentration endpoints are listed but never mapped to a question
1. Problem: line 262 lists Part 1 center shift, Part 1 scatter, Part 2 Shannon entropy, Part 2 Kullback-Leibler methods, Part 2 Kullback-Leibler square-root, Part 3 primary off-program count, and Part 3 stricter off-program count, then says they are "one pattern in several guises", so an outsider cannot say which endpoint answers which question or which one is primary.
2. Fix: replace the list at line 262 with a small table with three columns (endpoint, plain question it answers, primary or secondary), and add one sentence: "Scatter is the headline endpoint; Shannon entropy, Kullback-Leibler divergence, and off-program counts are the same concentration pattern measured in different ways, not four independent findings."

### B5. Lines 316-344: transport denominators shift without a map
1. Problem: the transport section moves between 252 units, 170 testable units, a 48-unit core comparison, a 79-unit expanded comparison, 144 retina units, 47 GTEx tissues, and 35 PBMC cell types, and the order statement at line 318 (core, rescue, reverse, organ sweep, PsychAD, PBMC) does not match the actual interleaving of PsychAD microglia inside the core paragraph at lines 320-322 and GTEx, retina, and retinal pigment epithelium inside the organ sweep.
2. Fix: insert a transport map table before line 320 with one row per block (core 48, expanded 79, reverse, heart, blood, retina, retinal pigment epithelium, GTEx Track 1, GTEx Track 2, PsychAD, PBMC), giving the denominator, the number passing both tests, and the direction for each, and rewrite line 318 to match the actual order.

### B6. Lines 380-392: the Bayesian section uses five unexplained technical terms to state its main qualification
1. Problem: lines 380-392 report posterior medians, 94 percent highest-density intervals, 95 percent credible intervals, region probabilities, a region of practical equivalence, and a Normal(0, 0.5) prior, without a plain sentence saying what a Bayesian refit adds or why two different interval widths appear, so an outsider cannot say what the section established beyond "agreed on direction".
2. Fix: insert a lead paragraph at line 382: "A Bayesian refit re-estimates each headline slope while stating a prior range for plausible slopes; the interval is the range holding most of the plausible values, the direction probability is the chance the true slope has the stated sign, and the region of practical equivalence is a pre-set band around zero treated as too small to matter. The 94 percent and 95 percent intervals are the same idea computed by two fitting methods," and define "conjugate refit" in that paragraph as a faster approximate refit of the same model.

### B7. Lines 302 and 322: Figure 4 appears before Figure 3
1. Problem: the interneuron figure is Figure 4 at line 302 and the transport figure is Figure 3 at line 322, so a reader following figure numbers in order meets Figure 4 first and cannot tell whether a figure is missing.
2. Fix: renumber the interneuron figure at line 302 to Figure 3 and the transport figure at line 322 to Figure 4, and update every in-text reference to Figure 3 and Figure 4 to match.

### B8. Lines 246 and 270: the two main figure captions do not stand alone
1. Problem: the Figure 1 caption at line 246 uses "the bar", "the refit family", and "five further rows" without defining any of them, and the Figure 2 caption at line 270 reports that raw clouds run the opposite way from adjusted slopes and that three Sim 2021 donors form a low cluster, without saying whether the figure shows raw or adjusted values.
2. Fix: rewrite both captions to define every term they use in place: in Figure 1, replace "the bar" with "q below 0.05 under the refit family on the 399-row mature scan"; in Figure 2, start with "Each point is a partial residual after removing model covariates, not a raw value; the line is the adjusted slope," and add one sentence naming the Sim 2021 donors as three donors from one 2021 study whose low values were not tested for influence on the slope.

## MAJOR issues

### M1. Line 42: three evolutionary theories are named with glosses too thin to use
1. Problem: antagonistic pleiotropy, disposable soma theory, and hyperfunction theory each get one clause, and mechanistic target of rapamycin (mTOR) is expanded but not explained, so an outsider cannot follow the Discussion return to these theories at line 462.
2. Fix: add one plain sentence per theory at line 42, for example: "Antagonistic pleiotropy means a gene that helps early in life harms late in life; disposable soma means the body spends limited energy on reproduction instead of repair; hyperfunction means growth programs keep running after development ends, including signals through mTOR, a protein that senses nutrients and drives cell growth."

### M2. Line 44: the DNA methylation clock is expanded but not explained
1. Problem: "90 cytosine-phosphate-guanine (CpG) sites" expands the abbreviation without saying what a methylation clock measures, so an outsider cannot tell why distinct mouse and human sites matter.
2. Fix: rewrite the sentence at line 44 as: "A mouse DNA methylation clock, a predictor built from chemical marks added to DNA at 90 CpG sites (places where a cytosine base sits next to a guanine base), detects calorie restriction in mice, but the sites that carry the signal differ between mouse and human."

### M3. Line 118: the mouse processing paragraph is a software and statistics list
1. Problem: line 118 uses DESeq2, SCANPY, log1p, counts per million scaling, t-based regression inference, and Benjamini-Hochberg false discovery rate correction in one paragraph with no plain gloss for any of them.
2. Fix: add one parenthetical gloss per term at line 118: "DESeq2 (a standard tool for testing gene-level change), SCANPY (a standard tool for handling single-cell data), log1p (a log scale that keeps zero counts at zero), counts per million (counts rescaled so each cell totals one million), t-based inference (a standard test that reports a p value), and Benjamini-Hochberg correction (an adjustment that limits false positives when many genes are tested)."

### M4. Lines 122-124: the signature derivation paragraph assumes five genetics terms
1. Problem: lines 122-124 use up-regulated, down-regulated, absolute log2 fold change of at least 0.25, one-to-one mouse-to-human orthologs, and Ensembl derived mapping without defining any of them, and "marrow" is never said to mean bone marrow.
2. Fix: rewrite lines 122-124 with in-place glosses: "up-regulated means higher in old mice, down-regulated means lower; log2 fold change is the size of that change on a doubling scale, so 0.25 is a small change; orthologs are the matching gene in the other species; Ensembl is the public gene catalogue that supplies the match; marrow means bone marrow, the blood-forming tissue inside bones."

### M5. Lines 140-142: the centered log ratio formula has no plain reading
1. Problem: lines 140-142 define a cell fraction as a centered log ratio slope and then give the formula CLR_ij with a pseudocount a, without saying in words what the number means or why a pseudocount is added.
2. Fix: insert before the formula at line 142: "A centered log ratio slope is the yearly change in one cell type's share relative to the average share of all types in that organ; a pseudocount of 0.5 is a small number added to every count so a type with zero observed cells does not break the log calculation."

### M6. Lines 152-154: off-program expression and Shannon entropy are defined by thresholds and a formula only
1. Problem: lines 152-154 define off-program genes by detection in at most 5 percent of mature target cells, at least 20 percent of another type, at least 20 counts per 10,000, and at least fivefold enrichment, and define Shannon entropy as minus the sum of p_g log(p_g), without saying why those cutoffs were chosen or what the entropy number means in words.
2. Fix: add two sentences at line 152: "Off-program means a cell turns on genes that belong to a different cell type; the cutoffs pick genes that are rare in the target type and common elsewhere, so a single stray read does not count. Shannon entropy is a single number for how evenly a cell spreads its RNA across genes: high means spread across many genes, low means focused on few genes."

### M7. Lines 162-168: four targeted tests are compressed into three paragraphs with four different verdict rules
1. Problem: lines 162-168 cover residual dispersion, gene coordination, categorical signature shape, and stress-score stratification, each with its own pass rule (finite-cell calibration, split-half noise floor, early-transition label, 90th percentile stress cells), so an outsider cannot say what any one test would have shown if it had passed.
2. Fix: split lines 162-168 into four short sub-paragraphs, each starting with one plain sentence: "Residual dispersion asks whether single genes vary more between cells after removing the average trend; gene coordination asks whether genes still rise and fall together; the categorical test asks whether change happens early or late; the stress test asks whether removing the most stressed cells erases the age signal."

### M8. Lines 158, 164, 186, and 262: four near-synonyms for distance from a reference are never disambiguated
1. Problem: drift, vector displacement, center shift, and scatter are used across lines 158, 164, 186, and 262 for related but different distance measures, and the paper never states in one place how they differ.
2. Fix: insert a four-row glossary box after line 186: "Drift is the whole-profile distance from a mature reference; center shift is the same distance averaged per donor; scatter is the distance between cells of the same type from the same donor; profile-cloud dispersion is the distance between donors."

### M9. Lines 196-216 (Table 3): the verdict-rules table cannot be applied by an outsider
1. Problem: Table 3 uses "BH-significant primary", "smallest detectable effect at most 0.50 SD", "split-half floor", "unresolved", and "powered non-detection" in every row without a legend or a worked example.
2. Fix: add a legend above Table 3 at line 196: "BH-significant means the result passes after correction for testing many types at once; smallest detectable effect is the smallest true effect the test had an 80 percent chance to catch, in standard deviation units; split-half floor is the smallest slope that cannot be explained by splitting each cell's molecules at random; unresolved means the test could neither pass nor rule out an effect."

### M10. Lines 250-252: study-aware combination and the propeller check are named but not explained at the point of use
1. Problem: lines 250-252 report study-aware fixed-effects models, meta-combination, and a propeller check without saying what each does, so an outsider cannot tell why the lung shifts fail them.
2. Fix: rewrite lines 250-252 with glosses: "Study-aware models estimate the slope inside each study separately and then combine them, so a shift driven by one study's sampling does not pass; the propeller check is a second, donor-level count model that tests the same share change in a different way."

### M11. Line 270 (Figure 2 caption): the caption reports adjusted slopes over raw clouds that run the other way
1. Problem: the Figure 2 caption at line 270 says raw Gao scatter rises with age while the adjusted slope falls, and raw Yazar pooled center shift falls while the adjusted slope rises, without explaining that the points are partial residuals, so an outsider reading the figure sees lines that contradict the clouds.
2. Fix: add to the caption at line 270: "Points are partial residuals after removing sex, batch, and depth, so the line is the adjusted slope and the raw cloud alone does not show it; read the slope, not the cloud direction."

### M12. Lines 286-292: the limited-signals section mixes three different kinds of result in one subsection
1. Problem: lines 286-292 present a surviving candidate (naive CD4 dispersion), a technical observation (captured library depth), and an inconclusive signal (naive CD8 scatter) in three paragraphs that all use slopes, q values, and sensitivity arms, so an outsider cannot say after reading which of the three is a biological claim.
2. Fix: add a header sentence to each paragraph at lines 286, 288, and 290 respectively: "This is a biological candidate with an unresolved mechanism"; "This is a technical measurement problem, not a biological claim"; "This is an inconclusive signal, not a finding."

### M13. Lines 300-302: the interneuron subtypes are named only by their marker proteins
1. Problem: lines 300-302 report VIP, parvalbumin, and somatostatin interneurons, a broad GABAergic comparator, TRMT11, hypergeometric expectations, and coexpression-block permutation without saying what any interneuron subtype is or what GABAergic means.
2. Fix: rewrite line 300 as: "VIP, parvalbumin, and somatostatin are three kinds of inhibitory brain cell, each named after a protein it makes; GABAergic means inhibitory, using the signal GABA; TRMT11 is a gene that adds a chemical mark to transfer RNA, with no known aging role; the hypergeometric expectation is the overlap expected by chance, and coexpression-block permutation shuffles donor labels while keeping groups of correlated genes together."

### M14. Lines 336-340: the PsychAD pooled rise flips sign between count scales without a plain explanation
1. Problem: lines 336-340 say ten positive slopes hold only on the log1p counts-per-million scale while nine have negative raw-count slopes, and the somatostatin unit is a direct sign conflict, without explaining in words why changing the count scale can flip a sign.
2. Fix: insert one sentence at line 338: "Counts per million rescales each cell to the same total, so a gene can rise as a share of the cell even while its raw count falls when the cell captures fewer molecules overall; the two scales answer different questions, and only microglia agree on both."

### M15. Lines 356-362: the concordance result is stated as a correlation against a null the reader cannot picture
1. Problem: lines 356-362 report observed mean Pearson 0.2105 against a null mean of 0.1465, an excess of 0.064 smaller than the smallest detectable effect of 0.068, and a null mean that is 70 percent of the observed correlation, without saying what the null permutes or why donor sharing inflates the correlation.
2. Fix: rewrite lines 358-360 as: "The test asks whether gene-level age effects point the same way in different cell types from the same donor. The null shuffles each donor's age label once and reuses it for all that donor's cell types, so cell types from the same donor still share a donor; most of the observed correlation comes from that sharing, not from age."

### M16. Lines 394-408: the disorder section reports three verdicts in three vocabularies
1. Problem: lines 394-408 cover residual dispersion ("does not prove", "inconclusive"), gene coordination ("zero prove, six do not prove, ten are inconclusive"), and survivor bias ("Part A", "Part B", "does not prove as a mechanism"), so an outsider cannot say in one sentence what the section established.
2. Fix: add a roadmap sentence at line 396: "This section tests three separate explanations for the concentration pattern: cells growing noisier gene by gene, genes losing coordination, and old samples looking narrower because fragile cells died; none of the three is supported in general, and the one surviving candidate is naive CD4 dispersion, reported above."

## MINOR issues

### N1. Line 38: "transcriptional program" and "transcriptomic" are never defined
1. Problem: line 38 uses "transcriptional program" as the paper's central phrase without saying that transcription means copying DNA into RNA messages.
2. Fix: add at line 38: "Transcription means copying DNA into RNA messages; a transcriptional program is the set of RNA messages a cell makes."

### N2. Line 54: the scope sentence lists six unmeasured things without glosses
1. Problem: line 54 lists chromatin, DNA damage, transposon RNA, protein, cell function, and actual cell turnover without a phrase saying what each is.
2. Fix: rewrite line 54 as: "Chromatin is the packaging of DNA, transposon RNA is RNA copied from jumping DNA elements, and cell turnover is how often cells are replaced; none of these was measured."

### N3. Lines 70-72: eleven 36-character dataset identifiers block the Datasets paragraph
1. Problem: lines 70-72 paste eleven UUIDs inline, so an outsider loses the cohort list inside the identifiers.
2. Fix: move the eleven identifiers at line 70 to Table 2 and keep only the cohort names in the text at line 70.

### N4. Line 74: the UMI definition is circular
1. Problem: line 74 defines depth as "total number of RNA molecules captured per cell, summarized as mean log unique molecular identifiers (UMIs)" without saying that a UMI is a barcode tag attached to each captured molecule so the same molecule is not counted twice.
2. Fix: rewrite line 74 as: "Depth is the number of RNA molecules captured per cell, counted by unique molecular identifiers (UMIs), barcode tags attached to each molecule so each is counted once."

### N5. Line 82 (Table 1): "doublet QC" is unexplained
1. Problem: Table 1 at line 82 lists "gene, UMI, and doublet QC" without saying a doublet is two cells captured together as if they were one.
2. Fix: change line 82 to "gene, UMI, and doublet QC (removal of droplets holding two cells at once)."

### N6. Line 106 (Table 2): "10x multiome" is unexplained
1. Problem: Table 2 at line 106 lists "10x 3' v3; 10x multiome" without saying multiome measures RNA and chromatin access in the same cell.
2. Fix: change line 106 to "10x 3' v3; 10x multiome (RNA plus chromatin access in the same cell)."

### N7. Line 124: "marrow" is never expanded to bone marrow
1. Problem: line 124 calls the signature "marrow-dominated" without ever saying marrow means bone marrow.
2. Fix: change the first use at line 124 to "bone marrow (marrow), the blood-forming tissue inside bones."

### N8. Line 128: "positive controls" is unexplained
1. Problem: line 128 says "Classic aging-gene positive controls were not significant" without saying a positive control is a known aging gene expected to show up if the method works.
2. Fix: rewrite line 128 as: "Classic aging-gene positive controls, known aging genes expected to appear if the method works, were not significant in the top derivation unit."

### N9. Line 134: the fluorescence-activated cell sorting arm is named without saying why it matters
1. Problem: line 134 says a sorting arm "was not run" and the sorting datasets hold ages 3, 18, and 24 months only, without saying sorting physically separates cell types before sequencing and would have given a cleaner test.
2. Fix: add at line 134: "Cell sorting physically separates cell types before sequencing, so it would have tested the signature in purified cells; it could not be run because no sorted 30-month mice exist."

### N10. Line 168: apoptosis, unfolded-protein response, and Hallmark sets are unexplained
1. Problem: line 168 builds a 271-gene score from "the apoptosis and unfolded-protein-response Hallmark sets" without saying apoptosis is programmed cell death, the unfolded-protein response is a cell stress reply to misfolded proteins, and Hallmark sets are curated public gene lists.
2. Fix: rewrite line 168 as: "Stress-score stratification used a fixed 271-gene score from curated public gene lists (Hallmark sets) for apoptosis, programmed cell death, and the unfolded-protein response, a cell stress reply to misfolded proteins."

### N11. Line 176: Levene tests and bootstrap intervals are unexplained
1. Problem: line 176 compares variance with "median-centered Levene tests and bootstrap variance-ratio intervals from 2,000 donor resamples" without saying a Levene test compares spreads and a bootstrap interval comes from resampling donors.
2. Fix: rewrite line 176 as: "Between-donor variance was compared with a Levene test, a test for whether two groups spread differently, and bootstrap intervals, ranges built by resampling donors 2,000 times."

### N12. Line 182: pseudobulk, precision-weighted least squares, and empirical-Bayes shrinkage are unexplained
1. Problem: line 182 sums "donor pseudobulk counts", fits with "precision-weighted least squares", and runs a sensitivity with "empirical-Bayes shrinkage" without defining any term.
2. Fix: rewrite line 182 as: "Donor pseudobulk counts are counts summed across all cells of one type from one donor; precision-weighted least squares gives more weight to donors with more cells; empirical-Bayes shrinkage pulls noisy gene estimates toward the average."

### N13. Line 186: AnAge values are cited without saying what AnAge is
1. Problem: line 186 scales effects "using AnAge values of 4 years for mouse, 122.5 years for human, and 31 years for naked mole rat" without saying AnAge is a public database of maximum lifespans.
2. Fix: rewrite line 186 as: "Effects were scaled to 10 percent of species maximum lifespan using AnAge, a public database of maximum lifespans, with 4 years for mouse, 122.5 years for human, and 31 years for naked mole rat."

### N14. Line 190: log1p is used without a gloss
1. Problem: line 190 normalizes counts "to counts per 10,000, and log1p transformed" without saying log1p is a log scale that keeps zero counts at zero.
2. Fix: change line 190 to "normalized to counts per 10,000, and log1p transformed (a log scale that keeps zero counts at zero)."

### N15. Line 220: the Bayesian Methods paragraph names its prior and likelihood without glosses
1. Problem: line 220 uses "Student-t likelihoods, donor intercepts, and unit residual scales, with a Normal(0, 0.5) prior on the global age slope" without saying a likelihood is the assumed shape of the data, a prior is the range of slopes assumed plausible before seeing data, and Normal(0, 0.5) means slopes near zero are assumed most plausible.
2. Fix: rewrite line 220 as: "Bayesian refits assume a Student-t shape for the data (a shape that tolerates outliers), give each donor its own baseline, and start from a Normal(0, 0.5) prior, meaning slopes near zero are assumed most plausible before seeing data."

### N16. Line 224: a cubic B spline with three degrees of freedom is unexplained
1. Problem: line 224 tests nonlinear age effects with "a cubic B spline with three degrees of freedom" without saying a spline is a flexible curved fit and degrees of freedom set how bendy it is.
2. Fix: rewrite line 224 as: "Nonlinear age effects used a cubic B spline, a flexible curved fit, with three degrees of freedom, a setting that limits how bendy the curve can be."

### N17. Line 228: "correction family" is used before it is defined
1. Problem: line 228 says "Each analysis family was defined before correction" without saying a family is the set of tests corrected together, so the size of the family sets how strict the correction is.
2. Fix: add at line 228: "A correction family is the set of tests corrected together; a larger family makes each single test harder to pass."

### N18. Lines 244-248: cortical layer and projection codes are unexplained
1. Problem: lines 244-248 report L6b glutamatergic neurons, L5 extratelencephalic neurons, and layer 2/3 neurons without saying L6b and L5 are layers of the cortex (the outer brain sheet) and extratelencephalic means projecting outside the cortex.
2. Fix: add glosses at first use in lines 244-248: "L6b glutamatergic neurons (excitatory cells in the deepest cortex layer), L5 extratelencephalic neurons (cells in layer 5 that project outside the cortex)."

### N19. Line 264: CD14-positive monocytes and effector memory CD8 cells are unexplained
1. Problem: line 264 lists "CD14-positive monocyte scatter" and "effector memory CD8 center shift" without saying CD14 is a surface marker on monocytes (a white blood cell type) and effector memory CD8 cells are experienced killer T cells.
2. Fix: rewrite line 264 as: "CD14-positive monocytes (a white blood cell type marked by the surface protein CD14) and effector memory CD8 cells (experienced killer T cells)."

### N20. Lines 292 and 342: PBMC is never expanded, and Zhu and Thomson are never introduced
1. Problem: line 292 uses "Thomson PBMC" and line 342 uses "Zhu conventional dendritic cells" without expanding PBMC as peripheral blood mononuclear cells (white blood cells from blood) or saying Zhu and Thomson are the first authors of the two blood cohorts.
2. Fix: expand at line 292 to "Thomson peripheral blood mononuclear cells (PBMC, white blood cells isolated from blood)" and at line 342 to "Zhu (GSE213516) and Thomson (GSE214546), the first authors of the two blood cohorts."

## Counts

1. BLOCKER: 8 issues (B1 to B8).
2. MAJOR: 16 issues (M1 to M16).
3. MINOR: 20 issues (N1 to N20).
4. Total: 44 issues.
5. Abstract check: 1 direct tension flagged (line 34 versus lines 332-334) and 1 wording ambiguity flagged (line 34 versus line 266); neither is counted in the 44 above, and no change to the fixed abstract is proposed.

## What an outsider can and cannot say after reading

1. An outsider can say the headline in one sentence: from maturity to old age, a few human cell-type shares shift, and single cells focus their RNA on fewer genes instead of growing noisier.
2. An outsider cannot say, without the fixes above, how the fraction counts relate to each other, what "pass both tests" means, which concentration endpoint is primary, how transport denominators relate, or what the Bayesian refits add beyond direction agreement.
3. The three fixes that matter most are B3 (the fraction counting chain), B2 (the two-test definition at the start of transport), and B5 (the transport map table), because those three passages carry the paper's two headline claims and its central negative claim.
