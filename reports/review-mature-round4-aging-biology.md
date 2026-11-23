# Round 4 review: aging biology, theory, and literature

Manuscript of record: commit dd6f503, tag mature-to-old-donor-scatterplots-2026-11-22.
File reviewed: report/plos-one/final-plos-one.tex, all 4 main figures, all 6 SI figures.
Reviewer role: hostile aging-biology, theory, and literature reviewer.
Date: 2026-11-22.

## Verdict: MAJOR REVISION

The data work is real and several headline numbers check out against the repo files.
The narrowing pattern in Yazar blood is real. The fraction scan is real.
The CD4 candidate numbers trace exactly to the mature Block 1 files.
But the paper hides three sensitivity failures that bear directly on its headline claims,
leads the Abstract with a hypothesis stated as fact, and carries figure label errors
that misstate what was measured. A reviewer who finds the Gao raw means or the
Yazar 3000-UMI rarefaction result first will not trust the rest. Fix the five
blockers, answer the majors, then resubmit.

Counts: 5 BLOCKER, 14 MAJOR, 9 MINOR.

## Method note

Every literature claim below was checked against the cited work or the repo file
named in the item. Claims marked CONFIRMED were verified by web search or file
read on 2026-11-22. Where the tex line number is given, it refers to the dd6f503
tex as read in this session.

## BLOCKER items

1. Gao scatter: the raw means are flat, and the paper never says so.
   Tex lines 249 to 250, 318 to 319, 402. The adjusted slope is minus 0.0723
   standard deviations per decade, p 0.041. The raw donor-by-cell-type means
   from results/followup-mature-fulldepth-gao-donor-part1.csv are YOUNG 30 to 39
   0.6029 (n 100), MID 0.6056, OLD 0.6048; donor-level young 0.6090 (n 12)
   versus old 0.6083 (n 41). The fall exists only after sex, depth, and batch
   adjustment. The tex presents the adjusted slope as the finding and uses it
   in the Discussion and Conclusions as one half of the narrowing evidence.
   The flat raw means appear nowhere in the tex or in
   reports/followup-mature-fulldepth.md. A paper whose headline is "cells get
   narrower, not messier" cannot hide that its Gao scatter means do not move.
   Fix: state the raw means, state that the Gao effect is adjustment-dependent,
   and downgrade the Gao scatter claim to a covariate-adjusted association.

2. Yazar Shannon entropy loses significance at 3000-UMI rarefaction, and the
   paper never says so. Tex lines 251 to 252, 318, 402. The stratified-global
   primary at rarefied_3000 gives slope minus 0.0268, p 0.281, from
   results/followup-mature-fulldepth-yazar-spread.csv. The value is listed in
   reports/followup-mature-fulldepth.md line 113 but never reaches the tex.
   The tex says only "Rarefaction at 1,500 molecules leaves 150 of 159 rows
   unchanged." The Yazar entropy fall is one of two pillars of the narrowing
   claim, and it dies at the stricter rarefaction depth. Fix: report the
   3000-UMI result in Results and carry the caveat into the Discussion and
   Conclusions.

3. The lung headline contradicts the paper's own study-aware test. Tex lines
   19, 231, 321, 403. The Abstract, Discussion, and Conclusions lead with
   "five human fraction shifts hold". The study-aware meta-combination in
   reports/followup-mature-composition.md gives club cell fixed effect minus
   0.0061, p 0.58 (line 112) and multiciliated epithelial cell fixed effect
   plus 0.0098, p 0.39 (line 117), with no study-specific slope significant.
   The tex discloses this once, at line 231, then leads with the opposite
   message in three other places. Both cannot lead. Fix: the Abstract and
   Conclusions must state that the two lung shifts do not survive
   study-aware combination, or the paper must show why the study-aware test
   is the wrong test. As written, the paper buries its own null.

4. The Abstract states the hypothesis as established fact and upgrades the
   interneuron finding. Tex lines 19 to 20 and 25 to 26. "Aging is local and
   plural" is the hypothesis under test, not a result, yet the Abstract opens
   with it as a finding. "Interneuron expression replicated in two cohorts"
   omits that the moderated check is unverified and that SEA-AD failed to
   replicate. The shared gene is TRMT11 only, from an unmoderated test, at
   about 4 counts per 10,000, with no established aging biology. Fix: rewrite
   the Abstract so the hypothesis is framed as tested, not proven, and state
   the interneuron result with its failures attached.

5. Figure S1 mislabels what was counted. The figure and its caption say "one
   test group is one tissue-by-species group" (scripts/build_si_figures.py
   lines 106 to 107). The P1m genes file
   (results/followup-mature-p1m-genes.csv, columns n_units_significant and
   units_list, for example "Marrow::hematopoietic precursor cell") counts
   tissue-by-cell-type units. The counts verify as tissue-by-cell-type:
   1212, 341, 163, 66, 48, 27, 24, 18, 9, 5, 1 sum to 1914; 702 in 2 or more,
   361 in 3 or more, 132 in 5 or more. The label is wrong about the unit of
   analysis in the paper's own signature-support figure. Fix the label and
   caption to tissue-by-cell-type units.

## MAJOR items

1. Gao Shannon entropy is specification-dependent, and the tex reports only
   the primary specification. Tex lines 251 to 252, 321, 403. Primary
   full-depth minus 0.1001, p 0.0539; depth-by-cell-type minus 0.1013,
   p 0.0344; spline-residual minus 0.1325, p 0.0288; rarefied_1500 p 0.0648;
   rarefied_3000 p 0.0727. Results calls it null, but Discussion line 321
   and Conclusions line 403 use it as corroboration ("falls in the same
   direction"). A null that flips to significant across defensible
   specifications cannot serve as corroboration. Fix: report the range and
   stop using the Gao Shannon result as supporting evidence.

2. Effect sizes are never translated into anything a reader can judge. Tex
   lines 249 to 252; Figure 2. The fitted lines in Figure 2 are nearly flat
   to the eye: Yazar scatter minus 0.0197 and Yazar Shannon minus 0.0298
   standard deviations per decade, Gao scatter minus 0.0723. The tex never
   states what a 0.02 to 0.10 standard-deviation-per-decade shift means in
   absolute terms or in biological terms. A 0.02 shift per decade is about
   0.1 standard deviations across a 50-year span. The paper also never
   addresses the Gao scatter outlier clusters near minus 3 to minus 4.5 at
   ages about 35 to 42. Fix: translate every headline slope into absolute
   change across the studied age span and explain or remove the outliers.

3. Figure S2 plots mouse per-month and human per-year slopes on one axis.
   The mouse limb-muscle macrophage slope minus 0.143 per month equals
   minus 1.72 per year but is plotted at minus 0.143. The footnote discloses
   the units, but the axis is uninterpretable as a single scale and the
   mouse rows look smaller than they are. Fix: convert to a common unit or
   split the figure by species.

4. Figure S5 mislabels its null. The legend and title say "random-gene
   null", but the null is the joint donor-permutation null (one permuted age
   per donor, shared across that donor's cell types; Methods). The error bars
   on the null bars are unlabeled. The nominal p 0.022 is emphasized while
   the q 0.066 appears only in the footnote. The primary three-test family
   gives q 0.066, which fails the 0.05 bar the paper uses elsewhere. Fix:
   relabel the null, label the error bars, and lead with q, not p.

5. The Oh et al. numbers in the Introduction do not match the cited paper.
   Tex lines 44 to 45. The abstract of Oh et al., Nature 2023, 624:164-172,
   reports 5,676 subjects, nearly 20 percent with strongly accelerated age
   in one organ, 1.7 percent multi-organ agers, and a 250 percent increased
   heart-failure risk for accelerated heart aging. It does not report an
   "organ-gap correlation r 0.29" or "heart failure about 2.4-fold per
   standard deviation". The 250 percent figure refers to extreme agers, not
   a per-standard-deviation hazard. Fix: verify both numbers against the
   full text and correct or delete them. CONFIRMED for the 5,676, 20 percent,
   and 1.7 percent figures; NOT VERIFIED for r 0.29 and 2.4-fold per SD.

6. No pathway-level uniformity test was run, yet the paper concludes against
   a shared program. Tex lines 53 to 56, 313 to 314, 320 to 322, 403 to 404.
   The four-criterion construct is the paper's own test, and criterion 3
   ("comparable magnitude after descriptive scaling by species lifespan") is
   descriptive scaling, not a test. No gene-set or pathway test of uniformity
   was run; the Dobner pan-endothelial gene-set test and the pathway-level
   uniformity test remain on the unresolved list. No mechanism was measured
   anywhere in the paper. Fix: downgrade the conclusion to "no shared
   signature was found at the gene level under these four criteria" and
   state plainly that pathway-level uniformity was not tested.

7. The cross-species headline rests on an unacknowledged design asymmetry.
   Tex lines 41 to 42, 313 to 314. Mouse mature is 18 months; human mature
   is 30 years. These are different life stages, and the Methods disclosure
   does not reach the headline. Cell types that share a name across mouse
   and human are not proven evolutionary matches, which the paper states in
   Methods but not where the cross-species conclusion is drawn. Fix: carry
   both caveats into the Results and Conclusions wherever a cross-species
   claim is made.

8. Typo: tex line 53 reads "so the P1 analysis neither confirms nor refutes
   that distinct claim". The context is Tyshkovskiy clocks and the mature
   transport analysis, so it should read P1m. Fix the typo.

9. Data Availability links to a stale tag. Tex lines 407, 411, 413 link to
   tag mature-to-old-rewrite-2026-11-21. The manuscript of record is dd6f503,
   tag mature-to-old-donor-scatterplots-2026-11-22. Fix all three links.

10. The paper cites Long et al. as a contradiction but never reconciles the
    two measures. Tex lines 302 to 303. Long et al., Cell Reports 2025,
    44:116612, report that individual-level noise rises with age, consistent
    across six cell types except dendritic cells, in the same OneK1K/Yazar
    blood the paper uses. CONFIRMED from the Long full text. The paper's
    Yazar scatter falls. Long measures gene-level CV2 noise aggregated across
    genes; the paper measures cell-level Hellinger scatter between cells.
    These may measure different things, or one may be wrong. The tex never
    says which. Fix: explain the measure difference and state whether the
    two results conflict or are compatible.

11. The "first view" the paper argues against is a strawman. Tex lines 47
    to 50, 313 to 314. The paper defines a first view in which one upstream
    failure drives the rest, so old cells of every type show the same broad
    pattern, then refutes it. No cited theory holds that view. The hallmarks
    framework, antagonistic pleiotropy, disposable soma, and hyperfunction
    are all compatible with heterogeneity, which the paper itself concedes
    at lines 313 to 314. Fix: drop the first-view framing or name a specific
    author who holds it. The theory section is otherwise fair: the hallmarks
    description, the pleiotropy and soma characterizations, and the
    hyperfunction description match the cited works.

12. The Gao E2 singleton-fraction metric rises with age and is never
    mentioned. Tex lines 253 to 254; reports/followup-mature-fulldepth.md
    line 127. Gao E2_stratum singleton_fraction_20cells gives slope 0.151,
    p 0.0116. The tex reports only the E2 stratum count, which falls. The
    singleton fraction is a distinct E2 metric, and its rise points the
    other way. Fix: report it and explain why the count metric governs the
    verdict, or drop the claim that all E2 evidence points one way.

13. The interneuron result is one unmoderated gene, and the paper still
    leads with it. Tex lines 25 to 26, 271 to 272, 282 to 283. TRMT11 is the
    only shared gene, from an unmoderated test, slopes 0.00205 to 0.00406
    log counts per 10,000 per year (about 2 to 4 percent per decade),
    expressed at about 4 counts per 10,000, with no established aging
    biology. The moderated mature intersection is unverified, and SEA-AD
    failed to replicate. The genome-wide correlations (VIP r 0.7523,
    parvalbumin 0.6299, somatostatin 0.6996) are real, but a correlation of
    effect sizes is not replication of a gene. Fix: demote the interneuron
    section to a replication-of-correlation result and remove TRMT11 from
    the Abstract.

14. The Abstract calls the MSSM depth decline a "candidate" alongside the
    CD4 result. Tex lines 19 to 20, 297 to 299. The paper itself states the
    depth decline is a technical observation, not a biological output:
    captured depth falls with age, sex-only adjustment, postmortem interval,
    brain bank, RNA quality, processing batch, and technical principal
    components all unavailable. A technical artifact is not a candidate
    aging signal. Fix: remove it from the Abstract's candidate sentence and
    keep it only as a caution about the brain data.

## MINOR items

1. Figure 4 panel-title capitalization is inconsistent ("VIP", "pvalb",
   "sst"). The SST y-axis runs to 0.012 while VIP and parvalbumin run to
   0.008, which visually inflates the SST panel. Fix both.

2. The Abstract names the moderated mature intersection (TRMT11) and then
   disclaims it as unverified in the same sentence. Tex lines 25 to 26.
   Do not name unverified results in the Abstract. Delete the clause.

3. The Tabula Sapiens 2.0 claim overstates the citation. Tex lines 290 to
   291. The paper says the senescence analysis "finds the same duality as
   this paper: senescent cells share common gene-expression modules
   alongside organ- and cell-specific patterns". The Tabula Sapiens 2.0
   abstract (Cell 2026, 189:7212-7227, DOI 10.1016/j.cell.2026.08.010)
   emphasizes heterogeneity and context-dependent senescence and says
   CDKN2A-positive phenotypes differ markedly by cell type and tissue
   context. The "common modules" half is not in the abstract. Fix: cite the
   specific figure or module, or soften the claim to heterogeneity only.
   Publication status CONFIRMED.

4. The Abstract's "five of fifteen fraction changes held" carries no family
   context. Tex lines 18 to 19. The 624-family total is 45 passes (17 core
   plus 28 exploratory-species), and the two promoted core tests are not
   established. A reader of the Abstract cannot tell what family the five
   survived. Fix: add one clause naming the family.

5. Figure 2 Gao scatter panel outlier clusters near minus 3 to minus 4.5
   at ages about 35 to 42 are unexplained in the tex or the caption. Fix:
   explain them or show they do not drive the slope.

6. The Abstract states "no uniform transcriptional response" without
   defining the four-criterion construct the claim rests on. Tex lines 19
   to 20. A reader cannot evaluate a verdict whose test is undefined in the
   Abstract. Fix: name the criteria briefly or soften the sentence.

7. The Jeffries citation does not note the age bins. Tex line 302. Jeffries
   et al., Nature 2025, 646:657-666, compare 7 elderly donors (70 plus)
   with 10 adult donors (15 to 49), so the "adult" group includes a
   15-year-old. The tex cites the IN-SST variability rise without noting
   this. Fix: add the bins. The 124-gene and 10-in-all-13-types figures are
   CONFIRMED verbatim from the Jeffries full text, as is the statement that
   no neuron subtype is lost in non-pathological aging.

8. Figure S6 uses E1 and E2 labels without spelling them out in the figure.
   This violates the paper's own no-invented-jargon rule, which requires
   every figure to stand alone with no internal abbreviations. Fix: spell
   out off-program count and stratum count in the figure and caption.

9. The tex never states whether De Man's entropy is computed the same way
   as the paper's Shannon entropy. Tex lines 300 to 301. De Man et al., Nat
   Commun 2026, 17:2095, report that transcriptional entropy rises with age
   in human lung and is an independent predictor of age. CONFIRMED. If De
   Man's entropy is a different quantity from the paper's per-cell Shannon
   entropy, the direct contradiction the paper concedes may be a measure
   mismatch instead. Fix: state the two definitions side by side.

## Literature checks that passed

1. Tyshkovskiy et al. [11]: CONFIRMED. Nature 2026, 654:173-188. CDKN1A
   and LGALS3 are named as shared features across four mammals; top shared
   pro-mortality changes include upregulation of Cdkn1a, Lgals3, Casp1,
   S100a8, and S100a4.
2. Long et al. [67]: CONFIRMED. Cell Reports 2025, 44:116612. Noise rises
   with age in blood, consistent across six cell types except dendritic
   cells, in 1.23 million cells from 981 individuals.
3. Hwang et al. [36]: CONFIRMED. Aging 2026, 18:575-592. Mouse 2,760 genes
   and human 613 genes change significantly; under 5 percent of enriched
   pathways are shared.
4. Buckley et al. [45]: CONFIRMED. Pan-mammalian clocks from 11,754
   samples, 59 tissue types, 185 mammalian species.
5. Sloan et al. [46]: CONFIRMED. Cell Genomics 2025. Loss of SST transcript
   and SST-positive neuron density; shared loss of oxidative-phosphorylation
   and cytoplasmic-translation genes across cortical types.
6. Kivimaki et al. [49]: CONFIRMED. Lancet Digit Health 2025, 7:e195-e204.
   6,235 participants, 30 diseases, 6 exclusively linked to their own
   organ's age gap.
7. Peters et al. [34]: CONFIRMED. Nat Commun 2015, 6:8570. Whole-blood
   meta-analysis in 14,983 individuals, 1,497 genes differentially
   expressed with age.
8. Petkovich et al. [35]: CONFIRMED. Cell Metab 2017, 25:954-960. Mouse
   blood clock from 90 CpG sites.
9. De Man et al. [41]: CONFIRMED. Nat Commun 2026, 17:2095. 32 human lungs
   plus 28 IPF control lungs, 199,400 cells, 60 donors, GEO GSE281219.
   General capillary 21.5 percent, alveolar type II 19.5 percent
   differentially expressed genes; transcriptional entropy rises with age.
10. Dobner et al. [42]: CONFIRMED. Nat Aging 2026, 6:2530-2547. Human Aging
    Endothelial Cell Atlas, 375,829 endothelial cells, 290 donors, 12
    tissues, 34 studies; rising inflammation, declining angiogenesis,
    extracellular-matrix remodeling, p21-linked senescence program.
11. Tabula Sapiens 2.0 [69]: CONFIRMED as published, Cell 2026,
    189:7212-7227, DOI 10.1016/j.cell.2026.08.010, PMID 42685689.
12. Jeffries et al. [40]: CONFIRMED in full. Nature 2025, 646:657-666.
    124 genes downregulated across multiple cell types; HSPA8, TUBA1A, and
    eight other genes downregulated in all 13 brain-cell types; IN-SST the
    one type with a significant rise in coefficient of variation, p 0.043;
    no neuron subtype lost in non-pathological aging.
13. CD4 candidate numbers: CONFIRMED against
    results/followup-mature-block1-residual-dispersion.csv and
    reports/followup-mature-block1.md lines 74 and 228 to 230. Primary
    0.00071 per decade, p 5.174e-08, q 8.279e-07, n 777; reweighted 0.00080,
    p 1.165e-06; unadjusted 0.00003, p 0.8142; rarefied 0.00008, p 0.7971;
    debiased 0.000313, p 0.000638 after 56.4 percent finite-cell calibration
    bias; depth-adjusted 0.000295, p 0.00106. The residual log cell-fraction
    coupling (beta minus 0.00111, p 4.48e-07) and the 12-of-30 imbalanced
    states after weighting are also on file.
14. Lung study-aware numbers: CONFIRMED against
    reports/followup-mature-composition.md lines 110 to 113 and 117. Club
    meta fixed effect minus 0.0061, p 0.58; multiciliated meta fixed effect
    plus 0.0098, p 0.39.
15. S3 Gao within-donor slope: CONFIRMED. Minus 0.062, p 0.085, in
    scripts/build_si_figures.py line 324.

## What would change the verdict

Fixing blockers 1 to 3 alone would move this to a minor-revision verdict.
Blocker 4 and the majors are framing and disclosure, not new data. No new
experiment is needed for any item above. The one analysis that would
strengthen the paper beyond this review is the pathway-level uniformity
test, which the paper has not run and which its conclusion assumes.
