# Mature-to-old rewrite of the Agescape manuscript

Date: 2026-11-21 (PT). Base: repo head 9644ade (Phase 2 complete, all results
committed). Previous manuscript state: optimized 0a022d8, 33 pages, abstract
287 words, 68 references, 4 figures. This rewrite supersedes it as the
submission candidate. Nothing was submitted anywhere.

## What changed

- `report/plos-one/final-plos-one.tex` rewritten as one coherent mature-to-old
  story, as if it was all one experiment. No process narration in the main
  text. Young-reference comparisons moved to Appendix A only.
- Sections replaced wholesale: Abstract, the signature-derivation Methods
  subsection, the whole Results section, the whole Discussion (including
  Predictions and Limitations), Conclusions, Data Availability. Appendix A
  rewritten fresh. Appendix B rows updated to mature 30-plus verdicts where
  rerun, with a PBMC row added.
- Figures: Fig 1 composition kept (caption updated to state the mature table
  is in the text). Fig 2 mouse drift moved to Appendix A with an honest
  young-reference caption. Fig 3 concordance regenerated from the mature CSVs
  (3 rows: mouse within, mouse cross, human within cohort). Fig 4 interneuron
  regenerated from the mature unmoderated gene table
  (`results/followup-mature-interneuron-genes.csv`).
- New figure builder: `scripts/build_mature_figures.py`. New figure data:
  `results/followup-mature-concordance-figure-data.csv`.
- Checker: `scripts/check_plos_manuscript.py` PASS. Abstract 249 words with
  p-values, q-values, and Bayesian results. Four figures. 68 cited references,
  cited set exactly equals listed set. Section order correct. No em or en
  dashes. No listed process residue.
- Compile: Tectonic, outside the repo. PDF builds, 32 pages, 877.74 KiB.
  Overfull hbox warnings are in the References list only (long DOI/URL
  entries), same as the prior version; no table overflow warnings.
- Fragments and splice script kept at
  `~/workspace/agescape-mature-rebuild/phase3-fragments/` on the VM (not in
  the repo).

## Issue-by-issue accounting (brief open questions 1 to 11, plus Item 9b)

1. P1m 14x size jump unexplained. Handled: Results signature subsection and
   Discussion/Limitations state the numbers (P1m 1,914 genes vs P1 139 genes;
   P0m 2,598) and state plainly that the cause is not established. No
   explanation is invented. Verified against
   `results/followup-mature-p1m-derivation.csv` by verifier B.
2. P1m sex confound. Handled: disclosed wherever P1m is introduced (Methods
   signature subsection, Results, Discussion, Limitations, Predictions). All
   four 30-month mice are male; the 18-month group is 2 female, 2 male. P1m
   is never presented as a clean age contrast.
3. Moderated interneuron intersection UNVERIFIED. Handled: only the
   unmoderated intersection is claimed (TRMT11 only, across VIP, parvalbumin,
   somatostatin). Methods, Results, and Limitations state the moderated check
   did not verify because a frozen recovery gate failed. The full-range
   six-gene moderated set (FBN1, GABRB1, IGF1R, PHF3, TRMT11, ZBTB16) is
   reported in Appendix A only, as context.
4. Mouse concordance p 0.0220 nominal. Handled: the q values are reported
   (0.0440 two-test family, 0.0460 six-row mouse family) and the p value is
   never given alone. Verified against
   `results/followup-mature-concordance-mouse.csv` (first row per scope,
   primary detection-filtered set): within observed 0.2105, null mean 0.1465,
   CI 0.1012 to 0.2055; cross observed 0.0604, null 0.0297, CI 0.0006 to
   0.0826; human observed 0.1793, null 0.2126, CI 0.1538 to 0.2861.
5. Propeller NOT_VERIFIED. Handled: Methods states the mature refit could not
   run the check because Rscript, speckle, and limma are unavailable in the
   analysis environment, and no mature propeller result is reported. Results
   composition and Limitations repeat that no propeller support is claimed.
6. Batch A runtime model_used label. Handled: the manuscript never cites the
   runtime label. Methods describes direct execution of the transport
   analyses. The label in the Batch A runtime CSV reads "Gemini 4 Neon
   (XHigh)" and is a hardcoded label, not the runner; it is not cited.
7. PsychAD IN_SST sign conflict. Handled: Results PsychAD paragraph discloses
   the scale dependence (10 positive slopes are log1p-CPM-only; 9 have
   negative raw-count slopes; IN_SST is a direct sign conflict, CPM adjusted
   +0.001201 q 0.000139 vs raw -0.001704 q 0.00940; only Micro holds on both
   scales). The pooled rise is not presented as scale-free.
8. PBMC 200-draw null miscalibration. Handled: no PBMC empirical-only pass is
   cited anywhere. PBMC transport rests on the composite test only (0 of 35
   both-pass). The miscalibration (26 of 35 rows at the 1/201 floor) is
   stated in Results.
9. Locked 624 numbers preserved. Handled: the 45-pass record, the two
   promoted tests (mouse lung vein endothelial, human kidney cortex
   collecting duct), and their sensitivity disclosures (core-only q,
   total-cell-adjusted q 0.0627 and 0.4231, propeller p 0.0528 and 0.0588)
   are preserved exactly as locked. The mature composition table is a
   separate declaration with its own correction families; its counts never
   enter the 624 family. Fly, worm, and zebrafish rows are never claimed as
   mature-restricted.
10. TS2 prose slip. RESOLVED by direct CSV check (verifier B also checked).
    The brief's framing conflated two cohorts. The Batch C report's "six" is
    the Lake kidney cohort: 6 original empirical-only rows per tier, verified
    in `results/followup-mature-organ-transports-c-lake.csv` (6 P1m, 6 P0m).
    The "four" is the TS2 remaining 41-group family: 4 original
    empirical-only groups (Bladder macrophage, Fat CD8-positive alpha-beta
    T cell, Fat fibroblast, Lymph_Node All_Cells), verified in
    `results/followup-mature-organ-transports-c-ts2.csv`. The fifth P1m
    original empirical-only row in that file is TS2_Blood_NK, which the
    report discusses separately (status DIED). Both counts are correct once
    the family is specified. The manuscript cites neither count.
11. Item 9b numbers. FILLED from the committed report
    `reports/followup-mature-between-subject.md` and verified against
    `results/followup-mature-between-subject-*.csv` by verifier A, with two
    corrections applied (see below). Results disorder subsection and
    Discussion/Limitations use the corrected values.

## Item 9b corrections applied (verifier A)

- CD4 total log variance ratio: the batch-arm q is 0.015992, not 0.007996.
  Only the depth-only arm is 0.007996. The 9b report text itself states it
  correctly; the brief's lumped value was wrong.
- Yazar per-cell-type distance: testable count is 6 of 28, not 6 of 19. The
  CSV has 28 non-pooled cell types for OLD vs MATURE: 6 TESTED, 22
  UNTESTABLE.
- All other 9b numbers confirmed: Yazar OLD vs MATURE raw ratio 1.542756,
  Levene p 3.461197e-05, CI 1.074845 to 2.623927; batch ratio 1.578469,
  p 0.136603, CI 1.062036 to 2.581521. Gao raw 1.048484, p 0.622569, CI
  0.387525 to 6.456563; batch 0.894063, p 0.823770, CI 0.482 to 2.42624.
  Yazar A1 -0.026112 p 3.549979e-06 n 925; A3 -0.077554 p 1.967207e-13 n 925.
  Gao A1 -0.007521 p 0.766955 n 104; A3 -0.061872 p 0.084606 n 105. Yazar
  gene-level variance: naive CD4 mean +0.277958 batch q 0.027986, +0.317544
  depth-only q 0.007996; naive CD8 mean +0.455961 batch q 0.007996,
  +0.621851 depth-only q 0.007996. Donors: Gao 109 (12 mature, 53 mid, 44
  old); Yazar 932 (59 mature, 193 mid, 680 old).

## Verification performed for this rewrite

- Verifier A (Item 9b): all numbers checked against CSVs; two corrections
  above, applied.
- Verifier B (P1m derivation, core transport, composition, drift, TS2):
  completed; only TS2 flagged, resolved above by direct CSV check.
- Verifier C (interneuron, Blocks 1, 2, 3, 6, 7, 8, concordance, Nov 17
  tests): completed; no corrections flagged in its report.
- Verifier D (organ transports A/B/C, full-depth, PsychAD, PBMC, Block 10):
  completed; all 8 items confirmed, no corrections.
- Own checks by the rewrite agent: Block 5 summary CSV (79 null-evaluated
  units; P1m both-pass 4 of 79; P0m both-pass 5 of 79; strict rescue 0 of
  48; all 4 P1m both-pass slopes negative; 4 of 5 P0m both-pass negative,
  only PsychAD IN_SST positive at +0.001226); concordance CSVs; interneuron
  replication CSV (VIP 4 vs 0.60 expected, pvalb 8 vs 1.00, sst 41 vs 11.66,
  block p 0.001; genome-wide r 0.7523, 0.6299, 0.6996); TS2 and Lake CSVs
  (see question 10); checker run; Tectonic compile.
- Reference [10] (Sikkema lung atlas) was uncited after the Results
  replacement; a [10] citation was added to the composition subsection.
  The checker then passed with cited set exactly equal to listed set.

## The story in short

Aging past maturity has no uniform transcriptional program. A mouse
signature derived from mature-to-old animals (1,914 genes) transports
almost nowhere in human cohorts: 2 of 48 core units pass both tests, both
microglia, both with negative slopes, and a broader unfiltered signature
rescues nothing (0 of 48). Cell-type shares shift in only 10 of 399
mature tests, all small. Whole-profile drift from a mature reference is
null in both species, as is gene coordination, and profile clouds widen
only in one blood cohort. Two narrow signals survive: microglial
expression falls with age in two brain cohorts, and interneuron age
effects replicate across two cohorts for one shared gene, TRMT11, with a
shared pattern of change in all three subtypes. Between-subject divergence
does not widen with age in heart or blood once pooling is accounted for,
though naive T cell gene-level variance widens in blood. Published
signatures mostly survive batch adjustment, so the artifact claim is
scoped, not broad. The nulls are bounded by power: most tests could only
detect effects of 0.7 to 1.6 standard deviations or more.

## Known limits carried into the manuscript

- P1m is sex-confounded (all four 30-month mice male) and its 14x size jump
  vs P1 is unexplained; both are disclosed.
- The mature mouse drift estimate rests on 14 mice; the Bayesian interval
  still includes effects up to +0.18 SD per decade.
- Propeller was not verified for the mature composition refit.
- The De Man human lung entropy replication stays blocked and stays a
  limitation.
- The 200-draw null miscalibration means no PBMC empirical-only pass is
  claimed; the 2,000-draw null is primary.
- The submission package at `~/workspace/agescape-final/plos-one-rev2/` was
  NOT refreshed from this commit in this task; it still reflects 0a022d8 and
  goes stale with this push. Refresh is a separate step.
