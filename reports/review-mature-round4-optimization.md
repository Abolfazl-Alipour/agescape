# Round 4 optimization: accounting of all three hostile reviews

Manuscript reviewed: `report/plos-one/final-plos-one.tex` at commit `dd6f503`,
tag `mature-to-old-donor-scatterplots-2026-11-22` (abstract 299 words, 4 main
figures, 6 SI figures, 69 cited references).
Reviews: top-tier (commit `9859d2b`, REJECT, 4 BLOCKER, 10 MAJOR, 8 MINOR),
aging biology (commit `40499ed`, MAJOR REVISION, 5 BLOCKER, 14 MAJOR, 9 MINOR),
statistician (commit `a4200e1`, MAJOR REVISION, 2 BLOCKER, 14 MAJOR, 9 MINOR).
75 items total across the three reviews, with overlaps noted below.

Optimization branches: `round4-tex` (commit `a7cd7ed`) and `round4-figs`
(commit `476c670`), merged to master as `01c77f7`. Coordinator fixes on top
of the merge: the narrowing-dependence sentence (top-tier MAJOR 6, tex line
239), the eligibility-floor disclosure (top-tier MAJOR 3, tex line 263), and
the Data Availability re-tag (35 occurrences, lines 407, 411, 413).

Method: every review item was verified against its source CSV, report, or
cited work before changing anything. Numbers were changed only where a review
proved the current number wrong against its source file (concordance family,
PBMC count, refit q values). Every other number is unchanged.

## Verdicts after optimization

The top-tier REJECT stands on its own terms: it requires the positive-control
package, which is new analysis (see UNRESOLVED 1). The two MAJOR REVISION
verdicts are addressed item by item below; every fixable item in them is
fixed, and the remainder is new analysis or data, listed as UNRESOLVED.

## FIXED (top-tier review)

- BLOCKER 2 (624 family): demote demand REJECTED (see REJECTED below).
  Disclosure strengthened in Abstract, Results, and Limitations: the family
  was chosen after both correction outcomes were available, false discovery
  rate control is therefore not guaranteed, the 45-pass count is not a valid
  false-discovery-rate-controlled result under a prespecified family, and both
  promoted tests fail every sensitivity check.
- BLOCKER 4 (straw-man construct): FIXED in framing. The first-view framing
  is dropped for neutral "two patterns" language. The conclusion is
  downgraded in Abstract, Discussion, Conclusions, and the Q&A summary to
  "no shared gene-level signature was found under these four criteria", with
  pathway-level uniformity stated as untested. The four criteria themselves
  are unchanged. The untested alternatives (pathway test, clocks) remain
  UNRESOLVED (items 3 and 4 below).
- MAJOR 3 (eligibility floors): FIXED as disclosure (tex line 263). The tex
  now states plainly that no rationale for the three-mouse and 50-cell floors
  appears in the record and no sensitivity analysis varies them. The
  sensitivity analysis itself is UNRESOLVED (item 7 below).
- MAJOR 5 (Figure 2): FIXED. Points are now partial residuals after removing
  the model covariates, so every plotted cloud matches its fitted line in
  sign (six slope checks pass against the committed adjusted slopes). Panels
  relabeled (a) through (e). Gao Shannon is marked "p 0.0539, not
  significant" in the figure. File renamed `fig2-narrowing-panels.pdf`.
- MAJOR 6 (narrowing not three independent findings): FIXED (coordinator,
  tex line 239). The tex now states that scatter, Shannon entropy,
  Kullback-Leibler divergence, and the off-program counts are not
  independent lines of evidence for narrowing: each measures whether a cell's
  expression concentrates on fewer genes, so the narrowing claim rests on one
  pattern in several guises, seen in two cohorts.
- MAJOR 7 (revision-response residue): FIXED. Methods Table 3 blocks renamed
  to descriptive names; "earlier verdicts" became "verdicts from their
  original age ranges"; the signature is defined in the Results introduction;
  "the P1 analysis" became "this analysis"; "Follow-up decision rules"
  became "Follow-up analysis rules"; "revision" dropped from Data
  Availability.
- MAJOR 9 (stale links): FIXED. All 35 Data Availability occurrences now
  point to the optimization tag `mature-to-old-round4-optimization-2026-11-22`.
- MAJOR 10 (S1/S5 mislabels): FIXED. S1: "one test group is one
  tissue-by-cell-type group". S5: "joint donor-permutation null", error bars
  labeled as the null 95 percent confidence interval, leading with q.
- MINOR 1 (Fig 1 title overlap): FIXED. Donor-panel titles wrapped to two
  lines; render read back, no collision.
- MINOR 2 (Fig 4 casing, TRMT11 label): FIXED. Titles VIP, PVALB, SST,
  Broad GABAergic neurons; one common axis scale (plus/minus 0.012);
  TRMT11 label moved clear of the ring cluster.
- MINOR 3 (S2 mixed units): FIXED. Two panels: mouse per month (114 rows),
  human per year (285 rows); rank numbers match their own panel axes.
- MINOR 4 (Fig 3 punchline): FIXED. The figure now carries the rescue
  annotation: the unfiltered signature passes both tests in 0 of the 29
  units where the cleaned signature passes neither.
- MINOR 5 (S4 table): FIXED. Rebuilt as a proper grouped bar figure with
  Gao zeros shown (Harries 0/0, Peters 0/0, Sood 3 before and 0 after,
  de Magalhaes 0/0) and the Sood cohort split stated (3 Gao heart, 1 Yazar
  blood).
- MINOR 6 (Abstract caveats): FIXED. Heart Shannon null, interneuron
  caveats (moderated check unverified, third cohort failed), and GTEx whole
  blood positive are all in the Abstract.
- MINOR 7 (body jargon): FIXED. After the purge: P1m 0, P0m 0, both-pass 0,
  E1/E2 0, MDE80 0, multihot 0, EN_L3_5_IT_3 0, IN_SST 0, All_Cells 0 in
  prose. M9/M10 defined in Data Availability; LACA and HCATA expanded;
  "multihot" became a cohort-indicator score; cell-type codes expanded.
- MINOR 8 (small nits): FIXED. Fig 2 panel (b) y-label consistent with
  panels (a) and (c); S3 Gao footnote numbers traceable in the main text;
  S5 error bars defined in figure and caption.

## FIXED (aging-biology review)

- BLOCKER 1 (Gao flat raw means): FIXED. Results states the raw means
  (young 0.6029, old 0.6048) and that the Gao fall is a covariate-adjusted
  association, not a raw-mean fall; carried into Discussion and Conclusions.
- BLOCKER 2 (Yazar Shannon 3000-UMI rarefaction): FIXED. Results reports
  p 0.281 at 3,000-molecule rarefaction; carried into Discussion and
  Conclusions.
- BLOCKER 3 (lung headline vs study-aware test): FIXED. Results reports
  both the meta-combination (club minus 0.0061, p 0.58; multiciliated plus
  0.0098, p 0.39) and the study fixed-effects slopes (club minus 0.0426,
  p 0.0155; multiciliated minus 0.0261, p 0.0518) with per-study
  heterogeneity noted. The "five hold" headline is qualified in Abstract
  and Conclusions.
- BLOCKER 4 (Abstract hypothesis as fact, interneuron upgrade): FIXED.
  The Abstract frames the hypothesis as tested, not proven; TRMT11 is
  removed from the Abstract; the interneuron caveats are attached.
- BLOCKER 5 (S1 mislabel): FIXED (see top-tier MAJOR 10).
- MAJOR 1 (Gao Shannon specification range): FIXED. Results reports the
  range (primary minus 0.1001, p 0.0539; depth-by-cell-type minus 0.1013,
  p 0.0344; spline-residual p 0.0288; rarefied 1,500 p 0.0648; rarefied
  3,000 p 0.0727) and Gao Shannon is no longer used as corroboration.
- MAJOR 2 (effect sizes): FIXED. Headline slopes translated into absolute
  change across the studied spans (Gao scatter about 0.51 SD, Yazar scatter
  about 0.15, Yazar Shannon about 0.23, Gao Shannon about 0.71, stricter
  off-program count about 0.90 and 0.26). The three Sim 2021 donors
  (GSM4742854_Adult1 age 35, GSM4742855_Adult2 age 42, GSM4742856_Adult3
  age 42) are named as the low cluster near minus 3 to minus 4.5; the tex
  states explicitly that whether they drive the slope was NOT tested,
  because no leave-these-donors-out refit was run.
- MAJOR 3 (S2 units): FIXED (see top-tier MINOR 3).
- MAJOR 4 (S5 null label): FIXED (see top-tier MAJOR 10).
- MAJOR 5 (Oh numbers): NOT UPHELD. The reviewer's mismatch claim does not
  hold: it checked the abstract only. Verified by the coordinator against
  the Oh et al. Nature 2023 full text: mean pairwise Pearson r 0.29
  ("low-to-moderate correlation between the age gaps of different organs
  (mean pairwise Pearson r = 0.29)"); hazard ratio 2.37 per standard
  deviation of extra heart age, where 1 SD equals 4.1 years ("every 4.1
  years of additional heart age (one standard deviation) conferred an
  almost 2.5-fold increased risk of heart failure"). No number changed;
  the tex now states those exact wordings.
- MAJOR 6 (pathway-level uniformity): FIXED as downgrade (see top-tier
  BLOCKER 4). The test itself is UNRESOLVED (item 3 below).
- MAJOR 7 (cross-species asymmetry): FIXED. Both caveats (mouse mature 18
  months vs human mature 30 years, different life stages; shared names are
  not proven evolutionary matches) are carried into Results (transport,
  concordance) and Conclusions wherever a cross-species claim is made.
- MAJOR 8 (P1 typo): FIXED. "The P1 analysis" became "this analysis".
- MAJOR 9 (stale tag): FIXED (see top-tier MAJOR 9).
- MAJOR 10 (Long reconciliation): FIXED. The tex states the measure
  difference (Long: gene-level CV2 noise aggregated across genes; paper:
  cell-level Hellinger scatter between cells) and that the two results are
  compatible, not conflicting.
- MAJOR 11 (first-view strawman): FIXED (see top-tier BLOCKER 4).
- MAJOR 12 (Gao E2 singleton fraction): FIXED. Results reports the
  singleton-fraction rise (slope 0.151, p 0.0116) and explains why the
  count metric governs the verdict (prespecified primary endpoint; the
  singleton fraction is a secondary ratio metric).
- MAJOR 13 (interneuron demotion): FIXED. Heading, Discussion,
  Conclusions, and Fig 4 caption now say replication of correlation, not a
  replicated gene program; TRMT11 out of the Abstract; moderated check
  unverified and third cohort failed. All numbers kept (genome-wide
  correlations VIP r 0.7523, parvalbumin 0.6299, somatostatin 0.6996).
- MAJOR 14 (MSSM depth as candidate): FIXED. The depth decline is out of
  the Abstract candidate sentence; one biological candidate (naive CD4)
  remains; depth is kept only as a caution about the brain data.
- MINOR 1 (Fig 4 casing/axis): FIXED (see top-tier MINOR 2).
- MINOR 2 (Abstract moderated clause): FIXED. TRMT11 removed from the
  Abstract; the unverified result is not named there.
- MINOR 3 (Tabula Sapiens overstatement): FIXED. The "same duality" claim
  is softened to heterogeneity only.
- MINOR 4 (Abstract family clause): FIXED. The "five of fifteen" claim
  names the 431-test family.
- MINOR 5 (Fig 2 outliers): FIXED (see MAJOR 2 above).
- MINOR 6 (Abstract construct undefined): FIXED. The four criteria are
  named briefly in the Abstract.
- MINOR 7 (Jeffries age bins): FIXED. The bins (7 elderly donors 70-plus
  vs 10 adult donors 15 to 49) are in Introduction and Discussion.
- MINOR 8 (S6 E1/E2): FIXED. The figure spells out "off-program count"
  and "off-program stratum count"; no E1/E2 text in the figure.
- MINOR 9 (De Man entropy definitions): FIXED. De Man's entropy
  definition (Kullback-Leibler divergence) is stated beside the paper's
  per-cell Shannon entropy definition.

## FIXED (statistician review)

- BLOCKER 1 (concordance family): FIXED. Verified from
  `results/followup-mature-concordance-bh-check.csv`: the file computes
  exactly four families and no three-test family exists. Primary is now
  the mature two-test family (mouse within-tissue and human within-cohort),
  q 0.0440, PASSES 0.05. The six-row family (q 0.0460) is a sensitivity,
  not primary, because its three gene sets are near-duplicates of one
  hypothesis. Fixed in Abstract, Methods, Results, Discussion,
  Conclusions, S5 caption, and Q&A Q5. The aging-biology review accepted
  the three-test family as primary without checking whether it exists in
  the record; the statistician is the authority here.
- BLOCKER 2 (PBMC 0 of 35 vs 0 of 36): FIXED. Primary is 0 of 35 (22 Zhu
  plus 13 Thomson, pooled 35-test family). Zhu CD4-Treg is absent because
  only 2 of its 8 donors are mature, below the 3-donor floor. The alternate
  per-cohort implementation gives 0 of 36; the verdict is 0 either way.
  Fixed in Results, Appendix B, and Q&A Q2, with the unit set and BH
  family stated.
- MAJOR 1 (624 selection): REJECTED (see REJECTED below). Disclosure
  strengthened (see top-tier BLOCKER 2).
- MAJOR 2 (Figure 2 raw vs adjusted): FIXED (see top-tier MAJOR 5). The
  caption discloses that the points are partial residuals and that the raw
  clouds run the other way in the Gao scatter and Yazar pooled
  center-shift panels (raw Gao scatter plus 0.0637 vs adjusted minus
  0.0723; raw Yazar pooled center shift minus 0.0855 vs adjusted plus
  0.0400).
- MAJOR 3 (lung fixed-effects slopes): FIXED (see aging-biology BLOCKER 3).
- MAJOR 4 (scan family incoherence): FIXED. One family defines "clear the
  bar": the refit family, q below 0.05, 10 rows. The 431 family (8 rows)
  and the within-file family (14 rows) are reported as sensitivities, used
  consistently in Results, the Fig 1 caption, and the S2 caption.
- MAJOR 5 (propeller scope mix): FIXED. The tex now quotes the mature
  propeller p for lung vein endothelial (0.6413, q 0.9967) alongside the
  full-range p (0.0528), with the scope stated for both values; kidney
  collecting duct mature p 0.0588.
- MAJOR 6 (Gao Shannon): FIXED (see aging-biology MAJOR 1).
- MAJOR 7 (S2 rank conflict): FIXED (see top-tier MINOR 3). Highlight
  numbers match the rank axis within each panel.
- MAJOR 8 (P1m leave-one-donor-out): UNRESOLVED as new work (item 9
  below). Disclosure FIXED: Methods, Limitations, and Results state
  explicitly that no leave-one-donor-out check was run on the cleaned
  signature, so transport conclusions rest on a signature whose
  donor-robustness is untested.
- MAJOR 9 (Yazar All_Cells family of one): FIXED. The tex states the
  pooled Yazar All_Cells gene-fraction signal is uncorrected (family of
  one, q equals p 0.0030) and would not survive inclusion in the cell-type
  family; it is not presented as a corrected pass.
- MAJOR 10 (spleen per-depth family): FIXED. The drift family is defined
  explicitly; per-depth BH results are reported alongside the pooled
  28-test result (spleen at 3,000 molecules survives the per-depth family,
  q 0.0375; none survives the pooled family).
- MAJOR 11 (CD8 gate citation): FIXED. References 70 (Rudd et al., PNAS
  2011) and 71 (Egorov et al., Frontiers in Immunology 2018) added with
  titles and DOIs verified by web search and cited at the CD8 gate
  sentence; Q&A Q8 fixed the same way.
- MAJOR 12 (lung-vs-brain slope): FIXED. The sentence is replaced with
  the per-tissue slopes and p-values (brain minus 0.00265, p 0.799; lung
  minus 0.01367, p 0.893; kidney plus 0.11393, p 0.312).
- MAJOR 13 (stale tag): FIXED (see top-tier MAJOR 9).
- MAJOR 14 (S1 label): FIXED (see top-tier MAJOR 10).
- MINOR 1 (Fig 2 panel b caption): FIXED. "Each point is one donor by one
  cell type, except panel (c), which is one point per donor." (Panel
  letters updated for the rebuilt figure.)
- MINOR 2 (Fig 4 axis ranges): FIXED (see top-tier MINOR 2).
- MINOR 3 (S5 error bars): FIXED (see top-tier MAJOR 10).
- MINOR 4 (S4 Gao/Sood): FIXED (see top-tier MINOR 5).
- MINOR 5 (between-subject epsilon): FIXED. The tex states the epsilon
  used (1e-08 for the quoted values) and the Yazar A1 depth sensitivity
  (full depth minus 0.0261, p 3.55e-06; at 1,500 molecules minus 0.0819;
  at 500 minus 0.1192; at 3,000 minus 0.0600, p 0.0176).
- MINOR 6 (P1 typo): FIXED (see aging-biology MAJOR 8).
- MINOR 7 (S3 title): FIXED. Retitled to cover both rows: "Donor-to-donor
  distances and within-donor spread shrink with age in blood".
- MINOR 8 (debiased CD4 CSV): FIXED. `results/followup-mature-block1-debiased-slope.csv`
  committed, transcribing the debiased-slope computation from
  `reports/followup-mature-block1.md` section 5.7 with a provenance header
  (primary 0.00071, p 5.174e-08; debiased 0.000313, p 0.000638; 56.4
  percent calibration bias; plus the sensitivity rows). Transcription, not
  new analysis.
- MINOR 9 (survivor-bias Spearman): FIXED. Spearman rho plus 0.1429,
  p 0.1866, permutation p 0.1956.

## Q&A fixes

`~/workspace/goals/agescape-manuscript-submission/files/agescape-questions-and-answers.md`
(edited in place, not in the repo): Q2 PBMC, Q3 propeller/lung/refit q
values/624 disclosure, Q5 concordance family, Q6 Gao Shannon range plus
Yazar 3000-UMI rarefaction plus plain-word off-program forms, Q7 A1/A3
plain forms with epsilon and depth sensitivity, Q8 CD8 gate with Rudd/Egorov
plain-word citations, Q9 both-FDR plain forms, Q13, Q14 one-candidate
rewrite, the one-paragraph summary, and the source line. The Q&A keeps
P1m/P0m short forms under its upfront glossary definition.

## Reviewer numbers not upheld (verified by the coordinator or the tex agent)

1. Oh numbers (aging-biology MAJOR 5): the reviewer's mismatch claim does
   not hold; see aging-biology MAJOR 5 above. Verified by the coordinator
   against the Nature full text.
2. Refit q values 0.0086/0.0203 quoted in the tex were untraceable in any
   committed file and were replaced with the CSV `q_refit_family` values
   (headline five: 0.0109 x4, 0.0258; scan five: 0.0258, 0.0169, 0.0306,
   0.0499, 0.0499). Verified by the coordinator from
   `results/followup-mature-composition-human.csv` (club 0.0109,
   multiciliated 0.0258).
3. CD8 epsilon: the file value at epsilon 0.001 is 0.430 for naive CD8, not
   the review's 0.402; q values unchanged across epsilon.
4. Spleen pooled 28-test family: smallest q recomputed as 0.1063 by proper
   step-up BH, not the review's 0.150; conclusion unchanged (none survives
   the pooled family), so the q is not quoted.
5. No IN-SST variability-rise claim was found anywhere in the tex; the
   Jeffries sentences are the 124-gene sentences, which now carry the age
   bins.
6. The three-test concordance family (accepted by the aging-biology review)
   does not exist in the analysis record; see statistician BLOCKER 1.

## UNRESOLVED (needs new analysis or data; none attempted)

1. Positive-control package for the transport test: held-out-mouse
   validation of the mature signature, a progeroid or accelerated-aging
   model, and re-derivation of the signature with both sexes at old age.
   Needed: a held-out mouse cohort, a progeroid model dataset, and female
   30-month mice. This is the single change the top-tier review says would
   most improve the paper's chances.
2. GTEx whole-blood composition adjustment or deconvolution. Needed: a
   deconvolution or composition-adjusted re-analysis of the GTEx whole
   blood positive transport (plus 0.0119 per decade bin). The paper states
   it as an untested bulk-composition candidate.
3. Pathway-level uniformity test. Needed: a gene-set or pathway-level test
   of uniformity across cell types, organs, and species.
4. Clock benchmarks: Tyshkovskiy chronological and mortality clocks,
   Buckley per-cell-type clocks, Dobner pan-endothelial gene-set test.
   Needed: the clock models or gene sets run against the same cohorts.
5. De Man Baylor-lung replication. Needed: author-released barcode-level
   labels for the 32 Baylor lungs in GSE281219, or the owner's explicit
   permission for a re-annotation sensitivity. Stays a limitation.
6. Independent replication of the five fraction shifts or the narrowing
   slopes in a new cohort. Needed: an independent cohort with the same
   design (mature-to-old, donor-level).
7. Eligibility-floor sensitivity. Needed: re-running the signature
   derivation with the 3-mouse and 50-cell floors varied.
8. Investigation of the 14-fold cleaned-vs-young signature size gap beyond
   stating the cause is not established. Needed: a power-matched
   comparison, age-span decomposition, or method-matched control.
9. P1m leave-one-donor-out check. Needed: drop each of the four old mice
   in turn, re-derive the signature, and report the overlap with the full
   1,914-gene signature.
10. Mouse arm power. Design fact, not fixable by editing: drift on 14 mice
    aged 18 to 30 months, derivation on 4 old males, fraction refit on 7
    to 13 donors per tissue. The tex states the nine mouse fraction losses
    are power losses, not nulls.
11. Whether the three Sim 2021 donors drive the Gao scatter slope. Needed:
    a leave-these-donors-out refit. The tex states explicitly that this was
    not tested.

## REJECTED

1. Demote the 624-test family to a sensitivity analysis and report the
   431-family 15 passes as primary (statistician MAJOR 1, top-tier
   BLOCKER 2; round 3 R3-M3). REJECTED per the owner's locked decision
   (2026-11-20): "Keep 624 primary and justify it, disclosing the
   sensitivity failures." The disclosure is strengthened in Abstract,
   Results, and Limitations instead: the family was chosen after both
   correction outcomes were available, false discovery rate control is
   therefore not guaranteed, the 45-pass count is not a valid
   false-discovery-rate-controlled result under a prespecified family, and
   both promoted tests fail every sensitivity check. Do not re-litigate.

## Verification (coordinator)

- Checker `scripts/check_plos_manuscript.py`: PASS on the merged tree
  (abstract 299 words; 4 main figures; 71 cited references; 6 SI figures;
  section order correct; no forbidden dashes; no process residue).
- Tectonic compile outside the repo: PDF 39 pages, letter size, no errors
  (compiled in /tmp/compile-test by the tex agent; re-verified by the
  coordinator on the merged tree).
- Coordinator independently verified: the concordance two-test family and
  q 0.0440 in Abstract, Results, Discussion, Conclusions, and the S5
  caption; the PBMC 0 of 35 count and CD4-Treg explanation; the Fig 2
  rename and its partial-residual slopes against the committed adjusted
  slopes (six checks pass); the refit q replacement values from the
  composition CSVs; the Oh r 0.29 and HR 2.37 per SD (4.1 years) against
  the Nature full text; the jargon grep (P1m 0, P0m 0, both-pass 0, MDE80
  0 in prose); the Fig 2 render read back (all five clouds match their
  fitted lines in sign; Gao Shannon marked not significant).
- Nothing submitted anywhere.
