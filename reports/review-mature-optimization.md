# Agescape manuscript optimization: issue-by-issue accounting

Parent commit: c5dac33 (round 1 review commit, head on origin/master at the
time of this pass). Base manuscript: the mature-to-old rewrite at 683a74d.
This pass edits `report/plos-one/final-plos-one.tex` and rebuilds the four
figures from committed result tables. Nothing was submitted anywhere.

Verification: `scripts/check_plos_manuscript.py` PASS (abstract 300 words,
four figures, 69 cited references, section order correct, no forbidden
dashes, no listed process residue). Tectonic compile outside the repo:
`final-plos-one.pdf`, 34 pages (was 33). Figures built by
`scripts/build_optimization_figures.py` from committed CSVs, with asserted
values: Fig 1 mature composition forest (9 mouse DIES, 5 human SURVIVES,
1 UNTESTABLE, from `results/mature-only-composition-verdicts.csv`);
Fig 2 mature drift forest (mouse global minus 0.0137, human global minus
0.0071, six per-organ rows, from `results/followup-mature-drift-mouse.csv`
and `results/followup-mature-drift-human.csv`); Fig 3 transport verdict
bars (P1m 2/13/4/29, P0m 3/15/5/25, from
`results/followup-mature-p1m-block5-summary.csv`); Fig 4 interneuron
slopes with TRMT11 labeled (from
`results/followup-mature-interneuron-genes.csv`). All four
`\includegraphics` calls now reference the PDF (vector) files. The
orphaned young-reference drift figure and concordance figure were deleted
(`fig2-mouse-drift`, `fig3-concordance`).

Counts: 47 BLOCKER and MAJOR issues across the three rounds. FIXED 36,
UNRESOLVED NEEDING NEW ANALYSIS 10, REJECTED WITH REASON 1. MINOR issues:
7 fixed, 1 process note needing no action.

Dispositions: FIXED means the tex now states the corrected wording or
number, verified against a CSV read in this pass. UNRESOLVED NEEDING NEW
ANALYSIS means the review demands a new analysis or new data; no result was
invented, and the tex now states honestly what was not tested. REJECTED
WITH REASON means Reza's locked decision or the files force it.

## Round 1: statistician (MAJOR REVISION, 3 BLOCKER, 14 MAJOR, 2 MINOR)

### R1-B1: concordance "nominal only" while q 0.0440 and q 0.0460 pass
FIXED. Recomputed both BH corrections from
`results/followup-mature-concordance-mouse.csv` and `-human.csv`. The
primary family is now declared as the three-test family (mouse
within-tissue p 0.0220, mouse cross-tissue p 0.1029, human p 0.8152, all on
the primary detection-filtered gene set): mouse within-tissue q 0.066, so
"nominal only" is now correct. The two-test family (q 0.044) and six-row
family (q 0.046) are reported as sensitivities with reasons: the two-test
family excludes the cross-tissue test without justification, and the
six-row family counts near-duplicate sensitivity gene sets as separate
tests. New Results wording: "Mouse within-tissue concordance is nominal
only: observed mean Pearson 0.2105 against a null mean of 0.1465
(p = 0.0220, q = 0.066 under the primary three-test family...)".

### R1-B2: Methods says families were defined before correction; line 143 admits the 624 choice was deferred
FIXED. The Methods sentence now reads: "Each analysis family was defined
before correction, except the primary composition family, whose choice was
deferred until after both correction outcomes were available, as stated
above." The 624 family stays primary per Reza's locked decision
(2026-11-20); the timing is disclosed in the abstract ("chosen after both
correction outcomes were available"), Results, Discussion, and Limitations.
"locked" was replaced by "declared" everywhere.

### R1-B3: ectopic headline quotes the E2_stratum bounded-depth rows
FIXED. Results now quote the primary E1 full-depth endpoint: "The primary
off-program set (E1) is null in both cohorts (Gao minus 0.063, p = 0.468;
Yazar minus 0.019, p = 0.073). The bounded-depth stratum variant (E2)
falls in both cohorts (Gao minus 0.127, p = 0.031; Yazar minus 0.033,
p = 0.0031), opposite the prediction." Verified in
`results/followup-mature-fulldepth-*-ectopic.csv`.

### R1-M1: mature composition q values cannot be reproduced from the mature CSVs
FIXED. The q values printed for the five survivors (0.0118, 0.0279) are
`mature_q_global_bh_orig_family` in
`results/mature-only-composition-verdicts.csv`. The tex now names the
family: "The q values printed here use the original 431-test family, the
most conservative of the three families computed; the refit-family q values
are 0.0086 for four rows and 0.0203 for multiciliated cells, and the
within-file q values are 0.00091 to 0.0102."

### R1-M2: marrow presented as a concordance carrier but fails BH across seven tissues
FIXED. Results now read: "The signal is carried by limb muscle alone
(observed 0.428 against a null mean of 0.190, p = 0.001, q = 0.007 across
seven tissues); marrow (0.292 against 0.189, p = 0.019) fails correction
across tissues (q = 0.066) and is not a carrier." Verified by recomputing
BH across the seven per-tissue p values in
`results/followup-mature-concordance-mouse.csv`.

### R1-M3: concordance excess 0.0640 smaller than MDE80 0.0677; joint null mean 70 percent of observed
FIXED. Results now state: "The excess over the null (0.064) is smaller than
the minimum detectable effect at 80 percent power (0.068), and the joint
null mean is 70 percent of the observed correlation, so most of the
correlation is donor sharing across cell types, not age concordance."
Verified in `results/followup-mature-concordance-mde80.csv`.

### R1-M4: naive CD4 headline slope uncalibrated; finite-N calibration explains 56.4 percent
FIXED. Results, abstract, and Conclusions now carry the debiased slope:
"Finite-cell calibration accounts for 56 percent of that slope: the
debiased slope is plus 0.000313 (p = 0.000638), and 0.000295 after depth
adjustment (p = 0.00106)." The unadjusted arm (p = 0.814) and rarefied arm
(p = 0.797) are stated. Verified in `reports/followup-mature-block1.md`
(the debiased slope is report-only; no committed CSV contains it, and the
tex says so by citing the report in Data Availability, not by inventing a
CSV).

### R1-M5: two BH-significant Block 1 gene-fraction signals never reported
FIXED. Reported in the "No general rise" section: "in Gao, myeloid cells
show 0.6658 of genes higher in old donors (q = 0.0270), plausibly
ambient-RNA artifact; in Yazar, pooled All_Cells show 0.5487 (q = 0.0030,
single-test family). Neither changes the dispersion verdict." Verified in
`results/followup-mature-block1-gene-fractions.csv`.

### R1-M6: composite q values uncalibrated
FIXED. Methods now states: "Composite q values are uncalibrated screening
statistics from a parametric model; empirical q values from the
matched-random null are the calibrated evidence numbers." Empirical q now
leads in every quoted case: Yazar All_Cells (empirical q = 0.00425;
composite q = 4.03e-10), GTEx whole blood (empirical q = 0.0015; composite
q = 1.29e-22), Harries (empirical q = 0.004248; composite q = 7.511e-47),
Peters (empirical q = 0.0006535; composite q = 2.599e-33).

### R1-M7: 200-draw null cannot pass by construction
FIXED. Methods, PsychAD, and PBMC sections now state that at 200 draws no
unit can pass the empirical test at q below 0.05. PsychAD: "the smallest q
in that file is 0.0597, above 0.05". PBMC: "at 200 draws the smallest
achievable q exceeds 0.05, so no row can pass". Verified in
`results/followup-mature-psychad-m10-null200.csv` (0 pass at 200 draws,
4 at 2,000; minimum q at 200 draws 0.0597) and
`reports/followup-mature-pbmc.md` (26 of 35 rows at the 1/201 floor).

### R1-M8: naive CD8 scatter undisclosed batch dominance and jackknife bias
FIXED. The narrow-signals section now states: "Batch explains 82 percent
of the outcome variance (batch partial R2 0.817) and the jackknife bias is
28 percent of the slope." Verified in
`results/followup-mature-fulldepth-yazar-cell-drift.csv`.

### R1-M9: spleen omitted from the per-organ drift list
FIXED. Drift Results now list spleen: "spleen plus 0.453 (p = 0.057;
p = 0.0076 at 1,500 molecules and p = 0.0054 at 3,000 molecules)".
Figure 2 plots the spleen row. Verified in
`results/followup-mature-drift-mouse.csv` (per_tissue variant). Note: the
mature drift file has no Heart per-tissue row, so Figure 2 shows the six
available per-organ rows, not seven.

### R1-M10: fourteen full-depth primary endpoints never enumerated
FIXED. The full-depth section now enumerates them: "The fourteen primary
endpoints are, per cohort: Part 1 center shift, Part 1 scatter, Part 2
Shannon entropy, Part 2 KL methods, Part 2 KL square-root, Part 3 E1
ectopic count, and Part 3 E2 stratum ectopic count. The verdict rule is a
positive-direction significant slope with adequate power; the stratified
global rows carry no BH q because the rule is per endpoint, not a family."
Per-type BH-significant falls added: "Yazar natural killer scatter
(q = 1.33e-07), CD14-positive monocyte scatter (q = 0.00360), effector
memory CD8 center shift (q = 0.00215), and natural killer center shift
(q = 0.00215)." All four verified on the primary model in
`results/followup-mature-fulldepth-yazar-cell-drift.csv`.

### R1-M11: Discussion 83 percent / 96 percent numbers are full-range, use changing test universes, no table citation
FIXED. Discussion now reads: "In the full age range, not the 30-plus
restriction, about 83 percent of gene hits did not persist after pool
adjustment in the 981-donor blood cohort (10,123 of 160,239 tests to 1,741
of 134,036) and concordance fell by about half (0.2430 to 0.1024), while
over 96 percent did not persist after cohort adjustment in the 99-donor
cohort (10,021 to 391). The test universe shrinks between arms, so part of
the attenuation is donor and gene loss rather than confounding removal."
Source: `reports/followup-hcata-pmc12511636.md` lines 61 and 89, linked in
Data Availability.

### R1-M12: RPE survival rests on 13 donors and a single-unit family
FIXED. Results now read: "Retinal pigment epithelium is the one earlier
both-pass transport that survives under both tiers, in the negative
direction (effect minus 0.0211, q = 0.0443, empirical q = 0.0005), on 13
donors in a single-unit family."

### R1-M13: Gao Block 7 hits run opposite the prediction; 14 of 96 predicted changes did not attenuate
FIXED. Results now read: "with 96 predicted pass-to-fail gene changes, 14
of which did not attenuate after removal (central memory CD4 4 of 37,
naive CD4 9 of 56, naive CD8 1 of 3), while Gao gives two hits with zero
predicted genes each, both negative in direction (lymphoid minus 0.0292,
q = 0.0318; pericyte minus 0.0242, q = 0.0378), opposite the
survivor-bias prediction." Verified in
`results/followup-mature-block7-*.csv` via the Block 7 report values.

### R1-M14: secondary-five q values printed as 0.050 are actually 0.04989
FIXED. Both printed as 0.0499. The first three q values are now labeled
431-family q (0.028, 0.018, 0.033). Verified in
`results/followup-mature-composition-secondary-checks.csv`.

### R1-m1: interaction bootstrap p values unmentioned
FIXED. Appendix A now reads: "The mouse minus human interaction was minus
0.0663 unadjusted (p = 0.0515; bootstrap p = 0.072) and minus 0.0447 after
depth adjustment (p = 0.0051; bootstrap p = 0.03)." Verified in
`results/reviewer-rev2-mouse-human-interaction.csv`.

### R1-m2: Table 1 says 479 vector donors; Appendix A says 480 donors
FIXED. Appendix A now reads: "human plus 0.1891 against a null mean of
plus 0.2076 (p = 0.7003) across 131 units and 480 donors (the vector
analysis in Table 1 uses 479 donors; one additional donor meets
concordance eligibility)".

## Round 2: aging biology / theory / literature (MAJOR REVISION, 2 BLOCKER, 15 MAJOR, 2 MINOR)

### R2-B1: naive CD8 scatter dies under log cell-count adjustment, flips without pool effects, dies under 20-cell subsampling
FIXED. The scatter signal is demoted to inconclusive and removed from the
abstract's survived list, the Results section heading, the Discussion
pattern paragraph, and the Conclusions. Results now quote all arms:
primary 0.0498 (p = 5.92e-05, q = 4.74e-04); plus log total cells 0.0080
(p = 0.526, q = 0.727); no pool minus 0.0312 (p = 0.131); 20-cell
subsampling minus 0.0016 (p = 0.887); rarefied 0.0885 (p = 5.5e-06);
spline 0.0925 (p = 6.9e-05). The Methods binding rule was extended to
per-type scatter: "a primary pass that fails the log-total-cells arm is
inconclusive." The naive CD8 gate contamination note and the standing
donor-cloud spreading and share-fall results are kept. Prediction 3 now
requires prespecified cell-count adjustment for any scatter claim.

### R2-B2: naive CD4 dispersion primary is real but 56.4 percent is finite-N calibration bias; Block 1 report calls it not proof of a biological dispersion rise
FIXED. See R1-M4. Results, abstract, and Conclusions carry the debiased
slope (0.000313, p = 0.000638) and the mechanistic caveat: "The signal
survives the frozen rule, but its mechanism is unresolved, so it is
presented as a candidate, not as proof of a biological dispersion rise."

### R2-M1: title and abstract state a claim no current theory makes
FIXED. Title is now "No uniform transcriptional response from maturity to
old age in mouse and human under a stringent four-criterion construct".
Abstract opens: "Aging theories posit shared upstream mechanisms, yet
organs and cell types age differently."

### R2-M2: criterion 1 (same direction in EVERY cell type) makes the negative near-tautological; shared cores fail only on "every"
FIXED. The abstract closing now reads: "No uniform transcriptional
response from maturity to old age was detected under the study's stringent
four-criterion construct, which does not exclude shared pathway-level
cores." The Introduction already states the construct is deliberately more
stringent than shared hallmarks, modular conserved signatures, or
transportable clocks. Criterion 1 was not softened to a prevalence
threshold because the four-criterion construct is the paper's stated
scope; the compatibility statement is the chosen remedy.

### R2-M3: no pathway-level uniformity test reported
UNRESOLVED NEEDING NEW ANALYSIS. Exact analysis required: take the
per-type differential-expression lists, run pathway enrichment per cell
type, and test whether enriched pathways overlap across cell types more
than a permutation null allows, separately in mouse and human. Text fix
done: Introduction, Discussion, and Limitations state "No pathway-level
uniformity test was run, so the gene-level null does not rule out
pathway-level uniformity."

### R2-M4: Dobner 2026 cited only for specificity; its headline is a conserved pan-endothelial program
UNRESOLVED NEEDING NEW ANALYSIS. Exact analysis required: score the Dobner
pan-endothelial gene set in the paper's endothelial units with the same
directional composite and matched-random rules used for P1m transport.
Text fix done: the Introduction and Discussion now cite Dobner for its
headline (rising inflammation, declining angiogenesis, extracellular-matrix
remodeling, p21-linked senescence program) and state the Dobner gene set
was not tested for transport here.

### R2-M5: Jeffries "124 genes across all cell types" is wrong
FIXED. Introduction and Discussion now read: "124 genes fall across
multiple cell types, with 10 falling in all 13 types". The mature-to-old
scope limit is stated: "That shared signal spans infancy to old age, so
the mature-to-old window used here cannot test it." Verified by web search
in the review pass.

### R2-M6: Tyshkovskiy 2026 never tested; CDKN1A absent from P1m, LGALS3 flipped to DOWN
UNRESOLVED NEEDING NEW ANALYSIS. Exact analysis required: benchmark the
Tyshkovskiy chronological and mortality clocks (and the universal gene
set/module direction test) on the same mouse and human cohorts. Text fix
done: Introduction and Discussion state "Two genes that work names as
shared features behave differently here: CDKN1A is absent from P1m, and
LGALS3 reverses direction, up in the young-to-old signature and down in
the mature signature" and "We did not benchmark their chronological or
mortality clocks". Verified in
`results/followup-mature-p1m-genes.csv` and
`results/followup-mature-p1m-p1-overlap.csv` (Cdkn1a lost from P1, was UP;
Lgals3 in P1m as DOWN, mean log2FC minus 2.71, 1 unit, common sign flip).

### R2-M7: Sloan 2026 (same DLPFC tissue) never engaged
FIXED. Sloan is engaged in Results (MSSM depth paragraph: "Sloan and
colleagues, in the same tissue, report somatostatin-positive neuron
density falling with age [46]; that composition change offers an
alternative explanation for the somatostatin depth result that was not
tested") and in Discussion (shared oxidative phosphorylation and
cytoplasmic translation decline is a pathway-level signal of the kind this
paper did not test).

### R2-M8: Yang 2026 three-phase lifespan structure unreconciled with Block 8
FIXED. Discussion now reads: "Yang and colleagues report three major
phases of molecular change across the lifespan in the same PsychAD atlas
[9]; those phases include development, where age explains the most
variance, so they are compatible with finding no mature-to-old step change
here. The mature-to-old window discards the developmental variance
component where Yang finds the largest age effect."

### R2-M9: Tabula Sapiens 2.0 cited as a preprint although published in Cell Oct 2026
FIXED. New reference 69: "The Tabula Sapiens Consortium (2026). Tabula
Sapiens 2.0: A comprehensive transcriptomic atlas of human cell types.
Cell 189, 7212-7227. https://doi.org/10.1016/j.cell.2026.08.010". Methods
now cites [69] with no preprint version. Results add: "Nine donors cannot
test a shared senescence program, so this is a power statement, not a
substantive null. The published Tabula Sapiens 2.0 senescence analysis
finds the same duality as this paper: senescent cells share common
gene-expression modules alongside organ- and cell-specific patterns [69]."

### R2-M10: Buckley clocks generalize across species and were never benchmarked
UNRESOLVED NEEDING NEW ANALYSIS. Exact analysis required: benchmark the
Buckley per-cell-type aging clocks on the same mouse and human cohorts for
cross-species transport. Text fix done: abstract and Introduction state
transport failure covers only this fixed signature.

### R2-M11: MSSM depth decline is a sex-only model on four nested units; abstract drops the qualifier
FIXED as a text fix; the covariate refit is UNRESOLVED NEEDING NEW DATA.
Abstract now calls it "captured-library-depth decline in nested MSSM
neuron units, a technical observation with unresolved explanations".
Results and Conclusions carry the same qualifier plus nesting and the
unavailable-covariate list. The refit cannot be run on the current input:
the MSSM input contains no postmortem interval, brain-bank, RNA-quality,
processing-batch, or technical-principal-component fields.

### R2-M12: microglia-only transport is marginal, sign-flipping, never decomposed
UNRESOLVED NEEDING NEW ANALYSIS. Exact analysis required: score the
disease-associated microglia (DAM) gene set in the microglia units with
the same directional composite and matched-random rules, and decompose the
P1m composite by gene in the two passing units. Text fix done: Results
state the PsychAD microglia sign flip (plus 0.002381 full-range to minus
0.001786 mature), the undecomposed composite ("which genes drive the
negative microglia slopes is unknown"), and the marker-panel disagreement
(rising microglial C1QA in MSSM, slope plus 0.00522, q = 0.000280, while
the P1m composite falls). MSSM microglia q 0.0497 is already in the tex.

### R2-M13: GTEx whole blood, the one positive transport, is bulk and most parsimoniously composition
UNRESOLVED NEEDING NEW ANALYSIS. Exact analysis required: composition
adjustment or deconvolution of GTEx whole blood (adjust the P1m score for
estimated cell-type fractions, or score sorted-cell blood data), to test
whether the positive transport survives. Text fix done: Results state
"Bulk whole blood conflates composition with within-cell change, and this
positive transport was not tested for composition confounding, so it stands
as a bulk-composition candidate, not as within-cell transport." The same
statement is in Limitations.

### R2-M14: TRMT11 is statistically real and biologically near-empty
FIXED. Results now read: "TRMT11 is a replicated statistical association
of unknown significance: a tRNA methyltransferase with age slopes of
0.00205 to 0.00406 log counts per 10,000 per year, about 2 to 4 percent per
decade, expressed at about 4 counts per 10,000, with no established aging
biology. It is not presented as a mechanistic lead." Prediction 1
reframed the same way. Slopes and expression verified in
`results/followup-mature-interneuron-genes.csv`. Figure 4 labels TRMT11 in
each subtype panel.

### R2-M15: "locked" in the abstract contradicts the Methods on 624-family selection timing
FIXED. "locked" replaced by "declared" in the abstract, Discussion, and
Conclusions. The abstract carries the timing clause: "The 624-test
composition family, chosen after both correction outcomes were available,
yielded 45 passes".

### R2-MINOR-1: Long "six immune types" needs a source check
FIXED. The count was removed: "Long and colleagues reported elevated
individual-level noise with age in the same 981-donor OneK1K resource
[67]."

### R2-MINOR-2: "five of fifteen" in the abstract vs "17 core passes" in Results
FIXED. Results now reconcile them: "The fifteen are the passes under the
original 431-test family; the 624-test family counts seventeen core
passes, the two extra being the promoted tests below." Abstract now reads
"five of the fifteen full-range fraction changes".

## Round 3: top-tier journal reviewer (REJECT at top-tier standard, 4 BLOCKER, 9 MAJOR, 4 MINOR)

### R3-B1: central negative rests on a confounded, unstable probe with no positive control
UNRESOLVED NEEDING NEW ANALYSIS. Exact analyses required: (1) re-derive
P1m with both sexes represented at old age; (2) held-out-mouse validation
across derivation methods; (3) a progeroid or accelerated-aging positive
control. Text fix done: the probe limitations are stated in Methods
(three-organ signature, sex confound), Results, Discussion ("Transport
failure therefore remains a property of an unstable and sex-confounded
probe as well as possible species divergence"), Limitations, and
Prediction 4.

### R3-B2: no discovery of top-tier scale (three thin, disconnected narrow signals)
UNRESOLVED NEEDING NEW ANALYSIS. Exact work required: replicate at least
one narrow signal in an independent cohort with a prespecified analysis,
or produce a new discovery of top-tier scale; no text fix can supply it.
Text fix done: the signals are now two candidates with caveats and one
inconclusive finding, not three survived signals.

### R3-B3: naive CD4 headline number is the un-debiased slope
FIXED. See R1-M4 and R2-B2. Debiased slope 0.000313 (p = 0.000638) is in
Results, abstract, and Conclusions.

### R3-B4: abstract omits the direction of the only transport passes
FIXED. Abstract now reads: "only microglia passed both transport tests,
both in the negative direction". Results, Discussion, and Conclusions
carry the same direction statement.

### R3-M1: Figure 1 does not show the paper's result; omits the two promoted tests by its own caption
FIXED. Figure 1 is now the mature composition forest plot (15 rows: nine
mouse DIES, five human SURVIVES, bladder untestable noted in the caption),
built from `results/mature-only-composition-verdicts.csv`. The caption
agrees with the text and states the 624-family passes and promoted tests
are reported in the text.

### R3-M2: the two central results (drift null, transport failure) have no main-text figure; the only drift figure sits in Appendix A and shows the disowned young-reference analysis
FIXED. New Figure 2 (mature drift: mouse and human globals plus six
per-organ rows) and new Figure 3 (core 48-unit transport verdicts) are in
the main text. The young-reference drift figure was removed from
Appendix A.

### R3-M3: the 624-test family is post-hoc scope-stretching; make the mature 399-row table primary with its own figure and demote 624 to a sensitivity, or drop it
REJECTED WITH REASON. Reza's locked decision (2026-11-20, confirmed in
the optimization brief) keeps the 624 family primary; re-litigating it
would disobey the user's explicit order. Maximum fix applied instead: the
mature 399-row table is now Figure 1, the primary visual; the 624 family
stays primary with its selection timing disclosed in the abstract,
Methods, Results, Discussion, and Limitations, and the two promoted tests
disclosed as sensitivity-fragile and unreplicated.

### R3-M4: the 83 percent / 96 percent Discussion numbers have no traceable source
FIXED. See R1-M11. Traced to `reports/followup-hcata-pmc12511636.md`
lines 61 and 89, denominators given, full-range scope stated, report
linked in Data Availability.

### R3-M5: the paper does not test the literature claim it positions itself against (never benchmarks Tyshkovskiy clocks; four-criterion construct stricter than its target)
UNRESOLVED NEEDING NEW ANALYSIS. Same as R2-M6: benchmark the
Tyshkovskiy chronological and mortality clocks on the same cohorts. Text
fix done: title and abstract are scoped to the four-criterion construct,
and the Introduction states the clocks were not benchmarked.

### R3-M6: interneuron replication rests on one gene, unmoderated, with a failed replication cohort (SEA-AD); moderated check unverified; Figure 4 does not label TRMT11
FIXED. Results state all four facts in the main text: TRMT11 only,
unmoderated test, SEA-AD did not replicate (15 donors without dementia),
moderated mature intersection unverified. Discussion and Limitations carry
the same. Figure 4 now labels TRMT11 in each subtype panel. TRMT11 is
reframed as a replicated association of unknown significance (see R2-M14).

### R3-M7: a third of the transport family is untestable because source files are missing (82 of 252)
FIXED. The abstract now states "In the 170 testable of 252 units, only
microglia passed both transport tests". Results and Limitations state the
82 untestable units (71 absent source files, 11 too few 30-plus donors).
Restoring the four missing cohorts (Kidney Normal Adult, RADC, SEA-AD,
Kidney SN) would require the absent source files and is the only way to
run those 71 units; that restoration is not part of this pass.

### R3-M8: P1m is a three-organ signature, not a mouse signature
FIXED. Stated in Methods ("P1m is therefore a three-organ signature, not
a whole-mouse signature, and every transport claim in this paper is a claim
about that three-organ probe"), Results, Discussion, and Limitations.
Only 36 of 151 units were eligible for derivation (11 kidney, 12 lung,
13 marrow).

### R3-M9: MSSM depth decline should not be a headline signal
FIXED. Demoted to a technical observation in the abstract ("a technical
observation with unresolved explanations"), Results ("It is one technical
observation, not a biological signal, and it is reported as such"),
Discussion, and Conclusions. Nesting (broad GABAergic unit contains the
other three) and sex-only adjustment are stated.

### R3-MINOR-1: mature composition q values do not name their family
FIXED. See R1-M1.

### R3-MINOR-2: figures are raster, not vector
FIXED. All four `\includegraphics` calls now reference the PDF (vector)
files; the PNG pairs are retained alongside.

### R3-MINOR-3: condense the concordance subsection and move Figure 3 to the supplement
FIXED. The concordance subsection is now one paragraph in the main text
and the concordance figure was removed from the tex (the figure file was
deleted; the underlying data remain in
`results/followup-mature-concordance-*.csv`).

### R3-MINOR-4: process note, not a manuscript defect (683a74d patch restored the naive CD4 signal)
No action. The naive CD4 signal is present and verified against
`results/followup-mature-block1-residual-dispersion.csv` (primary slope
0.0007139, p 5.17e-08, q 8.28e-07).

## Notes and limits of this pass

- The debiased naive CD4 slope (0.000313, p = 0.000638) exists only in
`reports/followup-mature-block1.md`, not in any committed CSV. The tex
does not claim a CSV source for it.
- The mature drift file has no Heart per-tissue row, so Figure 2 shows six
per-organ rows. No heart value was invented.
- The abstract is exactly 300 words under the checker's count; any future
abstract edit must cut words first.
- The optimization report and review reports are linked in Data
Availability. The submission package in
`~/workspace/agescape-final/plos-one-rev2/` is still stale (matches the
superseded optimized commit 0a022d8) and must be refreshed from the final
post-optimization commit before any submission step.
