# Mature PBMC rerun at 30-plus

## Scope and design

This is the mature-to-old rerun of the PBMC analysis in
`reports/followup-pbmc-analysis.md`. Two public blood cohorts, analyzed
separately and never pooled: GSE213516 (Zhu) and GSE214546 (Thomson).
The design is mature donors aged 30 to 39 versus old donors aged at
least 60, fitted as a strict two-group contrast with sex as a
covariate. Donors under 30 are excluded entirely. This is a narrower
baseline than the original young group (age at most 40), so any signal
that was driven by donors under 30 should shrink or die here.

The signature is P1m, the mature mouse signature derived from exactly
18 versus 30 months: 1,914 genes, 1,724 human-mapped. Detected in the
PBMC matrices: 1,007 up genes and 709 down genes, 1,716 of 1,724, for
99.5 percent coverage. The source file is
`results/followup-mature-p1m-genes.csv`.

Endpoints, all run per cohort: P1m donor-score transport with a
200-draw expression-matched null; CLR composition (centered log-ratio
of cell-type shares, the H2 fraction test); Hellinger donor distance
from a cross-fitted 30 to 39 reference (the drift test); and targeted
pseudobulk expression of six genes in CD8 naive T cells. The P1m
matched null uses 200 draws, not 2,000, because the frozen PBMC
pipeline sets NULL_ITERATIONS to 200. Its p-value floor is 1/201, or
0.004975.

A cell type is eligible when it has at least 6 mature donors with at
least 3 per group. The donor is the unit in every test. No
cell-level test is used anywhere.

## Cohorts and donor counts

| Cohort | Mature 30 to 39 | Old 60-plus | Total at 30-plus | Eligible types | Hellinger eligible |
| --- | --- | --- | --- | --- | --- |
| GSE213516 (Zhu) | 3 (ages 30, 36, 36) | 7 (ages 60, 60, 65, 66, 66, 74, 77) | 10 | 22 of 28 | 23 of 28 |
| GSE214546 (Thomson) | 11 (ages 30 to 35) | 9 (ages 60 to 63) | 20 | 13 of 13 | 13 of 13 |

Source: `results/followup-mature-pbmc-donor-counts.csv` (2 rows) and
`results/followup-mature-pbmc-eligibility-all.csv` (82 rows). Zhu is
severely underpowered: only 3 baseline donors, so every Zhu null below
is a weak null. The Thomson batch model has 15 parameters on 20
donors with 5 residual degrees of freedom; it is identifiable but
thin. All 9 covariate models are identifiable in both cohorts
(`results/followup-mature-pbmc-covariate-identifiability.csv`).

Runtime: Zhu 29.6 seconds, Thomson 37.9 seconds, total 67.6 seconds,
peak memory 1,209 and 1,330 MB
(`results/followup-mature-pbmc-runtime.csv`). Both cohorts completed.
The runtime file lists model_used as "Gemini 4 Neon XHigh", but that
string is hardcoded in `scripts/run_mature_pbmc.py`; the run was
executed directly with no delegated model, the same as the Batch A
reruns. Do not cite a model for this run.

## Verdicts

Overall: MIXED. Two of the four original PBMC claims partially
survive at 30-plus; two die.

| Endpoint | Verdict at 30-plus | Original claim |
| --- | --- | --- |
| P1m transport, composite test | DOES NOT PROVE: 0 of 35 pass both tests | Thomson P1: 11 of 13 both-pass |
| CLR composition | PARTIAL: CD8 naive fall survives in Thomson only | CD8 naive fall in both cohorts |
| Hellinger drift | NEGATIVE displacement in Zhu only, 4 types; Thomson null | Zhu 3 types negative; Thomson null |
| Targeted genes in CD8 naive T | DOES NOT PROVE: 0 of 12 pass | Zhu 3 of 6 pass; Thomson 0 of 6 |

## P1m transport: 0 of 35 pass both tests

Primary model, fully adjusted. No cell type passes the composite test
in either cohort. The smallest composite q is 0.6299 (Zhu cDC) and
0.8657 (Thomson, 7 types tied). All 11 original Thomson both-pass
rows died. This independently confirms the Batch C finding that every
original Thomson PBMC both-pass row dies under P1m.

Source: `results/followup-mature-pbmc-p1m-models.csv`, primary rows
only. BH family is per cohort, 22 tests in Zhu and 13 in Thomson.

| Cohort | Cell type | n | Effect | p | q | Random p | Random q |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GSE213516 | CD14-MC | 10 | +0.0209 | 0.729 | 0.857 | 0.00498 | 0.00576 |
| GSE213516 | CD14-MC-PPBP | 10 | -0.0109 | 0.881 | 0.970 | 0.00498 | 0.00576 |
| GSE213516 | CD16-MC | 10 | +0.1038 | 0.0859 | 0.630 | 0.00498 | 0.00576 |
| GSE213516 | CD4-Naive | 10 | +0.0813 | 0.276 | 0.640 | 0.00498 | 0.00576 |
| GSE213516 | CD4-Tem | 10 | +0.0602 | 0.454 | 0.769 | 0.00498 | 0.00576 |
| GSE213516 | CD4-Tm | 10 | +0.0788 | 0.320 | 0.640 | 0.00498 | 0.00576 |
| GSE213516 | CD8-CTL | 10 | +0.0545 | 0.507 | 0.797 | 0.00498 | 0.00576 |
| GSE213516 | CD8-Naive | 10 | +0.0825 | 0.268 | 0.640 | 0.00498 | 0.00576 |
| GSE213516 | CD8-Tem | 10 | +0.0751 | 0.381 | 0.698 | 0.00498 | 0.00576 |
| GSE213516 | Intermed-MC | 10 | +0.0039 | 0.973 | 0.973 | 0.00498 | 0.00576 |
| GSE213516 | Megakaryocytes | 10 | -0.0022 | 0.960 | 0.973 | 0.811 | 0.811 |
| GSE213516 | Memory-B-CRIP1 | 10 | +0.0902 | 0.200 | 0.640 | 0.00498 | 0.00576 |
| GSE213516 | Memory-B-IL32 | 10 | +0.1043 | 0.312 | 0.640 | 0.00498 | 0.00576 |
| GSE213516 | NK-FCER1G | 10 | +0.0396 | 0.617 | 0.803 | 0.00498 | 0.00576 |
| GSE213516 | NK-GZMH | 10 | +0.0283 | 0.740 | 0.857 | 0.00995 | 0.0109 |
| GSE213516 | NK-GZMK | 10 | +0.0379 | 0.619 | 0.803 | 0.00498 | 0.00576 |
| GSE213516 | NK-IL7R | 10 | +0.0439 | 0.621 | 0.803 | 0.0149 | 0.0156 |
| GSE213516 | NK-S100A2 | 9 | +0.1797 | 0.119 | 0.640 | 0.00498 | 0.00576 |
| GSE213516 | NK-SELL | 9 | +0.1284 | 0.213 | 0.640 | 0.00498 | 0.00576 |
| GSE213516 | Naive-B | 10 | +0.0718 | 0.287 | 0.640 | 0.00498 | 0.00576 |
| GSE213516 | cDC | 10 | +0.1429 | 0.0351 | 0.630 | 0.00498 | 0.00576 |
| GSE213516 | pDC | 8 | +0.1424 | 0.0774 | 0.630 | 0.00498 | 0.00576 |
| GSE214546 | B cell progenitor | 20 | +0.0164 | 0.443 | 0.866 | 0.00498 | 0.00924 |
| GSE214546 | CD14+ Monocytes | 20 | -0.0050 | 0.774 | 0.895 | 0.0398 | 0.0575 |
| GSE214546 | CD16+ Monocytes | 20 | +0.0035 | 0.826 | 0.895 | 0.358 | 0.388 |
| GSE214546 | CD4 Memory | 20 | +0.0097 | 0.629 | 0.866 | 0.0299 | 0.0485 |
| GSE214546 | CD4 Naive | 20 | +0.0244 | 0.178 | 0.866 | 0.00498 | 0.00924 |
| GSE214546 | CD8 Naive | 20 | +0.0166 | 0.403 | 0.866 | 0.00498 | 0.00924 |
| GSE214546 | CD8 effector | 20 | +0.0103 | 0.666 | 0.866 | 0.0547 | 0.0647 |
| GSE214546 | Dendritic cell | 20 | -0.0344 | 0.221 | 0.866 | 0.00498 | 0.00924 |
| GSE214546 | Double negative T cell | 20 | +0.0108 | 0.576 | 0.866 | 0.00498 | 0.00924 |
| GSE214546 | NK cell | 20 | +0.0014 | 0.950 | 0.950 | 0.741 | 0.741 |
| GSE214546 | Platelets | 19 | -0.0334 | 0.125 | 0.866 | 0.00498 | 0.00924 |
| GSE214546 | pDC | 20 | +0.0111 | 0.551 | 0.866 | 0.0498 | 0.0647 |
| GSE214546 | pre-B cell | 20 | +0.0114 | 0.532 | 0.866 | 0.00498 | 0.00924 |

Effects are old-minus-mature differences in the P1m donor score
(mean log expression of up genes minus down genes). No depth, CMV,
batch, or pool variant passes either: 0 composite passes in all 74
non-primary P1m rows as well.

### The 200-draw matched null is uninformative here

Do not cite the empirical-only passes as evidence. 29 of 35 primary
rows pass the empirical test alone (21 of 22 in Zhu, 8 of 13 in
Thomson), but 26 of 35 sit exactly at the 200-draw floor
(random p 0.004975), and the observed effects sit a median of 5.06
null standard deviations from the null mean (range 0.16 to 11.62).
The null is far tighter than the P1m score's own sampling
distribution: for Zhu cDC, the regression standard error is 0.0548
while the random-effect standard deviation is 0.0151, a 3.6-fold gap.
Random 1,716-gene sets average out donor noise; the co-expressed P1m
genes do not. The original 139-gene P1 null did not show this
pattern at this scale: 4 floor hits in Zhu and 11 in Thomson, with
Thomson still yielding 11 genuine both-pass rows. For P1m in these
tiny cohorts, the empirical test passes nearly everything and
therefore tests nearly nothing. Only the composite verdict counts,
and it is 0 of 35.

The random null also has nonzero means in five Zhu monocyte and
dendritic types (CD14-MC -0.0516, CD14-MC-PPBP -0.0675, Intermed-MC
-0.0395, CD16-MC -0.0363, cDC -0.0249), a transcriptome-wide age
shift that the test correctly centers on. It does not rescue the
calibration problem above.

## CLR composition: CD8 naive fall survives in Thomson only

Primary model, exploratory all-label BH family (22 tests in Zhu, 13
in Thomson). Source:
`results/followup-mature-pbmc-clr-models.csv`, endpoint
CLR_composition, model primary.

| Cohort | Cell type | n | CLR effect | p | q |
| --- | --- | --- | --- | --- | --- |
| GSE214546 | CD8 Naive | 20 | -1.319 | 0.000444 | 0.00578 |
| GSE214546 | CD4 Memory | 20 | +0.335 | 0.00601 | 0.0391 |
| GSE214546 | CD16+ Monocytes | 20 | +0.596 | 0.0242 | 0.105 |
| GSE214546 | Double negative T cell | 20 | -0.528 | 0.0369 | 0.120 |
| GSE213516 | Megakaryocytes | 10 | +0.635 | 0.0262 | 0.455 |
| GSE213516 | CD8-Naive | 10 | -0.954 | 0.196 | 0.479 |
| GSE213516 | CD16-MC | 10 | +0.583 | 0.0943 | 0.455 |

Full 35-row table is in the CSV. No other type passes in either
cohort. Zhu best q is 0.455.

Targeted two-label family (CD8 naive and non-classical monocytes,
the same family the original report used):
`results/followup-mature-pbmc-clr-targeted.csv`.

| Cohort | Cell type | Expected | CLR effect | p | Targeted q | Exploratory q |
| --- | --- | --- | --- | --- | --- | --- |
| GSE213516 | CD8-Naive | decrease | -0.954 | 0.196 | 0.196 | 0.479 |
| GSE213516 | CD16-MC | increase | +0.583 | 0.0943 | 0.189 | 0.455 |
| GSE214546 | CD8 Naive | decrease | -1.319 | 0.000444 | 0.000889 | 0.00578 |
| GSE214546 | CD16+ Monocytes | increase | +0.596 | 0.0242 | 0.0242 | 0.105 |
| Both | MAIT | decrease | untestable | - | - | - |

MAIT is untestable in both cohorts: the released label is absent or
not harmonizable, as in the original run.

Robustness of the two Thomson passes, exploratory q:

| Cell type | Primary | log total cells | CMV | Batch | Pool | Raw fraction |
| --- | --- | --- | --- | --- | --- | --- |
| CD8 Naive | 0.00578 | 0.00957 | 0.00693 | 0.283 | 0.00851 | 0.00162 |
| CD4 Memory | 0.0391 | 0.0571 | 0.00693 | 0.170 | 0.0760 | 0.206 |

CD8 Naive survives every sensitivity except batch, the same pattern
as the original run (original batch q 0.139). CD4 Memory passes only
the primary model and fails log total cells, batch, pool, and raw
fraction. Its original strict effect was nearly identical (+0.312,
exploratory q 0.0507) to the mature effect (+0.335, q 0.0391), so
this is the same effect crossing an arbitrary line, not a new
finding. Do not present CD4 Memory as a mature discovery.

Comparison with the original: CD8 naive CLR fell in both cohorts
originally (Zhu -1.606, targeted q 0.0150; Thomson -1.424, targeted
q 0.0000647). At 30-plus the Zhu effect shrinks to -0.954 and loses
significance (p 0.196, MDE80 2.18, so Zhu cannot detect even the
original effect size), while Thomson keeps the full effect (-1.319,
MDE80 0.903). The monocyte rise keeps its direction in both cohorts
but passes no exploratory family at 30-plus.

## Hellinger drift: negative displacement in Zhu, null in Thomson

Primary model is depth-adjusted (the Hellinger primary variant adds
mean log UMI as a covariate). Source:
`results/followup-mature-pbmc-hellinger-models.csv`, model primary,
33 fit rows (20 Zhu, 13 Thomson). All four Zhu passers have negative
effects: donor profiles move closer to the mature 30 to 39 reference
with age, opposite the drift prediction.

| Cohort | Cell type | n | Effect | p | q |
| --- | --- | --- | --- | --- | --- |
| GSE213516 | CD4-Tm | 10 | -0.0126 | 0.00227 | 0.0256 |
| GSE213516 | NK-GZMK | 10 | -0.0137 | 0.00256 | 0.0256 |
| GSE213516 | cDC | 10 | -0.00653 | 0.00593 | 0.0315 |
| GSE213516 | CD4-Naive | 10 | -0.00993 | 0.00629 | 0.0315 |

Thomson: 0 of 13 pass; smallest q is 0.149 (Platelets, n 9).
Untestable: Zhu CD4-Treg, NK-S100A2, and pDC (fewer than 3 donors
per group in primary and unadjusted models); Thomson Platelets in the
batch-depth model only (rank or residual degrees of freedom).

Comparison with the original strict depth-adjusted family: 3 Zhu
types passed originally (CD8-Tem -0.0083, NK-GZMK -0.0086, cDC
-0.0055), all negative. At 30-plus, NK-GZMK and cDC survive, CD8-Tem
dies (effect -0.0102, same direction and larger, but p 0.082,
q 0.103: lost on precision, not on sign), and CD4-Tm and CD4-Naive
are new passers. CD4-Treg passed only the original exact-age design
and is untestable at 30-plus. The direction of the surviving signal
is unchanged: negative, against drift.

## Targeted genes in CD8 naive T cells: 0 of 12 pass

Six genes, primary model, BH family of 6 per cohort. Source:
`results/followup-mature-pbmc-pseudobulk-targeted.csv`, model primary.
Expected direction is increase for all six.

| Cohort | Gene | n | Effect | p | q | Original effect |
| --- | --- | --- | --- | --- | --- | --- |
| GSE213516 | CCL5 | 10 | +0.856 | 0.385 | 0.462 | +1.823 |
| GSE213516 | CLIC1 | 10 | +0.763 | 0.0980 | 0.301 | +0.954 |
| GSE213516 | S100A4 | 10 | +0.630 | 0.116 | 0.301 | +0.989 |
| GSE213516 | PRF1 | 10 | +0.601 | 0.194 | 0.301 | +0.698 |
| GSE213516 | PTGER2 | 10 | +0.804 | 0.201 | 0.301 | +1.029 |
| GSE213516 | GZMA | 10 | -0.0250 | 0.977 | 0.977 | +0.850 |
| GSE214546 | CCL5 | 20 | +1.122 | 0.0281 | 0.0844 | +1.059 |
| GSE214546 | CLIC1 | 20 | +0.316 | 0.0608 | 0.110 | +0.231 |
| GSE214546 | S100A4 | 20 | +0.317 | 0.0731 | 0.110 | +0.258 |
| GSE214546 | PRF1 | 20 | +0.0620 | 0.753 | 0.787 | -0.0265 |
| GSE214546 | PTGER2 | 20 | +0.661 | 0.0274 | 0.0844 | +0.589 |
| GSE214546 | GZMA | 20 | -0.0924 | 0.787 | 0.787 | -0.166 |

Original: Zhu 3 of 6 passed BH (CLIC1, S100A4, CCL5); Thomson 0
of 6. At 30-plus: 0 of 12 pass in either cohort. Best q is 0.0844
(Thomson CCL5 and PTGER2). The direction pattern is mostly
preserved: all six Zhu effects stay positive except GZMA, which
flips to -0.0250, and Thomson is 5 positive with only GZMA negative.
PRF1 flips sign in Thomson, from -0.0265 originally to +0.0620 at
30-plus, but both are null. Every Zhu effect attenuates
versus the original. The gene-level replication dies at 30-plus in
both cohorts.

## Agreement with Yazar blood at 30-plus

Yazar is the other blood cohort, and PBMC is an independent check
on it at 30-plus. The two agree on the negative claims and partly
on the positive ones.

- No positive mature transport in either blood dataset. PBMC: 0 of
35 P1m composite passes. Yazar: the only P1m both-pass units at
30-plus are negative-direction (All_Cells effect -0.0151, CD14-low
CD16-positive monocyte -0.0154). Neither dataset supports a
positive mature signature signal in blood.
- No positive drift in either blood dataset. PBMC Hellinger: Zhu
shows significant negative displacement, Thomson null. Yazar drift
at 30-plus is null in the mature refits, and the full-depth redo
finds no positive center shift within cell types.
- CD8 naive cells change with age in both datasets, but the
quantity differs. PBMC Thomson: the CD8 naive share falls
(CLR -1.319, q 0.00578). Yazar: naive CD8 within-cell scatter
rises (SD slope 0.0498 per decade, p 5.92e-05, q 4.74e-04, 344
donors). Both point at naive CD8 T cells as the one blood type
that moves with age past 30, in different measurements.
- Naive CD4 does not replicate across the two. Yazar naive CD4
dispersion is the sole Block 1 PROVES at 30-plus (slope 0.00071
per decade, p 5.17e-08, q 8.28e-07, 777 donors). PBMC naive CD4
is null on every endpoint: CLR fraction (Thomson +0.141, p 0.311;
Zhu -0.0536, p 0.926), Hellinger (Thomson -0.00120, p 0.509; Zhu
CD4-Naive is a negative-displacement passer, opposite drift),
and P1m transport (Thomson +0.0244, q 0.866; Zhu +0.0813,
q 0.640). The Yazar naive CD4 dispersion signal has no PBMC
counterpart in fraction, drift, or transport.

## Limits

- Zhu has 3 baseline donors. Every Zhu null is weak, and the Zhu
CD8 naive composition null cannot detect even the original effect
(MDE80 2.18 versus original effect -1.61). Only the Thomson
results carry weight here.
- The P1m matched null is miscalibrated for a 1,716-gene set in
cohorts this small; see the null section above. The 200-draw run
is also weaker than the 2,000-draw primary null used elsewhere in
the rebuild. A 2,000-draw rerun would not fix the calibration gap,
which comes from the gene-set correlation structure, not the draw
count.
- CD8 naive composition fails the batch sensitivity in Thomson
(batch q 0.283), as it did originally (batch q 0.139). The
composition claim is not batch-robust.
- CD4 Memory passes only the primary CLR model. It is a threshold
crossing of the original effect, not a new finding.
- Zhu genome-build adjudication is unresolved from public records
(GEO says hg19, the primary Methods say GRCh38, the released
files carry no assembly identity). Zhu gene results are
symbol-level only, as in the original run.
- The run used the frozen original PBMC pipeline with no code
changes except the mature donor filter and the P1m gene set. The
model_used column in the runtime file is hardcoded and names no
real model for this run.

## Files

- Report: `reports/followup-mature-pbmc.md` (this file)
- Script: `scripts/run_mature_pbmc.py` (545 lines; mature donor
filter plus P1m gene set on the frozen original pipeline)
- Combined results: `results/followup-mature-pbmc-p1m-models.csv`
(109 rows), `results/followup-mature-pbmc-clr-models.csv` (144
rows), `results/followup-mature-pbmc-clr-targeted.csv` (6 rows),
`results/followup-mature-pbmc-hellinger-models.csv` (111 rows),
`results/followup-mature-pbmc-pseudobulk-targeted.csv` (42 rows),
`results/followup-mature-pbmc-donor-counts.csv` (2 rows),
`results/followup-mature-pbmc-covariate-identifiability.csv`
(9 rows), `results/followup-mature-pbmc-eligibility-all.csv`
(82 rows), `results/followup-mature-pbmc-runtime.csv` (2 rows).
There is no combined donor-scores file; the per-cohort
p1m-donor-scores files hold 216 rows (Zhu) and 259 rows (Thomson).
- Per-cohort files: `results/followup-mature-pbmc-{GSE213516,
GSE214546}-{donor-counts, eligibility, covariate-identifiability,
p1m-donor-scores, p1m-models, pseudobulk-targeted, clr-models,
clr-targeted, composition-donor-counts, hellinger-donors,
hellinger-models, hellinger-eligibility}.csv` (24 files)
- Test: `tests/test_followup_mature_pbmc.py` (asserts every
headline number in this report from the CSVs; all pass)
- Original comparison sources: `reports/followup-pbmc-analysis.md`,
`results/followup-pbmc-{GSE213516, GSE214546}-{p1-models,
clr-models, hellinger-models}.csv`

## Verification

Every BH q in this report was recomputed independently from the
p-values with statsmodels multipletests: 438 q values checked across
the six combined result files, 0 mismatches. The P1m, CLR,
Hellinger, and targeted-gene headline numbers were read directly
from the CSV rows cited in each table. The original-comparison
numbers come from `reports/followup-pbmc-analysis.md` and the
original result files listed above.
