# Reviewer response Phase 2: decisive analyses on existing data

## Scope and provenance

Phase 2 started from Phase 1 commits `ab6750aba98621eae9b9986add1a66fa8eedef39` and `c3e1269`, with the latter adding the propeller implementation. No new single-cell dataset was collected. The manuscript was not edited. Exact commands and protocols are in `notes/reviewer-phase2-methods.md`.

## 1. Off-program expression and entropy at full depth

The depth sweep used fixed young-reference marker definitions and animal-level means. Global models standardized outcomes within cell type and clustered standard errors by animal.

| Depth | Model | Off-program effect per 15 months | 95% CI | p value | Entropy effect per 15 months | 95% CI | p value |
|---|---|---:|---:|---:|---:|---:|---:|
| 500 UMIs | Age + sex | +0.2537 | -0.2237 to +0.7311 | 0.2751 | -0.0716 | -0.4992 to +0.3559 | 0.7260 |
| 1,500 UMIs | Age + sex | +0.2438 | -0.2349 to +0.7225 | 0.2948 | -0.0789 | -0.5135 to +0.3557 | 0.7041 |
| 3,000 UMIs | Age + sex | +0.1726 | -0.3506 to +0.6958 | 0.4928 | +0.0080 | -0.4144 to +0.4304 | 0.9683 |
| Full depth | Age + sex | +0.1743 | -0.3557 to +0.7044 | 0.4940 | -0.0906 | -0.5517 to +0.3704 | 0.6812 |
| Full depth | Age + sex + library size | +0.1590 | -0.3358 to +0.6539 | 0.5038 | -0.0602 | -0.3522 to +0.2318 | 0.6668 |

Source rows: `results/reviewer-phase2-depth-global.csv:2-17`. The 1,500-UMI point estimates exactly reproduce the original global effects. The new t-based p-values differ slightly from the earlier normal-approximation output; the conclusions do not change.

No depth produced a significant positive global off-program association. The prespecified withdrawal criterion was therefore not met. The full-depth null supports the previous null interpretation, but the confidence interval remains wide and equivalence was not established in Phase 1. The result does not prove biological stability or exclude depth-sensitive local effects. Cell-type results are in `results/reviewer-phase2-depth-cell-types.csv`, and animal-level aggregate inputs are in `results/reviewer-phase2-depth-donor-results.csv`.

The separate direct cell-type average expression-vector analysis was rerun without downsampling. In mouse lung, the unadjusted Hellinger global effect was +0.2914 SD per 15 months, 95% CI -0.1203 to +0.7030, p = 0.1653. After mean-log-UMI adjustment it was +0.2626, 95% CI +0.0806 to +0.4445, p = 0.00468. Source: `results/h3-expression-vector-global.csv:2-3`. This endpoint measures displacement from young profiles, not off-program expression or entropy.

## 2. Adult-only and nonlinear human models

Every public human vector and composition donor was already at least 20 years old. The minimum age across both analyses was 20. Adult-only refitting therefore uses the same donor sets, and no full-lifespan-to-adult donor change exists for these public analyses. Age distributions are committed in `results/reviewer-phase2-human-age-distributions.csv:2-19`. PsychAD developmental donors belong to the separate signature-transport analysis, not these public vector and composition inputs.

The nonlinear comparison used a fixed three-degree-of-freedom cubic B-spline. Its F test asks whether the nonlinear fit improves on the linear fit; it is not a test of a positive age slope.

- Vector analysis: MSSM had two linear-FDR hits and no nonlinear-FDR hits. Kidney_SCR had no linear-FDR hits and two nonlinear-FDR hits. All other vector cohorts had zero hits under both corrections.
- Composition analysis: MSSM had six linear-FDR hits and one nonlinear-FDR hit. Integrated lung had two linear-FDR hits and no nonlinear-FDR hits. Kidney Living Donor had one linear-FDR hit and no nonlinear-FDR hit. Other composition cohorts had zero hits under both corrections.

Sources: `results/reviewer-phase2-human-vector-nonlinear.csv` and `results/reviewer-phase2-human-composition-nonlinear.csv`. Nonlinearity and cohort-specific age coverage can change individual conclusions. These results should not be summarized as confirmation of linear age effects.

## 3. Annotation stability

No per-cell annotation-confidence probability was available in the TMS or public human inputs. Doublet scores were not relabeled as annotation confidence. Excluded-cell and unassigned-cell counts by age are in `results/reviewer-phase2-annotation-exclusion-by-age.csv`.

Age-related exclusion differences were detected in integrated human lung, Spearman rho = +0.4334 with q = 0.0251 across 48 ages, and mouse pancreas, rho = +1.0 across four age groups. Sources: `results/reviewer-phase2-annotation-exclusion-age-tests.csv:10` and `:22`. The mouse pancreas result has only four age groups and should not be generalized.

The fixed canonical-marker analysis found eight FDR-significant changes: seven in human and one in mouse. Sources: `results/reviewer-phase2-marker-stability-summary.csv:2-3`. Examples include increased microglial C1QA in MSSM, slope +0.00522 log CP10K per year with q = 0.000280, and decreased mouse cardiac-fibroblast COL1A1, slope -0.0397 with q = 0.00617. Sources: `results/reviewer-phase2-marker-stability.csv:55` and `:8`. Full results include 187 tests in `results/reviewer-phase2-marker-stability.csv`.

These findings do not demonstrate identity loss. Marker expression can change without loss of cell identity, and a targeted marker panel cannot validate every annotation. However, age-related marker changes and differential cell recovery could contribute to composition associations. Identity loss therefore cannot be excluded as an explanation for any affected composition result.

## 4. Clock-overlap resolution

GSE199979 contains 339 predominantly female BXD mouse liver profiles on the Horvath mammalian methylation array. The GEO record and original mouse study were checked, together with the Universal Clock publication and its complete publisher supplementary workbook.

The workbook's 70 sheets contain no occurrence of `GSE199979`, the tested identifier `GSM6000839`, or the tested sample title `H1228_B6D2F1_CD_F`. Sample-level membership in clock development therefore remains unresolved. This is not evidence of no overlap. Source: `results/reviewer-phase2-clock-overlap-audit.csv:2`.

GSE184211 remains documented as included in the consortium resource used to develop the clocks (`results/h1-training-overlap-audit.csv:2`). Both clock evaluations remain benchmark reproductions, not independent validation.

## Deliverables

The report and all supplementary tables are committed with exact reproduction commands in `notes/reviewer-phase2-methods.md`. Local verification gave 65 passing tests and one Bayesian-only skip in the core environment; all 146 app-data anchors passed. Compilation and `git diff --check` passed. No manuscript edit, submission, publication, or external contact occurred.
