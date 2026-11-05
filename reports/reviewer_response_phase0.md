# Reviewer response Phase 0: audit and verification

Date: 2026-11-05

## Scope and frozen state

The instructed path `/home/reza/workspace/agescape-repo` did not exist. Reza approved use of `/home/reza/projects/agescape`. The draft path under `/home/reza/workspace/your_files/` also did not exist, and no `agescape-paper-draft.txt` was found under `/home/reza`. The numerical audit below covers every draft Results value enumerated in the Phase 0 instructions. It is not an independent parse of the unavailable full draft.

Starting and current remote head before this audit: `7076f6704198d809a3fd15cf71a8077354785d41`, `Test proliferation status as aging-effect moderator`. The working tree was clean. This is exactly the stated draft build commit, not a newer head.

No manuscript file was edited. No Phase 1 analysis was started. No submission, publication, external message, or other external action was performed. Read-only GitHub API requests were used for the required CI audit.

## 1. Draft-to-repo numerical verification

Every enumerated value matches the committed aggregate output to the stated rounding. No numerical discrepancy was found in this enumerated set.

| Draft claim | Repository value | Source |
|---|---|---|
| Mouse drift +0.1028, p = 0.00000398 | +0.1027967050, p = 3.982884103e-6 | `results/h3-cross-species-vector-global.csv:114` |
| Human vector +0.0365, p = 0.2218 | +0.0364830449, p = 0.2218424928 | `results/h3-cross-species-vector-global.csv:110` |
| Fly depth-adjusted vector -0.0467, p = 0.7958 | -0.0466848017, p = 0.7958307338 | `results/extended-afca-vector-global.csv:7` |
| Worm day-1 reference +1.6755 | +1.6754876901 | `results/extended-worm-vector-global-platform-audit.csv:2` |
| Worm day-6 reference +0.5319 | +0.5319289961 | `results/extended-worm-vector-global-platform-audit-reference6.csv:2` |
| NMR -0.4744, exact p = 0.2 | -0.4744304248, p = 0.2 | `results/h3-cross-species-vector-global.csv:122` |
| MSSM VIP FDR genes 294, both-cohort FDR 24 | 294 and 24 | `results/h3-human-neuronal-gene-replication-summary.csv:3` |
| MSSM pvalb FDR genes 333, both-cohort FDR 20 | 333 and 20 | `results/h3-human-neuronal-gene-replication-summary.csv:4` |
| MSSM sst FDR genes 123, both-cohort FDR 62 | 123 and 62 | `results/h3-human-neuronal-gene-replication-summary.csv:5` |
| Signature transport 3 of 252 | Three rows survive; the family table has 252 data rows | `results/human-transport-global-survivors.csv:2-4`; `results/human-transport-global-multiplicity.csv:2-253` |
| Off-program global +0.244, p = 0.278 | +0.2438136459, p = 0.2776838469 | `results/h3-frequentist-global.csv:5` |
| Entropy global -0.079, p = 0.699 | -0.0789334810, p = 0.6986741274 | `results/h3-frequentist-global.csv:2` |
| Proliferation cell-fraction contrast -0.0033 | -0.0033022381 | `results/proliferation-moderator-tests.csv:43` |
| Proliferation expression-vector contrast -0.0051 | -0.0050508004 | `results/proliferation-moderator-tests.csv:15` |
| Human clock predictive R-squared 0.6383 | 0.6383234725 | `results/h1-epigenetic-clock-human.csv:2` |
| Mouse clock predictive R-squared 0.3169 | 0.3169095012 | `results/h1-epigenetic-clock-mouse.csv:2` |

The donor-overlap interpretation does require correction, as detailed next. The former audit table's `exact_identifier_overlap = 0` was produced by a hard-coded zero rather than an identifier test. That value was not one of the enumerated draft Results numbers.

## 2. Donor-overlap resolution

Exact identity can be established for a substantial subset. Both sources use the same `Donor_N` identifier format. Matching identifier, age, sex, and brain bank gives:

| Public cohort | Public donors | PsychAD donors in the same bank | Exact identifier matches with concordant age and sex | Source |
|---|---:|---:|---:|---|
| MSSM | 214 | 112 | 103 | `results/m9-m10-donor-overlap-audit.csv:2` |
| HBCC | 134 | 172 | 131 | `results/m9-m10-donor-overlap-audit.csv:3` |

The corrected audit computes identifier intersections in `scripts/audit_m9_m10_donor_overlap.py:11-18` and applies them at `scripts/audit_m9_m10_donor_overlap.py:38`. The previous script assigned zero without performing this test. Only aggregate counts were exported; no controlled identifiers were added to Git.

PsychAD analyses using MSSM-bank donors overlap the public MSSM cohort. PsychAD analyses using HBCC-bank donors overlap the public HBCC cohort. M9-versus-M10 comparisons involving these sources must therefore be described as analyses of overlapping donors, not independent replication. Unmatched donors remain unmatched, and RADC and SEA-AD have no corresponding M10 bank arm in this audit (`results/m9-m10-donor-overlap-audit.csv:4-5`).

This confirms the draft's existing non-independence precaution rather than changing it into an independent-replication claim. The MSSM-versus-HBCC interneuron comparison is between different source banks and is not reclassified by this M9-versus-M10 audit.

## 3. Spearman 0.000 diagnosis

The correlation is a genuine exact zero. It is not a coding error and not a degenerate test.

- Sample size: seven pairs.
- Mouse slope variance: 0.0026644287.
- Human slope variance: 0.0001709051.
- Neither arm is constant and there are no tied ranks.
- Mouse ranks: 2, 5, 6, 3, 1, 7, 4.
- Human ranks: 3, 5, 6, 2, 7, 4, 1.
- Squared rank differences sum to 56, so `rho = 1 - 6 * 56 / (7 * (49 - 1)) = 0` exactly.
- SciPy `spearmanr` returns rho = 0.0 and p = 1.0.

The seven pairs are alveolar macrophage, classical monocyte, intermediate mouse monocyte versus non-classical human monocyte, B cell, NK cell versus natural killer cell, CD8-positive alpha-beta T cell, and vein endothelial cell. Their slopes and pair order are committed at `results/h2-composition-cross-species-homologs.csv:2-8`. The calculation is performed at `scripts/run_h2_composition.py:164-172`, and its output is `results/h2-composition-cross-species-test.csv:2`.

The diagnosis is an exact null rank correlation in a very small seven-pair comparison. It is not evidence of equivalence or absence of association.

## 4. Resampling gap inventory

### What exists

- The entropy and off-program pipeline restricts cells to at least 1,500 UMIs and downsamples every included cell without replacement to exactly 1,500 UMIs: `scripts/run_h3_transcriptional_disorder.py:197-227`.
- Entropy and off-program expression are calculated from that downsampled matrix: `scripts/run_h3_transcriptional_disorder.py:222-226` and `scripts/run_h3_transcriptional_disorder.py:248-254`.
- Existing sensitivity arms balance at 20 cells per animal, exclude pooled identifiers, apply stricter gene exclusions, and use detected off-program gene counts: `notes/transcriptional-disorder-methods.md:17-19`.
- Cell-to-cell variability has a separate comparison of multinomial resampling and without-replacement subsampling at 1,500 UMIs: `scripts/run_d2_noise_sensitivity.py:85-105`. This is not an entropy or off-program depth sweep.
- A separate direct expression-vector analysis exists without UMI downsampling. It computes donor-mean Hellinger and raw-count Euclidean distances and has depth-adjusted models: `scripts/run_h3_expression_vector.py:88-124`, `scripts/run_h3_expression_vector.py:151`, and `scripts/run_h3_expression_vector.py:169-184`. The cross-species version is documented at `notes/cross-species-expression-vector-methods.md:9-17`.

### What does not exist for entropy and off-program expression

- No 500-UMI, 3,000-UMI, or other depth sweep. Repository-wide downsampling calls in these analyses use only 1,500 UMIs.
- No entropy or off-program arm without downsampling.
- No library-size covariate in the donor regression. The model contains age and, when identifiable, sex only: `scripts/run_h3_transcriptional_disorder.py:142-165`.
- No direct cell-type average expression-vector test inside the entropy/off-program pipeline.

The final point must be scoped correctly. A direct whole-expression-vector analysis does exist elsewhere in the repository, including without downsampling. It is not an arm of the off-program or entropy null analysis.

## 5. CI status

Current master head and commit `7076f67` are the same commit, so their CI result is identical.

- Workflow: `CI`, run 17, run ID 41802134886.
- Head: `7076f6704198d809a3fd15cf71a8077354785d41`.
- Status: completed.
- Conclusion: failure.
- Job: `audit-and-test`, check-run ID 126480211264.
- Job start: 2026-11-04T08:46:54Z.
- Job completion: 2026-11-04T08:47:03Z.
- Duration: 9 seconds.

Sources: [workflow run API](https://api.github.com/repos/Abolfazl-Alipour/agescape/actions/runs?head_sha=7076f6704198d809a3fd15cf71a8077354785d41) and [check-runs API](https://api.github.com/repos/Abolfazl-Alipour/agescape/commits/7076f6704198d809a3fd15cf71a8077354785d41/check-runs). The workflow defines the job at `.github/workflows/ci.yml:9-36`. No cause is inferred from the failure notice alone.

## Audit commands and local verification

All commands ran from `/home/reza/projects/agescape`. The initial audit ran at commit `7076f6704198d809a3fd15cf71a8077354785d41`. The donor-overlap correction was then run from the working tree based on that commit.

```bash
git status --short --branch
git rev-parse HEAD
git log --oneline -5
git remote get-url origin
/home/reza/miniconda3/envs/agescape/bin/python scripts/audit_m9_m10_donor_overlap.py
/home/reza/miniconda3/envs/agescape/bin/python -m pytest -q tests/test_donor_overlap.py
/home/reza/miniconda3/envs/agescape/bin/python -m pytest -q
/home/reza/miniconda3/envs/agescape/bin/python scripts/check_app_data.py
git diff --check
```

The Spearman reproduction used:

```bash
/home/reza/miniconda3/envs/agescape/bin/python -c "import pandas as pd; from scipy.stats import spearmanr; d = pd.read_csv('results/h2-composition-cross-species-homologs.csv'); print(spearmanr(d.mouse_slope_per_month, d.human_slope_per_year))"
```

Local results after the audit correction: 54 tests passed, one Bayesian-only test skipped in the core environment, all 146 app-data anchors passed, and `git diff --check` passed. These local checks do not change or explain the red GitHub Actions run.

## Stop assessment

No enumerated draft Results number differed from the repository by more than rounding. Donor identity resolution confirms that M9/M10 source overlap is real and substantial; it does not justify a new independent-replication claim. Phase 0 is complete. Work stops here pending Reza's explicit approval of Phase 1.
