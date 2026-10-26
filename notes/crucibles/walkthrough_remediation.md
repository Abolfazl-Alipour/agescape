# Walkthrough: Scientific Remediation and Ground-Truth Verification

## Overview of Remediation

This walkthrough documents the comprehensive correction and empirical re-execution of all analyses across the repository (`https://github.com/Abolfazl-Alipour/agescape`), resolving all issues identified in external peer review:

1. **Universal Clock 2 Link Function Correction:** Replaced bare linear sums with Lu et al.'s (2023) double-exponential inverse link function $\exp(-\exp(-\hat{\eta}))$.
2. **True 1,500-UMI Single-Cell Rarefaction:** Eradicated double-normalization and implemented exact multinomial downsampling directly from raw counts in `data/interim/m1-qc-cleaned.h5ad`.
3. **Committed Marker Pipeline & Multiplicity Control:** Codified `scripts/run_fibroblast_markers.py` and output `results/fibroblast-markers-audit.csv` across all 11 candidate genes with Benjamini-Hochberg FDR correction.
4. **Two-Sided Compositional Remodeling & Adult Stratum:** Replaced the one-sided test with two-sided exact binomial tests and evaluated the adult primary stratum ($\ge 20$ yr) to eliminate pediatric confounding on human B cells.
5. **Naked Mole-Rat Evolutionary Negative Control:** Restored the NMR comparative arm demonstrating that the mouse aging signature remains completely flat in negligible-senescence NMR skin ($p = 0.925$; $p_{\text{TOST}} = 0.027$).
6. **Clean-Clone CI & Test Verification:** Updated `requirements-ci.txt`, added data manifests, hardened `scripts/check_app_data.py`, and verified that all 16 automated tests pass cleanly.

---

## Metric Comparison: Before vs. Corrected

| Domain / Analysis | Previous Erroneous / Biased Result | Corrected Ground-Truth Result | Verification Script & Output |
|---|---|---|---|
| **Universal Clock 2 (Mouse)** | $\text{MAE} = 1.02\text{ yr}$ ($44\%$ lifespan error) | **$\text{MAE} = 0.36\text{ yr}$ ($4.3\text{ mo}$)**, $r = 0.923, R^2 = 85.1\%$ | [`run_h1_epigenetic_drift.py`](../../scripts/run_h1_epigenetic_drift.py)<br>[`h1-epigenetic-clock-mouse.csv`](../../results/h1-epigenetic-clock-mouse.csv) |
| **Universal Clock 2 (Human)** | $\text{MAE} = 25.3\text{ yr}$ (omitted from text) | **$\text{MAE} = 5.06\text{ yr}$**, $r = 0.801, R^2 = 64.2\%$ | [`run_h1_epigenetic_drift.py`](../../scripts/run_h1_epigenetic_drift.py)<br>[`h1-epigenetic-clock-human.csv`](../../results/h1-epigenetic-clock-human.csv) |
| **Mouse Lung Fibroblast Noise** | $q = 0.0404$ (double-normalized, un-rarefied) | **Unadjusted $p = 0.0443$, FDR $q = 0.2713$** (1,500-UMI rarefaction) | [`run_h3_entropy.py`](../../scripts/run_h3_entropy.py)<br>[`h3-entropy-mouse-lung.csv`](../../results/h3-entropy-mouse-lung.csv) |
| **Classical Monocyte Noise** | Labeled "un-aged" based on $p = 0.8925$ | **Noise flat ($p = 0.8923$), but mean DE prominent (9 genes, median $|\log_2\text{FC}| = 1.99$)** | [`run_h3_entropy.py`](../../scripts/run_h3_entropy.py)<br>[`m2-aging-signal-ranked.csv`](../../results/m2-aging-signal-ranked.csv) |
| **Fibroblast *Cdkn1a* (p21)** | Uncommitted scratch session ($\beta = +0.1097$) | **Committed pipeline:** $\beta = +0.1080/\text{mo}, p = 0.0029, \mathbf{q = 0.0322}, R^2 = 0.777$ | [`run_fibroblast_markers.py`](../../scripts/run_fibroblast_markers.py)<br>[`fibroblast-markers-audit.csv`](../../results/fibroblast-markers-audit.csv) |
| **Fibroblast *Sirt6*** | Uncorrected headline claim ($\beta = -0.0220$) | **Committed pipeline:** $\beta = -0.0136/\text{mo}, p = 0.0245, \mathbf{q = 0.1350}$ (exploratory trend) | [`run_fibroblast_markers.py`](../../scripts/run_fibroblast_markers.py)<br>[`fibroblast-markers-audit.csv`](../../results/fibroblast-markers-audit.csv) |
| **Human Lung B-Cell Trajectory** | Collapse $q = 2.9 \times 10^{-5}$ (pediatric confound) | **Adult stratum ($\ge 20$ yr, $N=117$):** $\beta = -0.0027/\text{yr}, \mathbf{p = 0.7601}$ (flat) | [`run_h2_composition.py`](../../scripts/run_h2_composition.py)<br>[`h2-composition-human-lung.csv`](../../results/h2-composition-human-lung.csv) |
| **Cross-Species Composition** | Claimed "divergence $p = 0.967$" (one-sided error) | **Two-sided binomial $p = 1.000$, Spearman $\rho = -0.286$** (unproven concordance) | [`run_h2_composition.py`](../../scripts/run_h2_composition.py)<br>[`h2-composition-cross-species-homologs.csv`](../../results/h2-composition-cross-species-homologs.csv) |
| **Naked Mole-Rat Contrast** | Omitted from flagship manuscript | **Restored:** $\beta = +0.0035, p = 0.925$; TOST equivalence $p_{\text{TOST}} = 0.027$ | [`report/conserved-aging-programs.md`](../../report/conserved-aging-programs.md)<br>[`m6-statistical-summary.csv`](../../results/m6-statistical-summary.csv) |
| **Human OPC Cross-Species Transport** | Suppressed in manuscript | **Restored:** Replicates across 4/4 independent cohorts ($q = 0.021$) | [`m9-regressions.csv`](../../results/m9-regressions.csv) |
| **Human Microglia Transport** | Suppressed in manuscript | **Restored:** Replicates across 3/4 cohorts + M10 adult stratum ($p = 0.0019$) | [`m9-regressions.csv`](../../results/m9-regressions.csv)<br>[`m10-cohort-summary.csv`](../../results/m10-cohort-summary.csv) |

---

## Verification and Execution Log

### 1. Verification of H1 (Pan-Mammalian Epigenetic Clock)
```bash
/home/reza/miniconda3/envs/agescape/bin/python scripts/run_h1_epigenetic_drift.py
```
* **Output:**
  - Human Multi-Tissue ($N = 661$): Pearson $r = 0.801$ ($p = 5.16 \times 10^{-149}$), $R^2 = 0.642$, $\text{MAE} = 5.06\text{ years}$.
  - Mouse Liver Lifespan ($N = 339$): Pearson $r = 0.923$ ($p = 1.99 \times 10^{-141}$), $R^2 = 0.851$, $\text{MAE} = 0.36\text{ years}$ ($4.3\text{ months}$).
  - Publication figure generated at [`figures/h1-pan-mammalian-epigenetic-clock.png`](../../figures/h1-pan-mammalian-epigenetic-clock.png).

### 2. Verification of H3 (Single-Cell Rarefied Entropy)
```bash
/home/reza/miniconda3/envs/agescape/bin/python scripts/run_h3_entropy.py
```
* **Output:**
  - 24,458 qualified lung cells downsampled to exactly 1,500 UMIs.
  - Complete 9-cell-type table generated with DM and multivariate dispersion.
  - Lung fibroblasts: $\beta = +0.00109/\text{mo}, p = 0.0443, q = 0.2713$.
  - Classical monocytes: $\beta = +0.00008/\text{mo}, p = 0.8923, q = 0.8923$.
  - Publication figure generated at [`figures/h3-transcriptional-noise-vs-age.png`](../../figures/h3-transcriptional-noise-vs-age.png).

### 3. Verification of Fibroblast Marker Audit
```bash
/home/reza/miniconda3/envs/agescape/bin/python scripts/run_fibroblast_markers.py
```
* **Output:**
  - Full audit of 11 markers across $N = 10$ donors ($\ge 20$ cells).
  - *Cdkn1a* (p21): $\beta = +0.1080/\text{mo}, p = 0.0029, q = 0.0322, R^2 = 0.777$.
  - *Sirt6*: $\beta = -0.0136/\text{mo}, p = 0.0245, q = 0.1350, R^2 = 0.540$.
  - Publication figure generated at [`figures/fibroblast-candidate-markers.png`](../../figures/fibroblast-candidate-markers.png).

### 4. Verification of H2 (Compositional Remodeling)
```bash
/home/reza/miniconda3/envs/agescape/bin/python scripts/run_h2_composition.py
```
* **Output:**
  - Mouse lung: 6 cell types pass FDR $q < 0.05$ (T cells up, fibroblasts down).
  - Human adult lung: 4 cell types pass FDR $q < 0.05$ (club cells down, AT2 cells up).
  - Cross-species concordance: 3/7 concordant pairs, two-sided exact binomial $p = 1.000$, Spearman $\rho = -0.286$.
  - Publication figure generated at [`figures/h2-composition-shifts-mouse-vs-human.png`](../../figures/h2-composition-shifts-mouse-vs-human.png).

### 5. Automated Test Suite
```bash
/home/reza/miniconda3/envs/agescape/bin/pytest tests/ -v
```
* **Result:** `16 passed in 2.07s`.
