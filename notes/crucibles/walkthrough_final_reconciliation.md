# Walkthrough: Final Reconciliation, Numerical Auditing, and Manuscript Harmonization

This walkthrough documents the comprehensive reconciliation of the flagship manuscript [`report/manuscript-aging-mechanisms.md`](../../report/manuscript-aging-mechanisms.md), the project [`README.md`](../../README.md), and the underlying computational pipeline against all critical external peer reviews.

Every quantitative claim, sample size, statistical model specification, and unit of measurement has been audited against committed, reproducible CSV outputs. Zero synthetic or mock data remain in any pipeline, and all 16 automated tests pass cleanly.

---

## 1. Summary of Changes and Reconciliations

### 1.1. Exact Binding of Whole-Tissue Cortex Regression Statistics
* **Issue Raised:** Prose previously quoted unadjusted $+0.00049$ and adjusted $+0.00015$ from an earlier iteration.
* **Ground Truth Output:** In [`results/m10-cohort-summary.csv`](../../results/m10-cohort-summary.csv), the SEA-AD DLPFC whole-tissue pseudobulk regressions are:
  - Unadjusted OLS: slope $\beta = +0.000318/\text{yr}$, standard error $= 0.000125$, $p = 0.0118$
  - Adjusted OLS (controlling for Brain Bank, Sex, PMI, and Sequencing Depth): slope $\beta = +0.000027/\text{yr}$, standard error $= 0.000156$, $p = 0.8626$
* **Resolution:** Both the abstract, Section 2.1, Table 2.1, Section 3.2, and `README.md` have been updated to quote $\beta = +0.000027/\text{yr}$ ($p = 0.8626$) and $\beta = +0.000318/\text{yr}$ ($p = 0.0118$).

### 1.2. Disclosing Glial vs. Neuronal Divergence in Human Cortex
* **Issue Raised:** Earlier drafts emphasized glial signal without disclosing the prominent negative trajectories in excitatory neurons across the DLPFC atlas.
* **Ground Truth Output:** In [`results/m10-celltype-regressions.csv`](../../results/m10-celltype-regressions.csv), the mouse 257-gene aging signature demonstrates divergent cell-type-private dynamics:
  - **Glial Upregulation:**
    - Microglia: $\beta = +0.00142/\text{yr}$, unadjusted $p = 0.00199$, FDR $q = 0.0384$
    - Oligodendrocyte Precursor Cells (OPC): $\beta = +0.00057/\text{yr}$, unadjusted $p = 0.0486$, FDR $q = 0.444$
  - **Neuronal Downregulation:**
    - Excitatory intratelencephalic projection neurons (`EN_L3_5_IT_3`): $\beta = -0.00140/\text{yr}$, standard error $= 0.000267$, unadjusted $p = 3.79 \times 10^{-7}$, FDR $q = 1.02 \times 10^{-5}$
    - Excitatory deep-layer neurons (`EN_L6_CT`): $\beta = -0.00115/\text{yr}$, unadjusted $p = 0.00316$, FDR $q = 0.0427$
* **Resolution:** Section 2.1, Section 3.2, and Section 4 explicitly report this divergence. Rather than claiming a uniform cortical aging signature, the manuscript shows that whole-tissue averaging attenuates to null because positive glial trajectories and negative neuronal trajectories cancel each other out.

### 1.3. Mouse Donor Sample Sizes in Lung Composition Analysis
* **Issue Raised:** Section 2.5 and Mermaid diagrams had previously stated $N = 23$ mice for mouse lung composition remodeling.
* **Ground Truth Output:** In *Tabula Muris Senis* droplet (`data/raw/tms_droplet_raw.h5ad`), 23 mice are profiled across all organs combined, but exactly **16 mice** have lung tissue (`results/h2-composition-mouse-lung.csv`).
* **Resolution:** Corrected all occurrences of mouse lung donor counts to $N = 16$ mice across Section 2.5, the Abstract, and the Mermaid architecture diagram.

### 1.4. Naked Mole-Rat (NMR) Units, Statistical Rigor, and Equivalence Calibration
* **Issue Raised:** NMR skin slope was previously mislabeled as $+0.0035/\text{mo}$, and the TOST equivalence test was criticized for using an unscaled mouse chronological margin without acknowledging degrees of freedom with $n = 5$ animals.
* **Ground Truth Output:** In [`results/m6-statistical-summary.csv`](../../results/m6-statistical-summary.csv):
  - Model: OLS on normalized lifespan units (`age / max_lifespan`, where max lifespan is 4 yr for mouse, 37 yr for NMR)
  - NMR Skin (*H. glaber*, $n = 5$): slope $\beta = +0.0058$, standard error $= 0.0631$, $95\%\text{ CI } [-0.093, +0.100]$, unadjusted $p = 0.9332$
  - Equivalence margin: under the mouse chronological progression margin ($\delta_{\min} = 0.2023$), TOST rejection of $\beta \ge 0.2023$ yields $p_{\text{TOST}} = 0.0275$. Under a zero-margin threshold ($\delta = 0$), $p_{\text{TOST}} = 0.329$.
* **Resolution:** Section 2.2 and Table 2.2 now:
  - Accurately label units as per unit of normalized maximum lifespan.
  - Formally state that with $n = 5$ animals (3 degrees of freedom in OLS), the 95% CI $[-0.093, +0.100]$ successfully excludes mouse-rate chronological progression ($0.202/\text{unit}$), but is too wide to exclude lifespan-scaled drift ($\sim 0.022/\text{unit}$).
  - Demote causal and conclusive rhetoric ("proves" $\to$ "consistent with non-senescence kinetics").
  - Disclose the tissue-mismatch limitation (mouse lung/kidney discovery vs. NMR skin testing).

### 1.5. Full Disclosure of Mouse Model Specifications in Table 2.2
* **Issue Raised:** Table 2.2 previously omitted the unpooled OLS and late-adult plateau specifications from `m6-statistical-summary.csv`.
* **Resolution:** Expanded Table 2.2 to report all 6 specifications:
  1. *H. sapiens* DLPFC Whole Cortex (Adjusted OLS, $N = 284$): $\beta = +0.000027/\text{yr}, p = 0.8626$
  2. *H. sapiens* DLPFC Microglia (Adjusted OLS, $N = 284$): $\beta = +0.00142/\text{yr}, p = 0.00199$
  3. *M. musculus* Whole-Body (REML Cluster-Robust, $N = 23$): $\beta = +0.2053, p = 0.00064$
  4. *M. musculus* Whole-Body (Weighted Linear Mixed-Effects, $N = 23$): $\beta = +0.2023, p = 0.00295$
  5. *M. musculus* Whole-Body (Unpooled OLS, $N = 23$): $\beta = +0.0465, p = 0.1511$
  6. *M. musculus* Late-Adult Plateau (21m–30m OLS, $N = 10$): $\beta = -0.0056, p = 0.9560$
  7. *H. glaber* Skin (OLS Normalized Lifespan, $n = 5$): $\beta = +0.0058, p = 0.9332$

### 1.6. Universal Clock 2 Double-Exponential Inverse Link Function
* **Issue Raised:** Universal Clock 2 evaluation previously used raw linear scores, resulting in poor human MAE.
* **Ground Truth Output:** Implemented the official Lu et al. (2023) double-exponential inverse link function:
  $$\text{Predicted Relative Age} = \exp\left(-\exp\left(-\hat{\eta}\right)\right)$$
  $$\text{Predicted Chronological Age} = (\text{Predicted Relative Age} \cdot (\text{Max Lifespan} + 0.08)) - 0.08$$
* **Verification in `results/h1-universal-clock2-benchmark.csv`:**
  - Mouse Liver ($N = 339$): Pearson $r = 0.923$ ($p = 1.99 \times 10^{-141}$), $R^2 = 85.1\%$, $\text{MAE} = 0.36\text{ years}$ ($4.3\text{ months}$)
  - Human Multi-Tissue ($N = 661$, 8 organs): Pearson $r = 0.801$ ($p = 5.16 \times 10^{-149}$), $R^2 = 64.2\%$, $\text{MAE} = 5.06\text{ years}$

### 1.7. Absolute Link Sanitation and References Section
* **Issue Raised:** Local paths (`file:///home/reza/...`) and missing bibliographic citations.
* **Resolution:**
  - Replaced all local links in `report/manuscript-aging-mechanisms.md` and `README.md` with repo-relative paths (`results/...`, `scripts/...`, `figures/...`).
  - Added a comprehensive References section with 10 primary peer-reviewed citations (Horvath 2013, Lu 2023, Tabula Muris Senis 2020, Tabula Sapiens 2022, SEA-AD 2024, Enge 2017, Ruby 2018, Aitchison 1986, Benjamini-Hochberg 1995, Schuirmann 1987).
  - Updated correspondence email to `reza@agescape.org`.

---

## 2. Quantitative Verification Matrix

| Claim / Statistic | Source Script | Committed Result File | Verified Value | Status |
| :--- | :--- | :--- | :--- | :--- |
| **SEA-AD Whole Cortex (Unadjusted)** | `run_m10_pipeline.py` | `results/m10-cohort-summary.csv` | $\beta = +0.000318/\text{yr}, p = 0.0118$ | **Exact Match** |
| **SEA-AD Whole Cortex (Adjusted)** | `run_m10_pipeline.py` | `results/m10-cohort-summary.csv` | $\beta = +0.000027/\text{yr}, p = 0.8626$ | **Exact Match** |
| **SEA-AD Excitatory Neurons (`EN_L3_5_IT_3`)** | `run_m10_pipeline.py` | `results/m10-celltype-regressions.csv` | $\beta = -0.00140/\text{yr}, p = 3.79 \times 10^{-7}$ | **Exact Match** |
| **SEA-AD Microglia (Adjusted)** | `run_m10_pipeline.py` | `results/m10-celltype-regressions.csv` | $\beta = +0.00142/\text{yr}, p = 0.00199$ | **Exact Match** |
| **Mouse Whole-Body Slope (REML)** | `run_m6_pipeline.py` | `results/m6-statistical-summary.csv` | $\beta = +0.2053, p = 0.00064$ | **Exact Match** |
| **Mouse Late-Adult Plateau (21–30m)** | `run_m6_pipeline.py` | `results/m6-statistical-summary.csv` | $\beta = -0.0056, p = 0.9560$ | **Exact Match** |
| **NMR Skin Slope (Norm. Lifespan)** | `run_m6_pipeline.py` | `results/m6-statistical-summary.csv` | $\beta = +0.0058, p = 0.9332$ | **Exact Match** |
| **NMR TOST Equivalence ($0.2023$)** | `run_m6_pipeline.py` | `results/m6-statistical-summary.csv` | $p_{\text{TOST}} = 0.0275$ | **Exact Match** |
| **Mouse Lung Downsampled Noise (Fibroblast)** | `run_h3_entropy.py` | `results/h3-transcriptional-entropy-results.csv` | $\beta = +0.00109/\text{mo}, p = 0.0443$ | **Exact Match** |
| **Mouse Lung Fibroblast *Cdkn1a* (p21)** | `run_fibroblast_markers.py` | `results/fibroblast-markers-audit.csv` | $\beta = +0.1080/\text{mo}, q = 0.0322$ | **Exact Match** |
| **Mouse Lung Donors in Composition** | `run_h2_composition.py` | `results/h2-composition-mouse-lung.csv` | $N = 16$ mice | **Exact Match** |
| **Universal Clock 2 Mouse Liver MAE** | `run_h1_epigenetic_drift.py` | `results/h1-universal-clock2-benchmark.csv` | $\text{MAE} = 0.36\text{ yr}, r = 0.923$ | **Exact Match** |
| **Universal Clock 2 Human MAE** | `run_h1_epigenetic_drift.py` | `results/h1-universal-clock2-benchmark.csv` | $\text{MAE} = 5.06\text{ yr}, r = 0.801$ | **Exact Match** |

---

## 3. Automated Test Verification Results

All 16 unit, integration, and security tests pass:

```
============================== test session starts ==============================
platform linux -- Python 3.11.17, pytest-9.1.1, pluggy-1.7.0
rootdir: /home/reza/projects/agescape
collected 16 items

tests/test_clean_clone_and_numbers.py::test_anchor_lock_against_committed_results PASSED [  6%]
tests/test_clean_clone_and_numbers.py::test_app_data_check_script PASSED [ 12%]
tests/test_clean_clone_and_numbers.py::test_app_data_sizes_and_count PASSED [ 18%]
tests/test_hypotheses_2_and_3.py::test_h2_composition_outputs_exist PASSED [ 25%]
tests/test_hypotheses_2_and_3.py::test_h2_composition_simplex_invariant PASSED [ 31%]
tests/test_hypotheses_2_and_3.py::test_h2_composition_results_content PASSED [ 37%]
tests/test_hypotheses_2_and_3.py::test_h3_entropy_outputs_exist_and_valid PASSED [ 43%]
tests/test_hypotheses_2_and_3.py::test_fibroblast_markers_audit PASSED   [ 50%]
tests/test_hypotheses_2_and_3.py::test_h2_and_h3_figures_exist PASSED    [ 56%]
tests/test_hypotheses_2_and_3.py::test_h1_epigenetic_outputs_exist_and_accurate PASSED [ 62%]
tests/test_security_and_dua.py::test_no_controlled_data_or_tokens_in_git PASSED [ 68%]
tests/test_security_and_dua.py::test_zero_innerhtml_outerhtml_in_app PASSED [ 75%]
tests/test_security_and_dua.py::test_zero_inline_styles_in_index_html PASSED [ 81%]
tests/test_security_and_dua.py::test_no_third_party_fonts_in_index_html PASSED [ 87%]
tests/test_security_and_dua.py::test_csp_headers_snippet_exists PASSED   [ 93%]
tests/test_security_and_dua.py::test_no_local_user_paths_in_tracked_manifests PASSED [100%]

============================== 16 passed in 2.14s ==============================
```

Additionally, `python scripts/check_app_data.py` executes completely green, verifying all 146 contract anchors.
