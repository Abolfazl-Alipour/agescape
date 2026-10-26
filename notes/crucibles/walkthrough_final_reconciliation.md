# Walkthrough: Final Reconciliation, Numerical Auditing, and Manuscript Harmonization

This walkthrough documents the programmatic reconciliation of the flagship manuscript [`report/manuscript-aging-mechanisms.md`](../../report/manuscript-aging-mechanisms.md), the project [`README.md`](../../README.md), and the underlying computational pipeline against all critical external peer reviews.

Every quantitative claim, sample size, statistical model specification, and unit of measurement has been audited directly against committed, reproducible CSV outputs.

---

## 1. Summary of Reconciliations and Audits

### 1.1. Exact Binding of Whole-Tissue Cortex Regression Statistics
* **Ground Truth Output:** In [`results/m10-cohort-summary.csv`](../../results/m10-cohort-summary.csv), the PsychAD DLPFC (284 control donors, ages 0.2–97 years, SD 24.0y) whole-tissue pseudobulk regressions are:
  - Unadjusted OLS: slope $\beta = +0.000318/\text{yr}$, standard error $= 0.000125$, $p = 0.0118$
  - Adjusted OLS (controlling for Brain Bank, Sex, PMI, and Technical PCs 1–3): slope $\beta = +0.000027/\text{yr}$, standard error $= 0.000155$, $p = 0.8626$
* **Resolution:** Abstract, Section 2.1, and `README.md` cite $\beta = +0.000027/\text{yr}$ ($p = 0.8626$) and $\beta = +0.000318/\text{yr}$ ($p = 0.0118$).

### 1.2. Disclosing Stratified Glial and Excitatory Neuronal Trajectories in Human Cortex
* **Ground Truth Output:** In [`results/m10-confounder-audit.csv`](../../results/m10-confounder-audit.csv) and [`results/m10-celltype-stratified-summary.csv`](../../results/m10-celltype-stratified-summary.csv):
  - **Glial Upregulation:**
    - Microglia: Pooled adjusted $\beta = +0.00142/\text{yr}, p = 0.0020$, within-family FDR $q = 0.0184$; Adult stratum ($\ge 20$ yr, $N = 231$) $\beta = +0.00200/\text{yr}, p = 0.0019$, FDR $q = 0.0518$.
    - Oligodendrocyte Precursor Cells (OPC): Pooled adjusted $\beta = +0.00065/\text{yr}, p = 0.0486$, within-family FDR $q = 0.2019$; Adult stratum $\beta = +0.00129/\text{yr}, p = 0.0049$, FDR $q = 0.0657$.
  - **Excitatory Neuronal Attenuation in Adults:**
    - Intratelencephalic projection neurons (`EN_L3_5_IT_3`): Pooled adjusted $\beta = -0.00140/\text{yr}, p = 3.79 \times 10^{-7}$, within-family FDR $q = 1.02 \times 10^{-5}$; Adult stratum ($\ge 20$ yr, $N = 225$) $\beta = -0.00060/\text{yr}, p = 0.0945$, within-family FDR $q = 0.4251$.
    - Upper-layer neurons (`EN_L2_3_IT`): Pooled adjusted $\beta = -0.00058/\text{yr}, p = 0.0180$; Adult stratum $p = 0.1487$.
* **Resolution:** Section 2.1 and the Abstract explicitly disclose the adult stratum for neurons: in adults, glia rise robustly and selectively, while excitatory neurons display a non-significant negative trend, demonstrating that pooled lifespan neuronal significance is substantially developmental.

### 1.3. Driver-Gene Decomposition of Glial Transport (D1)
* **Ground Truth Output:** In [`results/d1-glial-driver-decomposition.csv`](../../results/d1-glial-driver-decomposition.csv) and [`results/d1-glial-model-summary.csv`](../../results/d1-glial-model-summary.csv):
  - Microglia: Top drivers are *MAFB* ($p = 1.26 \times 10^{-7}$), *CD74* ($p = 5.56 \times 10^{-23}$), *CDKN1A* ($p = 4.71 \times 10^{-5}$), and *SPP1* ($p = 1.76 \times 10^{-4}$). The age slope survives removing a pre-declared 22-gene generic inflammatory/stress module ($\beta = +0.001081/\text{yr}, p = 0.0039$) and strengthens under the audited `P1_Clean` tier ($\beta = +0.002069/\text{yr}, p = 0.0015$).
  - OPC: Top drivers are *LPL* ($p = 2.69 \times 10^{-20}$), *AZGP1* ($p = 2.59 \times 10^{-26}$), and *CEBPD* ($p = 2.87 \times 10^{-12}$). The slope attenuates without generic inflammatory genes ($\beta = +0.000341/\text{yr}, p = 0.1483$), but strengthens under `P1_Clean` ($\beta = +0.001028/\text{yr}, p = 0.0343$).
  - Whole Cortex: Remains completely null under `P1_Clean` (Adjusted $\beta = -0.000209/\text{yr}, p = 0.4004$).

### 1.4. Naked Mole-Rat (NMR) Units, Statistical Rigor, and Equivalence Calibration
* **Ground Truth Output:** In [`results/m6-statistical-summary.csv`](../../results/m6-statistical-summary.csv):
  - Model: OLS on normalized lifespan units ($\text{slope}_{\text{norm}} = \text{slope}_{\text{year}} \times \text{max\_lifespan}$, where max lifespan is 4 yr for mouse, 37 yr for NMR).
  - NMR Skin (*H. glaber*, $n = 5$): slope $\beta = +0.0058$, standard error $= 0.0641$, $95\%\text{ CI } [-0.1981, +0.2098]$, $p = 0.9332$.
  - TOST Equivalence: under the same-normalized-lifespan margin ($\delta_{\min} = 0.2023$), the 90% CI is $[-0.1450, +0.1567]$, which sits inside $[-0.2023, +0.2023]$ ($p_{\text{TOST}} = 0.0274$). The 95% CI includes $\pm 0.2023$; standard TOST at $\alpha = 0.05$ uses the 90% CI.
  - Direct Mouse-vs-NMR Contrast: slope difference $= +0.1285$, standard error $= 0.0816$, $95\%\text{ CI } [-0.0314, +0.2885]$, $p = 0.1153$ (non-significant).
  - NMR Aged-vs-Young Exact Permutation: $p = 1.0000$.
  - Deletion of Zero-Margin TOST: The degenerate empty-interval zero-margin test has been completely removed across the manuscript, report, and all documentation. The defensible statistics (slope +0.0058, SE 0.0641, p = 0.9332, 90% CI [-0.1450, +0.1567], 95% CI [-0.1981, +0.2098], TOST delta = 0.2023, p = 0.0274) are preserved.
  - Mouse Late-Adult Interaction Test: In mice, late-adult subset slope is flat ($-0.0056, p = 0.9560$); a direct age-by-period interaction test yields $\beta_{\text{int}} = -0.1555, p = 0.563$, establishing that a true slope plateau is not statistically confirmed.

### 1.5. Noise Sensitivity Under Without-Replacement Subsampling (D2)
* **Ground Truth Output:** In [`results/d2-noise-sensitivity.csv`](../../results/d2-noise-sensitivity.csv):
  - Lung Fibroblasts: While multinomial downsampling showed a nominal unadjusted increase in Distance-to-Median noise ($\beta = +0.00109/\text{mo}, p = 0.0443$, FDR $q = 0.2713$), exact without-replacement hypergeometric subsampling to 1,500 UMIs/cell attenuates to null ($\beta = -0.000288/\text{mo}, p = 0.4748$, FDR $q = 0.5877$).
  - Classical Monocytes: Display flat noise metrics across both schemes (without-replacement $\beta = +0.000069/\text{mo}, p = 0.5142$).
  - Conclusion: Unadjusted single-cell transcriptional noise increases reflect sampling variance rather than biological entropy accumulation.

### 1.6. Human Composition Study Sensitivity (D3)
* **Ground Truth Output:** In [`results/d3-human-composition-study-sensitivity.csv`](../../results/d3-human-composition-study-sensitivity.csv) and [`results/d3-human-composition-loso.csv`](../../results/d3-human-composition-loso.csv):
  - Leave-One-Study-Out (LOSO) across all major study centers confirms club cell decline remains significant across all dropped studies ($p \in [0.0002, 0.0202]$, slopes $-0.0780$ to $-0.0446$).
  - However, full study dummy adjustment attenuates significance ($p = 0.2415$) due to structural collinearity between donor age and study center in `Lung_Integrated`.

### 1.7. Benchmark Reproduction of Pan-Mammalian Epigenetic Clock 2
* **Ground Truth Output:** In [`results/h1-epigenetic-clock-human.csv`](../../results/h1-epigenetic-clock-human.csv) and [`results/h1-epigenetic-clock-mouse.csv`](../../results/h1-epigenetic-clock-mouse.csv):
  - Double-exponential inverse link function:
    $$\text{Predicted Relative Age} = \exp\left(-\exp\left(-\hat{\eta}\right)\right)$$
    $$\text{Predicted Chronological Age} = (\text{Predicted Relative Age} \cdot (\text{Max Lifespan} + \text{Gestation})) - \text{Gestation}$$
    where $\text{Gestation} = 0.76\text{ yr}$ for human and $0.05\text{ yr}$ for mouse; $\text{MaxLifespan} = 122.5\text{ yr}$ for human and $4.0\text{ yr}$ for mouse.
  - Mouse Liver Lifespan (truncated, 724 / 816 CpGs): $r = 0.923$, $R^2 = 85.1\%$, $\text{MAE} = 0.36\text{ years}$ ($4.3\text{ months}$).
  - Human Multi-Tissue (consortium samples, 11 tissues): $r = 0.801$, $R^2 = 64.2\%$, $\text{MAE} = 5.06\text{ years}$ across all 661 samples (all 11 tissues accounted for).

---

## 2. Quantitative Verification Matrix

| Claim / Statistic | Source Script | Committed Result File | Verified Value | Status |
| :--- | :--- | :--- | :--- | :--- |
| **PsychAD Whole Cortex (Unadjusted)** | `run_m10_pipeline.py` | `results/m10-cohort-summary.csv` | $\beta = +0.000318/\text{yr}, p = 0.0118$ | **Exact Match** |
| **PsychAD Whole Cortex (Adjusted)** | `run_m10_pipeline.py` | `results/m10-cohort-summary.csv` | $\beta = +0.000027/\text{yr}, p = 0.8626$ | **Exact Match** |
| **PsychAD Excitatory Neurons (`EN_L3_5_IT_3` Pooled)** | `run_m10_pipeline.py` | `results/m10-celltype-stratified-summary.csv` | $\beta = -0.00140/\text{yr}, p = 3.79 \times 10^{-7}$ | **Exact Match** |
| **PsychAD Excitatory Neurons (`EN_L3_5_IT_3` Adult)** | `run_m10_pipeline.py` | `results/m10-celltype-stratified-summary.csv` | $\beta = -0.00060/\text{yr}, p = 0.0945$ | **Exact Match** |
| **PsychAD Microglia (Adjusted Pooled)** | `run_m10_pipeline.py` | `results/m10-celltype-stratified-summary.csv` | $\beta = +0.00142/\text{yr}, p = 0.0020$ | **Exact Match** |
| **PsychAD Microglia (Adult Stratum)** | `run_m10_pipeline.py` | `results/m10-celltype-stratified-summary.csv` | $\beta = +0.00200/\text{yr}, p = 0.0019$ | **Exact Match** |
| **PsychAD OPC (Adjusted Pooled)** | `run_m10_pipeline.py` | `results/m10-celltype-stratified-summary.csv` | $\beta = +0.00065/\text{yr}, p = 0.0486$ | **Exact Match** |
| **PsychAD OPC (Adult Stratum)** | `run_m10_pipeline.py` | `results/m10-celltype-stratified-summary.csv` | $\beta = +0.00129/\text{yr}, p = 0.0049$ | **Exact Match** |
| **Mouse Whole-Body Slope (REML)** | `run_m6_pipeline.py` | `results/m6-statistical-summary.csv` | $\beta = +0.2053, p = 0.00064$ | **Exact Match** |
| **Mouse Late-Adult Subset (21–30m)** | `run_m6_pipeline.py` | `results/m6-statistical-summary.csv` | $\beta = -0.0056, p = 0.9560$ | **Exact Match** |
| **NMR Skin Slope (Norm. Lifespan)** | `run_m6_pipeline.py` | `results/m6-statistical-summary.csv` | $\beta = +0.0058, \text{SE} = 0.0641, p = 0.9332$ | **Exact Match** |
| **NMR TOST Equivalence ($\delta=0.2023$)** | `run_m6_pipeline.py` | `results/m6-statistical-summary.csv` | $90\%\text{ CI } [-0.1450, +0.1567], p = 0.0274$ | **Exact Match** |
| **Direct Mouse-vs-NMR Contrast** | `run_m6_pipeline.py` | `results/m6-statistical-summary.csv` | $\Delta = +0.1285, p = 0.1153$ | **Exact Match** |
| **Mouse Lung Downsampled Noise (Without-Replacement)** | `run_d2_noise_sensitivity.py` | `results/d2-noise-sensitivity.csv` | $\beta = -0.000288/\text{mo}, p = 0.4748$ | **Exact Match** |
| **Mouse Lung Fibroblast *Cdkn1a* (p21)** | `run_fibroblast_markers.py` | `results/fibroblast-markers-audit.csv` | $\beta = +0.1080/\text{mo}, q = 0.0322$ | **Exact Match** |
| **Mouse Lung Donors in Composition** | `run_h2_composition.py` | `results/h2-composition-mouse-lung.csv` | $N = 16$ mice | **Exact Match** |
| **Universal Clock 2 Mouse Liver MAE** | `run_h1_epigenetic_drift.py` | `results/h1-epigenetic-clock-mouse.csv` | $\text{MAE} = 0.36\text{ yr}, r = 0.923$ | **Exact Match** |
| **Universal Clock 2 Human MAE (11 tissues)** | `run_h1_epigenetic_drift.py` | `results/h1-epigenetic-clock-human.csv` | $\text{MAE} = 5.06\text{ yr}, r = 0.801$ | **Exact Match** |

---

## 3. Automated Test Verification Results

All 16 unit, integration, and security tests pass cleanly:

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
