# Pre-Registered Age Bins and Demographic Contrast Specification

**Binding specification for all downstream differential expression and conservation analyses.**
**Established:** Milestone M0 (before any DE computation exists).

---

## 1. Mouse Age-Bin Design (Tabula Muris Senis)

### 1.1 Verified Demographics in Droplet Dataset
From file metadata inspection of `data/raw/tms_droplet_raw.h5ad` (245,389 cells across 16 tissues and 20 mice):
- **1m:** 25,980 cells — **Male only**
- **3m:** 45,602 cells — **Mixed (Male + Female)**
- **18m:** 44,645 cells — **Mixed (Male + Female)**
- **21m:** 35,828 cells — **Female only**
- **24m:** 37,660 cells — **Male only**
- **30m:** 55,674 cells — **Male only**

### 1.2 Primary Contrast: 3m vs 18m
- **Design:** `~ age_group + sex` via pseudobulk `pydeseq2` per (mouse × tissue × cell type).
- **Rationale:** 3m vs 18m is the **only age pair with both sexes represented in both groups**. This eliminates the catastrophic confounding between age and sex present in other contrasts.
- **Sex chromosome safeguard:** All DE gene lists must be checked for sex-chromosome enrichment (Y-chromosome genes must not drive aging rankings).
- **Eligibility criteria:** Minimum ≥3 mice per age group and ≥50 cells per group for a given tissue × cell type unit. Ineligible units are formally listed in results, never dropped silently.

### 1.3 Sensitivity and Robustness Contrasts (Pre-registered)
All sensitivity runs must be reported side-by-side with the primary contrast; none may be pooled silently:
1. **Extended Bins:** Young = {1m, 3m}, Old = {18m, 21m, 24m, 30m}.
   - *Pre-registered Hypothesis:* Wider age gap strengthens effect sizes and increases detected DE genes.
   - *Crucial Caveat:* 30m mice represent a **survival-selected cohort** (hardy survivors). Any signal strengthening cannot be ascribed to chronological aging alone.
2. **Male-Only Contrast:** 18m vs 24m (male-only).
   - Serves to test whether age progression continues without sex imbalance confounds.
3. **Continuous Age Model (Bridge):** Pseudobulk `~ age_months + sex`.
   - Methodological bridge to replicate and compare against the Tabula Muris Senis paper's continuous linear model.

---

## 2. Human Age Handling (SenNet & Tabula Sapiens)

### 2.1 Principle: Continuous Covariate Only
- **Human age is NEVER dichotomized into arbitrary young/old bins.**
- Human aging in single-cell cohorts spans heterogeneous decades; grouping donors into binary bins creates severe threshold artifacts and throws away continuous statistical power.

### 2.2 Model Specification
- **Design:** Pseudobulk per (donor × cell type): `~ age_years + sex`.
- **Primary cohort:** SenNet scRNA-seq datasets meeting selection rule (largest age span with ≥6 donors ≥60 years, 10x platform).
- **Replication cohort:** Tabula Sapiens v1 (10x-only subset, backed mode per tissue).
- **Cohorts are strictly analyzed separately, NEVER merged**, to prevent batch × cohort × age confounding.

### 2.3 Inclusion and Power Thresholds
- To qualify for DE testing, a human cell type must possess:
  - **≥ 6 donors** with cell representation.
  - **Donor age standard deviation (SD) ≥ 12 years**.
- Cell types failing these inclusion criteria are placed into a flagged `low-powered` partition. They are **never** interpreted as biologically "divergent" from mouse.
- Power calculations (minimum detectable |log2FC| at FDR 0.1 given donor n) are computed and exported to `results/m4-human-power.csv`.

---

## 3. Log2 Fold Change Orientation Standard

Across all analyses in this repository:
$$\log_2 \text{FC} = \log_2 \left(\frac{\text{Expression}_{\text{Old}}}{\text{Expression}_{\text{Young}}}\right)$$
- **Positive $\log_2 \text{FC}$:** Increased expression with age (upregulated with age).
- **Negative $\log_2 \text{FC}$:** Decreased expression with age (downregulated with age).
- Notebook 02 will programmatically assert this orientation using known aging biomarkers (e.g. *Cdkn2a*).
