# 10-Round Brutally Honest Scientific Crucible Review
## Evaluation of Implementation Plan for Testing Three Alternative Aging Hypotheses:
1. **Hypothesis 2 (Phase 1, Easiest): Cell-Type Compositional Remodeling**
2. **Hypothesis 3 (Phase 2): Cellular Entropy and Transcriptional Noise**
3. **Hypothesis 1 (Phase 3): Epigenetic Drift (Data Acquisition Protocol)**

---

### Executive Review Summary
The project proposes to pivot from a discredited "universal 200-gene mean expression clock" to three fundamental mechanisms of aging:
1. Macro-level tissue remodeling (compositional shifts).
2. Micro-level stochastic loss of coordination (transcriptional entropy/noise).
3. Epigenetic drift (DNA methylation erosion).

This document subjects the proposed methodologies to a gauntlet of 10 brutally honest scientific review rounds across 5 specialized perspectives.

---

## Rounds 1 & 2: Compositional Data Analysis & Simplex Geometry Gauntlet
**Reviewer:** Mathematical Statistician specialized in Aitchison Simplex Geometry & Microbiome/Single-Cell Composition.

### Round 1: The Negative Correlative Bias Trap
* **The Critique:** In single-cell RNA-seq, the total number of recovered cells per donor is an arbitrary sequencing/loading constraint. If cell counts are converted into simple percentages ($p_k = n_k / \sum n_j$), the data lies on a closed simplex ($\sum p_k = 1$). Under simplex closure, if one dominant cell type expands (e.g., classical monocytes infiltrate the lung during chronic inflammation), all other cell types are mathematically forced to decrease in percentage, even if their absolute in vivo count remained identical! 
* **The Gauntlet Rule:** Standard OLS regressions of raw proportions against age are mathematically invalid and generate spurious negative correlations. The plan MUST employ:
  1. **Hierarchical Dirichlet-Multinomial Bayesian modeling (scCODA; Büttner et al., *Nature Comm* 2021)** using spike-and-slab logit-normal priors to shrink false compositional shifts to zero.
  2. Or **Centred Log-Ratio (CLR)** / Isometric Log-Ratio (ILR) transformations with a pre-specified, biologically invariant reference cell type.
* **Verdict:** REVISE. Ban raw proportion linear regression. Mandate scCODA / CLR with explicit reference cell selection.

### Round 2: Dissociation Bias vs. In Vivo Remodeling
* **The Critique:** Single-cell tissue dissociation (enzymatic digestion) does not harvest cells with equal probability. Older, fibrotic tissues with dense collagen require harsher digestion, preferentially lysing fragile cells (e.g., alveolar type 1 cells or neurons) while hardier cells (fibroblasts, monocytes) survive. A shift in captured cell percentages between a 3-month and 24-month mouse lung may reflect differential survival in the dissociation tube rather than an in vivo population change.
* **The Gauntlet Rule:** 
  1. The plan must cross-validate mouse droplet results (`tms_droplet_raw.h5ad`) against FACS/microfluidic controls.
  2. In humans, single-nucleus RNA-seq (snRNA-seq, where cell membranes are lysed and only intact nuclei are captured, eliminating dissociation survival bias) must be compared against single-cell (scRNA-seq) where available (e.g., Kidney sn vs scr, Brain snRNA-seq).
* **Verdict:** ACCEPT with mandatory snRNA-seq concordance check.

---

## Rounds 3 & 4: Single-Cell Biophysics & Information Theory Gauntlet
**Reviewer:** Biophysicist & Stochastic Gene Expression Auditor.

### Round 3: The Poisson Sampling Noise Illusion
* **The Critique:** In single-cell RNA-seq, over 80% of observed zero values and high variances are Poisson/negative-binomial sampling noise driven by technical capture efficiency (typically 5–15% mRNA capture). Older cells frequently show altered total RNA content, altered cell diameter, or higher fractions of mitochondrial reads (a hallmark of membrane permeabilization). If an old cell has 1,200 detected UMIs while a young cell has 3,500 UMIs, the old cell will mathematically appear "noisier" purely due to Poisson sampling variance!
* **The Gauntlet Rule:** 
  1. Raw pairwise Euclidean distances or uncorrected Shannon entropy are banned.
  2. The plan MUST perform **rarefaction / downsampling**: every single cell must be downsampled to an identical UMI threshold (e.g., exactly 1,500 UMIs per cell) before computing any cell-to-cell variability metric. Cells below the threshold are discarded.
* **Verdict:** MANDATORY. All entropy/noise metrics must run on UMI-rarefied count matrices.

### Round 4: Decoupling Variance from Mean Expression (The DM Method)
* **The Critique:** In transcriptomics, gene variance is an intrinsic non-linear function of mean abundance ($CV^2 \propto 1/\mu + \alpha$). If aging alters the mean expression of certain genes, their variance will change automatically as a mathematical consequence of the mean-variance relationship. Claiming "cells are becoming noisier" is only biologically meaningful if variance increases *beyond what is expected for that gene's expression level*.
* **The Gauntlet Rule:** 
  1. The plan must implement the **Distance to Median (DM)** method established by Martinez-Jimenez et al. (*Science*, 2017).
  2. Fit a running median through the $CV^2$ vs mean expression curve across all genes.
  3. Compute each gene's DM as its vertical distance from the running median baseline.
  4. Aggregate DM across cell types per donor to evaluate whether aged cells exhibit positive excess variance.
* **Verdict:** MANDATORY. Use DM and multivariate dispersion (`betadisper`) on downsampled cells.

---

## Rounds 5 & 6: Cross-Species Homology & Cell Ontology Gauntlet
**Reviewer:** Mammalian Evolutionary Biologist & Anatomical Ontologist.

### Round 5: The Cell Label Semantic Trap
* **The Critique:** Mouse annotators in *Tabula Muris Senis* used free-text labels like `"classical monocyte"` and `"bronchial smooth muscle cell"`, while human annotators in *Tabula Sapiens* and CELLxGENE used labels like `"CD14-positive monocyte"` and `"vascular associated smooth muscle cell"`. Simple string matching will fail or compare non-homologous cell types.
* **The Gauntlet Rule:** 
  1. The plan must establish an explicit, pre-mapped **Cell Ontology (CL)** concordance dictionary between mouse TMS and human TS.
  2. Only homologous cell pairs with shared CL terms (e.g., CL:0000583 for alveolar type 2, CL:0000576 for monocytes, CL:0000077 for mesothelial cells) can enter the cross-species conservation test.
* **Verdict:** MANDATORY. Build `data/reference/celltype-ontology-mapping.csv` before running cross-species comparisons.

### Round 6: Inbred Vivarium Mice vs. Outbred Human Lifestyle
* **The Critique:** C57BL/6 mice live in sterile, pathogen-free barrier facilities with identical diets, circadian cycles, and genetics. Their aging is dominated by specific inbred strain pathologies (e.g., histiocytic sarcomas, glomerulonephritis, lung alveolar proteinosis). Humans live 80 years with lifelong viral infections (EBV, CMV), divergent diets, smoking, and air pollution. If mouse lungs show monocyte infiltration while human lungs do not, that does NOT invalidate human aging biology—it proves that laboratory mouse vivarium aging cannot be assumed to reflect human environmental aging.
* **The Gauntlet Rule:** 
  - Cross-species results must be analyzed both **within-species first** (does mouse lung remodel with age? does human lung remodel with age?) and only then compared for **directional concordance**. Do not force a joint multi-species regression model that assumes identical baseline rates.
* **Verdict:** ACCEPT. Decouple within-species testing from cross-species concordance.

---

## Rounds 7 & 8: Statistical Power, Donor Replication & Confounding Gauntlet
**Reviewer:** Biostatistical Replicability & Pseudoreplication Auditor.

### Round 7: The N = Cells Pseudoreplication Death Penalty
* **The Critique:** If an investigator has 16 mice and 24,000 cells, running a statistical test treating the 24,000 cells as independent $N=24,000$ degrees of freedom is fraudulent pseudoreplication. Cells from the same mouse share a common microenvironment, batch, and genetics. $p$-values will be inflated by orders of magnitude ($p = 10^{-50}$ on pure noise).
* **The Gauntlet Rule:** 
  1. **Strict Donor-Level Degrees of Freedom:** All hypothesis testing MUST be executed at the donor/animal level ($N = 16$ mice in TMS lung, $N = 137$ human lung donors, $N = 284$ human brain donors).
  2. Metrics (compositional shares $C_{i,k}$, transcriptional dispersion $DM_{i,k}$) are calculated per donor $i$ and cell type $k$.
  3. Regression models regress donor-level metrics on donor age and sex: $\text{Metric}_{i,k} \sim \text{Age}_i + \text{Sex}_i$.
* **Verdict:** MANDATORY HARD RULE. Zero cell-level regression.

### Round 8: Minimum Donor Cell-Count Floor
* **The Critique:** If donor A only has 4 cells of type X captured, estimating cell-to-cell variance or distance-to-median from 4 cells is wildly unstable and will create extreme outliers that bias the donor-level regression.
* **The Gauntlet Rule:** 
  - Minimum floor: A donor must have at least **$n \ge 30$ cells** for a given cell type to compute dispersion metrics. If a donor has $< 30$ cells, that cell type is marked `NA` for that donor.
  - Furthermore, a cell type must be represented across at least **$n \ge 8$ donors** to enter the age regression.
* **Verdict:** ACCEPT. Enforce $n_{\text{cells}} \ge 30$ per donor and $n_{\text{donors}} \ge 8$ per cell type.

---

## Rounds 9 & 10: Falsification Crucible & Pre-Registered Decision Boundaries
**Reviewer:** Senior Scientific Falsification Referee & Meta-Science Auditor.

### Round 9: Explicit Falsification Criteria for Hypothesis 2 (Composition Shifts)
* **The Question:** Exactly what numerical result will count as **falsifying** the compositional remodeling hypothesis?
* **Pre-Registered Falsification Protocol:**
  1. Compute scCODA credible effects / CLR regressions across all eligible cell types in Mouse Lung and Human Lung.
  2. **Falsification Boundary 1 (Within-Species Null):** If zero cell types show FDR-adjusted $q < 0.05$ with an annual effect size $|\beta| > 0.01$, Hypothesis 2 is rejected within that organ.
  3. **Falsification Boundary 2 (Cross-Species Divergence):** If cell types show significant age shifts in mouse, but their homologous human cell types show either null effects ($q > 0.10$) or opposing signs, the hypothesis of **Conserved Mammalian Compositional Remodeling** is officially FALSIFIED.
* **Verdict:** ACCEPT. No moving goalposts.

### Round 10: Explicit Falsification Criteria for Hypothesis 3 (Cellular Entropy / Noise)
* **The Question:** Exactly what numerical result will count as **falsifying** the cellular entropy hypothesis?
* **Pre-Registered Falsification Protocol:**
  1. For each homologous cell type (e.g. Alveolar Type 2, Capillary Endothelial, Monocyte/Macrophage), calculate UMI-rarefied Distance-to-Median (DM) and PCA distance-to-median per donor.
  2. Regress donor DM against donor age: $\text{DM}_{i,k} \sim \text{Age}_i + \text{Sex}_i$.
  3. **Falsification Boundary 1 (Flat Dispersion):** If the slope $\beta_{\text{age}}$ is not significantly positive ($p \ge 0.05$, or FDR $q \ge 0.10$ across cell types), or if the slope is negative, the hypothesis that aging increases transcriptional noise is FALSIFIED for that cell type.
  4. **Falsification Boundary 2 (Technical Confounder Dominance):** If the correlation between dispersion and age is eliminated when controlling for sequencing depth / UMI count, the signal is declared a technical artifact, and the biological entropy hypothesis is FALSIFIED.
* **Verdict:** ACCEPT. Boundaries locked before code execution.

---

## Crucible Synthesis & Implementation Plan Directive
All 10 review rounds agree on a clear path forward:
1. **Phase 1 (Hypothesis 2, Easiest):** Analyze cell-type composition shifts using Aitchison simplex geometry (CLR / scCODA), donor-level counts, and strict CL ontology mapping across Mouse TMS and Human TS/Census.
2. **Phase 2 (Hypothesis 3, Noise/Entropy):** Analyze cell-to-cell variability using UMI-rarefied single cells, Distance-to-Median (DM) abundance-adjusted variance, and donor-level replication.
3. **Phase 3 (Hypothesis 1, Epigenetic Drift):** Formulate a concrete public data acquisition plan (GEO methylation arrays) with verified accession IDs and processing pipeline specs.
