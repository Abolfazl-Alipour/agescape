# Implementation Plan: Testing Three Alternative Aging Hypotheses
## Moving Beyond the Mean-Expression Gene Clock: Composition, Noise, and Epigenetics

---

### Overview & Objectives
Our previous findings proved that a static 200-gene mean-expression signature derived from mouse does not serve as a universal cross-species aging clock. 

To determine what *does* change with age, we investigate three alternative hypotheses:
1. **Phase 1 (Low-Hanging Fruit / Easiest): Cell-Type Composition Shifts** — Does organ aging manifest as tissue remodeling (fibrosis, immune infiltration, parenchymal loss) across donors?
2. **Phase 2 (Stochastic Degradation): Cellular Entropy & Transcriptional Noise** — Does aging manifest as an explosion of within-cell-type cell-to-cell variance (loss of transcriptional coordination) rather than a shift in average expression?
3. **Phase 3 (Data Acquisition Protocol): Epigenetic Drift** — A concrete plan to download and analyze public DNA methylation arrays to test epigenetic clock erosion.

All methods incorporate the strict constraints established by the **10-Round Brutally Honest Scientific Crucible Review** ([crucible_review_three_hypotheses.md](crucible_review_three_hypotheses.md)).

---

## Phase 1: Cell-Type Compositional Remodeling (Hypothesis 2)

### 1. Data Inputs
* **Mouse Arm:** `data/raw/tms_droplet_raw.h5ad`
  - Focus organ: **Lung** (16 mice: 1m, 3m, 18m, 21m, 30m; both sexes; 24,540 cells) and **Kidney** (18 mice, 21,647 cells).
* **Human Arm:**
  - **Lung:** `data/raw/TS_Lung.h5ad` and `results/m9-donor-scores.csv` (`Lung_Integrated` cohort, 137 control donors with numeric ages 2–81).
  - **Brain:** `data/external/dlpfc-atlas/SupplementaryData1_Jan2026.csv` / `lifespan_pseudobulk.rds` (284 donors, age 0.2–97).

### 2. Mathematical Methodology (Crucible Rounds 1 & 2)
1. **Donor-by-CellType Matrix:** For each donor $i$, compute the count vector of captured cells across cell types $k \in \{1, \dots, K\}$: $\mathbf{n}_i = (n_{i1}, \dots, n_{iK})$.
2. **Aitchison Centred Log-Ratio (CLR) Transformation:**
   To overcome the closed-simplex negative correlative bias:
   $$\text{CLR}(p_{ik}) = \ln\left(\frac{p_{ik}}{g(\mathbf{p}_i)}\right), \quad \text{where } g(\mathbf{p}_i) = \left(\prod_{j=1}^K p_{ij}\right)^{1/K}$$
   with a Laplace pseudocount of 0.5 for zero cell counts.
3. **Statistical Modeling:**
   Regress CLR coordinates against age with sex adjustment:
   $$\text{CLR}(p_{ik}) = \beta_{0,k} + \beta_{\text{age},k} \cdot \text{Age}_i + \beta_{\text{sex},k} \cdot \text{Sex}_i + \epsilon_{ik}$$
   Degrees of freedom = number of donors ($N$).
4. **Multiple Testing Correction:** Benjamini-Hochberg FDR correction across the $K$ cell types in each organ.

### 3. Pre-Registered Falsification Criteria (Crucible Round 9)
* **Within-Species Falsification:** If zero cell types show significant age associations ($q < 0.05$ and $|\beta_{\text{age}}| > 0.005$) in mouse or human, tissue composition remodeling is declared a false hypothesis in that organ.
* **Cross-Species Conservation Test:** For homologous cell types mapped via Cell Ontology (e.g., monocytes/macrophages, fibroblasts, epithelial AT2):
  - If mouse and human show discordant signs (e.g., monocytes expand with age in mouse but contract or stay flat in human), **Conserved Mammalian Compositional Remodeling is formally FALSIFIED**.

### 4. Implementation Steps
* Create `scripts/run_h2_composition.py`.
* Output: `results/h2-composition-mouse-lung.csv`, `results/h2-composition-human-lung.csv`, `results/h2-composition-cross-species-concordance.csv`.
* Generate visualization: `figures/h2-composition-shifts-mouse-vs-human.png`.

---

## Phase 2: Cellular Entropy & Transcriptional Noise (Hypothesis 3)

### 1. Data Inputs
* **Mouse:** `data/raw/tms_droplet_raw.h5ad` (Lung: alveolar type 2 cells, classical monocytes, capillary endothelial cells).
* **Human:** `data/raw/TS_Lung.h5ad` (pulmonary alveolar type 2 cells, macrophages, capillary endothelial cells).

### 2. Biophysical & Mathematical Methodology (Crucible Rounds 3, 4, 7, 8)
1. **Homogeneous Cell-Type Subsetting:** Restrict analysis to clean, well-defined cell ontology types (no mixed broad classes).
2. **UMI Rarefaction (Downsampling Control):**
   - For all cells belonging to target cell type $k$, downsample raw counts to an identical library size of **$U = 1,500$ UMIs per cell** using multinomial sampling.
   - Discard cells with $< 1,500$ UMIs.
   - Enforce donor filter: Donor must retain at least **$n \ge 30$ rarefied cells**.
3. **Distance-to-Median (DM) Metric (Martinez-Jimenez et al., *Science* 2017):**
   - Within each donor $i$, compute gene mean expression $\mu_{ig}$ and squared coefficient of variation $CV^2_{ig}$.
   - Fit a robust running median through the $\ln(CV^2)$ vs $\ln(\mu)$ relationship across all expressed genes to establish technical baseline variability.
   - For each gene $g$, calculate the vertical distance: $\text{DM}_{ig} = \ln(CV^2_{ig}) - \text{baseline}(\mu_{ig})$.
   - The donor's transcriptional noise metric is the median DM across genes: $\overline{\text{DM}}_i$.
4. **Multivariate Dispersion (`betadisper`):**
   - Compute the centroid of cells in PCA space for each donor.
   - Measure each cell's Euclidean distance to the donor centroid. Compute mean distance $D_i$.
5. **Donor-Level Regression:**
   $$\overline{\text{DM}}_i = \gamma_0 + \gamma_{\text{age}} \cdot \text{Age}_i + \gamma_{\text{sex}} \cdot \text{Sex}_i + \eta_i$$

### 3. Pre-Registered Falsification Criteria (Crucible Round 10)
* **Entropy Falsification:** If $\gamma_{\text{age}} \le 0$ or $p \ge 0.05$ after UMI downsampling, the hypothesis that "aging increases cell-to-cell transcriptional noise/entropy" is **FALSIFIED** for that cell type.
* **Cross-Species Conservation Test:** If aged mouse cells become noisy ($\gamma > 0$) but aged human cells of the homologous type remain tightly regulated ($\gamma \le 0$), **Universal Mammalian Transcriptional Entropy is FALSIFIED**.

### 4. Implementation Steps
* Create `scripts/run_h3_entropy.py`.
* Output: `results/h3-entropy-mouse.csv`, `results/h3-entropy-human.csv`, `results/h3-entropy-comparison.csv`.
* Generate visualization: `figures/h3-transcriptional-noise-vs-age.png`.

---

## Phase 3: Epigenetic Drift & DNA Methylation (Hypothesis 1)

### 1. Rationale & Data Sourcing Plan
Because single-cell RNA-seq captures only transient transcriptional output, testing physical epigenetic drift requires DNA methylation profiling.

* **Target Public Dataset:**
  - **Mammalian Methylation Consortium:** GEO accession `GSE223748` / `GSE190666` (Infinium Mammalian Methylation 320k array; includes matched mouse and human tissues profiled on a shared orthologous CpG platform).
  - **Human Reference Cohort:** GEO `GSE40279` (Hannum et al., 656 human donors, whole blood, Illumina 450K) or `GSE115867` (Human brain DLPFC).

### 2. Analysis Protocol (When Data is Fetched)
1. Extract $\beta$-values ($\beta = M / (M + U + 100)$) for shared orthologous CpGs.
2. Compute **Shannon Epigenetic Entropy**:
   $$H_i = -\frac{1}{N_{\text{CpG}}} \sum_{j=1}^{N_{\text{CpG}}} \left[ \beta_{ij} \log_2(\beta_{ij}) + (1 - \beta_{ij}) \log_2(1 - \beta_{ij}) \right]$$
3. Test whether global epigenetic entropy $H_i$ increases with age at conserved rates across species.
4. Download script: `scripts/fetch_epigenetic_data.py`.

---

## Deliverables & Milestones Summary

| Phase | Hypothesis | Target Script | Primary Metric | Primary Output CSV | Status |
|---|---|---|---|---|---|
| **Phase 1** | H2: Composition Shifts | `scripts/run_h2_composition.py` | Aitchison CLR $\beta_{\text{age}}$ | `results/h2-composition-cross-species.csv` | Ready to execute |
| **Phase 2** | H3: Cellular Entropy | `scripts/run_h3_entropy.py` | Rarefied Distance-to-Median (DM) | `results/h3-entropy-comparison.csv` | Ready to execute |
| **Phase 3** | H1: Epigenetic Drift | `scripts/fetch_epigenetic_data.py` | CpG Shannon Entropy | `results/h1-epigenetic-entropy.csv` | Planned / Spec ready |

---

### Verification and Test Integration
* Add `tests/test_hypotheses_2_and_3.py` to the automated pytest suite.
* Ensure zero hardcoded paths, real data only, strict donor-level degrees of freedom, and automated verification of all generated results.
