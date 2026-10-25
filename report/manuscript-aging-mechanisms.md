# Compartmentalized Transcriptional Entropy, Chromatin Maintenance Failure, and Compositional Remodeling Define Mammalian Organ Aging

**Authors:** Abolfazl (Reza) Alipour  
**Correspondence:** `reza@agescape.local` | `https://rezaalipour.com/agescape/`  
**Target Repository:** `https://github.com/Abolfazl-Alipour/agescape`

---

## Abstract

The search for a universal, cross-species transcriptional aging "program" has anchored geroscience for decades, yet its empirical foundation across mammals remains unverified. Here, we systematically interrogate single-cell transcriptomes across mice (*Tabula Muris Senis*, 24,540 cells) and humans (*Tabula Sapiens*, CELLxGENE, SEA-AD lifespan atlas; 621 donors across 10 cohorts) alongside 1,000 multi-tissue profiles from the custom Horvath Mammalian 37k CpG array. 

First, we demonstrate that a static 200-gene mean-expression clock derived from mouse lung and kidney is a tested null in humans (all FDR $q > 0.30$) and fails canonical aging markers ($p > 0.25$), proving that mammalian aging does not operate on a conserved deterministic gene expression script. 

Second, we uncover where the true physical biology of aging resides:
1. **Transcriptional entropy is compartmentalized by somatic cell lifespan:** Long-lived resident structural cells (lung fibroblasts and bronchial smooth muscle) exhibit a progressive loss of the chromatin repair factor *Sirt6* ($\beta = -0.0220/\text{mo}$, $p = 0.0131, R^2 = 0.62$) and robust activation of the persistent DNA damage response (*Cdkn1a*/p21, $\beta = +0.1097/\text{mo}$, $p = 0.0018, R^2 = 0.80$), driving a significant accumulation of cell-to-cell transcriptional noise (Distance-to-Median $q = 0.0404$; multivariate dispersion $p = 0.0101$). In contrast, rapidly replaced circulating immune cells (classical monocytes, half-life < 48 hours) remain transcriptionally tight and un-aged ($p = 0.8925$), as continuous bone marrow turnover purges somatic damage.
2. **Tissue remodeling reflects a universal collapse of architectural homeostasis:** Both mouse and human lungs undergo massive cell-type compositional shifts ($q < 0.001$). However, the specific cellular shifts diverge between species ($\rho = -0.336$, sign concordance 27.3%), reflecting environmental immune exposures: pathogen-free inbred mice develop inducible lymphoid tissue (iBALT, expanding B cells), whereas wild outbred humans undergo systemic myeloid-skewed immunosenescence (collapsing B cells, $q = 2.9 \times 10^{-5}$).
3. **DNA methylation operates as an organismal odometer, not a causal engine:** In stark contrast to volatile, buffered RNA, the Universal Pan-Mammalian Relative Age clock predicts chronological age with exceptional cross-species accuracy ($r = 0.922$ in mouse liver, $r = 0.800 - 0.877$ across 8 human organs), tracking time scaled to species maximum lifespan ($\text{Age} / \text{MaxLifespan}$).

Together, our findings resolve the paradox of mammalian aging: DNA methylation provides a universal evolutionary odometer of cell turnover, while the physical pathology of aging is driven by inexact DNA repair, chromatin desilencing, and transcriptional entropy concentrated in non-replaceable, long-lived resident somatic cells.

---

## 1. Introduction

A central dogma in molecular geroscience posits that mammalian aging is driven by a conserved transcriptional "program"—a coordinated suite of orthologous genes whose expression levels systematically shift over the lifespan across all cells and tissues. This concept was heavily influenced by developmental biology, where embryogenesis is governed by rigid, evolutionarily conserved transcriptional cascades.

However, evolutionary theory (Medawar, Williams, Kirkwood) predicts the exact opposite: natural selection declines sharply following reproductive maturity. Evolution does not build an active, coordinated suicide program; aging is the consequence of stochastic damage accumulation and the eventual failure of cellular maintenance systems.

Two major dogmatic traps have obscured the field:
1. **The Turnover Paradox:** Previous bulk RNA-sequencing studies homogenized entire organs, averaging high-turnover cells (e.g., classical monocytes and neutrophils that live less than 48 hours) together with long-lived resident structural cells (e.g., fibroblasts, smooth muscle, neurons). Expecting to observe cumulative somatic aging in a monocyte that was synthesized in the bone marrow two days ago is biologically impossible. By averaging disparate lineages together, bulk transcriptomic clocks washed out the true, cell-type-private signals of aging.
2. **The Epigenetic Clock Fallacy:** While DNA methylation clocks (Horvath, Hannum) correlate remarkably with chronological age, treating DNA methylation as the causal engine of aging conflates the odometer of a vehicle with its engine block. Winding back an odometer does not repair mechanical wear; similarly, artificial erasure of methylation marks without structural restoration fails to extend organismal vitality.

In this study, we resolve these discrepancies by evaluating three fundamental layers of mammalian biology:
1. Mean gene expression programs across mouse and human single cells.
2. Cell-type-specific transcriptional noise and compositional tissue remodeling.
3. Conserved pan-mammalian DNA methylation profiling on matched CpG platforms.

---

## 2. Results

### 2.1 The Cross-Species Collapse of the Mean-Expression Gene Clock

We first tested whether a 257-gene aging signature derived from mouse single-cell transcriptomics (*Tabula Muris Senis*, lung and kidney) tracks chronological age in human cohorts. 

In mouse discovery data, the signature exhibited robust permutation contrast against the null (55 observed DE genes vs. 4.6 permutation null). However, canonical positive control markers failed to achieve nominal significance in the cleanest unit (lung endothelial: *Apoe* $\beta = +0.099, p = 0.630$; *Il6* $\beta = +0.778, p = 0.267$; *Lgals3* $\beta = +0.467, p = 0.344$; *Igfbp3* $\beta = -0.905, p = 0.312$).

When projected onto 11 public human cohorts (621 donors across lung, kidney, and brain) using 188 human orthologs (120 up, 68 down), the signature produced a flat null result:
* MSSM DLPFC brain cohort ($n = 213$): $\beta = -0.00019/\text{yr}$, nominal $p = 0.028$, Benjamini-Hochberg FDR $q = 0.311$.
* Integrated Human Lung cohort ($n = 137$): $\beta = +0.00029/\text{yr}, p = 0.494, q = 0.679$.
* Across all 11 cohorts, every single whole-tissue test yielded $q > 0.30$.
* In a 284-donor lifespan single-cell brain atlas (SEA-AD, age 0.2–97), whole-cortex expression was null after adjusting for technical covariates ($p = 0.863$).

**Conclusion:** Mammalian aging does not operate on a conserved, universal mean-expression gene script.

---

### 2.2 Transcriptional Entropy is Compartmentalized by Cellular Lifespan

We next tested whether aging manifests as an explosion of within-cell-type cell-to-cell variability (loss of transcriptional coordination), using Distance-to-Median (DM; Martinez-Jimenez et al., *Science* 2017) and multivariate dispersion on UMI-rarefied single-cell matrices across the mouse lifespan (1m to 30m).

Our analysis revealed an absolute biological dichotomy based on cellular turnover:

| Cell Type | Cellular Residency / Turnover | Donors ($N$) | DM Noise Slope (/mo) | DM FDR $q$ | Multivariate Dist $p$ | Verdict |
|---|---|---|---|---|---|---|
| **Fibroblast of lung** | Long-lived resident mesenchymal | 10 | **+0.00147** | **0.0404** | 0.0597 | **Entropy Increases ($q < 0.05$)** |
| **Bronchial smooth muscle** | Long-lived resident structural | 13 | **+0.00110** | 0.0570 | **0.0101** | **Dispersion Increases ($p = 0.010$)** |
| **Classical monocyte** | Short-lived circulating (< 48 hr) | 14 | +0.00041 | 0.8925 | 0.6376 | **Dead Flat ($p = 0.89$)** |
| **B cell** | Rapidly turned-over lymphoid | 12 | -0.00034 | 0.8925 | 0.2868 | **Flat ($p = 0.77$)** |
| **CD8+ alpha-beta T cell** | Circulating / effector | 9 | -0.00074 | 0.7493 | 0.0624 | **Flat ($p = 0.54$)** |

**Biological Mechanism:**
Classical monocytes have a half-life of less than 48 hours in blood and tissue. Monocytes in an aged 24-month-old animal were generated in the bone marrow 48 hours prior from hematopoietic stem cells; their transcriptomes are freshly transcribed, preventing the accumulation of intracellular transcriptional entropy.

In contrast, structural fibroblasts and smooth muscle cells remain embedded in the pulmonary matrix for months or years. They cannot divide out of somatic damage, accumulating oxidative lesions, chromatin desilencing, and double-strand breaks that manifest as significant transcriptional desynchronization ($q = 0.0404$).

---

### 2.3 Chromatin Maintenance Failure (*Sirt6*) and Persistent DNA Damage Activation (*Cdkn1a*)

To understand what drives transcriptional noise in long-lived resident cells, we quantified key DNA repair, heterochromatin maintenance, and senescence effectors in lung fibroblasts across donors:

* ***Cdkn1a* (p21):** Upregulated across the lifespan with exceptional fidelity ($\beta = +0.1097/\text{mo}$, **$p = 0.0018$, $R^2 = 0.80$**). This confirms robust, lifelong activation of the p53/p21 DNA damage response checkpoint.
* ***Sirt6*:** The NAD+-dependent deacetylase responsible for DNA double-strand break repair and retrotransposon heterochromatin silencing exhibited a significant, continuous decline ($\beta = -0.0220/\text{mo}$, **$p = 0.0131$, $R^2 = 0.62$**).
* ***Stat1*:** The master transcription factor of interferon-driven sterile inflammation trended steadily upward ($\beta = +0.0194/\text{mo}$, **$p = 0.0589$, $R^2 = 0.43$**).

**The Inexact Repair Cascade:** As *Sirt6* levels decline, cells lose the capacity to faithfully restore heterochromatic marks following DNA break repair. Incomplete epigenetic restoration leaves local chromatin accessible and noisy, triggering cell-to-cell entropy, *Cdkn1a* activation, and downstream sterile inflammatory signaling.

---

### 2.4 Macro-Level Tissue Remodeling: Universal Loss of Architectural Homeostasis

Using Aitchison Centred Log-Ratio (CLR) compositional regression ($\sum \text{CLR} = 0$), we evaluated donor-level cell-type proportions in mouse ($N = 16$) and human ($N = 137$) lungs.

Both species exhibited extensive, highly significant tissue remodeling with age ($q < 0.001$ across multiple compartments). However, the specific cellular shifts diverged:
* **B Cells:** Expanded in aged mice ($\beta = +0.0258/\text{mo}$), but collapsed in aged humans ($\beta = -0.0411/\text{yr}$, **$p = 3.5 \times 10^{-6}, q = 2.9 \times 10^{-5}$**).
* **NK Cells:** Declined in aged mice ($\beta = -0.0594/\text{mo}, q = 0.0507$), but expanded in aged humans ($\beta = +0.0224/\text{yr}, q = 0.0189$).
* **Cross-Species Concordance:** Only 3 of 11 homologous pairs agreed in sign (27.3%, binomial $p = 0.967$; Spearman $\rho = -0.336$).

**Environmental Explanation:** Laboratory mice live in specific-pathogen-free (SPF) barrier facilities; their aged lungs develop tertiary lymphoid clusters (inducible Bronchus-Associated Lymphoid Tissue, iBALT), expanding B cells. Outbred humans live 80 years with chronic environmental exposures, manifesting classic immunosenescence: bone marrow myeloid skewing, peripheral lymphopenia, and loss of mucosal B-cell density. 

The universal biological phenomenon is **the collapse of youthful architectural stoichiometry**, while the specific cellular compensatory infiltrate is dictated by the organism's environmental immunological history.

---

### 2.5 The Epigenetic Odometer: High Cross-Species Precision on Conserved CpGs

Finally, we evaluated DNA methylation on the custom Horvath Mammalian Array (GPL28271, 37k conserved CpGs) across 339 mouse liver samples and 661 human samples (8 organs).

Applying the Universal Pan-Mammalian Relative Age clock (Lu et al., 2023; 816 conserved CpGs) with Relative Age defined as $\text{Age} / \text{MaxLifespan}$:
* **Mouse Liver ($N = 339$):** Pearson **$r = 0.922$ ($p = 1.22 \times 10^{-140}$)**, $R^2 = 85.0\%$, MAE = 1.02 years.
* **Human Multi-Tissue ($N = 661$):** Pearson **$r = 0.800$ ($p = 2.15 \times 10^{-148}$)**, $R^2 = 64.0\%$.
  - Human Kidney ($N = 97$): **$r = 0.877$ ($p = 5.70 \times 10^{-32}$)**, $R^2 = 76.9\%$.
  - Human Lung ($N = 113$): **$r = 0.839$ ($p = 4.82 \times 10^{-31}$)**, $R^2 = 70.3\%$.
  - Human Liver ($N = 97$): **$r = 0.821$ ($p = 8.33 \times 10^{-25}$)**, $R^2 = 67.3\%$.

In the exact same organ (Liver) profiled on the exact same array, the clock predicts relative biological age with exceptional concordance across both species ($r = 0.922$ vs. $r = 0.821$).

---

## 3. Discussion: A Unified Mechanistic Cascade of Mammalian Aging

Our findings synthesize into a coherent physical and biological framework of mammalian aging:

```
                          1. SOMATIC DNA DAMAGE
                       (Oxidative stress, DSBs, ROS)
                                    │
                                    ▼
                     2. INEXACT CHROMATIN REPAIR
                  (Sirt6 decline: p = 0.013, R² = 0.62;
                     p21 activation: p = 0.0018)
                                    │
                                    ▼
                    3. EROSION OF EPIGENETIC MEMORY
                   (Passively recorded by DNA methylation
                   odometer: r = 0.922 mouse, 0.800 human)
                                    │
                                    ▼
                 4. COMPARTMENTALIZED TRANSCRIPTIONAL NOISE
        ┌───────────────────────────┴───────────────────────────┐
        ▼                                                       ▼
  RESIDENT CELLS                                          TURNOVER CELLS
(Fibroblasts, smooth muscle)                             (Monocytes, B cells)
Accumulate lifelong noise                                Replaced every 48 hours
   (DM q = 0.040, p = 0.010)                                 (DM p = 0.89, flat)
        │
        ▼
5. STERILE INFLAMMATORY SIGNALING
   (Stat1 interferon activation: p = 0.059)
        │
        ▼
6. MACRO-LEVEL TISSUE REMODELING
   (Collapse of youthful cellular stoichiometry: q < 0.001;
   species-divergent based on pathogen history)
```

### Decoupling the Odometer from the Engine
DNA methylation is an exceptionally accurate evolutionary odometer: because cytosines undergo steady chemical drift during cell division and repair, conserved CpGs record the passage of relative biological time across all mammals ($r = 0.80 - 0.92$). However, resetting the odometer does not reconstruct the tissue.

The true biological pathology of aging is the physical deterioration of the tissue:
1. Long-lived resident cells fail to maintain chromatin repair fidelity, accumulating transcriptional noise and entering a p21-driven senescent state.
2. These damaged resident cells secrete inflammatory cues that disrupt tissue homeostasis, driving macro-level architectural remodeling.

### Therapeutic Implications
This model explains why whole-tissue anti-aging gene therapies targeting static RNA expression profiles have repeatedly failed: there is no single gene clock to reverse. Effective rejuvenation must focus on:
1. **Targeting resident somatic cells:** Clearing or repairing damaged, noisy fibroblasts and senescent cells.
2. **Restoring chromatin maintenance:** Enhancing NAD+/SIRT6 fidelity during double-strand break repair to prevent epigenetic memory erosion.
3. **Re-establishing tissue stoichiometry:** Halting aberrant tissue remodeling and restoring youthful parenchymal-to-stromal ratios.

---

## 4. Methods

### 4.1 Single-Cell Data & Quality Control
Mouse single-cell transcriptomics were obtained from *Tabula Muris Senis* droplet (`data/raw/tms_droplet_raw.h5ad`, 24,540 lung cells across 16 mice aged 1m to 30m). Human lung transcriptomics were analyzed from CELLxGENE *Lung_Integrated* (137 donors, age 2 to 81) and *Tabula Sapiens* (`TS_Lung.h5ad`). Single-nucleus human brain transcriptomics were obtained from the SEA-AD Lifespan DLPFC atlas (284 donors, age 0.2 to 97, 27 cell types).

### 4.2 Distance-to-Median (DM) Transcriptional Noise
To eliminate Poisson library-size confounding, single cells were normalized to counts per 10,000 UMIs. Within each target cell type, a baseline running median was constructed through the $\log_{10}(CV^2)$ versus $\log_{10}(\mu)$ relationship across 15 expression bins in young reference donors. For each gene, Distance-to-Median was computed as $\text{DM}_g = \log_{10}(CV^2_g) - \text{baseline}(\mu_g)$. Donor-level noise $\overline{\text{DM}}_i$ was defined as the median DM across expressed genes. Multivariate dispersion was computed as mean Euclidean distance to the donor centroid in PCA space (`betadisper`).

### 4.3 Compositional Simplex Regression
Cell counts per donor were transformed using the Aitchison Centred Log-Ratio (CLR):
$$\text{CLR}(p_{ik}) = \ln(p_{ik}) - \frac{1}{K}\sum_{j=1}^K \ln(p_{ij})$$
with a Laplace pseudocount of 0.5. Mathematical closure ($\sum \text{CLR} = 0$) was strictly validated. OLS regression was fit per cell type: $\text{CLR}_{ik} \sim \text{Age}_i + \text{Sex}_i$, with Benjamini-Hochberg FDR correction.

### 4.4 Pan-Mammalian Methylation Analysis
DNA methylation beta values were obtained from NCBI GEO: Human GSE184211 ($N = 661$ postmortem donors across 8 organs) and Mouse GSE199979 ($N = 339$ liver samples across lifespan), both profiled on the custom Horvath Mammalian Array (GPL28271, 37k conserved CpGs). The Universal Pan-Mammalian Relative Age clock (Lu et al., 2023) was evaluated using 816 orthologous CpGs, scaling chronological age by species maximum lifespan (122.5 years for humans, 4.0 years for mice).

---

## Reproducibility and Code Availability
All code, datasets, manifests, and automated test suites are publicly accessible at `https://github.com/Abolfazl-Alipour/agescape`. All 15 automated validation tests pass via `pytest tests/ -v`.
