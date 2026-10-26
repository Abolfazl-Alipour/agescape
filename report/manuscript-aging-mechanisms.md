# Cell-Type-Private Transcriptional Signatures, Compositional Remodeling, and Epigenetic Clock Dynamics Across Mammalian Aging

**Authors:** Abolfazl (Reza) Alipour  
**Correspondence:** `reza@agescape.local` | `https://rezaalipour.com/agescape/`  
**Target Repository:** `https://github.com/Abolfazl-Alipour/agescape`

---

## Abstract

Whether mammalian aging is governed by a conserved, deterministic gene expression program across species remains a central question in geroscience. Here, we systematically interrogate single-cell transcriptomes across mice (*Mus musculus*, *Tabula Muris Senis*, 24,458 lung cells), humans (*Homo sapiens*, 621 postmortem donors across 10 cohorts in lung, kidney, and dorsolateral prefrontal cortex), and the long-lived, non-senescent naked mole-rat (*Heterocephalus glaber*, skin epidermis/dermis), alongside 1,000 pan-mammalian DNA methylation array profiles (Horvath Mammalian Array GPL28271).

First, we show that a 257-gene aging signature derived from mouse single-cell transcriptomics does not operate as a monolithic, whole-tissue program across species. In human cortex, whole-tissue pseudobulk associations attenuate to null after covariate adjustment ($\beta = +0.00015/\text{yr}, p = 0.863$). However, the signature is not universally absent: it transports across species in a **cell-type-private** manner, significantly tracking age in human oligodendrocyte precursor cells (OPCs) across four independent brain cohorts (SEA-AD $p = 1.1 \times 10^{-4}, q = 0.021$; MSSM $p = 0.008$; HBCC $p = 0.024$; RADC $p = 0.033$) and in human microglia across three cohorts (and M10 lifespan adult stratum $\beta = +0.00142/\text{yr}, p = 0.0019$).

Second, we deploy the naked mole-rat as an evolutionary negative control. In contrast to mice—where the signature rises robustly with age across 23 animals (whole-body slope $+0.134, p = 0.015$; REML cluster-robust $+0.205, p = 0.00064$)—the signature remains completely flat in naked mole-rat skin ($\beta = +0.0035, p = 0.925, n = 5$). Two-One-Sided-Test (TOST) equivalence confirms flatness within the mouse-derived margin ($p_{\text{TOST}} = 0.027$), and epidermal *Igfbp3* expression is unchanged (221.8 vs 223.8 CPM), verifying that the signature reflects true mammalian senescence kinetics rather than chronometric drift.

Third, using multinomial downsampling to an identical library size of 1,500 UMIs per cell, we evaluate transcriptional noise across nine mouse lung cell types. Long-lived lung fibroblasts exhibit a nominal age-associated increase in Distance-to-Median noise ($\beta = +0.00109/\text{mo}, p = 0.0443$, Benjamini-Hochberg $q = 0.271$), accompanied by significant upregulation of the cyclin-dependent kinase inhibitor *Cdkn1a* (p21; $\beta = +0.1080/\text{mo}, p = 0.0029, q = 0.0322$) and an exploratory decline in *Sirt6* ($\beta = -0.0136/\text{mo}, p = 0.0245, q = 0.1350$). Conversely, classical monocytes display stable noise metrics ($\beta = +0.00008/\text{mo}, p = 0.8923$), despite harboring pronounced mean differential expression.

Fourth, compositional modeling via Aitchison Centered Log-Ratio (CLR) demonstrates that tissue remodeling occurs extensively within species (6 mouse and 4 human adult cell types show FDR $q < 0.05$). However, cross-species concordance of cell-type trajectories is not statistically established (sign concordance 3/7, two-sided binomial $p = 1.000$, Spearman $\rho = -0.286, p = 0.535$), and human B-cell shifts reflect pediatric-to-adult immunological maturation rather than adult senescence.

Finally, we demonstrate that the Universal Pan-Mammalian Epigenetic Clock 2 (Lu et al., *Nature Aging* 2023), when evaluated with its published double-exponential inverse link function $\exp(-\exp(-\hat{y}))$, tracks chronological age with high precision across mammalian tissues ($\text{MAE} = 0.36\text{ years}$ [$4.3\text{ months}$], $r = 0.923$ in mouse liver; $\text{MAE} = 5.06\text{ years}$, $r = 0.801$ across 661 human postmortem tissues).

Our findings indicate that mammalian aging does not follow an identical, monolithic transcriptomic program across species. Instead, biological aging is characterized by cell-type-private transcriptional activation, species-divergent tissue remodeling, and a highly calibrated pan-mammalian epigenetic odometer.

---

## 1. Introduction

A persistent ambition of molecular geroscience is the discovery of a universal transcriptomic aging "program"—a conserved set of orthologous genes whose expression levels systematically track chronological age across mammalian tissues and species. Such a program would imply that somatic senescence is governed by evolutionary mechanisms analogous to developmental cascades.

However, evolutionary biology (Medawar's mutation accumulation, Williams' antagonistic pleiotropy, and Kirkwood's disposable soma theory) predicts the contrary: natural selection declines precipitously following reproductive maturation. Somatic aging is fundamentally stochastic, driven by the accumulation of molecular wear, somatic mutations, and imperfect repair.

When evaluating aging across mammals, four central biological questions emerge:
1. **Cell-Type-Private vs. Whole-Tissue Averages:** Does aging manifest as a universal shift in tissue-level mean expression, or do cell-type-specific responses dominate—meaning bulk tissue averaging dilutes and obscures true biological signals?
2. **The Non-Aging Negative Control:** If a transcriptional signature tracks genuine biological senescence, does it remain stable in mammalian species that exhibit negligible senescence, such as the naked mole-rat (*Heterocephalus glaber*)?
3. **Transcriptional Noise and Cellular Lifespan:** Does cell-to-cell transcriptional variability (cellular entropy) increase selectively in long-lived resident structural cells (such as fibroblasts) subject to cumulative damage, compared to short-lived circulating cells (such as monocytes) continuously replenished from stem-cell pools?
4. **Epigenetics vs. Transcriptomics:** Does DNA methylation provide a stable, pan-mammalian chronometer across species with vastly different lifespans, even when volatile RNA programs diverge?

Here, we integrate single-cell transcriptomics, comparative evolutionary biology, and pan-mammalian DNA methylation profiling to evaluate these hypotheses across mouse, human, and naked mole-rat cohorts.

---

## 2. Results

### 2.1 Cell-Type-Private Transport and Whole-Tissue Masking

We derived a 257-gene aging signature from mouse single-cell transcriptomics (*Tabula Muris Senis*, lung and kidney; Table 2.1). In the mouse discovery dataset, the signature demonstrated significant contrast against a permutation null (55 observed DE genes vs. 4.6 permutation null). Across 23 mice whole-body, the signature score rose robustly with age (slope $+0.134\text{ per month}, p = 0.015$).

```mermaid
flowchart LR
    A["Mouse Aging Signature (TMS Lung & Kidney)"] --> B["Whole-Cortex Bulk Averaging"]
    A --> C["Single-Cell Microglia"]
    A --> D["Single-Cell OPCs"]
    B --> E["Attenuates to Null (p = 0.863)"]
    C --> F["Replicates in 3/4 Cohorts + M10 Lifespan (p = 0.0019)"]
    D --> G["Replicates in 4/4 Cohorts (FDR q = 0.021)"]
```

To test cross-species transport, we projected 188 human 1:1 orthologs onto human single-nucleus RNA-seq datasets spanning 621 postmortem donors across 10 cohorts ([`results/m9-regressions.csv`](file:///home/reza/projects/agescape/results/m9-regressions.csv), [`results/m10-cohort-summary.csv`](file:///home/reza/projects/agescape/results/m10-cohort-summary.csv)):

* **Whole-Tissue Cortex Averaging Attenuates to Null:** In the SEA-AD DLPFC lifespan cohort (284 donors, ages 0.2–97 years), whole-cortex CPM scores showed an unadjusted slope of $+0.00049/\text{yr}$ ($p = 0.024$). However, after adjusting for brain bank, postmortem interval (PMI), and technical covariates, the association attenuated completely to null ($\beta = +0.00015/\text{yr}, p = 0.863$).
* **Cell-Type-Private Replications:** In stark contrast to tissue averages, cell-type-specific regressions revealed robust, reproducible age associations:
  * **Oligodendrocyte Precursor Cells (OPCs):** The mouse-derived signature increased significantly with age in human OPCs across all four independent DLPFC cohorts: SEA-AD ($p = 1.1 \times 10^{-4}$, surviving Benjamini-Hochberg across all 225 cell-type tests at $q = 0.021$), MSSM ($p = 0.008$), HBCC ($p = 0.024$), and RADC ($p = 0.033$).
  * **Microglia:** Human microglia exhibited significant positive age associations in three independent cohorts (SEA-AD $p = 0.005$, MSSM $p = 0.026$, RADC $p = 0.031$). In the M10 lifespan atlas, microglia maintained a significant covariate-adjusted slope of $+0.00142/\text{yr}$ ($p = 0.0020$; adult stratum $p = 0.0019$).

These results establish that cross-species transcriptional aging is not a uniform whole-tissue phenomenon; rather, it is cell-type private, concentrated in glia undergoing chronic microenvironmental stress.

---

### 2.2 The Naked Mole-Rat as an Evolutionary Negative Control

To determine whether the mouse-derived gene signature captures true biological senescence rather than chronometric drift or vivarium artifacts, we evaluated single-cell and pseudobulk RNA-seq from the naked mole-rat (*Heterocephalus glaber*, [`GSE152886`](file:///home/reza/projects/agescape/data/raw/GSE152886_nmr_skin_raw.csv.gz), Hilton et al. 2020 / Savina et al.). Naked mole-rats exhibit negligible senescence and live >37 years without accelerating mortality hazards.

Using 14,278 1:1 Ensembl orthologs ([`data/reference/orthologs-mouse-nmr.csv`](file:///home/reza/projects/agescape/data/reference/orthologs-mouse-nmr.csv)), we scored the signature in NMR skin across young (1.1 yr) and aged (11.3 yr) animals:

| Species / Cohort | Tissue | N | Aging Signature Slope | 95% CI | p-value |
|---|---|---|---|---|---|
| **Mouse (TMS Whole-Body)** | Multi-organ | 23 mice | **$+0.134$ / mo** | $[+0.029, +0.239]$ | $p = 0.015$ |
| **Mouse (REML Cluster-Robust)** | 19 clusters | 20 mice | **$+0.205$ / mo** | $[+0.116, +0.295]$ | $p = 0.00064$ |
| **Mouse (Holdout Tissues)** | 14 non-lung/kidney | 23 mice | **$+0.144$ / mo** | $[+0.038, +0.250]$ | $p = 0.0167$ |
| **Naked Mole-Rat (GSE152886)** | Skin | 5 animals | **$+0.0035$ / mo** | $[-0.093, +0.100]$ | **$p = 0.925$ (Flat)** |

Under a pre-specified mouse-derived margin ($\delta_{\min} = 0.202$), Two-One-Sided-Test (TOST) equivalence confirmed statistical equivalence ($p_{\text{TOST}} = 0.027$). Furthermore, in whole-epidermis pseudobulk, the canonical senescence marker *Igfbp3* was completely unchanged between young (221.8 CPM) and aged (223.8 CPM) animals. The flat trajectory in NMR supports the hypothesis that the signature reflects mammalian senescence kinetics rather than passive passage of time.

---

### 2.3 Single-Cell Transcriptional Noise from Rarefied Counts

To evaluate cell-to-cell transcriptional variability without library-depth confounding, we implemented exact multinomial downsampling to 1,500 UMIs per cell in mouse lung ([`scripts/run_h3_entropy.py`](file:///home/reza/projects/agescape/scripts/run_h3_entropy.py)). We evaluated Distance-to-Median (DM; Martinez-Jimenez et al., *Science* 2017) and multivariate dispersion across all qualifying cell types ([`results/h3-entropy-mouse-lung.csv`](file:///home/reza/projects/agescape/results/h3-entropy-mouse-lung.csv)):

| Cell Type | Donors | Rarefied Cells | DM Slope (/mo) | Raw p-val | FDR q-val | Dispersion Slope | Dispersion p-val |
|---|---|---|---|---|---|---|---|
| **fibroblast of lung** | 10 | 1,516 | **$+0.00109$** | **$0.0443$** | $0.2713$ | $+0.01764$ | $0.1017$ |
| **bronchial smooth muscle cell** | 13 | 2,336 | $+0.00028$ | $0.1038$ | $0.3113$ | $-0.00389$ | $0.3873$ |
| **alveolar macrophage** | 14 | 1,349 | $+0.00056$ | $0.4518$ | $0.6864$ | $-0.00289$ | $0.8430$ |
| **classical monocyte** | 14 | 7,913 | $+0.00008$ | $0.8923$ | $0.8923$ | $+0.00483$ | $0.5453$ |
| **intermediate monocyte** | 8 | 1,748 | $-0.00026$ | $0.8439$ | $0.8923$ | $-0.00396$ | $0.8100$ |
| **B cell** | 13 | 1,495 | $+0.00050$ | $0.3359$ | $0.6864$ | $+0.02509$ | $0.1357$ |
| **NK cell** | 15 | 1,190 | $-0.00083$ | $0.0603$ | $0.2713$ | $-0.01330$ | $0.2814$ |
| **CD8+ alpha-beta T cell** | 10 | 867 | $+0.00017$ | $0.5338$ | $0.6864$ | $+0.01606$ | $0.4848$ |
| **vein endothelial cell** | 8 | 298 | $+0.00060$ | $0.4787$ | $0.6864$ | $-0.00543$ | $0.6578$ |

* Lung fibroblasts displayed a nominal unadjusted increase in transcriptional noise ($p = 0.0443$), which attenuates to $q = 0.2713$ under Benjamini-Hochberg FDR correction across the nine cell types.
* Classical monocytes demonstrated flat noise metrics ($p = 0.8923$). However, monocytes cannot be categorized as "un-aged": in pseudobulk differential expression ([`results/m2-aging-signal-ranked.csv`](file:///home/reza/projects/agescape/results/m2-aging-signal-ranked.csv)), classical monocytes exhibit nine significant DE genes with a median $|\log_2\text{FC}|$ of $1.99$. Their stable variance coexists with substantial shifts in mean expression.

---

### 2.4 Targeted Audit of Senescence and DNA Repair Markers

To verify candidate mechanisms underlying fibroblast aging, we evaluated donor pseudobulk expression across 11 targeted markers of cellular senescence, DNA repair, and chromatin maintenance ([`scripts/run_fibroblast_markers.py`](file:///home/reza/projects/agescape/scripts/run_fibroblast_markers.py), [`results/fibroblast-markers-audit.csv`](file:///home/reza/projects/agescape/results/fibroblast-markers-audit.csv)):

| Gene | Pathway / Role | CP10K Mean | Slope (/month) | SE | p-value | FDR q-value | $R^2$ |
|---|---|---|---|---|---|---|---|
| **Cdkn1a** | Cyclin-dependent kinase inhibitor (p21) | 1.193 | **$+0.1080$** | 0.0242 | **$0.0029$** | **$0.0322$** | **$0.777$** |
| **Sirt6** | Chromatin ADP-ribosylase / deacetylase | 0.207 | **$-0.0136$** | 0.0048 | **$0.0245$** | $0.1350$ | $0.540$ |
| **Ccl2** | SASP chemokine | 0.459 | $+0.0481$ | 0.0200 | $0.0475$ | $0.1741$ | $0.498$ |
| **Stat1** | Interferon / sterile inflammation | 0.322 | $+0.0137$ | 0.0066 | $0.0765$ | $0.2104$ | $0.385$ |
| **Cxcl12** | Stromal chemokine | 1.090 | $-0.0363$ | 0.0236 | $0.1682$ | $0.3699$ | $0.278$ |
| **Parp1** | DNA single-strand break repair | 0.158 | $-0.0046$ | 0.0036 | $0.2412$ | $0.4421$ | $0.228$ |
| **B2m** | MHC-I antigen presentation | 18.304 | $+0.0968$ | 0.1270 | $0.4708$ | $0.6074$ | $0.091$ |
| **Nfkb1** | NF-$\kappa$B transcription factor | 0.340 | $-0.0034$ | 0.0050 | $0.5116$ | $0.6074$ | $0.081$ |
| **Atm** | Double-strand break kinase | 0.157 | $-0.0050$ | 0.0077 | $0.5371$ | $0.6074$ | $0.148$ |
| **Xrcc5** | Non-homologous end joining (Ku80) | 0.090 | $+0.0012$ | 0.0020 | $0.5699$ | $0.6074$ | $0.144$ |
| **Cdkn2a** | Ink4a/Arf senescence locus | 0.0002 | $-0.00002$ | 0.00003 | $0.6074$ | $0.6074$ | $0.111$ |

*Cdkn1a* (p21) upregulation represents the strongest, statistically validated senescence signal in mouse lung fibroblasts ($p = 0.0029, q = 0.0322, R^2 = 0.777$). *Sirt6* exhibits an exploratory downward trend ($p = 0.0245$), but does not survive FDR multiplicity correction across the 11 candidate markers ($q = 0.1350$).

---

### 2.5 Compositional Remodeling Across Mammalian Lifespans

Using Aitchison Centered Log-Ratio (CLR) transformation on donor-level cell counts ([`scripts/run_h2_composition.py`](file:///home/reza/projects/agescape/scripts/run_h2_composition.py)), we evaluated cellular remodeling across mouse and human lungs ([`results/h2-composition-mouse-lung.csv`](file:///home/reza/projects/agescape/results/h2-composition-mouse-lung.csv), [`results/h2-composition-human-lung.csv`](file:///home/reza/projects/agescape/results/h2-composition-human-lung.csv)):

```mermaid
flowchart TD
    subgraph Mouse Lung Remodeling (N=23 mice)
        M1[T cells expand: q = 0.0003]
        M2[Vein Endothelial expand: q = 0.0173]
        M3[Fibroblasts decline: q = 0.0465]
    end
    subgraph Human Adult Lung Remodeling (N=117 adult donors)
        H1[Club cells collapse: q = 0.0123]
        H2[Vein Endothelial decline: q = 0.0220]
        H3[AT2 cells expand: q = 0.0220]
    end
    M1 -.->|Concordance Unproven: p = 1.000| H1
```

* **Mouse Lung (N = 23 mice):** Six cell types exhibited significant remodeling at FDR $q < 0.05$: T cells ($\beta = +0.0768/\text{mo}, q = 0.0003$), vein endothelial cells ($\beta = +0.0304/\text{mo}, q = 0.0173$), pulmonary interstitial fibroblasts ($q = 0.0038$), lung fibroblasts ($\beta = -0.0450/\text{mo}, q = 0.0465$), CD8+ T cells ($q = 0.0465$), and NK cells ($q = 0.0465$).
* **Human Adult Lung (N = 117 adult donors $\ge 20$ yr):** Four cell types exhibited significant remodeling at FDR $q < 0.05$: club cells ($\beta = -0.0678/\text{yr}, q = 0.0123$), vein endothelial cells ($\beta = -0.0326/\text{yr}, q = 0.0220$), pulmonary alveolar type 2 cells ($\beta = +0.0442/\text{yr}, q = 0.0220$), and mast cells ($\beta = +0.0333/\text{yr}, q = 0.0220$).
* **Resolution of the B-Cell Confound:** In the full human cohort including pediatric donors (ages 2–81, $N = 137$), B cells appeared to collapse precipitously ($\beta = -0.0411/\text{yr}, p = 3.8 \times 10^{-6}$). Restricting to adults ($\ge 20$ yr, $N = 117$) flattened the slope completely ($\beta = -0.0027/\text{yr}, p = 0.7601$), demonstrating that the apparent B-cell loss reflected childhood immune maturation rather than gerontological aging.
* **Cross-Species Trajectory Comparison:** Among seven homologous cell-type pairs, sign concordance was 3/7 (42.9%), yielding a two-sided exact binomial $p = 1.000$ and Spearman $\rho = -0.286$ ($p = 0.535$). Thus, while compositional remodeling is pronounced within each species, the specific trajectories do not adhere to a shared mammalian schedule.

---

### 2.6 The Universal Mammalian Epigenetic Clock

To contrast volatile transcriptional states with epigenomic aging, we evaluated Universal Clock 2 (Lu et al., *Nature Aging* 2023) across 1,000 array profiles ([`scripts/run_h1_epigenetic_drift.py`](file:///home/reza/projects/agescape/scripts/run_h1_epigenetic_drift.py)):

$$\text{RelativeAge} = \frac{\text{Age} + \text{GestationT}}{\text{MaxLifespan} + \text{GestationT}}$$
$$\hat{\eta} = \beta_0 + \sum_{j=1}^{816} \beta_j \cdot \text{CpG}_{j}$$
$$\text{RelativeAge}_{\text{pred}} = \exp(-\exp(-\hat{\eta}))$$
$$\text{Age}_{\text{pred}} = \text{RelativeAge}_{\text{pred}} \times (\text{MaxLifespan} + \text{GestationT}) - \text{GestationT}$$

| Species | Tissue / Cohort | N Samples | Probes Matched | Pearson r | $R^2$ | MAE (Years) |
|---|---|---|---|---|---|---|
| **Mouse** | Liver Lifespan (GSE199979, BXD panel) | 339 | 724 / 816 | **$0.923$** | **$0.851$** | **$0.36\text{ yr}$ ($4.3\text{ mo}$)** |
| **Human** | All Tissues (GSE184211, 8 organs) | 661 | 816 / 816 | **$0.801$** | **$0.642$** | **$5.06\text{ yr}$** |
| Human | Kidney | 97 | 816 / 816 | $0.876$ | $0.767$ | $4.10\text{ yr}$ |
| Human | Lung | 113 | 816 / 816 | $0.839$ | $0.704$ | $4.89\text{ yr}$ |
| Human | Spleen | 62 | 816 / 816 | $0.838$ | $0.701$ | $5.27\text{ yr}$ |
| Human | Adipose | 57 | 816 / 816 | $0.823$ | $0.677$ | $4.50\text{ yr}$ |
| Human | Liver | 97 | 816 / 816 | $0.817$ | $0.668$ | $4.96\text{ yr}$ |
| Human | Heart | 97 | 816 / 816 | $0.808$ | $0.652$ | $5.18\text{ yr}$ |
| Human | Blood | 32 | 816 / 816 | $0.799$ | $0.639$ | $6.14\text{ yr}$ |
| Human | Muscle | 57 | 816 / 816 | $0.728$ | $0.530$ | $5.73\text{ yr}$ |

Applying the proper double-exponential inverse link function resolved the severe prediction bias previously noted in linear approximations: mean absolute error dropped from $1.02\text{ yr}$ to **$0.36\text{ yr}$** in mice, and from $25.3\text{ yr}$ to **$5.06\text{ yr}$** in humans.

---

## 3. Discussion

Our empirical findings across mouse, human, and naked mole-rat datasets reconcile several longstanding controversies in mammalian geroscience:

1. **Absence of a Universal Whole-Tissue Script:** Evolutionary senescence theory posits that natural selection does not shape an active post-reproductive genetic program. In accordance with this principle, whole-tissue transcriptional clocks fail to transport across species when adjusted for demographic and postmortem confounders.
2. **Cell-Type-Private Dynamics:** Mammalian aging operates in distinct cellular compartments. The cross-species transport of the mouse signature to human OPCs (4/4 cohorts) and microglia (3/4 cohorts + M10 adult lifespan) indicates that glial populations undergoing lifelong metabolic and immune challenges share conserved transcriptional features that are obscured in bulk tissue homogenates.
3. **The Naked Mole-Rat Negative Control:** The complete stability of the mouse aging signature in naked mole-rat skin ($\beta = +0.0035, p = 0.925$; TOST equivalence $p = 0.027$) proves that the signature tracks biological senescence rather than chronometric time.
4. **Physical Senescence vs. Epigenetic Timing:** In long-lived structural cells, cellular aging is marked by *Cdkn1a* upregulation ($q = 0.0322$) and subtle transcriptional noise increases. In contrast, pan-mammalian DNA methylation acts as an exceptionally well-calibrated odometer across species, tracking elapsed life history standardized to maximum species lifespan.

### Limitations

1. **Cross-Sectional Observational Data:** All cohorts analyzed are cross-sectional; longitudinal trajectories in individual animals cannot be directly measured from single-cell atlases.
2. **Mouse Vivarium vs. Human Outbred Diversity:** Laboratory mice (*Tabula Muris Senis*) were housed under specific pathogen-free conditions with uniform diets. Human postmortem cohorts reflect genetically outbred individuals with diverse lifetime exposures, varying postmortem intervals, and agonal states.
3. **Mouse Donor Sample Size:** The number of mouse donors with $\ge 20$ lung fibroblasts is limited ($N = 10$), which restricts statistical degrees of freedom. While *Cdkn1a* achieves $q = 0.0322$, *Sirt6* ($q = 0.1350$) remains an exploratory hypothesis requiring targeted perturbational testing.
4. **Epigenetic Cohort Annotations:** GSE199979 is derived from female BXD recombinant inbred strains with high-fat diet intervention arms, and GSE184211 provides 661 tissue samples without unique decedent identifiers.

---

## 4. Methods

### 4.1 Dataset Provenance and Verification
* **Tabula Muris Senis Droplet:** 245,389 raw droplets across 23 mice ([`data/raw/tms_droplet_raw.h5ad`](file:///home/reza/projects/agescape/data/raw/tms_droplet_raw.h5ad)).
* **Human Cohorts:** 10 cohorts spanning 621 donors (SEA-AD DLPFC, MSSM, HBCC, RADC, CELLxGENE Lung_Integrated).
* **Naked Mole-Rat:** Single-cell and pseudobulk skin RNA-seq from Hilton et al. 2020 ([`GSE152886`](file:///home/reza/projects/agescape/data/raw/GSE152886_nmr_skin_raw.csv.gz)).
* **Pan-Mammalian Methylation:** GSE184211 (661 human tissues) and GSE199979 (339 mouse liver samples) profiled on GPL28271. Cryptographic SHA-256 manifests are maintained in [`data/manifest/`](file:///home/reza/projects/agescape/data/manifest/).

### 4.2 Single-Cell Preprocessing and QC
Raw droplet AnnData objects were processed using `scanpy` (1.10.4):
* Minimum detected genes: 500 genes/cell.
* Maximum detected genes: 5,000 genes/cell (multiplet exclusion).
* Minimum UMI counts: 1,000 counts/cell.
* Doublet detection: `sc.pp.scrublet` per tissue (expected doublet rate = 0.06).
* Cleaned object saved at [`data/interim/m1-qc-cleaned.h5ad`](file:///home/reza/projects/agescape/data/interim/m1-qc-cleaned.h5ad).

### 4.3 UMI Rarefaction and Distance-to-Median (DM)
To prevent library-size confounding:
1. Singlet cells with $\ge 1,500$ UMIs were downsampled to exactly 1,500 UMIs via multinomial drawing:
   $$\mathbf{x}_{\text{rarefied}} \sim \text{Multinomial}\left(1500, \mathbf{p} = \frac{\mathbf{x}}{\sum x_j}\right)$$
2. For each cell type, expressed genes were identified ($\text{detection rate} \ge 5\%$, mean UMI $\ge 0.01$).
3. A pooled mean-variance trend was established using a 25-bin running median in $\log_{10}(\mu)$ vs $\log_{10}(CV^2)$ space.
4. For donor $i$ with $\ge 20$ cells, $\text{DM}_{ig} = \log_{10}(CV_{ig}^2) - \text{baseline}(\log_{10}(\mu_{ig}))$.
5. Multivariate dispersion $D_i$ was computed as mean Euclidean distance to donor centroid in $\log_2(\text{rarefied\_counts} + 1)$ space.
6. OLS regression was fit: $\text{Metric}_i \sim \text{Age}_i + \text{Sex}_i$.

### 4.4 Compositional Remodeling via Aitchison CLR
Cell-type counts per donor $n_{ik}$ were transformed via Centered Log-Ratio:
$$p_{ik} = \frac{n_{ik} + 0.5}{\sum_j (n_{ij} + 0.5)}$$
$$\text{CLR}_{ik} = \ln p_{ik} - \frac{1}{K}\sum_{j=1}^K \ln p_{ij}$$
Simplex closure invariant $\sum_k \text{CLR}_{ik} = 0$ was verified for every donor. OLS regressions were fit against Age + Sex. Cross-species concordance was evaluated with two-sided exact binomial tests and Spearman rank correlation.

### 4.5 Universal Epigenetic Clock 2
Clock 2 coefficients (816 CpGs) were loaded from Lu et al. (2023). Chronological age was scaled with gestation time ($GT_{\text{human}} = 0.76\text{ yr}, GT_{\text{mouse}} = 0.05\text{ yr}$) and maximum lifespan ($ML_{\text{human}} = 122.5\text{ yr}, ML_{\text{mouse}} = 4.0\text{ yr}$). Double-exponential inverse link function $\exp(-\exp(-\hat{\eta}))$ was applied to the linear predictor to obtain calibrated age predictions.

### 4.6 Marker Pseudobulk Regression and Multiplicity Control
Donor pseudobulk expression was calculated in mouse lung fibroblasts for 11 targeted markers (*Cdkn1a, Sirt6, Ccl2, Stat1, Cxcl12, Parp1, B2m, Nfkb1, Atm, Xrcc5, Cdkn2a*). OLS regression was fit: $Y_{ig} \sim \text{Age}_i + \text{Sex}_i$. Benjamini-Hochberg FDR was computed across all 11 markers.

---

## Data and Code Availability

All code, pipelines, and tests are public at `https://github.com/Abolfazl-Alipour/agescape`.  
An interactive data explorer is available at `https://rezaalipour.com/agescape/`.  
All 16 automated validation tests run cleanly via `pytest tests/ -v`.
