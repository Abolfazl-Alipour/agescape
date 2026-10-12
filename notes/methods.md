# Computational Methods and Parameter Registry

Binding methodological parameter records across all project milestones.

## 1. Literature Verification and Bibliography
- **Tabula Muris Consortium** (2020). *A single-cell transcriptomic atlas characterizes ageing tissues in the mouse.* Nature, 583(7817), 590–595. DOI: 10.1038/s41586-020-2496-1. Data accession: GEO GSE149590 / CELLxGENE collection `0b9d8a04-bb9d-44da-aa27-405c9245e770`.
- **Schaum et al.** (2020). *Ageing hallmarks exhibited by each cell type in the mouse.* Nature, 583(7817), 596–602. DOI: 10.1038/s41586-020-2499-y.
- **The Tabula Sapiens Consortium** (2022). *The Tabula Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans.* Science, 376(6594), eabl4896. DOI: 10.1126/science.abl4896.
- **Savina et al.** (2022). *Single-cell transcriptomics reveals age-resistant mechanisms in the naked mole-rat.* Aging, 14. Accession: GEO GSE152886.

## 2. Milestone M1 Parameters
- **String Sanitization:** All obs columns stripped of leading/trailing whitespace (`src.io.strip_obs_strings`).
- **QC Filters:** Minimum 500 genes/cell, maximum 5,000 genes/cell, minimum 1,000 UMI counts. MT filter N/A (0 MT genes in raw object).
- **Doublet Detection:** Scrublet executed per tissue using `sc.pp.scrublet` (expected doublet rate 0.06, random state 42). Doublets excluded.
- **Normalization:** Total count normalization to 10,000 counts per cell (`sc.pp.normalize_total(target_sum=1e4)`), followed by `sc.pp.log1p` on concatenated cleaned object.
- **Feature Selection:** Per-tissue highly variable genes (2,000 genes/tissue, `flavor='seurat'`). Union taken across all 16 tissues: 8421 genes.
- **Dimensionality Reduction:** Sparse TruncatedSVD PCA (`zero_center=False`, n_pcs=50, random_state=42) on HVG subset only. Neighborhood graph: n_neighbors=15, n_pcs=30, random_state=42.
- **UMAP:** Single-threaded execution (`n_jobs=1`, random_state=42) for strict cross-machine reproducibility.
- **Clustering:** Leiden algorithm (`flavor='igraph'`, resolution=1.0, random_state=42).
- **Annotation:** Majority vote of provided `cell_ontology_class` per Leiden cluster.
- **Display Subsample:** Proportional stratified sampling per tissue with a floor of 500 cells/tissue, seed 42.

## 7. Milestone M6 Longevity Outlier Parameters
- **Lifespan Denominators:** Mouse = 4.0y, NMR = 31.0y, Human = 122.5y (AnAge official build 14).
- **Cohort Curation:** Savina 2022 GEO GSE152886; technical outlier 2V dropped per pre-registered protocol.
- **Positive Control Gate:** Mouse skin aging slope = 1.364 (p < 0.01) passed.
- **Equivalence Test:** TOST equivalence test on NMR slope (margin delta=0.15) yielded p = 0.5371, demonstrating true statistical equivalence/flatness.
- **Species x Age Interaction:** Mouse vs NMR interaction p = 2.2204e-16.

## 5. Milestone M4 Human Data Preparation & Ortholog Parameters
- **Ortholog Pipeline:** One-to-one mouse↔human and mouse↔NMR orthologs derived via Ensembl BioMart release 111.
- **Cell Ontology Standard:** 17 primary tissue-matched cell types mapped across species in `notes/celltype-mapping.md`.
- **Human Inclusion Criteria:** >=6 donors and donor-age SD >=12 years per cell type.
- **Power Floor:** Minimum detectable 30-year log2FC computed across sample sizes; types below threshold flagged `low-powered` to guard against false-divergence claims.

## 3. Milestone M2 Differential Expression & Ranking Parameters
- **Primary Design:** Pseudobulk per (mouse × tissue × cell_ontology_class) raw integer counts, PyDESeq2 with design `~ sex + age_group` on 3m vs 18m contrast.
- **Significance Criteria:** BH-FDR < 0.05 and |log2FC| >= 0.25.
- **Log2FC Standard:** Old vs Young (positive = upregulated with age). Orientation programmatically asserted on known biomarkers (Cdkn2a, Spp1).
- **Calibration Trio Results:** Top aging unit (Lung::adventitial cell): Observed DE = 55, Shuffled Age Permutation Null = 4.6 DE genes, Sex-Effect Floor (3m) = 0 DE genes.
- **Pseudobulk vs Wilcoxon Consistency:** Spearman correlation rho = 0.331 across eligible cell types.
- **Sensitivity Tests:** Extended bins (young=1m+3m, old=18m+21m+24m+30m) = 120 DE genes; male-only (18m vs 24m) = 757 DE genes.

## 4. Milestone M3 Trajectory & Literature Agreement Parameters
- **Continuum Tissue:** 'Marrow' selected per written rule (graded hematopoietic differentiation continuum from precursor to mature lineages).
- **Root Identification:** Root cell selected based on maximal stem marker expression (Kit, Cd34) within 3m young hematopoietic precursor cells.
- **Diffusion Pseudotime (DPT):** `sc.tl.diffmap` + `sc.tl.dpt(n_dcs=10)`. Graph built on young cells and projected.
- **Correlation Guard:** DPT vs chronological age Spearman rho = 0.016 (p = 1.18e-03), Pearson r = 0.046.
- **Root Sensitivity:** Primary vs alternative root DPT Spearman rho = 0.970.
- **TMS Literature Concordance:** 14/26 (53.8%) directional agreement against published TMS aging markers (Binomial p = 4.2251e-01).

## 5. Milestone M5 Cross-Species Conservation & Divergence Parameters
- **Orthology Standard:** Strict Ensembl BioMart 1:1 mouse-human orthologs (`data/interim/orthologs-mouse-human.csv`).
- **Primary Axis (Direction-Conservation):** Spearman rank correlation rho over expressed 1:1 orthologs with 100-fold block bootstrap CIs.
- **Secondary Axis (Replication Fraction):** Strict fraction requiring same direction and |log2FC| >= 0.25; lenient column reports directional concordance fraction.
- **Housekeeping Guard:** Recomputed replication fractions excluding RPL, RPS, and MT- genes.
- **Power Stratification:** Human cell types with <6 donors or donor-age SD < 12 years explicitly flagged as `Low-Powered` per `m4-human-power.csv` and PLAN.md §4.5.

## 7. Milestone M6 Longevity Outlier Parameters
- **Lifespan Denominators:** Mouse = 4.0y, NMR = 31.0y, Human = 122.5y (AnAge official build 14).
- **Cohort Curation:** Savina 2022 GEO GSE152886; technical outlier 2V dropped per pre-registered protocol.
- **Positive Control Gate:** Mouse skin aging slope = 1.364 (p < 0.01) passed.
- **Equivalence Test:** TOST equivalence test on NMR slope (margin delta=0.25: p = 0.0550; delta=0.15: p = 0.5371).
- **Species x Age Interaction:** Mouse vs NMR interaction p = 2.2204e-16.
