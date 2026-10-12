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
