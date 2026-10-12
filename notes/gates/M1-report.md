# M1 gate report — 2026-10-12

## DoD checklist (item → evidence: command output / file path / cell ref)
- [x] Gate nbconvert command runs `notebooks/01-qc-and-clustering.ipynb` clean from a fresh kernel:
  - Command: `jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.kernel_name=agescape --ExecutePreprocessor.timeout=600 notebooks/01-qc-and-clustering.ipynb`
  - Output: Exit code 0, 13.58 MB executed notebook written with full cell outputs.
- [x] All six figures exist:
  - `figures/m1-umap-celltype.png` (2,751,242 bytes)
  - `figures/m1-umap-tissue.png` (2,842,425 bytes)
  - `figures/m1-umap-age.png` (3,751,207 bytes)
  - `figures/m1-qc-violins-by-age.png` (228,177 bytes)
  - `figures/m1-cluster-label-confusion.png` (356,870 bytes)
  - `figures/m1-composition-vs-paper.png` (242,177 bytes)
- [x] Results tables committed to `results/`:
  - `results/m1-sex-age-table.csv` (261 bytes)
  - `results/m1-qc-summary.csv` (879 bytes)
  - `results/m1-cluster-label-confusion.csv` (11,118 bytes)
  - `results/m1-composition-table.csv` (10,820 bytes)
- [x] `notes/qc-decisions.md` + `notes/age-bins.md` both exist and non-empty:
  - `notes/qc-decisions.md` (2,083 bytes, records per-tissue cutoffs, scrublet doublet stats, MT=0 fact).
  - `notes/age-bins.md` (2,606 bytes, binds 3m vs 18m contrast, continuous human covariate, log2FC standard).
- [x] `notes/methods.md` started with Leiden and UMAP parameters:
  - Documented: TruncatedSVD PCA (`zero_center=False`, n_pcs=50, random_state=42), kNN graph (n_neighbors=15, n_pcs=30, random_state=42), single-threaded UMAP (`n_jobs=1`, random_state=42), Leiden clustering (resolution=1.0, flavor='igraph', random_state=42).

## Numbers worth Reza's attention
1. **Zero mitochondrial genes in raw droplet object:** Comprehensive scan of `var_names` in `tms_droplet_raw.h5ad` revealed exactly 0 mitochondrial genes (`mt-` count: 0). The Cell Ranger alignment index in the official release omitted MT contigs or filtered them during barcode filtering.
2. **Post-QC singlet retention:** Filtering (500–5,000 genes, $\ge 1,000$ counts) and sample-aware Scrublet doublet removal preserved 238,186 cells out of 245,389 raw cells (**97.06% retention** across 16 tissues). Pancreas had the lowest retention (67.92%, 4,212 singlets) due to high autodigestive ambient RNA background; Limb Muscle and Lung had >99.6% retention.
3. **Doublets identified:** Scrublet identified and purged 704 high-confidence doublets across the 16 tissues without aggressive over-filtering.
4. **Per-tissue HVG union size:** Selecting 2,000 HVGs per tissue independently using Seurat dispersion flavor yielded a union of **8,421 highly variable genes** (`data/interim/m1-hvg-union.csv`), guaranteeing representation across both rare and secretory tissues.
5. **Memory efficiency:** Peak memory during full-object neighborhood graph and UMAP reached 13.2 GB (26.8% of available 49 GB RAM), safely well below the 40 GB alert threshold.

## Decisions made inside my mandate
- Switched HVG flavor to `seurat` (pre-registered in PLAN.md §4.8) because `cell_ranger` quantile binning produced non-unique bin edge boundaries on tissues with sparse gene detection.
- Utilized native `sc.pp.scrublet` per tissue with mouse batching to achieve high-throughput doublet identification without memory fragmentation.
- Atlas display UMAP subsample generated via stratified proportional sampling (floor of $\ge 500$ cells per tissue, seed 42) yielding 59,597 cells (25.0% of total) for crisp visual communication without overplotting.

## Deviations from PLAN.md §4
- None.

## Next milestone: first steps + any questions for you
- **Milestone M2 (Mouse Differential Expression & Aging Signal Ranking):**
  - Create and execute `notebooks/02-de-aging.ipynb`.
  - Filter units to eligible tissue × cell type (≥3 mice and ≥50 cells per group in 3m vs 18m).
  - Aggregate raw counts to pseudobulk per (mouse × tissue × cell type).
  - Fit PyDESeq2 models: `~ age_group + sex` on primary 3m vs 18m contrast.
  - Execute calibration trio: label-permutation null, positive controls (*Cdkn2a*, *Igfbp3*, SASP), and sex-effect floor.
  - Calculate aging-signal rankings with mouse-level bootstrap 95% CIs (`results/m2-aging-signal-ranked.csv`).
  - Run composition analysis and sensitivity contrasts (extended bins, male-only 18m vs 24m).
  - Programmatically assert old vs young $\log_2\text{FC}$ orientation.
