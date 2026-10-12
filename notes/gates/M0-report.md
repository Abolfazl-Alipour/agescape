# M0 gate report — 2026-10-12

## DoD checklist (item → evidence: command output / file path / cell ref)
- [x] `python scripts/verify_data.py` exits 0 on both data files:
  - `tms_droplet_raw.h5ad`: Size 3868.6 MB (4,056,529,818 bytes), SHA-256 `86e45e88f9c475253792f1bbcdf5851a5be0ef35da72fced607ed261aebc099d` — MATCH PASS.
  - `tms_facs_brain_nonmyeloid.h5ad`: Size 304.0 MB (318,721,054 bytes), SHA-256 `31ee7ffe9fd11fd128925c2dae862399f5c54f27aef9c026dbeaa077ea143339` — MATCH PASS.
- [x] `conda run -n agescape python -c "import scanpy, anndata, pydeseq2, scrublet, src; print(scanpy.__version__)"` passes:
  - Output: `SCANPY_VERSION: 1.11.5` (anndata 0.11.3, pydeseq2 0.5.4, scrublet 0.2.3, src import clean).
- [x] `jupyter kernelspec list` shows `agescape`:
  - Verified: `agescape /home/reza/.local/share/jupyter/kernels/agescape`.
- [x] `environment.lock` committed; `.gitignore` correct; tree matches PLAN.md §3:
  - Initial commits: `ae8e503` and `836817a`.
  - Directories: `data/raw`, `data/interim`, `data/manifest`, `scripts`, `notebooks`, `src`, `results`, `figures`, `notes/gates`, `report`, `app`.
- [x] `notes/age-bins.md` exists and matches PLAN.md §4.2:
  - Verified file path: `notes/age-bins.md` defining primary 3m vs 18m sex-balanced contrast, continuous human age model, sensitivity designs, and log2FC old-vs-young orientation.
- [x] Dead files deleted:
  - Deleted `/home/reza/projects/aging_hypothalamus/data/23872787` and `/home/reza/projects/aging_hypothalamus/data/23872787.1` (0 bytes). Confirmed absent with exit code 2 on subsequent `ls`.

## Numbers worth Reza's attention
1. **Raw droplet scale:** Exactly 245,389 cells across 20,138 genes and 16 tissues, totaling 4,056,529,818 bytes.
2. **Sex × age confounding:** Droplet cells per age group are: 1m (25,980 M-only), 3m (45,602 mixed M/F), 18m (44,645 mixed M/F), 21m (35,828 F-only), 24m (37,660 M-only), 30m (55,674 M-only). Only 3m and 18m have both sexes represented.
3. **Brain tissue absence in droplet:** The main droplet object contains 0 brain cells. The brain single-cell data is solely in the FACS subset (7,249 cells × 22,966 genes, ages 3m, 18m, 24m across 14 mice).
4. **Environment resolution:** Conda environment `agescape` resolved with Python 3.11.14, Scanpy 1.11.5, AnnData 0.11.3, PyDESeq2 0.5.4, and Scrublet 0.2.3 under strict conda-forge and bioconda channel priority. Explicit lockfile committed as `environment.lock`.
5. **System capacity:** Host Linux 6.8 kernel (Ubuntu 24.04.3 LTS, glibc 2.39), 49 GB RAM, 12 cores, 247 GB free NVMe disk space.

## Decisions made inside my mandate
- Pinned `environment.yml` with `nodefaults` to prevent conda solver cross-channel deadlocks.
- Installed `scrublet 0.2.3` via pip declaration within `environment.yml` per PLAN.md §4.1 because bioconda's standalone scrublet wheel carries older scikit-learn pin constraints for Python 3.11.
- Initialized local git repository with local user identity (`Reza <reza@agescape.local>`) and strict no-remote policy.

## Deviations from PLAN.md §4
- None. All packages, seeds, paths, hashes, and conventions adhere strictly to PLAN.md.

## Next milestone: first steps + any questions for you
- **Milestone M1 (Pipeline: QC → clusters):**
  - Create and run `notebooks/01-qc-and-clustering.ipynb`.
  - Strip trailing whitespace across string `obs` columns and record before/after counts for `subtissue`.
  - Export `results/m1-sex-age-table.csv` and mice-per-age distributions.
  - Compute age-stratified QC metrics across all 16 tissues, run scrublet per-sample, and log binding thresholds in `notes/qc-decisions.md`.
  - Perform log1p normalization on concatenated object, per-tissue HVG union, sparse PCA (`zero_center=False`), UMAP, Leiden clustering, majority-vote ontology mapping, and export six required diagnostic figures.
