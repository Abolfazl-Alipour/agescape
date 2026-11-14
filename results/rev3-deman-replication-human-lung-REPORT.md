# Item 5 Blocker Report: De Man et al. 2026 Human Lung Replication (`GSE281219` Baylor Cohort)

- **Date**: 2026-11-14
- **Status**: **BLOCKED** for the 32 Baylor donor lungs (`GSE281219`); **IPF Cell Atlas 28-control arm already completed** in [reports/followup-gse281219.md](file:///home/reza/projects/agescape-revision/reports/followup-gse281219.md) and [followup-de-man-entropy-models.csv](file:///home/reza/projects/agescape-revision/results/followup-de-man-entropy-models.csv).
- **Directive followed**: Per project rules ("or a blocker report if no labeled Baylor object exists" and "Never silently substitute a smaller analysis or a different question"), no custom de novo cell annotation was substituted in place of De Man et al.'s published cell-type assignments.

---

## 1. Exact Blocker for the 32 Baylor Lungs (`GSE281219`)

Live verification on **2026-11-14** across NCBI GEO, the CELLxGENE Curation API, De Man et al. 2026 (*Nature Communications*, `https://doi.org/10.1038/s41467-026-68810-9`) full-text XML, and the authors' public GitHub repository (`rubendeman/Lung-Aging-Cell-Atlas`, commit `0c11d44ce4103451327fadd106e1e6872edefa50`) establishes that **no publicly downloadable barcode-to-cell-type annotation table or labeled integrated object exists for the 32 Baylor lungs**:

1. **NCBI GEO `GSE281219` FTP directory (`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE281nnn/GSE281219/suppl/filelist.txt`)**: Lists only `GSE281219_RAW.tar` (`1,589,473,280` bytes; last modified `2025-12-31`), containing exactly `96` raw 10x CellRanger files (`32` samples `GSM8615073` to `GSM8615104` x `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`). No integrated `.RData`, `.rds`, `.h5ad`, or cell-level metadata table is present on GEO.
2. **Article Data Availability (`fullTextXML.xml`)**: Points exclusively to GEO `GSE281219` (raw Baylor scRNA-seq), `GSE136831` (IPF Cell Atlas), `GSE47460` (LTRC bulk), and `https://github.com/rubendeman/Lung-Aging-Cell-Atlas`. No Zenodo, Figshare, Dryad, or Synapse deposit is listed.
3. **CELLxGENE Curation API (`https://api.cellxgene.cziscience.com/curation/v1/collections`)**: Returns `0` matching collections for `GSE281219`, `s41467-026-68810-9`, `De Man`, or `Lung Aging Cell Atlas`.
4. **Authors' GitHub repository (`rubendeman/Lung-Aging-Cell-Atlas`, `60` tracked files)**:
   - [scRNAseq/add_meta.R](https://github.com/rubendeman/Lung-Aging-Cell-Atlas/blob/0c11d44ce4103451327fadd106e1e6872edefa50/scRNAseq/add_meta.R) reads a pre-built integrated Seurat object (`immune.combined` with `predicted.id`), attaches donor metadata from `demo2_17.csv`, and saves `agingseurat.RData`.
   - [scRNAseq/entropy.R](https://github.com/rubendeman/Lung-Aging-Cell-Atlas/blob/0c11d44ce4103451327fadd106e1e6872edefa50/scRNAseq/entropy.R) loads the saved integrated object (`immune.combined`) and subsets directly on `predicted.id %in% c('AT1', 'AT2', 'gCap', 'Aerocyte')`. Neither `agingseurat.RData` nor any barcode label export is in the repository.
5. **Local staging verified**: `/home/reza/projects/agescape-revision-data/de-man/GSE281219_RAW.tar` (`1,589,473,280` bytes, `sha256` `8c876b037c0ed46006f8fbedd79dedd778234bbc6bda39f0aad02c12aa8dd07d` recorded in [de_man_sources.json](file:///home/reza/projects/agescape-revision/data/manifest/de_man_sources.json)) is already downloaded on disk, along with `41467_2026_68810_MOESM3_ESM.csv` (`60` subjects: `32` Baylor + `28` IPF controls). Only the cell-type labels (`predicted.id` / `agingseurat.RData`) for Baylor cells are missing.

---

## 2. Summary of the Already-Completed 28 IPF Control Arm (`GSE136831`)

While the Baylor half cannot run without author-released cell labels, all **28 IPF Cell Atlas control subjects (`38` libraries, `96,303` control cells; `92,049` passing QC)** from `GSE136831` were previously analyzed at the donor level in [followup-gse281219.md](file:///home/reza/projects/agescape-revision/reports/followup-gse281219.md), [followup-de-man-entropy-models.csv](file:///home/reza/projects/agescape-revision/results/followup-de-man-entropy-models.csv), and [followup-de-man-entropy-robustness.csv](file:///home/reza/projects/agescape-revision/results/followup-de-man-entropy-robustness.csv):

- **Primary seed (`20261113`), donor-level OLS adjusted for sex and `log(actual_umi)`** (`BH` across the `4` target cell types per method): **zero** findings reach `q < 0.05`.
  - `published_KL_binomial` (`sqrt`-count `D_KL`): `AT1` (`n=19`, slope `-0.000382/yr, p=0.599, q=0.637`), `AT2` (`n=18`, `-0.000559/yr, p=0.521, q=0.637`), `Aerocyte` (`n=17`, `-0.000672/yr, p=0.637, q=0.637`), `gCap` (`n=22`, `-0.002005/yr, p=0.0687, q=0.275`; `6/20` random seeds at nominal `p < 0.05`, median `p = 0.0704`).
  - `shannon_1500` (`1,500`-UMI multinomial rarefaction): `AT1` (`n=24`, `+0.000629/yr, p=0.759, q=0.759`), `AT2` (`n=26`, `-0.000869/yr, p=0.599, q=0.759`), `Aerocyte` (`n=17`, `+0.001752/yr, p=0.641, q=0.759`), `gCap` (`n=22`, `+0.004459/yr, p=0.206, q=0.759`).
  - `shannon_full`: `AT1` (`+0.004575/yr, p=0.215, q=0.308`), `AT2` (`-0.001983/yr, p=0.387, q=0.387`), `Aerocyte` (`+0.005159/yr, p=0.231, q=0.308`), `gCap` (`+0.005619/yr, p=0.153, q=0.308`).

---

## 3. What Is Needed From Reza to Unblock the 32 Baylor Lungs

To run the Baylor arm without inventing cell labels, we need either:
1. Author-released barcode annotations (`agingseurat.RData`, `.rds`, `.h5ad`, or per-barcode CSV mapping `cell_barcode` + `orig.ident` to `predicted.id`), **or**
2. Explicit permission from Reza to run an independent reference-mapping step (for example Azimuth / scArches against HLCA, matching De Man et al.'s `predicted.id` pipeline) labeled clearly as a re-annotated sensitivity analysis rather than a replica of published cell labels.
