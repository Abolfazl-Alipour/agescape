# Extended dataset inventory and analysis decisions

Last verified: 2026-11-02.

## Datasets used in the current analysis

### Aging Fly Cell Atlas

- Primary source: [GEO GSE218661](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218661).
- Released objects: head and body H5AD files with cell annotations, sample identifiers, ages, sex, and total counts.
- Ages: 30, 50, and 70 days from AFCA, plus 5-day cells from the earlier Fly Cell Atlas.
- Primary design: 36 AFCA head libraries and 36 AFCA body libraries, representing 36 age-by-sex pools measured in both organs. The primary age contrast is 30 to 70 days. Thirty days is the young-adult expression reference.
- The 5-day cells are not part of the primary test because age and study of origin are perfectly confounded.
- The released matrix is log1p normalized to 10,000 per cell. Applying `expm1` recovers each cell's expression proportions up to a common scale. Hellinger distance uses those proportions. Raw-count Euclidean distance is not reported for this dataset.
- Composition means composition among annotated nuclei that survived the authors' processing, not absolute tissue cell counts.

### Gao and colleagues whole-worm atlas

- Primary source: [GEO GSE229022](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229022).
- Wild-type N2 worms fed OP50 are analyzed at days 1, 6, 12, and 14. Other genotypes, food conditions, and early-day variants are excluded from the normal-aging arm.
- Each library contains nuclei from a pool of approximately 2,000 worms. A pool is the inferential unit. Cells within a pool are not independent replicates.
- Raw counts and author cell annotations come from the integrated Seurat object released in the same GEO series.

### Roux and colleagues whole-worm atlas

- Primary sources: [GEO GSE208154](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208154), the [paper's data statement](https://www.sciencedirect.com/science/article/pii/S2211124723009130), and the [annotated H5AD release](https://storage.googleapis.com/worm_public/ad_worm_aging.h5ad).
- Ages are days 1, 3, 5, 8, 11, and 15. GEO labels the libraries as technical replicates. The day-15 object contains 219 cells after the authors' processing.
- This arm provides descriptive trajectories and technical-replicate consistency. It does not provide biological-replicate significance tests or an independent Bayesian estimate.

### Zebrafish testis atlas

- Primary sources: [GEO GSE275361](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE275361), the [paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11607696/), and the [authors' analysis repository](https://github.com/asposato/aging_zebrafish_testis).
- Ten individual male fish are analyzed at 5, 12, 20, 22, and 27 months. Technical channels from the same fish are aggregated before testing.
- The all-age analysis uses three comparable compartments: germ cells, immune cells, and blood plus other somatic cells. The finer germ-cell-stage analysis is restricted to the 5-to-22-month object, where comparable fine labels exist. Six fish contribute to that fine-grained test.
- Broad germ-cell expression displacement can reflect altered germ-cell-stage composition. It is not automatically cell-intrinsic drift.
- Incomplete spermatogenesis at 27 months is evaluated as a separate covariate. Exact age permutations are used for the six-fish fine-grained analysis.

## Verified datasets that remain limited

### Rat multi-tissue aging and caloric restriction

- Primary source: [GEO GSE137869](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137869).
- The accession and raw-count matrices are verified. One library is released per tissue, sex, and condition. These libraries do not provide biological replication for animal-level significance tests.
- Any reanalysis with newly generated clusters must be labeled descriptive. It must not reproduce the published p-values by treating cells as animals.

### Cynomolgus cardiopulmonary atlas

- The [paper's data availability record](https://www.nature.com/articles/s41422-020-00412-6) gives `CRA002577` for lung and `CRA002689` for heart.
- These accessions resolve the identifier question, but the count matrices and per-animal metadata have not yet been ingested. No result is claimed.

### Rhesus macaque brain lifespan preprint

- The [preprint](https://www.biorxiv.org/content/10.1101/2025.10.31.685880v2.full-text) reports 55 animals and 5.3 million nuclei. Its data statement says that data will be publicly available. A usable accession has not been verified.
- This dataset remains pending and is not treated as available.

## Datasets not suitable for H2 or H3

Bowhead whale, tortoise, Greenland shark, long-lived bats, ocean quahog, rockfish, and sea urchin do not currently provide the required single-cell age series. Genome-level, bulk-tissue, and single-animal data cannot estimate within-species cell-type fractions or per-cell expression displacement with age.

Hydra has a single-cell atlas but no young-versus-old age contrast. Adult-only naked mole-rat immune atlases likewise cannot test aging. These resources can support species maps or hypotheses, but they are not evidence for H2 or H3.

## Common analysis rules

1. The fish, mouse, worm pool, fly pool, or human donor is the inferential unit, never an individual cell.
2. H2 uses centered log-ratio transformed cell counts. A significant result means a change in relative abundance among recovered cells.
3. H3 uses Hellinger distance between each cell's complete expression-proportion vector and a cross-fitted young-adult reference.
4. UMI depth is included in sensitivity models. Raw-count Euclidean distances are not compared across assays.
5. Missing-data, total-cell-number, phenotype-adjusted, and exact-permutation checks are reported where applicable.
6. Bayesian intervals are exploratory when sampling diagnostics, small sample sizes, or model sensitivity limit interpretation. Frequentist sample-level tests remain the reference.
