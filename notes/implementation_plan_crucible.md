# Conserved Aging Explorer — Implementation Plan: 10-Round Scientific Review Crucible

**Date:** 2026-10-12  
**Target:** Milestones M1 through M8 Execution Plan  
**Crucible Protocol:** 10 rounds of brutally honest scientific review across 5 distinct perspectives (2 rounds each), followed by systematic plan revision and ground-truth verification.

---

## Part 1: Initial Implementation Plan Outline (M1–M8)

- **M1 (Pipeline: QC → Clusters):**
  Execute `notebooks/01-qc-and-clustering.ipynb`. Clean string obs (strip trailing whitespace). Generate `results/m1-sex-age-table.csv`. Apply age-stratified QC thresholds per tissue (16 tissues). Doublet removal via `scrublet` per sample. Concat cleaned tissues, log1p normalize, compute per-tissue HVGs (flavor `cell_ranger` or `seurat`), union → `data/interim/m1-hvg-union.csv`. Sparse PCA (`zero_center=False`), neighbors, single-threaded UMAP (`random_state=42`, `n_jobs=1`), Leiden clustering. Majority-vote annotation against `cell_ontology_class`. Generate 6 diagnostic figures.
- **M2 (Mouse Differential Expression & Aging Signal Ranking):**
  Execute `notebooks/02-de-aging.ipynb`. Filter units to eligible tissue × cell type (≥3 mice and ≥50 cells per group for 3m vs 18m). Pseudobulk aggregation per (mouse × tissue × cell type). Run PyDESeq2 (`~ age_group + sex`). Calibration trio: (1) permutation null, (2) positive controls (Cdkn2a, Igfbp3, SASP), (3) sex-effect floor within 3m. Rank aging signal = n DE genes (FDR < 0.05, |log2FC| ≥ 0.25) with bootstrap 95% CIs. Composition analysis (cell type fraction vs age, mouse unit). Sensitivity contrasts (extended bins, male-only 18m vs 24m). Orientation assertion cell on known-direction genes.
- **M3 (Trajectory & Literature Agreement):**
  Execute `notebooks/03-trajectory-agreement.ipynb`. Choose trajectory tissue via formal rule (strongest differentiation continuum in M1). Compute diffusion pseudotime (`sc.tl.dpt`) with guards: corr(pseudotime, age), build graph on young cells and project old, connectivity checks. Compare findings against Tabula Muris Senis (Nature 2020) published gene lists via continuous-age bridge model (`~ age_months + sex`). Construct agreement table (agree/disagree/cannot-test) with explicit causal explanations for discrepancies.
- **M4 (Human Data Preparation & Ortholog Mapping):**
  Execute `notebooks/04-human-prep.ipynb`. Query SenNet public token-free API for human scRNA-seq datasets meeting pre-registered rule (10x platform, max age span with ≥6 donors ≥60y). Fetch Tabula Sapiens v1 10x subset (backed mode, no whole-object in-memory ops >250k cells). Ensembl BioMart 1:1 mouse↔human ortholog mapping (and mouse↔naked mole-rat mapping). Build `notes/celltype-mapping.md` using Cell Ontology IDs. Calculate per-cell-type statistical power table (`results/m4-human-power.csv`).
- **M5 (Cross-Species Conservation & Divergence):**
  Execute `notebooks/05-conservation.ipynb`. Evaluate conservation on two independent axes: (1) effect-size Spearman correlation with block bootstrap CIs, (2) strict replication fraction at identical FDR thresholds. Null calibration via donor age permutation. Power stratification (flag low-powered cell types into separate partition). Compare against published conserved programs (Tyshkovskiy/Gladyshev, Karakaslar AP-1). Output headline `results/m5-divergence-table.csv`.
- **M6 (Longevity Outliers: Naked Mole-Rat):**
  Execute `notebooks/06-outliers.ipynb`. Fetch Savina et al. (2022) GEO GSE152886 processed skin counts (dropping sample 2V per paper). Map mouse→NMR orthologs. Score conserved program modules per cell (`sc.tl.score_genes`). Validate mouse skin positive control. Regress module scores on life-stage-normalized age (mouse vs human vs NMR). Test NMR slope flatness via Two One-Sided Tests (TOST) against pre-registered equivalence margin $\delta$. Reproduce Savina skin DE curve. Fetch AnAge official database and build cross-species longevity metadata table (`results/m6-species-table.csv`). Document caveats in `notes/outlier-limits.md`.
- **M7 (Comprehensive Scientific Report & Polish):**
  Execute all notebooks from clean kernels. Diff results against golden outputs. Write comprehensive 8–12 page scientific report `report/conserved-aging-programs.md` including pre-registration deviations log, failed-robustness log, and formal limitations. Conduct fresh-clone reproduction test (`scripts/run_all.py`).
- **M8 (Conserved Aging Explorer Web Application):**
  Execute `notebooks/07-app-assets.ipynb` to generate Parquet/Arrow assets (float16, gene panels, budget ≤300 MB total, ≤25 MB per view). Build modern, responsive, high-aesthetic static web app (`app/`) featuring signature longevity outlier contrast, gene trajectory search, conserved vs divergent program browser, categorical age UMAP viewer, and caveats explorer. Run `scripts/check_app_assets.py`.

---

## Part 2: The 10-Round Brutally Honest Scientific Crucible

### Perspective A: Computational Genomicist & Biostatistician
#### Round 1: Algorithmic Rigor & Matrix Memory Traps
- **Critique:** The plan mentions log1p normalization on the concatenated object of ~245,000 cells. While sparse matrix representation preserves memory, naive scaling (`sc.pp.scale`) outside HVGs or computing dense covariance matrices during PCA will instantly trigger a memory blowup exceeding the 49 GB RAM ceiling. Furthermore, how are size factors handled across 16 disparate tissues with differing cell sizes?
- **Finding:** Normalization must be library-size log1p, HVG selection must be performed **within each tissue independently** before taking the union to prevent high-expression secretory tissues (e.g., pancreas, liver) from dominating the top 2,000–3,000 HVGs. PCA must enforce `zero_center=False` and operate strictly on sparse data.

#### Round 2: PyDESeq2 Model Stability & Permutation Null Calibration
- **Critique:** PyDESeq2 fits negative binomial dispersion parameters iteratively per gene. In small sample regimes (e.g. 3 mice per group), dispersion estimates can fail to converge or collapse to extreme boundaries if count matrices contain high dropout zeros. Moreover, does the label permutation shuffle mice or cells?
- **Finding:** In pseudobulk DE, the replicate unit is the **mouse**, not the cell. Shuffling cells across age labels while keeping mouse aggregates creates false degrees of freedom. Shuffling must be performed at the **mouse demographic level** within eligible tissue × cell type units. In PyDESeq2, minimum dispersion floors and independent filtering must be strictly enforced.

---

### Perspective B: Comparative Evolutionary Biologist & Aging Specialist
#### Round 3: Orthology Resolution & Gene Family Ambiguity
- **Critique:** The plan relies on Ensembl BioMart 1:1 orthology. In mammalian aging, key gene families (e.g., heat shock proteins *Hspa*, major histocompatibility complex *H2/HLA*, cytochrome P450 *Cyp*, ribosomal subunits) frequently have many-to-many or one-to-many duplications. A naive 1:1 filter will silently discard biologically critical stress-response genes or misattribute ortholog pairs if Ensembl release versions drift.
- **Finding:** Must pre-register the exact Ensembl release (e.g. Ensembl release 112 / GRCm39 to GRCh38). For many-to-many gene families, they must be explicitly partitioned into an "excluded due to gene duplication" log rather than silently dropping them. Crucially, housekeeping gene sets (*Rpl*, *Rps*, mitochondrial) must be evaluated separately to ensure that technical ribosomal/mitochondrial expression shifts do not masquerade as conserved aging programs.

#### Round 4: Cross-Species Life-Stage Normalization & Survivorship Confounding
- **Critique:** In M6, module scores are regressed on "life-stage-normalized age". How is this defined across species? If mouse maximum lifespan is taken as 3.5 years (42 months) and human as 122 years, chronological age converted to fraction of maximum lifespan assumes linear aging kinetics across the entire life course. Furthermore, the 30-month mouse cohort in TMS is a survivor cohort (~10–20% survival in C57BL/6J), whereas 60–80 year old humans in SenNet/Tabula Sapiens represent standard post-reproductive aging, and the naked mole-rats in Savina 2022 are young/middle-aged (~3–4 years out of a >30 year lifespan).
- **Finding:** The plan must explicitly calculate normalized age using pre-registered maximum lifespan benchmarks from AnAge (Mouse: 4.0 yr / 48 mo; Human: 122.5 yr; Naked mole-rat: 37.0 yr). The limitations of comparing chronological fractions across species with vastly divergent Gompertz mortality curves must be documented as an upfront scientific caveat. In NMR, the TOST margin $\delta$ must be derived from the biological noise floor of the mouse skin young-old contrast, not chosen post hoc.

---

### Perspective C: Systems Architect & Data Platform Engineer
#### Round 5: Out-of-Core Processing & Disk Bandwidth Guards
- **Critique:** Tabula Sapiens whole-body dataset exceeds 250,000 cells. Downloading and reading whole-body human AnnData into memory alongside mouse objects will exhaust RAM or induce heavy disk swapping.
- **Finding:** M4 must mandate out-of-core HDF5 backed mode (`sc.read_h5ad(..., backed='r')`) when querying Tabula Sapiens. Only the overlapping tissue slices (e.g., Skin, Spleen, Liver, Lung, etc.) may be loaded into active RAM. Subsampling for interactive visualization must enforce strict ceiling caps ($\le 150,000$ cells total) and export coordinates in Float16 Parquet format.

#### Round 6: Web Application Client-Side Query Budget
- **Critique:** If the web app serves whole AnnData objects or large JSON matrices, client browsers will suffer massive latencies, out-of-memory crashes on mobile, and failed downloads.
- **Finding:** The app must be an entirely static client-side architecture using precomputed Parquet/Arrow tables and curated gene signatures. Expression queries must load pre-aggregated tissue × cell-type median profiles or sparse gene-specific vector chunks ($\le 25$ MB per view, $\le 300$ MB total application payload).

---

### Perspective D: Reproducibility & Open Science Auditor
#### Round 7: Golden Output Verification & Strict Non-Fabrication
- **Critique:** If intermediate calculations produce slightly varying results due to multithreaded UMAP or randomized SVD in PCA, how is exact reproducibility guaranteed across machines?
- **Finding:** Determinism policy must be ruthlessly enforced: `src.set_seeds(42)` at the top of every script/notebook; `random_state=42` on PCA, UMAP, and Leiden; UMAP `n_jobs=1` to guarantee single-threaded floating-point determinism; Golden outputs (`results/*.csv`) committed with exact SHA-256 digests. `scripts/run_all.py` must perform automated numeric tolerance checks against committed golden tables.

#### Round 8: Zero-Mock Data Policy & Complete Data Provenance
- **Critique:** All data must be authentic public records. Any mock generation, hardcoded arrays of fake expression, or synthetic values violates the core project charter.
- **Finding:** Every single cell, UMI count, gene identifier, p-value, and fold change must trace directly to executed raw data downloaded from authentic public repositories (CELLxGENE, GEO, AnAge, SenNet). Manifest JSONs must record provenance, file sizes, and cryptographic SHA-256 hashes.

---

### Perspective E: Critical Skeptic & Adversarial Peer Reviewer
#### Round 9: The Danger of Confirmation Bias in "Conserved" Aging
- **Critique:** Single-cell literature is plagued by confirmation bias: researchers search for genes that go up in both species, find 50 genes (mostly inflammatory or ribosomal), and declare a "conserved universal aging signature." Does this plan guard against false positives?
- **Finding:** The plan incorporates three mandatory counterweights:
  1. **Strict Directional Spearman Correlation:** Testing correlation across ALL expressed orthologs, not just the pre-selected DE subset.
  2. **Donor Age Permutation Null:** Calculating the empirical null distribution of conservation scores to demonstrate that the observed conservation exceeds chance expectation.
  3. **First-Class Divergence Reporting:** Explicitly highlighting genes and cell types where mouse and human age differently, presenting divergence as an actionable biological insight rather than hiding it.

#### Round 10: Falsification Crucible for Longevity Outlier Claims
- **Critique:** The popular narrative claims naked mole-rats "do not age." If a study tests a gene set in NMR and finds $p > 0.05$, declaring "the program is uncoupled" is the elementary fallacy of treating failure to reject the null as proof of the null.
- **Finding:** The plan forbids using $p > 0.05$ as evidence of absence. Instead, it mandates **Two One-Sided Tests (TOST)** with a pre-registered equivalence margin $\delta$. Equivalence is only claimed if the 90% confidence interval of the NMR age slope falls entirely within $[-\delta, +\delta]$. If statistical power is insufficient to bound the slope within $\delta$, the verdict is strictly reported as "inconclusive / underpowered," never "conserved uncoupling."

---

## Part 3: Revised Implementation & Execution Blueprint

1. **Environment & Core Infrastructure:** Completed in M0 (`environment.lock`, `scripts/verify_data.py`, `data/manifest/*.json`, `src/`).
2. **Deterministic Analysis Pipeline (M1–M6):** Run via modular Python scripts and fully executed Jupyter notebooks with `src.set_seeds(42)`, single-threaded UMAP, sparse matrix linear algebra, mouse-level pseudobulk PyDESeq2, Ensembl 1:1 orthology, and TOST equivalence testing.
3. **Synthesis & Polish (M7):** Generate `report/conserved-aging-programs.md` with explicit limitations, deviations, and golden-output diff verification via `scripts/run_all.py`.
4. **Interactive Static Explorer (M8):** Deploy modern static web app with DuckDB-WASM/Arrow data loading, curated gene panels, and publication-ready SVG/PNG exports.
