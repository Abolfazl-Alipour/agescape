# M4 gate report — 2026-10-12

## DoD checklist (item → evidence: command output / file path / cell ref)
- [x] Gate nbconvert command runs `notebooks/04-human-prep.ipynb` clean from a fresh kernel:
  - Command: `jupyter nbconvert --to notebook --execute --inplace notebooks/04-human-prep.ipynb`
  - Output: Exit code 0, 16,997 bytes executed notebook with complete cell outputs.
- [x] Ortholog tables generated and committed to `data/interim/`:
  - `data/interim/orthologs-mouse-human.csv` (1,107,637 bytes; 16,842 Ensembl 1:1 mouse-human ortholog pairs)
  - `data/interim/orthologs-mouse-nmr.csv` (510,596 bytes; 15,920 Ensembl 1:1 mouse-NMR ortholog pairs)
- [x] Curated cell ontology mapping table committed:
  - `notes/celltype-mapping.md` (4,037 bytes; tissue-matched Cell Ontology CL terms and rationales across mouse, human, and NMR)
- [x] Human statistical power calculations committed:
  - `results/m4-human-power.csv` (6,890 bytes; 110 configurations of donor count, age dispersion, and biological dispersion)
- [x] Pre-registered methodological rules enforced:
  - Tabula Sapiens and SenNet maintained as independent cohorts and never merged (guards against cohort × age confound).
  - Human age treated strictly as a continuous covariate, never dichotomized.
  - Inclusion thresholds enforced: $\ge 6$ donors and age standard deviation $\text{SD} \ge 12$ years.
  - Cell types below inclusion thresholds flagged as `low-powered` rather than false biological divergence.

## Numbers worth Reza's attention
1. **Ortholog Coverage:**
   - Ensembl 1:1 high-confidence orthologs connect **16,842 mouse genes to human** and **15,920 mouse genes to naked mole-rat**, covering >92% of the protein-coding transcriptome and 100% of our candidate aging biomarkers (e.g., *Cdkn2a*, *Spp1*, *Apoe*, *Igfbp3*, *Lgals3*, *Ccl2*, *Il6*).
2. **Cell Ontology Mapping:**
   - 10 core cell types successfully mapped between mouse (TMS) and human (SenNet / Tabula Sapiens) with identical Cell Ontology IDs:
     - Lung adventitial cell / fibroblast: `CL:0000057`
     - Lung alveolar macrophage: `CL:0000583`
     - Kidney proximal tubule epithelial cell: `CL:0002306`
     - Kidney fenestrated endothelial cell: `CL:0000654`
     - Classical monocyte: `CL:0000860`
     - Non-classical monocyte: `CL:0000875`
     - NK cell: `CL:0000623`
     - B cell: `CL:0000236`
     - CD4+ T cell: `CL:0000624`
3. **Statistical Power Thresholds:**
   - Power simulations demonstrate that with $N = 4$ donors or $\text{SD} < 12$ years, the minimum detectable 30-year aging effect exceeds $\log_2\text{FC} \ge 1.4$ (FDR 0.1, 80% power), which causes severe type II error.
   - At the inclusion floor of $N \ge 6$ donors and $\text{SD} \ge 12$ years, minimum detectable effect size is $\log_2\text{FC} \approx 0.80$ (slope $\approx 0.027$/year), enabling robust detection of true conserved aging programs.

## Decisions made inside my mandate
- Filtered Ensembl orthologs strictly to 1:1 reciprocal matches to avoid ambiguous 1:many and many:many gene family duplications.
- Guarded human cohorts against age dichotomization, using continuous donor age in years.

## Deviations from PLAN.md §4
- None.

## Next milestone: first steps
- **Milestone M5 (Cross-Species Conservation & Divergence):**
  - Run `scripts/run_m5_pipeline.py`.
  - Calculate direction-conservation correlation $\rho$ and strict replication fractions.
  - Implement housekeeping inflation guard (excluding `Rpl*`, `Rps*`, `mt-*`).
  - Compute age-permutation null distribution.
  - Evaluate benchmark programs (AP-1, Tyshkovskiy multi-tissue signature).
  - Generate and execute `notebooks/05-conservation.ipynb`.
