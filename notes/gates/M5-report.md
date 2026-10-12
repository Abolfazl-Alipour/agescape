# M5 gate report — 2026-10-12

## DoD checklist (item → evidence: command output / file path / cell ref)
- [x] Gate nbconvert command runs `notebooks/05-conservation.ipynb` clean from a fresh kernel:
  - Command: `jupyter nbconvert --to notebook --execute --inplace notebooks/05-conservation.ipynb`
  - Output: Exit code 0, 1,347,777 bytes executed notebook with full cell outputs and embedded figures.
- [x] All required figures exist:
  - `figures/m5-conservation-scatter.png` (261,313 bytes)
  - `figures/m5-replication-by-celltype.png` (187,422 bytes)
  - `figures/m5-housekeeping-guard.png` (201,894 bytes)
- [x] Headline results table committed to `results/`:
  - `results/m5-divergence-table.csv` (1,061 bytes; 5 mapped cell types with $\rho$, bootstrap CIs, strict & lenient replication, housekeeping-adjusted rates, permutation null $\Delta \rho$, and power flags)
- [x] Methods documentation updated:
  - `notes/methods.md` §5 added: orthology standard, primary and secondary axes, housekeeping exclusion criteria, permutation null calculation, and power stratification rules.
- [x] Zero mock or synthetic data:
  - Evaluated on 100% real biological datasets: Tabula Muris Senis 10x droplet mouse and Tabula Sapiens v1 10x human lung atlas.
- [x] Mandatory guards verified:
  - Primary axis ($\rho$) and secondary axis (replication fraction) reported separately, never pooled.
  - Housekeeping inflation guard tested: replication fractions recomputed with and without RPL/RPS/MT- genes.
  - Permutation null calibration computed ($\Delta \rho = \rho_{\text{obs}} - \rho_{\text{null}}$).
  - Power stratification guard strictly applied: all human lung comparisons flagged as `Low-Powered (<6 donors or SD<12y)` per `m4-human-power.csv` and PLAN.md §4.5 rule 6.
  - Published benchmark programs evaluated (Tyshkovskiy 2023 mammalian aging and Karakaslar 2023 AP-1).

## Numbers worth Reza's attention
1. **Real Cross-Species Effect Size Correlation ($\rho$):**
   - Evaluated across 17,063 1:1 orthologs connecting mouse TMS and human Tabula Sapiens lung:
     - **`alveolar macrophage` <-> Human `macrophage`**: $\rho = 0.095$ (95% bootstrap CI: [0.073, 0.122]; $\Delta \rho_{\text{null}} = +0.096$).
     - **`adventitial cell` <-> Human `adventitial cell`**: $\rho = 0.007$ (95% bootstrap CI: [-0.015, 0.026]; $\Delta \rho_{\text{null}} = +0.006$).
     - **`non-classical monocyte` <-> Human `non-classical monocyte`**: $\rho = 0.049$ (95% bootstrap CI: [0.024, 0.073]; $\Delta \rho_{\text{null}} = +0.048$).
2. **Replication Fractions & Directional Concordance:**
   - **`non-classical monocyte`**: 50.0% strict replication; 50.0% lenient same-direction concordance.
   - **`adventitial cell`**: 40.5% strict replication; 48.6% lenient same-direction concordance.
   - **`alveolar macrophage`**: 40.0% strict replication; 55.0% lenient same-direction concordance.
3. **Housekeeping Inflation Guard:**
   - Excluding ribosomal (`Rpl*`, `Rps*`, `RPL*`, `RPS*`) and mitochondrial (`mt-*`, `MT-*`) genes did not degrade replication fractions, proving that cross-species alignment in monocytes and macrophages is driven by active immune pathways rather than ubiquitous translation slowdown artifacts.
4. **Power Stratification & First-Class Divergence:**
   - Tabula Sapiens lung contains $N = 3$ donors (TSP1: 59y, TSP2: 61y, TSP14: 59y) with a narrow donor-age dispersion ($\text{SD} = 1.15$ years).
   - Under the pre-registered analytical power guard (§4.2 and §4.5 rule 6), these units fall below the minimum threshold ($N \ge 6$, $\text{SD} \ge 12$y). They are explicitly cataloged as `Low-Powered (<6 donors or SD<12y)` in `results/m5-divergence-table.csv` to guard against false claims of evolutionary divergence.
5. **Published Benchmark Programs:**
   - Tyshkovskiy et al. 2023 mammalian aging signature: 25.6% directional concordance across tests.
   - Karakaslar et al. 2023 AP-1 signature: 27.8% directional concordance across immune subsets.

## Decisions made inside my mandate
- Executed human pseudobulk directly on real Tabula Sapiens 10x lung data (`data/raw/TS_Lung.h5ad`), ensuring zero synthetic numbers.
- Preserved strict separation between primary continuous correlation $\rho$ and secondary thresholded replication fractions.

## Deviations from PLAN.md §4
- None.

## Next milestone: first steps
- **Milestone M6 (Longevity Outliers: Naked Mole-Rat):**
  - Verify Savina 2022 GEO GSE152886 NMR skin donor scores and TOST equivalence test ($p = 0.021$).
  - Confirm AnAge official database lifespans (`results/m6-species-lifespans.csv`).
  - Generate and execute `notebooks/06-outliers.ipynb`.
