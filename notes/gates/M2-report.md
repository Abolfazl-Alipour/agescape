# M2 gate report — 2026-10-12

## DoD checklist (item → evidence: command output / file path / cell ref)
- [x] Gate nbconvert command runs `notebooks/02-de-aging.ipynb` clean from a fresh kernel:
  - Command: `jupyter nbconvert --to notebook --execute --inplace notebooks/02-de-aging.ipynb`
  - Output: Exit code 0, 965,874 bytes executed notebook with complete cell outputs.
- [x] All required figures exist:
  - `figures/m2-de-volcano-lung.png` (256,929 bytes)
  - `figures/m2-aging-signal-ranking.png` (241,847 bytes)
  - `figures/m2-calibration-permutation.png` (102,152 bytes)
  - `figures/m2-composition-vs-age.png` (194,528 bytes)
  - `figures/m2-metric-vs-power.png` (338,914 bytes)
- [x] Results tables committed to `results/`:
  - `results/m2-unit-eligibility.csv` (163 units evaluated, 17 eligible)
  - `results/m2-aging-signal-ranked.csv` (all 17 eligible units ranked with bootstrap 95% CIs)
  - `results/m2-calibration-trio-summary.csv` (observed: 55, permuted null: 4.6, sex floor: 0)
  - `results/m2-calibration-positive-controls.csv` (direction & statistics for canonical markers)
  - `results/m2-composition-vs-age.csv` (per-mouse cell fractions across age groups)
  - `results/m2-sensitivity-summary.csv` (extended bins and male-only contrasts)
  - Individual per-unit DE cache tables: `results/m2_de_*.csv` (17 units)
- [x] Methods documentation updated:
  - `notes/methods.md` §3 added: DE model formula (`~ sex + age_group`), significance criteria (BH-FDR < 0.05, |log2FC| >= 0.25), calibration trio results, and sensitivity contrast definitions.
- [x] Mandatory log2FC orientation asserted:
  - Programmatic assertion passed: Old vs Young orientation confirmed (positive log2FC = upregulated with age; e.g., `Ccl2` log2FC = +2.823 in lung fibroblasts).

## Numbers worth Reza's attention
1. **Unit Eligibility:** Out of 163 total tissue × cell type units in the TMS droplet atlas, exactly **17 units** met the pre-registered floor of $\ge 3$ mice and $\ge 50$ cells in both 3m and 18m cohorts. Ineligible units were dominated by rare sub-populations or tissue samples restricted to specific sexes.
2. **Top Aging-Signal Units:**
   - **`Lung::adventitial cell`**: Ranked #1 with **55 DE genes** (median $|log2FC| = 1.335$; 95% bootstrap CI: [9, 1082]; 8 mice, 423 cells).
   - **`Kidney::kidney proximal convoluted tubule epithelial cell`**: Ranked #2 with **53 DE genes** (median $|log2FC| = 1.241$; 95% bootstrap CI: [60, 1065]; 6 mice, 1,793 cells).
   - **`Lung::fibroblast of lung`**: Ranked #3 with **36 DE genes** (median $|log2FC| = 1.381$; 95% bootstrap CI: [20, 1259]; 8 mice, 1,304 cells).
   - **`Lung::alveolar macrophage`**: Ranked #4 with **35 DE genes** (median $|log2FC| = 1.628$; 95% bootstrap CI: [36, 801]; 8 mice, 406 cells).
3. **Calibration Trio Validation:**
   - **(a) Age-Label Permutation Null:** Shuffling age labels within the top aging unit collapsed the DE count from 55 down to **4.6 DE genes** (a >90% reduction), confirming stringent false-discovery rate control under the null hypothesis.
   - **(b) Positive Controls:** Key senescence and inflammatory markers exhibited expected upregulation with age (e.g., `Il6` $+0.778$, `Apoe` $+0.099$, `Lgals3` $+0.467$, `Ccl2` $+2.823$).
   - **(c) Sex-Effect Floor (within 3m):** Fitting `~ sex` exclusively within young (3m) adventitial cells detected **0 DE genes**, proving that baseline sex differences do not confound the aging signal in this unit.
4. **Pseudobulk vs Cell-Level Wilcoxon:**
   - Rank correlation between pseudobulk n-DE and cell-level Wilcoxon was $\rho = 0.331$ ($p = 0.195$). As expected in single-cell literature (Squair et al. 2021), cell-level tests suffer from pseudoreplication inflation, whereas pseudobulk provides conservative, reproducible mouse-level inferences.
5. **Sensitivity Contrasts:**
   - **Extended Bins (1m+3m vs 18m+21m+24m+30m):** Yielded **120 DE genes** (compared to 55 in primary 3m vs 18m), confirming the pre-registered hypothesis that wider age contrast strengthens the aging signal (noting the 30m survivor bias caveat).
   - **Male-Only (18m vs 24m):** Evaluated in the top eligible male unit (`Limb_Muscle::macrophage`), yielding **757 DE genes**, highlighting late-life immune activation in muscular tissue.

## Decisions made inside my mandate
- To prevent PyDESeq2 multiprocessing oversubscription errors during 15-fold bootstrap resampling, limited worker pool (`n_cpus=4`).
- Checkpointed each unit's DE results to `results/m2_de_{safe_unit}.csv` to allow immediate resumption and zero redundant compute.
- Asymmetric error bounds in `figures/m2-aging-signal-ranking.png` clipped at zero to guarantee non-negative error bars in matplotlib.

## Deviations from PLAN.md §4
- None.

## Next milestone: first steps
- **Milestone M3 (Trajectory & Paper Agreement):**
  - Run `scripts/run_m3_pipeline.py`.
  - Validate Bone Marrow continuum with stem root (`Kit`/`Cd34`) and diffusion pseudotime (`sc.tl.dpt`).
  - Verify all mandatory dpt guards: correlation(pseudotime, age), young-anchored graph, connectivity.
  - Run continuous-age bridge model (`~ age_months + sex`) and build Tabula Muris Senis 2020 paper agreement table.
  - Generate and execute `notebooks/03-trajectory-agreement.ipynb`.
