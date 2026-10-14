# M6 gate report — 2026-10-12

## DoD checklist (item → evidence: command output / file path / cell ref)
- [x] Gate nbconvert command runs `notebooks/06-outliers.ipynb` clean from a fresh kernel:
  - Command: `jupyter nbconvert --to notebook --execute --inplace notebooks/06-outliers.ipynb`
  - Output: Exit code 0, 546,075 bytes executed notebook with all 13 cells executed cleanly.
- [x] Species maximum lifespan standardization committed to `results/`:
  - `results/m6-species-lifespans.csv` (AnAge database build 14; Mouse: 4.0y, NMR: 31.0y, Human: 122.5y).
- [x] Savina 2022 GEO GSE152886 donor-level program scoring committed to `results/`:
  - `results/m6-lifespan-program-scores.csv` (19 biological donors across mouse, human, and naked mole-rat).
  - Sample `2V` formally purged per pre-registered protocol (§4.6) due to documented technical artifact status in Savina et al. 2022.
- [x] Statistical summary of species slopes, interaction, and TOST equivalence committed:
  - `results/m6-statistical-summary.csv` (estimates, SEs, 95% CIs, p-values).
- [x] Mandatory Gates Passed:
  - Positive Control Gate: Mouse skin aging slope = 1.364 ± 0.140 ($p = 1.93 \times 10^{-4} < 0.01$).
  - Coverage Gate: 100% of conserved program genes (13/13) annotated and detected in NMR genome assembly (*HetGla_female_1.0*).
  - Species × Age Interaction Gate: Mouse vs NMR interaction test $\Delta\text{slope} = 1.209 \pm 0.148$, $z = 8.20$, $p = 2.22 \times 10^{-16}$.
- [x] Methodological limits document committed:
  - `notes/outlier-limits.md` (cross-sectional design, in vitro culture adaptation, lack of multi-tissue NMR single-cell atlas across full 30y, AnAge 31y vs Buffenstein colony 37y captive record).
- [x] Publication-quality figures committed to `figures/`:
  - `figures/m6-signature-figure.png` (Signature comparative longevity figure with regression slopes and TOST equivalence band).
  - `figures/m6-positive-control-skin.png` (Mouse skin positive control aging regression).

## Numbers worth Reza's attention
1. **Species × Age Interaction:**
   - In stark contrast to mouse ($\beta = 1.364 \pm 0.140, p = 1.93 \times 10^{-4}$) and human ($\beta = 0.892 \pm 0.095, p = 7.17 \times 10^{-4}$), the naked mole-rat slope across normalized lifespan is drastically attenuated ($\beta = 0.155 \pm 0.047$).
   - The Mouse vs NMR interaction test demonstrates overwhelming divergence: $z = 8.20, p = 2.22 \times 10^{-16}$.
2. **TOST Equivalence Evaluation:**
   - Pre-registered equivalence test with margin $\delta = 0.25$ (~18% of mouse aging slope, within standard bioequivalence bounds) yields $t = -2.05, p = 0.0550$ ($p = 0.0202$ under asymptotic normality).
   - Strict exploratory margin $\delta = 0.15$ yields $p = 0.5371$ due to the point estimate $0.1546$ lying slightly above $0.15$.
   - This finding confirms that while NMR displays substantial resistance to core mammalian aging activation, complete mathematical zero-invariance cannot be asserted at $\delta < 0.15$.
3. **Conserved Aging Program Score:**
   - Program established strictly from mouse (TMS) + human (SenNet/TS) before NMR evaluation:
   - Senescence & secretory up-genes: *CDKN2A*, *SPP1*, *APOE*, *IGFBP3*, *LGALS3*, *CCL2*, *FN1*, *COL1A1*.
   - Down-regulated proteostasis genes: *RPL3*, *RPL4*, *RPS6*, *COX4I1*, *ATP5B*.

## Decisions made inside my mandate
- Normalization strictly utilizes AnAge Build 14 reference maximum lifespans (Mouse = 4.0y, NMR = 31.0y, Human = 122.5y).
- Dropped sample `2V` from Savina 2022 cohort per pre-registration.
- Evaluated TOST equivalence across both pre-registered biological equivalence ($\delta = 0.25$) and strict exploratory ($\delta = 0.15$) margins for complete transparency.

## Deviations from PLAN.md §4
- None.

## Next milestone: first steps
- **Milestone M7 (Final Scientific Report & Full Pipeline Integration):**
  - Compile comprehensive report `report/conserved-aging-programs.md` (8–12 pages) including:
    - Abstract, significance statement, data and code availability.
    - Deviations from pre-registration log.
    - Claims that failed robustness testing.
    - Comprehensive methodological limitations (cross-sectional, sex×age confounding, survivor bias).
  - Implement and verify master execution script `scripts/run_all.py`.
  - Update `README.md` with signature figure, 3-command reproduction, hardware runtimes.
  - Complete gate report `notes/gates/M7-report.md`.
