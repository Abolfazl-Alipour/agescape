# Agescape App — Plot Data Contract (Step 0 output)

**Date:** 2026-10-18. **Spec:** `docs/webapp-spec-v2.md` (v2.0). **Status:** Step 0 complete — all Section-10 anchors verified against committed files; two spec-value corrections noted below (files win). Build proceeds on these values.

## Step 0 audit — every Section-10 anchor vs. the files

| Anchor (spec) | File value | Verdict |
|---|---|---|
| 163 units evaluated, 17 powered | `m2-unit-eligibility.csv` 163 rows; `m2-aging-signal-ranked.csv` 17 rows | ✓ |
| Powered confined to Lung & Kidney | ranked units → tissues {Kidney, Lung} only | ✓ |
| Ranking head 55, 53, 36, 35 | Lung::adventitial 55, Kidney::PCT 53, Lung::fibroblast 36, Lung::alveolar macrophage 35 | ✓ |
| Calibration 55 / 4.6 / 0 | `m2-calibration-trio-summary.csv`: observed 55, permutation mean 4.6, sex floor 0 | ✓ |
| Agreement 13 agree / 4 partial / 0 disagree | `m3-agreement-table.csv` verdict counts | ✓ |
| Kidney PCT rho 0.866 | 0.866, n_shared 13,019 | ✓ |
| Program 167 up / 94 down / 312 gene-unit entries | `m6-program-genes.csv` 167/94; `m5-gene-level-programs.csv` 312 rows | ✓ |
| Mouse whole-body slope 0.235 [0.067, 0.403] p=0.008, n=23 | 0.23496 [0.06727, 0.40265] p=0.0083, n=23 | ✓ |
| NMR slope 0.017 [-0.257, 0.290] p=0.858, n=5 | 0.01667 [-0.25667, 0.290] p=0.8585, n=5 | ✓ |
| Slope difference 0.218 [-0.013, 0.449] p=0.064 | 0.21828 [-0.01263, 0.44920] p=0.0639 | ✓ |
| TOST δ=0.0587, p=0.329 not established | δ=0.05874, p=0.32894 | ✓ |
| Exact permutation p=0.818 | 0.81818 | ✓ |
| Coverage 1.00 / 0.73 / 0.31 (mouse/human/NMR) | `m6-gates.csv` 1.0 / 0.73180 / 0.31414 | ✓ |
| Power anchors: 1.474 @ (4 donors, SD 8, disp 0.1); 0.59 @ (4, 20, disp 0.1) | `m4-human-power.csv` `min_detectable_30y_log2fc`: 1.474 and 0.59 | ✓ |
| Savina 29 @ stated rule; 9 & 6 @ 1.45; 560 @ loose rule | `m6-savina-reproduction.csv`: @1.0 → 29 (p<0.05), 20 (FDR); @1.45 → 9 (p), 6 (FDR); @0.25 → 811 (p), 560 (FDR) | ✓ |
| Lifespans 4 / 31 / 122.5 y | `m6-species-lifespans.csv` mouse 4.0, NMR 31, human 122.5 | ✓ |

**Corrections to the spec (files win; Reza informed per Step-0 rule):**
1. Spec Table Ch4 says mouse coverage "1.00, human 0.73, NMR 0.31" — verified; no change needed. (No numeric disagreement found anywhere else; the audit confirmed the spec's numbers.)
2. `m4-human-power.csv` contains BOTH `min_detectable_slope_per_year` and `min_detectable_30y_log2fc`; the Power Map reads only the 30y column (spec-correct), slope-per-year column is not displayed.
3. Power grid axes: donors {4,6,8,10,12,16,20,30,50} × SD {8,12,16,20} × dispersion {0.1,0.2,0.4} = 108 tiles per dispersion panel (spec said "as present in the grid" — present = 108).

## Conditional export verdicts (Section 13 decision #1)

- **C1 permutation replicates: NOT AVAILABLE.** No committed file retains per-replicate DE counts (pipeline stored only the mean 4.6 ± 2.1). Producing them requires an analysis re-run — **deferred to Reza** per spec default. → **G2 ships in fallback dot-plot form.**
- **C2 per-gene DE for Lung::adventitial cell: EXISTS** — `results/m2_de_Lung__adventitial_cell.csv` (10,150 genes; index, baseMean, log2FC, lfcSE, stat, pvalue, FDR, CIs). Pure transcription. → **G8 volcano unlocked.**
- **C3 paired effect sizes for Kidney::PCT: PRODUCIBLE BY PURE EXPORT** — our **continuous bridge slope per month** (`results/m3-continuous-bridge-model.csv`, the M3 comparator that `scripts/run_m3_agreement.py` prefers when >200 bridge genes exist) joined to the paper's published age coefficient (`data/raw/tms_nature2020_moesm8_dge_results.xlsx`, sheet `DGE_result.tissue_cell.age_coef`, column `Kidney.kidney proximal convoluted tubule epithelial cell`; bh_p from the sibling sheet). Both committed; the join is a reshape, no new analysis. → **G3 scatter unlocked as chapter lead.**
  - **Correction (2026-10-19, file wins):** this contract originally named the m2 dichotomous log2FC file as the "ours" side. The literal m2_de join gives n=11,836, rho=0.309 and does NOT reproduce the pinned anchor (n=13,019, rho=0.866). The anchor is defined by the bridge-slope join; the export prints the literal join to stdout as a diagnostic only. Export field renamed `ours_log2fc` → `ours_bridge_slope` to match what the values actually are.

## Stack & licence decisions (Step 0, per spec authority)

- **No Node.js on the build machine** (verified) → Astro impossible. **Equivalent settled: zero-build static site** — one `index.html` + ES-module JS + CSS, all charts **ECharts 5.5.1 vendored locally** at `app/vendor/echarts.min.js` (Apache-2.0; licence retained in file header). No runtime external calls.
- **Licences (spec working assumptions adopted):** code MIT; report/figures CC BY 4.0. Licence lines go in the Methods overlay and a top-level LICENCE.md.
- **Analytics: NONE** (spec default — cookieless option not exercised).
- **Base path:** configurable, default relative paths (works at `/`, `/agescape/`, or `file://`).

## Required exports (schemas pinned by this contract)

All written by `scripts/export_app_data.py` → `app/data/`, deterministic, JSON only, total < 5 MB. Freeze date recorded per export: **2026-10-14** (data) / export date 2026-10-18.

1. **`atlas.json`** ← m2-unit-eligibility.csv (d394b440aba5) + m2-aging-signal-ranked.csv (fb000922b087)
   `{units: [{unit, tissue, cell_ontology_class, eligible, powered, n_cells_3m, n_cells_18m, n_mice_3m, n_mice_18m}], ranked: [{unit, n_de_genes, median_abs_log2fc, n_de_ci_low, n_de_ci_high, tissue}]}` — 163 units, 17 ranked. Powers G0 (163 tiles, per-tissue "X of Y" counts) + G1 (bars + whiskers).
2. **`calibration.json`** ← m2-calibration-trio-summary.csv (ac3e65b60721) + m2-calibration-positive-controls.csv (0c78f4ae25f6)
   `{observed: 55, perm_mean: 4.6, sex_floor: 0, unit: "Lung::adventitial cell", positive_controls: [...]}` — powers G2 fallback dot-plot. (Histogram form blocked on C1 — deferred.)
3. **`agreement.json`** ← m3-agreement-table.csv (743ae1f1f262) + C3 join
   `{lollipop: [{unit, spearman_rho, verdict}], counts: {agree: 13, partial: 4, disagree: 0}, scatter: {unit: "Kidney::kidney proximal convoluted tubule epithelial cell", n: 13019, rho: 0.866, points: [{gene, ours_bridge_slope, paper_coef}]}}` — scatter points only for genes present in BOTH sources; `ours_bridge_slope` is the M3 continuous bridge slope per month (NOT the dichotomous m2 log2FC; see C3 correction above); rho asserted = 0.866 (recomputed at export for verification only). Powers G3 both layers.
4. **`design-grid.json`** ← m1-sex-age-table.csv (832e990f8803)
   `{age_groups: ["1m","3m","18m","21m","24m","30m"], cells: {female: [...6], male: [...6]}, mice: {female: [...6], male: [...6]}}` — powers G4; hatch rule: female-only columns (21m), male-only (1m, 24m, 30m); outline 3m & 18m as "the only fair fight".
5. **`program.json`** ← m5-mouse-program-catalog.csv (d5ab9ace5f3c) + m5-gene-level-programs.csv (5251165002a1) + m6-program-genes.csv (941766a85207)
   `{units: [{unit, n_up, n_down, frac_ribo_mito, frac_ribo_mito_background}], genes: [{gene, unit, direction, log2fc, fdr}], summary: {up: 167, down: 94, entries: 312}, matrix_top40: {rows: 40 (ranked by units-present, full set in `genes`), columns: 17 units}}` — powers G6, G7, unit-panel gene tables. **Stated sampling rule:** G7 displays top 40 genes by unit-recurrence; the full 312-row matrix is one tap down (drawer download). Not silent dropping.
6. **`signature.json`** ← m6-lifespan-program-scores.csv (dfff50c3d1d4) + m6-statistical-summary.csv (176a8dcab93f) + m6-gates.csv (1c1260534996) + m6-species-lifespans.csv (ad28535a0957)
   `{arms: {mouse_wb: {points: [{donor, age_years, age_normalized, score, n_cells}], n: 23, slope, ci, p}, mouse_skin: {...n:3, underpowered: true}, human_skin: {...n:4, low_powered: true}, human_fibroblast: {...}, nmr: {..., tested_span_years: [1.1, 11.3]}}, stats: {diff: {...}, tost: {delta: 0.05874, p: 0.32894, established: false}, perm_p: 0.81818}, coverage: {mouse: 1.0, human: 0.73180, nmr: 0.31414}, lifespans: {mouse: 4.0, nmr: 31.0, human: 122.5}}` — powers G9, G10, G11.
7. **`savina.json`** ← m6-savina-reproduction.csv (cfdab4a0f3a1)
   `{curve: [{log2fc_cutoff, n_de_pval05, n_de_fdr05}], markers: {paper_rule: {cutoff: 1.0, n_pval: 29}, converted_rule: {cutoff: 1.45, n_pval: 9, n_fdr: 6}, loose_rule: {cutoff: 0.25, n_fdr: 560}}, caveats: ["Igfbp3 absent from source file", "QC 10,251 vs paper 10,232"]}` — powers G12.
8. **`power-grid.json`** <- m4-human-power.csv (03dc469a5a05)
   `{donors: [4,6,8,10,12,16,20,30,50], sds: [8,12,16,20], dispersions: [0.1,0.2,0.4], tiles: [{n_donors, sd, dispersion, min_30y_log2fc, power_status, eligible_cohort}], ts_marker: {donors: 4, sd: 8, label: "Tabula Sapiens skin: 4 donors"}, anchors: {"4_8_low": 1.474, "4_20_low": 0.59}}` : powers G13. Readout may display only grid values.
9. **`provenance.json`** ← computed at export
   `{freeze_date: "2026-10-14", export_date, sources: {per-export: {file, sha256 (full), rows}}, assets: {per-file: {bytes}}}` — feeds every provenance drawer. SHA-256 short forms in drawers = first 12 hex chars of these.

## Figure assets (for compare-toggles & G5 & fallback images)

Copied verbatim from `figures/` (no re-rendering): `m3-dpt-vs-age.png` (17e50141a309), `m6-signature-figure.png` (d37d0ef831c4), `m6-positive-control-skin.png` (9371628bf60e), `m6-savina-reproduction.png` (ddc8ccee4034).

## Hard gates

- `scripts/check_app_data.py` asserts every anchor in the table above against the exports and exits non-zero on any miss. Run after every export; the page build consumes only checked exports.
- No export may contain a value not derivable from the listed sources. The C3 scatter is the only multi-source export: an inner join on gene symbol between the M3 continuous bridge slope (per month) and the paper's age coefficient, mirroring `scripts/run_m3_agreement.py`.
- Word/encoding budget enforcement (Section 5 text ≤700 words, ≥70% graph area) is a Step-5 checklist item, not an export concern.
