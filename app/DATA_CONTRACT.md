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
| Kidney PCT rho 0.569 | 0.568639, n_shared 12,713 | ✓ |
| Program 166 up / 91 down / 312 gene-unit entries | `m6-program-genes.csv` 166/91; `m5-gene-level-programs.csv` 312 rows | ✓ |
| Mouse whole-body slope 0.1344 [0.0293, 0.2395] p=0.0147, n=23 | 0.13437 [0.02926, 0.23949] p=0.014700, n=23 | ✓ |
| NMR slope 0.0058 [-0.1981, 0.2098] p=0.9332, n=5 | 0.005836 [-0.198122, 0.209795] p=0.933179, n=5 | ✓ |
| Slope difference 0.1285 [-0.0314, 0.2885] p=0.1153 | 0.128536 [-0.031443, 0.288515] p=0.115309 | ✓ |
| TOST δ=0.2023, p=0.02748 at a post hoc margin exceeding the mouse estimate | δ=0.2023, p=0.027479 | ✓ |
| Exact permutation p=1.000 | 1.000 | ✓ |
| Coverage 1.00 / 0.8489 / 0.6619 (mouse/human/NMR) | `m6-gates.csv` 1.0 / 0.84892 / 0.66187 | ✓ |
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

- **No Node.js on the build machine** (verified) → Astro impossible. **Equivalent settled: zero-build static site** — one `index.html` + ES-module JS + CSS, all charts **ECharts 5.5.1 vendored locally** at `app/vendor/echarts.min.js` (Apache-2.0; licence retained in file header). Data and charts make no runtime external calls; Google Fonts currently loads its stylesheet and font files externally.
- **Licences (spec working assumptions adopted):** code MIT; report/figures CC BY 4.0. Licence lines go in the Methods overlay and a top-level LICENCE.md.
- **Analytics: NONE** (spec default — cookieless option not exercised).
- **Base path:** configurable relative paths, working at `/` and `/agescape/`. Serve over HTTP because browser `fetch()` does not reliably permit `file://` JSON loading.

## Required exports (schemas pinned by this contract)

All written by `scripts/export_app_data.py` → `app/data/`, deterministic, JSON only, total < 5 MB. Freeze date recorded per export: **2026-10-14** (data) / export date 2026-10-18.

1. **`atlas.json`** ← m2-unit-eligibility.csv (d394b440aba5) + m2-aging-signal-ranked.csv (fb000922b087)
   `{units: [{unit, tissue, cell_ontology_class, eligible, powered, n_cells_3m, n_cells_18m, n_mice_3m, n_mice_18m}], ranked: [{unit, n_de_genes, median_abs_log2fc, n_de_ci_low, n_de_ci_high, tissue}]}` — 163 units, 17 ranked. Powers G0 (all 163 tiles, per-tissue "X of Y" counts) + G1 (ranked bars). Exported bootstrap fields are retained but not drawn as confidence intervals because 15 resamples are insufficient for stable tail estimates.
2. **`calibration.json`** ← m2-calibration-trio-summary.csv (ac3e65b60721) + m2-calibration-positive-controls.csv (0c78f4ae25f6)
   `{observed: 55, perm_mean: 4.6, sex_floor: 0, unit: "Lung::adventitial cell", positive_controls: [...]}` — powers G2 fallback dot-plot. (Histogram form blocked on C1 — deferred.)
3. **`agreement.json`** ← m3-agreement-table.csv (743ae1f1f262) + C3 join
   `{lollipop: [{unit, spearman_rho, sign_concordance, verdict}], counts: {agree: 13, partial: 4, disagree: 0}, scatter: {unit: "Kidney::kidney proximal convoluted tubule epithelial cell", n: 13019, rho: 0.866, points: [{gene, ours_bridge_slope, paper_coef}]}}` — scatter points only for genes present in BOTH sources; `ours_bridge_slope` is the M3 continuous bridge slope per month (NOT the dichotomous m2 log2FC; see C3 correction above); rho asserted = 0.866 (recomputed at export for verification only). Powers G3 both layers.
4. **`design-grid.json`** ← m1-sex-age-table.csv (832e990f8803)
   `{age_groups: ["1m","3m","18m","21m","24m","30m"], cells: {female: [...6], male: [...6]}, mice: {female: [...6], male: [...6]}}` — powers G4; hatch rule: female-only columns (21m), male-only (1m, 24m, 30m); outline 3m & 18m as "the only fair fight".
5. **`program.json`** ← m5-mouse-program-catalog.csv (d5ab9ace5f3c) + m5-gene-level-programs.csv (5251165002a1) + m6-program-genes.csv (941766a85207)
   `{units: [{unit, n_up, n_down, frac_ribo_mito, frac_ribo_mito_background}], genes: [{gene, unit, direction, log2fc, fdr}], summary: {up: 166, down: 91, entries: 312}, matrix_top40: {rows: 40 (ranked by units-present, full set in `genes`), columns: 17 units}}` — powers G6, G7, unit-panel gene tables. **Stated sampling rule:** G7 displays top 40 genes by unit-recurrence; the full 312-row matrix is one tap down (drawer download). Not silent dropping.
   - **Correction (2026-10-19, file wins):** the m5 catalog stores m2-DE-filename slugs (`Lung::adventitial_cell`); all `unit` fields and matrix columns in this export are translated back to display names (`Lung::adventitial cell`) so program.json joins to every other export (and to the app's panel lookups) on identical strings. Translation rule: filename stem form `re.sub(r"[^A-Za-z0-9_]+", "_", unit.replace("::", "__"))`, keyed back with `__` restored to `::`; asserted 1:1 against the 17 committed `m2_de_*.csv` stems at export.
6. **`signature.json`** ← m6-lifespan-program-scores.csv (dfff50c3d1d4) + m6-statistical-summary.csv (176a8dcab93f) + m6-gates.csv (1c1260534996) + m6-species-lifespans.csv (ad28535a0957)
   `{arms: {mouse_wb: {points: [{donor, age_years, age_normalized, score, n_cells}], n: 23, slope, ci, p}, mouse_skin: {...n:3, underpowered: true}, human_skin: {...n:4, low_powered: true}, human_fibroblast: {...}, nmr: {..., tested_span_years: [1.1, 11.3]}}, stats: {diff: {...}, tost: {delta: 0.2023, p: 0.02748, established: true}, perm_p: 1.0}, coverage: {mouse: 1.0, human: 0.84892, nmr: 0.66187}, lifespans: {mouse: 4.0, nmr: 31.0, human: 122.5}}` — powers G9, G10, G11. The large post hoc TOST margin exceeds the mouse whole-body estimate and cannot exclude a mouse-sized effect.
7. **`savina.json`** ← m6-savina-reproduction.csv (cfdab4a0f3a1)
   `{curve: [{log2fc_cutoff, n_de_pval05, n_de_fdr05}], markers: {paper_rule: {cutoff: 1.0, n_pval: 29}, converted_rule: {cutoff: 1.45, n_pval: 9, n_fdr: 6}, loose_rule: {cutoff: 0.25, n_fdr: 560}}, caveats: ["Igfbp3 absent from source file", "QC 10,251 vs paper 10,232"]}` — powers G12.
8. **`power-grid.json`** <- m4-human-power.csv (03dc469a5a05)
   `{donors: [4,6,8,10,12,16,20,30,50], sds: [8,12,16,20], dispersions: [0.1,0.2,0.4], tiles: [{n_donors, sd, dispersion, min_30y_log2fc, power_status, eligible_cohort}], ts_marker: {donors: 4, sd: 8, label: "Tabula Sapiens skin: 4 donors"}, anchors: {"4_8_low": 1.474, "4_20_low": 0.59}}` : powers G13. Readout may display only grid values.
9. **`provenance.json`** ← computed at export
   `{freeze_date: "2026-10-14", export_date, sources: {per-export: {file, sha256 (full), rows}}, assets: {per-file: {bytes}}}` — feeds every provenance drawer. SHA-256 short forms in drawers = first 12 hex chars of these.

## Figure assets (for compare-toggles & G5 & fallback images)

Copied verbatim from `figures/` (no re-rendering): `m3-dpt-vs-age.png` (17e50141a309), `m6-signature-figure.png` (d37d0ef831c4), `m6-positive-control-skin.png` (9371628bf60e), `m6-savina-reproduction.png` (ddc8ccee4034).

10. **`anage.json`** <- data/raw/anage_dataset.zip (official AnAge Build 14 export, e3ddb66e32e9)
   `[{n, latin, cls, long, kg, q, sc}]` — 4,141 species with recorded maximum longevity; `sc` flags the three species included in this explorer, not a global census of public single-cell data. Powers G14. Added 2026-10-21; wording corrected 2026-10-24.
11. **`humancohorts.json`** <- m9-cohort-summary.csv + m9-regressions.csv (scripts/run_m9_pipeline.py)
   `{cohorts: [{cohort, group, n_donors, age_min/max/sd, slope_per_year, se, ci_lo, ci_hi, p, q_bh, rho, verdict}] (11 cohorts, 621 scored control donors), celltypes: [top 40 eligible rows with global q values], total_control_donors: 621, program_human_genes: 118, program_tier: P1_Clean}` — powers G15. MSSM has slope -0.000404, p=0.2403, and BH q=0.9402 across 11 cohort tests. Primary scores use per-cell library-size normalization. Updated 2026-10-31.
12. **`dlpfc.json`** <- m10-cohort-summary.csv + m10-celltype-regressions.csv (scripts/run_m10_pipeline.py)
   `{cohort: {name: "DLPFC_Lifespan_Atlas", source: "Synapse syn52396927 (AD Knowledge Portal DUA)", n_donors: 284, age_min: 0.2, age_max: 97.0, age_sd: 24.0, n_cell_types: 27, whole_cortex: {slope_per_year: +0.000186, p: 0.3494, slope_covar_adj: -0.000209, p_covar_adj: 0.4004, ...}}, program: {up_slots: 73, down_slots: 40, measured_unique: 113, coverage: 0.9576, tier: P1_Clean}, n_concordant: 9, celltypes: [27 regression rows with global q values]}` — powers G16. `concordant` means nominal p<0.05 with the same sign in both schemes, not primary global significance. The whole-cortex summary averages only observed cell types per donor; missing cell types are not zeros. **DUA rule:** regression aggregates only; per-donor scores never leave `data/` (gitignored). Updated 2026-10-31.

## Hard gates

- `scripts/check_app_data.py` asserts every anchor in the table above against the exports and exits non-zero on any miss. Run after every export; the page build consumes only checked exports.
- No export may contain a value not derivable from the listed sources. The C3 scatter is the only multi-source export: an inner join on gene symbol between the M3 continuous bridge slope (per month) and the paper's age coefficient, mirroring `scripts/run_m3_agreement.py`.
- Word/encoding budget enforcement (Section 5 text ≤700 words, ≥70% graph area) is a Step-5 checklist item, not an export concern.
