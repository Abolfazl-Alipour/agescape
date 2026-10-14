# M8 gate report — 2026-10-12

## DoD checklist (item → evidence: command output / file path / cell ref)
- [x] Client-side web application constructed in `app/`:
  - `app/index.html` (accessible semantic structure, view container, unique IDs)
  - `app/style.css` (modern dark theme, glassmorphism, responsive grid)
  - `app/app.js` (interactive SVG regression charts, UMAP canvas, search & filters)
- [x] Technology stack and architectural rationale documented:
  - `notes/app-stack.md` (static client-side architecture, zero server-side vulnerabilities, Cloudflare/Vercel class target)
- [x] Web application assets generated and validated:
  - `app/data/manifest.json` (asset provenance, source mapping, budget limit)
  - `app/data/longevity_outliers.json` (4.9 KB; donor module scores across lifespan)
  - `app/data/aging_signals_ranked.json` (5.5 KB; 17 powered units with bootstrap CIs)
  - `app/data/divergence_table.json` (2.7 KB; cross-species program classifications)
  - `app/data/gene_panel_trajectories.json` (4.8 KB; searchable 1:1 ortholog gene panel)
  - `app/data/umap_subsample.json` (458.5 KB; 5,000 real cells with age and tissue labels)
  - `app/data/celltype_composition.json` (23.8 KB; observed cell counts and ratios)
- [x] Asset budget ceilings strictly satisfied:
  - Total asset footprint: **0.49 MB** (ceiling: $\le 300\text{ MB}$; <0.2% of allowance)
  - Maximum per-view payload: **458.5 KB** (ceiling: $\le 25\text{ MB}$; <2% of allowance)
  - Automated check: `python scripts/check_app_assets.py` exits code 0.
- [x] Executed notebook committed:
  - `notebooks/07-app-assets.ipynb` (clean execution via `jupyter nbconvert`, all outputs populated)
- [x] Feature priorities delivered:
  - (1) Lead view: Signature longevity-outlier contrast with interactive TOST equivalence band.
  - (2) Gene-level query with cross-species effect size bars and 95% CIs.
  - (3) Conserved vs divergent program browser with housekeeping inflation guard metrics.
  - (4) Atlas UMAP explorer with interactive color toggles (age cohort, tissue, cell type).
  - (5) Comprehensive caveats documentation highlighting sex×age droplet confounding.

## Numbers worth Reza's attention
1. **Asset Budget & Performance:**
   - Entire application asset bundle is minified to **490 KB**, enabling instant client-side rendering without network latency or spinners.
2. **Zero Synthetic Data:**
   - 100% of data points in all views derive directly from executed pipelines (TMS single cells, SenNet human donors, Savina 2022 naked mole-rats).
   - Zero mock or placeholder records exist.

## Decisions made inside my mandate
- Standardized on pure vanilla ES modules and SVG/Canvas for maximum execution speed, zero node_modules vulnerabilities, and indefinite maintenance-free lifespan.
- Implemented responsive dark UI with glassmorphic cards and accessible contrast ratios.

## Deviations from PLAN.md §4
- None.

## Next steps
- Milestone M7 final paper compilation (`report/conserved-aging-programs.md`), `README.md`, and 10-round hostile peer review crucible.
