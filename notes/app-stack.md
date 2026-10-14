# Web Application Technology Stack & Architecture Decision

**Milestone:** M8
**Component:** Conserved Aging Explorer Web Application (`app/`)
**Deployment Target:** any static host (Cloudflare Pages / Netlify / Vercel class), access-controlled staging first

---

## 1. Stack Selection & Rationale

### Choice: static vanilla JS + plain JSON assets, no build step
- **No build tooling.** `app/` is three hand-written files (`index.html`, `app.js`, `style.css`)
  plus generated `data/*.json` and `figures/*.png`. Nothing to compile, bundle, or npm-install;
  the repo convention against ad-hoc package installs makes a build toolchain a liability, not an asset.
- **No runtime dependencies.** No CDN links, no frameworks, no analytics; the only network activity is
  same-origin `fetch()` of the local JSON assets. Deployment is "copy the directory".
- **JSON assets are accepted for small curated tables.** Every asset is produced by
  `notebooks/07-app-assets.ipynb` directly from `results/*.csv`, with `app/data/manifest.json`
  recording source path + SHA-256 per asset and `scripts/check_app_assets.py` failing on any
  staleness. Total `app/data` payload is ~0.27 MB — the largest single asset is the 5,000-point
  UMAP subsample at ~132 KB.
- **Why not Parquet/Arrow + duckdb-wasm (the original M8 sketch)?** The pre-registered hosting
  budgets (≤25 MB per view, ≤300 MB total) were sized for a future with per-cell expression
  matrices. The data actually shipped at this milestone is ~270 KB across all views. Parsing
  ~270 KB of JSON costs low single-digit milliseconds in the browser; columnar formats and an
  in-browser query engine would add a dependency, a wasm payload, and generator complexity for
  zero user-visible benefit. **Honest tradeoff:** this decision is revisited the moment a view
  needs >10 MB of tabular payload — at that scale JSON parse time and memory become real and
  Parquet/Arrow is the right tool. The manifest + budget-audit machinery is format-agnostic, so
  swapping the asset encoder later does not change the verification story.
- **Serving note.** `fetch()` is blocked under `file://` in most browsers; the app shows an
  explicit banner telling the user to serve `app/` over HTTP (e.g. `python3 -m http.server`).

## 2. What ships in each view (priority order)

1. **Signature view (lead):** the longevity-outlier contrast — `figures/m6-signature-figure.png`
   prominent, the full M6 statistics table (slopes, CIs, honest p-values including the *failed*
   TOST), M6 gates (positive control, per-species program-gene coverage), AnAge denominators,
   an interactive per-sample scatter with the CSV slopes overlaid, and the Savina reproduction
   curve + control genes.
2. **Gene-level browser:** search any gene from `m5-gene-level-programs.csv`; shows per-unit
   log2FC/FDR, program membership (up/down/conserved), and the unit's catalog stats including
   `frac_ribo_mito`.
3. **Mouse aging-signal ranking:** SVG bar chart of the 17 powered units (of 163 evaluated) from
   `m2-aging-signal-ranked.csv` with bootstrap CIs, plus the calibration trio and sensitivity runs.
4. **Conservation / prior work:** M3 verdicts (13 agree / 4 partial / 0 disagree), the M5
   divergence table, the M4 human power grid, and the human-arm-untestable statement as a
   first-class result.
5. **Atlas UMAP:** seeded (seed 42) per-tissue proportional 5,000-cell subsample of
   `obsm['X_umap']`, colorable by age group / tissue / sex / cell type, with the sex×age
   confound table from `m1-sex-age-table.csv` rendered beside it.

Every view carries a provenance footnote naming its source CSV(s), and a persistent **Caveats**
button (header, all views) opens the limitations panel: human arm untestable, NMR 31.4% coverage,
TOST not established, group-mean NMR ages, human n=4, Igfbp3 absent from the NMR GEO matrix,
sex×age confounding, underpowered sub-cohorts, cross-sectional/survivorship limits, correlation≠causation.

## 3. Verification

- Generator: `notebooks/07-app-assets.ipynb` (kernelspec `agescape`), executed via the gate
  nbconvert command. It wipes stale assets, writes the JSON + manifest, asserts the golden
  numbers, and ends with an outputs-written cell (path + SHA-256 + size per file).
- Checker: `python scripts/check_app_assets.py` — manifest exists; every asset ≤1 MB; every
  source file exists with matching SHA-256; `app/data` total ≤5 MB; unlisted files under
  `app/data` are flagged. Exit 0/1.
