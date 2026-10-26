# Agescape Takeover — Final Implementation Plan v3

Companion files:
- [crucible_review_10_rounds.md](crucible_review_10_rounds.md): the 10-round review log.
- [module_specs_v2_base.md](module_specs_v2_base.md): the base module specs (Modules 1–6).

**Precedence rule:** where the v3 amendments in §3 conflict with the base specs, the amendments win.

## 1. Verdict on the three pasted reviews
- **Review 1 (external audit): correct, and the code is worse than it said.** Our audit confirmed all 8 flaws and found deeper causes underneath them (§2).
- **Review 2 (generic hallmarks/knowledge-graph prompt): reject.** It is an unrelated prompt. The only part kept is a small offline overlap table against published gene sets (Zhang 2021 GAGs, SenMayo mouse panel). Tyshkovskiy 2026 and Karakaslar 2023 are cited as "not incorporated".
- **Review 3 (closeout packet): fully integrated.**
  - Items A1–A4 are implemented.
  - Items B1–B12 are disclosed without spin.
  - Item C (verified fixes) is protected by regression gates.

## 2. Issues found beyond the pasted reviews (all verified on disk)
1. **Active DUA leak:** the tracked [m10 whole-score figure](../../figures/m10-dlpfc-whole-score-vs-age.png) plots all 284 controlled-access donors. It includes 9 donors aged 90+ and a 0.2-year infant.
2. **Fabricated NMR ortholog table:** it is an uppercased copy of the human table. Igfbp3 is present in GSE152886 as `ENSHGLG00000018612`, so "absent from GEO" is false.
   - Correct framing: global Igfbp3 shows ratio 1.019 (Wilcoxon p = 0.325), and animal-level pseudobulk CPM is 221.8 vs 223.8. Savina's signal is a batch-integrated BSC2-subcluster, cell-level result.
3. **M2:**
   - 3m male Lung pools are `immune-endo-depleted` and 2-mouse pooled.
   - `independent_filter=True` flips the #1 unit: under uniform BH, Kidney PCT overtakes Lung adventitial.
   - 1-cell pseudobulks exist.
   - Surviving units: 17 → 12 (≥10 cells/mouse) → 9 (≥20 cells/mouse); Lung adventitial is eliminated at ≥10.
   - Valid joint permutations: 15 group-size-preserving, 6 strict-margin.
4. **M3:** the code omits `sex` although the prose claims it adjusts for it. The Marrow DPT young-mask is dead code. The library-size denominator is computed after gene filtering.
5. **M5:** naive Cochran's Q is deflated because all 17 units come from 9 mice. Median Q/df is 0.59.
6. **M6:**
   - The mouse score is subcomposition-normalized (`Xsel.sum`).
   - The program contains Ins1/Ins2/Rn45s soup and bidirectional genes.
   - TOST bugs: the CI is slope ± δ, the p-value uses z instead of Welch t, and the permutation p is double-counted.
   - TOST has zero power at the proposed margins: δ_min = t(0.95,3)·SE = 0.202.
7. **M9:** no CP10K normalization, an unpinned Census version, and pediatric donors included.
8. **M10:** Brain_bank × Age × PMI confounding, and the 7/7 headline hits do not survive adjustment.
9. **Repo and reproducibility:**
   - `.gitignore` untracks 5 of 7 manifests.
   - `run_m1` raises `NameError` on a clean clone.
   - `finish_m1.py` wipes `notes/methods.md`.
   - `notebooks/07` deletes `app/data/*.json`.
   - The environment lock is missing scikit-image.
10. **Web app:** 22+3 inline styles collapse ECharts under CSP, plus 5 innerHTML sinks and 12 tooltip formatter HTML sinks. echarts has 2 `cssText` sites. Google Fonts are loaded without SRI.

## 3. Binding v3 amendments (override the base specs)

### 3.1 Execution DAG (replaces the "parallel Wave 1")
Workers run in these waves, in order:
1. **Wave 1a — Worker A** (sequential).
2. **Wave 1b — Worker B**, after 1a.
3. **Wave 1c — Worker C**, after 1b.
4. **Wave 2 — D1, D2 and D3 in parallel.** They write disjoint `m6-*`, `m9-*` and `m10-*` files.
5. **Orchestrator gate:** the parent writes `ANCHOR_LOCK.json`.
6. **Wave 3a — Worker F.**
7. **Wave 3b — E1, E2 and E3, sequential.** E3 runs the full CI sequence over the finished repo.
8. **Wave 4 — parent.** Final verification and the walkthrough crucible.

Worker scopes:
- **A:** `.gitignore`, `git rm` of the leaked figure, manifests, `data/reference/` and its manifest, `run_m1` and `set_seeds` fixes, `run_m4_power`, env files, `CITATION.cff`, path scrubs, `build_reference_orthologs.py`.
- **B:** M2/M3 fixes, `run_m2_audit.py`.
- **C:** M5 fixes, `run_m5_heterogeneity.py`.
- **D1 / D2 / D3:** M6 / M9 / M10.
- **F:** reports, notebooks 01–06, README/OVERVIEW/HANDOFF.
- **E1:** exporter, contract, `check_app_data`, `run_all`.
- **E2:** CSP, fonts, vendor patch, `app.js`, `index.html`, sync script.
- **E3:** `check_app_assets`, tests, CI.

### 3.2 Infrastructure
- **Never mutate the running `agescape` env** to regenerate the lock; build in a scratch prefix.
- **Reference data:**
  - Build it from the Ensembl 106 FTP: the mouse GRCm39 GTF plus the Compara 106 `protein_default` and `ncrna_default` mouse homologies.
  - The BioMart archive timed out live, so do not use it.
  - Keep the reference manifest at `data/reference/MANIFEST.json`, not in `data/manifest/`.
  - `fetch_data.py` and `verify_data.py` stay 1:1 with `data/raw/`; `verify_data.py` also validates `data/reference/MANIFEST.json`.
  - Commit `zhang2021-330-gags.csv`, `zhang2021-59-droplet-gags.csv` and `saul2022-senmayo-mouse.csv` (117 genes), each with its SHA-256 and source URL (Zhang supplement 3; Saul MOESM4 or MSigDB MM16098).
  - `run_all` verifies hashes by default and re-fetches only with `--refresh-references`.
- **Chromosome annotation:** `mouse-gene-chromosomes.csv` carries `chromosome_name` and `gene_biotype`. The sex-gene source is the empirical adult sex DE (`~ age + sex`, padj < 0.05, |log2FC| ≥ 1, all adult samples), plus X/Y genes. The Ransick table was unobtainable (HTTP 403), so it is cited only.
- **Wave labels:** read wave from `adata.obs["cell"]` (MACA vs 10X_P prefix), never from `obs_names`.
- **Thread variables:** set `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `NUMBA_NUM_THREADS` and `PYTHONHASHSEED` at the top of every script and notebook cell 1, before `import numpy`.
- **Census environment:**
  - M9 and M10 resolve their interpreter from `AGESCAPE_CENSUS_PYTHON`.
  - Freeze `requirements-census.lock`.
  - Tracked files contain no `/tmp` or `/home/reza`.
  - Scrub `KICKOFF.md`, `PLAN.md`, `docs/build-summary-v2.md`, `notes/gates/M0-report.md`, both manifests and both check scripts.
- **M1 rerun:** write via `.tmp` plus `os.replace`; emit `m1-animal-tissue-design.csv` without rewriting the cached h5ad needlessly.
- **Notebooks:**
  - Repoint notebooks 01, 04 and 05 to `data/reference`.
  - Fix the unsorted glob in 02.
  - Print relative paths only.
  - Make notebook 07 a thin read-only verifier. It must not call `p.unlink` and must not need `m1-processed.h5ad`.

### 3.3 Statistics
- **M2 design:**
  - Use `proto_tissue` covariate levels (Lung::depleted, Liver npc/hep, Pancreas endo/exo, bare tissue); drop strata with fewer than 2 ages.
  - Assert `matrix_rank == p`.
  - Emit gene-level `m2-s1-uniform-bh-genes.csv`.
  - Floor survival is locked: 17, 12, 9.
  - Global linked-animal permutation: 9 linked groups, M = 15 (and 6 strict); p_min = 1/15 and 1/6.
  - Report exact p = cnt / M.
- **M2 composition:**
  - Use `var_weights=total_cells` with `fit(scale="X2", use_t=True)`.
  - Apply Jeffreys shrinkage when a group is all-zero.
  - Note that undepleted Lung is female-identified.
- **M3:**
  - Fit `log1p(CP10K) ~ age_months + sex` (scale 1e4, not CPM).
  - Compute library size before gene filtering.
  - Fix or state the Marrow DPT young-mask honestly.
  - Replace the hardcoded "13 of 17" with a value read from `ANCHOR_LOCK`.
- **M5 Q\*:**
  - Estimate the correlation matrix from the Spearman → `2·sin(πρ/6)` mapping with no truncation, then eigen-floor it.
  - Report two versions: raw-SE `Q*_corr`, and empirical-null scaled `Q*_empnull` with s_u capped at 1.
  - Decompose organs by nested GLS: `Q_between = Q_total − Q_within_joint ≥ 0`.
  - Report marginal Lung/Kidney Q as diagnostics only, never subtracted.
  - Overlap tests use the 330 FACS GAGs and 59 droplet GAGs, plus SenMayo mouse (117).
- **Program flags:**
  - `^Gm\d+$` (digits to end of string only, preserving `Gm2a`).
  - `-ps\d*$`.
  - `^[0-9A-Z]+Rik$`.
  - Ensembl `gene_biotype != protein_coding`.
  - Soup flags are applied per (unit, gene) before collapsing. Examples: Mgp/Col3a1 in depleted-Lung immune units, and Ins1/Ins2.
- **M6 Spec (c):**
  - Estimate with REML point estimates plus a hand-written numpy CR2 + Satterthwaite df (about 25 lines; statsmodels 0.15 has neither). Do not install any new package.
  - Clusters are physical animals, with `3-M-7/8`, `3-M-8`, `3-M-8/9`, `3-M-9` merged.
  - Add (d1) 1m/21m/24m/30m, (d2) 21m–30m, and (d3) 14-organ holdout over adult ages.
  - Read NMR by stable ENSHGLG IDs and keep the symbol in `var["symbol"]` (never an empty string).
  - Report the TOST minimum margin δ_min, the true 90% CI, Welch-t, and exact permutation p.
  - Add a 5-animal pseudobulk NMR check.
- **M9:**
  - Request `raw_sum` and `cell_type_ontology_term_id` explicitly in `obs_column_names`.
  - Pin `census_version="2025-11-08"`.
  - Use `data/reference/m9-celltype-mapping.csv` keyed by CL ids (the 17-unit table in the Round 4 report; `Kidney::kidney cortex artery cell` is recorded as unmatched).
  - Require ≥5 genes for unit scoring.
  - Apply BH only on valid p-values.
  - Flag cohorts with age SD < 12 y as LOW-POWERED.
  - Adult primary (≥20 y), with all-ages secondary.
- **M10:**
  - Use CPM as primary.
  - The `Developmental` stratum fits without `Brain_bank` (it is 100% HBCC).
  - Assert `keep.sum() ≥ 10`, `dof ≥ 5` and full rank, with `.ravel()` on all colData vectors.
  - Output `m10-confounder-audit.csv` wide (27 rows).
  - No per-donor plots. Compute donor scores in RAM only. Delete `data/interim/m10-donor-whole-scores.csv`.

### 3.4 Web app
- **echarts and CSP:**
  - Patch the 2 `cssText` sites in `echarts.min.js` to CSSOM property writes, and document the SHA-256 in `app/vendor/README.md`.
  - Patch richText only when `option.tooltip` already exists, with `confine: true` and ≤5 lines.
  - Use UTF-8 characters only. Ban HTML entities as well as tags.
  - Add `app/favicon.svg` (no `data:` URI).
  - Use `classList.toggle("modal-open")` in place of `document.body.style`.
- **Provenance:**
  - Remove `ranking.json`, which leaves 12 files.
  - Fetch `provenance.json` in `init()` and hydrate the drawers with DOM nodes.
  - The manifest whitelist excludes free-text `notes`.
  - The provenance timestamp comes from git or `SOURCE_DATE_EPOCH`.
- **UI contract:** every chart shows the remediated primary as the main mark and the original pre-registered result as a muted baseline.
  - G1: S1 bars with S0 diamonds.
  - G9: plot partial residuals so the points match the Spec (c) slope.
  - G10: use fixed and calibrated TOST bands.
  - G12: remove the false "Igfbp3 absent" string.
  - G15: stays at the bottom of ch4, retitled.
  - G16: a paired unadjusted-vs-adjusted forest plot (7 → 0 under BH).
  - Cap copy at 35 words per support block, using `<details>` for caveats.
  - Retitle around the audit framing, subject to your decision Q2.
- **Tests:**
  - `check_app_assets.py` becomes a real asset and CSP verifier with a vendor size budget.
  - DUA tests are file-specific (row counts, no donor columns or ID tokens, no `0.2` age_min in M10 files, self-excluded patterns, and no hardcoded controlled IDs).
  - Add a local link checker.
  - CI pins `actions/checkout@v4` and `actions/setup-python@v5`.
  - CI order: `py_compile` → check scripts on the untouched checkout → unit tests → export plus `git diff --exit-code`.
- **Headless verification:**
  - Use the cached Chromium at `~/.cache/ms-playwright/chromium-1248/chrome-linux64/chrome`, driven over stdlib CDP.
  - The `chrome-headless-shell` build is missing 9 shared libraries, so do not use it.
  - Serve with `python -m http.server --bind 127.0.0.1`.
  - Assert 0 CSP violations, 0 console errors and canvas height > 0 on all 20 charts.

### 3.5 ANCHOR_LOCK.json keys (parent-generated, read-only for workers)
`n_raw_cells`, `n_post_qc_cells`, `n_mouse_pools`, `n_physical_mice`, `top_unit_s1`, `top_unit_s1_nde`, `lung_adv_s0_nde`, `lung_adv_s1_nde`, `n_floor10_units`, `n_floor20_units`, `n_p0_entries`, `n_p0_unique`, `n_p1_clean_unique`, `q_star_empnull_median`, `m3_agree_count`, `m6_mouse_specc_slope`, `m6_nmr_slope`, `m6_delta_min`, `nmr_igfbp3_umis`, `m9_n_cohorts`, `m9_n_fdr_sig`, `m10_n_unadj_nom`, `m10_n_adj_fdr_sig`. Workers must not invent numbers or revert to old values.

## 4. Security verification plan
The app is a static client-side viewer with no backend, auth, database or uploads. Those skill sections are N/A (`TODO(security)`: they become mandatory if an API is ever added).

- **Applied controls:**
  - Zero innerHTML, tag-free tooltips, and CSV-injection neutralization that preserves negative numbers.
  - Strict CSP with local fonts.
  - The DUA purge, with a `git filter-repo` option that needs your sign-off.
  - Path-scoped Cloudflare `_headers` snippet for `/agescape*` only. It never touches `/*`, which would break `rezaalipour.com`.
- **Dependencies:** no new packages. Existing locked ones need no new scan.
- **Verification:** run the security scanner skill on all new files, then produce the audit report in the walkthrough.

## 4b. Decisions made and added module
- **Scope:** Track 1 + Track 2 (full rebuild).
- **Title:** keep "Conserved Aging Explorer". Add the subtitle "within mouse; cross-species conservation not established". Do not use "audit/transport" wording anywhere.
- **History rewrite:** not done. Delete the leaked figure going forward and flag the old commits in `HANDOFF.md`.
- **Terminology:** use "aging signature" everywhere in text, not "program". "Program" survives only in code identifiers (`P0..P3`, `m6-program-genes.csv`), with the mapping noted in `HANDOFF.md`. Add a lint test that fails if "program" appears in the visible app text.
- **Execution:** the parent agent runs the waves sequentially in this session.

**Module 4b — pre-specified cross-species concordance tests** (exploratory, locked before any result is computed; reported whether null or positive):
1. **Gene-level concordance.** For each mouse unit, compare per-gene mouse log2FC with the human per-gene age slope in the mapped homologous cell type (`data/reference/m9-celltype-mapping.csv`), on 1:1 orthologs, using `P1_Clean` and all clean DE genes.
   - Metrics: sign agreement and Spearman rho.
   - Null: a gene-label permutation, plus a sex-trimmed sensitivity run.
   - Multiple testing: BH across units × cohorts.
   - Output: `results/m9-gene-level-concordance.csv`.
2. **Published-gene overlap.** A Fisher test with sign concordance of `P1_Clean` against the 330 FACS GAGs and 59 droplet GAGs (Zhang 2021) and the SenMayo mouse panel (117 genes). The same test is run against the human cohort slopes.
3. **Harmonized 3-species gene set.** Score mouse, human and naked mole-rat on the same shared 1:1 ortholog genes. Use the Hedges-corrected standardized slope and do not compare raw magnitudes.
4. **Reporting rule.** Any positive result is labelled exploratory and is never moved into the primary hero claim.

## 5. Open decisions for you (resolved)
1. **Scope:** Track 1 + Track 2 (full rebuild; recommended) or Track 1 only (hygiene, DUA purge, Review 3 fixes A1–A4)?
2. **Rename:** retitle "Conserved Aging Explorer" to an audit/transport framing (recommended)? If so, which wording?
3. **History rewrite:** run `git filter-repo` to purge the leaked figure from public history (force-push)? If not, it stays in old commits.
4. **Execution:** the Single-Cell Biology subagent reported quota exhaustion (429, resets in about 167h). If `flash` workers hit the same limit, I execute the waves myself in this session, sequentially. The plan is identical either way.

After execution I write `walkthrough.md` and run it through the same 10-round crucible.
