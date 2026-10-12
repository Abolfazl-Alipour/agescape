# Conserved Aging Explorer — Master Plan

**Project:** Cross-species single-cell aging program (mouse ↔ human ↔ longevity outliers) + public web app.
**Owner:** Reza. **Executor:** coding agent, operating from `HANDOFF.md`.
**Date:** 2026-10-11. **Status:** v2.0 — reviewed (10 passes, see `notes/reviews.md`).

---

## 0. The one-paragraph version

Reproduce and extend a rigorous cross-species single-cell aging analysis on public data
(Tabula Muris Senis mouse → SenNet/Tabula Sapiens human → naked mole-rat outlier contrast),
with every number traced to an executed notebook cell, every methodological choice
pre-registered and justified, divergence reported as a first-class result, and a static
precomputed web app (Conserved Aging Explorer) as the public centerpiece. Phases are
gated by Reza's sign-off. Nothing public without his explicit go-ahead.

---

## 1. Standing rules (non-negotiable, inherited from Reza)

1. **Never fabricate.** Every claim traces to a cell in a notebook that runs. No invented numbers.
2. **Cite only what was actually read.** A paper is cited only after the executor has opened its
   full text and recorded the exact citation + accession numbers in `notes/methods.md`.
   The bibliography in §10 lists papers *to read*; bibliographic details are verified at read time.
3. **Eval discipline.** Every quantitative claim carries a robustness check (across tissues, age
   definitions, parameters, species). Show the check.
4. **Repo hygiene.** No credentials ever — not in code, notebooks, committed outputs, or shell
   history. `.gitignore` from day one. No data files committed. Any required token lives in a
   gitignored `.env`, read via `os.environ`, never printed.
5. **Gated phases.** Reza signs off each milestone before the next starts (protocol in `HANDOFF.md`).
6. **Nothing public without explicit go-ahead.** Repo, report, app, staging — everything stops
   at Reza's check. Staging links are access-controlled (see M8).

---

## 2. Verified local data inventory (checked 2026-10-11, this machine)

| File | Size | SHA-256 | Status | Content (verified by reading HDF5 structure) |
|---|---|---|---|---|
| `/home/reza/projects/aging_hypothalamus/tabula-muris-senis-droplet-official-raw-obj.h5ad` | 3.8 GB | `86e45e88f9c475253792f1bbcdf5851a5be0ef35da72fced607ed261aebc099d` | **Intact** | Tabula Muris Senis (TMS) **droplet** raw counts: 245,389 cells × 20,138 genes, 16 tissues, ~124 cell types, 20 mice, 6 age groups. No PCA/UMAP precomputed — the pipeline computes them. |
| `/home/reza/projects/aging_hypothalamus/data/tms_brain.h5ad` | 304 MB | `31ee7ffe9fd11fd128925c2dae862399f5c54f27aef9c026dbeaa077ea143339` | **Intact** | TMS **FACS** `Brain_Non-Myeloid` subset: 7,249 cells × 22,966 genes, ages 3m/18m/24m, 14 mice, 17 cell types. |
| `/home/reza/projects/dry_lab/data/brain_aging.h5ad` | 304 MB | `31ee7ffe…ea143339` | **Intact** | Byte-identical duplicate of the file above. Backup only. |
| `/home/reza/projects/aging_hypothalamus/data/23872787{,.1}` | 0 B | — | **Dead** | Empty failed downloads. Executor deletes them in M0 (they are outside this repo; record the deletion in the M0 report). |

**Machine:** 49 GB RAM, 12 cores, 247 GB free disk, conda at `/home/reza/miniconda3`, internet
reachable (cellxgene.cziscience.com confirmed 200 OK). OS: record exact distro + glibc in M0 manifest.

### 2.1 Verified facts that shape the analysis (read from file metadata, not the paper)

- **Age groups (droplet):** 1m (25,980 cells), 3m (45,602), 18m (44,645), 21m (35,828),
  24m (37,660), 30m (55,674).
- **Sex balance per age group (verified):** 1m = male-only; 3m = mixed; 18m = mixed;
  21m = female-only; 24m = male-only; 30m = male-only.
  → **Age and sex are partially confounded.** Primary contrast must be sex-balanced (§4.2);
  every pooled analysis prints the sex×age table next to it.
- **The droplet object has no brain** (TMS brain was FACS-only). Local `tms_brain.h5ad` is
  `Brain_Non-Myeloid` only. Full-FACS-object download is an optional M3 extension (§8.3).
- **Known data-cleaning bug:** `subtissue` values contain trailing spaces (`'Cerebellum '` vs
  `'Cerebellum'`). All string obs columns are stripped in M1 with before/after value counts logged.

---

## 3. Repository structure and global conventions (from day one)

```
agescape/
├── PLAN.md  HANDOFF.md  README.md  environment.yml  environment.lock
├── .gitignore
├── scripts/            # fetch_data.py, verify_data.py, check_app_assets.py, run_all.py
├── data/
│   ├── raw/            # immutable inputs (gitignored); .partial/ for in-flight downloads
│   ├── interim/        # pipeline checkpoints (gitignored; prune aggressively)
│   └── manifest/       # one JSON per file: url, ids, version, size, sha256, dates
├── notebooks/          # 01-…07-, executed top-to-bottom, committed WITH outputs
├── src/                # importable functions only; each has a smoke check in Notebook 01
├── results/            # golden-output CSVs + SHA-256s (committed; small)
├── figures/
├── notes/              # qc-decisions.md, age-bins.md, methods.md, celltype-mapping.md, …
│   └── gates/          # M<n>-report.md gate reports to Reza
├── report/
└── app/                # static site; app/data/ precomputed assets + manifest.json
```

**Conventions (binding):**
- Notebooks never import other notebooks; reusable logic lives in `src/`.
- Notebooks carry kernelspec `agescape`; gates execute via
  `jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.kernel_name=agescape --ExecutePreprocessor.timeout=600`.
- **Determinism policy:** `src.set_seeds()` (numpy + python `random`, seed=42) called at the top of
  every notebook; every stochastic scanpy call takes an explicit `random_state`; Leiden
  `flavor`/`n_iterations` recorded in `notes/methods.md`; all subsamples/bootstraps documented seeds;
  UMAP run single-threaded (`n_jobs=1`) for reproducibility.
- **Outputs cell:** every notebook ends with a cell listing every file it wrote, each file's
  SHA-256, the git commit hash, and the data-manifest digest.
- **Golden outputs:** key tables (DE results, rankings, conservation scores) are written to
  `results/*.csv`; reference copies + SHA-256s are committed. Fresh-clone reproduction must
  reproduce them (exact match where deterministic, documented tolerance elsewhere).
- **log2FC orientation:** every fold-change in the project is **old vs young** (positive = up with
  age). Notebook 02 asserts this convention with a known-direction gene check.
- Notebooks are committed fully executed. `figures/` is written only by notebooks.

---

## 4. Pre-registered methodological decisions

### 4.1 Environment

Conda env `agescape`, Python 3.11, channels `[conda-forge, bioconda]`, `channel_priority: strict`.
**All** dependencies are installed at M0 (including "optional" ones — no mid-project ad-hoc installs):
`scanpy>=1.10, anndata>=0.10, leidenalg, python-igraph, scikit-learn, matplotlib, seaborn, jupyter,
nbconvert, ipykernel, requests, h5py, pyarrow, openpyxl, scikit-misc, pydeseq2, scrublet` (scrublet =
python doublet detection, run per sample; pip-declared in the lock if conda-forge lacks it). Executors record the resolved
lock with `conda list --explicit > environment.lock` at M0; the README reproduces from `environment.lock`.
OS (distro, glibc) recorded in the manifest. Rationale for py3.11: scanpy ecosystem maturity; document if changed.

### 4.2 Age and sex handling (written to `notes/age-bins.md` before any DE run)

- **Mouse primary contrast: 3m vs 18m** — the only pair with both sexes in both groups.
  The notebook prints the sex×age cell-count table and tests 3m-vs-18m sex ratios; DE lists are
  checked for sex-chromosome enrichment.
- **Sensitivity runs (all reported, none pooled silently):** (a) extended bins young=1m+3m,
  old=18m+21m+24m+30m (hypothesis: wider gap strengthens signal — stated before computing);
  (b) 18m+24m male-only run; (c) continuous-age model (§4.3 bridge).
- 30m is a **survival-selected cohort** — the Limitations section must say so; extended-bin
  strengthening cannot be attributed to age alone.
- **Human: age is a continuous covariate, never dichotomized.** Inclusion thresholds per
  cell type: ≥6 donors and donor-age SD ≥12 years; per-type power note (min detectable |log2FC|
  at FDR 0.1 given n donors) printed in Notebook 04 and stored in `results/`.

### 4.3 Differential expression design

- **Unit: tissue × cell type** (never pooled across tissues). Authority labels: TMS provided
  `cell_ontology_class` (use `free_annotation` only as a tie-break, logged).
- **Primary: pseudobulk per (mouse × tissue × cell type)**, sum counts, `pydeseq2`,
  design `~ age_group + sex` (primary contrast 3m vs 18m); independent filtering + BH-FDR;
  report log2FC + CI. Eligibility: ≥3 mice per group and ≥50 cells per group; ineligible
  tissue×type combos are listed, never silently dropped.
- **Secondary: cell-level Wilcoxon** (`sc.tl.rank_genes_groups`, FDR, log2FC, pct-expression)
  on the same units with the same eligibility floors, for comparability with atlas-style analyses.
- **Diagnostics on the top-ranked types:** age×sex interaction test; per-mouse effect-size
  consistency (fraction of DE genes same-direction in each mouse); leave-one-mouse-out refit;
  cells-per-mouse × age-group table.
- **Calibration trio (M2 gate artifacts):** (a) label-permutation null — shuffle age labels,
  expect ≈0 DE genes at the same thresholds (validates metric + FDR path); (b) positive controls —
  known TMS aging markers (e.g., Cdkn2a, Igfbp3, SASP genes) must appear with correct direction;
  (c) sex-effect floor — DE by sex within 3m only (sex is balanced there), reported as the
  noise floor for sex-linked artifacts.
- **Continuous-age bridge model** (for the M3 agreement check): pseudobulk `~ age_months + sex`.

### 4.4 Aging-signal ranking (defined before computing)

Per tissue×cell type, from the **pseudobulk primary results**: `Aging Signal = n DE genes at
BH-FDR < 0.05 with |log2FC| ≥ 0.25`, reported with **median |log2FC|** of significant genes and
**bootstrap-over-mice 95% CI on the ranking**. The full ranked table (all eligible types) goes to
`results/` and the report supplement — no top-10-only cherry-pick. Robustness: Spearman rank
correlation between pseudobulk and cell-level rankings across all eligible types; scatter of
metric vs n_mice and vs n_cells (power check); recomputation under extended bins. The claim
"cell type X ages fastest" is always accompanied by its CI and never made from a type with <3 mice
per group.

### 4.5 Conservation analysis, mouse → human (defined before computing)

1. **Orthologs:** one-to-one mouse↔human via Ensembl BioMart; record Ensembl release + query date
   in `notes/methods.md`; cache in `data/interim/`. Same pipeline later maps mouse→NMR (M6).
2. **Cell-type mapping:** a written `notes/celltype-mapping.md` table (Cell Ontology ID-based,
   curated, with rationale; one-to-many collapse rules; marker-based spot checks). Unmatched types
   are listed with disposition **excluded, reported** — never scored as zero, never dropped silently.
   Conservation is computed on **tissue-matched** pairs only.
3. **Human DE:** pseudobulk per (donor × cell type), `~ age_years + sex`. Tabula Sapiens and SenNet
   are **separate cohorts, never merged** (cohort×age confound); SenNet is the primary human aging
   contrast (chosen by the M4 selection rule), TS v1 (10x-only subset) is the replication cohort.
4. **Two axes, reported separately, never pooled:**
   - *Direction-conservation (primary):* Spearman correlation of effect sizes (mouse log2FC vs
     human log2FC) over all expressed 1:1 orthologs in the cell type, with per-mouse block bootstrap CI.
   - *Replication fraction (secondary):* fraction of mouse DE genes (FDR<0.05, |log2FC|≥0.25)
     whose human ortholog is same-direction significant at the **same FDR** (no asymmetric
     thresholds). A lenient same-direction-only fraction is reported as a separate column.
5. **Null calibration:** permutation null — shuffle human donor ages, recompute both axes,
   report observed-minus-null. Replication fractions are also reported **with and without
   RPL/RPS/MT- genes** (housekeeping inflation guard).
6. **Power stratification:** human cell types below the §4.2 inclusion thresholds are flagged
   `low-powered` and shown in a separate table section; their numbers are never read as
   "divergence."
7. **Divergence is a first-class result**, interpreted as: cross-species-identical programs are
   the better drug-target candidates. The divergence table (program | cell type | conserved /
   divergent | evidence | power flag) is a headline artifact.

### 4.6 Signature figure definition (the M6 falsifiable test)

- A **conserved program** = the set of mouse DE genes (§4.3) whose 1:1 human orthologs pass the
  strict replication rule (§4.5.4) in a matched cell type. Programs are defined on mouse + human
  **before** any NMR look (no circular "conserved because human rises").
- Per program: per-cell module score (`sc.tl.score_genes`, method fixed in `notes/methods.md`),
  regressed on **life-stage-normalized age** (fraction of species lifespan; conversion table from
  AnAge + documented sources, cited only after reading). Headline statistic: **species × age
  interaction** on program score, plus per-species slope + 95% CI.
- **NMR flatness = equivalence, not absence of significance:** two one-sided tests (TOST) with a
  pre-registered margin δ on the NMR slope; reported as slope + CI.
- **Positive control gate:** the same scoring pipeline must recover the known mouse-skin age
  effect before any NMR claim is made. **Coverage gate:** program-gene detection rate in NMR and
  the NMR age distribution are printed next to every NMR claim.

### 4.7 Pseudotime (M3)

Diffusion pseudotime (`sc.tl.diffmap` + `sc.tl.dpt`). Tissue chosen by written rule (tissue with the
strongest graded differentiation pattern in the M1 UMAP; choice recorded in the notebook). Root =
cluster with highest expression of a marker set named in the notebook; root-choice sensitivity run.
**Mandatory guards:** report corr(pseudotime, age); build the graph on young cells and project old
cells; check graph connectivity per age group. RNA velocity stays out of scope.

### 4.8 QC (thresholds proposed in Notebook 01, binding before DE)

Per tissue: min/max genes, mt% ceiling, doublet scoring/removal per sample (`scrublet`);
**age-stratified** n_genes/UMI/mt% distribution plots (old tissue dissociates worse — a fixed cutoff
is itself age-biased; thresholds must not preferentially drop old cells). All values + rationale in
`notes/qc-decisions.md` *before* Notebook 02 runs. String obs columns stripped (§2.1).

---

## 5. Milestones (binary acceptance + Reza gate; operational detail in `HANDOFF.md`)

- **M0 — Scaffold & verification.** Skeleton, `.gitignore`, env from §4.1, lock exported, kernelspec
  registered; copy + verify local data into `data/raw/`; `scripts/fetch_data.py` + `verify_data.py`
  working; manifest JSONs written (including CELLxGENE collection/dataset IDs for the local TMS file);
  delete the two dead 0 B files; OS recorded; golden-output convention in place.
  *DoD: verify_data.py exits 0 on all raw files; `python -c "import scanpy, src"` in env passes;
  environment.lock committed.*
- **M1 — Pipeline (TMS droplet).** Notebook 01: per-tissue QC (§4.8) → concat → normalize(log1p) →
  HVG (per tissue, union; flavor `seurat`/`cell_ranger` — **never `seurat_v3`**, never `pp.scale`
  outside the HVG subset) → PCA → neighbors → UMAP → Leiden → majority-vote annotation vs provided
  labels + confusion matrix; cluster-composition table vs paper; atlas UMAP on a per-tissue
  proportional subsample (≥500-cell floor, seed 42); checkpoints per stage.
  *DoD: 01 executes clean via nbconvert gate command; listed figures exist; composition table +
  confusion matrix in results/; qc-decisions.md + age-bins.md written.*
- **M2 — Mouse DE.** Notebook 02: §4.3 design on all tissue×type units; §4.4 ranking; calibration
  trio; composition analysis (fraction-vs-age per tissue×type, mouse = unit); pseudotime guards (§4.7)
  begin (tissue choice). *DoD: results/de-aging-ranked.csv with CIs; calibration trio artifacts;
  permutation null ≈ 0 DE; positive controls correct direction; methods.md updated.*
- **M3 — Trajectory & agreement.** Notebook 03: dpt with guards; agreement table vs TMS paper claims
  tested against the paper's **published gene lists** (continuous-age bridge model as the methodological
  bridge — the paper used a continuous-age MAST-style model, not dichotomous Wilcoxon; exact claim
  list in `notes/methods.md` before computing). *DoD: dpt-vs-age figure + corr; agreement table with
  agree/disagree/cannot-test + stated causes; every "disagree" has a checked explanation.*
- **M4 — Human data & orthologs.** Notebook 04: fetch SenNet public scRNA-seq dataset(s) by the
  pre-registered rule (public, token-free; largest donor-age span with ≥6 donors ≥60y, 10x); fetch
  TS v1 whole-body h5ad (permanent CELLxGENE dataset URL, `curl -C -`, local SHA-256 — never ETag;
  figshare release-1.0 per-tissue zips as documented fallback); per-tissue, backed-mode processing
  only — **no whole-object in-memory ops for >250K cells**; human DE cohorts built per §4.5.3;
  BioMart 1:1 ortholog map with release stamp; `notes/celltype-mapping.md` drafted; power table per
  human cell type. *DoD: verify_data.py covers new files; human cohorts load per-tissue; ortholog
  cache + release recorded; mapping table exists; power table in results/.*
- **M5 — Conservation.** Notebook 05: §4.5 axes + null + housekeeping guard + power strata;
  mouse-vs-human effect-size figures; validation vs the published conserved-aging programs (read
  full texts first — §10 bibliography; test = gene-set overlap + sign concordance, explicit
  cannot-test rows); divergence table. *DoD: conservation + divergence tables with stats in
  results/; each published program has an explicit agree/disagree/cannot-test row with evidence.*
- **M6 — Longevity outliers.** Notebook 06: Savina NMR skin (processed counts CSV from GEO —
  accessions read from the paper; mirror the paper's sample filtering, e.g. drop sample 2V, and say
  so); §4.6 signature figure with interaction test + TOST; Hilton NMR-vs-mouse immune data as an
  **identity contrast, explicitly not an aging contrast** (young NMRs only — honest framing or drop);
  AnAge official zipped-TSV export (build + date recorded) + species table with availability flags
  (single-cell / genome / metadata-only); `notes/outlier-limits.md` (thin outlier tissue coverage;
  no multi-tissue NMR atlas; NMR "aged" ≈ middle-aged by NMR standards — stated plainly).
  *DoD: signature figure with slope+CI per species + interaction p + TOST result; Savina
  reproduction row (DE-count-vs-cutoff curve); species table with flags; limitations note.*
- **M7 — Report & polish.** README (signature figure above the fold; sex×age honesty moment on the
  first screen; 3-command repro; per-notebook runtimes on this hardware; "20 minutes? read this"
  section); notebooks 01–06 re-executed via the gate command; `report/conserved-aging-programs.md`
  (8–12 pp) with: abstract, significance statement, data/code availability, reproducibility
  statement, **deviations-from-pre-registration log**, **claims that failed robustness** section,
  Limitations (correlation≠causation; divergence; TS age range; outlier coverage; survivorship;
  harvest-time uncertainty; dropout), AI-assistance disclosure; Zenodo-archived release of the
  gated commit. **Fresh-clone test:** clone → fetch_data.py → env from lock → execute all →
  diff golden outputs. *DoD: all prior DoDs green from a fresh clone; zero TODOs; every figure
  numbered and cited in the report; deviations log exists. Reza's explicit go-ahead gates any
  public step.*
- **M8 — Web app.** Static site; stack choice + rationale in `notes/app-stack.md` (host: Netlify/
  Vercel/Cloudflare class — **not GitHub Pages**; private-repo builds; deploy via GitHub App OAuth,
  never CLI tokens; staging behind access control or unlisted URL per Reza's M7 decision).
  Asset spec: Apache Parquet/Arrow only (no JSON arrays), float16 expression, curated gene panels,
  per-view lazy loading, per-view subsample cap ≤150K cells (seeded), budget caps ≤300 MB total
  app/data and ≤25 MB per view; `notebooks/07-app-assets.ipynb` generates all assets and writes
  `app/data/manifest.json` (asset → notebook → cell → commit → source-notebook hashes);
  `scripts/check_app_assets.py` fails on staleness. Feature priority: **(1) signature
  longevity-outlier contrast view** (lead the demo), **(2) gene-level query** (mouse/human/NMR age
  trajectories + conservation status + CIs), (3) conserved-vs-divergent program panels, (4)
  tissue/cell-type browser, (5) side-by-side mouse+human UMAPs colored by **age group** (categorical
  — no per-cell "transcriptomic age" pseudo-clock; if a within-species relative score is shown it is
  labeled as such), each with a composition toggle, (6) trajectory view. Shareable per-view URLs;
  SVG/PNG export with a methods footnote; provenance tooltips (notebook + cell); Caveats page
  linked from every view (sex×age confound lives there). *DoD: staging link live + access-controlled;
  check_app_assets.py passes; Reza click-through sign-off; public launch only on explicit go-ahead.*

---

## 6. Failure playbook (summary; full version in `HANDOFF.md`)

- **OOM guardrails (49 GB):** never `pp.scale` outside the HVG subset (`zero_center=False` to stay
  sparse); never HVG flavor `seurat_v3`; keep X float32; avoid `.raw` retention (or store backed);
  per-tissue processing; expected steady peak 15–25 GB — checkpoint `.h5ad` per stage; prune
  checkpoints (keep stage boundaries only); `df -h` check at M4.
- **Checksum mismatch / truncated download:** resume with `curl -C -`; rename into place only after
  local SHA-256 matches the manifest; never interpret a `.partial` file.
- **Human data >RAM:** per-tissue slices via `read_h5ad(backed='r')`, one tissue at a time; or
  `cellxgene-census` (record census release date). Never claim a result computed on a subset
  without saying so.
- **Ensembl/BioMart down:** retry with backoff; fallback Ensembl REST; record release; never
  hand-curate orthologs silently.
- **Result disagrees with a paper:** do not force agreement; log claim/result/why; check confounds
  first (sex×age, composition, QC thresholds, continuous vs dichotomous model).
- **Gate feedback:** conditional feedback = not approved; protocol in `HANDOFF.md`.
- **Stuck >30 min on one error, or any demand for credentials, or any result that would violate
  rule 1: STOP and report to Reza** (do not route around).

## 7. Scope guards (out of scope, from Reza)

No new biological claims beyond the data. No wet-lab. No single-cell for genome-only species
(AnAge metadata only). No GenAge/DrugAge/CellAge/LongevityMap scope creep. No bots; one app.
RNA velocity excluded. Tabula Sapiens v2 excluded by default (53 GB, adds cells not old donors —
v1 + SenNet carry the human arm); adoption requires a §8 note.

## 8. Changes an executor may not make silently

Any change to a number or choice printed in §4 (age contrast, thresholds, FDR levels, metric
definitions, subsample sizes, seeds policy, cohort choices) requires: a written justification +
robustness check in `notes/` referencing the original choice, **then** proceed; if the change is
not robustly justifiable, halt and report to Reza. Same protocol for: swapping DE method; changing
primary age contrast; adopting TS v2; downloading the full TMS FACS object (allowed only as an M3
extension with batch-effect handling documented, and only when M1–M2 are green); keeping any
"transcriptomic age" score (must be labeled within-species relative, never years).

## 9. Definition of done for the whole program

An interviewer at an aging company can open the repo and the app and conclude "this person has
done the actual work": clean code, documented decisions, honest limitations (including a visible
deviations log and failed-robustness section), figures that communicate, signature result that is
falsifiable and survived its own robustness checks, and a product a scientist would use.

## 10. Reading list (read before the milestone that needs them; record exact citations in notes/methods.md at read time — bibliographic details below are leads to verify, not to cite)

- **Tabula Muris Consortium, Nature 2020** (TMS paper) — before M1. Verify its DE method and
  published aging gene lists (supplementary tables / MSigDB `TABULA_MURIS_SENIS_*_AGEING` sets)
  before M3.
- **Schaum et al., Nature 2020** (cross-species life-stage single-cell predecessor) — before M5;
  its life-stage normalization informs §4.6.
- **Tyshkovskiy / Gladyshev et al., Nature (2025–2026, verify exact issue at read)** — cross-mammalian
  conserved aging signatures; before M5. The report needs a "relationship to prior work" paragraph
  (this project's differentiation: single-cell cell-type resolution + longevity-outlier contrast).
- **AP-1 inflammaging in immune aging, human + mouse (Karakaslar et al., Aging Cell 2023 — verify at
  read)** — before M5; validation restricted to immune cell types.
- **Savina et al., Aging 2022** (NMR skin; GEO GSE152886 — verify accession in the paper text) — before M6.
- **Hilton et al., PLoS Biology 2019** (NMR immune; GEO GSE132642 — verify) — before M6; remember:
  young NMRs only → identity contrast, not aging.
- **Tabula Sapiens v1, Science 2022 (+ v2 bioRxiv 2024 note)** — before M4.
- **AnAge database** (official zipped TSV export; record build date) — before M6.
- **SenNet consortium data portal docs** (public API, token-free queries) — before M4.

---

## 11. Review log

10 passes: 5 independent reviewer agents (reproducibility; scRNA-seq methods; portfolio/hiring
bar; literal-executor ambiguity; data/web engineering) + 5 targeted passes by the planner
(consistency; failure modes; scope/time; handoff completeness; final adversarial read).
Every finding → fix is recorded in `notes/reviews.md`. Plan version history: v1.0 draft →
v2.0 all 10 passes applied.
