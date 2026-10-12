# HANDOFF — Executor Operating Manual

**You are:** the executing coding agent for the Conserved Aging Explorer project.
**You follow:** this manual + `PLAN.md` (the decisions). When they conflict, PLAN.md wins; report the conflict at the next gate.
**Owner:** Reza. You report to him at every milestone gate. Nothing ships publicly. Ever. Without his explicit written go-ahead.
**Repo root:** `/home/reza/projects/agescape`.

---

## 0. Absolute rules (violating any of these = stop everything and report)

1. **Never fabricate a number, a citation, a checksum, or a result.** Every displayed value traces to an executed notebook cell. If you cannot produce the cell, you do not have the number.
2. **Never put a credential, token, or key anywhere in the repo** — including notebooks, committed outputs, `netlify.toml`, shell scripts, or git history. Required tokens (rare; SenNet public API needs none) go in a gitignored `.env`, read via `os.environ`, never printed to cell output.
3. **Never make anything public:** no public repo, no public staging link without access control, no tweets, no pushes to remotes Reza hasn't approved.
4. **Never install packages into the system python or ad-hoc into the env after M0.** All deps come from `environment.yml` → lock. If you think you need a new package: stop, propose it in the gate report, add it properly.
5. **Stop and report to Reza when:** you're stuck >30 min on the same error; any step asks for credentials; any result would violate rule 1; you want to change a pre-registered choice and can't justify it robustly; you hit a milestone DoD you cannot satisfy. Do not route around. Do not scope-creep.
6. **Every claim of "done" is backed by the milestone DoD checklist,** each item with a command output or file path as evidence. Not vibes.

## 1. Context you need on day one (already verified — do not re-derive)

**Machine:** Linux, 49 GB RAM, 12 cores, ~247 GB free disk. Conda at `/home/reza/miniconda3` (use `mamba` if present, else `conda`).

**Data already on disk (verified intact 2026-10-11):**

| Source path (copy FROM, do not modify originals) | Destination | Size | SHA-256 |
|---|---|---|---|
| `/home/reza/projects/aging_hypothalamus/tabula-muris-senis-droplet-official-raw-obj.h5ad` | `data/raw/tms_droplet_raw.h5ad` | 3.8 GB | `86e45e88f9c475253792f1bbcdf5851a5be0ef35da72fced607ed261aebc099d` |
| `/home/reza/projects/aging_hypothalamus/data/tms_brain.h5ad` | `data/raw/tms_facs_brain_nonmyeloid.h5ad` | 304 MB | `31ee7ffe9fd11fd128925c2dae862399f5c54f27aef9c026dbeaa077ea143339` |

Delete the two 0-byte dead files (`/home/reza/projects/aging_hypothalamus/data/23872787`, `23872787.1`) and record the deletion in the M0 gate report.

**Verified TMS droplet facts (from the file itself):** 245,389 cells × 20,138 genes; ages 1m/3m/18m/21m/24m/30m; sex by age group: 1m M-only, 3m mixed, 18m mixed, 21m F-only, 24m M-only, 30m M-only. **Primary contrast = 3m vs 18m** (only sex-balanced pair). The droplet object has **no brain tissue** and no PCA/UMAP — you compute everything. `subtissue` strings contain trailing spaces — strip all string obs columns in M1 and log before/after value counts.

**Environment:** create from `environment.yml` (Python 3.11, conda-forge + bioconda, strict priority; full package list in PLAN.md §4.1). After M0 succeeds: `conda list --explicit > environment.lock` and commit it. Reproduction always uses the lock.

**Repo conventions (binding, PLAN.md §3):** kernelspec `agescape` on every notebook; gate execution command:
```
jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.kernel_name=agescape --ExecutePreprocessor.timeout=600 <notebook>
```
Determinism: `src.set_seeds()` (seed 42) top of every notebook; explicit `random_state` on all stochastic scanpy calls; UMAP `n_jobs=1`. Every notebook ends with an "outputs written" cell: file path + SHA-256 + git commit + manifest digest for each output. Key tables also land in `results/*.csv` (golden outputs; reference copies committed). log2FC convention everywhere: **old vs young, positive = up with age** — Notebook 02 must assert this with a known-direction gene.

## 2. Gate protocol (how you talk to Reza)

- Finishing a milestone = writing `notes/gates/M<n>-report.md` using the template in §7, then **halting**.
- Resume only when Reza's reply contains the literal phrase `M<n> approved`.
- Conditional feedback ("looks good but fix X") = **not approved**. Implement X, re-run affected checks, re-report.
- While waiting at a gate: you may do non-public prep that commits nothing to the next milestone's analysis (e.g., read the next phase's papers, draft skeletons in `notes/drafts/`). No analysis runs past a gate.
- Reza gate = the weekly check-in (standing rule 5). If a milestone takes >1 week, send a brief mid-milestone progress note in the same format.

## 3. Milestone M0 — Scaffold & verification (target: days 1–2)

**Steps:**
1. Create the repo tree per PLAN.md §3 (incl. `scripts/`, `notes/gates/`). Write `.gitignore` (below). `git init` if not already a repo; initial commit.
2. Copy the two data files to `data/raw/` (paths in §1). `rm` the two 0-byte files. Record OS: `lsb_release -a`, `ldd --version` → manifest.
3. Write `environment.yml` per PLAN.md §4.1. Create env: `conda env create -f environment.yml -n agescape` (or mamba). Register kernel: `conda run -n agescape python -m ipykernel install --user --name agescape`.
4. Smoke test: `conda run -n agescape python -c "import scanpy, anndata, pydeseq2, src; print(scanpy.__version__)"` — fix the import list until clean. Then `conda list --explicit > environment.lock`; commit.
5. Write `scripts/verify_data.py`: walks `data/manifest/*.json` (schema: `{filename, url, collection_id, dataset_version_id, schema_version, size_bytes, sha256, downloaded, notes}`), recomputes SHA-256, prints a table, exits non-zero on mismatch/missing.
6. Write `scripts/fetch_data.py`: reads the same manifest; downloads missing files with `curl -C -` to `data/raw/.partial/`, renames into place only on SHA-256 match; `--only <name>` support. (Local copies are pre-seeded; the script exists for the M7 fresh-clone test and future data.)
7. Write manifest JSONs. For the TMS droplet file, record: local origin path, apparent source (CZ CELLxGENE Tabula Muris Senis droplet raw object), and — when you confirm it in M4's first CELLxGENE query — add the collection/dataset UUID + version. Note in the manifest that UUID confirmation happens at M4.
8. Write `src/__init__.py`, `src/set_seeds.py`, `src/io.py` (load + strip-strings + manifest helpers), `src/plotting.py` (style shared by all notebooks).
9. Write `notes/age-bins.md` from PLAN.md §4.2 (mouse 3m/18m primary; sensitivity runs; human continuous) — **before any DE exists anywhere**. This file is binding.
10. Commit.

**.gitignore:**
```
data/
!data/manifest/
*.h5ad
.env
.ipynb_checkpoints/
__pycache__/
*.pyc
.DS_Store
report/_build/
app/node_modules/
app/dist/
```

**DoD (all must pass):** `python scripts/verify_data.py` exits 0 on both files; smoke import passes; `environment.lock` committed; kernelspec visible (`jupyter kernelspec list | grep agescape`); age-bins.md exists; dead files deleted.

## 4. Milestone M1 — Pipeline: QC → clusters (target: week 1)

**Notebook `notebooks/01-qc-and-clustering.ipynb`.** Steps, in order:
1. Load `tms_droplet_raw.h5ad`. Strip string obs columns; log before/after value counts for `subtissue` etc.
2. Print + save the sex×age cell/mouse table (`results/m1-sex-age-table.csv`). Print cells-per-mouse × age-group.
3. **Per tissue** (16 tissues): age-stratified QC distribution plots (n_genes, n_counts, mt%); propose thresholds per PLAN.md §4.8; apply: min/max genes, mt% ceiling, doublet scoring/removal (per-sample). Log every threshold + rationale in `notes/qc-decisions.md` **in this notebook's first run, before Notebook 02 exists**.
4. Concatenate cleaned tissues. Normalize per cell (log1p) **on the concatenated object**. HVGs **per tissue** (flavor `seurat` or `cell_ranger` — never `seurat_v3`), union → `data/interim/m1-hvg-union.csv`. **Never `sc.pp.scale` outside the HVG subset; PCA uses `zero_center=False` on sparse HVG data.**
5. PCA → neighbors → UMAP (`random_state=42`, `n_jobs=1`) → Leiden (flavor + n_iterations recorded in methods) on the **full object**, HVGs only. Checkpoint `.h5ad` per stage to `data/interim/` (stage boundaries only; record SHA-256 in the outputs cell).
6. Annotation: majority vote of provided `cell_ontology_class` per Leiden cluster; confusion matrix cluster × label; composition table vs the paper's reported proportions (from the TMS paper — read it now, this is your first §10 read).
7. Atlas display UMAP: per-tissue proportional subsample, ≥500-cell/tissue floor, seed 42; composition table of the subsample printed next to the full-object table.
8. Memory watch: after neighbors, print `psutil.virtual_memory()`; if >40 GB used, checkpoint and switch remaining steps to per-tissue. (Guardrails in PLAN.md §6.)

**Figures (minimum):** `figures/m1-umap-celltype.png`, `m1-umap-tissue.png`, `m1-umap-age.png`, `m1-qc-violins-by-age.png`, `m1-cluster-label-confusion.png`, `m1-composition-vs-paper.png`.

**DoD:** gate nbconvert command runs 01 clean from a fresh kernel; all six figures exist; `results/m1-sex-age-table.csv` + confusion matrix + composition table committed to `results/`; qc-decisions.md + age-bins.md (from M0) both non-empty; methods.md started (Leiden/UMAP params). Reza gate.

## 5. Milestone M2 — Mouse DE + aging-signal ranking (target: weeks 2–3)

**Notebook `notebooks/02-de-aging.ipynb`.** Per PLAN.md §4.3–4.4:
1. Eligibility per tissue×cell type: ≥3 mice/group and ≥50 cells/group (3m vs 18m). Write the eligibility table to `results/`; ineligible combos listed.
2. **Pseudobulk** per (mouse × tissue × cell type) → pydeseq2, `~ age_group + sex`. Extract log2FC (old vs young) + CI + FDR. **Calibration trio, in order:** (a) label-permutation null (shuffled age labels within eligible units) → expect ≈0 DE at FDR<0.05 & |log2FC|≥0.25; (b) positive controls (Cdkn2a, Igfbp3, SASP set — build the list from the TMS paper you read in M1) must replicate direction; (c) sex-effect floor: DE by sex within 3m only.
3. Cell-level Wilcoxon on the same units (rank_genes_groups) — secondary.
4. **Ranking** per §4.4: pseudobulk n-DE + median |log2FC| + bootstrap-over-mice CIs; full table → `results/m2-aging-signal-ranked.csv` + report supplement copy. Spearman correlation pseudobulk-vs-cell-level rankings; metric-vs-n_mice and metric-vs-n_cells scatters.
5. Diagnostics on top types: age×sex interaction; per-mouse consistency; leave-one-mouse-out; sex-chromosome enrichment in DE lists; sex-ratio test.
6. **Composition analysis:** fraction-vs-age per tissue×cell type (mouse = unit; fractions + test + caveat about dissociation bias from the TMS paper).
7. Sensitivity: extended bins run (young=1m+3m, old=18m+21m+24m+30m) — hypothesis stated before computing (wider gap → stronger signal; 30m survivorship caveat printed beside it); 18m+24m male-only run.
8. **Assertion cell:** log2FC orientation check on a known up-with-age gene.
9. Update methods.md (DE design, thresholds, calibration results).

**Figures:** `m2-de-volcano-<toptissue>.png` (2–3 exemplar types), `m2-aging-signal-ranking.png` (CIs), `m2-calibration-permutation.png`, `m2-composition-vs-age.png`, `m2-metric-vs-power.png`.

**DoD:** ranked CSV with CIs for all eligible types; calibration trio artifacts (permutation ≈0 DE; positive controls correct; sex floor reported); sensitivity tables; methods.md updated. Reza gate.

## 6. Milestone M3 — Trajectory + paper agreement (target: week 4)

**Notebook `notebooks/03-trajectory-agreement.ipynb`.**
1. Choose pseudotime tissue by the §4.7 written rule; record the choice + justification. dpt with all §4.7 guards: corr(pseudotime, age) reported; graph built on young cells, old projected; connectivity check; root = named marker set; root-choice sensitivity.
2. **Agreement table vs the TMS paper:** first re-read the paper's methods + extract its published aging gene lists (supplementary tables / MSigDB TABULA_MURIS_SENIS aging sets — record exactly where you got them). Build the table: claim | our result | agree/disagree/cannot-test | why it might differ. Claims to test (pre-write this list in methods.md before computing): paper's reported DE direction/counts for its headline cell types; Cdkn2a/composition claims; any AP-1 note. Use the continuous-age bridge model (`~ age_months + sex`) as the methodological bridge — the paper used a continuous-age model, so a pure 3m-vs-18m disagreement is method news, not biology news; say which.
3. Figures: `m3-dpt-vs-age.png`, `m3-agreement-table.png` (or formatted markdown in report/).

**DoD:** dpt guards all present (numbers, not just plots); agreement table complete with explicit cannot-test rows; every disagree has a checked explanation. Reza gate.

## 7. Milestone M4 — Human data + orthologs (target: weeks 5–6)

**Notebook `notebooks/04-human-prep.ipynb`.**
1. **SenNet (primary human aging contrast):** query the SenNet **public** search API (token-free; docs read first per §10). Selection rule (pre-registered): public scRNA-seq (10x) datasets ranked by (donor-age span, then n donors ≥60y, then total n); pick the top 1–2 that also overlap a mouse tissue; record the exact query JSON + date + result count in the manifest and methods.md. **If nothing qualifies: stop, document the search, and escalate to Reza at the gate with options — do not silently relax the rule.** Download processed h5ad/AnnData per dataset (resumable; `.partial` + SHA-256). Record donor age + cohort per donor in the manifest.
2. **Tabula Sapiens v1 (replication cohort):** download the whole-body h5ad from its permanent CELLxGENE dataset URL (`curl -C - --retry 10`; local SHA-256, never ETag; figshare release-1.0 per-tissue zips are the documented fallback if the collection now serves v2 — record which path you used). Tissue subset = the tissues overlapping the mouse panel; **10x-only subset** (platform confound); per-tissue processing only, backed mode, **no whole-object in-memory ops** (>250K cells). `df -h` before and after.
3. Human cohorts built per §4.5.3: pseudobulk per (donor × cell type), `~ age_years + sex`; age continuous. Eligibility + per-type power table (min detectable |log2FC| at FDR 0.1 given n donors) → `results/m4-human-power.csv`.
4. **Orthologs:** Ensembl BioMart 1:1 mouse↔human; record release + query date; cache `data/interim/orthologs-mouse-human.csv`. Retry/backoff; REST fallback. Also fetch the mouse↔NMR mapping now (same pipeline) to de-risk M6.
5. Draft `notes/celltype-mapping.md`: CL-ID-based mouse→human table, one-to-many collapse rules, marker spot-checks, unmatched types with disposition.
6. Confirm + record the CELLxGENE collection/dataset UUIDs for the TMS droplet file in its manifest entry (promise from M0).

**DoD:** verify_data.py passes on all new raw files; human cohorts load per-tissue within RAM; ortholog cache + release recorded; mapping table exists; power table in results/; df -h logged. Reza gate.

## 8. Milestone M5 — Conservation (target: weeks 7–8)

**Notebook `notebooks/05-conservation.ipynb`.** Per PLAN.md §4.5 — read the §10 papers FIRST (Schaum 2020; Tyshkovskiy/Gladyshev — record exact citation; Karakaslar 2023 for AP-1), record citations in methods.md, THEN compute.
1. Two axes per tissue-matched cell-type pair: effect-size Spearman (primary; per-mouse block bootstrap CI) + strict replication fraction (same FDR both species) + separate lenient column. Denominator rule: mapped genes expressed ≥5% of cells in the human type — logged.
2. Permutation null (shuffled human donor ages); report observed−null. Housekeeping guard: with/without RPL/RPS/MT-.
3. Power strata: low-powered human types in a separate flagged section; never read as divergence.
4. Mouse-vs-human effect-size scatter figures per top pair.
5. **Validation vs published programs:** gene-set overlap + sign concordance vs Tyshkovskiy programs + AP-1 (immune types only) — explicit agree/disagree/cannot-test row each, with evidence. Add the "relationship to prior work" positioning paragraph (draft in `notes/drafts/`) for the report.
6. **Divergence table** → `results/m5-divergence-table.csv`: program | cell type | conserved/divergent | evidence | power flag. This is a headline artifact.

**DoD:** conservation + divergence tables with stats in results/; every published program has an evidence-backed verdict row; low-powered section separated; methods.md updated. Reza gate.

## 9. Milestone M6 — Longevity outliers (target: weeks 9–10)

Order matters: **(1) Savina signature figure, (2) AnAge table, (3) Hilton (optional, drop first if squeezed).**
**Notebook `notebooks/06-outliers.ipynb`.**
1. Read Savina 2022 full text (§10); record exact citation + GEO accession(s) **as printed in the paper**. Fetch GSE processed counts (the single supplementary CSV — do NOT go to SRA; record rationale in manifest). Mirror the paper's filtering (e.g., drop sample 2V if the paper did); state every filter.
2. Mouse→NMR ortholog mapping (cached in M4); report per-program coverage; <50% coverage → print it beside the claim, never present the score as comparable.
3. Build conserved programs per §4.6 from M5 outputs (mouse DE → strict human replication → program gene set). Score per cell (`sc.tl.score_genes`, method in methods). **Positive-control gate first:** recover the mouse-skin age effect with the identical pipeline.
4. **Signature figure:** per program, score vs life-stage-normalized age in mouse / human / NMR: per-species slope + 95% CI, species×age interaction p, NMR slope tested by TOST against pre-registered δ (set δ in `notes/methods.md` before running, from the mouse-skin slope scale — justify the value). `figures/m6-signature-programs.png`.
5. Savina reproduction: per-cell-type DE young-vs-old in NMR skin; DE-count-vs-log2FC-cutoff curve vs the paper's claim form; agreement row.
6. AnAge: download the official zipped TSV (record build + date + URL in manifest); build `results/m6-species-table.csv` with availability flags (single-cell / genome / metadata-only) for bowhead whale, tortoise, quahog, rockfish, NMR, human, mouse + full table for the app. Genome-only species: metadata rows only — never fake trajectories.
7. Hilton (optional): processed matrices only; frame as **NMR-vs-mouse immune identity contrast — explicitly NOT an aging contrast** (young NMRs only). One clean figure or drop.
8. `notes/outlier-limits.md`: thin outlier tissue coverage (skin + immune only); no multi-tissue NMR atlas; NMR "aged" cohort ≈ middle-aged in NMR terms — state plainly; ortholog-coverage caveats.

**DoD:** signature figure with slope+CI+interaction+TOST; positive-control gate passed; Savina reproduction row; species table with flags; outlier-limits.md written. Reza gate.

## 10. Milestone M7 — Report & polish (target: weeks 11–12)

1. README per M7 spec in PLAN.md §5 (signature figure above the fold; sex×age honesty moment first screen; 3-command repro: `python scripts/fetch_data.py && conda create -n agescape --file environment.lock && bash scripts/run_all.sh`; per-notebook runtimes; "20 minutes? read this").
2. `scripts/run_all.py`: executes notebooks 01→06 via the gate command in order, then diffs `results/*.csv` against committed golden copies (exact for deterministic columns, documented tolerance for stochastic ones) → prints PASS/FAIL per table.
3. Re-execute everything from clean kernels; regenerate all figures; zero TODOs in repo (`grep -ri todo --include="*.py" --include="*.ipynb" .` must return only intentional hits).
4. Report `report/conserved-aging-programs.md` with the full M7 section list from PLAN.md (abstract, significance, data/code availability, reproducibility statement, **deviations-from-pre-registration log**, **claims that failed robustness**, Limitations, AI-assistance disclosure). Every figure numbered + cited. Limitations ≥ the listed items (correlation≠causation; divergence; TS age range; outlier coverage; survivorship; harvest-time; dropout).
5. **Fresh-clone test:** `git clone` to `/tmp/agescape-fresh` → fetch_data.py → env from lock → run_all.py → all green. This is the acceptance test.
6. Archive the gated commit to Zenodo (concept DOI); put DOI in README + report.
7. Gate report to Reza with the public-go/no-go checklist as its own section. **Public release only on his explicit go-ahead.**

**DoD:** fresh-clone test all green; report complete; deviations log non-empty (even if it says "no deviations" — which would itself be suspicious; check). Reza gate (the big one).

## 11. Milestone M8 — Web app (target: weeks 13–16)

1. Stack choice per M8 spec (PLAN.md §5): static site, Parquet/Arrow assets, budget caps. Write `notes/app-stack.md` with the decision + host limits. Host = Netlify/Vercel/Cloudflare class (NOT GitHub Pages). **Ask Reza at the M7 gate to pick the host + staging access mode** (Cloudflare Access / basic auth / unlisted URL) — come to him with the three options summarized, don't pick for him.
2. `notebooks/07-app-assets.ipynb`: generates every asset from results/ + interim objects per the asset spec (float16, curated gene panels, ≤150K cells/view seeded subsamples); writes `app/data/manifest.json` (asset → notebook → cell → commit → source-notebook SHA-256s). `scripts/check_app_assets.py` fails on any staleness.
3. Features in M8 priority order (PLAN.md §5). Caveats page from the M6/M7 limitations; sex×age confound visible there; provenance tooltips on every number; shareable per-view URLs; SVG/PNG export with methods footnote.
4. Deploy to **staging only**, access-controlled. Reza click-through. Fix what he finds. **Public launch = new explicit decision, recorded, then execute.**
5. Production deploy via the host's GitHub App (OAuth) integration — never CLI tokens in the repo.

**DoD:** staging live + access-controlled; check_app_assets passes; Reza's click-through sign-off; launch only per explicit instruction. Program complete.

## 12. Failure playbook (quick reference)

| Symptom | Do this |
|---|---|
| OOM / RAM >40 GB | Per-tissue mode; PCA/scale on HVG subset only; `zero_center=False`; drop `.raw`; checkpoint; never seurat_v3 HVG; never t-test rank_genes. |
| Truncated download | Resume `curl -C -`; verify SHA-256 before use; never touch `.partial` files directly. |
| Checksum mismatch on local TMS file | Do NOT edit the file. Re-download from CELLxGENE to `.partial/`, verify, replace, log both checksums in the manifest + gate report. |
| BioMart/Ensembl down | Backoff ×5; REST fallback; if still down, pause ortholog steps, continue elsewhere, record in gate report. |
| Leiden↔provided labels disagree badly | Resolution sweep within 0.5–2.0 only; pick best-confusion point; document. Never hand-rename clusters silently. |
| Human type fails power thresholds | Goes to low-powered section. Never averaged in, never called divergent. |
| NMR program coverage <50% | Report coverage next to the claim; downgrade the claim's strength in text; never hide it. |
| SenNet search returns nothing qualifying | Stop; document query + date; escalate at gate with options. |
| Result disagrees with any paper | Agreement table row: claim/result/why. Check confounds first (sex×age, composition, continuous-vs-dichotomous, QC). Never force agreement. |
| "I need a token for this API" | Stop. Check whether a public/token-free path exists (SenNet: yes). Otherwise ask Reza; token goes in `.env`, never the repo. |
| Temptation to change a §4 number | §8 protocol: written justification + robustness check in notes/, then proceed; if not justifiable, halt and report. |
| Anything asking for >30 min on one error | Stop and report. |

## 13. Gate report template (`notes/gates/M<n>-report.md`)

```markdown
# M<n> gate report — <date>
## DoD checklist (item → evidence: command output / file path / cell ref)
- [ ] …
## Numbers worth Reza's attention (3–8 bullets; include the surprising/honest ones)
## Decisions made inside my mandate (thresholds chosen, tissues picked, …)
## Deviations from PLAN.md §4 (with justifications; "none" if none)
## Next milestone: first steps + any questions for you
## Public-go/no-go section (M7/M8 only)
```
Mid-milestone weekly notes use the same template with "partial" in the title.

## 14. Final reminders

- Your output is judged by a skeptical aging-biology scientist with 20 minutes. Boring honesty beats exciting vagueness everywhere they conflict.
- The deviations log, the failed-robustness section, and the sex×age table are features, not admissions. Lead with them.
- When in doubt: halt, write it down, ask. That is always acceptable. Fabricating or quietly improvising never is.
