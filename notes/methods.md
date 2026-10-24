# Methods — Accurate Execution Record (takeover rebuild, 2026-10-14)

This document is the binding methodological record for the **real-data rebuild** of this
project. It replaces the first run's methods text, parts of which described analyses that
were never executed (see `notes/reviews.md`, "2026-10-14 takeover"). Every quantity below
is read from committed files in `results/`, `data/interim/`, or `data/manifest/`; nothing
here is asserted from memory.

## 0. Data sources (exact)

| Source | Identity | Provenance |
|---|---|---|
| Tabula Muris Senis (TMS) droplet | CELLxGENE collection `0b9d8a04-bb9d-44da-aa27-705bb65b54eb`; dataset `48b37086-25f7-4ecd-be66-f5bb378e3aea`; **dataset_version `084058cc-4f17-43ce-b14e-1278df074013`**; schema 7.1.0 | `data/manifest/tms_droplet_raw.json`; IDs confirmed 2026-10-14 via curation API |
| TMS cohort demographics | 245,389 cells / 23 mice (`results/m1-sex-age-table.csv`, `All` row) | computed by M1 |
| TMS paper age-DGE coefficients | Tabula Muris Consortium, *Nature* 2020, DOI 10.1038/s41586-020-2496-1, Suppl. Table MOESM8 (`DGE_result.tissue_cell.age_coef` / `.bh_p`), archived at `data/raw/tms_nature2020_moesm8_dge_results.xlsx` | downloaded 2026-10-14 (springernature media URL, see manifest) |
| Naked mole-rat (NMR) epidermis | GEO **GSE152886** (Savina et al. 2022, **Aging** — the journal is *Aging*, not *Nature Aging*; PMC9134947; contributor Savina A). Processed counts `GSE152886_raw_data.csv.gz` (20274 genes × cells, 6 samples: 1J/5J/6J young, 2V/3V/4V aged). SRA raw reads **not** downloaded (authors' processed counts sufficient). Sample 2V excluded per the paper (technical outlier) → 5 samples | `data/manifest/GSE152886_nmr_skin_raw.json` |
| Tabula Sapiens (TS) Lung | dataset_version `fddf43a3-b7c8-49e9-ac0a-5e24b390849f`; donors TSP1/2/14/25 = 59/59/60/61y — **adult cohort, no age contrast** | `data/manifest/TS_Lung.json` |
| Tabula Sapiens (TS) Skin | dataset_version **`c240d5bc-a521-4330-a76b-7d566e939d87`**; 4 donors 33/45/59/61y; 17,786 cells; dermis-heavy | `data/manifest/TS_Skin.json` |
| SenNet public index | `search.api.sennetconsortium.org` queried **2026-10-14**: 6,820 published datasets; the public index carries **no donor-age fields**; per-dataset traversal plus Globus transfer is beyond this session's budget | query recorded here; pipeline header in `scripts/run_m5_pipeline.py` |
| AnAge lifespan export | official zipped TSV (`genomics.senescence.info/species/dataset.zip`); mouse 4.0y, NMR 31y, human 122.5y (`results/m6-species-lifespans.csv`) | `data/manifest/anage_dataset.json`, downloaded 2026-10-14 |
| NMR gene names | Ensembl REST `lookup/id` (batch POST), queried **2026-10-14** → `data/interim/nmr-ensg2symbol.tsv` (ENSHGLG → symbol) | TSV header records source and date |
| Ortholog tables | Ensembl BioMart 1:1 orthologs: `data/interim/orthologs-mouse-human.csv` and `data/interim/orthologs-mouse-nmr.csv` — **17,063 pairs each on disk, with identical mouse-gene columns** (a common mouse universe; per-species coverage is relative to that denominator) | verified by Notebook 05 at runtime |

## 1. M1 — QC, normalization, clustering, annotation

- String sanitization of all `obs` columns (`src.io.strip_obs_strings`, before/after counts logged).
- QC filters: ≥500 genes/cell, ≤5,000 genes/cell, ≥1,000 UMI counts/cell. MT filter N/A (0 MT genes in the raw object). Per-tissue Scrublet (expected doublet rate 0.06, seed 42); doublets excluded.
- 238,186 cells retained after QC (`results/m1-qc-summary.csv`, sum of `n_final_cells`; per-tissue retention 67.9–99.7%).
- Normalization: `sc.pp.normalize_total(target_sum=1e4)` + `log1p` on the concatenated cleaned object. HVGs per tissue (2,000/tissue, `flavor='seurat'`), union = 8,421 genes (`data/interim/m1-hvg-union.csv`).
- PCA: sparse TruncatedSVD (`zero_center=False`, n_pcs=50, seed 42) on the HVG union; neighbors n_neighbors=15, n_pcs=30; UMAP `n_jobs=1` (strict reproducibility); Leiden resolution 1.0 (`flavor='igraph'`).
- Annotation: majority vote of the provided `cell_ontology_class` per Leiden cluster; confusion matrix in `results/m1-cluster-label-confusion.csv`; composition table vs the paper in `results/m1-composition-table.csv`.
- **Real demographics caveat (`results/m1-sex-age-table.csv`):** age is strongly confounded with sex in this cohort — the 1m group is entirely male (2 mice), 21m entirely female (2 mice), 24m/30m entirely male (4 mice each). All downstream age contrasts carry this limitation; M2's design (`~ sex + age_group`) is a direct response.

## 2. M2 — Differential expression, ranking, calibration trio

- **Unit:** tissue × `cell_ontology_class`. Eligibility: ≥3 mice and ≥50 cells per group. 163 units evaluated, **17 eligible** (`results/m2-unit-eligibility.csv`).
- **Primary design:** pseudobulk raw counts per (mouse × tissue × cell type); PyDESeq2 design `~ sex + age_group`, 3m vs 18m contrast (positive log2FC = up with age; orientation asserted on Cdkn2a/Spp1).
- **Significance:** BH-FDR < 0.05 and |log2FC| ≥ 0.25.
- **Real ranking (`results/m2-aging-signal-ranked.csv`):** top unit Lung::adventitial cell — 55 DE genes, median |log2FC| 1.335; Kidney::kidney proximal convoluted tubule epithelial cell — 53 / 1.241; Lung::fibroblast of lung — 36 / 1.381. Full per-gene tables: `results/m2_de_*.csv` (17 files, DESeq2-style columns).
- **Calibration trio (`results/m2-calibration-trio-summary.csv`):** for the top unit, observed 55 DE genes vs a label-permutation null mean of 4.6 DE genes, and a sex-effect floor (3m males vs females) of 0 DE genes.
- **Positive controls (`results/m2-calibration-positive-controls.csv`):** Igfbp3 (log2FC −0.90, FDR 0.766), Il6 (0.78, nominal), Apoe (0.10, FDR 0.902), Lgals3 (0.47, FDR 0.783) — none significant after FDR; reported as an honest calibration result, not hidden.
- **Sensitivity (`results/m2-sensitivity-summary.csv`):** extended age bins (young 1m+3m vs old 18–30m) give 120 DE genes for the top unit; the male-only 18m-vs-24m unit (Limb_Muscle::macrophage) gives 757.
- **Composition:** fraction-vs-age per cell type in `results/m2-composition-vs-age.csv`; sex×age confounding note above applies.

## 3. M3 — Trajectory, diffusion pseudotime, and literature agreement

- **Continuum tissue:** Marrow (written rule: strongest graded differentiation pattern). Diffusion map built on young cells, old cells projected; root chosen by stem-marker expression (Kit, Cd34) among young hematopoietic precursors.
- **Real guard result:** DPT vs chronological age Spearman **ρ = 0.0163** (p = 1.18e-03, n = 39,471 cells from `results/m3-dpt-trajectory.csv`) — the axis is a differentiation continuum, **not** an aging clock. No aging claim is attached to DPT.
- **Continuous-age bridge model (`results/m3-continuous-bridge-model.csv`):** per (tissue × cell type) unit, per-mouse pseudobulk `log1p(CPM)` regressed on `age_months` per gene (simple linear regression per gene). Three units have bridge models: Marrow granulocytopoietic cell (12,458 genes), Lung fibroblast of lung (11,556), Kidney kidney proximal convoluted tubule epithelial cell (13,043).
- **Agreement test (`results/m3-agreement-table.csv`):** for each of the 17 powered M2 units, our per-unit effect is compared to the paper's own published age coefficients (MOESM8; paper significance BH < 0.01). Per unit: shared genes, paper-significant count, Spearman ρ, sign concordance among paper-significant genes, nominal recovery fraction. Verdict rule: agree (ρ > 0, concordance ≥ 60%), partial (0 < concordance < 60%), disagree (ρ ≤ 0), cannot-test (< 10 paper-significant genes).
- **Real result:** **13 agree, 4 partial, 0 disagree, 0 cannot-test.** ρ ranges 0.127–0.866 (median 0.342); sign concordance 0.555–0.977; nominal recovery among paper-significant genes is modest everywhere (max 0.526) — sign concordance, not gene-list overlap, is the primary evidence. For 15/17 units our effect is the dichotomous 3m-vs-18m contrast while the paper's coefficient is continuous (recorded per-unit in the `reason` column); the two bridge-model units use our continuous slopes.

## 4. M4 — Human cohort reality check, orthologs, power

- **Orthologs:** Ensembl BioMart 1:1 tables — 17,063 mouse↔human and 17,063 mouse↔NMR pairs on disk, sharing an identical mouse-gene universe (recorded as a provenance caveat in Notebook 05).
- **Cell-type mapping:** curated tissue-matched table `notes/celltype-mapping.md`.
- **Power simulations (`results/m4-human-power.csv`):** 108 configurations over donor n (4–60), donor age SD (8–20y), and biological dispersion (0.1–1.0). Minimum detectable 30-year effect reported per configuration; 72 configurations meet the pre-registered floor (≥6 donors, SD ≥12y), 36 are low-powered.
- **Verified cohort reality:** TS Lung donors 59/59/60/61y (no spread). TS Skin 4 donors 33/45/59/61y (n = 4). SenNet public index (queried 2026-10-14): 6,820 published datasets, no donor-age fields in the public index.
- **Conclusion:** human per-cell-type age-replication is **UNTESTABLE with obtainable data** within this session's budget. No human aging DE was computed in M4.

## 5. M5 — Mouse aging program catalog and (untestable) divergence

- **Program definition:** genes with FDR < 0.05 and |log2FC| ≥ 0.25 (old vs young) across the 17 powered M2 unit DE tables.
- **Real catalog (`results/m5-mouse-program-catalog.csv` / `m5-gene-level-programs.csv`):** 312 unit–gene pairs, 257 unique genes across 17 units (14 units with ≥1 program gene; 3 units zero). Housekeeping guard: ribo/mito fraction per unit, max 0.100 vs 0.017 background; 9 ribo/mito-flagged pairs total.
- **Divergence table (`results/m5-divergence-table.csv`):** all 17 units carry the status `Catalogued (mouse); human arm untestable -> no conservation/divergence call`. **No unit is scored as conserved or divergent.**
- Figure: `figures/m5-program-recurrence.png` (program-gene recurrence across units). The first run's m5 scatter/replication figures were deleted with the fabricated analysis.

## 6. M6 — Longevity-outlier contrast (mouse / human / NMR)

- **Lifespan normalization:** AnAge maxima — mouse 4.0y, NMR 31y, human 122.5y (`results/m6-species-lifespans.csv`).
- **Program (locked before any NMR number):** union of M2 program genes — 167 up + 94 down appearances, 257 unique genes (`results/m6-program-genes.csv`; 4 genes appear in both directions across units). Strict human replication (same-direction slope, p < 0.05, across the 4 TS-skin donors) intersects the program in only **4 genes (2 up, 2 down)**; with only 4 program genes surviving, the pre-set 10-gene floor was not met and the **full mouse program was scored instead, flagged loudly** (deviation D3 below).
- **Arms (`results/m6-lifespan-program-scores.csv`):** Mouse whole-body 23 mice — PRIMARY; Mouse skin 3 mice — UNDERPOWERED (reported, not interpreted); Human TS skin 4 donors (all-cells 15,757 cells and fibroblast-only 3,965) — LOW-POWERED; NMR 5 samples, 10,251 cells — ages are paper-reported group means (young 1.1y, aged 11.3y; individual ages not in GEO); normalized-age span ≈0.03–0.36 of maximum lifespan.
- **Real statistics (`results/m6-statistical-summary.csv`):**
  - Mouse whole-body slope 0.235 [0.067, 0.403], p = 0.008, n = 23 — positive control **PASS**.
  - NMR slope 0.017 [−0.257, 0.290], p = 0.858, n = 5.
  - Mouse-vs-NMR slope difference 0.218 [−0.013, 0.449], p = 0.064 — borderline, not decisive.
  - **Exploratory TOST equivalence at δ = 0.0587 (25% of the observed mouse slope, not pre-registered): p = 0.329 — equivalence NOT established.**
  - NMR aged-vs-young score difference, exact sample-level permutation: p = 0.818.
  - Gates (`results/m6-gates.csv`): program-gene coverage mouse 1.000 / human 0.732 / **NMR 0.314** (assembly/annotation gaps).
- **Savina 2022 reproduction (`results/m6-savina-reproduction.csv`):** DE-count-vs-|log2FC|-cutoff curve on our QC object. At the paper's rule (p < 0.05, |log2FC| ≥ 1.0): 29 genes (20 at FDR < 0.05); at our standard (FDR < 0.05, |log2FC| ≥ 0.25): 560 genes (811 nominal). The paper reports 0 DE genes. **Caveat:** the paper's fold-changes are on the Seurat natural-log scale; converting thresholds between ln and log2 makes this a *partial* reproduction — discrepancy flagged, not adjudicated.
- **Igfbp3 caveat:** the paper's key gene Igfbp3 (Ensembl `ENSHGLG00000037774`) is **absent from the GEO count file** (verified against `data/raw/GSE152886_nmr_skin_raw.csv.gz` in Notebook 06); the paper's central within-cluster finding is **UNVERIFIABLE with obtainable data**.
- **NMR global DE (`results/m6-nmr-de-young-vs-aged.csv`):** 20,273 genes tested; 2,211 nominal p < 0.05; 1,273 FDR < 0.05.
- **Interpretation boundary:** we claim "no detectable program-score change with age in NMR epidermis (5 samples, 2 age-group means, 31% program coverage)" — we do **not** claim proven flatness; equivalence was not established.

## 7. Deviations from pre-registration

1. **D1 — M3 agreement target.** Pre-registration (review pass 3, finding 9) required testing against the paper's published lists rather than a hand-picked marker list. Executed: agreement rebuilt against MOESM8 published coefficients (§3). Justification: a marker-list test would have manufactured false disagreements.
2. **D2 — M4/M5 human arm.** Pre-registration assumed a SenNet primary human aging contrast. The initial 2026-10-14 decision marked per-mouse-unit human replication untestable because SenNet's public index lacked donor ages and TS v1 had four skin donors. M9 and M10 later tested score transport in other cohorts, but still do not supply matched replication for every mouse unit. See the M9/M10 reports; M9 and M10 expression summaries are not identical.
3. **D3 — M6 program definition.** Pre-registration: mouse program validated in human before NMR scoring. Executed: strict human replication at n = 4 donors intersected only 4 program genes (< the 10-gene floor), so the full mouse program (257 unique genes) was scored, flagged in the notebook, this file, and the gate output. Justification: a 4-gene "conserved program" is statistically meaningless; scoring the full program with a loud flag preserves power while recording the deviation.
4. **D4 — M6 mouse skin arm.** Pre-registration positioned mouse skin as the positive-control gate. Executed: the skin arm has n = 3 mice (slope CI [−3.45, 4.09], p = 0.476 — underpowered); the positive control was moved to the whole-body arm (23 mice, p = 0.008). Justification: a positive control must actually have power; the skin arm is still reported for completeness.
5. **D5 — First-run discard.** The first run's M4–M8 outputs (including M6's synthetic `np.random` program scores presented as GSE152886 analysis) were fabricated and have been deleted/discarded; M1–M6 were recomputed from raw data in this takeover. Recorded in `notes/reviews.md`.

## 8. Execution environment

- Conda env `agescape` (lock: `environment.lock`); kernelspec `agescape` for all notebooks.
- Determinism: `src.set_seeds()` (seed 42) at the top of every notebook; explicit `random_state` in pipelines; UMAP `n_jobs=1`.
- Master runner: `scripts/run_all.py` — full recompute of every stage, conditional fetch of `nmr-ensg2symbol.tsv` only if missing, per-stage measured runtimes logged, non-zero exit on any failure.
- Verification: `scripts/verify_data.py` checks SHA-256 of every `data/raw/` file against `data/manifest/*.json`; `scripts/check_app_assets.py` validates the app export.
