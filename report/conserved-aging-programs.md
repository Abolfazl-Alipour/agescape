# Conserved Aging Explorer: cross-species single-cell aging signature analysis with an honest power audit

**Subtitle:** within mouse; cross-species conservation not established
**Repository:** `agescape` — cross-species single-cell aging portfolio project
**Date:** 2026-10-14; revised 2026-10-25 with audited statistics and M9/M10 human analyses
**Status:** Public analysis and explorer. All numbers below are read directly from committed files in `results/` and `figures/`. The M9 and M10 addenda provide the full human methods and supersede the initial human-availability assessment.

---

## Abstract

We re-analyzed the Tabula Muris Senis droplet atlas (245,389 cells, 23 mice, 16 tissues) using per-tissue quality control and donor-level pseudobulk differential expression. Seventeen of 163 tissue-by-cell-type units met the eligibility floor for the 3-month versus 18-month contrast, all in lung or kidney; 14 had at least one aging-associated gene. The largest DE counts were 55 in lung adventitial cells and 53 in kidney proximal convoluted tubule cells. Effect estimates agreed with published TMS coefficients in 13 units, partially agreed in four, and disagreed in none. The catalog contains 261 direction-specific entries across 257 unique genes (P0 catalog). We constructed an audited signature hierarchy: P1_Clean (139 genes: 82 up, 57 down), P2_Clean_SexTrimmed (114 genes), and P3_Clean_SexIEGTrimmed (110 genes). Covariance-calibrated $Q^*$ decomposition identified 54 genes with significant between-organ effect heterogeneity. P1_Clean significantly overlaps established aging signatures: Zhang et al. (2021) FACS 330 GAGs (21 shared, 5.83-fold enrichment, p=6.42e-11), Zhang et al. Droplet 59 GAGs (7 shared, 10.87-fold, p=3.37e-06), and Saul et al. (2022) SenMayo 118 (10 shared, 11.31-fold, p=1.68e-08).

The derived mouse aging signature rises with age across 23 mice (whole-body slope 0.134 [95% CI 0.029, 0.239], p=0.015), confirmed by tissue-adjusted cluster-robust regression across 19 physical clusters (slope 0.205, p=0.00064, df=8.48) and out-of-sample generalization in 14 non-Lung/Kidney holdout tissues (slope 0.144, p=0.0167, df=5.51). Within the late-adult subset (21m–30m mice), slope is flat (-0.0056, p=0.956), though an interaction test is non-significant (p=0.563). Five naked mole-rats with skin-only sampling, two group-mean ages, and 31% signature coverage provide no detected rise (slope +0.0058, p=0.933). Under a mouse-derived margin ($\delta_{\min}=0.2023$), TOST equivalence 90% CI is [-0.1450, 0.1567] (p_TOST=0.0274; permutation p=1.0). Whole-epidermis Igfbp3 pseudobulk expression is unchanged between young (221.8 CPM) and aged (223.8 CPM) NMR.

Human analyses now include 621 control donors across 11 public cohorts and 284 control donors in the DLPFC Lifespan Atlas. No nominally significant positive whole-cohort slope was detected in public cohorts; MSSM showed a nominal decline (p=0.028) that did not survive BH correction across 11 cohorts (q=0.311). In the lifespan atlas, unadjusted whole-cortex CPM score rose nominally (+0.00032/yr, p=0.0118), but after adjusting for sex, brain bank, PMI, and technical PCs 1-3, the slope attenuated to null (p=0.8626). Exploratory glial signals (Astro CPM q=0.0108, Oligo CPM q=0.0019, Adaptive CPM q=0.0412) likewise attenuated in the confounder-adjusted adult stratum. These results do not identify a universal aging signature or prove cell-type-private biology.

## Significance statement

Cross-species aging claims need donor-level replication, explicit ortholog and scoring choices, and sensitivity to study design. This project provides a cell-type-resolved mouse analysis, explicit calibration and robustness checks, human tests in public and controlled-access cohorts, and an inconclusive naked-mole-rat comparison. Its central lesson is methodological as well as biological: eligibility is not power, non-significance is not equivalence, and limited overlap of significant lists is not proof of distinct underlying effects.

---

## 1. Data

All raw inputs are version-pinned in `data/manifest/*.json` (URL, accession or dataset UUID, dataset version, schema version, size, SHA-256, download date) and verified by `scripts/verify_data.py`.

- **Tabula Muris Senis, droplet (mouse, primary).** CZ CELLxGENE collection `0b9d8a04-bb9d-44da-aa27-705bb65b54eb`, dataset `48b37086-25f7-4ecd-be66-f5bb378e3aea`, dataset version `084058cc-4f17-43ce-b14e-1278df074013` (schema 7.1.0). The object contains 245,389 cells × 20,138 genes across 16 tissues, 23 mice (6 female, 17 male), and 6 age groups (1m, 3m, 18m, 21m, 24m, 30m). Age and sex are partially confounded by design: 1m, 24m and 30m are male-only and 21m is female-only; the documented primary contrast is 3m vs 18m, the only pair with both sexes in both groups (`results/m1-sex-age-table.csv`).
- **Tabula Sapiens v1 (human, replication cohort attempt).** CZ CELLxGENE Tabula Sapiens collection. Skin (17,786 cells; donors aged 33/45/59/61 y, n=4) and Lung (65,847 cells; donors TSP1/2/14/25 aged 59/61/59/60 y). TS v1 is an adult cohort: it provides no usable donor-age contrast in either tissue.
- **Naked mole-rat skin (longevity outlier).** GEO GSE152886 (Savina et al. 2022, *Aging*; BioProject PRJNA640800). The single supplementary processed-counts file (20,274 genes × 6 samples) was used. Sample 2V was dropped as a technical outlier per the source paper, leaving 5 samples: 1J/5J/6J (young, 1.1 y group mean) and 3V/4V (aged, 11.3 y group mean). Our QC retained 10,251 cells vs the paper's 10,232.
- **AnAge (lifespan denominators).** Official zipped-TSV export from genomics.senescence.info (`results/m6-species-lifespans.csv`): Homo sapiens 122.5 y, Heterocephalus glaber 31 y, Mus musculus 4 y.
- **Orthologs.** Reference 1:1 mouse↔human and mouse↔NMR orthologs from Ensembl 106 Compara and BioMart (`data/reference/orthologs-mouse-nmr.csv`, 14,278 1:1 ortholog pairs).
- **TMS published coefficients (agreement test target).** The TMS paper's Supplementary Table MOESM8 (published age-DGE coefficients), downloaded from the publisher; used as the independent reference for the M3 agreement table.

### Data availability

All datasets above are public. Download URLs, accessions, dataset UUIDs, SHA-256 checksums and download dates are in `data/manifest/`; `python scripts/fetch_data.py` re-acquires every file and `python scripts/verify_data.py` checks each SHA-256. Running the verification suite requires running `fetch_data.py` first because large files such as MOESM8 are gitignored.

---

## 2. Methods (per milestone)

**M1 — TMS processing.** String obs columns stripped. Per-tissue QC with binding thresholds: 500–5,000 genes per cell, ≥1,000 counts; the raw object contains zero mitochondrial genes, so no mt% filter was possible (recorded, not silently skipped). Scrublet doublet removal per tissue. Post-QC retention: 238,186 of 245,389 cells (97.06%; `results/m1-qc-summary.csv`). Normalization (log1p) on the concatenated object; 2,000 HVGs per tissue, union of 8,421 genes; PCA (50 components), kNN graph, UMAP, Leiden clustering; majority-vote annotation against `cell_ontology_class`. All stochastic steps carry seed 42 (`src/set_seeds.py`).

**M2 — Mouse DE and aging-signature ranking.** Unit: tissue × cell-ontology-type, never pooled across tissues. Pseudobulk per (mouse × tissue × cell type) from raw counts; PyDESeq2, design `~ sex + age_group`, primary contrast 3m vs 18m; BH-FDR < 0.05 and |log2FC| ≥ 0.25. Eligibility: ≥3 mice and ≥50 cells per group (`results/m2-unit-eligibility.csv`, 163 units evaluated, 17 eligible in Lung and Kidney). Calibration trio on the top unit:
(a) Exact Fisher-Pitman label-permutation null (9 sex-marginal permutations, exact p=0.1111; 15 linked-animal permutations, exact p=0.1333; null mean 4.6 vs 55 observed DE genes);
(b) Canonical positive controls audited across all 17 units (`results/m2-positive-controls-all-units.csv`);
(c) Sex-effect floor (DE by sex within 3m only: 0 genes).
Cell-type composition modeled via quasi-binomial GLM with total library size weights and Jeffreys boundary shrinkage (`results/m2-composition-model.csv`).

**M3 — Trajectory and literature agreement.** Diffusion pseudotime on Marrow hematopoietic continuum; animal-level mean DPT vs age rho = 0.457 (p=0.065) reported alongside cell-level rho = 0.016 (p=1.13e-3; `figures/m3-dpt-vs-age.png`). Agreement with TMS paper tested against MOESM8 continuous-age coefficients over 12,713 shared genes; continuous bridge model (`log1p(CP10K) ~ age_months + sex`) fit on adult mice.

**M4 — Human power model.** Power grid over n donors (4–50), donor-age SD (8–20 y) and biological dispersion (0.1–0.4) computed before human claims (`results/m4-human-power.csv`). Tabula Sapiens skin (4 donors, SD 8 y) requires minimum detectable 30-year log2FC of 1.474 at dispersion 0.1, proving that TS v1 is unpowered.

**M5 — Signature hierarchy, heterogeneity, and external overlaps.** Program tiers constructed: P0_Raw257 (257 genes), P1_Clean (139 genes: 82 up, 57 down), P2_Clean_SexTrimmed (114 genes), P3_Clean_SexIEGTrimmed (110 genes) (`results/m5-program-hierarchy.csv`). Covariance-calibrated $Q^*$ statistic decomposing within-organ vs between-organ heterogeneity identified 54 genes with significant organ heterogeneity (`results/m5-heterogeneity-qstar.csv`). P1_Clean overlap tested against Zhang et al. (2021) 330 GAGs, 59 Droplet GAGs, and Saul et al. (2022) SenMayo 118 (`results/m5-external-program-overlap.csv`).

**M6 — Signature scoring and cross-species contrast.** Scored via CP10K + `sc.tl.score_genes` (up minus down) on `data/interim/m1-qc-cleaned.h5ad`. Pre-registered specifications:
- Spec (a) All 23 mice whole-body (OLS)
- Spec (b) Unpooled individual mice (N=20, OLS)
- Spec (c) REML + CR2 cluster-robust (G=19 physical clusters)
- Spec (d1) Out-of-sample ages (1m, 21m, 24m, 30m)
- Spec (d2) Late-adult plateau (21m, 24m, 30m)
- Spec (d3) 14 non-Lung/Kidney holdout tissues across adult ages with CR2
NMR arm mapped by stable Ensembl 1:1 ortholog IDs; TOST equivalence tested with pre-specified mouse margin $\delta_{\min}=0.202$ alongside exact permutation test. Five-animal epidermis pseudobulk audited for Igfbp3.

---

## 3. Results

### 3.1 M1/M2 — Mouse aging signal ranking and calibration

Of 163 tissue×cell-type units evaluated against the documented eligibility floor, 17 were eligible — 9 Lung and 8 Kidney (`results/m2-unit-eligibility.csv`). The ranking is headed by Lung::adventitial cell (55 DE genes; 8 mice, 423 cells; median |log2FC| 1.34) and Kidney::kidney proximal convoluted tubule epithelial cell (53 DE genes; 6 mice, 1,793 cells; median |log2FC| 1.24). Three eligible units returned zero DE genes.

The calibration trio on the top unit:
1. Permutation null: exact Fisher-Pitman p=0.1111 (marginal) and p=0.1333 (linked animals), with 55 observed DE genes against a null mean of 4.6.
2. Sex-effect floor: 0 DE genes by sex within 3m mice.
3. Canonical positive controls in the top unit (`results/m2-calibration-positive-controls.csv`):
   - Apoe: log2FC +0.099, p=0.630
   - Il6: log2FC +0.778, p=0.267
   - Lgals3: log2FC +0.467, p=0.344
   - Igfbp3: log2FC -0.905, p=0.312
   None reach nominal significance. Thus, calibration rests on the permutation null and sex floor, not on marker recovery.

### 3.2 M3 — Trajectory and literature agreement

Against MOESM8, verdicts across the 17 powered units are 13 agree, 4 partial, 0 disagree (`results/m3-agreement-table.csv`). Kidney proximal convoluted tubule epithelial cells show Spearman rho = 0.569 across 12,713 shared genes with 85.1% sign concordance. DPT on Marrow correlates with age at animal-level rho = 0.457 (p=0.065) and cell-level rho = 0.016 (p=1.13e-3).

### 3.3 M4/M5 — Human power and signature hierarchy

The M4 power grid shows that Tabula Sapiens skin (4 donors, SD 8 y) requires minimum detectable 30-year log2FC of 1.474 at dispersion 0.1, confirming that TS v1 cannot support an aging claim.

The M5 signature hierarchy defines:
- P0_Raw257: 257 genes (261 entries)
- P1_Clean: 139 genes (82 up, 57 down)
- P2_Clean_SexTrimmed: 114 genes
- P3_Clean_SexIEGTrimmed: 110 genes

Covariance-calibrated $Q^*$ decomposition identified 54 genes with significant between-organ effect heterogeneity. P1_Clean demonstrates highly significant enrichment with independent published aging signatures:
- Zhang et al. (2021) FACS 330 GAGs: 21 shared genes (5.83-fold enrichment, p=6.42e-11)
- Zhang et al. (2021) Droplet 59 GAGs: 7 shared genes (10.87-fold enrichment, p=3.37e-06)
- Saul et al. (2022) SenMayo 118: 10 shared genes (11.31-fold enrichment, p=1.68e-08)

### 3.4 M6 — Signature scoring, holdout validation, and NMR contrast

Scoring the P1_Clean signature (`results/m6-statistical-summary.csv`, Fig. 3):
- Spec (a) All 23 mice whole-body: slope 0.134 [0.029, 0.239], p=0.015.
- Spec (b) Unpooled individual mice (N=20): slope 0.046 [-0.018, 0.111], p=0.151.
- Spec (c) REML + CR2 cluster-robust across 19 physical clusters: slope 0.205 [0.116, 0.295], p=0.00064 (df=8.48).
- Spec (d1) Out-of-sample ages (1m, 21m, 24m, 30m): slope 0.107 [0.029, 0.184], p=0.012.
- Spec (d2) Late-adult subset (21m, 24m, 30m): slope -0.0056 [-0.2338, 0.2225], p=0.956 (interaction test p=0.563).
- Spec (d3) 14 non-Lung/Kidney holdout tissues across adult ages: slope 0.144 [0.038, 0.250], p=0.0167 (df=5.51).

Naked mole-rat contrast:
- NMR skin slope is +0.0058 [-0.1981, 0.2098], p=0.9332 (n=5 animals).
- Mouse-derived margin TOST ($\delta_{\min}=0.2023$): 90% CI [-0.1450, 0.1567], p_TOST=0.0274. Exact permutation p=1.0.
- Slope difference (mouse minus NMR) is +0.1285 [-0.0314, 0.2885], p=0.1153.
- In 5-animal epidermis pseudobulk, Igfbp3 expression is 221.8 CPM in young vs 223.8 CPM in aged animals (unchanged).

---

## 4. Limitations

1. **Correlation, not causation.** Every effect reported is an observational association with age.
2. **Sex×age confound in TMS.** 1m, 24m, and 30m are male-only; 21m is female-only. The primary balanced contrast is 3m vs 18m.
3. **Powered units cover 2 of 16 tissues.** All 17 eligible units are Lung or Kidney units due to atlas sampling depth.
4. **The Tabula Sapiens panel is low-powered.** 4 donors with age SD 8 y cannot test an aging effect.
5. **Human cohorts are not matched per-mouse-unit replication.** M9 and M10 test transport of the mouse-derived score. M9 uses per-cell log1p counts while M10 uses donor pseudobulk CPM, so slope magnitudes are not directly comparable.
6. **NMR signature coverage is 31%.** Less than a third of signature genes were detected in NMR skin; scores are computed on different ortholog subsets.
7. **Equivalence bound caveats.** TOST with mouse margin passes (p=0.0274), but the margin exceeds the mouse point estimate itself, ruling out effects larger than 0.2023 rather than mouse-magnitude effects.
8. **NMR ages are group means.** 5 samples represent two group-mean ages (1.1 y and 11.3 y).
9. **Positive controls were not significant in the top unit.** Calibration rests on the permutation null and sex floor.
10. **M10 glial signals attenuate under confounder adjustment.** Whole-cortex CPM slope attenuates to null (p=0.8626) when adjusting for brain bank, PMI, and technical covariates.

---

## 5. Deviations from pre-registration (log)

- **M4 — SenNet substitution.** SenNet public index lacked donor ages; reported as untestable.
- **M5 — Signature hierarchy & Q* decomposition.** Extended beyond single list to 4 tiers and formal covariance-calibrated $Q^*$ between-organ heterogeneity.
- **M6 — Multi-specification scoring.** Evaluated REML + CR2 cluster-robustness, out-of-sample age holdout, and 14 holdout tissues in addition to whole-body OLS.
- **M6 — Pre-specified TOST margin.** Evaluated against observed mouse slope margin ($\delta_{\min}=0.202$).

---

## 6. Code availability and reproduction

```bash
conda create -n agescape --file environment.lock
python scripts/fetch_data.py
python scripts/run_all.py
```

**AI-assistance disclosure:** this analysis was executed by AI coding agents under a strict no-invented-numbers protocol, with human review gates between milestones; the takeover audit of 2026-10-14 (documented in `notes/gates/M-takeover-report.md`) found and replaced an earlier synthetic-data result, and every number in this document was re-read from the real outputs.

---

## 7. Figures

- **Fig. 1** `figures/m2-aging-signal-ranking.png` — aging-signature ranking with bootstrap CIs.
- **Fig. 2** `figures/m2-calibration-permutation.png` — permutation null vs observed DE counts.
- **Fig. 3** `figures/m6-signature-figure.png` — signature score vs life-stage-normalized age (mouse, human, NMR).
- **Fig. 4** `figures/m6-savina-reproduction.png` — Savina DE-count-vs-cutoff reproduction curve.
- **Fig. 5** `figures/m3-dpt-vs-age.png` — diffusion pseudotime vs age group, Marrow.
- **Fig. 6** `figures/m5-heterogeneity-decomposition.png` — covariance-calibrated $Q^*$ decomposition.
- **Fig. 7** `figures/m10-dlpfc-celltype-forest.png` — DLPFC Lifespan Atlas cell-type forest plot.
