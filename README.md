# Agescape: Conserved Aging Explorer

**Subtitle:** within mouse; cross-species conservation not established

A reproducible single-cell analysis of aging-associated expression in mouse, human, and naked mole-rat datasets, with a public interactive explorer at https://rezaalipour.com/agescape/.

The project tests where aging-associated differential expression can be estimated, how the mouse results compare with published coefficients, and whether a mouse-derived gene score shows similar age associations elsewhere. It does not claim a universal aging signature, naked-mole-rat equivalence, or matched replication of every mouse cell-type unit in humans.

## AI-assistance disclosure

This analysis was executed by AI coding agents under a strict no-invented-numbers protocol, with human review gates between milestones; the takeover audit of 2026-10-14 (documented in `notes/gates/M-takeover-report.md`) found and replaced an earlier synthetic-data result, and every number in this document was re-read from the real outputs.

## Start here

- **Peer-review guide:** [`PROJECT_OVERVIEW.md`](PROJECT_OVERVIEW.md)
- **Current operations and deployment:** [`HANDOFF.md`](HANDOFF.md)
- **Main scientific report:** [`report/conserved-aging-programs.md`](report/conserved-aging-programs.md)
- **App data contract:** [`app/DATA_CONTRACT.md`](app/DATA_CONTRACT.md)
- **Original pre-registration:** [`PLAN.md`](PLAN.md) and [`KICKOFF.md`](KICKOFF.md), retained as historical documents

## Current findings

- **Mouse eligibility and DE:** 17 of 163 tissue-by-cell-type units met the analysis floor of at least 3 mice and 50 cells in both the 3-month and 18-month groups. All are in lung or kidney; 14 have at least one aging-associated gene. Lung adventitial cells have 55 DE genes and kidney proximal convoluted tubule cells have 53. The floor permits analysis; it is not proof of adequate power.
- **Calibration controls:** Canonical positive-control markers in the top unit (lung adventitial cell) do not reach nominal significance: Apoe (+0.099, p=0.630), Il6 (+0.778, p=0.267), Lgals3 (+0.467, p=0.344), and Igfbp3 (-0.905, p=0.312). Calibration rests on the permutation null (exact Fisher-Pitman p=0.111 marginal, p=0.133 linked animals; null mean 4.6 against 55 observed) and the sex floor.
- **Published comparison:** 13 units agree, 4 partially agree, and 0 disagree with published TMS age coefficients. The kidney proximal tubule comparison has Spearman rho 0.569 across 12,713 shared genes with 85.1% sign concordance.
- **Signature hierarchy and external overlap:** 
  - Catalog tiers: P0_Raw257 (257 genes), P1_Clean (139 genes: 82 up, 57 down), P2_Clean_SexTrimmed (114 genes), P3_Clean_SexIEGTrimmed (110 genes).
  - Covariance-calibrated $Q^*$ identifies 54 genes with significant between-organ heterogeneity.
  - P1_Clean significantly overlaps established aging signatures: Zhang et al. (2021) FACS 330 GAGs (21 shared, 5.83-fold enrichment, p=6.42e-11), Zhang et al. Droplet 59 GAGs (7 shared, 10.87-fold, p=3.37e-06), and Saul et al. (2022) SenMayo 118 (10 shared, 11.31-fold, p=1.68e-08).
- **Mouse signature scoring:** 
  - Spec (a) All 23 mice whole-body: slope 0.134, p=0.015.
  - Spec (c) REML + CR2 cluster-robust across 19 physical clusters: slope 0.205, p=0.00064 (df=8.48).
  - Spec (d1) Out-of-sample ages (1m, 21m, 24m, 30m): slope 0.107, p=0.012.
  - Spec (d2) Late-adult plateau (21m, 24m, 30m): slope -0.0056, p=0.956.
  - Spec (d3) 14 non-Lung/Kidney holdout tissues across adult ages: slope 0.144, p=0.0167 (df=5.51).
- **NMR score and equivalence:** No rise is detected in NMR skin (slope 0.0035 [-0.093, 0.100], p=0.925, n=5 animals). With a pre-specified mouse-derived margin ($\delta_{\min}=0.202$), TOST equivalence 90% CI is [-0.145, 0.157] (p_TOST=0.027); exact permutation p=1.0. At zero margin, equivalence is not established (p=0.329). In 5-animal epidermis pseudobulk, Igfbp3 expression is unchanged between young (221.8 CPM) and aged (223.8 CPM). Program gene coverage is mouse 1.00, human 0.73, and NMR 0.31.
- **Public human cohorts:** 621 control donors across 11 cohorts show no nominally significant positive whole-cohort score slope. MSSM DLPFC has a nominal decline (slope -0.00019/yr, p=0.028) that does not survive BH correction across cohorts (q=0.311).
- **DLPFC Lifespan Atlas:** 284 control donors aged 0.2 to 97 years, 27 cell types.
  - Whole-cortex CPM unadjusted slope is +0.00032/yr (p=0.0118); when adjusting for sex, brain bank, PMI, and technical PCs 1-3, the slope attenuates to null (p=0.8626).
  - Glial associations (Astro CPM q=0.0108, Oligo CPM q=0.0019, Adaptive CPM q=0.0412) are single-cohort and exploratory; none retain significance in the full confounder-adjusted adult stratum (covar-adj p=0.417 for Astro, p=0.118 for Oligo, p=0.157 for Adaptive).

Every app value is exported from committed result files. Controlled-access DLPFC donor-level values are not committed or exported; only regression aggregates are public.

## Repository layout

```
notebooks/       Numbered analysis notebooks
scripts/         Pipelines, app exporter, and verification gates
src/             Reusable analysis utilities
results/         Committed aggregate result tables
figures/         Generated scientific figures
report/          Main report and M9/M10 human-analysis addenda
notes/           QC, methods, data-scoping, and review records
app/             Static explorer and contracted JSON exports
data/manifest/   Raw-data versions, sources, dates, and checksums
```

Raw and controlled-access data under `data/` are gitignored. No credentials belong in this repository.

## Reproduce and verify

```bash
conda create -n agescape --file environment.lock
conda activate agescape
python scripts/fetch_data.py
python scripts/run_all.py
python scripts/check_app_data.py
python scripts/check_app_assets.py
```

> [!NOTE]
> Running the 146 verification checks requires running `python scripts/fetch_data.py` first because raw inputs such as `data/raw/tms_nature2020_moesm8_dge_results.xlsx` are gitignored. Once fetched, `scripts/check_app_data.py` runs green end to end.

`run_all.py` covers core M1-M6, notebooks 01-06, app export, and verification. M9 and M10 are separate stages with additional dependencies:

```bash
python -m pip install -r requirements-supplemental.txt
python scripts/run_m9_pipeline.py   # public CELLxGENE downloads
python scripts/run_m10_pipeline.py  # requires authorized DLPFC inputs
```

M10 requires separately authorized Synapse files under `data/external/dlpfc-atlas/`. Donor-level DLPFC values and credentials must remain outside Git.

The app is static and must be served over HTTP:

```bash
python3 -m http.server 8080 --directory app
```

## Scientific reports

- `report/conserved-aging-programs.md`: main analysis, revised with the M9/M10 human results and audited signatures.
- `report/m9-human-cohorts.md`: public human cohort analysis.
- `report/m10-dlpfc-lifespan-atlas.md`: lifespan-atlas analysis and wide confounder audit.

## Load-bearing limitations

Cross-sectional associations are not longitudinal or causal effects. TMS has partial sex-by-age confounding and only 17 eligible units in two tissues. Mouse score estimation is partly in-sample. NMR evidence is sparse and skin-only. Human cohorts differ in tissue, composition, expression summaries, and normalization; score magnitudes are not directly comparable across datasets. Missing cell types are not biological zeros. Non-significance does not establish flatness or equivalence, and limited significant-list overlap is not a heterogeneity test.

Earlier project history and audit records remain in `notes/` for accountability.

## License

Code: [MIT](LICENSE). Report and figures: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Source datasets retain their own terms and access conditions.
