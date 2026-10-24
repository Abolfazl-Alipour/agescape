# Agescape: Conserved Aging Explorer

A reproducible single-cell analysis of aging-associated expression in mouse, human, and naked mole-rat datasets, with a public interactive explorer at https://rezaalipour.com/agescape/.

The project tests where aging-associated differential expression can be estimated, how the mouse results compare with published coefficients, and whether a mouse-derived gene score shows similar age associations elsewhere. It does not claim a universal aging program, naked-mole-rat equivalence, or matched replication of every mouse cell-type unit in humans.

## Current findings

- **Mouse eligibility and DE:** 17 of 163 tissue-by-cell-type units met the analysis floor of at least 3 mice and 50 cells in both the 3-month and 18-month groups. All are in lung or kidney; 14 have at least one aging-associated gene. Lung adventitial cells have 55 DE genes and kidney proximal convoluted tubule cells have 53. The floor permits analysis; it is not proof of adequate power.
- **Published comparison:** 13 units agree, 4 partially agree, and 0 disagree with published TMS age coefficients. The kidney proximal tubule comparison has Spearman rho 0.866 across 13,019 shared genes.
- **Program catalog:** 261 direction-specific entries across 257 unique genes: 167 up and 94 down. Significant DE lists overlap little, but this alone does not establish cell-type-private effects.
- **Mouse and NMR scores:** the partly in-sample mouse score rises (slope 0.235 [0.067, 0.403], p=0.008, n=23). No rise is detected in NMR skin (slope 0.017 [-0.257, 0.290], p=0.858, n=5). Five animals, two group-mean ages, skin-only sampling, and 31% coverage preclude an equivalence claim (TOST p=0.329). The TOST margin was derived from the observed mouse slope, not pre-registered.
- **Public human cohorts:** 621 control donors across 11 cohorts show no nominally significant positive whole-cohort score slope. MSSM DLPFC has a nominal decline (p=0.028) that does not survive BH correction across cohorts (q=0.311).
- **DLPFC Lifespan Atlas:** 284 control donors aged 0.2 to 97 years, 27 cell types. Seven have same-direction nominal signals in raw-sum and CPM scoring; six also have BH q<0.05 in both, but adult-only robustness is not universal. The corrected available-cell-type average detects no rise (slope +0.0001464/year, p=0.3619).

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
python scripts/fetch_data.py
python scripts/run_all.py
python scripts/export_app_data.py
python scripts/check_app_data.py
python scripts/check_app_assets.py
```

The M10 controlled-access pipeline requires separately authorized DLPFC inputs and an environment with `rdata`:

```bash
python scripts/run_m10_pipeline.py
```

The app is static and must be served over HTTP:

```bash
python3 -m http.server 8080 --directory app
```

## Scientific reports

- `report/conserved-aging-programs.md`: main analysis, revised with the M9/M10 human results.
- `report/m9-human-cohorts.md`: public human cohort analysis.
- `report/m10-dlpfc-lifespan-atlas.md`: lifespan-atlas analysis and the 2026-10-24 missing-cell-type correction.

## Load-bearing limitations

Cross-sectional associations are not longitudinal or causal effects. TMS has partial sex-by-age confounding and only 17 eligible units in two tissues. Mouse score estimation is partly in-sample. NMR evidence is sparse and skin-only. Human cohorts differ in tissue, composition, expression summaries, and normalization; score magnitudes are not directly comparable across datasets. Missing cell types are not biological zeros. Non-significance does not establish flatness or equivalence, and limited significant-list overlap is not a heterogeneity test.

Earlier project history and audit records remain in `notes/` for accountability.

## License

Code: MIT. Report and figures: CC BY 4.0. Source datasets retain their own terms and access conditions.
