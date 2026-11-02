# Agescape: Conserved Aging Explorer

**Subtitle:** within mouse; cross-species conservation not established

Agescape is a reproducible analysis of mouse, human, naked mole-rat, fruit fly, worm, and zebrafish aging data. It tests how a mouse-derived transcriptional signature behaves across cohorts and cell types, how cell-type fractions within organs change with age, how whole-expression vectors move relative to young references, how transcriptional noise depends on UMI sampling, and how a published mammalian methylation clock performs under explicit predictive metrics. The static explorer is no longer the primary scientific artifact.

## Current manuscript

The primary manuscript is [`report/manuscript-aging-mechanisms.md`](report/manuscript-aging-mechanisms.md). The project design and reviewer guide are in [`PROJECT_OVERVIEW.md`](PROJECT_OVERVIEW.md). Operational details are in [`HANDOFF.md`](HANDOFF.md). `PLAN.md` and `KICKOFF.md` are historical documents.

## Main results

- The primary analysis uses the filtered 139-gene `P1_Clean` signature. The 257-gene `P0_Raw257` score is a labeled sensitivity analysis in every human arm.
- In PsychAD cortex, the primary whole-cortex score was not associated with age before adjustment (slope +0.000186 per year, p = 0.3494) or after covariate adjustment (slope -0.000209, p = 0.4004).
- Three of 252 primary human cell-type tests survived the shared BH correction: MSSM L2/3 neurons, MSSM oligodendrocytes, and PsychAD EN_L3_5_IT_3 neurons. No OPC or microglial result survived that primary correction. Glial transport is a candidate finding.
- Across 13 adequately covered mouse organs, 8 had at least one age-associated cell-type fraction change at FDR 0.05. Three of 11 human cohort datasets had such a change. Naked mole-rat skin had no significant change, but five animals cannot establish stability. See `report/organ-cell-composition.md`.
- Mouse lung fibroblast Cdkn1a increased with age (slope +0.1080 per month, q = 0.0322). Within-cell entropy showed no consistent age-related increase. Off-program expression had no FDR-significant cell-type result, while Bayesian hierarchical analysis gave a weak positive trend with intervals including zero. See `report/transcriptional-disorder.md` and `report/transcriptional-disorder-bayesian.md`.
- The naked mole-rat comparison did not detect a slope in five animals. It cannot establish equivalence or senescence specificity.
- Fly, worm, and zebrafish show cell-fraction changes and cell-type-specific expression changes, but not a universal increase in whole-transcriptome displacement. Worm has the strongest positive vector result, with a smaller and less secure effect under a balanced day-6 reference. Zebrafish has stronger evidence for immune-fraction change than cell-intrinsic drift. See `report/extended-cross-species-aging.md`.
- Universal Clock 2 gave MAE 0.36 years in mouse liver and 5.06 years across 661 human tissue profiles. Predictive R-squared was 0.3169 and 0.6383, respectively. The evaluations are not presented as independent validation.

Every public number is generated from a committed aggregate result table. The manuscript does not use the app as its source of evidence.

## AI-assistance disclosure

AI coding agents assisted with analysis implementation and review under human direction. The 2026-10-14 takeover audit found and replaced an earlier synthetic-data result. Current verification checks compare app exports and selected analysis outputs with committed result tables. These checks do not replace scientific review or independent reproduction.

## Reproduce and verify

```bash
git clone https://github.com/Abolfazl-Alipour/agescape.git
cd agescape
conda create -n agescape --file environment.lock
conda activate agescape
python scripts/fetch_data.py
pytest -v
python scripts/check_app_data.py
python scripts/check_app_assets.py
```

The app checks require fetched public raw inputs because some source files are gitignored. Controlled-access DLPFC donor-level values are not required for the app-data gate.

Run the core pipeline:

```bash
python scripts/run_all.py
```

Run the hierarchical Bayesian transcriptional-disorder analysis in its separate environment:

```bash
conda env create -f environment-bayesian.yml
conda activate agescape-bayesian
python scripts/run_h3_bayesian_disorder.py
```

The Bayesian stage uses a pinned public diagnostic-workflow revision. Large posterior NetCDF files and diagnostic plots stay local and are gitignored. Aggregate posterior and diagnostic tables are committed.

Run the human public-cohort and composition stages with the supplemental dependencies:

```bash
python -m pip install -r requirements-supplemental.txt
python scripts/run_m9_pipeline.py
python scripts/run_human_transport_inference.py
python scripts/run_h2_composition.py
python scripts/run_d3_human_composition_sensitivity.py
python scripts/run_organ_composition.py
python scripts/export_app_data.py
```

The PsychAD stage requires separately authorized inputs under `data/external/dlpfc-atlas/`:

```bash
python scripts/run_m10_pipeline.py
python scripts/generate_m10_stratified_summary.py
python scripts/run_human_transport_inference.py
python scripts/audit_m9_m10_donor_overlap.py
```

Do not commit Synapse credentials, controlled RDS files, donor-level scores, or donor identifiers from controlled data.

Prepare and run the extended fly, worm, and zebrafish analyses:

```bash
python scripts/fetch_data.py --include-optional
python scripts/prepare_extended_data.py
python scripts/run_afca_aging.py
python scripts/run_gao_worm_aging.py
python scripts/run_gao_worm_aging.py --vector-only --reference-age 6
python scripts/audit_gao_worm_vector.py
python scripts/audit_gao_worm_vector.py --suffix=-reference6
python scripts/run_roux_worm_aging.py
python scripts/run_zebrafish_testis_aging.py
python scripts/audit_zebrafish_aging.py
python scripts/plot_extended_aging.py
```

Bayesian extended commands and input preparation are listed in `notes/extended-analysis-methods.md`. Optional manifests download large public objects, so use `python scripts/fetch_data.py --only <target-filename>` when only one arm is needed.

## Repository layout

```text
scripts/          Analysis pipelines, exporters, and verification gates
src/              Shared analysis utilities and signature definitions
results/          Committed aggregate result tables
figures/          Generated scientific figures
report/           Manuscript and analysis reports
notes/            QC, methods, audits, and review records
app/              Static explorer and contracted JSON exports
tests/            Numerical, security, and reproducibility tests
data/manifest/    Public raw-data versions, sources, dates, and checksums
```

## Data-use acknowledgment

The results published here are in whole or in part based on data obtained from the AD Knowledge Portal (NIH-NIA, U24AG061340). See [`notes/adkp-data-use-verification.md`](notes/adkp-data-use-verification.md) for the inspected terms and aggregate-output boundary.

## License

Code: [MIT](LICENSE). Report and figures: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Source datasets retain their own terms and access conditions.
