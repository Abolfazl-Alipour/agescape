# Current Project Handoff

Last updated: 2026-10-25. This is the current operating handoff. The original milestone manual has been replaced because the project now includes the deployed app and M9/M10 human analyses.

## Current state

- Public explorer: https://rezaalipour.com/agescape/
- Analysis repository: https://github.com/Abolfazl-Alipour/agescape
- Website repository: `Abolfazl-Alipour/personalwebsite`, published by Cloudflare from `main`.
- Core analysis, app exports, M9 public human cohorts, and M10 DLPFC aggregates are implemented.
- Latest verification: 146 numeric, semantic, and provenance anchors pass. Desktop and 360 px browser checks pass without JavaScript errors or page overflow.
- External usability testing has not been performed. Browser checks are builder verification, not an independent usability study.

## Read in this order

1. `README.md`: findings, layout, and reproduction commands.
2. `PROJECT_OVERVIEW.md`: scientific question, design, evidence, and limitations.
3. `report/conserved-aging-programs.md`: main scientific report.
4. `report/m9-human-cohorts.md` and `report/m10-dlpfc-lifespan-atlas.md`: later human analyses.
5. `app/DATA_CONTRACT.md`: app export schemas and figure/data constraints.
6. `PLAN.md` and `KICKOFF.md`: historical pre-registration and initial instructions only.

## Non-negotiable rules

1. Never invent data, statistics, citations, or availability claims. App numbers come from committed result tables through `scripts/export_app_data.py`.
2. Do not treat eligibility as adequate power, non-significance as equivalence, or limited significant-list overlap as a heterogeneity test.
3. Keep donor/animal replication distinct from cell-level associations.
4. Never commit credentials, raw controlled-access data, or per-donor DLPFC values. Only M10 regression aggregates may be public.
5. Preserve the distinction between raw-sum and CPM results, nominal and multiplicity-adjusted evidence, and scheme-specific versus sign-reversing results.
6. Run the export and verification gates after changing results, app data, the data contract, or relevant copy.
7. Do not deploy unverified changes. Reza has approved the current public repository and website deployment path; any new public destination still requires approval.

## Pipeline map

| Stage | Entry point | Main committed outputs |
|---|---|---|
| M1 mouse QC and atlas | `scripts/run_m1_pipeline.py` | `results/m1-*.csv` |
| M2 pseudobulk DE and ranking | `scripts/run_m2_pipeline.py` | `results/m2-*.csv` |
| M3 trajectory and TMS agreement | `scripts/run_m3_pipeline.py`, `scripts/run_m3_agreement.py` | `results/m3-*.csv` |
| M5 mouse program catalog | `scripts/run_m5_pipeline.py` | `results/m5-*.csv` |
| M6 mouse, human skin, and NMR scores | `scripts/run_m6_pipeline.py` | `results/m6-*.csv` |
| M9 public human cohorts | `scripts/run_m9_pipeline.py` | `results/m9-*.csv` |
| M10 DLPFC lifespan atlas | `scripts/run_m10_pipeline.py` | `results/m10-*.csv` |
| App export | `scripts/export_app_data.py` | `app/data/*.json` |
| Verification | `scripts/check_app_data.py`, `scripts/check_app_assets.py` | 146 anchors |

`scripts/run_all.py` covers the core M1-M6 pipeline, notebooks 01-06, app export, and verification. It does not run M9 or M10; those have separate inputs and supplemental dependencies.

## Reproduction and verification

Core environment and pipeline:

```bash
conda create -n agescape --file environment.lock
conda activate agescape
python scripts/fetch_data.py
python scripts/run_all.py
```

M9 and M10 additionally require the two packages in `requirements-supplemental.txt`. M9 downloads public CELLxGENE data. M10 requires authorized Synapse inputs under `data/external/dlpfc-atlas/`; do not ask for or commit another user's token.

Fast app verification without recomputing the science:

```bash
python scripts/export_app_data.py
python scripts/check_app_data.py
python scripts/check_app_assets.py
python3 -m http.server 8080 --directory app
```

Open `http://localhost:8080`. Browser `fetch()` means `file://` is not a supported launch path.

## App architecture and deployment

The app is static HTML, CSS, and vanilla JavaScript with vendored ECharts. There is no backend, account system, analytics, or runtime data analysis. `app/index.html` contains the story; `app/app.js` renders charts and interactions; `app/style.css` contains the visual system.

Deployment from the analysis repository:

```bash
bash scripts/sync_to_website.sh ~/projects/personalwebsite
cd ~/projects/personalwebsite
git add agescape
git commit -m "Describe the Agescape change"
git push origin HEAD
```

The sync script excludes internal provenance and data-contract files. Cloudflare publishes the website after the website-repository push. Verify the live page, including mobile layout and any changed interaction, before calling deployment complete.

## Important corrections to preserve

- The TMS design has partial sex-by-age confounding, not complete confounding. The primary 3m versus 18m model adjusts for sex.
- Five M2 permutations are supportive descriptive checks, not precise calibration. Fifteen bootstrap resamples are not stable confidence-interval evidence and are not drawn as such.
- The mouse score is partly in-sample. NMR has five animals, two group-mean ages, skin-only data, and 31% program coverage; equivalence is not established. The TOST margin was post-hoc, not pre-registered.
- M9 MSSM is nominal p=0.028 and BH q=0.311 across 11 cohort tests.
- M10 absent cell types must remain missing, never zero. The corrected available-cell-type average has slope +0.0001464/year and p=0.3619.
- M9 and M10 use different expression summaries; their slope magnitudes are not directly comparable.
- G15 uses one standard line series per confidence interval and one scatter estimate per cohort. Do not restore the duplicated custom rendering.

## Remaining scientific work

- Direct effect-heterogeneity tests are needed before claiming cell-type-private aging.
- A defensible cross-species comparison requires harmonized tissues, scoring, age coverage, and sufficient animals/donors.
- Additional human replication should use comparable expression summaries and pre-specified multiple-testing and robustness criteria.
- Independent usability testing remains outstanding: five computational-biology or hiring reviewers performing the timed tasks on their own phones.

## Repository hygiene check

Before committing, inspect `git status`, keep `data/` and credentials out of the commit, run `git diff --check`, and rerun the relevant gates. Historical audit records remain under `notes/`; they are evidence of corrections, not current instructions.
