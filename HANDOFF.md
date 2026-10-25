# Current Project Handoff

Last updated: 2026-10-25. This is the current operating handoff. The original milestone manual has been replaced because the project now includes the deployed app, M9/M10 human analyses, and the complete audited statistical pipeline.

## Current state

- Public explorer: https://rezaalipour.com/agescape/
- Analysis repository: https://github.com/Abolfazl-Alipour/agescape
- Website repository: `Abolfazl-Alipour/personalwebsite`, published by Cloudflare from `main`.
- Core analysis, app exports, M9 public human cohorts, and M10 DLPFC aggregates are implemented.
- Latest verification: 146 numeric, semantic, and provenance anchors pass (requires running `fetch_data.py` first for gitignored raw data).
- Desktop and 360 px browser checks pass without JavaScript errors or page overflow.
- Security and DUA status: `figures/m10-dlpfc-whole-score-vs-age.png` (per-donor scatter) was staged and deleted (`git rm`) to ensure strict DUA compliance without rewriting public git history.
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
| M2 pseudobulk DE and ranking | `scripts/run_m2_pipeline.py`, `scripts/run_m2_audit.py` | `results/m2-*.csv` |
| M3 trajectory and TMS agreement | `scripts/run_m3_pipeline.py`, `scripts/run_m3_agreement.py` | `results/m3-*.csv` |
| M4 human power model | `scripts/run_m4_power.py` | `results/m4-*.csv` |
| M5 signature hierarchy & heterogeneity | `scripts/run_m5_pipeline.py`, `scripts/run_m5_heterogeneity.py` | `results/m5-*.csv` |
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

M9 and M10 additionally require the packages in `requirements-supplemental.txt` / `requirements-census.lock`. M9 downloads public CELLxGENE data. M10 requires authorized Synapse inputs under `data/external/dlpfc-atlas/`; do not ask for or commit another user's token.

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
- Calibration rests on the permutation null (exact Fisher-Pitman p=0.111 / 0.133) and the sex floor. Positive controls in the top unit (Apoe, Il6, Lgals3, Igfbp3) do not reach nominal significance.
- Mouse signature: P0 catalog (257 genes), P1_Clean (139 genes: 82 up, 57 down), P2 (114 genes), P3 (110 genes). 54 genes show between-organ heterogeneity via covariance-calibrated $Q^*$.
- Mouse score: whole-body Spec (a) slope 0.134 (p=0.015); Spec (c) cluster-robust CR2 slope 0.205 (p=0.00064, df 8.48); out-of-sample Spec (d1) slope 0.107 (p=0.012); late-adult plateau Spec (d2) slope -0.0056 (p=0.956); holdout tissues Spec (d3) slope 0.144 (p=0.0167, df 5.51).
- NMR: slope 0.0035, TOST equivalence delta_min 0.202, 90% CI [-0.145, 0.157], p_TOST 0.027; slope difference 0.131 (p=0.202); permutation p=1.0. 5-animal pseudobulk Igfbp3 young 221.8 vs aged 223.8 CPM.
- M9 MSSM is nominal p=0.028 and BH q=0.311 across 11 cohort tests.
- M10 absent cell types must remain missing, never zero. Whole cortex unadjusted CPM slope is +0.00032/yr (p=0.0118); after confounder adjustment for brain bank, PMI, and technical PCs, it attenuates to null (p=0.8626).
- M9 and M10 use different expression summaries; their slope magnitudes are not directly comparable.

## Remaining scientific work

- Direct effect-heterogeneity tests across cross-species single-cell datasets.
- A defensible cross-species comparison requires harmonized tissues, scoring, age coverage, and sufficient animals/donors.
- Additional human replication should use comparable expression summaries and pre-specified multiple-testing and robustness criteria.
- Independent usability testing remains outstanding: five computational-biology or hiring reviewers performing the timed tasks on their own phones.

## Repository hygiene check

Before committing, inspect `git status`, keep `data/` and credentials out of the commit, run `git diff --check`, and rerun the relevant gates. Historical audit records remain under `notes/`; they are evidence of corrections, not current instructions.
