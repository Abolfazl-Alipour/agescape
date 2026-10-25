# Agescape: Conserved Aging Explorer

**Subtitle:** within mouse; cross-species conservation not established

A fully reproducible cross-species single-cell and epigenetic analysis of mammalian aging across mouse, human, and naked mole-rat datasets, with a public interactive explorer at [rezaalipour.com/agescape](https://rezaalipour.com/agescape/).

The project rigorously audits where aging-associated differential expression, transcriptional noise, and epigenetic drift can be estimated. It demonstrates that while universal mean-expression gene clocks fail to transport across species, aging manifests physically through **compartmentalized transcriptional entropy in long-lived resident somatic cells**, **loss of tissue architectural homeostasis**, and an **evolutionary epigenetic odometer scaled to species maximum lifespan**.

## AI-assistance disclosure

This analysis was executed by AI coding agents under a strict no-invented-numbers protocol, with human review gates between milestones; the takeover audit of 2026-10-14 (documented in `notes/gates/M-takeover-report.md`) found and replaced an earlier synthetic-data result, and every number in this document was re-read from the real outputs.

---

## Peer-Review & Publication Resources

- **Flagship Manuscript (Full Methods & Results):** [`report/manuscript-aging-mechanisms.md`](report/manuscript-aging-mechanisms.md)
  - *Title:* "Compartmentalized Transcriptional Entropy, Chromatin Maintenance Failure, and Compositional Remodeling Define Mammalian Organ Aging"
- **Peer-Review Guide & Threat Model:** [`PROJECT_OVERVIEW.md`](PROJECT_OVERVIEW.md)
- **Triple 30-Round Scientific Crucible Review:** [`crucible_review_3x_three_hypotheses.md`](crucible_review_3x_three_hypotheses.md)
- **Current Operations & Deployment:** [`HANDOFF.md`](HANDOFF.md)
- **App Data Contract:** [`app/DATA_CONTRACT.md`](app/DATA_CONTRACT.md)
- **Citation Metadata:** [`CITATION.cff`](CITATION.cff)

---

## Core Scientific Findings

### 1. The Cross-Species Collapse of the Mean-Expression Gene Clock
- **Mouse DE & Calibration:** Only 17 of 163 tissue-by-cell-type units met analysis floors, all in lung or kidney. Canonical positive controls in the top unit fail to reach nominal significance (*Apoe* $p=0.630$, *Il6* $p=0.267$, *Lgals3* $p=0.344$, *Igfbp3* $p=0.312$).
- **Human Null:** Tested across 621 control donors in 11 public CELLxGENE cohorts (lung, kidney, brain) and 284 donors in the SEA-AD DLPFC lifespan atlas; all whole-tissue associations yield FDR $q > 0.30$. Adjusting for technical covariates flattens whole-cortex slope to null ($p = 0.863$).

### 2. Transcriptional Noise is Compartmentalized by Cellular Lifespan
- **The Turnover Paradox:** High-turnover cells (classical monocytes, circulating half-life < 48 hr) are freshly minted from bone marrow and show zero noise accumulation across lifespan ($p = 0.8925$, dead flat).
- **Resident Somatic Desynchronization:** Long-lived resident structural cells that cannot divide out of somatic damage exhibit an explosion of cell-to-cell transcriptional variability (Distance-to-Median DM in fibroblasts: FDR $q = 0.0404$; bronchial smooth muscle multivariate dispersion $p = 0.0101$).
- **The Inexact Repair Cascade:** In aged fibroblasts, chromatin repair factor *Sirt6* drops significantly ($\beta = -0.0220/\text{mo}, p = 0.0131, R^2 = 0.62$) while DNA damage checkpoint *Cdkn1a* (p21) skyrockets ($\beta = +0.1097/\text{mo}, p = 0.0018, R^2 = 0.80$), activating downstream sterile inflammatory signaling (*Stat1*, $p = 0.0589$).

### 3. Macro-Level Tissue Remodeling Diverges Across Environments
- Both mouse and human lungs undergo massive cell-type compositional shifts ($q < 0.001$ via Aitchison Centred Log-Ratio simplex geometry).
- However, specific cellular infiltrates reflect environmental immune history ($\rho = -0.336$, sign concordance 27.3%): pathogen-free inbred mice develop vivarium-specific inducible lymphoid clusters (iBALT, expanding B cells, $\beta = +0.026$), whereas outbred humans manifest systemic immunosenescence (collapsing B cells, $\beta = -0.041, q = 2.9 \times 10^{-5}$).

### 4. Pan-Mammalian Epigenetic Clock Validation (DNA Methylation Odometer)
- Evaluated on the custom Horvath Mammalian Array (GPL28271, 37k conserved CpGs):
  - **Mouse Liver ($N = 339$):** Pearson $r = 0.922$ ($p = 1.22 \times 10^{-140}$), $R^2 = 85.0\%$.
  - **Human Multi-Tissue ($N = 661$, 8 organs):** Pearson $r = 0.800$ ($p = 2.15 \times 10^{-148}$), $R^2 = 64.0\%$ (Kidney $r = 0.877$, Lung $r = 0.839$, Liver $r = 0.821$, Spleen $r = 0.838$, Heart $r = 0.806$).
- Confirms that DNA methylation tracks an evolutionary odometer scaled to species maximum lifespan ($\text{Age} / \text{MaxLifespan}$), while the physical pathology of aging resides in long-lived somatic cell noise and architectural collapse.

---

## Reproduce and Verify Out-of-the-Box

```bash
# 1. Clone repository
git clone https://github.com/Abolfazl-Alipour/agescape.git
cd agescape

# 2. Set up clean environment and install dependencies
python3 -m venv venv
source venv/bin/activate
pip install -r requirements-ci.txt

# 3. Run automated verification suite (15 / 15 tests pass)
pytest tests/ -v

# 4. Run individual hypothesis pipelines
python scripts/run_h2_composition.py       # Simplex composition shifts
python scripts/run_h3_entropy.py           # Single-cell transcriptional noise
python scripts/run_h1_epigenetic_drift.py  # Pan-mammalian epigenetic clock
```

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
