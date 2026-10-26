# Agescape: Conserved Aging Explorer

**Subtitle:** within mouse; cross-species conservation not established

A fully reproducible cross-species single-cell and epigenetic analysis of mammalian aging across mouse, human, and naked mole-rat datasets, with a public interactive explorer at [rezaalipour.com/agescape](https://rezaalipour.com/agescape/).

The project rigorously audits where aging-associated differential expression, transcriptional noise, and epigenetic drift can be estimated. It demonstrates that while universal mean-expression gene clocks fail to transport across species, aging manifests physically through **compartmentalized transcriptional entropy in long-lived resident somatic cells**, **loss of tissue architectural homeostasis**, and an **evolutionary epigenetic odometer scaled to species maximum lifespan**.

## AI-assistance disclosure

This analysis was executed by AI coding agents under a strict no-invented-numbers protocol, with human review gates between milestones; the takeover audit of 2026-10-14 (documented in `notes/gates/M-takeover-report.md`) found and replaced an earlier synthetic-data result, and every number in this document was re-read from the real outputs.

---

## Peer-Review & Publication Resources

- **Flagship Manuscript (Full Methods & Results):** [`report/manuscript-aging-mechanisms.md`](report/manuscript-aging-mechanisms.md)
  - *Title:* "Cell-Type-Private Transcriptional Signatures, Compositional Remodeling, and Epigenetic Clock Dynamics Across Mammalian Aging"
- **Peer-Review Guide & Threat Model:** [`PROJECT_OVERVIEW.md`](PROJECT_OVERVIEW.md)
- **Scientific Walkthroughs & Crucible Reviews:** [`notes/crucibles/`](notes/crucibles/)
  - Final Reconciliation Walkthrough: [`notes/crucibles/walkthrough_final_reconciliation.md`](notes/crucibles/walkthrough_final_reconciliation.md)
  - Final 10-Round Crucible Review: [`notes/crucibles/crucible_walkthrough_review_final_reconciliation.md`](notes/crucibles/crucible_walkthrough_review_final_reconciliation.md)
  - Triple 30-Round Scientific Crucible Review: [`notes/crucibles/crucible_review_3x_three_hypotheses.md`](notes/crucibles/crucible_review_3x_three_hypotheses.md)
- **Current Operations & Deployment:** [`HANDOFF.md`](HANDOFF.md)
- **App Data Contract:** [`app/DATA_CONTRACT.md`](app/DATA_CONTRACT.md)
- **Citation Metadata:** [`CITATION.cff`](CITATION.cff)

---

## Core Scientific Findings

### 1. Cell-Type-Private Transcriptional Aging vs. Whole-Tissue Averaging
- **Whole-Tissue Cortex Null:** In the SEA-AD DLPFC lifespan atlas (284 donors, ages 0.2–97 years), whole-tissue CPM scores attenuate completely to null after adjusting for brain bank, PMI, and technical covariates ($\beta = +0.000027/\text{yr}, p = 0.8626$; unadjusted $\beta = +0.000318/\text{yr}, p = 0.0118$).
- **Cell-Type-Private Cross-Species Replications:** The mouse-derived aging signature transports across species specifically in glial populations:
  - **Oligodendrocyte Precursor Cells (OPCs):** Significant positive slope across all 4 independent human DLPFC cohorts (SEA-AD $p = 1.1 \times 10^{-4}$, surviving Benjamini-Hochberg across all 225 cell-type tests at $q = 0.021$; MSSM $p = 0.008$; HBCC $p = 0.024$; RADC $p = 0.033$).
  - **Microglia:** Significant positive slope across 3 human cohorts (SEA-AD $p = 0.005$, MSSM $p = 0.026$, RADC $p = 0.031$; M10 adult lifespan stratum $\beta = +0.00142/\text{yr}, p = 0.0019$).
- **Neuronal Divergence:** Major excitatory projection neurons exhibit highly significant negative slopes (e.g., `EN_L3_5_IT_3`, adjusted $\beta = -0.00140/\text{yr}, p = 3.79 \times 10^{-7}, q = 1.02 \times 10^{-5}$), demonstrating cell-type-private divergence rather than whole-tissue synchronized aging.

### 2. Naked Mole-Rat as an Evolutionary Negative Control
- **Comparative Senescence Contrast:** While the mouse signature rises across adult mice (whole-body OLS $+0.1344, p = 0.0147$; REML cluster-robust $+0.2053, p = 0.00064$), it remains flat in the negligible-senescence naked mole-rat (*Heterocephalus glaber*, GSE152886 skin: slope $+0.0058$ per unit of normalized maximum lifespan, $p = 0.9332, n = 5$).
- **Equivalence Analysis:** Two-One-Sided-Test (TOST) equivalence confirms flatness within the mouse-derived margin ($p_{\text{TOST}} = 0.0275$ under $\delta_{\min} = 0.202$), and epidermal *Igfbp3* expression is unchanged (221.8 vs 223.8 CPM). At $n = 5$ animals, the 95% CI $[-0.093, +0.100]$ rules out mouse-rate chronological progression but cannot exclude lifespan-scaled drift. Additionally, mouse whole-body slopes exhibit an unpooled OLS slope of $+0.0465$ ($p = 0.1511$) and a plateau in late adulthood (21–30m slope $-0.0056, p = 0.9560$).

### 3. Transcriptional Noise and Fibroblast Marker Audit
- **Rarefied 1,500-UMI Single-Cell Noise:** In mouse lung downsampled to 1,500 UMIs/cell, lung fibroblasts display a nominal unadjusted increase in Distance-to-Median noise ($\beta = +0.00109/\text{mo}, p = 0.0443$, FDR $q = 0.2713$), whereas classical monocytes display flat noise ($\beta = +0.00008/\text{mo}, p = 0.8923$) alongside significant shifts in mean expression.
- **Fibroblast Marker Audit:** In mouse lung fibroblasts, *Cdkn1a* (p21) is significantly upregulated ($\beta = +0.1080/\text{mo}, p = 0.0029, q = 0.0322, R^2 = 0.777$), while *Sirt6* exhibits an exploratory downward trend ($\beta = -0.0136/\text{mo}, p = 0.0245, q = 0.1350$).

### 4. Compositional Remodeling Across Mammals
- **Within-Species Remodeling:** Aitchison Centered Log-Ratio (CLR) reveals significant remodeling in mouse lung (6 cell types across $N = 16$ lung-profiled mice show FDR $q < 0.05$) and human adult lung (4 cell types across $N = 117$ adult donors show FDR $q < 0.05$, including club cell decline $q = 0.0123$ and type 2 alveolar cell expansion $q = 0.0220$).
- **Resolution of Pediatric Confound:** In humans, the apparent B-cell collapse was driven entirely by pediatric donors (ages 2–19); restricting to adult donors ($\ge 20$ yr, $N = 117$) yields a flat B-cell trajectory ($\beta = -0.0027/\text{yr}, p = 0.7601$).
- **Cross-Species Divergence:** Homologous trajectory concordance is not statistically established (sign concordance 3/7, two-sided binomial $p = 1.000$, Spearman $\rho = -0.286$).

### 5. Calibrated Pan-Mammalian Epigenetic Clock 2
- Evaluated with Lu et al.'s (2023) double-exponential inverse link function $\exp(-\exp(-\hat{\eta}))$:
  - **Mouse Liver ($N = 339$):** Pearson $r = 0.923$ ($p = 1.99 \times 10^{-141}$), $R^2 = 85.1\%$, $\text{MAE} = 0.36\text{ years}$ ($4.3\text{ months}$).
  - **Human Multi-Tissue ($N = 661$, 8 organs):** Pearson $r = 0.801$ ($p = 5.16 \times 10^{-149}$), $R^2 = 64.2\%$, $\text{MAE} = 5.06\text{ years}$ (Kidney $r = 0.876, \text{MAE} = 4.10\text{ yr}$; Lung $r = 0.839, \text{MAE} = 4.89\text{ yr}$; Liver $r = 0.817, \text{MAE} = 4.96\text{ yr}$).

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

# 3. Run automated verification suite (16 / 16 tests pass)
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
