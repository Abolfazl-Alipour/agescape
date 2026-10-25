# M10 Addendum: Human DLPFC Lifespan Atlas Analysis

Date: 2026-10-22; revised 2026-10-25 with wide confounder audit and CPM primary normalization. Pipeline: `scripts/run_m10_pipeline.py`.
Results: `results/m10-cohort-summary.csv`, `results/m10-celltype-regressions.csv`, `results/m10-confounder-audit.csv`.
Figure: `figures/m10-dlpfc-celltype-forest.png`.

## 1. What this is

M9 scored the mouse aging signature (257 unique genes, 121 up / 70 down human orthologs) in 11 public human cohorts with donor ages on CELLxGENE. Those cohorts differed in age spread and composition.

This addendum analyzes the **DLPFC Lifespan Atlas** (Synapse syn52396927 / syn62064518, controlled access under the AD Knowledge Portal DUA). It profiles **284 control donors aged 0.2 to 97 years** (SD 24.0) across 27 DLPFC cell types as donor-by-cell-type pseudobulk (`lifespan_pseudobulk.rds`).

Data access: token-authenticated Synapse download; raw files live in `data/external/dlpfc-atlas/` (gitignored, DUA-controlled). Only aggregate regression output is committed. Per-donor scatter plots and donor identifiers are excluded from Git and the app to ensure strict DUA compliance.

## 2. Methods

1. Primary normalization is CPM (`log1p(CPM)`). Robustness arm uses raw pseudobulk counts.
2. Covariate adjustment: OLS adjusting for exact Age, Sex, Brain Bank, PMI, and technical PCs 1-3.
3. Confounder audit: wide 27-row table (`results/m10-confounder-audit.csv`, 1 row per cell type) evaluating 7 strata:
   - unadjusted all donors
   - covariate-adjusted all donors (~ Age + Sex + Brain_bank + PMI + TechPC1-3)
   - adult donors (age >= 20y) covariate-adjusted
   - mature adult donors (age >= 40y) covariate-adjusted
   - developmental donors (age < 20y) covariate-adjusted (without Brain_bank, as 100% of donors <20y originate from HBCC)
   - HBCC brain bank only
   - MSSM brain bank only
4. Missing cell types are excluded from donor averages rather than treated as zeros.
5. Multiple testing correction: FDR Benjamini-Hochberg applied across cell types within each stratum.

## 3. Results

### Whole cortex summary

Across 284 donors (age 0.2 to 97.0 years):
- **Unadjusted CPM slope:** +0.00032 per year, 95% CI [0.00007, 0.00056], p=0.0118.
- **Covariate-adjusted slope:** +0.000027 per year, 95% CI [-0.00028, 0.00033], p=0.8626.

When controlling for technical PCs, PMI, and brain bank origin (developmental donors are 100% HBCC), the whole-cortex aging signature slope attenuates completely to null. This demonstrates that the nominal unadjusted rise was driven by technical and demographic confounding between brain banks across the lifespan.

### Cell-type results and confounder audit

Forest plot: `figures/m10-dlpfc-celltype-forest.png`.
Key findings across cell types:
- **Astro:** Unadjusted CPM slope +0.00092/yr (p=0.0040, q_bh=0.0108). In the covariate-adjusted model, p=0.417. In the adult stratum (age >= 20y), slope +0.00104/yr (p=0.0316).
- **Oligo:** Unadjusted CPM slope +0.00075/yr (p=0.00048, q_bh=0.0019). In the covariate-adjusted model, p=0.118. In the adult stratum, p=0.159.
- **Adaptive:** Unadjusted CPM slope -0.00153/yr (p=0.0198, q_bh=0.0412). In the covariate-adjusted model, p=0.157. In the adult stratum, p=0.393.
- **Micro:** Unadjusted CPM slope +0.00154/yr (p=2.25e-5, q_bh=0.0002). In the covariate-adjusted model, slope +0.00142/yr (p=0.0020). In the adult stratum, slope +0.00200/yr (p=0.0019).
- **OPC:** Unadjusted CPM slope +0.00107/yr (p=4.55e-5, q_bh=0.0003). In the covariate-adjusted model, slope +0.00065/yr (p=0.0486). In the adult stratum, slope +0.00129/yr (p=0.0049).

The wide confounder audit confirms that while Micro and OPC retain significant positive slopes in covariate-adjusted models, Astro and Oligo associations are sensitive to brain bank and technical PC adjustments, and Adaptive is not robust in adults.

## 4. Limitations

1. Controlled-access DLPFC data cannot be distributed; public results consist solely of committed aggregate tables and forest plots.
2. Brain bank and developmental age are completely confounded: 100% of donors aged <20y originate from HBCC, while MSSM has zero developmental donors.
3. Pseudobulk CPM normalization handles sequencing depth variations, but cell-level heterogeneity within pseudobulk remains unmeasured.
4. Associations in single cohorts are exploratory hypothesis-generating findings, requiring independent multi-cohort replication.
