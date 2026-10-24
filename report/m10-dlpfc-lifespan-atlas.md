# M10 Addendum: Human DLPFC Lifespan Atlas Analysis

Date: 2026-10-22
Corrected: 2026-10-24
Status: pipeline complete, corrected result files available, app copy aligned
Companion: `report/m9-human-cohorts.md` (a complementary public-cohort analysis, not a directly comparable scoring pipeline)

---

## 1. What this is

M9 scored the mouse aging program, comprising 261 direction-specific entries across 257 unique genes (167 up entries and 94 down entries, `results/m6-program-genes.csv`), in 11 included public human cohorts with donor ages on CELLxGENE. Those cohorts contained 10 to 213 scored donors and differed in age spread and composition. Several met the pre-registered eligibility floor of at least 6 donors and donor-age SD of at least 12 years. That floor permits analysis but does not establish adequate statistical power.

This addendum analyzes the **DLPFC Lifespan Atlas** (Synapse syn52396927, controlled access under the AD Knowledge Portal DUA; folder `DLPFC_lifespan_snRNAseq`, syn62064535). It profiles **284 control donors aged 0.2 to 97 years** (SD 24.0) across 27 DLPFC cell types as donor-by-cell-type pseudobulk (counts summed by the authors, `lifespan_pseudobulk.rds`). The atlas exceeds the eligibility floor, but neither the floor nor its donor count guarantees power for every cell type or effect size.

Data access: token-authenticated Synapse download on 2026-10-22; raw files live in `data/external/dlpfc-atlas/` (gitignored, DUA-controlled). Only aggregate regression output is committed. Nothing per-donor enters the repo, the app, or this report.

## 2. Methods (what M10 runs, `scripts/run_m10_pipeline.py`)

1. Read the authors' SingleCellExperiment via `rdata`: 27 cell-type count matrices (34,890 genes x 284 donors, identical donor order and gene universe across matrices), donor colData (exact age, sex, ethnicity, PMI, pH, brain bank).
2. Score each donor x cell type as `mean(log1p(pseudobulk count) over the 121 up-program human orthologs present) minus mean(log1p(pseudobulk count) over the 70 down-program orthologs present)` (`data/interim/orthologs-mouse-human.csv`). This uses the same ortholog program as M9 but is not an identical scoring recipe. M9 averages per-cell scores before donor summaries, whereas M10 scores author-supplied pseudobulk expression. Their coefficients are not directly comparable. CKB, ID2, and STMN1 map to both directions and enter both means; this algebraic treatment is not evidence that their biological effects simply cancel.
3. Regress score on exact age per cell type. **Four arms, all committed per cell type:**
   - primary: unadjusted OLS on the pseudobulk score, plus Spearman;
   - normalization robustness: the same score after per-column CPM rescaling before log1p;
   - covariate robustness: OLS adjusted for sex (197 male / 87 female);
   - lifespan robustness: adults only (age at least 18) on the primary recipe, since the cohort spans development.
4. Donors with a zero-count cell-type column (no cells of that type sampled) are excluded from that cell type's regression and counted. Rare types lose donors: endothelial L5_ET loses 153 of 284. Missing cell types are not biological zeros.
5. Whole-cortex summary: for each donor, take the mean over only that donor's observed cell-type scores, then regress on age. Between 8 and 27 cell types contribute per donor. Varying composition makes this a descriptive available-cell-type average, not a common complete-panel measurement for every donor. Opposing cell-type estimates can offset numerically, but do not simply net out biologically.

### Code correction record, 2026-10-24

The corrected whole-cortex implementation averages only cell types observed in each donor and excludes absent cell types rather than treating them as zeros. The report and app copy were aligned to the corrected result files on 2026-10-24. The earlier whole-cortex values (+0.00026 per year, p=0.186) were replaced by +0.0001464 per year, p=0.361928 from `results/m10-cohort-summary.csv`. The cell-type regression table was not changed by this documentation alignment.

## 3. Results

Whole cortex (n=284): slope +0.0001464 per year, 95% CI [-0.0001692, +0.0004620], p=0.3619; sex-adjusted slope +0.0001519, p=0.3473; Spearman rho +0.0394, p=0.5079. **This analysis provides no supported evidence of a whole-cortex rise. It does not establish flatness or equivalence.** Figure: `figures/m10-dlpfc-whole-score-vs-age.png`.

### Seven same-direction nominal hits

Seven cell types have the same slope sign and nominal p < 0.05 in both raw-sum and CPM schemes. This rule is not a claim that every hit survives FDR correction in both schemes or remains supported in adults only.

| Cell type | Raw slope per year | Raw p | CPM p | Raw BH q | CPM BH q | Adults-only p |
|---|---:|---:|---:|---:|---:|---:|
| Oligo | +0.00177 | 9.2e-9 | 6.5e-4 | 1.7e-7 | 0.00218 | 8.8e-5 |
| Astro | +0.00241 | 1.3e-8 | 0.00284 | 1.7e-7 | 0.00851 | 3.1e-4 |
| PVM | -0.00224 | 3.0e-5 | 0.00321 | 0.000211 | 0.00867 | 0.00983 |
| IN_LAMP5_LHX6 | -0.000893 | 3.7e-5 | 0.0142 | 0.000211 | 0.0319 | 0.0117 |
| PC | +0.00110 | 5.9e-4 | 1.2e-7 | 0.00226 | 3.3e-6 | 0.00656 |
| Adaptive | -0.000767 | 0.00193 | 0.0475 | 0.00651 | 0.0917 | 0.127 |
| IN_VIP | +0.000637 | 0.00264 | 2.8e-4 | 0.00791 | 0.00110 | 0.0375 |

Six of the seven have BH q < 0.05 in both normalization schemes. Adaptive does not (CPM q=0.0917) and also lacks nominal adults-only support (p=0.127). The adults-only p-values are a separate robustness check, not an FDR-adjusted adults-only confirmation. OPC rises in the CPM arm (+0.00102, p=1.2e-4) but not at nominal significance in the raw arm (p=0.089), so it is not one of the seven.

### Other normalization-dependent results

Several other rows do not show a same-direction nominal signal in both schemes. Endo changes sign (raw -0.00196, p=3.9e-5; CPM +0.00168, p=2.3e-4). EN_L2_3_IT is positive in raw sums (p=5.7e-4) and negative but not nominally significant in CPM (p=0.360). Micro is not nominally significant in raw sums (p=0.616) and is positive in CPM (p=4.6e-6). IN_SST changes sign between schemes. EN_L3_5_IT_3 and IN_ADARB2 have CPM-only nominal signals. No direction is claimed from these discordant or scheme-specific results. Intervals that include zero do not establish neuronal flatness or equivalence.

Full results and all robustness arms: `results/m10-celltype-regressions.csv`.

## 4. Comparison with M9

- **Glial direction:** M9's exploratory OPC and glial estimates are directionally consistent with the M10 Astro and Oligo results. Different cohorts, scoring recipes, and normalization dependence mean this is not a formal direct replication claim.
- **Cortical neurons:** M9's MSSM L2/3 IT decline (nominal p=0.0002) is not directionally consistent with M10 EN_L2_3_IT, which is positive in raw sums and negative but not nominally significant in CPM. The discrepancy is unresolved; it does not justify calling either result a cohort artifact or declaring neurons equivalent to no change.
- **Immune and vascular rows:** PVM, Adaptive, and PC are among the seven same-direction nominal hits, with the FDR and adults-only qualifications shown above. These are cell-type-specific associations, not a general immune-down or vascular-up conclusion.
- **Whole-cohort summaries:** M9's MSSM whole-cohort decline is nominal (p=0.028) and is not supported after Benjamini-Hochberg correction across the 11 included cohorts (q=0.311). M10's available-cell-type average has p=0.362. The two summaries do not support a whole-cortex rise, but their different scoring and composition preclude a direct replication contrast or a claim of no cortical change.

## 5. Honest limits

1. Raw-sum versus CPM disagreement is itself important: pseudobulk library size per donor x cell type can reflect cell number and transcriptional activity. Only seven cell types meet the same-direction nominal rule in both schemes, and not all seven pass every FDR or adults-only check.
2. Control donors only (no psychiatric or neurodegenerative diagnosis), but "control" here means no Axis-1 diagnosis; medical comorbidities in the death-record metadata are not adjusted for. Sex is adjusted for in a dedicated arm.
3. One brain region (DLPFC). Cell-type findings are prefrontal-cortex statements, not whole-human conclusions.
4. The mouse program was defined in lung and kidney. Its use in human cortex is a cross-species, cross-tissue mapping with partial ortholog coverage.
5. M9 and M10 use related but non-identical scoring and inference recipes. Their coefficients and cohort summaries are not directly comparable.
6. The whole-cortex average includes different numbers of observed cell types across donors. Missing cell types are excluded, not zero-filled, and varying composition limits interpretation.
7. Per-donor values are DUA-controlled and deliberately absent from the repo and app; committed aggregates are the regression coefficients, intervals, and test statistics in `results/`.

## 6. Numbers the app commits to

- 284 donors, age 0.2 to 97 (SD 24.0), 27 cell types.
- Whole-cortex available-cell-type average: slope +0.0001464 per year, 95% CI [-0.0001692, +0.0004620], p=0.3619; no supported rise and no equivalence claim.
- Seven cell types with same-direction nominal p < 0.05 in raw and CPM schemes: Oligo, Astro, PVM, IN_LAMP5_LHX6, PC, Adaptive, IN_VIP.
- Six of those seven have BH q < 0.05 in both schemes. Adaptive has CPM q=0.0917 and adults-only p=0.127, illustrating why the seven are not described as uniformly robust.
- Gene coverage in this atlas: 178 unique human orthologs measured of 191 program slots (114 up + 67 down; coverage 0.9319). CKB, ID2, and STMN1 enter both directional means; this is not interpreted as biological cancellation.
