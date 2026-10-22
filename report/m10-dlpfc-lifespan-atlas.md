# M10 Addendum — Gold-Standard Human Lifespan Replication (DLPFC Lifespan Atlas)

Date: 2026-10-22
Status: pipeline complete, results committed, app-integrated
Companion: `report/m9-human-cohorts.md` (the underpowered-cohort pass that this addendum supersedes as the human answer)

---

## 1. What this is

M9 scored the 261-gene mouse aging program (167 up / 94 down, `results/m6-program-genes.csv`) in every public human cohort with donor ages on CELLxGENE. Its verdict was honest but unsatisfying: cohorts of 4 to 213 donors, narrow age spreads, and a pre-registered power floor (at least 6 donors, donor-age SD at least 12 years, M4) that only MSSM DLPFC barely cleared.

This addendum is the test M9 was a placeholder for. The **DLPFC Lifespan Atlas** (Synapse syn52396927, controlled access under the AD Knowledge Portal DUA; folder `DLPFC_lifespan_snRNAseq`, syn62064535) profiles **284 control donors aged 0.2 to 97 years** (SD 24.0) across 27 DLPFC cell types as donor-by-cell-type pseudobulk (counts summed by the authors, `lifespan_pseudobulk.rds`). It clears the M4 floor by an order of magnitude: with SD 24, the minimum detectable 30-year effect at the floor grid point (284 donors, dispersion 0.2) is far below any slope reported below.

Data access: token-authenticated Synapse download on 2026-10-22; raw files live in `data/external/dlpfc-atlas/` (gitignored, DUA-controlled). Only aggregate regression output is committed. Nothing per-donor enters the repo, the app, or this report.

## 2. Methods (what M10 runs, `scripts/run_m10_pipeline.py`)

1. Read the authors' SingleCellExperiment via `rdata`: 27 cell-type count matrices (34,890 genes x 284 donors, identical donor order and gene universe across matrices), donor colData (exact age, sex, ethnicity, PMI, pH, brain bank).
2. Score each donor x cell type exactly as M9: `mean(log1p(count) over the 121 up-program human orthologs present) minus mean(log1p(count) over the 70 down-program orthologs present)` (`data/interim/orthologs-mouse-human.csv`; 3 genes map to both directions in human space, CKB, ID2, STMN1, and net out in the up-minus-down score, matching M6 and M9).
3. Regress score on exact age per cell type. **Four arms, all committed per cell type:**
   - primary: M9-identical recipe, unadjusted OLS, plus Spearman;
   - normalization robustness: same score after per-column CPM rescale before log1p;
   - covariate robustness: OLS adjusted for sex (197 male / 87 female);
   - lifespan robustness: adults only (age at least 18) on the primary recipe, since the cohort spans development.
4. Donors with a zero-count cell-type column (no cells of that type sampled) are excluded from that cell type's regression and counted (rare types lose donors: endothelial L5_ET loses 153 of 284).
5. Whole-cortex headline: per-donor mean of the 27 cell-type scores, regressed on age (M9 cohort-summary style).

## 3. Results

Whole cortex (n=284): slope +0.00026 per year, 95% CI [-0.00012, +0.00063], p=0.186; sex-adjusted p=0.178; Spearman rho +0.080. **With full-lifespan power, the mouse program does not rise cortex-wide in humans.** Figure: `figures/m10-dlpfc-whole-score-vs-age.png`.

Cell-type level (`results/m10-celltype-regressions.csv`, sorted by primary p):

- **Rises with age in both normalization schemes (p < 0.05 in raw and CPM, same sign; the defensible core):** Astro +0.00241/yr (p=1.3e-8), Oligo +0.00177 (p=9.2e-9), PC +0.00110 (p=5.9e-4), IN_VIP +0.00064 (p=2.6e-3); plus OPC in the CPM arm (+0.00102, p=1.2e-4; raw p=0.089).
- **Falls with age in both schemes:** PVM -0.00224/yr (p=3.0e-5), IN_LAMP5_LHX6 -0.00089 (p=3.7e-5), Adaptive immune -0.00077 (p=1.9e-3).
- **Normalization-dependent (direction flips between raw-sum and CPM; not interpretable as aging direction):** Endo (raw -0.00196 p=3.9e-5; CPM +0.00168 p=2.3e-4), EN_L2_3_IT (raw + p=5.7e-4; CPM flat/negative), Micro (raw flat p=0.62; CPM +0.00168 p=4.6e-6), IN_SST (raw - p=0.025; CPM + p=0.014), EN_L3_5_IT_3 (CPM - p=8.4e-5; raw marginal), IN_ADARB2 (CPM + p=7.6e-5; raw flat).
- **Multiple testing:** BH-FDR q < 0.05 (primary model): Astro, Oligo, PVM, IN_LAMP5_LHX6, Endo, EN_L2_3_IT, PC, Adaptive, IN_VIP. Bonferroni p < 0.00185: Oligo, Astro, PVM, IN_LAMP5_LHX6, Endo, EN_L2_3_IT, PC (note Endo and EN_L2_3_IT flip under CPM, so their raw-only significance is treated as not established).
- **Adults-only arm:** conclusions unchanged for the core set (Astro p=3.1e-4, Oligo p=8.8e-5, PVM p=9.8e-3); development is not driving the glial rise.

## 4. What replicates, what dies

Against M9 (`report/m9-human-cohorts.md`):

- **Replicates: glia up.** M9's OPC rise (SEA-AD p=0.0001) and glial/secretory pattern find their cleanest form here: astrocytes and oligodendrocytes rise with age across the entire human lifespan, in both normalization schemes, at Bonferroni significance. This is now the strongest cross-species statement the project can make: the mouse lung-kidney aging program's glial component tracks age in human prefrontal cortex.
- **Dies: cortical neurons down.** M9's MSSM L2/3 IT decline (p=0.0002) does not replicate: in the lifespan atlas EN_L2_3_IT slopes positive on the raw recipe (p=5.7e-4) and flat on CPM. No excitatory neuronal subclass shows a consistent fall. The neuron-down arm of M9 is treated as a cohort artifact.
- **Immune/vascular: downward.** PVM macrophages, perivascular/adaptive immune populations fall with age here; this direction was not visible in M9's cohorts and is reported as new, with the CPM concordance noted.
- **Whole-cortex: flat here.** M9's MSSM whole-cohort decline (p=0.028) does not replicate either (p=0.186). The defensible whole-tissue statement across both milestones: no program-wide cortical drift; the signal is cell-type-private.

## 5. Honest limits

1. The raw-sum versus CPM disagreement is itself a result: pseudobulk library size per donor x cell type carries biological signal (cell number, transcriptional activity) that can masquerade as age direction. Only the 7 concordant cell types are claimed as directional. The app encodes this: concordant units get solid markers, normalization-dependent units get hollow ones.
2. Control donors only (no psychiatric or neurodegenerative diagnosis), but "control" here means no Axis-1 diagnosis; medical comorbidities (cardiovascular disease, diabetes, and so on in the death-record metadata) are part of normal human aging and are not adjusted for. Sex is adjusted for in a dedicated arm.
3. One brain region (DLPFC). The glia-up finding is a prefrontal-cortex statement, not a whole-human one.
4. The mouse program was defined in lung and kidney; its glial rise in human cortex is a cross-species, cross-tissue mapping and should be read as such.
5. Per-donor values are DUA-controlled and deliberately absent from the repo and app; committed aggregates are the regression coefficients, intervals, and test statistics in `results/`.

## 6. Numbers the app commits to (gate-anchored)

- 284 donors, age 0.2 to 97 (SD 24.0), 27 cell types.
- Whole cortex: slope +0.00026/yr, p=0.186, CI spans zero.
- Oligo +0.00177 (p=9.2e-9), Astro +0.00241 (p=1.3e-8), PVM -0.00224 (p=3.0e-5): Bonferroni-level, both schemes.
- 7 cell types concordant across normalization schemes.
- Gene coverage in this atlas: 178 unique human orthologs measured of 191 program slots (114 up + 67 down; coverage 0.9319), the 3 direction-conflict genes (CKB, ID2, STMN1) netting out as in M6 and M9.
