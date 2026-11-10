# Mature-only composition refit: report

## Question
Do the cell-fraction shifts reported in the Agescape manuscript survive when measured only from maturity to old age, rather than across the full age range including young animals and donors?

## Definitions (cutoffs are NOT equivalent across species)
- Mouse mature = 18 months and older (the paper's mature reference for the drift analysis; middle age in human terms, not equivalent to human 30).
- Human mature = 30 years and older (Reza's maturity definition; the brain completes development around age 30).

## Method
Exact reuse of scripts/run_organ_composition.py (pseudocount 0.5 CLR, OLS of CLR on age plus sex with age-only fallback, cell types kept only when present in at least 6 donors, BH q-values within each organ/cohort stratum). The only change from the original is the donor age filter. No eligibility rule was relaxed. The combined mouse file data/interim/m1-qc-cleaned.h5ad no longer exists, so mouse obs were read from the 16 per-tissue m1-qc-cleaned-*.h5ad files. Validation with --validate: the no-filter refit reproduces the original file. Human rows reproduce exactly (max effect difference 1e-16). Mouse rows reproduce to within 0.0005 on CLR slopes (21 of 136 mouse rows differ by more than 1e-6, all in small-n tissues); the 15 original survivors reproduce to within 1.6 percent on p-values. The small mouse differences are consistent with the per-tissue files having been re-cleaned after the combined file used in the original run was built. Human refit covers lung and brain cohorts only; kidney cohorts are out of scope and contain no original global survivor. NMR was not refit.

## Verdicts (15 original global survivors)
SURVIVES: 5. DIES: 9. UNTESTABLE: 1.
A survivor SURVIVES when its mature-only q-value under the ORIGINAL 431-row global BH family (q_global_bh_orig_family) is below 0.05. q_global_bh_refit is the same correction over the smaller refit-only family, shown for transparency.

| Species | Organ | Cell type | Original effect | Original q | Mature effect | Mature p | Mature q (orig family) | Verdict |
|---|---|---|---|---|---|---|---|---|
| Mouse | Bladder | bladder cell | -0.0613 | 0.0375 | - | - | - | **UNTESTABLE** |
| Mouse | Kidney | kidney capillary endothelial cell | -0.1103 | 0.0044 | -0.1567 | 5.36e-02 | 0.4012 | **DIES** |
| Mouse | Kidney | kidney mesangial cell | -0.0898 | 0.0044 | -0.0845 | 2.29e-01 | 0.7243 | **DIES** |
| Mouse | Limb_Muscle | skeletal muscle cell | +0.0532 | 0.0169 | +0.0475 | 1.82e-01 | 0.6648 | **DIES** |
| Mouse | Liver | B cell | +0.0847 | 0.0185 | +0.2202 | 6.58e-02 | 0.4399 | **DIES** |
| Mouse | Lung | T cell | +0.0619 | 0.0044 | +0.0552 | 2.30e-02 | 0.3138 | **DIES** |
| Mouse | Lung | pulmonary interstitial fibroblast | -0.1292 | 0.0130 | -0.1046 | 1.71e-02 | 0.3057 | **DIES** |
| Mouse | Marrow | late pro-B cell | +0.0834 | 0.0292 | +0.0328 | 4.91e-01 | 0.8676 | **DIES** |
| Mouse | Marrow | precursor B cell | -0.0520 | 0.0290 | -0.0454 | 2.89e-01 | 0.7534 | **DIES** |
| Mouse | Spleen | plasma cell | +0.0758 | 0.0053 | +0.0793 | 1.94e-01 | 0.6810 | **DIES** |
| Human | lung | club cell | -0.0602 | 0.0053 | -0.0804 | 8.48e-05 | 0.0118 | **SURVIVES** |
| Human | lung | multiciliated epithelial cell | -0.0341 | 0.0248 | -0.0470 | 4.52e-04 | 0.0279 | **SURVIVES** |
| Human | brain | L6b glutamatergic neuron of the primary motor cortex | +0.0160 | 0.0185 | +0.0198 | 7.30e-05 | 0.0118 | **SURVIVES** |
| Human | brain | oligodendrocyte | +0.0160 | 0.0064 | +0.0173 | 1.09e-04 | 0.0118 | **SURVIVES** |
| Human | brain | oligodendrocyte precursor cell | -0.0111 | 0.0053 | -0.0119 | 9.55e-05 | 0.0118 | **SURVIVES** |

## Reading the DIES verdicts (power, not disappearance)
All 9 DIES verdicts are mouse. Every one keeps the original sign (9 of 9), and the median mature-only slope is 0.89x the original slope. The standard errors inflate by a median factor of 3.28x (range 1.30 to 5.53x), because the mature-only subset has a 12-month age span (18 to 30 months) instead of 29 months (1 to 30 months) and 7 to 13 donors per tissue instead of 8 to 17. A slope standard error grows as the age spread shrinks, so the same effect size no longer clears the same false-discovery bar. Only marrow late pro-B cell shrinks substantially (0.39x its original slope); the other eight keep 0.81x to 2.60x of their original slope. In other words the mature-only mouse data do NOT show the effects vanishing. They show the effects persisting in point estimate but losing statistical significance because the subset is underpowered. Failure to clear q < 0.05 here does not prove the shifts are maturation artifacts; it shows the mature-only mouse data cannot confirm or refute them at this bar.

## Human lung caveat
SURVIVES here means only that the finding survives the 30+ age restriction under the original pooled method. It does NOT override the manuscript's study-aware reanalysis: human lung club and multiciliated cells do not hold within individual studies (study fixed effects and meta-analysis), and that result stands. The mature-only slopes for club (-0.080 vs -0.060 originally) and multiciliated (-0.047 vs -0.034) are slightly larger than the originals, so the age restriction does not weaken them; the study-confound concern is separate and unchanged.

## New mature-only findings (not among the original 15)
Three cell types that were NOT global survivors in the original analysis clear the original-family global bar (q < 0.05) in the mature-only refit. These are secondary findings from this refit, not pre-specified survivors.

- Human, lung, conventional dendritic cell: mature-only slope +0.0229 (SE 0.0064, p 0.00061, q_orig_family 0.033, n 91 donors); original slope +0.0141 (p 0.0053, q 0.110).
- Mouse, Limb_Muscle, macrophage: mature-only slope -0.1429 (SE 0.0265, p 0.00044, q_orig_family 0.028, n 12 donors); original slope +0.0076 (p 0.65, q 0.870).
- Mouse, Limb_Muscle, mesenchymal stem cell: mature-only slope -0.1005 (SE 0.0168, p 0.00021, q_orig_family 0.018, n 12 donors); original slope +0.0054 (p 0.66, q 0.870).

The two mouse Limb_Muscle types are not noise: their age trajectories are hump-shaped, rising from 1 to 18 months and then falling steadily from 18 to 30 months, so the full-range linear fit is flat (p about 0.65) while the mature-only slope is strongly negative. Leave-one-donor-out mature-only slopes stay between -0.08 and -0.12 for both types, so no single mouse drives them. The human lung conventional dendritic cell rises steadily across the 30+ age bands (mean CLR -1.54 at 30-39 to -0.51 at 70-81); its full-range slope was already positive and nominally significant (p 0.005), and the mature-only slope is steeper (+0.023 vs +0.014).

## Untestable cases
- Mouse Bladder bladder cell: tissue has only 4 mature donors; original rule requires a cell type present in at least 6 donors. Original n donors 8.

## Coverage of the mature-only refit
- Mouse Bladder TMS: 4 mature donors, 0 cell types tested.
- Mouse Fat TMS: 6 mature donors, 5 cell types tested.
- Mouse Heart_and_Aorta TMS: 9 mature donors, 8 cell types tested.
- Mouse Kidney TMS: 11 mature donors, 19 cell types tested.
- Mouse Large_Intestine TMS: 3 mature donors, 0 cell types tested.
- Mouse Limb_Muscle TMS: 12 mature donors, 9 cell types tested.
- Mouse Liver TMS: 7 mature donors, 6 cell types tested.
- Mouse Lung TMS: 10 mature donors, 29 cell types tested.
- Mouse Mammary_Gland TMS: 4 mature donors, 0 cell types tested.
- Mouse Marrow TMS: 13 mature donors, 18 cell types tested.
- Mouse Pancreas TMS: 8 mature donors, 8 cell types tested.
- Mouse Skin TMS: 3 mature donors, 0 cell types tested.
- Mouse Spleen TMS: 10 mature donors, 12 cell types tested.
- Mouse Thymus TMS: 5 mature donors, 0 cell types tested.
- Mouse Tongue TMS: 4 mature donors, 0 cell types tested.
- Human brain HBCC_DLPFC: 64 mature donors, 24 cell types tested.
- Human lung ImmuneAging_BAL: 8 mature donors, 31 cell types tested.
- Human lung Lung_CellRef: 13 mature donors, 29 cell types tested.
- Human lung Lung_Integrated: 91 mature donors, 50 cell types tested.
- Human brain MSSM_DLPFC: 209 mature donors, 25 cell types tested.
- Human brain RADC_DLPFC: 18 mature donors, 23 cell types tested.
- Human brain SEAAD_DLPFC: 20 mature donors, 18 cell types tested.

## Files
- scripts/run_mature_only_composition.py
- results/mature-only-composition-mouse.csv (one row per tested mouse cell type)
- results/mature-only-composition-human.csv (one row per tested human cell type)
- results/mature-only-composition-verdicts.csv (one row per original survivor)
