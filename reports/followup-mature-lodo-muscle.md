# Mouse limb-muscle mature-only LODO recompute

## Conflicting source claims

- Old report, `results/mature-only-composition-REPORT.md:48`: "Leave-one-donor-out mature-only slopes stay between -0.08 and -0.12 for both types" (mouse Limb_Muscle macrophage and mesenchymal stem cell).
- Current secondary file, `results/followup-mature-composition-secondary-checks.csv:2` (macrophage): "lodo_slope_min -0.162752438839134" and "lodo_slope_max -0.1186923753674386".
- Current secondary file, `results/followup-mature-composition-secondary-checks.csv:3` (mesenchymal stem cell): "lodo_slope_min -0.1120873064725101" and "lodo_slope_max -0.0835458673102327".

## Method and full-slope validation

The source was `results/followup-mature-composition-donor-clr.csv`, restricted to Mouse, Limb_Muscle, TMS. The subset contains 12 donors at ages 18, 21, 24, and 30 months, with no donors younger than 18 months. CLR values were used exactly as stored; CLR was not recomputed after donor removal, matching `scripts/run_mature_composition_complete.py:335-342`.

For each cell type, OLS regressed `clr_value` on age in months plus `sex_male`, with an intercept, following `fit_row` at `scripts/run_mature_composition_complete.py:60-89`. The age-only fallback applies if the age-plus-sex design is rank-deficient, and fits require at least six donors. All full and 11-donor LODO designs were full rank with age plus sex, so no fallback was needed. Slopes are CLR units per month despite the inherited per-year label (`scripts/run_mature_composition_complete.py:411`). Each LODO fit dropped one of the 12 donors and used the other 11.

Full-donor validation passed the required absolute tolerance of 1e-6 before LODO fitting:

| Cell type | Published full slope | Reproduced full slope | Reproduced minus published |
| --- | ---: | ---: | ---: |
| macrophage | -0.1429451154703087 | -0.14294511547030875 | -5.551115123125783e-17 |
| mesenchymal stem cell | -0.1004571758486359 | -0.1004571758486359 | 0.0 |

## Recomputed LODO ranges

| Cell type | Minimum slope | Held-out donor at minimum | Maximum slope | Held-out donor at maximum |
| --- | ---: | --- | ---: | --- |
| macrophage | -0.16275243883913407 | 30-M-5 (30 months) | -0.1186923753674386 | 30-M-4 (30 months) |
| mesenchymal stem cell | -0.11208730647251011 | 18-M-53 (18 months) | -0.0835458673102327 | 18-M-52 (18 months) |

All 12 refits per cell type are in `results/followup-mature-lodo-muscle.csv`, including age-coefficient p-values. The ranges agree with the current secondary file to floating-point round-off.

## Verdict

- **Macrophage:** the current secondary file is correct. The old report is incorrect: the recomputed range extends well outside -0.12 to -0.08 at both ends.
- **Mesenchymal stem cell:** both claims are numerically correct. Its recomputed range lies within -0.12 to -0.08 and agrees with the current secondary file.
- Overall, the old sentence is incorrect because its "for both types" assertion fails for macrophage; the current secondary file is correct for both cell types.

Model used: `kimi-code/k3.1-pro`.
