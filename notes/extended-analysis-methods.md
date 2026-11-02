# Extended analysis methods

## Estimands

H2 asks whether the relative abundance of annotated cell types changes with age. For sample `i` and cell type `j`, the outcome is the centered log ratio:

`CLR_ij = log((count_ij + 0.5) / sum_j(count_ij + 0.5)) - mean_j(log((count_ij + 0.5) / sum_j(count_ij + 0.5)))`.

The primary model regresses one cell type's CLR values on chronological age, with sex when it varies. Sample, fish, or pool is the inferential unit. Benjamini-Hochberg FDR is applied within each dataset and analysis family.

H3 asks whether cells move relative to a young-adult expression reference. For each cell, raw counts or recovered count proportions are converted to a probability vector `p`. The reference is the pooled young profile `q` for the same annotated cell type. The endpoint is mean Hellinger distance:

`H(p,q) = sqrt(sum_g((sqrt(p_g) - sqrt(q_g))^2) / 2)`.

Young reference cells from the evaluated sample are removed before calculating that sample's distances. A sample must have at least 20 cells in the cell type, and its cross-fitted reference must contain at least 50 cells. Unit models require at least six samples, except descriptive analyses explicitly labeled otherwise.

Depth-adjusted models add mean log UMI per sample. Additional covariates are sequencing platform for Gao worm, incomplete spermatogenesis for zebrafish, and total recovered cells for composition sensitivity analyses.

## Dataset-specific decisions

### AFCA

The primary age window is 30 to 70 days, with day 30 as reference. Five-day cells from FCA are excluded because study and age are inseparable. The released log1p matrix sums to 10,000 after `expm1`, so `expm1(X)` recovers expression proportions. No raw-count Euclidean endpoint is reported.

Head and body measurements from the same age-by-sex pool share a donor cluster in global models. FDR is calculated separately for head and body. Sensitivity analyses use observed records only and total-nuclei adjustment.

### Gao worm

Only wild-type N2 fed OP50 is included. Samples are `N2D1`, `N2D1R`, `N2OPD1`, and the day-6, day-12, and day-14 original, R, and R2 libraries. Day 1 is the primary reference and day 6 is a balanced-reference sensitivity analysis.

`N2D6R`, `N2D6R2`, `N2D12R`, `N2D12R2`, `N2D14R`, and `N2D14R2` are NextSeq 2000. Other included libraries are NovaSeq 6000, based on GEO platform metadata. Platform is included as a covariate. Permutations reassign ages within platform strata. Global permutations reassign ages at the donor level and apply the same assignment across cell types.

Whole-animal and somatic-only analyses exclude different reproductive compartments only in the sensitivity arm. `Germline`, `Embroynic_cells`, and `Sperms` are removed from the somatic arm.

### Roux worm

Raw counts and `annotate_name` labels come from the authors' H5AD. Ages are days 1, 3, 5, 8, 11, and 15. GEO labels libraries as technical replicates. Slopes and R-squared values are descriptive. No biological-replicate p-value, FDR claim, or Bayesian posterior is produced.

### Zebrafish

Technical channels are aggregated by fish. The broad analysis includes all ten fish and three compartments. The fine analysis includes six fish aged 5 to 22 months and the authors' fine cell labels.

The broad expression analysis is not interpreted as purely cell intrinsic because germ-cell-stage composition can change. Incomplete spermatogenesis is encoded for fish `mo27_fish1`, `mo27_fish2`, and `mo27_fish4`; `mo27_fish3` is complete. All earlier fish are complete.

The fine expression analysis enumerates all 720 age assignments per eligible cell type. Leave-one-fish-out slopes are reported. FDR uses the exact permutation p-values across the 12 eligible fine cell types.

## Bayesian models

Bayesian models use standardized sample-level outcomes, unit intercepts, unit age slopes, donor intercepts, unit residual scales, and Student-t residuals. Nuisance covariates are sex, sequencing platform, or total-cell depth as appropriate. Global effects use a Normal prior centered at zero. Vector models evaluate prior SDs 0.25, 0.5, and 1.0 for depth-adjusted effects. Composition models use prior SD 0.5.

Posterior means, 94% highest-density intervals, and probabilities of a positive effect are reported. Models with small cohorts, reference sensitivity, or influential-observation warnings are exploratory. Frequentist sample-level estimates remain the reference.

## Reproduction

Prepare inputs:

```bash
python scripts/fetch_data.py --include-optional
python scripts/prepare_extended_data.py
```

Run frequentist and descriptive analyses:

```bash
python scripts/run_afca_aging.py
python scripts/run_gao_worm_aging.py
python scripts/run_gao_worm_aging.py --vector-only --reference-age 6
python scripts/audit_gao_worm_vector.py
python scripts/audit_gao_worm_vector.py --suffix=-reference6
python scripts/run_roux_worm_aging.py
python scripts/run_zebrafish_testis_aging.py
python scripts/audit_zebrafish_aging.py
```

Prepare Bayesian composition inputs in the core environment:

```bash
python scripts/run_extended_bayesian_composition.py 'Fruit fly' --prepare-only
python scripts/run_extended_bayesian_composition.py 'C. elegans' --prepare-only
python scripts/run_extended_bayesian_composition.py 'Zebrafish' --prepare-only
```

Fit Bayesian models in the Bayesian environment:

```bash
python scripts/run_extended_bayesian_vector.py 'Fruit fly'
python scripts/run_extended_bayesian_vector.py 'C. elegans'
python scripts/run_extended_bayesian_vector.py 'Zebrafish'
python scripts/run_extended_bayesian_composition.py 'Fruit fly' --use-existing
python scripts/run_extended_bayesian_composition.py 'C. elegans' --use-existing
python scripts/run_extended_bayesian_composition.py 'Zebrafish' --use-existing
```

Generate the summary figure:

```bash
python scripts/plot_extended_aging.py
```
