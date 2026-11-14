# Item 1 Report: Cell-Type Fractions in Fly, Worm, and Zebrafish

Date: 2026-11-14
Primary output table: [rev3-composition-fly-worm-zebrafish.csv](file:///home/reza/projects/agescape-revision/results/rev3-composition-fly-worm-zebrafish.csv)
Combined global BH table (624 tests): [rev3-composition-global-bh-624.csv](file:///home/reza/projects/agescape-revision/results/rev3-composition-global-bh-624.csv)
Script: [run_rev3_composition_fly_worm_zebrafish.py](file:///home/reza/projects/agescape-revision/scripts/run_rev3_composition_fly_worm_zebrafish.py)

## 1. Objective and Exact Pipeline Specification

This report executes the exact centered log-ratio (CLR) compositional pipeline used for Mouse and Human in [run_organ_composition.py](file:///home/reza/projects/agescape-revision/scripts/run_organ_composition.py) cleanly from the raw and interim single-cell objects for:

1. **Fruit fly** (*Drosophila melanogaster*, Aging Fly Cell Atlas, `AFCA`): 72 total pools (36 Head pools: 119,833 AFCA annotated cells across 90 eligible cell types; 36 Body pools: 175,710 AFCA annotated cells across 72 eligible cell types; ages 5, 30, 50, and 70 days across both sexes).
2. **Worm** (*Caenorhabditis elegans*, Gao et al. 2024, `Gao_GSE229022`): 12 wild-type N2 hermaphrodite pools (81,196 singleton-filtered annotated cells across 16 eligible cell types; ages 1, 6, 12, and 14 days).
3. **Zebrafish** (*Danio rerio*, Testis, `GSE275361`, primary 6-fish cohort): 6 male fish (`5mo_A`, `5mo_B`, `12mo_A`, `12mo_B`, `mo20`, `22mo`; 36,523 cells total, 36,439 cells across 15 eligible `fine_cell_type`s present in all 6 fish; `liver` with 84 cells in 1 fish is excluded by the `present >= 6` rule).

Pipeline rules applied identically to the Mouse and Human specification:
- **Unit of inference**: donor, animal, pool, or fish (`N = 36`, `36`, `12`, `6`). No cell-level p-values are computed.
- **Eligibility**: cell types detected (`count > 0`) in at least 6 pools/fish (`n_present >= 6`). Primary runs do **not** impose an additional `total_cells >= 50` cutoff (a `total_cells >= 50` filter is run explicitly as a sensitivity check).
- **Primary transformation and regression**: `CLR(counts, pseudocount=0.5) ~ 1 + age + sex_male`, automatically dropping `sex_male` when constant/rank-deficient (Worm hermaphrodites and Zebrafish males).
- **Sensitivities**:
  - Pseudocounts `0.1` and `1.0`.
  - `observed_records_only` (subsetting regression to donors/pools with `count > 0`), `total_cells_adjusted` (`+ log_total_cells`), `platform_adjusted` and `platform_and_total_adjusted` (`+ platform_nextseq` in Worm), `somatic_only` (13 somatic Worm cell types), and `total_cells_ge_50`.
  - **Propeller + limma empirical Bayes check**: logit-transformed pseudocount proportions `log(p / (1 - p))` where `p = (counts + 0.5) / rowSums(counts + 0.5)`, followed by OLS and Smyth (2004) empirical-Bayes variance moderation (`fitFDist`) in Python under both the primary `age + sex` design and the covariate-adjusted `age + sex + log_total_cells (+ platform_nextseq)` design matching [run_reviewer_phase1_propeller.R](file:///home/reza/projects/agescape-revision/scripts/run_reviewer_phase1_propeller.R).

## 2. Explicit Multiple-Testing Families (Stated Before Correction)

1. **Family A (Within-stratum BH)**: Benjamini-Hochberg FDR separately within each of the four primary strata (`pseudocount = 0.5`, `age + sex`):
   - Fly Head (`m = 90` cell types)
   - Fly Body (`m = 72` cell types)
   - Worm Whole animal (`m = 16` cell types)
   - Zebrafish Testis 6-fish fine (`m = 15` cell types)
   Total new primary tests: `90 + 72 + 16 + 15 = 193`.
2. **Family B (Combined Global BH including the existing 431 tests)**: Benjamini-Hochberg FDR across one global family of **`m = 624` primary tests**, uniting the `431` existing Mouse, Human, and Naked mole-rat primary tests in [reviewer-phase1-composition-global-bh.csv](file:///home/reza/projects/agescape-revision/results/reviewer-phase1-composition-global-bh.csv) with the `193` new primary Fly, Worm, and Zebrafish tests.

## 3. Primary Results and Multiple-Testing Summary

| Species | Organ | Cohort / Stratum | Pools or Fish | Eligible Cell Types | Within-Stratum Sig (`q_within < 0.05`) | Global-624 Sig (`q_global < 0.05`) | Min p-value | Min `q_within` | Min `q_global_624` |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Fruit fly | head | AFCA | 36 pools | 90 | 19 | 12 | 4.416e-09 | 3.974e-07 | 1.378e-06 |
| Fruit fly | body | AFCA | 36 pools | 72 | 17 | 14 | 1.702e-09 | 1.226e-07 | 1.062e-06 |
| C. elegans | Whole animal | Gao_GSE229022 | 12 pools | 16 | 2 | 2 | 1.143e-05 | 1.829e-04 | 1.189e-03 |
| Zebrafish | Testis | GSE275361 (6 fish, 5 to 22m, fine) | 6 fish | 15 | 0 | 0 | 7.147e-02 | 9.393e-01 | 2.687e-01 |
| **Combined New** | **4 strata** | **Fly + Worm + Zebrafish** | - | **193** | **38** | **28** | - | - | - |

### Impact on the Combined 624-Test Global BH Family

- Under the original `431`-test Mouse + Human + NMR global family (`results/reviewer-phase1-composition-global-bh.csv`), `15` tests passed `q < 0.05`.
- When the `193` Fly, Worm, and Zebrafish primary tests are combined with the `431` existing tests (`m = 624` total tests), **`45` tests pass global BH `q < 0.05`**:
  - `28` new Fly (`12` head, `14` body) and Worm (`2`: `ARSC`, `Spermatheca`) hits.
  - `17` existing Mouse and Human hits (all `15` prior global hits stay below `0.05`, and 2 borderline tests move from `q_431` just above `0.05` to `q_624 < 0.05` due to the higher BH rank density of strong Fly/Worm p-values: Mouse Lung `vein endothelial cell`, `p = 0.002720`, `q_431 = 0.068955 -> q_624 = 0.039041`; Human Kidney `Kidney_LivingDonor` `kidney cortex collecting duct epithelial cell`, `p = 0.002478`, `q_431 = 0.066741 -> q_624 = 0.037708`).

### Notable Significant Shifts in New Strata (Primary CLR `pseudocount = 0.5`)

- **Fly Head** (top hits passing `q_global_624 < 0.05`, slopes in CLR units/day):
  - `epithelial cell_head`: `+0.01366` (SE `0.00173`, `p = 4.42e-09`, `q_within = 3.97e-07`, `q_624 = 1.38e-06`)
  - `skeletal muscle of head`: `+0.01127` (SE `0.00150`, `p = 1.21e-08`, `q_within = 5.42e-07`, `q_624 = 2.51e-06`)
  - `photoreceptor cell R7`: `-0.01240` (SE `0.00265`, `p = 4.83e-05`, `q_within = 0.00131`, `q_624 = 0.00251`)
  - `pigment cell`: `+0.01199` (SE `0.00260`, `p = 5.84e-05`, `q_within = 0.00131`, `q_624 = 0.00260`)
- **Fly Body** (top hits passing `q_global_624 < 0.05`, slopes in CLR units/day):
  - `adult fat body_body`: `-0.01848` (SE `0.00225`, `p = 1.70e-09`, `q_within = 1.23e-07`, `q_624 = 1.06e-06`)
  - `prefollicle cell/stalk follicle cell`: `+0.01791` (SE `0.00321`, `p = 3.35e-06`, `q_within = 1.21e-04`, `q_624 = 5.23e-04`)
  - `stalk follicle cell`: `+0.02194` (SE `0.00423`, `p = 1.06e-05`, `q_within = 2.54e-04`, `q_624 = 1.19e-03`)
  - `adult Malpighian tubule principal cell of lower ureter`: `-0.02552` (SE `0.00562`, `p = 7.08e-05`, `q_within = 0.00102`, `q_624 = 0.00273`)
- **C. elegans Whole animal** (primary hits passing both within-stratum and global-624 BH, slopes in CLR units/day):
  - `ARSC` (anterior/rectal/arcade support cell cluster): `+0.21600` (SE `0.02691`, `p = 1.14e-05`, `q_within = 1.83e-04`, `q_624 = 0.00119`)
  - `Spermatheca`: `-0.13662` (SE `0.01899`, `p = 2.94e-05`, `q_within = 2.35e-04`, `q_624 = 0.00175`)
- **Zebrafish Testis 6-fish fine** (15 cell types, 5 to 22 months, slopes in CLR units/month):
  - Zero cell types pass either within-stratum or global BH. Smallest p-value is `immune` (`+0.07069`, SE `0.02901`, 95% CI `-0.00986 to +0.15125`, `p = 0.07147`, `q_within = 0.93929`, `q_global_624 = 0.26868`). With `n = 6` fish (`df_residual = 4`), intervals are wide and do not establish equivalence.

## 4. Sensitivity Checks: Pseudocounts (0.1, 0.5, 1.0) and Propeller/Limma

### 4.1 Pseudocount Sensitivity (`0.1`, `0.5`, `1.0`)

Across all four primary strata, CLR slope estimates are strongly invariant to pseudocount choice, and major conclusions remain unchanged:
- **Fly Head (`m = 90`)**: `17`, `19`, and `19` within-stratum hits at pseudocounts `0.1`, `0.5`, `1.0` respectively; all `17` hits at `0.1` are shared across all three pseudocounts; slope correlation vs `0.5` is `r = 0.9873` (`0.1`) and `r = 0.9975` (`1.0`).
- **Fly Body (`m = 72`)**: `17`, `17`, and `21` within-stratum hits at pseudocounts `0.1`, `0.5`, `1.0`; all `17` primary hits at `0.5` remain significant at `0.1` and `1.0`; slope correlation vs `0.5` is `r = 0.9694` (`0.1`) and `r = 0.9915` (`1.0`).
- **C. elegans (`m = 16`)**: exactly `2` within-stratum hits (`ARSC`, `Spermatheca`) at all three pseudocounts (`r = 0.9996` across all three).
- **Zebrafish 6-fish fine (`m = 15`)**: `0` within-stratum hits at all three pseudocounts (`min q_within >= 0.928`; `r = 0.9999`).

### 4.2 Propeller + Limma Empirical Bayes Check

We implemented logit-transformed proportions `(counts + 0.5) / rowSums(counts + 0.5)` with Smyth (2004) empirical-Bayes moderated t-tests (`propeller_logit_eBayes`) under both the unadjusted `age + sex` design and the covariate-adjusted `age + sex + log_total_cells (+ platform_nextseq)` design.

1. **Cross-check against committed R `speckle`/`limma` results** (`results/reviewer-phase1-composition-propeller.csv`, which used the `log_total_cells` design):
   - **C. elegans** (`16` cell types): `max |logit slope diff| = 6.94e-16`, `r(-log10 p) = 0.999991`, `100%` `q < 0.05` agreement (`5/16` significant in both Python and R).
   - **Fly Body** (`72` cell types): `max |logit slope diff| = 1.80e-15`, `r(-log10 p) = 0.999980`, `100%` `q < 0.05` agreement (`32/72` significant in both Python and R).
   - **Fly Head** (`90` cell types): `max |logit slope diff| = 6.00e-16`, `r(-log10 p) = 0.999982`, `100%` `q < 0.05` agreement (`25/90` significant in both Python and R).
   - **Zebrafish 10-fish broad reconciliation** (`3` compartments): `max |logit slope diff| = 9.71e-17`, `r(-log10 p) = 0.997482`, `100%` `q < 0.05` agreement (`2/3` significant in both Python and R).
2. **Primary `age + sex` design (`propeller_limma_primary_design`)**:
   - **Fly Head**: `19` significant (`q_within < 0.05`), **Fly Body**: `30` significant, **C. elegans**: `4` significant (`ARSC`, `Spermatheca`, `Germline`, `Embroynic_cells`), **Zebrafish 6-fish fine**: `0` significant (`min p = 0.02169`, `min q_within = 0.32535`).
   - Empirical-Bayes variance borrowing detects more shifts than unmoderated single-type CLR OLS in Fly and Worm, while confirming `0` FDR-adjusted shifts in the 6-fish Zebrafish fine cohort.

## 5. Reconciliation with Nov 2 Exploratory Files and Phase 1 Propeller Files

Running the `run_organ_composition.py` pipeline fresh explains every discrepancy across earlier files:

1. **AFCA (`results/extended-afca-composition-tests.csv`, Nov 2)**:
   - Nov 2 `all_samples` matched our fresh `primary_clr_0.5` exactly (`162/162` rows matched, `max |effect diff| = 0.00`, `max |q diff| = 0.00`). What Nov 2 omitted were pseudocount `0.1` and `1.0` checks, `propeller/limma` moderation, and global BH with the `431` existing tests. Applying an extra `total_cells >= 50` filter (used in Worm and Zebrafish Nov 2 scripts) reduces AFCA Head from `90` to `88` cell types (`19` hits at `0.5`) and Body from `72` to `63` cell types (`15` hits at `0.5`).
2. **Gao Worm (`results/extended-worm-composition-tests.csv`, Nov 2)**:
   - Nov 2 added `counts.sum(axis=0) >= 50`, which did not drop any of the `16` cell types passing `n_present >= 6`. Thus Nov 2 `all_annotated_cells` / `all_samples` (`16` tests, `2` significant: `ARSC`, `Spermatheca`) matches our fresh `primary_clr_0.5` identically (`max |effect diff| = 0.00`).
   - Why later Phase 1 files (`reviewer-phase1-composition-propeller.csv`, `reviewer-phase1-composition-total-adjusted.csv`) reported more Worm hits (`3` to `5`): they added `log_total_cells` and `platform_nextseq` (`HiSeq` vs `NextSeq` batch) to the baseline model, under which `Germline` (`platform_adjusted` `effect = -0.07005`, `q = 0.03345`), `Neuron`, and `Embroynic_cells` / `Gonadal_sheath_cell` became significant. Notice that `59,896 / 81,196` (`73.77%`) of worm nuclei are reproductive (`Germline`, `Embroynic_cells`, `Sperms`); restricting to `13` somatic cell types preserves `ARSC` and `Spermatheca` (`q_within = 1.06e-04`).
3. **Zebrafish Testis (`results/extended-zebrafish-composition-tests.csv`, Nov 2; vs `reviewer-phase1-*`)**:
   - In the 6-fish `5 to 22 month` cohort (`36,523` cells), all `15` retained fine cell types satisfy both `n_present >= 6` and `total_cells >= 50`, so Nov 2 `fine_5_to_22_months` (`all_samples`) matches our `primary_clr_0.5` identically (`15/15` matched, `0` significant).
   - **Critical discrepancy with `reviewer-phase1-*`**: [run_reviewer_phase1_composition.py](file:///home/reza/projects/agescape-revision/scripts/run_reviewer_phase1_composition.py) and `reviewer-phase1-composition-propeller.csv` dropped the 6-fish fine stratification and instead pooled the 6 `5 to 22 month` fish (`cell_type`) with 4 separate `27 month` fish (`cell.type.cat`, two of which had `incomplete_spermatogenesis`) into 3 broad compartments (`germ cells`, `immune`, `blood_and_somatic`, `N = 10`).
   - Within the unconfounded 6-fish cohort even at broad level (`broad_6fish_5_to_22_months`), `0` of `3` broad compartments pass BH (`immune` `p = 0.0266`, `q_within = 0.0799`). Only when the 4 separate 27-month fish are appended (`broad_10fish_all_ages_reconciliation`, `N = 10`) does `immune` pass BH (`+0.0926/month`, `p = 0.00678`, `q_within = 0.02033`; `spermatogenesis_adjusted` `p = 0.000725`, `q_within = 0.00217`; `total_cells_adjusted` also flags `germ cells` decline `q = 0.04449`).

## 6. Mandatory Zebrafish Testis Germline Caveat

**Caveat**: Zebrafish testis in `GSE275361` is overwhelmingly **germline**, not somatic organ parenchyma:
- In the primary 6-fish cohort (`5 to 22 months`, `36,523` total cells; `36,439` excluding the single-fish `liver` contaminant), **`35,511 / 36,523` (`97.23%`)** of all cells belong to the spermatogenic germline lineage (`spermatogonia 1 to 4`, `spermatocytes 1 to 3`, `early/middle/late round spermatids`, `elongating spermatids`, and `spermatozoa`), leaving only `928` (`360` `immune`, `336` `blood & SMC-like`, `232` `Sertoli & Leydig`) somatic/immune cells in the 15 eligible cell types.
- Even in the 10-fish combined broad dataset (`69,576` cells), **`51,755 / 69,576` (`74.39%`)** are `germ cells`.
- Germline maturation and age-related spermatogenic arrest reflect reproductive-tract biology and meiotic progression rather than somatic organ aging, and proportions within a 97%-germline tissue are strongly coupled by the constant-sum closure constraint. Therefore, zebrafish testis compositional findings (null in 6-fish fine, immune enrichment when 27-month spermatogenic-arrest fish are added) must not be interpreted as evidence for or against somatic organ aging.
