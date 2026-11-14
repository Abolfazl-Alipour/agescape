# Item 8 Report: Cross-Sectional Randomness and Variability Test in Mouse Tabula Muris Senis (`18m` vs `30m`)

Date: 2026-11-14
Script: [run_rev3_randomness_variability_test.py](file:///home/reza/projects/agescape-revision/scripts/run_rev3_randomness_variability_test.py)
Result tables:
- Per-stratum arm results (`148` rows across `5` pseudobulk arms and within-mouse LOO dispersion): [rev3-randomness-variability-test.csv](file:///home/reza/projects/agescape-revision/results/rev3-randomness-variability-test.csv)
- Coordinated `8`-mouse exact global permutation table (`binom(8, 4) = 70` partitions): [rev3-randomness-variability-joint-global.csv](file:///home/reza/projects/agescape-revision/results/rev3-randomness-variability-joint-global.csv)

---

## 1. Fundamental Study-Design Constraints: Cross-Sectional Comparison and `30m` Sex Confounding

**Important methodological and biological limitation**: Tabula Muris Senis (`TMS`) droplet data ([m1-qc-cleaned.h5ad](file:///home/reza/projects/agescape/data/interim/m1-qc-cleaned.h5ad)) is **strictly cross-sectional**, not longitudinal (`18m` and `30m` are distinct sacrificed mice; no pre/post trajectory exists for any individual animal), and the `18m` vs `30m` comparison is structurally unbalanced by sex:
1. **`18m` mice (`n = 4`) are sex-balanced (`2` females, `2` males)**: `18-F-50` (`female`, `11,922` raw singlet-filtered cells across tissues), `18-F-51` (`female`, `9,694`), `18-M-52` (`male`, `17,429`), and `18-M-53` (`male`, `4,852`).
2. **`30m` mice (`n = 4`) are 100% male (`0` females)**: `30-M-2` (`male`, `23,182`), `30-M-3` (`male`, `9,453`), `30-M-4` (`male`, `7,764`), and `30-M-5` (`male`, `12,395`).
3. **Consequence for variance (`s_30m^2` vs `s_18m^2`)**: Pooling `2` females and `2` males at `18m` injects **inter-sex variance** into `s_18m^2` even after removing `Xist`, `Tsix`, and Y-chromosome markers (`Ddx3y, Eif2s3y, Uty, Kdm5d`) and library depth, whereas `s_30m^2` reflects solely inter-male variance among `30m` survivors (survivorship bias). Therefore, mixed-sex `30m` vs `18m` comparisons can underestimate true within-sex aging variance or attribute autosomal sex differences to age (`eta^2`). We address this directly with a **`male_only` sensitivity arm** comparing the `2` males at `18m` against `3 to 4` males at `30m` (`20` strata in `Lung` and `Marrow`).

---

## 2. Study Cohort, Metrics, and Multiple-Testing Design

Starting from all non-doublet singlets with `n_counts >= 500` and known `cell_ontology_class` at `18m` and `30m`, `34` `(tissue, cell_ontology_class)` strata meet eligibility (`>= 10` cells per mouse in `>= 3` `18m` mice and `>= 3` `30m` mice): **`11` Kidney**, **`11` Lung**, and **`12` Marrow** (`39,646` cells total in the primary arm; `7,482` to `12,238` detected autosomal genes in `>= 50%` of donor pseudobulks after removing `mt-`, `Rp[ls]`, `Mrp[ls]`, and sex-linked transcripts).

To decouple systematic aging drift from stochastic transcriptional noise and cell-level heterogeneity, we compute six complementary metrics/arms at the **mouse level** (never cell-level p-values):
1. **`primary_10cells_full_residualized` (`34` strata)**: Depth-residualized `log1p(CP10K)` donor pseudobulks (`mean_log_umi` linearly regressed out per gene). Measures:
   - Mean age variance fraction (`mean_eta2 = mean_g [SS_between / SS_total]` across genes, measuring between-group shift vs within-group spread).
   - Finite-df digamma-debiased gene-wise log variance ratio (`mean_log_var_ratio_30_vs_18 = mean_g [log((var_30 + 1e-3) / (var_18 + 1e-3)) - (bias(n_30) - bias(n_18))]`), unweighted `log_total_var_ratio_30_vs_18 = log(mean(var_30) / mean(var_18))`, and `frac_genes_s30_gt_s18`.
   - Exact single-stratum permutation p-values over all `binom(n_d, n_30)` (`20` or `35`) animal-label partitions.
2. **`sensitivity_10cells_full_unresidualized` (`34` strata)**: Full-depth `log1p(CP10K)` without depth residualization.
3. **`sensitivity_10cells_rarefied_1500_pseudobulk` (`34` strata, `39,570` rarefied cells)**: Exact multivariate hypergeometric UMI subsampling to `1,500` filtered UMIs/cell prior to pseudobulking.
4. **Within-mouse cell-to-centroid Hellinger dispersion (`34` strata at `1,500` rarefied UMIs)**: Leave-one-out cell-to-donor-centroid Hellinger distance (`diff_30m_minus_18m_disp`), isolating within-animal cellular noise free of inter-mouse differences.
5. **`sensitivity_20cells_full_residualized` (`26` strata, `37,814` cells)**: Stricter `>= 20` cells/mouse cutoff.
6. **`sensitivity_10cells_male_only_residualized` (`20` strata across `Lung` `10` and `Marrow` `10`, `26,948` cells)**: Unconfounded comparison among male mice only (`n_18 = 2`, `n_30 = 3 to 4`).

**Explicit multiple-testing families**: Benjamini-Hochberg `q` values are computed within each analysis arm separately across eligible strata (`m = 34`, `26`, or `20`). Because a `3 vs 4` split only supports `35` permutations (`min one-sided p = 0.02857`, `min two-sided p = 0.05714`), we additionally test a **single global hypothesis free of cross-stratum multiplicity**: the **coordinated `8`-mouse exact permutation test** across all `binom(8, 4) = 70` global age partitions simultaneously ([rev3-randomness-variability-joint-global.csv](file:///home/reza/projects/agescape-revision/results/rev3-randomness-variability-joint-global.csv)).

---

## 3. Headline Results Across All Six Arms

| Analysis Arm | Strata (`Cells`) | Med `mean_eta2` (Min `q_bh`) | Mean `log_var_ratio` (Med) | Mean `log_total_var_ratio` (Med) | Mean `frac_s30_gt_s18` (Med) | Min `q_bh` LVR (`+` / `-` / `2s`) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `primary_10cells_full_residualized` | `34` (`39,646`) | `0.2005` (`0.2429`) | `-0.1233` (`-0.1786`) | `-0.0480` (`-0.1071`) | `0.4451` (`0.4623`) | `0.9714` / `0.3886` / `0.7286` |
| `sensitivity_10cells_full_unresidualized` | `34` (`39,646`) | `0.1917` (`0.2159`) | `-0.1810` (`-0.1405`) | `-0.2484` (`-0.3022`) | `0.4164` (`0.4133`) | `1.0000` / `0.2776` / `0.5299` |
| `sensitivity_10cells_rarefied_1500_pseudobulk` | `34` (`39,570`) | `0.1812` (`0.3238`) | `-0.2818` (`-0.2761`) | `-0.3216` (`-0.3404`) | `0.4088` (`0.4240`) | `0.9714` / `0.1943` / `0.3238` |
| `sensitivity_20cells_full_residualized` | `26` (`37,814`) | `0.1982` (`0.2476`) | `-0.0651` (`-0.0317`) | `+0.0410` (`+0.0945`) | `0.4664` (`0.4764`) | `0.9774` / `0.4333` / `0.7429` |
| `sensitivity_10cells_male_only_residualized` | `20` (`26,948`) | `0.2487` (`0.8000`) | `-0.5520` (`-0.5323`)* | `+0.2707` (`+0.3051`) | `0.6117` (`0.6291`) | `0.7143` / `0.6667` / `1.0000` |
| **Within-mouse `1500`-UMI LOO Hellinger disp** | `34` (`39,570`) | - | `diff = +0.0014` (`+0.0012`) | - | `19/34` `30m > 18m` | `min q_2s = 0.3886` (`0/34` hits) |

*Technical note on `male_only` `mean_log_var_ratio` (`-0.5520`) vs `log_total_var_ratio` (`+0.2707`) and `frac_s30_gt_s18` (`0.6117`): with only `n_18 = 2` males (`df = 1`), Gaussian log-variance debiasing subtracts `digamma(0.5) - log(0.5) = -1.2704` from `var_18`, which over-penalizes `log((var_30 + 1e-3) / (var_18 + 1e-3))` when near-zero gene variances are bounded below by `eps = 1e-3`. By the unregularized `log_total_var_ratio` (`+0.3051` median = `1.36x` total pseudobulk variance at `30m`) and raw gene share (`62.9%` median with `s_30m^2 > s_18m^2`), male-only `30m` animals actually trend higher in variance, whereas `n_18 = 2` limits permutation power (`binom(6, 4) = 15` permutations, `p >= 0.10` per stratum).*

### Coordinated `8`-Mouse Global Permutation Result (`binom(8, 4) = 70` Partitions)

| Global Arm (`n = 34` Strata) | Observed Mean `eta^2` | Joint Null `eta^2` (SD) | Exact `p` (`eta^2`) | Observed Mean `log_var_ratio` | Joint Null `lvr` (SD) | Exact `p` (`30m > 18m` / `2-sided`) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `primary_full_residualized_34` | `0.2107` | `0.1704` (`0.0155`) | **`0.02857`** (`2/70`) | `-0.1233` | `0.0000` (`0.1028`) | `0.8714` / `0.2857` |
| `rarefied_1500_pseudobulk` | `0.1867` | `0.1704` (`0.0091`) | `0.08571` (`6/70`) | `-0.2818` | `0.0000` (`0.2026`) | `0.8857` / `0.2571` |

---

## 4. Scientific Interpretation and Representative Strata

1. **No evidence of increased cross-sectional transcriptional randomness in `30m` vs `18m` TMS mice**:
   - **Inter-animal pseudobulk variance**: Across all `34` mixed-sex strata, `30m` animals show **lower** (not higher) average log variance than `18m` animals in both primary residualized (`mean_log_var_ratio = -0.1233`, `joint p (30m > 18m) = 0.8714`, two-sided `0.2857`) and rarefied `1,500`-UMI pseudobulk (`-0.2818`, `joint p = 0.8857` / `0.2571`), with `0 / 34` strata surviving BH (`q >= 0.7286` two-sided).
   - **Within-animal cellular dispersion (`1,500` UMI LOO Hellinger)**: Average `30m - 18m` shift in within-mouse cell-to-centroid dispersion is virtually zero (`+0.0014` mean, `+0.0012` median; `19/34` positive; `0 / 34` BH hits, `min q_2s = 0.3886`).
2. **Modest coordinated between-age mean divergence in full-depth residualized data weakens upon UMI rarefaction**:
   - Globally across the `34` strata, `mean_eta2 = 0.2107` beats all but `1` of the `70` animal-label swaps in unrarefied residualized data (`joint p = 0.02857`), led by `Marrow` strata (`granulocytopoietic cell` `mean_eta2 = 0.3855`, `granulocyte = 0.3378`, `promonocyte = 0.2879`, `macrophage = 0.2808`, `naive T cell = 0.2778`, `megakaryocyte-erythroid progenitor cell = 0.2704`). However, exact `1,500`-UMI rarefaction lowers the global mean `mean_eta2` from `0.2107` to `0.1867` (`joint p = 0.08571`), indicating sensitivity to sequencing-depth heterogeneity across `18m` and `30m` animals (as well as female/male differences mixed into `18m` vs `30m`).
3. **Why mixed-sex `18m` vs `30m` variance masks potential within-sex patterns**:
   - Dropping `18m` females in the `20`-stratum `male_only` sensitivity reverses the fraction of genes with `s_30m^2 > s_18m^2` from `44.5%` (mixed-sex) to **`61.2%` mean / `62.9%` median** and flips `log_total_var_ratio` from `-0.0480` to **`+0.2707` mean / `+0.3051` median (`1.36x` total variance)**. However, having only `n = 2` males at `18m` prevents any stratum or test from achieving significance (`min p = 0.10`, `min q_eta2 = 0.80`).
4. **Tissue specificity across primary strata (`34` strata)**:
   - Top `eta^2` strata (`p_eta2_exact = 0.02857`, `q_bh = 0.2429`): `Marrow promonocyte` (`eta^2 = 0.2879`, `log_tot_var_ratio = +0.9049`), `Marrow naive T cell` (`0.2778`, `+0.5697`), `Marrow megakaryocyte-erythroid progenitor cell` (`0.2704`, `+0.2883`). Largest positive `log_total_var_ratio`: `Marrow monocyte` (`+0.9589`, `frac_s30_gt_s18 = 0.643`), `Marrow promonocyte` (`+0.9049`), `Kidney proximal convoluted tubule epithelial cell` (`+0.5344`, `frac = 0.671`), `Lung B cell` (`+0.5049`, `frac = 0.578`). Largest negative ratios concentrate in `Marrow proerythroblast` (`-1.1616`), `Lung vein endothelial cell` (`-0.6806`), and `Kidney T cell` (`-0.5700`) / `proximal tubule` (`-0.4030`).
   - Rarefied within-mouse LOO dispersion extremes (`p_two = 0.05714`, `q_two = 0.3886`): higher in `30m` for `Marrow granulocyte` (`0.4696 -> 0.5155`, `+0.0459`), `Kidney distal convoluted tubule` (`+0.0122`), and `Kidney collecting duct principal cell` (`+0.0066`); lower in `30m` for `Marrow hematopoietic precursor cell` (`0.6507 -> 0.6295`, `-0.0212`) and `Kidney T cell` (`-0.0144`).

**Conclusion**: Cross-sectional `18m` vs `30m` mouse TMS data do **not** establish a general age-related increase in inter-individual transcriptional variance (`joint p = 0.2857` two-sided) or within-animal cellular Hellinger dispersion (`mean diff = +0.0014`, `0/34` BH hits). Any suggestive within-sex male variance elevation (`62.9%` of genes, `+0.305` median `log_total_var_ratio`) is confounded and severely underpowered (`n = 2` `18m` males vs `n = 4` `30m` males) and cannot separate stochastic aging randomness from survivorship or autosomal sex-dimorphism effects.
