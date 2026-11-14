# Item 4 Report: Mature-Only (`>= 18 Months`) Transcriptional Entropy in Mouse Lung

- **Date**: 2026-11-14
- **Script**: [run_rev3_entropy_mature_mouse_lung.py](file:///home/reza/projects/agescape-revision/scripts/run_rev3_entropy_mature_mouse_lung.py)
- **Primary result tables**:
  - [rev3-entropy-mature-only-mouse-lung.csv](file:///home/reza/projects/agescape-revision/results/rev3-entropy-mature-only-mouse-lung.csv) (`1,593` rows across depths `500`, `1500`, `3000`, `full`; 3 entropy definitions; 3 covariate specifications; 4 eligibility families)
  - [rev3-entropy-mature-only-mouse-lung-donors.csv](file:///home/reza/projects/agescape-revision/results/rev3-entropy-mature-only-mouse-lung-donors.csv) (`449` donor-by-cell-type-by-depth rows)

---

## 1. Question and Mature Lung Cohort Design

Does within-cell transcriptional entropy increase from maturity to old age in mouse lung when young (`1m` and `3m`) mice are excluded? And do the De Man et al. 2026 (*Nature Communications*, `s41467-026-68810-9`) *Methods* KL-divergence formula and public GitHub `sqrt`-count variant agree with standard Shannon entropy across depth-matched UMI rarefactions?

- **Input**: Singlet `Lung` cells from `/home/reza/projects/agescape/data/interim/m1-qc-cleaned.h5ad` with `age_months >= 18` and non-`unknown` `cell_ontology_class` (`18,379` cells total).
- **Mature mice (`n = 10`, `24m` has 0 lung mice)**:
  - `18m` (`n = 4`): `18-F-50` (`328` cells), `18-F-51` (`1,154`), `18-M-52` (`1,596`), `18-M-53` (`1,278`)
  - `21m` (`n = 2`, **female only**): `21-F-54` (`1,364`), `21-F-55` (`1,578`)
  - `30m` (`n = 4`, **male only**): `30-M-2` (`8,014` cells = `43.6%` of mature lung cells, shallowest depth), `30-M-3` (`1,629`), `30-M-4` (`934`), `30-M-5` (`504`)
- **Cell-type eligibility (`>= 20` cells per mouse)**: `15` cell types have `>= 3` mature mice; **`11`** meet the `n >= 6` mouse regression floor (`all_mature_eligible_n6`), including **`8`** of the `9` original manuscript targets (`manuscript_target_n6`). `fibroblast of lung` has only **`n = 5`** mature mice (`manuscript_target_include_fibroblast_n5`). Because `21m` is 100% female and `30m` is 100% male, `age_only` models are reported alongside `age_sex` and `age_sex_library` (`mean_log_umi_actual`) for sensitivity against sex-by-age confounding.
- **Metrics compared** (positive `shannon_entropy` slope = higher disorder; **negative** `deman_kl_*` slope = closer to uniform = higher disorder):
  1. `shannon_entropy`: `H(p) = - sum p_g ln(p_g)` after removing `mt-`, `Rpl/Rps`, `Mrpl/Mrps`, and sex-linked (`Xist, Tsix, Ddx3y, Eif2s3y, Uty, Kdm5d`) genes.
  2. `deman_kl_methods`: Exact De Man et al. 2026 *Methods* formula `D_KL(p || u) = ln(G) - H(p)` on linear counts over `G` nonzero genes.
  3. `deman_kl_sqrt`: De Man public GitHub `sc.ent.KL.count` variant (`sqrt(count)` before normalizing to `p_g`).
- **Multiple testing family and MDE**: Benjamini-Hochberg `q` values are computed across eligible cell types within each `(family, depth, metric, model)` slice (`m = 11` for `all_mature_eligible_n6` at `500/1500/full`, `m = 9` at `3000`; `m = 8` for `manuscript_target_n6`). Minimal 80%-power detectable effect is `MDE_80 = (t_{0.975, df} + t_{0.80, df}) * SE`.

---

## 2. Headline Global Results (`n = 10` Mature Mice, Donor-Clustered OLS)

Across all rarefied depths (`500`, `1500`, `3000` UMIs) and unrarefied `full` depth, **no global mature-only shift in transcriptional disorder reaches statistical significance under `age_sex` or `age_only`**, and 80%-power MDEs span **`1.41 to 2.27` standard deviations per 15 months** (`1.13 to 1.82` SD across the observed `12`-month `18m -> 30m` interval), confirming low statistical power rather than evidence of absence.

| Family (`n` types) | Depth | Metric | Model | `Effect_15mo_SD` | `SE` | `MDE80_15mo_SD` | `pvalue` |
|---|---:|---|---|---:|---:|---:|---:|
| `all_mature_eligible_n6` (11) | `500` | `shannon_entropy` | `age_sex` | `+0.2800` | `0.5058` | `1.5910` | `0.5934` |
| `all_mature_eligible_n6` (11) | `1500` | `shannon_entropy` | `age_sex` | `+0.1604` | `0.5425` | `1.7064` | `0.7742` |
| `all_mature_eligible_n6` (9) | `3000` | `shannon_entropy` | `age_sex` | `+0.0508` | `0.5164` | `1.6245` | `0.9239` |
| `all_mature_eligible_n6` (11) | `full` | `shannon_entropy` | `age_sex` | `-0.1318` | `0.5473` | `1.7217` | `0.8151` |
| `manuscript_target_n6` (8) | `1500` | `shannon_entropy` | `age_sex` | `+0.1380` | `0.4757` | `1.4963` | `0.7783` |
| `manuscript_target_n6` (8) | `1500` | `deman_kl_methods` | `age_sex` | `-0.1613` | `0.4733` | `1.4889` | `0.7411` |
| `manuscript_target_n6` (8) | `1500` | `deman_kl_sqrt` | `age_sex` | `-0.6464` | `0.6288` | `1.9778` | `0.3308` |
| `manuscript_target_n6` (8) | `1500` | `shannon_entropy` | `age_only` | `+0.0356` | `0.4831` | `1.5197` | `0.9429` |
| `manuscript_target_n6` (8) | `full` | `shannon_entropy` | `age_sex_library` | `+0.2969` | `0.2968` | `0.9336` | `0.3433` |
| `manuscript_target_n6` (8) | `full` | `deman_kl_sqrt` | `age_sex_library` | `-0.9894` | `0.3408` | `1.0721` | `0.0175` |

Only one global model across the entire grid (`manuscript_target_n6`, `full` depth, `deman_kl_sqrt`, `age_sex_library`) reaches nominal `p = 0.0175` (`-0.9894 SD`), and it disappears when rarefied (`1500` `p = 0.0825`, `3000` `p = 0.1060`), when `sqrt` is omitted (`deman_kl_methods` `p = 0.1221`, `shannon_entropy` `p = 0.3433`), or when `sex` is omitted (`deman_kl_sqrt age_only` `p = 0.6898`).

---

## 3. Per-Cell-Type Results (`1500`-UMI Primary Depth) and Sex / Depth Confounding

At rarefied `1500` UMIs in `all_mature_eligible_n6` (`age_sex` model, `m = 11`):

| Cell Type | `n` Mice | `shannon` Slope/mo (`SE`) | `12mo` MDE80 | `p` | `q_bh` | `age_only p` (`q_bh`) |
|---|---:|---:|---:|---:|---:|---:|
| `CD8-positive, alpha-beta T cell` | `8` | `+0.01724` (`0.00310`) | `0.1300` | `0.00260` | `0.0286` | `0.0665` (`0.2779`) |
| `bronchial smooth muscle cell` | `9` | `-0.00759` (`0.00210`) | `0.0844` | `0.01112` | `0.0611` | `0.0216` (`0.2380`) |
| `B cell` | `9` | `+0.01041` (`0.00337`) | `0.1357` | `0.02148` | `0.0788` | `0.2282` (`deman_kl`) |
| `vein endothelial cell` | `7` | `-0.00818` (`0.00293`) | `0.1307` | `0.04916` | `0.1352` | `0.6600` (`full`) |
| `alveolar macrophage` | `9` | `-0.00014` (`0.00629`) | `0.2532` | `0.9830` | `0.9979` | - |
| `classical monocyte` | `10` | `-0.00001` (`0.00360`) | `0.1409` | `0.9979` | `0.9979` | `0.6281` (`0.7568`) |
| `fibroblast of lung` (`n>=5` arm) | `5` | `+0.02384` (`0.01941`) | `1.2493` | `0.3444` | `0.6198` (`m=9`) | - |

Two structural artifacts explain why per-cell-type `q < 0.05` hits are unstable across specifications:
1. **`CD8-positive, alpha-beta T cell` (`+0.0172/mo`, `age_sex q = 0.0286` at `1500`; `0.0208` in `m = 8`)** goes in the **opposite** direction of `bronchial smooth muscle cell` (`-0.0076/mo`, `q = 0.0445` in `m = 8`) and depends directly on sex collinearity: at `18m`, male mice (`18-M-52/53`) have lower Shannon entropy (`5.7081` and `5.6926`) than females (`18-F-51` `5.8700`, `21-F-54/55` `5.8889` and `5.8668`) and `30m` males (`5.8665` to `5.9793`). With `21m` 100% female and `30m` 100% male, fitting `sex_male` absorbs `-0.18` baseline male difference into an apparent upward slope; in `age_only`, **`CD8+` has zero `q < 0.05` hits at any depth (`q = 0.087 to 0.278`)**.
2. **`age_sex_library` post-rarefaction collider effect in `classical monocyte`**: After equalizing all cells to `1,500` UMIs, `classical monocyte` shows zero age slope under `age_sex` (`p = 0.9979`) or `age_only` (`p = 0.6281`), yet flips to `deman_kl_sqrt` `q = 0.0207` (`p = 0.00188`) when adjusting for raw `mean_log_umi_actual`. Inspecting donor rows explains why: mouse `30-M-2` contributes `5,383` classical monocytes with `mean_umi_actual = 3,911` (less than half the `8,046 to 11,014` UMI depth of the other 9 mice) AND lower rarefied Shannon entropy (`5.907` vs `5.99 to 6.07`). Conditioning on pre-rarefaction UMI depth when all `30m` mice are male induces leverage on `30-M-2`.

---

## 4. Direct Formula Comparison (Shannon vs De Man Methods `D_KL` vs De Man GitHub `sqrt`)

Across `all_mature_eligible_n6` (`age_sex` and `age_only`), `shannon_entropy` and `deman_kl_methods` (`sqrt=False`) behave as expected mirror images (`D_KL = ln(G) - H(p)` plus retained MT/RP counts), whereas `deman_kl_sqrt` shifts weighting toward low-expressed genes by compressing high-count transcripts via `sqrt(count)`, amplifying sensitivity to cell-depth variance (turning `classical monocyte age_sex_library` at `1500` UMIs and `manuscript_target_n6 full age_sex_library` nominally stronger while weakening `bronchial smooth muscle cell`).

**Conclusion**: Restricting mouse lung strictly to `>= 18m` mice does NOT reveal a general within-cell disorder increase under donor-level inference. Both Shannon entropy and De Man's KL divergence remain statistically inconclusive globally (`MDE_80 ~ 1.5 to 1.7 SD` per `15m`) due to `n = 10` mature animals and age-by-sex confounding across `18m`, `21m`, and `30m`.
