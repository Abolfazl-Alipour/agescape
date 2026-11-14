# Agescape Overnight Revision Summary (2026-11-14)

- **Repository**: `/home/reza/projects/agescape-revision` (`master`, base `origin/master` = `4009829`)
- **Starting `HEAD`**: `4009829` -> **Last Item Commit (`Item 8`)**: `7fb9488` -> **Morning Summary**: `27a61fc` (followed by CI/hygiene/`PRJNA1058087` update on `master`)
- **Push Status**: **PUSHED to `origin/master`** (`27a61fc` pushed earlier, plus CI/hygiene/PRJNA1058087 follow-up commit synchronized cleanly on `master`)
- **Frozen Files Status**: Verified untouched (`git diff --exit-code -- report/plos-one/ report/manuscript-aging-mechanisms.md` exited `0`)

---

## Step 0. Workspace & To-Do Tracker Lock (`DONE`)

- **Commit**: `dae5a98` (`Item 0: lock to-do list tracker`)
- **Tracker**: [todo-progress.md](file:///home/reza/projects/agescape-revision/notes/todo-progress.md)
- **Summary**: Verified working tree at `/home/reza/projects/agescape-revision`, data root `/home/reza/projects/agescape-revision-data`, Python env `/home/reza/miniconda3/envs/agescape/bin/python`, read `notes/agent-handoff-2026-11-14.md`, `notes/followup-execution-plan.md`, and frozen [final-plos-one.tex](file:///home/reza/projects/agescape-revision/report/plos-one/final-plos-one.tex), and locked the 9-item revision tracker.

---

## Item 1. CLR + Limma/Propeller Composition & Global 624-Test BH (`DONE`)

- **Commit**: `0485256`
- **Report**: [rev3-composition-fly-worm-zebrafish-REPORT.md](file:///home/reza/projects/agescape-revision/results/rev3-composition-fly-worm-zebrafish-REPORT.md)
- **Tables**: [rev3-composition-fly-worm-zebrafish.csv](file:///home/reza/projects/agescape-revision/results/rev3-composition-fly-worm-zebrafish.csv), [rev3-composition-global-bh-624.csv](file:///home/reza/projects/agescape-revision/results/rev3-composition-global-bh-624.csv)
- **Headline**: Across `193` new primary CLR composition tests (`90` Fly Head `36` pools: `19` within-BH / `12` global-624 hits; `72` Fly Body `36` pools: `17` within / `14` global; `16` C. elegans N2 `12` pools: `2` within/global hits `ARSC +0.2160/d q_624=0.00119` and `Spermatheca -0.1366/d q_624=0.00175`; `15` Zebrafish `6`-fish fine testis: `0` hits, `min p=0.0715 q_within=0.9393`, `97.23%` germline `35,511/36,523`), `45/624` total global BH tests reach `q < 0.05` (`28` new + `17` existing mouse/human/NMR).

## Item 2. LACA Naked Mole-Rat Multi-Tissue Atlas Verification (`DONE`)

- **Commit**: `33da505`
- **Report**: [rev3-laca-nmr-verification-REPORT.md](file:///home/reza/projects/agescape-revision/results/rev3-laca-nmr-verification-REPORT.md)
- **Headline**: Live inspection established `LACA-2024-002` (`1.2M` cells) and `"NMR Full Atlas v2.1"` (`660K` cells) have `0` downloadable files on `download.php`, `SP239` endpoints return `403 Forbidden`, site JS contains developer comment `SP239 unpublished self-test` redirecting to human `data/SP121/umap`, `0` public multi-tissue NMR aging atlases exist on GEO/BioProject/PubMed/CELLxGENE, and direct NCBI/CNCB/PMC inspection of `PRJNA1058087` (`PRJCA014972` / `PMC10652311`) confirmed it only generated new bat (`Myotis pilosus`) spleen scRNA-seq while reusing Hilton et al. 2019 (`GSE132642`) for NMR (zero new or aged NMR single-cell data); no download was attempted.

## Item 3. Audit of Young-vs-Mature Comparisons in Manuscript (`DONE`, Audit Only)

- **Commit**: `b68077e`
- **Report**: [rev3-young-vs-mature-audit-REPORT.md](file:///home/reza/projects/agescape-revision/results/rev3-young-vs-mature-audit-REPORT.md)
- **Headline**: Audited all `11` manuscript sections against result tables with zero edits to frozen `report/plos-one/*`: verified `P1_Clean` is `100%` `3m vs 18m` (`17` lung/kidney units), `mouse_ref_18m` in `reviewer-rev2-hellinger-donors.csv` still included all `23` `1m to 30m` mice, `9/9` testable mouse mature-only CLR losses keep sign (`median 0.89x`) under `3.28x` SE inflation while `5/5` human mature shifts hold (`q=0.0118 to 0.0279`), Worm drift collapses from Day 1 (`+1.5956 p=4.87e-5`) to Day 6 (`+0.3853 p=0.1001`), and no `>= 18m` mouse lung entropy model previously existed.

## Item 4. Mature-Only (`>= 18m`) Shannon Entropy & De Man KL in Mouse Lung (`DONE`)

- **Commit**: `4d28640`
- **Report**: [rev3-entropy-mature-only-mouse-lung-REPORT.md](file:///home/reza/projects/agescape-revision/results/rev3-entropy-mature-only-mouse-lung-REPORT.md)
- **Tables**: [rev3-entropy-mature-only-mouse-lung.csv](file:///home/reza/projects/agescape-revision/results/rev3-entropy-mature-only-mouse-lung.csv) (`1,593` rows), [rev3-entropy-mature-only-mouse-lung-donors.csv](file:///home/reza/projects/agescape-revision/results/rev3-entropy-mature-only-mouse-lung-donors.csv) (`449` rows)
- **Headline**: Across `n=10` mature TMS lung mice (`18m` 2F+2M, `21m` 2F only, `30m` 4M only), `0/96` global donor-clustered models under `age_sex` or `age_only` and `0/108` rarefied models reach `p < 0.05` (`1500` UMI `manuscript_target_n6` `shannon +0.1380 SD p=0.7783 MDE_80=1.4963`; `deman_kl_methods -0.1613 p=0.7411`; `deman_kl_sqrt -0.6464 p=0.3308`; only unrarefied `full`-depth `deman_kl_sqrt` `age_sex_library` hits nominal `p=0.0150 to 0.0270` across the `4` eligibility families), and structural age-sex collinearity explains why per-cell-type `CD8+ T` `shannon_entropy` nominal `age_sex` hits at `1500` UMI (`all_mature_eligible_n6 q=0.0286`, `manuscript_target_n6 q=0.0208`) disappear under `age_only` (`q=0.2779` and `0.2021`).

## Item 5. De Man Human Lung Replication on `GSE281219` (`BLOCKED`)

- **Commit**: `aa9ac6f`
- **Report**: [rev3-deman-replication-human-lung-REPORT.md](file:///home/reza/projects/agescape-revision/results/rev3-deman-replication-human-lung-REPORT.md)
- **Headline**: `BLOCKED` for the `32` Baylor control lungs (`GSE281219`): raw `1.59 GB` `GSE281219_RAW.tar` (`1,589,473,280` bytes) is staged locally, but GEO, CELLxGENE, and `rubendeman/Lung-Aging-Cell-Atlas` (`add_meta.R`/`entropy.R`) do not publicly release `agingseurat.RData` or barcode-level `predicted.id` cell-type annotations; the `28` IPF control lungs (`GSE136831`) remain completed with `0` `q < 0.05` hits (`min q=0.275` gCap `published_KL_binomial`).

## Item 6. Mature-Only Bayesian Whole-Vector Drift (`DONE`)

- **Commit**: `e56c797`
- **Report**: [cross-species-bayesian-vector-mature-only.md](file:///home/reza/projects/agescape-revision/report/cross-species-bayesian-vector-mature-only.md)
- **Tables**: [rev3-bayesian-vector-mature-only-posterior.csv](file:///home/reza/projects/agescape-revision/results/rev3-bayesian-vector-mature-only-posterior.csv), [rev3-bayesian-vector-mature-only-diagnostics.csv](file:///home/reza/projects/agescape-revision/results/rev3-bayesian-vector-mature-only-diagnostics.csv), [rev3-bayesian-vector-mature-only-prior-predictive.csv](file:///home/reza/projects/agescape-revision/results/rev3-bayesian-vector-mature-only-prior-predictive.csv), [rev3-bayesian-vector-mature-only-frequentist.csv](file:///home/reza/projects/agescape-revision/results/rev3-bayesian-vector-mature-only-frequentist.csv)
- **Headline**: `10/10` mature-only hierarchical PyMC models converged (`0` divergences, `max_rhat <= 1.0070`); depth adjustment improves in-sample LOO by `+229.8` ELPD (caveat: `8 to 19` high Pareto-k observations per model) and shrinks primary `mature_ge_18m` (`14` mice, `51` units, `480` obs) drift slope to `+0.0389` (`94% HDI [-0.0998, +0.1773], P(>0)=0.7068`, vs OLS `-0.0137 p=0.8163`), `pseudobulk +0.0765 ([-0.2251, +0.3813])`, and unconfounded `male_ge_24m +0.2580 ([-0.3078, +0.8059])`, whereas excluding `18m` leaves sex-imbalanced `post_ref_gt_18m` (`n=10`, `21m` 2F vs `24m/30m` M) at `+0.2300 ([+0.0509, +0.4072], P=0.9819)`.

## Item 7. Paper-Wide Bayesian Synthesis Across Five Families (`DONE`)

- **Commit**: `532b8b3`
- **Report**: [cross-species-bayesian-full.md](file:///home/reza/projects/agescape-revision/report/cross-species-bayesian-full.md)
- **Tables**: [rev3-bayesian-full-summary.csv](file:///home/reza/projects/agescape-revision/results/rev3-bayesian-full-summary.csv), [rev3-bayesian-full-diagnostics.csv](file:///home/reza/projects/agescape-revision/results/rev3-bayesian-full-diagnostics.csv), [rev3-bayesian-full-prior-predictive.csv](file:///home/reza/projects/agescape-revision/results/rev3-bayesian-full-prior-predictive.csv)
- **Headline**: `50/50` paper-wide PyMC models across `5` families converged (`0` divergences, `max_rhat <= 1.0079`, `min_ess_bulk >= 725.5`, `29/29` models with `prior_sd >= 0.50` and `19/21` `prior_0.25` models within 99.8% prior predictive interval; `42/50` have `1 to 276` high-Pareto-k unit-leverage obs): confirms net `mu_species_std ~ 0` in CLR composition (`620` finite-SE rows of `624`, excluding `4` NMR `se=NaN` rank tests), zero transport (`mu_global_std = -0.0278 [-0.0874, +0.0364]`, random-null `z_excess = -0.0855 [-0.6999, +0.5928]`), species-specific concordance in Fly (`+0.1247, P=0.9947`) and Worm (`+0.1198, P=0.9721`) vs null Human (`+0.0188, P=0.7380`), null lung disorder (`20` models), and EIV interneuron slopes (`VIP 0.9109, pvalb 0.9272, sst 1.3479; P=1.0000`) plus `6` moderated shared genes (`+0.1156 [+0.0802, +0.1494]`, with 7th unmoderated `DYNC1I1` documented).

## Item 8. Cross-Sectional Randomness & Variability Test in Mouse TMS (`DONE`)

- **Commit**: `7fb9488`
- **Report**: [rev3-randomness-variability-test-REPORT.md](file:///home/reza/projects/agescape-revision/results/rev3-randomness-variability-test-REPORT.md)
- **Tables**: [rev3-randomness-variability-test.csv](file:///home/reza/projects/agescape-revision/results/rev3-randomness-variability-test.csv) (`148` rows), [rev3-randomness-variability-joint-global.csv](file:///home/reza/projects/agescape-revision/results/rev3-randomness-variability-joint-global.csv) (`2` rows, `70` exact perms)
- **Headline**: Across `34` primary TMS `18m` (2F+2M) vs `30m` (4M only) strata (`Kidney 11, Lung 11, Marrow 12`; `39,646` singlets), `0/34` strata survive BH (`min q_eta2=0.2429`, `min q_lvr_two=0.7286`, `min q_disp_two=0.3886`) and `70`-perm global `mean_log_var_ratio = -0.1233` (`p(30m > 18m)=0.8714`, 2-sided `0.2857`), within-mouse `1500`-UMI LOO Hellinger dispersion is near zero (`+0.0014`), global `mean_eta2 = 0.2107` (`p=0.02857`) drops after `1,500`-UMI rarefaction (`0.1867`, `p=0.08571`), and `male_only` (`20` strata, `2` vs `3 to 4` males) shifts variance upward (`62.9%` median genes `s_30m^2 > s_18m^2`, `+0.3051` median `log_total_var_ratio`) without power (`min q=0.6667 to 1.0`).

## Item 9. Zebrafish Testis Germline Caveat (`DO_NOT_DO`, Frozen)

- **Status**: Intentionally untouched per instructions (left frozen pending Reza confirmation).

---

## Blocked Items Summary

1. **Item 5 (`GSE281219` 32 Baylor Human Lungs)**: Public GEO/CELLxGENE/`rubendeman/Lung-Aging-Cell-Atlas` repositories do not contain `agingseurat.RData` or barcode-level `predicted.id` annotations; only raw matrices (`GSE281219_RAW.tar`) are public.
2. **Remote Push (`git push origin master`)**: Resolved; `master` is synchronized with `origin/master`.
3. **Item 9 (Manuscript Revisions) & Frozen Files**: Held untouched by explicit instruction.

---

## System, Test Suite & Crucible Status

- **New Revision Tests (`pytest -q tests/test_rev3_todo_items.py`)**: `8 passed`
- **Full Test Suite (`pytest -q`)**: `163 passed, 1 failed, 1 warning`
  - **Pre-existing failure**: `tests/test_human_vector_investigation.py::test_bayesian_hierarchy_maps_groups_to_units_before_observations` fails with `ModuleNotFoundError: No module named 'arviz_plots'` inside `scripts/run_h3_bayesian_disorder.py:11` (previously skipped when `pymc` was absent; triggered once `pymc` became available without installing `arviz_plots`). Left untouched outside scope and reported transparently.
  - **Pre-existing warning**: `tests/test_followup_muscle.py::test_kedlian_pipeline_synthetic_end_to_end` emits `ImplicitModificationWarning: Transforming to str index`.
- **Disk & Memory Status (`df -h /`, `free -h`)**: `190G` available (`67%` used of `598G` total on `/`, well above the `30 GB` reserve); `39Gi` available RAM of `41Gi` total (`1.7Gi` used).

---

## 10-Round Scientific Review Crucible (Executed After Completion)

| Round | Perspective | Brutally Honest Finding & Resolution |
| --- | --- | --- |
| 1 | Statistical Rigor (Exact Counts & BH Families) | Audited all CSV outputs directly against Item 1, 4, 6, 7, 8 claims. Caught why Item 7 models `620` of `624` composition tests (`4` NMR rank rows lacked finite SEs), disambiguated Item 4 `CD8+ T` `q=0.0286` (`all_mature_eligible_n6`) vs `q=0.0208` (`manuscript_target_n6`), and earlier fixed Item 8 report line 67 where Marrow `log_total_var_ratio` (`+0.9589`) had been phrased beside `mean_eta2 = 0.3855`. |
| 2 | Statistical Assumptions & Estimator Artifacts | Verified Item 1 CLR sum-to-zero dependence caveat alongside Smyth `propeller` comparison; confirmed Item 8 `male_only` negative debiased `mean_log_var_ratio` (`-0.5520`) vs positive `log_total_var_ratio` (`+0.2707`, `62.9%` genes) arises from `df=1` digamma subtraction (`-1.2704`) on `eps=1e-3`-floored genes, bounded by exact permutation limits (`binom(8,4)=70`, `binom(7,4)=35`, `binom(6,4)=15`). |
| 3 | Single-Cell Biology (Sex & Depth Confounders) | Audited mouse metadata in `rev3-entropy-mature-only-mouse-lung-donors.csv` and TMS `m1-qc-cleaned.h5ad`: confirmed `30m` is `100%` male (`30-M-2..5`) vs sex-balanced `18m` (`18-F-50, 18-F-51, 18-M-52, 18-M-53`) and `21m` females (`21-F-54, 21-F-55`), and depth rarefaction weakens Item 8 `eta^2` (`p=0.02857 -> 0.08571`). |
| 4 | Biological Scope & Upstream Data Absence | Confirmed zebrafish `6`-fish testis is `97.23%` (`35,511/36,523`) germline; verified Item 2 (`LACA SP239` unpublished placeholder) and Item 5 (`GSE281219` Baylor `predicted.id` / `agingseurat.RData` unavailable across GEO, CELLxGENE, and `rubendeman/Lung-Aging-Cell-Atlas`) are real public absences. |
| 5 | Bayesian Convergence & Pareto-k Honesty | Audited all `10` Item 6 and `50` Item 7 diagnostics: `60/60` models converged (`0` div, `max_rhat <= 1.0079`). Caught high Pareto-k warnings (`pareto_k_n_bad = 8 to 19` in all `10` Item 6 models; `1 to 276` in `42/50` Item 7 models) and added this caveat to ELPD claims; verified `29/29` models with `prior_sd >= 0.50` pass 99.8% prior predictive checks (only `2/21` narrow `prior_0.25` specs sit outside extreme Z tails). |
| 6 | Bayesian Contrasts & Conjugate EIV Validity | Audited all `10` Item 6 species slopes: surfaced positive `post_ref_gt_18m` (`+0.2300, [+0.0509, +0.4072], P=0.9819` on `10` sex-imbalanced post-reference mice) so it is not hidden, and confirmed Item 7 conjugate EIV + `interneuron_6shared` explicitly note excluded unmoderated 7th gene `DYNC1I1`. |
| 7 | Manuscript Consistency & Frozen Compliance | Verified zero edits to frozen `report/plos-one/*` and `report/manuscript-aging-mechanisms.md` (`git diff` exit `0`), kept Item 9 `DO_NOT_DO`, and verified strict separation of `P1_Clean` (139-gene signature) vs `P1` (PC1 score) with zero contradiction of verified core results `(a)` to `(f)`. |
| 8 | Cross-Item Numerical Reconciliation | Reconciled sample hierarchies across Items 3, 4, 6, 7, 8 (`n=8` primary `18m/30m` TMS mice + `n=2` `21m` females = `n=10` lung entropy cohort + `24m` mice = `n=14` mature-only vector cohort) and Item 1 (`624` BH tests) vs Item 7 (`620` finite-SE effects). |
| 9 | Reproducibility & Honest Test Reporting | Verified `8/8` new `tests/test_rev3_todo_items.py` pass and documented the `1` pre-existing full-suite failure (`test_bayesian_hierarchy_maps_groups_to_units_before_observations`, failing on missing `arviz_plots` once `pymc` became available; `163 passed, 1 failed, 1 warning`). |
| 10 | Provenance, Hyperlinks, Dash & Resource Audit | Verified all `23` clickable `file://` links resolve on disk, zero `\u2014` or `\u2013` characters exist across `17` inspected deliverables/scripts, local commits `dae5a98` through `7fb9488` match `git log`, `190G` disk free (`> 30 GB`), and `39Gi` available RAM. |
