# Agescape final remediation preview

## 1. Commits and working tree

Starting commit: `3d02af8e20c2a71aedd882614f56c5b730eba316`.

Ending commit: recorded in the final response after the report commit.

Working-tree status: recorded in the final response after the report commit.

No remote push was performed.

## 2. Task-by-task changes

### Task 1: Human composition rebuild

Changed `scripts/run_m9_pipeline.py`, `scripts/run_h2_composition.py`, and `scripts/run_d3_human_composition_sensitivity.py`.

The public Census extraction now writes complete donor-by-cell-type counts before the expression-score filter. The human analysis retains counts from 1 to 19 cells and distinguishes true zero observations from rows removed by a score threshold. Primary human composition changed from 117 donors and 25 cell types to 119 donors and 44 cell types.

Added complete-count, observed-record-only, and presence-versus-age outputs:

- `results/m9-donor-cell-counts-prethreshold.csv`
- `results/h2-composition-human-lung.csv`
- `results/h2-composition-human-lung-observed-only.csv`
- `results/h2-composition-human-lung-presence.csv`
- `results/h2-composition-cross-species-test.csv`

Human FDR-significant compositional types changed from four to two. Club cells changed from slope -0.067799 and q = 0.0123 to -0.060227 and q = 0.004423. Multiciliated epithelial cells entered at slope -0.034154 and q = 0.022513. Club-cell study adjustment changed from p = 0.2415 to p = 0.103282. Cross-species Spearman correlation changed from -0.286 with p = 0.535 to 0.000 with p = 1.000.

### Task 2: Abstract overclaims

Removed species-divergent remodeling and evolutionary-calibrated pan-mammalian claims. The abstract now states that neither compositional concordance nor divergence was demonstrated and reports measured clock errors with predictive R-squared.

### Task 3: Human transport

Changed `scripts/run_m9_pipeline.py`, `scripts/run_m10_pipeline.py`, and `scripts/run_human_transport_inference.py`.

M9 now uses per-cell library-size normalization. Both signature tiers and raw-count sensitivity scores are exported. One BH family covers 225 public-cohort tests and 27 primary PsychAD tests.

MSSM whole-cohort results changed from slope -0.000193, p = 0.02825, and q = 0.3107 to slope -0.000404, p = 0.240302, and q = 0.940170.

Primary PsychAD whole-cortex results changed from the P0 unadjusted slope +0.000318 with p = 0.01182 to the P1 slope +0.000186 with p = 0.349415. The adjusted result changed from +0.000027 with p = 0.862622 to -0.000209 with p = 0.400449.

Three global survivors remain, but their identities changed after normalization and tier correction: MSSM L2/3 neurons, MSSM oligodendrocytes, and PsychAD EN_L3_5_IT_3 neurons. No OPC or microglial result survives the primary family. PsychAD microglia has pooled adjusted slope +0.002069, p = 0.001496, and global q = 0.062818. Adult microglia has slope +0.003115, p = 0.000617, and secondary adult q = 0.016655. Glial transport is a candidate.

Added `results/m9-m10-donor-overlap-audit.csv`. Candidate age, sex, and bank matches exist for 160 of 214 public MSSM donors and 133 of 134 public HBCC donors. Unique-candidate counts are 52 and 35. These are not confirmed identity matches, so PsychAD is not called independent replication.

### Task 4: Signature definition

Added `src/aging_programs.py` and `results/signature-tier-summary.csv`. P1 is primary throughout, with P0 as sensitivity. All result rows identify the tier and normalization.

Verification found that the remediation brief misdescribed P1. P1 applies technical, biotype, ribosomal, mitochondrial, and ambient-RNA filters. It does not remove sex-dimorphic or immediate-early genes. The manuscript now states the actual filters and makes no AP-1-independence claim.

### Task 5: Clock reporting

Changed `scripts/run_h1_epigenetic_drift.py`.

Added per-sample prediction files and predictive R-squared. Mouse liver changed from squared Pearson correlation 0.851 to predictive R-squared 0.3169. Human multi-tissue changed from 0.642 to 0.6383. MAE remains 0.36 years in mouse and 5.06 years in human.

Mouse scoring uses 724 of 816 probes. `results/h1-training-overlap-audit.csv` records that GSE184211 belongs to the consortium development resource. Mouse sample-level training membership remains unresolved, so neither evaluation is called independent validation.

### Task 6: Positive controls

Restored Apoe, Il6, Lgals3, and Igfbp3 in manuscript Methods and Results using `results/m2-calibration-positive-controls.csv`. All four are non-significant. No unverified values were added.

### Task 7: Literature comparison

Added comparisons with Tyshkovskiy et al. 2026 and Karakaslar et al. 2023. The manuscript distinguishes bulk multi-species modules, immune AP-1 mechanisms, shared TMS inputs, and Agescape's single-signature transport estimand. Full-text sources are recorded in `notes/remediation-number-map.md`.

### Task 8: Data-use terms

Rendered and inspected the current AD Knowledge Portal Data Access terms. Aggregate regression outputs may be published with required acknowledgments. Individual-level controlled data and derived donor-level data may not be redistributed. The current aggregate-only approach complies. See `notes/adkp-data-use-verification.md`.

### Task 9: Housekeeping

- Confirmed 11 public cohorts and 621 scored donors after recomputation.
- Confirmed 661 human tissue rows and 339 mouse rows, totaling 1,000 methylation profiles.
- Retained the corrected 16-mouse lung values.
- Corrected all nine H3 cell counts and added the missing intermediate-monocyte sensitivity arm.
- Added predictive-metric, signature, overlap, and cohort summary tables.
- Updated README, project overview, app contract, app copy, and G16 rendering.
- Added support for flat and nested manifest schemas in `scripts/verify_data.py`.
- Fixed an M1 pandas categorical fill compatibility failure.
- Removed local absolute paths from M9 provenance logging.

## 3. Abstract number mapping

| Abstract claim or value | Committed source | Row | Column |
|---|---|---|---|
| 11 public human cohorts | `results/analysis-cohort-summary.csv` | `m9_cohorts` | `value` |
| 621 scored control donors | `results/analysis-cohort-summary.csv` | `m9_scored_control_donors` | `value` |
| 284 PsychAD control donors | `results/analysis-cohort-summary.csv` | `psychad_control_donors` | `value` |
| 1,000 methylation profiles | `results/analysis-cohort-summary.csv` | `methylation_array_profiles` | `value` |
| 139 primary genes, 257 sensitivity genes | `results/signature-tier-summary.csv` | `P1_Clean`, `P0_Raw257` | `n_mouse_genes` |
| Primary whole-cortex slope +0.000186, p = 0.3494 | `results/m10-cohort-summary.csv` | `DLPFC_Lifespan_Atlas` | `slope`, `p` |
| Adjusted whole-cortex slope -0.000209, p = 0.4004 | `results/m10-cohort-summary.csv` | `DLPFC_Lifespan_Atlas` | `slope_covar_adj`, `p_covar_adj` |
| 252 primary tests and 3 survivors | `results/analysis-cohort-summary.csv` | `primary_human_transport_tests`, `primary_human_transport_survivors` | `value` |
| MSSM L2/3 neuron slope -0.001303, q = 0.0208 | `results/human-transport-global-multiplicity.csv` | `transport_022` | `slope_per_year`, `q_global_bh` |
| MSSM oligodendrocyte slope +0.001152, q = 0.0346 | `results/human-transport-global-multiplicity.csv` | `transport_023` | `slope_per_year`, `q_global_bh` |
| PsychAD EN_L3_5_IT_3 slope -0.001383, q = 0.0208 | `results/human-transport-global-multiplicity.csv` | `transport_230` | `slope_per_year`, `q_global_bh` |
| Fibroblast Cdkn1a slope +0.1080, p = 0.0029, q = 0.0322 | `results/fibroblast-markers-audit.csv` | `Cdkn1a` | `slope_per_month`, `p_value`, `q_value_bh` |
| Human club-cell slope -0.0602, q = 0.0044 | `results/h2-composition-human-lung.csv` | `club cell` | `slope_per_year_or_month`, `q_bh` |
| Human multiciliated-cell slope -0.0342, q = 0.0225 | `results/h2-composition-human-lung.csv` | `multiciliated epithelial cell` | `slope_per_year_or_month`, `q_bh` |
| 44 human lung cell types | `results/analysis-cohort-summary.csv` | `human_lung_tested_cell_types` | `value` |
| Multiciliated detection caveat | `results/h2-composition-human-lung-presence.csv` | `multiciliated epithelial cell` | `presence_age_r`, `presence_age_p`, `q_bh` |
| Club-cell study-adjusted p = 0.1033 | `results/d3-human-composition-study-sensitivity.csv` | `club cell` | `study_adj_p` |
| 3 of 7 concordant, binomial p = 1.000, rho = 0.000 | `results/h2-composition-cross-species-test.csv` | first row | `n_sign_concordant`, `n_pairs`, `binomial_two_sided_p`, `spearman_rho` |
| Mouse clock MAE 0.36 years, predictive R-squared 0.3169 | `results/h1-epigenetic-clock-mouse.csv` | `Mouse`, `Liver` | `mae_years`, `predictive_r2` |
| Human clock MAE 5.06 years, predictive R-squared 0.6383, 661 profiles | `results/h1-epigenetic-clock-human.csv` | `Human`, `All Tissues` | `mae_years`, `predictive_r2`, `n_samples` |
| 724 of 816 mouse clock probes | `results/analysis-cohort-summary.csv` | `mouse_clock_probes_used`, `universal_clock_probes` | `value` |
| Five NMR animals and no detected slope | `results/m6-statistical-summary.csv` | `Naked mole-rat skin aging slope (n=5 animals)` | `n`, `slope`, `p` |

No abstract number remains without a committed CSV source. A second copy of this map is in `notes/remediation-number-map.md`.

## 4. Forbidden-phrase scan

A case-insensitive scan of `report/`, `README.md`, and `PROJECT_OVERVIEW.md` found zero instances of:

- confirms flatness
- plateaus completely
- proves
- highly calibrated
- odometer
- species-divergent tissue remodeling
- pre-registered
- external validation

The manuscript, README, and project overview contain no em dashes or en dashes.

## 5. Verification outputs

Final outputs are recorded below after the last export and test run.

### Full pipeline

M1 through M6 completed in the main run. Notebook 06 initially failed on an outdated column name and NMR interpretation. After that notebook was corrected, all six notebooks executed successfully. The remaining app, H1, H2, H3, and marker stages then completed in a resumed invocation. The full set of stages passed, but this was not one uninterrupted `run_all.py` invocation.

### pytest

```text
============================= test session starts ==============================
platform linux -- Python 3.11.17, pytest-9.1.1, pluggy-1.7.0 -- /home/reza/miniconda3/envs/agescape/bin/python
cachedir: .pytest_cache
rootdir: /home/reza/projects/agescape
plugins: platformdirs-4.13.0, zarr-3.1.6, anyio-4.15.1
collecting ... collected 20 items

tests/test_clean_clone_and_numbers.py::test_anchor_lock_against_committed_results PASSED [  5%]
tests/test_clean_clone_and_numbers.py::test_app_data_check_script PASSED [ 10%]
tests/test_clean_clone_and_numbers.py::test_app_data_sizes_and_count PASSED [ 15%]
tests/test_human_transport.py::test_primary_signature_tiers_and_global_multiplicity PASSED [ 20%]
tests/test_human_transport.py::test_sensitivity_arms_use_distinct_tiers_and_normalization PASSED [ 25%]
tests/test_hypotheses_2_and_3.py::test_h2_composition_outputs_exist PASSED [ 30%]
tests/test_hypotheses_2_and_3.py::test_h2_composition_simplex_invariant PASSED [ 35%]
tests/test_hypotheses_2_and_3.py::test_h2_composition_results_content PASSED [ 40%]
tests/test_hypotheses_2_and_3.py::test_h3_entropy_outputs_exist_and_valid PASSED [ 45%]
tests/test_hypotheses_2_and_3.py::test_fibroblast_markers_audit PASSED   [ 50%]
tests/test_hypotheses_2_and_3.py::test_h2_and_h3_figures_exist PASSED    [ 55%]
tests/test_hypotheses_2_and_3.py::test_h1_epigenetic_outputs_exist_and_accurate PASSED [ 60%]
tests/test_hypotheses_2_and_3.py::test_h2_prethreshold_counts_and_presence_outputs PASSED [ 65%]
tests/test_hypotheses_2_and_3.py::test_h1_predictive_r2_recomputes_from_per_sample_predictions PASSED [ 70%]
tests/test_security_and_dua.py::test_no_controlled_data_or_tokens_in_git PASSED [ 75%]
tests/test_security_and_dua.py::test_zero_innerhtml_outerhtml_in_app PASSED [ 80%]
tests/test_security_and_dua.py::test_zero_inline_styles_in_index_html PASSED [ 85%]
tests/test_security_and_dua.py::test_no_third_party_fonts_in_index_html PASSED [ 90%]
tests/test_security_and_dua.py::test_csp_headers_snippet_exists PASSED   [ 95%]
tests/test_security_and_dua.py::test_no_local_user_paths_in_tracked_manifests PASSED [100%]

============================== 20 passed in 2.12s ==============================
```

### App-data gate

```text
ANCHOR                                                                       EXPECTED                                             OBSERVED                                             STATUS
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
export exists: atlas.json                                                    True                                                 True                                                 PASS
export exists: ranking.json                                                  True                                                 True                                                 PASS
export exists: calibration.json                                              True                                                 True                                                 PASS
export exists: agreement.json                                                True                                                 True                                                 PASS
export exists: design-grid.json                                              True                                                 True                                                 PASS
export exists: program.json                                                  True                                                 True                                                 PASS
export exists: signature.json                                                True                                                 True                                                 PASS
export exists: savina.json                                                   True                                                 True                                                 PASS
export exists: power-grid.json                                               True                                                 True                                                 PASS
export exists: provenance.json                                               True                                                 True                                                 PASS
export exists: anage.json                                                    True                                                 True                                                 PASS
export exists: humancohorts.json                                             True                                                 True                                                 PASS
export exists: dlpfc.json                                                    True                                                 True                                                 PASS
atlas units == 163                                                           163                                                  163                                                  PASS
atlas ranked == 17                                                           17                                                   17                                                   PASS
atlas powered == 17                                                          17                                                   17                                                   PASS
powered confined to Lung & Kidney                                            {'Lung', 'Kidney'}                                   {'Lung', 'Kidney'}                                   PASS
per-tissue powered {Kidney:8, Lung:9}                                        {'Kidney': 8, 'Lung': 9}                             {'Lung': 9, 'Kidney': 8}                             PASS
sum of per-tissue powered == 17                                              17                                                   17                                                   PASS
every unit carries all schema-1 fields                                       True                                                 True                                                 PASS
ranking head counts 55/53/36/35                                              [55, 53, 36, 35]                                     [55, 53, 36, 35]                                     PASS
ranking head units                                                           ['Lung::adventitial cell', 'Kidney::kidney proximal convoluted tubule epithelial cell', 'Lung::fibroblast of lung', 'Lung::alveolar macrophage'] ['Lung::adventitial cell', 'Kidney::kidney proximal convoluted tubule epithelial cell', 'Lung::fibroblast of lung', 'Lung::alveolar macrophage'] PASS
calibration observed 55                                                      55                                                   55                                                   PASS
calibration perm_mean 4.6                                                    4.6                                                  4.6                                                  PASS
calibration sex_floor 0                                                      0                                                    0                                                    PASS
calibration unit                                                             'Lung::adventitial cell'                             'Lung::adventitial cell'                             PASS
positive_controls == 4 rows                                                  4                                                    4                                                    PASS
agreement counts 13/4/0                                                      {'agree': 13, 'partial': 4, 'disagree': 0}           {'agree': 13, 'disagree': 0, 'partial': 4}           PASS
lollipop == 17 rows                                                          17                                                   17                                                   PASS
scatter unit                                                                 'Kidney::kidney proximal convoluted tubule epithelial cell' 'Kidney::kidney proximal convoluted tubule epithelial cell' PASS
scatter n == 12713                                                           12713                                                12713                                                PASS
scatter points == 12713                                                      12713                                                12713                                                PASS
scatter rho field == 0.569                                                   0.569                                                0.569                                                PASS
recomputed rho rounds to 0.569                                               0.569                                                0.569                                                PASS
|recomputed rho - 0.569| < 0.01                                              True                                                 True                                                 PASS
drawn rho matches recompute (full precision)                                 True                                                 True                                                 PASS
program summary up 166 / down 91 / entries 312                               {'up': 166, 'down': 91, 'entries': 312}              {'down': 91, 'entries': 312, 'up': 166}              PASS
program genes == 312 entries                                                 312                                                  312                                                  PASS
program units == 17                                                          17                                                   17                                                   PASS
summary up == m6-program-genes program_up non-null                           166                                                  166                                                  PASS
summary down == m6-program-genes program_down non-null                       91                                                   91                                                   PASS
matrix_top40 rows == 40                                                      40                                                   40                                                   PASS
matrix_top40 columns == 17 units                                             17                                                   17                                                   PASS
matrix rows ranked by units-present (desc)                                   True                                                 True                                                 PASS
mouse_wb slope                                                               0.13437                                              0.13437                                              PASS
mouse_wb ci                                                                  [0.02926, 0.23949]                                   [0.02926, 0.23949]                                   PASS
mouse_wb p                                                                   0.0147                                               0.0147                                               PASS
mouse_wb n                                                                   23                                                   23                                                   PASS
mouse_skin slope                                                             0.13437                                              0.13437                                              PASS
mouse_skin ci                                                                [0.02926, 0.23949]                                   [0.02926, 0.23949]                                   PASS
mouse_skin p                                                                 0.0147                                               0.0147                                               PASS
mouse_skin n                                                                 3                                                    3                                                    PASS
mouse_skin underpowered                                                      True                                                 True                                                 PASS
human_skin slope                                                             -0.59792                                             -0.59792                                             PASS
human_skin ci                                                                [-2.81388, 1.61803]                                  [-2.81388, 1.61803]                                  PASS
human_skin p                                                                 0.3655                                               0.3655                                               PASS
human_skin n                                                                 4                                                    4                                                    PASS
human_skin low_powered                                                       True                                                 True                                                 PASS
human_fibroblast slope                                                       -0.59792                                             -0.59792                                             PASS
human_fibroblast ci                                                          [-2.81388, 1.61803]                                  [-2.81388, 1.61803]                                  PASS
human_fibroblast p                                                           0.3655                                               0.3655                                               PASS
human_fibroblast n                                                           4                                                    4                                                    PASS
nmr slope                                                                    0.00584                                              0.00584                                              PASS
nmr ci                                                                       [-0.19812, 0.20979]                                  [-0.19812, 0.20979]                                  PASS
nmr p                                                                        0.9332                                               0.9332                                               PASS
nmr n                                                                        5                                                    5                                                    PASS
nmr tested_span_years [1.1, 11.3]                                            [1.1, 11.3]                                          [1.1, 11.3]                                          PASS
mouse_wb points == 23                                                        23                                                   23                                                   PASS
slope-diff 0.12854                                                           0.12854                                              0.12854                                              PASS
slope-diff ci                                                                [-0.03144, 0.28851]                                  [-0.03144, 0.28851]                                  PASS
slope-diff p 0.1153                                                          0.1153                                               0.1153                                               PASS
TOST delta 0.202                                                             0.202                                                0.202                                                PASS
TOST p 0.02748                                                               0.02748                                              0.02748                                              PASS
TOST established == true                                                     True                                                 True                                                 PASS
permutation p 1.0                                                            1.0                                                  1.0                                                  PASS
coverage mouse 1.0                                                           1.0                                                  1.0                                                  PASS
coverage human 0.84892                                                       0.84892                                              0.84892                                              PASS
coverage nmr 0.66187                                                         0.66187                                              0.66187                                              PASS
coverage values map to the right species (vs m6-gates.csv)                   True                                                 True                                                 PASS
lifespans mouse 4.0 / nmr 31.0 / human 122.5                                 {'mouse': 4.0, 'nmr': 31.0, 'human': 122.5}          {'human': 122.5, 'mouse': 4.0, 'nmr': 31.0}          PASS
power donors axis                                                            [4, 6, 8, 10, 12, 16, 20, 30, 50]                    [4, 6, 8, 10, 12, 16, 20, 30, 50]                    PASS
power sds axis                                                               [8, 12, 16, 20]                                      [8, 12, 16, 20]                                      PASS
power dispersions axis                                                       [0.1, 0.2, 0.4]                                      [0.1, 0.2, 0.4]                                      PASS
power tiles == 108                                                           108                                                  108                                                  PASS
anchor 1.474 @ (4 donors, SD 8, disp 0.1)                                    1.474                                                1.474                                                PASS
anchor 0.59 @ (4 donors, SD 20, disp 0.1)                                    0.59                                                 0.59                                                 PASS
anchors dict mirrors tiles                                                   {'4_8_low': 1.474, '4_20_low': 0.59}                 {'4_20_low': 0.59, '4_8_low': 1.474}                 PASS
anage species == 4141                                                        4141                                                 4141                                                 PASS
anage mouse 4.0 / nmr 31.0 / human 122.5 / bowhead 211.0 / quahog 507.0      {'Mus musculus': 4.0, 'Heterocephalus glaber': 31.0, 'Homo sapiens': 122.5, 'Balaena mysticetus': 211.0, 'Arctica islandica': 507.0} {'Mus musculus': 4.0, 'Heterocephalus glaber': 31.0, 'Homo sapiens': 122.5, 'Balaena mysticetus': 211.0, 'Arctica islandica': 507.0} PASS
anage single-cell flags: exactly mouse/human/nmr                             True                                                 True                                                 PASS
humancohorts cohorts == 11                                                   11                                                   11                                                   PASS
humancohorts total scored donors == 621                                      621                                                  621                                                  PASS
M9 MSSM n_donors 213, verdict POWERED-FLOOR                                  (213, 'POWERED-FLOOR')                               (213, 'POWERED-FLOOR')                               PASS
M9 MSSM whole-cohort result does not survive BH                              True                                                 True                                                 PASS
M9 MSSM slope_per_year -0.00040                                              -0.0004                                              -0.0004                                              PASS
M9 MSSM BH q across 11 cohorts 0.940                                         0.94                                                 0.94                                                 PASS
M9 cell-type hits present (>=10 rows)                                        True                                                 True                                                 PASS
M9 cohorts span brain/kidney/lung                                            {'lung', 'brain', 'kidney'}                          {'kidney', 'brain', 'lung'}                          PASS
M10 n_donors 284, 27 cell types, age 0.2-97                                  (284, 27, 0.2, 97.0)                                 (284, 27, 0.2, 97.0)                                 PASS
M10 age_sd 24.0                                                              24.0                                                 24.0                                                 PASS
M10 whole cortex covar-adjusted p > 0.1                                      True                                                 True                                                 PASS
M10 whole cortex P1 CPM slope +0.00019                                       0.00019                                              0.00019                                              PASS
M10 celltypes == 27 rows                                                     27                                                   27                                                   PASS
M10 concordant == 9                                                          9                                                    9                                                    PASS
M10 Oligo P1 slope +0.00115                                                  0.00115                                              0.00115                                              PASS
M10 Oligo p < 1e-3, concordant up                                            True                                                 True                                                 PASS
M10 Astro nominal p in (0.01, 0.05)                                          True                                                 True                                                 PASS
M10 Astro concordant up                                                      True                                                 True                                                 PASS
M10 PVM P1 slope -0.00444                                                    -0.00444                                             -0.00444                                             PASS
M10 PVM concordant down                                                      True                                                 True                                                 PASS
M10 Endo direction reverses, not concordant                                  True                                                 True                                                 PASS
M10 P1 program coverage 0.9576, 113 measured orthologs                       True                                                 True                                                 PASS
M10 BH-FDR q<0.05 set is non-empty                                           True                                                 True                                                 PASS
savina curve == 39 rows                                                      39                                                   39                                                   PASS
savina paper_rule 29 @ 1.0 (p<0.05)                                          {'cutoff': 1.0, 'n_pval': 29}                        {'cutoff': 1.0, 'n_pval': 29}                        PASS
savina converted_rule 9 & 6 @ 1.45                                           {'cutoff': 1.45, 'n_pval': 9, 'n_fdr': 6}            {'cutoff': 1.45, 'n_fdr': 6, 'n_pval': 9}            PASS
savina loose_rule 560 @ 0.25 (FDR)                                           {'cutoff': 0.25, 'n_fdr': 560}                       {'cutoff': 0.25, 'n_fdr': 560}                       PASS
design-grid age_groups                                                       ['1m', '3m', '18m', '21m', '24m', '30m']             ['1m', '3m', '18m', '21m', '24m', '30m']             PASS
cells female total == m1 All row                                             81362                                                81362                                                PASS
cells male total == m1 All row                                               164027                                               164027                                               PASS
mice female total == m1 All row                                              6                                                    6                                                    PASS
mice male total == m1 All row                                                17                                                   17                                                   PASS
cells female per group == source                                             [0, 23693, 21841, 35828, 0, 0]                       [0, 23693, 21841, 35828, 0, 0]                       PASS
cells male per group == source                                               [25980, 21909, 22804, 0, 37660, 55674]               [25980, 21909, 22804, 0, 37660, 55674]               PASS
female-only columns {21m}                                                    ['21m']                                              ['21m']                                              PASS
male-only columns {1m,24m,30m}                                               ['1m', '24m', '30m']                                 ['1m', '24m', '30m']                                 PASS
both-sex columns {3m,18m}                                                    ['3m', '18m']                                        ['3m', '18m']                                        PASS
freeze_date 2026-10-14                                                       '2026-10-14'                                         '2026-10-14'                                         PASS
export_date present                                                          True                                                 True                                                 PASS
all provenance SHA-256 values match files on disk                            True                                                 True                                                 PASS
every export has provenance sources                                          12                                                   12                                                   PASS
figure assets recorded == 4                                                  4                                                    4                                                    PASS
figure asset bytes and SHA-256 match app/figures                             True                                                 True                                                 PASS
every file < 1 MB                                                            True                                                 True                                                 PASS
total app/data < 5 MB                                                        True                                                 True                                                 PASS
S1: Lung::adventitial cell (55) ranks above Kidney PCT (53) in ranking.json  True                                                 True                                                 PASS
S1: calibration object names Lung::adventitial cell                          'Lung::adventitial cell'                             'Lung::adventitial cell'                             PASS
S2: preset (donors 4, SD 8, disp 0.1) -> 1.474 & LOW-POWERED                 True                                                 True                                                 PASS
S2: preset (donors 4, SD 20, disp 0.1) -> 0.59 & LOW-POWERED                 True                                                 True                                                 PASS
S3: Chapter 5 donor requirement statements map to power-grid.json tuples     True                                                 True                                                 PASS
S4: Chapter 1 mouse floor (3 mice, 50 cells) with no human floor leak        True                                                 True                                                 PASS
S4: Chapter 5 human floor (6 donors, 12 years SD) with no mouse floor leak   True                                                 True                                                 PASS
B1: zero pseudotime alignment/tracks age claims near pseudotime              []                                                   []                                                   PASS
B2: zero complete concordance, young-vs-young, 21 other organs, Savina 2024  []                                                   []                                                   PASS
B3: zero unbacked enrich, pathway, regulator terms                           []                                                   []                                                   PASS
B4: zero usability PASS labels without external tester evidence              []                                                   []                                                   PASS
-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
PASS: all 146 anchors hold.
```

### Additional checks

- Python compilation of `scripts/` and `src/`.
- `git diff --check`.
- `scripts/check_app_assets.py`.
- Browser rendering at 1440 by 1000 and 390 by 844 pixels.
- G16 console error count: 0 at both sizes.
- G16 canvas count: 1 at both sizes.
- Page width equals viewport width at both sizes.

## 6. Open decisions and unresolved inputs

1. Site title wording about conservation remains unchanged for Reza.
2. AI-disclosure placement remains unchanged for Reza.
3. Exact M9/M10 donor overlap cannot be established from shared identifiers. Only age, sex, and bank match bounds are reported.
4. GSE199979 sample-level membership in clock training data remains unresolved. No independent-validation claim is made.
5. P1 does not remove sex-dimorphic or immediate-early genes. A future sex- and immediate-early-trimmed primary analysis would require a new decision and recomputation.

## 7. External actions

Nothing was submitted. No pull request was opened. No email was sent. No issue or external review was created. No remote push was performed. This report stops at final preview.
