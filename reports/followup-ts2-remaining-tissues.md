# Follow-up Item #15: Tabula Sapiens 2.0 Complete 28-Tissue Audit and Remaining >=6-Donor Tissues

## 1. Executive Summary and Primary-Source Scope

This report completes the audit of all 28 tissues in Tabula Sapiens 2.0 (CELLxGENE collection `e5f58829`, GEO `GSE306755`, bioRxiv `2024.12.03.626516v1`) and analyzes every tissue outside previously completed `Skin` (TS1 4-donor baseline) and `Blood` natural killer cells (`reports/followup-tabula-sapiens-2.md`) that meets the minimum single-assay inclusion floor (`>=6` adult age-verified donors on a single platform).

Key conclusions:
1. **Strict single-assay >=6-donor eligibility across 28 tissues.** Reconciling supplementary workbooks (`results/followup-ts2-organ-summary.csv`, `results/followup-ts2-donor-organ-crosswalk.csv`) with direct H5AD inspection (`results/followup-ts2-all28-tissues-audit.csv`) leaves **five analyzed tissues**, all restricted to `10x 3' v3`: `Blood` non-NK cells (`9` donors, `33 to 69 yr`, `4M/5F`, `72,776` non-NK cells across `6` non-NK cell types + `All_NonNK_Cells`), `Bladder` (`6` donors, `45 to 61 yr`, `2M/4F`, `64,214` cells across `8` cell types + `All_Cells`), `Fat` (`6` donors, `33 to 61 yr`, `3M/3F`, `93,195` cells across `2` cell types + `All_Cells`), `Lymph_Node` (`6` donors, `45 to 69 yr`, `2M/4F`, `98,431` cells across `10` cell types + `All_Cells`), and `Spleen` (`6` donors, `45 to 69 yr`, `2M/4F`, `55,454` cells across `10` cell types + `All_Cells`).
2. **Bone Marrow excluded after H5AD inspection.** Although `Bone_Marrow` (`4aedb645`) has `7` donors overall across assays, direct H5AD inspection showed only `5` donors on `10x 3' v3` (`TSP2`, `TSP14`, `TSP21`, `TSP25`, `TSP27`), `5` on `10x 5' v2`, and `4` on `Smart-seq2`, with `LinEnriched` and `LinDepleted` FACS sorts. It fails the single-assay `>=6`-donor threshold. Three other 6-donor tissues (`Heart`, `Muscle`, `Trachea`) likewise have only `5` age-verified `10x` donors, and the remaining `19` tissues (plus already-analyzed `Skin` at `n=4`) have `1 to 5` donors overall.
3. **Why Family A fractions and cross-cell-type concordance are excluded.** As documented in Tabula Sapiens (`PMC9812260`), organ cell suspensions were compartment-normalized using magnetic microbeads or FACS across epithelial, endothelial, immune, and stromal fractions, and `Blood` `.obs['notes']` explicitly mixes `Enriched` (`TSP10`: `3,167` cells), `Whole` (`TSP10`: `1,193` cells), `PBMC` (`TSP25`: `1,719` cells), and `SheelaPrep` (`TSP2`: `1,542` cells) alongside `None` (`9,142` in `TSP2`). Accordingly, **all TS2 tissues run Family B (donor-level gene OLS) and Family C (`P1_Clean` signature) ONLY; Family A CLR composition fractions and unbiased cross-cell-type concordance are omitted.**
4. **Shared-donor non-independence.** Only **`9` unique donors** (`TSP1`, `TSP2`, `TSP7`, `TSP8`, `TSP10`, `TSP14`, `TSP21`, `TSP25`, `TSP27`) contribute across all 5 analyzed tissues (`results/followup-ts2-remaining-donors.csv`). Five core donors (`TSP2`, `TSP14`, `TSP21`, `TSP25`, `TSP27`) appear in every tissue, and `Lymph_Node` and `Spleen` share the exact same 6 donors.
5. **Family B gene-level results are globally null.** Across `36` cell-type groups (`598,189` tests) and `5` pooled organ groups (`688,294` total tests), **zero genes pass cross-tissue or per-tissue BH FDR < 0.05**. Only `2` genes pass exploratory within-single-group BH (`Spleen mature NK T cell AIG1`, `q_within = 0.0292`, `q_tissue = 0.175`, `q_ts2_ct36 = 0.634`; `Spleen classical monocyte IGHGP`, `q_within = 0.0314`, `q_tissue = 0.175`, `q_ts2_ct36 = 0.634`), both on `N = 6`, `residual_df = 2`.
6. **Family C (`P1_Clean`) results are unstable, small-sample, and directionally mixed.** At `N = 6` (`residual_df = 2` primary, `df = 3` in the `no_log_lib` sensitivity), permutation nulls show low `empirical_null_pvalue` tails with conflicting composite directions across tissues (`+` in `Bladder macrophage` `q_ts2_all41_empirical = 0.00683` and `Fat CD8+ T` `q_ts2_all41_empirical = 0.00683` vs `-` in `Lymph_Node All_Cells` `q_ts2_all41_empirical = 0.00683`), and parametric composite effects shift substantially when dropping `log_lib`. Neither a directional age trend nor a consistent P1 signature replication can be claimed from `N = 6` across `16 to 28` yr age spans.

## 2. Complete 28-Tissue Audit and H5AD Integer Count Gate

All 28 tissues in `results/followup-ts2-all28-tissues-audit.csv`:

- **Analyzed (`5` tissues, `10x 3' v3`, Family B + C only):** `Blood` (non-NK, `n=9`), `Bladder` (`n=6`), `Fat` (`n=6`), `Lymph_Node` (`n=6`), `Spleen` (`n=6`). (`Blood` NK cells were previously analyzed separately on `7` `10x 3' v3` donors and excluded from `All_NonNK_Cells` and `Blood` cell types here.)
- **Excluded at single-assay `>=6` floor (`4` tissues with `6 to 7` total donors but `5` on `10x`):** `Bone_Marrow` (`7` total -> `5` `10x 3' v3`, `5` `10x 5' v2`, H5AD directly audited), `Heart` (`6` total -> `5` `10x`), `Muscle` (`6` total -> `5` `10x`), `Trachea` (`6` total -> `5` `10x`).
- **Already analyzed / below floor (`1` tissue):** `Skin` (`4` donors `33 to 61 yr`, previously covered as TS1 baseline).
- **Excluded `<6` donors overall (`18` tissues):** `Ear` (`2`), `Eye` (`4`), `Kidney` (`1`), `Large_Intestine` (`5`), `Liver` (`4`), `Lung` (`4`), `Mammary` (`2`), `Ovary` (`3`), `Pancreas` (`2` age-verified), `Prostate` (`3`), `Salivary_Gland` (`4`), `Small_Intestine` (`5`), `Stomach` (`2` age-verified), `Testis` (`1`), `Thymus` (`5`), `Tongue` (`5` total / `4` `10x`), `Uterus` (`2`), `Vasculature` (`5`).

All 5 analyzed H5ADs passed the full raw integer count gate (`raw/X` finite, nonnegative, integer-valued; `results/followup-ts2-remaining-metadata-summary.json`):

| Tissue | Analyzed Donors (`10x 3' v3`) | Age Span (yr) | Sex (M/F) | Eligible Cell Types + Pool | Total Nonzeros in `raw/X` | Nonfinite / Negative / Noninteger |
|---|---:|---|---:|---:|---:|---:|
| Blood (non-NK) | 9 | 33 to 69 | 4 / 5 | 6 + `All_NonNK_Cells` | 201,330,996 | 0 / 0 / 0 |
| Bladder | 6 | 45 to 61 | 2 / 4 | 8 + `All_Cells` | 266,974,746 | 0 / 0 / 0 |
| Fat | 6 | 33 to 61 | 3 / 3 | 2 + `All_Cells` | 295,933,576 | 0 / 0 / 0 |
| Lymph_Node | 6 | 45 to 69 | 2 / 4 | 10 + `All_Cells` | 322,710,450 | 0 / 0 / 0 |
| Spleen | 6 | 45 to 69 | 2 / 4 | 10 + `All_Cells` | 207,399,916 | 0 / 0 / 0 |

Cross-check of `.obs['development_stage']` and `.obs['sex']` inside all five H5ADs matched `results/followup-ts2-donor-summary.csv` with zero discrepancies (`scripts/run_followup_ts2_remaining.py`).

## 3. Multiple-Testing Families Pre-Declared

- **Family A (Composition CLR fractions):** Not run (`PMC9812260` compartment microbead/FACS normalization + `Blood` prep heterogeneity in `results/followup-ts2-remaining-donors.csv`).
- **Family B (Donor-level Gene OLS `log2CPM ~ age_per_decade + sex_male + log_lib`, min `>=6` donors with `>=20` cells/donor, `mean/frac` filter `CPM >= 1.0` in `>=50%` donors, `residual_df >= 2`):**
  - Pre-declared cross-tissue primary cell-type family: all **`598,189`** tested `(tissue, cell_type, gene)` pairs across the `36` individual cell-type groups (`q_bh_ts2_celltypes`).
  - Pre-declared all-group sensitivity family (including `5` pooled organ groups): **`688,294`** tests (`q_bh_ts2_all41`).
  - Secondary within-tissue cell-type families: `83,051` (`Blood` 6 cell types), `137,712` (`Bladder` 8 cell types), `35,138` (`Fat` 2 cell types), `177,672` (`Lymph_Node` 10 cell types), `164,616` (`Spleen` 10 cell types) (`q_bh_tissue_celltypes`). Exploratory single-group BH (`q_bh_within_group`) is also recorded for completeness.
  - Cross-cell-type concordance omitted for sorted/compartment-normalized organs (`is_unbiased_single_cell=False`).
- **Family C (`P1_Clean` 139-gene mouse signature, `2,000` random gene-set draws matched by sign cardinality and expressed pool):**
  - Evaluated per tissue (`7`, `9`, `3`, `11`, `11` groups -> `q_bh_tissue_composite`, `q_bh_tissue_empirical`) and jointly across all **`41`** TS2 groups (`q_bh_ts2_all41_composite`, `q_bh_ts2_all41_empirical`).
  - Pre-declared `df`-preservation sensitivity (`sensitivity_age_sex_no_log_lib`): drops `log_lib` so `N = 6` models move from `residual_df = 2` to `residual_df = 3` (`results/followup-ts2-remaining-sensitivities-p1.csv`).

## 4. Family B Results: Donor-Level Gene Age Tests

Full summary statistics are in `results/followup-ts2-remaining-genes-summary.csv` and top hits in `results/followup-ts2-remaining-top-genes.csv`:

| Tissue | Cell-Type Groups | Cell-Type Gene Tests | Pooled Group Tests | Min Parametric `p` | Min Per-Tissue Cell-Type `q` | Sig at `q_ts2_ct36 < 0.05` | Sig at `q_within < 0.05` |
|---|---:|---:|---:|---:|---:|---:|---:|
| Blood (non-NK) | 6 | 83,051 | 15,033 (`All_NonNK_Cells`) | 5.35e-5 (`classical monocyte PMM2`) | 0.99997 | 0 | 0 |
| Bladder | 8 | 137,712 | 18,532 (`All_Cells`) | 3.39e-6 (`CD4+ T TRIM52`) | 0.467 | 0 | 0 |
| Fat | 2 | 35,138 | 19,285 (`All_Cells`) | 4.89e-5 (`CD8+ T SNRNP40`) | 0.99998 | 0 | 0 |
| Lymph_Node | 10 | 177,672 | 18,877 (`All_Cells`) | 6.39e-6 (`mature NK T POLR2J3`) | 0.706 | 0 | 0 |
| Spleen | 10 | 164,616 | 18,378 (`All_Cells`) | 1.71e-6 (`mature NK T AIG1`) | 0.175 | 0 | 2 (`AIG1`, `IGHGP`) |
| **Total TS2 Remaining** | **36** | **598,189** | **90,105 (total 688,294)** | **1.71e-6** | **0.175** | **0** | **2 (Spleen only)** |

The two exploratory within-group hits in `Spleen` (`AIG1` in `mature NK T cell`: effect `-0.949` log2CPM/decade, `p = 1.71e-6`, `q_within = 0.0292`, `q_tissue = 0.175`, `q_ts2_ct36 = 0.634`; `IGHGP` in `classical monocyte`: effect `+2.523`, `p = 2.12e-6`, `q_within = 0.0314`, `q_tissue = 0.175`, `q_ts2_ct36 = 0.634`) sit on `N = 6` (`df = 2`) over a narrow `45 to 69 yr` span and fail per-tissue and cross-tissue FDR.

## 5. Family C Results: P1_Clean 139-Gene Signature

Full results across all 41 groups are in `results/followup-ts2-remaining-p1-signature.csv` and `results/followup-ts2-remaining-sensitivities-p1.csv`:

- **Parametric composite t-test (`q_bh_ts2_all41_composite`):** **0 groups** pass `< 0.05` across the 41-group family (minimum `q_bh_ts2_all41_composite = 0.224` in `Bladder macrophage`). Within tissue, only `Bladder macrophage` reaches `q_bh_tissue_composite = 0.0492` (`composite_effect = +0.569`, `CI95 = [+0.387, +0.751]`, `p = 0.00547`, `df = 2`), which vanishes in the `no_log_lib` (`df = 3`) sensitivity (`composite_effect = +0.314`, `p = 0.153`, `q_bh_tissue_composite = 0.459`).
- **Matched-null permutation (`q_bh_ts2_all41_empirical`):** Six groups pass `< 0.05` in the primary specification, with contradictory signs across organs and donors:
  - Positive composite (`UP - DOWN` rises with age): `Bladder macrophage` (`+0.569`, `emp_p = 0.00050`, `q_ts2_41_emp = 0.00683`), `Fat CD8+ T` (`+0.183`, `emp_p = 0.00050`, `q_ts2_41_emp = 0.00683`), `Fat fibroblast` (`+0.141`, `emp_p = 0.00150`, `q_ts2_41_emp = 0.0154`; vanishes when dropping `log_lib`: `composite = +0.0167`, `emp_p = 0.699`).
  - Negative composite (`UP - DOWN` falls with age): `Lymph_Node All_Cells` (`-0.197`, `emp_p = 0.00050`, `q_ts2_41_emp = 0.00683`; persists in `no_log_lib` at `-0.199`, `emp_p = 0.00050`, `q_ts2_41_emp = 0.0102`), while `Blood classical monocyte` (`-0.109`, `emp_p = 0.00900`, `q_ts2_41_emp = 0.0738`) and `Spleen B cell` (`-0.109`, `composite_p = 0.0220`, `emp_p = 0.0520`; `no_log_lib` `composite_p = 0.00537`, `emp_p = 0.0430`) trend negative without passing cross-tissue empirical FDR.
- **Interpretation guard:** In these `N = 6` (`df = 2` or `3`) subsets spanning `16 to 28 yr` with `log_lib` collinearity, variance compression of individual gene effects yields tight SEs and discrepant parametric vs random-draw null tails. Combined with opposite composite signs (`Bladder`/`Fat` positive vs `Lymph_Node`/`Blood`/`Spleen` negative) across overlapping donors, TS2 remaining tissues provide **no consistent replication** of the mouse `P1_Clean` signature.

## 6. Limitations

1. Every tissue outside `Blood` (`n=9`) is limited to `N = 6` donors (`45 to 61 yr`, `33 to 61 yr`, or `45 to 69 yr`), yielding `residual_df = 2` in the primary OLS and `df = 3` in the `no_log_lib` sensitivity.
2. Only `9` unique TS2 donors underlie all 5 tissues, with `5` donors shared across all five and identical donor sets between `Lymph_Node` and `Spleen`.
3. `PMC9812260` compartment microbead/FACS normalization and `Blood` preparation mixing prevent compositional inference (Family A and concordance omitted).
4. Enriched/mixed preparations in `Blood` (`TSP10` `Enriched`/`Whole`, `TSP25` `PBMC`, `TSP2` `SheelaPrep`) may still subtly affect pseudobulk expression variance or ambient contamination across donors.

## 7. Verification and File Provenance

- Extraction & raw integer count gate script: `scripts/extract_ts2_remaining_caches.py`
- Analysis script: `scripts/run_followup_ts2_remaining.py`
- Committed outputs:
  - `results/followup-ts2-all28-tissues-audit.csv`
  - `results/followup-ts2-remaining-donors.csv`
  - `results/followup-ts2-remaining-qc.csv`
  - `results/followup-ts2-remaining-genes-summary.csv`
  - `results/followup-ts2-remaining-top-genes.csv`
  - `results/followup-ts2-remaining-p1-signature.csv`
  - `results/followup-ts2-remaining-sensitivities-p1.csv`
  - `results/followup-ts2-remaining-metadata-summary.json`
