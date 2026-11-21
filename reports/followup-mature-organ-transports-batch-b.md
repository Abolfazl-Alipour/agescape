# Item 5 Batch B: HRCA retina mature (30-plus) transport rerun

Status: COMPLETE. All 9 HRCA suites rerun at 30-plus with mature P1m and P0m. Nothing was submitted. No manuscript file was touched.

## Run record

- Cohort: Human Retina Cell Atlas (HRCA), Li et al. Nat Genet 2026.
- Script: scripts/run_followup_mature_organ_transports_b_hrca.py (860 lines, sha256 35e0dd67235a77127b8fd4901661e7c8dc8611a8dff597dcfc8b82eec8569b35). The exact verified bytes were transferred to beast and the hash was checked before running. The script was executed directly on beast by Muse Spark (direct execution, no coding-agent CLI), under Reza's explicit authorization.
- Protocol: frozen Batch A rules. Pseudobulk log2(CPM+1), technical-gene removal, CPM at least 1 in at least 50 percent of mature donors, original covariate models, 2000-draw expression-decile-matched null with seed 20261105, 1000-permutation concordance with seed 20261105. Donors aged at least 30.0 years were restricted BEFORE pseudobulk aggregation and detection filtering in the snRNA suites. Donor is the inferential unit, never the cell.
- Runtime: 201.5 seconds total, peak RSS 1199.85 MB. All 8 original input files SHA256-verified unchanged before and after.
- Outputs: results/followup-mature-organ-transports-b-hrca.csv (412 rows), results/followup-mature-organ-transports-b-hrca-bh-check.csv (539 rows), results/followup-mature-organ-transports-b-hrca-runtime.csv (37 rows).

## Suites and donors

| Suite | Donors original to mature | Age range mature |
|---|---|---|
| Sn_Adult_MacularOnly | 86 to 80 (6 excluded) | 30 to 90 |
| Sn_Adult_PeripheralOnly | 50 to 48 (2 excluded) | 30 to 90 |
| Sn_Adult_Unenriched | 86 to 80 (6 excluded) | 30 to 90 |
| Sn_Adult_ExactAge | 81 to 75 (6 excluded) | 30 to 88 |
| Sn_ChenAncestry_Adult | 51 to 45 (6 excluded) | 32 to 90 |
| Sn_AllAges_Unenriched | 93 to 80 (13 excluded) | 30 to 90 |
| Sc15_StudyAdj | 15 to 15 (0 excluded) | 42 to 89 |
| Sc9_External_StudyAdj | 9 to 9 (0 excluded) | 42 to 89 |
| Sc5_External_NoDisease | 5 to 5 (0 excluded) | 42 to 86 |

The scRNA suites needed no age restriction: all their donors were already 42-plus. The snRNA suites are the ones rerun at 30-plus.

## Multiple-testing families declared before correction

Per suite, BH correction was applied within each of these families separately, never pooled across suites:

- Fractions: one family per suite over the tested cell types (9, 9, 9, 9, 9, 9, 8, 6, 4 types).
- Genes global: one family per suite over all tested genes (145,470; 145,807; 148,051; 148,001; 146,881; 148,051; 114,446; 89,559; 67,446 tests; 1,153,712 tests in total).
- Transport composite: one family per suite and tier (P1m and P0m separately) over the tested cell types.
- Transport empirical: one family per suite and tier over the tested cell types.
- Concordance: one 1000-permutation test per suite, no BH family (single test each).

## Headline: no mature transport pass anywhere in retina

Transport was tested in 144 rows (72 P1m, 72 P0m) across the 9 suites. Both-pass count is 0 of 144 under either signature in every suite. Of the 144 rows, 66 are EMPIRICAL_ONLY (the matched-null test passes but the composite test does not) and 78 are FAIL_BOTH.

All 4 original P1m both-pass units DIED. All four were negative-direction originally, and all four now sit at EMPIRICAL_ONLY with the composite effect flipped to positive:

| Suite | Cell type | Original effect | Mature effect | Mature composite q | Mature empirical q | Mature class |
|---|---|---|---|---|---|---|
| Sn_Adult_MacularOnly | MG (Muller glia) | -0.0492 | +0.0114 | 0.303 | 0.0101 | EMPIRICAL_ONLY |
| Sn_Adult_MacularOnly | RGC (retinal ganglion cell) | -0.0379 | +0.0149 | 0.118 | 0.0022 | EMPIRICAL_ONLY |
| Sn_ChenAncestry_Adult | RGC | -0.0487 | +0.0118 | 0.518 | 0.0090 | EMPIRICAL_ONLY |
| Sn_AllAges_Unenriched | RGC | -0.0450 | +0.0129 | 0.245 | 0.0015 | EMPIRICAL_ONLY |

No NEW both-pass unit appears under either signature. The mature retina result is a clean null for transport.

## Fractions: one signal survives

The AC (amacrine cell) fraction signal survives the mature restriction in 3 snRNA suites:

| Suite | Mature q | Original q |
|---|---|---|
| Sn_Adult_Unenriched | 0.00143 | 0.00236 |
| Sn_Adult_ExactAge | 0.000505 | 0.00139 |
| Sn_AllAges_Unenriched | 0.00143 | 0.00612 |

No other fraction is significant at mature q below 0.05 in any suite. One fraction is UNTESTABLE: RPE (retinal pigment epithelium) in Sn_Adult_MacularOnly, with 13 donors present against a 15-donor floor.

## Genes: near-total collapse of the original signal

Across all 9 suites, 1,153,712 gene tests were run and only 3 genes are significant at global BH q below 0.05:

| Suite | Cell type | Gene | Mature q | Original q |
|---|---|---|---|---|
| Sn_Adult_MacularOnly | RGC | MAPK6 | 0.0471 | 0.0434 |
| Sn_Adult_MacularOnly | RGC | ADORA1 | 0.0471 | not recorded (no original record) |
| Sn_ChenAncestry_Adult | RGC | MAP4K3-DT | 0.00832 | 0.0252 |

At the cell-type level, only RGC is significant: Sn_Adult_MacularOnly (2 significant genes, q 0.0471) and Sn_ChenAncestry_Adult (1 significant gene, q 0.00832). The original gene signal collapses almost everywhere: Sn_Adult_MacularOnly 74 significant genes originally to 2 mature; Sn_AllAges_Unenriched 1500 to 0; Sn_Adult_Unenriched 4 to 0; Sn_Adult_ExactAge 4 to 0; Sn_ChenAncestry_Adult 5 to 1.

## Concordance: null in all 9 suites

Mouse-human gene-effect concordance is null in every suite (mean r 0.054 to 0.301, all two-sided p above 0.15).

## Cross-check: the two Unenriched suites agree exactly

Sn_AllAges_Unenriched and Sn_Adult_Unenriched have identical 30-plus donor sets by construction. Their mature transport rows are identical in keys and in values (composite q, empirical q, class), and their AC fraction q values are identical (0.00143). The rerun is internally consistent.

## Verdict accounting

Across all 412 rows: NULL 189, DIED 142, SURVIVED 7, UNTESTABLE 1, blank 73. The 7 SURVIVED rows are the 3 AC fractions, the 2 RGC gene cell-type rows, and the 2 RGC gene hits (MAPK6, MAP4K3-DT) listed above. The 73 blank verdicts are by design, not missing data: 72 are the P0m transport rows, which carry no verdict because no original P0 baseline exists (the original organ-10 HRCA run used P1 only; the verdict rule is defined on the P1m arm), and 1 is the ADORA1 gene hit, which has no original record to compare against.

## BH verification

- Fractions and transport families were independently recomputed from the stored raw p-values with a fresh BH implementation (not imported from the run script): 360 of 360 checked values match, 0 mismatches.
- The gene families were verified at the code level: the script computes q over the full per-suite gene p-value vector (up to 148,051 tests) and asserts the recomputed q matches the stored q within 1e-12 at runtime; the run passed that assert in every suite. The BH-check CSV lists only the per-cell-type gene rows, not the full vector, so a subset recompute over the listed rows is the wrong family and is not used.
- The script's own explicit BH versus statsmodels comparison gives max absolute difference 1.11e-16 across all 539 BH-check rows.

## Notes and limits

- The scRNA suites are small (15, 9, and 5 donors). Their transport nulls carry little power; a null there is weak evidence.
- The P0m transport verdicts are blank by design (see above). P0m both-pass counts are still reported: 0 of 72.
- The ADORA1 gene hit is mature-significant with no original record; it is reported as a fact, not as a SURVIVED verdict.
- Nothing under report/ was touched. No manuscript file was touched. Nothing was submitted.
