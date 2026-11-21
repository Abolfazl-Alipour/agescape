# Agescape mature-to-old rebuild: Item 5 Batch C report

TS2 plus the remaining small mature cohorts: Yan retina, RPE bulk, ovary
postmenopausal, PBMC, Lake kidney, Perez muscle. Transport rerun with the
mature signatures P1m and P0m at 30-plus, seed 20261105, 2000-draw
expression-matched null. Script:
`scripts/run_followup_mature_organ_transports_c.py`.

## Headline

237 rows total across 7 cohort CSVs. Both-pass (composite q below 0.05 AND
empirical q below 0.05) under the mature signatures: 3 rows out of 228
tested transport rows.

- RPE bulk: the one original both-pass transport SURVIVED under both P1m and
  P0m (negative direction, same sign as the original).
- Lake kidney scCv3 endothelial cell: one NEW both-pass under P0m only.
- Everything else: 0 both-pass under both tiers.

Every original both-pass transport in this batch DIED except the RPE one:
all 11 original Thomson PBMC both-pass rows died under both tiers, and the
Perez muscle gene signal died (0 of 18,058 genes significant, original 8 of
20,343). The empirical-only pattern is widespread (see the interpretation
note below): many rows pass the matched-null test while the parametric
composite test does not. Empirical-only is weak evidence on its own; the
headline verdict rests on the composite test.

## Per-cohort results

### TS2: `results/followup-mature-organ-transports-c-ts2.csv` (84 rows)

- Remaining tissues (Blood non-NK 7, Bladder 9, Fat 3, Lymph_Node 11,
  Spleen 11; one pooled 41-group BH family per tier): 0 both-pass under
  P1m and 0 under P0m. Empirical-only 18 of 41 under P1m, 19 of 41 under
  P0m. No original composite pass existed, so no verdict of DIED applies
  (all verdicts blank, as the protocol prescribes when no original q is
  below 0.05).
- Of the 4 original empirical-only groups (41-family, verified against
  `results/followup-ts2-remaining-p1-signature.csv`): Bladder macrophage
  and Fat CD8-positive alpha-beta T cell stay EMPIRICAL_ONLY under both
  tiers; Fat fibroblast and Lymph_Node All_Cells die to FAIL_BOTH under
  both tiers.
- Blood NK (single-endpoint family, n 7 donors, df 2): 0 both-pass under
  both tiers. The original empirical-only status DIED: original empirical
  q 0.00995, mature empirical q 0.1749 under P1m and 0.4688 under P0m,
  class EMPIRICAL_ONLY to FAIL_BOTH. Hellinger and off-program NK tests
  stay UNTESTABLE (the 30-to-39 reference has too few NK donors); no
  vector result was fabricated.

### Yan retina: `results/followup-mature-organ-transports-c-yan.csv` (14 rows)

- 7 groups, one 7-group BH family per tier: 0 both-pass under P1m and 0
  under P0m. Empirical-only 2 of 7 under each tier (Muller glia and RGC,
  both negative effects: P1m -0.0647 and -0.0533 per decade). No original
  q below 0.05 existed, so all verdicts are blank.

### RPE bulk: `results/followup-mature-organ-transports-c-rpe.csv` (2 rows)

- Single-endpoint family, 13 donors, batch-adjusted primary arm.
  Both-pass 1 of 1 under both tiers. Verdict SURVIVED under both tiers.
- P1m: original effect -0.0969 per decade, mature effect -0.0211 per
  decade, mature q 0.0443, mature empirical q 0.0005. P0m: mature effect
  -0.0194 per decade, mature q 0.0439, mature empirical q 0.0005.
  Detection: 933 up and 651 down P1m genes detected; 1063 up and 791 down
  P0m genes detected.
- Donor table came from a live GEO fetch (the frozen fallback file was
  not needed).

### Ovary postmenopausal: `results/followup-mature-organ-transports-c-ovary.csv` (18 rows)

- 9 cell types, one 9-group BH family per tier: 0 both-pass under P1m
  and 0 under P0m. Empirical-only 8 of 9 under each tier (P1m except
  Schwann, empirical q 0.0530; P0m except Macrophage, empirical q 0.1659).
- The 2 original empirical-only rows (ArtEnd and Stromal 1) stay
  EMPIRICAL_ONLY under both tiers.
- All 9 mature effects are negative under both tiers (-0.005 to -0.025
  per decade), opposite in sign to the mostly positive original P1
  effects. The parametric composite test does not pass anywhere
  (mature q 0.36 to 0.64), so this is a matched-null-only pattern, not a
  both-pass transport.
- One row was independently recomputed by hand from the cache: ArtEnd
  P1m effect -0.015253, SE 0.015690, p 0.3426, df 20, matching the CSV
  exactly. Detection 930 to 987 up and 635 to 689 down P1m genes per
  cell type.

### PBMC: `results/followup-mature-organ-transports-c-pbmc.csv` (72 rows)

- Binary mature design 30 to 39 vs 60-plus per the audit, not a
  continuous model. Per-cohort BH families (23 Zhu, 13 Thomson) per tier
  per test type.
- Both-pass 0 of 36 under P1m and 0 of 36 under P0m.
- Thomson (GSE214546): all 11 original both-pass rows DIED under both
  tiers (P1m mature q 0.8657 to 0.9503; mature effects collapse from
  about +0.08 to a mean of +0.0033, range -0.034 to +0.024). Platelets
  and pre-B cell had no original q below 0.05, so their verdicts are
  blank. Empirical-only 10 of 13 under P1m, 7 of 13 under P0m.
- Zhu (GSE213516): 0 original both-pass. The 4 original empirical-only
  rows (Memory-B-IL32, NK-SELL, cDC, pDC) stay EMPIRICAL_ONLY under both
  tiers. Empirical-only 22 of 23 under P1m, 23 of 23 under P0m.
- Cross-check: the other agent's independent mature PBMC P1m
  implementation (read-only file
  `results/followup-mature-pbmc-GSE214546-p1m-models.csv`) matches
  exactly for the rows compared, for example B cell progenitor effect
  0.016392, q 0.865666, same design and same donor counts.

### Lake kidney: `results/followup-mature-organ-transports-c-lake.csv` (44 rows)

- Per-assay BH families (9 scCv3, 10 snCv3) per tier per test type.
  Effects are PER YEAR, not per decade, for both original and mature.
- Both-pass: P1m 0, P0m 1. The one both-pass is scCv3 endothelial cell
  under P0m, verdict NEW: original q 0.5040, mature q 0.0402, mature
  empirical q 0.000643, mature effect +0.002081 per year (SE 0.000623,
  null SD 0.000458), original effect +0.007743 per year.
- Empirical-only: scCv3 9 of 9 under P1m, 8 of 9 under P0m; snCv3 7 of
  10 under P1m, 8 of 10 under P0m.
- Of the 6 original empirical-only rows, 5 stay EMPIRICAL_ONLY under
  P1m and the snCv3 loop of Henle thick ascending limb row dies to
  FAIL_BOTH under both tiers; under P0m the scCv3 endothelial row
  upgrades to BOTH_PASS and the other 4 stay EMPIRICAL_ONLY.
- 3 snCv3 cell types are UNTESTABLE under both tiers (T cell, kidney
  connecting tubule epithelial cell, mononuclear phagocyte): 5 donors
  after the 30-plus restriction, below the 6-donor floor. The snCv3
  primary included one 29-year-old donor who was excluded.

### Perez muscle: `results/followup-mature-organ-transports-c-muscle.csv` (3 rows)

- Continuous 30-plus whole-nucleus gene model, 11 donors aged 70 to 90,
  design intercept, age per 10 years, mean log UMI, residual df 8.
  18,058 genes tested, 0 significant at q below 0.05 (original: 8 of
  20,343). Top gene AC016588.2, effect +0.7319, p 4.51e-05, q 0.814.
  Verdict DIED.
- Transport P1m and P0m are UNTESTABLE (blocked, no released
  annotations), as the audit prescribed. No whole-nucleus transport
  claim is made.

## Multiple-testing families (declared before BH)

- TS2 remaining: one pooled 41-group family per tier, composite and
  empirical separately. TS2 NK: single-endpoint family per tier.
- Yan: one 7-group family per tier. RPE: single-endpoint family per
  tier. Ovary: one 9-group family per tier.
- PBMC: per-cohort families (23 Zhu, 13 Thomson) per tier.
- Lake: per-assay families (9 scCv3, 10 snCv3) per tier.
- Perez muscle genes: one 18,058-gene family.
- Every family was corrected with Batch A `bh_checked` (explicit BH
  cross-checked against statsmodels). Independent BH recomputation:
  `results/followup-mature-organ-transports-c-bh-check.csv`, 456 rows,
  max absolute diff 1.110e-16, and every cohort CSV mature q was
  cross-checked against a BH-check row with matching raw p.

## Interpretation note

The empirical-only pattern is systematic in this batch: 18 of 41 TS2
groups, 22 of 23 Zhu cell types, 8 of 9 ovary types, and most Lake rows
pass the matched-null test while the parametric composite test passes
almost nowhere. This follows from the mature signatures being large
(P1m 1,724 human-mapped genes, P0m 2,049, versus 118 in P1): averaging
many genes shrinks the matched-null spread, so modest age slopes look
large against the null even when donor-level noise (the composite SE)
does not support them. Empirical-only rows are therefore reported as
EMPIRICAL_ONLY, not as transport, and no both-pass claim rests on them.
The RPE survival and the Lake endothelial P0m pass are the only
composite-supported transports in this batch.

## Deviations and method notes

- PBMC used the binary mature design 30 to 39 vs 60-plus prescribed by
  the audit (section 9), preserving the original strict-group estimand,
  instead of a continuous age model.
- Lake effects are per year, as in the original, not per decade.
- PBMC program detection used a detectable mask (unique symbols only);
  the original excluded all P1 genes from the null pool instead. The
  null pool here keeps every matrix row. This is a stated deviation from
  the original PBMC null construction.
- RPE donor table was fetched live from GEO; the frozen donor file was
  not needed. RPE counts came from the downloaded GEO tar.
- TS2 CELLxGENE h5ad files were deleted after cache extraction, as the
  original pipeline does (delete_after is True); the extracted caches
  are retained under `tmp_organs/ts2_rem/`.
- Yan and RPE direct gene models were already mature; only transport
  was rerun, as tasked. Ovary direct genes and CLR were already mature;
  only transport was rerun.

## Bugs found and fixed during this run

- Lake runner crashed on `cache["sex"][order].astype(str).str.lower()`
  (numpy arrays have no `.str` accessor). Fixed to wrap in a pandas
  Series first. The Lake cohort was rerun after the fix; no other cohort
  used that pattern.
- BH-check CSV lost P1m rows for ovary, PBMC, Lake, and TS2 NK: those
  runners used tier-less family keys, and the incremental append deduped
  on (cohort, family, unit), so P0m rows overwrote P1m rows (326 rows
  instead of 456). The cohort CSVs were never affected (their q values
  were computed in memory before the append). Fix: family keys now
  carry the tier, and finalize rebuilds the BH-check CSV from the
  cohort CSVs with a fresh BH recomputation that asserts equality with
  every CSV q value. The rebuilt file has 456 rows and max diff
  1.110e-16.
- Pre-existing repo bug, not fixed here (out of scope for this batch):
  `scripts/extract_ts2_remaining_caches.py` calls
  `download_with_guard` with an argparse Namespace that lacks
  `connect_timeout`, while `scripts/resource_guard.py` reads
  `args.connect_timeout`, raising AttributeError. The four TS2 tissue
  downloads were done with plain curl plus size verification instead.

## Runtime

`results/followup-mature-organ-transports-c-runtime.csv`: 26 rows, 455.5
seconds total across all cohorts. Peak RSS 8,820 MB (Yan retina load),
under the 12 GiB cap. TS2 took 149.8 seconds, Yan 59.6, PBMC 58.7,
Lake 42.0, ovary 35.3, RPE 5.3, muscle 5.1.

## Files

- `scripts/run_followup_mature_organ_transports_c.py`
- `results/followup-mature-organ-transports-c-ts2.csv` (84 rows)
- `results/followup-mature-organ-transports-c-yan.csv` (14 rows)
- `results/followup-mature-organ-transports-c-rpe.csv` (2 rows)
- `results/followup-mature-organ-transports-c-ovary.csv` (18 rows)
- `results/followup-mature-organ-transports-c-pbmc.csv` (72 rows)
- `results/followup-mature-organ-transports-c-lake.csv` (44 rows)
- `results/followup-mature-organ-transports-c-muscle.csv` (3 rows)
- `results/followup-mature-organ-transports-c-bh-check.csv` (456 rows)
- `results/followup-mature-organ-transports-c-runtime.csv` (26 rows)
- `reports/followup-mature-organ-transports-batch-c.md` (this file)

Nothing under `report/` and no manuscript file was touched. Nothing was
submitted anywhere.
