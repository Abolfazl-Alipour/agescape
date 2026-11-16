# Full-Depth Redo: Cell-Level Drift, Spread, and Ectopic Identity (Gao Heart + Yazar Blood)

Combined report for the full-depth task: donor-level Part A plus cell-level
Parts 1/2/3. Every verdict is separate per dataset and per part. Every number
below names its source file. For the broader overnight program this task sits
in, see [followup-overnight-summary.md](followup-overnight-summary.md), which
is not duplicated here.

## 1. Design (frozen for this task)

- **Rule 0:** the primary analysis uses full-depth raw counts converted to
  within-cell proportions. `>= 1500` UMI (after the gene mask) is an
  eligibility filter only. 1500/3000-UMI rarefaction appears only as a labeled
  sensitivity, never as the primary.
- **Donor is the unit** for every p-value. Slopes are per decade of age;
  `slope_sd_per_decade` is in SD of the outcome per decade.
- **Technical null:** 5-seed split-half in every part. The split-half SNR
  (`sh_snr_slope`) and the 80%-power minimum detectable effect
  (`mde80_sd_per_decade`) are reported next to every effect.
- **BH families:** `apply_slice_bh` adjusts the per-cell-type non-pooled rows
  within each `[dataset, part, metric, gene_set, depth_mode, model_spec]`
  slice (e.g. Gao Part 1 center_shift full_depth primary: m = 12 cell types).
  Pooled rows carry their own singleton family (m = 1, q = p).
  Stratified-global primary rows are intentionally unadjusted (q = NaN by
  design); their verdict gates use the raw donor-level p-value together with
  the split-half SNR and MDE80 instead.
  Implementation:
  [run_followup_between_subject_gao_yazar.py](../scripts/run_followup_between_subject_gao_yazar.py),
  [run_followup_fulldepth_drift_entropy_identity.py](../scripts/run_followup_fulldepth_drift_entropy_identity.py).
- **Datasets:** Gao heart LV snRNA, 113 donors in the global model, 522,761
  nuclei at or above 1500 UMI; Yazar blood scRNA, 973 donors in the global
  model, 617,760 cells at or above 1500 UMI.
  Pipeline and Part 1/2 commits `ea79fce` (Gao) and `d2723b2` (Yazar); Part 3
  commit `adac77a`.

## 2. Verdicts (one per dataset per part; never collapsed)

| Dataset | Part A (donor level) | Part 1 center_shift | Part 1 scatter | Part 2 | Part 3 E1 | Part 3 E2 |
|---|---|---|---|---|---|---|
| Gao heart | DOES NOT PROVE drift or spread increase (A1 null; A3 slope negative) | DOES NOT PROVE (+0.0198, p 0.278) | DOES NOT PROVE increase (-0.0632, p 0.027, decreases) | DOES NOT PROVE more spread (entropy decreases) | DOES NOT PROVE (count p 0.564) | DOES NOT PROVE (stratum count p 0.658) |
| Yazar blood | DOES NOT PROVE drift or spread increase (A1/A3 slopes negative) | DOES NOT PROVE (-0.0028, p 0.384) | DOES NOT PROVE increase (-0.0169, p 1.5e-4, decreases) | DOES NOT PROVE more spread (entropy decreases) | DOES NOT PROVE (count p 0.0825) | DOES NOT PROVE (stratum count p 0.456) |

Machine-readable verdicts with reasons:
[followup-fulldepth-qc.json](../results/followup-fulldepth-qc.json)
(`/verdicts/gao`, `/verdicts/yazar`). Example reasons: Gao Part 1
center_shift "Global center_shift fails PROVES (slope=0.0198, p=0.278,
MDE80=0.051 [high_power<=0.25])"; Gao Part 3 E2 "E2 fails randomness PROVES
(cnt p=0.658, SNR=7.30, MDE80=0.114<=0.5)".

## 3. Part A: donor-level drift, variance, dispersion

Source: [followup-between-subject-gao-models.csv](../results/followup-between-subject-gao-models.csv),
[followup-between-subject-yazar-models.csv](../results/followup-between-subject-yazar-models.csv),
[followup-between-subject-variance-gao.csv](../results/followup-between-subject-variance-gao.csv),
[followup-between-subject-variance-yazar.csv](../results/followup-between-subject-variance-yazar.csv).
Rows below are `model_spec = stratified_global_primary`.

### 3.1 Gao

- A1 hellinger_young_ref, full_1500: +0.0234 SD/dec, p 0.238, MDE80 0.056.
  Rarefied depths agree on null (1500: +0.0094, p 0.792; 500: -0.0206,
  p 0.610; 3000: +0.0105, p 0.828). A1 hellinger_mature_ref, full_1500:
  +0.0163, p 0.647, MDE80 0.101.
- A2 OLD_vs_YOUNG, batch_residualized arm: mean_eta2 0.027,
  p_eta2_perm 0.0075, but median log-var-ratio +0.093 with two-sided
  p 0.738; the depth_only arm is null (p_eta2_perm 0.394) and MID_vs_YOUNG
  is null in both arms. Inconsistent across arms and contrasts: no proof of
  a variance increase.
- A3 within_donor_loo_disp: 1500-depth -0.1378, p 0.0070, MDE80 0.142;
  full_1500 -0.0661, p 0.0224, MDE80 0.081. Dispersion goes DOWN with age,
  so A3 does not prove an increase either.

### 3.2 Yazar

- A1 hellinger_young_ref, full_1500: -0.0221, p 1.0e-6, MDE80 0.0126
  (cells sit CLOSER to the young reference with age; all rarefied depths
  agree on the negative sign). A1 hellinger_mature_ref, full_1500: -0.0273,
  p 8.2e-5.
- A2 OLD_vs_YOUNG, batch_residualized: mean_eta2 0.0064, p_eta2_perm 0.0005,
  but median log-var-ratio -0.17, two-sided p 0.0010: variance DECREASES in
  OLD. No proof of a variance increase.
- A3 within_donor_loo_disp: 1500-depth -0.1219, p 1.9e-18; full_1500
  -0.0705, p 2.5e-17. Dispersion decreases with age.

## 4. Full-depth Part 1: center_shift and scatter

Source: [followup-fulldepth-gao-cell-drift.csv](../results/followup-fulldepth-gao-cell-drift.csv),
[followup-fulldepth-yazar-cell-drift.csv](../results/followup-fulldepth-yazar-cell-drift.csv).
Rows below are stratified_global_primary, full_depth.

### 4.1 Gao

- center_shift: +0.0198 SD/dec, p 0.278, SNR 1.72, MDE80 0.051.
  DOES NOT PROVE.
- scatter: -0.0632, p 0.027, SNR 5.68, MDE80 0.080. Nominally significant
  but NEGATIVE: cells are less scattered around their donor centroid with
  age. DOES NOT PROVE an increase.
- Sensitivity signals reported honestly although the primary is null:
  equal_ref_size +0.0500, p 0.0039; pythagorean_centroid_h2 +0.0682,
  p 0.0419; spline_resid_full +0.0483, p 0.199. These alternative
  reference constructions move but the pre-registered LOO-reference primary
  does not support them.

### 4.2 Yazar

- center_shift: -0.0028, p 0.384, SNR 0.59, MDE80 0.009. DOES NOT PROVE.
- scatter: -0.0169, p 1.48e-4, SNR 3.55, MDE80 0.012. Decreases with age.
  DOES NOT PROVE an increase.
- Composition note: the pooled All_Cells center_shift slope is positive
  (raw +0.000422, p 4.4e-6; SD +0.0370) while the within-cell-type
  stratified slope is null (-0.0028, p 0.384). The pooled sign is a
  cell-type composition effect, not drift within cell types.

## 5. Full-depth Part 2: entropy / spread

Source: [followup-fulldepth-gao-spread.csv](../results/followup-fulldepth-gao-spread.csv),
[followup-fulldepth-yazar-spread.csv](../results/followup-fulldepth-yazar-spread.csv).
Rows below are stratified_global_primary, full_depth.

### 5.1 Gao

- shannon -0.0998, p 0.0128, SNR 8.40, MDE80 0.112;
  deman_kl_methods +0.0998, same p (exact mirror by construction);
  deman_kl_sqrt +0.0382, p 0.0183, SNR 3.07, MDE80 0.045.
  Cells become MORE concentrated (lower entropy) with age.
  DOES NOT PROVE more spread.

### 5.2 Yazar

- shannon -0.0222, p 0.0079, SNR 4.59, MDE80 0.023;
  deman_kl_methods +0.0222, same p; deman_kl_sqrt +0.0199, p 1.06e-4,
  SNR 4.12, MDE80 0.014. Same direction as Gao: lower entropy with age.
  DOES NOT PROVE more spread.
- Composition note again: pooled All_Cells shannon is positive (raw
  +0.0035, p 0.0068) while the stratified slope is negative (-0.0222,
  p 0.0079). Same composition reading as Part 1.

## 6. Full-depth Part 3: ectopic identity (E1 off-program, E2 silent)

Source: [followup-fulldepth-gao-ectopic.csv](../results/followup-fulldepth-gao-ectopic.csv),
[followup-fulldepth-yazar-ectopic.csv](../results/followup-fulldepth-yazar-ectopic.csv).
Reference sets: Gao 12 E1 cell types (median 13.5 markers), 15,930 genes
kept; Yazar 14 E1 cell types (median 7.5 markers), 35,307 genes kept
([followup-fulldepth-qc.json](../results/followup-fulldepth-qc.json)).
Rows below are stratified_global_primary unless noted.

### 6.1 Gao

- E1 mean_ectopic_count full_depth: +0.0307, p 0.564, SNR 1.19,
  MDE80 0.150. Sensitivities: rarefied_1500 +0.0121, p 0.806;
  spline_resid +0.0293, p 0.616; E1_stringent +0.0323, p 0.541.
  E1 jaccard_both_pos +0.0408, p 0.355, SNR 3.95.
  E1 singleton_20cells -0.0257, p 0.431; E1 excess_jaccard +0.0318,
  p 0.216. DOES NOT PROVE (randomness or otherwise).
- E2_global mean_ectopic_count full_depth: -0.0301, p 0.226, SNR 1.25
  (tagged `sensitivity_global_depth_vuln`: the global silent set is
  depth-sensitive by construction, see 6.3). E2_global rarefied_1500:
  -0.0635, p 0.108. E2_global jaccard: +0.0066, p 0.869, SNR 0.57.
- E2_stratum (bounded depth strata, depth-matched primary):
  count -0.0179, p 0.658, SNR 7.30, MDE80 0.114; jaccard +0.0182,
  p 0.683, SNR 7.91. DOES NOT PROVE.

### 6.2 Yazar

- E1 mean_ectopic_count full_depth: -0.0121, p 0.0825, SNR 1.26,
  MDE80 0.019. Sensitivities: rarefied_1500 +0.00066, p 0.918;
  spline_resid -0.0114, p 0.122; E1_stringent -0.0138, p 0.0465
  (nominally significant but NEGATIVE: fewer ectopic UMIs with age).
  E1 singleton_20cells +0.0219, p 0.0078 (positive, but the count gate
  fails, so no randomness proof). E1 jaccard -0.0057, p 0.544;
  E1 excess_jaccard +0.0022, p 0.820. DOES NOT PROVE.
- E2_global mean_ectopic_count full_depth: -0.0210, p 0.0056, SNR 2.23;
  rarefied_1500 -0.0292, p 0.0017 (negative: fewer silent-gene UMIs with
  age). E2_global jaccard +0.0196, p 0.0287, SNR 3.92.
- E2_stratum: count -0.0066, p 0.456, SNR 1.88, MDE80 0.025; jaccard
  +0.0281, p 0.0070, SNR 6.23; singleton_20cells -0.0232, p 0.0220.
  Mixed signs with the count gate failed: DOES NOT PROVE.

### 6.3 Caveats that stay attached to Part 3

- E2_global is depth-vulnerable: silent genes are defined by detection
  rates, and detection depends on depth. The depth-matched E2_stratum rows
  are the primary E2 evidence; the E2_global rows are labeled
  `sensitivity_global_depth_vuln` in `analysis_role`.
- Rarefaction bug (found and fixed before Commit C, `adac77a`):
  `multinomial_downsample` returns `matrix.copy()` with new data, so genes
  dropping to zero stay as explicit stored zeros. Part 3 binarizes by
  stored entries (`np.diff(indptr)`, `data = 1`), so every Part 3
  `rarefied_1500` row was identical to its full-depth row (e.g. Gao E2
  jaccard 0.006619 / p 0.869412 in both; Gao E2 count -0.030072 /
  p 0.225889 in both), and the seed-0..4 YOUNG-1500 detection references
  were affected the same way. Fix: `eliminate_zeros()` on all five
  downsample outputs inside this pipeline's script; the shared sampler was
  left untouched because committed results depend on it. Full-depth
  primaries never had explicit zeros (checked: Gao shard 0 nnz 95,115,289
  with 0 explicit zeros; Yazar shard 0 nnz 42,915,505 with 0 explicit
  zeros), and the value-weighted Part 1/2 rarefied paths were unaffected:
  regenerated Part 1/2 CSVs are byte-identical to the pushed versions.
  Post-fix the rarefied rows differ from full depth (e.g. Gao rarefied
  E2 jaccard +0.0117, p 0.795; Gao E2_global rarefied count -0.0635,
  p 0.108).

## 7. Robustness

[followup-fulldepth-robustness.csv](../results/followup-fulldepth-robustness.csv):
159 rows, 152 pass `passes_robustness`.

## 8. Judgment-call log (decision, reasoning, rejected alternative)

| # | Decision | Reasoning | Rejected alternative |
|---|---|---|---|
| 1 | Full-depth within-cell proportions as primary; 1500 UMI as eligibility only | User Rule 0; rarefaction discards most UMIs and adds sampling noise | 1500-UMI rarefied primary (old Part B design, not run) |
| 2 | Donor as the unit for every p-value | Cells within a donor are not independent | Cell-level N (would inflate significance) |
| 3 | BH on per-cell-type rows per slice; pooled singleton families; stratified-global rows unadjusted with p + SNR + MDE80 gates | The stratified global is one pre-registered contrast, not a family | Global BH over everything, or BH on the single primary |
| 4 | 5-seed split-half null with SNR and MDE80 next to every effect | Raw p-values alone cannot separate signal from sampling noise | p-value-only reporting |
| 5 | Verdicts separate per dataset and per part | Gao and Yazar disagree in places; parts test different claims | One collapsed verdict |
| 6 | Fix explicit zeros at this pipeline's call sites; leave the shared sampler untouched | Other committed results depend on the shared function | Editing the shared sampler and re-validating everything downstream |
| 7 | Report Gao equal_ref_size / pythagorean sensitivities although the primary is null | They are real outputs of pre-registered sensitivities; hiding them would be selective | Reporting only the null primary |
| 8 | Read Yazar pooled All_Cells positives as composition effects | Within-cell-type stratified slopes are null or negative in the same models | Claiming drift/spread from the pooled slope |
| 9 | E2_stratum as primary E2 evidence; E2_global kept as depth-vulnerable sensitivity | Silent-gene detection depends on depth by construction | Treating E2_global as depth-proof |
| 10 | No new datasets, access applications, Kedlian retry, or GTEx | Explicit user scope rule | Expanding scope to chase power |
| 11 | No edits to `report/plos-one/`, `.tex`, `.pdf`, or figure files | Manuscript freeze; verified before every push | Editing frozen files |

## 9. Reproduction

```
PY=<agescape-conda-env>/bin/python
$PY scripts/resource_guard.py run --path <revision-data-root> --min-free-gb 30 --memory-gb 12 -- \
    $PY scripts/run_followup_fulldepth_drift_entropy_identity.py --dataset all --stage parts12
$PY scripts/run_followup_fulldepth_drift_entropy_identity.py --dataset all --stage part3 --reuse-pass1
<agescape-conda-env>/bin/pytest -p no:cacheprovider tests/test_followup_fulldepth.py -v
```

Timings observed: Gao Pass 1 589 s + Pass 2 490 s; Yazar Pass 1 308 s +
Pass 2 ~380 s; Part 3 both datasets 1668 s (peak RSS 7.5 GB).
Tests: 9 passed ([test_followup_fulldepth.py](../tests/test_followup_fulldepth.py)),
including the explicit-zero cleanup regression test.
