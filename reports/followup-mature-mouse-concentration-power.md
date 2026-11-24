# Mouse concentration-scatter detection floor

## Source and purpose

This report files the repository source for the mouse 18-month-vs-30-month scatter detection floor of **0.148 standard deviations (SD) per decade**. It reproduces section 6, "Power for the mouse 18m-vs-30m comparison," of the November 23 feasibility audit at `/tmp/kimi-mouse-report.txt` on Beast. The audit is filed here because a file under `/tmp` may not survive a Beast reboot.

Run `python3 scripts/compute_mature_mouse_scatter_mde80.py` from the repository root. The script uses only Python's `csv` and `math` modules and reads the human anchors directly from:

- `results/followup-mature-fulldepth-yazar-cell-drift.csv`
- `results/followup-mature-fulldepth-gao-cell-drift.csv`

## Donor-level design

The mouse comparison has 4 donors at 18 months and 4 at 30 months, for 8 donors total:

- 18 months: `18-F-50`, `18-F-51`, `18-M-52`, `18-M-53` (2 female, 2 male).
- 30 months: `30-M-2`, `30-M-3`, `30-M-4`, `30-M-5` (all male).

Thus, age is confounded with sex in this comparison: the 30-month group is entirely male. More cells per mouse would not remove the donor-level power limitation.

## Human anchor and scaling

The Yazar full-depth primary scatter anchor is the row with `model_spec=stratified_global_primary`, `group=Stratified_CellType_Global`, and `is_pooled=True`:

- slope: `-0.01966524117551119` SD per decade;
- SE: `0.004909531960818628` SD per decade (`0.0049` rounded);
- `p = 6.685230715256798e-05`;
- 925 donors in the fit; 932 total mature Yazar donors used by the audit for scaling.

The audit scales the human donor-level SE to the 8-mouse design as

`SE_mouse = SE_Yazar * sqrt(932 / 8)`.

For 80% power at two-sided alpha 0.05, the normal-approximation multiplier is `z(0.975) + z(0.80) ≈ 1.96 + 0.84 = 2.8`, so `MDE80 = 2.8 * SE_mouse`.

The three requested calculations are:

| Variant | SE_mouse (SD/decade) | MDE80 arithmetic | MDE80 (SD/decade) | Rounded |
|---|---:|---|---:|---:|
| Rounded SE, n=932 | `0.0049 * sqrt(932/8) = 0.0528882312` | `2.8 * 0.0528882312` | `0.1480870474` | **0.148** |
| Exact SE, n=932 (headline) | `0.004909531960818628 * sqrt(932/8) = 0.0529911146` | `2.8 * 0.0529911146` | `0.1483751208` | **0.148** |
| Exact SE, n=925 | `0.004909531960818628 * sqrt(925/8) = 0.0527917385` | `2.8 * 0.0527917385` | `0.1478168679` | **0.148** |

The filed value is therefore **0.148 SD per decade**. Using 925 fitted donors instead of 932 total mature donors does not change the value rounded to three decimals.

## Human effects relative to the floor

The Gao full-depth primary scatter anchor has slope `-0.07228251880607399` SD per decade, `p = 0.041028425872898774`, and 105 donors. Using the audit's rounded effect sizes and rounded-SE MDE80:

- Gao: `0.0723 / 0.1480870474 = 0.4882263593`, or **48.8%** of the mouse MDE80.
- Yazar: `0.0197 / 0.1480870474 = 0.1330298655`, or **13.3%** of the mouse MDE80.

Neither a Gao-sized nor a Yazar-sized human scatter effect would be detectable with 80% power in the 4-vs-4 mouse comparison under this scaling approximation.

## Small-sample sensitivity, not the filed value

As a sensitivity only, replacing the normal multiplier with `t(0.975, df=6) + t(0.80, df=6) = 2.4469 + 0.9057 = 3.3526` gives

`3.3526 * 0.0529911146 = 0.1776580107` SD per decade, approximately **0.1777** (about **0.178**).

This is **not** the 0.148 value filed above. The normal-approximation floor is optimistic; the exact small-sample floor is larger. The sensitivity therefore strengthens the conclusion that neither human effect is detectable in the mouse comparison.
