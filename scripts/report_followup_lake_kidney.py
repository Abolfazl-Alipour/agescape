#!/usr/bin/env python3
"""Render the restricted Lake kidney report from committed aggregate CSVs."""

from pathlib import Path
import json

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
PREFIX = "followup-lake-kidney"
OUT = REPO / "reports" / "followup-lake-kidney-analysis.md"


def read(name):
    return pd.read_csv(RESULTS / f"{PREFIX}-{name}.csv")


def num(value, digits=5):
    if pd.isna(value):
        return "NA"
    return f"{value:.{digits}g}"


def primary_table(frame, endpoint):
    part = frame[(frame["donor_set"] == "primary_adult") & frame["model"].str.startswith("primary")]
    lines = []
    for assay, group in part.groupby("assay_release"):
        lines.append(f"\n**{assay}: {endpoint} primary donor models**\n")
        if endpoint == "P1 transport":
            lines.append("| Cell type | Donors | Slope/year (95% CI) | p | BH q | Matched-null p | Null BH q |")
            lines.append("|---|---:|---:|---:|---:|---:|---:|")
            for row in group.itertuples():
                lines.append(f"| {row.cell_type} | {row.n_donors} | {num(row.slope_per_year)} [{num(row.ci95_low)}, {num(row.ci95_high)}] | {num(row.pvalue)} | {num(row.q_bh)} | {num(row.random_pvalue)} | {num(row.random_q_bh)} |")
        else:
            lines.append("| Cell type | Donors | Slope/year (95% CI) | p | BH q |")
            lines.append("|---|---:|---:|---:|---:|")
            for row in group.itertuples():
                lines.append(f"| {row.cell_type} | {row.n_donors} | {num(row.slope_per_year)} [{num(row.ci95_low)}, {num(row.ci95_high)}] | {num(row.pvalue)} | {num(row.q_bh)} |")
    if not lines:
        lines.append("\nNo eligible primary models.\n")
    return "\n".join(lines)


def sensitivity_summary(p1, composition, hellinger):
    frames = [pd.concat([p1, composition, hellinger], ignore_index=True, sort=False)]
    lines = ["| Endpoint | Assay | Donor set | Model | Tests | Minimum BH q | BH q<0.05 |",
             "|---|---|---|---|---:|---:|---:|"]
    for keys, part in frames[0].groupby(["endpoint", "assay_release", "donor_set", "model"], dropna=False):
        lines.append(f"| {keys[0]} | {keys[1]} | {keys[2]} | {keys[3]} | {len(part)} | {num(part['q_bh'].min())} | {int((part['q_bh'] < 0.05).sum())} |")
    null = p1[(p1["donor_set"] == "primary_adult") & (p1["model"] == "primary_age_sex")]
    for assay, part in null.groupby("assay_release"):
        lines.append(f"| P1_transport matched null | {assay} | primary_adult | expression-decile-matched random signatures | {len(part)} | {num(part['random_q_bh'].min())} | {int((part['random_q_bh'] < 0.05).sum())} |")
    return "\n".join(lines)


def main():
    p1 = read("p1-models")
    composition = read("composition-models")
    hellinger_path = RESULTS / f"{PREFIX}-hellinger-models.csv"
    hellinger = pd.read_csv(hellinger_path) if hellinger_path.exists() else p1.iloc[0:0].copy()
    eligibility = read("celltype-eligibility")
    h_eligibility = read("hellinger-eligibility")
    crosswalk = read("donor-crosswalk")
    comparison = read("existing-comparison")

    verification = {}
    for assay in ["scCv3", "snCv3"]:
        verification[assay] = json.loads((RESULTS / f"{PREFIX}-{assay}-verification.json").read_text())

    verified_lines = []
    for assay, v in verification.items():
        donors = crosswalk[crosswalk["assay_release"] == assay]
        verified_lines.append(
            f"- **{assay}**: {v['n_obs']:,} rows, raw.X {v['raw_X_shape'][0]:,} x {v['raw_X_shape'][1]:,}; "
            f"{v['raw_X_nnz']:,} stored values all finite, nonnegative, and integer-valued. "
            f"{v['healthy_donors']} healthy donors ({v['exact_age_healthy_donors']} exact-age; "
            f"{v['decade_only_healthy_donors']} decade-only), {v['primary_adult_donors']} primary adults. "
            f"GEO reconciliation: {len(donors)} donors, {donors['geo_specimens'].sum()} specimens, "
            f"{donors['geo_libraries'].sum()} libraries. Assay={v['assay_labels']}; suspension={v['suspension_labels']}; "
            f"release-wide is_primary_data={v['is_primary_data_labels']}."
        )

    h_primary = hellinger[(hellinger["donor_set"] == "primary_adult")
                          & (hellinger["model"] == "primary_depth_adjusted")]
    h_table = primary_table(h_primary, "Hellinger") if not h_primary.empty else "\nNo eligible Hellinger models.\n"

    overlap = comparison[comparison["relationship"] == "overlapping_prior_release_not_independent"]
    compare_lines = ["| Endpoint | Assay/prior cohort | Cell type | Shared donor IDs | New slope | Prior slope | New q | Prior q |",
                     "|---|---|---|---:|---:|---:|---:|---:|"]
    for row in overlap.itertuples():
        compare_lines.append(
            f"| {row.endpoint} | {row.assay_release}/{row.existing_cohort} | {row.cell_type} | "
            f"{row.n_shared_donor_ids} | {num(row.new_slope_per_year)} | {num(row.old_slope_per_year)} | "
            f"{num(row.new_q_bh)} | {num(row.old_q_bh)} |"
        )
    if overlap.empty:
        compare_lines.append("| No exact-label compatible overlapping-release estimates | NA | NA | NA | NA | NA | NA | NA |")

    ineligible = []
    for assay, part in eligibility[~eligibility["eligible"]].groupby("assay_release"):
        ineligible.append(f"- {assay}: {len(part)} released cell types failed the six-donor/20-cell gate: " + ", ".join(part["cell_type"]) + ".")
    sn_h = h_eligibility[(h_eligibility["assay_release"] == "snCv3") & (~h_eligibility["hellinger_eligible"])]
    sn_reasons = sorted(set(sn_h["reason"]))

    text = f"""# Restricted public Lake kidney analysis

## Verdict

The restricted public, assay-specific analysis is complete. It uses donors—not cells or libraries—as inferential units and does not pool scCv3 with snCv3. This is not the blocked full-45-reference-donor analysis and not controlled raw-sequencing work.

No result below is called replication merely because it is non-significant. The prior `Kidney_SCR_v15` and `Kidney_SN_v15` analyses use overlapping Lake releases and are release/method comparisons, not independent replication. The complete estimates, sensitivity models, per-30-year effects, deterministic `slope*30/122.5` lifespan-scaled effects, confidence intervals, and BH families are in the aggregate CSVs listed below.

## Frozen design and resources

The design was frozen before expression download or outcomes in `data/manifest/lake_kidney_analysis_plan.json` at HEAD `a335764`. Primary models include healthy exact-age adults aged 20–120 years. snCv3 donor 3490 (14 years) is excluded from every primary model and included only in the declared all-exact-age P1/CLR sensitivity. Decade-only donors 18-142, 18-162, and 18-312 are never assigned numeric ages.

Eligible cell types require at least six donors, 20 cells per donor, and identifiable covariates. P1 uses established `P1_Clean`, per-cell log1p(CPM), donor means, and age+sex OLS with a 200-draw expression-decile-matched random-signature null (seed 20261105). CLR uses all released non-unknown cell types per assay and a 0.5-cell pseudocount. Hellinger uses the established 20–40-year reference, cross-fitting by donor, 1,500-UMI cell filter, technical-gene exclusion, and reference/donor floors. No all-gene donor-pseudobulk model was predeclared or added.

Downloads and processing ran sequentially through `scripts/resource_guard.py` with a 30-GB disk reserve, 3-GB download cap, and 6-GB memory per processing job. scCv3 finished before snCv3 started. No commit or push was made.

## Verified inputs and exclusions

{chr(10).join(verified_lines)}

Sources are the public [CELLxGENE collection](https://cellxgene.cziscience.com/collections/bcb61471-2a44-4d00-a0af-ff085512674c) per-assay H5ADs and the public GEO metadata recorded in `data/manifest/lake_kidney.json`. Each object's own disease, development stage, donor ID, sex, assay, suspension, and released `cell_type` annotations were used. Feature IDs are unique; duplicate feature symbols are retained as distinct released features, and P1 symbol matching includes every matching feature as in the established scorer. GEO supplies specimen/library reconciliation but not ages.

The release-wide `is_primary_data=False` flag was recorded rather than silently discarding every cell. Each per-assay object was analyzed independently, with no cross-assay concatenation. Disease cells, decade-only ages, unknown annotations, and cell types failing the frozen donor/cell floors were excluded from the relevant models.

## P1 transport

{primary_table(p1, "P1 transport")}

The conventional donor-model p-values and the matched-random-null empirical p-values are separate inferential summaries with separate assay-wise BH families. A small matched-null q does not make a non-significant conventional donor-model result a confirmed replication.

## Donor-level CLR composition

{primary_table(composition, "CLR composition")}

These are associations among recovered annotated cell fractions, not absolute tissue cell counts. Zero cell-type counts remain in each donor's composition through the declared 0.5 pseudocount. Total-cell-adjusted and all-exact-age sensitivities are in `results/followup-lake-kidney-composition-models.csv`.

## Hellinger displacement

{h_table}

snCv3 Hellinger is **untestable** under the frozen reference rule: {', '.join(sn_reasons)}. Only KRP446 (29 years) is an exact-age primary adult in the 20–40-year reference window, so the required two young reference donors are unavailable. This is an eligibility outcome, not a null result. scCv3 reference-floor details are in `results/followup-lake-kidney-hellinger-eligibility.csv`; donor-level distances are in `results/followup-lake-kidney-hellinger-donors.csv`.

## Sensitivity and multiplicity families

BH was applied separately by assay, endpoint, donor set, and labeled model. The table summarizes every fitted family; it does not pool assays or endpoints.

{sensitivity_summary(p1, composition, hellinger)}

The lifespan quantity is only the predeclared deterministic rescaling of each slope by 30/122.5, not an additional test. The all-exact-age snCv3 sensitivity is the only lifespan/age-range sensitivity beyond primary adult models and depth/total-cell adjustment.

## Comparison with existing kidney results

Only exact cell-type labels and compatible donor-level estimands were compared. For snCv3 P1, the all-exact-age model is used because the prior 10-donor `Kidney_SN_v15` model included the 14-year-old; the primary adult-only snCv3 model is not presented as the same donor set. CLR was not compared across releases because the released composition universes and QC differ.

{chr(10).join(compare_lines)}

Full compatible comparisons with the separate `Kidney_NormalAdult` and `Kidney_LivingDonor` studies are in `results/followup-lake-kidney-existing-comparison.csv`, labeled `separate_study`. They are not pooled with Lake and do not turn overlapping releases into independent evidence.

## Untestable, not attempted, and interpretation limits

{chr(10).join(ineligible)}
- snCv3 Hellinger: untestable because the young-reference floor fails, as detailed above.
- Full 45-reference-donor analysis: still blocked by incomplete public exact-age/common-mapping coverage; no ages were imputed.
- All-gene donor-pseudobulk models: not predeclared and not performed; there is therefore no gene-level multiplicity family or gene-level claim.
- Controlled KPMP/NIDDK/dbGaP raw sequencing: not requested, contacted, downloaded, or analyzed.
- SNARE2, Slide-seq2, Visium, muscle, and other follow-up cohorts: outside this restricted task and not touched.
- Non-significant P1, CLR, or Hellinger results are not equivalence, absence-of-effect, or replication claims. Small donor numbers and release differences limit inference.

## Reproducibility and delivered aggregates

- `scripts/run_followup_lake_kidney.py`: backed raw-count validation, donor aggregation, P1/null, CLR, Hellinger, comparison, and manifest update.
- `scripts/report_followup_lake_kidney.py`: renders this report from aggregate CSVs.
- `results/followup-lake-kidney-*-verification.json`: per-assay raw-count and metadata gates.
- `results/followup-lake-kidney-donor-crosswalk.csv`: exact ages, health, primary/sensitivity sets, GEO patient/specimen/library mapping.
- `results/followup-lake-kidney-celltype-eligibility.csv` and `-hellinger-eligibility.csv`: inclusion gates and untestable reasons.
- `results/followup-lake-kidney-p1-donor-scores.csv`, `-p1-models.csv`: donor scores, null summaries, effects/CIs, and both BH families.
- `results/followup-lake-kidney-composition-donor-counts.csv`, `-composition-models.csv`: zero-inclusive donor counts and CLR models.
- `results/followup-lake-kidney-hellinger-donors.csv`, `-hellinger-models.csv`: eligible donor distances and models.
- `results/followup-lake-kidney-existing-comparison.csv`: exact-label compatible prior/separate-study comparisons and shared donor counts.
- `tests/test_followup_lake_kidney.py`: helper and delivered-output validation.

Run with the Agescape environment:

```bash
python scripts/resource_guard.py run --path /home/reza/projects/agescape-revision-data/lake_kidney \
  --min-free-gb 30 --memory-gb 6 -- \
  python scripts/run_followup_lake_kidney.py run --assay scCv3
# Repeat for snCv3 only after scCv3 completes, then:
python scripts/resource_guard.py run --path /home/reza/projects/agescape-revision-data/lake_kidney \
  --min-free-gb 30 --memory-gb 6 -- \
  python scripts/run_followup_lake_kidney.py finalize
python scripts/report_followup_lake_kidney.py
```

## Validation record

- `python -m pytest tests/test_followup_lake_kidney.py -q`: **10 passed**.
- `python -m pytest -q`: **116 passed, 1 skipped** after the TS2 manifest portability fix in commit `1c527c0`. The Lake manifest passes the same path/security checks.
- `git diff --check`: passed. HEAD at the frozen design was `a335764`; no commit or push was performed and parent review is required.
"""
    OUT.write_text(text)
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
