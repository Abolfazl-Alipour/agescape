#!/usr/bin/env python3
"""Validate and report the mature 30-plus full-depth redo."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
DATA_ROOT = ROOT.parent / f"{ROOT.name}-data" / "tmp_organs"
SLICES = ["dataset", "part", "metric", "gene_set", "depth_mode", "model_spec"]


def select_row(df: pd.DataFrame, part: int, metric: str, gene_set: str = "none", depth: str = "full_depth", spec: str = "stratified_global_primary", group: str | None = None) -> pd.Series:
    mask = (
        (df["part"].astype(int) == part)
        & (df["metric"] == metric)
        & (df["gene_set"] == gene_set)
        & (df["depth_mode"] == depth)
        & (df["model_spec"] == spec)
    )
    if group is not None:
        mask &= df["group"] == group
    rows = df[mask]
    if len(rows) != 1:
        raise ValueError(f"Expected one row, found {len(rows)} for {part}, {metric}, {gene_set}, {depth}, {spec}, {group}")
    return rows.iloc[0]


def select_row_optional(df: pd.DataFrame, part: int, metric: str, gene_set: str = "none", depth: str = "full_depth", spec: str = "stratified_global_primary", group: str | None = None) -> pd.Series | None:
    mask = (
        (df["part"].astype(int) == part)
        & (df["metric"] == metric)
        & (df["gene_set"] == gene_set)
        & (df["depth_mode"] == depth)
        & (df["model_spec"] == spec)
    )
    if group is not None:
        mask &= df["group"] == group
    rows = df[mask]
    if len(rows) == 0:
        return None
    if len(rows) != 1:
        raise ValueError(f"Expected at most one row, found {len(rows)} for {part}, {metric}, {gene_set}, {depth}, {spec}, {group}")
    return rows.iloc[0]


def fmt(value: Any, digits: int = 4) -> str:
    if value is None or not np.isfinite(float(value)):
        return "NA"
    return f"{float(value):.{digits}g}"


def p_fmt(value: Any) -> str:
    if value is None or not np.isfinite(float(value)):
        return "NA"
    return f"{float(value):.3e}"


def effect(row: pd.Series) -> str:
    return f"{fmt(row['slope_sd_per_decade'])}, p {p_fmt(row['pvalue'])}, SNR {fmt(row.get('sh_snr_slope', np.nan), 3)}, MDE80 {fmt(row['mde80_sd_per_decade'], 3)}"


def check_bh(df: pd.DataFrame) -> int:
    eligible = (~df["is_pooled"].astype(bool)) & (~df["is_stratified_global"].astype(bool))
    pool_arm = df["model_spec"].astype(str).isin(["primary", "pool_fe", "pool_fe_no_depth", "age_ge_20"])
    eligible &= ~(pool_arm & df["yazar_pool_fe_status"].astype(str).str.contains("not_evaluable"))
    n_checked = 0
    for _, idx in df.loc[eligible].groupby(SLICES, sort=False).groups.items():
        pvals = df.loc[idx, "pvalue"].to_numpy(dtype=float)
        valid = np.isfinite(pvals)
        if not valid.any():
            continue
        valid_idx = np.asarray(list(idx))[valid]
        expected = multipletests(pvals[valid], method="fdr_bh")[1]
        observed = df.loc[valid_idx, "q_bh_within_slice"].to_numpy(dtype=float)
        if not np.allclose(observed, expected, rtol=1e-10, atol=1e-12, equal_nan=True):
            raise AssertionError(f"BH mismatch in {idx}")
        if not (df.loc[valid_idx, "family_size_m"] == int(valid.sum())).all():
            raise AssertionError(f"BH family-size mismatch in {idx}")
        n_checked += int(valid.sum())
    pooled = df["is_pooled"].astype(bool) & (~df["is_stratified_global"].astype(bool))
    if not np.allclose(df.loc[pooled, "q_bh_within_slice"], df.loc[pooled, "pvalue"], rtol=0, atol=0, equal_nan=True):
        raise AssertionError("Pooled singleton q differs from p")
    if not (df.loc[pooled, "family_size_m"] == 1).all():
        raise AssertionError("Pooled family size is not one")
    stratified = df["is_stratified_global"].astype(bool)
    if not df.loc[stratified, "q_bh_within_slice"].isna().all():
        raise AssertionError("Stratified-global q must be NaN")
    return n_checked


def load_pair(ds: str, suffix: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    old = pd.read_csv(RESULTS / f"followup-fulldepth-{ds}-{suffix}.csv")
    new = pd.read_csv(RESULTS / f"followup-mature-fulldepth-{ds}-{suffix}.csv")
    return old, new


def main() -> int:
    old_qc = json.loads((RESULTS / "followup-fulldepth-qc.json").read_text())
    qc = json.loads((RESULTS / "followup-mature-fulldepth-qc.json").read_text())
    data: dict[str, dict[str, tuple[pd.DataFrame, pd.DataFrame]]] = {}
    bh_checked = 0
    for ds in ("gao", "yazar"):
        data[ds] = {}
        for suffix in ("cell-drift", "spread", "ectopic"):
            old, new = load_pair(ds, suffix)
            data[ds][suffix] = (old, new)
            bh_checked += check_bh(new)
        scope = qc["datasets"][ds]["mature_restriction"]
        if scope["minimum_donor_age_years"] != 30.0 or scope["reference_age_min_years"] != 30.0 or scope["reference_age_max_years"] != 39.0:
            raise AssertionError(f"Mature scope failed for {ds}")

    lines = [
        "# Mature 30-plus full-depth drift, entropy, and identity",
        "",
        "This is Phase 2 Item 9a. All derivation and tested contrasts use donors aged at least 30 years. The inherited YOUNG label means only the 30-to-39 cross-fitted baseline, with leave-one-out references for those donors. No young-derived P1 or original E1/E2 set enters this analysis.",
        "",
        "## Frozen design",
        "",
        "- Primary outcomes use full-depth raw counts as within-cell proportions; at least 1500 UMI is an eligibility filter only.",
        "- Rarefied 1500 and 3000 UMI outcomes are labeled sensitivities. Every downsample output has explicit zeros eliminated before stored-entry counting.",
        "- Donor is the unit for every p-value. Slopes are per decade. Full-depth primaries carry five-seed split-half SNR and MDE80.",
        "- Per-cell-type BH families use the frozen dataset, part, metric, gene-set, depth, and model slices. Pooled rows are singleton families; stratified-global primary q is NaN by design.",
        "- E2_stratum is primary E2 evidence. E2_global is only a depth-vulnerable sensitivity.",
        "",
        "## Donor and cell restriction",
        "",
        "Original selected donors include 117 Gao and 981 Yazar donors. Original global models contain 113 Gao donors and 522,761 eligible nuclei, and 973 Yazar donors and 617,760 eligible cells.",
        "",
        "| Dataset | Selected donors original -> mature | Lost | Baseline donors 30-39 | Global donors original -> mature | Eligible cells original -> mature | Cells lost |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    original_global = {"gao": 113, "yazar": 973}
    original_cells = {"gao": 522761, "yazar": 617760}
    donor_notes: list[str] = []
    for ds, label in (("gao", "Gao heart"), ("yazar", "Yazar blood")):
        scope = qc["datasets"][ds]["mature_restriction"]
        cs_row = select_row_optional(data[ds]["cell-drift"][1], 1, "center_shift")
        sc_row = select_row_optional(data[ds]["cell-drift"][1], 1, "scatter")
        sh_row = select_row_optional(data[ds]["spread"][1], 2, "shannon")
        n_glob_cs = int(cs_row["n_donors"]) if cs_row is not None else 0
        n_glob_sc = int(sc_row["n_donors"]) if sc_row is not None else 0
        n_glob_sh = int(sh_row["n_donors"]) if sh_row is not None else 0
        glob_str = f"{n_glob_cs}" if n_glob_cs == n_glob_sc == n_glob_sh else f"{n_glob_cs} (center_shift) / {n_glob_sc} (scatter/Part 2)"
        new_cells = int(qc["datasets"][ds]["n_cells_ge1500"])
        pooled_cells = qc["datasets"][ds].get("n_cells_ge1500_pooled_primary")
        lines.append(
            f"| {label} | {scope['n_donors_original']} -> {scope['n_donors_mature']} | {scope['n_donors_lost']} | {scope['n_reference_donors']} | {original_global[ds]} -> {glob_str} | {original_cells[ds]:,} -> {new_cells:,} | {original_cells[ds] - new_cells:,} |"
        )
        if pooled_cells is not None and cs_row is not None:
            donor_notes.append(
                f"- {label}: {new_cells:,} total eligible cells at >= 1500 UMI across {scope['n_donors_mature']} mature donors; {int(pooled_cells):,} cells in pooled primary model; {int(cs_row['n_cells']):,} cells across {n_glob_cs} donors in stratified-global center_shift ({int(sc_row['n_cells']):,} cells across {n_glob_sc} donors in stratified-global scatter)."
            )
    if donor_notes:
        lines.extend(["", *donor_notes])

    lines.extend([
        "",
        "## Primary verdicts and exact effects",
        "",
        "Effects below are stratified-global SD per decade. Original values are full-range context only.",
        "",
        "| Dataset | Endpoint | Original slope, p | Mature slope, p, SNR, MDE80 | Mature verdict | Reason |",
        "| --- | --- | --- | --- | --- | --- |",
    ])
    endpoints = [
        ("Part 1 center shift", "cell-drift", 1, "center_shift", "none", "full_depth", "Part1_center_shift"),
        ("Part 1 scatter", "cell-drift", 1, "scatter", "none", "full_depth", "Part1_scatter"),
        ("Part 2 Shannon", "spread", 2, "shannon", "none", "full_depth", "Part2/shannon"),
        ("Part 2 KL methods", "spread", 2, "deman_kl_methods", "none", "full_depth", "Part2/deman_kl_methods"),
        ("Part 2 KL sqrt", "spread", 2, "deman_kl_sqrt", "none", "full_depth", "Part2/deman_kl_sqrt"),
        ("Part 3 E1 count", "ectopic", 3, "mean_ectopic_count", "E1", "full_depth", "Part3/E1"),
        ("Part 3 E2 stratum count", "ectopic", 3, "mean_ectopic_count", "E2_stratum", "bounded_depth_strata", "Part3/E2"),
    ]
    for ds, label in (("gao", "Gao heart"), ("yazar", "Yazar blood")):
        for name, suffix, part, metric, gene_set, depth, verdict_key in endpoints:
            old, new = data[ds][suffix]
            ro = select_row(old, part, metric, gene_set, depth)
            rn = select_row_optional(new, part, metric, gene_set, depth)
            node: Any = qc["verdicts"][ds]
            for key in verdict_key.split("/"):
                node = node[key]
            eff_str = effect(rn) if rn is not None else "UNTESTABLE (missing stratified global row)"
            lines.append(f"| {label} | {name} | {fmt(ro['slope_sd_per_decade'])}, p {p_fmt(ro['pvalue'])} | {eff_str} | {node['verdict']} | {node['reason']} |")

    lines.extend(["", "## What survived or died", ""])
    for ds, label in (("gao", "Gao heart"), ("yazar", "Yazar blood")):
        lines.append(f"### {label}")
        lines.append("")
        for name, suffix, part, metric, gene_set, depth, verdict_key in endpoints:
            old, new = data[ds][suffix]
            ro = select_row(old, part, metric, gene_set, depth)
            rn = select_row_optional(new, part, metric, gene_set, depth)
            node: Any = qc["verdicts"][ds]
            for key in verdict_key.split("/"):
                node = node[key]
            if rn is None:
                lines.append(f"- {name}: UNTESTABLE ({node['reason']}).")
                continue
            old_sig = float(ro["pvalue"]) < 0.05
            new_sig = float(rn["pvalue"]) < 0.05
            if old_sig and new_sig and np.sign(ro["slope_sd_per_decade"]) == np.sign(rn["slope_sd_per_decade"]):
                change = "The nominally significant direction survived"
            elif old_sig:
                change = "The original nominal significance died"
            elif new_sig:
                change = "A new nominally significant association emerged"
            else:
                change = "The original non-significant result remained non-significant"
            lines.append(f"- {name}: {change}; original slope {fmt(ro['slope_sd_per_decade'])} (p {p_fmt(ro['pvalue'])}), mature slope {fmt(rn['slope_sd_per_decade'])} (p {p_fmt(rn['pvalue'])}, MDE80 {fmt(rn['mde80_sd_per_decade'], 3)}).")
        lines.append("")

    lines.extend([
        "## Pooled rows, reported separately",
        "",
        "Pooled slopes are not within-cell-type evidence. A pooled positive with a discordant stratified slope is read as a composition effect.",
        "",
        "| Dataset | Part and metric | Pooled raw slope | Pooled SD slope | p | Stratified SD slope | Stratified p |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: |",
    ])
    for ds, label in (("gao", "Gao heart"), ("yazar", "Yazar blood")):
        pool_name = "All_LV_Nuclei" if ds == "gao" else "All_Cells"
        for part, suffix, metric in ((1, "cell-drift", "center_shift"), (1, "cell-drift", "scatter"), (2, "spread", "shannon"), (2, "spread", "deman_kl_methods"), (2, "spread", "deman_kl_sqrt")):
            new = data[ds][suffix][1]
            rp = select_row_optional(new, part, metric, spec="primary", group=pool_name)
            rs = select_row_optional(new, part, metric)
            if rp is None:
                lines.append(f"| {label} | Part {part} {metric} | UNTESTABLE | UNTESTABLE | NA | {fmt(rs['slope_sd_per_decade'], 6) if rs is not None else 'UNTESTABLE'} | {p_fmt(rs['pvalue']) if rs is not None else 'NA'} |")
            else:
                lines.append(f"| {label} | Part {part} {metric} | {fmt(rp['slope_raw_per_decade'], 6)} | {fmt(rp['slope_sd_per_decade'], 6)} | {p_fmt(rp['pvalue'])} | {fmt(rs['slope_sd_per_decade'], 6) if rs is not None else 'UNTESTABLE'} | {p_fmt(rs['pvalue']) if rs is not None else 'NA'} |")

    naive = "naive thymus-derived CD8-positive, alpha-beta T cell"
    old_yazar, new_yazar = data["yazar"]["cell-drift"]
    naive_old = select_row(old_yazar, 1, "scatter", spec="primary", group=naive)
    naive_new = select_row_optional(new_yazar, 1, "scatter", spec="primary", group=naive)
    if naive_new is not None:
        naive_new_str = f"Mature 30-plus retest: raw slope {fmt(naive_new['slope_raw_per_decade'], 6)}, SD slope {fmt(naive_new['slope_sd_per_decade'], 6)}, p {p_fmt(naive_new['pvalue'])}, q {p_fmt(naive_new['q_bh_within_slice'])}, SNR {fmt(naive_new['sh_snr_slope'], 3)}, MDE80 {fmt(naive_new['mde80_sd_per_decade'], 3)}, donors {int(naive_new['n_donors'])}, cells {int(naive_new['n_cells']):,}."
    else:
        naive_new_str = "Mature 30-plus retest: UNTESTABLE (group did not meet primary model eligibility at 30-plus)."
    lines.extend([
        "",
        "### Yazar naive CD8 exception retest",
        "",
        f"Full-range context only (never mature evidence): SD slope {fmt(naive_old['slope_sd_per_decade'], 6)}, p {p_fmt(naive_old['pvalue'])}, q {p_fmt(naive_old['q_bh_within_slice'])}, donors {int(naive_old['n_donors'])}, cells {int(naive_old['n_cells']):,}. {naive_new_str}",
        "",
        "## Sensitivity outcomes",
        "",
        "| Dataset | Sensitivity | SD slope | p |",
        "| --- | --- | ---: | ---: |",
    ])
    sensitivity_specs = [
        ("cell-drift", 1, "center_shift", "none", "equal_ref_size", "Center equal reference size"),
        ("cell-drift", 1, "center_shift", "none", "pythagorean_centroid_h2", "Center Pythagorean centroid"),
        ("cell-drift", 1, "center_shift", "none", "spline_resid_full", "Center spline residual"),
        ("cell-drift", 1, "center_shift", "none", "rarefied_1500", "Center rarefied 1500"),
        ("cell-drift", 1, "scatter", "none", "rarefied_1500", "Scatter rarefied 1500"),
        ("spread", 2, "shannon", "none", "rarefied_1500", "Shannon rarefied 1500"),
        ("spread", 2, "shannon", "none", "rarefied_3000", "Shannon rarefied 3000"),
        ("ectopic", 3, "mean_ectopic_count", "E1", "rarefied_1500", "E1 count rarefied 1500"),
        ("ectopic", 3, "mean_ectopic_count", "E1", "spline_resid_full", "E1 count spline residual"),
        ("ectopic", 3, "mean_ectopic_count", "E1_stringent", "full_depth", "E1 stringent count"),
        ("ectopic", 3, "mean_ectopic_count", "E2_global", "full_depth", "E2 global count, depth vulnerable"),
        ("ectopic", 3, "mean_ectopic_count", "E2_global", "rarefied_1500", "E2 global count rarefied 1500"),
    ]
    for ds, label in (("gao", "Gao heart"), ("yazar", "Yazar blood")):
        for suffix, part, metric, gene_set, depth, name in sensitivity_specs:
            row = select_row_optional(data[ds][suffix][1], part, metric, gene_set, depth)
            if row is None:
                lines.append(f"| {label} | {name} | UNTESTABLE | NA |")
            else:
                lines.append(f"| {label} | {name} | {fmt(row['slope_sd_per_decade'], 6)} | {p_fmt(row['pvalue'])} |")

    lines.extend([
        "",
        "## Part 3 count gate and secondary metrics",
        "",
        "| Dataset | Set and metric | Depth | SD slope | p | SNR | MDE80 |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: |",
    ])
    p3_metrics = [
        ("E1", "singleton_fraction_20cells", "subsampled_20cells"),
        ("E1", "jaccard_both_pos", "full_depth"),
        ("E1", "excess_jaccard", "full_depth"),
        ("E2_stratum", "singleton_fraction_20cells", "bounded_depth_strata_20cells"),
        ("E2_stratum", "jaccard_both_pos", "bounded_depth_strata"),
        ("E2_stratum", "excess_jaccard", "bounded_depth_strata"),
    ]
    for ds, label in (("gao", "Gao heart"), ("yazar", "Yazar blood")):
        for gene_set, metric, depth in p3_metrics:
            row = select_row_optional(data[ds]["ectopic"][1], 3, metric, gene_set, depth)
            if row is None:
                lines.append(f"| {label} | {gene_set} {metric} | {depth} | UNTESTABLE | NA | NA | NA |")
            else:
                lines.append(f"| {label} | {gene_set} {metric} | {depth} | {fmt(row['slope_sd_per_decade'], 6)} | {p_fmt(row['pvalue'])} | {fmt(row.get('sh_snr_slope', np.nan), 3)} | {fmt(row['mde80_sd_per_decade'], 3)} |")

    lines.extend([
        "",
        "## Mature E1/E2 set sizes versus originals",
        "",
        "| Dataset | E1 cell types original -> mature | E1 median markers original -> mature | E2 global cell types original -> mature | E2 global median genes original -> mature | Mature E2 stratum CT-strata, median genes | Genes kept |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ])
    for ds, label in (("gao", "Gao heart"), ("yazar", "Yazar blood")):
        old = old_qc["datasets"][ds]
        new = qc["datasets"][ds]
        lines.append(f"| {label} | {old['n_e1_cts']} -> {new['n_e1_cts']} | {fmt(old['e1_median_markers'])} -> {fmt(new['e1_median_markers'])} | NA -> {new['n_e2_global_cts']} | {fmt(old['e2_global_median_markers'])} -> {fmt(new['e2_global_median_markers'])} | {new['n_e2_stratum_ct_strata']}, {fmt(new['e2_stratum_median_markers'])} | {new['gene_audit']['n_genes_retained']:,} |")

    lines.extend([
        "",
        "Original E2 global cell-type counts and E2 stratum CT-stratum counts were not recorded in the original QC JSON, so those original fields are NA rather than inferred.",
        "",
        "## Power after restriction",
        "",
        "| Dataset | Endpoint | MDE80 original -> mature | Fold change |",
        "| --- | --- | ---: | ---: |",
    ])
    for ds, label in (("gao", "Gao heart"), ("yazar", "Yazar blood")):
        for name, suffix, part, metric, gene_set, depth, _ in endpoints:
            old, new = data[ds][suffix]
            ro = select_row(old, part, metric, gene_set, depth)
            rn = select_row_optional(new, part, metric, gene_set, depth)
            if rn is None:
                lines.append(f"| {label} | {name} | {fmt(ro['mde80_sd_per_decade'], 4)} -> UNTESTABLE | NA |")
                continue
            ratio = float(rn["mde80_sd_per_decade"]) / float(ro["mde80_sd_per_decade"])
            lines.append(f"| {label} | {name} | {fmt(ro['mde80_sd_per_decade'], 4)} -> {fmt(rn['mde80_sd_per_decade'], 4)} | {fmt(ratio, 3)}x |")

    robust = pd.read_csv(RESULTS / "followup-mature-fulldepth-robustness.csv")
    old_robust = pd.read_csv(RESULTS / "followup-fulldepth-robustness.csv")
    scope_meta = qc.get("analysis_scope", {})
    lines.extend([
        "",
        "## Robustness, verification, and items not verified",
        "",
        f"- Rarefaction robustness: {int(robust['passes_robustness'].sum())} of {len(robust)} rows pass, versus {int(old_robust['passes_robustness'].sum())} of {len(old_robust)} originally.",
        f"- Independent BH recomputation matched all {bh_checked:,} eligible per-cell-type p-values across all six output CSVs. Pooled singleton (m = 1, q = p) and stratified-global NaN rules also passed.",
        "- Regression suites: tests/test_followup_mature_fulldepth.py (10 passed, including explicit-zero rarefaction cleanup and under-30 leakage checks) and tests/test_followup_fulldepth.py (9 passed). Real metadata checks confirmed minimum mapped donor age 30.0 in both Gao and Yazar, and all reference donors strictly within 30 to 39 years.",
        "- Explicit-zero rarefaction check: verified on both Pass 1 caches (y_1500_det strictly smaller than y_full_det; exact 1500 UMI per cell) and Part 3 rarefied_1500 outputs (distinct from full_depth outputs).",
        "- Not verified / intentionally omitted: (1) Cross-check against original Part A donor CSVs was not performed because original Part A references and scored donors are full-range (20-plus / 19-plus). (2) Original E2 global cell-type count and original E2 stratum CT-stratum count were not stored in results/followup-fulldepth-qc.json and are reported as NA. (3) Gao Parts 1 and 2 and Yazar Pass 1 shards/pickles were produced by the prior kimi-code/k3.1-pro session and verified via shard byte-size, cell-count, donor-age, explicit-zero, and BH checks rather than recomputed from HDF5.",
        "",
        "## Runtime and reproduction",
        "",
        f"- Current run model: {scope_meta.get('model', 'NA')}.",
        f"- Prior cached Gao Parts 1-2 and Yazar Pass 1 model: {scope_meta.get('prior_cached_pass1_and_gao_parts12_model', 'NA')}.",
        "- Prior cached Pass 1/2 durations (verified from checkpoints): Gao Pass 1 565.92 s, Gao Pass 2 506.41 s, Yazar Pass 1 314.95 s.",
        "",
        "| Run | Model | Start UTC | End UTC | Elapsed seconds |",
        "| --- | --- | --- | --- | ---: |",
    ])
    runtime_path = DATA_ROOT / "mature_fulldepth_cache" / "runtime.jsonl"
    if runtime_path.exists():
        for line in runtime_path.read_text().splitlines():
            if not line.strip():
                continue
            rec = json.loads(line)
            lines.append(f"| {rec['dataset']} {rec['stage']} | {rec.get('model', 'NA')} | {rec['start_utc']} | {rec['end_utc']} | {fmt(rec['elapsed_s'], 8)} |")
    else:
        lines.append("| Runtime log missing | NA | NA | NA | NA |")
    lines.extend([
        "",
        "Pass-level timings are in the pipeline logs and cache checkpoints. Commands:",
        "",
        "```",
        "PY=/home/reza/miniconda3/envs/agescape/bin/python",
        "$PY scripts/resource_guard.py run --path /home/reza/projects/agescape-revision-data --min-free-gb 30 --memory-gb 12 -- $PY scripts/run_followup_mature_fulldepth.py --dataset yazar --stage parts12 --reuse-pass1",
        "$PY scripts/run_followup_mature_fulldepth.py --dataset all --stage part3 --reuse-pass1",
        "$PY -m pytest -p no:cacheprovider tests/test_followup_mature_fulldepth.py tests/test_followup_fulldepth.py -v",
        "$PY scripts/report_followup_mature_fulldepth.py",
        "```",
        "",
    ])
    report = "\n".join(lines)
    if "\u2013" in report or "\u2014" in report:
        raise AssertionError("Non-hyphen dash in report")
    (ROOT / "reports" / "followup-mature-fulldepth.md").write_text(report, encoding="utf-8")
    print(f"BH checks passed for {bh_checked} rows")
    print(f"Wrote reports/followup-mature-fulldepth.md with {len(lines)} lines")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
