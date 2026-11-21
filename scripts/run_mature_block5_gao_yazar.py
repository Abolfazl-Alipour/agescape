#!/usr/bin/env python3
"""Extend mature Block 5 to Gao and Yazar using the frozen original models."""

from __future__ import annotations

import gc
import os
import pickle
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_var] = "1"

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem
from scipy import stats
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_followup_11pt_block5_p0_transport import YAZAR_CACHE
from scripts.run_followup_hcata_pool_sensitivity import load_cache
from scripts.run_followup_organ_17_heart_gao import (
    H5AD_PATH,
    STUDY_DUMMY_MAP,
    build_donor_audit,
    build_donor_meta_table,
    stream_pseudobulk_counts,
)
from scripts.run_mature_p1m_transport import load_programs
from src.followup_organs import (
    finalize_pseudobulk_group,
    run_family_c_p1_signature,
    select_design_columns,
)

RESULTS = REPO / "results"
REPORTS = REPO / "reports"
SCRATCH = REPO / "scratch" / "mature-blocks-458"
GAO_CACHE = SCRATCH / "mature_block5_gao_30plus_inputs.pkl"
TRANSPORT_CSV = RESULTS / "followup-mature-block5-gao-yazar-transport.csv"
COMPARISON_CSV = RESULTS / "followup-mature-block5-p0m-comparison.csv"
COVERAGE_CSV = RESULTS / "followup-mature-block5-coverage.csv"
SUMMARY_CSV = RESULTS / "followup-mature-block5-summary.csv"
RUNTIME_CSV = RESULTS / "followup-mature-block5-runtime.csv"
REPORT = REPORTS / "followup-mature-block5-p0m-comparison.md"
CORE_TRANSPORT = RESULTS / "followup-mature-p1m-transport.csv"
CORE_NULL = RESULTS / "followup-mature-p1m-transport-null.csv"
CORE_COMPARISON = RESULTS / "followup-mature-p1m-block5-p0m-comparison.csv"
CORE_UNTESTABLE = RESULTS / "followup-mature-p1m-transport-untestable.csv"
ORIGINAL_TRANSPORT = RESULTS / "followup-11pt-block5-p0-transport.csv"
ORIGINAL_FAMILY = RESULTS / "human-transport-global-multiplicity.csv"
MODEL = "claude-opus-5-5-high"
TIERS = ("P0m", "P1m")
MIN_AGE = 30.0
ITERATIONS = 2000
SEED = 20261105
ALPHA = 0.05


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def wait_for_memory(stage: str) -> None:
    while True:
        values = {}
        with open("/proc/meminfo", encoding="ascii") as handle:
            for line in handle:
                key, value = line.split(":", 1)
                values[key] = value.strip()
        available_kib = int(values["MemAvailable"].split()[0])
        available_gib = available_kib / 1024**2
        if available_gib >= 20.0:
            print(f"{stage}: {available_gib:.2f} GiB available", flush=True)
            return
        print(f"{stage}: only {available_gib:.2f} GiB available; retrying in 60 seconds", flush=True)
        time.sleep(60)


def original_rosters() -> tuple[pd.DataFrame, dict[tuple[str, str], list[str]]]:
    original = pd.read_csv(ORIGINAL_TRANSPORT)
    followup_original = original[original["cohort"].isin(["Gao_Heart", "Yazar_Blood"])]
    rosters = {}
    for (cohort, arm), sub in followup_original.groupby(["cohort", "arm"], sort=True):
        units = sorted(sub["unit"].unique())
        rosters[(cohort, arm)] = units
    assert {key: len(value) for key, value in rosters.items()} == {
        ("Gao_Heart", "primary"): 14,
        ("Gao_Heart", "unadjusted"): 14,
        ("Yazar_Blood", "primary"): 17,
        ("Yazar_Blood", "unadjusted"): 20,
    }
    return original, rosters


def load_gao_inputs(roster: list[str]):
    if GAO_CACHE.exists():
        print(f"Gao: loading mature cache {GAO_CACHE}", flush=True)
        with GAO_CACHE.open("rb") as handle:
            cached = pickle.load(handle)
        assert cached["min_age"] == MIN_AGE
        assert cached["roster"] == roster
        assert cached["donor_meta"]["age_years"].min() >= MIN_AGE
        for meta in cached["metas"].values():
            assert len(meta) == 0 or meta["age_years"].min() >= MIN_AGE
        return cached
    wait_for_memory("Gao stream")
    print("Gao: reading metadata and streaming 30-plus donors", flush=True)
    with h5py.File(H5AD_PATH, "r") as handle:
        obs = read_elem(handle["obs"])
        var = read_elem(handle["var"])
    audit = build_donor_audit(obs)
    donor_meta_all = build_donor_meta_table(audit)
    donor_meta = donor_meta_all[donor_meta_all["age_years"] >= MIN_AGE].copy()
    included = set(donor_meta.index)
    keep = obs["donor_id"].astype(str).isin(included).to_numpy()
    raw, counts, gate = stream_pseudobulk_counts(
        H5AD_PATH,
        keep,
        obs["donor_id"].astype(str).to_numpy(),
        obs["final_cell_type"].astype(str).to_numpy(),
        var.index.astype(str).to_numpy(),
    )
    expr = {}
    metas = {}
    for group in roster:
        if group not in raw:
            continue
        raw_group = raw[group]
        donors = list(raw_group.columns)
        cell_count = counts.loc[donors, "Total_Cells"] if group == "All_LV_Nuclei" else counts.loc[donors, group]
        expr[group], metas[group], _ = finalize_pseudobulk_group(
            raw_group, donor_meta.loc[donors], cell_count
        )
    cached = {
        "min_age": MIN_AGE,
        "roster": roster,
        "expr": expr,
        "metas": metas,
        "donor_meta": donor_meta,
        "n_original_donors": len(donor_meta_all),
        "gate": gate,
    }
    SCRATCH.mkdir(parents=True, exist_ok=True)
    temporary = GAO_CACHE.with_suffix(".tmp")
    with temporary.open("wb") as handle:
        pickle.dump(cached, handle, protocol=pickle.HIGHEST_PROTOCOL)
    temporary.replace(GAO_CACHE)
    del raw, expr, metas, obs, var
    gc.collect()
    return cached


def pool_inputs_30plus(roster: list[str]):
    donor_meta_all, cell_counts, expr_all, metas_all = load_cache(
        YAZAR_CACHE, drop_count_cols=("peripheral blood mononuclear cell",)
    )
    all_pools = sorted(donor_meta_all["pool_number"].dropna().astype(int).unique())
    pool_cols = [f"pool_{pool}" for pool in all_pools]
    mature_ids = donor_meta_all.index[donor_meta_all["age_years"] >= MIN_AGE].astype(str).tolist()
    expr = {}
    metas = {}
    for group in roster:
        if group not in expr_all:
            continue
        donors = [donor for donor in expr_all[group].columns.astype(str) if donor in set(mature_ids)]
        expr[group] = expr_all[group][donors].copy()
        metas[group] = metas_all[group].loc[donors].copy()
        assert len(donors) == 0 or metas[group]["age_years"].min() >= MIN_AGE
    donor_meta = donor_meta_all.loc[mature_ids].copy()
    del expr_all, metas_all
    gc.collect()
    return donor_meta, cell_counts, expr, metas, pool_cols


def apply_yazar_pool_rule(expr, metas, pool_cols):
    fe_expr = {}
    fe_meta = {}
    exclusions = {}
    for group, meta in metas.items():
        pool_counts = meta["pool_number"].value_counts()
        ge2 = pool_counts[pool_counts >= 2]
        keep = meta.index[meta["pool_number"].isin(ge2.index)].tolist()
        estimated_df = len(keep) - (4 + max(0, len(ge2) - 1))
        if len(ge2) < 2 or estimated_df < 10:
            exclusions[group] = {
                "n_donors": len(meta),
                "n_donors_ge2_pools": len(keep),
                "n_pools_ge2": len(ge2),
                "estimated_residual_df": estimated_df,
            }
            continue
        sub = meta.loc[keep].copy()
        reference = int(ge2.index[0])
        all_pools = sorted(int(column.split("_")[1]) for column in pool_cols)
        for pool in all_pools:
            sub[f"pool_{pool}"] = (
                (sub["pool_number"].astype(int) == pool).astype(float)
                if pool != reference else 0.0
            )
        fe_expr[group] = expr[group][keep]
        fe_meta[group] = sub
    return fe_expr, fe_meta, exclusions


def empty_followup_row(cohort: str, unit: str, tier: str, arm: str, programs) -> dict:
    return {
        "cohort": cohort,
        "unit": unit,
        "gene_set": tier,
        "arm": arm,
        "source": "Gao" if cohort == "Gao_Heart" else "Yazar",
        "model": "batch_adjusted" if arm == "primary" else "unadjusted",
        "status": "UNTESTABLE",
        "reason": "",
        "composite_slope": np.nan,
        "composite_slope_per_decade": np.nan,
        "composite_slope_per_year": np.nan,
        "composite_se": np.nan,
        "composite_pvalue": np.nan,
        "composite_q": np.nan,
        "empirical_pvalue": np.nan,
        "empirical_q": np.nan,
        "composite_pass": False,
        "empirical_pass": pd.NA,
        "both_pass": pd.NA,
        "n_donors": np.nan,
        "n_donors_30plus_cohort": np.nan,
        "n_genes_mapped": len(set(programs[tier]["up"]) | set(programs[tier]["down"])),
        "n_up_mapped": len(programs[tier]["up"]),
        "n_down_mapped": len(programs[tier]["down"]),
        "n_genes_detected": np.nan,
        "n_up_detected": np.nan,
        "n_down_detected": np.nan,
        "mde80_per_decade": np.nan,
        "residual_df": np.nan,
        "null_iterations": np.nan,
        "random_seed": np.nan,
        "multiplicity_family": f"mature_{cohort}_{arm}_{tier}_30plus",
        "age_scale": "per_decade",
        "min_fitted_age": np.nan,
        "max_fitted_age": np.nan,
        "covariates": "",
        "omitted_covariates": "",
    }


def result_rows(result: pd.DataFrame, cohort: str, tier: str, arm: str, programs, cohort_donors: int) -> list[dict]:
    rows = []
    for record in result.itertuples(index=False):
        row = empty_followup_row(cohort, record.group, tier, arm, programs)
        residual_df = int(record.residual_df)
        mde80 = (stats.t.ppf(0.975, residual_df) + stats.t.ppf(0.80, residual_df)) * float(record.composite_se)
        row.update({
            "status": "TESTED",
            "reason": "",
            "composite_slope": float(record.composite_effect),
            "composite_slope_per_decade": float(record.composite_effect),
            "composite_slope_per_year": float(record.composite_effect) / 10.0,
            "composite_se": float(record.composite_se),
            "composite_pvalue": float(record.composite_pvalue),
            "empirical_pvalue": float(record.empirical_null_pvalue),
            "n_donors": int(record.n_donors),
            "n_donors_30plus_cohort": int(cohort_donors),
            "n_genes_detected": int(record.n_p1_total_detected),
            "n_up_detected": int(record.n_p1_up_detected),
            "n_down_detected": int(record.n_p1_down_detected),
            "mde80_per_decade": float(mde80),
            "residual_df": residual_df,
            "null_iterations": int(record.null_iterations),
            "random_seed": SEED,
            "covariates": record.covariates,
            "omitted_covariates": record.omitted_covariates,
        })
        rows.append(row)
    return rows


def generic_untestable_reason(group: str, tier: str, arm: str, expr, metas, programs, covariates, min_df: int) -> tuple[str, int, float, float]:
    if group not in expr or group not in metas:
        return "no 30-plus pseudobulk group after the 20-cell donor floor", 0, np.nan, np.nan
    meta = metas[group]
    common = [donor for donor in meta.index if donor in expr[group].columns]
    n_donors = len(common)
    ages = meta.loc[common, "age_years"].astype(float) if common else pd.Series(dtype=float)
    age_min = float(ages.min()) if common else np.nan
    age_max = float(ages.max()) if common else np.nan
    if n_donors < 4:
        return f"only {n_donors} donors at 30-plus after technical gates; minimum is 4", n_donors, age_min, age_max
    matrix = expr[group][common]
    minimum_required = max(3, int(np.ceil(0.5 * n_donors)))
    keep = (matrix.to_numpy(dtype=float) >= np.log2(2.0)).sum(axis=1) >= minimum_required
    detected = matrix.loc[keep]
    if detected.shape[0] < 200:
        return f"only {detected.shape[0]} genes in the mature expression pool; minimum is 200", n_donors, age_min, age_max
    up = detected.index.isin(set(programs[tier]["up"])).sum()
    down = detected.index.isin(set(programs[tier]["down"])).sum()
    if up == 0 or down == 0:
        return f"composite is inestimable with {up} detected up genes and {down} detected down genes", n_donors, age_min, age_max
    try:
        select_design_columns(
            meta.loc[common], age_col="age_per_decade", candidate_covariates=tuple(covariates), min_residual_df=min_df
        )
    except ValueError as error:
        return f"30-plus design is ineligible: {error}", n_donors, age_min, age_max
    return f"UNRESOLVED: {tier} {arm} result missing despite passing preflight gates", n_donors, age_min, age_max


def run_followup_cohort(cohort: str, rosters, programs, expr, metas, cohort_donors: int,
                        primary_expr=None, primary_metas=None, primary_exclusions=None, pool_cols=None) -> pd.DataFrame:
    rows = []
    for arm in ("primary", "unadjusted"):
        roster = rosters[(cohort, arm)]
        if arm == "primary":
            arm_expr = primary_expr if primary_expr is not None else expr
            arm_metas = primary_metas if primary_metas is not None else metas
            if cohort == "Gao_Heart":
                covariates = ("sex_male", "log_lib", *STUDY_DUMMY_MAP.values())
                min_df = 2
            else:
                covariates = ("sex_male", "log_lib", *pool_cols)
                min_df = 10
        else:
            arm_expr = expr
            arm_metas = metas
            covariates = ("sex_male", "log_lib")
            min_df = 2
        for tier in TIERS:
            print(f"{cohort} {arm} {tier}: running 2000-draw matched nulls", flush=True)
            result = run_family_c_p1_signature(
                arm_expr, arm_metas, p1_program=programs[tier], candidate_covariates=tuple(covariates),
                min_donors=4, min_residual_df=min_df, iterations=ITERATIONS, seed=SEED,
            )
            tier_rows = result_rows(result, cohort, tier, arm, programs, cohort_donors)
            by_unit = {row["unit"]: row for row in tier_rows}
            for unit in roster:
                if unit in by_unit:
                    row = by_unit[unit]
                    meta = arm_metas[unit]
                    donors = list(arm_expr[unit].columns)
                    row["min_fitted_age"] = float(meta.loc[donors, "age_years"].min())
                    row["max_fitted_age"] = float(meta.loc[donors, "age_years"].max())
                    rows.append(row)
                    continue
                row = empty_followup_row(cohort, unit, tier, arm, programs)
                row["n_donors_30plus_cohort"] = int(cohort_donors)
                if arm == "primary" and primary_exclusions and unit in primary_exclusions:
                    detail = primary_exclusions[unit]
                    reason = (
                        f"pool pruning retained {detail['n_donors_ge2_pools']} donors in "
                        f"{detail['n_pools_ge2']} pools with at least 2 donors; estimated residual df "
                        f"{detail['estimated_residual_df']} is below 10 or fewer than 2 pools remain"
                    )
                    n_donors = detail["n_donors"]
                    ages = metas[unit]["age_years"].astype(float) if unit in metas else pd.Series(dtype=float)
                else:
                    reason, n_donors, age_min, age_max = generic_untestable_reason(
                        unit, tier, arm, arm_expr, arm_metas, programs, covariates, min_df
                    )
                    ages = pd.Series([age_min, age_max], dtype=float).dropna()
                assert not reason.startswith("UNRESOLVED"), f"{cohort} {unit} {tier} {arm}: {reason}"
                row["reason"] = reason
                row["n_donors"] = int(n_donors)
                if len(ages):
                    row["min_fitted_age"] = np.nan
                    row["max_fitted_age"] = np.nan
                rows.append(row)
    output = pd.DataFrame(rows)
    for (_, _, _), indices in output.groupby(["cohort", "arm", "gene_set"], sort=True).groups.items():
        tested = output.loc[indices]
        tested = tested.index[tested["status"] == "TESTED"]
        if len(tested) == 0:
            continue
        output.loc[tested, "composite_q"] = multipletests(
            output.loc[tested, "composite_pvalue"], method="fdr_bh"
        )[1]
        output.loc[tested, "empirical_q"] = multipletests(
            output.loc[tested, "empirical_pvalue"], method="fdr_bh"
        )[1]
    output["composite_pass"] = output["composite_q"] < ALPHA
    output["empirical_pass"] = output["empirical_q"].lt(ALPHA).astype("boolean")
    output["both_pass"] = output["composite_pass"].astype("boolean") & output["empirical_pass"]
    return output


def build_coverage() -> pd.DataFrame:
    rows = []
    files = {
        "P0m": "followup-mature-p1m-p0m-genes.csv",
        "P1m": "followup-mature-p1m-genes.csv",
    }
    for tier, filename in files.items():
        frame = pd.read_csv(RESULTS / filename)
        frame["mapped"] = frame["human_ortholog"].notna()
        frame["biotype_class"] = np.where(
            frame["biotype"] == "protein_coding", "protein_coding", "noncoding_or_unknown"
        )
        for label, sub in [("all", frame), *list(frame.groupby("biotype_class", sort=True))]:
            mapped = sub[sub["mapped"]]
            rows.append({
                "gene_set": tier,
                "biotype_class": label,
                "n_mouse_genes": int(sub["mouse_gene_symbol"].nunique()),
                "n_mapped": int(mapped["human_ortholog"].nunique()),
                "n_unmapped": int((~sub["mapped"]).sum()),
                "mapping_coverage": float(sub["mapped"].mean()),
                "n_mouse_up": int((sub["direction"] == "UP").sum()),
                "n_mouse_down": int((sub["direction"] == "DOWN").sum()),
                "n_human_up": int((mapped["direction"] == "UP").sum()),
                "n_human_down": int((mapped["direction"] == "DOWN").sum()),
                "unmapped_genes": ";".join(sorted(sub.loc[~sub["mapped"], "mouse_gene_symbol"])),
            })
    return pd.DataFrame(rows)


def status_label(composite_pass, empirical_pass, has_empirical=True, tested=True) -> str:
    if not tested:
        return "UNTESTABLE"
    if not has_empirical:
        return "COMPOSITE_PASS_NO_EMPIRICAL" if composite_pass else "NO_EMPIRICAL"
    if composite_pass and empirical_pass:
        return "BOTH_PASS"
    if empirical_pass:
        return "EMPIRICAL_ONLY"
    if composite_pass:
        return "COMPOSITE_ONLY"
    return "FAIL_BOTH"


def tier_values(row: pd.Series, prefix: str) -> dict:
    tested = row["status"] == "TESTED"
    has_empirical = tested and pd.notna(row["empirical_pvalue"])
    composite_pass = bool(tested and row["composite_q"] < ALPHA)
    empirical_pass = bool(has_empirical and row["empirical_q"] < ALPHA)
    both = bool(composite_pass and empirical_pass)
    fails_both = bool(has_empirical and not composite_pass and not empirical_pass)
    return {
        f"{prefix}_status": status_label(composite_pass, empirical_pass, has_empirical, tested),
        f"{prefix}_reason": "" if tested else row["reason"],
        f"{prefix}_n_donors": row["n_donors"],
        f"{prefix}_age_min": row.get("age_min", row.get("min_fitted_age", np.nan)),
        f"{prefix}_age_max": row.get("age_max", row.get("max_fitted_age", np.nan)),
        f"{prefix}_n_genes_detected": row["n_genes_detected"],
        f"{prefix}_n_up_detected": row["n_up_detected"],
        f"{prefix}_n_down_detected": row["n_down_detected"],
        f"{prefix}_composite_slope": row["composite_slope"],
        f"{prefix}_composite_slope_per_year": row["composite_slope_per_year"],
        f"{prefix}_composite_se": row.get("composite_se", np.nan),
        f"{prefix}_composite_pvalue": row["composite_pvalue"],
        f"{prefix}_composite_q": row["composite_q"],
        f"{prefix}_composite_pass": composite_pass if tested else pd.NA,
        f"{prefix}_empirical_pvalue": row["empirical_pvalue"],
        f"{prefix}_empirical_q": row["empirical_q"],
        f"{prefix}_empirical_pass": empirical_pass if has_empirical else pd.NA,
        f"{prefix}_both_pass": both if has_empirical else pd.NA,
        f"{prefix}_fails_both_strict": fails_both,
        f"{prefix}_mde80_per_decade": row.get("mde80_per_decade", np.nan),
        f"{prefix}_null_iterations": row["null_iterations"],
        f"{prefix}_random_seed": row["random_seed"],
        f"{prefix}_multiplicity_family": row["multiplicity_family"],
    }


def build_comparison(followup: pd.DataFrame) -> pd.DataFrame:
    core_transport = pd.read_csv(CORE_TRANSPORT)
    core_null = pd.read_csv(CORE_NULL)
    core = core_transport.merge(
        core_null[["cohort", "unit", "gene_set", "empirical_pvalue", "empirical_q", "null_iterations", "random_seed"]],
        on=["cohort", "unit", "gene_set"], how="left", validate="one_to_one"
    )
    core["arm"] = "primary"
    core["reason"] = ""
    core["mde80_per_decade"] = np.nan
    core["composite_se"] = np.nan
    core_untestable = pd.read_csv(CORE_UNTESTABLE)
    family = pd.read_csv(ORIGINAL_FAMILY)[["cohort", "cell_type", "source"]].rename(columns={"cell_type": "unit"})
    followup_primary = followup[followup["arm"] == "primary"].copy()
    followup_primary["source"] = followup_primary["source"].map({"Gao": "Gao", "Yazar": "Yazar"})
    all_long = pd.concat([core, followup_primary], ignore_index=True, sort=False)
    all_long["age_min"] = all_long["age_min"].combine_first(all_long.get("min_fitted_age"))
    all_long["age_max"] = all_long["age_max"].combine_first(all_long.get("max_fitted_age"))
    rows = []
    for tier in TIERS:
        tier_core = core[core["gene_set"] == tier]
        tested_keys = set(zip(tier_core["cohort"], tier_core["unit"]))
        for record in family.itertuples(index=False):
            if (record.cohort, record.unit) in tested_keys:
                continue
            missing = core_untestable[
                (core_untestable["cohort"] == record.cohort) & (core_untestable["unit"] == record.unit)
            ]
            assert len(missing) == 1
            missing = missing.iloc[0]
            rows.append({
                "cohort": record.cohort,
                "unit": record.unit,
                "gene_set": tier,
                "source": record.source,
                "model": "not_fitted",
                "status": "UNTESTABLE",
                "reason": missing["reason"],
                "composite_slope": np.nan,
                "composite_slope_per_year": np.nan,
                "composite_se": np.nan,
                "composite_pvalue": np.nan,
                "composite_q": np.nan,
                "empirical_pvalue": np.nan,
                "empirical_q": np.nan,
                "n_donors": missing["n_donors_30plus"],
                "n_genes_detected": np.nan,
                "n_up_detected": np.nan,
                "n_down_detected": np.nan,
                "mde80_per_decade": np.nan,
                "null_iterations": np.nan,
                "random_seed": np.nan,
                "multiplicity_family": "mature_P1m_P0m_M9_M10_30plus",
                "age_min": missing["age_min_30plus"],
                "age_max": missing["age_max_30plus"],
            })
    all_long = pd.concat([all_long, pd.DataFrame(rows)], ignore_index=True, sort=False)
    output = []
    for (cohort, unit), sub in all_long.groupby(["cohort", "unit"], sort=True):
        assert len(sub) == 2 and set(sub["gene_set"]) == set(TIERS), f"bad pair {cohort} {unit}"
        p0 = sub[sub["gene_set"] == "P0m"].iloc[0]
        p1 = sub[sub["gene_set"] == "P1m"].iloc[0]
        p0_values = tier_values(p0, "p0m")
        p1_values = tier_values(p1, "p1m")
        row = {
            "cohort": cohort,
            "unit": unit,
            "source": p0["source"],
            "status": "TESTED" if p0["status"] == "TESTED" and p1["status"] == "TESTED" else "UNTESTABLE",
            **p0_values,
            **p1_values,
        }
        row["strict_p0m_rescue"] = bool(row["p0m_both_pass"] is True and row["p1m_fails_both_strict"])
        output.append(row)
    comparison = pd.DataFrame(output)
    expected = 252 + 14 + 17
    assert len(comparison) == expected
    assert int((comparison["status"] == "TESTED").sum()) == 170 + len(followup_primary[followup_primary["status"] == "TESTED"].query("gene_set == 'P0m'"))
    comparison["strict_p0m_rescue"] = comparison["p0m_both_pass"].fillna(False).astype(bool) & comparison["p1m_fails_both_strict"]
    return comparison


def verify_core(comparison: pd.DataFrame) -> dict:
    committed = pd.read_csv(CORE_COMPARISON)
    keys = ["cohort", "unit"]
    current = comparison[comparison["p0m_empirical_pvalue"].notna() & comparison["source"].isin(["M9", "M10"])].copy()
    assert len(current) == 48
    joined = committed.merge(current, on=keys, suffixes=("_committed", "_new"), validate="one_to_one")
    assert len(joined) == 48
    for tier in ("p0m", "p1m"):
        for column in ("composite_pvalue", "composite_q", "empirical_pvalue", "empirical_q"):
            assert np.allclose(joined[f"{tier}_{column}_committed"], joined[f"{tier}_{column}_new"], rtol=1e-12, atol=1e-15)
        assert (joined[f"{tier}_both_pass_committed"].astype(bool) == joined[f"{tier}_both_pass_new"].fillna(False).astype(bool)).all()
    anchors = {
        "core_units": len(committed),
        "p0m_both_pass": int(committed["p0m_both_pass"].sum()),
        "p1m_both_pass": int(committed["p1m_both_pass"].sum()),
        "p1m_fail_both": int(committed["p1m_fails_both_strict"].sum()),
        "strict_rescue": int(committed["strict_p0m_rescue"].sum()),
    }
    assert anchors == {"core_units": 48, "p0m_both_pass": 3, "p1m_both_pass": 2, "p1m_fail_both": 29, "strict_rescue": 0}
    current_anchors = {
        "core_units": len(current),
        "p0m_both_pass": int(current["p0m_both_pass"].fillna(False).sum()),
        "p1m_both_pass": int(current["p1m_both_pass"].fillna(False).sum()),
        "p1m_fail_both": int(current["p1m_fails_both_strict"].sum()),
        "strict_rescue": int(current["strict_p0m_rescue"].sum()),
    }
    assert current_anchors == anchors
    return anchors


def verify_followup(followup: pd.DataFrame) -> dict:
    discrepancies = []
    for (_, _, _), sub in followup[followup["status"] == "TESTED"].groupby(["cohort", "arm", "gene_set"]):
        for p_column, q_column in (("composite_pvalue", "composite_q"), ("empirical_pvalue", "empirical_q")):
            recomputed = multipletests(sub[p_column], method="fdr_bh")[1]
            discrepancies.extend(np.abs(recomputed - sub[q_column].to_numpy(float)).tolist())
    tested = followup[followup["status"] == "TESTED"]
    assert (tested["null_iterations"] == ITERATIONS).all()
    assert (tested["random_seed"] == SEED).all()
    assert (tested["min_fitted_age"] >= MIN_AGE).all()
    expected = {("Gao_Heart", "primary"): 14, ("Gao_Heart", "unadjusted"): 14,
                ("Yazar_Blood", "primary"): 17, ("Yazar_Blood", "unadjusted"): 20}
    for key, n_units in expected.items():
        sub = followup[(followup["cohort"] == key[0]) & (followup["arm"] == key[1])]
        assert len(sub) == 2 * n_units
        assert sub.groupby("gene_set").size().eq(n_units).all()
    return {
        "max_abs_bh_discrepancy": float(max(discrepancies, default=0.0)),
        "gao_min_fitted_age": float(tested.loc[tested["cohort"] == "Gao_Heart", "min_fitted_age"].min()),
        "yazar_min_fitted_age": float(tested.loc[tested["cohort"] == "Yazar_Blood", "min_fitted_age"].min()),
        "null_rows": len(tested),
    }


def summarize(comparison: pd.DataFrame) -> pd.DataFrame:
    rows = []
    n_rescue = int(comparison["strict_p0m_rescue"].sum())
    verdict = "PROVES" if n_rescue else "DOES NOT PROVE"
    for tier in ("p0m", "p1m"):
        tested = comparison[comparison["status"] == "TESTED"]
        null = tested[tested[f"{tier}_empirical_pvalue"].notna()]
        statuses = null[f"{tier}_status"]
        rows.append({
            "gene_set": {"p0m": "P0m", "p1m": "P1m"}[tier],
            "original_design_units": len(comparison),
            "primary_units": len(tested),
            "untestable_units": int((comparison["status"] == "UNTESTABLE").sum()),
            "nominal_composite_hits": int((tested[f"{tier}_composite_pvalue"] < ALPHA).sum()),
            "composite_passes": int(tested[f"{tier}_composite_pass"].fillna(False).sum()),
            "null_evaluated": len(null),
            "random_passes": int(null[f"{tier}_empirical_pass"].fillna(False).sum()),
            "both_pass": int(null[f"{tier}_both_pass"].fillna(False).sum()),
            "empirical_only": int((statuses == "EMPIRICAL_ONLY").sum()),
            "composite_only": int((statuses == "COMPOSITE_ONLY").sum()),
            "fail_both": int((statuses == "FAIL_BOTH").sum()),
            "p1m_fails_both_strict": int(comparison["p1m_fails_both_strict"].sum()),
            "strict_rescue_count": n_rescue,
            "strict_verdict": verdict,
        })
    return pd.DataFrame(rows)


def fmt(value, digits=6) -> str:
    if value is None or pd.isna(value):
        return "NA"
    if isinstance(value, (bool, np.bool_)):
        return "True" if value else "False"
    return f"{float(value):.{digits}g}"


def exact(value) -> str:
    if value is None or pd.isna(value):
        return "NA"
    return f"{float(value):.17g}"


def report_counts(comparison: pd.DataFrame, cohort: str, tier: str) -> dict:
    sub = comparison[(comparison["cohort"] == cohort) & (comparison["status"] == "TESTED")]
    null = sub[sub[f"{tier}_empirical_pvalue"].notna()]
    return {
        "units": len(sub),
        "nominal": int((sub[f"{tier}_composite_pvalue"] < ALPHA).sum()),
        "composite": int(sub[f"{tier}_composite_pass"].fillna(False).sum()),
        "null": len(null),
        "random": int(null[f"{tier}_empirical_pass"].fillna(False).sum()),
        "both": int(null[f"{tier}_both_pass"].fillna(False).sum()),
    }


def write_report(comparison: pd.DataFrame, followup: pd.DataFrame, coverage: pd.DataFrame,
                 summary: pd.DataFrame, original: pd.DataFrame, verification: dict,
                 core_anchors: dict, cohort_donors: dict) -> None:
    verdict = summary["strict_verdict"].iloc[0]
    fail_both = int(comparison["p1m_fails_both_strict"].sum())
    rescue = int(comparison["strict_p0m_rescue"].sum())
    null_eval = int(comparison["p0m_empirical_pvalue"].notna().sum())
    lines = [
        "# Mature Block 5: P0m versus P1m filter comparison",
        "",
        "## Spec",
        "",
        "This is the full mature-to-old extension of frozen Block 5. Mouse derivation used exactly 18-versus-30-month animals and every fitted human donor in this report is at least 30 years old. P0m and P1m are the committed human-mapped probe sets, loaded by the same load_programs function used by the mature core implementation. P1m is not a subset of P0m, so their coverage and transport results are evaluated independently.",
        "",
        "The committed 170-unit M9/M10 composite family, 48-unit MSSM/PsychAD matched-null comparison, and 82 core untestable dispositions from commit 3f5b57a are reused without refitting or recomputing core nulls. New work is restricted to Gao heart and Yazar blood. Both P0m and P1m primaries are computed fresh. Gao uses sex, depth, and the original six study indicators, at least 4 donors, residual df at least 2, and fresh 30-plus pseudobulk aggregation. Yazar uses the original cache, the same pool-pruning rule, residual df floor 10 for the primary and 2 for the unadjusted arm. Unadjusted arms retain sex and depth only.",
        "",
        f"Every tested Gao/Yazar unit and probe set has {ITERATIONS} expression-decile-matched null draws at seed {SEED}. Composite and empirical BH correction is separate within each cohort, arm, and probe set. Core composite and empirical corrections remain the committed 170-unit, 21-unit MSSM, and 27-unit PsychAD families. No empirical family was added for other M9 units.",
        "",
        "The strict verdict is evaluated across all null-evaluated units: the 48 core units plus tested Gao and Yazar units. P0m must pass both corrected tests where P1m fails both. A unit without an empirical result cannot qualify and is not counted as a P1m failure. UNTESTABLE records an eligibility or input failure and is not a biological null.",
        "",
        "## Mapping coverage",
        "",
        "| Probe set | Biotype class | Mouse genes | Mapped | Unmapped | Coverage | Mapped up | Mapped down |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in coverage.itertuples(index=False):
        lines.append(f"| {row.gene_set} | {row.biotype_class} | {row.n_mouse_genes} | {row.n_mapped} | {row.n_unmapped} | {100*row.mapping_coverage:.2f}% | {row.n_human_up} | {row.n_human_down} |")
    lines.extend([
        "",
        "## Global comparison",
        "",
        "Primary units are tested units. Original design units also include all explicitly untestable units. Empirical-only, composite-only, and fail-both counts are among null-evaluated units only.",
        "",
        "| Probe set | Original design | Primary tested | Untestable | Nominal composite | Composite pass | Null evaluated | Random pass | Both pass | Empirical only | Composite only | Fail both |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ])
    for row in summary.itertuples(index=False):
        lines.append(f"| {row.gene_set} | {row.original_design_units} | {row.primary_units} | {row.untestable_units} | {row.nominal_composite_hits} | {row.composite_passes} | {row.null_evaluated} | {row.random_passes} | {row.both_pass} | {row.empirical_only} | {row.composite_only} | {row.fail_both} |")
    lines.extend([
        "",
        "## Cohort comparison",
        "",
        "| Cohort | Probe set | Tested | Untestable | Nominal composite | Composite pass | Null evaluated | Random pass | Both pass |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ])
    for cohort in sorted(comparison["cohort"].unique()):
        for tier in ("p0m", "p1m"):
            counts = report_counts(comparison, cohort, tier)
            untestable = int(((comparison["cohort"] == cohort) & (comparison["status"] == "UNTESTABLE")).sum())
            tier_label = {"p0m": "P0m", "p1m": "P1m"}[tier]
            lines.append(f"| {cohort} | {tier_label} | {counts['units']} | {untestable} | {counts['nominal']} | {counts['composite']} | {counts['null']} | {counts['random']} | {counts['both']} |")
    interesting = comparison[
        comparison["p0m_composite_pass"].fillna(False) | comparison["p1m_composite_pass"].fillna(False)
        | (comparison["p0m_both_pass"].fillna(False) != comparison["p1m_both_pass"].fillna(False))
        | comparison["strict_p0m_rescue"]
    ]
    lines.extend([
        "",
        "## Side-by-side units with a corrected composite pass or verdict change",
        "",
        "Slopes are per year for M9/M10 and per decade for Gao/Yazar. Exact q values are shown.",
        "",
        "| Cohort | Unit | P0m slope | P0m composite q | P0m empirical q | P0m both | P1m slope | P1m composite q | P1m empirical q | P1m both |",
        "| --- | --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | --- |",
    ])
    for row in interesting.sort_values(["cohort", "unit"]).itertuples(index=False):
        lines.append(f"| {row.cohort} | {row.unit} | {fmt(row.p0m_composite_slope)} | {exact(row.p0m_composite_q)} | {exact(row.p0m_empirical_q)} | {fmt(row.p0m_both_pass)} | {fmt(row.p1m_composite_slope)} | {exact(row.p1m_composite_q)} | {exact(row.p1m_empirical_q)} | {fmt(row.p1m_both_pass)} |")
    disagreement = comparison[comparison["p0m_both_pass"].fillna(False) != comparison["p1m_both_pass"].fillna(False)]
    lines.extend([
        "",
        "## Both-pass disagreements",
        "",
        f"P0m and P1m disagree on both-pass status in {len(disagreement)} units. P0m gains status in {int((disagreement['p0m_both_pass'].fillna(False) & ~disagreement['p1m_both_pass'].fillna(False)).sum())} and loses it in {int((~disagreement['p0m_both_pass'].fillna(False) & disagreement['p1m_both_pass'].fillna(False)).sum())}.",
        "",
        "| Cohort | Unit | P0m composite q | P0m empirical q | P0m both | P1m composite q | P1m empirical q | P1m both |",
        "| --- | --- | ---: | ---: | --- | ---: | ---: | --- |",
    ])
    for row in disagreement.sort_values(["cohort", "unit"]).itertuples(index=False):
        lines.append(f"| {row.cohort} | {row.unit} | {exact(row.p0m_composite_q)} | {exact(row.p0m_empirical_q)} | {fmt(row.p0m_both_pass)} | {exact(row.p1m_composite_q)} | {exact(row.p1m_empirical_q)} | {fmt(row.p1m_both_pass)} |")
    lines.extend([
        "",
        "## Unadjusted sensitivity",
        "",
        "| Cohort | Probe set | Tested | Untestable | Nominal composite | Composite pass | Random pass | Both pass |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ])
    unadjusted = followup[followup["arm"] == "unadjusted"]
    for (cohort, tier), sub in unadjusted.groupby(["cohort", "gene_set"], sort=True):
        tested = sub[sub["status"] == "TESTED"]
        lines.append(f"| {cohort} | {tier} | {len(tested)} | {int((sub['status'] == 'UNTESTABLE').sum())} | {int((tested['composite_pvalue'] < ALPHA).sum())} | {int(tested['composite_pass'].sum())} | {int(tested['empirical_pass'].fillna(False).sum())} | {int(tested['both_pass'].fillna(False).sum())} |")
    lines.extend([
        "",
        "## MDE80 and power",
        "",
        "MDE80 per decade is (t_0.975 + t_0.80) times the composite standard error at the fitted residual df. The committed core transport CSV has no standard-error column, so no core MDE80 is fabricated.",
        "",
        "| Cohort | Probe set | Arm | Tested | MDE80 median | MDE80 minimum | MDE80 maximum |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: |",
    ])
    for (cohort, tier, arm), sub in followup[followup["status"] == "TESTED"].groupby(["cohort", "gene_set", "arm"], sort=True):
        values = sub["mde80_per_decade"].astype(float)
        lines.append(f"| {cohort} | {tier} | {arm} | {len(sub)} | {fmt(values.median())} | {fmt(values.min())} | {fmt(values.max())} |")
    lines.extend([
        "",
        "| Cohort | Original selected donors | Donors at 30-plus | Minimum fitted age | Original primary tested | Mature primary tested | Mature primary untestable |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ])
    for cohort in ("Gao_Heart", "Yazar_Blood"):
        info = cohort_donors[cohort]
        primary = followup[(followup["cohort"] == cohort) & (followup["arm"] == "primary") & (followup["gene_set"] == "P1m")]
        original_n = len(original[(original["cohort"] == cohort) & (original["arm"] == "primary") & (original["gene_set"] == "P1_Clean")])
        lines.append(f"| {cohort} | {info['original']} | {info['mature']} | {fmt(info['min_fitted_age'])} | {original_n} | {int((primary['status'] == 'TESTED').sum())} | {int((primary['status'] == 'UNTESTABLE').sum())} |")
    for cohort in ("Gao_Heart", "Yazar_Blood"):
        info = cohort_donors[cohort]
        removed = info["original"] - info["mature"]
        primary = followup[(followup["cohort"] == cohort) & (followup["arm"] == "primary") & (followup["gene_set"] == "P1m")]
        mde = primary.loc[primary["status"] == "TESTED", "mde80_per_decade"].astype(float)
        lines.extend([
            "",
            f"{cohort} loses {removed} of {info['original']} selected donors at the 30-plus cutoff. Primary MDE80 has median {fmt(mde.median()) if len(mde) else 'NA'} and range {fmt(mde.min()) if len(mde) else 'NA'} to {fmt(mde.max()) if len(mde) else 'NA'} per decade. {info['power_statement']}",
        ])
    lines.extend([
        "",
        "## Strict verdict",
        "",
        f"**{verdict}**. Across {null_eval} null-evaluated units, P1m fails both tests in {fail_both}. P0m passes both tests in {rescue} of those units. Core anchors are {core_anchors['p1m_fail_both']} P1m fail-both units and {core_anchors['strict_rescue']} rescues among {core_anchors['core_units']} units. The additional Gao/Yazar units contribute {fail_both-core_anchors['p1m_fail_both']} P1m fail-both units and {rescue-core_anchors['strict_rescue']} rescues.",
        "",
        "## Survived or died relative to original Block 5",
        "",
        "SURVIVED means the same cohort and unit passes both tests with the corresponding mature probe set. DIED means it no longer passes both; an untestable mature unit is labeled DIED_UNTESTABLE and is not treated as a null.",
        "",
    ])
    original_primary = original[original["arm"] == "primary"].copy()
    original_p0 = original_primary[original_primary["gene_set"] == "P0_Raw257"]
    original_p1 = original_primary[original_primary["gene_set"] == "P1_Clean"]
    assert len(original_primary[original_primary["gene_set"] == "P0_Raw257"]) == 283
    old_p0_both = original_p0[original_p0["both_pass"].fillna(False)]
    old_p1_both = original_p1[original_p1["both_pass"].fillna(False)]
    assert len(old_p0_both) == 10 and len(old_p1_both) == 7
    old_fail = original_p1[
        original_p1["empirical_pvalue"].notna() & ~original_p1["composite_pass"] & ~original_p1["empirical_pass"].fillna(True)
    ]
    old_rescue = old_fail.merge(
        original_p0[["cohort", "unit", "both_pass"]], on=["cohort", "unit"], suffixes=("_p1", "_p0")
    )
    assert len(old_fail) == 59 and int(old_rescue["both_pass_p0"].fillna(False).sum()) == 0
    lines.append(f"The original full-range design had {len(old_p0_both)} P0 both-pass and {len(old_p1_both)} P1 both-pass units among 283, with 0 strict rescues among {len(old_fail)} P1 fail-both units.")
    lines.extend([
        "",
        "| Probe set | Cohort | Unit | Original composite q | Original empirical q | Mature composite q | Mature empirical q | Disposition |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | --- |",
    ])
    for old_tier, new_prefix, old_both in (("P0", "p0m", old_p0_both), ("P1", "p1m", old_p1_both)):
        for old in old_both.sort_values(["cohort", "unit"]).itertuples(index=False):
            match = comparison[(comparison["cohort"] == old.cohort) & (comparison["unit"] == old.unit)]
            assert len(match) == 1
            new = match.iloc[0]
            if new["status"] == "UNTESTABLE":
                disposition = "DIED_UNTESTABLE"
            else:
                disposition = "SURVIVED" if bool(new[f"{new_prefix}_both_pass"]) else "DIED"
            lines.append(f"| {old_tier} | {old.cohort} | {old.unit} | {exact(old.composite_q)} | {exact(old.empirical_q)} | {exact(new[f'{new_prefix}_composite_q'])} | {exact(new[f'{new_prefix}_empirical_q'])} | {disposition} |")
    untestable_followup = followup[(followup["arm"] == "primary") & (followup["status"] == "UNTESTABLE")]
    lines.extend([
        "",
        "## Untestable Gao and Yazar primary units",
        "",
        "| Cohort | Unit | Probe set | Donors at 30-plus | Reason |",
        "| --- | --- | --- | ---: | --- |",
    ])
    if untestable_followup.empty:
        lines.append("| none | none | none | 0 | none |")
    else:
        for row in untestable_followup.sort_values(["cohort", "unit", "gene_set"]).itertuples(index=False):
            lines.append(f"| {row.cohort} | {row.unit} | {row.gene_set} | {row.n_donors} | {row.reason} |")
    lines.extend([
        "",
        "The comparison CSV separately retains all 82 committed core untestable units and their original reasons.",
        "",
        "## Verification",
        "",
        f"Independent statsmodels BH recomputation within every stated Gao/Yazar cohort, arm, and probe-set family gives maximum absolute q discrepancy {verification['max_abs_bh_discrepancy']:.17g}. All {verification['null_rows']} tested Gao/Yazar transport rows have {ITERATIONS} null iterations and seed {SEED}. Minimum fitted ages are {fmt(verification['gao_min_fitted_age'])} for Gao and {fmt(verification['yazar_min_fitted_age'])} for Yazar.",
        "",
        f"The rebuilt full comparison reproduces the committed core values exactly: P1m both-pass {core_anchors['p1m_both_pass']}, P0m both-pass {core_anchors['p0m_both_pass']}, P1m fail-both {core_anchors['p1m_fail_both']}, and strict rescue {core_anchors['strict_rescue']} in {core_anchors['core_units']} core units.",
        "",
        "## Files",
        "",
        "- `scripts/run_mature_block5_gao_yazar.py`: new implementation; original and mature core scripts are unchanged.",
        "- `results/followup-mature-block5-gao-yazar-transport.csv`: all Gao/Yazar primary and unadjusted rows, including untestable dispositions.",
        "- `results/followup-mature-block5-p0m-comparison.csv`: complete 283-unit side-by-side comparison.",
        "- `results/followup-mature-block5-coverage.csv`: independent P0m and P1m mapping coverage.",
        "- `results/followup-mature-block5-summary.csv`: global counts and strict verdict.",
        "- `results/followup-mature-block5-runtime.csv`: stage runtimes and null settings.",
        "- `scratch/mature-blocks-458/mature_block5_gao_30plus_inputs.pkl`: large mature Gao pseudobulk cache.",
        "",
    ])
    REPORT.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    overall_start = time.monotonic()
    runtimes = []

    def stage_start(name):
        return {"stage": name, "model": MODEL, "start_utc": utc_now(), "start": time.monotonic()}

    def stage_end(timer):
        runtimes.append({
            "model": timer["model"],
            "stage": timer["stage"],
            "start_utc": timer["start_utc"],
            "end_utc": utc_now(),
            "duration_seconds": time.monotonic() - timer["start"],
            "iterations": ITERATIONS,
            "seed": SEED,
        })

    timer = stage_start("load_inputs")
    original, rosters = original_rosters()
    programs = load_programs()
    coverage = build_coverage()
    stage_end(timer)

    timer = stage_start("gao_30plus_transport")
    gao = load_gao_inputs(rosters[("Gao_Heart", "primary")])
    gao_output = run_followup_cohort(
        "Gao_Heart", rosters, programs, gao["expr"], gao["metas"], len(gao["donor_meta"])
    )
    del gao
    gc.collect()
    stage_end(timer)

    wait_for_memory("Yazar cache load")
    timer = stage_start("yazar_30plus_transport")
    primary_roster = rosters[("Yazar_Blood", "primary")]
    unadjusted_roster = rosters[("Yazar_Blood", "unadjusted")]
    donor_meta, cell_counts, expr, metas, pool_cols = pool_inputs_30plus(unadjusted_roster)
    primary_expr_all = {key: value for key, value in expr.items() if key in primary_roster}
    primary_metas_all = {key: value for key, value in metas.items() if key in primary_roster}
    fe_expr, fe_meta, exclusions = apply_yazar_pool_rule(primary_expr_all, primary_metas_all, pool_cols)
    yazar_output = run_followup_cohort(
        "Yazar_Blood", rosters, programs, expr, metas, len(donor_meta),
        primary_expr=fe_expr, primary_metas=fe_meta, primary_exclusions=exclusions, pool_cols=pool_cols,
    )
    del donor_meta, cell_counts, expr, metas, primary_expr_all, primary_metas_all, fe_expr, fe_meta
    gc.collect()
    stage_end(timer)

    timer = stage_start("comparison_and_verification")
    followup = pd.concat([gao_output, yazar_output], ignore_index=True, sort=False)
    followup = followup.sort_values(["cohort", "arm", "unit", "gene_set"]).reset_index(drop=True)
    comparison = build_comparison(followup)
    summary = summarize(comparison)
    core_anchors = verify_core(comparison)
    verification = verify_followup(followup)
    null_eval = comparison[comparison["p0m_empirical_pvalue"].notna()]
    fail_both = int(null_eval["p1m_fails_both_strict"].sum())
    rescue = int(null_eval["strict_p0m_rescue"].sum())
    assert fail_both == int(comparison["p1m_fails_both_strict"].sum())
    assert rescue == int(comparison["strict_p0m_rescue"].sum())
    assert (summary["strict_verdict"] == ("PROVES" if rescue else "DOES NOT PROVE")).all()
    cohort_donors = {}
    for cohort, original_n, mature_n in (
        ("Gao_Heart", 117, int(gao_output["n_donors_30plus_cohort"].iloc[0])),
        ("Yazar_Blood", 981, int(yazar_output["n_donors_30plus_cohort"].iloc[0])),
    ):
        primary = followup[(followup["cohort"] == cohort) & (followup["arm"] == "primary")]
        tested = primary[primary["status"] == "TESTED"]
        untestable = primary[primary["status"] == "UNTESTABLE"]
        if cohort == "Gao_Heart":
            statement = "The smaller donor pool reduces precision and any lost units are untestable rather than null."
        else:
            statement = "Most donors remain, but pool pruning can still remove units and inflate MDE80 in sparse cell types."
        cohort_donors[cohort] = {
            "original": original_n,
            "mature": mature_n,
            "min_fitted_age": float(tested["min_fitted_age"].min()),
            "power_statement": statement,
        }
    stage_end(timer)

    timer = stage_start("write_outputs")
    followup.to_csv(TRANSPORT_CSV, index=False)
    comparison.to_csv(COMPARISON_CSV, index=False)
    coverage.to_csv(COVERAGE_CSV, index=False)
    summary.to_csv(SUMMARY_CSV, index=False)
    write_report(comparison, followup, coverage, summary, original, verification, core_anchors, cohort_donors)
    stage_end(timer)

    total_end = utc_now()
    runtimes.append({
        "model": MODEL,
        "stage": "total",
        "start_utc": runtimes[0]["start_utc"],
        "end_utc": total_end,
        "duration_seconds": time.monotonic() - overall_start,
        "iterations": ITERATIONS,
        "seed": SEED,
    })
    pd.DataFrame(runtimes).to_csv(RUNTIME_CSV, index=False)
    print(summary.to_string(index=False), flush=True)
    print(f"Strict verdict: {summary['strict_verdict'].iloc[0]}; fail-both {fail_both}; rescue {rescue}", flush=True)
    print(f"BH maximum absolute discrepancy: {verification['max_abs_bh_discrepancy']:.17g}", flush=True)
    print(f"Wrote {TRANSPORT_CSV}", flush=True)
    print(f"Wrote {COMPARISON_CSV}", flush=True)
    print(f"Wrote {REPORT}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
