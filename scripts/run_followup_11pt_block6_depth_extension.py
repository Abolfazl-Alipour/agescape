#!/usr/bin/env python3
"""Block 6 donor-level depth versus age models for Gao, Yazar, and MSSM."""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from scripts.run_rev2_de_interneuron import _age_from_stage, _cat_eq
from src.followup_organs import write_clean_csv, write_clean_markdown

RESULTS = REPO / "results"
REPORT = REPO / "reports" / "followup-11pt-block6-depth-extension.md"
MSSM_PATH = Path("/home/reza/workspace/rev2-mssm-null-scratch/mssm.h5ad")
SEED = 20261120
MIN_CELLS = 20
MIN_DONORS = 6
MSSM_UNITS = [
    "GABAergic neuron",
    "VIP GABAergic cortical interneuron",
    "pvalb GABAergic cortical interneuron",
    "sst GABAergic cortical interneuron",
]
EXPECTED_FAMILIES = {
    ("gao", "mean_log_umi_ge1500"): 13,
    ("yazar", "mean_log_umi_ge1500"): 18,
    ("gao", "mean_log_umi_all"): 13,
    ("yazar", "mean_log_umi_all"): 19,
}
OUTCOME_CELLS = {
    "mean_log_umi_ge1500": "n_cells_ge1500",
    "mean_log_umi_all": "n_cells_all",
}


def fit_ols(frame: pd.DataFrame, outcome: str, adjusted: bool, batch: bool) -> dict:
    work = frame.dropna(subset=[outcome, "age_per_decade"]).copy()
    n_cells = int(work["n_cells"].sum())
    base = {
        "n_donors": int(len(work)),
        "n_cells": n_cells,
        "age_min": float(work["age_years"].min()) if len(work) else np.nan,
        "age_max": float(work["age_years"].max()) if len(work) else np.nan,
        "age_median": float(work["age_years"].median()) if len(work) else np.nan,
    }
    if len(work) < MIN_DONORS:
        return {**base, "status": "UNTESTABLE", "status_reason": "fewer_than_6_donors"}
    if work["age_per_decade"].nunique() < 2:
        return {**base, "status": "UNTESTABLE", "status_reason": "age_has_fewer_than_2_values"}
    y = work[outcome].to_numpy(dtype=float)
    outcome_sd = float(np.std(y, ddof=1))
    if not np.isfinite(outcome_sd) or outcome_sd <= 0:
        return {**base, "status": "UNTESTABLE", "status_reason": "outcome_has_zero_variance"}

    design = pd.DataFrame(
        {"intercept": np.ones(len(work)), "age_per_decade": work["age_per_decade"].to_numpy(dtype=float)},
        index=work.index,
    )
    active = ["intercept", "age_per_decade"]
    omitted: list[str] = []

    def add_column(name: str, values: np.ndarray) -> None:
        candidate = design.copy()
        candidate[name] = values
        if np.linalg.matrix_rank(candidate.to_numpy(dtype=float)) > np.linalg.matrix_rank(design.to_numpy(dtype=float)):
            design[name] = values
            active.append(name)
        else:
            omitted.append(f"{name}(collinear)")

    if adjusted:
        sex = work["sex_male"].to_numpy(dtype=float)
        if np.isfinite(sex).all() and np.unique(sex).size > 1:
            add_column("sex_male", sex)
        else:
            omitted.append("sex_male(constant_or_na)")
        if batch:
            batches = pd.get_dummies(
                work["study_or_pool"].astype(str), prefix="batch", drop_first=True, dtype=float
            )
            for name in sorted(batches.columns):
                add_column(name, batches[name].to_numpy(dtype=float))

    residual_df = len(work) - np.linalg.matrix_rank(design.to_numpy(dtype=float))
    if residual_df < 1:
        return {
            **base,
            "status": "UNTESTABLE",
            "status_reason": "insufficient_residual_df",
            "covariates_active": ";".join(active),
            "covariates_omitted": ";".join(omitted) if omitted else "none",
        }

    fit = sm.OLS(y, design.astype(float)).fit()
    slope = float(fit.params["age_per_decade"])
    se = float(fit.bse["age_per_decade"])
    df_resid = int(fit.df_resid)
    t975 = float(stats.t.ppf(0.975, df_resid))
    t80 = float(stats.t.ppf(0.80, df_resid))
    mde80 = (t975 + t80) * se
    return {
        **base,
        "status": "FIT",
        "status_reason": "none",
        "slope_per_decade": slope,
        "se": se,
        "ci95_low": slope - t975 * se,
        "ci95_high": slope + t975 * se,
        "pvalue": float(fit.pvalues["age_per_decade"]),
        "residual_df": df_resid,
        "outcome_sd": outcome_sd,
        "slope_sd_per_decade": slope / outcome_sd,
        "mde80_per_decade": mde80,
        "mde80_sd_per_decade": mde80 / outcome_sd,
        "covariates_active": ";".join(active),
        "covariates_omitted": ";".join(omitted) if omitted else "none",
    }


def empty_fit(n_donors: int = 0) -> dict:
    return {
        "status": "UNTESTABLE",
        "status_reason": "fewer_than_6_donors",
        "n_donors": int(n_donors),
        "n_cells": 0,
        "age_min": np.nan,
        "age_max": np.nan,
        "age_median": np.nan,
    }


def run_part_a() -> pd.DataFrame:
    rows: list[dict] = []
    candidate_counts: dict[tuple[str, str], int] = {}
    for dataset in ("gao", "yazar"):
        source = RESULTS / f"followup-between-subject-{dataset}-donors.csv"
        donors = pd.read_csv(source)
        pooled = donors["is_pooled"].astype(str).str.lower().eq("true")
        donors = donors.loc[~pooled].copy()
        if donors.duplicated(["group", "donor_id"]).any():
            raise ValueError(f"Duplicate donor cell-type rows in {source}")
        cell_types = sorted(donors["group"].astype(str).unique())
        for outcome, cell_count in OUTCOME_CELLS.items():
            eligible = donors.loc[
                (donors[cell_count] >= MIN_CELLS)
                & donors[outcome].notna()
                & donors["age_per_decade"].notna()
                & donors["sex_male"].notna()
                & donors["study_or_pool"].notna()
            ].copy()
            counts = eligible.groupby("group")["donor_id"].nunique()
            candidates = {str(cell_type) for cell_type, n in counts.items() if n >= MIN_DONORS}
            candidate_counts[(dataset, outcome)] = len(candidates)
            for cell_type in cell_types:
                group = eligible.loc[eligible["group"].astype(str) == cell_type].copy()
                group["n_cells"] = group[cell_count]
                if len(group) < MIN_DONORS:
                    fit = empty_fit(len(group))
                else:
                    fit = fit_ols(group, outcome, adjusted=True, batch=True)
                rows.append({
                    "dataset": dataset,
                    "cell_type": cell_type,
                    "outcome": outcome,
                    "model": "primary" if outcome == "mean_log_umi_ge1500" else "sensitivity_adjusted",
                    "in_family": cell_type in candidates,
                    "family": f"{dataset}_{outcome}_adjusted",
                    **fit,
                })
                if len(group) < MIN_DONORS:
                    unadjusted = empty_fit(len(group))
                else:
                    unadjusted = fit_ols(group, outcome, adjusted=False, batch=False)
                rows.append({
                    "dataset": dataset,
                    "cell_type": cell_type,
                    "outcome": outcome,
                    "model": "unadjusted" if outcome == "mean_log_umi_ge1500" else "sensitivity_unadjusted",
                    "in_family": False,
                    "family": "none",
                    **unadjusted,
                })

    if candidate_counts != EXPECTED_FAMILIES:
        raise ValueError(f"Prefrozen family sizes changed: {candidate_counts} != {EXPECTED_FAMILIES}")
    out = pd.DataFrame(rows)
    out["q_bh"] = np.nan
    out["family_size"] = np.nan
    out["family_n_fitted"] = np.nan
    for (dataset, outcome), expected in EXPECTED_FAMILIES.items():
        mask = (
            out["dataset"].eq(dataset)
            & out["outcome"].eq(outcome)
            & out["model"].isin(["primary", "sensitivity_adjusted"])
            & out["in_family"]
        )
        out.loc[mask, "family_size"] = expected
        fitted = mask & out["status"].eq("FIT") & out["pvalue"].notna()
        out.loc[mask, "family_n_fitted"] = int(fitted.sum())
        if fitted.any():
            out.loc[fitted, "q_bh"] = multipletests(out.loc[fitted, "pvalue"], method="fdr_bh")[1]

    columns = [
        "dataset", "cell_type", "outcome", "model", "status", "status_reason",
        "slope_per_decade", "se", "ci95_low", "ci95_high", "pvalue", "q_bh",
        "n_donors", "n_cells", "mde80_per_decade", "outcome_sd",
        "slope_sd_per_decade", "mde80_sd_per_decade", "residual_df",
        "age_min", "age_median", "age_max", "family", "in_family",
        "family_size", "family_n_fitted", "covariates_active", "covariates_omitted",
    ]
    out = out.reindex(columns=columns).sort_values(
        ["dataset", "outcome", "cell_type", "model"]
    ).reset_index(drop=True)
    write_clean_csv(out, RESULTS / "followup-11pt-block6-depth-per-cell-type.csv")
    return out


def mssm_donor_units() -> dict[str, pd.DataFrame]:
    if not MSSM_PATH.exists():
        raise FileNotFoundError(f"Required MSSM file is missing: {MSSM_PATH}")
    backed = ad.read_h5ad(MSSM_PATH, backed="r")
    obs = backed.obs
    age = _age_from_stage(obs["development_stage"])
    common = _cat_eq(obs["disease"], "normal") & _cat_eq(obs["is_primary_data"], "true")
    common &= pd.to_numeric(obs["n_counts"], errors="coerce").ge(1500)
    common &= age.between(20, 120)
    units: dict[str, pd.DataFrame] = {}
    for unit in MSSM_UNITS:
        mask = common & _cat_eq(obs["cell_type"], unit)
        selected = pd.DataFrame({
            "donor_id": obs.loc[mask, "donor_id"].astype(str).to_numpy(),
            "age_years": age.loc[mask].to_numpy(dtype=float),
            "sex_male": (obs.loc[mask, "sex"].astype(str).str.lower().to_numpy() == "male").astype(float),
            "log_umi": np.log(obs.loc[mask, "n_counts"].to_numpy(dtype=float)),
        })
        grouped = selected.groupby("donor_id", sort=True).agg(
            n_cells=("log_umi", "size"),
            mean_log_umi=("log_umi", "mean"),
            age_years=("age_years", "first"),
            sex_male=("sex_male", "first"),
        ).reset_index()
        grouped = grouped.loc[grouped["n_cells"] >= MIN_CELLS].copy()
        grouped["age_per_decade"] = grouped["age_years"] / 10.0
        units[unit] = grouped
    backed.file.close()
    return units


def run_part_b() -> pd.DataFrame:
    rows: list[dict] = []
    for unit, donors in mssm_donor_units().items():
        for model in ("primary", "unadjusted", "mature_age30_plus"):
            work = donors.loc[donors["age_years"] >= 30].copy() if model == "mature_age30_plus" else donors.copy()
            work = work.rename(columns={"mean_log_umi": "mean_log_umi_full_depth"})
            fit = fit_ols(
                work,
                "mean_log_umi_full_depth",
                adjusted=model != "unadjusted",
                batch=False,
            )
            rows.append({
                "dataset": "MSSM_DLPFC",
                "cell_type": unit,
                "outcome": "mean_log_umi_full_depth",
                "model": model,
                "family": "mssm_neuronal_primary" if model == "primary" else "none",
                **fit,
            })
    out = pd.DataFrame(rows)
    out["q_bh"] = np.nan
    out["family_size"] = np.nan
    primary = out["model"].eq("primary")
    out.loc[primary, "family_size"] = len(MSSM_UNITS)
    fitted = primary & out["status"].eq("FIT") & out["pvalue"].notna()
    if fitted.any():
        out.loc[fitted, "q_bh"] = multipletests(out.loc[fitted, "pvalue"], method="fdr_bh")[1]
    columns = [
        "dataset", "cell_type", "outcome", "model", "status", "status_reason",
        "slope_per_decade", "se", "ci95_low", "ci95_high", "pvalue", "q_bh",
        "n_donors", "n_cells", "mde80_per_decade", "outcome_sd",
        "slope_sd_per_decade", "mde80_sd_per_decade", "residual_df",
        "age_min", "age_median", "age_max", "family", "family_size",
        "covariates_active", "covariates_omitted",
    ]
    out = out.reindex(columns=columns).sort_values(["cell_type", "model"]).reset_index(drop=True)
    write_clean_csv(out, RESULTS / "followup-11pt-block6-mssm-depth.csv")
    return out


def fmt(value: float, digits: int = 4) -> str:
    if pd.isna(value):
        return "NA"
    if value != 0 and abs(value) < 0.001:
        return f"{value:.2e}"
    return f"{value:.{digits}f}"


def markdown_table(headers: list[str], rows: list[list[str]]) -> str:
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join(["---"] * len(headers)) + "|"]
    lines.extend("| " + " | ".join(row) + " |" for row in rows)
    return "\n".join(lines)


def verdict_a(frame: pd.DataFrame, dataset: str) -> tuple[str, str]:
    primary = frame.loc[
        frame["dataset"].eq(dataset) & frame["model"].eq("primary") & frame["in_family"]
    ]
    hits = primary.loc[
        primary["status"].eq("FIT") & primary["slope_per_decade"].lt(0) & primary["q_bh"].lt(0.05)
    ]
    if len(hits):
        names = ", ".join(hits["cell_type"])
        return "PROVES", f"Negative batch-adjusted slopes pass BH correction in {names}."
    untestable = primary["status"].ne("FIT")
    fitted = primary.loc[primary["status"].eq("FIT")]
    adequately_powered = fitted["mde80_sd_per_decade"].le(0.50)
    if len(fitted) and not untestable.any() and adequately_powered.all():
        return "DOES NOT PROVE", "No adjusted primary slope passes BH for a decline, and every candidate has MDE80 at most 0.50 SD per decade."
    return "INCONCLUSIVE", "No adjusted primary slope passes BH for a decline, but at least one candidate is UNTESTABLE or has MDE80 above 0.50 SD per decade."


def verdict_b(frame: pd.DataFrame) -> tuple[str, str]:
    primary = frame.loc[frame["model"].eq("primary")]
    hits = primary.loc[
        primary["status"].eq("FIT") & primary["slope_per_decade"].lt(0) & primary["q_bh"].lt(0.05)
    ]
    if len(hits):
        return "PROVES", f"The direct formal test supports declining depth in {', '.join(hits['cell_type'])}."
    if (
        len(primary) == len(MSSM_UNITS)
        and primary["status"].eq("FIT").all()
        and primary["mde80_sd_per_decade"].le(0.50).all()
    ):
        return "DOES NOT PROVE", "The direct formal test does not support a decline, and all four units have MDE80 at most 0.50 SD per decade."
    return "INCONCLUSIVE", "The direct formal test does not support a decline, but at least one unit is UNTESTABLE or has MDE80 above 0.50 SD per decade."


def append_report(part_a: pd.DataFrame, part_b: pd.DataFrame) -> None:
    spec = REPORT.read_text(encoding="utf-8").split("## Results", 1)[0].rstrip()
    spec = spec.replace(
        "Pre-analysis specification frozen before model fitting. Results and verdicts are appended by the analysis script after the fits.",
        "Complete. The specification and correction families below were frozen before model fitting.",
    )
    sections = [spec, "", "## Results", "", "All slopes and MDE80 values are natural-log UMI per decade. FIT and UNTESTABLE describe model evaluability, not scientific verdicts."]

    primary = part_a.loc[part_a["model"].eq("primary") & part_a["in_family"]].copy()
    table_rows = []
    for row in primary.itertuples():
        table_rows.append([
            row.dataset, row.cell_type, str(row.n_donors), fmt(row.slope_per_decade),
            f"{fmt(row.ci95_low)} to {fmt(row.ci95_high)}", fmt(row.pvalue), fmt(row.q_bh),
            fmt(row.mde80_per_decade), fmt(row.mde80_sd_per_decade), row.status,
        ])
    sections.extend([
        "", "### Part A primary fits", "",
        markdown_table(
            ["Dataset", "Cell type", "N", "Slope", "95% CI", "p", "BH q", "MDE80", "MDE80 SD", "Status"],
            table_rows,
        ),
    ])

    below_floor = part_a.loc[
        part_a["model"].eq("primary") & ~part_a["in_family"] & part_a["status"].eq("UNTESTABLE")
    ]
    if len(below_floor):
        labels = ", ".join(f"{r.dataset}: {r.cell_type} (N={r.n_donors})" for r in below_floor.itertuples())
        sections.extend(["", "Cell types below the six-donor primary floor are UNTESTABLE: " + labels + "."])

    sensitivity_rows = []
    for (dataset, model), group in part_a.loc[part_a["model"].ne("primary")].groupby(["dataset", "model"]):
        fitted = group.loc[group["status"].eq("FIT")]
        if model == "sensitivity_adjusted":
            negative = int((fitted["slope_per_decade"].lt(0) & fitted["q_bh"].lt(0.05)).sum())
            positive = int((fitted["slope_per_decade"].gt(0) & fitted["q_bh"].lt(0.05)).sum())
            threshold = "BH q < 0.05"
        else:
            negative = int((fitted["slope_per_decade"].lt(0) & fitted["pvalue"].lt(0.05)).sum())
            positive = int((fitted["slope_per_decade"].gt(0) & fitted["pvalue"].lt(0.05)).sum())
            threshold = "nominal p < 0.05"
        sensitivity_rows.append([
            dataset, model, str(len(fitted)), str(int(group["status"].ne("FIT").sum())),
            str(negative), str(positive), threshold,
        ])
    sections.extend([
        "", "### Part A sensitivity summary", "",
        markdown_table(
            ["Dataset", "Model", "Fitted", "UNTESTABLE", "Negative hits", "Positive hits", "Threshold"],
            sensitivity_rows,
        ),
        "", "The adjusted `mean_log_umi_all` fits use their separate pre-stated BH families. Unadjusted fits are age-only and are labeled with nominal p-values. Full per-unit sensitivity estimates, CIs, and MDE80 values are in the Part A CSV.",
    ])

    part_b_rows = []
    for row in part_b.itertuples():
        part_b_rows.append([
            row.cell_type, row.model, str(row.n_donors), fmt(row.slope_per_decade),
            f"{fmt(row.ci95_low)} to {fmt(row.ci95_high)}", fmt(row.pvalue), fmt(row.q_bh),
            fmt(row.mde80_per_decade), fmt(row.mde80_sd_per_decade), row.status,
        ])
    sections.extend([
        "", "### Part B MSSM formal test", "",
        markdown_table(
            ["Cell type", "Model", "N", "Slope", "95% CI", "p", "BH q", "MDE80", "MDE80 SD", "Status"],
            part_b_rows,
        ),
        "", "The primary q-values are corrected across exactly four neuronal units. The age-only and mature 30 years and older fits are sensitivities and have no q-values.",
    ])

    sections.extend(["", "## Verdicts", ""])
    for dataset in ("gao", "yazar"):
        verdict, reason = verdict_a(part_a, dataset)
        sections.append(f"- Part A, {dataset}: {verdict}. {reason}")
    verdict, reason = verdict_b(part_b)
    sections.extend([
        f"- Part B, MSSM: {verdict}. {reason}",
        "",
        f"Does the MSSM neuron depth decline survive a direct formal test? {verdict}. {reason}",
        "",
        "The Part A verdicts concern per-cell-type depth slopes and do not repeat the pooled November 17 tests. The Part B verdict concerns the four pre-stated neuronal units only.",
        "", "## Output files", "",
        "- `results/followup-11pt-block6-depth-per-cell-type.csv`",
        "- `results/followup-11pt-block6-mssm-depth.csv`",
        "- `scripts/run_followup_11pt_block6_depth_extension.py`",
        "",
    ])
    write_clean_markdown("\n".join(sections), REPORT)


def main() -> None:
    np.random.seed(SEED)
    part_a = run_part_a()
    part_b = run_part_b()
    append_report(part_a, part_b)
    for dataset in ("gao", "yazar"):
        verdict, reason = verdict_a(part_a, dataset)
        print(f"Part A {dataset}: {verdict}. {reason}", flush=True)
    verdict, reason = verdict_b(part_b)
    print(f"Part B MSSM: {verdict}. {reason}", flush=True)
    print(f"Wrote {len(part_a)} Part A rows and {len(part_b)} Part B rows.", flush=True)


if __name__ == "__main__":
    main()
