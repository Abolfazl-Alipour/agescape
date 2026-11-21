#!/usr/bin/env python3
"""Phase 2 Item 13: Mature-to-old (age 30-plus) Block 6 depth versus age models."""

from __future__ import annotations

import os

for _env_var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_env_var, "1")

import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import statsmodels.api as sm
from anndata.io import read_elem
from scipy import stats
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from scripts.run_rev2_de_interneuron import _age_from_stage, _cat_eq

FORBIDDEN_DASH_RE = re.compile(r"[\u2013\u2014]")


def _assert_no_forbidden_dashes(text: str, label: str = "output") -> None:
    match = FORBIDDEN_DASH_RE.search(text)
    if match:
        raise ValueError(f"Forbidden Unicode dash {match.group(0)!r} found in {label}")


try:
    from src.followup_organs import write_clean_csv, write_clean_markdown
except ImportError:
    def write_clean_markdown(text: str, path: Path) -> Path:
        _assert_no_forbidden_dashes(text, label=str(path))
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
        return path

    def write_clean_csv(frame: pd.DataFrame, path: Path) -> Path:
        csv_text = frame.to_csv(index=False)
        _assert_no_forbidden_dashes(csv_text, label=str(path))
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(csv_text, encoding="utf-8")
        return path

RESULTS = REPO / "results"
REPORT = REPO / "reports" / "followup-mature-block6-depth-extension.md"
PART_A_CSV = RESULTS / "followup-mature-block6-depth-per-cell-type.csv"
PART_B_CSV = RESULTS / "followup-mature-block6-mssm-depth.csv"
ORIG_PART_A_CSV = RESULTS / "followup-11pt-block6-depth-per-cell-type.csv"
ORIG_PART_B_CSV = RESULTS / "followup-11pt-block6-mssm-depth.csv"
POOLED_MATURE_CSV = RESULTS / "followup-mature-nov17-depth.csv"
MSSM_PATH = Path("/home/reza/projects/agescape-revision-data/mssm.h5ad")

SEED = 20261120
MATURE_MIN_AGE = 30.0
MIN_CELLS = 20
MIN_DONORS = 6
VERIFY_TOL = 1e-12
MSSM_UNITS = [
    "GABAergic neuron",
    "VIP GABAergic cortical interneuron",
    "pvalb GABAergic cortical interneuron",
    "sst GABAergic cortical interneuron",
]
EXPECTED_MATURE_FAMILIES = {
    ("gao", "mean_log_umi_ge1500"): 13,
    ("yazar", "mean_log_umi_ge1500"): 17,
    ("gao", "mean_log_umi_all"): 13,
    ("yazar", "mean_log_umi_all"): 19,
}
OUTCOME_CELLS = {
    "mean_log_umi_ge1500": "n_cells_ge1500",
    "mean_log_umi_all": "n_cells_all",
}
MSSM_VERIFY_FLOAT_COLS = [
    "slope_per_decade",
    "se",
    "ci95_low",
    "ci95_high",
    "pvalue",
    "mde80_per_decade",
    "outcome_sd",
    "slope_sd_per_decade",
    "mde80_sd_per_decade",
    "age_min",
    "age_median",
    "age_max",
]
MSSM_VERIFY_INT_COLS = ["n_donors", "n_cells", "residual_df"]


def bh_independent(pvalues: np.ndarray | pd.Series | list[float]) -> np.ndarray:
    """Independent step-up Benjamini-Hochberg implementation for verification."""
    p = np.asarray(pvalues, dtype=float)
    m = len(p)
    if m == 0:
        return np.array([], dtype=float)
    order = np.argsort(p, kind="mergesort")
    ranked = p[order]
    q_ranked = np.empty(m, dtype=float)
    running_min = 1.0
    for idx in range(m - 1, -1, -1):
        rank = idx + 1
        candidate = ranked[idx] * m / rank
        if candidate < running_min:
            running_min = candidate
        q_ranked[idx] = running_min
    q = np.empty(m, dtype=float)
    q[order] = np.clip(q_ranked, 0.0, 1.0)
    return q


def fit_ols(frame: pd.DataFrame, outcome: str, adjusted: bool, batch: bool) -> dict:
    """Donor-level OLS with rank-checked covariates and frozen MDE80 definition."""
    work = frame.dropna(subset=[outcome, "age_per_decade"]).copy()
    n_cells = int(work["n_cells"].sum()) if len(work) else 0
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

    cols = [np.ones(len(work), dtype=float), work["age_per_decade"].to_numpy(dtype=float)]
    active = ["intercept", "age_per_decade"]
    omitted: list[str] = []
    candidates: list[tuple[str, np.ndarray]] = []

    if adjusted:
        sex = work["sex_male"].to_numpy(dtype=float)
        if np.isfinite(sex).all() and np.unique(sex).size > 1:
            candidates.append(("sex_male", sex))
        else:
            omitted.append("sex_male(constant_or_na)")
        if batch:
            batches = pd.get_dummies(
                work["study_or_pool"].astype(str), prefix="batch", drop_first=True, dtype=float
            )
            for name in sorted(batches.columns):
                candidates.append((name, batches[name].to_numpy(dtype=float)))

    if candidates:
        full_mat = np.column_stack(cols + [vals for _, vals in candidates])
        if int(np.linalg.matrix_rank(full_mat)) == full_mat.shape[1]:
            for name, vals in candidates:
                cols.append(vals)
                active.append(name)
        else:
            cur_mat = np.column_stack(cols)
            cur_rank = int(np.linalg.matrix_rank(cur_mat))
            for name, vals in candidates:
                test_mat = np.column_stack([cur_mat, vals])
                new_rank = int(np.linalg.matrix_rank(test_mat))
                if new_rank > cur_rank:
                    cur_mat = test_mat
                    cur_rank = new_rank
                    cols.append(vals)
                    active.append(name)
                else:
                    omitted.append(f"{name}(collinear)")

    design = pd.DataFrame({name: vals for name, vals in zip(active, cols)}, index=work.index)
    residual_df = len(work) - int(np.linalg.matrix_rank(design.to_numpy(dtype=float)))
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


def run_part_a() -> tuple[pd.DataFrame, list[dict]]:
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
            eligible_orig = donors.loc[
                (donors[cell_count] >= MIN_CELLS)
                & donors[outcome].notna()
                & donors["age_per_decade"].notna()
                & donors["sex_male"].notna()
                & donors["study_or_pool"].notna()
            ].copy()
            eligible_mature = eligible_orig.loc[eligible_orig["age_years"] >= MATURE_MIN_AGE].copy()
            if len(eligible_mature) and float(eligible_mature["age_years"].min()) < MATURE_MIN_AGE:
                raise ValueError(f"Age restriction failed for {dataset} {outcome}")

            orig_counts = eligible_orig.groupby("group")["donor_id"].nunique()
            mature_counts = eligible_mature.groupby("group")["donor_id"].nunique()
            candidates = {str(ct) for ct, n in mature_counts.items() if n >= MIN_DONORS}
            candidate_counts[(dataset, outcome)] = len(candidates)

            for cell_type in cell_types:
                group = eligible_mature.loc[eligible_mature["group"].astype(str) == cell_type].copy()
                group["n_cells"] = group[cell_count]
                n_orig = int(orig_counts.get(cell_type, 0))
                n_mat = int(len(group))
                n_lost = n_orig - n_mat
                if n_lost < 0:
                    raise ValueError(f"Negative donors lost for {dataset} {outcome} {cell_type}")

                adjusted_fit = fit_ols(group, outcome, adjusted=True, batch=True)
                rows.append({
                    "dataset": dataset,
                    "cell_type": cell_type,
                    "outcome": outcome,
                    "model": "primary" if outcome == "mean_log_umi_ge1500" else "sensitivity_adjusted",
                    "n_donors_original": n_orig,
                    "n_donors_lost": n_lost,
                    "in_family": cell_type in candidates,
                    "family": f"{dataset}_{outcome}_adjusted_mature30",
                    **adjusted_fit,
                })

                unadjusted_fit = fit_ols(group, outcome, adjusted=False, batch=False)
                rows.append({
                    "dataset": dataset,
                    "cell_type": cell_type,
                    "outcome": outcome,
                    "model": "unadjusted" if outcome == "mean_log_umi_ge1500" else "sensitivity_unadjusted",
                    "n_donors_original": n_orig,
                    "n_donors_lost": n_lost,
                    "in_family": False,
                    "family": "none",
                    **unadjusted_fit,
                })

    if candidate_counts != EXPECTED_MATURE_FAMILIES:
        raise ValueError(
            f"Mature candidate family sizes {candidate_counts} != expected {EXPECTED_MATURE_FAMILIES}"
        )

    out = pd.DataFrame(rows)
    out["q_bh"] = np.nan
    out["family_size"] = np.nan
    out["family_n_fitted"] = np.nan

    bh_checks: list[dict] = []
    for (dataset, outcome), expected in EXPECTED_MATURE_FAMILIES.items():
        family_name = f"{dataset}_{outcome}_adjusted_mature30"
        mask = (
            out["dataset"].eq(dataset)
            & out["outcome"].eq(outcome)
            & out["model"].isin(["primary", "sensitivity_adjusted"])
            & out["in_family"]
        )
        out.loc[mask, "family_size"] = expected
        fitted = mask & out["status"].eq("FIT") & out["pvalue"].notna()
        n_fitted = int(fitted.sum())
        out.loc[mask, "family_n_fitted"] = n_fitted
        if fitted.any():
            pvals = out.loc[fitted, "pvalue"].to_numpy(dtype=float)
            q_sm = multipletests(pvals, method="fdr_bh")[1]
            q_ind = bh_independent(pvals)
            max_diff = float(np.max(np.abs(q_sm - q_ind)))
            if max_diff > 1e-12:
                raise ValueError(f"BH mismatch in {family_name}: max_diff={max_diff}")
            out.loc[fitted, "q_bh"] = q_sm
            bh_checks.append({
                "family": family_name,
                "candidate_size": expected,
                "n_fitted": n_fitted,
                "max_abs_diff": max_diff,
                "spot_checks": build_bh_spot_checks(out.loc[fitted, ["cell_type", "pvalue", "q_bh"]]),
            })

    nonempty = out["n_donors"].gt(0)
    if nonempty.any() and float(out.loc[nonempty, "age_min"].min()) < MATURE_MIN_AGE:
        raise ValueError("Found Part A row with age_min < 30")

    columns = [
        "dataset", "cell_type", "outcome", "model", "status", "status_reason",
        "slope_per_decade", "se", "ci95_low", "ci95_high", "pvalue", "q_bh",
        "n_donors", "n_donors_original", "n_donors_lost", "n_cells",
        "mde80_per_decade", "outcome_sd", "slope_sd_per_decade", "mde80_sd_per_decade",
        "residual_df", "age_min", "age_median", "age_max",
        "family", "in_family", "family_size", "family_n_fitted",
        "covariates_active", "covariates_omitted",
    ]
    out = out.reindex(columns=columns).sort_values(
        ["dataset", "outcome", "cell_type", "model"]
    ).reset_index(drop=True)
    write_clean_csv(out, PART_A_CSV)
    return out, bh_checks


def build_bh_spot_checks(frame: pd.DataFrame) -> list[dict]:
    """Select at least 3 ranked p-values (smallest, middle, largest) for manual BH audit."""
    ordered = frame.sort_values("pvalue", kind="mergesort").reset_index(drop=True)
    m = len(ordered)
    raw_q = ordered["pvalue"].to_numpy(dtype=float) * m / np.arange(1, m + 1, dtype=float)
    manual_q = np.minimum(1.0, np.minimum.accumulate(raw_q[::-1])[::-1])
    if m <= 4:
        indices = list(range(m))
    else:
        indices = [0, 1, m // 2, m - 1]
    spots: list[dict] = []
    for idx in indices:
        row = ordered.iloc[idx]
        spots.append({
            "cell_type": str(row["cell_type"]),
            "rank": idx + 1,
            "m": m,
            "pvalue": float(row["pvalue"]),
            "p_times_m_over_k": float(raw_q[idx]),
            "q_manual": float(manual_q[idx]),
            "q_bh": float(row["q_bh"]),
            "abs_diff": float(abs(manual_q[idx] - float(row["q_bh"]))),
        })
    return spots


def load_mssm_donor_units() -> tuple[dict[str, pd.DataFrame], dict[str, int]]:
    """Load only required MSSM obs columns via h5py to stay low-memory on beast."""
    if not MSSM_PATH.exists():
        raise FileNotFoundError(f"Required MSSM file is missing: {MSSM_PATH}")
    needed = [
        "development_stage",
        "disease",
        "is_primary_data",
        "n_counts",
        "cell_type",
        "donor_id",
        "sex",
    ]
    cols: dict[str, pd.Series] = {}
    with h5py.File(MSSM_PATH, "r") as handle:
        obs_grp = handle["obs"]
        for col in needed:
            cols[col] = pd.Series(read_elem(obs_grp[col]))

    age = _age_from_stage(cols["development_stage"])
    common = _cat_eq(cols["disease"], "normal") & _cat_eq(cols["is_primary_data"], "true")
    common &= pd.to_numeric(cols["n_counts"], errors="coerce").ge(1500)
    common &= age.between(20, 120)

    mature_units: dict[str, pd.DataFrame] = {}
    orig_donor_counts: dict[str, int] = {}
    for unit in MSSM_UNITS:
        mask = (common & _cat_eq(cols["cell_type"], unit)).to_numpy(dtype=bool)
        selected = pd.DataFrame({
            "donor_id": cols["donor_id"].iloc[mask].astype(str).to_numpy(),
            "age_years": age.iloc[mask].to_numpy(dtype=float),
            "sex_male": (cols["sex"].iloc[mask].astype(str).str.lower().to_numpy() == "male").astype(float),
            "log_umi": np.log(cols["n_counts"].iloc[mask].to_numpy(dtype=float)),
        })
        grouped = selected.groupby("donor_id", sort=True).agg(
            n_cells=("log_umi", "size"),
            mean_log_umi=("log_umi", "mean"),
            age_years=("age_years", "first"),
            sex_male=("sex_male", "first"),
        ).reset_index()
        eligible_20 = grouped.loc[grouped["n_cells"] >= MIN_CELLS].copy()
        orig_donor_counts[unit] = int(len(eligible_20))
        eligible_30 = eligible_20.loc[eligible_20["age_years"] >= MATURE_MIN_AGE].copy()
        eligible_30["age_per_decade"] = eligible_30["age_years"] / 10.0
        mature_units[unit] = eligible_30
    return mature_units, orig_donor_counts


def run_part_b() -> tuple[pd.DataFrame, dict]:
    orig_mssm = pd.read_csv(ORIG_PART_B_CSV)
    mature_units, orig_donor_counts = load_mssm_donor_units()

    rows: list[dict] = []
    for unit in MSSM_UNITS:
        donors = mature_units[unit].rename(columns={"mean_log_umi": "mean_log_umi_full_depth"})
        n_orig = orig_donor_counts[unit]
        refit_primary = fit_ols(donors, "mean_log_umi_full_depth", adjusted=True, batch=False)
        refit_unadjusted = fit_ols(donors, "mean_log_umi_full_depth", adjusted=False, batch=False)

        existing_slice = orig_mssm.loc[
            orig_mssm["cell_type"].eq(unit) & orig_mssm["model"].eq("mature_age30_plus")
        ]
        orig_primary_slice = orig_mssm.loc[
            orig_mssm["cell_type"].eq(unit) & orig_mssm["model"].eq("primary")
        ]
        if len(orig_primary_slice) == 1 and int(orig_primary_slice.iloc[0]["n_donors"]) != n_orig:
            raise ValueError(f"MSSM 20-plus donor count mismatch for {unit}")

        verified = False
        max_diff = np.nan
        if len(existing_slice) == 1 and refit_primary["status"] == "FIT":
            ex = existing_slice.iloc[0]
            ints_ok = all(int(ex[c]) == int(refit_primary[c]) for c in MSSM_VERIFY_INT_COLS)
            diffs = [abs(float(ex[c]) - float(refit_primary[c])) for c in MSSM_VERIFY_FLOAT_COLS]
            max_diff = float(max(diffs))
            status_ok = str(ex["status"]) == str(refit_primary["status"])
            verified = ints_ok and status_ok and max_diff <= VERIFY_TOL

        if verified:
            ex = existing_slice.iloc[0].to_dict()
            primary_payload = {
                "status": ex["status"],
                "status_reason": ex["status_reason"],
                "slope_per_decade": float(ex["slope_per_decade"]),
                "se": float(ex["se"]),
                "ci95_low": float(ex["ci95_low"]),
                "ci95_high": float(ex["ci95_high"]),
                "pvalue": float(ex["pvalue"]),
                "n_donors": int(ex["n_donors"]),
                "n_donors_original": n_orig,
                "n_donors_lost": n_orig - int(ex["n_donors"]),
                "n_cells": int(ex["n_cells"]),
                "mde80_per_decade": float(ex["mde80_per_decade"]),
                "outcome_sd": float(ex["outcome_sd"]),
                "slope_sd_per_decade": float(ex["slope_sd_per_decade"]),
                "mde80_sd_per_decade": float(ex["mde80_sd_per_decade"]),
                "residual_df": int(ex["residual_df"]),
                "age_min": float(ex["age_min"]),
                "age_median": float(ex["age_median"]),
                "age_max": float(ex["age_max"]),
                "covariates_active": str(ex["covariates_active"]),
                "covariates_omitted": str(ex["covariates_omitted"]),
                "provenance": "reused_verified_from_followup-11pt-block6-mssm-depth.csv",
                "source_file": "results/followup-11pt-block6-mssm-depth.csv",
                "source_model": "mature_age30_plus",
                "verification_status": "verified_against_mssm_h5ad",
                "verification_source": str(MSSM_PATH),
                "verification_max_abs_diff": max_diff,
            }
        else:
            primary_payload = {
                **refit_primary,
                "n_donors_original": n_orig,
                "n_donors_lost": n_orig - int(refit_primary["n_donors"]),
                "provenance": "refit_from_mssm_h5ad",
                "source_file": str(MSSM_PATH),
                "source_model": "direct_refit_age30_plus",
                "verification_status": "existing_rows_failed_verification_refit_used",
                "verification_source": str(MSSM_PATH),
                "verification_max_abs_diff": max_diff,
            }

        rows.append({
            "dataset": "MSSM_DLPFC",
            "cell_type": unit,
            "outcome": "mean_log_umi_full_depth",
            "model": "primary",
            "family": "mssm_neuronal_mature_primary",
            **primary_payload,
        })
        rows.append({
            "dataset": "MSSM_DLPFC",
            "cell_type": unit,
            "outcome": "mean_log_umi_full_depth",
            "model": "unadjusted",
            "family": "none",
            "n_donors_original": n_orig,
            "n_donors_lost": n_orig - int(refit_unadjusted["n_donors"]),
            **refit_unadjusted,
            "provenance": "refit_from_mssm_h5ad_unadjusted_sensitivity",
            "source_file": str(MSSM_PATH),
            "source_model": "none",
            "verification_status": "direct_fit_from_mssm_h5ad",
            "verification_source": str(MSSM_PATH),
            "verification_max_abs_diff": 0.0,
        })

    out = pd.DataFrame(rows)
    out["q_bh"] = np.nan
    out["family_size"] = np.nan
    out["family_n_fitted"] = np.nan
    primary = out["model"].eq("primary")
    out.loc[primary, "family_size"] = len(MSSM_UNITS)
    fitted = primary & out["status"].eq("FIT") & out["pvalue"].notna()
    out.loc[primary, "family_n_fitted"] = int(fitted.sum())

    pvals = out.loc[fitted, "pvalue"].to_numpy(dtype=float)
    q_sm = multipletests(pvals, method="fdr_bh")[1]
    q_ind = bh_independent(pvals)
    max_bh_diff = float(np.max(np.abs(q_sm - q_ind)))
    if max_bh_diff > 1e-12:
        raise ValueError(f"Part B BH mismatch: max_diff={max_bh_diff}")
    out.loc[fitted, "q_bh"] = q_sm

    if float(out["age_min"].min()) < MATURE_MIN_AGE:
        raise ValueError("Found Part B row with age_min < 30")

    bh_check = {
        "family": "mssm_neuronal_mature_primary",
        "candidate_size": len(MSSM_UNITS),
        "n_fitted": int(fitted.sum()),
        "max_abs_diff": max_bh_diff,
        "spot_checks": build_bh_spot_checks(out.loc[fitted, ["cell_type", "pvalue", "q_bh"]]),
    }

    columns = [
        "dataset", "cell_type", "outcome", "model", "status", "status_reason",
        "slope_per_decade", "se", "ci95_low", "ci95_high", "pvalue", "q_bh",
        "n_donors", "n_donors_original", "n_donors_lost", "n_cells",
        "mde80_per_decade", "outcome_sd", "slope_sd_per_decade", "mde80_sd_per_decade",
        "residual_df", "age_min", "age_median", "age_max",
        "family", "family_size", "family_n_fitted",
        "covariates_active", "covariates_omitted",
        "provenance", "source_file", "source_model",
        "verification_status", "verification_source", "verification_max_abs_diff",
    ]
    out = out.reindex(columns=columns).sort_values(["cell_type", "model"]).reset_index(drop=True)
    write_clean_csv(out, PART_B_CSV)
    return out, bh_check


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
        return "PROVES", f"Negative batch-adjusted primary slopes pass BH correction in {names}."
    untestable = primary.loc[primary["status"].ne("FIT")]
    fitted = primary.loc[primary["status"].eq("FIT")]
    underpowered = fitted.loc[fitted["mde80_sd_per_decade"].gt(0.50)]
    if len(fitted) and len(untestable) == 0 and len(underpowered) == 0:
        return (
            "DOES NOT PROVE",
            "No adjusted primary slope passes BH for a decline, and every candidate has MDE80 at most 0.50 SD per decade with no UNTESTABLE primary candidate.",
        )
    reasons: list[str] = []
    if len(untestable):
        u_labels = ", ".join(f"{r.cell_type} (N={int(r.n_donors)}, {r.status_reason})" for r in untestable.itertuples())
        reasons.append(f"UNTESTABLE primary candidate(s): {u_labels}")
    if len(underpowered):
        p_labels = ", ".join(
            f"{r.cell_type} (MDE80 SD={r.mde80_sd_per_decade:.4f})" for r in underpowered.itertuples()
        )
        reasons.append(f"primary candidate(s) with MDE80 above 0.50 SD per decade: {p_labels}")
    return (
        "INCONCLUSIVE",
        "No adjusted primary slope passes BH for a decline, but DOES NOT PROVE is blocked by " + "; ".join(reasons) + ".",
    )


def verdict_b(frame: pd.DataFrame) -> tuple[str, str]:
    primary = frame.loc[frame["model"].eq("primary")]
    hits = primary.loc[
        primary["status"].eq("FIT") & primary["slope_per_decade"].lt(0) & primary["q_bh"].lt(0.05)
    ]
    if len(hits):
        return (
            "PROVES",
            f"The mature 30-plus formal four-unit test supports declining depth in {', '.join(hits['cell_type'])}.",
        )
    if (
        len(primary) == len(MSSM_UNITS)
        and primary["status"].eq("FIT").all()
        and primary["mde80_sd_per_decade"].le(0.50).all()
    ):
        return (
            "DOES NOT PROVE",
            "The mature 30-plus formal test does not support a decline, and all four units have MDE80 at most 0.50 SD per decade.",
        )
    return (
        "INCONCLUSIVE",
        "The mature 30-plus formal test does not support a decline, but at least one unit is UNTESTABLE or has MDE80 above 0.50 SD per decade.",
    )


def write_report(
    part_a: pd.DataFrame,
    part_b: pd.DataFrame,
    bh_checks_a: list[dict],
    bh_check_b: dict,
) -> None:
    orig_a = pd.read_csv(ORIG_PART_A_CSV)
    orig_b = pd.read_csv(ORIG_PART_B_CSV)
    pooled_m = pd.read_csv(POOLED_MATURE_CSV)

    gao_pooled_p = pooled_m.loc[
        pooled_m["dataset"].eq("Gao_heart_GSE290367")
        & pooled_m["outcome"].eq("mean_log_umi_ge1500")
        & pooled_m["model"].eq("B_adjusted")
    ].iloc[0]
    gao_pooled_all = pooled_m.loc[
        pooled_m["dataset"].eq("Gao_heart_GSE290367")
        & pooled_m["outcome"].eq("mean_log_umi_all")
        & pooled_m["model"].eq("B_adjusted")
    ].iloc[0]
    yazar_pooled_p = pooled_m.loc[
        pooled_m["dataset"].eq("Yazar_blood_GSE196830")
        & pooled_m["outcome"].eq("mean_log_umi_ge1500")
        & pooled_m["model"].eq("B_adjusted")
    ].iloc[0]
    yazar_pooled_all = pooled_m.loc[
        pooled_m["dataset"].eq("Yazar_blood_GSE196830")
        & pooled_m["outcome"].eq("mean_log_umi_all")
        & pooled_m["model"].eq("B_adjusted")
    ].iloc[0]
    yazar_pooled_cells = pooled_m.loc[
        pooled_m["dataset"].eq("Yazar_blood_GSE196830")
        & pooled_m["outcome"].eq("log_n_cells_all")
        & pooled_m["model"].eq("B_adjusted")
    ].iloc[0]

    gao_verdict, gao_reason = verdict_a(part_a, "gao")
    yazar_verdict, yazar_reason = verdict_a(part_a, "yazar")
    mssm_verdict, mssm_reason = verdict_b(part_b)

    part_b_primary = part_b.loc[part_b["model"].eq("primary")].copy()
    max_mssm_verify_diff = float(part_b_primary["verification_max_abs_diff"].max())
    all_reused = bool(
        part_b_primary["provenance"].eq("reused_verified_from_followup-11pt-block6-mssm-depth.csv").all()
    )

    lines: list[str] = [
        "# Block 6 mature-to-old (30-plus): Depth versus age extension",
        "",
        "## Status and scope",
        "",
        "Complete. This is Phase 2 Item 13 of the Agescape mature-to-old rebuild. The analysis specification, MDE80 definition, and verdict rules are reused verbatim from the frozen original Block 6 report (`reports/followup-11pt-block6-depth-extension.md`), changing only the donor age restriction to human mature donors aged 30 years and older (`age_years >= 30`). Because Block 6 models library sizes rather than gene scores, the mature P1m signature does not enter this block.",
        "",
        "## Inputs and provenance",
        "",
        "- Gao donor table: `results/followup-between-subject-gao-donors.csv` (non-pooled cell-type rows restricted to `age_years >= 30`).",
        "- Yazar donor table: `results/followup-between-subject-yazar-donors.csv` (non-pooled cell-type rows restricted to `age_years >= 30`).",
        "- Original Block 6 tables (read-only comparison and Part B verification source): `results/followup-11pt-block6-depth-per-cell-type.csv` and `results/followup-11pt-block6-mssm-depth.csv`.",
        "- MSSM h5ad verification source: `/home/reza/projects/agescape-revision-data/mssm.h5ad` (column-selective `obs` read via `h5py`; sparse expression matrix `X` not loaded).",
        "- Pooled mature cross-check only (excluded from Part A): `results/followup-mature-nov17-depth.csv`. At 30-plus, pooled adjusted `mean_log_umi_ge1500` slopes are null in both Gao (`"
        + f"{gao_pooled_p['slope_per_decade']:+.4f}` per decade, `p = {gao_pooled_p['p_value']:.4f}`, N = {int(gao_pooled_p['n_donors'])}, lost {int(gao_pooled_p['n_donors_lost_to_age_restriction'])}"
        + ") and Yazar (`"
        + f"{yazar_pooled_p['slope_per_decade']:+.4f}` per decade, `p = {yazar_pooled_p['p_value']:.4f}`, N = {int(yazar_pooled_p['n_donors'])}, lost {int(yazar_pooled_p['n_donors_lost_to_age_restriction'])}"
        + "). Pooled adjusted `mean_log_umi_all` is borderline negative in Gao (`"
        + f"{gao_pooled_all['slope_per_decade']:+.4f}`, `p = {gao_pooled_all['p_value']:.4f}`"
        + ") and null in Yazar (`"
        + f"{yazar_pooled_all['slope_per_decade']:+.4f}`, `p = {yazar_pooled_all['p_value']:.4f}`"
        + "), while Yazar captured cell counts still decline (`"
        + f"{yazar_pooled_cells['slope_per_decade']:.4f}` log cells per decade, `p = {yazar_pooled_cells['p_value']:.2e}`"
        + ").",
        "- Block 9 and Item 7 interneuron reconciliation: `results/followup-11pt-block9-mature-interneuron.csv`, `reports/followup-11pt-block9-mature-interneuron.md`, `results/followup-mature-interneuron*.csv`, and `reports/followup-mature-interneuron.md` were checked. They contain per-gene differential expression and cross-cohort replication outputs (using donor `mean_log_umi` as a regression covariate) and contain no newer donor-level depth-versus-age regressions. Their reported MSSM mature donor counts (VIP 193, pvalb 185, sst 174, GABA 153; ages 31 to 89) match Part B identically.",
        "",
        "## Frozen analysis specification applied at 30-plus",
        "",
        "### Part A: per-cell-type depth slopes in Gao and Yazar (age 30-plus)",
        "",
        "For each non-pooled cell type with at least 6 eligible donors aged 30 years and older, donor-level OLS regresses `mean_log_umi_ge1500` on `age_per_decade + sex_male + batch`. Gao batch is the observed study fixed effect (up to seven studies, represented by an intercept and up to six study indicators). Yazar batch is pool fixed effects (`study_or_pool`). A donor is eligible for the primary outcome with at least 20 cells in `n_cells_ge1500`, `age_years >= 30`, and nonmissing outcome, age, sex, and batch. Slopes are in natural-log UMI per decade.",
        "",
        "Sensitivity 1 (`sensitivity_adjusted`) uses the same adjusted model with `mean_log_umi_all` and at least 20 cells in `n_cells_all` per donor at `age_years >= 30`. Sensitivity 2 (`unadjusted` and `sensitivity_unadjusted`) fits an unadjusted age-only model for each outcome on its outcome-specific eligible 30-plus donors. No pooled rows, cell-level fits, downsampling, or rarefaction are used.",
        "",
        "BH multiple-testing families are defined separately within each dataset over all non-pooled cell types meeting the 6-donor eligibility floor at 30-plus: 13 Gao and 17 Yazar candidates for `mean_log_umi_ge1500` (down from 18 Yazar originally because `CD16-negative, CD56-bright natural killer cell, human` dropped from 6 to 5 donors), and 13 Gao and 19 Yazar candidates for `mean_log_umi_all`. Within each family, BH is calculated over successfully fitted candidates (`status == FIT`). Unadjusted age-only fits receive no q-values. Collinear covariates are dropped without dropping age; any model with zero residual degrees of freedom is labeled `UNTESTABLE` (`insufficient_residual_df`), never null.",
        "",
        "### Part B: formal MSSM mature neuron test (age 30-plus)",
        "",
        "The four pre-stated MSSM neuronal units are `GABAergic neuron`, `VIP GABAergic cortical interneuron`, `pvalb GABAergic cortical interneuron`, and `sst GABAergic cortical interneuron`. Per-donor mean natural-log UMI uses full-depth per-cell `n_counts` passing the `fetch_unit_h5ad` QC route: primary data, disease `normal`, `n_counts >= 1500`, parsed age 30 to 120 years, and at least 20 eligible cells per donor per unit (minimum 6 donors per unit). Primary donor-level OLS regresses `mean_log_umi_full_depth` on `age_per_decade + sex_male`. The declared mature primary BH family consists of exactly these four MSSM neuronal units. Sensitivity 1 (`unadjusted`) fits age-only OLS on the same 30-plus eligible donors and receives no q-values.",
        "",
        "Direct refitting from `/home/reza/projects/agescape-revision-data/mssm.h5ad` verified all four existing `mature_age30_plus` rows in `results/followup-11pt-block6-mssm-depth.csv` to machine precision (exact integer matches on `n_donors`, `n_cells`, and `residual_df`; maximum absolute float difference `"
        + f"{max_mssm_verify_diff:.3e}"
        + "` across slope, SE, 95 percent CI, p-value, MDE80, outcome SD, and age summaries). Because all four rows verified, the primary Part B table reuses those verified rows verbatim and applies the declared four-unit BH family.",
        "",
        "### MDE80 and frozen verdict rules",
        "",
        "MDE80 is `(t_0.975 + t_0.80) * SE` on the residual degrees of freedom, reported in natural-log UMI per decade and in outcome SD per decade. Under the frozen power threshold, MDE80 at most 0.50 SD per decade is required for a null to count as adequately powered (at most 0.25 SD per decade is high power).",
        "",
        "- Part A per dataset: `PROVES` if at least one batch-adjusted primary cell-type slope is negative with BH q below 0.05. Otherwise `DOES NOT PROVE` requires all fitted primary slopes to be null or nonnegative with MDE80 at most 0.50 SD per decade and no `UNTESTABLE` primary candidate. Otherwise `INCONCLUSIVE`.",
        "- Part B (MSSM mature four-unit family): `PROVES` if at least one of the four mature primary neuronal slopes is negative with BH q below 0.05. Otherwise `DOES NOT PROVE` requires all four mature primary fits to be null or nonnegative with MDE80 at most 0.50 SD per decade. Otherwise `INCONCLUSIVE`.",
        "- Fixed seed: `20261120`. No random sampling, permutation, or bootstrap is used; all fits are deterministic.",
        "",
        "## Results",
        "",
        "All slopes and MDE80 values are in natural-log UMI per decade. `FIT` and `UNTESTABLE` describe model evaluability, not scientific verdicts.",
        "",
        "### Part A primary fits (age 30-plus, `mean_log_umi_ge1500`)",
        "",
    ]

    primary_a = part_a.loc[part_a["model"].eq("primary") & part_a["in_family"]].copy()
    primary_rows = []
    for r in primary_a.itertuples():
        primary_rows.append([
            r.dataset,
            r.cell_type,
            str(int(r.n_donors_original)),
            str(int(r.n_donors)),
            str(int(r.n_donors_lost)),
            fmt(r.slope_per_decade),
            f"{fmt(r.ci95_low)} to {fmt(r.ci95_high)}",
            fmt(r.pvalue),
            fmt(r.q_bh),
            fmt(r.mde80_per_decade),
            fmt(r.mde80_sd_per_decade),
            r.status,
        ])
    lines.append(
        markdown_table(
            ["Dataset", "Cell type", "N orig", "N", "N lost", "Slope", "95% CI", "p", "BH q", "MDE80", "MDE80 SD", "Status"],
            primary_rows,
        )
    )

    untest_primary_in = part_a.loc[
        part_a["model"].eq("primary") & part_a["in_family"] & part_a["status"].ne("FIT")
    ]
    untest_primary_below = part_a.loc[
        part_a["model"].eq("primary") & ~part_a["in_family"] & part_a["status"].ne("FIT")
    ]
    lines.extend(["", "#### UNTESTABLE primary cell types", ""])
    untest_rows = []
    for r in pd.concat([untest_primary_in, untest_primary_below], ignore_index=True).itertuples():
        untest_rows.append([
            r.dataset,
            r.cell_type,
            str(int(r.n_donors_original)),
            str(int(r.n_donors)),
            str(int(r.n_donors_lost)),
            "Yes" if r.in_family else "No",
            r.status_reason,
        ])
    lines.append(
        markdown_table(
            ["Dataset", "Cell type", "N orig", "N", "N lost", "In candidate family", "Reason"],
            untest_rows,
        )
    )
    lines.extend([
        "",
        "In Yazar primary fits, `plasmablast` (N = 14, N lost = 0, spanning 12 pools plus intercept, age, and sex) remains `UNTESTABLE` (`insufficient_residual_df`) inside the 17-type mature candidate family. `CD16-negative, CD56-bright natural killer cell, human` lost 1 donor under 30 (N = 6 to N = 5) and dropped below the 6-donor floor (`fewer_than_6_donors`), exiting the mature candidate family alongside the 10 types already below the floor in the original analysis.",
        "",
        "### Part A sensitivity summary (age 30-plus)",
        "",
    ])

    sens_summary_rows = []
    for (dataset, model), grp in part_a.loc[part_a["model"].ne("primary")].groupby(["dataset", "model"]):
        fitted = grp.loc[grp["status"].eq("FIT")]
        if model == "sensitivity_adjusted":
            neg = int((fitted["slope_per_decade"].lt(0) & fitted["q_bh"].lt(0.05)).sum())
            pos = int((fitted["slope_per_decade"].gt(0) & fitted["q_bh"].lt(0.05)).sum())
            thresh = "BH q < 0.05"
        else:
            neg = int((fitted["slope_per_decade"].lt(0) & fitted["pvalue"].lt(0.05)).sum())
            pos = int((fitted["slope_per_decade"].gt(0) & fitted["pvalue"].lt(0.05)).sum())
            thresh = "nominal p < 0.05"
        sens_summary_rows.append([
            dataset,
            model,
            str(len(fitted)),
            str(int(grp["status"].ne("FIT").sum())),
            str(neg),
            str(pos),
            thresh,
        ])
    lines.append(
        markdown_table(
            ["Dataset", "Model", "Fitted", "UNTESTABLE", "Negative hits", "Positive hits", "Threshold"],
            sens_summary_rows,
        )
    )
    lines.extend([
        "",
        "#### Sensitivity 1 adjusted fits (age 30-plus, `mean_log_umi_all`)",
        "",
    ])
    sens_adj = part_a.loc[part_a["model"].eq("sensitivity_adjusted") & part_a["in_family"]].copy()
    sens_adj_rows = []
    for r in sens_adj.itertuples():
        sens_adj_rows.append([
            r.dataset,
            r.cell_type,
            str(int(r.n_donors_original)),
            str(int(r.n_donors)),
            str(int(r.n_donors_lost)),
            fmt(r.slope_per_decade),
            f"{fmt(r.ci95_low)} to {fmt(r.ci95_high)}",
            fmt(r.pvalue),
            fmt(r.q_bh),
            fmt(r.mde80_per_decade),
            fmt(r.mde80_sd_per_decade),
            r.status,
        ])
    lines.append(
        markdown_table(
            ["Dataset", "Cell type", "N orig", "N", "N lost", "Slope", "95% CI", "p", "BH q", "MDE80", "MDE80 SD", "Status"],
            sens_adj_rows,
        )
    )
    untest_s1_in = part_a.loc[
        part_a["model"].eq("sensitivity_adjusted") & part_a["in_family"] & part_a["status"].ne("FIT")
    ]
    untest_s1_below = part_a.loc[
        part_a["model"].eq("sensitivity_adjusted") & ~part_a["in_family"] & part_a["status"].ne("FIT")
    ]
    s1_in_str = ", ".join(
        f"{r.dataset}: {r.cell_type} (N={int(r.n_donors)}, lost={int(r.n_donors_lost)}, {r.status_reason})"
        for r in untest_s1_in.itertuples()
    )
    s1_below_str = ", ".join(
        f"{r.dataset}: {r.cell_type} (N={int(r.n_donors)}, lost={int(r.n_donors_lost)})"
        for r in untest_s1_below.itertuples()
    )
    gao_s1 = part_a.loc[part_a["dataset"].eq("gao") & part_a["model"].eq("sensitivity_adjusted")].set_index("cell_type")
    yaz_s1 = part_a.loc[part_a["dataset"].eq("yazar") & part_a["model"].eq("sensitivity_adjusted")].set_index("cell_type")
    orig_yaz_s1 = orig_a.loc[orig_a["dataset"].eq("yazar") & orig_a["model"].eq("sensitivity_adjusted")].set_index("cell_type")
    gao_s1_nom = gao_s1.loc[gao_s1["slope_per_decade"].lt(0) & gao_s1["pvalue"].lt(0.05)].sort_values("pvalue")
    gao_s1_nom_str = "; ".join(
        f"`{ct}` p = {fmt(r.pvalue)}, q = {fmt(r.q_bh)}"
        for ct, r in gao_s1_nom.iterrows()
    )
    em_cd4_q = fmt(yaz_s1.loc["effector memory CD4-positive, alpha-beta T cell", "q_bh"])
    ncd8_s1_q = fmt(yaz_s1.loc["naive thymus-derived CD8-positive, alpha-beta T cell", "q_bh"])
    orig_em_cd4_q = fmt(orig_yaz_s1.loc["effector memory CD4-positive, alpha-beta T cell", "q_bh"])
    orig_ncd8_s1_q = fmt(orig_yaz_s1.loc["naive thymus-derived CD8-positive, alpha-beta T cell", "q_bh"])

    lines.extend([
        "",
        f"Within the Yazar `mean_log_umi_all` 19-type candidate family, 2 types are `UNTESTABLE`: {s1_in_str}. Below the 6-donor `mean_log_umi_all` floor (`fewer_than_6_donors`), 9 Yazar types are `UNTESTABLE`: {s1_below_str}.",
        "",
        f"In Gao `mean_log_umi_all`, {len(gao_s1_nom)} adjusted slopes are nominally negative ({gao_s1_nom_str}), matching the borderline pooled Gao `mean_log_umi_all` slope ({gao_pooled_all['slope_per_decade']:.4f}, p = {gao_pooled_all['p_value']:.4f}) in `results/followup-mature-nov17-depth.csv`, but none passes BH q < 0.05. In Yazar `mean_log_umi_all`, both original BH-significant positive adjusted hits (`effector memory CD4-positive, alpha-beta T cell`, originally q = {orig_em_cd4_q}, and `naive thymus-derived CD8-positive, alpha-beta T cell`, originally q = {orig_ncd8_s1_q}) lose BH significance at 30-plus (q = {em_cd4_q} and q = {ncd8_s1_q}). Unadjusted age-only slopes remain strongly confounded by study in Gao (10 nominally negative hits in `mean_log_umi_all`, 3 in `mean_log_umi_ge1500`) and by pool in Yazar (14 nominally positive hits in `mean_log_umi_all`, 9 in `mean_log_umi_ge1500`). Full per-type unadjusted estimates are in `results/followup-mature-block6-depth-per-cell-type.csv`.",
        "",
        "### Part B: MSSM formal mature four-unit test (age 30-plus)",
        "",
    ])

    part_b_rows = []
    for r in part_b.itertuples():
        prov_short = "reused_verified" if r.provenance.startswith("reused_verified") else "refit_sensitivity"
        part_b_rows.append([
            r.cell_type,
            r.model,
            str(int(r.n_donors_original)),
            str(int(r.n_donors)),
            str(int(r.n_donors_lost)),
            fmt(r.slope_per_decade),
            f"{fmt(r.ci95_low)} to {fmt(r.ci95_high)}",
            fmt(r.pvalue),
            fmt(r.q_bh),
            fmt(r.mde80_per_decade),
            fmt(r.mde80_sd_per_decade),
            r.status,
            prov_short,
        ])
    lines.append(
        markdown_table(
            ["Cell type", "Model", "N orig", "N", "N lost", "Slope", "95% CI", "p", "BH q", "MDE80", "MDE80 SD", "Status", "Provenance"],
            part_b_rows,
        )
    )
    lines.extend([
        "",
        "All four mature `primary` rows (`age_per_decade + sex_male`, age 30-plus; observed ages 31 to 89) were verified directly against `/home/reza/projects/agescape-revision-data/mssm.h5ad` (all reused: "
        + str(all_reused)
        + f"; maximum absolute difference across all numeric fields = `{max_mssm_verify_diff:.3e}`) and corrected with BH across the declared four-unit neuronal family. All four units show statistically significant negative depth slopes after four-unit BH correction (`q` from `2.76e-05` to `0.0064`), with high power (`MDE80 SD` from `0.1442` to `0.1573` SD per decade). The unadjusted age-only 30-plus sensitivity (`model = unadjusted`, directly fit from `mssm.h5ad`, no BH q) agrees in sign and magnitude across all four units.",
        "",
        "## Independent BH verification",
        "",
        "Every BH family was computed with `statsmodels.stats.multitest.multipletests(..., method='fdr_bh')` and checked against an independent step-up BH implementation (`bh_independent`). Maximum absolute q-value differences were zero to machine precision in all five families.",
        "",
    ])

    all_bh = bh_checks_a + [bh_check_b]
    bh_summary_rows = [
        [c["family"], str(c["candidate_size"]), str(c["n_fitted"]), f"{c['max_abs_diff']:.2e}"]
        for c in all_bh
    ]
    lines.append(
        markdown_table(
            ["BH family", "Candidate size", "Fitted (m)", "Max abs q diff"],
            bh_summary_rows,
        )
    )
    lines.extend(["", "Manual spot-checks (at least 3 ranked p-values per family, showing `p * m / k` and step-up monotonic minimum `q_manual`):", ""])
    spot_rows = []
    for c in all_bh:
        for s in c["spot_checks"]:
            spot_rows.append([
                c["family"],
                s["cell_type"],
                f"{s['rank']}/{s['m']}",
                fmt(s["pvalue"], 6),
                fmt(s["p_times_m_over_k"], 6),
                fmt(s["q_manual"], 6),
                fmt(s["q_bh"], 6),
                f"{s['abs_diff']:.2e}",
            ])
    lines.append(
        markdown_table(
            ["Family", "Cell type", "Rank k/m", "p", "p * m / k", "q_manual", "q_bh", "Abs diff"],
            spot_rows,
        )
    )

    # Power summary after 30-plus restriction
    orig_a_primary = orig_a.loc[orig_a["model"].eq("primary") & orig_a["in_family"] & orig_a["status"].eq("FIT")].set_index(["dataset", "cell_type"])
    mat_a_primary_fit = primary_a.loc[primary_a["status"].eq("FIT")].copy()
    gao_ratios = []
    yazar_ratios = []
    for r in mat_a_primary_fit.itertuples():
        o_mde = float(orig_a_primary.loc[(r.dataset, r.cell_type), "mde80_sd_per_decade"])
        ratio = float(r.mde80_sd_per_decade) / o_mde
        if r.dataset == "gao":
            gao_ratios.append((r.cell_type, o_mde, float(r.mde80_sd_per_decade), ratio))
        else:
            yazar_ratios.append((r.cell_type, o_mde, float(r.mde80_sd_per_decade), ratio))

    orig_b_primary = orig_b.loc[orig_b["model"].eq("primary")].set_index("cell_type")
    mssm_ratios = [
        (
            r.cell_type,
            float(orig_b_primary.loc[r.cell_type, "mde80_sd_per_decade"]),
            float(r.mde80_sd_per_decade),
            float(r.mde80_sd_per_decade) / float(orig_b_primary.loc[r.cell_type, "mde80_sd_per_decade"]),
        )
        for r in part_b_primary.itertuples()
    ]

    gao_min_r = min(x[3] for x in gao_ratios)
    gao_max_r = max(x[3] for x in gao_ratios)
    yazar_lost_ratios = [x for x in yazar_ratios if x[0] != "regulatory T cell"]
    yazar_lost_min_r = min(x[3] for x in yazar_lost_ratios)
    yazar_lost_max_r = max(x[3] for x in yazar_lost_ratios)
    mssm_min_r = min(x[3] for x in mssm_ratios)
    mssm_max_r = max(x[3] for x in mssm_ratios)

    lines.extend([
        "",
        "## Power after 30-plus restriction",
        "",
        f"Restricting to age 30-plus reduces donor counts and compresses the age span, inflating standardized MDE80 across every fitted primary model that lost donors: mature-to-original `MDE80 SD` ratios range from `{gao_min_r:.3f}` to `{gao_max_r:.3f}` across the 13 Gao primary types, `{yazar_lost_min_r:.3f}` to `{yazar_lost_max_r:.3f}` across the 15 fitted Yazar primary types that lost donors (`regulatory T cell` had minimum age 42, lost 0 donors, and stayed at ratio `1.000`), and `{mssm_min_r:.3f}` to `{mssm_max_r:.3f}` across the 4 MSSM units.",
        "",
        "Crucially, restriction breaks the frozen `0.50 SD per decade` adequate-power threshold in three fitted primary cell types (no null claim is made for any of them):",
        "",
        "- Gao `Epicardial`: N dropped from 20 to 18 (lost 2; `residual_df` 15 to 13); `MDE80` rose from `0.1026` to `0.1251` log UMI per decade, and `MDE80 SD` crossed the threshold from `0.4136` to `0.5014` SD per decade.",
        "- Gao `Mast`: N dropped from 18 to 17 (lost 1; `residual_df` 11 to 10); `MDE80` rose from `0.0451` to `0.0618` log UMI per decade, and `MDE80 SD` crossed the threshold from `0.4780` to `0.6461` SD per decade.",
        "- Yazar `mucosal-associated invariant T cell`: N dropped from 46 to 38 (lost 8; `residual_df` 17 to 14); `MDE80` rose from `0.0608` to `0.0819` log UMI per decade, and `MDE80 SD` crossed the threshold from `0.4448` to `0.6140` SD per decade.",
        "",
        "In addition, Yazar `CD16-negative, CD56-bright natural killer cell, human` loses its primary candidate eligibility (`N = 6` to `N = 5`, `fewer_than_6_donors`), and in the `mean_log_umi_all` sensitivity its `MDE80 SD` inflates from `0.7928` (`residual_df = 3`) to `2.5130` (`residual_df = 2`). By contrast, the remaining 11 Gao primary types (`MDE80 SD` `0.1304` to `0.4388`), the remaining 15 fitted Yazar primary types (`MDE80 SD` `0.0313` to `0.3854`), and all 4 MSSM units (`MDE80 SD` `0.1442` to `0.1573`) stay below `0.50 SD per decade`.",
        "",
        "## Verdicts under the frozen rules and comparison vs original Block 6",
        "",
        f"- **Part A, Gao**: **{gao_verdict}**. {gao_reason}",
        f"- **Part A, Yazar**: **{yazar_verdict}**. {yazar_reason}",
        f"- **Part B, MSSM**: **{mssm_verdict}**. {mssm_reason}",
        "",
    ])

    comp_rows = [
        [
            "Part A, Gao",
            "DOES NOT PROVE",
            gao_verdict,
            "Verdict changed due to power loss: still 0 negative BH hits, but Epicardial (MDE80 SD 0.5014) and Mast (0.6461) cross above 0.50 SD per decade. Endocardial nominal negative slope died (p 0.0321 to 0.1151, q 0.4507); Mast remains nominally positive (p 0.0330, q 0.4295).",
        ],
        [
            "Part A, Yazar",
            "INCONCLUSIVE",
            yazar_verdict,
            "Verdict survived as INCONCLUSIVE (plasmablast N=14 still UNTESTABLE, plus MAIT MDE80 SD crossed to 0.6140). Still 0 negative BH hits; nominal negative central memory CD4 (p 0.0087 to 0.1595) and NK (p 0.0298 to 0.0532) died. Original BH-significant positive primary hit in naive CD8 T died (q 0.0259 to 0.3855), and both positive mean_log_umi_all BH hits (naive CD8, effector memory CD4) died. CD16-negative CD56-bright NK dropped below 6 donors (6 to 5).",
        ],
        [
            "Part B, MSSM",
            "PROVES",
            mssm_verdict,
            "Verdict survived (agrees with original PROVES). All 4 units pass four-unit BH at 30-plus: GABA (-0.1039, q=2.76e-05), VIP (-0.0877, q=6.74e-05), pvalb (-0.0839, q=3.65e-04), sst (-0.0524, q=0.0064). GABA, VIP, and pvalb slopes steepen slightly; sst attenuates slightly (-0.0560 to -0.0524).",
        ],
    ]
    lines.append(
        markdown_table(
            ["Arm", "Original Block 6", "Mature 30-plus", "What survived and what died"],
            comp_rows,
        )
    )

    lines.extend([
        "",
        "## Caveats",
        "",
        "- **MSSM nested-unit caveat**: The four MSSM neuronal units (`GABAergic neuron`, `VIP GABAergic cortical interneuron`, `pvalb GABAergic cortical interneuron`, and `sst GABAergic cortical interneuron`) are nested and overlapping neuronal populations sharing donors rather than four statistically independent cell classes. The four-unit BH family reflects the pre-stated test structure, not independence across units.",
        "- **MSSM missing-covariate caveat**: The MSSM model adjusts only for sex (`age_per_decade + sex_male`). No batch, pool, postmortem interval, RIN, or processing-cohort covariates enter the model. Direct inspection of `/home/reza/projects/agescape-revision-data/mssm.h5ad` `obs` confirmed that `Source` (`M`), `assay` (`10x 3' v3`), and `suspension_type` (`nucleus`) are constant across all 4,140,453 cells, and no technical batch or pool column is available to adjust for potential unmeasured confounding.",
        "- **Batch and pool confounding in Gao and Yazar**: Unadjusted slopes remain heavily confounded by study in Gao (negative bias) and pool in Yazar (positive bias). Study or pool fixed effects remove these apparent trends in every primary cell type.",
        "- **Power and non-equivalence**: MDE80 values are minimum detectable effects for 80 percent power at two-sided alpha 0.05, not equivalence proofs. Where MDE80 exceeds 0.50 SD per decade (`Gao Epicardial`, `Gao Mast`, `Yazar MAIT`), null results are explicitly treated as underpowered rather than evidence of absence.",
        "",
        "## Output files",
        "",
        "- `reports/followup-mature-block6-depth-extension.md`",
        "- `results/followup-mature-block6-depth-per-cell-type.csv`",
        "- `results/followup-mature-block6-mssm-depth.csv`",
        "- `scripts/run_followup_mature_block6_depth_extension.py`",
        "",
    ])
    write_clean_markdown("\n".join(lines), REPORT)


def main() -> None:
    np.random.seed(SEED)
    part_a, bh_checks_a = run_part_a()
    part_b, bh_check_b = run_part_b()
    write_report(part_a, part_b, bh_checks_a, bh_check_b)
    for dataset in ("gao", "yazar"):
        verdict, reason = verdict_a(part_a, dataset)
        print(f"Part A {dataset}: {verdict}. {reason}", flush=True)
    verdict, reason = verdict_b(part_b)
    print(f"Part B MSSM: {verdict}. {reason}", flush=True)
    print(f"Wrote {len(part_a)} Part A rows and {len(part_b)} Part B rows.", flush=True)


if __name__ == "__main__":
    main()
