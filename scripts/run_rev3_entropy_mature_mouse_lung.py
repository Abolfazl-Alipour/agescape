#!/usr/bin/env python3
"""Item 4: Mature-only (>= 18 months) transcriptional entropy in mouse lung.

Evaluates age associations on Tabula Muris Senis Lung cells from mice aged
>= 18 months (18m, 21m, 30m; 10 mice total) using:
  1. Shannon entropy H(p) after removing MT/RP/sex genes (Agescape protocol).
  2. De Man et al. 2026 published Methods KL divergence from uniform:
     D_KL(p || u) = ln(G) - H(p) using raw counts over expressed genes.
  3. De Man et al. 2026 public GitHub variant:
     D_KL(p || u) after square-root transformation of UMI counts.

Run under scripts/resource_guard.py with a 30 GB free reserve and 8 GB cap.
"""

from __future__ import annotations

import os
from pathlib import Path
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import scipy.stats as stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_h3_transcriptional_disorder import AGE_MAP, TARGET_CELL_TYPES, excluded_gene_mask

RESULTS = REPO / "results"
SOURCE_CANDIDATES = [
    REPO / "data" / "interim" / "m1-qc-cleaned.h5ad",
    Path("/home/reza/projects/agescape/data/interim/m1-qc-cleaned.h5ad"),
]
DEPTH_SPECS = (
    ("500", 500),
    ("1500", 1500),
    ("3000", 3000),
    ("full", 1500),
)
METRICS = (
    "shannon_entropy",
    "deman_kl_methods",
    "deman_kl_sqrt",
)


def resolve_input_h5ad() -> Path:
    for path in SOURCE_CANDIDATES:
        if path.exists():
            return path
    raise FileNotFoundError(f"Could not locate m1-qc-cleaned.h5ad in {SOURCE_CANDIDATES}")


def compute_shannon(matrix: sp.csr_matrix) -> np.ndarray:
    matrix = sp.csr_matrix(matrix, dtype=np.float64)
    totals = np.asarray(matrix.sum(axis=1)).ravel()
    wlogs = matrix.copy()
    positive = wlogs.data > 0
    wlogs.data = np.where(positive, wlogs.data * np.log(np.maximum(wlogs.data, 1e-300)), 0.0)
    moments = np.asarray(wlogs.sum(axis=1)).ravel()
    out = np.full(matrix.shape[0], np.nan, dtype=np.float64)
    valid = totals > 0
    out[valid] = np.log(totals[valid]) - moments[valid] / totals[valid]
    return out


def compute_deman_kl(matrix: sp.csr_matrix, sqrt_counts: bool = False) -> np.ndarray:
    """Compute D_KL(p || u) = ln(G) - H(p) per De Man et al. 2026."""
    mat = sp.csr_matrix(matrix, dtype=np.float64).copy()
    if sqrt_counts:
        mat.data = np.sqrt(mat.data)
    totals = np.asarray(mat.sum(axis=1)).ravel()
    out = np.full(mat.shape[0], np.nan, dtype=np.float64)
    valid = totals > 0
    n_genes = mat.shape[1]
    sub = mat[valid].copy()
    sub_totals = totals[valid]
    sub.data /= np.repeat(sub_totals, np.diff(sub.indptr))
    pos = sub.data > 0
    sub.data[pos] *= np.log(sub.data[pos] * n_genes)
    out[valid] = np.asarray(sub.sum(axis=1)).ravel()
    return out


def fit_donor_celltype_ols(
    donors_df: pd.DataFrame,
    metric: str,
    depth: str,
    model_name: str,
    family_name: str,
    min_donors: int = 6,
) -> list[dict]:
    rows: list[dict] = []
    for cell_type, part in donors_df.dropna(subset=[metric]).groupby("cell_type"):
        part = part.copy()
        if len(part) < min_donors:
            continue
        if model_name == "age_only":
            predictors = ["age_months"]
        elif model_name == "age_sex":
            predictors = ["age_months", "sex_male"]
        elif model_name == "age_sex_library":
            predictors = ["age_months", "sex_male", "mean_log_umi_actual"]
        else:
            raise ValueError(f"Unknown model {model_name}")

        used_predictors = list(predictors)
        X = sm.add_constant(part[used_predictors].astype(float), has_constant="add")
        actual_model = model_name
        if "sex_male" in used_predictors and (
            part["sex_male"].nunique() < 2 or np.linalg.matrix_rank(X.to_numpy()) < X.shape[1]
        ):
            used_predictors.remove("sex_male")
            X = sm.add_constant(part[used_predictors].astype(float), has_constant="add")
            actual_model = f"{model_name}_rank_deficient_sex"

        if np.linalg.matrix_rank(X.to_numpy()) < X.shape[1] or len(part) - X.shape[1] < 1:
            continue

        y = part[metric].astype(float)
        fit = sm.OLS(y, X).fit()
        df_resid = int(fit.df_resid)
        slope = float(fit.params["age_months"])
        se = float(fit.bse["age_months"])
        tcrit_975 = float(stats.t.ppf(0.975, df_resid))
        tcrit_80 = float(stats.t.ppf(0.80, df_resid))
        mde_per_month = (tcrit_975 + tcrit_80) * se
        ci_low = slope - tcrit_975 * se
        ci_high = slope + tcrit_975 * se
        metric_sd = float(y.std(ddof=1))
        rows.append(
            {
                "scope": "per_cell_type",
                "family": family_name,
                "depth": str(depth),
                "metric": metric,
                "model": actual_model,
                "cell_type": cell_type,
                "is_manuscript_target": cell_type in TARGET_CELL_TYPES,
                "n_donors": int(len(part)),
                "n_cell_types": 1,
                "n_donor_units": int(len(part)),
                "n_cells": int(part["n_cells"].sum()),
                "df_resid": df_resid,
                "age_min_months": float(part["age_months"].min()),
                "age_max_months": float(part["age_months"].max()),
                "slope_per_month": slope,
                "se_per_month": se,
                "ci95_low_per_month": ci_low,
                "ci95_high_per_month": ci_high,
                "slope_per_12mo": slope * 12.0,
                "se_per_12mo": se * 12.0,
                "ci95_low_per_12mo": ci_low * 12.0,
                "ci95_high_per_12mo": ci_high * 12.0,
                "effect_15mo_sd": (slope * 15.0 / metric_sd) if metric_sd > 0 else np.nan,
                "se_15mo_sd": (se * 15.0 / metric_sd) if metric_sd > 0 else np.nan,
                "mde80_per_month": mde_per_month,
                "mde80_per_12mo": mde_per_month * 12.0,
                "mde80_15mo_sd": (mde_per_month * 15.0 / metric_sd) if metric_sd > 0 else np.nan,
                "pvalue": float(fit.pvalues["age_months"]),
                "r2": float(fit.rsquared),
            }
        )

    if rows:
        qvals = multipletests([r["pvalue"] for r in rows], method="fdr_bh")[1]
        for r, q in zip(rows, qvals):
            r["q_bh_within_depth_family"] = float(q)
    return rows


def fit_global_clustered(
    donors_df: pd.DataFrame,
    metric: str,
    depth: str,
    model_name: str,
    family_name: str,
    min_donors_per_ct: int = 6,
) -> dict | None:
    data = donors_df.dropna(subset=[metric]).copy()
    parts = []
    for cell_type, part in data.groupby("cell_type"):
        if len(part) < min_donors_per_ct:
            continue
        sd = float(part[metric].std(ddof=1))
        if not np.isfinite(sd) or sd <= 0:
            continue
        part = part.copy()
        part["outcome"] = (part[metric] - part[metric].mean()) / sd
        part["age15"] = (part["age_months"] - 18.0) / 15.0
        part["sex_centered"] = part["sex_male"] - part["sex_male"].mean()
        part["depth_centered"] = part["mean_log_umi_actual"] - part["mean_log_umi_actual"].mean()
        parts.append(part)
    if not parts:
        return None
    df = pd.concat(parts, ignore_index=True)
    dummies = pd.get_dummies(df["cell_type"], drop_first=True, dtype=float)
    if model_name == "age_only":
        predictors = ["age15"]
    elif model_name == "age_sex":
        predictors = ["age15", "sex_centered"]
    elif model_name == "age_sex_library":
        predictors = ["age15", "sex_centered", "depth_centered"]
    else:
        raise ValueError(f"Unknown model {model_name}")

    X = sm.add_constant(pd.concat([df[predictors], dummies], axis=1), has_constant="add")
    fit = sm.OLS(df["outcome"], X).fit(
        cov_type="cluster", cov_kwds={"groups": df["donor"]}, use_t=True
    )
    df_cluster = int(df["donor"].nunique() - 1)
    est = float(fit.params["age15"])
    se = float(fit.bse["age15"])
    ci = fit.conf_int().loc["age15"]
    tcrit_975 = float(stats.t.ppf(0.975, df_cluster))
    tcrit_80 = float(stats.t.ppf(0.80, df_cluster))
    mde_15mo_sd = (tcrit_975 + tcrit_80) * se
    return {
        "scope": "global_donor_clustered",
        "family": family_name,
        "depth": str(depth),
        "metric": metric,
        "model": model_name,
        "cell_type": "ALL_ELIGIBLE",
        "is_manuscript_target": family_name == "manuscript_target_cell_types",
        "n_donors": int(df["donor"].nunique()),
        "n_cell_types": int(df["cell_type"].nunique()),
        "n_donor_units": int(len(df)),
        "n_cells": int(df["n_cells"].sum()),
        "df_resid": df_cluster,
        "age_min_months": float(df["age_months"].min()),
        "age_max_months": float(df["age_months"].max()),
        "slope_per_month": np.nan,
        "se_per_month": np.nan,
        "ci95_low_per_month": np.nan,
        "ci95_high_per_month": np.nan,
        "slope_per_12mo": est * (12.0 / 15.0),
        "se_per_12mo": se * (12.0 / 15.0),
        "ci95_low_per_12mo": float(ci.iloc[0]) * (12.0 / 15.0),
        "ci95_high_per_12mo": float(ci.iloc[1]) * (12.0 / 15.0),
        "effect_15mo_sd": est,
        "se_15mo_sd": se,
        "mde80_per_month": np.nan,
        "mde80_per_12mo": mde_15mo_sd * (12.0 / 15.0),
        "mde80_15mo_sd": mde_15mo_sd,
        "pvalue": float(fit.pvalues["age15"]),
        "r2": float(fit.rsquared),
        "q_bh_within_depth_family": np.nan,
    }


def main() -> None:
    RESULTS.mkdir(exist_ok=True)
    h5_path = resolve_input_h5ad()
    print(f"[item4] Reading mature mouse lung from {h5_path}", flush=True)
    atlas = ad.read_h5ad(h5_path, backed="r")
    obs = atlas.obs
    age_months = obs["age"].astype(str).map(AGE_MAP)
    mask = (
        (obs["tissue"] == "Lung")
        & (~obs["predicted_doublet"].to_numpy(dtype=bool))
        & (age_months >= 18.0).to_numpy(dtype=bool)
        & (obs["n_counts"] >= 500).to_numpy(dtype=bool)
        & (obs["cell_ontology_class"].astype(str) != "unknown").to_numpy(dtype=bool)
    )
    lung = atlas[mask].to_memory()
    lung.obs["age_months"] = lung.obs["age"].astype(str).map(AGE_MAP).astype(float)
    lung.obs["sex_male"] = (lung.obs["sex"].astype(str).str.lower() == "male").astype(float)
    lung.obs["mouse.id"] = lung.obs["mouse.id"].astype(str)
    lung.obs["actual_n_counts"] = np.asarray(lung.X.sum(axis=1)).ravel()
    print(
        f"[item4] Loaded {lung.n_obs} mature lung singlet cells across "
        f"{lung.obs['mouse.id'].nunique()} mice "
        f"(ages={sorted(lung.obs['age_months'].unique().tolist())})",
        flush=True,
    )

    gene_names = lung.var_names.astype(str).to_numpy()
    technical, _ = excluded_gene_mask(gene_names)

    donor_records: list[dict] = []
    for depth_label, min_umi in DEPTH_SPECS:
        selected = lung[lung.obs["actual_n_counts"] >= min_umi].copy()
        if depth_label == "full":
            matrix = sp.csr_matrix(selected.X)
        else:
            sc.pp.downsample_counts(
                selected, counts_per_cell=int(depth_label), replace=False, random_state=42
            )
            matrix = sp.csr_matrix(selected.X)
            totals_check = np.asarray(matrix.sum(axis=1)).ravel()
            assert np.all(totals_check == int(depth_label))

        shannon_mat = matrix[:, ~technical]
        nonzero_genes = np.asarray(matrix.sum(axis=0)).ravel() > 0
        deman_mat = matrix[:, nonzero_genes]
        selected.obs["shannon_entropy"] = compute_shannon(shannon_mat)
        selected.obs["deman_kl_methods"] = compute_deman_kl(deman_mat, sqrt_counts=False)
        selected.obs["deman_kl_sqrt"] = compute_deman_kl(deman_mat, sqrt_counts=True)
        selected.obs["deman_n_genes_G"] = int(deman_mat.shape[1])

        for cell_type, ct_obs in selected.obs.groupby("cell_ontology_class", observed=True):
            for donor, d_obs in ct_obs.groupby("mouse.id", observed=True):
                if len(d_obs) < 20:
                    continue
                donor_records.append(
                    {
                        "depth": depth_label,
                        "cell_type": str(cell_type),
                        "is_manuscript_target": str(cell_type) in TARGET_CELL_TYPES,
                        "donor": str(donor),
                        "age_months": float(d_obs["age_months"].iloc[0]),
                        "sex_male": float(d_obs["sex_male"].iloc[0]),
                        "n_cells": int(len(d_obs)),
                        "mean_umi_actual": float(d_obs["actual_n_counts"].mean()),
                        "mean_log_umi_actual": float(np.log(d_obs["actual_n_counts"]).mean()),
                        "deman_n_genes_G": int(d_obs["deman_n_genes_G"].iloc[0]),
                        "shannon_entropy": float(d_obs["shannon_entropy"].mean()),
                        "deman_kl_methods": float(d_obs["deman_kl_methods"].mean()),
                        "deman_kl_sqrt": float(d_obs["deman_kl_sqrt"].mean()),
                    }
                )
        print(f"[item4] Completed depth {depth_label}", flush=True)

    donors_df = pd.DataFrame(donor_records)
    donors_out = RESULTS / "rev3-entropy-mature-only-mouse-lung-donors.csv"
    donors_df.to_csv(donors_out, index=False)
    print(f"[item4] Wrote {donors_out} ({len(donors_df)} rows)", flush=True)

    result_rows: list[dict] = []
    for depth_label, _ in DEPTH_SPECS:
        part = donors_df[donors_df["depth"] == depth_label]
        manuscript_part = part[part["is_manuscript_target"]]
        for family_name, fam_df, min_d in (
            ("manuscript_target_n6", manuscript_part, 6),
            ("manuscript_target_include_fibroblast_n5", manuscript_part, 5),
            ("all_mature_eligible_n6", part, 6),
            ("all_mature_eligible_n5", part, 5),
        ):
            for metric in METRICS:
                for model_name in ("age_sex", "age_sex_library", "age_only"):
                    g = fit_global_clustered(
                        fam_df, metric, depth_label, model_name, family_name, min_donors_per_ct=min_d
                    )
                    if g is not None:
                        result_rows.append(g)
                    result_rows.extend(
                        fit_donor_celltype_ols(
                            fam_df, metric, depth_label, model_name, family_name, min_donors=min_d
                        )
                    )

    out_df = pd.DataFrame(result_rows)
    out_path = RESULTS / "rev3-entropy-mature-only-mouse-lung.csv"
    out_df.to_csv(out_path, index=False)
    print(f"[item4] Wrote {out_path} ({len(out_df)} rows)", flush=True)


if __name__ == "__main__":
    main()
