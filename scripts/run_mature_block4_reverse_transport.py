#!/usr/bin/env python3
"""Mature Block 4 reverse signature transport."""

from __future__ import annotations

import gc
import os
import sys
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_key, "1")

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_followup_11pt_block4_reverse_transport import (
    EXPECTED_TISSUES,
    GAO_H5AD,
    ITERATIONS,
    MIN_CELLS,
    MIN_DONORS,
    YAZAR_H5AD,
    derive_gao,
    derive_yazar,
    fit_gene_models,
    human_biotypes,
    make_signature,
    wait_for_memory,
)
from scripts.run_followup_between_subject_gao_yazar import build_human_gene_projection
from src.followup_muscle import NULL_SEED, fit_ols
from src.followup_organs import run_family_c_p1_signature, write_clean_csv, write_clean_markdown

RESULTS = REPO / "results"
REPORTS = REPO / "reports"
MOUSE_DIR = REPO / "data" / "interim"
GENE_OUT = RESULTS / "followup-mature-block4-reverse-signature-genes.csv"
TRANSPORT_OUT = RESULTS / "followup-mature-block4-reverse-transport.csv"
RUNTIME_OUT = RESULTS / "followup-mature-block4-runtime.csv"
REPORT_OUT = REPORTS / "followup-mature-block4-reverse-transport.md"
AGE_MAP = {"18m": 18.0, "21m": 21.0, "24m": 24.0, "30m": 30.0}
MODEL_NAME = "kimi-code/k3.1-pro"
PRIMARY_ARM = "mature_18plus_slope"
PAIRWISE_ARM = "mature_18_vs_30"
SCORING_ARMS = (PRIMARY_ARM, PAIRWISE_ARM)
GENE_COLUMNS = [
    "derivation_cohort",
    "derivation_model",
    "human_gene_symbol",
    "mean_age_slope_log2cpm_per_decade",
    "min_pvalue",
    "min_q_bh",
    "n_discovery_cell_types",
    "discovery_cell_types",
    "direction",
    "human_biotype",
    "passes_clean_filter",
    "mouse_gene_symbol",
    "mapping_status",
    "in_transport_signature",
    "discovery_stage",
]


class RuntimeLog:
    def __init__(self) -> None:
        self.rows: list[dict] = []
        self.started = datetime.now(timezone.utc)

    @contextmanager
    def stage(self, name: str):
        start = datetime.now(timezone.utc)
        try:
            yield
        finally:
            end = datetime.now(timezone.utc)
            self.rows.append({
                "model": MODEL_NAME,
                "stage": name,
                "start_utc": start.isoformat(),
                "end_utc": end.isoformat(),
                "duration_seconds": (end - start).total_seconds(),
                "iterations": ITERATIONS,
                "seed": NULL_SEED,
            })
            print(f"Stage {name}: {(end - start).total_seconds():.3f} seconds", flush=True)

    def finish(self) -> pd.DataFrame:
        end = datetime.now(timezone.utc)
        self.rows.append({
            "model": MODEL_NAME,
            "stage": "total",
            "start_utc": self.started.isoformat(),
            "end_utc": end.isoformat(),
            "duration_seconds": (end - self.started).total_seconds(),
            "iterations": ITERATIONS,
            "seed": NULL_SEED,
        })
        return pd.DataFrame(self.rows)


def restrict_human_to_30plus(
    expr: dict[str, pd.DataFrame],
    meta: dict[str, pd.DataFrame],
    cohort: str,
) -> dict[str, float | int]:
    if set(expr) != set(meta):
        raise ValueError(f"{cohort} expression and metadata groups differ")
    for group in sorted(expr):
        group_meta = meta[group]
        if "age_years" not in group_meta:
            raise ValueError(f"{cohort} cache cannot be subset by donor age")
        ages = pd.to_numeric(group_meta["age_years"], errors="raise")
        donors = group_meta.index[ages >= 30.0].astype(str)
        expr[group] = expr[group].loc[:, donors].copy()
        meta[group] = group_meta.loc[donors].copy()
        fitted_ages = pd.to_numeric(meta[group]["age_years"], errors="raise")
        if len(fitted_ages) and fitted_ages.min() < 30.0:
            raise AssertionError(f"{cohort} {group} retained a donor under 30")
    cached_donors = set().union(*(set(frame.index.astype(str)) for frame in meta.values()))
    fit_groups = [group for group in sorted(meta) if len(meta[group]) >= MIN_DONORS]
    fitted_donors = set().union(*(set(meta[group].index.astype(str)) for group in fit_groups))
    fitted_ages = np.concatenate([
        pd.to_numeric(meta[group].loc[sorted(fitted_donors.intersection(meta[group].index)), "age_years"]).to_numpy()
        for group in fit_groups
    ])
    audit = {
        "cohort": cohort,
        "n_donors_30plus_cached": len(cached_donors),
        "n_donors_fitted": len(fitted_donors),
        "n_cell_types_cached": len(meta),
        "n_cell_types_fitted": len(fit_groups),
        "min_fitted_age_years": float(fitted_ages.min()),
        "max_fitted_age_years": float(fitted_ages.max()),
        "cache_route": "reused_original_donor_mean_cache",
    }
    print(f"{cohort} 30-plus cache audit: {audit}", flush=True)
    return audit


def derive_signature(
    tests: pd.DataFrame,
    cohort: str,
    model: str,
    orthologs: pd.DataFrame,
    biotypes: dict[str, str],
):
    discovered = tests[(tests["q_bh"] < 0.05) & (tests["effect"].abs() >= 0.25)]
    if not discovered.empty:
        return make_signature(tests, cohort, model, orthologs, biotypes)
    frame = pd.DataFrame(columns=GENE_COLUMNS)
    program = {"up": [], "down": []}
    counts = {
        "discovery": 0,
        "clean": 0,
        "mapped": 0,
        "lost": 0,
        "coverage": 0.0,
        "up": 0,
        "down": 0,
    }
    print(f"{cohort} {model} signature counts: {counts}", flush=True)
    return frame, program, counts


def mouse_targets() -> dict[str, tuple[pd.DataFrame, pd.DataFrame, str, str]]:
    paths = sorted(MOUSE_DIR.glob("m1-qc-cleaned-*.h5ad"))
    if len(paths) != 16 or {path.stem.removeprefix("m1-qc-cleaned-") for path in paths} != set(EXPECTED_TISSUES):
        raise ValueError(f"Expected 16 restored TMS files, found {len(paths)}")
    targets = {}
    for path in paths:
        tissue = path.stem.removeprefix("m1-qc-cleaned-")
        wait_for_memory(f"loading {tissue}")
        adata = ad.read_h5ad(path)
        if adata.shape != (EXPECTED_TISSUES[tissue], 20138):
            raise ValueError(f"{path}: unexpected shape {adata.shape}")
        keep = (
            (~adata.obs["predicted_doublet"].astype(bool).to_numpy())
            & (adata.obs["cell_ontology_class"].astype(str) != "unknown").to_numpy()
            & adata.obs["age"].astype(str).isin(AGE_MAP).to_numpy()
        )
        adata = adata[keep].copy()
        obs = adata.obs
        for cell_type in sorted(obs["cell_ontology_class"].astype(str).unique()):
            cell_mask = obs["cell_ontology_class"].astype(str).eq(cell_type).to_numpy()
            cell_obs = obs.loc[cell_mask]
            counts = cell_obs.groupby("mouse.id", observed=True).size()
            donors = sorted(counts[counts >= MIN_CELLS].index.astype(str))
            if len(donors) < MIN_DONORS:
                continue
            positions = np.flatnonzero(cell_mask & obs["mouse.id"].astype(str).isin(donors).to_numpy())
            donor_values = obs["mouse.id"].astype(str).to_numpy()[positions]
            donor_index = {donor: index for index, donor in enumerate(donors)}
            codes = np.array([donor_index[donor] for donor in donor_values], dtype=int)
            indicator = sp.csr_matrix(
                (np.ones(len(positions)), (np.arange(len(positions)), codes)),
                shape=(len(positions), len(donors)),
            )
            projection = build_human_gene_projection(adata.var_names.astype(str).to_numpy())
            projected = (adata.X[positions] @ projection["proj"]).tocsr()
            cell_library = np.asarray(projected.sum(axis=1)).ravel()
            scale = 1e6 / np.where(cell_library > 0, cell_library, 1.0)
            projected.data *= np.repeat(scale, np.diff(projected.indptr))
            np.log2(projected.data + 1.0, out=projected.data)
            cell_counts = counts.loc[donors].to_numpy(dtype=float)
            mean_expression = (indicator.T @ projected).toarray().T / cell_counts[None, :]
            log2cpm = pd.DataFrame(
                mean_expression,
                index=projection["audit"]["retained_symbols"],
                columns=donors,
            )
            donor_obs = cell_obs[cell_obs["mouse.id"].astype(str).isin(donors)].copy()
            donor_obs["mouse.id"] = donor_obs["mouse.id"].astype(str)
            donor_meta = donor_obs.groupby("mouse.id", observed=True).agg(
                age_months=("age", lambda values: AGE_MAP[str(values.iloc[0])]),
                sex_male=("sex", lambda values: float(str(values.iloc[0]).lower() == "male")),
            ).loc[donors]
            donor_meta["n_cells"] = cell_counts.astype(int)
            donor_meta["lib_size"] = np.bincount(codes, weights=cell_library, minlength=len(donors))
            donor_meta["log_lib"] = np.log1p(donor_meta["lib_size"])
            if donor_meta["age_months"].min() < 18.0:
                raise AssertionError(f"{tissue} {cell_type} retained a mouse under 18 months")
            key = f"{tissue}::{cell_type}"
            targets[key] = (log2cpm, donor_meta, tissue, cell_type)
        print(f"Prepared {tissue}; cumulative 18-plus eligible units {len(targets)}", flush=True)
        del adata
        gc.collect()
    if not targets:
        raise RuntimeError("No eligible 18-plus TMS units")
    return targets


def pairwise_view(expr: pd.DataFrame, meta: pd.DataFrame):
    pair_meta = meta[meta["age_months"].isin([18.0, 30.0])].copy()
    pair_meta["age_group"] = (pair_meta["age_months"] == 30.0).astype(float)
    pair_expr = expr.loc[:, pair_meta.index]
    n18 = int((pair_meta["age_months"] == 18.0).sum())
    n30 = int((pair_meta["age_months"] == 30.0).sum())
    eligible = len(pair_meta) >= MIN_DONORS and n18 >= 3 and n30 >= 3
    return pair_expr, pair_meta, eligible, n18, n30


def empty_transport_row(
    cohort: str,
    model: str,
    scoring_arm: str,
    tissue: str,
    cell_type: str,
    meta: pd.DataFrame,
    status: str,
    status_reason: str,
) -> dict:
    ages = sorted(meta["age_months"].unique()) if len(meta) else []
    return {
        "derivation_cohort": cohort,
        "derivation_model": model,
        "scoring_direction": "human_to_mouse",
        "scoring_arm": scoring_arm,
        "tissue": tissue,
        "cell_type": cell_type,
        "status": status,
        "status_reason": status_reason,
        "n_donors": len(meta),
        "n_donors_18m": int((meta["age_months"] == 18.0).sum()) if len(meta) else 0,
        "n_donors_30m": int((meta["age_months"] == 30.0).sum()) if len(meta) else 0,
        "min_donor_age_months": float(meta["age_months"].min()) if len(meta) else np.nan,
        "max_donor_age_months": float(meta["age_months"].max()) if len(meta) else np.nan,
        "fitted_mouse_ages_months": ";".join(str(int(age)) for age in ages),
        "n_cells_total": int(meta["n_cells"].sum()),
        "n_genes_mapped": 0,
        "n_up_mapped": 0,
        "n_down_mapped": 0,
        "n_random_sets": ITERATIONS,
        "random_seed": NULL_SEED,
        "passes_both": False,
    }


def binary_orientation_check(
    expr: pd.DataFrame,
    meta: pd.DataFrame,
    program: dict[str, list[str]],
) -> dict[str, float | bool]:
    log2_thresh = float(np.log2(2.0))
    min_req = max(3, int(np.ceil(0.5 * len(meta))))
    keep = (expr.to_numpy(dtype=float) >= log2_thresh).sum(axis=1) >= min_req
    mat = expr.loc[keep]
    genes = mat.index.astype(str)
    up_idx = np.flatnonzero(genes.isin(set(program["up"])))
    down_idx = np.flatnonzero(genes.isin(set(program["down"])))
    values = mat.to_numpy(dtype=float)
    scores = values[up_idx, :].mean(axis=0) - values[down_idx, :].mean(axis=0)
    groups = meta["age_group"].to_numpy(dtype=float)
    design = np.column_stack([np.ones(len(meta)), groups])
    fitted = fit_ols(scores, design, effect_column=1)
    raw_difference = float(scores[groups == 1.0].mean() - scores[groups == 0.0].mean())
    unadjusted_slope = float(fitted["effect"][0])
    error = abs(unadjusted_slope - raw_difference)
    if error > 1e-10:
        raise AssertionError(f"Binary orientation check failed: {error}")
    return {
        "binary_mean_score_difference_30_minus_18": raw_difference,
        "binary_unadjusted_slope_30_minus_18": unadjusted_slope,
        "binary_orientation_abs_error": error,
    }


def score_targets(
    targets: dict[str, tuple[pd.DataFrame, pd.DataFrame, str, str]],
    programs: dict[tuple[str, str], dict[str, list[str]]],
    runtime: RuntimeLog,
) -> pd.DataFrame:
    rows = []
    for scoring_arm in SCORING_ARMS:
        for (cohort, model), program in sorted(programs.items()):
            stage = f"score_{cohort}_{model}_{scoring_arm}"
            with runtime.stage(stage):
                print(
                    f"Scoring {cohort} {model} {scoring_arm}: planned family has {len(targets)} units",
                    flush=True,
                )
                for key in sorted(targets):
                    full_expr, full_meta, tissue, cell_type = targets[key]
                    if scoring_arm == PRIMARY_ARM:
                        expr, meta = full_expr, full_meta
                        age_col = "age_months"
                        eligible = True
                    else:
                        expr, meta, eligible, n18, n30 = pairwise_view(full_expr, full_meta)
                        age_col = "age_group"
                    if not program["up"] or not program["down"]:
                        missing = "UP" if not program["up"] else "DOWN"
                        status = f"UNTESTABLE_MISSING_{missing}_DIRECTION"
                        rows.append(empty_transport_row(
                            cohort, model, scoring_arm, tissue, cell_type, meta, status,
                            f"Derived mapped signature has no {missing.lower()} genes",
                        ))
                        continue
                    if not eligible:
                        status = "UNTESTABLE_18_VS_30_DONOR_FLOOR"
                        rows.append(empty_transport_row(
                            cohort, model, scoring_arm, tissue, cell_type, meta, status,
                            f"Requires at least 6 total and 3 per side; observed {n18} at 18m and {n30} at 30m",
                        ))
                        continue
                    result = run_family_c_p1_signature(
                        {key: expr},
                        {key: meta},
                        p1_program=program,
                        age_col=age_col,
                        candidate_covariates=("sex_male", "log_lib"),
                        min_donors=MIN_DONORS,
                        min_residual_df=2,
                        iterations=ITERATIONS,
                        seed=NULL_SEED,
                    )
                    if result.empty:
                        rows.append(empty_transport_row(
                            cohort, model, scoring_arm, tissue, cell_type, meta,
                            "UNTESTABLE_SCORING_ELIGIBILITY",
                            "Shared scorer returned no fit because an expression, detected-direction, or design floor failed",
                        ))
                        continue
                    fitted = result.iloc[0].to_dict()
                    df = int(fitted["residual_df"])
                    effect = float(fitted["composite_effect"])
                    se = float(fitted["composite_se"])
                    mde80 = (stats.t.ppf(0.975, df) + stats.t.ppf(0.80, df)) * se
                    row = empty_transport_row(
                        cohort, model, scoring_arm, tissue, cell_type, meta, "TESTED", "Both scoring tests fit"
                    )
                    row.update({
                        "n_donors": int(fitted["n_donors"]),
                        "n_genes_mapped": int(fitted["n_p1_total_detected"]),
                        "n_up_mapped": int(fitted["n_p1_up_detected"]),
                        "n_down_mapped": int(fitted["n_p1_down_detected"]),
                        "n_expression_pool": int(fitted["n_expressed_pool"]),
                        "effect_units": "log2cpm_per_month" if scoring_arm == PRIMARY_ARM else "log2cpm_30_minus_18",
                        "composite_effect": effect,
                        "composite_slope_per_month": effect if scoring_arm == PRIMARY_ARM else np.nan,
                        "composite_slope_per_12_months": effect * 12.0 if scoring_arm == PRIMARY_ARM else np.nan,
                        "composite_slope_30_minus_18": effect if scoring_arm == PAIRWISE_ARM else np.nan,
                        "composite_se": se,
                        "composite_se_per_12_months": se * 12.0 if scoring_arm == PRIMARY_ARM else np.nan,
                        "composite_ci95_low": float(fitted["composite_ci95_low"]),
                        "composite_ci95_high": float(fitted["composite_ci95_high"]),
                        "composite_ci95_low_per_12_months": (
                            float(fitted["composite_ci95_low"]) * 12.0 if scoring_arm == PRIMARY_ARM else np.nan
                        ),
                        "composite_ci95_high_per_12_months": (
                            float(fitted["composite_ci95_high"]) * 12.0 if scoring_arm == PRIMARY_ARM else np.nan
                        ),
                        "composite_pvalue": float(fitted["composite_pvalue"]),
                        "empirical_pvalue": float(fitted["empirical_null_pvalue"]),
                        "null_mean_effect": float(fitted["null_mean_effect"]),
                        "null_sd_effect": float(fitted["null_sd_effect"]),
                        "mde80": mde80,
                        "mde80_per_month": mde80 if scoring_arm == PRIMARY_ARM else np.nan,
                        "mde80_per_12_months": mde80 * 12.0 if scoring_arm == PRIMARY_ARM else np.nan,
                        "mde80_30_minus_18": mde80 if scoring_arm == PAIRWISE_ARM else np.nan,
                        "residual_df": df,
                        "predictor": age_col,
                        "covariates": fitted["covariates"],
                        "omitted_covariates": fitted["omitted_covariates"],
                    })
                    if scoring_arm == PAIRWISE_ARM:
                        orientation = binary_orientation_check(expr, meta, program)
                        row.update(orientation)
                        row["binary_adjusted_sign_matches_raw_difference"] = bool(
                            np.sign(effect) == np.sign(orientation["binary_mean_score_difference_30_minus_18"])
                        )
                    rows.append(row)
                    print(
                        f"  {key} {scoring_arm}: effect {effect:.5g}, empirical p "
                        f"{fitted['empirical_null_pvalue']:.5g}",
                        flush=True,
                    )
    result = pd.DataFrame(rows)
    result["family_composite_size"] = 0
    result["family_empirical_size"] = 0
    result["composite_q"] = np.nan
    result["empirical_q"] = np.nan
    group_columns = ["derivation_cohort", "derivation_model", "scoring_arm"]
    for _, indices in result.groupby(group_columns, sort=True).groups.items():
        family = result.loc[indices]
        tested = family.index[family["status"] == "TESTED"]
        result.loc[tested, "family_composite_size"] = len(tested)
        result.loc[tested, "family_empirical_size"] = len(tested)
        if len(tested):
            result.loc[tested, "composite_q"] = multipletests(
                result.loc[tested, "composite_pvalue"], method="fdr_bh"
            )[1]
            result.loc[tested, "empirical_q"] = multipletests(
                result.loc[tested, "empirical_pvalue"], method="fdr_bh"
            )[1]
    result["passes_both"] = (
        (result["status"] == "TESTED")
        & (result["composite_q"] < 0.05)
        & (result["empirical_q"] < 0.05)
    )
    return result


def arms_disagree_on_every_unit(transport: pd.DataFrame, scoring_arm: str) -> tuple[bool, int]:
    adjusted = transport[
        (transport["derivation_model"] == "batch_adjusted")
        & (transport["scoring_arm"] == scoring_arm)
        & (transport["status"] == "TESTED")
    ]
    wide = adjusted.pivot(
        index=["tissue", "cell_type"], columns="derivation_cohort", values="composite_effect"
    ).dropna()
    if not {"Yazar", "Gao"}.issubset(wide.columns) or wide.empty:
        return False, 0
    disagreement = np.sign(wide["Yazar"]) * np.sign(wide["Gao"]) < 0
    return bool(disagreement.all()), len(wide)


def arm_verdict(
    cohort: str,
    scoring_arm: str,
    counts: dict,
    transport: pd.DataFrame,
    every_unit_disagrees: bool,
) -> str:
    if counts["coverage"] < 0.50 or counts["up"] == 0 or counts["down"] == 0:
        return "INCONCLUSIVE"
    arm = transport[
        (transport["derivation_cohort"] == cohort)
        & (transport["derivation_model"] == "batch_adjusted")
        & (transport["scoring_arm"] == scoring_arm)
    ]
    if arm["passes_both"].any():
        return "PROVES"
    if every_unit_disagrees or not (arm["status"] == "TESTED").any():
        return "INCONCLUSIVE"
    return "DOES NOT PROVE"


def verify_serialized(
    genes: pd.DataFrame,
    transport: pd.DataFrame,
    audits: dict[str, dict],
) -> dict[str, float | int]:
    discrepancies = {"composite": 0.0, "empirical": 0.0}
    group_columns = ["derivation_cohort", "derivation_model", "scoring_arm"]
    for _, indices in transport.groupby(group_columns, sort=True).groups.items():
        family = transport.loc[indices]
        tested = family[family["status"] == "TESTED"]
        if tested.empty:
            continue
        for test_type, p_column, q_column in (
            ("composite", "composite_pvalue", "composite_q"),
            ("empirical", "empirical_pvalue", "empirical_q"),
        ):
            recalculated = multipletests(tested[p_column].to_numpy(dtype=float), method="fdr_bh")[1]
            discrepancy = float(np.max(np.abs(recalculated - tested[q_column].to_numpy(dtype=float))))
            discrepancies[test_type] = max(discrepancies[test_type], discrepancy)
        if not (tested["family_composite_size"] == len(tested)).all():
            raise AssertionError("Composite family size mismatch")
        if not (tested["family_empirical_size"] == len(tested)).all():
            raise AssertionError("Empirical family size mismatch")
    if max(discrepancies.values()) > 1e-12:
        raise AssertionError(f"BH verification failed: {discrepancies}")
    tested = transport[transport["status"] == "TESTED"]
    if not (tested["n_random_sets"] == 2000).all() or not (tested["random_seed"] == 20261105).all():
        raise AssertionError("Matched-null iterations or seed mismatch")
    for row in tested.itertuples():
        ages = {float(value) for value in str(row.fitted_mouse_ages_months).split(";")}
        if min(ages) < 18.0:
            raise AssertionError("A tested mouse donor is under 18 months")
        if row.scoring_arm == PAIRWISE_ARM:
            if not ages.issubset({18.0, 30.0}) or ages != {18.0, 30.0}:
                raise AssertionError("Pairwise arm includes an age other than 18 or 30 months")
            if row.n_donors < 6 or row.n_donors_18m < 3 or row.n_donors_30m < 3:
                raise AssertionError("Pairwise donor floor failed for a tested row")
    for audit in audits.values():
        if audit["min_fitted_age_years"] < 30.0:
            raise AssertionError("A fitted human donor is under 30")
    if len(genes) and genes["derivation_min_donor_age_years"].min() < 30.0:
        raise AssertionError("Signature CSV contains an under-30 derivation minimum")
    pair = tested[tested["scoring_arm"] == PAIRWISE_ARM]
    max_orientation_error = float(pair["binary_orientation_abs_error"].max()) if len(pair) else np.nan
    if len(pair) and max_orientation_error > 1e-10:
        raise AssertionError("Pairwise orientation verification failed")
    return {
        "max_composite_q_discrepancy": discrepancies["composite"],
        "max_empirical_q_discrepancy": discrepancies["empirical"],
        "n_tested_null_verified": len(tested),
        "min_mouse_age_months": float(tested["min_donor_age_months"].min()),
        "max_binary_orientation_error": max_orientation_error,
        "n_pairwise_tested": len(pair),
        "n_pairwise_adjusted_sign_agreements": int(pair["binary_adjusted_sign_matches_raw_difference"].sum()),
    }


def counts_from_genes(genes: pd.DataFrame, key: tuple[str, str]) -> dict[str, float | int]:
    subset = genes[
        (genes["derivation_cohort"] == key[0]) & (genes["derivation_model"] == key[1])
    ] if len(genes) else genes
    clean = subset[subset["passes_clean_filter"]] if len(subset) else subset
    mapped = subset[subset["in_transport_signature"]] if len(subset) else subset
    return {
        "discovery": len(subset),
        "clean": len(clean),
        "mapped": len(mapped),
        "lost": len(clean) - len(mapped),
        "coverage": len(mapped) / len(clean) if len(clean) else 0.0,
        "up": int((mapped["direction"] == "up").sum()) if len(mapped) else 0,
        "down": int((mapped["direction"] == "down").sum()) if len(mapped) else 0,
    }


def fmt(value: float, digits: int = 6) -> str:
    return f"{value:.{digits}g}" if pd.notna(value) else "NA"


def exact(value: float) -> str:
    return f"{value:.17g}" if pd.notna(value) else "NA"


def build_report(
    genes: pd.DataFrame,
    transport: pd.DataFrame,
    counts: dict[tuple[str, str], dict],
    audits: dict[str, dict],
    verification: dict,
    runtime: RuntimeLog,
) -> tuple[str, dict, dict]:
    csv_counts = {key: counts_from_genes(genes, key) for key in counts}
    for key in counts:
        for field in ("discovery", "clean", "mapped", "lost", "up", "down"):
            if int(counts[key][field]) != int(csv_counts[key][field]):
                raise AssertionError(f"Report count mismatch for {key} {field}")
    disagreements = {arm: arms_disagree_on_every_unit(transport, arm) for arm in SCORING_ARMS}
    verdicts = {
        (cohort, arm): arm_verdict(
            cohort, arm, csv_counts[(cohort, "batch_adjusted")], transport, disagreements[arm][0]
        )
        for cohort in ("Yazar", "Gao")
        for arm in SCORING_ARMS
    }
    arm_labels = {PRIMARY_ARM: "18-plus slope", PAIRWISE_ARM: "18-vs-30"}
    lines = [
        "# Mature Block 4: Reverse signature transport",
        "",
        "## Spec",
        "",
        "This is a mature-to-old rerun of the frozen 11-point Block 4. Human derivation uses only donors aged at least 30 years before every model fit. Mouse scoring uses only donors aged at least 18 months. Both derivation and tested contrasts therefore satisfy the mature-to-old rule.",
        "",
        "The original derivation is otherwise unchanged. For each human cell type and donor with at least 20 cells, full-depth raw counts received the human technical gene mask, were converted to per-cell log2(CPM + 1), and were averaged within donor and gene. The cached donor means were subset to 30-plus donors before OLS. Models use age per decade plus sex plus batch. Yazar uses the 75 pool fixed-effect candidates and Gao uses the 7 study indicators, with the original full-rank and residual-df handling. At least 6 donors are required per slope. Discovery requires BH FDR below 0.05 within each derivation cell type and absolute log2 fold change at least 0.25 per decade. A gene discovered in multiple cell types appears once with direction from the sign of its mean discovery slope.",
        "",
        "The P1 clean filter retains protein-coding genes and excludes GM-prefixed, -PS, RIK, ribosomal, mitochondrial, and MRP genes. Mapping uses the Ensembl release 106 one-to-one table in the reverse human-to-mouse direction. No original cache was overwritten.",
        "",
        "Two mouse scoring arms are prespecified for every derivation signature. The primary 18-plus arm fits continuous age in months plus sex plus log library depth over 18, 21, 24, and 30 months where present. The secondary 18-vs-30 arm uses only 18-month and 30-month donors and fits a binary covariate coded 0 at 18 months and 1 at 30 months, plus sex and log library depth. The coefficient is the adjusted 30-minus-18 contrast. Both arms require at least 20 cells per donor, at least 6 donors, exclusion of predicted doublets and unknown cell types, and the original expression and design floors. The pairwise arm additionally requires at least 3 donors on each side after the cell floor.",
        "",
        "Composite scores are mean up-gene expression minus mean down-gene expression. Both arms reuse run_family_c_p1_signature and its matched random-set null unchanged: 2000 draws, direction-specific matching on gene count and donor-mean expression decile, seed 20261105, and a two-sided empirical p centered on the null mean. MDE80 is (t_0.975 + t_0.80) times composite SE using residual df.",
        "",
        "BH correction is separate within each derivation cohort, derivation model, scoring arm, and test type, across tested TMS units only. A unit passes only when composite and empirical q are both below 0.05. The unadjusted sensitivity removes pool or study indicators during derivation, repeats both scoring arms, and remains in separate families.",
        "",
        "The binary orientation was checked independently for every tested pairwise unit. An intercept-plus-age_group OLS with no covariates was fit to donor composite scores, and its coefficient was required to equal the raw mean score at 30 months minus the raw mean at 18 months to absolute tolerance 1e-10. The adjusted coefficient is also compared with that raw direction. Covariate adjustment can legitimately reverse a sign, so the exact equality assertion applies to the unadjusted orientation check.",
        "",
        "All 30-month mice are male, as verified in the mature P1m derivation. Age and sex are not completely aliased because the 18-month group contains both sexes, but the 18-vs-30 contrast retains this age-sex confound and limited sex generalizability.",
        "",
        "## Cache route and human donor coverage",
        "",
        "The existing per-donor caches under tmp_organs/block4_reverse_transport_cache were reused. Both contained age_years and per-donor expression, so no h5ad re-streaming and no new cache directory were needed. Cached donor counts below are unique 30-plus donors with an eligible cached donor-cell-type mean. Fitted donors are the unique donors in cell types retaining at least 6 donors before model fitting.",
        "",
        "| Cohort | 30-plus cached donors | Fitted donors | Cell types cached | Cell types fitted | Minimum fitted age | Maximum fitted age |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for cohort in ("Yazar", "Gao"):
        audit = audits[cohort]
        lines.append(
            f"| {cohort} | {audit['n_donors_30plus_cached']} | {audit['n_donors_fitted']} | "
            f"{audit['n_cell_types_cached']} | {audit['n_cell_types_fitted']} | "
            f"{fmt(audit['min_fitted_age_years'])} | {fmt(audit['max_fitted_age_years'])} |"
        )
    lines.extend([
        "",
        "## Derivation coverage and family sizes",
        "",
        "| Cohort | Model | Scoring arm | Discovery | Clean | Mapped | Lost | Coverage | Up | Down | Composite family | Empirical family | Passing both |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ])
    for key in sorted(counts):
        c = csv_counts[key]
        for scoring_arm in SCORING_ARMS:
            arm = transport[
                (transport["derivation_cohort"] == key[0])
                & (transport["derivation_model"] == key[1])
                & (transport["scoring_arm"] == scoring_arm)
            ]
            tested = arm[arm["status"] == "TESTED"]
            lines.append(
                f"| {key[0]} | {key[1]} | {arm_labels[scoring_arm]} | {c['discovery']} | {c['clean']} | "
                f"{c['mapped']} | {c['lost']} | {100*c['coverage']:.1f}% | {c['up']} | {c['down']} | "
                f"{len(tested)} | {len(tested)} | {int(tested['passes_both'].sum())} |"
            )
    lines.extend([
        "",
        "Discovery counts are unique genes after pooling cell-type discovery rows. Clean counts retain the P1 filter. Mapped counts and coverage are relative to the clean set. Per-gene mean slopes, minimum discovery p and q values, biotypes, directions, and mapping statuses are in the signature-genes CSV.",
        "",
        "## Asymmetry",
        "",
        "Human batch noise enters the DERIVATION in this block, so a reverse null is weaker evidence than the mouse-to-human null, where derivation came from one controlled experiment. A reverse null is therefore not equally strong evidence of symmetric transport failure.",
        "",
        "## Verdict rules and sign disagreement",
        "",
        "PROVES requires at least one TMS unit passing both corrected tests in that batch-adjusted cohort and scoring arm. DOES NOT PROVE requires no passing unit in that arm and per-unit MDE80 reporting. INCONCLUSIVE applies when under 50 percent of the clean set maps, when a mapped signature lacks an up or down direction, when cross-cohort adjusted slopes disagree in sign on every commonly tested unit, or when no unit is testable for another reason. Missing-direction units are UNTESTABLE, not nulls.",
        "",
        "The sign-disagreement rule was evaluated before verdict calculation separately in each scoring arm.",
        "",
        "| Scoring arm | Commonly tested units | Opposite signs on every shared unit |",
        "|---|---:|---|",
    ])
    for scoring_arm in SCORING_ARMS:
        every, shared = disagreements[scoring_arm]
        lines.append(f"| {arm_labels[scoring_arm]} | {shared} | {every} |")
    lines.extend([
        "",
        "| Derivation cohort | Scoring arm | Verdict | Tested units | Passing both |",
        "|---|---|---|---:|---:|",
    ])
    for cohort in ("Yazar", "Gao"):
        for scoring_arm in SCORING_ARMS:
            arm = transport[
                (transport["derivation_cohort"] == cohort)
                & (transport["derivation_model"] == "batch_adjusted")
                & (transport["scoring_arm"] == scoring_arm)
            ]
            tested = arm[arm["status"] == "TESTED"]
            lines.append(
                f"| {cohort} | {arm_labels[scoring_arm]} | {verdicts[(cohort, scoring_arm)]} | "
                f"{len(tested)} | {int(tested['passes_both'].sum())} |"
            )
    for cohort in ("Yazar", "Gao"):
        for scoring_arm in SCORING_ARMS:
            arm = transport[
                (transport["derivation_cohort"] == cohort)
                & (transport["derivation_model"] == "batch_adjusted")
                & (transport["scoring_arm"] == scoring_arm)
            ].copy()
            tested = arm[arm["status"] == "TESTED"]
            effect_name = "Slope per 12 months" if scoring_arm == PRIMARY_ARM else "30-minus-18 effect"
            mde_name = "MDE80 per 12 months" if scoring_arm == PRIMARY_ARM else "MDE80 30-minus-18"
            lines.extend([
                "",
                f"### {cohort} batch-adjusted {arm_labels[scoring_arm]} results",
                "",
                f"Passing both tests: {int(tested['passes_both'].sum())} of {len(tested)} tested units. Composite and empirical family sizes are both {len(tested)}.",
                "",
                f"| Tissue | Cell type | N donors | N genes | {effect_name} | SE | CI95 low | CI95 high | Composite p | Composite q | Empirical p | Empirical q | {mde_name} | Status |",
                "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|",
            ])
            for row in arm.sort_values(["tissue", "cell_type"]).itertuples():
                if row.status == "TESTED":
                    if scoring_arm == PRIMARY_ARM:
                        effect = row.composite_slope_per_12_months
                        se = row.composite_se_per_12_months
                        low = row.composite_ci95_low_per_12_months
                        high = row.composite_ci95_high_per_12_months
                        mde = row.mde80_per_12_months
                    else:
                        effect = row.composite_slope_30_minus_18
                        se = row.composite_se
                        low = row.composite_ci95_low
                        high = row.composite_ci95_high
                        mde = row.mde80_30_minus_18
                    lines.append(
                        f"| {row.tissue} | {row.cell_type} | {row.n_donors} | {row.n_genes_mapped} | "
                        f"{fmt(effect)} | {fmt(se)} | {fmt(low)} | {fmt(high)} | {fmt(row.composite_pvalue)} | "
                        f"{exact(row.composite_q)} | {fmt(row.empirical_pvalue)} | {exact(row.empirical_q)} | "
                        f"{fmt(mde)} | {row.status} |"
                    )
                else:
                    lines.append(
                        f"| {row.tissue} | {row.cell_type} | {row.n_donors} | 0 | NA | NA | NA | NA | "
                        f"NA | NA | NA | NA | NA | {row.status} |"
                    )
    lines.extend([
        "",
        "## MDE80 by arm",
        "",
        "| Cohort | Model | Scoring arm | Tested units | Median MDE80 | Minimum MDE80 | Maximum MDE80 | Units |",
        "|---|---|---|---:|---:|---:|---:|---|",
    ])
    for key in sorted(counts):
        for scoring_arm in SCORING_ARMS:
            mask = (
                (transport["derivation_cohort"] == key[0])
                & (transport["derivation_model"] == key[1])
                & (transport["scoring_arm"] == scoring_arm)
                & (transport["status"] == "TESTED")
            )
            tested = transport[mask]
            units = "log2CPM per month" if scoring_arm == PRIMARY_ARM else "log2CPM for 30 minus 18"
            lines.append(
                f"| {key[0]} | {key[1]} | {arm_labels[scoring_arm]} | {len(tested)} | "
                f"{fmt(tested['mde80'].median())} | {fmt(tested['mde80'].min())} | "
                f"{fmt(tested['mde80'].max())} | {units} |"
            )
    lines.extend([
        "",
        "These MDE80 values are the smallest composite effects with 80 percent power at the two-sided 0.05 level under the fitted t model. Large values show reduced resolution after age restriction. Non-passage is not equivalence or proof of no mature effect.",
        "",
        "## Unadjusted sensitivity",
        "",
        "| Cohort | Scoring arm | Tested units | Passing both | Minimum composite q | Minimum empirical q | Median MDE80 | MDE80 range |",
        "|---|---|---:|---:|---:|---:|---:|---|",
    ])
    for cohort in ("Yazar", "Gao"):
        for scoring_arm in SCORING_ARMS:
            arm = transport[
                (transport["derivation_cohort"] == cohort)
                & (transport["derivation_model"] == "unadjusted")
                & (transport["scoring_arm"] == scoring_arm)
            ]
            tested = arm[arm["status"] == "TESTED"]
            mde_range = (
                f"{fmt(tested['mde80'].min())} to {fmt(tested['mde80'].max())}" if len(tested) else "NA"
            )
            lines.append(
                f"| {cohort} | {arm_labels[scoring_arm]} | {len(tested)} | {int(tested['passes_both'].sum())} | "
                f"{exact(tested['composite_q'].min())} | {exact(tested['empirical_q'].min())} | "
                f"{fmt(tested['mde80'].median())} | {mde_range} |"
            )
    yazar_comparisons = []
    for scoring_arm in SCORING_ARMS:
        verdict = verdicts[("Yazar", scoring_arm)]
        if verdict == "DOES NOT PROVE":
            comparison = "SURVIVED"
            explanation = "No mature unit passes both corrected tests."
        elif verdict == "INCONCLUSIVE":
            comparison = "UNTESTABLE"
            explanation = "The mature arm cannot confirm the original non-passing result under the frozen verdict rules."
        else:
            comparison = "DIED"
            explanation = "At least one mature unit passes both corrected tests."
        yazar_comparisons.append((scoring_arm, comparison, explanation))
    gao_counts = csv_counts[("Gao", "batch_adjusted")]
    gao_missing = gao_counts["up"] == 0 or gao_counts["down"] == 0 or gao_counts["coverage"] < 0.50
    lines.extend([
        "",
        "## Comparison with original Block 4",
        "",
        "The original Yazar batch-adjusted derivation gave 14 discovered genes, 13 clean genes, 10 mapped genes, 6 up and 4 down, with 0 of 51 TMS units passing both corrected tests and a DOES NOT PROVE verdict. The original Gao batch-adjusted derivation gave only 1 down gene and no up gene, so that arm was UNTESTABLE_MISSING_UP_DIRECTION and INCONCLUSIVE, never a null.",
        "",
        "| Original result | Mature scoring arm | Comparison | Evidence |",
        "|---|---|---|---|",
    ])
    for scoring_arm, comparison, explanation in yazar_comparisons:
        tested = transport[
            (transport["derivation_cohort"] == "Yazar")
            & (transport["derivation_model"] == "batch_adjusted")
            & (transport["scoring_arm"] == scoring_arm)
            & (transport["status"] == "TESTED")
        ]
        lines.append(
            f"| Original Yazar 0 of 51 both-pass | {arm_labels[scoring_arm]} | {comparison} | "
            f"Mature result is {int(tested['passes_both'].sum())} of {len(tested)} both-pass; {explanation} |"
        )
    for scoring_arm in SCORING_ARMS:
        if gao_missing:
            comparison = "UNTESTABLE"
            explanation = (
                f"Mature Gao still lacks a fully estimable mapped bidirectional signature: {gao_counts['up']} up, "
                f"{gao_counts['down']} down, coverage {100*gao_counts['coverage']:.1f}%."
            )
        else:
            comparison = "DIED"
            explanation = (
                f"The original missing-up inestimability is removed: mature Gao has {gao_counts['up']} up and "
                f"{gao_counts['down']} down mapped genes. Transport verdict is {verdicts[('Gao', scoring_arm)]}."
            )
        lines.append(
            f"| Original Gao 1 down, 0 up, inestimable | {arm_labels[scoring_arm]} | {comparison} | {explanation} |"
        )
    lines.extend([
        "",
        "Mature family sizes differ from the original 51-unit family, so the comparison concerns the verdict and passage rule, not identical denominators or identical donor sets.",
        "",
        "## Verification",
        "",
        f"BH q-values were recomputed from raw p-values in the serialized transport CSV with statsmodels multipletests. Maximum absolute discrepancy: composite {exact(verification['max_composite_q_discrepancy'])}; empirical {exact(verification['max_empirical_q_discrepancy'])}. Passing counts and family sizes in this report were recomputed from that CSV.",
        "",
        f"Minimum fitted human ages: Yazar {fmt(audits['Yazar']['min_fitted_age_years'])} years; Gao {fmt(audits['Gao']['min_fitted_age_years'])} years. Minimum fitted mouse age across TESTED rows: {fmt(verification['min_mouse_age_months'])} months. Pairwise TESTED rows contain exactly ages 18 and 30 months and satisfy at least 6 total donors with at least 3 per side.",
        "",
        f"All {verification['n_tested_null_verified']} TESTED rows have n_random_sets = 2000 and random_seed = 20261105. The pairwise orientation check covered {verification['n_pairwise_tested']} units with maximum absolute unadjusted-coefficient error {exact(verification['max_binary_orientation_error'])}. Adjusted and raw mean-difference signs agreed in {verification['n_pairwise_adjusted_sign_agreements']} of {verification['n_pairwise_tested']} pairwise units.",
        "",
        "Forbidden Unicode em dash and en dash characters are rejected by the CSV and Markdown writers. A final source and output character scan is performed after the run.",
        "",
        "## Runtime",
        "",
        f"Model: {MODEL_NAME}. All matched nulls use {ITERATIONS} iterations and seed {NULL_SEED}. Exact stage start and end timestamps are in the runtime CSV.",
        "",
        "| Stage | Seconds |",
        "|---|---:|",
    ])
    for row in runtime.rows:
        lines.append(f"| {row['stage']} | {row['duration_seconds']:.3f} |")
    lines.extend([
        "",
        "## Files",
        "",
        "- scripts/run_mature_block4_reverse_transport.py",
        "- results/followup-mature-block4-reverse-signature-genes.csv",
        "- results/followup-mature-block4-reverse-transport.csv",
        "- results/followup-mature-block4-runtime.csv",
        "- reports/followup-mature-block4-reverse-transport.md",
        "",
    ])
    return "\n".join(lines), verdicts, disagreements


def main() -> None:
    runtime = RuntimeLog()
    wait_for_memory("Mature Block 4")
    for path in (YAZAR_H5AD, GAO_H5AD, REPO / "data" / "reference" / "orthologs-mouse-human.csv"):
        if not path.exists():
            raise FileNotFoundError(path)
    orthologs = pd.read_csv(REPO / "data" / "reference" / "orthologs-mouse-human.csv")
    orthologs = orthologs[orthologs["orthology_type"] == "ortholog_one2one"]
    orthologs = orthologs.drop_duplicates("human_gene_symbol", keep="last")
    with runtime.stage("load_human_biotypes"):
        biotypes = human_biotypes()

    all_genes = []
    programs = {}
    counts = {}
    audits = {}
    for cohort, loader in (("Yazar", derive_yazar), ("Gao", derive_gao)):
        with runtime.stage(f"load_and_restrict_{cohort}_cache"):
            print(f"Loading {cohort} cached derivation cohort", flush=True)
            expr, meta, model_covariates = loader()
            audits[cohort] = restrict_human_to_30plus(expr, meta, cohort)
        for model, covariates in model_covariates.items():
            with runtime.stage(f"derive_{cohort}_{model}"):
                print(f"Deriving {cohort} {model} mature signature", flush=True)
                tests = fit_gene_models(expr, meta, tuple(covariates))
                gene_frame, program, count = derive_signature(
                    tests, cohort, model, orthologs, biotypes
                )
                all_genes.append(gene_frame)
                programs[(cohort, model)] = program
                counts[(cohort, model)] = count
        del expr, meta
        gc.collect()

    genes = pd.concat(all_genes, ignore_index=True)
    if len(genes):
        genes["derivation_min_donor_age_years"] = genes["derivation_cohort"].map(
            {cohort: audit["min_fitted_age_years"] for cohort, audit in audits.items()}
        )
        genes["derivation_n_donors_fitted"] = genes["derivation_cohort"].map(
            {cohort: audit["n_donors_fitted"] for cohort, audit in audits.items()}
        )
    with runtime.stage("prepare_18plus_mouse_targets"):
        targets = mouse_targets()
    transport = score_targets(targets, programs, runtime)
    del targets
    gc.collect()
    with runtime.stage("write_gene_and_transport_csv"):
        write_clean_csv(genes, GENE_OUT)
        write_clean_csv(transport, TRANSPORT_OUT)
    with runtime.stage("verify_serialized_statistics_and_ages"):
        serialized_genes = pd.read_csv(GENE_OUT)
        serialized_transport = pd.read_csv(TRANSPORT_OUT)
        verification = verify_serialized(serialized_genes, serialized_transport, audits)
    with runtime.stage("build_report_from_serialized_csv"):
        report, verdicts, disagreements = build_report(
            serialized_genes, serialized_transport, counts, audits, verification, runtime
        )
    with runtime.stage("write_report"):
        write_clean_markdown(report, REPORT_OUT)
    with runtime.stage("write_runtime_csv"):
        runtime_frame = runtime.finish()
        write_clean_csv(runtime_frame, RUNTIME_OUT)
    print(f"Wrote {GENE_OUT} ({len(genes)} rows)", flush=True)
    print(f"Wrote {TRANSPORT_OUT} ({len(transport)} rows)", flush=True)
    print(f"Wrote {RUNTIME_OUT} ({len(runtime_frame)} rows)", flush=True)
    print(f"Wrote {REPORT_OUT}; verdicts {verdicts}; sign checks {disagreements}", flush=True)


if __name__ == "__main__":
    main()
