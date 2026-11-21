#!/usr/bin/env python3
"""Finalize mature mouse and human Hellinger drift and their interaction."""

from __future__ import annotations

import argparse
import gc
import os
import sys
import time
import warnings
from datetime import datetime, timezone
from pathlib import Path

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, sparse_dataset

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_cross_species_expression_vector import technical_gene_mask
from scripts.run_rev2_hellinger import (
    fit_global2,
    hellinger_distance,
    load_tissue_arrays,
    mouse_tissue_endpoints,
    mouse_tissue_names,
)

RESULTS = REPO / "results"
SCRATCH = REPO / "scratch" / "mature-drift"
DATA_ROOT = Path("/home/reza/projects/agescape-revision-data")
HUMAN_H5ADS = DATA_ROOT / "human-h5ads"
MODEL_USED = "kimi-code/k3.1-pro"
RUNTIME_CSV = RESULTS / "followup-mature-drift-runtime.csv"
MOUSE_CSV = RESULTS / "followup-mature-drift-mouse.csv"
HUMAN_CSV = RESULTS / "followup-mature-drift-human.csv"
INTERACTION_CSV = RESULTS / "followup-mature-drift-interaction.csv"
REUSED_CSV = RESULTS / "followup-mature-drift-reused.csv"
MOUSE_DONORS_CSV = RESULTS / "followup-mature-drift-mouse-donors.csv"
HUMAN_DONORS_CSV = RESULTS / "followup-mature-drift-human-donors.csv"
AGE_RE = r"(\d+)-year-old"
HUMAN_COHORTS = {
    "MSSM_DLPFC": ("37a17b78-4864-4a42-b67b-31c00962795a", "brain"),
    "HBCC_DLPFC": ("5c97eeeb-7e52-44b3-b010-b832b1f5424c", "brain"),
    "Kidney_LivingDonor": ("0bae7ebf-eb54-46a6-be9a-3461cecefa4c", "kidney"),
    "Kidney_SCR_v15": ("dea717d4-7bc0-4e46-950f-fd7e1cc8df7d", "kidney"),
    "Lung_CellRef": ("eb499fd8-7000-419f-8854-926b9b61b11f", "lung"),
    "Lung_Integrated": ("9f222629-9e39-47d0-b83f-e08d610c7479", "lung"),
    "ImmuneAging_BAL": ("1b350d0a-4535-4879-beb6-1142f3f94947", "lung"),
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def record_runtime(pipeline: str, start_utc: str, started: float, status: str, details: str = "") -> None:
    row = pd.DataFrame([{
        "pipeline": pipeline,
        "start_utc": start_utc,
        "end_utc": utc_now(),
        "duration_seconds": time.monotonic() - started,
        "model_used": MODEL_USED,
        "status": status,
        "details": details,
    }])
    if RUNTIME_CSV.exists():
        previous = pd.read_csv(RUNTIME_CSV)
        row = pd.concat([previous, row], ignore_index=True)
    row.to_csv(RUNTIME_CSV, index=False)


def run_stage(name: str, function) -> None:
    started = time.monotonic()
    start_utc = utc_now()
    print(f"{name} started {start_utc}; model {MODEL_USED}", flush=True)
    try:
        function()
    except Exception as exc:
        record_runtime(name, start_utc, started, "FAILED", f"{type(exc).__name__}: {exc}")
        raise
    record_runtime(name, start_utc, started, "completed")
    print(f"{name} completed after {time.monotonic() - started:.3f} seconds", flush=True)


def untestable_row(species: str, variant: str, endpoint: str, adjustment: str,
                   depth: str, tissue: str = "All tested tissues", cohort: str = "All cohorts",
                   reason: str = "", n_donors: int = 0, n_units: int = 0,
                   n_donor_units: int = 0, **extra) -> dict:
    return {
        "status": "UNTESTABLE",
        "provenance": "new mature restriction",
        "species": species,
        "variant": variant,
        "endpoint": endpoint,
        "adjustment": adjustment,
        "depth": depth,
        "reference": "18m_cross_fitted" if species == "Mouse" else "human_30_39_cross_fitted",
        "tissue": tissue,
        "cohort": cohort,
        "age_scale_years": 0.4 if species == "Mouse" else 12.25,
        "n_donors": int(n_donors),
        "n_units": int(n_units),
        "n_donor_units": int(n_donor_units),
        "slope": np.nan,
        "se": np.nan,
        "ci95_low": np.nan,
        "ci95_high": np.nan,
        "pvalue": np.nan,
        "notes": reason,
        **extra,
    }


def fitted_row(fit: dict, status: str, provenance: str, **extra) -> dict:
    row = dict(fit)
    row["slope"] = row.pop("effect_sd")
    row["status"] = status
    row["provenance"] = provenance
    row.update(extra)
    return row


def mature_full_mouse() -> pd.DataFrame:
    donors = pd.read_csv(RESULTS / "reviewer-rev2-hellinger-donors.csv")
    part = donors[(donors["table"] == "mouse_ref_18m") & (donors["age"] >= 1.5 - 1e-8)].copy()
    assert len(part) and part["age"].min() >= 1.5
    assert part["donor"].nunique() == 14
    return part


def reused_frequentist_rows(full: pd.DataFrame) -> list[dict]:
    source = RESULTS / "rev3-bayesian-vector-mature-only-frequentist.csv"
    frame = pd.read_csv(source)
    rows = []
    for endpoint in ("profile_hellinger_centroid_mean", "pseudobulk_hellinger"):
        for adjustment, model in (("unadjusted", "unadjusted"), ("linear_log_umi", "depth_adjusted")):
            spec = f"mature_ge_18m-{endpoint}-{model}-prior_0.50"
            match = frame[frame["spec_id"] == spec]
            assert len(match) == 1
            item = match.iloc[0]
            rows.append({
                "status": "REUSED",
                "provenance": str(source.relative_to(REPO)),
                "species": "Mouse",
                "variant": "global",
                "endpoint": endpoint,
                "adjustment": adjustment,
                "depth": "full",
                "reference": "18m_cross_fitted",
                "tissue": "All tested tissues",
                "cohort": "All cohorts",
                "age_scale_years": 0.4,
                "n_donors": int(item["n_donors"]),
                "n_units": int(item["n_units"]),
                "n_donor_units": int(item["n_obs"]),
                "slope": float(item["freq_effect_10pct_lifespan"]),
                "se": float(item["freq_se"]),
                "ci95_low": float(item["freq_ci95_low"]),
                "ci95_high": float(item["freq_ci95_high"]),
                "pvalue": float(item["freq_pvalue"]),
                "notes": f"Reused as-is from {spec}; no full-depth global refit.",
            })
    return rows


def fit_mouse_variants(part: pd.DataFrame, depth: str, rows: list[dict],
                       include_influence: bool = True) -> None:
    provenance = "new fixed-depth 18m reference" if depth != "full" else "new model on filtered committed 18m donor endpoints"
    endpoints = ("profile_hellinger_centroid_mean", "pseudobulk_hellinger")
    for endpoint in endpoints:
        for adjustment in ("unadjusted", "linear_log_umi", "inv_sqrt_umi", "spline_log_umi"):
            if depth == "full" and adjustment in ("unadjusted", "linear_log_umi"):
                continue
            if depth == "full" and adjustment == "inv_sqrt_umi":
                continue
            fit = fit_global2(
                part, endpoint, "Mouse", "All tested tissues", "All cohorts", adjustment,
                0.4, variant="global", depth_label=depth, reference="18m_cross_fitted",
                endpoint=endpoint, source="mature_restricted",
            )
            if fit is None:
                rows.append(untestable_row(
                    "Mouse", "global", endpoint, adjustment, depth,
                    reason="No estimable fit after unit, donor, and rank floors.",
                    n_donors=part["donor"].nunique(), n_donor_units=len(part),
                ))
            else:
                rows.append(fitted_row(fit, "NEW", provenance))
    for tissue in mouse_tissue_names():
        tissue_part = part[part["tissue"] == tissue]
        fit = fit_global2(
            tissue_part, "profile_hellinger_centroid_mean", "Mouse", tissue, "All cohorts",
            "linear_log_umi", 0.4, variant="per_tissue", depth_label=depth,
            reference="18m_cross_fitted", endpoint="profile_hellinger_centroid_mean",
            source="mature_restricted",
        ) if len(tissue_part) else None
        if fit is None:
            rows.append(untestable_row(
                "Mouse", "per_tissue", "profile_hellinger_centroid_mean", "linear_log_umi",
                depth, tissue=tissue,
                reason="No estimable per-organ fit after mature unit, donor, and rank floors.",
                n_donors=tissue_part["donor"].nunique(), n_units=tissue_part["cell_type"].nunique(),
                n_donor_units=len(tissue_part),
            ))
        else:
            rows.append(fitted_row(fit, "NEW", provenance))
    if not include_influence:
        return
    influence_specs = [("profile_hellinger_centroid_mean", "linear_log_umi")]
    if depth == "full":
        influence_specs = [
            ("profile_hellinger_centroid_mean", "unadjusted"),
            ("profile_hellinger_centroid_mean", "linear_log_umi"),
            ("pseudobulk_hellinger", "unadjusted"),
            ("pseudobulk_hellinger", "linear_log_umi"),
        ]
    for endpoint, adjustment in influence_specs:
        for omitted in sorted(part["donor"].unique()):
            remaining = part[part["donor"] != omitted]
            fit = fit_global2(
                remaining, endpoint, "Mouse", "All tested tissues", "All cohorts", adjustment,
                0.4, variant="leave_one_donor_out", depth_label=depth,
                reference="18m_cross_fitted", endpoint=endpoint, source="mature_restricted",
            )
            if fit is not None:
                rows.append(fitted_row(
                    fit, "NEW", provenance, omitted_donor=omitted,
                    notes="Regression influence refit; cross-fitted reference endpoints are held fixed.",
                ))


def mouse_pipeline() -> None:
    full = mature_full_mouse()
    rows = reused_frequentist_rows(full)
    fit_mouse_variants(full, "full", rows)
    donor_outputs = [full.assign(depth="full", mature_input_source="committed_mouse_ref_18m")]
    for depth in ("500", "1500", "3000"):
        all_rows = []
        for tissue in mouse_tissue_names():
            X, obs, var_names = load_tissue_arrays(tissue)
            mature = obs["age"].to_numpy(dtype=float) >= 1.5 - 1e-8
            if mature.any():
                computed = mouse_tissue_endpoints(
                    X[mature], obs.loc[mature], var_names, tissue, int(depth), "18m"
                )
                all_rows.extend(computed)
            del X, obs, var_names
            gc.collect()
            print(f"Mouse depth {depth} {tissue}: cumulative donor units {len(all_rows)}", flush=True)
        part = pd.DataFrame(all_rows)
        assert len(part) and part["age"].min() >= 1.5
        part["table"] = f"mouse_depth_{depth}_ref_18m_mature"
        donor_outputs.append(part.assign(depth=depth, mature_input_source="new_mature_cell_level"))
        fit_mouse_variants(part, depth, rows)
        checkpoint = pd.DataFrame(rows)
        checkpoint.to_csv(SCRATCH / f"mouse-through-{depth}.csv", index=False)
    mouse = pd.DataFrame(rows)
    mouse.to_csv(MOUSE_CSV, index=False)
    pd.concat(donor_outputs, ignore_index=True).to_csv(MOUSE_DONORS_CSV, index=False)
    primary = mouse[(mouse["variant"] == "global") & (mouse["depth"] == "full")]
    print(primary[["endpoint", "adjustment", "status", "slope", "pvalue", "n_donors"]].to_string(index=False), flush=True)


def read_human_obs(handle: h5py.File) -> dict[str, np.ndarray]:
    names = ("donor_id", "development_stage", "disease", "cell_type", "sex", "is_primary_data")
    return {name: np.asarray(read_elem(handle["obs"][name])) for name in names}


def matrix_rows(dataset, rows: np.ndarray, chunk_size: int = 2000):
    for start in range(0, len(rows), chunk_size):
        yield rows[start:start + chunk_size], dataset[rows[start:start + chunk_size]].tocsr()


def add_reference_counts(stats: dict, cell_type: str, donor: str, matrix, n_genes: int) -> None:
    counts = np.asarray(matrix.sum(axis=0)).ravel().astype(np.float64)
    cell_stats = stats.setdefault(cell_type, {"counts": np.zeros(n_genes), "n_cells": 0, "donors": {}})
    donor_stats = cell_stats["donors"].setdefault(
        donor, {"counts": np.zeros(n_genes), "n_cells": 0}
    )
    cell_stats["counts"] += counts
    cell_stats["n_cells"] += matrix.shape[0]
    donor_stats["counts"] += counts
    donor_stats["n_cells"] += matrix.shape[0]


def human_cohort_distances(cohort: str, original: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    SCRATCH.mkdir(parents=True, exist_ok=True)
    dataset_id, tissue = HUMAN_COHORTS[cohort]
    path = HUMAN_H5ADS / f"{dataset_id}.h5ad"
    old = original[(original["cohort"] == cohort) & (original["age"] >= 30)].copy()
    original_units = set(old["cell_type"].astype(str))
    scored_donor_ids = set(old["donor"].astype(str))
    assert original_units and scored_donor_ids
    with h5py.File(path, "r") as handle:
        columns = read_human_obs(handle)
        donor = columns["donor_id"].astype(str)
        cell_type = columns["cell_type"].astype(str)
        stage = pd.Series(columns["development_stage"].astype(str))
        age = stage.str.extract(AGE_RE, expand=False).astype(float).to_numpy()
        primary = columns["is_primary_data"].astype(bool)
        normal = columns["disease"].astype(str) == "normal"
        base = (
            primary & normal & np.isfinite(age) & (age >= 30) & (age <= 120)
            & np.isin(cell_type, list(original_units))
        )
        row_map = {}
        for row in np.flatnonzero(base):
            row_map.setdefault(str(donor[row]), []).append(int(row))
        row_map = {key: np.asarray(value, dtype=np.int64) for key, value in row_map.items()}
        reference_donors = sorted({
            str(donor[row]) for row in np.flatnonzero(base & (age >= 30) & (age <= 39))
        })
        scored_donors = sorted(scored_donor_ids & set(row_map))
        var_names = np.asarray(read_elem(handle["var"]["feature_name"])).astype(str)
        keep_genes = ~technical_gene_mask(var_names)
        matrix_dataset = sparse_dataset(handle["X"])
        stats = {}
        for index, ref_donor in enumerate(reference_donors, start=1):
            rows = row_map[ref_donor]
            for selected_rows, matrix in matrix_rows(matrix_dataset, rows):
                library = np.asarray(matrix.sum(axis=1)).ravel()
                qc = library >= 1500
                if not qc.any():
                    continue
                selected_rows = selected_rows[qc]
                matrix = matrix[qc][:, keep_genes]
                for ct in np.unique(cell_type[selected_rows]):
                    mask = cell_type[selected_rows] == str(ct)
                    add_reference_counts(stats, str(ct), ref_donor, matrix[mask], matrix.shape[1])
            if index % 5 == 0 or index == len(reference_donors):
                print(f"Human {cohort}: reference donors {index}/{len(reference_donors)}", flush=True)
        valid_reference_units = {
            ct for ct, value in stats.items()
            if sum(d["n_cells"] >= 20 for d in value["donors"].values()) >= 2
        }
        output = {}
        for index, scored_donor in enumerate(scored_donors, start=1):
            rows = row_map[scored_donor]
            accumulated = {}
            for selected_rows, matrix in matrix_rows(matrix_dataset, rows):
                library = np.asarray(matrix.sum(axis=1)).ravel()
                qc = library >= 1500
                if not qc.any():
                    continue
                selected_rows = selected_rows[qc]
                matrix = matrix[qc][:, keep_genes].astype(np.float64)
                totals = np.asarray(matrix.sum(axis=1)).ravel()
                if np.any(totals <= 0):
                    raise ValueError(f"{cohort} {scored_donor}: zero library after technical gene removal")
                for ct in np.unique(cell_type[selected_rows]):
                    ct = str(ct)
                    if ct not in valid_reference_units:
                        continue
                    mask = cell_type[selected_rows] == ct
                    donor_matrix = matrix[mask]
                    ref = stats[ct]
                    own = ref["donors"].get(scored_donor)
                    ref_counts = ref["counts"] - (own["counts"] if own else 0)
                    ref_n = ref["n_cells"] - (own["n_cells"] if own else 0)
                    if ref_n < 50 or ref_counts.sum() <= 0:
                        continue
                    distances = hellinger_distance(donor_matrix, ref_counts)
                    item = accumulated.setdefault(ct, {
                        "n_cells": 0, "distance_sum": 0.0, "log_umi_sum": 0.0,
                        "reference_cells": int(ref_n),
                        "reference_donors": sum(d["n_cells"] >= 20 for d in ref["donors"].values()),
                    })
                    item["n_cells"] += donor_matrix.shape[0]
                    item["distance_sum"] += float(distances.sum())
                    item["log_umi_sum"] += float(np.log(totals[mask]).sum())
            donor_age = old.loc[old["donor"].astype(str) == scored_donor, "age"]
            donor_sex = old.loc[old["donor"].astype(str) == scored_donor, "sex_male"]
            assert len(donor_age) and np.allclose(donor_age, donor_age.iloc[0])
            for ct, item in accumulated.items():
                if item["n_cells"] < 20:
                    continue
                key = (scored_donor, ct)
                output[key] = {
                    "species": "Human", "tissue": tissue, "cohort": cohort,
                    "cell_type": ct, "donor": scored_donor,
                    "age": float(donor_age.iloc[0]), "sex_male": float(donor_sex.iloc[0]),
                    "n_cells": item["n_cells"],
                    "mean_log_umi": item["log_umi_sum"] / item["n_cells"],
                    "profile_hellinger": item["distance_sum"] / item["n_cells"],
                    "n_reference_donors_for_unit": item["reference_donors"],
                    "n_reference_cells_after_loo": item["reference_cells"],
                    "reference_window_min": 30.0, "reference_window_max": 39.0,
                    "source_object": str(path),
                }
            if index % 10 == 0 or index == len(scored_donors):
                print(f"Human {cohort}: scored donors {index}/{len(scored_donors)}", flush=True)
        frame = pd.DataFrame(output.values())
        if len(frame):
            assert frame["age"].min() >= 30
            assert frame["n_reference_donors_for_unit"].min() >= 2
            assert frame["n_reference_cells_after_loo"].min() >= 50
        reference_used = {
            donor_id for ct in valid_reference_units for donor_id, value in stats[ct]["donors"].items()
            if value["n_cells"] >= 20
        }
        info = {
            "cohort": cohort,
            "n_source_scored_donors_30plus": len(scored_donors),
            "n_scored_donors": frame["donor"].nunique() if len(frame) else 0,
            "n_reference_donors": len(reference_used),
            "n_reference_units": len(valid_reference_units),
            "n_donor_units": len(frame),
            "age_min": frame["age"].min() if len(frame) else np.nan,
            "age_max": frame["age"].max() if len(frame) else np.nan,
            "source_object": str(path),
        }
    frame.to_csv(SCRATCH / f"human-{cohort}.csv", index=False)
    return frame, info


def add_human_fit(rows: list[dict], part: pd.DataFrame, adjustment: str, variant: str,
                  tissue: str = "All tested tissues", cohort: str = "All cohorts",
                  info: dict | None = None) -> None:
    fit = fit_global2(
        part, "profile_hellinger", "Human", tissue, cohort, adjustment, 12.25,
        variant=variant, depth_label="full", reference="human_30_39_cross_fitted",
        endpoint="profile_hellinger_centroid_mean", source="local_census_mature",
    ) if len(part) else None
    extra = {}
    if info:
        extra = {
            "n_reference_donors": info["n_reference_donors"],
            "n_source_scored_donors_30plus": info["n_source_scored_donors_30plus"],
            "age_min": info["age_min"], "age_max": info["age_max"],
        }
    if fit is None:
        reason = "No estimable mature fit. Check reference and scored donor floors."
        if info and info["n_reference_units"] == 0:
            reason = "No original vector unit has at least 2 reference donors aged 30-39 with at least 20 QC cells."
        rows.append(untestable_row(
            "Human", variant, "profile_hellinger_centroid_mean", adjustment, "full",
            tissue=tissue, cohort=cohort, reason=reason,
            n_donors=part["donor"].nunique() if len(part) else 0,
            n_units=part["cell_type"].nunique() if len(part) else 0,
            n_donor_units=len(part), **extra,
        ))
    else:
        rows.append(fitted_row(
            fit, "NEW", "new local Census 30-39 reference and 30-plus scores", **extra
        ))


def human_pipeline() -> None:
    original = pd.read_csv(RESULTS / "h3-cross-species-vector-donors.csv")
    original = original[original["species"] == "Human"].copy()
    frames = []
    infos = []
    for cohort in HUMAN_COHORTS:
        frame, info = human_cohort_distances(cohort, original)
        frames.append(frame)
        infos.append(info)
        print(f"Human {cohort}: {info}", flush=True)
    donors = pd.concat(frames, ignore_index=True)
    assert donors["age"].min() >= 30
    assert donors["reference_window_min"].min() >= 30
    assert donors["reference_window_max"].max() <= 39
    donors.to_csv(HUMAN_DONORS_CSV, index=False)
    info_by_cohort = {item["cohort"]: item for item in infos}
    rows = []
    for adjustment in ("unadjusted", "linear_log_umi", "inv_sqrt_umi", "linear_plus_inv_sqrt",
                       "quadratic_log_umi", "spline_log_umi"):
        add_human_fit(rows, donors, adjustment, "global")
    for cohort in HUMAN_COHORTS:
        part = donors[donors["cohort"] == cohort]
        for adjustment in ("unadjusted", "linear_log_umi"):
            add_human_fit(rows, part, adjustment, "per_cohort", cohort=cohort,
                          info=info_by_cohort[cohort])
    for tissue, part in donors.groupby("tissue"):
        for adjustment in ("unadjusted", "linear_log_umi", "inv_sqrt_umi", "spline_log_umi"):
            add_human_fit(rows, part, adjustment, "per_tissue", tissue=tissue)
    result = pd.DataFrame(rows)
    unit_frame = donors.copy()
    unit_frame["unit"] = (
        unit_frame["tissue"].astype(str) + " | " + unit_frame["cohort"].astype(str)
        + " | " + unit_frame["cell_type"].astype(str)
    )
    unit_stats = unit_frame.groupby("unit")["profile_hellinger"].agg(["size", "std"])
    valid_units = set(unit_stats[(unit_stats["size"] >= 6) & (unit_stats["std"] > 0)].index)
    global_donors = unit_frame[unit_frame["unit"].isin(valid_units)]
    donor_counts = global_donors.groupby("cohort")["donor"].nunique()
    reference_counts = {item["cohort"]: item["n_reference_donors"] for item in infos}
    donor_label = ";".join(f"{cohort}:{int(donor_counts.get(cohort, 0))}" for cohort in HUMAN_COHORTS)
    reference_label = ";".join(f"{cohort}:{reference_counts[cohort]}" for cohort in HUMAN_COHORTS)
    result.loc[result["variant"] == "global", "donors_by_cohort"] = donor_label
    result.loc[result["variant"] == "global", "reference_donors_by_cohort"] = reference_label
    result.to_csv(HUMAN_CSV, index=False)
    print(result[result["variant"].isin(["global", "per_cohort"])][
        ["variant", "cohort", "adjustment", "status", "slope", "pvalue", "n_donors"]
    ].to_string(index=False), flush=True)


def interaction_pipeline() -> None:
    from scripts.run_rev2_concordance import (
        cluster_bootstrap_interaction,
        pooled_interaction,
        pooled_prep,
    )
    mouse = mature_full_mouse().rename(columns={"profile_hellinger_centroid_mean": "profile_hellinger"})
    human = pd.read_csv(HUMAN_DONORS_CSV)
    columns = list(human.columns)
    human = human[columns]
    combined = pd.concat([mouse, human], ignore_index=True)
    assert combined.loc[combined["species"] == "Mouse", "age"].min() >= 1.5
    assert combined.loc[combined["species"] == "Human", "age"].min() >= 30
    rows = []
    for model in ("unadjusted", "depth_adjusted"):
        prepared = pooled_prep(combined, "profile_hellinger", model, "published")
        fit, X = pooled_interaction(prepared, model)
        boot = cluster_bootstrap_interaction(
            prepared, X, model, reps=1000, seed=20261121
        )
        interaction = fit["age_x_human"]
        mouse_fit = fit["age_scaled"]
        rows.append({
            "analysis": "primary_mature_published_unadjusted" if model == "unadjusted" else "sensitivity_mature_published_depth_adjusted",
            "status": "NEW",
            "model": model,
            "age_scale_mouse_years": 0.4,
            "age_scale_human_years": 12.25,
            "mouse_effect_sd": mouse_fit["effect"],
            "human_effect_sd": mouse_fit["effect"] + interaction["effect"],
            "interaction_human_minus_mouse": interaction["effect"],
            "interaction_se": interaction["se"],
            "interaction_ci95_low": interaction["ci_low"],
            "interaction_ci95_high": interaction["ci_high"],
            "interaction_p": interaction["p"],
            "interaction_bootstrap_ci95_low": float(np.quantile(boot, 0.025)),
            "interaction_bootstrap_ci95_high": float(np.quantile(boot, 0.975)),
            "interaction_bootstrap_p_two_sided": float(min(1.0, 2 * min((boot <= 0).mean(), (boot >= 0).mean()))),
            "interaction_bootstrap_reps": 1000,
            "n_mouse_donors": prepared.loc[prepared["species"] == "Mouse", "donor_cluster"].nunique(),
            "n_human_donors": prepared.loc[prepared["species"] == "Human", "donor_cluster"].nunique(),
            "n_donor_units": len(prepared),
            "n_clusters": fit["n_clusters"],
            "provenance": "new mature mouse and mature human donor arms only",
            "notes": "Pooled OLS with unit fixed effects, species-specific sex and depth terms, and donor-cluster-robust inference.",
        })
    result = pd.DataFrame(rows)
    result.to_csv(INTERACTION_CSV, index=False)
    print(result.to_string(index=False), flush=True)


def reuse_pipeline() -> None:
    rows = []
    posterior_path = RESULTS / "rev3-bayesian-vector-mature-only-posterior.csv"
    posterior = pd.read_csv(posterior_path)
    for item in posterior[posterior["scope"] == "species"].itertuples():
        rows.append({
            "arm": "mouse_full_depth_bayesian",
            "disposition": "REUSED_AS_IS",
            "file_path": str(posterior_path.relative_to(REPO)),
            "spec_id": item.spec_id,
            "endpoint": item.endpoint,
            "model": item.model,
            "effect": item.posterior_mean,
            "se": item.posterior_sd,
            "interval_low": item.hdi94_low,
            "interval_high": item.hdi94_high,
            "interval_type": "94% highest density interval",
            "pvalue": np.nan,
            "probability_positive": item.prob_positive,
            "n_donors": item.n_donors,
            "n_units": item.n_units,
            "n_obs": item.n_obs,
            "notes": "Posterior species row copied without a new Bayesian fit.",
        })
    frequentist_path = RESULTS / "rev3-bayesian-vector-mature-only-frequentist.csv"
    frequentist = pd.read_csv(frequentist_path)
    for item in frequentist.itertuples():
        rows.append({
            "arm": "mouse_full_depth_frequentist",
            "disposition": "REUSED_AS_IS",
            "file_path": str(frequentist_path.relative_to(REPO)),
            "spec_id": item.spec_id,
            "endpoint": item.spec_id.split("-")[1],
            "model": "depth_adjusted" if "depth_adjusted" in item.spec_id else "unadjusted",
            "effect": item.freq_effect_10pct_lifespan,
            "se": item.freq_se,
            "interval_low": item.freq_ci95_low,
            "interval_high": item.freq_ci95_high,
            "interval_type": "95% CI",
            "pvalue": item.freq_pvalue,
            "probability_positive": np.nan,
            "n_donors": item.n_donors,
            "n_units": item.n_units,
            "n_obs": item.n_obs,
            "notes": "Frequentist row copied without a new full-depth global fit.",
        })
    entropy_path = RESULTS / "rev3-entropy-mature-only-mouse-lung.csv"
    entropy = pd.read_csv(entropy_path)
    for index, item in enumerate(entropy[entropy["scope"] == "global_donor_clustered"].itertuples()):
        rows.append({
            "arm": "mouse_lung_mature_entropy",
            "disposition": "REUSED_AS_IS",
            "file_path": str(entropy_path.relative_to(REPO)),
            "spec_id": f"global_row_{index}",
            "endpoint": item.metric,
            "model": f"{item.family}|{item.depth}|{item.model}",
            "effect": item.effect_15mo_sd,
            "se": item.se_15mo_sd,
            "interval_low": item.ci95_low_per_12mo * 15 / 12,
            "interval_high": item.ci95_high_per_12mo * 15 / 12,
            "interval_type": "95% CI per 15 months rescaled from reused 12-month columns",
            "pvalue": item.pvalue,
            "probability_positive": np.nan,
            "n_donors": item.n_donors,
            "n_units": item.n_cell_types,
            "n_obs": item.n_donor_units,
            "notes": "Existing mature global entropy result reused; no entropy recomputation.",
        })
    pd.DataFrame(rows).to_csv(REUSED_CSV, index=False)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=("all", "reuse", "mouse", "human", "interaction"), default="all")
    args = parser.parse_args()
    RESULTS.mkdir(exist_ok=True)
    SCRATCH.mkdir(parents=True, exist_ok=True)
    stages = {
        "reuse": reuse_pipeline,
        "mouse": mouse_pipeline,
        "human": human_pipeline,
        "interaction": interaction_pipeline,
    }
    selected = list(stages) if args.stage == "all" else [args.stage]
    for name in selected:
        run_stage(name, stages[name])
    return 0


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=FutureWarning)
    raise SystemExit(main())
