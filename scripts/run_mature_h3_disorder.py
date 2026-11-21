#!/usr/bin/env python3
"""Compute mature mouse lung off-program disorder and 30-plus De Man entropy."""

from __future__ import annotations

import os
import sys
import time
import warnings
from datetime import datetime, timezone
from pathlib import Path

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_followup_de_man_entropy import (
    SEEDS,
    TARGETS,
    binomial_downsample,
    kl_uniform,
    model_rows,
    multinomial_downsample,
    shannon,
)
from scripts.run_h3_transcriptional_disorder import (
    AGE_MAP,
    TARGET_CELL_TYPES,
    offprogram_markers,
    reference_stats,
)
from scripts.run_rev3_entropy_mature_mouse_lung import (
    fit_donor_celltype_ols,
    fit_global_clustered,
)

RESULTS = REPO / "results"
DATA = Path("/home/reza/projects/agescape-revision-data/de-man")
OUTPUT = RESULTS / "followup-mature-h3-disorder.csv"
DONOR_OUTPUT = RESULTS / "followup-mature-h3-disorder-donors.csv"
MARKER_OUTPUT = RESULTS / "followup-mature-h3-marker-provenance.csv"
RUNTIME_CSV = RESULTS / "followup-mature-drift-runtime.csv"
MODEL_USED = "kimi-code/k3.1-pro"
OFFPROGRAM_METRICS = (
    "offprogram_fraction",
    "offprogram_fraction_stringent",
    "offprogram_detected_genes",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def record_runtime(start_utc: str, started: float, status: str, details: str = "") -> None:
    new = pd.DataFrame([{
        "pipeline": "mature_h3_disorder",
        "start_utc": start_utc,
        "end_utc": utc_now(),
        "duration_seconds": time.monotonic() - started,
        "model_used": MODEL_USED,
        "status": status,
        "details": details,
    }])
    if RUNTIME_CSV.exists():
        new = pd.concat([pd.read_csv(RUNTIME_CSV), new], ignore_index=True)
    new.to_csv(RUNTIME_CSV, index=False)


def offprogram_donors() -> tuple[pd.DataFrame, pd.DataFrame]:
    path = REPO / "data" / "interim" / "m1-qc-cleaned-Lung.h5ad"
    atlas = ad.read_h5ad(path, backed="r")
    age_months = atlas.obs["age"].astype(str).map(AGE_MAP)
    mask = (
        (~atlas.obs["predicted_doublet"].to_numpy(dtype=bool))
        & (age_months >= 18).to_numpy(dtype=bool)
        & (atlas.obs["n_counts"] >= 500).to_numpy(dtype=bool)
        & (atlas.obs["cell_ontology_class"].astype(str) != "unknown").to_numpy(dtype=bool)
    )
    lung = atlas[mask].to_memory()
    atlas.file.close()
    lung.obs["age_months"] = lung.obs["age"].astype(str).map(AGE_MAP).astype(float)
    lung.obs["sex_male"] = (lung.obs["sex"].astype(str).str.lower() == "male").astype(float)
    lung.obs["mouse.id"] = lung.obs["mouse.id"].astype(str)
    lung.obs["actual_n_counts"] = np.asarray(lung.X.sum(axis=1)).ravel()
    assert lung.obs["age_months"].min() >= 18
    assert lung.obs["mouse.id"].nunique() == 10
    reference = lung[
        (lung.obs["age_months"] == 18) & (lung.obs["actual_n_counts"] >= 1500)
    ].copy()
    stats = reference_stats(reference.X.tocsr(), reference.obs)
    gene_names = lung.var_names.astype(str).to_numpy()
    marker_cache = {}
    provenance_cache = {}
    stringent_cache = {}
    for cell_type in TARGET_CELL_TYPES:
        for donor in sorted(lung.obs["mouse.id"].unique()):
            if cell_type not in stats:
                marker_cache[(cell_type, donor)] = np.array([], dtype=int)
                stringent_cache[(cell_type, donor)] = np.array([], dtype=int)
                provenance_cache[(cell_type, donor)] = {}
                continue
            markers, provenance = offprogram_markers(
                stats, cell_type, donor, gene_names, stringent=False
            )
            stringent, _ = offprogram_markers(
                stats, cell_type, donor, gene_names, stringent=True
            )
            marker_cache[(cell_type, donor)] = markers
            stringent_cache[(cell_type, donor)] = stringent
            provenance_cache[(cell_type, donor)] = provenance
    marker_rows = []
    for (cell_type, donor), provenance in provenance_cache.items():
        for source in sorted(set(provenance.values())):
            marker_rows.append({
                "cell_type": cell_type,
                "excluded_donor": donor,
                "source_cell_type": source,
                "n_marker_genes": sum(value == source for value in provenance.values()),
                "n_all_marker_genes": len(marker_cache[(cell_type, donor)]),
                "n_stringent_marker_genes": len(stringent_cache[(cell_type, donor)]),
                "reference_months": 18,
            })
    donor_rows = []
    for depth, minimum in ((500, 500), (1500, 1500), (3000, 3000), ("full", 1500)):
        selected = lung[lung.obs["actual_n_counts"] >= minimum].copy()
        if depth == "full":
            matrix = selected.X.tocsr()
        else:
            sc.pp.downsample_counts(
                selected, counts_per_cell=depth, replace=False, random_state=42
            )
            matrix = selected.X.tocsr()
            assert np.all(np.asarray(matrix.sum(axis=1)).ravel() == depth)
        totals = np.asarray(matrix.sum(axis=1)).ravel()
        rng = np.random.default_rng(20261121 + (0 if depth == "full" else int(depth)))
        for cell_type in TARGET_CELL_TYPES:
            ct_mask = (selected.obs["cell_ontology_class"].astype(str) == cell_type).to_numpy()
            ct_obs = selected.obs.loc[ct_mask]
            ct_matrix = matrix[ct_mask]
            for donor in sorted(ct_obs["mouse.id"].unique()):
                donor_mask = (ct_obs["mouse.id"] == donor).to_numpy()
                if donor_mask.sum() < 20:
                    continue
                donor_obs = ct_obs.loc[donor_mask]
                donor_matrix = ct_matrix[donor_mask]
                markers = marker_cache[(cell_type, donor)]
                stringent = stringent_cache[(cell_type, donor)]
                donor_totals = totals[ct_mask][donor_mask]
                offprogram = (
                    np.asarray(donor_matrix[:, markers].sum(axis=1)).ravel() / donor_totals
                    if len(markers) else np.full(donor_matrix.shape[0], np.nan)
                )
                offprogram_stringent = (
                    np.asarray(donor_matrix[:, stringent].sum(axis=1)).ravel() / donor_totals
                    if len(stringent) else np.full(donor_matrix.shape[0], np.nan)
                )
                detected = (
                    np.asarray((donor_matrix[:, markers] > 0).sum(axis=1)).ravel()
                    if len(markers) else np.full(donor_matrix.shape[0], np.nan)
                )
                values = {
                    "offprogram_fraction": offprogram,
                    "offprogram_fraction_stringent": offprogram_stringent,
                    "offprogram_detected_genes": detected,
                }
                for sampling in ("all_cells", "balanced_20_cells"):
                    indices = np.arange(donor_matrix.shape[0])
                    if sampling == "balanced_20_cells":
                        indices = rng.choice(indices, size=20, replace=False)
                    row = {
                        "depth": depth,
                        "sampling": sampling,
                        "cell_type": cell_type,
                        "donor": donor,
                        "age_months": float(donor_obs["age_months"].iloc[0]),
                        "sex_male": float(donor_obs["sex_male"].iloc[0]),
                        "n_cells": len(indices),
                        "mean_log_umi_actual": float(np.log(donor_obs["actual_n_counts"]).mean()),
                        "n_offprogram_genes": len(markers),
                        "n_stringent_offprogram_genes": len(stringent),
                        "reference_months": 18,
                    }
                    for metric, metric_values in values.items():
                        chosen = metric_values[indices]
                        row[metric] = float(np.nanmean(chosen)) if np.isfinite(chosen).any() else np.nan
                    donor_rows.append(row)
        print(f"Mature off-program depth {depth} completed", flush=True)
    donors = pd.DataFrame(donor_rows)
    assert donors["age_months"].min() >= 18
    assert donors["reference_months"].eq(18).all()
    return donors, pd.DataFrame(marker_rows)


def common_row(item: dict, domain: str, status: str, provenance: str) -> dict:
    scope = item.get("scope", "")
    is_global = scope == "global_donor_clustered"
    effect = item.get("effect_15mo_sd") if is_global else item.get("slope_per_month")
    se = item.get("se_15mo_sd") if is_global else item.get("se_per_month")
    low = item.get("ci95_low_per_12mo") * 15 / 12 if is_global else item.get("ci95_low_per_month")
    high = item.get("ci95_high_per_12mo") * 15 / 12 if is_global else item.get("ci95_high_per_month")
    mde = item.get("mde80_15mo_sd") if is_global else item.get("mde80_per_month")
    return {
        "domain": domain,
        "scope": scope,
        "family": item.get("family", ""),
        "depth": item.get("depth", ""),
        "sampling": item.get("sampling", ""),
        "metric": item.get("metric", ""),
        "model": item.get("model", ""),
        "cell_type": item.get("cell_type", ""),
        "status": status,
        "n_donors": item.get("n_donors", np.nan),
        "n_cell_types": item.get("n_cell_types", np.nan),
        "n_donor_units": item.get("n_donor_units", np.nan),
        "n_cells": item.get("n_cells", np.nan),
        "age_min": item.get("age_min_months", np.nan),
        "age_max": item.get("age_max_months", np.nan),
        "effect": effect,
        "se": se,
        "ci95_low": low,
        "ci95_high": high,
        "pvalue": item.get("pvalue", np.nan),
        "q_bh": item.get("q_bh_within_depth_family", np.nan),
        "mde80": mde,
        "effect_units": "within_cell_type_SD_per_15_months" if is_global else "metric_units_per_month",
        "mde80_units": "within_cell_type_SD_per_15_months" if is_global else "metric_units_per_month",
        "provenance": provenance,
        "notes": item.get("notes", ""),
    }


def offprogram_models(donors: pd.DataFrame) -> list[dict]:
    rows = []
    for depth in (500, 1500, 3000, "full"):
        depth_part = donors[donors["depth"] == depth]
        for sampling in ("all_cells", "balanced_20_cells"):
            sampled = depth_part[depth_part["sampling"] == sampling]
            for metric in OFFPROGRAM_METRICS:
                for model in ("age_sex", "age_sex_library"):
                    family = f"mature_18m_markers_{sampling}"
                    global_fit = fit_global_clustered(
                        sampled, metric, depth, model, family, min_donors_per_ct=1
                    )
                    if global_fit:
                        global_fit["sampling"] = sampling
                        global_fit["notes"] = "New 18-month cross-fitted marker sets; no entropy calculation."
                        rows.append(common_row(
                            global_fit, "mouse_lung_offprogram", "NEW",
                            "new 18m reference markers and 18-plus mice",
                        ))
                    cell_fits = fit_donor_celltype_ols(
                        sampled, metric, depth, model, family, min_donors=6
                    )
                    for fit in cell_fits:
                        fit["sampling"] = sampling
                        fit["notes"] = "New 18-month cross-fitted marker sets; no entropy calculation."
                        rows.append(common_row(
                            fit, "mouse_lung_offprogram", "NEW",
                            "new 18m reference markers and 18-plus mice",
                        ))
                    fitted_types = {fit["cell_type"] for fit in cell_fits}
                    for cell_type in TARGET_CELL_TYPES:
                        if cell_type in fitted_types:
                            continue
                        part = sampled[sampled["cell_type"] == cell_type].dropna(subset=[metric])
                        rows.append(common_row({
                            "scope": "per_cell_type",
                            "family": family,
                            "depth": depth,
                            "sampling": sampling,
                            "metric": metric,
                            "model": model,
                            "cell_type": cell_type,
                            "n_donors": len(part),
                            "n_cell_types": 1,
                            "n_donor_units": len(part),
                            "n_cells": part["n_cells"].sum() if len(part) else 0,
                            "age_min_months": part["age_months"].min() if len(part) else np.nan,
                            "age_max_months": part["age_months"].max() if len(part) else np.nan,
                            "notes": "UNTESTABLE: fewer than 6 eligible mature donors, no markers, or an inestimable design.",
                        }, "mouse_lung_offprogram", "UNTESTABLE", "new mature restriction"))
    return rows


def reused_entropy_rows() -> list[dict]:
    path = RESULTS / "rev3-entropy-mature-only-mouse-lung.csv"
    frame = pd.read_csv(path)
    rows = []
    for item in frame[frame["scope"] == "global_donor_clustered"].to_dict("records"):
        item["sampling"] = "all_cells"
        item["notes"] = "Reused mature global entropy output as-is; no entropy rerun."
        rows.append(common_row(
            item, "mouse_lung_entropy", "REUSED",
            "results/rev3-entropy-mature-only-mouse-lung.csv",
        ))
    return rows


def de_man_per_cell() -> pd.DataFrame:
    cache = DATA / "target_matrices"
    required = [cache / f"{target}.npz" for target in TARGETS] + [cache / "selected_metadata.csv"]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"De Man cached matrices are unavailable: {missing}")
    metadata = pd.read_csv(
        DATA / "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz", sep="\t"
    )
    metadata["matrix_column"] = np.arange(len(metadata))
    metadata = metadata[metadata["Disease_Identity"] == "Control"].copy()
    selected = pd.read_csv(cache / "selected_metadata.csv")
    metadata = metadata.merge(selected, on="matrix_column", how="left")
    subjects = pd.read_csv(DATA / "41467_2026_68810_MOESM3_ESM.csv").rename(
        columns={"Sample ID": "sample_id", "Sex": "sex", "Age": "age"}
    )
    metadata = metadata.merge(
        subjects[["sample_id", "sex", "age"]],
        left_on="Subject_Identity", right_on="sample_id", how="left"
    )
    gene_info = pd.read_csv(DATA / "GSE136831_AllCells.GeneIDs.txt.gz", sep="\t")
    symbols = gene_info["HGNC_EnsemblAlt_GeneID"].astype(str).str.upper()
    technical = (
        symbols.str.match(r"^(MT-|RP[LS]|MRP[LS])")
        | symbols.isin(["XIST", "TSIX", "DDX3Y", "EIF2S3Y", "UTY", "KDM5D"])
    ).to_numpy()
    per_cell = []
    for target, label in TARGETS.items():
        obs = metadata[metadata["Manuscript_Identity"] == target].copy()
        matrix = sp.load_npz(cache / f"{target}.npz").tocsr()
        assert len(obs) == matrix.shape[0]
        keep = obs["qc_pass"].to_numpy(dtype=bool) & (obs["age"].to_numpy(dtype=float) >= 30)
        obs = obs.loc[keep].copy()
        matrix = matrix[keep]
        assert len(obs) and obs["age"].min() >= 30
        shannon_matrix = matrix[:, ~technical]
        matrix = matrix[:, np.asarray(matrix.sum(axis=0)).ravel() > 0]
        shannon_matrix = shannon_matrix[:, np.asarray(shannon_matrix.sum(axis=0)).ravel() > 0]
        full_shannon = shannon(shannon_matrix)
        depth_totals = np.asarray(shannon_matrix.sum(axis=1)).ravel()
        fixed_idx = np.flatnonzero(depth_totals >= 1500)
        for seed in SEEDS:
            rng = np.random.default_rng(seed)
            thinned = binomial_downsample(matrix, rng)
            young = obs["age"].to_numpy() <= 60
            n_sample = min(young.sum(), (~young).sum())
            if n_sample >= 10:
                old_idx = rng.choice(np.flatnonzero(~young), size=n_sample, replace=False)
                young_idx = rng.choice(np.flatnonzero(young), size=n_sample, replace=False)
                indices = np.concatenate([old_idx, young_idx])
                values = kl_uniform(thinned[indices], sqrt=True)
                for index, value in zip(indices, values):
                    item = obs.iloc[index]
                    per_cell.append({
                        "seed": seed, "method": "published_KL_binomial", "cell_type": label,
                        "sample_id": item["sample_id"], "age": item["age"], "sex": item["sex"],
                        "actual_umi": item["actual_umi"], "entropy": value,
                    })
            fixed_shannon = shannon(multinomial_downsample(shannon_matrix[fixed_idx], 1500, rng))
            for method, indices, values in (
                ("shannon_full", np.arange(len(obs)), full_shannon),
                ("shannon_1500", fixed_idx, fixed_shannon),
            ):
                for local_index, value in zip(indices, values):
                    item = obs.iloc[local_index]
                    per_cell.append({
                        "seed": seed, "method": method, "cell_type": label,
                        "sample_id": item["sample_id"], "age": item["age"], "sex": item["sex"],
                        "actual_umi": item["actual_umi"], "entropy": value,
                    })
        print(f"De Man 30-plus {label} completed", flush=True)
    frame = pd.DataFrame(per_cell)
    assert frame["age"].min() >= 30
    return frame


def de_man_models() -> list[dict]:
    per_cell = de_man_per_cell()
    rows = []
    for (seed, method), part in per_cell.groupby(["seed", "method"]):
        for fit in model_rows(part, seed, method):
            donor_ages = part.groupby("sample_id")["age"].first()
            fit.update({
                "age_min": donor_ages.min(), "age_max": donor_ages.max(),
            })
            rows.append(fit)
    models = pd.DataFrame(rows)
    primary = models[models["seed"] == SEEDS[0]].copy()
    primary["q_bh"] = np.nan
    for _, indices in primary.groupby("method").groups.items():
        primary.loc[indices, "q_bh"] = multipletests(primary.loc[indices, "p"], method="fdr_bh")[1]
    models = models.merge(
        primary[["seed", "method", "cell_type", "q_bh"]],
        on=["seed", "method", "cell_type"], how="left"
    )
    output = []
    for item in models.to_dict("records"):
        output.append({
            "domain": "de_man_30plus",
            "scope": "per_cell_type_seed_model",
            "family": "De_man_30plus_20_seeds",
            "depth": item["method"],
            "sampling": "age_balanced_cells" if item["method"] == "published_KL_binomial" else "all_eligible_cells",
            "metric": item["method"],
            "model": "age_sex_log_umi",
            "cell_type": item["cell_type"],
            "status": "NEW",
            "seed": item["seed"],
            "n_donors": item["n_donors"],
            "n_cell_types": 1,
            "n_donor_units": item["n_donors"],
            "n_cells": item["n_cells"],
            "age_min": item["age_min"],
            "age_max": item["age_max"],
            "effect": item["slope_per_year"],
            "se": item["se"],
            "ci95_low": item["ci95_low"],
            "ci95_high": item["ci95_high"],
            "pvalue": item["p"],
            "q_bh": item["q_bh"],
            "mde80": np.nan,
            "effect_units": "entropy_units_per_year",
            "mde80_units": "",
            "provenance": "new 30-plus refit from cached raw-count matrices",
            "notes": "Age restriction applied before rarefaction, nonzero gene selection, balancing, and regression.",
        })
    for method in sorted(models["method"].unique()):
        for cell_type in sorted(TARGETS.values()):
            part = models[(models["method"] == method) & (models["cell_type"] == cell_type)]
            if part.empty:
                output.append({
                    "domain": "de_man_30plus", "scope": "per_cell_type_seed_model",
                    "family": "De_man_30plus_20_seeds", "depth": method,
                    "metric": method, "model": "age_sex_log_umi", "cell_type": cell_type,
                    "status": "UNTESTABLE", "effect_units": "entropy_units_per_year",
                    "provenance": "new 30-plus refit from cached raw-count matrices",
                    "notes": "UNTESTABLE: fewer than 6 donors or fewer than 10 cells in an age-balanced group at 30-plus.",
                })
                continue
            output.append({
                "domain": "de_man_30plus",
                "scope": "seed_robustness",
                "family": "De_man_30plus_20_seeds",
                "depth": method,
                "sampling": "20_seeds",
                "metric": method,
                "model": "age_sex_log_umi",
                "cell_type": cell_type,
                "status": "NEW",
                "n_seeds": part["seed"].nunique(),
                "n_donors": part["n_donors"].median(),
                "n_cell_types": 1,
                "effect": part["slope_per_year"].median(),
                "pvalue": part["p"].median(),
                "q_bh": np.nan,
                "effect_units": "median_entropy_units_per_year",
                "provenance": "new 30-plus refit from cached raw-count matrices",
                "seed_slope_min": part["slope_per_year"].min(),
                "seed_slope_max": part["slope_per_year"].max(),
                "seed_p_min": part["p"].min(),
                "seed_p_max": part["p"].max(),
                "fraction_negative": (part["slope_per_year"] < 0).mean(),
                "fraction_positive": (part["slope_per_year"] > 0).mean(),
                "n_p_lt_005": int((part["p"] < 0.05).sum()),
                "notes": "Seed ranges are robustness summaries, not confidence intervals.",
            })
    return output


def main() -> int:
    started = time.monotonic()
    start_utc = utc_now()
    print(f"Mature H3 disorder started {start_utc}; model {MODEL_USED}", flush=True)
    try:
        donors, markers = offprogram_donors()
        donors.to_csv(DONOR_OUTPUT, index=False)
        markers.to_csv(MARKER_OUTPUT, index=False)
        rows = offprogram_models(donors)
        rows.extend(reused_entropy_rows())
        rows.extend(de_man_models())
        result = pd.DataFrame(rows)
        result.to_csv(OUTPUT, index=False)
    except Exception as exc:
        record_runtime(start_utc, started, "FAILED", f"{type(exc).__name__}: {exc}")
        raise
    record_runtime(start_utc, started, "completed")
    print(f"Mature H3 disorder completed after {time.monotonic() - started:.3f} seconds", flush=True)
    primary = result[
        (result["domain"] == "mouse_lung_offprogram")
        & (result["scope"] == "global_donor_clustered")
        & (result["metric"] == "offprogram_fraction")
    ]
    print(primary[["depth", "sampling", "model", "effect", "pvalue", "n_donors"]].to_string(index=False), flush=True)
    return 0


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=FutureWarning)
    raise SystemExit(main())
