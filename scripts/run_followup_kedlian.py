#!/usr/bin/env python3
"""Validate and analyze the public Kedlian raw-count H5AD by donor and assay."""

from __future__ import annotations

import argparse
import hashlib
import json
import resource
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from src.aging_programs import load_program_tiers  # noqa: E402
from src.followup_muscle import (  # noqa: E402
    MIN_CELLS, MIN_DONORS, MIN_PER_GROUP, MIN_REFERENCE_CELLS,
    MIN_RESIDUAL_DF, MIN_UMI_HELLINGER, NULL_ITERATIONS, NULL_SEED,
    add_bh, aggregate_block, assert_count_gate, build_design, clr_transform,
    fit_ols, hellinger_distances, known_cell_type, matched_null_effect_pvalue,
    random_signature_scores, technical_gene_mask,
)

DATA = REPO.parent / f"{REPO.name}-data" / "muscle"
H5AD = DATA / "kedlian" / "SKM_human_raw_cells2nuclei_2023-06-22.h5ad"
CACHE = DATA / "interim" / "kedlian_donor_aggregates.npz"
MANIFEST = REPO / "data" / "manifest" / "E-MTAB-13874.json"
RESULTS = REPO / "results"
PREFIX = "followup-muscle-kedlian"
EXPECTED_SIZE = 4025846616
CHUNK = 256
MODELS = {
    "primary_adult": ["primary_group_depth_adjusted", "unadjusted_group_sensitivity"],
    "all_donor": ["all_donor_group_depth_adjusted"],
}


def write_csv(frame: pd.DataFrame, name: str) -> None:
    path = RESULTS / name
    frame.to_csv(path, index=False)
    print(f"wrote {path}: {len(frame)} rows", flush=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def donor_metadata() -> pd.DataFrame:
    manifest = json.loads(MANIFEST.read_text())
    rows = []
    for donor in manifest["donors"]:
        for assay in ["scRNA", "snRNA"]:
            rows.append({
                "donor_id": donor["id"], "assay": assay,
                "sample_ids": donor[f"{assay}_sample_ids"],
                "age_band": donor["age_band"], "group": donor["group"],
                "old": int(donor["group"] == "old"),
                "sex": donor["sex"], "male": int(donor["sex"] == "M"),
                "donor_type": donor["donor_type"], "dcd": int(donor["donor_type"] == "DCD"),
                "primary_adult": donor["id"] != "470BR",
            })
    return pd.DataFrame(rows)


def inspect_object() -> None:
    if not H5AD.is_file() or H5AD.stat().st_size != EXPECTED_SIZE:
        raise ValueError("Kedlian H5AD is absent or has an unexpected size")
    adata = ad.read_h5ad(H5AD, backed="r")
    try:
        obs = adata.obs
        columns = {}
        for column in obs.columns:
            values = obs[column].astype(str)
            unique = sorted(values.unique())
            columns[column] = {
                "dtype": str(obs[column].dtype), "n_unique": int(values.nunique()),
                "values": unique if len(unique) <= 100 else unique[:100],
            }
        schema = {
            "n_obs": int(adata.n_obs), "n_vars": int(adata.n_vars),
            "obs_columns": columns,
            "var_columns": list(adata.var.columns.astype(str)),
            "var_names_unique": bool(adata.var_names.is_unique),
            "var_names_examples": adata.var_names.astype(str)[:20].tolist(),
            "layers": list(adata.layers.keys()),
            "counts_type": str(type(adata.layers["counts"])),
        }
    finally:
        adata.file.close()
    path = RESULTS / f"{PREFIX}-schema.json"
    path.write_text(json.dumps(schema, indent=2, sort_keys=True) + "\n")
    print(json.dumps(schema, indent=2, sort_keys=True), flush=True)


def load_metadata() -> tuple[ad.AnnData, pd.DataFrame, pd.DataFrame, dict]:
    adata = ad.read_h5ad(H5AD, backed="r")
    obs = adata.obs.copy()
    manifest = json.loads(MANIFEST.read_text())
    sample_assay = {}
    sample_donor = {}
    for row in manifest["library_to_ena"]:
        sample_assay[row["sample_id"]] = row["assay"]
        sample_donor[row["sample_id"]] = row["donor_id"]
    scores = {
        column: obs[column].astype(str).isin(sample_assay).mean()
        for column in obs.columns
    }
    sample_column = max(scores, key=scores.get)
    if scores[sample_column] != 1.0:
        adata.file.close()
        raise ValueError(f"No obs sample column maps every row to the 43 released libraries: {scores}")
    obs["sample_id"] = obs[sample_column].astype(str)
    obs["assay"] = obs["sample_id"].map(sample_assay)
    obs["donor_id"] = obs["sample_id"].map(sample_donor)
    donors = donor_metadata()
    donor_assay = donors.set_index(["donor_id", "assay"])
    for row in manifest["library_to_ena"]:
        expected_samples = set(donor_assay.loc[(row["donor_id"], row["assay"]), "sample_ids"])
        if row["sample_id"] not in expected_samples:
            adata.file.close()
            raise ValueError(f"Manifest sample mapping mismatch: {row['sample_id']}")
    unexpected = set(zip(obs["donor_id"], obs["assay"])) - set(donor_assay.index)
    if unexpected:
        adata.file.close()
        raise ValueError(f"Unexpected donor-assay mappings: {sorted(unexpected)}")
    meta = donors.set_index(["donor_id", "assay"])
    keys = pd.MultiIndex.from_arrays([obs["donor_id"], obs["assay"]])
    for column in ["age_band", "group", "old", "sex", "male", "donor_type", "dcd", "primary_adult"]:
        obs[column] = meta[column].reindex(keys).to_numpy()
    obs["annotation"] = ""
    for assay, field in [("scRNA", "cell"), ("snRNA", "nuclei")]:
        if field not in obs:
            adata.file.close()
            raise ValueError(f"Released annotation field is absent: {field}")
        mask = obs["assay"].eq(assay)
        obs.loc[mask, "annotation"] = obs.loc[mask, field].astype(str).str.strip()
    obs["known_annotation"] = known_cell_type(obs["annotation"])
    obs["_row"] = np.arange(len(obs))
    verification = {
        "sample_column": sample_column, "all_rows_mapped_to_released_libraries": True,
        "annotation_fields": {"scRNA": "cell", "snRNA": "nuclei"},
        "released_donor_metadata_source": "data/manifest/E-MTAB-13874.json",
        "numeric_ages_used": False,
    }
    return adata, obs, donors, verification


def eligibility(obs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    selected = obs[obs["known_annotation"]]
    for assay in ["scRNA", "snRNA"]:
        assay_obs = selected[selected["assay"].eq(assay)]
        for donor_set in ["primary_adult", "all_donor"]:
            part0 = assay_obs[assay_obs["primary_adult"].to_numpy(bool) | (donor_set == "all_donor")]
            counts = part0.groupby(["annotation", "donor_id"]).size().rename("n_cells").reset_index()
            for cell_type, part in counts.groupby("annotation"):
                eligible = part[part["n_cells"] >= MIN_CELLS]
                groups = eligible.groupby("donor_id")["n_cells"].first().reset_index()
                donor_groups = part0[["donor_id", "old"]].drop_duplicates().set_index("donor_id")["old"]
                groups["old"] = groups["donor_id"].map(donor_groups)
                n_young = int((groups["old"] == 0).sum())
                n_old = int((groups["old"] == 1).sum())
                passed = len(groups) >= MIN_DONORS and n_young >= MIN_PER_GROUP and n_old >= MIN_PER_GROUP
                rows.append({
                    "assay": assay, "donor_set": donor_set, "cell_type": cell_type,
                    "n_donors_observed": int(len(part)), "n_donors_ge20": int(len(groups)),
                    "n_young_ge20": n_young, "n_old_ge20": n_old,
                    "n_cells_total": int(part["n_cells"].sum()), "eligible": passed,
                    "reason": "pass" if passed else "failed_six_donor_two_per_group_20_cell_floors",
                })
    return pd.DataFrame(rows).sort_values(["assay", "donor_set", "cell_type"]).reset_index(drop=True)


def composition_counts(obs: pd.DataFrame) -> pd.DataFrame:
    selected = obs[obs["known_annotation"]].copy()
    counts = selected.groupby(["assay", "donor_id", "annotation"]).size().rename("n_cells").reset_index()
    rows = []
    for assay, part in counts.groupby("assay"):
        wide = part.pivot(index="donor_id", columns="annotation", values="n_cells").fillna(0).astype(int)
        meta = selected[selected["assay"].eq(assay)].groupby("donor_id").agg(
            old=("old", "first"), male=("male", "first"), dcd=("dcd", "first"),
            primary_adult=("primary_adult", "first"), age_band=("age_band", "first"),
            group=("group", "first"), sex=("sex", "first"), donor_type=("donor_type", "first"),
        )
        wide = wide.reindex(meta.index, fill_value=0)
        long = wide.stack().rename("n_cells").reset_index().rename(columns={"annotation": "cell_type"})
        long = long.merge(meta.reset_index(), on="donor_id", validate="many_to_one")
        long.insert(0, "assay", assay)
        long["total_cells"] = long.groupby("donor_id")["n_cells"].transform("sum")
        rows.append(long)
    return pd.concat(rows, ignore_index=True)


def build_cache() -> None:
    started = time.time()
    adata, obs, donors, verification = load_metadata()
    try:
        if adata.n_obs != len(obs) or adata.layers["counts"].dtype != np.dtype("float32"):
            raise ValueError("Unexpected Kedlian count-layer schema")
        eligible = eligibility(obs)
        write_csv(eligible, f"{PREFIX}-celltype-eligibility.csv")
        comp = composition_counts(obs)
        write_csv(comp, f"{PREFIX}-composition-donor-counts.csv")
        eligible_pairs = set(eligible.loc[eligible["eligible"], ["assay", "cell_type"]].itertuples(index=False, name=None))
        selected_mask = np.array([
            (assay, annotation) in eligible_pairs
            for assay, annotation in zip(obs["assay"], obs["annotation"])
        ])
        selected = obs[selected_mask].copy()
        units = selected.groupby(["assay", "annotation", "donor_id"], sort=True).agg(
            old=("old", "first"), male=("male", "first"), dcd=("dcd", "first"),
            primary_adult=("primary_adult", "first"), age_band=("age_band", "first"),
            group=("group", "first"), sex=("sex", "first"), donor_type=("donor_type", "first"),
        ).reset_index().rename(columns={"annotation": "cell_type"})
        unit_lookup = {(r.assay, r.cell_type, r.donor_id): i for i, r in enumerate(units.itertuples())}
        codes = np.full(len(obs), -1, dtype=np.int64)
        codes[selected["_row"].to_numpy()] = [
            unit_lookup[(a, c, d)] for a, c, d in zip(selected["assay"], selected["annotation"], selected["donor_id"])
        ]
        known = obs["known_annotation"].to_numpy()
        comp_donors = (obs.loc[known, ["assay", "donor_id"]].drop_duplicates()
                       .sort_values(["assay", "donor_id"]).reset_index(drop=True))
        comp_lookup = {(r.assay, r.donor_id): i for i, r in enumerate(comp_donors.itertuples())}
        comp_codes = np.full(len(obs), -1, dtype=np.int64)
        comp_codes[known] = [comp_lookup[(a, d)] for a, d in zip(obs.loc[known, "assay"], obs.loc[known, "donor_id"])]
        gene_ids = adata.var_names.astype(str).to_numpy(dtype=str)
        gene_names = (adata.var["feature_name"].astype(str).to_numpy(dtype=str)
                      if "feature_name" in adata.var else gene_ids.copy())
        if not adata.var_names.is_unique:
            raise ValueError("Kedlian feature IDs are not unique")
        n_units, n_genes, n_comp = len(units), len(gene_ids), len(comp_donors)
        raw_counts = np.zeros((n_units, n_genes), dtype=np.int64)
        log_sums = np.zeros((n_units, n_genes), dtype=np.float32)
        qc_counts = np.zeros((n_units, n_genes), dtype=np.int64)
        n_cells = np.zeros(n_units, dtype=np.int64)
        qc_n = np.zeros(n_units, dtype=np.int64)
        sum_log_umi = np.zeros(n_units, dtype=np.float64)
        sum_umi = np.zeros(n_units, dtype=np.float64)
        comp_n = np.zeros(n_comp, dtype=np.int64)
        comp_log = np.zeros(n_comp, dtype=np.float64)
        raw_nnz, raw_sum = 0, 0.0
        counts_layer = adata.layers["counts"]
        for start in range(0, adata.n_obs, CHUNK):
            stop = min(start + CHUNK, adata.n_obs)
            matrix = sp.csr_matrix(counts_layer[start:stop])
            values = matrix.data
            assert_count_gate(values)
            raw_nnz += values.size
            raw_sum += float(values.sum(dtype=np.float64))
            local_codes = codes[start:stop]
            chosen = np.flatnonzero(local_codes >= 0)
            if len(chosen):
                sub = matrix[chosen].tocsr().astype(np.float64)
                sub_codes = local_codes[chosen]
                libraries = np.asarray(sub.sum(axis=1)).ravel()
                if (libraries <= 0).any():
                    raise ValueError("Eligible donor-cell-type row has zero raw UMI")
                np.add.at(n_cells, sub_codes, 1)
                np.add.at(sum_log_umi, sub_codes, np.log(libraries))
                np.add.at(sum_umi, sub_codes, libraries)
                log_sub = sub.copy()
                log_sub.data *= np.repeat(1e6 / libraries, np.diff(log_sub.indptr))
                np.log1p(log_sub.data, out=log_sub.data)
                aggregate_block(log_sub, sub_codes, log_sums)
                aggregate_block(sub, sub_codes, raw_counts)
                qc = libraries >= MIN_UMI_HELLINGER
                if qc.any():
                    aggregate_block(sub[qc], sub_codes[qc], qc_counts)
                    np.add.at(qc_n, sub_codes[qc], 1)
            local_comp = comp_codes[start:stop]
            chosen_comp = np.flatnonzero(local_comp >= 0)
            if len(chosen_comp):
                sub = matrix[chosen_comp].tocsr().astype(np.float64)
                libraries = np.asarray(sub.sum(axis=1)).ravel()
                if (libraries <= 0).any():
                    raise ValueError("Annotated composition row has zero raw UMI")
                np.add.at(comp_n, local_comp[chosen_comp], 1)
                np.add.at(comp_log, local_comp[chosen_comp], np.log(libraries))
            if start % 10240 == 0:
                print(f"Kedlian counts rows {start}:{stop}; {raw_nnz:,} stored values checked", flush=True)
        if (n_cells <= 0).any() or (comp_n <= 0).any():
            raise ValueError("Empty aggregation unit")
        cache = {
            "gene_ids": gene_ids, "gene_names": gene_names,
            "assay": units["assay"].astype(str).to_numpy(dtype=str),
            "cell_type": units["cell_type"].astype(str).to_numpy(dtype=str),
            "donor_id": units["donor_id"].astype(str).to_numpy(dtype=str),
            "old": units["old"].to_numpy(np.int8), "male": units["male"].to_numpy(np.int8),
            "dcd": units["dcd"].to_numpy(np.int8),
            "primary_adult": units["primary_adult"].to_numpy(bool),
            "age_band": units["age_band"].astype(str).to_numpy(dtype=str),
            "group": units["group"].astype(str).to_numpy(dtype=str),
            "sex": units["sex"].astype(str).to_numpy(dtype=str),
            "donor_type": units["donor_type"].astype(str).to_numpy(dtype=str),
            "n_cells": n_cells, "qc_n": qc_n, "mean_log_umi": sum_log_umi / n_cells,
            "sum_umi": sum_umi, "raw_counts": raw_counts, "log_sums": log_sums,
            "qc_counts": qc_counts, "obs_unit": codes,
            "composition_assay": comp_donors["assay"].astype(str).to_numpy(dtype=str),
            "composition_donor_id": comp_donors["donor_id"].astype(str).to_numpy(dtype=str),
            "composition_mean_log_umi": comp_log / comp_n,
            "composition_n_cells": comp_n,
        }
        CACHE.parent.mkdir(parents=True, exist_ok=True)
        np.savez(CACHE, **cache)
        verification.update({
            "file_size_bytes": H5AD.stat().st_size, "n_obs": int(adata.n_obs),
            "n_vars": int(adata.n_vars), "count_layer": "layers/counts",
            "count_dtype": "float32", "stored_values_checked": int(raw_nnz),
            "stored_values_sum": raw_sum, "all_counts_finite_nonnegative_integer_valued": True,
            "feature_ids_unique": True, "duplicate_feature_name_rows": int(pd.Index(gene_names).duplicated(keep=False).sum()),
            "aggregation_units": n_units, "eligible_assay_celltype_pairs": len(eligible_pairs),
            "cache_bytes": CACHE.stat().st_size, "elapsed_seconds": time.time() - started,
            "peak_rss_bytes": int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024),
            "memory_limit_gb": 6, "minimum_free_gb": 30,
            "donor_counts": {
                assay: {
                    "all": int(donors[(donors["assay"] == assay) & donors["sample_ids"].map(bool)]["donor_id"].nunique()),
                    "primary": int(donors[(donors["assay"] == assay) & donors["sample_ids"].map(bool) & donors["primary_adult"]]["donor_id"].nunique()),
                } for assay in ["scRNA", "snRNA"]
            },
        })
        (RESULTS / f"{PREFIX}-verification.json").write_text(json.dumps(verification, indent=2, sort_keys=True) + "\n")
    finally:
        adata.file.close()
    print("Kedlian raw-count gate and donor cache completed; no outcomes fitted in this command.", flush=True)


def load_cache() -> dict:
    with np.load(CACHE, allow_pickle=False) as data:
        return {key: data[key] for key in data.files}


def unit_frame(cache: dict, indices: np.ndarray) -> pd.DataFrame:
    return pd.DataFrame({
        "assay": cache["assay"][indices].astype(str),
        "cell_type": cache["cell_type"][indices].astype(str),
        "donor_id": cache["donor_id"][indices].astype(str),
        "old": cache["old"][indices], "male": cache["male"][indices], "dcd": cache["dcd"][indices],
        "age_band": cache["age_band"][indices].astype(str), "group": cache["group"][indices].astype(str),
        "sex": cache["sex"][indices].astype(str), "donor_type": cache["donor_type"][indices].astype(str),
        "primary_adult": cache["primary_adult"][indices], "n_cells": cache["n_cells"][indices],
        "mean_log_umi": cache["mean_log_umi"][indices], "sum_umi": cache["sum_umi"][indices],
    })


def selected_units(cache: dict, eligible: pd.DataFrame, assay: str, donor_set: str,
                   cell_type: str, qc: bool = False) -> np.ndarray:
    row = eligible[(eligible["assay"] == assay) & (eligible["donor_set"] == donor_set)
                   & (eligible["cell_type"] == cell_type)]
    if row.empty or not bool(row["eligible"].iloc[0]):
        return np.array([], dtype=int)
    mask = (cache["assay"].astype(str) == assay) & (cache["cell_type"].astype(str) == cell_type)
    if donor_set == "primary_adult":
        mask &= cache["primary_adult"].astype(bool)
    counts = cache["qc_n"] if qc else cache["n_cells"]
    mask &= counts >= MIN_CELLS
    indices = np.flatnonzero(mask)
    return indices[np.argsort(cache["donor_id"][indices].astype(str))]


def fit_frame(frame: pd.DataFrame, outcome: np.ndarray, model: str) -> tuple[dict, np.ndarray]:
    depth = model != "unadjusted_group_sensitivity"
    design, names, omitted = build_design(frame, depth=depth)
    fitted = fit_ols(np.asarray(outcome), design)
    result = {key: float(value[0]) for key, value in fitted.items()}
    result.update({
        "n_donors": len(frame), "n_young": int((frame["old"] == 0).sum()),
        "n_old": int((frame["old"] == 1).sum()), "n_cells": int(frame["n_cells"].sum()),
        "design_columns": ";".join(names), "omitted_covariates": ";".join(omitted),
    })
    return result, design


def status_row(assay: str, donor_set: str, model: str, endpoint: str, cell_type: str,
               status: str, reason: str) -> dict:
    return {"assay": assay, "donor_set": donor_set, "model": model, "endpoint": endpoint,
            "cell_type": cell_type, "status": status, "reason": reason}


def run_p1(cache: dict, eligible: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[dict]]:
    program = load_program_tiers(REPO, "human")["P1_Clean"]
    genes = cache["gene_names"].astype(str)
    up = np.flatnonzero(np.isin(genes, program["up"]))
    down = np.flatnonzero(np.isin(genes, program["down"]))
    if not len(up) or not len(down):
        raise ValueError("P1 lacks detected up or down genes")
    donor_rows, model_rows, null_rows, statuses = [], [], [], []
    for assay in ["scRNA", "snRNA"]:
        for donor_set, models in MODELS.items():
            cell_types = eligible.loc[eligible["eligible"] & eligible["assay"].eq(assay)
                                      & eligible["donor_set"].eq(donor_set), "cell_type"].astype(str)
            for cell_type in sorted(cell_types):
                indices = selected_units(cache, eligible, assay, donor_set, cell_type)
                frame = unit_frame(cache, indices)
                matrix = (cache["log_sums"][indices] / cache["n_cells"][indices, None]).astype(float).T
                scores = matrix[up, :].mean(axis=0) - matrix[down, :].mean(axis=0)
                donor_frame = frame.copy()
                donor_frame.insert(1, "donor_set", donor_set)
                donor_frame["P1_score"] = scores
                donor_rows.append(donor_frame)
                for model in models:
                    try:
                        fitted, design = fit_frame(frame, scores, model)
                    except ValueError as exc:
                        statuses.append(status_row(assay, donor_set, model, "P1_transport", cell_type, "untestable", str(exc)))
                        continue
                    row = {"assay": assay, "donor_set": donor_set, "model": model,
                           "endpoint": "P1_transport", "cell_type": cell_type,
                           "p1_up_detected": len(up), "p1_down_detected": len(down), **fitted}
                    if donor_set == "primary_adult" and model == "primary_group_depth_adjusted":
                        random_scores = random_signature_scores(matrix, up, down)
                        null_effects = fit_ols(random_scores.T, design)["effect"]
                        observed = fitted["effect"]
                        row.update({
                            "random_iterations": NULL_ITERATIONS, "random_seed": NULL_SEED,
                            "random_mean_effect": float(null_effects.mean()),
                            "random_sd_effect": float(null_effects.std(ddof=1)),
                            "random_q025": float(np.quantile(null_effects, 0.025)),
                            "random_q975": float(np.quantile(null_effects, 0.975)),
                            "random_pvalue": matched_null_effect_pvalue(observed, null_effects),
                        })
                        for iteration, effect in enumerate(null_effects, 1):
                            null_rows.append({"assay": assay, "cell_type": cell_type,
                                              "iteration": iteration, "effect_old_minus_young": effect,
                                              "seed": NULL_SEED})
                    model_rows.append(row)
                    statuses.append(status_row(assay, donor_set, model, "P1_transport", cell_type, "completed", "pass"))
    models = add_bh(pd.DataFrame(model_rows), ["assay", "donor_set", "model"])
    models["random_q_bh"] = np.nan
    primary = models["donor_set"].eq("primary_adult") & models["model"].eq("primary_group_depth_adjusted")
    if primary.any():
        for assay in models.loc[primary, "assay"].unique():
            indices = models.index[primary & models["assay"].eq(assay)]
            from statsmodels.stats.multitest import multipletests
            models.loc[indices, "random_q_bh"] = multipletests(models.loc[indices, "random_pvalue"], method="fdr_bh")[1]
    return pd.concat(donor_rows, ignore_index=True), models, pd.DataFrame(null_rows), statuses


def run_composition(cache: dict, counts: pd.DataFrame, eligible: pd.DataFrame) -> tuple[pd.DataFrame, list[dict]]:
    depth = pd.DataFrame({
        "assay": cache["composition_assay"].astype(str),
        "donor_id": cache["composition_donor_id"].astype(str),
        "mean_log_umi": cache["composition_mean_log_umi"],
        "composition_cells": cache["composition_n_cells"],
    })
    counts = counts.merge(depth, on=["assay", "donor_id"], validate="many_to_one")
    model_rows, statuses = [], []
    for assay in ["scRNA", "snRNA"]:
        for donor_set, models in MODELS.items():
            part = counts[counts["assay"].eq(assay) & (counts["primary_adult"].to_numpy(bool) | (donor_set == "all_donor"))].copy()
            wide = part.pivot(index="donor_id", columns="cell_type", values="n_cells").fillna(0)
            meta = part.groupby("donor_id").agg(old=("old", "first"), male=("male", "first"),
                                                dcd=("dcd", "first"), total_cells=("total_cells", "first"),
                                                mean_log_umi=("mean_log_umi", "first"))
            clr = clr_transform(wide)
            cell_types = eligible.loc[eligible["eligible"] & eligible["assay"].eq(assay)
                                      & eligible["donor_set"].eq(donor_set), "cell_type"].astype(str)
            for cell_type in sorted(cell_types):
                frame = meta.copy()
                frame["n_cells"] = wide[cell_type]
                for model in models:
                    try:
                        fitted, _ = fit_frame(frame, clr[cell_type].to_numpy(), model)
                    except ValueError as exc:
                        statuses.append(status_row(assay, donor_set, model, "CLR_composition", cell_type, "untestable", str(exc)))
                        continue
                    model_rows.append({"assay": assay, "donor_set": donor_set, "model": model,
                                       "endpoint": "CLR_composition", "cell_type": cell_type,
                                       "composition_universe_n_cell_types": wide.shape[1], **fitted})
                    statuses.append(status_row(assay, donor_set, model, "CLR_composition", cell_type, "completed", "pass"))
    return add_bh(pd.DataFrame(model_rows), ["assay", "donor_set", "model"]), statuses


def run_gene_models(cache: dict, eligible: pd.DataFrame) -> tuple[dict[str, pd.DataFrame], list[dict]]:
    outputs, statuses = {}, []
    for assay in ["scRNA", "snRNA"]:
        rows = []
        for donor_set, models in MODELS.items():
            cell_types = eligible.loc[eligible["eligible"] & eligible["assay"].eq(assay)
                                      & eligible["donor_set"].eq(donor_set), "cell_type"].astype(str)
            for cell_type in sorted(cell_types):
                indices = selected_units(cache, eligible, assay, donor_set, cell_type)
                frame = unit_frame(cache, indices)
                counts = cache["raw_counts"][indices].astype(float)
                totals = counts.sum(axis=1)
                expression = np.log1p(counts * (1e6 / totals[:, None]))
                gene_filter = (counts.sum(axis=0) >= 10) & ((counts > 0).sum(axis=0) >= 2)
                for model in models:
                    try:
                        design, names, omitted = build_design(frame, depth=model != "unadjusted_group_sensitivity")
                    except ValueError as exc:
                        statuses.append(status_row(assay, donor_set, model, "donor_pseudobulk_gene", cell_type, "untestable", str(exc)))
                        continue
                    fitted = fit_ols(expression[:, gene_filter], design)
                    testable = np.isfinite(fitted["pvalue"]) & np.isfinite(fitted["se"]) & (fitted["se"] > 0)
                    gene_indices = np.flatnonzero(gene_filter)[testable]
                    for position, gene_index in enumerate(gene_indices):
                        rows.append({
                            "assay": assay, "donor_set": donor_set, "model": model,
                            "cell_type": cell_type, "gene_id": cache["gene_ids"][gene_index],
                            "gene": cache["gene_names"][gene_index],
                            "effect_old_minus_young": fitted["effect"][testable][position],
                            "se": fitted["se"][testable][position],
                            "ci95_low": fitted["ci95_low"][testable][position],
                            "ci95_high": fitted["ci95_high"][testable][position],
                            "pvalue": fitted["pvalue"][testable][position],
                            "r2": fitted["r2"][testable][position],
                            "residual_df": fitted["residual_df"][testable][position],
                            "n_donors": len(frame), "n_young": int((frame["old"] == 0).sum()),
                            "n_old": int((frame["old"] == 1).sum()), "n_cells": int(frame["n_cells"].sum()),
                            "design_columns": ";".join(names), "omitted_covariates": ";".join(omitted),
                            "effect_scale": "log1p_CPM_donor_summed_raw_counts",
                        })
                    statuses.append(status_row(assay, donor_set, model, "donor_pseudobulk_gene", cell_type, "completed",
                                               f"tested_{int(testable.sum())}_of_{int(gene_filter.sum())}_filtered_genes"))
                    print(f"Kedlian genes {assay} {donor_set} {model} {cell_type}: {int(testable.sum())}", flush=True)
        output = add_bh(pd.DataFrame(rows), ["assay", "donor_set", "model"])
        outputs[assay] = output
    return outputs, statuses


def prepare_hellinger(cache: dict, eligible: pd.DataFrame) -> tuple[list[dict], pd.DataFrame]:
    references, rows = [], []
    genes = cache["gene_names"].astype(str)
    technical = technical_gene_mask(genes)
    for assay in ["scRNA", "snRNA"]:
        for donor_set in MODELS:
            cell_types = eligible.loc[eligible["eligible"] & eligible["assay"].eq(assay)
                                      & eligible["donor_set"].eq(donor_set), "cell_type"].astype(str)
            for cell_type in sorted(cell_types):
                indices = selected_units(cache, eligible, assay, donor_set, cell_type, qc=True)
                base = {"assay": assay, "donor_set": donor_set, "cell_type": cell_type,
                        "n_donors_qc20": len(indices)}
                if not len(indices):
                    rows.append({**base, "hellinger_eligible": False, "reason": "failed_20_qc_cell_or_base_donor_floors"})
                    continue
                old = cache["old"][indices]
                young = indices[old == 0]
                base["n_young_reference_donors"] = len(young)
                if len(young) < 2:
                    rows.append({**base, "hellinger_eligible": False, "reason": "fewer_than_two_released_young_reference_donors"})
                    continue
                total_reference = cache["qc_counts"][young].sum(axis=0)
                total_reference[technical] = 0
                total_young_cells = int(cache["qc_n"][young].sum())
                scored = []
                for index in indices:
                    own = cache["qc_counts"][index].copy() if index in young else np.zeros(len(genes), dtype=np.int64)
                    own[technical] = 0
                    reference = total_reference - own
                    own_cells = int(cache["qc_n"][index]) if index in young else 0
                    reference_n = total_young_cells - own_cells
                    if reference_n < MIN_REFERENCE_CELLS or reference.sum() <= 0:
                        continue
                    scored.append((index, reference_n, reference / reference.sum()))
                scored_indices = np.array([x[0] for x in scored], dtype=int)
                n_young = int((cache["old"][scored_indices] == 0).sum()) if len(scored) else 0
                n_old = int((cache["old"][scored_indices] == 1).sum()) if len(scored) else 0
                passed = len(scored) >= MIN_DONORS and n_young >= MIN_PER_GROUP and n_old >= MIN_PER_GROUP
                rows.append({**base, "hellinger_eligible": passed,
                             "reason": "pass" if passed else "failed_post_crossfit_six_donor_two_per_group_reference_floors",
                             "n_donors_after_reference_floors": len(scored),
                             "n_young_after_reference_floors": n_young, "n_old_after_reference_floors": n_old,
                             "min_cross_fitted_reference_cells": min([x[1] for x in scored], default=0),
                             "analysis_label": "NONNUMERIC_AGE_GROUP_SENSITIVITY"})
                if passed:
                    for index, reference_n, reference in scored:
                        references.append({"assay": assay, "donor_set": donor_set, "cell_type": cell_type,
                                           "unit_index": int(index), "reference_n": int(reference_n),
                                           "reference": reference})
    return references, pd.DataFrame(rows)


def run_hellinger(cache: dict, eligible: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[dict]]:
    references, eligibility_rows = prepare_hellinger(cache, eligible)
    statuses = []
    if not references:
        return pd.DataFrame(), pd.DataFrame(), eligibility_rows, statuses
    adata = ad.read_h5ad(H5AD, backed="r")
    try:
        n_refs = len(references)
        sums = np.zeros(n_refs)
        numbers = np.zeros(n_refs, dtype=np.int64)
        log_umi = np.zeros(n_refs)
        actual_umi = np.zeros(n_refs)
        for donor_set in MODELS:
            ref_indices = [i for i, ref in enumerate(references) if ref["donor_set"] == donor_set]
            if not ref_indices:
                continue
            unit_to_local = np.full(len(cache["donor_id"]), -1, dtype=np.int64)
            for local, ref_index in enumerate(ref_indices):
                unit_to_local[references[ref_index]["unit_index"]] = local
            refs = [references[i] for i in ref_indices]
            local_sums = np.zeros(len(refs))
            local_numbers = np.zeros(len(refs), dtype=np.int64)
            local_log = np.zeros(len(refs))
            local_actual = np.zeros(len(refs))
            counts_layer = adata.layers["counts"]
            for start in range(0, adata.n_obs, CHUNK):
                stop = min(start + CHUNK, adata.n_obs)
                units = cache["obs_unit"][start:stop]
                local = np.full(stop - start, -1, dtype=np.int64)
                valid_units = units >= 0
                local[valid_units] = unit_to_local[units[valid_units]]
                chosen = np.flatnonzero(local >= 0)
                if not len(chosen):
                    continue
                matrix = sp.csr_matrix(counts_layer[start:stop])[chosen].astype(np.float64)
                actual = np.asarray(matrix.sum(axis=1)).ravel()
                keep_genes = ~technical_gene_mask(cache["gene_names"].astype(str))
                nontechnical = np.asarray(matrix[:, keep_genes].sum(axis=1)).ravel()
                valid = (actual >= MIN_UMI_HELLINGER) & (nontechnical > 0)
                if not valid.any():
                    continue
                matrix = matrix[valid]
                actual = actual[valid]
                nontechnical = nontechnical[valid]
                local_codes = local[chosen][valid]
                for code in np.unique(local_codes):
                    choose = local_codes == code
                    distances = hellinger_distances(matrix[choose], refs[code]["reference"], nontechnical[choose])
                    local_sums[code] += float(distances.sum())
                    local_numbers[code] += int(choose.sum())
                    local_log[code] += float(np.log(actual[choose]).sum())
                    local_actual[code] += float(actual[choose].sum())
            for local, ref_index in enumerate(ref_indices):
                sums[ref_index] = local_sums[local]
                numbers[ref_index] = local_numbers[local]
                log_umi[ref_index] = local_log[local]
                actual_umi[ref_index] = local_actual[local]
    finally:
        adata.file.close()
    donor_rows = []
    for i, ref in enumerate(references):
        index = ref["unit_index"]
        if numbers[i] != cache["qc_n"][index] or numbers[i] < MIN_CELLS:
            raise ValueError(f"Hellinger cell count mismatch for {ref}")
        donor_rows.append({
            "assay": ref["assay"], "donor_set": ref["donor_set"], "cell_type": ref["cell_type"],
            "donor_id": str(cache["donor_id"][index]), "age_band": str(cache["age_band"][index]),
            "group": str(cache["group"][index]), "old": int(cache["old"][index]),
            "male": int(cache["male"][index]), "dcd": int(cache["dcd"][index]),
            "n_cells": int(numbers[i]), "mean_log_umi": log_umi[i] / numbers[i],
            "mean_umi": actual_umi[i] / numbers[i], "cross_fitted_reference_cells": ref["reference_n"],
            "hellinger": sums[i] / numbers[i], "analysis_label": "NONNUMERIC_AGE_GROUP_SENSITIVITY",
        })
    donors = pd.DataFrame(donor_rows)
    model_rows = []
    for (assay, donor_set, cell_type), frame in donors.groupby(["assay", "donor_set", "cell_type"]):
        for model in MODELS[donor_set]:
            try:
                fitted, _ = fit_frame(frame, frame["hellinger"].to_numpy(), model)
            except ValueError as exc:
                statuses.append(status_row(assay, donor_set, model, "Hellinger_nonnumeric_age_sensitivity", cell_type, "untestable", str(exc)))
                continue
            model_rows.append({"assay": assay, "donor_set": donor_set, "model": model,
                               "endpoint": "Hellinger_nonnumeric_age_sensitivity", "cell_type": cell_type,
                               "analysis_label": "NONNUMERIC_AGE_GROUP_SENSITIVITY", **fitted})
            statuses.append(status_row(assay, donor_set, model, "Hellinger_nonnumeric_age_sensitivity", cell_type, "completed", "pass"))
    models = add_bh(pd.DataFrame(model_rows), ["assay", "donor_set", "model"])
    return donors, models, eligibility_rows, statuses


def analyze() -> None:
    started = time.time()
    cache = load_cache()
    eligible = pd.read_csv(RESULTS / f"{PREFIX}-celltype-eligibility.csv")
    composition = pd.read_csv(RESULTS / f"{PREFIX}-composition-donor-counts.csv")
    p1_donors, p1_models, p1_null, p1_status = run_p1(cache, eligible)
    write_csv(p1_donors, f"{PREFIX}-p1-donor-scores.csv")
    write_csv(p1_models, f"{PREFIX}-p1-models.csv")
    write_csv(p1_null, f"{PREFIX}-p1-null.csv")
    comp_models, comp_status = run_composition(cache, composition, eligible)
    write_csv(comp_models, f"{PREFIX}-composition-models.csv")
    gene_models, gene_status = run_gene_models(cache, eligible)
    for assay, frame in gene_models.items():
        write_csv(frame, f"{PREFIX}-{assay}-gene-models.csv")
    h_donors, h_models, h_eligible, h_status = run_hellinger(cache, eligible)
    write_csv(h_donors, f"{PREFIX}-hellinger-donors.csv")
    write_csv(h_models, f"{PREFIX}-hellinger-models.csv")
    write_csv(h_eligible, f"{PREFIX}-hellinger-eligibility.csv")
    write_csv(pd.DataFrame(p1_status + comp_status + gene_status + h_status), f"{PREFIX}-model-status.csv")
    resources = {
        "elapsed_seconds": time.time() - started,
        "peak_rss_bytes": int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024),
        "memory_limit_gb": 6, "minimum_free_gb": 30,
        "backed_Hellinger_reads": True,
    }
    (RESULTS / f"{PREFIX}-analysis-resource.json").write_text(json.dumps(resources, indent=2, sort_keys=True) + "\n")
    print("Kedlian eligible outcomes completed; assays and multiplicity families remain separate.", flush=True)


def finalize() -> None:
    manifest = json.loads(MANIFEST.read_text())
    verification = json.loads((RESULTS / f"{PREFIX}-verification.json").read_text())
    entry = manifest["files"]["human_raw_h5ad"]
    entry.update({
        "downloaded": True, "local_path": str(Path("..") / H5AD.relative_to(REPO.parent)),
        "size_bytes": H5AD.stat().st_size, "sha256": sha256(H5AD),
        "guarded_download": {"minimum_free_gb": 30, "maximum_file_gb": 5},
        "integer_valued_counts_verified": True, "expression_values_read": True,
        "count_validation": f"results/{PREFIX}-verification.json",
    })
    manifest["scope"].update({"metadata_and_access_only": False, "expression_matrices_downloaded": True,
                              "expression_values_read": True, "biological_outcomes_examined": True})
    manifest["analysis"] = {
        "completed": "2026-11-14", "frozen_plan": "data/manifest/muscle_analysis_plan.json",
        "script": "scripts/run_followup_kedlian.py", "report": "reports/followup-muscle-analysis.md",
        "assays_pooled": False, "inferential_unit": "donor", "numeric_ages_used": False,
        "primary_excluded_donor": "470BR", "all_donor_sensitivity_predeclared": True,
        "p1_random_null_iterations": NULL_ITERATIONS, "p1_random_seed": NULL_SEED,
        "gene_model_multiplicity": "assay-wide by donor set and model",
        "Hellinger": "NONNUMERIC_AGE_GROUP_SENSITIVITY with cross-fitted released-young reference",
        "peak_rss_bytes_observed_in_cache_job": verification["peak_rss_bytes"],
        "raw_fastq_downloaded": False,
    }
    manifest["verdict"].update({"expression_analysis_performed": True,
                                "raw_count_value_validation": "passed_all_layers_counts_float32_values_finite_nonnegative_integer_valued",
                                "young_old_group_analysis": "completed_eligible_assay_specific_donor_models"})
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["inspect", "build-cache", "analyze", "finalize"])
    args = parser.parse_args()
    RESULTS.mkdir(exist_ok=True)
    if args.command == "inspect":
        inspect_object()
    elif args.command == "build-cache":
        build_cache()
    elif args.command == "analyze":
        analyze()
    else:
        finalize()


if __name__ == "__main__":
    main()
