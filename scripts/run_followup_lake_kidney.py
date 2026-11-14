#!/usr/bin/env python3
"""Restricted public Lake kidney analysis using per-assay CELLxGENE H5ADs.

Run one assay at a time through scripts/resource_guard.py. The first pass
validates all of raw.X and checkpoints donor aggregates; outcome models run
only after that gate succeeds. Hellinger uses a second backed pass only for
cell types passing the frozen young-reference and cross-fitting floors.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats as stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from src.aging_programs import load_program_tiers  # noqa: E402

RESULTS = REPO / "results"
PREFIX = "followup-lake-kidney"
DEFAULT_DATA = Path("/home/reza/projects/agescape-revision-data/lake_kidney")
CHUNK = 500
MIN_DONORS = 6
MIN_CELLS = 20
MIN_REFERENCE_CELLS = 50
MIN_UMI_HELLINGER = 1500
NULL_ITERATIONS = 200
NULL_SEED = 20261105
HUMAN_MAX_LIFESPAN = 122.5
AGE_RE = re.compile(r"^(\d+)-year-old stage$")
UNKNOWN_LABELS = {"", "unknown", "unassigned", "nan", "none"}

RELEASES = {
    "scCv3": {
        "dataset_id": "32b9bdce-2481-4c85-ba1b-6ad5fcea844c",
        "version": "3ae68468-89f7-4941-a945-f7bae9b2cbea",
        "filename": "lake_scCv3_3ae68468.h5ad",
        "url": "https://datasets.cellxgene.cziscience.com/3ae68468-89f7-4941-a945-f7bae9b2cbea.h5ad",
        "cells": 107344,
        "assay": "10x 3' v3",
        "suspension": "cell",
        "geo_file": "GSE183276_scCv3_Metadata_03282022.txt.gz",
        "geo_accession": "GSE183276",
        "expected_geo": (18, 20, 20),
        "prior_cohort": "Kidney_SCR_v15",
    },
    "snCv3": {
        "dataset_id": "07854d9c-5375-4a9b-ac34-fa919d3c3686",
        "version": "5dc03d39-ac5f-49e1-8630-a64d2e6e293b",
        "filename": "lake_snCv3_5dc03d39.h5ad",
        "url": "https://datasets.cellxgene.cziscience.com/5dc03d39-ac5f-49e1-8630-a64d2e6e293b.h5ad",
        "cells": 172847,
        "assay": "10x 3' v3",
        "suspension": "nucleus",
        "geo_file": "GSE183277_snCv3_Metadata_03282022.txt.gz",
        "geo_accession": "GSE183277",
        "expected_geo": (13, 18, 21),
        "prior_cohort": "Kidney_SN_v15",
    },
}

SC_AGES = {
    "Sample1153-EO1": 35, "Sample1153-EO2": 35, "Sample1153-EO3": 30,
    "Sample1157-EO1": 49, "Sample1157-EO2": 41, "Sample1157-EO3": 47,
    "Sample1158-EO1": 48, "Sample1158-EO2": 56, "Sample1158-EO3": 35,
    "Sample1162-EO1": 39, "Sample1162-EO2": 52, "PRE018-1": 56,
    "PRE019": 53, "PRE027": 31, "PRE038": 42, "PRE055-1": 66,
    "PRE062-1": 40, "PRE98sc": 57,
}
SN_AGES = {
    "3490": 14, "3499": 57, "3504": 45, "3535": 76, "3593": 66,
    "3613": 59, "KRP446": 29, "KRP460": 63, "KRP461": 54, "KRP462": 69,
}
SN_DECADE_ONLY = {"18-142", "18-162", "18-312"}
SC_GEO_IDS = {
    "Sample1153-EO1": "1153EO1", "Sample1153-EO2": "1153EO2",
    "Sample1153-EO3": "1153EO3", "Sample1157-EO1": "1157EO1",
    "Sample1157-EO2": "1157EO2", "Sample1157-EO3": "1157EO3",
    "Sample1158-EO1": "1158EO1", "Sample1158-EO2": "1158EO2",
    "Sample1158-EO3": "1158EO3", "Sample1162-EO1": "1162EO1",
    "Sample1162-EO2": "1162EO2", "PRE018-1": "PRE018", "PRE019": "PRE019",
    "PRE027": "PRE027", "PRE038": "PRE038", "PRE055-1": "PRE055",
    "PRE062-1": "PRE062", "PRE98sc": "PRE098",
}


def expected_ages(assay: str) -> dict[str, int]:
    return dict(SC_AGES if assay == "scCv3" else SN_AGES)


def expected_healthy(assay: str) -> set[str]:
    return set(expected_ages(assay)) | (SN_DECADE_ONLY if assay == "snCv3" else set())


def geo_donor_id(assay: str, donor: str) -> str:
    return SC_GEO_IDS.get(donor, donor) if assay == "scCv3" else donor


def parse_exact_age(value: object) -> float:
    match = AGE_RE.fullmatch(str(value))
    return float(match.group(1)) if match else np.nan


def known_cell_type(values: pd.Series) -> np.ndarray:
    return ~values.astype(str).str.strip().str.lower().isin(UNKNOWN_LABELS).to_numpy()


def clr_transform(counts: pd.DataFrame, pseudocount: float = 0.5) -> pd.DataFrame:
    if pseudocount <= 0:
        raise ValueError("CLR pseudocount must be positive")
    adjusted = counts.astype(float) + pseudocount
    log_props = np.log(adjusted.div(adjusted.sum(axis=1), axis=0))
    result = log_props.sub(log_props.mean(axis=1), axis=0)
    if not np.isfinite(result.to_numpy()).all():
        raise ValueError("Non-finite CLR value")
    if not np.allclose(result.sum(axis=1), 0.0, atol=1e-10):
        raise AssertionError("CLR rows must sum to zero")
    return result


def technical_gene_mask(gene_names: np.ndarray) -> np.ndarray:
    genes = pd.Index(gene_names).astype(str).str.upper()
    return np.array([
        bool(re.match(r"^(MT-|RP[LS]|MRP[LS])", gene))
        or gene in {"XIST", "TSIX", "DDX3Y", "EIF2S3Y", "UTY", "KDM5D"}
        for gene in genes
    ])


def fit_donor_model(data: pd.DataFrame, outcome: str, adjustment: str = "none") -> dict | None:
    """One-row-per-donor OLS; age is always the second design column."""
    required = [outcome, "age", "sex"]
    frame = data.dropna(subset=required).copy()
    if len(frame) < MIN_DONORS or frame["age"].nunique() < 2:
        return None
    columns = [np.ones(len(frame)), frame["age"].to_numpy(float)]
    names = ["intercept", "age"]
    male = (frame["sex"].astype(str).str.lower() == "male").to_numpy(float)
    if len(np.unique(male)) == 2:
        columns.append(male)
        names.append("male")
    if adjustment == "mean_log_umi":
        columns.append(frame["mean_log_umi"].to_numpy(float))
        names.append("mean_log_umi")
    elif adjustment == "log_total_cells":
        columns.append(np.log(frame["total_cells"].to_numpy(float)))
        names.append("log_total_cells")
    elif adjustment != "none":
        raise ValueError(f"Unknown adjustment: {adjustment}")
    X = np.column_stack(columns)
    n, k = X.shape
    if n <= k + 1 or np.linalg.matrix_rank(X) < k or not np.isfinite(X).all():
        return None
    y = frame[outcome].to_numpy(float)
    fit = sm.OLS(y, X).fit()
    ci = fit.conf_int()[1]
    slope = float(fit.params[1])
    return {
        "n_donors": int(n), "n_cells": int(frame["n_cells"].sum()) if "n_cells" in frame else np.nan,
        "age_min": float(frame["age"].min()), "age_max": float(frame["age"].max()),
        "age_sd": float(frame["age"].std(ddof=1)), "sex_covariate": "male" in names,
        "adjustment": adjustment, "slope_per_year": slope, "se": float(fit.bse[1]),
        "ci95_low": float(ci[0]), "ci95_high": float(ci[1]), "pvalue": float(fit.pvalues[1]),
        "r2": float(fit.rsquared), "effect_per_30y": slope * 30.0,
        "norm_slope_30y_max_lifespan": slope * 30.0 / HUMAN_MAX_LIFESPAN,
    }


def add_bh(frame: pd.DataFrame, group_columns: list[str], p_column: str = "pvalue",
           q_column: str = "q_bh") -> pd.DataFrame:
    frame = frame.copy()
    frame[q_column] = np.nan
    if frame.empty:
        return frame
    for _, indices in frame.groupby(group_columns, dropna=False).groups.items():
        valid = frame.loc[indices, p_column].dropna().index
        if len(valid):
            frame.loc[valid, q_column] = multipletests(
                frame.loc[valid, p_column], method="fdr_bh"
            )[1]
    return frame


def random_signature_scores(matrix: np.ndarray, p1_up: np.ndarray, p1_down: np.ndarray,
                            iterations: int, rng: np.random.Generator) -> np.ndarray:
    """Expression-decile-matched random signatures; genes x donors matrix."""
    means = np.asarray(matrix.mean(axis=1)).ravel()
    deciles = pd.Series(pd.qcut(means, 10, labels=False, duplicates="drop")).fillna(-1).to_numpy()
    available = np.ones(matrix.shape[0], dtype=bool)
    available[np.concatenate([p1_up, p1_down])] = False
    scores = []
    for _ in range(iterations):
        directions = []
        for program_genes in (p1_up, p1_down):
            selected = []
            for decile, count in pd.Series(deciles[program_genes]).value_counts().items():
                pool = np.flatnonzero((deciles == decile) & available)
                if len(pool) < int(count):
                    raise ValueError("Insufficient expression-matched random-signature pool")
                selected.extend(rng.choice(pool, size=int(count), replace=False))
            directions.append(selected)
        scores.append(matrix[directions[0], :].mean(axis=0) - matrix[directions[1], :].mean(axis=0))
    return np.vstack(scores)


def fast_age_sex_slope(score: np.ndarray, age: np.ndarray, sex: np.ndarray) -> float:
    columns = [np.ones(len(age)), np.asarray(age, float)]
    male = (pd.Series(sex).astype(str).str.lower() == "male").to_numpy(float)
    if len(np.unique(male)) == 2:
        columns.append(male)
    X = np.column_stack(columns)
    return float(np.linalg.lstsq(X, np.asarray(score, float), rcond=None)[0][1])


def upsert_csv(path: Path, new: pd.DataFrame, key_columns: list[str]) -> None:
    if path.exists():
        old = pd.read_csv(path)
        if not old.empty:
            old = old[~old.set_index(key_columns).index.isin(new.set_index(key_columns).index)]
        combined = pd.concat([old, new], ignore_index=True, sort=False)
    else:
        combined = new
    combined.to_csv(path, index=False)


def read_release_obs(assay: str, path: Path):
    import anndata as ad

    adata = ad.read_h5ad(path, backed="r")
    required = ["donor_id", "development_stage", "disease", "assay", "suspension_type",
                "cell_type", "sex", "is_primary_data"]
    missing = [column for column in required if column not in adata.obs.columns]
    if missing:
        raise ValueError(f"Missing obs columns: {missing}")
    obs = adata.obs[required].copy()
    for column in required:
        obs[column] = obs[column].astype(str)
    obs["age"] = obs["development_stage"].map(parse_exact_age)
    obs["_row"] = np.arange(len(obs))
    return adata, obs


def geo_crosswalk(assay: str, data_dir: Path) -> pd.DataFrame:
    release = RELEASES[assay]
    columns = ["patient", "specimen", "library", "condition.l1", "sex", "assay"]
    geo = pd.read_csv(data_dir / release["geo_file"], sep="\t", usecols=columns, dtype=str)
    geo = geo[geo["condition.l1"] == "Ref"].copy()
    geo["cellxgene_donor_id"] = geo["patient"].map(
        {geo_donor_id(assay, donor): donor for donor in expected_healthy(assay)}
    )
    if geo["cellxgene_donor_id"].isna().any():
        raise ValueError(f"Unmapped GEO reference patients: {sorted(geo.loc[geo['cellxgene_donor_id'].isna(), 'patient'].unique())}")
    rows = []
    for donor, part in geo.groupby("cellxgene_donor_id"):
        sex = set(part["sex"].dropna().astype(str))
        rows.append({
            "cellxgene_donor_id": donor, "geo_patient": str(part["patient"].iloc[0]),
            "geo_specimens": int(part["specimen"].nunique()),
            "geo_libraries": int(part["library"].nunique()),
            "geo_sex_values": ";".join(sorted(sex)),
            "geo_assay_values": ";".join(sorted(part["assay"].dropna().astype(str).unique())),
        })
    crosswalk = pd.DataFrame(rows)
    actual = (len(crosswalk), int(crosswalk["geo_specimens"].sum()),
              int(crosswalk["geo_libraries"].sum()))
    if actual != release["expected_geo"]:
        raise ValueError(f"GEO healthy gate mismatch: {actual} != {release['expected_geo']}")
    return crosswalk


def verify_obs(assay: str, adata, obs: pd.DataFrame, data_dir: Path) -> tuple[pd.DataFrame, dict]:
    release = RELEASES[assay]
    if adata.n_obs != release["cells"]:
        raise ValueError(f"Cell count {adata.n_obs} != released {release['cells']}")
    if set(obs["assay"]) != {release["assay"]}:
        raise ValueError(f"Assay identity mismatch: {set(obs['assay'])}")
    if set(obs["suspension_type"]) != {release["suspension"]}:
        raise ValueError(f"Suspension identity mismatch: {set(obs['suspension_type'])}")
    normal = obs[obs["disease"] == "normal"].copy()
    healthy = set(normal["donor_id"])
    if healthy != expected_healthy(assay):
        raise ValueError(f"Healthy donors differ from frozen mapping: {sorted(healthy)}")
    for donor, part in normal.groupby("donor_id"):
        for column in ["development_stage", "disease", "sex", "assay", "suspension_type"]:
            if part[column].nunique() != 1:
                raise ValueError(f"Non-unique {column} for donor {donor}")
    ages = expected_ages(assay)
    observed_ages = normal.groupby("donor_id")["age"].first().to_dict()
    for donor, age in ages.items():
        if not np.isfinite(observed_ages.get(donor, np.nan)) or int(observed_ages[donor]) != age:
            raise ValueError(f"Age mismatch for {donor}: {observed_ages.get(donor)} != {age}")
    decade = healthy - set(ages)
    if assay == "snCv3" and decade != SN_DECADE_ONLY:
        raise ValueError(f"Decade-only mapping mismatch: {decade}")
    if normal.loc[normal["donor_id"].isin(decade), "age"].notna().any():
        raise ValueError("Decade-only donor unexpectedly parsed as exact age")
    crosswalk = geo_crosswalk(assay, data_dir)
    donors = normal.groupby("donor_id").agg(
        development_stage=("development_stage", "first"), age=("age", "first"),
        disease=("disease", "first"), sex=("sex", "first"), assay=("assay", "first"),
        suspension_type=("suspension_type", "first"), n_cells=("donor_id", "size"),
    ).reset_index().rename(columns={"donor_id": "cellxgene_donor_id"})
    donors = donors.merge(crosswalk, on="cellxgene_donor_id", how="left", validate="one_to_one")
    donors.insert(0, "assay_release", assay)
    donors["geo_patient_expected"] = donors["cellxgene_donor_id"].map(
        lambda donor: geo_donor_id(assay, donor)
    )
    if not (donors["geo_patient"] == donors["geo_patient_expected"]).all():
        raise ValueError("GEO donor crosswalk mismatch")
    donors["exact_age"] = donors["age"].notna()
    donors["primary_adult"] = donors["age"].between(20, 120) & donors["exact_age"]
    donors["sensitivity_all_exact_age"] = donors["exact_age"]
    expected_primary = 18 if assay == "scCv3" else 9
    if int(donors["primary_adult"].sum()) != expected_primary:
        raise ValueError(f"Primary adult donor count != {expected_primary}")
    verification = {
        "assay_release": assay, "dataset_id": release["dataset_id"],
        "dataset_version_id": release["version"], "source_url": release["url"],
        "n_obs": int(adata.n_obs), "n_vars_X": int(adata.n_vars),
        "assay_labels": sorted(obs["assay"].unique()),
        "suspension_labels": sorted(obs["suspension_type"].unique()),
        "is_primary_data_labels": sorted(obs["is_primary_data"].unique()),
        "released_annotation_labels": sorted(obs["cell_type"].unique()),
        "healthy_donors": len(healthy), "exact_age_healthy_donors": len(ages),
        "decade_only_healthy_donors": len(decade),
        "primary_adult_donors": int(donors["primary_adult"].sum()),
        "geo_reference_donors_specimens_libraries": list(release["expected_geo"]),
        "geo_age_field_used": False,
        "donor_age_source": "H5AD development_stage, verified against frozen public mapping",
    }
    return donors, verification


def eligibility_table(obs: pd.DataFrame) -> pd.DataFrame:
    selected = obs[obs["disease"].eq("normal") & obs["age"].between(20, 120)
                   & known_cell_type(obs["cell_type"])].copy()
    counts = selected.groupby(["cell_type", "donor_id"]).size().rename("n_cells").reset_index()
    rows = []
    for cell_type, part in counts.groupby("cell_type"):
        eligible = part[part["n_cells"] >= MIN_CELLS]
        rows.append({
            "cell_type": cell_type, "n_donors_observed": int(len(part)),
            "n_donors_ge20": int(len(eligible)), "n_cells_total": int(part["n_cells"].sum()),
            "eligible": len(eligible) >= MIN_DONORS,
            "reason": "pass" if len(eligible) >= MIN_DONORS else "fewer_than_six_donors_with_20_cells",
        })
    return pd.DataFrame(rows).sort_values("cell_type").reset_index(drop=True)


def aggregate_block(matrix: sp.csr_matrix, codes: np.ndarray, output: np.ndarray,
                    positions: np.ndarray | None = None) -> None:
    unique = np.unique(codes)
    remap = np.searchsorted(unique, codes)
    indicator = sp.csr_matrix((
        np.ones(len(codes), dtype=np.float64), (np.arange(len(codes)), remap)
    ), shape=(len(codes), len(unique)))
    block = (indicator.T @ matrix).toarray()
    target = unique if positions is None else positions[unique]
    if np.issubdtype(output.dtype, np.integer):
        if not np.equal(block, np.floor(block)).all():
            raise ValueError("Non-integer aggregate assigned to integer count cache")
        block = block.astype(output.dtype)
    output[target] += block


def build_cache(assay: str, adata, obs: pd.DataFrame, donors: pd.DataFrame,
                eligible: pd.DataFrame, cache_path: Path) -> dict:
    if adata.raw is None:
        raise ValueError("raw.X is absent")
    raw = adata.raw
    genes = (raw.var["feature_name"].astype(str).to_numpy()
             if "feature_name" in raw.var.columns else raw.var_names.astype(str).to_numpy())
    if not raw.var_names.is_unique:
        raise ValueError("raw.var feature IDs are not unique")
    eligible_types = set(eligible.loc[eligible["eligible"], "cell_type"])
    selected_mask = (obs["disease"].eq("normal") & obs["age"].notna()
                     & obs["cell_type"].isin(eligible_types)).to_numpy()
    selected = obs.loc[selected_mask].copy()
    units = selected.groupby(["donor_id", "cell_type"], sort=True).agg(
        age=("age", "first"), sex=("sex", "first"), primary=("age", lambda x: bool(20 <= x.iloc[0] <= 120)),
    ).reset_index()
    unit_lookup = {(row.donor_id, row.cell_type): i for i, row in enumerate(units.itertuples())}
    codes = np.full(len(obs), -1, dtype=np.int64)
    codes[selected["_row"].to_numpy()] = [unit_lookup[(d, c)] for d, c in zip(selected["donor_id"], selected["cell_type"])]
    n_units, n_genes = len(units), len(genes)
    log_sums = np.zeros((n_units, n_genes), dtype=np.float64)
    qc_counts = np.zeros((n_units, n_genes), dtype=np.int64)
    n_cells = np.zeros(n_units, dtype=np.int64)
    qc_n = np.zeros(n_units, dtype=np.int64)
    sum_log_umi = np.zeros(n_units, dtype=np.float64)
    sum_umi = np.zeros(n_units, dtype=np.float64)
    primary_unit = units["primary"].to_numpy(bool)
    raw_nnz = 0
    raw_sum = 0.0
    X = raw.X
    for start in range(0, adata.n_obs, CHUNK):
        stop = min(start + CHUNK, adata.n_obs)
        matrix = sp.csr_matrix(X[start:stop])
        values = matrix.data
        raw_nnz += values.size
        raw_sum += float(values.sum(dtype=np.float64))
        if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.floor(values)).all():
            raise ValueError(f"raw.X failed finite nonnegative integer check at rows {start}:{stop}")
        local_codes = codes[start:stop]
        chosen = np.flatnonzero(local_codes >= 0)
        if not len(chosen):
            continue
        sub = matrix[chosen].tocsr().astype(np.float64)
        sub_codes = local_codes[chosen]
        libraries = np.asarray(sub.sum(axis=1)).ravel()
        if (libraries <= 0).any():
            raise ValueError("Selected healthy exact-age cell has zero raw library size")
        np.add.at(n_cells, sub_codes, 1)
        np.add.at(sum_log_umi, sub_codes, np.log(libraries))
        np.add.at(sum_umi, sub_codes, libraries)
        log_sub = sub.copy()
        scale = (1e6 / libraries).astype(np.float64)
        log_sub.data *= np.repeat(scale, np.diff(log_sub.indptr))
        np.log1p(log_sub.data, out=log_sub.data)
        aggregate_block(log_sub, sub_codes, log_sums)
        qc = primary_unit[sub_codes] & (libraries >= MIN_UMI_HELLINGER)
        if qc.any():
            aggregate_block(sub[qc], sub_codes[qc], qc_counts)
            np.add.at(qc_n, sub_codes[qc], 1)
        if start % 10000 == 0:
            print(f"[{assay}] raw.X rows {start}:{stop}; {raw_nnz:,} stored values checked", flush=True)
    cache = {
        "genes": genes, "donor_id": units["donor_id"].astype(str).to_numpy(),
        "cell_type": units["cell_type"].astype(str).to_numpy(),
        "age": units["age"].to_numpy(float), "sex": units["sex"].astype(str).to_numpy(),
        "primary": primary_unit, "n_cells": n_cells, "qc_n": qc_n,
        "mean_log_umi": sum_log_umi / n_cells, "sum_umi": sum_umi,
        "log_sums": log_sums.astype(np.float32), "qc_counts": qc_counts,
    }
    np.savez(cache_path, **cache)
    verification = {
        "raw_X_present": True, "raw_X_shape": [int(adata.n_obs), int(n_genes)],
        "raw_X_nnz": int(raw_nnz), "raw_X_sum": raw_sum,
        "raw_X_all_stored_values_finite_nonnegative_integer": True,
        "raw_feature_ids_unique": True,
        "raw_unique_feature_names": len(set(genes)),
        "raw_duplicate_feature_name_rows": int(pd.Index(genes).duplicated(keep=False).sum()),
        "aggregation_units": n_units,
        "cache_path": str(cache_path), "cache_bytes": cache_path.stat().st_size,
    }
    return {"arrays": cache, "verification": verification}


def load_cache(path: Path) -> dict:
    with np.load(path, allow_pickle=False) as data:
        return {key: data[key] for key in data.files}


def run_p1(assay: str, cache: dict, eligible: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    program = load_program_tiers(REPO, "human")["P1_Clean"]
    genes = cache["genes"].astype(str)
    up = np.flatnonzero(np.isin(genes, program["up"]))
    down = np.flatnonzero(np.isin(genes, program["down"]))
    if not len(up) or not len(down):
        raise ValueError("P1 lacks detected up or down genes")
    donor_frames, model_rows = [], []
    rng = np.random.default_rng(NULL_SEED)
    for cell_type in sorted(eligible.loc[eligible["eligible"], "cell_type"]):
        indices = np.flatnonzero(cache["cell_type"].astype(str) == cell_type)
        frames = {}
        for donor_set, primary_only in [("primary_adult", True), ("all_exact_age", False)]:
            keep = indices[(cache["n_cells"][indices] >= MIN_CELLS)
                           & (cache["primary"][indices] | (not primary_only))]
            order = keep[np.argsort(cache["donor_id"][keep].astype(str))]
            matrix = (cache["log_sums"][order] / cache["n_cells"][order, None]).astype(float).T
            frame = pd.DataFrame({
                "donor_id": cache["donor_id"][order].astype(str), "cell_type": cell_type,
                "age": cache["age"][order], "sex": cache["sex"][order].astype(str),
                "n_cells": cache["n_cells"][order],
                "mean_log_umi": cache["mean_log_umi"][order],
                "score": matrix[up, :].mean(axis=0) - matrix[down, :].mean(axis=0),
            })
            frames[donor_set] = (frame, matrix, order)
            frame.insert(0, "assay_release", assay)
            frame["donor_set"] = donor_set
            donor_frames.append(frame)
            for model, adjustment in [("primary_age_sex", "none"), ("depth_sensitivity", "mean_log_umi")]:
                if donor_set == "all_exact_age" and model != "primary_age_sex":
                    continue
                fit = fit_donor_model(frame, "score", adjustment)
                if fit is None:
                    raise ValueError(f"Eligible P1 model unidentifiable: {assay} {cell_type} {donor_set} {model}")
                model_rows.append({
                    "assay_release": assay, "endpoint": "P1_transport", "cell_type": cell_type,
                    "donor_set": donor_set, "model": model, "p1_up_detected": len(up),
                    "p1_down_detected": len(down), "p1_coverage": (len(up) + len(down)) /
                    (len(program["up"]) + len(program["down"])), **fit,
                })
        primary_frame, matrix, _ = frames["primary_adult"]
        target = next(row for row in reversed(model_rows)
                      if row["cell_type"] == cell_type and row["donor_set"] == "primary_adult"
                      and row["model"] == "primary_age_sex")
        observed = target["slope_per_year"]
        random_scores = random_signature_scores(matrix, up, down, NULL_ITERATIONS, rng)
        random_slopes = np.array([
            fast_age_sex_slope(score, primary_frame["age"].to_numpy(), primary_frame["sex"].to_numpy())
            for score in random_scores
        ])
        null_mean = float(random_slopes.mean())
        empirical_p = (1 + int(np.sum(np.abs(random_slopes - null_mean) >= abs(observed - null_mean)))) / (NULL_ITERATIONS + 1)
        target.update({
            "random_iterations": NULL_ITERATIONS, "random_seed": NULL_SEED,
            "random_mean_slope": null_mean, "random_sd": float(random_slopes.std(ddof=1)),
            "random_q025": float(np.quantile(random_slopes, 0.025)),
            "random_q975": float(np.quantile(random_slopes, 0.975)),
            "random_pvalue": float(empirical_p),
        })
        print(f"[{assay}] P1 {cell_type}: n={len(primary_frame)} null p={empirical_p:.4f}", flush=True)
    models = pd.DataFrame(model_rows)
    models = add_bh(models, ["assay_release", "endpoint", "donor_set", "model"])
    primary = models["donor_set"].eq("primary_adult") & models["model"].eq("primary_age_sex")
    models["random_q_bh"] = np.nan
    if primary.any():
        models.loc[primary, "random_q_bh"] = multipletests(
            models.loc[primary, "random_pvalue"], method="fdr_bh"
        )[1]
    return pd.concat(donor_frames, ignore_index=True), models


def composition_counts(obs: pd.DataFrame, assay: str) -> pd.DataFrame:
    selected = obs[obs["disease"].eq("normal") & obs["age"].notna()
                   & known_cell_type(obs["cell_type"])].copy()
    counts = selected.groupby(["donor_id", "cell_type"]).size().rename("n_cells").reset_index()
    meta = selected.groupby("donor_id").agg(age=("age", "first"), sex=("sex", "first")).reset_index()
    all_types = sorted(selected["cell_type"].unique())
    wide = counts.pivot(index="donor_id", columns="cell_type", values="n_cells").fillna(0)
    wide = wide.reindex(index=meta["donor_id"], columns=all_types, fill_value=0).astype(int)
    long = wide.stack().rename("n_cells").reset_index()
    long = long.merge(meta, on="donor_id", validate="many_to_one")
    long.insert(0, "assay_release", assay)
    long["primary_adult"] = long["age"].between(20, 120)
    long["total_cells"] = long.groupby("donor_id")["n_cells"].transform("sum")
    return long


def run_composition(assay: str, counts: pd.DataFrame, eligible: pd.DataFrame) -> pd.DataFrame:
    rows = []
    eligible_types = sorted(eligible.loc[eligible["eligible"], "cell_type"])
    for donor_set, primary_only in [("primary_adult", True), ("all_exact_age", False)]:
        part = counts[counts["primary_adult"] | (not primary_only)].copy()
        wide = part.pivot(index="donor_id", columns="cell_type", values="n_cells").fillna(0)
        meta = part.groupby("donor_id").agg(age=("age", "first"), sex=("sex", "first"),
                                             total_cells=("total_cells", "first"))
        clr = clr_transform(wide)
        for cell_type in eligible_types:
            frame = meta.copy()
            frame["clr"] = clr[cell_type]
            frame["n_cells"] = wide[cell_type]
            for model, adjustment in [("primary_age_sex", "none"), ("total_cells_sensitivity", "log_total_cells")]:
                if donor_set == "all_exact_age" and model != "primary_age_sex":
                    continue
                fit = fit_donor_model(frame.reset_index(), "clr", adjustment)
                if fit is None:
                    raise ValueError(f"Eligible CLR model unidentifiable: {assay} {cell_type} {donor_set} {model}")
                rows.append({"assay_release": assay, "endpoint": "CLR_composition",
                             "cell_type": cell_type, "donor_set": donor_set, "model": model, **fit})
    models = add_bh(pd.DataFrame(rows), ["assay_release", "endpoint", "donor_set", "model"])
    return models


def hellinger_distances(matrix: sp.csr_matrix, reference_props: np.ndarray,
                        nontechnical_totals: np.ndarray) -> np.ndarray:
    root = matrix.copy().astype(np.float64)
    root.data = np.sqrt(root.data)
    affinity = np.asarray(root @ np.sqrt(reference_props)).ravel() / np.sqrt(nontechnical_totals)
    return np.sqrt(np.maximum(0.0, 1.0 - affinity))


def run_hellinger(assay: str, adata, obs: pd.DataFrame, cache: dict,
                  eligible: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    genes = cache["genes"].astype(str)
    technical = technical_gene_mask(genes)
    keep_genes = ~technical
    donor_rows, eligibility_rows = [], []
    primary = cache["primary"].astype(bool)
    for cell_type in sorted(eligible.loc[eligible["eligible"], "cell_type"]):
        indices = np.flatnonzero((cache["cell_type"].astype(str) == cell_type) & primary)
        indices = indices[cache["qc_n"][indices] >= MIN_CELLS]
        young = indices[cache["age"][indices] <= 40]
        base = {"assay_release": assay, "cell_type": cell_type,
                "n_donors_qc20": len(indices), "n_young_reference_donors": len(young)}
        if len(young) < 2:
            eligibility_rows.append({**base, "hellinger_eligible": False,
                                     "reason": "fewer_than_two_young_reference_donors_with_20_qc_cells"})
            continue
        total_reference = cache["qc_counts"][young].sum(axis=0)
        total_reference[technical] = 0
        total_young_cells = int(cache["qc_n"][young].sum())
        scored = []
        references = {}
        for index in indices:
            donor = str(cache["donor_id"][index])
            own = cache["qc_counts"][index].copy() if index in young else np.zeros(len(genes), dtype=np.int64)
            own[technical] = 0
            reference = total_reference - own
            own_cells = int(cache["qc_n"][index]) if index in young else 0
            reference_n = total_young_cells - own_cells
            if reference_n < MIN_REFERENCE_CELLS or reference.sum() <= 0:
                continue
            references[donor] = reference / reference.sum()
            scored.append((index, donor, reference_n))
        if len(scored) < MIN_DONORS:
            eligibility_rows.append({**base, "hellinger_eligible": False,
                                     "reason": "fewer_than_six_donors_after_cross_fitted_reference_floors",
                                     "n_donors_after_reference_floors": len(scored)})
            continue
        eligibility_rows.append({**base, "hellinger_eligible": True, "reason": "pass",
                                 "n_donors_after_reference_floors": len(scored),
                                 "min_cross_fitted_reference_cells": min(x[2] for x in scored)})
        scored_donors = {donor for _, donor, _ in scored}
        row_donor = obs["donor_id"].to_numpy(str)
        wanted = (obs["donor_id"].isin(scored_donors) & obs["cell_type"].eq(cell_type)
                  & obs["disease"].eq("normal") & obs["age"].between(20, 120)).to_numpy()
        global_rows = np.flatnonzero(wanted)
        accum = {donor: {"sum": 0.0, "n": 0, "log_umi": 0.0, "actual_umi": 0.0}
                 for donor in scored_donors}
        X = adata.raw.X
        for offset in range(0, len(global_rows), CHUNK):
            rows = global_rows[offset:offset + CHUNK]
            matrix = sp.csr_matrix(X[rows]).astype(np.float64)
            actual_lib = np.asarray(matrix.sum(axis=1)).ravel()
            nontechnical = np.asarray(matrix[:, keep_genes].sum(axis=1)).ravel()
            valid = (actual_lib >= MIN_UMI_HELLINGER) & (nontechnical > 0)
            if not valid.any():
                continue
            matrix = matrix[valid]
            actual_lib = actual_lib[valid]
            nontechnical = nontechnical[valid]
            donors_chunk = row_donor[rows][valid]
            for donor in np.unique(donors_chunk):
                choose = donors_chunk == donor
                distances = hellinger_distances(matrix[choose], references[str(donor)], nontechnical[choose])
                accum[str(donor)]["sum"] += float(distances.sum())
                accum[str(donor)]["n"] += int(choose.sum())
                accum[str(donor)]["log_umi"] += float(np.log(actual_lib[choose]).sum())
                accum[str(donor)]["actual_umi"] += float(actual_lib[choose].sum())
        for index, donor, reference_n in scored:
            values = accum[donor]
            if values["n"] != cache["qc_n"][index] or values["n"] < MIN_CELLS:
                raise ValueError(f"Hellinger cell count mismatch for {assay} {cell_type} {donor}")
            donor_rows.append({
                "assay_release": assay, "cell_type": cell_type, "donor_id": donor,
                "age": float(cache["age"][index]), "sex": str(cache["sex"][index]),
                "n_cells": int(values["n"]), "mean_log_umi": values["log_umi"] / values["n"],
                "mean_umi": values["actual_umi"] / values["n"],
                "cross_fitted_reference_cells": reference_n,
                "hellinger": values["sum"] / values["n"],
            })
        print(f"[{assay}] Hellinger {cell_type}: {len(scored)} donors", flush=True)
    donors = pd.DataFrame(donor_rows)
    model_rows = []
    if not donors.empty:
        for (ct,), part in donors.groupby(["cell_type"]):
            for model, adjustment in [("primary_depth_adjusted", "mean_log_umi"),
                                      ("unadjusted_sensitivity", "none")]:
                fit = fit_donor_model(part, "hellinger", adjustment)
                if fit is None:
                    raise ValueError(f"Eligible Hellinger model unidentifiable: {assay} {ct} {model}")
                model_rows.append({"assay_release": assay, "endpoint": "Hellinger",
                                   "cell_type": ct, "donor_set": "primary_adult", "model": model, **fit})
    models = add_bh(pd.DataFrame(model_rows), ["assay_release", "endpoint", "donor_set", "model"])
    return donors, models, pd.DataFrame(eligibility_rows)


def empty_hellinger_outputs(assay: str, eligible: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    reasons = pd.DataFrame({
        "assay_release": assay,
        "cell_type": sorted(eligible.loc[eligible["eligible"], "cell_type"]),
    })
    reasons["hellinger_eligible"] = False
    reasons["reason"] = "no_cell_type_passed_frozen_reference_and_donor_floors"
    return pd.DataFrame(), pd.DataFrame(), reasons


def run_assay(assay: str, data_dir: Path) -> None:
    release = RELEASES[assay]
    path = data_dir / "expression" / release["filename"]
    if not path.exists():
        raise FileNotFoundError(path)
    interim = data_dir / "interim"
    interim.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(exist_ok=True)
    adata, obs = read_release_obs(assay, path)
    try:
        donors, verification = verify_obs(assay, adata, obs, data_dir)
        eligible = eligibility_table(obs)
        eligible.insert(0, "assay_release", assay)
        cache_path = interim / f"{assay}_donor_aggregates.npz"
        built = build_cache(assay, adata, obs, donors, eligible, cache_path)
        verification.update(built["verification"])
        cache = built["arrays"]
        with open(RESULTS / f"{PREFIX}-{assay}-verification.json", "w") as handle:
            json.dump(verification, handle, indent=2, sort_keys=True)
            handle.write("\n")
        upsert_csv(RESULTS / f"{PREFIX}-donor-crosswalk.csv", donors,
                   ["assay_release", "cellxgene_donor_id"])
        upsert_csv(RESULTS / f"{PREFIX}-celltype-eligibility.csv", eligible,
                   ["assay_release", "cell_type"])
        comp_counts = composition_counts(obs, assay)
        upsert_csv(RESULTS / f"{PREFIX}-composition-donor-counts.csv", comp_counts,
                   ["assay_release", "donor_id", "cell_type"])
        p1_donors, p1_models = run_p1(assay, cache, eligible)
        upsert_csv(RESULTS / f"{PREFIX}-p1-donor-scores.csv", p1_donors,
                   ["assay_release", "donor_set", "donor_id", "cell_type"])
        upsert_csv(RESULTS / f"{PREFIX}-p1-models.csv", p1_models,
                   ["assay_release", "donor_set", "model", "cell_type"])
        comp_models = run_composition(assay, comp_counts, eligible)
        upsert_csv(RESULTS / f"{PREFIX}-composition-models.csv", comp_models,
                   ["assay_release", "donor_set", "model", "cell_type"])
        h_donors, h_models, h_eligible = run_hellinger(assay, adata, obs, cache, eligible)
        if not h_donors.empty:
            upsert_csv(RESULTS / f"{PREFIX}-hellinger-donors.csv", h_donors,
                       ["assay_release", "cell_type", "donor_id"])
        if not h_models.empty:
            upsert_csv(RESULTS / f"{PREFIX}-hellinger-models.csv", h_models,
                       ["assay_release", "model", "cell_type"])
        upsert_csv(RESULTS / f"{PREFIX}-hellinger-eligibility.csv", h_eligible,
                   ["assay_release", "cell_type"])
    finally:
        adata.file.close()
    print(f"Completed {assay}; results use donor-level inference only.", flush=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def compare_existing() -> pd.DataFrame:
    p1_path = RESULTS / f"{PREFIX}-p1-models.csv"
    if not p1_path.exists():
        raise FileNotFoundError("Run both assays before comparison")
    new = pd.read_csv(p1_path)
    new = new[new["model"].eq("primary_age_sex")].copy()
    old = pd.read_csv(RESULTS / "m9-regressions.csv")
    cohorts = ["Kidney_NormalAdult", "Kidney_LivingDonor", "Kidney_SCR_v15", "Kidney_SN_v15"]
    old = old[old["cohort"].isin(cohorts)].rename(columns={
        "n_donors": "old_n_donors", "slope_per_year": "old_slope_per_year",
        "p": "old_pvalue", "q_bh": "old_q_bh",
    })
    current_donors = pd.read_csv(RESULTS / f"{PREFIX}-donor-crosswalk.csv")
    old_donors = pd.read_csv(RESULTS / "m9-donor-scores.csv", usecols=["cohort", "donor_id"])
    old_donors = old_donors[old_donors["cohort"].isin(cohorts)].drop_duplicates()
    for column in ["n_donors", "slope_per_year", "pvalue", "q_bh"]:
        new[f"{column}_new"] = new[column]
    rows = []

    def shared_donors(assay, donor_set, cohort):
        eligible = current_donors[current_donors["assay_release"] == assay]
        if donor_set == "primary_adult":
            eligible = eligible[eligible["primary_adult"]]
        else:
            eligible = eligible[eligible["sensitivity_all_exact_age"]]
        previous = set(old_donors.loc[old_donors["cohort"] == cohort, "donor_id"].astype(str))
        return len(set(eligible["cellxgene_donor_id"].astype(str)) & previous)

    for assay in ["scCv3", "snCv3"]:
        donor_set = "primary_adult" if assay == "scCv3" else "all_exact_age"
        current = new[(new["assay_release"] == assay) & (new["donor_set"] == donor_set)]
        for cohort, prior in old.groupby("cohort"):
            merged = current.merge(prior, on="cell_type", suffixes=("_new", "_old"))
            relationship = ("overlapping_prior_release_not_independent"
                            if cohort == RELEASES[assay]["prior_cohort"] else "separate_study")
            for _, row in merged.iterrows():
                rows.append({
                    "endpoint": "P1_transport", "assay_release": assay, "existing_cohort": cohort,
                    "cell_type": row["cell_type"], "new_donor_set": donor_set,
                    "new_n_donors": row["n_donors_new"], "old_n_donors": row["old_n_donors"],
                    "n_shared_donor_ids": shared_donors(assay, donor_set, cohort),
                    "new_slope_per_year": row["slope_per_year_new"],
                    "old_slope_per_year": row["old_slope_per_year"],
                    "new_ci95_low": row["ci95_low"], "new_ci95_high": row["ci95_high"],
                    "new_pvalue": row["pvalue"], "old_pvalue": row["old_pvalue"],
                    "new_q_bh": row["q_bh_new"], "old_q_bh": row["old_q_bh"],
                    "sign_concordant": bool(row["slope_per_year_new"] * row["old_slope_per_year"] > 0),
                    "relationship": relationship,
                    "comparison_note": "Same exact cell_type and donor-level age+sex log1p(CPM) P1 estimand; not an independent replication for overlapping releases.",
                })
    h_path = RESULTS / f"{PREFIX}-hellinger-models.csv"
    if h_path.exists():
        hnew = pd.read_csv(h_path)
        hnew = hnew[hnew["model"].eq("primary_depth_adjusted")]
        for column in ["n_donors", "slope_per_year", "pvalue", "q_bh"]:
            hnew[f"{column}_new"] = hnew[column]
        hold = pd.read_csv(RESULTS / "h3-cross-species-vector-tests.csv")
        hold = hold[(hold["species"] == "Human") & hold["tissue"].eq("kidney")
                    & hold["metric"].eq("profile_hellinger") & hold["model"].eq("depth_adjusted")]
        hold = hold.rename(columns={
            "n_donors": "old_n_donors", "slope_per_year": "old_slope_per_year",
            "pvalue": "old_pvalue", "q_bh": "old_q_bh",
        })
        for assay in ["scCv3", "snCv3"]:
            current = hnew[hnew["assay_release"] == assay]
            for cohort, prior in hold.groupby("cohort"):
                merged = current.merge(prior, on="cell_type", suffixes=("_new", "_old"))
                relationship = ("overlapping_prior_release_not_independent"
                                if cohort == RELEASES[assay]["prior_cohort"] else "separate_study")
                for _, row in merged.iterrows():
                    rows.append({
                        "endpoint": "Hellinger", "assay_release": assay, "existing_cohort": cohort,
                        "cell_type": row["cell_type"], "new_donor_set": "primary_adult",
                        "new_n_donors": row["n_donors_new"], "old_n_donors": row["old_n_donors"],
                        "n_shared_donor_ids": shared_donors(assay, "primary_adult", cohort),
                        "new_slope_per_year": row["slope_per_year_new"],
                        "old_slope_per_year": row["old_slope_per_year"],
                        "new_ci95_low": row["ci95_low_new"], "new_ci95_high": row["ci95_high_new"],
                        "new_pvalue": row["pvalue_new"], "old_pvalue": row["old_pvalue"],
                        "new_q_bh": row["q_bh_new"], "old_q_bh": row["old_q_bh"],
                        "sign_concordant": bool(row["slope_per_year_new"] * row["old_slope_per_year"] > 0),
                        "relationship": relationship,
                        "comparison_note": "Same exact cell_type and donor-level cross-fitted young-reference depth-adjusted Hellinger estimand.",
                    })
    result = pd.DataFrame(rows)
    result.to_csv(RESULTS / f"{PREFIX}-existing-comparison.csv", index=False)
    return result


def update_manifest(data_dir: Path) -> None:
    path = REPO / "data" / "manifest" / "lake_kidney.json"
    manifest = json.loads(path.read_text())
    for assay, release in RELEASES.items():
        local = data_dir / "expression" / release["filename"]
        verification_path = RESULTS / f"{PREFIX}-{assay}-verification.json"
        if not local.exists() or not verification_path.exists():
            raise FileNotFoundError(f"Missing {assay} expression or verification")
        verification = json.loads(verification_path.read_text())
        entry = manifest["public_sc_sn_gate"][assay]
        entry["cellxgene_h5ad_url"] = release["url"]
        entry["cellxgene_collection_version_id"] = "379e9a9b-8436-496f-9370-042491bc27d6"
        entry["cellxgene_h5ad_downloaded"] = True
        entry["cellxgene_h5ad_size_bytes"] = local.stat().st_size
        entry["cellxgene_h5ad_sha256"] = sha256(local)
        entry["cellxgene_h5ad_local_path"] = str(
            Path("..") / "agescape-revision-data" / "lake_kidney" / "expression" / release["filename"]
        )
        entry["raw_integer_counts_verified"] = verification["raw_X_all_stored_values_finite_nonnegative_integer"]
        entry["released_annotations_verified"] = True
        entry["donor_age_disease_assay_mapping_verified"] = True
        entry["primary_adult_exact_age_donors"] = verification["primary_adult_donors"]
        entry["analysis_verification"] = f"results/{PREFIX}-{assay}-verification.json"
    manifest["local_metadata"]["directory"] = "../agescape-revision-data/lake_kidney"
    manifest["analysis"] = {
        "completed": "2026-11-14", "frozen_plan": "data/manifest/lake_kidney_analysis_plan.json",
        "script": "scripts/run_followup_lake_kidney.py", "report": "reports/followup-lake-kidney-analysis.md",
        "assays_pooled": False, "inferential_unit": "donor", "primary_adult_age_range": [20, 120],
        "snCv3_14_year_old": "excluded_primary_included_predeclared_all_exact_age_sensitivity",
        "p1_random_null_iterations": NULL_ITERATIONS, "p1_random_seed": NULL_SEED,
        "donor_pseudobulk_gene_models": "not_predeclared_not_performed",
        "controlled_raw_sequencing_pursued": False,
    }
    manifest["verdict"]["expression_analysis_performed"] = True
    manifest["verdict"]["restricted_assay_specific_signature_and_donor_CLR"] = "completed_public_per_assay_H5AD"
    manifest["verdict"]["Hellinger"] = "evaluated_under_frozen_floors_see_eligibility_and_report"
    path.write_text(json.dumps(manifest, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run")
    run.add_argument("--assay", choices=RELEASES, required=True)
    run.add_argument("--data-dir", type=Path, default=DEFAULT_DATA)
    final = sub.add_parser("finalize")
    final.add_argument("--data-dir", type=Path, default=DEFAULT_DATA)
    final.add_argument("--skip-manifest", action="store_true")
    args = parser.parse_args()
    if args.command == "run":
        run_assay(args.assay, args.data_dir)
    else:
        compare_existing()
        if not args.skip_manifest:
            update_manifest(args.data_dir)


if __name__ == "__main__":
    main()
