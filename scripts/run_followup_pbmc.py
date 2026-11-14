#!/usr/bin/env python3
"""Validate and analyze one public PBMC cohort without cross-cohort pooling."""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import re
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from src.aging_programs import load_program_tiers
from src.followup_pbmc import (
    add_bh, clr_transform, common_cell_type, design_matrix, fast_effect, fit_donor_model,
    fit_vectorized, hellinger_distances, known_cell_type, random_signature_scores,
    technical_gene_mask,
)

RESULTS = REPO / "results"
PREFIX = "followup-pbmc"
DEFAULT_DATA = REPO.parent / "agescape-revision-data" / "pbmc"
MIN_CELLS = 20
MIN_DONORS = 6
MIN_REFERENCE_CELLS = 50
MIN_UMI = 1500
NULL_ITERATIONS = 200
NULL_SEED = 20261114
TARGET_GENES = ["PTGER2", "CLIC1", "S100A4", "PRF1", "GZMA", "CCL5"]
TARGET_COMPOSITION = {
    "CD8 naive T": "decrease", "MAIT": "decrease",
    "non-classical monocytes": "increase",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def decode(values) -> np.ndarray:
    return np.array([value.decode() if isinstance(value, bytes) else str(value) for value in values], dtype=str)


def donor_files(data_root: Path, cohort: str, donor: str) -> dict[str, Path]:
    directory = data_root / "raw" / cohort / donor
    return {path.name: path for path in directory.iterdir() if path.is_file()}


def read_features(path: Path) -> pd.DataFrame:
    features = pd.read_csv(path, sep="\t", header=None, dtype=str)
    if features.shape[1] < 2:
        raise ValueError(f"Invalid feature file: {path}")
    features = features.iloc[:, :3].copy()
    features.columns = ["feature_id", "symbol", "feature_type"][:features.shape[1]]
    if not features["feature_id"].is_unique:
        raise ValueError(f"Non-unique feature IDs: {path}")
    return features


def validate_raw_matrix(matrix: sp.csr_matrix) -> dict:
    values = matrix.data
    if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.floor(values)).all():
        raise ValueError("Raw counts failed finite nonnegative integer validation")
    return {
        "shape_genes_by_cells": [int(matrix.shape[0]), int(matrix.shape[1])],
        "nnz": int(matrix.nnz), "count_sum": float(values.sum(dtype=np.float64)),
        "all_stored_counts_finite_nonnegative_integer": True,
    }


def load_zhu(data_root: Path, donor_meta: dict, annotations: pd.DataFrame) -> tuple:
    donor = donor_meta["donor"]
    files = donor_files(data_root, "GSE213516", donor)
    matrix_path = next(path for name, path in files.items() if "matrix.mtx.gz" in name)
    barcode_path = next(path for name, path in files.items() if "barcodes.tsv.gz" in name)
    feature_path = next(path for name, path in files.items() if "features.tsv.gz" in name)
    matrix = sp.csr_matrix(scipy.io.mmread(str(matrix_path)))
    validation = validate_raw_matrix(matrix)
    barcodes = pd.read_csv(barcode_path, header=None, sep="\t")[0].astype(str).to_numpy()
    features = read_features(feature_path)
    if len(barcodes) != matrix.shape[1] or len(features) != matrix.shape[0]:
        raise ValueError(f"Matrix/barcode/feature dimensions differ for {donor}")
    if not pd.Index(barcodes).is_unique:
        raise ValueError(f"Duplicate raw barcodes for {donor}")
    annotation_aliases = {"F40": "HR", "F65": "HYD", "F60": "XMZ", "F46": "ZY"}
    annotation_donor = annotation_aliases.get(donor, donor)
    part = annotations[annotations["sampleID"].eq(annotation_donor)].copy()
    if part.empty or not part["Cohort"].eq("SE").all():
        raise ValueError(f"Missing SE annotations for {donor}")
    if not part["X"].is_unique or not np.all(part["Age"].to_numpy(float) == float(donor_meta["age"])):
        raise ValueError(f"Annotation IDs/ages invalid for {donor}")
    stripped = pd.Index(barcodes).str.replace(r"-1$", "", regex=True)
    baseline_format = part["X"].str.startswith("baseline_").all() and part["X"].str.endswith("-" + donor).all()
    alias_format = donor in annotation_aliases
    if alias_format:
        expected = stripped + "-1-" + annotation_donor
        join_convention = "<barcode>-1-" + annotation_donor
    elif baseline_format:
        expected = "baseline_" + stripped + "-" + donor
        join_convention = "baseline_<barcode>-" + donor
    else:
        expected = donor + "_" + stripped
        join_convention = donor + "_<barcode>"
    lookup = pd.Series(np.arange(len(barcodes)), index=expected)
    if lookup.index.has_duplicates:
        raise ValueError(f"Duplicate normalized raw barcodes for {donor}")
    missing = ~part["X"].isin(lookup.index)
    unmatched_annotations = int(missing.sum())
    part = part.loc[~missing].copy()
    indices = lookup.loc[part["X"]].to_numpy()
    if len(set(indices)) != len(indices):
        raise ValueError(f"Non-one-to-one annotation join for {donor}")
    obs = pd.DataFrame({"cell_type": part["Celltype"].astype(str).to_numpy(), "barcode": part["X"].to_numpy()})
    selected = matrix[:, indices].T.tocsr()
    validation.update({
        "donor": donor, "gsm": donor_meta["gsm"], "matrix_file": matrix_path.name,
        "matrix_sha256": sha256(matrix_path), "raw_barcodes": len(barcodes),
        "annotation_rows": len(part) + unmatched_annotations, "annotation_rows_matched": len(indices),
        "annotation_coverage_of_raw": len(indices) / len(barcodes),
        "annotation_match_fraction": len(indices) / (len(part) + unmatched_annotations),
        "unmatched_annotation_rows_excluded": unmatched_annotations,
        "unmatched_raw_barcodes_excluded": len(barcodes) - len(indices),
        "barcode_join_one_to_one": True, "barcode_join_convention": join_convention,
        "feature_file": feature_path.name, "feature_sha256": sha256(feature_path),
    })
    return selected, obs, features, validation


def load_thomson(data_root: Path, donor_meta: dict, annotations=None) -> tuple:
    donor = donor_meta["donor"]
    files = donor_files(data_root, "GSE214546", donor)
    path = next(path for name, path in files.items() if name.endswith("_labeled.h5"))
    with h5py.File(path, "r") as handle:
        group = handle["matrix"]
        shape = tuple(int(x) for x in group["shape"][:])
        matrix = sp.csc_matrix((group["data"][:], group["indices"][:], group["indptr"][:]), shape=shape)
        barcodes = decode(group["barcodes"][:])
        observations = group["observations"]
        labels = decode(observations["seurat_pbmc_type"][:])
        sample_values = decode(observations["pbmc_sample_id"][:]) if "pbmc_sample_id" in observations else np.array([])
        pool_values = decode(observations["pool_id"][:]) if "pool_id" in observations else np.array([])
        features = pd.DataFrame({
            "feature_id": decode(group["features"]["id"][:]),
            "symbol": decode(group["features"]["name"][:]),
            "feature_type": decode(group["features"]["feature_type"][:]),
        })
    validation = validate_raw_matrix(matrix)
    if len(barcodes) != matrix.shape[1] or len(features) != matrix.shape[0] or len(labels) != matrix.shape[1]:
        raise ValueError(f"H5 matrix/metadata dimensions differ for {donor}")
    if not pd.Index(barcodes).is_unique or not features["feature_id"].is_unique:
        raise ValueError(f"Duplicate H5 barcodes/features for {donor}")
    obs = pd.DataFrame({"cell_type": labels, "barcode": barcodes})
    match = re.search(r"_(B\d+)-(P\d+)_", path.name)
    if not match:
        raise ValueError(f"Cannot parse batch/pool from {path.name}")
    validation.update({
        "donor": donor, "gsm": donor_meta["gsm"], "matrix_file": path.name,
        "matrix_sha256": sha256(path), "raw_barcodes": len(barcodes),
        "annotation_rows": len(labels), "annotation_rows_matched": len(labels),
        "annotation_coverage_of_raw": 1.0, "unmatched_raw_barcodes_excluded": 0,
        "barcode_join_one_to_one": True, "batch_from_filename": match.group(1),
        "pool_from_filename": match.group(2), "h5_pool_values": sorted(set(pool_values)),
        "h5_pbmc_sample_values": sorted(set(sample_values)),
        "feature_sha256": sha256(path),
    })
    return matrix.T.tocsr(), obs, features, validation


def aggregate_donor(matrix: sp.csr_matrix, obs: pd.DataFrame, donor_meta: dict,
                    batch: str = "", pool: str = "") -> dict:
    keep = known_cell_type(obs["cell_type"])
    if not keep.all():
        obs = obs.loc[keep].copy()
        matrix = matrix[keep].copy()
    libraries = np.asarray(matrix.sum(axis=1)).ravel()
    if (libraries <= 0).any():
        raise ValueError(f"Zero-library annotated cell for {donor_meta['donor']}")
    labels = sorted(obs["cell_type"].unique())
    result = {
        "donor": [], "cell_type": [], "age": [], "sex": [], "cmv": [], "batch": [], "pool": [],
        "n_cells": [], "mean_log_umi": [], "sum_umi": [], "log_sums": [], "counts": [],
        "qc_n": [], "qc_mean_log_umi": [], "qc_sum_umi": [], "qc_counts": [],
    }
    for label in labels:
        mask = (obs["cell_type"].to_numpy() == label)
        part = matrix[mask].astype(np.float64)
        lib = libraries[mask]
        normalized = part.copy()
        normalized.data *= np.repeat(1e6 / lib, np.diff(normalized.indptr))
        np.log1p(normalized.data, out=normalized.data)
        qc = lib >= MIN_UMI
        result["donor"].append(donor_meta["donor"])
        result["cell_type"].append(label)
        result["age"].append(float(donor_meta["age"]))
        result["sex"].append(donor_meta["sex"])
        result["cmv"].append(donor_meta.get("cmv", ""))
        result["batch"].append(batch)
        result["pool"].append(pool)
        result["n_cells"].append(len(lib))
        result["mean_log_umi"].append(float(np.log(lib).mean()))
        result["sum_umi"].append(float(lib.sum()))
        result["log_sums"].append(np.asarray(normalized.sum(axis=0)).ravel().astype(np.float32))
        result["counts"].append(np.asarray(part.sum(axis=0)).ravel().astype(np.uint32))
        result["qc_n"].append(int(qc.sum()))
        result["qc_mean_log_umi"].append(float(np.log(lib[qc]).mean()) if qc.any() else np.nan)
        result["qc_sum_umi"].append(float(lib[qc].sum()))
        result["qc_counts"].append(np.asarray(part[qc].sum(axis=0)).ravel().astype(np.uint32))
    for key in ["log_sums", "counts", "qc_counts"]:
        result[key] = np.vstack(result[key]).astype(np.float32 if key == "log_sums" else np.uint32)
    for key in ["n_cells", "qc_n"]:
        result[key] = np.asarray(result[key], dtype=np.int64)
    for key in ["age", "mean_log_umi", "sum_umi", "qc_mean_log_umi", "qc_sum_umi"]:
        result[key] = np.asarray(result[key], dtype=float)
    for key in ["donor", "cell_type", "sex", "cmv", "batch", "pool"]:
        result[key] = np.asarray(result[key], dtype=str)
    return result


def build_cache(cohort: str, data_root: Path, manifest: dict) -> tuple[dict, pd.DataFrame, dict]:
    interim = data_root / "interim" / cohort
    interim.mkdir(parents=True, exist_ok=True)
    cache_path = interim / "donor-aggregates.npz"
    verification_path = interim / "validation.json"
    donor_records = manifest[cohort]["donors"]
    annotations = None
    if cohort == "GSE213516":
        annotation_path = data_root / "raw" / cohort / manifest[cohort]["annotations"]["file"]
        annotations = pd.read_csv(annotation_path)
        annotations = annotations[annotations["Cohort"].eq("SE")].copy()
    if cache_path.exists() and verification_path.exists():
        with np.load(cache_path, allow_pickle=False) as data:
            cache = {key: data[key] for key in data.files}
        verification = json.loads(verification_path.read_text())
        features = pd.DataFrame({"feature_id": cache["feature_id"].astype(str), "symbol": cache["symbol"].astype(str)})
        return cache, features, verification
    loader = load_zhu if cohort == "GSE213516" else load_thomson
    units, verifications, features_reference = [], [], None
    for donor_meta in donor_records:
        print(f"[{cohort}] validating {donor_meta['donor']}", flush=True)
        matrix, obs, features, verification = loader(data_root, donor_meta, annotations)
        if features_reference is None:
            features_reference = features[["feature_id", "symbol"]].copy()
        elif not features_reference.equals(features[["feature_id", "symbol"]]):
            raise ValueError(f"Feature reference differs for {donor_meta['donor']}")
        batch, pool = (verification.get("batch_from_filename", ""), verification.get("pool_from_filename", ""))
        unit = aggregate_donor(matrix, obs, donor_meta, batch, pool)
        verification["released_labels"] = sorted(obs["cell_type"].unique())
        verifications.append(verification)
        units.append(unit)
        del matrix, obs, features, unit
        gc.collect()
    cache = {key: np.concatenate([unit[key] for unit in units], axis=0) for key in units[0]}
    cache["feature_id"] = features_reference["feature_id"].to_numpy(str)
    cache["symbol"] = features_reference["symbol"].to_numpy(str)
    np.savez(cache_path, **cache)
    label_inventory = []
    for label in sorted(set(cache["cell_type"])):
        mask = cache["cell_type"] == label
        label_inventory.append({
            "cell_type": label, "common_cell_type": common_cell_type(label),
            "donors_observed": int(mask.sum()), "cells": int(cache["n_cells"][mask].sum()),
        })
    pd.DataFrame(label_inventory).to_csv(RESULTS / f"{PREFIX}-{cohort}-label-inventory.csv", index=False)
    verification = {
        "cohort": cohort, "donors": verifications,
        "all_donor_features_identical": True,
        "unique_feature_ids": len(set(cache["feature_id"])),
        "unique_gene_symbols": len(set(cache["symbol"])),
        "duplicate_symbol_feature_rows": int(pd.Index(cache["symbol"]).duplicated(keep=False).sum()),
        "cache_path": str(cache_path), "cache_bytes": cache_path.stat().st_size,
    }
    if cohort == "GSE213516":
        verification["genome_reference_adjudication"] = {
            "geo_record": "hg19", "primary_methods": "GRCh38",
            "actual_feature_columns": ["Ensembl gene ID", "gene symbol", "feature type"],
            "matrix_metadata": "Cell Ranger 3.0.0; no reference genome identity in released feature/matrix records",
            "verdict": "unresolved_from_released_features",
            "restriction": "All Zhu gene-level results are symbol-level only; no coordinate, locus, or assembly-specific claim. Ambiguous duplicate symbols are excluded from P1/targeted tests.",
        }
        verification["annotation_sha256"] = sha256(annotation_path)
    verification_path.write_text(json.dumps(verification, indent=2) + "\n")
    features = features_reference
    return cache, features, verification


def donor_table(cache: dict) -> pd.DataFrame:
    columns = ["donor", "age", "sex", "cmv", "batch", "pool"]
    return pd.DataFrame({key: cache[key] for key in columns}).drop_duplicates("donor").reset_index(drop=True)


def select_donors(cache: dict, design: str, metadata: dict) -> pd.DataFrame:
    donors = donor_table(cache).copy()
    records = {record["donor"]: record for record in metadata["donors"]}
    donors["strict_group"] = donors["donor"].map(lambda x: records[x].get("strict_group"))
    donors["atlas_included"] = donors["donor"].map(lambda x: bool(records[x].get("atlas_included")))
    if design == "strict_group":
        donors = donors[donors["strict_group"].notna()].copy()
        donors["old"] = donors["strict_group"].eq("old").astype(float)
    elif design == "atlas_group":
        donors = donors[donors["atlas_included"]].copy()
        donors["old"] = donors["age"].ge(60).astype(float)
    elif design == "exact_age":
        donors["old"] = np.nan
    else:
        raise ValueError(design)
    return donors


def unit_frame(cache: dict, cell_type: str, design: str, metadata: dict,
               cell_floor: str = "n_cells") -> pd.DataFrame:
    mask = (cache["cell_type"] == cell_type) & (cache[cell_floor] >= MIN_CELLS)
    columns = ["donor", "cell_type", "age", "sex", "cmv", "batch", "pool", "n_cells",
               "mean_log_umi", "sum_umi", "qc_n", "qc_mean_log_umi", "qc_sum_umi"]
    frame = pd.DataFrame({key: cache[key][mask] for key in columns})
    donors = select_donors(cache, design, metadata)
    frame = frame.merge(donors[["donor", "old"]], on="donor", how="inner", validate="many_to_one")
    return frame


def base_row(cohort: str, endpoint: str, cell_type: str, design: str, model: str) -> dict:
    return {
        "cohort": cohort, "endpoint": endpoint, "cell_type": cell_type,
        "common_cell_type": common_cell_type(cell_type), "design": design, "model": model,
    }


def eligibility(cache: dict, metadata: dict) -> pd.DataFrame:
    rows = []
    for cell_type in sorted(set(cache["cell_type"])):
        for design in ["strict_group", "atlas_group", "exact_age"]:
            frame = unit_frame(cache, cell_type, design, metadata)
            term = "age" if design == "exact_age" else "old"
            group_ok = True if design == "exact_age" else min((frame[term] == 0).sum(), (frame[term] == 1).sum()) >= 3
            eligible = len(frame) >= MIN_DONORS and group_ok
            rows.append({
                **base_row(metadata["accession"] if "accession" in metadata else metadata.get("cohort", ""), "eligibility", cell_type, design, "primary"),
                "n_donors": len(frame), "n_cells": int(frame["n_cells"].sum()),
                "n_young": int((frame[term] == 0).sum()) if term == "old" else np.nan,
                "n_old": int((frame[term] == 1).sum()) if term == "old" else np.nan,
                "eligible": eligible, "reason": "pass" if eligible else "donor_or_group_floor",
            })
    return pd.DataFrame(rows)


def run_p1(cohort: str, cache: dict, metadata: dict, eligible: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    program = load_program_tiers(REPO, "human")["P1_Clean"]
    symbols = cache["symbol"].astype(str)
    unique_symbols = ~np.asarray(pd.Index(symbols).duplicated(keep=False))
    up = np.flatnonzero(np.isin(symbols, program["up"]) & unique_symbols)
    down = np.flatnonzero(np.isin(symbols, program["down"]) & unique_symbols)
    if not len(up) or not len(down):
        raise ValueError(f"P1 lacks unique detected genes for {cohort}")
    model_rows, donor_rows = [], []
    rng = np.random.default_rng(NULL_SEED)
    for cell_type in sorted(eligible.loc[eligible["eligible"], "cell_type"].unique()):
        for design in ["strict_group", "atlas_group", "exact_age"]:
            if not eligible.query("cell_type == @cell_type and design == @design")["eligible"].iloc[0]:
                continue
            mask = (cache["cell_type"] == cell_type) & (cache["n_cells"] >= MIN_CELLS)
            indices = np.flatnonzero(mask)
            frame = unit_frame(cache, cell_type, design, metadata)
            lookup = pd.Series(indices, index=cache["donor"][mask])
            order = lookup.loc[frame["donor"]].to_numpy()
            matrix = (cache["log_sums"][order] / cache["n_cells"][order, None]).astype(float).T
            frame["score"] = matrix[up].mean(axis=0) - matrix[down].mean(axis=0)
            frame["cohort"] = cohort
            frame["design"] = design
            frame["common_cell_type"] = common_cell_type(cell_type)
            donor_rows.append(frame)
            variants = ["primary"] + (["depth"] if design == "strict_group" else [])
            if cohort == "GSE214546" and design == "strict_group":
                variants += ["cmv", "batch", "pool"]
            for variant in variants:
                fit = fit_donor_model(frame, "score", design, variant, "P1_transport")
                model_rows.append({
                    **base_row(cohort, "P1_transport", cell_type, design, variant),
                    "p1_up_detected": len(up), "p1_down_detected": len(down),
                    "p1_coverage": (len(up) + len(down)) / (len(program["up"]) + len(program["down"])), **fit,
                })
            if design == "strict_group":
                observed = next(row["effect"] for row in reversed(model_rows)
                                if row["cell_type"] == cell_type and row["design"] == "strict_group" and row["model"] == "primary")
                all_p1 = np.flatnonzero(np.isin(symbols, program["up"] + program["down"]))
                random_scores = random_signature_scores(matrix, up, down, NULL_ITERATIONS, rng, exclude=all_p1)
                random_effects = np.array([fast_effect(score, frame, design) for score in random_scores])
                null_mean = float(random_effects.mean())
                empirical = (1 + int(np.sum(np.abs(random_effects - null_mean) >= abs(observed - null_mean)))) / (NULL_ITERATIONS + 1)
                target = next(row for row in reversed(model_rows)
                              if row["cell_type"] == cell_type and row["design"] == "strict_group" and row["model"] == "primary")
                target.update({
                    "random_iterations": NULL_ITERATIONS, "random_seed": NULL_SEED,
                    "random_mean_effect": null_mean, "random_sd": float(random_effects.std(ddof=1)),
                    "random_q025": float(np.quantile(random_effects, 0.025)),
                    "random_q975": float(np.quantile(random_effects, 0.975)),
                    "random_pvalue": float(empirical),
                })
                print(f"[{cohort}] P1 {cell_type}: n={len(frame)}, null p={empirical:.4f}", flush=True)
    models = pd.DataFrame(model_rows)
    models = add_bh(models, ["cohort", "endpoint", "design", "model"])
    models["random_q_bh"] = np.nan
    primary = models["design"].eq("strict_group") & models["model"].eq("primary") & models["status"].eq("fit")
    if primary.any():
        models.loc[primary, "random_q_bh"] = multipletests(models.loc[primary, "random_pvalue"], method="fdr_bh")[1]
    donors = pd.concat(donor_rows, ignore_index=True) if donor_rows else pd.DataFrame()
    return donors, models


def pseudobulk_outcomes(counts: np.ndarray) -> np.ndarray:
    totals = counts.sum(axis=1)
    return np.log2((counts.T + 0.5) / (totals + 1) * 1e6)


def run_pseudobulk(cohort: str, cache: dict, metadata: dict,
                   eligible: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    targeted_rows, exploratory_frames = [], []
    symbols = cache["symbol"].astype(str)
    unique_symbols = ~np.asarray(pd.Index(symbols).duplicated(keep=False))
    for cell_type in sorted(eligible.loc[eligible["eligible"], "cell_type"].unique()):
        for design in ["strict_group", "atlas_group", "exact_age"]:
            if not eligible.query("cell_type == @cell_type and design == @design")["eligible"].iloc[0]:
                continue
            mask = (cache["cell_type"] == cell_type) & (cache["n_cells"] >= MIN_CELLS)
            indices = np.flatnonzero(mask)
            frame = unit_frame(cache, cell_type, design, metadata)
            lookup = pd.Series(indices, index=cache["donor"][mask])
            order = lookup.loc[frame["donor"]].to_numpy()
            counts = cache["counts"][order].astype(np.float64)
            outcomes = pseudobulk_outcomes(counts)
            detected = (counts > 0).mean(axis=0) >= 0.20
            variants = ["primary"] + (["depth"] if design == "strict_group" else [])
            for variant in variants:
                fitted = fit_vectorized(outcomes[detected], frame, design, variant)
                if fitted is None:
                    continue
                selected = np.flatnonzero(detected)
                result = pd.DataFrame({
                    "cohort": cohort, "endpoint": "pseudobulk_exploratory", "cell_type": cell_type,
                    "common_cell_type": common_cell_type(cell_type), "design": design, "model": variant,
                    "feature_id": cache["feature_id"][selected].astype(str), "gene": symbols[selected],
                    "effect_term": "age_per_year" if design == "exact_age" else "old_vs_young",
                    "n_donors": len(frame), "n_cells": int(frame["n_cells"].sum()),
                    "n_young": int((frame["old"] == 0).sum()) if design != "exact_age" else np.nan,
                    "n_old": int((frame["old"] == 1).sum()) if design != "exact_age" else np.nan,
                    "status": "fit", **{key: fitted[key] for key in ["effect", "se", "ci95_low", "ci95_high", "pvalue"]},
                })
                result["q_bh"] = multipletests(result["pvalue"], method="fdr_bh")[1]
                exploratory_frames.append(result)
            targeted_variants = ["primary"] + (["depth"] if design == "strict_group" else [])
            if cohort == "GSE214546" and design == "strict_group":
                targeted_variants += ["cmv", "batch", "pool"]
            if common_cell_type(cell_type) == "CD8 naive T":
                for variant in targeted_variants:
                    for gene in TARGET_GENES:
                        selected = np.flatnonzero((symbols == gene) & unique_symbols)
                        base = {
                            **base_row(cohort, "pseudobulk_targeted", cell_type, design, variant),
                            "gene": gene, "expected_direction": "increase",
                        }
                        if not len(selected):
                            targeted_rows.append({**base, "status": "untestable", "reason": "gene_missing_or_ambiguous_symbol"})
                            continue
                        fitted = fit_vectorized(outcomes[selected], frame, design, variant)
                        if fitted is None:
                            targeted_rows.append({**base, "feature_id": cache["feature_id"][selected[0]],
                                                  "status": "untestable", "reason": "rank_or_group_floor"})
                            continue
                        targeted_rows.append({
                            **base, "feature_id": cache["feature_id"][selected[0]],
                            "status": "fit", "n_donors": len(frame), "n_cells": int(frame["n_cells"].sum()),
                            **{key: float(fitted[key][0]) for key in ["effect", "se", "ci95_low", "ci95_high", "pvalue"]},
                        })
    targeted = add_bh(pd.DataFrame(targeted_rows), ["cohort", "endpoint", "design", "model"])
    exploratory = pd.concat(exploratory_frames, ignore_index=True) if exploratory_frames else pd.DataFrame()
    return targeted, exploratory


def composition_counts(cache: dict, metadata: dict, design: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    donors = select_donors(cache, design, metadata)
    all_labels = sorted(set(cache["cell_type"]))
    unit_labels = pd.DataFrame({"donor": cache["donor"], "cell_type": cache["cell_type"],
                                "n_cells": cache["n_cells"]})
    wide = unit_labels.pivot(index="donor", columns="cell_type", values="n_cells").fillna(0)
    wide = wide.reindex(index=donors["donor"], columns=all_labels, fill_value=0).astype(int)
    meta = donors.set_index("donor", drop=False)
    meta["total_cells"] = wide.sum(axis=1)
    return wide, meta


def run_composition(cohort: str, cache: dict, metadata: dict,
                    eligible: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    exploratory_rows, targeted_rows, donor_frames = [], [], []
    for design in ["strict_group", "atlas_group", "exact_age"]:
        wide, meta = composition_counts(cache, metadata, design)
        clr = clr_transform(wide)
        donor_long = wide.stack().rename("n_cells").reset_index()
        donor_long["cohort"] = cohort
        donor_long["design"] = design
        donor_long["total_cells"] = donor_long["donor"].map(meta["total_cells"])
        donor_frames.append(donor_long)
        for cell_type in sorted(eligible.loc[eligible["eligible"], "cell_type"].unique()):
            if not eligible.query("cell_type == @cell_type and design == @design")["eligible"].iloc[0]:
                continue
            frame = meta.copy().reset_index(drop=True)
            frame["n_cells"] = wide[cell_type].to_numpy()
            frame["clr"] = clr[cell_type].to_numpy()
            frame["fraction"] = frame["n_cells"] / frame["total_cells"]
            variants = ["primary"] + (["log_total_cells"] if design == "strict_group" else [])
            if cohort == "GSE214546" and design == "strict_group":
                variants += ["cmv", "batch", "pool"]
            for variant in variants:
                fit = fit_donor_model(frame, "clr", design, variant, "CLR_composition")
                exploratory_rows.append({**base_row(cohort, "CLR_composition", cell_type, design, variant), **fit})
            if design == "strict_group":
                fit = fit_donor_model(frame, "fraction", design, "primary", "CLR_composition")
                exploratory_rows.append({**base_row(cohort, "raw_fraction_sensitivity", cell_type, design, "primary"), **fit})
        for common, direction in TARGET_COMPOSITION.items():
            matching = [label for label in wide.columns if common_cell_type(label) == common]
            base = base_row(cohort, "CLR_targeted", common, design, "primary")
            base["expected_direction"] = direction
            if not matching:
                targeted_rows.append({**base, "cell_type": common, "common_cell_type": common,
                                      "status": "untestable", "reason": "released_label_absent_or_not_harmonizable"})
                continue
            label = matching[0]
            pass_row = eligible.query("cell_type == @label and design == @design")["eligible"].iloc[0]
            if not pass_row:
                targeted_rows.append({**base, "cell_type": label, "status": "untestable",
                                      "reason": "donor_or_cell_floor"})
                continue
            frame = meta.copy().reset_index(drop=True)
            frame["n_cells"] = wide[label].to_numpy()
            frame["clr"] = clr[label].to_numpy()
            fit = fit_donor_model(frame, "clr", design, "primary", "CLR_composition")
            targeted_rows.append({**base, "cell_type": label, **fit})
    exploratory = add_bh(pd.DataFrame(exploratory_rows), ["cohort", "endpoint", "design", "model"])
    targeted = add_bh(pd.DataFrame(targeted_rows), ["cohort", "endpoint", "design", "model"])
    return exploratory, targeted, pd.concat(donor_frames, ignore_index=True)


def run_hellinger(cohort: str, cache: dict, metadata: dict, data_root: Path,
                  annotations: pd.DataFrame | None) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    technical = technical_gene_mask(cache["symbol"].astype(str))
    keep_genes = ~technical
    eligibility_rows, donor_rows, model_rows = [], [], []
    references_by_type = {}
    for cell_type in sorted(set(cache["cell_type"])):
        strict_young = np.flatnonzero((cache["cell_type"] == cell_type) & (cache["age"] <= 40)
                                      & (cache["qc_n"] >= MIN_CELLS))
        total_reference = cache["qc_counts"][strict_young].sum(axis=0).astype(np.int64)
        total_reference[technical] = 0
        total_young_cells = int(cache["qc_n"][strict_young].sum())
        references_by_type[cell_type] = {}
        for design in ["strict_group", "atlas_group", "exact_age"]:
            frame = unit_frame(cache, cell_type, design, metadata, cell_floor="qc_n")
            scored, minimum_reference = [], []
            for row in frame.itertuples():
                index = np.flatnonzero((cache["cell_type"] == cell_type) & (cache["donor"] == row.donor))[0]
                own = cache["qc_counts"][index].astype(np.int64) if cache["age"][index] <= 40 else np.zeros(len(technical), dtype=np.int64)
                own[technical] = 0
                reference = total_reference - own
                own_cells = int(cache["qc_n"][index]) if cache["age"][index] <= 40 else 0
                reference_n = total_young_cells - own_cells
                if len(strict_young) < 2 or reference_n < MIN_REFERENCE_CELLS or reference.sum() <= 0:
                    continue
                references_by_type[cell_type][(design, row.donor)] = reference / reference.sum()
                scored.append(row.donor)
                minimum_reference.append(reference_n)
            eligible = len(scored) >= MIN_DONORS
            eligibility_rows.append({
                **base_row(cohort, "Hellinger_eligibility", cell_type, design, "primary"),
                "n_donors_qc20": len(frame), "n_young_reference_donors": len(strict_young),
                "n_scored_donors": len(scored), "min_cross_fitted_reference_cells": min(minimum_reference) if minimum_reference else 0,
                "eligible": eligible,
                "reason": "pass" if eligible else "young_reference_or_scored_donor_floor",
            })
    needed = {(cell_type, donor) for cell_type, mapping in references_by_type.items() for _, donor in mapping}
    loader = load_zhu if cohort == "GSE213516" else load_thomson
    donor_lookup = {record["donor"]: record for record in metadata["donors"]}
    for donor in sorted({donor for _, donor in needed}):
        donor_meta = donor_lookup[donor]
        matrix, obs, features, _ = loader(data_root, donor_meta, annotations)
        libraries = np.asarray(matrix.sum(axis=1)).ravel()
        nontechnical = np.asarray(matrix[:, keep_genes].sum(axis=1)).ravel()
        for cell_type in sorted({ct for ct, d in needed if d == donor}):
            mask = (obs["cell_type"].to_numpy() == cell_type) & (libraries >= MIN_UMI) & (nontechnical > 0)
            if not mask.any():
                continue
            for design in ["strict_group", "atlas_group", "exact_age"]:
                key = (design, donor)
                if key not in references_by_type[cell_type]:
                    continue
                distances = hellinger_distances(matrix[mask], references_by_type[cell_type][key], nontechnical[mask])
                donor_rows.append({
                    "cohort": cohort, "cell_type": cell_type, "common_cell_type": common_cell_type(cell_type),
                    "design": design, "donor": donor, "age": float(donor_meta["age"]),
                    "sex": donor_meta["sex"], "cmv": donor_meta.get("cmv", ""),
                    "n_cells": int(mask.sum()), "mean_log_umi": float(np.log(libraries[mask]).mean()),
                    "qc_mean_log_umi": float(np.log(libraries[mask]).mean()),
                    "mean_umi": float(libraries[mask].mean()), "hellinger": float(distances.mean()),
                    "cross_fitted_reference_cells": int(references_by_type[cell_type].get((design, "reference_n_" + donor), 0)),
                })
                cache_index = np.flatnonzero((cache["cell_type"] == cell_type) & (cache["donor"] == donor))[0]
                own_cells = int(cache["qc_n"][cache_index]) if cache["age"][cache_index] <= 40 else 0
                total_young = np.flatnonzero((cache["cell_type"] == cell_type) & (cache["age"] <= 40) & (cache["qc_n"] >= MIN_CELLS))
                donor_rows[-1]["cross_fitted_reference_cells"] = int(cache["qc_n"][total_young].sum()) - own_cells
        del matrix, obs, features
        gc.collect()
    donors = pd.DataFrame(donor_rows)
    eligibility_frame = pd.DataFrame(eligibility_rows)
    if not donors.empty:
        cache_donors = donor_table(cache)
        for (cell_type, design), part in donors.groupby(["cell_type", "design"]):
            eligibility_row = eligibility_frame.query("cell_type == @cell_type and design == @design")
            if eligibility_row.empty or not eligibility_row["eligible"].iloc[0]:
                continue
            part = part.merge(cache_donors[["donor", "batch", "pool"]], on="donor", how="left", validate="many_to_one")
            selected = select_donors(cache, design, metadata)[["donor", "old"]]
            part = part.merge(selected, on="donor", how="left", validate="many_to_one")
            variants = ["primary"] + (["unadjusted"] if design == "strict_group" else [])
            if cohort == "GSE214546" and design == "strict_group":
                variants += ["cmv_depth", "batch_depth", "pool_depth"]
            for variant in variants:
                fit = fit_donor_model(part, "hellinger", design, variant, "Hellinger")
                model_rows.append({**base_row(cohort, "Hellinger", cell_type, design, variant), **fit})
    models = add_bh(pd.DataFrame(model_rows), ["cohort", "endpoint", "design", "model"])
    return donors, models, eligibility_frame


def identifiability(cohort: str, cache: dict, metadata: dict) -> pd.DataFrame:
    rows = []
    for design in ["strict_group", "atlas_group", "exact_age"]:
        donors = select_donors(cache, design, metadata)
        for variant in ["primary", "depth", "log_total_cells", "cmv", "batch", "pool"]:
            if cohort != "GSE214546" and variant in {"cmv", "batch", "pool"}:
                continue
            frame = donors.copy()
            frame["n_cells"] = 1
            frame["mean_log_umi"] = frame["donor"].map(pd.DataFrame({"donor": cache["donor"], "value": cache["mean_log_umi"]}).groupby("donor")["value"].mean())
            frame["qc_mean_log_umi"] = frame["mean_log_umi"]
            frame["total_cells"] = frame["donor"].map(pd.DataFrame({"donor": cache["donor"], "value": cache["n_cells"]}).groupby("donor")["value"].sum())
            X, names = design_matrix(frame, design, variant, "P1_transport")
            rows.append({
                "cohort": cohort, "design": design, "model": variant, "n_donors": len(frame),
                "n_parameters": X.shape[1], "rank": int(np.linalg.matrix_rank(X)),
                "residual_df": len(frame) - X.shape[1],
                "identifiable": np.linalg.matrix_rank(X) == X.shape[1] and len(frame) > X.shape[1] + 1,
                "covariates": ";".join(names),
                "sex_levels": ";".join(sorted(frame["sex"].unique())),
                "cmv_levels": ";".join(sorted(frame["cmv"].unique())),
                "batch_levels": ";".join(sorted(frame["batch"].unique())),
                "pool_levels": ";".join(sorted(frame["pool"].unique())),
            })
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", required=True, choices=["GSE213516", "GSE214546"])
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA)
    args = parser.parse_args()
    started = time.time()
    RESULTS.mkdir(exist_ok=True)
    full_manifest = json.loads((REPO / "data" / "manifest" / "pbmc_geo_cohorts.json").read_text())
    metadata = dict(full_manifest[args.cohort])
    metadata["accession"] = args.cohort
    cache, features, verification = build_cache(args.cohort, args.data_root, full_manifest)
    cache["batch"] = cache["batch"].astype(str)
    cache["pool"] = cache["pool"].astype(str)
    eligible = eligibility(cache, metadata)
    eligible.to_csv(RESULTS / f"{PREFIX}-{args.cohort}-eligibility.csv", index=False)
    identifiability(args.cohort, cache, metadata).to_csv(
        RESULTS / f"{PREFIX}-{args.cohort}-covariate-identifiability.csv", index=False)
    annotations = None
    if args.cohort == "GSE213516":
        annotations = pd.read_csv(args.data_root / "raw" / args.cohort / metadata["annotations"]["file"])
        annotations = annotations[annotations["Cohort"].eq("SE")].copy()
    p1_donors, p1_models = run_p1(args.cohort, cache, metadata, eligible)
    p1_donors.to_csv(RESULTS / f"{PREFIX}-{args.cohort}-p1-donor-scores.csv", index=False)
    p1_models.to_csv(RESULTS / f"{PREFIX}-{args.cohort}-p1-models.csv", index=False)
    targeted, exploratory = run_pseudobulk(args.cohort, cache, metadata, eligible)
    targeted.to_csv(RESULTS / f"{PREFIX}-{args.cohort}-pseudobulk-targeted.csv", index=False)
    exploratory.to_csv(RESULTS / f"{PREFIX}-{args.cohort}-pseudobulk-exploratory.csv.gz", index=False, compression="gzip")
    comp, comp_target, comp_donors = run_composition(args.cohort, cache, metadata, eligible)
    comp.to_csv(RESULTS / f"{PREFIX}-{args.cohort}-clr-models.csv", index=False)
    comp_target.to_csv(RESULTS / f"{PREFIX}-{args.cohort}-clr-targeted.csv", index=False)
    comp_donors.to_csv(RESULTS / f"{PREFIX}-{args.cohort}-composition-donor-counts.csv", index=False)
    h_donors, h_models, h_eligible = run_hellinger(args.cohort, cache, metadata, args.data_root, annotations)
    h_donors.to_csv(RESULTS / f"{PREFIX}-{args.cohort}-hellinger-donors.csv", index=False)
    h_models.to_csv(RESULTS / f"{PREFIX}-{args.cohort}-hellinger-models.csv", index=False)
    h_eligible.to_csv(RESULTS / f"{PREFIX}-{args.cohort}-hellinger-eligibility.csv", index=False)
    import resource
    verification["resource_use"] = {
        "elapsed_seconds": time.time() - started,
        "max_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "guard_memory_gb": 8, "guard_min_free_gb": 30, "guard_max_file_gb": 3,
        "parallel_heavy_jobs": 1,
    }
    (RESULTS / f"{PREFIX}-{args.cohort}-validation.json").write_text(json.dumps(verification, indent=2) + "\n")
    print(f"Completed {args.cohort} in {time.time() - started:.1f} seconds", flush=True)


if __name__ == "__main__":
    main()
