#!/usr/bin/env python3
"""Block 3 full profile-cloud donor dispersion in Gao heart and Yazar blood."""

from __future__ import annotations

import argparse
import gc
import hashlib
import h5py
import json
import math
import os
from pathlib import Path
import sys
import time

for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "4"

import numpy as np
import pandas as pd
from scipy import sparse as sp
from statsmodels.stats.multitest import multipletests

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.run_followup_between_subject_gao_yazar import (
    GAO_H5AD,
    YAZAR_CACHE,
    YAZAR_H5AD,
    build_human_gene_projection,
    read_csr_chunk,
    reconcile_gao_metadata,
    reconcile_yazar_metadata,
)
from scripts.run_h3_expression_vector import hellinger_distance
from src.followup_muscle import count_gate
from src.followup_organs import write_clean_csv, write_clean_markdown

RESULT_CSV = REPO_ROOT / "results" / "followup-11pt-block3-profile-cloud.csv"
REPORT_MD = REPO_ROOT / "reports" / "followup-11pt-block3-profile-cloud.md"
DATA_ROOT = REPO_ROOT.parent / f"{REPO_ROOT.name}-data" / "tmp_organs"
CACHE_ROOT = DATA_ROOT / "block3_profile_cloud_cache"
ANALYSIS_VERSION = "2026-11-20-v1"
MIN_CELLS = 20
MIN_DONORS_PER_BIN = 6
MIN_DEPTH = 1500
N_PERM = 2000
PERM_SEED = 20261118
SPLIT_SEEDS = tuple(20261116 + i for i in range(5))
CHUNK_ROWS = 8000
MDE80_Z = 2.801585218497961
ARMS = ("batch_residualized", "raw")
COMPARISONS = (("OLD_vs_YOUNG", "OLD", "YOUNG"), ("OLD_vs_MID", "OLD", "MID"), ("MID_vs_YOUNG", "MID", "YOUNG"))
BINS = ("YOUNG", "MID", "OLD")


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def stable_seed(base: int, *parts: object) -> int:
    payload = (str(base) + "|" + "|".join(str(p) for p in parts)).encode("utf-8")
    return int.from_bytes(hashlib.blake2b(payload, digest_size=8).digest(), "little") & 0x7FFFFFFFFFFFFFFF


def cache_signature(meta: dict[str, object], audit: dict[str, object]) -> dict[str, object]:
    path = Path(meta["h5_path"])
    stat = path.stat()
    return {
        "analysis_version": ANALYSIS_VERSION,
        "dataset": meta["dataset"],
        "input_path": str(path),
        "input_size": stat.st_size,
        "input_mtime_ns": stat.st_mtime_ns,
        "matrix_slot": meta["matrix_slot"],
        "n_obs": meta["n_obs"],
        "n_raw_genes": len(meta["gene_symbols"]),
        "n_clean_genes": audit["n_genes_retained"],
        "technical_excluded": audit["n_technical_excluded"],
        "min_depth": MIN_DEPTH,
    }


def cache_is_valid(path: Path, expected: dict[str, object]) -> bool:
    if not path.exists():
        return False
    try:
        with h5py.File(path, "r") as handle:
            if not handle.attrs.get("complete", False):
                return False
            saved = json.loads(handle.attrs["signature"])
            for key, value in expected.items():
                if key != "analysis_version" and saved.get(key) != value:
                    return False
        return True
    except (OSError, ValueError, KeyError, json.JSONDecodeError):
        return False


def validate_cache_counts(meta: dict[str, object], cache_path: Path) -> None:
    donor_table = pd.read_csv(REPO_ROOT / "results" / f"followup-between-subject-{meta['dataset']}-donors.csv")
    donor_table = donor_table[~donor_table["is_pooled"]]
    expected = np.zeros((len(meta["donor_ids"]), len(meta["cell_types"])), dtype=np.int64)
    donor_lookup = {donor: idx for idx, donor in enumerate(meta["donor_ids"])}
    ct_lookup = {ct: idx for idx, ct in enumerate(meta["cell_types"])}
    for row in donor_table.itertuples():
        expected[donor_lookup[str(row.donor_id)], ct_lookup[str(row.group)]] = int(row.n_cells_ge1500)
    observed = np.zeros_like(expected)
    with h5py.File(cache_path, "r") as handle:
        for ci in range(len(meta["cell_types"])):
            donors = handle[f"ct_{ci:03d}"]["donor_idx"][:]
            observed[:, ci] = np.bincount(donors, minlength=len(meta["donor_ids"]))
    if not np.array_equal(expected, observed):
        bad = np.argwhere(expected != observed)
        raise AssertionError(f"Cache donor cell counts differ at {bad[:10].tolist()}")
    log(f"Validated {meta['dataset']} cache against every donor-cell-type count")


def build_count_cache(meta: dict[str, object], projection: sp.csc_matrix, audit: dict[str, object]) -> Path:
    dataset = str(meta["dataset"])
    CACHE_ROOT.mkdir(parents=True, exist_ok=True)
    final_path = CACHE_ROOT / f"{dataset}_counts.h5"
    signature = cache_signature(meta, audit)
    temp_path = CACHE_ROOT / f"{dataset}_counts.tmp.h5"
    temp_path.unlink(missing_ok=True)
    n_raw_genes = len(meta["gene_symbols"])
    n_clean_genes = projection.shape[1]
    gate_total = {"nonfinite": 0, "negative": 0, "noninteger": 0}
    started = time.perf_counter()
    with h5py.File(meta["h5_path"], "r") as source, h5py.File(temp_path, "w") as target:
        group = source["layers"]["counts"] if meta["matrix_slot"] == "layers/counts" else source["X"]
        for ci, cell_type in enumerate(meta["cell_types"]):
            ct_group = target.create_group(f"ct_{ci:03d}")
            ct_group.attrs["cell_type"] = cell_type
            ct_group.create_dataset("data", shape=(0,), maxshape=(None,), dtype="float32", chunks=(1000000,))
            ct_group.create_dataset("indices", shape=(0,), maxshape=(None,), dtype="int32", chunks=(1000000,))
            ct_group.create_dataset("indptr", shape=(1,), maxshape=(None,), dtype="int64", chunks=(10000,))
            ct_group["indptr"][0] = 0
            ct_group.create_dataset("donor_idx", shape=(0,), maxshape=(None,), dtype="int32", chunks=(10000,))
            ct_group.create_dataset("umi", shape=(0,), maxshape=(None,), dtype="int32", chunks=(10000,))
        target.attrs["signature"] = json.dumps(signature, sort_keys=True)
        for row0 in range(0, meta["n_obs"], CHUNK_ROWS):
            row1 = min(row0 + CHUNK_ROWS, meta["n_obs"])
            local_keep = meta["keep_mask"][row0:row1]
            if not np.any(local_keep):
                continue
            raw = read_csr_chunk(group, row0, row1, n_raw_genes, meta["indptr"])[local_keep]
            gate = count_gate(raw.data)
            for key in gate_total:
                gate_total[key] += gate[key]
            if any(gate_total[key] for key in gate_total):
                raise ValueError(f"Raw count gate failed: {gate_total}")
            clean = (raw @ projection).tocsr()
            clean.eliminate_zeros()
            umi = np.rint(np.asarray(clean.sum(axis=1)).ravel()).astype(np.int64)
            eligible = umi >= MIN_DEPTH
            if not np.any(eligible):
                continue
            clean = clean[eligible]
            clean.data = np.rint(clean.data).astype(np.float32)
            clean.indices = clean.indices.astype(np.int32)
            donors = meta["cell_donor_idx"][row0:row1][local_keep][eligible]
            types = meta["cell_ct_idx"][row0:row1][local_keep][eligible]
            depths = umi[eligible]
            for ci in np.unique(types):
                selected = types == ci
                matrix = clean[selected]
                ct_group = target[f"ct_{int(ci):03d}"]
                old_cells = ct_group["donor_idx"].shape[0]
                old_nnz = ct_group["data"].shape[0]
                for name, values in (
                    ("data", matrix.data),
                    ("indices", matrix.indices),
                    ("donor_idx", donors[selected].astype(np.int32)),
                    ("umi", depths[selected].astype(np.int32)),
                ):
                    output = ct_group[name]
                    output.resize(output.shape[0] + len(values), axis=0)
                    output[-len(values):] = values
                ptr = matrix.indptr.astype(np.int64) + old_nnz
                ptr_out = ct_group["indptr"]
                ptr_out.resize(old_cells + matrix.shape[0] + 1, axis=0)
                ptr_out[old_cells + 1:] = ptr[1:]
            if (row0 // CHUNK_ROWS) % 20 == 0:
                log(f"{dataset} cache rows {row1}/{meta['n_obs']}")
        for ci in range(len(meta["cell_types"])):
            ptr = target[f"ct_{ci:03d}"]["indptr"]
            assert ptr[-1] == target[f"ct_{ci:03d}"]["data"].shape[0]
        target.attrs["gate"] = json.dumps(gate_total, sort_keys=True)
        target.attrs["complete"] = True
    os.replace(temp_path, final_path)
    log(f"Built {dataset} count cache in {(time.perf_counter() - started) / 60:.1f} minutes")
    return final_path


def resolve_count_cache(meta: dict[str, object], projection: sp.csc_matrix, audit: dict[str, object]) -> tuple[Path, dict[str, object]]:
    signature = cache_signature(meta, audit)
    own = CACHE_ROOT / f"{meta['dataset']}_counts.h5"
    prior = DATA_ROOT / "block2_gcl_cache" / f"{meta['dataset']}_counts.h5"
    for candidate in (own, prior):
        if cache_is_valid(candidate, signature):
            log(f"Reusing validated count cache {candidate}")
            validate_cache_counts(meta, candidate)
            return candidate, signature
    path = build_count_cache(meta, projection, audit)
    validate_cache_counts(meta, path)
    return path, signature


def aggregate_profiles(matrix: sp.csr_matrix, local_donors: np.ndarray, n_donors: int) -> np.ndarray:
    counts = np.bincount(local_donors, minlength=n_donors)
    if np.any(counts == 0):
        raise ValueError("Every donor must contribute cells")
    weights = 1.0 / counts[local_donors].astype(float)
    indicator = sp.csr_matrix(
        (weights, (local_donors, np.arange(len(local_donors)))),
        shape=(n_donors, len(local_donors)),
    )
    profiles = (indicator @ matrix).toarray().astype(np.float64, copy=False)
    if not np.allclose(profiles.sum(axis=1), 1.0, atol=1e-8, rtol=1e-8):
        raise AssertionError("Donor centroids are not proportions")
    return profiles


def batch_design(batch: np.ndarray) -> np.ndarray:
    levels = sorted(str(x) for x in np.unique(batch))
    return np.column_stack([np.ones(len(batch), dtype=float), *[(batch == level).astype(float) for level in levels[1:]]])


def batch_residualize(profiles: np.ndarray, batch: np.ndarray) -> np.ndarray:
    x = batch_design(batch)
    fitted = x @ np.linalg.lstsq(x, profiles, rcond=None)[0]
    adjusted = profiles - fitted + profiles.mean(axis=0, keepdims=True)
    adjusted = np.maximum(adjusted, 0.0)
    totals = adjusted.sum(axis=1, keepdims=True)
    if np.any(totals <= 0):
        raise ValueError("Batch residualization produced an empty profile")
    adjusted /= totals
    return adjusted


def pairwise_hellinger(profiles: np.ndarray) -> np.ndarray:
    totals = profiles.sum(axis=1)
    if np.any(profiles < 0) or np.any(totals <= 0):
        raise ValueError("Hellinger profiles must be nonnegative and nonempty")
    proportions = profiles / totals[:, None]
    roots = np.sqrt(proportions)
    affinity = roots @ roots.T
    distances = np.sqrt(np.maximum(0.0, 1.0 - affinity))
    distances = np.minimum(distances, 1.0)
    np.fill_diagonal(distances, 0.0)
    return distances


def mean_pairwise(distances: np.ndarray, mask: np.ndarray) -> float:
    selected = distances[np.ix_(mask, mask)]
    n = int(mask.sum())
    return float(selected.sum() / (n * (n - 1))) if n >= 2 else np.nan


def own_bin_distance(profiles: np.ndarray, mask: np.ndarray) -> float:
    if int(mask.sum()) == 0:
        return np.nan
    center = profiles[mask].mean(axis=0)
    return float(np.mean(hellinger_distance(profiles[mask], center)))


def permuted_difference(distances: np.ndarray, first_mask: np.ndarray, seed: int) -> dict[str, object]:
    n_first = int(first_mask.sum())
    n = len(first_mask)
    observed_first = mean_pairwise(distances, first_mask)
    observed_second = mean_pairwise(distances, ~first_mask)
    observed = observed_first - observed_second
    rng = np.random.default_rng(seed)
    labels = np.zeros((2 * N_PERM, n), dtype=np.float32)
    for perm in range(N_PERM):
        order = rng.permutation(n)
        labels[perm, order[:n_first]] = 1.0
        labels[N_PERM + perm] = 1.0 - labels[perm]
    products = labels @ distances.astype(np.float32)
    pair_sums = 0.5 * np.einsum("ij,ij->i", products, labels, dtype=np.float64)
    first_sums = pair_sums[:N_PERM]
    second_sums = pair_sums[N_PERM:]
    null = first_sums / (n_first * (n_first - 1)) - second_sums / ((n - n_first) * (n - n_first - 1))
    pvalue = float((1 + np.sum(np.abs(null) >= abs(observed) - 1e-12)) / (N_PERM + 1))
    null_sd = float(np.std(null, ddof=1))
    return {
        "difference": observed,
        "permutation_pvalue": pvalue,
        "permutation_null_mean": float(np.mean(null)),
        "permutation_null_sd": null_sd,
        "mde80": MDE80_Z * null_sd,
    }


def empty_comparison(
    dataset: str,
    cell_type: str,
    arm: str,
    comparison: str,
    bin1: str,
    bin2: str,
    n_genes: int,
    n_eligible: int,
    bins: np.ndarray,
    cell_counts: np.ndarray,
) -> dict[str, object]:
    n1 = int(np.sum(bins == bin1))
    n2 = int(np.sum(bins == bin2))
    return {
        "dataset": dataset,
        "cell_type": cell_type,
        "arm": arm,
        "comparison": comparison,
        "bin1": bin1,
        "bin2": bin2,
        "status": "EVALUATED" if min(n1, n2) >= MIN_DONORS_PER_BIN else "UNTESTABLE_BIN_FLOOR",
        "n_genes": n_genes,
        "n_donors_eligible": n_eligible,
        "n_donors": n1 + n2,
        "n_donors_bin1": n1,
        "n_donors_bin2": n2,
        "n_cells_bin1": int(cell_counts[bins == bin1].sum()),
        "n_cells_bin2": int(cell_counts[bins == bin2].sum()),
        "mean_pairwise_bin1": np.nan,
        "mean_pairwise_bin2": np.nan,
        "difference": np.nan,
        "mean_to_bin_centroid_bin1": np.nan,
        "mean_to_bin_centroid_bin2": np.nan,
        "secondary_difference": np.nan,
        "permutation_pvalue": np.nan,
        "q_bh": np.nan,
        "family": "old_young" if comparison == "OLD_vs_YOUNG" else "mid_secondary",
        "family_size": np.nan,
        "n_permutations": N_PERM if min(n1, n2) >= MIN_DONORS_PER_BIN else 0,
        "permutation_seed": stable_seed(PERM_SEED, dataset, cell_type, arm, comparison),
        "permutation_null_mean": np.nan,
        "permutation_null_sd": np.nan,
        "mde80": np.nan,
        "split_half_noise_floor": np.nan,
        "split_half_noise_floor_secondary": np.nan,
        "split_half_abs_change_bin1": np.nan,
        "split_half_abs_change_bin2": np.nan,
    }


def comparison_values(
    profiles: np.ndarray,
    bins: np.ndarray,
    row: dict[str, object],
    distances: np.ndarray | None = None,
    permute: bool = True,
) -> dict[str, object]:
    mask1 = bins == row["bin1"]
    mask2 = bins == row["bin2"]
    pair_mask = mask1 | mask2
    out = dict(row)
    if not np.any(pair_mask):
        out["status"] = "UNTESTABLE_BIN_FLOOR"
        out["n_permutations"] = 0
        return out
    sub_profiles = profiles[pair_mask]
    sub_mask1 = mask1[pair_mask]
    sub_distances = (
        pairwise_hellinger(sub_profiles)
        if distances is None
        else distances[np.ix_(pair_mask, pair_mask)]
    )
    out["mean_pairwise_bin1"] = mean_pairwise(sub_distances, sub_mask1)
    out["mean_pairwise_bin2"] = mean_pairwise(sub_distances, ~sub_mask1)
    out["difference"] = out["mean_pairwise_bin1"] - out["mean_pairwise_bin2"]
    out["mean_to_bin_centroid_bin1"] = own_bin_distance(sub_profiles, sub_mask1)
    out["mean_to_bin_centroid_bin2"] = own_bin_distance(sub_profiles, ~sub_mask1)
    out["secondary_difference"] = out["mean_to_bin_centroid_bin1"] - out["mean_to_bin_centroid_bin2"]
    if row["status"] == "EVALUATED" and permute:
        out.update(permuted_difference(sub_distances, sub_mask1, int(row["permutation_seed"])))
    return out


def process_cell_type(
    meta: dict[str, object],
    cache_path: Path,
    signature: dict[str, object],
    ct_idx: int,
    reuse_checkpoint: bool,
) -> list[dict[str, object]]:
    dataset = str(meta["dataset"])
    cell_type = str(meta["cell_types"][ct_idx])
    checkpoint = CACHE_ROOT / "checkpoints" / dataset / f"ct_{ct_idx:03d}.pkl"
    if reuse_checkpoint and checkpoint.exists():
        saved = pd.read_pickle(checkpoint)
        if saved.get("version") == ANALYSIS_VERSION and saved.get("signature") == signature:
            log(f"Reusing {dataset} {cell_type} Block 3 checkpoint")
            return saved["rows"]
    started = time.perf_counter()
    n_genes = int(signature["n_clean_genes"])
    with h5py.File(cache_path, "r") as handle:
        group = handle[f"ct_{ct_idx:03d}"]
        data = group["data"][:]
        indices = group["indices"][:]
        indptr = group["indptr"][:]
        donor_idx = group["donor_idx"][:]
        umi = group["umi"][:].astype(np.float64)
    counts_by_donor = np.bincount(donor_idx, minlength=len(meta["donor_ids"]))
    eligible_donors = np.flatnonzero(counts_by_donor >= MIN_CELLS)
    donor_meta = meta["donors_df"].iloc[eligible_donors].copy()
    keep_age = donor_meta["bin"].astype(str).isin(BINS).to_numpy()
    eligible_donors = eligible_donors[keep_age]
    donor_meta = donor_meta.iloc[np.flatnonzero(keep_age)].copy()
    bins = donor_meta["bin"].astype(str).to_numpy()
    batch = donor_meta["study" if dataset == "gao" else "pool_number"].astype(str).to_numpy()
    cell_counts = counts_by_donor[eligible_donors]
    rows = [
        empty_comparison(dataset, cell_type, arm, comp, bin1, bin2, n_genes, len(eligible_donors), bins, cell_counts)
        for arm in ARMS
        for comp, bin1, bin2 in COMPARISONS
    ]
    row_lookup = {(str(r["arm"]), str(r["comparison"])): r for r in rows}
    if len(eligible_donors) == 0:
        checkpoint.parent.mkdir(parents=True, exist_ok=True)
        pd.to_pickle({"version": ANALYSIS_VERSION, "signature": signature, "rows": rows}, checkpoint)
        return rows

    selected = np.isin(donor_idx, eligible_donors)
    matrix = sp.csr_matrix((data, indices, indptr), shape=(len(donor_idx), n_genes))[selected].tocsr()
    donor_idx = donor_idx[selected]
    umi = umi[selected]
    del data, indices, indptr, selected
    order = np.argsort(donor_idx, kind="stable")
    matrix = matrix[order].tocsr()
    donor_idx = donor_idx[order]
    umi = umi[order]
    matrix.data = matrix.data.astype(np.float64) * np.repeat(1.0 / umi, np.diff(matrix.indptr))
    local_donors = np.searchsorted(eligible_donors, donor_idx).astype(np.int32)
    boundaries = np.r_[0, np.cumsum(counts_by_donor[eligible_donors])]
    full = aggregate_profiles(matrix, local_donors, len(eligible_donors))
    full_profiles = {"raw": full, "batch_residualized": batch_residualize(full, batch)}
    full_distances = {arm: pairwise_hellinger(profiles) for arm, profiles in full_profiles.items()}
    for arm in ARMS:
        for comp, _, _ in COMPARISONS:
            key = (arm, comp)
            row_lookup[key] = comparison_values(
                full_profiles[arm], bins, row_lookup[key], distances=full_distances[arm], permute=True
            )
    del full_distances
    noise = {(arm, comp): [] for arm in ARMS for comp, _, _ in COMPARISONS}
    noise_secondary = {(arm, comp): [] for arm in ARMS for comp, _, _ in COMPARISONS}
    noise_bin = {(arm, comp, b): [] for arm in ARMS for comp, _, _ in COMPARISONS for b in BINS}

    for seed in SPLIT_SEEDS:
        half = np.empty(len(donor_idx), dtype=np.int8)
        for pos, global_donor in enumerate(eligible_donors):
            start, stop = boundaries[pos:pos + 2]
            n = stop - start
            rng = np.random.default_rng(stable_seed(seed, dataset, cell_type, int(global_donor)))
            permutation = rng.permutation(n)
            local_half = np.ones(n, dtype=np.int8)
            local_half[permutation[: n // 2]] = 0
            half[start:stop] = local_half
        half_profiles: dict[str, dict[str, np.ndarray]] = {"raw": {}, "batch_residualized": {}}
        for half_idx in (0, 1):
            mask = half == half_idx
            profiles = aggregate_profiles(matrix[mask], local_donors[mask], len(eligible_donors))
            half_profiles["raw"][half_idx] = profiles
            half_profiles["batch_residualized"][half_idx] = batch_residualize(profiles, batch)
        for arm in ARMS:
            half_values = {}
            for half_idx in (0, 1):
                profiles = half_profiles[arm][half_idx]
                distances = pairwise_hellinger(profiles)
                half_values[half_idx] = {}
                for comp, _, _ in COMPARISONS:
                    base = empty_comparison(dataset, cell_type, arm, comp, comp.split("_vs_")[0], comp.split("_vs_")[1], n_genes, len(eligible_donors), bins, cell_counts)
                    half_values[half_idx][comp] = comparison_values(profiles, bins, base, distances=distances, permute=False)
            for comp, _, _ in COMPARISONS:
                a = half_values[0][comp]
                b = half_values[1][comp]
                noise[(arm, comp)].append(abs(a["difference"] - b["difference"]))
                noise_secondary[(arm, comp)].append(abs(a["secondary_difference"] - b["secondary_difference"]))
                for field, bin_name in (("bin1", a["bin1"]), ("bin2", a["bin2"])):
                    noise_bin[(arm, comp, bin_name)].append(abs(a[f"mean_pairwise_{field}"] - b[f"mean_pairwise_{field}"]))
        del half_profiles, half_values
        gc.collect()

    for key, row in row_lookup.items():
        values = np.asarray(noise[key], dtype=float)
        secondary_values = np.asarray(noise_secondary[key], dtype=float)
        row["split_half_noise_floor"] = float(np.nanmean(values)) if np.isfinite(values).any() else np.nan
        row["split_half_noise_floor_secondary"] = float(np.nanmean(secondary_values)) if np.isfinite(secondary_values).any() else np.nan
        for field, bin_name in (("split_half_abs_change_bin1", row["bin1"]), ("split_half_abs_change_bin2", row["bin2"])):
            bin_values = np.asarray(noise_bin[(*key, bin_name)], dtype=float)
            row[field] = float(np.nanmean(bin_values)) if np.isfinite(bin_values).any() else np.nan
    rows = list(row_lookup.values())
    checkpoint.parent.mkdir(parents=True, exist_ok=True)
    pd.to_pickle({"version": ANALYSIS_VERSION, "signature": signature, "rows": rows}, checkpoint)
    log(f"Finished {dataset} {cell_type} in {(time.perf_counter() - started) / 60:.1f} minutes")
    del matrix, full, full_profiles
    gc.collect()
    return rows


def apply_bh(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out["q_bh"] = np.nan
    out["family_size"] = np.nan
    evaluated = out["status"] == "EVALUATED"
    for _, indices in out[evaluated].groupby(["dataset", "arm", "family"], sort=False).groups.items():
        idx = np.asarray(list(indices))
        pvalues = out.loc[idx, "permutation_pvalue"].to_numpy(dtype=float)
        valid = np.isfinite(pvalues)
        if not np.any(valid):
            continue
        idx = idx[valid]
        out.loc[idx, "q_bh"] = multipletests(pvalues[valid], method="fdr_bh")[1]
        out.loc[idx, "family_size"] = int(valid.sum())
    return out


def validate_results(results: pd.DataFrame, meta_by_dataset: dict[str, dict[str, object]]) -> None:
    expected = sum(len(meta["cell_types"]) * len(ARMS) * len(COMPARISONS) for meta in meta_by_dataset.values())
    assert len(results) == expected
    keys = ["dataset", "cell_type", "arm", "comparison"]
    assert not results.duplicated(keys).any()
    evaluated = results[results["status"] == "EVALUATED"]
    for column in ("difference", "permutation_pvalue", "q_bh", "mde80", "split_half_noise_floor"):
        assert np.isfinite(evaluated[column].to_numpy(dtype=float)).all(), column
    assert (evaluated["n_donors_bin1"] >= MIN_DONORS_PER_BIN).all()
    assert (evaluated["n_donors_bin2"] >= MIN_DONORS_PER_BIN).all()
    assert (evaluated["n_permutations"] == N_PERM).all()
    assert evaluated["permutation_pvalue"].between(0, 1).all()
    assert evaluated["q_bh"].between(0, 1).all()
    assert (evaluated["mde80"] >= 0).all()
    for dataset, meta in meta_by_dataset.items():
        subset = results[results["dataset"] == dataset]
        assert set(subset["cell_type"]) == set(meta["cell_types"])


def cell_type_verdict(primary: pd.Series, raw: pd.Series) -> str:
    if primary["status"] != "EVALUATED":
        return "UNTESTABLE"
    if raw["status"] == "EVALUATED" and primary["difference"] * raw["difference"] < 0:
        return "INCONCLUSIVE"
    if primary["difference"] > 0 and primary["q_bh"] < 0.05 and primary["difference"] > primary["split_half_noise_floor"]:
        return "PROVES"
    if (primary["q_bh"] >= 0.05 or primary["difference"] <= 0) and primary["mde80"] <= primary["split_half_noise_floor"]:
        return "DOES NOT PROVE"
    return "INCONCLUSIVE"


def fmt(value: object) -> str:
    if value is None or (isinstance(value, (float, np.floating)) and not np.isfinite(value)):
        return "NA"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.5g}"
    return str(value).replace("|", "/")


def markdown_table(frame: pd.DataFrame, columns: list[tuple[str, str]]) -> str:
    lines = ["| " + " | ".join(label for _, label in columns) + " |"]
    lines.append("| " + " | ".join("---" for _ in columns) + " |")
    for _, row in frame.iterrows():
        lines.append("| " + " | ".join(fmt(row[key]) for key, _ in columns) + " |")
    return "\n".join(lines)


def write_report(results: pd.DataFrame) -> None:
    spec = REPORT_MD.read_text(encoding="utf-8").split("## Results", 1)[0]
    spec = spec.replace("Status: analysis pending.", "Status: complete.")
    primary = results[(results["arm"] == "batch_residualized") & (results["comparison"] == "OLD_vs_YOUNG")].copy()
    raw = results[(results["arm"] == "raw") & (results["comparison"] == "OLD_vs_YOUNG")].copy()
    raw_lookup = raw.set_index(["dataset", "cell_type"])
    summary_rows = []
    for _, row in primary.iterrows():
        key = (row["dataset"], row["cell_type"])
        raw_row = raw_lookup.loc[key]
        summary_rows.append({
            "dataset": row["dataset"],
            "cell_type": row["cell_type"],
            "n_young": row["n_donors_bin2"],
            "n_old": row["n_donors_bin1"],
            "young": row["mean_pairwise_bin2"],
            "old": row["mean_pairwise_bin1"],
            "difference": row["difference"],
            "pvalue": row["permutation_pvalue"],
            "q": row["q_bh"],
            "noise": row["split_half_noise_floor"],
            "mde80": row["mde80"],
            "raw_difference": raw_row["difference"],
            "raw_q": raw_row["q_bh"],
            "raw_noise": raw_row["split_half_noise_floor"],
            "raw_mde80": raw_row["mde80"],
            "secondary_difference": row["secondary_difference"],
            "verdict": cell_type_verdict(row, raw_row),
        })
    summary = pd.DataFrame(summary_rows)
    dataset_verdicts = {}
    lines = [
        "## Results",
        "",
        "All pairwise effects are first-bin minus second-bin Hellinger distance differences. The batch-residualized arm is primary; raw is the sensitivity. NA denotes an UNTESTABLE or not applicable field.",
        "",
        "### OLD vs YOUNG primary results and raw sensitivity",
        "",
        markdown_table(summary, [
            ("dataset", "Dataset"), ("cell_type", "Cell type"), ("n_young", "N YOUNG"), ("n_old", "N OLD"),
            ("young", "Residualized YOUNG"), ("old", "Residualized OLD"), ("difference", "Residualized difference"),
            ("pvalue", "p"), ("q", "BH q"), ("noise", "Noise floor"), ("mde80", "MDE80"),
            ("raw_difference", "Raw difference"), ("raw_q", "Raw q"), ("raw_noise", "Raw noise"), ("raw_mde80", "Raw MDE80"),
            ("secondary_difference", "Own-bin secondary difference"), ("verdict", "Verdict"),
        ]),
        "",
        "### MID secondary comparisons",
        "",
        "The MID BH family combines OLD vs MID and MID vs YOUNG within each dataset and arm. Permutation p-values and q-values apply to mean pairwise distance, not the descriptive own-bin-centroid outcome.",
        "",
    ]
    mid = results[results["comparison"] != "OLD_vs_YOUNG"].copy()
    lines.append(markdown_table(mid, [
        ("dataset", "Dataset"), ("cell_type", "Cell type"), ("arm", "Arm"), ("comparison", "Comparison"),
        ("n_donors_bin1", "N bin 1"), ("n_donors_bin2", "N bin 2"),
        ("mean_pairwise_bin1", "Bin 1 pairwise"), ("mean_pairwise_bin2", "Bin 2 pairwise"),
        ("difference", "Difference"), ("permutation_pvalue", "p"), ("q_bh", "BH q"),
        ("split_half_noise_floor", "Noise floor"), ("mde80", "MDE80"),
        ("secondary_difference", "Own-bin secondary difference"), ("status", "Status"),
    ]))
    lines.extend(["", "### Verdicts", ""])
    for dataset in ("gao", "yazar"):
        subset = summary[summary["dataset"] == dataset]
        evaluated = subset[subset["verdict"] != "UNTESTABLE"]
        if (evaluated["verdict"] == "PROVES").any():
            verdict = "PROVES"
        elif not evaluated.empty and (evaluated["verdict"] == "DOES NOT PROVE").all():
            verdict = "DOES NOT PROVE"
        else:
            verdict = "INCONCLUSIVE"
        dataset_verdicts[dataset] = verdict
        proved = subset.loc[subset["verdict"] == "PROVES", "cell_type"].tolist()
        inconclusive = subset.loc[subset["verdict"] == "INCONCLUSIVE", "cell_type"].tolist()
        untestable = subset.loc[subset["verdict"] == "UNTESTABLE", "cell_type"].tolist()
        lines.append(f"- {dataset.capitalize()}: {verdict}. PROVES {len(proved)}, DOES NOT PROVE {int((subset['verdict'] == 'DOES NOT PROVE').sum())}, INCONCLUSIVE {len(inconclusive)}, UNTESTABLE {len(untestable)} cell types.")
        if proved:
            lines.append(f"  Proven units: {', '.join(proved)}.")
        if inconclusive:
            lines.append(f"  Inconclusive units: {', '.join(inconclusive)}.")
        if untestable:
            lines.append(f"  Untestable units: {', '.join(untestable)}.")
    lines.extend(["", "### Distinction from the Nov 17 distance-variance test", ""])
    lines.append("The Nov 17 test reduced each donor to one distance to a young reference and then compared the variance of those donor distances. Block 3 instead measures dispersion within each age bin directly from pairwise distances among donor expression centroids. It does not use a young reference and does not analyze the variance of one distance per donor. These are different estimands.")
    for dataset in ("gao", "yazar"):
        verdict = dataset_verdicts[dataset]
        if verdict == "PROVES":
            relation = "does not agree with the Nov 17 batch-adjusted null, although the estimands differ"
        elif verdict == "DOES NOT PROVE":
            relation = "agrees with the Nov 17 batch-adjusted null at the level of direction and inference"
        else:
            relation = "does not clearly agree or disagree with the Nov 17 batch-adjusted null because this test is inconclusive"
        lines.append(f"- {dataset.capitalize()}: this profile-cloud test {relation}. The Nov 17 headline was Gao OLD to YOUNG distance-variance ratio 0.79, p 0.85; Yazar raw ratio 1.62, p 1.9e-07, changing to pool-residualized ratio 1.25, p 0.43.")
    family_parts = []
    for (dataset, arm, family), group in results[results["status"] == "EVALUATED"].groupby(["dataset", "arm", "family"]):
        family_parts.append(f"{dataset} {arm} {family} m={len(group)}")
    lines.extend([
        "",
        "### Reproducibility and family sizes",
        "",
        f"Permutation base seed {PERM_SEED}; split-half seeds {', '.join(map(str, SPLIT_SEEDS))}. Exact comparison seeds are in the CSV. All 2000 permutations resample donor labels, never cells.",
        "Evaluated BH family sizes: " + "; ".join(family_parts) + ".",
        "No primary cells were downsampled. Full-depth within-cell proportions use only cells passing the 1500-count eligibility filter.",
        "",
    ])
    write_clean_markdown(spec + "\n".join(lines), REPORT_MD)


def self_test() -> None:
    rng = np.random.default_rng(123)
    profiles = rng.dirichlet(np.ones(6), size=8)
    distances = pairwise_hellinger(profiles)
    for i in range(len(profiles)):
        for j in range(i + 1, len(profiles)):
            expected = float(hellinger_distance(profiles[j:j + 1], profiles[i])[0])
            assert np.isclose(distances[i, j], expected, atol=1e-10), (i, j, distances[i, j], expected)
    counts = sp.csr_matrix(np.array([[1, 1], [3, 0]], dtype=float))
    normalized = counts.multiply(1.0 / np.asarray(counts.sum(axis=1))).tocsr()
    centroid = aggregate_profiles(normalized, np.array([0, 0], dtype=np.int32), 1)
    assert np.allclose(centroid, [[0.75, 0.25]])
    batch = np.array(["a", "a", "b", "b"])
    safe_profiles = np.array([[0.4, 0.6], [0.42, 0.58], [0.6, 0.4], [0.62, 0.38]])
    adjusted = batch_residualize(safe_profiles, batch)
    assert np.allclose(adjusted[batch == "a"].mean(axis=0), adjusted[batch == "b"].mean(axis=0))
    assert np.allclose(adjusted.sum(axis=1), 1.0)
    first = np.array([True, True, False, False])
    small_distances = pairwise_hellinger(profiles[:4])
    a = permuted_difference(small_distances, first, 777)
    b = permuted_difference(small_distances, first, 777)
    assert a == b
    log("Self test passes for Hellinger distances, equal-cell centroids, residualization, and permutations")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--no-reuse-checkpoints", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if args.report_only:
        write_report(pd.read_csv(RESULT_CSV))
        return
    if not REPORT_MD.exists() or "## Exact analysis specification" not in REPORT_MD.read_text(encoding="utf-8"):
        raise RuntimeError("The pre-analysis specification report must exist")
    self_test()
    rows = []
    meta_by_dataset = {}
    for dataset, path in (("gao", GAO_H5AD), ("yazar", YAZAR_H5AD)):
        meta = reconcile_gao_metadata(path) if dataset == "gao" else reconcile_yazar_metadata(path, YAZAR_CACHE)
        projection_info = build_human_gene_projection(meta["gene_symbols"])
        cache_path, signature = resolve_count_cache(meta, projection_info["proj"], projection_info["audit"])
        meta_by_dataset[dataset] = meta
        for ct_idx in range(len(meta["cell_types"])):
            rows.extend(process_cell_type(meta, cache_path, signature, ct_idx, not args.no_reuse_checkpoints))
        del projection_info
        gc.collect()
    results = apply_bh(pd.DataFrame(rows))
    validate_results(results, meta_by_dataset)
    preferred = [
        "dataset", "cell_type", "arm", "comparison", "bin1", "bin2", "status",
        "mean_pairwise_bin1", "mean_pairwise_bin2", "difference", "permutation_pvalue", "q_bh",
        "family", "family_size", "n_donors", "n_donors_bin1", "n_donors_bin2", "n_cells_bin1", "n_cells_bin2",
        "n_donors_eligible", "n_genes", "split_half_noise_floor", "mde80",
        "mean_to_bin_centroid_bin1", "mean_to_bin_centroid_bin2", "secondary_difference", "split_half_noise_floor_secondary",
        "split_half_abs_change_bin1", "split_half_abs_change_bin2", "permutation_null_mean", "permutation_null_sd",
        "n_permutations", "permutation_seed",
    ]
    results = results[preferred].sort_values(["dataset", "cell_type", "arm", "comparison"]).reset_index(drop=True)
    write_clean_csv(results, RESULT_CSV)
    write_report(results)
    log(f"Wrote {RESULT_CSV} and {REPORT_MD}")


if __name__ == "__main__":
    main()
