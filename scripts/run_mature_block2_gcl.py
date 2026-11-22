#!/usr/bin/env python3
"""Mature Block 2 donor-level gene coordination in Gao heart and Yazar blood.

Mature-to-old rerun of scripts/run_followup_11pt_block2_gcl.py.
Only change vs the frozen original: donors aged under 30 are excluded,
detected-gene sets and fixed gene pairs are rebuilt from mature (30-plus)
donors only, and the group comparison baseline YOUNG (20-39) is replaced by
MATURE (30-39). Every other frozen rule (eligibility floors, Spearman GCL,
primary/depth/rarefied/unadjusted models, 2000 permutations with seed
20261118, split-half floors, BH families, verdict rules, MDE80) is
unchanged. The original Block 2 count cache is reused: its per-cell
projected counts are donor-local, and detected-gene derivation, pair
sampling, GCL scoring, and all models in this script use mature donors only.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import gc
import hashlib
import h5py
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import sys
import time

for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"

import numpy as np
import pandas as pd
from scipy import sparse as sp
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.run_followup_between_subject_gao_yazar import (
    build_human_gene_projection,
    read_csr_chunk,
    reconcile_gao_metadata,
    reconcile_yazar_metadata,
)
from scripts.run_followup_de_man_entropy import multinomial_downsample
from src.followup_muscle import count_gate
from src.followup_organs import write_clean_csv, write_clean_markdown

RESULT_CSV = REPO_ROOT / "results" / "followup-mature-block2-gcl.csv"
REPORT_MD = REPO_ROOT / "reports" / "followup-mature-block2-gcl.md"
RUNTIME_CSV = REPO_ROOT / "results" / "followup-mature-block2-runtime.csv"
DATA_ROOT = REPO_ROOT.parent / f"{REPO_ROOT.name}-data" / "tmp_organs"
CACHE_ROOT = DATA_ROOT / "block2_gcl_cache"
CHECKPOINT_ROOT = DATA_ROOT / "block2_gcl_mature_cache"
ANALYSIS_VERSION = "2026-11-21-mature-v1"
# Count-cache contents are donor-local per-cell projections identical to the
# original 2026-11-20 run, so the cache signature keeps the original version
# string; the mature analysis version must not invalidate cache reuse.
CACHE_VERSION = "2026-11-20-v1"
MIN_FREE_GIB = 20.0
MIN_CELLS = 20
MIN_DONORS = 6
MIN_DEPTH = 1500
PAIR_SEED = 20261117
PERM_SEED = 20261118
RAREFY_SEED = 20261119
SPLIT_SEEDS = tuple(20261116 + i for i in range(5))
MAX_SAMPLED_PAIRS = 20000
N_PERM = 2000
CHUNK_ROWS = 8000

_WORK_FULL = None
_WORK_RARE = None
_WORK_LEFT = None
_WORK_RIGHT = None
_WORK_ROWS = None
_WORK_DATASET_CODE = None
_WORK_CT_CODE = None


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def mem_free_gib() -> float:
    fields = {}
    for line in Path("/proc/meminfo").read_text().splitlines():
        key, value = line.split(":", 1)
        fields[key] = float(value.strip().split()[0])
    return fields["MemFree"] / 1024.0 / 1024.0


def release_clean_file_cache() -> None:
    paths = [
        DATA_ROOT / "organ17_gao_heart" / "GSE290367_integrated_RNA_adata.h5ad",
        DATA_ROOT / "hcata" / "yazar_1e44db10.h5ad",
        *CACHE_ROOT.glob("*_counts.h5"),
    ]
    for path in paths:
        if not path.exists():
            continue
        fd = os.open(path, os.O_RDONLY)
        try:
            os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
        finally:
            os.close(fd)


def wait_for_memory(label: str) -> None:
    while True:
        free = mem_free_gib()
        if free >= MIN_FREE_GIB:
            log(f"Memory guard passes for {label}: {free:.2f} GiB free")
            return
        release_clean_file_cache()
        free = mem_free_gib()
        if free >= MIN_FREE_GIB:
            log(f"Memory guard passes for {label} after releasing clean file cache: {free:.2f} GiB free")
            return
        log(f"Memory guard waits for {label}: {free:.2f} GiB free, need {MIN_FREE_GIB:.2f} GiB")
        time.sleep(60)


def stable_seed(base: int, *parts: object) -> int:
    payload = (str(base) + "|" + "|".join(str(p) for p in parts)).encode("utf-8")
    return int.from_bytes(hashlib.blake2b(payload, digest_size=8).digest(), "little") & 0x7FFFFFFFFFFFFFFF


def fixed_pairs(n_genes: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    if n_genes < 2:
        return np.array([], dtype=np.int32), np.array([], dtype=np.int32)
    if n_genes <= 2000:
        left, right = np.triu_indices(n_genes, 1)
        return left.astype(np.int32), right.astype(np.int32)
    n_pairs = min(MAX_SAMPLED_PAIRS, n_genes * (n_genes - 1) // 2)
    rng = np.random.default_rng(seed)
    selected = set()
    while len(selected) < n_pairs:
        size = max(1024, 2 * (n_pairs - len(selected)))
        a = rng.integers(0, n_genes, size=size, dtype=np.int64)
        b = rng.integers(0, n_genes, size=size, dtype=np.int64)
        left = np.minimum(a, b)
        right = np.maximum(a, b)
        valid = left != right
        for pair_key in (left[valid] * n_genes + right[valid]).tolist():
            selected.add(pair_key)
            if len(selected) == n_pairs:
                break
    keys = np.array(sorted(selected), dtype=np.int64)
    left = keys // n_genes
    right = keys % n_genes
    return left.astype(np.int32), right.astype(np.int32)


def gcl_from_matrix(matrix: sp.csr_matrix, left: np.ndarray, right: np.ndarray) -> float:
    n_cells, n_score_genes = matrix.shape
    if n_cells < 2 or len(left) == 0:
        return float("nan")
    csc = matrix.tocsc()
    ranks = np.empty((n_cells, n_score_genes), dtype=np.float32)
    values = np.zeros(n_cells, dtype=np.float64)
    for gene in range(n_score_genes):
        values.fill(0.0)
        start, stop = csc.indptr[gene:gene + 2]
        values[csc.indices[start:stop]] = csc.data[start:stop]
        ranks[:, gene] = stats.rankdata(values, method="average").astype(np.float32)
    ranks -= ranks.mean(axis=0, dtype=np.float64).astype(np.float32)
    norms = np.sqrt(np.einsum("ij,ij->j", ranks, ranks, dtype=np.float64))
    pair_valid = (norms[left] > 0.0) & (norms[right] > 0.0)
    if not np.any(pair_valid):
        return float("nan")
    total = 0.0
    n_valid = 0
    batch_size = max(1, 1000000 // max(n_cells, 1))
    for start in range(0, len(left), batch_size):
        stop = min(start + batch_size, len(left))
        li = left[start:stop]
        ri = right[start:stop]
        valid = pair_valid[start:stop]
        if not np.any(valid):
            continue
        li = li[valid]
        ri = ri[valid]
        a = ranks[:, li]
        b = ranks[:, ri]
        numerator = np.einsum("ij,ij->j", a, b, dtype=np.float64)
        corr = numerator / (norms[li] * norms[ri])
        total += float(np.sum(corr))
        n_valid += len(corr)
    return total / n_valid if n_valid else float("nan")


def _worker_init() -> None:
    gc.freeze()


def _donor_gcl(task: tuple[int, str]) -> dict[str, object]:
    donor_idx, donor_id = task
    rows = _WORK_ROWS[donor_idx]
    full = _WORK_FULL[rows]
    rare = _WORK_RARE[rows]
    result = {
        "donor_id": donor_id,
        "donor_idx": donor_idx,
        "n_cells": len(rows),
        "gcl_full": gcl_from_matrix(full, _WORK_LEFT, _WORK_RIGHT),
        "gcl_rarefied_1500": gcl_from_matrix(rare, _WORK_LEFT, _WORK_RIGHT),
    }
    n_rows = len(rows)
    for seed_idx, seed in enumerate(SPLIT_SEEDS):
        rng = np.random.default_rng(stable_seed(seed, _WORK_DATASET_CODE, _WORK_CT_CODE, donor_idx))
        order = rng.permutation(n_rows)
        half_a = order[: n_rows // 2]
        half_b = order[n_rows // 2 :]
        gcl_a = gcl_from_matrix(full[half_a], _WORK_LEFT, _WORK_RIGHT)
        gcl_b = gcl_from_matrix(full[half_b], _WORK_LEFT, _WORK_RIGHT)
        result[f"split_abs_seed_{seed_idx}"] = abs(gcl_a - gcl_b)
    return result


def cache_signature(meta: dict[str, object], projection_audit: dict[str, object]) -> dict[str, object]:
    path = Path(meta["h5_path"])
    stat = path.stat()
    return {
        "analysis_version": CACHE_VERSION,
        "dataset": meta["dataset"],
        "input_path": str(path),
        "input_size": stat.st_size,
        "input_mtime_ns": stat.st_mtime_ns,
        "matrix_slot": meta["matrix_slot"],
        "n_obs": meta["n_obs"],
        "n_raw_genes": len(meta["gene_symbols"]),
        "n_clean_genes": projection_audit["n_genes_retained"],
        "technical_excluded": projection_audit["n_technical_excluded"],
        "min_depth": MIN_DEPTH,
    }


def build_count_cache(meta: dict[str, object], projection: sp.csc_matrix, audit: dict[str, object]) -> Path:
    dataset = str(meta["dataset"])
    CACHE_ROOT.mkdir(parents=True, exist_ok=True)
    final_path = CACHE_ROOT / f"{dataset}_counts.h5"
    signature = cache_signature(meta, audit)
    if final_path.exists():
        with h5py.File(final_path, "r") as handle:
            if handle.attrs.get("complete", False) and json.loads(handle.attrs["signature"]) == signature:
                log(f"Reusing validated {dataset} count cache")
                validate_cache_counts(meta, final_path)
                return final_path
    wait_for_memory(f"building {dataset} count cache")
    temp_path = CACHE_ROOT / f"{dataset}_counts.tmp.h5"
    temp_path.unlink(missing_ok=True)
    n_donors = len(meta["donor_ids"])
    n_cts = len(meta["cell_types"])
    keep_mask = meta["keep_mask"]
    cell_donor = meta["cell_donor_idx"]
    cell_ct = meta["cell_ct_idx"]
    indptr = meta["indptr"]
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
            local_keep = keep_mask[row0:row1]
            if not np.any(local_keep):
                continue
            raw = read_csr_chunk(group, row0, row1, n_raw_genes, indptr)[local_keep]
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
            donors = cell_donor[row0:row1][local_keep][eligible]
            types = cell_ct[row0:row1][local_keep][eligible]
            depths = umi[eligible]
            for ci in np.unique(types):
                selected = types == ci
                matrix = clean[selected]
                ct_group = target[f"ct_{int(ci):03d}"]
                old_cells = ct_group["donor_idx"].shape[0]
                new_cells = matrix.shape[0]
                old_nnz = ct_group["data"].shape[0]
                new_nnz = matrix.nnz
                for name, values in (
                    ("data", matrix.data),
                    ("indices", matrix.indices),
                    ("donor_idx", donors[selected].astype(np.int32)),
                    ("umi", depths[selected].astype(np.int32)),
                ):
                    dataset_out = ct_group[name]
                    dataset_out.resize(dataset_out.shape[0] + len(values), axis=0)
                    dataset_out[-len(values):] = values
                ptr = matrix.indptr.astype(np.int64) + old_nnz
                ptr_out = ct_group["indptr"]
                ptr_out.resize(old_cells + new_cells + 1, axis=0)
                ptr_out[old_cells + 1:] = ptr[1:]
            if (row0 // CHUNK_ROWS) % 20 == 0:
                log(f"{dataset} cache rows {row1}/{meta['n_obs']}")
        for ci in range(n_cts):
            ptr = target[f"ct_{ci:03d}"]["indptr"]
            assert ptr[-1] == target[f"ct_{ci:03d}"]["data"].shape[0]
        target.attrs["gate"] = json.dumps(gate_total, sort_keys=True)
        target.attrs["complete"] = True
    os.replace(temp_path, final_path)
    elapsed = time.perf_counter() - started
    log(f"Built {dataset} count cache in {elapsed / 60:.1f} minutes")
    validate_cache_counts(meta, final_path)
    return final_path


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


def checkpoint_path(dataset: str, ct_idx: int) -> Path:
    return CHECKPOINT_ROOT / "checkpoints" / dataset / f"ct_{ct_idx:03d}.pkl"


def mature_donor_positions(meta: dict[str, object]) -> np.ndarray:
    """Positions into meta donor_ids of donors aged 30-plus (mature only)."""
    ages = meta["donors_df"]["age_years"].astype(float).to_numpy()
    mature = np.flatnonzero(ages >= 30.0)
    log(
        f"Mature restriction {meta['dataset']}: "
        f"{len(mature)} of {len(ages)} donors at 30-plus "
        f"({int(((ages >= 30.0) & (ages <= 39.0)).sum())} MATURE, "
        f"{int(((ages >= 40.0) & (ages <= 59.0)).sum())} MID, "
        f"{int((ages >= 60.0).sum())} OLD)"
    )
    return mature


def process_cell_type(
    meta: dict[str, object],
    cache_path: Path,
    ct_idx: int,
    max_workers: int,
    reuse_checkpoints: bool,
) -> tuple[dict[str, object], pd.DataFrame]:
    dataset = str(meta["dataset"])
    cell_type = str(meta["cell_types"][ct_idx])
    checkpoint = checkpoint_path(dataset, ct_idx)
    if reuse_checkpoints and checkpoint.exists():
        saved = pd.read_pickle(checkpoint)
        if saved["version"] == ANALYSIS_VERSION:
            log(f"Reusing {dataset} {cell_type} checkpoint")
            return saved["info"], saved["donors"]
    started = time.perf_counter()
    wait_for_memory(f"{dataset} {cell_type} scoring")
    with h5py.File(cache_path, "r") as handle:
        group = handle[f"ct_{ct_idx:03d}"]
        n_clean_genes = int(json.loads(handle.attrs["signature"])["n_clean_genes"])
        data = group["data"][:]
        indices = group["indices"][:]
        indptr = group["indptr"][:]
        donor_idx = group["donor_idx"][:]
        umi = group["umi"][:].astype(np.float64)
    matrix = sp.csr_matrix((data, indices, indptr), shape=(len(donor_idx), n_clean_genes))
    del data, indices, indptr
    mature_positions = mature_donor_positions(meta)
    mature_keep = np.isin(donor_idx, mature_positions)
    matrix = matrix[mature_keep].tocsr()
    donor_idx = donor_idx[mature_keep]
    umi = umi[mature_keep]
    del mature_keep
    detected_counts = np.bincount(matrix.indices, minlength=n_clean_genes)
    detected = np.flatnonzero(detected_counts >= math.ceil(0.10 * len(donor_idx)))
    pair_seed = stable_seed(PAIR_SEED, dataset, cell_type)
    left, right = fixed_pairs(len(detected), pair_seed)
    score_local = np.unique(np.concatenate([left, right]))
    score_genes = detected[score_local]
    local_to_score = np.full(len(detected), -1, dtype=np.int32)
    local_to_score[score_local] = np.arange(len(score_local), dtype=np.int32)
    work_left = local_to_score[left]
    work_right = local_to_score[right]
    counts_by_donor = np.bincount(donor_idx, minlength=len(meta["donor_ids"]))
    eligible_donors = np.flatnonzero(counts_by_donor >= MIN_CELLS)
    info = {
        "dataset": dataset,
        "cell_type": cell_type,
        "n_cells_pooled": int(len(donor_idx)),
        "n_genes": int(len(detected)),
        "n_pairs": int(len(left)),
        "pair_seed": int(pair_seed),
        "n_donors_ge20": int(len(eligible_donors)),
        "split_half_noise_floor": np.nan,
    }
    if len(eligible_donors) < MIN_DONORS or len(left) == 0:
        info["status"] = "UNTESTABLE"
        donors_out = pd.DataFrame(columns=["donor_id", "donor_idx", "n_cells", "gcl_full", "gcl_rarefied_1500"])
        elapsed = time.perf_counter() - started
        info["elapsed_seconds"] = elapsed
        checkpoint.parent.mkdir(parents=True, exist_ok=True)
        pd.to_pickle({"version": ANALYSIS_VERSION, "info": info, "donors": donors_out}, checkpoint)
        log(f"{dataset} {cell_type}: UNTESTABLE with {len(eligible_donors)} eligible donors and {len(detected)} genes")
        return info, donors_out
    eligible_mask = np.isin(donor_idx, eligible_donors)
    matrix = matrix[eligible_mask]
    donor_idx = donor_idx[eligible_mask]
    umi = umi[eligible_mask]
    score_matrix = matrix[:, score_genes].tocsr()
    score_matrix = score_matrix.multiply(1.0 / umi[:, None]).tocsr()
    rare_seed = stable_seed(RAREFY_SEED, dataset, cell_type)
    rare_matrix = multinomial_downsample(matrix, MIN_DEPTH, np.random.default_rng(rare_seed))
    rare_matrix = rare_matrix[:, score_genes].tocsr()
    rare_matrix.data = (rare_matrix.data / MIN_DEPTH).astype(np.float32)
    rare_matrix.eliminate_zeros()
    order = np.argsort(donor_idx, kind="stable")
    sorted_donors = donor_idx[order]
    boundaries = np.r_[0, np.cumsum(counts_by_donor[eligible_donors])]
    rows_by_donor = {}
    for pos, di in enumerate(eligible_donors):
        rows_by_donor[int(di)] = order[boundaries[pos]:boundaries[pos + 1]]
    global _WORK_FULL, _WORK_RARE, _WORK_LEFT, _WORK_RIGHT, _WORK_ROWS, _WORK_DATASET_CODE, _WORK_CT_CODE
    _WORK_FULL = score_matrix
    _WORK_RARE = rare_matrix
    _WORK_LEFT = work_left
    _WORK_RIGHT = work_right
    _WORK_ROWS = rows_by_donor
    _WORK_DATASET_CODE = 1 if dataset == "gao" else 2
    _WORK_CT_CODE = ct_idx
    tasks = [(int(di), str(meta["donor_ids"][di])) for di in eligible_donors]
    log(f"{dataset} {cell_type}: {len(tasks)} donors, {len(detected)} genes, {len(left)} pairs, {len(donor_idx)} cells")
    results = []
    if max_workers == 1:
        for task in tasks:
            results.append(_donor_gcl(task))
    else:
        context = mp.get_context("fork")
        with ProcessPoolExecutor(max_workers=max_workers, mp_context=context, initializer=_worker_init) as executor:
            for result in executor.map(_donor_gcl, tasks):
                results.append(result)
    donors_out = pd.DataFrame(results).sort_values("donor_id").reset_index(drop=True)
    split_cols = [f"split_abs_seed_{i}" for i in range(len(SPLIT_SEEDS))]
    info["split_half_noise_floor"] = float(np.nanmean(donors_out[split_cols].to_numpy(dtype=float)))
    info["rarefy_seed"] = int(rare_seed)
    info["status"] = "EVALUATED"
    info["elapsed_seconds"] = time.perf_counter() - started
    checkpoint.parent.mkdir(parents=True, exist_ok=True)
    pd.to_pickle({"version": ANALYSIS_VERSION, "info": info, "donors": donors_out}, checkpoint)
    _WORK_FULL = None
    _WORK_RARE = None
    _WORK_ROWS = None
    gc.collect()
    log(f"Finished {dataset} {cell_type} in {info['elapsed_seconds'] / 60:.1f} minutes; noise floor {info['split_half_noise_floor']:.6g}")
    return info, donors_out


def design_matrix(frame: pd.DataFrame, predictors: list[str], batch: bool) -> tuple[np.ndarray, list[str], list[str]]:
    columns = [np.ones(len(frame), dtype=float)]
    names = ["intercept"]
    omitted = []
    for predictor in predictors:
        values = frame[predictor].to_numpy(dtype=float)
        if not np.isfinite(values).all() or np.unique(values).size < 2:
            omitted.append(predictor)
            continue
        candidate = np.column_stack([*columns, values])
        if np.linalg.matrix_rank(candidate) == candidate.shape[1]:
            columns.append(values)
            names.append(predictor)
        else:
            omitted.append(predictor)
    if batch:
        batch_values = frame["study_or_pool"].astype(str)
        levels = sorted(batch_values.unique())
        for level in levels[1:]:
            values = (batch_values == level).to_numpy(dtype=float)
            candidate = np.column_stack([*columns, values])
            if np.linalg.matrix_rank(candidate) == candidate.shape[1]:
                columns.append(values)
                names.append(f"batch_{level}")
            else:
                omitted.append(f"batch_{level}")
    return np.column_stack(columns), names, omitted


def fit_effect(frame: pd.DataFrame, outcome: str, predictor: str, covariates: list[str], batch: bool) -> dict[str, object]:
    work = frame.dropna(subset=[outcome, predictor]).copy()
    base = {"status": "UNTESTABLE", "n_donors": int(len(work))}
    if len(work) < MIN_DONORS or work[predictor].nunique() < 2:
        return base
    predictors = [predictor, *covariates]
    x, names, omitted = design_matrix(work, predictors, batch)
    if predictor not in names or len(work) - x.shape[1] < 1:
        return {**base, "covariates_omitted": ";".join(omitted)}
    y = work[outcome].to_numpy(dtype=float)
    fit = sm.OLS(y, x).fit()
    idx = names.index(predictor)
    estimate = float(fit.params[idx])
    se = float(fit.bse[idx])
    df_resid = int(fit.df_resid)
    t95 = float(stats.t.ppf(0.975, df_resid))
    t80 = float(stats.t.ppf(0.80, df_resid))
    return {
        "status": "EVALUATED",
        "n_donors": int(len(work)),
        "effect": estimate,
        "se": se,
        "ci95_low": estimate - t95 * se,
        "ci95_high": estimate + t95 * se,
        "pvalue": float(fit.pvalues[idx]),
        "mde80": (t95 + t80) * se,
        "residual_df": df_resid,
        "covariates": ";".join(names),
        "covariates_omitted": ";".join(omitted),
    }


def group_comparison(frame: pd.DataFrame, seed: int) -> dict[str, object]:
    work = frame[frame["bin"].isin(["MATURE", "OLD"])].dropna(subset=["gcl_full"]).copy()
    mature = work[work["bin"] == "MATURE"]
    old = work[work["bin"] == "OLD"]
    base = {
        "status": "UNTESTABLE",
        "n_donors": int(len(work)),
        "n_donors_mature": int(len(mature)),
        "n_donors_old": int(len(old)),
    }
    if len(mature) < MIN_DONORS or len(old) < MIN_DONORS:
        return base
    x_batch, batch_names, omitted = design_matrix(work, [], True)
    if len(work) - x_batch.shape[1] < 2:
        return base
    y = work["gcl_full"].to_numpy(dtype=float)
    residual = y - x_batch @ np.linalg.lstsq(x_batch, y, rcond=None)[0]
    work["gcl_residual"] = residual + y.mean()
    work["old_indicator"] = (work["bin"] == "OLD").astype(int)
    fit = fit_effect(work, "gcl_residual", "old_indicator", [], False)
    if fit["status"] != "EVALUATED":
        return {**base, "covariates_omitted": ";".join(omitted)}
    is_old = work["old_indicator"].to_numpy(dtype=bool)
    rng = np.random.default_rng(seed)
    observed = float(residual[is_old].mean() - residual[~is_old].mean())
    null = np.empty(N_PERM, dtype=float)
    for perm in range(N_PERM):
        permuted = rng.permutation(is_old)
        null[perm] = residual[permuted].mean() - residual[~permuted].mean()
    pvalue = float((1 + np.sum(np.abs(null) >= abs(observed) - 1e-15)) / (N_PERM + 1))
    return {
        **fit,
        "effect": observed,
        "pvalue": pvalue,
        "mean_mature": float(work.loc[~is_old, "gcl_residual"].mean()),
        "mean_old": float(work.loc[is_old, "gcl_residual"].mean()),
        "n_donors_mature": int(len(mature)),
        "n_donors_old": int(len(old)),
        "n_permutations": N_PERM,
        "permutation_seed": int(seed),
        "batch_covariates": ";".join(batch_names),
        "covariates_omitted": ";".join(omitted),
    }


def empty_model_row(info: dict[str, object], model: str, n_donors: int) -> dict[str, object]:
    return {
        "dataset": info["dataset"],
        "cell_type": info["cell_type"],
        "model": model,
        "status": "UNTESTABLE",
        "effect_type": "old_minus_mature" if model == "group_comparison" else "slope_per_decade",
        "effect": np.nan,
        "slope_per_decade": np.nan,
        "group_difference_old_minus_mature": np.nan,
        "se": np.nan,
        "ci95_low": np.nan,
        "ci95_high": np.nan,
        "pvalue": np.nan,
        "q_bh": np.nan,
        "n_donors": int(n_donors),
        "n_genes": info["n_genes"],
        "n_pairs": info["n_pairs"],
        "pair_seed": info["pair_seed"],
        "split_half_noise_floor": info["split_half_noise_floor"],
        "mde80": np.nan,
    }


def relabel_mature_bins(frame: pd.DataFrame) -> pd.DataFrame:
    """Relabel donor-table bins to the mature scheme and drop under-30 rows.

    The donor table carries age_per_decade and the original YOUNG/MID/OLD
    bins. Donors 30-39 become MATURE; donors under 30 are dropped (they
    cannot appear after the mature-only inner merge, but the drop is kept
    as a guard so no young donor can leak into a group comparison).
    """
    out = frame.copy()
    ages = out["age_per_decade"].astype(float).to_numpy() * 10.0
    out["bin"] = np.where(
        ages < 30.0,
        "EXCLUDED",
        np.where((ages >= 30.0) & (ages < 40.0), "MATURE", out["bin"].astype(str).to_numpy()),
    )
    return out[out["bin"] != "EXCLUDED"].copy()


def model_rows(info: dict[str, object], donors: pd.DataFrame, donor_table: pd.DataFrame) -> list[dict[str, object]]:
    cell_type = str(info["cell_type"])
    metadata = donor_table[donor_table["group"] == cell_type].copy()
    metadata["donor_id"] = metadata["donor_id"].astype(str)
    merged = metadata.merge(donors, on="donor_id", how="inner", validate="one_to_one")
    merged = relabel_mature_bins(merged)
    eligible_metadata = metadata[metadata["n_cells_ge1500"] >= MIN_CELLS].copy()
    eligible_metadata = relabel_mature_bins(eligible_metadata)
    n_eligible = len(merged) if not donors.empty else int(info["n_donors_ge20"])
    rows = []
    specifications = (
        ("primary", "gcl_full", ["sex_male"], True),
        ("depth_adjusted", "gcl_full", ["sex_male", "mean_log_umi_ge1500"], True),
        ("rarefied_1500", "gcl_rarefied_1500", ["sex_male"], True),
        ("unadjusted", "gcl_full", ["sex_male"], False),
    )
    for model, outcome, covariates, batch in specifications:
        fit = fit_effect(merged, outcome, "age_per_decade", covariates, batch)
        row = empty_model_row(info, model, n_eligible)
        row.update(fit)
        if donors.empty:
            row["n_donors"] = n_eligible
        if fit["status"] == "EVALUATED":
            row["slope_per_decade"] = fit["effect"]
        rows.append(row)
    group_seed = stable_seed(PERM_SEED, info["dataset"], cell_type)
    fit = group_comparison(merged, group_seed)
    group_metadata = eligible_metadata[eligible_metadata["bin"].isin(["MATURE", "OLD"])]
    row = empty_model_row(info, "group_comparison", len(group_metadata))
    row.update(fit)
    row["n_donors"] = len(group_metadata)
    row["n_donors_mature"] = int((group_metadata["bin"] == "MATURE").sum())
    row["n_donors_old"] = int((group_metadata["bin"] == "OLD").sum())
    if fit["status"] == "EVALUATED":
        row["group_difference_old_minus_mature"] = fit["effect"]
    rows.append(row)
    return rows


def apply_bh(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out["q_bh"] = np.nan
    out["family_size"] = np.nan
    for (_, _), indices in out[out["status"] == "EVALUATED"].groupby(["dataset", "model"]).groups.items():
        pvalues = out.loc[indices, "pvalue"].to_numpy(dtype=float)
        valid = np.isfinite(pvalues)
        if not np.any(valid):
            continue
        idx = np.asarray(list(indices))[valid]
        out.loc[idx, "q_bh"] = multipletests(pvalues[valid], method="fdr_bh")[1]
        out.loc[idx, "family_size"] = int(valid.sum())
    return out


def validate_results(results: pd.DataFrame, infos: pd.DataFrame) -> None:
    models = {"primary", "group_comparison", "depth_adjusted", "rarefied_1500", "unadjusted"}
    assert len(results) == 5 * len(infos)
    for info in infos.to_dict("records"):
        subset = results[(results["dataset"] == info["dataset"]) & (results["cell_type"] == info["cell_type"])]
        assert len(subset) == 5 and set(subset["model"]) == models
        assert (subset["n_genes"] == info["n_genes"]).all()
        assert (subset["n_pairs"] == info["n_pairs"]).all()
    evaluated = results[results["status"] == "EVALUATED"]
    for column in ("effect", "se", "ci95_low", "ci95_high", "pvalue", "q_bh", "mde80", "split_half_noise_floor"):
        assert np.isfinite(evaluated[column].to_numpy(dtype=float)).all(), column
    assert (evaluated["n_donors"] >= MIN_DONORS).all()
    assert evaluated["pvalue"].between(0, 1).all()
    assert evaluated["q_bh"].between(0, 1).all()
    assert (evaluated["mde80"] >= 0).all()
    groups = evaluated[evaluated["model"] == "group_comparison"]
    assert (groups["n_donors_mature"] >= MIN_DONORS).all()
    assert (groups["n_donors_old"] >= MIN_DONORS).all()
    assert (groups["n_permutations"] == N_PERM).all()


def verdict_for_cell_type(primary: pd.Series, group: pd.Series) -> str:
    if primary["status"] != "EVALUATED":
        return ""
    if group["status"] != "EVALUATED":
        return "INCONCLUSIVE"
    slope = float(primary["effect"])
    difference = float(group["effect"])
    floor = float(primary["split_half_noise_floor"])
    if np.sign(slope) != np.sign(difference):
        return "INCONCLUSIVE"
    if float(primary["q_bh"]) < 0.05 and abs(slope) > floor:
        return "PROVES"
    if float(primary["q_bh"]) >= 0.05 and float(primary["mde80"]) <= floor:
        return "DOES NOT PROVE"
    return "INCONCLUSIVE"


def dataset_verdict(frame: pd.DataFrame) -> str:
    evaluated = frame[frame["verdict"] != "UNTESTABLE"]
    if evaluated.empty:
        return "INCONCLUSIVE"
    if (evaluated["verdict"] == "PROVES").any():
        failed_significant = evaluated[
            (evaluated["verdict"] == "INCONCLUSIVE") & (evaluated["primary_q"] < 0.05)
        ]
        return "PROVES" if failed_significant.empty else "INCONCLUSIVE"
    if (evaluated["verdict"] == "DOES NOT PROVE").all():
        return "DOES NOT PROVE"
    return "INCONCLUSIVE"


def fmt(value: object) -> str:
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return "NA"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.4g}"
    return str(value).replace("|", "/")


def markdown_table(frame: pd.DataFrame, columns: list[tuple[str, str]]) -> str:
    lines = ["| " + " | ".join(label for _, label in columns) + " |"]
    lines.append("| " + " | ".join("---" for _ in columns) + " |")
    for _, row in frame.iterrows():
        lines.append("| " + " | ".join(fmt(row[key]) for key, _ in columns) + " |")
    return "\n".join(lines)


def write_report(results: pd.DataFrame, infos: pd.DataFrame) -> None:
    spec = REPORT_MD.read_text(encoding="utf-8").split("## Results", 1)[0]
    spec = spec.replace("Status: analysis pending.", "Status: complete.")
    primary = results[results["model"] == "primary"].set_index(["dataset", "cell_type"])
    group = results[results["model"] == "group_comparison"].set_index(["dataset", "cell_type"])
    summary_rows = []
    for info in infos.to_dict("records"):
        key = (info["dataset"], info["cell_type"])
        p = primary.loc[key]
        g = group.loc[key]
        verdict = verdict_for_cell_type(p, g)
        direction = "coordination rises" if np.isfinite(p["effect"]) and p["effect"] > 0 else "coordination falls"
        summary_rows.append({
            **info,
            "n_donors": p["n_donors"],
            "slope": p["effect"],
            "primary_q": p["q_bh"],
            "mde80": p["mde80"],
            "group_diff": g["effect"],
            "group_q": g["q_bh"],
            "verdict": verdict if verdict else "UNTESTABLE",
            "direction": direction if verdict else "UNTESTABLE",
        })
    summary = pd.DataFrame(summary_rows)
    lines = [
        "## Results",
        "",
        "All effects are changes in mean pairwise Spearman correlation. Primary slopes are per decade. Group effects are OLD minus MATURE after batch residualization. NA denotes an UNTESTABLE or not applicable field.",
        "",
        markdown_table(summary, [
            ("dataset", "Dataset"), ("cell_type", "Cell type"), ("n_genes", "N genes"),
            ("n_pairs", "N pairs"), ("n_donors", "N donors"), ("slope", "Primary slope"),
            ("primary_q", "Primary q"), ("split_half_noise_floor", "Noise floor"),
            ("mde80", "MDE80"), ("group_diff", "OLD-MATURE"), ("group_q", "Group q"),
            ("verdict", "Verdict"), ("direction", "Direction"),
        ]),
        "",
        "### Sensitivity models",
        "",
    ]
    sensitivity_rows = []
    for info in infos.to_dict("records"):
        key = (info["dataset"], info["cell_type"])
        subset = results[(results["dataset"] == key[0]) & (results["cell_type"] == key[1])].set_index("model")
        sensitivity_rows.append({
            "dataset": key[0], "cell_type": key[1],
            "depth_effect": subset.loc["depth_adjusted", "effect"],
            "depth_q": subset.loc["depth_adjusted", "q_bh"],
            "unadjusted_effect": subset.loc["unadjusted", "effect"],
            "unadjusted_q": subset.loc["unadjusted", "q_bh"],
            "rare_effect": subset.loc["rarefied_1500", "effect"],
            "rare_q": subset.loc["rarefied_1500", "q_bh"],
        })
    sensitivity = pd.DataFrame(sensitivity_rows)
    lines.extend([
        markdown_table(sensitivity, [
            ("dataset", "Dataset"), ("cell_type", "Cell type"),
            ("depth_effect", "Depth adjusted slope"), ("depth_q", "Depth q"),
            ("unadjusted_effect", "Unadjusted slope"), ("unadjusted_q", "Unadjusted q"),
            ("rare_effect", "1500 rarefied slope"), ("rare_q", "1500 q"),
        ]),
        "",
        "### Verdicts",
        "",
    ])
    dataset_verdicts = {}
    for dataset in ("gao", "yazar"):
        subset = summary[summary["dataset"] == dataset]
        verdict = dataset_verdict(subset)
        dataset_verdicts[dataset] = verdict
        proved = subset[subset["verdict"] == "PROVES"]
        nulls = subset[subset["verdict"] == "DOES NOT PROVE"]
        inconclusive = subset[subset["verdict"] == "INCONCLUSIVE"]
        untestable = subset[subset["verdict"] == "UNTESTABLE"]
        direction = "coordination rises" if (proved["slope"] > 0).all() and not proved.empty else "coordination falls" if not proved.empty else "no proven direction"
        lines.append(f"- {dataset.capitalize()}: {verdict}; {direction}. PROVES {len(proved)}, DOES NOT PROVE {len(nulls)}, INCONCLUSIVE {len(inconclusive)}, UNTESTABLE {len(untestable)} cell types.")
        if not proved.empty:
            lines.append(f"  Proven units: {', '.join(proved['cell_type'])}.")
        if not inconclusive.empty:
            lines.append(f"  Inconclusive units: {', '.join(inconclusive['cell_type'])}.")
        if not untestable.empty:
            lines.append(f"  Untestable units: {', '.join(untestable['cell_type'])}.")
    overall = "PROVES" if all(v == "PROVES" for v in dataset_verdicts.values()) else "DOES NOT PROVE" if all(v == "DOES NOT PROVE" for v in dataset_verdicts.values()) else "INCONCLUSIVE"
    family_parts = []
    for dataset in ("gao", "yazar"):
        n_primary = len(results[(results["dataset"] == dataset) & (results["model"] == "primary") & (results["status"] == "EVALUATED")])
        n_group = len(results[(results["dataset"] == dataset) & (results["model"] == "group_comparison") & (results["status"] == "EVALUATED")])
        family_parts.append(f"{dataset} primary {n_primary}, group {n_group}")
    family_summary = ", ".join(family_parts)
    lines.extend([
        f"- Commit-label summary: {overall}. This does not pool donor data or replace the separate Gao and Yazar verdicts.",
        "- A rise would differ from the Levy 2020 prediction of falling coordination. No significant rise is described as confirmation of that prediction.",
        "",
        "### Reproducibility",
        "",
        f"Pair base seed {PAIR_SEED}, permutation base seed {PERM_SEED}, rarefaction base seed {RAREFY_SEED}, and split-half seeds {', '.join(map(str, SPLIT_SEEDS))}. Exact derived per-cell-type pair and permutation seeds are in the result CSV and analysis checkpoints.",
        f"Primary and group comparison family sizes: {family_summary}.",
        "No primary cells were downsampled. The only downsampling is the explicitly labeled 1500 rarefaction sensitivity.",
        "",
    ])
    write_clean_markdown(spec + "\n".join(lines), REPORT_MD)


def self_test() -> None:
    rng = np.random.default_rng(123)
    dense = rng.normal(size=(30, 4))
    left, right = fixed_pairs(4, 1)
    observed = gcl_from_matrix(sp.csr_matrix(dense), left, right)
    expected = np.mean([stats.spearmanr(dense[:, i], dense[:, j]).statistic for i, j in zip(left, right)])
    assert np.isclose(observed, expected, atol=1e-6), (observed, expected)
    sampled_left, sampled_right = fixed_pairs(2001, 12345)
    assert len(sampled_left) == MAX_SAMPLED_PAIRS
    assert len(np.unique(sampled_left * 2001 + sampled_right)) == MAX_SAMPLED_PAIRS
    a1, b1 = fixed_pairs(2001, 12345)
    assert np.array_equal(sampled_left, a1) and np.array_equal(sampled_right, b1)
    log("Self test passes against scipy.stats.spearmanr")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", choices=("gao", "yazar", "both"), default="both")
    parser.add_argument("--max-workers", type=int, default=2)
    parser.add_argument("--no-reuse-checkpoints", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if not REPORT_MD.exists() or "## Exact analysis specification" not in REPORT_MD.read_text(encoding="utf-8"):
        raise RuntimeError("The pre-analysis specification report must exist")
    wait_for_memory("Block 2 input loading")
    started_all = time.time()
    runtime_rows = []
    infos = []
    donors_by_key = {}
    datasets = ("gao", "yazar") if args.dataset == "both" else (args.dataset,)
    for dataset in datasets:
        started = time.time()
        if dataset == "gao":
            wait_for_memory("Gao h5ad metadata and count streaming")
            meta = reconcile_gao_metadata(DATA_ROOT / "organ17_gao_heart" / "GSE290367_integrated_RNA_adata.h5ad")
        else:
            wait_for_memory("Yazar h5ad metadata and count streaming")
            meta = reconcile_yazar_metadata(DATA_ROOT / "hcata" / "yazar_1e44db10.h5ad", DATA_ROOT / "hcata" / "yazar_cache")
        projection_info = build_human_gene_projection(meta["gene_symbols"])
        cache_path = build_count_cache(meta, projection_info["proj"], projection_info["audit"])
        for ct_idx in range(len(meta["cell_types"])):
            info, donors = process_cell_type(meta, cache_path, ct_idx, args.max_workers, not args.no_reuse_checkpoints)
            infos.append(info)
            donors_by_key[(dataset, str(meta["cell_types"][ct_idx]))] = donors
        del meta, projection_info
        gc.collect()
        runtime_rows.append({
            "stage": dataset,
            "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(started)),
            "finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(time.time())),
            "wall_minutes": (time.time() - started) / 60.0,
            "model": "direct-runner-no-delegated-model",
        })
    donor_tables = {
        dataset: pd.read_csv(REPO_ROOT / "results" / f"followup-between-subject-{dataset}-donors.csv")
        for dataset in datasets
    }
    for dataset in datasets:
        donor_tables[dataset] = donor_tables[dataset][~donor_tables[dataset]["is_pooled"]].copy()
    rows = []
    for info in infos:
        key = (str(info["dataset"]), str(info["cell_type"]))
        rows.extend(model_rows(info, donors_by_key[key], donor_tables[key[0]]))
    results = apply_bh(pd.DataFrame(rows))
    validate_results(results, pd.DataFrame(infos))
    results["age_window"] = "30-plus"
    results["mature_baseline"] = "30-39"
    runtime_rows.append({
        "stage": "all",
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(started_all)),
        "finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(time.time())),
        "wall_minutes": (time.time() - started_all) / 60.0,
        "model": "direct-runner-no-delegated-model",
    })
    write_clean_csv(pd.DataFrame(runtime_rows), RUNTIME_CSV)
    preferred = [
        "dataset", "cell_type", "model", "status", "effect_type", "effect",
        "slope_per_decade", "group_difference_old_minus_mature", "se", "ci95_low", "ci95_high",
        "pvalue", "q_bh", "family_size", "n_donors", "n_donors_mature", "n_donors_old",
        "n_genes", "n_pairs", "pair_seed", "split_half_noise_floor", "mde80", "residual_df",
        "n_permutations", "permutation_seed", "mean_mature", "mean_old", "covariates",
        "covariates_omitted", "batch_covariates",
    ]
    for column in preferred:
        if column not in results:
            results[column] = np.nan
    remaining = [column for column in results.columns if column not in preferred]
    results = results[[*preferred, *remaining]].sort_values(["dataset", "cell_type", "model"]).reset_index(drop=True)
    write_clean_csv(results, RESULT_CSV)
    write_report(results, pd.DataFrame(infos))
    log(f"Wrote {RESULT_CSV} and {REPORT_MD}")


if __name__ == "__main__":
    main()
