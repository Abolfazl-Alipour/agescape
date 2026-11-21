#!/usr/bin/env python3
"""Recompute November 17 donor distance variance with a 30-39 reference."""

from datetime import datetime, timezone
import os
from pathlib import Path
import resource
import sys
import time

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[name] = "1"

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import h5py
import numpy as np
import pandas as pd
from scipy import sparse as sp
from scipy import stats

from scripts.run_cross_species_expression_vector import reference_without_donor
from scripts.run_h3_expression_vector import hellinger_distance
from scripts.run_followup_between_subject_gao_yazar import (
    GAO_H5AD,
    YAZAR_CACHE,
    YAZAR_H5AD,
    build_human_gene_projection,
    read_csr_chunk,
    reconcile_gao_metadata,
    reconcile_yazar_metadata,
)
from src.followup_muscle import count_gate

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
OUTPUT = RESULTS / "followup-mature-nov17-donor-spread.csv"
DESCRIPTIVES = RESULTS / "followup-mature-nov17-donor-spread-descriptives.csv"
BY_CELL_TYPE = RESULTS / "followup-mature-nov17-donor-spread-by-cell-type.csv"
DISTANCES = RESULTS / "followup-mature-nov17-donor-spread-distances.csv"
RUNTIME = RESULTS / "followup-mature-nov17-runtime.csv"
SEED = 7
N_BOOT = 2000
CHUNK_ROWS = 8000
MIN_CELLS_PER_DONOR = 20
MIN_REFERENCE_DONORS = 2
MIN_LOO_REF_CELLS = 50
NEW_DISTANCE = "mature_30_39_ref_full_1500"
RETAINED_DISTANCE = "retained_mid_40_59_ref_full_1500"
rng = np.random.default_rng(SEED)


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def append_runtime(pipeline, start_utc, started, status, details=""):
    row = pd.DataFrame([{
        "pipeline": pipeline,
        "start_utc": start_utc,
        "end_utc": utc_now(),
        "duration_seconds": time.monotonic() - started,
        "model_used": "kimi-code/k3.1-pro",
        "status": status,
        "details": details,
    }])
    if RUNTIME.exists():
        row = pd.concat([pd.read_csv(RUNTIME), row], ignore_index=True)
    row.to_csv(RUNTIME, index=False)


def describe(series):
    values = series.dropna().to_numpy(dtype=float)
    if len(values) == 0:
        return None
    q1, median, q3 = np.quantile(values, [0.25, 0.50, 0.75])
    return {
        "n": int(len(values)),
        "mean": float(values.mean()),
        "sd": float(values.std(ddof=1)) if len(values) > 1 else np.nan,
        "q1": float(q1),
        "median": float(median),
        "q3": float(q3),
        "iqr": float(q3 - q1),
    }


def compare(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    stat, pvalue = stats.levene(a, b, center="median")
    var_a = float(np.var(a, ddof=1))
    var_b = float(np.var(b, ddof=1))
    boots = np.empty(N_BOOT)
    for index in range(N_BOOT):
        resampled_a = rng.choice(a, size=len(a), replace=True)
        resampled_b = rng.choice(b, size=len(b), replace=True)
        boots[index] = np.var(resampled_a, ddof=1) / np.var(resampled_b, ddof=1)
    low, high = np.percentile(boots, [2.5, 97.5])
    return {
        "levene_stat": float(stat),
        "levene_p": float(pvalue),
        "var_a": var_a,
        "var_b": var_b,
        "var_ratio": var_a / var_b,
        "ci_low": float(low),
        "ci_high": float(high),
    }


def pool_residualize(frame, variable, output):
    work = frame.loc[frame[variable].notna()].copy()
    pool_sizes = work.groupby("study_or_pool")["donor_id"].transform("size")
    singleton_donors = int((pool_sizes == 1).sum())
    work = work.loc[pool_sizes > 1].copy()
    if work.empty:
        raise ValueError(f"No non-singleton pools for {variable}")
    grand_mean = work[variable].mean()
    pool_means = work.groupby("study_or_pool")[variable].transform("mean")
    work[output] = work[variable] - pool_means + grand_mean
    return work, singleton_donors, int(work["study_or_pool"].nunique())


def add_comparisons(frame, variable, distance, dataset, rows, descriptions, adjustment="raw"):
    for bin_name in ["MATURE_30_39", "MID_40_59", "OLD_60_PLUS"]:
        summary = describe(frame.loc[frame["mature_bin"] == bin_name, variable])
        if summary is not None:
            descriptions.append({
                "dataset": dataset,
                "distance": distance,
                "adjustment": adjustment,
                "variable": variable,
                "bin": bin_name,
                **summary,
            })
    comparisons = [
        ("OLD_vs_MATURE_30_39", "OLD_60_PLUS", "MATURE_30_39"),
        ("MID_vs_MATURE_30_39", "MID_40_59", "MATURE_30_39"),
        ("OLD_vs_MID", "OLD_60_PLUS", "MID_40_59"),
    ]
    if distance == RETAINED_DISTANCE:
        comparisons = [comparisons[-1]]
    for comparison, group_a, group_b in comparisons:
        a = frame.loc[frame["mature_bin"] == group_a, variable].dropna().to_numpy(dtype=float)
        b = frame.loc[frame["mature_bin"] == group_b, variable].dropna().to_numpy(dtype=float)
        if len(a) < 2 or len(b) < 2:
            rows.append({
                "dataset": dataset,
                "distance": distance,
                "adjustment": adjustment,
                "variable": variable,
                "comparison": comparison,
                "group_a": group_a,
                "group_b": group_b,
                "n_a": len(a),
                "n_b": len(b),
                "status": "UNTESTABLE",
                "reason": "Fewer than two donors with a distance in at least one group",
            })
            continue
        rows.append({
            "dataset": dataset,
            "distance": distance,
            "adjustment": adjustment,
            "variable": variable,
            "comparison": comparison,
            "group_a": group_a,
            "group_b": group_b,
            "n_a": len(a),
            "n_b": len(b),
            "sd_a": float(np.std(a, ddof=1)),
            "sd_b": float(np.std(b, ddof=1)),
            "status": "TESTED",
            "reason": "",
            "method": "Levene median-centered (Brown-Forsythe); ratio bootstrap percentile CI, 2000 resamples, seed 7",
            **compare(a, b),
        })


def build_reference_stats(meta, counts, ref_counts, ref_donor_indices, n_groups, n_genes):
    donor_ids = meta["donor_ids"]
    donor_ages = meta["donors_df"]["age_years"].to_numpy(dtype=float)
    ref_row_by_donor = {donor: index for index, donor in enumerate(ref_donor_indices)}
    references = {}
    audit = []
    for group_index in range(n_groups):
        eligible = [
            donor for donor in ref_donor_indices
            if counts[donor, group_index] >= MIN_CELLS_PER_DONOR
        ]
        if len(eligible) < MIN_REFERENCE_DONORS:
            audit.append({"group_index": group_index, "n_reference_donors": len(eligible), "eligible": False})
            continue
        per_donor = {}
        total_counts = np.zeros(n_genes, dtype=np.int64)
        total_cells = 0
        for donor in eligible:
            vector = ref_counts[ref_row_by_donor[donor], group_index].copy()
            n_cells = int(counts[donor, group_index])
            per_donor[donor_ids[donor]] = {"counts": vector, "n_cells": n_cells}
            total_counts += vector
            total_cells += n_cells
        minimum_loo = min(total_cells - per_donor[donor_ids[donor]]["n_cells"] for donor in eligible)
        eligible_reference = minimum_loo >= MIN_LOO_REF_CELLS and total_counts.sum() > 0
        audit.append({
            "group_index": group_index,
            "n_reference_donors": len(eligible),
            "min_loo_reference_cells": int(minimum_loo),
            "eligible": bool(eligible_reference),
        })
        if eligible_reference:
            references[group_index] = {
                "counts": total_counts,
                "n_cells": total_cells,
                "donors": per_donor,
                "n_reference_donors": len(eligible),
                "min_loo_reference_cells": int(minimum_loo),
                "reference_age_min": float(donor_ages[eligible].min()),
                "reference_age_max": float(donor_ages[eligible].max()),
            }
    return references, audit


def recompute_distances(meta):
    start_utc = utc_now()
    started = time.monotonic()
    dataset = meta["dataset"]
    try:
        donor_ids = meta["donor_ids"]
        donor_ages = meta["donors_df"]["age_years"].to_numpy(dtype=float)
        mature_donors = np.flatnonzero(donor_ages >= 30)
        ref_donor_indices = np.flatnonzero((donor_ages >= 30) & (donor_ages <= 39))
        mature_mask = np.isin(meta["cell_donor_idx"], mature_donors)
        groups = [*meta["cell_types"], meta["pooled_name"]]
        n_donors = len(donor_ids)
        n_groups = len(groups)
        n_raw_vars = len(meta["gene_symbols"])
        projection_data = build_human_gene_projection(meta["gene_symbols"])
        projection = projection_data["proj"]
        n_genes = projection.shape[1]
        ref_row_by_donor = {donor: index for index, donor in enumerate(ref_donor_indices)}
        ref_counts = np.zeros((len(ref_donor_indices), n_groups, n_genes), dtype=np.int64)
        counts = np.zeros((n_donors, n_groups), dtype=np.int64)
        with h5py.File(meta["h5_path"], "r") as handle:
            if meta["matrix_slot"] == "layers/counts":
                matrix_group = handle["layers"]["counts"]
            elif meta["matrix_slot"] == "raw/X":
                matrix_group = handle["raw"]["X"]
            else:
                matrix_group = handle["X"]
            for start in range(0, meta["n_obs"], CHUNK_ROWS):
                stop = min(start + CHUNK_ROWS, meta["n_obs"])
                local_keep = mature_mask[start:stop]
                if not local_keep.any():
                    continue
                chunk = read_csr_chunk(matrix_group, start, stop, n_raw_vars, meta["indptr"])[local_keep]
                gate = count_gate(chunk.data)
                if gate["nonfinite"] or gate["negative"] or gate["noninteger"]:
                    raise ValueError(f"Count gate failed: {gate}")
                clean = (chunk @ projection).tocsr()
                umi = np.rint(np.asarray(clean.sum(axis=1)).ravel()).astype(np.int64)
                keep = umi >= 1500
                if not keep.any():
                    continue
                clean = clean[keep]
                donor_chunk = meta["cell_donor_idx"][start:stop][local_keep][keep]
                ct_chunk = meta["cell_ct_idx"][start:stop][local_keep][keep]
                pairs = np.concatenate([
                    donor_chunk * n_groups + ct_chunk,
                    donor_chunk * n_groups + n_groups - 1,
                ])
                counts_flat = np.bincount(pairs, minlength=n_donors * n_groups)
                counts += counts_flat.reshape(n_donors, n_groups)
                ref_local = np.isin(donor_chunk, ref_donor_indices)
                if ref_local.any():
                    ref_d = np.concatenate([donor_chunk[ref_local], donor_chunk[ref_local]])
                    ref_g = np.concatenate([ct_chunk[ref_local], np.full(ref_local.sum(), n_groups - 1)])
                    ref_pairs = ref_d * n_groups + ref_g
                    ref_matrix = sp.vstack([clean[ref_local], clean[ref_local]], format="csr")
                    unique_pairs, inverse = np.unique(ref_pairs, return_inverse=True)
                    indicator = sp.csr_matrix(
                        (np.ones(len(inverse)), (inverse, np.arange(len(inverse)))),
                        shape=(len(unique_pairs), len(inverse)),
                    )
                    sums = (indicator @ ref_matrix).tocsr()
                    for index, pair in enumerate(unique_pairs):
                        donor = int(pair // n_groups)
                        group = int(pair % n_groups)
                        begin, end = int(sums.indptr[index]), int(sums.indptr[index + 1])
                        ref_counts[ref_row_by_donor[donor], group, sums.indices[begin:end]] += np.rint(sums.data[begin:end]).astype(np.int64)
                if start % (CHUNK_ROWS * 25) == 0:
                    print(f"{dataset} reference pass at row {stop}/{meta['n_obs']}", flush=True)
        references, reference_audit = build_reference_stats(
            meta, counts, ref_counts, ref_donor_indices, n_groups, n_genes
        )
        del ref_counts
        eligible = counts >= MIN_CELLS_PER_DONOR
        for group_index in range(n_groups):
            if group_index not in references:
                eligible[:, group_index] = False
        distance_sums = np.zeros((n_donors, n_groups), dtype=np.float64)
        distance_counts = np.zeros((n_donors, n_groups), dtype=np.int64)
        with h5py.File(meta["h5_path"], "r") as handle:
            if meta["matrix_slot"] == "layers/counts":
                matrix_group = handle["layers"]["counts"]
            elif meta["matrix_slot"] == "raw/X":
                matrix_group = handle["raw"]["X"]
            else:
                matrix_group = handle["X"]
            for start in range(0, meta["n_obs"], CHUNK_ROWS):
                stop = min(start + CHUNK_ROWS, meta["n_obs"])
                local_keep = mature_mask[start:stop]
                if not local_keep.any():
                    continue
                chunk = read_csr_chunk(matrix_group, start, stop, n_raw_vars, meta["indptr"])[local_keep]
                clean = (chunk @ projection).tocsr()
                umi = np.rint(np.asarray(clean.sum(axis=1)).ravel()).astype(np.int64)
                keep = umi >= 1500
                if not keep.any():
                    continue
                clean = clean[keep]
                donor_chunk = meta["cell_donor_idx"][start:stop][local_keep][keep]
                ct_chunk = meta["cell_ct_idx"][start:stop][local_keep][keep]
                pairs = np.concatenate([
                    donor_chunk * n_groups + ct_chunk,
                    donor_chunk * n_groups + n_groups - 1,
                ])
                matrix = sp.vstack([clean, clean], format="csr")
                unique_pairs = np.unique(pairs)
                unique_pairs = unique_pairs[eligible.flat[unique_pairs]]
                for pair in unique_pairs:
                    donor = int(pair // n_groups)
                    group = int(pair % n_groups)
                    selected = pairs == pair
                    reference = reference_without_donor(references[group], donor_ids[donor])
                    if reference["n_cells"] < MIN_LOO_REF_CELLS or reference["counts"].sum() <= 0:
                        continue
                    distances = hellinger_distance(matrix[selected], reference["counts"])
                    distance_sums[donor, group] += float(distances.sum())
                    distance_counts[donor, group] += len(distances)
                if start % (CHUNK_ROWS * 25) == 0:
                    print(f"{dataset} scoring pass at row {stop}/{meta['n_obs']}", flush=True)
        output_rows = []
        for donor in mature_donors:
            for group in range(n_groups):
                reference = references.get(group)
                distance = np.nan
                if distance_counts[donor, group] > 0:
                    distance = distance_sums[donor, group] / distance_counts[donor, group]
                output_rows.append({
                    "dataset": "Gao heart" if dataset == "gao" else "Yazar blood",
                    "donor_id": donor_ids[donor],
                    "age_years": donor_ages[donor],
                    "mature_bin": "MATURE_30_39" if donor_ages[donor] <= 39 else ("MID_40_59" if donor_ages[donor] <= 59 else "OLD_60_PLUS"),
                    "group": groups[group],
                    "is_pooled": group == n_groups - 1,
                    "n_cells_ge1500": int(counts[donor, group]),
                    "distance_mature_30_39": distance,
                    "n_reference_donors": reference["n_reference_donors"] if reference else 0,
                    "min_loo_reference_cells": reference["min_loo_reference_cells"] if reference else 0,
                    "reference_age_min": reference["reference_age_min"] if reference else np.nan,
                    "reference_age_max": reference["reference_age_max"] if reference else np.nan,
                })
        result = pd.DataFrame(output_rows)
        assert result["age_years"].min() >= 30
        valid_refs = result.loc[result["distance_mature_30_39"].notna()]
        assert valid_refs["reference_age_min"].min() >= 30
        assert valid_refs["reference_age_max"].max() <= 39
        committed = pd.read_csv(RESULTS / f"followup-between-subject-{dataset}-donors.csv")
        committed = committed.loc[committed["age_years"] >= 30, ["donor_id", "group", "n_cells_ge1500"]]
        check = committed.merge(result, on=["donor_id", "group"], suffixes=("_committed", "_computed"), validate="one_to_one")
        if len(check) != len(committed) or not np.array_equal(check["n_cells_ge1500_committed"], check["n_cells_ge1500_computed"]):
            raise ValueError(f"{dataset}: recomputed mature cell counts do not match committed donor table")
        details = (
            f"mature donors {len(mature_donors)}; reference donors {len(ref_donor_indices)}; "
            f"retained genes {n_genes}; max RSS MB {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024:.1f}"
        )
        append_runtime(f"mature_nov17_distances_{dataset}", start_utc, started, "completed", details)
        print(details, flush=True)
        return result
    except Exception as exc:
        append_runtime(
            f"mature_nov17_distances_{dataset}", start_utc, started, "FAILED", f"{type(exc).__name__}: {exc}"
        )
        raise


def main():
    started = time.monotonic()
    start_utc = utc_now()
    try:
        distance_frames = [
            recompute_distances(reconcile_gao_metadata(GAO_H5AD)),
            recompute_distances(reconcile_yazar_metadata(YAZAR_H5AD, YAZAR_CACHE)),
        ]
        distances = pd.concat(distance_frames, ignore_index=True)
        distances.to_csv(DISTANCES, index=False)
        rows = []
        descriptions = []
        supplement = []
        for dataset, short_name, pooled_group in [
            ("Gao heart", "gao", "All_LV_Nuclei"),
            ("Yazar blood", "yazar", "All_Cells"),
        ]:
            original = pd.read_csv(RESULTS / f"followup-between-subject-{short_name}-donors.csv")
            original = original.loc[original["age_years"] >= 30].copy()
            original["mature_bin"] = np.select(
                [original["age_years"] <= 39, original["age_years"] <= 59],
                ["MATURE_30_39", "MID_40_59"], default="OLD_60_PLUS"
            )
            recomputed = distances.loc[distances["dataset"] == dataset].copy()
            pooled = original.loc[original["group"] == pooled_group].merge(
                recomputed.loc[recomputed["is_pooled"], ["donor_id", "distance_mature_30_39"]],
                on="donor_id", validate="one_to_one"
            )
            add_comparisons(pooled, "distance_mature_30_39", NEW_DISTANCE, dataset, rows, descriptions)
            add_comparisons(original.loc[original["group"] == pooled_group], "a1_dist_mature_ref_mid_loo", RETAINED_DISTANCE, dataset, rows, descriptions)
            if dataset == "Yazar blood":
                for variable, distance in [
                    ("distance_mature_30_39", NEW_DISTANCE),
                    ("a1_dist_mature_ref_mid_loo", RETAINED_DISTANCE),
                ]:
                    residual, n_singleton, n_pools = pool_residualize(pooled, variable, "pool_residualized")
                    note = f"Excluded {n_singleton} singleton-pool donors; {n_pools} pools retained"
                    before = len(rows)
                    add_comparisons(residual, "pool_residualized", distance, dataset, rows, descriptions, "pool_residualized")
                    for row in rows[before:]:
                        row["method"] = row.get("method", "") + "; " + note
            for group, group_frame in original.loc[original["group"] != pooled_group].groupby("group"):
                cell_distances = recomputed.loc[recomputed["group"] == group, ["donor_id", "distance_mature_30_39"]]
                group_frame = group_frame.merge(cell_distances, on="donor_id", how="left", validate="one_to_one")
                raw_a = group_frame.loc[group_frame["mature_bin"] == "OLD_60_PLUS", "distance_mature_30_39"].dropna().to_numpy(dtype=float)
                raw_b = group_frame.loc[group_frame["mature_bin"] == "MATURE_30_39", "distance_mature_30_39"].dropna().to_numpy(dtype=float)
                if len(raw_a) < 20 or len(raw_b) < 20:
                    continue
                versions = [("raw", group_frame, "distance_mature_30_39", 0)]
                if dataset == "Yazar blood":
                    residual, n_singleton, _ = pool_residualize(group_frame, "distance_mature_30_39", "pool_residualized")
                    versions.append(("pool_residualized", residual, "pool_residualized", n_singleton))
                for adjustment, version_frame, variable, n_singleton in versions:
                    a = version_frame.loc[version_frame["mature_bin"] == "OLD_60_PLUS", variable].dropna().to_numpy(dtype=float)
                    b = version_frame.loc[version_frame["mature_bin"] == "MATURE_30_39", variable].dropna().to_numpy(dtype=float)
                    if len(a) < 20 or len(b) < 20:
                        continue
                    supplement.append({
                        "dataset": dataset,
                        "cell_type_group": group,
                        "distance": NEW_DISTANCE,
                        "adjustment": adjustment,
                        "n_old": len(a),
                        "n_mature_30_39": len(b),
                        "sd_old": float(np.std(a, ddof=1)),
                        "sd_mature_30_39": float(np.std(b, ddof=1)),
                        "n_singleton_pool_donors_excluded": n_singleton,
                        **compare(a, b),
                    })
        pd.DataFrame(rows).to_csv(OUTPUT, index=False)
        pd.DataFrame(descriptions).to_csv(DESCRIPTIVES, index=False)
        pd.DataFrame(supplement).to_csv(BY_CELL_TYPE, index=False)
        append_runtime("mature_nov17_donor_spread", start_utc, started, "completed", f"comparisons {len(rows)}; supplement rows {len(supplement)}")
        print(pd.DataFrame(rows).to_string(index=False))
        print(pd.DataFrame(supplement).to_string(index=False))
    except Exception as exc:
        append_runtime("mature_nov17_donor_spread", start_utc, started, "FAILED", f"{type(exc).__name__}: {exc}")
        raise


if __name__ == "__main__":
    main()
