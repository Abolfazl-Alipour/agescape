#!/usr/bin/env python3
"""Validate and analyze Jin GSE202601 donor-level ovary endpoints."""

from __future__ import annotations

import argparse
import json
import resource
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from src.aging_programs import load_program_tiers  # noqa: E402
from src.followup_muscle import (  # noqa: E402
    NULL_ITERATIONS, NULL_SEED, add_bh, aggregate_block, assert_count_gate,
    clr_transform, hellinger_distances, matched_null_effect_pvalue,
    random_signature_scores, technical_gene_mask,
)
from src.followup_ovary import (  # noqa: E402
    JIN_MODELS, build_jin_design, fit_gene_models, fit_one, sha256,
    status_row, write_csv,
)

DATA = REPO.parent / f"{REPO.name}-data" / "ovary" / "GSE202601"
CACHE = DATA / "interim" / "jin_donor_aggregates.npz"
COHORTS = REPO / "data" / "manifest" / "ovary_cohorts.json"
METADATA = REPO.parent / f"{REPO.name}-data" / "ovary-metadata" / "GSE202601" / "GSE202601_human_ovary_snRNA-seq_metadata.txt.gz"
RESULTS = REPO / "results"
PREFIX = "followup-ovary-jin"
COHORT = "GSE202601"
MIN_CELLS = 20
MIN_REFERENCE_CELLS = 50
MIN_UMI = 1500


def manifest_donors():
    return json.loads(COHORTS.read_text())[COHORT]["donors"]


def donor_path(donor):
    return DATA / donor["rna_filtered_count_matrix"]["url"].rsplit("/", 1)[-1]


def load_metadata():
    meta = pd.read_csv(METADATA, sep="\t", index_col=0)
    if len(meta) != 42568 or not meta.index.is_unique:
        raise ValueError("Unexpected Jin snRNA metadata dimensions or duplicate barcodes")
    if set(meta.columns) != {"age", "group", "celltype1", "celltype2"}:
        raise ValueError("Unexpected Jin metadata columns")
    meta["celltype2"] = meta["celltype2"].astype(str).str.strip()
    if set(meta["celltype2"]) != {"BEC", "EpiC", "GC", "IC", "LEC", "SC", "SMC", "TC"}:
        raise ValueError("Released Jin broad annotations changed")
    for donor in manifest_donors():
        label = f"ovary{donor['snrna_metadata_age_years']}"
        part = meta[meta["age"].eq(label)]
        if len(part) != donor["snrna_metadata_nuclei"]:
            raise ValueError(f"Jin metadata nucleus mismatch for {donor['donor']}")
        expected_group = "young" if donor["group"] == "young" else "old"
        if set(part["group"]) != {expected_group}:
            raise ValueError(f"Jin group mismatch for {donor['donor']}")
    return meta


def read_h5(donor, donor_number, meta, validate_counts=True):
    path = donor_path(donor)
    expected_size = donor["rna_filtered_count_matrix"]["size_bytes"]
    if not path.is_file() or path.stat().st_size != expected_size:
        raise ValueError(f"Missing or unexpected Jin H5 size: {path}")
    label = f"ovary{donor['snrna_metadata_age_years']}"
    annotations = meta[meta["age"].eq(label)]
    with h5py.File(path, "r") as handle:
        group = handle["matrix"]
        shape = tuple(int(x) for x in group["shape"][:])
        data = group["data"]
        total_raw_umi = 0
        if validate_counts:
            for start in range(0, len(data), 5_000_000):
                values = data[start:start + 5_000_000]
                assert_count_gate(values)
                total_raw_umi += int(values.sum(dtype=np.int64))
        matrix = sp.csc_matrix((data[:], group["indices"][:], group["indptr"][:]), shape=shape)
        barcodes = pd.Index(group["barcodes"][:].astype(str))
        gene_ids = group["features"]["id"][:].astype(str)
        gene_names = group["features"]["name"][:].astype(str)
        feature_types = group["features"]["feature_type"][:].astype(str)
    if not barcodes.is_unique or not pd.Index(gene_ids).is_unique:
        raise ValueError(f"Duplicate Jin barcodes or feature IDs for {donor['donor']}")
    if set(feature_types) != {"Gene Expression"} or shape != (len(gene_ids), len(barcodes)):
        raise ValueError(f"Unexpected Jin feature schema for {donor['donor']}")
    if not barcodes.str.endswith("-1").all():
        raise ValueError(f"Unexpected Jin raw barcode suffix for {donor['donor']}")
    converted = pd.Index(barcodes.str[:-2] + f".{donor_number}")
    missing = annotations.index.difference(converted)
    if len(missing) or converted.duplicated().any():
        raise ValueError(f"Jin annotated barcode join failed for {donor['donor']}: {len(missing)} missing")
    selected = converted.get_indexer(annotations.index)
    annotated = matrix[:, selected].T.tocsr()
    depths = np.asarray(annotated.sum(axis=1)).ravel()
    if (depths <= 0).any():
        raise ValueError(f"Zero-UMI annotated Jin nucleus for {donor['donor']}")
    summary = {
        "donor": donor["donor"],
        "rna_gsm": donor["rna_gsm"],
        "snrna_metadata_age_label": label,
        "geo_exact_age_years": donor["geo_exact_age_years"],
        "snrna_metadata_age_years": donor["snrna_metadata_age_years"],
        "exact_age_status": donor["exact_age_status"],
        "group": donor["group"],
        "old": int(donor["group"] != "young"),
        "raw_h5_barcodes": len(barcodes),
        "annotated_nuclei": len(annotations),
        "unannotated_h5_barcodes": len(barcodes) - len(annotations),
        "stored_count_values": matrix.nnz,
        "total_raw_umi_all_h5_barcodes": total_raw_umi if validate_counts else np.nan,
        "total_raw_umi_annotated": int(depths.sum()),
        "mean_umi_annotated": float(depths.mean()),
        "median_umi_annotated": float(np.median(depths)),
        "mean_log_umi_annotated": float(np.log(depths).mean()),
        "all_counts_finite_nonnegative_integer_valued": validate_counts,
        "all_annotations_matched_one_to_one": True,
    }
    return summary, annotated, annotations, gene_ids, gene_names


def build_cache():
    started = time.time()
    meta = load_metadata()
    donors = manifest_donors()
    cell_types = sorted(meta["celltype2"].unique())
    units = pd.MultiIndex.from_product(
        [[d["donor"] for d in donors], cell_types], names=["donor", "cell_type"]
    ).to_frame(index=False)
    donor_info = {d["donor"]: d for d in donors}
    units["old"] = units["donor"].map(lambda d: int(donor_info[d]["group"] != "young"))
    units["group"] = units["donor"].map(lambda d: donor_info[d]["group"])
    unit_lookup = {(r.donor, r.cell_type): i for i, r in enumerate(units.itertuples())}
    donor_lookup = {d["donor"]: i for i, d in enumerate(donors)}
    raw_counts = log_sums = qc_counts = None
    n_units = len(units)
    n_cells = np.zeros(n_units, dtype=np.int64)
    qc_n = np.zeros(n_units, dtype=np.int64)
    sum_log_umi = np.zeros(n_units)
    sum_umi = np.zeros(n_units)
    qc_sum_log_umi = np.zeros(n_units)
    qc_sum_umi = np.zeros(n_units)
    composition_n = np.zeros(len(donors), dtype=np.int64)
    composition_log = np.zeros(len(donors))
    composition_umi = np.zeros(len(donors))
    summaries = []
    expected_ids = expected_names = None
    technical = None
    stored_values = 0
    for donor_number, donor in enumerate(donors, 1):
        summary, matrix, annotations, gene_ids, gene_names = read_h5(donor, donor_number, meta)
        if expected_ids is None:
            expected_ids, expected_names = gene_ids, gene_names
            technical = technical_gene_mask(expected_names)
            raw_counts = np.zeros((n_units, len(gene_ids)), dtype=np.int64)
            log_sums = np.zeros((n_units, len(gene_ids)), dtype=np.float32)
            qc_counts = np.zeros((n_units, len(gene_ids)), dtype=np.int64)
        elif not np.array_equal(gene_ids, expected_ids) or not np.array_equal(gene_names, expected_names):
            raise ValueError(f"Jin feature order differs for {donor['donor']}")
        codes = np.asarray([unit_lookup[(donor["donor"], cell)] for cell in annotations["celltype2"]])
        depths = np.asarray(matrix.sum(axis=1)).ravel()
        nontechnical = np.asarray(matrix[:, ~technical].sum(axis=1)).ravel()
        np.add.at(n_cells, codes, 1)
        np.add.at(sum_log_umi, codes, np.log(depths))
        np.add.at(sum_umi, codes, depths)
        donor_index = donor_lookup[donor["donor"]]
        composition_n[donor_index] = matrix.shape[0]
        composition_log[donor_index] = np.log(depths).sum()
        composition_umi[donor_index] = depths.sum()
        log_matrix = matrix.astype(np.float64).tocsr()
        log_matrix.data *= np.repeat(1e6 / depths, np.diff(log_matrix.indptr))
        np.log1p(log_matrix.data, out=log_matrix.data)
        aggregate_block(log_matrix, codes, log_sums)
        aggregate_block(matrix.astype(np.float64), codes, raw_counts)
        qc = (depths >= MIN_UMI) & (nontechnical > 0)
        if qc.any():
            aggregate_block(matrix[qc].astype(np.float64), codes[qc], qc_counts)
            np.add.at(qc_n, codes[qc], 1)
            np.add.at(qc_sum_log_umi, codes[qc], np.log(depths[qc]))
            np.add.at(qc_sum_umi, codes[qc], depths[qc])
        summaries.append(summary)
        stored_values += summary["stored_count_values"]
        print(f"Jin {donor['donor']}: {matrix.shape[0]} annotated nuclei; raw gate passed", flush=True)
    if (n_cells <= 0).any() or (composition_n <= 0).any():
        raise ValueError("Unexpected empty Jin donor-cell-type unit")
    cache = {
        "gene_ids": expected_ids,
        "gene_names": expected_names,
        "donor": units["donor"].astype(str).to_numpy(dtype=str),
        "cell_type": units["cell_type"].astype(str).to_numpy(dtype=str),
        "group": units["group"].astype(str).to_numpy(dtype=str),
        "old": units["old"].to_numpy(np.int8),
        "n_cells": n_cells,
        "mean_log_umi": sum_log_umi / n_cells,
        "mean_umi": sum_umi / n_cells,
        "raw_counts": raw_counts,
        "log_sums": log_sums,
        "qc_n": qc_n,
        "qc_mean_log_umi": np.divide(qc_sum_log_umi, qc_n, out=np.full(n_units, np.nan), where=qc_n > 0),
        "qc_mean_umi": np.divide(qc_sum_umi, qc_n, out=np.full(n_units, np.nan), where=qc_n > 0),
        "qc_counts": qc_counts,
        "composition_donor": np.asarray([d["donor"] for d in donors]),
        "composition_n_cells": composition_n,
        "composition_mean_log_umi": composition_log / composition_n,
        "composition_mean_umi": composition_umi / composition_n,
    }
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    np.savez(CACHE, **cache)
    donor_frame = pd.DataFrame(summaries)
    write_csv(donor_frame, RESULTS / f"{PREFIX}-donor-qc.csv")
    counts = pd.DataFrame({
        "donor": cache["donor"], "cell_type": cache["cell_type"],
        "group": cache["group"], "n_cells": n_cells, "qc_n_cells": qc_n,
        "mean_log_umi": cache["mean_log_umi"],
    })
    write_csv(counts, RESULTS / f"{PREFIX}-celltype-donor-counts.csv")
    eligibility_rows = []
    for cell_type, part in counts.groupby("cell_type"):
        eligible = part["n_cells"] >= MIN_CELLS
        n_young = int((eligible & part["old"].eq(0)).sum()) if "old" in part else int((eligible & part["group"].eq("young")).sum())
        n_old = int((eligible & part["group"].ne("young")).sum())
        passed = n_young == 4 and n_old == 4
        eligibility_rows.append({
            "cell_type": cell_type,
            "n_donors_observed": len(part),
            "n_donors_ge20": int(eligible.sum()),
            "n_young_ge20": n_young,
            "n_reproductively_aged_ge20": n_old,
            "min_cells_per_donor": int(part["n_cells"].min()),
            "eligible": passed,
            "reason": "pass" if passed else "failed_all_four_donors_per_group_20_nuclei_floor",
        })
    eligibility = pd.DataFrame(eligibility_rows).sort_values("cell_type").reset_index(drop=True)
    write_csv(eligibility, RESULTS / f"{PREFIX}-celltype-eligibility.csv")
    verification = {
        "cohort": COHORT,
        "donors": len(donors),
        "young_donors": 4,
        "reproductively_aged_donors": 4,
        "metadata_rows": len(meta),
        "annotated_nuclei": int(donor_frame["annotated_nuclei"].sum()),
        "raw_h5_barcodes": int(donor_frame["raw_h5_barcodes"].sum()),
        "unannotated_h5_barcodes": int(donor_frame["unannotated_h5_barcodes"].sum()),
        "stored_count_values_checked": int(stored_values),
        "all_counts_finite_nonnegative_integer_valued": True,
        "feature_order_identical_across_donors": True,
        "all_annotations_matched_one_to_one": True,
        "ovary5_numeric_age_used": False,
        "continuous_age_model": False,
        "cache_bytes": CACHE.stat().st_size,
        "elapsed_seconds": time.time() - started,
        "peak_rss_bytes": int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024),
        "memory_limit_gb": 8,
        "minimum_free_gb": 30,
        "outcomes_examined": False,
    }
    (RESULTS / f"{PREFIX}-verification.json").write_text(json.dumps(verification, indent=2, sort_keys=True) + "\n")
    print("Jin raw-count and annotation gates passed; no outcomes fitted in build-cache.", flush=True)


def load_cache():
    with np.load(CACHE, allow_pickle=False) as data:
        return {key: data[key] for key in data.files}


def unit_frame(cache, indices, qc=False):
    count_key = "qc_n" if qc else "n_cells"
    depth_key = "qc_mean_log_umi" if qc else "mean_log_umi"
    return pd.DataFrame({
        "donor": cache["donor"][indices].astype(str),
        "cell_type": cache["cell_type"][indices].astype(str),
        "group": cache["group"][indices].astype(str),
        "old": cache["old"][indices],
        "n_cells": cache[count_key][indices],
        "mean_log_umi": cache[depth_key][indices],
    })


def selected_units(cache, eligibility, cell_type, qc=False):
    row = eligibility[eligibility["cell_type"].eq(cell_type)]
    if row.empty or not bool(row["eligible"].iloc[0]):
        return np.array([], dtype=int)
    counts = cache["qc_n"] if qc else cache["n_cells"]
    indices = np.flatnonzero((cache["cell_type"].astype(str) == cell_type) & (counts >= MIN_CELLS))
    return indices[np.argsort(cache["donor"][indices].astype(str))]


def run_p1(cache, eligibility):
    program = load_program_tiers(REPO, "human")["P1_Clean"]
    gene_names = cache["gene_names"].astype(str)
    donor_rows, model_rows, null_rows, statuses = [], [], [], []
    for cell_type in eligibility["cell_type"]:
        indices = selected_units(cache, eligibility, cell_type)
        if not len(indices):
            for model in JIN_MODELS:
                statuses.append(status_row(COHORT, model, "P1_transport", cell_type, "untestable", "failed_base_donor_nuclei_floors"))
            continue
        frame = unit_frame(cache, indices)
        matrix = (cache["log_sums"][indices] / cache["n_cells"][indices, None]).astype(float).T
        detected = cache["raw_counts"][indices].sum(axis=0) > 0
        up = np.flatnonzero(np.isin(gene_names, program["up"]) & detected)
        down = np.flatnonzero(np.isin(gene_names, program["down"]) & detected)
        if not len(up) or not len(down):
            for model in JIN_MODELS:
                statuses.append(status_row(COHORT, model, "P1_transport", cell_type, "untestable", "no_detected_P1_up_or_down_program_genes"))
            continue
        scores = matrix[up].mean(axis=0) - matrix[down].mean(axis=0)
        donor_frame = frame.copy()
        donor_frame["P1_score"] = scores
        donor_frame["p1_up_detected"] = len(up)
        donor_frame["p1_down_detected"] = len(down)
        donor_rows.append(donor_frame)
        for model in JIN_MODELS:
            fitted, design = fit_one(frame, scores, model, COHORT)
            row = {"model": model, "endpoint": "P1_transport", "cell_type": cell_type,
                   "p1_up_detected": len(up), "p1_down_detected": len(down), **fitted}
            if model == "primary_group_depth_adjusted":
                try:
                    random_scores = random_signature_scores(matrix, up, down)
                    null_effects = []
                    for random_score in random_scores:
                        random_design, _ = build_jin_design(frame, model)
                        from src.followup_muscle import fit_ols
                        null_effects.append(fit_ols(random_score, random_design)["effect"][0])
                    null_effects = np.asarray(null_effects)
                    row.update({
                        "random_iterations": NULL_ITERATIONS,
                        "random_seed": NULL_SEED,
                        "random_mean_effect": float(null_effects.mean()),
                        "random_sd_effect": float(null_effects.std(ddof=1)),
                        "random_q025": float(np.quantile(null_effects, 0.025)),
                        "random_q975": float(np.quantile(null_effects, 0.975)),
                        "random_pvalue": matched_null_effect_pvalue(fitted["effect"], null_effects),
                    })
                    for iteration, effect in enumerate(null_effects, 1):
                        null_rows.append({"cell_type": cell_type, "iteration": iteration,
                                          "effect_reproductively_aged_minus_young": effect, "seed": NULL_SEED})
                except ValueError as exc:
                    row["random_null_status"] = f"untestable: {exc}"
            model_rows.append(row)
            statuses.append(status_row(COHORT, model, "P1_transport", cell_type, "completed", "pass"))
    models = add_bh(pd.DataFrame(model_rows), ["model"])
    models["random_q_bh"] = np.nan
    primary = models["model"].eq("primary_group_depth_adjusted") & models["random_pvalue"].notna() if "random_pvalue" in models else pd.Series(False, index=models.index)
    if primary.any():
        models.loc[primary, "random_q_bh"] = multipletests(models.loc[primary, "random_pvalue"], method="fdr_bh")[1]
    return pd.concat(donor_rows, ignore_index=True) if donor_rows else pd.DataFrame(), models, pd.DataFrame(null_rows), statuses


def run_composition(cache, eligibility):
    counts = pd.DataFrame({
        "donor": cache["donor"].astype(str), "cell_type": cache["cell_type"].astype(str),
        "n_cells": cache["n_cells"],
    })
    wide = counts.pivot(index="donor", columns="cell_type", values="n_cells")
    depth = pd.DataFrame({
        "donor": cache["composition_donor"].astype(str),
        "mean_log_umi": cache["composition_mean_log_umi"],
        "total_cells": cache["composition_n_cells"],
    }).set_index("donor").reindex(wide.index)
    donor_info = {d["donor"]: d for d in manifest_donors()}
    depth["old"] = [int(donor_info[d]["group"] != "young") for d in depth.index]
    clr = clr_transform(wide)
    rows, statuses = [], []
    for cell_type in eligibility["cell_type"]:
        eligible = bool(eligibility.loc[eligibility["cell_type"].eq(cell_type), "eligible"].iloc[0])
        for model in JIN_MODELS:
            if not eligible:
                statuses.append(status_row(COHORT, model, "CLR_composition", cell_type, "untestable", "failed_base_donor_nuclei_floors"))
                continue
            frame = depth.copy()
            frame["n_cells"] = wide[cell_type]
            fitted, _ = fit_one(frame.reset_index(drop=True), clr[cell_type].to_numpy(), model, COHORT)
            rows.append({"model": model, "endpoint": "CLR_composition", "cell_type": cell_type,
                         "composition_universe_n_cell_types": wide.shape[1], **fitted})
            statuses.append(status_row(COHORT, model, "CLR_composition", cell_type, "completed", "pass"))
    return add_bh(pd.DataFrame(rows), ["model"]), statuses


def run_genes(cache, eligibility):
    rows, statuses = [], []
    for cell_type in eligibility["cell_type"]:
        indices = selected_units(cache, eligibility, cell_type)
        if not len(indices):
            for model in JIN_MODELS:
                statuses.append(status_row(COHORT, model, "donor_pseudobulk_gene", cell_type, "untestable", "failed_base_donor_nuclei_floors"))
            continue
        frame = unit_frame(cache, indices)
        counts = cache["raw_counts"][indices].astype(float)
        totals = counts.sum(axis=1)
        expression = np.log1p(counts * (1e6 / totals[:, None]))
        gene_filter = (counts.sum(axis=0) >= 10) & ((counts > 0).sum(axis=0) >= 2)
        for model in JIN_MODELS:
            fitted, _ = fit_gene_models(frame, counts, expression[:, gene_filter], model, COHORT)
            testable = np.isfinite(fitted["pvalue"]) & np.isfinite(fitted["se"]) & (fitted["se"] > 0)
            gene_indices = np.flatnonzero(gene_filter)[testable]
            for position, gene_index in enumerate(gene_indices):
                rows.append({
                    "model": model, "endpoint": "donor_pseudobulk_gene", "cell_type": cell_type,
                    "gene_id": cache["gene_ids"][gene_index], "gene": cache["gene_names"][gene_index],
                    "effect": fitted["effect"][testable][position], "effect_definition": "reproductively_aged_minus_young",
                    "se": fitted["se"][testable][position], "ci95_low": fitted["ci95_low"][testable][position],
                    "ci95_high": fitted["ci95_high"][testable][position], "pvalue": fitted["pvalue"][testable][position],
                    "r2": fitted["r2"][testable][position], "residual_df": fitted["residual_df"][testable][position],
                    "n_donors": len(frame), "n_young": 4, "n_reproductively_aged": 4,
                    "n_cells": int(frame["n_cells"].sum()), "design_columns": fitted["design_columns"][0],
                    "effect_scale": "log1p_CPM_donor_summed_raw_counts",
                })
            statuses.append(status_row(COHORT, model, "donor_pseudobulk_gene", cell_type, "completed",
                                       f"tested_{int(testable.sum())}_of_{int(gene_filter.sum())}_filtered_genes"))
            print(f"Jin genes {cell_type} {model}: {int(testable.sum())}", flush=True)
    return add_bh(pd.DataFrame(rows), ["model"]), statuses


def prepare_hellinger(cache, eligibility):
    references, rows = [], []
    technical = technical_gene_mask(cache["gene_names"].astype(str))
    for cell_type in eligibility["cell_type"]:
        indices = selected_units(cache, eligibility, cell_type, qc=True)
        base = {"cell_type": cell_type, "n_donors_qc20": len(indices)}
        if not len(indices):
            rows.append({**base, "hellinger_eligible": False,
                         "reason": "failed_base_or_20_qualifying_nuclei_per_donor_floors"})
            continue
        young = indices[cache["old"][indices] == 0]
        total_reference = cache["qc_counts"][young].sum(axis=0)
        total_reference[technical] = 0
        total_young_cells = int(cache["qc_n"][young].sum())
        scored = []
        for index in indices:
            own = cache["qc_counts"][index].copy() if index in young else np.zeros(cache["qc_counts"].shape[1], dtype=np.int64)
            own[technical] = 0
            reference = total_reference - own
            own_cells = int(cache["qc_n"][index]) if index in young else 0
            reference_n = total_young_cells - own_cells
            if reference_n >= MIN_REFERENCE_CELLS and reference.sum() > 0:
                scored.append((index, reference_n, reference / reference.sum()))
        n_young = int((cache["old"][[x[0] for x in scored]] == 0).sum()) if scored else 0
        n_old = int((cache["old"][[x[0] for x in scored]] == 1).sum()) if scored else 0
        passed = n_young == 4 and n_old == 4
        rows.append({**base, "hellinger_eligible": passed,
                     "reason": "pass" if passed else "failed_all_four_per_group_or_50_cross_fitted_reference_nuclei",
                     "n_donors_after_reference_floors": len(scored),
                     "n_young_after_reference_floors": n_young,
                     "n_reproductively_aged_after_reference_floors": n_old,
                     "min_cross_fitted_reference_cells": min([x[1] for x in scored], default=0)})
        if passed:
            for index, reference_n, reference in scored:
                references.append({"cell_type": cell_type, "unit_index": int(index),
                                   "reference_n": int(reference_n), "reference": reference})
    return references, pd.DataFrame(rows)


def run_hellinger(cache, eligibility):
    references, eligibility_rows = prepare_hellinger(cache, eligibility)
    meta = load_metadata()
    donor_lookup = {d["donor"]: (i, d) for i, d in enumerate(manifest_donors(), 1)}
    technical = technical_gene_mask(cache["gene_names"].astype(str))
    distances_by_unit = {}
    for cell_type in sorted({r["cell_type"] for r in references}):
        refs = [r for r in references if r["cell_type"] == cell_type]
        by_donor = {str(cache["donor"][r["unit_index"]]): r for r in refs}
        for donor_name, ref in by_donor.items():
            donor_number, donor = donor_lookup[donor_name]
            _, matrix, annotations, _, _ = read_h5(donor, donor_number, meta, validate_counts=False)
            keep = annotations["celltype2"].eq(cell_type).to_numpy()
            matrix = matrix[keep].tocsr()
            depths = np.asarray(matrix.sum(axis=1)).ravel()
            nontechnical = np.asarray(matrix[:, ~technical].sum(axis=1)).ravel()
            valid = (depths >= MIN_UMI) & (nontechnical > 0)
            matrix = matrix[valid]
            depths = depths[valid]
            nontechnical = nontechnical[valid]
            if len(depths) != cache["qc_n"][ref["unit_index"]]:
                raise ValueError(f"Jin Hellinger QC count mismatch for {donor_name} {cell_type}")
            distances = hellinger_distances(matrix, ref["reference"], nontechnical)
            distances_by_unit[ref["unit_index"]] = {
                "hellinger": float(distances.mean()), "n_cells": len(distances),
                "mean_log_umi": float(np.log(depths).mean()), "mean_umi": float(depths.mean()),
                "cross_fitted_reference_cells": ref["reference_n"],
            }
            print(f"Jin Hellinger {cell_type} {donor_name}: {len(distances)} nuclei", flush=True)
    donor_rows = []
    for ref in references:
        index = ref["unit_index"]
        donor_rows.append({
            "donor": str(cache["donor"][index]), "cell_type": ref["cell_type"],
            "group": str(cache["group"][index]), "old": int(cache["old"][index]),
            **distances_by_unit[index],
        })
    donors = pd.DataFrame(donor_rows)
    model_rows, statuses = [], []
    for cell_type in eligibility["cell_type"]:
        part = donors[donors["cell_type"].eq(cell_type)] if not donors.empty else pd.DataFrame()
        for model in JIN_MODELS:
            if part.empty:
                reason = eligibility_rows.loc[eligibility_rows["cell_type"].eq(cell_type), "reason"].iloc[0]
                statuses.append(status_row(COHORT, model, "Hellinger_displacement", cell_type, "untestable", reason))
                continue
            fitted, _ = fit_one(part, part["hellinger"].to_numpy(), model, COHORT)
            model_rows.append({"model": model, "endpoint": "Hellinger_displacement", "cell_type": cell_type, **fitted})
            statuses.append(status_row(COHORT, model, "Hellinger_displacement", cell_type, "completed", "pass"))
    return donors, add_bh(pd.DataFrame(model_rows), ["model"]), eligibility_rows, statuses


def analyze():
    started = time.time()
    cache = load_cache()
    eligibility = pd.read_csv(RESULTS / f"{PREFIX}-celltype-eligibility.csv")
    p1_donors, p1_models, p1_null, p1_status = run_p1(cache, eligibility)
    write_csv(p1_donors, RESULTS / f"{PREFIX}-p1-donor-scores.csv")
    write_csv(p1_models, RESULTS / f"{PREFIX}-p1-models.csv")
    write_csv(p1_null, RESULTS / f"{PREFIX}-p1-null.csv")
    comp_models, comp_status = run_composition(cache, eligibility)
    write_csv(comp_models, RESULTS / f"{PREFIX}-composition-models.csv")
    gene_models, gene_status = run_genes(cache, eligibility)
    write_csv(gene_models, RESULTS / f"{PREFIX}-gene-models.csv")
    h_donors, h_models, h_eligibility, h_status = run_hellinger(cache, eligibility)
    write_csv(h_donors, RESULTS / f"{PREFIX}-hellinger-donors.csv")
    write_csv(h_models, RESULTS / f"{PREFIX}-hellinger-models.csv")
    write_csv(h_eligibility, RESULTS / f"{PREFIX}-hellinger-eligibility.csv")
    linkage_status = [
        status_row(COHORT, "not_applicable", "same_nucleus_RNA_ATAC_coupling", "all", "blocked", "separate_RNA_ATAC_assays_and_BioSamples", False),
        status_row(COHORT, "not_applicable", "donor_matched_RNA_ATAC_linkage", "all", "blocked", "no_valid_released_peak_gene_links_and_no_validated_EC_BEC_LEC_harmonization", False),
        status_row(COHORT, "not_applicable", "continuous_chronological_age", "all", "blocked", "not_planned_and_ovary5_47_49_conflict", False),
    ]
    write_csv(pd.DataFrame(p1_status + comp_status + gene_status + h_status + linkage_status),
              RESULTS / f"{PREFIX}-model-status.csv")
    resources = {
        "elapsed_seconds": time.time() - started,
        "peak_rss_bytes": int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024),
        "memory_limit_gb": 8,
        "minimum_free_gb": 30,
        "one_donor_at_a_time": True,
    }
    (RESULTS / f"{PREFIX}-analysis-resource.json").write_text(json.dumps(resources, indent=2, sort_keys=True) + "\n")


def finalize():
    manifest = json.loads(COHORTS.read_text())
    download = json.loads((RESULTS / f"followup-ovary-{COHORT}-download.json").read_text())
    by_donor = {row["donor"]: row for row in download["files"]}
    for donor in manifest[COHORT]["donors"]:
        local = donor_path(donor)
        donor["rna_filtered_count_matrix"].update({
            "downloaded": True,
            "local_path": str(Path("..") / local.relative_to(REPO.parent)),
            "sha256": by_donor[donor["donor"]]["sha256"],
            "count_values_validated": True,
        })
    verification = json.loads((RESULTS / f"{PREFIX}-verification.json").read_text())
    resources = json.loads((RESULTS / f"{PREFIX}-analysis-resource.json").read_text())
    manifest["scope"].update({
        "metadata_and_public_access_only": False,
        "expression_matrices_downloaded": True,
        "biological_outcomes_analyzed": True,
        "heavy_jobs_run": True,
        "expression_data_root": "../agescape-revision-data/ovary",
        "note": "Metadata remain under ../agescape-revision-data/ovary-metadata. Step 7 Jin per-donor H5 files and the GSE329417 final atlas were acquired sequentially under ../agescape-revision-data/ovary; raw-count gates and eligible donor-level outcomes are complete. No RAW tar was downloaded.",
    })
    manifest[COHORT]["raw_counts"]["validation_not_done"] = False
    manifest[COHORT]["raw_counts"]["validation"] = f"results/{PREFIX}-verification.json"
    manifest[COHORT]["analysis"] = {
        "completed": "2026-11-14",
        "frozen_plan": "data/manifest/ovary_jin_analysis_plan.json",
        "script": "scripts/run_followup_ovary_jin.py",
        "report": "reports/followup-ovary-analysis.md",
        "inferential_unit": "donor",
        "continuous_age_used": False,
        "ovary5_conflict_preserved": True,
        "annotated_nuclei": verification["annotated_nuclei"],
        "all_raw_counts_validated": True,
        "same_nucleus_coupling": "blocked",
        "donor_matched_linkage": "blocked_no_valid_peak_gene_links_or_label_harmonization",
        "peak_rss_bytes_build_cache": verification["peak_rss_bytes"],
        "peak_rss_bytes_analysis": resources["peak_rss_bytes"],
    }
    COHORTS.write_text(json.dumps(manifest, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["build-cache", "analyze", "finalize"])
    args = parser.parse_args()
    RESULTS.mkdir(exist_ok=True)
    if args.command == "build-cache":
        build_cache()
    elif args.command == "analyze":
        analyze()
    else:
        finalize()


if __name__ == "__main__":
    main()
