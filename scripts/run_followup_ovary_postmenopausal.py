#!/usr/bin/env python3
"""Validate and analyze GSE329417 within-postmenopausal donor endpoints."""

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
    clr_transform, fit_ols, known_cell_type, matched_null_effect_pvalue,
    random_signature_scores,
)
from src.followup_ovary import (  # noqa: E402
    POST_MODELS, build_postmenopausal_design, fit_gene_models, fit_one,
    sha256, status_row, write_csv,
)

DATA = REPO.parent / f"{REPO.name}-data" / "ovary" / "GSE329417"
ATLAS = DATA / "GSE329417_human_ovary_snucRNAseq_final_annotated_atlas.h5ad"
CACHE = DATA / "interim" / "postmenopausal_donor_aggregates.npz"
METADATA = REPO.parent / f"{REPO.name}-data" / "ovary-metadata" / "GSE329417"
COHORTS = REPO / "data" / "manifest" / "ovary_cohorts.json"
RESULTS = REPO / "results"
PREFIX = "followup-ovary-postmenopausal"
COHORT = "GSE329417"
CHUNK = 256
MIN_CELLS = 20
MIN_DONORS = 20
EXPECTED_NUCLEI = 439011
DONOR_CANDIDATES = ["SenNet_ID", "donor_id", "Donor_ID", "donor"]
SAMPLE_CANDIDATES = ["unique_sample_ID", "processing_sample_ID", "sample_id", "Sample_ID"]
ANNOTATION_CANDIDATES = [
    "cell_type", "celltype", "broad_cell_type", "coarse_cell_type",
    "cell_type_coarse", "cell_type_coarse1", "CellType", "annotation", "cell_type_l1",
]
GENE_NAME_CANDIDATES = ["feature_name", "gene_name", "Gene", "gene_symbol", "symbol"]


def cohort_manifest():
    return json.loads(COHORTS.read_text())[COHORT]


def atlas_sample_identifier(sample):
    entry = sample["cellbender_h5ad"]
    filename = entry["url"].rsplit("/", 1)[-1]
    prefix = f"{sample['geo_library_gsm']}_"
    suffix = "_cellbender.h5ad"
    if not filename.startswith(prefix) or not filename.endswith(suffix):
        raise ValueError(f"Unexpected released CellBender filename: {filename}")
    stem, plex = filename[len(prefix):-len(suffix)].rsplit("_p", 1)
    if int(plex) != int(sample["plex"]):
        raise ValueError(f"CellBender filename/processing plex mismatch: {filename}")
    return f"{stem}_{plex}"


def clinical_donors():
    path = METADATA / "GSE329417_human_ovary_snucRNAseq_public_donor_metadata.tsv.gz"
    donors = pd.read_csv(path, sep="\t").rename(columns={
        "SenNet_ID": "donor", "Age": "age_years", "Para": "parity_live_births",
        "NonOvarian_Cancer_State": "non_ovarian_cancer_state",
        "Tissue_Pieces_Analyzed": "expected_tissue_pieces",
    })
    if len(donors) != 28 or donors["donor"].nunique() != 28:
        raise ValueError("GSE329417 clinical table is not 28 unique donors")
    donors["non_ovarian_abnormal"] = donors["non_ovarian_cancer_state"].ne("Benign").astype(int)
    if set(donors["non_ovarian_cancer_state"]) != {"Benign", "Pre-malignant", "Malignant"}:
        raise ValueError("Unexpected non-ovarian cancer state labels")
    if donors[["age_years", "BMI", "parity_live_births", "non_ovarian_abnormal"]].isna().any().any():
        raise ValueError("Frozen clinical covariates contain missing values")
    expected = {d["donor"]: d for d in cohort_manifest()["donors"]}
    for row in donors.itertuples():
        entry = expected[row.donor]
        if int(row.age_years) != entry["exact_age_years"] or float(row.BMI) != float(entry["BMI"]):
            raise ValueError(f"Clinical/manifest mismatch for {row.donor}")
        if int(row.parity_live_births) != int(entry["parity_live_births"]):
            raise ValueError(f"Parity/manifest mismatch for {row.donor}")
        if row.non_ovarian_cancer_state != entry["non_ovarian_cancer_state"]:
            raise ValueError(f"Cancer-state/manifest mismatch for {row.donor}")
    return donors.sort_values("donor").reset_index(drop=True)


def resolve_mapping_field(obs, candidates, valid_values, expected_levels, label):
    valid_values = set(valid_values)
    scores = {}
    for column in obs.columns:
        values = obs[column].astype(str).str.strip()
        scores[column] = float(values.isin(valid_values).mean())
    for column in candidates:
        if column in scores and scores[column] == 1.0:
            values = obs[column].astype(str).str.strip()
            if values.nunique() == expected_levels:
                return column, scores
    complete = [column for column, score in scores.items() if score == 1.0
                and obs[column].astype(str).str.strip().nunique() == expected_levels]
    if len(complete) != 1:
        raise ValueError(f"No unique exact atlas {label} mapping field: candidates={complete}, best={sorted(scores.items(), key=lambda x: -x[1])[:8]}")
    return complete[0], scores


def resolve_annotation_field(obs):
    candidates = []
    for column in ANNOTATION_CANDIDATES:
        if column not in obs:
            continue
        values = obs[column].astype(str).str.strip()
        known = values[known_cell_type(values)]
        if known.nunique() == 21:
            candidates.append(column)
    if len(candidates) != 1:
        raise ValueError(f"Expected one released 21-level coarse annotation field; found {candidates}")
    return candidates[0]


def read_obs_column(obs_group, column):
    item = obs_group[column]
    if isinstance(item, h5py.Group) and "categories" in item:
        categories = item["categories"][:].astype(str)
        codes = item["codes"][:]
        values = categories[np.maximum(codes, 0)]
        values[codes < 0] = "nan"
        return values
    return item[:]


def read_csr_chunk(group, start, stop, n_vars):
    indptr = group["indptr"][start:stop + 1].astype(np.int64)
    first, last = int(indptr[0]), int(indptr[-1])
    data = group["data"][first:last]
    indices = group["indices"][first:last]
    return sp.csr_matrix((data, indices, indptr - first), shape=(stop - start, n_vars))


def load_atlas_metadata():
    handle = h5py.File(ATLAS, "r")
    try:
        n_obs = len(handle["X/indptr"]) - 1
        obs_names = handle["obs/_index"][:].astype(str)
        if n_obs != EXPECTED_NUCLEI or not pd.Index(obs_names).is_unique:
            raise ValueError("Atlas nucleus count or observation-ID uniqueness gate failed")
        selected = set(DONOR_CANDIDATES + SAMPLE_CANDIDATES + ANNOTATION_CANDIDATES + [
            "Age", "BMI", "Para", "NonOvarian_Cancer_State", "plex",
        ])
        obs = pd.DataFrame({
            column: read_obs_column(handle["obs"], column)
            for column in sorted(selected.intersection(handle["obs"].keys()))
        }, index=obs_names)
        donors = clinical_donors()
        samples = cohort_manifest()["samples"]
        donor_field, _ = resolve_mapping_field(obs, DONOR_CANDIDATES, donors["donor"], 28, "donor")
        atlas_to_processing = {atlas_sample_identifier(s): s["processing_sample_id"] for s in samples}
        if len(atlas_to_processing) != 64:
            raise ValueError("Released CellBender filenames do not produce 64 unique atlas sample IDs")
        sample_field, _ = resolve_mapping_field(obs, SAMPLE_CANDIDATES, atlas_to_processing, 64, "sample")
        annotation_field = resolve_annotation_field(obs)
        obs["donor"] = obs[donor_field].astype(str).str.strip()
        obs["atlas_sample_id"] = obs[sample_field].astype(str).str.strip()
        obs["processing_sample_id"] = obs["atlas_sample_id"].map(atlas_to_processing)
        sample_map = {s["processing_sample_id"]: s for s in samples}
        obs["mapped_donor"] = obs["processing_sample_id"].map(lambda x: sample_map[x]["donor"])
        obs["released_plex"] = obs["plex"].to_numpy(int)
        obs["plex"] = obs["processing_sample_id"].map(lambda x: int(sample_map[x]["plex"])).astype(int)
        if not obs["donor"].equals(obs["mapped_donor"]):
            raise ValueError("Atlas donor and sample-to-donor mappings disagree")
        if not np.array_equal(obs["released_plex"], obs["plex"]):
            raise ValueError("Atlas and processing-log plex mappings disagree")
        obs["cell_type"] = obs[annotation_field].astype(str).str.strip()
        obs["known_annotation"] = known_cell_type(obs["cell_type"])
        if obs.loc[obs["known_annotation"], "cell_type"].nunique() != 21:
            raise ValueError("Released coarse annotation coverage is not 21 labels")
        obs["released_BMI"] = obs["BMI"].to_numpy(float)
        clinical = donors.set_index("donor")
        for column in ["age_years", "BMI", "parity_live_births", "non_ovarian_cancer_state", "non_ovarian_abnormal"]:
            obs[column] = obs["donor"].map(clinical[column])
        if obs[["age_years", "BMI", "parity_live_births", "non_ovarian_abnormal"]].isna().any().any():
            raise ValueError("Atlas clinical join is incomplete")
        if not np.array_equal(obs["Age"].to_numpy(int), obs["age_years"].to_numpy(int)):
            raise ValueError("Atlas and public exact ages disagree")
        if not np.allclose(obs["released_BMI"].to_numpy(float), obs["BMI"].to_numpy(float)):
            raise ValueError("Atlas and public BMI values disagree")
        if not np.array_equal(obs["Para"].to_numpy(int), obs["parity_live_births"].to_numpy(int)):
            raise ValueError("Atlas and public parity values disagree")
        if not obs["NonOvarian_Cancer_State"].astype(str).equals(obs["non_ovarian_cancer_state"].astype(str)):
            raise ValueError("Atlas and public non-ovarian cancer states disagree")
        for donor, part in obs.groupby("donor"):
            if part["age_years"].nunique() != 1:
                raise ValueError(f"Non-unique exact age for {donor}")
        obs["_row"] = np.arange(len(obs))
        join = {
            "donor_field": donor_field,
            "sample_field": sample_field,
            "coarse_annotation_field": annotation_field,
            "sample_crosswalk": "released per-sample CellBender filename stem plus its _p1-_p4 suffix; exact one-to-one map to processing_sample_ID",
            "all_observations_exactly_mapped": True,
            "all_64_samples_present": obs["processing_sample_id"].nunique() == 64,
            "all_28_donors_present": obs["donor"].nunique() == 28,
            "technical_replicates_combined_within_donor": ["3_4", "16_4"],
        }
        if not join["all_64_samples_present"] or not join["all_28_donors_present"]:
            raise ValueError("Atlas donor or sample coverage is incomplete")
        return handle, obs, donors, join
    except Exception:
        handle.close()
        raise


def unit_frame_from_obs(obs):
    counts = obs[obs["known_annotation"]].groupby(["cell_type", "donor"]).size().rename("n_cells").reset_index()
    rows = []
    for cell_type, part in counts.groupby("cell_type"):
        eligible = part[part["n_cells"] >= MIN_CELLS]
        rows.append({
            "cell_type": cell_type,
            "n_donors_observed": len(part),
            "n_donors_ge20": len(eligible),
            "n_cells_total": int(part["n_cells"].sum()),
            "eligible": len(eligible) >= MIN_DONORS,
            "reason": "pass" if len(eligible) >= MIN_DONORS else "failed_20_donors_with_20_nuclei_floor",
        })
    return pd.DataFrame(rows).sort_values("cell_type").reset_index(drop=True)


def build_cache():
    started = time.time()
    handle, obs, donors, join = load_atlas_metadata()
    try:
        eligibility = unit_frame_from_obs(obs)
        write_csv(eligibility, RESULTS / f"{PREFIX}-celltype-eligibility.csv")
        cell_types = sorted(obs.loc[obs["known_annotation"], "cell_type"].unique())
        donor_ids = donors["donor"].tolist()
        units = pd.MultiIndex.from_product([donor_ids, cell_types], names=["donor", "cell_type"]).to_frame(index=False)
        clinical = donors.set_index("donor")
        for column in ["age_years", "BMI", "parity_live_births", "non_ovarian_abnormal"]:
            units[column] = units["donor"].map(clinical[column])
        unit_lookup = {(r.donor, r.cell_type): i for i, r in enumerate(units.itertuples())}
        donor_lookup = {donor: i for i, donor in enumerate(donor_ids)}
        known = obs["known_annotation"].to_numpy(bool)
        codes = np.full(len(obs), -1, dtype=np.int64)
        codes[known] = [unit_lookup[(d, c)] for d, c in zip(obs.loc[known, "donor"], obs.loc[known, "cell_type"])]
        donor_codes = np.full(len(obs), -1, dtype=np.int64)
        donor_codes[known] = obs.loc[known, "donor"].map(donor_lookup).to_numpy()
        plex_codes = obs["plex"].to_numpy(int) - 1
        gene_names = handle["var/_index"][:].astype(str)
        gene_ids = (handle["var/gene_ids"][:].astype(str)
                    if "gene_ids" in handle["var"] else gene_names.copy())
        if not pd.Index(gene_ids).is_unique:
            raise ValueError("Atlas feature IDs are not unique")
        gene_name_field = "var/_index"
        if not (pd.Index(gene_names).str.len() > 0).all():
            raise ValueError("No released gene symbols are available for P1 mapping")
        n_obs = len(obs)
        n_units, n_genes, n_donors = len(units), len(gene_ids), len(donor_ids)
        raw_counts = np.zeros((n_units, n_genes), dtype=np.int64)
        log_sums = np.zeros((n_units, n_genes), dtype=np.float32)
        n_cells = np.zeros(n_units, dtype=np.int64)
        sum_log_umi = np.zeros(n_units)
        sum_umi = np.zeros(n_units)
        unit_plex = np.zeros((n_units, 4), dtype=np.int64)
        composition_n = np.zeros(n_donors, dtype=np.int64)
        composition_log = np.zeros(n_donors)
        composition_umi = np.zeros(n_donors)
        composition_plex = np.zeros((n_donors, 4), dtype=np.int64)
        if "raw_counts" not in handle["layers"] or "log1p_norm" not in handle["layers"] or "X" not in handle:
            raise ValueError("Atlas raw_counts, log1p_norm, or X is absent")
        raw_layer = handle["layers/raw_counts"]
        normalized_layer = handle["layers/log1p_norm"]
        raw_nnz = 0
        raw_sum = 0.0
        for start in range(0, n_obs, CHUNK):
            stop = min(start + CHUNK, n_obs)
            matrix = read_csr_chunk(raw_layer, start, stop, n_genes)
            assert_count_gate(matrix.data)
            x_matrix = read_csr_chunk(handle["X"], start, stop, n_genes)
            normalized = read_csr_chunk(normalized_layer, start, stop, n_genes)
            if x_matrix.shape != normalized.shape or (x_matrix != normalized).nnz:
                raise ValueError(f"Atlas X and log1p_norm differ in rows {start}:{stop}")
            raw_nnz += matrix.nnz
            raw_sum += float(matrix.sum(dtype=np.float64))
            local = codes[start:stop]
            chosen = np.flatnonzero(local >= 0)
            if len(chosen):
                sub = matrix[chosen].astype(np.float64)
                sub_codes = local[chosen]
                libraries = np.asarray(sub.sum(axis=1)).ravel()
                if (libraries <= 0).any():
                    raise ValueError("Zero-UMI annotated atlas nucleus")
                np.add.at(n_cells, sub_codes, 1)
                np.add.at(sum_log_umi, sub_codes, np.log(libraries))
                np.add.at(sum_umi, sub_codes, libraries)
                np.add.at(unit_plex, (sub_codes, plex_codes[start:stop][chosen]), 1)
                log_sub = sub.copy()
                log_sub.data *= np.repeat(1e6 / libraries, np.diff(log_sub.indptr))
                np.log1p(log_sub.data, out=log_sub.data)
                aggregate_block(log_sub, sub_codes, log_sums)
                aggregate_block(sub, sub_codes, raw_counts)
                local_donors = donor_codes[start:stop][chosen]
                np.add.at(composition_n, local_donors, 1)
                np.add.at(composition_log, local_donors, np.log(libraries))
                np.add.at(composition_umi, local_donors, libraries)
                np.add.at(composition_plex, (local_donors, plex_codes[start:stop][chosen]), 1)
            if start % 10240 == 0:
                print(f"atlas rows {start}:{stop}; {raw_nnz:,} corrected values checked", flush=True)
        cache = {
            "gene_ids": gene_ids, "gene_names": gene_names,
            "donor": units["donor"].astype(str).to_numpy(dtype=str),
            "cell_type": units["cell_type"].astype(str).to_numpy(dtype=str),
            "age_years": units["age_years"].to_numpy(float),
            "BMI": units["BMI"].to_numpy(float),
            "parity_live_births": units["parity_live_births"].to_numpy(float),
            "non_ovarian_abnormal": units["non_ovarian_abnormal"].to_numpy(np.int8),
            "n_cells": n_cells,
            "mean_log_umi": np.divide(sum_log_umi, n_cells, out=np.full(n_units, np.nan), where=n_cells > 0),
            "mean_umi": np.divide(sum_umi, n_cells, out=np.full(n_units, np.nan), where=n_cells > 0),
            "raw_counts": raw_counts, "log_sums": log_sums, "unit_plex_counts": unit_plex,
            "composition_donor": np.asarray(donor_ids),
            "composition_n_cells": composition_n,
            "composition_mean_log_umi": composition_log / composition_n,
            "composition_mean_umi": composition_umi / composition_n,
            "composition_plex_counts": composition_plex,
        }
        CACHE.parent.mkdir(parents=True, exist_ok=True)
        np.savez(CACHE, **cache)
        counts = pd.DataFrame({
            "donor": cache["donor"], "cell_type": cache["cell_type"], "n_cells": n_cells,
            "age_years": cache["age_years"], "mean_log_umi": cache["mean_log_umi"],
            "plex_1_n": unit_plex[:, 0], "plex_2_n": unit_plex[:, 1],
            "plex_3_n": unit_plex[:, 2], "plex_4_n": unit_plex[:, 3],
        })
        write_csv(counts, RESULTS / f"{PREFIX}-celltype-donor-counts.csv")
        sample_qc = obs.groupby(["processing_sample_id", "donor", "plex"]).size().rename("atlas_nuclei").reset_index()
        expected_samples = pd.DataFrame([{
            "processing_sample_id": s["processing_sample_id"],
            "technical_replicate_group": s["technical_replicate_group"],
        } for s in cohort_manifest()["samples"]])
        sample_qc = sample_qc.merge(expected_samples, on="processing_sample_id", validate="one_to_one")
        write_csv(sample_qc, RESULTS / f"{PREFIX}-sample-qc.csv")
        donor_qc = donors.copy()
        observed = obs.groupby("donor").agg(atlas_nuclei=("donor", "size"), observed_samples=("processing_sample_id", "nunique")).reset_index()
        donor_qc = donor_qc.merge(observed, on="donor", validate="one_to_one")
        donor_qc["composition_n_cells"] = donor_qc["donor"].map(dict(zip(donor_ids, composition_n)))
        donor_qc["mean_log_umi"] = donor_qc["donor"].map(dict(zip(donor_ids, cache["composition_mean_log_umi"])))
        if not donor_qc["observed_samples"].eq(donor_qc["expected_tissue_pieces"]).all():
            raise ValueError("Atlas tissue-piece counts disagree with donor metadata")
        write_csv(donor_qc, RESULTS / f"{PREFIX}-donor-qc.csv")
        identifiability = build_identifiability(cache, eligibility)
        write_csv(identifiability, RESULTS / f"{PREFIX}-covariate-identifiability.csv")
        verification = {
            "cohort": COHORT,
            "atlas_nuclei": int(n_obs),
            "atlas_genes": int(n_genes),
            "donors": len(donor_ids),
            "samples": int(obs["processing_sample_id"].nunique()),
            "coarse_cell_types": len(cell_types),
            "unknown_annotation_nuclei_excluded": int((~known).sum()),
            "stored_corrected_count_values_checked": int(raw_nnz),
            "corrected_count_sum": raw_sum,
            "all_counts_finite_nonnegative_integer_valued": True,
            "x_identical_to_log1p_norm_in_every_chunk": True,
            "x_used_as_raw_counts": False,
            "uncorrected_cell_ranger_counts_available": False,
            "uncorrected_sensitivity_performed": False,
            "gene_name_field": gene_name_field or "var_names",
            "feature_ids_unique": True,
            "duplicate_gene_name_rows": int(pd.Index(gene_names).duplicated(keep=False).sum()),
            "join": join,
            "cache_bytes": CACHE.stat().st_size,
            "elapsed_seconds": time.time() - started,
            "peak_rss_bytes": int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024),
            "memory_limit_gb": 8,
            "minimum_free_gb": 30,
            "outcomes_examined": False,
        }
        (RESULTS / f"{PREFIX}-verification.json").write_text(json.dumps(verification, indent=2, sort_keys=True) + "\n")
    finally:
        handle.close()
    print("GSE329417 corrected-count and atlas-join gates passed; no outcomes fitted in build-cache.", flush=True)


def load_cache():
    with np.load(CACHE, allow_pickle=False) as data:
        return {key: data[key] for key in data.files}


def add_plex_fractions(frame, counts, denominator=None):
    if denominator is None:
        denominator = frame["n_cells"].to_numpy(float)
    for plex in range(1, 5):
        frame[f"plex_fraction_{plex}"] = counts[:, plex - 1] / denominator
    return frame


def unit_frame(cache, indices):
    frame = pd.DataFrame({
        "donor": cache["donor"][indices].astype(str),
        "cell_type": cache["cell_type"][indices].astype(str),
        "age_years": cache["age_years"][indices],
        "BMI": cache["BMI"][indices],
        "parity_live_births": cache["parity_live_births"][indices],
        "non_ovarian_abnormal": cache["non_ovarian_abnormal"][indices],
        "n_cells": cache["n_cells"][indices],
        "mean_log_umi": cache["mean_log_umi"][indices],
    })
    return add_plex_fractions(frame, cache["unit_plex_counts"][indices])


def selected_units(cache, eligibility, cell_type):
    row = eligibility[eligibility["cell_type"].eq(cell_type)]
    if row.empty or not bool(row["eligible"].iloc[0]):
        return np.array([], dtype=int)
    indices = np.flatnonzero((cache["cell_type"].astype(str) == cell_type) & (cache["n_cells"] >= MIN_CELLS))
    if len(indices) < MIN_DONORS:
        return np.array([], dtype=int)
    return indices[np.argsort(cache["donor"][indices].astype(str))]


def composition_frame(cache, frame):
    lookup = {donor: i for i, donor in enumerate(cache["composition_donor"].astype(str))}
    indices = np.asarray([lookup[d] for d in frame["donor"]])
    result = frame[["donor", "cell_type", "age_years", "BMI", "parity_live_births", "non_ovarian_abnormal", "n_cells"]].copy()
    result["mean_log_umi"] = cache["composition_mean_log_umi"][indices]
    return add_plex_fractions(result, cache["composition_plex_counts"][indices], cache["composition_n_cells"][indices])


def build_identifiability(cache, eligibility):
    rows = []
    for cell_type in eligibility["cell_type"]:
        indices = selected_units(cache, eligibility, cell_type)
        if not len(indices):
            for endpoint_class in ["expression", "composition"]:
                for model in POST_MODELS:
                    rows.append({"cell_type": cell_type, "endpoint_class": endpoint_class, "model": model,
                                 "identifiable": False, "n_donors": 0, "reason": "failed_base_donor_nuclei_floors"})
            continue
        expression = unit_frame(cache, indices)
        for endpoint_class in ["expression", "composition"]:
            frame = expression if endpoint_class == "expression" else composition_frame(cache, expression)
            for model in POST_MODELS:
                try:
                    design, names = build_postmenopausal_design(frame, model)
                    rows.append({"cell_type": cell_type, "endpoint_class": endpoint_class, "model": model,
                                 "identifiable": True, "n_donors": len(frame), "n_design_columns": design.shape[1],
                                 "residual_df": len(frame) - design.shape[1], "design_columns": ";".join(names),
                                 "reason": "pass"})
                except ValueError as exc:
                    rows.append({"cell_type": cell_type, "endpoint_class": endpoint_class, "model": model,
                                 "identifiable": False, "n_donors": len(frame), "reason": str(exc)})
    return pd.DataFrame(rows)


def identifiable(identifiability, cell_type, endpoint_class, model):
    rows = identifiability[(identifiability["cell_type"] == cell_type)
                           & (identifiability["endpoint_class"] == endpoint_class)
                           & (identifiability["model"] == model)]
    return not rows.empty and bool(rows["identifiable"].iloc[0])


def run_p1(cache, eligibility, identifiability):
    program = load_program_tiers(REPO, "human")["P1_Clean"]
    gene_names = cache["gene_names"].astype(str)
    donor_rows, model_rows, null_rows, statuses = [], [], [], []
    for cell_type in eligibility["cell_type"]:
        indices = selected_units(cache, eligibility, cell_type)
        if not len(indices):
            for model in POST_MODELS:
                statuses.append(status_row(COHORT, model, "P1_transport", cell_type, "untestable", "failed_base_donor_nuclei_floors"))
            continue
        frame = unit_frame(cache, indices)
        matrix = (cache["log_sums"][indices] / cache["n_cells"][indices, None]).astype(float).T
        detected = cache["raw_counts"][indices].sum(axis=0) > 0
        up = np.flatnonzero(np.isin(gene_names, program["up"]) & detected)
        down = np.flatnonzero(np.isin(gene_names, program["down"]) & detected)
        if not len(up) or not len(down):
            for model in POST_MODELS:
                statuses.append(status_row(COHORT, model, "P1_transport", cell_type, "untestable", "no_detected_P1_up_or_down_program_genes"))
            continue
        scores = matrix[up].mean(axis=0) - matrix[down].mean(axis=0)
        donor_frame = frame.copy()
        donor_frame["P1_score"] = scores
        donor_frame["p1_up_detected"] = len(up)
        donor_frame["p1_down_detected"] = len(down)
        donor_rows.append(donor_frame)
        for model in POST_MODELS:
            if not identifiable(identifiability, cell_type, "expression", model):
                statuses.append(status_row(COHORT, model, "P1_transport", cell_type, "untestable", "frozen_covariate_design_not_identifiable"))
                continue
            fitted, design = fit_one(frame, scores, model, COHORT)
            row = {"model": model, "endpoint": "P1_transport", "cell_type": cell_type,
                   "p1_up_detected": len(up), "p1_down_detected": len(down), **fitted}
            if model == "primary_depth_plex_clinical":
                try:
                    random_scores = random_signature_scores(matrix, up, down)
                    null_effects = fit_ols(random_scores.T, design)["effect"]
                    row.update({
                        "random_iterations": NULL_ITERATIONS, "random_seed": NULL_SEED,
                        "random_mean_effect": float(null_effects.mean()),
                        "random_sd_effect": float(null_effects.std(ddof=1)),
                        "random_q025": float(np.quantile(null_effects, 0.025)),
                        "random_q975": float(np.quantile(null_effects, 0.975)),
                        "random_pvalue": matched_null_effect_pvalue(fitted["effect"], null_effects),
                    })
                    for iteration, effect in enumerate(null_effects, 1):
                        null_rows.append({"cell_type": cell_type, "iteration": iteration,
                                          "effect_per_10_exact_years": effect, "seed": NULL_SEED})
                except ValueError as exc:
                    row["random_null_status"] = f"untestable: {exc}"
            model_rows.append(row)
            statuses.append(status_row(COHORT, model, "P1_transport", cell_type, "completed", "pass"))
    models = add_bh(pd.DataFrame(model_rows), ["model"])
    models["random_q_bh"] = np.nan
    primary = models["model"].eq("primary_depth_plex_clinical") & models["random_pvalue"].notna() if "random_pvalue" in models else pd.Series(False, index=models.index)
    if primary.any():
        models.loc[primary, "random_q_bh"] = multipletests(models.loc[primary, "random_pvalue"], method="fdr_bh")[1]
    return pd.concat(donor_rows, ignore_index=True) if donor_rows else pd.DataFrame(), models, pd.DataFrame(null_rows), statuses


def run_composition(cache, eligibility, identifiability):
    rows, statuses = [], []
    for cell_type in eligibility["cell_type"]:
        indices = selected_units(cache, eligibility, cell_type)
        if not len(indices):
            for model in POST_MODELS:
                statuses.append(status_row(COHORT, model, "CLR_composition", cell_type, "untestable", "failed_base_donor_nuclei_floors"))
            continue
        expression = unit_frame(cache, indices)
        frame = composition_frame(cache, expression)
        all_counts = pd.DataFrame({
            "donor": cache["donor"].astype(str), "cell_type": cache["cell_type"].astype(str),
            "n_cells": cache["n_cells"],
        })
        wide = all_counts.pivot(index="donor", columns="cell_type", values="n_cells")
        clr = clr_transform(wide)
        outcome = clr.loc[frame["donor"], cell_type].to_numpy()
        for model in POST_MODELS:
            if not identifiable(identifiability, cell_type, "composition", model):
                statuses.append(status_row(COHORT, model, "CLR_composition", cell_type, "untestable", "frozen_covariate_design_not_identifiable"))
                continue
            fitted, _ = fit_one(frame, outcome, model, COHORT)
            rows.append({"model": model, "endpoint": "CLR_composition", "cell_type": cell_type,
                         "composition_universe_n_cell_types": wide.shape[1], **fitted})
            statuses.append(status_row(COHORT, model, "CLR_composition", cell_type, "completed", "pass"))
    return add_bh(pd.DataFrame(rows), ["model"]), statuses


def run_genes(cache, eligibility, identifiability):
    statuses = []
    for model in POST_MODELS:
        rows = []
        for cell_type in eligibility["cell_type"]:
            indices = selected_units(cache, eligibility, cell_type)
            if not len(indices):
                statuses.append(status_row(COHORT, model, "donor_pseudobulk_gene", cell_type, "untestable", "failed_base_donor_nuclei_floors"))
                continue
            if not identifiable(identifiability, cell_type, "expression", model):
                statuses.append(status_row(COHORT, model, "donor_pseudobulk_gene", cell_type, "untestable", "frozen_covariate_design_not_identifiable"))
                continue
            frame = unit_frame(cache, indices)
            counts = cache["raw_counts"][indices].astype(float)
            totals = counts.sum(axis=1)
            expression = np.log1p(counts * (1e6 / totals[:, None]))
            gene_filter = (counts.sum(axis=0) >= 10) & ((counts > 0).sum(axis=0) >= 2)
            fitted, _ = fit_gene_models(frame, counts, expression[:, gene_filter], model, COHORT)
            testable = np.isfinite(fitted["pvalue"]) & np.isfinite(fitted["se"]) & (fitted["se"] > 0)
            gene_indices = np.flatnonzero(gene_filter)[testable]
            for position, gene_index in enumerate(gene_indices):
                rows.append({
                    "model": model, "endpoint": "donor_pseudobulk_gene", "cell_type": cell_type,
                    "gene_id": cache["gene_ids"][gene_index], "gene": cache["gene_names"][gene_index],
                    "effect": fitted["effect"][testable][position], "effect_definition": "within_postmenopausal_per_10_exact_years",
                    "se": fitted["se"][testable][position], "ci95_low": fitted["ci95_low"][testable][position],
                    "ci95_high": fitted["ci95_high"][testable][position], "pvalue": fitted["pvalue"][testable][position],
                    "r2": fitted["r2"][testable][position], "residual_df": fitted["residual_df"][testable][position],
                    "n_donors": len(frame), "n_cells": int(frame["n_cells"].sum()),
                    "age_min": float(frame["age_years"].min()), "age_max": float(frame["age_years"].max()),
                    "design_columns": fitted["design_columns"][0],
                    "effect_scale": "log1p_CPM_donor_summed_CellBender_corrected_counts",
                })
            statuses.append(status_row(COHORT, model, "donor_pseudobulk_gene", cell_type, "completed",
                                       f"tested_{int(testable.sum())}_of_{int(gene_filter.sum())}_filtered_genes"))
            print(f"postmenopausal genes {cell_type} {model}: {int(testable.sum())}", flush=True)
        output = add_bh(pd.DataFrame(rows), ["model"])
        write_csv(output, RESULTS / f"{PREFIX}-gene-models-{model}.csv")
    return statuses


def analyze():
    started = time.time()
    cache = load_cache()
    eligibility = pd.read_csv(RESULTS / f"{PREFIX}-celltype-eligibility.csv")
    identifiability = pd.read_csv(RESULTS / f"{PREFIX}-covariate-identifiability.csv")
    p1_donors, p1_models, p1_null, p1_status = run_p1(cache, eligibility, identifiability)
    write_csv(p1_donors, RESULTS / f"{PREFIX}-p1-donor-scores.csv")
    write_csv(p1_models, RESULTS / f"{PREFIX}-p1-models.csv")
    write_csv(p1_null, RESULTS / f"{PREFIX}-p1-null.csv")
    comp_models, comp_status = run_composition(cache, eligibility, identifiability)
    write_csv(comp_models, RESULTS / f"{PREFIX}-composition-models.csv")
    gene_status = run_genes(cache, eligibility, identifiability)
    blocked = []
    for cell_type in eligibility["cell_type"]:
        blocked.append(status_row(COHORT, "not_applicable", "Hellinger_displacement", cell_type,
                                  "untestable", "no_20_to_40_year_reference_in_postmenopausal_cohort", False))
    blocked.extend([
        status_row(COHORT, "not_applicable", "uncorrected_Cell_Ranger_sensitivity", "all",
                   "blocked", "uncorrected_Cell_Ranger_counts_not_separately_released", False),
        status_row(COHORT, "not_applicable", "premenopausal_or_Jin_pooled_analysis", "all",
                   "blocked", "blocked_by_design", False),
    ])
    write_csv(pd.DataFrame(p1_status + comp_status + gene_status + blocked), RESULTS / f"{PREFIX}-model-status.csv")
    resources = {
        "elapsed_seconds": time.time() - started,
        "peak_rss_bytes": int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024),
        "memory_limit_gb": 8,
        "minimum_free_gb": 30,
        "donor_aggregated_models_only": True,
    }
    (RESULTS / f"{PREFIX}-analysis-resource.json").write_text(json.dumps(resources, indent=2, sort_keys=True) + "\n")


def finalize():
    manifest = json.loads(COHORTS.read_text())
    download = json.loads((RESULTS / f"followup-ovary-{COHORT}-download.json").read_text())
    decompression = json.loads((RESULTS / f"{PREFIX}-decompression.json").read_text())
    verification = json.loads((RESULTS / f"{PREFIX}-verification.json").read_text())
    resources = json.loads((RESULTS / f"{PREFIX}-analysis-resource.json").read_text())
    atlas_entry = next(row for row in manifest[COHORT]["public_processed_objects"]
                       if row["url"].endswith("final_annotated_atlas.h5ad.gz"))
    atlas_entry.update({
        "downloaded": True,
        "content": "Final post-QC annotated atlas; CellBender-corrected raw_counts and log-normalized layers, public obs/var, annotations and embeddings; downloaded and validated in Step 7",
        "local_path": str(Path("..") / Path(download["files"][0]["local_path"]).relative_to(REPO.parent)),
        "sha256": download["files"][0]["sha256"],
        "decompressed_local_path": str(Path("..") / ATLAS.relative_to(REPO.parent)),
        "decompressed_size_bytes": decompression["decompressed_size_bytes"],
        "decompressed_sha256": decompression["decompressed_sha256"],
        "count_values_validated": True,
    })
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
    manifest[COHORT]["raw_counts"]["uncorrected_cell_ranger_counts_unavailable"] = True
    manifest[COHORT]["annotations"]["validation_not_done"] = False
    manifest[COHORT]["annotations"]["validated_field"] = verification["join"]["coarse_annotation_field"]
    manifest[COHORT]["analysis"] = {
        "completed": "2026-11-14",
        "publication_status": "preprint",
        "frozen_plan": "data/manifest/ovary_postmenopausal_analysis_plan.json",
        "script": "scripts/run_followup_ovary_postmenopausal.py",
        "report": "reports/followup-ovary-analysis.md",
        "inferential_unit": "donor",
        "effect": "exact_age_per_10_years_within_postmenopausal_cohort",
        "technical_replicates_combined": True,
        "x_identical_to_log1p_norm": True,
        "x_used_as_raw_counts": False,
        "uncorrected_sensitivity": "blocked_counts_unavailable",
        "Hellinger": "untestable_no_20_to_40_reference",
        "cross_cohort_pooling": False,
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
