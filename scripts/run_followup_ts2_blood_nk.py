#!/usr/bin/env python3
"""Run the frozen exploratory TS2 blood NK donor-level analysis.

Uses direct HDF5 group reads and streamed raw.X CSR slices. The seven-donor
family, model, P1 null, and endpoint floors were fixed before outcomes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
from statsmodels.stats.multitest import multipletests
from anndata.io import read_elem, sparse_dataset

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from scripts.prepare_ts2_blood_nk_gate import (  # noqa: E402
    EXPECTED_CANDIDATES, EXPECTED_DONORS, GIB, build_donor_eligibility,
    build_endpoint_status, donor_sort_key,
)
from src.aging_programs import load_program_tiers  # noqa: E402

RESULTS = REPO / "results"
DATA = REPO.parent / f"{REPO.name}-data" / "ts2"
H5AD = DATA / "blood.h5ad"
CACHE = DATA / "blood_nk_pseudobulk.npz"
MANIFEST = REPO / "data" / "manifest" / "tabula_sapiens_2_0_blood.json"
EXPECTED_SIZE = 2568875081
EXPECTED_SHA256 = "667680ecdb017d3b4f3393ad4b668ff82345a9297aa38d7ab2b2830df0d1aa52"
DONORS = sorted(EXPECTED_CANDIDATES, key=donor_sort_key)
SEED = 20261105
NULL_ITERATIONS = 200
CHUNK_CELLS = 128


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def nk_mask(obs):
    return (
        (obs["method"].astype(str) == "10X")
        & (obs["cell_type"].astype(str) == "natural killer cell")
    ).to_numpy()


def validate_raw_values(raw_group):
    data = raw_group["data"]
    bad = {"negative": 0, "noninteger": 0, "nonfinite": 0}
    for start in range(0, len(data), 5_000_000):
        values = data[start:start + 5_000_000]
        bad["negative"] += int((values < 0).sum())
        bad["noninteger"] += int((values != np.floor(values)).sum())
        bad["nonfinite"] += int((~np.isfinite(values)).sum())
    assert not any(bad.values()), f"raw.X count validation failed: {bad}"
    return {"n_stored_values": len(data), **bad}


def load_and_validate_metadata():
    assert H5AD.is_file()
    assert H5AD.stat().st_size == EXPECTED_SIZE
    with h5py.File(H5AD, "r") as handle:
        assert handle["X"].attrs["encoding-type"] == "csr_matrix"
        assert handle["raw/X"].attrs["encoding-type"] == "csr_matrix"
        assert tuple(handle["raw/X"].attrs["shape"]) == (85233, 60606)
        raw_check = validate_raw_values(handle["raw/X"])
        obs = read_elem(handle["obs"])
        var = read_elem(handle["raw/var"])
    assert obs.shape == (85233, 51)
    assert var.shape[0] == 60606
    assert obs.index.is_unique and obs["observation_joinid"].is_unique
    assert var.index.is_unique
    rows = np.flatnonzero(nk_mask(obs))
    nk = obs.iloc[rows].copy()
    nk["donor"] = nk["donor_id"].astype(str)
    nk["age"] = nk["development_stage"].astype(str).str.extract(r"(\d+)-year-old").astype(int)
    assert len(nk) == 2985 and set(nk["donor"]) == EXPECTED_DONORS
    for column, expected in {
        "tissue": "blood", "tissue_in_publication": "Blood", "compartment": "Immune",
        "disease": "normal", "method": "10X", "suspension_type": "cell",
    }.items():
        assert set(nk[column].astype(str)) == {expected}
    assert not nk["is_primary_data"].any()

    coverage = pd.read_csv(RESULTS / "followup-ts2-celltype-coverage.csv")
    expected = coverage[
        coverage["tissue"].eq("Blood") & coverage["compartment"].eq("Immune")
        & coverage["cell_type"].eq("natural killer cell") & coverage["method"].eq("10X")
    ].set_index("donor_id")
    observed = nk.groupby("donor", observed=True).size()
    observed.index = observed.index.astype(str)
    pd.testing.assert_series_equal(
        observed.sort_index(), expected["n_cells_supplementary"].sort_index(),
        check_names=False, check_index_type=False,
    )
    for donor, row in expected.iterrows():
        cells = nk[nk["donor"] == donor]
        assert set(cells["age"]) == {int(row["age_years"])}
        assert set(cells["sex"].astype(str).str.upper().str[0]) == {row["sex"]}
    return obs, var, rows, nk, expected, raw_check


def build_pseudobulk(rows, nk, var, rebuild=False):
    if CACHE.exists() and not rebuild:
        cached = np.load(CACHE, allow_pickle=False)
        assert cached["donors"].tolist() == DONORS
        assert cached["gene_ids"].tolist() == var.index.astype(str).tolist()
        return {key: cached[key] for key in cached.files}
    eligible = np.flatnonzero(nk["donor"].isin(DONORS).to_numpy())
    source_rows = rows[eligible]
    donor_values = nk.iloc[eligible]["donor"].to_numpy()
    n_genes = var.shape[0]
    sums = np.zeros((len(DONORS), n_genes), dtype=np.float64)
    depths = np.zeros(len(eligible), dtype=np.float64)
    donor_pos = {donor: index for index, donor in enumerate(DONORS)}
    with h5py.File(H5AD, "r") as handle:
        raw = sparse_dataset(handle["raw/X"])
        for start in range(0, len(eligible), CHUNK_CELLS):
            stop = min(start + CHUNK_CELLS, len(eligible))
            matrix = sp.csr_matrix(raw[source_rows[start:stop]], dtype=np.float64)
            library = np.asarray(matrix.sum(axis=1)).ravel()
            assert (library > 0).all()
            depths[start:stop] = library
            scale = 1e6 / library
            matrix.data *= np.repeat(scale, np.diff(matrix.indptr))
            np.log1p(matrix.data, out=matrix.data)
            for donor in set(donor_values[start:stop]):
                local = donor_values[start:stop] == donor
                sums[donor_pos[donor]] += np.asarray(matrix[local].sum(axis=0)).ravel()
    counts = np.array([(donor_values == donor).sum() for donor in DONORS])
    pseudobulk = sums / counts[:, None]
    donor_depth = {
        donor: float(np.log(depths[donor_values == donor]).mean()) for donor in DONORS
    }
    result = {
        "donors": np.asarray(DONORS),
        "gene_ids": var.index.astype(str).to_numpy(dtype=str),
        "gene_names": var["feature_name"].astype(str).to_numpy(dtype=str),
        "pseudobulk": pseudobulk,
        "n_cells": counts,
        "mean_log_umi": np.asarray([donor_depth[donor] for donor in DONORS]),
        "minimum_cell_umi": np.asarray([
            depths[donor_values == donor].min() for donor in DONORS
        ]),
        "maximum_cell_umi": np.asarray([
            depths[donor_values == donor].max() for donor in DONORS
        ]),
    }
    np.savez_compressed(CACHE, **result)
    return result


def donor_frame(nk, expression):
    rows = []
    for index, donor in enumerate(DONORS):
        cells = nk[nk["donor"] == donor]
        assays = set(cells["assay"].astype(str))
        rows.append({
            "donor_id": donor,
            "age_years": int(cells["age"].iloc[0]),
            "sex": cells["sex"].astype(str).iloc[0],
            "male": int(cells["sex"].astype(str).iloc[0] == "male"),
            "assay_labels": ";".join(sorted(assays)),
            "assay_5prime_present": int("10x 5' v2" in assays),
            "n_cells": int(expression["n_cells"][index]),
            "mean_log_umi": float(expression["mean_log_umi"][index]),
            "minimum_cell_umi": float(expression["minimum_cell_umi"][index]),
            "maximum_cell_umi": float(expression["maximum_cell_umi"][index]),
            "disease": "normal",
            "is_primary_data_all_false": True,
            "within_object_observation_ids_unique": True,
        })
    frame = pd.DataFrame(rows)
    assert frame["n_cells"].sum() == 2960
    return frame


def design_matrix(donors):
    matrix = np.column_stack([
        np.ones(len(donors)), donors["age_years"].to_numpy() / 10,
        donors["male"].to_numpy(), donors["assay_5prime_present"].to_numpy(),
        donors["mean_log_umi"].to_numpy(),
    ])
    assert np.linalg.matrix_rank(matrix) == matrix.shape[1]
    assert len(matrix) - matrix.shape[1] > 0
    return matrix


def fit_age_models(outcomes, design):
    """Fit all genes simultaneously; outcomes is donors x genes."""
    outcomes = np.asarray(outcomes, dtype=np.float64)
    inverse = np.linalg.inv(design.T @ design)
    coefficients = inverse @ design.T @ outcomes
    residual = outcomes - design @ coefficients
    residual_df = len(design) - design.shape[1]
    residual_ss = np.sum(residual ** 2, axis=0)
    variance = residual_ss / residual_df
    effect = coefficients[1]
    se = np.sqrt(variance * inverse[1, 1])
    with np.errstate(divide="ignore", invalid="ignore"):
        statistic = effect / se
    pvalue = 2 * stats.t.sf(np.abs(statistic), residual_df)
    critical = stats.t.ppf(0.975, residual_df)
    total_ss = np.sum((outcomes - outcomes.mean(axis=0)) ** 2, axis=0)
    return {
        "effect_age_10y": effect,
        "se": se,
        "ci95_low": effect - critical * se,
        "ci95_high": effect + critical * se,
        "pvalue": pvalue,
        "residual_df": np.full(outcomes.shape[1], residual_df),
        "r2": 1 - residual_ss / total_ss,
        "residual_ss": residual_ss,
    }


def prepare_p1(expression):
    program = load_program_tiers(REPO, "human")["P1_Clean"]
    assert len(program["up"]) == 76 and len(program["down"]) == 42
    names = expression["gene_names"].astype(str)
    lookup = {gene: index for index, gene in enumerate(names)}
    coverage_rows = []
    detected = {}
    matrix = expression["pseudobulk"].T
    for direction in ("up", "down"):
        indices = []
        for gene in program[direction]:
            present = gene in lookup
            index = lookup.get(gene, -1)
            total = float(matrix[index].sum()) if present else 0.0
            if present:
                indices.append(index)
            coverage_rows.append({
                "direction": direction, "gene": gene, "present_in_raw_var": present,
                "pseudobulk_sum_across_donors": total,
            })
        detected[direction] = np.asarray(indices, dtype=int)
    means = matrix.mean(axis=1)
    deciles = pd.qcut(means, 10, labels=False, duplicates="drop")
    deciles = pd.Series(deciles).fillna(-1).to_numpy()
    available = np.ones(len(names), dtype=bool)
    available[np.concatenate([detected["up"], detected["down"]])] = False
    pool_requirements = []
    valid = True
    for direction in ("up", "down"):
        valid &= len(detected[direction]) > 0
        for decile, count in pd.Series(deciles[detected[direction]]).value_counts().items():
            pool = int(((deciles == decile) & available).sum())
            valid &= pool >= int(count)
            pool_requirements.append({
                "direction": direction, "decile": int(decile),
                "signature_genes_needed": int(count), "available_pool_genes": pool,
                "passes": pool >= int(count),
            })
    return {
        "program": program, "detected": detected, "deciles": deciles,
        "available": available, "valid": bool(valid),
        "coverage": pd.DataFrame(coverage_rows),
        "pools": pd.DataFrame(pool_requirements),
    }


def p1_scores(matrix, indices):
    return np.asarray(matrix[indices].mean(axis=0)).ravel()


def run_p1(expression, design, setup):
    matrix = expression["pseudobulk"].T
    up = setup["detected"]["up"]
    down = setup["detected"]["down"]
    observed_donor = p1_scores(matrix, up) - p1_scores(matrix, down)
    observed = fit_age_models(observed_donor[:, None], design)
    rng = np.random.default_rng(SEED)
    null_scores = []
    null_effects = []
    deciles = setup["deciles"]
    available = setup["available"]
    for iteration in range(NULL_ITERATIONS):
        selected = {}
        for direction, genes in (("up", up), ("down", down)):
            chosen = []
            for decile, count in pd.Series(deciles[genes]).value_counts().items():
                pool = np.flatnonzero((deciles == decile) & available)
                chosen.extend(rng.choice(pool, size=int(count), replace=False))
            selected[direction] = np.asarray(chosen, dtype=int)
        score = p1_scores(matrix, selected["up"]) - p1_scores(matrix, selected["down"])
        null_scores.append(score)
        null_effects.append(float(fit_age_models(score[:, None], design)["effect_age_10y"][0]))
    null_effects = np.asarray(null_effects)
    observed_effect = float(observed["effect_age_10y"][0])
    centered_null = null_effects - null_effects.mean()
    empirical_p = (
        1 + int((np.abs(centered_null) >= abs(observed_effect - null_effects.mean())).sum())
    ) / (NULL_ITERATIONS + 1)
    model = {key: float(value[0]) for key, value in observed.items()}
    model.update({
        "endpoint": "P1_transport", "n_donors": len(design),
        "n_up_genes": len(up), "n_down_genes": len(down),
        "null_iterations": NULL_ITERATIONS, "null_seed": SEED,
        "null_mean_effect_age_10y": float(null_effects.mean()),
        "null_sd_effect_age_10y": float(null_effects.std(ddof=1)),
        "empirical_pvalue": empirical_p, "q_bh_single_endpoint": empirical_p,
    })
    null = pd.DataFrame({
        "iteration": np.arange(1, NULL_ITERATIONS + 1),
        "effect_age_10y": null_effects,
        "centered_effect_age_10y": centered_null,
        "donor_scores": list(np.asarray(null_scores)),
    })
    return observed_donor, model, null


def young_reference_floors(nk):
    eligible = nk[nk["donor"].isin(DONORS)]
    young = eligible[eligible["age"].between(20, 40)]
    young_donors = sorted(young["donor"].unique())
    remaining = {
        donor: int((young["donor"] != donor).sum()) for donor in DONORS
    }
    return {
        "young_donors": young_donors,
        "n_young_donors": len(young_donors),
        "young_cells": int(len(young)),
        "minimum_remaining_cells": min(remaining.values()),
        "passes": len(young_donors) >= 2 and min(remaining.values()) >= 50,
    }


def write_csv(frame, filename):
    path = RESULTS / filename
    frame.to_csv(path, index=False)
    print(f"wrote {path.relative_to(REPO)}: {len(frame)} rows", flush=True)


def update_static_outputs(asset_hash, donors, young, p1, statuses):
    asset_path = RESULTS / "followup-ts2-blood-nk-asset.csv"
    asset = pd.read_csv(asset_path)
    asset["maximum_file_gib"] = 3
    asset["within_file_cap"] = True
    asset["H5AD_verification"] = "verified_size_sha256_raw_counts_and_metadata"
    asset["local_H5AD"] = str(H5AD)
    asset["sha256"] = asset_hash
    write_csv(asset, asset_path.name)

    eligibility = build_donor_eligibility(pd.read_csv(RESULTS / "followup-ts2-celltype-coverage.csv"))
    eligibility["H5AD_metadata_verification"] = "verified_against_blood_H5AD"
    eligibility["raw_counts_verified"] = True
    eligibility["primary_data_deduplication_verified"] = True
    eligibility["is_primary_data_all_false"] = True
    eligibility["within_object_observation_ids_unique"] = True
    eligibility["deduplication_scope"] = "single_blood_object_no_cross_object_concatenation"
    eligibility["expression_eligible"] = eligibility["supplementary_candidate"]
    eligibility["status"] = "H5AD_verified_seven_donor_exploratory"
    eligibility.loc[~eligibility["supplementary_candidate"], "status"] = (
        "excluded_H5AD_verified_20_cell_floor"
    )
    eligibility["blocker"] = ""
    eligibility["eligible_expression_endpoints"] = "donor_pseudobulk_gene_models"
    if p1["valid"]:
        eligibility.loc[eligibility["supplementary_candidate"], "eligible_expression_endpoints"] += ";P1_transport"
    eligibility.loc[~eligibility["supplementary_candidate"], "eligible_expression_endpoints"] = ""
    write_csv(eligibility, "followup-ts2-blood-nk-donor-eligibility.csv")

    endpoints = build_endpoint_status()
    endpoints["donor_family_status"] = "seven_donor_exploratory_H5AD_verified"
    for index, row in endpoints.iterrows():
        status, reason, examined = statuses[row["endpoint"]]
        endpoints.loc[index, ["status", "reason", "outcomes_examined"]] = [status, reason, examined]
    write_csv(endpoints, "followup-ts2-blood-nk-endpoint-status.csv")

    manifest = json.loads(MANIFEST.read_text())
    manifest["access_gate"].pop("download_command_if_policy_changes", None)
    manifest["access_gate"].update({
        "maximum_file_gib": 3, "within_file_cap": True,
        "download_permitted": True, "local_H5AD": str(H5AD),
        "local_filesize_bytes": EXPECTED_SIZE, "sha256": asset_hash,
        "guarded_download_command": (
            "python scripts/resource_guard.py download --url "
            "https://datasets.cellxgene.cziscience.com/1bf55482-ad49-446a-87c6-f3442bebcef1.h5ad "
            "--output /home/reza/projects/agescape-revision-data/ts2/blood.h5ad "
            "--min-free-gb 30 --max-file-gb 3"
        ),
    })
    manifest["verification"].update({
        "raw_counts": True, "primary_data_and_deduplication": True,
        "per_cell_donor_age_sex_assay_disease": True,
        "P1_signature_floors": bool(p1["valid"]),
        "off_program_floors": bool(young["passes"]),
        "Hellinger_reference_floors": bool(young["passes"]),
        "is_primary_data_all_false": True,
        "deduplication_scope": "single_blood_object_unique_observation_and_join_ids",
    })
    manifest["donor_gate"]["candidate_ages_are_H5AD_verified"] = True
    manifest["expression_analysis"] = {
        "status": "completed_eligible_endpoints",
        "analysis_label": "exploratory",
        "design": "age_per_10y_plus_sex_assay5_presence_mean_log_raw_UMI_full_rank",
        "residual_df": 2,
        "endpoint_statuses": {key: value[0] for key, value in statuses.items()},
        "P1_null_iterations": NULL_ITERATIONS,
        "P1_null_seed": SEED,
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rebuild-pseudobulk", action="store_true")
    args = parser.parse_args()
    asset_hash = sha256(H5AD)
    assert asset_hash == EXPECTED_SHA256
    obs, var, rows, nk, expected, raw_check = load_and_validate_metadata()
    expression = build_pseudobulk(rows, nk, var, rebuild=args.rebuild_pseudobulk)
    donors = donor_frame(nk, expression)
    donors["age_years"] = donors["age_years"].astype(int)
    design = design_matrix(donors)
    young = young_reference_floors(nk)
    p1 = prepare_p1(expression)
    print(f"P1 floors pass: {p1['valid']}; young reference floors pass: {young['passes']}", flush=True)

    gates = pd.DataFrame([
        ("file_size_sha256", True, EXPECTED_SIZE, "exact approved asset and SHA-256"),
        ("raw_X_integer_counts", True, raw_check["n_stored_values"], "all stored values finite nonnegative integers"),
        ("observation_and_join_ids_unique", True, len(obs), "no duplicate within the single blood object"),
        ("is_primary_data", True, "all_False", "verified; no cross-object concatenation or primary-only claim"),
        ("NK_donor_counts_age_sex", True, len(nk), "exact supplementary agreement for nine donors"),
        ("NK_disease_assay_suspension", True, "normal;10x_3prime_v3_and_5prime_v2;cell", "per-cell H5AD metadata"),
        ("seven_donor_20_cell_floor", True, len(DONORS), "2,960 eligible cells; TSP8/TSP27 excluded"),
        ("continuous_age_design", True, "rank_5_of_5_residual_df_2", "age/sex/assay5/depth identifiable"),
        ("P1_signature_and_null_pools", p1["valid"], f"{len(p1['detected']['up'])}_up;{len(p1['detected']['down'])}_down", "nonempty directions and sufficient matched pools"),
        ("young_reference", young["passes"], f"{young['n_young_donors']}_donors;{young['young_cells']}_cells;min_remaining_{young['minimum_remaining_cells']}", "at least two young donors and 50 remaining cells for every donor"),
    ], columns=["gate", "passes", "observed", "requirement_or_scope"])
    write_csv(gates, "followup-ts2-blood-nk-h5ad-gates.csv")
    write_csv(donors, "followup-ts2-blood-nk-donor-summary.csv")
    write_csv(p1["coverage"], "followup-ts2-blood-nk-p1-gene-coverage.csv")
    write_csv(p1["pools"], "followup-ts2-blood-nk-p1-null-pools.csv")

    statuses = {
        "donor_pseudobulk_gene_models": ("completed", "Full-rank seven-donor exploratory model; BH across successfully fitted variable genes.", True),
        "P1_transport": (
            "completed" if p1["valid"] else "untestable",
            "Mapped P1 directions and expression-matched null pools passed." if p1["valid"] else "P1 signature/null-pool floor failed.",
            bool(p1["valid"]),
        ),
        "off_program_expression": (
            "untestable",
            f"Only {young['n_young_donors']} young donor (20-40 years); cross-fitted reference has zero cells for TSP10. No off-program markers or outcomes computed.",
            False,
        ),
        "Hellinger_displacement": (
            "untestable",
            f"Established human rule requires at least two young donors; observed {young['n_young_donors']}. No Hellinger distances computed.",
            False,
        ),
    }

    gene_counts = {"tested": 0, "q_lt_0.05": 0}
    if statuses["donor_pseudobulk_gene_models"][2]:
        matrix = expression["pseudobulk"]
        variable = np.ptp(matrix, axis=0) > 0
        fitted = fit_age_models(matrix[:, variable], design)
        genes = pd.DataFrame({
            "ensembl_id": expression["gene_ids"][variable],
            "gene": expression["gene_names"][variable],
            **{key: value for key, value in fitted.items()},
        })
        genes["model"] = "OLS_age10_sex_assay5_mean_log_raw_UMI"
        genes["n_donors"] = len(DONORS)
        genes["q_bh"] = np.nan
        testable = np.isfinite(genes["pvalue"]) & np.isfinite(genes["se"]) & (genes["se"] > 0)
        genes["test_status"] = "tested"
        genes.loc[~testable, "test_status"] = "untestable_zero_standard_error"
        genes.loc[testable, "q_bh"] = multipletests(genes.loc[testable, "pvalue"], method="fdr_bh")[1]
        genes = genes.sort_values(["q_bh", "pvalue", "ensembl_id"], na_position="last")
        write_csv(genes, "followup-ts2-blood-nk-gene-models.csv")
        gene_counts = {"tested": int(testable.sum()), "q_lt_0.05": int((genes["q_bh"] < 0.05).sum())}

    if p1["valid"]:
        donor_scores, p1_model, null = run_p1(expression, design, p1)
        donors["P1_score"] = donor_scores
        write_csv(donors, "followup-ts2-blood-nk-donor-summary.csv")
        write_csv(pd.DataFrame([p1_model]), "followup-ts2-blood-nk-p1-model.csv")
        null.drop(columns="donor_scores").to_csv(RESULTS / "followup-ts2-blood-nk-p1-null.csv", index=False)
        np.savez_compressed(DATA / "blood_nk_p1_null_donor_scores.npz", scores=np.asarray(null["donor_scores"].tolist()))
        print(f"P1 effect per 10y {p1_model['effect_age_10y']:+.8g}; empirical p={p1_model['empirical_pvalue']:.6g}", flush=True)
    print(f"Gene models tested {gene_counts['tested']}; BH q<0.05 {gene_counts['q_lt_0.05']}", flush=True)

    update_static_outputs(asset_hash, donors, young, p1, statuses)
    manifest = json.loads(MANIFEST.read_text())
    manifest["expression_analysis"]["gene_models_tested"] = gene_counts["tested"]
    manifest["expression_analysis"]["gene_models_q_lt_0_05"] = gene_counts["q_lt_0.05"]
    MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print("Completed all eligible endpoints; Hellinger and off-program remain untestable.", flush=True)


if __name__ == "__main__":
    main()
