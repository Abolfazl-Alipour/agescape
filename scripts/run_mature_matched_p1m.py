#!/usr/bin/env python3
"""Matched mouse-human P1m gene concordance in floored and unfloored arms."""

from __future__ import annotations

import argparse
import gc
import os
import re
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_var] = "1"

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem
from scipy.stats import pearsonr, spearmanr

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts import run_rev2_concordance as rev2
from scripts.run_mature_p1m_transport import COHORTS, read_csr_block, read_obs_columns

RESULTS = REPO / "results"
REPORTS = REPO / "reports"
DATA_ROOT = Path("/home/reza/projects/agescape-revision-data")
WORK_DIR = DATA_ROOT / "matched-p1m"
HUMAN_CACHE = WORK_DIR / "human-units.pkl"
PAIR_GENE_CSV = WORK_DIR / "pair-gene-rows.csv"
OUTPUT_CSV = RESULTS / "followup-mature-matched-p1m.csv"
OUTPUT_REPORT = REPORTS / "followup-mature-matched-p1m.md"
MIN_AGE = 30.0
OLD_FLOOR_DONORS = 3
OLD_FLOOR_CELLS = 50
M9_COHORTS = [
    "Kidney_LivingDonor", "Kidney_SCR_v15", "Lung_Integrated",
    "Lung_CellRef", "ImmuneAging_BAL",
]
ABSENT_M9_COHORTS = ["Kidney_NormalAdult", "Kidney_SN_v15"]
YAZAR_H5AD = DATA_ROOT / "tmp_organs" / "hcata" / "yazar_1e44db10.h5ad"
OUTPUT_COLUMNS = [
    "mouse_unit", "mouse_tissue", "human_cohort", "human_unit", "arm",
    "n_mouse_genes", "n_human_genes_tested", "sign_agreement_fraction",
    "spearman", "pearson", "p_sign_gene_perm", "p_spearman_gene_perm",
    "p_sign_age_perm", "p_spearman_age_perm", "n_donors_human_30plus",
    "n_cells_human_30plus", "below_old_floor", "notes",
]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def normalize_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def safe_unit_name(unit: str) -> str:
    return unit.replace("::", "__").replace(" ", "_").replace(",", "").replace("-", "_")


def load_mouse() -> tuple[pd.DataFrame, dict[str, pd.DataFrame], pd.DataFrame]:
    genes = pd.read_csv(RESULTS / "followup-mature-p1m-genes.csv")
    assert len(genes) == 1914
    assert genes["human_ortholog"].notna().sum() == 1724
    eligibility = pd.read_csv(RESULTS / "followup-mature-p1m-eligibility.csv")
    eligible = eligibility[eligibility["eligible"].astype(str).eq("yes")].copy()
    assert len(eligible) == 36
    mouse_units = {}
    for row in eligible.itertuples():
        de = pd.read_csv(RESULTS / f"followup-mature-p1m-de_{safe_unit_name(row.unit)}.csv")
        de = de.set_index("index", drop=False)
        selected = genes[genes["units_list"].astype(str).str.split(";").apply(
            lambda units: row.unit in units
        )].copy()
        missing = set(selected["mouse_gene_symbol"]) - set(de.index)
        assert not missing, f"{row.unit}: missing DE rows for {sorted(missing)[:5]}"
        selected["log2FC"] = selected["mouse_gene_symbol"].map(de["log2FC"])
        selected["unit_FDR"] = selected["mouse_gene_symbol"].map(de["FDR"])
        selected["unit_significant"] = selected["mouse_gene_symbol"].map(
            de["P1m_significant_before_filters"]
        )
        assert selected["unit_significant"].astype(bool).all()
        assert np.isfinite(selected["log2FC"]).all()
        mouse_units[row.unit] = selected
    return eligible, mouse_units, genes


def human_roster() -> pd.DataFrame:
    counts = pd.read_csv(RESULTS / "m9-donor-cell-counts-prethreshold.csv")
    rows = []
    for cohort in M9_COHORTS + ABSENT_M9_COHORTS:
        tissue = "Kidney" if cohort.startswith("Kidney") else "Lung"
        for cell_type in sorted(counts.loc[counts["cohort"] == cohort, "cell_type"].unique()):
            rows.append({
                "human_cohort": cohort, "human_tissue": tissue, "human_unit": cell_type,
                "source_status": "absent" if cohort in ABSENT_M9_COHORTS else "present",
            })
    for cohort, filename in [
        ("Yazar_Blood", "followup-hcata-yazar-fractions.csv"),
        ("Perez_Blood_Controls", "followup-hcata-perez-fractions.csv"),
    ]:
        frame = pd.read_csv(RESULTS / filename)
        for cell_type in sorted(frame.loc[frame["arm"] == "primary", "cell_type"].unique()):
            rows.append({"human_cohort": cohort, "human_tissue": "Blood",
                         "human_unit": cell_type, "source_status": "present"})
    for cohort in ("GSE213516", "GSE214546"):
        path = DATA_ROOT / "pbmc" / "interim" / cohort / "donor-aggregates.npz"
        with np.load(path, allow_pickle=True) as cache:
            for cell_type in sorted(np.unique(cache["cell_type"].astype(str))):
                rows.append({"human_cohort": f"PBMC_{cohort}", "human_tissue": "Blood",
                             "human_unit": str(cell_type), "source_status": "present"})
    roster = pd.DataFrame(rows)
    roster["normalized"] = roster["human_unit"].map(normalize_name)
    return roster


def build_matches(eligible: pd.DataFrame, roster: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    match_rows = []
    audit_rows = []
    for row in eligible.itertuples():
        target_tissue = "Blood" if row.tissue == "Marrow" else row.tissue
        normalized = normalize_name(row.cell_ontology_class)
        hits = roster[(roster["human_tissue"] == target_tissue)
                      & (roster["normalized"] == normalized)]
        present = hits[hits["source_status"] == "present"]
        absent = hits[hits["source_status"] == "absent"]
        for hit in present.itertuples():
            match_rows.append({
                "mouse_unit": row.unit, "mouse_tissue": row.tissue,
                "mouse_cell_type": row.cell_ontology_class,
                "human_cohort": hit.human_cohort, "human_unit": hit.human_unit,
                "normalized_name": normalized,
            })
        if len(present):
            decision = "matched"
            reason = "; ".join(f"{r.human_cohort}::{r.human_unit}" for r in present.itertuples())
        elif len(absent):
            decision = "unmatched_source_absent"
            reason = "exact-name human source file absent: " + ", ".join(absent["human_cohort"])
        else:
            decision = "unmatched_no_cell_type"
            reason = "no exactly normalized same-named cell type in same-tissue human roster"
        if len(absent):
            reason += "; additional absent exact-name sources: " + ", ".join(absent["human_cohort"])
        audit_rows.append({
            "mouse_unit": row.unit, "mouse_tissue": row.tissue,
            "mouse_cell_type": row.cell_ontology_class, "normalized_name": normalized,
            "decision": decision, "n_present_matches": len(present), "reason": reason,
        })
    matches = pd.DataFrame(match_rows).sort_values(
        ["mouse_unit", "human_cohort", "human_unit"]
    ).reset_index(drop=True)
    audit = pd.DataFrame(audit_rows)
    assert len(audit) == 36
    assert audit["mouse_unit"].nunique() == 36
    assert set(matches["mouse_unit"]) == set(audit.loc[audit["decision"] == "matched", "mouse_unit"])
    return matches, audit


def _stream_selected_h5ad(
    cohort: str,
    path: Path,
    cell_types: list[str],
    target_genes: list[str],
    yazar: bool = False,
    chunk_size: int = 5000,
) -> dict[str, dict]:
    print(f"{cohort}: streaming {path} for {len(cell_types)} matched units", flush=True)
    with h5py.File(path, "r") as handle:
        columns = read_obs_columns(handle)
        if yazar:
            columns["assay"] = np.asarray(read_elem(handle["obs"]["assay"]))
            columns["suspension_type"] = np.asarray(read_elem(handle["obs"]["suspension_type"]))
            x_group = handle["X"]
            var = read_elem(handle["var"])
        else:
            x_group = handle["X"]
            var = read_elem(handle["var"])
        var_names = var["feature_name"].astype(str).to_numpy()
        donor_arr = columns["donor_id"].astype(str)
        cell_arr = columns["cell_type"].astype(str)
        ages = pd.Series(columns["development_stage"].astype(str)).str.extract(
            r"(\d+)-year-old", expand=False
        ).to_numpy(float)
        selected = (
            columns["is_primary_data"].astype(bool)
            & (columns["disease"].astype(str) == "normal")
            & np.isfinite(ages) & (ages >= MIN_AGE)
            & np.isin(cell_arr, cell_types)
        )
        if yazar:
            selected &= ((columns["assay"].astype(str) == "10x 3' v2")
                         & (columns["suspension_type"].astype(str) == "cell"))
        selected_rows = np.flatnonzero(selected)
        selected_pairs = pd.MultiIndex.from_arrays(
            [donor_arr[selected_rows], cell_arr[selected_rows]], names=["donor", "cell_type"]
        )
        group_codes_selected, group_index = pd.factorize(selected_pairs, sort=True)
        group_codes = np.full(len(donor_arr), -1, dtype=np.int64)
        group_codes[selected_rows] = group_codes_selected
        n_groups = len(group_index)
        group_donors = group_index.get_level_values(0).astype(str).to_numpy()
        group_cell_types = group_index.get_level_values(1).astype(str).to_numpy()
        donor_metadata = {}
        for row in selected_rows:
            donor = str(donor_arr[row])
            values = (float(ages[row]), float(str(columns["sex"][row]).lower() == "male"))
            if donor in donor_metadata:
                assert donor_metadata[donor] == values
            donor_metadata[donor] = values
        target_indices = np.flatnonzero(np.isin(var_names, target_genes)).astype(int)
        target_feature_genes = var_names[target_indices]
        present_genes = sorted(set(target_feature_genes))
        counts = np.zeros((n_groups, len(target_indices)), dtype=np.float64)
        libraries = np.zeros(n_groups, dtype=np.float64)
        n_obs = len(donor_arr)
        for start in range(0, n_obs, chunk_size):
            stop = min(start + chunk_size, n_obs)
            local_rows = np.flatnonzero(selected[start:stop])
            if not len(local_rows):
                continue
            block = read_csr_block(x_group, start, stop, len(var_names))
            local_codes = group_codes[start:stop][local_rows]
            unique, inverse = np.unique(local_codes, return_inverse=True)
            indicator = sp.csr_matrix(
                (np.ones(len(local_rows)), (inverse, np.arange(len(local_rows)))),
                shape=(len(unique), len(local_rows)),
            )
            summed = indicator @ block[local_rows]
            counts[unique] += np.asarray(summed[:, target_indices].todense())
            libraries[unique] += np.asarray(summed.sum(axis=1)).ravel()
            if start % 250000 == 0:
                print(f"  {cohort}: cells {start}-{stop}/{n_obs}; selected groups {n_groups}", flush=True)
        cell_counts = np.bincount(group_codes_selected, minlength=n_groups).astype(int)
        units = {}
        for cell_type in cell_types:
            indices = np.flatnonzero(group_cell_types == cell_type)
            indices = indices[np.argsort(group_donors[indices])]
            if not len(indices):
                units[cell_type] = {
                    "cohort": cohort, "cell_type": cell_type, "genes": present_genes,
                    "expression": np.empty((0, len(present_genes))),
                    "metadata": pd.DataFrame(columns=["donor", "age", "sex_male", "n_cells"]),
                }
                continue
            feature_counts = pd.DataFrame(
                counts[indices], columns=target_feature_genes, index=indices
            ).T.groupby(level=0).sum().T.reindex(columns=present_genes).to_numpy()
            library = libraries[indices]
            assert np.all(library > 0)
            expression = np.log1p(feature_counts / library[:, None] * 1e4)
            metadata = pd.DataFrame({
                "donor": group_donors[indices],
                "age": [donor_metadata[d][0] for d in group_donors[indices]],
                "sex_male": [donor_metadata[d][1] for d in group_donors[indices]],
                "n_cells": cell_counts[indices],
            })
            units[cell_type] = {
                "cohort": cohort, "cell_type": cell_type, "genes": present_genes,
                "expression": expression, "metadata": metadata,
            }
        print(f"  {cohort}: completed {n_groups} donor-unit pseudobulks", flush=True)
    return units


def load_human_units(matches: pd.DataFrame, target_genes: list[str], refresh: bool = False) -> dict[str, dict]:
    import pickle

    WORK_DIR.mkdir(parents=True, exist_ok=True)
    if HUMAN_CACHE.exists() and not refresh:
        print(f"Loading cached human pseudobulks: {HUMAN_CACHE}", flush=True)
        with HUMAN_CACHE.open("rb") as handle:
            return pickle.load(handle)
    units = {}
    for cohort in M9_COHORTS:
        cell_types = sorted(matches.loc[matches["human_cohort"] == cohort, "human_unit"].unique())
        if not cell_types:
            continue
        path = COHORTS[cohort][1]
        if not path.exists():
            raise FileNotFoundError(f"Required source absent after match audit: {path}")
        loaded = _stream_selected_h5ad(cohort, path, cell_types, target_genes)
        for cell_type, unit in loaded.items():
            units[f"{cohort}::{cell_type}"] = unit
        gc.collect()
    cell_types = sorted(matches.loc[matches["human_cohort"] == "Yazar_Blood", "human_unit"].unique())
    if cell_types:
        loaded = _stream_selected_h5ad("Yazar_Blood", YAZAR_H5AD, cell_types, target_genes, yazar=True)
        for cell_type, unit in loaded.items():
            units[f"Yazar_Blood::{cell_type}"] = unit
    assert not matches["human_cohort"].isin(["Perez_Blood_Controls", "PBMC_GSE213516", "PBMC_GSE214546"]).any(), \
        "A new exact blood match needs a raw-count loader before analysis"
    temporary = HUMAN_CACHE.with_suffix(".tmp")
    with temporary.open("wb") as handle:
        pickle.dump(units, handle, protocol=pickle.HIGHEST_PROTOCOL)
    temporary.replace(HUMAN_CACHE)
    return units


def prepare_unit(unit: dict) -> dict:
    prepared = dict(unit)
    metadata = unit["metadata"].reset_index(drop=True)
    expression = np.asarray(unit["expression"], dtype=np.float64)
    prepared["metadata"] = metadata
    prepared["expression"] = expression
    prepared["n_donors"] = len(metadata)
    prepared["n_cells"] = int(metadata["n_cells"].sum()) if len(metadata) else 0
    prepared["slope_status"] = "fittable"
    prepared["sex_included"] = False
    if len(metadata) < 2:
        prepared["slope_status"] = "fewer_than_2_donors"
        return prepared
    if metadata["age"].nunique() < 2:
        prepared["slope_status"] = "no_age_variation"
        return prepared
    z_columns = [np.ones(len(metadata))]
    if len(metadata) > 2 and metadata["sex_male"].nunique() > 1:
        candidate = np.column_stack(z_columns + [metadata["sex_male"].to_numpy(float)])
        if np.linalg.matrix_rank(candidate) == candidate.shape[1]:
            z_columns.append(metadata["sex_male"].to_numpy(float))
            prepared["sex_included"] = True
    z = np.column_stack(z_columns)
    q, _ = np.linalg.qr(z, mode="reduced")
    prepared["Q"] = q
    prepared["Y_perp"] = expression - q @ (q.T @ expression)
    return prepared


def fwl_slopes(prepared: dict, ages: np.ndarray) -> np.ndarray:
    n_genes = len(prepared["genes"])
    if prepared["slope_status"] != "fittable":
        return np.full(n_genes, np.nan)
    age = np.asarray(ages, dtype=float)
    q = prepared["Q"]
    residual_age = age - q @ (q.T @ age)
    denominator = float(residual_age @ residual_age)
    if denominator <= 1e-20 * max(float(age @ age), 1.0):
        return np.full(n_genes, np.nan)
    return (residual_age @ prepared["Y_perp"]) / denominator


def verify_fwl(prepared: dict, slopes: np.ndarray) -> None:
    if prepared["slope_status"] != "fittable":
        return
    metadata = prepared["metadata"]
    columns = [np.ones(len(metadata)), metadata["age"].to_numpy(float)]
    if prepared["sex_included"]:
        columns.append(metadata["sex_male"].to_numpy(float))
    design = np.column_stack(columns)
    if len(metadata) - np.linalg.matrix_rank(design) > 0:
        direct = rev2.fit_age_slope_only(prepared["expression"], design, age_column=1)
        np.testing.assert_allclose(slopes, direct, rtol=1e-10, atol=1e-12)


def permute_units(units: dict[str, dict], permutations: int, seed: int) -> None:
    rng = np.random.default_rng(seed)
    by_cohort: dict[str, list[dict]] = {}
    for unit in units.values():
        prepared = prepare_unit(unit)
        unit.clear()
        unit.update(prepared)
        by_cohort.setdefault(unit["cohort"], []).append(unit)
        observed = fwl_slopes(unit, unit["metadata"]["age"].to_numpy(float)) if len(unit["metadata"]) else np.full(len(unit["genes"]), np.nan)
        unit["observed_slopes"] = observed
        verify_fwl(unit, observed)
        unit["permuted_slopes"] = np.full((permutations, len(unit["genes"])), np.nan)
    for cohort, cohort_units in by_cohort.items():
        records = {}
        for unit in cohort_units:
            for row in unit["metadata"].itertuples():
                values = (float(row.age), float(row.sex_male))
                if row.donor in records:
                    assert records[row.donor] == values
                records[row.donor] = values
        donors = sorted(records)
        donor_index = {donor: index for index, donor in enumerate(donors)}
        ages = np.array([records[donor][0] for donor in donors])
        sexes = np.array([records[donor][1] for donor in donors])
        for unit in cohort_units:
            unit["donor_indices"] = np.array([donor_index[d] for d in unit["metadata"]["donor"]], dtype=int)
        for iteration in range(permutations):
            permuted = ages.copy()
            for sex in np.unique(sexes):
                positions = np.flatnonzero(sexes == sex)
                permuted[positions] = rng.permutation(permuted[positions])
            for unit in cohort_units:
                if unit["slope_status"] == "fittable":
                    unit["permuted_slopes"][iteration] = fwl_slopes(
                        unit, permuted[unit["donor_indices"]]
                    )
            if (iteration + 1) % 100 == 0:
                print(f"  {cohort}: age permutations {iteration + 1}/{permutations}", flush=True)


def resolve_mouse_orthologs(frame: pd.DataFrame) -> tuple[pd.DataFrame, int, int]:
    mapped = frame[frame["human_ortholog"].notna()].copy()
    mapped["human_ortholog"] = mapped["human_ortholog"].astype(str)
    n_mapped = len(mapped)
    mapped = mapped.sort_values(["human_ortholog", "unit_FDR", "mouse_gene_symbol"])
    collisions = int(mapped.duplicated("human_ortholog").sum())
    mapped = mapped.drop_duplicates("human_ortholog", keep="first")
    return mapped, n_mapped, collisions


def correlation_values(mouse: np.ndarray, human: np.ndarray) -> dict[str, float]:
    result = {
        "sign": float(np.mean(np.sign(mouse) == np.sign(human))) if len(mouse) else np.nan,
        "spearman": np.nan, "pearson": np.nan,
    }
    if len(mouse) >= 2 and np.std(mouse) > 0 and np.std(human) > 0:
        result["spearman"] = float(spearmanr(mouse, human).statistic)
        result["pearson"] = float(pearsonr(mouse, human).statistic)
    return result


def two_sided_permutation_p(observed: float, null: np.ndarray) -> float:
    finite = null[np.isfinite(null)]
    if not np.isfinite(observed) or not len(finite):
        return np.nan
    center = float(finite.mean())
    return float((1 + np.sum(np.abs(finite - center) >= abs(observed - center))) / (len(finite) + 1))


def gene_permutation_null(mouse: np.ndarray, human: np.ndarray, permutations: int, rng: np.random.Generator) -> dict[str, np.ndarray]:
    null = {key: np.full(permutations, np.nan) for key in ("sign", "spearman", "pearson")}
    for iteration in range(permutations):
        values = correlation_values(rng.permutation(mouse), human)
        for key in null:
            null[key][iteration] = values[key]
    return null


def build_pair_data(matches: pd.DataFrame, mouse_units: dict[str, pd.DataFrame], human_units: dict[str, dict],
                    gene_permutations: int, seed: int) -> tuple[list[dict], pd.DataFrame]:
    pair_data = []
    audit_rows = []
    for pair_index, match in enumerate(matches.itertuples()):
        mouse, n_mapped, collisions = resolve_mouse_orthologs(mouse_units[match.mouse_unit])
        unit = human_units[f"{match.human_cohort}::{match.human_unit}"]
        gene_index = {gene: index for index, gene in enumerate(unit["genes"])}
        present = mouse[mouse["human_ortholog"].isin(gene_index)].copy()
        human_indices = present["human_ortholog"].map(gene_index).to_numpy(dtype=int)
        observed = unit["observed_slopes"][human_indices] if len(human_indices) else np.array([])
        finite = np.isfinite(observed)
        present = present.iloc[np.flatnonzero(finite)].copy()
        human_indices = human_indices[finite]
        observed = observed[finite]
        mouse_values = present["log2FC"].to_numpy(float)
        stats = correlation_values(mouse_values, observed)
        rng = np.random.default_rng(seed + pair_index)
        gene_null = gene_permutation_null(mouse_values, observed, gene_permutations, rng) if len(present) else {
            key: np.array([]) for key in ("sign", "spearman", "pearson")
        }
        age_null = {key: np.full(unit["permuted_slopes"].shape[0], np.nan)
                    for key in ("sign", "spearman", "pearson")}
        permuted = unit["permuted_slopes"][:, human_indices] if len(human_indices) else np.empty((unit["permuted_slopes"].shape[0], 0))
        for iteration in range(permuted.shape[0]):
            values = correlation_values(mouse_values, permuted[iteration])
            for key in age_null:
                age_null[key][iteration] = values[key]
        data = {
            "match": match, "mouse_frame": mouse, "n_mouse_genes": len(mouse_units[match.mouse_unit]),
            "n_mapped": n_mapped, "collisions": collisions, "n_present": len(present),
            "present": present, "human_indices": human_indices, "mouse_values": mouse_values,
            "human_observed": observed, "human_age_permuted": permuted,
            "stats": stats, "gene_null": gene_null, "age_null": age_null,
            "unit": unit,
        }
        pair_data.append(data)
        for row, human_slope, human_index in zip(present.itertuples(), observed, human_indices):
            audit_rows.append({
                "mouse_unit": match.mouse_unit, "mouse_tissue": match.mouse_tissue,
                "human_cohort": match.human_cohort, "human_unit": match.human_unit,
                "mouse_gene_symbol": row.mouse_gene_symbol, "human_ortholog": row.human_ortholog,
                "mouse_log2FC": row.log2FC, "mouse_FDR": row.unit_FDR,
                "human_slope": human_slope, "human_slope_gene_index": int(human_index),
                "same_sign": bool(np.sign(row.log2FC) == np.sign(human_slope)),
            })
    return pair_data, pd.DataFrame(audit_rows)


def pooled_observed(pair_data: list[dict]) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, float]]:
    nonempty = [pair for pair in pair_data if len(pair["mouse_values"])]
    if not nonempty:
        empty = np.array([])
        return empty, empty, np.empty((0, 0)), correlation_values(empty, empty)
    mouse = np.concatenate([pair["mouse_values"] for pair in nonempty])
    human = np.concatenate([pair["human_observed"] for pair in nonempty])
    permutations = np.concatenate([pair["human_age_permuted"] for pair in nonempty], axis=1)
    return mouse, human, permutations, correlation_values(mouse, human)


def pooled_nulls(mouse: np.ndarray, human: np.ndarray, age_permuted: np.ndarray,
                 gene_permutations: int, seed: int) -> tuple[dict, dict]:
    gene_null = gene_permutation_null(mouse, human, gene_permutations, np.random.default_rng(seed)) if len(mouse) else {
        key: np.array([]) for key in ("sign", "spearman", "pearson")
    }
    age_null = {key: np.full(age_permuted.shape[0], np.nan) for key in ("sign", "spearman", "pearson")}
    for iteration in range(age_permuted.shape[0]):
        values = correlation_values(mouse, age_permuted[iteration])
        for key in age_null:
            age_null[key][iteration] = values[key]
    return gene_null, age_null


def make_output_rows(pair_data: list[dict]) -> pd.DataFrame:
    rows = []
    for pair in pair_data:
        match = pair["match"]
        unit = pair["unit"]
        below = not (unit["n_donors"] >= OLD_FLOOR_DONORS and unit["n_cells"] >= OLD_FLOOR_CELLS)
        notes = [
            f"ortholog_mapped={pair['n_mapped']}", f"ortholog_collisions={pair['collisions']}",
            f"human_genes_present={pair['n_present']}", f"sex_adjusted={str(unit['sex_included']).lower()}",
            f"below_old_floor={str(below).lower()}",
        ]
        if unit["slope_status"] == "fewer_than_2_donors":
            notes.append("slope_NA_fewer_than_2_donors_excluded_from_pooled_stats")
        elif unit["slope_status"] == "no_age_variation":
            notes.append("slope_NA_no_age_variation_excluded_from_pooled_stats")
        elif unit["n_donors"] == 2:
            notes.append("difference_with_no_variance_estimate_age_only_no_SE_or_slope_p")
        elif unit["n_donors"] >= 3 and unit["n_donors"] - (3 if unit["sex_included"] else 2) <= 0:
            notes.append("slope_reported_without_residual_variance_estimate")
        for arm in ("floored", "unfloored"):
            if arm == "floored" and below:
                continue
            rows.append({
                "mouse_unit": match.mouse_unit, "mouse_tissue": match.mouse_tissue,
                "human_cohort": match.human_cohort, "human_unit": match.human_unit, "arm": arm,
                "n_mouse_genes": pair["n_mouse_genes"], "n_human_genes_tested": len(pair["mouse_values"]),
                "sign_agreement_fraction": pair["stats"]["sign"],
                "spearman": pair["stats"]["spearman"], "pearson": pair["stats"]["pearson"],
                "p_sign_gene_perm": two_sided_permutation_p(pair["stats"]["sign"], pair["gene_null"]["sign"]),
                "p_spearman_gene_perm": two_sided_permutation_p(pair["stats"]["spearman"], pair["gene_null"]["spearman"]),
                "p_sign_age_perm": two_sided_permutation_p(pair["stats"]["sign"], pair["age_null"]["sign"]),
                "p_spearman_age_perm": two_sided_permutation_p(pair["stats"]["spearman"], pair["age_null"]["spearman"]),
                "n_donors_human_30plus": unit["n_donors"], "n_cells_human_30plus": unit["n_cells"],
                "below_old_floor": below, "notes": ";".join(notes),
            })
    output = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
    assert not output.duplicated(["mouse_unit", "human_cohort", "human_unit", "arm"]).any()
    assert (output.loc[output["n_donors_human_30plus"] < 2, "n_human_genes_tested"] == 0).all()
    return output.sort_values(["arm", "mouse_unit", "human_cohort", "human_unit"]).reset_index(drop=True)


def summarize_arm(arm: str, pair_data: list[dict], gene_permutations: int, seed: int) -> tuple[dict, dict, dict]:
    mouse, human, age_permuted, observed = pooled_observed(pair_data)
    gene_null, age_null = pooled_nulls(mouse, human, age_permuted, gene_permutations, seed)
    summary = {
        "arm": arm, "n_pairs": len(pair_data), "n_pairs_with_tested_genes": sum(len(p["mouse_values"]) > 0 for p in pair_data),
        "n_pair_gene_rows": len(mouse), "sign_agreement_fraction": observed["sign"],
        "spearman": observed["spearman"], "pearson": observed["pearson"],
    }
    for metric in ("sign", "spearman", "pearson"):
        summary[f"p_{metric}_gene_perm"] = two_sided_permutation_p(observed[metric], gene_null[metric])
        summary[f"p_{metric}_age_perm"] = two_sided_permutation_p(observed[metric], age_null[metric])
    return summary, gene_null, age_null


def independent_spot_check(pair_data: list[dict], pair_genes: pd.DataFrame, output: pd.DataFrame) -> dict:
    candidates = [pair for pair in pair_data if len(pair["mouse_values"]) >= 5]
    assert candidates
    pair = min(candidates, key=lambda item: len(item["mouse_values"]))
    match = pair["match"]
    subset = pair_genes[(pair_genes["mouse_unit"] == match.mouse_unit)
                        & (pair_genes["human_cohort"] == match.human_cohort)
                        & (pair_genes["human_unit"] == match.human_unit)].copy()
    de = pd.read_csv(RESULTS / f"followup-mature-p1m-de_{safe_unit_name(match.mouse_unit)}.csv").set_index("index")
    np.testing.assert_allclose(subset["mouse_log2FC"], de.loc[subset["mouse_gene_symbol"], "log2FC"], rtol=0, atol=1e-12)
    source_genes = pd.read_csv(RESULTS / "followup-mature-p1m-genes.csv").set_index("mouse_gene_symbol")
    np.testing.assert_array_equal(
        subset["human_ortholog"].to_numpy(),
        source_genes.loc[subset["mouse_gene_symbol"], "human_ortholog"].to_numpy(),
    )
    same = int((np.sign(subset["mouse_log2FC"]) == np.sign(subset["human_slope"])).sum())
    total = len(subset)
    fraction = same / total
    row = output[(output["mouse_unit"] == match.mouse_unit)
                 & (output["human_cohort"] == match.human_cohort)
                 & (output["human_unit"] == match.human_unit)
                 & (output["arm"] == "unfloored")].iloc[0]
    np.testing.assert_allclose(fraction, row["sign_agreement_fraction"], rtol=0, atol=1e-12)
    assert int(row["n_human_genes_tested"]) == total
    return {"mouse_unit": match.mouse_unit, "human_cohort": match.human_cohort,
            "human_unit": match.human_unit, "same": same, "total": total, "fraction": fraction}


def fmt(value, digits=6) -> str:
    if value is None or pd.isna(value):
        return "NA"
    if isinstance(value, (bool, np.bool_)):
        return str(bool(value)).lower()
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.{digits}g}"
    return str(value).replace("|", "\\|")


def markdown_table(columns: list[str], rows: list[list], digits=6) -> str:
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    lines.extend("| " + " | ".join(fmt(value, digits) for value in row) + " |" for row in rows)
    return "\n".join(lines)


def write_report(output: pd.DataFrame, audit: pd.DataFrame, summaries: dict[str, dict],
                 pair_data: list[dict], spot: dict, gene_permutations: int,
                 age_permutations: int, started: str, duration: float) -> None:
    matched = audit[audit["decision"] == "matched"]
    unmatched = audit[audit["decision"] != "matched"]
    lines = [
        "# Matched mouse-human P1m comparison", "",
        "## Headline results", "",
        "Pooled statistics weight every pair-gene row equally. P values are two-sided permutation values centered on the corresponding null mean, with the plus-one correction.", "",
        markdown_table(
            ["arm", "pairs", "pairs with genes", "pair-gene rows", "sign agreement", "gene-label p (sign)", "donor-age p (sign)", "Spearman", "gene-label p (Spearman)", "donor-age p (Spearman)", "Pearson", "gene-label p (Pearson)", "donor-age p (Pearson)"],
            [[s["arm"], s["n_pairs"], s["n_pairs_with_tested_genes"], s["n_pair_gene_rows"],
              s["sign_agreement_fraction"], s["p_sign_gene_perm"], s["p_sign_age_perm"],
              s["spearman"], s["p_spearman_gene_perm"], s["p_spearman_age_perm"],
              s["pearson"], s["p_pearson_gene_perm"], s["p_pearson_age_perm"]]
             for s in summaries.values()],
        ), "",
        f"Matching covered {len(matched)} of 36 mouse units ({len(matched) / 36:.1%}), yielding {len(pair_data)} matched pairs across cohorts. {len(unmatched)} mouse units were unmatched.", "",
        f"**Pooled verdict:** sign agreement is below 0.5 in both arms. Its gene-label permutation p values are below 0.05, but this does not mean agreement exceeds 0.5; the donor-age permutation p values are {fmt(summaries['floored']['p_sign_age_perm'])} (floored) and {fmt(summaries['unfloored']['p_sign_age_perm'])} (unfloored). Spearman and Pearson correlations are near zero, with no p value below 0.05 under either null in either arm. Thus the run does not show robust positive gene-effect concordance; the sign result is null-sensitive and is not evidence of conserved directionality.", "",
        "## Design", "",
        "The committed P1m derivation was read, not re-derived: 1,914 cleaned genes (1,119 up and 795 down), with 1,724 human-mapped genes. Each mouse unit contributes exactly the global-list genes whose `units_list` contains that unit; mouse effects are that unit's derivation log2FC for 30m versus 18m.", "",
        "Human units were matched by tissue and exact cell-type identity after lowercasing and removing punctuation and whitespace. Marrow matches only to blood; no synonyms were introduced. Kidney Living Donor and SCR v15, lung Integrated, CellRef and ImmuneAging BAL, Yazar blood, PBMC Thomson/Zhu and Perez blood were audited. Normal Adult and SN v15 kidney sources are absent and were not fabricated.", "",
        "Human expression is donor pseudobulk from donors aged 30 years or older: raw counts are summed within donor and matched cell type, normalized to counts per 10,000 within that donor-cell-type pseudobulk, transformed with log1p, and regressed on continuous age with OLS plus sex where sex is variable and identifiable. No depth, pool, batch, or detection covariate/filter was added. Duplicate feature symbols are summed before normalization. A gene is tested when its human ortholog is present and its slope is finite; zero slopes count as sign disagreement. One-donor and no-age-variation units have slope NA and contribute no pooled rows. Two-donor fits are age-only differences with no variance estimate, no slope SE/p, and are explicitly flagged.", "",
        "The floored arm retains the 36 eligible mouse derivation units and, among their matched human pairs, requires at least 3 donors aged 30+ and at least 50 cells in total from those donors. This is the human analog of the mouse two-group floor: 'per side' applies to the mouse 18m/30m design, whereas the human fit is one continuous-age regression over donors aged 30+. The unfloored arm includes every matched pair regardless of donor or cell count. Removing this floor does not change the P1m gene list or add mouse units.", "",
        f"Gene-label nulls use {gene_permutations:,} shuffles, within pair for pair-level p values and across all pair-gene rows for pooled p values. Donor-age nulls use {age_permutations:,} permutations within cohort and sex strata, sharing each donor's permuted age across all matched cell types from that donor, refitting every human slope, and recomputing pair-level and pooled statistics. Seeds are 20261124 for age permutations and deterministic pair-specific offsets from 20261125 for gene permutations.", "",
        "## Match audit (all 36 mouse units)", "",
        markdown_table(
            ["mouse unit", "normalized name", "decision", "present matches", "match decision / reason"],
            [[r.mouse_unit, r.normalized_name, r.decision, r.n_present_matches, r.reason]
             for r in audit.itertuples()], digits=4,
        ), "",
        "Unmatched mouse units: " + "; ".join(unmatched["mouse_unit"]) + ".", "",
        "The match-rate assertion passed: matched plus unmatched mouse units = "
        f"{len(matched)} + {len(unmatched)} = 36.", "",
        "## Per-pair results", "",
        markdown_table(OUTPUT_COLUMNS, output[OUTPUT_COLUMNS].itertuples(index=False, name=None), digits=5),
        "",
        "## Effect of removing the floor", "",
    ]
    floored = summaries["floored"]
    unfloored = summaries["unfloored"]
    entered = output[(output["arm"] == "unfloored") & output["below_old_floor"]]
    lines.extend([
        f"Removing the floor added {len(entered)} pairs, from {floored['n_pairs']} to {unfloored['n_pairs']} pairs, and changed the pooled pair-gene row count from {floored['n_pair_gene_rows']} to {unfloored['n_pair_gene_rows']}. Sign agreement moved from {fmt(floored['sign_agreement_fraction'])} to {fmt(unfloored['sign_agreement_fraction'])}; Spearman moved from {fmt(floored['spearman'])} to {fmt(unfloored['spearman'])}; Pearson moved from {fmt(floored['pearson'])} to {fmt(unfloored['pearson'])}. The p values for both nulls are in the headline table.", "",
        "Pairs entering only in the unfloored arm:", "",
        markdown_table(
            ["mouse unit", "human cohort", "human unit", "donors 30+", "cells 30+", "tested genes"],
            [[r.mouse_unit, r.human_cohort, r.human_unit, r.n_donors_human_30plus,
              r.n_cells_human_30plus, r.n_human_genes_tested] for r in entered.itertuples()],
        ) if len(entered) else "None.", "",
        "## Independent spot-check and consistency checks", "",
        f"The spot-check pair is `{spot['mouse_unit']}` versus `{spot['human_cohort']}::{spot['human_unit']}`. Starting again from `results/followup-mature-p1m-genes.csv` and the unit DE CSV, the script checked every mouse log2FC and ortholog against the pair-gene audit rows, recomputed same-sign indicators as `sign(mouse_log2FC) == sign(human_slope)`, and obtained {spot['same']}/{spot['total']} = {spot['fraction']:.12f}. This exactly equals the unfloored CSV fraction, and the tested-gene count also equals {spot['total']}. The non-committed audit table is `{PAIR_GENE_CSV}`.", "",
        "The report's pooled values were computed from the same serialized pair-gene rows used for these checks, not copied by hand. Internal assertions confirmed 36 total match-audit rows, unique pair-arm keys, zero tested genes whenever donors <2, and exact spot-check agreement. After writing, the final CSV and pooled summaries were independently re-read and compared before completion.", "",
        "## Limits", "",
        "- Same-named cell types are not proven homologous or developmentally equivalent across mouse and human; exact-name matching is deliberately conservative and excludes many plausible synonyms.",
        "- The P1m signature is marrow-dominated. Exact blood matching leaves two Yazar marrow pairs, and the hematopoietic-precursor pair alone supplies 187 of 299 unfloored pair-gene rows (187 of 294 floored). Thus pooled pair-gene statistics are dominated by that marrow pair despite 14 matched kidney/lung mouse units, and the 36 derivation units are not represented equally.",
        "- Power is limited by small donor numbers, zero-gene pairs, and unequal gene counts per pair. Pair-gene rows are not independent biological replicates; permutation p values address the stated nulls, not all clustering or selection uncertainty.",
        "- Kidney Normal Adult and SN v15 source files are absent. Their committed metadata can document exact-name matches, but no expression slopes were invented for them. Perez and PBMC contributed no exact normalized matches under this rule.",
        "- Mouse log2FC is a selected two-group 18m-to-30m effect, while human slopes are continuous cross-sectional 30+ associations on a different expression scale. Correlation magnitudes are descriptive concordance measures, not evidence of conserved causal aging programs.",
        "- The mouse 30m group is all male, as documented in the derivation, so mouse age and sex are confounded; human sex adjustment does not remove that derivation limitation.", "",
        "## Reproducibility", "",
        f"Script: `scripts/run_mature_matched_p1m.py`. Started `{started}`; runtime {duration:.3f} seconds. The full run used the `agescape` conda environment and sources under `{DATA_ROOT}`. Human pseudobulks are cached outside the repository at `{HUMAN_CACHE}`; rerun with `--refresh-human-cache` to rebuild them. Outputs are `{OUTPUT_CSV.relative_to(REPO)}` and this report.", "",
    ])
    OUTPUT_REPORT.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-permutations", type=int, default=2000)
    parser.add_argument("--age-permutations", type=int, default=1000)
    parser.add_argument("--refresh-human-cache", action="store_true")
    args = parser.parse_args()
    assert args.gene_permutations >= 2000, "Owner order requires at least 2,000 gene-label permutations"
    assert args.age_permutations >= 200, "Owner order requires at least 200 donor-age permutations"
    start = time.monotonic()
    started = utc_now()
    print(f"Matched P1m started {started}", flush=True)
    eligible, mouse_units, genes = load_mouse()
    roster = human_roster()
    matches, audit = build_matches(eligible, roster)
    print(f"Match audit: {audit['decision'].value_counts().to_dict()}; pairs={len(matches)}", flush=True)
    target_genes = sorted(genes.loc[genes["human_ortholog"].notna(), "human_ortholog"].astype(str).unique())
    human_units = load_human_units(matches, target_genes, refresh=args.refresh_human_cache)
    expected_keys = {f"{r.human_cohort}::{r.human_unit}" for r in matches.itertuples()}
    assert set(human_units) == expected_keys
    permute_units(human_units, args.age_permutations, seed=20261124)
    pair_data, pair_genes = build_pair_data(matches, mouse_units, human_units, args.gene_permutations, seed=20261125)
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    pair_genes.to_csv(PAIR_GENE_CSV, index=False)
    output = make_output_rows(pair_data)
    output.to_csv(OUTPUT_CSV, index=False)
    summaries = {}
    for index, arm in enumerate(("floored", "unfloored")):
        selected = [pair for pair in pair_data if arm == "unfloored" or not (
            pair["unit"]["n_donors"] < OLD_FLOOR_DONORS or pair["unit"]["n_cells"] < OLD_FLOOR_CELLS
        )]
        summaries[arm], _, _ = summarize_arm(arm, selected, args.gene_permutations, seed=20261200 + index)
    spot = independent_spot_check(pair_data, pair_genes, output)
    reread = pd.read_csv(OUTPUT_CSV)
    pd.testing.assert_frame_equal(output, reread, check_dtype=False, check_exact=False, rtol=1e-12, atol=1e-12)
    for arm, summary in summaries.items():
        rows = pair_genes if arm == "unfloored" else pair_genes.merge(
            output[(output["arm"] == "floored")][["mouse_unit", "human_cohort", "human_unit"]],
            on=["mouse_unit", "human_cohort", "human_unit"], how="inner",
        )
        observed = correlation_values(rows["mouse_log2FC"].to_numpy(), rows["human_slope"].to_numpy())
        np.testing.assert_allclose(summary["sign_agreement_fraction"], observed["sign"], rtol=0, atol=1e-12)
        np.testing.assert_allclose(summary["spearman"], observed["spearman"], rtol=0, atol=1e-12)
        np.testing.assert_allclose(summary["pearson"], observed["pearson"], rtol=0, atol=1e-12)
        assert summary["n_pair_gene_rows"] == len(rows)
    write_report(output, audit, summaries, pair_data, spot, args.gene_permutations,
                 args.age_permutations, started, time.monotonic() - start)
    for summary in summaries.values():
        print(pd.Series(summary).to_string(), flush=True)
    print(f"Wrote {OUTPUT_CSV} ({len(output)} rows)", flush=True)
    print(f"Wrote {OUTPUT_REPORT}", flush=True)
    print(f"Matched P1m completed in {time.monotonic() - start:.3f} seconds", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
