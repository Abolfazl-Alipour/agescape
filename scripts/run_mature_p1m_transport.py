#!/usr/bin/env python3
"""Run the original M9 and M10 composite transport with P1m and P0m at 30-plus."""

from __future__ import annotations

import gc
import hashlib
import os
import re
import sys
import time
import warnings
from datetime import datetime, timezone
from pathlib import Path

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_var] = "1"

import h5py
import numpy as np
import pandas as pd
import rdata
import scipy.sparse as sp
from anndata.io import read_elem
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_m10_pipeline import ols as m10_ols
from scripts.run_m9_pipeline import ols as m9_ols

RESULTS = REPO / "results"
DATA_ROOT = Path("/home/reza/projects/agescape-revision-data")
HUMAN_H5ADS = DATA_ROOT / "human-h5ads"
MSSM_H5AD = DATA_ROOT / "mssm.h5ad"
RDS = REPO / "data" / "external" / "dlpfc-atlas" / "lifespan_pseudobulk.rds"
RDS_MD5 = "aef628e7bc2e5134b790d98c6f57ec48"
TRANSPORT_CSV = RESULTS / "followup-mature-p1m-transport.csv"
UNTESTABLE_CSV = RESULTS / "followup-mature-p1m-transport-untestable.csv"
SCORES_CSV = RESULTS / "followup-mature-p1m-m9-donor-scores.csv"
DONOR_SUMMARY_CSV = RESULTS / "followup-mature-p1m-donor-summary.csv"
RUNTIME_CSV = RESULTS / "followup-mature-p1m-transport-runtime.csv"
MODEL = "kimi-code/k3.1-pro"
MIN_AGE = 30.0
MIN_CELLS_PER_DONOR_CT = 20
MIN_DONORS_M9 = 4
AGE_RE = re.compile(r"(\d+)-year-old")
TIERS = ("P0m", "P1m")
COHORTS = {
    "MSSM_DLPFC": ("37a17b78-4864-4a42-b67b-31c00962795a", MSSM_H5AD),
    "HBCC_DLPFC": ("5c97eeeb-7e52-44b3-b010-b832b1f5424c", HUMAN_H5ADS / "5c97eeeb-7e52-44b3-b010-b832b1f5424c.h5ad"),
    "RADC_DLPFC": ("5e57cd50-8e42-42d6-940d-5c1660d06864", HUMAN_H5ADS / "5e57cd50-8e42-42d6-940d-5c1660d06864.h5ad"),
    "SEAAD_DLPFC": ("6f7fd0f1-a2ed-4ff1-80d3-33dde731cbc3", HUMAN_H5ADS / "6f7fd0f1-a2ed-4ff1-80d3-33dde731cbc3.h5ad"),
    "Kidney_NormalAdult": ("0b4a15a7-4e9e-4555-9733-2423e5c66469", HUMAN_H5ADS / "0b4a15a7-4e9e-4555-9733-2423e5c66469.h5ad"),
    "Kidney_LivingDonor": ("0bae7ebf-eb54-46a6-be9a-3461cecefa4c", HUMAN_H5ADS / "0bae7ebf-eb54-46a6-be9a-3461cecefa4c.h5ad"),
    "Kidney_SN_v15": ("a12ccb9b-4fbe-457d-8590-ac78053259ef", HUMAN_H5ADS / "a12ccb9b-4fbe-457d-8590-ac78053259ef.h5ad"),
    "Kidney_SCR_v15": ("dea717d4-7bc0-4e46-950f-fd7e1cc8df7d", HUMAN_H5ADS / "dea717d4-7bc0-4e46-950f-fd7e1cc8df7d.h5ad"),
    "Lung_CellRef": ("eb499fd8-7000-419f-8854-926b9b61b11f", HUMAN_H5ADS / "eb499fd8-7000-419f-8854-926b9b61b11f.h5ad"),
    "Lung_Integrated": ("9f222629-9e39-47d0-b83f-e08d610c7479", HUMAN_H5ADS / "9f222629-9e39-47d0-b83f-e08d610c7479.h5ad"),
    "ImmuneAging_BAL": ("1b350d0a-4535-4879-beb6-1142f3f94947", HUMAN_H5ADS / "1b350d0a-4535-4879-beb6-1142f3f94947.h5ad"),
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_programs() -> dict[str, dict[str, list[str]]]:
    expected = {
        "P0m": ("followup-mature-p1m-p0m-genes.csv", 2598, 1474, 1124, 2049, 1157, 892),
        "P1m": ("followup-mature-p1m-genes.csv", 1914, 1119, 795, 1724, 1013, 711),
    }
    programs = {}
    for tier, (filename, n_mouse, n_up, n_down, n_mapped, n_human_up, n_human_down) in expected.items():
        frame = pd.read_csv(RESULTS / filename)
        assert len(frame) == n_mouse
        assert int((frame["direction"] == "UP").sum()) == n_up
        assert int((frame["direction"] == "DOWN").sum()) == n_down
        mapped = frame[frame["human_ortholog"].notna()].copy()
        mapped["human_ortholog"] = mapped["human_ortholog"].astype(str)
        assert len(mapped) == mapped["human_ortholog"].nunique() == n_mapped
        up = mapped.loc[mapped["direction"] == "UP", "human_ortholog"].tolist()
        down = mapped.loc[mapped["direction"] == "DOWN", "human_ortholog"].tolist()
        assert (len(up), len(down)) == (n_human_up, n_human_down)
        assert not set(up) & set(down)
        programs[tier] = {"up": up, "down": down}
    return programs


def original_family() -> pd.DataFrame:
    family = pd.read_csv(RESULTS / "human-transport-global-multiplicity.csv")
    assert len(family) == 252
    assert not family.duplicated(["cohort", "cell_type"]).any()
    assert (family["source"] == "M9").sum() == 225
    assert (family["source"] == "M10").sum() == 27
    return family[["cohort", "cell_type", "source"]].rename(columns={"cell_type": "unit"})


def read_obs_columns(handle: h5py.File) -> dict[str, np.ndarray]:
    names = ("donor_id", "development_stage", "disease", "cell_type", "sex", "is_primary_data")
    return {name: np.asarray(read_elem(handle["obs"][name])) for name in names}


def read_csr_block(group: h5py.Group, start: int, stop: int, n_genes: int) -> sp.csr_matrix:
    assert group.attrs["encoding-type"] == "csr_matrix"
    indptr = group["indptr"][start:stop + 1]
    lo = int(indptr[0])
    hi = int(indptr[-1])
    return sp.csr_matrix(
        (group["data"][lo:hi], group["indices"][lo:hi], indptr - lo),
        shape=(stop - start, n_genes),
    )


def donor_metadata(columns: dict[str, np.ndarray]) -> pd.DataFrame:
    age = pd.Series(columns["development_stage"].astype(str)).str.extract(AGE_RE)[0].to_numpy(float)
    primary = columns["is_primary_data"].astype(bool)
    normal = columns["disease"].astype(str) == "normal"
    keep = primary & normal & np.isfinite(age) & (age >= MIN_AGE)
    frame = pd.DataFrame({
        "donor_id": columns["donor_id"].astype(str)[keep],
        "age": age[keep],
        "sex": columns["sex"].astype(str)[keep],
    })
    donors = frame.groupby("donor_id", sort=True).agg(age=("age", "first"), sex=("sex", "first"))
    assert len(donors) == 0 or donors["age"].min() >= MIN_AGE
    return donors


def score_m9_cohort(cohort: str, path: Path, programs: dict[str, dict[str, list[str]]],
                    chunk_size: int = 5000) -> tuple[pd.DataFrame, dict[str, object]]:
    print(f"{cohort}: reading {path}", flush=True)
    with h5py.File(path, "r") as handle:
        columns = read_obs_columns(handle)
        donors = donor_metadata(columns)
        donor_arr = columns["donor_id"].astype(str)
        cell_type_arr = columns["cell_type"].astype(str)
        primary = columns["is_primary_data"].astype(bool)
        donor_codes, donor_levels = pd.factorize(donor_arr, sort=False)
        cell_codes, cell_levels = pd.factorize(cell_type_arr, sort=False)
        n_cell_levels = len(cell_levels)
        selected = primary & np.isin(donor_arr, donors.index.to_numpy()) & (cell_type_arr != "unknown")
        var_names = np.asarray(read_elem(handle["var"]["feature_name"])).astype(str)
        gene_sets = {
            tier: {
                direction: np.flatnonzero(np.isin(var_names, programs[tier][direction])).astype(int)
                for direction in ("up", "down")
            }
            for tier in TIERS
        }
        gene_indices = np.array(sorted({
            int(i) for tier in TIERS for direction in ("up", "down") for i in gene_sets[tier][direction]
        }), dtype=int)
        local_lookup = {gene: i for i, gene in enumerate(gene_indices)}
        local_sets = {
            tier: {
                direction: np.array([local_lookup[i] for i in gene_sets[tier][direction]], dtype=int)
                for direction in ("up", "down")
            }
            for tier in TIERS
        }
        n_obs = len(donor_arr)
        pair_sums = np.zeros((len(donor_levels) * n_cell_levels, len(TIERS)), dtype=np.float64)
        pair_counts = np.zeros(len(donor_levels) * n_cell_levels, dtype=np.int64)
        for start in range(0, n_obs, chunk_size):
            stop = min(start + chunk_size, n_obs)
            keep = selected[start:stop]
            if not keep.any():
                continue
            block = read_csr_block(handle["X"], start, stop, len(var_names))
            library = np.asarray(block.sum(axis=1)).ravel()
            matrix = block[:, gene_indices].astype(np.float64)
            scale = 1e6 / np.where(library > 0, library, 1.0)
            matrix.data *= np.repeat(scale, np.diff(matrix.indptr))
            np.log1p(matrix.data, out=matrix.data)
            pair_code = (donor_codes[start:stop] * n_cell_levels + cell_codes[start:stop])[keep]
            unique, inverse, counts = np.unique(pair_code, return_inverse=True, return_counts=True)
            pair_counts[unique] += counts
            for tier_index, tier in enumerate(TIERS):
                up = local_sets[tier]["up"]
                down = local_sets[tier]["down"]
                up_mean = np.asarray(matrix[keep][:, up].mean(axis=1)).ravel() if len(up) else np.zeros(keep.sum())
                down_mean = np.asarray(matrix[keep][:, down].mean(axis=1)).ravel() if len(down) else np.zeros(keep.sum())
                values = up_mean - down_mean
                pair_sums[unique, tier_index] += np.bincount(
                    inverse, weights=values, minlength=len(unique)
                )
            del block, matrix
            if start % 500000 == 0:
                print(f"  {cohort}: cells {start} to {stop} of {n_obs}", flush=True)
        present = np.flatnonzero(pair_counts)
        scores = pd.DataFrame({
            "donor_id": donor_levels[present // n_cell_levels].astype(str),
            "cell_type": cell_levels[present % n_cell_levels].astype(str),
            "n_cells": pair_counts[present],
        })
        for tier_index, tier in enumerate(TIERS):
            scores[f"score_{tier}"] = pair_sums[present, tier_index] / pair_counts[present]
        scores["age"] = scores["donor_id"].map(donors["age"])
        scores["sex"] = scores["donor_id"].map(donors["sex"])
        scores["cohort"] = cohort
        assert scores["age"].min() >= MIN_AGE
        for tier in TIERS:
            scores[f"n_up_detected_{tier}"] = len(gene_sets[tier]["up"])
            scores[f"n_down_detected_{tier}"] = len(gene_sets[tier]["down"])
        scores = scores[scores["n_cells"] >= MIN_CELLS_PER_DONOR_CT].copy()
        info = {
            "cohort": cohort,
            "n_donors_30plus": len(donors),
            "age_min_30plus": float(donors["age"].min()),
            "age_max_30plus": float(donors["age"].max()),
            "n_score_rows": len(scores),
        }
    validate_m9_cell_counts(cohort, scores)
    print(f"  {cohort}: {len(donors)} donors at 30-plus, {len(scores)} donor-cell-type rows", flush=True)
    del columns
    gc.collect()
    return scores, info


def validate_m9_cell_counts(cohort: str, scores: pd.DataFrame) -> None:
    original = pd.read_csv(RESULTS / "m9-donor-scores.csv")
    original = original[(original["cohort"] == cohort) & (original["age"] >= MIN_AGE)].copy()
    original["donor_id"] = original["donor_id"].astype(str)
    current = scores[["donor_id", "cell_type", "n_cells"]]
    joined = original[["donor_id", "cell_type", "n_cells"]].merge(
        current, on=["donor_id", "cell_type"], how="outer", suffixes=("_original", "_new"), indicator=True
    )
    assert (joined["_merge"] == "both").all(), f"{cohort}: local donor-cell-type pairs differ from M9"
    assert (joined["n_cells_original"] == joined["n_cells_new"]).all(), f"{cohort}: local cell counts differ from M9"


def fit_m9_unit(sub: pd.DataFrame, score_column: str) -> tuple[dict[str, float], bool]:
    donors = sub.groupby("donor_id", observed=True).agg(
        score=(score_column, "mean"), age=("age", "first"), sex=("sex", "first")
    ).dropna(subset=["age"])
    male = (donors["sex"].astype(str).str.lower() == "male").astype(float).to_numpy()
    columns = [np.ones(len(donors)), donors["age"].to_numpy(float)]
    sex_included = len(np.unique(male)) > 1
    if sex_included:
        columns.append(male)
    return m9_ols(donors["score"].to_numpy(), np.column_stack(columns)), sex_included


def transport_row(cohort: str, unit: str, tier: str, source: str, model: str,
                  programs: dict[str, dict[str, list[str]]]) -> dict[str, object]:
    return {
        "cohort": cohort,
        "unit": unit,
        "gene_set": tier,
        "source": source,
        "model": model,
        "status": "TESTED",
        "composite_slope": np.nan,
        "composite_slope_per_year": np.nan,
        "composite_pvalue": np.nan,
        "composite_q": np.nan,
        "composite_pass": False,
        "n_donors": np.nan,
        "n_donors_30plus_cohort": np.nan,
        "age_min": np.nan,
        "age_max": np.nan,
        "age_sd": np.nan,
        "sex_included": pd.NA,
        "n_genes_mapped": len(set(programs[tier]["up"]) | set(programs[tier]["down"])),
        "n_up_mapped": len(programs[tier]["up"]),
        "n_down_mapped": len(programs[tier]["down"]),
        "n_genes_detected": np.nan,
        "n_up_detected": np.nan,
        "n_down_detected": np.nan,
        "multiplicity_family": "mature_P1m_P0m_M9_M10_30plus",
    }


def untestable_row(cohort: str, unit: str, source: str, reason: str,
                   n_donors: int, n_cohort_donors: int, age_min: float = np.nan,
                   age_max: float = np.nan) -> dict[str, object]:
    dataset_id = COHORTS.get(cohort, ("", ""))[0]
    source_object = str(COHORTS[cohort][1]) if cohort in COHORTS else str(RDS)
    return {
        "cohort": cohort,
        "unit": unit,
        "source": source,
        "dataset_id": dataset_id,
        "source_object": source_object,
        "status": "UNTESTABLE",
        "reason": reason,
        "n_donors_30plus": int(n_donors),
        "n_donors_30plus_cohort": int(n_cohort_donors),
        "age_min_30plus": age_min,
        "age_max_30plus": age_max,
    }


def historical_unit_meta(cohort: str, unit: str, prethreshold: pd.DataFrame) -> tuple[int, float, float]:
    sub = prethreshold[
        (prethreshold["cohort"] == cohort) & (prethreshold["cell_type"] == unit)
        & (prethreshold["age"] >= MIN_AGE) & (prethreshold["n_cells"] >= MIN_CELLS_PER_DONOR_CT)
    ]
    donors = sub.drop_duplicates("donor_id")
    return len(donors), float(donors["age"].min()), float(donors["age"].max())


def run_m9(family: pd.DataFrame, programs: dict[str, dict[str, list[str]]]):
    all_scores = []
    infos = []
    rows = []
    untestable = []
    prethreshold = pd.read_csv(RESULTS / "m9-donor-cell-counts-prethreshold.csv")
    cohort_donor_counts = prethreshold[prethreshold["age"] >= MIN_AGE].groupby("cohort")["donor_id"].nunique()
    for cohort, (dataset_id, path) in COHORTS.items():
        units = family.loc[family["cohort"] == cohort, "unit"].tolist()
        if not path.exists():
            print(f"{cohort}: UNTESTABLE, source object absent: {path}", flush=True)
            for unit in units:
                n_donors, age_min, age_max = historical_unit_meta(cohort, unit, prethreshold)
                untestable.append(untestable_row(
                    cohort, unit, "M9", f"source object absent: {dataset_id}", n_donors,
                    int(cohort_donor_counts.get(cohort, 0)), age_min, age_max
                ))
            infos.append({
                "cohort": cohort,
                "n_donors_30plus": int(cohort_donor_counts.get(cohort, 0)),
                "age_min_30plus": np.nan,
                "age_max_30plus": np.nan,
                "n_score_rows": 0,
            })
            continue
        scores, info = score_m9_cohort(cohort, path, programs)
        all_scores.append(scores)
        infos.append(info)
        for unit in units:
            sub = scores[(scores["cell_type"] == unit) & (scores["n_cells"] >= MIN_CELLS_PER_DONOR_CT)].copy()
            n_donors = sub["donor_id"].nunique()
            if n_donors < MIN_DONORS_M9:
                untestable.append(untestable_row(
                    cohort, unit, "M9",
                    f"only {n_donors} donors at 30-plus with at least {MIN_CELLS_PER_DONOR_CT} cells; M9 requires {MIN_DONORS_M9}",
                    n_donors, info["n_donors_30plus"], info["age_min_30plus"], info["age_max_30plus"]
                ))
                continue
            donor_ages = sub.groupby("donor_id")["age"].first()
            for tier in TIERS:
                fit, sex_included = fit_m9_unit(sub, f"score_{tier}")
                row = transport_row(cohort, unit, tier, "M9", "age_plus_sex", programs)
                row.update({
                    "composite_slope": fit["slope"],
                    "composite_slope_per_year": fit["slope"],
                    "composite_pvalue": fit["p"],
                    "n_donors": n_donors,
                    "n_donors_30plus_cohort": info["n_donors_30plus"],
                    "age_min": donor_ages.min(),
                    "age_max": donor_ages.max(),
                    "age_sd": donor_ages.std(),
                    "sex_included": sex_included,
                    "n_genes_detected": int(sub[f"n_up_detected_{tier}"].iloc[0] + sub[f"n_down_detected_{tier}"].iloc[0]),
                    "n_up_detected": int(sub[f"n_up_detected_{tier}"].iloc[0]),
                    "n_down_detected": int(sub[f"n_down_detected_{tier}"].iloc[0]),
                })
                rows.append(row)
    scores_out = pd.concat(all_scores, ignore_index=True)
    return scores_out, infos, rows, untestable


def psychad_design(metadata: dict[str, np.ndarray], keep: np.ndarray) -> np.ndarray:
    return np.column_stack([
        np.ones(keep.sum()), metadata["Age"][keep].astype(float),
        (metadata["Sex"][keep] == "Male").astype(float),
        (metadata["Brain_bank"][keep] == "HBCC").astype(float),
        metadata["PMI"][keep].astype(float), metadata["TechPC1"][keep].astype(float),
        metadata["TechPC2"][keep].astype(float), metadata["TechPC3"][keep].astype(float),
    ])


def run_m10(family: pd.DataFrame, programs: dict[str, dict[str, list[str]]]):
    digest = md5sum(RDS)
    assert digest == RDS_MD5, f"PsychAD MD5 mismatch: {digest}"
    print(f"PsychAD: MD5 verified {digest}", flush=True)
    obj = rdata.read_rds(str(RDS))
    cd = obj.colData
    metadata = {
        name: np.asarray(cd.listData[name]).ravel()
        for name in ("Age", "Sex", "Brain_bank", "PMI", "TechPC1", "TechPC2", "TechPC3")
    }
    mature = metadata["Age"].astype(float) >= MIN_AGE
    assert int(mature.sum()) > 0
    rows = []
    untestable = []
    for unit in family.loc[family["source"] == "M10", "unit"].tolist():
        source = obj.assays.data.listData[unit]
        genes = np.array([str(gene) for gene in source.Dimnames[0]])
        matrix = sp.csc_matrix(
            (np.asarray(source.x), np.asarray(source.i), np.asarray(source.p)), shape=tuple(source.Dim)
        )
        library = np.asarray(matrix.sum(axis=0)).ravel()
        keep = mature & (library > 0)
        n_donors = int(keep.sum())
        design = psychad_design(metadata, keep)
        if n_donors <= design.shape[1] or np.linalg.matrix_rank(design) < design.shape[1]:
            untestable.append(untestable_row(
                "DLPFC_Lifespan_Atlas", unit, "M10",
                f"30-plus fully adjusted design is ineligible: {n_donors} donors, {design.shape[1]} columns, rank {np.linalg.matrix_rank(design)}",
                n_donors, int(mature.sum())
            ))
            continue
        normalized = matrix[:, keep] @ sp.diags(1e6 / library[keep])
        normalized = normalized.tocsc()
        normalized.data = np.log1p(normalized.data)
        lookup = {gene: i for i, gene in enumerate(genes)}
        for tier in TIERS:
            up = np.array([lookup[g] for g in programs[tier]["up"] if g in lookup], dtype=int)
            down = np.array([lookup[g] for g in programs[tier]["down"] if g in lookup], dtype=int)
            up_mean = np.asarray(normalized[up, :].mean(axis=0)).ravel() if len(up) else np.zeros(n_donors)
            down_mean = np.asarray(normalized[down, :].mean(axis=0)).ravel() if len(down) else np.zeros(n_donors)
            fit = m10_ols(up_mean - down_mean, design)
            row = transport_row("DLPFC_Lifespan_Atlas", unit, tier, "M10", "fully_adjusted", programs)
            ages = metadata["Age"][keep].astype(float)
            row.update({
                "composite_slope": fit["slope"],
                "composite_slope_per_year": fit["slope"],
                "composite_pvalue": fit["p"],
                "n_donors": n_donors,
                "n_donors_30plus_cohort": int(mature.sum()),
                "age_min": ages.min(),
                "age_max": ages.max(),
                "age_sd": ages.std(ddof=1),
                "sex_included": True,
                "n_genes_detected": len(up) + len(down),
                "n_up_detected": len(up),
                "n_down_detected": len(down),
            })
            rows.append(row)
    info = {
        "cohort": "DLPFC_Lifespan_Atlas",
        "n_donors_30plus": int(mature.sum()),
        "age_min_30plus": float(metadata["Age"][mature].astype(float).min()),
        "age_max_30plus": float(metadata["Age"][mature].astype(float).max()),
        "n_score_rows": 0,
    }
    return rows, untestable, info, digest


def main() -> int:
    start = time.monotonic()
    start_utc = utc_now()
    print(f"Mature transport started {start_utc}; model {MODEL}", flush=True)
    programs = load_programs()
    family = original_family()
    m9_scores, infos, m9_rows, m9_untestable = run_m9(family, programs)
    m10_rows, m10_untestable, m10_info, digest = run_m10(family, programs)
    transport = pd.DataFrame(m9_rows + m10_rows)
    for tier in TIERS:
        mask = transport["gene_set"] == tier
        transport.loc[mask, "composite_q"] = multipletests(
            transport.loc[mask, "composite_pvalue"], method="fdr_bh"
        )[1]
    transport["composite_pass"] = transport["composite_q"] < 0.05
    transport = transport.sort_values(["source", "cohort", "unit", "gene_set"]).reset_index(drop=True)
    untestable = pd.DataFrame(m9_untestable + m10_untestable).sort_values(["source", "cohort", "unit"]).reset_index(drop=True)
    tested_keys = set(zip(transport["cohort"], transport["unit"]))
    untestable_keys = set(zip(untestable["cohort"], untestable["unit"]))
    original_keys = set(zip(family["cohort"], family["unit"]))
    assert tested_keys.isdisjoint(untestable_keys)
    assert tested_keys | untestable_keys == original_keys
    assert transport.groupby("gene_set").size().nunique() == 1
    assert transport["age_min"].min() >= MIN_AGE
    m9_scores.to_csv(SCORES_CSV, index=False)
    transport.to_csv(TRANSPORT_CSV, index=False)
    untestable.to_csv(UNTESTABLE_CSV, index=False)
    infos = pd.DataFrame(infos + [m10_info])
    original_counts = family.groupby("cohort").size()
    tested_counts = transport[transport["gene_set"] == "P1m"].groupby("cohort").size()
    untestable_counts = untestable.groupby("cohort").size()
    infos["n_original_units"] = infos["cohort"].map(original_counts)
    infos["n_tested_units"] = infos["cohort"].map(tested_counts).fillna(0).astype(int)
    infos["n_untestable_units"] = infos["cohort"].map(untestable_counts).fillna(0).astype(int)
    infos.to_csv(DONOR_SUMMARY_CSV, index=False)
    end_utc = utc_now()
    duration = time.monotonic() - start
    pd.DataFrame([{
        "model": MODEL,
        "start_utc": start_utc,
        "end_utc": end_utc,
        "duration_seconds": duration,
        "n_original_units": len(family),
        "n_tested_units_per_tier": len(tested_keys),
        "n_untestable_units": len(untestable),
        "psychad_md5": digest,
    }]).to_csv(RUNTIME_CSV, index=False)
    print(f"Transport family: {len(tested_keys)} tested and {len(untestable)} untestable per tier", flush=True)
    print(transport[transport["composite_pass"]][["cohort", "unit", "gene_set", "composite_pvalue", "composite_q"]].to_string(index=False), flush=True)
    print(f"Mature transport ended {end_utc}; duration {duration:.3f} seconds", flush=True)
    return 0


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    raise SystemExit(main())
