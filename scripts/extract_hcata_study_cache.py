#!/usr/bin/env python3
"""Audit .obs, verify raw integer counts, and extract compact pseudobulk caches for HCATA studies."""

from __future__ import annotations

import argparse
import gc
import json
from pathlib import Path
import pickle
import sys
from typing import Any

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem
from scipy import sparse as sp

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.followup_muscle import count_gate

MIN_DONORS = 6
MIN_CELLS_PER_DONOR = 20
CHUNK_ROWS = 25_000
MAX_GROUP_BATCH = 6  # keeps dense (batch, n_donors, n_vars) float64 well below 2 GB even at N=981


def summarize_obs_columns(obs: pd.DataFrame) -> dict[str, Any]:
    out: dict[str, Any] = {"n_obs": int(len(obs)), "columns": list(obs.columns)}
    for col in obs.columns:
        if col in {"observation_joinid"}:
            continue
        nunq = int(obs[col].nunique(dropna=False))
        if nunq <= 50:
            vc = obs[col].astype(str).value_counts(dropna=False)
            out[f"col_{col}"] = {str(k): int(v) for k, v in vc.items() if int(v) > 0}
        else:
            out[f"col_{col}_nunique"] = nunq
    return out


def study_cell_mask_and_meta(study: str, obs: pd.DataFrame) -> tuple[pd.Series, pd.DataFrame, dict[str, Any]]:
    """Return boolean inclusion mask, per-donor metadata DataFrame, and audit dictionary."""
    base_audit = summarize_obs_columns(obs)
    if study == "guilliams":
        mask = (
            (obs["disease"].astype(str) == "normal")
            & (obs["assay"].astype(str) == "10x 3' v3")
            & (obs["suspension_type"].astype(str) == "cell")
            & (obs["development_stage"].astype(str).str.endswith("-year-old stage"))
        )
    elif study in {"yazar", "perez"}:
        assay_target = "10x 3' v2"
        mask = (
            (obs["disease"].astype(str) == "normal")
            & (obs["assay"].astype(str) == assay_target)
            & (obs["suspension_type"].astype(str) == "cell")
            & (obs["development_stage"].astype(str).str.endswith("-year-old stage"))
        )
    else:
        raise ValueError(f"Unknown study {study!r}")

    sub = obs.loc[mask].copy()
    sub["_donor"] = sub["donor_id"].astype(str)
    sub["_ct"] = sub["cell_type"].astype(str)
    sub["_stage"] = sub["development_stage"].astype(str)
    sub["_sex"] = sub["sex"].astype(str)

    donor_rows = []
    for did, dsub in sub.groupby("_donor", observed=True):
        stages = sorted(dsub["_stage"].unique().tolist())
        sexes = sorted(dsub["_sex"].unique().tolist())
        assert len(stages) == 1, f"{study} {did}: inconsistent stages {stages}"
        assert len(sexes) == 1, f"{study} {did}: inconsistent sexes {sexes}"
        age_years = float(int(stages[0].split("-")[0]))
        sex_val = sexes[0]
        assert sex_val in {"male", "female"}, f"{study} {did}: unexpected sex {sex_val!r}"
        rec: dict[str, Any] = {
            "donor_id": did,
            "age_years": age_years,
            "age_per_decade": age_years / 10.0,
            "sex": sex_val,
            "sex_male": int(sex_val == "male"),
            "n_cells": int(len(dsub)),
        }
        for extra in [
            "reported_diseases",
            "diet",
            "suspension_enrichment_factors",
            "typeSample",
            "self_reported_ethnicity",
            "pool_number",
            "Processing_Cohort",
        ]:
            if extra in dsub.columns:
                vals = sorted(set(dsub[extra].astype(str)))
                rec[extra] = ";".join(vals)
        donor_rows.append(rec)

    donor_meta = pd.DataFrame(donor_rows).sort_values(["age_years", "donor_id"]).set_index("donor_id")
    base_audit["included_cells"] = int(len(sub))
    base_audit["included_donors"] = int(len(donor_meta))
    base_audit["age_min"] = float(donor_meta["age_years"].min())
    base_audit["age_max"] = float(donor_meta["age_years"].max())
    base_audit["n_male"] = int(donor_meta["sex_male"].sum())
    base_audit["n_female"] = int((1 - donor_meta["sex_male"]).sum())
    return mask, donor_meta, base_audit


def extract_study(study: str, h5_path: Path, cache_dir: Path) -> None:
    cache_dir.mkdir(parents=True, exist_ok=True)
    summary_path = cache_dir / "obs_audit.json"
    pb_manifest = cache_dir / "pb_manifest.json"
    if summary_path.exists() and pb_manifest.exists():
        print(f"[{study}] Cache already exists at {cache_dir}", flush=True)
        return

    with h5py.File(h5_path, "r") as f:
        obs = read_elem(f["obs"])
        if "raw" in f and "var" in f["raw"]:
            raw_grp = f["raw/X"]
            var = read_elem(f["raw/var"])
            raw_source = "raw/X"
        else:
            raw_grp = f["X"]
            var = read_elem(f["var"])
            raw_source = "X"

        enc = raw_grp.attrs.get("encoding-type", "")
        if isinstance(enc, bytes):
            enc = enc.decode()
        shape = tuple(int(x) for x in raw_grp.attrs["shape"])
        n_obs, n_vars = shape
        assert n_obs == len(obs), f"{study}: obs len {len(obs)} != {n_obs}"
        assert n_vars == len(var), f"{study}: var len {len(var)} != {n_vars}"
        assert enc == "csr_matrix", f"{study}: unexpected encoding {enc!r}"

        # Sequential full count_gate over every nonzero element
        data_ds = raw_grp["data"]
        total_nnz = int(data_ds.shape[0])
        gate = {"nonfinite": 0, "negative": 0, "noninteger": 0, "total_nnz": total_nnz, "n_obs": n_obs, "raw_source": raw_source}
        for s in range(0, total_nnz, 5_000_000):
            e = min(s + 5_000_000, total_nnz)
            g = count_gate(data_ds[s:e])
            gate["nonfinite"] += g["nonfinite"]
            gate["negative"] += g["negative"]
            gate["noninteger"] += g["noninteger"]
        print(f"[{study}] count_gate: {gate}", flush=True)
        if gate["nonfinite"] or gate["negative"] or gate["noninteger"]:
            raise ValueError(f"{study} failed count_gate: {gate}")

        mask, donor_meta, audit = study_cell_mask_and_meta(study, obs)
        audit["gate"] = gate

        sub_obs = obs.loc[mask]
        ct_counts = sub_obs.groupby(["donor_id", "cell_type"], observed=True).size().unstack(fill_value=0)
        ct_counts.index = ct_counts.index.astype(str)
        ct_counts.columns = ct_counts.columns.astype(str)
        ct_counts = ct_counts.reindex(donor_meta.index, fill_value=0)
        ct_counts.to_csv(cache_dir / "cell_counts.csv")
        donor_meta.to_csv(cache_dir / "donor_meta.csv")

        # Determine eligible cell types + All_Cells (+ Guilliams CD45+-only sensitivity pool)
        elig_cts = [
            ct for ct in sorted(ct_counts.columns)
            if int((ct_counts[ct] >= MIN_CELLS_PER_DONOR).sum()) >= MIN_DONORS
        ]
        groups_to_extract = list(elig_cts) + ["All_Cells"]
        cd45_bool = None
        if study == "guilliams":
            cd45_bool = (
                mask.to_numpy()
                & (obs["suspension_enrichment_factors"].astype(str).to_numpy() == "CD45+")
            )
            groups_to_extract.append("All_CD45pos_Cells")

        if "feature_name" in var.columns:
            gene_names = var["feature_name"].astype(str).to_numpy()
        else:
            gene_names = var.index.astype(str).to_numpy()

        donors_arr = obs["donor_id"].astype(str).to_numpy()
        ct_arr = obs["cell_type"].astype(str).to_numpy()
        keep_bool = mask.to_numpy()
        unique_donors = list(donor_meta.index)
        donor_idx_map = {d: i for i, d in enumerate(unique_donors)}
        donor_code_arr = np.array([donor_idx_map.get(d, -1) for d in donors_arr], dtype=np.int32)

        indptr_ds = raw_grp["indptr"]
        indices_ds = raw_grp["indices"]

        pb_files: dict[str, str] = {}
        # Stream in batches of <= MAX_GROUP_BATCH groups to keep dense buffers < 2 GB
        for b0 in range(0, len(groups_to_extract), MAX_GROUP_BATCH):
            batch_groups = groups_to_extract[b0 : b0 + MAX_GROUP_BATCH]
            print(f"[{study}] Streaming batch {b0//MAX_GROUP_BATCH + 1}: {batch_groups}", flush=True)
            acc = {g: np.zeros((len(unique_donors), n_vars), dtype=np.float64) for g in batch_groups}

            for r0 in range(0, n_obs, CHUNK_ROWS):
                r1 = min(r0 + CHUNK_ROWS, n_obs)
                cmask = keep_bool[r0:r1]
                if not np.any(cmask):
                    continue
                ptrs = indptr_ds[r0 : r1 + 1]
                p0 = int(ptrs[0])
                p1 = int(ptrs[-1])
                mat = sp.csr_matrix(
                    (data_ds[p0:p1].astype(np.float64, copy=False), indices_ds[p0:p1], ptrs - p0),
                    shape=(r1 - r0, n_vars),
                )
                sub_rows = np.flatnonzero(cmask)
                mat_sub = mat[sub_rows]
                d_sub = donor_code_arr[r0:r1][sub_rows]
                c_sub = ct_arr[r0:r1][sub_rows]

                for grp in batch_groups:
                    if grp == "All_Cells":
                        rows_g = np.arange(len(sub_rows))
                    elif grp == "All_CD45pos_Cells":
                        cd_sub = cd45_bool[r0:r1][sub_rows]
                        rows_g = np.flatnonzero(cd_sub)
                    else:
                        rows_g = np.flatnonzero(c_sub == grp)
                    if len(rows_g) == 0:
                        continue
                    dg = d_sub[rows_g]
                    mg = mat_sub[rows_g]
                    # Sparse binary indicator (n_donors_in_chunk x len(rows_g)) @ mg
                    uniq_d, inv_d = np.unique(dg, return_inverse=True)
                    ind = sp.csr_matrix(
                        (np.ones(len(rows_g), dtype=np.float64), (inv_d, np.arange(len(rows_g)))),
                        shape=(len(uniq_d), len(rows_g)),
                    )
                    summed = (ind @ mg).toarray()
                    acc[grp][uniq_d, :] += summed

            for grp in batch_groups:
                if grp == "All_Cells":
                    elig_donors = unique_donors
                elif grp == "All_CD45pos_Cells":
                    cd_counts = obs.loc[cd45_bool].groupby("donor_id", observed=True).size()
                    elig_donors = [d for d in unique_donors if int(cd_counts.get(d, 0)) >= MIN_CELLS_PER_DONOR]
                else:
                    elig_donors = [d for d in unique_donors if int(ct_counts.loc[d, grp]) >= MIN_CELLS_PER_DONOR]
                d_indices = [donor_idx_map[d] for d in elig_donors]
                mat_g = acc[grp][d_indices, :].T  # (n_vars, n_elig_donors)
                df_g = pd.DataFrame(mat_g, index=gene_names, columns=elig_donors)
                safe_name = (
                    grp.replace("/", "_")
                    .replace(" ", "_")
                    .replace(",", "")
                    .replace("+", "pos")
                )
                fname = f"pb_{safe_name}.pkl"
                with open(cache_dir / fname, "wb") as pf:
                    pickle.dump(df_g, pf, protocol=pickle.HIGHEST_PROTOCOL)
                pb_files[grp] = fname
            del acc
            gc.collect()

        if study == "guilliams":
            cd45_ser = obs.loc[cd45_bool].groupby("donor_id", observed=True).size()
            cd45_ser.index = cd45_ser.index.astype(str)
            cd45_ser.to_csv(cache_dir / "cd45pos_counts.csv", header=["n_cells"])

    audit["groups_extracted"] = list(pb_files.keys())
    audit["pb_files"] = pb_files
    summary_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    pb_manifest.write_text(json.dumps(pb_files, indent=2), encoding="utf-8")
    print(f"[{study}] Extraction complete -> {cache_dir} ({len(pb_files)} groups)", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", required=True, choices=["guilliams", "yazar", "perez"])
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--cache-dir", required=True, type=Path)
    args = parser.parse_args()
    extract_study(args.study, args.h5ad, args.cache_dir)


if __name__ == "__main__":
    main()
