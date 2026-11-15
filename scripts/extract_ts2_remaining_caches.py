#!/usr/bin/env python3
"""Download, audit .obs sorting/assays, verify integer raw.X, and cache pseudobulks for TS2 tissues >= 6 donors."""

from __future__ import annotations

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

from scripts import resource_guard
from src.followup_muscle import count_gate
from src.followup_organs import download_with_guard

CACHE_ROOT = Path("/home/reza/projects/agescape-revision-data/tmp_organs/ts2_rem")
EXISTING_BLOOD = Path("/home/reza/projects/agescape-revision-data/ts2/blood.h5ad")

TISSUE_ASSETS = [
    ("Blood", None, EXISTING_BLOOD, False),
    ("Bladder", "https://datasets.cellxgene.cziscience.com/8ffc335d-36fa-461e-be47-cd6503833420.h5ad", CACHE_ROOT / "bladder.h5ad", True),
    ("Fat", "https://datasets.cellxgene.cziscience.com/8115812f-9acf-4f62-b15f-bde18086a1c8.h5ad", CACHE_ROOT / "fat.h5ad", True),
    ("Lymph_Node", "https://datasets.cellxgene.cziscience.com/6dd28080-1d36-4905-8b44-a9954a4438a2.h5ad", CACHE_ROOT / "lymph_node.h5ad", True),
    ("Spleen", "https://datasets.cellxgene.cziscience.com/5498457f-6a18-4766-b9a3-a2c8c4662703.h5ad", CACHE_ROOT / "spleen.h5ad", True),
]

TARGET_ASSAY = "10x 3' v3"
CHUNK_ROWS = 20000


def verify_and_aggregate_h5ad(tissue: str, h5_path: Path) -> dict[str, Any]:
    """Audit .obs metadata, check raw.X count gate over all rows, and cache 10x 3' v3 raw pseudobulks."""
    t_dir = CACHE_ROOT / f"{tissue}_cache"
    t_dir.mkdir(parents=True, exist_ok=True)
    summary_file = t_dir / "obs_audit.json"
    pb_file = t_dir / "raw_pb_by_group.pkl"
    if summary_file.exists() and pb_file.exists():
        return json.loads(summary_file.read_text(encoding="utf-8"))

    with h5py.File(h5_path, "r") as f:
        obs = read_elem(f["obs"])
        if "raw" in f and "var" in f["raw"]:
            raw_grp = f["raw/X"]
            var = read_elem(f["raw/var"])
        else:
            raw_grp = f["X"]
            var = read_elem(f["var"])

        enc = raw_grp.attrs.get("encoding-type", "")
        if isinstance(enc, bytes):
            enc = enc.decode()
        shape = tuple(int(x) for x in raw_grp.attrs["shape"])
        n_obs, n_vars = shape
        assert n_obs == len(obs), f"{tissue}: obs mismatch {len(obs)} vs {n_obs}"
        assert n_vars == len(var), f"{tissue}: var mismatch {len(var)} vs {n_vars}"
        assert enc == "csr_matrix", f"{tissue}: unexpected raw encoding {enc!r}"

        # Check count gate sequentially across all stored values in raw/X
        data_ds = raw_grp["data"]
        total_nnz = int(data_ds.shape[0])
        gate = {"nonfinite": 0, "negative": 0, "noninteger": 0, "total_nnz": total_nnz, "n_obs": n_obs}
        for s in range(0, total_nnz, 5_000_000):
            e = min(s + 5_000_000, total_nnz)
            vals = data_ds[s:e]
            g = count_gate(vals)
            gate["nonfinite"] += g["nonfinite"]
            gate["negative"] += g["negative"]
            gate["noninteger"] += g["noninteger"]
        if gate["nonfinite"] or gate["negative"] or gate["noninteger"]:
            raise ValueError(f"{tissue} failed raw count_gate: {gate}")

        # Inspect obs tables for 10x 3' v3 and all assays
        assay_x_donor = (
            obs.groupby(["donor_id", "assay"], observed=True).size().unstack(fill_value=0)
        )
        tenx3_mask = (obs["assay"].astype(str) == TARGET_ASSAY) & (obs["disease"].astype(str) == "normal")
        tenx3 = obs.loc[tenx3_mask].copy()
        comp_x_donor = (
            tenx3.groupby(["donor_id", "compartment"], observed=True).size().unstack(fill_value=0)
        )
        notes_col = tenx3["notes"].astype(str) if "notes" in tenx3.columns else pd.Series("MISSING", index=tenx3.index)
        notes_x_donor = (
            tenx3.assign(_notes=notes_col)
            .groupby(["donor_id", "_notes"], observed=True)
            .size()
            .unstack(fill_value=0)
        )

        ct_counts = (
            tenx3.groupby(["donor_id", "cell_type"], observed=True).size().unstack(fill_value=0)
        )

        # Build donor meta from H5AD .obs
        donor_rows = []
        for did, sub in tenx3.groupby("donor_id", observed=True):
            ages = sorted(set(sub["development_stage"].astype(str)))
            sexes = sorted(set(sub["sex"].astype(str)))
            pos = sorted(set(sub["anatomical_position"].astype(str))) if "anatomical_position" in sub.columns else []
            susp = sorted(set(sub["suspension_type"].astype(str)))
            tip = sorted(set(sub["tissue_in_publication"].astype(str)))
            n_notes = dict(sub["notes"].astype(str).value_counts()) if "notes" in sub.columns else {}
            n_comp = dict(sub["compartment"].astype(str).value_counts())
            donor_rows.append({
                "donor_id": str(did),
                "development_stage": ";".join(ages),
                "sex": ";".join(sexes),
                "n_cells_10x3v3": int(len(sub)),
                "anatomical_positions": ";".join(pos),
                "suspension_types": ";".join(susp),
                "tissue_in_publication": ";".join(tip),
                "compartment_counts": {str(k): int(v) for k, v in n_comp.items()},
                "notes_counts": {str(k): int(v) for k, v in n_notes.items()},
            })

        # Gene symbols from var['feature_name'] (fallback index)
        if "feature_name" in var.columns:
            gene_names = var["feature_name"].astype(str).to_numpy()
        else:
            gene_names = var.index.astype(str).to_numpy()

        donors_arr = obs["donor_id"].astype(str).to_numpy()
        ct_arr = obs["cell_type"].astype(str).to_numpy()
        unique_donors = sorted(tenx3["donor_id"].astype(str).unique())
        unique_cts = sorted(tenx3["cell_type"].astype(str).unique())

        pair_keys = []
        donor_index = {d: i for i, d in enumerate(unique_donors)}
        ct_index = {c: i for i, c in enumerate(unique_cts)}

        # Allocate pseudobulk accumulators for every (cell_type, donor) plus All_Cells (+ All_NonNK_Cells for Blood)
        # Keep in dict of float64 matrices indexed by (group -> (n_donors, n_vars))
        group_names = unique_cts + (["All_NonNK_Cells"] if tissue == "Blood" else ["All_Cells"])
        acc = {g: np.zeros((len(unique_donors), n_vars), dtype=np.float64) for g in group_names}

        indptr_ds = raw_grp["indptr"]
        indices_ds = raw_grp["indices"]
        keep_bool = tenx3_mask.to_numpy()

        for r0 in range(0, n_obs, CHUNK_ROWS):
            r1 = min(r0 + CHUNK_ROWS, n_obs)
            chunk_mask = keep_bool[r0:r1]
            if not np.any(chunk_mask):
                continue
            ptrs = indptr_ds[r0 : r1 + 1]
            p0 = int(ptrs[0])
            p1 = int(ptrs[-1])
            chunk_indptr = ptrs - p0
            chunk_indices = indices_ds[p0:p1]
            chunk_data = data_ds[p0:p1].astype(np.float64, copy=False)
            mat = sp.csr_matrix((chunk_data, chunk_indices, chunk_indptr), shape=(r1 - r0, n_vars))

            sub_rows = np.flatnonzero(chunk_mask)
            mat_sub = mat[sub_rows]
            d_sub = donors_arr[r0:r1][sub_rows]
            c_sub = ct_arr[r0:r1][sub_rows]

            # All_Cells or All_NonNK_Cells
            all_grp = "All_NonNK_Cells" if tissue == "Blood" else "All_Cells"
            for d in set(d_sub):
                d_mask = d_sub == d
                if tissue == "Blood":
                    d_mask_all = d_mask & (c_sub != "natural killer cell")
                else:
                    d_mask_all = d_mask
                if np.any(d_mask_all):
                    acc[all_grp][donor_index[d]] += np.asarray(mat_sub[d_mask_all].sum(axis=0)).ravel()

            # Per cell type
            for c in set(c_sub):
                c_mask = c_sub == c
                for d in set(d_sub[c_mask]):
                    cd_mask = c_mask & (d_sub == d)
                    acc[c][donor_index[d]] += np.asarray(mat_sub[cd_mask].sum(axis=0)).ravel()

        raw_pb_by_group = {}
        for g, arr in acc.items():
            # transpose to genes x donors dataframe
            raw_pb_by_group[g] = pd.DataFrame(arr.T, index=gene_names, columns=unique_donors)

        with open(pb_file, "wb") as pf:
            pickle.dump(raw_pb_by_group, pf, protocol=pickle.HIGHEST_PROTOCOL)
        ct_counts.to_csv(t_dir / "cell_counts_10x3v3.csv")

        audit = {
            "tissue": tissue,
            "file_bytes": int(h5_path.stat().st_size),
            "obs_shape": [int(n_obs), int(obs.shape[1])],
            "n_vars": int(n_vars),
            "gate": gate,
            "is_primary_data_true_count": int(obs["is_primary_data"].sum()) if "is_primary_data" in obs.columns else -1,
            "unique_diseases": sorted(obs["disease"].astype(str).unique().tolist()),
            "unique_tissues": sorted(obs["tissue"].astype(str).unique().tolist()),
            "unique_tissue_in_pub": sorted(obs["tissue_in_publication"].astype(str).unique().tolist()),
            "unique_suspension_types": sorted(obs["suspension_type"].astype(str).unique().tolist()),
            "assay_x_donor": assay_x_donor.to_dict(),
            "comp_x_donor_10x3v3": comp_x_donor.to_dict(),
            "notes_x_donor_10x3v3": notes_x_donor.to_dict(),
            "donors_10x3v3": donor_rows,
        }
        summary_file.write_text(json.dumps(audit, indent=2), encoding="utf-8")
        return audit


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--tissue", choices=[t[0] for t in TISSUE_ASSETS], default=None)
    args = parser.parse_args()
    CACHE_ROOT.mkdir(parents=True, exist_ok=True)
    targets = [t for t in TISSUE_ASSETS if args.tissue is None or t[0] == args.tissue]
    for tissue, url, path, delete_after in targets:
        print(f"\n=== PROCESSING {tissue} ===", flush=True)
        if not path.exists():
            if url is None:
                raise FileNotFoundError(path)
            download_with_guard(url, path, max_file_gb=5.0)
        audit = verify_and_aggregate_h5ad(tissue, path)
        if delete_after and path.exists():
            path.unlink()
            print(f"Deleted temporary H5AD: {path}", flush=True)

        print(f"{tissue}: shape={audit['obs_shape']}, gate={audit['gate']}, "
              f"10x3v3 donors={len(audit['donors_10x3v3'])}, suspensions={audit['unique_suspension_types']}, "
              f"tissues_in_pub={audit['unique_tissue_in_pub']}")
        for d in audit["donors_10x3v3"]:
            print(f"  {d['donor_id']:6s} | {d['development_stage']:20s} | {d['sex']:6s} | n={d['n_cells_10x3v3']:6d} | "
                  f"comp={d['compartment_counts']} | notes={d['notes_counts']}")


if __name__ == "__main__":
    main()
