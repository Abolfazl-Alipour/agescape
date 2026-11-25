#!/usr/bin/env python3
"""Retry input build: append the 11 missing units to the Stage 1 SAMap inputs.

Stage 1 (commit 9aeadde) recovered only 5 of 28 name-matched pairs; 10 pairs
were untestable because their clusters were absent from the SAMap input
(the input roster was the 311-unit concordance Level 0 set). This script
adds exactly the 11 units those 10 pairs need and nothing else:

  Human, M9 cohorts (8 units, same QC as build_samap_inputs.build_human:
  primary data, disease normal, age 30-120, library sum >= 1500):
    Kidney_LivingDonor B cell, Kidney_SCR_v15 B cell,
    Kidney_LivingDonor kidney distal convoluted tubule epithelial cell,
    ImmuneAging_BAL B cell, Lung_CellRef B cell,
    Lung_Integrated T cell, Lung_CellRef classical monocyte,
    Lung_CellRef non-classical monocyte (4 cells, 4 donors after QC).
  Human, Yazar blood (2 units, same mature window and QC, plus the Yazar
  base mask: disease normal, assay 10x 3' v2, suspension cell, doublet
  label "peripheral blood mononuclear cell" excluded):
    hematopoietic precursor cell (1,718 cells, 670 donors at 30+),
    naive B cell (61,914 cells, 931 donors at 30+; subsampled to 3,000).
  Mouse, TMS lung (1 unit, same QC as build_samap_inputs.build_mouse:
  not doublet, n_counts >= 1500, not unknown, age >= 18 months):
    T cell (225 cells, 10 donors at 18m+).

The Stage 1 input files are NOT modified. The existing cells are carried
over verbatim (no resampling), so Stage 1 results stay comparable; only the
new units are subsampled, with a separate RNG (RETRY_SEED). Outputs:
  data/samap/inputs/human-retry.h5ad, mouse-retry.h5ad,
  data/samap/inputs/manifest-retry.csv (added units only).

Zebrafish and naked mole rat are NOT added (dropped by Reza's order
Nov 25 01:08 PT). No match-rule or gate change is made here.
"""
from __future__ import annotations

import gc
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from scripts.build_samap_inputs import (  # noqa: E402
    AGE_MAP_MONTHS,
    COHORTS,
    MAX_PER_UNIT,
    subsample_indices,
    standardize,
)

DATA = REPO.parent / "agescape-revision-data"
IN_DIR = REPO / "data" / "samap" / "inputs"
YAZAR_H5AD = DATA / "tmp_organs" / "hcata" / "yazar_1e44db10.h5ad"
YAZAR_DOUBLET_LABEL = "peripheral blood mononuclear cell"

RETRY_SEED = 20261125 + 100
MAX_PER_UNIT_RETRY = MAX_PER_UNIT  # 3000, same cap as Stage 1

# (cohort, unit) for the 8 missing M9 human units
HUMAN_EXTRA = {
    "Kidney_LivingDonor": {"B cell", "kidney distal convoluted tubule epithelial cell"},
    "Kidney_SCR_v15": {"B cell"},
    "ImmuneAging_BAL": {"B cell"},
    "Lung_CellRef": {"B cell", "classical monocyte", "non-classical monocyte"},
    "Lung_Integrated": {"T cell"},
}
YAZAR_EXTRA = {"hematopoietic precursor cell", "naive B cell"}


def _to_int32(part):
    """Round float counts to int32 in place (sparse-safe)."""
    x = part.X.tocsr(copy=False) if sp.issparse(part.X) else sp.csr_matrix(part.X)
    if x.dtype != np.int32:
        x.data = np.round(x.data).astype(np.int32)
        x = x.astype(np.int32)
    part.X = x
    return part


def build_human_extra(rng):
    """The 8 missing M9 units, mirroring build_samap_inputs.build_human."""
    import anndata as ad
    chunks, kept = [], 0
    for uuid, (cohort, tissue) in sorted(COHORTS.items(), key=lambda kv: kv[1][0]):
        wanted = HUMAN_EXTRA.get(cohort, set())
        if not wanted:
            continue
        path = DATA / "human-h5ads" / f"{uuid}.h5ad"
        atlas = ad.read_h5ad(path, backed="r")
        obs = atlas.obs
        age = obs["development_stage"].astype(str).str.extract(
            r"(\d+)-year-old", expand=False).astype(float).to_numpy()
        base = ((obs["is_primary_data"].astype(str).str.lower().to_numpy() == "true")
                & (obs["disease"].astype(str).to_numpy() == "normal")
                & (obs["cell_type"].astype(str).isin(wanted).to_numpy())
                & (age >= 30) & (age <= 120))
        rows = np.where(base)[0]
        lib = np.empty(len(rows), dtype=np.float64)
        for start in range(0, len(rows), 10000):
            chunk = rows[start:start + 10000]
            lib[start:start + len(chunk)] = np.asarray(atlas.X[chunk].sum(axis=1)).ravel()
        rows = rows[lib >= 1500]
        unit_labels = obs["cell_type"].astype(str).to_numpy()[rows]
        sel = subsample_indices(unit_labels, rng)
        rows = rows[sel]
        part = atlas[rows].to_memory()
        atlas.file.close()
        part.var["symbol"] = part.var["feature_name"].astype(str)
        part.obs["unit"] = part.obs["cell_type"].astype(str)
        part.obs["donor"] = part.obs["donor_id"].astype(str)
        part.obs["age"] = part.obs["development_stage"].astype(str).str.extract(
            r"(\d+)-year-old", expand=False).astype(float)
        part.obs["sex_male"] = (part.obs["sex"].astype(str).str.lower() == "male").astype(float)
        part = _to_int32(part)
        chunks.append(standardize(part, "Human", tissue, cohort))
        kept += len(rows)
        print(f"  Human {cohort}: {len(rows)} cells kept "
              f"({sorted(set(unit_labels[sel]))})", flush=True)
        del part
        gc.collect()
    return chunks, kept


def build_yazar(rng):
    """The 2 missing Yazar blood units (same mature window and QC)."""
    import anndata as ad
    atlas = ad.read_h5ad(YAZAR_H5AD, backed="r")
    obs = atlas.obs
    age = obs["development_stage"].astype(str).str.extract(
        r"(\d+)-year-old", expand=False).astype(float).to_numpy()
    base = ((obs["disease"].astype(str).to_numpy() == "normal")
            & (obs["assay"].astype(str).to_numpy() == "10x 3' v2")
            & (obs["suspension_type"].astype(str).to_numpy() == "cell")
            & (obs["development_stage"].astype(str).str.endswith("-year-old stage"))
            & (obs["is_primary_data"].astype(str).str.lower().to_numpy() == "true")
            & (obs["cell_type"].astype(str).isin(YAZAR_EXTRA).to_numpy())
            & (obs["cell_type"].astype(str).to_numpy() != YAZAR_DOUBLET_LABEL)
            & (age >= 30) & (age <= 120))
    rows = np.where(base)[0]
    print(f"  Yazar pre-library-filter rows: {len(rows)}", flush=True)
    lib = np.empty(len(rows), dtype=np.float64)
    for start in range(0, len(rows), 10000):
        chunk = rows[start:start + 10000]
        lib[start:start + len(chunk)] = np.asarray(atlas.X[chunk].sum(axis=1)).ravel()
    rows = rows[lib >= 1500]
    print(f"  Yazar post-library-filter rows: {len(rows)}", flush=True)
    unit_labels = obs["cell_type"].astype(str).to_numpy()[rows]
    sel = subsample_indices(unit_labels, rng)
    rows = rows[sel]
    part = atlas[rows].to_memory()
    atlas.file.close()
    part.var["symbol"] = part.var["feature_name"].astype(str)
    part.obs["unit"] = part.obs["cell_type"].astype(str)
    part.obs["donor"] = part.obs["donor_id"].astype(str)
    part.obs["age"] = part.obs["development_stage"].astype(str).str.extract(
        r"(\d+)-year-old", expand=False).astype(float)
    part.obs["sex_male"] = (part.obs["sex"].astype(str).str.lower() == "male").astype(float)
    part = _to_int32(part)
    print(f"  Yazar: {len(rows)} cells kept "
          f"({sorted(set(unit_labels[sel]))})", flush=True)
    return [standardize(part, "Human", "blood", "Yazar_Blood")], len(rows)


def build_mouse_extra(rng):
    """Mouse lung T cell, mirroring build_samap_inputs.build_mouse."""
    import anndata as ad
    atlas = ad.read_h5ad(REPO / "data" / "interim" / "m1-qc-cleaned-Lung.h5ad", backed="r")
    obs = atlas.obs
    age = obs["age"].astype(str).map(AGE_MAP_MONTHS).to_numpy()
    keep = ((~obs["predicted_doublet"].astype(bool).to_numpy())
            & (obs["n_counts"].to_numpy() >= 1500)
            & (obs["cell_ontology_class"].astype(str).to_numpy() != "unknown")
            & (age >= 18)
            & (obs["cell_ontology_class"].astype(str).to_numpy() == "T cell"))
    rows = np.where(keep)[0]
    unit_labels = obs["cell_ontology_class"].astype(str).to_numpy()[rows]
    sel = subsample_indices(unit_labels, rng)
    rows = rows[sel]
    part = atlas[rows].to_memory()
    atlas.file.close()
    part.obs["unit"] = part.obs["cell_ontology_class"].astype(str)
    part.obs["donor"] = part.obs["mouse.id"].astype(str)
    part.obs["age"] = part.obs["age"].astype(str).map(AGE_MAP_MONTHS).astype(float)
    part.obs["sex_male"] = (part.obs["sex"].astype(str).str.lower() == "male").astype(float)
    part = _to_int32(part)
    print(f"  Mouse Lung T cell: {len(rows)} cells kept", flush=True)
    return [standardize(part, "Mouse", "Lung", "TMS")], len(rows)


def append_species(slug, old_path, new_chunks, manifest_rows, species):
    import anndata as ad
    old = ad.read_h5ad(old_path)
    n_old = old.n_obs
    # verify old clusters are carried over verbatim (counts per cluster)
    old_counts = old.obs.groupby(["tissue", "cohort", "unit"], observed=True).size()
    chunks = [old] + new_chunks
    out = ad.concat(chunks, join="outer", merge="same") if len(chunks) > 1 else chunks[0]
    out.obs_names_make_unique()
    if not sp.issparse(out.X):
        out.X = sp.csr_matrix(out.X)
    out.X = out.X.astype(np.int32)
    # new clusters only, for the manifest
    new_clusters = set()
    for ch in new_chunks:
        for (t, c, u) in ch.obs[["tissue", "cohort", "unit"]].drop_duplicates().itertuples(index=False):
            new_clusters.add((t, c, u))
    for (tissue, cohort, unit), group in out.obs.groupby(
            ["tissue", "cohort", "unit"], observed=True):
        if (tissue, cohort, unit) in new_clusters:
            manifest_rows.append({
                "species": species, "tissue": tissue, "cohort": cohort,
                "unit": unit, "n_cells_kept": int(len(group)),
                "n_donors": int(group["donor"].nunique()),
            })
    # verbatim check: old clusters keep their exact cell counts
    new_counts = out.obs.groupby(["tissue", "cohort", "unit"], observed=True).size()
    for key, n in old_counts.items():
        assert new_counts[key] == n, f"old cluster {key} changed: {n} -> {new_counts[key]}"
    out_path = IN_DIR / f"{slug}-retry.h5ad"
    out.write_h5ad(out_path, compression="gzip")
    print(f"{species}: wrote {out_path} ({out.n_obs} cells x {out.n_vars} genes, "
          f"{out.obs['unit'].nunique()} units; {n_old} carried over, "
          f"{out.n_obs - n_old} added)", flush=True)
    return out.n_obs


def main():
    rng_h = np.random.default_rng(RETRY_SEED + 1)
    rng_m = np.random.default_rng(RETRY_SEED)
    manifest_rows = []

    print("Building human extra units ...", flush=True)
    h_chunks, _ = build_human_extra(rng_h)
    print("Building Yazar units ...", flush=True)
    y_chunks, _ = build_yazar(rng_h)
    print("Appending to human input ...", flush=True)
    append_species("human", IN_DIR / "human.h5ad", h_chunks + y_chunks,
                   manifest_rows, "Human")

    print("Building mouse extra units ...", flush=True)
    m_chunks, _ = build_mouse_extra(rng_m)
    print("Appending to mouse input ...", flush=True)
    append_species("mouse", IN_DIR / "mouse.h5ad", m_chunks,
                   manifest_rows, "Mouse")

    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(IN_DIR / "manifest-retry.csv", index=False)
    print(f"Wrote {IN_DIR / 'manifest-retry.csv'} ({len(manifest)} added units)",
          flush=True)
    print(manifest.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
