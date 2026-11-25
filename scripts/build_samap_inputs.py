#!/usr/bin/env python3
"""Build SAMap inputs: per-species mature-only single-cell AnnData.

For each of the six species, loads the same source data and unit labels as
the concordance pipeline, applies the same mature windows and QC filters,
subsamples to at most 3000 cells per unit (seed 20261125), and writes one
AnnData per species to data/samap/inputs/<slug>.h5ad.

Output obs columns (same for every species): species, tissue, cohort, unit,
donor, age (native unit: mouse months, human years, fly days, worm days,
zebrafish months, NMR years), sex_male. X is raw counts (fly counts are
recovered from log1p with expm1, as the concordance pipeline does).
Var index keeps the source gene identifier (mouse/fly/zebrafish symbols,
human ENSG, worm sequence names/symbols, NMR ENSHGLG); var['symbol'] holds
the gene symbol where it differs from the index.

Memory rule: never load a full human cohort file. Pass 1 counts cells per
unit from obs (plus chunked library sums where the pipeline filters on them);
pass 2 samples row indices per unit, then loads exactly the sampled rows.
"""
from __future__ import annotations

import argparse
import gc
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
DATA = REPO.parent / "agescape-revision-data"
OLD_REPO = Path("/home/reza/projects/agescape")
RESULTS = REPO / "results"
OUT = REPO / "data" / "samap" / "inputs"

SEED = 20261125
MAX_PER_UNIT = 3000

AGE_MAP_MONTHS = {"1m": 1.0, "3m": 3.0, "18m": 18.0, "21m": 21.0, "24m": 24.0, "30m": 30.0}

COHORTS = {  # Census dataset UUID -> (cohort, tissue)
    "37a17b78-4864-4a42-b67b-31c00962795a": ("MSSM_DLPFC", "brain"),
    "5c97eeeb-7e52-44b3-b010-b832b1f5424c": ("HBCC_DLPFC", "brain"),
    "0bae7ebf-eb54-46a6-be9a-3461cecefa4c": ("Kidney_LivingDonor", "kidney"),
    "dea717d4-7bc0-4e46-950f-fd7e1cc8df7d": ("Kidney_SCR_v15", "kidney"),
    "eb499fd8-7000-419f-8854-926b9b61b11f": ("Lung_CellRef", "lung"),
    "9f222629-9e39-47d0-b83f-e08d610c7479": ("Lung_Integrated", "lung"),
    "1b350d0a-4535-4879-beb6-1142f3f94947": ("ImmuneAging_BAL", "lung"),
}

MOUSE_TISSUES = ["Fat", "Heart_and_Aorta", "Kidney", "Limb_Muscle", "Lung", "Marrow", "Spleen"]
SLUGS = {"Mouse": "mouse", "Human": "human", "Fruit fly": "fly",
         "C. elegans": "worm", "Zebrafish": "zebrafish", "Naked mole-rat": "nmr"}


def roster():
    level0 = pd.read_csv(RESULTS / "followup-mature-hierarchy-level0.csv")
    out = {}
    for species, group in level0.groupby("species"):
        units = {}
        for row in group.itertuples():
            units.setdefault(row.tissue_cohort, set()).add(row.cell_type)
        out[species] = units
    return out


def subsample_indices(unit_labels, rng):
    """Per-unit uniform subsample of row positions, at most MAX_PER_UNIT each."""
    labels = pd.Series(unit_labels)
    keep = []
    for unit, idx in labels.groupby(labels).groups.items():
        idx = np.asarray(sorted(idx))
        if len(idx) > MAX_PER_UNIT:
            idx = np.sort(rng.choice(idx, size=MAX_PER_UNIT, replace=False))
        keep.append(idx)
    return np.sort(np.concatenate(keep)) if keep else np.array([], dtype=int)


def standardize(adata, species, tissue, cohort):
    """Rename obs to the common schema and drop everything else."""
    import anndata as ad
    obs = pd.DataFrame({
        "species": species,
        "tissue": tissue,
        "cohort": cohort,
        "unit": adata.obs["unit"].astype(str).to_numpy(),
        "donor": adata.obs["donor"].astype(str).to_numpy(),
        "age": adata.obs["age"].astype(float).to_numpy(),
        "sex_male": adata.obs["sex_male"].astype(float).to_numpy(),
    }, index=adata.obs_names)
    var = pd.DataFrame(index=adata.var_names)
    if "symbol" in adata.var.columns:
        var["symbol"] = adata.var["symbol"].astype(str).to_numpy()
    out = ad.AnnData(X=adata.X, obs=obs, var=var)
    return out


def write_species(chunks, species, manifest_rows):
    import anndata as ad
    adata = ad.concat(chunks, join="outer", merge="same") if len(chunks) > 1 else chunks[0]
    adata.obs_names_make_unique()
    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)
    adata.X = adata.X.astype(np.int32)
    slug = SLUGS[species]
    path = OUT / f"{slug}.h5ad"
    adata.write_h5ad(path, compression="gzip")
    print(f"{species}: wrote {path} ({adata.n_obs} cells x {adata.n_vars} genes, "
          f"{adata.obs['unit'].nunique()} units)", flush=True)
    for (tissue, cohort, unit), group in adata.obs.groupby(
            ["tissue", "cohort", "unit"], observed=True):
        manifest_rows.append({
            "species": species, "tissue": tissue, "cohort": cohort, "unit": unit,
            "n_cells_kept": int(len(group)),
            "n_donors": int(group["donor"].nunique()),
        })
    return adata.n_obs


# ---------------------------------------------------------------- mouse
def build_mouse(roster_map):
    import anndata as ad
    rng = np.random.default_rng(SEED)
    chunks, kept = [], 0
    for tissue in MOUSE_TISSUES:
        key = f"{tissue} | TMS"
        wanted = roster_map["Mouse"].get(key, set())
        if not wanted:
            continue
        path = REPO / "data" / "interim" / f"m1-qc-cleaned-{tissue}.h5ad"
        atlas = ad.read_h5ad(path, backed="r")
        obs = atlas.obs
        age = obs["age"].astype(str).map(AGE_MAP_MONTHS).to_numpy()
        keep = ((~obs["predicted_doublet"].astype(bool).to_numpy())
                & (obs["n_counts"].to_numpy() >= 1500)
                & (obs["cell_ontology_class"].astype(str).to_numpy() != "unknown")
                & (age >= 18)
                & (obs["cell_ontology_class"].astype(str).isin(wanted).to_numpy()))
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
        chunks.append(standardize(part, "Mouse", tissue, "TMS"))
        kept += len(rows)
        print(f"  Mouse {tissue}: {len(rows)} cells kept", flush=True)
        del part
        gc.collect()
    return chunks, kept


# ---------------------------------------------------------------- human
def build_human(roster_map):
    import anndata as ad
    rng = np.random.default_rng(SEED + 1)
    chunks, kept = [], 0
    for uuid, (cohort, tissue) in sorted(COHORTS.items(), key=lambda kv: kv[1][0]):
        key = f"{tissue} | {cohort}"
        wanted = roster_map["Human"].get(key, set())
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
        # chunked library sums (raw_sum >= 1500), as the concordance loader does
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
        chunks.append(standardize(part, "Human", tissue, cohort))
        kept += len(rows)
        print(f"  Human {cohort}: {len(rows)} cells kept", flush=True)
        del part
        gc.collect()
    return chunks, kept


# ---------------------------------------------------------------- fly
def build_fly(roster_map):
    import anndata as ad
    from scripts.run_afca_aging import afca_obs
    rng = np.random.default_rng(SEED + 2)
    chunks, kept = [], 0
    for tissue in ("head", "body"):
        key = f"{tissue} | AFCA"
        wanted = roster_map["Fruit fly"].get(key, set())
        if not wanted:
            continue
        path = REPO / "data" / "raw" / "afca" / f"afca-{tissue}.h5ad"
        atlas = ad.read_h5ad(path, backed="r")
        obs = afca_obs(atlas)
        keep = ((obs["dataset"].astype(str).to_numpy() == "AFCA")
                & (obs["cell_type"].astype(str).to_numpy() != "unannotated")
                & (obs["age"].to_numpy() >= 30)
                & (obs["cell_type"].astype(str).isin(wanted).to_numpy()))
        rows = np.where(keep)[0]
        unit_labels = obs["cell_type"].astype(str).to_numpy()[rows]
        sel = subsample_indices(unit_labels, rng)
        rows = rows[sel]
        part = atlas[rows].to_memory()
        atlas.file.close()
        # AFCA X is log1p-normalized; recover counts exactly as the pipeline does
        x = part.X
        if sp.issparse(x):
            x = x.tocsr(copy=True)
            x.data = np.expm1(x.data)
        else:
            x = np.expm1(x)
        part.X = sp.csr_matrix(np.round(x).astype(np.int32))
        part.obs = obs.loc[obs.index[rows]].copy()
        part.obs["unit"] = part.obs["cell_type"].astype(str)
        part.obs["donor"] = part.obs["sample"].astype(str)
        part.obs["age"] = part.obs["age"].astype(float)
        part.obs["sex_male"] = part.obs["sex_male"].astype(float)
        chunks.append(standardize(part, "Fruit fly", tissue, "AFCA"))
        kept += len(rows)
        print(f"  Fly {tissue}: {len(rows)} cells kept", flush=True)
        del part
        gc.collect()
    return chunks, kept


# ---------------------------------------------------------------- worm
def build_worm(roster_map):
    import anndata as ad
    from scripts.run_gao_worm_aging import normal_obs
    from scripts.run_cross_species_expression_vector import technical_gene_mask
    rng = np.random.default_rng(SEED + 3)
    wanted = roster_map["C. elegans"].get("Whole animal | Gao_GSE229022", set())
    atlas = ad.read_h5ad(REPO / "data" / "interim" / "worm-gao-integrated.h5ad", backed="r")
    obs = normal_obs(atlas)
    keep = (obs["age"].to_numpy() >= 6) & (obs["cell_type"].astype(str).isin(wanted).to_numpy())
    idx = np.where(keep)[0]
    unit_labels = obs["cell_type"].astype(str).to_numpy()[idx]
    sel = subsample_indices(unit_labels, rng)
    idx = idx[sel]
    rows = obs.index[idx]
    part = atlas[rows].to_memory()
    atlas.file.close()
    gene_keep = ~technical_gene_mask(part.var_names.astype(str))
    part = part[:, gene_keep].copy()
    part.obs = obs.loc[rows].copy()
    part.obs["unit"] = part.obs["cell_type"].astype(str)
    part.obs["donor"] = part.obs["donor"].astype(str)
    part.obs["age"] = part.obs["age"].astype(float)
    part.obs["sex_male"] = part.obs["sex_male"].astype(float)
    return [standardize(part, "C. elegans", "Whole animal", "Gao_GSE229022")], len(idx)


# ---------------------------------------------------------------- zebrafish
def build_zebrafish(roster_map):
    import anndata as ad
    from scripts.run_zebrafish_testis_aging import standard_obs
    from scripts.run_cross_species_expression_vector import technical_gene_mask
    rng = np.random.default_rng(SEED + 4)
    wanted = roster_map["Zebrafish"].get("Testis fine cell types | GSE275361", set())
    atlas = ad.read_h5ad(REPO / "data" / "interim" / "zebrafish-testis-5-22mo.h5ad", backed="r")
    obs = standard_obs(atlas, "5-22mo")
    keep = (obs["fine_cell_type"].astype(str).to_numpy() != "unknown") & \
        (obs["fine_cell_type"].astype(str).isin(wanted).to_numpy())
    idx = np.where(keep)[0]
    unit_labels = obs["fine_cell_type"].astype(str).to_numpy()[idx]
    sel = subsample_indices(unit_labels, rng)
    idx = idx[sel]
    part = atlas[idx].to_memory()
    atlas.file.close()
    gene_keep = ~technical_gene_mask(part.var_names.astype(str))
    part = part[:, gene_keep].copy()
    part.obs = obs.iloc[idx].copy()
    part.obs["unit"] = part.obs["fine_cell_type"].astype(str)
    part.obs["donor"] = part.obs["donor"].astype(str)
    part.obs["age"] = part.obs["age"].astype(float)
    part.obs["sex_male"] = part.obs["sex_male"].astype(float)
    return [standardize(part, "Zebrafish", "Testis fine cell types", "GSE275361")], len(idx)


# ---------------------------------------------------------------- NMR
def build_nmr(roster_map):
    import anndata as ad
    import scanpy as sc
    from scripts.run_cross_species_expression_vector import (
        technical_gene_mask, NMR_AGES)
    rng = np.random.default_rng(SEED + 5)
    wanted = roster_map["Naked mole-rat"].get("Skin | GSE152886", set())
    raw = pd.read_csv(OLD_REPO / "data" / "raw" / "GSE152886_nmr_skin_raw.csv.gz",
                      sep="\t", index_col=0)
    obs = pd.DataFrame({"barcode": raw.columns})
    obs["sample"] = obs["barcode"].str.split("_").str[0]
    obs = obs[obs["sample"] != "2V"].copy()
    obs["age"] = obs["sample"].map(NMR_AGES)
    X = sp.csr_matrix(raw[obs["barcode"]].T.values)
    adata = sc.AnnData(X=X, obs=obs.set_index("barcode"))
    adata.var_names = raw.index.astype(str)
    adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
    adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    adata = adata[(adata.obs["n_counts"] >= 2000) & (adata.obs["n_genes"] >= 700)].copy()
    # replicate the annotation clustering exactly (same params, same seed)
    raw_counts = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    markers = {
        "immune": ["ENSHGLG00000017278", "ENSHGLG00000015756", "ENSHGLG00000005722"],
        "basal": ["ENSHGLG00000000726", "ENSHGLG00000016277", "ENSHGLG00000017682"],
        "corneum": ["ENSHGLG00000010577", "ENSHGLG00000001162"],
    }
    marker_genes = {gene for genes in markers.values() for gene in genes}
    retain = adata.var["highly_variable"] | adata.var_names.isin(marker_genes)
    sub = adata[:, retain].copy()
    sc.pp.scale(sub, max_value=10)
    sc.pp.pca(sub, n_comps=50, random_state=42)
    sc.pp.neighbors(sub, n_neighbors=15, n_pcs=30, random_state=42)
    sc.tl.leiden(sub, resolution=0.4, random_state=42)
    annotations = []
    for cluster in sorted(sub.obs["leiden"].unique()):
        cells = sub[sub.obs["leiden"] == cluster]
        scores = {}
        for name, genes in markers.items():
            present = [g for g in genes if g in cells.var_names]
            scores[name] = float(np.asarray(cells[:, present].X.mean()).ravel()[0]) if present else 0.0
        if scores["immune"] > 0.1 and scores["immune"] >= max(scores["basal"], scores["corneum"]):
            compartment = "immune"
        elif scores["corneum"] > scores["basal"]:
            compartment = "corneum"
        elif scores["basal"] > 0:
            compartment = "basal"
        else:
            compartment = "intermediate"
        annotations.append({"cluster": cluster, "compartment": compartment})
    annotation_map = pd.DataFrame(annotations).set_index("cluster")["compartment"]
    adata.obs["cell_type"] = sub.obs["leiden"].astype(str).map(annotation_map)
    adata.X = raw_counts
    symbols = pd.read_csv(
        OLD_REPO / "data" / "interim" / "nmr-ensg2symbol.tsv", sep="\t", skiprows=1,
        names=["ensembl_id", "symbol"]).set_index("ensembl_id")["symbol"]
    gene_symbols = adata.var_names.map(symbols).fillna("").astype(str)
    keep_genes = ~technical_gene_mask(gene_symbols)
    adata = adata[:, keep_genes].copy()
    adata.var["symbol"] = gene_symbols.to_numpy()[keep_genes]
    # keep only roster units; subsample
    keep_cells = adata.obs["cell_type"].astype(str).isin(wanted).to_numpy()
    adata = adata[keep_cells].copy()
    unit_labels = adata.obs["cell_type"].astype(str).to_numpy()
    sel = subsample_indices(unit_labels, rng)
    adata = adata[sel].copy()
    adata.obs["unit"] = adata.obs["cell_type"].astype(str)
    adata.obs["donor"] = adata.obs["sample"].astype(str)
    adata.obs["age"] = adata.obs["age"].astype(float)
    adata.obs["sex_male"] = 0.0
    return [standardize(adata, "Naked mole-rat", "Skin", "GSE152886")], adata.n_obs


BUILDERS = {"Mouse": build_mouse, "Human": build_human, "Fruit fly": build_fly,
            "C. elegans": build_worm, "Zebrafish": build_zebrafish,
            "Naked mole-rat": build_nmr}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", choices=list(BUILDERS), nargs="+", default=list(BUILDERS))
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    roster_map = roster()
    manifest_rows = []
    for species in args.species:
        print(f"Building {species} ...", flush=True)
        chunks, n_kept = BUILDERS[species](roster_map)
        write_species(chunks, species, manifest_rows)
        del chunks
        gc.collect()
    manifest = pd.DataFrame(manifest_rows)
    manifest_path = OUT / "manifest.csv"
    if manifest_path.exists():
        old = pd.read_csv(manifest_path)
        manifest = pd.concat([old[~old["species"].isin(args.species)], manifest],
                             ignore_index=True)
    manifest.to_csv(manifest_path, index=False)
    print(f"Wrote {manifest_path} ({len(manifest)} units)", flush=True)


if __name__ == "__main__":
    main()
