#!/usr/bin/env python3
"""Reproduce De Man KL entropy and Agescape Shannon entropy in IPF controls."""

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
DATA = REPO.parent / f"{REPO.name}-data" / "de-man"
RESULTS = REPO / "results"
TARGETS = {"ATI": "AT1", "ATII": "AT2", "VE_Capillary_A": "Aerocyte", "VE_Capillary_B": "gCap"}
SEEDS = [20261113 + i for i in range(20)]


def load_target_matrices():
    cache = DATA / "target_matrices"
    cache.mkdir(exist_ok=True)
    metadata = pd.read_csv(DATA / "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz", sep="\t")
    barcodes = pd.read_csv(DATA / "GSE136831_AllCells.cellBarcodes.txt.gz", sep="\t", header=None)[0].to_numpy()
    assert len(barcodes) == len(metadata)
    metadata["matrix_column"] = np.arange(len(metadata))
    metadata = metadata[metadata["Disease_Identity"] == "Control"].copy()
    gene_info = pd.read_csv(DATA / "GSE136831_AllCells.GeneIDs.txt.gz", sep="\t")
    mito_genes = gene_info["HGNC_EnsemblAlt_GeneID"].astype(str).str.upper().str.startswith("MT-").to_numpy()
    if all((cache / f"{target}.npz").exists() for target in TARGETS) and (cache / "selected_metadata.csv").exists():
        selected = pd.read_csv(cache / "selected_metadata.csv")
        metadata = metadata.merge(selected, on="matrix_column", how="left")
        return metadata, {target: sp.load_npz(cache / f"{target}.npz") for target in TARGETS}

    selected = metadata[metadata["Manuscript_Identity"].isin(TARGETS)].copy()
    selected["local_column"] = np.arange(len(selected))
    local_to_global = selected["matrix_column"].to_numpy()
    entries = {target: ([], [], []) for target in TARGETS}
    totals = np.zeros(len(selected), dtype=np.float64)
    mito = np.zeros(len(selected), dtype=np.float64)
    selected_type = selected["Manuscript_Identity"].to_numpy()
    target_columns = selected.groupby("Manuscript_Identity").cumcount().to_numpy()
    with gzip.open(DATA / "GSE136831_RawCounts_Sparse.mtx.gz", "rt") as handle:
        header = handle.readline()
        dimensions = handle.readline()
        n_genes, n_cells, _ = map(int, dimensions.split())
        assert n_cells == len(barcodes)
        chunks = pd.read_csv(handle, sep=" ", names=["gene", "cell", "count"], chunksize=5_000_000)
        for chunk in chunks:
            mask = np.isin(chunk["cell"].to_numpy(), local_to_global)
            if not mask.any():
                continue
            part = chunk.loc[mask]
            local_cells = np.searchsorted(local_to_global, part["cell"].to_numpy())
            genes = part["gene"].to_numpy()
            counts = part["count"].to_numpy()
            totals += np.bincount(local_cells, weights=counts, minlength=len(selected))
            is_mito = mito_genes[genes - 1]
            if is_mito.any():
                mito += np.bincount(local_cells[is_mito], weights=counts[is_mito], minlength=len(selected))
            for index, target in enumerate(TARGETS):
                type_mask = selected_type[local_cells] == target
                if type_mask.any():
                    entries[target][0].append(genes[type_mask])
                    entries[target][1].append(target_columns[local_cells][type_mask])
                    entries[target][2].append(counts[type_mask])
    for target, (gene_parts, cell_parts, count_parts) in entries.items():
        genes = np.concatenate(gene_parts)
        cells = np.concatenate(cell_parts)
        counts = np.concatenate(count_parts)
        n_target_cells = int((selected_type == target).sum())
        matrix = sp.coo_matrix((counts, (cells, genes - 1)), shape=(n_target_cells, n_genes)).tocsr()
        sp.save_npz(cache / f"{target}.npz", matrix)
    selected["actual_umi"] = totals
    selected["mito_umi"] = mito
    selected["qc_pass"] = selected["nGene"].between(750, 7500) & (selected["actual_umi"] > 0) & (
        selected["mito_umi"] / selected["actual_umi"] < 0.10)
    selected[["matrix_column", "actual_umi", "mito_umi", "qc_pass"]].to_csv(
        cache / "selected_metadata.csv", index=False)
    metadata = metadata.merge(selected[["matrix_column", "actual_umi", "mito_umi", "qc_pass"]],
                              on="matrix_column", how="left")
    return metadata, {target: sp.load_npz(cache / f"{target}.npz") for target in TARGETS}


def binomial_downsample(matrix, rng):
    matrix = matrix.tocsr().astype(np.int64)
    totals = np.asarray(matrix.sum(axis=1)).ravel()
    factor = totals.min() / totals
    result = matrix.copy().astype(float)
    for row in range(matrix.shape[0]):
        start, end = matrix.indptr[row:row + 2]
        result.data[start:end] = rng.binomial(matrix.data[start:end], factor[row])
    return result


def multinomial_downsample(matrix, depth, rng):
    matrix = matrix.tocsr().astype(np.int64)
    totals = np.asarray(matrix.sum(axis=1)).ravel()
    assert (totals >= depth).all()
    result = matrix.copy().astype(float)
    for row in range(matrix.shape[0]):
        start, end = matrix.indptr[row:row + 2]
        result.data[start:end] = rng.multinomial(depth, matrix.data[start:end] / totals[row])
    return result


def kl_uniform(matrix, sqrt=False):
    matrix = matrix.tocsr().astype(float)
    if sqrt:
        matrix.data = np.sqrt(matrix.data)
    totals = np.asarray(matrix.sum(axis=1)).ravel()
    keep = totals > 0
    proportions = matrix[keep].copy()
    proportions.data /= np.repeat(totals[keep], np.diff(proportions.indptr))
    n_genes = matrix.shape[1]
    positive = proportions.data > 0
    proportions.data[positive] *= np.log(proportions.data[positive] * n_genes)
    return np.asarray(proportions.sum(axis=1)).ravel()


def shannon(matrix):
    matrix = matrix.tocsr().astype(float)
    totals = np.asarray(matrix.sum(axis=1)).ravel()
    keep = totals > 0
    proportions = matrix[keep].copy()
    proportions.data /= np.repeat(totals[keep], np.diff(proportions.indptr))
    positive = proportions.data > 0
    proportions.data[positive] *= np.log(proportions.data[positive])
    return -np.asarray(proportions.sum(axis=1)).ravel()


def model_rows(per_cell, seed, method):
    rows = []
    for cell_type, part in per_cell.groupby("cell_type"):
        donor = part.groupby("sample_id").agg(
            entropy=("entropy", "mean"), age=("age", "first"), sex=("sex", "first"),
            actual_umi=("actual_umi", "mean"), n_cells=("entropy", "size")).reset_index()
        if len(donor) < 6:
            continue
        donor["log_umi"] = np.log(donor["actual_umi"])
        donor["male"] = (donor["sex"].str.lower() == "male").astype(float)
        X = sm.add_constant(donor[["age", "male", "log_umi"]].astype(float), has_constant="add")
        fit = sm.OLS(donor["entropy"], X).fit()
        ci = fit.conf_int().loc["age"]
        rows.append({"seed": seed, "method": method, "cell_type": cell_type,
                     "n_donors": len(donor), "n_cells": int(donor["n_cells"].sum()),
                     "age_sd": donor["age"].std(), "slope_per_year": fit.params["age"],
                     "se": fit.bse["age"], "ci95_low": ci.iloc[0], "ci95_high": ci.iloc[1],
                     "p": fit.pvalues["age"]})
    return rows


def main():
    metadata, matrices = load_target_matrices()
    subjects = pd.read_csv(DATA / "41467_2026_68810_MOESM3_ESM.csv").rename(
        columns={"Sample ID": "sample_id", "Sex": "sex", "Age": "age"})
    metadata = metadata.merge(subjects[["sample_id", "sex", "age"]],
                              left_on="Subject_Identity", right_on="sample_id", how="left")
    assert metadata.loc[metadata["Manuscript_Identity"].isin(TARGETS), "qc_pass"].notna().all()
    per_cell = []
    gene_info = pd.read_csv(DATA / "GSE136831_AllCells.GeneIDs.txt.gz", sep="\t")
    symbols = gene_info["HGNC_EnsemblAlt_GeneID"].astype(str).str.upper()
    technical = (symbols.str.match(r"^(MT-|RP[LS]|MRP[LS])")
                 | symbols.isin(["XIST", "TSIX", "DDX3Y", "EIF2S3Y", "UTY", "KDM5D"])).to_numpy()
    for target, label in TARGETS.items():
        obs = metadata[metadata["Manuscript_Identity"] == target].copy()
        assert len(obs) == matrices[target].shape[0]
        keep = obs["qc_pass"].to_numpy(dtype=bool)
        obs = obs.loc[keep].copy()
        matrix = matrices[target][keep]
        shannon_matrix = matrix[:, ~technical]
        matrix = matrix[:, np.asarray(matrix.sum(axis=0)).ravel() > 0]
        shannon_matrix = shannon_matrix[:, np.asarray(shannon_matrix.sum(axis=0)).ravel() > 0]
        obs["local_column"] = np.arange(len(obs))
        for seed in SEEDS:
            rng = np.random.default_rng(seed)
            thinned = binomial_downsample(matrix, rng)
            young = obs["age"].to_numpy() <= 60
            n_sample = min(young.sum(), (~young).sum())
            if n_sample < 10:
                continue
            old_idx = rng.choice(np.flatnonzero(~young), size=n_sample, replace=False)
            young_idx = rng.choice(np.flatnonzero(young), size=n_sample, replace=False)
            selected = np.concatenate([old_idx, young_idx])
            kl = kl_uniform(thinned[selected], sqrt=True)
            for index, value in zip(selected, kl):
                row = obs.iloc[index]
                per_cell.append({"seed": seed, "method": "published_KL_binomial", "cell_type": label,
                                 "sample_id": row["sample_id"], "age": row["age"], "sex": row["sex"],
                                 "actual_umi": row["actual_umi"], "entropy": value})

            full_shannon = shannon(shannon_matrix)
            depth_totals = np.asarray(shannon_matrix.sum(axis=1)).ravel()
            fixed_idx = np.flatnonzero(depth_totals >= 1500)
            fixed_shannon = shannon(multinomial_downsample(shannon_matrix[fixed_idx], 1500, rng))
            for method, indices, values in (("shannon_full", np.arange(len(obs)), full_shannon),
                                            ("shannon_1500", fixed_idx, fixed_shannon)):
                for local_index, value in zip(indices, values):
                    row = obs.iloc[local_index]
                    per_cell.append({"seed": seed, "method": method, "cell_type": label,
                                     "sample_id": row["sample_id"], "age": row["age"], "sex": row["sex"],
                                     "actual_umi": row["actual_umi"], "entropy": value})
        print(f"Finished {label}", flush=True)
    per_cell = pd.DataFrame(per_cell)
    per_cell.to_csv(RESULTS / "followup-de-man-entropy-per-cell.csv.gz", index=False,
                    compression="gzip")
    rows = []
    for (seed, method), part in per_cell.groupby(["seed", "method"]):
        rows.extend(model_rows(part, seed, method))
    models = pd.DataFrame(rows)
    primary = models[models["seed"] == SEEDS[0]].copy()
    primary["q_bh"] = np.nan
    for _, indices in primary.groupby("method").groups.items():
        primary.loc[indices, "q_bh"] = multipletests(primary.loc[indices, "p"], method="fdr_bh")[1]
    models = models.merge(primary[["seed", "method", "cell_type", "q_bh"]],
                          on=["seed", "method", "cell_type"], how="left")
    models.to_csv(RESULTS / "followup-de-man-entropy-models.csv", index=False)
    print(primary.to_string(index=False))


if __name__ == "__main__":
    main()
