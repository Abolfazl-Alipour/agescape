#!/usr/bin/env python3
"""Test within-cell transcriptional entropy and off-program expression in mouse lung."""

from pathlib import Path
import re

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
FIGURES = REPO / "figures"
AGE_MAP = {"1m": 1.0, "3m": 3.0, "18m": 18.0, "21m": 21.0, "24m": 24.0, "30m": 30.0}
TARGET_CELL_TYPES = [
    "fibroblast of lung",
    "bronchial smooth muscle cell",
    "alveolar macrophage",
    "classical monocyte",
    "intermediate monocyte",
    "B cell",
    "NK cell",
    "CD8-positive, alpha-beta T cell",
    "vein endothelial cell",
]
SEX_GENES = {"Xist", "Tsix", "Ddx3y", "Eif2s3y", "Uty", "Kdm5d"}
AMBIENT_PRONE_GENES = {"Jchain", "Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2"}
IMMEDIATE_EARLY_GENES = {
    "Fos", "Fosb", "Jun", "Junb", "Jund", "Egr1", "Egr2", "Egr3",
    "Ier2", "Ier3", "Zfp36", "Atf3", "Dusp1",
}


def excluded_gene_mask(gene_names):
    genes = pd.Index(gene_names).astype(str)
    technical = np.array([
        bool(re.match(r"^(mt-|Rp[ls]|Mrp[ls])", gene)) or gene in SEX_GENES
        for gene in genes
    ])
    stringent = np.array([
        bool(re.match(r"^(Ig[hk]|Hb[ab])", gene))
        or gene in AMBIENT_PRONE_GENES | IMMEDIATE_EARLY_GENES
        for gene in genes
    ])
    return technical, stringent


def shannon_entropy(matrix):
    matrix = sp.csr_matrix(matrix, dtype=float)
    totals = np.asarray(matrix.sum(axis=1)).ravel()
    weighted_logs = matrix.copy()
    weighted_logs.data = weighted_logs.data * np.log(weighted_logs.data)
    moments = np.asarray(weighted_logs.sum(axis=1)).ravel()
    entropy = np.full(matrix.shape[0], np.nan)
    valid = totals > 0
    entropy[valid] = np.log(totals[valid]) - moments[valid] / totals[valid]
    return entropy


def reference_stats(matrix, obs):
    stats = {}
    for cell_type in obs["cell_ontology_class"].unique():
        ct_mask = (obs["cell_ontology_class"] == cell_type).to_numpy()
        ct_matrix = matrix[ct_mask]
        ct_obs = obs.loc[ct_mask]
        per_donor = {}
        for donor in ct_obs["mouse.id"].unique():
            donor_mask = (ct_obs["mouse.id"] == donor).to_numpy()
            donor_matrix = ct_matrix[donor_mask]
            per_donor[str(donor)] = {
                "counts": np.asarray(donor_matrix.sum(axis=0)).ravel(),
                "detected": np.asarray((donor_matrix > 0).sum(axis=0)).ravel(),
                "library": float(donor_matrix.sum()),
                "n_cells": int(donor_matrix.shape[0]),
            }
        stats[str(cell_type)] = {
            "counts": np.asarray(ct_matrix.sum(axis=0)).ravel(),
            "detected": np.asarray((ct_matrix > 0).sum(axis=0)).ravel(),
            "library": float(ct_matrix.sum()),
            "n_cells": int(ct_matrix.shape[0]),
            "donors": per_donor,
        }
    return stats


def without_donor(stats, donor):
    result = {key: stats[key].copy() if isinstance(stats[key], np.ndarray) else stats[key]
              for key in ("counts", "detected", "library", "n_cells")}
    if donor in stats["donors"]:
        donor_stats = stats["donors"][donor]
        result["counts"] -= donor_stats["counts"]
        result["detected"] -= donor_stats["detected"]
        result["library"] -= donor_stats["library"]
        result["n_cells"] -= donor_stats["n_cells"]
    return result


def offprogram_markers(stats, target_cell_type, excluded_donor, gene_names, stringent=False):
    target = without_donor(stats[target_cell_type], excluded_donor)
    if target["n_cells"] < 50 or target["library"] <= 0:
        return np.array([], dtype=int), {}
    target_detection = target["detected"] / target["n_cells"]
    target_cp10k = target["counts"] / target["library"] * 1e4
    other_detection = np.zeros(len(gene_names))
    other_cp10k = np.zeros(len(gene_names))
    other_source = np.array([""] * len(gene_names), dtype=object)
    for cell_type in stats:
        if cell_type == target_cell_type:
            continue
        other = without_donor(stats[cell_type], excluded_donor)
        if other["n_cells"] < 50 or other["library"] <= 0:
            continue
        detection = other["detected"] / other["n_cells"]
        cp10k = other["counts"] / other["library"] * 1e4
        replace = cp10k > other_cp10k
        other_detection[replace] = detection[replace]
        other_cp10k[replace] = cp10k[replace]
        other_source[replace] = cell_type
    technical, stringent_mask = excluded_gene_mask(gene_names)
    marker_mask = (
        (target_detection <= 0.05)
        & (other_detection >= 0.20)
        & (other_cp10k >= 20)
        & (other_cp10k >= 5 * (target_cp10k + 1))
        & ~technical
    )
    if stringent:
        marker_mask &= ~stringent_mask
    indices = np.flatnonzero(marker_mask)
    provenance = {
        gene_names[index]: other_source[index]
        for index in indices
    }
    return indices, provenance


def donor_regression(records, metric, analysis, cell_type):
    if len(records) < 6:
        return None
    X = sm.add_constant(records[["age_months", "sex_male"]].astype(float), has_constant="add")
    model = "age_plus_sex"
    if np.linalg.matrix_rank(X.to_numpy()) < X.shape[1]:
        X = sm.add_constant(records[["age_months"]].astype(float), has_constant="add")
        model = "age_only_rank_deficient_sex"
    fit = sm.OLS(records[metric].astype(float), X).fit()
    ci = fit.conf_int().loc["age_months"]
    return {
        "cell_type": cell_type,
        "analysis": analysis,
        "metric": metric,
        "model": model,
        "n_donors": len(records),
        "n_cells": int(records["n_cells"].sum()),
        "slope_per_month": fit.params["age_months"],
        "se": fit.bse["age_months"],
        "ci95_low": ci.iloc[0],
        "ci95_high": ci.iloc[1],
        "pvalue": fit.pvalues["age_months"],
        "r2": fit.rsquared,
    }


def make_figure(results):
    primary = results[results["analysis"] == "all_cells_all_donor_ids"]
    metrics = [("entropy", "Shannon entropy slope per month"),
               ("offprogram_fraction", "Off-program UMI-fraction slope per month")]
    fig, axes = plt.subplots(1, 2, figsize=(15, 6), sharey=True)
    cell_types = TARGET_CELL_TYPES[::-1]
    positions = np.arange(len(cell_types))
    for ax, (metric, xlabel) in zip(axes, metrics):
        part = primary[primary["metric"] == metric].set_index("cell_type")
        for position, cell_type in zip(positions, cell_types):
            if cell_type not in part.index:
                continue
            row = part.loc[cell_type]
            color = "#0072B2" if row["q_bh"] >= 0.05 else "#D55E00"
            ax.plot([row["ci95_low"], row["ci95_high"]], [position, position], color=color, lw=2)
            ax.scatter(row["slope_per_month"], position, color=color, s=45, zorder=3)
        ax.axvline(0, color="black", lw=0.8, ls="--")
        ax.set_xlabel(xlabel)
        ax.grid(axis="x", alpha=0.25)
    axes[0].set_yticks(positions, cell_types)
    fig.suptitle("Age associations with within-cell transcriptional disorder in mouse lung")
    fig.tight_layout()
    fig.savefig(FIGURES / "h3-transcriptional-disorder.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    RESULTS.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    adata = ad.read_h5ad(REPO / "data" / "interim" / "m1-qc-cleaned.h5ad", backed="r")
    lung_mask = (
        (adata.obs["tissue"] == "Lung")
        & ~adata.obs["predicted_doublet"]
        & (adata.obs["n_counts"] >= 1500)
        & (adata.obs["cell_ontology_class"].astype(str) != "unknown")
    )
    lung = adata[lung_mask].to_memory()
    lung.obs["age_months"] = lung.obs["age"].astype(str).map(AGE_MAP)
    lung.obs["sex_male"] = (lung.obs["sex"].astype(str).str.lower() == "male").astype(float)
    lung.obs["mouse.id"] = lung.obs["mouse.id"].astype(str)
    print(f"Loaded {lung.n_obs} eligible lung singlet cells", flush=True)

    young_mask = lung.obs["age_months"] <= 3
    young = lung[young_mask].copy()
    stats = reference_stats(young.X.tocsr(), young.obs)
    print(f"Built young reference profiles for {len(stats)} cell types", flush=True)

    downsampled = lung.copy()
    sc.pp.downsample_counts(downsampled, counts_per_cell=1500, replace=False, random_state=42)
    downsampled.X = downsampled.X.tocsr()
    totals = np.asarray(downsampled.X.sum(axis=1)).ravel()
    assert np.all(totals == 1500)
    gene_names = downsampled.var_names.astype(str).to_numpy()
    technical, stringent = excluded_gene_mask(gene_names)
    entropy_matrix = downsampled.X[:, ~technical]
    lung.obs["entropy"] = shannon_entropy(entropy_matrix)
    stringent_entropy_matrix = downsampled.X[:, ~(technical | stringent)]
    lung.obs["entropy_stringent"] = shannon_entropy(stringent_entropy_matrix)
    lung.X = downsampled.X
    print("Completed without-replacement UMI rarefaction and entropy calculation", flush=True)

    rng = np.random.default_rng(20261101)
    donor_rows = []
    marker_rows = []
    for cell_type in TARGET_CELL_TYPES:
        ct_mask = (lung.obs["cell_ontology_class"] == cell_type).to_numpy()
        ct_obs = lung.obs.loc[ct_mask]
        ct_matrix = lung.X[ct_mask]
        for donor in ct_obs["mouse.id"].unique():
            donor_mask = (ct_obs["mouse.id"] == donor).to_numpy()
            if donor_mask.sum() < 20:
                continue
            donor_obs = ct_obs.loc[donor_mask]
            donor_matrix = ct_matrix[donor_mask]
            markers, provenance = offprogram_markers(
                stats, cell_type, str(donor), gene_names, stringent=False
            )
            stringent_markers, _ = offprogram_markers(
                stats, cell_type, str(donor), gene_names, stringent=True
            )
            offprogram = (
                np.asarray(donor_matrix[:, markers].sum(axis=1)).ravel() / 1500
                if len(markers) else np.full(donor_matrix.shape[0], np.nan)
            )
            offprogram_stringent = (
                np.asarray(donor_matrix[:, stringent_markers].sum(axis=1)).ravel() / 1500
                if len(stringent_markers) else np.full(donor_matrix.shape[0], np.nan)
            )
            detected = (
                np.asarray((donor_matrix[:, markers] > 0).sum(axis=1)).ravel()
                if len(markers) else np.full(donor_matrix.shape[0], np.nan)
            )
            for source_cell_type in sorted(set(provenance.values())):
                marker_rows.append({
                    "cell_type": cell_type,
                    "excluded_donor": donor,
                    "source_cell_type": source_cell_type,
                    "n_marker_genes": sum(source == source_cell_type for source in provenance.values()),
                    "n_all_marker_genes": len(markers),
                    "n_stringent_marker_genes": len(stringent_markers),
                })
            base = {
                "cell_type": cell_type,
                "donor": donor,
                "age_months": float(donor_obs["age_months"].iloc[0]),
                "sex_male": float(donor_obs["sex_male"].iloc[0]),
                "n_offprogram_genes": len(markers),
                "n_stringent_offprogram_genes": len(stringent_markers),
            }
            values = {
                "entropy": donor_obs["entropy"].to_numpy(),
                "entropy_stringent": donor_obs["entropy_stringent"].to_numpy(),
                "offprogram_fraction": offprogram,
                "offprogram_fraction_stringent": offprogram_stringent,
                "offprogram_detected_genes": detected,
            }
            for sampling in ("all_cells", "balanced_20_cells"):
                selected = np.arange(donor_matrix.shape[0])
                if sampling == "balanced_20_cells":
                    selected = rng.choice(selected, size=20, replace=False)
                row = base | {"sampling": sampling, "n_cells": len(selected)}
                for metric, metric_values in values.items():
                    selected_values = metric_values[selected]
                    row[metric] = (
                        float(np.nanmean(selected_values))
                        if np.isfinite(selected_values).any() else np.nan
                    )
                donor_rows.append(row)
        print(f"Completed {cell_type}", flush=True)

    donors = pd.DataFrame(donor_rows)
    donors.to_csv(RESULTS / "h3-transcriptional-disorder-donors.csv", index=False)
    pd.DataFrame(marker_rows).to_csv(RESULTS / "h3-offprogram-marker-provenance.csv", index=False)
    pooled = set(pd.read_csv(RESULTS / "m1-animal-tissue-design.csv")
                 .query("is_pooled_id")["mouse.id"].astype(str))
    results = []
    metrics = [
        "entropy", "entropy_stringent", "offprogram_fraction",
        "offprogram_fraction_stringent", "offprogram_detected_genes",
    ]
    for cell_type in TARGET_CELL_TYPES:
        ct_donors = donors[donors["cell_type"] == cell_type]
        for sampling in ("all_cells", "balanced_20_cells"):
            sampled = ct_donors[ct_donors["sampling"] == sampling]
            for donor_set, selected in (("all_donor_ids", sampled),
                                        ("exclude_pooled_ids", sampled[~sampled["donor"].isin(pooled)])):
                analysis = f"{sampling}_{donor_set}"
                for metric in metrics:
                    records = selected.dropna(subset=[metric])
                    row = donor_regression(records, metric, analysis, cell_type)
                    if row:
                        results.append(row)
    results = pd.DataFrame(results)
    results["q_bh"] = np.nan
    primary_metrics = ["entropy", "offprogram_fraction"]
    for (analysis, metric), indices in results.groupby(["analysis", "metric"]).groups.items():
        qvalues = multipletests(results.loc[indices, "pvalue"], method="fdr_bh")[1]
        results.loc[indices, "q_bh"] = qvalues
    results["endpoint_role"] = np.where(
        results["metric"].isin(primary_metrics), "primary", "sensitivity"
    )
    results.to_csv(RESULTS / "h3-transcriptional-disorder.csv", index=False)
    make_figure(results)
    primary = results[(results["analysis"] == "all_cells_all_donor_ids")
                      & results["metric"].isin(primary_metrics)]
    print(primary[["cell_type", "metric", "n_donors", "slope_per_month",
                   "ci95_low", "ci95_high", "pvalue", "q_bh"]].to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
