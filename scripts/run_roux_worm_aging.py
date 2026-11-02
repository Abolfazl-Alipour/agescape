#!/usr/bin/env python3

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import statsmodels.api as sm

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.extended_atlas_utils import RESULTS, clr_table
from scripts.run_cross_species_expression_vector import (
    build_reference, donor_distances, technical_gene_mask,
)


def metadata(adata):
    obs = adata.obs.copy()
    obs["donor"] = obs.index.astype(str).str.split("-").str[-1]
    obs["age"] = obs["timepoint"].astype(str).str.removeprefix("d").astype(float)
    obs["sex_male"] = 0.0
    obs["cell_type"] = obs["annotate_name"].astype(str)
    return obs


def descriptive_slopes(data, outcome, species, tissue, cohort):
    rows = []
    for cell_type, part in data.groupby("cell_type"):
        if len(part) < 6:
            continue
        for adjusted in (False, True):
            predictors = ["age"] + (["mean_log_umi"] if adjusted else [])
            X = sm.add_constant(part[predictors], has_constant="add")
            if np.linalg.matrix_rank(X.to_numpy()) < X.shape[1]:
                continue
            fit = sm.OLS(part[outcome], X).fit()
            rows.append({
                "species": species, "tissue": tissue, "cohort": cohort, "cell_type": cell_type,
                "model": "depth_adjusted" if adjusted else "unadjusted",
                "effect_per_day": fit.params["age"], "r2": fit.rsquared,
                "n_libraries": len(part), "inference": "descriptive_only_technical_replicates",
            })
    return rows


def composition_arm(adata):
    obs = metadata(adata)
    counts = (obs.groupby(["donor", "cell_type"], observed=True).size()
              .unstack(fill_value=0).sort_index())
    counts = counts.loc[:, ((counts > 0).sum(axis=0) >= 6) & (counts.sum(axis=0) >= 50)]
    meta = obs.groupby("donor").agg(age=("age", "first"), n_counts=("n_counts", "sum"))
    meta["mean_log_umi"] = np.log(meta["n_counts"])
    clr = clr_table(counts)
    data = clr.stack().rename("clr").reset_index()
    data = data.merge(meta[["age", "mean_log_umi"]], left_on="donor", right_index=True)
    counts.stack().rename("n_cells").reset_index().to_csv(
        RESULTS / "extended-roux-composition-counts.csv", index=False
    )
    return descriptive_slopes(data, "clr", "C. elegans", "Whole animal", "Roux_GSE208154")


def vector_arm(adata):
    obs = metadata(adata)
    keep = ~technical_gene_mask(adata.var["gene_names"].fillna("").astype(str))
    rows = []
    checkpoint = RESULTS / "extended-roux-vector-libraries.csv"
    if checkpoint.exists():
        previous = pd.read_csv(checkpoint)
        rows = previous.to_dict("records")
        completed = set(previous["cell_type"])
    else:
        completed = set()
    for cell_type in sorted(obs["cell_type"].unique()):
        if cell_type in completed:
            continue
        part = adata[(obs["cell_type"] == cell_type).to_numpy(), keep].to_memory()
        part.obs = metadata(part)
        young = part[part.obs["age"] == 1].copy()
        if young.obs["donor"].nunique() < 1 or young.n_obs < 50:
            continue
        reference = build_reference(young.X, young.obs, "cell_type", "donor")
        before = len(rows)
        rows.extend(donor_distances(
            part.X, part.obs, reference, "C. elegans", "Whole animal", "Roux_GSE208154",
            "cell_type", "donor", "age", "sex_male"
        ))
        print(f"Roux {cell_type}: {len(rows) - before} library rows", flush=True)
        pd.DataFrame(rows).to_csv(checkpoint, index=False)
    return rows


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--composition-only", action="store_true")
    parser.add_argument("--vector-only", action="store_true")
    args = parser.parse_args()
    adata = ad.read_h5ad(REPO / "data" / "raw" / "worm-roux" / "ad_worm_aging.h5ad", backed="r")
    if not args.vector_only:
        rows = composition_arm(adata)
        pd.DataFrame(rows).to_csv(RESULTS / "extended-roux-composition-effects.csv", index=False)
    if not args.composition_only:
        donors = pd.DataFrame(vector_arm(adata))
        donors.to_csv(RESULTS / "extended-roux-vector-libraries.csv", index=False)
        rows = descriptive_slopes(
            donors, "profile_hellinger", "C. elegans", "Whole animal", "Roux_GSE208154"
        )
        pd.DataFrame(rows).to_csv(RESULTS / "extended-roux-vector-effects.csv", index=False)


if __name__ == "__main__":
    main()
