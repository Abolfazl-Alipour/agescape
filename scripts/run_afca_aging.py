#!/usr/bin/env python3

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.extended_atlas_utils import (
    RESULTS, afca_sample_id, composition_counts, composition_tests, write_vector_results,
)
from scripts.run_cross_species_expression_vector import (
    build_reference, donor_distances, technical_gene_mask,
)


def afca_obs(adata):
    obs = adata.obs.copy()
    obs["sample"] = afca_sample_id(obs.index)
    obs["age"] = obs["age"].astype(float)
    obs["sex_male"] = (obs["sex"].astype(str) == "male").astype(float)
    obs["cell_type"] = obs["afca_annotation"].astype(str)
    return obs


def composition_arm(adata, organ):
    obs = afca_obs(adata)
    obs = obs[(obs["dataset"] == "AFCA") & (obs["cell_type"] != "unannotated")].copy()
    counts = composition_counts(obs, "sample", "cell_type")
    metadata = obs.groupby("sample").agg(age=("age", "first"), sex_male=("sex_male", "first"))
    rows = composition_tests(counts, metadata, "Fruit fly", organ, "AFCA", "day")
    rows += composition_tests(counts, metadata, "Fruit fly", organ, "AFCA", "day",
                              sensitivity="observed_records_only", observed_only=True)
    rows += composition_tests(counts, metadata, "Fruit fly", organ, "AFCA", "day",
                              sensitivity="total_nuclei_adjusted", adjust_total=True)
    design = (obs.groupby(["sample", "age", "sex"], observed=True).size()
              .rename("n_annotated_cells").reset_index())
    return rows, design


def vector_arm(adata, organ):
    obs = afca_obs(adata)
    keep_gene = ~technical_gene_mask(adata.var_names) & ~adata.var_names.astype(str).str.startswith("mt:")
    checkpoint = RESULTS / f"extended-afca-vector-donors-{organ}.csv"
    if checkpoint.exists():
        previous = pd.read_csv(checkpoint)
        rows = previous.to_dict("records")
        completed = set(previous["cell_type"])
    else:
        rows = []
        completed = set()
    cell_types = sorted(set(obs.loc[
        (obs["dataset"] == "AFCA") & (obs["cell_type"] != "unannotated"), "cell_type"
    ]))
    for cell_type in cell_types:
        if cell_type in completed:
            continue
        mask = (
            (obs["dataset"] == "AFCA")
            & (obs["cell_type"] == cell_type)
        ).to_numpy()
        part = adata[mask, keep_gene].to_memory()
        part.obs = afca_obs(part)
        part.X = sp.csr_matrix(part.X)
        part.X.data = np.expm1(part.X.data)
        young = part[part.obs["age"] == 30].copy()
        before = len(rows)
        if young.obs["sample"].nunique() >= 2:
            reference = build_reference(young.X, young.obs, "cell_type", "sample")
            rows.extend(donor_distances(
                part.X, part.obs, reference, "Fruit fly", organ, "AFCA", "cell_type",
                "sample", "age", "sex_male"
            ))
        print(f"AFCA {organ} {cell_type}: {len(rows) - before} sample rows", flush=True)
        pd.DataFrame(rows).to_csv(checkpoint, index=False)
    return rows


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--composition-only", action="store_true")
    parser.add_argument("--vector-only", action="store_true")
    args = parser.parse_args()
    RESULTS.mkdir(exist_ok=True)
    all_composition = []
    all_design = []
    all_donors = []
    for organ in ("head", "body"):
        adata = ad.read_h5ad(REPO / "data" / "raw" / "afca" / f"afca-{organ}.h5ad", backed="r")
        if not args.vector_only:
            rows, design = composition_arm(adata, organ)
            all_composition.extend(rows)
            all_design.append(design)
        if not args.composition_only:
            all_donors.extend(vector_arm(adata, organ))
    if not args.vector_only:
        composition = pd.DataFrame(all_composition)
        composition.to_csv(RESULTS / "extended-afca-composition-tests.csv", index=False)
        pd.concat(all_design, ignore_index=True).to_csv(
            RESULTS / "extended-afca-sample-design.csv", index=False
        )
        primary = composition[composition["sensitivity"] == "all_samples"]
        print(primary.groupby("organ").agg(
            n_cell_types=("cell_type", "nunique"),
            n_significant=("q_bh", lambda x: int((x < 0.05).sum())),
            min_q=("q_bh", "min"),
        ).to_string(), flush=True)
    if not args.composition_only:
        donors = pd.DataFrame(all_donors)
        donors["donor"] = (donors["donor"].str.replace("_head_", "_", regex=False)
                           .str.replace("_body_", "_", regex=False))
        donors.to_csv(RESULTS / "extended-afca-vector-donors.csv", index=False)
        _, global_results = write_vector_results(
            donors, {"Fruit fly": 40.0}, "extended-afca-vector", {"Fruit fly": "day"}
        )
        print(global_results.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
