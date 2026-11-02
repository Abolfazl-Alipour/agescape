#!/usr/bin/env python3

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.extended_atlas_utils import (
    RESULTS, composition_counts, composition_tests, write_vector_results,
)
from scripts.run_cross_species_expression_vector import (
    build_reference, donor_distances, technical_gene_mask,
)

NORMAL_SAMPLES = {
    "N2D1", "N2D1R", "N2OPD1",
    "N2D6", "N2D6R", "N2D6R2",
    "N2D12", "N2D12R", "N2D12R2",
    "N2D14", "N2D14R", "N2D14R2",
}
REPRODUCTIVE = {"Germline", "Embroynic_cells", "Sperms"}


def standard_obs(adata):
    obs = adata.obs.copy()
    obs["donor"] = obs["orig.ident"].astype(str)
    obs["age"] = obs["donor"].str.extract(r"D(\d+)")[0].astype(float)
    obs["sex_male"] = 0.0
    nextseq = {"N2D6R", "N2D6R2", "N2D12R", "N2D12R2", "N2D14R", "N2D14R2"}
    obs["platform_nextseq"] = obs["donor"].isin(nextseq).astype(float)
    obs["cell_type"] = obs["tissue"].astype(str)
    return obs


def normal_obs(adata):
    obs = standard_obs(adata)
    return obs[(obs["donor"].isin(NORMAL_SAMPLES)) & (obs["singlelet"] == True)
               & (obs["cell_type"] != "singleton")].copy()


def composition_arm(adata):
    obs = normal_obs(adata)
    rows = []
    for analysis, selected in (
        ("all_annotated_cells", obs),
        ("somatic_only", obs[~obs["cell_type"].isin(REPRODUCTIVE)]),
    ):
        counts = composition_counts(selected, "donor", "cell_type")
        counts = counts.loc[:, counts.sum(axis=0) >= 50]
        metadata = selected.groupby("donor").agg(
            age=("age", "first"), sex_male=("sex_male", "first"),
            platform_nextseq=("platform_nextseq", "first"),
        )
        for sensitivity, observed_only, adjust_total, platform_adjusted in (
            ("all_samples", False, False, False),
            ("observed_records_only", True, False, False),
            ("total_nuclei_adjusted", False, True, False),
            ("platform_adjusted", False, False, True),
        ):
            current = composition_tests(
                counts, metadata, "C. elegans", "Whole animal", "Gao_GSE229022", "day",
                sensitivity=sensitivity, observed_only=observed_only, adjust_total=adjust_total,
                extra_covariates=["platform_nextseq"] if platform_adjusted else None,
            )
            for row in current:
                row["analysis"] = analysis
            rows.extend(current)
    design = (obs.groupby(["donor", "age"], observed=True).size()
              .rename("n_annotated_cells").reset_index())
    return rows, design


def vector_arm(adata, reference_age=1):
    obs = normal_obs(adata)
    full_obs = standard_obs(adata)
    keep = ~technical_gene_mask(adata.var_names)
    checkpoint = RESULTS / (
        "extended-worm-vector-donors.csv" if reference_age == 1
        else f"extended-worm-vector-donors-reference{reference_age}.csv"
    )
    if checkpoint.exists():
        previous = pd.read_csv(checkpoint)
        rows = previous.to_dict("records")
        completed = set(previous["cell_type"])
    else:
        rows = []
        completed = set()
    for cell_type in sorted(obs["cell_type"].unique()):
        if cell_type in completed:
            continue
        mask = (
            (full_obs["cell_type"] == cell_type)
            & full_obs["donor"].isin(NORMAL_SAMPLES)
            & (full_obs["singlelet"] == True)
        ).to_numpy()
        part = adata[mask, keep].to_memory()
        part.obs = normal_obs(part)
        young = part[part.obs["age"] == reference_age].copy()
        before = len(rows)
        if young.obs["donor"].nunique() >= 2 and young.n_obs >= 50:
            reference = build_reference(young.X, young.obs, "cell_type", "donor")
            rows.extend(donor_distances(
                part.X, part.obs, reference, "C. elegans", "Whole animal", "Gao_GSE229022",
                "cell_type", "donor", "age", "sex_male"
            ))
        print(f"Gao worm {cell_type}: {len(rows) - before} pool rows", flush=True)
        pd.DataFrame(rows).to_csv(checkpoint, index=False)
    return rows


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--composition-only", action="store_true")
    parser.add_argument("--vector-only", action="store_true")
    parser.add_argument("--reference-age", type=int, default=1)
    args = parser.parse_args()
    adata = ad.read_h5ad(REPO / "data" / "interim" / "worm-gao-integrated.h5ad", backed="r")
    if not args.vector_only:
        rows, design = composition_arm(adata)
        pd.DataFrame(rows).to_csv(RESULTS / "extended-worm-composition-tests.csv", index=False)
        design.to_csv(RESULTS / "extended-worm-sample-design.csv", index=False)
    if not args.composition_only:
        suffix = "" if args.reference_age == 1 else f"-reference{args.reference_age}"
        donors = pd.DataFrame(vector_arm(adata, args.reference_age))
        nextseq = {"N2D6R", "N2D6R2", "N2D12R", "N2D12R2", "N2D14R", "N2D14R2"}
        donors["platform_nextseq"] = donors["donor"].isin(nextseq).astype(float)
        donors.to_csv(RESULTS / f"extended-worm-vector-donors{suffix}.csv", index=False)
        _, global_results = write_vector_results(
            donors, {"C. elegans": 14.0}, f"extended-worm-vector{suffix}", {"C. elegans": "day"}
        )
        somatic = donors[~donors["cell_type"].isin(REPRODUCTIVE)].copy()
        somatic.to_csv(RESULTS / f"extended-worm-vector-donors-somatic{suffix}.csv", index=False)
        write_vector_results(
            somatic, {"C. elegans": 14.0}, f"extended-worm-vector-somatic{suffix}", {"C. elegans": "day"}
        )
        print(global_results.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
