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

DONOR_MAP = {
    "5mo_A": "mo5_fish1", "5mo_B": "mo5_fish2",
    "12mo_A": "mo12_fish1", "12mo_B": "mo12_fish2",
    "mo20": "mo20_fish1", "22mo": "mo22_fish1",
    "27mo_A": "mo27_fish1", "27mo_B": "mo27_fish2",
    "27mo_C": "mo27_fish3", "27mo_D": "mo27_fish4",
}


def standard_obs(adata, atlas):
    obs = adata.obs.copy()
    obs["age"] = obs["age"].astype(str).str.extract(r"(\d+)")[0].astype(float)
    obs["donor"] = obs["sub.sample"].astype(str).map(DONOR_MAP)
    if obs["donor"].isna().any():
        raise ValueError("Unmapped zebrafish donor")
    obs["sex_male"] = 1.0
    obs["incomplete_spermatogenesis"] = (obs["sperm.status"].astype(str) == "incomplete").astype(float)
    if atlas == "27mo":
        fine = obs["cell.type.cat"].astype(str)
        broad = fine.map({
            "germ cells": "germ cells", "lymphocytes": "immune",
            "macrophages": "immune", "leukocytes": "immune", "neutrophils": "immune",
            "blood": "blood_and_somatic", "Sertoli, Leydig, SMC-like": "blood_and_somatic",
            "epithelial cells": "blood_and_somatic", "unknown": "unknown",
        })
    else:
        fine = obs["cell_type"].astype(str)
        broad = pd.Series(
            np.where(fine == "immune", "immune",
                     np.where(fine.str.startswith(("spermat", "elongating", "early round",
                                                   "middle round", "late round")),
                              "germ cells", "blood_and_somatic")),
            index=obs.index,
        )
    obs["fine_cell_type"] = fine
    obs["broad_cell_type"] = broad.astype(str)
    return obs


def composition_arm(adatas):
    obs = pd.concat([standard_obs(a, name) for name, a in adatas.items()])
    rows = []
    designs = []
    for analysis, cell_column, selected in (
        ("broad_all_ages", "broad_cell_type", obs),
        ("fine_5_to_22_months", "fine_cell_type", obs[obs["age"] < 27]),
    ):
        selected = selected[selected[cell_column] != "unknown"]
        counts = composition_counts(selected, "donor", cell_column)
        counts = counts.loc[:, counts.sum(axis=0) >= 50]
        metadata = selected.groupby("donor").agg(
            age=("age", "first"), sex_male=("sex_male", "first"),
            incomplete_spermatogenesis=("incomplete_spermatogenesis", "first"),
        )
        for sensitivity, observed_only, adjust_total in (
            ("all_samples", False, False),
            ("observed_records_only", True, False),
            ("total_cells_adjusted", False, True),
        ):
            current = composition_tests(
                counts, metadata, "Zebrafish", "Testis", "GSE275361", "month",
                sensitivity=sensitivity, observed_only=observed_only, adjust_total=adjust_total,
            )
            if analysis == "broad_all_ages" and sensitivity == "all_samples":
                current += composition_tests(
                    counts, metadata, "Zebrafish", "Testis", "GSE275361", "month",
                    sensitivity="spermatogenesis_adjusted",
                    extra_covariates=["incomplete_spermatogenesis"],
                )
            for row in current:
                row["analysis"] = analysis
            rows.extend(current)
        design = (selected.groupby(["donor", "age"], observed=True).size()
                  .rename("n_annotated_cells").reset_index())
        design["analysis"] = analysis
        designs.append(design)
    return rows, designs


def vector_arm(adatas, analysis):
    cell_column = "broad_cell_type" if analysis == "broad" else "fine_cell_type"
    selected_atlases = adatas if analysis == "broad" else {"5-22mo": adatas["5-22mo"]}
    checkpoint = RESULTS / f"extended-zebrafish-vector-donors-{analysis}.csv"
    if checkpoint.exists():
        previous = pd.read_csv(checkpoint)
        rows = previous.to_dict("records")
        completed = set(previous["cell_type"])
    else:
        rows = []
        completed = set()
    obs_by_atlas = {name: standard_obs(a, name) for name, a in selected_atlases.items()}
    cell_types = sorted(set().union(*[
        set(obs.loc[(obs[cell_column] != "unknown"), cell_column]) for obs in obs_by_atlas.values()
    ]))
    for cell_type in cell_types:
        if cell_type in completed:
            continue
        parts = []
        for name, atlas in selected_atlases.items():
            obs = obs_by_atlas[name]
            keep = ~technical_gene_mask(atlas.var_names)
            mask = (obs[cell_column] == cell_type).to_numpy()
            if mask.any():
                part = atlas[mask, keep].to_memory()
                part.obs = standard_obs(part, name)
                parts.append(part)
        if not parts:
            continue
        part = parts[0] if len(parts) == 1 else ad.concat(parts, join="inner")
        part.obs["cell_type"] = cell_type
        young = part[part.obs["age"] == 5].copy()
        before = len(rows)
        if young.obs["donor"].nunique() >= 2 and young.n_obs >= 50:
            reference = build_reference(young.X, young.obs, "cell_type", "donor")
            tissue = "Testis broad compartments" if analysis == "broad" else "Testis fine cell types"
            rows.extend(donor_distances(
                part.X, part.obs, reference, "Zebrafish", tissue, "GSE275361",
                "cell_type", "donor", "age", "sex_male"
            ))
        print(f"Zebrafish {analysis} {cell_type}: {len(rows) - before} fish rows", flush=True)
        pd.DataFrame(rows).to_csv(checkpoint, index=False)
    return rows


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--composition-only", action="store_true")
    parser.add_argument("--vector-only", action="store_true")
    args = parser.parse_args()
    adatas = {
        "5-22mo": ad.read_h5ad(REPO / "data" / "interim" / "zebrafish-testis-5-22mo.h5ad", backed="r"),
        "27mo": ad.read_h5ad(REPO / "data" / "interim" / "zebrafish-testis-27mo.h5ad", backed="r"),
    }
    if not args.vector_only:
        rows, designs = composition_arm(adatas)
        pd.DataFrame(rows).to_csv(RESULTS / "extended-zebrafish-composition-tests.csv", index=False)
        pd.concat(designs, ignore_index=True).to_csv(
            RESULTS / "extended-zebrafish-sample-design.csv", index=False
        )
    if not args.composition_only:
        for analysis in ("broad", "fine"):
            donors = pd.DataFrame(vector_arm(adatas, analysis))
            donors.to_csv(RESULTS / f"extended-zebrafish-vector-donors-{analysis}.csv", index=False)
            if analysis == "broad":
                donors.to_csv(RESULTS / "extended-zebrafish-vector-donors.csv", index=False)
            _, global_results = write_vector_results(
                donors, {"Zebrafish": 22.0}, f"extended-zebrafish-vector-{analysis}",
                {"Zebrafish": "month"},
            )
            print(global_results.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
