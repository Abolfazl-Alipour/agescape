#!/usr/bin/env python3
"""Joint donor-level age permutation for human, fly, worm, and zebrafish.

One age assignment is applied to every cell type contributed by a donor or
pool. Permutations stay within cohort, sex, and platform strata. The primary
statistic is mean pairwise Pearson correlation of donor-pseudobulk age slopes,
using the Rev2 detection filter. Species summaries weight strata by pair count.
"""

from __future__ import annotations

import argparse
import gc
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from scripts import run_reviewer_phase1_gene_concordance as phase1  # noqa: E402
from scripts import run_rev2_concordance as rev2  # noqa: E402

RESULTS = REPO / "results"
CACHE = REPO.parent / f"{REPO.name}-data" / "concordance-cache"
SPECIES = ["Human", "Fruit fly", "C. elegans", "Zebrafish"]
SLUGS = {"Human": "human", "Fruit fly": "fly", "C. elegans": "worm", "Zebrafish": "zebrafish"}


def human_strata():
    import cellxgene_census
    from scripts.run_human_vector_gene_drivers import pseudobulk
    from scripts.run_m9_pipeline import COHORTS
    from scripts.run_rev2_transport_null import tiledb_config_from_env

    tests = pd.read_csv(RESULTS / "h3-cross-species-vector-tests.csv")
    eligible = tests[(tests["species"] == "Human")
                     & (tests["metric"] == "profile_hellinger")
                     & (tests["model"] == "depth_adjusted")][
        ["cohort", "tissue", "cell_type"]].drop_duplicates()
    census = cellxgene_census.open_soma(
        census_version="2025-11-08", tiledb_config=tiledb_config_from_env())
    try:
        for cohort, selected in eligible.groupby("cohort"):
            dataset_id = next(key for key, value in COHORTS.items() if value[0] == cohort)
            cell_filter = ",".join(f"'{value}'" for value in selected["cell_type"])
            base_filter = (f"dataset_id == '{dataset_id}' and is_primary_data == True "
                           f"and disease == 'normal' and raw_sum >= 1500 "
                           f"and cell_type in [{cell_filter}]")
            obs = cellxgene_census.get_obs(
                census, "Homo sapiens", value_filter=base_filter,
                column_names=["donor_id", "development_stage"])
            obs["age"] = obs["development_stage"].astype(str).str.extract(
                r"(\d+)-year-old").astype(float)
            donors = sorted(obs.loc[obs["age"].between(20, 120), "donor_id"].astype(str).unique())
            collected = {}
            genes = None
            for start in range(0, len(donors), 5):
                batch = donors[start:start + 5]
                donor_filter = ",".join(f"'{donor}'" for donor in batch)
                adata = cellxgene_census.get_anndata(
                    census, "Homo sapiens",
                    obs_value_filter=f"{base_filter} and donor_id in [{donor_filter}]",
                    obs_column_names=["donor_id", "development_stage", "sex", "cell_type"],
                    var_column_names=["feature_name"])
                for cell_type in selected["cell_type"]:
                    part = adata[adata.obs["cell_type"].astype(str) == cell_type]
                    if part.n_obs == 0 or not (
                            part.obs.groupby("donor_id", observed=True).size() >= 20).any():
                        continue
                    expression, metadata, gene_names = pseudobulk(part)
                    genes = gene_names
                    if len(metadata):
                        collected.setdefault(cell_type, []).append((expression, metadata))
                print(f"Human {cohort}: donors {start + 1}-{min(start + 5, len(donors))} "
                      f"of {len(donors)}", flush=True)
                del adata
                gc.collect()
            units = []
            for row in selected.itertuples():
                pieces = collected.get(row.cell_type, [])
                if not pieces:
                    continue
                expression = np.vstack([piece[0] for piece in pieces])
                metadata = pd.concat([piece[1] for piece in pieces], ignore_index=True)
                if len(metadata) >= 6:
                    units.append({"species": "Human", "tissue": row.tissue,
                                  "cohort": cohort, "cell_type": row.cell_type,
                                  "expression": expression, "metadata": metadata})
            yield selected["tissue"].iloc[0], cohort, genes, units
    finally:
        census.close()


LOADERS = {"Human": human_strata, "Fruit fly": phase1.fly_strata,
           "C. elegans": phase1.worm_strata, "Zebrafish": phase1.zebrafish_strata}


def load_strata(species: str, refresh: bool = False):
    directory = CACHE / SLUGS[species]
    directory.mkdir(parents=True, exist_ok=True)
    for index, (tissue, cohort, genes, units) in enumerate(LOADERS[species]()):
        path = directory / f"stratum-{index}.pkl"
        if path.exists() and not refresh:
            with path.open("rb") as handle:
                yield pickle.load(handle)
            continue
        payload = (tissue, cohort, genes, units)
        with path.open("wb") as handle:
            pickle.dump(payload, handle, protocol=pickle.HIGHEST_PROTOCOL)
        yield payload
        gc.collect()


def cached_strata(species: str):
    directory = CACHE / SLUGS[species]
    paths = sorted(directory.glob("stratum-*.pkl"),
                   key=lambda path: int(path.stem.split("-")[-1]))
    if not paths:
        raise FileNotFoundError(f"No cached strata for {species}; run without --analysis-only")
    for path in paths:
        with path.open("rb") as handle:
            yield pickle.load(handle)


def detection_mask(units, n_genes):
    fractions = np.vstack([(unit["expression"] > 0).mean(axis=0) for unit in units])
    return (fractions >= 0.5).sum(axis=0) >= np.ceil(0.5 * len(units)), n_genes


def donor_pool(strata):
    records = {}
    for stratum in strata:
        cohort = stratum["cohort"]
        for unit in stratum["units"]:
            for row in unit["metadata"].itertuples():
                key = f"{cohort}|{row.donor}"
                values = (float(row.age), float(row.sex_male),
                          float(getattr(row, "platform_nextseq", 0.0)))
                if key in records:
                    assert records[key] == values, f"Inconsistent metadata for {key}"
                records[key] = values
    donors = sorted(records)
    return donors, np.array([records[d][0] for d in donors]), np.array(
        [records[d][1:] for d in donors])


def permute_donor_ages(ages, strata, rng):
    permuted = ages.copy()
    frame = pd.DataFrame(strata, columns=["sex_male", "platform_nextseq"])
    frame["cohort"] = [key.split("|", 1)[0] for key in CURRENT_DONORS]
    group_columns = ["cohort"] + [c for c in ("sex_male", "platform_nextseq")
                                  if frame[c].nunique() > 1]
    for _, indices in frame.groupby(group_columns, sort=False).groups.items():
        positions = np.asarray(indices, dtype=int)
        permuted[positions] = rng.permutation(permuted[positions])
    return permuted


def analyze_species(species, strata, permutations, seed):
    global CURRENT_DONORS
    prepared = []
    for stratum_index, (tissue, cohort, genes, units) in enumerate(strata):
        units = [unit for unit in units if len(unit["metadata"]) >= 6]
        if len(units) < 2:
            continue
        mask, _ = detection_mask(units, len(genes))
        prepared.append({"stratum_index": stratum_index, "tissue": tissue, "cohort": cohort,
                         "genes": genes, "units": units, "mask": mask})
    if not prepared:
        raise ValueError(f"No evaluable strata for {species}")
    donors, donor_age, donor_strata = donor_pool(prepared)
    CURRENT_DONORS = donors
    donor_index = {donor: index for index, donor in enumerate(donors)}
    for stratum in prepared:
        for unit in stratum["units"]:
            unit["expr64"] = np.asarray(unit["expression"], dtype=np.float64)
            unit["donor_idx"] = np.array([
                donor_index[f"{stratum['cohort']}|{donor}"]
                for donor in unit["metadata"]["donor"]])
    rng = np.random.default_rng(seed)
    nulls = {stratum["stratum_index"]: [] for stratum in prepared}
    stratum_rows = []
    observed_values = {}
    for iteration in range(permutations + 1):
        if iteration == 0:
            ages = donor_age
        else:
            ages = permute_donor_ages(donor_age, donor_strata, rng)
        for stratum in prepared:
            slopes = []
            for unit in stratum["units"]:
                age = ages[unit["donor_idx"]]
                slopes.append(rev2.fit_age_slope_only(
                    unit["expr64"], rev2.unit_design(unit["metadata"], age=age)))
            values = rev2.pairwise_pearson(np.vstack(slopes)[:, stratum["mask"]])
            value = float(np.mean(values))
            if iteration == 0:
                observed_values[stratum["stratum_index"]] = value
            else:
                nulls[stratum["stratum_index"]].append(value)
        if iteration % 100 == 0:
            print(f"{species}: permutation {iteration}/{permutations}", flush=True)
    for stratum in prepared:
        index = stratum["stratum_index"]
        observed = observed_values[index]
        null = np.asarray(nulls[index])
        stratum_rows.append({
            "species": species, "tissue": stratum["tissue"], "cohort": stratum["cohort"],
            "n_units": len(stratum["units"]),
            "n_pairs": len(stratum["units"]) * (len(stratum["units"]) - 1) // 2,
            "n_donors_in_stratum": len({donor for unit in stratum["units"]
                                        for donor in unit["metadata"]["donor"]}),
            "n_genes": int(stratum["mask"].sum()),
            "observed_mean_pearson": observed,
            "joint_null_mean": float(null.mean()),
            "joint_null_sd": float(null.std(ddof=1)),
            "joint_null_ci95_low": float(np.quantile(null, 0.025)),
            "joint_null_ci95_high": float(np.quantile(null, 0.975)),
            "joint_p_greater": float((1 + (null >= observed).sum()) / (len(null) + 1)),
            "n_permutations": permutations,
        })
    weights = np.array([row["n_pairs"] for row in stratum_rows], dtype=float)
    observed = float(np.average([row["observed_mean_pearson"] for row in stratum_rows],
                                weights=weights))
    null = np.average(np.vstack([nulls[stratum["stratum_index"]]
                                 for stratum in prepared]), axis=0, weights=weights)
    summary = {
        "species": species, "scope": "within_tissue_or_cohort",
        "gene_set": "primary_detection_filtered", "observed_mean_pearson": observed,
        "n_pairs": int(weights.sum()), "n_units": sum(row["n_units"] for row in stratum_rows),
        "n_donors": len(donors), "n_strata": len(stratum_rows),
        "detection_rule": "detected in >=50% of donor pseudobulks in >=50% of stratum units",
        "joint_null_mean": float(null.mean()), "joint_null_sd": float(null.std(ddof=1)),
        "joint_null_ci95_low": float(np.quantile(null, 0.025)),
        "joint_null_ci95_high": float(np.quantile(null, 0.975)),
        "joint_p_greater": float((1 + (null >= observed).sum()) / (len(null) + 1)),
        "n_permutations_joint": permutations, "random_seed": seed,
        "permutation_scheme": "one age assignment per donor shared across cell types; "
                              "within cohort, sex, and platform strata",
        "status": "redone_joint",
    }
    return summary, stratum_rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", choices=SPECIES, required=True)
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--refresh-cache", action="store_true")
    parser.add_argument("--analysis-only", action="store_true")
    args = parser.parse_args()
    seed = args.seed if args.seed is not None else 20261108 + sum(map(ord, args.species))
    source = cached_strata(args.species) if args.analysis_only else load_strata(
        args.species, refresh=args.refresh_cache)
    strata = list(source)
    expected = {"Human": 7, "Fruit fly": 2, "C. elegans": 1, "Zebrafish": 1}[args.species]
    if len(strata) != expected:
        raise ValueError(f"{args.species}: expected {expected} strata, found {len(strata)}")
    summary, rows = analyze_species(args.species, strata, args.permutations, seed)
    slug = SLUGS[args.species]
    pd.DataFrame([summary]).to_csv(RESULTS / f"reviewer-rev2-joint-concordance-{slug}.csv",
                                   index=False)
    pd.DataFrame(rows).to_csv(RESULTS / f"reviewer-rev2-joint-concordance-{slug}-strata.csv",
                              index=False)
    print(pd.Series(summary).to_string(), flush=True)


if __name__ == "__main__":
    main()
