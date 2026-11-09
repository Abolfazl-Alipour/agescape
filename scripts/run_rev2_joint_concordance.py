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
    import anndata as ad
    from scripts.run_human_vector_gene_drivers import pseudobulk
    from scripts.run_m9_pipeline import COHORTS

    tests = pd.read_csv(RESULTS / "h3-cross-species-vector-tests.csv")
    eligible = tests[(tests["species"] == "Human")
                     & (tests["metric"] == "profile_hellinger")
                     & (tests["model"] == "depth_adjusted")][
        ["cohort", "tissue", "cell_type"]].drop_duplicates()
    directory = CACHE.parent / "human-h5ads"
    cohort_order = sorted(eligible["cohort"].unique())
    for cohort, selected in eligible.groupby("cohort"):
        cache_path = CACHE / "human" / f"stratum-{cohort_order.index(cohort)}.pkl"
        if cache_path.exists():
            with cache_path.open("rb") as handle:
                yield pickle.load(handle)
            continue
        dataset_id = next(key for key, value in COHORTS.items() if value[0] == cohort)
        path = directory / f"{dataset_id}.h5ad"
        if not path.exists():
            raise FileNotFoundError(f"Missing Census 2025-11-08 source file: {path}")
        atlas = ad.read_h5ad(path, backed="r")
        try:
            columns = ["donor_id", "development_stage", "disease", "cell_type", "sex",
                       "is_primary_data"]
            obs = atlas.obs[columns].copy()
            obs["row"] = np.arange(len(obs))
            obs["age"] = obs["development_stage"].astype(str).str.extract(
                r"(\d+)-year-old", expand=False).astype(float)
            keep = (obs["is_primary_data"].astype(str).str.lower().eq("true")
                    & obs["disease"].astype(str).eq("normal")
                    & obs["cell_type"].astype(str).isin(set(selected["cell_type"]))
                    & obs["age"].between(20, 120))
            obs = obs.loc[keep].copy().sort_values("row")
            rows = obs["row"].to_numpy()
            library = np.empty(len(rows), dtype=np.float64)
            for start in range(0, len(rows), 10000):
                chunk = rows[start:start + 10000]
                library[start:start + len(chunk)] = np.asarray(
                    atlas.X[chunk].sum(axis=1)).ravel()
            obs = obs.loc[library >= 1500].copy()
            donors = sorted(obs["donor_id"].astype(str).unique())
            collected = {}
            genes = None
            for start in range(0, len(donors), 1):
                batch = donors[start:start + 1]
                batch_rows = np.sort(obs.loc[obs["donor_id"].isin(batch), "row"].to_numpy())
                adata = atlas[batch_rows].to_memory()
                adata.uns = {}
                for cell_type in selected["cell_type"]:
                    part = adata[adata.obs["cell_type"].astype(str) == cell_type]
                    if part.n_obs == 0 or not (
                            part.obs.groupby("donor_id", observed=True).size() >= 20).any():
                        continue
                    expression, metadata, gene_names = pseudobulk(part)
                    genes = gene_names
                    if len(metadata):
                        collected.setdefault(cell_type, []).append((expression, metadata))
                print(f"Human {cohort}: donors {start + 1}-{min(start + 1, len(donors))} "
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
            atlas.file.close()


LOADERS = {"Human": human_strata, "Fruit fly": phase1.fly_strata,
           "C. elegans": phase1.worm_strata, "Zebrafish": phase1.zebrafish_strata}


def load_strata(species: str, refresh: bool = False):
    directory = CACHE / SLUGS[species]
    directory.mkdir(parents=True, exist_ok=True)
    if refresh:
        for path in directory.glob("stratum-*.pkl"):
            path.unlink()
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
    return summary, stratum_rows, null


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
    expected = {"Human": 7, "Fruit fly": 2, "C. elegans": 1, "Zebrafish": 1}[args.species]
    rows = []
    nulls = []
    n_donors = 0
    summary = None
    for index, stratum in enumerate(source):
        stratum_summary, stratum_rows, null = analyze_species(
            args.species, [stratum], args.permutations, seed + index)
        rows.extend(stratum_rows)
        nulls.append(null)
        n_donors += stratum_summary["n_donors"]
        summary = stratum_summary
        del stratum, stratum_summary, null
        gc.collect()
    if len(rows) != expected:
        raise ValueError(f"{args.species}: expected {expected} strata, found {len(rows)}")
    weights = np.array([row["n_pairs"] for row in rows], dtype=float)
    observed = float(np.average([row["observed_mean_pearson"] for row in rows], weights=weights))
    null = np.average(np.vstack(nulls), axis=0, weights=weights)
    summary.update({
        "observed_mean_pearson": observed, "n_pairs": int(weights.sum()),
        "n_units": sum(row["n_units"] for row in rows), "n_donors": n_donors,
        "n_strata": len(rows), "joint_null_mean": float(null.mean()),
        "joint_null_sd": float(null.std(ddof=1)),
        "joint_null_ci95_low": float(np.quantile(null, 0.025)),
        "joint_null_ci95_high": float(np.quantile(null, 0.975)),
        "joint_p_greater": float((1 + (null >= observed).sum()) / (len(null) + 1)),
        "random_seed": seed,
    })
    slug = SLUGS[args.species]
    pd.DataFrame([summary]).to_csv(RESULTS / f"reviewer-rev2-joint-concordance-{slug}.csv",
                                   index=False)
    pd.DataFrame(rows).to_csv(RESULTS / f"reviewer-rev2-joint-concordance-{slug}-strata.csv",
                              index=False)
    print(pd.Series(summary).to_string(), flush=True)


if __name__ == "__main__":
    main()
