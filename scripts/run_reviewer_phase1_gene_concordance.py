#!/usr/bin/env python3

import argparse
import itertools
from pathlib import Path
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import statsmodels.api as sm
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_cross_species_expression_vector import technical_gene_mask
from src.reviewer_stats import pairwise_pearson, vectorized_ols

RESULTS = REPO / "results"
INTERIM = REPO / "data" / "interim"
AGE_MAP = {"1m": 1.0, "3m": 3.0, "18m": 18.0, "21m": 21.0, "24m": 24.0, "30m": 30.0}


def aggregate_unit(matrix, obs, cell_type, species, tissue, cohort, cell_column, donor_column,
                   recovered_normalized=False):
    selected = (obs[cell_column].astype(str) == cell_type).to_numpy()
    part_obs = obs[selected].copy()
    part = matrix[selected]
    if sp.issparse(part):
        part = part.tocsr()
    if recovered_normalized:
        part = part.copy()
        if sp.issparse(part):
            part.data = np.expm1(part.data)
        else:
            part = np.expm1(part)
    counts = []
    metadata = []
    donor_ids = part_obs[donor_column].astype(str)
    for donor in sorted(donor_ids.unique()):
        mask = (donor_ids == donor).to_numpy()
        if mask.sum() < 20:
            continue
        donor_matrix = part[mask]
        count = np.asarray(donor_matrix.sum(axis=0)).ravel()
        if count.sum() <= 0:
            continue
        totals = np.asarray(donor_matrix.sum(axis=1)).ravel()
        metadata.append({
            "donor": donor, "age": float(part_obs.loc[mask, "age"].iloc[0]),
            "sex_male": float(part_obs.loc[mask, "sex_male"].iloc[0]),
            "mean_log_umi": float(np.log(np.maximum(totals, 1)).mean()),
            "n_cells": int(mask.sum()),
        })
        counts.append(count)
    if len(counts) < 5:
        return None
    counts = np.vstack(counts)
    expression = np.log1p(counts / counts.sum(axis=1, keepdims=True) * 1e4)
    return {
        "species": species, "tissue": tissue, "cohort": cohort, "cell_type": cell_type,
        "expression": expression, "metadata": pd.DataFrame(metadata),
    }


def mouse_strata():
    atlas = ad.read_h5ad(INTERIM / "m1-qc-cleaned.h5ad", backed="r")
    for tissue in sorted(atlas.obs["tissue"].unique()):
        mask = ((atlas.obs["tissue"] == tissue) & ~atlas.obs["predicted_doublet"]
                & (atlas.obs["n_counts"] >= 1500)
                & (atlas.obs["cell_ontology_class"].astype(str) != "unknown"))
        part = atlas[mask].to_memory()
        part.obs["age"] = part.obs["age"].astype(str).map(AGE_MAP)
        part.obs["sex_male"] = (part.obs["sex"].astype(str).str.lower() == "male").astype(float)
        keep = ~technical_gene_mask(part.var_names)
        part = part[:, keep]
        units = [aggregate_unit(part.X, part.obs, cell_type, "Mouse", tissue, "TMS",
                                "cell_ontology_class", "mouse.id")
                 for cell_type in sorted(part.obs["cell_ontology_class"].unique())]
        yield tissue, "TMS", part.var_names.to_numpy(), [unit for unit in units if unit]


def fly_strata():
    from scripts.run_afca_aging import afca_obs
    for tissue in ("head", "body"):
        atlas = ad.read_h5ad(REPO / "data" / "raw" / "afca" / f"afca-{tissue}.h5ad", backed="r")
        obs = afca_obs(atlas)
        mask = ((obs["dataset"] == "AFCA") & (obs["cell_type"] != "unannotated")
                & (obs["age"] >= 30))
        part = atlas[mask.to_numpy()].to_memory()
        part.obs = obs.loc[mask].copy()
        units = [aggregate_unit(part.X, part.obs, cell_type, "Fruit fly", tissue, "AFCA",
                                "cell_type", "sample", recovered_normalized=True)
                 for cell_type in sorted(part.obs["cell_type"].unique())]
        yield tissue, "AFCA", part.var_names.to_numpy(), [unit for unit in units if unit]


def worm_strata():
    from scripts.run_gao_worm_aging import normal_obs
    atlas = ad.read_h5ad(INTERIM / "worm-gao-integrated.h5ad", backed="r")
    obs = normal_obs(atlas)
    part = atlas[obs.index].to_memory()
    part.obs = obs
    keep = ~technical_gene_mask(part.var_names)
    part = part[:, keep]
    units = [aggregate_unit(part.X, part.obs, cell_type, "C. elegans", "Whole animal", "Gao_GSE229022",
                            "cell_type", "donor")
             for cell_type in sorted(part.obs["cell_type"].unique())]
    yield "Whole animal", "Gao_GSE229022", part.var_names.to_numpy(), [unit for unit in units if unit]


def zebrafish_strata():
    from scripts.run_zebrafish_testis_aging import standard_obs
    atlas = ad.read_h5ad(INTERIM / "zebrafish-testis-5-22mo.h5ad", backed="r")
    obs = standard_obs(atlas, "5-22mo")
    part = atlas.to_memory()
    part.obs = obs
    keep = ~technical_gene_mask(part.var_names)
    part = part[:, keep]
    units = [aggregate_unit(part.X, part.obs, cell_type, "Zebrafish", "Testis fine cell types",
                            "GSE275361", "fine_cell_type", "donor")
             for cell_type in sorted(part.obs["fine_cell_type"].unique()) if cell_type != "unknown"]
    yield "Testis fine cell types", "GSE275361", part.var_names.to_numpy(), [unit for unit in units if unit]


def nmr_strata():
    from scripts.run_cross_species_expression_vector import nmr_annotated_data
    annotated = nmr_annotated_data()
    raw = pd.read_csv(REPO / "data" / "raw" / "GSE152886_nmr_skin_raw.csv.gz", sep="\t", index_col=0)
    counts = raw[annotated.obs_names].T.to_numpy()
    obs = annotated.obs.copy()
    obs["sex_male"] = 0.0
    keep = ~technical_gene_mask(raw.index.astype(str))
    units = [aggregate_unit(counts[:, keep], obs, cell_type, "Naked mole-rat", "Skin", "GSE152886",
                            "cell_type", "sample")
             for cell_type in sorted(obs["cell_type"].unique())]
    yield "Skin", "GSE152886", raw.index.to_numpy()[keep], [unit for unit in units if unit]


def human_strata():
    import cellxgene_census
    from scripts.run_human_vector_gene_drivers import pseudobulk
    from scripts.run_m9_pipeline import COHORTS

    tests = pd.read_csv(RESULTS / "h3-cross-species-vector-tests.csv")
    eligible = tests[(tests["species"] == "Human") & (tests["metric"] == "profile_hellinger")
                     & (tests["model"] == "depth_adjusted")][["cohort", "tissue", "cell_type"]].drop_duplicates()
    census = cellxgene_census.open_soma(census_version="2025-11-08")
    for cohort, selected in eligible.groupby("cohort"):
        dataset_id = next(key for key, value in COHORTS.items() if value[0] == cohort)
        cell_filter = ",".join(f"'{value}'" for value in selected["cell_type"])
        obs = cellxgene_census.get_obs(
            census, "Homo sapiens",
            value_filter=(
                f"dataset_id == '{dataset_id}' and is_primary_data == True and disease == 'normal' "
                f"and raw_sum >= 1500 and cell_type in [{cell_filter}]"
            ),
            column_names=["donor_id", "development_stage"],
        )
        obs["age"] = obs["development_stage"].astype(str).str.extract(r"(\d+)-year-old").astype(float)
        donors = sorted(obs.loc[obs["age"].between(20, 120), "donor_id"].astype(str).unique())
        collected = {}
        genes = None
        for start in range(0, len(donors), 25):
            batch = donors[start:start + 25]
            donor_filter = ",".join(f"'{donor}'" for donor in batch)
            adata = cellxgene_census.get_anndata(
                census, "Homo sapiens",
                obs_value_filter=(
                    f"dataset_id == '{dataset_id}' and is_primary_data == True and disease == 'normal' "
                    f"and raw_sum >= 1500 and cell_type in [{cell_filter}] and donor_id in [{donor_filter}]"
                ),
                obs_column_names=["donor_id", "development_stage", "sex"],
                var_column_names=["feature_name"],
            )
            for cell_type in selected["cell_type"]:
                part = adata[adata.obs["cell_type"].astype(str) == cell_type]
                if part.n_obs == 0 or not (part.obs.groupby("donor_id", observed=True).size() >= 20).any():
                    continue
                expression, metadata, gene_names = pseudobulk(part)
                genes = gene_names
                if len(metadata):
                    collected.setdefault(cell_type, []).append((expression, metadata))
            print(f"Human {cohort} | donors {start + 1} to {min(start + 25, len(donors))} of {len(donors)}", flush=True)
        units = []
        for row in selected.itertuples():
            pieces = collected.get(row.cell_type, [])
            if not pieces:
                continue
            expression = np.vstack([piece[0] for piece in pieces])
            metadata = pd.concat([piece[1] for piece in pieces], ignore_index=True)
            if len(metadata) >= 6:
                units.append({"species": "Human", "tissue": row.tissue, "cohort": cohort,
                              "cell_type": row.cell_type, "expression": expression, "metadata": metadata})
        yield selected["tissue"].iloc[0], cohort, genes, units


def design(metadata, age=None):
    data = metadata.copy()
    data["age_model"] = data["age"] if age is None else age
    columns = ["age_model"]
    if data["sex_male"].nunique() > 1:
        columns.append("sex_male")
    if "mean_log_umi" in data and data["mean_log_umi"].nunique() > 1:
        columns.append("mean_log_umi")
    if "platform_nextseq" in data and data["platform_nextseq"].nunique() > 1:
        columns.append("platform_nextseq")
    return sm.add_constant(data[columns].astype(float), has_constant="add").to_numpy()


def permuted_ages(metadata, rng):
    age = metadata["age"].to_numpy().copy()
    strata = []
    for column in ("sex_male", "platform_nextseq"):
        if column in metadata and metadata[column].nunique() > 1:
            strata.append(column)
    if not strata:
        return rng.permutation(age)
    for _, indices in metadata.groupby(strata).groups.items():
        positions = metadata.index.get_indexer(indices)
        age[positions] = rng.permutation(age[positions])
    return age


def unit_slopes(unit, age=None):
    X = design(unit["metadata"], age=age)
    slopes, se, pvalues = vectorized_ols(unit["expression"], X, age_column=1)
    return slopes, se, pvalues


def correlations(slopes):
    matrix = np.vstack(slopes)
    pearson = pairwise_pearson(matrix)
    ranked = np.vstack([rankdata(row) for row in matrix])
    spearman = pairwise_pearson(ranked)
    return pearson, spearman


def analyze_stratum(species, tissue, cohort, genes, units, permutations, rng):
    units = [unit for unit in units if len(unit["metadata"]) >= (5 if species == "Naked mole-rat" else 6)]
    if len(units) < 2:
        return [], [], []
    fitted = [unit_slopes(unit) for unit in units]
    observed_pearson, observed_spearman = correlations([item[0] for item in fitted])
    pair_rows = []
    index = 0
    for first, second in itertools.combinations(range(len(units)), 2):
        pair_rows.append({
            "species": species, "tissue": tissue, "cohort": cohort,
            "cell_type_1": units[first]["cell_type"], "cell_type_2": units[second]["cell_type"],
            "pearson": observed_pearson[index], "spearman": observed_spearman[index],
        })
        index += 1
    permutation_rows = []
    assignments = [None]
    if species == "Naked mole-rat":
        sample_ages = sorted(units[0]["metadata"]["age"].unique())
        assignments = []
        for aged_indices in itertools.combinations(range(5), 2):
            age = np.array([sample_ages[1] if i in aged_indices else sample_ages[0] for i in range(5)])
            assignments.append(age)
        iterator = enumerate(assignments)
    else:
        iterator = ((iteration, None) for iteration in range(permutations + 1))
    for iteration, fixed_age in iterator:
        slopes = []
        for unit in units:
            if iteration == 0:
                slopes.append(unit_slopes(unit)[0])
            elif species == "Naked mole-rat":
                age_map = dict(zip(sorted(unit["metadata"]["donor"]), fixed_age))
                age = unit["metadata"]["donor"].map(age_map).to_numpy()
                slopes.append(unit_slopes(unit, age=age)[0])
            else:
                slopes.append(unit_slopes(unit, age=permuted_ages(unit["metadata"], rng))[0])
        pearson, spearman = correlations(slopes)
        permutation_rows.append({
            "species": species, "tissue": tissue, "cohort": cohort, "iteration": iteration,
            "n_units": len(units), "n_pairs": len(pearson),
            "median_pearson": np.median(pearson), "mean_pearson": np.mean(pearson),
            "median_spearman": np.median(spearman), "mean_spearman": np.mean(spearman),
        })
    unit_rows = []
    for unit, (slopes, se, pvalues) in zip(units, fitted):
        unit_rows.append({
            "species": species, "tissue": tissue, "cohort": cohort, "cell_type": unit["cell_type"],
            "n_donors": len(unit["metadata"]), "n_genes": len(genes),
            "mean_age_slope": np.mean(slopes), "median_age_slope": np.median(slopes),
            "fraction_positive": np.mean(slopes > 0),
        })
        if species == "Fruit fly":
            qvalues = multipletests(pvalues, method="fdr_bh")[1]
            significant = qvalues < 0.05
            fly_rows = pd.DataFrame({
                "species": species, "tissue": tissue, "cohort": cohort, "cell_type": unit["cell_type"],
                "gene": genes[significant], "slope": slopes[significant], "se": se[significant],
                "pvalue": pvalues[significant], "q_bh": qvalues[significant],
                "n_donors": len(unit["metadata"]),
            })
            fly_rows.to_csv(RESULTS / "reviewer-phase1-fly-gene-effects.csv", mode="a",
                            header=not (RESULTS / "reviewer-phase1-fly-gene-effects.csv").exists(), index=False)
    return pair_rows, permutation_rows, unit_rows


def summarize(permutations):
    rows = []
    for species, part in permutations.groupby("species"):
        aggregated = part.groupby("iteration").apply(lambda frame: pd.Series({
            metric: np.average(frame[metric], weights=frame["n_pairs"])
            for metric in ("median_pearson", "mean_pearson", "median_spearman", "mean_spearman")
        }), include_groups=False)
        observed = aggregated.loc[0]
        null = aggregated.drop(index=0)
        for metric in ("median_pearson", "mean_pearson", "median_spearman", "mean_spearman"):
            value = observed[metric]
            values = null[metric]
            stratum_values = part[part["iteration"] == 0][metric]
            rows.append({
                "species": species, "metric": metric, "observed": value,
                "pair_or_strata_ci95_low": stratum_values.quantile(0.025) if len(stratum_values) > 1 else np.nan,
                "pair_or_strata_ci95_high": stratum_values.quantile(0.975) if len(stratum_values) > 1 else np.nan,
                "null_mean": values.mean(), "null_ci95_low": values.quantile(0.025),
                "null_ci95_high": values.quantile(0.975),
                "permutation_p_greater": (1 + (values >= value).sum()) / (len(values) + 1),
                "permutation_p_two_sided": (1 + (np.abs(values - values.mean()) >= abs(value - values.mean())).sum()) / (len(values) + 1),
                "n_permutations": len(null),
            })
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--skip-human", action="store_true")
    parser.add_argument("--species", choices=["Mouse", "Human", "Fruit fly", "C. elegans", "Zebrafish", "Naked mole-rat"])
    parser.add_argument("--output-suffix", default="")
    args = parser.parse_args()
    output = RESULTS / f"reviewer-phase1-fly-gene-effects{args.output_suffix}.csv"
    if output.exists() and args.species in (None, "Fruit fly"):
        output.unlink()
    loaders = {"Mouse": mouse_strata, "Human": human_strata, "Fruit fly": fly_strata,
               "C. elegans": worm_strata, "Zebrafish": zebrafish_strata, "Naked mole-rat": nmr_strata}
    pairs = []
    permuted = []
    unit_rows = []
    for species, loader in loaders.items():
        if args.species and species != args.species:
            continue
        if species == "Human" and args.skip_human:
            continue
        rng = np.random.default_rng(20261105 + sum(map(ord, species)))
        for tissue, cohort, genes, units in loader():
            print(f"Analyzing {species} | {tissue} | {cohort}", flush=True)
            p, perm, u = analyze_stratum(species, tissue, cohort, genes, units, args.permutations, rng)
            pairs.extend(p)
            permuted.extend(perm)
            unit_rows.extend(u)
            pd.DataFrame(pairs).to_csv(RESULTS / f"reviewer-phase1-universality-pairs{args.output_suffix}.csv", index=False)
            pd.DataFrame(permuted).to_csv(RESULTS / f"reviewer-phase1-universality-permutations{args.output_suffix}.csv", index=False)
            pd.DataFrame(unit_rows).to_csv(RESULTS / f"reviewer-phase1-universality-units{args.output_suffix}.csv", index=False)
    result = pd.DataFrame(permuted)
    if not result.empty:
        summarize(result).to_csv(RESULTS / f"reviewer-phase1-universality-summary{args.output_suffix}.csv", index=False)


if __name__ == "__main__":
    main()
