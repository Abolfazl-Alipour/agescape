#!/usr/bin/env python3
"""Test whole-expression-vector displacement independently in mouse, human, and NMR."""

import itertools
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_h3_expression_vector import hellinger_distance, raw_euclidean_distance
from scripts.run_m9_pipeline import COHORTS

RESULTS = REPO / "results"
MOUSE_AGES = {"1m": 1.0, "3m": 3.0, "18m": 18.0, "21m": 21.0, "24m": 24.0, "30m": 30.0}
NMR_AGES = {"1J": 1.1, "5J": 1.1, "6J": 1.1, "2V": 11.3, "3V": 11.3, "4V": 11.3}
NMR_GROUPS = {"1J": "young", "5J": "young", "6J": "young", "2V": "aged", "3V": "aged", "4V": "aged"}
AGE_RE = re.compile(r"(\d+)-year-old")


def technical_gene_mask(gene_names):
    genes = pd.Index(gene_names).astype(str).str.upper()
    return np.array([
        bool(re.match(r"^(MT-|RP[LS]|MRP[LS])", gene))
        or gene in {"XIST", "TSIX", "DDX3Y", "EIF2S3Y", "UTY", "KDM5D"}
        for gene in genes
    ])


def build_reference(matrix, obs, cell_column, donor_column):
    stats = {}
    for cell_type in obs[cell_column].unique():
        ct_mask = (obs[cell_column] == cell_type).to_numpy()
        ct_matrix = matrix[ct_mask]
        ct_obs = obs.loc[ct_mask]
        per_donor = {}
        for donor in ct_obs[donor_column].unique():
            donor_mask = (ct_obs[donor_column] == donor).to_numpy()
            donor_matrix = ct_matrix[donor_mask]
            per_donor[str(donor)] = {
                "counts": np.asarray(donor_matrix.sum(axis=0)).ravel(),
                "n_cells": int(donor_matrix.shape[0]),
            }
        stats[str(cell_type)] = {
            "counts": np.asarray(ct_matrix.sum(axis=0)).ravel(),
            "n_cells": int(ct_matrix.shape[0]),
            "donors": per_donor,
        }
    return stats


def reference_without_donor(stats, donor):
    result = {"counts": stats["counts"].copy(), "n_cells": stats["n_cells"]}
    if donor in stats["donors"]:
        result["counts"] -= stats["donors"][donor]["counts"]
        result["n_cells"] -= stats["donors"][donor]["n_cells"]
    return result


def donor_distances(matrix, obs, reference, species, tissue, cohort, cell_column, donor_column,
                    age_column, sex_column, eligible_units=None):
    rows = []
    for cell_type in obs[cell_column].unique():
        cell_type = str(cell_type)
        if eligible_units is not None and cell_type not in eligible_units:
            continue
        if cell_type not in reference:
            continue
        ct_mask = (obs[cell_column].astype(str) == cell_type).to_numpy()
        ct_obs = obs.loc[ct_mask]
        ct_matrix = matrix[ct_mask]
        for donor in ct_obs[donor_column].unique():
            donor_mask = (ct_obs[donor_column] == donor).to_numpy()
            if donor_mask.sum() < 20:
                continue
            donor_obs = ct_obs.loc[donor_mask]
            donor_matrix = ct_matrix[donor_mask]
            ref = reference_without_donor(reference[cell_type], str(donor))
            if ref["n_cells"] < 50 or ref["counts"].sum() <= 0:
                continue
            totals = np.asarray(donor_matrix.sum(axis=1)).ravel()
            rows.append({
                "species": species,
                "tissue": tissue,
                "cohort": cohort,
                "cell_type": cell_type,
                "donor": str(donor),
                "age": float(donor_obs[age_column].iloc[0]),
                "sex_male": float(donor_obs[sex_column].iloc[0]),
                "n_cells": int(donor_matrix.shape[0]),
                "mean_log_umi": float(np.log(totals).mean()),
                "profile_hellinger": float(hellinger_distance(
                    donor_matrix, ref["counts"]
                ).mean()),
                "raw_count_euclidean": float(raw_euclidean_distance(
                    donor_matrix, ref["counts"] / ref["n_cells"]
                ).mean()),
            })
    return rows


def fit_unit(records, metric, species, tissue, cohort, cell_type, model):
    records = records.dropna(subset=[metric]).copy()
    if len(records) < 6:
        return None
    predictors = ["age", "sex_male"]
    if model == "depth_adjusted":
        predictors.append("mean_log_umi")
    X = sm.add_constant(records[predictors].astype(float), has_constant="add")
    if np.linalg.matrix_rank(X.to_numpy()) < X.shape[1]:
        predictors.remove("sex_male")
        X = sm.add_constant(records[predictors].astype(float), has_constant="add")
    if np.linalg.matrix_rank(X.to_numpy()) < X.shape[1]:
        return None
    fit = sm.OLS(records[metric].astype(float), X).fit()
    ci = fit.conf_int().loc["age"]
    return {
        "species": species, "tissue": tissue, "cohort": cohort, "cell_type": cell_type,
        "metric": metric, "model": model, "n_donors": len(records),
        "n_cells": int(records["n_cells"].sum()), "slope_per_year": fit.params["age"],
        "se": fit.bse["age"], "ci95_low": ci.iloc[0], "ci95_high": ci.iloc[1],
        "pvalue": fit.pvalues["age"], "r2": fit.rsquared,
    }


def fit_global(records, metric, species, tissue, cohort, model, age_scale, donor_set="all_donor_ids"):
    data = records.dropna(subset=[metric]).copy()
    if data.empty:
        return None
    data["unit"] = data["tissue"] + " | " + data["cohort"] + " | " + data["cell_type"]
    standardized = []
    for _, part in data.groupby("unit"):
        part = part.copy()
        if len(part) < 6:
            continue
        scale = part[metric].std(ddof=1)
        if not np.isfinite(scale) or scale <= 0:
            continue
        part["outcome"] = (part[metric] - part[metric].mean()) / scale
        part["age_scaled"] = part["age"] / age_scale
        part["sex_centered"] = part["sex_male"] - part["sex_male"].mean()
        part["depth_centered"] = part["mean_log_umi"] - part["mean_log_umi"].mean()
        standardized.append(part)
    if not standardized:
        return None
    data = pd.concat(standardized, ignore_index=True)
    data["donor_cluster"] = data["cohort"] + " | " + data["donor"]
    if data["donor_cluster"].nunique() < 6:
        return None
    unit_dummies = pd.get_dummies(data["unit"], drop_first=True, dtype=float)
    predictors = ["age_scaled", "sex_centered"]
    if model == "depth_adjusted":
        predictors.append("depth_centered")
    predictors = [name for name in predictors if np.ptp(data[name]) > 1e-12]
    X = pd.concat([data[predictors], unit_dummies], axis=1)
    X = sm.add_constant(X, has_constant="add")
    if np.linalg.matrix_rank(X.to_numpy()) < X.shape[1]:
        return None
    fit = sm.OLS(data["outcome"], X).fit(
        cov_type="cluster", cov_kwds={"groups": data["donor_cluster"]}, use_t=True
    )
    ci = fit.conf_int().loc["age_scaled"]
    return {
        "species": species, "tissue": tissue, "cohort": cohort, "metric": metric, "model": model,
        "donor_set": donor_set, "age_scale_years": age_scale,
        "n_donors": data["donor_cluster"].nunique(),
        "n_units": data["unit"].nunique(), "n_donor_units": len(data),
        "effect_sd": fit.params["age_scaled"], "se": fit.bse["age_scaled"],
        "ci95_low": ci.iloc[0], "ci95_high": ci.iloc[1], "pvalue": fit.pvalues["age_scaled"],
    }


def mouse_analysis():
    adata = ad.read_h5ad(REPO / "data" / "interim" / "m1-qc-cleaned.h5ad", backed="r")
    donor_rows = []
    coverage = []
    for tissue in sorted(adata.obs["tissue"].unique()):
        mask = (
            (adata.obs["tissue"] == tissue)
            & ~adata.obs["predicted_doublet"]
            & (adata.obs["n_counts"] >= 1500)
            & (adata.obs["cell_ontology_class"].astype(str) != "unknown")
        )
        part = adata[mask].to_memory()
        part.obs["age"] = part.obs["age"].astype(str).map(MOUSE_AGES) / 12.0
        part.obs["sex_male"] = (part.obs["sex"].astype(str).str.lower() == "male").astype(float)
        part.obs["mouse.id"] = part.obs["mouse.id"].astype(str)
        keep = ~technical_gene_mask(part.var_names)
        part = part[:, keep].copy()
        young = part[part.obs["age"] <= 0.25].copy()
        reference = build_reference(young.X, young.obs, "cell_ontology_class", "mouse.id")
        before = len(donor_rows)
        donor_rows.extend(donor_distances(
            part.X, part.obs, reference, "Mouse", tissue, "TMS", "cell_ontology_class",
            "mouse.id", "age", "sex_male"
        ))
        coverage.append({"species": "Mouse", "tissue": tissue, "cohort": "TMS",
                         "n_donor_units": len(donor_rows) - before})
        print(f"Mouse {tissue}: {len(donor_rows) - before} donor-cell-type rows", flush=True)
    return donor_rows, coverage


def human_eligible_units(obs):
    counts = (obs.groupby(["cell_type", "donor_id"], observed=True).size().rename("n_cells").reset_index())
    donor_ages = obs.groupby("donor_id", observed=True)["age"].first()
    counts["age"] = counts["donor_id"].map(donor_ages)
    eligible = []
    for cell_type, part in counts.groupby("cell_type"):
        part = part[part["n_cells"] >= 20]
        young = part[part["age"].between(20, 40)]
        if part["donor_id"].nunique() >= 6 and young["donor_id"].nunique() >= 2:
            eligible.append(str(cell_type))
    return set(eligible)


def human_analysis():
    import cellxgene_census

    census = cellxgene_census.open_soma(census_version="2025-11-08")
    donor_rows = []
    coverage = []
    for dataset_id, (cohort, tissue) in COHORTS.items():
        obs = cellxgene_census.get_obs(
            census, "Homo sapiens",
            value_filter=f"dataset_id == '{dataset_id}'",
            column_names=["dataset_id", "donor_id", "development_stage", "disease", "cell_type",
                          "raw_sum", "sex", "is_primary_data"],
        )
        obs["donor_id"] = obs["donor_id"].astype(str)
        obs["age"] = obs["development_stage"].astype(str).str.extract(AGE_RE).astype(float)
        obs = obs[
            (obs["is_primary_data"] == True)
            & (obs["disease"] == "normal")
            & obs["age"].between(20, 120)
            & (obs["raw_sum"] >= 1500)
            & (obs["cell_type"].astype(str) != "unknown")
        ].copy()
        if obs.empty:
            coverage.append({"species": "Human", "tissue": tissue, "cohort": cohort, "n_donor_units": 0})
            continue
        obs["sex_male"] = (obs["sex"].astype(str).str.lower() == "male").astype(float)
        eligible_units = human_eligible_units(obs)
        obs = obs[obs["cell_type"].astype(str).isin(eligible_units)].copy()
        young_donors = sorted(obs.loc[obs["age"] <= 40, "donor_id"].unique())
        donors = sorted(obs["donor_id"].unique())
        if not eligible_units or len(young_donors) < 2:
            coverage.append({"species": "Human", "tissue": tissue, "cohort": cohort, "n_donor_units": 0})
            print(f"Human {cohort}: no eligible units", flush=True)
            continue
        donor_filter = ",".join(f"'{donor}'" for donor in young_donors)
        reference_data = cellxgene_census.get_anndata(
            census, "Homo sapiens",
            obs_value_filter=(
                f"dataset_id == '{dataset_id}' and is_primary_data == True and disease == 'normal' "
                f"and raw_sum >= 1500 and donor_id in [{donor_filter}]"
            ),
            obs_column_names=["donor_id", "cell_type", "development_stage", "sex"],
            var_column_names=["feature_name"],
        )
        reference_data.obs["donor_id"] = reference_data.obs["donor_id"].astype(str)
        keep = ~technical_gene_mask(reference_data.var["feature_name"])
        reference_data = reference_data[:, keep].copy()
        reference = build_reference(reference_data.X, reference_data.obs, "cell_type", "donor_id")
        before = len(donor_rows)
        batch_size = 50
        for start in range(0, len(donors), batch_size):
            batch = donors[start:start + batch_size]
            donor_filter = ",".join(f"'{donor}'" for donor in batch)
            part = cellxgene_census.get_anndata(
                census, "Homo sapiens",
                obs_value_filter=(
                    f"dataset_id == '{dataset_id}' and is_primary_data == True and disease == 'normal' "
                    f"and raw_sum >= 1500 and donor_id in [{donor_filter}]"
                ),
                obs_column_names=["donor_id", "cell_type", "development_stage", "sex"],
                var_column_names=["feature_name"],
            )
            part.obs["donor_id"] = part.obs["donor_id"].astype(str)
            part.obs["age"] = part.obs["development_stage"].astype(str).str.extract(AGE_RE).astype(float)
            part.obs["sex_male"] = (part.obs["sex"].astype(str).str.lower() == "male").astype(float)
            part = part[:, keep].copy()
            donor_rows.extend(donor_distances(
                part.X, part.obs, reference, "Human", tissue, cohort, "cell_type", "donor_id",
                "age", "sex_male", eligible_units=eligible_units
            ))
            print(f"Human {cohort}: donors {start + 1} to {min(start + batch_size, len(donors))} of {len(donors)}", flush=True)
        coverage.append({"species": "Human", "tissue": tissue, "cohort": cohort,
                         "n_donor_units": len(donor_rows) - before})
    return donor_rows, coverage


def nmr_annotated_data():
    raw = pd.read_csv(REPO / "data" / "raw" / "GSE152886_nmr_skin_raw.csv.gz", sep="\t", index_col=0)
    obs = pd.DataFrame({"barcode": raw.columns})
    obs["sample"] = obs["barcode"].str.split("_").str[0]
    obs = obs[obs["sample"] != "2V"].copy()
    obs["age"] = obs["sample"].map(NMR_AGES)
    obs["group"] = obs["sample"].map(NMR_GROUPS)
    X = sp.csr_matrix(raw[obs["barcode"]].T.values)
    adata = sc.AnnData(X=X, obs=obs.set_index("barcode"))
    adata.var_names = raw.index.astype(str)
    adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
    adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    adata = adata[(adata.obs["n_counts"] >= 2000) & (adata.obs["n_genes"] >= 700)].copy()
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
    adata = adata[:, retain].copy()
    sc.pp.scale(adata, max_value=10)
    sc.pp.pca(adata, n_comps=50, random_state=42)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=42)
    sc.tl.leiden(adata, resolution=0.4, random_state=42)
    annotations = []
    for cluster in sorted(adata.obs["leiden"].unique()):
        cells = adata[adata.obs["leiden"] == cluster]
        scores = {}
        for name, genes in markers.items():
            present = [gene for gene in genes if gene in cells.var_names]
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
    adata.obs["cell_type"] = adata.obs["leiden"].astype(str).map(annotation_map)
    adata.X = raw_counts[:, retain.to_numpy()]
    symbols = pd.read_csv(
        REPO / "data" / "interim" / "nmr-ensg2symbol.tsv", sep="\t", skiprows=1,
        names=["ensembl_id", "symbol"],
    ).set_index("ensembl_id")["symbol"]
    gene_symbols = adata.var_names.map(symbols).fillna("").astype(str)
    keep = ~technical_gene_mask(gene_symbols)
    return adata[:, keep].copy()


def nmr_analysis():
    adata = nmr_annotated_data()
    adata.obs["sex_male"] = 0.0
    young = adata[adata.obs["group"] == "young"].copy()
    reference = build_reference(young.X, young.obs, "cell_type", "sample")
    rows = donor_distances(
        adata.X, adata.obs, reference, "Naked mole-rat", "Skin", "GSE152886", "cell_type",
        "sample", "age", "sex_male"
    )
    print(f"NMR skin: {len(rows)} donor-cell-type rows", flush=True)
    return rows, [{"species": "Naked mole-rat", "tissue": "Skin", "cohort": "GSE152886",
                   "n_donor_units": len(rows)}]


def nmr_tests(donors):
    rows = []
    global_rows = []
    scales = {}
    for cell_type, part in donors.groupby("cell_type"):
        if len(part) != 5:
            continue
        for metric in ("profile_hellinger", "raw_count_euclidean"):
            part = part.copy()
            aged = part["group"] == "aged"
            observed = part.loc[aged, metric].mean() - part.loc[~aged, metric].mean()
            permutations = []
            for indices in itertools.combinations(range(len(part)), 2):
                mask = np.zeros(len(part), dtype=bool)
                mask[list(indices)] = True
                permutations.append(part.loc[mask, metric].mean() - part.loc[~mask, metric].mean())
            pvalue = float(np.mean(np.abs(permutations) >= abs(observed) - 1e-12))
            pooled_sd = np.sqrt((part.loc[aged, metric].var(ddof=1) + part.loc[~aged, metric].var(ddof=1)) / 2)
            scales[(cell_type, metric)] = pooled_sd
            rows.append({"species": "Naked mole-rat", "tissue": "Skin", "cohort": "GSE152886",
                         "cell_type": cell_type, "metric": metric, "model": "exact_two_aged_of_five",
                         "n_donors": len(part), "n_cells": int(part["n_cells"].sum()),
                         "aged_minus_young": observed, "standardized_effect": observed / pooled_sd,
                         "pvalue": pvalue})
    tests = pd.DataFrame(rows)
    tests["q_bh"] = np.nan
    for metric, indices in tests.groupby("metric").groups.items():
        tests.loc[indices, "q_bh"] = multipletests(tests.loc[indices, "pvalue"], method="fdr_bh")[1]
    sample_ids = sorted(donors["donor"].unique())
    for metric, part in tests.groupby("metric"):
        observed = part["standardized_effect"].mean()
        permuted = []
        for indices in itertools.combinations(range(5), 2):
            aged_samples = {sample_ids[index] for index in indices}
            effects = []
            for cell_type in part["cell_type"]:
                ct = donors[donors["cell_type"] == cell_type]
                aged = ct["donor"].isin(aged_samples)
                effect = (ct.loc[aged, metric].mean() - ct.loc[~aged, metric].mean()) / scales[(cell_type, metric)]
                effects.append(effect)
            permuted.append(np.mean(effects))
        pvalue = float(np.mean(np.abs(permuted) >= abs(observed) - 1e-12))
        global_rows.append({"species": "Naked mole-rat", "tissue": "Skin", "cohort": "GSE152886",
                            "metric": metric, "model": "exact_two_aged_of_five", "n_donors": 5,
                            "n_units": len(part), "age_scale_years": 3.1,
                            "effect_sd": observed * 3.1 / 10.2, "pvalue": pvalue})
    return tests, pd.DataFrame(global_rows)


def summarize(donors):
    continuous = donors[donors["species"] != "Naked mole-rat"].copy()
    result_rows = []
    global_rows = []
    for (species, tissue, cohort, cell_type), part in continuous.groupby(
        ["species", "tissue", "cohort", "cell_type"]
    ):
        for metric in ("profile_hellinger", "raw_count_euclidean"):
            for model in ("unadjusted", "depth_adjusted"):
                row = fit_unit(part, metric, species, tissue, cohort, cell_type, model)
                if row:
                    result_rows.append(row)
    for (species, tissue, cohort), part in continuous.groupby(["species", "tissue", "cohort"]):
        age_scale = 0.4 if species == "Mouse" else 12.25
        donor_sets = [("all_donor_ids", part)]
        if species == "Mouse":
            pooled = set(pd.read_csv(RESULTS / "m1-animal-tissue-design.csv")
                         .query("is_pooled_id")["mouse.id"].astype(str))
            donor_sets.append(("exclude_pooled_ids", part[~part["donor"].isin(pooled)]))
        for donor_set, selected in donor_sets:
            for metric in ("profile_hellinger", "raw_count_euclidean"):
                for model in ("unadjusted", "depth_adjusted"):
                    row = fit_global(selected, metric, species, tissue, cohort, model, age_scale, donor_set)
                    if row:
                        global_rows.append(row)
    for species, part in continuous.groupby("species"):
        age_scale = 0.4 if species == "Mouse" else 12.25
        donor_sets = [("all_donor_ids", part)]
        if species == "Mouse":
            pooled = set(pd.read_csv(RESULTS / "m1-animal-tissue-design.csv")
                         .query("is_pooled_id")["mouse.id"].astype(str))
            donor_sets.append(("exclude_pooled_ids", part[~part["donor"].isin(pooled)]))
        for donor_set, selected in donor_sets:
            for metric in ("profile_hellinger", "raw_count_euclidean"):
                for model in ("unadjusted", "depth_adjusted"):
                    row = fit_global(selected, metric, species, "All tested tissues", "All cohorts",
                                     model, age_scale, donor_set)
                    if row:
                        global_rows.append(row)
    results = pd.DataFrame(result_rows)
    results["q_bh"] = np.nan
    for (_, _, _, _), indices in results.groupby(["cohort", "metric", "model", "species"]).groups.items():
        results.loc[indices, "q_bh"] = multipletests(results.loc[indices, "pvalue"], method="fdr_bh")[1]
    nmr_result, nmr_global = nmr_tests(donors[donors["species"] == "Naked mole-rat"])
    nmr_global["donor_set"] = "all_donor_ids"
    results = pd.concat([results, nmr_result], ignore_index=True)
    globals_out = pd.concat([pd.DataFrame(global_rows), nmr_global], ignore_index=True)
    results.to_csv(RESULTS / "h3-cross-species-vector-tests.csv", index=False)
    globals_out.to_csv(RESULTS / "h3-cross-species-vector-global.csv", index=False)
    print(globals_out.to_string(index=False), flush=True)


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--summarize-only", action="store_true")
    args = parser.parse_args()
    if args.summarize_only:
        summarize(pd.read_csv(RESULTS / "h3-cross-species-vector-donors.csv"))
        return
    all_donors = []
    coverage = []
    for analysis in (mouse_analysis, human_analysis, nmr_analysis):
        rows, covered = analysis()
        all_donors.extend(rows)
        coverage.extend(covered)
    donors = pd.DataFrame(all_donors)
    donors["group"] = np.where(donors["species"] == "Naked mole-rat",
                               np.where(donors["age"] > 5, "aged", "young"), "continuous")
    donors.to_csv(RESULTS / "h3-cross-species-vector-donors.csv", index=False)
    pd.DataFrame(coverage).to_csv(RESULTS / "h3-cross-species-vector-coverage.csv", index=False)
    summarize(donors)


if __name__ == "__main__":
    main()
