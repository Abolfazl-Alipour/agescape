#!/usr/bin/env python3
"""Test age-related changes in cell-type fractions within organs in three species."""

import itertools
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
RESULTS = REPO / "results"
AGE_MAP = {"1m": 1.0, "3m": 3.0, "18m": 18.0, "21m": 21.0, "24m": 24.0, "30m": 30.0}
NMR_AGES = {"1J": 1.1, "5J": 1.1, "6J": 1.1, "2V": 11.3, "3V": 11.3, "4V": 11.3}
NMR_GROUPS = {"1J": "young", "5J": "young", "6J": "young", "2V": "aged", "3V": "aged", "4V": "aged"}


def clr_table(counts, pseudocount=0.5):
    adjusted = counts.astype(float) + pseudocount
    proportions = adjusted.div(adjusted.sum(axis=1), axis=0)
    log_proportions = np.log(proportions)
    return log_proportions.sub(log_proportions.mean(axis=1), axis=0)


def count_table(obs, donor_column, cell_column, eligible_donors=None):
    donors = obs[donor_column].astype(str).unique() if eligible_donors is None else list(eligible_donors)
    counts = (obs.groupby([donor_column, cell_column], observed=True).size()
              .unstack(fill_value=0).reindex(donors, fill_value=0))
    valid = (counts > 0).sum(axis=0)
    return counts.loc[:, valid >= 6]


def regress_fractions(counts, metadata, species, organ, cohort, sensitivity, observed_only=False):
    clr = clr_table(counts)
    rows = []
    for cell_type in counts.columns:
        if observed_only:
            donors = counts.index[counts[cell_type] > 0]
        else:
            donors = counts.index
        if len(donors) < 6:
            continue
        X = sm.add_constant(metadata.loc[donors, ["age", "sex"]])
        model = "age_plus_sex"
        if np.linalg.matrix_rank(X.to_numpy()) < X.shape[1]:
            X = sm.add_constant(metadata.loc[donors, ["age"]])
            model = "age_only_rank_deficient_sex"
        fit = sm.OLS(clr.loc[donors, cell_type], X).fit()
        rows.append({
            "species": species, "organ": organ, "cohort": cohort,
            "cell_type": cell_type, "sensitivity": sensitivity, "model": model,
            "effect_type": "CLR_slope_per_year", "effect": fit.params["age"],
            "se": fit.bse["age"], "pvalue": fit.pvalues["age"],
            "n_donors": len(donors), "n_present": int((counts[cell_type] > 0).sum()),
            "total_cells": int(counts[cell_type].sum()),
        })
    if not rows:
        return []
    qvalues = multipletests([row["pvalue"] for row in rows], method="fdr_bh")[1]
    for row, qvalue in zip(rows, qvalues):
        row["q_bh"] = qvalue
    return rows


def mouse_tests():
    adata = ad.read_h5ad(REPO / "data" / "interim" / "m1-qc-cleaned.h5ad", backed="r")
    obs = adata.obs[(~adata.obs["predicted_doublet"])].copy()
    obs["age"] = obs["age"].astype(str).map(AGE_MAP)
    obs["sex"] = (obs["sex"].astype(str).str.lower() == "male").astype(float)
    pooled = set(pd.read_csv(RESULTS / "m1-animal-tissue-design.csv")
                 .query("is_pooled_id")["mouse.id"].astype(str))
    rows = []
    coverage = []
    for tissue in sorted(obs["tissue"].unique()):
        part = obs[(obs["tissue"] == tissue) & (obs["cell_ontology_class"].astype(str) != "unknown")].copy()
        coverage.append({"species": "Mouse", "organ": tissue, "cohort": "TMS",
                         "n_donors": part["mouse.id"].nunique()})
        counts = count_table(part, "mouse.id", "cell_ontology_class")
        metadata = part.groupby("mouse.id").agg(age=("age", "first"), sex=("sex", "first"))
        rows.extend(regress_fractions(counts, metadata, "Mouse", tissue, "TMS", "all_donor_ids"))
        unpooled = counts.index[~counts.index.isin(pooled)]
        rows.extend(regress_fractions(counts.loc[unpooled], metadata, "Mouse", tissue, "TMS", "exclude_pooled_ids"))
    return rows, coverage


def human_tests():
    source = pd.read_csv(RESULTS / "m9-donor-cell-counts-prethreshold.csv")
    rows = []
    coverage = []
    for cohort in sorted(source["cohort"].unique()):
        part = source[(source["cohort"] == cohort) & (source["age"] >= 20)].copy()
        organ = part["group"].iloc[0]
        coverage.append({"species": "Human", "organ": organ, "cohort": cohort,
                         "n_donors": part["donor_id"].nunique()})
        counts = (part.set_index(["donor_id", "cell_type"])["n_cells"]
                  .unstack(fill_value=0).sort_index())
        counts = counts.loc[:, (counts > 0).sum(axis=0) >= 6]
        metadata = part.groupby("donor_id").agg(age=("age", "first"), sex=("sex", "first"))
        metadata["sex"] = (metadata["sex"].astype(str).str.lower() == "male").astype(float)
        rows.extend(regress_fractions(counts, metadata, "Human", organ, cohort, "all_adult_donors"))
        rows.extend(regress_fractions(counts, metadata, "Human", organ, cohort, "observed_records_only", observed_only=True))
    return rows, coverage


def nmr_annotations(resolution):
    raw = pd.read_csv(REPO / "data" / "raw" / "GSE152886_nmr_skin_raw.csv.gz", sep="\t", index_col=0)
    obs = pd.DataFrame({"barcode": raw.columns})
    obs["sample"] = obs["barcode"].str.split("_").str[0]
    obs = obs[obs["sample"] != "2V"].copy()
    obs["age_years"] = obs["sample"].map(NMR_AGES)
    obs["group"] = obs["sample"].map(NMR_GROUPS)
    X = sp.csr_matrix(raw[obs["barcode"]].T.values)
    adata = sc.AnnData(X=X, obs=obs.set_index("barcode"))
    adata.var_names = raw.index.astype(str)
    adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
    adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    adata = adata[(adata.obs["n_counts"] >= 2000) & (adata.obs["n_genes"] >= 700)].copy()
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
    sc.tl.leiden(adata, resolution=resolution, random_state=42)

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
        annotations.append({
            "resolution": resolution, "cluster": cluster, "compartment": compartment,
            "n_cells": cells.n_obs, "immune_marker_mean": scores["immune"],
            "basal_marker_mean": scores["basal"], "corneum_marker_mean": scores["corneum"],
        })
    annotation_df = pd.DataFrame(annotations)
    adata.obs["compartment"] = adata.obs["leiden"].astype(str).map(
        annotation_df.set_index("cluster")["compartment"]
    )
    counts = (adata.obs.groupby(["sample", "compartment"], observed=True).size()
              .unstack(fill_value=0).reindex(columns=["basal", "intermediate", "corneum", "immune"], fill_value=0))
    return counts, annotation_df


def nmr_tests():
    all_rows = []
    annotation_tables = []
    for resolution in (0.4, 0.6):
        counts, annotations = nmr_annotations(resolution)
        annotation_tables.append(annotations)
        fractions = counts.div(counts.sum(axis=1), axis=0)
        fractions.to_csv(RESULTS / f"nmr-skin-cell-fractions-resolution-{resolution:.1f}.csv")
        aged = fractions.index[fractions.index.map(NMR_GROUPS) == "aged"]
        young = fractions.index[fractions.index.map(NMR_GROUPS) == "young"]
        for compartment in fractions.columns:
            observed = fractions.loc[aged, compartment].mean() - fractions.loc[young, compartment].mean()
            permutations = []
            for indices in itertools.combinations(range(5), 2):
                aged_mask = np.zeros(fractions.shape[0], dtype=bool)
                aged_mask[list(indices)] = True
                permutations.append(
                    fractions.loc[aged_mask, compartment].mean() - fractions.loc[~aged_mask, compartment].mean()
                )
            pvalue = float(np.mean(np.abs(permutations) >= abs(observed) - 1e-12))
            all_rows.append({
                "species": "Naked mole-rat", "organ": "Skin", "cohort": "GSE152886",
                "cell_type": compartment, "sensitivity": f"annotation_resolution_{resolution:.1f}",
                "effect_type": "aged_minus_young_fraction", "effect": observed,
                "se": np.nan, "pvalue": pvalue, "n_donors": 5,
                "n_present": int((counts[compartment] > 0).sum()),
                "total_cells": int(counts[compartment].sum()),
            })
        qvalues = multipletests([row["pvalue"] for row in all_rows[-4:]], method="fdr_bh")[1]
        for row, qvalue in zip(all_rows[-4:], qvalues):
            row["q_bh"] = qvalue
    pd.concat(annotation_tables, ignore_index=True).to_csv(RESULTS / "nmr-skin-cluster-annotation.csv", index=False)
    return all_rows


def main():
    RESULTS.mkdir(exist_ok=True)
    mouse_rows, mouse_coverage = mouse_tests()
    human_rows, human_coverage = human_tests()
    nmr_coverage = [{"species": "Naked mole-rat", "organ": "Skin", "cohort": "GSE152886", "n_donors": 5}]
    results = pd.DataFrame(mouse_rows + human_rows + nmr_tests())
    results.to_csv(RESULTS / "organ-composition-age-tests.csv", index=False)
    primary = results[results["sensitivity"].isin(["all_donor_ids", "all_adult_donors", "annotation_resolution_0.4"])]
    statistics = (primary.groupby(["species", "organ", "cohort"], observed=True)
                  .agg(n_cell_types=("cell_type", "nunique"),
                       n_significant=("q_bh", lambda x: int((x < 0.05).sum())), min_q=("q_bh", "min"))
                  .reset_index())
    summary = (pd.DataFrame(mouse_coverage + human_coverage + nmr_coverage)
               .merge(statistics, on=["species", "organ", "cohort"], how="left"))
    summary[["n_cell_types", "n_significant"]] = summary[["n_cell_types", "n_significant"]].fillna(0).astype(int)
    summary["answer"] = "insufficient donor coverage"
    tested = summary["n_cell_types"] > 0
    summary.loc[tested, "answer"] = np.where(summary.loc[tested, "n_significant"] > 0,
                                             "age-associated fraction change detected",
                                             "no fraction change detected at FDR 0.05")
    summary.to_csv(RESULTS / "organ-composition-summary.csv", index=False)
    print(summary.to_string(index=False))
    print("\nNMR tests:")
    print(results[results["species"] == "Naked mole-rat"].to_string(index=False))


if __name__ == "__main__":
    main()
