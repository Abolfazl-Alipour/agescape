#!/usr/bin/env python3

from pathlib import Path
import sys

import anndata as ad
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

RESULTS = REPO / "results"
AGE_MAP = {"1m": 1.0, "3m": 3.0, "18m": 18.0, "21m": 21.0, "24m": 24.0, "30m": 30.0}


def marker_for_cell_type(cell_type):
    label = cell_type.lower()
    rules = [
        (("interneuron",), {"vip": "VIP", "pvalb": "PVALB", "sst": "SST"}),
        (("goblet",), "MUC5AC"),
        (("t cell", "thymocyte"), "CD3E"), (("b cell",), "CD79A"),
        (("nk cell", "natural killer"), "NKG7"), (("macrophage",), "C1QA"),
        (("monocyte",), "LY6C2"), (("endothelial", "endocardial", "fenestrated"), "PECAM1"),
        (("fibroblast", "stromal"), "COL1A1"), (("hepatocyte",), "ALB"),
        (("podocyte",), "NPHS2"), (("cardiomyocyte",), "MYH6"),
        (("neuron",), "SNAP25"), (("oligodendrocyte precursor",), "PDGFRA"),
        (("oligodendrocyte",), "MBP"), (("astrocyte",), "GFAP"),
        (("microglial",), "C1QA"), (("pericyte",), "RGS5"),
        (("pulmonary alveolar type 1",), "AGER"), (("pulmonary alveolar type 2", "type ii pneumocyte"), "SFTPC"),
        (("basal cell",), "KRT5"), (("satellite cell",), "PAX7"),
        (("smooth muscle",), "ACTA2"), (("ciliated",), "FOXJ1"),
        (("goblet",), "MUC5AC"), (("club cell",), "SCGB1A1"),
        (("pancreatic b cell",), "INS"), (("acinar",), "PRSS1"),
    ]
    for patterns, marker in rules:
        if any(pattern in label for pattern in patterns):
            if isinstance(marker, dict):
                return next((value for key, value in marker.items() if key in label), None)
            return marker
    return None


def marker_regression(records, species, cohort, cell_type, marker):
    if len(records) < 6 or records["marker_log_cp10k"].nunique() < 2:
        return None
    predictors = ["age", "sex_male", "mean_log_umi"]
    X = sm.add_constant(records[predictors].astype(float), has_constant="add")
    if np.linalg.matrix_rank(X.to_numpy()) < X.shape[1]:
        predictors.remove("sex_male")
        X = sm.add_constant(records[predictors].astype(float), has_constant="add")
    fit = sm.OLS(records["marker_log_cp10k"], X).fit()
    ci = fit.conf_int().loc["age"]
    return {"species": species, "cohort": cohort, "cell_type": cell_type, "marker": marker,
            "n_donors": len(records), "slope_log_cp10k_per_year": fit.params["age"],
            "se": fit.bse["age"], "ci95_low": ci.iloc[0], "ci95_high": ci.iloc[1],
            "pvalue": fit.pvalues["age"], "mean_marker_positive_fraction": records["marker_positive_fraction"].mean()}


def mouse_analysis():
    raw = ad.read_h5ad(REPO / "data" / "raw" / "tms_droplet_raw.h5ad", backed="r")
    clean = ad.read_h5ad(REPO / "data" / "interim" / "m1-qc-cleaned.h5ad", backed="r")
    clean_ids = set(clean.obs_names)
    excluded_rows = []
    for (tissue, age), part in raw.obs.groupby(["tissue", "age"], observed=True):
        unknown = part["cell_ontology_class"].astype(str).str.lower().isin(["unknown", "nan", "none", ""])
        excluded = ~part.index.isin(clean_ids)
        excluded_rows.append({"species": "Mouse", "cohort": "TMS", "tissue": tissue, "age": age,
                              "n_raw_cells": len(part), "n_excluded_qc": int(excluded.sum()),
                              "n_unassigned_raw": int(unknown.sum())})
    marker_rows = []
    for tissue in sorted(raw.obs["tissue"].unique()):
        raw_part = raw[raw.obs["tissue"] == tissue].to_memory()
        clean_part = raw_part[raw_part.obs_names.isin(clean_ids)].copy()
        clean_part.obs["age_numeric"] = clean_part.obs["age"].astype(str).map(AGE_MAP)
        clean_part.obs["sex_male"] = (clean_part.obs["sex"].astype(str).str.lower() == "male").astype(float)
        totals = np.asarray(clean_part.X.sum(axis=1)).ravel()
        gene_lookup = {gene.upper(): index for index, gene in enumerate(clean_part.var_names)}
        for cell_type in sorted(clean_part.obs["cell_ontology_class"].unique()):
            marker = marker_for_cell_type(str(cell_type))
            if marker is None or marker not in gene_lookup:
                continue
            ct_mask = (clean_part.obs["cell_ontology_class"] == cell_type).to_numpy()
            ct_obs = clean_part.obs[ct_mask]
            expression = np.asarray(clean_part.X[ct_mask, gene_lookup[marker]].toarray()).ravel()
            ct_totals = totals[ct_mask]
            records = []
            for donor in ct_obs["mouse.id"].unique():
                mask = (ct_obs["mouse.id"] == donor).to_numpy()
                if mask.sum() < 20:
                    continue
                records.append({"age": ct_obs.loc[mask, "age_numeric"].iloc[0],
                                "sex_male": ct_obs.loc[mask, "sex_male"].iloc[0],
                                "mean_log_umi": np.log(np.maximum(ct_totals[mask], 1)).mean(),
                                "marker_log_cp10k": np.log1p(expression[mask].sum() / ct_totals[mask].sum() * 1e4),
                                "marker_positive_fraction": np.mean(expression[mask] > 0)})
            row = marker_regression(pd.DataFrame(records), "Mouse", "TMS", str(cell_type), marker)
            if row:
                row["tissue"] = tissue
                marker_rows.append(row)
    return excluded_rows, marker_rows


def human_analysis():
    import cellxgene_census
    from scripts.run_m9_pipeline import COHORTS
    census = cellxgene_census.open_soma(census_version="2025-11-08")
    excluded_rows = []
    marker_rows = []
    for dataset_id, (cohort, tissue) in COHORTS.items():
        obs = cellxgene_census.get_obs(
            census, "Homo sapiens", value_filter=f"dataset_id == '{dataset_id}'",
            column_names=["donor_id", "development_stage", "disease", "cell_type", "is_primary_data", "raw_sum"],
        )
        obs["age"] = obs["development_stage"].astype(str).str.extract(r"(\d+)-year-old").astype(float)
        selected = obs[(obs["is_primary_data"] == True) & (obs["disease"] == "normal")
                       & obs["age"].between(20, 120)].copy()
        selected["unassigned"] = selected["cell_type"].astype(str).str.lower().isin(["unknown", "nan", "none", ""])
        selected["excluded"] = selected["raw_sum"] < 1500
        for age, part in selected.groupby("age"):
            excluded_rows.append({"species": "Human", "cohort": cohort, "tissue": tissue, "age": age,
                                  "n_raw_cells": len(part), "n_excluded_qc": int(part["excluded"].sum()),
                                  "n_unassigned_raw": int(part["unassigned"].sum())})
        cell_types = selected.loc[~selected["excluded"] & ~selected["unassigned"], "cell_type"].astype(str).unique()
        markers = {cell_type: marker_for_cell_type(cell_type) for cell_type in cell_types}
        markers = {cell_type: marker for cell_type, marker in markers.items() if marker}
        if not markers:
            continue
        marker_list = sorted(set(markers.values()))
        var_filter = ",".join(f"'{gene}'" for gene in marker_list)
        cell_filter = ",".join(f"'{cell_type}'" for cell_type in markers)
        adata = cellxgene_census.get_anndata(
            census, "Homo sapiens",
            obs_value_filter=(f"dataset_id == '{dataset_id}' and is_primary_data == True and disease == 'normal' "
                              f"and raw_sum >= 1500 and cell_type in [{cell_filter}]"),
            obs_column_names=["donor_id", "development_stage", "sex", "raw_sum"],
            var_value_filter=f"feature_name in [{var_filter}]",
            var_column_names=["feature_name"],
        )
        adata.obs["age"] = adata.obs["development_stage"].astype(str).str.extract(r"(\d+)-year-old").astype(float)
        adata = adata[adata.obs["age"].between(20, 120)].copy()
        gene_lookup = {gene: index for index, gene in enumerate(adata.var["feature_name"].astype(str))}
        for cell_type, marker in markers.items():
            if marker not in gene_lookup:
                continue
            mask = (adata.obs["cell_type"].astype(str) == cell_type).to_numpy()
            part = adata[mask]
            expression = np.asarray(part.X[:, gene_lookup[marker]].toarray()).ravel()
            records = []
            for donor in part.obs["donor_id"].unique():
                donor_mask = (part.obs["donor_id"] == donor).to_numpy()
                if donor_mask.sum() < 20:
                    continue
                library = part.obs.loc[donor_mask, "raw_sum"].to_numpy(dtype=float)
                records.append({"age": part.obs.loc[donor_mask, "age"].iloc[0],
                                "sex_male": float(part.obs.loc[donor_mask, "sex"].astype(str).str.lower().iloc[0] == "male"),
                                "mean_log_umi": np.log(np.maximum(library, 1)).mean(),
                                "marker_log_cp10k": np.log1p(expression[donor_mask].sum() / library.sum() * 1e4),
                                "marker_positive_fraction": np.mean(expression[donor_mask] > 0)})
            row = marker_regression(pd.DataFrame(records), "Human", cohort, cell_type, marker)
            if row:
                row["tissue"] = tissue
                marker_rows.append(row)
        print(f"Human annotation audit: {cohort}", flush=True)
    return excluded_rows, marker_rows


def main():
    mouse_excluded, mouse_markers = mouse_analysis()
    human_excluded, human_markers = human_analysis()
    excluded = pd.DataFrame(mouse_excluded + human_excluded)
    excluded["fraction_excluded_qc"] = excluded["n_excluded_qc"] / excluded["n_raw_cells"]
    excluded["fraction_unassigned_raw"] = excluded["n_unassigned_raw"] / excluded["n_raw_cells"]
    excluded.to_csv(RESULTS / "reviewer-phase2-annotation-exclusion-by-age.csv", index=False)
    markers = pd.DataFrame(mouse_markers + human_markers)
    if not markers.empty:
        markers["q_bh"] = np.nan
        valid = markers["pvalue"].notna()
        markers.loc[valid, "q_bh"] = multipletests(markers.loc[valid, "pvalue"], method="fdr_bh")[1]
    markers.to_csv(RESULTS / "reviewer-phase2-marker-stability.csv", index=False)


if __name__ == "__main__":
    main()
