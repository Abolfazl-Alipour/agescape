#!/usr/bin/env python3
"""Replicate MSSM neuronal depth-adjusted gene associations in HBCC."""

import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_human_vector_gene_drivers import gene_regressions, pseudobulk
from scripts.run_m9_pipeline import COHORTS

RESULTS = REPO / "results"
CELL_TYPES = [
    "GABAergic neuron", "VIP GABAergic cortical interneuron",
    "pvalb GABAergic cortical interneuron", "sst GABAergic cortical interneuron",
]


def main():
    import cellxgene_census

    dataset_id = next(key for key, value in COHORTS.items() if value[0] == "HBCC_DLPFC")
    census = cellxgene_census.open_soma(census_version="2025-11-08")
    rows = []
    for cell_type in CELL_TYPES:
        adata = cellxgene_census.get_anndata(
            census, "Homo sapiens",
            obs_value_filter=(
                f"dataset_id == '{dataset_id}' and is_primary_data == True and disease == 'normal' "
                f"and raw_sum >= 1500 and cell_type == '{cell_type}'"
            ),
            obs_column_names=["donor_id", "development_stage", "sex"],
            var_column_names=["feature_name"],
        )
        expression, metadata, gene_names = pseudobulk(adata)
        print(f"HBCC | {cell_type}: {len(metadata)} donors", flush=True)
        rows.extend(gene_regressions(expression, metadata, gene_names, "HBCC_DLPFC", cell_type))
    results = pd.DataFrame(rows)
    results.to_csv(RESULTS / "h3-human-hbcc-gene-replication-all.csv", index=False)
    results[results["model"] == "age_plus_sex_depth"].to_csv(
        RESULTS / "h3-human-hbcc-gene-replication-depth.csv", index=False
    )
    mssm = pd.read_csv(RESULTS / "h3-human-vector-gene-drivers-all.csv")
    mssm = mssm[mssm["model"] == "age_plus_sex_depth"].copy()
    hbcc = results[results["model"] == "age_plus_sex_depth"].copy()
    merged = mssm.merge(
        hbcc[["cell_type", "gene", "slope_log_cp10k_per_year", "pvalue", "q_bh"]],
        on=["cell_type", "gene"], suffixes=("_mssm", "_hbcc")
    )
    merged["same_direction"] = (
        merged["slope_log_cp10k_per_year_mssm"] * merged["slope_log_cp10k_per_year_hbcc"] > 0
    )
    merged["mssm_fdr"] = merged["q_bh_mssm"] < 0.05
    merged["hbcc_fdr"] = merged["q_bh_hbcc"] < 0.05
    merged.to_csv(RESULTS / "h3-human-neuronal-gene-replication-all.csv", index=False)
    merged[merged["mssm_fdr"]].to_csv(RESULTS / "h3-human-neuronal-gene-replication-significant.csv", index=False)
    summary = []
    for cell_type, part in merged.groupby("cell_type"):
        discovery = part[part["mssm_fdr"]]
        summary.append({
            "cell_type": cell_type,
            "mssm_fdr_genes": len(discovery),
            "same_direction_in_hbcc": int(discovery["same_direction"].sum()),
            "hbcc_nominal_same_direction": int((discovery["same_direction"] & (discovery["pvalue_hbcc"] < 0.05)).sum()),
            "both_fdr": int((discovery["same_direction"] & discovery["hbcc_fdr"]).sum()),
        })
    summary = pd.DataFrame(summary)
    summary.to_csv(RESULTS / "h3-human-neuronal-gene-replication-summary.csv", index=False)
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
