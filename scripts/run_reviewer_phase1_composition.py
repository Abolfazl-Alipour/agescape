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

from scripts.extended_atlas_utils import clr_table, composition_counts
from scripts.run_organ_composition import AGE_MAP, count_table

RESULTS = REPO / "results"


def core_inputs():
    adata = ad.read_h5ad(REPO / "data" / "interim" / "m1-qc-cleaned.h5ad", backed="r")
    obs = adata.obs[~adata.obs["predicted_doublet"]].copy()
    obs["age_numeric"] = obs["age"].astype(str).map(AGE_MAP)
    obs["sex_male"] = (obs["sex"].astype(str).str.lower() == "male").astype(float)
    for tissue in sorted(obs["tissue"].unique()):
        part = obs[(obs["tissue"] == tissue) & (obs["cell_ontology_class"].astype(str) != "unknown")]
        counts = count_table(part, "mouse.id", "cell_ontology_class")
        metadata = part.groupby("mouse.id").agg(age=("age_numeric", "first"), sex_male=("sex_male", "first"))
        yield "Mouse", tissue, "TMS", counts, metadata
    source = pd.read_csv(RESULTS / "m9-donor-cell-counts-prethreshold.csv")
    for cohort in sorted(source["cohort"].unique()):
        part = source[(source["cohort"] == cohort) & (source["age"] >= 20)].copy()
        counts = part.set_index(["donor_id", "cell_type"])["n_cells"].unstack(fill_value=0).sort_index()
        counts = counts.loc[:, (counts > 0).sum(axis=0) >= 6]
        metadata = part.groupby("donor_id").agg(age=("age", "first"), sex_male=("sex", "first"))
        metadata["sex_male"] = (metadata["sex_male"].astype(str).str.lower() == "male").astype(float)
        yield "Human", part["group"].iloc[0], cohort, counts, metadata
    from scripts.run_organ_composition import nmr_annotations
    counts, _ = nmr_annotations(0.4)
    metadata = pd.DataFrame(index=counts.index)
    metadata["age"] = metadata.index.map({"1J": 1.1, "5J": 1.1, "6J": 1.1, "3V": 11.3, "4V": 11.3})
    metadata["sex_male"] = 0.0
    yield "Naked mole-rat", "Skin", "GSE152886", counts, metadata


def extended_inputs():
    from scripts.run_afca_aging import afca_obs
    for organ in ("head", "body"):
        atlas = ad.read_h5ad(REPO / "data" / "raw" / "afca" / f"afca-{organ}.h5ad", backed="r")
        obs = afca_obs(atlas)
        obs = obs[(obs["dataset"] == "AFCA") & (obs["cell_type"] != "unannotated")]
        counts = composition_counts(obs, "sample", "cell_type")
        metadata = obs.groupby("sample").agg(age=("age", "first"), sex_male=("sex_male", "first"))
        yield "Fruit fly", organ, "AFCA", counts, metadata
    from scripts.run_gao_worm_aging import normal_obs
    atlas = ad.read_h5ad(REPO / "data" / "interim" / "worm-gao-integrated.h5ad", backed="r")
    obs = normal_obs(atlas)
    counts = composition_counts(obs, "donor", "cell_type")
    metadata = obs.groupby("donor").agg(age=("age", "first"), sex_male=("sex_male", "first"),
                                         platform_nextseq=("platform_nextseq", "first"))
    yield "C. elegans", "Whole animal", "Gao_GSE229022", counts, metadata
    from scripts.run_zebrafish_testis_aging import standard_obs
    atlases = {
        "5-22mo": ad.read_h5ad(REPO / "data" / "interim" / "zebrafish-testis-5-22mo.h5ad", backed="r"),
        "27mo": ad.read_h5ad(REPO / "data" / "interim" / "zebrafish-testis-27mo.h5ad", backed="r"),
    }
    obs = pd.concat([standard_obs(a, name) for name, a in atlases.items()])
    obs = obs[obs["broad_cell_type"] != "unknown"]
    counts = composition_counts(obs, "donor", "broad_cell_type")
    metadata = obs.groupby("donor").agg(age=("age", "first"), sex_male=("sex_male", "first"))
    yield "Zebrafish", "Testis", "GSE275361", counts, metadata


def fit_models(counts, metadata, species, organ, cohort, pseudocount):
    metadata = metadata.loc[counts.index].copy()
    clr = clr_table(counts, pseudocount=pseudocount)
    total = counts.sum(axis=1)
    rows = []
    for cell_type in counts.columns:
        y = clr[cell_type]
        predictors = pd.DataFrame({"age": metadata["age"], "log_total_cells": np.log(total)})
        if metadata["sex_male"].nunique() > 1:
            predictors["sex_male"] = metadata["sex_male"]
        if "platform_nextseq" in metadata and metadata["platform_nextseq"].nunique() > 1:
            predictors["platform_nextseq"] = metadata["platform_nextseq"]
        X = sm.add_constant(predictors.astype(float), has_constant="add")
        if len(y) <= X.shape[1] + 1 or np.linalg.matrix_rank(X.to_numpy()) < X.shape[1]:
            continue
        fit = sm.OLS(y, X).fit()
        ci = fit.conf_int().loc["age"]
        rows.append({
            "species": species, "organ": organ, "cohort": cohort, "cell_type": cell_type,
            "pseudocount": pseudocount, "model": "age_sex_platform_log_total_cells",
            "effect": fit.params["age"], "se": fit.bse["age"], "ci95_low": ci.iloc[0],
            "ci95_high": ci.iloc[1], "pvalue": fit.pvalues["age"], "n_donors": len(y),
            "total_cells": int(counts[cell_type].sum()),
        })
    return rows


def main():
    rows = []
    for inputs in (core_inputs, extended_inputs):
        for species, organ, cohort, counts, metadata in inputs():
            print(f"{species} | {organ} | {cohort}", flush=True)
            for pseudocount in (0.1, 0.5, 1.0):
                rows.extend(fit_models(counts, metadata, species, organ, cohort, pseudocount))
    result = pd.DataFrame(rows)
    result["q_bh"] = np.nan
    for _, indices in result.groupby(["species", "organ", "cohort", "pseudocount"]).groups.items():
        result.loc[indices, "q_bh"] = multipletests(result.loc[indices, "pvalue"], method="fdr_bh")[1]
    primary = result[result["pseudocount"] == 0.5].copy()
    primary["q_global_bh"] = multipletests(primary["pvalue"], method="fdr_bh")[1]
    result["q_global_bh"] = np.nan
    result.loc[primary.index, "q_global_bh"] = primary["q_global_bh"]
    result.to_csv(RESULTS / "reviewer-phase1-composition-total-adjusted.csv", index=False)
    print(result[result["pseudocount"] == 0.5].groupby(["species", "organ", "cohort"]).agg(
        n=("cell_type", "size"), n_q=("q_bh", lambda x: int((x < 0.05).sum())),
        n_global=("q_global_bh", lambda x: int((x < 0.05).sum()))
    ).to_string(), flush=True)


if __name__ == "__main__":
    main()
