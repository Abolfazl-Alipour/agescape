#!/usr/bin/env python3

import itertools
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import requests
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
INTERIM = REPO / "data" / "interim"


def fetch_go(dataset, filename):
    INTERIM.mkdir(exist_ok=True)
    path = INTERIM / filename
    service = "FlyEnrichr" if dataset == "dmelanogaster_gene_ensembl" else "Enrichr"
    library = "GO_Biological_Process_2018" if service == "FlyEnrichr" else "GO_Biological_Process_2023"
    url = (f"https://maayanlab.cloud/{service}/geneSetLibrary?mode=text"
           f"&libraryName={library}")
    if not path.exists():
        response = requests.get(url, timeout=180)
        response.raise_for_status()
        path.write_text(response.text)
    rows = []
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        term = fields[0]
        match = pd.Series([term]).str.extract(r"\((GO:\d+)\)").iloc[0, 0]
        go_id = term if pd.isna(match) else match
        rows.extend((gene.upper(), go_id, term, "biological_process")
                    for gene in fields[2:] if gene)
    return pd.DataFrame(rows, columns=["gene", "go_id", "term", "namespace"])


def enrichment(gene_lists, go, analysis, universe_genes):
    go = go[go["gene"].isin(set(universe_genes))]
    universe = set(go["gene"])
    term_sizes = go.groupby(["go_id", "term"])["gene"].nunique().reset_index(name="n_universe")
    eligible = term_sizes[(term_sizes["n_universe"] >= 5) & (term_sizes["n_universe"] <= 500)]
    term_gene_sets = go.groupby(["go_id", "term"])["gene"].agg(set).to_dict()
    rows = []
    for list_name, genes in gene_lists.items():
        selected = set(genes) & universe
        for term in eligible.itertuples():
            term_genes = term_gene_sets[(term.go_id, term.term)]
            overlap = len(selected & term_genes)
            if overlap == 0:
                continue
            rows.append({
                "analysis": analysis, "list_name": list_name, "go_id": term.go_id, "term": term.term,
                "n_list": len(selected), "n_universe": len(universe), "n_term": term.n_universe,
                "overlap": overlap, "expected": len(selected) * term.n_universe / len(universe),
                "fold_enrichment": overlap / (len(selected) * term.n_universe / len(universe)),
                "pvalue": hypergeom.sf(overlap - 1, len(universe), term.n_universe, len(selected)),
                "genes": ";".join(sorted(selected & term_genes)),
            })
    result = pd.DataFrame(rows)
    if not result.empty:
        result["q_bh"] = multipletests(result["pvalue"], method="fdr_bh")[1]
    return result


def fly_overlap(go):
    path = RESULTS / "reviewer-phase1-fly-gene-effects.csv"
    effects = pd.read_csv(path)
    vector = pd.read_csv(RESULTS / "extended-afca-vector-tests.csv")
    hits = vector[(vector["model"] == "depth_adjusted") & (vector["q_bh"] < 0.05)].copy()
    hit_units = set(zip(hits["tissue"], hits["cell_type"]))
    effects["unit_key"] = list(zip(effects["tissue"], effects["cell_type"]))
    selected = effects[effects["unit_key"].isin(hit_units)]
    gene_lists = {f"{tissue} | {cell_type}": set(part["gene"])
                  for (tissue, cell_type), part in selected.groupby(["tissue", "cell_type"])}
    rows = []
    for first, second in itertools.combinations(sorted(hit_units), 2):
        a = gene_lists.get(" | ".join(first), set())
        b = gene_lists.get(" | ".join(second), set())
        overlap = len(a & b)
        rows.append({
            "unit_1": " | ".join(first), "unit_2": " | ".join(second),
            "n_1": len(a), "n_2": len(b), "overlap": overlap,
            "jaccard": overlap / len(a | b) if a | b else np.nan,
            "expected_overlap": len(a) * len(b) / 15992,
            "hypergeometric_p": hypergeom.sf(overlap - 1, 15992, len(b), len(a)) if overlap else np.nan,
        })
    pairs = pd.DataFrame(rows)
    if pairs["overlap"].notna().any():
        valid = pairs["overlap"] > 0
        pairs.loc[valid, "q_bh"] = multipletests(pairs.loc[valid, "hypergeometric_p"], method="fdr_bh")[1]
    pairs.to_csv(RESULTS / "reviewer-phase1-fly-hit-gene-overlap.csv", index=False)
    import anndata as ad
    atlas = ad.read_h5ad(REPO / "data" / "raw" / "afca" / "afca-head.h5ad", backed="r")
    return enrichment(gene_lists, go, "fly_vector_hit_gene_sets", atlas.var_names)


def main():
    human_go = fetch_go("hsapiens_gene_ensembl", "reviewer-human-go.tsv")
    fly_go = fetch_go("dmelanogaster_gene_ensembl", "reviewer-fly-go.tsv")
    significant = pd.read_csv(RESULTS / "h3-human-neuronal-gene-replication-significant.csv")
    significant = significant[significant["mssm_fdr"]]
    gene_lists = {name: set(part["gene"]) for name, part in significant.groupby("cell_type")}
    shared = pd.read_csv(RESULTS / "h3-human-mssm-shared-neuronal-genes.csv")
    gene_lists["shared_all_three"] = set(shared[shared["n_fdr_subtypes"] == 3]["gene"])
    all_genes = pd.read_csv(RESULTS / "h3-human-neuronal-gene-replication-all.csv", usecols=["gene"])
    interneuron = enrichment(gene_lists, human_go, "interneuron_mssm_fdr", all_genes["gene"].unique())
    interneuron.to_csv(RESULTS / "reviewer-phase1-interneuron-go-enrichment.csv", index=False)
    fly = fly_overlap(fly_go)
    fly.to_csv(RESULTS / "reviewer-phase1-fly-go-enrichment.csv", index=False)


if __name__ == "__main__":
    main()
