#!/usr/bin/env python3
"""Run SAMap for one species pair and apply the pre-stated match rule.

Match rule (all must hold, fixed before seeing results):
  1. alignment score >= 0.30
  2. reciprocal best hit: each cluster's highest-scoring partner is the other
  3. >= 40 enriched gene pairs (GenePairFinder, align_thr 0.1)
  4. score in the top 1% of all pairs in that comparison

Clusters are the concordance units (obs['cluster'] = tissue | cohort | unit).
Outputs (results/samap-<id1><id2>/):
  - alignment_scores.csv (full pairwise alignment matrix)
  - gene_pairs.csv (enriched gene pairs per mapped pair)
  - gene_pair_counts.csv (n enriched pairs per pair)
  - matches.csv (pairs passing the full rule)
  - run_summary.json (cells, genes matched, convergence, counts)
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
PREP_DIR = REPO / "data" / "samap" / "prepped"
MAP_DIR = REPO / "data" / "samap" / "maps"
RESULTS = REPO / "results"

SLUG = {"hu": "human", "ms": "mouse", "dm": "fly",
        "ce": "worm", "dr": "zebrafish", "hg": "nmr"}
ALIGN_THR = 0.30
MIN_GENE_PAIRS = 40
TOP_PCT = 1.0


def load_names(sid):
    table = pd.read_csv(MAP_DIR / f"names_{sid}.csv", dtype=str)
    return list(zip(table["protein_id"], table["gene_id"]))


def apply_match_rule(scores, pair_counts):
    """scores: DataFrame indexed/columned by 'sid_cluster'. Returns matches df."""
    id1, id2 = None, None
    # split index labels into (sid, cluster)
    rows = []
    for a in scores.index:
        for b in scores.columns:
            sid_a, cl_a = a.split("_", 1)
            sid_b, cl_b = b.split("_", 1)
            if sid_a == sid_b:
                continue
            rows.append((sid_a, cl_a, sid_b, cl_b, float(scores.loc[a, b])))
    df = pd.DataFrame(rows, columns=["sid1", "cluster1", "sid2", "cluster2", "score"])
    if df.empty:
        return df
    # reciprocal best hit within the cross-species block, per direction
    best1 = df.loc[df.groupby(["sid1", "cluster1"])["score"].idxmax()]
    best2 = df.loc[df.groupby(["sid2", "cluster2"])["score"].idxmax()]
    best_pairs = set(map(tuple, best1[["sid1", "cluster1", "sid2", "cluster2"]].values)) & \
        set(map(tuple, best2[["sid1", "cluster1", "sid2", "cluster2"]].values))
    df["reciprocal_best"] = [tuple(r) in best_pairs for r in
                             df[["sid1", "cluster1", "sid2", "cluster2"]].values]
    # top 1% of scores in this comparison
    thr_top = float(np.quantile(df["score"], 1 - TOP_PCT / 100.0))
    df["top1pct"] = df["score"] >= thr_top
    df["top1pct_thr"] = thr_top
    # enriched gene pairs (keyed "sid_cluster;sid_cluster" in GenePairFinder output)
    def n_pairs(row):
        key1 = f"{row['sid1']}_{row['cluster1']};{row['sid2']}_{row['cluster2']}"
        key2 = f"{row['sid2']}_{row['cluster2']};{row['sid1']}_{row['cluster1']}"
        return int(pair_counts.get(key1, pair_counts.get(key2, 0)))
    df["n_gene_pairs"] = df.apply(n_pairs, axis=1)
    df["passes"] = ((df["score"] >= ALIGN_THR) & df["reciprocal_best"]
                    & (df["n_gene_pairs"] >= MIN_GENE_PAIRS) & df["top1pct"])
    # one row per unordered pair (both orientations carry the same score)
    df["pair_key"] = [tuple(sorted([(a, b), (c, d)])) for a, b, c, d in
                      df[["sid1", "cluster1", "sid2", "cluster2"]].values]
    df = df.drop_duplicates("pair_key").drop(columns="pair_key")
    return df


def main():
    global ALIGN_THR
    parser = argparse.ArgumentParser()
    parser.add_argument("id1", choices=list(SLUG))
    parser.add_argument("id2", choices=list(SLUG))
    parser.add_argument("--outdir", default=None)
    parser.add_argument("--align-thr", type=float, default=ALIGN_THR)
    parser.add_argument("--no-gene-pairs", action="store_true",
                        help="skip GenePairFinder (alignment scores only)")
    parser.add_argument("--override", action="append", default=[],
                        metavar="SID:PATH",
                        help="use a different prepped file for a species, e.g. dm:path/to/fly-perm.h5ad")
    args = parser.parse_args()
    ALIGN_THR = args.align_thr
    outdir = Path(args.outdir) if args.outdir else RESULTS / f"samap-{args.id1}{args.id2}"
    outdir.mkdir(parents=True, exist_ok=True)

    from samap import SAMAP, GenePairFinder, get_mapping_scores

    names = {sid: load_names(sid) for sid in (args.id1, args.id2)}
    sams = {sid: str(PREP_DIR / f"{SLUG[sid]}.h5ad") for sid in (args.id1, args.id2)}
    for item in args.override:
        sid, path = item.split(":", 1)
        sams[sid] = path
    keys = {args.id1: "cluster", args.id2: "cluster"}
    print(f"Running SAMap {args.id1} vs {args.id2} ...", flush=True)
    sm = SAMAP(sams=sams, f_maps=str(MAP_DIR) + "/", names=names, keys=keys,
               save_processed=False)
    print(f"  genes matched: " + ", ".join(
        f"{sid}={len(sm.gns_dict[sid])}" for sid in (args.id1, args.id2)), flush=True)
    sm.run(umap=False)  # UMAP is only a visualization projection; the mapping
    # scores and gene pairs do not use it. Skipping it saves substantial
    # time and memory on 100k+ cell inputs.
    print("  mapping converged.", flush=True)
    _, scores = get_mapping_scores(sm, keys)
    scores.to_csv(outdir / "alignment_scores.csv")
    pair_counts = {}
    if not args.no_gene_pairs:
        gpf = GenePairFinder(sm, keys=keys)
        gp = gpf.find_all(align_thr=0.1)
        gp.to_csv(outdir / "gene_pairs.csv")
        for col in gp.columns:
            if col.endswith("_pval1") or col.endswith("_pval2"):
                continue
            pair_counts[col] = int(gp[col].notna().sum())
        pd.Series(pair_counts, name="n_gene_pairs").to_csv(
            outdir / "gene_pair_counts.csv")
    matches = apply_match_rule(scores, pair_counts)
    matches.to_csv(outdir / "matches.csv", index=False)
    n_pass = int(matches["passes"].sum()) if len(matches) else 0
    print(f"  pairs scored: {len(matches)}, passing full rule: {n_pass}", flush=True)
    summary = {
        "id1": args.id1, "id2": args.id2,
        "n_cells": {sid: int(sm.sams[sid].adata.n_obs) for sid in (args.id1, args.id2)},
        "n_clusters": {sid: int(sm.sams[sid].adata.obs["cluster"].nunique())
                       for sid in (args.id1, args.id2)},
        "n_genes_matched": {sid: int(len(sm.gns_dict[sid])) for sid in (args.id1, args.id2)},
        "n_pairs_scored": int(len(matches)),
        "n_passing": n_pass,
        "align_thr": ALIGN_THR, "min_gene_pairs": MIN_GENE_PAIRS, "top_pct": TOP_PCT,
    }
    (outdir / "run_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
