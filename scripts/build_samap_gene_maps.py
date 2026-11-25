#!/usr/bin/env python3
"""Harmonize gene IDs between SAMap inputs and proteome FASTAs.

For each species, parses its proteome FASTA (all Ensembl-style headers with
gene: and gene_symbol: fields), renames the input AnnData var index to the
stable gene ID (ENSG / ENSMUSG / FBgn / WBGene / ENSDARG / ENSHGLG), and
writes:
  - data/samap/prepped/<slug>.h5ad  (SAMap-ready AnnData, var = gene IDs)
  - data/samap/maps/names_<sid>.csv (FASTA header -> gene ID, for SAMap `names`)
  - data/samap/gene_map_counts.csv (matched-ID counts per species)

A var gene counts as matched if its gene ID appears in the FASTA (i.e. the
protein exists). Unmatched var genes are dropped from the prepped file and
listed in the counts table. Duplicate var names after renaming keep the copy
with the highest total counts; the drop is reported.
"""
from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
IN_DIR = REPO / "data" / "samap" / "inputs"
PREP_DIR = REPO / "data" / "samap" / "prepped"
MAP_DIR = REPO / "data" / "samap" / "maps"
FASTA_DIR = REPO / "data" / "samap" / "proteomes"

# sid -> (input slug, fasta file, var id kind)
SPECIES = {
    "hu": ("human", "Homo_sapiens.GRCh38.pep.all.fa.gz", "ensg"),
    "ms": ("mouse", "Mus_musculus.GRCm39.pep.all.fa.gz", "symbol"),
    "dm": ("fly", "Drosophila_melanogaster.BDGP6.54.pep.all.fa.gz", "symbol"),
    "ce": ("worm", "Caenorhabditis_elegans.WBcel235.pep.all.fa.gz", "symbol_or_seq"),
    "dr": ("zebrafish", "Danio_rerio.GRCz11.pep.all.fa.gz", "symbol"),
    "hg": ("nmr", "Heterocephalus_glaber_female.Naked_mole-rat_maternal.pep.all.fa.gz", "ensg"),
}

GENE_RE = re.compile(r"gene:(\S+)")
SYMBOL_RE = re.compile(r"gene_symbol:(\S+)")


def strip_version(gene_id):
    return gene_id.split(".")[0]


def parse_fasta(path, split_colon=False):
    """Return (names_df, proteins). names_df: protein_id, gene_id, symbol.

    split_colon=True only for fly (dm): its gene_symbol field is
    'His2A:CG33832' (symbol:FlyBase ID). For other species the symbol is
    used whole; zebrafish symbols themselves contain colons
    (si:ch211-..., zgc:...), so splitting would corrupt them."""
    rows = []
    with gzip.open(path, "rt") as handle:
        header = None
        for line in handle:
            if line.startswith(">"):
                header = line[1:].strip()
                first = header.split()[0]
                gene = GENE_RE.search(header)
                symbol = SYMBOL_RE.search(header)
                sym = symbol.group(1) if symbol else ""
                if split_colon and sym:
                    sym = sym.split(":")[0]
                rows.append({
                    "protein_id": first,
                    "gene_id": strip_version(gene.group(1)) if gene else "",
                    "symbol": sym,
                })
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", choices=list(SPECIES), nargs="+", default=list(SPECIES))
    args = parser.parse_args()
    PREP_DIR.mkdir(parents=True, exist_ok=True)
    MAP_DIR.mkdir(parents=True, exist_ok=True)
    import anndata as ad
    count_rows = []
    for sid in args.species:
        slug, fasta_name, kind = SPECIES[sid]
        fasta = parse_fasta(FASTA_DIR / fasta_name, split_colon=(sid == "dm"))
        print(f"{sid} ({slug}): {len(fasta)} proteins in FASTA", flush=True)
        # gene -> symbol lookup (first symbol seen per gene) and symbol -> gene
        gene_symbol = (fasta[fasta["symbol"] != ""]
                       .drop_duplicates("gene_id").set_index("gene_id")["symbol"])
        symbol_gene = (fasta[fasta["symbol"] != ""]
                       .drop_duplicates("symbol").set_index("symbol")["gene_id"])
        # sequence-name lookup for worm: protein id minus isoform suffix
        seq_gene = {}
        if kind == "symbol_or_seq":
            for row in fasta.itertuples():
                seq = row.protein_id.rsplit(".", 1)[0] if "." in row.protein_id else row.protein_id
                seq_gene.setdefault(seq, row.gene_id)
        # names map: FASTA header -> gene ID (what SAMap `names` needs)
        names = fasta[["protein_id", "gene_id"]].copy()
        names = names[names["gene_id"] != ""]
        names.to_csv(MAP_DIR / f"names_{sid}.csv", index=False)
        print(f"  names map: {len(names)} protein->gene rows", flush=True)

        adata = ad.read_h5ad(IN_DIR / f"{slug}.h5ad")
        # cluster key = the concordance unit (tissue/cohort + cell type)
        adata.obs["cluster"] = (adata.obs["tissue"].astype(str) + " | "
                                + adata.obs["cohort"].astype(str) + " | "
                                + adata.obs["unit"].astype(str))
        var = adata.var_names.astype(str)
        symbols = (adata.var["symbol"].astype(str).to_numpy()
                   if "symbol" in adata.var.columns else np.array([""] * len(var)))
        gene_ids = np.array([""] * len(var), dtype=object)
        if kind == "ensg":
            gene_ids = np.array([strip_version(v) for v in var], dtype=object)
        elif kind == "symbol":
            gene_ids = np.array([symbol_gene.get(v, "") for v in var], dtype=object)
        elif kind == "symbol_or_seq":
            for i, (v, s) in enumerate(zip(var, symbols)):
                gene_ids[i] = symbol_gene.get(v, "") or seq_gene.get(v, "") or symbol_gene.get(s, "")
        fasta_genes = set(fasta["gene_id"]) - {""}
        matched = np.array([g in fasta_genes for g in gene_ids])
        print(f"  var genes: {len(var)}, mapped to a gene ID: {int((gene_ids != '').sum())}, "
              f"matched in FASTA: {int(matched.sum())} "
              f"({100 * matched.mean():.1f}%)", flush=True)
        # drop unmatched, dedup by gene id keeping highest total counts
        keep = np.where(matched)[0]
        totals = np.asarray(adata.X[:, keep].sum(axis=0)).ravel() if len(keep) else np.array([])
        sub = adata[:, keep].copy()
        gene_keep = gene_ids[keep]
        order = np.argsort(-totals, kind="stable")
        seen, dedup_idx = set(), []
        for pos in order:
            g = gene_keep[pos]
            if g not in seen:
                seen.add(g)
                dedup_idx.append(pos)
        dedup_idx = np.sort(np.array(dedup_idx))
        dropped_dup = len(keep) - len(dedup_idx)
        sub = sub[:, dedup_idx].copy()
        sub.var_names = gene_keep[dedup_idx]
        sub.var["symbol"] = [gene_symbol.get(g, "") for g in sub.var_names]
        sub.write_h5ad(PREP_DIR / f"{slug}.h5ad", compression="gzip")
        print(f"  prepped: {sub.n_obs} cells x {sub.n_vars} genes "
              f"(dropped {int((~matched).sum())} unmatched, {dropped_dup} duplicate-symbol)",
              flush=True)
        count_rows.append({
            "sid": sid, "species": slug, "fasta_proteins": int(len(fasta)),
            "fasta_genes": int(len(fasta_genes)),
            "var_genes": int(len(var)),
            "mapped_to_gene_id": int((gene_ids != "").sum()),
            "matched_in_fasta": int(matched.sum()),
            "pct_matched": round(100 * float(matched.mean()), 2),
            "prepped_genes": int(sub.n_vars),
            "prepped_cells": int(sub.n_obs),
            "dropped_duplicate_symbol": int(dropped_dup),
        })
    counts = pd.DataFrame(count_rows)
    counts_path = REPO / "data" / "samap" / "gene_map_counts.csv"
    if counts_path.exists():
        old = pd.read_csv(counts_path)
        counts = pd.concat([old[~old["sid"].isin(args.species)], counts], ignore_index=True)
    counts.to_csv(counts_path, index=False)
    print(counts.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
