#!/usr/bin/env python3
"""Reciprocal BLASTP for SAMap, following the official map_genes.sh:
makeblastdb per proteome, then blastp both directions with
-outfmt 6 -max_hsps 1 -evalue 1e-6.

SPEED CHANGE (Nov 25, 2026, documented, thresholds untouched): the full
proteomes are NOT blasted. The first full-isoform hu->ms run was measured at
~1,000 queries/min: 245,535 human proteins would take ~4 hours for one
direction and ~13 hours for all remaining arms, outside the 1-2 day budget.
Instead each species contributes ONE representative protein per gene present
in its prepped expression file: the longest protein for that gene. SAMap
collapses protein hits to gene level by max bit score (_coarsen_blast_graph),
so the gene-level graph it consumes is the same up to rare isoform-specific
edge cases (a shorter isoform aligning better than the longest); the
validation gate (Stage 1 must recover >= 14 of the 28 name-matched pairs)
empirically tests whether the graph is good enough, and the full-isoform run
remains available as a fallback if it fails. E-value threshold, match rule,
and gates are unchanged.

SCOPE CHANGE (Reza, Nov 25 01:08 PT): zebrafish (dr) and naked mole rat (hg)
are removed from the paper; their arms are dropped. Default pair list covers
only hu, ms, dm, ce. Their proteomes and prepped files stay on disk untouched.

Outputs data/samap/maps/<id1><id2>/<id1>_to_<id2>.txt and the reverse file,
the layout SAMap's _calculate_blast_graph reads. Species IDs are the 2-char
codes SAMap expects (map_genes.sh truncates to 2 chars).

Default runs all 6 unordered pairs among hu/ms/dm/ce. Pass --pairs hu-ms
to run one pair (Stage 1 validation needs only hu-ms first).
"""
from __future__ import annotations

import argparse
import gzip
import os
import re
import subprocess
import sys
from itertools import combinations
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
FASTA_DIR = REPO / "data" / "samap" / "proteomes"
PREP_DIR = REPO / "data" / "samap" / "prepped"
MAP_DIR = REPO / "data" / "samap" / "maps"

FASTA = {
    "hu": "Homo_sapiens.GRCh38.pep.all.fa.gz",
    "ms": "Mus_musculus.GRCm39.pep.all.fa.gz",
    "dm": "Drosophila_melanogaster.BDGP6.54.pep.all.fa.gz",
    "ce": "Caenorhabditis_elegans.WBcel235.pep.all.fa.gz",
    # dropped Nov 25 01:08 PT at Reza's order; proteomes stay on disk, never blasted:
    "dr": "Danio_rerio.GRCz11.pep.all.fa.gz",
    "hg": "Heterocephalus_glaber_female.Naked_mole-rat_maternal.pep.all.fa.gz",
}
SLUG = {"hu": "human", "ms": "mouse", "dm": "fly", "ce": "worm",
        "dr": "zebrafish", "hg": "nmr"}
IDS = ["hu", "ms", "dm", "ce"]  # active species only (dr/hg dropped)
GENE_RE = re.compile(r"gene:(\S+)")


def gunzip_if_needed(gz_path, out_path):
    import shutil
    if out_path.exists():
        return out_path
    with gzip.open(gz_path, "rb") as src, open(out_path, "wb") as dst:
        shutil.copyfileobj(src, dst)
    return out_path


def expressed_genes(sid):
    """Gene IDs present in the prepped expression file for this species."""
    import anndata as ad
    path = PREP_DIR / f"{SLUG[sid]}.h5ad"
    adata = ad.read_h5ad(path, backed="r")
    genes = set(adata.var_names.astype(str))
    adata.file.close()
    return genes


def build_representative_fasta(sid):
    """One FASTA with the longest protein per expressed gene.

    Returns the path. Rebuilt only if missing; the gene set and proteome are
    fixed for this run, so a cached file is valid.
    """
    out = MAP_DIR / f"rep_{sid}.fa"
    if out.exists() and out.stat().st_size > 0:
        print(f"  rep {sid}: exists ({out.stat().st_size} bytes)", flush=True)
        return out
    genes = expressed_genes(sid)
    fa = gunzip_if_needed(FASTA_DIR / FASTA[sid],
                          FASTA_DIR / FASTA[sid].replace(".fa.gz", ".fa"))
    best = {}  # gene_id -> (length, header, sequence)
    header, seq_parts = None, []
    n_proteins = 0

    def flush():
        nonlocal header, seq_parts
        if header is None:
            return
        seq = "".join(seq_parts)
        match = GENE_RE.search(header)
        gene = match.group(1).split(".")[0] if match else ""
        if gene and gene in genes:
            length = len(seq)
            if gene not in best or length > best[gene][0]:
                best[gene] = (length, header, seq)
        header, seq_parts = None, []

    with open(fa) as handle:
        for line in handle:
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                seq_parts = []
                n_proteins += 1
            else:
                seq_parts.append(line.strip())
    flush()
    with open(out, "w") as handle:
        for gene in sorted(best):
            _, header, seq = best[gene]
            handle.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                handle.write(seq[i:i + 60] + "\n")
    print(f"  rep {sid}: {n_proteins} proteins scanned, "
          f"{len(best)} representative genes written (of {len(genes)} expressed)",
          flush=True)
    return out


def run(cmd, log):
    print(" ".join(str(c) for c in cmd), flush=True)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log.write(result.stdout + result.stderr)
        raise RuntimeError(f"command failed: {' '.join(str(c) for c in cmd)}\n"
                           f"{result.stdout[-2000:]}{result.stderr[-2000:]}")
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", nargs="+", default=None,
                        help="e.g. hu-ms dm-hu (default: all 6 unordered pairs)")
    parser.add_argument("--threads", type=int, default=min(16, os.cpu_count() or 8))
    args = parser.parse_args()
    if args.pairs:
        pairs = [tuple(p.split("-")) for p in args.pairs]
        for pair in pairs:
            if any(sid not in IDS for sid in pair):
                sys.exit(f"pair {pair} includes a dropped species (dr/hg); "
                         f"active species are {IDS}")
    else:
        pairs = list(combinations(IDS, 2))
    MAP_DIR.mkdir(parents=True, exist_ok=True)
    rep_paths = {}
    for sid in IDS:
        rep_paths[sid] = build_representative_fasta(sid)
    for id1, id2 in pairs:
        outdir = MAP_DIR / f"{id1}{id2}"
        outdir.mkdir(parents=True, exist_ok=True)
        for sid in (id1, id2):
            db = rep_paths[sid]
            if not Path(str(db) + ".phr").exists():
                run(["makeblastdb", "-in", str(db), "-dbtype", "prot"],
                    log=open(outdir / "makeblastdb.log", "a"))
        for query, db in ((id1, id2), (id2, id1)):
            out = outdir / f"{query}_to_{db}.txt"
            if out.exists() and out.stat().st_size > 0:
                print(f"  skip {out.name} (exists)", flush=True)
                continue
            run(["blastp", "-query", str(rep_paths[query]), "-db", str(rep_paths[db]),
                 "-outfmt", "6", "-out", str(out), "-num_threads", str(args.threads),
                 "-max_hsps", "1", "-evalue", "1e-6"],
                log=open(outdir / "blast.log", "a"))
            print(f"  wrote {out} ({out.stat().st_size} bytes)", flush=True)
    print("BLAST done.", flush=True)


if __name__ == "__main__":
    main()
