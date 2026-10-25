#!/usr/bin/env python3
"""
build_reference_orthologs.py - builds the committed reference tables in data/reference/
from Ensembl Release 106 (the release whose naked mole-rat gene IDs, ENSHGLG00000000001-
ENSHGLG00000031858, match GSE152886).

Sources (downloaded to a scratch dir, never committed):
  gtf   https://ftp.ensembl.org/pub/release-106/gtf/mus_musculus/Mus_musculus.GRCm39.106.gtf.gz
  homol https://ftp.ensembl.org/pub/release-106/tsv/ensembl-compara/homologies/mus_musculus/
        Compara.106.protein_default.homologies.tsv.gz

Outputs (data/reference/):
  orthologs-mouse-nmr.csv         mouse -> Heterocephalus glaber (female, HetGla_female_1.0), 1:1 only
  mouse-gene-chromosomes.csv      mouse gene id, symbol, chromosome, gene_biotype
  MANIFEST.json                   SHA-256 of every reference file and source URL
Existing orthologs-mouse-human.csv, nmr-ensg2symbol.tsv, m1-hvg-union.csv are copied from
data/interim/ by --copy-existing.
"""
import argparse
import csv
import gzip
import hashlib
import json
import re
import shutil
from pathlib import Path

GTF_URL = "https://ftp.ensembl.org/pub/release-106/gtf/mus_musculus/Mus_musculus.GRCm39.106.gtf.gz"
HOM_URL = ("https://ftp.ensembl.org/pub/release-106/tsv/ensembl-compara/homologies/mus_musculus/"
           "Compara.106.protein_default.homologies.tsv.gz")
NMR_SPECIES = "heterocephalus_glaber_female"

REPO = Path(__file__).resolve().parent.parent
REF = REPO / "data" / "reference"


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_gtf(gtf: Path):
    rx = {k: re.compile(k + r' "([^"]+)"') for k in ("gene_id", "gene_name", "gene_biotype")}
    genes = {}
    with gzip.open(gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if c[2] != "gene":
                continue
            attrs = {k: (m.group(1) if (m := r.search(c[8])) else "") for k, r in rx.items()}
            genes[attrs["gene_id"]] = (attrs["gene_name"] or attrs["gene_id"], c[0], attrs["gene_biotype"])
    return genes


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scratch", default=str(REPO / "scratch_dl"))
    ap.add_argument("--copy-existing", action="store_true")
    args = ap.parse_args()
    scratch = Path(args.scratch)
    gtf = scratch / "Mus_musculus.GRCm39.106.gtf.gz"
    hom = scratch / "Compara.106.protein_default.homologies.tsv.gz"
    REF.mkdir(parents=True, exist_ok=True)

    genes = parse_gtf(gtf)
    with open(REF / "mouse-gene-chromosomes.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["mouse_ensembl_id", "mouse_gene_symbol", "chromosome_name", "gene_biotype"])
        for gid in sorted(genes):
            sym, chrom, bt = genes[gid]
            w.writerow([gid, sym, chrom, bt])

    n = 0
    with gzip.open(hom, "rt") as f, open(REF / "orthologs-mouse-nmr.csv", "w", newline="") as out:
        w = csv.writer(out)
        w.writerow(["mouse_ensembl_id", "mouse_gene_symbol", "nmr_ensembl_id", "orthology_type"])
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            if r["homology_species"] != NMR_SPECIES or r["homology_type"] != "ortholog_one2one":
                continue
            gid = r["gene_stable_id"]
            w.writerow([gid, genes.get(gid, (gid,))[0], r["homology_gene_stable_id"], "ortholog_one2one"])
            n += 1
    print(f"1:1 mouse-NMR orthologs: {n}")

    if args.copy_existing:
        for name in ("orthologs-mouse-human.csv", "nmr-ensg2symbol.tsv", "m1-hvg-union.csv"):
            src = REPO / "data" / "interim" / name
            if src.exists():
                shutil.copy2(src, REF / name)

    files = {}
    for p in sorted(REF.iterdir()):
        if p.name != "MANIFEST.json" and p.is_file():
            files[p.name] = {"sha256": sha256(p), "size_bytes": p.stat().st_size}
    manifest = {"ensembl_release": 106, "sources": {"gtf": GTF_URL, "homologies": HOM_URL,
                "nmr_species": NMR_SPECIES}, "files": files}
    (REF / "MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("wrote", REF / "MANIFEST.json")


if __name__ == "__main__":
    main()
