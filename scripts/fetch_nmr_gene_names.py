#!/usr/bin/env python3
"""Fetch ENSHGLG (Heterocephalus glaber) Ensembl gene ID -> display_name mapping
from Ensembl REST (batch POST), for genes present in GSE152886. Cached TSV.

Source: https://rest.ensembl.org/lookup/id (Ensembl release live at query date;
recorded in output header). One POST per 1000 ids. Retries with backoff.
"""
import gzip, json, sys, time, urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
RAW = REPO / "data/raw/GSE152886_nmr_skin_raw.csv.gz"
OUT = REPO / "data/interim/nmr-ensg2symbol.tsv"

# gene ids = first column (skip header)
ids = []
with gzip.open(RAW, "rt") as f:
    header = f.readline()
    for line in f:
        ids.append(line.split("\t", 1)[0].strip('"'))
print(f"{len(ids)} gene ids in GSE152886", flush=True)

mapping = {}
B = 1000
for i in range(0, len(ids), B):
    chunk = ids[i:i+B]
    body = json.dumps({"ids": chunk}).encode()
    req = urllib.request.Request(
        "https://rest.ensembl.org/lookup/id",
        data=body, headers={"Content-Type": "application/json"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                res = json.loads(r.read().decode())
            break
        except Exception as e:
            print(f"  chunk {i}: attempt {attempt} failed: {e}", flush=True)
            time.sleep(5 * (attempt + 1))
    else:
        sys.exit(f"chunk {i} failed permanently")
    for gid, rec in res.items():
        mapping[gid] = (rec or {}).get("display_name", "")
    print(f"  mapped {len(mapping)}/{len(ids)}", flush=True)
    time.sleep(0.4)

with open(OUT, "w") as f:
    f.write("# ensembl_gene_id\tdisplay_name\t# source=Ensembl REST lookup/id, queried 2026-10-14\n")
    for gid in ids:
        f.write(f"{gid}\t{mapping.get(gid, '')}\n")
n_named = sum(1 for v in mapping.values() if v)
print(f"DONE: {len(mapping)} mapped, {n_named} with display_name -> {OUT}")
