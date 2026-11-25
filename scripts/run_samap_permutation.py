#!/usr/bin/env python3
"""Label-permutation negative control for a SAMap pair run.

Shuffles the cluster (concordance-unit) labels of one species' prepped file,
re-runs scripts/run_samap_pair.py on the permuted file, and reports how many
pairs pass the full match rule. Under the null (no real homology signal) the
expectation is zero passing pairs.

Usage:
  python scripts/run_samap_permutation.py ms dm --permute dm --seed 20261125
  python scripts/run_samap_permutation.py hu dm --permute dm --seed 20261125

The permuted file is written to data/samap/prepped/<slug>-perm.h5ad and the
results to results/samap-<id1><id2>-perm/. The original prepped files are
never modified.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
PREP_DIR = REPO / "data" / "samap" / "prepped"
SLUG = {"hu": "human", "ms": "mouse", "dm": "fly",
        "ce": "worm", "dr": "zebrafish", "hg": "nmr"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("id1", choices=list(SLUG))
    parser.add_argument("id2", choices=list(SLUG))
    parser.add_argument("--permute", required=True, choices=list(SLUG),
                        help="species whose cluster labels get shuffled")
    parser.add_argument("--seed", type=int, default=20261125)
    args = parser.parse_args()

    import anndata as ad
    src = PREP_DIR / f"{SLUG[args.permute]}.h5ad"
    dst = PREP_DIR / f"{SLUG[args.permute]}-perm.h5ad"
    print(f"Permuting cluster labels in {src} (seed {args.seed}) ...", flush=True)
    adata = ad.read_h5ad(src)
    rng = np.random.default_rng(args.seed)
    labels = adata.obs["cluster"].to_numpy()
    adata.obs["cluster"] = rng.permutation(labels)
    adata.write_h5ad(dst, compression="gzip")
    print(f"  wrote {dst}", flush=True)

    outdir = REPO / "results" / f"samap-{args.id1}{args.id2}-perm"
    cmd = [sys.executable, str(REPO / "scripts" / "run_samap_pair.py"),
           args.id1, args.id2, "--outdir", str(outdir),
           "--override", f"{args.permute}:{dst}"]
    print("Running:", " ".join(cmd), flush=True)
    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        print("PERMUTATION RUN FAILED", flush=True)
        sys.exit(result.returncode)
    import json
    summary = json.loads((outdir / "run_summary.json").read_text())
    print(f"Permutation control {args.id1}-{args.id2} (permuted {args.permute}): "
          f"{summary['n_passing']} passing pairs "
          f"(expect 0 under the null)", flush=True)


if __name__ == "__main__":
    main()
