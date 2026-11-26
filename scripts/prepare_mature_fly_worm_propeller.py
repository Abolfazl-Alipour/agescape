#!/usr/bin/env python3
"""Build mature propeller inputs for fruit fly and worm (pool x cell-type counts + pool metadata).

Mirrors scripts/prepare_mature_propeller.py, but for the two non-mammalian
strata of the pooled mature fraction scan: AFCA fly body, AFCA fly head, and
Gao worm (GSE229022). The inferential unit is the pool, not the animal.
Mature only: fly 30 days and older, worm 6 days and older.

Layout matches data/interim/mature-propeller: one subdir per cohort with
counts.csv and metadata.csv, plus manifest.csv at the root. Metadata always
carries age and sex_male. Worm metadata also carries platform_nextseq: the
day 6+ worm data spans two sequencing platforms, and the propeller R script
adds platform_nextseq to the design whenever it varies. Worm sex is constant
(all pools female/hermaphrodite), so the R script drops sex_male for worm.
Fly has no platform column in the AFCA obs.
"""

import hashlib
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

OUTPUT = REPO / "data" / "interim" / "mature-propeller-fly-worm"


def write_cohort(root, species, organ, cohort, counts, metadata):
    key = hashlib.sha256(f"{species}|{organ}|{cohort}".encode()).hexdigest()[:16]
    target = root / key
    target.mkdir(parents=True, exist_ok=True)
    counts.to_csv(target / "counts.csv")
    metadata.loc[counts.index].to_csv(target / "metadata.csv")
    return {"key": key, "species": species, "organ": organ, "cohort": cohort,
            "n_pools": len(counts), "n_cell_types": counts.shape[1]}


def main():
    import anndata as ad

    from scripts.run_gao_worm_aging import normal_obs as worm_normal_obs
    from scripts.run_mature_all_species_fractions import load_strata

    OUTPUT.mkdir(parents=True, exist_ok=True)
    manifest = []
    for species, organ, cohort, counts, metadata in load_strata():
        if species not in ("Fruit fly", "C. elegans"):
            continue
        meta = metadata.copy()
        meta["sex_male"] = meta["sex"].astype(float)
        if species == "C. elegans":
            obs = worm_normal_obs(ad.read_h5ad(
                REPO / "data" / "interim" / "worm-gao-integrated.h5ad", backed="r"))
            plat = obs[obs["age"] >= 6].groupby("donor")["platform_nextseq"].first()
            meta["platform_nextseq"] = plat.loc[meta.index].astype(float)
        meta = meta[[c for c in ("age", "sex_male", "platform_nextseq") if c in meta.columns]]
        counts = counts.astype(int)
        assert set(counts.index) == set(meta.index), f"pool mismatch {species} {organ} {cohort}"
        manifest.append(write_cohort(OUTPUT, species, organ, cohort, counts, meta))

    pd.DataFrame(manifest).to_csv(OUTPUT / "manifest.csv", index=False)
    print(f"wrote {len(manifest)} cohorts to {OUTPUT}", flush=True)
    for row in manifest:
        print(row, flush=True)


if __name__ == "__main__":
    main()
