#!/usr/bin/env python3

import hashlib
from pathlib import Path
import sys

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_reviewer_phase1_composition import core_inputs, extended_inputs

OUTPUT = REPO / "data" / "interim" / "reviewer-phase1-propeller"


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    manifest = []
    for loaders in (core_inputs, extended_inputs):
        for species, organ, cohort, counts, metadata in loaders():
            key = hashlib.sha256(f"{species}|{organ}|{cohort}".encode()).hexdigest()[:16]
            target = OUTPUT / key
            target.mkdir(exist_ok=True)
            counts.to_csv(target / "counts.csv")
            metadata.loc[counts.index].to_csv(target / "metadata.csv")
            manifest.append({"key": key, "species": species, "organ": organ, "cohort": cohort,
                             "n_donors": len(counts), "n_cell_types": counts.shape[1]})
    pd.DataFrame(manifest).to_csv(OUTPUT / "manifest.csv", index=False)


if __name__ == "__main__":
    main()
