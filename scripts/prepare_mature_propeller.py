#!/usr/bin/env python3
"""Build mature propeller inputs (donor x cell-type counts + donor metadata).

Mature only: mice 18 months and older, human donors 30 years and older.
Layout mirrors scripts/prepare_reviewer_phase1_propeller.py: one subdir per
cohort with counts.csv and metadata.csv, plus manifest.csv at the root.

HRCA retina suites are NOT built here: the HRCA aggregate files tripped an
execution-risk check and need the user's explicit go-ahead before loading.
"""

import hashlib
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

OUTPUT = REPO / "data" / "interim" / "mature-propeller"


def write_cohort(root, species, organ, cohort, counts, metadata):
    key = hashlib.sha256(f"{species}|{organ}|{cohort}".encode()).hexdigest()[:16]
    target = root / key
    target.mkdir(parents=True, exist_ok=True)
    counts.to_csv(target / "counts.csv")
    metadata.loc[counts.index].to_csv(target / "metadata.csv")
    return {"key": key, "species": species, "organ": organ, "cohort": cohort,
            "n_donors": len(counts), "n_cell_types": counts.shape[1]}


def main():
    from scripts.run_mature_composition_complete import human_inputs, mouse_inputs

    OUTPUT.mkdir(parents=True, exist_ok=True)
    manifest = []
    for species, organ, cohort, counts, metadata in list(mouse_inputs()) + list(human_inputs()):
        meta = metadata.copy()
        meta["sex_male"] = meta["sex"].astype(float)
        meta = meta[["age", "sex_male"]]
        counts = counts.astype(int)
        assert set(counts.index) == set(meta.index), f"donor mismatch {species} {organ} {cohort}"
        manifest.append(write_cohort(OUTPUT, species, organ, cohort, counts, meta))

    results = REPO / "results"
    for cohort in ("GSE214546", "GSE213516"):
        comp = pd.read_csv(results / f"followup-mature-pbmc-{cohort}-composition-donor-counts.csv")
        hell = pd.read_csv(results / f"followup-mature-pbmc-{cohort}-hellinger-donors.csv")
        wide = comp.pivot(index="donor", columns="cell_type", values="n_cells").fillna(0).astype(int)
        dmeta = hell.groupby("donor").agg(age=("age", "first"), sex=("sex", "first"))
        assert set(wide.index) == set(dmeta.index), f"PBMC donor mismatch {cohort}"
        dmeta["sex_male"] = (dmeta["sex"].astype(str).str.lower() == "male").astype(float)
        dmeta = dmeta[["age", "sex_male"]]
        wide = wide.reindex(index=dmeta.index).fillna(0).astype(int)
        manifest.append(write_cohort(OUTPUT, "Human", "blood", f"PBMC_{cohort}", wide, dmeta))

    pd.DataFrame(manifest).to_csv(OUTPUT / "manifest.csv", index=False)
    print(f"wrote {len(manifest)} cohorts to {OUTPUT}", flush=True)
    for row in manifest:
        print(row, flush=True)


if __name__ == "__main__":
    main()
