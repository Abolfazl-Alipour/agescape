#!/usr/bin/env python3
"""Export per-donor centered-log-ratio fraction values (export only).

Calls the frozen mature composition pipeline's own input loaders and CLR
function (scripts/run_mature_composition_complete.py) and writes one row per
donor per cell type. No regression is rerun here and no analysis changes:
the same donor counts after the same at-least-6-donor presence filter, the
same pseudocount 0.5 CLR transform. The committed regression slopes are
refit from this file by the figure builder and must match the committed
composition tables.
"""

import sys
from pathlib import Path

import pandas as pd

SCRIPTS = Path(__file__).resolve().parent
REPO = SCRIPTS.parent
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import run_mature_composition_complete as comp  # noqa: E402


def main():
    rows = []
    for loader in (comp.mouse_inputs, comp.human_inputs):
        for species, organ, cohort, all_counts, metadata in loader():
            present = (all_counts > 0).sum(axis=0)
            counts = all_counts.loc[:, present >= 6].copy()
            clr = comp.clr_table(counts)
            totals = counts.sum(axis=1)
            for cell_type in counts.columns:
                frame = pd.DataFrame({
                    "species": species,
                    "organ": organ,
                    "cohort": cohort,
                    "donor_id": counts.index.astype(str),
                    "age": metadata.loc[counts.index, "age"].to_numpy(dtype=float),
                    "sex_male": metadata.loc[counts.index, "sex"].to_numpy(dtype=float),
                    "cell_type": cell_type,
                    "clr_value": clr[cell_type].to_numpy(dtype=float),
                    "raw_count": counts[cell_type].to_numpy(dtype=float),
                    "total_cells": totals.to_numpy(dtype=float),
                })
                rows.append(frame)
    out = pd.concat(rows, ignore_index=True)
    out_path = REPO / "results" / "followup-mature-composition-donor-clr.csv"
    out.to_csv(out_path, index=False)
    print(f"wrote {out_path} rows={len(out)} donors={out['donor_id'].nunique()}")


if __name__ == "__main__":
    main()
