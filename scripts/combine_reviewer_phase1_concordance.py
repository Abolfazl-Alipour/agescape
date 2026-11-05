#!/usr/bin/env python3

from pathlib import Path
import sys

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_reviewer_phase1_gene_concordance import summarize

RESULTS = REPO / "results"


def main():
    for stem in ("pairs", "permutations", "units"):
        frames = []
        for suffix in ("", "-human", "-worm", "-zebrafish", "-nmr"):
            path = RESULTS / f"reviewer-phase1-universality-{stem}{suffix}.csv"
            if path.exists():
                frame = pd.read_csv(path)
                frame["source_suffix"] = suffix
                frames.append(frame)
        pd.concat(frames, ignore_index=True).to_csv(
            RESULTS / f"reviewer-phase1-universality-{stem}.csv", index=False
        )
    permutations = pd.read_csv(RESULTS / "reviewer-phase1-universality-permutations.csv")
    summarize(permutations).to_csv(RESULTS / "reviewer-phase1-universality-summary.csv", index=False)


if __name__ == "__main__":
    main()
