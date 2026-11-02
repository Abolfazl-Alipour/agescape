#!/usr/bin/env python3

import argparse
from pathlib import Path
import sys

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_cross_species_bayesian_vector import fit_model

RESULTS = REPO / "results"
INPUTS = {
    "Fruit fly": "extended-afca-vector-donors.csv",
    "C. elegans": "extended-worm-vector-donors.csv",
    "Zebrafish": "extended-zebrafish-vector-donors.csv",
}
PREFIXES = {"Fruit fly": "afca", "C. elegans": "worm", "Zebrafish": "zebrafish"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("species", choices=INPUTS)
    parser.add_argument("--skip-diagnostics", action="store_true")
    args = parser.parse_args()
    donors = pd.read_csv(RESULTS / INPUTS[args.species])
    rows = []
    priors = []
    diagnostics = []
    for depth_adjusted, prior_sd in ((False, 0.5), (True, 0.25), (True, 0.5), (True, 1.0)):
        print(f"{args.species}: depth_adjusted={depth_adjusted}, prior={prior_sd}", flush=True)
        posterior, prior, diagnostic = fit_model(
            donors, args.species, depth_adjusted, prior_sd,
            diagnostics=not args.skip_diagnostics, draws=2000,
            fixed_group_effects=args.species == "Fruit fly",
            covariate_column="platform_nextseq" if args.species == "C. elegans" else "sex_male",
        )
        rows.extend(posterior)
        priors.append(prior)
        diagnostics.append(diagnostic)
    prefix = PREFIXES[args.species]
    pd.DataFrame(rows).to_csv(RESULTS / f"extended-{prefix}-bayesian-posterior.csv", index=False)
    pd.DataFrame(priors).to_csv(RESULTS / f"extended-{prefix}-bayesian-prior-predictive.csv", index=False)
    pd.DataFrame(diagnostics).to_csv(RESULTS / f"extended-{prefix}-bayesian-diagnostics.csv", index=False)


if __name__ == "__main__":
    main()
