#!/usr/bin/env python3

from pathlib import Path
import sys

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_reviewer_phase1_core import cluster_nmr_effect
from src.reviewer_stats import random_effects_meta

RESULTS = Path(__file__).resolve().parent.parent / "results"


def main():
    core = pd.read_csv(RESULTS / "h3-cross-species-vector-global.csv")
    fly = pd.read_csv(RESULTS / "extended-afca-vector-global.csv")
    fish = pd.read_csv(RESULTS / "extended-zebrafish-vector-fine-global.csv")
    rows = []
    for reference in (1, 6):
        worm = pd.read_csv(
            RESULTS / ("extended-worm-vector-global-platform-audit.csv" if reference == 1
                       else "extended-worm-vector-global-platform-audit-reference6.csv")
        )
        mouse = core[(core["species"] == "Mouse") & (core["tissue"] == "All tested tissues")
                     & (core["metric"] == "profile_hellinger") & (core["model"] == "unadjusted")
                     & (core["donor_set"] == "all_donor_ids")].iloc[0]
        human = core[(core["species"] == "Human") & (core["tissue"] == "All tested tissues")
                     & (core["metric"] == "profile_hellinger") & (core["model"] == "unadjusted")].iloc[0]
        nmr_effect, nmr_se, _ = cluster_nmr_effect()
        fly_row = fly[(fly["tissue"] == "All tested tissues") & (fly["model"] == "depth_adjusted")].iloc[0]
        worm_row = worm[worm["analysis"] == "all_annotated_cells"].iloc[0]
        fish_row = fish[(fish["tissue"] == "All tested tissues") & (fish["model"] == "depth_adjusted")].iloc[0]
        inputs = [
            ("Mouse", mouse["effect_sd"], mouse["se"], "results/h3-cross-species-vector-global.csv:114"),
            ("Human", human["effect_sd"], human["se"], "results/h3-cross-species-vector-global.csv:110"),
            ("Naked mole-rat", nmr_effect, nmr_se, "results/h3-cross-species-vector-donors.csv; clustered refit"),
            ("Fruit fly", fly_row["effect_sd"], fly_row["se"], "results/extended-afca-vector-global.csv:7"),
            ("C. elegans", worm_row["effect_sd_per_14_days"], worm_row["se"],
             f"results/extended-worm-vector-global-platform-audit{'-reference6' if reference == 6 else ''}.csv:2"),
            ("Zebrafish", fish_row["effect_sd"], fish_row["se"],
             "results/extended-zebrafish-vector-fine-global.csv; depth-adjusted all-tissues row"),
        ]
        for species, effect, se, source in inputs:
            rows.append({"worm_reference_day": reference, "scope": "species_input", "species": species,
                         "effect": effect, "se": se, "source": source})
        meta = random_effects_meta([row[1] for row in inputs], [row[2] for row in inputs])
        rows.append({"worm_reference_day": reference, "scope": "random_effects_meta", "species": "All six",
                     **meta, "source": "six species rows in this table"})
    pd.DataFrame(rows).to_csv(RESULTS / "reviewer-phase1-species-random-effects.csv", index=False)


if __name__ == "__main__":
    main()
