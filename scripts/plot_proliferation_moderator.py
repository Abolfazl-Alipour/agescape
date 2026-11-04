#!/usr/bin/env python3

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from src.plotting import set_pub_style

RESULTS = REPO / "results"
FIGURES = REPO / "figures"
SPECIES = ["Mouse", "Human", "Fruit fly", "C. elegans", "Zebrafish"]


def panel(ax, data, title, xlabel):
    data = data.set_index("scope").reindex(SPECIES).reset_index()
    y = np.arange(len(data))[::-1]
    for yi, row in zip(y, data.itertuples()):
        ax.plot([row.ci95_low, row.ci95_high], [yi, yi], color="#555555", lw=1.2)
        ax.scatter(row.post_minus_proliferative, yi, s=38, color="#0072B2", zorder=3)
    ax.axvline(0, color="#999999", lw=0.8, ls="--")
    ax.set_yticks(y, data["scope"])
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.spines[["top", "right"]].set_visible(False)


def main():
    set_pub_style()
    tests = pd.read_csv(RESULTS / "proliferation-moderator-tests.csv")
    primary = tests[
        (((tests["analysis"] == "expression_vector") & (tests["sensitivity"] == "depth_adjusted"))
         | ((tests["analysis"] == "cell_fraction") & (tests["sensitivity"] == "all_donors")))
        & tests["scope"].isin(SPECIES)
    ]
    fig, axes = plt.subplots(2, 2, figsize=(13, 8))
    for row, (analysis, label) in enumerate([
        ("expression_vector", "Expression-vector displacement"),
        ("cell_fraction", "Cell-type fraction change"),
    ]):
        for column, (response, endpoint) in enumerate([
            ("effect_sd", "Signed effect"), ("absolute_effect_sd", "Absolute effect"),
        ]):
            selected = primary[(primary["analysis"] == analysis) & (primary["response"] == response)]
            panel(
                axes[row, column], selected, f"{label}: {endpoint}",
                "Post-mitotic minus proliferative-capable, SD per 10% lifespan",
            )
    fig.suptitle("Proliferation class does not show a consistent aging-effect difference")
    fig.tight_layout()
    FIGURES.mkdir(exist_ok=True)
    fig.savefig(FIGURES / "proliferation-moderator.png")
    fig.savefig(FIGURES / "proliferation-moderator.pdf")


if __name__ == "__main__":
    main()
