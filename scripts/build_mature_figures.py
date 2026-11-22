#!/usr/bin/env python3
"""Build mature-to-old concordance (Fig 3) and interneuron (Fig 4) figures
from the committed mature result tables. Same visual style as the Rev2
builders: teal circles = observed, gray squares + bars = null mean + 95% CI."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
FIGURES = REPO / "report" / "plos-one" / "figures-v2"


def concordance_figure():
    mouse = pd.read_csv(RESULTS / "followup-mature-concordance-mouse.csv")
    human = pd.read_csv(RESULTS / "followup-mature-concordance-human.csv").iloc[0]
    rows = []
    for scope, label in (("within_tissue", "Mouse, within tissue"),
                         ("cross_tissue", "Mouse, cross tissue")):
        row = mouse[(mouse["scope"] == scope)
                    & (mouse["gene_set"] == "primary_detection_filtered")].iloc[0]
        rows.append({"species_scope": label, "observed": row["observed_mean_pearson"],
                     "null_mean": row["joint_null_mean"],
                     "null_ci95_low": row["joint_null_ci95_low"],
                     "null_ci95_high": row["joint_null_ci95_high"],
                     "p_value": row["joint_p_greater"]})
    rows.append({"species_scope": "Human, within cohort",
                 "observed": human["observed_mean_pearson"],
                 "null_mean": human["joint_null_mean"],
                 "null_ci95_low": human["joint_null_ci95_low"],
                 "null_ci95_high": human["joint_null_ci95_high"],
                 "p_value": human["joint_p_greater"]})
    data = pd.DataFrame(rows)
    data.to_csv(RESULTS / "followup-mature-concordance-figure-data.csv", index=False)

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(9, 5.6))
    y = np.arange(len(data))[::-1]
    ax.scatter(data["observed"], y, color="#147d82", s=58, zorder=3, label="Observed")
    ax.scatter(data["null_mean"], y, color="#737373", marker="s", s=38, zorder=3,
               label="Joint null mean")
    for yi, row in zip(y, data.itertuples()):
        ax.plot([row.null_ci95_low, row.null_ci95_high], [yi, yi], color="#737373", lw=1.4)
    ax.set_yticks(y, data["species_scope"])
    ax.set_xlabel("Mean pairwise correlation of gene-level age effects")
    ax.axvline(0, color="#cccccc", lw=0.8, zorder=0)
    ax.legend(loc="lower right", frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig3-concordance.png", dpi=300)
    fig.savefig(FIGURES / "fig3-concordance.pdf")
    plt.close(fig)


def interneuron_figure():
    genes = pd.read_csv(RESULTS / "followup-mature-interneuron-genes.csv")
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(8.5, 7.2))
    for ax, short in zip(axes.ravel(), ("VIP", "pvalb", "sst", "GABA")):
        part = genes[genes["short_type"] == short]
        both = part["tier"] == "both_FDR"
        base_color = "#a0a0a0" if short == "GABA" else "#4478a6"
        ax.scatter(part.loc[~both, "slope_mssm"], part.loc[~both, "slope_hbcc"],
                   s=10, color=base_color, alpha=0.55, linewidths=0)
        ax.scatter(part.loc[both, "slope_mssm"], part.loc[both, "slope_hbcc"],
                   s=22, facecolors="none", edgecolors="#147d82", linewidths=0.9)
        ax.axhline(0, color="#dddddd", lw=0.7)
        ax.axvline(0, color="#dddddd", lw=0.7)
        lim = max(np.abs(part[["slope_mssm", "slope_hbcc"]].to_numpy()).max(), 1e-6) * 1.05
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.locator_params(axis="both", nbins=4)
        ax.set_title("Broad GABAergic neurons" if short == "GABA" else f"{short} interneurons")
        ax.set_xlabel("MSSM age slope")
        ax.set_ylabel("HBCC age slope")
    fig.suptitle("Mature MSSM FDR genes; rings mark same-direction HBCC FDR replication")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig4-interneuron.png", dpi=300)
    fig.savefig(FIGURES / "fig4-interneuron.pdf")
    plt.close(fig)


if __name__ == "__main__":
    FIGURES.mkdir(parents=True, exist_ok=True)
    concordance_figure()
    interneuron_figure()
    print("fig3 and fig4 rebuilt from mature tables")
