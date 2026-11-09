#!/usr/bin/env python3
"""Build Rev2 concordance and interneuron figures from committed result tables."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
FIGURES = REPO / "report" / "plos-one" / "figures-v2"


def concordance_figure():
    old = pd.read_csv(RESULTS / "reviewer-rev2-concordance.csv")
    rows = []
    for scope, label in (("within_tissue", "Mouse, within tissue"),
                         ("cross_tissue", "Mouse, cross tissue")):
        row = old[(old["species"] == "Mouse") & (old["scope"] == scope)
                  & (old["gene_set"] == "primary_detection_filtered")].iloc[0]
        rows.append({"species_scope": label, "observed": row["observed_mean_pearson"],
                     "null_mean": row["joint_null_mean"],
                     "null_ci95_low": row["joint_null_ci95_low"],
                     "null_ci95_high": row["joint_null_ci95_high"],
                     "p_value": row["joint_p_greater"], "status": "redone_joint"})
    row = old[(old["species"] == "Naked mole-rat") & (old["scope"] == "within_tissue")].iloc[0]
    rows.append({"species_scope": "Naked mole rat, within tissue",
                 "observed": row["observed_mean_pearson"], "null_mean": row["joint_null_mean"],
                 "null_ci95_low": row["joint_null_ci95_low"],
                 "null_ci95_high": row["joint_null_ci95_high"],
                 "p_value": row["joint_p_greater"], "status": "already_joint"})
    for slug, label in (("human", "Human, within cohort"), ("fly", "Fruit fly, within tissue"),
                        ("worm", "C. elegans, within tissue"),
                        ("zebrafish", "Zebrafish, within tissue")):
        row = pd.read_csv(RESULTS / f"reviewer-rev2-joint-concordance-{slug}.csv").iloc[0]
        rows.append({"species_scope": label, "observed": row["observed_mean_pearson"],
                     "null_mean": row["joint_null_mean"],
                     "null_ci95_low": row["joint_null_ci95_low"],
                     "null_ci95_high": row["joint_null_ci95_high"],
                     "p_value": row["joint_p_greater"], "status": row["status"]})
    data = pd.DataFrame(rows)
    data.to_csv(RESULTS / "reviewer-rev2-concordance-figure-data.csv", index=False)

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(9, 5.6))
    y = np.arange(len(data))[::-1]
    ax.scatter(data["observed"], y, color="#147d82", s=58, zorder=3,
               label="Observed mean pairwise Pearson correlation")
    ax.scatter(data["null_mean"], y, color="#737373", marker="s", s=38, zorder=3,
               label="Joint donor-permutation null mean")
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
    genes = pd.read_csv(RESULTS / "reviewer-rev2-interneuron-genes.csv")
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.8), sharex=False, sharey=False)
    for ax, short in zip(axes, ("VIP", "pvalb", "sst")):
        part = genes[genes["short_type"] == short]
        both = part["tier"] == "both_FDR"
        ax.scatter(part.loc[~both, "slope_mssm"], part.loc[~both, "slope_hbcc"],
                   s=10, color="#a0a0a0", alpha=0.65, linewidths=0)
        ax.scatter(part.loc[both, "slope_mssm"], part.loc[both, "slope_hbcc"],
                   s=22, facecolors="none", edgecolors="#147d82", linewidths=0.9)
        ax.axhline(0, color="#dddddd", lw=0.7)
        ax.axvline(0, color="#dddddd", lw=0.7)
        lim = max(np.abs(part[["slope_mssm", "slope_hbcc"]].to_numpy()).max(), 1e-6) * 1.05
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_title(f"{short} interneurons")
        ax.set_xlabel("MSSM age slope")
        ax.set_ylabel("HBCC age slope")
    fig.suptitle("Age-effect slopes for MSSM FDR genes; rings mark same-direction HBCC FDR replication")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig4-interneuron.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    FIGURES.mkdir(parents=True, exist_ok=True)
    concordance_figure()
    interneuron_figure()
