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


def save(fig, name):
    fig.savefig(FIGURES / f"reviewer-phase3-{name}.png")
    fig.savefig(FIGURES / f"reviewer-phase3-{name}.pdf")
    plt.close(fig)


def mouse_drift():
    data = pd.read_csv(RESULTS / "reviewer-phase1-mouse-drift-unit-effects.csv")
    data["unit"] = data["tissue"] + " | " + data["cell_type"]
    data = data.sort_values("slope_per_year")
    data.to_csv(RESULTS / "reviewer-phase3-mouse-drift-figure-data.csv", index=False)
    fig, ax = plt.subplots(figsize=(10, 14))
    y = np.arange(len(data))
    ax.errorbar(data["slope_per_year"], y,
                xerr=[data["slope_per_year"] - data["ci95_low"], data["ci95_high"] - data["slope_per_year"]],
                fmt="o", ms=3, color="#0072B2", ecolor="#B0B0B0", elinewidth=0.6)
    ax.axvline(0, color="#555555", ls="--", lw=0.8)
    ax.set_yticks(y, data["unit"], fontsize=6)
    ax.set_xlabel("Age slope in Hellinger distance per year, with 95% CI")
    ax.set_title("Mouse expression-vector displacement by tissue and cell type")
    fig.tight_layout()
    save(fig, "mouse-vector-drift")


def interneurons():
    data = pd.read_csv(RESULTS / "reviewer-phase1-interneuron-dotplot-data.csv")
    data["short_type"] = data["cell_type"].map({
        "VIP GABAergic cortical interneuron": "VIP",
        "pvalb GABAergic cortical interneuron": "pvalb",
        "sst GABAergic cortical interneuron": "sst",
    })
    data = data.dropna(subset=["short_type"])
    gene_order = data.groupby("gene")["q_bh_mssm"].min().sort_values().index
    data["gene"] = pd.Categorical(data["gene"], categories=gene_order, ordered=True)
    data.to_csv(RESULTS / "reviewer-phase3-interneuron-figure-data.csv", index=False)
    fig, ax = plt.subplots(figsize=(8, 14))
    scatter = ax.scatter(data["short_type"], data["gene"], c=data["slope_log_cp10k_per_year_mssm"],
                         s=np.maximum(12, -np.log10(data["q_bh_mssm"]) * 20), cmap="RdBu_r",
                         vmin=-max(abs(data["slope_log_cp10k_per_year_mssm"])),
                         vmax=max(abs(data["slope_log_cp10k_per_year_mssm"])))
    ax.set_ylabel("Gene")
    ax.set_xlabel("Interneuron subtype")
    ax.set_title("MSSM age-associated interneuron gene changes")
    fig.colorbar(scatter, ax=ax, label="MSSM slope, log CP10K per year")
    ax.tick_params(axis="y", labelsize=6)
    fig.tight_layout()
    save(fig, "interneuron-gene-changes")


def composition():
    data = pd.read_csv(RESULTS / "reviewer-phase1-composition-total-adjusted.csv")
    data = data[(data["pseudocount"] == 0.5) & (data["q_bh"] < 0.05)
                & (data["species"].isin(["Mouse", "Human"]))].copy()
    data["unit"] = data["species"] + " | " + data["organ"] + " | " + data["cohort"] + " | " + data["cell_type"]
    data = data.sort_values("effect")
    data.to_csv(RESULTS / "reviewer-phase3-composition-figure-data.csv", index=False)
    fig, ax = plt.subplots(figsize=(11, 9))
    y = np.arange(len(data))
    colors = data["species"].map({"Mouse": "#0072B2", "Human": "#D55E00"})
    ax.errorbar(data["effect"], y,
                xerr=[data["effect"] - data["ci95_low"], data["ci95_high"] - data["effect"]],
                fmt="none", ecolor="#B0B0B0", elinewidth=0.8)
    ax.scatter(data["effect"], y, c=colors, s=30)
    ax.axvline(0, color="#555555", ls="--", lw=0.8)
    ax.set_yticks(y, data["unit"], fontsize=6)
    ax.set_xlabel("Total-cell-adjusted CLR age slope, with 95% CI")
    ax.set_title("FDR-significant mouse and human cell-fraction changes")
    fig.tight_layout()
    save(fig, "composition-changes")


def universality():
    summary = pd.read_csv(RESULTS / "reviewer-phase1-universality-summary.csv")
    data = summary[summary["metric"] == "mean_pearson"].copy()
    order = ["Mouse", "Human", "Naked mole-rat", "Fruit fly", "C. elegans", "Zebrafish"]
    data["species"] = pd.Categorical(data["species"], order, ordered=True)
    data = data.sort_values("species")
    data.to_csv(RESULTS / "reviewer-phase3-universality-figure-data.csv", index=False)
    fig, ax = plt.subplots(figsize=(9, 5))
    y = np.arange(len(data))
    ax.errorbar(data["observed"], y, xerr=[np.zeros(len(data)), np.zeros(len(data))], fmt="o",
                color="#0072B2", label="Observed mean pairwise concordance")
    ax.errorbar(data["null_mean"], y + 0.15,
                xerr=[data["null_mean"] - data["null_ci95_low"], data["null_ci95_high"] - data["null_mean"]],
                fmt="s", color="#D55E00", ecolor="#D55E00", label="Age-permutation null mean and 95% interval")
    ax.axvline(0, color="#555555", ls="--", lw=0.8)
    ax.set_yticks(y, data["species"])
    ax.set_xlabel("Mean pairwise gene-effect Pearson correlation")
    ax.set_title("Gene-effect concordance across cell types")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    save(fig, "universality-concordance")


def main():
    set_pub_style()
    mouse_drift()
    interneurons()
    composition()
    universality()


if __name__ == "__main__":
    main()
