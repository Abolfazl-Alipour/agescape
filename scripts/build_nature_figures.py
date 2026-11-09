#!/usr/bin/env python3
"""Rebuild the four Rev2 manuscript figures in journal style.

Nature-style presentation: Liberation Serif throughout (the metric-identical
open clone of Times New Roman; Times itself is not installed on this machine),
built at print size so type is 9-11 pt on the page, thick black axes, large
outward tick marks, no in-figure titles, subtitles, side panels, footnotes
or pale gridlines. The figure captions in report/plos-one/final-plos-one.tex
carry the descriptions; nothing in the figures says anything the captions
and Results text do not already say.

This script is the single owner of figures-v2/fig1-fig4. It supersedes the
figure-building parts of scripts/build_rev2_manuscript_figures.py (fig3, fig4)
and ~/workspace/agescape-final/figures-v2/build_v2.py (fig1, fig2). Data are
NOT changed: same CSVs, same rows, same values, same colours.

Outputs: report/plos-one/figures-v2/fig1-composition.png/.pdf,
fig2-mouse-drift.png/.pdf, fig3-concordance.png/.pdf, fig4-interneuron.png/.pdf.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
FIGURES = REPO / "report" / "plos-one" / "figures-v2"

TEAL = "#147d82"
AMBER = "#c07f1f"
GREY = "#737373"
LIGHT_GREY = "#999999"

# Times New Roman is not installed; Liberation Serif is its metric-identical
# open clone (same glyph widths, same look on the page).
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Liberation Serif", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 11,
    "axes.titlesize": 11,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "axes.linewidth": 1.1,
    "axes.edgecolor": "black",
    "axes.labelcolor": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "text.color": "black",
    "xtick.major.size": 6,
    "ytick.major.size": 6,
    "xtick.major.width": 1.1,
    "ytick.major.width": 1.1,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "savefig.facecolor": "white",
})


def save(fig, name):
    fig.savefig(FIGURES / f"{name}.png", dpi=300)
    fig.savefig(FIGURES / f"{name}.pdf")
    plt.close(fig)
    print("saved", name)


def clean(ax):
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(length=6, width=1.1)


# ---------------------------------------------------------------- fig 1
def fig1():
    df = pd.read_csv(RESULTS / "reviewer-phase3-composition-figure-data.csv")
    df = df.sort_values("effect").reset_index(drop=True)
    short = {"hepatocyte": "Hepatocyte", "pulmonary interstitial fibroblast": "Interstitial fibroblast",
             "kidney capillary endothelial cell": "Capillary endothelial", "kidney mesangial cell": "Mesangial cell",
             "endothelial cell of hepatic sinusoid": "Sinusoid endothelial", "club cell": "Club cell",
             "bladder cell": "Bladder cell", "precursor B cell": "Precursor B cell",
             "hepatic stellate cell": "Hepatic stellate", "multiciliated epithelial cell": "Multiciliated epithelial",
             "smooth muscle cell": "Smooth muscle", "oligodendrocyte precursor cell": "Oligodendrocyte precursor",
             "sst GABAergic cortical interneuron": "SST interneuron", "GABAergic neuron": "GABAergic neuron",
             "L6b glutamatergic neuron of the primary motor cortex": "L6b glutamatergic neuron",
             "oligodendrocyte": "Oligodendrocyte", "pulmonary alveolar type 2 cell": "Alveolar type 2",
             "non-classical monocyte": "Non-classical monocyte", "T cell": "T cell", "NK cell": "NK cell",
             "plasma cell": "Plasma cell", "myeloid leukocyte": "Myeloid leukocyte", "B cell": "B cell",
             "epithelial cell of proximal tubule": "Proximal tubule epithelial", "late pro-B cell": "Late pro-B cell"}
    fig, ax = plt.subplots(figsize=(7.0, 6.6))
    fig.subplots_adjust(left=0.40, right=0.97, top=0.94, bottom=0.16)
    y = np.arange(len(df))
    for i, r in df.iterrows():
        c = TEAL if r.species == "Mouse" else AMBER
        ax.plot([r.ci95_low, r.ci95_high], [i, i], color=c, lw=2.0, alpha=0.6,
                solid_capstyle="round", zorder=1)
        ax.plot([r.effect], [i], "o", ms=7, color=c, mec="white", mew=1.2, zorder=3)
    ax.axvline(0, color="black", lw=1.1)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{short.get(t, t)} ({o})" for t, o in zip(df.cell_type, df.organ)],
                       fontsize=9.5)
    ax.set_xlabel("CLR age slope (share of one cell type, per year)")
    ax.set_xlim(-0.32, 0.19)
    ax.invert_yaxis()
    clean(ax)
    ax.text(0.02, 1.02, "Falling share", transform=ax.transAxes, fontsize=10)
    ax.text(0.98, 1.02, "Rising share", transform=ax.transAxes, fontsize=10, ha="right")
    ax.legend(handles=[Line2D([], [], marker="o", color="none", markerfacecolor=TEAL,
                               markersize=8, label="Mouse"),
                       Line2D([], [], marker="o", color="none", markerfacecolor=AMBER,
                               markersize=8, label="Human")],
              loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, frameon=False)
    save(fig, "fig1-composition")


# ---------------------------------------------------------------- fig 2
def fig2():
    df = pd.read_csv(RESULTS / "reviewer-phase3-mouse-drift-figure-data.csv")
    df = df.sort_values("slope_per_year").reset_index(drop=True)
    sig = df[df.q_bh < 0.05]
    short = {"bladder cell": "bladder", "mesenchymal stem cell": "mesenchymal stem",
             "podocyte": "podocyte", "fibroblast of lung": "fibroblast",
             "fenestrated cell": "fenestrated"}
    fig, ax = plt.subplots(figsize=(7.0, 6.6))
    fig.subplots_adjust(left=0.11, right=0.97, top=0.96, bottom=0.16)
    y = np.arange(len(df))
    for i, r in df.iterrows():
        s = r.q_bh < 0.05
        c = TEAL if s else GREY
        ax.plot([r.ci95_low, r.ci95_high], [i, i], color=c, lw=1.6 if s else 1.0,
                alpha=0.9 if s else 0.5, solid_capstyle="round", zorder=1)
        ax.plot([r.slope_per_year], [i], "o", ms=7 if s else 4.5, color=c,
                mec="white", mew=1.1, zorder=3)
    for i, r in sig.iterrows():
        if r.cell_type == "podocyte":
            # Podocyte sits directly below mesenchymal stem; its label would
            # collide, so it gets a short leader line to clear space.
            ax.plot([r.ci95_high, 0.073], [i, i], color=LIGHT_GREY, lw=0.8, zorder=1)
            ax.text(0.075, i, short[r.cell_type], va="center", fontsize=9.5)
        else:
            ax.text(r.ci95_high + 0.004, i, short[r.cell_type], va="center",
                    fontsize=9.5)
    ax.axvline(0, color="black", lw=1.1)
    ax.set_yticks([])
    ax.set_ylabel("62 units, ranked by slope")
    ax.set_xlabel("Hellinger slope per year (distance from young reference)")
    ax.set_xlim(-0.045, 0.135)
    clean(ax)
    ax.legend(handles=[Line2D([], [], marker="o", color="none", markerfacecolor=TEAL,
                               markersize=8, label="Passing correction (5 units)"),
                       Line2D([], [], marker="o", color="none", markerfacecolor=GREY,
                               markersize=6, label="Other units")],
              loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, frameon=False)
    save(fig, "fig2-mouse-drift")


# ---------------------------------------------------------------- fig 3
def fig3():
    data = pd.read_csv(RESULTS / "reviewer-rev2-concordance-figure-data.csv")
    order = ["Mouse, within tissue", "Mouse, cross tissue", "Naked mole rat, within tissue",
             "Human, within cohort", "Fruit fly, within tissue", "C. elegans, within tissue",
             "Zebrafish, within tissue"]
    data["species_scope"] = pd.Categorical(data["species_scope"], order, ordered=True)
    data = data.sort_values("species_scope").reset_index(drop=True)

    def ptext(p):
        # Three decimals throughout, matching the Results text exactly
        # (mouse 0.044, cross-tissue 0.153, human 0.700, zebrafish 0.361).
        return f"p = {p:.3f}"

    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    fig.subplots_adjust(left=0.36, right=0.97, top=0.96, bottom=0.20)
    y = np.arange(len(data))[::-1]
    ax.scatter(data["observed"], y, color=TEAL, s=81, zorder=3, label="Observed")
    ax.scatter(data["null_mean"], y, color=GREY, marker="s", s=49, zorder=3,
               label="Joint null mean")
    for yi, row in zip(y, data.itertuples()):
        ax.plot([row.null_ci95_low, row.null_ci95_high], [yi, yi], color=GREY, lw=1.6)
    ax.set_yticks(y, [f"{row.species_scope}, {ptext(row.p_value)}"
                      for row in data.itertuples()], fontsize=9.5)
    ax.set_xlabel("Mean pairwise correlation of gene-level age effects")
    ax.axvline(0, color=LIGHT_GREY, lw=0.9, zorder=0)
    clean(ax)
    ax.legend(loc="lower right", frameon=False)
    save(fig, "fig3-concordance")


# ---------------------------------------------------------------- fig 4
def fig4():
    genes = pd.read_csv(RESULTS / "reviewer-rev2-interneuron-genes-moderated.csv")
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 6.2))
    fig.subplots_adjust(left=0.11, right=0.97, top=0.92, bottom=0.10,
                        wspace=0.28, hspace=0.55)
    panels = (("VIP", "VIP interneurons"), ("pvalb", "pvalb interneurons"),
              ("sst", "sst interneurons"), ("GABA", "Broad GABAergic neurons"))
    for ax, (short, title), letter in zip(axes.ravel(), panels, "abcd"):
        part = genes[genes["short_type"] == short]
        both = part["tier"] == "both_FDR"
        base_color = "#a0a0a0" if short == "GABA" else "#4478a6"
        ax.scatter(part.loc[~both, "slope_mssm"], part.loc[~both, "slope_hbcc"],
                   s=12, color=base_color, alpha=0.6, linewidths=0)
        ax.scatter(part.loc[both, "slope_mssm"], part.loc[both, "slope_hbcc"],
                   s=30, facecolors="none", edgecolors=TEAL, linewidths=1.0)
        ax.axhline(0, color=LIGHT_GREY, lw=0.8)
        ax.axvline(0, color=LIGHT_GREY, lw=0.8)
        lim = max(np.abs(part[["slope_mssm", "slope_hbcc"]].to_numpy()).max(), 1e-6) * 1.05
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.locator_params(axis="both", nbins=4)
        ax.set_title(title, fontsize=11, pad=6)
        ax.set_xlabel("MSSM age slope", fontsize=10)
        ax.set_ylabel("HBCC age slope", fontsize=10)
        ax.tick_params(labelsize=9)
        clean(ax)
        ax.text(-0.10, 1.06, f"({letter})", transform=ax.transAxes,
                fontsize=12, fontweight="bold", va="bottom", ha="right")
    save(fig, "fig4-interneuron")


if __name__ == "__main__":
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig1()
    fig2()
    fig3()
    fig4()
    print("all done")
