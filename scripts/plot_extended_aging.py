#!/usr/bin/env python3

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
FIGURES = REPO / "figures"


def interval(axis, label, mean, low, high, position, color, text):
    axis.errorbar(mean, position, xerr=[[mean - low], [high - mean]], fmt="o", color=color,
                  capsize=4, markersize=7)
    axis.text(mean, position + 0.22, text, ha="center", va="bottom", fontsize=9)
    return label


def main():
    FIGURES.mkdir(exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    blue = "#40566B"
    red = "#C24A36"
    gray = "#8B8F94"

    ax = axes[0, 0]
    fly = pd.read_csv(RESULTS / "extended-afca-vector-global.csv").query(
        "model == 'depth_adjusted' and tissue == 'All tested tissues'"
    ).iloc[0]
    worm1 = pd.read_csv(RESULTS / "extended-worm-vector-global-platform-audit.csv").query(
        "analysis == 'all_annotated_cells'"
    ).iloc[0]
    worm6 = pd.read_csv(RESULTS / "extended-worm-vector-global-platform-audit-reference6.csv").query(
        "analysis == 'all_annotated_cells'"
    ).iloc[0]
    fish = pd.read_csv(RESULTS / "extended-zebrafish-vector-broad-global.csv").query(
        "model == 'depth_adjusted' and tissue == 'Testis broad compartments'"
    ).iloc[0]
    entries = [
        ("Fly, 30 to 70 days", fly.effect_sd, fly.ci95_low, fly.ci95_high, blue, "depth adjusted"),
        ("Worm, day-1 reference", worm1.effect_sd_per_14_days,
         worm1.effect_sd_per_14_days - 1.96 * worm1.se,
         worm1.effect_sd_per_14_days + 1.96 * worm1.se, red,
         f"platform permutation p = {worm1.platform_stratified_permutation_p:.4f}"),
        ("Worm, day-6 reference", worm6.effect_sd_per_14_days,
         worm6.effect_sd_per_14_days - 1.96 * worm6.se,
         worm6.effect_sd_per_14_days + 1.96 * worm6.se, red,
         f"platform permutation p = {worm6.platform_stratified_permutation_p:.4f}"),
        ("Zebrafish, 5 to 27 months", fish.effect_sd, fish.ci95_low, fish.ci95_high, blue,
         "broad compartments"),
    ]
    for i, (label, mean, low, high, color, text) in enumerate(entries[::-1]):
        interval(ax, label, mean, low, high, i, color, text)
    ax.axvline(0, color="black", linewidth=1, linestyle="--")
    ax.set_yticks(range(len(entries)), [entry[0] for entry in entries[::-1]])
    ax.set_ylim(-0.6, 3.65)
    ax.set_xlabel("Standardized age effect over the dataset age window")
    ax.set_title("A. Whole-vector effects depend on species and reference")
    ax.grid(axis="x", alpha=0.25)

    ax = axes[0, 1]
    fly_c = pd.read_csv(RESULTS / "extended-afca-composition-tests.csv").query("sensitivity == 'all_samples'")
    worm_c = pd.read_csv(RESULTS / "extended-worm-composition-tests.csv").query(
        "sensitivity == 'platform_adjusted' and analysis == 'all_annotated_cells'"
    )
    fish_c = pd.read_csv(RESULTS / "extended-zebrafish-composition-tests.csv")
    counts = [
        fly_c.q_bh.lt(0.05).sum(), len(fly_c),
        worm_c.q_bh.lt(0.05).sum(), len(worm_c),
        fish_c.query("analysis == 'broad_all_ages' and sensitivity == 'all_samples'").q_bh.lt(0.05).sum(), 3,
        fish_c.query("analysis == 'fine_5_to_22_months' and sensitivity == 'all_samples'").q_bh.lt(0.05).sum(), 15,
    ]
    labels = ["Fly", "Worm", "Zebrafish broad", "Zebrafish fine"]
    significant = counts[::2]
    totals = counts[1::2]
    ax.bar(labels, significant, color=[blue, red, blue, blue], width=0.65)
    for i, (sig, total) in enumerate(zip(significant, totals)):
        ax.text(i, sig + 0.7, f"{sig}/{total}", ha="center", fontsize=11)
    ax.set_ylabel("Cell types with FDR-significant fraction change")
    ax.set_title("B. Cell-fraction changes occur in all three new datasets")
    ax.tick_params(axis="x", rotation=15)
    ax.set_ylim(0, max(significant) * 1.18)
    ax.grid(axis="y", alpha=0.25)

    ax = axes[1, 0]
    parametric = pd.read_csv(RESULTS / "extended-zebrafish-vector-fine-tests.csv").query(
        "model == 'depth_adjusted'"
    )[["cell_type", "q_bh"]]
    permutation = pd.read_csv(RESULTS / "extended-zebrafish-vector-permutation.csv")[
        ["cell_type", "permutation_q_bh"]
    ]
    merged = parametric.merge(permutation, on="cell_type")
    ax.scatter(merged.q_bh, merged.permutation_q_bh, s=55, color=blue)
    limits = [0.008, 1.1]
    ax.plot(limits, limits, color="black", linestyle="--", linewidth=1)
    ax.axhline(0.05, color=red, linewidth=1)
    ax.axvline(0.05, color=red, linewidth=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_xlabel("Parametric depth-adjusted q")
    ax.set_ylabel("Exact age-permutation q")
    ax.set_title("C. Zebrafish fine-grained hits do not survive exact FDR")
    ax.grid(alpha=0.25)

    ax = axes[1, 1]
    bayes_entries = []
    for label, prefix, model in (
        ("Fly", "afca", "depth_adjusted_fixed_group"),
        ("Worm, day-1 reference", "worm", "depth_adjusted_platform_nextseq"),
        ("Zebrafish broad", "zebrafish", "depth_adjusted"),
    ):
        data = pd.read_csv(RESULTS / f"extended-{prefix}-bayesian-posterior.csv").query(
            "scope == 'species' and model == @model and prior_sd == 0.5"
        ).iloc[0]
        bayes_entries.append((label, data.posterior_mean, data.hdi94_low, data.hdi94_high))
    for i, (label, mean, low, high) in enumerate(bayes_entries[::-1]):
        ax.errorbar(mean, i, xerr=[[mean - low], [high - mean]], fmt="o", color=red if "Worm" in label else blue,
                    capsize=4, markersize=7)
    ax.axvline(0, color="black", linewidth=1, linestyle="--")
    ax.set_yticks(range(len(bayes_entries)), [entry[0] for entry in bayes_entries[::-1]])
    ax.set_xlabel("Posterior mean and 94% HDI")
    ax.set_title("D. Bayesian global effects remain uncertain")
    ax.grid(axis="x", alpha=0.25)

    fig.tight_layout()
    fig.savefig(FIGURES / "extended-cross-species-summary.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
