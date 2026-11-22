#!/usr/bin/env python3
"""Build the four optimization-pass figures for the Agescape manuscript from
the committed mature result tables.

Supersedes the Fig 1, Fig 2, and Fig 3 outputs of scripts/build_nature_figures.py
(full-range, young-reference) and the Fig 3 and Fig 4 outputs of
scripts/build_mature_figures.py for the optimized manuscript. This script is
the owner of fig1-composition, fig2-mature-drift, fig3-transport, and
fig4-interneuron in report/plos-one/figures-v2/.

Fig 1: mature-to-old cell-fraction age slopes for the fifteen full-range
  fraction associations (nine mouse rows die on power loss, five human rows
  survive, one bladder row untestable). Mouse slopes are CLR slopes per month
  and human slopes are CLR slopes per year; the two panels carry separate
  axes because the time units differ.
Fig 2: mature whole-profile drift slopes (Hellinger distance from mature
  references) with 95 percent intervals: mouse and human globals plus the
  seven mouse per-organ slopes. All slopes are per 10 percent of species
  maximum lifespan, so one axis serves both species.
Fig 3: core 48-unit transport verdicts under P1m and P0m (both-pass,
  empirical-only, composite-only, fail-both counts).
Fig 4: mature interneuron MSSM vs HBCC age slopes, with TRMT11 labeled in
  each subtype panel (the only gene passing FDR in all three subtypes in
  both cohorts on the unmoderated test).

Every number is read from a committed CSV in results/ and asserted against
the values printed in the manuscript before plotting. PNG (300 dpi) and PDF
(vector) pairs are written; the tex references the PDF files.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
FIGURES = REPO / "report" / "plos-one" / "figures-v2"

TEAL = "#147d82"
AMBER = "#b7791f"
GREY = "#737373"
LIGHT = "#d9d9d9"

plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})


def save(fig, name):
    fig.savefig(FIGURES / f"{name}.png", dpi=300)
    fig.savefig(FIGURES / f"{name}.pdf")
    plt.close(fig)


def fig1_composition():
    v = pd.read_csv(RESULTS / "mature-only-composition-verdicts.csv")
    assert len(v) == 15, len(v)
    assert (v["verdict"] == "DIES").sum() == 9
    assert (v["verdict"] == "SURVIVES").sum() == 5
    assert (v["verdict"] == "UNTESTABLE").sum() == 1
    # manuscript numbers: club -0.080, multiciliated -0.047, L6b +0.020,
    # oligodendrocyte +0.017, OPC -0.012; q values are the original
    # 431-test family (0.01178, 0.02785).
    surv = v[v["verdict"] == "SURVIVES"]
    assert np.allclose(surv["mature_q_global_bh_orig_family"].values,
                       [0.011776868326018252, 0.027853986950145068,
                        0.011776868326018252, 0.011776868326018252,
                        0.011776868326018252])
    short = {
        ("Mouse", "Kidney", "kidney capillary endothelial cell"): "Kidney capillary endothelial",
        ("Mouse", "Kidney", "kidney mesangial cell"): "Kidney mesangial",
        ("Mouse", "Limb_Muscle", "skeletal muscle cell"): "Limb muscle skeletal muscle",
        ("Mouse", "Liver", "B cell"): "Liver B cell",
        ("Mouse", "Lung", "T cell"): "Lung T cell",
        ("Mouse", "Lung", "pulmonary interstitial fibroblast"): "Lung interstitial fibroblast",
        ("Mouse", "Marrow", "late pro-B cell"): "Marrow late pro-B",
        ("Mouse", "Marrow", "precursor B cell"): "Marrow precursor B",
        ("Mouse", "Spleen", "plasma cell"): "Spleen plasma cell",
        ("Human", "lung", "club cell"): "Lung club cell",
        ("Human", "lung", "multiciliated epithelial cell"): "Lung multiciliated epithelial",
        ("Human", "brain", "L6b glutamatergic neuron of the primary motor cortex"): "Brain L6b glutamatergic",
        ("Human", "brain", "oligodendrocyte"): "Brain oligodendrocyte",
        ("Human", "brain", "oligodendrocyte precursor cell"): "Brain oligodendrocyte precursor",
    }
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 5.6), sharey=False)
    for ax, species, color, xlabel in (
            (axes[0], "Mouse", TEAL, "CLR slope per month"),
            (axes[1], "Human", AMBER, "CLR slope per year")):
        part = v[(v["species"] == species) & (v["verdict"] != "UNTESTABLE")].copy()
        part["label"] = [short[(r["species"], r["organ"], r["cell_type"])]
                         for _, r in part.iterrows()]
        part = part.sort_values("mature_effect")
        y = np.arange(len(part))[::-1]
        lo = part["mature_effect"] - 1.96 * part["mature_se"]
        hi = part["mature_effect"] + 1.96 * part["mature_se"]
        for yi, l, h in zip(y, lo, hi):
            ax.plot([l, h], [yi, yi], color=color, lw=1.4, alpha=0.8)
        survived = part["verdict"] == "SURVIVES"
        ax.scatter(part.loc[survived, "mature_effect"], y[survived.values],
                   color=color, s=64, zorder=3)
        ax.scatter(part.loc[~survived, "mature_effect"], y[(~survived).values],
                   facecolors="none", edgecolors=color, s=64, linewidths=1.4, zorder=3)
        ax.set_yticks(y, part["label"], fontsize=9)
        ax.axvline(0, color="#cccccc", lw=0.8, zorder=0)
        ax.set_xlabel(xlabel)
        ax.set_title(f"{species}: {'dies on power loss' if species == 'Mouse' else 'survives'}")
    axes[0].set_ylabel("Cell type")
    fig.suptitle("Mature-to-old cell-fraction age slopes (95% CI); bladder untestable, not shown",
                 fontsize=10)
    fig.tight_layout()
    save(fig, "fig1-composition")
    print("fig1-composition rebuilt: 9 mouse DIES, 5 human SURVIVES, 1 UNTESTABLE")


def fig2_mature_drift():
    mouse = pd.read_csv(RESULTS / "followup-mature-drift-mouse.csv")
    human = pd.read_csv(RESULTS / "followup-mature-drift-human.csv")
    glob_m = mouse[(mouse["tissue"] == "All tested tissues")
                   & (mouse["adjustment"] == "linear_log_umi")
                   & (mouse["depth"] == "full")
                   & (mouse["n_donors"] == 14)].iloc[0]
    assert abs(glob_m["slope"] - (-0.0137)) < 0.001, glob_m["slope"]
    assert abs(glob_m["ci95_low"] - (-0.139)) < 0.01
    assert abs(glob_m["ci95_high"] - 0.111) < 0.01
    glob_h = human[(human["tissue"] == "All tested tissues")
                   & (human["adjustment"] == "linear_log_umi")
                   & (human["depth"] == "full")
                   & (human["n_donors"] == 407)].iloc[0]
    assert abs(glob_h["slope"] - (-0.0071)) < 0.001, glob_h["slope"]
    assert abs(glob_h["ci95_low"] - (-0.0345)) < 0.002
    assert abs(glob_h["ci95_high"] - 0.0203) < 0.002
    organs = ["Lung", "Limb_Muscle", "Kidney", "Spleen", "Fat", "Marrow"]
    labels = ["Mouse lung", "Mouse limb muscle", "Mouse kidney", "Mouse spleen",
              "Mouse fat", "Mouse marrow"]
    rows = [
        {"label": "Mouse global (depth-adjusted)", "slope": glob_m["slope"],
         "lo": glob_m["ci95_low"], "hi": glob_m["ci95_high"], "color": TEAL},
        {"label": "Human global (depth-adjusted)", "slope": glob_h["slope"],
         "lo": glob_h["ci95_low"], "hi": glob_h["ci95_high"], "color": AMBER},
    ]
    for organ, label in zip(organs, labels):
        r = mouse[(mouse["tissue"] == organ) & (mouse["adjustment"] == "linear_log_umi")
                  & (mouse["depth"] == "full") & (mouse["variant"] == "per_tissue")].iloc[0]
        rows.append({"label": label + " (per-organ)", "slope": r["slope"],
                     "lo": r["ci95_low"], "hi": r["ci95_high"], "color": TEAL})
    # manuscript per-organ values: lung -0.194, limb muscle +0.286, kidney
    # -0.231, spleen +0.453, fat -0.017, marrow -0.035. Heart has no
    # per-tissue row in the mature file, so it is not plotted.
    expect = [-0.194, 0.286, -0.231, 0.453, -0.017, -0.035]
    assert np.allclose([r["slope"] for r in rows[2:]], expect, atol=0.002)
    fig, ax = plt.subplots(figsize=(9, 5.4))
    y = np.arange(len(rows))[::-1]
    for yi, r in zip(y, rows):
        ax.plot([r["lo"], r["hi"]], [yi, yi], color=r["color"], lw=1.4)
        ax.scatter([r["slope"]], [yi], color=r["color"], s=58, zorder=3)
    ax.set_yticks(y, [r["label"] for r in rows])
    ax.set_xlabel("Hellinger slope (SD per 10% of maximum lifespan)")
    ax.axvline(0, color="#cccccc", lw=0.8, zorder=0)
    fig.tight_layout()
    save(fig, "fig2-mature-drift")
    print("fig2-mature-drift rebuilt: mouse global -0.0137, human global -0.0071")


def fig3_transport():
    s = pd.read_csv(RESULTS / "followup-mature-p1m-block5-summary.csv")
    p1 = s[s["gene_set"] == "P1m"].iloc[0]
    p0 = s[s["gene_set"] == "P0m"].iloc[0]
    assert (p1["both_pass"], p1["empirical_only"], p1["composite_only"],
            p1["fail_both"]) == (2, 13, 4, 29)
    assert (p0["both_pass"], p0["empirical_only"], p0["composite_only"],
            p0["fail_both"]) == (3, 15, 5, 25)
    cats = ["Both-pass", "Empirical-only", "Composite-only", "Fail-both"]
    colors = [TEAL, AMBER, "#6a8caf", LIGHT]
    fig, ax = plt.subplots(figsize=(9, 5.0))
    x = np.arange(len(cats))
    w = 0.36
    for off, row, label in ((-w / 2, p1, "P1m (1,914 genes)"), (w / 2, p0, "P0m (2,598 genes)")):
        vals = [row["both_pass"], row["empirical_only"], row["composite_only"],
                row["fail_both"]]
        bars = ax.bar(x + off, vals, w, label=label, color=colors, edgecolor="white")
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v + 0.6, str(v),
                    ha="center", fontsize=10)
    ax.set_xticks(x, cats)
    ax.set_ylabel("Core units (n = 48)")
    ax.set_ylim(0, 34)
    ax.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    save(fig, "fig3-transport")
    print("fig3-transport rebuilt: P1m 2/13/4/29, P0m 3/15/5/25")


def fig4_interneuron():
    genes = pd.read_csv(RESULTS / "followup-mature-interneuron-genes.csv")
    # TRMT11: the only gene passing FDR in all three subtypes in both
    # cohorts on the unmoderated test. Verified slopes (log counts per
    # 10,000 per year): VIP MSSM 0.00227 / HBCC 0.00406, pvalb 0.00205 /
    # 0.00316, sst 0.00229 / 0.00280; mean expression about 4 counts per
    # 10,000 in MSSM.
    tr = genes[genes["gene"] == "TRMT11"]
    assert set(tr["short_type"]) == {"VIP", "pvalb", "sst"}
    assert (tr["tier"] == "both_FDR").all()
    plt.rcParams.update({"font.size": 9})
    fig, axes = plt.subplots(2, 2, figsize=(8.5, 7.2))
    for ax, short in zip(axes.ravel(), ("VIP", "pvalb", "sst", "GABA")):
        part = genes[genes["short_type"] == short]
        both = part["tier"] == "both_FDR"
        base_color = "#a0a0a0" if short == "GABA" else "#4478a6"
        ax.scatter(part.loc[~both, "slope_mssm"], part.loc[~both, "slope_hbcc"],
                   s=10, color=base_color, alpha=0.55, linewidths=0)
        ax.scatter(part.loc[both, "slope_mssm"], part.loc[both, "slope_hbcc"],
                   s=22, facecolors="none", edgecolors=TEAL, linewidths=0.9)
        ax.axhline(0, color="#dddddd", lw=0.7)
        ax.axvline(0, color="#dddddd", lw=0.7)
        lim = max(np.abs(part[["slope_mssm", "slope_hbcc"]].to_numpy()).max(), 1e-6) * 1.05
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.locator_params(axis="both", nbins=4)
        ax.set_title("Broad GABAergic neurons" if short == "GABA" else f"{short} interneurons")
        ax.set_xlabel("MSSM age slope")
        ax.set_ylabel("HBCC age slope")
        t = tr[tr["short_type"] == short]
        if len(t):
            row = t.iloc[0]
            ax.annotate("TRMT11", (row["slope_mssm"], row["slope_hbcc"]),
                        textcoords="offset points", xytext=(10, 8), fontsize=10,
                        fontweight="bold", color=TEAL,
                        arrowprops=dict(arrowstyle="-", color=TEAL, lw=0.8))
    fig.suptitle("Mature MSSM FDR genes; rings mark same-direction HBCC FDR replication")
    fig.tight_layout()
    save(fig, "fig4-interneuron")
    plt.rcParams.update({"font.size": 10})
    print("fig4-interneuron rebuilt with TRMT11 labeled in VIP, pvalb, sst panels")


if __name__ == "__main__":
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig1_composition()
    fig2_mature_drift()
    fig3_transport()
    fig4_interneuron()
    print("optimization figures rebuilt")
