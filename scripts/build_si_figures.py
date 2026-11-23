#!/usr/bin/env python3
"""Build Agescape Supporting Information figures S1 to S6.

The transport-across-cohorts figure is NOT built: its expected numbers do
not verify against the sources (see verify_dropped_transport(), which
prints the exact mismatches and skips the figure). The six built figures
are numbered S1 to S6 with no gap.

Every figure reads its numbers from the committed source files, asserts
them against the expected values, and only then plots. Style matches the
main figures: Liberation Serif, slate blue / amber palette, spines
top/right off, explicit ticks, 300 dpi PNG plus PDF.

No-jargon rule: no bare codes (P1m, P0m, E1, E2, A1, A3, Block 5,
Block 10, MDE80, 624 family, 431 family, DIES, SURVIVES, both-pass,
empirical-only, composite-only, fail-both) and no em/en dashes anywhere
in figure text. Matplotlib is set to ASCII minus signs.
"""
import csv
import math
import os
from collections import Counter

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = os.path.expanduser("~/workspace/agescape-repo")
OUT = os.path.join(REPO, "report", "plos-one", "figures-si")
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Liberation Serif", "DejaVu Serif"],
        "font.size": 11,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "axes.linewidth": 1.3,
        "axes.unicode_minus": False,
        "svg.fonttype": "none",
    }
)

MOUSE = "#2e4a62"   # mouse dark slate blue
HUMAN = "#9a6b1e"   # human deep amber/brown
STEEL = "#6a8caf"   # muted steel blue
GREY = "#d9d9d9"    # light grey
GRID = "#cccccc"    # grid
GREEN = "#2e7d32"
RED = "#b03030"


def style_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(length=4, width=1.2)
    ax.grid(True, axis="x", color=GRID, linewidth=0.7, alpha=0.7)
    ax.set_axisbelow(True)


def save(fig, stem):
    fig.savefig(os.path.join(OUT, stem + ".png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(OUT, stem + ".pdf"), bbox_inches="tight")
    plt.close(fig)
    print("wrote", stem + ".png/.pdf")


def rows(path):
    with open(os.path.join(REPO, path), newline="") as fh:
        return list(csv.DictReader(fh))


# ---------------------------------------------------------------- S1
def build_s1():
    data = rows("results/followup-mature-p1m-genes.csv")
    assert len(data) == 1914, len(data)
    n_units = [int(r["n_units_significant"]) for r in data]
    dist = Counter(n_units)
    assert dist[1] == 1212, dist[1]
    assert sum(1 for n in n_units if n >= 2) == 702
    assert sum(1 for n in n_units if n >= 3) == 361
    assert sum(1 for n in n_units if n >= 5) == 132
    assert max(n_units) == 11
    both = sum(1 for r in data if int(r["n_units_up"]) > 0 and int(r["n_units_down"]) > 0)
    assert both == 57, both

    xs = list(range(1, 12))
    ys = [dist[x] for x in xs]
    fig, ax = plt.subplots(figsize=(9, 5.6))
    style_ax(ax)
    ax.bar(xs, ys, color=MOUSE, width=0.72, zorder=3)
    for x, y in zip(xs, ys):
        ax.text(x, y + 28, f"{y:,}", ha="center", va="bottom", fontsize=10)
    ax.set_title(
        "Genes in the cleaned mature-to-old signature (1,914 genes),\n"
        "by number of supporting test groups",
        fontsize=12,
    )
    ax.set_xlabel(
        "Supporting test groups per gene\n"
        "(one test group is one tissue-by-species group)"
    )
    ax.set_ylabel("Genes")
    ax.set_xticks(xs)
    ax.set_yticks([0, 400, 800, 1200])
    ax.set_ylim(0, 1360)
    ax.set_xlim(0.3, 11.7)
    ax.text(
        0.015,
        -0.30,
        "702 genes have support in 2 or more groups; 361 in 3 or more; 132 in 5 or more.\n"
        "57 genes have support in both directions (up in some groups, down in others).",
        transform=ax.transAxes,
        fontsize=10,
        va="top",
        ha="left",
    )
    save(fig, "s1-p1m-multi-unit-support")


# ---------------------------------------------------------------- S2
def build_s2():
    mouse = rows("results/followup-mature-composition-mouse.csv")
    human = rows("results/followup-mature-composition-human.csv")
    assert len(mouse) == 114, len(mouse)
    assert len(human) == 285, len(human)
    all_rows = [(r, "mouse") for r in mouse] + [(r, "human") for r in human]
    assert len(all_rows) == 399

    def find(cell, cohort, organ):
        for r, sp in all_rows:
            if r["cell_type"] == cell and r["cohort"] == cohort and r["organ"] == organ:
                return float(r["effect"]), sp
        raise AssertionError(f"not found: {cell} {cohort} {organ}")

    # (short key label, cell_type, cohort, organ, expected slope)
    targets = [
        ("Human lung club cell", "club cell", "Lung_Integrated", "lung", -0.080407),
        ("Human lung multiciliated epithelial cell", "multiciliated epithelial cell",
         "Lung_Integrated", "lung", -0.047002),
        ("Mount Sinai brain samples (MSSM) L6b glutamatergic neuron",
         "L6b glutamatergic neuron of the primary motor cortex", "MSSM_DLPFC", "brain", 0.019765),
        ("Mount Sinai brain samples (MSSM) oligodendrocyte", "oligodendrocyte",
         "MSSM_DLPFC", "brain", 0.017269),
        ("Mount Sinai brain samples (MSSM) oligodendrocyte precursor cell",
         "oligodendrocyte precursor cell", "MSSM_DLPFC", "brain", -0.011896),
        ("Mouse limb-muscle macrophage", "macrophage", "TMS", "Limb_Muscle", -0.142945),
        ("Mouse limb-muscle mesenchymal stem cell", "mesenchymal stem cell",
         "TMS", "Limb_Muscle", -0.100457),
        ("Human lung conventional dendritic cell", "conventional dendritic cell",
         "Lung_Integrated", "lung", 0.022897),
        ("Human lung multiciliated columnar cell of tracheobronchial tree",
         "multiciliated columnar cell of tracheobronchial tree", "Lung_Integrated",
         "lung", -0.038079),
        ("Seattle Alzheimer brain atlas (SEA-AD) L5 extratelencephalic projecting "
         "glutamatergic cortical neuron",
         "L5 extratelencephalic projecting glutamatergic cortical neuron",
         "SEAAD_DLPFC", "brain", -0.021299),
    ]
    # Same-scan rows live in the secondary-checks file, not the main files.
    sec = rows("results/followup-mature-composition-secondary-checks.csv")
    assert len(sec) == 5, len(sec)
    assert all(r["status"] == "SCAN_FINDING" for r in sec)
    sec_lookup = {(r["cell_type"], r["cohort"]): float(r["effect"]) for r in sec}

    found = {}
    for label, cell, cohort, organ, expected in targets:
        try:
            val, _ = find(cell, cohort, organ)
        except AssertionError:
            val = sec_lookup[(cell, cohort)]
        assert abs(val - expected) < 5e-7, (label, val, expected)
        found[label] = val
    # Order the numbered key by slope, most negative first.
    ordered = sorted(found.items(), key=lambda kv: kv[1])

    slopes = sorted(float(r["effect"]) for r, _ in all_rows)
    n = len(slopes)
    highlight = {round(v, 6): i for i, (lab, v) in enumerate(ordered, start=1)}

    fig = plt.figure(figsize=(11, 8.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.1, 1.35], hspace=0.55)
    ax = fig.add_subplot(gs[0])
    style_ax(ax)
    ax.scatter(slopes, range(1, n + 1), s=14, color="#8a8a8a", zorder=2,
               label="Other 389 rows")
    for lab, val in ordered:
        rank = slopes.index(val) + 1
        num = highlight[round(val, 6)]
        col = MOUSE if "Mouse" in lab or "mouse" in lab.lower()[:6] else HUMAN
        ax.scatter([val], [rank], s=90, color=col, edgecolor="white",
                   linewidth=1.0, zorder=4)
        ax.text(val + 0.0035, rank, str(num), fontsize=10, fontweight="bold",
                color=col, va="center", ha="left", zorder=5)
    ax.axvline(0, color="black", linewidth=1.0, zorder=1)
    ax.set_title("All 399 cell-type fraction slopes, with the 10 that clear the bar highlighted",
                 fontsize=12)
    ax.set_xlabel("CLR slope (centered log ratio slope)")
    ax.set_ylabel("Rows sorted by slope (1 = most negative)")
    ax.set_xticks([-0.15, -0.10, -0.05, 0, 0.05])
    ax.set_xlim(-0.165, 0.075)
    ax.set_yticks([1, 100, 200, 300, 399])
    ax.set_ylim(0, 415)
    ax.legend(loc="lower right", frameon=False)
    ax.text(0.015, -0.24,
            "Mouse slopes are per month and human slopes are per year.\n"
            "MSSM: Mount Sinai brain samples. SEA-AD: Seattle Alzheimer brain atlas.",
            transform=ax.transAxes, fontsize=10, va="top", ha="left")

    axk = fig.add_subplot(gs[1])
    axk.axis("off")
    axk.set_xlim(0, 11)
    axk.set_ylim(0, 10)
    axk.text(0, 9.6, "The 10 highlighted rows:", fontsize=11, fontweight="bold",
             va="top", ha="left")
    # Short key labels; MSSM and SEA-AD are expanded in the note above the key.
    short = {
        "Mouse limb-muscle macrophage": "Mouse limb-muscle macrophage",
        "Mouse limb-muscle mesenchymal stem cell": "Mouse limb-muscle mesenchymal stem cell",
        "Human lung club cell": "Human lung club cell",
        "Human lung multiciliated epithelial cell": "Human lung multiciliated epithelial cell",
        "Human lung multiciliated columnar cell of tracheobronchial tree":
            "Human lung multiciliated columnar cell of tracheobronchial tree",
        "Seattle Alzheimer brain atlas (SEA-AD) L5 extratelencephalic projecting "
        "glutamatergic cortical neuron":
            "SEA-AD L5 extratelencephalic projecting glutamatergic cortical neuron",
        "Mount Sinai brain samples (MSSM) oligodendrocyte precursor cell":
            "MSSM oligodendrocyte precursor cell",
        "Mount Sinai brain samples (MSSM) oligodendrocyte": "MSSM oligodendrocyte",
        "Mount Sinai brain samples (MSSM) L6b glutamatergic neuron":
            "MSSM L6b glutamatergic neuron",
        "Human lung conventional dendritic cell": "Human lung conventional dendritic cell",
    }
    import textwrap
    left = ordered[:5]
    right = ordered[5:]
    for j, col_items in enumerate([left, right]):
        x0 = 0.1 + j * 5.6
        y = 8.2
        for lab, val in col_items:
            num = highlight[round(val, 6)]
            sign = "+" if val > 0 else ""
            text = f"{num}. {short[lab]}, {sign}{val:.3f}"
            lines = textwrap.wrap(text, width=66)
            for k, ln in enumerate(lines):
                axk.text(x0 + (0.35 if k > 0 else 0), y, ln,
                         fontsize=9.5, va="top", ha="left")
                y -= 0.8
            y -= 0.75
    save(fig, "s2-fraction-scan")


# ---------------------------------------------------------------- S3 (DROP)
def verify_dropped_transport():
    """Expected S3 numbers do not verify. Print mismatches and build nothing."""
    print("Transport figure DROPPED:")
    print("  1. PBMC: task expects 0 of 35; "
          "reports/followup-mature-organ-transports-batch-c.md says "
          "0 of 36 under the cleaned signature (23 Zhu + 13 Thomson).")
    print("  2. Retinal pigment epithelium: task expects no passing units; "
          "batch C report says RPE bulk 1 of 1 passes both tests and SURVIVED "
          "under both signatures (cleaned effect -0.0211 per decade, q 0.0443, "
          "random-signature q 0.0005).")
    print("  3. GTEx: task expects the whole-blood positive bulk result in "
          "Track 1; the GTEx report puts the Whole Blood both-pass result in "
          "Track 2 (effect +0.0119, q 1.289e-22), while Track 1 both-pass "
          "tissues are 5, all negative direction.")
    print("  Verified parts (not plotted): core 48 units, 2 pass both tests, "
          "both negative (PsychAD microglia, Mount Sinai microglia); "
          "retina 0 of 144; liver 0 of 8.")


# ---------------------------------------------------------------- S4
def build_s3():
    data = rows("results/followup-mature-between-subject-yazar-models.csv")
    sel = {}
    for r in data:
        if (r["is_stratified_global"] == "True" and r["depth"] == "full_1500"
                and r["model_spec"] == "stratified_global_primary"):
            sel[r["part"]] = r
    a1, a3 = sel["A1"], sel["A3"]
    assert a1["n_donors"] == "925" and a3["n_donors"] == "925"
    assert abs(float(a1["slope_sd_per_decade"]) - (-0.026112)) < 5e-7
    assert abs(float(a1["pvalue"]) - 3.55e-06) < 1e-7
    assert abs(float(a3["slope_sd_per_decade"]) - (-0.077554)) < 5e-7
    assert abs(float(a3["pvalue"]) - 1.97e-13) < 1e-14

    pts = [
        ("Distance between donors\n"
         "(Hellinger distance: a distance\nbetween two expression profiles)",
         float(a1["slope_sd_per_decade"]), float(a1["ci95_low_sd"]),
         float(a1["ci95_high_sd"]), "p 3.6e-06"),
        ("Spread of expression\nwithin a donor",
         float(a3["slope_sd_per_decade"]), float(a3["ci95_low_sd"]),
         float(a3["ci95_high_sd"]), "p 2.0e-13"),
    ]
    fig, ax = plt.subplots(figsize=(9, 4.6))
    style_ax(ax)
    ys = [1, 0]
    for (lab, slope, lo, hi, ptxt), y in zip(pts, ys):
        ax.errorbar([slope], [y], xerr=[[slope - lo], [hi - slope]],
                    fmt="o", color=MOUSE, markersize=9, elinewidth=2.5,
                    capsize=6, capthick=2.5, zorder=3)
        ax.text(slope - 0.004, y + 0.30,
                f"slope {slope:.3f}, {ptxt}", fontsize=10, ha="right",
                va="bottom", color=MOUSE)
    ax.axvline(0, color="black", linewidth=1.0)
    ax.set_yticks(ys)
    ax.set_yticklabels([p[0] for p in pts])
    ax.set_xticks([-0.10, -0.08, -0.06, -0.04, -0.02, 0])
    ax.set_xlim(-0.115, 0.015)
    ax.set_ylim(-0.6, 1.6)
    ax.set_xlabel("Change per decade of age (standard deviations)")
    ax.set_title("Donor-to-donor distances shrink with age in blood",
                 fontsize=12)
    ax.text(0.0, -0.34,
            "Yazar blood cohort, 925 donors, full-depth data. Slopes are change per decade.\n"
            "Gao heart cohort: no clear change (slopes -0.008, p 0.77 and -0.062, p 0.085).",
            transform=ax.transAxes, fontsize=10, va="top", ha="left")
    save(fig, "s3-between-subject-shrinkage")


# ---------------------------------------------------------------- S5
def build_s4():
    # Verdicts from the report table; counts cross-checked against the CSV.
    data = rows("results/followup-mature-block10-benchmark-results.csv")
    cnt = Counter()
    for r in data:
        if r["both_pass"] == "True":
            cnt[(r["signature_id"], r["arm"])] += 1
    assert cnt[("Harries2011_6gene", "unadjusted")] == 7
    assert cnt[("Harries2011_6gene", "batch_adjusted")] == 5
    assert cnt[("Peters2015_1497gene", "unadjusted")] == 16
    assert cnt[("Peters2015_1497gene", "batch_adjusted")] == 16
    assert cnt[("Sood2015_150probe", "unadjusted")] == 4
    assert cnt[("Sood2015_150probe", "batch_adjusted")] == 0
    assert cnt[("deMagalhaes2009_73gene", "unadjusted")] == 0
    assert cnt[("deMagalhaes2009_73gene", "batch_adjusted")] == 0
    assert len(data) == 246

    sigs = [
        ("Harries 2011, 6 genes", GREEN,
         "Holds up: passes both tests in 5 blood cell groups after batch adjustment (7 before)."),
        ("Peters 2015, 1,497 genes", GREEN,
         "Holds up: passes both tests in 16 blood cell groups, with and without batch adjustment."),
        ("Sood 2015, 150 probes", RED,
         "Lost with batch adjustment: 4 groups pass before adjustment, 0 after."),
        ("de Magalhaes 2009, 73 genes", "#8a8a8a",
         "No verdict: 0 groups pass both tests in either cohort, before or after adjustment."),
    ]
    fig, ax = plt.subplots(figsize=(10, 5.4))
    ax.axis("off")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.text(0, 9.7, "Four published aging signatures, tested in Gao heart and Yazar blood",
            fontsize=12, fontweight="bold", va="top", ha="left")
    ax.text(0, 8.9,
            "Passes both tests means the age signal beats both the composite test and the "
            "random-signature test.",
            fontsize=10, va="top", ha="left")
    for i, (name, col, verdict) in enumerate(sigs):
        y = 7.4 - i * 1.7
        ax.scatter([0.25], [y], s=170, color=col, zorder=3)
        ax.text(0.55, y + 0.28, name, fontsize=11, fontweight="bold",
                va="center", ha="left")
        ax.text(0.55, y - 0.32, verdict, fontsize=10, va="center", ha="left")
        ax.plot([0.1, 9.9], [y - 0.85, y - 0.85], color=GRID, linewidth=0.8)
    ax.text(0, 0.15,
            "Two cohorts: Gao heart (13 groups) and Yazar blood (17 groups). "
            "Four signatures cannot prove a claim about all published signatures.",
            fontsize=10, va="bottom", ha="left")
    save(fig, "s4-published-signatures")


# ---------------------------------------------------------------- S6
def bh_q(pvals):
    """Benjamini-Hochberg q values, monotone, for a list of p values."""
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    q = [0.0] * m
    running = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        running = min(running, pvals[i] * m / rank)
        q[i] = running
    return q


def build_s5():
    data = rows("results/followup-mature-concordance-figure-data.csv")
    assert len(data) == 3, len(data)
    d = {r["species_scope"]: r for r in data}
    mw = d["Mouse, within tissue"]
    mc = d["Mouse, cross tissue"]
    hu = d["Human, within cohort"]
    assert abs(float(mw["p_value"]) - 0.0219780219780219) < 1e-9
    assert abs(float(hu["p_value"]) - 0.8151848151848152) < 1e-9
    qs = bh_q([float(mw["p_value"]), float(mc["p_value"]), float(hu["p_value"])])
    assert abs(qs[0] - 0.065934) < 1e-4, qs[0]
    assert round(qs[0], 3) == 0.066

    groups = [
        ("Mouse\nwithin tissue", float(mw["observed"]), float(mw["null_mean"]),
         float(mw["null_ci95_low"]), float(mw["null_ci95_high"]), "p 0.022", MOUSE),
        ("Mouse\ncross tissue", float(mc["observed"]), float(mc["null_mean"]),
         float(mc["null_ci95_low"]), float(mc["null_ci95_high"]), "p 0.103", MOUSE),
        ("Human\nwithin cohort", float(hu["observed"]), float(hu["null_mean"]),
         float(hu["null_ci95_low"]), float(hu["null_ci95_high"]), "p 0.815", HUMAN),
    ]
    fig, ax = plt.subplots(figsize=(9, 5.8))
    style_ax(ax)
    x = np.arange(3)
    w = 0.36
    ax.bar(x - w / 2, [g[1] for g in groups], w, color=[g[6] for g in groups],
           label="Observed", zorder=3)
    ax.bar(x + w / 2, [g[2] for g in groups], w, color="#bdbdbd",
           label="Random-gene null mean", zorder=3)
    for xi, g in zip(x, groups):
        ax.errorbar([xi + w / 2], [g[2]], yerr=[[g[2] - g[3]], [g[4] - g[2]]],
                    fmt="none", ecolor="black", elinewidth=2.0, capsize=6,
                    capthick=2.0, zorder=4)
        ax.text(xi, max(g[1], g[4]) + 0.018, g[5], ha="center", va="bottom",
                fontsize=10)
    ax.set_xticks(x)
    ax.set_xticklabels([g[0] for g in groups])
    ax.set_yticks([0, 0.1, 0.2, 0.3])
    ax.set_ylim(0, 0.36)
    ax.set_ylabel("Mean correlation of age effects (Pearson)")
    ax.set_title("Age-signal agreement across cell types, observed versus random-gene null",
                 fontsize=12)
    ax.legend(loc="upper right", frameon=False)
    ax.text(0.0, -0.30,
            "Mouse within-tissue: p 0.022, q 0.066 (false discovery rate correction over the 3 rows "
            "shown).\nNominal only, not significant after correction. Human: p 0.815, null.\n"
            "Mouse: 59 groups, 14 mice. Human: 127 groups, 416 donors, 7 cohorts.",
            transform=ax.transAxes, fontsize=10, va="top", ha="left")
    save(fig, "s5-concordance")


# ---------------------------------------------------------------- S7
def build_s6():
    def get(path, gene_set, role):
        for r in rows(path):
            if (r["is_stratified_global"] == "True" and r["metric"] == "mean_ectopic_count"
                    and r["gene_set"] == gene_set and r["analysis_role"] == role):
                return r
        raise AssertionError(f"not found: {path} {gene_set} {role}")

    e1g = get("results/followup-mature-fulldepth-gao-ectopic.csv", "E1", "primary_full_depth")
    e1y = get("results/followup-mature-fulldepth-yazar-ectopic.csv", "E1", "primary_full_depth")
    e2g = get("results/followup-mature-fulldepth-gao-ectopic.csv", "E2_stratum",
              "primary_depth_matched_e2")
    e2y = get("results/followup-mature-fulldepth-yazar-ectopic.csv", "E2_stratum",
              "primary_depth_matched_e2")
    for r, slope, p in [(e1g, -0.063, 0.468), (e1y, -0.019, 0.073),
                        (e2g, -0.127, 0.031), (e2y, -0.033, 0.0031)]:
        assert abs(float(r["slope_sd_per_decade"]) - slope) < 5e-4, (r, slope)
        assert abs(float(r["pvalue"]) - p) < 5e-4, (r, p)

    pts = [
        ("Gao heart cohort\nprimary off-program count",
         e1g, "slope -0.063, p 0.47", MOUSE),
        ("Gao heart cohort\nstricter off-program count",
         e2g, "slope -0.127, p 0.031", MOUSE),
        ("Yazar blood cohort\nprimary off-program count",
         e1y, "slope -0.019, p 0.073", HUMAN),
        ("Yazar blood cohort\nstricter off-program count",
         e2y, "slope -0.033, p 0.0031", HUMAN),
    ]
    fig, ax = plt.subplots(figsize=(9, 5.6))
    style_ax(ax)
    ys = [3, 2, 1, 0]
    for (lab, r, txt, col), y in zip(pts, ys):
        slope = float(r["slope_sd_per_decade"])
        lo = float(r["ci95_low_sd"])
        hi = float(r["ci95_high_sd"])
        ax.errorbar([slope], [y], xerr=[[slope - lo], [hi - slope]],
                    fmt="o", color=col, markersize=9, elinewidth=2.5,
                    capsize=6, capthick=2.5, zorder=3)
        ax.text(slope - 0.006, y + 0.30, txt, fontsize=10, ha="right",
                va="bottom", color=col)
    ax.axvline(0, color="black", linewidth=1.0)
    ax.set_yticks(ys)
    ax.set_yticklabels([p[0] for p in pts])
    ax.set_xticks([-0.25, -0.20, -0.15, -0.10, -0.05, 0, 0.05, 0.10])
    ax.set_xlim(-0.28, 0.14)
    ax.set_ylim(-0.6, 3.6)
    ax.set_xlabel("Change per decade of age (standard deviations)")
    ax.set_title("Stricter off-program count falls with age; primary count shows no clear change",
                 fontsize=12)
    ax.text(0.0, -0.30,
            "Off-program: genes a cell expresses that belong to other cell types.\n"
            "Primary count: Gao 104 donors, Yazar 925 donors. Stricter count: Gao 103 donors, "
            "Yazar 919 donors.",
            transform=ax.transAxes, fontsize=10, va="top", ha="left")
    save(fig, "s6-ectopic-e1-e2")


if __name__ == "__main__":
    build_s1()
    build_s2()
    verify_dropped_transport()
    build_s3()
    build_s4()
    build_s5()
    build_s6()
    print("done")
