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
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = str(Path(__file__).resolve().parent.parent)
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
        "(one test group is one tissue-by-cell-type group)"
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
    pooled = rows("results/mature-all-species-composition.csv")
    assert len(mouse) == 114, len(mouse)
    assert len(human) == 285, len(human)
    assert len(pooled) == 577, len(pooled)
    assert sum(1 for r in pooled if float(r["q_pooled_mature_bh"]) < 0.05) == 39

    pooled_by_key = {
        (r["species"], r["organ"], r["cohort"], r["cell_type"]): r
        for r in pooled
    }
    assert len(pooled_by_key) == 577, len(pooled_by_key)
    for r in mouse + human:
        key = (r["species"], r["organ"], r["cohort"], r["cell_type"])
        assert key in pooled_by_key, key
        assert abs(float(r["effect"]) - float(pooled_by_key[key]["effect"])) < 1e-12, key

    def find(df, cell, cohort, organ):
        hits = [r for r in df if r["cell_type"] == cell and r["cohort"] == cohort
                and r["organ"] == organ]
        assert len(hits) == 1, (cell, cohort, organ, len(hits))
        return hits[0]

    # (label, cell_type, cohort, organ, expected slope, expected pooled q, expected rank, species)
    targets = [
        ("Mouse limb-muscle macrophage", "macrophage", "TMS", "Limb_Muscle",
         -0.142945, 0.01243, 8, "mouse"),
        ("Mouse limb-muscle mesenchymal stem cell", "mesenchymal stem cell",
         "TMS", "Limb_Muscle", -0.100457, 0.00812, 15, "mouse"),
        ("Human lung club cell", "club cell", "Lung_Integrated", "lung",
         -0.080407, 0.00424, 6, "human"),
        ("Human lung multiciliated epithelial cell", "multiciliated epithelial cell",
         "Lung_Integrated", "lung", -0.047002, 0.01243, 21, "human"),
        ("Human lung multiciliated columnar cell of tracheobronchial tree",
         "multiciliated columnar cell of tracheobronchial tree", "Lung_Integrated",
         "lung", -0.038079, 0.02488, 26, "human"),
        ("Seattle Alzheimer brain atlas (SEA-AD) L5 extratelencephalic projecting "
         "glutamatergic cortical neuron",
         "L5 extratelencephalic projecting glutamatergic cortical neuron",
         "SEAAD_DLPFC", "brain", -0.021299, 0.02488, 51, "human"),
        ("Mount Sinai brain samples (MSSM) sst GABAergic cortical interneuron",
         "sst GABAergic cortical interneuron", "MSSM_DLPFC", "brain",
         -0.012646, 0.04678, 77, "human"),
        ("Mount Sinai brain samples (MSSM) oligodendrocyte precursor cell",
         "oligodendrocyte precursor cell", "MSSM_DLPFC", "brain",
         -0.011896, 0.00424, 81, "human"),
        ("Mount Sinai brain samples (MSSM) GABAergic neuron",
         "GABAergic neuron", "MSSM_DLPFC", "brain",
         -0.011759, 0.02944, 82, "human"),
        ("Mount Sinai brain samples (MSSM) oligodendrocyte", "oligodendrocyte",
         "MSSM_DLPFC", "brain", 0.017269, 0.00450, 219, "human"),
        ("Mount Sinai brain samples (MSSM) L6b glutamatergic neuron",
         "L6b glutamatergic neuron of the primary motor cortex", "MSSM_DLPFC",
         "brain", 0.019765, 0.00421, 227, "human"),
        ("Human lung conventional dendritic cell", "conventional dendritic cell",
         "Lung_Integrated", "lung", 0.022897, 0.01449, 237, "human"),
        ("Human kidney cortex collecting duct epithelial cell",
         "kidney cortex collecting duct epithelial cell", "Kidney_LivingDonor",
         "kidney", 0.080346, 0.04085, 282, "human"),
    ]
    # Same-scan rows live in the secondary-checks file: all five of them are
    # among the 13 targets below.
    sec = rows("results/followup-mature-composition-secondary-checks.csv")
    assert len(sec) == 5, len(sec)
    assert all(r["status"] == "SCAN_FINDING" for r in sec)
    sec_keys = {(r["cell_type"], r["cohort"]) for r in sec}

    # The primary scan family is the pooled mature 577-test family
    # (q_pooled_mature_bh below 0.05; 39 passes overall, 13 mammalian:
    # 2 mouse and 11 human).
    passing = set()
    for r in mouse + human:
        pr = pooled_by_key[(r["species"], r["organ"], r["cohort"], r["cell_type"])]
        if float(pr["q_pooled_mature_bh"]) < 0.05:
            passing.add((r["cell_type"], r["cohort"], r["organ"]))
    found = {}
    for label, cell, cohort, organ, expected, exp_q, _exp_rank, species in targets:
        df = mouse if species == "mouse" else human
        r = find(df, cell, cohort, organ)
        pr = pooled_by_key[(r["species"], organ, cohort, cell)]
        q_val = float(pr["q_pooled_mature_bh"])
        assert abs(float(r["effect"]) - expected) < 5e-7, (label, r["effect"], expected)
        assert q_val < 0.05, (label, q_val)
        assert abs(q_val - exp_q) < 5e-6, (label, q_val, exp_q)
        found[label] = (float(r["effect"]), species, q_val)
    assert passing == {(c, co, o) for (_, c, co, o, _, _, _, _) in targets}, passing
    assert sec_keys <= {(c, co) for (_, c, co, _, _, _, _, _) in targets}

    # Rank within each species panel, most negative first (1 = most negative).
    # The printed number next to each highlighted dot IS its rank on that
    # panel's axis.
    ranks = {}
    for species, df in (("mouse", mouse), ("human", human)):
        ordered = sorted(df, key=lambda r: float(r["effect"]))
        for i, r in enumerate(ordered, start=1):
            ranks[(r["cell_type"], r["cohort"], r["organ"])] = i
    numbers = {}
    for label, cell, cohort, organ, _expected, _exp_q, exp_rank, _species in targets:
        rk = ranks[(cell, cohort, organ)]
        assert rk == exp_rank, (label, rk, exp_rank)
        numbers[label] = rk

    n_mouse_pass = sum(1 for t in targets if t[7] == "mouse")
    n_human_pass = sum(1 for t in targets if t[7] == "human")
    assert (n_mouse_pass, n_human_pass) == (2, 11)
    n_other_mouse = len(mouse) - n_mouse_pass
    n_other_human = len(human) - n_human_pass
    assert (n_other_mouse, n_other_human) == (112, 274)

    print("S2 highlighted rows (pooled 577-test family, 13 mammalian passes):")
    for label, cell, cohort, organ, _expected, _exp_q, _exp_rank, species in targets:
        val, sp, qv = found[label]
        print(f"  {sp:5s} rank={numbers[label]:3d} | {organ:11s} | {cohort:18s} | {cell} | "
              f"effect={val:+.6f} | pooled_q={qv:.10f}")

    # Explicit per-rank offsets (in points) so clustered labels (21/26,
    # 77/81/82, 219/227/237, and 282) never collide with one another,
    # the x=0 axis line, the legend, or the curve.
    # Tuple: (dx_pt, dy_pt, ha, va, draw_leader)
    label_offsets = {
        ("mouse", 8): (-10, 7, "right", "center", True),
        ("mouse", 15): (8, 1, "left", "center", False),
        ("human", 6): (7, -1, "left", "center", False),
        ("human", 21): (-10, 8, "right", "center", True),
        ("human", 26): (8, -2, "left", "center", False),
        ("human", 51): (8, 0, "left", "center", False),
        ("human", 77): (-14, -12, "right", "center", True),
        ("human", 81): (-20, 2, "right", "center", True),
        ("human", 82): (-14, 17, "right", "center", True),
        ("human", 219): (14, -10, "left", "center", True),
        ("human", 227): (18, 1, "left", "center", True),
        ("human", 237): (14, 10, "left", "center", True),
        ("human", 282): (0, -14, "center", "top", True),
    }

    fig = plt.figure(figsize=(11.8, 9.8))
    gs = fig.add_gridspec(2, 1, height_ratios=[2.9, 2.05], hspace=0.58)
    gs_top = gs[0].subgridspec(1, 2, wspace=0.28)
    ax_m = fig.add_subplot(gs_top[0])
    ax_h = fig.add_subplot(gs_top[1])
    panels = [
        (ax_m, mouse, "mouse", "Mouse rows sorted by slope (1 = most negative)",
         "CLR slope per month", [-0.50, -0.25, 0.00, 0.25], (-0.60, 0.30),
         [1, 30, 60, 90, 114], (0, 122), f"Other {n_other_mouse} mouse rows"),
        (ax_h, human, "human", "Human rows sorted by slope (1 = most negative)",
         "CLR slope per year", [-0.10, -0.05, 0.00, 0.05, 0.10], (-0.13, 0.12),
         [1, 100, 200, 285], (0, 300), f"Other {n_other_human} human rows"),
    ]
    panel_annots = {"mouse": [], "human": []}
    for ax, df, species, ylabel, xlabel, xticks, xlim, yticks, ylim, other_label in panels:
        style_ax(ax)
        slopes = sorted(float(r["effect"]) for r in df)
        n = len(df)
        ax.scatter(slopes, range(1, n + 1), s=14, color="#8a8a8a", zorder=2,
                   label=other_label)
        for label, cell, cohort, organ, _expected, _exp_q, _exp_rank, sp in targets:
            if sp != species:
                continue
            val, _, _ = found[label]
            rank = numbers[label]
            col = MOUSE if species == "mouse" else HUMAN
            ax.scatter([val], [rank], s=90, color=col, edgecolor="white",
                       linewidth=1.0, zorder=4)
            dx_pt, dy_pt, ha, va, leader = label_offsets[(species, rank)]
            arrowprops = (dict(arrowstyle="-", color=col, lw=0.8, shrinkA=0, shrinkB=4)
                          if leader else None)
            ann = ax.annotate(
                str(rank),
                xy=(val, rank),
                xytext=(dx_pt, dy_pt),
                textcoords="offset points",
                fontsize=9.5,
                fontweight="bold",
                color=col,
                va=va,
                ha=ha,
                arrowprops=arrowprops,
                zorder=5,
            )
            panel_annots[species].append((rank, ann))
        ax.axvline(0, color="black", linewidth=1.0, zorder=1)
        ax.set_title("Mouse: 114 rows" if species == "mouse" else "Human: 285 rows",
                     fontsize=12)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xticks(xticks)
        ax.set_xlim(*xlim)
        ax.set_yticks(yticks)
        ax.set_ylim(*ylim)
        ax.legend(loc="lower right", frameon=False)
    fig.suptitle("All 399 cell-type fraction slopes, split by species; the 13 that clear the bar are highlighted",
                 fontsize=12, y=0.98)
    footnote = fig.text(
        0.5, 0.395,
        "The 13 highlighted rows clear the bar: q below 0.05 in the pooled 577-test family. "
        "Numbers are each row's rank in its own panel (1 = most negative).\n"
        "Mouse slopes are per month; human slopes are per year. "
        "MSSM: Mount Sinai brain samples. SEA-AD: Seattle Alzheimer brain atlas.",
        ha="center", va="top", fontsize=10,
    )

    axk = fig.add_subplot(gs[1])
    axk.axis("off")
    axk.set_xlim(0, 13.05)
    axk.set_ylim(0, 10)
    key_header = axk.text(
        0, 9.85, "The 13 highlighted rows (number = rank in its panel):",
        fontsize=11, fontweight="bold", va="top", ha="left",
    )
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
        "Mount Sinai brain samples (MSSM) sst GABAergic cortical interneuron":
            "MSSM sst GABAergic cortical interneuron",
        "Mount Sinai brain samples (MSSM) oligodendrocyte precursor cell":
            "MSSM oligodendrocyte precursor cell",
        "Mount Sinai brain samples (MSSM) GABAergic neuron":
            "MSSM GABAergic neuron",
        "Mount Sinai brain samples (MSSM) oligodendrocyte": "MSSM oligodendrocyte",
        "Mount Sinai brain samples (MSSM) L6b glutamatergic neuron":
            "MSSM L6b glutamatergic neuron",
        "Human lung conventional dendritic cell": "Human lung conventional dendritic cell",
        "Human kidney cortex collecting duct epithelial cell":
            "Kidney cortex collecting duct epithelial cell",
    }
    import textwrap
    # Key order: mouse ranks first, then human ranks, each most negative first.
    # 3-column layout (5, 4, 4 entries) so all 13 fit legibly with no overlap.
    keyed = sorted(found.items(), key=lambda kv: (kv[1][1] != "mouse", numbers[kv[0]]))
    col_slices = [keyed[:5], keyed[5:9], keyed[9:]]
    col_xs = [0.0, 4.35, 8.70]
    key_texts = []
    for x0, col_items in zip(col_xs, col_slices):
        y = 8.70
        for lab, (val, species, _qv) in col_items:
            num = numbers[lab]
            sign = "+" if val > 0 else ""
            unit = "per month" if species == "mouse" else "per year"
            text = f"{num}. {short[lab]}, {sign}{val:.3f} {unit}"
            lines = textwrap.wrap(text, width=39)
            for k, ln in enumerate(lines):
                t = axk.text(x0 + (0.24 if k > 0 else 0), y, ln,
                             fontsize=8.8, va="top", ha="left")
                key_texts.append((text, t))
                y -= 0.66
            y -= 0.28

    # Verify rendered text bounding boxes: no overlaps and no overflow.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    # 1. Check rank label text bounding boxes do not overlap within either panel.
    for sp, ann_list in panel_annots.items():
        boxes = [(rk, matplotlib.text.Text.get_window_extent(ann, renderer))
                 for rk, ann in ann_list]
        for i in range(len(boxes)):
            for j in range(i + 1, len(boxes)):
                assert not boxes[i][1].overlaps(boxes[j][1]), (
                    f"Overlapping rank labels in {sp}: {boxes[i][0]} and {boxes[j][0]}"
                )
    # 2. Check footnote sits cleanly above the key header and below the x-labels.
    fn_box = footnote.get_window_extent(renderer)
    hdr_box = key_header.get_window_extent(renderer)
    xlab_y0 = min(ax_m.xaxis.label.get_window_extent(renderer).y0,
                  ax_h.xaxis.label.get_window_extent(renderer).y0)
    assert fn_box.y0 > hdr_box.y1 + 2.0, (fn_box, hdr_box)
    assert fn_box.y1 < xlab_y0 - 2.0, (fn_box, xlab_y0)
    # 3. Check all key text boxes are non-overlapping and stay within axk.
    axk_box = axk.get_window_extent(renderer)
    kboxes = [(lbl, t.get_window_extent(renderer)) for lbl, t in key_texts]
    for lbl, b in kboxes:
        assert b.x0 >= axk_box.x0 - 2.0 and b.x1 <= axk_box.x1 + 2.0, (lbl, b, axk_box)
        assert b.y0 >= axk_box.y0 - 2.0 and b.y1 <= hdr_box.y0 + 1.0, (lbl, b, axk_box, hdr_box)
    for i in range(len(kboxes)):
        for j in range(i + 1, len(kboxes)):
            assert not kboxes[i][1].overlaps(kboxes[j][1]), (
                f"Overlapping key texts: {kboxes[i][0]} vs {kboxes[j][0]}"
            )
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
    # Gao heart cohort, same rows: no clear change on either estimand.
    gdata = rows("results/followup-mature-between-subject-gao-models.csv")
    gsel = {}
    for r in gdata:
        if (r["is_stratified_global"] == "True" and r["depth"] == "full_1500"
                and r["model_spec"] == "stratified_global_primary"):
            gsel[r["part"]] = r
    assert abs(float(gsel["A1"]["slope_sd_per_decade"]) - (-0.007521)) < 5e-7
    assert abs(float(gsel["A1"]["pvalue"]) - 0.766955) < 1e-5
    assert abs(float(gsel["A3"]["slope_sd_per_decade"]) - (-0.061872)) < 5e-7
    assert abs(float(gsel["A3"]["pvalue"]) - 0.084606) < 1e-5

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
    ax.set_title("Donor-to-donor distances and within-donor spread shrink with age in blood",
                 fontsize=12)
    ax.text(0.0, -0.34,
            "Yazar blood cohort, 925 donors, full-depth data. Slopes are change per decade.\n"
            "Gao heart cohort: no clear change (slopes -0.008, p 0.77 and -0.062, p 0.085).",
            transform=ax.transAxes, fontsize=10, va="top", ha="left")
    save(fig, "s3-between-subject-shrinkage")


# ---------------------------------------------------------------- S5
def build_s4():
    # Both-test pass counts per signature, cohort, and arm, from the
    # benchmark CSV; denominators are the tested cell groups per cohort.
    data = rows("results/followup-mature-block10-benchmark-results.csv")
    assert len(data) == 246
    counts = {}
    units = {}
    for r in data:
        key = (r["signature_id"], r["cohort"], r["arm"])
        counts[key] = counts.get(key, 0) + (1 if r["both_pass"] == "True" else 0)
        units.setdefault((r["signature_id"], r["cohort"]), set()).add(r["unit"])
    expected = {
        ("Harries2011_6gene", "Gao_Heart", "unadjusted"): 0,
        ("Harries2011_6gene", "Gao_Heart", "batch_adjusted"): 0,
        ("Harries2011_6gene", "Yazar_Blood", "unadjusted"): 7,
        ("Harries2011_6gene", "Yazar_Blood", "batch_adjusted"): 5,
        ("Peters2015_1497gene", "Gao_Heart", "unadjusted"): 0,
        ("Peters2015_1497gene", "Gao_Heart", "batch_adjusted"): 0,
        ("Peters2015_1497gene", "Yazar_Blood", "unadjusted"): 16,
        ("Peters2015_1497gene", "Yazar_Blood", "batch_adjusted"): 16,
        ("Sood2015_150probe", "Gao_Heart", "unadjusted"): 3,
        ("Sood2015_150probe", "Gao_Heart", "batch_adjusted"): 0,
        ("Sood2015_150probe", "Yazar_Blood", "unadjusted"): 1,
        ("Sood2015_150probe", "Yazar_Blood", "batch_adjusted"): 0,
        ("deMagalhaes2009_73gene", "Gao_Heart", "unadjusted"): 0,
        ("deMagalhaes2009_73gene", "Gao_Heart", "batch_adjusted"): 0,
        ("deMagalhaes2009_73gene", "Yazar_Blood", "unadjusted"): 0,
        ("deMagalhaes2009_73gene", "Yazar_Blood", "batch_adjusted"): 0,
    }
    for key, val in expected.items():
        assert counts.get(key, 0) == val, (key, counts.get(key, 0), val)
    # Denominators: Harries tested 13 Gao heart groups; the other three
    # signatures tested 14. Yazar blood: 17 groups for all four.
    assert len(units[("Harries2011_6gene", "Gao_Heart")]) == 13
    for sig in ("Peters2015_1497gene", "Sood2015_150probe", "deMagalhaes2009_73gene"):
        assert len(units[(sig, "Gao_Heart")]) == 14, sig
    for sig in ("Harries2011_6gene", "Peters2015_1497gene", "Sood2015_150probe",
                "deMagalhaes2009_73gene"):
        assert len(units[(sig, "Yazar_Blood")]) == 17, sig

    sigs = [
        ("Harries2011_6gene", "Harries 2011\n6 genes"),
        ("Peters2015_1497gene", "Peters 2015\n1,497 genes"),
        ("Sood2015_150probe", "Sood 2015\n150 probes"),
        ("deMagalhaes2009_73gene", "de Magalhaes 2009\n73 genes"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 5.8), sharey=True)
    x = np.arange(len(sigs))
    w = 0.36
    for ax, cohort, title in (
            (axes[0], "Gao_Heart", "Gao heart"),
            (axes[1], "Yazar_Blood", "Yazar blood")):
        style_ax(ax)
        before = [counts[(sig, cohort, "unadjusted")] for sig, _ in sigs]
        after = [counts[(sig, cohort, "batch_adjusted")] for sig, _ in sigs]
        b1 = ax.bar(x - w / 2, before, w, color="#bdbdbd", label="Before batch adjustment",
                    zorder=3)
        b2 = ax.bar(x + w / 2, after, w, color=MOUSE, label="After batch adjustment",
                    zorder=3)
        for bars in (b1, b2):
            for b in bars:
                ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.4,
                        str(int(b.get_height())), ha="center", fontsize=11)
        ax.set_title(title, fontsize=12)
        ax.set_xticks(x)
        ax.set_xticklabels([lab for _, lab in sigs], fontsize=10)
        ax.set_yticks([0, 5, 10, 15])
        ax.set_ylim(0, 18.5)
        ax.grid(True, axis="y", color=GRID, linewidth=0.7, alpha=0.7)
        ax.set_axisbelow(True)
    axes[0].set_ylabel("Cell groups passing both tests", fontsize=11)
    axes[0].legend(loc="upper left", frameon=False)
    fig.suptitle("Four published aging signatures, tested in Gao heart and Yazar blood",
                 fontsize=12, y=0.98)
    fig.text(0.5, 0.02,
             "Passes both tests means the age signal beats both the composite test and the random-signature test.\n"
             "Harries 2011 tested 13 Gao heart groups; the other three signatures tested 14. Yazar blood: 17 groups for all four signatures.\n"
             "Sood 2015: the 4 groups that pass before adjustment are 3 Gao heart and 1 Yazar blood. "
             "Harries and Peters still pass after batch adjustment; Sood passes 0 after adjustment; de Magalhaes passes 0 in either arm.\n"
             "Four signatures cannot prove a claim about all published signatures.",
             ha="center", va="bottom", fontsize=10)
    fig.tight_layout(rect=[0, 0.14, 1, 0.94])
    save(fig, "s4-published-signatures")



def build_s5():
    data = rows("results/followup-mature-concordance-figure-data.csv")
    assert len(data) == 3, len(data)
    d = {r["species_scope"]: r for r in data}
    mw = d["Mouse, within tissue"]
    mc = d["Mouse, cross tissue"]
    hu = d["Human, within cohort"]
    assert abs(float(mw["p_value"]) - 0.0219780219780219) < 1e-9
    assert abs(float(hu["p_value"]) - 0.8151848151848152) < 1e-9
    assert abs(float(mw["observed"]) - 0.210467) < 1e-5
    assert abs(float(mw["null_ci95_low"]) - 0.101241) < 1e-5
    assert abs(float(mw["null_ci95_high"]) - 0.205516) < 1e-5
    assert float(mw["observed"]) > float(mw["null_ci95_high"])
    # q values use the two-test primary family (mouse within-tissue and
    # human), read from the BH check file. The three-test family q 0.066
    # does not exist in the analysis record and must not appear here.
    bh = rows("results/followup-mature-concordance-bh-check.csv")
    two = {r["test"]: float(r["q_bh"]) for r in bh
           if r["family"] == "mature_two_test_family"}
    assert len(two) == 2, two
    q_mouse = two["Mouse within_tissue primary_detection_filtered"]
    q_human = two["Human within_tissue_or_cohort primary_detection_filtered"]
    assert abs(q_mouse - 0.0439560439560438) < 1e-9, q_mouse
    assert abs(q_human - 0.8151848151848152) < 1e-9, q_human

    groups = [
        ("Mouse\nwithin tissue", float(mw["observed"]), float(mw["null_mean"]),
         float(mw["null_ci95_low"]), float(mw["null_ci95_high"]),
         f"q {q_mouse:.3f} (p 0.022)", MOUSE),
        ("Mouse\ncross tissue", float(mc["observed"]), float(mc["null_mean"]),
         float(mc["null_ci95_low"]), float(mc["null_ci95_high"]), "p 0.103", MOUSE),
        ("Human\nwithin cohort", float(hu["observed"]), float(hu["null_mean"]),
         float(hu["null_ci95_low"]), float(hu["null_ci95_high"]),
         f"q {q_human:.3f} (p 0.815)", MOUSE),  # observed bars all dark blue per legend
    ]
    fig, ax = plt.subplots(figsize=(9, 5.8))
    style_ax(ax)
    x = np.arange(3)
    w = 0.36
    ax.bar(x - w / 2, [g[1] for g in groups], w, color=[g[6] for g in groups],
           label="Observed", zorder=3)
    ax.bar(x + w / 2, [g[2] for g in groups], w, color="#bdbdbd",
           label="Joint donor-permutation null mean", zorder=3)
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
    ax.set_title("Age-signal agreement across cell types, observed versus joint donor-permutation null",
                 fontsize=12)
    ax.legend(loc="upper right", frameon=False)
    ax.text(0.0, -0.30,
            "Error bars: null 95% confidence interval. Mouse within-tissue: q 0.044, passes the 0.05 bar "
            "(false discovery rate correction over the two-test primary family: mouse within-tissue and human; p 0.022).\n"
            "Observed 0.2105 sits just above the null 95% interval (0.101 to 0.206). "
            "Mouse cross-tissue is a separate scope and carries no q in the primary family (p 0.103). Human: q 0.815, null.\n"
            "Mouse: 59 groups, 14 mice. Human: 127 groups, 416 donors, 7 cohorts.",
            transform=ax.transAxes, fontsize=10, va="top", ha="left")
    save(fig, "s5-concordance")



def build_s6():
    def get(path, gene_set, role):
        for r in rows(path):
            if (r["is_stratified_global"] == "True" and r["metric"] == "mean_ectopic_count"
                    and r["gene_set"] == gene_set and r["analysis_role"] == role):
                return r
        raise AssertionError(f"not found: {path} {gene_set} {role}")

    e2g = get("results/followup-mature-fulldepth-gao-ectopic.csv", "E2_stratum",
              "primary_depth_matched_e2")
    e2y = get("results/followup-mature-fulldepth-yazar-ectopic.csv", "E2_stratum",
              "primary_depth_matched_e2")
    for r, slope, p in [(e2g, -0.127, 0.031), (e2y, -0.033, 0.0031)]:
        assert abs(float(r["slope_sd_per_decade"]) - slope) < 5e-4, (r, slope)
        assert abs(float(r["pvalue"]) - p) < 5e-4, (r, p)

    pts = [
        ("Gao heart cohort\nstricter off-program count",
         e2g, "slope -0.127, p 0.031", MOUSE),
        ("Yazar blood cohort\nstricter off-program count",
         e2y, "slope -0.033, p 0.0031", HUMAN),
    ]
    fig, ax = plt.subplots(figsize=(9, 5.6))
    style_ax(ax)
    ys = [1, 0]
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
    ax.set_ylim(-0.6, 1.6)
    ax.set_xlabel("Change per decade of age (standard deviations)")
    ax.set_title("Stricter off-program count falls with age in both cohorts", fontsize=12)
    ax.text(0.0, -0.30,
            "Stricter count: genes silent in young cells of the same type, matched by depth.\n"
            "Gao 103 donors, Yazar 919 donors.",
            transform=ax.transAxes, fontsize=10, va="top", ha="left")
    save(fig, "s6-ectopic-e2")


if __name__ == "__main__":
    build_s1()
    build_s2()
    verify_dropped_transport()
    build_s3()
    build_s4()
    build_s5()
    build_s6()
    print("done")
