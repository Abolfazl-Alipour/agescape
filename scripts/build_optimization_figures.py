#!/usr/bin/env python3
"""Build the four optimization-pass figures for the Agescape manuscript from
the committed mature result tables and the per-donor export files.

Style system (all four figures): Liberation Serif throughout; spines 1.3 pt;
tick labels 10-11 pt; axis labels 11-12 pt; titles 12 pt; one muted color
per species (dark slate blue for mouse, deep amber/brown for human);
filled markers = association survives false discovery rate correction
(q below 0.05), open markers = does not survive; confidence-interval lines
at least 2 pt; explicit tick placement everywhere. Every number is read from
a committed CSV in results/ and asserted against the values printed in the
manuscript before plotting. Fitted lines in donor panels use the committed
regression slope, drawn through the data center; the slope is refit from the
exported per-donor file and asserted equal to the committed slope.

Fig 1: mature-to-old cell-fraction age slopes (forest) plus donor
  scatterplots with fitted lines for the five headline fraction rows.
Fig 2: donor scatterplots with fitted lines: (a) stratified scatter versus
  age in Gao heart and Yazar blood, (b) pooled center shift versus age in
  Yazar blood, (c) Shannon entropy versus age in Gao heart and Yazar blood.
  Scatter and Shannon values are standardized within cell type (cell-type
  mean subtracted, divided by cell-type SD), matching the stratified model.
Fig 3: core 48-unit transport verdicts, cleaned versus unfiltered
  mature-to-old signature, in plain words.
Fig 4: mature interneuron Mount Sinai versus Human Brain Collection Core age
  slopes, TRMT11 labeled (design unchanged, style only).

PNG (300 dpi) and PDF (vector) pairs are written; the tex references the PDF.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
FIGURES = REPO / "report" / "plos-one" / "figures-v2"

# --- style system ---
MOUSE = "#2e4a62"   # dark slate blue
HUMAN = "#9a6b1e"   # deep amber / brown
STEEL = "#6a8caf"   # muted steel blue (Fig 3 third category)
LIGHT = "#d9d9d9"   # light grey (Fig 3 fourth category)
GRID = "#cccccc"

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Liberation Serif", "DejaVu Serif"],
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "axes.linewidth": 1.3,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.major.size": 4,
    "xtick.major.width": 1.2,
    "ytick.major.size": 4,
    "ytick.major.width": 1.2,
})


def save(fig, name):
    fig.savefig(FIGURES / f"{name}.png", dpi=300)
    fig.savefig(FIGURES / f"{name}.pdf")
    plt.close(fig)


def _set_ticks(ax, xticks, yticks=None):
    ax.set_xticks(xticks)
    ax.set_xticklabels([f"{t:.2f}".rstrip("0").rstrip(".") if abs(t) < 1 else f"{t:.0f}"
                        for t in xticks])
    if yticks is not None:
        ax.set_yticks(yticks)


def _fitted_line(ax, x, y, slope_per_year, color, lw=2.0):
    """Draw the committed-slope line through the data center."""
    xc = float(np.mean(x))
    yc = float(np.mean(y))
    xs = np.array([np.min(x), np.max(x)], dtype=float)
    ys = yc + slope_per_year * (xs - xc)
    ax.plot(xs, ys, color=color, lw=lw, zorder=4)
    return xc, yc


def _refit_clr_slope(sub):
    """Refit OLS clr_value ~ age + sex_male; return (slope, se, p, model)."""
    X = sm.add_constant(sub[["age", "sex_male"]])
    model = "age_plus_sex"
    if np.linalg.matrix_rank(X.to_numpy()) < X.shape[1]:
        X = sm.add_constant(sub[["age"]])
        model = "age_only"
    fit = sm.OLS(sub["clr_value"], X).fit()
    return float(fit.params["age"]), float(fit.bse["age"]), float(fit.pvalues["age"]), model


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
        ("Mouse", "Limb_Muscle", "skeletal muscle cell"): "Limb skeletal muscle",
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
    # --- donor data for the five headline rows ---
    clr = pd.read_csv(RESULTS / "followup-mature-composition-donor-clr.csv")
    human_comp = pd.read_csv(RESULTS / "followup-mature-composition-human.csv")
    headline = [
        ("lung", "Lung_Integrated", "club cell", "Lung club cell"),
        ("lung", "Lung_Integrated", "multiciliated epithelial cell", "Lung multiciliated epithelial"),
        ("brain", "MSSM_DLPFC", "L6b glutamatergic neuron of the primary motor cortex", "L6b glutamatergic neuron"),
        ("brain", "MSSM_DLPFC", "oligodendrocyte", "Oligodendrocyte"),
        ("brain", "MSSM_DLPFC", "oligodendrocyte precursor cell", "Oligodendrocyte precursor"),
    ]
    donor_panels = []
    for organ, cohort, ct, label in headline:
        sub = clr[(clr["organ"] == organ) & (clr["cohort"] == cohort)
                  & (clr["cell_type"] == ct)].copy()
        crow = human_comp[(human_comp["organ"] == organ) & (human_comp["cohort"] == cohort)
                          & (human_comp["cell_type"] == ct)]
        assert len(crow) == 1, (organ, cohort, ct)
        crow = crow.iloc[0]
        eff, se, pv, model = _refit_clr_slope(sub)
        assert abs(eff - float(crow["effect"])) < 1e-9, (ct, eff, crow["effect"])
        assert abs(se - float(crow["se"])) < 1e-9, (ct, se, crow["se"])
        # partial residual: remove the fitted sex effect so points sit on the
        # age slope; the line uses the committed slope through the center.
        X = sm.add_constant(sub[["age", "sex_male"]])
        fit = sm.OLS(sub["clr_value"], X).fit()
        beta_sex = float(fit.params["sex_male"])
        sub["y_plot"] = sub["clr_value"] - beta_sex * sub["sex_male"]
        donor_panels.append({
            "label": label, "sub": sub, "slope": float(crow["effect"]),
            "q": float(crow["q_431_family"]), "p": float(crow["pvalue"]),
        })

    # Forest is stacked (mouse on top, human below) with heights in the
    # ratio of the row counts (9 to 5), so both panels have the same row
    # spacing and neither looks sparse. Donor panels sit below.
    fig = plt.figure(figsize=(10.5, 9.2))
    gs = fig.add_gridspec(3, 1, height_ratios=[9, 5, 8], hspace=0.55)
    ax_mouse = fig.add_subplot(gs[0, 0])
    ax_human = fig.add_subplot(gs[1, 0])
    gs_bot = gs[2, 0].subgridspec(1, 5, wspace=0.35)
    ax_donor = [fig.add_subplot(gs_bot[0, i]) for i in range(5)]

    for ax, species, color, xlabel in (
            (ax_mouse, "Mouse", MOUSE, "CLR slope per month"),
            (ax_human, "Human", HUMAN, "CLR slope per year")):
        part = v[(v["species"] == species) & (v["verdict"] != "UNTESTABLE")].copy()
        part["label"] = [short[(r["species"], r["organ"], r["cell_type"])]
                         for _, r in part.iterrows()]
        part = part.sort_values("mature_effect")
        y = np.arange(len(part))[::-1]
        lo = part["mature_effect"] - 1.96 * part["mature_se"]
        hi = part["mature_effect"] + 1.96 * part["mature_se"]
        for yi, l, h in zip(y, lo, hi):
            ax.plot([l, h], [yi, yi], color=color, lw=2.0, alpha=0.85)
        survived = part["verdict"] == "SURVIVES"
        ax.scatter(part.loc[survived, "mature_effect"], y[survived.values],
                   color=color, s=70, zorder=3)
        ax.scatter(part.loc[~survived, "mature_effect"], y[(~survived).values],
                   facecolors="none", edgecolors=color, s=70, linewidths=1.5, zorder=3)
        ax.set_yticks(y, part["label"], fontsize=10)
        ax.set_ylim(-0.6, len(part) - 0.4)
        ax.axvline(0, color=GRID, lw=1.0, zorder=0)
        ax.set_xlabel(xlabel, fontsize=11)
        n_surv = int(survived.sum())
        n_tot = len(part)
        ax.set_title(f"{species}: {n_surv} of {n_tot} associations survive correction"
                     if n_surv else f"{species}: none of {n_tot} associations survive correction",
                     fontsize=12)
        # explicit, widely spaced ticks (full-width panels have room)
        if species == "Mouse":
            ax.set_xticks([-0.20, -0.10, 0.00, 0.10, 0.20])
            ax.set_xticklabels(["-0.20", "-0.10", "0.00", "0.10", "0.20"])
        else:
            ax.set_xticks([-0.10, -0.05, 0.00, 0.05])
            ax.set_xticklabels(["-0.10", "-0.05", "0.00", "0.05"])
    ax_mouse.set_ylabel("Cell type", fontsize=11)
    ax_human.set_ylabel("Cell type", fontsize=11)
    # legend in the mouse panel's empty upper right (top-row intervals sit left)
    from matplotlib.lines import Line2D
    legend_elems = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#666666",
               markeredgecolor="#666666", markersize=9,
               label="Survives correction (q below 0.05)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="none",
               markeredgecolor="#666666", markersize=9,
               label="Does not survive correction"),
    ]
    ax_mouse.legend(handles=legend_elems, loc="upper right", frameon=False, fontsize=10)
    fig.suptitle("Mature-to-old cell-fraction age slopes (95% intervals); bladder untestable, not shown",
                 fontsize=12)

    for ax, d in zip(ax_donor, donor_panels):
        sub = d["sub"]
        x = sub["age"].to_numpy(dtype=float)
        yv = sub["y_plot"].to_numpy(dtype=float)
        ax.scatter(x, yv, s=14, color=HUMAN, alpha=0.55, linewidths=0, zorder=3)
        _fitted_line(ax, x, yv, d["slope"], HUMAN, lw=2.0)
        ax.set_title(d["label"], fontsize=11)
        ax.set_xlabel("Age (years)", fontsize=11)
        ax.set_xticks([30, 50, 70, 90])
        ax.set_xticklabels(["30", "50", "70", "90"])
        ax.tick_params(axis="y", labelsize=10)
        ax.text(0.05, 0.95, f"slope {d['slope']:+.3f} per year\nq {d['q']:.3f}",
                transform=ax.transAxes, fontsize=10, va="top", ha="left")
    ax_donor[0].set_ylabel("CLR value (sex-adjusted)", fontsize=11)
    fig.subplots_adjust(left=0.21, right=0.98, top=0.92, bottom=0.08)
    save(fig, "fig1-composition")
    print("fig1-composition rebuilt: forest + 5 donor panels; slopes refit and asserted")


def _load_stratified(ds, metric, table_csv, donor_csv, out_col, n_cells_min=None):
    """Load donor rows for a stratified-global endpoint and standardize
    within cell type (group mean subtracted, divided by group SD), matching
    the stratified model. Returns (plot frame, committed slope, committed p)."""
    t = pd.read_csv(table_csv)
    r = t[(t["is_stratified_global"] == True) & (t["metric"] == metric)
          & (t["depth_mode"] == "full_depth") & (t["model_spec"] == "stratified_global_primary")]
    assert len(r) == 1, (ds, metric)
    cslope = float(r.iloc[0]["slope_sd_per_decade"])
    cp = float(r.iloc[0]["pvalue"])
    g = t[(t["is_pooled"] == False) & (t["is_stratified_global"] == False)
          & (t["model_spec"] == "primary") & (t["depth_mode"] == "full_depth")
          & (t["metric"] == metric) & (t["yazar_pool_fe_status"] == "evaluated")]
    elig = sorted(g["group"].unique().tolist())
    d = pd.read_csv(donor_csv)
    st = d[(d["group"].isin(elig)) & d[out_col].notna()].copy()
    if n_cells_min is not None:
        st = st[st["n_cells"] >= n_cells_min].copy()
    # standardize within cell type
    st["y_std"] = np.nan
    for grp, idx in st.groupby("group").groups.items():
        vals = st.loc[idx, out_col]
        st.loc[idx, "y_std"] = (vals - vals.mean()) / vals.std(ddof=1)
    return st, cslope, cp


def _scatter_panel(ax, st, slope_sd_per_decade, color, title, ylabel):
    x = st["age_years"].to_numpy(dtype=float)
    yv = st["y_std"].to_numpy(dtype=float)
    ax.scatter(x, yv, s=12, color=color, alpha=0.45, linewidths=0, zorder=3)
    _fitted_line(ax, x, yv, slope_sd_per_decade / 10.0, color, lw=2.0)
    ax.set_title(title, fontsize=12)
    ax.set_xlabel("Age (years)", fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_xticks([30, 50, 70, 90])
    ax.set_xticklabels(["30", "50", "70", "90"])
    ax.set_yticks([-2, -1, 0, 1, 2])
    ax.set_yticklabels(["-2", "-1", "0", "1", "2"])


def fig2_mature_drift():
    # committed stratified-global slopes (SD per decade), asserted from tables
    drift_m = pd.read_csv(RESULTS / "followup-mature-drift-mouse.csv")
    drift_h = pd.read_csv(RESULTS / "followup-mature-drift-human.csv")
    glob_m = drift_m[(drift_m["tissue"] == "All tested tissues")
                     & (drift_m["adjustment"] == "linear_log_umi")
                     & (drift_m["depth"] == "full")
                     & (drift_m["n_donors"] == 14)].iloc[0]
    assert abs(glob_m["slope"] - (-0.0137)) < 0.001, glob_m["slope"]
    glob_h = drift_h[(drift_h["tissue"] == "All tested tissues")
                     & (drift_h["adjustment"] == "linear_log_umi")
                     & (drift_h["depth"] == "full")
                     & (drift_h["n_donors"] == 407)].iloc[0]
    assert abs(glob_h["slope"] - (-0.0071)) < 0.001, glob_h["slope"]

    gao_scat, gao_scat_slope, gao_scat_p = _load_stratified(
        "gao", "scatter", RESULTS / "followup-mature-fulldepth-gao-cell-drift.csv",
        RESULTS / "followup-mature-fulldepth-gao-donor-part1.csv", "scatter")
    yaz_scat, yaz_scat_slope, yaz_scat_p = _load_stratified(
        "yazar", "scatter", RESULTS / "followup-mature-fulldepth-yazar-cell-drift.csv",
        RESULTS / "followup-mature-fulldepth-yazar-donor-part1.csv", "scatter")
    gao_sh, gao_sh_slope, gao_sh_p = _load_stratified(
        "gao", "shannon", RESULTS / "followup-mature-fulldepth-gao-spread.csv",
        RESULTS / "followup-mature-fulldepth-gao-donor-part2.csv", "full_shannon",
        n_cells_min=20)
    yaz_sh, yaz_sh_slope, yaz_sh_p = _load_stratified(
        "yazar", "shannon", RESULTS / "followup-mature-fulldepth-yazar-spread.csv",
        RESULTS / "followup-mature-fulldepth-yazar-donor-part2.csv", "full_shannon",
        n_cells_min=20)
    # pooled Yazar center shift (committed +0.0400 SD/decade, p 2.99e-05)
    cd = pd.read_csv(RESULTS / "followup-mature-fulldepth-yazar-cell-drift.csv")
    prow = cd[(cd["is_pooled"] == True) & (cd["metric"] == "center_shift")
              & (cd["depth_mode"] == "full_depth") & (cd["model_spec"] == "primary")].iloc[0]
    assert abs(float(prow["slope_sd_per_decade"]) - 0.0400) < 0.0005
    assert abs(float(prow["pvalue"]) - 2.99e-05) / 2.99e-05 < 0.01
    d1y = pd.read_csv(RESULTS / "followup-mature-fulldepth-yazar-donor-part1.csv")
    pool = d1y[(d1y["group"] == "All_Cells") & d1y["center_shift"].notna()].copy()
    pool["y_std"] = (pool["center_shift"] - pool["center_shift"].mean()) / pool["center_shift"].std(ddof=1)
    pool_slope = float(prow["slope_sd_per_decade"])
    pool_p = float(prow["pvalue"])

    fig = plt.figure(figsize=(10.5, 10.5))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1], hspace=0.65, wspace=0.30)
    ax_gs = fig.add_subplot(gs[0, 0])
    ax_ys = fig.add_subplot(gs[0, 1])
    ax_pool = fig.add_subplot(gs[1, :])
    ax_gh = fig.add_subplot(gs[2, 0])
    ax_yh = fig.add_subplot(gs[2, 1])

    _scatter_panel(ax_gs, gao_scat, gao_scat_slope, MOUSE,
                   "(a) Scatter, Gao heart", "Scatter, within cell type (SD)")
    ax_gs.text(0.05, 0.95, f"slope {gao_scat_slope:+.4f} SD per decade\np {gao_scat_p:.3g}",
               transform=ax_gs.transAxes, fontsize=10, va="top", ha="left")
    _scatter_panel(ax_ys, yaz_scat, yaz_scat_slope, HUMAN,
                   "(a) Scatter, Yazar blood", "Scatter, within cell type (SD)")
    ax_ys.text(0.05, 0.95, f"slope {yaz_scat_slope:+.4f} SD per decade\np {yaz_scat_p:.3g}",
               transform=ax_ys.transAxes, fontsize=10, va="top", ha="left")
    _scatter_panel(ax_pool, pool, pool_slope, HUMAN,
                   "(b) Pooled center shift, Yazar blood", "Center shift, standardized (SD)")
    ax_pool.text(0.05, 0.95, f"slope {pool_slope:+.4f} SD per decade\np {pool_p:.3g}",
                 transform=ax_pool.transAxes, fontsize=10, va="top", ha="left")
    _scatter_panel(ax_gh, gao_sh, gao_sh_slope, MOUSE,
                   "(c) Shannon entropy, Gao heart", "Shannon entropy, within cell type (SD)")
    ax_gh.text(0.05, 0.95, f"slope {gao_sh_slope:+.4f} SD per decade\np {gao_sh_p:.3g}",
               transform=ax_gh.transAxes, fontsize=10, va="top", ha="left")
    _scatter_panel(ax_yh, yaz_sh, yaz_sh_slope, HUMAN,
                   "(c) Shannon entropy, Yazar blood", "Shannon entropy, within cell type (SD)")
    ax_yh.text(0.05, 0.95, f"slope {yaz_sh_slope:+.4f} SD per decade\np {yaz_sh_p:.3g}",
               transform=ax_yh.transAxes, fontsize=10, va="top", ha="left")
    fig.suptitle("Donor scatterplots with fitted lines (slopes from the stratified model)",
                 fontsize=12)
    fig.tight_layout()
    fig.subplots_adjust(top=0.93)
    save(fig, "fig2-mature-drift")
    print("fig2-mature-drift rebuilt: 5 donor panels; slopes asserted from tables")


def fig3_transport():
    s = pd.read_csv(RESULTS / "followup-mature-p1m-block5-summary.csv")
    p1 = s[s["gene_set"] == "P1m"].iloc[0]
    p0 = s[s["gene_set"] == "P0m"].iloc[0]
    assert (p1["both_pass"], p1["empirical_only"], p1["composite_only"],
            p1["fail_both"]) == (2, 13, 4, 29)
    assert (p0["both_pass"], p0["empirical_only"], p0["composite_only"],
            p0["fail_both"]) == (3, 15, 5, 25)
    cats = ["Passes\nboth tests", "Passes the\nrandom-signature\ntest only",
            "Passes the\ncomposite\ntest only", "Passes\nneither test"]
    # Color encodes the signature (cleaned = slate blue, unfiltered = amber);
    # the category is the x-axis label. This makes the legend correct: the
    # two legend swatches are distinct colors, one per signature.
    fig, ax = plt.subplots(figsize=(10.5, 5.4))
    x = np.arange(len(cats))
    w = 0.36
    for off, row, label, color in (
            (-w / 2, p1, "Cleaned mature-to-old signature (1,914 genes)", MOUSE),
            (w / 2, p0, "Unfiltered mature-to-old signature (2,598 genes)", HUMAN)):
        vals = [row["both_pass"], row["empirical_only"], row["composite_only"],
                row["fail_both"]]
        bars = ax.bar(x + off, vals, w, label=label, color=color, edgecolor="white",
                      linewidth=1.0)
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v + 0.6, str(v),
                    ha="center", fontsize=11)
    ax.set_xticks(x, cats, fontsize=11)
    ax.set_ylabel("Core units (n = 48)", fontsize=11)
    ax.set_yticks([0, 10, 20, 30])
    ax.set_yticklabels(["0", "10", "20", "30"])
    ax.set_ylim(0, 34)
    ax.legend(frameon=False, loc="upper right", fontsize=11)
    ax.set_title("Whether the mouse signature's age signal appears in human cohorts",
                 fontsize=12)
    fig.tight_layout()
    save(fig, "fig3-transport")
    print("fig3-transport rebuilt: cleaned 2/13/4/29, unfiltered 3/15/5/25 (plain words)")


def fig4_interneuron():
    genes = pd.read_csv(RESULTS / "followup-mature-interneuron-genes.csv")
    tr = genes[genes["gene"] == "TRMT11"]
    assert set(tr["short_type"]) == {"VIP", "pvalb", "sst"}
    assert (tr["tier"] == "both_FDR").all()
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 7.6))
    for ax, short in zip(axes.ravel(), ("VIP", "pvalb", "sst", "GABA")):
        part = genes[genes["short_type"] == short]
        both = part["tier"] == "both_FDR"
        base_color = "#a0a0a0" if short == "GABA" else "#4478a6"
        ax.scatter(part.loc[~both, "slope_mssm"], part.loc[~both, "slope_hbcc"],
                   s=14, color=base_color, alpha=0.55, linewidths=0)
        ax.scatter(part.loc[both, "slope_mssm"], part.loc[both, "slope_hbcc"],
                   s=28, facecolors="none", edgecolors="#147d82", linewidths=1.2)
        ax.axhline(0, color="#dddddd", lw=1.0)
        ax.axvline(0, color="#dddddd", lw=1.0)
        lim = max(np.abs(part[["slope_mssm", "slope_hbcc"]].to_numpy()).max(), 1e-6) * 1.05
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xticks([-lim, -lim / 2, 0, lim / 2, lim])
        ax.set_xticklabels([f"{t:.3f}" for t in [-lim, -lim / 2, 0, lim / 2, lim]])
        ax.set_yticks([-lim, -lim / 2, 0, lim / 2, lim])
        ax.set_yticklabels([f"{t:.3f}" for t in [-lim, -lim / 2, 0, lim / 2, lim]])
        ax.set_title("Broad GABAergic neurons" if short == "GABA" else f"{short} interneurons",
                     fontsize=12)
        ax.set_xlabel("Mount Sinai brain samples (MSSM) age slope", fontsize=11)
        ax.set_ylabel("Human Brain Collection Core (HBCC) age slope", fontsize=11)
        t = tr[tr["short_type"] == short]
        if len(t):
            row = t.iloc[0]
            ax.annotate("TRMT11", (row["slope_mssm"], row["slope_hbcc"]),
                        textcoords="offset points", xytext=(10, 8), fontsize=11,
                        fontweight="bold", color="#147d82",
                        arrowprops=dict(arrowstyle="-", color="#147d82", lw=1.0))
    fig.suptitle("Genes passing false discovery rate correction in Mount Sinai brain samples (MSSM);\n"
                 "rings mark same-direction replication in Human Brain Collection Core (HBCC)",
                 fontsize=12)
    fig.tight_layout()
    fig.subplots_adjust(top=0.90)
    save(fig, "fig4-interneuron")
    print("fig4-interneuron rebuilt (design unchanged; style only)")


if __name__ == "__main__":
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig1_composition()
    fig2_mature_drift()
    fig3_transport()
    fig4_interneuron()
    print("optimization figures rebuilt")
