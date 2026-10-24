#!/usr/bin/env python3
"""
run_m9_pipeline.py — Human lifespan program-scoring on public CELLxGENE cohorts.

Run with the census environment (has cellxgene-census; no scanpy needed):
    /tmp/censusenv/bin/python scripts/run_m9_pipeline.py

M9 question: does the 261-gene mouse aging program (167 up / 94 down, results/
m6-program-genes.csv, the exact set scored in the app) track donor age in human
single-cell data with real per-donor ages?

Cohorts (CELLxGENE census stable, all public, no auth): control donors
(disease == "normal") with numeric ages from HsapDv development_stage.
Brain: MSSM/HBCC/RADC/SEA-AD DLPFC. Kidney: normal adult / living donor /
adult kidney v1.5 (snRNA + scRNA). Lung: CellRef 1.0 / integrated atlas
(normal donors) / CZI immune-aging BAL.

Outputs (repo-rooted):
  results/m9-donor-scores.csv       one row per cohort x donor x cell type
  results/m9-regressions.csv        one row per cohort x cell type: OLS score~age(+sex)
  results/m9-cohort-summary.csv     per-cohort headline (whole-donor score regression)
  results/m9-provenance.log         cohort composition + coverage decisions
  figures/m9-*.png                  score-vs-age panels

Honesty rules (same as other milestones): no invented numbers; coverage reported
per cohort; genes absent from a cohort's measurement are excluded there and the
coverage is recorded; 3 genes map to both directions in human space (CKB, ID2,
STMN1). They remain in unequally sized up/down means and do not simply cancel.
"""
from __future__ import annotations

import json
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats as stats

import cellxgene_census
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
FIGURES = REPO / "figures"
RESULTS.mkdir(exist_ok=True)
FIGURES.mkdir(exist_ok=True)

PROV = []


def prov(msg: str) -> None:
    line = f"[m9] {msg}"
    print(line, flush=True)
    PROV.append(line)


HUMAN_MAX_LIFESPAN = 122.5  # AnAge Build 14 (results/m6-species-lifespans.csv)

# cohorts: dataset_id -> (label, group)
COHORTS = {
    "37a17b78-4864-4a42-b67b-31c00962795a": ("MSSM_DLPFC", "brain"),
    "5c97eeeb-7e52-44b3-b010-b832b1f5424c": ("HBCC_DLPFC", "brain"),
    "5e57cd50-8e42-42d6-940d-5c1660d06864": ("RADC_DLPFC", "brain"),
    "6f7fd0f1-a2ed-4ff1-80d3-33dde731cbc3": ("SEAAD_DLPFC", "brain"),
    "0b4a15a7-4e9e-4555-9733-2423e5c66469": ("Kidney_NormalAdult", "kidney"),
    "0bae7ebf-eb54-46a6-be9a-3461cecefa4c": ("Kidney_LivingDonor", "kidney"),
    "a12ccb9b-4fbe-457d-8590-ac78053259ef": ("Kidney_SN_v15", "kidney"),
    "dea717d4-7bc0-4e46-950f-fd7e1cc8df7d": ("Kidney_SCR_v15", "kidney"),
    "eb499fd8-7000-419f-8854-926b9b61b11f": ("Lung_CellRef", "lung"),
    "9f222629-9e39-47d0-b83f-e08d610c7479": ("Lung_Integrated", "lung"),
    "1b350d0a-4535-4879-beb6-1142f3f94947": ("ImmuneAging_BAL", "lung"),
}

MIN_CELLS_PER_DONOR_CT = 20   # min cells in a donor x cell-type unit to score it
MIN_DONORS_FOR_REGRESSION = 6  # pre-registered human floor (M4)
MIN_AGE_SD = 12.0              # pre-registered human floor (M4)

AGE_RE = re.compile(r"(\d+)-year-old")


# ---------------------------------------------------------------- program genes

def load_program():
    pg = pd.read_csv(RESULTS / "m6-program-genes.csv")
    up_m = pg["program_up"].dropna().tolist()
    dn_m = pg["program_down"].dropna().tolist()
    orth = pd.read_csv(REPO / "data" / "interim" / "orthologs-mouse-human.csv")
    m2h = dict(zip(orth["mouse_gene_symbol"], orth["human_gene_symbol"]))
    up = sorted({m2h[g] for g in up_m if g in m2h and isinstance(m2h[g], str)})
    dn = sorted({m2h[g] for g in dn_m if g in m2h and isinstance(m2h[g], str)})
    prov(f"program: {len(up_m)} up / {len(dn_m)} down mouse genes -> "
         f"{len(up)} up / {len(dn)} down human orthologs")
    prov(f"human-space direction conflicts retained in unequal up/down means: "
         f"{sorted(set(up) & set(dn))}")
    return up, dn


# ---------------------------------------------------------------- stats helper

def ols(y, X):
    """OLS with intercept included in X. Returns slope of X[:,1], SE, p, CI, plus Spearman."""
    y = np.asarray(y, float)
    X = np.asarray(X, float)
    n, k = X.shape
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = n - k
    s2 = (resid @ resid) / max(dof, 1)
    cov = s2 * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    tstat = beta[1] / se[1]
    p = 2 * stats.t.sf(abs(tstat), dof)
    tcrit = stats.t.ppf(0.975, dof)
    return {
        "slope": beta[1], "se": se[1], "p": p,
        "ci_lo": beta[1] - tcrit * se[1], "ci_hi": beta[1] + tcrit * se[1],
        "rho": stats.spearmanr(X[:, 1], y).statistic,
        "rho_p": stats.spearmanr(X[:, 1], y).pvalue,
    }


def power_verdict(n, sd):
    """M4-floor verdict + nearest MDES from the committed power grid (dispersion 0.2)."""
    pg = pd.read_csv(RESULTS / "m4-human-power.csv")
    mid = pg[pg["biological_dispersion"].round(6) == 0.2]
    tiles = mid.assign(dist=(mid["n_donors"] - n).abs() + (mid["donor_age_sd_years"] - sd).abs() * 3)
    tile = tiles.sort_values("dist").iloc[0]
    floor = (n >= MIN_DONORS_FOR_REGRESSION) and (sd >= MIN_AGE_SD)
    return ("POWERED-FLOOR" if floor else "LOW-POWERED",
            float(tile["min_detectable_30y_log2fc"]), int(tile["n_donors"]), float(tile["donor_age_sd_years"]))


# ---------------------------------------------------------------- main

def main():
    up, dn = load_program()
    symbols = sorted(set(up) | set(dn))
    var_filter = "feature_name in [" + ",".join(f"'{s}'" for s in symbols) + "]"

    census = cellxgene_census.open_soma(census_version="stable")
    ids = ",".join(f"'{k}'" for k in COHORTS)

    # ---- control-donor selection from metadata (reproducible, no /tmp state)
    obs = cellxgene_census.get_obs(
        census, "Homo sapiens",
        value_filter=f"dataset_id in [{ids}]",
        column_names=["dataset_id", "donor_id", "development_stage", "disease",
                      "cell_type", "sex", "is_primary_data"],
    )
    obs["dataset_id"] = obs["dataset_id"].astype(str)
    obs["donor_id"] = obs["donor_id"].astype(str)
    obs["age"] = obs["development_stage"].astype(str).str.extract(AGE_RE).astype(float)
    obs = obs[obs["is_primary_data"] == True]  # noqa: E712

    donor_meta = (obs.groupby(["dataset_id", "donor_id"], observed=True)
                    .agg(age=("age", "first"),
                         disease=("disease", lambda s: str(s.iloc[0])),
                         sex=("sex", lambda s: str(s.iloc[0])),
                         cells=("donor_id", "size"))
                    .reset_index())
    controls = donor_meta[(donor_meta["disease"] == "normal") & donor_meta["age"].notna()].copy()
    prov(f"donors total {donor_meta.shape[0]}; control donors with numeric age: {controls.shape[0]}")

    score_rows = []
    for ds_id, (label, group) in COHORTS.items():
        cds = controls[controls["dataset_id"] == ds_id]
        if len(cds) < 4:
            prov(f"[{label}] only {len(cds)} control donors with age; SKIPPED")
            continue
        donor_list = ",".join(f"'{d}'" for d in cds["donor_id"])
        prov(f"[{label}] downloading {len(cds)} control donors "
             f"(ages {cds['age'].min():.0f}-{cds['age'].max():.0f}) ...")
        adata = cellxgene_census.get_anndata(
            census, "Homo sapiens",
            obs_value_filter=f"dataset_id == '{ds_id}' and is_primary_data == True "
                             f"and donor_id in [{donor_list}]",
            var_value_filter=var_filter,
        )
        X = adata.X
        if sp.issparse(X):
            X = X.tocsr()
        logX = X.copy()
        logX.data = np.log1p(logX.data)
        names = list(adata.var["feature_name"])
        up_i = [names.index(g) for g in up if g in names]
        dn_i = [names.index(g) for g in dn if g in names]
        coverage = len(set(up_i) | set(dn_i)) / (len(up) + len(dn))
        prov(f"[{label}] cells={adata.n_obs:,} genes_measured={adata.n_vars} "
             f"program_coverage={coverage:.3f}")

        up_mean = np.asarray(logX[:, up_i].mean(axis=1)).ravel() if up_i else np.zeros(adata.n_obs)
        dn_mean = np.asarray(logX[:, dn_i].mean(axis=1)).ravel() if dn_i else np.zeros(adata.n_obs)
        cell_score = up_mean - dn_mean

        sdf = pd.DataFrame({
            "donor_id": adata.obs["donor_id"].astype(str).values,
            "cell_type": adata.obs["cell_type"].astype(str).values,
            "score": cell_score,
        })
        sdf = sdf[sdf["cell_type"] != "unknown"]
        agg = (sdf.groupby(["donor_id", "cell_type"], observed=True)
                  .agg(score=("score", "mean"), n_cells=("score", "size"))
                  .reset_index())
        agg = agg[agg["n_cells"] >= MIN_CELLS_PER_DONOR_CT]
        meta = cds.set_index("donor_id")
        agg["age"] = agg["donor_id"].map(meta["age"])
        agg["sex"] = agg["donor_id"].map(meta["sex"])
        agg["cohort"] = label
        agg["group"] = group
        agg["coverage"] = round(coverage, 4)
        score_rows.append(agg)
        prov(f"[{label}] scored donor x cell-type units: {agg.shape[0]}")

        # per-donor whole-cohort score (all cells pooled)
        whole = sdf.groupby("donor_id", observed=True).agg(score=("score", "mean"),
                                                           n_cells=("score", "size")).reset_index()
        whole["age"] = whole["donor_id"].map(meta["age"])
        whole["sex"] = whole["donor_id"].map(meta["sex"])
        _ = whole  # used in regressions below via summary assembly

        del adata, logX, X
        import gc; gc.collect()

    scores = pd.concat(score_rows, ignore_index=True)
    scores.to_csv(RESULTS / "m9-donor-scores.csv", index=False)
    prov(f"total donor x cell-type units: {scores.shape[0]}")

    # ---- whole-donor summary per cohort
    summary_rows = []
    for ds_id, (label, group) in COHORTS.items():
        sub = scores[scores["cohort"] == label]
        if sub.empty:
            continue
        whole = (sub.groupby("donor_id", observed=True)
                    .apply(lambda d: np.average(d["score"], weights=d["n_cells"]))
                    .rename("score").reset_index())
        meta = controls[controls["dataset_id"] == ds_id].set_index("donor_id")
        whole["age"] = whole["donor_id"].map(meta["age"])
        whole["sex"] = whole["donor_id"].map(meta["sex"]).astype(str)
        whole = whole.dropna(subset=["age"])
        if len(whole) < 4:
            continue
        male = (whole["sex"].str.lower() == "male").astype(float).values
        Xd = np.column_stack([np.ones(len(whole)), whole["age"].values, male])
        r = ols(whole["score"].values, Xd)
        verdict, mdes, tn, tsd = power_verdict(len(whole), whole["age"].std())
        summary_rows.append({
            "cohort": label, "group": group, "n_donors": len(whole),
            "age_min": whole["age"].min(), "age_max": whole["age"].max(),
            "age_sd": round(whole["age"].std(), 2),
            "slope_per_year": r["slope"], "se": r["se"],
            "ci_lo": r["ci_lo"], "ci_hi": r["ci_hi"], "p": r["p"],
            "rho": r["rho"], "rho_p": r["rho_p"],
            "norm_slope_30y": r["slope"] * 30 / HUMAN_MAX_LIFESPAN,
            "verdict": verdict, "mdes30y_disp02": mdes,
        })
        fig, ax = plt.subplots(figsize=(5.2, 4))
        ax.scatter(whole["age"], whole["score"], s=18, alpha=0.6, color="#0369a1")
        xs = np.linspace(whole["age"].min(), whole["age"].max(), 50)
        male_c = male.mean()
        ax.plot(xs, r["slope"] * xs + (whole["score"].mean() - r["slope"] * whole["age"].mean()),
                color="#dc2626", lw=2)
        ax.set_xlabel("donor age (years)")
        ax.set_ylabel("mouse aging-program score")
        ax.set_title(f"{label}: n={len(whole)}, slope={r['slope']:.4f}/yr, p={r['p']:.2g}")
        fig.tight_layout()
        fig.savefig(FIGURES / f"m9-{label.lower()}-score-vs-age.png", dpi=150)
        plt.close(fig)
    summary = pd.DataFrame(summary_rows)
    from statsmodels.stats.multitest import multipletests
    summary["q_bh"] = multipletests(summary["p"], method="fdr_bh")[1]
    summary.to_csv(RESULTS / "m9-cohort-summary.csv", index=False)

    # ---- per-cell-type regressions
    group_of = {label: group for _, (label, group) in COHORTS.items()}
    reg_rows = []
    for (label, ct), sub in scores.groupby(["cohort", "cell_type"], observed=True):
        d = sub.groupby("donor_id", observed=True).agg(score=("score", "mean"),
                                                       age=("age", "first"),
                                                       sex=("sex", "first")).dropna(subset=["age"])
        if len(d) < 4:
            continue
        male = (d["sex"].astype(str).str.lower() == "male").astype(float).values
        Xd = np.column_stack([np.ones(len(d)), d["age"].values, male])
        try:
            r = ols(d["score"].values, Xd)
        except np.linalg.LinAlgError:
            continue
        verdict, mdes, tn, tsd = power_verdict(len(d), d["age"].std())
        reg_rows.append({
            "cohort": label, "group": group_of.get(label, ""), "cell_type": ct, "n_donors": len(d),
            "age_sd": round(d["age"].std(), 2),
            "slope_per_year": r["slope"], "se": r["se"], "p": r["p"],
            "ci_lo": r["ci_lo"], "ci_hi": r["ci_hi"], "rho": r["rho"],
            "verdict": verdict, "mdes30y_disp02": mdes,
        })
    regs = pd.DataFrame(reg_rows).sort_values(["group", "cohort", "p"])
    regs.to_csv(RESULTS / "m9-regressions.csv", index=False)

    (RESULTS / "m9-provenance.log").write_text("\n".join(PROV) + "\n")
    print("\n== cohort summary ==")
    print(summary.to_string(index=False))
    print("\ntop cell-type regressions by significance:")
    print(regs.sort_values("p").head(18).to_string(index=False))
    print("\n[m9] done")


if __name__ == "__main__":
    sys.exit(main())
