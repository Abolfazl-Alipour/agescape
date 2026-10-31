#!/usr/bin/env python3
"""Score public CELLxGENE cohorts with both mouse-derived aging-signature tiers."""

from __future__ import annotations

import argparse
import gc
import re
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from src.aging_programs import load_program_tiers

RESULTS = REPO / "results"
FIGURES = REPO / "figures"
RESULTS.mkdir(exist_ok=True)
FIGURES.mkdir(exist_ok=True)
PROV = []
HUMAN_MAX_LIFESPAN = 122.5
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
MIN_CELLS_PER_DONOR_CT = 20
MIN_DONORS_FOR_REGRESSION = 6
MIN_AGE_SD = 12.0
AGE_RE = re.compile(r"(\d+)-year-old")


def prov(msg: str) -> None:
    line = f"[m9] {msg}"
    print(line, flush=True)
    PROV.append(line)


def ols(y, X):
    y = np.asarray(y, float)
    X = np.asarray(X, float)
    n, k = X.shape
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = n - k
    s2 = (resid @ resid) / max(dof, 1)
    cov = s2 * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.maximum(np.diag(cov), 1e-300))
    tstat = beta[1] / se[1]
    p = 2 * stats.t.sf(abs(tstat), max(dof, 1))
    tcrit = stats.t.ppf(0.975, max(dof, 1))
    return {
        "slope": beta[1], "se": se[1], "p": p,
        "ci_lo": beta[1] - tcrit * se[1], "ci_hi": beta[1] + tcrit * se[1],
        "rho": stats.spearmanr(X[:, 1], y).statistic,
        "rho_p": stats.spearmanr(X[:, 1], y).pvalue,
    }


def power_verdict(n, sd):
    pg = pd.read_csv(RESULTS / "m4-human-power.csv")
    mid = pg[pg["biological_dispersion"].round(6) == 0.2]
    tiles = mid.assign(dist=(mid["n_donors"] - n).abs() + (mid["donor_age_sd_years"] - sd).abs() * 3)
    tile = tiles.sort_values("dist").iloc[0]
    floor = (n >= MIN_DONORS_FOR_REGRESSION) and (sd >= MIN_AGE_SD)
    return (
        "POWERED-FLOOR" if floor else "LOW-POWERED",
        float(tile["min_detectable_30y_log2fc"]),
        int(tile["n_donors"]),
        float(tile["donor_age_sd_years"]),
    )


def score_cells(matrix, names, up, dn, library_sizes=None):
    X = matrix.tocsr() if sp.issparse(matrix) else sp.csr_matrix(matrix)
    X = X.astype(np.float64)
    if library_sizes is not None and X.shape[0] != len(library_sizes) and X.shape[1] == len(library_sizes):
        X = X.T.tocsr()
    if library_sizes is not None and X.shape[0] != len(library_sizes):
        raise ValueError("Expression matrix and library-size dimensions do not match")
    if library_sizes is not None:
        scale = 1e6 / np.where(library_sizes > 0, library_sizes, 1.0)
        X = sp.diags(scale) @ X
    X = X.tocsr()
    X.data = np.log1p(X.data)
    up_i = [i for i, gene in enumerate(names) if gene in set(up)]
    dn_i = [i for i, gene in enumerate(names) if gene in set(dn)]
    up_mean = np.asarray(X[:, up_i].mean(axis=1)).ravel() if up_i else np.zeros(X.shape[0])
    dn_mean = np.asarray(X[:, dn_i].mean(axis=1)).ravel() if dn_i else np.zeros(X.shape[0])
    coverage = len(set(up_i) | set(dn_i)) / max(len(set(up) | set(dn)), 1)
    return up_mean - dn_mean, coverage


def extract_scores(programs):
    import cellxgene_census

    prov("Opening CELLxGENE Census version 2025-11-08")
    census = cellxgene_census.open_soma(census_version="2025-11-08")
    ids = ",".join(f"'{key}'" for key in COHORTS)
    obs = cellxgene_census.get_obs(
        census, "Homo sapiens",
        value_filter=f"dataset_id in [{ids}]",
        column_names=["dataset_id", "donor_id", "development_stage", "disease", "cell_type",
                      "cell_type_ontology_term_id", "raw_sum", "sex", "is_primary_data"],
    )
    obs["dataset_id"] = obs["dataset_id"].astype(str)
    obs["donor_id"] = obs["donor_id"].astype(str)
    obs["age"] = obs["development_stage"].astype(str).str.extract(AGE_RE).astype(float)
    obs = obs[(obs["is_primary_data"] == True) & (obs["disease"] == "normal") & obs["age"].notna()].copy()

    donor_meta = (obs.groupby(["dataset_id", "donor_id"], observed=True)
                  .agg(age=("age", "first"), sex=("sex", lambda s: str(s.iloc[0])),
                       disease=("disease", lambda s: str(s.iloc[0])), cells=("donor_id", "size"))
                  .reset_index())
    prov(f"Normal control donors with numeric age: {len(donor_meta)}")

    prethreshold = (obs[obs["cell_type"].astype(str) != "unknown"]
                    .groupby(["dataset_id", "donor_id", "cell_type"], observed=True)
                    .agg(n_cells=("donor_id", "size"), raw_sum=("raw_sum", "sum"),
                         age=("age", "first"), sex=("sex", lambda s: str(s.iloc[0])))
                    .reset_index())
    prethreshold["cohort"] = prethreshold["dataset_id"].map(lambda x: COHORTS[x][0])
    prethreshold["group"] = prethreshold["dataset_id"].map(lambda x: COHORTS[x][1])
    prethreshold["census_version"] = "2025-11-08"
    prethreshold = prethreshold[["cohort", "group", "dataset_id", "donor_id", "cell_type", "n_cells",
                                 "raw_sum", "age", "sex", "census_version"]]
    count_path = RESULTS / "m9-donor-cell-counts-prethreshold.csv"
    prethreshold.to_csv(count_path, index=False)
    prov(f"Saved {len(prethreshold)} unfiltered donor-cell-type count rows to results/{count_path.name}")

    symbols = sorted({gene for program in programs.values() for values in program.values() for gene in values})
    var_filter = "feature_name in [" + ",".join(f"'{gene}'" for gene in symbols) + "]"
    score_rows = []
    for ds_id, (label, group) in COHORTS.items():
        cds = donor_meta[donor_meta["dataset_id"] == ds_id]
        if len(cds) < 4:
            prov(f"[{label}] only {len(cds)} control donors with age; skipped")
            continue
        donor_list = ",".join(f"'{donor}'" for donor in cds["donor_id"])
        prov(f"[{label}] downloading {len(cds)} control donors")
        adata = cellxgene_census.get_anndata(
            census, "Homo sapiens",
            obs_value_filter=f"dataset_id == '{ds_id}' and is_primary_data == True and donor_id in [{donor_list}]",
            var_value_filter=var_filter,
            obs_column_names=["donor_id", "cell_type", "raw_sum"],
            var_column_names=["feature_name"],
        )
        names = list(adata.var["feature_name"])
        raw_sum = np.asarray(adata.obs["raw_sum"], float)
        cell_scores = {}
        coverages = {}
        for tier, program in programs.items():
            cell_scores[tier], coverages[tier] = score_cells(adata.X, names, program["up"], program["down"], raw_sum)
            prov(f"[{label}] {tier}: {coverages[tier]:.3f} coverage")
        raw_score, _ = score_cells(adata.X, names, programs["P1_Clean"]["up"], programs["P1_Clean"]["down"])
        sdf = pd.DataFrame({
            "donor_id": adata.obs["donor_id"].astype(str).values,
            "cell_type": adata.obs["cell_type"].astype(str).values,
            "score": cell_scores["P1_Clean"],
            "score_p0": cell_scores["P0_Raw257"],
            "score_raw_counts": raw_score,
        })
        sdf = sdf[sdf["cell_type"] != "unknown"]
        agg = (sdf.groupby(["donor_id", "cell_type"], observed=True)
               .agg(score=("score", "mean"), score_p0=("score_p0", "mean"),
                    score_raw_counts=("score_raw_counts", "mean"), n_cells=("score", "size"))
               .reset_index())
        agg = agg[agg["n_cells"] >= MIN_CELLS_PER_DONOR_CT].copy()
        meta = cds.set_index("donor_id")
        agg["age"] = agg["donor_id"].map(meta["age"])
        agg["sex"] = agg["donor_id"].map(meta["sex"])
        agg["cohort"] = label
        agg["group"] = group
        agg["coverage"] = round(coverages["P1_Clean"], 4)
        agg["coverage_p0"] = round(coverages["P0_Raw257"], 4)
        score_rows.append(agg)
        del adata, sdf, agg
        gc.collect()
    scores = pd.concat(score_rows, ignore_index=True)
    scores.to_csv(RESULTS / "m9-donor-scores.csv", index=False)
    prov(f"Total thresholded donor-cell-type score rows: {len(scores)}")
    return scores


def fit_scores(scores, score_column, tier, normalization):
    summary_rows = []
    reg_rows = []
    for label, group in COHORTS.values():
        sub = scores[scores["cohort"] == label].copy()
        if sub.empty:
            continue
        whole = (sub.groupby("donor_id", observed=True)
                 .apply(lambda d: np.average(d[score_column], weights=d["n_cells"]), include_groups=False)
                 .rename("score").reset_index())
        meta = sub.groupby("donor_id", observed=True).agg(age=("age", "first"), sex=("sex", "first"))
        whole["age"] = whole["donor_id"].map(meta["age"])
        whole["sex"] = whole["donor_id"].map(meta["sex"]).astype(str)
        whole = whole.dropna(subset=["age"])
        if len(whole) >= 4:
            male = (whole["sex"].str.lower() == "male").astype(float).values
            columns = [np.ones(len(whole)), whole["age"].values]
            if len(np.unique(male)) > 1:
                columns.append(male)
            result = ols(whole["score"].values, np.column_stack(columns))
            verdict, mdes, _, _ = power_verdict(len(whole), whole["age"].std())
            summary_rows.append({
                "cohort": label, "group": group, "n_donors": len(whole),
                "age_min": whole["age"].min(), "age_max": whole["age"].max(),
                "age_sd": round(whole["age"].std(), 2), "slope_per_year": result["slope"],
                "se": result["se"], "ci_lo": result["ci_lo"], "ci_hi": result["ci_hi"],
                "p": result["p"], "rho": result["rho"], "rho_p": result["rho_p"],
                "norm_slope_30y": result["slope"] * 30 / HUMAN_MAX_LIFESPAN,
                "verdict": verdict, "mdes30y_disp02": mdes,
                "program_tier": tier, "normalization": normalization,
            })
            if normalization == "log1p_CPM_per_cell":
                fig, ax = plt.subplots(figsize=(5.2, 4))
                ax.scatter(whole["age"], whole["score"], s=18, alpha=0.6, color="#0369a1")
                xs = np.linspace(whole["age"].min(), whole["age"].max(), 50)
                intercept = whole["score"].mean() - result["slope"] * whole["age"].mean()
                ax.plot(xs, result["slope"] * xs + intercept, color="#dc2626", lw=2)
                ax.set_xlabel("Donor age (years)")
                ax.set_ylabel(f"{tier} aging signature score")
                ax.set_title(f"{label}: n={len(whole)}, slope={result['slope']:.4f}/yr, p={result['p']:.2g}")
                fig.tight_layout()
                fig.savefig(FIGURES / f"m9-{label.lower()}-score-vs-age.png", dpi=150)
                plt.close(fig)

    for (label, cell_type), sub in scores.groupby(["cohort", "cell_type"], observed=True):
        donors = (sub.groupby("donor_id", observed=True)
                  .agg(score=(score_column, "mean"), age=("age", "first"), sex=("sex", "first"))
                  .dropna(subset=["age"]))
        if len(donors) < 4:
            continue
        male = (donors["sex"].astype(str).str.lower() == "male").astype(float).values
        columns = [np.ones(len(donors)), donors["age"].values]
        if len(np.unique(male)) > 1:
            columns.append(male)
        try:
            result = ols(donors["score"].values, np.column_stack(columns))
        except np.linalg.LinAlgError:
            continue
        verdict, mdes, _, _ = power_verdict(len(donors), donors["age"].std())
        reg_rows.append({
            "cohort": label, "group": dict(COHORTS.values()).get(label, ""), "cell_type": cell_type,
            "n_donors": len(donors), "age_sd": round(donors["age"].std(), 2),
            "slope_per_year": result["slope"], "se": result["se"], "p": result["p"],
            "ci_lo": result["ci_lo"], "ci_hi": result["ci_hi"], "rho": result["rho"],
            "verdict": verdict, "mdes30y_disp02": mdes,
            "program_tier": tier, "normalization": normalization,
        })

    summary = pd.DataFrame(summary_rows)
    summary["q_bh"] = multipletests(summary["p"], method="fdr_bh")[1]
    regs = pd.DataFrame(reg_rows).sort_values(["group", "cohort", "p"]).reset_index(drop=True)
    regs["q_bh"] = multipletests(regs["p"], method="fdr_bh")[1]
    return summary, regs


def run_gene_level_concordance():
    mapping_path = REPO / "data" / "reference" / "m9-celltype-mapping.csv"
    if not mapping_path.exists():
        return
    rows = []
    for _, row in pd.read_csv(mapping_path).iterrows():
        if row["disposition"] != "Matched":
            continue
        slug = f"{row['mouse_tissue']}__{row['mouse_cell_type']}".replace(" ", "_").replace("-", "_").replace(",", "")
        de_path = RESULTS / f"m2_de_{slug}.csv"
        if not de_path.exists():
            continue
        de = pd.read_csv(de_path)
        rows.append({
            "mouse_tissue": row["mouse_tissue"], "mouse_cell_type": row["mouse_cell_type"],
            "human_mapped_cell_type": row["human_mapped_cell_type"],
            "mouse_de_genes_fdr05": int((de["FDR"] < 0.05).sum()) if "FDR" in de else 0,
            "status": "Homologous unit matched",
        })
    pd.DataFrame(rows).to_csv(RESULTS / "m9-gene-level-concordance.csv", index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--use-cache", action="store_true")
    args = parser.parse_args()
    programs = load_program_tiers(REPO, "human")
    for tier, program in programs.items():
        prov(f"{tier}: {len(program['up'])} up and {len(program['down'])} down human orthologs")
    scores_path = RESULTS / "m9-donor-scores.csv"
    required = {"score", "score_p0", "score_raw_counts"}
    if args.use_cache and scores_path.exists() and required.issubset(pd.read_csv(scores_path, nrows=0).columns):
        scores = pd.read_csv(scores_path)
        prov(f"Loaded cached scores from results/{scores_path.name}")
    else:
        scores = extract_scores(programs)

    outputs = [
        ("score", "P1_Clean", "log1p_CPM_per_cell", "m9-cohort-summary.csv", "m9-regressions.csv"),
        ("score_p0", "P0_Raw257", "log1p_CPM_per_cell", "m9-cohort-summary-p0.csv", "m9-regressions-p0.csv"),
        ("score_raw_counts", "P1_Clean", "log1p_raw_counts", "m9-cohort-summary-raw-counts.csv", "m9-regressions-raw-counts.csv"),
    ]
    for score_column, tier, normalization, summary_name, reg_name in outputs:
        summary, regs = fit_scores(scores, score_column, tier, normalization)
        summary.to_csv(RESULTS / summary_name, index=False)
        regs.to_csv(RESULTS / reg_name, index=False)
        prov(f"Saved {summary_name} and {reg_name}")
    run_gene_level_concordance()
    (RESULTS / "m9-provenance.log").write_text("\n".join(PROV) + "\n")


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    sys.exit(main())
