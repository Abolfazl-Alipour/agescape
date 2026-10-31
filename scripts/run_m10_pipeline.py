#!/usr/bin/env python3
"""Evaluate both aging-signature tiers in the controlled DLPFC lifespan atlas."""

from __future__ import annotations

import hashlib
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rdata
import scipy.sparse as sp
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from src.aging_programs import load_program_tiers

RESULTS = REPO / "results"
FIGURES = REPO / "figures"
EXTERNAL = REPO / "data" / "external" / "dlpfc-atlas"
RESULTS.mkdir(exist_ok=True)
FIGURES.mkdir(exist_ok=True)
RDS = EXTERNAL / "lifespan_pseudobulk.rds"
MD5 = "aef628e7bc2e5134b790d98c6f57ec48"


def ols(y, X):
    y = np.asarray(y, float).ravel()
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


def score_matrix(mat, genes, up, dn, normalize=True):
    gidx = {gene: i for i, gene in enumerate(genes)}
    up_i = [gidx[gene] for gene in up if gene in gidx]
    dn_i = [gidx[gene] for gene in dn if gene in gidx]
    X = mat.astype(np.float64)
    if normalize:
        library = np.asarray(X.sum(axis=0)).ravel()
        X = X @ sp.diags(1e6 / np.where(library > 0, library, 1.0))
    X = X.tocsc()
    X.data = np.log1p(X.data)
    up_mean = np.asarray(X[up_i, :].mean(axis=0)).ravel() if up_i else np.zeros(X.shape[1])
    dn_mean = np.asarray(X[dn_i, :].mean(axis=0)).ravel() if dn_i else np.zeros(X.shape[1])
    coverage = len(set(up_i) | set(dn_i)) / max(len(set(up) | set(dn)), 1)
    return up_mean - dn_mean, coverage, len(up_i), len(dn_i)


def analyze_tier(obj, donor_ids, metadata, tier, program):
    age, sex, brain_bank, pmi, pc1, pc2, pc3 = metadata
    donor_scores = {donor: [] for donor in donor_ids}
    confounder_rows = []
    regression_rows = []
    for cell_type in sorted(obj.assays.data.listData.keys()):
        matrix = obj.assays.data.listData[cell_type]
        genes = [str(gene) for gene in matrix.Dimnames[0]]
        assert [str(donor) for donor in matrix.Dimnames[1]] == donor_ids
        mat = sp.csc_matrix((np.asarray(matrix.x), np.asarray(matrix.i), np.asarray(matrix.p)), shape=tuple(matrix.Dim))
        keep = np.asarray(mat.sum(axis=0)).ravel() > 0
        score, coverage, n_up, n_down = score_matrix(mat, genes, program["up"], program["down"])
        raw_score, _, _, _ = score_matrix(mat, genes, program["up"], program["down"], normalize=False)
        for donor, value, present in zip(donor_ids, score, keep):
            if present:
                donor_scores[donor].append(value)

        y = score[keep]
        a = age[keep]
        male = (sex[keep] == "Male").astype(float)
        hbcc = (brain_bank[keep] == "HBCC").astype(float)
        pmi_k = pmi[keep]
        tech1, tech2, tech3 = pc1[keep], pc2[keep], pc3[keep]
        unadjusted = ols(y, np.column_stack([np.ones(len(y)), a]))
        covar = ols(y, np.column_stack([np.ones(len(y)), a, male, hbcc, pmi_k, tech1, tech2, tech3]))

        def fit_mask(mask, include_bank=True):
            if mask.sum() < 10:
                return None
            columns = [np.ones(mask.sum()), a[mask], male[mask]]
            if include_bank:
                columns.append(hbcc[mask])
            columns.extend([pmi_k[mask], tech1[mask], tech2[mask], tech3[mask]])
            return ols(y[mask], np.column_stack(columns))

        adult = fit_mask(a >= 20)
        age40 = fit_mask(a >= 40)
        developmental = fit_mask(a < 20, include_bank=False)
        hbcc_only = fit_mask(hbcc == 1, include_bank=False)
        mssm_only = fit_mask(hbcc == 0, include_bank=False)

        def values(fit, prefix):
            return {
                f"slope_{prefix}": fit["slope"] if fit else np.nan,
                f"se_{prefix}": fit["se"] if fit else np.nan,
                f"p_{prefix}": fit["p"] if fit else np.nan,
            }

        confounder_rows.append({
            "cell_type": cell_type, "n_all": len(y), **values(unadjusted, "all"),
            "n_covar_adj": len(y), **values(covar, "covar_adj"),
            "n_adult_ge20": int((a >= 20).sum()), **values(adult, "adult_ge20"),
            "n_age_ge40": int((a >= 40).sum()), **values(age40, "age_ge40"),
            "n_dev_lt20": int((a < 20).sum()), **values(developmental, "dev_lt20"),
            "n_hbcc": int((hbcc == 1).sum()), **values(hbcc_only, "hbcc"),
            "n_mssm": int((hbcc == 0).sum()), **values(mssm_only, "mssm"),
            "program_tier": tier, "normalization": "log1p_CPM_pseudobulk",
        })
        raw_unadjusted = ols(raw_score[keep], np.column_stack([np.ones(len(y)), a]))
        regression_rows.append({
            "cell_type": cell_type, "n_donors": int(keep.sum()), "n_zero_donors": int((~keep).sum()),
            "coverage": round(coverage, 4), "n_up_genes": n_up, "n_down_genes": n_down,
            "slope": unadjusted["slope"], "se": unadjusted["se"], "p": unadjusted["p"],
            "ci_lo": unadjusted["ci_lo"], "ci_hi": unadjusted["ci_hi"],
            "rho": unadjusted["rho"], "rho_p": unadjusted["rho_p"],
            "slope_raw_counts": raw_unadjusted["slope"], "p_raw_counts": raw_unadjusted["p"],
            "slope_covar_adj": covar["slope"], "p_covar_adj": covar["p"],
            "n_adults": int((a >= 20).sum()),
            "slope_adult": adult["slope"] if adult else np.nan,
            "p_adult": adult["p"] if adult else np.nan,
            "program_tier": tier, "normalization": "log1p_CPM_pseudobulk",
        })

    confounders = pd.DataFrame(confounder_rows)
    for column in ["p_all", "p_covar_adj", "p_adult_ge20", "p_age_ge40", "p_dev_lt20", "p_hbcc", "p_mssm"]:
        valid = confounders[column].notna()
        qcolumn = column.replace("p_", "q_")
        confounders[qcolumn] = np.nan
        confounders.loc[valid, qcolumn] = multipletests(confounders.loc[valid, column], method="fdr_bh")[1]
    regressions = pd.DataFrame(regression_rows).sort_values("p").reset_index(drop=True)
    regressions["q_bh"] = multipletests(regressions["p"], method="fdr_bh")[1]

    whole_scores = np.array([np.mean(donor_scores[d]) if donor_scores[d] else np.nan for d in donor_ids])
    valid = ~np.isnan(whole_scores)
    y = whole_scores[valid]
    a = age[valid]
    male = (sex[valid] == "Male").astype(float)
    hbcc = (brain_bank[valid] == "HBCC").astype(float)
    unadjusted = ols(y, np.column_stack([np.ones(len(y)), a]))
    covar = ols(y, np.column_stack([
        np.ones(len(y)), a, male, hbcc, pmi[valid], pc1[valid], pc2[valid], pc3[valid],
    ]))
    summary = pd.DataFrame([{
        "cohort": "DLPFC_Lifespan_Atlas", "n_donors": int(valid.sum()),
        "age_min": a.min(), "age_max": a.max(), "age_sd": a.std(),
        "n_cell_types": len(obj.assays.data.listData), "slope": unadjusted["slope"],
        "se": unadjusted["se"], "p": unadjusted["p"], "ci_lo": unadjusted["ci_lo"],
        "ci_hi": unadjusted["ci_hi"], "rho": unadjusted["rho"], "rho_p": unadjusted["rho_p"],
        "slope_covar_adj": covar["slope"], "p_covar_adj": covar["p"],
        "ci_lo_covar_adj": covar["ci_lo"], "ci_hi_covar_adj": covar["ci_hi"],
        "program_tier": tier, "normalization": "log1p_CPM_pseudobulk",
    }])
    return confounders, regressions, summary


def main():
    warnings.filterwarnings("ignore")
    digest = hashlib.md5(RDS.read_bytes()).hexdigest()
    assert digest == MD5, f"RDS md5 mismatch: {digest}"
    obj = rdata.read_rds(str(RDS))
    cd = obj.colData
    donor_ids = [str(donor) for donor in cd.rownames]
    metadata = tuple(np.asarray(cd.listData[name]).astype(dtype).ravel() for name, dtype in (
        ("Age", float), ("Sex", str), ("Brain_bank", str), ("PMI", float),
        ("TechPC1", float), ("TechPC2", float), ("TechPC3", float),
    ))
    programs = load_program_tiers(REPO, "human")
    for tier, suffix in (("P1_Clean", ""), ("P0_Raw257", "-p0")):
        confounders, regressions, summary = analyze_tier(obj, donor_ids, metadata, tier, programs[tier])
        confounders.to_csv(RESULTS / f"m10-confounder-audit{suffix}.csv", index=False)
        regressions.to_csv(RESULTS / f"m10-celltype-regressions{suffix}.csv", index=False)
        summary.to_csv(RESULTS / f"m10-cohort-summary{suffix}.csv", index=False)
        print(tier)
        print(summary.to_string(index=False))
        if tier == "P1_Clean":
            ordered = regressions.sort_values("slope")
            fig, ax = plt.subplots(figsize=(6.5, 7.5))
            ax.errorbar(ordered["slope"], ordered["cell_type"],
                        xerr=[ordered["slope"] - ordered["ci_lo"], ordered["ci_hi"] - ordered["slope"]],
                        fmt="o", color="#0369a1", ecolor="#9ec3dd", ms=4)
            ax.axvline(0, color="#444", lw=0.8, linestyle="--")
            ax.set_xlabel("P1_Clean score slope per year (95% CI)")
            ax.set_title("DLPFC lifespan association by cell type")
            fig.tight_layout()
            fig.savefig(FIGURES / "m10-dlpfc-celltype-forest.png", dpi=150)
            plt.close(fig)
    (RESULTS / "m10-provenance.log").write_text(
        f"source: syn62064518/lifespan_pseudobulk.rds md5 {digest}\n"
        "primary: P1_Clean log1p(CPM); sensitivity: P0_Raw257 and raw-count scores\n"
        "only aggregate regression outputs committed\n"
    )


if __name__ == "__main__":
    main()
