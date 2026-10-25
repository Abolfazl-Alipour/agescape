#!/usr/bin/env python3
"""
run_m10_pipeline.py — Gold-standard human lifespan replication on the DLPFC Lifespan Atlas
(controlled-access, AD Knowledge Portal DUA; source: syn52396927 / lifespan_pseudobulk.rds).

Run with the census environment:
    /tmp/censusenv/bin/python scripts/run_m10_pipeline.py

Rigorous Rebuild Highlights:
1. Primary normalization is CPM: log1p(CPM) within donor x cell type.
2. PsychAD confounder adjustment: ~ Age + Sex + Brain_bank + PMI + TechPC1-3.
3. Developmental stratum (<20y) fits WITHOUT Brain_bank (100% HBCC, singular otherwise).
4. Outputs wide results/m10-confounder-audit.csv (27 rows, 1 per cell type).
5. Strict DUA hygiene: NO per-donor scatter plots or donor files committed.
   Compute scores in RAM; generates only aggregate forest plot.
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats as stats

import rdata
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
warnings.filterwarnings("ignore", category=UserWarning, module="rdata")

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
FIGURES = REPO / "figures"
REF_DIR = REPO / "data" / "reference"
EXTERNAL = REPO / "data" / "external" / "dlpfc-atlas"
for d in (RESULTS, FIGURES):
    d.mkdir(parents=True, exist_ok=True)

RDS = EXTERNAL / "lifespan_pseudobulk.rds"
MD5 = "aef628e7bc2e5134b790d98c6f57ec48"

PROV = []


def prov(msg: str) -> None:
    line = f"[m10] {msg}"
    print(line, flush=True)
    PROV.append(line)


def load_program():
    pg = pd.read_csv(RESULTS / "m6-program-genes.csv")
    up_m = pg["program_up"].dropna().tolist()
    dn_m = pg["program_down"].dropna().tolist()
    orth = pd.read_csv(REF_DIR / "orthologs-mouse-human.csv")
    m2h = dict(zip(orth["mouse_gene_symbol"], orth["human_gene_symbol"]))
    up = sorted({m2h[g] for g in up_m if g in m2h and isinstance(m2h[g], str)})
    dn = sorted({m2h[g] for g in dn_m if g in m2h and isinstance(m2h[g], str)})
    prov(f"program: {len(up_m)} up / {len(dn_m)} down mouse genes -> {len(up)} up / {len(dn)} down human orthologs")
    return up, dn


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
        "dof": dof, "n": n
    }


def score_matrix(mat, genes, up, dn, norm="cpm"):
    gidx = {g: i for i, g in enumerate(genes)}
    up_i = [gidx[g] for g in up if g in gidx]
    dn_i = [gidx[g] for g in dn if g in gidx]
    X = mat.astype(np.float64)
    if norm == "cpm":
        lib = np.asarray(X.sum(axis=0)).ravel()
        lib = np.where(lib > 0, lib, 1.0)
        X = X @ sp.diags(1e6 / lib)
    Xl = X.copy()
    Xl.data = np.log1p(Xl.data)
    up_mean = np.asarray(Xl[up_i, :].mean(axis=0)).ravel() if up_i else np.zeros(X.shape[1])
    dn_mean = np.asarray(Xl[dn_i, :].mean(axis=0)).ravel() if dn_i else np.zeros(X.shape[1])
    coverage = len(set(up_i) | set(dn_i)) / (len(up) + len(dn))
    return up_mean - dn_mean, coverage, len(up_i), len(dn_i)


def main():
    import hashlib
    md5 = hashlib.md5(RDS.read_bytes()).hexdigest()
    assert md5 == MD5, f"RDS md5 mismatch: {md5}"
    prov(f"source: syn62064518/lifespan_pseudobulk.rds md5 verified {md5}")
    prov("DUA: AD Knowledge Portal controlled access; aggregates only committed")

    up, dn = load_program()
    prov("scoring recipe: primary CPM-normalized log1p(CPM); robustness raw counts")

    obj = rdata.read_rds(str(RDS))
    cd = obj.colData
    donor_ids = list(cd.rownames)
    
    # Ravel all colData fields
    age = np.asarray(cd.listData["Age"], float).ravel()
    sex = np.asarray(cd.listData["Sex"]).astype(str).ravel()
    brain_bank = np.asarray(cd.listData["Brain_bank"]).astype(str).ravel()
    pmi = np.asarray(cd.listData["PMI"], float).ravel()
    tech_pc1 = np.asarray(cd.listData["TechPC1"], float).ravel()
    tech_pc2 = np.asarray(cd.listData["TechPC2"], float).ravel()
    tech_pc3 = np.asarray(cd.listData["TechPC3"], float).ravel()
    
    n_donors = len(donor_ids)
    prov(f"cohort: {n_donors} control donors, age {age.min():.1f}-{age.max():.1f} (SD {age.std():.1f})")

    cts = sorted(obj.assays.data.listData.keys())
    prov(f"cell types: {len(cts)}")

    # We will build wide results for confounder audit and regression table
    confounder_rows = []
    regression_rows = []
    donor_scores_cpm = {d: [] for d in donor_ids}

    for ct in cts:
        m = obj.assays.data.listData[ct]
        genes = [str(g) for g in m.Dimnames[0]]
        assert [str(d) for d in m.Dimnames[1]] == donor_ids, "donor order drift"
        mat = sp.csc_matrix((np.asarray(m.x), np.asarray(m.i), np.asarray(m.p)), shape=tuple(m.Dim))

        colsum = np.asarray(mat.sum(axis=0)).ravel()
        keep = colsum > 0
        n_zero = int((~keep).sum())

        # CPM primary, raw counts secondary
        sc_cpm, cov, nu, nd = score_matrix(mat, genes, up, dn, norm="cpm")
        sc_raw, _, _, _ = score_matrix(mat, genes, up, dn, norm=None)

        for d, s, present in zip(donor_ids, sc_cpm, keep):
            if present:
                donor_scores_cpm[d].append(s)

        # Base vectors for present donors
        y = sc_cpm[keep]
        a = age[keep]
        s_male = (sex[keep] == "Male").astype(float)
        bb_hbcc = (brain_bank[keep] == "HBCC").astype(float)
        pmi_k = pmi[keep]
        pc1 = tech_pc1[keep]
        pc2 = tech_pc2[keep]
        pc3 = tech_pc3[keep]

        # 1. Unadjusted All Donors
        X_unadj = np.column_stack([np.ones(len(y)), a])
        r_unadj = ols(y, X_unadj)

        # 2. Covariate Adjusted All Donors (~ Age + Sex + Brain_bank + PMI + TechPC1-3)
        X_covar = np.column_stack([np.ones(len(y)), a, s_male, bb_hbcc, pmi_k, pc1, pc2, pc3])
        r_covar = ols(y, X_covar)

        # 3. Adult (>=20y) Covariate Adjusted
        ad_mask = a >= 20.0
        X_adult = np.column_stack([np.ones(ad_mask.sum()), a[ad_mask], s_male[ad_mask], bb_hbcc[ad_mask], pmi_k[ad_mask], pc1[ad_mask], pc2[ad_mask], pc3[ad_mask]])
        r_adult = ols(y[ad_mask], X_adult) if ad_mask.sum() >= 10 else None

        # 4. Age >= 40y Covariate Adjusted
        ge40_mask = a >= 40.0
        X_ge40 = np.column_stack([np.ones(ge40_mask.sum()), a[ge40_mask], s_male[ge40_mask], bb_hbcc[ge40_mask], pmi_k[ge40_mask], pc1[ge40_mask], pc2[ge40_mask], pc3[ge40_mask]])
        r_ge40 = ols(y[ge40_mask], X_ge40) if ge40_mask.sum() >= 10 else None

        # 5. Developmental (<20y) — WITHOUT Brain_bank because 100% HBCC!
        dev_mask = a < 20.0
        X_dev = np.column_stack([np.ones(dev_mask.sum()), a[dev_mask], s_male[dev_mask], pmi_k[dev_mask], pc1[dev_mask], pc2[dev_mask], pc3[dev_mask]])
        r_dev = ols(y[dev_mask], X_dev) if dev_mask.sum() >= 10 else None

        # 6. HBCC Only
        hbcc_mask = bb_hbcc == 1.0
        X_hbcc = np.column_stack([np.ones(hbcc_mask.sum()), a[hbcc_mask], s_male[hbcc_mask], pmi_k[hbcc_mask], pc1[hbcc_mask], pc2[hbcc_mask], pc3[hbcc_mask]])
        r_hbcc = ols(y[hbcc_mask], X_hbcc) if hbcc_mask.sum() >= 10 else None

        # 7. MSSM Only
        mssm_mask = bb_hbcc == 0.0
        X_mssm = np.column_stack([np.ones(mssm_mask.sum()), a[mssm_mask], s_male[mssm_mask], pmi_k[mssm_mask], pc1[mssm_mask], pc2[mssm_mask], pc3[mssm_mask]])
        r_mssm = ols(y[mssm_mask], X_mssm) if mssm_mask.sum() >= 10 else None

        # Wide confounder audit row
        confounder_rows.append({
            "cell_type": ct,
            "n_all": len(y),
            "slope_all": r_unadj["slope"], "se_all": r_unadj["se"], "p_all": r_unadj["p"],
            "n_covar_adj": len(y),
            "slope_covar_adj": r_covar["slope"], "se_covar_adj": r_covar["se"], "p_covar_adj": r_covar["p"],
            "n_adult_ge20": int(ad_mask.sum()),
            "slope_adult_ge20": r_adult["slope"] if r_adult else np.nan, "se_adult_ge20": r_adult["se"] if r_adult else np.nan, "p_adult_ge20": r_adult["p"] if r_adult else np.nan,
            "n_age_ge40": int(ge40_mask.sum()),
            "slope_age_ge40": r_ge40["slope"] if r_ge40 else np.nan, "se_age_ge40": r_ge40["se"] if r_ge40 else np.nan, "p_age_ge40": r_ge40["p"] if r_ge40 else np.nan,
            "n_dev_lt20": int(dev_mask.sum()),
            "slope_dev_lt20": r_dev["slope"] if r_dev else np.nan, "se_dev_lt20": r_dev["se"] if r_dev else np.nan, "p_dev_lt20": r_dev["p"] if r_dev else np.nan,
            "n_hbcc": int(hbcc_mask.sum()),
            "slope_hbcc": r_hbcc["slope"] if r_hbcc else np.nan, "se_hbcc": r_hbcc["se"] if r_hbcc else np.nan, "p_hbcc": r_hbcc["p"] if r_hbcc else np.nan,
            "n_mssm": int(mssm_mask.sum()),
            "slope_mssm": r_mssm["slope"] if r_mssm else np.nan, "se_mssm": r_mssm["se"] if r_mssm else np.nan, "p_mssm": r_mssm["p"] if r_mssm else np.nan,
        })

        # Regression summary row (CPM primary, raw counts secondary)
        r_raw_unadj = ols(sc_raw[keep], X_unadj)
        r_raw_sex = ols(sc_raw[keep], np.column_stack([np.ones(len(y)), a, s_male]))
        regression_rows.append({
            "cell_type": ct, "n_donors": int(keep.sum()), "n_zero_donors": n_zero,
            "coverage": round(cov, 4), "n_up_genes": nu, "n_down_genes": nd,
            "slope": r_unadj["slope"], "se": r_unadj["se"], "p": r_unadj["p"],  # primary is CPM
            "ci_lo": r_unadj["ci_lo"], "ci_hi": r_unadj["ci_hi"],
            "rho": r_unadj["rho"], "rho_p": r_unadj["rho_p"],
            "slope_raw_counts": r_raw_unadj["slope"], "p_raw_counts": r_raw_unadj["p"],
            "slope_covar_adj": r_covar["slope"], "p_covar_adj": r_covar["p"],
            "n_adults": int(ad_mask.sum()),
            "slope_adult": r_adult["slope"] if r_adult else np.nan,
            "p_adult": r_adult["p"] if r_adult else np.nan,
        })

    # Save wide confounder audit
    confound_df = pd.DataFrame(confounder_rows)
    from statsmodels.stats.multitest import multipletests
    for col in ["p_all", "p_covar_adj", "p_adult_ge20", "p_age_ge40", "p_dev_lt20", "p_hbcc", "p_mssm"]:
        valid = confound_df[col].notna()
        qcol = col.replace("p_", "q_")
        confound_df[qcol] = np.nan
        if valid.sum() > 0:
            confound_df.loc[valid, qcol] = multipletests(confound_df.loc[valid, col], method="fdr_bh")[1]

    confound_path = RESULTS / "m10-confounder-audit.csv"
    confound_df.to_csv(confound_path, index=False)
    prov(f"Saved wide {confound_path} ({len(confound_df)} rows, 1 per cell type)")

    # Save celltype regressions table
    reg_df = pd.DataFrame(regression_rows).sort_values("p").reset_index(drop=True)
    valid_p = reg_df["p"].notna()
    reg_df["q_bh"] = np.nan
    if valid_p.sum() > 0:
        reg_df.loc[valid_p, "q_bh"] = multipletests(reg_df.loc[valid_p, "p"], method="fdr_bh")[1]
    reg_df.to_csv(RESULTS / "m10-celltype-regressions.csv", index=False)
    prov(f"Saved {RESULTS / 'm10-celltype-regressions.csv'}")

    # Whole cortex summary (computed in RAM, no per-donor file saved)
    whole_scores = [np.mean(donor_scores_cpm[d]) if len(donor_scores_cpm[d]) > 0 else np.nan for d in donor_ids]
    valid_donors = ~np.isnan(whole_scores)
    y_wc = np.array(whole_scores)[valid_donors]
    a_wc = age[valid_donors]
    s_wc = (sex[valid_donors] == "Male").astype(float)
    bb_wc = (brain_bank[valid_donors] == "HBCC").astype(float)
    pmi_wc = pmi[valid_donors]
    pc1_wc = tech_pc1[valid_donors]
    pc2_wc = tech_pc2[valid_donors]
    pc3_wc = tech_pc3[valid_donors]

    r_wc_unadj = ols(y_wc, np.column_stack([np.ones(len(y_wc)), a_wc]))
    r_wc_covar = ols(y_wc, np.column_stack([np.ones(len(y_wc)), a_wc, s_wc, bb_wc, pmi_wc, pc1_wc, pc2_wc, pc3_wc]))

    summary_df = pd.DataFrame([{
        "cohort": "DLPFC_Lifespan_Atlas", "n_donors": int(valid_donors.sum()),
        "age_min": a_wc.min(), "age_max": a_wc.max(), "age_sd": a_wc.std(),
        "n_cell_types": len(cts),
        "slope": r_wc_unadj["slope"], "se": r_wc_unadj["se"], "p": r_wc_unadj["p"],
        "ci_lo": r_wc_unadj["ci_lo"], "ci_hi": r_wc_unadj["ci_hi"],
        "rho": r_wc_unadj["rho"], "rho_p": r_wc_unadj["rho_p"],
        "slope_covar_adj": r_wc_covar["slope"], "p_covar_adj": r_wc_covar["p"],
        "ci_lo_covar_adj": r_wc_covar["ci_lo"], "ci_hi_covar_adj": r_wc_covar["ci_hi"],
    }])
    summary_df.to_csv(RESULTS / "m10-cohort-summary.csv", index=False)
    prov(f"Saved {RESULTS / 'm10-cohort-summary.csv'}")

    # Forest plot of cell types
    fig, ax = plt.subplots(figsize=(6.5, 7.5))
    d = reg_df.sort_values("slope")
    ax.errorbar(d["slope"], d["cell_type"],
                xerr=[d["slope"] - d["ci_lo"], d["ci_hi"] - d["slope"]],
                fmt="o", color="#0369a1", ecolor="#9ec3dd", ms=4)
    ax.axvline(0, color="#444", lw=0.8, linestyle="--")
    ax.set_xlabel("P1 Aging Signature Score Slope per Year (95% CI)")
    ax.set_title("Lifespan Aging Signature vs Human Age by DLPFC Cell Type")
    fig.tight_layout()
    fig.savefig(FIGURES / "m10-dlpfc-celltype-forest.png", dpi=150)
    plt.close(fig)
    prov(f"Saved {FIGURES / 'm10-dlpfc-celltype-forest.png'}")

    (RESULTS / "m10-provenance.log").write_text("\n".join(PROV) + "\n")
    print(f"\n[m10] Whole cortex CPM slope: {r_wc_unadj['slope']:+.5f}/yr (p={r_wc_unadj['p']:.4g}); covar-adj p={r_wc_covar['p']:.4g}")
    print("[m10] Done successfully.")


if __name__ == "__main__":
    sys.exit(main())
