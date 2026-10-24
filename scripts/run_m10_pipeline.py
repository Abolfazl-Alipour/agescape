#!/usr/bin/env python3
"""
run_m10_pipeline.py — Gold-standard human lifespan replication of the mouse
aging program on the DLPFC Lifespan Atlas (controlled-access, AD Knowledge
Portal DUA; source: Synapse syn52396927, folder syn62064535,
lifespan_pseudobulk.rds + SingleCellExperiment colData).

Run with the census environment (has rdata; no scanpy needed):
    /tmp/censusenv/bin/python scripts/run_m10_pipeline.py

M10 question: does the 261-gene mouse aging program (167 up / 94 down,
results/m6-program-genes.csv) track chronological age across the FULL human
lifespan (0.2 to 97 years, 284 control donors, 27 DLPFC cell types,
donor x cell-type pseudobulk already summed by the authors)?

This is the cohort M9 could only approximate: M9's public CELLxGENE cohorts
were 4 to 213 donors each with narrow age spreads; this atlas is 284 donors
with age SD 24 years, which clears the pre-registered M4 human floor
(>=6 donors, donor-age SD >=12) by an order of magnitude.

The score uses the same directional contrast and ortholog set as M9, but the
input summaries differ. M9 averages per-cell log1p counts; M10 applies log1p
to donor-level summed counts. Cross-cohort slope magnitudes are therefore not
directly comparable: score = mean(log1p(count) over up-program genes present)
minus mean(log1p(count) over down-program genes present).

Robustness (eval discipline): (a) CPM-normalized scores (log1p of counts per
million within donor x cell type) re-run for every regression; (b) sex-
adjusted OLS (score ~ age + sex) alongside the primary unadjusted model that
matches M9; (c) Spearman rank correlation alongside OLS; (d) donors with a
zero-count cell-type column (no cells of that type sampled) are excluded from
that cell type's regression and counted.

DUA hygiene: per-donor scores are controlled data. Only aggregate regression
output is committed to results/; per-donor scores go to data/interim/
(gitignored). Nothing per-donor ever enters the repo, the app, or the report.

Outputs:
  results/m10-celltype-regressions.csv   one row per cell type: OLS + Spearman,
                                          primary and CPM-normalized, sex-adj
  results/m10-cohort-summary.csv         whole-cortex headline (per-donor mean score)
  results/m10-provenance.log             cohort composition + recipe decisions
  figures/m10-dlpfc-whole-score-vs-age.png
  figures/m10-dlpfc-celltype-forest.png
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
INTERIM = REPO / "data" / "interim"
EXTERNAL = REPO / "data" / "external" / "dlpfc-atlas"
for d in (RESULTS, FIGURES, INTERIM):
    d.mkdir(parents=True, exist_ok=True)

RDS = EXTERNAL / "lifespan_pseudobulk.rds"
MD5 = "aef628e7bc2e5134b790d98c6f57ec48"  # synapseclient-verified at download

PROV = []


def prov(msg: str) -> None:
    line = f"[m10] {msg}"
    print(line, flush=True)
    PROV.append(line)


# ---------------------------------------------------------------- program genes

def load_program():
    """Identical recipe to M9: m6 program genes -> human orthologs."""
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


# ---------------------------------------------------------------- stats helpers

def ols(y, X):
    """OLS with intercept in X. Slope of X[:,1], SE, p, 95% CI, plus Spearman.
    (Copied from run_m9_pipeline.py so M10 uses byte-identical inference.)"""
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


def score_matrix(mat, genes, up, dn, norm=None):
    """M9 score per donor column: mean(log1p up) - mean(log1p down).

    norm=None: primary pseudobulk recipe (log1p of summed counts).
    norm='cpm': robustness arm, counts rescaled to per-million within the
    donor x cell-type column before log1p.
    Returns (score, coverage, n_up_used, n_dn_used)."""
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


# ---------------------------------------------------------------- main pipeline

def main():
    import hashlib
    md5 = hashlib.md5(RDS.read_bytes()).hexdigest()
    assert md5 == MD5, f"RDS md5 mismatch: {md5}"
    prov(f"source: syn62064518/lifespan_pseudobulk.rds md5 verified {md5}")
    prov("DUA: AD Knowledge Portal controlled access; aggregates only committed")

    up, dn = load_program()
    prov("scoring recipe: donor-pseudobulk up_mean(log1p) - dn_mean(log1p); "
         "robustness arm log1p(CPM); input summaries differ from M9")

    obj = rdata.read_rds(str(RDS))
    cd = obj.colData
    donor_ids = list(cd.rownames)
    age = np.asarray(cd.listData["Age"], float)
    sex = np.asarray(cd.listData["Sex"]).astype(str)
    n_donors = len(donor_ids)
    prov(f"cohort: {n_donors} control donors, age {age.min():.1f}-{age.max():.1f} "
         f"(SD {age.std():.1f}), sexes: {dict(pd.Series(sex).value_counts())}")

    cts = sorted(obj.assays.data.listData.keys())
    prov(f"cell types: {len(cts)}")

    rows = []
    donor_scores = {d: [] for d in donor_ids}
    for ct in cts:
        m = obj.assays.data.listData[ct]
        genes = [str(g) for g in m.Dimnames[0]]
        assert [str(d) for d in m.Dimnames[1]] == donor_ids, "donor order drift"
        mat = sp.csc_matrix((np.asarray(m.x), np.asarray(m.i), np.asarray(m.p)),
                            shape=tuple(m.Dim))

        colsum = np.asarray(mat.sum(axis=0)).ravel()
        keep = colsum > 0
        n_zero = int((~keep).sum())
        if n_zero:
            prov(f"[{ct}] excluded {n_zero} donors with zero cells of this type")

        sc_raw, cov, nu, nd = score_matrix(mat, genes, up, dn)
        sc_cpm, _, _, _ = score_matrix(mat, genes, up, dn, norm="cpm")

        for d, s, present in zip(donor_ids, sc_raw, keep):
            if present:
                donor_scores[d].append(s)

        d_age, d_sex = age[keep], sex[keep]
        r_raw = ols(sc_raw[keep], np.column_stack([np.ones(keep.sum()), d_age]))
        r_cpm = ols(sc_cpm[keep], np.column_stack([np.ones(keep.sum()), d_age]))
        Xs = np.column_stack([np.ones(keep.sum()), d_age, (d_sex == "Female").astype(float)])
        r_sex = ols(sc_raw[keep], Xs)
        ad = d_age >= 18
        r_adult = ols(sc_raw[keep][ad],
                      np.column_stack([np.ones(ad.sum()), d_age[ad]])) if ad.sum() > 10 else None

        rows.append({
            "cell_type": ct, "n_donors": int(keep.sum()), "n_zero_donors": n_zero,
            "coverage": round(cov, 4), "n_up_genes": nu, "n_down_genes": nd,
            "slope": r_raw["slope"], "se": r_raw["se"], "p": r_raw["p"],
            "ci_lo": r_raw["ci_lo"], "ci_hi": r_raw["ci_hi"],
            "rho": r_raw["rho"], "rho_p": r_raw["rho_p"],
            "slope_cpm": r_cpm["slope"], "p_cpm": r_cpm["p"],
            "ci_lo_cpm": r_cpm["ci_lo"], "ci_hi_cpm": r_cpm["ci_hi"],
            "slope_sex_adj": r_sex["slope"], "p_sex_adj": r_sex["p"],
            "ci_lo_sex_adj": r_sex["ci_lo"], "ci_hi_sex_adj": r_sex["ci_hi"],
            "n_adults": int(ad.sum()),
            "slope_adult": (r_adult["slope"] if r_adult else np.nan),
            "p_adult": (r_adult["p"] if r_adult else np.nan),
            "ci_lo_adult": (r_adult["ci_lo"] if r_adult else np.nan),
            "ci_hi_adult": (r_adult["ci_hi"] if r_adult else np.nan),
        })
        adult_p = f"{r_adult['p']:.4g}" if r_adult else "nan"
        prov(f"[{ct}] slope {r_raw['slope']:+.5f}/yr "
             f"CI [{r_raw['ci_lo']:+.5f},{r_raw['ci_hi']:+.5f}] p={r_raw['p']:.4g} "
             f"n={keep.sum()} (cpm p={r_cpm['p']:.4g}, sex-adj p={r_sex['p']:.4g}, "
             f"adult p={adult_p})")

    ct_df = pd.DataFrame(rows).sort_values("p").reset_index(drop=True)
    # Benjamini-Hochberg FDR across the 27 cell-type tests, both models
    from statsmodels.stats.multitest import multipletests
    ct_df["q_bh"] = multipletests(ct_df["p"], method="fdr_bh")[1]
    ct_df["q_bh_cpm"] = multipletests(ct_df["p_cpm"], method="fdr_bh")[1]
    ct_df.to_csv(RESULTS / "m10-celltype-regressions.csv", index=False)

    # Average only cell types observed in each donor; absent cell types are missing,
    # not biological zeros. The number contributing therefore varies by donor.
    whole = pd.DataFrame({
        "donor_id": donor_ids, "age": age, "sex": sex,
        "n_cell_types_observed": [len(donor_scores[d]) for d in donor_ids],
        "score": [np.mean(donor_scores[d]) for d in donor_ids],
    })
    prov(f"whole cortex: {whole['n_cell_types_observed'].min()} to "
         f"{whole['n_cell_types_observed'].max()} observed cell types averaged per donor")
    whole.to_csv(INTERIM / "m10-donor-whole-scores.csv", index=False)  # gitignored (DUA)
    r_all = ols(whole["score"].values,
                np.column_stack([np.ones(n_donors), whole["age"].values]))
    Xs_all = np.column_stack([np.ones(n_donors), whole["age"].values,
                              (whole["sex"] == "Female").astype(float)])
    r_all_sex = ols(whole["score"].values, Xs_all)
    summary = pd.DataFrame([{
        "cohort": "DLPFC_Lifespan_Atlas", "n_donors": n_donors,
        "age_min": age.min(), "age_max": age.max(), "age_sd": age.std(),
        "n_cell_types": len(cts),
        "slope": r_all["slope"], "se": r_all["se"], "p": r_all["p"],
        "ci_lo": r_all["ci_lo"], "ci_hi": r_all["ci_hi"],
        "rho": r_all["rho"], "rho_p": r_all["rho_p"],
        "slope_sex_adj": r_all_sex["slope"], "p_sex_adj": r_all_sex["p"],
        "ci_lo_sex_adj": r_all_sex["ci_lo"], "ci_hi_sex_adj": r_all_sex["ci_hi"],
    }])
    summary.to_csv(RESULTS / "m10-cohort-summary.csv", index=False)
    prov(f"WHOLE CORTEX slope {r_all['slope']:+.5f}/yr "
         f"CI [{r_all['ci_lo']:+.5f},{r_all['ci_hi']:+.5f}] p={r_all['p']:.4g} "
         f"(sex-adj p={r_all_sex['p']:.4g}, rho {r_all['rho']:+.3f})")

    # ---- figures
    fig, ax = plt.subplots(figsize=(6, 4.2))
    colors = np.where(whole["sex"] == "Female", "#7e57c2", "#0369a1")
    ax.scatter(whole["age"], whole["score"], s=16, alpha=0.55, c=colors)
    xs = np.linspace(age.min(), age.max(), 100)
    ax.plot(xs, r_all["slope"] * xs + (whole["score"].mean() - r_all["slope"] * age.mean()),
            color="#c23b22", lw=2)
    for lo, hi in [(r_all["ci_lo"], r_all["ci_hi"])]:
        pass
    ax.set_xlabel("donor age (years)")
    ax.set_ylabel("mouse aging-program score")
    ax.set_title(f"DLPFC Lifespan Atlas (n={n_donors}): slope {r_all['slope']:+.4f}/yr, "
                 f"p={r_all['p']:.4g}")
    fig.tight_layout()
    fig.savefig(FIGURES / "m10-dlpfc-whole-score-vs-age.png", dpi=150)
    plt.close(fig)

    f2, a2 = plt.subplots(figsize=(6, 7))
    d = ct_df.sort_values("slope")
    a2.errorbar(d["slope"], d["cell_type"],
                xerr=[d["slope"] - d["ci_lo"], d["ci_hi"] - d["slope"]],
                fmt="o", color="#0369a1", ecolor="#9ec3dd", ms=4)
    a2.axvline(0, color="#444", lw=0.8)
    a2.set_xlabel("program score slope per year (95% CI)")
    a2.set_title("Mouse aging program vs human age, by DLPFC cell type")
    f2.tight_layout()
    f2.savefig(FIGURES / "m10-dlpfc-celltype-forest.png", dpi=150)
    plt.close(f2)

    (RESULTS / "m10-provenance.log").write_text("\n".join(PROV) + "\n")
    print("\nsignificant cell types (uncorrected p<0.05, primary model):",
          (ct_df["p"] < 0.05).sum(), "of", len(ct_df))
    print(ct_df.head(8)[["cell_type", "slope", "p", "slope_cpm", "p_cpm",
                         "slope_sex_adj", "p_sex_adj"]].to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())
