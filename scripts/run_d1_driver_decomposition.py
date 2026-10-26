#!/usr/bin/env python3
"""
scripts/run_d1_driver_decomposition.py

D1: Driver-Gene Decomposition of Glial Transport Signal in Human DLPFC (PsychAD, N=284).
1. Computes per-gene contributions to the signature's age slope in Microglia and OPCs.
2. Recomputes score age associations after removing a pre-declared generic inflammatory/stress subset.
3. Evaluates P1_Clean tier in Microglia, OPCs, and Whole Cortex.

Outputs:
  - results/d1-glial-driver-decomposition.csv
  - results/d1-glial-model-summary.csv
"""
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats as stats
import rdata

REPO = Path(__file__).resolve().parent.parent

INFLAMMATORY_STRESS_GENES = {
    "CCL2", "ISG15", "ITGAM", "STAT3", "NFKBIA", "NFKBIZ", "SOCS2", "SOCS3",
    "CD74", "LGALS3", "SPP1", "TNFRSF12A", "CDKN1A", "GADD45A", "GADD45B",
    "IER3", "FOS", "JUN", "JUNB", "JUND", "CXCL10", "CCL11"
}

def load_programs():
    pg = pd.read_csv(REPO / "results" / "m6-program-genes.csv")
    up_m = pg["program_up"].dropna().tolist()
    dn_m = pg["program_down"].dropna().tolist()
    orth = pd.read_csv(REPO / "data" / "reference" / "orthologs-mouse-human.csv")
    m2h = dict(zip(orth["mouse_gene_symbol"], orth["human_gene_symbol"]))
    up_p0 = sorted({m2h[g] for g in up_m if g in m2h and isinstance(m2h[g], str)})
    dn_p0 = sorted({m2h[g] for g in dn_m if g in m2h and isinstance(m2h[g], str)})

    hier = pd.read_csv(REPO / "results" / "m5-program-hierarchy.csv")
    p1 = hier[hier["tier"] == "P1_Clean"]
    up_p1_m = p1[p1["direction"] == "UP"]["gene"].tolist()
    dn_p1_m = p1[p1["direction"] == "DOWN"]["gene"].tolist()
    up_p1 = sorted({m2h[g] for g in up_p1_m if g in m2h and isinstance(m2h[g], str)})
    dn_p1 = sorted({m2h[g] for g in dn_p1_m if g in m2h and isinstance(m2h[g], str)})

    up_p0_no_inflam = [g for g in up_p0 if g not in INFLAMMATORY_STRESS_GENES]
    dn_p0_no_inflam = [g for g in dn_p0 if g not in INFLAMMATORY_STRESS_GENES]

    return {
        "P0_Raw257": (up_p0, dn_p0),
        "P0_Minus_Inflammatory": (up_p0_no_inflam, dn_p0_no_inflam),
        "P1_Clean": (up_p1, dn_p1)
    }

def fit_unadj(y, keep, age):
    y_sub = y[keep]
    a = age[keep]
    X = np.column_stack([np.ones(len(y_sub)), a])
    beta, *_ = np.linalg.lstsq(X, y_sub, rcond=None)
    resid = y_sub - X @ beta
    dof = len(y_sub) - X.shape[1]
    s2 = (resid @ resid) / max(dof, 1)
    cov = s2 * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))[1]
    tstat = beta[1] / se
    p = 2 * stats.t.sf(abs(tstat), dof)
    return beta[1], se, p

def fit_covar_adj(y, keep, age, sex, brain_bank, pmi, tech_pc1, tech_pc2, tech_pc3):
    y_sub = y[keep]
    a = age[keep]
    s_male = (sex[keep] == "Male").astype(float)
    bb_hbcc = (brain_bank[keep] == "HBCC").astype(float)
    pmi_k = pmi[keep]
    pc1 = tech_pc1[keep]
    pc2 = tech_pc2[keep]
    pc3 = tech_pc3[keep]

    X = np.column_stack([np.ones(len(y_sub)), a, s_male, bb_hbcc, pmi_k, pc1, pc2, pc3])
    beta, *_ = np.linalg.lstsq(X, y_sub, rcond=None)
    resid = y_sub - X @ beta
    dof = len(y_sub) - X.shape[1]
    s2 = (resid @ resid) / max(dof, 1)
    cov = s2 * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))[1]
    tstat = beta[1] / se
    p = 2 * stats.t.sf(abs(tstat), dof)
    return beta[1], se, p

def main():
    print("[D1 Driver Decomposition] Loading lifespan_pseudobulk.rds ...")
    rds_path = REPO / "data" / "external" / "dlpfc-atlas" / "lifespan_pseudobulk.rds"
    obj = rdata.read_rds(str(rds_path))
    cd = obj.colData
    donor_ids = list(cd.rownames)
    age = np.asarray(cd.listData["Age"], float).ravel()
    sex = np.asarray(cd.listData["Sex"]).astype(str).ravel()
    brain_bank = np.asarray(cd.listData["Brain_bank"]).astype(str).ravel()
    pmi = np.asarray(cd.listData["PMI"], float).ravel()
    tech_pc1 = np.asarray(cd.listData["TechPC1"], float).ravel()
    tech_pc2 = np.asarray(cd.listData["TechPC2"], float).ravel()
    tech_pc3 = np.asarray(cd.listData["TechPC3"], float).ravel()

    programs = load_programs()
    up_p0, dn_p0 = programs["P0_Raw257"]

    # 1. Per-gene driver analysis for Microglia and OPC
    gene_records = []
    model_records = []

    for ct in ["Micro", "OPC"]:
        m = obj.assays.data.listData[ct]
        genes = [str(g) for g in m.Dimnames[0]]
        gidx = {g: i for i, g in enumerate(genes)}
        mat = sp.csc_matrix((np.asarray(m.x), np.asarray(m.i), np.asarray(m.p)), shape=tuple(m.Dim))
        keep = np.asarray(mat.sum(axis=0)).ravel() > 0

        # Normalization
        lib = np.asarray(mat.sum(axis=0)).ravel()
        lib = np.where(lib > 0, lib, 1.0)
        X = mat.astype(np.float64) @ sp.diags(1e6 / lib)
        Xl = X.copy()
        Xl.data = np.log1p(Xl.data)

        # Decompose P0 genes
        up_matched = [g for g in up_p0 if g in gidx]
        dn_matched = [g for g in dn_p0 if g in gidx]
        n_up = len(up_matched)
        n_dn = len(dn_matched)

        for g in up_matched:
            expr = np.asarray(Xl[gidx[g], :].toarray()).ravel()
            b, se, p = fit_covar_adj(expr, keep, age, sex, brain_bank, pmi, tech_pc1, tech_pc2, tech_pc3)
            contrib = b / n_up
            gene_records.append({
                "cell_type": "Microglia" if ct == "Micro" else "OPC",
                "gene": g,
                "direction": "UP",
                "is_inflammatory_stress": g in INFLAMMATORY_STRESS_GENES,
                "mean_cpm": np.expm1(np.mean(expr[keep])),
                "gene_slope": b,
                "gene_se": se,
                "gene_p": p,
                "score_slope_contribution": contrib
            })

        for g in dn_matched:
            expr = np.asarray(Xl[gidx[g], :].toarray()).ravel()
            b, se, p = fit_covar_adj(expr, keep, age, sex, brain_bank, pmi, tech_pc1, tech_pc2, tech_pc3)
            contrib = -b / n_dn
            gene_records.append({
                "cell_type": "Microglia" if ct == "Micro" else "OPC",
                "gene": g,
                "direction": "DOWN",
                "is_inflammatory_stress": g in INFLAMMATORY_STRESS_GENES,
                "mean_cpm": np.expm1(np.mean(expr[keep])),
                "gene_slope": b,
                "gene_se": se,
                "gene_p": p,
                "score_slope_contribution": contrib
            })

        # Evaluate Models
        for prog_name, (up_g, dn_g) in programs.items():
            up_i = [gidx[g] for g in up_g if g in gidx]
            dn_i = [gidx[g] for g in dn_g if g in gidx]
            up_m_score = np.asarray(Xl[up_i, :].mean(axis=0)).ravel() if up_i else np.zeros(Xl.shape[1])
            dn_m_score = np.asarray(Xl[dn_i, :].mean(axis=0)).ravel() if dn_i else np.zeros(Xl.shape[1])
            score = up_m_score - dn_m_score
            b, se, p = fit_covar_adj(score, keep, age, sex, brain_bank, pmi, tech_pc1, tech_pc2, tech_pc3)
            model_records.append({
                "compartment": "Microglia" if ct == "Micro" else "OPC",
                "program_tier": prog_name,
                "n_up_mapped": len(up_i),
                "n_down_mapped": len(dn_i),
                "slope": b,
                "se": se,
                "p": p
            })

    # Evaluate Whole Cortex across the 3 program tiers
    cts = sorted(obj.assays.data.listData.keys())
    for prog_name, (up_g, dn_g) in programs.items():
        donor_scores = {d: [] for d in donor_ids}
        for ct in cts:
            m = obj.assays.data.listData[ct]
            genes = [str(g) for g in m.Dimnames[0]]
            gidx = {g: i for i, g in enumerate(genes)}
            mat = sp.csc_matrix((np.asarray(m.x), np.asarray(m.i), np.asarray(m.p)), shape=tuple(m.Dim))
            keep = np.asarray(mat.sum(axis=0)).ravel() > 0
            lib = np.asarray(mat.sum(axis=0)).ravel()
            lib = np.where(lib > 0, lib, 1.0)
            X = mat.astype(np.float64) @ sp.diags(1e6 / lib)
            Xl = X.copy()
            Xl.data = np.log1p(Xl.data)

            up_i = [gidx[g] for g in up_g if g in gidx]
            dn_i = [gidx[g] for g in dn_g if g in gidx]
            up_m_score = np.asarray(Xl[up_i, :].mean(axis=0)).ravel() if up_i else np.zeros(Xl.shape[1])
            dn_m_score = np.asarray(Xl[dn_i, :].mean(axis=0)).ravel() if dn_i else np.zeros(Xl.shape[1])
            sc = up_m_score - dn_m_score

            for d, s, present in zip(donor_ids, sc, keep):
                if present:
                    donor_scores[d].append(s)

        whole_scores = np.array([np.mean(donor_scores[d]) for d in donor_ids])
        b_u, se_u, p_u = fit_unadj(whole_scores, np.ones(len(donor_ids), bool), age)
        b_c, se_c, p_c = fit_covar_adj(whole_scores, np.ones(len(donor_ids), bool), age, sex, brain_bank, pmi, tech_pc1, tech_pc2, tech_pc3)
        model_records.append({
            "compartment": "Whole Cortex (Unadjusted)",
            "program_tier": prog_name,
            "n_up_mapped": len(up_g),
            "n_down_mapped": len(dn_g),
            "slope": b_u,
            "se": se_u,
            "p": p_u
        })
        model_records.append({
            "compartment": "Whole Cortex (Adjusted)",
            "program_tier": prog_name,
            "n_up_mapped": len(up_g),
            "n_down_mapped": len(dn_g),
            "slope": b_c,
            "se": se_c,
            "p": p_c
        })

    # Save
    gene_df = pd.DataFrame(gene_records)
    gene_df.to_csv(REPO / "results" / "d1-glial-driver-decomposition.csv", index=False)

    model_df = pd.DataFrame(model_records)
    model_df.to_csv(REPO / "results" / "d1-glial-model-summary.csv", index=False)

    print(f"Saved D1 driver decomposition ({len(gene_df)} rows) to results/d1-glial-driver-decomposition.csv")
    print(f"Saved D1 model summary to results/d1-glial-model-summary.csv:")
    print(model_df.to_string())

if __name__ == "__main__":
    main()
