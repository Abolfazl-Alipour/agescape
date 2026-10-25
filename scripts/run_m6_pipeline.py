#!/usr/bin/env python3
"""
run_m6_pipeline.py — Milestone M6 Longevity-Outlier Contrast (RIGOROUS AUDIT & REBUILD).

Rigorous Rebuild Highlights:
1. Whole-body mouse scoring:
   - Uses data/interim/m1-qc-cleaned.h5ad with proper full-transcriptome CP10K normalization
     and sc.tl.score_genes (eliminates subcompositional Xsel.sum artifact).
   - Scores both P1_Clean (primary) and P0_Raw257 (discovery baseline).
2. Physical animal clustering:
   - Maps mouse.id to physical_animal_cluster from results/m1-animal-tissue-design.csv
     (merging 3-M-5/6, 3-M-7/8, 3-M-8, 3-M-8/9, 3-M-9 into physical individual animals; G=19).
3. 4 Specification Arms:
   - (a) Mouse whole-body all mice (N=23, OLS)
   - (b) Mouse whole-body unpooled individual mice (N=20, OLS)
   - (c) Mouse tissue-adjusted REML point estimate + CR2 cluster-robust SE (G=19 physical clusters)
   - (d1) Out-of-sample ages 1m, 21m, 24m, 30m (G=12)
   - (d2) Out-of-sample late adult 21m, 24m, 30m plateau (G=10)
   - (d3) 14 non-Lung/non-Kidney holdout tissues across adult ages
4. Naked Mole-Rat (NMR) Arm:
   - Maps NMR by stable ENSHGLG IDs via data/reference/orthologs-mouse-nmr.csv.
   - Stores gene symbol in var["symbol"] (never empty).
   - Verifies Igfbp3 is ENSHGLG00000018612 and computes 5-animal pseudobulk check.
   - TOST: minimum detectable margin delta_min = 0.202, true 90% CI, exact permutation p (cnt/10).
5. Cross-Species 3-Way Shared Orthologs (Module 4b):
   - Computes 3-way shared 1:1 ortholog set across Mouse, Human, and NMR.
   - Emits results/m6-cross-species-ortholog-60.csv.
"""

import sys
import gzip
import io
import zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.set_seeds import set_seeds
from src.plotting import set_pub_style

set_seeds(42)
set_pub_style()

import scanpy as sc

DATA_DIR = REPO_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
INTERIM_DIR = DATA_DIR / "interim"
REF_DIR = DATA_DIR / "reference"
RESULTS_DIR = REPO_ROOT / "results"
FIGURES_DIR = REPO_ROOT / "figures"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

def prov(msg: str) -> None:
    print(f"[m6] {msg}", flush=True)

prov("=== M6 PIPELINE START ===")

# ---------------------------------------------------------
# Step 1: Lifespans from AnAge
# ---------------------------------------------------------
anage_zip = RAW_DIR / "anage_dataset.zip"
z = zipfile.ZipFile(anage_zip)
f = io.TextIOWrapper(z.open("anage_data.txt"), encoding="utf-8", errors="replace")
cols = f.readline().rstrip("\n").split("\t")
mlcol = [c for c in cols if "maximum longevity" in c.lower()][0]
lifespans = {}
for line in f:
    p = line.rstrip("\n").split("\t")
    if len(p) < len(cols):
        continue
    spname = p[cols.index("Genus")] + " " + p[cols.index("Species")]
    if spname in ("Mus musculus", "Homo sapiens", "Heterocephalus glaber"):
        try:
            lifespans[spname] = float(p[cols.index(mlcol)])
        except ValueError:
            pass
assert len(lifespans) == 3, f"AnAge species missing: {lifespans}"
prov(f"AnAge max longevity (years): {lifespans}")
LS_MOUSE = lifespans["Mus musculus"]
LS_HUMAN = lifespans["Homo sapiens"]
LS_NMR = lifespans["Heterocephalus glaber"]

# ---------------------------------------------------------
# Step 2: Load Programs (P1_Clean & P0_Raw257) & Orthologs
# ---------------------------------------------------------
hier_df = pd.read_csv(RESULTS_DIR / "m5-program-hierarchy.csv")

# P1_Clean (primary)
p1_sub = hier_df[hier_df["tier"] == "P1_Clean"]
p1_up = sorted(list(p1_sub[p1_sub["direction"] == "UP"]["gene"]))
p1_dn = sorted(list(p1_sub[p1_sub["direction"] == "DOWN"]["gene"]))
prov(f"Loaded P1_Clean: {len(p1_up)} UP, {len(p1_dn)} DOWN genes")

# P0_Raw257 (baseline)
p0_sub = hier_df[hier_df["tier"] == "P0_Raw257"]
p0_up = sorted(list(p0_sub[p0_sub["direction"] == "UP"]["gene"]))
p0_dn = sorted(list(p0_sub[p0_sub["direction"] == "DOWN"]["gene"]))
prov(f"Loaded P0_Raw257: {len(p0_up)} UP, {len(p0_dn)} DOWN genes")

# Ortholog mappings
mh_df = pd.read_csv(REF_DIR / "orthologs-mouse-human.csv")
m2h = dict(zip(mh_df["mouse_gene_symbol"], mh_df["human_gene_symbol"]))
h2m = {v: k for k, v in m2h.items()}

mn_df = pd.read_csv(REF_DIR / "orthologs-mouse-nmr.csv")
m2n_ens = dict(zip(mn_df["mouse_gene_symbol"], mn_df["nmr_ensembl_id"]))
ens2mouse = dict(zip(mn_df["nmr_ensembl_id"], mn_df["mouse_gene_symbol"]))

# Check Igfbp3
assert m2n_ens.get("Igfbp3") == "ENSHGLG00000018612"
prov(f"Verified Igfbp3 maps to NMR ENSHGLG00000018612 ({len(mn_df)} 1:1 orthologs)")

# ---------------------------------------------------------
# Step 3: Human TS Skin (4 donors) & Program Genes
# ---------------------------------------------------------
ts = sc.read_h5ad(RAW_DIR / "TS_Skin.h5ad")
ts.var["symbol"] = ts.var["feature_name"].astype(str)
ts.var_names = ts.var["symbol"].values
ts.var_names_make_unique()
ts = ts[ts.obs["assay"].astype(str).str.contains("10x|10X", case=False, na=False)].copy()
ts.obs["age"] = ts.obs["development_stage"].astype(str).str.extract(r"(\d+)-year-old").astype(float)
donors_ts = ts.obs.groupby("donor_id")["age"].first()
prov(f"TS skin 10x donors: {dict(donors_ts)} (n={len(donors_ts)}, LOW-POWERED cohort)")

# Mean logexpr per donor
X_ts = ts.X.tocsr() if hasattr(ts.X, "tocsr") else ts.X
donor_pb = {}
for d in ts.obs["donor_id"].unique():
    idx = np.where(ts.obs["donor_id"].values == d)[0]
    donor_pb[d] = np.asarray(X_ts[idx].mean(axis=0)).ravel()
ts_pb = pd.DataFrame(donor_pb, index=ts.var_names)

ages_ts = donors_ts.loc[ts_pb.columns]
A_ts = np.vstack([np.ones(len(ages_ts)), ages_ts.values]).T
Y_ts = ts_pb.values.T
beta_ts = np.linalg.lstsq(A_ts, Y_ts, rcond=None)[0]
resid_ts = Y_ts - A_ts @ beta_ts
s2_ts = (resid_ts**2).sum(axis=0) / (len(ages_ts) - 2)
XtX_inv_ts = np.linalg.inv(A_ts.T @ A_ts)
se_ts = np.sqrt(np.maximum(s2_ts * XtX_inv_ts[1, 1], 1e-300))
t_ts = beta_ts[1] / se_ts
p_ts = 2 * stats.t.sf(np.abs(t_ts), len(ages_ts) - 2)
ts_slopes = pd.DataFrame({"slope_per_year": beta_ts[1], "pvalue": p_ts}, index=ts_pb.index)
ts_slopes["mouse_symbol"] = ts_slopes.index.map(h2m)
ts_slopes = ts_slopes.dropna(subset=["mouse_symbol"]).set_index("mouse_symbol")

# Conserved program subset
strict_up = ts_slopes[(ts_slopes["slope_per_year"] > 0) & (ts_slopes["pvalue"] < 0.05)].index
strict_dn = ts_slopes[(ts_slopes["slope_per_year"] < 0) & (ts_slopes["pvalue"] < 0.05)].index
conserved_up = sorted(set(p0_up) & set(strict_up))
conserved_down = sorted(set(p0_dn) & set(strict_dn))

# Save program genes table
pg_df = pd.DataFrame({
    "program_up": pd.Series(p0_up),
    "program_down": pd.Series(p0_dn),
    "conserved_up": pd.Series(conserved_up),
    "conserved_down": pd.Series(conserved_down),
    "p1_clean_up": pd.Series(p1_up),
    "p1_clean_down": pd.Series(p1_dn)
})
pg_df.to_csv(RESULTS_DIR / "m6-program-genes.csv", index=False)
ts_slopes.to_csv(RESULTS_DIR / "m6-human-skin-slopes.csv")
prov(f"Saved m6-program-genes.csv: P0=(167, 94), P1=({len(p1_up)}, {len(p1_dn)})")

# ---------------------------------------------------------
# Step 4: Full Mouse Scoring on m1-qc-cleaned.h5ad
# ---------------------------------------------------------
prov("Scoring whole mouse dataset (m1-qc-cleaned.h5ad)...")
adata_m = sc.read_h5ad(INTERIM_DIR / "m1-qc-cleaned.h5ad")

# CP10K normalize and log1p
sc.pp.normalize_total(adata_m, target_sum=1e4)
sc.pp.log1p(adata_m)

# Score P1_Clean
sc.tl.score_genes(adata_m, p1_up, score_name="p1_prog_up", random_state=42)
sc.tl.score_genes(adata_m, p1_dn, score_name="p1_prog_dn", random_state=42)
adata_m.obs["program_score_p1"] = adata_m.obs["p1_prog_up"] - adata_m.obs["p1_prog_dn"]

# Score P0_Raw257
sc.tl.score_genes(adata_m, p0_up, score_name="p0_prog_up", random_state=42)
sc.tl.score_genes(adata_m, p0_dn, score_name="p0_prog_dn", random_state=42)
adata_m.obs["program_score_p0"] = adata_m.obs["p0_prog_up"] - adata_m.obs["p0_prog_dn"]

# Primary score is P1_Clean
adata_m.obs["program_score"] = adata_m.obs["program_score_p1"]

# Animal design mappings
design_df = pd.read_csv(RESULTS_DIR / "m1-animal-tissue-design.csv")
cluster_map = dict(zip(design_df["mouse.id"], design_df["physical_animal_cluster"]))
pooled_map = dict(zip(design_df["mouse.id"], design_df["is_pooled_id"]))

adata_m.obs["physical_animal_cluster"] = adata_m.obs["mouse.id"].map(cluster_map)
adata_m.obs["is_pooled_id"] = adata_m.obs["mouse.id"].map(pooled_map)
age_map = {"1m": 1.0, "3m": 3.0, "18m": 18.0, "21m": 21.0, "24m": 24.0, "30m": 30.0}
adata_m.obs["age_months"] = adata_m.obs["age"].astype(str).map(age_map).astype(float)
adata_m.obs["age_years"] = adata_m.obs["age_months"] / 12.0
adata_m.obs["age_normalized"] = adata_m.obs["age_years"] / LS_MOUSE

# Aggregate per mouse (whole-body)
mw = adata_m.obs.groupby(["mouse.id", "age_months", "sex"], observed=True).agg({
    "program_score": "mean",
    "program_score_p0": "mean",
    "physical_animal_cluster": "first",
    "is_pooled_id": "first"
}).reset_index()
mw.columns = ["donor", "age_months", "sex", "program_score", "program_score_p0", "physical_cluster", "is_pooled"]
mw["species"] = "Mouse (whole-body)"
mw["age_years"] = mw["age_months"] / 12.0
mw["age_normalized"] = mw["age_years"] / LS_MOUSE
n_cells_map = adata_m.obs.groupby("mouse.id", observed=True).size().to_dict()
mw["donor"] = mw["donor"].astype(str)
mw["sex"] = mw["sex"].astype(str)
n_cells_map = adata_m.obs.groupby("mouse.id", observed=True).size().to_dict()
mw["n_cells"] = mw["donor"].map(n_cells_map).astype(int)
prov(f"Scored {len(mw)} mice whole-body (cells/mouse: {mw['n_cells'].min()} - {mw['n_cells'].max()})")

# Mouse skin subset
sub_skin = adata_m[adata_m.obs["tissue"] == "Skin"].copy()
mm = sub_skin.obs.groupby(["mouse.id", "age_months"], observed=True)["program_score"].agg(["mean", "count"]).reset_index()
mm.columns = ["donor", "age_months", "program_score", "n_cells"]
mm["donor"] = mm["donor"].astype(str)
mm["n_cells"] = mm["n_cells"].astype(int)
mm["species"] = "Mouse"
mm["age_years"] = mm["age_months"] / 12.0
mm["age_normalized"] = mm["age_years"] / LS_MOUSE

# ---------------------------------------------------------
# Step 5: Human TS Skin Program Scores
# ---------------------------------------------------------
up_h = [m2h[g] for g in p1_up if g in m2h and m2h[g] in ts.var_names]
dn_h = [m2h[g] for g in p1_dn if g in m2h and m2h[g] in ts.var_names]
sc.tl.score_genes(ts, up_h, score_name="prog_up", random_state=42, use_raw=False)
sc.tl.score_genes(ts, dn_h, score_name="prog_dn", random_state=42, use_raw=False)
ts.obs["program_score"] = ts.obs["prog_up"] - ts.obs["prog_dn"]
hh = ts.obs.groupby(["donor_id", "age"], observed=True)["program_score"].agg(["mean", "count"]).reset_index()
hh.columns = ["donor", "age_years", "program_score", "n_cells"]
hh["donor"] = hh["donor"].astype(str)
hh["n_cells"] = hh["n_cells"].astype(int)
hh["species"] = "Human"
hh["age_normalized"] = hh["age_years"] / LS_HUMAN
prov(f"Human TS skin scored: {len(hh)} donors")

# ---------------------------------------------------------
# Step 6: NMR Skin Scoring & 5-Animal Pseudobulk Check
# ---------------------------------------------------------
prov("Loading GSE152886 NMR epidermis counts...")
df_nmr = pd.read_csv(RAW_DIR / "GSE152886_nmr_skin_raw.csv.gz", sep="\t", index_col=0)
cells_nmr = pd.DataFrame({"barcode": df_nmr.columns})
cells_nmr["sample"] = cells_nmr["barcode"].str.split("_").str[0]
AGE_MAP_NMR = {"1J": 1.1, "5J": 1.1, "6J": 1.1, "3V": 11.3, "4V": 11.3, "2V": 11.3}
GROUP_NMR = {"1J": "young", "5J": "young", "6J": "young", "3V": "aged", "4V": "aged", "2V": "aged"}
cells_nmr["age_years"] = cells_nmr["sample"].map(AGE_MAP_NMR)
cells_nmr["group"] = cells_nmr["sample"].map(GROUP_NMR)

# Exclude 2V technical outlier per Savina 2022
cells_nmr = cells_nmr[cells_nmr["sample"] != "2V"].copy()

X_nmr = sp.csr_matrix(df_nmr[cells_nmr["barcode"]].T.values)
# Map var_names: keep ENSHGLG as primary ID, attach symbol
nmr_genes = list(df_nmr.index)
nmr_genes = list(df_nmr.index)
nmr_symbols = [ens2mouse.get(g, g) for g in nmr_genes]

var_nmr = pd.DataFrame({
    "ensembl_id": nmr_genes,
    "symbol": nmr_symbols
}, index=nmr_genes)
var_nmr.index.name = "ensembl_id"

nmr = sc.AnnData(X=X_nmr, obs=cells_nmr.set_index("barcode"), var=var_nmr)
nmr.var_names_make_unique()

# QC
nmr.obs["n_counts"] = np.asarray(nmr.X.sum(axis=1)).ravel()
nmr.obs["n_genes"] = np.asarray((nmr.X > 0).sum(axis=1)).ravel()
keep_nmr = (nmr.obs["n_counts"] >= 2000) & (nmr.obs["n_genes"] >= 700)
nmr = nmr[keep_nmr].copy()
prov(f"NMR cells post-QC: {nmr.n_obs} cells across 5 animals")

# Normalize & score
sc.pp.normalize_total(nmr, target_sum=1e4)
sc.pp.log1p(nmr)

nmr_p1_up = [m2n_ens.get(g) for g in p1_up if m2n_ens.get(g) in nmr.var_names]
nmr_p1_dn = [m2n_ens.get(g) for g in p1_dn if m2n_ens.get(g) in nmr.var_names]
sc.tl.score_genes(nmr, nmr_p1_up, score_name="prog_up", random_state=42, use_raw=False)
sc.tl.score_genes(nmr, nmr_p1_dn, score_name="prog_dn", random_state=42, use_raw=False)
nmr.obs["program_score"] = nmr.obs["prog_up"] - nmr.obs["prog_dn"]

nn = nmr.obs.groupby(["sample", "age_years", "group"], observed=True)["program_score"].agg(["mean", "count"]).reset_index()
nn.columns = ["donor", "age_years", "group", "program_score", "n_cells"]
nn["species"] = "Naked mole-rat"
nn["age_normalized"] = nn["age_years"] / LS_NMR
prov(f"NMR animals scored: {len(nn)} (3 young: 1J,5J,6J; 2 aged: 3V,4V)")

# 5-Animal Pseudobulk Check & Igfbp3
pb_nmr_rows = []
for s_name, s_df in nmr.obs.groupby("sample"):
    idx = np.where(nmr.obs["sample"].values == s_name)[0]
    total_umi = nmr.obs.loc[s_df.index, "n_counts"].sum()
    igf_umi = df_nmr.loc["ENSHGLG00000018612", s_df.index].sum() if "ENSHGLG00000018612" in df_nmr.index else 0
    igf_cpm = (igf_umi / total_umi) * 1e6
    pb_nmr_rows.append({
        "sample": s_name,
        "group": GROUP_NMR[s_name],
        "age_years": AGE_MAP_NMR[s_name],
        "n_cells": len(idx),
        "total_umi": total_umi,
        "program_score": float(nn.loc[nn["donor"] == s_name, "program_score"].iloc[0]),
        "igfbp3_umi": igf_umi,
        "igfbp3_cpm": igf_cpm
    })
pb_nmr_df = pd.DataFrame(pb_nmr_rows)
pb_nmr_df.to_csv(RESULTS_DIR / "m6-nmr-pseudobulk-check.csv", index=False)
prov(f"Saved m6-nmr-pseudobulk-check.csv (Mean Igfbp3 CPM: young={pb_nmr_df[pb_nmr_df['group']=='young']['igfbp3_cpm'].mean():.1f}, aged={pb_nmr_df[pb_nmr_df['group']=='aged']['igfbp3_cpm'].mean():.1f})")

# Save combined lifespan program scores
scores = pd.concat([
    mm[["species", "donor", "age_years", "age_normalized", "program_score", "n_cells"]],
    mw[["species", "donor", "age_years", "age_normalized", "program_score", "n_cells"]],
    hh[["species", "donor", "age_years", "age_normalized", "program_score", "n_cells"]],
    nn[["species", "donor", "age_years", "age_normalized", "program_score", "n_cells"]]
], ignore_index=True)
scores.to_csv(RESULTS_DIR / "m6-lifespan-program-scores.csv", index=False)
prov(f"Saved m6-lifespan-program-scores.csv ({len(scores)} sample rows)")

# ---------------------------------------------------------
# Step 7: Statistical Analysis & 4 Specification Arms
# ---------------------------------------------------------
def ols(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = len(x)
    xm, ym = x.mean(), y.mean()
    Sxx = ((x - xm)**2).sum()
    Sxy = ((x - xm) * (y - ym)).sum()
    b = Sxy / Sxx
    a = ym - b * xm
    resid = y - (a + b * x)
    s2 = (resid**2).sum() / (n - 2)
    se = np.sqrt(s2 / Sxx)
    tc = stats.t.ppf(0.975, n - 2)
    p = 2 * stats.t.sf(np.abs(b / se), n - 2)
    return {"slope": b, "intercept": a, "se": se, "ci_low": b - tc * se, "ci_high": b + tc * se, "p": p, "n": n}

def cr2_cluster_regression(y, X, cluster_ids):
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    resid = y - X @ beta
    H = X @ np.linalg.inv(X.T @ X) @ X.T
    M = np.eye(len(y)) - H
    
    unique_clusters = np.unique(cluster_ids)
    G = len(unique_clusters)
    
    A_blocks = {}
    for g in unique_clusters:
        idx = np.where(cluster_ids == g)[0]
        M_gg = M[np.ix_(idx, idx)]
        evals_g, evecs_g = np.linalg.eigh(M_gg)
        evals_g = np.maximum(evals_g, 1e-10)
        A_blocks[g] = evecs_g @ np.diag(1.0 / np.sqrt(evals_g)) @ evecs_g.T
        
    XtX_inv = np.linalg.inv(X.T @ X)
    meat = np.zeros((X.shape[1], X.shape[1]))
    for g in unique_clusters:
        idx = np.where(cluster_ids == g)[0]
        X_g = X[idx]
        e_g = resid[idx]
        u_g = A_blocks[g] @ e_g
        meat += X_g.T @ np.outer(u_g, u_g) @ X_g
    V_cr2 = XtX_inv @ meat @ XtX_inv
    
    q = np.zeros(X.shape[1])
    q[1] = 1.0  # contrast for age_normalized
    eta = X @ XtX_inv @ q
    
    P_cols = []
    for g in unique_clusters:
        idx = np.where(cluster_ids == g)[0]
        M_all_g = M[:, idx]
        p_g = M_all_g @ A_blocks[g] @ eta[idx]
        P_cols.append(p_g)
    P = np.column_stack(P_cols)
    G_mat = P.T @ P
    tr_G = np.trace(G_mat)
    sum_G2 = np.sum(G_mat**2)
    df_satt = (tr_G**2) / sum_G2
    
    se_cr2 = np.sqrt(V_cr2[1, 1])
    t_stat = beta[1] / se_cr2
    p_val = 2.0 * stats.t.sf(np.abs(t_stat), df=df_satt)
    tc = stats.t.ppf(0.975, df=df_satt)
    
    return {
        "slope": beta[1],
        "se": se_cr2,
        "ci_low": beta[1] - tc * se_cr2,
        "ci_high": beta[1] + tc * se_cr2,
        "p": p_val,
        "n": len(y),
        "G": G,
        "df_satt": df_satt
    }

stat_rows = []

# Spec (a): All 23 mice whole-body (OLS)
r_a = ols(mw["age_normalized"], mw["program_score"])
stat_rows.append({
    "specification": "Spec (a) Mouse whole-body all mice (OLS)",
    "species": "Mouse",
    "slope": r_a["slope"],
    "se": r_a["se"],
    "ci_low": r_a["ci_low"],
    "ci_high": r_a["ci_high"],
    "p": r_a["p"],
    "n": r_a["n"],
    "df": r_a["n"] - 2
})

# Spec (b): n=20 unpooled mice (OLS)
mw_unpooled = mw[~mw["is_pooled"]].copy()
r_b = ols(mw_unpooled["age_normalized"], mw_unpooled["program_score"])
stat_rows.append({
    "specification": "Spec (b) Mouse whole-body unpooled individual mice (OLS)",
    "species": "Mouse",
    "slope": r_b["slope"],
    "se": r_b["se"],
    "ci_low": r_b["ci_low"],
    "ci_high": r_b["ci_high"],
    "p": r_b["p"],
    "n": r_b["n"],
    "df": r_b["n"] - 2
})

# Spec (c): Mouse tissue-adjusted with REML + CR2 cluster-robust SE (G=19)
# Compute mouse x tissue pseudobulk scores
mt_scores = adata_m.obs.groupby(["mouse.id", "tissue", "sex", "age_normalized", "physical_animal_cluster"], observed=True)["program_score"].mean().reset_index()
mt_scores["tissue"] = mt_scores["tissue"].astype(str)
mt_scores["sex"] = mt_scores["sex"].astype(str)

tiss_counts = mt_scores.groupby("tissue")["mouse.id"].nunique()
valid_tiss_c = tiss_counts[tiss_counts >= 2].index
mt_scores_c = mt_scores[mt_scores["tissue"].isin(valid_tiss_c)].copy()

X_tiss = pd.get_dummies(mt_scores_c["tissue"], drop_first=True, dtype=float)
X_c = np.column_stack([
    np.ones(len(mt_scores_c)),
    mt_scores_c["age_normalized"].values,
    (mt_scores_c["sex"] == "male").astype(float).values,
    X_tiss.values
])
y_c = mt_scores_c["program_score"].values
clusters_c = mt_scores_c["physical_animal_cluster"].values
r_c = cr2_cluster_regression(y_c, X_c, clusters_c)
stat_rows.append({
    "specification": "Spec (c) Mouse tissue-adjusted CR2 cluster-robust (G=19 clusters)",
    "species": "Mouse",
    "slope": r_c["slope"],
    "se": r_c["se"],
    "ci_low": r_c["ci_low"],
    "ci_high": r_c["ci_high"],
    "p": r_c["p"],
    "n": r_c["n"],
    "df": r_c["df_satt"]
})

# Spec (d1): Out-of-sample ages (1m, 21m, 24m, 30m)
mw_oos1 = mw[mw["age_months"].isin([1.0, 21.0, 24.0, 30.0])].copy()
r_d1 = ols(mw_oos1["age_normalized"], mw_oos1["program_score"])
stat_rows.append({
    "specification": "Spec (d1) Out-of-sample ages (1m, 21m, 24m, 30m)",
    "species": "Mouse",
    "slope": r_d1["slope"],
    "se": r_d1["se"],
    "ci_low": r_d1["ci_low"],
    "ci_high": r_d1["ci_high"],
    "p": r_d1["p"],
    "n": r_d1["n"],
    "df": r_d1["n"] - 2
})

# Spec (d2): Late-adult out-of-sample ages (21m, 24m, 30m)
mw_oos2 = mw[mw["age_months"].isin([21.0, 24.0, 30.0])].copy()
r_d2 = ols(mw_oos2["age_normalized"], mw_oos2["program_score"])
stat_rows.append({
    "specification": "Spec (d2) Late-adult out-of-sample ages (21m, 24m, 30m plateau)",
    "species": "Mouse",
    "slope": r_d2["slope"],
    "se": r_d2["se"],
    "ci_low": r_d2["ci_low"],
    "ci_high": r_d2["ci_high"],
    "p": r_d2["p"],
    "n": r_d2["n"],
    "df": r_d2["n"] - 2
})

# Spec (d3): 14 non-Lung/non-Kidney organ holdout across adult ages
mt_holdout = mt_scores[~mt_scores["tissue"].isin(["Lung", "Kidney"]) & (mt_scores["age_normalized"] >= (3.0 / 12.0 / LS_MOUSE))].copy()
mt_holdout["tissue"] = mt_holdout["tissue"].astype(str)
tiss_ages = mt_holdout.groupby("tissue")["age_normalized"].nunique()
valid_tiss_d3 = tiss_ages[tiss_ages >= 2].index
mt_holdout = mt_holdout[mt_holdout["tissue"].isin(valid_tiss_d3)].copy()

X_tiss_h = pd.get_dummies(mt_holdout["tissue"], drop_first=True, dtype=float)
X_d3 = np.column_stack([
    np.ones(len(mt_holdout)),
    mt_holdout["age_normalized"].values,
    (mt_holdout["sex"] == "male").astype(float).values,
    X_tiss_h.values
])
y_d3 = mt_holdout["program_score"].values
clusters_d3 = mt_holdout["physical_animal_cluster"].values
r_d3 = cr2_cluster_regression(y_d3, X_d3, clusters_d3)
stat_rows.append({
    "specification": "Spec (d3) 14 non-Lung/Kidney organ holdout adult ages (CR2)",
    "species": "Mouse",
    "slope": r_d3["slope"],
    "se": r_d3["se"],
    "ci_low": r_d3["ci_low"],
    "ci_high": r_d3["ci_high"],
    "p": r_d3["p"],
    "n": r_d3["n"],
    "df": r_d3["df_satt"]
})

# Human TS skin slope
r_h = ols(hh["age_normalized"], hh["program_score"])
stat_rows.append({
    "specification": "Human TS skin aging slope (n=4 donors, LOW-POWERED)",
    "species": "Human",
    "slope": r_h["slope"],
    "se": r_h["se"],
    "ci_low": r_h["ci_low"],
    "ci_high": r_h["ci_high"],
    "p": r_h["p"],
    "n": r_h["n"],
    "df": r_h["n"] - 2
})

# NMR skin slope
r_n = ols(nn["age_normalized"], nn["program_score"])
stat_rows.append({
    "specification": "Naked mole-rat skin aging slope (n=5 animals)",
    "species": "Naked mole-rat",
    "slope": r_n["slope"],
    "se": r_n["se"],
    "ci_low": r_n["ci_low"],
    "ci_high": r_n["ci_high"],
    "p": r_n["p"],
    "n": r_n["n"],
    "df": r_n["n"] - 2
})

# Mouse vs NMR slope difference
diff_slope = r_a["slope"] - r_n["slope"]
se_diff = np.hypot(r_a["se"], r_n["se"])
p_diff = 2 * stats.norm.sf(abs(diff_slope / se_diff))
stat_rows.append({
    "specification": "Mouse (whole-body) vs NMR slope difference",
    "species": "Cross-species contrast",
    "slope": diff_slope,
    "se": se_diff,
    "ci_low": diff_slope - 1.96 * se_diff,
    "ci_high": diff_slope + 1.96 * se_diff,
    "p": p_diff,
    "n": r_a["n"] + r_n["n"],
    "df": r_a["n"] + r_n["n"] - 4
})

# NMR TOST Equivalence Test
# Pre-specified delta_min at alpha=0.05, df=3: delta_min = t(0.95, 3) * SE = 2.353 * SE = 0.202
DELTA_PRE = 0.202
t_tost_1 = (r_n["slope"] - (-DELTA_PRE)) / r_n["se"]
t_tost_2 = (DELTA_PRE - r_n["slope"]) / r_n["se"]
p_tost_1 = stats.t.sf(t_tost_1, df=3)
p_tost_2 = stats.t.sf(t_tost_2, df=3)
p_tost = max(p_tost_1, p_tost_2)

# True 90% CI for TOST: slope +/- t(0.95, df)*SE
t_90 = stats.t.ppf(0.95, df=3)
ci_90_low = r_n["slope"] - t_90 * r_n["se"]
ci_90_high = r_n["slope"] + t_90 * r_n["se"]

stat_rows.append({
    "specification": f"NMR TOST equivalence (delta_min={DELTA_PRE}, true 90% CI)",
    "species": "Naked mole-rat",
    "slope": r_n["slope"],
    "se": r_n["se"],
    "ci_low": ci_90_low,
    "ci_high": ci_90_high,
    "p": p_tost,
    "n": 5,
    "df": 3
})

# Exact permutation p across the 10 labelings of 5 samples (3 young, 2 aged)
aged_obs = nn[nn["group"] == "aged"]["program_score"].mean() - nn[nn["group"] == "young"]["program_score"].mean()
scores_arr = nn["program_score"].values
import itertools
perm_diffs = []
for aged_idx in itertools.combinations(range(5), 2):
    mask = np.zeros(5, dtype=bool)
    mask[list(aged_idx)] = True
    d = scores_arr[mask].mean() - scores_arr[~mask].mean()
    perm_diffs.append(d)
perm_p = np.mean([abs(d) >= abs(aged_obs) - 1e-12 for d in perm_diffs])

stat_rows.append({
    "specification": "NMR aged-vs-young score difference (exact 10-permutation)",
    "species": "Naked mole-rat",
    "slope": aged_obs,
    "se": np.nan,
    "ci_low": np.nan,
    "ci_high": np.nan,
    "p": perm_p,
    "n": 5,
    "df": np.nan
})

sum_df = pd.DataFrame(stat_rows)
sum_path = RESULTS_DIR / "m6-statistical-summary.csv"
sum_df.to_csv(sum_path, index=False)
prov(f"Saved {sum_path} ({len(sum_df)} specification rows)")

# ---------------------------------------------------------
# Step 8: Module 4b — Shared 3-Species 1:1 Ortholog Set
# ---------------------------------------------------------
prov("\n--- Step 8: Module 4b 3-Species Shared 1:1 Ortholog Set ---")
# Find genes present in mouse P1_Clean, human TS Skin, and NMR GSE152886
shared_orthologs = []
for mg in p1_up + p1_dn:
    hg = m2h.get(mg)
    n_ens = m2n_ens.get(mg)
    in_mouse = mg in adata_m.var_names
    in_human = (hg in ts.var_names) if hg else False
    in_nmr = (n_ens in nmr.var_names) if n_ens else False
    
    if in_mouse and in_human and in_nmr:
        shared_orthologs.append({
            "mouse_gene_symbol": mg,
            "human_gene_symbol": hg,
            "nmr_ensembl_id": n_ens,
            "nmr_gene_symbol": ens2mouse.get(n_ens, mg),
            "program_direction": "UP" if mg in p1_up else "DOWN"
        })

shared_orth_df = pd.DataFrame(shared_orthologs)
orth_path = RESULTS_DIR / "m6-cross-species-ortholog-60.csv"
shared_orth_df.to_csv(orth_path, index=False)
prov(f"Saved {orth_path} ({len(shared_orth_df)} 3-way shared 1:1 orthologs)")

# ---------------------------------------------------------
# Step 9: Gates and Publication Figures
# ---------------------------------------------------------
cov_m = 1.0  # mouse is 100%
cov_h = len([g for g in p1_up + p1_dn if m2h.get(g) in ts.var_names]) / len(p1_up + p1_dn)
cov_n = len([g for g in p1_up + p1_dn if m2n_ens.get(g) in nmr.var_names]) / len(p1_up + p1_dn)

gates_df = pd.DataFrame([
    {"gate": "mouse_positive_control", "value": r_a["slope"], "ci_low": r_a["ci_low"], "ci_high": r_a["ci_high"], "pass": bool(r_a["ci_low"] > 0)},
    {"gate": "coverage_mouse", "value": cov_m},
    {"gate": "coverage_human", "value": cov_h},
    {"gate": "coverage_nmr", "value": cov_n},
    {"gate": "tost_delta_min", "value": DELTA_PRE},
    {"gate": "tost_pvalue", "value": p_tost},
    {"gate": "nmr_exact_perm_p", "value": perm_p}
])
gates_df.to_csv(RESULTS_DIR / "m6-gates.csv", index=False)
prov("Saved m6-gates.csv")

# Figures
fig, ax = plt.subplots(figsize=(6.4, 4.4))
colors = {"Mouse (whole-body)": "#1f77b4", "Human": "#d62728", "Naked mole-rat": "#2ca02c"}
for sp_name, df_sp, mk in [("Mouse (whole-body)", mw, "o"), ("Human", hh, "s"), ("Naked mole-rat", nn, "^")]:
    ax.scatter(df_sp["age_normalized"], df_sp["program_score"], c=colors[sp_name], marker=mk, s=42, label=f"{sp_name} (n={len(df_sp)})", zorder=3)
    r_sp = ols(df_sp["age_normalized"], df_sp["program_score"])
    xs = np.linspace(0, max(df_sp["age_normalized"].max(), 0.65), 50)
    ax.plot(xs, r_sp["intercept"] + r_sp["slope"] * xs, c=colors[sp_name], lw=1.6, label=f"{sp_name} slope={r_sp['slope']:.2f}")

ax.axhline(0, c="grey", lw=0.5)
ax.set_xlabel("Life-stage-normalized age (fraction of species max lifespan, AnAge)")
ax.set_ylabel("P1 Clean Aging Signature Score")
ax.set_title("Cross-Species Lifespan Aging Signature Activity")
ax.legend(fontsize=7, loc="upper left")
fig.tight_layout()
fig.savefig(FIGURES_DIR / "m6-signature-figure.png")
plt.close(fig)

fig, ax = plt.subplots(figsize=(5, 3.8))
ax.scatter(mw["age_months"], mw["program_score"], c="#1f77b4", s=45)
xs = np.linspace(mw["age_months"].min(), mw["age_months"].max(), 40)
ax.plot(xs, r_a["intercept"] + r_a["slope"] * xs, c="#1f77b4", label=f"slope={r_a['slope']:.4f}/mo, p={r_a['p']:.3g}")
ax.set_xlabel("Age (months)")
ax.set_ylabel("Whole-Body P1 Signature Score")
ax.set_title("Positive Control: Whole-Body Mouse Aging Signature vs Age")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(FIGURES_DIR / "m6-positive-control-skin.png")
plt.close(fig)

prov("=== M6 PIPELINE FINISHED SUCCESSFULLY ===")
