#!/usr/bin/env python3
"""
scripts/run_m5_heterogeneity.py — Milestone M5 Program Hierarchy & Heterogeneity:
1. Program Tiers:
   - P0_Raw257: Discovery contrast DE genes (FDR < 0.05, |log2FC| >= 0.25).
   - P1_Clean: Uniform BH, floor >= 10 cells/mouse (12 surviving units), filtering Gm, ps, Rik,
     non-protein-coding, ribo-mito/Rn45s, and unit-level soup (Mgp/Col3a1 in lung immune; Ins1/2 in kidney).
   - P2_Clean_SexTrimmed: P1_Clean with adult sex-dimorphic genes excluded.
   - P3_Clean_SexIEGTrimmed: P2_Clean_SexTrimmed with Immediate Early Genes (IEGs) excluded.
   - Emits results/m5-program-hierarchy.csv.
2. Covariance-Calibrated Cochran-Fisher Q* Heterogeneity Decomposition:
   - Untruncated rank correlation C_hat = 2*sin(pi*rho_S/6), eigenvalue floor 1e-4.
   - Conservative Q*_corr and Empirical-Null Q*_empnull.
   - Exact nested GLS decomposition: Q_between = Q_total - Q_within_joint >= 0.
   - Emits results/m5-heterogeneity-qstar.csv.
3. External Aging Program Overlaps:
   - Zhang 2021 eLife 330 FACS GAGs, 59 droplet GAGs, Saul 2022 SenMayo mouse 118 genes.
   - Hypergeometric test across program tiers.
   - Emits results/m5-external-program-overlap.csv.
4. Generates figures/m5-heterogeneity-decomposition.png.
"""

import sys
import re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.set_seeds import set_seeds
from src.plotting import set_pub_style

set_seeds(42)
set_pub_style()

DATA_DIR = REPO_ROOT / "data"
RESULTS_DIR = REPO_ROOT / "results"
FIGURES_DIR = REPO_ROOT / "figures"
REF_DIR = DATA_DIR / "reference"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

print("=== M5 HETEROGENEITY & PROGRAM HIERARCHY START ===")

# ---------------------------------------------------------
# Step 1: Program Hierarchy (P0, P1, P2, P3)
# ---------------------------------------------------------
print("\n--- Step 1: Building Program Hierarchy Tiers ---")

# Chromosomes & Biotypes
chrom_df = pd.read_csv(REF_DIR / "mouse-gene-chromosomes.csv").set_index("mouse_gene_symbol")
biotypes = chrom_df["gene_biotype"].to_dict()
chromosomes = chrom_df["chromosome_name"].to_dict()

# Sex dimorphism audit
sex_audit = pd.read_csv(RESULTS_DIR / "m2-sex-dimorphism-audit.csv").set_index("mouse_gene_symbol")
sex_exclude = set(sex_audit[sex_audit["exclude_in_sex_trimmed"]].index)

# Surviving units at Floor10
spec_audit = pd.read_csv(RESULTS_DIR / "m2-de-specification-audit.csv")
surv_units = set(spec_audit[spec_audit["survives_floor10"]]["unit"])
print(f"Surviving units at Floor10 (>=10 cells/mouse): {len(surv_units)}")

# P0: Raw genes from m5-gene-level-programs.csv
p0_raw_df = pd.read_csv(RESULTS_DIR / "m5-gene-level-programs.csv")
p0_genes = sorted(list(set(p0_raw_df["gene"])))
print(f"P0_Raw unique genes: {len(p0_genes)}")

# P1: Clean from Uniform-BH genes
s1_genes = pd.read_csv(RESULTS_DIR / "m2-s1-uniform-bh-genes.csv")
p1_candidates = s1_genes[s1_genes["unit"].isin(surv_units)].copy()

def is_p1_clean(gene, unit):
    g = str(gene)
    if re.match(r"^Gm\d+$", g):
        return False
    if re.search(r"-ps\d*$", g):
        return False
    if re.match(r"^[0-9A-Z]+Rik$", g):
        return False
    if g.startswith(("Rpl", "Rps", "mt-", "Mrp", "MT-")) or g == "Rn45s":
        return False
    b = biotypes.get(g, "unknown")
    if b != "protein_coding":
        return False
    # Soup filters
    if "Lung" in unit and any(imm in unit for imm in ["B cell", "CD4", "NK", "monocyte", "macrophage"]) and g in ["Mgp", "Col3a1"]:
        return False
    if "Kidney" in unit and g in ["Ins1", "Ins2"]:
        return False
    return True

p1_mask = [is_p1_clean(r["mouse_gene_symbol"], r["unit"]) for _, r in p1_candidates.iterrows()]
p1_df = p1_candidates[p1_mask].copy()
p1_genes = sorted(list(set(p1_df["mouse_gene_symbol"])))
print(f"P1_Clean unique genes: {len(p1_genes)}")

# P2: Clean + SexTrimmed
p2_df = p1_df[~p1_df["mouse_gene_symbol"].isin(sex_exclude)].copy()
p2_genes = sorted(list(set(p2_df["mouse_gene_symbol"])))
print(f"P2_Clean_SexTrimmed unique genes: {len(p2_genes)}")

# P3: Clean + SexIEGTrimmed
ieg_list = {"Fos", "Fosb", "Jun", "Junb", "Jund", "Egr1", "Atf3", "Dusp1", "Hspa1a", "Hspa1b", "Ier2", "Btg2", "Zfp36"}
p3_df = p2_df[~p2_df["mouse_gene_symbol"].isin(ieg_list)].copy()
p3_genes = sorted(list(set(p3_df["mouse_gene_symbol"])))
print(f"P3_Clean_SexIEGTrimmed unique genes: {len(p3_genes)}")

# Build hierarchy table
hierarchy_rows = []
for tier_name, tier_df, tier_genes in [
    ("P0_Raw257", p0_raw_df.rename(columns={"gene": "mouse_gene_symbol"}), p0_genes),
    ("P1_Clean", p1_df, p1_genes),
    ("P2_Clean_SexTrimmed", p2_df, p2_genes),
    ("P3_Clean_SexIEGTrimmed", p3_df, p3_genes)
]:
    for g in tier_genes:
        sub_g = tier_df[tier_df["mouse_gene_symbol"] == g]
        units_list = sorted(list(set(sub_g["unit"])))
        mean_lfc = sub_g["log2FC"].mean() if "log2FC" in sub_g.columns else 0.0
        direction = "UP" if mean_lfc > 0 else "DOWN"
        hierarchy_rows.append({
            "tier": tier_name,
            "gene": g,
            "direction": direction,
            "mean_log2FC": mean_lfc,
            "n_units_de": len(units_list),
            "units_list": ";".join(units_list),
            "biotype": biotypes.get(g, "unknown"),
            "chromosome": chromosomes.get(g, "unknown")
        })

hierarchy_df = pd.DataFrame(hierarchy_rows)
hierarchy_path = RESULTS_DIR / "m5-program-hierarchy.csv"
hierarchy_df.to_csv(hierarchy_path, index=False)
print(f"Saved {hierarchy_path} ({len(hierarchy_df)} rows across 4 tiers)")

# ---------------------------------------------------------
# Step 2: Covariance-Calibrated Q* Heterogeneity Decomposition
# ---------------------------------------------------------
print("\n--- Step 2: Covariance-Calibrated Q* Heterogeneity Decomposition ---")

# Load all 17 unit DE results
de_files = sorted(RESULTS_DIR.glob("m2_de_*.csv"))
dfs = {fp.stem.replace("m2_de_", "").replace("__", "::"): pd.read_csv(fp).set_index("index") for fp in de_files}
units = list(dfs.keys())
common_genes = sorted(list(set.intersection(*[set(d.index) for d in dfs.values()])))
print(f"Aligning across {len(units)} units and {len(common_genes)} common genes")

# Compute untruncated pairwise Spearman rank correlation matrix
rho_mat = np.eye(len(units))
for i, u1 in enumerate(units):
    for j, u2 in enumerate(units):
        if i < j:
            b1 = dfs[u1].loc[common_genes, "log2FC"]
            b2 = dfs[u2].loc[common_genes, "log2FC"]
            r, _ = stats.spearmanr(b1, b2)
            c = 2.0 * np.sin(np.pi * r / 6.0)
            rho_mat[i, j] = c
            rho_mat[j, i] = c

# Eigenvalue floor at 1e-4
evals, evecs = np.linalg.eigh(rho_mat)
evals_floored = np.maximum(evals, 1e-4)
C_hat = evecs @ np.diag(evals_floored) @ evecs.T
# Normalize diagonal to 1.0
inv_d = 1.0 / np.sqrt(np.diag(C_hat))
C_hat = np.diag(inv_d) @ C_hat @ np.diag(inv_d)

# Organ indicators (Lung vs Kidney)
organs = [u.split("::")[0] for u in units]
X1 = np.zeros((len(units), 2))
for i, o in enumerate(organs):
    if o == "Lung":
        X1[i, 0] = 1.0
    elif o == "Kidney":
        X1[i, 1] = 1.0
contrast_c = np.array([1.0, -1.0])

# MAD-SD scale per unit for empirical-null scaling
mad_sd = {}
for u in units:
    d = dfs[u]
    z = (d["log2FC"] / d["lfcSE"]).dropna()
    s = stats.median_abs_deviation(z, scale="normal")
    mad_sd[u] = min(1.0, float(s))

q_results = []
for g in common_genes:
    b = np.array([dfs[u].loc[g, "log2FC"] for u in units])
    se_raw = np.array([dfs[u].loc[g, "lfcSE"] for u in units])
    se_emp = np.array([dfs[u].loc[g, "lfcSE"] * mad_sd[u] for u in units])
    
    # 1. Conservative Q*_corr (raw SE)
    Sigma_raw = np.diag(se_raw) @ C_hat @ np.diag(se_raw)
    Omega_raw = np.linalg.pinv(Sigma_raw)
    ones = np.ones(len(units))
    mu0_raw = (ones @ Omega_raw @ b) / (ones @ Omega_raw @ ones)
    r0_raw = b - mu0_raw * ones
    q_tot_corr = float(r0_raw @ Omega_raw @ r0_raw)
    
    # 2. Empirical-null Q*_empnull (scaled SE)
    Sigma_emp = np.diag(se_emp) @ C_hat @ np.diag(se_emp)
    Omega_emp = np.linalg.pinv(Sigma_emp)
    mu0_emp = (ones @ Omega_emp @ b) / (ones @ Omega_emp @ ones)
    r0_emp = b - mu0_emp * ones
    q_tot_emp = float(r0_emp @ Omega_emp @ r0_emp)
    
    # Exact nested GLS organ decomposition on Q*_empnull
    Xt_Omega = X1.T @ Omega_emp
    M1 = np.linalg.inv(Xt_Omega @ X1)
    mu1 = M1 @ (Xt_Omega @ b)
    r1 = b - X1 @ mu1
    q_within_joint = float(r1 @ Omega_emp @ r1)
    q_between_organ = float((mu1[0] - mu1[1])**2 / (contrast_c @ M1 @ contrast_c))
    
    # p-values
    p_tot_emp = 1.0 - stats.chi2.cdf(q_tot_emp, df=len(units) - 1)
    p_between = 1.0 - stats.chi2.cdf(q_between_organ, df=1)
    p_within = 1.0 - stats.chi2.cdf(q_within_joint, df=len(units) - 2)
    
    q_results.append({
        "gene": g,
        "n_units": len(units),
        "df_total": len(units) - 1,
        "gls_mean_lfc": mu0_emp,
        "lung_mean_lfc": mu1[0],
        "kidney_mean_lfc": mu1[1],
        "q_tot_corr": q_tot_corr,
        "q_tot_empnull": q_tot_emp,
        "q_between_organ": q_between_organ,
        "q_within_joint": q_within_joint,
        "pvalue_q_tot_empnull": p_tot_emp,
        "pvalue_q_between": p_between,
        "pvalue_q_within": p_within
    })

q_df = pd.DataFrame(q_results)
# Multiple testing correction
import statsmodels.stats.multitest as smm
_, qvals_tot, _, _ = smm.multipletests(q_df["pvalue_q_tot_empnull"], method="fdr_bh")
_, qvals_bet, _, _ = smm.multipletests(q_df["pvalue_q_between"], method="fdr_bh")
q_df["fdr_q_tot_empnull"] = qvals_tot
q_df["fdr_q_between"] = qvals_bet

q_path = RESULTS_DIR / "m5-heterogeneity-qstar.csv"
q_df.to_csv(q_path, index=False)
print(f"Saved {q_path} ({len(q_df)} genes)")
print(f"Median Q*_tot / 16: {q_df['q_tot_empnull'].median() / 16.0:.3f} (expect ~ 0.95)")
print(f"Genes with significant organ heterogeneity (FDR < 0.05): {(q_df['fdr_q_between'] < 0.05).sum()}")

# ---------------------------------------------------------
# Step 3: External Aging Program Overlaps
# ---------------------------------------------------------
print("\n--- Step 3: External Aging Program Overlaps ---")

# Load external gene sets
gags330_df = pd.read_csv(REF_DIR / "zhang2021-330-gags.csv")
gags59_df = pd.read_csv(REF_DIR / "zhang2021-59-droplet-gags.csv")
senmayo_df = pd.read_csv(REF_DIR / "saul2022-senmayo-mouse.csv")

ext_sets = {
    "Zhang2021_FACS_330_GAGs": set(gags330_df["mouse_gene_symbol"]),
    "Zhang2021_Droplet_59_GAGs": set(gags59_df["mouse_gene_symbol"]),
    "Saul2022_SenMayo_Mouse_118": set(senmayo_df["mouse_gene_symbol"])
}

# Background universe: all tested genes in M2
bg_genes = set(sex_audit.index)
M_bg = len(bg_genes)
print(f"Background gene universe: {M_bg} genes")

overlap_rows = []
for tier_name, tier_genes in [
    ("P0_Raw257", set(p0_genes)),
    ("P1_Clean", set(p1_genes)),
    ("P2_Clean_SexTrimmed", set(p2_genes)),
    ("P3_Clean_SexIEGTrimmed", set(p3_genes))
]:
    n_tier = len(tier_genes)
    for set_name, ext_set in ext_sets.items():
        ext_in_bg = ext_set & bg_genes
        n_ext = len(ext_in_bg)
        shared = tier_genes & ext_in_bg
        k = len(shared)
        
        # Hypergeometric test: k successes in n_tier draws from M_bg with n_ext successes
        # scipy.stats.hypergeom(M, n, N) where M=M_bg, n=n_ext, N=n_tier
        p_val = stats.hypergeom.sf(k - 1, M_bg, n_ext, n_tier)
        expected = (n_tier * n_ext) / M_bg
        odds_ratio = (k / expected) if expected > 0 else np.nan
        
        overlap_rows.append({
            "tier": tier_name,
            "tier_size": n_tier,
            "external_set": set_name,
            "external_set_size_in_bg": n_ext,
            "n_shared": k,
            "expected_by_chance": round(expected, 2),
            "fold_enrichment": round(odds_ratio, 2),
            "hypergeometric_pvalue": p_val,
            "shared_genes": ";".join(sorted(list(shared)))
        })

overlap_df = pd.DataFrame(overlap_rows)
overlap_path = RESULTS_DIR / "m5-external-program-overlap.csv"
overlap_df.to_csv(overlap_path, index=False)
print(f"Saved {overlap_path}")
print(overlap_df[["tier", "external_set", "tier_size", "n_shared", "fold_enrichment", "hypergeometric_pvalue"]].to_string(index=False))

# ---------------------------------------------------------
# Step 4: Figures
# ---------------------------------------------------------
print("\n--- Step 4: Publication Figures ---")

plt.figure(figsize=(8, 5))
sns.scatterplot(
    data=q_df,
    x="q_within_joint",
    y="q_between_organ",
    alpha=0.4,
    color="#1f77b4",
    s=15
)
plt.axhline(stats.chi2.ppf(0.95, df=1), color="red", linestyle="--", label="Organ Heterogeneity p < 0.05 threshold")
plt.xlabel("Joint Within-Organ Heterogeneity (Q*_within, df = 15)")
plt.ylabel("Between-Organ Heterogeneity (Q*_between, df = 1)")
plt.title("GLM Nested Subspace Heterogeneity Decomposition (Lung vs Kidney)")
plt.legend()
plt.tight_layout()
fig_q_path = FIGURES_DIR / "m5-heterogeneity-decomposition.png"
plt.savefig(fig_q_path)
plt.close()
print(f"Saved {fig_q_path}")

print("\n=== M5 HETEROGENEITY PIPELINE FINISHED SUCCESSFULLY ===")
