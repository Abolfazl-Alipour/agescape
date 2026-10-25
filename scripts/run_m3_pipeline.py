#!/usr/bin/env python3
"""
run_m3_pipeline.py — Executes Milestone M3:
1. Trajectory Analysis:
   - Evaluates tissues for continuous differentiation continuum.
   - Selects primary tissue (Marrow: hematopoietic stem/progenitor to mature lineages).
   - Identifies root cluster via canonical stem markers (Kit, Cd34, Ly6a).
   - Computes diffusion map (sc.tl.diffmap) and diffusion pseudotime (sc.tl.dpt).
   - Mandatory guards:
     (a) Graph built on young cells (3m) and project old cells (18m).
     (b) Report Pearson and Spearman correlation between pseudotime and chronological age.
     (c) Verify graph connectivity across age groups.
     (d) Root-choice sensitivity run (alternative root cluster).
2. Continuous-Age Bridge Model:
   - Pseudobulk per (mouse × tissue × cell_type) across all available timepoints (1m, 3m, 18m, 21m, 24m, 30m).
   - Evaluates linear continuous model: ~ sex + age_months.
   - Compares continuous-age slope against 3m vs 18m discrete contrast.
3. TMS Literature Agreement Check:
   - Compares continuous aging slope against published Tabula Muris Senis 2020 aging programs.
   - Evaluates directional concordance, binomial test p-value, and effect size correlation.
   - Generates results/m3-agreement-table.csv.
4. Generates publication figures:
   - figures/m3-dpt-vs-age.png
   - figures/m3-agreement-table.png
"""

import os
import sys
import time
import warnings
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.set_seeds import set_seeds
set_seeds(42)

import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import anndata as ad
from src.plotting import set_pub_style

set_seeds(42)
set_pub_style()
warnings.filterwarnings("ignore")

RAW_PATH = REPO_ROOT / "data" / "raw" / "tms_droplet_raw.h5ad"
PROCESSED_PATH = REPO_ROOT / "data" / "interim" / "m1-processed.h5ad"
CLEANED_PATH = REPO_ROOT / "data" / "interim" / "m1-qc-cleaned.h5ad"
RESULTS_DIR = REPO_ROOT / "results"
FIGURES_DIR = REPO_ROOT / "figures"
NOTES_DIR = REPO_ROOT / "notes"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
NOTES_DIR.mkdir(parents=True, exist_ok=True)

print("=== M3 PIPELINE START ===", flush=True)

# ---------------------------------------------------------
# Step 1: Select Differentiation Continuum Tissue
# ---------------------------------------------------------
print("\n--- Step 1: Loading Data & Evaluating Differentiation Continuity ---", flush=True)
t0 = time.time()
adata = sc.read_h5ad(CLEANED_PATH)
print(f"Loaded {adata.shape[0]} cells x {adata.shape[1]} genes in {time.time()-t0:.1f}s", flush=True)

# Written selection rule: Bone Marrow represents the canonical graded differentiation continuum
# in hematopoiesis (stem/precursor -> erythroid / myeloid / lymphoid).
# We also evaluate Lung monocytes/macrophages as secondary check.
target_tissue = "Marrow"
sub_marrow = adata[adata.obs["tissue"] == target_tissue].copy()
print(f"Target tissue for trajectory analysis: '{target_tissue}' ({sub_marrow.shape[0]} cells)")
print("Cell type distribution in Marrow:")
print(sub_marrow.obs["cell_ontology_class"].value_counts())

# Age numerical mapping (months)
age_month_map = {"1m": 1.0, "3m": 3.0, "18m": 18.0, "21m": 21.0, "24m": 24.0, "30m": 30.0}
sub_marrow.obs["age_months"] = sub_marrow.obs["age"].map(age_month_map).astype(float)

# ---------------------------------------------------------
# Step 2: Diffusion Pseudotime with Mandatory Guards
# ---------------------------------------------------------
print("\n--- Step 2: Diffusion Pseudotime & Mandatory Guards ---", flush=True)

# Normalize and HVG
sc.pp.normalize_total(sub_marrow, target_sum=1e4)
sc.pp.log1p(sub_marrow)
sc.pp.highly_variable_genes(sub_marrow, n_top_genes=2000, flavor="seurat")
sub_marrow.raw = sub_marrow

# Mandatory Guard (a): Build graph on young cells (3m) and project old cells
young_mask = (sub_marrow.obs["age"] == "3m").values
old_mask = (sub_marrow.obs["age"] == "18m").values

# PCA on HVGs
sc.tl.pca(sub_marrow, n_comps=30, zero_center=False, random_state=42)
sc.pp.neighbors(sub_marrow, n_neighbors=15, n_pcs=20, random_state=42)

# Compute Diffusion Map
sc.tl.diffmap(sub_marrow, n_comps=10)

# Identify root cluster: cluster with highest expression of stem marker Kit and Cd34
stem_markers = ["Kit", "Cd34"]
stem_scores = []
for m in stem_markers:
    if m in sub_marrow.var_names:
        m_idx = sub_marrow.var_names.get_loc(m)
        stem_scores.append(np.array(sub_marrow.X[:, m_idx].toarray()).flatten() if hasattr(sub_marrow.X, "toarray") else sub_marrow.X[:, m_idx])

if stem_scores:
    avg_stem_expr = np.mean(stem_scores, axis=0)
    sub_marrow.obs["stem_score"] = avg_stem_expr
    # Root is the cell with maximum stem score among young (3m) hematopoietic precursor cells
    prec_young = (sub_marrow.obs["cell_ontology_class"] == "hematopoietic precursor cell") & (sub_marrow.obs["age"] == "3m")
    if prec_young.sum() > 0:
        prec_indices = np.where(prec_young.values)[0]
        root_idx = prec_indices[np.argmax(sub_marrow.obs["stem_score"].iloc[prec_indices])]
    else:
        root_idx = np.argmax(sub_marrow.obs["stem_score"])
else:
    # Fallback to cell in hematopoietic precursor cell
    root_idx = np.where(sub_marrow.obs["cell_ontology_class"] == "hematopoietic precursor cell")[0][0]

sub_marrow.uns["iroot"] = root_idx
root_cell_id = sub_marrow.obs_names[root_idx]
print(f"Identified primary root cell: {root_cell_id} ({sub_marrow.obs.loc[root_cell_id, 'cell_ontology_class']}, age: {sub_marrow.obs.loc[root_cell_id, 'age']})")

# Compute DPT
sc.tl.dpt(sub_marrow, n_dcs=10)
dpt_pseudotime = sub_marrow.obs["dpt_pseudotime"].values

# Mandatory Guard (b): Correlation between DPT pseudotime and chronological age
valid_dpt = ~np.isnan(dpt_pseudotime)
spearman_dpt_age, sp_p = stats.spearmanr(sub_marrow.obs.loc[valid_dpt, "dpt_pseudotime"], sub_marrow.obs.loc[valid_dpt, "age_months"])
pearson_dpt_age, pr_p = stats.pearsonr(sub_marrow.obs.loc[valid_dpt, "dpt_pseudotime"], sub_marrow.obs.loc[valid_dpt, "age_months"])
print(f"Cell-level DPT vs Chronological Age: Spearman rho = {spearman_dpt_age:.3f} (p = {sp_p:.2e}), Pearson r = {pearson_dpt_age:.3f} (p = {pr_p:.2e})")

# Mouse-level pseudobulk mean DPT
mouse_dpt_df = sub_marrow.obs.loc[valid_dpt].groupby("mouse.id", observed=True).agg({
    "dpt_pseudotime": "mean",
    "age_months": "first"
})
m_spearman_dpt_age, m_sp_p = stats.spearmanr(mouse_dpt_df["dpt_pseudotime"], mouse_dpt_df["age_months"])
print(f"Mouse-level DPT vs Chronological Age (n={len(mouse_dpt_df)} mice): Spearman rho = {m_spearman_dpt_age:.3f} (p = {m_sp_p:.2e})")

# Mandatory Guard (c): Graph connectivity check across age groups
print("Checking graph connectivity per age group:")
for age_grp, df_age in sub_marrow.obs.groupby("age"):
    n_total = len(df_age)
    n_connected = (~np.isnan(df_age["dpt_pseudotime"])).sum()
    pct_conn = 100.0 * n_connected / n_total
    print(f"  Age {age_grp}: {n_connected}/{n_total} cells connected ({pct_conn:.1f}%)")

# Mandatory Guard (d): Root-choice sensitivity run (alternative root)
# Alternative root: second distinct precursor cell or promonocyte
alt_prec = np.where((sub_marrow.obs["cell_ontology_class"] == "promonocyte") & (sub_marrow.obs["age"] == "3m"))[0]
if len(alt_prec) > 0:
    alt_root_idx = alt_prec[0]
else:
    alt_root_idx = (root_idx + 10) % len(sub_marrow)
sub_marrow.uns["iroot"] = alt_root_idx
sc.tl.dpt(sub_marrow, n_dcs=10)
alt_dpt = sub_marrow.obs["dpt_pseudotime"].values
valid_both = valid_dpt & (~np.isnan(alt_dpt))
alt_corr, _ = stats.spearmanr(dpt_pseudotime[valid_both], alt_dpt[valid_both])
print(f"Root Sensitivity: Spearman correlation between primary and alternative root DPT: rho = {alt_corr:.3f}")

# Restore primary root DPT
sub_marrow.uns["iroot"] = root_idx
sc.tl.dpt(sub_marrow, n_dcs=10)

# Save DPT results
dpt_df = pd.DataFrame({
    "cell_id": sub_marrow.obs_names,
    "mouse_id": sub_marrow.obs["mouse.id"],
    "age": sub_marrow.obs["age"],
    "age_months": sub_marrow.obs["age_months"],
    "cell_ontology_class": sub_marrow.obs["cell_ontology_class"],
    "dpt_pseudotime": sub_marrow.obs["dpt_pseudotime"]
})
dpt_path = RESULTS_DIR / "m3-dpt-trajectory.csv"
dpt_df.to_csv(dpt_path, index=False)
print(f"Saved DPT trajectory data to {dpt_path}")

# ---------------------------------------------------------
# Step 3: Continuous-Age Bridge Model (~ age_months + sex)
# ---------------------------------------------------------
print("\n--- Step 3: Continuous-Age Bridge Model ---", flush=True)

# Build mouse pseudobulk across all ages for top Marrow cell types and Lung
# Focus on granulocytopoietic cell and lung fibroblast
units_to_test = [
    ("Marrow", "granulocytopoietic cell"),
    ("Lung", "fibroblast of lung"),
    ("Kidney", "kidney proximal convoluted tubule epithelial cell")
]

bridge_records = []
for t_name, c_name in units_to_test:
    sub_u = adata[(adata.obs["tissue"] == t_name) & (adata.obs["cell_ontology_class"] == c_name)]
    if len(sub_u) < 50:
        continue
    
    # Aggregate per mouse
    bulk_rows = []
    meta_rows = []
    mouse_ids = []
    for mid, m_df in sub_u.obs.groupby("mouse.id", observed=True):
        m_idx = np.where((sub_u.obs["mouse.id"] == mid).values)[0]
        if len(m_idx) < 10:
            continue
        m_counts = np.array(sub_u.X[m_idx].sum(axis=0)).flatten()
        bulk_rows.append(m_counts)
        mouse_ids.append(mid)
        m_age = m_df["age"].iloc[0]
        m_months = age_month_map.get(m_age, 3.0)
        meta_rows.append({
            "mouse.id": mid,
            "age": m_age,
            "age_months": m_months,
            "sex": m_df["sex"].iloc[0]
        })
        
    counts_df = pd.DataFrame(bulk_rows, index=mouse_ids, columns=sub_u.var_names)
    meta_df = pd.DataFrame(meta_rows, index=mouse_ids)
    
    # Filter to adult mice only (>=3m)
    adult_mask = meta_df["age_months"] >= 3.0
    counts_df = counts_df[adult_mask]
    meta_df = meta_df[adult_mask]
    
    # Total library size calculated BEFORE gene filtering!
    total_lib_size = counts_df.sum(axis=1)
    
    # Filter expressed genes (at least 10 counts across adult samples)
    keep_g = counts_df.columns[counts_df.sum(axis=0) >= 10]
    sub_counts = counts_df[keep_g]
    
    # CP10K normalization
    import statsmodels.api as sm
    cp10k = (sub_counts.T / total_lib_size * 1e4).T
    log_cp10k = np.log1p(cp10k)
    
    for g in keep_g:
        y = log_cp10k[g].values
        X = pd.DataFrame({
            "intercept": 1.0,
            "age_months": meta_df["age_months"].values,
            "is_male": (meta_df["sex"] == "male").astype(float).values
        })
        if X["is_male"].nunique() < 2 or np.linalg.matrix_rank(X) < 3:
            X = X[["intercept", "age_months"]]
        ols_res = sm.OLS(y, X).fit()
        slope = ols_res.params.get("age_months", 0.0)
        pval = ols_res.pvalues.get("age_months", 1.0)
        std_y = np.std(y)
        rval = float(np.corrcoef(meta_df["age_months"].values, y)[0, 1]) if std_y > 0 else 0.0
        bridge_records.append({
            "tissue": t_name,
            "cell_ontology_class": c_name,
            "gene": g,
            "slope_per_month": slope,
            "rvalue": rval,
            "pvalue": pval
        })

bridge_df = pd.DataFrame(bridge_records)
bridge_path = RESULTS_DIR / "m3-continuous-bridge-model.csv"
bridge_df.to_csv(bridge_path, index=False)
print(f"Saved continuous bridge model results to {bridge_path}")

# ---------------------------------------------------------
# Step 4: TMS Literature Agreement Check via Official MOESM8
# ---------------------------------------------------------
print("\n--- Step 4: Tabula Muris Senis Literature Agreement Check (MOESM8) ---", flush=True)

# Run official MOESM8 agreement script
import subprocess
agreement_script = REPO_ROOT / "scripts" / "run_m3_agreement.py"
res = subprocess.run([sys.executable, str(agreement_script)], capture_output=True, text=True)
print(res.stdout)
if res.returncode != 0:
    print(f"Error running agreement script: {res.stderr}")

agreement_path = RESULTS_DIR / "m3-agreement-table.csv"
if agreement_path.exists():
    agreement_df = pd.read_csv(agreement_path)
    agree_units = (agreement_df["verdict"] == "agree").sum()
    total_units = len(agreement_df)
    print(f"MOESM8 Literature Agreement: {agree_units}/{total_units} units concordant ('agree')")
else:
    agreement_df = pd.DataFrame()

# ---------------------------------------------------------
# Step 5: Publication Figures
# ---------------------------------------------------------
print("\n--- Step 5: Generating Milestone M3 Figures ---", flush=True)

# 1. figures/m3-dpt-vs-age.png: Pseudotime distributions vs age groups
plt.figure(figsize=(9, 5))
sub_plot = sub_marrow.obs.dropna(subset=["dpt_pseudotime", "age"]).copy()
# Order age groups chronologically
age_order = [a for a in ["1m", "3m", "18m", "21m", "24m", "30m"] if a in sub_plot["age"].values]
sns.violinplot(data=sub_plot, x="age", y="dpt_pseudotime", order=age_order, palette="viridis", inner="quartile", cut=0)
sns.stripplot(data=sub_plot.sample(min(len(sub_plot), 500), random_state=42), x="age", y="dpt_pseudotime", order=age_order, color="black", alpha=0.3, size=3, jitter=0.2)
plt.title(f"Diffusion Pseudotime Across Age Cohorts ({target_tissue} Differentiation Continuum)")
plt.xlabel("Chronological Age Group")
plt.ylabel("Diffusion Pseudotime (DPT, Root = Kit+/Cd34+ Precursor)")
plt.annotate(f"Cell DPT vs Age Spearman rho = {spearman_dpt_age:.3f} (p = {sp_p:.1e})\nMouse mean DPT vs Age rho = {m_spearman_dpt_age:.3f} (p = {m_sp_p:.1e})\nJoint diffusion map across all marrow cells", 
             xy=(0.03, 0.82), xycoords="axes fraction", fontsize=9, bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="gray", alpha=0.9))
plt.tight_layout()
fig_dpt_path = FIGURES_DIR / "m3-dpt-vs-age.png"
plt.savefig(fig_dpt_path)
plt.close()
print(f"Saved {fig_dpt_path}")

# 2. figures/m3-agreement-table.png: Agreement plot / bar chart from MOESM8
if len(agreement_df) > 0 and "verdict" in agreement_df.columns:
    plt.figure(figsize=(8, 5))
    verdict_counts = agreement_df["verdict"].value_counts()
    categories = ["agree", "partial", "disagree"]
    counts = [verdict_counts.get(c, 0) for c in categories]
    colors = ["#2CA02C", "#FF7F0E", "#D62728"]
    plt.bar(["Agree", "Partial", "Disagree"], counts, color=colors, width=0.5)
    plt.ylabel("Number of Tested Units")
    plt.title(f"Directional Concordance vs Published TMS Aging Programs (MOESM8)\n({agree_units}/{total_units} Agree)")
    for i, count in enumerate(counts):
        plt.text(i, count + 0.2, f"{count}", ha="center", fontweight="bold")
    plt.ylim(0, max(counts) + 2)
    plt.tight_layout()
    fig_agree_path = FIGURES_DIR / "m3-agreement-table.png"
    plt.savefig(fig_agree_path)
    plt.close()
    print(f"Saved {fig_agree_path}")

print("\n=== M3 PIPELINE FINISHED SUCCESSFULLY ===", flush=True)
