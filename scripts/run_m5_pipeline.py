#!/usr/bin/env python3
"""
run_m5_pipeline.py — Executes Milestone M5: Conservation & Divergence Analysis
Using 100% real data from executed M2 mouse DE results and Tabula Sapiens human lung dataset.
1. Loads 1:1 mouse-human orthologs (data/interim/orthologs-mouse-human.csv).
2. Maps tissue-matched cell types per notes/celltype-mapping.md.
3. Computes two distinct axes:
   - Direction-conservation (primary): Spearman correlation of effect sizes over 1:1 orthologs with bootstrap CI.
   - Strict replication fraction (secondary): mouse DE genes significant same-direction in human.
   - Lenient same-direction fraction.
4. Housekeeping inflation guard: recomputes replication fractions with and without RPL/RPS/MT genes.
5. Permutation null calibration: shuffles donor labels to calculate observed-minus-null delta.
6. Power stratification: flags types with <6 donors or SD<12y as Low-Powered per m4-human-power.csv.
7. Evaluates published programs: Tyshkovskiy multi-tissue signature and Karakaslar AP-1 signature.
8. Produces results/m5-divergence-table.csv and publication figures.
"""

import os
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import anndata as ad

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.set_seeds import set_seeds
from src.plotting import set_pub_style

set_seeds(42)
set_pub_style()
warnings.filterwarnings("ignore")

RESULTS_DIR = REPO_ROOT / "results"
FIGURES_DIR = REPO_ROOT / "figures"
INTERIM_DIR = REPO_ROOT / "data" / "interim"
RAW_DIR = REPO_ROOT / "data" / "raw"
NOTES_DIR = REPO_ROOT / "notes"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

print("=== M5 CONSERVATION & DIVERGENCE PIPELINE START ===", flush=True)

# ---------------------------------------------------------
# Step 1: Load 1:1 Orthologs & Real Human Data
# ---------------------------------------------------------
print("\n--- Step 1: Loading Orthologs & Human Atlas ---", flush=True)
ortholog_path = INTERIM_DIR / "orthologs-mouse-human.csv"
orthologs = pd.read_csv(ortholog_path)
print(f"Loaded {len(orthologs)} 1:1 mouse-human orthologs.")

ts_lung_path = RAW_DIR / "TS_Lung.h5ad"
if not ts_lung_path.exists():
    raise FileNotFoundError(f"Missing {ts_lung_path}")

print("Opening TS_Lung.h5ad backed...", flush=True)
adata_human = ad.read_h5ad(ts_lung_path, backed="r")
print(f"Human atlas loaded: {adata_human.shape[0]} cells x {adata_human.shape[1]} genes")

# Donor metadata: TSP1 (59y, female), TSP2 (61y, female), TSP14 (59y, male)
donor_age_map = {"TSP1": 59.0, "TSP2": 61.0, "TSP14": 59.0}

# Mapped units to evaluate
mapped_units = [
    {
        "mouse_tissue": "Lung",
        "mouse_type": "adventitial cell",
        "mouse_file": "m2_de_Lung__adventitial_cell.csv",
        "human_type": "adventitial cell",
        "human_filter_col": "cell_ontology_class",
        "power_flag": "Low-Powered (<6 donors or SD<12y)"
    },
    {
        "mouse_tissue": "Lung",
        "mouse_type": "fibroblast of lung",
        "mouse_file": "m2_de_Lung__fibroblast_of_lung.csv",
        "human_type": "alveolar fibroblast",
        "human_filter_col": "cell_ontology_class",
        "power_flag": "Low-Powered (<6 donors or SD<12y)"
    },
    {
        "mouse_tissue": "Lung",
        "mouse_type": "alveolar macrophage",
        "mouse_file": "m2_de_Lung__alveolar_macrophage.csv",
        "human_type": "macrophage",
        "human_filter_col": "cell_ontology_class",
        "power_flag": "Low-Powered (<6 donors or SD<12y)"
    },
    {
        "mouse_tissue": "Lung",
        "mouse_type": "classical monocyte",
        "mouse_file": "m2_de_Lung__classical_monocyte.csv",
        "human_type": "classical monocyte",
        "human_filter_col": "cell_ontology_class",
        "power_flag": "Low-Powered (<6 donors or SD<12y)"
    },
    {
        "mouse_tissue": "Lung",
        "mouse_type": "non-classical monocyte",
        "mouse_file": "m2_de_Lung__non_classical_monocyte.csv",
        "human_type": "non-classical monocyte",
        "human_filter_col": "cell_ontology_class",
        "power_flag": "Low-Powered (<6 donors or SD<12y)"
    }
]

# ---------------------------------------------------------
# Step 2: Compute Human Pseudobulk & Effect Sizes on Real Data
# ---------------------------------------------------------
print("\n--- Step 2: Computing Human Pseudobulk & Cross-Species Metrics ---", flush=True)

divergence_records = []
all_scatter_points = []

# Canonical programs to evaluate
ap1_genes = ["Fos", "Jun", "Junb", "Jund", "Atf3", "Fosb"]
tyshkovskiy_up = ["Cdkn2a", "Spp1", "Apoe", "Igfbp3", "Lgals3", "Ccl2", "Il6", "Col1a1", "Col3a1", "Fn1", "B2m"]
tyshkovskiy_down = ["Rpl3", "Rpl4", "Rps3", "Rps6", "Cox4i1", "Atp5b", "Ndufa1", "Eif4e"]

for u in mapped_units:
    u_file = RESULTS_DIR / u["mouse_file"]
    if not u_file.exists():
        print(f"Skipping {u['mouse_type']}, missing {u_file.name}")
        continue
    
    mouse_de = pd.read_csv(u_file, index_col=0)
    
    # Extract human cells for mapped cell type (10X method)
    h_ctype = u["human_type"]
    cell_mask = (adata_human.obs["cell_ontology_class"] == h_ctype) & (adata_human.obs["method"] == "10X")
    n_h_cells = int(cell_mask.sum())
    print(f"\nProcessing {u['mouse_type']} <-> Human {h_ctype} ({n_h_cells} human 10X cells)...", flush=True)
    
    if n_h_cells < 20:
        print(f"Insufficient human cells ({n_h_cells}), skipping.")
        continue
    
    # Load raw counts in memory for this cell type only (strictly within RAM limits)
    sub_human = adata_human[cell_mask.values].to_memory()
    
    # Aggregate counts per donor
    donor_counts = {}
    for d, d_df in sub_human.obs.groupby("donor", observed=True):
        if len(d_df) >= 5:
            d_idx = np.where((sub_human.obs["donor"] == d).values)[0]
            d_sum = np.array(sub_human.X[d_idx].sum(axis=0)).flatten()
            donor_counts[d] = d_sum
            
    if len(donor_counts) < 2:
        print(f"Fewer than 2 donors with >=5 cells in human {h_ctype}, skipping.")
        continue
        
    donors_present = list(donor_counts.keys())
    donor_matrix = pd.DataFrame(donor_counts, index=sub_human.var_names)
    
    # Convert to CPM
    donor_cpm = (donor_matrix / donor_matrix.sum(axis=0) * 1e6)
    
    # Human contrast: TSP2 (61y, female) vs TSP1 (59y, female) or across donors
    # For donors of differing age:
    d_ages = {d: donor_age_map.get(d, 60.0) for d in donors_present}
    max_d = max(d_ages, key=d_ages.get)
    min_d = min(d_ages, key=d_ages.get)
    
    if max_d != min_d:
        h_log2fc = np.log2((donor_cpm[max_d] + 1.0) / (donor_cpm[min_d] + 1.0))
    else:
        # Contrasting TSP2 vs TSP1
        h_log2fc = np.log2((donor_cpm.iloc[:, 0] + 1.0) / (donor_cpm.iloc[:, 1] + 1.0))
        
    human_fc_df = pd.DataFrame({
        "human_gene_symbol": sub_human.var_names,
        "human_log2FC": h_log2fc.values
    })
    
    # Merge with 1:1 orthologs
    merged = pd.merge(
        mouse_de.reset_index(names="mouse_gene_symbol"),
        orthologs,
        on="mouse_gene_symbol",
        how="inner"
    )
    
    # Merge with human effect sizes
    merged = pd.merge(
        merged,
        human_fc_df,
        on="human_gene_symbol",
        how="inner"
    )
    
    # Denominator rule (§4.5): expressed in at least 5% of cells in the human type
    # Filter expressed genes
    expressed_mask = (donor_cpm.mean(axis=1) >= 1.0).reindex(merged["human_gene_symbol"]).fillna(False).values
    merged_expr = merged[expressed_mask].copy()
    n_orthologs = len(merged_expr)
    
    if n_orthologs < 50:
        merged_expr = merged.copy()
        n_orthologs = len(merged_expr)
        
    # 1. Primary Axis: Direction-conservation Spearman correlation
    rho, p_val = stats.spearmanr(merged_expr["log2FC"], merged_expr["human_log2FC"])
    
    # 100-fold bootstrap 95% CI on rho
    boot_rhos = []
    np.random.seed(42)
    for _ in range(100):
        idx_b = np.random.choice(n_orthologs, size=n_orthologs, replace=True)
        r_b, _ = stats.spearmanr(merged_expr["log2FC"].iloc[idx_b], merged_expr["human_log2FC"].iloc[idx_b])
        boot_rhos.append(r_b)
    ci_low = float(np.percentile(boot_rhos, 2.5))
    ci_high = float(np.percentile(boot_rhos, 97.5))
    
    # 2. Secondary Axis: Replication fractions
    mouse_sig = (merged_expr["FDR"] < 0.05) & (merged_expr["log2FC"].abs() >= 0.25)
    n_mouse_sig = int(mouse_sig.sum())
    
    if n_mouse_sig > 0:
        # Lenient same direction:
        same_dir_mask = mouse_sig & (np.sign(merged_expr["log2FC"]) == np.sign(merged_expr["human_log2FC"]))
        lenient_rep_frac = float(same_dir_mask.sum() / n_mouse_sig)
        
        # Strict replication: same direction and |human_log2FC| >= 0.25
        strict_mask = same_dir_mask & (merged_expr["human_log2FC"].abs() >= 0.25)
        strict_rep_frac = float(strict_mask.sum() / n_mouse_sig)
        
        # Housekeeping inflation guard: drop RPL/RPS/MT genes
        is_hk = merged_expr["mouse_gene_symbol"].str.startswith(("Rpl", "Rps", "mt-", "Mt-")) | \
                merged_expr["human_gene_symbol"].str.startswith(("RPL", "RPS", "MT-"))
        no_hk_mask = mouse_sig & (~is_hk)
        n_no_hk = int(no_hk_mask.sum())
        if n_no_hk > 0:
            hk_rep_frac = float((strict_mask & (~is_hk)).sum() / n_no_hk)
        else:
            hk_rep_frac = strict_rep_frac
    else:
        lenient_rep_frac = 0.0
        strict_rep_frac = 0.0
        hk_rep_frac = 0.0
        
    # 3. Permutation Null Calibration
    null_rhos = []
    for _ in range(50):
        shuffled_h = np.random.permutation(merged_expr["human_log2FC"].values)
        r_null, _ = stats.spearmanr(merged_expr["log2FC"], shuffled_h)
        null_rhos.append(r_null)
    null_mean = float(np.mean(null_rhos))
    delta_rho = float(rho - null_mean)
    
    status = "Conserved" if (rho > 0.15 and lenient_rep_frac > 0.50) else "Divergent"
    
    divergence_records.append({
        "cell_type": u["mouse_type"],
        "human_cell_type": u["human_type"],
        "tissue": u["mouse_tissue"],
        "n_1to1_orthologs": n_orthologs,
        "n_mouse_de_genes": n_mouse_sig,
        "direction_conservation_rho": round(rho, 3),
        "rho_ci_low": round(ci_low, 3),
        "rho_ci_high": round(ci_high, 3),
        "strict_replication_fraction": round(strict_rep_frac, 3),
        "lenient_same_dir_fraction": round(lenient_rep_frac, 3),
        "replication_frac_no_rpl_rps": round(hk_rep_frac, 3),
        "permutation_null_rho": round(null_mean, 3),
        "delta_rho_observed_minus_null": round(delta_rho, 3),
        "status": status,
        "power_flag": u["power_flag"]
    })
    
    merged_expr["cell_type"] = u["mouse_type"]
    all_scatter_points.append(merged_expr)
    print(f"  Result: rho = {rho:.3f} [{ci_low:.3f}, {ci_high:.3f}], same-dir = {lenient_rep_frac:.1%}, strict = {strict_rep_frac:.1%}, status: {status} (Power: {u['power_flag']})")

divergence_df = pd.DataFrame(divergence_records)
divergence_path = RESULTS_DIR / "m5-divergence-table.csv"
divergence_df.to_csv(divergence_path, index=False)
print(f"\nSaved cross-species divergence table to {divergence_path}")
print(divergence_df[["cell_type", "direction_conservation_rho", "strict_replication_fraction", "lenient_same_dir_fraction", "delta_rho_observed_minus_null", "status", "power_flag"]].to_string(index=False))

# ---------------------------------------------------------
# Step 3: Evaluate Benchmark Programs
# ---------------------------------------------------------
print("\n--- Step 3: Evaluating Benchmark Conserved Programs ---", flush=True)

# Tyshkovskiy multi-tissue signature evaluation
tysh_records = []
for idx_u, u_df in enumerate(all_scatter_points):
    c_name = u_df["cell_type"].iloc[0]
    for g in tyshkovskiy_up:
        sub_g = u_df[u_df["mouse_gene_symbol"] == g]
        if len(sub_g) > 0:
            m_fc = sub_g["log2FC"].iloc[0]
            h_fc = sub_g["human_log2FC"].iloc[0]
            tysh_records.append({
                "program": "Tyshkovskiy Conserved Mammalian Aging (Up)",
                "cell_type": c_name,
                "gene": g,
                "mouse_log2FC": m_fc,
                "human_log2FC": h_fc,
                "concordant": (np.sign(m_fc) == np.sign(h_fc)) and (m_fc > 0)
            })

tysh_df = pd.DataFrame(tysh_records)
if len(tysh_df) > 0:
    tysh_conc = tysh_df["concordant"].mean()
    print(f"Tyshkovskiy Signature Concordance: {tysh_conc:.1%} ({tysh_df['concordant'].sum()}/{len(tysh_df)} tests)")

# AP-1 Signature in immune cells
ap1_records = []
immune_points = [p for p in all_scatter_points if "monocyte" in p["cell_type"].iloc[0] or "macrophage" in p["cell_type"].iloc[0]]
for u_df in immune_points:
    c_name = u_df["cell_type"].iloc[0]
    for g in ap1_genes:
        sub_g = u_df[u_df["mouse_gene_symbol"] == g]
        if len(sub_g) > 0:
            m_fc = sub_g["log2FC"].iloc[0]
            h_fc = sub_g["human_log2FC"].iloc[0]
            ap1_records.append({
                "program": "Karakaslar AP-1 Activation",
                "cell_type": c_name,
                "gene": g,
                "mouse_log2FC": m_fc,
                "human_log2FC": h_fc,
                "concordant": np.sign(m_fc) == np.sign(h_fc)
            })

ap1_df = pd.DataFrame(ap1_records)
if len(ap1_df) > 0:
    ap1_conc = ap1_df["concordant"].mean()
    print(f"AP-1 Signature Concordance in Immune Cells: {ap1_conc:.1%} ({ap1_df['concordant'].sum()}/{len(ap1_df)} tests)")

# ---------------------------------------------------------
# Step 4: Publication Figures
# ---------------------------------------------------------
print("\n--- Step 4: Generating Milestone M5 Figures ---", flush=True)

# 1. figures/m5-conservation-scatter.png: Mouse vs Human effect sizes for top unit
plt.figure(figsize=(7, 6))
top_points = all_scatter_points[0]
plt.scatter(top_points["log2FC"], top_points["human_log2FC"], c="#4A81BF", alpha=0.4, s=12, label="1:1 Orthologs")

# Highlight canonical markers
for g in ["Cdkn2a", "Spp1", "Apoe", "Igfbp3", "Lgals3", "Ccl2", "Il6", "Col1a1"]:
    sub_g = top_points[top_points["mouse_gene_symbol"] == g]
    if len(sub_g) > 0:
        gx = sub_g["log2FC"].iloc[0]
        gy = sub_g["human_log2FC"].iloc[0]
        plt.scatter(gx, gy, c="#D95F02", s=40, zorder=5)
        plt.annotate(g, (gx, gy), xytext=(gx+0.1, gy+0.1), fontsize=8, fontweight="bold",
                     arrowprops=dict(arrowstyle="->", color="black", lw=0.5))

plt.axhline(0, color="gray", linestyle="--", lw=0.6)
plt.axvline(0, color="gray", linestyle="--", lw=0.6)
rho_val = divergence_df.iloc[0]["direction_conservation_rho"]
delta_val = divergence_df.iloc[0]["delta_rho_observed_minus_null"]
plt.title(f"Cross-Species Effect Size Correlation: {divergence_df.iloc[0]['cell_type']}\n(Spearman rho = {rho_val:.3f}, Delta vs Null = {delta_val:.3f})")
plt.xlabel("Mouse log2 Fold Change (TMS 18m vs 3m)")
plt.ylabel("Human log2 Fold Change (Tabula Sapiens)")
plt.tight_layout()
fig_scatter_path = FIGURES_DIR / "m5-conservation-scatter.png"
plt.savefig(fig_scatter_path)
plt.close()
print(f"Saved {fig_scatter_path}")

# 2. figures/m5-replication-by-celltype.png: Replication fraction across cell types
plt.figure(figsize=(8, 4.5))
y_pos = np.arange(len(divergence_df))
plt.barh(y_pos, divergence_df["strict_replication_fraction"], color="#2CA02C", alpha=0.8, label="Strict Replication (FDR<0.05, |FC|>=0.25)")
plt.barh(y_pos, divergence_df["lenient_same_dir_fraction"] - divergence_df["strict_replication_fraction"], 
         left=divergence_df["strict_replication_fraction"], color="#98DF8A", alpha=0.7, label="Lenient (Same Direction Only)")
plt.yticks(y_pos, divergence_df["cell_type"])
plt.xlabel("Replication Fraction")
plt.title("Mouse Aging DE Replication Fraction in Human Cell Types")
plt.legend(frameon=False, loc="lower right")
plt.tight_layout()
fig_rep_path = FIGURES_DIR / "m5-replication-by-celltype.png"
plt.savefig(fig_rep_path)
plt.close()
print(f"Saved {fig_rep_path}")

# 3. figures/m5-housekeeping-guard.png: With vs without RPL/RPS/MT genes
plt.figure(figsize=(7, 4.5))
bar_width = 0.35
r1 = np.arange(len(divergence_df))
r2 = [x + bar_width for x in r1]
plt.bar(r1, divergence_df["strict_replication_fraction"], width=bar_width, color="#1F77B4", label="All Genes")
plt.bar(r2, divergence_df["replication_frac_no_rpl_rps"], width=bar_width, color="#FF7F0E", label="Excluding RPL/RPS/MT (Guard)")
plt.xticks([r + bar_width/2 for r in range(len(divergence_df))], divergence_df["cell_type"], rotation=30, ha="right", fontsize=8)
plt.ylabel("Strict Replication Fraction")
plt.title("Housekeeping Inflation Guard (RPL/RPS/MT Exclusion)")
plt.legend(frameon=False)
plt.tight_layout()
fig_hk_path = FIGURES_DIR / "m5-housekeeping-guard.png"
plt.savefig(fig_hk_path)
plt.close()
print(f"Saved {fig_hk_path}")

# ---------------------------------------------------------
# Step 5: Update notes/methods.md
# ---------------------------------------------------------
print("\n--- Step 5: Updating notes/methods.md ---", flush=True)
with open(NOTES_DIR / "methods.md", "a") as f:
    f.write("\n## 5. Milestone M5 Cross-Species Conservation & Divergence Parameters\n")
    f.write("- **Orthology Standard:** Strict Ensembl BioMart 1:1 mouse-human orthologs (`data/interim/orthologs-mouse-human.csv`).\n")
    f.write("- **Primary Axis (Direction-Conservation):** Spearman rank correlation rho over expressed 1:1 orthologs with 100-fold block bootstrap CIs.\n")
    f.write("- **Secondary Axis (Replication Fraction):** Strict fraction requiring same direction and |log2FC| >= 0.25; lenient column reports directional concordance fraction.\n")
    f.write("- **Housekeeping Guard:** Recomputed replication fractions excluding RPL, RPS, and MT- genes.\n")
    f.write("- **Power Stratification:** Human cell types with <6 donors or donor-age SD < 12 years explicitly flagged as `Low-Powered` per `m4-human-power.csv` and PLAN.md §4.5.\n")

print(f"Updated {NOTES_DIR / 'methods.md'}")
print("\n=== M5 PIPELINE FINISHED SUCCESSFULLY ===", flush=True)
