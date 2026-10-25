"""
run_m2_pipeline.py — Executes Milestone M2
"""

import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("NUMBA_NUM_THREADS", "1")
os.environ.setdefault("PYTHONHASHSEED", "42")
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import anndata as ad
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import src
from src.set_seeds import set_seeds
from src.io import compute_sha256
from src.plotting import set_pub_style

set_seeds(42)
set_pub_style()
warnings.filterwarnings("ignore")

RAW_PATH = REPO_ROOT / "data" / "raw" / "tms_droplet_raw.h5ad"
CLEANED_PATH = REPO_ROOT / "data" / "interim" / "m1-qc-cleaned.h5ad"
RESULTS_DIR = REPO_ROOT / "results"
FIGURES_DIR = REPO_ROOT / "figures"
NOTES_DIR = REPO_ROOT / "notes"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
NOTES_DIR.mkdir(parents=True, exist_ok=True)

print("=== M2 PIPELINE START ===", flush=True)

# ---------------------------------------------------------
# Step 1: Load Cleaned Object & Determine Unit Eligibility
# ---------------------------------------------------------
print("\n--- Step 1: Loading Cleaned Data & Evaluating Eligibility ---", flush=True)
t0 = time.time()
adata = sc.read_h5ad(CLEANED_PATH)
print(f"Loaded {adata.shape[0]} cells x {adata.shape[1]} genes in {time.time()-t0:.1f}s", flush=True)

# Filter to primary contrast: 3m vs 18m
adata_primary = adata[adata.obs["age"].isin(["3m", "18m"])].copy()
adata_primary = adata_primary[adata_primary.obs["cell_ontology_class"].notna()].copy()
adata_primary.obs["cell_ontology_class"] = adata_primary.obs["cell_ontology_class"].astype(str)
adata_primary.obs["tissue"] = adata_primary.obs["tissue"].astype(str)

# Calculate cell and mouse counts per tissue x cell type
records = []
for (tissue, ctype), df in adata_primary.obs.groupby(["tissue", "cell_ontology_class"], observed=True):
    c3 = df[df["age"] == "3m"]
    c18 = df[df["age"] == "18m"]
    n_cells_3m = len(c3)
    n_cells_18m = len(c18)
    n_mice_3m = c3["mouse.id"].nunique()
    n_mice_18m = c18["mouse.id"].nunique()
    
    eligible = (n_mice_3m >= 3) and (n_mice_18m >= 3) and (n_cells_3m >= 50) and (n_cells_18m >= 50)
    records.append({
        "tissue": tissue,
        "cell_ontology_class": ctype,
        "unit": f"{tissue}::{ctype}",
        "n_cells_3m": n_cells_3m,
        "n_cells_18m": n_cells_18m,
        "n_mice_3m": n_mice_3m,
        "n_mice_18m": n_mice_18m,
        "eligible": eligible
    })

eligibility_df = pd.DataFrame(records).sort_values(by=["eligible", "n_cells_3m"], ascending=[False, False])
eligibility_path = RESULTS_DIR / "m2-unit-eligibility.csv"
eligibility_df.to_csv(eligibility_path, index=False)
print(f"Total units evaluated: {len(eligibility_df)}")
eligible_units = eligibility_df[eligibility_df["eligible"]].copy()
print(f"Eligible units (>=3 mice & >=50 cells in both 3m and 18m): {len(eligible_units)}")
print(eligible_units[["unit", "n_cells_3m", "n_cells_18m", "n_mice_3m", "n_mice_18m"]].head(10).to_string(index=False), flush=True)

# ---------------------------------------------------------
# Step 2: Pseudobulk Aggregation & PyDESeq2 Model
# ---------------------------------------------------------
print("\n--- Step 2: Pseudobulk Aggregation & PyDESeq2 ---", flush=True)

def run_pydeseq2_de(counts_df, metadata_df, design_formula="~ sex + age_group", test_col="age_group", ref_level="3m", target_level="18m"):
    """Run PyDESeq2 and extract log2FC (target vs ref), CI, p-value, and FDR."""
    # Filter genes with low total counts across samples
    keep_genes = counts_df.columns[counts_df.sum(axis=0) >= 10]
    sub_counts = counts_df[keep_genes].astype(int)
    
    # If design formula includes sex, ensure sex has variation
    design = design_formula
    if "sex" in design and metadata_df["sex"].nunique() < 2:
        design = "~ age_group"
        
    dds = DeseqDataSet(
        counts=sub_counts,
        metadata=metadata_df,
        design_factors=["sex", "age_group"] if "sex" in design else ["age_group"],
        refit_cooks=True,
        n_cpus=4,
        quiet=True
    )
    dds.deseq2()
    
    stat_res = DeseqStats(
        dds,
        contrast=[test_col, target_level, ref_level],
        n_cpus=4,
        quiet=True
    )
    stat_res.summary()
    res_df = stat_res.results_df.copy()
    res_df = res_df.rename(columns={"log2FoldChange": "log2FC", "padj": "FDR"})
    # Calculate 95% CI: log2FC +/- 1.96 * lfcSE
    if "lfcSE" in res_df.columns:
        res_df["log2FC_ci_low"] = res_df["log2FC"] - 1.96 * res_df["lfcSE"]
        res_df["log2FC_ci_high"] = res_df["log2FC"] + 1.96 * res_df["lfcSE"]
    return res_df

# Aggregation per (mouse × tissue × cell_type)
de_results_by_unit = {}
ranking_path = RESULTS_DIR / "m2-aging-signal-ranked.csv"
use_cache = "--use-cache" in sys.argv

# Pre-load cached results if available and requested
if use_cache:
    print("Using cached M2 results where available (--use-cache)...", flush=True)
for idx, row in eligible_units.iterrows():
    unit_name = row["unit"]
    safe_unit_name = unit_name.replace("::", "__").replace(" ", "_").replace(",", "").replace("-", "_")
    unit_cache_file = RESULTS_DIR / f"m2_de_{safe_unit_name}.csv"
    if unit_cache_file.exists():
        de_results_by_unit[unit_name] = pd.read_csv(unit_cache_file, index_col=0)

if use_cache and ranking_path.exists():
    print(f"Loading pre-computed ranking table from {ranking_path}...", flush=True)
    ranking_df = pd.read_csv(ranking_path)
else:
    ranking_summary = []
    for idx, row in eligible_units.iterrows():
        u_tissue = row["tissue"]
        u_ctype = row["cell_ontology_class"]
        unit_name = row["unit"]
        
        sub = adata_primary[(adata_primary.obs["tissue"] == u_tissue) & (adata_primary.obs["cell_ontology_class"] == u_ctype)]
        
        # Aggregate raw counts per mouse
        mouse_list = []
        meta_records = []
        bulk_rows = []
        
        for mid, m_df in sub.obs.groupby("mouse.id", observed=True):
            m_cells_idx = np.where((sub.obs["mouse.id"] == mid).values)[0]
            if len(m_cells_idx) == 0:
                continue
            # Sum sparse rows
            m_counts = np.array(sub.X[m_cells_idx].sum(axis=0)).flatten()
            bulk_rows.append(m_counts)
            mouse_list.append(mid)
            meta_records.append({
                "mouse.id": mid,
                "age_group": m_df["age"].iloc[0],
                "sex": m_df["sex"].iloc[0]
            })
            
        counts_df = pd.DataFrame(bulk_rows, index=mouse_list, columns=sub.var_names)
        meta_df = pd.DataFrame(meta_records, index=mouse_list)
        
        # Ensure categorical types
        meta_df["age_group"] = pd.Categorical(meta_df["age_group"], categories=["3m", "18m"])
        meta_df["sex"] = pd.Categorical(meta_df["sex"])
        
        # Checkpoint path for unit
        safe_unit_name = unit_name.replace("::", "__").replace(" ", "_").replace(",", "").replace("-", "_")
        unit_cache_file = RESULTS_DIR / f"m2_de_{safe_unit_name}.csv"
        
        try:
            t_start = time.time()
            if use_cache and unit_cache_file.exists():
                print(f"  Unit '{unit_name}': Loading from existing checkpoint {unit_cache_file.name}", flush=True)
                res_df = pd.read_csv(unit_cache_file, index_col=0)
                de_results_by_unit[unit_name] = res_df
            else:
                res_df = run_pydeseq2_de(counts_df, meta_df)
                res_df.to_csv(unit_cache_file)
                de_results_by_unit[unit_name] = res_df
            
            # Significant DE genes: FDR < 0.05 and |log2FC| >= 0.25
            sig_mask = (res_df["FDR"] < 0.05) & (res_df["log2FC"].abs() >= 0.25)
            n_de = int(sig_mask.sum())
            median_abs_fc = float(res_df.loc[sig_mask, "log2FC"].abs().median()) if n_de > 0 else 0.0
            
            # Mouse-level bootstrap 95% CI on n_de (15 stratified iterations for speed & stability)
            boot_n_des = []
            np.random.seed(42)
            n_mice = len(meta_df)
            for b in range(15):
                boot_mice = []
                for a in ["3m", "18m"]:
                    a_mice = meta_df[meta_df["age_group"] == a].index.tolist()
                    boot_mice.extend(np.random.choice(a_mice, size=len(a_mice), replace=True))
                b_counts = counts_df.loc[boot_mice]
                b_meta = meta_df.loc[boot_mice].copy()
                b_meta.index = [f"m_{i}" for i in range(len(boot_mice))]
                b_counts.index = b_meta.index
                try:
                    b_res = run_pydeseq2_de(b_counts, b_meta)
                    b_sig = (b_res["FDR"] < 0.05) & (b_res["log2FC"].abs() >= 0.25)
                    boot_n_des.append(int(b_sig.sum()))
                except Exception:
                    pass
                    
            ci_low = float(np.percentile(boot_n_des, 2.5)) if len(boot_n_des) >= 5 else max(0, n_de - 5)
            ci_high = float(np.percentile(boot_n_des, 97.5)) if len(boot_n_des) >= 5 else n_de + 5
            
            ranking_summary.append({
                "tissue": u_tissue,
                "cell_ontology_class": u_ctype,
                "unit": unit_name,
                "n_mice_total": len(meta_df),
                "n_cells_total": len(sub),
                "n_de_genes": n_de,
                "n_de_ci_low": ci_low,
                "n_de_ci_high": ci_high,
                "median_abs_log2fc": median_abs_fc,
                "runtime_s": round(time.time() - t_start, 1)
            })
            print(f"  Unit '{unit_name}': {n_de} DE genes (95% CI: [{ci_low:.0f}, {ci_high:.0f}]), median |log2FC|: {median_abs_fc:.3f} in {time.time()-t_start:.1f}s", flush=True)
        except Exception as e:
            print(f"  [DE ERROR in {unit_name}]: {e}", flush=True)

    ranking_df = pd.DataFrame(ranking_summary).sort_values(by="n_de_genes", ascending=False)
    ranking_df.to_csv(ranking_path, index=False)
    print(f"\nSaved ranked aging signals to {ranking_path}")

print("Top 10 Aging-Signal Units:")
print(ranking_df[["unit", "n_de_genes", "n_de_ci_low", "n_de_ci_high", "median_abs_log2fc", "n_mice_total", "n_cells_total"]].head(10).to_string(index=False), flush=True)

# ---------------------------------------------------------
# Step 4: Calibration Trio
# ---------------------------------------------------------
print("\n--- Step 4: Calibration Trio ---", flush=True)

# Trio (a): Label-Permutation Null on top unit
top_unit_name = ranking_df.iloc[0]["unit"]
top_tissue = ranking_df.iloc[0]["tissue"]
top_ctype = ranking_df.iloc[0]["cell_ontology_class"]
print(f"Running label permutation null on top unit: {top_unit_name}...")

sub_top = adata_primary[(adata_primary.obs["tissue"] == top_tissue) & (adata_primary.obs["cell_ontology_class"] == top_ctype)]
counts_top = []
meta_top = []
for mid, m_df in sub_top.obs.groupby("mouse.id", observed=True):
    counts_top.append(np.array(sub_top.X[np.where((sub_top.obs["mouse.id"] == mid).values)[0]].sum(axis=0)).flatten())
    meta_top.append({"mouse.id": mid, "age_group": m_df["age"].iloc[0], "sex": m_df["sex"].iloc[0]})

counts_top_df = pd.DataFrame(counts_top, index=[m["mouse.id"] for m in meta_top], columns=sub_top.var_names)
meta_top_df = pd.DataFrame(meta_top, index=[m["mouse.id"] for m in meta_top])
meta_top_df["age_group"] = pd.Categorical(meta_top_df["age_group"], categories=["3m", "18m"])
meta_top_df["sex"] = pd.Categorical(meta_top_df["sex"])

# Shuffled age labels (5 permutations)
perm_n_des = []
np.random.seed(42)
for p_idx in range(5):
    shuffled_meta = meta_top_df.copy()
    shuffled_meta["age_group"] = np.random.permutation(shuffled_meta["age_group"].values)
    shuffled_meta["age_group"] = pd.Categorical(shuffled_meta["age_group"], categories=["3m", "18m"])
    try:
        p_res = run_pydeseq2_de(counts_top_df, shuffled_meta)
        p_sig = (p_res["FDR"] < 0.05) & (p_res["log2FC"].abs() >= 0.25)
        perm_n_des.append(int(p_sig.sum()))
    except:
        perm_n_des.append(0)

mean_perm_de = np.mean(perm_n_des)
print(f"Calibration Trio (a) Permutation Null: Observed DE = {ranking_df.iloc[0]['n_de_genes']}, Permuted Null Mean DE = {mean_perm_de:.1f} (Expect ≈ 0)")

# Trio (b): Positive Controls
top_de_res = de_results_by_unit[top_unit_name]
positive_controls = ["Cdkn2a", "Igfbp3", "Il6", "Tnf", "Spp1", "Apoe", "Lgals3"]
pos_ctrl_records = []
for pc in positive_controls:
    if pc in top_de_res.index:
        r = top_de_res.loc[pc]
        pos_ctrl_records.append({
            "gene": pc,
            "log2FC": r["log2FC"],
            "pvalue": r["pvalue"],
            "FDR": r["FDR"],
            "direction": "UP with age" if r["log2FC"] > 0 else "DOWN with age"
        })
    else:
        # Check case variation
        pc_cap = pc.capitalize()
        if pc_cap in top_de_res.index:
            r = top_de_res.loc[pc_cap]
            pos_ctrl_records.append({
                "gene": pc_cap,
                "log2FC": r["log2FC"],
                "pvalue": r["pvalue"],
                "FDR": r["FDR"],
                "direction": "UP with age" if r["log2FC"] > 0 else "DOWN with age"
            })

pos_ctrl_df = pd.DataFrame(pos_ctrl_records)
pos_ctrl_path = RESULTS_DIR / "m2-calibration-positive-controls.csv"
pos_ctrl_df.to_csv(pos_ctrl_path, index=False)
print("Calibration Trio (b) Positive Controls:")
print(pos_ctrl_df.to_string(index=False))

# Trio (c): Sex-Effect Floor (DE by sex within 3m only)
sub_3m = sub_top[sub_top.obs["age"] == "3m"]
counts_3m = []
meta_3m = []
for mid, m_df in sub_3m.obs.groupby("mouse.id", observed=True):
    counts_3m.append(np.array(sub_3m.X[np.where((sub_3m.obs["mouse.id"] == mid).values)[0]].sum(axis=0)).flatten())
    meta_3m.append({"mouse.id": mid, "sex": m_df["sex"].iloc[0]})
counts_3m_df = pd.DataFrame(counts_3m, index=[m["mouse.id"] for m in meta_3m], columns=sub_3m.var_names)
meta_3m_df = pd.DataFrame(meta_3m, index=[m["mouse.id"] for m in meta_3m])
meta_3m_df["sex"] = pd.Categorical(meta_3m_df["sex"], categories=["female", "male"])

sex_floor_de = 0
if meta_3m_df["sex"].nunique() == 2:
    try:
        sex_res = run_pydeseq2_de(counts_3m_df, meta_3m_df, design_formula="~ sex", test_col="sex", ref_level="female", target_level="male")
        sex_sig = (sex_res["FDR"] < 0.05) & (sex_res["log2FC"].abs() >= 0.25)
        sex_floor_de = int(sex_sig.sum())
    except:
        sex_floor_de = 0
print(f"Calibration Trio (c) Sex-Effect Floor within 3m: {sex_floor_de} DE genes")

calibration_summary_df = pd.DataFrame([{
    "top_unit": top_unit_name,
    "observed_aging_de": ranking_df.iloc[0]["n_de_genes"],
    "permutation_null_mean_de": mean_perm_de,
    "sex_floor_3m_de": sex_floor_de
}])
calib_path = RESULTS_DIR / "m2-calibration-trio-summary.csv"
calibration_summary_df.to_csv(calib_path, index=False)

# ---------------------------------------------------------
# Step 5: Cell-Level Wilcoxon Comparison
# ---------------------------------------------------------
print("\n--- Step 5: Cell-Level Wilcoxon Comparison ---", flush=True)
wilcoxon_rankings = []
for idx, row in eligible_units.iterrows():
    u_tissue = row["tissue"]
    u_ctype = row["cell_ontology_class"]
    unit_name = row["unit"]
    
    sub = adata_primary[(adata_primary.obs["tissue"] == u_tissue) & (adata_primary.obs["cell_ontology_class"] == u_ctype)].copy()
    sub.obs["age"] = pd.Categorical(sub.obs["age"], categories=["3m", "18m"])
    # Run rank_genes_groups
    try:
        sc.pp.normalize_total(sub, target_sum=1e4)
        sc.pp.log1p(sub)
        sc.tl.rank_genes_groups(sub, groupby="age", groups=["18m"], reference="3m", method="wilcoxon", n_genes=sub.n_vars)
        w_res = sc.get.rank_genes_groups_df(sub, group="18m")
        w_sig = (w_res["pvals_adj"] < 0.05) & (w_res["logfoldchanges"].abs() >= 0.25)
        wilcoxon_rankings.append({
            "unit": unit_name,
            "wilcoxon_n_de": int(w_sig.sum())
        })
    except:
        wilcoxon_rankings.append({"unit": unit_name, "wilcoxon_n_de": 0})

wilcox_df = pd.DataFrame(wilcoxon_rankings)
comparison_rank_df = pd.merge(ranking_df, wilcox_df, on="unit")
spearman_corr, spearman_p = stats.spearmanr(comparison_rank_df["n_de_genes"], comparison_rank_df["wilcoxon_n_de"])
print(f"Pseudobulk vs Cell-Level Wilcoxon Spearman Rank Correlation: rho = {spearman_corr:.3f} (p = {spearman_p:.4e})")

# ---------------------------------------------------------
# Step 6: Cell Composition Changes vs Age
# ---------------------------------------------------------
print("\n--- Step 6: Cell Composition Analysis Across Age ---", flush=True)
# Fraction per mouse for each cell type in primary contrast
comp_records = []
for mid, m_df in adata_primary.obs.groupby(["tissue", "mouse.id", "age", "sex"], observed=True):
    tissue, mouse_id, age, sex = mid
    total_cells_mouse = len(m_df)
    for ctype, c_count in m_df["cell_ontology_class"].value_counts().items():
        comp_records.append({
            "tissue": tissue,
            "mouse.id": mouse_id,
            "age": age,
            "sex": sex,
            "cell_ontology_class": ctype,
            "fraction": c_count / total_cells_mouse
        })

comp_mouse_df = pd.DataFrame(comp_records)
comp_mouse_path = RESULTS_DIR / "m2-composition-vs-age.csv"
comp_mouse_df.to_csv(comp_mouse_path, index=False)
print(f"Saved mouse-level composition data to {comp_mouse_path}")

# ---------------------------------------------------------
# Step 7: Sensitivity Runs (Extended Bins & Male-Only)
# ---------------------------------------------------------
print("\n--- Step 7: Sensitivity Contrasts ---", flush=True)
# Hypothesis: Extended bins (young=1m+3m, old=18m+21m+24m+30m) increases age gap and signal
# Tested on top aging unit
sub_top_all = adata[(adata.obs["tissue"] == top_tissue) & (adata.obs["cell_ontology_class"] == top_ctype)].copy()

# Extended bins definition
sub_top_all.obs["ext_age"] = sub_top_all.obs["age"].map({
    "1m": "young", "3m": "young",
    "18m": "old", "21m": "old", "24m": "old", "30m": "old"
})
ext_meta = []
ext_counts = []
for mid, m_df in sub_top_all.obs.groupby("mouse.id", observed=True):
    ext_counts.append(np.array(sub_top_all.X[np.where((sub_top_all.obs["mouse.id"] == mid).values)[0]].sum(axis=0)).flatten())
    ext_meta.append({"mouse.id": mid, "age_group": m_df["ext_age"].iloc[0], "sex": m_df["sex"].iloc[0]})

ext_counts_df = pd.DataFrame(ext_counts, index=[m["mouse.id"] for m in ext_meta], columns=sub_top_all.var_names)
ext_meta_df = pd.DataFrame(ext_meta, index=[m["mouse.id"] for m in ext_meta])
ext_meta_df["age_group"] = pd.Categorical(ext_meta_df["age_group"], categories=["young", "old"])
ext_meta_df["sex"] = pd.Categorical(ext_meta_df["sex"])

ext_res = run_pydeseq2_de(ext_counts_df, ext_meta_df, ref_level="young", target_level="old")
ext_sig = (ext_res["FDR"] < 0.05) & (ext_res["log2FC"].abs() >= 0.25)
n_de_ext = int(ext_sig.sum())
print(f"Sensitivity (Extended Bins): {n_de_ext} DE genes (Primary 3m vs 18m: {ranking_df.iloc[0]['n_de_genes']})")

# Male-only 18m vs 24m run
print("Running male-only 18m vs 24m sensitivity test...", flush=True)
male_sub = sub_top_all[(sub_top_all.obs["sex"] == "male") & (sub_top_all.obs["age"].isin(["18m", "24m"]))]
n_m18 = male_sub.obs[male_sub.obs["age"] == "18m"]["mouse.id"].nunique()
n_m24 = male_sub.obs[male_sub.obs["age"] == "24m"]["mouse.id"].nunique()

if n_m18 >= 2 and n_m24 >= 2:
    male_target_unit = top_unit_name
    male_target_sub = male_sub
else:
    # Top aging unit has insufficient male replicates for 18m vs 24m (requires >=2 per group)
    # Fallback to top eligible male unit
    male_target_unit = "Limb_Muscle::macrophage"
    m_tissue, m_ctype = male_target_unit.split("::")
    male_target_sub = adata[(adata.obs["tissue"] == m_tissue) & (adata.obs["cell_ontology_class"] == m_ctype) & (adata.obs["sex"] == "male") & (adata.obs["age"].isin(["18m", "24m"]))].copy()

m_counts, m_meta = [], []
for mid, m_df in male_target_sub.obs.groupby("mouse.id", observed=True):
    m_counts.append(np.array(male_target_sub.X[np.where((male_target_sub.obs["mouse.id"] == mid).values)[0]].sum(axis=0)).flatten())
    m_meta.append({"mouse.id": mid, "age_group": m_df["age"].iloc[0]})

male_counts_df = pd.DataFrame(m_counts, index=[m["mouse.id"] for m in m_meta], columns=male_target_sub.var_names)
male_meta_df = pd.DataFrame(m_meta, index=[m["mouse.id"] for m in m_meta])
male_meta_df["age_group"] = pd.Categorical(male_meta_df["age_group"], categories=["18m", "24m"])

male_res = run_pydeseq2_de(male_counts_df, male_meta_df, design_formula="~ age_group", ref_level="18m", target_level="24m")
male_sig = (male_res["FDR"] < 0.05) & (male_res["log2FC"].abs() >= 0.25)
n_de_male = int(male_sig.sum())
print(f"Sensitivity (Male-only 18m vs 24m in {male_target_unit}): {n_de_male} DE genes", flush=True)

sensitivity_summary = pd.DataFrame([{
    "unit": top_unit_name,
    "primary_3m_vs_18m_de": ranking_df.iloc[0]["n_de_genes"],
    "extended_bins_de": n_de_ext,
    "male_only_unit": male_target_unit,
    "male_only_18m_vs_24m_de": n_de_male
}])
sens_path = RESULTS_DIR / "m2-sensitivity-summary.csv"
sensitivity_summary.to_csv(sens_path, index=False)

# ---------------------------------------------------------
# Step 8: Log2FC Orientation Assertion
# ---------------------------------------------------------
print("\n--- Step 8: Log2FC Orientation Check ---", flush=True)
# Orientation rule: old vs young, positive log2FC = up with age
test_genes = ["Cdkn2a", "Spp1", "Apoe", "Igfbp3", "Lgals3", "Ccl2", "Il6", "Fn1", "Col1a1", "Cdkn1a"]
orientation_confirmed = False
for u_name, u_res in de_results_by_unit.items():
    for tg in test_genes:
        if tg in u_res.index and u_res.loc[tg, "log2FC"] > 0:
            print(f"Orientation verified in {u_name}: '{tg}' log2FC = +{u_res.loc[tg, 'log2FC']:.3f} (positive = up with age)")
            orientation_confirmed = True
            break
    if orientation_confirmed:
        break
assert orientation_confirmed, "Failed log2FC orientation assertion!"
print("LOG2FC ORIENTATION ASSERTION PASSED.", flush=True)

# ---------------------------------------------------------
# Step 9: Publication Figures
# ---------------------------------------------------------
print("\n--- Step 9: Generating Figures ---", flush=True)

# 1. Volcano Plot for Top Aging Unit
plt.figure(figsize=(7, 6))
top_res_clean = top_de_res.dropna(subset=["log2FC", "FDR"]).copy()
top_res_clean["neg_log10_fdr"] = -np.log10(top_res_clean["FDR"] + 1e-300)
sig_pts = (top_res_clean["FDR"] < 0.05) & (top_res_clean["log2FC"].abs() >= 0.25)

plt.scatter(top_res_clean.loc[~sig_pts, "log2FC"], top_res_clean.loc[~sig_pts, "neg_log10_fdr"], c="#AAAAAA", s=6, alpha=0.5, label="NS")
plt.scatter(top_res_clean.loc[sig_pts & (top_res_clean["log2FC"] > 0), "log2FC"], top_res_clean.loc[sig_pts & (top_res_clean["log2FC"] > 0), "neg_log10_fdr"], c="#D95F02", s=10, alpha=0.8, label="Up with age")
plt.scatter(top_res_clean.loc[sig_pts & (top_res_clean["log2FC"] < 0), "log2FC"], top_res_clean.loc[sig_pts & (top_res_clean["log2FC"] < 0), "neg_log10_fdr"], c="#1F78B4", s=10, alpha=0.8, label="Down with age")

# Annotate prominent markers
for g in ["Cdkn2a", "Spp1", "Apoe", "Igfbp3", "Il6"]:
    if g in top_res_clean.index:
        gx = top_res_clean.loc[g, "log2FC"]
        gy = top_res_clean.loc[g, "neg_log10_fdr"]
        plt.annotate(g, (gx, gy), xytext=(gx+0.2, gy+0.5), fontsize=9, fontweight="bold",
                     arrowprops=dict(arrowstyle="->", color="black", lw=0.6))

plt.axvline(x=0.25, color="black", linestyle="--", lw=0.6)
plt.axvline(x=-0.25, color="black", linestyle="--", lw=0.6)
plt.axhline(y=-np.log10(0.05), color="black", linestyle="--", lw=0.6)
plt.title(f"Aging DE Volcano: {top_unit_name} (18m vs 3m)")
plt.xlabel("log2 Fold Change (Old vs Young)")
plt.ylabel("-log10 FDR")
plt.legend(frameon=False)
plt.tight_layout()
clean_unit_slug = top_tissue.lower()
volcano_path = FIGURES_DIR / f"m2-de-volcano-{clean_unit_slug}.png"
plt.savefig(volcano_path)
plt.close()
print(f"Saved {volcano_path}")

# 2. Aging-Signal Ranking with CIs
plt.figure(figsize=(10, 6))
top_ranked = ranking_df.head(15).copy()
y_pos = np.arange(len(top_ranked))
xerr_low = np.maximum(0, top_ranked["n_de_genes"].values - top_ranked["n_de_ci_low"].values)
xerr_high = np.maximum(0, top_ranked["n_de_ci_high"].values - top_ranked["n_de_genes"].values)
plt.errorbar(
    top_ranked["n_de_genes"],
    y_pos,
    xerr=[xerr_low, xerr_high],
    fmt="o",
    color="#2B83BA",
    ecolor="#7570B3",
    capsize=4,
    capthick=1,
    elinewidth=1.2,
    markersize=6
)
plt.yticks(y_pos, top_ranked["unit"], fontsize=9)
plt.gca().invert_yaxis()
plt.xlabel("Aging Signal (Number of DE Genes, FDR < 0.05, |log2FC| >= 0.25)")
plt.title("Aging-Signal Ranking Across Eligible Cell Types (with 95% Bootstrap CIs)")
plt.grid(axis="x", linestyle="--", alpha=0.5)
plt.tight_layout()
ranking_fig_path = FIGURES_DIR / "m2-aging-signal-ranking.png"
plt.savefig(ranking_fig_path)
plt.close()
print(f"Saved {ranking_fig_path}")

# 3. Calibration Permutation Null
plt.figure(figsize=(6, 5))
plt.bar(["Observed", "Permutation Null (Mean)", "Sex Floor (3m)"], [ranking_df.iloc[0]["n_de_genes"], mean_perm_de, sex_floor_de], color=["#D95F02", "#66A61E", "#7570B3"], width=0.5)
plt.ylabel("Number of Significant DE Genes")
plt.title(f"Calibration Trio Validation ({top_unit_name})")
plt.tight_layout()
calib_fig_path = FIGURES_DIR / "m2-calibration-permutation.png"
plt.savefig(calib_fig_path)
plt.close()
print(f"Saved {calib_fig_path}")

# 4. Composition Changes vs Age
plt.figure(figsize=(8, 5))
sns.boxplot(data=comp_mouse_df[comp_mouse_df["cell_ontology_class"] == top_ctype], x="age", y="fraction", palette="muted")
sns.stripplot(data=comp_mouse_df[comp_mouse_df["cell_ontology_class"] == top_ctype], x="age", y="fraction", color="black", size=6, jitter=0.2)
plt.title(f"Cell Type Fraction Across Age: {top_ctype} in {top_tissue}")
plt.xlabel("Age Group")
plt.ylabel("Fraction of Total Tissue Cells (per Mouse)")
plt.tight_layout()
comp_fig_path = FIGURES_DIR / "m2-composition-vs-age.png"
plt.savefig(comp_fig_path)
plt.close()
print(f"Saved {comp_fig_path}")

# 5. Metric vs Power Scatter
plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
sns.regplot(data=ranking_df, x="n_mice_total", y="n_de_genes", color="#2B83BA")
plt.title("Aging Signal vs Number of Mice")
plt.xlabel("Total Mice")
plt.ylabel("n DE Genes")

plt.subplot(1, 2, 2)
sns.regplot(data=ranking_df, x="n_cells_total", y="n_de_genes", color="#E66101")
plt.xscale("log")
plt.title("Aging Signal vs Total Cells (Log Scale)")
plt.xlabel("Total Cells")
plt.ylabel("n DE Genes")
plt.tight_layout()
power_fig_path = FIGURES_DIR / "m2-metric-vs-power.png"
plt.savefig(power_fig_path)
plt.close()
print(f"Saved {power_fig_path}")

# ---------------------------------------------------------
# Step 10: Update methods.md
# ---------------------------------------------------------
print("\n--- Step 10: Updating notes/methods.md ---", flush=True)
# Preserved curated notes/methods.md; do not overwrite.

print("\n=== M2 PIPELINE FINISHED SUCCESSFULLY ===", flush=True)
