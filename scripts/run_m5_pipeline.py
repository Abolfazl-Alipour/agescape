#!/usr/bin/env python3
"""
run_m5_pipeline.py — Milestone M5 (HONEST REBUILD): mouse aging programs + the
human-replication reality check.

WHY REBUILT: the first version scored mouse DE against a "human DE" contrast between
TS lung donors aged 59y and 61y (verified: TS Lung new-schema CELLxGENE obs — donors
TSP1/2/14/25 = 59/61/59/60y). That is not an aging contrast. TS v1 is an adult cohort:
TS Lung SD ~1y; TS Skin (4 donors) spans 33-61y but n=4. SenNet public search was
queried 2026-10-14 (search.api.sennetconsortium.org: 6,820 published datasets; no
donor-age fields in the public index; per-dataset traversal + Globus transfer beyond
session budget — documented in notes/methods.md). CONCLUSION: human per-cell-type
age-replication is UNTESTABLE with the human data obtainable here. This pipeline
therefore reports: (1) the real mouse program catalog, (2) housekeeping composition
of mouse DE lists, (3) a divergence table whose status column says exactly what is
and is not testable. No human number is presented as evidence of conservation.
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parent.parent
RESULTS, FIGURES = REPO/"results", REPO/"figures"
FDR_T, FC_T = 0.05, 0.25
HK = lambda g: str(g).startswith(("Rpl","Rps","mt-","Mrp")) or str(g).startswith("MT-")

de_files = sorted(RESULTS.glob("m2_de_*.csv"))
rows, gene_rows = [], []
bg = pd.read_csv(de_files[0])
all_genes_bg = set(bg["index"])
hk_bg = np.mean([HK(g) for g in all_genes_bg])
for fp in de_files:
    unit = fp.stem.replace("m2_de_","").replace("__","::")
    df = pd.read_csv(fp)
    sig = df[(df["FDR"] < FDR_T) & (df["log2FC"].abs() >= FC_T)]
    up = sig[sig["log2FC"] > 0]; dn = sig[sig["log2FC"] < 0]
    hk_frac = np.mean([HK(g) for g in sig["index"]]) if len(sig) else 0.0
    rows.append({"unit": unit, "n_de": len(sig), "n_up": len(up), "n_down": len(dn),
                 "median_abs_log2fc": sig["log2FC"].abs().median() if len(sig) else np.nan,
                 "frac_ribo_mito": round(float(hk_frac), 3),
                 "frac_ribo_mito_background": round(float(hk_bg), 3),
                 "human_replication": "UNTESTABLE (TS adult cohort, no age contrast; SenNet documented)"})
    for _, r in sig.iterrows():
        gene_rows.append({"unit": unit, "gene": r["index"], "log2FC": r["log2FC"],
                          "FDR": r["FDR"], "is_ribo_mito": HK(r["index"])})

cat = pd.DataFrame(rows).sort_values("n_de", ascending=False)
cat.to_csv(RESULTS/"m5-mouse-program-catalog.csv", index=False)
gd = pd.DataFrame(gene_rows)
gd.to_csv(RESULTS/"m5-gene-level-programs.csv", index=False)
cat["status"] = "Catalogued (mouse); human arm untestable -> no conservation/divergence call"
cat[["unit","n_de","n_up","n_down","median_abs_log2fc","frac_ribo_mito","status"]]\
  .to_csv(RESULTS/"m5-divergence-table.csv", index=False)
print(f"[m5] catalogued {len(cat)} units, {len(gd)} program genes -> results/m5-*.csv")

# figures: ribo/mito inflation + top genes recurrence
fig, ax = plt.subplots(figsize=(5.2, 3.4))
x = np.arange(len(cat))
ax.bar(x, cat["frac_ribo_mito"], color="#d62728", label="DE genes (FDR<0.05, |log2FC|>=0.25)")
ax.axhline(hk_bg, color="k", ls="--", lw=1, label=f"background fraction ({hk_bg:.2f})")
ax.set_xticks(x); ax.set_xticklabels(cat["unit"], rotation=75, ha="right", fontsize=6)
ax.set_ylabel("fraction ribosomal / mitochondrial"); ax.legend(fontsize=7)
ax.set_title("Housekeeping composition of mouse aging DE lists")
fig.tight_layout(); fig.savefig(FIGURES/"m5-housekeeping-guard.png"); plt.close(fig)

top = gd["gene"].value_counts().head(15)
fig, ax = plt.subplots(figsize=(5.2, 3.4))
ax.barh(range(len(top))[::-1], top.values, color="#1f77b4")
ax.set_yticks(range(len(top))[::-1]); ax.set_yticklabels(top.index, fontsize=7)
ax.set_xlabel("number of powered units with this gene DE")
ax.set_title("Most recurrent genes across mouse aging programs (17 units)")
fig.tight_layout(); fig.savefig(FIGURES/"m5-program-recurrence.png"); plt.close(fig)
print("[m5] figures written")
