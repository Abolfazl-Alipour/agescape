#!/usr/bin/env python3
"""
scripts/run_h1_epigenetic_drift.py

Phase 3: Empirical Execution of Epigenetic Drift & Universal Epigenetic Aging
across Mammals on the Custom Horvath Mammalian Array (GPL28271, 37k conserved CpGs).

Uses the exact published Universal Clock 2 (Lu et al., Nature Aging 2023)
double-exponential inverse link transformation:
  RelativeAge = (Age + GestationT) / (MaxLifespan + GestationT)
  LinearPredictor (eta) = Intercept + sum(beta_j * CpG_j)
  RelativeAge_pred = exp(-exp(-eta))
  Age_pred = RelativeAge_pred * (MaxLifespan + GestationT) - GestationT

Evaluates:
  1. Human Multi-Tissue Aging (GSE184211, N=661 tissue samples across 8 organs).
  2. Mouse Lifespan Aging (GSE199979, N=339 liver samples from BXD recombinant panel).
  3. Universal Relative Age Prediction.

Outputs:
  - results/h1-epigenetic-clock-human.csv
  - results/h1-epigenetic-clock-mouse.csv
  - results/h1-epigenetic-cross-species-summary.csv
  - figures/h1-pan-mammalian-epigenetic-clock.png
"""

import sys
import gzip
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import statsmodels.api as sm
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

# Published parameters from Lu et al. (Nature Aging, 2023) / AnAge benchmark
MAX_LIFESPAN_HUMAN = 122.5
GESTATION_HUMAN = 0.76  # ~280 days / 365.25 yr

MAX_LIFESPAN_MOUSE = 4.0
GESTATION_MOUSE = 0.05  # ~19 days / 365.25 yr


def trafo_relative_age(age_years, max_lifespan, gestation):
    """Computes bounded relative age incorporating gestation time."""
    return (age_years + gestation) / (max_lifespan + gestation)


def inverse_trafo_age(linear_pred, max_lifespan, gestation):
    """
    Inverse link transformation for Universal Clock 2 (Lu et al. 2023):
    eta = -ln(-ln(RelativeAge))
    RelativeAge_pred = exp(-exp(-eta))
    Age_pred = RelativeAge_pred * (MaxLifespan + GestationT) - GestationT
    """
    pred_rel = np.exp(-np.exp(-linear_pred))
    pred_years = pred_rel * (max_lifespan + gestation) - gestation
    return pred_rel, pred_years


def main():
    print("[H1 Epigenetics] Starting Phase 3: Pan-Mammalian Epigenetic Clock & Drift...")
    results_dir = REPO_ROOT / "results"
    figures_dir = REPO_ROOT / "figures"
    results_dir.mkdir(exist_ok=True)
    figures_dir.mkdir(exist_ok=True)

    # -------------------------------------------------------------
    # 1. Load Universal Pan-Mammalian Epigenetic Clock (Lu et al. 2023)
    # -------------------------------------------------------------
    clock_file = REPO_ROOT / "data" / "reference" / "Universal_Coef.CombAnn.csv"
    if not clock_file.exists():
        raise FileNotFoundError(f"Missing {clock_file}")

    clock_df = pd.read_csv(clock_file, index_col=0)
    intercept = float(clock_df.loc["(Intercept)", "Coef.Universal2_RelativeAge"])
    coefs = clock_df["Coef.Universal2_RelativeAge"].dropna().drop("(Intercept)", errors="ignore")
    print(f"[H1 Epigenetics] Loaded Universal Clock: {len(coefs)} CpGs, Intercept: {intercept:.4f}")

    # -------------------------------------------------------------
    # 2. Human Multi-Tissue Analysis (GSE184211, N=661 tissue samples)
    # -------------------------------------------------------------
    print("\n[H1 Epigenetics] Processing Human GSE184211 (661 tissue samples across 8 tissues)...")
    h_meta = {}
    with gzip.open(REPO_ROOT / "data" / "raw" / "GSE184211_series_matrix.txt.gz", "rt") as f:
        for line in f:
            if line.startswith("!Sample_title"):
                h_meta["title"] = [x.strip('\"') for x in line.strip().split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                h_meta["gsm"] = [x.strip('\"') for x in line.strip().split("\t")[1:]]
            elif line.startswith("!Sample_source_name_ch1"):
                h_meta["tissue"] = [x.strip('\"') for x in line.strip().split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                parts = [x.strip('\"') for x in line.strip().split("\t")[1:]]
                if len(parts) > 0 and ":" in parts[0]:
                    k = parts[0].split(":")[0].strip().lower()
                    h_meta[k] = [p.split(":", 1)[1].strip() if ":" in p else p for p in parts]

    h_df = pd.DataFrame(h_meta)
    h_df["age"] = pd.to_numeric(h_df["age"])
    h_df["sex_bin"] = (h_df["sex"].str.lower() == "male").astype(float)
    h_df["actual_rel_age"] = trafo_relative_age(h_df["age"], MAX_LIFESPAN_HUMAN, GESTATION_HUMAN)

    # Load Human Beta Matrix
    beta_df = pd.read_csv(REPO_ROOT / "data" / "raw" / "GSE184211_datBetaNormalized.csv.gz", index_col=0)
    beta_df.index = beta_df.index.str.strip('\"')
    beta_df.columns = beta_df.columns.str.strip('\"')

    common_h_cpgs = list(set(coefs.index).intersection(set(beta_df.index)))
    print(f"  Matched {len(common_h_cpgs)} / {len(coefs)} Universal CpGs in Human dataset.")

    sub_h_beta = beta_df.loc[common_h_cpgs]
    sub_h_coefs = coefs.loc[common_h_cpgs].values
    linear_pred_h = intercept + np.dot(sub_h_beta.values.T, sub_h_coefs)
    h_df["pred_rel_age"], h_df["pred_age_years"] = inverse_trafo_age(
        linear_pred_h, MAX_LIFESPAN_HUMAN, GESTATION_HUMAN
    )

    # Evaluate Human Performance by Tissue
    h_tissue_rows = []
    tissues = ["All Tissues", "Lung", "Kidney", "Liver", "Heart", "Muscle", "Blood", "Adipose", "Spleen"]

    for t in tissues:
        sub = h_df if t == "All Tissues" else h_df[h_df["tissue"] == t]
        sub = sub.dropna(subset=["actual_rel_age", "pred_rel_age"])
        r, p = pearsonr(sub["actual_rel_age"], sub["pred_rel_age"])
        rho, rho_p = spearmanr(sub["actual_rel_age"], sub["pred_rel_age"])
        mae_rel = np.mean(np.abs(sub["actual_rel_age"] - sub["pred_rel_age"]))
        mae_years = np.mean(np.abs(sub["age"] - sub["pred_age_years"]))

        h_tissue_rows.append({
            "species": "Human",
            "tissue": t,
            "n_samples": len(sub),
            "age_min": sub["age"].min(),
            "age_max": sub["age"].max(),
            "pearson_r": r,
            "pearson_p": p,
            "spearman_rho": rho,
            "r2": r**2,
            "mae_relative_age": mae_rel,
            "mae_years": mae_years
        })

    h_res_df = pd.DataFrame(h_tissue_rows)
    h_res_df.to_csv(results_dir / "h1-epigenetic-clock-human.csv", index=False)
    print(f"[H1 Epigenetics] Saved human tissue results to results/h1-epigenetic-clock-human.csv")
    print(h_res_df[["tissue", "n_samples", "pearson_r", "pearson_p", "r2", "mae_years"]].to_string(index=False))

    # -------------------------------------------------------------
    # 3. Mouse Lifespan Analysis (GSE199979, N=339 Liver, BXD panel)
    # -------------------------------------------------------------
    print("\n[H1 Epigenetics] Processing Mouse GSE199979 (339 liver samples across lifespan)...")
    m_meta = {}
    data_lines = []
    probes = []
    with gzip.open(REPO_ROOT / "data" / "raw" / "GSE199979_series_matrix.txt.gz", "rt") as f:
        for line in f:
            if line.startswith("!Sample_title"):
                m_meta["title"] = [x.strip('\"') for x in line.strip().split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                m_meta["gsm"] = [x.strip('\"') for x in line.strip().split("\t")[1:]]
            elif line.startswith("!Sample_source_name_ch1"):
                m_meta["tissue"] = [x.strip('\"') for x in line.strip().split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                parts = [x.strip('\"') for x in line.strip().split("\t")[1:]]
                if len(parts) > 0 and ":" in parts[0]:
                    k = parts[0].split(":")[0].strip().lower()
                    m_meta[k] = [p.split(":", 1)[1].strip() if ":" in p else p for p in parts]
            elif line.startswith("!series_matrix_table_begin"):
                header = f.readline().strip().split("\t")
                for row_line in f:
                    if row_line.startswith("!series_matrix_table_end"):
                        break
                    parts = row_line.strip().split("\t")
                    cg = parts[0].strip('\"')
                    if cg in coefs.index:
                        probes.append(cg)
                        data_lines.append([float(x) if x not in ("null", "NA", "") else np.nan for x in parts[1:]])
                break

    m_df = pd.DataFrame(m_meta)
    mat = np.array(data_lines)  # probes x samples
    sub_m_coefs = coefs.loc[probes].values
    print(f"  Matched {len(probes)} / {len(coefs)} Universal CpGs in Mouse dataset.")

    linear_pred_m = intercept + np.dot(mat.T, sub_m_coefs)
    m_df["pred_rel_age"], m_df["pred_age_years"] = inverse_trafo_age(
        linear_pred_m, MAX_LIFESPAN_MOUSE, GESTATION_MOUSE
    )

    # Convert age (days) to years and relative age
    m_df["age_days"] = pd.to_numeric(m_df["age (days)"])
    m_df["age_years"] = m_df["age_days"] / 365.25
    m_df["actual_rel_age"] = trafo_relative_age(m_df["age_years"], MAX_LIFESPAN_MOUSE, GESTATION_MOUSE)

    m_valid = m_df.dropna(subset=["actual_rel_age", "pred_rel_age"])
    r_m, p_m = pearsonr(m_valid["actual_rel_age"], m_valid["pred_rel_age"])
    rho_m, rho_p_m = spearmanr(m_valid["actual_rel_age"], m_valid["pred_rel_age"])
    mae_m_rel = np.mean(np.abs(m_valid["actual_rel_age"] - m_valid["pred_rel_age"]))
    mae_m_years = np.mean(np.abs(m_valid["age_years"] - m_valid["pred_age_years"]))

    m_res_df = pd.DataFrame([{
        "species": "Mouse",
        "tissue": "Liver",
        "n_samples": len(m_valid),
        "age_min": m_valid["age_years"].min(),
        "age_max": m_valid["age_years"].max(),
        "pearson_r": r_m,
        "pearson_p": p_m,
        "spearman_rho": rho_m,
        "r2": r_m**2,
        "mae_relative_age": mae_m_rel,
        "mae_years": mae_m_years
    }])
    m_res_df.to_csv(results_dir / "h1-epigenetic-clock-mouse.csv", index=False)
    print(f"[H1 Epigenetics] Saved mouse liver results to results/h1-epigenetic-clock-mouse.csv")
    print(m_res_df[["species", "tissue", "n_samples", "pearson_r", "pearson_p", "r2", "mae_years"]].to_string(index=False))

    # -------------------------------------------------------------
    # 4. Cross-Species Liver-to-Liver Benchmark Summary
    # -------------------------------------------------------------
    cross_df = pd.concat([
        m_res_df,
        h_res_df[h_res_df["tissue"] == "Liver"],
        h_res_df[h_res_df["tissue"] == "Lung"],
        h_res_df[h_res_df["tissue"] == "Kidney"],
        h_res_df[h_res_df["tissue"] == "All Tissues"]
    ], ignore_index=True)
    cross_df.to_csv(results_dir / "h1-epigenetic-cross-species-summary.csv", index=False)
    print(f"\n[H1 Epigenetics] Saved cross-species summary to results/h1-epigenetic-cross-species-summary.csv")

    # -------------------------------------------------------------
    # 5. Generate Publication-Quality Figure
    # -------------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    # Panel A: Mouse Liver (N=339)
    axes[0].scatter(m_valid["actual_rel_age"], m_valid["pred_rel_age"], color="crimson", alpha=0.6, s=40, edgecolors="none")
    lims_a = [0, max(m_valid["actual_rel_age"].max(), m_valid["pred_rel_age"].max()) + 0.05]
    axes[0].plot(lims_a, lims_a, color="black", linestyle="--", linewidth=1.2, label="Identity (y = x)")
    m_fit = np.poly1d(np.polyfit(m_valid["actual_rel_age"], m_valid["pred_rel_age"], 1))
    axes[0].plot(lims_a, m_fit(lims_a), color="darkred", linewidth=1.8, label=f"Fit (r = {r_m:.3f})")
    axes[0].set_xlabel("Actual Relative Age (with gestation)")
    axes[0].set_ylabel("Predicted Relative Age (Universal Clock 2)")
    axes[0].set_title(f"A. Mouse Liver (N = {len(m_valid)})\nPearson r = {r_m:.3f}, p = {p_m:.1e}, MAE = {mae_m_years:.2f} yr ({mae_m_years*12:.1f} mo)")
    axes[0].legend(frameon=True, fontsize=9)

    # Panel B: Human Liver (N=97)
    h_liv = h_df[h_df["tissue"] == "Liver"].dropna(subset=["actual_rel_age", "pred_rel_age"])
    r_hl, p_hl = pearsonr(h_liv["actual_rel_age"], h_liv["pred_rel_age"])
    mae_hl = np.mean(np.abs(h_liv["age"] - h_liv["pred_age_years"]))
    axes[1].scatter(h_liv["actual_rel_age"], h_liv["pred_rel_age"], color="royalblue", alpha=0.7, s=50, edgecolors="black")
    lims_b = [0, max(h_liv["actual_rel_age"].max(), h_liv["pred_rel_age"].max()) + 0.05]
    axes[1].plot(lims_b, lims_b, color="black", linestyle="--", linewidth=1.2, label="Identity (y = x)")
    hl_fit = np.poly1d(np.polyfit(h_liv["actual_rel_age"], h_liv["pred_rel_age"], 1))
    axes[1].plot(lims_b, hl_fit(lims_b), color="darkblue", linewidth=1.8, label=f"Fit (r = {r_hl:.3f})")
    axes[1].set_xlabel("Actual Relative Age (with gestation)")
    axes[1].set_ylabel("Predicted Relative Age (Universal Clock 2)")
    axes[1].set_title(f"B. Human Liver (N = {len(h_liv)})\nPearson r = {r_hl:.3f}, p = {p_hl:.1e}, MAE = {mae_hl:.2f} yr")
    axes[1].legend(frameon=True, fontsize=9)

    # Panel C: Human Multi-Tissue (N=661 across 8 tissues)
    tissue_colors = {
        "Lung": "#e41a1c",
        "Kidney": "#377eb8",
        "Liver": "#4daf4a",
        "Heart": "#984ea3",
        "Muscle": "#ff7f00",
        "Blood": "#a65628",
        "Adipose": "#f781bf",
        "Spleen": "#999999"
    }
    for t_name, color in tissue_colors.items():
        sub_t = h_df[h_df["tissue"] == t_name].dropna(subset=["actual_rel_age", "pred_rel_age"])
        if len(sub_t) > 0:
            axes[2].scatter(sub_t["actual_rel_age"], sub_t["pred_rel_age"], label=f"{t_name} ({len(sub_t)})", color=color, alpha=0.6, s=35)

    lims_c = [0, 0.85]
    axes[2].plot(lims_c, lims_c, color="black", linestyle="--", linewidth=1.2, label="Identity")
    r_all = h_res_df[h_res_df["tissue"] == "All Tissues"].iloc[0]["pearson_r"]
    p_all = h_res_df[h_res_df["tissue"] == "All Tissues"].iloc[0]["pearson_p"]
    mae_all = h_res_df[h_res_df["tissue"] == "All Tissues"].iloc[0]["mae_years"]
    axes[2].set_xlabel("Actual Relative Age (with gestation)")
    axes[2].set_ylabel("Predicted Relative Age (Universal Clock 2)")
    axes[2].set_title(f"C. Human Multi-Tissue Cohort (N = 661)\nPearson r = {r_all:.3f}, p = {p_all:.1e}, MAE = {mae_all:.2f} yr")
    axes[2].legend(frameon=True, fontsize=8, ncol=2)

    plt.tight_layout()
    fig_path = figures_dir / "h1-pan-mammalian-epigenetic-clock.png"
    plt.savefig(fig_path, dpi=200)
    plt.close()
    print(f"[H1 Epigenetics] Saved publication figure to {fig_path}")


if __name__ == "__main__":
    main()
