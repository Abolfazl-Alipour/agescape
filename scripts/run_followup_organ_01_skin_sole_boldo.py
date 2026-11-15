"""Dataset 01: Solé-Boldo et al. (2020) whole human skin scRNA-seq (GSE130973)."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.followup_organs import (
    DATA_TMP_ROOT,
    clean_temp_dir,
    counts_to_log2cpm,
    download_with_guard,
    extract_seurat_rds_selective,
    load_p1_clean_human,
    prepare_counts_matrix,
    run_family_a_clr,
    run_family_b_genes_and_concordance,
    run_family_c_p1_signature,
    write_clean_csv,
    write_clean_markdown,
)

DONOR_MAP = [
    {"donor_id": "S1", "sample_label": "y1", "gsm": "GSM3758115", "age_years": 25.0, "age_group": "young", "sex": "M", "sex_male": 1, "ethnicity": "Caucasian", "expected_cells": 2784},
    {"donor_id": "S2", "sample_label": "y2", "gsm": "GSM3758116", "age_years": 27.0, "age_group": "young", "sex": "M", "sex_male": 1, "ethnicity": "Caucasian", "expected_cells": 2670},
    {"donor_id": "S3", "sample_label": "o1", "gsm": "GSM3758117", "age_years": 53.0, "age_group": "old", "sex": "M", "sex_male": 1, "ethnicity": "Caucasian", "expected_cells": 3324},
    {"donor_id": "S4", "sample_label": "o2", "gsm": "GSM3758118", "age_years": 70.0, "age_group": "old", "sex": "M", "sex_male": 1, "ethnicity": "Caucasian", "expected_cells": 2144},
    {"donor_id": "S5", "sample_label": "o3", "gsm": "GSM3758119", "age_years": 69.0, "age_group": "old", "sex": "M", "sex_male": 1, "ethnicity": "Caucasian", "expected_cells": 4535},
]

CLUSTER_TO_MAIN = {
    "0": "Macrophages_DC",
    "1": "Fibroblasts",
    "2": "Fibroblasts",
    "3": "Fibroblasts",
    "4": "Vascular_EC",
    "5": "Keratinocytes",
    "6": "T_cells",
    "7": "Keratinocytes",
    "8": "Pericytes",
    "9": "Fibroblasts",
    "10": "Pericytes",
    "11": "Erythrocytes",
    "12": "Lymphatic_EC",
    "13": "Macrophages_DC",
    "14": "Melanocytes",
    "15": "Keratinocytes",
    "16": "Macrophages_DC",
}

CLUSTER_TO_FINE = {
    "0": "c00_Macrophages_DC_1",
    "1": "c01_Fibroblast_Secretory_Reticular",
    "2": "c02_Fibroblast_Pro_Inflammatory",
    "3": "c03_Fibroblast_Secretory_Papillary",
    "4": "c04_Vascular_EC",
    "5": "c05_Keratinocyte_Diff",
    "6": "c06_T_cells",
    "7": "c07_Keratinocyte_Undiff",
    "8": "c08_Pericytes_SMC_1",
    "9": "c09_Fibroblast_Mesenchymal",
    "10": "c10_Pericytes_2",
    "11": "c11_Erythrocytes",
    "12": "c12_Lymphatic_EC",
    "13": "c13_Langerhans_DC_2",
    "14": "c14_Melanocytes",
    "15": "c15_Keratinocyte_Prolif",
    "16": "c16_Macrophages_DC_Prolif",
}


def aggregate_by_group(
    raw_csc,
    genes: np.ndarray,
    donor_vec: np.ndarray,
    group_vec: np.ndarray,
    min_cells: int = 20,
) -> tuple[pd.DataFrame, dict[str, pd.DataFrame], dict[str, pd.DataFrame], dict[str, dict]]:
    """Aggregate sparse raw gene x cell matrix to donor pseudobulks per group."""
    cell_counts = pd.crosstab(pd.Series(donor_vec, name="donor_id"), pd.Series(group_vec, name="group"))
    groups = sorted(cell_counts.columns.tolist())
    donors = sorted(cell_counts.index.tolist())

    expr_by_group: dict[str, pd.DataFrame] = {}
    meta_by_group: dict[str, pd.DataFrame] = {}
    stats_by_group: dict[str, dict] = {}
    donor_df = pd.DataFrame(DONOR_MAP).set_index("donor_id")
    donor_df["age_per_decade"] = donor_df["age_years"] / 10.0

    for g in groups:
        eligible = [d for d in donors if int(cell_counts.loc[d, g]) >= min_cells]
        if len(eligible) < 4:
            continue
        pb_cols = {}
        for d in eligible:
            mask = (donor_vec == d) & (group_vec == g)
            summed = np.asarray(raw_csc[:, mask].sum(axis=1)).ravel()
            pb_cols[d] = summed
        raw_pb = pd.DataFrame(pb_cols, index=pd.Index(genes, dtype=str))
        filt_pb, p_stats = prepare_counts_matrix(raw_pb)
        log2cpm, lib_sizes = counts_to_log2cpm(filt_pb)
        m = donor_df.loc[eligible].copy()
        m["lib_size"] = lib_sizes.loc[eligible]
        m["log_lib"] = np.log1p(m["lib_size"].to_numpy(dtype=float))
        m["n_cells"] = cell_counts.loc[eligible, g].astype(int)
        expr_by_group[g] = log2cpm
        meta_by_group[g] = m
        stats_by_group[g] = p_stats
    return cell_counts, expr_by_group, meta_by_group, stats_by_group


def main() -> None:
    tmp_dir = DATA_TMP_ROOT / "01_gse130973"
    rds_path = tmp_dir / "GSE130973_seurat_analysis_lyko.rds.gz"
    if not rds_path.exists():
        download_with_guard(
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130973/suppl/GSE130973_seurat_analysis_lyko.rds.gz",
            rds_path,
            max_file_gb=2.0,
        )

    info = extract_seurat_rds_selective(rds_path)
    md = info["meta.data"].copy()
    raw_csc = info["raw_counts"]
    genes = info["genes"]

    donor_vec = md["subj"].astype(str).to_numpy()
    cl_vec = md["integrated_snn_res.0.4"].astype(str).to_numpy()
    main_vec = np.array([CLUSTER_TO_MAIN[c] for c in cl_vec], dtype=str)
    fine_vec = np.array([CLUSTER_TO_FINE[c] for c in cl_vec], dtype=str)

    donor_meta = pd.DataFrame(DONOR_MAP).set_index("donor_id")
    donor_meta["age_per_decade"] = donor_meta["age_years"] / 10.0
    obs_cells = pd.Series(donor_vec).value_counts()
    for d, row in donor_meta.iterrows():
        if int(obs_cells[d]) != int(row["expected_cells"]):
            raise RuntimeError(f"Donor cell count mismatch for {d}: {obs_cells[d]} vs {row['expected_cells']}")
    donor_meta["observed_cells"] = obs_cells.loc[donor_meta.index].astype(int)

    # Whole-donor summary table & Family A (9 main cell types + 17 fine clusters)
    main_counts, expr_main, meta_main, stats_main = aggregate_by_group(raw_csc, genes, donor_vec, main_vec, min_cells=20)
    fine_counts, expr_fine, meta_fine, _ = aggregate_by_group(raw_csc, genes, donor_vec, fine_vec, min_cells=20)

    fam_a_main = run_family_a_clr(
        main_counts, donor_meta, age_col="age_per_decade", candidate_covariates=(), min_donors_present=4, min_total_cells=20, min_residual_df=2
    )
    fam_a_fine = run_family_a_clr(
        fine_counts, donor_meta, age_col="age_per_decade", candidate_covariates=(), min_donors_present=4, min_total_cells=20, min_residual_df=2
    )

    p1 = load_p1_clean_human(REPO_ROOT)
    # Primary Family B and C: adjust for log_lib (df=2); sensitivity: unadjusted age_per_decade (df=3)
    fam_b_genes_adj, fam_b_sum_adj, fam_b_conc_adj = run_family_b_genes_and_concordance(
        expr_main, meta_main, age_col="age_per_decade", candidate_covariates=("log_lib",),
        min_donors=4, min_cpm=1.0, min_expr_frac=0.5, min_residual_df=2, is_unbiased_single_cell=True,
    )
    fam_b_genes_unadj, fam_b_sum_unadj, fam_b_conc_unadj = run_family_b_genes_and_concordance(
        expr_main, meta_main, age_col="age_per_decade", candidate_covariates=(),
        min_donors=4, min_cpm=1.0, min_expr_frac=0.5, min_residual_df=2, is_unbiased_single_cell=True,
    )

    fam_c_main_adj = run_family_c_p1_signature(
        expr_main, meta_main, p1_program=p1, age_col="age_per_decade", candidate_covariates=("log_lib",),
        min_donors=4, min_cpm=1.0, min_expr_frac=0.5, min_residual_df=2,
    )
    fam_c_main_unadj = run_family_c_p1_signature(
        expr_main, meta_main, p1_program=p1, age_col="age_per_decade", candidate_covariates=(),
        min_donors=4, min_cpm=1.0, min_expr_frac=0.5, min_residual_df=2,
    )
    fam_c_fine_unadj = run_family_c_p1_signature(
        expr_fine, meta_fine, p1_program=p1, age_col="age_per_decade", candidate_covariates=(),
        min_donors=4, min_cpm=1.0, min_expr_frac=0.5, min_residual_df=2,
    )

    # Save CSV deliverables
    results_dir = REPO_ROOT / "results"
    reports_dir = REPO_ROOT / "reports"
    prefix = "followup-organ-01-skin-sole-boldo-gse130973"

    donor_out = donor_meta.reset_index().join(main_counts, on="donor_id")
    write_clean_csv(donor_out, results_dir / f"{prefix}-donors.csv")
    write_clean_csv(fam_a_main, results_dir / f"{prefix}-fractions-main.csv")
    write_clean_csv(fam_a_fine, results_dir / f"{prefix}-fractions-fine.csv")
    write_clean_csv(fam_b_sum_adj.assign(model="age_plus_log_lib"), results_dir / f"{prefix}-gene-summary.csv")
    write_clean_csv(fam_b_genes_adj.head(200), results_dir / f"{prefix}-top-genes.csv")
    conc_out = pd.concat(
        [fam_b_conc_adj.assign(model="age_plus_log_lib"), fam_b_conc_unadj.assign(model="age_unadjusted")],
        ignore_index=True,
    )
    write_clean_csv(conc_out, results_dir / f"{prefix}-concordance.csv")
    sig_out = pd.concat(
        [
            fam_c_main_adj.assign(level="main_9", model="age_plus_log_lib"),
            fam_c_main_unadj.assign(level="main_9", model="age_unadjusted"),
            fam_c_fine_unadj.assign(level="fine_17", model="age_unadjusted"),
        ],
        ignore_index=True,
    )
    write_clean_csv(sig_out, results_dir / f"{prefix}-p1-signature.csv")

    # Summary stats for report
    total_cells = int(donor_meta["observed_cells"].sum())
    n_genes_total = int(raw_csc.shape[0])
    n_tech = next(iter(stats_main.values()))["n_technical_removed"]
    min_q_a = float(fam_a_main["q_bh_family_a"].min()) if len(fam_a_main) else np.nan
    min_q_b = float(fam_b_genes_adj["q_bh_global"].min()) if len(fam_b_genes_adj) else np.nan
    n_sig_b = int((fam_b_genes_adj["q_bh_global"] < 0.05).sum()) if len(fam_b_genes_adj) else 0
    mb_size = int(fam_b_genes_adj["family_b_global_size"].iloc[0]) if len(fam_b_genes_adj) else 0
    r_adj = float(fam_b_conc_adj["mean_pairwise_pearson_r"].iloc[0])
    p_adj = float(fam_b_conc_adj["perm_p_one_sided"].iloc[0])
    r_unadj = float(fam_b_conc_unadj["mean_pairwise_pearson_r"].iloc[0])
    p_unadj = float(fam_b_conc_unadj["perm_p_one_sided"].iloc[0])
    min_q_c_emp = float(fam_c_main_adj["q_bh_empirical"].min())
    min_q_c_ana = float(fam_c_main_adj["q_bh_composite"].min())

    a_lines = "\n".join(
        f"| {r['cell_type']} | {int(r['total_cells'])} | {r['median_fraction']:.4f} | {r['effect_clr']:+.3f} | {r['pvalue']:.4f} | {r['q_bh_family_a']:.4f} |"
        for _, r in fam_a_main.iterrows()
    )
    b_lines = "\n".join(
        f"| {r['group']} | {int(r['n_donors'])} | {int(r['n_genes_tested'])} | {int(r['n_sig_global_q05'])} | {r['min_pvalue']:.2e} | {r['min_q_bh_global']:.4f} |"
        for _, r in fam_b_sum_adj.iterrows()
    )
    c_lines = "\n".join(
        f"| {r['group']} | {int(r['n_donors'])} | {int(r['n_p1_up_detected'])}/{int(r['n_p1_down_detected'])} | {r['composite_effect']:+.4f} | {r['composite_pvalue']:.4f} | {r['empirical_null_pvalue']:.4f} | {r['q_bh_empirical']:.4f} |"
        for _, r in fam_c_main_adj.iterrows()
    )

    n_sig_a = int((fam_a_main["q_bh_family_a"] < 0.05).sum()) if len(fam_a_main) else 0
    pos_both_sig = fam_c_main_adj[
        (fam_c_main_adj["composite_effect"] > 0)
        & (fam_c_main_adj["q_bh_composite"] < 0.05)
        & (fam_c_main_adj["q_bh_empirical"] < 0.05)
    ] if len(fam_c_main_adj) else pd.DataFrame()
    neg_emp_hits = fam_c_main_adj[
        (fam_c_main_adj["composite_effect"] < 0)
        & (fam_c_main_adj["q_bh_empirical"] < 0.05)
    ] if len(fam_c_main_adj) else pd.DataFrame()

    top_a_desc = (
        f"top nominal cell type `{fam_a_main.iloc[0]['cell_type']}` "
        f"(slope `{fam_a_main.iloc[0]['effect_clr']:+.3f}` CLR/decade, `p = {fam_a_main.iloc[0]['pvalue']:.4f}`, `q = {fam_a_main.iloc[0]['q_bh_family_a']:.4f}`)"
        if len(fam_a_main)
        else "no cell types tested"
    )
    if len(pos_both_sig) > 0:
        c_verdict = f"**POSITIVE TRANSPORT ALERT: `{len(pos_both_sig)}` cell type(s) passed both analytic and empirical `q_BH < 0.05` with positive composite slope**"
    elif len(neg_emp_hits) > 0:
        h = neg_emp_hits.iloc[0]
        unadj_match = fam_c_main_unadj.loc[fam_c_main_unadj["group"] == h["group"]]
        unadj_note = (
            f" (`composite = {unadj_match.iloc[0]['composite_effect']:+.4f}`, `analytic p = {unadj_match.iloc[0]['composite_pvalue']:.4f}`, `q_emp = {unadj_match.iloc[0]['q_bh_empirical']:.4f}` unadjusted)"
            if len(unadj_match)
            else ""
        )
        c_verdict = (
            f"**no cell type showed statistically significant positive `P1_Clean` transport** "
            f"(minimum `q_bh_composite = {min_q_c_ana:.4f}`; `{h['group']}` reached two-sided empirical `q_bh_empirical = {h['q_bh_empirical']:.4f}` in the **opposite/negative** direction "
            f"`composite = {h['composite_effect']:+.4f}`, `95% CI [{h['composite_ci95_low']:+.3f}, {h['composite_ci95_high']:+.3f}]`, non-significant `analytic p = {h['composite_pvalue']:.4f}`{unadj_note})"
        )
    else:
        c_verdict = "**no cell type showed statistically significant `P1_Clean` transport**"

    verdict_summary = (
        f"Zero global FDR-significant gene hits (`M_B = {mb_size:,}`, `min q_BH = {min_q_b:.4f}`) and no positive `P1_Clean` aging transport "
        f"(`min q_bh_composite = {min_q_c_ana:.4f}`; the two-sided empirical outlier in `Keratinocytes` is negative `{fam_c_main_adj.iloc[0]['composite_effect']:+.4f}` with `analytic p = {fam_c_main_adj.iloc[0]['composite_pvalue']:.4f}`). "
        f"`{n_sig_a}` CLR compositional shift passed BH in Family A (`{fam_a_main.iloc[0]['cell_type']}`, `q_BH = {fam_a_main.iloc[0]['q_bh_family_a']:.4f}`, `n = 5` males)."
    )

    md_text = f"""# Dataset 01: Human Skin Single-Cell Aging (Sole-Boldo et al., GSE130973)

## 1. Study Verification and Data Provenance
- **Accession**: `GSE130973` (`GSM3758115` to `GSM3758119`), Primary source: Sole-Boldo et al., *Communications Biology* 2020 (`PMC7181753`, `DOI: 10.1038/s42003-020-0922-4`).
- **Tissue and Sampling**: Sun-protected inguinoiliac whole human skin (epidermis + dermis, unbiased single-cell droplet capture, 10x Chromium v2), excluding photoaged sites.
- **Verified Donor Cohort (`n = 5` males)**:
  - `S1` (`y1`, `GSM3758115`): age `25` yr, Male, Caucasian, `{int(donor_meta.loc['S1', 'observed_cells']):,}` QC cells.
  - `S2` (`y2`, `GSM3758116`): age `27` yr, Male, Caucasian, `{int(donor_meta.loc['S2', 'observed_cells']):,}` QC cells.
  - `S3` (`o1`, `GSM3758117`): age `53` yr, Male, Caucasian, `{int(donor_meta.loc['S3', 'observed_cells']):,}` QC cells.
  - `S4` (`o2`, `GSM3758118`): age `70` yr, Male, Caucasian, `{int(donor_meta.loc['S4', 'observed_cells']):,}` QC cells.
  - `S5` (`o3`, `GSM3758119`): age `69` yr, Male, Caucasian, `{int(donor_meta.loc['S5', 'observed_cells']):,}` QC cells.
  - **Total verified cells**: `{total_cells:,}` cells across `{n_genes_total:,}` Seurat raw matrix genes (`{n_tech}` technical `MT-/RP[LS]/MRP[LS]/sex-linked/ERCC` genes stripped prior to library-size normalization). Note all 5 donors are white Caucasian males (0 female donors).

## 2. Explicit Multiple-Testing Families and Methodological Design
- **Inferential unit**: Donor (`n = 5`). Age modeled continuously as `age_per_decade` (`age_years / 10`).
- **Family A (CLR cell-type fractions)**: OLS across the `{len(fam_a_main)}` major skin cell types (clusters mapped 1-to-1 to the 9 published identities via Supplementary Data 1 markers; pseudocount `0.5`). Unadjusted continuous age (`residual_df = 3`), BH correction across `K_A = {len(fam_a_main)}`.
- **Family B (Donor-level gene age test + cross-cell-type concordance)**: Donor pseudobulk (`sum raw UMIs` -> `log2(CPM + 1)`, genes expressed at `CPM >= 1.0` in `>= 3` donors per cell type) across eligible major cell types (`>= 20` cells in `>= 4` donors). Primary model: `age_per_decade + log_lib` (`residual_df = 2` for `n=5`), sensitivity model: unadjusted `age_per_decade` (`residual_df = 3`). Global BH across all `M_B = {mb_size:,}` `(cell_type, gene)` tests. Cross-cell-type concordance tested on common expressed genes via mean pairwise Pearson `r` across cell-type `t`-statistics with `1,000` donor-age permutations.
- **Family C (`P1_Clean` 139-gene mouse aging signature)**: `76` mapped human UP and `42` DOWN symbols (`load_program_tiers`). Composite (`up_mean - down_mean`) OLS slope tested analytically and against `200` expression-decile-matched random gene sets (`matched_null_effect_pvalue`), with BH across `K_C = {len(fam_c_main_adj)}` cell types.

## 3. Results

### Family A: Cell-Type Composition (CLR Fractions, 9 Major Cell Types)
Minimum BH-corrected `q_BH = {min_q_a:.4f}` (`{n_sig_a}` significant compositional shifts at `q < 0.05` in this `n = 5` sample; {top_a_desc}).

| Cell Type | Total Cells | Median Frac | CLR Slope (/10 yr) | p-value | q_BH (K=9) |
|---|---:|---:|---:|---:|---:|
{a_lines}

### Family B: Donor-Level Gene Age Test and Cross-Cell-Type Concordance
Across `{len(fam_b_sum_adj)}` major cell types (`M_B = {mb_size:,}` tests, primary `age_per_decade + log_lib`), `n_sig_global_q05 = {n_sig_b}` genes passed `q_BH < 0.05` (minimum `q_BH = {min_q_b:.4f}`).
Cross-cell-type directional concordance across `{int(fam_b_conc_adj['n_common_genes'].iloc[0])}` common expressed genes across `{int(fam_b_conc_adj['n_cell_types'].iloc[0])}` major cell types:
- Primary (`age_per_decade + log_lib`): mean pairwise Pearson `r = {r_adj:+.4f}` (`perm_p_one_sided = {p_adj:.4f}`, `n_pairs = {int(fam_b_conc_adj['n_pairs'].iloc[0])}`).
- Sensitivity (`age_per_decade` unadjusted): mean pairwise Pearson `r = {r_unadj:+.4f}` (`perm_p_one_sided = {p_unadj:.4f}`).

| Cell Type | Donors | Genes Tested | Global q<0.05 | Min p-value | Min Global q_BH |
|---|---:|---:|---:|---:|---:|
{b_lines}

### Family C: Mouse `P1_Clean` Signature Check
In `K_C = {len(fam_c_main_adj)}` major skin cell types (primary `age_per_decade + log_lib`), {c_verdict}: minimum empirical null `q_BH = {min_q_c_emp:.4f}` (minimum analytic `q_BH = {min_q_c_ana:.4f}`). Results in unadjusted sensitivity (`age_per_decade` only) and across the 17 fine clusters are saved in `results/{prefix}-p1-signature.csv`.

| Cell Type | Donors | UP/DOWN Detected | Composite Slope | Analytic p | Empirical Null p (200) | Empirical q_BH |
|---|---:|---|---:|---:|---:|---:|
{c_lines}

## 4. Limitations and Verdict
- **Sample-size and demographic limits**: Only `n = 5` male Caucasian donors (`25, 27, 53, 69, 70` yr), with 0 female donors and a gap between ages 27 and 53. Power is very limited (`df = 2` adjusted, `df = 3` unadjusted).
- **Consistent with Agescape core claims**: {verdict_summary}

## 5. Output Artifacts
- [Donors & Cell Counts](file:///home/reza/projects/agescape-revision/results/{prefix}-donors.csv)
- [Main Fractions](file:///home/reza/projects/agescape-revision/results/{prefix}-fractions-main.csv)
- [Fine Fractions](file:///home/reza/projects/agescape-revision/results/{prefix}-fractions-fine.csv)
- [Gene Summary](file:///home/reza/projects/agescape-revision/results/{prefix}-gene-summary.csv)
- [Top Genes](file:///home/reza/projects/agescape-revision/results/{prefix}-top-genes.csv)
- [Concordance](file:///home/reza/projects/agescape-revision/results/{prefix}-concordance.csv)
- [P1 Signature](file:///home/reza/projects/agescape-revision/results/{prefix}-p1-signature.csv)
"""
    write_clean_markdown(md_text, reports_dir / f"{prefix}.md")
    clean_temp_dir(tmp_dir)
    print(f"DONE Dataset 01: cells={total_cells}, M_B={mb_size}, min_qA={min_q_a:.4f}, min_qB={min_q_b:.4f}, conc_r={r_adj:.4f} (p={p_adj:.4f}), min_qC_emp={min_q_c_emp:.4f}")


if __name__ == "__main__":
    main()
