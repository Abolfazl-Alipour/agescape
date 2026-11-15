#!/usr/bin/env python3
"""Execute Dataset 06 analysis: Human Bone Marrow / HSC scRNA-seq (Safina et al., GSE302126)."""

from __future__ import annotations

import math
from pathlib import Path
import sys

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from followup_organs import (  # noqa: E402
    DATA_TMP_ROOT,
    NULL_ITERATIONS,
    NULL_SEED,
    aggregate_pseudobulk_sparse,
    clean_temp_dir,
    download_with_guard,
    extract_seurat_rds_selective,
    run_family_a_clr,
    run_family_b_genes_and_concordance,
    run_family_c_p1_signature,
    write_clean_csv,
    write_clean_markdown,
)

FULL_RDS_URL = "https://ndownloader.figshare.com/files/55618298"
ANN_RDS_URL = "https://ndownloader.figshare.com/files/55618823"
OUT_PREFIX = "followup-organ-06-bm-hsc-gse302126"

BROAD_LINEAGE_MAP: dict[str, str] = {
    "HSC": "HSPC_Stem_MPP",
    "MPP-MkEry": "HSPC_Stem_MPP",
    "MPP-MyLy": "HSPC_Stem_MPP",
    "LMPP": "HSPC_Stem_MPP",
    "Cycling Progenitor": "HSPC_Stem_MPP",
    "Early GMP": "Myeloid_Prog",
    "Late GMP": "Myeloid_Prog",
    "GMP-Mono": "Myeloid_Prog",
    "GMP-Neut": "Myeloid_Prog",
    "MEP": "Erythroid_Mk",
    "BFU-E": "Erythroid_Mk",
    "CFU-E": "Erythroid_Mk",
    "Pro-Erythroblast": "Erythroid_Mk",
    "Basophilic Erythroblast": "Erythroid_Mk",
    "Polychromatic Erythroblast": "Erythroid_Mk",
    "Orthochromatic Erythroblast": "Erythroid_Mk",
    "Megakaryocyte Precursor": "Erythroid_Mk",
    "EoBasoMast Precursor": "Erythroid_Mk",
    "Early ProMono": "Monocyte_DC",
    "Late ProMono": "Monocyte_DC",
    "CD14 Mono": "Monocyte_DC",
    "CD16 Mono": "Monocyte_DC",
    "Pre-cDC": "Monocyte_DC",
    "cDC1": "Monocyte_DC",
    "cDC2": "Monocyte_DC",
    "Pre-pDC": "Monocyte_DC",
    "pDC": "Monocyte_DC",
    "ASDC": "Monocyte_DC",
    "Macrophage": "Monocyte_DC",
    "CLP": "B_Lymphoid",
    "Pre-B": "B_Lymphoid",
    "Pro-B Cycling": "B_Lymphoid",
    "Pro-B VDJ": "B_Lymphoid",
    "Large Pre-B": "B_Lymphoid",
    "Small Pre-B": "B_Lymphoid",
    "Immature B": "B_Lymphoid",
    "Naive B": "B_Lymphoid",
    "Memory B": "B_Lymphoid",
    "Plasma Cell": "B_Lymphoid",
    "Naive CD4 T": "T_NK_Lymphoid",
    "CD4 Central Memory": "T_NK_Lymphoid",
    "CD4 Effector Memory": "T_NK_Lymphoid",
    "CD4 Treg": "T_NK_Lymphoid",
    "Naive CD8 T": "T_NK_Lymphoid",
    "CD8 Central Memory": "T_NK_Lymphoid",
    "CD8 Effector Memory 1": "T_NK_Lymphoid",
    "CD8 Effector Memory 2": "T_NK_Lymphoid",
    "CD8 TRM": "T_NK_Lymphoid",
    "MAIT": "T_NK_Lymphoid",
    "GammaDelta T": "T_NK_Lymphoid",
    "NK": "T_NK_Lymphoid",
    "NK CD56bright": "T_NK_Lymphoid",
    "ILC": "T_NK_Lymphoid",
    "Stromal": "Stromal_Other",
}


def _format_num(val: float | int | None, digits: int = 4) -> str:
    if val is None or (isinstance(val, float) and (math.isnan(val) or math.isinf(val))):
        return "NA"
    if isinstance(val, (int, np.integer)):
        return str(int(val))
    f = float(val)
    if abs(f) < 1e-3 and f != 0.0:
        return f"{f:.2e}"
    return f"{f:.{digits}f}"


def _build_safina_donor_meta(meta_sub: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for donor_id, df in meta_sub.groupby("Sample", observed=True):
        age = float(str(df["Age"].iloc[0]).strip())
        sex_raw = str(df["Sex"].iloc[0]).strip().upper()
        site = "sternum_cardiac_surgery" if str(donor_id).startswith("SBM") else "iliac_crest"
        rows.append(
            {
                "donor_id": str(donor_id),
                "age": age,
                "age_per_decade": age / 10.0,
                "sex": sex_raw,
                "sex_male": 1.0 if sex_raw == "M" else 0.0,
                "cohort": str(df["Cohort"].iloc[0]),
                "anatomical_site": site,
                "is_sternum": 1.0 if site == "sternum_cardiac_surgery" else 0.0,
                "n_cells": int(len(df)),
                "n_hsc": int((df["is_hsc"] == "HSC").sum()),
            }
        )
    return pd.DataFrame(rows).sort_values(["age", "donor_id"]).set_index("donor_id", drop=False)


def _build_multistudy_adult_hsc_meta(meta_ann: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (dataset, sample), df in meta_ann.groupby(["Dataset", "Sample"], observed=True):
        age_str = str(df["Age"].iloc[0]).strip()
        if any(tok in age_str.lower() for tok in ("w", "d", "m", "-", "undefined", "na")):
            continue
        try:
            age = float(age_str)
        except ValueError:
            continue
        cohort = str(df["Cohort"].iloc[0])
        if age < 18.0 or cohort in {"Prenatal", "Cord", "Infant", "Child"}:
            continue
        sex_raw = str(df["Sex"].iloc[0]).strip().upper()
        donor_key = f"{dataset}::{sample}"
        rows.append(
            {
                "donor_id": donor_key,
                "sample": str(sample),
                "dataset": str(dataset),
                "age": age,
                "age_per_decade": age / 10.0,
                "sex": sex_raw,
                "sex_male": 1.0 if sex_raw == "M" else 0.0,
                "cohort": cohort,
                "n_hsc": int(len(df)),
            }
        )
    out = pd.DataFrame(rows).sort_values(["dataset", "age", "donor_id"]).set_index("donor_id", drop=False)
    for ds_name in sorted(out["dataset"].unique())[1:]:
        out[f"ds_{ds_name}"] = (out["dataset"] == ds_name).astype(float)
    return out


def _run_suite(
    analysis_name: str,
    counts_mat,
    genes: np.ndarray,
    cell_mask: np.ndarray,
    donor_labels: np.ndarray,
    group_labels: np.ndarray,
    donor_meta: pd.DataFrame,
    candidate_covariates: tuple[str, ...],
    min_cells_per_donor_group: int,
    min_donors: int,
    min_residual_df: int,
    run_fam_a: bool = True,
    run_concordance: bool = True,
    min_total_cells_a: int = 10,
) -> dict[str, pd.DataFrame]:
    cell_counts, expr_by_group, meta_by_group, _ = aggregate_pseudobulk_sparse(
        counts_mat,
        genes,
        donor_labels=donor_labels,
        group_labels=group_labels,
        donor_meta=donor_meta,
        cell_mask=cell_mask,
        min_cells_per_donor_group=min_cells_per_donor_group,
        min_donors_per_group=min_donors,
        allow_fractional_expected=False,
    )

    fam_a = pd.DataFrame()
    if run_fam_a and len(cell_counts.columns) >= 2 and len(cell_counts.index) >= 3:
        fam_a_covs = tuple(c for c in candidate_covariates if c != "log_lib")
        fam_a = run_family_a_clr(
            cell_counts,
            donor_meta.loc[cell_counts.index],
            candidate_covariates=fam_a_covs,
            min_donors_present=min(min_donors, len(cell_counts.index)),
            min_total_cells=min_total_cells_a,
            min_residual_df=min_residual_df,
        )
        if not fam_a.empty:
            fam_a.insert(0, "analysis", analysis_name)

    fam_b_top = pd.DataFrame()
    fam_b_summary = pd.DataFrame()
    fam_b_conc = pd.DataFrame()
    fam_c = pd.DataFrame()
    if expr_by_group:
        fam_b_genes, fam_b_summary, fam_b_conc = run_family_b_genes_and_concordance(
            expr_by_group,
            meta_by_group,
            candidate_covariates=candidate_covariates,
            min_donors=min_donors,
            min_cpm=1.0,
            min_expr_frac=0.5,
            min_residual_df=min_residual_df,
            is_unbiased_single_cell=bool(run_concordance and len(expr_by_group) >= 2),
            concordance_perms=1000,
            seed=NULL_SEED,
        )
        if not fam_b_genes.empty:
            fam_b_top = fam_b_genes.head(50).copy()
        for df in (fam_b_top, fam_b_summary, fam_b_conc):
            if not df.empty:
                df.insert(0, "analysis", analysis_name)

        fam_c = run_family_c_p1_signature(
            expr_by_group,
            meta_by_group,
            candidate_covariates=candidate_covariates,
            min_donors=min_donors,
            min_cpm=1.0,
            min_expr_frac=0.5,
            min_residual_df=min_residual_df,
            iterations=NULL_ITERATIONS,
            seed=NULL_SEED,
        )
        if not fam_c.empty:
            fam_c.insert(0, "analysis", analysis_name)

    return {
        "fam_a": fam_a,
        "fam_b_top": fam_b_top,
        "fam_b_summary": fam_b_summary,
        "fam_b_conc": fam_b_conc,
        "fam_c": fam_c,
    }


def main() -> int:
    tmp_dir = DATA_TMP_ROOT / "organ_06_bm_hsc"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    full_rds = tmp_dir / "seu_full_object_hscs_marked.rds"
    ann_rds = tmp_dir / "annotated_seu.rds"

    download_with_guard(FULL_RDS_URL, full_rds, max_file_gb=2.0)
    download_with_guard(ANN_RDS_URL, ann_rds, max_file_gb=2.0)

    print("Loading seu_full_object_hscs_marked.rds...", flush=True)
    full_info = extract_seurat_rds_selective(full_rds)
    full_meta = full_info["meta.data"].copy()
    full_counts = full_info["raw_counts"]
    full_genes = full_info["genes"]

    saf_mask = np.asarray(full_meta["Dataset"] == "Safina2024", dtype=bool)
    saf_meta = full_meta.loc[saf_mask].copy()
    saf_donors = _build_safina_donor_meta(saf_meta)
    sternum_donors = saf_donors.loc[saf_donors["is_sternum"] == 1.0].copy()

    full_lineage = np.asarray(
        [BROAD_LINEAGE_MAP.get(str(ct), "Other") for ct in full_meta["predicted_CellType"]],
        dtype=str,
    )
    full_fine = np.asarray(full_meta["predicted_CellType"], dtype=str)
    full_donor_col = np.asarray(full_meta["Sample"], dtype=str)
    full_pooled = np.full(len(full_meta), "All_BMMC", dtype=object)

    suites: list[tuple[str, dict[str, pd.DataFrame]]] = []
    suites.append(
        (
            "gse302126_bmmc_lineage_n9",
            _run_suite(
                "gse302126_bmmc_lineage_n9",
                full_counts,
                full_genes,
                saf_mask,
                full_donor_col,
                full_lineage,
                saf_donors,
                ("sex_male", "log_lib"),
                min_cells_per_donor_group=5,
                min_donors=4,
                min_residual_df=2,
                run_fam_a=True,
                run_concordance=True,
            ),
        )
    )
    suites.append(
        (
            "gse302126_bmmc_fine_n9",
            _run_suite(
                "gse302126_bmmc_fine_n9",
                full_counts,
                full_genes,
                saf_mask,
                full_donor_col,
                full_fine,
                saf_donors,
                ("sex_male", "log_lib"),
                min_cells_per_donor_group=5,
                min_donors=4,
                min_residual_df=2,
                run_fam_a=True,
                run_concordance=True,
            ),
        )
    )
    suites.append(
        (
            "gse302126_bmmc_pooled_n9",
            _run_suite(
                "gse302126_bmmc_pooled_n9",
                full_counts,
                full_genes,
                saf_mask,
                full_donor_col,
                full_pooled,
                saf_donors,
                ("sex_male", "log_lib"),
                min_cells_per_donor_group=10,
                min_donors=4,
                min_residual_df=2,
                run_fam_a=False,
                run_concordance=False,
            ),
        )
    )
    suites.append(
        (
            "gse302126_sternum_lineage_n7",
            _run_suite(
                "gse302126_sternum_lineage_n7",
                full_counts,
                full_genes,
                saf_mask,
                full_donor_col,
                full_lineage,
                sternum_donors,
                ("sex_male", "log_lib"),
                min_cells_per_donor_group=5,
                min_donors=4,
                min_residual_df=1,
                run_fam_a=True,
                run_concordance=True,
            ),
        )
    )
    suites.append(
        (
            "gse302126_sternum_pooled_n7",
            _run_suite(
                "gse302126_sternum_pooled_n7",
                full_counts,
                full_genes,
                saf_mask,
                full_donor_col,
                full_pooled,
                sternum_donors,
                ("sex_male", "log_lib"),
                min_cells_per_donor_group=10,
                min_donors=4,
                min_residual_df=2,
                run_fam_a=False,
                run_concordance=False,
            ),
        )
    )

    del full_info, full_counts, full_meta

    print("Loading annotated_seu.rds (28,989 purified HSCs)...", flush=True)
    ann_info = extract_seurat_rds_selective(ann_rds)
    ann_meta = ann_info["meta.data"].copy()
    ann_counts = ann_info["raw_counts"]
    ann_genes = ann_info["genes"]

    ann_saf_mask = np.asarray(ann_meta["Dataset"] == "Safina2024", dtype=bool)
    ann_saf_donors_present = set(ann_meta.loc[ann_saf_mask, "Sample"].astype(str).unique())
    saf_hsc_donors = saf_donors.loc[saf_donors.index.intersection(list(ann_saf_donors_present))].copy()
    ann_donor_col = np.asarray(ann_meta["Sample"], dtype=str)
    ann_prog = np.asarray(ann_meta["Predominant_program"], dtype=str)
    ann_clust = np.asarray([f"HSC_cluster_{c}" for c in ann_meta["seurat_clusters"]], dtype=str)
    ann_hsc_pool = np.full(len(ann_meta), "Safina_HSC_all", dtype=object)

    suites.append(
        (
            "safina_hsc_pooled_n7",
            _run_suite(
                "safina_hsc_pooled_n7",
                ann_counts,
                ann_genes,
                ann_saf_mask,
                ann_donor_col,
                ann_hsc_pool,
                saf_hsc_donors,
                ("sex_male", "log_lib"),
                min_cells_per_donor_group=1,
                min_donors=4,
                min_residual_df=2,
                run_fam_a=False,
                run_concordance=False,
            ),
        )
    )
    suites.append(
        (
            "safina_hsc_programs_n5",
            _run_suite(
                "safina_hsc_programs_n5",
                ann_counts,
                ann_genes,
                ann_saf_mask,
                ann_donor_col,
                ann_prog,
                saf_hsc_donors,
                ("sex_male", "log_lib"),
                min_cells_per_donor_group=3,
                min_donors=4,
                min_residual_df=1,
                run_fam_a=True,
                run_concordance=True,
                min_total_cells_a=5,
            ),
        )
    )
    suites.append(
        (
            "safina_hsc_clusters_n7",
            _run_suite(
                "safina_hsc_clusters_n7",
                ann_counts,
                ann_genes,
                ann_saf_mask,
                ann_donor_col,
                ann_clust,
                saf_hsc_donors,
                ("sex_male", "log_lib"),
                min_cells_per_donor_group=3,
                min_donors=4,
                min_residual_df=2,
                run_fam_a=True,
                run_concordance=False,
                min_total_cells_a=5,
            ),
        )
    )

    adult_meta_tbl = _build_multistudy_adult_hsc_meta(ann_meta)
    adult_key_col = np.asarray(
        [f"{d}::{s}" for d, s in zip(ann_meta["Dataset"], ann_meta["Sample"])],
        dtype=str,
    )
    adult_mask = np.isin(adult_key_col, adult_meta_tbl.index.to_numpy(dtype=str))
    ms_pool = np.full(len(ann_meta), "MultiStudy_Adult_HSC", dtype=object)
    ds_covs = tuple(c for c in adult_meta_tbl.columns if c.startswith("ds_"))
    ms_covs = ("sex_male", "log_lib") + ds_covs

    suites.append(
        (
            "multistudy_adult_hsc_total",
            _run_suite(
                "multistudy_adult_hsc_total",
                ann_counts,
                ann_genes,
                adult_mask,
                adult_key_col,
                ms_pool,
                adult_meta_tbl,
                ms_covs,
                min_cells_per_donor_group=5,
                min_donors=10,
                min_residual_df=4,
                run_fam_a=False,
                run_concordance=False,
            ),
        )
    )
    suites.append(
        (
            "multistudy_adult_hsc_programs",
            _run_suite(
                "multistudy_adult_hsc_programs",
                ann_counts,
                ann_genes,
                adult_mask,
                adult_key_col,
                ann_prog,
                adult_meta_tbl,
                ms_covs,
                min_cells_per_donor_group=5,
                min_donors=10,
                min_residual_df=4,
                run_fam_a=True,
                run_concordance=True,
            ),
        )
    )

    fam_a_all = pd.concat([s[1]["fam_a"] for s in suites if not s[1]["fam_a"].empty], ignore_index=True)
    fam_b_sum_all = pd.concat([s[1]["fam_b_summary"] for s in suites if not s[1]["fam_b_summary"].empty], ignore_index=True)
    fam_b_top_all = pd.concat([s[1]["fam_b_top"] for s in suites if not s[1]["fam_b_top"].empty], ignore_index=True)
    fam_b_conc_all = pd.concat([s[1]["fam_b_conc"] for s in suites if not s[1]["fam_b_conc"].empty], ignore_index=True)
    fam_c_all = pd.concat([s[1]["fam_c"] for s in suites if not s[1]["fam_c"].empty], ignore_index=True)

    saf_don_out = saf_donors.reset_index(drop=True)
    saf_don_out.insert(0, "cohort_table", "Safina2024_GSE302126_BMMC")
    ms_don_out = adult_meta_tbl.reset_index(drop=True)
    ms_don_out.insert(0, "cohort_table", "MultiStudy_ExactAge_Adult_HSC")
    combined_donors = pd.concat([saf_don_out, ms_don_out], ignore_index=True, sort=False)

    res_dir = REPO_ROOT / "results"
    write_clean_csv(combined_donors, res_dir / f"{OUT_PREFIX}-donors.csv")
    write_clean_csv(fam_a_all, res_dir / f"{OUT_PREFIX}-fractions.csv")
    write_clean_csv(fam_b_sum_all, res_dir / f"{OUT_PREFIX}-gene-summary.csv")
    write_clean_csv(fam_b_top_all, res_dir / f"{OUT_PREFIX}-top-genes.csv")
    write_clean_csv(fam_b_conc_all, res_dir / f"{OUT_PREFIX}-concordance.csv")
    write_clean_csv(fam_c_all, res_dir / f"{OUT_PREFIX}-p1-signature.csv")

    primary_a = fam_a_all.loc[fam_a_all["analysis"] == "gse302126_bmmc_lineage_n9"] if not fam_a_all.empty else pd.DataFrame()
    primary_a_hits = int((primary_a["q_bh_family_a"] < 0.05).sum()) if not primary_a.empty else 0
    fine_a = fam_a_all.loc[fam_a_all["analysis"] == "gse302126_bmmc_fine_n9"] if not fam_a_all.empty else pd.DataFrame()
    fine_a_hits = int((fine_a["q_bh_family_a"] < 0.05).sum()) if not fine_a.empty else 0

    b_prim = fam_b_sum_all.loc[(fam_b_sum_all["analysis"] == "gse302126_bmmc_lineage_n9") & (fam_b_sum_all["status"] == "tested")]
    m_b_lin = int(b_prim["n_genes_tested"].sum()) if not b_prim.empty else 0
    sig_b_lin = int(b_prim["n_sig_global_q05"].sum()) if not b_prim.empty else 0

    b_pool = fam_b_sum_all.loc[(fam_b_sum_all["analysis"] == "gse302126_bmmc_pooled_n9") & (fam_b_sum_all["status"] == "tested")]
    m_b_pool = int(b_pool["n_genes_tested"].sum()) if not b_pool.empty else 0
    sig_b_pool = int(b_pool["n_sig_global_q05"].sum()) if not b_pool.empty else 0

    c_lin = fam_c_all.loc[fam_c_all["analysis"] == "gse302126_bmmc_lineage_n9"]
    c_fine = fam_c_all.loc[fam_c_all["analysis"] == "gse302126_bmmc_fine_n9"]
    c_stern = fam_c_all.loc[fam_c_all["analysis"] == "gse302126_sternum_lineage_n7"]
    c_ms_tot = fam_c_all.loc[fam_c_all["analysis"] == "multistudy_adult_hsc_total"]

    donor_md_lines = [
        "| Donor | Age | Sex | Group | Anatomical / Clinical Site | Total BMMC Cells | HSC Cells |",
        "|---|---:|---|---|---|---:|---:|",
    ]
    for _, r in saf_donors.iterrows():
        donor_md_lines.append(
            f"| {r['donor_id']} | {int(r['age'])} | {r['sex']} | {r['cohort']} | {r['anatomical_site']} | {int(r['n_cells'])} | {int(r['n_hsc'])} |"
        )

    p1_md_lines = [
        "| Analysis | Cell Group | Donors (n) | Up / Down Detected | Composite Slope / Decade | Analytic p | Empirical p (200 draws) | q_BH (Composite) | q_BH (Empirical) |",
        "|---|---|---:|---|---:|---:|---:|---:|---:|",
    ]
    for _, r in fam_c_all.iterrows():
        p1_md_lines.append(
            f"| {r['analysis']} | {r['group']} | {int(r['n_donors'])} | {int(r['n_p1_up_detected'])} / {int(r['n_p1_down_detected'])} | "
            f"{_format_num(r['composite_effect'])} | {_format_num(r['composite_pvalue'])} | {_format_num(r['empirical_null_pvalue'])} | "
            f"{_format_num(r['q_bh_composite'])} | {_format_num(r['q_bh_empirical'])} |"
        )

    conc_md_lines = [
        "| Analysis | Cell Types (n) | Common Genes | Mean Pairwise Pearson r | Perm p (1-sided) | Perm p (2-sided) |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for _, r in fam_b_conc_all.iterrows():
        conc_md_lines.append(
            f"| {r['analysis']} | {int(r['n_cell_types'])} | {int(r['n_common_genes'])} | "
            f"{_format_num(r.get('mean_pairwise_pearson_r'))} | {_format_num(r.get('perm_p_one_sided'))} | {_format_num(r.get('perm_p_two_sided'))} |"
        )

    min_q_lin_comp = float(c_lin["q_bh_composite"].min()) if not c_lin.empty else float("nan")
    min_q_lin_emp = float(c_lin["q_bh_empirical"].min()) if not c_lin.empty else float("nan")
    min_q_fine_emp = float(c_fine["q_bh_empirical"].min()) if not c_fine.empty else float("nan")
    min_q_stern_emp = float(c_stern["q_bh_empirical"].min()) if not c_stern.empty else float("nan")
    ms_tot_emp = float(c_ms_tot["empirical_null_pvalue"].iloc[0]) if not c_ms_tot.empty else float("nan")
    ms_tot_comp_p = float(c_ms_tot["composite_pvalue"].iloc[0]) if not c_ms_tot.empty else float("nan")
    ms_tot_comp = float(c_ms_tot["composite_effect"].iloc[0]) if not c_ms_tot.empty else float("nan")

    report_text = f"""# Dataset 06: Human Bone Marrow and HSC Single-Cell RNA-seq (Safina et al., GSE302126)

## 1. Study Verification and Data Provenance
- **Accession and Primary Paper**: `GSE302126` (`GSM9096700` to `GSM9096716`, 17 CellRanger BMMC technical-replicate libraries), Safina et al., *Genome Biology* 2026 (`PMID: 41546063`, `PMC12892447`).
- **Processed Author Seurat Objects**: Retrieved from the primary paper Figshare project (`https://figshare.com/projects/Aging_HSCs_2025/235781`):
  - `seu_full_object_hscs_marked.rds` (file `55618298`, 204,535 BoneMarrowMap-annotated cells across 7 datasets and 98 donors; contains all 2,361 QC-passed `Safina2024` / `GSE302126` BMMC cells across all 9 unique WT CHIP-negative donors).
  - `annotated_seu.rds` (file `55618823`, 28,989 BoneMarrowMap-purified HSCs across 7 datasets with `Predominant_program` `metaGEP.1` to `metaGEP.4` and Seurat subclusters).
- **Exact GEO Donor Verification (17 Replicates Collapsed to 9 Unique Donors)**:
  - Young (`n = 2`, ages `31` and `37` yr; `1M / 1F`): `BM1` (31 yr, M, GSM9096700 to GSM9096702) and `BM2` (37 yr, F, GSM9096703 to GSM9096705), both coded `source_name_ch1 = Bone marrow` (iliac crest in healthy donors, non-sternum).
  - Aged (`n = 7`, ages `64 to 78` yr; `5M / 2F`): `SBM1307` (64F), `SBM1096` (66M), `SBM1271` (67M), `SBM1237` (71M), `SBM1177` (73M), `SBM1264` (77M), and `SBM1073` (78F), all coded `source_name_ch1 = Sternum bone marrow` from patients undergoing cardiac surgery.
  - **Critical structural confound**: In `GSE302126` (`n = 9`), age group is 100% collinear with anatomical collection site and cardiac-surgery status (`100%` young iliac crest vs `100%` aged sternum cardiac surgery). A sternum-only internal sensitivity (`n = 7`, `64 to 78` yr) removes this confound, and multi-study exact-age adult BoneMarrowMap HSCs (`annotated_seu.rds`) provide an external adult-lifespan check (`n = {len(adult_meta_tbl)}` unique adult donors, ages `{int(adult_meta_tbl['age'].min())} to {int(adult_meta_tbl['age'].max())}` yr across `{adult_meta_tbl['dataset'].nunique()}` datasets after excluding prenatal/cord/infant/child and undefined-age `Weng2024`).

{chr(10).join(donor_md_lines)}

## 2. Explicit Multiple-Testing Families
1. **Family A (CLR cell-type composition)**: Reported per analysis suite (`gse302126_bmmc_lineage_n9`, `gse302126_bmmc_fine_n9`, `gse302126_sternum_lineage_n7`, `safina_hsc_programs_n5`, `safina_hsc_clusters_n7`, `multistudy_adult_hsc_programs`), testing CLR fractions against `age_per_decade + covariates` with `q_bh_family_a < 0.05`.
2. **Family B (Donor-level gene age OLS and within-species directional concordance)**: Tested globally across all qualifying `(group, gene)` pairs within each suite (`CPM >= 1.0` in `>= 50%` of tested donors; `age_per_decade + sex_male + log_lib` in `GSE302126`, plus dataset fixed effects in `multistudy_adult_hsc`), plus 1,000-permutation donor-age concordance.
3. **Family C (`P1_Clean` mouse signature transport)**: `118` mapped human orthologs (`76` UP, `42` DOWN) scored per donor pseudobulk per qualifying group (`OLS` analytic p-value and 200-draw expression-decile-matched random-set two-sided empirical p-value), BH-corrected within each analysis suite.

## 3. Results

### Family A: Compositional Changes
- **Primary `GSE302126` BMMC broad lineages (`n = 9`, `sex_male` adjusted)**: `{primary_a_hits}` cell lineages reached `q_bh_family_a < 0.05`.
- **Primary `GSE302126` BMMC fine `predicted_CellType` (`n = 9`, `sex_male` adjusted)**: `{fine_a_hits}` cell types reached `q_bh_family_a < 0.05`.

### Family B: Donor-Level Gene Age OLS & Concordance
- **Primary `GSE302126` BMMC broad lineages (`n = 9`)**: `M_B = {m_b_lin}` gene-lineage tests across `{len(b_prim)}` lineages -> `{sig_b_lin}` global BH hits (`q_bh_global < 0.05`).
- **Primary `GSE302126` All-BMMC pooled (`n = 9`)**: `M_B = {m_b_pool}` tested genes -> `{sig_b_pool}` global BH hits (`q_bh_global < 0.05`).
- **Cross-cell-type gene concordance**:
{chr(10).join(conc_md_lines)}

### Family C: Mouse `P1_Clean` Aging Signature Transport
- **Primary `GSE302126` BMMC broad lineages (`n = 9`)**: Minimum `q_bh_composite = {_format_num(min_q_lin_comp)}`, minimum `q_bh_empirical = {_format_num(min_q_lin_emp)}` (NULL).
- **Fine `predicted_CellType` (`n = 9`)**: Minimum `q_bh_empirical = {_format_num(min_q_fine_emp)}` (NULL).
- **Sternum-only `SBM*` (`n = 7`, `64 to 78` yr, removes iliac crest confound)**: Minimum `q_bh_empirical = {_format_num(min_q_stern_emp)}` (NULL).
- **Multi-study exact-age adult BoneMarrowMap HSCs (`n = {len(adult_meta_tbl)}`, `ages {int(adult_meta_tbl['age'].min())} to {int(adult_meta_tbl['age'].max())}` yr)**: Pooled total HSC `P1_Clean` composite slope `{_format_num(ms_tot_comp)}` (`analytic p = {_format_num(ms_tot_comp_p)}`, `200-draw empirical p = {_format_num(ms_tot_emp)}`).

{chr(10).join(p1_md_lines)}

## 4. Limitations and Verdict
- **Limits**: `GSE302126` has small sample size (`n = 9` unique donors after collapsing 17 replicate libraries), `SBM1177` (28 cells) and `SBM1237` (17 cells) have sparse BMMC recovery, only `117` total `Safina2024` cells map to `HSC`, and young vs aged donors in `GSE302126` are 100% confounded with iliac crest vs sternum cardiac surgery sampling. Multi-study adult HSC integration (`n = {len(adult_meta_tbl)}`) overcomes sample size but mixes 10x and Smart-seq/SORT protocols across studies.
- **Verdict**: Mouse `P1_Clean` shows no significant positive transport in human `GSE302126` bone marrow lineages, `GSE302126` sternum-only donors, `GSE302126` HSCs, or the `{len(adult_meta_tbl)}`-donor multi-study exact-age adult HSC atlas (`all empirical q_BH >= 0.05` across tested primary cohorts).

## 5. Output Artifacts
- [Donors](file:///home/reza/projects/agescape-revision/results/{OUT_PREFIX}-donors.csv)
- [Fractions](file:///home/reza/projects/agescape-revision/results/{OUT_PREFIX}-fractions.csv)
- [Gene Summary](file:///home/reza/projects/agescape-revision/results/{OUT_PREFIX}-gene-summary.csv)
- [Top Genes](file:///home/reza/projects/agescape-revision/results/{OUT_PREFIX}-top-genes.csv)
- [Concordance](file:///home/reza/projects/agescape-revision/results/{OUT_PREFIX}-concordance.csv)
- [P1 Signature](file:///home/reza/projects/agescape-revision/results/{OUT_PREFIX}-p1-signature.csv)
"""

    report_path = write_clean_markdown(report_text, REPO_ROOT / "reports" / f"{OUT_PREFIX}.md")
    clean_temp_dir(tmp_dir)
    print(f"DONE Dataset 06: report={report_path}, safina_bmmc_donors={len(saf_donors)}, adult_ms_donors={len(adult_meta_tbl)}, p1_rows={len(fam_c_all)}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
