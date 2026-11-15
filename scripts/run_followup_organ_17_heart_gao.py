#!/usr/bin/env python3
"""Organ 17 (Item #14): Gao et al. 2026 integrated human heart snRNA-seq atlas (GSE290367)."""

from __future__ import annotations

import json
from pathlib import Path
import sys
from typing import Any

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem
from scipy import sparse as sp
from statsmodels.stats.multitest import multipletests

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.followup_muscle import NULL_SEED, count_gate
from src.followup_organs import (
    finalize_pseudobulk_group,
    load_p1_clean_human,
    run_family_a_clr,
    run_family_b_genes_and_concordance,
    run_family_c_p1_signature,
    write_clean_csv,
)

H5AD_PATH = Path(
    "/home/reza/projects/agescape-revision-data/tmp_organs/organ17_gao_heart/GSE290367_integrated_RNA_adata.h5ad"
)
RESULTS_DIR = REPO_ROOT / "results"
PREFIX = "followup-organ-17-heart-gao"

P1_ITERATIONS = 2000
CONCORDANCE_PERMS = 1000
MIN_DONOR_CELLS = 20
CHUNK_ROWS = 25000

MIDPOINT_STUDIES = {"Litvinukova 2020", "Kanemaru 2023", "Reichart 2022"}
DUP_PENN_DONORS = {
    "Penn:K1549",
    "Penn:K1561",
    "Penn:K1622",
    "Penn:K1702",
    "Penn:K1718",
    "Penn:K1739",
}
STUDY_DUMMY_MAP = {
    "Chaffin 2022": "study_chaffin_2022",
    "Koenig 2022": "study_koenig_2022",
    "Penn": "study_penn",
    "Simonson 2023": "study_simonson_2023",
    "Kuppe 2022": "study_kuppe_2022",
    "Sim 2021": "study_sim_2021",
}


def build_donor_audit(obs: pd.DataFrame) -> pd.DataFrame:
    """Create a complete 299-donor disposition table across all studies in GSE290367."""
    work = obs[["donor_id", "study", "age", "age_status", "disease", "sex", "technology", "region"]].copy()
    for col in work.columns:
        work[col] = work[col].astype(str)
    work["age_num"] = pd.to_numeric(work["age"], errors="coerce")

    grouped = (
        work.groupby(["donor_id", "study"], sort=True, observed=True)
        .agg(
            age_raw=("age", "first"),
            age_num=("age_num", "first"),
            age_status=("age_status", lambda s: ";".join(sorted(set(s)))),
            disease=("disease", lambda s: ";".join(sorted(set(s)))),
            sex=("sex", lambda s: ";".join(sorted(set(s)))),
            technology=("technology", lambda s: ";".join(sorted(set(s)))),
            regions=("region", lambda s: ";".join(sorted(set(s)))),
            n_nuclei_total=("donor_id", "size"),
        )
        .reset_index()
    )

    reasons: list[str] = []
    included: list[bool] = []
    for row in grouped.itertuples(index=False):
        if row.study in MIDPOINT_STUDIES:
            included.append(False)
            reasons.append("excluded_midpoint_age_brackets")
        elif row.age_status != "postnatal":
            included.append(False)
            reasons.append("excluded_fetal")
        elif row.disease != "ND":
            included.append(False)
            reasons.append(f"excluded_disease_{row.disease}")
        elif not np.isfinite(row.age_num) or row.age_num < 18.0:
            included.append(False)
            reasons.append("excluded_pediatric_under_18")
        elif row.donor_id in DUP_PENN_DONORS:
            included.append(False)
            reasons.append("excluded_duplicate_penn_donor_kept_10x_copy")
        elif row.sex not in {"male", "female"}:
            included.append(False)
            reasons.append("excluded_unknown_sex")
        else:
            included.append(True)
            reasons.append("included_primary_exact_adult_lv")

    grouped["included_primary"] = included
    grouped["disposition_reason"] = reasons
    return grouped.sort_values(["included_primary", "study", "age_num", "donor_id"], ascending=[False, True, True, True]).reset_index(drop=True)


def stream_pseudobulk_counts(
    h5_path: Path,
    keep_mask: np.ndarray,
    donor_ids: np.ndarray,
    cell_types: np.ndarray,
    gene_names: np.ndarray,
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame, dict[str, int]]:
    """Stream layers/counts from HDF5 in sequential chunks and aggregate raw integer donor x gene sums."""
    unique_donors = sorted(set(donor_ids[keep_mask].tolist()))
    unique_cts = sorted(set(cell_types[keep_mask].tolist()))
    donor_to_idx = {d: i for i, d in enumerate(unique_donors)}
    ct_to_idx = {c: i for i, c in enumerate(unique_cts)}
    n_donors = len(unique_donors)
    n_cts = len(unique_cts)
    n_pairs = n_donors * n_cts
    n_vars = len(gene_names)

    pair_sums = np.zeros((n_pairs, n_vars), dtype=np.float64)
    cell_counts = np.zeros((n_donors, n_cts), dtype=np.int64)

    keep_donor_idx = np.full(len(keep_mask), -1, dtype=np.int64)
    keep_ct_idx = np.full(len(keep_mask), -1, dtype=np.int64)
    kept_pos = np.flatnonzero(keep_mask)
    keep_donor_idx[kept_pos] = np.array([donor_to_idx[d] for d in donor_ids[kept_pos]], dtype=np.int64)
    keep_ct_idx[kept_pos] = np.array([ct_to_idx[c] for c in cell_types[kept_pos]], dtype=np.int64)

    gate_totals = {"nonfinite": 0, "negative": 0, "noninteger": 0, "total_nnz": 0, "cells_checked": 0}
    n_obs = len(keep_mask)

    with h5py.File(h5_path, "r") as handle:
        counts_grp = handle["layers"]["counts"]
        indptr = counts_grp["indptr"][:]
        data_ds = counts_grp["data"]
        indices_ds = counts_grp["indices"]

        for r0 in range(0, n_obs, CHUNK_ROWS):
            r1 = min(r0 + CHUNK_ROWS, n_obs)
            local_keep = keep_mask[r0:r1]
            if not np.any(local_keep):
                continue
            lo, hi = int(indptr[r0]), int(indptr[r1])
            local_indptr = indptr[r0 : r1 + 1] - lo
            blk_data = data_ds[lo:hi]
            blk_indices = indices_ds[lo:hi]
            blk_csr = sp.csr_matrix((blk_data, blk_indices, local_indptr), shape=(r1 - r0, n_vars))
            sub_csr = blk_csr[local_keep]

            gate = count_gate(sub_csr.data)
            gate_totals["nonfinite"] += gate["nonfinite"]
            gate_totals["negative"] += gate["negative"]
            gate_totals["noninteger"] += gate["noninteger"]
            gate_totals["total_nnz"] += int(len(sub_csr.data))
            gate_totals["cells_checked"] += int(sub_csr.shape[0])
            if gate["nonfinite"] or gate["negative"] or gate["noninteger"]:
                raise ValueError(f"layers/counts failed integer check in rows {r0}:{r1}: {gate}")

            d_local = keep_donor_idx[r0:r1][local_keep]
            c_local = keep_ct_idx[r0:r1][local_keep]
            np.add.at(cell_counts, (d_local, c_local), 1)
            pair_ids = d_local * n_cts + c_local
            ind_mat = sp.csr_matrix(
                (np.ones(len(pair_ids), dtype=np.float64), (pair_ids, np.arange(len(pair_ids)))),
                shape=(n_pairs, len(pair_ids)),
            )
            blk_sum = ind_mat @ sub_csr
            pair_sums += blk_sum.toarray()

    counts_df = pd.DataFrame(cell_counts, index=unique_donors, columns=unique_cts)
    raw_pb_by_group: dict[str, pd.DataFrame] = {}
    for c_idx, ct_name in enumerate(unique_cts):
        rows = np.arange(n_donors) * n_cts + c_idx
        sub_counts = counts_df[ct_name]
        elig_donors = sub_counts.index[sub_counts >= MIN_DONOR_CELLS].tolist()
        if not elig_donors:
            continue
        elig_rows = [donor_to_idx[d] * n_cts + c_idx for d in elig_donors]
        mat = pair_sums[elig_rows, :].T
        raw_pb_by_group[ct_name] = pd.DataFrame(mat, index=gene_names, columns=elig_donors)

    total_cells = counts_df.sum(axis=1)
    elig_all = total_cells.index[total_cells >= MIN_DONOR_CELLS].tolist()
    all_sum = pair_sums.reshape(n_donors, n_cts, n_vars).sum(axis=1)
    all_rows = [donor_to_idx[d] for d in elig_all]
    raw_pb_by_group["All_LV_Nuclei"] = pd.DataFrame(all_sum[all_rows, :].T, index=gene_names, columns=elig_all)
    counts_df["Total_Cells"] = total_cells
    return raw_pb_by_group, counts_df, gate_totals


def build_donor_meta_table(audit_df: pd.DataFrame) -> pd.DataFrame:
    """Build donor-level covariate frame with age_per_decade, sex_male, and study dummies."""
    inc = audit_df[audit_df["included_primary"]].copy().set_index("donor_id")
    inc["age_years"] = inc["age_num"].astype(float)
    inc["age_per_decade"] = inc["age_years"] / 10.0
    inc["sex_male"] = (inc["sex"] == "male").astype(float)
    for study_name, col_name in STUDY_DUMMY_MAP.items():
        inc[col_name] = (inc["study"] == study_name).astype(float)
    return inc


def run_suite_for_subset(
    label: str,
    donor_subset: list[str],
    raw_pb_by_group: dict[str, pd.DataFrame],
    cell_counts_df: pd.DataFrame,
    donor_meta_all: pd.DataFrame,
    p1_program: dict[str, list[str]],
    adjust_study: bool,
    celltype_names: list[str],
) -> dict[str, Any]:
    """Run Family A, Family B (genes + celltype concordance), and Family C for a defined donor subset."""
    sub_meta = donor_meta_all.loc[donor_subset].copy()
    sub_counts = cell_counts_df.loc[donor_subset, celltype_names].copy()

    cov_a: list[str] = ["sex_male"]
    cov_bc: list[str] = ["sex_male", "log_lib"]
    if adjust_study:
        active_studies = set(sub_meta["study"].tolist())
        for st_name, col_name in STUDY_DUMMY_MAP.items():
            if st_name in active_studies:
                cov_a.append(col_name)
                cov_bc.append(col_name)

    fam_a = run_family_a_clr(
        sub_counts,
        sub_meta,
        age_col="age_per_decade",
        candidate_covariates=tuple(cov_a),
        min_donors_present=3,
        min_total_cells=MIN_DONOR_CELLS,
        min_residual_df=2,
    )
    fam_a.insert(0, "subset", label)

    expr_ct: dict[str, pd.DataFrame] = {}
    meta_ct: dict[str, pd.DataFrame] = {}
    qc_rows: list[dict[str, Any]] = []
    for grp_name, full_raw in raw_pb_by_group.items():
        cols = [d for d in donor_subset if d in full_raw.columns]
        if len(cols) < 4:
            continue
        raw_sub = full_raw[cols]
        c_ser = (
            cell_counts_df.loc[cols, "Total_Cells"]
            if grp_name == "All_LV_Nuclei"
            else cell_counts_df.loc[cols, grp_name]
        )
        log2cpm, grp_meta, stats = finalize_pseudobulk_group(raw_sub, sub_meta.loc[cols], c_ser)
        expr_ct[grp_name] = log2cpm
        meta_ct[grp_name] = grp_meta
        qc_rows.append(
            {
                "subset": label,
                "group": grp_name,
                "n_donors": len(cols),
                "total_cells_in_group": int(c_ser.sum()),
                "median_cells_per_donor": float(c_ser.median()),
                **{k: v for k, v in stats.items() if k != "count_gate"},
            }
        )

    expr_only_ct = {k: v for k, v in expr_ct.items() if k != "All_LV_Nuclei"}
    meta_only_ct = {k: v for k, v in meta_ct.items() if k != "All_LV_Nuclei"}

    genes_ct, summary_ct, conc_df = run_family_b_genes_and_concordance(
        expr_only_ct,
        meta_only_ct,
        age_col="age_per_decade",
        candidate_covariates=tuple(cov_bc),
        min_donors=4,
        min_cpm=1.0,
        min_expr_frac=0.5,
        min_residual_df=2,
        is_unbiased_single_cell=True,
        concordance_perms=CONCORDANCE_PERMS,
        seed=NULL_SEED,
    )

    parts_genes = [genes_ct] if not genes_ct.empty else []
    parts_summary = [summary_ct] if not summary_ct.empty else []
    if "All_LV_Nuclei" in expr_ct:
        genes_all, summary_all, _ = run_family_b_genes_and_concordance(
            {"All_LV_Nuclei": expr_ct["All_LV_Nuclei"]},
            {"All_LV_Nuclei": meta_ct["All_LV_Nuclei"]},
            age_col="age_per_decade",
            candidate_covariates=tuple(cov_bc),
            min_donors=4,
            min_cpm=1.0,
            min_expr_frac=0.5,
            min_residual_df=2,
            is_unbiased_single_cell=False,
            concordance_perms=CONCORDANCE_PERMS,
            seed=NULL_SEED,
        )
        if not genes_all.empty:
            parts_genes.append(genes_all)
        if not summary_all.empty:
            parts_summary.append(summary_all)

    genes_df = pd.concat(parts_genes, ignore_index=True) if parts_genes else pd.DataFrame()
    if not genes_df.empty:
        genes_df["family_b_celltypes13_size"] = int(len(genes_ct))
        genes_df["family_b_all14_size"] = int(len(genes_df))
        genes_df["q_bh_all14"] = np.nan
        valid_all = genes_df["pvalue"].notna()
        if bool(valid_all.any()):
            genes_df.loc[valid_all, "q_bh_all14"] = multipletests(genes_df.loc[valid_all, "pvalue"], method="fdr_bh")[1]
        genes_df["q_bh_celltypes13"] = np.nan
        valid_ct = (genes_df["group"] != "All_LV_Nuclei") & genes_df["pvalue"].notna()
        if bool(valid_ct.any()):
            genes_df.loc[valid_ct, "q_bh_celltypes13"] = multipletests(genes_df.loc[valid_ct, "pvalue"], method="fdr_bh")[1]
        genes_df.insert(0, "subset", label)
    summary_df = pd.concat(parts_summary, ignore_index=True) if parts_summary else pd.DataFrame()
    if not summary_df.empty:
        if not genes_df.empty:
            all14_q = genes_df.groupby("group")["q_bh_all14"].apply(lambda s: int((s < 0.05).sum()))
            summary_df["n_sig_q05_all14"] = summary_df["group"].map(all14_q).fillna(0).astype(int)
        summary_df.insert(0, "subset", label)
    if not conc_df.empty:
        conc_df.insert(0, "subset", label)

    fam_c = run_family_c_p1_signature(
        expr_ct,
        meta_ct,
        p1_program=p1_program,
        age_col="age_per_decade",
        candidate_covariates=tuple(cov_bc),
        min_donors=4,
        min_cpm=1.0,
        min_expr_frac=0.5,
        min_residual_df=2,
        iterations=P1_ITERATIONS,
        seed=NULL_SEED,
    )
    if not fam_c.empty:
        fam_c.insert(0, "subset", label)

    return {
        "fam_a": fam_a,
        "genes": genes_df,
        "summary": summary_df,
        "conc": conc_df,
        "fam_c": fam_c,
        "qc": pd.DataFrame(qc_rows),
    }


def main() -> None:
    if not H5AD_PATH.exists():
        raise FileNotFoundError(f"Missing unpacked H5AD: {H5AD_PATH}")

    print("1. Reading .obs and .var from GSE290367 H5AD...", flush=True)
    with h5py.File(H5AD_PATH, "r") as handle:
        obs = read_elem(handle["obs"])
        var = read_elem(handle["var"])
    gene_names = var.index.astype(str).to_numpy()
    p1_program = load_p1_clean_human(REPO_ROOT)

    print("2. Building full 299-donor disposition audit...", flush=True)
    audit_df = build_donor_audit(obs)
    donor_meta_all = build_donor_meta_table(audit_df)
    inc_set = set(donor_meta_all.index.tolist())
    obs_donors = obs["donor_id"].astype(str).to_numpy()
    obs_cts = obs["final_cell_type"].astype(str).to_numpy()
    keep_mask = np.array([d in inc_set for d in obs_donors], dtype=bool)
    print(f"Primary cohort: {len(inc_set)} unique adult LV donors, {int(keep_mask.sum())} nuclei.", flush=True)

    print("3. Streaming raw integer counts & verifying count_gate across all kept cells...", flush=True)
    raw_pb_by_group, cell_counts_df, gate_totals = stream_pseudobulk_counts(
        H5AD_PATH, keep_mask, obs_donors, obs_cts, gene_names
    )
    print(f"Count gate PASSED: {gate_totals}", flush=True)

    celltype_names = sorted([c for c in cell_counts_df.columns if c != "Total_Cells"])
    merged_counts = donor_meta_all[["study", "age_years", "sex", "technology", "regions"]].join(cell_counts_df, how="right")
    merged_counts.index.name = "donor_id"

    all_donors = donor_meta_all.index.tolist()
    print("4. Running Primary Pooled 7-study model (N=117, study-adjusted)...", flush=True)
    primary = run_suite_for_subset(
        "pooled_7study_primary",
        all_donors,
        raw_pb_by_group,
        cell_counts_df,
        donor_meta_all,
        p1_program,
        adjust_study=True,
        celltype_names=celltype_names,
    )

    sens_specs = [
        ("sensitivity_10x_only_no_penn", donor_meta_all.index[donor_meta_all["study"] != "Penn"].tolist(), True),
        ("sensitivity_studies_ge8", donor_meta_all.index[~donor_meta_all["study"].isin(["Kuppe 2022", "Sim 2021"])].tolist(), True),
        ("study_ENCODE_v4_Snyder", donor_meta_all.index[donor_meta_all["study"] == "ENCODE v4 (Snyder)"].tolist(), False),
        ("study_Koenig_2022", donor_meta_all.index[donor_meta_all["study"] == "Koenig 2022"].tolist(), False),
        ("study_Chaffin_2022", donor_meta_all.index[donor_meta_all["study"] == "Chaffin 2022"].tolist(), False),
        ("study_Penn_unique9", donor_meta_all.index[donor_meta_all["study"] == "Penn"].tolist(), False),
        ("study_Simonson_2023", donor_meta_all.index[donor_meta_all["study"] == "Simonson 2023"].tolist(), False),
    ]

    sens_results: dict[str, dict[str, Any]] = {}
    for label, sub_donors, adj_study in sens_specs:
        print(f"5. Running {label} (N={len(sub_donors)}, adjust_study={adj_study})...", flush=True)
        sens_results[label] = run_suite_for_subset(
            label, sub_donors, raw_pb_by_group, cell_counts_df, donor_meta_all, p1_program, adj_study, celltype_names
        )

    print("6. Exporting verified CSV artifacts...", flush=True)
    write_clean_csv(audit_df, RESULTS_DIR / f"{PREFIX}-donors.csv")
    write_clean_csv(merged_counts.reset_index(), RESULTS_DIR / f"{PREFIX}-cellcounts.csv")
    write_clean_csv(primary["qc"], RESULTS_DIR / f"{PREFIX}-qc.csv")
    write_clean_csv(primary["fam_a"], RESULTS_DIR / f"{PREFIX}-fractions.csv")
    write_clean_csv(primary["summary"], RESULTS_DIR / f"{PREFIX}-genes-summary.csv")
    top_genes = (
        primary["genes"]
        .sort_values(["pvalue", "group", "gene"])
        .groupby("group", sort=False)
        .head(30)
        .reset_index(drop=True)
    )
    sig_genes = primary["genes"][primary["genes"]["q_bh_all14"] < 0.05].sort_values(["q_bh_all14", "pvalue"])
    top_export = pd.concat([sig_genes, top_genes], ignore_index=True).drop_duplicates(subset=["group", "gene"]).sort_values(["pvalue", "group", "gene"])
    write_clean_csv(top_export, RESULTS_DIR / f"{PREFIX}-top-genes.csv")
    write_clean_csv(primary["conc"], RESULTS_DIR / f"{PREFIX}-concordance.csv")
    write_clean_csv(primary["fam_c"], RESULTS_DIR / f"{PREFIX}-p1-signature.csv")

    sens_a = pd.concat([r["fam_a"] for r in sens_results.values() if not r["fam_a"].empty], ignore_index=True)
    sens_sum = pd.concat([r["summary"] for r in sens_results.values() if not r["summary"].empty], ignore_index=True)
    sens_conc = pd.concat([r["conc"] for r in sens_results.values() if not r["conc"].empty], ignore_index=True)
    sens_c = pd.concat([r["fam_c"] for r in sens_results.values() if not r["fam_c"].empty], ignore_index=True)
    write_clean_csv(sens_a, RESULTS_DIR / f"{PREFIX}-sensitivities-fractions.csv")
    write_clean_csv(sens_sum, RESULTS_DIR / f"{PREFIX}-sensitivities-genes-summary.csv")
    write_clean_csv(sens_conc, RESULTS_DIR / f"{PREFIX}-sensitivities-concordance.csv")
    write_clean_csv(sens_c, RESULTS_DIR / f"{PREFIX}-sensitivities-p1-signature.csv")

    ct_genes_mask = primary["genes"]["group"] != "All_LV_Nuclei"
    conc_q_cols = [c for c in primary["conc"].columns if c.startswith("q_bh")]
    conc_q_col = conc_q_cols[0] if conc_q_cols else None
    meta_summary = {
        "gate_totals": gate_totals,
        "audit_counts": audit_df["disposition_reason"].value_counts().to_dict(),
        "primary_donors": len(all_donors),
        "primary_nuclei": int(keep_mask.sum()),
        "family_a_size": int(len(primary["fam_a"])),
        "family_b_genes_all14": int(len(primary["genes"])),
        "family_b_genes_celltypes13": int(ct_genes_mask.sum()),
        "family_b_sig_all14_q05": int((primary["genes"]["q_bh_all14"] < 0.05).sum()),
        "family_b_sig_celltypes13_q05": int((primary["genes"].loc[ct_genes_mask, "q_bh_celltypes13"] < 0.05).sum()),
        "family_b_concordance_pairs": int(len(primary["conc"])),
        "family_b_concordance_sig_q05": int((primary["conc"][conc_q_col] < 0.05).sum()) if conc_q_col else 0,
        "family_c_size": int(len(primary["fam_c"])),
    }
    (RESULTS_DIR / f"{PREFIX}-metadata-summary.json").write_text(json.dumps(meta_summary, indent=2) + "\n", encoding="utf-8")

    print("\n=== PRIMARY SUMMARY ===", flush=True)
    print(json.dumps(meta_summary, indent=2), flush=True)
    print("\n=== FAMILY A (FRACTIONS) ===", flush=True)
    print(primary["fam_a"].to_string(index=False), flush=True)
    print("\n=== FAMILY B GENE SUMMARY ===", flush=True)
    print(primary["summary"].to_string(index=False), flush=True)
    print("\n=== FAMILY B CONCORDANCE TOP 15 ===", flush=True)
    print(primary["conc"].head(15).to_string(index=False), flush=True)
    print("\n=== FAMILY C (P1 SIGNATURE, 2000 perms) ===", flush=True)
    print(primary["fam_c"].to_string(index=False), flush=True)
    print("\n=== ENCODE v4 STANDALONE FAMILY C ===", flush=True)
    enc_c = sens_results["study_ENCODE_v4_Snyder"]["fam_c"]
    print(enc_c.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
