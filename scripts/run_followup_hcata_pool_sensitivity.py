"""Pool and processing-cohort batch sensitivity checks for HCATA Yazar and Perez caches."""

from __future__ import annotations

import json
import os
from pathlib import Path
import pickle
import sys
import time

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_var] = "1"

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.followup_organs import (  # noqa: E402
    add_bh,
    finalize_pseudobulk_group,
    load_p1_clean_human,
    run_family_a_clr,
    run_family_b_genes_and_concordance,
    run_family_c_p1_signature,
    write_clean_csv,
)

DATA_ROOT = Path("/home/reza/projects/agescape-revision-data")
HCATA_TMP = DATA_ROOT / "tmp_organs" / "hcata"
RESULTS_DIR = REPO_ROOT / "results"
P1_ITERATIONS = 2000
CONCORDANCE_PERMS = 1000
SEED = 20261105
YAZAR_DOUBLETS_LABEL = "peripheral blood mononuclear cell"


def load_cache(cache_dir: Path, drop_count_cols: tuple[str, ...] = ()) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, pd.DataFrame], dict[str, pd.DataFrame]]:
    """Load cached donor metadata, cell counts, and finalized log2CPM pseudobulk matrices."""
    donor_meta = pd.read_csv(cache_dir / "donor_meta.csv")
    donor_meta["donor_id"] = donor_meta["donor_id"].astype(str)
    donor_meta = donor_meta.set_index("donor_id", drop=False)

    cell_counts = pd.read_csv(cache_dir / "cell_counts.csv", index_col=0)
    cell_counts.index = cell_counts.index.astype(str)
    for col in drop_count_cols:
        if col in cell_counts.columns:
            cell_counts = cell_counts.drop(columns=[col])

    with open(cache_dir / "pb_manifest.json", encoding="utf-8") as fh:
        manifest: dict[str, str] = json.load(fh)

    expr_by_group: dict[str, pd.DataFrame] = {}
    meta_by_group: dict[str, pd.DataFrame] = {}
    for group, filename in manifest.items():
        with open(cache_dir / filename, "rb") as pf:
            raw_pb: pd.DataFrame = pickle.load(pf)
        raw_pb.columns = raw_pb.columns.astype(str)
        if group == "All_Cells":
            group_counts = cell_counts.loc[raw_pb.columns].sum(axis=1)
        else:
            group_counts = cell_counts.loc[raw_pb.columns, group]
        log2cpm, m_df, _ = finalize_pseudobulk_group(raw_pb, donor_meta, group_counts)
        expr_by_group[group] = log2cpm
        meta_by_group[group] = m_df
    return donor_meta, cell_counts, expr_by_group, meta_by_group


def add_yazar_pool_dummies(df: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """Attach pool_number fixed-effect indicator columns ordered smallest-to-largest (largest omitted naturally or as ref)."""
    out = df.copy()
    vc = out["pool_number"].value_counts()
    pools_sorted = [p for p, cnt in sorted(vc.items(), key=lambda kv: (-int(kv[1]), int(kv[0])))]
    ref_pool = pools_sorted[0]
    dummy_cols: list[str] = []
    for p in pools_sorted[1:]:
        col = f"pool_{int(p)}"
        out[col] = (out["pool_number"] == p).astype(float)
        dummy_cols.append(col)
    out.attrs["ref_pool"] = int(ref_pool)
    return out, dummy_cols


def run_yazar_pool_fe(p1_program: dict[str, list[str]]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Run Yazar OneK1K PBMC adjusted for sex_male, log_lib, and 75-pool fixed effects."""
    donor_meta, cell_counts, expr_by_group, meta_by_group = load_cache(
        HCATA_TMP / "yazar_cache", drop_count_cols=(YAZAR_DOUBLETS_LABEL,)
    )
    donor_clr, clr_pool_cols = add_yazar_pool_dummies(donor_meta)
    frac_df = run_family_a_clr(
        cell_counts,
        donor_clr,
        age_col="age_per_decade",
        candidate_covariates=("sex_male", *clr_pool_cols),
        min_donors_present=3,
        min_total_cells=20,
        min_residual_df=10,
    )
    frac_df["study"] = "yazar"
    frac_df["arm"] = "pool_fe"

    # For within-pool gene & P1 models, restrict each group to donors in pools with >= 2 donors in that group
    # and require >= 10 residual df after intercept + age + sex + log_lib + (n_pools - 1) dummies.
    all_pools = sorted(donor_meta["pool_number"].dropna().astype(int).unique().tolist())
    all_pool_cols = [f"pool_{p}" for p in all_pools]

    fe_expr: dict[str, pd.DataFrame] = {}
    fe_meta: dict[str, pd.DataFrame] = {}
    excluded_rows: list[dict[str, object]] = []

    for group, m_df in meta_by_group.items():
        vc = m_df["pool_number"].value_counts()
        ge2_pools = vc[vc >= 2]
        keep_donors = m_df.index[m_df["pool_number"].isin(ge2_pools.index)].tolist()
        n_ge2 = len(keep_donors)
        n_pools_ge2 = len(ge2_pools)
        # Exact residual df estimate: intercept + age + sex + log_lib + (n_pools_ge2 - 1) dummies
        est_df = n_ge2 - (4 + max(0, n_pools_ge2 - 1))
        if n_pools_ge2 < 2 or est_df < 10:
            excluded_rows.append(
                {
                    "study": "yazar",
                    "arm": "pool_fe",
                    "group": group,
                    "n_donors_original": len(m_df),
                    "n_donors_ge2_pools": n_ge2,
                    "n_pools_ge2": n_pools_ge2,
                    "estimated_residual_df": est_df,
                    "status": "not_evaluable_within_pool_df",
                }
            )
            continue
        sub_meta = m_df.loc[keep_donors].copy()
        # Keep the largest pool in this group zeroed in reference by omitting it via collinear check or explicit dummy set
        ref_p = int(ge2_pools.index[0])
        for p in all_pools:
            sub_meta[f"pool_{p}"] = (sub_meta["pool_number"].astype(int) == p).astype(float) if p != ref_p else 0.0
        fe_meta[group] = sub_meta
        fe_expr[group] = expr_by_group[group][keep_donors]
        excluded_rows.append(
            {
                "study": "yazar",
                "arm": "pool_fe",
                "group": group,
                "n_donors_original": len(m_df),
                "n_donors_ge2_pools": n_ge2,
                "n_pools_ge2": n_pools_ge2,
                "estimated_residual_df": est_df,
                "status": "evaluated",
            }
        )

    cov_tuple = ("sex_male", "log_lib", *all_pool_cols)
    cell_groups = {g: fe_expr[g] for g in fe_expr if g != "All_Cells"}
    cell_metas = {g: fe_meta[g] for g in fe_meta if g != "All_Cells"}

    genes_all, gene_sum_cell, conc_df = run_family_b_genes_and_concordance(
        cell_groups,
        cell_metas,
        age_col="age_per_decade",
        candidate_covariates=cov_tuple,
        min_donors=4,
        min_residual_df=10,
        is_unbiased_single_cell=True,
        concordance_perms=CONCORDANCE_PERMS,
        seed=SEED,
    )
    genes_pool, gene_sum_pool, _ = run_family_b_genes_and_concordance(
        {"All_Cells": fe_expr["All_Cells"]},
        {"All_Cells": fe_meta["All_Cells"]},
        age_col="age_per_decade",
        candidate_covariates=cov_tuple,
        min_donors=4,
        min_residual_df=10,
        is_unbiased_single_cell=True,
        concordance_perms=0,
        seed=SEED,
    )
    gene_sum_cell["is_pooled"] = False
    gene_sum_pool["is_pooled"] = True
    gene_summary = pd.concat([gene_sum_cell, gene_sum_pool], ignore_index=True)
    gene_summary["study"] = "yazar"
    gene_summary["arm"] = "pool_fe"

    # Drop verbose pool dummy names in csv summary to keep concise, recording count instead
    gene_summary["n_pool_dummies_active"] = gene_summary["covariates"].astype(str).str.count("pool_")
    frac_df["n_pool_dummies_active"] = frac_df["covariates"].astype(str).str.count("pool_")

    p1_cell = run_family_c_p1_signature(
        cell_groups,
        cell_metas,
        p1_program=p1_program,
        age_col="age_per_decade",
        candidate_covariates=cov_tuple,
        min_donors=4,
        min_residual_df=10,
        iterations=P1_ITERATIONS,
        seed=SEED,
    )
    p1_cell["is_pooled"] = False
    p1_pool = run_family_c_p1_signature(
        {"All_Cells": fe_expr["All_Cells"]},
        {"All_Cells": fe_meta["All_Cells"]},
        p1_program=p1_program,
        age_col="age_per_decade",
        candidate_covariates=cov_tuple,
        min_donors=4,
        min_residual_df=10,
        iterations=P1_ITERATIONS,
        seed=SEED,
    )
    p1_pool["is_pooled"] = True

    p1_all_groups = pd.concat([p1_cell, p1_pool], ignore_index=True)
    p1_all_groups = add_bh(p1_all_groups, ["composite_pvalue"], p_column="composite_pvalue", q_column="q_bh_all_composite")
    p1_all_groups = add_bh(p1_all_groups, ["empirical_null_pvalue"], p_column="empirical_null_pvalue", q_column="q_bh_all_empirical")
    p1_all_groups["family_c_all_groups_size"] = len(p1_all_groups)
    p1_all_groups["study"] = "yazar"
    p1_all_groups["arm"] = "pool_fe"
    p1_all_groups["n_pool_dummies_active"] = p1_all_groups["covariates"].astype(str).str.count("pool_")

    conc_df["study"] = "yazar"
    conc_df["arm"] = "pool_fe"
    elig_df = pd.DataFrame(excluded_rows)

    # Return also compact top genes (top 50 cell + top 20 pooled) for inspection
    if not genes_all.empty:
        genes_all = genes_all.sort_values("pvalue").head(50).copy()
        genes_all["study"] = "yazar"
        genes_all["arm"] = "pool_fe"
    return frac_df, gene_summary, conc_df, p1_all_groups, elig_df


def add_perez_cohort_columns(df: pd.DataFrame) -> tuple[pd.DataFrame, list[str], list[str]]:
    """Attach both multi-hot Processing_Cohort indicators and exact level dummies (vs '4.0')."""
    out = df.copy()
    tokens = out["Processing_Cohort"].fillna("").astype(str).apply(lambda s: [t.strip() for t in s.split(";") if t.strip()])
    for k in ("1.0", "2.0", "3.0", "4.0"):
        key = k.replace(".", "_")
        out[f"cohort_{key}"] = tokens.apply(lambda ts, kk=k: float(kk in ts))
    multihot_cols = ["cohort_1_0", "cohort_2_0", "cohort_3_0", "cohort_4_0"]

    levels = ["1.0", "1.0;2.0", "1.0;2.0;3.0", "3.0", "3.0;4.0", "2.0;3.0"]
    level_cols: list[str] = []
    for lv in levels:
        cname = "lvl_" + lv.replace(".", "").replace(";", "_")
        out[cname] = (out["Processing_Cohort"].astype(str) == lv).astype(float)
        level_cols.append(cname)
    return out, multihot_cols, level_cols


def run_perez_cohort_sensitivity(p1_program: dict[str, list[str]]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Run Perez PBMC adjusted for Processing_Cohort (multi-hot indicators AND exact levels)."""
    donor_meta, cell_counts, expr_by_group, meta_by_group = load_cache(HCATA_TMP / "perez_cache")
    donor_clr, multihot_cols, level_cols = add_perez_cohort_columns(donor_meta)

    frac_list = []
    gsum_list = []
    conc_list = []
    p1_list = []

    arms = [
        ("cohort_multihot", False, multihot_cols),
        ("cohort_level_ge2", True, level_cols),
    ]
    for arm_name, drop_singleton_lvl, batch_cols in arms:
        if drop_singleton_lvl:
            keep_ids = donor_clr.index[donor_clr["Processing_Cohort"].astype(str) != "2.0"].tolist()
            d_clr = donor_clr.loc[keep_ids]
            c_clr = cell_counts.loc[keep_ids]
        else:
            d_clr = donor_clr
            c_clr = cell_counts

        fr = run_family_a_clr(
            c_clr,
            d_clr,
            age_col="age_per_decade",
            candidate_covariates=("sex_male", *batch_cols),
            min_donors_present=3,
            min_total_cells=20,
            min_residual_df=10,
        )
        fr["study"] = "perez"
        fr["arm"] = arm_name
        frac_list.append(fr)

        arm_expr: dict[str, pd.DataFrame] = {}
        arm_meta: dict[str, pd.DataFrame] = {}
        for g, m_df in meta_by_group.items():
            m_aug, _, _ = add_perez_cohort_columns(m_df)
            if drop_singleton_lvl:
                m_aug = m_aug.loc[m_aug["Processing_Cohort"].astype(str) != "2.0"]
            arm_meta[g] = m_aug
            arm_expr[g] = expr_by_group[g][list(m_aug.index)]

        cov_tuple = ("sex_male", "log_lib", *batch_cols)
        cell_g = {g: arm_expr[g] for g in arm_expr if g != "All_Cells"}
        cell_m = {g: arm_meta[g] for g in arm_meta if g != "All_Cells"}
        _, gsum_cell, conc = run_family_b_genes_and_concordance(
            cell_g,
            cell_m,
            age_col="age_per_decade",
            candidate_covariates=cov_tuple,
            min_donors=4,
            min_residual_df=10,
            is_unbiased_single_cell=True,
            concordance_perms=CONCORDANCE_PERMS,
            seed=SEED,
        )
        _, gsum_pool, _ = run_family_b_genes_and_concordance(
            {"All_Cells": arm_expr["All_Cells"]},
            {"All_Cells": arm_meta["All_Cells"]},
            age_col="age_per_decade",
            candidate_covariates=cov_tuple,
            min_donors=4,
            min_residual_df=10,
            is_unbiased_single_cell=True,
            concordance_perms=0,
            seed=SEED,
        )
        gsum_cell["is_pooled"] = False
        gsum_pool["is_pooled"] = True
        gsum = pd.concat([gsum_cell, gsum_pool], ignore_index=True)
        gsum["study"] = "perez"
        gsum["arm"] = arm_name
        gsum_list.append(gsum)

        conc["study"] = "perez"
        conc["arm"] = arm_name
        conc_list.append(conc)

        p1_c = run_family_c_p1_signature(
            cell_g,
            cell_m,
            p1_program=p1_program,
            age_col="age_per_decade",
            candidate_covariates=cov_tuple,
            min_donors=4,
            min_residual_df=10,
            iterations=P1_ITERATIONS,
            seed=SEED,
        )
        p1_c["is_pooled"] = False
        p1_p = run_family_c_p1_signature(
            {"All_Cells": arm_expr["All_Cells"]},
            {"All_Cells": arm_meta["All_Cells"]},
            p1_program=p1_program,
            age_col="age_per_decade",
            candidate_covariates=cov_tuple,
            min_donors=4,
            min_residual_df=10,
            iterations=P1_ITERATIONS,
            seed=SEED,
        )
        p1_p["is_pooled"] = True
        p1_df = pd.concat([p1_c, p1_p], ignore_index=True)
        p1_df = add_bh(p1_df, ["composite_pvalue"], p_column="composite_pvalue", q_column="q_bh_all_composite")
        p1_df = add_bh(p1_df, ["empirical_null_pvalue"], p_column="empirical_null_pvalue", q_column="q_bh_all_empirical")
        p1_df["family_c_all_groups_size"] = len(p1_df)
        p1_df["study"] = "perez"
        p1_df["arm"] = arm_name
        p1_list.append(p1_df)

    return (
        pd.concat(frac_list, ignore_index=True),
        pd.concat(gsum_list, ignore_index=True),
        pd.concat(conc_list, ignore_index=True),
        pd.concat(p1_list, ignore_index=True),
    )


def main() -> int:
    t0 = time.monotonic()
    p1_program = load_p1_clean_human(REPO_ROOT)

    y_frac, y_gsum, y_conc, y_p1, y_elig = run_yazar_pool_fe(p1_program)
    p_frac, p_gsum, p_conc, p_p1 = run_perez_cohort_sensitivity(p1_program)

    # Condense verbose pool covariate strings in output CSVs while recording exact active dummy count
    for df in (y_frac, y_gsum, y_p1):
        df["covariates"] = df.apply(
            lambda r: ";".join([c for c in str(r["covariates"]).split(";") if not c.startswith("pool_")])
            + f";pool_fe({int(r['n_pool_dummies_active'])}_dummies)",
            axis=1,
        )
        df["omitted_covariates"] = df["omitted_covariates"].apply(
            lambda s: ";".join([c for c in str(s).split(";") if not c.startswith("pool_")]) or "ref_or_unused_pools_only"
        )

    write_clean_csv(pd.concat([y_frac, p_frac], ignore_index=True), RESULTS_DIR / "followup-hcata-pool-sensitivity-fractions.csv")
    write_clean_csv(pd.concat([y_gsum, p_gsum], ignore_index=True), RESULTS_DIR / "followup-hcata-pool-sensitivity-genes-summary.csv")
    write_clean_csv(pd.concat([y_conc, p_conc], ignore_index=True), RESULTS_DIR / "followup-hcata-pool-sensitivity-concordance.csv")
    write_clean_csv(pd.concat([y_p1, p_p1], ignore_index=True), RESULTS_DIR / "followup-hcata-pool-sensitivity-p1-signature.csv")
    write_clean_csv(y_elig, RESULTS_DIR / "followup-hcata-yazar-pool-fe-eligibility.csv")

    elapsed = time.monotonic() - t0
    print(f"HCATA pool/cohort sensitivity complete in {elapsed:.1f}s")
    print("Yazar pool_fe P1 passing q_bh_empirical < 0.05 (cell types):",
          int(((~y_p1["is_pooled"]) & (y_p1["q_bh_empirical"] < 0.05)).sum()),
          "of", int((~y_p1["is_pooled"]).sum()))
    print("Yazar pool_fe CLR passing q < 0.05:", int((y_frac["q_bh_family_a"] < 0.05).sum()), "of", len(y_frac))
    print("Yazar pool_fe concordance:", y_conc[["mean_pairwise_pearson_r", "median_pairwise_pearson_r", "perm_p_two_sided"]].to_dict("records"))
    print("Perez cohort sensitivities P1 passing:",
          p_p1.groupby("arm").apply(lambda d: int(((~d["is_pooled"]) & (d["q_bh_empirical"] < 0.05)).sum()), include_groups=False).to_dict())
    print("Perez cohort concordance:", p_conc[["arm", "mean_pairwise_pearson_r", "median_pairwise_pearson_r", "perm_p_two_sided"]].to_dict("records"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
