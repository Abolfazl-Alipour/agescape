#!/usr/bin/env python3
"""Item #17: Run HCATA newly identified open studies (Guilliams liver, Yazar OneK1K PBMC, Perez control PBMC)."""

from __future__ import annotations

import gc
import json
from pathlib import Path
import pickle
import sys
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats
from statsmodels.stats.multitest import multipletests

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.followup_muscle import NULL_SEED
from src.followup_organs import (
    finalize_pseudobulk_group,
    load_p1_clean_human,
    run_family_a_clr,
    run_family_b_genes_and_concordance,
    run_family_c_p1_signature,
    write_clean_csv,
)

CACHE_ROOT = REPO_ROOT.parent / "agescape-revision-data" / "tmp_organs" / "hcata"
RESULTS_DIR = REPO_ROOT / "results"
MIN_DONORS = 6
MIN_CELLS_PER_DONOR = 20
P1_ITERATIONS = 2000
CONC_PERMS = 1000


def load_study_inputs(
    study: str,
    trim_for_large_n: bool = True,
) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame, dict[str, pd.DataFrame], dict[str, pd.DataFrame], pd.DataFrame]:
    s_dir = CACHE_ROOT / f"{study}_cache"
    audit = json.loads((s_dir / "obs_audit.json").read_text(encoding="utf-8"))
    donor_meta = pd.read_csv(s_dir / "donor_meta.csv", index_col=0)
    donor_meta.index = donor_meta.index.astype(str)
    cell_counts = pd.read_csv(s_dir / "cell_counts.csv", index_col=0)
    cell_counts.index = cell_counts.index.astype(str)
    pb_files = json.loads((s_dir / "pb_manifest.json").read_text(encoding="utf-8"))

    cd45_counts = None
    cd45_path = s_dir / "cd45pos_counts.csv"
    if cd45_path.exists():
        cd45_counts = pd.read_csv(cd45_path, index_col=0)["n_cells"]
        cd45_counts.index = cd45_counts.index.astype(str)

    expr_by_group: dict[str, pd.DataFrame] = {}
    meta_by_group: dict[str, pd.DataFrame] = {}
    qc_rows: list[dict[str, Any]] = []
    thresh = float(np.log2(2.0))  # log2(1.0 + 1)

    for grp, fname in pb_files.items():
        with open(s_dir / fname, "rb") as pf:
            raw_df: pd.DataFrame = pickle.load(pf)
        raw_df.columns = raw_df.columns.astype(str)
        if grp == "All_Cells":
            c_ser = donor_meta["n_cells"].astype(int)
        elif grp == "All_CD45pos_Cells" and cd45_counts is not None:
            c_ser = cd45_counts
        else:
            c_ser = cell_counts[grp].astype(int)

        elig_donors = [d for d in donor_meta.index if d in raw_df.columns and int(c_ser.get(d, 0)) >= MIN_CELLS_PER_DONOR]
        if len(elig_donors) < MIN_DONORS:
            continue
        log2cpm, grp_meta, stats = finalize_pseudobulk_group(
            raw_df[elig_donors], donor_meta.loc[elig_donors], c_ser.loc[elig_donors]
        )
        # Retain >=35%-expressed genes for large cohorts (N >= 50) to cap peak memory well below 4 GB
        # while remaining strictly looser than the >=50% filter inside run_family_b / run_family_c.
        if trim_for_large_n and len(elig_donors) >= 50:
            min_keep = max(3, int(np.floor(0.35 * len(elig_donors))))
            keep_mask = (log2cpm.to_numpy(dtype=float) >= thresh).sum(axis=1) >= min_keep
            log2cpm = log2cpm.loc[keep_mask].copy()

        expr_by_group[grp] = log2cpm
        meta_by_group[grp] = grp_meta
        qc_rows.append(
            {
                "study": study,
                "group": grp,
                "is_pooled": grp.startswith("All_"),
                "n_donors": int(len(elig_donors)),
                "age_min": float(grp_meta["age_years"].min()),
                "age_max": float(grp_meta["age_years"].max()),
                "n_male": int(grp_meta["sex_male"].sum()),
                "n_female": int((1 - grp_meta["sex_male"]).sum()),
                "total_cells_in_group": int(c_ser.loc[elig_donors].sum()),
                "median_cells_per_donor": float(c_ser.loc[elig_donors].median()),
                **{k: v for k, v in stats.items() if k != "count_gate"},
            }
        )
        del raw_df
        gc.collect()

    return audit, donor_meta, cell_counts, expr_by_group, meta_by_group, pd.DataFrame(qc_rows)


def extract_top_genes(genes_df: pd.DataFrame, max_per_group: int = 30) -> pd.DataFrame:
    if genes_df.empty:
        return genes_df
    top = (
        genes_df.sort_values(["pvalue", "group", "gene"], ascending=[True, True, True])
        .groupby("group", sort=False)
        .head(max_per_group)
        .reset_index(drop=True)
    )
    return top


def summarize_study_genes(
    study: str,
    expr_by_group: dict[str, pd.DataFrame],
    meta_by_group: dict[str, pd.DataFrame],
    is_unbiased: bool,
    covariates: tuple[str, ...] = ("sex_male", "log_lib"),
    run_concordance: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    ct_expr = {k: v for k, v in expr_by_group.items() if not k.startswith("All_")}
    ct_meta = {k: meta_by_group[k] for k in ct_expr}
    pool_expr = {k: v for k, v in expr_by_group.items() if k.startswith("All_")}
    pool_meta = {k: meta_by_group[k] for k in pool_expr}

    genes_ct, sum_ct, conc_df = run_family_b_genes_and_concordance(
        ct_expr,
        ct_meta,
        age_col="age_per_decade",
        candidate_covariates=covariates,
        min_donors=MIN_DONORS,
        min_cpm=1.0,
        min_expr_frac=0.5,
        min_residual_df=2,
        is_unbiased_single_cell=is_unbiased and run_concordance,
        concordance_perms=CONC_PERMS,
        seed=NULL_SEED,
    )
    genes_pool, sum_pool, _ = run_family_b_genes_and_concordance(
        pool_expr,
        pool_meta,
        age_col="age_per_decade",
        candidate_covariates=covariates,
        min_donors=MIN_DONORS,
        min_cpm=1.0,
        min_expr_frac=0.5,
        min_residual_df=2,
        is_unbiased_single_cell=False,
        seed=NULL_SEED,
    )
    if not genes_ct.empty:
        genes_ct["is_pooled"] = False
    if not genes_pool.empty:
        genes_pool["is_pooled"] = True
        # Keep cell-type global BH strictly separate from pooled runs; add combined BH for completeness
        genes_pool["q_bh_pool_only"] = genes_pool["q_bh_global"]
        genes_pool["q_bh_global"] = np.nan
    genes_all = pd.concat([genes_ct, genes_pool], ignore_index=True)
    if not genes_all.empty:
        genes_all.insert(0, "study", study)
        genes_all["family_b_all_including_pooled_size"] = int(len(genes_all))
        genes_all["q_bh_all_including_pooled"] = multipletests(genes_all["pvalue"], method="fdr_bh")[1]

    if not sum_ct.empty:
        sum_ct["is_pooled"] = False
    if not sum_pool.empty:
        sum_pool["is_pooled"] = True
    sum_all = pd.concat([sum_ct, sum_pool], ignore_index=True)
    if not sum_all.empty:
        sum_all.insert(0, "study", study)
    if not conc_df.empty:
        conc_df.insert(0, "study", study)
    return genes_all, sum_all, conc_df, genes_ct


def run_p1_with_both_families(
    study: str,
    expr_by_group: dict[str, pd.DataFrame],
    meta_by_group: dict[str, pd.DataFrame],
    p1_program: dict[str, list[str]],
    covariates: tuple[str, ...] = ("sex_male", "log_lib"),
) -> pd.DataFrame:
    ct_expr = {k: v for k, v in expr_by_group.items() if not k.startswith("All_")}
    ct_meta = {k: meta_by_group[k] for k in ct_expr}
    pool_expr = {k: v for k, v in expr_by_group.items() if k.startswith("All_")}
    pool_meta = {k: meta_by_group[k] for k in pool_expr}

    p1_ct = run_family_c_p1_signature(
        ct_expr,
        ct_meta,
        p1_program=p1_program,
        age_col="age_per_decade",
        candidate_covariates=covariates,
        min_donors=MIN_DONORS,
        min_cpm=1.0,
        min_expr_frac=0.5,
        min_residual_df=2,
        iterations=P1_ITERATIONS,
        seed=NULL_SEED,
    )
    if not p1_ct.empty:
        p1_ct["is_pooled"] = False
    p1_pool = run_family_c_p1_signature(
        pool_expr,
        pool_meta,
        p1_program=p1_program,
        age_col="age_per_decade",
        candidate_covariates=covariates,
        min_donors=MIN_DONORS,
        min_cpm=1.0,
        min_expr_frac=0.5,
        min_residual_df=2,
        iterations=P1_ITERATIONS,
        seed=NULL_SEED,
    )
    if not p1_pool.empty:
        p1_pool["is_pooled"] = True
        p1_pool["q_bh_composite"] = np.nan
        p1_pool["q_bh_empirical"] = np.nan

    combined = pd.concat([p1_ct, p1_pool], ignore_index=True)
    if not combined.empty:
        combined.insert(0, "study", study)
        combined["family_c_all_groups_size"] = int(len(combined))
        combined["q_bh_all_composite"] = multipletests(combined["composite_pvalue"], method="fdr_bh")[1]
        combined["q_bh_all_empirical"] = multipletests(combined["empirical_null_pvalue"], method="fdr_bh")[1]
    return combined


def subset_expr_meta(
    expr_by_group: dict[str, pd.DataFrame],
    meta_by_group: dict[str, pd.DataFrame],
    keep_donors: set[str],
) -> tuple[dict[str, pd.DataFrame], dict[str, pd.DataFrame]]:
    e_sub: dict[str, pd.DataFrame] = {}
    m_sub: dict[str, pd.DataFrame] = {}
    for g, mat in expr_by_group.items():
        cols = [d for d in mat.columns if d in keep_donors]
        if len(cols) >= MIN_DONORS:
            e_sub[g] = mat[cols]
            m_sub[g] = meta_by_group[g].loc[cols]
    return e_sub, m_sub


def run_guilliams(p1_program: dict[str, list[str]]) -> dict[str, Any]:
    audit, donor_meta, _, expr_by_group, meta_by_group, qc_df = load_study_inputs("guilliams", trim_for_large_n=False)
    donor_meta["is_citeseq"] = donor_meta["typeSample"].astype(str).str.contains("citeSeq").astype(float)
    for g in meta_by_group:
        meta_by_group[g]["is_citeseq"] = donor_meta.loc[meta_by_group[g].index, "is_citeseq"]

    genes_all, sum_all, _, genes_ct = summarize_study_genes(
        "guilliams", expr_by_group, meta_by_group, is_unbiased=False
    )
    p1_df = run_p1_with_both_families("guilliams", expr_by_group, meta_by_group, p1_program)

    # Sensitivities for Guilliams FACS/surgical cohort
    sens_genes_rows = []
    sens_p1_frames = []
    no_malig = set(donor_meta.index) - {"H16", "H21"}
    no_tumor = set(donor_meta.index) - {"H13", "H16", "H21"}
    female_ids = set(donor_meta[donor_meta["sex_male"] == 0].index)

    specs = [
        ("no_malignant_resection_n8", no_malig, ("sex_male", "log_lib")),
        ("no_any_tumor_n7", no_tumor, ("log_lib",)),
        ("female_only_n8", female_ids, ("log_lib",)),
        ("adjust_citeseq_n10", set(donor_meta.index), ("sex_male", "log_lib", "is_citeseq")),
    ]
    for sname, dset, covs in specs:
        e_sub, m_sub = subset_expr_meta(expr_by_group, meta_by_group, dset)
        g_s, sum_s, _, g_ct_s = summarize_study_genes(
            f"guilliams_{sname}", e_sub, m_sub, is_unbiased=False, covariates=covs
        )
        p1_s = run_p1_with_both_families(f"guilliams_{sname}", e_sub, m_sub, p1_program, covariates=covs)
        if not sum_s.empty:
            sum_s.insert(0, "sensitivity", sname)
            sens_genes_rows.append(sum_s)
        if not p1_s.empty:
            p1_s.insert(0, "sensitivity", sname)
            sens_p1_frames.append(p1_s)

    write_clean_csv(donor_meta.reset_index(), RESULTS_DIR / "followup-hcata-guilliams-donors.csv")
    write_clean_csv(qc_df, RESULTS_DIR / "followup-hcata-guilliams-qc.csv")
    write_clean_csv(sum_all, RESULTS_DIR / "followup-hcata-guilliams-genes-summary.csv")
    write_clean_csv(extract_top_genes(genes_all, 30), RESULTS_DIR / "followup-hcata-guilliams-top-genes.csv")
    write_clean_csv(p1_df, RESULTS_DIR / "followup-hcata-guilliams-p1-signature.csv")
    if sens_genes_rows:
        write_clean_csv(pd.concat(sens_genes_rows, ignore_index=True), RESULTS_DIR / "followup-hcata-guilliams-sensitivities-genes.csv")
    if sens_p1_frames:
        write_clean_csv(pd.concat(sens_p1_frames, ignore_index=True), RESULTS_DIR / "followup-hcata-guilliams-sensitivities-p1.csv")

    return {
        "study": "guilliams",
        "donors": int(len(donor_meta)),
        "age_min": float(donor_meta["age_years"].min()),
        "age_max": float(donor_meta["age_years"].max()),
        "n_male": int(donor_meta["sex_male"].sum()),
        "n_female": int((1 - donor_meta["sex_male"]).sum()),
        "sampling": "FACS_sorted_CD45pos_surgical",
        "gate": audit["gate"],
        "n_ct_tests": int(len(genes_ct)),
        "n_ct_sig_global_q05": int((genes_ct["q_bh_global"] < 0.05).sum()) if not genes_ct.empty else 0,
        "min_ct_q": float(genes_ct["q_bh_global"].min()) if not genes_ct.empty else np.nan,
        "p1_ct_sig_comp_q05": int((p1_df[~p1_df["is_pooled"]]["q_bh_composite"] < 0.05).sum()) if not p1_df.empty else 0,
        "p1_ct_sig_emp_q05": int((p1_df[~p1_df["is_pooled"]]["q_bh_empirical"] < 0.05).sum()) if not p1_df.empty else 0,
    }


def check_pool_age_confounding(donor_meta: pd.DataFrame, col: str) -> dict[str, Any]:
    if col not in donor_meta.columns:
        return {"column": col, "present": False}
    pools = donor_meta[col].astype(str)
    valid = ~pools.str.contains(";")
    sub = donor_meta.loc[valid]
    if len(sub) < 10 or pools.loc[valid].nunique() < 2:
        return {"column": col, "present": True, "single_pool_donors": int(valid.sum()), "n_pools": int(pools.loc[valid].nunique())}
    groups = [g["age_years"].to_numpy(dtype=float) for _, g in sub.groupby(pools.loc[valid])]
    groups = [g for g in groups if len(g) >= 2]
    fstat, pval = scipy_stats.f_oneway(*groups)
    # eta^2 = SS_between / SS_total
    y = np.concatenate(groups)
    gm = y.mean()
    ss_tot = float(np.sum((y - gm) ** 2))
    ss_bet = float(sum(len(g) * (g.mean() - gm) ** 2 for g in groups))
    return {
        "column": col,
        "present": True,
        "single_pool_donors": int(valid.sum()),
        "n_pools": int(len(groups)),
        "anova_r2": ss_bet / ss_tot if ss_tot > 0 else 0.0,
        "anova_pvalue": float(pval),
    }


def run_unbiased_pbmc_study(
    study: str,
    p1_program: dict[str, list[str]],
    sens_name: str,
    sens_keep_fn,
    sens_covs: tuple[str, ...],
) -> dict[str, Any]:
    audit, donor_meta, cell_counts, expr_by_group, meta_by_group, qc_df = load_study_inputs(study, trim_for_large_n=True)
    # In Yazar OneK1K, Azimuth 'Doublet' (131 cells) maps to 'peripheral blood mononuclear cell'; exclude explicit doublets from Family A
    fam_a_counts = (
        cell_counts.drop(columns=["peripheral blood mononuclear cell"], errors="ignore")
        if study == "yazar"
        else cell_counts
    )
    # Primary Family A (CLR on cell types present in >=6 donors and >=100 cells)
    fam_a = run_family_a_clr(
        fam_a_counts,
        donor_meta,
        age_col="age_per_decade",
        candidate_covariates=("sex_male",),
        min_donors_present=6,
        min_total_cells=100,
        min_residual_df=2,
    )
    if not fam_a.empty:
        fam_a.insert(0, "arm", "primary")
        fam_a.insert(0, "study", study)

    genes_all, sum_all, conc_df, genes_ct = summarize_study_genes(
        study, expr_by_group, meta_by_group, is_unbiased=True, covariates=("sex_male", "log_lib"), run_concordance=True
    )
    if not conc_df.empty:
        conc_df.insert(1, "arm", "primary")
    p1_df = run_p1_with_both_families(study, expr_by_group, meta_by_group, p1_program, covariates=("sex_male", "log_lib"))
    if not p1_df.empty:
        p1_df.insert(1, "arm", "primary")

    # Sensitivity arm
    keep_ids = set(sens_keep_fn(donor_meta))
    dmeta_s = donor_meta.loc[donor_meta.index.isin(keep_ids)]
    fam_a_s = run_family_a_clr(
        fam_a_counts.loc[dmeta_s.index],
        dmeta_s,
        age_col="age_per_decade",
        candidate_covariates=("sex_male",) if "sex_male" in sens_covs else (),
        min_donors_present=6,
        min_total_cells=100,
        min_residual_df=2,
    )
    if not fam_a_s.empty:
        fam_a_s.insert(0, "arm", sens_name)
        fam_a_s.insert(0, "study", study)

    e_sub, m_sub = subset_expr_meta(expr_by_group, meta_by_group, keep_ids)
    _, sum_s, conc_s, genes_ct_s = summarize_study_genes(
        study, e_sub, m_sub, is_unbiased=True, covariates=sens_covs, run_concordance=True
    )
    if not sum_s.empty:
        sum_s.insert(1, "arm", sens_name)
    if not conc_s.empty:
        conc_s.insert(1, "arm", sens_name)
    p1_s = run_p1_with_both_families(study, e_sub, m_sub, p1_program, covariates=sens_covs)
    if not p1_s.empty:
        p1_s.insert(1, "arm", sens_name)

    frac_combined = pd.concat([fam_a, fam_a_s], ignore_index=True)
    conc_combined = pd.concat([conc_df, conc_s], ignore_index=True)
    p1_combined = pd.concat([p1_df, p1_s], ignore_index=True)
    sum_sens = sum_s if not sum_s.empty else pd.DataFrame()

    write_clean_csv(donor_meta.reset_index(), RESULTS_DIR / f"followup-hcata-{study}-donors.csv")
    write_clean_csv(qc_df, RESULTS_DIR / f"followup-hcata-{study}-qc.csv")
    write_clean_csv(frac_combined, RESULTS_DIR / f"followup-hcata-{study}-fractions.csv")
    write_clean_csv(sum_all, RESULTS_DIR / f"followup-hcata-{study}-genes-summary.csv")
    write_clean_csv(extract_top_genes(genes_all, 30), RESULTS_DIR / f"followup-hcata-{study}-top-genes.csv")
    write_clean_csv(conc_combined, RESULTS_DIR / f"followup-hcata-{study}-concordance.csv")
    write_clean_csv(p1_combined, RESULTS_DIR / f"followup-hcata-{study}-p1-signature.csv")
    if not sum_sens.empty:
        write_clean_csv(sum_sens, RESULTS_DIR / f"followup-hcata-{study}-sensitivity-genes-summary.csv")

    pool_info = check_pool_age_confounding(donor_meta, "pool_number" if study == "yazar" else "Processing_Cohort")

    out = {
        "study": study,
        "donors": int(len(donor_meta)),
        "age_min": float(donor_meta["age_years"].min()),
        "age_max": float(donor_meta["age_years"].max()),
        "n_male": int(donor_meta["sex_male"].sum()),
        "n_female": int((1 - donor_meta["sex_male"]).sum()),
        "sampling": "unbiased_PBMC",
        "gate": audit["gate"],
        "pool_age_check": pool_info,
        "fam_a_tested": int(len(fam_a)),
        "fam_a_sig_q05": int((fam_a["q_bh_family_a"] < 0.05).sum()) if not fam_a.empty else 0,
        "n_ct_tests": int(len(genes_ct)),
        "n_ct_sig_global_q05": int((genes_ct["q_bh_global"] < 0.05).sum()) if not genes_ct.empty else 0,
        "min_ct_q": float(genes_ct["q_bh_global"].min()) if not genes_ct.empty else np.nan,
        "conc_mean_r": float(conc_df["mean_pairwise_pearson_r"].iloc[0]) if not conc_df.empty else np.nan,
        "conc_perm_p_two": float(conc_df["perm_p_two_sided"].iloc[0]) if not conc_df.empty else np.nan,
        "p1_ct_sig_comp_q05": int((p1_df[~p1_df["is_pooled"]]["q_bh_composite"] < 0.05).sum()) if not p1_df.empty else 0,
        "p1_ct_sig_emp_q05": int((p1_df[~p1_df["is_pooled"]]["q_bh_empirical"] < 0.05).sum()) if not p1_df.empty else 0,
        "sens_name": sens_name,
        "sens_donors": int(len(dmeta_s)),
        "sens_fam_a_sig_q05": int((fam_a_s["q_bh_family_a"] < 0.05).sum()) if not fam_a_s.empty else 0,
        "sens_n_ct_sig_global_q05": int((genes_ct_s["q_bh_global"] < 0.05).sum()) if not genes_ct_s.empty else 0,
        "sens_conc_mean_r": float(conc_s["mean_pairwise_pearson_r"].iloc[0]) if not conc_s.empty else np.nan,
        "sens_conc_perm_p_two": float(conc_s["perm_p_two_sided"].iloc[0]) if not conc_s.empty else np.nan,
    }
    del expr_by_group, e_sub, genes_all, genes_ct
    gc.collect()
    return out


def main() -> None:
    p1_program = load_p1_clean_human(REPO_ROOT)
    summaries = []

    print("=== Running Guilliams 2022 Liver ===", flush=True)
    summaries.append(run_guilliams(p1_program))

    print("=== Running Yazar 2022 OneK1K PBMC ===", flush=True)
    summaries.append(
        run_unbiased_pbmc_study(
            "yazar",
            p1_program,
            sens_name="age_ge_20_n976",
            sens_keep_fn=lambda dm: dm.index[dm["age_years"] >= 20.0],
            sens_covs=("sex_male", "log_lib"),
        )
    )

    print("=== Running Perez 2022 Control PBMC ===", flush=True)
    summaries.append(
        run_unbiased_pbmc_study(
            "perez",
            p1_program,
            sens_name="female_only_n97",
            sens_keep_fn=lambda dm: dm.index[dm["sex_male"] == 0],
            sens_covs=("log_lib",),
        )
    )

    ov_rows = []
    for s in summaries:
        row = {k: v for k, v in s.items() if k not in {"gate", "pool_age_check"}}
        ov_rows.append(row)
    write_clean_csv(pd.DataFrame(ov_rows), RESULTS_DIR / "followup-hcata-studies-overview.csv")
    (RESULTS_DIR / "followup-hcata-metadata-summary.json").write_text(
        json.dumps(summaries, indent=2), encoding="utf-8"
    )
    print("DONE HCATA Studies:", json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
