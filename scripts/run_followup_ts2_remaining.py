#!/usr/bin/env python3
"""Item #15: Complete 28-tissue TS2 audit and Family B + C analyses across remaining >=6 single-assay tissues."""

from __future__ import annotations

import json
from pathlib import Path
import pickle
import sys
from typing import Any

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.followup_muscle import NULL_SEED
from src.followup_organs import (
    finalize_pseudobulk_group,
    load_p1_clean_human,
    run_family_b_genes_and_concordance,
    run_family_c_p1_signature,
    write_clean_csv,
)

CACHE_ROOT = Path("/home/reza/projects/agescape-revision-data/tmp_organs/ts2_rem")
RESULTS_DIR = REPO_ROOT / "results"
PREFIX = "followup-ts2-remaining"

TARGET_TISSUES = ["Blood", "Bladder", "Fat", "Lymph_Node", "Spleen"]
MIN_DONORS = 6
MIN_CELLS_PER_DONOR = 20
P1_ITERATIONS = 2000

CZ_TISSUE_MAP = {
    "Bladder": ("8ffc335d-36fa-461e-be47-cd6503833420", 3317048769),
    "Blood": ("1bf55482-ad49-446a-87c6-f3442bebcef1", 2568875081),
    "Bone_Marrow": ("4aedb645-7487-4ab0-ad95-ee4db1b8eccd", 1030692363),
    "Ear": ("eef916d5-41c0-434f-a50c-617aa1413987", 181131420),
    "Eye": ("e9b8ea4b-7901-4992-94c7-8e85ccd06fb5", 1384715251),
    "Fat": ("8115812f-9acf-4f62-b15f-bde18086a1c8", 3689061644),
    "Heart": ("82197577-7914-45c0-9fde-6597a1da920f", 1160037664),
    "Kidney": ("76b6cf23-56db-4eac-8ee8-d84ce64c0395", 451232856),
    "Large_Intestine": ("81781bd2-7a71-43a8-ab31-2a9009c1aea6", 1020224367),
    "Liver": ("a4b3a49e-062b-4e3f-8915-02f40607651f", 837059369),
    "Lung": ("fddf43a3-b7c8-49e9-ac0a-5e24b390849f", 3199502267),
    "Lymph_Node": ("6dd28080-1d36-4905-8b44-a9954a4438a2", 4086623268),
    "Mammary": ("33b25f68-a9d6-4643-9283-653e937ae608", 1326569473),
    "Muscle": ("4919a450-62fb-4741-b750-5591443edee8", 2051218332),
    "Ovary": ("6ad78e84-6715-458d-95b4-69e8509cd785", 1817624603),
    "Pancreas": ("ba587dee-e050-4a31-bcb8-05e5418f9086", 740671625),
    "Prostate": ("4a289bf1-97c2-4616-8e1c-585088d28c3e", 904569252),
    "Salivary_Gland": ("a9f49247-9866-422e-a9c3-d315f6e99369", 1806149308),
    "Skin": ("c240d5bc-a521-4330-a76b-7d566e939d87", 504731465),
    "Small_Intestine": ("6db4dfe6-d421-4287-a9dc-22eaa59be966", 1806934696),
    "Spleen": ("5498457f-6a18-4766-b9a3-a2c8c4662703", 2654711797),
    "Stomach": ("3e8bd2b2-fe53-4b95-abb2-ce971b4f1c14", 1072465874),
    "Testis": ("09e3ae14-f355-4ed7-ad81-787b587c8a16", 386411220),
    "Thymus": ("b58baee9-3e90-4491-aa16-eb5818a6db27", 1715061877),
    "Tongue": ("07cf9c5f-6e79-4422-a4de-4427f0f4a76d", 2005012044),
    "Trachea": ("03500d0b-e9ae-42e0-9ae0-576a551099dc", 965745201),
    "Uterus": ("ec18d0d1-3b1b-4413-b011-0078c8d69bd3", 894132831),
    "Vasculature": ("5d54a606-3cd9-4164-b862-a3c51e622790", 2180756396),
}


def build_all28_audit(tissue_gates: dict[str, Any]) -> pd.DataFrame:
    """Reconcile all 28 Tabula Sapiens 2.0 tissues across supplementary workbooks and H5AD checks."""
    org = pd.read_csv(RESULTS_DIR / "followup-ts2-organ-summary.csv")
    cw = pd.read_csv(RESULTS_DIR / "followup-ts2-donor-organ-crosswalk.csv")
    rows = []
    for row in org.itertuples(index=False):
        t = str(row.tissue)
        sub_cw = cw[cw["tissue"] == t]
        n_10x_supp = int((sub_cw["n_cells_10x_supplementary"] > 0).sum())
        n_ss_supp = int((sub_cw["n_cells_smartseq_supplementary"] > 0).sum())
        n_10x_age_verified = int(
            ((sub_cw["n_cells_10x_supplementary"] > 0) & sub_cw["age_years"].notna()).sum()
        )
        ver_id, sz_bytes = CZ_TISSUE_MAP[t]
        n_10x3v3_h5ad = None
        if t == "Bone_Marrow":
            n_10x3v3_h5ad = 5
            status = "excluded_fewer_than_6_single_assay_donors"
            note = "H5AD checked: 7 total donors, 5 on 10x 3' v3 (<6), 5 on 10x 5' v2 (<6), plus LinEnriched/LinDepleted FACS preps"
        elif t == "Blood":
            n_10x3v3_h5ad = 9
            status = "analyzed_non_nk_sorted_immune"
            note = "9 donors on 10x 3' v3 (33-69y); Blood NK already analyzed; remaining 6 non-NK cell types + All_NonNK_Cells analyzed Family B+C only (Enriched/PBMC/Whole/SheelaPrep mixed preps)"
        elif t in {"Bladder", "Fat", "Lymph_Node", "Spleen"}:
            n_10x3v3_h5ad = int(len(tissue_gates[t]["donors_10x3v3"]))
            status = "analyzed_remaining_compartment_sorted"
            note = (
                f"{n_10x3v3_h5ad} donors on 10x 3' v3 ({row.age_min:.0f}-{row.age_max:.0f}y); "
                "Family B+C only (no Family A/concordance: TS microbead/FACS compartment normalization PMC9812260)"
            )
        elif t == "Skin":
            status = "already_analyzed_ts1_and_below_6_donors"
            note = "4 donors (33-61y < 6 floor); prior TS_Skin (4-donor TS1) already analyzed"
        elif t in {"Heart", "Muscle", "Trachea"}:
            status = "excluded_fewer_than_6_10x_donors"
            note = f"6 total donors overall, but only {n_10x_age_verified} donors in 10x (<6 single-assay floor)"
        else:
            status = "excluded_fewer_than_6_donors"
            note = f"{row.n_donors_age_verified} age-verified donors ({n_10x_age_verified} in 10x < 6 floor)"

        rows.append(
            {
                "tissue": t,
                "cz_dataset_version_id": ver_id,
                "cz_h5ad_bytes": sz_bytes,
                "n_donors_supplementary": int(row.n_donors_supplementary),
                "n_donors_age_verified": int(row.n_donors_age_verified),
                "n_donors_10x_age_verified_supp": n_10x_age_verified,
                "n_donors_smartseq_supp": n_ss_supp,
                "n_donors_10x3v3_h5ad_verified": n_10x3v3_h5ad if n_10x3v3_h5ad is not None else "",
                "n_cells_supplementary": int(row.n_cells_supplementary),
                "age_min": row.age_min,
                "age_max": row.age_max,
                "analysis_status": status,
                "disposition_note": note,
            }
        )
    out = pd.DataFrame(rows)
    return out


def load_curated_donor_table() -> pd.DataFrame:
    df = pd.read_csv(RESULTS_DIR / "followup-ts2-donor-summary.csv")
    return df.set_index("donor_id")


def run_single_tissue(tissue: str, curated_donors: pd.DataFrame, p1_program: dict[str, list[str]]) -> dict[str, Any]:
    t_dir = CACHE_ROOT / f"{tissue}_cache"
    audit = json.loads((t_dir / "obs_audit.json").read_text(encoding="utf-8"))
    cell_counts = pd.read_csv(t_dir / "cell_counts_10x3v3.csv", index_col=0)
    with open(t_dir / "raw_pb_by_group.pkl", "rb") as pf:
        raw_pb_by_group: dict[str, pd.DataFrame] = pickle.load(pf)

    # Reconcile H5AD .obs age & sex against curated supplementary table
    donor_rows = []
    for rec in audit["donors_10x3v3"]:
        did = rec["donor_id"]
        c_row = curated_donors.loc[did]
        age_h5ad = int(str(rec["development_stage"]).split("-")[0])
        age_curated = int(c_row["age_years"])
        assert age_h5ad == age_curated, f"{tissue} {did}: H5AD age {age_h5ad} != curated {age_curated}"
        sex_h5ad = str(rec["sex"])
        sex_curated_char = "M" if sex_h5ad == "male" else "F"
        assert sex_curated_char == str(c_row["sex"]), f"{tissue} {did}: sex mismatch"
        donor_rows.append(
            {
                "tissue": tissue,
                "donor_id": did,
                "age_years": float(age_curated),
                "age_per_decade": float(age_curated) / 10.0,
                "sex": sex_h5ad,
                "sex_male": int(sex_h5ad == "male"),
                "n_cells_10x3v3": int(rec["n_cells_10x3v3"]),
                "compartments": ";".join(f"{k}:{v}" for k, v in sorted(rec["compartment_counts"].items())),
                "notes": ";".join(f"{k}:{v}" for k, v in sorted(rec["notes_counts"].items())),
            }
        )
    donor_meta = pd.DataFrame(donor_rows).set_index("donor_id")

    # Select eligible groups (>= MIN_DONORS with >= MIN_CELLS_PER_DONOR); exclude Blood NK
    all_grp_name = "All_NonNK_Cells" if tissue == "Blood" else "All_Cells"
    if tissue == "Blood" and "natural killer cell" in cell_counts.columns:
        non_nk_counts = cell_counts.drop(columns=["natural killer cell"]).sum(axis=1)
    else:
        non_nk_counts = cell_counts.sum(axis=1)

    eligible_groups: list[str] = []
    for ct in sorted(cell_counts.columns):
        if tissue == "Blood" and ct == "natural killer cell":
            continue
        if int((cell_counts[ct] >= MIN_CELLS_PER_DONOR).sum()) >= MIN_DONORS:
            eligible_groups.append(ct)
    if int((non_nk_counts >= MIN_CELLS_PER_DONOR).sum()) >= MIN_DONORS:
        eligible_groups.append(all_grp_name)

    expr_ct: dict[str, pd.DataFrame] = {}
    meta_ct: dict[str, pd.DataFrame] = {}
    qc_rows: list[dict[str, Any]] = []

    for grp in eligible_groups:
        c_ser = non_nk_counts if grp == all_grp_name else cell_counts[grp]
        elig_donors = [d for d in donor_meta.index if int(c_ser.loc[d]) >= MIN_CELLS_PER_DONOR]
        raw_sub = raw_pb_by_group[grp][elig_donors]
        log2cpm, grp_meta, stats = finalize_pseudobulk_group(raw_sub, donor_meta.loc[elig_donors], c_ser.loc[elig_donors])
        expr_ct[grp] = log2cpm
        meta_ct[grp] = grp_meta
        qc_rows.append(
            {
                "tissue": tissue,
                "group": grp,
                "is_pooled_tissue": grp == all_grp_name,
                "n_donors": len(elig_donors),
                "age_min": float(grp_meta["age_years"].min()),
                "age_max": float(grp_meta["age_years"].max()),
                "n_male": int(grp_meta["sex_male"].sum()),
                "n_female": int((1 - grp_meta["sex_male"]).sum()),
                "total_cells_in_group": int(c_ser.loc[elig_donors].sum()),
                "median_cells_per_donor": float(c_ser.loc[elig_donors].median()),
                **{k: v for k, v in stats.items() if k != "count_gate"},
            }
        )

    expr_only_ct = {k: v for k, v in expr_ct.items() if k != all_grp_name}
    meta_only_ct = {k: v for k, v in meta_ct.items() if k != all_grp_name}

    # Family B: sorted/compartment-normalized -> gene OLS only, NO concordance (is_unbiased_single_cell=False)
    genes_ct, summary_ct, _ = run_family_b_genes_and_concordance(
        expr_only_ct,
        meta_only_ct,
        age_col="age_per_decade",
        candidate_covariates=("sex_male", "log_lib"),
        min_donors=MIN_DONORS,
        min_cpm=1.0,
        min_expr_frac=0.5,
        min_residual_df=2,
        is_unbiased_single_cell=False,
        seed=NULL_SEED,
    )
    genes_pool, summary_pool, _ = run_family_b_genes_and_concordance(
        {all_grp_name: expr_ct[all_grp_name]},
        {all_grp_name: meta_ct[all_grp_name]},
        age_col="age_per_decade",
        candidate_covariates=("sex_male", "log_lib"),
        min_donors=MIN_DONORS,
        min_cpm=1.0,
        min_expr_frac=0.5,
        min_residual_df=2,
        is_unbiased_single_cell=False,
        seed=NULL_SEED,
    )

    if not genes_ct.empty:
        genes_ct["q_bh_tissue_celltypes"] = multipletests(genes_ct["pvalue"], method="fdr_bh")[1]
        genes_ct["family_b_tissue_celltypes_size"] = int(len(genes_ct))
    else:
        genes_ct["q_bh_tissue_celltypes"] = np.nan
        genes_ct["family_b_tissue_celltypes_size"] = 0
    if not genes_pool.empty:
        genes_pool["q_bh_tissue_celltypes"] = np.nan
        genes_pool["family_b_tissue_celltypes_size"] = 0

    genes_all = pd.concat([genes_ct, genes_pool], ignore_index=True)
    genes_all.insert(0, "tissue", tissue)
    summary_all = pd.concat([summary_ct, summary_pool], ignore_index=True)
    summary_all.insert(0, "tissue", tissue)
    summary_all["is_pooled_tissue"] = summary_all["group"] == all_grp_name

    # Family C: P1_Clean 2,000 iterations
    fam_c = run_family_c_p1_signature(
        expr_ct,
        meta_ct,
        p1_program=p1_program,
        age_col="age_per_decade",
        candidate_covariates=("sex_male", "log_lib"),
        min_donors=MIN_DONORS,
        min_cpm=1.0,
        min_expr_frac=0.5,
        min_residual_df=2,
        iterations=P1_ITERATIONS,
        seed=NULL_SEED,
    )
    fam_c.insert(0, "tissue", tissue)
    fam_c = fam_c.rename(columns={"q_bh_composite": "q_bh_tissue_composite", "q_bh_empirical": "q_bh_tissue_empirical"})

    # Sensitivity: drop log_lib to preserve residual_df >= 3 for n=6 tissues
    fam_c_nolib = run_family_c_p1_signature(
        expr_ct,
        meta_ct,
        p1_program=p1_program,
        age_col="age_per_decade",
        candidate_covariates=("sex_male",),
        min_donors=MIN_DONORS,
        min_cpm=1.0,
        min_expr_frac=0.5,
        min_residual_df=2,
        iterations=P1_ITERATIONS,
        seed=NULL_SEED,
    )
    fam_c_nolib.insert(0, "tissue", tissue)
    fam_c_nolib.insert(0, "model", "sensitivity_age_sex_no_log_lib")

    return {
        "audit": audit,
        "donor_meta": donor_meta.reset_index(),
        "qc": pd.DataFrame(qc_rows),
        "genes": genes_all,
        "summary": summary_all,
        "fam_c": fam_c,
        "fam_c_nolib": fam_c_nolib,
        "all_grp_name": all_grp_name,
    }


def main() -> None:
    curated_donors = load_curated_donor_table()
    p1_program = load_p1_clean_human(REPO_ROOT)

    results_by_tissue = {}
    tissue_gates = {}
    for t in TARGET_TISSUES:
        print(f"Running TS2 tissue {t}...", flush=True)
        res = run_single_tissue(t, curated_donors, p1_program)
        results_by_tissue[t] = res
        tissue_gates[t] = res["audit"]

    all28 = build_all28_audit(tissue_gates)
    donors_df = pd.concat([results_by_tissue[t]["donor_meta"] for t in TARGET_TISSUES], ignore_index=True)
    qc_df = pd.concat([results_by_tissue[t]["qc"] for t in TARGET_TISSUES], ignore_index=True)
    genes_df = pd.concat([results_by_tissue[t]["genes"] for t in TARGET_TISSUES], ignore_index=True)
    summary_df = pd.concat([results_by_tissue[t]["summary"] for t in TARGET_TISSUES], ignore_index=True)
    fam_c_df = pd.concat([results_by_tissue[t]["fam_c"] for t in TARGET_TISSUES], ignore_index=True)
    fam_c_sens_df = pd.concat([results_by_tissue[t]["fam_c_nolib"] for t in TARGET_TISSUES], ignore_index=True)

    # Combined multi-tissue BH alongside per-tissue BH
    ct_mask = ~genes_df["group"].isin(["All_NonNK_Cells", "All_Cells"])
    genes_df["family_b_ts2_celltypes_size"] = int(ct_mask.sum())
    genes_df["q_bh_ts2_celltypes"] = np.nan
    genes_df.loc[ct_mask, "q_bh_ts2_celltypes"] = multipletests(genes_df.loc[ct_mask, "pvalue"], method="fdr_bh")[1]
    genes_df["family_b_ts2_all41_size"] = int(len(genes_df))
    genes_df["q_bh_ts2_all41"] = multipletests(genes_df["pvalue"], method="fdr_bh")[1]

    # Attach global hit counts into summary_df
    sig_tiss = genes_df[ct_mask & (genes_df["q_bh_tissue_celltypes"] < 0.05)].groupby(["tissue", "group"]).size()
    sig_ct36 = genes_df[ct_mask & (genes_df["q_bh_ts2_celltypes"] < 0.05)].groupby(["tissue", "group"]).size()
    sig_all41 = genes_df[genes_df["q_bh_ts2_all41"] < 0.05].groupby(["tissue", "group"]).size()
    summary_df["n_sig_tissue_celltypes_q05"] = [int(sig_tiss.get((r.tissue, r.group), 0)) for r in summary_df.itertuples(index=False)]
    summary_df["n_sig_ts2_celltypes36_q05"] = [int(sig_ct36.get((r.tissue, r.group), 0)) for r in summary_df.itertuples(index=False)]
    summary_df["n_sig_ts2_all41_q05"] = [int(sig_all41.get((r.tissue, r.group), 0)) for r in summary_df.itertuples(index=False)]

    # Family C cross-tissue BH (all 41)
    fam_c_df["family_c_ts2_all41_size"] = int(len(fam_c_df))
    fam_c_df["q_bh_ts2_all41_composite"] = multipletests(fam_c_df["composite_pvalue"], method="fdr_bh")[1]
    fam_c_df["q_bh_ts2_all41_empirical"] = multipletests(fam_c_df["empirical_null_pvalue"], method="fdr_bh")[1]

    fam_c_sens_df["family_c_ts2_all41_size"] = int(len(fam_c_sens_df))
    fam_c_sens_df["q_bh_ts2_all41_composite"] = multipletests(fam_c_sens_df["composite_pvalue"], method="fdr_bh")[1]
    fam_c_sens_df["q_bh_ts2_all41_empirical"] = multipletests(fam_c_sens_df["empirical_null_pvalue"], method="fdr_bh")[1]

    # Top genes per group (top 20 + any global/tissue FDR hits)
    sig_mask = (
        (genes_df["q_bh_within_group"] < 0.05)
        | (genes_df["q_bh_tissue_celltypes"] < 0.05)
        | (genes_df["q_bh_ts2_celltypes"] < 0.05)
        | (genes_df["q_bh_ts2_all41"] < 0.05)
    )
    top_n = (
        genes_df.sort_values("pvalue", ascending=True)
        .groupby(["tissue", "group"], sort=False)
        .head(20)
    )
    top_genes_df = pd.concat([genes_df[sig_mask], top_n], ignore_index=True).drop_duplicates(
        subset=["tissue", "group", "gene"]
    ).sort_values(["pvalue", "tissue", "group"]).reset_index(drop=True)

    write_clean_csv(all28, RESULTS_DIR / "followup-ts2-all28-tissues-audit.csv")
    write_clean_csv(donors_df, RESULTS_DIR / f"{PREFIX}-donors.csv")
    write_clean_csv(qc_df, RESULTS_DIR / f"{PREFIX}-qc.csv")
    write_clean_csv(summary_df, RESULTS_DIR / f"{PREFIX}-genes-summary.csv")
    write_clean_csv(top_genes_df, RESULTS_DIR / f"{PREFIX}-top-genes.csv")
    write_clean_csv(fam_c_df, RESULTS_DIR / f"{PREFIX}-p1-signature.csv")
    write_clean_csv(fam_c_sens_df, RESULTS_DIR / f"{PREFIX}-sensitivities-p1.csv")

    meta_summary = {
        "tissues_analyzed": TARGET_TISSUES,
        "unique_donors_across_5_tissues": sorted(donors_df["donor_id"].unique().tolist()),
        "n_unique_donors_across_5_tissues": int(donors_df["donor_id"].nunique()),
        "n_groups_total": int(len(qc_df)),
        "n_celltype_groups": int(ct_mask.loc[genes_df.drop_duplicates(["tissue","group"]).index].sum()) if False else int((~qc_df["is_pooled_tissue"]).sum()),
        "n_gene_tests_celltypes36": int(ct_mask.sum()),
        "n_gene_tests_all41": int(len(genes_df)),
        "n_sig_within_q05": int((genes_df["q_bh_within_group"] < 0.05).sum()),
        "n_sig_tissue_celltypes_q05": int((ct_mask & (genes_df["q_bh_tissue_celltypes"] < 0.05)).sum()),
        "n_sig_ts2_celltypes36_q05": int((ct_mask & (genes_df["q_bh_ts2_celltypes"] < 0.05)).sum()),
        "n_sig_ts2_all41_q05": int((genes_df["q_bh_ts2_all41"] < 0.05).sum()),
        "gates": {t: tissue_gates[t]["gate"] for t in TARGET_TISSUES},
    }
    (RESULTS_DIR / f"{PREFIX}-metadata-summary.json").write_text(json.dumps(meta_summary, indent=2), encoding="utf-8")
    print("DONE TS2 Remaining:", json.dumps(meta_summary, indent=2))


if __name__ == "__main__":
    main()
