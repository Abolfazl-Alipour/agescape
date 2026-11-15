"""Donor-level aging analysis of naked mole-rat splenic T cells (Lin et al. 2024, GSE214390)."""
from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.aging_programs import load_program_tiers
from src.followup_organs import (
    NULL_SEED,
    counts_to_log2cpm,
    run_family_b_genes_and_concordance,
    run_family_c_p1_signature,
)
from src.followup_muscle import count_gate

FAMILY_C_ITERATIONS = 2000
RAW_DIR = REPO_ROOT.parent / "agescape-revision-data" / "tmp_organs" / "nmr_gse214390"
RESULTS_DIR = REPO_ROOT / "results"
REPORTS_DIR = REPO_ROOT / "reports"

# Verified 1-to-1 sample map across GEO GSM6605086..092, Data S2 (MOESM5), and Data S21 (MOESM24).
NMR_SAMPLES: tuple[dict[str, Any], ...] = (
    {
        "sample_id": "nmr.2y_1",
        "gsm": "GSM6605086",
        "file_tag": "GSM6605086_nmr_spleen.T.cells_2y_2",
        "barcode_suffix": "nmr.2y_2",
    },
    {
        "sample_id": "nmr.2y_2",
        "gsm": "GSM6605087",
        "file_tag": "GSM6605087_nmr_spleen.T.cells_2y_3",
        "barcode_suffix": "nmr.2y_3",
    },
    {
        "sample_id": "nmr.2y_3",
        "gsm": "GSM6605088",
        "file_tag": "GSM6605088_nmr_spleen.T.cells_2y_4",
        "barcode_suffix": "nmr.2y_4",
    },
    {
        "sample_id": "nmr.24y_1",
        "gsm": "GSM6605089",
        "file_tag": "GSM6605089_nmr_spleen.T.cells_24y_1",
        "barcode_suffix": "nmr.24y_1",
    },
    {
        "sample_id": "nmr.24y_2",
        "gsm": "GSM6605090",
        "file_tag": "GSM6605090_nmr_spleen.T.cells_24y_2",
        "barcode_suffix": "nmr.24y_2",
    },
    {
        "sample_id": "nmr.24y_3",
        "gsm": "GSM6605091",
        "file_tag": "GSM6605091_nmr_spleen.T.cells_24y_3",
        "barcode_suffix": "nmr.24y_3",
    },
    {
        "sample_id": "nmr.24y_4",
        "gsm": "GSM6605092",
        "file_tag": "GSM6605092_nmr_spleen.T.cells_24y_4",
        "barcode_suffix": "nmr.24y_4",
    },
)

PRIMARY_SUBSETS: tuple[str, ...] = ("CD4.T.cell", "CD8.T.cell", "non.cytotoxic.gdT.cell")
EXCLUDED_SUBSETS: tuple[str, ...] = ("cytotoxic.gdT.cell", "ILC")
NMR_MT_CONTIG_IDS: frozenset[str] = frozenset(f"ENSHGLG000000000{i:02d}" for i in range(1, 39))
NMR_TECH_RE = re.compile(r"^(MT-|RP[LS]|MRP[LS]|ERCC-)")
NMR_SEX_SET: frozenset[str] = frozenset({"XIST", "TSIX", "DDX3Y", "EIF2S3Y", "UTY", "KDM5D"})
MIN_CELLS_PER_DONOR = 20


def load_gene_annotation() -> tuple[pd.DataFrame, dict[str, str], dict[str, dict[str, list[str]]]]:
    """Load 10x NMR features, symbol lookup, and mouse P1_Clean ortholog mappings."""
    genes_path = RAW_DIR / f"{NMR_SAMPLES[0]['file_tag']}_genes.tsv.gz"
    g_df = pd.read_csv(genes_path, sep="\t", header=None, names=["gene_id", "raw_symbol", "feature_type"], dtype=str)
    ens_df = pd.read_csv(
        REPO_ROOT / "data" / "reference" / "nmr-ensg2symbol.tsv",
        sep="\t",
        comment="#",
        header=None,
        names=["gene_id", "ensembl_symbol"],
        dtype=str,
    )
    ens_map = {
        str(gid).strip(): str(sym).strip()
        for gid, sym in zip(ens_df["gene_id"], ens_df["ensembl_symbol"].fillna(""))
        if str(sym).strip()
    }

    best_syms = []
    for gid, rsym in zip(g_df["gene_id"], g_df["raw_symbol"]):
        s = str(rsym).strip()
        best_syms.append(ens_map.get(gid, s) if (s == gid or not s) else s)
    g_df["best_symbol"] = best_syms
    g_df["is_technical"] = [
        (gid in NMR_MT_CONTIG_IDS)
        or bool(NMR_TECH_RE.match(sym.upper()))
        or (sym.upper() in NMR_SEX_SET)
        for gid, sym in zip(g_df["gene_id"], g_df["best_symbol"])
    ]

    ortho = pd.read_csv(REPO_ROOT / "data" / "reference" / "orthologs-mouse-nmr.csv", dtype=str)
    ortho_121 = ortho[ortho["orthology_type"] == "ortholog_one2one"].copy()
    nmr_to_mouse = (
        ortho_121.groupby("nmr_ensembl_id")["mouse_gene_symbol"]
        .apply(lambda s: ";".join(sorted(set(s.dropna().astype(str)))))
        .to_dict()
    )

    p1_mouse = load_program_tiers(REPO_ROOT, "mouse")["P1_Clean"]
    valid_ids = set(g_df["gene_id"])
    sym_10x_map = {str(sym).strip().upper(): str(gid) for gid, sym in zip(g_df["gene_id"], g_df["raw_symbol"])}
    sym_dups = set(g_df["raw_symbol"].str.upper().value_counts()[lambda s: s > 1].index)

    p1_primary: dict[str, list[str]] = {"up": [], "down": []}
    p1_extended: dict[str, list[str]] = {"up": [], "down": []}
    for direction in ("up", "down"):
        m_set = set(p1_mouse[direction])
        sub = ortho_121[ortho_121["mouse_gene_symbol"].isin(m_set) & ortho_121["nmr_ensembl_id"].isin(valid_ids)]
        mapped_ids = sorted(set(sub["nmr_ensembl_id"].astype(str)))
        p1_primary[direction] = mapped_ids
        ext_ids = set(mapped_ids)
        mapped_mouse = set(sub["mouse_gene_symbol"].astype(str))
        for mg in sorted(m_set - mapped_mouse):
            key = mg.upper()
            if key in sym_10x_map and key not in sym_dups:
                ext_ids.add(sym_10x_map[key])
        p1_extended[direction] = sorted(ext_ids)

    g_df["mouse_ortholog_121"] = g_df["gene_id"].map(nmr_to_mouse).fillna("")
    p1_tag_map: dict[str, str] = {}
    for gid in p1_primary["up"]:
        p1_tag_map[gid] = "P1_Clean_UP_ortholog121"
    for gid in p1_primary["down"]:
        p1_tag_map[gid] = "P1_Clean_DOWN_ortholog121"
    for gid in set(p1_extended["up"]) - set(p1_primary["up"]):
        p1_tag_map[gid] = "P1_Clean_UP_symbol_match"
    for gid in set(p1_extended["down"]) - set(p1_primary["down"]):
        p1_tag_map[gid] = "P1_Clean_DOWN_symbol_match"
    g_df["p1_clean_role"] = g_df["gene_id"].map(p1_tag_map).fillna("none")

    return g_df, {"primary": p1_primary, "extended": p1_extended}, nmr_to_mouse


def aggregate_nmr_samples(g_df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, pd.DataFrame], pd.DataFrame]:
    """Read raw 10x counts and aggregate donor pseudobulks stratified by Data S2 cell_subset."""
    df21 = pd.read_excel(RAW_DIR / "41467_2024_47264_MOESM24_ESM.xlsx", sheet_name="DataS21", header=3)
    s21_nmr = df21[(df21["species"] == "NMR") & (df21["tissue"] == "spleen.T.cells")].set_index("sample")
    df2 = pd.read_excel(RAW_DIR / "41467_2024_47264_MOESM5_ESM.xlsx", sheet_name="DataS2", header=3)
    df2["raw_barcode"] = df2["cell_barcode"].astype(str).str.split("_", n=1).str[0]
    df2["suffix"] = df2["cell_barcode"].astype(str).str.split("_", n=1).str[1]

    gene_ids = g_df["gene_id"].tolist()
    all_groups = list(PRIMARY_SUBSETS) + list(EXCLUDED_SUBSETS) + ["All_Splenic_T_Cells", "CD3plus_T_Cells_excl_ILC"]
    pb_counts: dict[str, dict[str, np.ndarray]] = {g: {} for g in all_groups}
    count_rows: list[dict[str, Any]] = []

    for spec in NMR_SAMPLES:
        sid = spec["sample_id"]
        tag = spec["file_tag"]
        with gzip.open(RAW_DIR / f"{tag}_barcodes.tsv.gz", "rt") as bf:
            bc_list = [line.strip() for line in bf]
        bc_to_col = {bc: i for i, bc in enumerate(bc_list)}

        with gzip.open(RAW_DIR / f"{tag}_matrix.mtx.gz", "rb") as mf:
            mat = scipy.io.mmread(mf)
        if not sp.issparse(mat):
            mat = sp.csc_matrix(mat)
        else:
            mat = mat.tocsc()
        if mat.shape != (len(gene_ids), len(bc_list)):
            raise ValueError(f"Shape mismatch in {sid}: {mat.shape} vs ({len(gene_ids)}, {len(bc_list)})")
        if not np.array_equal(mat.data, np.floor(mat.data)):
            raise ValueError(f"Non-integer counts in {sid}")

        sub = df2[df2["sample"] == sid].copy()
        if set(sub["suffix"]) != {spec["barcode_suffix"]}:
            raise ValueError(f"Suffix mismatch for {sid}: {set(sub['suffix'])}")
        col_idx = sub["raw_barcode"].map(bc_to_col)
        if col_idx.isna().any():
            raise ValueError(f"Unmatched barcodes in {sid}: {int(col_idx.isna().sum())}")
        sub["col_idx"] = col_idx.astype(int)

        meta_row = s21_nmr.loc[sid]
        rec: dict[str, Any] = {
            "sample_id": sid,
            "gsm": spec["gsm"],
            "s1_fig_label": str(meta_row["S1 Fig label"]),
            "colony": str(meta_row["colony.or.strain"]),
            "chip_id": str(meta_row["chip_ID"]),
            "sex": str(meta_row["sex"]).lower(),
            "age_months": float(meta_row["age.months"]),
            "age_years": float(meta_row["age.months"]) / 12.0,
            "age_group": "old_24y" if float(meta_row["age.months"]) > 100 else "young_2y",
            "old_vs_young": 1.0 if float(meta_row["age.months"]) > 100 else 0.0,
            "raw_10x_barcodes": int(len(bc_list)),
            "s21_n_cells": int(meta_row["n.cells"]),
            "qc_cells_total": int(len(sub)),
        }

        for ct in PRIMARY_SUBSETS + EXCLUDED_SUBSETS:
            idx = sub.loc[sub["cell_subset"] == ct, "col_idx"].to_numpy(dtype=int)
            rec[f"cells_{ct}"] = int(len(idx))
            if len(idx) > 0:
                pb_counts[ct][sid] = np.asarray(mat[:, idx].sum(axis=1)).ravel().astype(np.float64)

        idx_all = sub["col_idx"].to_numpy(dtype=int)
        pb_counts["All_Splenic_T_Cells"][sid] = np.asarray(mat[:, idx_all].sum(axis=1)).ravel().astype(np.float64)
        idx_cd3 = sub.loc[sub["cell_subset"] != "ILC", "col_idx"].to_numpy(dtype=int)
        rec["cells_CD3plus_excl_ILC"] = int(len(idx_cd3))
        pb_counts["CD3plus_T_Cells_excl_ILC"][sid] = np.asarray(mat[:, idx_cd3].sum(axis=1)).ravel().astype(np.float64)
        count_rows.append(rec)

    donor_df = pd.DataFrame(count_rows).set_index("sample_id")
    raw_pb_by_group = {
        g: pd.DataFrame(pb_counts[g], index=gene_ids)
        for g in all_groups
        if len(pb_counts[g]) > 0
    }
    return donor_df, raw_pb_by_group, df2


def build_group_inputs(
    g_df: pd.DataFrame,
    donor_df: pd.DataFrame,
    raw_pb_by_group: dict[str, pd.DataFrame],
    group_names: tuple[str, ...],
) -> tuple[dict[str, pd.DataFrame], dict[str, pd.DataFrame], list[dict[str, Any]]]:
    """Validate raw counts, mask NMR technical genes, normalize log2CPM, and attach log_lib."""
    keep_non_tech = ~g_df["is_technical"].to_numpy(dtype=bool)
    expr_by_group: dict[str, pd.DataFrame] = {}
    meta_by_group: dict[str, pd.DataFrame] = {}
    qc_records: list[dict[str, Any]] = []

    for g in group_names:
        cell_col = (
            "qc_cells_total"
            if g == "All_Splenic_T_Cells"
            else ("cells_CD3plus_excl_ILC" if g == "CD3plus_T_Cells_excl_ILC" else f"cells_{g}")
        )
        eligible = donor_df.index[donor_df[cell_col] >= MIN_CELLS_PER_DONOR].tolist()
        young_n = int((donor_df.loc[eligible, "old_vs_young"] == 0.0).sum())
        old_n = int((donor_df.loc[eligible, "old_vs_young"] == 1.0).sum())
        if len(eligible) < 4 or young_n < 2 or old_n < 2:
            qc_records.append(
                {
                    "group": g,
                    "eligible": False,
                    "skip_reason": f"low_cells_per_age_group(young={young_n},old={old_n})",
                    "n_donors": int(len(eligible)),
                    "young_n": young_n,
                    "old_n": old_n,
                    "total_cells": int(donor_df.loc[eligible, cell_col].sum()) if eligible else 0,
                }
            )
            continue

        raw_df = raw_pb_by_group[g][eligible]
        gate = count_gate(raw_df.to_numpy(dtype=np.float64))
        if gate["nonfinite"] > 0 or gate["negative"] > 0 or gate["noninteger"] > 0:
            raise ValueError(f"Invalid pseudobulk counts in {g}: {gate}")

        filt_df = raw_df.loc[keep_non_tech].copy()
        log2cpm, lib_sizes = counts_to_log2cpm(filt_df)
        m = donor_df.loc[eligible].copy()
        m["n_cells"] = m[cell_col].astype(int)
        m["lib_size_nontech"] = lib_sizes.loc[eligible].astype(float)
        m["log_lib"] = np.log10(np.maximum(m["lib_size_nontech"], 1.0))

        expr_by_group[g] = log2cpm
        meta_by_group[g] = m
        qc_records.append(
            {
                "group": g,
                "eligible": True,
                "skip_reason": "none",
                "n_donors": int(len(eligible)),
                "young_n": young_n,
                "old_n": old_n,
                "total_cells": int(m["n_cells"].sum()),
                "min_lib_size": float(lib_sizes.min()),
                "median_lib_size": float(lib_sizes.median()),
                "max_lib_size": float(lib_sizes.max()),
            }
        )
    return expr_by_group, meta_by_group, qc_records


def annotate_top_genes(genes_df: pd.DataFrame, g_df: pd.DataFrame, top_n: int = 50) -> pd.DataFrame:
    """Attach NMR symbol, mouse ortholog, P1_Clean membership and keep all global hits + top N per group."""
    ann = g_df.set_index("gene_id")[["raw_symbol", "best_symbol", "mouse_ortholog_121", "p1_clean_role"]]
    merged = genes_df.join(ann, on="gene")
    parts = []
    for _, sub in merged.groupby("group", sort=False):
        s = sub.sort_values(["pvalue", "gene"], ascending=[True, True])
        hits = s[s["q_bh_global"] < 0.05] if "q_bh_global" in s.columns else s.iloc[0:0]
        if len(hits) > top_n:
            parts.append(hits)
        else:
            parts.append(pd.concat([hits, s.drop(hits.index).head(top_n - len(hits))], axis=0))
    return pd.concat(parts, axis=0).reset_index(drop=True) if parts else merged


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    g_df, p1_maps, _ = load_gene_annotation()
    donor_df, raw_pb_by_group, _ = aggregate_nmr_samples(g_df)

    # Descriptive library-size & cell counts across all cells per donor (for audit table only).
    expr_primary, meta_primary, qc_primary = build_group_inputs(
        g_df, donor_df, raw_pb_by_group, PRIMARY_SUBSETS + EXCLUDED_SUBSETS
    )
    expr_pool, meta_pool, qc_pool = build_group_inputs(
        g_df, donor_df, raw_pb_by_group, ("All_Splenic_T_Cells", "CD3plus_T_Cells_excl_ILC")
    )

    donor_out = donor_df.copy()
    donor_out["all_cells_lib_size_nontech"] = meta_pool["All_Splenic_T_Cells"]["lib_size_nontech"]
    donor_out["all_cells_log_lib"] = meta_pool["All_Splenic_T_Cells"]["log_lib"]
    for ct in PRIMARY_SUBSETS + EXCLUDED_SUBSETS:
        donor_out[f"frac_{ct}"] = donor_out[f"cells_{ct}"] / donor_out["qc_cells_total"]
    donor_out.to_csv(RESULTS_DIR / "followup-nmr-spleen-t-gse214390-donors.csv")

    genes_prim, summary_prim, _ = run_family_b_genes_and_concordance(
        expr_primary,
        meta_primary,
        age_col="old_vs_young",
        candidate_covariates=("log_lib",),
        min_donors=4,
        min_residual_df=3,
        is_unbiased_single_cell=False,
    )
    summary_prim["analysis_tier"] = "primary_cell_subset"

    # Run pooled sensitivity separately so pooled tests do not double-count into primary 3-subset global BH.
    pool_summaries = []
    pool_gene_frames = []
    for pool_name in ("All_Splenic_T_Cells", "CD3plus_T_Cells_excl_ILC"):
        gp, sp_df, _ = run_family_b_genes_and_concordance(
            {pool_name: expr_pool[pool_name]},
            {pool_name: meta_pool[pool_name]},
            age_col="old_vs_young",
            candidate_covariates=("log_lib",),
            min_donors=4,
            min_residual_df=3,
            is_unbiased_single_cell=False,
        )
        sp_df["analysis_tier"] = "pooled_sensitivity"
        sp_df["n_sig_global_q05"] = np.nan
        sp_df["min_q_bh_global"] = np.nan
        sp_df["family_b_global_size"] = np.nan
        gp["q_bh_global"] = np.nan
        gp["family_b_global_size"] = np.nan
        pool_summaries.append(sp_df)
        pool_gene_frames.append(gp)

    # Summary CSV combining skipped subsets + tested subsets + pooled sensitivity.
    qc_map = {r["group"]: r for r in qc_primary + qc_pool}
    skip_rows = [
        {
            "group": r["group"],
            "status": f"skipped_{r['skip_reason']}",
            "n_donors": r["n_donors"],
            "n_genes_tested": 0,
            "analysis_tier": "excluded_sparse_subset",
        }
        for r in qc_primary
        if not r["eligible"]
    ]
    summary_all = pd.concat([summary_prim] + pool_summaries + [pd.DataFrame(skip_rows)], ignore_index=True)
    summary_all["total_cells"] = summary_all["group"].map(lambda g: qc_map.get(g, {}).get("total_cells", 0))
    summary_all["young_donors"] = summary_all["group"].map(lambda g: qc_map.get(g, {}).get("young_n", 0))
    summary_all["old_donors"] = summary_all["group"].map(lambda g: qc_map.get(g, {}).get("old_n", 0))
    summary_all["family_b_global_tests"] = int(len(genes_prim))
    summary_all.to_csv(RESULTS_DIR / "followup-nmr-spleen-t-gse214390-summary.csv", index=False)

    top_prim = annotate_top_genes(genes_prim, g_df, top_n=50)
    top_prim["analysis_tier"] = "primary_cell_subset"
    top_pools = []
    for gf in pool_gene_frames:
        tp = annotate_top_genes(gf, g_df, top_n=50)
        tp["analysis_tier"] = "pooled_sensitivity"
        tp["q_bh_global"] = np.nan
        tp["family_b_global_size"] = np.nan
        top_pools.append(tp)
    top_all = pd.concat([top_prim] + top_pools, ignore_index=True)
    top_all.to_csv(RESULTS_DIR / "followup-nmr-spleen-t-gse214390-top-genes.csv", index=False)

    # Save detected P1_Clean gene stats across primary & pooled subsets for biological audit.
    ann_idx = g_df.set_index("gene_id")[["raw_symbol", "best_symbol", "mouse_ortholog_121", "p1_clean_role"]]
    p1_genes_df = pd.concat([genes_prim] + pool_gene_frames, ignore_index=True).join(ann_idx, on="gene")
    p1_genes_df = p1_genes_df[p1_genes_df["p1_clean_role"] != "none"].sort_values(["group", "pvalue"])
    p1_genes_df.to_csv(RESULTS_DIR / "followup-nmr-spleen-t-gse214390-p1-genes.csv", index=False)

    # Family C: P1_Clean mouse signature across primary subsets and pooled sensitivity (both ortholog sets).
    p1_frames = []
    for variant, p1_prog in (("ortholog_one2one_105", p1_maps["primary"]), ("extended_symbol_113", p1_maps["extended"])):
        c_prim = run_family_c_p1_signature(
            expr_primary,
            meta_primary,
            p1_program=p1_prog,
            age_col="old_vs_young",
            candidate_covariates=("log_lib",),
            min_donors=4,
            min_residual_df=3,
            iterations=FAMILY_C_ITERATIONS,
            seed=NULL_SEED,
        )
        c_prim["mapping_variant"] = variant
        c_prim["analysis_tier"] = "primary_cell_subset"
        p1_frames.append(c_prim)

        for pool_name in ("All_Splenic_T_Cells", "CD3plus_T_Cells_excl_ILC"):
            c_pool = run_family_c_p1_signature(
                {pool_name: expr_pool[pool_name]},
                {pool_name: meta_pool[pool_name]},
                p1_program=p1_prog,
                age_col="old_vs_young",
                candidate_covariates=("log_lib",),
                min_donors=4,
                min_residual_df=3,
                iterations=FAMILY_C_ITERATIONS,
                seed=NULL_SEED,
            )
            c_pool["mapping_variant"] = variant
            c_pool["analysis_tier"] = "pooled_sensitivity"
            p1_frames.append(c_pool)

    p1_out = pd.concat(p1_frames, ignore_index=True)
    p1_out.to_csv(RESULTS_DIR / "followup-nmr-spleen-t-gse214390-p1-signature.csv", index=False)

    # Animal-outlier and colony sensitivity checks on P1_Clean (ortholog_one2one_105).
    sens_specs = [
        ("drop_nmr.24y_4_gdT_outlier", {"nmr.24y_4"}),
        ("drop_nmr.24y_1_ILC_outlier", {"nmr.24y_1"}),
        ("drop_nmr.2y_1_lowcell", {"nmr.2y_1"}),
        ("old_colony5000_only", {"nmr.24y_4"}),
    ]
    sens_frames = []
    combined_expr = {**expr_primary, **expr_pool}
    combined_meta = {**meta_primary, **meta_pool}
    for sname, drop_ids in sens_specs:
        sub_expr = {
            g: df[[c for c in df.columns if c not in drop_ids]]
            for g, df in combined_expr.items()
        }
        sub_meta = {
            g: m.loc[[idx for idx in m.index if idx not in drop_ids]]
            for g, m in combined_meta.items()
        }
        cs = run_family_c_p1_signature(
            sub_expr,
            sub_meta,
            p1_program=p1_maps["primary"],
            age_col="old_vs_young",
            candidate_covariates=("log_lib",),
            min_donors=4,
            min_residual_df=2,
            iterations=FAMILY_C_ITERATIONS,
            seed=NULL_SEED,
        )
        cs["sensitivity_check"] = sname
        sens_frames.append(cs)
    sens_out = pd.concat(sens_frames, ignore_index=True)
    sens_out.to_csv(RESULTS_DIR / "followup-nmr-spleen-t-gse214390-sensitivity.csv", index=False)

    print("=== Summary ===")
    print(summary_all.to_string(index=False))
    print("\n=== Family C P1_Clean (2000 perms) ===")
    print(
        p1_out[
            [
                "analysis_tier",
                "mapping_variant",
                "group",
                "n_donors",
                "n_p1_up_detected",
                "n_p1_down_detected",
                "composite_effect",
                "composite_pvalue",
                "empirical_null_pvalue",
                "q_bh_empirical",
            ]
        ].to_string(index=False)
    )
    print("\n=== Outlier Sensitivity (ortholog_one2one_105) ===")
    print(
        sens_out[
            [
                "sensitivity_check",
                "group",
                "n_donors",
                "residual_df",
                "composite_effect",
                "composite_pvalue",
                "empirical_null_pvalue",
                "q_bh_empirical",
            ]
        ].to_string(index=False)
    )
    print("\n=== Global BH < 0.05 hits (primary subsets) ===")
    prim_hits = genes_prim[genes_prim["q_bh_global"] < 0.05]
    print("Count:", len(prim_hits))
    if len(prim_hits) > 0:
        print(top_prim[top_prim["q_bh_global"] < 0.05].head(20)[["group", "gene", "best_symbol", "effect", "pvalue", "q_bh_global", "p1_clean_role"]].to_string(index=False))


if __name__ == "__main__":
    main()
