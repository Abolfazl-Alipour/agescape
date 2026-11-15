#!/usr/bin/env python3
"""Run GTEx v8 open-portal bulk RNA-seq age analyses (Track 1: all 49 >=50-donor SMTSD; Track 2: 6-tissue PyDESeq2 + 5 frozen sensitivities)."""

from __future__ import annotations

import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import argparse
import base64
from dataclasses import dataclass
import gc
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
from typing import Any

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import nb_nll
from scipy import stats as scipy_stats
from statsmodels.stats.multitest import multipletests

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.followup_muscle import (
    NULL_ITERATIONS,
    NULL_SEED,
    fit_ols,
    matched_null_effect_pvalue,
    random_signature_scores,
)
from src.followup_organs import (
    DATA_TMP_ROOT,
    counts_to_log2cpm,
    extended_technical_mask,
    load_p1_clean_human,
    run_family_b_genes_and_concordance,
    run_family_c_p1_signature,
    write_clean_csv,
)

GTEX_META_DIR = REPO_ROOT.parent / "agescape-revision-data" / "gtex-metadata"
MANIFEST_PATH = REPO_ROOT / "data" / "manifest" / "gtex_v8.json"
RESULTS_DIR = REPO_ROOT / "results"
EXPECTED_GCT_MD5_B64 = "HW7Jcxikd9f9grpyFDkl7g=="
EXPECTED_GCT_BYTES = 917_091_994

AGE_BIN_TO_ORDER = {
    "20-29": 1.0,
    "30-39": 2.0,
    "40-49": 3.0,
    "50-59": 4.0,
    "60-69": 5.0,
    "70-79": 6.0,
}

ISOLATED_CELL_LINES = {
    "Cells - Cultured fibroblasts",
    "Cells - EBV-transformed lymphocytes",
}

TRACK2_PRIMARY_TISSUES = [
    "Whole Blood",
    "Muscle - Skeletal",
    "Lung",
    "Kidney - Cortex",
    "Ovary",
    "Testis",
]

P1_ITERATIONS = 2000
MIN_DONORS_TRACK1 = 50
MIN_CPM = 1.0
MIN_EXPR_FRAC = 0.5
ENSEMBL_VER_RE = re.compile(r"\.\d+(?=(_PAR_Y)?$)")


@dataclass
class CovariateDesign:
    meta: pd.DataFrame
    numeric_cols: list[str]
    categorical_cols: list[str]
    all_cov_cols: list[str]
    cat_dummy_cols: list[str]
    n_samples_tissue: int
    n_samples_complete: int
    dropped_collinear: list[str]


def strip_ensembl_version_keep_par(gene_ids: np.ndarray | list[str]) -> list[str]:
    """Strip GENCODE .version while keeping _PAR_Y distinct so zero Ensembl ID collisions occur."""
    cleaned = [ENSEMBL_VER_RE.sub("", str(g).strip()) for g in gene_ids]
    if len(set(cleaned)) != len(cleaned):
        raise ValueError("Ensembl ID collision detected after version stripping")
    return cleaned


def verify_gct_file(gct_path: Path, expected_md5_b64: str = EXPECTED_GCT_MD5_B64) -> dict[str, Any]:
    """Verify GCT size and base64 MD5 checksum."""
    size = int(gct_path.stat().st_size)
    md5 = hashlib.md5()
    with gct_path.open("rb") as fh:
        while True:
            chunk = fh.read(8 * 1024 * 1024)
            if not chunk:
                break
            md5.update(chunk)
    actual_b64 = base64.b64encode(md5.digest()).decode("ascii")
    if actual_b64 != expected_md5_b64:
        raise ValueError(f"MD5 mismatch for {gct_path}: got {actual_b64}, expected {expected_md5_b64}")
    return {"path": str(gct_path), "size_bytes": size, "md5_base64": actual_b64, "md5_verified": True}


def load_gtex_sample_subject_metadata() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load GTEx v8 open-access subject phenotypes + RNASEQ sample attributes with 0 missing donor matches."""
    subj_path = GTEX_META_DIR / "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"
    samp_path = GTEX_META_DIR / "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
    subj = pd.read_csv(subj_path, sep="\t", dtype={"SUBJID": str, "AGE": str})
    samp = pd.read_csv(samp_path, sep="\t", low_memory=False)
    rnaseq = samp.loc[samp["SMAFRZE"].astype(str) == "RNASEQ"].copy()
    rnaseq["SUBJID"] = rnaseq["SAMPID"].astype(str).str.split("-").str[:2].str.join("-")
    merged = rnaseq.merge(subj, on="SUBJID", how="left", validate="many_to_one")
    if merged["AGE"].isna().any() or merged["SEX"].isna().any():
        raise ValueError("Unmatched RNASEQ sample to subject phenotype")
    merged["age_bin"] = merged["AGE"].astype(str).str.strip()
    merged["age_bin_order"] = merged["age_bin"].map(AGE_BIN_TO_ORDER)
    if merged["age_bin_order"].isna().any():
        raise ValueError(f"Unmapped AGE bins: {merged.loc[merged['age_bin_order'].isna(), 'AGE'].unique()}")
    merged["sex_male"] = (merged["SEX"].astype(int) == 1).astype(float)
    merged["SMTSISCH_hr"] = pd.to_numeric(merged["SMTSISCH"], errors="coerce") / 60.0
    merged["SMRIN"] = pd.to_numeric(merged["SMRIN"], errors="coerce")
    merged["SMEXNCRT"] = pd.to_numeric(merged["SMEXNCRT"], errors="coerce")
    merged["SMATSSCR"] = pd.to_numeric(merged["SMATSSCR"], errors="coerce")
    merged["DTHHRDY"] = pd.to_numeric(merged["DTHHRDY"], errors="coerce")
    merged["SMRDLGTH"] = pd.to_numeric(merged["SMRDLGTH"], errors="coerce")
    return merged.set_index("SAMPID", drop=False), subj


def stream_gct_to_memmap(
    gct_path: Path,
    valid_sample_ids: set[str],
    work_dir: Path,
) -> tuple[Path, tuple[int, int], pd.DataFrame, list[str], dict[str, Any]]:
    """Single-pass gzip -dc parse into int32 disk memmap + full non-negative integer count gate."""
    work_dir.mkdir(parents=True, exist_ok=True)
    mmap_path = work_dir / "counts_i32.mmap"
    gene_meta_path = work_dir / "genes_meta.csv"
    header_json_path = work_dir / "gct_header.json"
    if mmap_path.exists() and gene_meta_path.exists() and header_json_path.exists():
        info = json.loads(header_json_path.read_text(encoding="utf-8"))
        genes_df = pd.read_csv(gene_meta_path)
        return mmap_path, tuple(info["shape"]), genes_df, info["sample_ids"], info["gate"]

    proc = subprocess.Popen(
        ["gzip", "-dc", str(gct_path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        bufsize=16 * 1024 * 1024,
    )
    assert proc.stdout is not None
    try:
        v_line = proc.stdout.readline().decode("utf-8").strip()
        if not v_line.startswith("#1.2"):
            raise ValueError(f"Unexpected GCT version header: {v_line!r}")
        dims = proc.stdout.readline().decode("utf-8").strip().split("\t")
        n_genes_hdr, n_samples_hdr = int(dims[0]), int(dims[1])
        col_line = proc.stdout.readline().decode("utf-8").rstrip("\r\n").split("\t")
        if col_line[:2] != ["Name", "Description"]:
            raise ValueError(f"Unexpected GCT columns prefix: {col_line[:5]}")
        sample_ids = col_line[2:]
        if len(sample_ids) != n_samples_hdr:
            raise ValueError(f"Header sample count {n_samples_hdr} != {len(sample_ids)}")
        if set(sample_ids) != valid_sample_ids:
            raise ValueError(
                f"GCT sample set mismatch vs RNASEQ metadata: "
                f"missing={len(valid_sample_ids - set(sample_ids))}, extra={len(set(sample_ids) - valid_sample_ids)}"
            )

        mm = np.memmap(mmap_path, dtype=np.int32, mode="w+", shape=(n_genes_hdr, n_samples_hdr))
        raw_ids: list[str] = []
        descriptions: list[str] = []
        noninteger_lines = 0
        negative_values = 0
        total_nonzero = 0
        int32_max = np.iinfo(np.int32).max

        row_idx = 0
        for raw_line in proc.stdout:
            p1 = raw_line.find(b"\t")
            p2 = raw_line.find(b"\t", p1 + 1)
            gid = raw_line[:p1].decode("utf-8")
            sym = raw_line[p1 + 1 : p2].decode("utf-8")
            rest = raw_line[p2 + 1 :]
            if b"." in rest or b"e" in rest or b"E" in rest:
                noninteger_lines += 1
            arr = np.fromstring(rest, dtype=np.int64, sep="\t")
            if arr.size != n_samples_hdr:
                raise ValueError(f"Row {row_idx} ({gid}) parsed {arr.size} counts vs {n_samples_hdr}")
            neg_ct = int((arr < 0).sum())
            if neg_ct:
                negative_values += neg_ct
            if int(arr.max()) > int32_max:
                raise OverflowError(f"Count exceeds int32 max at row {row_idx} ({gid})")
            mm[row_idx, :] = arr.astype(np.int32, copy=False)
            total_nonzero += int((arr > 0).sum())
            raw_ids.append(gid)
            descriptions.append(sym)
            row_idx += 1
            if row_idx % 10000 == 0:
                print(f"[GCT stream] {row_idx}/{n_genes_hdr} genes parsed...", flush=True)

        if row_idx != n_genes_hdr:
            raise ValueError(f"Parsed {row_idx} gene rows vs header {n_genes_hdr}")
        mm.flush()
        del mm

        rc = proc.wait()
        stderr_txt = proc.stderr.read().decode("utf-8", errors="replace") if proc.stderr else ""
        if rc != 0:
            raise RuntimeError(f"gzip -dc failed rc={rc}: {stderr_txt}")
    finally:
        if proc.poll() is None:
            proc.kill()

    gate = {
        "n_genes": n_genes_hdr,
        "n_samples": n_samples_hdr,
        "nonfinite": 0,
        "negative": negative_values,
        "noninteger_lines": noninteger_lines,
        "total_nonzero_entries": total_nonzero,
    }
    if noninteger_lines or negative_values:
        raise ValueError(f"GCT failed count gate: {gate}")

    clean_ensembl = strip_ensembl_version_keep_par(raw_ids)
    genes_df = pd.DataFrame(
        {
            "Name": raw_ids,
            "ensembl_id": clean_ensembl,
            "Description": descriptions,
        }
    )
    genes_df.to_csv(gene_meta_path, index=False)
    header_json_path.write_text(
        json.dumps(
            {
                "shape": [n_genes_hdr, n_samples_hdr],
                "sample_ids": sample_ids,
                "gate": gate,
            }
        ),
        encoding="utf-8",
    )
    print(f"[GCT stream] Complete: shape=({n_genes_hdr}, {n_samples_hdr}), gate={gate}", flush=True)
    return mmap_path, (n_genes_hdr, n_samples_hdr), genes_df, sample_ids, gate


def build_clean_covariate_design(
    meta_sub: pd.DataFrame,
    numeric_candidates: list[str],
    categorical_candidates: list[str],
    age_col: str = "age_bin_order",
    min_level_count: int = 3,
    extra_age_cols: list[str] | None = None,
) -> CovariateDesign:
    """Filter complete cases and build full-rank numeric + categorical design columns."""
    req_cols = [age_col] + [c for c in numeric_candidates if c in meta_sub.columns] + [
        c for c in categorical_candidates if c in meta_sub.columns
    ]
    req_cols = list(dict.fromkeys(req_cols))
    n_tissue = int(len(meta_sub))
    work = meta_sub.dropna(subset=req_cols).copy()
    n_complete = int(len(work))

    kept_num: list[str] = []
    for c in numeric_candidates:
        if c not in work.columns:
            continue
        if work[c].nunique(dropna=True) >= 2 and float(work[c].std()) > 1e-12:
            work[c] = work[c].astype(float)
            kept_num.append(c)

    cat_dummy_cols: list[str] = []
    kept_cats: list[str] = []
    for c in categorical_candidates:
        if c not in work.columns:
            continue
        ser = work[c].astype(str)
        # Merge or drop rare levels? Only encode levels with >= min_level_count?
        # Wait: collapsing ultra-rare batch/category labels into 'OTHER' preserves samples cleanly while keeping rank stable;
        # if only 1 level remains with >= min_level_count, drop category cleanly.
        vc = ser.value_counts()
        major = vc[vc >= min_level_count].index.tolist()
        if len(major) < 2:
            continue
        # Ensure clean categorical names (most frequent level first as reference)
        ordered_levels = vc.index.tolist()
        ref = ordered_levels[0]
        kept_cats.append(c)
        work[c] = pd.Categorical(ser, categories=ordered_levels)
        for lvl in ordered_levels[1:]:
            if int(vc[lvl]) < min_level_count:
                # pool tiny (<3) isolated batch/site into reference/other? Better: only create dummy if >= min_level_count,
                # so rare levels don't absorb degrees of freedom; if needed for exactness, PyDESeq2 categorical handles >=2 levels cleanly.
                continue
            safe_lvl = re.sub(r"[^A-Za-z0-9_]+", "_", str(lvl)).strip("_")
            dcol = f"cov_{c}_{safe_lvl}"
            work[dcol] = (ser == lvl).astype(float)
            cat_dummy_cols.append(dcol)

    # Full-rank verification alongside age column(s) + intercept
    check_age_cols = [age_col] if extra_age_cols is None else list(extra_age_cols)
    base_cols = list(check_age_cols) + kept_num
    accepted_covs = list(kept_num)
    dropped_collinear: list[str] = []

    def _is_full_rank(cols: list[str]) -> bool:
        X = np.column_stack([np.ones(len(work), dtype=float), work[cols].to_numpy(dtype=float)])
        return np.linalg.matrix_rank(X) == X.shape[1]

    if not _is_full_rank(base_cols):
        # Drop any numeric column collinear with age/other covariates
        accepted_covs = []
        for c in kept_num:
            if _is_full_rank(check_age_cols + accepted_covs + [c]):
                accepted_covs.append(c)
            else:
                dropped_collinear.append(c)

    final_dummies: list[str] = []
    for dcol in cat_dummy_cols:
        if _is_full_rank(check_age_cols + accepted_covs + final_dummies + [dcol]):
            final_dummies.append(dcol)
        else:
            dropped_collinear.append(dcol)

    all_covs = accepted_covs + final_dummies
    return CovariateDesign(
        meta=work,
        numeric_cols=accepted_covs,
        categorical_cols=kept_cats,
        all_cov_cols=all_covs,
        cat_dummy_cols=final_dummies,
        n_samples_tissue=n_tissue,
        n_samples_complete=n_complete,
        dropped_collinear=dropped_collinear,
    )


def run_track1_all_tissues(
    mm: np.memmap,
    genes_df: pd.DataFrame,
    sample_id_to_col: dict[str, int],
    merged_meta: pd.DataFrame,
    p1_program: dict[str, list[str]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Track 1: all >=50-donor SMTSD tissues (47 primary separated from 2 isolated cell lines), OLS log2CPM + P1_Clean."""
    tech_mask = extended_technical_mask(genes_df["Description"].astype(str).to_numpy())
    nontech_idx = np.flatnonzero(~tech_mask)
    nontech_syms = genes_df["Description"].astype(str).to_numpy()[nontech_idx]
    # Group duplicate gene symbols once as CSR summation matrix -> avoids per-tissue pandas groupby overhead!
    uniq_syms, sym_inv = np.unique(nontech_syms, return_inverse=True)
    # Exclude empty/nan symbol names cleanly
    valid_sym = np.array([s != "" and s.lower() not in {"nan", "na"} for s in uniq_syms], dtype=bool)
    from scipy import sparse as sp
    sym_agg = sp.csr_matrix(
        (np.ones(len(nontech_idx), dtype=np.float64), (sym_inv, np.arange(len(nontech_idx)))),
        shape=(len(uniq_syms), len(nontech_idx)),
    )[valid_sym, :]
    clean_symbols = uniq_syms[valid_sym]
    sym_index = pd.Index(clean_symbols, dtype=str)

    tissue_vc = merged_meta["SMTSD"].value_counts()
    audit_rows: list[dict[str, Any]] = []
    primary_genes_list: list[pd.DataFrame] = []
    cellline_genes_list: list[pd.DataFrame] = []
    primary_p1_list: list[pd.DataFrame] = []
    cellline_p1_list: list[pd.DataFrame] = []
    tissue_summary_rows: list[dict[str, Any]] = []

    for smtsd in sorted(tissue_vc.index):
        t_meta = merged_meta.loc[merged_meta["SMTSD"] == smtsd].copy()
        if not t_meta["SUBJID"].is_unique:
            raise ValueError(f"Duplicate donor in {smtsd}")
        n_tot = int(len(t_meta))
        cohort_class = "excluded_lt50" if n_tot < MIN_DONORS_TRACK1 else (
            "isolated_cell_line" if smtsd in ISOLATED_CELL_LINES else "primary_tissue"
        )

        if n_tot < MIN_DONORS_TRACK1:
            audit_rows.append(
                {
                    "SMTSD": smtsd,
                    "cohort_class": cohort_class,
                    "n_samples": n_tot,
                    "n_complete": 0,
                    "n_male": int(t_meta["sex_male"].sum()),
                    "n_female": int((1 - t_meta["sex_male"]).sum()),
                    "age_bins_present": int(t_meta["age_bin"].nunique()),
                    "included_in_analysis": False,
                    "exclusion_reason": f"n_donors={n_tot} < {MIN_DONORS_TRACK1}",
                }
            )
            continue

        sids = list(t_meta.index.astype(str))
        col_indices = [sample_id_to_col[sid] for sid in sids]
        sub_raw = np.asarray(mm[nontech_idx, :][:, col_indices], dtype=np.float64)
        # Sum duplicate gene symbols -> (n_clean_symbols, n_samples)
        sym_counts = sym_agg @ sub_raw
        counts_df = pd.DataFrame(sym_counts, index=sym_index, columns=sids)
        log2cpm, lib_sizes = counts_to_log2cpm(counts_df)
        t_meta["lib_size"] = lib_sizes.loc[sids].to_numpy(dtype=float)
        t_meta["log_lib"] = np.log1p(t_meta["lib_size"].to_numpy(dtype=float))

        cov_design = build_clean_covariate_design(
            t_meta,
            numeric_candidates=["sex_male", "log_lib", "SMRIN", "SMTSISCH_hr", "SMEXNCRT"],
            categorical_candidates=["SMCENTER", "SMNABTCHT"],
            age_col="age_bin_order",
        )
        m_comp = cov_design.meta
        e_comp = log2cpm[list(m_comp.index)]

        audit_rows.append(
            {
                "SMTSD": smtsd,
                "cohort_class": cohort_class,
                "n_samples": n_tot,
                "n_complete": cov_design.n_samples_complete,
                "n_male": int(m_comp["sex_male"].sum()),
                "n_female": int((1 - m_comp["sex_male"]).sum()),
                "age_bins_present": int(m_comp["age_bin"].nunique()),
                "min_age_bin": sorted(m_comp["age_bin"].unique())[0],
                "max_age_bin": sorted(m_comp["age_bin"].unique())[-1],
                "covariates_used": ";".join(cov_design.all_cov_cols),
                "included_in_analysis": True,
                "exclusion_reason": "",
            }
        )

        # Bulk-only OLS (Family B + Family C, NO Family A, NO concordance)
        g_df, s_df, conc_df = run_family_b_genes_and_concordance(
            {smtsd: e_comp},
            {smtsd: m_comp},
            age_col="age_bin_order",
            candidate_covariates=tuple(cov_design.all_cov_cols),
            min_donors=MIN_DONORS_TRACK1,
            min_cpm=MIN_CPM,
            min_expr_frac=MIN_EXPR_FRAC,
            min_residual_df=10,
            is_unbiased_single_cell=False,
        )
        if not conc_df.empty:
            raise AssertionError(f"Bulk concordance must be empty for {smtsd}")
        p1_df = run_family_c_p1_signature(
            {smtsd: e_comp},
            {smtsd: m_comp},
            p1_program=p1_program,
            age_col="age_bin_order",
            candidate_covariates=tuple(cov_design.all_cov_cols),
            min_donors=MIN_DONORS_TRACK1,
            min_cpm=MIN_CPM,
            min_expr_frac=MIN_EXPR_FRAC,
            min_residual_df=10,
            iterations=P1_ITERATIONS,
            seed=NULL_SEED,
        )

        g_df["cohort_class"] = cohort_class
        p1_df["cohort_class"] = cohort_class
        if cohort_class == "primary_tissue":
            primary_genes_list.append(g_df)
            primary_p1_list.append(p1_df)
        else:
            cellline_genes_list.append(g_df)
            cellline_p1_list.append(p1_df)

        s_row = s_df.iloc[0].to_dict() if not s_df.empty else {"group": smtsd}
        s_row["cohort_class"] = cohort_class
        s_row["covariates_count"] = int(len(cov_design.all_cov_cols))
        s_row["p1_composite_effect"] = float(p1_df["composite_effect"].iloc[0]) if not p1_df.empty else np.nan
        s_row["p1_composite_t"] = (
            float(p1_df["composite_effect"].iloc[0] / p1_df["composite_se"].iloc[0])
            if (not p1_df.empty and float(p1_df["composite_se"].iloc[0]) > 0)
            else np.nan
        )
        s_row["p1_composite_p"] = float(p1_df["composite_pvalue"].iloc[0]) if not p1_df.empty else np.nan
        s_row["p1_emp_p"] = float(p1_df["empirical_null_pvalue"].iloc[0]) if not p1_df.empty else np.nan
        tissue_summary_rows.append(s_row)
        print(
            f"[Track1 {cohort_class}] {smtsd}: n={cov_design.n_samples_complete}/{n_tot}, "
            f"genes={s_row.get('n_genes_tested', 0)}, within_q<0.05={s_row.get('n_sig_within_q05', 0)}, "
            f"P1_eff={s_row['p1_composite_effect']:.4f} (emp_p={s_row['p1_emp_p']:.4g})",
            flush=True,
        )

    # Apply explicit separate BH corrections: primary_tissue (47 tissues) vs isolated_cell_line (2 lines)
    def _add_global_family_bh(glist: list[pd.DataFrame], p1list: list[pd.DataFrame], fname: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        gcat = pd.concat(glist, ignore_index=True)
        gcat["bh_family_name"] = fname
        gcat["bh_family_size"] = int(len(gcat))
        gcat["q_bh_family"] = multipletests(gcat["pvalue"], method="fdr_bh")[1]
        pcat = pd.concat(p1list, ignore_index=True)
        pcat["bh_family_name"] = fname
        pcat["family_c_size"] = int(len(pcat))
        pcat["q_bh_composite"] = multipletests(pcat["composite_pvalue"], method="fdr_bh")[1]
        pcat["q_bh_empirical"] = multipletests(pcat["empirical_null_pvalue"], method="fdr_bh")[1]
        return gcat, pcat

    g_prim, p1_prim = _add_global_family_bh(primary_genes_list, primary_p1_list, "track1_47_primary_tissues")
    g_cell, p1_cell = _add_global_family_bh(cellline_genes_list, cellline_p1_list, "track1_2_isolated_cell_lines")
    g_all = pd.concat([g_prim, g_cell], ignore_index=True)
    p1_all = pd.concat([p1_prim, p1_cell], ignore_index=True)

    # Attach global family q<0.05 counts back onto summary rows
    sig_fam = g_all.groupby("group")["q_bh_family"].apply(lambda s: int((s < 0.05).sum())).to_dict()
    sum_df = pd.DataFrame(tissue_summary_rows)
    sum_df["sig_q_family_0_05"] = sum_df["group"].map(sig_fam).fillna(0).astype(int)
    p1_qmap = p1_all.set_index("group")[["q_bh_composite", "q_bh_empirical"]].to_dict("index")
    sum_df["p1_q_bh_composite"] = sum_df["group"].map(lambda g: p1_qmap.get(g, {}).get("q_bh_composite", np.nan))
    sum_df["p1_q_bh_empirical"] = sum_df["group"].map(lambda g: p1_qmap.get(g, {}).get("q_bh_empirical", np.nan))

    top_genes = (
        g_all.sort_values(["pvalue", "group", "gene"], ascending=[True, True, True])
        .groupby("group", sort=False)
        .head(25)
        .reset_index(drop=True)
    )
    return pd.DataFrame(audit_rows), sum_df.sort_values(["cohort_class", "sig_q_family_0_05"], ascending=[True, False]), top_genes, p1_all


def fit_pydeseq2_continuous(
    counts_genes_by_samples: pd.DataFrame,
    meta_df: pd.DataFrame,
    cov_design: CovariateDesign,
    age_col: str = "age_bin_order",
    n_cpus: int = 12,
) -> pd.DataFrame:
    """Fit PyDESeq2 NB-GLM (~ age_bin_order + covariates) on filtered non-technical integer counts."""
    work_meta = cov_design.meta.copy()
    # Rename covariates cleanly for patsy formula if needed; use exact safe column names
    design_terms = [age_col] + cov_design.numeric_cols + cov_design.cat_dummy_cols
    # Confirm zero missing values & exact sample order
    work_meta = work_meta.loc[:, list(dict.fromkeys(design_terms))]
    # Center non-age numeric covariates? PyDESeq2 fits GLM with intercept; standardizing numeric covariates improves IRLS convergence cleanly without changing age slope or Wald p!
    for c in cov_design.numeric_cols:
        sd = float(work_meta[c].std())
        if sd > 1e-12:
            work_meta[c] = (work_meta[c] - float(work_meta[c].mean())) / sd
    formula = "~ " + " + ".join(design_terms)

    c_df = counts_genes_by_samples[list(work_meta.index)].T  # samples x genes
    dds = DeseqDataSet(
        counts=c_df,
        metadata=work_meta,
        design=formula,
        fit_type="parametric",
        refit_cooks=True,
        n_cpus=n_cpus,
        quiet=True,
    )
    dds.deseq2()
    dm_cols = list(dds.obsm["design_matrix"].columns)
    cvec = np.zeros(len(dm_cols), dtype=float)
    cvec[dm_cols.index(age_col)] = 1.0
    ds = DeseqStats(
        dds,
        contrast=cvec,
        cooks_filter=False,
        independent_filter=False,
        quiet=True,
        n_cpus=n_cpus,
    )
    ds.summary()
    res = ds.results_df.copy()
    res.index.name = "ensembl_id"
    res = res.reset_index()
    res["residual_df_approx"] = int(len(work_meta) - len(dm_cols))
    del dds, ds
    gc.collect()
    return res


def fit_pydeseq2_categorical_lrt(
    counts_genes_by_samples: pd.DataFrame,
    cov_design: CovariateDesign,
    n_cpus: int = 12,
) -> pd.DataFrame:
    """Fit PyDESeq2 categorical age_bin (~ age_cat + covariates) and compute joint NB-GLM LRT across bins."""
    work_meta = cov_design.meta.copy()
    bins_present = sorted(work_meta["age_bin"].astype(str).unique(), key=lambda b: AGE_BIN_TO_ORDER[b])
    work_meta["age_cat"] = pd.Categorical(work_meta["age_bin"].astype(str), categories=bins_present, ordered=False)
    for c in cov_design.numeric_cols:
        sd = float(work_meta[c].std())
        if sd > 1e-12:
            work_meta[c] = (work_meta[c] - float(work_meta[c].mean())) / sd

    cov_terms = cov_design.numeric_cols + cov_design.cat_dummy_cols
    full_formula = "~ age_cat" + ((" + " + " + ".join(cov_terms)) if cov_terms else "")
    c_df = counts_genes_by_samples[list(work_meta.index)].T

    dds = DeseqDataSet(
        counts=c_df,
        metadata=work_meta,
        design=full_formula,
        fit_type="parametric",
        refit_cooks=True,
        n_cpus=n_cpus,
        quiet=True,
    )
    dds.deseq2()
    dm_full = dds.obsm["design_matrix"]
    if np.linalg.matrix_rank(dm_full.to_numpy(dtype=float)) < dm_full.shape[1]:
        raise np.linalg.LinAlgError(f"Rank-deficient categorical design matrix: cols={list(dm_full.columns)}")
    age_dm_cols = [col for col in dm_full.columns if str(col).startswith("age_cat")]
    red_dm = dm_full.drop(columns=age_dm_cols)
    df_diff = int(len(age_dm_cols))

    nz_mask = dds.var["non_zero"].to_numpy(dtype=bool)
    nz_genes = dds.var_names[nz_mask]
    X_nz = np.asarray(dds.X[:, nz_mask], dtype=np.float64)
    disp_nz = dds.var.loc[nz_mask, "dispersions"].to_numpy(dtype=np.float64)
    sf = dds.obs["size_factors"].to_numpy(dtype=np.float64)
    # _mu_hat has shape (n_samples, n_genes_total) or (n_samples, n_nz)? Check exact width safely!
    mu_raw = dds.layers["_mu_hat"]
    mu_full_nz = np.asarray(mu_raw[:, nz_mask] if mu_raw.shape[1] == dds.n_vars else mu_raw, dtype=np.float64)

    _, mu_red_nz, _, _ = dds.inference.irls(
        counts=X_nz,
        size_factors=sf,
        design_matrix=red_dm.to_numpy(dtype=np.float64),
        disp=disp_nz,
        min_mu=dds.min_mu,
        beta_tol=dds.beta_tol,
    )
    ll_full = -nb_nll(X_nz, mu_full_nz, disp_nz)
    ll_red = -nb_nll(X_nz, mu_red_nz, disp_nz)
    lrt_stat = 2.0 * np.maximum(0.0, ll_full - ll_red)
    lrt_p = scipy_stats.chi2.sf(lrt_stat, df=df_diff)
    base_mean = np.asarray((X_nz / sf[:, None]).mean(axis=0), dtype=float)

    out = pd.DataFrame(
        {
            "ensembl_id": nz_genes.astype(str),
            "baseMean": base_mean,
            "lrt_stat": lrt_stat,
            "df_diff": df_diff,
            "pvalue": lrt_p,
        }
    )
    valid = np.isfinite(out["pvalue"].to_numpy(dtype=float))
    out["padj_within"] = np.nan
    if valid.any():
        out.loc[valid, "padj_within"] = multipletests(out.loc[valid, "pvalue"], method="fdr_bh")[1]
    del dds
    gc.collect()
    return out


def compute_p1_on_nb_results(
    nb_res: pd.DataFrame,
    gene_desc_map: dict[str, str],
    p1_program: dict[str, list[str]],
    iterations: int = P1_ITERATIONS,
) -> dict[str, Any]:
    """Evaluate P1_Clean UP-minus-DOWN effect, Stouffer z, and 2000-iter expression-matched null on NB-GLM results."""
    work = nb_res.copy()
    work["symbol"] = work["ensembl_id"].map(gene_desc_map).fillna("")
    work = work.loc[(work["symbol"] != "") & work["log2FoldChange"].notna() & work["lfcSE"].notna() & (work["lfcSE"] > 0)].copy()
    # Keep highest baseMean Ensembl ID per gene symbol so no duplicate symbol is double-counted
    work = work.sort_values("baseMean", ascending=False).drop_duplicates(subset=["symbol"], keep="first")
    z_all = (work["log2FoldChange"] / work["lfcSE"]).to_numpy(dtype=float)
    b_all = work["log2FoldChange"].to_numpy(dtype=float)
    mean_expr = np.log2(np.maximum(work["baseMean"].to_numpy(dtype=float), 0.0) + 1.0)
    sym_to_idx = {s: i for i, s in enumerate(work["symbol"].tolist())}

    up_idx = np.asarray([sym_to_idx[g] for g in p1_program["up"] if g in sym_to_idx], dtype=int)
    down_idx = np.asarray([sym_to_idx[g] for g in p1_program["down"] if g in sym_to_idx], dtype=int)
    eff = float(np.mean(b_all[up_idx]) - np.mean(b_all[down_idx]))
    z_signed = np.concatenate([z_all[up_idx], -z_all[down_idx]])
    comp_z = float(z_signed.sum() / np.sqrt(len(z_signed)))
    comp_p = float(2.0 * scipy_stats.norm.sf(abs(comp_z)))

    rng = np.random.default_rng(NULL_SEED)
    deciles = pd.Series(pd.qcut(mean_expr, 10, labels=False, duplicates="drop")).fillna(-1).to_numpy()
    available = np.ones(len(b_all), dtype=bool)
    available[np.concatenate([up_idx, down_idx])] = False
    up_dec_counts = pd.Series(deciles[up_idx]).value_counts().to_dict()
    down_dec_counts = pd.Series(deciles[down_idx]).value_counts().to_dict()
    pools_by_dec = {dec: np.flatnonzero((deciles == dec) & available) for dec in set(up_dec_counts) | set(down_dec_counts)}
    null_draws = np.empty(iterations, dtype=float)
    for it in range(iterations):
        r_up: list[np.ndarray] = []
        r_dn: list[np.ndarray] = []
        for dec, cnt in up_dec_counts.items():
            r_up.append(rng.choice(pools_by_dec[dec], size=int(cnt), replace=False))
        for dec, cnt in down_dec_counts.items():
            r_dn.append(rng.choice(pools_by_dec[dec], size=int(cnt), replace=False))
        null_draws[it] = float(b_all[np.concatenate(r_up)].mean() - b_all[np.concatenate(r_dn)].mean())
    emp_p = matched_null_effect_pvalue(eff, null_draws)
    return {
        "nb_p1_n_up": int(len(up_idx)),
        "nb_p1_n_down": int(len(down_idx)),
        "nb_p1_effect_lfc": eff,
        "nb_p1_composite_z": comp_z,
        "nb_p1_composite_p": comp_p,
        "nb_p1_empirical_p": emp_p,
    }


def run_track2_manifest_nb_glm(
    mm: np.memmap,
    genes_df: pd.DataFrame,
    sample_id_to_col: dict[str, int],
    merged_meta: pd.DataFrame,
    p1_program: dict[str, list[str]],
    n_cpus: int = 12,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Track 2: 6 primary tissues PyDESeq2 NB-GLM + all 5 frozen sensitivity arms from data/manifest/gtex_v8.json."""
    tech_mask = extended_technical_mask(genes_df["Description"].astype(str).to_numpy())
    nontech_mask = ~tech_mask
    nontech_idx = np.flatnonzero(nontech_mask)
    nontech_ens = genes_df["ensembl_id"].astype(str).to_numpy()[nontech_idx]
    ens_to_sym = dict(zip(genes_df["ensembl_id"].astype(str), genes_df["Description"].astype(str)))

    primary_gene_frames: list[pd.DataFrame] = []
    primary_summary_rows: list[dict[str, Any]] = []
    p1_nb_rows: list[dict[str, Any]] = []
    sens_rows: list[dict[str, Any]] = []

    def _prepare_counts_for_samples(sids: list[str]) -> pd.DataFrame:
        cols = [sample_id_to_col[sid] for sid in sids]
        mat = np.asarray(mm[nontech_idx, :][:, cols], dtype=np.int32)
        min_samp = max(5, int(np.ceil(0.20 * len(sids))))
        keep_g = (mat >= 10).sum(axis=1) >= min_samp
        return pd.DataFrame(mat[keep_g, :], index=pd.Index(nontech_ens[keep_g], dtype=str), columns=sids)

    def _summarize_wald(df: pd.DataFrame, tissue: str, arm: str, n_samp: int, rdf: int) -> dict[str, Any]:
        valid = df["pvalue"].notna().to_numpy()
        q_in = np.full(len(df), np.nan)
        if valid.any():
            q_in[valid] = multipletests(df.loc[valid, "pvalue"], method="fdr_bh")[1]
        df["q_bh_within_tissue"] = q_in
        sig = df.loc[df["q_bh_within_tissue"] < 0.05]
        return {
            "tissue": tissue,
            "arm": arm,
            "n_donors": n_samp,
            "residual_df_approx": rdf,
            "genes_tested": int(valid.sum()),
            "sig_q05_within": int(len(sig)),
            "n_up_q05": int((sig["log2FoldChange"] > 0).sum()),
            "n_down_q05": int((sig["log2FoldChange"] < 0).sum()),
            "min_pvalue": float(df["pvalue"].min()),
        }

    prim_by_tissue: dict[str, pd.DataFrame] = {}

    for tissue in TRACK2_PRIMARY_TISSUES:
        t_meta = merged_meta.loc[merged_meta["SMTSD"] == tissue].copy()
        num_cands = ["sex_male", "SMRIN", "SMTSISCH_hr", "SMEXNCRT"]
        cat_cands = ["SMCENTER", "SMNABTCHT"]
        cov_d = build_clean_covariate_design(t_meta, num_cands, cat_cands, age_col="age_bin_order")
        c_mat = _prepare_counts_for_samples(list(cov_d.meta.index))

        print(f"[Track2 Primary] Fitting PyDESeq2 {tissue}: n={cov_d.n_samples_complete}, genes={len(c_mat)}", flush=True)
        res_p = fit_pydeseq2_continuous(c_mat, cov_d.meta, cov_d, age_col="age_bin_order", n_cpus=n_cpus)
        res_p.insert(0, "tissue", tissue)
        res_p["symbol"] = res_p["ensembl_id"].map(ens_to_sym)
        s_row = _summarize_wald(
            res_p, tissue, "primary_continuous", cov_d.n_samples_complete, int(res_p["residual_df_approx"].iloc[0])
        )
        p1_nb = compute_p1_on_nb_results(res_p, ens_to_sym, p1_program, iterations=P1_ITERATIONS)
        p1_nb["tissue"] = tissue
        p1_nb_rows.append(p1_nb)
        primary_gene_frames.append(res_p)
        primary_summary_rows.append({**s_row, **p1_nb})
        prim_by_tissue[tissue] = res_p

        # Sensitivity 1: categorical_age LRT
        n_bins = int(cov_d.meta["age_bin"].nunique())
        c_cat = _prepare_counts_for_samples(list(cov_d.meta.index))
        print(f"[Track2 Sens 1] categorical_age LRT {tissue}: bins={n_bins}", flush=True)
        lrt_df = fit_pydeseq2_categorical_lrt(c_cat, cov_d, n_cpus=n_cpus)
        m_lrt = res_p.merge(lrt_df[["ensembl_id", "lrt_stat", "df_diff", "pvalue", "padj_within"]], on="ensembl_id", suffixes=("_prim", "_lrt"))
        both_sig = int(((m_lrt["q_bh_within_tissue"] < 0.05) & (m_lrt["padj_within"] < 0.05)).sum())
        sens_rows.append(
            {
                "tissue": tissue,
                "sensitivity": "categorical_age_lrt",
                "n_donors": cov_d.n_samples_complete,
                "genes_tested": int(lrt_df["pvalue"].notna().sum()),
                "sig_q05_within": int((lrt_df["padj_within"] < 0.05).sum()),
                "overlap_with_primary_q05": both_sig,
                "spearman_stat_vs_primary": float(
                    scipy_stats.spearmanr(m_lrt["stat"].abs(), m_lrt["lrt_stat"], nan_policy="omit").statistic
                ),
                "notes": f"df_diff={int(lrt_df['df_diff'].iloc[0])}",
            }
        )

        # Sensitivity 2: oldest_bin_tail (exclude 70-79)
        t_no70 = t_meta.loc[t_meta["age_bin"] != "70-79"].copy()
        cov_no70 = build_clean_covariate_design(t_no70, num_cands, cat_cands, age_col="age_bin_order")
        c_no70 = _prepare_counts_for_samples(list(cov_no70.meta.index))
        print(f"[Track2 Sens 2] oldest_bin_tail {tissue}: n={cov_no70.n_samples_complete}", flush=True)
        res_no70 = fit_pydeseq2_continuous(c_no70, cov_no70.meta, cov_no70, age_col="age_bin_order", n_cpus=n_cpus)
        s_no70 = _summarize_wald(res_no70, tissue, "oldest_bin_tail_excl_70_79", cov_no70.n_samples_complete, int(res_no70["residual_df_approx"].iloc[0]))
        m_no70 = res_p.merge(res_no70[["ensembl_id", "log2FoldChange", "stat", "q_bh_within_tissue"]], on="ensembl_id", suffixes=("_prim", "_sens"))
        sens_rows.append(
            {
                "tissue": tissue,
                "sensitivity": "oldest_bin_tail_excl_70_79",
                "n_donors": s_no70["n_donors"],
                "genes_tested": s_no70["genes_tested"],
                "sig_q05_within": s_no70["sig_q05_within"],
                "overlap_with_primary_q05": int(((m_no70["q_bh_within_tissue_prim"] < 0.05) & (m_no70["q_bh_within_tissue_sens"] < 0.05)).sum()),
                "spearman_stat_vs_primary": float(scipy_stats.spearmanr(m_no70["log2FoldChange_prim"], m_no70["log2FoldChange_sens"], nan_policy="omit").statistic),
                "notes": f"dropped {cov_d.n_samples_complete - s_no70['n_donors']} donors in 70-79",
            }
        )

        # Sensitivity 3: extended_sample_quality_and_death (+ SMATSSCR if non-blood, + DTHHRDY)
        ext_cats = list(cat_cands) + (["DTHHRDY"] if tissue == "Whole Blood" else ["SMATSSCR", "DTHHRDY"])
        cov_ext = build_clean_covariate_design(t_meta, num_cands, ext_cats, age_col="age_bin_order")
        c_ext = _prepare_counts_for_samples(list(cov_ext.meta.index))
        print(f"[Track2 Sens 3] extended_quality_death {tissue}: n={cov_ext.n_samples_complete}, covs={len(cov_ext.all_cov_cols)}", flush=True)
        res_ext = fit_pydeseq2_continuous(c_ext, cov_ext.meta, cov_ext, age_col="age_bin_order", n_cpus=n_cpus)
        s_ext = _summarize_wald(res_ext, tissue, "extended_quality_and_death", cov_ext.n_samples_complete, int(res_ext["residual_df_approx"].iloc[0]))
        m_ext = res_p.merge(res_ext[["ensembl_id", "log2FoldChange", "q_bh_within_tissue"]], on="ensembl_id", suffixes=("_prim", "_sens"))
        sens_rows.append(
            {
                "tissue": tissue,
                "sensitivity": "extended_quality_and_death",
                "n_donors": s_ext["n_donors"],
                "genes_tested": s_ext["genes_tested"],
                "sig_q05_within": s_ext["sig_q05_within"],
                "overlap_with_primary_q05": int(((m_ext["q_bh_within_tissue_prim"] < 0.05) & (m_ext["q_bh_within_tissue_sens"] < 0.05)).sum()),
                "spearman_stat_vs_primary": float(scipy_stats.spearmanr(m_ext["log2FoldChange_prim"], m_ext["log2FoldChange_sens"], nan_policy="omit").statistic),
                "notes": f"added {','.join(c for c in ext_cats if c not in cat_cands)}",
            }
        )

        # Sensitivity 4: blood_signed_ischemic_time (Whole Blood SMTSISCH >= 0)
        if int((t_meta["SMTSISCH"] < 0).sum()) > 0:
            t_pos = t_meta.loc[t_meta["SMTSISCH"] >= 0].copy()
            cov_pos = build_clean_covariate_design(t_pos, num_cands, cat_cands, age_col="age_bin_order")
            c_pos = _prepare_counts_for_samples(list(cov_pos.meta.index))
            print(f"[Track2 Sens 4] blood_signed_ischemic_time {tissue}: n={cov_pos.n_samples_complete}", flush=True)
            res_pos = fit_pydeseq2_continuous(c_pos, cov_pos.meta, cov_pos, age_col="age_bin_order", n_cpus=n_cpus)
            s_pos = _summarize_wald(res_pos, tissue, "blood_signed_ischemic_time", cov_pos.n_samples_complete, int(res_pos["residual_df_approx"].iloc[0]))
            m_pos = res_p.merge(res_pos[["ensembl_id", "log2FoldChange", "q_bh_within_tissue"]], on="ensembl_id", suffixes=("_prim", "_sens"))
            sens_rows.append(
                {
                    "tissue": tissue,
                    "sensitivity": "blood_signed_ischemic_time",
                    "n_donors": s_pos["n_donors"],
                    "genes_tested": s_pos["genes_tested"],
                    "sig_q05_within": s_pos["sig_q05_within"],
                    "overlap_with_primary_q05": int(((m_pos["q_bh_within_tissue_prim"] < 0.05) & (m_pos["q_bh_within_tissue_sens"] < 0.05)).sum()),
                    "spearman_stat_vs_primary": float(scipy_stats.spearmanr(m_pos["log2FoldChange_prim"], m_pos["log2FoldChange_sens"], nan_policy="omit").statistic),
                    "notes": "restricted to SMTSISCH >= 0",
                }
            )

        # Sensitivity 5: read_length_76bp (only refit if non-76bp samples present; otherwise note identical to primary)
        non76 = int((cov_d.meta["SMRDLGTH"] != 76).sum())
        if non76 > 0:
            t_76 = t_meta.loc[t_meta["SMRDLGTH"] == 76].copy()
            cov_76 = build_clean_covariate_design(t_76, num_cands, cat_cands, age_col="age_bin_order")
            c_76 = _prepare_counts_for_samples(list(cov_76.meta.index))
            print(f"[Track2 Sens 5] read_length_76bp {tissue}: n={cov_76.n_samples_complete} (removed {non76})", flush=True)
            res_76 = fit_pydeseq2_continuous(c_76, cov_76.meta, cov_76, age_col="age_bin_order", n_cpus=n_cpus)
            s_76 = _summarize_wald(res_76, tissue, "read_length_76bp", cov_76.n_samples_complete, int(res_76["residual_df_approx"].iloc[0]))
            m_76 = res_p.merge(res_76[["ensembl_id", "log2FoldChange", "q_bh_within_tissue"]], on="ensembl_id", suffixes=("_prim", "_sens"))
            sens_rows.append(
                {
                    "tissue": tissue,
                    "sensitivity": "read_length_76bp",
                    "n_donors": s_76["n_donors"],
                    "genes_tested": s_76["genes_tested"],
                    "sig_q05_within": s_76["sig_q05_within"],
                    "overlap_with_primary_q05": int(((m_76["q_bh_within_tissue_prim"] < 0.05) & (m_76["q_bh_within_tissue_sens"] < 0.05)).sum()),
                    "spearman_stat_vs_primary": float(scipy_stats.spearmanr(m_76["log2FoldChange_prim"], m_76["log2FoldChange_sens"], nan_policy="omit").statistic),
                    "notes": f"dropped {non76} non-76bp (101bp) samples",
                }
            )
        else:
            sens_rows.append(
                {
                    "tissue": tissue,
                    "sensitivity": "read_length_76bp",
                    "n_donors": cov_d.n_samples_complete,
                    "genes_tested": s_row["genes_tested"],
                    "sig_q05_within": s_row["sig_q05_within"],
                    "overlap_with_primary_q05": s_row["sig_q05_within"],
                    "spearman_stat_vs_primary": 1.0,
                    "notes": "100% 76bp already (identical to primary)",
                }
            )

    all_prim = pd.concat(primary_gene_frames, ignore_index=True)
    valid_g = all_prim["pvalue"].notna().to_numpy()
    all_prim["family_b_track2_size"] = int(valid_g.sum())
    all_prim["q_bh_track2_global"] = np.nan
    all_prim.loc[valid_g, "q_bh_track2_global"] = multipletests(all_prim.loc[valid_g, "pvalue"], method="fdr_bh")[1]

    sig_glob = all_prim.groupby("tissue")["q_bh_track2_global"].apply(lambda s: int((s < 0.05).sum())).to_dict()
    sum_t2_df = pd.DataFrame(primary_summary_rows)
    sum_t2_df["sig_q05_global_6tissues"] = sum_t2_df["tissue"].map(sig_glob).fillna(0).astype(int)

    p1_nb_df = pd.DataFrame(p1_nb_rows)
    p1_nb_df["family_size"] = int(len(p1_nb_df))
    p1_nb_df["q_bh_composite"] = multipletests(p1_nb_df["nb_p1_composite_p"], method="fdr_bh")[1]
    p1_nb_df["q_bh_empirical"] = multipletests(p1_nb_df["nb_p1_empirical_p"], method="fdr_bh")[1]

    top_t2 = (
        all_prim.sort_values(["pvalue", "tissue", "ensembl_id"], ascending=[True, True, True])
        .groupby("tissue", sort=False)
        .head(40)
        .reset_index(drop=True)
    )
    return sum_t2_df, top_t2, pd.DataFrame(sens_rows), p1_nb_df


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gct", required=True, type=Path)
    parser.add_argument("--work-dir", default=DATA_TMP_ROOT / "gtex_v8", type=Path)
    parser.add_argument("--n-cpus", default=12, type=int)
    args = parser.parse_args()

    file_audit = verify_gct_file(args.gct)
    file_audit["path"] = Path(file_audit["path"]).name
    print(f"[GTEx v8] Verified GCT: {file_audit}", flush=True)

    merged_meta, _ = load_gtex_sample_subject_metadata()
    mmap_path, shape, genes_df, sample_ids, gate = stream_gct_to_memmap(
        args.gct, set(merged_meta.index.astype(str)), args.work_dir
    )
    sample_id_to_col = {sid: idx for idx, sid in enumerate(sample_ids)}
    mm = np.memmap(mmap_path, dtype=np.int32, mode="r", shape=shape)
    p1_program = load_p1_clean_human(REPO_ROOT)

    print("=== GTEx Track 1: All 49 >=50-donor SMTSD tissues ===", flush=True)
    t1_audit, t1_sum, t1_top, t1_p1 = run_track1_all_tissues(
        mm, genes_df, sample_id_to_col, merged_meta, p1_program
    )
    write_clean_csv(t1_audit, RESULTS_DIR / "followup-gtex-v8-tissue-audit.csv")
    write_clean_csv(t1_sum, RESULTS_DIR / "followup-gtex-v8-track1-genes-summary.csv")
    write_clean_csv(t1_top, RESULTS_DIR / "followup-gtex-v8-track1-top-genes.csv")
    write_clean_csv(t1_p1, RESULTS_DIR / "followup-gtex-v8-track1-p1-signature.csv")

    print("=== GTEx Track 2: 6-tissue PyDESeq2 NB-GLM + 5 sensitivities ===", flush=True)
    t2_sum, t2_top, t2_sens, t2_p1 = run_track2_manifest_nb_glm(
        mm, genes_df, sample_id_to_col, merged_meta, p1_program, n_cpus=args.n_cpus
    )
    write_clean_csv(t2_sum, RESULTS_DIR / "followup-gtex-v8-track2-nbglm-summary.csv")
    write_clean_csv(t2_top, RESULTS_DIR / "followup-gtex-v8-track2-nbglm-top-genes.csv")
    write_clean_csv(t2_sens, RESULTS_DIR / "followup-gtex-v8-track2-sensitivities.csv")
    write_clean_csv(t2_p1, RESULTS_DIR / "followup-gtex-v8-track2-p1-nbglm.csv")

    summary_json = {
        "file_audit": file_audit,
        "count_gate": gate,
        "track1": {
            "primary_tissues": int((t1_audit["cohort_class"] == "primary_tissue").sum()),
            "isolated_cell_lines": int((t1_audit["cohort_class"] == "isolated_cell_line").sum()),
            "excluded_lt50": int((t1_audit["cohort_class"] == "excluded_lt50").sum()),
            "primary_p1_emp_q05": int((t1_p1.loc[t1_p1["cohort_class"] == "primary_tissue", "q_bh_empirical"] < 0.05).sum()),
            "cellline_p1_emp_q05": int((t1_p1.loc[t1_p1["cohort_class"] == "isolated_cell_line", "q_bh_empirical"] < 0.05).sum()),
        },
        "track2": {
            "tissues": TRACK2_PRIMARY_TISSUES,
            "total_nb_tests_global": int(t2_sum["genes_tested"].sum()),
            "total_sig_q05_global": int(t2_sum["sig_q05_global_6tissues"].sum()),
            "nb_p1_emp_q05": int((t2_p1["q_bh_empirical"] < 0.05).sum()),
        },
    }
    (RESULTS_DIR / "followup-gtex-v8-summary.json").write_text(
        json.dumps(summary_json, indent=2), encoding="utf-8"
    )
    print("DONE GTEx v8:", json.dumps(summary_json, indent=2), flush=True)


if __name__ == "__main__":
    main()
