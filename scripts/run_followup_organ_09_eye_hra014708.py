"""Dataset 09: Retinal aging cohort (Ophthalmology Science 2026, PMC12892068, PMID 41684508).

Primary source verification:
- Paper: Ophthalmology Science (2026), DOI 10.1016/j.xops.2025.101062, PMID 41684508, PMCID PMC12892068.
- Data Availability cites two accessions:
  1) GSA-Human HRA014708 (controlled access at CNCB-NGDC; HTTP browse page returns anti-bot
     SecurityCheck HTML interstitial; not bypassed and no DAC application submitted).
  2) Figshare doi:10.6084/m9.figshare.28757048 (public article with 18 Cell Ranger v8.0.1
     10x count matrix triplets for sample_1 to sample_18 plus README.docx and Script.py,
     omitting the author-clustered /RPE.h5ad object).
- Table 1 of PMC12892068 maps the 18 eye samples (Sample 1..18) to 12 unique donors (D1..D12):
  6 bilateral cadaver donors (D7..D12, L+R eyes collapsed to donor unit) and 6 unilateral donors (D1..D6).
  8 samples (from 4 donors D9, D10, D11, D12) have Diabetic Retinopathy (DR); 3 samples (2 donors D1, D8)
  have DM without DR; 7 samples (6 donors D2..D7) are Non-DM.
- Streams Figshare 10x matrices one sample at a time (<1.5 GB RAM, <300 MB disk), verifies sex-marker
  concordance against Table 1, applies the paper's exact QC (>=200 genes, >500 UMI, <15% mito),
  aggregates raw integer counts to unique donors x cell classes (via Figure 1C canonical markers)
  plus exact unpartitioned whole-retina donor pseudobulk, and runs disease-adjusted and sensitivity models.
"""

from __future__ import annotations

import gc
import gzip
import io
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.followup_organs import (
    DATA_TMP_ROOT,
    clean_temp_dir,
    counts_to_log2cpm,
    finalize_pseudobulk_group,
    load_p1_clean_human,
    run_family_a_clr,
    run_family_b_genes_and_concordance,
    run_family_c_p1_signature,
    write_clean_csv,
    write_clean_markdown,
)

NULL_ITERATIONS = 200
NULL_SEED = 202609

# Verified verbatim from PMC12892068 Table 1
TABLE1_SAMPLES = pd.DataFrame(
    [
        (1, "sample_1", "D1", "B01", "Live", "L", "M", 76.0, "DM", "Corneal staphyloma", "L.H.", "Old"),
        (2, "sample_2", "D2", "B02", "Live", "L", "M", 55.0, "Non-DM", "Phthisis bulbi", "L.H.", "Young"),
        (3, "sample_3", "D3", "B03", "Live", "R", "M", 71.0, "Non-DM", "Corneal ulcer", "Y.H.", "Old"),
        (4, "sample_4", "D4", "B04", "Live", "L", "F", 50.0, "Non-DM", "Phthisis bulbi", "Y.H.", "Young"),
        (5, "sample_5", "D5", "B04", "Live", "R", "F", 68.0, "Non-DM", "Corneal perforation", "Y.H.", "Old"),
        (6, "sample_6", "D6", "B05", "Cadaver", "L", "M", 34.0, "Non-DM", "Cerebral hemorrhage", "Y.H.", "Young"),
        (7, "sample_7", "D7", "B06", "Cadaver", "L", "F", 92.0, "Non-DM", "Natural death", "Y.H.", "Old"),
        (8, "sample_8", "D7", "B06", "Cadaver", "R", "F", 92.0, "Non-DM", "Natural death", "Y.H.", "Old"),
        (9, "sample_9", "D8", "B07", "Cadaver", "L", "F", 55.0, "DM", "Traumatic brain injury", "Y.H.", "Young"),
        (10, "sample_10", "D8", "B07", "Cadaver", "R", "F", 55.0, "DM", "Traumatic brain injury", "Y.H.", "Young"),
        (11, "sample_11", "D9", "B08", "Cadaver", "L", "M", 55.0, "DR", "Brainstem hemorrhage", "Y.H.", "Young"),
        (12, "sample_12", "D9", "B08", "Cadaver", "R", "M", 55.0, "DR", "Brainstem hemorrhage", "Y.H.", "Young"),
        (13, "sample_13", "D10", "B09", "Cadaver", "L", "M", 35.0, "DR", "Traumatic shock", "Y.H.", "Young"),
        (14, "sample_14", "D10", "B09", "Cadaver", "R", "M", 35.0, "DR", "Traumatic shock", "Y.H.", "Young"),
        (15, "sample_15", "D11", "B10", "Cadaver", "L", "M", 84.0, "DR", "Prostate cancer", "Y.H.", "Old"),
        (16, "sample_16", "D11", "B10", "Cadaver", "R", "M", 84.0, "DR", "Prostate cancer", "Y.H.", "Old"),
        (17, "sample_17", "D12", "B11", "Cadaver", "L", "M", 69.0, "DR", "Lung cancer", "Y.H.", "Old"),
        (18, "sample_18", "D12", "B11", "Cadaver", "R", "M", 69.0, "DR", "Lung cancer", "Y.H.", "Old"),
    ],
    columns=[
        "sample_num",
        "figshare_folder",
        "donor_id",
        "batch",
        "tissue_source",
        "eye",
        "sex",
        "age",
        "diabetic_status",
        "indication_or_cod",
        "hospital",
        "group_paper",
    ],
)

# Canonical retinal markers informed by Figure 1C & retinal consensus
CANONICAL_MARKERS = {
    "Rod": ["RHO", "PDE6G", "GNAT1", "NR2E3"],
    "Cone": ["ARR3", "OPN1LW", "OPN1MW", "GNAT2"],
    "BC": ["GRM6", "VSX2", "TRPM1", "PCP2"],
    "MGC": ["DKK3", "RLBP1", "SLC1A3", "GLUL"],
    "Microglia": ["TYROBP", "C1QA", "C1QB", "CSF1R"],
    "AC": ["SLC32A1", "GAD1", "GAD2", "TFAP2A"],
    "HC": ["TNR", "ONECUT1", "ONECUT2", "LHX1"],
    "TC": ["CD3D", "CD3E", "TRAC", "CD52"],
    "Astrocyte": ["GFAP", "AQP4", "S100B"],
    "RGC": ["RBPMS", "SNCG", "THY1", "NEFL"],
    "RPE": ["RPE65", "BEST1", "LRAT"],
}


def _curl_bytes(url: str) -> bytes:
    return subprocess.check_output(["curl", "-sSL", "--retry", "3", "--retry-delay", "2", url])


def _aggregate_figshare_samples() -> tuple[
    np.ndarray,
    dict[tuple[str, str], np.ndarray],
    dict[str, np.ndarray],
    dict[tuple[str, str], int],
    pd.DataFrame,
]:
    art = json.loads(_curl_bytes("https://api.figshare.com/v2/articles/28757048").decode("utf-8"))
    folders = art["folder_structure"]
    by_folder: dict[str, dict[str, dict]] = {}
    for fl in art["files"]:
        sname = folders.get(str(fl["id"]), "")
        if sname:
            by_folder.setdefault(sname, {})[fl["name"]] = fl

    unique_genes: np.ndarray | None = None
    proj_mat: sp.csr_matrix | None = None
    gene_to_idx: dict[str, list[int]] = {}
    mt_idx: np.ndarray | None = None

    pb_ct: dict[tuple[str, str], np.ndarray] = {}
    pb_whole: dict[str, np.ndarray] = {}
    cell_counts: dict[tuple[str, str], int] = {}
    sample_qc_rows = []

    for _, row in TABLE1_SAMPLES.iterrows():
        sfolder = str(row["figshare_folder"])
        donor_id = str(row["donor_id"])
        files = by_folder[sfolder]

        feat_lines = gzip.decompress(_curl_bytes(files["features.tsv.gz"]["download_url"])).decode("utf-8", "replace").splitlines()
        raw_syms = np.array([ln.split("\t")[1].strip() for ln in feat_lines if ln.strip()])

        if unique_genes is None:
            unique_genes, inv = np.unique(raw_syms, return_inverse=True)
            proj_mat = sp.csr_matrix(
                (np.ones(len(raw_syms), dtype=np.float64), (inv, np.arange(len(raw_syms)))),
                shape=(len(unique_genes), len(raw_syms)),
            )
            for idx_u, gsym in enumerate(unique_genes):
                gene_to_idx[str(gsym)] = [idx_u]
            mt_idx = np.where(np.char.startswith(unique_genes.astype(str), "MT-"))[0]
        elif len(raw_syms) != proj_mat.shape[1]:
            raise RuntimeError(f"Feature length mismatch in {sfolder}: {len(raw_syms)} vs {proj_mat.shape[1]}")

        mtx_bytes = gzip.decompress(_curl_bytes(files["matrix.mtx.gz"]["download_url"]))
        mat = scipy.io.mmread(io.BytesIO(mtx_bytes)).tocsc()
        del mtx_bytes

        if mat.nnz > 0 and not np.all(np.equal(np.mod(mat.data, 1.0), 0.0)):
            raise ValueError(f"Non-integer counts detected in {sfolder}")

        # Project raw features to unique symbols (cells x unique_genes -> sparse u_genes x cells)
        mat_u = (proj_mat @ mat.astype(np.float64)).tocsc()
        del mat

        umi = np.asarray(mat_u.sum(axis=0)).ravel()
        n_genes_det = np.diff(mat_u.indptr)
        mt_umi = np.asarray(mat_u[mt_idx, :].sum(axis=0)).ravel()
        pct_mt = 100.0 * mt_umi / np.maximum(umi, 1.0)

        qc = (n_genes_det >= 200) & (umi > 500.0) & (pct_mt < 15.0)
        n_raw = mat_u.shape[1]
        n_qc = int(qc.sum())

        mat_q = mat_u[:, qc].tocsr()
        umi_q = umi[qc]
        del mat_u

        def _score_gene(sym: str) -> np.ndarray:
            idxs = gene_to_idx.get(sym, [])
            if not idxs:
                return np.zeros(n_qc, dtype=np.float64)
            v = np.asarray(mat_q[idxs[0], :].todense()).ravel()
            return np.log1p(1e4 * v / umi_q)

        xist_mean = float(_score_gene("XIST").mean()) if n_qc > 0 else 0.0
        y_mean = float(0.5 * (_score_gene("RPS4Y1").mean() + _score_gene("DDX3Y").mean())) if n_qc > 0 else 0.0
        inferred_sex = "M" if y_mean > xist_mean else "F"

        scores_df = pd.DataFrame({ct: np.mean([_score_gene(g) for g in gs], axis=0) for ct, gs in CANONICAL_MARKERS.items()})
        mx = scores_df.max(axis=1)
        labels = scores_df.idxmax(axis=1).where(mx >= 0.5, "Unassigned").to_numpy(dtype=object)

        # Aggregate QC-filtered whole retina and per marker class into donor accumulators
        whole_sum = np.asarray(mat_q.sum(axis=1)).ravel()
        if donor_id in pb_whole:
            pb_whole[donor_id] += whole_sum
        else:
            pb_whole[donor_id] = whole_sum.copy()

        for ct_val in np.unique(labels):
            mask_ct = labels == ct_val
            key = (donor_id, str(ct_val))
            ct_sum = np.asarray(mat_q[:, mask_ct].sum(axis=1)).ravel()
            cell_counts[key] = cell_counts.get(key, 0) + int(mask_ct.sum())
            if key in pb_ct:
                pb_ct[key] += ct_sum
            else:
                pb_ct[key] = ct_sum.copy()

        sample_qc_rows.append(
            {
                "sample_num": int(row["sample_num"]),
                "figshare_folder": sfolder,
                "donor_id": donor_id,
                "sex_table1": str(row["sex"]),
                "xist_mean": xist_mean,
                "rps4y1_ddx3y_mean": y_mean,
                "inferred_sex": inferred_sex,
                "sex_match": bool(inferred_sex == str(row["sex"])),
                "n_raw_barcodes": int(n_raw),
                "n_qc_cells": n_qc,
            }
        )
        print(
            f"  processed {sfolder} ({donor_id} {row['sex']} {int(row['age'])}y {row['diabetic_status']}): "
            f"raw={n_raw} qc={n_qc} sex_check={inferred_sex}/{row['sex']}",
            flush=True,
        )
        del mat_q, scores_df
        gc.collect()

    sample_qc_df = pd.DataFrame(sample_qc_rows)
    # Any sample whose expression-inferred sex disagrees with Table 1 is recorded, not dropped.
    # Primary models keep Table 1 sex; a sensitivity suite uses expression-verified donor sex.
    bad = sample_qc_df[~sample_qc_df["sex_match"]]
    if len(bad):
        print(f"WARNING: expression-inferred sex disagrees with Table 1 for:\n{bad}", flush=True)
    # Donor-level consistency: all eyes of one donor must agree on inferred sex
    per_donor = sample_qc_df.groupby("donor_id")["inferred_sex"].nunique()
    if bool((per_donor > 1).any()):
        raise RuntimeError(f"Eyes of one donor disagree on inferred sex: {per_donor[per_donor > 1]}")

    return unique_genes.astype(object), pb_ct, pb_whole, cell_counts, sample_qc_df


def _aggregate_cached() -> tuple:
    cache_dir = DATA_TMP_ROOT / "organ09_retina"
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache = cache_dir / "aggregates.npz"
    qc_csv = cache_dir / "samples_qc.csv"
    if cache.exists() and qc_csv.exists():
        z = np.load(cache, allow_pickle=True)
        genes = z["genes"]
        pb_ct = {tuple(k.split("||")): v for k, v in zip(z["ct_keys"], z["ct_vals"])}
        pb_whole = {k: v for k, v in zip(z["wh_keys"], z["wh_vals"])}
        cell_counts = {tuple(k.split("||")): int(v) for k, v in zip(z["cc_keys"], z["cc_vals"])}
        return genes, pb_ct, pb_whole, cell_counts, pd.read_csv(qc_csv)
    genes, pb_ct, pb_whole, cell_counts, sample_qc_df = _aggregate_figshare_samples()
    np.savez(
        cache,
        genes=genes,
        ct_keys=np.array(["||".join(k) for k in pb_ct.keys()], dtype=object),
        ct_vals=np.stack(list(pb_ct.values())),
        wh_keys=np.array(list(pb_whole.keys()), dtype=object),
        wh_vals=np.stack(list(pb_whole.values())),
        cc_keys=np.array(["||".join(k) for k in cell_counts.keys()], dtype=object),
        cc_vals=np.array(list(cell_counts.values()), dtype=np.int64),
    )
    sample_qc_df.to_csv(qc_csv, index=False)
    return genes, pb_ct, pb_whole, cell_counts, sample_qc_df


def _build_donor_meta(sample_qc_df: pd.DataFrame, cell_counts_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for donor_id, grp in TABLE1_SAMPLES.groupby("donor_id", sort=False):
        r0 = grp.iloc[0]
        sq = sample_qc_df[sample_qc_df["donor_id"] == donor_id]
        rows.append(
            {
                "donor_id": str(donor_id),
                "n_eyes_sampled": int(len(grp)),
                "sample_nums": ";".join(str(x) for x in grp["sample_num"].tolist()),
                "batch": str(r0["batch"]),
                "age": float(r0["age"]),
                "age_per_decade": float(r0["age"]) / 10.0,
                "sex": str(r0["sex"]),
                "sex_male": 1.0 if str(r0["sex"]) == "M" else 0.0,
                "sex_expr": str(sq["inferred_sex"].iloc[0]),
                "sex_expr_male": 1.0 if str(sq["inferred_sex"].iloc[0]) == "M" else 0.0,
                "xist_mean": float(sq["xist_mean"].mean()),
                "rps4y1_ddx3y_mean": float(sq["rps4y1_ddx3y_mean"].mean()),
                "tissue_source": str(r0["tissue_source"]),
                "is_live": 1.0 if str(r0["tissue_source"]) == "Live" else 0.0,
                "diabetic_status": str(r0["diabetic_status"]),
                "is_dr": 1.0 if str(r0["diabetic_status"]) == "DR" else 0.0,
                "is_dm_no_dr": 1.0 if str(r0["diabetic_status"]) == "DM" else 0.0,
                "is_any_dm": 1.0 if str(r0["diabetic_status"]) in {"DM", "DR"} else 0.0,
                "indication_or_cod": str(r0["indication_or_cod"]),
                "n_raw_barcodes": int(sq["n_raw_barcodes"].sum()),
                "n_qc_cells": int(sq["n_qc_cells"].sum()),
            }
        )
    df = pd.DataFrame(rows).sort_values(["age", "donor_id"]).set_index("donor_id", drop=False)
    return df.join(cell_counts_df)


def _run_suite(
    analysis_name: str,
    genes: np.ndarray,
    donor_subset: pd.DataFrame,
    pb_ct: dict[tuple[str, str], np.ndarray],
    pb_whole: dict[str, np.ndarray],
    cell_counts_df: pd.DataFrame,
    candidate_covariates: tuple[str, ...],
    min_cells_per_donor_group: int,
    min_donors: int,
    min_residual_df: int,
    run_cell_types: bool = True,
    whole_only: bool = False,
) -> dict[str, pd.DataFrame]:
    donors = donor_subset.index.tolist()
    cc_sub = cell_counts_df.loc[donors].copy()
    valid_cts = [c for c in cc_sub.columns if c != "Unassigned"]

    fam_a = pd.DataFrame()
    if run_cell_types and not whole_only and len(valid_cts) >= 2 and len(donors) >= 3:
        covs_a = tuple(c for c in candidate_covariates if c != "log_lib")
        fam_a = run_family_a_clr(
            cc_sub[valid_cts],
            donor_subset,
            age_col="age_per_decade",
            candidate_covariates=covs_a,
            min_donors_present=min(min_donors, len(donors)),
            min_total_cells=20,
            min_residual_df=min_residual_df,
        )
        if not fam_a.empty:
            fam_a.insert(0, "analysis", analysis_name)

    expr_by_group: dict[str, pd.DataFrame] = {}
    meta_by_group: dict[str, pd.DataFrame] = {}

    if whole_only:
        raw_df = pd.DataFrame({d: pb_whole[d] for d in donors}, index=genes)
        c_series = donor_subset["n_qc_cells"].loc[donors]
        expr_g, meta_g, _ = finalize_pseudobulk_group(
            raw_df,
            donor_subset.loc[donors],
            c_series,
            allow_fractional_expected=False,
        )
        expr_by_group["All_Retina_QC"] = expr_g
        meta_by_group["All_Retina_QC"] = meta_g
    else:
        for ct in valid_cts:
            eligible = [d for d in donors if int(cc_sub.loc[d, ct]) >= min_cells_per_donor_group]
            if len(eligible) < min_donors:
                continue
            raw_df = pd.DataFrame({d: pb_ct[(d, ct)] for d in eligible}, index=genes)
            c_series = cc_sub.loc[eligible, ct]
            expr_g, meta_g, _ = finalize_pseudobulk_group(
                raw_df,
                donor_subset.loc[eligible],
                c_series,
                allow_fractional_expected=False,
            )
            expr_by_group[ct] = expr_g
            meta_by_group[ct] = meta_g

    p1 = load_p1_clean_human(REPO_ROOT)
    fam_b_top = pd.DataFrame()
    fam_b_summary = pd.DataFrame()
    fam_b_conc = pd.DataFrame()
    fam_c = pd.DataFrame()
    if expr_by_group:
        fam_b_genes, fam_b_summary, fam_b_conc = run_family_b_genes_and_concordance(
            expr_by_group,
            meta_by_group,
            age_col="age_per_decade",
            candidate_covariates=candidate_covariates,
            min_donors=min_donors,
            min_cpm=1.0,
            min_expr_frac=0.5,
            min_residual_df=min_residual_df,
            is_unbiased_single_cell=bool((not whole_only) and len(expr_by_group) >= 2),
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
            p1_program=p1,
            age_col="age_per_decade",
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


def main() -> None:
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
    genes, pb_ct, pb_whole, cell_counts, sample_qc_df = _aggregate_cached()

    all_donors = sorted(TABLE1_SAMPLES["donor_id"].unique())
    all_cts = sorted({ct for _, ct in cell_counts.keys()})
    cell_counts_df = pd.DataFrame(0, index=pd.Index(all_donors, name="donor_id"), columns=all_cts, dtype=int)
    for (d, ct), n in cell_counts.items():
        cell_counts_df.loc[d, ct] = n

    donor_meta = _build_donor_meta(sample_qc_df, cell_counts_df)
    gidx = {str(g): i for i, g in enumerate(genes)}
    y_genes = ["RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY"]
    donor_meta["xist_umi_qc"] = [int(pb_whole[d][gidx["XIST"]]) for d in donor_meta.index]
    donor_meta["tsix_umi_qc"] = [int(pb_whole[d][gidx["TSIX"]]) if "TSIX" in gidx else -1 for d in donor_meta.index]
    donor_meta["y5_umi_qc"] = [int(sum(pb_whole[d][gidx[g]] for g in y_genes if g in gidx)) for d in donor_meta.index]
    donor_meta["sex_table1_vs_expr_match"] = donor_meta["sex"] == donor_meta["sex_expr"]

    non_dm_meta = donor_meta[donor_meta["diabetic_status"] == "Non-DM"].copy()
    cadaver_meta = donor_meta[donor_meta["tissue_source"] == "Cadaver"].copy()

    suites = [
        _run_suite(
            "All12_DiseaseAdj",
            genes,
            donor_meta,
            pb_ct,
            pb_whole,
            cell_counts_df,
            candidate_covariates=("is_dr", "is_dm_no_dr", "sex_male", "log_lib"),
            min_cells_per_donor_group=10,
            min_donors=6,
            min_residual_df=2,
            run_cell_types=True,
            whole_only=False,
        ),
        _run_suite(
            "All12_DiseaseAdj_ExprSex",
            genes,
            donor_meta,
            pb_ct,
            pb_whole,
            cell_counts_df,
            candidate_covariates=("is_dr", "is_dm_no_dr", "sex_expr_male", "log_lib"),
            min_cells_per_donor_group=10,
            min_donors=6,
            min_residual_df=2,
            run_cell_types=True,
            whole_only=False,
        ),
        _run_suite(
            "All12_DiseaseLiveAdj",
            genes,
            donor_meta,
            pb_ct,
            pb_whole,
            cell_counts_df,
            candidate_covariates=("is_dr", "is_dm_no_dr", "is_live", "sex_male", "log_lib"),
            min_cells_per_donor_group=10,
            min_donors=6,
            min_residual_df=2,
            run_cell_types=True,
            whole_only=False,
        ),
        _run_suite(
            "All12_WholeRetina_DiseaseAdj",
            genes,
            donor_meta,
            pb_ct,
            pb_whole,
            cell_counts_df,
            candidate_covariates=("is_dr", "is_dm_no_dr", "is_live", "sex_male", "log_lib"),
            min_cells_per_donor_group=10,
            min_donors=6,
            min_residual_df=2,
            run_cell_types=False,
            whole_only=True,
        ),
        _run_suite(
            "NonDM6_AgeOnly",
            genes,
            non_dm_meta,
            pb_ct,
            pb_whole,
            cell_counts_df,
            candidate_covariates=("log_lib",),
            min_cells_per_donor_group=10,
            min_donors=5,
            min_residual_df=2,
            run_cell_types=True,
            whole_only=False,
        ),
        _run_suite(
            "Cadaver7_DiseaseAdj",
            genes,
            cadaver_meta,
            pb_ct,
            pb_whole,
            cell_counts_df,
            candidate_covariates=("is_dr", "is_dm_no_dr", "log_lib"),
            min_cells_per_donor_group=10,
            min_donors=5,
            min_residual_df=2,
            run_cell_types=True,
            whole_only=False,
        ),
    ]

    fam_a_df = pd.concat([s["fam_a"] for s in suites if not s["fam_a"].empty], ignore_index=True)
    fam_b_sum_df = pd.concat([s["fam_b_summary"] for s in suites if not s["fam_b_summary"].empty], ignore_index=True)
    fam_b_top_df = pd.concat([s["fam_b_top"] for s in suites if not s["fam_b_top"].empty], ignore_index=True)
    fam_b_conc_df = pd.concat([s["fam_b_conc"] for s in suites if not s["fam_b_conc"].empty], ignore_index=True)
    fam_c_df = pd.concat([s["fam_c"] for s in suites if not s["fam_c"].empty], ignore_index=True)

    access_df = pd.DataFrame(
        [
            {
                "accession": "HRA014708",
                "repository": "CNCB-NGDC GSA-Human",
                "url": "https://ngdc.cncb.ac.cn/gsa-human/browse/HRA014708",
                "access_status": "Controlled access; HTTP browse gated by SecurityCheck HTML interstitial (not bypassed, no application submitted)",
                "usable_raw_counts_via_figshare": f"Yes (doi:10.6084/m9.figshare.28757048, 18 Cell Ranger v8.0.1 count matrix triplets, {int(sample_qc_df['n_raw_barcodes'].sum())} raw barcodes, {int(sample_qc_df['n_qc_cells'].sum())} QC cells across 12 unique Table 1 donors D1 to D12)",
                "missing_artifact": "Integrated scVI/Leiden-annotated /RPE.h5ad referenced in Figshare Script.py is absent",
            }
        ]
    )

    results_dir = REPO_ROOT / "results"
    reports_dir = REPO_ROOT / "reports"
    prefix = "followup-organ-09-eye-hra014708"

    write_clean_csv(donor_meta.reset_index(drop=True), results_dir / f"{prefix}-donors.csv")
    write_clean_csv(sample_qc_df, results_dir / f"{prefix}-samples-qc.csv")
    write_clean_csv(access_df, results_dir / f"{prefix}-access-check.csv")
    write_clean_csv(fam_a_df, results_dir / f"{prefix}-fractions.csv")
    write_clean_csv(fam_b_sum_df, results_dir / f"{prefix}-gene-summary.csv")
    write_clean_csv(fam_b_top_df, results_dir / f"{prefix}-top-genes.csv")
    write_clean_csv(fam_b_conc_df, results_dir / f"{prefix}-concordance.csv")
    write_clean_csv(fam_c_df, results_dir / f"{prefix}-p1-signature.csv")

    donor_tbl = "\n".join(
        f"| {r['donor_id']} | {int(r['age'])} | {r['sex']} | {r['diabetic_status']} | {r['tissue_source']} | {int(r['n_eyes_sampled'])} ({r['sample_nums']}) | {int(r['n_qc_cells'])} | {r['indication_or_cod']} |"
        for _, r in donor_meta.iterrows()
    )
    fa_prim = fam_a_df[fam_a_df["analysis"] == "All12_DiseaseAdj"]
    fa_tbl = "\n".join(
        f"| {r['cell_type']} | {int(r['n_donors_present'])} | {int(r['total_cells'])} | {r['median_fraction']:.4f} | {r['effect_clr']:+.4f} | {r['pvalue']:.4f} | {r['q_bh_family_a']:.4f} |"
        for _, r in fa_prim.iterrows()
    )
    fb_prim = fam_b_sum_df[fam_b_sum_df["analysis"] == "All12_DiseaseAdj"]
    fb_tbl = "\n".join(
        f"| {r['group']} | {int(r['n_donors'])} | {int(r['n_genes_tested'])} | {int(r['n_sig_within_q05'])} | {int(r['n_sig_global_q05'])} | {r['min_pvalue']:.2e} | {r['min_q_bh_global']:.4f} |"
        for _, r in fb_prim.iterrows()
    )
    conc_tbl = "\n".join(
        f"| {r['analysis']} | {int(r['n_cell_types'])} | {int(r['n_common_genes'])} | {r['mean_pairwise_pearson_r']:+.4f} | {r['perm_p_one_sided']:.4f} | {r['perm_p_two_sided']:.4f} |"
        for _, r in fam_b_conc_df.iterrows()
    )
    fc_prim = fam_c_df[fam_c_df["analysis"].isin(["All12_DiseaseAdj", "All12_WholeRetina_DiseaseAdj", "NonDM6_AgeOnly"])]
    fc_tbl = "\n".join(
        f"| {r['analysis']} | {r['group']} | {int(r['n_donors'])} | {int(r['n_p1_up_detected'])}/{int(r['n_p1_down_detected'])} | {r['composite_effect']:+.4f} | {r['composite_pvalue']:.4f} | {r['empirical_null_pvalue']:.4f} | {r['q_bh_composite']:.4f} | {r['q_bh_empirical']:.4f} |"
        for _, r in fc_prim.iterrows()
    )

    n_a_sig = int((fa_prim["q_bh_family_a"] < 0.05).sum()) if len(fa_prim) else 0
    mb_prim = int(fb_prim["n_genes_tested"].sum())
    n_b_sig = int(fb_prim["n_sig_global_q05"].sum())

    def _hits(name: str) -> str:
        sub = fam_b_sum_df[fam_b_sum_df["analysis"] == name]
        if sub.empty:
            return "not run (no group met donor thresholds)"
        return f"{int(sub['n_sig_global_q05'].sum())} of {int(sub['n_genes_tested'].sum())} (min q_BH = {float(sub['min_q_bh_global'].min()):.4f})"

    def _conc(name: str) -> str:
        sub = fam_b_conc_df[fam_b_conc_df["analysis"] == name] if len(fam_b_conc_df) else fam_b_conc_df
        if sub is None or len(sub) == 0:
            return "not computed"
        r = sub.iloc[0]
        return f"r = {float(r['mean_pairwise_pearson_r']):+.4f}, one-sided permutation p = {float(r['perm_p_one_sided']):.4f}"

    def _fa_hits(name: str) -> str:
        sub = fam_a_df[fam_a_df["analysis"] == name] if len(fam_a_df) else fam_a_df
        if sub is None or len(sub) == 0:
            return "not run"
        return f"{int((sub['q_bh_family_a'] < 0.05).sum())} of {len(sub)}"

    min_q_c_comp = float(fam_c_df["q_bh_composite"].min()) if len(fam_c_df) else np.nan
    min_q_c_emp = float(fam_c_df["q_bh_empirical"].min()) if len(fam_c_df) else np.nan
    n_c_both = int(((fam_c_df["q_bh_composite"] < 0.05) & (fam_c_df["q_bh_empirical"] < 0.05)).sum()) if len(fam_c_df) else 0
    n_c_pos_both = int(((fam_c_df["composite_effect"] > 0) & (fam_c_df["q_bh_composite"] < 0.05) & (fam_c_df["q_bh_empirical"] < 0.05)).sum()) if len(fam_c_df) else 0
    tot_raw = int(sample_qc_df["n_raw_barcodes"].sum())
    tot_qc = int(sample_qc_df["n_qc_cells"].sum())
    n_match = int(sample_qc_df["sex_match"].sum())
    n_samp = int(len(sample_qc_df))
    mism = donor_meta[~donor_meta["sex_table1_vs_expr_match"]]
    n_m_t1 = int((donor_meta["sex"] == "M").sum())
    n_f_t1 = int((donor_meta["sex"] == "F").sum())
    n_m_ex = int((donor_meta["sex_expr"] == "M").sum())
    n_f_ex = int((donor_meta["sex_expr"] == "F").sum())
    mism_txt = "; ".join(
        f"donor {r['donor_id']} (sample(s) {r['sample_nums']}, {int(r['age'])} yr): Table 1 sex {r['sex']}, expression sex {r['sex_expr']} (QC-cell pseudobulk XIST = {int(r['xist_umi_qc'])} UMIs, TSIX = {int(r['tsix_umi_qc'])} UMIs, RPS4Y1 + DDX3Y + UTY + KDM5D + EIF1AY = {int(r['y5_umi_qc'])} UMIs)"
        for _, r in mism.iterrows()
    ) or "none"

    donor_tbl = "\n".join(
        f"| {r['donor_id']} | {int(r['age'])} | {r['sex']} | {r['sex_expr']} | {int(r['xist_umi_qc'])} | {int(r['y5_umi_qc'])} | {r['diabetic_status']} | {r['tissue_source']} | {int(r['n_eyes_sampled'])} ({r['sample_nums']}) | {int(r['n_qc_cells'])} | {r['indication_or_cod']} |"
        for _, r in donor_meta.iterrows()
    )
    fb_all_tbl = "\n".join(
        f"| {r['analysis']} | {r['group']} | {int(r['n_donors'])} | {int(r['n_genes_tested'])} | {int(r['residual_df'])} | {int(r['n_sig_within_q05'])} | {int(r['n_sig_global_q05'])} | {r['min_pvalue']:.2e} | {r['min_q_bh_global']:.4f} |"
        for _, r in fam_b_sum_df.iterrows()
        if str(r.get("status", "ok")) == "ok" or pd.notna(r.get("n_genes_tested"))
    )
    fc_all_tbl = "\n".join(
        f"| {r['analysis']} | {r['group']} | {int(r['n_donors'])} | {int(r['n_p1_up_detected'])}/{int(r['n_p1_down_detected'])} | {r['composite_effect']:+.4f} | {r['composite_pvalue']:.4f} | {r['empirical_null_pvalue']:.4f} | {r['q_bh_composite']:.4f} | {r['q_bh_empirical']:.4f} |"
        for _, r in fam_c_df.iterrows()
    )
    suite_names = ["All12_DiseaseAdj", "All12_DiseaseAdj_ExprSex", "All12_DiseaseLiveAdj", "All12_WholeRetina_DiseaseAdj", "NonDM6_AgeOnly", "Cadaver7_DiseaseAdj"]
    suite_tbl = "\n".join(
        f"| {s} | {_fa_hits(s)} | {_hits(s)} | {_conc(s)} |" for s in suite_names
    )
    ages = donor_meta["age"].astype(int)

    notes = []
    conc_sig = fam_b_conc_df[fam_b_conc_df["perm_p_one_sided"] < 0.05] if len(fam_b_conc_df) else fam_b_conc_df
    if len(conc_sig):
        n_conc = int(len(fam_b_conc_df))
        notes.append(
            "- Cross-cell-type concordance reached one-sided permutation p < 0.05 in "
            + ", ".join(f"`{r['analysis']}` (r = {r['mean_pairwise_pearson_r']:+.4f}, p = {r['perm_p_one_sided']:.4f})" for _, r in conc_sig.iterrows())
            + f". These p values are not corrected across the {n_conc} concordance tests in this report, and the two donor subsets that remove the disease or surgical strata ("
            + ", ".join(f"`{r['analysis']}` p = {r['perm_p_one_sided']:.4f}" for _, r in fam_b_conc_df[fam_b_conc_df['perm_p_one_sided'] >= 0.05].iterrows())
            + ") did not reach p < 0.05. Shared donor-level factors that the covariates capture only coarsely (diabetic retinopathy, enucleation for ocular disease, cadaver handling) and ambient RNA from the dominant rod population can align age effects across marker-scored classes. This is a flag for review against the paper's human concordance result (p = 0.70 in the surveyed human datasets), not a replication of concordance: it comes from one 12-donor cohort with disease confounding."
        )
    hits = fam_b_top_df[fam_b_top_df["q_bh_global"] < 0.05] if len(fam_b_top_df) else fam_b_top_df
    for _, r in hits.iterrows():
        notes.append(
            f"- Single Family B hit: `{r['analysis']}` {r['group']} {r['gene']} (effect {r['effect']:+.4f} log2CPM per decade, p = {r['pvalue']:.2e}, q_BH global = {r['q_bh_global']:.4f}, n = {int(r['n_donors'])}, residual df = {int(r['residual_df'])}). With residual df = {int(r['residual_df'])} this is treated as unstable and is not replicated in the 12-donor suites."
        )
    emp_only = fam_c_df[(fam_c_df["q_bh_empirical"] < 0.05) & (fam_c_df["q_bh_composite"] >= 0.05)] if len(fam_c_df) else fam_c_df
    if len(emp_only):
        notes.append(
            "- `P1_Clean` groups significant on the empirical null only (analytic q_BH >= 0.05): "
            + "; ".join(f"`{r['analysis']}` {r['group']} (composite {r['composite_effect']:+.4f}, q_BH empirical = {r['q_bh_empirical']:.4f}, q_BH analytic = {r['q_bh_composite']:.4f})" for _, r in emp_only.iterrows())
            + ". All of these are in the negative direction (opposite to mouse-signature transport) and none passes both criteria."
            if bool((emp_only["composite_effect"] < 0).all())
            else "- `P1_Clean` groups significant on the empirical null only (analytic q_BH >= 0.05): "
            + "; ".join(f"`{r['analysis']}` {r['group']} (composite {r['composite_effect']:+.4f}, q_BH empirical = {r['q_bh_empirical']:.4f}, q_BH analytic = {r['q_bh_composite']:.4f})" for _, r in emp_only.iterrows())
            + ". None passes both criteria."
        )
    notable_txt = "\n".join(notes) if notes else "- None."

    report_md = f"""# Dataset 09: Human Retinal Aging Cohort (Ophthalmology Science 2026, PMC12892068, HRA014708 / Figshare 28757048)

## Verdict
Across the 12 unique Table 1 donors (18 retinal samples, ages {ages.min()} to {ages.max()} yr; Table 1 sex {n_m_t1}M / {n_f_t1}F, expression-verified sex {n_m_ex}M / {n_f_ex}F; 4 DR donors, 2 DM-without-DR donors, 6 non-diabetic donors), with paired eyes collapsed to the donor and diabetic status adjusted, the primary suite `All12_DiseaseAdj` gave {n_a_sig} of {len(fa_prim)} cell-class fraction shifts at q_BH < 0.05, {n_b_sig} of {mb_prim} global gene-level age hits at q_BH < 0.05, cross-cell-type concordance {_conc('All12_DiseaseAdj')}, and {n_c_both} of {len(fam_c_df)} `P1_Clean` signature tests (all suites) at q_BH < 0.05 on both the analytic and the 200-draw matched null ({n_c_pos_both} in the positive, transport-supporting direction; min q_BH analytic = {min_q_c_comp:.4f}, min q_BH empirical = {min_q_c_emp:.4f}). With n = 12 and heavy disease and surgical confounding, these are low-power results and do not establish absence of age effects.

## Notable Results and Caveats
{notable_txt}

## Primary Source and Dual-Accession Audit
- Publication: Ophthalmology Science (2026), DOI `10.1016/j.xops.2025.101062`, PMID `41684508`, PMCID `PMC12892068`.
- GSA-Human `HRA014708`: controlled access. The browse page returned an anti-bot SecurityCheck interstitial; it was not bypassed and no access application was submitted.
- Figshare `doi:10.6084/m9.figshare.28757048`: openly hosts 18 Cell Ranger v8.0.1 count matrix triplets (`sample_1` to `sample_18`) plus `README.docx` and `Script.py`. The author-clustered `/RPE.h5ad` referenced by `Script.py` is not deposited. Counts are integer UMIs: {tot_raw:,} raw barcodes, {tot_qc:,} cells passing the paper's QC (at least 200 genes, more than 500 UMIs, under 15% mitochondrial UMIs).
- Cell classes were assigned by canonical marker scoring informed by Figure 1C (mean log1p(CP10k) of the marker set at least 0.5; otherwise Unassigned), because the author labels are not deposited. An unpartitioned whole-retina donor pseudobulk (`All12_WholeRetina_DiseaseAdj`) is reported as a label-free control.

## Metadata Discrepancies Found (reported, not resolved by assumption)
1. Donor sex: expression-inferred sex (mean log1p(CP10k) of XIST versus mean of RPS4Y1 and DDX3Y) matched Table 1 on {n_match} of {n_samp} samples. Mismatch: {mism_txt}. Primary models use Table 1 sex as printed; `All12_DiseaseAdj_ExprSex` repeats the primary model with expression-verified sex.
2. Cohort counts: the Figshare description states 18 retinal samples from 13 Chinese donors (8 male, 5 female, ages 34 to 92 years) with 7 DR, 3 diabetics without retinopathy, and 8 non-diabetic samples, whereas Table 1 of PMC12892068 lists 18 samples from 12 donors (D1 to D12; Table 1 sex {n_m_t1}M / {n_f_t1}F) with 8 DR, 3 DM, and 7 non-DM samples. The source of the difference cannot be determined from the deposited files. This analysis uses the Table 1 donor mapping (12 donor units); the cohort-level counts here therefore follow Table 1, and the sample-level sex check above is the only independent check applied.

## Donor Summary (n = 12 donors across 18 eyes)
| Donor | Age (yr) | Sex (Table 1) | Sex (expression) | XIST UMIs (QC cells) | Y-gene UMIs (5 genes) | Diabetic Status | Source | Eyes (Sample #) | QC Cells | Surgical Indication or Cause of Death |
|---|---:|:---:|:---:|---:|---:|---|---|---|---:|---|
{donor_tbl}

## Multiple-Testing Families
- Family A (per suite): all cell classes passing presence thresholds (present in at least the suite's minimum donor count, at least 20 total cells), BH across classes within the suite. Model CLR fraction ~ age_per_decade + suite covariates (log_lib excluded).
- Family B (per suite): all gene by cell-class tests in the suite, BH applied globally across every gene in every class of that suite (q_BH global); within-class BH also reported. Genes kept if CPM at least 1 in at least half of donors. Concordance: mean pairwise Pearson r of age effects across classes, 1,000 donor-label permutations.
- Family C (per suite): one `P1_Clean` composite test per cell class (118 human orthologs, 76 up, 42 down), BH across classes within the suite for the analytic p and separately for the 200-draw expression-matched random-set empirical p.

## Suite Overview
| Suite | Family A q<0.05 | Family B global q<0.05 | Concordance |
|---|---|---|---|
{suite_tbl}

Suite definitions: `All12_DiseaseAdj` covariates is_dr, is_dm_no_dr, sex_male (Table 1), log_lib; `All12_DiseaseAdj_ExprSex` same with expression-verified sex; `All12_DiseaseLiveAdj` adds is_live; `All12_WholeRetina_DiseaseAdj` unpartitioned whole-retina pseudobulk with is_dr, is_dm_no_dr, is_live, sex_male, log_lib; `NonDM6_AgeOnly` the 6 non-diabetic donors with log_lib only; `Cadaver7_DiseaseAdj` the cadaver donors with is_dr, is_dm_no_dr, log_lib. Cell classes require at least 10 cells per donor.

## Family A: Primary Suite Detail (`All12_DiseaseAdj`)
| Cell Class | Donors Present | QC Cells | Median Fraction | Effect (CLR/dec) | p-value | q_BH (Family A) |
|---|---:|---:|---:|---:|---:|---:|
{fa_tbl}

## Family B: All Suites
| Suite | Group | n_donors | Genes Tested | Residual df | Within q<0.05 | Global q<0.05 | Min p | Min q_BH Global |
|---|---|---:|---:|---:|---:|---:|---:|---:|
{fb_all_tbl}

## Family C: P1_Clean Mouse Signature, All Suites
| Suite | Group | n_donors | P1 Up/Down | Composite Effect | Analytic p | Empirical p (200 draws) | q_BH (Analytic) | q_BH (Empirical) |
|---|---|---:|:---:|---:|---:|---:|---:|---:|
{fc_all_tbl}

## Key Limitations
1. Missing author cell labels and controlled raw reads; marker scoring approximates cell classes.
2. Confounding: 4 of 12 donors have DR, 2 have DM without DR, and all 5 living donors had eyes enucleated for ocular disease (phthisis bulbi, corneal perforation, ulcer, or staphyloma). Disease is adjusted as covariates, but with n = 12 the adjustment is coarse.
3. Sample size: residual df is at most 6 in the primary suite and lower in sensitivity suites and in classes with fewer donors.
4. The Table 1 sex entry for one donor disagrees with expression, and the Figshare cohort counts disagree with Table 1 (see above).

## Output Files
- `results/{prefix}-donors.csv`
- `results/{prefix}-samples-qc.csv`
- `results/{prefix}-access-check.csv`
- `results/{prefix}-fractions.csv`
- `results/{prefix}-gene-summary.csv`
- `results/{prefix}-top-genes.csv`
- `results/{prefix}-concordance.csv`
- `results/{prefix}-p1-signature.csv`
"""
    write_clean_markdown(report_md, reports_dir / f"{prefix}.md")

    blocker_md = f"""# Dataset 09 Access Audit: Ophthalmology Science 2026 (HRA014708 / Figshare 28757048)

## Status
Partial block only. GSA-Human `HRA014708` raw reads are controlled access (not applied for, not bypassed), and the author-clustered `/RPE.h5ad` is not deposited on Figshare. The donor-level analysis was run from the 18 open Figshare Cell Ranger count matrices; see [`reports/{prefix}.md`]({prefix}.md).

## Access Verification Details
1. CNCB-NGDC GSA-Human `HRA014708` (`https://ngdc.cncb.ac.cn/gsa-human/browse/HRA014708`): HTTP request returned an anti-bot SecurityCheck interstitial; access is governed by GSA-Human data access approval. No bypass and no application.
2. Figshare `doi:10.6084/m9.figshare.28757048` (`https://api.figshare.com/v2/articles/28757048`): 18 count matrix triplets ({tot_raw:,} raw barcodes, {tot_qc:,} QC cells, integer UMIs), `README.docx`, `Script.py`. Expression sex matched Table 1 on {n_match} of {n_samp} samples (mismatch: {mism_txt}).
3. The Figshare description's cohort counts (13 donors, 8 male, 5 female; 7 DR, 3 DM, 8 non-DM samples) differ from Table 1 (12 donors; 8 DR, 3 DM, 7 non-DM samples). Not resolvable from deposited files; Table 1 mapping used.
"""
    write_clean_markdown(blocker_md, reports_dir / f"{prefix}-BLOCKED.md")
    if os.environ.get("ORGAN09_KEEP_CACHE") != "1":
        clean_temp_dir(DATA_TMP_ROOT / "organ09_retina")
    print(f"Organ 09 done: raw={tot_raw}, qc={tot_qc}, FamA_sig={n_a_sig}, MB={mb_prim}, FamB_sig={n_b_sig}, min_qC_comp={min_q_c_comp:.4f}, min_qC_emp={min_q_c_emp:.4f}")


if __name__ == "__main__":
    main()
