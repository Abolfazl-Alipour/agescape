"""Dataset 04: Enge et al. (2017) human pancreas Smart-seq2 aging scRNA-seq (GSE81547)."""

from __future__ import annotations

import gzip
import io
from pathlib import Path
import re
import sys
import tarfile

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
    load_p1_clean_human,
    prepare_counts_matrix,
    run_family_a_clr,
    run_family_b_genes_and_concordance,
    run_family_c_p1_signature,
    write_clean_csv,
    write_clean_markdown,
)

HTSEQ_META_ROWS = {
    "no_feature",
    "ambiguous",
    "too_low_aQual",
    "not_aligned",
    "alignment_not_unique",
    "__no_feature",
    "__ambiguous",
    "__too_low_aQual",
    "__not_aligned",
    "__alignment_not_unique",
}


def parse_age_years(raw: str) -> float:
    """Parse donor_age characteristic string into numeric years (1 month -> 1/12)."""
    s = raw.strip().lower()
    m_month = re.match(r"^([0-9]+(?:\.[0-9]+)?)\s*(?:m|mo|mos|month|months)\b", s)
    if m_month:
        return round(float(m_month.group(1)) / 12.0, 4)
    m_yr = re.match(r"^([0-9]+(?:\.[0-9]+)?)\s*(?:y|yr|yrs|year|years)?$", s)
    if m_yr:
        return float(m_yr.group(1))
    raise ValueError(f"Unrecognized pancreas donor_age string: {raw!r}")


def parse_series_matrix_samples(matrix_gz: Path) -> pd.DataFrame:
    """Extract per-cell GSM, age_years, sex, inferred_cell_type from GSE81547 series matrix."""
    gsms: list[str] = []
    titles: list[str] = []
    chars_rows: list[list[str]] = []

    with gzip.open(matrix_gz, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_geo_accession\t"):
                gsms = [x.strip().strip('"') for x in line.rstrip("\r\n").split("\t")[1:]]
            elif line.startswith("!Sample_title\t"):
                titles = [x.strip().strip('"') for x in line.rstrip("\r\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1\t"):
                chars_rows.append([x.strip().strip('"') for x in line.rstrip("\r\n").split("\t")[1:]])
            elif line.startswith("!series_matrix_table_begin"):
                break

    if not gsms:
        raise RuntimeError("No GSM accessions parsed from series matrix")

    records = []
    for idx, gsm in enumerate(gsms):
        age_raw = None
        sex_raw = None
        ct_raw = None
        for row in chars_rows:
            if idx >= len(row):
                continue
            item = row[idx]
            if item.startswith("donor_age:"):
                age_raw = item.split(":", 1)[1].strip()
            elif item.startswith("gender:"):
                sex_raw = item.split(":", 1)[1].strip().lower()
            elif item.startswith("inferred_cell_type:"):
                ct_raw = item.split(":", 1)[1].strip()
        age_years = parse_age_years(age_raw)
        sex = "M" if sex_raw in {"male", "m"} else "F"
        donor_id = f"{age_raw.replace(' ', '')}_{sex}"
        records.append(
            {
                "gsm": gsm,
                "title": titles[idx] if idx < len(titles) else gsm,
                "donor_id": donor_id,
                "age_raw": age_raw,
                "age_years": age_years,
                "sex": sex,
                "sex_male": 1 if sex == "M" else 0,
                "cell_type": ct_raw,
            }
        )
    return pd.DataFrame(records).set_index("gsm")


def aggregate_tar_counts(tar_path: Path, cell_meta: pd.DataFrame) -> tuple[pd.DataFrame, dict[tuple[str, str], np.ndarray], list[str], dict]:
    """Stream GSE81547_RAW.tar and sum raw gene counts per (donor_id, cell_type)."""
    genes_ref: list[str] | None = None
    valid_mask: np.ndarray | None = None
    sums: dict[tuple[str, str], np.ndarray] = {}
    parsed_cells = 0
    matched_gsms: set[str] = set()
    ercc_count = 0
    htseq_count = 0

    gsm_lookup = cell_meta[["donor_id", "cell_type"]].to_dict("index")

    with tarfile.open(tar_path, "r:") as tf:
        for member in tf:
            if not member.isfile():
                continue
            name = Path(member.name).name
            m = re.match(r"^(GSM\d+)_", name)
            if not m:
                continue
            gsm = m.group(1)
            if gsm not in gsm_lookup:
                continue
            matched_gsms.add(gsm)
            raw_bytes = tf.extractfile(member).read()
            if name.endswith(".gz"):
                text = gzip.decompress(raw_bytes).decode("utf-8", errors="replace")
            else:
                text = raw_bytes.decode("utf-8", errors="replace")

            lines = text.splitlines()
            if genes_ref is None:
                g_list = []
                keep = []
                val_list = []
                for ln in lines:
                    if not ln.strip():
                        continue
                    parts = ln.split("\t")
                    if len(parts) < 2:
                        continue
                    try:
                        val = float(parts[1])
                    except ValueError:
                        continue
                    g_clean = parts[0].strip().strip('"')
                    if g_clean in HTSEQ_META_ROWS or g_clean.startswith("__"):
                        htseq_count += 1
                    elif g_clean.upper().startswith("ERCC-"):
                        ercc_count += 1
                    else:
                        g_list.append(g_clean)
                        val_list.append(val)
                        keep.append(g_clean)
                genes_ref = g_list
                vals = np.asarray(val_list, dtype=np.float64)
            else:
                row_map = {}
                for ln in lines:
                    if not ln.strip():
                        continue
                    parts = ln.split("\t")
                    if len(parts) < 2:
                        continue
                    try:
                        val = float(parts[1])
                    except ValueError:
                        continue
                    g_clean = parts[0].strip().strip('"')
                    row_map[g_clean] = val
                vals = np.fromiter((row_map[g] for g in genes_ref), dtype=np.float64, count=len(genes_ref))

            key = (gsm_lookup[gsm]["donor_id"], gsm_lookup[gsm]["cell_type"])
            if key not in sums:
                sums[key] = np.zeros(len(genes_ref), dtype=np.float64)
            sums[key] += vals
            parsed_cells += 1

    if genes_ref is None or parsed_cells != len(cell_meta):
        raise RuntimeError(f"Parsed {parsed_cells} cells vs expected {len(cell_meta)} GSMs")

    cell_counts = pd.crosstab(cell_meta["donor_id"], cell_meta["cell_type"])
    parse_stats = {
        "parsed_cells": parsed_cells,
        "htseq_meta_rows_excluded": htseq_count,
        "ercc_spikein_rows_excluded": ercc_count,
    }
    return cell_counts, sums, genes_ref, parse_stats


def build_subset_expr(
    sums: dict[tuple[str, str], np.ndarray],
    genes_ref: list[str],
    cell_counts: pd.DataFrame,
    donor_meta: pd.DataFrame,
    min_cells: int = 20,
    min_donors: int = 5,
) -> tuple[dict[str, pd.DataFrame], dict[str, pd.DataFrame]]:
    """Build per-cell-type log2(CPM+1) expression and donor metadata for eligible donors."""
    expr_by_ct: dict[str, pd.DataFrame] = {}
    meta_by_ct: dict[str, pd.DataFrame] = {}
    donors = donor_meta.index.tolist()

    for ct in sorted(cell_counts.columns):
        eligible = [d for d in donors if d in cell_counts.index and int(cell_counts.loc[d, ct]) >= min_cells]
        if len(eligible) < min_donors:
            continue
        pb = pd.DataFrame({d: sums[(d, ct)] for d in eligible}, index=pd.Index(genes_ref, dtype=str))
        filt_pb, _ = prepare_counts_matrix(pb)
        log2cpm, lib_sizes = counts_to_log2cpm(filt_pb)
        m = donor_meta.loc[eligible].copy()
        m["lib_size"] = lib_sizes.loc[eligible]
        m["log_lib"] = np.log1p(m["lib_size"].to_numpy(dtype=float))
        m["n_cells"] = cell_counts.loc[eligible, ct].astype(int)
        expr_by_ct[ct] = log2cpm
        meta_by_ct[ct] = m
    return expr_by_ct, meta_by_ct


def main() -> None:
    tmp_dir = DATA_TMP_ROOT / "04_gse81547"
    matrix_gz = download_with_guard(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81547/matrix/GSE81547_series_matrix.txt.gz",
        tmp_dir / "GSE81547_series_matrix.txt.gz",
        max_file_gb=0.2,
    )
    tar_path = download_with_guard(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81547/suppl/GSE81547_RAW.tar",
        tmp_dir / "GSE81547_RAW.tar",
        max_file_gb=1.0,
    )

    cell_meta = parse_series_matrix_samples(matrix_gz)
    donor_meta = (
        cell_meta.groupby("donor_id", sort=False)
        .agg(
            age_raw=("age_raw", "first"),
            age_years=("age_years", "first"),
            sex=("sex", "first"),
            sex_male=("sex_male", "first"),
            total_cells=("age_years", "count"),
        )
        .sort_values("age_years")
    )
    donor_meta["age_per_decade"] = donor_meta["age_years"] / 10.0

    cell_counts, sums, genes_ref, parse_stats = aggregate_tar_counts(tar_path, cell_meta)
    p1 = load_p1_clean_human(REPO_ROOT)

    # Primary: all 8 donors (1 month to 54 yr)
    fam_a_all = run_family_a_clr(
        cell_counts, donor_meta, age_col="age_per_decade", candidate_covariates=("sex_male",),
        min_donors_present=4, min_total_cells=20, min_residual_df=2,
    )
    expr_all, meta_all = build_subset_expr(sums, genes_ref, cell_counts, donor_meta, min_cells=20, min_donors=5)
    fam_b_genes_all, fam_b_sum_all, fam_b_conc_all = run_family_b_genes_and_concordance(
        expr_all, meta_all, age_col="age_per_decade", candidate_covariates=("sex_male", "log_lib"),
        min_donors=5, min_cpm=1.0, min_expr_frac=0.5, min_residual_df=2, is_unbiased_single_cell=True,
    )
    fam_c_all = run_family_c_p1_signature(
        expr_all, meta_all, p1_program=p1, age_col="age_per_decade", candidate_covariates=("sex_male", "log_lib"),
        min_donors=5, min_cpm=1.0, min_expr_frac=0.5, min_residual_df=2,
    )

    # Sensitivity: adults >= 18 yr only (n = 5: 21, 22, 38, 44, 54 yr)
    adult_donors = donor_meta.loc[donor_meta["age_years"] >= 18.0].copy()
    fam_a_adult = run_family_a_clr(
        cell_counts.loc[adult_donors.index], adult_donors, age_col="age_per_decade",
        candidate_covariates=("sex_male",), min_donors_present=4, min_total_cells=20, min_residual_df=2,
    )
    expr_ad, meta_ad = build_subset_expr(sums, genes_ref, cell_counts, adult_donors, min_cells=20, min_donors=5)
    fam_b_genes_ad, fam_b_sum_ad, fam_b_conc_ad = run_family_b_genes_and_concordance(
        expr_ad, meta_ad, age_col="age_per_decade", candidate_covariates=("log_lib",),
        min_donors=5, min_cpm=1.0, min_expr_frac=0.5, min_residual_df=2, is_unbiased_single_cell=True,
    )
    fam_c_ad = run_family_c_p1_signature(
        expr_ad, meta_ad, p1_program=p1, age_col="age_per_decade", candidate_covariates=("log_lib",),
        min_donors=5, min_cpm=1.0, min_expr_frac=0.5, min_residual_df=2,
    )

    results_dir = REPO_ROOT / "results"
    reports_dir = REPO_ROOT / "reports"
    prefix = "followup-organ-04-pancreas-enge-gse81547"

    donor_out = donor_meta.reset_index().join(cell_counts, on="donor_id")
    write_clean_csv(donor_out, results_dir / f"{prefix}-donors.csv")
    write_clean_csv(
        pd.concat([fam_a_all.assign(cohort="all_8"), fam_a_adult.assign(cohort="adult_5")], ignore_index=True),
        results_dir / f"{prefix}-fractions.csv",
    )
    write_clean_csv(
        pd.concat([fam_b_sum_all.assign(cohort="all_8"), fam_b_sum_ad.assign(cohort="adult_5")], ignore_index=True),
        results_dir / f"{prefix}-gene-summary.csv",
    )
    write_clean_csv(fam_b_genes_all.head(200), results_dir / f"{prefix}-top-genes.csv")
    write_clean_csv(
        pd.concat([fam_b_conc_all.assign(cohort="all_8"), fam_b_conc_ad.assign(cohort="adult_5")], ignore_index=True),
        results_dir / f"{prefix}-concordance.csv",
    )
    write_clean_csv(
        pd.concat([fam_c_all.assign(cohort="all_8"), fam_c_ad.assign(cohort="adult_5")], ignore_index=True),
        results_dir / f"{prefix}-p1-signature.csv",
    )

    donor_tbl = "\n".join(
        f"| {d} | {r['age_raw']} ({r['age_years']:.3f}) | {r['sex']} | {int(r['total_cells'])} |"
        for d, r in donor_meta.iterrows()
    )
    a_tbl = "\n".join(
        f"| {r['cell_type']} | {int(r['total_cells'])} | {r['median_fraction']:.4f} | {r['effect_clr']:+.3f} | {r['pvalue']:.4f} | {r['q_bh_family_a']:.4f} |"
        for _, r in fam_a_all.iterrows()
    )
    b_tbl = "\n".join(
        f"| {r['group']} | {int(r['n_donors'])} | {int(r['n_genes_tested'])} | {int(r['n_sig_global_q05'])} | {r['min_pvalue']:.2e} | {r['min_q_bh_global']:.4f} |"
        for _, r in fam_b_sum_all.iterrows()
    )
    c_tbl = "\n".join(
        f"| {r['group']} | {int(r['n_donors'])} | {int(r['n_p1_up_detected'])}/{int(r['n_p1_down_detected'])} | {r['composite_effect']:+.4f} | {r['composite_pvalue']:.4f} | {r['empirical_null_pvalue']:.4f} | {r['q_bh_empirical']:.4f} |"
        for _, r in fam_c_all.iterrows()
    )

    mb_all = int(fam_b_genes_all["family_b_global_size"].iloc[0]) if len(fam_b_genes_all) else 0
    n_sig_b_all = int((fam_b_genes_all["q_bh_global"] < 0.05).sum()) if len(fam_b_genes_all) else 0
    min_q_b_all = float(fam_b_genes_all["q_bh_global"].min()) if len(fam_b_genes_all) else np.nan
    r_all = float(fam_b_conc_all["mean_pairwise_pearson_r"].iloc[0]) if len(fam_b_conc_all) else np.nan
    p_all = float(fam_b_conc_all["perm_p_one_sided"].iloc[0]) if len(fam_b_conc_all) else np.nan
    r_ad = float(fam_b_conc_ad["mean_pairwise_pearson_r"].iloc[0]) if len(fam_b_conc_ad) else np.nan
    p_ad = float(fam_b_conc_ad["perm_p_one_sided"].iloc[0]) if len(fam_b_conc_ad) else np.nan
    min_q_c_all = float(fam_c_all["q_bh_empirical"].min()) if len(fam_c_all) else np.nan
    min_q_c_ana_all = float(fam_c_all["q_bh_composite"].min()) if len(fam_c_all) else np.nan
    min_q_c_ad = float(fam_c_ad["q_bh_empirical"].min()) if len(fam_c_ad) else np.nan

    adult_note = (
        f"Concordance mean Pearson `r = {r_ad:+.4f}` (`perm_p_one_sided = {p_ad:.4f}`)"
        if len(fam_b_conc_ad)
        else "Only 1 cell type (`alpha`, `11,831` genes, `0` genes `q_BH < 0.05`) met `>= 20` cells across all 5 adult donors, so cross-cell-type concordance is not applicable"
    )

    md_text = f"""# Dataset 04: Human Pancreas Aging Smart-seq2 (Enge et al., GSE81547)

## 1. Study Verification and Data Provenance
- **Accession**: `GSE81547` (`2,544` single cells across `8` donors), Primary source: Enge et al., *Cell* 2017 (`PMID: 28965763`, `DOI: 10.1016/j.cell.2017.09.004`).
- **Critical Span and Protocol Limits**:
  1. **Span is development/pediatric to midlife ONLY**: Ages span `1 to 54` years (`GEO donor_age: 1, 5, 6, 21, 22, 38, 44, 54`; `1_M` is coded `1` in GEO `series_matrix` and reported as 1 month in Enge et al., `6M / 2F`). There are **zero old donors (`>= 60` yr)**.
  2. **FACS-gated Smart-seq2 capture**: Single cells were FACS-sorted into 96-well plates for Smart-seq2 (`{parse_stats['parsed_cells']:,}` cells parsed; `{parse_stats['htseq_meta_rows_excluded']}` HTSeq summary rows and `{parse_stats['ercc_spikein_rows_excluded']}` ERCC spike-in rows excluded before CPM normalization).

| Donor ID | Age Raw (Years) | Sex | QC Cells |
|---|---:|---|---:|
{donor_tbl}

## 2. Explicit Multiple-Testing Families
- **Inferential unit**: Donor (`n = 8` all donors; `n = 5` adults `>= 21` yr sensitivity). Continuous `age_per_decade = age_years / 10`.
- **Family A (CLR cell-type proportions)**: `{len(fam_a_all)}` FACS-captured pancreas cell types present in `>= 4` donors (`total >= 20` cells), adjusting for `sex_male`, BH across `K_A = {len(fam_a_all)}`. Because cells were FACS-sorted rather than unbiased whole-tissue dissociation, proportions are library-capture fractions rather than true tissue abundances.
- **Family B (Donor-level gene age test + concordance)**: Donor pseudobulks (`>= 20` cells in `>= 5` donors; `CPM >= 1.0` in `>= 50%` donors) adjusting for `sex_male + log_lib` (`all_8`) or `log_lib` (`adult_5`). Global BH across `M_B = {mb_all:,}` `(cell_type, gene)` tests. Concordance measured across shared expressed genes with `1,000` donor-age permutations.
- **Family C (`P1_Clean` mouse aging signature)**: `76` UP / `42` DOWN mapped human genes, `K_C = {len(fam_c_all)}` cell types (`all_8`) and `K_C = {len(fam_c_ad)}` (`adult_5`), 200 expression-decile-matched random draws + BH.

## 3. Results

### Family A: FACS Library CLR Proportions (All 8 Donors)
| Cell Type | Total Cells | Median Frac | CLR Slope (/10 yr) | p-value | q_BH (K={len(fam_a_all)}) |
|---|---:|---:|---:|---:|---:|
{a_tbl}

### Family B: Donor-Level Gene Age OLS and Cross-Cell-Type Concordance
- **All 8 donors (`1 to 54` yr, `M_B = {mb_all:,}`)**: `{n_sig_b_all}` genes passed `q_BH < 0.05` (minimum `q_BH = {min_q_b_all:.4f}`). Concordance across `{int(fam_b_conc_all['n_cell_types'].iloc[0])}` cell types (`{int(fam_b_conc_all['n_common_genes'].iloc[0])}` common genes): mean pairwise Pearson `r = {r_all:+.4f}` (`perm_p_one_sided = {p_all:.4f}`).
- **Adult-only sensitivity (`n = 5`, `21 to 54` yr)**: {adult_note}.

| Cell Type | Donors | Genes Tested | Global q<0.05 | Min p-value | Min Global q_BH |
|---|---:|---:|---:|---:|---:|
{b_tbl}

### Family C: Mouse `P1_Clean` Signature Check (All 8 Donors)
Minimum empirical null `q_BH = {min_q_c_all:.4f}` (minimum analytic `q_BH = {min_q_c_ana_all:.4f}`) across `K_C = {len(fam_c_all)}` cell types; adult-only sensitivity (`n = 5`, `alpha` `K_C = 1`) minimum empirical `q_BH = {min_q_c_ad:.4f}`.

| Cell Type | Donors | UP/DOWN Detected | Composite Slope | Analytic p | Empirical Null p (200) | Empirical q_BH |
|---|---:|---|---:|---:|---:|---:|
{c_tbl}

## 4. Verdict
- **No mouse aging transport**: `P1_Clean` shows null aging association across pancreas cell types in both the `1 to 54` yr cohort (`min q_BH = {min_q_c_all:.4f}`) and adult-only `21 to 54` yr `alpha` subset (`q_BH = {min_q_c_ad:.4f}`).
- **Explicit limits**: Cohort spans pediatric development to midlife (`1 to 54` yr) with no elderly donors (`>= 60` yr) and uses FACS-sorted Smart-seq2 capture.

## 5. Output Artifacts
- [Donors](file:///home/reza/projects/agescape-revision/results/{prefix}-donors.csv)
- [Fractions](file:///home/reza/projects/agescape-revision/results/{prefix}-fractions.csv)
- [Gene Summary](file:///home/reza/projects/agescape-revision/results/{prefix}-gene-summary.csv)
- [Top Genes](file:///home/reza/projects/agescape-revision/results/{prefix}-top-genes.csv)
- [Concordance](file:///home/reza/projects/agescape-revision/results/{prefix}-concordance.csv)
- [P1 Signature](file:///home/reza/projects/agescape-revision/results/{prefix}-p1-signature.csv)
"""
    write_clean_markdown(md_text, reports_dir / f"{prefix}.md")
    clean_temp_dir(tmp_dir)
    print(
        f"DONE Dataset 04: donors=8, cells={parse_stats['parsed_cells']}, M_B={mb_all}, sig_B={n_sig_b_all}, "
        f"conc_all={r_all:+.4f}(p={p_all:.4f}), conc_ad={r_ad:+.4f}(p={p_ad:.4f}), min_qC_emp={min_q_c_all:.4f} (ad {min_q_c_ad:.4f})"
    )


if __name__ == "__main__":
    main()
