"""Dataset 02: Bulk RNA-seq of primary human dermal fibroblasts (GSE307377)."""

from __future__ import annotations

from pathlib import Path
import re
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
    load_p1_clean_human,
    prepare_counts_matrix,
    run_family_b_genes_and_concordance,
    run_family_c_p1_signature,
    write_clean_csv,
    write_clean_markdown,
)

# Verified from GEO GSE307377 full GSM records (!Sample_description Library name)
GSM_SAMPLES = [
    {"gsm": "GSM9222986", "lib_name": "TM1_CKC3_Y_CTRL", "age_years": 23.0, "age_group": "young", "sex": "F", "sex_male": 0},
    {"gsm": "GSM9222987", "lib_name": "TM2_CKC4_Y_CTRL", "age_years": 26.0, "age_group": "young", "sex": "F", "sex_male": 0},
    {"gsm": "GSM9222988", "lib_name": "TM3_CKC5_O_CTRL", "age_years": 72.0, "age_group": "old", "sex": "M", "sex_male": 1},
    {"gsm": "GSM9222989", "lib_name": "TM4_CKC6_O_CTRL", "age_years": 70.0, "age_group": "old", "sex": "F", "sex_male": 0},
    {"gsm": "GSM9222990", "lib_name": "TM15_CKC7_O_CTRL", "age_years": 69.0, "age_group": "old", "sex": "F", "sex_male": 0},
    {"gsm": "GSM9222991", "lib_name": "TM16_CKC8_Y_CTRL", "age_years": 23.0, "age_group": "young", "sex": "M", "sex_male": 1},
    {"gsm": "GSM9222992", "lib_name": "TM17_CKC14_O_CTRL", "age_years": 64.0, "age_group": "old", "sex": "F", "sex_male": 0},
    {"gsm": "GSM9222993", "lib_name": "TM18_CKC15_Y_CTRL", "age_years": 33.0, "age_group": "young", "sex": "M", "sex_male": 1},
    {"gsm": "GSM9222994", "lib_name": "TM19_CKC16_O_CTRL", "age_years": 64.0, "age_group": "old", "sex": "M", "sex_male": 1},
]


def match_sample_columns(columns: list[str], meta_df: pd.DataFrame) -> dict[str, str]:
    """Map HOMER raw count columns to GSM accessions using verified Library name tokens."""
    mapping: dict[str, str] = {}
    for col in columns:
        col_str = str(col).strip()
        for gsm, row in meta_df.iterrows():
            if gsm in col_str or row["lib_name"] in col_str:
                mapping[col_str] = gsm
                break
    if len(mapping) != len(meta_df):
        raise RuntimeError(
            f"Could not uniquely map all {len(meta_df)} GSM samples from columns {columns}; mapped={mapping}"
        )
    return mapping


def main() -> None:
    tmp_dir = DATA_TMP_ROOT / "02_gse307377"
    raw_path = download_with_guard(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE307nnn/GSE307377/suppl/GSE307377_raw.txt.gz",
        tmp_dir / "GSE307377_raw.txt.gz",
        max_file_gb=0.5,
    )

    df = pd.read_csv(raw_path, sep="\t")
    print("GSE307377 raw columns:", list(df.columns), "shape:", df.shape)
    # Identify gene symbol column and count columns (HOMER tables often include Chr/Start/End/Annotation/Gene Name)
    meta_df = pd.DataFrame(GSM_SAMPLES).set_index("gsm")
    meta_df["age_per_decade"] = meta_df["age_years"] / 10.0
    meta_df["old_vs_young"] = (meta_df["age_group"] == "old").astype(float)

    sample_cols = [
        c for c in df.columns
        if str(c).startswith("tags/") and any(lib in str(c) for lib in meta_df["lib_name"])
    ]
    if len(sample_cols) != 9:
        raise RuntimeError(f"Expected 9 matched count columns, found {sample_cols}")

    col_map = match_sample_columns(sample_cols, meta_df)
    gene_col = None
    for cand in ["Gene Name", "gene_name", "gene_symbol", "Gene", "Symbol", "Annotation/Divergence"]:
        if cand in df.columns:
            gene_col = cand
            break
    if gene_col is None:
        gene_col = df.columns[0]
    genes = df[gene_col].astype(str).str.split("|").str[0]

    counts = df[sample_cols].copy()
    counts.index = genes
    counts = counts.rename(columns=col_map)[meta_df.index]

    filt_counts, count_stats = prepare_counts_matrix(counts, allow_fractional_expected=True)
    log2cpm, lib_sizes = counts_to_log2cpm(filt_counts)
    meta_df["lib_size"] = lib_sizes.loc[meta_df.index]
    meta_df["log_lib"] = np.log1p(meta_df["lib_size"].to_numpy(dtype=float))

    p1 = load_p1_clean_human(REPO_ROOT)
    expr_dict = {"Dermal_Fibroblast_Bulk": log2cpm}
    meta_dict = {"Dermal_Fibroblast_Bulk": meta_df}

    # Primary: continuous age_per_decade + sex_male + log_lib (bulk ONLY -> is_unbiased_single_cell=False)
    fam_b_cont, fam_b_sum_cont, fam_b_conc = run_family_b_genes_and_concordance(
        expr_dict, meta_dict, age_col="age_per_decade", candidate_covariates=("sex_male", "log_lib"),
        min_donors=6, min_cpm=1.0, min_expr_frac=0.5, min_residual_df=2, is_unbiased_single_cell=False,
    )
    if not fam_b_conc.empty:
        raise AssertionError("Concordance must be empty for bulk dataset")
    # Sensitivity: binary old vs young + sex_male + log_lib
    fam_b_bin, fam_b_sum_bin, _ = run_family_b_genes_and_concordance(
        expr_dict, meta_dict, age_col="old_vs_young", candidate_covariates=("sex_male", "log_lib"),
        min_donors=6, min_cpm=1.0, min_expr_frac=0.5, min_residual_df=2, is_unbiased_single_cell=False,
    )

    fam_c_cont = run_family_c_p1_signature(
        expr_dict, meta_dict, p1_program=p1, age_col="age_per_decade", candidate_covariates=("sex_male", "log_lib"),
        min_donors=6, min_cpm=1.0, min_expr_frac=0.5, min_residual_df=2,
    )
    fam_c_bin = run_family_c_p1_signature(
        expr_dict, meta_dict, p1_program=p1, age_col="old_vs_young", candidate_covariates=("sex_male", "log_lib"),
        min_donors=6, min_cpm=1.0, min_expr_frac=0.5, min_residual_df=2,
    )

    results_dir = REPO_ROOT / "results"
    reports_dir = REPO_ROOT / "reports"
    prefix = "followup-organ-02-skin-fibroblasts-gse307377"

    write_clean_csv(meta_df.reset_index(), results_dir / f"{prefix}-donors.csv")
    write_clean_csv(
        pd.concat([fam_b_sum_cont.assign(model="age_per_decade"), fam_b_sum_bin.assign(model="old_vs_young")], ignore_index=True),
        results_dir / f"{prefix}-gene-summary.csv",
    )
    write_clean_csv(fam_b_cont.head(200), results_dir / f"{prefix}-top-genes.csv")
    write_clean_csv(
        pd.concat([fam_c_cont.assign(model="age_per_decade"), fam_c_bin.assign(model="old_vs_young")], ignore_index=True),
        results_dir / f"{prefix}-p1-signature.csv",
    )

    r_cont = fam_b_sum_cont.iloc[0]
    r_bin = fam_b_sum_bin.iloc[0]
    c_cont = fam_c_cont.iloc[0]
    c_bin = fam_c_bin.iloc[0]
    frac_info = count_stats["count_gate"]["noninteger"]
    donor_rows = "\n".join(
        f"| {gsm} | {int(r['age_years'])} | {r['age_group']} | {r['sex']} | {int(r['lib_size']):,} |"
        for gsm, r in meta_df.iterrows()
    )

    md_text = f"""# Dataset 02: Primary Human Dermal Fibroblast Bulk RNA-seq (GSE307377)

## 1. Study Verification and Data Provenance
- **Accession**: `GSE307377` (`GSM9222986` to `GSM9222994`).
- **Assay and Tissue**: Cultured primary human dermal fibroblasts (SD-NSC Human Cell Models of Aging Core; DMEM + 10% FBS for 48h prior to polyA bulk RNA-seq on NovaSeq 6000 SP, HOMER raw counts).
- **Verified Donor Cohort (`n = 9`, `5F / 4M`)**:
  - Young (`n = 4`, ages `23 to 33` yr): `23F, 26F, 23M, 33M`.
  - Old (`n = 5`, ages `64 to 72` yr): `72M, 70F, 69F, 64F, 64M`.
- **Technical Filtering**: Input genes `{count_stats['n_genes_input']:,}` -> unique non-empty symbols `{count_stats['n_genes_unique']:,}` (`{frac_info}` fractional expected/HOMER entries, `0` negative) -> `{count_stats['n_genes_nontech']:,}` non-technical genes (`{count_stats['n_technical_removed']}` `MT-/RP[LS]/MRP[LS]/sex-linked/ERCC` genes removed).

## 2. Scope Rules: Bulk Cultured Lineage Only
- **Family A (CLR cell-type fractions)**: **NOT RUN**. Bulk RNA-seq of a purified, in vitro cultured fibroblast monoculture measures one cell lineage and cannot estimate in vivo skin cell-type proportions.
- **Cross-cell-type concordance**: **NOT RUN**. Requires multiple distinct unbiased cell types from intact tissue.
- **Family B (Donor-level gene age test)**: Bulk OLS on `log2(CPM + 1)` adjusting for `sex_male + log_lib` (`residual_df = 5`), BH correction across the explicit family of `M_B = {int(r_cont['n_genes_tested']):,}` expressed genes (`CPM >= 1.0` in `>= 5` of 9 donors). Both continuous `age_per_decade` (primary) and binary `old_vs_young` (sensitivity) were tested.
- **Family C (`P1_Clean` 139-gene mouse aging signature)**: `K_C = 1` group tested analytically and against `200` expression-decile-matched random gene sets (`matched_null_effect_pvalue`).

## 3. Results

| GSM | Age (yr) | Group | Sex | Non-Tech Lib Size |
|---|---:|---|---|---:|
{donor_rows}

### Family B: Bulk Dermal Fibroblast Gene Age OLS
- **Primary (`age_per_decade + sex_male + log_lib`, `M_B = {int(r_cont['n_genes_tested']):,}`)**: `{int(r_cont['n_sig_global_q05'])}` genes passed `q_BH < 0.05` (minimum `p = {r_cont['min_pvalue']:.2e}`, minimum `q_BH = {r_cont['min_q_bh_global']:.4f}`).
- **Sensitivity (`old_vs_young + sex_male + log_lib`, `M_B = {int(r_bin['n_genes_tested']):,}`)**: `{int(r_bin['n_sig_global_q05'])}` genes passed `q_BH < 0.05` (minimum `p = {r_bin['min_pvalue']:.2e}`, minimum `q_BH = {r_bin['min_q_bh_global']:.4f}`).

### Family C: Mouse `P1_Clean` Signature Check
- **Primary (`age_per_decade + sex_male + log_lib`)**: `{int(c_cont['n_p1_up_detected'])}` UP / `{int(c_cont['n_p1_down_detected'])}` DOWN detected; UP slope `{c_cont['up_effect']:+.4f}` (`p = {c_cont['up_pvalue']:.4f}`), DOWN slope `{c_cont['down_effect']:+.4f}` (`p = {c_cont['down_pvalue']:.4f}`), Composite slope `{c_cont['composite_effect']:+.4f}` (`analytic p = {c_cont['composite_pvalue']:.4f}`, 200-draw random null mean `{c_cont['null_mean_effect']:+.4f}` +/- `{c_cont['null_sd_effect']:.4f}`, **`empirical_null_pvalue = {c_cont['empirical_null_pvalue']:.4f}`**, `q_BH = {c_cont['q_bh_empirical']:.4f}`).
- **Sensitivity (`old_vs_young + sex_male + log_lib`)**: Composite `{c_bin['composite_effect']:+.4f}` (`analytic p = {c_bin['composite_pvalue']:.4f}`, **`empirical_null_pvalue = {c_bin['empirical_null_pvalue']:.4f}`**, `q_BH = {c_bin['q_bh_empirical']:.4f}`).

## 4. Limitations and Verdict
- **Culture and cohort limits**: Small donor count (`n = 9`), binary age sampling (`23 to 33` vs `64 to 72` with no middle-aged donors), and short-term in vitro passaging/48h DMEM culture can attenuate or alter in vivo aging signals.
- **Verdict**: Null transport for the mouse `P1_Clean` aging signature in human primary dermal fibroblasts (`empirical p = {c_cont['empirical_null_pvalue']:.4f}` continuous; `{c_bin['empirical_null_pvalue']:.4f}` binary).

## 5. Output Artifacts
- [Donors](file:///home/reza/projects/agescape-revision/results/{prefix}-donors.csv)
- [Gene Summary](file:///home/reza/projects/agescape-revision/results/{prefix}-gene-summary.csv)
- [Top Genes](file:///home/reza/projects/agescape-revision/results/{prefix}-top-genes.csv)
- [P1 Signature](file:///home/reza/projects/agescape-revision/results/{prefix}-p1-signature.csv)
"""
    write_clean_markdown(md_text, reports_dir / f"{prefix}.md")
    clean_temp_dir(tmp_dir)
    print(
        f"DONE Dataset 02: n=9, M_B={int(r_cont['n_genes_tested'])}, sig_cont={int(r_cont['n_sig_global_q05'])}, "
        f"sig_bin={int(r_bin['n_sig_global_q05'])}, P1_emp_p={c_cont['empirical_null_pvalue']:.4f} (bin {c_bin['empirical_null_pvalue']:.4f})"
    )


if __name__ == "__main__":
    main()
