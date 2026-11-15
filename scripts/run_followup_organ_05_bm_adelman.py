"""Dataset 05: Adelman et al. (2019) bulk FACS-sorted human BM HSCe RNA-seq (GSE104406)."""

from __future__ import annotations

import gzip
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
    run_family_b_genes_and_concordance,
    run_family_c_p1_signature,
    write_clean_csv,
    write_clean_markdown,
)

# Verified verbatim from GEO GSE104406 GSM2797287..GSM2797306
ADELMAN_GSMS = [
    {"gsm": "GSM2797287", "title": "Young_1_RNA", "age_raw": "18", "exact_age": 18.0, "age_group": "young", "sex": "M", "sex_male": 1},
    {"gsm": "GSM2797288", "title": "Young_12_RNA", "age_raw": "25", "exact_age": 25.0, "age_group": "young", "sex": "M", "sex_male": 1},
    {"gsm": "GSM2797289", "title": "Young_16_RNA", "age_raw": "27", "exact_age": 27.0, "age_group": "young", "sex": "F", "sex_male": 0},
    {"gsm": "GSM2797290", "title": "Young_17_RNA", "age_raw": "28", "exact_age": 28.0, "age_group": "young", "sex": "M", "sex_male": 1},
    {"gsm": "GSM2797291", "title": "Young_25_RNA", "age_raw": "<30", "exact_age": np.nan, "age_group": "young", "sex": "M", "sex_male": 1},
    {"gsm": "GSM2797292", "title": "Young_4_RNA", "age_raw": "21", "exact_age": 21.0, "age_group": "young", "sex": "F", "sex_male": 0},
    {"gsm": "GSM2797293", "title": "Young_5_RNA", "age_raw": "22", "exact_age": 22.0, "age_group": "young", "sex": "M", "sex_male": 1},
    {"gsm": "GSM2797294", "title": "Young_6_RNA", "age_raw": "22", "exact_age": 22.0, "age_group": "young", "sex": "M", "sex_male": 1},
    {"gsm": "GSM2797295", "title": "Young_7_RNA", "age_raw": "22", "exact_age": 22.0, "age_group": "young", "sex": "F", "sex_male": 0},
    {"gsm": "GSM2797296", "title": "Young_9_RNA", "age_raw": "22", "exact_age": 22.0, "age_group": "young", "sex": "M", "sex_male": 1},
    {"gsm": "GSM2797297", "title": "Old_11_RNA", "age_raw": "67", "exact_age": 67.0, "age_group": "old", "sex": "F", "sex_male": 0},
    {"gsm": "GSM2797298", "title": "Old_18_RNA", "age_raw": "70", "exact_age": 70.0, "age_group": "old", "sex": "F", "sex_male": 0},
    {"gsm": "GSM2797299", "title": "Old_19_RNA", "age_raw": "70", "exact_age": 70.0, "age_group": "old", "sex": "F", "sex_male": 0},
    {"gsm": "GSM2797300", "title": "Old_21_RNA", "age_raw": "71", "exact_age": 71.0, "age_group": "old", "sex": "M", "sex_male": 1},
    {"gsm": "GSM2797301", "title": "Old_22_RNA", "age_raw": "71", "exact_age": 71.0, "age_group": "old", "sex": "F", "sex_male": 0},
    {"gsm": "GSM2797302", "title": "Old_23_RNA", "age_raw": "72", "exact_age": 72.0, "age_group": "old", "sex": "F", "sex_male": 0},
    {"gsm": "GSM2797303", "title": "Old_25_RNA", "age_raw": "75", "exact_age": 75.0, "age_group": "old", "sex": "F", "sex_male": 0},
    {"gsm": "GSM2797304", "title": "Old_28_RNA", "age_raw": "74", "exact_age": 74.0, "age_group": "old", "sex": "F", "sex_male": 0},
    {"gsm": "GSM2797305", "title": "Old_20_RNA", "age_raw": "75", "exact_age": 75.0, "age_group": "old", "sex": "F", "sex_male": 0},
    {"gsm": "GSM2797306", "title": "Old_4_RNA", "age_raw": "65", "exact_age": 65.0, "age_group": "old", "sex": "F", "sex_male": 0},
]

HTSEQ_BAD = {"no_feature", "ambiguous", "too_low_aQual", "not_aligned", "alignment_not_unique"}


def load_ensembl_to_symbol_map() -> dict[str, str]:
    """Build Ensembl ID -> human gene symbol dictionary from reference ortholog table."""
    ortho = pd.read_csv(REPO_ROOT / "data/reference/orthologs-mouse-human.csv")
    mapping = {}
    for _, r in ortho.dropna(subset=["human_ensembl_id", "human_gene_symbol"]).iterrows():
        eid = str(r["human_ensembl_id"]).split(".")[0].strip()
        sym = str(r["human_gene_symbol"]).strip()
        if eid and sym:
            mapping[eid] = sym
    return mapping


def parse_raw_tar(tar_path: Path, gsm_ids: list[str]) -> tuple[pd.DataFrame, str]:
    """Read 20 per-GSM count files from GSE104406_RAW.tar into genes x GSM raw matrix."""
    series_by_gsm: dict[str, pd.Series] = {}
    with tarfile.open(tar_path, "r:") as tf:
        for member in tf:
            if not member.isfile():
                continue
            fname = Path(member.name).name
            m = re.match(r"^(GSM\d+)_", fname)
            if not m or m.group(1) not in gsm_ids:
                continue
            gsm = m.group(1)
            raw_bytes = tf.extractfile(member).read()
            text = gzip.decompress(raw_bytes).decode("utf-8", errors="replace") if fname.endswith(".gz") else raw_bytes.decode("utf-8", errors="replace")
            genes = []
            vals = []
            for ln in text.splitlines():
                if not ln.strip():
                    continue
                parts = ln.split("\t")
                if len(parts) < 2:
                    continue
                try:
                    val = float(parts[-1])
                except ValueError:
                    continue
                g = parts[0].strip().strip('"')
                if g.startswith("__") or g in HTSEQ_BAD:
                    continue
                genes.append(g)
                vals.append(val)
            series_by_gsm[gsm] = pd.Series(vals, index=pd.Index(genes, dtype=str))

    if len(series_by_gsm) != len(gsm_ids):
        raise RuntimeError(f"Expected {len(gsm_ids)} GSM count files, found {len(series_by_gsm)}")

    raw_df = pd.DataFrame(series_by_gsm)[gsm_ids]
    id_style = "ensembl" if raw_df.index.astype(str).str.startswith("ENSG").mean() > 0.5 else "symbol"
    if id_style == "ensembl":
        ens_map = load_ensembl_to_symbol_map()
        raw_df = raw_df.copy()
        clean_eids = raw_df.index.astype(str).str.split(".").str[0]
        raw_df.index = pd.Index([ens_map.get(e, e) for e in clean_eids], dtype=str)
    return raw_df, id_style


def main() -> None:
    tmp_dir = DATA_TMP_ROOT / "05_gse104406"
    tar_path = download_with_guard(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104406/suppl/GSE104406_RAW.tar",
        tmp_dir / "GSE104406_RAW.tar",
        max_file_gb=0.5,
    )

    meta_df = pd.DataFrame(ADELMAN_GSMS).set_index("gsm")
    raw_df, id_style = parse_raw_tar(tar_path, meta_df.index.tolist())
    filt_counts, c_stats = prepare_counts_matrix(raw_df, allow_fractional_expected=True)
    log2cpm, lib_sizes = counts_to_log2cpm(filt_counts)

    meta_df["lib_size"] = lib_sizes.loc[meta_df.index]
    meta_df["log_lib"] = np.log1p(meta_df["lib_size"].to_numpy(dtype=float))
    meta_df["old_vs_young"] = (meta_df["age_group"] == "old").astype(float)
    meta_df["age_20_imputed25"] = meta_df["exact_age"].fillna(25.0) / 10.0

    p1 = load_p1_clean_human(REPO_ROOT)

    # 1. Primary (n=19 exact ages)
    m19 = meta_df.loc[meta_df["exact_age"].notna()].copy()
    m19["age_per_decade"] = m19["exact_age"] / 10.0
    e19 = log2cpm[m19.index]
    b19, bs19, conc19 = run_family_b_genes_and_concordance(
        {"HSCe_Sorted_Bulk": e19}, {"HSCe_Sorted_Bulk": m19},
        age_col="age_per_decade", candidate_covariates=("sex_male", "log_lib"),
        min_donors=10, min_cpm=1.0, min_expr_frac=0.5, min_residual_df=2, is_unbiased_single_cell=False,
    )
    if not conc19.empty:
        raise AssertionError("Concordance must be empty for sorted bulk HSCe dataset")
    c19 = run_family_c_p1_signature(
        {"HSCe_Sorted_Bulk": e19}, {"HSCe_Sorted_Bulk": m19},
        p1_program=p1, age_col="age_per_decade", candidate_covariates=("sex_male", "log_lib"),
        min_donors=10, min_cpm=1.0, min_expr_frac=0.5, min_residual_df=2,
    )

    # 2. Sensitivity (n=20, Young_25_RNA coded as 25)
    b20_cont, bs20_cont, _ = run_family_b_genes_and_concordance(
        {"HSCe_Sorted_Bulk": log2cpm}, {"HSCe_Sorted_Bulk": meta_df},
        age_col="age_20_imputed25", candidate_covariates=("sex_male", "log_lib"),
        min_donors=10, min_cpm=1.0, min_expr_frac=0.5, min_residual_df=2, is_unbiased_single_cell=False,
    )
    c20_cont = run_family_c_p1_signature(
        {"HSCe_Sorted_Bulk": log2cpm}, {"HSCe_Sorted_Bulk": meta_df},
        p1_program=p1, age_col="age_20_imputed25", candidate_covariates=("sex_male", "log_lib"),
        min_donors=10, min_cpm=1.0, min_expr_frac=0.5, min_residual_df=2,
    )

    # 3. Sensitivity (n=20 binary Aged vs Young)
    b20_bin, bs20_bin, _ = run_family_b_genes_and_concordance(
        {"HSCe_Sorted_Bulk": log2cpm}, {"HSCe_Sorted_Bulk": meta_df},
        age_col="old_vs_young", candidate_covariates=("sex_male", "log_lib"),
        min_donors=10, min_cpm=1.0, min_expr_frac=0.5, min_residual_df=2, is_unbiased_single_cell=False,
    )
    c20_bin = run_family_c_p1_signature(
        {"HSCe_Sorted_Bulk": log2cpm}, {"HSCe_Sorted_Bulk": meta_df},
        p1_program=p1, age_col="old_vs_young", candidate_covariates=("sex_male", "log_lib"),
        min_donors=10, min_cpm=1.0, min_expr_frac=0.5, min_residual_df=2,
    )

    results_dir = REPO_ROOT / "results"
    reports_dir = REPO_ROOT / "reports"
    prefix = "followup-organ-05-bm-adelman-gse104406"

    write_clean_csv(meta_df.reset_index(), results_dir / f"{prefix}-donors.csv")
    write_clean_csv(
        pd.concat(
            [
                bs19.assign(model="exact_age_n19"),
                bs20_cont.assign(model="imputed25_n20"),
                bs20_bin.assign(model="binary_old_vs_young_n20"),
            ],
            ignore_index=True,
        ),
        results_dir / f"{prefix}-gene-summary.csv",
    )
    write_clean_csv(b19.head(200), results_dir / f"{prefix}-top-genes.csv")
    write_clean_csv(
        pd.concat(
            [
                c19.assign(model="exact_age_n19"),
                c20_cont.assign(model="imputed25_n20"),
                c20_bin.assign(model="binary_old_vs_young_n20"),
            ],
            ignore_index=True,
        ),
        results_dir / f"{prefix}-p1-signature.csv",
    )

    s19 = bs19.iloc[0]
    s20c = bs20_cont.iloc[0]
    s20b = bs20_bin.iloc[0]
    p19 = c19.iloc[0]
    p20c = c20_cont.iloc[0]
    p20b = c20_bin.iloc[0]
    donor_tbl = "\n".join(
        f"| {gsm} | {r['title']} | {r['age_raw']} | {r['age_group']} | {r['sex']} | {int(r['lib_size']):,} |"
        for gsm, r in meta_df.iterrows()
    )

    md_text = f"""# Dataset 05: Human Bone Marrow Sorted HSCe Bulk RNA-seq (Adelman et al., GSE104406)

## 1. Study Verification and Data Provenance
- **Accession**: `GSE104406` (`GSM2797287` to `GSM2797306`), Primary source: Adelman et al., *Cancer Discovery* 2019 (`PMID: 31085557`).
- **Assay and Cell Population**: Bulk RNA-seq on FACS-sorted human bone marrow hematopoietic stem-cell-enriched populations (`Lin- CD34+ CD38-`), `n = 20` healthy donors without known hematological malignancy.
- **Age and Sex Breakdown**:
  - Young (`n = 10`, ages `18 to 28` plus `Young_25_RNA` labeled `<30`): `7M / 3F`.
  - Aged (`n = 10`, ages `65 to 75`): `1M / 9F`.
  - **Confounding note**: Strong age-by-sex imbalance (`Young` 70% male vs `Aged` 10% male), requiring explicit `sex_male` adjustment alongside `log_lib`.
- **Feature ID normalization**: Raw tar count ID style `{id_style}`; `{c_stats['n_genes_input']:,}` rows -> `{c_stats['n_genes_unique']:,}` unique symbols -> `{c_stats['n_genes_nontech']:,}` non-technical genes (`{c_stats['n_technical_removed']}` removed).

| GSM | Title | Age Raw | Group | Sex | Non-Tech Lib Size |
|---|---|---|---|---|---:|
{donor_tbl}

## 2. Scope Rules: FACS-Sorted Bulk HSCe Only
- **Family A (CLR cell-type fractions)**: **NOT RUN**. Cells were FACS-sorted (`Lin- CD34+ CD38-`) and sequenced in bulk, so proportions of bone marrow cell types cannot be inferred.
- **Cross-cell-type concordance**: **NOT RUN**. Single sorted population in bulk RNA-seq.
- **Family B (Donor-level gene age OLS)**: Adjusting for `sex_male + log_lib` (`CPM >= 1.0` in `>= 50%` of donors):
  - **Primary (`n = 19` exact ages)**: `M_B = {int(s19['n_genes_tested']):,}` genes (`residual_df = {int(s19['residual_df'])}`).
  - **Sensitivity 1 (`n = 20`, `Young_25_RNA` coded `25`)**: `M_B = {int(s20c['n_genes_tested']):,}` genes.
  - **Sensitivity 2 (`n = 20`, binary `old_vs_young`)**: `M_B = {int(s20b['n_genes_tested']):,}` genes.
- **Family C (`P1_Clean` signature)**: `K_C = 1`, 200 expression-decile-matched random draws + analytic OLS.

## 3. Results

### Family B: Bulk Sorted HSCe Gene Age OLS
- **Primary (`n = 19` exact ages)**: **`{int(s19['n_sig_global_q05'])}`** genes passed `q_BH < 0.05` (minimum `p = {s19['min_pvalue']:.2e}`, minimum `q_BH = {s19['min_q_bh_global']:.4f}`).
- **Sensitivity (`n = 20` continuous, `Young_25_RNA = 25`)**: `{int(s20c['n_sig_global_q05'])}` genes passed `q_BH < 0.05` (min `q_BH = {s20c['min_q_bh_global']:.4f}`).
- **Sensitivity (`n = 20` binary `old_vs_young`)**: `{int(s20b['n_sig_global_q05'])}` genes passed `q_BH < 0.05` (min `q_BH = {s20b['min_q_bh_global']:.4f}`).

### Family C: Mouse `P1_Clean` Signature Check
- **Primary (`n = 19` exact ages)**: `{int(p19['n_p1_up_detected'])}` UP / `{int(p19['n_p1_down_detected'])}` DOWN detected; UP slope `{p19['up_effect']:+.4f}` (`p = {p19['up_pvalue']:.4f}`), DOWN slope `{p19['down_effect']:+.4f}` (`p = {p19['down_pvalue']:.4f}`), Composite `{p19['composite_effect']:+.4f}` (`analytic p = {p19['composite_pvalue']:.4f}`, 200-draw random null mean `{p19['null_mean_effect']:+.4f}` +/- `{p19['null_sd_effect']:.4f}`, **`empirical_null_pvalue = {p19['empirical_null_pvalue']:.4f}`**, `q_BH = {p19['q_bh_empirical']:.4f}`).
- **Sensitivity (`n = 20` continuous)**: Composite `{p20c['composite_effect']:+.4f}` (`analytic p = {p20c['composite_pvalue']:.4f}`, **`empirical_null_pvalue = {p20c['empirical_null_pvalue']:.4f}`**).
- **Sensitivity (`n = 20` binary)**: Composite `{p20b['composite_effect']:+.4f}` (`analytic p = {p20b['composite_pvalue']:.4f}`, **`empirical_null_pvalue = {p20b['empirical_null_pvalue']:.4f}`**).

## 4. Limitations and Verdict
- **Limits**: FACS-sorted (`Lin- CD34+ CD38-`) bulk RNA-seq cannot deconvolve HSC subpopulations or evaluate tissue composition, and sex is strongly imbalanced across age (`7M/3F` young vs `1M/9F` aged).
- **Verdict**: Null transport for the mouse `P1_Clean` aging signature in sorted human bone marrow HSCe (`empirical p = {p19['empirical_null_pvalue']:.4f}` primary exact-age `n = 19`; `{p20b['empirical_null_pvalue']:.4f}` binary `n = 20`).

## 5. Output Artifacts
- [Donors](file:///home/reza/projects/agescape-revision/results/{prefix}-donors.csv)
- [Gene Summary](file:///home/reza/projects/agescape-revision/results/{prefix}-gene-summary.csv)
- [Top Genes](file:///home/reza/projects/agescape-revision/results/{prefix}-top-genes.csv)
- [P1 Signature](file:///home/reza/projects/agescape-revision/results/{prefix}-p1-signature.csv)
"""
    write_clean_markdown(md_text, reports_dir / f"{prefix}.md")
    clean_temp_dir(tmp_dir)
    print(
        f"DONE Dataset 05: id_style={id_style}, M_B={int(s19['n_genes_tested'])}, sig_19={int(s19['n_sig_global_q05'])}, "
        f"P1_emp_p19={p19['empirical_null_pvalue']:.4f} (bin20 {p20b['empirical_null_pvalue']:.4f})"
    )


if __name__ == "__main__":
    main()
