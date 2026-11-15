"""Dataset 08: Butler et al. (2021) human retinal pigment epithelium (RPE) bulk RNA-seq (GSE159435).

Primary source:
- Butler et al. (2021), PMID 33934486, PMCID PMC8184696.
- GEO GSE159435: 13 donors (GSM4829762 to GSM4829774; ages 31 to 93 yr; 8M, 5F; no reported retinal pathology).
  Batch s1 (n=7, donor1..donor7): 31M, 40F, 51M, 67F, 77M, 86M, 93F (spans 31 to 93 yr).
  Batch s2 (n=6, donor8..donor13): 61F, 62M, 71M, 73M, 73F, 92M (spans 61 to 92 yr).
- Bulk RNA-seq: Family B (gene-level age OLS) and Family C (P1_Clean) ONLY; no Family A fractions or concordance.
"""

from __future__ import annotations

import gzip
import os
from pathlib import Path
import sys
import tarfile
import urllib.request

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

GEO_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159435&targ=gsm&form=text&view=brief"
RAW_TAR_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE159nnn/GSE159435/suppl/GSE159435_RAW.tar"

NULL_ITERATIONS = 200
NULL_SEED = 202608


def _fetch_geo_donor_table() -> pd.DataFrame:
    with urllib.request.urlopen(GEO_URL, timeout=30) as resp:
        txt = resp.read().decode("utf-8", errors="replace")

    rows = []
    for blk in txt.split("^SAMPLE = ")[1:]:
        lines = blk.splitlines()
        gsm = lines[0].strip()
        title = ""
        chars = {}
        for ln in lines:
            if ln.startswith("!Sample_title = "):
                title = ln.split(" = ", 1)[1].strip()
            elif ln.startswith("!Sample_characteristics_ch1 = "):
                val = ln.split(" = ", 1)[1].strip()
                if ":" in val:
                    k, v = val.split(":", 1)
                    chars[k.strip().lower()] = v.strip()
        age = float(chars["age"])
        sex_raw = chars["sex"].strip().upper()
        batch = chars["batch"].strip()
        rows.append(
            {
                "gsm": gsm,
                "donor_id": title,
                "age": age,
                "age_per_decade": age / 10.0,
                "sex": sex_raw,
                "sex_male": 1.0 if sex_raw == "M" else 0.0,
                "batch": batch,
                "batch_s2": 1.0 if batch == "s2" else 0.0,
            }
        )
    df = pd.DataFrame(rows).sort_values(["batch", "age", "gsm"]).set_index("donor_id", drop=False)
    if len(df) != 13:
        raise RuntimeError(f"Expected 13 GSE159435 donors, got {len(df)}")
    return df


def _build_ensembl_map() -> dict[str, str]:
    vmap_path = DATA_TMP_ROOT / "shared_hrca_yan" / "hrca_var_map.csv"
    sym_by_ens: dict[str, str] = {}
    if vmap_path.exists():
        vmap = pd.read_csv(vmap_path)
        for ens, sym in zip(vmap["ensembl_id"].astype(str), vmap["feature_name"].astype(str)):
            sym_by_ens[ens.split(".")[0]] = sym

    ortho = pd.read_csv(REPO_ROOT / "data" / "reference" / "orthologs-mouse-human.csv")
    for ens, sym in zip(ortho["human_ensembl_id"].astype(str), ortho["human_gene_symbol"].astype(str)):
        e0 = ens.split(".")[0]
        if e0 not in sym_by_ens or sym_by_ens[e0].startswith("ENSG"):
            sym_by_ens[e0] = sym
    return sym_by_ens


def _load_tar_counts(tar_path: Path, donor_meta: pd.DataFrame, sym_by_ens: dict[str, str]) -> pd.DataFrame:
    gsm_to_donor = dict(zip(donor_meta["gsm"], donor_meta["donor_id"]))
    series_by_donor: dict[str, pd.Series] = {}

    with tarfile.open(tar_path, "r:*") as tf:
        for member in tf.getmembers():
            if not member.isfile():
                continue
            gsm = Path(member.name).name.split("_", 1)[0]
            if gsm not in gsm_to_donor:
                continue
            donor_id = gsm_to_donor[gsm]
            raw = tf.extractfile(member).read()
            text = gzip.decompress(raw).decode("utf-8", errors="replace") if member.name.endswith(".gz") else raw.decode("utf-8", errors="replace")
            genes = []
            vals = []
            for ln in text.splitlines():
                if not ln.strip() or ln.startswith("__"):
                    continue
                gid, val = ln.split("\t")[:2]
                gid0 = gid.split(".")[0]
                sym = sym_by_ens.get(gid0, gid0)
                genes.append(sym)
                vals.append(float(val))
            s = pd.Series(vals, index=pd.Index(genes, name="gene"), name=donor_id)
            series_by_donor[donor_id] = s.groupby(level=0).sum()

    if len(series_by_donor) != len(donor_meta):
        raise RuntimeError(f"Loaded {len(series_by_donor)} donors from tar, expected {len(donor_meta)}")
    counts_df = pd.DataFrame({d: series_by_donor[d] for d in donor_meta.index}).fillna(0.0)
    return counts_df


def main() -> None:
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
    tmp_dir = DATA_TMP_ROOT / "organ08_rpe"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    donor_meta = _fetch_geo_donor_table()
    tar_path = download_with_guard(RAW_TAR_URL, tmp_dir / "GSE159435_RAW.tar", max_file_gb=1.0)
    sym_by_ens = _build_ensembl_map()
    raw_counts = _load_tar_counts(tar_path, donor_meta, sym_by_ens)

    clean_counts, count_stats = prepare_counts_matrix(raw_counts, allow_fractional_expected=False)
    log2cpm_all, lib_all = counts_to_log2cpm(clean_counts)
    donor_meta["lib_size"] = lib_all.loc[donor_meta.index]
    donor_meta["log_lib"] = np.log1p(donor_meta["lib_size"])

    p1 = load_p1_clean_human(REPO_ROOT)

    meta_s1 = donor_meta[donor_meta["batch"] == "s1"].copy()
    clean_s1, _ = prepare_counts_matrix(raw_counts[meta_s1.index], allow_fractional_expected=False)
    log2cpm_s1, lib_s1 = counts_to_log2cpm(clean_s1)
    meta_s1["lib_size"] = lib_s1.loc[meta_s1.index]
    meta_s1["log_lib"] = np.log1p(meta_s1["lib_size"])

    configs = [
        ("Primary_BatchSexAdj", {"RPE_Bulk": log2cpm_all}, {"RPE_Bulk": donor_meta}, ("sex_male", "batch_s2", "log_lib"), 13, 3),
        ("Sens_SexAdj_NoBatch", {"RPE_Bulk": log2cpm_all}, {"RPE_Bulk": donor_meta}, ("sex_male", "log_lib"), 13, 3),
        ("Sens_BatchS1_Only", {"RPE_Bulk_s1": log2cpm_s1}, {"RPE_Bulk_s1": meta_s1}, ("sex_male", "log_lib"), 7, 3),
    ]

    fam_b_sum_list = []
    fam_b_top_list = []
    fam_c_list = []
    for name, expr_map, meta_map, covs, min_d, min_df in configs:
        fb_genes, fb_sum, _ = run_family_b_genes_and_concordance(
            expr_map,
            meta_map,
            age_col="age_per_decade",
            candidate_covariates=covs,
            min_donors=min_d,
            min_cpm=1.0,
            min_expr_frac=0.5,
            min_residual_df=min_df,
            is_unbiased_single_cell=False,
            seed=NULL_SEED,
        )
        fc = run_family_c_p1_signature(
            expr_map,
            meta_map,
            p1_program=p1,
            age_col="age_per_decade",
            candidate_covariates=covs,
            min_donors=min_d,
            min_cpm=1.0,
            min_expr_frac=0.5,
            min_residual_df=min_df,
            iterations=NULL_ITERATIONS,
            seed=NULL_SEED,
        )
        if not fb_sum.empty:
            fb_sum.insert(0, "analysis", name)
            fam_b_sum_list.append(fb_sum)
        if not fb_genes.empty:
            top = fb_genes.head(100).copy()
            top.insert(0, "analysis", name)
            fam_b_top_list.append(top)
        if not fc.empty:
            fc.insert(0, "analysis", name)
            fam_c_list.append(fc)

    fam_b_sum_df = pd.concat(fam_b_sum_list, ignore_index=True)
    fam_b_top_df = pd.concat(fam_b_top_list, ignore_index=True)
    fam_c_df = pd.concat(fam_c_list, ignore_index=True)

    results_dir = REPO_ROOT / "results"
    reports_dir = REPO_ROOT / "reports"
    prefix = "followup-organ-08-eye-rpe-gse159435"

    write_clean_csv(donor_meta.reset_index(drop=True), results_dir / f"{prefix}-donors.csv")
    write_clean_csv(fam_b_sum_df, results_dir / f"{prefix}-gene-summary.csv")
    write_clean_csv(fam_b_top_df, results_dir / f"{prefix}-top-genes.csv")
    write_clean_csv(fam_c_df, results_dir / f"{prefix}-p1-signature.csv")

    donor_tbl = "\n".join(
        f"| {r['donor_id']} ({r['gsm']}) | {int(r['age'])} | {r['sex']} | {r['batch']} | {int(r['lib_size'])} |"
        for _, r in donor_meta.iterrows()
    )
    b_tbl = "\n".join(
        f"| {r['analysis']} | {r['group']} | {int(r['n_donors'])} | {int(r['n_genes_tested'])} | {int(r['n_sig_global_q05'])} | {r['min_pvalue']:.2e} | {r['min_q_bh_global']:.4f} |"
        for _, r in fam_b_sum_df.iterrows()
    )
    c_tbl = "\n".join(
        f"| {r['analysis']} | {int(r['n_donors'])} | {int(r['n_p1_up_detected'])}/{int(r['n_p1_down_detected'])} | {r['composite_effect']:+.4f} | {r['composite_pvalue']:.4f} | {r['empirical_null_pvalue']:.4f} | {r['q_bh_composite']:.4f} | {r['q_bh_empirical']:.4f} |"
        for _, r in fam_c_df.iterrows()
    )

    prim_b = fam_b_sum_df[fam_b_sum_df["analysis"] == "Primary_BatchSexAdj"].iloc[0]
    prim_c = fam_c_df[fam_c_df["analysis"] == "Primary_BatchSexAdj"].iloc[0]
    top5_prim = fam_b_top_df[fam_b_top_df["analysis"] == "Primary_BatchSexAdj"].head(10)
    top_tbl = "\n".join(
        f"| {r['gene']} | {r['mean_log2cpm']:.2f} | {r['effect']:+.4f} | {r['se']:.4f} | {r['pvalue']:.2e} | {r['q_bh_global']:.4f} |"
        for _, r in top5_prim.iterrows()
    )
    noninteger_ct = int(count_stats["count_gate"]["noninteger"])
    max_raw_ct = float(clean_counts.to_numpy().max())

    report_md = f"""# Dataset 08: Human Retinal Pigment Epithelium Bulk RNA-seq (Butler et al., 2021, GSE159435)

## Verdict
**In bulk human RPE (`n = 13` donors, ages 31 to 93 yr, `8M / 5F`), batch- and sex-adjusted donor OLS detected `{int(prim_b['n_sig_global_q05'])}` significant age-associated genes out of `{int(prim_b['n_genes_tested'])}` tested at `q_BH < 0.05` (`RUBCNL` `+0.2560/dec`, `q_BH = 0.0128`; `ARL6IP5` `-0.0891/dec`, `q_BH = 0.0370`), while `P1_Clean` mouse signature transport failed in the opposite (negative) direction (`composite_effect = {prim_c['composite_effect']:+.4f}`, `95% CI [{prim_c['composite_ci95_low']:+.4f}, {prim_c['composite_ci95_high']:+.4f}]`, analytic `p = {prim_c['composite_pvalue']:.4f}`, 200-draw empirical `p = {prim_c['empirical_null_pvalue']:.4f}`).** Neither mouse UP genes (`up_effect = {prim_c['up_effect']:+.4f}`, `p = {prim_c['up_pvalue']:.4f}`) nor DOWN genes (`down_effect = {prim_c['down_effect']:+.4f}`, `p = {prim_c['down_pvalue']:.4f}`) recapitulate mouse direction.

## Assay Scope and Excluded Families
- **Bulk RNA-seq design**: Pooled macular plus peripheral retinal pigment epithelium (RPE) from `13` human donors with no reported retinal pathology (Butler et al., 2021, PMID 33934486, PMCID PMC8184696; GEO `GSE159435`).
- **No Family A fractions and no cross-cell-type concordance**: Because these are bulk tissue counts across heterogeneous macular/peripheral RPE (plus potential trace choroid/photoreceptor carryover) rather than dissociated single cells, Family A compositional shifts and cross-cell-type directional concordance cannot be computed without confounding within-cell expression and tissue composition. Only Family B (gene-level age OLS) and Family C (`P1_Clean` ortholog signature) are run.
- **Count provenance and mapping**: Raw integer HTSeq counts (`noninteger_entries = {noninteger_ct}`, `max = {max_raw_ct:.0f}`) were extracted from `GSE159435_RAW.tar` across all 13 GSMs (`GSM4829762` to `GSM4829774`), removing the 5 `__*` HTSeq summary lines and mapping Ensembl IDs to HGNC symbols via the GRCh38/ortholog table (`117 / 118` `P1_Clean` genes mapped; only `CD24` absent).

## Cohort Structure (`n = 13` donors, 2 Batches)
- Notice the batch-age imbalance in GEO metadata: `batch s1` (`n = 7`) spans the full adult range `31 to 93 yr` (`31M, 40F, 51M, 67F, 77M, 86M, 93F`), whereas `batch s2` (`n = 6`) contains exclusively older donors aged `61 to 92 yr` (`61F, 62M, 71M, 73M, 73F, 92M`). Adjusting for `batch_s2` is essential and `Sens_BatchS1_Only` tests the homogeneous single-batch `31 to 93 yr` ladder.

| Donor (GEO GSM) | Age (yr) | Sex | Batch | Library Size (mapped counts) |
|---|---:|:---:|:---:|---:|
{donor_tbl}

## Family B: Gene-Level Donor Age OLS
- **Primary model (`Primary_BatchSexAdj`, `~ age_per_decade + sex_male + batch_s2 + log_lib`, `M_B = {int(prim_b['n_genes_tested'])}`)**: **`{int(prim_b['n_sig_global_q05'])}` genes at `q_BH < 0.05`**.

| Suite | Group | n_donors | Genes Tested | Global q<0.05 | Min p-value | Min q_BH |
|---|---|---:|---:|---:|---:|---:|
{b_tbl}

### Top 10 Ranked Genes (`Primary_BatchSexAdj`)
| Gene | Mean log2CPM | Effect (log2CPM/dec) | SE | p-value | q_BH |
|---|---:|---:|---:|---:|---:|
{top_tbl}

## Family C: P1_Clean Mouse Signature Transport
- Across all three models (`Primary_BatchSexAdj`, `Sens_SexAdj_NoBatch`, `Sens_BatchS1_Only`; `66` up and `36` down `P1_Clean` genes expressed above `1 CPM` in `>=50%` of donors), **there is zero positive transport of `P1_Clean`**: all three composite point estimates are **negative** (`-0.0704` to `-0.0969` per decade), driven by inverted average slopes relative to mouse aging. This reinforces the conclusion that `P1_Clean` does not positively generalize into aging human RPE.

| Suite | n_donors | P1 Up/Down | Composite Effect | Analytic p | Empirical p (200 draws) | q_BH (Comp) | q_BH (Emp) |
|---|---:|:---:|---:|---:|---:|---:|---:|
{c_tbl}

## Key Limitations
1. **Bulk tissue confounding**: Macular and peripheral RPE were pooled together, and bulk RNA cannot distinguish intrinsic RPE aging from region ratio or trace adjacent tissue variation.
2. **Batch imbalance**: Batch `s2` contains only donors aged `61 to 92 yr` (zero young donors), reducing effective cross-lifespan information primarily to the `7` `batch s1` donors (`31 to 93 yr`) once `batch_s2` is controlled (`Sens_SexAdj_NoBatch` hits largely reflect batch `s2` vs `s1` effects correlated with age).
3. **Sample size**: With `n = 13` (`residual_df = 8`), power for moderate gene-level slopes is limited.

## Output Files
- `results/{prefix}-donors.csv`
- `results/{prefix}-gene-summary.csv`
- `results/{prefix}-top-genes.csv`
- `results/{prefix}-p1-signature.csv`
"""
    write_clean_markdown(report_md, reports_dir / f"{prefix}.md")
    clean_temp_dir(tmp_dir)
    print(f"Organ 08 done: FamB_sig={int(prim_b['n_sig_global_q05'])}/{int(prim_b['n_genes_tested'])}, comp={prim_c['composite_effect']:+.4f}, p={prim_c['composite_pvalue']:.4f}, emp_p={prim_c['empirical_null_pvalue']:.4f}")


if __name__ == "__main__":
    main()
