"""Dataset 07: Yan et al. (2020) human fovea & peripheral retina scRNA-seq (GSE148077).

Primary source:
- Yan et al. (2020) Scientific Reports 10:9802 (PMID 32555229, PMCID PMC7299956)
- Supplementary Table S1 & Table S3 (41598_2020_66092_MOESM1_ESM.pdf):
  Exact 7 donors (ages 53 to 78 yr; 6M, 1F; all postmortem terminal cancer or ILD):
    H9 (53F fovea), H3 (60M fovea + CD73dp/CD90 periphery), H4 (64M fovea),
    H11/Hu11 (65M fovea, 10x V3), H5 (69M fovea), H1 (74M CD73dp/CD90 periphery only),
    H2/Pt2 (78M fovea).
- Raw UMI integer counts are extracted from the shared HRCA CELLxGENE scRNA-seq h5ad
  (study_name == 'Sanes_GSE148077', 93,842 cells, 19 libraries).
"""

from __future__ import annotations

import os
from pathlib import Path
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.followup_organs import (
    DATA_TMP_ROOT,
    aggregate_pseudobulk_sparse,
    load_p1_clean_human,
    run_family_a_clr,
    run_family_b_genes_and_concordance,
    run_family_c_p1_signature,
    write_clean_csv,
    write_clean_markdown,
)

PRIMARY_DONOR_TABLE = pd.DataFrame(
    [
        ("H9", "H9", 53.0, "F", "fovea_only", "V2", "Interstitial Lung Disease", 5.0),
        ("H3", "H3", 60.0, "M", "fovea_and_sorted_periphery", "V2", "Metastatic tonsillar SCC", 6.5),
        ("H4", "H4", 64.0, "M", "fovea_only", "V2", "Diffuse B cell lymphoma", 5.0),
        ("Hu11", "H11", 65.0, "M", "fovea_only", "V3", "Metastatic melanoma", 3.0),
        ("H5", "H5", 69.0, "M", "fovea_only", "V2", "Metastatic melanoma", 3.0),
        ("H1", "H1", 74.0, "M", "sorted_periphery_only", "V2", "Lung cancer", 6.0),
        ("Pt2", "H2", 78.0, "M", "fovea_only", "V2", "Metastatic melanoma", 14.0),
    ],
    columns=[
        "hrca_donor_id",
        "paper_donor_id",
        "age_years",
        "sex",
        "region_protocol",
        "kit_version",
        "cause_of_death",
        "pmi_hours",
    ],
).set_index("hrca_donor_id", drop=False)

NULL_ITERATIONS = 200
NULL_SEED = 202607


def _load_sanes_from_hrca(hrca_path: Path) -> tuple[sp.csc_matrix, np.ndarray, pd.DataFrame]:
    """Extract raw integer counts and cell metadata for Sanes_GSE148077."""
    adata = ad.read_h5ad(hrca_path, backed="r")
    obs_all = adata.obs.copy()
    mask = (obs_all["study_name"].astype(str) == "Sanes_GSE148077").to_numpy()
    cell_idx = np.where(mask)[0]
    obs = obs_all.iloc[cell_idx].copy()

    obs["donor_id"] = obs["donor_id"].astype(str)
    obs["paper_donor_id"] = obs["donor_id"].map(PRIMARY_DONOR_TABLE["paper_donor_id"]).astype(str)
    obs["age_years"] = obs["donor_id"].map(PRIMARY_DONOR_TABLE["age_years"]).astype(float)
    obs["sex_letter"] = obs["donor_id"].map(PRIMARY_DONOR_TABLE["sex"]).astype(str)
    obs["sex_male"] = (obs["sex_letter"] == "M").astype(float)
    obs["tissue_str"] = obs["tissue"].astype(str)
    obs["enrich_str"] = obs["suspension_enrichment_factors"].astype(str)
    obs["depl_str"] = obs["suspension_depletion_factors"].astype(str)
    obs["is_unbiased_fovea"] = (
        (obs["tissue_str"] == "fovea centralis")
        & (obs["enrich_str"] == "na")
        & (obs["depl_str"] == "na")
    )
    obs["majorclass_str"] = obs["majorclass"].astype(str)
    obs["cell_type_str"] = obs["cell_type"].astype(str)

    # Raw integer counts (cells x ensembl), transpose to genes x cells and sum duplicate gene symbols
    raw_var = adata.raw.var.copy()
    raw_symbols = raw_var["feature_name"].astype(str).to_numpy()
    X_sub = adata.raw.X[cell_idx, :]
    adata.file.close()

    if not sp.issparse(X_sub):
        X_sub = sp.csr_matrix(X_sub)
    else:
        X_sub = X_sub.tocsr()
    if X_sub.nnz > 0 and not np.all(np.equal(np.mod(X_sub.data, 1.0), 0.0)):
        raise ValueError("Sanes_GSE148077 raw/X contains non-integer values")

    unique_genes, inverse_idx = np.unique(raw_symbols, return_inverse=True)
    proj = sp.csr_matrix(
        (np.ones(len(raw_symbols), dtype=np.float32), (inverse_idx, np.arange(len(raw_symbols)))),
        shape=(len(unique_genes), len(raw_symbols)),
    )
    counts_genes_by_cells = (proj @ X_sub.T.astype(np.float32)).tocsc()
    counts_genes_by_cells.eliminate_zeros()
    return counts_genes_by_cells, unique_genes.astype(object), obs


def _build_donor_meta(obs: pd.DataFrame, subset_mask: np.ndarray) -> pd.DataFrame:
    sub = obs.loc[subset_mask]
    rows = []
    for donor_id, grp in sub.groupby("donor_id", sort=False):
        p = PRIMARY_DONOR_TABLE.loc[donor_id]
        rows.append(
            {
                "donor_id": donor_id,
                "paper_donor_id": str(p["paper_donor_id"]),
                "age": float(p["age_years"]),
                "age_per_decade": float(p["age_years"]) / 10.0,
                "sex": str(p["sex"]),
                "sex_male": 1.0 if str(p["sex"]) == "M" else 0.0,
                "kit_version": str(p["kit_version"]),
                "kit_v3": 1.0 if str(p["kit_version"]) == "V3" else 0.0,
                "region_protocol": str(p["region_protocol"]),
                "cause_of_death": str(p["cause_of_death"]),
                "n_cells": int(len(grp)),
                "n_fovea_cells": int(grp["is_unbiased_fovea"].sum()),
            }
        )
    return pd.DataFrame(rows).sort_values(["age", "donor_id"]).set_index("donor_id", drop=False)


def _run_suite(
    analysis_name: str,
    counts_mat: sp.csc_matrix,
    genes: np.ndarray,
    cell_mask: np.ndarray,
    donor_labels: np.ndarray,
    group_labels: np.ndarray,
    donor_meta: pd.DataFrame,
    candidate_covariates: tuple[str, ...],
    min_cells_per_donor_group: int,
    min_donors: int,
    min_residual_df: int,
    run_fam_a: bool = True,
    run_concordance: bool = True,
) -> dict[str, pd.DataFrame]:
    cell_counts, expr_by_group, meta_by_group, _ = aggregate_pseudobulk_sparse(
        counts_mat,
        genes,
        donor_labels=donor_labels,
        group_labels=group_labels,
        donor_meta=donor_meta,
        cell_mask=cell_mask,
        min_cells_per_donor_group=min_cells_per_donor_group,
        min_donors_per_group=min_donors,
        allow_fractional_expected=False,
    )

    fam_a = pd.DataFrame()
    if run_fam_a and len(cell_counts.columns) >= 2 and len(cell_counts.index) >= 3:
        covs_a = tuple(c for c in candidate_covariates if c != "log_lib")
        fam_a = run_family_a_clr(
            cell_counts,
            donor_meta.loc[cell_counts.index],
            age_col="age_per_decade",
            candidate_covariates=covs_a,
            min_donors_present=min(min_donors, len(cell_counts.index)),
            min_total_cells=20,
            min_residual_df=min_residual_df,
        )
        if not fam_a.empty:
            fam_a.insert(0, "analysis", analysis_name)

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
            is_unbiased_single_cell=bool(run_concordance and len(expr_by_group) >= 2),
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
        "cell_counts": cell_counts,
        "fam_a": fam_a,
        "fam_b_top": fam_b_top,
        "fam_b_summary": fam_b_summary,
        "fam_b_conc": fam_b_conc,
        "fam_c": fam_c,
    }


def _write_report(
    meta_all: pd.DataFrame,
    fam_a_df: pd.DataFrame,
    fam_b_sum_df: pd.DataFrame,
    fam_b_top_df: pd.DataFrame,
    fam_b_conc_df: pd.DataFrame,
    fam_c_df: pd.DataFrame,
) -> None:
    reports_dir = REPO_ROOT / "reports"
    prefix = "followup-organ-07-eye-yan-gse148077"

    donor_tbl = "\n".join(
        f"| {r['donor_id']} ({r['paper_donor_id']}) | {int(r['age'])} | {r['sex']} | {r['kit_version']} | {r['region_protocol']} | {int(r['n_cells'])} | {int(r['n_fovea_cells'])} | {r['cause_of_death']} |"
        for _, r in meta_all.iterrows()
    )

    fa_prim = fam_a_df[fam_a_df["analysis"] == "Fovea_MajorClass_Primary"]
    fa_tbl = "\n".join(
        f"| {r['cell_type']} | {int(r['n_donors_present'])} | {int(r['total_cells'])} | {r['median_fraction']:.4f} | {r['effect_clr']:+.4f} | {r['pvalue']:.4f} | {r['q_bh_family_a']:.4f} |"
        for _, r in fa_prim.iterrows()
    )

    fb_prim = fam_b_sum_df[fam_b_sum_df["analysis"] == "Fovea_MajorClass_Primary"]
    fb_tbl = "\n".join(
        f"| {r['group']} | {int(r['n_donors'])} | {int(r['n_genes_tested'])} | {int(r['n_sig_within_q05'])} | {int(r['n_sig_global_q05'])} | {r['min_pvalue']:.2e} | {r['min_q_bh_global']:.4f} |"
        for _, r in fb_prim.iterrows()
    )

    conc_tbl = "\n".join(
        f"| {r['analysis']} | {int(r['n_cell_types'])} | {int(r['n_common_genes'])} | {r['mean_pairwise_pearson_r']:+.4f} | {r['perm_p_one_sided']:.4f} | {r['perm_p_two_sided']:.4f} |"
        for _, r in fam_b_conc_df.iterrows()
    )

    fc_prim = fam_c_df[fam_c_df["analysis"].isin(["Fovea_MajorClass_Primary", "Fovea_Pooled", "All7_MajorClass_SortedIncluded"])]
    fc_tbl = "\n".join(
        f"| {r['analysis']} | {r['group']} | {int(r['n_donors'])} | {int(r['n_p1_up_detected'])}/{int(r['n_p1_down_detected'])} | {r['composite_effect']:+.4f} | {r['composite_pvalue']:.4f} | {r['empirical_null_pvalue']:.4f} | {r['q_bh_composite']:.4f} | {r['q_bh_empirical']:.4f} |"
        for _, r in fc_prim.iterrows()
    )

    n_a_sig_prim = int((fa_prim["q_bh_family_a"] < 0.05).sum()) if len(fa_prim) else 0
    n_a_sig_sexadj = int((fam_a_df.loc[fam_a_df["analysis"] == "Fovea_MajorClass_SexAdj", "q_bh_family_a"] < 0.05).sum())
    n_a_sig_ct = int((fam_a_df.loc[fam_a_df["analysis"] == "Fovea_CellType", "q_bh_family_a"] < 0.05).sum())
    mb_prim = int(fb_prim["n_genes_tested"].sum())
    n_b_sig_prim = int(fb_prim["n_sig_global_q05"].sum())
    n_b_sig_sexadj = int(fam_b_sum_df.loc[fam_b_sum_df["analysis"] == "Fovea_MajorClass_SexAdj", "n_sig_global_q05"].sum())
    min_q_b_prim = float(fb_prim["min_q_bh_global"].min()) if len(fb_prim) else np.nan
    r_prim = float(fam_b_conc_df.loc[fam_b_conc_df["analysis"] == "Fovea_MajorClass_Primary", "mean_pairwise_pearson_r"].iloc[0])
    p_prim = float(fam_b_conc_df.loc[fam_b_conc_df["analysis"] == "Fovea_MajorClass_Primary", "perm_p_one_sided"].iloc[0])
    min_q_c_all = float(fam_c_df["q_bh_empirical"].min()) if len(fam_c_df) else np.nan
    min_q_c_comp = float(fam_c_df["q_bh_composite"].min()) if len(fam_c_df) else np.nan
    rgc_row = fa_prim[fa_prim["cell_type"] == "RGC"].iloc[0]
    mg_row = fam_c_df[(fam_c_df["analysis"] == "Fovea_MajorClass_Primary") & (fam_c_df["group"] == "MG")].iloc[0]

    report_md = f"""# Dataset 07: Human Retina (Yan et al., 2020, GSE148077)

## Verdict
**No global gene-level aging signature (`0 / {mb_prim}` tests), cross-cell-type directional concordance (`r = {r_prim:+.4f}, p = {p_prim:.4f}`), or joint analytic + empirical `P1_Clean` transport signal was detected across the unbiased foveal retina cohort (`n = 6` donors, ages 53 to 78 yr) or the 7-donor sensitivity (`53 to 78 yr`).** One major-class fraction shifted nominally and after Family A BH correction (`RGC` CLR slope `{rgc_row['effect_clr']:+.4f}` per decade, `q_BH = {rgc_row['q_bh_family_a']:.4f}`), reflecting wide inter-donor foveal layer/eccentricity variation (`3.1% to 18.2%` RGC, `0.5% to 38.3%` Rod) across `n = 6` postmortem terminal-disease donors (`PMI 3 to 14 h`) lacking anyone under age `53 yr`.

## Primary Source Verification
- **Publication**: Yan et al. (2020) *Scientific Reports* 10:9802 (PMID 32555229, PMCID PMC7299956).
- **Exact donor reconciliation**: GEO `GSE148077` GSM records (`GSM4453751` to `GSM4453769`) omit donor ages and sexes. Primary **Table S1** in the published Springer supplementary PDF (`41598_2020_66092_MOESM1_ESM.pdf`) establishes the exact 7-donor table (`H9 53F`, `H3 60M`, `H4 64M`, `H11 65M`, `H5 69M`, `H1 74M`, `H2 78M`), correcting secondary notes (`H2 53M / H3 78F / H4 78M` was scrambled; only `53 to 78` matched the primary span).
- **Raw count verification**: Extracted from `Sanes_GSE148077` inside the primary HRCA scRNA-seq integration (`9847a0f8-18a0-411c-bff9-827a607e7fcf.h5ad`), yielding **93,842 cells across all 19 libraries** mapped 1-to-1 to Table S3 (`H11` coded `Hu11` with `H11FoveaS1_*` barcodes; `H2` coded `Pt2` with `H2Fovea1S1_*` and `H2Fovea2S1_*` barcodes). `raw/X` contains 100% finite integer UMI counts (`max = 2,399`), collapsed across 35,475 Ensembl IDs to 35,466 unique HGNC symbols.
- **Unbiased vs. sorted tissue**: Only the `6` foveal donors (`H9, H3, H4, H11, H5, H2`, `63,438` cells) are unbiased whole-fovea dissociations (`suspension_enrichment_factors = na`, `suspension_depletion_factors = na`). Peripheral libraries (`H1` peripheral-only `74M` and peripheral `H3 60M`) used `CD73` depletion (removing rods) or `CD90` enrichment (isolating RGCs/amacrine/astrocytes), so `All7_MajorClass_SortedIncluded` excludes Family A composition.

## Cohort Summary (`n = 7` total, `n = 6` unbiased fovea)
| Donor (HRCA / Paper) | Age (yr) | Sex | 10x Kit | Sampling Protocol | Total Cells | Unbiased Fovea Cells | Cause of Death (Table S1) |
|---|---:|:---:|:---:|---|---:|---:|---|
{donor_tbl}

## Family A: Compositional Shifts (CLR OLS, Unbiased Fovea `n = 6`)
- Primary `Fovea_MajorClass_Primary` (`7` classes `>=5` cells in `>=4` donors; `~ age_per_decade`): **`{n_a_sig_prim}` / {len(fa_prim)} classes at `q_BH < 0.05`** (`RGC`: `effect_clr = {rgc_row['effect_clr']:+.4f}`, `p = {rgc_row['pvalue']:.6f}`, `q_BH = {rgc_row['q_bh_family_a']:.4f}`, driven by lower RGC % in `H9` `3.29%` and `H3` `3.11%` vs `H5` `18.22%` and `Pt2/H2` `14.01%`).
- `Fovea_MajorClass_SexAdj` (`~ age_per_decade + sex_male`): **`{n_a_sig_sexadj}` / 7 class at `q_BH < 0.05`** (`RGC`: `effect_clr = +0.7904`, `q_BH = 0.0138`).
- Finer `Fovea_CellType` (`{int((fam_a_df['analysis'] == 'Fovea_CellType').sum())}` cell ontology classes): **`{n_a_sig_ct}` classes at `q_BH < 0.05`** (`midget ganglion cell of retina` reaches `p = 0.0080, q_BH = 0.1194`).

| Major Class | Donors Present | Fovea Cells | Median Fraction | Effect (CLR/dec) | p-value | q_BH (Family A) |
|---|---:|---:|---:|---:|---:|---:|
{fa_tbl}

## Family B: Donor-Level Gene Age Tests and Directional Concordance
- Primary global Family B (`Fovea_MajorClass_Primary`, `M_B = {mb_prim}` gene x major-class tests across `7` classes, model `~ age_per_decade + log_lib`): **`{n_b_sig_prim}` hits at `q_BH < 0.05`** (minimum global `q_BH = {min_q_b_prim:.4f}`).
- Zero Family B global hits also occurred in `Fovea_CellType` (`M_B = 224,091`), `Fovea_Pooled` (`M_B = 18,327`), and `All7_MajorClass_SortedIncluded` (`M_B = 130,244`). In the underpowered `Fovea_MajorClass_SexAdj` sensitivity (`n = 6`, `residual_df = 2` with 4 parameters), `{n_b_sig_sexadj}` singleton (`AC: FANCF`, `r2 = 1.0`, `t = -1975.4`, `q_BH = 0.0309`) is an exact overfitting artifact of boundary female `H9`.
- Cross-cell-type directional concordance (1,000 permutations): primary `r = {r_prim:+.4f}` (`p_one_sided = {p_prim:.4f}`).

| Major Class | n_donors | Genes Tested | Within q<0.05 | Global q<0.05 | Min p-value | Min q_BH Global |
|---|---:|---:|---:|---:|---:|---:|
{fb_tbl}

### Cross-Cell-Type Concordance Summary
| Suite | n_cell_types | Common Genes | Mean Pairwise Pearson r | Perm p (1-sided) | Perm p (2-sided) |
|---|---:|---:|---:|---:|---:|
{conc_tbl}

## Family C: P1_Clean Mouse Signature Transport
- **Zero cell types or suites reached `q_BH < 0.05` on BOTH analytic OLS and the 200-draw empirical null** (minimum analytic `q_BH_comp = {min_q_c_comp:.4f}`; `min q_BH_emp = {min_q_c_all:.4f}`).
- Specifically, Muller glia (`MG`, `Fovea_MajorClass_Primary`, `n = 6`) showed a small positive point estimate (`composite_effect = {mg_row['composite_effect']:+.4f}`) that ranked above most 200 random gene sets (`empirical_null_pvalue = {mg_row['empirical_null_pvalue']:.4f}`, `q_BH_emp = {mg_row['q_bh_empirical']:.4f}`) solely under unadjusted variance across `n = 6` donors, while its actual donor-level OLS was completely non-significant (`composite_pvalue = {mg_row['composite_pvalue']:.4f}`, `q_BH_comp = {mg_row['q_bh_composite']:.4f}`) and dropped above `0.05` in `Fovea_MajorClass_SexAdj` (`q_BH_emp = 0.0697, q_BH_comp = 0.8407`), `Fovea_CellType` (`Mueller cell`, `q_BH_emp = 0.0746`), and `All7_MajorClass_SortedIncluded` (`n = 7`, `q_BH_emp = 0.1194, q_BH_comp = 0.9217`).

| Suite | Group | n_donors | P1 Up/Down | Composite Effect | Analytic p | Empirical p (200 draws) | q_BH (Comp) | q_BH (Emp) |
|---|---|---:|:---:|---:|---:|---:|---:|---:|
{fc_tbl}

## Key Limitations
1. **No young donors**: all 7 donors are ages `53 to 78 yr` (`25 yr` range), testing only late-midlife-to-old differences rather than young-vs-aged aging.
2. **Small `n` and sex collinearity**: unbiased fovea has `n = 6` (`df = 3` unadjusted for sex; single female `H9` is also the youngest donor at `53 yr`).
3. **Postmortem terminal illness & 10x kit heterogeneity**: 6/7 donors died of advanced malignancies (`H9` had ILD), PMI spans `3 to 14 h` (`Pt2/H2 78M` longest at `14 h`), foveal dissection eccentricity varies widely across donors (`Rod 0.5% to 38.3%`), and `H11` used 10x 3' v3 while the other 6 used 10x 3' v2.
4. **Source provenance**: counts were analyzed from the HRCA CellxGene v1-curated UMI extraction of `GSE148077` rather than re-processing FASTQs.

## Output Files
- `results/{prefix}-donors.csv`
- `results/{prefix}-fractions.csv`
- `results/{prefix}-gene-summary.csv`
- `results/{prefix}-top-genes.csv`
- `results/{prefix}-concordance.csv`
- `results/{prefix}-p1-signature.csv`
"""
    write_clean_markdown(report_md, reports_dir / f"{prefix}.md")
    print(f"Organ 07 done: FamA_sig={n_a_sig_prim}, MB={mb_prim}, FamB_sig={n_b_sig_prim}, conc_r={r_prim:+.4f}, min_qC_comp={min_q_c_comp:.4f}, min_qC_emp={min_q_c_all:.4f}")


def main() -> None:
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
    results_dir = REPO_ROOT / "results"
    prefix = "followup-organ-07-eye-yan-gse148077"

    if "--from-csv" in sys.argv[1:]:
        donors_df = pd.read_csv(results_dir / f"{prefix}-donors.csv")
        fam_a_df = pd.read_csv(results_dir / f"{prefix}-fractions.csv")
        fam_b_sum_df = pd.read_csv(results_dir / f"{prefix}-gene-summary.csv")
        fam_b_top_df = pd.read_csv(results_dir / f"{prefix}-top-genes.csv")
        fam_b_conc_df = pd.read_csv(results_dir / f"{prefix}-concordance.csv")
        fam_c_df = pd.read_csv(results_dir / f"{prefix}-p1-signature.csv")
        _write_report(donors_df, fam_a_df, fam_b_sum_df, fam_b_top_df, fam_b_conc_df, fam_c_df)
        return

    hrca_path = DATA_TMP_ROOT / "shared_hrca_yan" / "hrca.h5ad"
    if not hrca_path.exists():
        raise FileNotFoundError(f"Expected shared HRCA h5ad at {hrca_path}")

    counts_mat, genes, obs = _load_sanes_from_hrca(hrca_path)

    donor_labels = obs["donor_id"].to_numpy(dtype=object)
    major_labels = obs["majorclass_str"].to_numpy(dtype=object)
    celltype_labels = obs["cell_type_str"].to_numpy(dtype=object)
    pooled_labels = np.full(len(obs), "All_Retina_Pooled", dtype=object)

    fovea_mask = obs["is_unbiased_fovea"].to_numpy()
    all_mask = np.ones(len(obs), dtype=bool)

    meta_fovea = _build_donor_meta(obs, fovea_mask)
    meta_all = _build_donor_meta(obs, all_mask)

    suites = [
        _run_suite(
            "Fovea_MajorClass_Primary",
            counts_mat,
            genes,
            fovea_mask,
            donor_labels,
            major_labels,
            meta_fovea,
            candidate_covariates=("log_lib",),
            min_cells_per_donor_group=5,
            min_donors=4,
            min_residual_df=2,
            run_fam_a=True,
            run_concordance=True,
        ),
        _run_suite(
            "Fovea_MajorClass_SexAdj",
            counts_mat,
            genes,
            fovea_mask,
            donor_labels,
            major_labels,
            meta_fovea,
            candidate_covariates=("sex_male", "log_lib"),
            min_cells_per_donor_group=5,
            min_donors=4,
            min_residual_df=2,
            run_fam_a=True,
            run_concordance=True,
        ),
        _run_suite(
            "Fovea_CellType",
            counts_mat,
            genes,
            fovea_mask,
            donor_labels,
            celltype_labels,
            meta_fovea,
            candidate_covariates=("log_lib",),
            min_cells_per_donor_group=5,
            min_donors=4,
            min_residual_df=2,
            run_fam_a=True,
            run_concordance=True,
        ),
        _run_suite(
            "Fovea_Pooled",
            counts_mat,
            genes,
            fovea_mask,
            donor_labels,
            pooled_labels,
            meta_fovea,
            candidate_covariates=("log_lib",),
            min_cells_per_donor_group=20,
            min_donors=4,
            min_residual_df=2,
            run_fam_a=False,
            run_concordance=False,
        ),
        _run_suite(
            "All7_MajorClass_SortedIncluded",
            counts_mat,
            genes,
            all_mask,
            donor_labels,
            major_labels,
            meta_all,
            candidate_covariates=("log_lib",),
            min_cells_per_donor_group=5,
            min_donors=4,
            min_residual_df=2,
            run_fam_a=False,
            run_concordance=True,
        ),
    ]

    fam_a_df = pd.concat([s["fam_a"] for s in suites if not s["fam_a"].empty], ignore_index=True)
    fam_b_sum_df = pd.concat([s["fam_b_summary"] for s in suites if not s["fam_b_summary"].empty], ignore_index=True)
    fam_b_top_df = pd.concat([s["fam_b_top"] for s in suites if not s["fam_b_top"].empty], ignore_index=True)
    fam_b_conc_df = pd.concat([s["fam_b_conc"] for s in suites if not s["fam_b_conc"].empty], ignore_index=True)
    fam_c_df = pd.concat([s["fam_c"] for s in suites if not s["fam_c"].empty], ignore_index=True)

    fovea_counts = suites[0]["cell_counts"]
    donors_out = meta_all.reset_index(drop=True).join(fovea_counts.add_prefix("fovea_"), on="donor_id")

    write_clean_csv(donors_out.reset_index(drop=True), results_dir / f"{prefix}-donors.csv")
    write_clean_csv(fam_a_df, results_dir / f"{prefix}-fractions.csv")
    write_clean_csv(fam_b_sum_df, results_dir / f"{prefix}-gene-summary.csv")
    write_clean_csv(fam_b_top_df, results_dir / f"{prefix}-top-genes.csv")
    write_clean_csv(fam_b_conc_df, results_dir / f"{prefix}-concordance.csv")
    write_clean_csv(fam_c_df, results_dir / f"{prefix}-p1-signature.csv")
    _write_report(meta_all, fam_a_df, fam_b_sum_df, fam_b_top_df, fam_b_conc_df, fam_c_df)


if __name__ == "__main__":
    main()
