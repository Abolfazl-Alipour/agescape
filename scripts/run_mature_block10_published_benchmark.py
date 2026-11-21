#!/usr/bin/env python3
"""Mature (30-plus) sensitivity rerun of the frozen Block 10 published-signature benchmark.

Only change vs the frozen original: scoring is restricted to donors aged 30-plus
in both cohorts (Gao heart, Yazar blood) in both arms, and matched nulls are
recomputed on the 30-plus donors with 2000 draws. The four published gene lists,
directions, coverage rule, eligibility floors, and all 32 multiplicity families
are unchanged. The four external sets were NOT rederived from mature donors, so
this is a mature-TESTED sensitivity only, not a mature-derived benchmark.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import os
import sys
from pathlib import Path

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_var] = "1"

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem
from scipy import stats
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_followup_hcata_pool_sensitivity import load_cache
from scripts.run_followup_organ_17_heart_gao import (
    H5AD_PATH,
    STUDY_DUMMY_MAP,
    build_donor_audit,
    build_donor_meta_table,
    stream_pseudobulk_counts,
)
from src.followup_organs import (
    finalize_pseudobulk_group,
    run_family_c_p1_signature,
    select_design_columns,
    write_clean_csv,
    write_clean_markdown,
)

RESULTS = REPO / "results"
REGISTRY_CSV = RESULTS / "followup-mature-block10-benchmark-signatures.csv"
RESULTS_CSV = RESULTS / "followup-mature-block10-benchmark-results.csv"
REPORT = REPO / "reports" / "followup-mature-block10.md"
ORIGINAL_RESULTS_CSV = RESULTS / "followup-11pt-block10-benchmark-results.csv"
DONOR_COUNTS_CSV = RESULTS / "followup-mature-block10-donor-counts.csv"
SOURCE_DIR = Path("/home/reza/projects/agescape-revision-data/tmp_organs/block10_signatures")
YAZAR_CACHE = Path("/home/reza/projects/agescape-revision-data/tmp_organs/hcata/yazar_cache")
ITERATIONS = 2000
SEED = 20261105
MIN_DONORS = 6
MIN_CELLS = 20
MIN_AGE_YEARS = 30.0
COHORTS = ("Gao_Heart", "Yazar_Blood")
ARMS = ("unadjusted", "batch_adjusted")
SIGNATURE_ORDER = (
    "Harries2011_6gene",
    "deMagalhaes2009_73gene",
    "Sood2015_150probe",
    "Peters2015_1497gene",
)

SOURCE_URLS = {
    "Harries2011_6gene": "https://pmc.ncbi.nlm.nih.gov/articles/PMC3173580/",
    "deMagalhaes2009_73gene_up": "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=DEMAGALHAES_AGING_UP&fileType=gmt",
    "deMagalhaes2009_73gene_down": "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=DEMAGALHAES_AGING_DN&fileType=gmt",
    "Sood2015_150probe": "https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-015-0750-x/MediaObjects/13059_2015_750_MOESM1_ESM.xlsx",
    "Peters2015_1497gene": "https://static-content.springer.com/esm/art%3A10.1038%2Fncomms9570/MediaObjects/41467_2015_BFncomms9570_MOESM436_ESM.xlsx",
}

CITATIONS = {
    "Harries2011_6gene": (
        "Harries LW, Hernandez D, Henley W, Wood AR, Holly AC, Bradley-Smith RM, "
        "Yaghootkar H, Dutta A, Murray A, Frayling TM, Guralnik JM, Bandinelli S, "
        "Singleton A, Ferrucci L, Melzer D. 2011. Human aging is characterized by focused "
        "changes in gene expression and deregulation of alternative splicing. Aging Cell "
        "10(5):868-878. DOI: 10.1111/j.1474-9726.2011.00726.x."
    ),
    "deMagalhaes2009_73gene": (
        "de Magalhaes JP, Curado J, Church GM. 2009. Meta-analysis of age-related gene "
        "expression profiles identifies common signatures of aging. Bioinformatics "
        "25(7):875-881. DOI: 10.1093/bioinformatics/btp073."
    ),
    "Sood2015_150probe": (
        "Sood S, Gallagher IJ, Lunnon K, Rullman E, Keohane A, Crossland H, Phillips BE, "
        "Cederholm T, Jensen T, van Loon LJC, Lannfelt L, Kraus WE, Atherton PJ, Howard R, "
        "Gustafsson T, Hodges A, Timmons JA. "
        "2015. A novel multi-tissue RNA diagnostic of healthy ageing relates to cognitive "
        "health status. Genome Biology 16:185. DOI: 10.1186/s13059-015-0750-x."
    ),
    "Peters2015_1497gene": (
        "Peters MJ, Joehanes R, Pilling LC, Schurmann C, Conneely KN, Powell J, Reinmaa E, "
        "et al. 2015. The transcriptional landscape of age in human peripheral blood. "
        "Nature Communications 6:8570. DOI: 10.1038/ncomms9570."
    ),
}

SOURCE_DESCRIPTIONS = {
    "Harries2011_6gene": (
        "Main article section 'Six of the most deregulated 25 probes provide good discrimination "
        "between young and old', with directions from Table 1a and Table 1b. Exact model genes: "
        "LRRN3, CD248, CCR6, GRAP, VAMP5, CD27."
    ),
    "deMagalhaes2009_73gene": (
        "MSigDB C2 CGP gene sets DEMAGALHAES_AGING_UP and DEMAGALHAES_AGING_DN, contributed "
        "by the paper author and sourced to the paper's published meta-signature. The 56 up and "
        "17 down source identifiers map to 54 and 16 current unique human symbols in MSigDB."
    ),
    "Sood2015_150probe": (
        "Additional file 1, worksheet 'Top 150', columns 'HU133 Probeset_ID', 'Gene Symbol', "
        "and 'Ratio of Y:0 muscle'. The 150 probe sets comprise 148 unique symbols because "
        "RUNX1 and FAM20C each occur twice."
    ),
    "Peters2015_1497gene": (
        "Supplementary Data 1, worksheet 'Data.1'. Select rows with sign-DISC = y and "
        "sign-REPL = y, then use NEW-Gene-ID and Direction. This yields exactly 1497 unique "
        "genes, 600 positive and 897 negative, matching the paper."
    ),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_gmt(path: Path) -> list[str]:
    fields = path.read_text(encoding="utf-8").strip().split("\t")
    genes = [str(g).strip() for g in fields[2:] if str(g).strip()]
    assert genes and len(genes) == len(set(genes))
    return genes


def load_frozen_signatures() -> dict[str, dict[str, object]]:
    required = [
        SOURCE_DIR / "demagalhaes_aging_up.gmt",
        SOURCE_DIR / "demagalhaes_aging_dn.gmt",
        SOURCE_DIR / "sood_2015_additional_file_1.xlsx",
        SOURCE_DIR / "peters_2015_supplementary_data_1.xlsx",
    ]
    for path in required:
        if not path.exists():
            raise FileNotFoundError(f"Missing frozen published source file: {path}")

    signatures: dict[str, dict[str, object]] = {
        "Harries2011_6gene": {
            "up": ["VAMP5"],
            "down": ["LRRN3", "CD248", "CCR6", "GRAP", "CD27"],
            "n_published": 6,
            "n_published_up": 1,
            "n_published_down": 5,
            "published_unit": "genes",
            "species": "human",
            "source_url": SOURCE_URLS["Harries2011_6gene"],
            "source_paths": [],
        }
    }

    up = read_gmt(SOURCE_DIR / "demagalhaes_aging_up.gmt")
    down = read_gmt(SOURCE_DIR / "demagalhaes_aging_dn.gmt")
    assert (len(up), len(down)) == (54, 16)
    signatures["deMagalhaes2009_73gene"] = {
        "up": up,
        "down": down,
        "n_published": 73,
        "n_published_up": 56,
        "n_published_down": 17,
        "published_unit": "source identifiers",
        "species": "human",
        "source_url": SOURCE_URLS["deMagalhaes2009_73gene_up"] + " ; " + SOURCE_URLS["deMagalhaes2009_73gene_down"],
        "source_paths": [SOURCE_DIR / "demagalhaes_aging_up.gmt", SOURCE_DIR / "demagalhaes_aging_dn.gmt"],
    }

    sood_path = SOURCE_DIR / "sood_2015_additional_file_1.xlsx"
    sood = pd.read_excel(sood_path, sheet_name="Top 150", header=1)
    assert len(sood) == 150
    assert sood["Gene Symbol"].nunique() == 148
    assert sood["Ratio of Y:0 muscle"].value_counts().to_dict() == {"down": 138, "up": 12}
    sood_up = sood.loc[sood["Ratio of Y:0 muscle"] == "up", "Gene Symbol"].astype(str).drop_duplicates().tolist()
    sood_down = sood.loc[sood["Ratio of Y:0 muscle"] == "down", "Gene Symbol"].astype(str).drop_duplicates().tolist()
    assert (len(sood_up), len(sood_down)) == (12, 136)
    signatures["Sood2015_150probe"] = {
        "up": sood_up,
        "down": sood_down,
        "n_published": 150,
        "n_published_up": 12,
        "n_published_down": 138,
        "published_unit": "probe sets",
        "species": "human",
        "source_url": SOURCE_URLS["Sood2015_150probe"],
        "source_paths": [sood_path],
    }

    peters_path = SOURCE_DIR / "peters_2015_supplementary_data_1.xlsx"
    peters = pd.read_excel(peters_path, sheet_name="Data.1", header=2)
    selected = peters[(peters["sign-DISC"] == "y") & (peters["sign-REPL"] == "y")].copy()
    assert len(selected) == 1497
    assert selected["NEW-Gene-ID"].nunique() == 1497
    assert selected["Direction"].value_counts().to_dict() == {"-": 897, "+": 600}
    assert (selected["Direction"] == selected["META:Direction"]).all()
    peters_up = selected.loc[selected["Direction"] == "+", "NEW-Gene-ID"].astype(str).tolist()
    peters_down = selected.loc[selected["Direction"] == "-", "NEW-Gene-ID"].astype(str).tolist()
    signatures["Peters2015_1497gene"] = {
        "up": peters_up,
        "down": peters_down,
        "n_published": 1497,
        "n_published_up": 600,
        "n_published_down": 897,
        "published_unit": "genes",
        "species": "human",
        "source_url": SOURCE_URLS["Peters2015_1497gene"],
        "source_paths": [peters_path],
    }

    assert tuple(signatures) == SIGNATURE_ORDER
    for signature_id, program in signatures.items():
        up_list = list(program["up"])
        down_list = list(program["down"])
        assert len(up_list) == len(set(up_list))
        assert len(down_list) == len(set(down_list))
        assert not set(up_list).intersection(down_list)
        program["n_unique"] = len(up_list) + len(down_list)
        program["source_sha256"] = ";".join(sha256(path) for path in program["source_paths"])
    return signatures


def filter_cell_floor(
    expr: dict[str, pd.DataFrame], metas: dict[str, pd.DataFrame]
) -> tuple[dict[str, pd.DataFrame], dict[str, pd.DataFrame]]:
    filtered_expr = {}
    filtered_meta = {}
    for group, meta in metas.items():
        keep = meta.index[meta["n_cells"].astype(int) >= MIN_CELLS].tolist()
        if len(keep) >= MIN_DONORS:
            filtered_meta[group] = meta.loc[keep].copy()
            filtered_expr[group] = expr[group][keep]
    return filtered_expr, filtered_meta


def load_gao_inputs():
    print("Gao: reading metadata and streaming full-depth counts", flush=True)
    with h5py.File(H5AD_PATH, "r") as handle:
        obs = read_elem(handle["obs"])
        var = read_elem(handle["var"])
    audit = build_donor_audit(obs)
    donor_meta = build_donor_meta_table(audit)
    included = set(donor_meta.index)
    keep = obs["donor_id"].astype(str).isin(included).to_numpy()
    raw, counts, _ = stream_pseudobulk_counts(
        H5AD_PATH,
        keep,
        obs["donor_id"].astype(str).to_numpy(),
        obs["final_cell_type"].astype(str).to_numpy(),
        var.index.astype(str).to_numpy(),
    )
    expr = {}
    metas = {}
    for group, raw_group in raw.items():
        donors = list(raw_group.columns)
        cell_count = counts.loc[donors, "Total_Cells"] if group == "All_LV_Nuclei" else counts.loc[donors, group]
        expr[group], metas[group], _ = finalize_pseudobulk_group(raw_group, donor_meta.loc[donors], cell_count)
    expr, metas = filter_cell_floor(expr, metas)
    covariates = ("sex_male", "log_lib", *STUDY_DUMMY_MAP.values())
    del raw, counts, obs, var
    gc.collect()
    return {
        "unadjusted": (expr, metas, ("sex_male", "log_lib"), 2),
        "batch_adjusted": (expr, metas, covariates, 2),
    }


def load_yazar_inputs():
    print("Yazar: loading cached full-depth pseudobulk", flush=True)
    donor_meta, _, expr, metas = load_cache(
        YAZAR_CACHE, drop_count_cols=("peripheral blood mononuclear cell",)
    )
    expr, metas = filter_cell_floor(expr, metas)
    all_pools = sorted(donor_meta["pool_number"].dropna().astype(int).unique())
    all_pool_cols = [f"pool_{p}" for p in all_pools]
    fe_expr = {}
    fe_meta = {}
    for group, meta in metas.items():
        pool_counts = meta["pool_number"].value_counts()
        ge2 = pool_counts[pool_counts >= 2]
        keep = meta.index[meta["pool_number"].isin(ge2.index)].tolist()
        estimated_df = len(keep) - (4 + max(0, len(ge2) - 1))
        if len(keep) < MIN_DONORS or len(ge2) < 2 or estimated_df < 10:
            continue
        sub = meta.loc[keep].copy()
        reference = int(ge2.index[0])
        for pool in all_pools:
            sub[f"pool_{pool}"] = (
                (sub["pool_number"].astype(int) == pool).astype(float)
                if pool != reference
                else 0.0
            )
        fe_expr[group] = expr[group][keep]
        fe_meta[group] = sub
    return {
        "unadjusted": (expr, metas, ("sex_male", "log_lib"), 2),
        "batch_adjusted": (fe_expr, fe_meta, ("sex_male", "log_lib", *all_pool_cols), 10),
    }


def restrict_to_mature(cohort_inputs):
    """Restrict every arm to donors aged 30-plus; return restricted inputs, audit rows, totals.

    Only donor filtering changes: the same expr/metas objects from the frozen
    loaders are subset to donors with age_years >= 30.0 in both cohorts and
    both arms. Pool and study fixed-effect columns built by the loaders are
    kept; columns that become constant after restriction are dropped later by
    select_design_columns inside eligibility checks. Matched nulls are
    recomputed on the restricted donors by run_family_c_p1_signature.
    """
    rows = []
    totals = {}
    restricted = {}
    for cohort, arms in cohort_inputs.items():
        restricted[cohort] = {}
        for arm, (expr, metas, covariates, min_residual_df) in arms.items():
            new_expr = {}
            new_metas = {}
            for group, meta in metas.items():
                ages = meta["age_years"].astype(float)
                keep = meta.index[ages >= MIN_AGE_YEARS].tolist()
                assert ages.loc[keep].ge(MIN_AGE_YEARS).all()
                rows.append({
                    "cohort": cohort,
                    "arm": arm,
                    "unit": group,
                    "n_donors_full_range": len(meta),
                    "n_donors_30plus": len(keep),
                    "n_donors_excluded_under_30": len(meta) - len(keep),
                    "min_age_years_30plus": float(ages.loc[keep].min()) if keep else float("nan"),
                    "max_age_years_30plus": float(ages.loc[keep].max()) if keep else float("nan"),
                })
                new_metas[group] = meta.loc[keep].copy()
                new_expr[group] = expr[group].loc[:, keep]
            restricted[cohort][arm] = (new_expr, new_metas, covariates, min_residual_df)
        donors_before = set()
        donors_after = set()
        for arm in ARMS:
            for meta in cohort_inputs[cohort][arm][1].values():
                donors_before.update(meta.index)
            for meta in restricted[cohort][arm][1].values():
                donors_after.update(meta.index)
        totals[cohort] = (len(donors_before), len(donors_after))
        print(f"{cohort}: {len(donors_after)} of {len(donors_before)} unique donors at 30-plus", flush=True)
    return restricted, pd.DataFrame(rows), totals


def expressed_pool(log2cpm: pd.DataFrame) -> pd.DataFrame:
    threshold = float(np.log2(1.0 + 1.0))
    min_req = max(3, int(np.ceil(0.5 * log2cpm.shape[1])))
    keep = (log2cpm.to_numpy(dtype=float) >= threshold).sum(axis=1) >= min_req
    return log2cpm.loc[keep]


def eligible_units(
    expr: dict[str, pd.DataFrame],
    metas: dict[str, pd.DataFrame],
    covariates: tuple[str, ...],
    min_residual_df: int,
    program: dict[str, object],
) -> tuple[set[str], dict[str, set[str]]]:
    up_set = set(program["up"])
    down_set = set(program["down"])
    eligible = set()
    scored: dict[str, set[str]] = {}
    for group, log2cpm in expr.items():
        meta = metas[group]
        common = [d for d in meta.index if d in log2cpm.columns]
        if len(common) < MIN_DONORS:
            continue
        mat = expressed_pool(log2cpm[common])
        if mat.shape[0] < 200:
            continue
        genes = set(mat.index.astype(str))
        if not genes.intersection(up_set) or not genes.intersection(down_set):
            continue
        try:
            select_design_columns(
                meta.loc[common],
                age_col="age_per_decade",
                candidate_covariates=covariates,
                min_residual_df=min_residual_df,
            )
        except ValueError:
            continue
        eligible.add(group)
        scored[group] = genes.intersection(up_set | down_set)
    return eligible, scored


def build_preflight(signatures, cohort_inputs) -> tuple[pd.DataFrame, dict]:
    rows = []
    details = {}
    for signature_id in SIGNATURE_ORDER:
        program = signatures[signature_id]
        details[signature_id] = {}
        for cohort in COHORTS:
            arm_eligible = {}
            arm_scored = {}
            for arm in ARMS:
                arm_eligible[arm], arm_scored[arm] = eligible_units(
                    *cohort_inputs[cohort][arm], program
                )
            units = sorted(arm_eligible["unadjusted"] & arm_eligible["batch_adjusted"])
            all_genes = set()
            for arm in ARMS:
                for group in units:
                    all_genes.update(arm_scored[arm].get(group, set()))
            cohort_gene_universe = set()
            for arm in ARMS:
                expr = cohort_inputs[cohort][arm][0]
                for frame in expr.values():
                    cohort_gene_universe.update(frame.index.astype(str))
            published_symbols = set(program["up"]) | set(program["down"])
            detected = published_symbols & cohort_gene_universe
            coverage = len(detected) / int(program["n_unique"])
            status = "TESTABLE" if coverage >= 0.5 else "UNTESTABLE_COVERAGE"
            if status != "TESTABLE":
                units = []
            details[signature_id][cohort] = {
                "units": units,
                "coverage": coverage,
                "detected": detected,
                "scored": all_genes,
                "status": status,
            }
            rows.append({
                "signature_id": signature_id,
                "cohort": cohort,
                "n_detected": len(detected),
                "n_scored": len(all_genes) if status == "TESTABLE" else 0,
                "coverage_fraction": coverage,
                "cohort_status": status,
                "family_size": len(units),
                "units": ";".join(units),
            })
    return pd.DataFrame(rows), details


def registry_frame(signatures, preflight: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for signature_id in SIGNATURE_ORDER:
        program = signatures[signature_id]
        row = {
            "signature_id": signature_id,
            "citation": CITATIONS[signature_id],
            "year": np.nan,
            "journal": "",
            "doi": CITATIONS[signature_id].split("DOI: ")[1].rstrip("."),
            "source_url": program["source_url"],
            "source_file_or_table": SOURCE_DESCRIPTIONS[signature_id],
            "source_sha256": program["source_sha256"],
            "retrieval_status": "RETRIEVED_AND_FROZEN_BEFORE_INFERENCE",
            "species": program["species"],
            "ortholog_mapping": "not_applicable_human_symbols_used_as_published",
            "published_unit": program["published_unit"],
            "n_genes_published": program["n_published"],
            "n_published_up": program["n_published_up"],
            "n_published_down": program["n_published_down"],
            "n_unique_symbols_used": program["n_unique"],
            "n_up_used": len(program["up"]),
            "n_down_used": len(program["down"]),
            "n_genes_mapped": program["n_unique"],
            "direction_definition": "composite = mean log2CPM of published age-up genes minus mean log2CPM of published age-down genes",
            "genes_up": ";".join(program["up"]),
            "genes_down": ";".join(program["down"]),
        }
        journals = {
            "Harries2011_6gene": "Aging Cell",
            "deMagalhaes2009_73gene": "Bioinformatics",
            "Sood2015_150probe": "Genome Biology",
            "Peters2015_1497gene": "Nature Communications",
        }
        years = {
            "Harries2011_6gene": 2011,
            "deMagalhaes2009_73gene": 2009,
            "Sood2015_150probe": 2015,
            "Peters2015_1497gene": 2015,
        }
        row["journal"] = journals[signature_id]
        row["year"] = years[signature_id]
        for cohort in COHORTS:
            p = preflight[(preflight["signature_id"] == signature_id) & (preflight["cohort"] == cohort)].iloc[0]
            prefix = "gao" if cohort == "Gao_Heart" else "yazar"
            row[f"{prefix}_n_detected"] = int(p["n_detected"])
            row[f"{prefix}_n_scored"] = int(p["n_scored"])
            row[f"{prefix}_coverage_fraction"] = float(p["coverage_fraction"])
            row[f"{prefix}_status"] = p["cohort_status"]
            for arm in ARMS:
                row[f"{prefix}_{arm}_composite_family_size"] = int(p["family_size"])
                row[f"{prefix}_{arm}_empirical_family_size"] = int(p["family_size"])
        rows.append(row)
    return pd.DataFrame(rows)


def markdown_gene_list(label: str, genes: list[str]) -> str:
    return f"**{label} ({len(genes)}):** " + ", ".join(genes)


def frozen_report_text(registry: pd.DataFrame, preflight: pd.DataFrame, totals: dict) -> str:
    lines = [
        "# Block 10 mature rerun: Published-signature benchmark at 30-plus",
        "",
        "## Pre-analysis freeze",
        "",
        "This registry, the exact gene lists, the analysis specification, the coverage rule, and all multiplicity families were frozen before any composite slope or matched-null model was fit. Four published human signatures were selected before running. No candidate required substitution because all four lists and directions were retrievable. No P1 or P0 genes were used as a benchmark signature.",
        "",
        "## Mature sensitivity: what changed and what did not",
        "",
        "Scoring is restricted to donors aged 30-plus in both primary cohorts (Gao heart, Yazar blood) and in both arms. Matched nulls are recomputed on the 30-plus donors with 2000 draws and the same fixed seed 20261105. Every other frozen rule is unchanged: the four published gene lists and directions, the 50 percent coverage rule, the eligibility floors (at least 20 cells per donor per unit and at least 6 donors per unit), the Block 5 composite and empirical tests, and BH correction within each signature, cohort, and arm across the frozen units, separately for composite and empirical tests.",
        "",
        "Limitation, stated explicitly: the four external published signatures were NOT rederived from mature donors. Their gene lists and directions come from the original publications, whose derivations used younger donors. This rerun is a mature-TESTED sensitivity only: it asks whether the published signatures still pass when scored on donors aged 30-plus. It is not a mature-derived benchmark. Full compliance with a mature-to-old derivation would require rederiving each signature from its source study restricted to mature donors, which is outside this sensitivity.",
        "",
        f"Unique donors scored: Gao heart {totals['Gao_Heart'][1]} of {totals['Gao_Heart'][0]} full-range donors at 30-plus; Yazar blood {totals['Yazar_Blood'][1]} of {totals['Yazar_Blood'][0]} full-range donors at 30-plus. Donor restriction audit: `results/followup-mature-block10-donor-counts.csv` records, for every cohort, arm, and unit, the donor count at full range and at 30-plus. The unit is the donor, never the cell.",
        "",
        "### Frozen registry and citations",
        "",
    ]
    for signature_id in SIGNATURE_ORDER:
        row = registry[registry["signature_id"] == signature_id].iloc[0]
        lines.extend([
            f"### {signature_id}",
            "",
            CITATIONS[signature_id],
            "",
            f"Source: {SOURCE_DESCRIPTIONS[signature_id]}",
            "",
            f"Retrieval URL: {row['source_url']}",
            "",
            f"Published content: {row['n_genes_published']} {row['published_unit']}, {row['n_published_up']} up and {row['n_published_down']} down. Scoring uses {row['n_unique_symbols_used']} unique human symbols, {row['n_up_used']} up and {row['n_down_used']} down.",
            "",
            markdown_gene_list("Exact up list used", list(row["genes_up"].split(";"))),
            "",
            markdown_gene_list("Exact down list used", list(row["genes_down"].split(";"))),
            "",
        ])
    lines.extend([
        "## Exact analysis specification frozen before running",
        "",
        "The unit is the donor, never the cell. The analysis reuses `run_family_c_p1_signature` from the Block 5 pipeline without changing its scoring or random-set code. Full-depth raw counts are aggregated to donor pseudobulk, the human technical mask is applied before normalization, and expression is log2(CPM + 1). There is no downsampling or rarefaction. The mask excludes uppercase MT-, RPL/RPS, MRPL/MRPS, XIST, TSIX, DDX3Y, EIF2S3Y, UTY, and KDM5D, plus ERCC features through the Block 5 extended mask.",
        "",
        "Only donors aged 30-plus are scored, in both cohorts and both arms. For each signature, cohort, arm, and cell-type unit, the composite score is mean published age-up log2CPM minus mean published age-down log2CPM. The OLS predictor is age_per_decade. Both arms include sex_male and log_lib. The batch-adjusted arm adds Gao study fixed effects for seven studies (six indicator columns) or Yazar pool fixed effects, using the Block 5 pool construction and minimum residual df of 10. The unadjusted arm removes only batch covariates. Gao uses minimum residual df 2. The expression filter is the Block 5 filter: at least 1 CPM in at least half of donors, at least 200 genes in the expression pool, and at least one detected up and down signature gene.",
        "",
        "Shared eligibility floors are at least 20 cells per donor per unit and at least 6 donors per unit. Units below either floor are excluded rather than treated as nulls. For each signature and cohort, only units eligible in both arms enter the batch contrast, so arm comparisons contain the same units. Yazar uses pool fixed effects and excludes the peripheral blood mononuclear cell doublet label as in Block 5.",
        "",
        "The matched random-set comparison uses exactly 2000 iterations, matched separately on up and down set sizes and expression deciles, with fixed seed 20261105 for every fit. Composite empirical p-values use the Block 5 centered two-sided matched-null statistic. A unit passes only when composite BH q < 0.05 and empirical BH q < 0.05. MDE80 is (t(0.975, df) + t(0.80, df)) times the composite standard error and is reported for every tested unit, including every null.",
        "",
        "Coverage denominator is the number of unique published human symbols used. N detected is the number present in the cohort post-mask gene universe. N scored is the number detected in at least one eligible unit expression pool. Per-unit scored counts are also reported. Coverage below 50 percent makes that signature UNTESTABLE in that cohort, not null. All four signatures are human, so mouse one-to-one ortholog mapping is not applicable; n_genes_mapped equals the unique symbols used.",
        "",
        "The optional paper human cohorts were not added. The complete MSSM/HBCC brain, kidney, and lung extension does not have one already cached batch-contrast pipeline with these exact covariates across all those cohorts. No partial secondary arm or new download was substituted. The two required primary cohorts use only existing local data.",
        "",
        "A split-half technical null is not part of this transport design; the required design-specific null is the 2000-iteration expression-matched random-set comparison. Age bins are not used because the prespecified model is continuous age per decade.",
        "",
        "## Multiplicity families frozen before correction",
        "",
        "BH correction is within each signature, cohort, and arm across the frozen units, separately for composite and empirical tests. There are 32 families: 4 signatures x 2 cohorts x 2 arms x 2 test types. An UNTESTABLE coverage slice has family size zero and no q-values. Family sizes are identical for composite and empirical tests within a slice. Family sizes at 30-plus are in the table below and in the registry CSV; a family can shrink where units become ineligible at 30-plus.",
        "",
        "| Signature | Cohort | Coverage | Detected | Scored | Status | Units per arm and test family |",
        "| --- | --- | ---: | ---: | ---: | --- | ---: |",
    ])
    for r in preflight.itertuples(index=False):
        lines.append(
            f"| {r.signature_id} | {r.cohort} | {100*r.coverage_fraction:.2f}% | {r.n_detected} | {r.n_scored} | {r.cohort_status} | {r.family_size} |"
        )
    lines.extend([
        "",
        "## Verdict rule frozen before running",
        "",
        "A signature SURVIVES if at least one unit passes both tests after batch adjustment in either testable primary cohort. It is LOST_WITH_BATCH only if it has at least one unadjusted both-test pass in each primary cohort and zero adjusted both-test passes in both cohorts. Other outcomes are recorded as NO_MAJORITY_CATEGORY or PARTIALLY_UNTESTABLE and are not mislabeled as batch losses. PROVES requires a majority of the four benchmarked signatures to be LOST_WITH_BATCH. DOES NOT PROVE requires a majority to SURVIVE. INCONCLUSIVE applies if most signatures are coverage-untestable, results split evenly, or neither majority rule is met. The verdict names every surviving and non-surviving signature. Four signatures cannot prove a claim about ALL published signatures; the verdict is about this frozen benchmarked set only.",
        "",
        "## Results",
        "",
        "Inferential results are appended only after this frozen pre-analysis section has been written.",
        "",
    ])
    return "\n".join(lines)


def result_row(signature_id: str, cohort: str, arm: str, source_row: pd.Series) -> dict[str, object]:
    se = float(source_row["composite_se"])
    df = int(source_row["residual_df"])
    mde80 = float((stats.t.ppf(0.975, df) + stats.t.ppf(0.80, df)) * se)
    covariates = str(source_row["covariates"])
    if cohort == "Yazar_Blood":
        n_pool = covariates.count("pool_")
        covariates = ";".join(c for c in covariates.split(";") if not c.startswith("pool_"))
        if arm == "batch_adjusted":
            covariates += f";pool_fe({n_pool}_dummies)"
    return {
        "signature_id": signature_id,
        "cohort": cohort,
        "arm": arm,
        "unit": source_row["group"],
        "status": "TESTED",
        "composite_slope": float(source_row["composite_effect"]),
        "composite_slope_per_decade": float(source_row["composite_effect"]),
        "composite_se": se,
        "composite_ci95_low": float(source_row["composite_ci95_low"]),
        "composite_ci95_high": float(source_row["composite_ci95_high"]),
        "composite_pvalue": float(source_row["composite_pvalue"]),
        "composite_q": np.nan,
        "empirical_pvalue": float(source_row["empirical_null_pvalue"]),
        "empirical_q": np.nan,
        "composite_pass": False,
        "empirical_pass": False,
        "both_pass": False,
        "n_donors": int(source_row["n_donors"]),
        "n_genes_published": np.nan,
        "n_genes_mapped": np.nan,
        "n_genes_detected": int(source_row["n_p1_total_detected"]),
        "n_up_detected": int(source_row["n_p1_up_detected"]),
        "n_down_detected": int(source_row["n_p1_down_detected"]),
        "n_genes_scored": int(source_row["n_p1_total_detected"]),
        "mde80": mde80,
        "null_mean_effect": float(source_row["null_mean_effect"]),
        "null_sd_effect": float(source_row["null_sd_effect"]),
        "null_iterations": int(source_row["null_iterations"]),
        "random_seed": SEED,
        "residual_df": df,
        "covariates": covariates,
        "omitted_covariates": source_row["omitted_covariates"],
        "multiplicity_family": f"{signature_id}_{cohort}_{arm}",
        "composite_family_size": np.nan,
        "empirical_family_size": np.nan,
    }


def untestable_rows(registry: pd.DataFrame, preflight: pd.DataFrame) -> list[dict[str, object]]:
    rows = []
    for r in preflight[preflight["cohort_status"] != "TESTABLE"].itertuples(index=False):
        for arm in ARMS:
            rows.append({
                "signature_id": r.signature_id,
                "cohort": r.cohort,
                "arm": arm,
                "unit": "COHORT_LEVEL",
                "status": "UNTESTABLE_COVERAGE",
                "composite_slope": np.nan,
                "composite_slope_per_decade": np.nan,
                "composite_se": np.nan,
                "composite_ci95_low": np.nan,
                "composite_ci95_high": np.nan,
                "composite_pvalue": np.nan,
                "composite_q": np.nan,
                "empirical_pvalue": np.nan,
                "empirical_q": np.nan,
                "composite_pass": False,
                "empirical_pass": False,
                "both_pass": False,
                "n_donors": np.nan,
                "n_genes_published": np.nan,
                "n_genes_mapped": np.nan,
                "n_genes_detected": r.n_detected,
                "n_up_detected": np.nan,
                "n_down_detected": np.nan,
                "n_genes_scored": 0,
                "mde80": np.nan,
                "null_mean_effect": np.nan,
                "null_sd_effect": np.nan,
                "null_iterations": 0,
                "random_seed": SEED,
                "residual_df": np.nan,
                "covariates": "not_fit_coverage_below_50_percent",
                "omitted_covariates": "not_applicable",
                "multiplicity_family": f"{r.signature_id}_{r.cohort}_{arm}",
                "composite_family_size": 0,
                "empirical_family_size": 0,
            })
    return rows


def run_inference(signatures, cohort_inputs, registry: pd.DataFrame, preflight: pd.DataFrame, details: dict) -> pd.DataFrame:
    rows = untestable_rows(registry, preflight)
    for signature_id in SIGNATURE_ORDER:
        program = {"up": list(signatures[signature_id]["up"]), "down": list(signatures[signature_id]["down"])}
        for cohort in COHORTS:
            info = details[signature_id][cohort]
            if info["status"] != "TESTABLE":
                continue
            units = set(info["units"])
            for arm in ARMS:
                expr, metas, covariates, min_residual_df = cohort_inputs[cohort][arm]
                expr_sub = {group: expr[group] for group in sorted(units)}
                meta_sub = {group: metas[group] for group in sorted(units)}
                print(f"Running {signature_id} {cohort} {arm}: {len(units)} units", flush=True)
                fitted = run_family_c_p1_signature(
                    expr_sub,
                    meta_sub,
                    p1_program=program,
                    age_col="age_per_decade",
                    candidate_covariates=covariates,
                    min_donors=MIN_DONORS,
                    min_residual_df=min_residual_df,
                    iterations=ITERATIONS,
                    seed=SEED,
                )
                assert len(fitted) == len(units), (signature_id, cohort, arm, len(fitted), len(units))
                assert set(fitted["group"]) == units
                for _, source_row in fitted.iterrows():
                    rows.append(result_row(signature_id, cohort, arm, source_row))
    out = pd.DataFrame(rows)
    published_n = registry.set_index("signature_id")["n_genes_published"].to_dict()
    mapped_n = registry.set_index("signature_id")["n_unique_symbols_used"].to_dict()
    out["n_genes_published"] = out["signature_id"].map(published_n).astype(int)
    out["n_genes_mapped"] = out["signature_id"].map(mapped_n).astype(int)
    for (_, _, _), indices in out[out["status"] == "TESTED"].groupby(["signature_id", "cohort", "arm"]).groups.items():
        n = len(indices)
        out.loc[indices, "composite_q"] = multipletests(out.loc[indices, "composite_pvalue"], method="fdr_bh")[1]
        out.loc[indices, "empirical_q"] = multipletests(out.loc[indices, "empirical_pvalue"], method="fdr_bh")[1]
        out.loc[indices, "composite_family_size"] = n
        out.loc[indices, "empirical_family_size"] = n
    out["composite_pass"] = out["composite_q"].lt(0.05).fillna(False)
    out["empirical_pass"] = out["empirical_q"].lt(0.05).fillna(False)
    out["both_pass"] = out["composite_pass"] & out["empirical_pass"]
    out = out.sort_values(["signature_id", "cohort", "arm", "unit"]).reset_index(drop=True)
    return out


def fmt(value) -> str:
    if pd.isna(value):
        return "NA"
    return f"{float(value):.4g}"


def append_results(report_text: str, results: pd.DataFrame, preflight: pd.DataFrame, original: pd.DataFrame) -> tuple[str, str]:
    lines = [
        "The pre-analysis section above was written before model fitting. This section contains the completed batch contrast.",
        "",
        "### Batch-contrast summary",
        "",
        "Both pass means composite BH q < 0.05 and empirical BH q < 0.05 in the frozen slice. N is the number of tested units.",
        "",
        "| Signature | Cohort | Arm | N units | Composite passes | Empirical passes | Both pass |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: |",
    ]
    for signature_id in SIGNATURE_ORDER:
        for cohort in COHORTS:
            for arm in ARMS:
                sub = results[(results["signature_id"] == signature_id) & (results["cohort"] == cohort) & (results["arm"] == arm)]
                tested = sub[sub["status"] == "TESTED"]
                lines.append(
                    f"| {signature_id} | {cohort} | {arm} | {len(tested)} | "
                    f"{int(tested['composite_pass'].sum())} | {int(tested['empirical_pass'].sum())} | {int(tested['both_pass'].sum())} |"
                )
    classifications = {}
    lines.extend([
        "",
        "### Per-signature outcome",
        "",
        "| Signature | Unadjusted both-pass units, Gao | Adjusted, Gao | Unadjusted both-pass units, Yazar | Adjusted, Yazar | Classification |",
        "| --- | ---: | ---: | ---: | ---: | --- |",
    ])
    for signature_id in SIGNATURE_ORDER:
        counts = {}
        statuses = {}
        for cohort in COHORTS:
            sub = results[(results["signature_id"] == signature_id) & (results["cohort"] == cohort)]
            statuses[cohort] = "UNTESTABLE" if (sub["status"] == "UNTESTABLE_COVERAGE").any() else "TESTABLE"
            for arm in ARMS:
                counts[(cohort, arm)] = int(sub[sub["arm"] == arm]["both_pass"].sum())
        if "UNTESTABLE" in statuses.values() and not any(counts[(c, "batch_adjusted")] > 0 for c in COHORTS):
            classification = "PARTIALLY_UNTESTABLE"
        elif any(counts[(c, "batch_adjusted")] > 0 for c in COHORTS):
            classification = "SURVIVES"
        elif all(statuses[c] == "TESTABLE" and counts[(c, "unadjusted")] > 0 and counts[(c, "batch_adjusted")] == 0 for c in COHORTS):
            classification = "LOST_WITH_BATCH"
        else:
            classification = "NO_MAJORITY_CATEGORY"
        classifications[signature_id] = classification
        lines.append(
            f"| {signature_id} | {counts[('Gao_Heart', 'unadjusted')]} | {counts[('Gao_Heart', 'batch_adjusted')]} | "
            f"{counts[('Yazar_Blood', 'unadjusted')]} | {counts[('Yazar_Blood', 'batch_adjusted')]} | {classification} |"
        )
    lines.extend([
        "",
        "### Mature vs original: survived vs died vs untestable",
        "",
        "Original counts are both-pass counts from the full-range benchmark (`results/followup-11pt-block10-benchmark-results.csv`, 246 tested rows). Mature counts are both-pass counts at 30-plus from this rerun. Died means the signature SURVIVED originally and has zero adjusted both-pass units at 30-plus in both cohorts while remaining testable. Untestable means coverage below 50 percent at 30-plus.",
        "",
        "| Signature | Original Gao unadj | Original Gao adj | Original Yazar unadj | Original Yazar adj | Original class | Mature Gao unadj | Mature Gao adj | Mature Yazar unadj | Mature Yazar adj | Mature class | Change |",
        "| --- | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | --- | --- |",
    ])
    for signature_id in SIGNATURE_ORDER:
        orig_counts = {}
        orig_statuses = {}
        for cohort in COHORTS:
            osub = original[(original["signature_id"] == signature_id) & (original["cohort"] == cohort)]
            orig_statuses[cohort] = "UNTESTABLE" if (osub["status"] == "UNTESTABLE_COVERAGE").any() else "TESTABLE"
            for arm in ARMS:
                orig_counts[(cohort, arm)] = int(osub[osub["arm"] == arm]["both_pass"].sum())
        if "UNTESTABLE" in orig_statuses.values() and not any(orig_counts[(c, "batch_adjusted")] > 0 for c in COHORTS):
            orig_class = "PARTIALLY_UNTESTABLE"
        elif any(orig_counts[(c, "batch_adjusted")] > 0 for c in COHORTS):
            orig_class = "SURVIVES"
        elif all(orig_statuses[c] == "TESTABLE" and orig_counts[(c, "unadjusted")] > 0 and orig_counts[(c, "batch_adjusted")] == 0 for c in COHORTS):
            orig_class = "LOST_WITH_BATCH"
        else:
            orig_class = "NO_MAJORITY_CATEGORY"
        mat_counts = {}
        for cohort in COHORTS:
            for arm in ARMS:
                mat_counts[(cohort, arm)] = int(results[(results["signature_id"] == signature_id) & (results["cohort"] == cohort) & (results["arm"] == arm)]["both_pass"].sum())
        mat_class = classifications[signature_id]
        if orig_class == mat_class:
            change = "unchanged"
        elif orig_class == "SURVIVES" and mat_class != "SURVIVES":
            change = "no longer survives at 30-plus"
        elif orig_class != "SURVIVES" and mat_class == "SURVIVES":
            change = "now survives at 30-plus"
        else:
            change = "category changed"
        lines.append(
            f"| {signature_id} | {orig_counts[('Gao_Heart', 'unadjusted')]} | {orig_counts[('Gao_Heart', 'batch_adjusted')]} | "
            f"{orig_counts[('Yazar_Blood', 'unadjusted')]} | {orig_counts[('Yazar_Blood', 'batch_adjusted')]} | {orig_class} | "
            f"{mat_counts[('Gao_Heart', 'unadjusted')]} | {mat_counts[('Gao_Heart', 'batch_adjusted')]} | "
            f"{mat_counts[('Yazar_Blood', 'unadjusted')]} | {mat_counts[('Yazar_Blood', 'batch_adjusted')]} | {mat_class} | {change} |"
        )
    n_lost = sum(v == "LOST_WITH_BATCH" for v in classifications.values())
    n_survive = sum(v == "SURVIVES" for v in classifications.values())
    n_untestable = sum(v == "PARTIALLY_UNTESTABLE" for v in classifications.values())
    if n_lost > len(SIGNATURE_ORDER) / 2:
        verdict = "PROVES"
        reason = "A majority of benchmarked signatures lost all passing units after adjustment in both primary cohorts."
    elif n_survive > len(SIGNATURE_ORDER) / 2:
        verdict = "DOES NOT PROVE"
        reason = "A majority of benchmarked signatures kept at least one passing unit after adjustment in a primary cohort."
    else:
        verdict = "INCONCLUSIVE"
        reason = "Neither majority rule was met for the four frozen signatures."
    survivors = [s for s, v in classifications.items() if v == "SURVIVES"]
    nonsurvivors = [s for s, v in classifications.items() if v != "SURVIVES"]
    lines.extend([
        "",
        "### Units passing both tests after batch adjustment",
        "",
        "| Signature | Cohort | Unit | Slope per decade | Composite q | Empirical q | N donors | MDE80 |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |",
    ])
    passing = results[(results["arm"] == "batch_adjusted") & results["both_pass"]]
    if passing.empty:
        lines.append("| None | None | None | NA | NA | NA | NA | NA |")
    for r in passing.itertuples(index=False):
        lines.append(
            f"| {r.signature_id} | {r.cohort} | {r.unit} | {fmt(r.composite_slope)} | "
            f"{fmt(r.composite_q)} | {fmt(r.empirical_q)} | {r.n_donors} | {fmt(r.mde80)} |"
        )
    lines.extend([
        "",
        "### Complete per-unit results and MDE80",
        "",
        "Every tested unit is shown so each null has its MDE80 alongside it. UNTESTABLE coverage slices have no model or MDE80 and remain listed in the frozen coverage table.",
        "",
        "| Signature | Cohort | Arm | Unit | Slope per decade | Composite q | Empirical q | Both pass | N donors | MDE80 |",
        "| --- | --- | --- | --- | ---: | ---: | ---: | --- | ---: | ---: |",
    ])
    for r in results[results["status"] == "TESTED"].itertuples(index=False):
        lines.append(
            f"| {r.signature_id} | {r.cohort} | {r.arm} | {r.unit} | {fmt(r.composite_slope)} | "
            f"{fmt(r.composite_q)} | {fmt(r.empirical_q)} | {bool(r.both_pass)} | {r.n_donors} | {fmt(r.mde80)} |"
        )
    lines.extend([
        "",
        "## Verdict",
        "",
        f"**{verdict}**. {reason}",
        "",
        "Surviving signatures: " + (", ".join(survivors) if survivors else "none") + ".",
        "",
        "Non-surviving signatures and their categories: " + "; ".join(f"{s} = {classifications[s]}" for s in nonsurvivors) + ".",
        "",
        f"Counts: {n_survive} SURVIVES, {n_lost} LOST_WITH_BATCH, {n_untestable} PARTIALLY_UNTESTABLE, and {len(SIGNATURE_ORDER) - n_survive - n_lost - n_untestable} NO_MAJORITY_CATEGORY.",
        "",
        "These four signatures cannot prove a claim about ALL published signatures. This verdict is only about the frozen benchmarked set and its tested Gao and Yazar units. Absence of a passing unit is not proof of no biological age association; per-unit MDE80 values are supplied in the result CSV so nulls can be judged for power.",
        "",
        "## Reproducibility and outputs",
        "",
        "- `results/followup-mature-block10-benchmark-signatures.csv`: frozen citations, sources, exact up and down lists, coverage at 30-plus, and family sizes at 30-plus.",
        "- `results/followup-mature-block10-benchmark-results.csv`: one row per signature, cohort, arm, and tested unit at 30-plus, plus cohort-level UNTESTABLE rows where required.",
        "- `results/followup-mature-block10-donor-counts.csv`: donor counts at full range and at 30-plus for every cohort, arm, and unit.",
        "- `scripts/run_mature_block10_published_benchmark.py`: mature restriction, preflight freeze, Block 5 transport, matched nulls recomputed at 30-plus, BH correction, and report generation.",
        "- Fixed seed: 20261105. Random iterations: 2000 per tested signature, cohort, arm, and unit, recomputed on the 30-plus donors.",
        "- Original full-range benchmark for comparison: `results/followup-11pt-block10-benchmark-results.csv` and `reports/followup-11pt-block10-benchmark.md`.",
        "",
    ])
    text = report_text.replace(
        "Inferential results are appended only after this frozen pre-analysis section has been written.\n",
        "\n".join(lines),
    )
    return text, verdict


def validate_outputs(registry: pd.DataFrame, results: pd.DataFrame, preflight: pd.DataFrame) -> None:
    assert len(registry) == 4
    assert set(registry["signature_id"]) == set(SIGNATURE_ORDER)
    assert (registry["retrieval_status"] == "RETRIEVED_AND_FROZEN_BEFORE_INFERENCE").all()
    tested = results[results["status"] == "TESTED"]
    assert (tested["n_donors"] >= MIN_DONORS).all()
    assert (tested["null_iterations"] == ITERATIONS).all()
    assert (tested["random_seed"] == SEED).all()
    assert np.isfinite(tested[["composite_slope", "composite_pvalue", "empirical_pvalue", "mde80"]]).all().all()
    for r in preflight.itertuples(index=False):
        sub = tested[(tested["signature_id"] == r.signature_id) & (tested["cohort"] == r.cohort)]
        for arm in ARMS:
            arm_sub = sub[sub["arm"] == arm]
            assert len(arm_sub) == r.family_size
            if r.family_size:
                assert (arm_sub["composite_family_size"] == r.family_size).all()
                assert (arm_sub["empirical_family_size"] == r.family_size).all()
                assert np.allclose(
                    arm_sub["composite_q"],
                    multipletests(arm_sub["composite_pvalue"], method="fdr_bh")[1],
                )
                assert np.allclose(
                    arm_sub["empirical_q"],
                    multipletests(arm_sub["empirical_pvalue"], method="fdr_bh")[1],
                )
    assert not results.duplicated(["signature_id", "cohort", "arm", "unit"]).any()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--preflight", action="store_true", help="Freeze registry and family sizes without fitting models")
    args = parser.parse_args()
    signatures = load_frozen_signatures()
    cohort_inputs = {
        "Gao_Heart": load_gao_inputs(),
        "Yazar_Blood": load_yazar_inputs(),
    }
    cohort_inputs, donor_counts, donor_totals = restrict_to_mature(cohort_inputs)
    write_clean_csv(donor_counts, DONOR_COUNTS_CSV)
    print(f"Wrote {DONOR_COUNTS_CSV} ({len(donor_counts)} rows)", flush=True)
    preflight, details = build_preflight(signatures, cohort_inputs)
    registry = registry_frame(signatures, preflight)
    frozen_text = frozen_report_text(registry, preflight, donor_totals)
    write_clean_csv(registry, REGISTRY_CSV)
    write_clean_markdown(frozen_text, REPORT)
    print(f"Preflight freeze written to {REPORT}", flush=True)
    if args.preflight:
        write_clean_csv(pd.DataFrame(columns=["signature_id"]), RESULTS_CSV)
        return 0
    results = run_inference(signatures, cohort_inputs, registry, preflight, details)
    validate_outputs(registry, results, preflight)
    write_clean_csv(results, RESULTS_CSV)
    original = pd.read_csv(ORIGINAL_RESULTS_CSV)
    final_text, verdict = append_results(frozen_text, results, preflight, original)
    write_clean_markdown(final_text, REPORT)
    print(f"Mature Block 10 verdict: {verdict}", flush=True)
    print(f"Wrote {REGISTRY_CSV} ({len(registry)} rows)", flush=True)
    print(f"Wrote {RESULTS_CSV} ({len(results)} rows)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
