#!/usr/bin/env python3
"""Step 9: SEA-AD MTG donor-pseudobulk interneuron gene replication.

Only SEA-AD expression is read. Existing MSSM/HBCC aggregate effects and Rev2
membership tables are inputs; no existing cohort expression is refitted. Raw
per-donor caches stay outside the repository. Public result CSVs contain only
gene-level or cohort/subtype aggregates, not donor-level expression values.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from src.aging_programs import load_program_tiers
from scripts.run_rev2_de_interneuron import (
    coexpression_blocks,
    filter_genes,
    lodo_sensitivity,
    ols_arm,
    voom_wls_counts,
)

RESULTS = REPO / "results"
MANIFEST = REPO / "data" / "manifest"
PLAN_JSON = MANIFEST / "sea_ad_interneuron_analysis_plan.json"
PLAN_MD = REPO / "notes" / "followup-sea-ad-interneurons-plan.md"
WORKBOOK_NAME = "sea-ad_cohort_donor_metadata_072524.xlsx"
WORKBOOK_URL = (
    "https://cdn.prod.website-files.com/689cfbd308fa7373b604d290/"
    "68debdfdd1b8e9f8fd64dab0_sea-ad_cohort_donor_metadata_072524.xlsx"
)
OBJECT_BASE = (
    "https://sea-ad-single-cell-profiling.s3.us-west-2.amazonaws.com/"
    "MTG/RNAseq/previous_objects/donor_objects"
)
RELEASE = "2024-02-13"
SUBCLASS = {"VIP": "Vip", "pvalb": "Pvalb", "sst": "Sst"}
CELL_TYPE = {
    "VIP": "VIP GABAergic cortical interneuron",
    "pvalb": "pvalb GABAergic cortical interneuron",
    "sst": "sst GABAergic cortical interneuron",
}
PRIMARY_COLUMNS = ("age", "sex_male", "mean_log_umi", "multiome_fraction", "adnc_low")
FAMILY_ORDER = (
    "psychad_hbcc_both_fdr",
    "psychad_discovery_all",
    "psychad_without_hbcc_both_fdr",
    "mouse_p1_all",
    "mouse_p1_psychad_overlap",
    "mouse_p1_outside_psychad",
)


def log(*args):
    print(*args, flush=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_csv(frame: pd.DataFrame, path: Path):
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False)
    temporary.replace(path)
    log(f"wrote {path.relative_to(REPO)}: {len(frame)} rows")


def load_plan() -> dict:
    plan = json.loads(PLAN_JSON.read_text())
    assert sha256(PLAN_MD) == plan["frozen_plan_sha256"], "Frozen plan hash changed"
    assert sha256(REPO / plan["targeted_set_source"]["path"]) == plan["targeted_set_source"]["sha256"]
    assert sha256(REPO / plan["ortholog_source"]["path"]) == plan["ortholog_source"]["sha256"]
    return plan


def eligible_donors(workbook: Path) -> pd.DataFrame:
    donors = pd.read_excel(workbook, sheet_name="SEA-AD_Cohort_Metadata")
    required = {
        "Donor ID", "Age at Death", "Sex", "Cognitive Status",
        "Consensus Clinical Dx (choice=Control)",
        "Overall AD neuropathological Change",
    }
    assert required.issubset(donors.columns)
    assert donors["Donor ID"].is_unique
    selected = donors[
        (donors["Cognitive Status"] == "No dementia")
        & (donors["Consensus Clinical Dx (choice=Control)"] == "Checked")
        & donors["Overall AD neuropathological Change"].isin(["Not AD", "Low"])
    ].copy()
    assert len(selected) == 15, "The frozen metadata gate no longer yields 15 donors"
    assert selected["Age at Death"].between(72, 102).all()
    return selected.sort_values("Donor ID").reset_index(drop=True)


def object_name(donor: str) -> str:
    return f"{donor}_SEAAD_MTG_RNAseq_final-nuclei.{RELEASE}.h5ad"


def acquire(data_root: Path):
    metadata_dir = data_root / "metadata"
    object_dir = data_root / "donor-h5ad"
    metadata_dir.mkdir(parents=True, exist_ok=True)
    object_dir.mkdir(parents=True, exist_ok=True)
    guard = [sys.executable, str(REPO / "scripts" / "resource_guard.py"), "download"]
    workbook = metadata_dir / WORKBOOK_NAME
    subprocess.run(
        guard + ["--url", WORKBOOK_URL, "--output", str(workbook),
                 "--min-free-gb", "30", "--max-file-gb", "0.1"],
        check=True,
    )
    for donor in eligible_donors(workbook)["Donor ID"]:
        name = object_name(donor)
        subprocess.run(
            guard + ["--url", f"{OBJECT_BASE}/{name}", "--output", str(object_dir / name),
                     "--min-free-gb", "30", "--max-file-gb", "2"],
            check=True,
        )


def sparse_counts(matrix) -> sp.csr_matrix:
    matrix = matrix.tocsr() if sp.issparse(matrix) else sp.csr_matrix(matrix)
    matrix = matrix.astype(np.float64)
    assert np.isfinite(matrix.data).all()
    assert (matrix.data >= 0).all()
    assert np.allclose(matrix.data, np.round(matrix.data)), "UMIs layer is not raw integer counts"
    return matrix


def summarize_counts(matrix: sp.csr_matrix):
    library = np.asarray(matrix.sum(axis=1)).ravel()
    return (
        np.asarray(matrix.sum(axis=0)).ravel().astype(np.int64),
        float(np.log(library).mean()),
        float(library.sum()),
    )


def build_donor_cache(path: Path, donor_row: pd.Series, output: Path):
    import anndata as ad

    donor = str(donor_row["Donor ID"])
    adata = ad.read_h5ad(path, backed="r")
    obs = adata.obs
    assert adata.n_vars == 36601
    assert adata.var_names.is_unique
    genes = np.asarray(adata.var_names.astype(str), dtype=str)
    assert obs["Donor ID"].astype(str).eq(donor).all()
    assert set(obs["Brain Region"].astype(str)).issubset({"Human MTG", "Human MTG All Layers"})
    assert obs["Age at Death"].nunique() == 1
    assert float(obs["Age at Death"].iloc[0]) == float(donor_row["Age at Death"])
    assert obs["Sex"].astype(str).eq(str(donor_row["Sex"])).all()
    assert obs["Cognitive Status"].astype(str).eq("No dementia").all()
    assert obs["Overall AD neuropathological Change"].astype(str).eq(
        str(donor_row["Overall AD neuropathological Change"])
    ).all()
    assert set(obs["method"].astype(str)).issubset({"10Xv3.1", "10xMulti"})
    qc = (
        obs["Used in analysis"].astype(bool)
        & obs["rna_amplification_pass_fail"].astype(str).eq("Pass")
        & obs["library_prep_pass_fail"].astype(str).eq("Pass")
        & ~obs["experiment_component_failed"].astype(bool)
        & pd.to_numeric(obs["Number of UMIs"], errors="coerce").ge(1500)
    )
    multiome = obs["method"].astype(str).str.lower().str.contains("multi").to_numpy()
    arrays = {"genes": genes, "donor": np.array(donor)}
    qc_row = {
        "donor": donor,
        "n_released_nuclei": int(len(obs)),
        "n_qc_pass_nuclei": int(qc.sum()),
        "n_multiome_nuclei": int(multiome.sum()),
        "all_used_in_analysis": bool(obs["Used in analysis"].all()),
        "all_preparation_flags_pass": bool(
            obs["rna_amplification_pass_fail"].astype(str).eq("Pass").all()
            and obs["library_prep_pass_fail"].astype(str).eq("Pass").all()
            and not obs["experiment_component_failed"].any()
        ),
    }
    for short, label in SUBCLASS.items():
        subclass = obs["Subclass"].astype(str).eq(label)
        if subclass.any():
            assert obs.loc[subclass, "Class"].astype(str).eq("Neuronal: GABAergic").all()
        selected = (subclass & qc).to_numpy()
        singleome = selected & ~multiome
        qc_row[f"{short}_nuclei"] = int(selected.sum())
        qc_row[f"{short}_singleome_nuclei"] = int(singleome.sum())
        for arm, mask in (("all", selected), ("singleome", singleome)):
            if mask.sum() == 0:
                raise ValueError(f"No {arm} nuclei for {donor} {short}")
            matrix = sparse_counts(adata.layers["UMIs"][np.flatnonzero(mask)])
            observed_library = np.asarray(matrix.sum(axis=1)).ravel()
            released_library = pd.to_numeric(
                obs.loc[mask, "Number of UMIs"], errors="raise"
            ).to_numpy(dtype=float)
            assert np.allclose(observed_library, released_library)
            counts, mean_log_umi, library = summarize_counts(matrix)
            prefix = short if arm == "all" else f"{short}_singleome"
            arrays[f"{prefix}_counts"] = counts
            arrays[f"{prefix}_n_cells"] = np.array(int(mask.sum()), dtype=np.int64)
            arrays[f"{prefix}_mean_log_umi"] = np.array(mean_log_umi)
            arrays[f"{prefix}_library"] = np.array(library)
            del matrix
        arrays[f"{short}_multiome_fraction"] = np.array(
            float(multiome[selected].mean()) if selected.sum() else 0.0
        )
    arrays["age"] = np.array(float(donor_row["Age at Death"]))
    arrays["sex_male"] = np.array(float(donor_row["Sex"] == "Male"))
    arrays["adnc_low"] = np.array(
        float(donor_row["Overall AD neuropathological Change"] == "Low")
    )
    temporary = output.with_suffix(".tmp.npz")
    np.savez_compressed(temporary, **arrays)
    temporary.replace(output)
    del adata, obs
    gc.collect()
    return qc_row


def build_metadata(data_root: Path, force: bool = False):
    workbook = data_root / "metadata" / WORKBOOK_NAME
    donors = eligible_donors(workbook)
    interim = data_root / "interim" / "donor-pseudobulk"
    interim.mkdir(parents=True, exist_ok=True)
    qc_rows = []
    external_root = Path("../agescape-revision-data/sea-ad")
    source_files = [{
        "kind": "donor_metadata",
        "url": WORKBOOK_URL,
        "local_path": str(external_root / "metadata" / WORKBOOK_NAME),
        "size_bytes": workbook.stat().st_size,
        "sha256": sha256(workbook),
    }]
    for _, row in donors.iterrows():
        donor = row["Donor ID"]
        path = data_root / "donor-h5ad" / object_name(donor)
        assert path.exists(), f"Missing donor object: {path}"
        cache = interim / f"{donor}.npz"
        if force or not cache.exists():
            qc_row = build_donor_cache(path, row, cache)
        else:
            qc_row = {
                "donor": donor,
                "n_released_nuclei": None,
                "n_qc_pass_nuclei": None,
                "n_multiome_nuclei": None,
                "all_used_in_analysis": None,
                "all_preparation_flags_pass": None,
            }
            cached = np.load(cache, allow_pickle=False)
            for short in SUBCLASS:
                qc_row[f"{short}_nuclei"] = int(cached[f"{short}_n_cells"])
                qc_row[f"{short}_singleome_nuclei"] = int(cached[f"{short}_singleome_n_cells"])
        qc_rows.append(qc_row)
        source_files.append({
            "kind": "donor_h5ad",
            "donor": donor,
            "url": f"{OBJECT_BASE}/{object_name(donor)}",
            "local_path": str(external_root / "donor-h5ad" / object_name(donor)),
            "size_bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
        log(f"metadata gate {donor}: complete")
    qc = pd.DataFrame(qc_rows)
    qc.to_csv(data_root / "interim" / "donor-qc.csv", index=False)
    eligibility = []
    for short in SUBCLASS:
        nuclei = qc[f"{short}_nuclei"].astype(int)
        singleome = qc[f"{short}_singleome_nuclei"].astype(int)
        ages = donors.loc[nuclei >= 20, "Age at Death"].astype(float)
        sex = donors.loc[nuclei >= 20, "Sex"]
        eligibility.append({
            "subtype": short,
            "metadata_eligible_donors": len(donors),
            "eligible_donors_20_nuclei": int((nuclei >= 20).sum()),
            "eligible_donors_50_nuclei": int((nuclei >= 50).sum()),
            "singleome_eligible_donors_20_nuclei": int((singleome >= 20).sum()),
            "excluded_donors_below_20": int((nuclei < 20).sum()),
            "total_eligible_nuclei": int(nuclei[nuclei >= 20].sum()),
            "median_nuclei_per_donor": float(nuclei[nuclei >= 20].median()),
            "minimum_nuclei_per_donor": int(nuclei[nuclei >= 20].min()),
            "maximum_nuclei_per_donor": int(nuclei[nuclei >= 20].max()),
            "age_min": float(ages.min()),
            "age_max": float(ages.max()),
            "age_sd": float(ages.std()),
            "female_donors": int((sex == "Female").sum()),
            "male_donors": int((sex == "Male").sum()),
            "not_ad_donors": int((donors["Overall AD neuropathological Change"] == "Not AD").sum()),
            "low_adnc_donors": int((donors["Overall AD neuropathological Change"] == "Low").sum()),
            "multiome_contributing_donors": int(sum(
                np.load(interim / f"{d}.npz", allow_pickle=False)[f"{short}_multiome_fraction"] > 0
                for d in donors["Donor ID"]
            )),
            "annotation_mapping": f"{SUBCLASS[short]}; Class=Neuronal: GABAergic",
            "raw_count_layer": "UMIs",
            "qc_rule": "released final nuclei, preparation flags pass, Number of UMIs >=1500",
        })
    atomic_csv(pd.DataFrame(eligibility), RESULTS / "followup-sea-ad-eligibility.csv")
    manifest = {
        "accession": "SEA-AD MTG RNAseq final nuclei 2024-02-13",
        "paper": {
            "title": "Integrated multimodal cell atlas of Alzheimer's disease",
            "doi": "10.1038/s41593-024-01774-5",
        },
        "data_page": "https://brain-map.org/consortia/sea-ad/our-data",
        "bucket_inventory": "https://sea-ad-single-cell-profiling.s3.us-west-2.amazonaws.com/?list-type=2&prefix=MTG/RNAseq/previous_objects/",
        "created_from_frozen_plan": "data/manifest/sea_ad_interneuron_analysis_plan.json",
        "release": RELEASE,
        "region": "MTG",
        "assay": {"singleome": "10Xv3.1", "multiome": "10xMulti"},
        "raw_count_layer": "UMIs",
        "normalized_X_used": False,
        "donor_filter": "No dementia; consensus Control; ADNC Not AD or Low",
        "donor_count": len(donors),
        "annotation_mapping": SUBCLASS,
        "sources": source_files,
        "resource_guard": {
            "guard": "scripts/resource_guard.py for acquisition and expression jobs",
            "minimum_free_gb": 30,
            "analysis_memory_gb": 8,
            "maximum_single_donor_file_gb": 2,
            "parallel_heavy_jobs": 1,
            "external_data_root": str(external_root),
        },
        "public_outputs_contain_donor_expression": False,
    }
    (MANIFEST / "sea-ad-mtg-2024-02-13.json").write_text(json.dumps(manifest, indent=2) + "\n")
    log(f"wrote data/manifest/sea-ad-mtg-2024-02-13.json with {len(source_files)} verified sources")


def load_counts(data_root: Path, short: str, arm: str = "all", floor: int = 20,
                not_ad_only: bool = False):
    interim = data_root / "interim" / "donor-pseudobulk"
    rows, matrices = [], []
    genes = None
    for path in sorted(interim.glob("*.npz")):
        cached = np.load(path, allow_pickle=False)
        prefix = short if arm == "all" else f"{short}_singleome"
        n_cells = int(cached[f"{prefix}_n_cells"])
        if n_cells < floor:
            continue
        if not_ad_only and int(cached["adnc_low"]) == 1:
            continue
        if genes is None:
            genes = cached["genes"].astype(str)
        else:
            assert np.array_equal(genes, cached["genes"].astype(str))
        matrices.append(cached[f"{prefix}_counts"])
        rows.append({
            "donor": str(cached["donor"]),
            "age": float(cached["age"]),
            "sex_male": float(cached["sex_male"]),
            "mean_log_umi": float(cached[f"{prefix}_mean_log_umi"]),
            "library": float(cached[f"{prefix}_library"]),
            "n_cells": n_cells,
            "multiome_fraction": 0.0 if arm == "singleome" else float(cached[f"{short}_multiome_fraction"]),
            "adnc_low": float(cached["adnc_low"]),
        })
    if not rows:
        return None, pd.DataFrame(), np.array([])
    return np.vstack(matrices), pd.DataFrame(rows), genes


def checked_design(metadata: pd.DataFrame, columns):
    columns = [c for c in columns if metadata[c].nunique() > 1]
    matrix = np.column_stack([np.ones(len(metadata))] + [metadata[c].astype(float) for c in columns])
    rank = np.linalg.matrix_rank(matrix)
    if rank < matrix.shape[1]:
        raise ValueError("Age design is rank-deficient; frozen rule marks this arm untestable")
    dof = len(metadata) - rank
    if len(metadata) < 10 or dof < 5:
        raise ValueError("Donor or residual-df floor failed")
    return columns, dof


def load_discovery(source: Path, short: str) -> pd.DataFrame:
    columns = [
        "cell_type", "model", "gene", "slope_log_cp10k_per_year_mssm",
        "pvalue_mssm", "q_bh_mssm", "slope_log_cp10k_per_year_hbcc",
        "pvalue_hbcc", "q_bh_hbcc",
    ]
    frame = pd.read_csv(source, usecols=columns)
    frame = frame[(frame["cell_type"] == CELL_TYPE[short]) & (frame["model"] == "age_plus_sex_depth")].copy()
    assert not frame["gene"].duplicated().any()
    return frame.drop(columns=["cell_type", "model"]).rename(columns={
        "slope_log_cp10k_per_year_mssm": "mssm_slope_original_ols",
        "pvalue_mssm": "mssm_p_original_ols",
        "q_bh_mssm": "mssm_q_original_ols",
        "slope_log_cp10k_per_year_hbcc": "hbcc_slope_original_ols",
        "pvalue_hbcc": "hbcc_p_original_ols",
        "q_bh_hbcc": "hbcc_q_original_ols",
    })


def load_targets(short: str):
    frame = pd.read_csv(RESULTS / "reviewer-rev2-interneuron-genes.csv")
    frame = frame[frame["cell_type"] == CELL_TYPE[short]].copy()
    assert len(frame) == {"VIP": 202, "pvalb": 143, "sst": 375}[short]
    assert not frame["gene"].duplicated().any()
    both = set(frame.loc[frame["tier"] == "both_FDR", "gene"])
    assert len(both) == {"VIP": 76, "pvalb": 66, "sst": 249}[short]
    directions = dict(zip(frame["gene"], np.sign(frame["slope_mssm"])))
    p1 = load_program_tiers(REPO, target="human")["P1_Clean"]
    p1_direction = {**{g: 1.0 for g in p1["up"]}, **{g: -1.0 for g in p1["down"]}}
    assert len(p1_direction) == 118
    discovery = set(frame["gene"])
    return frame, both, discovery, directions, p1_direction


def fit_wls(counts, metadata, genes, columns=PRIMARY_COLUMNS):
    active, dof = checked_design(metadata, list(columns))
    fitted = voom_wls_counts(counts, metadata, genes, columns=active)
    fitted.attrs["active_columns"] = active
    fitted.attrs["dof"] = dof
    return fitted


def family_membership(genes: pd.Index, short: str):
    target, both, discovery, directions, p1_direction = load_targets(short)
    rev2 = target.set_index("gene")
    membership = {
        "psychad_hbcc_both_fdr": genes.isin(both),
        "psychad_discovery_all": genes.isin(discovery),
        "psychad_without_hbcc_both_fdr": genes.isin(discovery - both),
        "mouse_p1_all": genes.isin(p1_direction),
        "mouse_p1_psychad_overlap": genes.isin(set(p1_direction) & discovery),
        "mouse_p1_outside_psychad": genes.isin(set(p1_direction) - discovery),
    }
    membership = {key: np.asarray(value, dtype=bool) for key, value in membership.items()}
    source_sign = {}
    source_size = {}
    for family in FAMILY_ORDER:
        if family.startswith("mouse_p1"):
            selected = set(p1_direction)
            if family == "mouse_p1_psychad_overlap":
                selected &= discovery
            elif family == "mouse_p1_outside_psychad":
                selected -= discovery
            source_sign[family] = np.array([p1_direction.get(g, 0.0) for g in genes])
            source_size[family] = len(selected)
        else:
            selected = discovery
            if family == "psychad_hbcc_both_fdr":
                selected = both
            elif family == "psychad_without_hbcc_both_fdr":
                selected = discovery - both
            source_sign[family] = np.array([directions.get(g, 0.0) for g in genes])
            source_size[family] = len(selected)
    return membership, source_sign, source_size, rev2


def correlation_rows(genes: pd.DataFrame, x_column: str, y_column: str):
    shared = genes[[x_column, y_column]].dropna()
    if len(shared) < 3:
        return len(shared), np.nan, np.nan
    return (
        len(shared),
        float(stats.pearsonr(shared[x_column], shared[y_column]).statistic),
        float(stats.spearmanr(shared[x_column], shared[y_column]).statistic),
    )


def targeted_summary(frame: pd.DataFrame, family: str, membership: np.ndarray,
                     signs: np.ndarray, total: int, slope="slope",
                     pvalue="pvalue", qvalue="q_bh"):
    selected = frame.loc[membership].copy()
    direction = signs[membership]
    same = selected[slope].to_numpy() * direction > 0
    return {
        "family": family,
        "n_source_genes": int(total),
        "n_tested": int(len(selected)),
        "n_missing_or_filtered": int(total - len(selected)),
        "same_direction": int(same.sum()),
        "same_direction_fraction": float(same.mean()) if len(same) else np.nan,
        "same_direction_nominal": int((same & (selected[pvalue].to_numpy() < 0.05)).sum()),
        "same_direction_seaad_fdr": int((same & (selected[qvalue].to_numpy() < 0.05)).sum()),
    }


def permute_within_blocks(values, blocks, rng):
    output = values.copy()
    for block in np.unique(blocks):
        idx = np.flatnonzero(blocks == block)
        output[idx] = values[rng.permutation(idx)]
    return output


def bootstrap_indices(blocks, rng):
    selected = rng.choice(np.unique(blocks), size=len(np.unique(blocks)), replace=True)
    return np.concatenate([np.flatnonzero(blocks == block) for block in selected])


def analyze_subtype(data_root: Path, discovery_source: Path, short: str,
                    permutations: int, bootstraps: int):
    counts, metadata, genes = load_counts(data_root, short)
    wls = fit_wls(counts, metadata, genes)
    ols = ols_arm(counts, metadata, genes, columns=wls.attrs["active_columns"])
    frame = wls.merge(ols, on="gene", how="left")
    discovery = load_discovery(discovery_source, short)
    frame = frame.merge(discovery, on="gene", how="left")
    membership, source_sign, source_size, rev2 = family_membership(frame["gene"], short)
    frame["in_psychad_hbcc_both_fdr"] = membership["psychad_hbcc_both_fdr"]
    frame["in_psychad_discovery"] = membership["psychad_discovery_all"]
    frame["in_mouse_p1"] = membership["mouse_p1_all"]
    for column in ("slope_mssm", "se_mssm", "ci95_low_mssm", "ci95_high_mssm",
                   "pvalue_mssm", "q_bh_mssm", "slope_hbcc", "q_bh_hbcc"):
        if column in rev2.columns:
            frame[f"rev2_{column}"] = frame["gene"].map(rev2[column])
    frame["subtype"] = short
    frame["model"] = "primary_voom_wls"
    frame["n_donors"] = len(metadata)
    frame["resid_dof"] = wls.attrs["dof"]
    frame["design"] = " + ".join(["intercept"] + wls.attrs["active_columns"])
    frame["slope_per_10_years"] = 10 * frame["slope"]
    frame["ci95_low_per_10_years"] = 10 * frame["ci95_low"]
    frame["ci95_high_per_10_years"] = 10 * frame["ci95_high"]

    shared = frame["mssm_slope_original_ols"].notna().to_numpy()
    shared_frame = frame.loc[shared].reset_index(drop=True)
    shared_genes = shared_frame["gene"]
    shared_membership = {key: value[shared] for key, value in membership.items()}
    shared_sign = {key: value[shared] for key, value in source_sign.items()}
    n_shared, pearson, spearman = correlation_rows(
        shared_frame, "mssm_slope_original_ols", "slope_ols"
    )
    _, pearson_wls, spearman_wls = correlation_rows(
        shared_frame, "mssm_slope_original_ols", "slope"
    )
    shared_positions = np.flatnonzero(shared)
    keep, _ = filter_genes(counts, genes)
    blocks = coexpression_blocks(counts[:, keep][:, shared_positions], shared_genes.to_numpy())
    rng = np.random.default_rng(20261108)
    null_correlation = np.empty(permutations)
    null_same = {family: np.empty(permutations) for family in FAMILY_ORDER}
    null_fdr = {family: np.empty(permutations) for family in FAMILY_ORDER}
    seaad_sign = np.sign(shared_frame["slope"].to_numpy())
    seaad_fdr = shared_frame["q_bh"].to_numpy() < 0.05
    seaad_ols = shared_frame["slope_ols"].to_numpy()
    mssm_ols = shared_frame["mssm_slope_original_ols"].to_numpy()
    for permutation in range(permutations):
        perm_slope = permute_within_blocks(seaad_ols, blocks, rng)
        perm_sign = permute_within_blocks(seaad_sign, blocks, rng)
        perm_fdr = permute_within_blocks(seaad_fdr, blocks, rng)
        null_correlation[permutation] = stats.pearsonr(mssm_ols, perm_slope).statistic
        for family in FAMILY_ORDER:
            selected = shared_membership[family]
            same = perm_sign[selected] * shared_sign[family][selected] > 0
            null_same[family][permutation] = same.sum()
            null_fdr[family][permutation] = (same & perm_fdr[selected]).sum()

    boot = np.random.default_rng(20261208)
    boot_correlation = np.empty(bootstraps)
    boot_spearman = np.empty(bootstraps)
    boot_same = {family: np.empty(bootstraps) for family in FAMILY_ORDER}
    observed_sign = np.sign(shared_frame["slope_ols"].to_numpy())
    for index in range(bootstraps):
        idx = bootstrap_indices(blocks, boot)
        boot_correlation[index] = stats.pearsonr(mssm_ols[idx], seaad_ols[idx]).statistic
        boot_spearman[index] = stats.spearmanr(mssm_ols[idx], seaad_ols[idx]).statistic
        for family in FAMILY_ORDER:
            selected = shared_membership[family]
            family_idx = np.flatnonzero(selected)
            # Resample blocks within the target's shared-gene universe.
            target_blocks = blocks[family_idx]
            if len(family_idx):
                target_idx = bootstrap_indices(target_blocks, boot)
                same = observed_sign[family_idx][target_idx] * shared_sign[family][family_idx][target_idx] > 0
                boot_same[family][index] = same.mean() if len(same) else np.nan
            else:
                boot_same[family][index] = np.nan
    genome = {
        "subtype": short,
        "n_donors": len(metadata),
        "n_shared_genes": n_shared,
        "pearson_ols": pearson,
        "pearson_ols_ci95_low": float(np.nanpercentile(boot_correlation, 2.5)),
        "pearson_ols_ci95_high": float(np.nanpercentile(boot_correlation, 97.5)),
        "spearman_ols": spearman,
        "spearman_ols_ci95_low": float(np.nanpercentile(boot_spearman, 2.5)),
        "spearman_ols_ci95_high": float(np.nanpercentile(boot_spearman, 97.5)),
        "pearson_wls_secondary": pearson_wls,
        "spearman_wls_secondary": spearman_wls,
        "block_perm_null_mean": float(null_correlation.mean()),
        "block_perm_observed_minus_null": pearson - float(null_correlation.mean()),
        "block_perm_p_greater": float((1 + (null_correlation >= pearson).sum()) / (permutations + 1)),
        "block_perm_n": permutations,
        "block_bootstrap_n": bootstraps,
        "n_coexpression_blocks": int(len(np.unique(blocks))),
        "discovery_model": "original MSSM age_plus_sex_depth OLS",
    }
    targeted = []
    for family in FAMILY_ORDER:
        row = targeted_summary(
            shared_frame, family, shared_membership[family], shared_sign[family],
            source_size[family],
        )
        observed_same = row["same_direction"]
        observed_fdr = row["same_direction_seaad_fdr"]
        row.update({
            "subtype": short,
            "block_perm_same_null_mean": float(null_same[family].mean()),
            "block_perm_same_observed_minus_null": observed_same - float(null_same[family].mean()),
            "block_perm_same_p_greater": float((1 + (null_same[family] >= observed_same).sum()) / (permutations + 1)),
            "block_perm_fdr_null_mean": float(null_fdr[family].mean()),
            "block_perm_fdr_observed_minus_null": observed_fdr - float(null_fdr[family].mean()),
            "block_perm_fdr_p_greater": float((1 + (null_fdr[family] >= observed_fdr).sum()) / (permutations + 1)),
            "same_direction_fraction_ci95_low": float(np.nanpercentile(boot_same[family], 2.5)),
            "same_direction_fraction_ci95_high": float(np.nanpercentile(boot_same[family], 97.5)),
            "block_perm_n": permutations,
            "block_bootstrap_n": bootstraps,
        })
        targeted.append(row)

    sensitivity = []
    arms = [
        ("primary_wls", counts, metadata, wls, wls.attrs["active_columns"]),
        ("ols_same_design", counts, metadata, ols, wls.attrs["active_columns"]),
    ]
    no_path = [c for c in PRIMARY_COLUMNS if c != "adnc_low"]
    arms.append(("wls_without_adnc", counts, metadata,
                 fit_wls(counts, metadata, genes, no_path), no_path))
    single_counts, single_meta, single_genes = load_counts(data_root, short, arm="singleome")
    single_columns = ("age", "sex_male", "mean_log_umi", "adnc_low")
    arms.append(("singleome_wls", single_counts, single_meta,
                 fit_wls(single_counts, single_meta, single_genes, single_columns), single_columns))
    floor_counts, floor_meta, floor_genes = load_counts(data_root, short, floor=50)
    arms.append(("floor_50_nuclei_wls", floor_counts, floor_meta,
                 fit_wls(floor_counts, floor_meta, floor_genes), PRIMARY_COLUMNS))
    for name, arm_counts, arm_meta, fitted, active in arms:
        arm_frame = fitted if "q_bh" in fitted else fitted
        merged = arm_frame.merge(discovery[["gene", "mssm_slope_original_ols"]], on="gene", how="left")
        slope_column = "slope" if "slope" in fitted else "slope_ols"
        p_column = "pvalue" if "pvalue" in fitted else "pvalue_ols"
        q_column = "q_bh" if "q_bh" in fitted else "q_bh_ols"
        n, r, rho = correlation_rows(merged, "mssm_slope_original_ols", slope_column)
        arm_membership, arm_sign, arm_size, _ = family_membership(merged["gene"], short)
        both_row = targeted_summary(
            merged, "psychad_hbcc_both_fdr",
            arm_membership["psychad_hbcc_both_fdr"],
            arm_sign["psychad_hbcc_both_fdr"],
            arm_size["psychad_hbcc_both_fdr"],
            slope=slope_column, pvalue=p_column, qvalue=q_column,
        )
        sensitivity.append({
            "subtype": short, "arm": name, "status": "tested",
            "n_donors": len(arm_meta), "design": " + ".join(["intercept"] + list(active)),
            "n_shared_genes": n, "pearson_vs_original_mssm_ols": r,
            "spearman_vs_original_mssm_ols": rho,
            **{f"both_fdr_{key}": value for key, value in both_row.items() if key != "family"},
        })
    not_ad_counts, not_ad_meta, not_ad_genes = load_counts(data_root, short, not_ad_only=True)
    try:
        not_ad_fit = fit_wls(not_ad_counts, not_ad_meta, not_ad_genes)
        not_ad_frame = not_ad_fit.merge(
            discovery[["gene", "mssm_slope_original_ols"]], on="gene", how="left"
        )
        n, r, rho = correlation_rows(not_ad_frame, "mssm_slope_original_ols", "slope")
        sensitivity.append({
            "subtype": short, "arm": "not_ad_only_wls", "status": "tested",
            "n_donors": len(not_ad_meta), "n_shared_genes": n,
            "pearson_vs_original_mssm_ols": r, "spearman_vs_original_mssm_ols": rho,
        })
    except (ValueError, AssertionError) as exc:
        sensitivity.append({
            "subtype": short, "arm": "not_ad_only_wls", "status": "untestable",
            "n_donors": len(not_ad_meta), "reason": str(exc),
        })

    lodo_rows = []
    lodo_gene_frames = []
    both_genes = sorted(load_targets(short)[1])
    lodo_genes = lodo_sensitivity(counts, metadata, genes, both_genes, columns=wls.attrs["active_columns"])
    if len(lodo_genes):
        lodo_genes["subtype"] = short
        lodo_gene_frames.append(lodo_genes)
    for drop in range(len(metadata)):
        keep_donor = np.arange(len(metadata)) != drop
        try:
            dropped = fit_wls(counts[keep_donor], metadata.loc[keep_donor].reset_index(drop=True), genes)
            merged = dropped.merge(discovery[["gene", "mssm_slope_original_ols"]], on="gene", how="left")
            n, r, _ = correlation_rows(merged, "mssm_slope_original_ols", "slope")
            drop_membership, drop_sign, drop_size, _ = family_membership(merged["gene"], short)
            row = targeted_summary(
                merged, "psychad_hbcc_both_fdr", drop_membership["psychad_hbcc_both_fdr"],
                drop_sign["psychad_hbcc_both_fdr"], drop_size["psychad_hbcc_both_fdr"],
            )
            lodo_rows.append({"drop_index": drop, "pearson": r, **row})
        except ValueError:
            continue
    lodo_summary = {
        "subtype": short,
        "n_successful_donor_drops": len(lodo_rows),
        "pearson_min": min((r["pearson"] for r in lodo_rows), default=np.nan),
        "pearson_max": max((r["pearson"] for r in lodo_rows), default=np.nan),
        "pearson_sign_flips": int(sum(np.sign(r["pearson"]) != np.sign(pearson_wls) for r in lodo_rows)),
        "both_fdr_same_direction_fraction_min": min((r["same_direction_fraction"] for r in lodo_rows), default=np.nan),
        "both_fdr_same_direction_fraction_max": max((r["same_direction_fraction"] for r in lodo_rows), default=np.nan),
        "both_fdr_genes_evaluated": len(lodo_genes),
        "both_fdr_genes_with_any_slope_sign_flip": int((lodo_genes["n_sign_flip"] > 0).sum()) if len(lodo_genes) else 0,
        "both_fdr_genes_losing_nominal_in_any_drop": int((lodo_genes["n_lose_nominal"] > 0).sum()) if len(lodo_genes) else 0,
    }
    log(f"analysis {short}: {len(metadata)} donors, {len(frame)} genes, {n_shared} shared")
    return frame, genome, targeted, sensitivity, lodo_summary, lodo_gene_frames


def analyze(data_root: Path, discovery_source: Path, permutations: int, bootstraps: int):
    plan = load_plan()
    assert sha256(discovery_source) == plan["genome_wide_discovery_source"]["sha256"]
    frames, genomes, targeted, sensitivities, lodo, lodo_genes = [], [], [], [], [], []
    for short in SUBCLASS:
        outputs = analyze_subtype(data_root, discovery_source, short, permutations, bootstraps)
        frames.append(outputs[0])
        genomes.append(outputs[1])
        targeted.extend(outputs[2])
        sensitivities.extend(outputs[3])
        lodo.append(outputs[4])
        lodo_genes.extend(outputs[5])
    genome = pd.DataFrame(genomes)
    genome["block_perm_q_bh_across_3_subtypes"] = multipletests(
        genome["block_perm_p_greater"], method="fdr_bh"
    )[1]
    targeted = pd.DataFrame(targeted)
    primary = targeted["family"] == "psychad_hbcc_both_fdr"
    targeted.loc[primary, "block_perm_same_q_bh_across_3_subtypes"] = multipletests(
        targeted.loc[primary, "block_perm_same_p_greater"], method="fdr_bh"
    )[1]
    atomic_csv(pd.concat(frames, ignore_index=True), RESULTS / "followup-sea-ad-gene-effects.csv")
    atomic_csv(genome, RESULTS / "followup-sea-ad-genome-wide.csv")
    atomic_csv(targeted, RESULTS / "followup-sea-ad-replication-summary.csv")
    atomic_csv(pd.DataFrame(sensitivities), RESULTS / "followup-sea-ad-sensitivities.csv")
    atomic_csv(pd.DataFrame(lodo), RESULTS / "followup-sea-ad-lodo-summary.csv")
    atomic_csv(pd.concat(lodo_genes, ignore_index=True), RESULTS / "followup-sea-ad-lodo-genes.csv")
    output_manifest = {
        "analysis": "Step 9 SEA-AD interneuron extension",
        "frozen_plan_sha256": plan["frozen_plan_sha256"],
        "script_sha256": sha256(Path(__file__)),
        "source_manifest_sha256": sha256(MANIFEST / "sea-ad-mtg-2024-02-13.json"),
        "report_sha256": (
            sha256(REPO / "reports" / "followup-sea-ad-interneurons.md")
            if (REPO / "reports" / "followup-sea-ad-interneurons.md").exists() else None
        ),
        "permutations": permutations,
        "block_bootstraps": bootstraps,
        "outputs": {},
    }
    for name in (
        "followup-sea-ad-eligibility.csv", "followup-sea-ad-gene-effects.csv",
        "followup-sea-ad-genome-wide.csv", "followup-sea-ad-replication-summary.csv",
        "followup-sea-ad-sensitivities.csv", "followup-sea-ad-lodo-summary.csv",
        "followup-sea-ad-lodo-genes.csv",
    ):
        output_manifest["outputs"][name] = {"sha256": sha256(RESULTS / name)}
    (MANIFEST / "sea_ad_interneuron_analysis_outputs.json").write_text(
        json.dumps(output_manifest, indent=2) + "\n"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=["acquire", "metadata", "analyze", "all"], default="all")
    parser.add_argument("--data-root", type=Path,
                        default=REPO.parent / "agescape-revision-data" / "sea-ad")
    parser.add_argument("--discovery-source", type=Path,
                        default=REPO.parent / "agescape" / "results" /
                        "h3-human-neuronal-gene-replication-all.csv")
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--bootstraps", type=int, default=1000)
    parser.add_argument("--force-metadata", action="store_true")
    args = parser.parse_args()
    assert args.permutations > 0 and args.bootstraps > 0
    if args.stage == "acquire":
        acquire(args.data_root)
        return
    if args.stage in {"metadata", "all"}:
        load_plan()
        build_metadata(args.data_root, force=args.force_metadata)
    if args.stage in {"analyze", "all"}:
        analyze(args.data_root, args.discovery_source, args.permutations, args.bootstraps)


if __name__ == "__main__":
    main()
