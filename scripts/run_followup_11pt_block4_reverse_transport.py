#!/usr/bin/env python3
"""Block 4: derive human aging signatures and transport them to TMS mouse units."""

from __future__ import annotations

import gc
import os
import pickle
import re
import sys
import time
from pathlib import Path

for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_key, "1")

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata.io import read_elem
from scipy import stats
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_followup_between_subject_gao_yazar import (
    build_human_gene_projection,
    read_csr_chunk,
    read_obs_column,
    read_yazar_var_and_slot,
)
from scripts.run_followup_hcata_pool_sensitivity import add_yazar_pool_dummies
from scripts.run_followup_organ_17_heart_gao import (
    H5AD_PATH as GAO_H5AD,
    STUDY_DUMMY_MAP,
    build_donor_audit,
    build_donor_meta_table,
)
from src.followup_muscle import NULL_SEED, fit_ols
from src.followup_organs import (
    run_family_c_p1_signature,
    select_design_columns,
    write_clean_csv,
    write_clean_markdown,
)

RESULTS = REPO / "results"
REPORTS = REPO / "reports"
DATA_ROOT = REPO.parent / f"{REPO.name}-data" / "tmp_organs"
YAZAR_H5AD = DATA_ROOT / "hcata" / "yazar_1e44db10.h5ad"
YAZAR_CACHE = DATA_ROOT / "hcata" / "yazar_cache"
CACHE = DATA_ROOT / "block4_reverse_transport_cache"
MOUSE_DIR = REPO / "data" / "interim"
GENE_OUT = RESULTS / "followup-11pt-block4-reverse-signature-genes.csv"
TRANSPORT_OUT = RESULTS / "followup-11pt-block4-reverse-transport.csv"
REPORT_OUT = REPORTS / "followup-11pt-block4-reverse-transport.md"
AGE_MAP = {"1m": 1.0, "3m": 3.0, "18m": 18.0, "21m": 21.0, "24m": 24.0, "30m": 30.0}
EXPECTED_TISSUES = {
    "Marrow": 39471,
    "Spleen": 34989,
    "Limb_Muscle": 28760,
    "Lung": 24458,
    "Kidney": 20922,
    "Tongue": 20031,
    "Mammary_Gland": 12064,
    "Thymus": 9173,
    "Bladder": 8669,
    "Heart_and_Aorta": 8335,
    "Trachea": 7898,
    "Liver": 7084,
    "Fat": 6383,
    "Pancreas": 4212,
    "Skin": 3962,
    "Large_Intestine": 1775,
}
ITERATIONS = 2000
MIN_CELLS = 20
MIN_DONORS = 6


def available_gb() -> float:
    with open("/proc/meminfo", encoding="utf-8") as handle:
        fields = {line.split(":")[0]: line for line in handle}
    return float(fields["MemAvailable"].split()[1]) / 1024**2


def wait_for_memory(label: str) -> None:
    while available_gb() < 20.0:
        print(f"Waiting for 20 GiB available before {label}; currently {available_gb():.1f} GiB", flush=True)
        time.sleep(60)


def fit_gene_models(
    expr_by_group: dict[str, pd.DataFrame],
    meta_by_group: dict[str, pd.DataFrame],
    covariates: tuple[str, ...],
) -> pd.DataFrame:
    frames = []
    for group in sorted(expr_by_group):
        expr = expr_by_group[group]
        meta = meta_by_group[group]
        donors = [donor for donor in meta.index if donor in expr.columns]
        if len(donors) < MIN_DONORS:
            continue
        design, names, omitted = select_design_columns(
            meta.loc[donors],
            age_col="age_per_decade",
            candidate_covariates=covariates,
            min_residual_df=2,
        )
        outcomes = expr[donors].to_numpy(dtype=np.float64).T
        fitted = fit_ols(outcomes, design, effect_column=1)
        frame = pd.DataFrame({
            "group": group,
            "gene": expr.index.astype(str),
            "n_donors": len(donors),
            "effect": fitted["effect"],
            "se": fitted["se"],
            "pvalue": fitted["pvalue"],
            "residual_df": fitted["residual_df"].astype(int),
            "covariates": ";".join(names[2:]) if len(names) > 2 else "none",
            "omitted_covariates": ";".join(omitted) if omitted else "none",
        })
        valid = frame["pvalue"].notna()
        frame["q_bh"] = np.nan
        if valid.any():
            frame.loc[valid, "q_bh"] = multipletests(frame.loc[valid, "pvalue"], method="fdr_bh")[1]
        frames.append(frame)
        print(f"  {group}: {len(donors)} donors, {valid.sum()} genes tested", flush=True)
    if not frames:
        raise RuntimeError("No eligible human cell-type gene models")
    return pd.concat(frames, ignore_index=True)


def stream_mean_log2cpm(
    h5_path: Path,
    slot: str,
    obs_donors: np.ndarray,
    obs_groups: np.ndarray,
    keep: np.ndarray,
    gene_names: np.ndarray,
    donor_meta: pd.DataFrame,
) -> tuple[dict[str, pd.DataFrame], dict[str, pd.DataFrame]]:
    donors = np.asarray(obs_donors, dtype=str)
    groups = np.asarray(obs_groups, dtype=str)
    keep = np.asarray(keep, dtype=bool) & np.isin(donors, donor_meta.index.astype(str))
    pair_counts = (
        pd.DataFrame({"donor": donors[keep], "group": groups[keep]})
        .groupby(["group", "donor"], sort=True).size()
    )
    eligible = pair_counts[pair_counts >= MIN_CELLS]
    pairs = list(eligible.index)
    pair_lookup = {pair: index for index, pair in enumerate(pairs)}
    cell_pairs = np.full(len(keep), -1, dtype=np.int64)
    for position in np.flatnonzero(keep):
        pair = (groups[position], donors[position])
        if pair in pair_lookup:
            cell_pairs[position] = pair_lookup[pair]

    projection = build_human_gene_projection(gene_names)
    retained_genes = projection["audit"]["retained_symbols"]
    pair_sums = np.zeros((len(pairs), len(retained_genes)), dtype=np.float64)
    library = np.zeros(len(pairs), dtype=np.float64)
    with h5py.File(h5_path, "r") as handle:
        matrix_group = handle[slot]
        indptr = matrix_group["indptr"][:]
        for start in range(0, len(keep), 4000):
            stop = min(start + 4000, len(keep))
            pair_ids = cell_pairs[start:stop]
            local_keep = pair_ids >= 0
            if not local_keep.any():
                continue
            block = read_csr_chunk(matrix_group, start, stop, len(gene_names), indptr)[local_keep]
            projected = (block @ projection["proj"]).tocsr()
            local_library = np.asarray(projected.sum(axis=1)).ravel()
            scale = 1e6 / np.where(local_library > 0, local_library, 1.0)
            projected.data *= np.repeat(scale, np.diff(projected.indptr))
            np.log2(projected.data + 1.0, out=projected.data)
            selected_pairs = pair_ids[local_keep]
            local_pairs, inverse = np.unique(selected_pairs, return_inverse=True)
            indicator = sp.csr_matrix(
                (np.ones(len(selected_pairs)), (np.arange(len(selected_pairs)), inverse)),
                shape=(len(selected_pairs), len(local_pairs)),
            )
            pair_sums[local_pairs] += (indicator.T @ projected).toarray()
            np.add.at(library, selected_pairs, local_library)
            if stop % 40000 == 0 or stop == len(keep):
                print(f"  streamed {stop} of {len(keep)} cells", flush=True)

    expr = {}
    meta = {}
    for group in sorted({pair[0] for pair in pairs}):
        group_donors = sorted(pair[1] for pair in pairs if pair[0] == group)
        if len(group_donors) < MIN_DONORS:
            continue
        positions = [pair_lookup[(group, donor)] for donor in group_donors]
        cell_counts = np.array([eligible.loc[(group, donor)] for donor in group_donors], dtype=float)
        expr[group] = pd.DataFrame(
            (pair_sums[positions] / cell_counts[:, None]).T,
            index=retained_genes,
            columns=group_donors,
        )
        group_meta = donor_meta.loc[group_donors].copy()
        group_meta["n_cells"] = cell_counts.astype(int)
        group_meta["lib_size"] = library[positions]
        group_meta["log_lib"] = np.log1p(group_meta["lib_size"])
        meta[group] = group_meta
    return expr, meta


def load_or_save_cache(cache_path: Path, builder) -> tuple[dict[str, pd.DataFrame], dict[str, pd.DataFrame]]:
    CACHE.mkdir(parents=True, exist_ok=True)
    if cache_path.exists():
        with open(cache_path, "rb") as handle:
            return pickle.load(handle)
    result = builder()
    temporary = cache_path.with_suffix(".tmp.pkl")
    with open(temporary, "wb") as handle:
        pickle.dump(result, handle, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(temporary, cache_path)
    return result


def derive_yazar() -> tuple[dict[str, pd.DataFrame], dict[str, pd.DataFrame], dict[str, list[str]]]:
    donor_meta = pd.read_csv(YAZAR_CACHE / "donor_meta.csv")
    donor_meta["donor_id"] = donor_meta["donor_id"].astype(str)
    donor_meta = donor_meta.set_index("donor_id")
    if len(donor_meta) != 981 or donor_meta["pool_number"].nunique() != 75:
        raise ValueError("Yazar cache does not contain 981 donors and 75 pools")
    donor_meta, pool_columns = add_yazar_pool_dummies(donor_meta)

    def build():
        wait_for_memory("Yazar mean-logCPM streaming")
        with h5py.File(YAZAR_H5AD, "r") as handle:
            slot, genes = read_yazar_var_and_slot(handle)
            donors = read_obs_column(handle["obs"], "donor_id")
            groups = read_obs_column(handle["obs"], "cell_type")
        keep = groups != "peripheral blood mononuclear cell"
        print(f"Streaming {keep.sum()} Yazar cells from {YAZAR_H5AD}", flush=True)
        return stream_mean_log2cpm(YAZAR_H5AD, slot, donors, groups, keep, genes, donor_meta)

    expr, meta = load_or_save_cache(CACHE / "yazar_mean_log2cpm.pkl", build)
    return expr, meta, {"batch_adjusted": ["sex_male", *pool_columns], "unadjusted": ["sex_male"]}


def derive_gao() -> tuple[dict[str, pd.DataFrame], dict[str, pd.DataFrame], dict[str, list[str]]]:
    def build():
        wait_for_memory("Gao mean-logCPM streaming")
        print("Reading Gao obs and var", flush=True)
        with h5py.File(GAO_H5AD, "r") as handle:
            obs = read_elem(handle["obs"])
            var = read_elem(handle["var"])
        genes = var.index.astype(str).to_numpy()
        audit = build_donor_audit(obs)
        donor_meta = build_donor_meta_table(audit)
        if len(donor_meta) != 117 or donor_meta["study"].nunique() != 7:
            raise ValueError("Gao cohort does not contain 117 donors and 7 studies")
        included = set(donor_meta.index)
        obs_donors = obs["donor_id"].astype(str).to_numpy()
        obs_types = obs["final_cell_type"].astype(str).to_numpy()
        keep = np.array([donor in included for donor in obs_donors], dtype=bool)
        print(f"Streaming {keep.sum()} Gao nuclei from {GAO_H5AD}", flush=True)
        return stream_mean_log2cpm(
            GAO_H5AD, "layers/counts", obs_donors, obs_types, keep, genes, donor_meta
        )

    expr, meta = load_or_save_cache(CACHE / "gao_mean_log2cpm.pkl", build)
    return expr, meta, {
        "batch_adjusted": ["sex_male", *STUDY_DUMMY_MAP.values()],
        "unadjusted": ["sex_male"],
    }


def human_biotypes() -> dict[str, str]:
    adata = ad.read_h5ad(YAZAR_H5AD, backed="r")
    var = adata.var[["feature_name", "feature_type"]].copy()
    adata.file.close()
    var["feature_name"] = var["feature_name"].astype(str)
    var["feature_type"] = var["feature_type"].astype(str)
    return var.drop_duplicates("feature_name").set_index("feature_name")["feature_type"].to_dict()


def clean_human_gene(gene: str, biotype: str) -> bool:
    symbol = str(gene).upper()
    if re.match(r"^GM\d+$", symbol):
        return False
    if re.search(r"-PS\d*$", symbol):
        return False
    if re.match(r"^[0-9A-Z]+RIK$", symbol):
        return False
    if symbol.startswith(("RPL", "RPS", "MT-", "MRP")) or symbol == "RN45S":
        return False
    return biotype == "protein_coding"


def make_signature(
    tests: pd.DataFrame,
    cohort: str,
    model: str,
    orthologs: pd.DataFrame,
    biotypes: dict[str, str],
) -> tuple[pd.DataFrame, dict[str, list[str]], dict[str, float | int]]:
    discovered = tests[(tests["q_bh"] < 0.05) & (tests["effect"].abs() >= 0.25)].copy()
    if discovered.empty:
        raise RuntimeError(f"{cohort} {model}: discovery set is empty")
    aggregated = discovered.groupby("gene", sort=True).agg(
        mean_age_slope_log2cpm_per_decade=("effect", "mean"),
        min_pvalue=("pvalue", "min"),
        min_q_bh=("q_bh", "min"),
        n_discovery_cell_types=("group", "nunique"),
        discovery_cell_types=("group", lambda values: ";".join(sorted(set(values)))),
    ).reset_index().rename(columns={"gene": "human_gene_symbol"})
    aggregated["direction"] = np.where(aggregated["mean_age_slope_log2cpm_per_decade"] > 0, "up", "down")
    aggregated["human_biotype"] = aggregated["human_gene_symbol"].map(biotypes).fillna("unknown")
    aggregated["passes_clean_filter"] = [
        clean_human_gene(gene, biotype) for gene, biotype in zip(
            aggregated["human_gene_symbol"], aggregated["human_biotype"], strict=True
        )
    ]
    mapping = orthologs.set_index("human_gene_symbol")["mouse_gene_symbol"].to_dict()
    aggregated["mouse_gene_symbol"] = aggregated["human_gene_symbol"].map(mapping)
    aggregated["mapping_status"] = np.where(
        aggregated["mouse_gene_symbol"].notna(), "mapped_one_to_one", "lost_no_one_to_one"
    )
    aggregated["in_transport_signature"] = (
        aggregated["passes_clean_filter"] & aggregated["mouse_gene_symbol"].notna()
    )
    aggregated.insert(0, "derivation_model", model)
    aggregated.insert(0, "derivation_cohort", cohort)
    aggregated["discovery_stage"] = "discovery"
    clean = aggregated[aggregated["passes_clean_filter"]]
    mapped = aggregated[aggregated["in_transport_signature"]]
    program = {
        "up": mapped.loc[mapped["direction"] == "up", "mouse_gene_symbol"].tolist(),
        "down": mapped.loc[mapped["direction"] == "down", "mouse_gene_symbol"].tolist(),
    }
    counts = {
        "discovery": int(len(aggregated)),
        "clean": int(len(clean)),
        "mapped": int(len(mapped)),
        "lost": int(len(clean) - len(mapped)),
        "coverage": float(len(mapped) / len(clean)) if len(clean) else 0.0,
        "up": int(len(program["up"])),
        "down": int(len(program["down"])),
    }
    print(f"{cohort} {model} signature counts: {counts}", flush=True)
    return aggregated, program, counts


def mouse_targets() -> dict[str, tuple[pd.DataFrame, pd.DataFrame, str, str]]:
    paths = sorted(MOUSE_DIR.glob("m1-qc-cleaned-*.h5ad"))
    if len(paths) != 16 or {path.stem.removeprefix("m1-qc-cleaned-") for path in paths} != set(EXPECTED_TISSUES):
        raise ValueError(f"Expected 16 restored TMS files, found {len(paths)}")
    targets = {}
    for path in paths:
        tissue = path.stem.removeprefix("m1-qc-cleaned-")
        wait_for_memory(f"loading {tissue}")
        adata = ad.read_h5ad(path)
        if adata.shape != (EXPECTED_TISSUES[tissue], 20138):
            raise ValueError(f"{path}: unexpected shape {adata.shape}")
        keep = (
            (~adata.obs["predicted_doublet"].astype(bool).to_numpy())
            & (adata.obs["cell_ontology_class"].astype(str) != "unknown").to_numpy()
            & adata.obs["age"].astype(str).isin(AGE_MAP).to_numpy()
        )
        adata = adata[keep].copy()
        obs = adata.obs
        for cell_type in sorted(obs["cell_ontology_class"].astype(str).unique()):
            cell_mask = obs["cell_ontology_class"].astype(str).eq(cell_type).to_numpy()
            cell_obs = obs.loc[cell_mask]
            counts = cell_obs.groupby("mouse.id", observed=True).size()
            donors = sorted(counts[counts >= MIN_CELLS].index.astype(str))
            if len(donors) < MIN_DONORS:
                continue
            positions = np.flatnonzero(cell_mask & obs["mouse.id"].astype(str).isin(donors).to_numpy())
            donor_values = obs["mouse.id"].astype(str).to_numpy()[positions]
            donor_index = {donor: index for index, donor in enumerate(donors)}
            codes = np.array([donor_index[donor] for donor in donor_values], dtype=int)
            indicator = sp.csr_matrix(
                (np.ones(len(positions)), (np.arange(len(positions)), codes)),
                shape=(len(positions), len(donors)),
            )
            projection = build_human_gene_projection(adata.var_names.astype(str).to_numpy())
            projected = (adata.X[positions] @ projection["proj"]).tocsr()
            cell_library = np.asarray(projected.sum(axis=1)).ravel()
            scale = 1e6 / np.where(cell_library > 0, cell_library, 1.0)
            projected.data *= np.repeat(scale, np.diff(projected.indptr))
            np.log2(projected.data + 1.0, out=projected.data)
            cell_counts = counts.loc[donors].to_numpy(dtype=float)
            mean_expression = (indicator.T @ projected).toarray().T / cell_counts[None, :]
            log2cpm = pd.DataFrame(
                mean_expression,
                index=projection["audit"]["retained_symbols"],
                columns=donors,
            )
            donor_obs = cell_obs[cell_obs["mouse.id"].astype(str).isin(donors)].copy()
            donor_obs["mouse.id"] = donor_obs["mouse.id"].astype(str)
            donor_meta = donor_obs.groupby("mouse.id", observed=True).agg(
                age_months=("age", lambda values: AGE_MAP[str(values.iloc[0])]),
                sex_male=("sex", lambda values: float(str(values.iloc[0]).lower() == "male")),
            ).loc[donors]
            donor_meta["n_cells"] = cell_counts.astype(int)
            donor_meta["lib_size"] = np.bincount(codes, weights=cell_library, minlength=len(donors))
            donor_meta["log_lib"] = np.log1p(donor_meta["lib_size"])
            key = f"{tissue}::{cell_type}"
            targets[key] = (log2cpm, donor_meta, tissue, cell_type)
        print(f"Prepared {tissue}; cumulative eligible units {len(targets)}", flush=True)
        del adata
        gc.collect()
    if not targets:
        raise RuntimeError("No eligible TMS units")
    return targets


def empty_transport_row(
    cohort: str,
    model: str,
    tissue: str,
    cell_type: str,
    meta: pd.DataFrame,
    status: str = "UNTESTABLE_SCORING_ELIGIBILITY",
) -> dict:
    return {
        "derivation_cohort": cohort,
        "derivation_model": model,
        "scoring_direction": "human_to_mouse",
        "tissue": tissue,
        "cell_type": cell_type,
        "status": status,
        "n_donors": len(meta),
        "n_cells_total": int(meta["n_cells"].sum()),
        "n_genes_mapped": 0,
        "n_up_mapped": 0,
        "n_down_mapped": 0,
        "n_random_sets": ITERATIONS,
        "random_seed": NULL_SEED,
    }


def score_targets(
    targets: dict[str, tuple[pd.DataFrame, pd.DataFrame, str, str]],
    programs: dict[tuple[str, str], dict[str, list[str]]],
) -> pd.DataFrame:
    rows = []
    for (cohort, model), program in programs.items():
        print(f"Scoring {cohort} {model}: planned TMS family has {len(targets)} eligible units", flush=True)
        for key in sorted(targets):
            expr, meta, tissue, cell_type = targets[key]
            if not program["up"] or not program["down"]:
                missing = "UP" if not program["up"] else "DOWN"
                rows.append(empty_transport_row(
                    cohort, model, tissue, cell_type, meta, f"UNTESTABLE_MISSING_{missing}_DIRECTION"
                ))
                continue
            result = run_family_c_p1_signature(
                {key: expr},
                {key: meta},
                p1_program=program,
                age_col="age_months",
                candidate_covariates=("sex_male", "log_lib"),
                min_donors=MIN_DONORS,
                min_residual_df=2,
                iterations=ITERATIONS,
                seed=NULL_SEED,
            )
            if result.empty:
                rows.append(empty_transport_row(cohort, model, tissue, cell_type, meta))
                continue
            row = result.iloc[0].to_dict()
            df = int(row["residual_df"])
            mde80 = (stats.t.ppf(0.975, df) + stats.t.ppf(0.80, df)) * float(row["composite_se"])
            rows.append({
                "derivation_cohort": cohort,
                "derivation_model": model,
                "scoring_direction": "human_to_mouse",
                "tissue": tissue,
                "cell_type": cell_type,
                "status": "TESTED",
                "n_donors": int(row["n_donors"]),
                "n_cells_total": int(meta["n_cells"].sum()),
                "n_genes_mapped": int(row["n_p1_total_detected"]),
                "n_up_mapped": int(row["n_p1_up_detected"]),
                "n_down_mapped": int(row["n_p1_down_detected"]),
                "n_expression_pool": int(row["n_expressed_pool"]),
                "composite_slope_per_month": float(row["composite_effect"]),
                "composite_slope_per_12_months": float(row["composite_effect"]) * 12.0,
                "composite_se": float(row["composite_se"]),
                "composite_ci95_low": float(row["composite_ci95_low"]),
                "composite_ci95_high": float(row["composite_ci95_high"]),
                "composite_pvalue": float(row["composite_pvalue"]),
                "empirical_pvalue": float(row["empirical_null_pvalue"]),
                "null_mean_slope": float(row["null_mean_effect"]),
                "null_sd_slope": float(row["null_sd_effect"]),
                "mde80_per_month": mde80,
                "mde80_per_12_months": mde80 * 12.0,
                "residual_df": df,
                "covariates": row["covariates"],
                "omitted_covariates": row["omitted_covariates"],
                "n_random_sets": ITERATIONS,
                "random_seed": NULL_SEED,
            })
            print(f"  {key}: slope {row['composite_effect']:.5g}, empirical p {row['empirical_null_pvalue']:.5g}", flush=True)
    result = pd.DataFrame(rows)
    result["family_composite_size"] = 0
    result["family_empirical_size"] = 0
    result["composite_q"] = np.nan
    result["empirical_q"] = np.nan
    for (_, _), indices in result.groupby(["derivation_cohort", "derivation_model"], sort=True).groups.items():
        valid = result.loc[indices]
        tested = valid.index[valid["status"] == "TESTED"]
        result.loc[tested, "family_composite_size"] = len(tested)
        result.loc[tested, "family_empirical_size"] = len(tested)
        if len(tested):
            result.loc[tested, "composite_q"] = multipletests(
                result.loc[tested, "composite_pvalue"], method="fdr_bh"
            )[1]
            result.loc[tested, "empirical_q"] = multipletests(
                result.loc[tested, "empirical_pvalue"], method="fdr_bh"
            )[1]
    result["passes_both"] = (
        (result["status"] == "TESTED")
        & (result["composite_q"] < 0.05)
        & (result["empirical_q"] < 0.05)
    )
    return result


def arms_disagree_on_every_unit(transport: pd.DataFrame) -> tuple[bool, int]:
    adjusted = transport[(transport["derivation_model"] == "batch_adjusted") & (transport["status"] == "TESTED")]
    wide = adjusted.pivot(index=["tissue", "cell_type"], columns="derivation_cohort", values="composite_slope_per_month").dropna()
    if not {"Yazar", "Gao"}.issubset(wide.columns):
        return False, 0
    wide = wide.dropna()
    if wide.empty:
        return False, 0
    disagreement = np.sign(wide["Yazar"]) * np.sign(wide["Gao"]) < 0
    return bool(disagreement.all()), len(wide)


def arm_verdict(cohort: str, counts: dict, transport: pd.DataFrame, every_unit_disagrees: bool) -> str:
    if counts["coverage"] < 0.50 or counts["up"] == 0 or counts["down"] == 0:
        return "INCONCLUSIVE"
    arm = transport[
        (transport["derivation_cohort"] == cohort)
        & (transport["derivation_model"] == "batch_adjusted")
    ]
    if arm["passes_both"].any():
        return "PROVES"
    if every_unit_disagrees:
        return "INCONCLUSIVE"
    return "DOES NOT PROVE"


def fmt(value: float) -> str:
    return f"{value:.4g}" if pd.notna(value) else "NA"


def build_report(
    genes: pd.DataFrame,
    transport: pd.DataFrame,
    counts: dict[tuple[str, str], dict],
    verdicts: dict[str, str],
    every_unit_disagrees: bool,
    shared_units: int,
) -> str:
    lines = [
        "# Block 4: Reverse signature transport",
        "",
        "## Spec",
        "",
        "The primary derivation cohort is Yazar blood and the sensitivity derivation cohort is Gao heart. For each human cell type and donor with at least 20 cells, full-depth raw counts receive the human technical gene mask and are converted to per-cell log2(CPM + 1), then averaged across cells within each donor and gene. Donor mean log2CPM expression is fit by OLS on age per decade plus sex plus batch. Yazar uses 75 pool fixed effects and Gao uses 7 study indicators. At least 6 donors are required for every slope. Discovery requires BH FDR below 0.05 within each derivation cell type and absolute log2 fold change at least 0.25, mirroring P0. A gene discovered in multiple cell types is included once and its direction is the sign of its mean discovery slope, as in the P0 hierarchy.",
        "",
        "The P1 clean filter is adapted to human symbols: protein-coding only, excluding GM-prefixed, -PS, RIK, ribosomal, mitochondrial, and MRP genes. Clean genes are split by age-slope sign and mapped to mouse with the Ensembl release 106 one-to-one table used by scripts/build_reference_orthologs.py, in the reverse human-to-mouse direction.",
        "",
        "Each TMS target unit requires at least 20 cells per donor and at least 6 donors. Predicted doublets and unknown cell types are excluded. TMS expression is likewise per-cell log2(CPM + 1) averaged within each donor. The composite score is mean up-gene expression minus mean down-gene expression. OLS uses mouse age in months plus sex plus log library depth. The matched random-set comparison reuses run_family_c_p1_signature and random_signature_scores unchanged: 2000 random sets, matched separately by direction on gene count and donor-mean expression decile, seed 20261105, and a two-sided empirical p centered on the null mean. MDE80 is (t_0.975 + t_0.80) times the composite standard error using residual degrees of freedom.",
        "",
        "The shared-rule unadjusted sensitivity repeats derivation without pool or study indicators and repeats TMS scoring with the resulting signature. These sensitivity families are separate from the primary batch-adjusted families.",
        "",
        "Before correction, families are defined within each derivation cohort, derivation model, and human-to-mouse scoring direction across tested TMS units, separately for composite and empirical p-values. The family sizes are recorded below and in the CSV. A unit passes only when both q-values are below 0.05.",
        "",
        "The human technical mask excludes uppercase MT-, RPL, RPS, MRPL, MRPS, XIST, TSIX, DDX3Y, EIF2S3Y, UTY, and KDM5D genes before library normalization. No downsampling or rarefaction is used.",
        "",
        "## Input restoration",
        "",
        "The 16 per-tissue TMS files were restored by mechanically splitting the combined m1-qc-cleaned.h5ad per Reza's explicit instruction on Nov 20, 2026. Split verification passed for exactly 16 files, 238186 total cells, and 20138 genes per file. Every cell retained its tissue label, every gene list was identical to the combined file, and the split obs indices formed a complete duplicate-free union of the combined index.",
        "",
        "## Derivation coverage and family sizes",
        "",
        "| Cohort | Model | Discovery | Clean | Mapped | Lost | Coverage | Up | Down | Composite family | Empirical family |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for key in sorted(counts):
        c = counts[key]
        arm = transport[(transport["derivation_cohort"] == key[0]) & (transport["derivation_model"] == key[1])]
        tested = arm[arm["status"] == "TESTED"]
        lines.append(
            f"| {key[0]} | {key[1]} | {c['discovery']} | {c['clean']} | {c['mapped']} | {c['lost']} | "
            f"{100*c['coverage']:.1f}% | {c['up']} | {c['down']} | {len(tested)} | {len(tested)} |"
        )
    lines.extend([
        "",
        "Discovery counts are unique genes after pooling cell-type discovery rows. Clean counts retain the P1 filter. Mapped counts and coverage are relative to the clean set. Full per-gene slopes, p-values, q-values, biotypes, directions, and mapping statuses are in the signature-genes CSV.",
        "",
        "## Asymmetry",
        "",
        "Human batch noise enters the DERIVATION in this block, so a reverse null is weaker evidence than the mouse-to-human null, where derivation came from one controlled experiment. A reverse null is therefore not equally strong evidence of symmetric transport failure.",
        "",
        "## Verdict rules and results",
        "",
        "PROVES requires at least one TMS unit passing both corrected tests in that batch-adjusted arm. DOES NOT PROVE requires no unit passing both tests in either derivation arm, with MDE80 reported per unit. INCONCLUSIVE applies when under 50 percent of the clean set maps, when a mapped signature has no up or no down genes so its composite is not estimable, or when the Yazar and Gao adjusted composite slopes disagree in sign on every commonly tested unit. Units in a missing-direction arm are UNTESTABLE, not nulls. The sign-disagreement condition was evaluated before verdict calculation; it covered "
        f"{shared_units} shared units and was {every_unit_disagrees}.",
        "",
        f"- Yazar primary arm: {verdicts['Yazar']}.",
        f"- Gao sensitivity arm: {verdicts['Gao']}.",
        "",
    ])
    for cohort in ("Yazar", "Gao"):
        arm = transport[
            (transport["derivation_cohort"] == cohort)
            & (transport["derivation_model"] == "batch_adjusted")
        ].copy()
        passing = arm[arm["passes_both"]]
        lines.extend([
            f"### {cohort} batch-adjusted TMS results",
            "",
            f"Passing both tests: {len(passing)} of {int((arm['status'] == 'TESTED').sum())} tested units.",
            "",
            "| Tissue | Cell type | N donors | N genes | Slope per 12 months | Composite q | Empirical q | MDE80 per 12 months | Status |",
            "|---|---|---:|---:|---:|---:|---:|---:|---|",
        ])
        for row in arm.sort_values(["tissue", "cell_type"]).itertuples():
            if row.status == "TESTED":
                lines.append(
                    f"| {row.tissue} | {row.cell_type} | {row.n_donors} | {row.n_genes_mapped} | "
                    f"{fmt(row.composite_slope_per_12_months)} | {fmt(row.composite_q)} | {fmt(row.empirical_q)} | "
                    f"{fmt(row.mde80_per_12_months)} | {row.status} |"
                )
            else:
                lines.append(
                    f"| {row.tissue} | {row.cell_type} | {row.n_donors} | 0 | NA | NA | NA | NA | {row.status} |"
                )
        lines.append("")
    lines.extend([
        "## Unadjusted sensitivity",
        "",
        "| Cohort | Tested units | Passing both | Minimum composite q | Minimum empirical q |",
        "|---|---:|---:|---:|---:|",
    ])
    for cohort in ("Yazar", "Gao"):
        arm = transport[(transport["derivation_cohort"] == cohort) & (transport["derivation_model"] == "unadjusted")]
        tested = arm[arm["status"] == "TESTED"]
        lines.append(
            f"| {cohort} | {len(tested)} | {int(tested['passes_both'].sum())} | "
            f"{fmt(tested['composite_q'].min())} | {fmt(tested['empirical_q'].min())} |"
        )
    lines.extend([
        "",
        "## Files",
        "",
        "- scripts/split_m1_qc_cleaned_per_tissue.py",
        "- scripts/run_followup_11pt_block4_reverse_transport.py",
        "- results/followup-11pt-block4-reverse-signature-genes.csv",
        "- results/followup-11pt-block4-reverse-transport.csv",
        "",
    ])
    return "\n".join(lines)


def main() -> None:
    wait_for_memory("Block 4")
    for path in (YAZAR_H5AD, GAO_H5AD, REPO / "data" / "reference" / "orthologs-mouse-human.csv"):
        if not path.exists():
            raise FileNotFoundError(path)
    orthologs = pd.read_csv(REPO / "data" / "reference" / "orthologs-mouse-human.csv")
    orthologs = orthologs[orthologs["orthology_type"] == "ortholog_one2one"]
    orthologs = orthologs.drop_duplicates("human_gene_symbol", keep="last")
    biotypes = human_biotypes()

    all_genes = []
    programs = {}
    counts = {}
    for cohort, loader in (("Yazar", derive_yazar), ("Gao", derive_gao)):
        print(f"Loading {cohort} derivation cohort", flush=True)
        expr, meta, model_covariates = loader()
        for model, covariates in model_covariates.items():
            print(f"Deriving {cohort} {model} signature", flush=True)
            tests = fit_gene_models(expr, meta, tuple(covariates))
            gene_frame, program, count = make_signature(tests, cohort, model, orthologs, biotypes)
            all_genes.append(gene_frame)
            programs[(cohort, model)] = program
            counts[(cohort, model)] = count
        del expr, meta
        gc.collect()

    genes = pd.concat(all_genes, ignore_index=True)
    targets = mouse_targets()
    transport = score_targets(targets, programs)
    every_unit_disagrees, shared_units = arms_disagree_on_every_unit(transport)
    verdicts = {
        cohort: arm_verdict(cohort, counts[(cohort, "batch_adjusted")], transport, every_unit_disagrees)
        for cohort in ("Yazar", "Gao")
    }
    report = build_report(genes, transport, counts, verdicts, every_unit_disagrees, shared_units)
    write_clean_csv(genes, GENE_OUT)
    write_clean_csv(transport, TRANSPORT_OUT)
    write_clean_markdown(report, REPORT_OUT)
    print(f"Wrote {GENE_OUT} ({len(genes)} rows)", flush=True)
    print(f"Wrote {TRANSPORT_OUT} ({len(transport)} rows)", flush=True)
    print(f"Wrote {REPORT_OUT}; verdicts {verdicts}", flush=True)


if __name__ == "__main__":
    main()
