#!/usr/bin/env python3
"""Mature concordance hierarchy: Level 0 per-unit age-slope significance and
the missing Level 2 whole-species tests (human cross-cohort, fly head vs body).

Level 0 refits the observed per-gene OLS age models of the existing
concordance pipeline (no permutations), applies Benjamini-Hochberg within each
unit over the gene set the pipeline used for that unit's slope vector, and
verifies the reconstruction by recomputing the committed Level 1 observed
statistics from the same slopes. Level 2 computes the mean Pearson correlation
over unit pairs whose cohort (human) or organ (fly) differs, with a donor-level
age-permutation null. Existing within-tissue logic and outputs are not touched.

Modes:
  level0 --species <name>   per-unit OLS + BH, fidelity check  -> scratch CSVs
  level2 --species Human|Fruit fly [--permutations 1000]        -> scratch CSVs + pairs CSV
  assemble                  final results/followup-mature-hierarchy-level{0,2}.csv
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
from scripts import run_rev2_concordance as rev2  # noqa: E402
from scripts import run_rev2_joint_concordance as joint  # noqa: E402
from scripts import run_reviewer_phase1_gene_concordance as phase1  # noqa: E402

RESULTS = REPO / "results"
SCRATCH = REPO / "scratch" / "hierarchy"
NMR_RAW_CANDIDATES = [
    REPO / "data" / "raw" / "GSE152886_nmr_skin_raw.csv.gz",
    REPO.parent / "agescape" / "data" / "raw" / "GSE152886_nmr_skin_raw.csv.gz",
]
LEVEL0_COLUMNS = [
    "species", "tissue_cohort", "cell_type", "n_donors_or_pools",
    "n_genes_tested", "n_significant_q05", "frac_significant",
]
LEVEL2_COLUMNS = [
    "species", "scope", "observed_mean_pearson", "joint_null_mean",
    "joint_null_sd", "joint_p_greater", "n_pairs", "n_units",
    "n_donors_or_pools", "n_permutations",
]
SLUGS = {
    "Mouse": "mouse", "Human": "human", "Fruit fly": "fly",
    "C. elegans": "worm", "Zebrafish": "zebrafish", "Naked mole-rat": "nmr",
}
EXPECTED_UNITS = {
    "Mouse": 59, "Human": 127, "Fruit fly": 99, "C. elegans": 10, "Zebrafish": 14,
    "Naked mole-rat": 2,
}
EXPECTED_LEVEL0_TOTAL = sum(EXPECTED_UNITS.values())  # 311: 59 + 127 + 99 + 10 + 14 + 2
# The 7 unidentifiable fly units stay in the table as NA rows (n_significant_q05
# and frac_significant NA, n_genes_tested from the stratum mask), so the final
# Level 0 CSV has exactly 311 rows, one per unit.
# Gene set behind each species' committed Level 1 slope vectors.
PIPELINE_GENE_SET = {
    "Mouse": "species_primary",        # rev2.analyse_mouse: one mask over all 59 units
    "Human": "stratum_primary",        # joint.analyze_species: per-stratum mask
    "Fruit fly": "stratum_primary",
    "C. elegans": "stratum_primary",
    "Zebrafish": "stratum_primary",
    "Naked mole-rat": "unfiltered",     # phase 1 NMR run used every gene after the technical mask
}
FIDELITY_TOLERANCE = 1e-6
# Mouse mean tolerance. The mouse strata are rebuilt from the raw TMS file
# (no cache exists), and the rebuild no longer reproduces the original
# scrublet doublet calls bit-for-bit: Fat 20 vs 15 committed, Kidney 77 vs
# 75, the other 14 tissues exact (results/m1-qc-summary.csv). That fixed
# offset shifts the within-tissue mean by 0.000667 and the cross-tissue mean
# by 0.000015, while units (59), primary genes (10633) and within pairs
# (260) match exactly. 1e-3 sits between the rebuild noise (<1e-3) and a
# real error such as the wrong gene set (>7e-3: primary 0.2105 vs lenient
# 0.2039 vs unfiltered 0.2028 in results/followup-mature-concordance-mouse.csv),
# so the mean checks use it and every structural check keeps 1e-6.
MOUSE_MEAN_FIDELITY_TOLERANCE = 1e-3


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def bh_qvalues(pvalues):
    pvalues = np.asarray(pvalues, dtype=float)
    if pvalues.size == 0:
        return np.array([], dtype=float)
    if not np.isfinite(pvalues).all():
        raise ValueError("BH correction requires finite p-values")
    return multipletests(pvalues, method="fdr_bh")[1]


def bh_count(pvalues, alpha=0.05):
    qvalues = bh_qvalues(pvalues)
    return int((qvalues < alpha).sum()), qvalues


def detection_fractions(unit):
    values = unit["counts"] if "counts" in unit else unit["expression"]
    return (np.asarray(values) > 0).mean(axis=0)


def primary_mask(units):
    """Rev2 primary rule: detected in >=50% of donor pseudobulks in >=50% of units."""
    fractions = np.vstack([detection_fractions(unit) for unit in units])
    return (fractions >= 0.5).sum(axis=0) >= np.ceil(0.5 * len(units))


def clean_unit_design(metadata, age=None):
    """Pipeline design (rev2.unit_design) with float32-jitter columns dropped.

    rev2.unit_design keeps mean_log_umi whenever it has more than one unique
    value. The fly caches store recovered log1p-CP10K values, so mean_log_umi
    there is log(10000) plus float32 rounding noise (range under 1e-4): a
    near-constant column, nearly collinear with the intercept, which
    destabilises pinv(X.T @ X). Drop it when its range is below 1e-4. Column
    order matches rev2.unit_design exactly (const, age, sex?, umi?), so
    age_column=1 still holds.
    """
    import statsmodels.api as sm

    data = metadata.copy()
    data["age_model"] = data["age"] if age is None else age
    columns = ["age_model"]
    if data["sex_male"].nunique() > 1:
        columns.append("sex_male")
    if "mean_log_umi" in data and data["mean_log_umi"].nunique() > 1:
        if np.ptp(np.asarray(data["mean_log_umi"], dtype=float)) >= 1e-4:
            columns.append("mean_log_umi")
    return sm.add_constant(data[columns].astype(float), has_constant="add").to_numpy()


def safe_vectorized_ols(expression, design, age_column=1):
    """OLS slopes, SEs and p-values via Frisch-Waugh-Lovell (QR), no pinv(X.T @ X).

    rev2.vectorized_ols computes the age variance as pinv(design.T @ design),
    which squares the condition number: for ill-conditioned designs (a
    near-constant jitter column, near-collinear age/sex) pinv truncates
    singular values inconsistently and the SEs come out wrong. FWL instead
    residualises age and every gene on the other columns with one QR of the
    reduced design, then age_variance = 1 / ||r_age||^2 exactly.

    Returns (slopes, se, pvalues, identifiable). When the age slope is not
    identified (rank-deficient design, zero residual degrees of freedom, a
    single age value, or age residualised to ~0) returns NaN slopes/SEs and
    pvalues of 1 with identifiable=False; the caller must not turn those
    into hits. Genes whose expression is constant or whose fit is exact to
    floating-point dust (rss <= 1e-15 * tss) get se=0, t=0, pvalue=1: an
    interpolated point is not evidence.
    """
    from scipy import stats

    expression = np.asarray(expression, dtype=float)
    design = np.asarray(design, dtype=float)
    n, p = design.shape
    rank = np.linalg.matrix_rank(design)
    age = design[:, age_column]
    if rank < p or n - rank <= 0 or np.ptp(age) == 0:
        nan = np.full(expression.shape[1], np.nan)
        return nan, nan, np.ones(expression.shape[1]), False
    others = np.delete(design, age_column, axis=1)
    q, _ = np.linalg.qr(others)
    r_age = age - q @ (q.T @ age)
    age_resid_ss = float(r_age @ r_age)
    if age_resid_ss <= 1e-10 * float(np.sum((age - age.mean()) ** 2)):
        nan = np.full(expression.shape[1], np.nan)
        return nan, nan, np.ones(expression.shape[1]), False
    r_y = expression - q @ (q.T @ expression)
    slopes = (r_age @ r_y) / age_resid_ss
    df = n - rank
    fitted = np.outer(r_age, slopes)
    rss = np.sum((r_y - fitted) ** 2, axis=0)
    tss = np.sum((expression - expression.mean(axis=0)) ** 2, axis=0)
    residual_variance = rss / df
    se = np.sqrt(residual_variance / age_resid_ss)
    se = np.where((tss == 0) | (rss <= 1e-15 * np.maximum(tss, 1e-300)), 0.0, se)
    statistics = np.divide(slopes, se, out=np.zeros_like(slopes), where=se > 0)
    pvalues = 2 * stats.t.sf(np.abs(statistics), df)
    pvalues = np.where(se == 0, 1.0, pvalues)
    return slopes, se, pvalues, True


def clean_fit_age_slope_only(expression, metadata, ages):
    """Pipeline age slope only, on the clean design (jitter UMI dropped).

    Same slope as rev2.fit_age_slope_only on a well-conditioned design, but
    uses pinv of the design itself (never pinv(X.T @ X)) and the clean
    design, so the fly jitter column cannot leak into Level 2 slopes.
    """
    design = clean_unit_design(metadata, age=np.asarray(ages, dtype=float))
    return np.linalg.pinv(design)[1] @ np.asarray(expression, dtype=float)


def unit_ols(unit):
    """Observed fit of the concordance pipeline's model, safe OLS version."""
    if "platform_nextseq" in unit["metadata"]:
        raise ValueError("platform_nextseq present: rev2.unit_design and phase1.design differ")
    expression = np.asarray(unit["expression"], dtype=np.float64)
    slopes, se, pvalues, identifiable = safe_vectorized_ols(
        expression, clean_unit_design(unit["metadata"]), age_column=1)
    return slopes, se, pvalues


def age_identifiable(unit):
    """True when the age slope is identified: the clean design must be full
    rank, have residual degrees of freedom, a varying age column, and age
    residualised on the other columns must keep non-trivial variance.
    Units failing this (single age value, or age collinear with sex / mean
    log UMI) carry no Level 0 p-value and no Level 2 slope. Fidelity checks
    still use every pipeline unit."""
    design = clean_unit_design(unit["metadata"])
    n, p = design.shape
    if np.linalg.matrix_rank(design) < p or n - np.linalg.matrix_rank(design) <= 0:
        return False
    age = design[:, 1]
    if np.ptp(age) == 0:
        return False
    others = np.delete(design, 1, axis=1)
    q, _ = np.linalg.qr(others)
    r_age = age - q @ (q.T @ age)
    return bool(float(r_age @ r_age) > 1e-10 * float(np.sum((age - age.mean()) ** 2)))


def mean_pearson(slopes, mask):
    values = rev2.pairwise_pearson(np.asarray(slopes)[:, mask])
    if len(values) != len(slopes) * (len(slopes) - 1) // 2:
        raise ValueError("Zero-variance slope vector: pair order would be incomplete")
    return values


# ---------------------------------------------------------------------------
# Data sources
# ---------------------------------------------------------------------------

def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 22), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolve_nmr_raw():
    manifest = json.loads((REPO / "data" / "manifest" / "GSE152886_nmr_skin_raw.json").read_text())
    for path in NMR_RAW_CANDIDATES:
        if path.exists():
            digest = sha256(path)
            if digest != manifest["sha256"]:
                raise ValueError(f"sha256 mismatch for {path}: {digest}")
            print(f"NMR raw counts: {path} (sha256 matches manifest)", flush=True)
            return path
    raise FileNotFoundError("No naked mole-rat raw skin count file found: "
                            + ", ".join(map(str, NMR_RAW_CANDIDATES)))


@contextmanager
def cross_species_repo_root(root):
    """nmr_annotated_data() reads <REPO>/data/raw and <REPO>/data/interim; point
    it at the checkout that holds the sha-verified raw file and symbol table."""
    from scripts import run_cross_species_expression_vector as cross

    original = cross.REPO
    cross.REPO = root
    try:
        yield cross
    finally:
        cross.REPO = original


def nmr_strata():
    """Exact replica of phase1.nmr_strata with the verified raw-file location."""
    from scripts.run_cross_species_expression_vector import technical_gene_mask

    path = resolve_nmr_raw()
    root = path.parents[2]
    with cross_species_repo_root(root) as cross:
        annotated = cross.nmr_annotated_data()
    raw = pd.read_csv(path, sep="\t", index_col=0)
    counts = raw[annotated.obs_names].T.to_numpy()
    obs = annotated.obs.copy()
    obs["sex_male"] = 0.0
    keep = ~technical_gene_mask(raw.index.astype(str))
    units = [phase1.aggregate_unit(counts[:, keep], obs, cell_type, "Naked mole-rat", "Skin",
                                   "GSE152886", "cell_type", "sample")
             for cell_type in sorted(obs["cell_type"].unique())]
    yield "Skin", "GSE152886", raw.index.to_numpy()[keep], [unit for unit in units if unit]


def filter_unit_min_age(unit, min_age):
    """Drop donor rows below min_age from a cached pseudobulk unit.

    The worm cache spans days 1, 6, 12, 14; the committed mature worm run is
    6-plus days, so day-1 rows are dropped from metadata, expression, and
    counts (when present) before any mask or fit is computed. Pseudobulks are
    per donor, so dropping donor rows is exactly the pre-aggregation filter.
    """
    keep = np.asarray(unit["metadata"]["age"]) >= min_age
    if not keep.any():
        return None
    unit = dict(unit)
    unit["metadata"] = unit["metadata"].loc[keep].reset_index(drop=True)
    unit["expression"] = np.asarray(unit["expression"])[keep]
    if "counts" in unit:
        unit["counts"] = np.asarray(unit["counts"])[keep]
    return unit


def level0_source(species):
    if species == "Mouse":
        for tissue, genes, units, qc in rev2.mouse_strata(use_scrublet=True, mature_only=True):
            print(f"Mouse QC: {qc}", flush=True)
            yield tissue, "TMS", genes, units
    elif species == "Naked mole-rat":
        yield from nmr_strata()
    else:
        # joint.cached_strata takes no age argument: the human cache is
        # already 30-plus and the fly cache already 30-plus. Only worm
        # needs an explicit mature filter (6-plus days).
        for tissue, cohort, genes, units in joint.cached_strata(species):
            if species == "C. elegans":
                units = [u for u in (filter_unit_min_age(unit, 6) for unit in units) if u]
            yield tissue, cohort, genes, units


# ---------------------------------------------------------------------------
# Level 0
# ---------------------------------------------------------------------------

def level0_rows_for_unit(species, tissue, cohort, unit, pvalues, se, masks):
    rows = []
    for gene_set, mask in masks.items():
        tested = pvalues[mask]
        n_significant, qvalues = bh_count(tested)
        rows.append({
            "species": species,
            "tissue_cohort": f"{tissue} | {cohort}",
            "cell_type": unit["cell_type"],
            "n_donors_or_pools": len(unit["metadata"]),
            "n_genes_tested": int(mask.sum()),
            "n_significant_q05": n_significant,
            "frac_significant": n_significant / int(mask.sum()),
            "gene_set": gene_set,
            "n_zero_se_genes": int((se[mask] == 0).sum()),
            "min_q": float(qvalues.min()),
        })
    return rows


def pipeline_ols(unit):
    """The pipeline's own fit (rev2.vectorized_ols on rev2.unit_design).

    Used ONLY for the Level 1 fidelity check: the committed Level 1 numbers
    were produced with this exact fit, so reproducing them requires it. The
    safe OLS is used for every Level 0 inference; the two agree to ~1e-12 on
    well-conditioned designs and differ only where pinv(X.T @ X) is unstable.
    """
    expression = np.asarray(unit["expression"], dtype=np.float64)
    return rev2.vectorized_ols(expression, rev2.unit_design(unit["metadata"]), age_column=1)


def collect_level0(species, strata, min_donors):
    """Fit every unit once (safe OLS for inference, pipeline OLS for fidelity);
    return per-unit fits plus all candidate masks."""
    prepared = []
    for tissue, cohort, genes, units in strata:
        units = [unit for unit in units if len(unit["metadata"]) >= min_donors]
        if not units:
            continue
        fits = [unit_ols(unit) for unit in units]
        pipeline_fits = [pipeline_ols(unit) for unit in units]
        prepared.append({"tissue": tissue, "cohort": cohort, "genes": np.asarray(genes),
                         "units": units, "fits": fits, "pipeline_fits": pipeline_fits,
                         "identifiable": [age_identifiable(unit) for unit in units],
                         "stratum_primary": primary_mask(units)})
        print(f"Level 0 {species} | {tissue} | {cohort}: {len(units)} units fitted", flush=True)
    same_genes = all(np.array_equal(item["genes"], prepared[0]["genes"]) for item in prepared)
    species_mask = (primary_mask([unit for item in prepared for unit in item["units"]])
                    if same_genes else None)
    return prepared, species_mask


def level0_table(species, prepared, species_mask):
    rows, excluded = [], []
    for item in prepared:
        masks = {"stratum_primary": item["stratum_primary"],
                 "unfiltered": np.ones(len(item["genes"]), dtype=bool)}
        if species_mask is not None:
            masks["species_primary"] = species_mask
        for unit, (_, se, pvalues), ok in zip(item["units"], item["fits"],
                                              item["identifiable"]):
            if not ok:
                # Keep the unit in the table as an NA row so the final CSV has
                # exactly one row per unit (311 total). n_genes_tested is the
                # pipeline gene-set size for its stratum; the counts are NA
                # because no age slope is identified.
                pipeline_set = PIPELINE_GENE_SET[species]
                pipeline_mask = masks[pipeline_set]
                excluded.append({"species": species,
                                 "tissue_cohort": f"{item['tissue']} | {item['cohort']}",
                                 "cell_type": unit["cell_type"],
                                 "n_donors_or_pools": len(unit["metadata"]),
                                 "reason": "age slope not identified (design rank unchanged when age dropped)"})
                rows.append({
                    "species": species,
                    "tissue_cohort": f"{item['tissue']} | {item['cohort']}",
                    "cell_type": unit["cell_type"],
                    "n_donors_or_pools": len(unit["metadata"]),
                    "n_genes_tested": int(pipeline_mask.sum()),
                    "n_significant_q05": pd.NA,
                    "frac_significant": pd.NA,
                    "gene_set": pipeline_set,
                    "n_zero_se_genes": pd.NA,
                    "min_q": pd.NA,
                })
                continue
            rows.extend(level0_rows_for_unit(species, item["tissue"], item["cohort"], unit,
                                             pvalues, se, masks))
    frame = pd.DataFrame(rows)
    frame["pipeline_gene_set"] = frame["gene_set"] == PIPELINE_GENE_SET[species]
    frame["n_significant_q05"] = frame["n_significant_q05"].astype("Int64")
    return frame, pd.DataFrame(excluded)



def fidelity_checks(species, prepared, species_mask):
    """Recompute committed Level 1 observed statistics from the Level 0 slopes."""
    rows = []

    def add(label, observed, expected, source, tolerance=FIDELITY_TOLERANCE):
        rows.append({"species": species, "check": label, "recomputed": observed,
                     "committed": expected, "abs_diff": abs(observed - expected),
                     "source": source, "passed": abs(observed - expected) <= tolerance})

    if species == "Mouse":
        committed = pd.read_csv(RESULTS / "followup-mature-concordance-mouse.csv")
        committed = committed[committed["gene_set"] == "primary_detection_filtered"].set_index("scope")
        slopes = np.vstack([fit[0] for item in prepared for fit in item["pipeline_fits"]])
        tissues = np.array([item["tissue"] for item in prepared for _ in item["units"]])
        values = mean_pearson(slopes, species_mask)
        first, second = np.triu_indices(len(slopes), 1)
        same = tissues[first] == tissues[second]
        source = "results/followup-mature-concordance-mouse.csv"
        add("within_tissue_mean", float(values[same].mean()),
            float(committed.loc["within_tissue", "observed_mean_pearson"]), source,
            tolerance=MOUSE_MEAN_FIDELITY_TOLERANCE)
        add("cross_tissue_mean", float(values[~same].mean()),
            float(committed.loc["cross_tissue", "observed_mean_pearson"]), source,
            tolerance=MOUSE_MEAN_FIDELITY_TOLERANCE)
        add("n_units", len(slopes), float(committed.loc["within_tissue", "n_units"]), source)
        add("n_genes_primary", int(species_mask.sum()),
            float(committed.loc["within_tissue", "n_genes"]), source)
        add("n_within_pairs", int(same.sum()), float(committed.loc["within_tissue", "n_pairs"]), source)
    elif species == "Naked mole-rat":
        source = "results/reviewer-phase1-universality-permutations-nmr.csv"
        committed = pd.read_csv(RESULTS / "reviewer-phase1-universality-permutations-nmr.csv")
        observed_row = committed[committed["iteration"] == 0].iloc[0]
        item = prepared[0]
        slopes = np.vstack([fit[0] for fit in item["pipeline_fits"]])
        values = mean_pearson(slopes, np.ones(slopes.shape[1], dtype=bool))
        add("within_skin_mean_unfiltered", float(values.mean()),
            float(observed_row["mean_pearson"]), source)
        add("n_units", len(slopes), float(observed_row["n_units"]), source)
        units = pd.read_csv(RESULTS / "reviewer-phase1-universality-units.csv")
        nmr_genes = units.loc[units["species"] == "Naked mole-rat", "n_genes"].unique()
        add("n_genes_unfiltered", slopes.shape[1], float(nmr_genes[0]),
            "results/reviewer-phase1-universality-units.csv")
    else:
        path = {
            "Human": "followup-mature-concordance-human-strata.csv",
            "Fruit fly": "followup-mature-concordance-fly-strata.csv",
            "C. elegans": "followup-mature-concordance-worm-strata.csv",
            "Zebrafish": "reviewer-rev2-joint-concordance-zebrafish-strata.csv",
        }[species]
        committed = pd.read_csv(RESULTS / path)
        for item in prepared:
            row = committed[(committed["tissue"] == item["tissue"])
                            & (committed["cohort"] == item["cohort"])]
            if len(row) != 1:
                raise ValueError(f"No unique committed stratum row for {item['tissue']} {item['cohort']}")
            row = row.iloc[0]
            slopes = np.vstack([fit[0] for fit in item["pipeline_fits"]])
            values = mean_pearson(slopes, item["stratum_primary"])
            label = f"{item['tissue']} | {item['cohort']}"
            add(f"{label} mean", float(values.mean()), float(row["observed_mean_pearson"]),
                f"results/{path}")
            add(f"{label} n_units", len(slopes), float(row["n_units"]), f"results/{path}")
            add(f"{label} n_genes", int(item["stratum_primary"].sum()), float(row["n_genes"]),
                f"results/{path}")
            donors = {donor for unit in item["units"] for donor in unit["metadata"]["donor"]}
            add(f"{label} n_donors", len(donors), float(row["n_donors_in_stratum"]), f"results/{path}")
    return pd.DataFrame(rows)


def run_level0(species, scratch_dir):
    minimum = 5 if species == "Naked mole-rat" else 6
    prepared, species_mask = collect_level0(species, list(level0_source(species)), minimum)
    table, excluded = level0_table(species, prepared, species_mask)
    fidelity = fidelity_checks(species, prepared, species_mask)
    slug = SLUGS[species]
    table.to_csv(scratch_dir / f"level0-{slug}.csv", index=False)
    excluded.to_csv(scratch_dir / f"level0-{slug}-excluded.csv", index=False)
    if len(excluded):
        print(f"{species}: excluded {len(excluded)} units with no identifiable age slope:",
              flush=True)
        print(excluded.to_string(index=False), flush=True)
    fidelity.to_csv(scratch_dir / f"level0-{slug}-fidelity.csv", index=False)
    print(fidelity.to_string(index=False), flush=True)
    primary = table[table["pipeline_gene_set"]]
    n_with_sig = int(((primary["n_significant_q05"] > 0).fillna(False)).sum())
    print(f"{species}: {len(primary)} units; units with >=1 significant gene = "
          f"{n_with_sig}; median significant = "
          f"{primary['n_significant_q05'].median()}", flush=True)
    if len(primary) != EXPECTED_UNITS[species]:
        raise ValueError(f"{species}: expected {EXPECTED_UNITS[species]} units, got {len(primary)}")
    if not fidelity["passed"].all():
        raise ValueError(f"{species}: Level 1 fidelity check failed")
    print(f"Wrote {scratch_dir / f'level0-{slug}.csv'} ({len(table)} rows)", flush=True)


# ---------------------------------------------------------------------------
# Level 2
# ---------------------------------------------------------------------------

def fly_pair_id(donor):
    pair = re.sub(r"_(?:head|body)_", "_", str(donor), count=1)
    if pair == str(donor):
        raise ValueError(f"Cannot derive head/body pair identity from {donor}")
    return pair


def cross_pair_mask(group_labels):
    """Boolean mask over np.triu_indices(n, 1) order: True where groups differ."""
    labels = np.asarray(group_labels)
    first, second = np.triu_indices(len(labels), 1)
    return labels[first] != labels[second]


def common_genes(gene_lists):
    indexes = [pd.Index(np.asarray(genes).astype(str)) for genes in gene_lists]
    for index in indexes:
        if index.has_duplicates:
            raise ValueError("Duplicate gene names prevent alignment")
    common = indexes[0]
    for index in indexes[1:]:
        common = common.intersection(index, sort=False)
    if common.empty:
        raise ValueError("No common genes across strata")
    common = pd.Index(sorted(common))
    return common, [index.get_indexer(common) for index in indexes]


TECHNICAL_GENE_RE = re.compile(r"^(MT[-:]|RP[LS]|MRP[LS])", re.IGNORECASE)


def technical_gene_names(genes):
    """Genes matching the technical set (mitochondrial, ribosomal, mitochondrial
    ribosomal): names starting with MT-/MT:, RPL/RPS, or MRPL/MRPS."""
    return np.array([bool(TECHNICAL_GENE_RE.match(str(gene))) for gene in genes])


def prepare_cross(strata, species, min_donors=6):
    genes, indexers = common_genes([stratum[2] for stratum in strata])
    records = []
    per_group_masks = []
    for (tissue, cohort, _, units), indexer in zip(strata, indexers):
        units = [unit for unit in units if len(unit["metadata"]) >= min_donors]
        n_before = len(units)
        units = [unit for unit in units if age_identifiable(unit)]
        if len(units) != n_before:
            print(f"Level 2 {species} | {tissue} | {cohort}: dropped "
                  f"{n_before - len(units)} units with no identifiable age slope", flush=True)
        aligned = []
        for unit in units:
            expression = np.asarray(unit["expression"], dtype=np.float64)[:, indexer]
            aligned.append({"species": species, "tissue": tissue, "cohort": cohort,
                            "cell_type": unit["cell_type"], "metadata": unit["metadata"].copy(),
                            "expression": expression})
        per_group_masks.append(primary_mask(aligned))
        records.extend(aligned)
    masks = {"primary": primary_mask(records),
             "all_strata_primary": np.logical_and.reduce(per_group_masks)}
    if species == "Fruit fly":
        technical = technical_gene_names(genes)
        masks["technical_masked"] = masks["primary"] & ~technical
        print(f"Level 2 {species}: technical genes masked: {int(technical.sum())} of "
              f"{len(genes)} common genes", flush=True)
    return records, genes, masks


def donor_key(record, raw_donor, scheme):
    if scheme == "paired_shared":
        return fly_pair_id(raw_donor)
    if scheme == "independent_by_tissue":
        return f"{record['tissue']}|{raw_donor}"
    return f"{record['cohort']}|{raw_donor}"


def donor_table(records, scheme):
    values = {}
    for record in records:
        for row in record["metadata"].itertuples():
            key = donor_key(record, row.donor, scheme)
            entry = (float(row.age), float(row.sex_male), record["cohort"],
                     record["tissue"] if scheme == "independent_by_tissue" else "")
            if key in values and values[key] != entry:
                raise ValueError(f"Inconsistent donor metadata for {key}: {values[key]} vs {entry}")
            values[key] = entry
    keys = sorted(values)
    return pd.DataFrame({
        "donor": keys,
        "age": [values[key][0] for key in keys],
        "sex_male": [values[key][1] for key in keys],
        "cohort": [values[key][2] for key in keys],
        "tissue": [values[key][3] for key in keys],
    })


def permutation_columns(donors, scheme):
    columns = []
    if scheme == "independent_by_cohort":
        columns.append("cohort")
    elif scheme == "independent_by_tissue":
        columns.append("tissue")
    if donors["sex_male"].nunique() > 1:
        columns.append("sex_male")
    return columns


def permute_ages(donors, columns, rng):
    permuted = donors["age"].to_numpy(dtype=float, copy=True)
    if not columns:
        return rng.permutation(permuted)
    for _, positions in sorted(donors.groupby(columns, sort=True).indices.items()):
        permuted[positions] = rng.permutation(permuted[positions])
    return permuted


def cross_null(records, masks, cross, scheme, permutations, seed, report_every=100):
    """Observed and permuted mean cross-group Pearson for each gene-set mask.

    Slopes use the Frisch-Waugh-Lovell form on the clean design (jitter UMI
    dropped): slope = (r_age . Y) / ||r_age||^2, where r_age is the age
    column residualised on the other design columns. This is exactly the OLS
    slope and never forms pinv(X.T @ X). The non-age columns are fixed across
    permutations, so their QR is precomputed once per record. Per iteration,
    a unit whose permuted age column is unidentifiable (age residualised to
    ~0, or no residual degrees of freedom) is masked out of that iteration's
    mean and recorded; the observed iteration must have every prepared unit
    identifiable. Per-pair null draws for the primary mask are returned as a
    full-triu float32 matrix (NaN where a pair was masked) so decompositions
    can be computed from the same draws.
    """
    donors = donor_table(records, scheme)
    columns = permutation_columns(donors, scheme)
    index = {key: position for position, key in enumerate(donors["donor"])}
    union = np.logical_or.reduce(list(masks.values()))
    within_union = {name: mask[union] for name, mask in masks.items()}
    n_records = len(records)
    full_first, full_second = np.triu_indices(n_records, 1)
    n_full_pairs = len(full_first)
    full_pos = np.full((n_records, n_records), -1, dtype=int)
    full_pos[full_first, full_second] = np.arange(n_full_pairs)
    cross_pos = np.where(cross)[0]
    for record in records:
        record["expr_union"] = record["expression"][:, union]
        record["donor_idx"] = np.array([index[donor_key(record, donor, scheme)]
                                        for donor in record["metadata"]["donor"]])
    prepped = []
    for record in records:
        design0 = clean_unit_design(record["metadata"])
        others = np.delete(design0, 1, axis=1)
        q, _ = np.linalg.qr(others)
        prepped.append({"q": q, "rank_others": int(np.linalg.matrix_rank(others)),
                        "n": len(record["donor_idx"])})
    rng = np.random.default_rng(seed)
    observed, observed_values, observed_within = {}, {}, {}
    observed_pairs = None
    nulls = {name: [] for name in masks}
    null_pairs = []
    masked_iterations = 0
    for iteration in range(permutations + 1):
        ages = (donors["age"].to_numpy(dtype=float) if iteration == 0
                else permute_ages(donors, columns, rng))
        slopes = np.full((n_records, int(union.sum())), np.nan)
        ident = np.zeros(n_records, dtype=bool)
        for ri, record in enumerate(records):
            age_vec = ages[record["donor_idx"]]
            prep = prepped[ri]
            if np.ptp(age_vec) == 0:
                continue
            r_age = age_vec - prep["q"] @ (prep["q"].T @ age_vec)
            rss_age = float(r_age @ r_age)
            if rss_age <= 1e-10 * float(np.sum((age_vec - age_vec.mean()) ** 2)):
                continue
            if prep["n"] - (prep["rank_others"] + 1) <= 0:
                continue
            ident[ri] = True
            slopes[ri] = (r_age @ record["expr_union"]) / rss_age
        if iteration == 0 and not ident.all():
            raise ValueError("Observed iteration has unidentifiable units; "
                             "prepare_cross should have dropped them")
        if not ident.all():
            masked_iterations += 1
        idx = np.where(ident)[0]
        sub_first, sub_second = np.triu_indices(len(idx), 1)
        full_i, full_j = idx[sub_first], idx[sub_second]
        pos = full_pos[full_i, full_j]
        cross_sub = cross[pos]
        for name, mask in within_union.items():
            values = mean_pearson(slopes[idx], mask)
            value = float(values[cross_sub].mean())
            if iteration == 0:
                observed[name] = value
                observed_values[name] = values[cross]
                observed_within[name] = float(values[~cross].mean()) if (~cross).any() else np.nan
                if name == "primary":
                    full_values = np.full(n_full_pairs, np.nan)
                    full_values[pos] = values
                    observed_pairs = full_values
            else:
                nulls[name].append(value)
                if name == "primary":
                    full_values = np.full(n_full_pairs, np.nan, dtype=np.float32)
                    full_values[pos] = values.astype(np.float32)
                    null_pairs.append(full_values)
        if iteration % report_every == 0:
            print(f"  {scheme}: permutation {iteration}/{permutations}", flush=True)
    return {"donors": donors, "columns": columns, "observed": observed,
            "observed_values": observed_values, "observed_within": observed_within,
            "nulls": {name: np.asarray(values) for name, values in nulls.items()},
            "observed_pairs": observed_pairs,
            "null_pairs": np.asarray(null_pairs, dtype=np.float32),
            "cross_pos": cross_pos, "masked_iterations": masked_iterations}


def summarize_null(observed, null):
    return {
        "observed_mean_pearson": observed,
        "joint_null_mean": float(null.mean()),
        "joint_null_sd": float(null.std(ddof=1)),
        "joint_null_ci95_low": float(np.quantile(null, 0.025)),
        "joint_null_ci95_high": float(np.quantile(null, 0.975)),
        "joint_p_greater": float((1 + (null >= observed).sum()) / (len(null) + 1)),
    }


def fly_pool_overlap(strata):
    by_tissue = {tissue: {str(donor) for unit in units for donor in unit["metadata"]["donor"]}
                 for tissue, _, _, units in strata}
    head, body = by_tissue["head"], by_tissue["body"]
    head_meta = {}
    body_meta = {}
    for tissue, _, _, units in strata:
        target = head_meta if tissue == "head" else body_meta
        for unit in units:
            for row in unit["metadata"].itertuples():
                target[fly_pair_id(row.donor)] = (row.age, row.sex_male)
    shared = set(head_meta) & set(body_meta)
    return {
        "n_head_pools": len(head), "n_body_pools": len(body),
        "n_identical_raw_ids": len(head & body),
        "n_inferred_pairs": len(shared),
        "pairs_match_age_sex": all(head_meta[key] == body_meta[key] for key in shared),
        "complete_one_to_one": set(head_meta) == set(body_meta) and len(shared) == len(head) == len(body),
    }


LEVEL2_CONFIG = {
    "Human": {"group": "cohort", "scope": "cross_cohort", "expected_pairs": 6675,
              "schemes": [("independent_by_cohort", 20261613 + 7, "primary")]},
    "Fruit fly": {"group": "tissue", "scope": "cross_stratum_head_body", "expected_pairs": 2016,
                  "schemes": [("paired_shared", 20261993 + 2, "primary"),
                              ("independent_by_tissue", 20261993 + 3, "sensitivity")]},
}
SCHEME_TEXT = {
    "independent_by_cohort": ("one permuted age per donor shared across its cell types; ages permuted "
                              "independently within each cohort (and sex where variable); donors keyed "
                              "cohort|donor"),
    "paired_shared": ("one permuted age per inferred head/body pool pair, shared across both organs and "
                      "all cell types; permuted within sex"),
    "independent_by_tissue": ("ages permuted independently within each organ (head, body) and sex; one "
                              "age per pool shared across that organ's cell types"),
}


HUMAN_ORGAN = {
    "HBCC_DLPFC": "Brain", "MSSM_DLPFC": "Brain",
    "Lung_CellRef": "Lung", "Lung_Integrated": "Lung",
    "Kidney_LivingDonor": "Kidney", "Kidney_SCR_v15": "Kidney",
    "ImmuneAging_BAL": "BAL",
}


def subset_null_summary(observed_values, null_matrix, selector):
    """summarize_null over a subset of pairs. observed_values and null_matrix
    share pair order; selector is boolean over pairs. Masked (NaN) pairs are
    ignored per iteration."""
    observed = float(np.nanmean(observed_values[selector]))
    null = np.nanmean(null_matrix[:, selector], axis=1)
    return summarize_null(observed, null)


def group_weighted_null_summary(observed_values, null_matrix, group_labels):
    """Mean over groups of the within-group mean, per iteration. Used for
    equal cohort-pair weighting: each of the 21 cohort pairs counts 1/21,
    so an 8-donor cohort cannot dominate by pair count. Group labels may be
    tuples (cohort pairs); they are flattened to "|"-joined strings first,
    because comparing an object array of tuples with == broadcasts the
    tuple instead of matching elementwise (and pd.Series misreads each
    tuple as a 2D row)."""
    keys = np.array(["|".join(g) if isinstance(g, (tuple, list)) else str(g)
                     for g in group_labels])
    unique = np.unique(keys)
    observed = float(np.mean([np.nanmean(observed_values[keys == g])
                              for g in unique]))
    null = np.array([np.mean([np.nanmean(null_matrix[i][keys == g])
                              for g in unique]) for i in range(null_matrix.shape[0])])
    return summarize_null(observed, null)


def human_level2_decompositions(pair_frame, null_pairs, cross_pos):
    """Human cross-cohort decompositions from the same null draws.

    pair_frame rows are the cross-cohort pairs in cross_pos order;
    null_pairs is the (n_permutations, n_full_pairs) primary-mask matrix.
    Returns (detail_rows, cohort_pair_rows). All rows carry the same
    permutation null as the primary human row; they are decompositions of
    it, not independent tests.
    """
    observed = np.asarray(pair_frame["pearson"], dtype=float)
    null_cross = np.asarray(null_pairs[:, cross_pos], dtype=float)
    cohort_1 = pair_frame["cohort_1"].to_numpy()
    cohort_2 = pair_frame["cohort_2"].to_numpy()
    organ_1 = np.array([HUMAN_ORGAN[c] for c in cohort_1])
    organ_2 = np.array([HUMAN_ORGAN[c] for c in cohort_2])
    pair_key = np.array(["|".join(sorted(p)) for p in zip(cohort_1, cohort_2)])
    detail_rows = []

    def add(scope, summary, n_pairs):
        row = {"scope": scope}
        row.update(summary)
        row["n_pairs"] = int(n_pairs)
        detail_rows.append(row)

    add("cross_cohort_equal_cohort_pair_weighted",
        group_weighted_null_summary(observed, null_cross, pair_key), len(observed))
    strict = (organ_1 != organ_2) & (organ_1 != "BAL") & (organ_2 != "BAL")
    add("cross_cohort_strict_cross_organ", subset_null_summary(observed, null_cross, strict),
        int(strict.sum()))
    same_organ = (organ_1 == organ_2) & (cohort_1 != cohort_2)
    add("cross_cohort_same_organ", subset_null_summary(observed, null_cross, same_organ),
        int(same_organ.sum()))
    same_type = same_organ & (pair_frame["cell_type_1"].to_numpy()
                              == pair_frame["cell_type_2"].to_numpy())
    add("cross_cohort_same_organ_same_cell_type",
        subset_null_summary(observed, null_cross, same_type), int(same_type.sum()))
    add("cross_cohort_same_organ_different_cell_type",
        subset_null_summary(observed, null_cross, same_organ & ~same_type),
        int((same_organ & ~same_type).sum()))
    bal = (cohort_1 == "ImmuneAging_BAL") | (cohort_2 == "ImmuneAging_BAL")
    add("cross_cohort_bal_involved", subset_null_summary(observed, null_cross, bal),
        int(bal.sum()))
    for cohort in sorted(set(cohort_1) | set(cohort_2)):
        keep = (cohort_1 != cohort) & (cohort_2 != cohort)
        add(f"cross_cohort_leave_{cohort}_out",
            subset_null_summary(observed, null_cross, keep), int(keep.sum()))

    cohort_pair_rows = []
    for key in sorted(set(pair_key.tolist())):
        selector = pair_key == key
        summary = subset_null_summary(observed, null_cross, selector)
        cohort_a, cohort_b = key.split("|")
        row = {"cohort_1": cohort_a, "cohort_2": cohort_b}
        row.update(summary)
        row["n_pairs"] = int(selector.sum())
        cohort_pair_rows.append(row)
    return detail_rows, cohort_pair_rows


def run_level2(species, permutations, scratch_dir):
    config = LEVEL2_CONFIG[species]
    strata = list(joint.cached_strata(species))
    overlap = None
    if species == "Fruit fly":
        overlap = fly_pool_overlap(strata)
        print(f"Fly pool overlap: {overlap}", flush=True)
        if not (overlap["complete_one_to_one"] and overlap["pairs_match_age_sex"]):
            raise ValueError("Head/body counterpart labels are not a complete 1:1 age/sex match")
    records, genes, masks = prepare_cross(strata, species)
    labels = np.array([record[config["group"]] for record in records])
    cross = cross_pair_mask(labels)
    n_units = len(records)
    print(f"{species}: {n_units} units, {int(cross.sum())} cross pairs, {int((~cross).sum())} "
          f"within pairs, common genes {len(genes)}, primary genes {int(masks['primary'].sum())}, "
          f"all-strata genes {int(masks['all_strata_primary'].sum())}", flush=True)
    if int(cross.sum()) != config["expected_pairs"]:
        raise ValueError(f"{species}: expected {config['expected_pairs']} cross pairs, got {int(cross.sum())}")
    n_libraries = len({f"{record['cohort']}|{record['tissue']}|{donor}" for record in records
                       for donor in record["metadata"]["donor"]})
    details, pair_frame = [], None
    primary_result = None
    for scheme, seed, role in config["schemes"]:
        print(f"{species}: running {scheme} null, seed {seed}", flush=True)
        result = cross_null(records, masks, cross, scheme, permutations, seed)
        if role == "primary":
            primary_result = result
        for gene_set in masks:
            row = {"species": species, "scope": config["scope"], "null_scheme": scheme,
                   "null_role": role, "gene_set": gene_set}
            row.update(summarize_null(result["observed"][gene_set], result["nulls"][gene_set]))
            row.update({
                "n_pairs": int(cross.sum()), "n_units": n_units,
                "n_donors_or_pools": n_libraries,
                "n_exchangeable_units": len(result["donors"]),
                "n_permutations": permutations,
                "n_common_genes": len(genes), "n_genes_tested": int(masks[gene_set].sum()),
                "observed_within_group_mean_same_mask": result["observed_within"][gene_set],
                "n_within_group_pairs": int((~cross).sum()),
                "random_seed": seed, "permutation_columns": "+".join(result["columns"]),
                "permutation_scheme": SCHEME_TEXT[scheme],
                "detection_rule": ("detected in >=50% of donor pseudobulks in >=50% of all "
                                   f"{n_units} units (common genes)" if gene_set == "primary" else
                                   "passes the per-stratum primary rule in every stratum (common genes)"),
                "fly_pool_overlap": json.dumps(overlap) if overlap else "",
            })
            details.append(row)
            print(pd.Series(row).to_string(), flush=True)
        if pair_frame is None:
            first, second = np.triu_indices(n_units, 1)
            first, second = first[cross], second[cross]
            pair_frame = pd.DataFrame({
                "species": species, "scope": config["scope"],
                "tissue_1": [records[i]["tissue"] for i in first],
                "cohort_1": [records[i]["cohort"] for i in first],
                "cell_type_1": [records[i]["cell_type"] for i in first],
                "tissue_2": [records[i]["tissue"] for i in second],
                "cohort_2": [records[i]["cohort"] for i in second],
                "cell_type_2": [records[i]["cell_type"] for i in second],
                "pearson": result["observed_values"]["primary"],
            })
    if species == "Human":
        decomp_rows, cohort_pair_rows = human_level2_decompositions(
            pair_frame, primary_result["null_pairs"], primary_result["cross_pos"])
        for decomp in decomp_rows:
            row = {"species": species, "scope": decomp["scope"],
                   "null_scheme": config["schemes"][0][0], "null_role": "decomposition",
                   "gene_set": "primary"}
            row.update({k: decomp[k] for k in
                        ("observed_mean_pearson", "joint_null_mean", "joint_null_sd",
                         "joint_null_ci95_low", "joint_null_ci95_high", "joint_p_greater")})
            row.update({
                "n_pairs": decomp["n_pairs"], "n_units": n_units,
                "n_donors_or_pools": n_libraries,
                "n_exchangeable_units": len(primary_result["donors"]),
                "n_permutations": permutations,
                "n_common_genes": len(genes),
                "n_genes_tested": int(masks["primary"].sum()),
                "observed_within_group_mean_same_mask": np.nan,
                "n_within_group_pairs": int((~cross).sum()),
                "random_seed": config["schemes"][0][1],
                "permutation_columns": "+".join(primary_result["columns"]),
                "permutation_scheme": SCHEME_TEXT[config["schemes"][0][0]],
                "detection_rule": ("detected in >=50% of donor pseudobulks in >=50% of all "
                                   f"{n_units} units (common genes)"),
                "fly_pool_overlap": "",
            })
            details.append(row)
            print(pd.Series(row).to_string(), flush=True)
        cohort_path = RESULTS / "followup-mature-hierarchy-level2-human-cohort-pairs.csv"
        pd.DataFrame(cohort_pair_rows).to_csv(cohort_path, index=False)
        print(f"Wrote {cohort_path} ({len(cohort_pair_rows)} cohort pairs)", flush=True)
    slug = SLUGS[species]
    pd.DataFrame(details).to_csv(scratch_dir / f"level2-{slug}-details.csv", index=False)
    pair_path = RESULTS / f"followup-mature-hierarchy-level2-{slug}-pairs.csv"
    pair_frame.to_csv(pair_path, index=False)
    print(f"Wrote {scratch_dir / f'level2-{slug}-details.csv'} and {pair_path} "
          f"({len(pair_frame)} pairs)", flush=True)


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------

def single_stratum_rows():
    rows = []
    worm = pd.read_csv(RESULTS / "followup-mature-concordance-worm.csv").iloc[0]
    rows.append({"species": "C. elegans", "scope": "single_stratum_level1_is_level2",
                 "observed_mean_pearson": worm["observed_mean_pearson"],
                 "joint_null_mean": worm["joint_null_mean"], "joint_null_sd": worm["joint_null_sd"],
                 "joint_p_greater": worm["joint_p_greater"], "n_pairs": worm["n_pairs"],
                 "n_units": worm["n_units"], "n_donors_or_pools": worm["n_donors"],
                 "n_permutations": worm["n_permutations_joint"]})
    fish = pd.read_csv(RESULTS / "reviewer-rev2-joint-concordance-zebrafish.csv").iloc[0]
    rows.append({"species": "Zebrafish", "scope": "single_stratum_level1_is_level2",
                 "observed_mean_pearson": fish["observed_mean_pearson"],
                 "joint_null_mean": fish["joint_null_mean"], "joint_null_sd": fish["joint_null_sd"],
                 "joint_p_greater": fish["joint_p_greater"], "n_pairs": fish["n_pairs"],
                 "n_units": fish["n_units"], "n_donors_or_pools": fish["n_donors"],
                 "n_permutations": fish["n_permutations_joint"]})
    rev = pd.read_csv(RESULTS / "reviewer-rev2-concordance.csv")
    nmr = rev[rev["species"] == "Naked mole-rat"].iloc[0]
    perms = pd.read_csv(RESULTS / "reviewer-phase1-universality-permutations-nmr.csv")
    rows.append({"species": "Naked mole-rat", "scope": "single_stratum_level1_is_level2",
                 "observed_mean_pearson": nmr["observed_mean_pearson"],
                 "joint_null_mean": nmr["joint_null_mean"], "joint_null_sd": nmr["joint_null_sd"],
                 "joint_p_greater": nmr["joint_p_greater"], "n_pairs": nmr["n_pairs"],
                 "n_units": int(perms["n_units"].iloc[0]), "n_donors_or_pools": nmr["n_donors"],
                 "n_permutations": int((perms["iteration"] != 0).sum())})
    return rows


def assemble(scratch_dir):
    level0 = pd.concat([pd.read_csv(scratch_dir / f"level0-{slug}.csv") for slug in SLUGS.values()],
                       ignore_index=True)
    for slug in SLUGS.values():
        fidelity = pd.read_csv(scratch_dir / f"level0-{slug}-fidelity.csv")
        if not fidelity["passed"].all():
            raise ValueError(f"Fidelity check failed for {slug}")
    primary = level0[level0["pipeline_gene_set"]].copy()
    counts = primary.groupby("species").size().to_dict()
    if counts != EXPECTED_UNITS or len(primary) != EXPECTED_LEVEL0_TOTAL:
        raise ValueError(f"Unexpected Level 0 unit counts: {counts}")
    identified = primary[primary["n_significant_q05"].notna()]
    if not np.allclose(identified["frac_significant"],
                       identified["n_significant_q05"] / identified["n_genes_tested"],
                       rtol=0, atol=1e-15):
        raise ValueError("Level 0 fractions do not match counts")
    order = list(SLUGS)
    primary["species_order"] = primary["species"].map(order.index)
    primary = primary.sort_values(["species_order", "tissue_cohort", "cell_type"])
    level0_path = RESULTS / "followup-mature-hierarchy-level0.csv"
    primary[LEVEL0_COLUMNS].to_csv(level0_path, index=False)
    sensitivity = level0[~level0["pipeline_gene_set"]].drop(columns="pipeline_gene_set")
    sensitivity.to_csv(RESULTS / "followup-mature-hierarchy-level0-gene-set-sensitivity.csv", index=False)

    mouse = pd.read_csv(RESULTS / "followup-mature-concordance-mouse.csv")
    row = mouse[(mouse["scope"] == "cross_tissue")
                & (mouse["gene_set"] == "primary_detection_filtered")].iloc[0]
    rows = [{"species": "Mouse", "scope": "cross_tissue",
             "observed_mean_pearson": row["observed_mean_pearson"],
             "joint_null_mean": row["joint_null_mean"], "joint_null_sd": row["joint_null_sd"],
             "joint_p_greater": row["joint_p_greater"], "n_pairs": row["n_pairs"],
             "n_units": row["n_units"], "n_donors_or_pools": row["n_donors"],
             "n_permutations": row["n_permutations_joint"]}]
    details = pd.concat([pd.read_csv(scratch_dir / f"level2-{slug}-details.csv")
                         for slug in ("human", "fly")], ignore_index=True)
    for record in details.itertuples():
        # level2.csv holds strictly the 6 primary rows (mouse, human, fly,
        # worm, zebrafish, naked mole-rat). Sensitivity and decomposition rows
        # stay in level2-details.csv only.
        if record.gene_set != "primary" or record.null_role != "primary":
            continue
        rows.append({"species": record.species, "scope": record.scope,
                     "observed_mean_pearson": record.observed_mean_pearson,
                     "joint_null_mean": record.joint_null_mean, "joint_null_sd": record.joint_null_sd,
                     "joint_p_greater": record.joint_p_greater, "n_pairs": record.n_pairs,
                     "n_units": record.n_units, "n_donors_or_pools": record.n_donors_or_pools,
                     "n_permutations": record.n_permutations})
    rows.extend(single_stratum_rows())
    level2 = pd.DataFrame(rows, columns=LEVEL2_COLUMNS)
    if len(level2) != 6:
        raise ValueError(f"level2.csv must have exactly 6 primary rows, got {len(level2)}")
    level2_path = RESULTS / "followup-mature-hierarchy-level2.csv"
    level2.to_csv(level2_path, index=False)
    details.to_csv(RESULTS / "followup-mature-hierarchy-level2-details.csv", index=False)
    fidelity = pd.concat([pd.read_csv(scratch_dir / f"level0-{slug}-fidelity.csv")
                          for slug in SLUGS.values()], ignore_index=True)
    fidelity.to_csv(RESULTS / "followup-mature-hierarchy-level0-fidelity.csv", index=False)
    print(f"Wrote {level0_path} ({len(primary)} rows) and {level2_path} ({len(level2)} rows)", flush=True)
    print(level2.to_string(index=False), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("mode", choices=["level0", "level2", "assemble"])
    parser.add_argument("--species", choices=list(SLUGS))
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--scratch-dir", type=Path, default=SCRATCH)
    args = parser.parse_args()
    args.scratch_dir.mkdir(parents=True, exist_ok=True)
    if args.mode == "assemble":
        assemble(args.scratch_dir)
        return
    if args.species is None:
        parser.error("--species is required for level0 and level2")
    if args.mode == "level0":
        run_level0(args.species, args.scratch_dir)
    else:
        if args.species not in LEVEL2_CONFIG:
            parser.error("Level 2 is new only for Human and Fruit fly")
        run_level2(args.species, args.permutations, args.scratch_dir)


if __name__ == "__main__":
    main()
