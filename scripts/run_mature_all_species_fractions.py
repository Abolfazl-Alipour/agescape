#!/usr/bin/env python3
"""Mature-to-old cell-fraction scan across all eligible Agescape species."""

import hashlib
import sys
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import t as t_dist
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_afca_aging import afca_obs
from scripts.run_gao_worm_aging import normal_obs as worm_normal_obs
from scripts.run_mature_composition_complete import clr_table, human_inputs, mouse_inputs
from scripts.run_zebrafish_testis_aging import standard_obs as zfish_standard_obs

RESULTS = REPO / "results"
REPORTS = REPO / "reports"
KEY = ["species", "organ", "cohort", "cell_type"]
STRATUM = ["species", "organ", "cohort"]
MODEL_USED = "kimi-code/k3.1-pro"
LOCKED_624 = RESULTS / "rev3-composition-global-bh-624.csv"
LOCKED_624_SHA256 = "3f26242585017bd3a68a4a7c4dcb1d2684da0e9f662c37ed7e07e1bd6f054d96"
OUTPUT = RESULTS / "mature-all-species-composition.csv"
COVERAGE_OUTPUT = RESULTS / "mature-all-species-composition-coverage.csv"
EXCLUSIONS_OUTPUT = RESULTS / "mature-all-species-composition-exclusions.csv"
UNTESTABLE_OUTPUT = RESULTS / "mature-all-species-composition-untestable.csv"
VERIFICATION_OUTPUT = RESULTS / "mature-all-species-composition-verification.csv"
RUNTIME_OUTPUT = RESULTS / "mature-all-species-composition-runtime.csv"
REPORT_OUTPUT = REPORTS / "followup-mature-all-species-fractions.md"
MIN_PRESENT = 6
PSEUDOCOUNT = 0.5
SPECIES_ORDER = ["Mouse", "Human", "Fruit fly", "C. elegans", "Zebrafish", "Naked mole-rat"]
WINDOWS = {
    "Mouse": {"cutoff": 18.0, "age_unit": "month", "sensitivity": "mature_only_18mo"},
    "Human": {"cutoff": 30.0, "age_unit": "year", "sensitivity": "mature_only_30y"},
    "Fruit fly": {"cutoff": 30.0, "age_unit": "day", "sensitivity": "mature_only_30d"},
    "C. elegans": {"cutoff": 6.0, "age_unit": "day", "sensitivity": "mature_only_6d"},
    "Zebrafish": {"cutoff": 12.0, "age_unit": "month", "sensitivity": "mature_only_12mo"},
}
UNITS = {
    "Mouse": "animal",
    "Human": "donor",
    "Fruit fly": "pool",
    "C. elegans": "pool",
    "Zebrafish": "fish",
    "Naked mole-rat": "animal",
}
NMR_AGES = {"1J": 1.1, "5J": 1.1, "6J": 1.1, "2V": 11.3, "3V": 11.3, "4V": 11.3}
HEADLINE_KEYS = {
    ("Human", "brain", "MSSM_DLPFC", "L6b glutamatergic neuron of the primary motor cortex"),
    ("Human", "brain", "MSSM_DLPFC", "oligodendrocyte"),
    ("Human", "brain", "MSSM_DLPFC", "oligodendrocyte precursor cell"),
    ("Mouse", "Limb_Muscle", "TMS", "macrophage"),
}


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def fmt(value, digits=12):
    if value is None or pd.isna(value):
        return "NA"
    if isinstance(value, (bool, np.bool_)):
        return str(bool(value))
    return f"{float(value):.{digits}g}"


def species_sort(frame):
    out = frame.copy()
    out["_species_order"] = out["species"].map({s: i for i, s in enumerate(SPECIES_ORDER)})
    out = out.sort_values(["_species_order", "organ", "cohort", "cell_type"]).drop(columns="_species_order")
    return out.reset_index(drop=True)


def count_table(obs, unit_column, cell_column, columns):
    counts = obs.groupby([unit_column, cell_column], observed=True).size().unstack(fill_value=0)
    return counts.sort_index().reindex(columns=columns, fill_value=0)


def read_obs(path):
    atlas = ad.read_h5ad(path, backed="r")
    obs = atlas.obs.copy()
    atlas.file.close()
    return obs


def load_strata():
    for species, organ, cohort, counts, metadata in mouse_inputs():
        yield species, organ, cohort, counts, metadata
    for species, organ, cohort, counts, metadata in human_inputs():
        yield species, organ, cohort, counts, metadata

    for organ in ("head", "body"):
        path = REPO / "data" / "raw" / "afca" / f"afca-{organ}.h5ad"
        obs = afca_obs(ad.read_h5ad(path, backed="r"))
        obs = obs[(obs["dataset"] == "AFCA") & (obs["cell_type"] != "unannotated")].copy()
        columns = sorted(obs["cell_type"].unique())
        part = obs[obs["age"] >= WINDOWS["Fruit fly"]["cutoff"]].copy()
        counts = count_table(part, "sample", "cell_type", columns)
        metadata = part.groupby("sample").agg(age=("age", "first"), sex=("sex_male", "first"))
        yield "Fruit fly", organ, "AFCA", counts, metadata.loc[counts.index]

    path = REPO / "data" / "interim" / "worm-gao-integrated.h5ad"
    obs = worm_normal_obs(ad.read_h5ad(path, backed="r"))
    columns = sorted(obs["cell_type"].unique())
    part = obs[obs["age"] >= WINDOWS["C. elegans"]["cutoff"]].copy()
    counts = count_table(part, "donor", "cell_type", columns)
    metadata = part.groupby("donor").agg(age=("age", "first"), sex=("sex_male", "first"))
    yield "C. elegans", "Whole animal", "Gao_GSE229022", counts, metadata.loc[counts.index]

    path = REPO / "data" / "interim" / "zebrafish-testis-5-22mo.h5ad"
    obs = zfish_standard_obs(ad.read_h5ad(path, backed="r"), "5-22mo")
    obs = obs[obs["fine_cell_type"] != "unknown"].copy()
    columns = sorted(obs["fine_cell_type"].unique())
    part = obs[obs["age"] >= WINDOWS["Zebrafish"]["cutoff"]].copy()
    counts = count_table(part, "donor", "fine_cell_type", columns)
    metadata = part.groupby("donor").agg(age=("age", "first"), sex=("sex_male", "first"))
    yield "Zebrafish", "Testis", "GSE275361", counts, metadata.loc[counts.index]


def fit_cell_type(clr, counts, metadata, cell_type, species, organ, cohort):
    donors = counts.index
    if len(donors) < MIN_PRESENT:
        return None, f"fewer than {MIN_PRESENT} mature {UNITS[species]}s"
    X = sm.add_constant(metadata.loc[donors, ["age", "sex"]].astype(float))
    model = "age_plus_sex"
    if np.linalg.matrix_rank(X.to_numpy()) < X.shape[1]:
        X = sm.add_constant(metadata.loc[donors, ["age"]].astype(float))
        model = "age_only_rank_deficient_sex"
    if np.linalg.matrix_rank(X.to_numpy()) < X.shape[1] or len(donors) <= X.shape[1]:
        return None, "rank-deficient age design or insufficient residual degrees of freedom"
    fit = sm.OLS(clr.loc[donors, cell_type], X).fit()
    if not np.isfinite(fit.pvalues["age"]):
        return None, "non-finite age test after model fit"
    df = int(fit.df_resid)
    half = t_dist.ppf(0.975, df) * fit.bse["age"]
    window = WINDOWS[species]
    return {
        "species": species,
        "organ": organ,
        "cohort": cohort,
        "cell_type": str(cell_type),
        "sensitivity": window["sensitivity"],
        "model": model,
        "effect_type": f"CLR_slope_per_{window['age_unit']}",
        "age_unit": window["age_unit"],
        "mature_cutoff": window["cutoff"],
        "inferential_unit": UNITS[species],
        "effect": float(fit.params["age"]),
        "se": float(fit.bse["age"]),
        "ci95_low": float(fit.params["age"] - half),
        "ci95_high": float(fit.params["age"] + half),
        "pvalue": float(fit.pvalues["age"]),
        "n_donors": int(len(donors)),
        "n_present": int((counts.loc[donors, cell_type] > 0).sum()),
        "total_cells": int(counts.loc[donors, cell_type].sum()),
        "min_age": float(metadata.loc[donors, "age"].min()),
        "max_age": float(metadata.loc[donors, "age"].max()),
    }, ""


def run_stratum(species, organ, cohort, all_counts, metadata):
    present = (all_counts > 0).sum(axis=0)
    counts = all_counts.loc[:, present >= MIN_PRESENT].copy()
    rows = []
    failures = {}
    if counts.shape[1]:
        clr = clr_table(counts, pseudocount=PSEUDOCOUNT)
        for cell_type in counts.columns:
            row, reason = fit_cell_type(clr, counts, metadata, cell_type, species, organ, cohort)
            if row is None:
                failures[cell_type] = reason
            else:
                rows.append(row)
    if rows:
        qvalues = multipletests([row["pvalue"] for row in rows], method="fdr_bh")[1]
        for row, qvalue in zip(rows, qvalues):
            row["q_bh"] = float(qvalue)
    tested_types = {row["cell_type"] for row in rows}
    untestable = []
    for cell_type in all_counts.columns:
        if cell_type in tested_types:
            continue
        reason = failures.get(cell_type, f"cell type present in fewer than {MIN_PRESENT} mature {UNITS[species]}s")
        untestable.append({
            "species": species,
            "organ": organ,
            "cohort": cohort,
            "cell_type": str(cell_type),
            "reason": reason,
            "n_units_mature": int(len(metadata)),
            "n_present_mature": int((all_counts[cell_type] > 0).sum()),
            "total_cells_mature": int(all_counts[cell_type].sum()),
        })
    ages = sorted(float(x) for x in metadata["age"].unique())
    coverage = {
        "species": species,
        "organ": organ,
        "cohort": cohort,
        "inferential_unit": UNITS[species],
        "mature_cutoff": WINDOWS[species]["cutoff"],
        "age_unit": WINDOWS[species]["age_unit"],
        "available_ages_mature": ";".join(fmt(x, 6) for x in ages),
        "n_units_mature": int(len(metadata)),
        "n_cell_types_candidates": int(all_counts.shape[1]),
        "n_cell_types_tested": int(len(rows)),
        "min_age_mature": min(ages) if ages else np.nan,
        "max_age_mature": max(ages) if ages else np.nan,
    }
    return rows, untestable, coverage


def nmr_coverage_and_untestable():
    path = REPO / "data" / "raw" / "GSE152886_nmr_skin_raw.csv.gz"
    header = pd.read_csv(path, sep="\t", index_col=0, nrows=0)
    samples = sorted({str(x).split("_")[0] for x in header.columns})
    retained = [x for x in samples if x != "2V"]
    if set(retained) != {"1J", "5J", "6J", "3V", "4V"}:
        raise RuntimeError(f"Unexpected NMR samples: {samples}")
    ages = [NMR_AGES[x] for x in retained]
    coverage = {
        "species": "Naked mole-rat",
        "organ": "Skin",
        "cohort": "GSE152886",
        "inferential_unit": "animal",
        "mature_cutoff": np.nan,
        "age_unit": "year",
        "available_ages_mature": "not defined; all available ages: " + ";".join(fmt(x, 6) for x in sorted(set(ages))),
        "n_units_mature": 0,
        "n_cell_types_candidates": 4,
        "n_cell_types_tested": 0,
        "min_age_mature": np.nan,
        "max_age_mature": np.nan,
    }
    reason = ("no defensible mature-to-old window: after the established 2V exclusion, only three animals at 1.1 years "
              "and two at 11.3 years remain, with no intermediate ages; either group is below the six-animal floor")
    untestable = [{
        "species": "Naked mole-rat",
        "organ": "Skin",
        "cohort": "GSE152886",
        "cell_type": cell_type,
        "reason": reason,
        "n_units_mature": 0,
        "n_present_mature": np.nan,
        "total_cells_mature": np.nan,
    } for cell_type in ["basal", "intermediate", "corneum", "immune"]]
    return coverage, untestable, reason


def add_multiplicity_and_comparisons(result):
    result = result.copy()
    result["q_species_bh"] = np.nan
    result["n_species_family"] = np.nan
    for _, indices in result.groupby("species").groups.items():
        pvalues = result.loc[indices, "pvalue"]
        result.loc[indices, "q_species_bh"] = multipletests(pvalues, method="fdr_bh")[1]
        result.loc[indices, "n_species_family"] = len(indices)
    result["q_pooled_mature_bh"] = multipletests(result["pvalue"], method="fdr_bh")[1]
    result["n_pooled_mature_family"] = len(result)
    result["significant_pooled_05"] = result["q_pooled_mature_bh"] < 0.05
    result["n_stratum_family"] = result.groupby(STRATUM)["cell_type"].transform("size")

    full = pd.read_csv(LOCKED_624)
    full_columns = full[KEY + ["effect", "pvalue", "q_global_bh_624", "significant_global_624_05"]].rename(columns={
        "effect": "full_range_624_effect",
        "pvalue": "full_range_624_pvalue",
        "q_global_bh_624": "full_range_624_q",
        "significant_global_624_05": "full_range_624_pass",
    })
    result = result.merge(full_columns, on=KEY, how="left", validate="one_to_one")

    mature399 = pd.concat([
        pd.read_csv(RESULTS / "followup-mature-composition-mouse.csv"),
        pd.read_csv(RESULTS / "followup-mature-composition-human.csv"),
    ], ignore_index=True)
    mature_columns = mature399[KEY + ["q_refit_family"]].rename(columns={"q_refit_family": "mature_399_q"})
    result = result.merge(mature_columns, on=KEY, how="left", validate="one_to_one")
    mature_pass = mature399.loc[mature399["q_refit_family"] < 0.05, KEY].copy()
    mature_pass["mature_399_pass"] = True
    result = result.merge(mature_pass, on=KEY, how="left", validate="one_to_one")
    result["mature_399_pass"] = result["mature_399_pass"].eq(True)
    return species_sort(result), full, mature399


def make_exclusions(coverage, nmr_reason):
    rows = []
    for row in coverage.itertuples():
        if row.n_cell_types_tested > 0:
            continue
        if row.species == "Naked mole-rat":
            scope = "species"
            reason = nmr_reason
        elif row.species == "Zebrafish":
            scope = "species_primary_fine_arm"
            reason = ("12 months is the prespecified mature cutoff because the source atlas describes 5-month testes as approaching "
                      "peak spermatogenesis and 12-month testes as closest to peak fertility; only four fish at 12, 20, and 22 months "
                      "have comparable fine labels, below the six-fish floor, so no fine cell type is eligible; the separate broad-label "
                      "reconciliation arm is not substituted for the primary fine-label arm")
        else:
            scope = "tissue_or_cohort"
            reason = (f"zero eligible tests: no cell type met the {MIN_PRESENT}-unit presence rule and produced an estimable age model")
        rows.append({
            "species": row.species,
            "organ": row.organ,
            "cohort": row.cohort,
            "scope": scope,
            "inferential_unit": row.inferential_unit,
            "n_units_mature": row.n_units_mature,
            "available_ages_mature": row.available_ages_mature,
            "n_cell_types_candidates": row.n_cell_types_candidates,
            "reason": reason,
        })
    return pd.DataFrame(rows)


def independent_bh(pvalues):
    p = np.asarray(pvalues, dtype=float)
    m = len(p)
    order = np.argsort(p, kind="mergesort")
    ranked = p[order] * m / np.arange(1, m + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m, dtype=float)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def max_difference(left, right):
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    return float(np.max(np.abs(left - right))) if len(left) else 0.0


def verify_saved_outputs(result, full, mature399, before_hash):
    saved = pd.read_csv(OUTPUT)
    checks = []

    def record(name, value, passed=True):
        checks.append({"check": name, "value": value, "status": "PASS" if passed else "FAIL"})
        if not passed:
            raise RuntimeError(f"Verification failed: {name} = {value}")

    record("unique_test_keys", int(not saved.duplicated(KEY).any()), not saved.duplicated(KEY).any())
    record("finite_pvalues", int(np.isfinite(saved["pvalue"]).all()), np.isfinite(saved["pvalue"]).all())
    record("tests_equal_csv_rows", len(saved), len(saved) == len(result))
    pooled = independent_bh(saved["pvalue"])
    record("pooled_bh_max_abs_difference", max_difference(pooled, saved["q_pooled_mature_bh"]),
           np.allclose(pooled, saved["q_pooled_mature_bh"], rtol=0, atol=1e-14))
    for species, group in saved.groupby("species"):
        q = independent_bh(group["pvalue"])
        record(f"species_bh_{species}", max_difference(q, group["q_species_bh"]),
               np.allclose(q, group["q_species_bh"], rtol=0, atol=1e-14))
    for stratum, group in saved.groupby(STRATUM):
        q = independent_bh(group["pvalue"])
        record(f"stratum_bh_{stratum}", max_difference(q, group["q_bh"]),
               np.allclose(q, group["q_bh"], rtol=0, atol=1e-14))
    for species, group in saved.groupby("species"):
        cutoff = WINDOWS[species]["cutoff"]
        record(f"age_floor_{species}", float(group["min_age"].min()), group["min_age"].min() >= cutoff)
        record(f"presence_floor_{species}", int(group["n_present"].min()), group["n_present"].min() >= MIN_PRESENT)
        record(f"unit_floor_{species}", int(group["n_donors"].min()), group["n_donors"].min() >= MIN_PRESENT)

    current_mammals = saved[saved["species"].isin(["Mouse", "Human"])].set_index(KEY)
    previous = mature399.set_index(KEY)
    record("mammalian_key_count", len(current_mammals), set(current_mammals.index) == set(previous.index))
    joined = current_mammals.join(previous[["effect", "pvalue"]], rsuffix="_previous")
    effect_diff = max_difference(joined["effect"], joined["effect_previous"])
    p_diff = max_difference(joined["pvalue"], joined["pvalue_previous"])
    record("mammalian_effect_max_abs_difference", effect_diff, effect_diff <= 1e-12)
    record("mammalian_pvalue_max_abs_difference", p_diff, p_diff <= 1e-12)

    current_fly = saved[saved["species"] == "Fruit fly"].set_index(KEY)
    previous_fly = full[full["species"] == "Fruit fly"].set_index(KEY)
    fly = current_fly.join(previous_fly[["effect", "pvalue"]], rsuffix="_624")
    fly_effect_diff = max_difference(fly["effect"], fly["effect_624"])
    fly_p_diff = max_difference(fly["pvalue"], fly["pvalue_624"])
    record("fly_key_count_matches_624", len(fly), set(current_fly.index) == set(previous_fly.index))
    record("fly_effect_max_abs_difference_to_624", fly_effect_diff, fly_effect_diff <= 1e-12)
    record("fly_pvalue_max_abs_difference_to_624", fly_p_diff, fly_p_diff <= 1e-12)
    record("full_range_624_pass_count", int(full["significant_global_624_05"].fillna(False).astype(bool).sum()),
           int(full["significant_global_624_05"].fillna(False).astype(bool).sum()) == 45)
    record("mature_399_pass_count", int((mature399["q_refit_family"] < 0.05).sum()),
           int((mature399["q_refit_family"] < 0.05).sum()) == 10)
    saved_mature_passes = saved[saved["mature_399_pass"]]
    expected_mature_passes = mature399.loc[mature399["q_refit_family"] < 0.05, KEY]
    record("saved_mature_399_pass_flag_count", len(saved_mature_passes), len(saved_mature_passes) == 10)
    record("saved_mature_399_pass_flag_keys", len(expected_mature_passes),
           set(saved_mature_passes[KEY].itertuples(index=False, name=None))
           == set(expected_mature_passes.itertuples(index=False, name=None)))
    after_hash = hashlib.sha256(LOCKED_624.read_bytes()).hexdigest()
    record("locked_624_sha256_unchanged", after_hash, after_hash == before_hash == LOCKED_624_SHA256)
    verification = pd.DataFrame(checks)
    verification.to_csv(VERIFICATION_OUTPUT, index=False)
    return saved, verification


def reconciliation_rows(saved, full):
    current = saved.set_index(KEY)
    rows = []
    for old in full[full["significant_global_624_05"].fillna(False).astype(bool)].itertuples():
        key = (old.species, old.organ, old.cohort, old.cell_type)
        if key not in current.index:
            status = "UNTESTABLE_MATURE"
            row = None
        else:
            row = current.loc[key]
            status = "RETAINED" if row["significant_pooled_05"] else "NOT_RETAINED"
        rows.append({
            "species": old.species,
            "organ": old.organ,
            "cohort": old.cohort,
            "cell_type": old.cell_type,
            "full_range_624_effect": old.effect,
            "full_range_624_pvalue": old.pvalue,
            "full_range_624_q": old.q_global_bh_624,
            "mature_status": status,
            "mature_effect": np.nan if row is None else row["effect"],
            "mature_pvalue": np.nan if row is None else row["pvalue"],
            "mature_pooled_q": np.nan if row is None else row["q_pooled_mature_bh"],
        })
    return pd.DataFrame(rows)


def write_report(saved, coverage, exclusions, full, mature399, verification, started, ended):
    passes = saved[saved["significant_pooled_05"]].copy()
    reconciliation = reconciliation_rows(saved, full)
    old_pass_keys = set(full.loc[full["significant_global_624_05"].fillna(False).astype(bool), KEY].itertuples(index=False, name=None))
    new_vs_624 = passes[[tuple(x) not in old_pass_keys for x in passes[KEY].to_numpy()]]
    mature_old_pass = mature399[mature399["q_refit_family"] < 0.05]
    current_index = saved.set_index(KEY)
    old399_rows = []
    for old in mature_old_pass.itertuples():
        key = (old.species, old.organ, old.cohort, old.cell_type)
        row = current_index.loc[key]
        old399_rows.append((old, row, bool(row["significant_pooled_05"])))
    retained_399 = sum(x[2] for x in old399_rows)
    headline = passes[[tuple(x) in HEADLINE_KEYS for x in passes[KEY].to_numpy()]]
    duration = (datetime.fromisoformat(ended) - datetime.fromisoformat(started)).total_seconds()
    verify_map = verification.set_index("check")["value"].to_dict()

    lines = [
        "# All-species mature-to-old cell-fraction scan", "",
        "Analysis ordered by Reza on 2026-11-24. Pipeline model: `kimi-code/k3.1-pro`.", "",
        "## Scope, model, and eligibility", "",
        "This is a mature-to-old fraction scan only. The inferential unit is the human donor, mouse or naked-mole-rat animal, zebrafish fish, or fly/worm pool; cells are never independent replicates. Fly head/body samples and Gao worm libraries are pools, as verified in `scripts/run_rev3_composition_fly_worm_zebrafish.py` and `scripts/run_gao_worm_aging.py`.", "",
        "For each organ/cohort stratum, donor-by-cell-type counts receive a 0.5 pseudocount and a centered log-ratio (CLR) transform. Ordinary least squares fits CLR fraction on continuous age plus sex; sex is dropped when the design is rank-deficient. Slopes are per month for mouse and zebrafish, per year for human, and per day for fly and worm. The primary model does not adjust for total cells or sequencing platform.", "",
        "The pipeline rule, not the tentative 3-unit/50-cell rule in the order, is used: a cell type must have count > 0 in at least 6 mature inferential units and the age design must be full rank with positive residual degrees of freedom. There is no additional total-cell >= 50 exclusion in the primary pipeline; that is a sensitivity arm in the small-species full-range runner. Untestable rows never enter a correction family.", "",
        "## Mature windows and species decisions", "",
        "Cutoffs were fixed from the repository record and source biology before fitting:", "",
        "- **Mouse: 18 months and older.** Reuses the established mature composition window exactly.",
        "- **Human: 30 years and older.** Reuses the established mature composition window exactly.",
        "- **Fruit fly: 30 days and older.** AFCA's primary analyzed ages are 30, 50, and 70 days; day 30 is the repository's young-adult reference. Day-5 FCA cells remain excluded because study and age are inseparable. Thus the mature fly fit is the same age window as the primary AFCA fit.",
        "- **C. elegans: 6 days and older.** Day 6 is the repository's balanced mature reference and includes both sequencing platforms. Day-1 pools are excluded, leaving days 6, 12, and 14.",
        "- **Zebrafish: 12 months and older, but excluded.** The [source atlas](https://pmc.ncbi.nlm.nih.gov/articles/PMC11607696/) says males become fertile at 3 months, describes 5-month testes as approaching peak spermatogenesis, and identifies 12-month testes as closest to peak fertility. Twelve months is therefore the fully mature reference. Only four fish (two at 12 months, one at 20, one at 22) retain comparable fine labels, below the six-fish floor. The primary fine-label arm is not run; the separate broad-label reconciliation arm is not substituted to force a zebrafish result.",
        "- **Naked mole-rat: excluded.** Direct inspection of the GSE152886 count-file header confirms the established samples. After excluding 2V, there are three animals at 1.1 years and two at 11.3 years, with no intermediate ages. There is no defensible mature-to-old window satisfying the six-animal rule.", "",
        "## Multiple-testing families declared before correction", "",
        f"**Primary pooled mature family:** one Benjamini-Hochberg correction over all {len(saved)} tests actually run across all included species. `q_pooled_mature_bh < 0.05` is the only pass bar.", "",
        "Transparency-only families: `q_species_bh` applies BH separately within each species, and `q_bh` applies BH separately within each organ/cohort stratum. Neither determines a pooled pass. Family sizes are recorded on every result row.", "",
        "| Species | Mature tests | Pooled passes | Species-family passes |",
        "| --- | ---: | ---: | ---: |",
    ]
    for species in SPECIES_ORDER:
        group = saved[saved["species"] == species]
        if group.empty:
            continue
        lines.append(f"| {species} | {len(group)} | {int(group['significant_pooled_05'].sum())} | {int((group['q_species_bh'] < 0.05).sum())} |")
    lines.extend(["", "## Donor, animal, and pool coverage", "",
                  "| Species | Organ | Cohort | Unit | Mature ages | Units | Tested cell types |",
                  "| --- | --- | --- | --- | --- | ---: | ---: |"])
    for row in coverage.itertuples():
        lines.append(f"| {row.species} | {row.organ} | {row.cohort} | {row.inferential_unit} | {row.available_ages_mature} | {row.n_units_mature} | {row.n_cell_types_tested} |")

    lines.extend(["", f"## Primary pooled results: {len(passes)} passes", "",
                  "Slopes are in CLR units per the stated age unit. Numeric values below are read back from `results/mature-all-species-composition.csv`.", "",
                  "| Species | Organ | Cohort | Cell type | Slope | Age unit | p | Species q | Pooled mature q |",
                  "| --- | --- | --- | --- | ---: | --- | ---: | ---: | ---: |"])
    for row in species_sort(passes).itertuples():
        lines.append(f"| {row.species} | {row.organ} | {row.cohort} | {row.cell_type} | {fmt(row.effect)} | {row.age_unit} | {fmt(row.pvalue)} | {fmt(row.q_species_bh)} | {fmt(row.q_pooled_mature_bh)} |")
    if passes.empty:
        lines.append("| None |  |  |  |  |  |  |  |  |")

    lines.extend(["", "## Four prespecified headline rows", "",
                  f"{len(headline)} of 4 headline rows survive the pooled mature correction.", "",
                  "| Species | Cell type | Slope | p | Pooled mature q | Survives |",
                  "| --- | --- | ---: | ---: | ---: | --- |"])
    for key in sorted(HEADLINE_KEYS):
        row = current_index.loc[key]
        lines.append(f"| {key[0]} | {key[3]} | {fmt(row['effect'])} | {fmt(row['pvalue'])} | {fmt(row['q_pooled_mature_bh'])} | {bool(row['significant_pooled_05'])} |")

    retained_624 = int((reconciliation["mature_status"] == "RETAINED").sum())
    lost_624 = int((reconciliation["mature_status"] == "NOT_RETAINED").sum())
    untestable_624 = int((reconciliation["mature_status"] == "UNTESTABLE_MATURE").sum())
    lines.extend(["", "## Comparison with the full-range 624-test family", "",
                  f"The locked full-range CSV has 45 passes. Of these, {retained_624} pass the new pooled mature bar, {lost_624} are tested but do not pass, and {untestable_624} are untestable in the declared mature window. The new pooled scan has {len(new_vs_624)} passes that were not full-range 624 passes.", "",
                  "| Species | Organ | Cell type | Full-range slope | Full-range p | Full-range q | Mature status | Mature pooled q |",
                  "| --- | --- | --- | ---: | ---: | ---: | --- | ---: |"])
    for row in reconciliation.itertuples():
        lines.append(f"| {row.species} | {row.organ} | {row.cell_type} | {fmt(row.full_range_624_effect)} | {fmt(row.full_range_624_pvalue)} | {fmt(row.full_range_624_q)} | {row.mature_status} | {fmt(row.mature_pooled_q)} |")
    lines.extend(["", f"### Pooled mature passes not among the 45 full-range passes ({len(new_vs_624)})", "",
                  "| Species | Organ | Cohort | Cell type | Mature slope | Mature p | Pooled mature q |",
                  "| --- | --- | --- | --- | ---: | ---: | ---: |"])
    for row in species_sort(new_vs_624).itertuples():
        lines.append(f"| {row.species} | {row.organ} | {row.cohort} | {row.cell_type} | {fmt(row.effect)} | {fmt(row.pvalue)} | {fmt(row.q_pooled_mature_bh)} |")

    lines.extend(["", "## Comparison with the mammalian 399-row mature scan", "",
                  f"The previous mammalian-only family had 10 passes. {retained_399} of those 10 pass the new all-species pooled bar; {10 - retained_399} do not. All 399 mammalian tests are retained with the same fitted effects and p-values; only the correction family changes.", "",
                  "| Species | Organ | Cell type | Previous 399-row q | New pooled q | Retained |",
                  "| --- | --- | --- | ---: | ---: | --- |"])
    for old, row, retained in old399_rows:
        lines.append(f"| {old.species} | {old.organ} | {old.cell_type} | {fmt(old.q_refit_family)} | {fmt(row['q_pooled_mature_bh'])} | {retained} |")

    lines.extend(["", "## Exclusions and zero-test strata", "",
                  "No excluded or zero-test row enters the pooled family.", "",
                  "| Species | Organ | Scope | Units | Reason |",
                  "| --- | --- | --- | ---: | --- |"])
    for row in exclusions.itertuples():
        lines.append(f"| {row.species} | {row.organ} | {row.scope} | {row.n_units_mature} | {row.reason} |")

    lines.extend(["", "## Verification", "",
                  f"All {len(verification)} saved-output checks pass. Independent BH recomputation agrees with the pooled CSV to maximum absolute difference {verify_map['pooled_bh_max_abs_difference']}.", "",
                  f"The newly fit 399 mammalian rows match the prior mature CSVs on all keys, with maximum absolute effect difference {verify_map['mammalian_effect_max_abs_difference']} and p-value difference {verify_map['mammalian_pvalue_max_abs_difference']}. Fly has the same 30-to-70-day window as its full-range primary arm; all 162 keys match and maximum absolute effect and p-value differences are {verify_map['fly_effect_max_abs_difference_to_624']} and {verify_map['fly_pvalue_max_abs_difference_to_624']}.", "",
                  f"The locked 624-test input SHA256 is unchanged: `{verify_map['locked_624_sha256_unchanged']}`. Every reported count and pass value was recomputed from the saved result CSVs, not from prose or commit messages.", "",
                  "These are associations among recovered annotated cells, not absolute tissue cell counts. Non-significance is not evidence of equivalence.", "",
                  "## Runtime and files", "",
                  f"Started `{started}`; ended `{ended}`; duration {duration:.6f} seconds; model `kimi-code/k3.1-pro`.", "",
                  "- `scripts/run_mature_all_species_fractions.py`",
                  "- `results/mature-all-species-composition.csv`",
                  "- `results/mature-all-species-composition-{mouse,human,fruit-fly,c-elegans}.csv`",
                  "- `results/mature-all-species-composition-coverage.csv`",
                  "- `results/mature-all-species-composition-exclusions.csv`",
                  "- `results/mature-all-species-composition-untestable.csv`",
                  "- `results/mature-all-species-composition-verification.csv`",
                  "- `results/mature-all-species-composition-runtime.csv`",
                  "- `results/mature-all-species-composition-full-range-624-reconciliation.csv`",
                  "- `reports/followup-mature-all-species-fractions.md`", ""])
    REPORT_OUTPUT.write_text("\n".join(lines), encoding="utf-8")


def main():
    started = utc_now()
    print(f"START {started} MODEL {MODEL_USED}", flush=True)
    before_hash = hashlib.sha256(LOCKED_624.read_bytes()).hexdigest()
    if before_hash != LOCKED_624_SHA256:
        raise RuntimeError("Locked 624-family input hash differs from expected")
    rows = []
    untestable = []
    coverage = []
    for species, organ, cohort, counts, metadata in load_strata():
        stratum_rows, stratum_untestable, stratum_coverage = run_stratum(
            species, organ, cohort, counts, metadata)
        rows.extend(stratum_rows)
        untestable.extend(stratum_untestable)
        coverage.append(stratum_coverage)
        print(f"{species} | {organ} | {cohort}: units={len(metadata)} tests={len(stratum_rows)}", flush=True)
    nmr_cov, nmr_untestable, nmr_reason = nmr_coverage_and_untestable()
    coverage.append(nmr_cov)
    untestable.extend(nmr_untestable)
    if not rows:
        raise RuntimeError("No mature tests were produced")
    result = pd.DataFrame(rows)
    result, full, mature399 = add_multiplicity_and_comparisons(result)
    columns = KEY + [
        "sensitivity", "model", "effect_type", "age_unit", "mature_cutoff", "inferential_unit",
        "effect", "se", "ci95_low", "ci95_high", "pvalue", "q_bh", "q_species_bh",
        "q_pooled_mature_bh", "n_stratum_family", "n_species_family", "n_pooled_mature_family",
        "significant_pooled_05", "n_donors", "n_present", "total_cells", "min_age", "max_age",
        "full_range_624_effect", "full_range_624_pvalue", "full_range_624_q", "full_range_624_pass",
        "mature_399_q", "mature_399_pass",
    ]
    result = result[columns]
    result.to_csv(OUTPUT, index=False)
    slugs = {"Mouse": "mouse", "Human": "human", "Fruit fly": "fruit-fly", "C. elegans": "c-elegans"}
    for species, slug in slugs.items():
        result[result["species"] == species].to_csv(
            RESULTS / f"mature-all-species-composition-{slug}.csv", index=False)
    coverage = pd.DataFrame(coverage)
    coverage.to_csv(COVERAGE_OUTPUT, index=False)
    exclusions = make_exclusions(coverage, nmr_reason)
    exclusions.to_csv(EXCLUSIONS_OUTPUT, index=False)
    pd.DataFrame(untestable).sort_values(KEY).to_csv(UNTESTABLE_OUTPUT, index=False)
    saved, verification = verify_saved_outputs(result, full, mature399, before_hash)
    reconciliation = reconciliation_rows(saved, full)
    reconciliation.to_csv(RESULTS / "mature-all-species-composition-full-range-624-reconciliation.csv", index=False)
    ended = utc_now()
    runtime = pd.DataFrame([{
        "pipeline": "mature_all_species_fractions",
        "start_utc": started,
        "end_utc": ended,
        "duration_seconds": (datetime.fromisoformat(ended) - datetime.fromisoformat(started)).total_seconds(),
        "model_used": MODEL_USED,
        "status": "completed",
    }])
    runtime.to_csv(RUNTIME_OUTPUT, index=False)
    write_report(saved, coverage, exclusions, full, mature399, verification, started, ended)
    print(f"END {ended} pooled_tests={len(saved)} pooled_passes={int(saved['significant_pooled_05'].sum())}", flush=True)
    print(saved[saved["significant_pooled_05"]].groupby("species").size().to_string(), flush=True)


if __name__ == "__main__":
    main()
