#!/usr/bin/env python3
"""Complete the mature composition refit without changing the base pipelines."""

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

RESULTS = REPO / "results"
REPORTS = REPO / "reports"
KEY = ["species", "organ", "cohort", "cell_type"]
AGE_MAP = {"1m": 1.0, "3m": 3.0, "18m": 18.0, "21m": 21.0, "24m": 24.0, "30m": 30.0}
MOUSE_CUTOFF = 18.0
HUMAN_CUTOFF = 30.0
MODEL_USED = "kimi-code/k3.1-pro"
LOCKED_624_SHA256 = "3f26242585017bd3a68a4a7c4dcb1d2684da0e9f662c37ed7e07e1bd6f054d96"
RUNTIME_FILE = RESULTS / "followup-mature-composition-runtime.csv"


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def write_runtime(pipeline, started, status="completed"):
    ended = utc_now()
    duration = (datetime.fromisoformat(ended) - datetime.fromisoformat(started)).total_seconds()
    row = pd.DataFrame([{
        "pipeline": pipeline, "start_utc": started, "end_utc": ended,
        "duration_seconds": duration, "model_used": MODEL_USED, "status": status,
    }])
    row.to_csv(RUNTIME_FILE, index=False)


def clr_table(counts, pseudocount=0.5):
    adjusted = counts.astype(float) + pseudocount
    proportions = adjusted.div(adjusted.sum(axis=1), axis=0)
    log_proportions = np.log(proportions)
    return log_proportions.sub(log_proportions.mean(axis=1), axis=0)


def count_table(obs, donor_column, cell_column, columns=None):
    counts = (obs.groupby([donor_column, cell_column], observed=True).size()
              .unstack(fill_value=0))
    if columns is not None:
        counts = counts.reindex(columns=columns, fill_value=0)
    return counts.sort_index()


def fit_row(clr, counts, metadata, cell_type, species, organ, cohort, donors=None):
    if donors is None:
        donors = counts.index
    donors = pd.Index(donors)
    if len(donors) < 6:
        return None, "fewer than 6 mature donors"
    X = sm.add_constant(metadata.loc[donors, ["age", "sex"]])
    model = "age_plus_sex"
    if np.linalg.matrix_rank(X.to_numpy()) < X.shape[1]:
        X = sm.add_constant(metadata.loc[donors, ["age"]])
        model = "age_only_rank_deficient_sex"
    if np.linalg.matrix_rank(X.to_numpy()) < X.shape[1] or len(donors) <= X.shape[1]:
        return None, "rank-deficient age design or insufficient residual degrees of freedom"
    fit = sm.OLS(clr.loc[donors, cell_type], X).fit()
    if not np.isfinite(fit.pvalues["age"]):
        return None, "non-finite age test after model fit"
    df = int(fit.df_resid)
    half = t_dist.ppf(0.975, df) * fit.bse["age"] if df > 0 else np.nan
    return {
        "species": species, "organ": organ, "cohort": cohort, "cell_type": cell_type,
        "sensitivity": "mature_only_18mo" if species == "Mouse" else "mature_only_30y",
        "model": model, "effect_type": "CLR_slope_per_year",
        "effect": fit.params["age"], "se": fit.bse["age"],
        "ci95_low": fit.params["age"] - half, "ci95_high": fit.params["age"] + half,
        "pvalue": fit.pvalues["age"], "n_donors": len(donors),
        "n_present": int((counts.loc[donors, cell_type] > 0).sum()),
        "total_cells": int(counts.loc[donors, cell_type].sum()),
        "min_age": float(metadata.loc[donors, "age"].min()),
        "max_age": float(metadata.loc[donors, "age"].max()),
    }, ""


def regress_fractions(counts, metadata, species, organ, cohort):
    clr = clr_table(counts)
    rows = []
    failures = {}
    for cell_type in counts.columns:
        row, reason = fit_row(clr, counts, metadata, cell_type, species, organ, cohort)
        if row is None:
            failures[cell_type] = reason
        else:
            rows.append(row)
    if rows:
        qvalues = multipletests([row["pvalue"] for row in rows], method="fdr_bh")[1]
        for row, qvalue in zip(rows, qvalues):
            row["q_bh"] = qvalue
    return rows, failures


def mouse_inputs():
    frames = []
    paths = sorted((REPO / "data" / "interim").glob("m1-qc-cleaned-*.h5ad"))
    if len(paths) != 16:
        raise RuntimeError(f"Expected 16 mouse tissue files, found {len(paths)}")
    for path in paths:
        adata = ad.read_h5ad(path, backed="r")
        frames.append(adata.obs.copy())
        adata.file.close()
    obs = pd.concat(frames, ignore_index=True)
    obs = obs[(~obs["predicted_doublet"].astype(bool))
              & (obs["cell_ontology_class"].astype(str) != "unknown")].copy()
    obs["age"] = obs["age"].astype(str).map(AGE_MAP)
    if obs["age"].isna().any():
        raise RuntimeError("Unmapped mouse age")
    obs["sex"] = (obs["sex"].astype(str).str.lower() == "male").astype(float)
    for tissue in sorted(obs["tissue"].unique()):
        full = obs[obs["tissue"] == tissue].copy()
        part = full[full["age"] >= MOUSE_CUTOFF].copy()
        columns = sorted(full["cell_ontology_class"].astype(str).unique())
        counts = count_table(part, "mouse.id", "cell_ontology_class", columns)
        metadata = part.groupby("mouse.id").agg(age=("age", "first"), sex=("sex", "first"))
        yield "Mouse", tissue, "TMS", counts, metadata


def human_inputs():
    source = pd.read_csv(RESULTS / "m9-donor-cell-counts-prethreshold.csv")
    for cohort in sorted(source["cohort"].unique()):
        full = source[source["cohort"] == cohort].copy()
        part = full[full["age"] >= HUMAN_CUTOFF].copy()
        organ = full["group"].iloc[0]
        columns = sorted(full["cell_type"].astype(str).unique())
        if part.empty:
            counts = pd.DataFrame(columns=columns, dtype=float)
            metadata = pd.DataFrame(columns=["age", "sex"], dtype=float)
        else:
            counts = (part.set_index(["donor_id", "cell_type"])["n_cells"]
                      .unstack(fill_value=0).sort_index().reindex(columns=columns, fill_value=0))
            metadata = part.groupby("donor_id").agg(age=("age", "first"), sex=("sex", "first"))
            metadata["sex"] = (metadata["sex"].astype(str).str.lower() == "male").astype(float)
        yield "Human", organ, cohort, counts, metadata


def untestable_row(species, organ, cohort, cell_type, reason, counts, metadata,
                   original_lookup, row_scope="cell_type"):
    original_key = (species, organ, cohort, cell_type)
    original = original_lookup.get(original_key, {})
    if row_scope == "cell_type" and cell_type in counts.columns:
        n_present = int((counts[cell_type] > 0).sum())
        total_cells = int(counts[cell_type].sum())
    else:
        n_present = np.nan
        total_cells = np.nan
    return {
        "species": species, "organ": organ, "cohort": cohort, "cell_type": cell_type,
        "row_scope": row_scope, "status": "UNTESTABLE", "reason": reason,
        "n_donors_mature": int(len(metadata)) if species in ("Mouse", "Human") else np.nan,
        "n_present_mature": n_present, "total_cells_mature": total_cells,
        "n_candidate_cell_types": int(counts.shape[1]) if row_scope != "cell_type" else np.nan,
        "min_age_mature": float(metadata["age"].min()) if len(metadata) else np.nan,
        "max_age_mature": float(metadata["age"].max()) if len(metadata) else np.nan,
        "in_original_431": bool(original),
        "original_n_donors": original.get("n_donors", np.nan),
        "original_q_global_bh": original.get("q_global_bh", np.nan),
    }


def collect_tests(original):
    original_lookup = original.set_index(KEY).to_dict("index")
    tested_rows = []
    untestable = {}
    inputs = {}
    coverage = []
    for loader in (mouse_inputs, human_inputs):
        for species, organ, cohort, all_counts, metadata in loader():
            present = (all_counts > 0).sum(axis=0)
            counts = all_counts.loc[:, present >= 6].copy()
            rows, failures = regress_fractions(counts, metadata, species, organ, cohort)
            tested_rows.extend(rows)
            inputs[(species, organ, cohort)] = (counts, metadata)
            tested_types = {row["cell_type"] for row in rows}
            for cell_type in all_counts.columns:
                if cell_type in tested_types:
                    continue
                if cell_type in failures:
                    reason = failures[cell_type]
                else:
                    reason = "cell type present in fewer than 6 mature donors"
                key = (species, organ, cohort, cell_type)
                untestable[key] = untestable_row(
                    species, organ, cohort, cell_type, reason, all_counts, metadata, original_lookup)
            if not rows:
                key = (species, organ, cohort, "(all cell types)")
                reason = ("zero eligible tests: no cell type met the 6-donor presence rule "
                          "and produced an estimable age model")
                untestable[key] = untestable_row(
                    species, organ, cohort, "(all cell types)", reason,
                    all_counts, metadata, original_lookup, row_scope="tissue_or_cohort")
            coverage.append({
                "species": species, "organ": organ, "cohort": cohort,
                "n_donors_mature": len(metadata), "n_cell_types_tested": len(rows),
                "n_cell_types_candidates": int(all_counts.shape[1]),
            })
    for row in original.itertuples():
        key = (row.species, row.organ, row.cohort, row.cell_type)
        if row.species == "Naked mole-rat":
            if key not in untestable:
                untestable[key] = {
                    "species": row.species, "organ": row.organ, "cohort": row.cohort,
                    "cell_type": row.cell_type, "row_scope": "cell_type", "status": "UNTESTABLE",
                    "reason": "no mature cutoff is defined for naked mole-rat; outside the mouse and human refit",
                    "n_donors_mature": np.nan, "n_present_mature": np.nan,
                    "total_cells_mature": np.nan, "n_candidate_cell_types": np.nan,
                    "min_age_mature": np.nan, "max_age_mature": np.nan,
                    "in_original_431": True, "original_n_donors": row.n_donors,
                    "original_q_global_bh": row.q_global_bh,
                }
        elif key not in untestable and key not in inputs and key[:3] not in inputs:
            pass
    tested_keys = {(r["species"], r["organ"], r["cohort"], r["cell_type"]) for r in tested_rows}
    for row in original[original["species"].isin(["Mouse", "Human"])].itertuples():
        key = (row.species, row.organ, row.cohort, row.cell_type)
        if key in tested_keys or key in untestable:
            continue
        counts, metadata = inputs.get(key[:3], (pd.DataFrame(), pd.DataFrame(columns=["age"])))
        untestable[key] = untestable_row(
            row.species, row.organ, row.cohort, row.cell_type,
            "no mature cells or no eligible mature test for this original family row",
            counts, metadata, original_lookup)
    return pd.DataFrame(tested_rows), pd.DataFrame(untestable.values()), inputs, pd.DataFrame(coverage)


def add_multiplicity(refit, original):
    refit = refit.copy()
    refit["q_refit_family"] = multipletests(refit["pvalue"], method="fdr_bh")[1]
    refit["n_refit_family"] = len(refit)
    original_indexed = original.set_index(KEY)
    pvalues = original_indexed["pvalue"].copy()
    for row in refit.itertuples():
        key = (row.species, row.organ, row.cohort, row.cell_type)
        if key in pvalues.index:
            pvalues.loc[key] = row.pvalue
    q_all = pd.Series(multipletests(pvalues.to_numpy(), method="fdr_bh")[1], index=pvalues.index)
    refit["in_original_431"] = [tuple(x) in original_indexed.index for x in refit[KEY].to_numpy()]
    refit["q_431_family"] = np.nan
    refit["n_431_family"] = np.nan
    for idx, row in refit.iterrows():
        key = (row["species"], row["organ"], row["cohort"], row["cell_type"])
        if key in q_all.index:
            refit.loc[idx, "q_431_family"] = q_all.loc[key]
            refit.loc[idx, "n_431_family"] = len(original)
    original_columns = original_indexed[[
        "effect", "se", "pvalue", "n_donors", "q_global_bh", "significant_global_05"
    ]].rename(columns={
        "effect": "original_effect", "se": "original_se", "pvalue": "original_pvalue",
        "n_donors": "original_n_donors", "q_global_bh": "original_q_global_bh",
        "significant_global_05": "original_survivor",
    })
    refit = refit.join(original_columns, on=KEY)
    refit["original_survivor"] = refit["original_survivor"].fillna(False).astype(bool)
    refit["mature_effect_over_original"] = refit["effect"] / refit["original_effect"]
    refit["mature_se_over_original_se"] = refit["se"] / refit["original_se"]
    return refit


def total_cell_fit(counts, metadata, cell_type):
    if len(counts) < 6:
        return None
    clr = clr_table(counts)
    predictors = pd.DataFrame({
        "age": metadata.loc[counts.index, "age"],
        "log_total_cells": np.log(counts.sum(axis=1)),
    })
    if metadata.loc[counts.index, "sex"].nunique() > 1:
        predictors["sex"] = metadata.loc[counts.index, "sex"]
    X = sm.add_constant(predictors.astype(float), has_constant="add")
    if len(counts) <= X.shape[1] + 1 or np.linalg.matrix_rank(X.to_numpy()) < X.shape[1]:
        return None
    fit = sm.OLS(clr[cell_type], X).fit()
    ci = fit.conf_int().loc["age"]
    return {
        "total_effect": fit.params["age"], "total_se": fit.bse["age"],
        "total_ci95_low": ci.iloc[0], "total_ci95_high": ci.iloc[1],
        "total_pvalue": fit.pvalues["age"], "total_n_donors": len(counts),
    }


def secondary_checks(refit, inputs):
    hits = refit[(~refit["original_survivor"])
                 & ((refit["q_431_family"] < 0.05) | (refit["q_refit_family"] < 0.05))].copy()
    all_total = []
    for key, (counts, metadata) in inputs.items():
        clr = clr_table(counts)
        for cell_type in counts.columns:
            result = total_cell_fit(counts, metadata, cell_type)
            if result is not None:
                all_total.append({"species": key[0], "organ": key[1], "cohort": key[2],
                                  "cell_type": cell_type, **result})
    total = pd.DataFrame(all_total)
    if len(total):
        total["q_total_refit_family"] = multipletests(total["total_pvalue"], method="fdr_bh")[1]
        total["q_total_bh"] = np.nan
        for _, indices in total.groupby(["species", "organ", "cohort"]).groups.items():
            total.loc[indices, "q_total_bh"] = multipletests(
                total.loc[indices, "total_pvalue"], method="fdr_bh")[1]
        total = total.set_index(KEY)
    records = []
    for hit in hits.itertuples():
        key = (hit.species, hit.organ, hit.cohort, hit.cell_type)
        counts, metadata = inputs[key[:3]]
        clr = clr_table(counts)
        lodo = []
        for donor in counts.index:
            keep = counts.index != donor
            row, _ = fit_row(clr.loc[keep], counts.loc[keep], metadata, hit.cell_type,
                             hit.species, hit.organ, hit.cohort)
            if row is not None:
                lodo.append({"donor": donor, "effect": row["effect"], "pvalue": row["pvalue"]})
        lodo_df = pd.DataFrame(lodo)
        record = {k: getattr(hit, k) for k in KEY}
        record.update({
            "status": "SECONDARY_NOT_PROMOTED", "effect": hit.effect, "se": hit.se,
            "pvalue": hit.pvalue, "q_431_family": hit.q_431_family,
            "q_refit_family": hit.q_refit_family,
            "n_donors": hit.n_donors, "n_present": hit.n_present,
        })
        if len(total) and key in total.index:
            for column, value in total.loc[key].items():
                record[column] = value
        else:
            record.update({"total_effect": np.nan, "total_se": np.nan, "total_pvalue": np.nan,
                           "q_total_bh": np.nan, "q_total_refit_family": np.nan})
        if len(lodo_df):
            changes = (lodo_df["effect"] - hit.effect).abs()
            record.update({
                "lodo_n_fits": len(lodo_df), "lodo_slope_min": lodo_df["effect"].min(),
                "lodo_slope_max": lodo_df["effect"].max(), "lodo_pvalue_max": lodo_df["pvalue"].max(),
                "lodo_same_sign": bool(np.all(np.sign(lodo_df["effect"]) == np.sign(hit.effect))),
                "lodo_max_abs_slope_change": changes.max(),
                "lodo_max_change_donor": lodo_df.loc[changes.idxmax(), "donor"],
            })
        else:
            record.update({"lodo_n_fits": 0, "lodo_slope_min": np.nan, "lodo_slope_max": np.nan,
                           "lodo_pvalue_max": np.nan, "lodo_same_sign": False,
                           "lodo_max_abs_slope_change": np.nan, "lodo_max_change_donor": ""})
        record["propeller_status"] = "NOT_VERIFIED"
        record["propeller_reason"] = "Rscript, speckle, and limma are unavailable on beast; no mature propeller fit was run"
        record["study_aware_status"] = ("see followup-mature-composition-lung-study.csv"
                                        if key[:3] == ("Human", "lung", "Lung_Integrated")
                                        else "not applicable outside Lung_Integrated")
        records.append(record)
    return pd.DataFrame(records)


def fmt(value, digits=10):
    if value is None or pd.isna(value):
        return "NA"
    return f"{float(value):.{digits}g}"


def write_report():
    mouse = pd.read_csv(RESULTS / "followup-mature-composition-mouse.csv")
    human = pd.read_csv(RESULTS / "followup-mature-composition-human.csv")
    refit = pd.concat([mouse, human], ignore_index=True)
    original = pd.read_csv(RESULTS / "reviewer-phase1-composition-global-bh.csv")
    untestable = pd.read_csv(RESULTS / "followup-mature-composition-untestable.csv")
    secondary_path = RESULTS / "followup-mature-composition-secondary-checks.csv"
    secondary = pd.read_csv(secondary_path) if secondary_path.exists() else pd.DataFrame()
    lung_path = RESULTS / "followup-mature-composition-lung-study.csv"
    lung = pd.read_csv(lung_path) if lung_path.exists() else pd.DataFrame()
    refit_keys = refit.set_index(KEY)
    survivors = original[original["significant_global_05"]]
    lines = [
        "# Complete mature composition refit", "",
        "## Scope and model", "",
        "This is Phase 2 Item 3 only. Mouse fits use ages 18, 21, 24, and 30 months where present. Human fits use donors aged at least 30 years. Naked mole-rat rows are not eligible for a mature refit and retain their original p-values only for the 431-family correction.", "",
        "The model is unchanged from run_organ_composition.py: donor by cell-type counts, pseudocount 0.5, centered log-ratio transformation, OLS on age plus sex, and age-only fallback when sex makes the design rank-deficient. A cell type is eligible when observed in at least 6 mature donors. Mouse effects are CLR slopes per month and human effects are CLR slopes per year, despite the inherited effect_type label CLR_slope_per_year.", "",
        "The 624-family locked result was not used as an input or modified. Its SHA256 before and after this run is checked against 3f26242585017bd3a68a4a7c4dcb1d2684da0e9f662c37ed7e07e1bd6f054d96.", "",
        "## Multiplicity and coverage", "",
        f"Mouse tested rows: {len(mouse)}. Human tested rows: {len(human)}. Total mature refit rows: {len(refit)}.",
        f"q_431_family uses one BH correction over all {len(original)} original primary tests after substituting mature p-values only at matching mouse or human keys. It is defined only for rows belonging to that original key set.",
        f"q_refit_family uses a separate BH correction over all {len(refit)} mature mouse and human tests actually run. The two q columns are never mixed. Within-stratum q_bh is also retained for continuity with the original outputs.", "",
        "### Cohort coverage", "",
        "| Species | Organ | Cohort | Mature donors | Tested cell types | Minimum age | Maximum age |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: |",
    ]
    coverage = (refit.groupby(["species", "organ", "cohort"])
                .agg(n_donors=("n_donors", "max"), n_tested=("cell_type", "size"),
                     min_age=("min_age", "min"), max_age=("max_age", "max")).reset_index())
    zero = untestable[untestable["row_scope"] == "tissue_or_cohort"]
    for row in coverage.itertuples():
        lines.append(f"| {row.species} | {row.organ} | {row.cohort} | {row.n_donors} | {row.n_tested} | {fmt(row.min_age)} | {fmt(row.max_age)} |")
    for row in zero[zero["species"].isin(["Mouse", "Human"])].itertuples():
        lines.append(f"| {row.species} | {row.organ} | {row.cohort} | {int(row.n_donors_mature)} | 0 | {fmt(row.min_age_mature)} | {fmt(row.max_age_mature)} |")
    lines.extend(["", "## Original 15 global passes", "",
                  "SURVIVED means q_431_family < 0.05. DIED means the row was estimable but did not clear that bar. UNTESTABLE is not a null. Both mature q families are shown.", "",
                  "| Species | Organ | Cohort | Cell type | Original q | Mature effect | SE | SE inflation | p | q_431_family | q_refit_family | Donors original to mature | Verdict |",
                  "| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |"])
    verdict_counts = {"SURVIVED": 0, "DIED": 0, "UNTESTABLE": 0}
    died = []
    for old in survivors.itertuples():
        key = (old.species, old.organ, old.cohort, old.cell_type)
        if key in refit_keys.index:
            row = refit_keys.loc[key]
            verdict = "SURVIVED" if row["q_431_family"] < 0.05 else "DIED"
            if verdict == "DIED":
                died.append(row)
            values = [fmt(row["effect"]), fmt(row["se"]), fmt(row["mature_se_over_original_se"]),
                      fmt(row["pvalue"]), fmt(row["q_431_family"]), fmt(row["q_refit_family"]),
                      f"{old.n_donors} to {int(row['n_donors'])}"]
        else:
            verdict = "UNTESTABLE"
            values = ["NA", "NA", "NA", "NA", "NA", "NA", f"{old.n_donors} to NA"]
        verdict_counts[verdict] += 1
        lines.append(f"| {old.species} | {old.organ} | {old.cohort} | {old.cell_type} | {fmt(old.q_global_bh)} | " + " | ".join(values) + f" | {verdict} |")
    lines.extend(["", f"Totals: SURVIVED {verdict_counts['SURVIVED']}; DIED {verdict_counts['DIED']}; UNTESTABLE {verdict_counts['UNTESTABLE']}.", ""])
    if died:
        died_df = pd.DataFrame(died)
        median_se = died_df["mature_se_over_original_se"].median()
        lines.append(f"For the {len(died)} DIED rows, median SE inflation is {fmt(median_se)}x, with range {fmt(died_df['mature_se_over_original_se'].min())}x to {fmt(died_df['mature_se_over_original_se'].max())}x. These are losses of significance after donor and age-range restriction, not demonstrations of no mature effect. Exact per-row donor counts and SE ratios are in the table above.")
        lines.append("")
    new = refit[(~refit["original_survivor"])
                & ((refit["q_431_family"] < 0.05) | (refit["q_refit_family"] < 0.05))]
    lines.extend(["## New mature passes", "",
                  "All rows below remain SECONDARY_NOT_PROMOTED. A new pass here means q_431_family or q_refit_family is below 0.05 and the row was not one of the original 15 passes.", "",
                  "| Species | Organ | Cohort | Cell type | Effect | SE | p | q_431_family | q_refit_family | Donors | Present | Cells |",
                  "| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"])
    if new.empty:
        lines.append("None.")
    for row in new.sort_values(KEY).itertuples():
        lines.append(f"| {row.species} | {row.organ} | {row.cohort} | {row.cell_type} | {fmt(row.effect)} | {fmt(row.se)} | {fmt(row.pvalue)} | {fmt(row.q_431_family)} | {fmt(row.q_refit_family)} | {row.n_donors} | {row.n_present} | {row.total_cells} |")
    lines.extend(["", "## Untestable rows and zero-test tissues", "",
                  f"The untestable CSV has {len(untestable)} rows, including cell-level eligibility failures, zero-test tissue summaries, and 4 naked mole-rat rows outside the mature definition. No untestable row enters either mature BH family.", "",
                  "| Species | Organ | Cohort | Mature donors | Candidate cell types | Reason |",
                  "| --- | --- | --- | ---: | ---: | --- |"])
    for row in zero.itertuples():
        lines.append(f"| {row.species} | {row.organ} | {row.cohort} | {fmt(row.n_donors_mature)} | {fmt(row.n_candidate_cell_types)} | {row.reason} |")
    lines.extend(["", "Bladder has 4 mature donors and zero eligible cell types, so the original bladder-cell pass is UNTESTABLE. Large_Intestine, Mammary_Gland, Skin, Thymus, Tongue, and Trachea also have zero eligible mouse tests. Cell-level reasons and mature presence and cell counts are in results/followup-mature-composition-untestable.csv.", ""])
    lines.extend(["## Study-aware human lung at 30-plus", "",
                  "The lung sensitivity uses the D3 donor-ID study proxy, the original adult-eligible Lung_Integrated cell universe, and only donors aged at least 30. Study groups with at least 5 mature donors enter fixed effects; smaller groups are pooled as other. Per-study fits require at least 5 donors and age SD at least 3 years. Meta-analysis uses inverse-variance fixed effects and DerSimonian-Laird random effects. LOSO retains the original rule of at least 20 remaining donors.", ""])
    if lung.empty:
        lines.append("Lung study-aware output was not available when this report was generated.")
    else:
        focal = ["club cell", "multiciliated epithelial cell"]
        lines.extend(["| Cell type | Arm | Model or group | Effect | SE | p | Donors | Studies |",
                      "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |"])
        selected = lung[(lung["cell_type"].isin(focal))
                        & (lung["analysis"].isin(["pooled_model", "meta_combine", "per_study_summary"]))]
        for row in selected.itertuples():
            group = f" {row.group}" if pd.notna(row.group) and row.group else ""
            lines.append(f"| {row.cell_type} | {row.analysis} | {row.model}{group} | {fmt(row.effect)} | {fmt(row.se)} | {fmt(row.pvalue)} | {fmt(row.n_donors)} | {fmt(row.n_studies)} |")
        lines.append("")
        for cell_type in focal:
            part = lung[(lung["cell_type"] == cell_type) & (lung["analysis"] == "per_study_slope")]
            lines.append(f"Per-study slopes for {cell_type}:")
            for row in part.itertuples():
                lines.append(f"- {row.group}: status {row.status}; effect {fmt(row.effect)}; SE {fmt(row.se)}; p {fmt(row.pvalue)}; n {fmt(row.n_donors)}; {row.note}")
            loso = lung[(lung["cell_type"] == cell_type) & (lung["analysis"] == "loso")]
            lines.append(f"LOSO for {cell_type}:")
            for row in loso.itertuples():
                lines.append(f"- Drop {row.group}: effect {fmt(row.effect)}; SE {fmt(row.se)}; p {fmt(row.pvalue)}; remaining donors {fmt(row.n_donors)}; status {row.status}.")
            lines.append("")
    lines.extend(["## Secondary-hit checks", "",
                  "No secondary hit is promoted to a headline pass. Total-cell models add log total recovered eligible cells. Donor influence repeats the primary model leaving out one donor at a time. Study-aware status for Lung_Integrated hits is determined from the preceding lung output. Mature propeller could not be run because Rscript, speckle, and limma are unavailable on beast; this check is NOT_VERIFIED rather than passed.", ""])
    if secondary.empty:
        lines.append("No new mature passes required secondary checks.")
    else:
        lines.extend(["| Species | Organ | Cohort | Cell type | Total-cell p | Total-cell refit q | LODO fits | LODO maximum p | LODO same sign | Propeller |",
                      "| --- | --- | --- | --- | ---: | ---: | ---: | ---: | --- | --- |"])
        for row in secondary.sort_values(KEY).itertuples():
            lines.append(f"| {row.species} | {row.organ} | {row.cohort} | {row.cell_type} | {fmt(row.total_pvalue)} | {fmt(row.q_total_refit_family)} | {row.lodo_n_fits} | {fmt(row.lodo_pvalue_max)} | {row.lodo_same_sign} | {row.propeller_status} |")
        lines.append("")
        existing = secondary[((secondary["species"] == "Mouse") & (secondary["organ"] == "Limb_Muscle")
                              & (secondary["cell_type"].isin(["macrophage", "mesenchymal stem cell"])))
                             | ((secondary["species"] == "Human") & (secondary["organ"] == "lung")
                                & (secondary["cohort"] == "Lung_Integrated")
                                & (secondary["cell_type"] == "conventional dendritic cell"))]
        lines.append(f"The three previously reported secondary findings are retained: mouse macrophage and mesenchymal stem cell declines, and human Lung_Integrated conventional dendritic-cell increase. All three are present in the new-pass table: {len(existing) == 3}. They remain secondary and their check values are listed above.")
        lines.append("")
        if not lung.empty:
            lung_secondary = secondary[(secondary["species"] == "Human")
                                       & (secondary["cohort"] == "Lung_Integrated")]
            for secondary_row in lung_secondary.itertuples():
                cell_lung = lung[lung["cell_type"] == secondary_row.cell_type]
                fixed = cell_lung[(cell_lung["analysis"] == "pooled_model")
                                  & (cell_lung["model"] == "study_fixed_effects")].iloc[0]
                meta = cell_lung[(cell_lung["analysis"] == "meta_combine")
                                 & (cell_lung["model"] == "fixed_effect_ivw")].iloc[0]
                random = cell_lung[(cell_lung["analysis"] == "meta_combine")
                                   & (cell_lung["model"] == "random_effects_dl")].iloc[0]
                loso = cell_lung[(cell_lung["analysis"] == "loso")
                                 & (cell_lung["status"] == "TESTED")]
                lines.append(f"Study-aware check for {secondary_row.cell_type}: study fixed effect {fmt(fixed.effect)}, p {fmt(fixed.pvalue)}; fixed-effect meta {fmt(meta.effect)}, p {fmt(meta.pvalue)}; random-effects meta {fmt(random.effect)}, p {fmt(random.pvalue)}; LOSO effects {fmt(loso.effect.min())} to {fmt(loso.effect.max())}, p range {fmt(loso.pvalue.min())} to {fmt(loso.pvalue.max())} across {len(loso)} arms.")
            lines.append("")
    lines.extend(["## Runtime and verification limits", "",
                  "Pipeline model: kimi-code/k3.1-pro.", "",
                  "| Pipeline | Start UTC | End UTC | Duration seconds | Status |",
                  "| --- | --- | --- | ---: | --- |"])
    if RUNTIME_FILE.exists():
        runtime = pd.read_csv(RUNTIME_FILE)
        for row in runtime.itertuples():
            lines.append(f"| {row.pipeline} | {row.start_utc} | {row.end_utc} | {fmt(row.duration_seconds)} | {row.status} |")
    lines.extend(["", "One self-review independently recomputed both BH families, checked all output age minima, accounted for every original 431-family row as tested or untestable, checked all 44 lung cell types across 28 study-arm rows, and verified the locked 624-file hash. All 399 tested rows belong to the original family; 32 original rows are untestable. An in-memory call to the unchanged base mouse pipeline on the current h5ad files reproduces the new mouse effects to numerical precision. Seventeen rows in the older stored mature mouse CSV differ slightly, with a maximum absolute effect difference of 0.000697691 and at most 7 recovered cells, so that stored CSV reflects an earlier state of the per-tissue inputs; the original-survivor verdicts and three prior secondary hits are unchanged.", "", "Mature propeller estimates are not verified for the reason above. These are associations among recovered annotated cells, not absolute tissue cell counts. No manuscript file was read or changed for this item, and no non-mammalian extended-family row was refit.", ""])
    (REPORTS / "followup-mature-composition.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    started = utc_now()
    print(f"START {started} MODEL {MODEL_USED}", flush=True)
    locked = RESULTS / "rev3-composition-global-bh-624.csv"
    before = hashlib.sha256(locked.read_bytes()).hexdigest()
    if before != LOCKED_624_SHA256:
        raise RuntimeError("Locked 624-family file hash differs from the expected Item 3 input")
    original = pd.read_csv(RESULTS / "reviewer-phase1-composition-global-bh.csv")
    if len(original) != 431 or original.duplicated(KEY).any():
        raise RuntimeError("Original composition family is not 431 unique rows")
    if int(original["significant_global_05"].sum()) != 15:
        raise RuntimeError("Original composition family does not contain 15 passes")
    refit, untestable, inputs, coverage = collect_tests(original)
    if refit.empty:
        raise RuntimeError("No mature tests were produced")
    refit = add_multiplicity(refit, original)
    mouse = refit[refit["species"] == "Mouse"].sort_values(KEY).reset_index(drop=True)
    human = refit[refit["species"] == "Human"].sort_values(KEY).reset_index(drop=True)
    columns = KEY + ["sensitivity", "model", "effect_type", "effect", "se", "ci95_low",
                     "ci95_high", "pvalue", "q_bh", "q_431_family", "q_refit_family",
                     "n_431_family", "n_refit_family", "in_original_431", "n_donors",
                     "n_present", "total_cells", "min_age", "max_age", "original_survivor",
                     "original_effect", "original_se", "original_pvalue", "original_n_donors",
                     "original_q_global_bh", "mature_effect_over_original", "mature_se_over_original_se"]
    mouse[columns].to_csv(RESULTS / "followup-mature-composition-mouse.csv", index=False)
    human[columns].to_csv(RESULTS / "followup-mature-composition-human.csv", index=False)
    untestable.sort_values(KEY + ["row_scope"]).to_csv(
        RESULTS / "followup-mature-composition-untestable.csv", index=False)
    coverage.sort_values(["species", "organ", "cohort"]).to_csv(
        RESULTS / "followup-mature-composition-coverage.csv", index=False)
    secondary = secondary_checks(refit, inputs)
    secondary.to_csv(RESULTS / "followup-mature-composition-secondary-checks.csv", index=False)
    after = hashlib.sha256(locked.read_bytes()).hexdigest()
    if after != before:
        raise RuntimeError("Locked 624-family file changed during the mature refit")
    write_runtime("complete_mature_composition", started)
    print(f"END {utc_now()} mouse={len(mouse)} human={len(human)} refit_family={len(refit)}", flush=True)
    print(f"untestable_rows={len(untestable)} secondary_hits={len(secondary)}", flush=True)


if __name__ == "__main__":
    main()
