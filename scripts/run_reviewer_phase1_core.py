#!/usr/bin/env python3

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import binomtest, hypergeom, pearsonr, spearmanr, t
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_cross_species_expression_vector import fit_global

RESULTS = REPO / "results"
PREFIX = "reviewer-phase1-"


def write(frame, name):
    frame.to_csv(RESULTS / f"{PREFIX}{name}.csv", index=False)


def mouse_drift():
    donors = pd.read_csv(RESULTS / "h3-cross-species-vector-donors.csv")
    mouse = donors[donors["species"] == "Mouse"].copy()
    baseline = fit_global(mouse, "profile_hellinger", "Mouse", "All tested tissues", "All cohorts",
                          "unadjusted", 0.4, "all_donor_ids")
    rows = []
    for cell_type in sorted(mouse["cell_type"].unique()):
        selected = mouse[mouse["cell_type"] != cell_type]
        fit = fit_global(selected, "profile_hellinger", "Mouse", "All tested tissues", "All cohorts",
                         "unadjusted", 0.4, "all_donor_ids")
        if fit:
            rows.append({
                "cell_type": cell_type, "n_units_removed": len(
                    mouse.loc[mouse["cell_type"] == cell_type, ["tissue", "cell_type"]].drop_duplicates()
                ),
                "effect_sd": fit["effect_sd"], "se": fit["se"], "pvalue": fit["pvalue"],
                "delta_from_full": fit["effect_sd"] - baseline["effect_sd"],
            })
    write(pd.DataFrame(rows).sort_values("delta_from_full"), "mouse-drift-leave-one-cell-type-out")
    tests = pd.read_csv(RESULTS / "h3-cross-species-vector-tests.csv")
    units = tests[(tests["species"] == "Mouse") & (tests["metric"] == "profile_hellinger")
                  & (tests["model"] == "unadjusted")].copy()
    summary = pd.DataFrame([{
        "n_units": len(units), "n_positive": int((units["slope_per_year"] > 0).sum()),
        "fraction_positive": float((units["slope_per_year"] > 0).mean()),
        "sign_test_p": binomtest(int((units["slope_per_year"] > 0).sum()), len(units), 0.5).pvalue,
        "full_effect_sd": baseline["effect_sd"], "loo_min_effect_sd": min(row["effect_sd"] for row in rows),
        "loo_max_effect_sd": max(row["effect_sd"] for row in rows),
        "max_absolute_loo_delta": max(abs(row["delta_from_full"]) for row in rows),
    }])
    write(summary, "mouse-drift-summary")
    write(units, "mouse-drift-unit-effects")


def cluster_nmr_effect():
    donors = pd.read_csv(RESULTS / "h3-cross-species-vector-donors.csv")
    data = donors[donors["species"] == "Naked mole-rat"].copy()
    data["unit"] = data["cell_type"]
    data["age_scaled"] = data["age"] / 3.1
    standardized = []
    for _, part in data.groupby("unit"):
        part = part.copy()
        part["outcome"] = ((part["profile_hellinger"] - part["profile_hellinger"].mean())
                           / part["profile_hellinger"].std(ddof=1))
        standardized.append(part)
    data = pd.concat(standardized)
    X = sm.add_constant(pd.concat([
        data[["age_scaled"]], pd.get_dummies(data["unit"], drop_first=True, dtype=float)
    ], axis=1))
    fit = sm.OLS(data["outcome"], X).fit(
        cov_type="cluster", cov_kwds={"groups": data["donor"]}, use_t=True
    )
    return fit.params["age_scaled"], fit.bse["age_scaled"], 4


def null_row(claim, effect, se, df, source, source_line, margin):
    ci = t.interval(0.95, df, loc=effect, scale=se)
    lower_t = (effect + margin) / se
    upper_t = (effect - margin) / se
    tost_p = max(t.sf(lower_t, df), t.cdf(upper_t, df))
    mde = (t.ppf(0.975, df) + t.ppf(0.8, df)) * se
    return {
        "claim": claim, "effect": effect, "se": se, "df": df,
        "ci95_low": ci[0], "ci95_high": ci[1], "equivalence_margin_sd": margin,
        "tost_p": tost_p, "equivalent_at_0.05": tost_p < 0.05,
        "mde80_two_sided_0.05": mde, "source": source, "source_line": source_line,
    }


def nulls(margin):
    rows = []
    vector = pd.read_csv(RESULTS / "h3-cross-species-vector-global.csv")
    row = vector[(vector["species"] == "Human") & (vector["tissue"] == "All tested tissues")
                 & (vector["metric"] == "profile_hellinger") & (vector["model"] == "unadjusted")].iloc[0]
    rows.append(null_row("Human whole-vector displacement", row["effect_sd"], row["se"],
                         int(row["n_donors"] - 1), "results/h3-cross-species-vector-global.csv", 110, margin))
    disorder = pd.read_csv(RESULTS / "h3-frequentist-global.csv")
    for claim, metric, line in [("Entropy", "entropy", 2), ("Off-program expression", "offprogram_fraction", 5)]:
        row = disorder[(disorder["metric"] == metric)
                       & (disorder["analysis"] == "all_cells_all_donor_ids")].iloc[0]
        df = int(row["n_donors"] - 1)
        se = (row["ci95_high_sd15"] - row["ci95_low_sd15"]) / (2 * t.ppf(0.975, df))
        rows.append(null_row(claim, row["effect_15mo_sd"], se, df,
                             "results/h3-frequentist-global.csv", line, margin))
    proliferation = pd.read_csv(RESULTS / "proliferation-moderator-tests.csv")
    for claim, analysis, response, line in [
        ("Proliferation expression magnitude", "expression_vector", "absolute_effect_sd", 15),
        ("Proliferation cell-fraction magnitude", "cell_fraction", "absolute_effect_sd", 43),
    ]:
        row = proliferation[(proliferation["analysis"] == analysis) & (proliferation["response"] == response)
                            & (proliferation["scope"] == "All species combined")
                            & (proliferation["sensitivity"].isin(["depth_adjusted", "all_donors"]))].iloc[0]
        df = int(row["n_units"] - row["n_strata"] - 1)
        rows.append(null_row(claim, row["post_minus_proliferative"], row["se"], df,
                             "results/proliferation-moderator-tests.csv", line, margin))
    effect, se, df = cluster_nmr_effect()
    rows.append(null_row("NMR expression-vector displacement", effect, se, df,
                         "results/h3-cross-species-vector-donors.csv", "2-6788; model refit", margin))
    fish = pd.read_csv(RESULTS / "extended-zebrafish-vector-fine-global.csv")
    row = fish[(fish["model"] == "depth_adjusted") & (fish["tissue"] == "All tested tissues")].iloc[0]
    rows.append(null_row("Zebrafish fine expression-vector displacement", row["effect_sd"], row["se"],
                         int(row["n_donors"] - 1), "results/extended-zebrafish-vector-fine-global.csv",
                         "depth-adjusted all-tissues row", margin))
    result = pd.DataFrame(rows)
    write(result, "nulls-intervals-power-equivalence")
    return result


def composition_multiplicity():
    core = pd.read_csv(RESULTS / "organ-composition-age-tests.csv")
    primary = core[core["sensitivity"].isin(["all_donor_ids", "all_adult_donors", "annotation_resolution_0.4"])].copy()
    primary["q_global_bh"] = multipletests(primary["pvalue"], method="fdr_bh")[1]
    primary["significant_global_05"] = primary["q_global_bh"] < 0.05
    write(primary, "composition-global-bh")
    summary = primary.groupby(["species", "organ", "cohort"]).agg(
        n_cell_types=("cell_type", "nunique"),
        n_significant_original=("q_bh", lambda values: int((values < 0.05).sum())),
        n_significant_global=("significant_global_05", "sum"),
        min_q_global=("q_global_bh", "min"),
    ).reset_index()
    write(summary, "composition-global-bh-by-organ")
    signature = pd.read_csv(RESULTS / "human-transport-global-multiplicity.csv")
    signature["q_overlap_aware_recheck"] = multipletests(signature["p"], method="fdr_bh")[1]
    write(signature, "signature-global-bh-recheck")


def interneurons():
    significant = pd.read_csv(RESULTS / "h3-human-neuronal-gene-replication-significant.csv")
    all_genes = pd.read_csv(RESULTS / "h3-human-neuronal-gene-replication-all.csv")
    rows = []
    for cell_type, part in all_genes.groupby("cell_type"):
        discovery = part[part["mssm_fdr"]]
        n = len(part)
        mssm = int(part["mssm_fdr"].sum())
        hbcc = int(part["hbcc_fdr"].sum())
        both = int((part["mssm_fdr"] & part["hbcc_fdr"] & part["same_direction"]).sum())
        for subset, selected in [("all_shared_genes", part), ("mssm_fdr_genes", discovery)]:
            if len(selected) < 3:
                continue
            rows.append({
                "cell_type": cell_type, "subset": subset, "n_genes": len(selected),
                "pearson_effect_r": pearsonr(selected["slope_log_cp10k_per_year_mssm"],
                                               selected["slope_log_cp10k_per_year_hbcc"]).statistic,
                "spearman_effect_rho": spearmanr(selected["slope_log_cp10k_per_year_mssm"],
                                                   selected["slope_log_cp10k_per_year_hbcc"]).statistic,
                "same_direction": int(selected["same_direction"].sum()),
                "direction_concordance": selected["same_direction"].mean(),
                "direction_binomial_p": binomtest(int(selected["same_direction"].sum()),
                                                   len(selected), 0.5).pvalue,
                "n_universe": n, "n_mssm_fdr": mssm, "n_hbcc_fdr": hbcc,
                "both_fdr_same_direction": both,
                "expected_both_fdr_under_independence": mssm * hbcc / n,
                "chance_overlap_hypergeometric_p": hypergeom.sf(both - 1, n, hbcc, mssm),
            })
    write(pd.DataFrame(rows), "interneuron-replication-statistics")
    subtypes = ["VIP GABAergic cortical interneuron", "pvalb GABAergic cortical interneuron",
                "sst GABAergic cortical interneuron"]
    sets = {name: set(significant[(significant["cell_type"] == name)
                                  & significant["mssm_fdr"]]["gene"]) for name in subtypes}
    universe = len(all_genes[all_genes["cell_type"] == subtypes[0]])
    rows = []
    for i, first in enumerate(subtypes):
        for second in subtypes[i + 1:]:
            overlap = len(sets[first] & sets[second])
            rows.append({
                "cell_type_1": first, "cell_type_2": second, "n_1": len(sets[first]),
                "n_2": len(sets[second]), "overlap": overlap,
                "jaccard": overlap / len(sets[first] | sets[second]),
                "expected_overlap_under_independence": len(sets[first]) * len(sets[second]) / universe,
                "hypergeometric_p": hypergeom.sf(overlap - 1, universe, len(sets[second]), len(sets[first])),
                "n_universe": universe,
            })
    triple = set.intersection(*sets.values())
    rows.append({"cell_type_1": "VIP", "cell_type_2": "pvalb and sst", "n_1": len(sets[subtypes[0]]),
                 "n_2": len(sets[subtypes[1]] | sets[subtypes[2]]), "overlap": len(triple),
                 "jaccard": np.nan, "expected_overlap_under_independence": np.nan,
                 "hypergeometric_p": np.nan, "n_universe": universe})
    write(pd.DataFrame(rows), "interneuron-subtype-overlap")
    top = significant[significant["mssm_fdr"] & significant["cell_type"].isin(subtypes)].copy()
    top["rank"] = top.groupby("cell_type")["q_bh_mssm"].rank(method="first")
    shared = set(pd.read_csv(RESULTS / "h3-human-mssm-shared-neuronal-genes.csv")["gene"])
    top = top[(top["rank"] <= 10) | top["gene"].isin(shared)]
    write(top.sort_values(["cell_type", "q_bh_mssm"]), "interneuron-dotplot-data")


def inventory():
    rows = [
        ("Mouse aging DE", "M2 within unit and gene", "BH within cell type", 0.05, "primary", "results/m2_de_*.csv"),
        ("Mouse aging ranking", "M2 cell-type aging signal", "BH and permutation audit", 0.05, "primary", "results/m2-aging-signal-ranked.csv"),
        ("Human signature transport", "252 M9/M10 cell-type tests", "BH across 252", 0.05, "primary", "results/human-transport-global-multiplicity.csv"),
        ("PsychAD adult signature transport", "27 adult cell types", "BH across 27", 0.05, "secondary", "results/m10-celltype-stratified-summary.csv"),
        ("Core organ composition", "Mouse, human, NMR cell types", "BH within organ or cohort; global recheck added", 0.05, "primary", "results/organ-composition-age-tests.csv"),
        ("Core expression vectors", "Cell types within tissue or cohort", "BH within tissue or cohort", 0.05, "primary", "results/h3-cross-species-vector-tests.csv"),
        ("Entropy and off-program", "Nine cell types per metric", "BH within metric", 0.05, "primary", "results/h3-transcriptional-disorder.csv"),
        ("Entropy and off-program sensitivity", "Alternative cells, genes, and endpoints", "BH within sensitivity family", 0.05, "secondary", "results/h3-transcriptional-disorder.csv"),
        ("Bayesian disorder and vector models", "Posterior effects", "94% HDI; no p-value correction", np.nan, "exploratory", "results/h3-bayesian-posterior.csv"),
        ("Fly composition and vectors", "Cell types within head or body", "BH within tissue and endpoint", 0.05, "primary", "results/extended-afca-*.csv"),
        ("Gao worm composition and vectors", "Cell groups", "BH plus platform-stratified permutation", 0.05, "primary", "results/extended-worm-*.csv"),
        ("Roux worm", "Technical replicate series", "No inferential correction", np.nan, "exploratory", "results/extended-roux-*.csv"),
        ("Zebrafish composition and vectors", "Broad and fine cell types", "BH and exact permutation", 0.05, "primary", "results/extended-zebrafish-*.csv"),
        ("Proliferation moderator", "Species-specific signed and absolute contrasts", "BH across five species per response", 0.05, "exploratory", "results/proliferation-moderator-tests.csv"),
        ("Methylation clock", "Predictive performance", "Not applicable", np.nan, "benchmark", "results/h1-epigenetic-clock-*.csv"),
        ("Reviewer Phase 1", "All new reviewer-response tests", "Specified per output family", 0.05, "exploratory", "results/reviewer-phase1-*.csv"),
    ]
    write(pd.DataFrame(rows, columns=["hypothesis", "family", "correction", "threshold", "status", "source"]),
          "multiplicity-inventory")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--equivalence-margin", type=float, nargs="+", default=[0.1, 0.25])
    args = parser.parse_args()
    mouse_drift()
    null_results = [nulls(margin) for margin in args.equivalence_margin]
    write(pd.concat(null_results, ignore_index=True), "nulls-intervals-power-equivalence")
    composition_multiplicity()
    interneurons()
    inventory()


if __name__ == "__main__":
    main()
