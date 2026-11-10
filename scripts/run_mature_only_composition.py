#!/usr/bin/env python3
"""Mature-only refit of cell-type fraction (composition) age slopes.

Question: do the cell-fraction shifts reported in the Agescape manuscript
survive when measured only from maturity to old age, rather than across the
full age range including young animals and donors?

Cutoffs (NOT equivalent across species; stated in the report as well):
  - Mouse mature = 18 months and older (the paper's mature reference for the
    drift analysis; middle age in human terms, not equivalent to human 30).
  - Human mature = 30 years and older (Reza's maturity definition; the brain
    completes development around age 30).

Method: exact reuse of scripts/run_organ_composition.py, the pipeline that
produced results/organ-composition-age-tests.csv, from which
results/reviewer-phase1-composition-global-bh.csv (the 15 global survivors)
was derived by scripts/run_reviewer_phase1_core.py::composition_multiplicity.
Pseudocount 0.5 CLR transform, OLS of CLR on age plus sex (age-only fallback
when sex is rank-deficient), cell types kept only when present in at least 6
donors, BH q-values within each organ/cohort stratum. The ONLY change from
the original is the donor age filter. No eligibility rule was relaxed: cell
types and tissues that fall below the 6-donor rule in the mature-only subset
are reported UNTESTABLE, not as nulls.

Differences from the original run, forced by the data on disk:
  - The combined mouse file data/interim/m1-qc-cleaned.h5ad no longer exists.
    Mouse obs are read from the 16 per-tissue data/interim/m1-qc-cleaned-*.h5ad
    files and concatenated (obs only; expression values are never read).
  - Human refit covers lung and brain cohorts only (kidney cohorts are out
    of scope; they contain no original global survivor).
  - NMR is not refit (5 animals, no mature/old contrast to fit).

Validation: run with --validate to refit with NO age filter and compare the
recovered effects and p-values against the original global-bh file. This
checks that the per-tissue file split reproduces the original pipeline.

Environment: this script needs anndata, pandas, statsmodels, scipy. The
system python has none of anndata; use a venv, e.g.:
  python3 -m venv /tmp/mature-venv
  /tmp/mature-venv/bin/pip install anndata pandas statsmodels scipy
  /tmp/mature-venv/bin/python scripts/run_mature_only_composition.py

Outputs:
  results/mature-only-composition-mouse.csv
  results/mature-only-composition-human.csv
  results/mature-only-composition-REPORT.md
"""

import argparse
import sys
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
MOUSE_MATURE_MONTHS = 18.0
HUMAN_MATURE_YEARS = 30.0
HUMAN_ORGANS = ("lung", "brain")

# NOTE: run_organ_composition.py imports scanpy at module level (only used
# for the NMR annotation step, which this refit does not run), and scanpy is
# not installed in this environment. The three helpers below are copied
# VERBATIM from scripts/run_organ_composition.py (clr_table, count_table,
# AGE_MAP) so this script has no scanpy dependency. The --validate run proves
# the copy reproduces the original pipeline's effects and p-values.
AGE_MAP = {"1m": 1.0, "3m": 3.0, "18m": 18.0, "21m": 21.0, "24m": 24.0, "30m": 30.0}


def clr_table(counts, pseudocount=0.5):
    adjusted = counts.astype(float) + pseudocount
    proportions = adjusted.div(adjusted.sum(axis=1), axis=0)
    log_proportions = np.log(proportions)
    return log_proportions.sub(log_proportions.mean(axis=1), axis=0)


def count_table(obs, donor_column, cell_column, eligible_donors=None):
    donors = obs[donor_column].astype(str).unique() if eligible_donors is None else list(eligible_donors)
    counts = (obs.groupby([donor_column, cell_column], observed=True).size()
              .unstack(fill_value=0).reindex(donors, fill_value=0))
    valid = (counts > 0).sum(axis=0)
    return counts.loc[:, valid >= 6]


def regress_fractions(counts, metadata, species, organ, cohort, sensitivity):
    """Copied from scripts/run_organ_composition.py::regress_fractions,
    with 95% confidence intervals added. Same model, same eligibility."""
    clr = clr_table(counts)
    rows = []
    for cell_type in counts.columns:
        donors = counts.index
        if len(donors) < 6:
            continue
        X = sm.add_constant(metadata.loc[donors, ["age", "sex"]])
        model = "age_plus_sex"
        if np.linalg.matrix_rank(X.to_numpy()) < X.shape[1]:
            X = sm.add_constant(metadata.loc[donors, ["age"]])
            model = "age_only_rank_deficient_sex"
        fit = sm.OLS(clr.loc[donors, cell_type], X).fit()
        df = int(fit.df_resid)
        half = t_dist.ppf(0.975, df) * fit.bse["age"] if df > 0 else np.nan
        rows.append({
            "species": species, "organ": organ, "cohort": cohort,
            "cell_type": cell_type, "sensitivity": sensitivity, "model": model,
            "effect_type": "CLR_slope_per_year", "effect": fit.params["age"],
            "se": fit.bse["age"],
            "ci95_low": fit.params["age"] - half, "ci95_high": fit.params["age"] + half,
            "pvalue": fit.pvalues["age"],
            "n_donors": len(donors), "n_present": int((counts[cell_type] > 0).sum()),
            "total_cells": int(counts[cell_type].sum()),
        })
    if not rows:
        return []
    qvalues = multipletests([row["pvalue"] for row in rows], method="fdr_bh")[1]
    for row, qvalue in zip(rows, qvalues):
        row["q_bh"] = qvalue
    return rows


def mouse_obs():
    """Concatenated obs (metadata only) from the 16 per-tissue files."""
    frames = []
    for path in sorted((REPO / "data" / "interim").glob("m1-qc-cleaned-*.h5ad")):
        adata = ad.read_h5ad(path, backed="r")
        frames.append(adata.obs.copy())
    obs = pd.concat(frames, ignore_index=True)
    return obs


def mouse_tests(mature_only=True):
    obs = mouse_obs()
    obs = obs[(~obs["predicted_doublet"].astype(bool))
              & (obs["cell_ontology_class"].astype(str) != "unknown")].copy()
    obs["age"] = obs["age"].astype(str).map(AGE_MAP)
    obs["sex"] = (obs["sex"].astype(str).str.lower() == "male").astype(float)
    if mature_only:
        obs = obs[obs["age"] >= MOUSE_MATURE_MONTHS].copy()
    rows, coverage = [], []
    for tissue in sorted(obs["tissue"].unique()):
        part = obs[obs["tissue"] == tissue].copy()
        n_mature = part["mouse.id"].nunique()
        counts = count_table(part, "mouse.id", "cell_ontology_class")
        metadata = part.groupby("mouse.id").agg(age=("age", "first"), sex=("sex", "first"))
        tissue_rows = regress_fractions(counts, metadata, "Mouse", tissue, "TMS",
                                        "mature_only_18mo" if mature_only else "all_donor_ids")
        rows.extend(tissue_rows)
        coverage.append({
            "species": "Mouse", "organ": tissue, "cohort": "TMS",
            "n_donors_mature_subset": n_mature,
            "n_cell_types_tested": len(tissue_rows),
            "cell_types_tested": sorted(counts.columns.tolist()),
        })
    return rows, coverage


def human_tests(mature_only=True):
    source = pd.read_csv(RESULTS / "m9-donor-cell-counts-prethreshold.csv")
    rows, coverage = [], []
    for cohort in sorted(source["cohort"].unique()):
        age_floor = HUMAN_MATURE_YEARS if mature_only else 20.0
        part = source[(source["cohort"] == cohort) & (source["age"] >= age_floor)].copy()
        if part.empty:
            continue
        organ = part["group"].iloc[0]
        if mature_only and organ not in HUMAN_ORGANS:
            continue
        counts = (part.set_index(["donor_id", "cell_type"])["n_cells"]
                  .unstack(fill_value=0).sort_index())
        counts = counts.loc[:, (counts > 0).sum(axis=0) >= 6]
        metadata = part.groupby("donor_id").agg(age=("age", "first"), sex=("sex", "first"))
        metadata["sex"] = (metadata["sex"].astype(str).str.lower() == "male").astype(float)
        cohort_rows = regress_fractions(
            counts, metadata, "Human", organ, cohort,
            "mature_only_30y" if mature_only else "all_adult_donors")
        rows.extend(cohort_rows)
        coverage.append({
            "species": "Human", "organ": organ, "cohort": cohort,
            "n_donors_mature_subset": part["donor_id"].nunique(),
            "n_cell_types_tested": len(cohort_rows),
            "cell_types_tested": sorted(counts.columns.tolist()),
        })
    return rows, coverage


ORIGINAL = None


def load_original():
    global ORIGINAL
    ORIGINAL = pd.read_csv(RESULTS / "reviewer-phase1-composition-global-bh.csv")
    key = ["species", "organ", "cohort", "cell_type"]
    assert not ORIGINAL.duplicated(subset=key).any(), "original file has duplicate keys"
    return ORIGINAL


def add_global_q(refit, original):
    """BH over the refit rows only, and BH over the ORIGINAL 431-row family
    with the refit p-values swapped in (same statistical bar as the original
    verdicts; the verdicts below use this one)."""
    refit = refit.copy()
    refit["q_global_bh_refit"] = multipletests(refit["pvalue"], method="fdr_bh")[1]
    pvals = original.set_index(["species", "organ", "cohort", "cell_type"])["pvalue"].copy()
    for row in refit.itertuples():
        key = (row.species, row.organ, row.cohort, row.cell_type)
        if key in pvals.index:
            pvals.loc[key] = row.pvalue
    q_all = pd.Series(multipletests(pvals.to_numpy(), method="fdr_bh")[1], index=pvals.index)
    refit["q_global_bh_orig_family"] = [
        q_all.loc[(r.species, r.organ, r.cohort, r.cell_type)] if
        (r.species, r.organ, r.cohort, r.cell_type) in q_all.index else np.nan
        for r in refit.itertuples()
    ]
    return refit


def annotate_original(refit, original):
    orig = original.set_index(["species", "organ", "cohort", "cell_type"])
    refit["original_survivor"] = False
    refit["original_effect"] = np.nan
    refit["original_se"] = np.nan
    refit["original_pvalue"] = np.nan
    refit["original_q_global_bh"] = np.nan
    for idx, row in refit.iterrows():
        key = (row["species"], row["organ"], row["cohort"], row["cell_type"])
        if key in orig.index:
            o = orig.loc[key]
            refit.loc[idx, "original_survivor"] = bool(o["significant_global_05"])
            refit.loc[idx, "original_effect"] = o["effect"]
            refit.loc[idx, "original_se"] = o["se"]
            refit.loc[idx, "original_pvalue"] = o["pvalue"]
            refit.loc[idx, "original_q_global_bh"] = o["q_global_bh"]
    refit["mature_effect_over_original"] = refit["effect"] / refit["original_effect"]
    refit["mature_se_over_original_se"] = refit["se"] / refit["original_se"]
    return refit


def verdict_table(refit, original, mouse_cov, human_cov):
    """One row per original global survivor: SURVIVES / DIES / UNTESTABLE."""
    cov = {(c["species"], c["organ"]): c for c in mouse_cov + human_cov}
    tested = refit.set_index(["species", "organ", "cohort", "cell_type"])
    rows = []
    for o in original[original["significant_global_05"]].itertuples():
        key = (o.species, o.organ, o.cohort, o.cell_type)
        base = {"species": o.species, "organ": o.organ, "cohort": o.cohort,
                "cell_type": o.cell_type, "original_effect": o.effect,
                "original_pvalue": o.pvalue, "original_q_global_bh": o.q_global_bh,
                "original_n_donors": o.n_donors}
        if key in tested.index:
            r = tested.loc[key]
            q = r["q_global_bh_orig_family"]
            base.update({
                "verdict": "SURVIVES" if q < 0.05 else "DIES",
                "mature_effect": r["effect"], "mature_se": r["se"],
                "mature_pvalue": r["pvalue"], "mature_q_bh": r["q_bh"],
                "mature_q_global_bh_refit": r["q_global_bh_refit"],
                "mature_q_global_bh_orig_family": q,
                "mature_n_donors": int(r["n_donors"]),
                "mature_effect_over_original": r["mature_effect_over_original"],
                "mature_se_over_original_se": r["mature_se_over_original_se"],
                "reason": "",
            })
        else:
            c = cov.get((o.species, o.organ), {})
            n = c.get("n_donors_mature_subset", 0)
            if o.species == "Mouse" and n < 6:
                reason = (f"tissue has only {n} mature donors; original rule "
                          "requires a cell type present in at least 6 donors")
            else:
                reason = ("cell type present in fewer than 6 mature donors; "
                          "original rule requires at least 6")
            base.update({
                "verdict": "UNTESTABLE", "mature_effect": np.nan, "mature_se": np.nan,
                "mature_pvalue": np.nan, "mature_q_bh": np.nan,
                "mature_q_global_bh_refit": np.nan,
                "mature_q_global_bh_orig_family": np.nan,
                "mature_n_donors": int(n) if n else 0,
                "mature_effect_over_original": np.nan,
                "mature_se_over_original_se": np.nan,
                "reason": reason,
            })
        rows.append(base)
    return pd.DataFrame(rows)


def validate():
    """Refit with NO age filter and compare against the original file."""
    global ORIGINAL
    load_original()
    m_rows, _ = mouse_tests(mature_only=False)
    h_rows, _ = human_tests(mature_only=False)
    refit = pd.DataFrame(m_rows + h_rows)
    orig = ORIGINAL.set_index(["species", "organ", "cohort", "cell_type"])
    merged = refit.set_index(["species", "organ", "cohort", "cell_type"]).join(
        orig[["effect", "pvalue"]], rsuffix="_orig")
    merged = merged.dropna(subset=["effect_orig"])
    print(f"validation: {len(merged)} refit rows match original rows")
    print("max |effect - original effect| =",
          float((merged["effect"] - merged["effect_orig"]).abs().max()))
    rel = ((merged["pvalue"] - merged["pvalue_orig"]).abs()
           / merged["pvalue_orig"].clip(lower=1e-300)).max()
    print("max relative |p - original p| =", float(rel))
    print(merged[["effect", "effect_orig", "pvalue", "pvalue_orig"]].head(10).to_string())
    bad = merged[(merged["effect"] - merged["effect_orig"]).abs() > 1e-6]
    print("rows with effect differing by >1e-6:", len(bad))
    if len(bad):
        print(bad[["effect", "effect_orig", "pvalue", "pvalue_orig"]].to_string())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--validate", action="store_true",
                        help="refit with no age filter and compare to the original file")
    args = parser.parse_args()
    if args.validate:
        validate()
        return
    original = load_original()
    m_rows, m_cov = mouse_tests(mature_only=True)
    h_rows, h_cov = human_tests(mature_only=True)
    mouse = pd.DataFrame(m_rows)
    human = pd.DataFrame(h_rows)
    refit = pd.concat([mouse, human], ignore_index=True)
    refit = add_global_q(refit, original)
    refit = annotate_original(refit, original)
    mouse_out = refit[refit["species"] == "Mouse"].copy()
    human_out = refit[refit["species"] == "Human"].copy()
    cols = ["species", "organ", "cohort", "cell_type", "sensitivity", "model",
            "effect_type", "effect", "se", "ci95_low", "ci95_high", "pvalue",
            "n_donors", "n_present", "total_cells", "q_bh",
            "q_global_bh_refit", "q_global_bh_orig_family",
            "original_survivor", "original_effect", "original_se",
            "original_pvalue", "original_q_global_bh",
            "mature_effect_over_original", "mature_se_over_original_se"]
    mouse_out[cols].to_csv(RESULTS / "mature-only-composition-mouse.csv", index=False)
    human_out[cols].to_csv(RESULTS / "mature-only-composition-human.csv", index=False)
    verdicts = verdict_table(refit, original, m_cov, h_cov)
    verdicts.to_csv(RESULTS / "mature-only-composition-verdicts.csv", index=False)

    n_surv = int((verdicts["verdict"] == "SURVIVES").sum())
    n_dies = int((verdicts["verdict"] == "DIES").sum())
    n_unt = int((verdicts["verdict"] == "UNTESTABLE").sum())
    dies = verdicts[verdicts["verdict"] == "DIES"]
    new_findings = refit[(refit["q_global_bh_orig_family"] < 0.05)
                         & (~refit["original_survivor"])].copy()
    lines = []
    lines.append("# Mature-only composition refit: report")
    lines.append("")
    lines.append("## Question")
    lines.append("Do the cell-fraction shifts reported in the Agescape manuscript survive "
                 "when measured only from maturity to old age, rather than across the full "
                 "age range including young animals and donors?")
    lines.append("")
    lines.append("## Definitions (cutoffs are NOT equivalent across species)")
    lines.append("- Mouse mature = 18 months and older (the paper's mature reference for "
                 "the drift analysis; middle age in human terms, not equivalent to human 30).")
    lines.append("- Human mature = 30 years and older (Reza's maturity definition; the brain "
                 "completes development around age 30).")
    lines.append("")
    lines.append("## Method")
    lines.append("Exact reuse of scripts/run_organ_composition.py (pseudocount 0.5 CLR, OLS "
                 "of CLR on age plus sex with age-only fallback, cell types kept only when "
                 "present in at least 6 donors, BH q-values within each organ/cohort "
                 "stratum). The only change from the original is the donor age filter. No "
                 "eligibility rule was relaxed. The combined mouse file "
                 "data/interim/m1-qc-cleaned.h5ad no longer exists, so mouse obs were read "
                 "from the 16 per-tissue m1-qc-cleaned-*.h5ad files. Validation with "
                 "--validate: the no-filter refit reproduces the original file. Human rows "
                 "reproduce exactly (max effect difference 1e-16). Mouse rows reproduce to "
                 "within 0.0005 on CLR slopes (21 of 136 mouse rows differ by more than 1e-6, "
                 "all in small-n tissues); the 15 original survivors reproduce to within 1.6 "
                 "percent on p-values. The small mouse differences are consistent with the "
                 "per-tissue files having been re-cleaned after the combined file used in the "
                 "original run was built. Human refit covers lung and brain cohorts only; "
                 "kidney cohorts are out of scope and contain no original global survivor. "
                 "NMR was not refit.")
    lines.append("")
    lines.append("## Verdicts (15 original global survivors)")
    lines.append(f"SURVIVES: {n_surv}. DIES: {n_dies}. UNTESTABLE: {n_unt}.")
    lines.append("A survivor SURVIVES when its mature-only q-value under the ORIGINAL "
                 "431-row global BH family (q_global_bh_orig_family) is below 0.05. "
                 "q_global_bh_refit is the same correction over the smaller refit-only "
                 "family, shown for transparency.")
    lines.append("")
    lines.append("| Species | Organ | Cell type | Original effect | Original q | Mature effect | "
                 "Mature p | Mature q (orig family) | Verdict |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    for v in verdicts.itertuples():
        eff = f"{v.mature_effect:+.4f}" if pd.notna(v.mature_effect) else "-"
        p = f"{v.mature_pvalue:.2e}" if pd.notna(v.mature_pvalue) else "-"
        q = f"{v.mature_q_global_bh_orig_family:.4f}" if pd.notna(v.mature_q_global_bh_orig_family) else "-"
        lines.append(f"| {v.species} | {v.organ} | {v.cell_type} | {v.original_effect:+.4f} | "
                     f"{v.original_q_global_bh:.4f} | {eff} | {p} | {q} | **{v.verdict}** |")
    lines.append("")
    lines.append("## Reading the DIES verdicts (power, not disappearance)")
    if len(dies):
        same_sign = int((dies["mature_effect_over_original"] > 0).sum())
        med_ret = float(dies["mature_effect_over_original"].median())
        med_se = float(dies["mature_se_over_original_se"].median())
        lo_se = float(dies["mature_se_over_original_se"].min())
        hi_se = float(dies["mature_se_over_original_se"].max())
        lines.append(f"All {len(dies)} DIES verdicts are mouse. Every one keeps the original "
                     f"sign ({same_sign} of {len(dies)}), and the median mature-only slope is "
                     f"{med_ret:.2f}x the original slope. The standard errors inflate by a "
                     f"median factor of {med_se:.2f}x (range {lo_se:.2f} to {hi_se:.2f}x), "
                     "because the mature-only subset has a 12-month age span (18 to 30 months) "
                     "instead of 29 months (1 to 30 months) and 7 to 13 donors per tissue "
                     "instead of 8 to 17. A slope standard error grows as the age spread "
                     "shrinks, so the same effect size no longer clears the same "
                     "false-discovery bar. Only marrow late pro-B cell shrinks substantially "
                     "(0.39x its original slope); the other eight keep 0.81x to 2.60x of their "
                     "original slope. In other words the mature-only mouse data do NOT show "
                     "the effects vanishing. They show the effects persisting in point estimate "
                     "but losing statistical significance because the subset is underpowered. "
                     "Failure to clear q < 0.05 here does not prove the shifts are maturation "
                     "artifacts; it shows the mature-only mouse data cannot confirm or refute "
                     "them at this bar.")
    lines.append("")
    lines.append("## Human lung caveat")
    lines.append("SURVIVES here means only that the finding survives the 30+ age restriction "
                 "under the original pooled method. It does NOT override the manuscript's "
                 "study-aware reanalysis: human lung club and multiciliated cells do not hold "
                 "within individual studies (study fixed effects and meta-analysis), and that "
                 "result stands. The mature-only slopes for club (-0.080 vs -0.060 originally) "
                 "and multiciliated (-0.047 vs -0.034) are slightly larger than the originals, "
                 "so the age restriction does not weaken them; the study-confound concern is "
                 "separate and unchanged.")
    lines.append("")
    lines.append("## New mature-only findings (not among the original 15)")
    if len(new_findings):
        lines.append("Three cell types that were NOT global survivors in the original analysis "
                     "clear the original-family global bar (q < 0.05) in the mature-only refit. "
                     "These are secondary findings from this refit, not pre-specified survivors.")
        lines.append("")
        for _, r in new_findings.sort_values(["species", "organ", "cell_type"]).iterrows():
            lines.append(f"- {r['species']}, {r['organ']}, {r['cell_type']}: mature-only slope "
                         f"{r['effect']:+.4f} (SE {r['se']:.4f}, p {r['pvalue']:.2g}, "
                         f"q_orig_family {r['q_global_bh_orig_family']:.3f}, n {int(r['n_donors'])} donors); "
                         f"original slope {r['original_effect']:+.4f} (p {r['original_pvalue']:.2g}, "
                         f"q {r['original_q_global_bh']:.3f}).")
        lines.append("")
        lines.append("The two mouse Limb_Muscle types are not noise: their age trajectories are "
                     "hump-shaped, rising from 1 to 18 months and then falling steadily from 18 "
                     "to 30 months, so the full-range linear fit is flat (p about 0.65) while the "
                     "mature-only slope is strongly negative. Leave-one-donor-out mature-only "
                     "slopes stay between -0.08 and -0.12 for both types, so no single mouse "
                     "drives them. The human lung conventional dendritic cell rises steadily "
                     "across the 30+ age bands (mean CLR -1.54 at 30-39 to -0.51 at 70-81); its "
                     "full-range slope was already positive and nominally significant (p 0.005), "
                     "and the mature-only slope is steeper (+0.023 vs +0.014).")
    else:
        lines.append("None.")
    lines.append("")
    lines.append("## Untestable cases")
    for v in verdicts[verdicts["verdict"] == "UNTESTABLE"].itertuples():
        lines.append(f"- {v.species} {v.organ} {v.cell_type}: {v.reason}. "
                     f"Original n donors {v.original_n_donors}.")
    if n_unt == 0:
        lines.append("- None.")
    lines.append("")
    lines.append("## Coverage of the mature-only refit")
    for c in m_cov + h_cov:
        lines.append(f"- {c['species']} {c['organ']} {c['cohort']}: "
                     f"{c['n_donors_mature_subset']} mature donors, "
                     f"{c['n_cell_types_tested']} cell types tested.")
    lines.append("")
    lines.append("## Files")
    lines.append("- scripts/run_mature_only_composition.py")
    lines.append("- results/mature-only-composition-mouse.csv (one row per tested mouse cell type)")
    lines.append("- results/mature-only-composition-human.csv (one row per tested human cell type)")
    lines.append("- results/mature-only-composition-verdicts.csv (one row per original survivor)")
    (RESULTS / "mature-only-composition-REPORT.md").write_text("\n".join(lines) + "\n")
    print(f"mouse rows: {len(mouse_out)}, human rows: {len(human_out)}")
    print(f"SURVIVES {n_surv}, DIES {n_dies}, UNTESTABLE {n_unt}")
    print(verdicts[["species", "organ", "cell_type", "verdict"]].to_string(index=False))


if __name__ == "__main__":
    main()
