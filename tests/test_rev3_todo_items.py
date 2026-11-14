"""Tests for Rev3 overnight to-do list scripts and deliverables."""

from pathlib import Path

import numpy as np
import pandas as pd

from scripts.run_rev3_composition_fly_worm_zebrafish import (
    clr_table,
    fit_f_dist_smyth,
    run_clr_regressions,
    run_propeller_limma,
)
from scripts.run_rev3_entropy_mature_mouse_lung import compute_deman_kl, compute_shannon

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"


def test_clr_zero_sum_and_propeller_synthetic():
    counts = pd.DataFrame(
        {
            "cell_a": [100, 120, 150, 200, 220, 250, 300, 320],
            "cell_b": [400, 380, 350, 300, 280, 250, 200, 180],
            "cell_c": [50, 52, 48, 51, 49, 50, 52, 48],
        },
        index=[f"d{i}" for i in range(8)],
    )
    meta = pd.DataFrame(
        {
            "age": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0],
            "sex_male": [0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0],
        },
        index=counts.index,
    )
    clr = clr_table(counts, pseudocount=0.5)
    assert np.allclose(clr.sum(axis=1).to_numpy(), 0.0, atol=1e-12)

    s02, df0 = fit_f_dist_smyth(np.array([0.05, 0.06, 0.07]), df1=5.0)
    assert s02 > 0
    assert df0 > 0

    clr_rows = run_clr_regressions(counts, meta, "Fruit fly", "head", "AFCA", "all", "primary", 0.5)
    assert len(clr_rows) == 3
    assert np.isclose(sum(r["effect"] for r in clr_rows), 0.0, atol=1e-12)

    prop_rows = run_propeller_limma(counts, meta, "Fruit fly", "head", "AFCA", "all", "propeller", 0.5)
    assert len(prop_rows) == 3
    cell_a = [r for r in prop_rows if r["cell_type"] == "cell_a"][0]
    assert cell_a["effect"] > 0
    assert cell_a["q_within_bh"] < 0.05


def test_item1_outputs_and_global_bh_family():
    comp = pd.read_csv(RESULTS / "rev3-composition-fly-worm-zebrafish.csv")
    primary = comp[comp["is_primary_test"]]
    assert len(primary) == 193
    assert set(comp.loc[comp["is_primary_stratum"], "pseudocount"].dropna().unique()) >= {0.1, 0.5, 1.0}
    assert (comp["sensitivity"] == "propeller_limma_primary_design").sum() >= 193

    g624 = pd.read_csv(RESULTS / "rev3-composition-global-bh-624.csv")
    assert len(g624) == 624
    assert int((g624["source_family"] == "existing_431_mouse_human_nmr").sum()) == 431
    assert int((g624["source_family"] == "new_193_fly_worm_zebrafish").sum()) == 193


def test_item2_and_item3_reports_clean():
    for fname in [
        "rev3-laca-nmr-verification-REPORT.md",
        "rev3-young-vs-mature-audit-REPORT.md",
    ]:
        text = (RESULTS / fname).read_text(encoding="utf-8")
        assert "\u2014" not in text
        assert "\u2013" not in text
        assert len(text) > 500
    assert "AUDIT ONLY" in (RESULTS / "rev3-young-vs-mature-audit-REPORT.md").read_text(encoding="utf-8")


def test_item4_entropy_math_and_outputs():
    import scipy.sparse as sp

    mat = sp.csr_matrix(np.array([[10, 10, 10, 10], [37, 1, 1, 1]], dtype=float))
    h = compute_shannon(mat)
    kl = compute_deman_kl(mat, sqrt_counts=False)
    kl_sqrt = compute_deman_kl(mat, sqrt_counts=True)
    assert np.isclose(h[0], np.log(4.0))
    assert np.isclose(kl[0], 0.0, atol=1e-12)
    assert np.allclose(h + kl, np.log(4.0), atol=1e-12)
    assert kl_sqrt[1] < kl[1]

    out = pd.read_csv(RESULTS / "rev3-entropy-mature-only-mouse-lung.csv")
    don = pd.read_csv(RESULTS / "rev3-entropy-mature-only-mouse-lung-donors.csv")
    assert don["age_months"].min() >= 18.0
    assert don["donor"].nunique() == 10
    assert set(out["metric"].unique()) == {"shannon_entropy", "deman_kl_methods", "deman_kl_sqrt"}
    assert set(out["depth"].astype(str).unique()) == {"500", "1500", "3000", "full"}
    assert (out["mde80_15mo_sd"] > 0).all()
    txt = (RESULTS / "rev3-entropy-mature-only-mouse-lung-REPORT.md").read_text(encoding="utf-8")
    assert "\u2014" not in txt and "\u2013" not in txt


def test_item5_blocker_report_clean():
    txt = (RESULTS / "rev3-deman-replication-human-lung-REPORT.md").read_text(encoding="utf-8")
    assert "BLOCKED" in txt
    assert "GSE281219_RAW.tar" in txt
    assert "agingseurat.RData" in txt
    assert "\u2014" not in txt and "\u2013" not in txt


def test_item6_bayesian_vector_mature_only_deliverables():
    diag = pd.read_csv(RESULTS / "rev3-bayesian-vector-mature-only-diagnostics.csv")
    post = pd.read_csv(RESULTS / "rev3-bayesian-vector-mature-only-posterior.csv")
    prior = pd.read_csv(RESULTS / "rev3-bayesian-vector-mature-only-prior-predictive.csv")
    freq = pd.read_csv(RESULTS / "rev3-bayesian-vector-mature-only-frequentist.csv")
    assert len(diag) == 10 and bool(diag["convergence_ok"].all())
    assert int(diag["divergences"].sum()) == 0
    assert float(diag["max_rhat"].max()) <= 1.01
    assert bool(prior["observed_within_prior_998"].all())
    assert len(freq) == 10 and len(post) == 507
    txt = (REPO / "report" / "cross-species-bayesian-vector-mature-only.md").read_text(encoding="utf-8")
    assert "0.0389" in txt and "0.7068" in txt and "male_ge_24m" in txt
    assert "\u2014" not in txt and "\u2013" not in txt


def test_item7_bayesian_full_paper_deliverables():
    diag = pd.read_csv(RESULTS / "rev3-bayesian-full-diagnostics.csv")
    summ = pd.read_csv(RESULTS / "rev3-bayesian-full-summary.csv")
    prior = pd.read_csv(RESULTS / "rev3-bayesian-full-prior-predictive.csv")
    assert len(diag) == 50 and bool(diag["convergence_ok"].all())
    assert int(diag["divergences"].sum()) == 0
    assert float(diag["max_rhat"].max()) <= 1.01
    assert len(prior) == 50 and int(prior["observed_within_prior_998"].sum()) == 48
    assert set(summ["family"].unique()) == {
        "composition",
        "transport",
        "concordance",
        "disorder_entropy",
        "interneurons",
    }
    assert len(summ) == 138
    txt = (REPO / "report" / "cross-species-bayesian-full.md").read_text(encoding="utf-8")
    for tok in ["Naked mole-rat", "DYNC1I1", "FBN1", "0.9272", "0.1156"]:
        assert tok in txt
    assert "\u2014" not in txt and "\u2013" not in txt


def test_item8_randomness_variability_deliverables():
    df = pd.read_csv(RESULTS / "rev3-randomness-variability-test.csv")
    jg = pd.read_csv(RESULTS / "rev3-randomness-variability-joint-global.csv")
    assert len(df) == 148
    assert set(df["arm"].unique()) == {
        "primary_10cells_full_residualized",
        "sensitivity_10cells_full_unresidualized",
        "sensitivity_10cells_rarefied_1500_pseudobulk",
        "sensitivity_20cells_full_residualized",
        "sensitivity_10cells_male_only_residualized",
    }
    primary = df[df["arm"] == "primary_10cells_full_residualized"]
    assert len(primary) == 34
    assert set(primary["tissue"].unique()) == {"Kidney", "Lung", "Marrow"}
    for qcol in [
        "q_bh_eta2_within_arm",
        "q_bh_lvr_gt0_within_arm",
        "q_bh_lvr_lt0_within_arm",
        "q_bh_lvr_two_within_arm",
        "q_bh_tot_lvr_two_within_arm",
        "q_bh_disp_gt0_within_arm",
        "q_bh_disp_two_within_arm",
    ]:
        valid_q = df[qcol].dropna()
        assert len(valid_q) > 0 and float(valid_q.min()) >= 0.05
    assert len(jg) == 2
    assert (jg["n_global_mouse_perms"] == 70).all()
    assert (jg["n_strata"] == 34).all()
    txt = (RESULTS / "rev3-randomness-variability-test-REPORT.md").read_text(encoding="utf-8")
    for tok in [
        "cross-sectional",
        "18-F-50",
        "30-M-2",
        "-0.1233",
        "0.02857",
        "sensitivity_10cells_male_only_residualized",
    ]:
        assert tok in txt
    assert "\u2014" not in txt and "\u2013" not in txt
