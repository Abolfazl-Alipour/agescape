"""Tests for the mature PBMC rerun at 30-plus.

Every headline number in reports/followup-mature-pbmc.md is asserted
here from the result CSVs. Run from the repo root:
    python -m pytest tests/test_followup_mature_pbmc.py -q
"""

from pathlib import Path

import pandas as pd
import pytest
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
PREFIX = "followup-mature-pbmc"


def load(name):
    return pd.read_csv(RESULTS / f"{PREFIX}-{name}")


def row(frame, **match):
    subset = frame
    for key, value in match.items():
        subset = subset[subset[key] == value]
    assert len(subset) == 1, f"expected one row for {match}, got {len(subset)}"
    return subset.iloc[0]


def check_bh(frame, groups, pcol="pvalue", qcol="q_bh"):
    checked = 0
    for _, indices in frame.groupby(groups, dropna=False).groups.items():
        valid = frame.loc[indices, pcol].dropna().index
        if not len(valid):
            continue
        expected = multipletests(frame.loc[valid, pcol], method="fdr_bh")[1]
        for idx, value in zip(valid, expected):
            assert frame.loc[idx, qcol] == pytest.approx(value, abs=1e-9)
            checked += 1
    return checked


def test_bh_recomputation_zero_mismatches():
    total = 0
    total += check_bh(load("p1m-models.csv"), ["cohort", "endpoint", "design", "model"])
    total += check_bh(
        load("p1m-models.csv"), ["cohort", "endpoint", "design", "model"],
        pcol="random_pvalue", qcol="random_q_bh",
    )
    total += check_bh(load("clr-models.csv"), ["cohort", "endpoint", "design", "model"])
    total += check_bh(load("clr-targeted.csv"), ["cohort", "endpoint", "design", "model"])
    total += check_bh(load("hellinger-models.csv"), ["cohort", "endpoint", "design", "model"])
    total += check_bh(load("pseudobulk-targeted.csv"), ["cohort", "endpoint", "design", "model"])
    assert total == 438


def test_donor_counts_and_eligibility():
    counts = load("donor-counts.csv")
    assert len(counts) == 2
    zhu = row(counts, cohort="GSE213516")
    assert (zhu["n_total_selected"], zhu["n_mature_30_39"], zhu["n_old_60_plus"]) == (10, 3, 7)
    assert zhu["status"] == "TESTABLE_AT_COHORT_LEVEL"
    thomson = row(counts, cohort="GSE214546")
    assert (thomson["n_total_selected"], thomson["n_mature_30_39"], thomson["n_old_60_plus"]) == (20, 11, 9)
    assert thomson["status"] == "TESTABLE_AT_COHORT_LEVEL"
    eligibility = load("eligibility-all.csv")
    for cohort, expected_total, expected_eligible in (("GSE213516", 28, 22), ("GSE214546", 13, 13)):
        subset = eligibility[(eligibility["cohort"] == cohort) & (eligibility["endpoint"] == "eligibility")]
        assert len(subset) == expected_total
        assert int(subset["eligible"].sum()) == expected_eligible
    runtime = load("runtime.csv")
    assert set(runtime["status"]) == {"completed"}
    assert len(runtime) == 2


def test_p1m_transport_zero_both_pass():
    models = load("p1m-models.csv")
    primary = models[(models["model"] == "primary") & (models["status"] == "fit")]
    assert len(primary) == 35
    assert (primary["cohort"] == "GSE213516").sum() == 22
    assert (primary["cohort"] == "GSE214546").sum() == 13
    assert ((primary["q_bh"] < 0.05) & (primary["random_q_bh"] < 0.05)).sum() == 0
    assert (primary["q_bh"] < 0.05).sum() == 0
    assert primary["q_bh"].min() == pytest.approx(0.6299204155120731, rel=1e-9)
    # P1m gene coverage: 1007 up + 709 down detected of 1724 human-mapped.
    assert set(primary["p1m_up_detected"]) == {1007}
    assert set(primary["p1m_down_detected"]) == {709}
    assert primary["p1m_coverage"].iloc[0] == pytest.approx(1716 / 1724, rel=1e-9)
    assert set(primary["random_iterations"].dropna()) == {200.0}
    # Closest composite call: Zhu cDC.
    cdc = row(primary, cohort="GSE213516", cell_type="cDC")
    assert cdc["effect"] == pytest.approx(0.1428616610214657, rel=1e-9)
    assert cdc["pvalue"] == pytest.approx(0.0350833624288907, rel=1e-9)
    assert cdc["q_bh"] == pytest.approx(0.6299204155120731, rel=1e-9)


def test_p1m_matched_null_floor_hits():
    models = load("p1m-models.csv")
    primary = models[(models["model"] == "primary") & (models["status"] == "fit")].copy()
    floor = 1 / 201
    assert (primary["random_pvalue"].sub(floor).abs() < 1e-9).sum() == 26
    ratio = (primary["effect"] - primary["random_mean_effect"]).abs() / primary["random_sd"]
    assert ratio.median() == pytest.approx(5.06, abs=0.01)
    assert ratio.min() == pytest.approx(0.16, abs=0.01)
    assert ratio.max() == pytest.approx(11.62, abs=0.01)
    # Empirical-only passes are 21 of 22 in Zhu and 8 of 13 in Thomson.
    for cohort, expected in (("GSE213516", 21), ("GSE214546", 8)):
        subset = primary[primary["cohort"] == cohort]
        empirical_only = ((subset["random_q_bh"] < 0.05) & (subset["q_bh"] >= 0.05)).sum()
        assert empirical_only == expected
    # Null SD far below the regression SE for Zhu cDC.
    cdc = row(primary, cohort="GSE213516", cell_type="cDC")
    assert cdc["se"] / cdc["random_sd"] == pytest.approx(3.63, abs=0.01)


def test_clr_composition_passes():
    models = load("clr-models.csv")
    primary = models[
        (models["model"] == "primary")
        & (models["endpoint"] == "CLR_composition")
        & (models["status"] == "fit")
    ]
    assert len(primary) == 35
    cd8 = row(primary, cohort="GSE214546", cell_type="CD8 Naive")
    assert cd8["effect"] == pytest.approx(-1.318945126196282, rel=1e-9)
    assert cd8["pvalue"] == pytest.approx(0.0004444880274658, rel=1e-9)
    assert cd8["q_bh"] == pytest.approx(0.0057783443570562, rel=1e-9)
    memory = row(primary, cohort="GSE214546", cell_type="CD4 Memory")
    assert memory["effect"] == pytest.approx(0.3347, abs=0.001)
    assert memory["q_bh"] == pytest.approx(0.0390604260672042, rel=1e-9)
    zhu = primary[primary["cohort"] == "GSE213516"]
    assert zhu["q_bh"].min() == pytest.approx(0.4546437084394136, rel=1e-9)
    zhu_cd8 = row(primary, cohort="GSE213516", cell_type="CD8-Naive")
    assert zhu_cd8["effect"] == pytest.approx(-0.9542303694165808, rel=1e-9)
    assert zhu_cd8["q_bh"] == pytest.approx(0.4792433744144554, rel=1e-9)
    assert zhu_cd8["mde80"] == pytest.approx(2.1772377392374995, rel=1e-9)
    # CD4 Memory passes only the primary model.
    thomson_memory = models[
        (models["cohort"] == "GSE214546")
        & (models["cell_type"] == "CD4 Memory")
        & (models["endpoint"] == "CLR_composition")
    ]
    assert row(thomson_memory, model="log_total_cells")["q_bh"] > 0.05
    assert row(thomson_memory, model="batch")["q_bh"] > 0.05
    assert row(thomson_memory, model="pool")["q_bh"] > 0.05
    raw = models[
        (models["cohort"] == "GSE214546")
        & (models["cell_type"] == "CD4 Memory")
        & (models["endpoint"] == "raw_fraction_sensitivity")
    ]
    assert row(raw, model="primary")["q_bh"] > 0.05
    # CD8 Naive fails batch but passes raw fraction.
    thomson_cd8 = models[
        (models["cohort"] == "GSE214546")
        & (models["cell_type"] == "CD8 Naive")
        & (models["endpoint"] == "CLR_composition")
    ]
    assert row(thomson_cd8, model="batch")["q_bh"] > 0.05
    raw_cd8 = models[
        (models["cohort"] == "GSE214546")
        & (models["cell_type"] == "CD8 Naive")
        & (models["endpoint"] == "raw_fraction_sensitivity")
    ]
    assert row(raw_cd8, model="primary")["q_bh"] == pytest.approx(0.0016170432431603, rel=1e-9)


def test_clr_targeted_family():
    targeted = load("clr-targeted.csv")
    assert len(targeted) == 6
    cd8 = row(targeted, cohort="GSE214546", cell_type="CD8 Naive")
    assert cd8["q_bh"] == pytest.approx(0.0008889760549317, rel=1e-9)
    mono = row(targeted, cohort="GSE214546", cell_type="CD16+ Monocytes")
    assert mono["q_bh"] == pytest.approx(0.0242254838234206, rel=1e-9)
    assert row(targeted, cohort="GSE213516", cell_type="CD8-Naive")["q_bh"] > 0.05
    mait = targeted[targeted["cell_type"] == "MAIT"]
    assert len(mait) == 2
    assert set(mait["status"]) == {"untestable"}
    assert set(mait["reason"]) == {"released_label_absent_or_not_harmonizable"}


def test_hellinger_negative_displacement_zhu():
    models = load("hellinger-models.csv")
    primary = models[(models["model"] == "primary") & (models["status"] == "fit")]
    assert len(primary) == 33
    zhu = primary[primary["cohort"] == "GSE213516"]
    assert len(zhu) == 20
    assert (zhu["q_bh"] < 0.05).sum() == 4
    expected = {
        "CD4-Tm": (-0.0126205096017646, 0.0256393901923262),
        "NK-GZMK": (-0.0136932843237572, 0.0256393901923262),
        "cDC": (-0.0065254372243776, 0.0314525390599707),
        "CD4-Naive": (-0.0099291450457932, 0.0314525390599707),
    }
    for cell_type, (effect, q) in expected.items():
        found = row(zhu, cell_type=cell_type)
        assert found["effect"] == pytest.approx(effect, rel=1e-9)
        assert found["q_bh"] == pytest.approx(q, rel=1e-9)
        assert found["effect"] < 0
    thomson = primary[primary["cohort"] == "GSE214546"]
    assert len(thomson) == 13
    assert (thomson["q_bh"] < 0.05).sum() == 0
    # CD8-Tem dies on precision, same direction as the original pass.
    cd8tem = row(zhu, cell_type="CD8-Tem")
    assert cd8tem["effect"] == pytest.approx(-0.0102066951155472, rel=1e-9)
    assert cd8tem["q_bh"] == pytest.approx(0.1030916958662582, rel=1e-9)
    # CD4-Treg is untestable at 30-plus.
    treg = models[(models["cell_type"] == "CD4-Treg") & (models["model"] == "primary")]
    assert len(treg) == 1
    assert treg.iloc[0]["status"] == "untestable"


def test_targeted_genes_zero_pass():
    genes = load("pseudobulk-targeted.csv")
    primary = genes[(genes["model"] == "primary") & (genes["status"] == "fit")]
    assert len(primary) == 12
    assert (primary["q_bh"] < 0.05).sum() == 0
    assert primary["q_bh"].min() == pytest.approx(0.0843783738076173, rel=1e-9)
    zhu = primary[primary["cohort"] == "GSE213516"]
    assert (zhu["effect"] > 0).sum() == 5
    assert row(zhu, gene="GZMA")["effect"] < 0
    thomson = primary[primary["cohort"] == "GSE214546"]
    assert (thomson["effect"] > 0).sum() == 5
    assert row(thomson, gene="PRF1")["effect"] > 0
    assert row(thomson, gene="GZMA")["effect"] < 0
    assert row(thomson, gene="CCL5")["effect"] == pytest.approx(1.1224, abs=0.001)
    assert row(zhu, gene="CCL5")["effect"] == pytest.approx(0.8560, abs=0.001)


def test_runtime_model_column_is_hardcoded():
    runtime = load("runtime.csv")
    assert set(runtime["model_used"]) == {"Gemini 4 Neon XHigh"}
    script = (REPO / "scripts" / "run_mature_pbmc.py").read_text()
    assert '"Gemini 4 Neon XHigh"' in script
