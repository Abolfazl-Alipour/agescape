import json
import subprocess
import sys
from pathlib import Path
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent


def test_anchor_lock_against_committed_results():
    lock_file = REPO_ROOT / "tests" / "fixtures" / "ANCHOR_LOCK.json"
    if not lock_file.exists():
        lock_file = REPO_ROOT / "scratch" / "ANCHOR_LOCK.json"
    assert lock_file.exists(), "ANCHOR_LOCK.json is missing!"
    lock = json.loads(lock_file.read_text())

    # 1. M2 positive controls
    pctrl = pd.read_csv(REPO_ROOT / "results" / "m2-calibration-positive-controls.csv").set_index("gene")
    for gene, vals in lock["m2_positive_controls"].items():
        assert round(float(pctrl.loc[gene, "log2FC"]), 3) == vals["log2fc"]
        assert round(float(pctrl.loc[gene, "pvalue"]), 3) == vals["pvalue"]

    # 2. M4 power
    pow_df = pd.read_csv(REPO_ROOT / "results" / "m4-human-power.csv")
    for _, row in pow_df.iterrows():
        n = int(row["n_donors"])
        sd = int(row["donor_age_sd_years"])
        disp = row["biological_dispersion"]
        k = f"n{n}_sd{sd}_disp{disp}"
        if k in lock["m4_power"]:
            assert round(float(row["min_detectable_30y_log2fc"]), 3) == lock["m4_power"][k]

    # 3. M5 hierarchy
    hier = pd.read_csv(REPO_ROOT / "results" / "m5-program-hierarchy.csv")
    counts = hier["tier"].value_counts().to_dict()
    assert counts == lock["m5_hierarchy"]

    # 4. M5 overlaps
    ov = pd.read_csv(REPO_ROOT / "results" / "m5-external-program-overlap.csv")
    p1_ov = ov[ov["tier"] == "P1_Clean"].set_index("external_set")
    for s, vals in lock["m5_external_overlaps"].items():
        assert int(p1_ov.loc[s, "n_shared"]) == vals["shared"]
        assert round(float(p1_ov.loc[s, "fold_enrichment"]), 2) == vals["fold_enrichment"]

    # 5. M6 stats
    m6_stat = pd.read_csv(REPO_ROOT / "results" / "m6-statistical-summary.csv").set_index("specification")
    for spec, vals in lock["m6_stats"].items():
        if spec in m6_stat.index:
            obs_slope = round(float(m6_stat.loc[spec, "slope"]), 4) if pd.notna(m6_stat.loc[spec, "slope"]) else None
            assert obs_slope == vals["slope"], f"Slope mismatch for {spec}: {obs_slope} vs {vals['slope']}"

    # 6. M9 cohorts
    m9 = pd.read_csv(REPO_ROOT / "results" / "m9-cohort-summary.csv")
    assert int(m9["n_donors"].sum()) == lock["m9_cohorts"]["total_donors"]

    # 7. M10 lifespan
    m10_wc = pd.read_csv(REPO_ROOT / "results" / "m10-cohort-summary.csv").iloc[0]
    assert int(m10_wc["n_donors"]) == lock["m10_lifespan"]["whole_cortex"]["n_donors"]
    assert round(float(m10_wc["slope"]), 5) == lock["m10_lifespan"]["whole_cortex"]["slope"]


def test_app_data_check_script():
    cmd = [sys.executable, "scripts/check_app_data.py"]
    res = subprocess.run(cmd, cwd=REPO_ROOT, capture_output=True, text=True)
    assert res.returncode == 0, f"check_app_data.py failed:\n{res.stdout}\n{res.stderr}"
    assert "PASS: all 146 anchors hold." in res.stdout


def test_app_data_sizes_and_count():
    data_dir = REPO_ROOT / "app" / "data"
    json_files = list(data_dir.glob("*.json"))
    assert len(json_files) >= 12, f"Expected at least 12 json files, got {len(json_files)}"
    total_size = sum(f.stat().st_size for f in json_files)
    assert total_size < 5_000_000, f"Total app data {total_size} exceeds 5MB limit"
    for f in json_files:
        assert f.stat().st_size < 1_000_000, f"File {f} exceeds 1MB limit ({f.stat().st_size} bytes)"
