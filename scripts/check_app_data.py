#!/usr/bin/env python3
"""
check_app_data.py — Step 1 fidelity gate for the Agescape webapp exports.

Asserts every anchor in the app/DATA_CONTRACT.md audit table against the JSON
exports in app/data/ (plus the derived rules: design-grid totals, per-tissue
powered counts, file sizes, provenance checksums). Prints a PASS/FAIL table and
exits 0 only when every anchor holds. Run after every export:

    /home/reza/miniconda3/bin/conda run -n agescape python scripts/check_app_data.py
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "app" / "data"
RESULTS = REPO / "results"

results: list[tuple[str, str, str, bool]] = []


def check(anchor: str, expected, observed) -> bool:
    ok = expected == observed
    results.append((anchor, repr(expected), repr(observed), ok))
    return ok


def f(x):
    return float(x)


def fchk(anchor: str, expected: float, observed: float, ndigits: int) -> bool:
    """Float anchor at the contract's stated precision (files store more)."""
    try:
        return check(anchor, expected, round(f(observed), ndigits))
    except (TypeError, ValueError):
        return check(anchor, expected, repr(observed))


def main() -> int:
    j = {p.stem: json.loads(p.read_text()) for p in sorted(DATA.glob("*.json"))}
    for name in ["atlas", "calibration", "agreement", "design-grid", "program",
                 "signature", "savina", "power-grid", "provenance"]:
        check(f"export exists: {name}.json", True, name in j)

    # ---- atlas: 163 units / 17 powered, confined to Lung & Kidney
    atlas = j["atlas"]
    check("atlas units == 163", 163, len(atlas["units"]))
    check("atlas ranked == 17", 17, len(atlas["ranked"]))
    powered = [u for u in atlas["units"] if u["powered"]]
    check("atlas powered == 17", 17, len(powered))
    check("powered confined to Lung & Kidney", {"Lung", "Kidney"},
          {u["tissue"] for u in powered})
    per_tissue = {}
    for u in powered:
        per_tissue[u["tissue"]] = per_tissue.get(u["tissue"], 0) + 1
    check("per-tissue powered {Kidney:8, Lung:9}", {"Kidney": 8, "Lung": 9}, per_tissue)
    check("sum of per-tissue powered == 17", 17, sum(per_tissue.values()))
    check("every unit carries all schema-1 fields", True,
          all(set(u) == {"unit", "tissue", "cell_ontology_class", "eligible", "powered",
                         "n_cells_3m", "n_cells_18m", "n_mice_3m", "n_mice_18m"}
              for u in atlas["units"]))

    # ---- ranking head 55/53/36/35 in order
    head = atlas["ranked"][:4]
    check("ranking head counts 55/53/36/35", [55, 53, 36, 35],
          [r["n_de_genes"] for r in head])
    check("ranking head units", [
        "Lung::adventitial cell",
        "Kidney::kidney proximal convoluted tubule epithelial cell",
        "Lung::fibroblast of lung",
        "Lung::alveolar macrophage",
    ], [r["unit"] for r in head])

    # ---- calibration 55 / 4.6 / 0
    cal = j["calibration"]
    check("calibration observed 55", 55, cal["observed"])
    check("calibration perm_mean 4.6", 4.6, cal["perm_mean"])
    check("calibration sex_floor 0", 0, cal["sex_floor"])
    check("calibration unit", "Lung::adventitial cell", cal["unit"])
    check("positive_controls == 4 rows", 4, len(cal["positive_controls"]))

    # ---- agreement 13/4/0 + rho 0.866 (n 13,019) + recompute over points
    agr = j["agreement"]
    check("agreement counts 13/4/0", {"agree": 13, "partial": 4, "disagree": 0}, agr["counts"])
    check("lollipop == 17 rows", 17, len(agr["lollipop"]))
    sc = agr["scatter"]
    check("scatter unit", "Kidney::kidney proximal convoluted tubule epithelial cell", sc["unit"])
    check("scatter n == 13019", 13019, sc["n"])
    check("scatter points == 13019", 13019, len(sc["points"]))
    rho_re = spearmanr([p["ours_log2fc"] for p in sc["points"]],
                       [p["paper_coef"] for p in sc["points"]]).statistic
    check("scatter rho field == 0.866", 0.866, sc["rho"])
    check("recomputed rho rounds to 0.866", 0.866, round(f(rho_re), 3))
    check("|recomputed rho - 0.866| < 0.01", True, abs(f(rho_re) - 0.866) < 0.01)
    check("drawn rho matches recompute (full precision)", True,
          abs(f(sc["rho_full"]) - f(rho_re)) < 1e-9)

    # ---- program 167 / 94 / 312 + matrix rules
    prog = j["program"]
    check("program summary up 167 / down 94 / entries 312",
          {"up": 167, "down": 94, "entries": 312}, prog["summary"])
    check("program genes == 312 entries", 312, len(prog["genes"]))
    check("program units == 17", 17, len(prog["units"]))
    pg_src = pd.read_csv(RESULTS / "m6-program-genes.csv")
    check("summary up == m6-program-genes program_up non-null", 167,
          int(pg_src["program_up"].notna().sum()))
    check("summary down == m6-program-genes program_down non-null", 94,
          int(pg_src["program_down"].notna().sum()))
    m40 = prog["matrix_top40"]
    check("matrix_top40 rows == 40", 40, len(m40["rows"]))
    check("matrix_top40 columns == 17 units", 17, len(m40["columns"]))
    check("matrix rows ranked by units-present (desc)",
          True, all(r["n_units"] >= s["n_units"] for r, s in zip(m40["rows"], m40["rows"][1:])))

    # ---- signature: every m6 slope/CI/p, TOST, permutation, coverage, lifespans
    sig = j["signature"]
    arms = sig["arms"]

    def arm_check(name, slope, ci, p, n, slp_d=5, ci_d=(5, 5), p_d=4):
        a = arms[name]
        fchk(f"{name} slope", slope, a["slope"], slp_d)
        check(f"{name} ci", [ci[0], ci[1]],
              [round(f(a["ci"][0]), ci_d[0]), round(f(a["ci"][1]), ci_d[1])])
        fchk(f"{name} p", p, a["p"], p_d)
        check(f"{name} n", n, a["n"])

    arm_check("mouse_wb", 0.23496, [0.06727, 0.40265], 0.0083, 23)
    arm_check("mouse_skin", 0.31947, [-3.44814, 4.08708], 0.4763, 3)
    check("mouse_skin underpowered", True, arms["mouse_skin"].get("underpowered") is True)
    arm_check("human_skin", -0.21297, [-1.37366, 0.94772], 0.5126, 4)
    check("human_skin low_powered", True, arms["human_skin"].get("low_powered") is True)
    arm_check("human_fibroblast", 0.07099, [-0.82588, 0.96787], 0.7659, 4)
    arm_check("nmr", 0.01667, [-0.25667, 0.29], 0.8585, 5, ci_d=(5, 3), p_d=4)
    check("nmr tested_span_years [1.1, 11.3]", [1.1, 11.3], arms["nmr"]["tested_span_years"])
    check("mouse_wb points == 23", 23, len(arms["mouse_wb"]["points"]))

    diff = sig["stats"]["diff"]
    fchk("slope-diff 0.21828", 0.21828, diff["slope"], 5)
    check("slope-diff ci", [-0.01263, 0.4492],
          [round(f(diff["ci"][0]), 5), round(f(diff["ci"][1]), 4)])
    fchk("slope-diff p 0.0639", 0.0639, diff["p"], 4)
    tost = sig["stats"]["tost"]
    fchk("TOST delta 0.05874", 0.05874, tost["delta"], 5)
    fchk("TOST p 0.32894", 0.32894, tost["p"], 5)
    check("TOST established == false", False, tost["established"])
    fchk("permutation p 0.81818", 0.81818, sig["stats"]["perm_p"], 5)

    cov = sig["coverage"]
    check("coverage mouse 1.0", 1.0, cov["mouse"])
    check("coverage human 0.73180", 0.73180, round(cov["human"], 5))
    check("coverage nmr 0.31414", 0.31414, round(cov["nmr"], 5))
    gates_src = pd.read_csv(RESULTS / "m6-gates.csv")
    gv = {r["gate"]: r["value"] for _, r in gates_src.iterrows()}
    check("coverage values map to the right species (vs m6-gates.csv)", True,
          cov["mouse"] == f(gv["coverage_mouse"])
          and round(cov["human"], 6) == round(f(gv["coverage_human"]), 6)
          and round(cov["nmr"], 6) == round(f(gv["coverage_nmr"]), 6))
    check("lifespans mouse 4.0 / nmr 31.0 / human 122.5",
          {"mouse": 4.0, "nmr": 31.0, "human": 122.5}, sig["lifespans"])

    # ---- power grid anchors on the correct tiles
    pg = j["power-grid"]
    check("power donors axis", [4, 6, 8, 10, 12, 16, 20, 30, 50], pg["donors"])
    check("power sds axis", [8, 12, 16, 20], pg["sds"])
    check("power dispersions axis", [0.1, 0.2, 0.4], pg["dispersions"])
    check("power tiles == 108", 108, len(pg["tiles"]))
    tile = {(t["n_donors"], t["sd"], t["dispersion"]): t["min_30y_log2fc"] for t in pg["tiles"]}
    check("anchor 1.474 @ (4 donors, SD 8, disp 0.1)", 1.474, tile[(4, 8, 0.1)])
    check("anchor 0.59 @ (4 donors, SD 20, disp 0.1)", 0.59, tile[(4, 20, 0.1)])
    check("anchors dict mirrors tiles", {"4_8_low": 1.474, "4_20_low": 0.59}, pg["anchors"])

    # ---- savina markers
    sav = j["savina"]
    check("savina curve == 39 rows", 39, len(sav["curve"]))
    mrk = sav["markers"]
    check("savina paper_rule 29 @ 1.0 (p<0.05)", {"cutoff": 1.0, "n_pval": 29}, mrk["paper_rule"])
    check("savina converted_rule 9 & 6 @ 1.45",
          {"cutoff": 1.45, "n_pval": 9, "n_fdr": 6}, mrk["converted_rule"])
    check("savina loose_rule 560 @ 0.25 (FDR)", {"cutoff": 0.25, "n_fdr": 560}, mrk["loose_rule"])

    # ---- design grid totals match m1-sex-age-table.csv (+ hatch rules)
    dg = j["design-grid"]
    src = pd.read_csv(RESULTS / "m1-sex-age-table.csv")
    allrow = src[src["age_group"] == "All"].iloc[0]
    check("design-grid age_groups", ["1m", "3m", "18m", "21m", "24m", "30m"], dg["age_groups"])
    check("cells female total == m1 All row", int(allrow["cells_female"]), sum(dg["cells"]["female"]))
    check("cells male total == m1 All row", int(allrow["cells_male"]), sum(dg["cells"]["male"]))
    check("mice female total == m1 All row", int(allrow["mice_female"]), sum(dg["mice"]["female"]))
    check("mice male total == m1 All row", int(allrow["mice_male"]), sum(dg["mice"]["male"]))
    d = src[src["age_group"] != "All"].reset_index(drop=True)
    check("cells female per group == source", list(d["cells_female"]), dg["cells"]["female"])
    check("cells male per group == source", list(d["cells_male"]), dg["cells"]["male"])
    fem = [g for g, v in zip(dg["age_groups"], dg["cells"]["female"]) if v > 0 and
           dg["cells"]["male"][dg["age_groups"].index(g)] == 0]
    mal = [g for g, v in zip(dg["age_groups"], dg["cells"]["male"]) if v > 0 and
           dg["cells"]["female"][dg["age_groups"].index(g)] == 0]
    both = [g for g in dg["age_groups"] if dg["cells"]["female"][dg["age_groups"].index(g)] > 0
            and dg["cells"]["male"][dg["age_groups"].index(g)] > 0]
    check("female-only columns {21m}", ["21m"], fem)
    check("male-only columns {1m,24m,30m}", ["1m", "24m", "30m"], mal)
    check("both-sex columns {3m,18m}", ["3m", "18m"], both)

    # ---- provenance: freeze date, checksums match files on disk, sizes
    pv = j["provenance"]
    check("freeze_date 2026-10-14", "2026-10-14", pv["freeze_date"])
    check("export_date present", True, bool(pv.get("export_date")))
    ok_sha, n_sources = True, 0
    for export, entries in pv["sources"].items():
        for e in entries:
            n_sources += 1
            actual = hashlib.sha256((REPO / e["file"]).read_bytes()).hexdigest()
            if actual != e["sha256"]:
                ok_sha = False
    check("all provenance SHA-256 values match files on disk", True, ok_sha)
    check("every export has provenance sources", 9, len(pv["sources"]))
    check("figure assets recorded == 4", 4, len(pv["assets"]))
    fig_ok = all((REPO / "app" / "figures" / name).stat().st_size == meta["bytes"]
                 for name, meta in pv["assets"].items())
    check("figure asset bytes match app/figures", True, fig_ok)

    # ---- size budget
    sizes = {p.name: p.stat().st_size for p in DATA.glob("*.json")}
    check("every file < 1 MB", True, all(sz < 1_000_000 for sz in sizes.values()))
    check("total app/data < 5 MB", True, sum(sizes.values()) < 5_000_000)

    # ---- report
    width_a = max(len(r[0]) for r in results)
    print(f"{'ANCHOR'.ljust(width_a)}  {'EXPECTED'.ljust(52)} {'OBSERVED'.ljust(52)} STATUS")
    print("-" * (width_a + 112))
    n_fail = 0
    for anchor, exp, obs, ok in results:
        print(f"{anchor.ljust(width_a)}  {exp.ljust(52)} {obs.ljust(52)} {'PASS' if ok else 'FAIL'}")
        if not ok:
            n_fail += 1
    print("-" * (width_a + 112))
    if n_fail:
        print(f"FAIL: {n_fail} of {len(results)} anchors missed:")
        for anchor, exp, obs, ok in results:
            if not ok:
                print(f"  - {anchor}: expected {exp}, observed {obs}")
        return 1
    print(f"PASS: all {len(results)} anchors hold.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
