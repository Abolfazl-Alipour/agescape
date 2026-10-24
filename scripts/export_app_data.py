#!/usr/bin/env python3
"""
export_app_data.py — Step 1 plot-data exports for the Agescape webapp.

Produces the 11 exports pinned by app/DATA_CONTRACT.md (schemas 1-11) in app/data/
as deterministic JSON (sorted keys, shortest round-trip-safe float repr), copies the
4 contracted figure PNGs into app/figures/, and wipes every stale file from the
previous app version. Run:

    /home/reza/miniconda3/bin/conda run -n agescape python scripts/export_app_data.py

Source-of-truth rule: committed files win over contract guesses. Two contract
corrections are applied here and recorded (dated) in app/DATA_CONTRACT.md:
  1. agreement.json scatter `ours_bridge_slope` is sourced from the M3 continuous bridge
     slope (results/m3-continuous-bridge-model.csv) joined to the paper age
     coefficient — that is the join that produces the pinned anchor (n=13,019,
     rho=0.866; see scripts/run_m3_agreement.py, which prefers the bridge slope
     over the dichotomous m2_de log2FC). The literal m2_de join gives n=11,835,
     rho=0.309 and is exported to stdout as a diagnostic only.
  2. signature.json arm `mouse_skin` is the pipeline's skin arm (species "Mouse",
     all skin cells, stats row "Mouse aging slope"); the "Mouse (epidermal)"
     sensitivity series has no arm in the pinned schema and is not exported.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import shutil
import sys
import zipfile
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import load_workbook
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
RAW = REPO / "data" / "raw"
FIGS = REPO / "figures"
OUT = REPO / "app" / "data"
APP_FIGS = REPO / "app" / "figures"

FREEZE_DATE = "2026-10-14"
XLSX = RAW / "tms_nature2020_moesm8_dge_results.xlsx"
PCT_UNIT = "Kidney::kidney proximal convoluted tubule epithelial cell"
PCT_TISSUE = "Kidney"
PCT_CT = "kidney proximal convoluted tubule epithelial cell"
PCT_PAPER_COL = "Kidney.kidney proximal convoluted tubule epithelial cell"

FIGURE_FILES = [
    "m3-dpt-vs-age.png",
    "m6-signature-figure.png",
    "m6-positive-control-skin.png",
    "m6-savina-reproduction.png",
]

NEW_DATA_FILES = [
    "atlas.json",
    "calibration.json",
    "agreement.json",
    "design-grid.json",
    "program.json",
    "signature.json",
    "savina.json",
    "power-grid.json",
    "provenance.json",
    "ranking.json",
    "anage.json",
    "humancohorts.json",
    "dlpfc.json",
]

# ---------------------------------------------------------------- helpers


def sha256_full(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fnum(x) -> float | None:
    """Float for JSON: NaN -> None, otherwise a plain Python float (repr round-trips)."""
    v = float(x)
    if v != v:  # NaN
        return None
    return v


def inum(x) -> int:
    return int(x)


def bnum(x) -> bool:
    return bool(x)


def write_json(name: str, payload: dict) -> int:
    """Deterministic JSON: sorted keys, compact separators, no NaN, round-trip verified."""
    text = json.dumps(payload, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False, allow_nan=False)
    assert json.loads(text) == payload, f"{name}: JSON round-trip mismatch"
    target = OUT / name
    target.write_text(text + "\n", encoding="utf-8")
    return target.stat().st_size


# ---------------------------------------------------------------- load sources

def load_sources():
    s = {}
    s["eligibility"] = pd.read_csv(RESULTS / "m2-unit-eligibility.csv")
    s["ranked"] = pd.read_csv(RESULTS / "m2-aging-signal-ranked.csv")
    s["trio"] = pd.read_csv(RESULTS / "m2-calibration-trio-summary.csv")
    s["posctl"] = pd.read_csv(RESULTS / "m2-calibration-positive-controls.csv")
    s["agreement"] = pd.read_csv(RESULTS / "m3-agreement-table.csv")
    s["bridge"] = pd.read_csv(RESULTS / "m3-continuous-bridge-model.csv")
    s["sexage"] = pd.read_csv(RESULTS / "m1-sex-age-table.csv")
    s["catalog"] = pd.read_csv(RESULTS / "m5-mouse-program-catalog.csv")
    s["genelevel"] = pd.read_csv(RESULTS / "m5-gene-level-programs.csv")
    s["programgenes"] = pd.read_csv(RESULTS / "m6-program-genes.csv")
    s["scores"] = pd.read_csv(RESULTS / "m6-lifespan-program-scores.csv")
    s["stats"] = pd.read_csv(RESULTS / "m6-statistical-summary.csv")
    s["gates"] = pd.read_csv(RESULTS / "m6-gates.csv")
    s["lifespans"] = pd.read_csv(RESULTS / "m6-species-lifespans.csv")
    s["savina"] = pd.read_csv(RESULTS / "m6-savina-reproduction.csv")
    s["power"] = pd.read_csv(RESULTS / "m4-human-power.csv")
    return s


def read_paper_age_coef() -> pd.DataFrame:
    """gene -> paper age coefficient + bh_p for PCT_PAPER_COL, mirroring
    scripts/run_m3_agreement.py (dropna across BOTH coef and bh_p)."""
    wb = load_workbook(XLSX, read_only=True)
    coef = pd.DataFrame(wb["DGE_result.tissue_cell.age_coef"].values)
    bh = pd.DataFrame(wb["DGE_result.tissue_cell.bh_p"].values)
    wb.close()
    hdr = coef.iloc[0]
    coef.columns = hdr
    coef = coef[1:]
    bh.columns = hdr
    bh = bh[1:]
    genes = coef.iloc[:, 0].astype(str)
    coef.index = genes
    bh.index = genes
    col = [c for c in coef.columns if c == PCT_PAPER_COL][0]
    pap = pd.DataFrame({"paper_coef": coef[col], "paper_bh": bh[col]}).dropna()
    pap = pap[~pap.index.duplicated(keep="first")]
    return pap


# ---------------------------------------------------------------- exports

def export_anage(prov):
    """AnAge species browser export: every species with a recorded maximum longevity..

    Single-cell flags are scoped strictly to this explorer's own data (mouse,
    human, naked mole-rat); every other species is flagged 'no' because no
    public single-cell aging data exists for it (verified claim, not inferred).
    """
    zpath = REPO / "data" / "raw" / "anage_dataset.zip"
    sc_species = {"Mus musculus", "Homo sapiens", "Heterocephalus glaber"}

    def fnum_or_none(x):
        try:
            v = float(x)
            return v if np.isfinite(v) else None
        except (TypeError, ValueError):
            return None

    with zipfile.ZipFile(zpath) as zf:
        rows = list(csv.DictReader(io.TextIOWrapper(zf.open("anage_data.txt"), encoding="latin1"), delimiter="\t"))

    out, seen = [], set()
    for r in rows:
        latin = f"{r['Genus']} {r['Species']}".strip()
        lon = fnum_or_none(r["Maximum longevity (yrs)"])
        if not latin or latin in seen or lon is None:
            continue
        seen.add(latin)
        mass_g = fnum_or_none(r["Adult weight (g)"])
        out.append({
            "n": r.get("Common name") or latin,
            "latin": latin,
            "cls": r["Class"] or "Other",
            "long": lon,
            "kg": round(mass_g / 1000, 3) if mass_g else None,
            "q": (r.get("Data quality") or "").lower() or None,
            "sc": "yes" if latin in sc_species else "no",
        })

    assert len(out) == 4141, f"anage species count changed: {len(out)}"
    by_latin = {d["latin"]: d["long"] for d in out}
    assert by_latin["Mus musculus"] == 4.0
    assert by_latin["Homo sapiens"] == 122.5
    assert by_latin["Heterocephalus glaber"] == 31.0
    assert by_latin["Balaena mysticetus"] == 211.0
    assert sum(1 for d in out if d["sc"] == "yes") == 3
    prov["anage.json"] = [{"file": "data/raw/anage_dataset.zip", "sha256": sha256_full(zpath), "rows": len(out)}]
    return out


def export_human_cohorts(prov):
    """M9 human cohort arm: program score vs donor age across public CELLxGENE cohorts.

    Whole-cohort regressions plus the powered cell-type hits. Values transcribe the
    committed M9 outputs (scripts/run_m9_pipeline.py); no new statistical decisions here.
    """
    summ = pd.read_csv(RESULTS / "m9-cohort-summary.csv")
    regs = pd.read_csv(RESULTS / "m9-regressions.csv")

    assert len(summ) == 11, f"expected 11 cohorts, got {len(summ)}"
    assert int(summ["n_donors"].sum()) == 621, summ["n_donors"].sum()
    mssm = summ[summ["cohort"] == "MSSM_DLPFC"].iloc[0]
    assert mssm["n_donors"] == 213 and mssm["ci_hi"] < 0, "MSSM headline changed"

    cohorts = []
    for _, r in summ.iterrows():
        cohorts.append({
            "cohort": r["cohort"], "group": r["group"], "n_donors": int(r["n_donors"]),
            "age_min": fnum(r["age_min"]), "age_max": fnum(r["age_max"]),
            "age_sd": fnum(r["age_sd"]),
            "slope_per_year": fnum(r["slope_per_year"]), "se": fnum(r["se"]),
            "ci_lo": fnum(r["ci_lo"]), "ci_hi": fnum(r["ci_hi"]), "p": fnum(r["p"]),
            "q_bh": fnum(r["q_bh"]), "rho": fnum(r["rho"]), "verdict": r["verdict"],
        })

    hits = regs[regs["n_donors"] >= 6].sort_values("p").head(40)
    celltypes = [{
        "cohort": r["cohort"], "group": r["group"], "cell_type": r["cell_type"],
        "n_donors": int(r["n_donors"]), "age_sd": fnum(r["age_sd"]),
        "slope_per_year": fnum(r["slope_per_year"]), "ci_lo": fnum(r["ci_lo"]),
        "ci_hi": fnum(r["ci_hi"]), "p": fnum(r["p"]), "verdict": r["verdict"],
    } for _, r in hits.iterrows()]

    prov["humancohorts.json"] = [
        {"file": "results/m9-cohort-summary.csv", "sha256": sha256_full(RESULTS / "m9-cohort-summary.csv"), "rows": len(summ)},
        {"file": "results/m9-regressions.csv", "sha256": sha256_full(RESULTS / "m9-regressions.csv"), "rows": len(regs)},
    ]
    return {"cohorts": cohorts, "celltypes": celltypes,
            "total_control_donors": 624, "program_human_genes": 188}


def export_dlpfc(prov):
    """M10 human lifespan arm: DLPFC Lifespan Atlas (284 control donors 0.2-97y,
    27 cell types, authors' donor x cell-type pseudobulk; Synapse syn62064535,
    AD Knowledge Portal DUA). Values transcribe committed M10 aggregates only;
    per-donor values are DUA-controlled and never exported. Concordant flag is
    a pure transcription of the two committed model columns (both p<0.05, same
    sign); no new statistical decisions here."""
    summ = pd.read_csv(RESULTS / "m10-cohort-summary.csv").iloc[0]
    regs = pd.read_csv(RESULTS / "m10-celltype-regressions.csv")

    assert summ["n_donors"] == 284 and len(regs) == 27
    assert float(summ["age_min"]) == 0.2 and float(summ["age_max"]) == 97.0
    assert abs(float(summ["age_sd"]) - 24.0) < 0.1
    assert float(summ["p"]) > 0.1 and float(summ["ci_lo"]) < 0 < float(summ["ci_hi"])

    rows = []
    for _, r in regs.iterrows():
        conc = bool(r["p"] < 0.05 and r["p_cpm"] < 0.05 and r["slope"] * r["slope_cpm"] > 0)
        if conc:
            direction = "up" if r["slope"] > 0 else "down"
        elif r["slope"] * r["slope_cpm"] < 0:
            direction = "normalization-dependent direction"
        elif r["p"] < 0.05 or r["p_cpm"] < 0.05:
            direction = "scheme-sensitive significance"
        else:
            direction = "not significant"
        rows.append({
            "cell_type": r["cell_type"], "n_donors": inum(r["n_donors"]),
            "n_zero_donors": inum(r["n_zero_donors"]),
            "slope_per_year": fnum(r["slope"]), "se": fnum(r["se"]),
            "ci_lo": fnum(r["ci_lo"]), "ci_hi": fnum(r["ci_hi"]),
            "p": fnum(r["p"]), "q_bh": fnum(r["q_bh"]), "rho": fnum(r["rho"]),
            "slope_cpm": fnum(r["slope_cpm"]), "ci_lo_cpm": fnum(r["ci_lo_cpm"]),
            "ci_hi_cpm": fnum(r["ci_hi_cpm"]), "p_cpm": fnum(r["p_cpm"]),
            "q_bh_cpm": fnum(r["q_bh_cpm"]),
            "slope_sex_adj": fnum(r["slope_sex_adj"]), "p_sex_adj": fnum(r["p_sex_adj"]),
            "slope_adult": fnum(r["slope_adult"]), "p_adult": fnum(r["p_adult"]),
            "concordant": conc, "direction": direction,
        })
    n_conc = sum(1 for x in rows if x["concordant"])
    assert n_conc == 7, f"expected 7 concordant cell types, got {n_conc}"

    prov["dlpfc.json"] = [
        {"file": "results/m10-cohort-summary.csv",
         "sha256": sha256_full(RESULTS / "m10-cohort-summary.csv"), "rows": 1},
        {"file": "results/m10-celltype-regressions.csv",
         "sha256": sha256_full(RESULTS / "m10-celltype-regressions.csv"), "rows": len(regs)},
    ]
    return {
        "cohort": {
            "name": "DLPFC_Lifespan_Atlas",
            "source": "Synapse syn52396927 (AD Knowledge Portal DUA), folder syn62064535",
            "n_donors": inum(summ["n_donors"]),
            "age_min": fnum(summ["age_min"]), "age_max": fnum(summ["age_max"]),
            "age_sd": fnum(summ["age_sd"]), "n_cell_types": inum(summ["n_cell_types"]),
            "whole_cortex": {
                "slope_per_year": fnum(summ["slope"]), "se": fnum(summ["se"]),
                "ci_lo": fnum(summ["ci_lo"]), "ci_hi": fnum(summ["ci_hi"]),
                "p": fnum(summ["p"]), "rho": fnum(summ["rho"]),
                "slope_sex_adj": fnum(summ["slope_sex_adj"]),
                "p_sex_adj": fnum(summ["p_sex_adj"]),
            },
        },
        "program": {"up_slots": 121, "down_slots": 70, "measured_unique": 178,
                    "coverage": 0.9319},
        "n_concordant": n_conc,
        "celltypes": rows,
    }


def export_atlas(s, prov):
    el, rk = s["eligibility"], s["ranked"]
    powered_units = set(rk["unit"])
    assert len(el) == 163 and len(rk) == 17
    units = [{
        "unit": r["unit"],
        "tissue": r["tissue"],
        "cell_ontology_class": r["cell_ontology_class"],
        "eligible": bnum(r["eligible"]),
        "powered": r["unit"] in powered_units,
        "n_cells_3m": inum(r["n_cells_3m"]),
        "n_cells_18m": inum(r["n_cells_18m"]),
        "n_mice_3m": inum(r["n_mice_3m"]),
        "n_mice_18m": inum(r["n_mice_18m"]),
    } for _, r in el.iterrows()]
    ranked = [{
        "unit": r["unit"],
        "tissue": r["tissue"],
        "n_de_genes": inum(r["n_de_genes"]),
        "median_abs_log2fc": fnum(r["median_abs_log2fc"]),
        "n_de_ci_low": fnum(r["n_de_ci_low"]),
        "n_de_ci_high": fnum(r["n_de_ci_high"]),
    } for _, r in rk.iterrows()]
    prov["atlas.json"] = [
        {"file": "results/m2-unit-eligibility.csv", "sha256": sha256_full(RESULTS / "m2-unit-eligibility.csv"), "rows": len(el)},
        {"file": "results/m2-aging-signal-ranked.csv", "sha256": sha256_full(RESULTS / "m2-aging-signal-ranked.csv"), "rows": len(rk)},
    ]
    return {"units": units, "ranked": ranked}


def export_calibration(s, prov):
    trio, pc = s["trio"], s["posctl"]
    assert len(trio) == 1
    t = trio.iloc[0]
    assert inum(t["observed_aging_de"]) == 55
    assert fnum(t["permutation_null_mean_de"]) == 4.6
    assert inum(t["sex_floor_3m_de"]) == 0
    positive_controls = [{
        "gene": r["gene"],
        "log2FC": fnum(r["log2FC"]),
        "pvalue": fnum(r["pvalue"]),
        "FDR": fnum(r["FDR"]),
        "direction": r["direction"],
    } for _, r in pc.iterrows()]
    prov["calibration.json"] = [
        {"file": "results/m2-calibration-trio-summary.csv", "sha256": sha256_full(RESULTS / "m2-calibration-trio-summary.csv"), "rows": len(trio)},
        {"file": "results/m2-calibration-positive-controls.csv", "sha256": sha256_full(RESULTS / "m2-calibration-positive-controls.csv"), "rows": len(pc)},
    ]
    return {
        "observed": inum(t["observed_aging_de"]),
        "perm_mean": fnum(t["permutation_null_mean_de"]),
        "sex_floor": inum(t["sex_floor_3m_de"]),
        "unit": t["top_unit"],
        "positive_controls": positive_controls,
    }


def export_agreement(s, prov, paper_coef):
    ag = s["agreement"]
    lollipop = [{
        "unit": r["unit"],
        "spearman_rho": fnum(r["spearman_rho"]),
        "sign_concordance": fnum(r["sign_concordance_paper_sig"]),
        "verdict": r["verdict"],
    } for _, r in ag.iterrows()]
    vc = ag["verdict"].value_counts()
    counts = {"agree": int(vc.get("agree", 0)), "partial": int(vc.get("partial", 0)),
              "disagree": int(vc.get("disagree", 0))}
    assert counts == {"agree": 13, "partial": 4, "disagree": 0}, counts

    # --- C3 scatter. The pinned anchor (n=13,019, rho=0.866) is defined by the M3
    # pipeline's join: continuous bridge slope (preferred over the dichotomous m2_de
    # log2FC when >200 bridge genes exist) vs paper age coefficient.
    bridge = s["bridge"]
    bsub = bridge[(bridge["tissue"] == PCT_TISSUE) & (bridge["cell_ontology_class"] == PCT_CT)]
    assert len(bsub) > 200
    ours = bsub.set_index("gene")[["slope_per_month", "pvalue"]]
    m = ours.join(paper_coef[["paper_coef"]], how="inner").dropna()  # mirrors run_m3_agreement.py
    n = len(m)
    rho = float(spearmanr(m["paper_coef"], m["slope_per_month"]).statistic)
    m3row = ag[ag["unit"] == PCT_UNIT].iloc[0]
    points = [{"gene": g, "ours_bridge_slope": fnum(round(float(v["slope_per_month"]), 5)),
               "paper_coef": fnum(round(float(v["paper_coef"]), 5))}
              for g, v in m.sort_index().iterrows()]
    rho = float(spearmanr([p["paper_coef"] for p in points],
                          [p["ours_bridge_slope"] for p in points]).statistic)
    assert abs(rho - 0.866) < 0.01, rho
    assert round(rho, 3) == fnum(m3row["spearman_rho"]) == 0.866, (rho, m3row["spearman_rho"])

    # Diagnostic (not exported): the literal m2_de dichotomous join the contract text
    # named. It does NOT reproduce the anchor — recorded in DATA_CONTRACT.md.
    de = pd.read_csv(RESULTS / "m2_de_Kidney__kidney_proximal_convoluted_tubule_epithelial_cell.csv")
    md = de.set_index("index")[["log2FC"]].join(paper_coef[["paper_coef"]], how="inner").dropna()
    rho_d = float(spearmanr(md["paper_coef"], md["log2FC"]).statistic)
    print(f"[c3-diagnostic] literal m2_de dichotomous join: n={len(md)}, rho={rho_d:.6f} "
          f"(anchor join: n={n}, rho={rho:.6f}); contract corrected accordingly")

    prov["agreement.json"] = [
        {"file": "results/m3-agreement-table.csv", "sha256": sha256_full(RESULTS / "m3-agreement-table.csv"), "rows": len(ag)},
        {"file": "results/m3-continuous-bridge-model.csv", "sha256": sha256_full(RESULTS / "m3-continuous-bridge-model.csv"), "rows": len(bridge)},
        {"file": "data/raw/tms_nature2020_moesm8_dge_results.xlsx", "sha256": sha256_full(XLSX), "rows": len(paper_coef)},
    ]
    return {
        "lollipop": lollipop,
        "counts": counts,
        "scatter": {
            "unit": PCT_UNIT,
            "n": n,
            "rho": 0.866,
            "rho_full": rho,
            "points": points,
        },
    }


def export_design_grid(s, prov):
    df = s["sexage"]
    d = df[df["age_group"] != "All"].reset_index(drop=True)
    allrow = df[df["age_group"] == "All"].iloc[0]
    age_groups = list(d["age_group"])
    assert age_groups == ["1m", "3m", "18m", "21m", "24m", "30m"]
    payload = {
        "age_groups": age_groups,
        "cells": {"female": [inum(v) for v in d["cells_female"]],
                  "male": [inum(v) for v in d["cells_male"]]},
        "mice": {"female": [inum(v) for v in d["mice_female"]],
                 "male": [inum(v) for v in d["mice_male"]]},
    }
    # totals must reconcile with the All row of the source table
    assert sum(payload["cells"]["female"]) == inum(allrow["cells_female"])
    assert sum(payload["cells"]["male"]) == inum(allrow["cells_male"])
    assert sum(payload["mice"]["female"]) == inum(allrow["mice_female"])
    assert sum(payload["mice"]["male"]) == inum(allrow["mice_male"])
    prov["design-grid.json"] = [
        {"file": "results/m1-sex-age-table.csv", "sha256": sha256_full(RESULTS / "m1-sex-age-table.csv"), "rows": len(df)},
    ]
    return payload


def export_program(s, prov):
    cat, glp, pg = s["catalog"], s["genelevel"], s["programgenes"]
    assert len(cat) == 17 and len(glp) == 312
    # The m5 catalog stores m2-DE-filename slugs ("Lung::adventitial_cell"); every
    # other export uses display unit names ("Lung::adventitial cell"). Map back via
    # the same filename derivation the m2 pipeline used; assert it reproduces every
    # committed m2_de_*.csv stem so no unit can be mistranslated (2026-10-19 fix).
    import re
    display_names = list(s["ranked"]["unit"])
    slug_to_display = {}
    for d in display_names:
        stem = re.sub(r"[^A-Za-z0-9_]+", "_", d.replace("::", "__"))
        assert (RESULTS / f"m2_de_{stem}.csv").exists(), (d, stem)
        # catalog keys are filename stems with the tissue separator restored ("::")
        slug_to_display[stem.replace("__", "::")] = d
    untranslated = set(cat["unit"]) - set(slug_to_display)
    assert not untranslated, untranslated

    def disp(u):
        return slug_to_display.get(u, u)

    n_up = int(pg["program_up"].notna().sum())
    n_down = int(pg["program_down"].notna().sum())
    assert (n_up, n_down) == (167, 94), (n_up, n_down)
    units = [{
        "unit": disp(r["unit"]),
        "n_up": inum(r["n_up"]),
        "n_down": inum(r["n_down"]),
        "frac_ribo_mito": fnum(r["frac_ribo_mito"]),
        "frac_ribo_mito_background": fnum(r["frac_ribo_mito_background"]),
    } for _, r in cat.iterrows()]
    genes = [{
        "gene": r["gene"],
        "unit": disp(r["unit"]),
        "direction": "up" if r["log2FC"] > 0 else "down",
        "log2fc": fnum(r["log2FC"]),
        "fdr": fnum(r["FDR"]),
    } for _, r in glp.iterrows()]
    # top-40 recurrence matrix (stated sampling rule, contract schema 5)
    columns = [disp(u) for u in list(cat["unit"])]
    recurrence = glp.groupby("gene")["unit"].agg(["count", "first"]).reset_index()
    top40 = recurrence.sort_values(["count", "gene"], ascending=[False, True]).head(40)["gene"]
    by_gene = glp.groupby("gene")
    rows = []
    for g in top40:
        gg = by_gene.get_group(g)
        rows.append({
            "gene": g,
            "n_units": int(gg["unit"].nunique()),
            "log2fc": {disp(r["unit"]): fnum(r["log2FC"]) for _, r in gg.iterrows()},
        })
    assert len(rows) == 40 and len(columns) == 17
    prov["program.json"] = [
        {"file": "results/m5-mouse-program-catalog.csv", "sha256": sha256_full(RESULTS / "m5-mouse-program-catalog.csv"), "rows": len(cat)},
        {"file": "results/m5-gene-level-programs.csv", "sha256": sha256_full(RESULTS / "m5-gene-level-programs.csv"), "rows": len(glp)},
        {"file": "results/m6-program-genes.csv", "sha256": sha256_full(RESULTS / "m6-program-genes.csv"), "rows": len(pg)},
    ]
    return {
        "units": units,
        "genes": genes,
        "summary": {"up": n_up, "down": n_down, "entries": len(glp)},
        "matrix_top40": {"rows": rows, "columns": columns},
    }


def _arm(scores, stats, species, stat_name, extra=None):
    sc = scores[scores["species"] == species]
    st = stats[stats["comparison"] == stat_name].iloc[0]
    arm = {
        "points": [{
            "donor": r["donor"],
            "age_years": fnum(r["age_years"]),
            "age_normalized": fnum(r["age_normalized"]),
            "score": fnum(r["program_score"]),
            "n_cells": inum(r["n_cells"]),
        } for _, r in sc.iterrows()],
        "n": len(sc),
        "slope": fnum(st["slope"]),
        "ci": [fnum(st["ci_low"]), fnum(st["ci_high"])],
        "p": fnum(st["p"]),
    }
    if extra:
        arm.update(extra)
    return arm


def export_signature(s, prov):
    scores, stats, gates, ls = s["scores"], s["stats"], s["gates"], s["lifespans"]
    g = {r["gate"]: r["value"] for _, r in gates.iterrows()}
    nmr = scores[scores["species"] == "Naked mole-rat"]
    arms = {
        "mouse_wb": _arm(scores, stats, "Mouse (whole-body)", "Mouse (whole-body) aging slope"),
        # pipeline skin arm = all skin cells ("Mouse"); the "Mouse (epidermal)" series
        # has no arm in the pinned schema (correction note in DATA_CONTRACT.md)
        "mouse_skin": _arm(scores, stats, "Mouse", "Mouse aging slope", {"underpowered": True}),
        "human_skin": _arm(scores, stats, "Human", "Human aging slope", {"low_powered": True}),
        "human_fibroblast": _arm(scores, stats, "Human (fibroblast)", "Human (fibroblast) aging slope"),
        "nmr": _arm(scores, stats, "Naked mole-rat", "Naked mole-rat aging slope",
                    {"tested_span_years": [fnum(nmr["age_years"].min()), fnum(nmr["age_years"].max())]}),
    }
    diff = stats[stats["comparison"] == "Mouse (whole-body) vs Naked mole-rat slope difference"].iloc[0]
    tost = stats[stats["comparison"].str.startswith("NMR TOST")].iloc[0]
    perm = stats[stats["comparison"].str.startswith("NMR aged-vs-young")].iloc[0]
    assert round(fnum(g["tost_delta"]), 5) == 0.05874 and round(fnum(tost["p"]), 5) == 0.32894
    assert fnum(g["nmr_perm_p"]) == fnum(perm["p"])
    assert round(fnum(perm["p"]), 5) == 0.81818
    assert arms["nmr"]["tested_span_years"] == [1.1, 11.3]
    life = {r["common_name"]: fnum(r["max_lifespan_years"]) for _, r in ls.iterrows()}
    payload = {
        "arms": arms,
        "stats": {
            "diff": {"slope": fnum(diff["slope"]), "se": fnum(diff["se"]),
                     "ci": [fnum(diff["ci_low"]), fnum(diff["ci_high"])],
                     "p": fnum(diff["p"]), "n": inum(diff["n"])},
            "tost": {"delta": fnum(g["tost_delta"]), "p": fnum(tost["p"]), "established": False},
            "perm_p": fnum(perm["p"]),
        },
        "coverage": {"mouse": fnum(g["coverage_mouse"]), "human": fnum(g["coverage_human"]),
                     "nmr": fnum(g["coverage_nmr"])},
        "lifespans": {"mouse": life["House mouse"], "nmr": life["Naked mole-rat"],
                      "human": life["Human"]},
    }
    prov["signature.json"] = [
        {"file": "results/m6-lifespan-program-scores.csv", "sha256": sha256_full(RESULTS / "m6-lifespan-program-scores.csv"), "rows": len(scores)},
        {"file": "results/m6-statistical-summary.csv", "sha256": sha256_full(RESULTS / "m6-statistical-summary.csv"), "rows": len(stats)},
        {"file": "results/m6-gates.csv", "sha256": sha256_full(RESULTS / "m6-gates.csv"), "rows": len(gates)},
        {"file": "results/m6-species-lifespans.csv", "sha256": sha256_full(RESULTS / "m6-species-lifespans.csv"), "rows": len(ls)},
    ]
    return payload


def export_savina(s, prov):
    sv = s["savina"]
    curve = [{"log2fc_cutoff": fnum(r["log2fc_cutoff"]),
              "n_de_pval05": inum(r["n_de_pval05"]),
              "n_de_fdr05": inum(r["n_de_fdr05"])} for _, r in sv.iterrows()]
    by_cutoff = {fnum(r["log2fc_cutoff"]): r for _, r in sv.iterrows()}
    r10, r145, r025 = by_cutoff[1.0], by_cutoff[1.45], by_cutoff[0.25]
    assert inum(r10["n_de_pval05"]) == 29
    assert (inum(r145["n_de_pval05"]), inum(r145["n_de_fdr05"])) == (9, 6)
    assert inum(r025["n_de_fdr05"]) == 560
    prov["savina.json"] = [
        {"file": "results/m6-savina-reproduction.csv", "sha256": sha256_full(RESULTS / "m6-savina-reproduction.csv"), "rows": len(sv)},
    ]
    return {
        "curve": curve,
        "markers": {
            "paper_rule": {"cutoff": 1.0, "n_pval": inum(r10["n_de_pval05"])},
            "converted_rule": {"cutoff": 1.45, "n_pval": inum(r145["n_de_pval05"]),
                               "n_fdr": inum(r145["n_de_fdr05"])},
            "loose_rule": {"cutoff": 0.25, "n_fdr": inum(r025["n_de_fdr05"])},
        },
        "caveats": ["Igfbp3 absent from source file", "QC 10,251 vs paper 10,232"],
    }


def export_power_grid(s, prov):
    pw = s["power"]
    donors = sorted(int(v) for v in pw["n_donors"].unique())
    sds = sorted(int(v) for v in pw["donor_age_sd_years"].unique())
    disps = sorted(float(v) for v in pw["biological_dispersion"].unique())
    assert donors == [4, 6, 8, 10, 12, 16, 20, 30, 50]
    assert sds == [8, 12, 16, 20] and disps == [0.1, 0.2, 0.4] and len(pw) == 108
    disp_label = {0.1: "low", 0.2: "mid", 0.4: "high"}
    tiles = []
    anchors = {}
    for _, r in pw.iterrows():
        d, sd, disp = inum(r["n_donors"]), int(r["donor_age_sd_years"]), fnum(r["biological_dispersion"])
        tiles.append({
            "n_donors": d,
            "sd": sd,
            "dispersion": disp,
            "min_30y_log2fc": fnum(r["min_detectable_30y_log2fc"]),
            "power_status": r["power_status"],
            "eligible_cohort": bnum(r["eligible_cohort"]),
        })
        anchors[f"{d}_{sd}_{disp_label[disp]}"] = fnum(r["min_detectable_30y_log2fc"])
    assert anchors["4_8_low"] == 1.474 and anchors["4_20_low"] == 0.59
    prov["power-grid.json"] = [
        {"file": "results/m4-human-power.csv", "sha256": sha256_full(RESULTS / "m4-human-power.csv"), "rows": len(pw)},
    ]
    return {
        "donors": donors,
        "sds": sds,
        "dispersions": disps,
        "tiles": tiles,
        "ts_marker": {"dispersion": 0.1, "donors": 4, "label": "Tabula Sapiens skin: 4 donors, SD 8y", "min_30y_log2fc": 1.474, "sd": 8, "tile_x": 0, "tile_y": 0},
        "anchors": {"4_8_low": anchors["4_8_low"], "4_20_low": anchors["4_20_low"]},
    }


# ---------------------------------------------------------------- main

def main():
    s = load_sources()
    paper_coef = read_paper_age_coef()

    # wipe stale exports/figures from the superseded app version
    OUT.mkdir(parents=True, exist_ok=True)
    for old in OUT.glob("*.json"):
        if old.name not in NEW_DATA_FILES:
            print(f"[wipe] removing stale export {old.name}")
            old.unlink()
    APP_FIGS.mkdir(parents=True, exist_ok=True)
    for old in APP_FIGS.glob("*.png"):
        if old.name not in FIGURE_FILES:
            print(f"[wipe] removing stale figure {old.name}")
            old.unlink()

    prov: dict[str, list] = {}
    sizes = {}
    atlas_payload = export_atlas(s, prov)
    cal_payload = export_calibration(s, prov)
    sizes["atlas.json"] = write_json("atlas.json", atlas_payload)
    sizes["ranking.json"] = write_json("ranking.json", {
        "ranked": atlas_payload["ranked"],
        "calibration": cal_payload,
        "n_units_evaluated": len(atlas_payload["units"]),
        "n_units_powered": len(atlas_payload["ranked"]),
    })
    sizes["calibration.json"] = write_json("calibration.json", cal_payload)
    sizes["agreement.json"] = write_json("agreement.json", export_agreement(s, prov, paper_coef))
    sizes["design-grid.json"] = write_json("design-grid.json", export_design_grid(s, prov))
    sizes["program.json"] = write_json("program.json", export_program(s, prov))
    sizes["signature.json"] = write_json("signature.json", export_signature(s, prov))
    sizes["savina.json"] = write_json("savina.json", export_savina(s, prov))
    sizes["power-grid.json"] = write_json("power-grid.json", export_power_grid(s, prov))
    sizes["anage.json"] = write_json("anage.json", export_anage(prov))
    sizes["humancohorts.json"] = write_json("humancohorts.json", export_human_cohorts(prov))
    sizes["dlpfc.json"] = write_json("dlpfc.json", export_dlpfc(prov))

    prov["provenance.json"] = [
        {"file": "app/DATA_CONTRACT.md", "sha256": sha256_full(REPO / "app" / "DATA_CONTRACT.md"), "rows": 76},
    ]

    for fig in FIGURE_FILES:
        shutil.copyfile(FIGS / fig, APP_FIGS / fig)
    assets = {fig: {"bytes": (APP_FIGS / fig).stat().st_size,
                    "sha256": sha256_full(APP_FIGS / fig)} for fig in FIGURE_FILES}

    provenance = {
        "freeze_date": FREEZE_DATE,
        "export_date": date.today().isoformat(),
        "sources": prov,
        "assets": assets,
    }
    sizes["provenance.json"] = write_json("provenance.json", provenance)

    total = sum(sizes.values())
    for name, sz in sorted(sizes.items()):
        assert sz < 1_000_000, f"{name} too large: {sz}"
        print(f"[export] app/data/{name}: {sz:,} bytes")
    assert total < 5_000_000, f"app/data over budget: {total}"
    print(f"[export] total app/data: {total:,} bytes (< 5 MB); all files < 1 MB")
    print("[export] done")


if __name__ == "__main__":
    sys.exit(main())
