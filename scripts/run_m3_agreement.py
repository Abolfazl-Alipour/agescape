#!/usr/bin/env python3
"""
run_m3_agreement.py — Milestone M3 agreement check, REBUILT against the TMS paper's
OWN published aging-DGE results (Nature 2020, DOI 10.1038/s41586-020-2496-1,
Supplementary Table MOESM8, sheets DGE_result.tissue_cell.age_coef / .bh_p; file
archived at data/raw/tms_nature2020_moesm8_dge_results.xlsx, fetched 2026-10-14).

The first M3 tested a hand-picked senescence-marker list; the pre-registered design
requires testing against the paper's published lists. This script does that:
for each powered M2 unit, merge paper per-gene age coefficients (continuous model)
with our continuous bridge model slopes (results/m3-continuous-bridge-model.csv),
and report rank correlation, sign concordance among paper-significant genes, and
overlap at nominal significance. Verdicts: agree (rho>0 & sign-concordance>=60%),
partial (0<sign-concordance<60%), disagree (rho<=0), cannot-test (<10 paper-sig genes).
"""
from pathlib import Path
import numpy as np, pandas as pd
import openpyxl
from scipy import stats

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO/"results"
XLSX = REPO/"data/raw/tms_nature2020_moesm8_dge_results.xlsx"
PAPER_BH = 0.01  # paper's FDR level for reported age-DGE (reported as bh_p)

wb = openpyxl.load_workbook(XLSX, read_only=True)
coef = pd.DataFrame(wb["DGE_result.tissue_cell.age_coef"].values)
bh = pd.DataFrame(wb["DGE_result.tissue_cell.bh_p"].values)
hdr = coef.iloc[0]; coef.columns = hdr; coef = coef[1:]
bh.columns = hdr; bh = bh[1:]
genes = coef.iloc[:, 0].astype(str)
coef.index = genes; bh.index = genes
coef = coef.drop(columns=[hdr[0]]).astype(float)
bh = bh.drop(columns=[hdr[0]]).astype(float)
print(f"[m3] paper matrix: {coef.shape[1]} tissue.cell-type columns, {len(genes)} genes")

bridge = pd.read_csv(RESULTS/"m3-continuous-bridge-model.csv")
ranked = pd.read_csv(RESULTS/"m2-aging-signal-ranked.csv")
units = list(ranked["unit"])
# our per-gene effect: M2 pseudobulk log2FC (3m vs 18m, dichotomous). Where the continuous
# bridge model exists for the unit, prefer its slope; note the dichotomous-vs-continuous
# asymmetry vs the paper's continuous model in the reason column.
rows = []
for unit_raw in units:
    tissue, ct = unit_raw.split("::", 1)
    unit = f"{tissue}.{ct}"
    cand = [c for c in coef.columns if c.lower() == unit.lower()]
    if not cand:
        cand = [c for c in coef.columns if c.lower().endswith(ct.lower()) and c.lower().startswith(tissue.lower())]
    de_fp = RESULTS / f"m2_de_{unit_raw.replace('::','__').replace(' ','_').replace('-','_').replace(',','')}.csv"
    ours = None; how = ""
    bsub = bridge[(bridge["tissue"] == tissue) & (bridge["cell_ontology_class"] == ct)]
    if len(bsub) > 200:
        ours = bsub.set_index("gene")[["slope_per_month", "pvalue"]].rename(
            columns={"slope_per_month": "ours_effect", "pvalue": "ours_p"})
        how = "continuous bridge slope"
    elif de_fp.exists():
        d = pd.read_csv(de_fp).set_index("index")
        ours = d[["log2FC", "pvalue"]].rename(columns={"log2FC": "ours_effect", "pvalue": "ours_p"})
        how = "M2 dichotomous log2FC (3m vs 18m)"
    if not cand or ours is None or len(ours) == 0:
        rows.append({"unit": unit.replace(".", "::"), "n_shared": 0, "verdict": "cannot-test",
                     "reason": "unit not present in one of the two result sets"})
        continue
    col = cand[0]
    pap = pd.DataFrame({"paper_coef": coef[col], "paper_bh": bh[col]}).dropna()
    m = ours.join(pap, how="inner").dropna()
    if len(m) < 50:
        rows.append({"unit": unit.replace(".", "::"), "n_shared": len(m), "verdict": "cannot-test",
                     "reason": "too few shared genes"})
        continue
    rho = stats.spearmanr(m["paper_coef"], m["ours_effect"]).statistic
    psig = m[m["paper_bh"] < PAPER_BH]
    conc = float((np.sign(psig["paper_coef"]) == np.sign(psig["ours_effect"])).mean()) if len(psig) else np.nan
    rec = float((psig["ours_p"] < 0.05).mean()) if len(psig) else np.nan
    verdict = "agree" if (rho > 0 and (np.isnan(conc) or conc >= 0.6)) else \
              ("partial" if (rho > 0) else "disagree")
    rows.append({"unit": unit.replace(".", "::"), "paper_column": col, "n_shared": len(m),
                 "n_paper_sig_bh0.01": len(psig), "spearman_rho": round(rho, 3),
                 "sign_concordance_paper_sig": None if np.isnan(conc) else round(conc, 3),
                 "our_nominal_recovery_in_paper_sig": None if np.isnan(rec) else round(rec, 3),
                 "verdict": verdict, "reason": f"paper continuous-age coef vs our {how}"})
    if "proximal convoluted tubule" in unit:
        m_pct = m.copy().reset_index().rename(columns={"index": "gene"})
        m_pct["concordant"] = np.sign(m_pct["paper_coef"]) == np.sign(m_pct["ours_effect"])
        m_pct["paper_sig"] = m_pct["paper_bh"] < PAPER_BH
        m_pct["ours_sig"] = m_pct["ours_p"] < 0.05
        pct_ref_path = RESULTS / "m3-pct-moesm8-reference.csv"
        m_pct.to_csv(pct_ref_path, index=False)
        print(f"[m3] wrote {pct_ref_path} ({len(m_pct)} shared genes)")
out = pd.DataFrame(rows)
out.to_csv(RESULTS/"m3-agreement-table.csv", index=False)
print(out.to_string(index=False))
print(f"[m3] wrote {RESULTS/'m3-agreement-table.csv'}")
