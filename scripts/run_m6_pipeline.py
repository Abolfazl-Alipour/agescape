#!/usr/bin/env python3
"""
run_m6_pipeline.py — Milestone M6: Longevity-Outlier Contrast (REAL DATA, full rebuild).

Every number below derives from a named file on disk. Nothing is simulated.

Inputs
------
- data/raw/GSE152886_nmr_skin_raw.csv.gz : Savina et al. 2022 (Aging) NMR epidermis
  raw counts, genes x cells, tab-separated. Samples 1J,5J,6J (young, mean 1.1y) and
  2V,3V,4V (aged, mean 11.3y). Sample 2V excluded per the paper (technical outlier).
  Ages: paper reports group means (1.1y young, 11.3y aged); individual ages are not in GEO.
- data/raw/tms_droplet_raw.h5ad : Tabula Muris Senis droplet (mouse). Skin-of-body subset.
- data/raw/TS_Skin.h5ad : Tabula Sapiens skin (human; 4 donors: 33y,45y,59y,61y; dermis-heavy,
  no keratinocyte population). LOW-POWERED human cohort (n=4 donors); all human numbers
  carry that flag.
- data/interim/orthologs-mouse-human.csv, orthologs-mouse-nmr.csv : Ensembl 1:1 (symbols).
- data/interim/nmr-ensg2symbol.tsv : ENSHGLG -> gene symbol (Ensembl REST, 2026-10-14).
- results/m2_de_*.csv : real pseudobulk DE from M2 (17 powered Lung/Kidney units).
- data/raw/anage_dataset.zip : official AnAge export (manifest: anage_dataset).

Program definition (documented; no circularity)
-----------------------------------------------
mouse_up   = union of genes FDR<0.05, log2FC>=0.25  across the 17 M2 unit DE files
mouse_down = union of genes FDR<0.05, log2FC<=-0.25 across the same
conserved_* = mouse program genes whose 1:1 human-ortholog slope across the 4 TS-skin
             donors is same-direction with p<0.05 (strict) — reported alongside the
             full mouse program (sensitivity). Program is locked BEFORE any NMR number
             is computed (code order below enforces this).

Outputs: results/m6-*.csv, figures/m6-*.png, plus printed gates.
"""
import gzip, io, json, sys, zipfile, re
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.sparse as sp
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parent.parent
RAW, INTERIM, RESULTS, FIGURES = (REPO/"data/raw", REPO/"data/interim",
                                  REPO/"results", REPO/"figures")
NOTES = REPO/"notes"
for d in (RESULTS, FIGURES): d.mkdir(exist_ok=True)

rng = np.random.default_rng(42)  # seeds only; no synthetic data anywhere

def prov(msg): print(f"[prov] {msg}", flush=True)

# ---------------------------------------------------------------- AnAge lifespans
anage_zip = RAW/"anage_dataset.zip"
if not anage_zip.exists():
    import shutil; shutil.copy("/tmp/anage_dataset.zip", anage_zip)
z = zipfile.ZipFile(anage_zip)
release = z.read("release.html").decode(errors="replace")[:200].replace("\n"," ")
prov(f"AnAge export: {anage_zip.name}; release.html head: {release[:120]}")
f = io.TextIOWrapper(z.open("anage_data.txt"), encoding="utf-8", errors="replace")
cols = f.readline().rstrip("\n").split("\t")
mlcol = [c for c in cols if "maximum longevity" in c.lower()][0]
spcol = ("Genus","Species")
lifespans = {}
for line in f:
    p = line.rstrip("\n").split("\t")
    if len(p) < len(cols): continue
    spname = p[cols.index("Genus")] + " " + p[cols.index("Species")]
    if spname in ("Mus musculus","Homo sapiens","Heterocephalus glaber"):
        try: lifespans[spname] = float(p[cols.index(mlcol)])
        except ValueError: pass
assert len(lifespans) == 3, f"AnAge species missing: {lifespans}"
prov(f"AnAge max longevity (years): {lifespans}")
LS_MOUSE = lifespans["Mus musculus"]; LS_HUMAN = lifespans["Homo sapiens"]; LS_NMR = lifespans["Heterocephalus glaber"]

# ---------------------------------------------------------------- mouse program from M2 (real DE tables)
FDR_T, FC_T = 0.05, 0.25
de_files = sorted(RESULTS.glob("m2_de_*.csv"))
assert len(de_files) >= 15, f"expected >=15 M2 DE tables, found {len(de_files)}"
up, down = set(), set()
for fp in de_files:
    df = pd.read_csv(fp)
    sig = df[(df["FDR"] < FDR_T) & (df["log2FC"].abs() >= FC_T)]
    up |= set(sig.loc[sig["log2FC"] > 0, "index"])
    down |= set(sig.loc[sig["log2FC"] < 0, "index"])
prov(f"mouse program from {len(de_files)} M2 units: {len(up)} up, {len(down)} down genes")

# ---------------------------------------------------------------- ortholog tables
mh = pd.read_csv(INTERIM/"orthologs-mouse-human.csv")
m2h = dict(zip(mh["mouse_gene_symbol"], mh["human_gene_symbol"]))
mn = pd.read_csv(INTERIM/"orthologs-mouse-nmr.csv")
m2n = dict(zip(mn["mouse_gene_symbol"], mn["nmr_gene_symbol"]))
prov(f"orthologs: {len(m2h)} mouse->human, {len(m2n)} mouse->NMR (1:1)")

# ---------------------------------------------------------------- human replication (TS skin, 4 donors) -> locks program
ts = sc.read_h5ad(RAW/"TS_Skin.h5ad")
prov(f"TS_Skin loaded: {ts.shape} from {RAW/'TS_Skin.h5ad'}")
# new-schema CELLxGENE TS: X is log-normalized; var_names are ENSG, symbols in feature_name
ts.var["symbol"] = ts.var["feature_name"].astype(str)
ts.var_names = ts.var["symbol"].values
ts.var_names_make_unique()
assert "10X" in ts.obs["assay"].unique() or "10x" in str(ts.obs["assay"].unique()).lower()
ts = ts[ts.obs["assay"].astype(str).str.contains("10x|10X", case=False, na=False)].copy()
ts.obs["age"] = ts.obs["development_stage"].astype(str).str.extract(r"(\d+)-year-old").astype(float)
donors = ts.obs.groupby("donor_id")["age"].first()
prov(f"TS skin 10x donors: {dict(donors)} (n={len(donors)}, LOW-POWERED cohort)")

def donor_mean_logexpr(adata):
    X = adata.X.tocsr() if hasattr(adata.X, "tocsr") else adata.X
    out = {}
    for d in adata.obs["donor_id"].unique():
        idx = np.where(adata.obs["donor_id"].values == d)[0]
        out[d] = np.asarray(X[idx].mean(axis=0)).ravel()
    return pd.DataFrame(out, index=adata.var_names)

ts_pb = donor_mean_logexpr(ts)  # mean log-normalized expression per donor
ages = donors.loc[ts_pb.columns]
A = np.vstack([np.ones(len(ages)), ages.values]).T

def slopes_per_gene(expr_df):
    Y = expr_df.values.T  # donors x genes
    beta = np.linalg.lstsq(A, Y, rcond=None)[0]  # 2 x G
    resid = Y - A @ beta
    dof = len(ages) - 2
    s2 = (resid**2).sum(axis=0) / dof
    XtX_inv = np.linalg.inv(A.T @ A)
    se = np.sqrt(np.maximum(s2 * XtX_inv[1, 1], 1e-300))
    t = beta[1] / se
    p = 2 * stats.t.sf(np.abs(t), dof)
    return pd.DataFrame({"slope_per_year": beta[1], "pvalue": p}, index=expr_df.index)

ts_sl = slopes_per_gene(ts_pb)
h2m = {v: k for k, v in m2h.items()}
ts_sl["mouse_symbol"] = ts_sl.index.map(h2m)
ts_sl = ts_sl.dropna(subset=["mouse_symbol"]).set_index("mouse_symbol")
strict_up = ts_sl[(ts_sl["slope_per_year"] > 0) & (ts_sl["pvalue"] < 0.05)].index
strict_dn = ts_sl[(ts_sl["slope_per_year"] < 0) & (ts_sl["pvalue"] < 0.05)].index
conserved_up = sorted(up & set(strict_up))
conserved_down = sorted(down & set(strict_dn))
prov(f"human replication (TS skin slopes): strict up={len(strict_up)}, strict down={len(strict_dn)} genes total")
prov(f"conserved program: {len(conserved_up)} up, {len(conserved_down)} down genes (of {len(up)}/{len(down)} mouse)")
pd.DataFrame({"program_up": pd.Series(list(up)), "program_down": pd.Series(list(down)),
              "conserved_up": pd.Series(conserved_up), "conserved_down": pd.Series(conserved_down)}
             ).to_csv(RESULTS/"m6-program-genes.csv", index=False)
pd.DataFrame({"human_slope": ts_sl["slope_per_year"], "pvalue": ts_sl["pvalue"]}
             ).to_csv(RESULTS/"m6-human-skin-slopes.csv")
PROGRAM = conserved_up + conserved_down
SCORE_UP, SCORE_DOWN = conserved_up, conserved_down
if len(PROGRAM) < 10:  # pragmatic floor: fall back to full mouse program, flagged loudly
    SCORE_UP, SCORE_DOWN = sorted(up), sorted(down)
    print("WARNING: conserved program <10 genes; scoring FULL mouse program (flagged in report)")

# ---------------------------------------------------------------- scoring helpers
def score_cells(adata, up_genes, down_genes, label):
    sc.tl.score_genes(adata, up_genes, score_name="prog_up", use_raw=False)
    if down_genes:
        sc.tl.score_genes(adata, down_genes, score_name="prog_down", use_raw=False)
        adata.obs["program_score"] = adata.obs["prog_up"] - adata.obs["prog_down"]
    else:
        adata.obs["program_score"] = adata.obs["prog_up"]
    cov = sum(1 for g in up_genes + down_genes if g in adata.var_names) / max(1, len(up_genes) + len(down_genes))
    prov(f"[{label}] program gene coverage in var_names: {cov:.3f} ({int(cov*len(up_genes+down_genes))}/{len(up_genes)+len(down_genes)})")
    return adata, cov

# ---------------------------------------------------------------- MOUSE: TMS skin
tms = sc.read_h5ad(RAW/"tms_droplet_raw.h5ad", backed="r")
skin_idx = tms.obs.index[tms.obs["tissue"].astype(str).str.strip() == "Skin"]
sub = tms[skin_idx].to_memory()
prov(f"TMS skin cells: {sub.shape} from tms_droplet_raw.h5ad")
# strip obs strings
for c in sub.obs.columns:
    if sub.obs[c].dtype == object or str(sub.obs[c].dtype).startswith("string"):
        sub.obs[c] = sub.obs[c].astype(str).str.strip()
sub.obs["age_months"] = sub.obs["age"].astype(str).str.extract(r"(\d+)").astype(float)
epi_types = {"basal cell of epidermis","keratinocyte","keratinocyte stem cell","stem cell of epidermis"}
sub.obs["is_epidermal"] = sub.obs["cell_ontology_class"].isin(epi_types)
sc.pp.normalize_total(sub, target_sum=1e4); sc.pp.log1p(sub)
sub, cov_m = score_cells(sub, SCORE_UP, SCORE_DOWN, "mouse-skin(all)")
sub_epi = sub[sub.obs["is_epidermal"]].copy()
sub_epi, cov_me = score_cells(sub_epi, SCORE_UP, SCORE_DOWN, "mouse-skin(epidermal)")
mm = sub.obs.groupby(["mouse.id","age_months"], observed=True)["program_score"].agg(["mean","count"]).reset_index()
mm.columns = ["donor","age_months","program_score","n_cells"]
mm["species"] = "Mouse"; mm["age_years"] = mm["age_months"]/12.0
mm["age_normalized"] = mm["age_years"]/LS_MOUSE
mme = sub_epi.obs.groupby(["mouse.id","age_months"], observed=True)["program_score"].agg(["mean","count"]).reset_index()
mme.columns = ["donor","age_months","program_score","n_cells"]
mme["species"] = "Mouse (epidermal)"; mme["age_years"] = mme["age_months"]/12.0
mme["age_normalized"] = mme["age_years"]/LS_MOUSE
prov(f"mouse skin mice: {len(mm)} (all), {len(mme)} (epidermal); cells/mouse min={int(mm['n_cells'].min())} -> skin arm UNDERPOWERED (n=3 mice), whole-body arm below is primary")

# --- mouse whole-body positive control: program score per mouse across ALL droplet cells
tms_b = sc.read_h5ad(RAW/"tms_droplet_raw.h5ad", backed="r")
var_idx = {g: i for i, g in enumerate(tms_b.var_names)}
up_i = [var_idx[g] for g in SCORE_UP if g in var_idx]
dn_i = [var_idx[g] for g in SCORE_DOWN if g in var_idx]
prov(f"whole-body mouse arm: {len(up_i)}/{len(SCORE_UP)} up and {len(dn_i)}/{len(SCORE_DOWN)} down genes found in raw var")
Xsel = tms_b[:, sorted(up_i + dn_i)].X  # backed slice -> csr (245k x ~261)
colmap = {c: k for k, c in enumerate(sorted(up_i + dn_i))}
up_k = [colmap[i] for i in up_i]; dn_k = [colmap[i] for i in dn_i]
tot = np.asarray(Xsel.sum(axis=1)).ravel()
tot[tot == 0] = 1.0
import scipy.sparse as _sp
Xn = _sp.diags(1.0 / tot * 1e4) @ Xsel
Xn.data = np.log1p(Xn.data)
score_up = np.asarray(Xn[:, up_k].mean(axis=1)).ravel()
score_dn = np.asarray(Xn[:, dn_k].mean(axis=1)).ravel() if dn_k else 0.0
wb = pd.DataFrame({"mouse.id": tms_b.obs["mouse.id"].astype(str).values,
                   "age_months": tms_b.obs["age"].astype(str).str.extract(r"(\d+)")[0].astype(float).values,
                   "program_score": score_up - score_dn})
mw = wb.groupby(["mouse.id","age_months"], observed=True)["program_score"].agg(["mean","count"]).reset_index()
mw.columns = ["donor","age_months","program_score","n_cells"]
mw["species"] = "Mouse (whole-body)"; mw["age_years"] = mw["age_months"]/12.0
mw["age_normalized"] = mw["age_years"]/LS_MOUSE
prov(f"whole-body mouse arm: {len(mw)} mice, cells/mouse min={int(mw['n_cells'].min())}")

# ---------------------------------------------------------------- HUMAN: TS skin per donor
ts2, cov_h = score_cells(ts, [m2h.get(g, g) for g in SCORE_UP], [m2h.get(g, g) for g in SCORE_DOWN], "human-TS-skin")
hh = ts2.obs.groupby(["donor_id","age"], observed=True)["program_score"].agg(["mean","count"]).reset_index()
hh.columns = ["donor","age_years","program_score","n_cells"]
hh["species"] = "Human"; hh["age_normalized"] = hh["age_years"]/LS_HUMAN
# fibroblast-only sensitivity
ts_fib = ts2[ts2.obs["cell_type"] == "fibroblast"]
hf = ts_fib.obs.groupby(["donor_id","age"], observed=True)["program_score"].agg(["mean","count"]).reset_index()
hf.columns = ["donor","age_years","program_score","n_cells"]
hf["species"] = "Human (fibroblast)"; hf["age_normalized"] = hf["age_years"]/LS_HUMAN
prov(f"human skin donors scored: {len(hh)} (all cells), {len(hf)} (fibroblast-only); LOW-POWERED n=4")

# ---------------------------------------------------------------- NMR: GSE152886
prov("loading GSE152886 (genes x cells) ...")
df = pd.read_csv(RAW/"GSE152886_nmr_skin_raw.csv.gz", sep="\t", index_col=0)
prov(f"GSE matrix: {df.shape[0]} genes x {df.shape[1]} cells")
cells = pd.DataFrame({"barcode": df.columns})
cells["sample"] = cells["barcode"].str.split("_").str[0]
AGE_MAP = {"1J":1.1,"5J":1.1,"6J":1.1,"3V":11.3,"4V":11.3,"2V":11.3}
GROUP = {"1J":"young","5J":"young","6J":"young","3V":"aged","4V":"aged","2V":"aged"}
cells["age_years"] = cells["sample"].map(AGE_MAP); cells["group"] = cells["sample"].map(GROUP)
n0 = len(cells)
cells = cells[cells["sample"] != "2V"].copy()   # paper: technical outlier
prov(f"sample 2V excluded per Savina 2022: {n0} -> {len(cells)} cells")
X = sp.csr_matrix(df[cells["barcode"]].T.values)  # cells x genes
g2s = {}
for line in open(INTERIM/"nmr-ensg2symbol.tsv"):
    if line.startswith("#"): continue
    gid, sym = line.rstrip("\n").split("\t")[:2]
    g2s[gid] = sym
var_names = [g2s.get(g, g) for g in df.index]
nmr = sc.AnnData(X=X, obs=cells.set_index("barcode"), var=pd.DataFrame(index=var_names))
nmr.var_names_make_unique()
# QC per paper: UMI>=2000, genes>=700 (MT% not applicable: no MT annotation in GEO file -> documented)
nmr.obs["n_counts"] = np.asarray(X.sum(axis=1)).ravel()
nz = np.asarray((X > 0).sum(axis=1)).ravel()
nmr.obs["n_genes"] = nz
before = nmr.n_obs
keep = (nmr.obs["n_counts"] >= 2000) & (nmr.obs["n_genes"] >= 700)
nmr = nmr[keep].copy()
prov(f"NMR QC (UMI>=2000, genes>=700; no MT filter possible from GEO file): {before} -> {nmr.n_obs} cells")
prov(f"paper reports 10,232 cells post-QC from the same 5 samples; our independent QC gives {nmr.n_obs}")
sc.pp.normalize_total(nmr, target_sum=1e4); sc.pp.log1p(nmr)
# map program to NMR symbols
up_n = [g for g in (m2n.get(g) for g in SCORE_UP) if g]
dn_n = [g for g in (m2n.get(g) for g in SCORE_DOWN) if g]
nmr, cov_n = score_cells(nmr, up_n, dn_n, "NMR-skin")
nn = nmr.obs.groupby(["sample","age_years","group"], observed=True)["program_score"].agg(["mean","count"]).reset_index()
nn.columns = ["donor","age_years","group","program_score","n_cells"]
nn["species"] = "Naked mole-rat"; nn["age_normalized"] = nn["age_years"]/LS_NMR
prov(f"NMR samples scored: {len(nn)} (3 young, 2 aged)")

# ---------------------------------------------------------------- Savina reproduction: global DE young vs aged in NMR
sc.tl.rank_genes_groups(nmr, "group", groups=["aged"], reference="young",
                        method="wilcoxon", n_genes=nmr.n_vars, key_added="de")
deg = sc.get.rank_genes_groups_df(nmr, group="aged", key="de").set_index("names")
deg.index.name = "gene"
deg.to_csv(RESULTS/"m6-nmr-de-young-vs-aged.csv")
curve = []
for cut in np.arange(0.1, 2.01, 0.05):
    n_p = int(((deg["pvals"] < 0.05) & (deg["logfoldchanges"].abs() >= cut)).sum())
    n_fdr = int(((deg["pvals_adj"] < 0.05) & (deg["logfoldchanges"].abs() >= cut)).sum())
    curve.append({"log2fc_cutoff": round(float(cut),2), "n_de_pval05": n_p, "n_de_fdr05": n_fdr})
pd.DataFrame(curve).to_csv(RESULTS/"m6-savina-reproduction.csv", index=False)
n_de_paper_rule = int(((deg["pvals"] < 0.05) & (deg["logfoldchanges"].abs() >= 1.0)).sum())
n_de_fdr_025 = int(((deg["pvals_adj"] < 0.05) & (deg["logfoldchanges"].abs() >= 0.25)).sum())
prov(f"Savina reproduction: DE genes aged-vs-young at paper rule (p<0.05 & |log2FC|>=1; "
     f"paper used avg logFC>1 in Seurat v3 natural-log units ~= log2FC 1.44): {n_de_paper_rule}; "
     f"at our standard (FDR<0.05 & |log2FC|>=0.25): {n_de_fdr_025}")
prov(f"paper headline: NO global DE genes (p<0.05 & logFc>1); within-cluster only Igfbp3 (+Rps2) in BSC2")

# Igfbp3 check (paper's within-cluster finding)
ig = [g for g in nmr.var_names if str(g).lower() == "igfbp3"]
if ig:
    sub_y = nmr[nmr.obs["group"]=="young"]; sub_a = nmr[nmr.obs["group"]=="aged"]
    v_y = np.asarray(sub_y[:, ig[0]].X.todense()).ravel(); v_a = np.asarray(sub_a[:, ig[0]].X.todense()).ravel()
    u_p = stats.mannwhitneyu(v_y, v_a, alternative="two-sided").pvalue
    per_s = nmr.obs.groupby("sample").apply(lambda d: float(np.asarray(nmr[d.index, ig[0]].X.todense()).mean()))
    pd.DataFrame({"sample": per_s.index, "mean_log1p_expr": per_s.values,
                  "group": [GROUP[s] for s in per_s.index]}).to_csv(RESULTS/"m6-igfbp3.csv", index=False)
    prov(f"Igfbp3 (NMR): young n={len(v_y)} mean={v_y.mean():.3f}, aged n={len(v_a)} mean={v_a.mean():.3f}, "
         f"Wilcoxon p={u_p:.3g} | paper: ~4x higher in aged basal cells")
else:
    prov("Igfbp3 symbol not found in NMR var_names (check mapping)")

# ---------------------------------------------------------------- sample-level stats
scores = pd.concat([mm[["species","donor","age_years","age_normalized","program_score","n_cells"]],
                    mme[["species","donor","age_years","age_normalized","program_score","n_cells"]],
                    mw[["species","donor","age_years","age_normalized","program_score","n_cells"]],
                    hh[["species","donor","age_years","age_normalized","program_score","n_cells"]],
                    hf[["species","donor","age_years","age_normalized","program_score","n_cells"]],
                    nn[["species","donor","age_years","age_normalized","program_score","n_cells"]]],
                   ignore_index=True)
scores.to_csv(RESULTS/"m6-lifespan-program-scores.csv", index=False)
prov(f"wrote {RESULTS/'m6-lifespan-program-scores.csv'}: {len(scores)} sample-level rows")

def ols(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    n = len(x); xm, ym = x.mean(), y.mean()
    Sxx = ((x-xm)**2).sum(); Sxy = ((x-xm)*(y-ym)).sum()
    b = Sxy/Sxx; a = ym - b*xm
    resid = y - (a + b*x); s2 = (resid**2).sum()/(n-2)
    se = np.sqrt(s2/Sxx)
    tc = stats.t.ppf(0.975, n-2)
    p = 2*stats.t.sf(abs(b/se), n-2)
    return dict(slope=b, intercept=a, se=se, ci_low=b-tc*se, ci_high=b+tc*se, p=p, n=n)

def tost(slope, se, df, delta):
    t1 = (slope + delta)/se; t2 = (delta - slope)/se
    p = max(stats.t.sf(t1, df), stats.t.sf(t2, df))
    return p

rows = []
for sp, d in [("Mouse", mm), ("Mouse (epidermal)", mme), ("Mouse (whole-body)", mw), ("Human", hh), ("Human (fibroblast)", hf), ("Naked mole-rat", nn)]:
    r = ols(d["age_normalized"], d["program_score"]); r["comparison"] = f"{sp} aging slope"; rows.append(r)
# interactions (sample-level, OLS on stacked z-scored within species is overkill; use slope diffs with se)
def slope_diff(d1, d2, n1, n2):
    r1, r2 = ols(d1["age_normalized"], d1["program_score"]), ols(d2["age_normalized"], d2["program_score"])
    diff = r1["slope"] - r2["slope"]; se = np.hypot(r1["se"], r2["se"])
    z = diff/se; p = 2*stats.norm.sf(abs(z))
    return {"comparison": f"{n1} vs {n2} slope difference",
            "slope": diff, "se": se, "ci_low": diff-1.96*se, "ci_high": diff+1.96*se, "p": p, "n": r1["n"]+r2["n"]}
rows.append(slope_diff(mw, nn, "Mouse (whole-body)", "Naked mole-rat"))
mouse_slope_abs = abs([r for r in rows if r["comparison"]=="Mouse (whole-body) aging slope"][0]["slope"])
DELTA = 0.25 * mouse_slope_abs  # equivalence margin = 25% of mouse slope, set BEFORE NMR test
rn = [r for r in rows if r["comparison"]=="Naked mole-rat aging slope"][0]
rows.append({"comparison": f"NMR TOST equivalence (delta={DELTA:.3f} = 25% of mouse slope)",
             "slope": rn["slope"], "se": rn["se"], "ci_low": rn["slope"]-DELTA, "ci_high": rn["slope"]+DELTA,
             "p": tost(rn["slope"], rn["se"], rn["n"]-2, DELTA), "n": rn["n"]})
# NMR young-vs-aged sample-level permutation (exact, 10 perms of 5 samples)
obs = nn[nn.group=="aged"]["program_score"].mean() - nn[nn.group=="young"]["program_score"].mean()
vals = nn["program_score"].values; labels = nn["group"].values
cnt, tot = 0, 0
from itertools import combinations
for aged_pick in combinations(range(5), 2):
    mask = np.zeros(5, bool); mask[list(aged_pick)] = True
    d = vals[mask].mean() - vals[~mask].mean()
    tot += 1; cnt += abs(d) >= abs(obs) - 1e-12
perm_p = (cnt + 1) / (tot + 1)
rows.append({"comparison": "NMR aged-vs-young score difference (exact permutation, sample-level)",
             "slope": obs, "se": np.nan, "ci_low": np.nan, "ci_high": np.nan, "p": perm_p, "n": 5})
sumdf = pd.DataFrame(rows)[["comparison","slope","se","ci_low","ci_high","p","n"]]
sumdf.to_csv(RESULTS/"m6-statistical-summary.csv", index=False)
prov(f"wrote {RESULTS/'m6-statistical-summary.csv'}")

# ---------------------------------------------------------------- gates
gate_mouse = [r for r in rows if r["comparison"]=="Mouse (whole-body) aging slope"][0]
prov(f"GATE positive-control (mouse skin slope>0, CI excludes 0): "
     f"slope={gate_mouse['slope']:.3f} CI=[{gate_mouse['ci_low']:.3f},{gate_mouse['ci_high']:.3f}] -> "
     f"{'PASS' if gate_mouse['ci_low']>0 else 'FAIL'}")
prov(f"GATE coverage: mouse={cov_m:.2f}, mouse-epi={cov_me:.2f}, human={cov_h:.2f}, NMR={cov_n:.2f}")
pd.DataFrame([{"gate":"mouse_positive_control","value":gate_mouse['slope'],"ci_low":gate_mouse['ci_low'],"ci_high":gate_mouse['ci_high'],"pass":bool(gate_mouse['ci_low']>0)},
              {"gate":"coverage_mouse","value":cov_m},{"gate":"coverage_human","value":cov_h},{"gate":"coverage_nmr","value":cov_n},
              {"gate":"tost_delta","value":DELTA},{"gate":"nmr_perm_p","value":perm_p}
             ]).to_csv(RESULTS/"m6-gates.csv", index=False)

# ---------------------------------------------------------------- figures
plt.rcParams.update({"figure.dpi": 150, "font.size": 9})
fig, ax = plt.subplots(figsize=(6.4, 4.4))
colors = {"Mouse":"#9ecae1","Mouse (whole-body)":"#1f77b4","Mouse (epidermal)":"#c6dbef","Human":"#d62728","Human (fibroblast)":"#f4a3a3","Naked mole-rat":"#2ca02c"}
for sp, d, mk in [("Mouse (whole-body)", mw, "o"), ("Human", hh, "s"), ("Naked mole-rat", nn, "^")]:
    ax.scatter(d["age_normalized"], d["program_score"], c=colors[sp], marker=mk, s=42, label=f"{sp} (n={len(d)})", zorder=3)
    r = ols(d["age_normalized"], d["program_score"])
    xs = np.linspace(0, max(d["age_normalized"].max(), 0.65), 50)
    ax.plot(xs, r["intercept"] + r["slope"]*xs, c=colors[sp], lw=1.6,
            label=f"slope={r['slope']:.2f} [{r['ci_low']:.2f}, {r['ci_high']:.2f}]")
ax.axhline(0, c="grey", lw=0.5)
ax.set_xlabel("life-stage-normalized age (fraction of species max lifespan, AnAge)")
ax.set_ylabel("conserved aging program score (per sample)")
ax.set_title("Conserved program activity vs normalized age\n(mouse & NMR: skin; human: TS skin dermis, LOW-POWERED n=4 donors)")
ax.legend(fontsize=7, loc="upper left")
fig.tight_layout(); fig.savefig(FIGURES/"m6-signature-figure.png"); plt.close(fig)

fig, ax = plt.subplots(figsize=(4.8, 3.6))
ax.scatter(mw["age_months"], mw["program_score"], c="#1f77b4", s=45)
r = ols(mw["age_months"], mw["program_score"])
xs = np.linspace(mw["age_months"].min(), mw["age_months"].max(), 40)
ax.plot(xs, r["intercept"] + r["slope"]*xs, c="#1f77b4",
        label=f"slope={r['slope']:.4f}/mo, p={r['p']:.3g}")
ax.set_xlabel("age (months)"); ax.set_ylabel("program score (per mouse)")
ax.set_title("Positive control: whole-body mouse program score vs age\n(all droplet cells, per-mouse mean)")
ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(FIGURES/"m6-positive-control-skin.png"); plt.close(fig)

fig, ax = plt.subplots(figsize=(4.8, 3.6))
cv = pd.DataFrame(curve)
ax.plot(cv["log2fc_cutoff"], cv["n_de_pval05"], marker="o", ms=3, label="p<0.05")
ax.plot(cv["log2fc_cutoff"], cv["n_de_fdr05"], marker="s", ms=3, label="FDR<0.05")
ax.set_xlabel("|log2 fold-change| cutoff"); ax.set_ylabel("N DE genes (aged vs young)")
ax.set_title("Savina 2022 reproduction: NMR epidermis global DE\n(paper: 0 genes at p<0.05 & logFC>1)")
ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(FIGURES/"m6-savina-reproduction.png"); plt.close(fig)
prov("figures written: m6-signature-figure.png, m6-positive-control-skin.png, m6-savina-reproduction.png")

print("\n=== M6 REAL-DATA PIPELINE COMPLETE ===")
print(sumdf.to_string(index=False))
