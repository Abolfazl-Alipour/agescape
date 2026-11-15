"""Dataset 10: Human Retina Cell Atlas (HRCA), Li et al. Nat Genet 2026.

Primary source verification:
- Paper: Nature Genetics (2026), DOI 10.1038/s41588-025-02454-1, PMID 41578023, PMCID PMC13050529.
- CELLxGENE collection 4c6eaf5c-6d57-4c76-b1e9-60df8c655f1e ("Single cell atlas of the human retina"),
  queried through the CELLxGENE curation API. The collection holds 8 datasets; two are all-cell objects:
    * be41a86a-b606-4b1c-8055-32f334898775 scRNA-seq all cells (265,767 cells, 20 donors).
    * d6505c89-c43d-4c28-8c4f-7351a5fd5528 snRNA-seq all cells (3,177,310 nuclei, 104 donor entries).
  The other six are cell-class subsets of these and are not re-analyzed.
- Donor age, sex, region, enrichment, and disease come from the deposited per-cell metadata
  (donor_age, development_stage, sex, tissue, suspension_enriched_cell_types, reported_diseases).

Unbiased sampling rule: only cells or nuclei with suspension_enriched_cell_types == "na" (and, for
scRNA-seq, suspension_depleted_cell_types == "na") are used. NeuN-enriched ("neuron") nuclei are excluded
from all analyses. Pooled multi-donor libraries and donors without an exact age are excluded.

Raw integer counts are read from raw/X by streaming row chunks with h5py (no full matrix in memory),
summed to donor x region-group x major class keys, and then to donor-level pseudobulk.
"""

from __future__ import annotations

import argparse
import gc
import os
from pathlib import Path
import sys

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from anndata.io import read_elem  # noqa: E402

from src.followup_organs import (  # noqa: E402
    DATA_TMP_ROOT,
    download_with_guard,
    finalize_pseudobulk_group,
    load_p1_clean_human,
    run_family_a_clr,
    run_family_b_genes_and_concordance,
    run_family_c_p1_signature,
    write_clean_csv,
    write_clean_markdown,
)

NULL_ITERATIONS = 2000
NULL_SEED = 202610
SC_PATH = DATA_TMP_ROOT / "shared_hrca_yan" / "hrca.h5ad"
SN_URL = "https://datasets.cellxgene.cziscience.com/25f1558d-c1f8-49f8-addc-b778d8cbd4aa.h5ad"
SN_PATH = DATA_TMP_ROOT / "organ10_hrca" / "hrca_snrna_all.h5ad"
SN_OBS_PARQUET = DATA_TMP_ROOT / "hrca_sn_obs.parquet"
CACHE_DIR = DATA_TMP_ROOT / "organ10_hrca"
CLASSES = ["Rod", "Cone", "BC", "AC", "HC", "RGC", "MG", "Astrocyte", "Microglia", "RPE"]
MACULAR = {"macula lutea", "macula lutea proper", "fovea centralis"}
PREFIX = "followup-organ-10-eye-hrca"


def _read_obs(h: h5py.File, cols: list[str]) -> pd.DataFrame:
    return pd.DataFrame({c: np.asarray(read_elem(h["obs"][c])).astype(str) for c in cols})


def _raw_symbols(h: h5py.File) -> np.ndarray:
    var = h["raw"]["var"]
    if "feature_name" in var:
        return np.asarray(read_elem(var["feature_name"])).astype(str)
    return np.asarray(read_elem(h["var"]["feature_name"])).astype(str)


def _stream_sum(h: h5py.File, rows: np.ndarray, keys: np.ndarray, n_keys: int, chunk: int = 40000) -> np.ndarray:
    """Sum raw/X rows (CSR) into dense [n_keys, n_features] float64, reading row chunks."""
    X = h["raw"]["X"]
    enc = X.attrs.get("encoding-type", "")
    if isinstance(enc, bytes):
        enc = enc.decode()
    if enc != "csr_matrix":
        raise RuntimeError(f"raw/X encoding {enc!r} is not csr_matrix")
    n_feat = int(X.attrs["shape"][1])
    indptr = X["indptr"][:]
    out = np.zeros((n_keys, n_feat), dtype=np.float64)
    order = np.argsort(rows)
    rows = rows[order]
    keys = keys[order]
    i = 0
    n_done = 0
    while i < len(rows):
        r0 = rows[i]
        j = np.searchsorted(rows, r0 + chunk, side="left")
        r1 = rows[j - 1] + 1
        a, b = int(indptr[r0]), int(indptr[r1])
        data = X["data"][a:b].astype(np.float64)
        if data.size and not np.all(np.mod(data, 1.0) == 0.0):
            raise ValueError("raw/X contains non-integer values")
        idx = X["indices"][a:b]
        ip = indptr[r0 : r1 + 1] - a
        block = sp.csr_matrix((data, idx, ip), shape=(r1 - r0, n_feat))
        local = rows[i:j] - r0
        sel = block[local, :]
        ind = sp.csr_matrix(
            (np.ones(len(local)), (keys[i:j], np.arange(len(local)))), shape=(n_keys, len(local))
        )
        out += (ind @ sel).toarray()
        n_done += len(local)
        i = j
        del data, idx, block, sel
        gc.collect()
        print(f"    streamed {n_done:,} / {len(rows):,}", flush=True)
    return out


def _collapse_symbols(mat: np.ndarray, symbols: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    uniq, inv = np.unique(symbols, return_inverse=True)
    proj = sp.csr_matrix((np.ones(len(symbols)), (np.arange(len(symbols)), inv)), shape=(len(symbols), len(uniq)))
    return np.asarray(mat @ proj), uniq.astype(object)


def _parse_age(s: str) -> float:
    s = s.strip()
    if s.startswith("pooled") or s in {"", "nan", "unknown"}:
        return np.nan
    if s.startswith(">89"):
        return 90.0
    try:
        return float(s.replace(" years", ""))
    except ValueError:
        return np.nan


def _repair_sn_study_from_obs(kdf: pd.DataFrame) -> pd.DataFrame:
    """Repair cached per-key study_name using (donor_id, regiongroup) from hrca_sn_obs.parquet."""
    if not SN_OBS_PARQUET.exists():
        return kdf
    obs = pd.read_parquet(
        SN_OBS_PARQUET,
        columns=["donor_id", "tissue", "study_name", "suspension_enriched_cell_types"],
    )
    unb = obs[obs["suspension_enriched_cell_types"].astype(str) == "na"].copy()
    unb["regiongroup"] = np.where(
        unb["tissue"].astype(str).isin(MACULAR),
        "macular",
        np.where(unb["tissue"].astype(str) == "peripheral region of retina", "peripheral", "other"),
    )
    study_map = (
        unb.groupby(["donor_id", "regiongroup"], observed=True)["study_name"]
        .agg(lambda s: "+".join(sorted(set(s.astype(str)))))
        .to_dict()
    )
    kdf = kdf.copy()
    kdf["study"] = [
        study_map.get((str(d), str(r)), str(st))
        for d, r, st in zip(kdf["donor_id"], kdf["regiongroup"], kdf["study"])
    ]
    return kdf


def _aggregate(path: Path, assay_label: str) -> dict:
    """Return per (donor, regiongroup, class) sums and cell counts, plus donor metadata, for unbiased cells."""
    cache = CACHE_DIR / f"{assay_label}_agg.npz"
    meta_csv = CACHE_DIR / f"{assay_label}_cells_meta.csv"
    if cache.exists() and meta_csv.exists():
        z = np.load(cache, allow_pickle=True)
        kdf = pd.read_csv(meta_csv)
        if assay_label == "sn":
            fixed = _repair_sn_study_from_obs(kdf)
            if not fixed["study"].equals(kdf["study"]):
                fixed.to_csv(meta_csv, index=False)
            kdf = fixed
        return {"sums": z["sums"], "genes": z["genes"], "keys": kdf}
    with h5py.File(path, "r") as h:
        cols = ["donor_id", "donor_age", "development_stage", "sex", "tissue", "study_name",
                "suspension_enriched_cell_types", "majorclass", "reported_diseases"]
        if "suspension_depleted_cell_types" in h["obs"]:
            cols.append("suspension_depleted_cell_types")
        obs = _read_obs(h, cols)
        if "suspension_depleted_cell_types" not in obs:
            obs["suspension_depleted_cell_types"] = "na"
        obs["age"] = obs["donor_age"].map(_parse_age)
        keep = (
            (obs["suspension_enriched_cell_types"] == "na")
            & (obs["suspension_depleted_cell_types"] == "na")
            & obs["age"].notna()
            & ~obs["donor_id"].str.startswith("pooled")
            & obs["majorclass"].isin(CLASSES)
        )
        if assay_label == "sc":
            # Roska_EGAS00001004561 donors carry only "adult stage" and are already excluded by the age filter.
            pass
        obs["regiongroup"] = np.where(obs["tissue"].isin(MACULAR), "macular", np.where(obs["tissue"] == "peripheral region of retina", "peripheral", "other"))
        sub = obs[keep].copy()
        sub["key"] = sub["donor_id"] + "||" + sub["regiongroup"] + "||" + sub["majorclass"]
        key_codes, key_names = pd.factorize(sub["key"])
        rows = np.where(keep.to_numpy())[0]
        print(f"  {assay_label}: {len(rows):,} unbiased cells, {sub['donor_id'].nunique()} donors, {len(key_names)} keys", flush=True)
        sums = _stream_sum(h, rows, key_codes.astype(np.int64), len(key_names))
        symbols = _raw_symbols(h)
    sums, genes = _collapse_symbols(sums, symbols)
    kdf = pd.DataFrame({"key": key_names})
    parts = kdf["key"].str.split(r"\|\|", expand=True)
    kdf["donor_id"], kdf["regiongroup"], kdf["majorclass"] = parts[0], parts[1], parts[2]
    kdf["n_cells"] = np.bincount(key_codes, minlength=len(key_names))
    dmeta = sub.groupby("donor_id").agg(
        age=("age", "first"), donor_age_raw=("donor_age", "first"),
        sex=("sex", "first"), disease=("reported_diseases", "first"),
    )
    kdf = kdf.join(dmeta, on="donor_id")
    smeta = sub.groupby(["donor_id", "regiongroup"])["study_name"].agg(
        lambda s: "+".join(sorted(set(s.astype(str))))
    ).rename("study")
    kdf = kdf.join(smeta, on=["donor_id", "regiongroup"])
    np.savez(cache, sums=sums, genes=genes)
    kdf.to_csv(meta_csv, index=False)
    return {"sums": sums, "genes": genes, "keys": kdf}


def _donor_meta(keys: pd.DataFrame, donors: list[str], studies_ref: str) -> pd.DataFrame:
    present = [d_id for d_id in donors if d_id in set(keys["donor_id"])]
    d = keys.groupby("donor_id").agg(
        study=("study", lambda s: "+".join(sorted(set(s.astype(str))))),
        age=("age", "first"),
        donor_age_raw=("donor_age_raw", "first"),
        sex=("sex", "first"),
        disease=("disease", "first"),
        n_cells=("n_cells", "sum"),
    ).loc[present]
    per = keys[keys["regiongroup"] == "peripheral"].groupby("donor_id")["n_cells"].sum()
    d["frac_peripheral"] = (per.reindex(d.index).fillna(0) / d["n_cells"]).astype(float)
    d["age_per_decade"] = d["age"] / 10.0
    d["sex_male"] = (d["sex"] == "male").astype(float)
    d["disease_any"] = (d["disease"].astype(str) != "none").astype(float)
    d["age_topcoded"] = d["donor_age_raw"].astype(str).str.startswith(">89")
    ref = studies_ref if studies_ref in set(d["study"]) else str(d["study"].value_counts().idxmax())
    for s in sorted(d["study"].unique()):
        if s != ref:
            d[f"study_{s}"] = (d["study"] == s).astype(float)
    d.index.name = None
    d.insert(0, "donor_id", d.index)
    return d.sort_values(["age", "donor_id"])


def _run_suite(name, agg, donors, regions, covs, min_cells, min_donors, min_df, studies_ref):
    keys = agg["keys"]
    k = keys[keys["donor_id"].isin(donors) & keys["regiongroup"].isin(regions)]
    meta = _donor_meta(k, donors, studies_ref)
    study_covs = tuple(c for c in meta.columns if c.startswith("study_"))
    covs_full = tuple(c if c != "STUDY" else None for c in covs)
    cov_list: list[str] = []
    for c in covs_full:
        if c is None:
            cov_list.extend(study_covs)
        else:
            cov_list.append(c)
    covariates = tuple(cov_list)

    cc = k.pivot_table(index="donor_id", columns="majorclass", values="n_cells", aggfunc="sum", fill_value=0)
    cc = cc.reindex(index=meta.index, columns=[c for c in CLASSES if c in cc.columns], fill_value=0)
    fam_a = run_family_a_clr(
        cc, meta, age_col="age_per_decade",
        candidate_covariates=tuple(c for c in covariates if c != "log_lib"),
        min_donors_present=min(min_donors, len(meta)), min_total_cells=20, min_residual_df=min_df,
    )
    if not fam_a.empty:
        fam_a.insert(0, "analysis", name)

    genes = agg["genes"]
    expr_by, meta_by = {}, {}
    for ct in cc.columns:
        eligible = [d for d in meta.index if int(cc.loc[d, ct]) >= min_cells]
        if len(eligible) < min_donors:
            continue
        cols = {}
        for d in eligible:
            idx = k.index[(k["donor_id"] == d) & (k["majorclass"] == ct)]
            cols[d] = agg["sums"][idx].sum(axis=0)
        raw_df = pd.DataFrame(cols, index=genes)
        e, m, _ = finalize_pseudobulk_group(raw_df, meta.loc[eligible], cc.loc[eligible, ct])
        expr_by[ct], meta_by[ct] = e, m

    p1 = load_p1_clean_human(REPO_ROOT)
    if len(expr_by) > 0:
        fb_genes, fb_sum, fb_conc = run_family_b_genes_and_concordance(
            expr_by, meta_by, age_col="age_per_decade", candidate_covariates=covariates,
            min_donors=min_donors, min_cpm=1.0, min_expr_frac=0.5, min_residual_df=min_df,
            is_unbiased_single_cell=len(expr_by) >= 2, concordance_perms=1000, seed=NULL_SEED,
        )
        fam_c = run_family_c_p1_signature(
            expr_by, meta_by, p1_program=p1, age_col="age_per_decade", candidate_covariates=covariates,
            min_donors=min_donors, min_cpm=1.0, min_expr_frac=0.5, min_residual_df=min_df,
            iterations=NULL_ITERATIONS, seed=NULL_SEED,
        )
    else:
        fb_genes, fb_sum, fb_conc, fam_c = pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    top = fb_genes.head(50).copy() if not fb_genes.empty else pd.DataFrame()
    for df in (top, fb_sum, fb_conc, fam_c):
        if not df.empty:
            df.insert(0, "analysis", name)
    meta_out = meta.copy()
    meta_out.insert(0, "analysis", name)
    print(f"  suite {name}: n={len(meta)} groups={len(expr_by)} FamB_global_hits={int(fb_sum['n_sig_global_q05'].sum()) if not fb_sum.empty else 0}", flush=True)
    return {"fam_a": fam_a, "top": top, "sum": fb_sum, "conc": fb_conc, "fam_c": fam_c, "meta": meta_out}


def _clean_text(s: Any) -> str:
    return str(s).replace("\u2014", " - ").replace("\u2013", "-")


def _build_donor_audit() -> pd.DataFrame:
    """Build donor-level inclusion and exclusion audit across snRNA (104) and scRNA (20) objects."""
    rows: list[dict[str, Any]] = []
    if SN_OBS_PARQUET.exists():
        sn_obs = pd.read_parquet(
            SN_OBS_PARQUET,
            columns=[
                "donor_id", "donor_age", "sex", "tissue", "study_name",
                "suspension_enriched_cell_types", "reported_diseases",
            ],
        )
        sn_obs["suspension_depleted_cell_types"] = "na"
        assay_frames = [("snRNA", sn_obs)]
    else:
        assay_frames = []
    if SC_PATH.exists():
        with h5py.File(SC_PATH, "r") as h:
            cols = [
                "donor_id", "donor_age", "sex", "tissue", "study_name",
                "suspension_enriched_cell_types", "suspension_depleted_cell_types", "reported_diseases",
            ]
            sc_obs = _read_obs(h, [c for c in cols if c in h["obs"]])
            for c in cols:
                if c not in sc_obs:
                    sc_obs[c] = "na"
            assay_frames.append(("scRNA", sc_obs))

    for assay_label, obs in assay_frames:
        obs = obs.copy()
        obs["regiongroup"] = np.where(
            obs["tissue"].astype(str).isin(MACULAR),
            "macular",
            np.where(obs["tissue"].astype(str) == "peripheral region of retina", "peripheral", "other"),
        )
        obs["is_unbiased"] = (
            (obs["suspension_enriched_cell_types"].astype(str) == "na")
            & (obs["suspension_depleted_cell_types"].astype(str) == "na")
        )
        for donor_id, grp in obs.groupby("donor_id", observed=True):
            did = _clean_text(donor_id)
            age_raw = _clean_text(grp["donor_age"].iloc[0])
            age_num = _parse_age(age_raw)
            sex = _clean_text(grp["sex"].iloc[0])
            disease = _clean_text(grp["reported_diseases"].iloc[0])
            studies = "+".join(sorted(set(_clean_text(x) for x in grp["study_name"])))
            unb = grp[grp["is_unbiased"]]
            unb_studies = "+".join(sorted(set(_clean_text(x) for x in unb["study_name"]))) if not unb.empty else ""
            regions = "+".join(sorted(set(_clean_text(x) for x in unb["regiongroup"]))) if not unb.empty else "none"
            n_total = int(len(grp))
            n_unb = int(len(unb))
            n_biased = n_total - n_unb
            if did.startswith("pooled"):
                status = "excluded_pooled"
                note = "multi-donor pool; donor unit ambiguous"
            elif n_unb == 0:
                if assay_label == "snRNA":
                    status = "excluded_neun_enriched_only"
                    note = "all nuclei NeuN-enriched (neuron)"
                else:
                    status = "excluded_rod_depleted_rgc_enriched"
                    note = "no unbiased na/na cells"
            elif not np.isfinite(age_num):
                status = "excluded_age_nonexact"
                note = f"non-exact donor_age ({age_raw})"
            elif assay_label == "snRNA":
                if age_num < 18.0:
                    status = "included_sn_pediatric_sensitivity"
                    note = "pediatric <18 yr; excluded from Sn_Adult suites, included in Sn_AllAges"
                elif age_raw.startswith(">89"):
                    status = "included_sn_adult_topcoded"
                    note = "topcoded >89 yr coded as 90; excluded in Sn_Adult_ExactAge"
                else:
                    status = "included_sn_adult_exact"
                    note = "exact-age disease-free unbiased adult donor"
            else:
                if "Sanes_GSE148077" in studies:
                    status = "included_sc_sanes_organ07_dup"
                    note = f"duplicates Organ 07 Sanes GSE148077 fovea; disease={disease}"
                elif disease != "none":
                    status = "included_sc_external_diseased"
                    note = f"external scRNA donor with ocular/systemic disease={disease}"
                else:
                    status = "included_sc_external_healthy"
                    note = "external disease-free scRNA donor"
            rows.append({
                "assay": assay_label,
                "donor_id": did,
                "studies_all": studies,
                "studies_unbiased": unb_studies,
                "regions_unbiased": regions,
                "donor_age_raw": age_raw,
                "age_numeric": age_num if np.isfinite(age_num) else "",
                "sex": sex,
                "reported_diseases": disease,
                "total_cells": n_total,
                "unbiased_cells": n_unb,
                "enriched_or_depleted_cells": n_biased,
                "status": status,
                "note": note,
            })
    return pd.DataFrame(rows).sort_values(["assay", "status", "donor_id"])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["aggregate-sc", "aggregate-sn", "analyze"], required=True)
    args = ap.parse_args()
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    if args.stage == "aggregate-sc":
        _aggregate(SC_PATH, "sc")
        return
    if args.stage == "aggregate-sn":
        if not SN_PATH.exists():
            download_with_guard(SN_URL, SN_PATH, max_file_gb=40.0)
        _aggregate(SN_PATH, "sn")
        return

    sc = _aggregate(SC_PATH, "sc")
    sn = _aggregate(SN_PATH, "sn")
    sck = sc["keys"]
    # scRNA: Sanes fovea only (peripheral Sanes libraries were rod-depleted / RGC-enriched and are removed by
    # the enrichment filter; keep macular region only for Sanes to match Organ 07), all regions for others.
    sck = sck[~((sck["study"] == "Sanes_GSE148077") & (sck["regiongroup"] != "macular"))].copy()
    sc = {**sc, "keys": sck}
    sc_all = sorted(sck["donor_id"].unique())
    sc_ext = sorted(sck.loc[sck["study"] != "Sanes_GSE148077", "donor_id"].unique())
    sc_ext_nodis = sorted(
        sck.loc[(sck["study"] != "Sanes_GSE148077") & (sck["disease"] == "none"), "donor_id"].unique()
    )

    snk = sn["keys"]
    dages = snk.groupby("donor_id")["age"].first()
    dtop = snk.groupby("donor_id")["donor_age_raw"].first().str.startswith(">89")
    sn_adult = sorted(dages.index[dages >= 18])
    sn_all = sorted(dages.index)
    sn_adult_exact = sorted(dages.index[(dages >= 18) & ~dtop])
    sn_anc_donors = set(snk.loc[snk["study"] == "Chen_ancestry", "donor_id"].unique())
    sn_anc_adult = sorted(d for d in sn_adult if d in sn_anc_donors)

    suites = [
        _run_suite("Sn_Adult_MacularOnly", sn, sn_adult, ["macular"],
                   ("STUDY", "sex_male", "log_lib"), 10, 15, 5, "Chen_ancestry"),
        _run_suite("Sn_Adult_PeripheralOnly", sn, sn_adult, ["peripheral"],
                   ("STUDY", "sex_male", "log_lib"), 10, 15, 5, "Chen_ancestry"),
        _run_suite("Sn_Adult_Unenriched", sn, sn_adult, ["macular", "peripheral"],
                   ("STUDY", "sex_male", "frac_peripheral", "log_lib"), 10, 20, 5, "Chen_ancestry"),
        _run_suite("Sn_Adult_ExactAge", sn, sn_adult_exact, ["macular", "peripheral"],
                   ("STUDY", "sex_male", "frac_peripheral", "log_lib"), 10, 20, 5, "Chen_ancestry"),
        _run_suite("Sn_ChenAncestry_Adult", sn, sn_anc_adult, ["macular", "peripheral"],
                   ("sex_male", "frac_peripheral", "log_lib"), 10, 15, 5, "Chen_ancestry"),
        _run_suite("Sn_AllAges_Unenriched", sn, sn_all, ["macular", "peripheral"],
                   ("STUDY", "sex_male", "frac_peripheral", "log_lib"), 10, 20, 5, "Chen_ancestry"),
        _run_suite("Sc15_StudyAdj", sc, sc_all, ["macular", "peripheral", "other"],
                   ("STUDY", "disease_any", "sex_male", "log_lib"), 10, 6, 3, "Sanes_GSE148077"),
        _run_suite("Sc9_External_StudyAdj", sc, sc_ext, ["macular", "peripheral", "other"],
                   ("STUDY", "disease_any", "sex_male", "log_lib"), 10, 6, 3, "Hafler_GSE137537"),
        _run_suite("Sc5_External_NoDisease", sc, sc_ext_nodis, ["macular", "peripheral", "other"],
                   ("STUDY", "sex_male", "log_lib"), 10, 5, 2, "Hafler_GSE137537"),
    ]

    def cat(key):
        frames = [s[key] for s in suites if not s[key].empty]
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    out = {k: cat(k) for k in ["fam_a", "top", "sum", "conc", "fam_c", "meta"]}
    rd = REPO_ROOT / "results"
    write_clean_csv(out["meta"], rd / f"{PREFIX}-donors.csv")
    write_clean_csv(out["fam_a"], rd / f"{PREFIX}-fractions.csv")
    write_clean_csv(out["sum"], rd / f"{PREFIX}-gene-summary.csv")
    write_clean_csv(out["top"], rd / f"{PREFIX}-top-genes.csv")
    write_clean_csv(out["conc"], rd / f"{PREFIX}-concordance.csv")
    write_clean_csv(out["fam_c"], rd / f"{PREFIX}-p1-signature.csv")
    audit = _build_donor_audit()
    if not audit.empty:
        write_clean_csv(audit, rd / f"{PREFIX}-donor-audit.csv")
    print("analysis done", flush=True)


if __name__ == "__main__":
    main()
