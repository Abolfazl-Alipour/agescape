#!/usr/bin/env python3
"""11pt Block 7: within-type death/stress scoring and narrowing localization."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import time
from typing import Any

for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_key, "1")

import numpy as np
import pandas as pd
from scipy import sparse as sp
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.run_followup_between_subject_gao_yazar import _prune_covariate_matrix
from src.followup_organs import (
    counts_to_log2cpm,
    run_family_b_genes_and_concordance,
    write_clean_csv,
    write_clean_markdown,
)

DATA_ROOT = REPO_ROOT.parent / f"{REPO_ROOT.name}-data" / "tmp_organs"
BLOCK1_CACHE = DATA_ROOT / "block1_cache"
CACHE_ROOT = DATA_ROOT / "block7_cache"
RESULTS_DIR = REPO_ROOT / "results"
REPORTS_DIR = REPO_ROOT / "reports"
REPORT_PATH = REPORTS_DIR / "followup-11pt-block7-survivor-within-type.md"
DEATH_PATH = RESULTS_DIR / "followup-11pt-block7-death-stress.csv"
REMOVAL_PATH = RESULTS_DIR / "followup-11pt-block7-removal-sensitivity.csv"
NARROWING_PATH = RESULTS_DIR / "followup-11pt-block7-narrowing-localization.csv"

MASTER_SEED = 20261117
DATASET_CODE = {"gao": 1, "yazar": 2}
MIN_CELLS = 20
MIN_DONORS = 6
MIN_RESID_DF = {"gao": 5, "yazar": 10}
BATCH_REF = {"gao": "ENCODE v4 (Snyder)", "yazar": None}
MDE80_SMALL_SD = 0.50
GENE_SOURCE = (
    "MSigDB HALLMARK_APOPTOSIS M5902 and HALLMARK_UNFOLDED_PROTEIN_RESPONSE M5922; "
    "PMID 26771021; DOI 10.1016/j.cels.2015.12.004"
)

APOPTOSIS = """ADD1 AIFM3 ANKH ANXA1 APP ATF3 AVPR1A BAX BCAP31 BCL10 BCL2L1 BCL2L10 BCL2L11 BCL2L2 BGN BID BIK BIRC3 BMF BMP2 BNIP3L BRCA1 BTG2 BTG3 CASP1 CASP2 CASP3 CASP4 CASP6 CASP7 CASP8 CASP9 CAV1 CCNA1 CCND1 CCND2 CD14 CD2 CD38 CD44 CD69 CDC25B CDK2 CDKN1A CDKN1B CFLAR CLU CREBBP CTH CTNNB1 CYLD DAP DAP3 DCN DDIT3 DFFA DIABLO DNAJA1 DNAJC3 DNM1L DPYD EBP EGR3 EMP1 ENO2 ERBB2 ERBB3 EREG ETF1 F2 F2R FAS FASLG FDXR FEZ1 GADD45A GADD45B GCH1 GNA15 GPX1 GPX3 GPX4 GSN GSR GSTM1 GUCY2D H1-0 HGF HMGB2 HMOX1 HSPB1 IER3 IFITM3 IFNB1 IFNGR1 IGF2R IGFBP6 IL18 IL1A IL1B IL6 IRF1 ISG20 JUN KRT18 LEF1 LGALS3 LMNA LUM MADD MCL1 MGMT MMP2 NEDD9 NEFH PAK1 PDCD4 PDGFRB PEA15 PLAT PLCB2 PLPPR4 PMAIP1 PPP2R5B PPP3R1 PPT1 PRF1 PSEN1 PSEN2 PTK2 RARA RELA RETSAT RHOB RHOT2 RNASEL ROCK1 SAT1 SATB1 SC5D SLC20A1 SMAD7 SOD1 SOD2 SPTAN1 SQSTM1 TAP1 TGFB2 TGFBR3 TIMP1 TIMP2 TIMP3 TNF TNFRSF12A TNFSF10 TOP2A TSPO TXNIP VDAC2 WEE1 XIAP""".split()
UPR = """ALDH18A1 ARFGAP1 ASNS ATF3 ATF4 ATF6 ATP6V0D1 BAG3 BANF1 CALR CCL2 CEBPB CEBPG CHAC1 CKS1B CNOT2 CNOT4 CNOT6 CXXC1 DCP1A DCP2 DCTN1 DDIT4 DDX10 DKC1 DNAJA4 DNAJB9 DNAJC3 EDC4 EDEM1 EEF2 EIF2AK3 EIF2S1 EIF4A1 EIF4A2 EIF4A3 EIF4E EIF4EBP1 EIF4G1 ERN1 ERO1A EXOC2 EXOSC1 EXOSC10 EXOSC2 EXOSC4 EXOSC5 EXOSC9 FKBP14 FUS GEMIN4 GOSR2 H2AX HERPUD1 HSP90B1 HSPA5 HSPA9 HYOU1 IARS1 IFIT1 IGFBP1 IMP3 KDELR3 KHSRP KIF5B LSM1 LSM4 MTHFD2 MTREX NABP1 NFYA NFYB NHP2 NOLC1 NOP14 NOP56 NPM1 PAIP1 PARN PDIA5 PDIA6 POP4 PREB PSAT1 RPS14 RRP9 SDAD1 SEC11A SEC31A SERP1 SHC1 SKIC3 SLC1A4 SLC30A5 SLC7A5 SPCS1 SPCS3 SRPRA SRPRB SSR1 STC2 TARS1 TATDN2 TSPYL2 TUBB2A VEGFA WFS1 WIPI1 XBP1 XPOT YIF1A YWHAZ ZBTB17""".split()
SOURCE_UNION = sorted(set(APOPTOSIS) | set(UPR))
MASKED_LIST_GENES = ["RPS14"]
SCORING_GENES = sorted(set(SOURCE_UNION) - set(MASKED_LIST_GENES))

assert len(APOPTOSIS) == 161
assert len(UPR) == 113
assert len(SOURCE_UNION) == 272
assert len(SCORING_GENES) == 271


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def frozen_spec() -> str:
    apoptosis = ", ".join(APOPTOSIS)
    upr = ", ".join(UPR)
    union = ", ".join(SCORING_GENES)
    return f"""# 11pt Follow-Up Block 7: Survivor-Bias Within-Type Tests

## Frozen protocol and gene list

This protocol and gene list were frozen before any Block 7 score, model, or result was computed. The source reference lookup did not download expression data. No gene symbol, threshold, model, family, or verdict rule may be edited after seeing results.

### Inputs and exact source files

- Gao h5ad: `/home/reza/projects/agescape-revision-data/tmp_organs/organ17_gao_heart/GSE290367_integrated_RNA_adata.h5ad`.
- Yazar h5ad: `/home/reza/projects/agescape-revision-data/tmp_organs/hcata/yazar_1e44db10.h5ad`, with `YAZAR_CACHE` and exclusion of rows labeled `peripheral blood mononuclear cell`.
- Donor tables: `results/followup-between-subject-gao-donors.csv` and `results/followup-between-subject-yazar-donors.csv`.
- Full-depth Part 2 spread inputs, identified by listing `results/`: `results/followup-fulldepth-gao-spread.csv` and `results/followup-fulldepth-yazar-spread.csv`.
- Paper-cohort fraction input: `results/reviewer-phase1-composition-global-bh.csv`, 431 tests, `effect` is CLR slope per year. It is recorded as required, but is not joined to Gao or Yazar because it contains no Gao or Yazar rows and no matching full-depth Part 2 spread slopes.
- Gao fraction slopes: `results/followup-organ-17-heart-gao-fractions.csv`, `subset=pooled_7study_primary`, `effect_clr` per decade.
- Yazar fraction slopes: `results/followup-hcata-pool-sensitivity-fractions.csv`, `study=yazar`, `arm=pool_fe`, `effect_clr` per decade.
- The existing Block 1 donor-sharded cache is reused. It contains raw integer counts after duplicate-symbol collapse and the human technical mask, restricted to post-mask UMI >= 1500 cells. Thus this block does not reread the large h5ad matrices and does not load them into memory.

### Frozen death/stress gene list and citation

Source: Liberzon et al. 2015, Cell Systems, The Molecular Signatures Database Hallmark Gene Set Collection, DOI `10.1016/j.cels.2015.12.004`, PMID `26771021`. Official MSigDB sets: `HALLMARK_APOPTOSIS`, systematic ID `M5902`, 161 genes, https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_APOPTOSIS.html; and `HALLMARK_UNFOLDED_PROTEIN_RESPONSE`, systematic ID `M5922`, 113 genes, https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_UNFOLDED_PROTEIN_RESPONSE.html. Their union is 272 genes because `ATF3` and `DNAJC3` occur in both.

Apoptosis source genes: {apoptosis}.

Unfolded protein response source genes: {upr}.

The shared human technical mask excludes uppercase `MT-`, `RP[L/S]`, `MRP[L/S]`, `XIST`, `TSIX`, `DDX3Y`, `EIF2S3Y`, `UTY`, and `KDM5D`. `RPS14` is the one source-list gene excluded by this mask. It is explicitly recorded here and is not silently dropped. The frozen scoring list is the remaining 271 unique symbols. Symbols absent from a dataset contribute zero counts and remain in the fixed denominator.

Final scoring symbols: {union}.

### Part A specification

- Per cell, program score is the mean within-cell proportion across the fixed 271 scoring genes: the sum of their post-mask raw counts divided by `271 * post-mask UMI`. Primary analyses use full-depth counts and no downsampling. Cells with post-mask UMI < 1500 were excluded as an eligibility filter only.
- For each dataset and cell type, donor-cell eligibility is at least 20 eligible cells. The high-score threshold is the 90th percentile, computed with numpy linear quantile interpolation, among cells from eligible donor-cell-type pairs pooled across donors. High means `score >= threshold`, including ties. The threshold and tie convention are fixed here. The same full-depth high-score set is removed in test 2.
- Test 1 outcome is each donor's fraction of high-scoring cells among that donor's eligible cells in the cell type. Primary OLS is outcome on `age_per_decade + sex_male + batch`: Gao 7-study fixed effects and Yazar 75-pool fixed effects. At least 6 donors and the repository rank and residual-df gates are required. The unadjusted sensitivity retains age and sex but removes batch. No unadjusted fallback is substituted for a failed batch-adjusted model.
- The labeled 1500 sensitivity repeats test 1 only. For each cell, aggregate source-list counts after multinomial rarefaction to 1500 UMIs have the exact marginal distribution `Binomial(original_list_count, 1500 / original_UMI)`. Independent fixed-seed draws across cells therefore simulate the exact aggregate-list score without changing any primary data. The 90th-percentile rule is recomputed on rarefied scores within each cell type. Seeds derive from `SeedSequence([20261117, 7, dataset_code, donor_idx])`, with cells sorted by original h5ad row index.
- Test 2 is run only for cell types with primary test 1 BH q < 0.05. Donor pseudobulk gene models mirror the follow-up gene tests: log2CPM, age per decade, sex, log library size, and the same study or pool fixed effects; genes pass the original expression filter of CPM >= 1 in at least half of donors, with a minimum of 3 donors. Both arms require at least 20 cells per donor after their respective cell filters and at least 6 common donors. Unremoved and removed fits use the same common donor set. The gene family is the original `q_bh_within_group` family: BH separately within each dataset, cell type, and removal arm across that arm's expressed genes. Pass means q < 0.05. The descriptive predicted change is an unremoved pass becoming a removed fail while the absolute age slope becomes smaller. Slope change is removed minus unremoved. All other pass/fail changes are also counted.
- Removing captured high-scoring cells is a lower bound. Dying cells are under-captured, so removal cannot remove all dying or stressed cells.
- Test 1 BH families are stated before correction as each `(dataset, test)` slice across evaluable non-pooled cell types, separately for primary full depth, unadjusted full depth, and rarefied 1500. Test 2 uses the gene-level family stated above. No BH correction is applied across slope changes.
- MDE80 is `(t_0.975 + t_0.80) * SE`, using residual df, on the raw outcome scale and standardized by the donor-outcome SD. A small MDE80 for test 1 is predefined as at most 0.50 outcome SD per decade, matching the repository convention.

### Part B specification

- Spread slopes are read, not recomputed, from the two exact full-depth CSVs above. Required rows are `metric=shannon`, `depth_mode=full_depth`, `analysis_role=primary_full_depth`, `model_spec=primary`, non-pooled, and non-stratified-global. These are the batch-adjusted and depth-adjusted per-cell-type primary rows needed for a cell-type join. A `stratified_global_primary` row is one cross-type summary and cannot provide per-cell-type slopes.
- Gao and Yazar are joined separately on exact cell-type keys to their batch-adjusted CLR fraction slopes above. Fraction slope negative means losing share. Spread slope is Shannon change per decade, so a decrease is negative.
- Prediction frozen before running: survivor bias predicts a positive association between fraction slope and spread slope, with more negative spread slopes concentrated in fraction-losing types.
- For each dataset, fit weighted least squares of spread slope on fraction slope across cell types with inverse variance weights `1 / SE_spread^2`. Report the association slope, 95 percent CI, two-sided p, and MDE80. Also report mean spread slope and cell-type counts for fraction-losing versus fraction-gaining types. Zero fraction slope is included in the gaining group.
- There is one test per dataset and no multiplicity correction across the two datasets. The datasets are never pooled.

### Verdict rules

- Part A, per dataset: PROVES only if test 1 passes BH correction and test 2 produces at least one predicted pass-to-fail change with absolute-slope attenuation. DOES NOT PROVE if no test 1 result passes BH and every evaluable primary null has MDE80 at most 0.50 SD per decade. INCONCLUSIVE if a null has larger MDE80, or test 1 is significant but removal changes nothing in the predicted direction, or the required removal comparison is not evaluable. The overall Part A summary is PROVES if either dataset meets PROVES, DOES NOT PROVE only if both datasets meet DOES NOT PROVE, and otherwise INCONCLUSIVE.
- Part B: PROVES only if the cross-type association is positive and significant in at least one dataset and the other dataset is not significantly opposite. DOES NOT PROVE if the association is null or opposite in both datasets. A significant positive result in one dataset with a significant opposite result in the other is INCONCLUSIVE.
- Neither Part A nor Part B alone proves survivor bias as a mechanism. The Nov 17 cross-type correlation is not rerun.

### Outputs and reproducibility

- `results/followup-11pt-block7-death-stress.csv`
- `results/followup-11pt-block7-removal-sensitivity.csv`
- `results/followup-11pt-block7-narrowing-localization.csv`
- `reports/followup-11pt-block7-survivor-within-type.md`
- `scripts/run_followup_11pt_block7_survivor_within_type.py`

All randomness uses fixed seed `20261117` and deterministic derived seeds. Donor is the unit in every regression. No permutation or bootstrap is used. No split-half technical null is required by this design. No primary analysis uses downsampling. Exactly one deep review round follows completion of the outputs.
"""


def verify_frozen_spec() -> str:
    frozen_path = CACHE_ROOT / "frozen_spec.md"
    if not frozen_path.exists():
        raise FileNotFoundError(f"Run --freeze-spec first: {frozen_path}")
    text = frozen_path.read_text(encoding="utf-8")
    if text != frozen_spec():
        raise RuntimeError("The frozen specification or gene list differs from the analysis script")
    if not REPORT_PATH.exists():
        raise FileNotFoundError(f"Missing frozen report: {REPORT_PATH}")
    report = REPORT_PATH.read_text(encoding="utf-8")
    if report != text and not report.startswith(text + "\n## Results\n"):
        raise RuntimeError("The frozen protocol portion of the report was edited")
    return text


def freeze() -> None:
    CACHE_ROOT.mkdir(parents=True, exist_ok=True)
    text = frozen_spec()
    frozen_path = CACHE_ROOT / "frozen_spec.md"
    if frozen_path.exists() and frozen_path.read_text(encoding="utf-8") != text:
        raise RuntimeError("A different Block 7 specification is already frozen")
    write_clean_markdown(text, frozen_path)
    if REPORT_PATH.exists():
        report = REPORT_PATH.read_text(encoding="utf-8")
        if report != text and not report.startswith(text + "\n## Results\n"):
            raise RuntimeError("An existing Block 7 report has a different frozen protocol")
    else:
        write_clean_markdown(text, REPORT_PATH)
    digest = hashlib.sha256(text.encode("utf-8")).hexdigest()
    (CACHE_ROOT / "frozen_spec_sha256.txt").write_text(digest + "\n", encoding="utf-8")
    log(f"Frozen Block 7 protocol and 271-gene scoring list: sha256 {digest}")


def load_info(ds: str) -> dict[str, Any]:
    path = BLOCK1_CACHE / ds / "CACHE_DONE.json"
    if not path.exists():
        raise FileNotFoundError(f"Missing required cached input: {path}")
    info = json.loads(path.read_text(encoding="utf-8"))
    if not (BLOCK1_CACHE / ds / "donors_meta.csv").exists():
        raise FileNotFoundError(f"Missing donor metadata in {path.parent}")
    return info


def donor_metadata(ds: str, info: dict[str, Any]) -> pd.DataFrame:
    donors = pd.read_csv(BLOCK1_CACHE / ds / "donors_meta.csv")
    donors["donor_id"] = donors["donor_id"].astype(str)
    donors["donor_idx"] = np.arange(len(donors), dtype=np.int32)
    donors["study_or_pool"] = (
        donors["study"].astype(str) if ds == "gao" else donors["pool_number"].astype(int).astype(str)
    )
    if list(donors["donor_id"]) != list(info["donor_ids"]):
        raise AssertionError(f"Donor order mismatch in {ds} cache")
    return donors


def fit_donor_ols(
    frame: pd.DataFrame,
    outcome: str,
    ds: str,
    adjusted: bool,
) -> dict[str, Any]:
    work = frame.dropna(subset=[outcome, "age_per_decade"]).copy()
    base: dict[str, Any] = {
        "status": "EVALUATED",
        "n_donors": int(len(work)),
        "n_cells": int(work["n_cells"].sum()) if "n_cells" in work else 0,
        "residual_df": np.nan,
        "slope_per_decade": np.nan,
        "se": np.nan,
        "ci95_low": np.nan,
        "ci95_high": np.nan,
        "pvalue": np.nan,
        "outcome_sd": np.nan,
        "mde80": np.nan,
        "mde80_sd_per_decade": np.nan,
        "covariates_active": "none",
        "covariates_omitted": "none",
    }
    if len(work) < MIN_DONORS or work["age_per_decade"].nunique() < 2:
        base["status"] = "UNTESTABLE_DONOR_FLOOR"
        return base
    batch_col = "study_or_pool" if adjusted else None
    min_df = MIN_RESID_DF[ds] if adjusted else 5
    sub, active, omitted, status = _prune_covariate_matrix(
        work,
        nonbatch_covs=["sex_male"],
        batch_col=batch_col,
        min_batch_donors=2,
        min_residual_df=min_df,
        ref_level=BATCH_REF[ds] if adjusted else None,
        require_is_comp_collinear_check="age_per_decade",
        atomic_batch_block=True,
    )
    if status != "evaluated":
        base["status"] = f"UNTESTABLE_{status.upper()}"
        base["covariates_omitted"] = ";".join(omitted)
        return base
    y = sub[outcome].to_numpy(dtype=float)
    sd = float(np.std(y, ddof=1))
    base["n_donors"] = int(len(sub))
    base["n_cells"] = int(sub["n_cells"].sum())
    base["outcome_sd"] = sd
    if not np.isfinite(sd) or sd <= 0:
        base["status"] = "UNTESTABLE_ZERO_VAR"
        return base
    X = sm.add_constant(sub[["age_per_decade", *active]].astype(float), has_constant="add")
    fit = sm.OLS(y, X).fit()
    slope = float(fit.params["age_per_decade"])
    se = float(fit.bse["age_per_decade"])
    df = int(fit.df_resid)
    t975 = float(stats.t.ppf(0.975, df))
    t80 = float(stats.t.ppf(0.80, df))
    mde = (t975 + t80) * se
    base.update(
        residual_df=df,
        slope_per_decade=slope,
        se=se,
        ci95_low=slope - t975 * se,
        ci95_high=slope + t975 * se,
        pvalue=float(fit.pvalues["age_per_decade"]),
        mde80=mde,
        mde80_sd_per_decade=mde / sd,
        covariates_active=";".join(active) if active else "none",
        covariates_omitted=";".join(omitted) if omitted else "none",
    )
    return base


def run_scores(ds: str) -> None:
    verify_frozen_spec()
    info = load_info(ds)
    donors = donor_metadata(ds, info)
    out_dir = CACHE_ROOT / ds / "scores"
    out_dir.mkdir(parents=True, exist_ok=True)
    symbols = list(info["retained_symbols"])
    symbol_to_idx = {gene: i for i, gene in enumerate(symbols)}
    score_idx = np.array([symbol_to_idx[g] for g in SCORING_GENES if g in symbol_to_idx], dtype=np.int32)
    absent = [g for g in SCORING_GENES if g not in symbol_to_idx]
    n_genes = len(symbols)
    n_shards = int(info["n_shards"])
    for shard in range(n_shards):
        out_path = out_dir / f"scores_shard{shard:02d}.npz"
        if out_path.exists():
            continue
        files = sorted((BLOCK1_CACHE / ds).glob(f"shard{shard:02d}_chunk*.npz"))
        if not files:
            continue
        mats: list[sp.csr_matrix] = []
        meta_parts: list[np.ndarray] = []
        for file in files:
            z = np.load(file)
            n_rows = len(z["donor"])
            mat = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=(n_rows, n_genes))
            mats.append(mat)
            meta_parts.append(np.column_stack([z["donor"], z["ct"], z["umi"], z["obs_idx"]]))
        mat = sp.vstack(mats, format="csr")
        meta = np.vstack(meta_parts)
        order = np.argsort(meta[:, 3], kind="mergesort")
        mat = mat[order]
        meta = meta[order]
        del mats, meta_parts
        donor = meta[:, 0].astype(np.int32)
        ct = meta[:, 1].astype(np.int32)
        umi = meta[:, 2].astype(np.int64)
        obs_idx = meta[:, 3].astype(np.int32)
        list_counts = np.rint(np.asarray(mat[:, score_idx].sum(axis=1)).ravel()).astype(np.int64)
        full_score = list_counts.astype(np.float64) / (float(len(SCORING_GENES)) * umi.astype(np.float64))
        rare_score = np.empty(len(meta), dtype=np.float64)
        for d in np.unique(donor):
            rows = np.flatnonzero(donor == d)
            rng = np.random.default_rng(np.random.SeedSequence([MASTER_SEED, 7, DATASET_CODE[ds], int(d)]))
            probability = 1500.0 / umi[rows].astype(np.float64)
            rare_counts = rng.binomial(list_counts[rows], probability)
            rare_score[rows] = rare_counts.astype(np.float64) / float(len(SCORING_GENES) * 1500)
        np.savez_compressed(
            out_path,
            donor=donor,
            ct=ct,
            umi=umi,
            obs_idx=obs_idx,
            list_counts=list_counts,
            full_score=full_score,
            rare_score=rare_score,
        )
        log(f"{ds} shard {shard + 1}/{n_shards}: scored {len(meta)} cells")
    audit = {
        "dataset": ds,
        "source_union_genes": len(SOURCE_UNION),
        "masked_source_genes": MASKED_LIST_GENES,
        "scoring_genes": len(SCORING_GENES),
        "present_scoring_genes": int(len(score_idx)),
        "absent_scoring_genes": absent,
        "n_shards": n_shards,
    }
    (out_dir / "SCORES_DONE.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    log(f"{ds}: full-depth and rarefied score cache complete")


def load_scores(ds: str) -> dict[str, np.ndarray]:
    score_dir = CACHE_ROOT / ds / "scores"
    if not (score_dir / "SCORES_DONE.json").exists():
        raise FileNotFoundError(f"Missing completed score cache for {ds}")
    files = sorted(score_dir.glob("scores_shard*.npz"))
    pieces = {key: [] for key in ("donor", "ct", "umi", "obs_idx", "list_counts", "full_score", "rare_score")}
    for file in files:
        z = np.load(file)
        for key in pieces:
            pieces[key].append(z[key])
    return {key: np.concatenate(value) for key, value in pieces.items()}


def apply_test1_bh(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out["q_bh"] = np.nan
    out["family_size"] = np.nan
    eligible = (out["status"] == "EVALUATED") & np.isfinite(out["pvalue"].astype(float))
    for _, indices in out.loc[eligible].groupby(["dataset", "test"], sort=False).groups.items():
        idx = list(indices)
        out.loc[idx, "family_size"] = len(idx)
        out.loc[idx, "q_bh"] = multipletests(out.loc[idx, "pvalue"].astype(float), method="fdr_bh")[1]
    return out


def run_test1() -> pd.DataFrame:
    verify_frozen_spec()
    rows: list[dict[str, Any]] = []
    for ds in ("gao", "yazar"):
        info = load_info(ds)
        donors = donor_metadata(ds, info)
        scores = load_scores(ds)
        for ci, cell_type in enumerate(info["cell_types"]):
            ct_mask = scores["ct"] == ci
            ct_donor = scores["donor"][ct_mask]
            ct_full = scores["full_score"][ct_mask]
            ct_rare = scores["rare_score"][ct_mask]
            all_donors, all_counts = np.unique(ct_donor, return_counts=True)
            count_map = dict(zip(all_donors.tolist(), all_counts.tolist()))
            eligible_donors = np.array(sorted(d for d, n in count_map.items() if n >= MIN_CELLS), dtype=np.int32)
            eligible_cell_mask = np.isin(ct_donor, eligible_donors)
            ed = ct_donor[eligible_cell_mask]
            ef = ct_full[eligible_cell_mask]
            er = ct_rare[eligible_cell_mask]
            base = {
                "dataset": ds,
                "cell_type": cell_type,
                "n_donors_with_cells": int(len(all_donors)),
                "n_eligible_cells": int(len(ed)),
                "n_source_genes": len(SOURCE_UNION),
                "n_scoring_genes": len(SCORING_GENES),
                "masked_list_genes": ";".join(MASKED_LIST_GENES),
                "gene_list_source": GENE_SOURCE,
            }
            donor_rows: dict[str, pd.DataFrame] = {}
            for arm, values in (("full", ef), ("rare", er)):
                if len(ed) == 0:
                    donor_frame = pd.DataFrame()
                    threshold = np.nan
                    n_high = 0
                else:
                    threshold = float(np.quantile(values, 0.90, method="linear"))
                    high = values >= threshold
                    unique_d, inverse = np.unique(ed, return_inverse=True)
                    n_cells = np.bincount(inverse)
                    n_high_by_d = np.bincount(inverse, weights=high.astype(float)).astype(int)
                    dmeta = donors.set_index("donor_idx").loc[unique_d].reset_index()
                    donor_frame = dmeta.copy()
                    donor_frame["n_cells"] = n_cells
                    donor_frame["n_high"] = n_high_by_d
                    donor_frame["fraction_high"] = n_high_by_d / n_cells
                    donor_frame["score_threshold"] = threshold
                    n_high = int(high.sum())
                donor_rows[arm] = donor_frame
                for model, adjusted in (("primary", True), ("unadjusted", False)):
                    if arm == "rare" and model == "unadjusted":
                        continue
                    test = f"test1_{model}_{arm}" if arm == "full" else "test1_rarefied_1500"
                    if donor_frame.empty:
                        fit = {"status": "UNTESTABLE_DONOR_FLOOR", "n_donors": 0, "n_cells": 0}
                    else:
                        fit = fit_donor_ols(donor_frame, "fraction_high", ds, adjusted)
                    rows.append(
                        {
                            **base,
                            "test": test,
                            "model": model,
                            "depth": "full_depth" if arm == "full" else "rarefied_1500",
                            "score_threshold": threshold,
                            "n_high_cells": n_high,
                            **fit,
                        }
                    )
            if not donor_rows["full"].empty:
                tmp = donor_rows["full"][["donor_idx", "n_cells", "n_high"]].copy()
                tmp.to_csv(CACHE_ROOT / ds / f"donor_high_{ci:02d}.csv", index=False)
    result = apply_test1_bh(pd.DataFrame(rows))
    preferred = [
        "dataset", "cell_type", "test", "model", "depth", "status", "slope_per_decade", "se",
        "ci95_low", "ci95_high", "pvalue", "q_bh", "family_size", "n_donors", "n_cells",
        "n_donors_with_cells", "n_eligible_cells", "n_high_cells", "score_threshold", "outcome_sd",
        "mde80", "mde80_sd_per_decade", "residual_df", "covariates_active", "covariates_omitted",
        "n_source_genes", "n_scoring_genes", "masked_list_genes", "gene_list_source",
    ]
    for col in preferred:
        if col not in result:
            result[col] = np.nan
    result = result[preferred].sort_values(["dataset", "test", "cell_type"]).reset_index(drop=True)
    write_clean_csv(result, DEATH_PATH)
    for ds in ("gao", "yazar"):
        sub = result[(result["dataset"] == ds) & (result["test"] == "test1_primary_full")]
        log(f"{ds}: test 1 evaluated={int((sub['status'] == 'EVALUATED').sum())}, BH passing={int((sub['q_bh'] < 0.05).sum())}")
    return result


GAO_STUDY_COLS = (
    "study_chaffin_2022",
    "study_koenig_2022",
    "study_penn",
    "study_simonson_2023",
    "study_kuppe_2022",
    "study_sim_2021",
)
GENE_MIN_RESID_DF = {"gao": 2, "yazar": 10}
REMOVAL_COLUMNS = [
    "dataset", "cell_type", "row_type", "status", "gene", "n_donors", "n_cells_unremoved",
    "n_cells_removed", "n_high_cells_removed", "score_threshold", "n_genes_unremoved",
    "n_genes_removed", "n_common_genes", "n_pass_unremoved", "n_pass_removed",
    "n_pass_to_fail_predicted", "n_pass_to_fail_not_attenuated", "n_fail_to_pass",
    "n_pass_to_pass", "n_fail_to_fail", "slope_unremoved", "se_unremoved",
    "ci95_low_unremoved", "ci95_high_unremoved", "pvalue_unremoved", "q_bh_unremoved",
    "mde80_unremoved", "slope_removed", "se_removed", "ci95_low_removed",
    "ci95_high_removed", "pvalue_removed", "q_bh_removed", "mde80_removed",
    "slope_change_removed_minus_unremoved", "absolute_slope_attenuated", "pass_fail_change",
    "predicted_pass_to_fail", "covariates_unremoved", "covariates_removed",
    "omitted_covariates_unremoved", "omitted_covariates_removed", "gene_list_source",
]


def empty_removal() -> pd.DataFrame:
    return pd.DataFrame(columns=REMOVAL_COLUMNS)


def load_shard_for_removal(
    ds: str,
    shard: int,
    info: dict[str, Any],
) -> tuple[sp.csr_matrix, np.ndarray, np.ndarray, np.ndarray] | None:
    files = sorted((BLOCK1_CACHE / ds).glob(f"shard{shard:02d}_chunk*.npz"))
    if not files:
        return None
    mats: list[sp.csr_matrix] = []
    meta_parts: list[np.ndarray] = []
    n_genes = len(info["retained_symbols"])
    for file in files:
        z = np.load(file)
        n_rows = len(z["donor"])
        mats.append(sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=(n_rows, n_genes)))
        meta_parts.append(np.column_stack([z["donor"], z["ct"], z["obs_idx"]]))
    mat = sp.vstack(mats, format="csr")
    meta = np.vstack(meta_parts)
    order = np.argsort(meta[:, 2], kind="mergesort")
    mat = mat[order]
    meta = meta[order]
    donor = meta[:, 0].astype(np.int32)
    ct = meta[:, 1].astype(np.int32)
    obs_idx = meta[:, 2].astype(np.int32)
    score_path = CACHE_ROOT / ds / "scores" / f"scores_shard{shard:02d}.npz"
    if not score_path.exists():
        raise FileNotFoundError(f"Missing score shard: {score_path}")
    score = np.load(score_path)
    if not np.array_equal(obs_idx, score["obs_idx"]):
        raise AssertionError(f"Cell order mismatch between count and score shard: {ds} shard {shard}")
    if not np.array_equal(donor, score["donor"]) or not np.array_equal(ct, score["ct"]):
        raise AssertionError(f"Metadata mismatch between count and score shard: {ds} shard {shard}")
    return mat, donor, ct, score["full_score"]


def prepare_removal_specs(
    ds: str,
    hits: pd.DataFrame,
    info: dict[str, Any],
    donors: pd.DataFrame,
) -> list[dict[str, Any]]:
    specs: list[dict[str, Any]] = []
    donor_lookup = donors.set_index("donor_idx")
    for _, hit in hits.iterrows():
        cell_type = str(hit["cell_type"])
        ci = info["cell_types"].index(cell_type)
        donor_high_path = CACHE_ROOT / ds / f"donor_high_{ci:02d}.csv"
        if not donor_high_path.exists():
            raise FileNotFoundError(f"Missing test 1 donor fractions: {donor_high_path}")
        high = pd.read_csv(donor_high_path)
        high["n_removed"] = high["n_cells"] - high["n_high"]
        eligible = high[(high["n_cells"] >= MIN_CELLS) & (high["n_removed"] >= MIN_CELLS)].copy()
        common = np.sort(eligible["donor_idx"].to_numpy(dtype=np.int32))
        if ds == "yazar" and len(common):
            pools = donor_lookup.loc[common, "pool_number"].astype(int)
            pool_counts = pools.value_counts()
            common = np.sort(pools[pools.map(pool_counts) >= 2].index.to_numpy(dtype=np.int32))
        eligible = eligible.set_index("donor_idx").loc[common] if len(common) else eligible.iloc[0:0]
        base = {column: np.nan for column in REMOVAL_COLUMNS}
        base.update(
            dataset=ds,
            cell_type=cell_type,
            row_type="summary",
            status="EVALUATED" if len(common) >= MIN_DONORS else "UNTESTABLE_DONOR_FLOOR",
            gene="ALL_GENES",
            n_donors=int(len(common)),
            n_cells_unremoved=int(eligible["n_cells"].sum()) if len(common) else 0,
            n_cells_removed=int(eligible["n_removed"].sum()) if len(common) else 0,
            n_high_cells_removed=int(eligible["n_high"].sum()) if len(common) else 0,
            score_threshold=float(hit["score_threshold"]),
            n_genes_unremoved=0,
            n_genes_removed=0,
            n_common_genes=0,
            n_pass_unremoved=0,
            n_pass_removed=0,
            n_pass_to_fail_predicted=0,
            n_pass_to_fail_not_attenuated=0,
            n_fail_to_pass=0,
            n_pass_to_pass=0,
            n_fail_to_fail=0,
            gene_list_source=GENE_SOURCE,
        )
        specs.append(
            {
                "dataset": ds,
                "cell_type": cell_type,
                "ci": ci,
                "threshold": float(hit["score_threshold"]),
                "common": common,
                "summary": base,
                "unremoved_blocks": [],
                "removed_blocks": [],
            }
        )
    return specs


def aggregate_removal_blocks(ds: str, specs: list[dict[str, Any]], info: dict[str, Any]) -> None:
    active = [spec for spec in specs if spec["summary"]["status"] == "EVALUATED"]
    if not active:
        return
    n_donors = len(info["donor_ids"])
    n_shards = int(info["n_shards"])
    for shard in range(n_shards):
        loaded = load_shard_for_removal(ds, shard, info)
        if loaded is None:
            continue
        mat, donor, ct, full_score = loaded
        for spec in active:
            local = spec["common"][spec["common"] % n_shards == shard]
            if len(local) == 0:
                continue
            rows = np.flatnonzero((ct == spec["ci"]) & np.isin(donor, local))
            if len(rows) == 0:
                raise AssertionError(f"No cells for {ds} {spec['cell_type']} in expected shard {shard}")
            lookup = np.full(n_donors, -1, dtype=np.int32)
            lookup[local] = np.arange(len(local), dtype=np.int32)
            columns = lookup[donor[rows]]
            if (columns < 0).any():
                raise AssertionError("Donor mapping failure in removal aggregation")
            indicator = sp.csr_matrix(
                (np.ones(len(rows), dtype=np.float64), (np.arange(len(rows)), columns)),
                shape=(len(rows), len(local)),
            )
            ct_mat = mat[rows]
            unremoved = (indicator.T @ ct_mat).T.tocsr()
            high = full_score[rows] >= spec["threshold"]
            removed = (indicator[~high].T @ ct_mat[~high]).T.tocsr()
            spec["unremoved_blocks"].append((local, unremoved))
            spec["removed_blocks"].append((local, removed))
        del mat, donor, ct, full_score
        log(f"{ds}: aggregated removal counts for shard {shard + 1}/{n_shards}")


def combine_count_blocks(
    blocks: list[tuple[np.ndarray, sp.csr_matrix]],
    common: np.ndarray,
    n_genes: int,
) -> np.ndarray:
    if not blocks:
        raise ValueError("No pseudobulk blocks for an evaluable removal cell type")
    donor_order = np.concatenate([block[0] for block in blocks])
    matrix = sp.hstack([block[1] for block in blocks], format="csr")
    order = np.argsort(donor_order, kind="mergesort")
    donor_order = donor_order[order]
    matrix = matrix[:, order]
    if not np.array_equal(donor_order, common):
        raise AssertionError("Aggregated donor columns do not match the common donor set")
    if matrix.shape != (n_genes, len(common)):
        raise AssertionError(f"Unexpected pseudobulk shape {matrix.shape}")
    return matrix.toarray().astype(np.float64, copy=False)


def removal_meta(
    ds: str,
    donors: pd.DataFrame,
    common: np.ndarray,
) -> tuple[pd.DataFrame, tuple[str, ...]]:
    meta = donors.set_index("donor_idx").loc[common].copy()
    meta.index = meta["donor_id"].astype(str)
    if ds == "gao":
        return meta, ("sex_male", "log_lib", *GAO_STUDY_COLS)
    pools = meta["pool_number"].astype(int)
    counts = pools.value_counts()
    ordered = sorted(counts.index.tolist(), key=lambda p: (-int(counts[p]), int(p)))
    reference = ordered[0]
    all_pools = sorted(donors["pool_number"].dropna().astype(int).unique().tolist())
    pool_cols: list[str] = []
    for pool in all_pools:
        column = f"pool_{pool}"
        meta[column] = 0.0 if pool == reference else (pools == pool).astype(float)
        pool_cols.append(column)
    return meta, ("sex_male", "log_lib", *pool_cols)


def fit_removal_arm(
    counts: np.ndarray,
    symbols: list[str],
    donor_ids: list[str],
    meta: pd.DataFrame,
    covariates: tuple[str, ...],
    cell_type: str,
    ds: str,
) -> pd.DataFrame:
    raw = pd.DataFrame(counts, index=symbols, columns=donor_ids)
    log2cpm, library = counts_to_log2cpm(raw)
    arm_meta = meta.copy()
    arm_meta["lib_size"] = library.loc[donor_ids].to_numpy(dtype=float)
    arm_meta["log_lib"] = np.log1p(arm_meta["lib_size"].to_numpy(dtype=float))
    genes, _, _ = run_family_b_genes_and_concordance(
        {cell_type: log2cpm},
        {cell_type: arm_meta},
        age_col="age_per_decade",
        candidate_covariates=covariates,
        min_donors=MIN_DONORS,
        min_cpm=1.0,
        min_expr_frac=0.5,
        min_residual_df=GENE_MIN_RESID_DF[ds],
        is_unbiased_single_cell=True,
        concordance_perms=0,
        seed=MASTER_SEED,
    )
    return genes


def rename_gene_arm(frame: pd.DataFrame, arm: str) -> pd.DataFrame:
    columns = {
        "gene": "gene",
        "effect": f"slope_{arm}",
        "se": f"se_{arm}",
        "ci95_low": f"ci95_low_{arm}",
        "ci95_high": f"ci95_high_{arm}",
        "pvalue": f"pvalue_{arm}",
        "q_bh_within_group": f"q_bh_{arm}",
        "covariates": f"covariates_{arm}",
        "omitted_covariates": f"omitted_covariates_{arm}",
        "residual_df": f"residual_df_{arm}",
    }
    if frame.empty:
        return pd.DataFrame(columns=list(columns.values()))
    out = frame[list(columns)].rename(columns=columns).copy()
    df = frame["residual_df"].to_numpy(dtype=float)
    critical = stats.t.ppf(0.975, df)
    t80 = stats.t.ppf(0.80, df)
    out[f"mde80_{arm}"] = (critical + t80) * out[f"se_{arm}"].to_numpy(dtype=float)
    return out


def compare_removal_arms(
    spec: dict[str, Any],
    unremoved: pd.DataFrame,
    removed: pd.DataFrame,
) -> pd.DataFrame:
    left = rename_gene_arm(unremoved, "unremoved")
    right = rename_gene_arm(removed, "removed")
    genes = left.merge(right, on="gene", how="outer", validate="one_to_one")
    both = genes["slope_unremoved"].notna() & genes["slope_removed"].notna()
    genes["row_type"] = "gene"
    genes["status"] = np.select(
        [both, genes["slope_unremoved"].notna(), genes["slope_removed"].notna()],
        ["EVALUATED_BOTH", "UNREMOVED_ONLY", "REMOVED_ONLY"],
        default="NOT_EXPRESSED_EITHER",
    )
    genes["slope_change_removed_minus_unremoved"] = genes["slope_removed"] - genes["slope_unremoved"]
    genes["absolute_slope_attenuated"] = np.where(
        both,
        genes["slope_removed"].abs() < genes["slope_unremoved"].abs(),
        np.nan,
    )
    unremoved_pass = genes["q_bh_unremoved"] < 0.05
    removed_pass = genes["q_bh_removed"] < 0.05
    pass_to_fail = both & unremoved_pass & ~removed_pass
    fail_to_pass = both & ~unremoved_pass & removed_pass
    predicted = pass_to_fail & (genes["absolute_slope_attenuated"] == True)  # noqa: E712
    genes["pass_fail_change"] = "NOT_COMPARABLE"
    genes.loc[both & ~unremoved_pass & ~removed_pass, "pass_fail_change"] = "FAIL_TO_FAIL"
    genes.loc[both & unremoved_pass & removed_pass, "pass_fail_change"] = "PASS_TO_PASS"
    genes.loc[pass_to_fail & ~predicted, "pass_fail_change"] = "PASS_TO_FAIL_NOT_ATTENUATED"
    genes.loc[predicted, "pass_fail_change"] = "PASS_TO_FAIL_PREDICTED"
    genes.loc[fail_to_pass, "pass_fail_change"] = "FAIL_TO_PASS"
    genes["predicted_pass_to_fail"] = np.where(both, predicted.astype(float), np.nan)
    summary = spec["summary"].copy()
    summary.update(
        status="EVALUATED" if int(both.sum()) > 0 else "UNTESTABLE_NO_COMMON_GENES",
        n_genes_unremoved=int(len(unremoved)),
        n_genes_removed=int(len(removed)),
        n_common_genes=int(both.sum()),
        n_pass_unremoved=int((both & unremoved_pass).sum()),
        n_pass_removed=int((both & removed_pass).sum()),
        n_pass_to_fail_predicted=int(predicted.sum()),
        n_pass_to_fail_not_attenuated=int((pass_to_fail & ~predicted).sum()),
        n_fail_to_pass=int(fail_to_pass.sum()),
        n_pass_to_pass=int((both & unremoved_pass & removed_pass).sum()),
        n_fail_to_fail=int((both & ~unremoved_pass & ~removed_pass).sum()),
    )
    for column in REMOVAL_COLUMNS:
        if column not in genes:
            genes[column] = np.nan
    shared_columns = {
        "dataset", "cell_type", "n_donors", "n_cells_unremoved", "n_cells_removed",
        "n_high_cells_removed", "score_threshold", "n_genes_unremoved", "n_genes_removed",
        "n_common_genes", "n_pass_unremoved", "n_pass_removed", "n_pass_to_fail_predicted",
        "n_pass_to_fail_not_attenuated", "n_fail_to_pass", "n_pass_to_pass", "n_fail_to_fail",
        "gene_list_source",
    }
    for column in shared_columns:
        genes[column] = summary[column]
    summary_row = pd.DataFrame([summary])
    out = pd.concat([summary_row, genes], ignore_index=True, sort=False)
    for column in REMOVAL_COLUMNS:
        if column not in out:
            out[column] = np.nan
    return out[REMOVAL_COLUMNS]


def run_test2(death: pd.DataFrame | None = None) -> pd.DataFrame:
    verify_frozen_spec()
    if death is None:
        if not DEATH_PATH.exists():
            raise FileNotFoundError(f"Run test 1 first: {DEATH_PATH}")
        death = pd.read_csv(DEATH_PATH, float_precision="round_trip")
    pieces: list[pd.DataFrame] = []
    for ds in ("gao", "yazar"):
        info = load_info(ds)
        donors = donor_metadata(ds, info)
        hits = death[
            (death["dataset"] == ds)
            & (death["test"] == "test1_primary_full")
            & (death["status"] == "EVALUATED")
            & (death["q_bh"] < 0.05)
        ].copy().sort_values("cell_type")
        log(f"{ds}: test 2 qualifying cell types={len(hits)}")
        if hits.empty:
            continue
        specs = prepare_removal_specs(ds, hits, info, donors)
        wait_for_memory(f"{ds} removal aggregation")
        aggregate_removal_blocks(ds, specs, info)
        symbols = list(info["retained_symbols"])
        for spec in specs:
            if spec["summary"]["status"] != "EVALUATED":
                pieces.append(pd.DataFrame([spec["summary"]]))
                continue
            common = spec["common"]
            donor_ids = donors.set_index("donor_idx").loc[common, "donor_id"].astype(str).tolist()
            unremoved_counts = combine_count_blocks(spec["unremoved_blocks"], common, len(symbols))
            removed_counts = combine_count_blocks(spec["removed_blocks"], common, len(symbols))
            meta, covariates = removal_meta(ds, donors, common)
            try:
                unremoved = fit_removal_arm(
                    unremoved_counts, symbols, donor_ids, meta, covariates, spec["cell_type"], ds
                )
                removed = fit_removal_arm(
                    removed_counts, symbols, donor_ids, meta, covariates, spec["cell_type"], ds
                )
            except ValueError as exc:
                summary = spec["summary"].copy()
                summary["status"] = f"UNTESTABLE_MODEL_{exc}"
                pieces.append(pd.DataFrame([summary]))
                continue
            compared = compare_removal_arms(spec, unremoved, removed)
            pieces.append(compared)
            summary = compared.iloc[0]
            log(
                f"{ds} {spec['cell_type']}: common genes={int(summary['n_common_genes'])}, "
                f"predicted pass-to-fail={int(summary['n_pass_to_fail_predicted'])}"
            )
            spec["unremoved_blocks"].clear()
            spec["removed_blocks"].clear()
    result = pd.concat(pieces, ignore_index=True, sort=False) if pieces else empty_removal()
    for column in REMOVAL_COLUMNS:
        if column not in result:
            result[column] = np.nan
    result = result[REMOVAL_COLUMNS]
    if not result.empty:
        result["_row_order"] = result["row_type"].map({"summary": 0, "gene": 1}).fillna(2)
        result = result.sort_values(["dataset", "cell_type", "_row_order", "gene"]).drop(columns="_row_order")
    write_clean_csv(result.reset_index(drop=True), REMOVAL_PATH)
    return result.reset_index(drop=True)


def bool_column(frame: pd.DataFrame, column: str) -> pd.Series:
    return frame[column].astype(str).str.lower().isin(["true", "1"])


def run_narrowing() -> pd.DataFrame:
    verify_frozen_spec()
    paper = pd.read_csv(RESULTS_DIR / "reviewer-phase1-composition-global-bh.csv")
    if len(paper) != 431:
        raise ValueError(f"Expected 431 paper composition tests, found {len(paper)}")
    cohorts = paper["cohort"].astype(str).str.lower()
    if cohorts.isin(["gao", "yazar"]).any():
        raise ValueError("Paper composition input unexpectedly contains Gao or Yazar rows")
    rows: list[dict[str, Any]] = []
    for ds in ("gao", "yazar"):
        spread_file = RESULTS_DIR / f"followup-fulldepth-{ds}-spread.csv"
        fraction_file = RESULTS_DIR / (
            "followup-organ-17-heart-gao-fractions.csv"
            if ds == "gao"
            else "followup-hcata-pool-sensitivity-fractions.csv"
        )
        spread = pd.read_csv(spread_file)
        mask = (
            spread["metric"].eq("shannon")
            & spread["depth_mode"].eq("full_depth")
            & spread["analysis_role"].eq("primary_full_depth")
            & spread["model_spec"].eq("primary")
            & ~bool_column(spread, "is_pooled")
            & ~bool_column(spread, "is_stratified_global")
        )
        spread = spread.loc[mask].copy()
        if spread["group"].duplicated().any():
            raise ValueError(f"Duplicate spread cell-type keys in {spread_file}")
        fraction = pd.read_csv(fraction_file)
        if ds == "gao":
            fraction = fraction[fraction["subset"].eq("pooled_7study_primary")].copy()
        else:
            fraction = fraction[fraction["study"].eq("yazar") & fraction["arm"].eq("pool_fe")].copy()
        if fraction["cell_type"].duplicated().any():
            raise ValueError(f"Duplicate fraction cell-type keys in {fraction_file}")
        spread = spread.rename(columns={"group": "cell_type"})
        joined = spread.merge(fraction, on="cell_type", how="inner", suffixes=("_spread", "_fraction"))
        for _, item in joined.sort_values("cell_type").iterrows():
            fraction_slope = float(item["effect_clr"])
            rows.append(
                {
                    "dataset": ds,
                    "row_type": "cell_type",
                    "cell_type": item["cell_type"],
                    "status": "MATCHED",
                    "spread_slope_per_decade": float(item["slope_raw_per_decade"]),
                    "spread_se": float(item["se_raw"]),
                    "spread_pvalue": float(item["pvalue_spread"]),
                    "spread_q_bh": float(item["q_bh_within_slice"]),
                    "spread_n_donors": int(item["n_donors_spread"]),
                    "fraction_slope_clr_per_decade": fraction_slope,
                    "fraction_se": float(item["se"]),
                    "fraction_pvalue": float(item["pvalue_fraction"]),
                    "fraction_q_bh": float(item["q_bh_family_a"]),
                    "fraction_n_donors": int(item["n_donors_fraction"]),
                    "fraction_direction": "losing" if fraction_slope < 0 else "gaining",
                    "spread_source": str(spread_file.relative_to(REPO_ROOT)),
                    "fraction_source": str(fraction_file.relative_to(REPO_ROOT)),
                }
            )
        model = joined[np.isfinite(joined["slope_raw_per_decade"]) & np.isfinite(joined["effect_clr"])].copy()
        model = model[np.isfinite(model["se_raw"]) & (model["se_raw"] > 0)]
        base: dict[str, Any] = {
            "dataset": ds,
            "row_type": "cross_type_model",
            "cell_type": "ALL_MATCHED_CELL_TYPES",
            "status": "EVALUATED",
            "n_cell_types": int(len(model)),
            "n_fraction_losing": int((model["effect_clr"] < 0).sum()),
            "n_fraction_gaining": int((model["effect_clr"] >= 0).sum()),
            "mean_spread_slope_losing": float(model.loc[model["effect_clr"] < 0, "slope_raw_per_decade"].mean()),
            "mean_spread_slope_gaining": float(model.loc[model["effect_clr"] >= 0, "slope_raw_per_decade"].mean()),
            "spread_source": str(spread_file.relative_to(REPO_ROOT)),
            "fraction_source": str(fraction_file.relative_to(REPO_ROOT)),
            "weighting": "inverse_variance_1_over_spread_se_squared",
            "family_size": 1,
            "q_bh": np.nan,
        }
        if len(model) < 3 or model["effect_clr"].nunique() < 2:
            base["status"] = "UNTESTABLE_CELL_TYPE_FLOOR"
        else:
            y = model["slope_raw_per_decade"].to_numpy(dtype=float)
            x = model["effect_clr"].to_numpy(dtype=float)
            X = sm.add_constant(pd.DataFrame({"fraction_slope": x}), has_constant="add")
            weights = 1.0 / np.square(model["se_raw"].to_numpy(dtype=float))
            fit = sm.WLS(y, X, weights=weights).fit()
            slope = float(fit.params["fraction_slope"])
            se = float(fit.bse["fraction_slope"])
            df = int(fit.df_resid)
            t975 = float(stats.t.ppf(0.975, df))
            t80 = float(stats.t.ppf(0.80, df))
            mde = (t975 + t80) * se
            base.update(
                association_slope=slope,
                association_se=se,
                ci95_low=slope - t975 * se,
                ci95_high=slope + t975 * se,
                pvalue=float(fit.pvalues["fraction_slope"]),
                residual_df=df,
                mde80=mde,
                spread_slope_sd=float(np.std(y, ddof=1)),
                mde80_spread_sd_units=mde / float(np.std(y, ddof=1)),
                weight_sum=float(weights.sum()),
            )
        rows.append(base)
        log(f"{ds}: narrowing join={len(joined)}, weighted model cell types={len(model)}")
    result = pd.DataFrame(rows)
    preferred = [
        "dataset", "row_type", "cell_type", "status", "association_slope", "association_se",
        "ci95_low", "ci95_high", "pvalue", "q_bh", "family_size", "mde80",
        "mde80_spread_sd_units", "residual_df", "n_cell_types", "n_fraction_losing",
        "n_fraction_gaining", "mean_spread_slope_losing", "mean_spread_slope_gaining",
        "spread_slope_sd", "weight_sum", "weighting", "spread_slope_per_decade", "spread_se",
        "spread_pvalue", "spread_q_bh", "spread_n_donors", "fraction_slope_clr_per_decade",
        "fraction_se", "fraction_pvalue", "fraction_q_bh", "fraction_n_donors",
        "fraction_direction", "spread_source", "fraction_source",
    ]
    for column in preferred:
        if column not in result:
            result[column] = np.nan
    result = result[preferred]
    result["_order"] = result["row_type"].map({"cell_type": 0, "cross_type_model": 1})
    result = result.sort_values(["dataset", "_order", "cell_type"]).drop(columns="_order").reset_index(drop=True)
    write_clean_csv(result, NARROWING_PATH)
    return result


def available_memory_gib() -> float:
    for line in Path("/proc/meminfo").read_text(encoding="utf-8").splitlines():
        if line.startswith("MemAvailable:"):
            return float(line.split()[1]) / (1024.0 ** 2)
    raise RuntimeError("MemAvailable is missing from /proc/meminfo")


def wait_for_memory(label: str) -> None:
    while True:
        available = available_memory_gib()
        if available >= 20.0:
            log(f"Memory guard for {label}: {available:.1f} Gi available")
            return
        log(f"Memory guard for {label}: only {available:.1f} Gi available; waiting 60 seconds")
        time.sleep(60)


def part_a_verdicts(death: pd.DataFrame, removal: pd.DataFrame) -> dict[str, str]:
    verdicts: dict[str, str] = {}
    for ds in ("gao", "yazar"):
        primary = death[(death["dataset"] == ds) & (death["test"] == "test1_primary_full")]
        evaluated = primary[primary["status"] == "EVALUATED"]
        hits = evaluated[evaluated["q_bh"] < 0.05]
        if len(hits):
            predicted = 0
            if not removal.empty:
                summaries = removal[(removal["dataset"] == ds) & (removal["row_type"] == "summary")]
                predicted = int(summaries["n_pass_to_fail_predicted"].fillna(0).sum())
            verdicts[ds] = "PROVES" if predicted >= 1 else "INCONCLUSIVE"
        elif len(evaluated) and (evaluated["mde80_sd_per_decade"] <= MDE80_SMALL_SD).all():
            verdicts[ds] = "DOES NOT PROVE"
        else:
            verdicts[ds] = "INCONCLUSIVE"
    if "PROVES" in verdicts.values():
        verdicts["overall"] = "PROVES"
    elif all(value == "DOES NOT PROVE" for value in verdicts.values()):
        verdicts["overall"] = "DOES NOT PROVE"
    else:
        verdicts["overall"] = "INCONCLUSIVE"
    return verdicts


def part_b_verdict(narrowing: pd.DataFrame) -> str:
    models = narrowing[narrowing["row_type"] == "cross_type_model"]
    if len(models) != 2 or not (models["status"] == "EVALUATED").all():
        return "INCONCLUSIVE"
    positive = (models["association_slope"] > 0) & (models["pvalue"] < 0.05)
    opposite = (models["association_slope"] < 0) & (models["pvalue"] < 0.05)
    if positive.any() and not opposite.any():
        return "PROVES"
    if positive.any() and opposite.any():
        return "INCONCLUSIVE"
    null_or_opposite = (models["pvalue"] >= 0.05) | (models["association_slope"] < 0)
    return "DOES NOT PROVE" if null_or_opposite.all() else "INCONCLUSIVE"


def fmt(value: Any, digits: int = 4, signed: bool = False) -> str:
    if value is None or pd.isna(value):
        return "NA"
    return f"{float(value):+.{digits}f}" if signed else f"{float(value):.{digits}f}"


def fmt_p(value: Any) -> str:
    if value is None or pd.isna(value):
        return "NA"
    return f"{float(value):.3e}"


def markdown_table(headers: list[str], rows: list[list[Any]]) -> str:
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    lines.extend("| " + " | ".join(str(value) for value in row) + " |" for row in rows)
    return "\n".join(lines)


def generate_report() -> dict[str, str]:
    frozen = verify_frozen_spec()
    death = pd.read_csv(DEATH_PATH)
    removal = pd.read_csv(REMOVAL_PATH)
    narrowing = pd.read_csv(NARROWING_PATH)
    paper = pd.read_csv(RESULTS_DIR / "reviewer-phase1-composition-global-bh.csv")
    if len(paper) != 431:
        raise ValueError(f"Expected 431 paper composition tests, found {len(paper)}")
    verdict_a = part_a_verdicts(death, removal)
    verdict_b = part_b_verdict(narrowing)
    lines: list[str] = [
        "",
        "## Results",
        "",
        "All three result CSVs were generated from the frozen protocol above. Primary analyses use full-depth counts. Donor is the unit in every regression.",
        "",
        "### Part A test 1: high death/stress-score cell fractions",
        "",
    ]
    for ds in ("gao", "yazar"):
        primary = death[(death["dataset"] == ds) & (death["test"] == "test1_primary_full")].sort_values("cell_type")
        evaluated = primary[primary["status"] == "EVALUATED"]
        hits = evaluated[evaluated["q_bh"] < 0.05]
        lines.extend(
            [
                f"#### {ds.capitalize()}",
                "",
                f"Evaluated {len(evaluated)} of {len(primary)} cell types. BH family size {int(evaluated['family_size'].max()) if len(evaluated) else 0}; {len(hits)} passed q < 0.05.",
                "",
                markdown_table(
                    ["Cell type", "Status", "N donors", "Slope per decade", "95 percent CI", "p", "BH q", "MDE80 SD"],
                    [
                        [
                            row["cell_type"], row["status"], int(row["n_donors"]),
                            fmt(row["slope_per_decade"], 6, True),
                            f"{fmt(row['ci95_low'], 6)} to {fmt(row['ci95_high'], 6)}",
                            fmt_p(row["pvalue"]), fmt_p(row["q_bh"]), fmt(row["mde80_sd_per_decade"], 3),
                        ]
                        for _, row in primary.iterrows()
                    ],
                ),
                "",
            ]
        )
        sensitivity = death[(death["dataset"] == ds) & (death["test"] != "test1_primary_full")]
        for test, group in sensitivity.groupby("test"):
            evaluated_s = group[group["status"] == "EVALUATED"]
            lines.append(
                f"{ds} {test}: {len(evaluated_s)} evaluated, {int((evaluated_s['q_bh'] < 0.05).sum())} BH hits."
            )
        lines.append("")
    lines.extend(["### Part A test 2: removal sensitivity", ""])
    summaries = removal[removal["row_type"] == "summary"] if not removal.empty else removal
    if summaries.empty:
        lines.append("No cell type passed primary test 1 BH correction, so test 2 was not run. The removal CSV contains the complete required schema and no gene rows.")
    else:
        lines.append(
            markdown_table(
                ["Dataset", "Cell type", "Status", "N donors", "Common genes", "Unremoved pass", "Removed pass", "Predicted pass to fail", "Other pass to fail", "Fail to pass"],
                [
                    [
                        row["dataset"], row["cell_type"], row["status"], int(row["n_donors"]),
                        int(row["n_common_genes"]), int(row["n_pass_unremoved"]), int(row["n_pass_removed"]),
                        int(row["n_pass_to_fail_predicted"]), int(row["n_pass_to_fail_not_attenuated"]),
                        int(row["n_fail_to_pass"]),
                    ]
                    for _, row in summaries.iterrows()
                ],
            )
        )
    lines.extend(
        [
            "",
            "Pass/fail change counts use genes evaluated in both arms. Slope change is removed minus unremoved. Removal is a lower bound because dying cells are under-captured.",
            "",
            "### Part B: narrowing localization",
            "",
        ]
    )
    models = narrowing[narrowing["row_type"] == "cross_type_model"]
    lines.append(
        markdown_table(
            ["Dataset", "N types", "Association slope", "95 percent CI", "p", "MDE80", "Losing N", "Gaining N", "Mean spread losing", "Mean spread gaining"],
            [
                [
                    row["dataset"], int(row["n_cell_types"]), fmt(row["association_slope"], 5, True),
                    f"{fmt(row['ci95_low'], 5)} to {fmt(row['ci95_high'], 5)}", fmt_p(row["pvalue"]),
                    fmt(row["mde80"], 5), int(row["n_fraction_losing"]), int(row["n_fraction_gaining"]),
                    fmt(row["mean_spread_slope_losing"], 5, True), fmt(row["mean_spread_slope_gaining"], 5, True),
                ]
                for _, row in models.iterrows()
            ],
        )
    )
    lines.extend(
        [
            "",
            "The association slope is the change in Shannon spread slope per unit CLR fraction slope per decade. Fraction slope negative means losing share. The prediction is a positive association. There is one test per dataset and no multiplicity correction across datasets.",
            "",
            "## Verdicts",
            "",
            f"- Part A Gao: {verdict_a['gao']}.",
            f"- Part A Yazar: {verdict_a['yazar']}.",
            f"- Part A overall: {verdict_a['overall']}.",
            f"- Part B: {verdict_b}.",
            "",
            "For Part A, PROVES requires both a BH-passing test 1 result and at least one predicted pass-to-fail gene change with absolute-slope attenuation. A no-hit DOES NOT PROVE verdict requires every evaluated primary null to have MDE80 at most 0.50 outcome SD per decade. Larger MDE80 or a required but unevaluable removal comparison is INCONCLUSIVE.",
            "",
            "Neither Part A nor Part B alone proves survivor bias as a mechanism. These tests address within-type stress-cell dependence and cross-type localization only. The Nov 17 cross-type correlation was not rerun.",
            "",
            "## Reproducibility and input audit",
            "",
            f"The paper composition file contains {len(paper)} tests as required. It is not part of either cross-type model because it has no Gao or Yazar cohort rows. Gao joins 13 cell types and Yazar joins 17 cell types. Exact spread and fraction filenames, model rows, mask, gene source, fixed seeds, families, and eligibility rules are recorded in the frozen protocol.",
            "",
        ]
    )
    write_clean_markdown(frozen + "\n".join(lines), REPORT_PATH)
    log(f"Report complete: Part A {verdict_a['overall']}; Part B {verdict_b}")
    return {"part_a_gao": verdict_a["gao"], "part_a_yazar": verdict_a["yazar"], "part_a_overall": verdict_a["overall"], "part_b": verdict_b}


def validate_outputs() -> None:
    for path in (DEATH_PATH, REMOVAL_PATH, NARROWING_PATH, REPORT_PATH):
        if not path.exists():
            raise FileNotFoundError(path)
        text = path.read_text(encoding="utf-8")
        if re.search(r"[\u2013\u2014]", text):
            raise ValueError(f"Forbidden dash in {path}")
    death = pd.read_csv(DEATH_PATH)
    evaluated = death[death["status"] == "EVALUATED"]
    if not np.isfinite(evaluated[["slope_per_decade", "se", "pvalue", "q_bh", "mde80"]].to_numpy(dtype=float)).all():
        raise ValueError("Non-finite evaluated test 1 result")
    if not evaluated["pvalue"].between(0, 1).all() or not evaluated["q_bh"].between(0, 1).all():
        raise ValueError("Invalid test 1 p or q value")
    if (evaluated["n_donors"] < MIN_DONORS).any():
        raise ValueError("Evaluated test 1 row below donor floor")
    narrowing = pd.read_csv(NARROWING_PATH)
    models = narrowing[narrowing["row_type"] == "cross_type_model"]
    if len(models) != 2 or len(narrowing[narrowing["row_type"] == "cell_type"]) != 30:
        raise ValueError("Unexpected narrowing output row counts")
    removal = pd.read_csv(REMOVAL_PATH)
    if not removal.empty:
        common = removal[removal["status"] == "EVALUATED_BOTH"]
        if len(common) and not np.isfinite(common[["slope_unremoved", "slope_removed", "slope_change_removed_minus_unremoved"]].to_numpy(dtype=float)).all():
            raise ValueError("Non-finite paired removal result")
    if "Final scoring symbols:" not in REPORT_PATH.read_text(encoding="utf-8") or "## Verdicts" not in REPORT_PATH.read_text(encoding="utf-8"):
        raise ValueError("Report is missing the frozen list or verdicts")
    log("Output validation passed")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze-spec", action="store_true", help="Freeze the report protocol and exit")
    parser.add_argument(
        "--stage",
        choices=["all", "freeze", "scores", "test1", "test2", "narrowing", "report", "validate"],
        default="all",
    )
    parser.add_argument("--dataset", choices=["both", "gao", "yazar"], default="both")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.freeze_spec or args.stage == "freeze":
        freeze()
        return
    datasets = ("gao", "yazar") if args.dataset == "both" else (args.dataset,)
    if args.stage == "all":
        if not (CACHE_ROOT / "frozen_spec.md").exists():
            freeze()
        else:
            verify_frozen_spec()
        wait_for_memory("Block 7 score cache")
        for ds in datasets:
            run_scores(ds)
        death = run_test1()
        run_test2(death)
        run_narrowing()
        generate_report()
        validate_outputs()
        return
    verify_frozen_spec()
    if args.stage == "scores":
        wait_for_memory("Block 7 score cache")
        for ds in datasets:
            run_scores(ds)
    elif args.stage == "test1":
        run_test1()
    elif args.stage == "test2":
        run_test2()
    elif args.stage == "narrowing":
        run_narrowing()
    elif args.stage == "report":
        generate_report()
    elif args.stage == "validate":
        validate_outputs()


if __name__ == "__main__":
    main()

