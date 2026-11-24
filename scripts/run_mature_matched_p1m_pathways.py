#!/usr/bin/env python3
"""GO Biological Process direction agreement across matched mouse-human P1m pairs."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_var] = "1"

import numpy as np
import pandas as pd
import scipy
import scipy.sparse as sp

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts import run_mature_matched_p1m as matched_test

RESULTS = REPO / "results"
REPORTS = REPO / "reports"
WORK_DIR = matched_test.WORK_DIR
PAIR_GENE_CSV = matched_test.PAIR_GENE_CSV
GO_CSV = WORK_DIR / "mouse-go-bp.csv"
GO_METADATA_CSV = WORK_DIR / "mouse-go-bp-metadata.csv"
OUTPUT_CSV = RESULTS / "followup-mature-matched-p1m-pathways.csv"
SUMMARY_CSV = RESULTS / "followup-mature-matched-p1m-pathways-summary.csv"
OUTPUT_REPORT = REPORTS / "followup-mature-matched-p1m-pathways.md"
PAIR_KEYS = ["mouse_unit", "human_cohort", "human_unit"]
MIN_GENES = 5
SEED = 20261126


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def load_matches() -> pd.DataFrame:
    frame = pd.read_csv(RESULTS / "followup-mature-matched-p1m.csv")
    frame = frame[frame["arm"].eq("unfloored")].copy()
    assert len(frame) == 28
    assert not frame.duplicated(PAIR_KEYS).any()
    assert frame["n_human_genes_tested"].sum() == 299
    return frame.sort_values(PAIR_KEYS).reset_index(drop=True)


def reconstruct_pair_genes(
    matches: pd.DataFrame,
    mouse_units: dict[str, pd.DataFrame],
    human_units: dict[str, dict],
) -> pd.DataFrame:
    expected_keys = {f"{row.human_cohort}::{row.human_unit}" for row in matches.itertuples()}
    assert set(human_units) == expected_keys
    rows = []
    for match in matches.itertuples():
        mouse, _, _ = matched_test.resolve_mouse_orthologs(mouse_units[match.mouse_unit])
        key = f"{match.human_cohort}::{match.human_unit}"
        unit = matched_test.prepare_unit(human_units[key])
        slopes = matched_test.fwl_slopes(unit, unit["metadata"]["age"].to_numpy(float))
        matched_test.verify_fwl(unit, slopes)
        gene_index = {gene: index for index, gene in enumerate(unit["genes"])}
        present = mouse[mouse["human_ortholog"].isin(gene_index)].copy()
        human_indices = present["human_ortholog"].map(gene_index).to_numpy(dtype=int)
        observed = slopes[human_indices] if len(human_indices) else np.array([])
        finite = np.isfinite(observed)
        present = present.iloc[np.flatnonzero(finite)].copy()
        human_indices = human_indices[finite]
        observed = observed[finite]
        for row, human_slope, human_index in zip(present.itertuples(), observed, human_indices):
            rows.append({
                "mouse_unit": match.mouse_unit,
                "mouse_tissue": match.mouse_tissue,
                "human_cohort": match.human_cohort,
                "human_unit": match.human_unit,
                "mouse_gene_symbol": row.mouse_gene_symbol,
                "human_ortholog": row.human_ortholog,
                "mouse_log2FC": row.log2FC,
                "mouse_FDR": row.unit_FDR,
                "human_slope": human_slope,
                "human_slope_gene_index": int(human_index),
                "same_sign": bool(np.sign(row.log2FC) == np.sign(human_slope)),
            })
    columns = [
        "mouse_unit", "mouse_tissue", "human_cohort", "human_unit",
        "mouse_gene_symbol", "human_ortholog", "mouse_log2FC", "mouse_FDR",
        "human_slope", "human_slope_gene_index", "same_sign",
    ]
    return pd.DataFrame(rows, columns=columns).sort_values(
        PAIR_KEYS + ["mouse_gene_symbol"]
    ).reset_index(drop=True)


def load_pair_genes(
    matches: pd.DataFrame,
    mouse_units: dict[str, pd.DataFrame],
    genes: pd.DataFrame,
    refresh: bool,
) -> pd.DataFrame:
    target_genes = sorted(genes.loc[genes["human_ortholog"].notna(), "human_ortholog"].astype(str).unique())
    print("Loading owner-approved human-units.pkl cache (Reza, 2026-11-24 09:10 PT).", flush=True)
    human_units = matched_test.load_human_units(matches, target_genes, refresh=refresh)
    reconstructed = reconstruct_pair_genes(matches, mouse_units, human_units)
    if PAIR_GENE_CSV.exists() and not refresh:
        cached = pd.read_csv(PAIR_GENE_CSV).sort_values(
            PAIR_KEYS + ["mouse_gene_symbol"]
        ).reset_index(drop=True)
        pd.testing.assert_frame_equal(
            reconstructed, cached, check_dtype=False, check_exact=False, rtol=1e-12, atol=1e-12
        )
        pair_genes = cached
        print(f"Validated and reused {len(cached)} cached pair-gene rows.", flush=True)
    else:
        pair_genes = reconstructed
        WORK_DIR.mkdir(parents=True, exist_ok=True)
        pair_genes.to_csv(PAIR_GENE_CSV, index=False)
        print(f"Rebuilt {len(pair_genes)} pair-gene rows.", flush=True)
    expected_nonzero = set(map(tuple, matches.loc[matches["n_human_genes_tested"].gt(0), PAIR_KEYS].to_numpy()))
    actual_nonzero = set(map(tuple, pair_genes[PAIR_KEYS].drop_duplicates().to_numpy()))
    assert actual_nonzero == expected_nonzero
    assert len(pair_genes) == 299
    counts = pair_genes.groupby(PAIR_KEYS).size()
    expected_counts = matches.set_index(PAIR_KEYS).loc[counts.index, "n_human_genes_tested"]
    pd.testing.assert_series_equal(counts, expected_counts, check_names=False, check_dtype=False)
    return pair_genes


def build_go_mapping(genes: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    gene_path = WORK_DIR / "pathway-mouse-symbols.txt"
    gene_path.write_text("\n".join(genes["mouse_gene_symbol"].astype(str)) + "\n", encoding="utf-8")
    r_code = r'''
args <- commandArgs(trailingOnly=TRUE)
.libPaths(c("/home/reza/R/library", .libPaths()))
symbols <- readLines(args[1])
suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(GO.db))
suppressPackageStartupMessages(library(AnnotationDbi))
mapped <- AnnotationDbi::select(org.Mm.eg.db, keys=symbols, keytype="SYMBOL",
                                columns=c("GOALL", "ONTOLOGYALL"))
mapped <- mapped[!is.na(mapped$GOALL) & mapped$ONTOLOGYALL == "BP", c("SYMBOL", "GOALL")]
mapped <- unique(mapped)
terms <- unname(AnnotationDbi::Term(GO.db::GOTERM[mapped$GOALL]))
out <- data.frame(mouse_gene_symbol=mapped$SYMBOL, function_id=mapped$GOALL,
                  function_name=terms, stringsAsFactors=FALSE)
out <- out[order(out$mouse_gene_symbol, out$function_id), ]
write.csv(out, args[2], row.names=FALSE, na="")
mm <- AnnotationDbi::metadata(org.Mm.eg.db)
gg <- AnnotationDbi::metadata(GO.db)
meta <- rbind(
  data.frame(key=c("R_version", "org.Mm.eg.db_package_version", "GO.db_package_version",
                   "AnnotationDbi_package_version", "BiocManager_package_version"),
             value=c(R.version.string, as.character(packageVersion("org.Mm.eg.db")),
                     as.character(packageVersion("GO.db")), as.character(packageVersion("AnnotationDbi")),
                     as.character(packageVersion("BiocManager")))),
  data.frame(key=paste0("org.Mm.eg.db:", mm$name), value=mm$value),
  data.frame(key=paste0("GO.db:", gg$name), value=gg$value)
)
write.csv(meta, args[3], row.names=FALSE, na="")
'''
    subprocess.run(
        ["Rscript", "--vanilla", "-e", r_code, str(gene_path), str(GO_CSV), str(GO_METADATA_CSV)],
        check=True, text=True,
    )
    annotation = pd.read_csv(GO_CSV)
    metadata = pd.read_csv(GO_METADATA_CSV)
    assert annotation["function_id"].str.startswith("GO:").all()
    assert annotation["function_name"].notna().all()
    assert not annotation.duplicated(["mouse_gene_symbol", "function_id"]).any()
    annotation = annotation[annotation["mouse_gene_symbol"].isin(genes["mouse_gene_symbol"])].copy()
    annotation.to_csv(GO_CSV, index=False)
    return annotation, metadata


def dense(matrix) -> np.ndarray:
    return np.asarray(matrix.todense()) if sp.issparse(matrix) else np.asarray(matrix)


def agreement_counts(
    mouse_mean: np.ndarray,
    human_mean: np.ndarray,
    mouse_floor: float | np.ndarray,
    human_floor: float | np.ndarray,
) -> dict[str, float]:
    mouse_near = np.abs(mouse_mean) < mouse_floor
    human_near = np.abs(human_mean) < human_floor
    unclear = mouse_near | human_near
    same = np.sign(mouse_mean) == np.sign(human_mean)
    agree = int(np.sum(~unclear & same))
    disagree = int(np.sum(~unclear & ~same))
    n_unclear = int(np.sum(unclear))
    testable = len(mouse_mean)
    decisive = agree + disagree
    return {
        "n_testable_rows": testable,
        "n_agree": agree,
        "n_disagree": disagree,
        "n_unclear": n_unclear,
        "agree_per_testable": agree / testable if testable else np.nan,
        "directional_agreement": agree / decisive if decisive else np.nan,
    }


def noise_floor_lookup(samples: dict[int, list[np.ndarray]], max_count: int) -> np.ndarray:
    quantiles = {count: float(np.quantile(np.concatenate(values), 0.95))
                 for count, values in samples.items()}
    available = np.asarray(sorted(quantiles))
    lookup = np.full(max_count + 1, np.nan)
    for count in range(MIN_GENES, max_count + 1):
        nearest = available[np.argmin(np.abs(available - count))]
        lookup[count] = quantiles[int(nearest)]
    return lookup


def permutation_analysis(
    membership: sp.csr_matrix,
    mouse_values: np.ndarray,
    human_values: np.ndarray,
    permutations: int,
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    mouse_mask = np.isfinite(mouse_values).astype(float)
    human_mask = np.isfinite(human_values).astype(float)
    features = np.column_stack([
        mouse_mask, np.nan_to_num(mouse_values), human_mask, np.nan_to_num(human_values)
    ])
    rng = np.random.default_rng(SEED)
    null_rows = []
    mouse_samples: dict[int, list[np.ndarray]] = {}
    human_samples: dict[int, list[np.ndarray]] = {}
    transpose = membership.T.tocsr()
    for iteration in range(permutations):
        order = rng.permutation(membership.shape[0])
        aggregate = dense(transpose @ features[order])
        cm, sm, ch, sh = np.split(aggregate, 4, axis=1)
        mouse_mean = np.divide(sm, cm, out=np.full_like(sm, np.nan), where=cm > 0)
        human_mean = np.divide(sh, ch, out=np.full_like(sh, np.nan), where=ch > 0)
        mouse_ok = cm >= MIN_GENES
        human_ok = ch >= MIN_GENES
        for count in np.unique(cm[mouse_ok]).astype(int):
            selected = mouse_ok & (cm == count)
            mouse_samples.setdefault(count, []).append(np.abs(mouse_mean[selected]))
        for count in np.unique(ch[human_ok]).astype(int):
            selected = human_ok & (ch == count)
            human_samples.setdefault(count, []).append(np.abs(human_mean[selected]))
        both = mouse_ok & human_ok
        function_index, pair_index = np.nonzero(both)
        null_rows.append((
            function_index.astype(np.int32), pair_index.astype(np.int32),
            cm[both].astype(np.int32), ch[both].astype(np.int32),
            mouse_mean[both], human_mean[both],
        ))
        if (iteration + 1) % 100 == 0:
            print(f"GO-label permutations {iteration + 1}/{permutations}", flush=True)
    max_count = max(np.max(cm), np.max(ch), mouse_values.shape[0])
    mouse_floor = noise_floor_lookup(mouse_samples, max_count)
    human_floor = noise_floor_lookup(human_samples, max_count)
    null = {
        "all": {key: [] for key in ("agree_per_testable", "directional_agreement")},
        "repeatable": {key: [] for key in ("agree_per_testable", "directional_agreement")},
    }
    for function_index, pair_index, cm_both, ch_both, mouse_mean, human_mean in null_rows:
        mouse_thresholds = mouse_floor[cm_both]
        human_thresholds = human_floor[ch_both]
        counts = agreement_counts(mouse_mean, human_mean, mouse_thresholds, human_thresholds)
        for key in null["all"]:
            null["all"][key].append(counts[key])
        n_pairs_per_function = np.bincount(function_index, minlength=membership.shape[1])
        repeatable = n_pairs_per_function >= 3
        keep = repeatable[function_index]
        repeat_counts = agreement_counts(
            mouse_mean[keep], human_mean[keep], mouse_thresholds[keep], human_thresholds[keep]
        )
        for key in null["repeatable"]:
            null["repeatable"][key].append(repeat_counts[key])
    assert len(null_rows) == permutations
    return mouse_floor, human_floor, {
        label: {key: np.asarray(values, dtype=float) for key, values in metrics.items()}
        for label, metrics in null.items()
    }


def build_observed(
    matches: pd.DataFrame,
    mouse_units: dict[str, pd.DataFrame],
    pair_genes: pd.DataFrame,
    genes: pd.DataFrame,
    annotation: pd.DataFrame,
) -> tuple[pd.DataFrame, sp.csr_matrix, np.ndarray, np.ndarray]:
    gene_names = genes["mouse_gene_symbol"].astype(str).tolist()
    gene_index = {gene: index for index, gene in enumerate(gene_names)}
    functions = annotation[["function_id", "function_name"]].drop_duplicates().sort_values(
        "function_id"
    ).reset_index(drop=True)
    function_index = {term: index for index, term in enumerate(functions["function_id"])}
    rows = annotation["mouse_gene_symbol"].map(gene_index).to_numpy(dtype=int)
    columns = annotation["function_id"].map(function_index).to_numpy(dtype=int)
    membership = sp.csr_matrix(
        (np.ones(len(annotation)), (rows, columns)), shape=(len(gene_names), len(functions))
    )
    n_pairs = len(matches)
    mouse_values = np.full((len(gene_names), n_pairs), np.nan)
    human_values = np.full((len(gene_names), n_pairs), np.nan)
    for pair_index, match in enumerate(matches.itertuples()):
        unit = mouse_units[match.mouse_unit]
        indices = unit["mouse_gene_symbol"].map(gene_index).to_numpy(dtype=int)
        mouse_values[indices, pair_index] = unit["log2FC"].to_numpy(float)
    for pair_index, match in enumerate(matches.itertuples()):
        key = (match.mouse_unit, match.human_cohort, match.human_unit)
        subset = pair_genes[
            pair_genes["mouse_unit"].eq(key[0])
            & pair_genes["human_cohort"].eq(key[1])
            & pair_genes["human_unit"].eq(key[2])
        ]
        indices = subset["mouse_gene_symbol"].map(gene_index).to_numpy(dtype=int)
        assert len(np.unique(indices)) == len(indices)
        human_values[indices, pair_index] = subset["human_slope"].to_numpy(float)
    return functions, membership, mouse_values, human_values


def aggregate_observed(
    matches: pd.DataFrame,
    functions: pd.DataFrame,
    membership: sp.csr_matrix,
    mouse_values: np.ndarray,
    human_values: np.ndarray,
    mouse_floor: float | np.ndarray,
    human_floor: float | np.ndarray,
) -> pd.DataFrame:
    mouse_mask = np.isfinite(mouse_values).astype(float)
    human_mask = np.isfinite(human_values).astype(float)
    transpose = membership.T.tocsr()
    cm = dense(transpose @ mouse_mask)
    ch = dense(transpose @ human_mask)
    sm = dense(transpose @ np.nan_to_num(mouse_values))
    sh = dense(transpose @ np.nan_to_num(human_values))
    mouse_up = dense(transpose @ (mouse_values > 0).astype(float))
    mouse_down = dense(transpose @ (mouse_values < 0).astype(float))
    human_up = dense(transpose @ (human_values > 0).astype(float))
    human_down = dense(transpose @ (human_values < 0).astype(float))
    n_functions, n_pairs = cm.shape
    function_ids = functions["function_id"].to_numpy()
    pair_frame = matches.reset_index(drop=True)
    fi = np.repeat(np.arange(n_functions), n_pairs)
    pi = np.tile(np.arange(n_pairs), n_functions)
    cm_flat = cm.ravel()
    ch_flat = ch.ravel()
    mouse_mean = np.divide(sm.ravel(), cm_flat, out=np.full(cm.size, np.nan), where=cm_flat > 0)
    human_mean = np.divide(sh.ravel(), ch_flat, out=np.full(ch.size, np.nan), where=ch_flat > 0)
    testable = (cm_flat >= MIN_GENES) & (ch_flat >= MIN_GENES)
    mouse_thresholds = np.full(len(fi), np.nan)
    human_thresholds = np.full(len(fi), np.nan)
    mouse_thresholds[testable] = mouse_floor[cm_flat[testable].astype(int)]
    human_thresholds[testable] = human_floor[ch_flat[testable].astype(int)]
    verdict = np.full(len(fi), "UNTESTABLE", dtype=object)
    direction = np.full(len(fi), "", dtype=object)
    unclear = testable & ((np.abs(mouse_mean) < mouse_thresholds) | (np.abs(human_mean) < human_thresholds))
    verdict[unclear] = "UNCLEAR"
    decisive = testable & ~unclear
    agree = decisive & (np.sign(mouse_mean) == np.sign(human_mean))
    disagree = decisive & (np.sign(mouse_mean) != np.sign(human_mean))
    verdict[agree] = "AGREE"
    verdict[disagree] = "DISAGREE"
    direction[agree & (mouse_mean > 0)] = "up/up"
    direction[agree & (mouse_mean < 0)] = "down/down"
    direction[disagree & (mouse_mean > 0)] = "up/down"
    direction[disagree & (mouse_mean < 0)] = "down/up"
    mouse_zero = cm_flat - mouse_up.ravel() - mouse_down.ravel()
    human_zero = ch_flat - human_up.ravel() - human_down.ravel()
    output = pd.DataFrame({
        "function_id": function_ids[fi],
        "function_name": functions["function_name"].to_numpy()[fi],
        "mouse_unit": pair_frame["mouse_unit"].to_numpy()[pi],
        "mouse_tissue": pair_frame["mouse_tissue"].to_numpy()[pi],
        "human_cohort": pair_frame["human_cohort"].to_numpy()[pi],
        "human_unit": pair_frame["human_unit"].to_numpy()[pi],
        "n_mouse_genes": cm_flat.astype(int),
        "n_human_genes": ch_flat.astype(int),
        "mouse_n_up": mouse_up.ravel().astype(int),
        "mouse_n_down": mouse_down.ravel().astype(int),
        "mouse_n_zero": mouse_zero.astype(int),
        "mouse_share_up": np.divide(mouse_up.ravel(), cm_flat, out=np.full(cm.size, np.nan), where=cm_flat > 0),
        "mouse_share_down": np.divide(mouse_down.ravel(), cm_flat, out=np.full(cm.size, np.nan), where=cm_flat > 0),
        "mouse_mean_log2FC": mouse_mean,
        "human_n_positive": human_up.ravel().astype(int),
        "human_n_negative": human_down.ravel().astype(int),
        "human_n_zero": human_zero.astype(int),
        "human_share_positive": np.divide(human_up.ravel(), ch_flat, out=np.full(ch.size, np.nan), where=ch_flat > 0),
        "human_share_negative": np.divide(human_down.ravel(), ch_flat, out=np.full(ch.size, np.nan), where=ch_flat > 0),
        "human_mean_slope": human_mean,
        "testable": testable,
        "mouse_noise_floor": mouse_thresholds,
        "human_noise_floor": human_thresholds,
        "verdict": verdict,
        "direction": direction,
    })
    assert len(output) == n_functions * 28
    assert not output.duplicated(["function_id"] + PAIR_KEYS).any()
    assert output.loc[~output["testable"], "verdict"].eq("UNTESTABLE").all()
    assert output.loc[output["n_mouse_genes"].eq(0), "mouse_mean_log2FC"].isna().all()
    assert output.loc[output["n_human_genes"].eq(0), "human_mean_slope"].isna().all()
    return output.sort_values(["function_id"] + PAIR_KEYS).reset_index(drop=True)


def summarize_observed(output: pd.DataFrame, label: str) -> tuple[dict, pd.DataFrame]:
    testable = output[output["testable"]].copy()
    pair_counts = testable.groupby("function_id")["mouse_unit"].size()
    if label == "repeatable":
        repeatable_ids = pair_counts[pair_counts >= 3].index
        testable = testable[testable["function_id"].isin(repeatable_ids)].copy()
    counts = agreement_counts(
        testable["mouse_mean_log2FC"].to_numpy(), testable["human_mean_slope"].to_numpy(),
        -np.inf, -np.inf,
    )
    unclear = testable["verdict"].eq("UNCLEAR")
    same = np.sign(testable["mouse_mean_log2FC"]) == np.sign(testable["human_mean_slope"])
    counts["n_agree"] = int((~unclear & same).sum())
    counts["n_disagree"] = int((~unclear & ~same).sum())
    counts["n_unclear"] = int(unclear.sum())
    counts["n_testable_rows"] = len(testable)
    decisive = counts["n_agree"] + counts["n_disagree"]
    counts["agree_per_testable"] = counts["n_agree"] / len(testable) if len(testable) else np.nan
    counts["directional_agreement"] = counts["n_agree"] / decisive if decisive else np.nan
    counts.update({"set": label, "n_functions": testable["function_id"].nunique()})
    return counts, testable


def make_summary(
    output: pd.DataFrame,
    null: dict[str, dict[str, np.ndarray]],
    mouse_floor: np.ndarray,
    human_floor: np.ndarray,
    permutations: int,
) -> pd.DataFrame:
    rows = []
    for label in ("all", "repeatable"):
        observed, testable = summarize_observed(output, label)
        row = dict(observed)
        row.update({
            "mouse_noise_floor_min": testable["mouse_noise_floor"].min(),
            "mouse_noise_floor_max": testable["mouse_noise_floor"].max(),
            "human_noise_floor_min": testable["human_noise_floor"].min(),
            "human_noise_floor_max": testable["human_noise_floor"].max(),
            "mouse_noise_floor_n5": mouse_floor[5],
            "human_noise_floor_n5": human_floor[5],
            "n_permutations": permutations,
            "seed": SEED,
        })
        for metric in ("agree_per_testable", "directional_agreement"):
            values = null[label][metric]
            finite = values[np.isfinite(values)]
            assert len(finite) > 0
            row[f"n_null_finite_{metric}"] = len(finite)
            row[f"null_mean_{metric}"] = float(finite.mean())
            row[f"null_q025_{metric}"] = float(np.quantile(finite, 0.025))
            row[f"null_q975_{metric}"] = float(np.quantile(finite, 0.975))
            row[f"p_upper_{metric}"] = float(
                (1 + np.sum(finite >= observed[metric])) / (len(finite) + 1)
            ) if np.isfinite(observed[metric]) else np.nan
        rows.append(row)
    columns = [
        "set", "n_functions", "n_testable_rows", "n_agree", "n_disagree", "n_unclear",
        "agree_per_testable", "null_mean_agree_per_testable", "null_q025_agree_per_testable",
        "null_q975_agree_per_testable", "p_upper_agree_per_testable", "n_null_finite_agree_per_testable",
        "directional_agreement", "null_mean_directional_agreement", "null_q025_directional_agreement",
        "null_q975_directional_agreement", "p_upper_directional_agreement", "n_null_finite_directional_agreement",
        "mouse_noise_floor_min", "mouse_noise_floor_max", "human_noise_floor_min", "human_noise_floor_max",
        "mouse_noise_floor_n5", "human_noise_floor_n5", "n_permutations", "seed",
    ]
    return pd.DataFrame(rows, columns=columns)


def pair_summary(matches: pd.DataFrame, output: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for match in matches.itertuples():
        subset = output[
            output["mouse_unit"].eq(match.mouse_unit)
            & output["human_cohort"].eq(match.human_cohort)
            & output["human_unit"].eq(match.human_unit)
        ]
        testable = subset[subset["testable"]]
        agree = int(testable["verdict"].eq("AGREE").sum())
        disagree = int(testable["verdict"].eq("DISAGREE").sum())
        unclear = int(testable["verdict"].eq("UNCLEAR").sum())
        rows.append({
            "mouse_unit": match.mouse_unit,
            "human_cohort": match.human_cohort,
            "human_unit": match.human_unit,
            "n_pair_genes": int(match.n_human_genes_tested),
            "n_testable_functions": len(testable),
            "n_agree": agree,
            "n_disagree": disagree,
            "n_unclear": unclear,
            "agree_per_testable": agree / len(testable) if len(testable) else np.nan,
            "directional_agreement": agree / (agree + disagree) if agree + disagree else np.nan,
        })
    result = pd.DataFrame(rows)
    assert len(result) == 28
    return result


def independent_spot_check(
    output: pd.DataFrame,
    pair_genes: pd.DataFrame,
    genes: pd.DataFrame,
    annotation: pd.DataFrame,
) -> dict:
    candidates = output[output["verdict"].isin(["AGREE", "DISAGREE"])].copy()
    if not len(candidates):
        candidates = output[output["testable"]].copy()
    assert len(candidates)
    candidates["total_genes"] = candidates["n_mouse_genes"] + candidates["n_human_genes"]
    row = candidates.sort_values(["total_genes", "function_id"] + PAIR_KEYS).iloc[0]
    source_genes = genes[genes["units_list"].astype(str).str.split(";").apply(
        lambda units: row.mouse_unit in units
    )].copy()
    de_path = RESULTS / f"followup-mature-p1m-de_{matched_test.safe_unit_name(row.mouse_unit)}.csv"
    de = pd.read_csv(de_path).set_index("index")
    source_genes["source_log2FC"] = source_genes["mouse_gene_symbol"].map(de["log2FC"])
    function_genes = set(annotation.loc[annotation["function_id"].eq(row.function_id), "mouse_gene_symbol"])
    mouse = source_genes[source_genes["mouse_gene_symbol"].isin(function_genes)]
    human = pair_genes[
        pair_genes["mouse_unit"].eq(row.mouse_unit)
        & pair_genes["human_cohort"].eq(row.human_cohort)
        & pair_genes["human_unit"].eq(row.human_unit)
        & pair_genes["mouse_gene_symbol"].isin(function_genes)
    ]
    assert len(mouse) == row.n_mouse_genes
    assert len(human) == row.n_human_genes
    np.testing.assert_allclose(mouse["source_log2FC"].mean(), row.mouse_mean_log2FC, rtol=0, atol=1e-12)
    np.testing.assert_allclose(human["human_slope"].mean(), row.human_mean_slope, rtol=0, atol=1e-12)
    assert int((mouse["source_log2FC"] > 0).sum()) == row.mouse_n_up
    assert int((human["human_slope"] > 0).sum()) == row.human_n_positive
    return {
        "function_id": row.function_id,
        "function_name": row.function_name,
        "mouse_unit": row.mouse_unit,
        "human_cohort": row.human_cohort,
        "human_unit": row.human_unit,
        "n_mouse": len(mouse),
        "n_human": len(human),
        "mouse_mean": float(mouse["source_log2FC"].mean()),
        "human_mean": float(human["human_slope"].mean()),
        "verdict": row.verdict,
    }


def fmt(value, digits=6) -> str:
    if value is None or pd.isna(value):
        return "NA"
    if isinstance(value, (bool, np.bool_)):
        return str(bool(value)).lower()
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.{digits}g}"
    return str(value).replace("|", "\\|")


def markdown_table(columns: list[str], rows: list[list], digits=6) -> str:
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
    lines.extend("| " + " | ".join(fmt(value, digits) for value in row) + " |" for row in rows)
    return "\n".join(lines)


def function_tables(output: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    testable = output[output["testable"]].copy()
    rows = []
    for function_id, frame in testable.groupby("function_id"):
        agree = frame[frame["verdict"].eq("AGREE")]
        disagree = frame[frame["verdict"].eq("DISAGREE")]
        directions = frame["direction"].value_counts()
        rows.append({
            "function_id": function_id,
            "function_name": frame["function_name"].iloc[0],
            "n_tested_pairs": len(frame),
            "n_agree": len(agree),
            "n_disagree": len(disagree),
            "n_unclear": int(frame["verdict"].eq("UNCLEAR").sum()),
            "agree_per_testable": len(agree) / len(frame),
            "directions": "; ".join(f"{key}:{value}" for key, value in sorted(directions.items()) if key),
        })
    functions = pd.DataFrame(rows)
    repeatable = functions[(functions["n_tested_pairs"] >= 3) & (functions["n_agree"] >= 3)].sort_values(
        ["n_agree", "n_tested_pairs", "function_id"], ascending=[False, False, True]
    )
    disagreeing = functions[functions["n_disagree"] > 0].sort_values(
        ["n_disagree", "n_tested_pairs", "function_id"], ascending=[False, False, True]
    ).head(15)
    return repeatable, disagreeing


def metadata_value(metadata: pd.DataFrame, key: str) -> str:
    selected = metadata.loc[metadata["key"].eq(key), "value"]
    return str(selected.iloc[0]) if len(selected) else "NA"


def write_report(
    output: pd.DataFrame,
    summary: pd.DataFrame,
    pairs: pd.DataFrame,
    metadata: pd.DataFrame,
    annotation: pd.DataFrame,
    genes: pd.DataFrame,
    spot: dict,
    permutations: int,
    started: str,
    duration: float,
) -> None:
    repeatable, disagreeing = function_tables(output)
    lines = [
        "# Pathway-level uniformity test, matched mouse-human P1m pairs", "",
        "## Numbers", "",
        "The requested row-weighted fraction is AGREE divided by all testable function-pair rows; UNCLEAR rows remain in that denominator. Directional agreement, AGREE/(AGREE+DISAGREE), excludes UNCLEAR and is shown alongside it. Null p values are upper-tail, plus-one corrected.", "",
        markdown_table(
            ["set", "functions", "testable rows", "AGREE", "DISAGREE", "UNCLEAR", "AGREE/testable", "null mean", "null 95% range", "p", "directional agreement", "null mean", "p"],
            [[r.set, r.n_functions, r.n_testable_rows, r.n_agree, r.n_disagree, r.n_unclear,
              r.agree_per_testable, r.null_mean_agree_per_testable,
              f"{fmt(r.null_q025_agree_per_testable)}–{fmt(r.null_q975_agree_per_testable)}",
              r.p_upper_agree_per_testable, r.directional_agreement,
              r.null_mean_directional_agreement, r.p_upper_directional_agreement]
             for r in summary.itertuples()],
        ), "",
        f"Noise floors are gene-count matched: for a row with n genes on a side, its threshold is the 95th percentile of absolute permuted means among all side-testable rows with that n (nearest available n if absent). Across observed testable rows, mouse thresholds range {fmt(summary['mouse_noise_floor_min'].iloc[0])}–{fmt(summary['mouse_noise_floor_max'].iloc[0])} log2FC and human thresholds range {fmt(summary['human_noise_floor_min'].iloc[0])}–{fmt(summary['human_noise_floor_max'].iloc[0])} slope units. At n=5 the thresholds are {fmt(summary['mouse_noise_floor_n5'].iloc[0])} and {fmt(summary['human_noise_floor_n5'].iloc[0])}, respectively. Every row's exact thresholds are in the main CSV; estimates use {permutations:,} shuffles.", "",
        "### Repeatable functions with at least three AGREE pairs", "",
        "The repeatable analysis contains functions testable in at least three pairs. The table below further requires at least three AGREE pairs; directions are mouse/human.", "",
        markdown_table(
            ["GO id", "name", "tested pairs", "AGREE", "DISAGREE", "UNCLEAR", "AGREE/testable", "directions"],
            [[r.function_id, r.function_name, r.n_tested_pairs, r.n_agree, r.n_disagree,
              r.n_unclear, r.agree_per_testable, r.directions] for r in repeatable.itertuples()],
        ) if len(repeatable) else "None.", "",
        "### Top disagreeing functions", "",
        markdown_table(
            ["GO id", "name", "tested pairs", "AGREE", "DISAGREE", "UNCLEAR", "directions"],
            [[r.function_id, r.function_name, r.n_tested_pairs, r.n_agree, r.n_disagree,
              r.n_unclear, r.directions] for r in disagreeing.itertuples()],
        ) if len(disagreeing) else "None.", "",
        "### Per-pair summary", "",
        markdown_table(
            ["mouse unit", "human cohort", "human unit", "pair genes", "testable functions", "AGREE", "DISAGREE", "UNCLEAR", "AGREE/testable", "directional agreement"],
            [[r.mouse_unit, r.human_cohort, r.human_unit, r.n_pair_genes, r.n_testable_functions,
              r.n_agree, r.n_disagree, r.n_unclear, r.agree_per_testable, r.directional_agreement]
             for r in pairs.itertuples()],
        ), "",
        "## Interpretation", "",
    ]
    overall = summary.iloc[0]
    relation = "above" if overall.agree_per_testable > overall.null_mean_agree_per_testable else "below or equal to"
    lines.extend([
        f"Observed AGREE/testable is {overall.agree_per_testable:.6f}, {relation} the permutation mean of {overall.null_mean_agree_per_testable:.6f}; the upper-tail p value is {overall.p_upper_agree_per_testable:.6g}. The corresponding repeatable-set values are {summary.iloc[1].agree_per_testable:.6f}, {summary.iloc[1].null_mean_agree_per_testable:.6f}, and {summary.iloc[1].p_upper_agree_per_testable:.6g}. These are descriptive function-level direction comparisons, not evidence that GO annotations themselves are causal aging mechanisms.", "",
        "## Design and grouping", "",
        f"The committed P1m derivation supplied {len(genes):,} genes and 36 eligible mouse units; only the 28 matched pairs in `results/followup-mature-matched-p1m.csv` were analyzed, without floored/unfloored arms. Mouse effects are per-unit 30-month versus 18-month log2FC values from the existing DE CSVs. The P1m list was not re-derived. The cached 299 pair-gene rows cover every one of the 22 pairs with tested human genes; the other six pairs have zero tested genes. All 299 rows were checked against freshly computed slopes from the owner-approved human cache.", "",
        f"GO grouping used direct database mapping through `org.Mm.eg.db` `GOALL`, filtered to `ONTOLOGYALL == BP`, with names from `GO.db`; this includes direct and inherited GO Biological Process annotations. clusterProfiler was not used. The mapping contains {annotation['mouse_gene_symbol'].nunique():,} annotated mouse genes, {annotation['function_id'].nunique():,} GO terms, and {len(annotation):,} gene-term memberships. A gene can belong to many terms. The CSV is the complete {annotation['function_id'].nunique():,} terms × 28 pairs grid; absent directions are NA, never zero-filled.", "",
        "Human orthologs and collision handling are exactly those of the matched test: use `human_ortholog` from the P1m gene CSV, sort within each mouse unit by human ortholog, unit FDR, and mouse symbol, and keep the first mouse gene per human ortholog. Human slopes are continuous-age OLS fits for donors aged 30+ on donor-cell-type log1p counts-per-10,000 pseudobulks, with sex included where variable and identifiable. Mouse n counts all annotated P1m DE genes in that unit; human n counts annotated, ortholog-mapped genes present with finite slopes in that human unit.", "",
        f"A row is testable only with at least {MIN_GENES} genes on each side. Otherwise it is UNTESTABLE. For testable rows, either absolute mean below its noise floor makes the row UNCLEAR; otherwise equal mean signs give AGREE and opposite signs DISAGREE. Up/down and positive/negative counts and shares are recorded in the main CSV.", "",
        f"The null uses {permutations:,} deterministic shuffles (seed {SEED}) of gene labels relative to the complete GO membership table. Each shuffle preserves every gene's number of GO memberships, every GO term's global size, and the linkage between each gene's mouse and human measurements; the same shuffle is applied across all 28 pairs. Counts, means, noise-floor classifications, and both agreement fractions are recomputed. The repeatable null reselects functions testable in at least three pairs within each shuffle. Noise floors are fixed at the gene-count-matched pooled null 95th percentiles for observed and null classifications.", "",
        "## Verification", "",
        f"The match check found {len(pairs)} unique matched pairs. The independent spot-check selected `{spot['function_id']}` ({spot['function_name']}) in `{spot['mouse_unit']}` versus `{spot['human_cohort']}::{spot['human_unit']}`. Recomputing from the unit DE CSV, P1m gene list, GO mapping, and pair-gene slopes gave mouse n={spot['n_mouse']}, mean={spot['mouse_mean']:.12g}; human n={spot['n_human']}, mean={spot['human_mean']:.12g}; verdict={spot['verdict']}. Counts, means, and direction counts exactly match the serialized main CSV.", "",
        "After writing, the script re-read both result CSVs and recomputed the pooled and per-pair tables used in this report from those files. Assertions cover 28-pair completeness, unique function-pair keys, NA means at zero counts, cached pair counts, and exact observed/null summary consistency.", "",
        "## Limits", "",
        "- GO terms overlap extensively and are not independent tests. Row weighting gives more influence to terms and pairs with many annotations, and the pooled fractions are not gene- or pair-independent estimates.",
        "- The 95th-percentile noise floors adapt to gene count but pool across functions and pairs; they do not adapt to function identity, tissue, effect distribution, or donor number. UNCLEAR is a threshold classification, not a confidence interval or a multiple-testing-adjusted pathway p value.",
        "- Permutation tests exchangeability of gene-to-function labels while preserving the measured values and annotation multiplicity. It does not model donor sampling, cell-type mismatches, mouse gene-selection uncertainty, GO annotation error, or dependence among terms.",
        "- Only 22 pairs have any tested human genes, and gene coverage is highly unequal. Six zero-gene pairs cannot contribute testable functions. The repeatable set reduces sparsity but does not remove unequal tissue, cohort, donor, or gene-count representation.",
        "- Mouse effects are selected 18m-to-30m two-group log2FC values; human effects are cross-sectional continuous-age slopes on a different scale. Equal signs are not equivalent effect sizes or proof of conserved causal programs. The mouse 30m group is all male, so mouse age and sex remain confounded.",
        "- Same-named cell types are not proven homologous across species. Existing ortholog and exact-name matching choices propagate into every function result.", "",
        "## Reproducibility and manuscript boundary", "",
        f"Ran with Python {sys.version.split()[0]}, pandas {pd.__version__}, NumPy {np.__version__}, SciPy {scipy.__version__}; {metadata_value(metadata, 'R_version')}; org.Mm.eg.db {metadata_value(metadata, 'org.Mm.eg.db_package_version')}; GO.db {metadata_value(metadata, 'GO.db_package_version')}; AnnotationDbi {metadata_value(metadata, 'AnnotationDbi_package_version')}. Database source metadata follows.", "",
        markdown_table(["metadata key", "value"], [[r.key, r.value] for r in metadata.itertuples()], digits=4), "",
        f"Script: `scripts/run_mature_matched_p1m_pathways.py`. Use `--refresh-human-cache` to rebuild the external cache and pair-gene rows. Started `{started}`; runtime {duration:.3f} seconds. Main data: `{OUTPUT_CSV.relative_to(REPO)}`; pooled data: `{SUMMARY_CSV.relative_to(REPO)}`. Human cache and GO mapping are under `{WORK_DIR}`.", "",
        "No manuscript or submission file was changed. Any TeX change is only a proposal for later owner review: if these descriptive results are cited, report both the observed and null fractions, the noise-floor rule, and the limits above rather than claiming pathway conservation from agreement alone.", "",
    ])
    OUTPUT_REPORT.write_text("\n".join(lines), encoding="utf-8")


def verify_serialized(
    output: pd.DataFrame,
    summary: pd.DataFrame,
    matches: pd.DataFrame,
    permutations: int,
) -> pd.DataFrame:
    reread = pd.read_csv(OUTPUT_CSV, low_memory=False)
    reread["direction"] = reread["direction"].fillna("")
    pd.testing.assert_frame_equal(output, reread, check_dtype=False, check_exact=False, rtol=1e-12, atol=1e-12)
    reread_summary = pd.read_csv(SUMMARY_CSV)
    pd.testing.assert_frame_equal(summary, reread_summary, check_dtype=False, check_exact=False, rtol=1e-12, atol=1e-12)
    for label in ("all", "repeatable"):
        observed, _ = summarize_observed(reread, label)
        stored = reread_summary.set_index("set").loc[label]
        for key in ("n_testable_rows", "n_agree", "n_disagree", "n_unclear", "agree_per_testable", "directional_agreement"):
            np.testing.assert_allclose(observed[key], stored[key], rtol=0, atol=1e-12)
        assert int(stored["n_permutations"]) == permutations
    return pair_summary(matches, reread)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--permutations", type=int, default=2000)
    parser.add_argument("--refresh-human-cache", action="store_true")
    args = parser.parse_args()
    assert args.permutations >= 2000, "Owner order requires at least 2,000 GO-label permutations"
    start = time.monotonic()
    started = utc_now()
    print(f"Matched P1m pathway analysis started {started}", flush=True)
    matches = load_matches()
    _, mouse_units, genes = matched_test.load_mouse()
    pair_genes = load_pair_genes(matches, mouse_units, genes, args.refresh_human_cache)
    annotation, metadata = build_go_mapping(genes)
    print(f"GO BP mapping: {len(annotation)} edges, {annotation['function_id'].nunique()} terms", flush=True)
    functions, membership, mouse_values, human_values = build_observed(
        matches, mouse_units, pair_genes, genes, annotation
    )
    mouse_floor, human_floor, null = permutation_analysis(
        membership, mouse_values, human_values, args.permutations
    )
    print(f"n=5 noise floors: mouse={mouse_floor[5]:.12g}; human={human_floor[5]:.12g}", flush=True)
    output = aggregate_observed(
        matches, functions, membership, mouse_values, human_values, mouse_floor, human_floor
    )
    summary = make_summary(output, null, mouse_floor, human_floor, args.permutations)
    output.to_csv(OUTPUT_CSV, index=False)
    summary.to_csv(SUMMARY_CSV, index=False)
    spot = independent_spot_check(output, pair_genes, genes, annotation)
    pairs = verify_serialized(output, summary, matches, args.permutations)
    write_report(
        output, summary, pairs, metadata, annotation, genes, spot,
        args.permutations, started, time.monotonic() - start,
    )
    print(summary.to_string(index=False), flush=True)
    print("Spot check: " + json.dumps(spot), flush=True)
    print(f"Wrote {OUTPUT_CSV} ({len(output)} rows)", flush=True)
    print(f"Wrote {SUMMARY_CSV} ({len(summary)} rows)", flush=True)
    print(f"Wrote {OUTPUT_REPORT}", flush=True)
    print(f"Pathway analysis completed in {time.monotonic() - start:.3f} seconds", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
