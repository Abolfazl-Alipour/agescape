#!/usr/bin/env python3
"""Aggregate PBMC donor-level models and draw separate-cohort forest plots."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"
FIGURES = REPO / "figures"
PREFIX = "followup-pbmc"
COHORTS = ["GSE213516", "GSE214546"]
DATA = REPO.parent / "agescape-revision-data" / "pbmc"


def read(cohort: str, suffix: str) -> pd.DataFrame:
    return pd.read_csv(RESULTS / f"{PREFIX}-{cohort}-{suffix}")


def concatenate(suffix: str, output: str) -> pd.DataFrame:
    frame = pd.concat([read(cohort, suffix) for cohort in COHORTS], ignore_index=True, sort=False)
    frame.to_csv(RESULTS / output, index=False)
    return frame


def forest(frame: pd.DataFrame, path: Path, title: str, expected: dict | None = None) -> None:
    frame = frame.copy()
    frame = frame[frame["status"].eq("fit") & frame["common_cell_type"].notna()].copy()
    if frame.empty:
        raise ValueError(f"No fitted harmonized rows for {path.name}")
    if "plot_label" in frame:
        frame["label"] = frame["plot_label"] + " | " + frame["cohort"]
    else:
        frame["label"] = frame["common_cell_type"] + " | " + frame["cohort"]
    frame = frame.sort_values("label", ascending=False).reset_index(drop=True)
    height = max(3.5, 0.42 * len(frame) + 1.8)
    fig, ax = plt.subplots(figsize=(10.5, height))
    colors = {"GSE213516": "#0369a1", "GSE214546": "#b45309"}
    y = np.arange(len(frame))
    for index, row in frame.iterrows():
        color = colors[row["cohort"]]
        ax.plot([row["ci95_low"], row["ci95_high"]], [index, index], color=color, lw=1.8)
        ax.scatter(row["effect"], index, color=color, s=38, zorder=3)
        ax.text(row["ci95_high"], index + 0.22, f"q={row['q_bh']:.3g}", fontsize=7, va="bottom")
    ax.axvline(0, color="black", lw=0.8, ls="--")
    ax.set_yticks(y, frame["label"], fontsize=8)
    ax.set_xlabel("Effect (95% CI); old minus young")
    ax.set_title(title)
    ax.grid(axis="x", alpha=0.2)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def targeted_forest(pseudobulk: pd.DataFrame, composition: pd.DataFrame) -> None:
    pseudo = pseudobulk[pseudobulk["design"].eq("strict_group") & pseudobulk["model"].eq("primary")].copy()
    pseudo["common_cell_type"] = pseudo["gene"]
    pseudo["plot_label"] = pseudo["gene"]
    comp = composition[composition["design"].eq("strict_group") & composition["model"].eq("primary")].copy()
    comp["plot_label"] = "CLR: " + comp["common_cell_type"]
    combined = pd.concat([pseudo, comp], ignore_index=True, sort=False)
    forest(combined, FIGURES / f"{PREFIX}-targeted-forest.png",
           "Targeted PBMC replication families (strict primary)")


def shared_primary(frame: pd.DataFrame, output: Path, title: str) -> None:
    part = frame[frame["design"].eq("strict_group") & frame["model"].eq("primary")].copy()
    shared = set(part.groupby("common_cell_type")["cohort"].nunique().loc[lambda x: x == 2].index) - {""}
    part = part[part["common_cell_type"].isin(shared)]
    forest(part, output, title)


def main() -> None:
    FIGURES.mkdir(exist_ok=True)
    p1 = concatenate("p1-models.csv", f"{PREFIX}-p1-models.csv")
    pseudo = concatenate("pseudobulk-targeted.csv", f"{PREFIX}-pseudobulk-targeted.csv")
    comp = concatenate("clr-targeted.csv", f"{PREFIX}-clr-targeted.csv")
    hell = concatenate("hellinger-models.csv", f"{PREFIX}-hellinger-models.csv")
    exploratory_clr = concatenate("clr-models.csv", f"{PREFIX}-clr-models.csv")
    labels = concatenate("label-inventory.csv", f"{PREFIX}-label-harmonization.csv")
    eligibility = pd.concat([
        pd.concat([read(c, "eligibility.csv") for c in COHORTS], ignore_index=True),
        pd.concat([read(c, "hellinger-eligibility.csv") for c in COHORTS], ignore_index=True),
    ], ignore_index=True, sort=False)
    eligibility.to_csv(RESULTS / f"{PREFIX}-eligibility-all.csv", index=False)
    identifiability = pd.concat([read(c, "covariate-identifiability.csv") for c in COHORTS], ignore_index=True)
    identifiability.to_csv(RESULTS / f"{PREFIX}-covariate-identifiability.csv", index=False)
    strict = pd.concat([
        p1[p1["design"].eq("strict_group") & p1["model"].eq("primary")],
        pseudo[pseudo["design"].eq("strict_group") & pseudo["model"].eq("primary")],
        comp[comp["design"].eq("strict_group") & comp["model"].eq("primary")],
        exploratory_clr[exploratory_clr["design"].eq("strict_group") & exploratory_clr["model"].eq("primary")],
        hell[hell["design"].eq("strict_group") & hell["model"].eq("primary")],
    ], ignore_index=True, sort=False)
    strict.to_csv(RESULTS / f"{PREFIX}-replication-strict.csv", index=False)
    sensitivities = pd.concat([
        p1[~p1["design"].eq("strict_group")],
        pseudo[~pseudo["design"].eq("strict_group")],
        comp[~comp["design"].eq("strict_group")],
        exploratory_clr[~exploratory_clr["design"].eq("strict_group")],
        hell[~hell["design"].eq("strict_group")],
    ], ignore_index=True, sort=False)
    sensitivities.to_csv(RESULTS / f"{PREFIX}-replication-sensitivities.csv", index=False)
    targeted_forest(pseudo, comp)
    shared_primary(p1, FIGURES / f"{PREFIX}-p1-forest.png", "P1 transport, shared released-label mappings (strict primary)")
    shared_primary(hell, FIGURES / f"{PREFIX}-hellinger-forest.png", "Cross-fitted Hellinger, shared released-label mappings (strict depth-adjusted)")
    survivor_frames = []
    for cohort in COHORTS:
        path = RESULTS / f"{PREFIX}-{cohort}-pseudobulk-exploratory.csv.gz"
        for chunk in pd.read_csv(path, chunksize=200000):
            survivors = chunk[chunk["q_bh"] < 0.05].copy()
            if not survivors.empty:
                survivor_frames.append(survivors)
    survivors = pd.concat(survivor_frames, ignore_index=True) if survivor_frames else pd.DataFrame()
    survivors.to_csv(RESULTS / f"{PREFIX}-pseudobulk-exploratory-survivors.csv.gz", index=False, compression="gzip")
    resources, validation_summary = {}, {}
    for cohort in COHORTS:
        verification = json.loads((RESULTS / f"{PREFIX}-{cohort}-validation.json").read_text())
        inventory = json.loads((DATA / "interim" / cohort / "download-inventory.json").read_text())
        resources[cohort] = {
            **verification["resource_use"],
            "download_files": len(inventory["files"]),
            "download_bytes": sum(file["size_bytes"] for file in inventory["files"]),
            "donor_aggregate_cache_bytes": verification["cache_bytes"],
        }
        validation_summary[cohort] = {
            "donors": len(verification["donors"]),
            "raw_cells": sum(donor["raw_barcodes"] for donor in verification["donors"]),
            "annotation_rows": sum(donor["annotation_rows"] for donor in verification["donors"]),
            "annotation_rows_matched": sum(donor["annotation_rows_matched"] for donor in verification["donors"]),
            "unmatched_annotation_rows_excluded": sum(donor.get("unmatched_annotation_rows_excluded", 0) for donor in verification["donors"]),
            "all_counts_finite_nonnegative_integer": all(donor["all_stored_counts_finite_nonnegative_integer"] for donor in verification["donors"]),
            "all_barcode_joins_one_to_one": all(donor["barcode_join_one_to_one"] for donor in verification["donors"]),
            "duplicate_symbol_feature_rows": verification["duplicate_symbol_feature_rows"],
        }
        if cohort == "GSE213516":
            validation_summary[cohort]["genome_reference_adjudication"] = verification["genome_reference_adjudication"]
    (RESULTS / f"{PREFIX}-resource-use.json").write_text(json.dumps(resources, indent=2) + "\n")
    (RESULTS / f"{PREFIX}-validation-summary.json").write_text(json.dumps(validation_summary, indent=2) + "\n")
    print(f"Wrote replication tables, three forest plots, {len(survivors)} exploratory survivor rows", flush=True)


if __name__ == "__main__":
    main()
