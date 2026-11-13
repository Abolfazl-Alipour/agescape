#!/usr/bin/env python3
"""Build verified De Man donor, sample, and cell-type eligibility tables."""

import gzip
import hashlib
import json
import tarfile
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
DATA = REPO.parent / f"{REPO.name}-data" / "de-man"
RESULTS = REPO / "results"
MANIFEST = REPO / "data" / "manifest"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 22), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    subjects = pd.read_csv(DATA / "41467_2026_68810_MOESM3_ESM.csv").rename(
        columns={"Sample ID": "sample_id", "Dataset": "source", "Sex": "sex", "Age": "age",
                 "Smoking": "smoking"})
    subjects["source"] = subjects["source"].map({"Dataset 1": "Baylor", "Dataset 2": "IPF_Cell_Atlas"})
    assert subjects["sample_id"].is_unique
    assert len(subjects) == 60
    assert (subjects["source"] == "Baylor").sum() == 32
    assert (subjects["source"] == "IPF_Cell_Atlas").sum() == 28
    assert subjects["age"].between(11, 80).all()

    ipf = pd.read_csv(DATA / "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz", sep="\t")
    ipf = ipf[ipf["Disease_Identity"] == "Control"].copy()
    ipf["qc_gene_pass"] = ipf["nGene"].between(750, 7500)
    ipf["subject"] = ipf["Subject_Identity"].astype(str)
    assert set(ipf["subject"]) == set(subjects.loc[subjects["source"] == "IPF_Cell_Atlas", "sample_id"])
    assert ipf["CellBarcode_Identity"].is_unique

    cell_counts = (ipf.groupby(["subject", "Manuscript_Identity"], observed=True)
                   .agg(n_cells=("CellBarcode_Identity", "size"),
                        n_qc_gene_pass=("qc_gene_pass", "sum")).reset_index())
    cell_counts.to_csv(RESULTS / "followup-de-man-celltype-counts.csv", index=False)

    subject_counts = (ipf.groupby("subject", observed=True)
                      .agg(n_cells=("CellBarcode_Identity", "size"),
                           n_qc_gene_pass=("qc_gene_pass", "sum"),
                           n_libraries=("Library_Identity", "nunique"))
                      .reset_index().rename(columns={"subject": "sample_id"}))
    crosswalk = subjects.merge(subject_counts, on="sample_id", how="left")
    crosswalk[["n_cells", "n_qc_gene_pass", "n_libraries"]] = crosswalk[
        ["n_cells", "n_qc_gene_pass", "n_libraries"]].fillna(0).astype(int)
    crosswalk["annotation_status"] = crosswalk["source"].map(
        {"Baylor": "not_released_in_GEO_raw_tar", "IPF_Cell_Atlas": "released"})
    crosswalk.to_csv(RESULTS / "followup-de-man-donor-crosswalk.csv", index=False)

    summary = crosswalk.groupby("source").agg(
        n_subjects=("sample_id", "size"), age_min=("age", "min"), age_max=("age", "max"),
        age_sd=("age", "std"), n_cells=("n_cells", "sum"),
        n_qc_gene_pass=("n_qc_gene_pass", "sum"), n_libraries=("n_libraries", "sum"))
    summary["age_11_to_80"] = True
    summary["annotation_note"] = [
        "Raw matrices downloaded; cell annotations and integrated object not present in GEO raw tar",
        "Released control annotations present; libraries pooled by subject for inference",
    ]
    summary.to_csv(RESULTS / "followup-de-man-eligibility-summary.csv")

    files = {
        "GSE281219_RAW.tar": DATA / "GSE281219_RAW.tar",
        "GSE136831_RawCounts_Sparse.mtx.gz": DATA / "GSE136831_RawCounts_Sparse.mtx.gz",
        "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz": DATA / "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz",
        "41467_2026_68810_MOESM3_ESM.csv": DATA / "41467_2026_68810_MOESM3_ESM.csv",
    }
    provenance = {name: {"size_bytes": path.stat().st_size, "sha256": sha256(path)}
                  for name, path in files.items()}
    (MANIFEST / "de_man_sources.json").write_text(json.dumps({
        "downloaded": "2026-11-13",
        "paper": "https://doi.org/10.1038/s41467-026-68810-9",
        "github_tree": "0c11d44ce4103451327fadd106e1e6872edefa50",
        "subjects": {"Baylor": 32, "IPF_Cell_Atlas_subjects": 28, "IPF_Cell_Atlas_libraries": 38},
        "files": provenance,
        "baylor_annotation_status": "not_released_in_GEO_raw_tar",
    }, indent=2) + "\n")
    print(summary.to_string())


if __name__ == "__main__":
    main()
