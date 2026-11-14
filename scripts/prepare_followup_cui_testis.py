#!/usr/bin/env python3
"""Build the Cui testis donor map, public inventory, and blocked endpoint table."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path
import re

from openpyxl import load_workbook

REPO = Path(__file__).resolve().parent.parent
DATA = REPO.parent / "agescape-revision-data" / "cui-testis"
METADATA = DATA / "metadata"
RESULTS = REPO / "results"
PREFIX = "followup-cui-testis"
EXTERNAL_DONORS = [1, 2, 3, 26, 27, 28, 29, 31, 32, 33]
SOMATIC = [
    "Sertoli", "TPC", "Leydig", "SMC", "Endothelial", "Macrophage", "Lymphocyte", "B",
]
GERM = [
    "Undifferentiated spermatogonia", "Differentiating spermatogonia",
    "Early primary spermatocytes", "Late primary spermatocytes",
    "Round spermatids", "Elongating spermatids",
]
CAVEAT = (
    "Changing germ-cell-stage composition can mimic broad expression displacement; "
    "testis hits/nulls cannot support general whole-body aging claims."
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def external_path(path) -> str:
    path = Path(path)
    if not path.is_absolute():
        return str(path)
    return str(Path("..") / path.relative_to(REPO.parent))


def parse_soft(path: Path) -> tuple[dict, dict]:
    donors: dict[int, dict] = {}
    series: dict[str, list[str] | str] = {}
    current = None
    with gzip.open(path, "rt") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("!Series_overall_design = "):
                series["overall_design"] = line.split(" = ", 1)[1]
            elif line.startswith("!Series_summary = "):
                series["summary"] = line.split(" = ", 1)[1]
            elif line.startswith("^SAMPLE = "):
                current = {"gsm": line.split(" = ", 1)[1], "files": {}}
            elif current is not None and line.startswith("!Sample_title = Donor"):
                current["donor_number"] = int(line.rsplit("Donor", 1)[1])
                donors[current["donor_number"]] = current
            elif current is not None and line.startswith("!Sample_characteristics_ch1 = age (years): "):
                current["geo_age"] = int(line.rsplit(" ", 1)[1])
            elif current is not None and line.startswith("!Sample_data_processing = Assembly: "):
                current["assembly"] = line.rsplit(" ", 1)[1]
            elif current is not None and line.startswith("!Sample_supplementary_file_"):
                url = line.split(" = ", 1)[1].replace("ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov")
                filename = url.rsplit("/", 1)[-1]
                kind = next(kind for kind in ["matrix", "barcodes", "features"] if f"_{kind}" in filename)
                current["files"][kind] = {"url": url, "file": filename}
    if len(donors) != 25:
        raise ValueError(f"Expected 25 GEO donors, found {len(donors)}")
    return donors, series


def parse_filelist(path: Path) -> dict[str, int]:
    sizes = {}
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["Name"].endswith((".mtx.gz", ".tsv.gz")):
                sizes[row["Name"]] = int(row["Size"])
    if len(sizes) != 75:
        raise ValueError(f"Expected 75 per-donor files, found {len(sizes)}")
    return sizes


def primary_donor_table(path: Path) -> list[dict]:
    sheet = load_workbook(path, read_only=True, data_only=True)["F1a"]
    rows = []
    for donor_number, values in enumerate(list(sheet.iter_rows(values_only=True))[1:], start=1):
        if values[1] is None:
            continue
        age = int(str(values[1]).removesuffix("yr"))
        rows.append({
            "donor_number": donor_number,
            "donor": f"Donor{donor_number}",
            "primary_age": age,
            "primary_decade": f"{age // 10 * 10}s",
            "bmi": float(values[2]),
        })
    if len(rows) != 35 or [row["donor_number"] for row in rows] != list(range(1, 36)):
        raise ValueError("Fig. 1a source is not the expected ordered 35-donor table")
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows for {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    soft = METADATA / "geo" / "GSE254315_family.soft.gz"
    filelist = METADATA / "geo" / "filelist.txt"
    fig1 = METADATA / "supplements" / "43587_2025_824_MOESM4_ESM.xlsx"
    ed1 = METADATA / "supplements" / "43587_2025_824_MOESM10_ESM.xlsx"
    plan_path = REPO / "data" / "manifest" / "cui_testis_analysis_plan.json"
    geo, series = parse_soft(soft)
    sizes = parse_filelist(filelist)
    primary = primary_donor_table(fig1)
    included = [row for row in primary if row["donor_number"] in geo]
    exact_matches = sum(row["primary_age"] == geo[row["donor_number"]]["geo_age"] for row in included)
    if exact_matches != 24:
        raise ValueError(f"Expected 24 exact primary/GEO age matches, found {exact_matches}")
    discrepancies = [row["donor_number"] for row in included if row["primary_age"] != geo[row["donor_number"]]["geo_age"]]
    if discrepancies != [8] or geo[8]["geo_age"] != 39 or primary[7]["primary_age"] != 40:
        raise ValueError(f"Unexpected age discrepancy: {discrepancies}")

    manifest_path = REPO / "data" / "manifest" / "GSE254315.json"
    previous = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    previous_files = {entry["file"]: entry for entry in previous.get("files", [])}
    files = []
    donor_rows = []
    manifest_donors = []
    for row in primary:
        number = row["donor_number"]
        is_included = number in geo
        record = {**row, "included_in_gse254315_step": is_included}
        if is_included:
            current = geo[number]
            library = re.search(r"_([^_]+)_matrix", current["files"]["matrix"]["file"]).group(1)
            record.update({
                "gsm": current["gsm"], "library": library,
                "geo_age": current["geo_age"], "geo_decade": f"{current['geo_age'] // 10 * 10}s",
                "assembly": current["assembly"],
                "age_discrepancy": number == 8,
            })
            for kind in ["matrix", "barcodes", "features"]:
                source = current["files"][kind]
                old = previous_files.get(source["file"], {})
                files.append({
                    **source, "donor": row["donor"], "gsm": current["gsm"], "kind": kind,
                    "size_bytes": sizes[source["file"]],
                    "local_path": external_path(old["local_path"]) if old.get("local_path") else "",
                    "sha256": old.get("sha256", ""),
                    "download_verified": old.get("download_verified", False),
                })
            manifest_donors.append(record)
        donor_rows.append({
            "donor": row["donor"], "scope": "GSE254315" if is_included else "external_atlas_only",
            "included_in_inferential_n": int(is_included), "gsm": record.get("gsm", ""),
            "library": record.get("library", ""), "primary_age": row["primary_age"],
            "geo_age": record.get("geo_age", ""), "primary_decade": row["primary_decade"],
            "geo_decade": record.get("geo_decade", ""), "bmi": row["bmi"],
            "age_discrepancy": number == 8,
            "mapping_note": "separate donor library" if is_included else "reused external atlas donor; excluded",
        })

    plan = json.loads(plan_path.read_text())
    manifest = {
        "accession": "GSE254315",
        "paper": {
            "title": "Single-cell transcriptomic atlas of the human testis across the reproductive lifespan",
            "doi": "10.1038/s43587-025-00824-2", "pmcid": "PMC12003174",
            "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12003174/",
        },
        "created_from_frozen_plan": str(plan_path.relative_to(REPO)),
        "analysis_plan_sha256": sha256(plan_path),
        "series_summary": series.get("summary", ""),
        "series_overall_design": series.get("overall_design", ""),
        "mapping_resolution": {
            "gse254315_independent_donor_libraries": 25,
            "integrated_paper_donors": 35,
            "external_atlas_only_donors": EXTERNAL_DONORS,
            "external_source_accessions_stated_by_paper": ["GSE182786", "GSE215754"],
            "external_donors_excluded": True,
            "primary_geo_exact_age_matches": exact_matches,
            "primary_geo_age_discrepancy": {"donor": "Donor8", "primary_age": 40, "geo_age": 39},
            "pooled_or_demultiplexed_libraries": False,
            "maximum_independent_n": 25,
            "mapping_status": "resolved_for_all_25_GSE254315_donors",
        },
        "annotation_gate": {
            "geo_files": "Raw counts, barcodes, and features only.",
            "primary_fig1_source": "Donor order, exact age, and BMI; no barcode labels.",
            "primary_extended_data_fig1_source": "Age-aggregated cell-label totals; no donor-and-barcode crosswalk.",
            "released_code_repository": "https://github.com/JingtaoLab/Testicular_Aging; notebooks refer to local metadata/objects but release no annotation files.",
            "public_barcode_to_cell_type_crosswalk": False,
            "status": "blocked_for_all_cell_resolved_expression_and_composition_models",
        },
        "sources": [
            {"kind": "geo_soft", "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/GSE254315/soft/GSE254315_family.soft.gz", "path": external_path(soft), "sha256": sha256(soft)},
            {"kind": "geo_filelist", "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE254nnn/GSE254315/suppl/filelist.txt", "path": external_path(filelist), "sha256": sha256(filelist)},
            {"kind": "primary_fig1_source", "url": "https://www.ebi.ac.uk/biostudies/files/S-EPMC12003174/43587_2025_824_MOESM4_ESM.xlsx", "path": external_path(fig1), "sha256": sha256(fig1)},
            {"kind": "primary_extended_data_fig1_source", "url": "https://www.ebi.ac.uk/biostudies/files/S-EPMC12003174/43587_2025_824_MOESM10_ESM.xlsx", "path": external_path(ed1), "sha256": sha256(ed1)},
        ],
        "donors": manifest_donors,
        "files": files,
        "biological_caveat": CAVEAT,
        "resource_guard": plan["resources"],
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    write_csv(RESULTS / f"{PREFIX}-donor-map.csv", donor_rows)
    write_csv(RESULTS / f"{PREFIX}-file-inventory.csv", files)

    status_rows = []
    for compartment, labels in [("somatic", SOMATIC), ("germ_cell_stage", GERM)]:
        for label in labels:
            for endpoint in ["P1", "donor_pseudobulk", "CLR", "Hellinger"]:
                status_rows.append({
                    "compartment": compartment, "cell_type": label, "endpoint": endpoint,
                    "status": "blocked", "reason": "no_released_donor_barcode_cell_type_crosswalk",
                    "n_donors": 0, "effect": "", "ci95_low": "", "ci95_high": "",
                    "pvalue": "", "q_bh": "", "bh_family": "not_run",
                    "caveat": CAVEAT,
                })
    write_csv(RESULTS / f"{PREFIX}-endpoint-status.csv", status_rows)

    sheet = load_workbook(ed1, read_only=True, data_only=True).active
    released = []
    for values in list(sheet.iter_rows(values_only=True))[1:]:
        if values[1] is None:
            continue
        released.append({
            "source_row": values[0], "released_label": values[1], "released_age": values[2],
            "released_cells": values[3], "age_total_cells": values[4], "fraction": values[5],
            "donor_crosswalk_available": False,
            "use": "descriptive age-aggregated released totals only; not independent donor observations",
        })
    for item in released:
        item["caveat"] = CAVEAT
    write_csv(RESULTS / f"{PREFIX}-released-age-cell-totals.csv", released)

    fig6 = METADATA / "supplements" / "43587_2025_824_MOESM9_ESM.xlsx"
    if fig6.exists():
        sheet = load_workbook(fig6, read_only=True, data_only=True)["F6a"]
        by_age_bmi = {(row["primary_age"], row["bmi"]): row for row in primary}
        if len(by_age_bmi) != 35:
            raise ValueError("Primary age/BMI pairs are not unique donor keys")
        released_freq = {}
        for values in list(sheet.iter_rows(values_only=True))[1:]:
            if values[0] is None:
                continue
            key = (int(str(values[0]).removesuffix("yr")), float(values[1]))
            donor = by_age_bmi[key]
            if donor["donor_number"] in released_freq:
                raise ValueError("Duplicate Fig. 6a age/BMI donor key")
            released_freq[donor["donor_number"]] = values[2]
        if set(released_freq) != set(range(1, 36)) - {34}:
            raise ValueError("Unexpected Fig. 6a donor coverage")
        spermatogenesis_rows = []
        for donor in primary:
            if donor["donor_number"] not in geo:
                continue
            missing = donor["donor_number"] not in released_freq
            spermatogenesis_rows.append({
                "donor": donor["donor"], "gsm": geo[donor["donor_number"]]["gsm"],
                "primary_age": donor["primary_age"], "bmi": donor["bmi"],
                "released_elongating_spermatid_freq": released_freq.get(donor["donor_number"], ""),
                "released_freq_status": "missing_from_Fig6a" if missing else "released",
                "incomplete_spermatogenesis_status": "not_released",
                "donor_pas_sperm_count_per_tubule": "not_released",
                "inference": "descriptive only; no spermatogenesis status inferred from composition",
                "caveat": CAVEAT,
            })
        write_csv(RESULTS / f"{PREFIX}-spermatogenesis-metadata.csv", spermatogenesis_rows)
        manifest["sources"].append({
            "kind": "primary_fig6_source",
            "url": "https://www.ebi.ac.uk/biostudies/files/S-EPMC12003174/43587_2025_824_MOESM9_ESM.xlsx",
            "path": external_path(fig6), "sha256": sha256(fig6),
        })
        manifest["spermatogenesis_gate"] = {
            "released_donor_elongating_spermatid_freq": True,
            "fig6a_missing_donor": "Donor34",
            "released_donor_incomplete_spermatogenesis_status": False,
            "released_donor_pas_sperm_count_per_tubule": False,
            "rule": "Do not infer incomplete spermatogenesis from germ-cell composition; adjustment and complete-only sensitivities are untestable.",
        }
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Resolved 25 GSE donors plus {len(EXTERNAL_DONORS)} excluded external atlas donors")
    print(f"Wrote {manifest_path} and {len(status_rows)} explicitly blocked endpoint rows")


if __name__ == "__main__":
    main()
