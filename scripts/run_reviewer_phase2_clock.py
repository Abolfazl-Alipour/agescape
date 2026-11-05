#!/usr/bin/env python3

from pathlib import Path

import pandas as pd
from openpyxl import load_workbook

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "results"


def main():
    path = REPO / "data" / "interim" / "universal-clock-supplementary-data-1-14.xlsx"
    workbook = load_workbook(path, read_only=True, data_only=True)
    needles = ["GSE199979", "GSM6000839", "H1228_B6D2F1_CD_F"]
    counts = {needle: 0 for needle in needles}
    for sheet in workbook.sheetnames:
        for row in workbook[sheet].iter_rows(values_only=True):
            for value in row:
                text = str(value)
                for needle in needles:
                    counts[needle] += needle in text
    rows = [{
        "dataset": "GSE199979", "evaluation_samples": 339,
        "source_checked": "https://static-content.springer.com/esm/art%3A10.1038%2Fs43587-023-00462-6/MediaObjects/43587_2023_462_MOESM3_ESM.xlsx",
        "supplementary_sheets_checked": len(workbook.sheetnames),
        "accession_mentions": counts["GSE199979"], "tested_sample_id_mentions": counts["GSM6000839"],
        "tested_sample_title_mentions": counts["H1228_B6D2F1_CD_F"],
        "training_overlap_status": "sample-level training membership unresolved; benchmark reproduction, not independent validation",
    }]
    pd.DataFrame(rows).to_csv(RESULTS / "reviewer-phase2-clock-overlap-audit.csv", index=False)


if __name__ == "__main__":
    main()
