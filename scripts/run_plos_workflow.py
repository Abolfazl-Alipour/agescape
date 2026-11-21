#!/usr/bin/env python3
"""Verify or execute the ordered PLOS ONE workflow manifest."""

import argparse
import csv
import shlex
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
MANIFEST = REPO / "report" / "plos-one" / "workflow-manifest.tsv"
FIELDS = ["order", "stage", "command", "inputs", "outputs", "seeds", "environment"]


def read_manifest():
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != FIELDS:
            raise ValueError(f"Manifest columns must be {FIELDS}")
        rows = list(reader)
    for expected, row in enumerate(rows, 1):
        if int(row["order"]) != expected:
            raise ValueError(f"Nonsequential order at manifest row {expected}")
        if any(not row[field].strip() for field in FIELDS):
            raise ValueError(f"Empty manifest field at stage {expected}")
    return rows


def commands_for(row):
    tokens = shlex.split(row["command"])
    groups = []
    current = []
    for token in tokens:
        if token == "&&":
            if not current:
                raise ValueError(f"Empty command in stage {row['stage']}")
            groups.append(current)
            current = []
        else:
            current.append(token.replace("{python}", sys.executable))
    if not current:
        raise ValueError(f"Empty command in stage {row['stage']}")
    groups.append(current)
    for command in groups:
        for part in command[1:]:
            if part.endswith(".py") and not (REPO / part).is_file():
                raise FileNotFoundError(f"Missing analysis script: {part}")
    return groups


def verify_outputs(rows):
    checked = 0
    for row in rows:
        for output in row["outputs"].split(";"):
            path = REPO / output.strip()
            if not path.exists():
                raise FileNotFoundError(f"Missing output for {row['stage']}: {output.strip()}")
            checked += 1
    return checked


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execute", action="store_true", help="Run every stage in order after verification")
    args = parser.parse_args()
    rows = read_manifest()
    commands = [commands_for(row) for row in rows]
    checked = verify_outputs(rows)
    print(f"PASS: {len(rows)} ordered stages, {checked} existing output paths, and all analysis scripts verified")
    if args.execute:
        for row, groups in zip(rows, commands):
            print(f"START {row['order']}: {row['stage']}", flush=True)
            for command in groups:
                subprocess.run(command, cwd=REPO, check=True)
            print(f"DONE {row['order']}: {row['stage']}", flush=True)
        print("PASS: all workflow stages completed")


if __name__ == "__main__":
    main()
