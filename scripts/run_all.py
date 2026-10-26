#!/usr/bin/env python3
"""
scripts/run_all.py — master end-to-end runner for the (rebuilt, real-data) project.

Executes every stage in order with full recompute: no stage skips its work because
an output already exists. Each stage's measured wall time is logged honestly at the
end. Any non-zero exit aborts the run with a non-zero exit code.

Stage order:
  1.  scripts/verify_data.py
  2.  scripts/run_m1_pipeline.py
  3.  scripts/run_m2_pipeline.py
  4.  scripts/run_m3_pipeline.py
  5.  scripts/run_m3_agreement.py
  6.  scripts/run_m5_pipeline.py
  7.  scripts/run_m6_pipeline.py
  8.  scripts/fetch_nmr_gene_names.py   (ONLY if data/interim/nmr-ensg2symbol.tsv is missing)
  9.  nbconvert --execute notebooks 01..06 (kernelspec agescape, timeout 600)
  10. scripts/check_app_assets.py
"""

import subprocess
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
PYTHON_BIN = sys.executable
NMR_SYMBOL_TSV = REPO_ROOT / "data" / "interim" / "nmr-ensg2symbol.tsv"

NOTEBOOKS = [
    "notebooks/01-qc-and-clustering.ipynb",
    "notebooks/02-de-aging.ipynb",
    "notebooks/03-trajectory-agreement.ipynb",
    "notebooks/04-human-prep.ipynb",
    "notebooks/05-conservation.ipynb",
    "notebooks/06-outliers.ipynb",
]


def run_step(step_name, cmd):
    """Run one stage; print start/end and measured runtime. Abort the run on failure."""
    print(f"\n{'=' * 72}")
    print(f"START: {step_name}")
    print(f"  cmd: {' '.join(str(c) for c in cmd)}")
    print(f"{'=' * 72}", flush=True)
    t0 = time.time()
    res = subprocess.run([str(c) for c in cmd], cwd=REPO_ROOT)
    dt = time.time() - t0
    if res.returncode != 0:
        print(f"\nFAILED: {step_name} (exit {res.returncode}) after {dt:.1f}s", flush=True)
        print("Aborting run_all.py; no later stage was executed.", flush=True)
        sys.exit(res.returncode if res.returncode != 0 else 1)
    print(f"\nPASSED: {step_name} in {dt:.1f}s", flush=True)
    return dt


def nbconvert_cmd(notebook):
    return [
        PYTHON_BIN, "-m", "jupyter", "nbconvert",
        "--to", "notebook", "--execute", "--inplace",
        "--ExecutePreprocessor.kernel_name=agescape",
        "--ExecutePreprocessor.timeout=600",
        notebook,
    ]


def main():
    runtimes = {}

    print("=" * 72)
    print("FULL PIPELINE RUNNER — real-data rebuild (takeover 2026-10-14)")
    print(f"Repo root: {REPO_ROOT}")
    print(f"Python:    {PYTHON_BIN}")
    print("Every stage recomputes from scratch; nothing is skipped for existing outputs.")
    print("=" * 72, flush=True)

    runtimes["verify_data"] = run_step(
        "Data verification (SHA-256 vs manifests)",
        [PYTHON_BIN, "scripts/verify_data.py"],
    )
    runtimes["m1_pipeline"] = run_step(
        "M1 pipeline (QC, clustering, annotation)",
        [PYTHON_BIN, "scripts/run_m1_pipeline.py"],
    )
    runtimes["m2_pipeline"] = run_step(
        "M2 pipeline (pseudobulk DE, ranking, calibration trio)",
        [PYTHON_BIN, "scripts/run_m2_pipeline.py"],
    )
    runtimes["m3_pipeline"] = run_step(
        "M3 pipeline (DPT trajectory, bridge model)",
        [PYTHON_BIN, "scripts/run_m3_pipeline.py"],
    )
    runtimes["m3_agreement"] = run_step(
        "M3 agreement vs TMS paper published coefficients",
        [PYTHON_BIN, "scripts/run_m3_agreement.py"],
    )
    runtimes["m5_pipeline"] = run_step(
        "M5 pipeline (mouse program catalog, divergence status)",
        [PYTHON_BIN, "scripts/run_m5_pipeline.py"],
    )
    runtimes["m6_pipeline"] = run_step(
        "M6 pipeline (lifespan program scores, NMR contrast)",
        [PYTHON_BIN, "scripts/run_m6_pipeline.py"],
    )

    if NMR_SYMBOL_TSV.exists():
        print(f"\nSKIP (conditional stage): nmr-ensg2symbol.tsv already present at "
              f"{NMR_SYMBOL_TSV}", flush=True)
        runtimes["fetch_nmr_gene_names"] = 0.0
    else:
        runtimes["fetch_nmr_gene_names"] = run_step(
            "Fetch NMR gene symbols from Ensembl REST (cached TSV missing)",
            [PYTHON_BIN, "scripts/fetch_nmr_gene_names.py"],
        )

    for nb in NOTEBOOKS:
        name = Path(nb).stem
        runtimes[f"nb_{name}"] = run_step(f"Execute {nb}", nbconvert_cmd(nb))

    runtimes["export_app_data"] = run_step(
        "Export app data (Spec v2)",
        [PYTHON_BIN, "scripts/export_app_data.py"],
    )

    runtimes["check_app_assets"] = run_step(
        "Check app assets (146 anchors)",
        [PYTHON_BIN, "scripts/check_app_assets.py"],
    )

    runtimes["h1_epigenetics"] = run_step(
        "H1: Pan-mammalian epigenetic clock & drift (Universal Clock 2)",
        [PYTHON_BIN, "scripts/run_h1_epigenetic_drift.py"],
    )

    runtimes["h2_composition"] = run_step(
        "H2: Compositional remodeling (Aitchison CLR)",
        [PYTHON_BIN, "scripts/run_h2_composition.py"],
    )

    runtimes["h3_entropy"] = run_step(
        "H3: Single-cell transcriptional noise & entropy (1,500-UMI rarefaction)",
        [PYTHON_BIN, "scripts/run_h3_entropy.py"],
    )

    runtimes["fibroblast_markers"] = run_step(
        "Fibroblast candidate markers audit (11 markers, FDR BH)",
        [PYTHON_BIN, "scripts/run_fibroblast_markers.py"],
    )

    print("\n" + "=" * 72)
    print("ALL STAGES PASSED. Measured stage runtimes (wall seconds, this machine):")
    for stage, dt in runtimes.items():
        print(f"  {stage:<28} {dt:>10.1f}s")
    print("=" * 72, flush=True)


if __name__ == "__main__":
    main()
