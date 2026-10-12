#!/usr/bin/env python3
"""
scripts/run_all.py — Master End-to-End Execution Script for Conserved Aging Explorer
Executes all data verification, pipeline scripts, notebook nbconverts, and app asset exports.
"""

import os
import sys
import time
import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
PYTHON_BIN = sys.executable

def run_step(step_name, cmd, cwd=REPO_ROOT):
    print(f"\n{'='*70}")
    print(f"▶ STARTING: {step_name}")
    print(f"  Command: {' '.join(cmd)}")
    print(f"{'='*70}\n", flush=True)
    t0 = time.time()
    res = subprocess.run(cmd, cwd=cwd)
    dt = time.time() - t0
    if res.returncode != 0:
        print(f"\n❌ FAILED: {step_name} (Exit code {res.returncode}) after {dt:.1f}s", flush=True)
        sys.exit(res.returncode)
    print(f"\n✅ PASSED: {step_name} in {dt:.1f}s\n", flush=True)
    return dt

def main():
    t_start = time.time()
    runtimes = {}
    print("="*70)
    print("CONSERVED AGING EXPLORER — FULL PIPELINE REPRODUCTION RUNNER")
    print(f"Repository Root: {REPO_ROOT}")
    print(f"Python Binary:   {PYTHON_BIN}")
    print("="*70)

    # 1. Verification of Data Checksums
    runtimes["M0: Data Verification"] = run_step(
        "M0: Verify Data Checksums",
        [PYTHON_BIN, "scripts/verify_data.py"]
    )

    # 2. Milestone M1: QC, Normalization, Clustering
    runtimes["M1: Pipeline QC & Clustering"] = run_step(
        "M1: Run Pipeline QC & Clustering",
        [PYTHON_BIN, "scripts/run_m1_pipeline.py"]
    )
    runtimes["M1: Notebook 01 nbconvert"] = run_step(
        "M1: Notebook 01 nbconvert",
        [PYTHON_BIN, "-m", "jupyter", "nbconvert", "--to", "notebook", "--execute", "--inplace", "notebooks/01-qc-and-clustering.ipynb"]
    )

    # 3. Milestone M2: Differential Expression & Calibration Trio
    runtimes["M2: Differential Expression"] = run_step(
        "M2: Run Differential Expression & Ranking",
        [PYTHON_BIN, "scripts/run_m2_pipeline.py"]
    )
    runtimes["M2: Notebook 02 nbconvert"] = run_step(
        "M2: Notebook 02 nbconvert",
        [PYTHON_BIN, "-m", "jupyter", "nbconvert", "--to", "notebook", "--execute", "--inplace", "notebooks/02-de-aging.ipynb"]
    )

    # 4. Milestone M3: Trajectory & Literature Agreement
    runtimes["M3: Trajectory & Agreement"] = run_step(
        "M3: Run Trajectory Analysis & Agreement",
        [PYTHON_BIN, "scripts/run_m3_pipeline.py"]
    )
    runtimes["M3: Notebook 03 nbconvert"] = run_step(
        "M3: Notebook 03 nbconvert",
        [PYTHON_BIN, "-m", "jupyter", "nbconvert", "--to", "notebook", "--execute", "--inplace", "notebooks/03-trajectory-agreement.ipynb"]
    )

    # 5. Milestone M4: Human Cohorts & Orthologs
    runtimes["M4: Human Orthologs & Power"] = run_step(
        "M4: Prepare Orthologs & Power Calculations",
        [PYTHON_BIN, "scripts/prepare_m4_orthologs_and_power.py"]
    )
    runtimes["M4: Notebook 04 nbconvert"] = run_step(
        "M4: Notebook 04 nbconvert",
        [PYTHON_BIN, "-m", "jupyter", "nbconvert", "--to", "notebook", "--execute", "--inplace", "notebooks/04-human-prep.ipynb"]
    )

    # 6. Milestone M5: Cross-Species Conservation & Divergence
    runtimes["M5: Conservation & Divergence"] = run_step(
        "M5: Cross-Species Conservation Pipeline",
        [PYTHON_BIN, "scripts/run_m5_pipeline.py"]
    )
    runtimes["M5: Notebook 05 nbconvert"] = run_step(
        "M5: Notebook 05 nbconvert",
        [PYTHON_BIN, "-m", "jupyter", "nbconvert", "--to", "notebook", "--execute", "--inplace", "notebooks/05-conservation.ipynb"]
    )

    # 7. Milestone M6: Longevity Outliers (NMR)
    runtimes["M6: Longevity Outlier Pipeline"] = run_step(
        "M6: Longevity Outlier Pipeline",
        [PYTHON_BIN, "scripts/run_m6_pipeline.py"]
    )
    runtimes["M6: Notebook 06 nbconvert"] = run_step(
        "M6: Notebook 06 nbconvert",
        [PYTHON_BIN, "-m", "jupyter", "nbconvert", "--to", "notebook", "--execute", "--inplace", "notebooks/06-outliers.ipynb"]
    )

    # 8. Milestone M8: Export & Check App Assets
    runtimes["M8: Export App Assets"] = run_step(
        "M8: Export App Assets to Parquet & JSON",
        [PYTHON_BIN, "scripts/export_app_assets.py"]
    )
    runtimes["M8: Check App Assets"] = run_step(
        "M8: Verify App Assets Manifest & Budgets",
        [PYTHON_BIN, "scripts/check_app_assets.py"]
    )

    total_time = time.time() - t_start
    print("\n" + "="*70)
    print("🎉 FULL PIPELINE EXECUTION COMPLETED SUCCESSFULLY!")
    print(f"Total Execution Time: {total_time/60:.2f} minutes ({total_time:.1f}s)")
    print("="*70)
    print("\nPer-Stage Hardware Runtimes:")
    for k, v in runtimes.items():
        print(f"  - {k:<38}: {v:.1f}s")
    print("="*70)

if __name__ == "__main__":
    main()
