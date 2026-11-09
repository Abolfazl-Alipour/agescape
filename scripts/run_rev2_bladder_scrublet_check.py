#!/usr/bin/env python3
"""Bladder scrublet sensitivity for the Rev2 concordance redo.

Question: does omitting scrublet (the --no-scrublet deviation in
run_rev2_concordance.py, forced by an OOM kill in this shared VM) move the
Bladder unit age slopes? Bladder is the smallest tissue (8,702 cells after
QC, 33 committed doublets = 0.38% of cells), so if slopes barely move here,
the deviation is bounded everywhere it matters.

Runs in three separate processes because holding two Bladder loads at once
OOM-killed the combined version in this shared VM:
  --side on      Bladder WITH scrublet (M1 settings, rate 0.06) -> saves slopes
  --side off     Bladder WITHOUT scrublet -> saves slopes
  --side compare reads both slope files and writes the CSV

Slopes are the full-design age slopes (const + age + sex + mean_log_umi)
for the 4 Phase 1 Bladder units, on all genes after the technical mask.
Writes results/reviewer-rev2-bladder-scrublet-check.csv (compare step only).
CPU-only; no new data.
"""

import argparse
import gc
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib.util
_spec = importlib.util.spec_from_file_location(
    "run_rev2_concordance", Path(__file__).resolve().parent / "run_rev2_concordance.py")
rc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rc)

RESULTS = rc.RESULTS
SCRATCH = Path.home() / "workspace" / "rev2-scratch"


def load_bladder(use_scrublet):
    """Return {cell_type: unit}, genes, qc info for Bladder (one load only)."""
    gen = rc.mouse_strata(use_scrublet=use_scrublet, tissues_only=["Bladder"])
    tissue, genes, units, info = next(gen)
    return {u["cell_type"]: u for u in units}, genes, info


def save_slopes(side):
    use_scrublet = side == "on"
    units, genes, info = load_bladder(use_scrublet)
    phase1 = pd.read_csv(RESULTS / "reviewer-phase1-universality-units.csv")
    bladder_units = sorted({r.cell_type for r in phase1[
        (phase1["species"] == "Mouse") & (phase1["tissue"] == "Bladder")].itertuples()})
    out = {"cell_types": np.array(bladder_units, dtype=object),
           "n_donors": {}, "slopes": {}, "n_genes": {}}
    for ct in bladder_units:
        unit = units.get(ct)
        if unit is None:
            continue
        # Full gene vector (all genes after the technical mask): using a
        # per-unit detection mask here would misalign the two sides, since
        # the mask itself depends on which cells scrublet removed.
        out["n_donors"][ct] = len(unit["metadata"])
        out["n_genes"][ct] = unit["expression"].shape[1]
        out["slopes"][ct] = rc.unit_slopes(unit)
        del unit
    gc.collect()
    path = SCRATCH / f"bladder_slopes_scrublet_{side}.npz"
    np.savez(path, cell_types=out["cell_types"],
             n_donors=np.array([out["n_donors"].get(ct, 0) for ct in bladder_units]),
             n_genes=np.array([out["n_genes"].get(ct, 0)
                               for ct in bladder_units]),
             **{f"slope_{ct}": out["slopes"][ct] for ct in out["slopes"]})
    print(f"Bladder scrublet {side.upper()}: n_after_qc={info['n_after_qc']} "
          f"doublets={info['n_doublets']} final={info['n_final']}; "
          f"saved {len(out['slopes'])} unit slope vectors to {path}", flush=True)


def compare():
    from scipy.stats import spearmanr
    on = np.load(SCRATCH / "bladder_slopes_scrublet_on.npz", allow_pickle=True)
    off = np.load(SCRATCH / "bladder_slopes_scrublet_off.npz", allow_pickle=True)
    rows = []
    for i, ct in enumerate(on["cell_types"]):
        key = f"slope_{ct}"
        if key not in on or key not in off:
            rows.append({"tissue": "Bladder", "cell_type": str(ct),
                         "n_donors_scrublet_on": int(on["n_donors"][i]),
                         "n_donors_scrublet_off": int(off["n_donors"][i]),
                         "n_genes": int(on["n_genes"][i]),
                         "slope_pearson": float("nan"), "slope_spearman": float("nan"),
                         "mean_abs_slope_diff": float("nan"),
                         "max_abs_slope_diff": float("nan"),
                         "slope_sd_scrublet_on": float("nan"),
                         "note": "unit missing on one side; not comparable"})
            continue
        s_on, s_off = on[key], off[key]
        diff = np.abs(s_on - s_off)
        rows.append({
            "tissue": "Bladder", "cell_type": str(ct),
            "n_donors_scrublet_on": int(on["n_donors"][i]),
            "n_donors_scrublet_off": int(off["n_donors"][i]),
            "n_genes": int(on["n_genes"][i]),
            "slope_pearson": float(np.corrcoef(s_on, s_off)[0, 1]),
            "slope_spearman": float(spearmanr(s_on, s_off).statistic),
            "mean_abs_slope_diff": float(diff.mean()),
            "max_abs_slope_diff": float(diff.max()),
            "slope_sd_scrublet_on": float(np.std(s_on)),
            "note": "Slopes on all genes after the technical mask; diff in log1p(CP10K) per month.",
        })
        print(f"  {ct}: Pearson={rows[-1]['slope_pearson']:.4f} "
              f"Spearman={rows[-1]['slope_spearman']:.4f} "
              f"mean|diff|={rows[-1]['mean_abs_slope_diff']:.5f} "
              f"max|diff|={rows[-1]['max_abs_slope_diff']:.5f} "
              f"(slope SD {rows[-1]['slope_sd_scrublet_on']:.4f})", flush=True)
    frame = pd.DataFrame(rows)
    rc.emit(frame, "reviewer-rev2-bladder-scrublet-check.csv")
    rc.robust_write(frame, RESULTS / "reviewer-rev2-bladder-scrublet-check.csv")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--side", choices=["on", "off", "compare"], required=True)
    args = parser.parse_args()
    if args.side == "compare":
        compare()
    else:
        save_slopes(args.side)


if __name__ == "__main__":
    main()
