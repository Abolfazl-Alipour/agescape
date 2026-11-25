"""Re-file of the 605-test family recompute (originally /tmp/recompute_605.py, lost).

Recomputes the 605-test full-range family from the committed 624-test file:
drop the 15 zebrafish and 4 naked mole rat rows, then Benjamini-Hochberg
correction over the remaining 605 p-values. Also recomputes the 427-test
core family (mouse + human rows) the same way.

Checks, all against committed files:
  1. This BH implementation reproduces the committed q_global_bh_624 column
     over all 624 p-values (max abs diff).
  2. It reproduces the committed q_global_bh_431 column over the 431 core
     rows (mouse + human + naked mole rat).
  3. The recomputed 605 rows and q_global_bh_605 column match the committed
     results/rev3-composition-global-bh-605.csv exactly (row keys and q
     values), with 45 passes and the same 45 rows.
  4. The 427-family BH has 15 passes, the same rows as the 431-family's
     15 passes.

Usage: python3 scripts/recompute_605_family.py [--write PATH]
  --write writes the recomputed 605 CSV to PATH (default: checks only).
Exit code 0 only if every check passes.
"""

import argparse
import csv
import sys

REPO_624 = "results/rev3-composition-global-bh-624.csv"
REPO_605 = "results/rev3-composition-global-bh-605.csv"
DROP_SPECIES = {"Zebrafish", "Naked mole-rat"}
CORE_SPECIES = {"Mouse", "Human", "Naked mole-rat"}
CORE_427_SPECIES = {"Mouse", "Human"}
TOL = 1e-12

KEY_COLS = ["species", "organ", "cohort", "cell_type", "sensitivity",
            "model", "effect_type"]


def bh_qvalues(pvalues):
    """Standard Benjamini-Hochberg q-values (as in statsmodels fdrcorrection).

    q_i = p_(i) * n / rank_i, then a running minimum from the largest p
    downward so q is non-decreasing with p, clipped at 1.
    """
    n = len(pvalues)
    order = sorted(range(n), key=lambda i: pvalues[i])
    q = [0.0] * n
    running_min = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        val = pvalues[i] * n / rank
        if val < running_min:
            running_min = val
        q[i] = min(running_min, 1.0)
    return q


def row_key(row):
    return tuple(row[c] for c in KEY_COLS)


def load(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def max_diff(computed, committed):
    return max(abs(a - b) for a, b in zip(computed, committed))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", default=None,
                    help="also write the recomputed 605 CSV to this path")
    args = ap.parse_args()

    rows624 = load(REPO_624)
    assert len(rows624) == 624, f"expected 624 rows, got {len(rows624)}"
    p624 = [float(r["pvalue"]) for r in rows624]

    failures = []

    # Check 1: BH reproduces the committed 624 q column.
    q624 = bh_qvalues(p624)
    d1 = max_diff(q624, [float(r["q_global_bh_624"]) for r in rows624])
    print(f"check 1: BH vs committed q_global_bh_624, max abs diff {d1:.3e}")
    if d1 > TOL:
        failures.append("check 1")

    # Check 2: BH reproduces the committed 431 q column on core rows.
    core_idx = [i for i, r in enumerate(rows624)
                if r["species"] in CORE_SPECIES]
    assert len(core_idx) == 431, f"expected 431 core rows, got {len(core_idx)}"
    q431 = bh_qvalues([p624[i] for i in core_idx])
    committed431 = []
    for i in core_idx:
        val = rows624[i]["q_global_bh_431"]
        committed431.append(float(val) if val.strip() else None)
    pairs = [(a, b) for a, b in zip(q431, committed431) if b is not None]
    d2 = max_diff([a for a, _ in pairs], [b for _, b in pairs])
    print(f"check 2: BH vs committed q_global_bh_431 ({len(pairs)} rows), "
          f"max abs diff {d2:.3e}")
    if d2 > TOL:
        failures.append("check 2")

    # Recompute the 605 family: drop zebrafish and naked mole rat rows.
    keep_idx = [i for i, r in enumerate(rows624)
                if r["species"] not in DROP_SPECIES]
    dropped = [rows624[i]["species"] for i in range(624) if i not in keep_idx]
    assert len(keep_idx) == 605, f"expected 605 rows, got {len(keep_idx)}"
    assert dropped.count("Zebrafish") == 15, dropped
    assert dropped.count("Naked mole-rat") == 4, dropped
    q605 = bh_qvalues([p624[i] for i in keep_idx])
    recomputed = [dict(rows624[i]) for i in keep_idx]
    for row, q in zip(recomputed, q605):
        row["q_global_bh_605"] = repr(q)

    # Check 3: recomputed 605 matches the committed 605 CSV.
    rows605 = load(REPO_605)
    assert len(rows605) == 605, f"expected 605 rows, got {len(rows605)}"
    committed_by_key = {row_key(r): r for r in rows605}
    assert len(committed_by_key) == 605, "duplicate keys in committed 605 CSV"
    key_mismatch = 0
    q_diffs = []
    for row in recomputed:
        committed = committed_by_key.get(row_key(row))
        if committed is None:
            key_mismatch += 1
            continue
        q_diffs.append(abs(float(row["q_global_bh_605"])
                           - float(committed["q_global_bh_605"])))
    d3 = max(q_diffs) if q_diffs else float("inf")
    passes_recomputed = {row_key(r) for r in recomputed
                         if float(r["q_global_bh_605"]) < 0.05}
    passes_committed = {row_key(r) for r in rows605
                        if float(r["q_global_bh_605"]) < 0.05}
    print(f"check 3: recomputed vs committed 605 CSV, key mismatches "
          f"{key_mismatch}, max q diff {d3:.3e}, passes "
          f"{len(passes_recomputed)} vs {len(passes_committed)}, "
          f"same rows: {passes_recomputed == passes_committed}")
    if key_mismatch or d3 > TOL or passes_recomputed != passes_committed:
        failures.append("check 3")

    # Check 4: 427-family BH, 15 passes, same rows as the 431 family's.
    idx427 = [i for i, r in enumerate(rows624)
              if r["species"] in CORE_427_SPECIES]
    assert len(idx427) == 427, f"expected 427 rows, got {len(idx427)}"
    q427 = bh_qvalues([p624[i] for i in idx427])
    passes427 = {row_key(rows624[i]) for i, q in zip(idx427, q427)
                 if q < 0.05}
    passes431 = {row_key(rows624[i]) for i in core_idx
                 if rows624[i]["q_global_bh_431"].strip()
                 and float(rows624[i]["q_global_bh_431"]) < 0.05}
    print(f"check 4: 427-family passes {len(passes427)}, 431-family passes "
          f"{len(passes431)}, same rows: {passes427 == passes431}")
    if len(passes427) != 15 or passes427 != passes431:
        failures.append("check 4")

    if args.write:
        with open(args.write, "w", newline="") as f:
            fieldnames = list(rows624[0].keys()) + ["q_global_bh_605"]
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(recomputed)
        print(f"wrote {args.write}")

    if failures:
        print("FAIL: " + ", ".join(failures))
        return 1
    print("ALL CHECKS PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
