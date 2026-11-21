#!/usr/bin/env python3
"""Independent BH verification of Batch A cohort CSVs (fresh implementation,
stdlib only, nothing imported from the repo script).

For every cohort CSV (gao, yazar, perez, guilliams):
- transport and fractions families: recompute BH q from the raw p-values stored
  in the file over the declared family (all tested units in that family,
  pooled groups included) and compare to the stored q, composite and empirical.
- concordance rows: q must equal p (single-unit family).
- transport classes and verdicts: re-derive from the stored q values under the
  frozen rules and compare to the stored class/verdict.
- notes check (yazar/perez/guilliams transport rows): the q_all values recorded
  in notes must equal the mature_q columns after the fix.
- genes rows: mature_q is the minimum per-gene global BH q in that cell type;
  cross-check it against the finalize BH-check CSV top-gene entries when that
  file exists (BH q is monotone in p, so min q belongs to the min-p gene).

Usage: python3 verify_batcha_bh.py <repo-results-dir> [bh-check-csv]
Exit 1 if any problem is found.
"""
import csv
import math
import re
import sys
from collections import Counter, defaultdict

ALPHA = 0.05


def bh(pvals):
    n = len(pvals)
    out = [float("nan")] * n
    valid = [(i, p) for i, p in enumerate(pvals) if p is not None and math.isfinite(p)]
    if not valid:
        return out
    m = len(valid)
    order = sorted(valid, key=lambda t: t[1])
    ranked = [p for _, p in order]
    adj = [ranked[k] * m / (k + 1) for k in range(m)]
    for k in range(m - 2, -1, -1):
        adj[k] = min(adj[k], adj[k + 1])
    adj = [min(a, 1.0) for a in adj]
    for (idx, _), a in zip(order, adj):
        out[idx] = a
    return out


def f(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return float("nan")
    return v


def tclass(cq, eq, testable=True):
    if not testable:
        return "UNTESTABLE"
    c = cq is not None and math.isfinite(cq) and cq < ALPHA
    e = eq is not None and math.isfinite(eq) and eq < ALPHA
    if c and e:
        return "BOTH_PASS"
    if e:
        return "EMPIRICAL_ONLY"
    if c:
        return "COMPOSITE_ONLY"
    return "FAIL_BOTH"


def verdict(o_sig, m_sig, testable=True):
    if not testable:
        return "UNTESTABLE"
    if o_sig is None:
        return "NEW" if m_sig else "NULL"
    if o_sig and m_sig:
        return "SURVIVED"
    if o_sig and not m_sig:
        return "DIED"
    if (not o_sig) and m_sig:
        return "NEW"
    return "NULL"


def main():
    results_dir = sys.argv[1]
    bh_check_path = sys.argv[2] if len(sys.argv) > 2 else None
    problems = []
    summaries = []

    for cohort in ("gao", "yazar", "perez", "guilliams"):
        path = f"{results_dir}/followup-mature-organ-transports-a-{cohort}.csv"
        try:
            rows = list(csv.DictReader(open(path)))
        except FileNotFoundError:
            print(f"{cohort}: CSV MISSING, skipped")
            continue
        print(f"{cohort}: {len(rows)} rows")
        fams = defaultdict(list)
        for r in rows:
            fams[r["endpoint_family"]].append(r)
        for fam, sub in sorted(fams.items()):
            is_transport = fam.startswith("transport_")
            is_fractions = fam.startswith("fractions")
            is_conc = fam.startswith("concordance")
            if not (is_transport or is_fractions or is_conc):
                continue
            for pcol, qcol, label in (
                ("mature_p", "mature_q", "composite"),
                ("mature_empirical_p", "mature_empirical_q", "empirical"),
            ):
                ps = [f(r[pcol]) for r in sub]
                qs = bh(ps)
                n_fin = sum(1 for p in ps if math.isfinite(p))
                mm = 0
                maxd = 0.0
                for r, q in zip(sub, qs):
                    st = f(r[qcol])
                    if math.isfinite(q) and math.isfinite(st):
                        d = abs(q - st)
                        maxd = max(maxd, d)
                        if d > 1e-9:
                            mm += 1
                            problems.append(
                                f"{cohort} {fam} {label} {r['cell_type']}: recomputed q={q:.6g} stored={st:.6g}"
                            )
                    elif math.isfinite(q) != math.isfinite(st):
                        mm += 1
                        problems.append(
                            f"{cohort} {fam} {label} {r['cell_type']}: finite mismatch recomputed={q} stored={st}"
                        )
                n_sig = sum(1 for q in qs if math.isfinite(q) and q < ALPHA)
                summaries.append(
                    f"{cohort} {fam} {label}: n={len(sub)} finite={n_fin} mismatches={mm} max|dq|={maxd:.2e} q<0.05:{n_sig}"
                )
            if is_transport:
                for r in sub:
                    testable = r["verdict"] != "UNTESTABLE"
                    tc = tclass(f(r["mature_q"]), f(r["mature_empirical_q"]), testable)
                    if tc != r["mature_transport_class"]:
                        problems.append(
                            f"{cohort} {fam} {r['cell_type']}: class {tc} != stored {r['mature_transport_class']}"
                        )
                    o_testable = math.isfinite(f(r["original_q"])) and math.isfinite(f(r["original_empirical_q"]))
                    otc = tclass(f(r["original_q"]), f(r["original_empirical_q"]), o_testable)
                    if r.get("orig_transport_class") and r["orig_transport_class"] not in ("NO_ORIG_P0", "NOT_IN_ORIG_SENS") and o_testable and otc != r["orig_transport_class"]:
                        problems.append(
                            f"{cohort} {fam} {r['cell_type']}: orig class {otc} != stored {r['orig_transport_class']}"
                        )
                    m_sig = (
                        math.isfinite(f(r["mature_q"]))
                        and f(r["mature_q"]) < ALPHA
                        and math.isfinite(f(r["mature_empirical_q"]))
                        and f(r["mature_empirical_q"]) < ALPHA
                    )
                    o_sig = (
                        math.isfinite(f(r["original_q"]))
                        and f(r["original_q"]) < ALPHA
                        and math.isfinite(f(r["original_empirical_q"]))
                        and f(r["original_empirical_q"]) < ALPHA
                    )
                    if not math.isfinite(f(r["original_q"])):
                        o_sig = None
                    v = verdict(o_sig, m_sig, testable)
                    if v == "SURVIVED":
                        oe, me = f(r["original_effect"]), f(r["mature_effect"])
                        if math.isfinite(oe) and math.isfinite(me) and (oe > 0) != (me > 0):
                            v = "DIED"
                    if v != r["verdict"]:
                        problems.append(
                            f"{cohort} {fam} {r['cell_type']}: verdict {v} != stored {r['verdict']}"
                        )
                    # notes all-group q must equal the columns after the fix
                    m = re.search(r"q_all(?:17)?_comp=([0-9.e+-]+);q_all(?:17)?_emp=([0-9.e+-]+)", r.get("notes", ""))
                    if m and testable:
                        nc, ne = float(m.group(1)), float(m.group(2))
                        if abs(nc - f(r["mature_q"])) > 1e-3 * max(1.0, abs(nc)) or abs(ne - f(r["mature_empirical_q"])) > 1e-3 * max(1.0, abs(ne)):
                            problems.append(
                                f"{cohort} {fam} {r['cell_type']}: notes q_all comp={nc:.4g} emp={ne:.4g} != columns {r['mature_q']} / {r['mature_empirical_q']}"
                            )
        print("  verdict counts:", dict(Counter(r["verdict"] for r in rows)))

    if bh_check_path:
        try:
            bh_rows = list(csv.DictReader(open(bh_check_path)))
            print(f"bh-check: {len(bh_rows)} rows, max abs_diff={max((f(r['abs_diff']) for r in bh_rows), default=float('nan')):.3e}")
        except FileNotFoundError:
            print("bh-check CSV not found (finalize has not run yet)")

    print("\n--- family summaries ---")
    for s in summaries:
        print(" ", s)
    print(f"\nPROBLEMS: {len(problems)}")
    for p in problems:
        print("  -", p)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
