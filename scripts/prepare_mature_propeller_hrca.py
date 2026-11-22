#!/usr/bin/env python3
"""Build mature HRCA propeller inputs (donor x cell-type counts + donor metadata).

Mature only: human donors 30 years and older, per the SUITES definitions in
scripts/run_followup_mature_organ_transports_b_hrca.py. Loads the HRCA
aggregate caches via _aggregate from scripts/run_followup_organ_10_eye_hrca.py
(user authorized loading sn_agg.npz and sc_agg.npz on beast). Layout mirrors
scripts/prepare_mature_propeller.py: one subdir per suite with counts.csv and
metadata.csv, plus manifest.csv at the root. Cell types with zero cells in a
suite are dropped before writing.
"""

import hashlib
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

OUTPUT = REPO / "data" / "interim" / "mature-propeller-hrca"


def write_cohort(root, species, organ, cohort, counts, metadata):
    key = hashlib.sha256(f"{species}|{organ}|{cohort}".encode()).hexdigest()[:16]
    target = root / key
    target.mkdir(parents=True, exist_ok=True)
    counts.to_csv(target / "counts.csv")
    metadata.loc[counts.index].to_csv(target / "metadata.csv")
    return {"key": key, "species": species, "organ": organ, "cohort": cohort,
            "n_donors": len(counts), "n_cell_types": counts.shape[1]}


def main():
    from scripts.run_followup_mature_organ_transports_b_hrca import SUITES, build_donor_lists
    from scripts.run_followup_organ_10_eye_hrca import CLASSES, SC_PATH, SN_PATH, _aggregate, _donor_meta

    OUTPUT.mkdir(parents=True, exist_ok=True)
    sn = _aggregate(SN_PATH, "sn")
    sc = _aggregate(SC_PATH, "sc")
    sck = sc["keys"]
    sck = sck[~((sck["study"] == "Sanes_GSE148077") & (sck["regiongroup"] != "macular"))].copy()
    sc = {**sc, "keys": sck}
    donor_lists = build_donor_lists(sn["keys"], sck)

    manifest = []
    for suite in SUITES:
        agg = sn if suite["assay"] == "sn" else sc
        keys = agg["keys"]
        donors = donor_lists[suite["donors_key"]]
        regions = list(suite["regions"])
        k = keys[keys["donor_id"].isin(donors) & keys["regiongroup"].isin(regions)].copy()
        meta = _donor_meta(k, donors, str(suite["studies_ref"]))
        assert float(meta["age"].min()) >= 30.0, f"{suite['name']}: under-30 donor leaked"
        cc = k.pivot_table(index="donor_id", columns="majorclass", values="n_cells",
                           aggfunc="sum", fill_value=0)
        cc = cc.reindex(index=meta.index, columns=[c for c in CLASSES if c in cc.columns],
                        fill_value=0)
        cc = cc.loc[:, cc.sum(axis=0) > 0]
        meta2 = meta[["age", "sex"]].copy()
        meta2["sex_male"] = (meta2["sex"] == "male").astype(float)
        meta2 = meta2[["age", "sex_male"]]
        assert set(cc.index) == set(meta2.index), f"{suite['name']}: donor mismatch"
        cohort = "HRCA_" + suite["name"]
        row = write_cohort(OUTPUT, "Human", "retina", cohort, cc.astype(int), meta2)
        row["min_donors_suite"] = suite["min_donors"]
        manifest.append(row)
        print(f"{cohort}: donors={len(cc)} types={cc.shape[1]} "
              f"age {meta['age'].min():.0f}-{meta['age'].max():.0f}", flush=True)

    pd.DataFrame(manifest).to_csv(OUTPUT / "manifest.csv", index=False)
    print(f"wrote {len(manifest)} suites to {OUTPUT}", flush=True)


if __name__ == "__main__":
    main()
