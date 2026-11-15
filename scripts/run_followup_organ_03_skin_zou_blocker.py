"""Dataset 03: Zou et al. (2021) eyelid skin scRNA-seq (HRA000395) access-check report."""

from __future__ import annotations

from pathlib import Path
import sys

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.followup_organs import write_clean_csv, write_clean_markdown

ZOU_DONORS = [
    {"donor_label": "Y1", "age_group": "Young", "age_years": 18, "sex": "F", "tissue": "Upper eyelid skin", "accession": "HRA000395", "bioproject": "PRJCA003032", "access_status": "Controlled_DAC_SecurityCheck"},
    {"donor_label": "Y2", "age_group": "Young", "age_years": 22, "sex": "F", "tissue": "Upper eyelid skin", "accession": "HRA000395", "bioproject": "PRJCA003032", "access_status": "Controlled_DAC_SecurityCheck"},
    {"donor_label": "Y3", "age_group": "Young", "age_years": 23, "sex": "F", "tissue": "Upper eyelid skin", "accession": "HRA000395", "bioproject": "PRJCA003032", "access_status": "Controlled_DAC_SecurityCheck"},
    {"donor_label": "M1", "age_group": "Middle", "age_years": 44, "sex": "F", "tissue": "Upper eyelid skin", "accession": "HRA000395", "bioproject": "PRJCA003032", "access_status": "Controlled_DAC_SecurityCheck"},
    {"donor_label": "M2", "age_group": "Middle", "age_years": 47, "sex": "F", "tissue": "Upper eyelid skin", "accession": "HRA000395", "bioproject": "PRJCA003032", "access_status": "Controlled_DAC_SecurityCheck"},
    {"donor_label": "M3", "age_group": "Middle", "age_years": 48, "sex": "F", "tissue": "Upper eyelid skin", "accession": "HRA000395", "bioproject": "PRJCA003032", "access_status": "Controlled_DAC_SecurityCheck"},
    {"donor_label": "O1", "age_group": "Old", "age_years": 70, "sex": "F", "tissue": "Upper eyelid skin", "accession": "HRA000395", "bioproject": "PRJCA003032", "access_status": "Controlled_DAC_SecurityCheck"},
    {"donor_label": "O2", "age_group": "Old", "age_years": 73, "sex": "F", "tissue": "Upper eyelid skin", "accession": "HRA000395", "bioproject": "PRJCA003032", "access_status": "Controlled_DAC_SecurityCheck"},
    {"donor_label": "O3", "age_group": "Old", "age_years": 76, "sex": "F", "tissue": "Upper eyelid skin", "accession": "HRA000395", "bioproject": "PRJCA003032", "access_status": "Controlled_DAC_SecurityCheck"},
]


def main() -> None:
    results_dir = REPO_ROOT / "results"
    reports_dir = REPO_ROOT / "reports"
    prefix = "followup-organ-03-skin-zou-hra000395"

    df = pd.DataFrame(ZOU_DONORS)
    write_clean_csv(df, results_dir / f"{prefix}-metadata.csv")

    donor_rows = "\n".join(
        f"| {r['donor_label']} | {r['age_group']} | {r['age_years']} | {r['sex']} | {r['tissue']} |"
        for _, r in df.iterrows()
    )

    md_text = f"""# Dataset 03: Human Eyelid Skin Aging Access Check (Zou et al., HRA000395)

## 1. Access-Check Disposition
- **Status**: **BLOCKED (Controlled Access, Access Check Only; Do Not Apply)**.
- **Primary Accession**: `HRA000395` (`PRJCA003032`) in NGDC/CNCB GSA-Human (`https://ngdc.cncb.ac.cn/gsa-human/browse/HRA000395`).
- **Publication**: Zou et al., *Developmental Cell* 2021;56(3):383-397.e8 (`PMID: 33238152`, `DOI: 10.1016/j.devcel.2020.11.002`).
- **Observed Access Terms**: Anonymous HTTPS inspection of `https://ngdc.cncb.ac.cn/gsa-human/browse/HRA000395` returned an automated `SecurityCheck` gate (`Content-Type: text/html`, `<title>SecurityCheck</title>`, `Content-Length: 1014`). GSA-Human individual human sequencing datasets require registered investigator authentication, institutional email verification, and formal Data Access Committee (DAC) approval prior to data release. Per project rules, no DAC application was submitted and no gate bypass was attempted.

## 2. Verification of Donor Ages and Sex
- **Do exact ages exist?** **Yes**. Exact chronological ages are documented in the primary publication (`n = 9` healthy female upper-eyelid donors; ages `18, 22, 23` in Young, `44, 47, 48` in Middle, and `70, 73, 76` in Old; note 100% female cohort with upper eyelid sun exposure differences).

| Donor Slot | Age Group | Age (yr) | Sex | Tissue |
|---|---|---:|---|---|
{donor_rows}

## 3. Analyses Run
- **None run** (no expression or cell-count matrices downloaded). Documented blocker recorded cleanly.

## 4. Output Artifacts
- [Verified Metadata CSV](file:///home/reza/projects/agescape-revision/results/{prefix}-metadata.csv)
"""
    write_clean_markdown(md_text, reports_dir / f"{prefix}.md")
    print("DONE Dataset 03 (BLOCKED access check documented)")


if __name__ == "__main__":
    main()
