from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent


def test_phase3_figures_have_data_and_captions():
    names = ["mouse-vector-drift", "interneuron-gene-changes", "composition-changes", "universality-concordance"]
    data_files = ["mouse-drift", "interneuron", "composition", "universality"]
    for name, data_name in zip(names, data_files):
        for extension in ("png", "pdf"):
            path = REPO / "figures" / f"reviewer-phase3-{name}.{extension}"
            assert path.exists() and path.stat().st_size > 1000
        data = pd.read_csv(REPO / "results" / f"reviewer-phase3-{data_name}-figure-data.csv")
        assert not data.empty
    captions = (REPO / "reports" / "reviewer_phase3_figure_captions.md").read_text()
    assert captions.count("## ") == 4


def test_reference_list_excludes_unverified_matches():
    metadata = pd.read_csv(REPO / "results" / "reviewer-phase3-reference-metadata.csv")
    verified = metadata[metadata["verification_status"] == "verified"]
    assert len(verified) == 32
    assert verified["title_similarity"].min() >= 0.85
    assert verified["doi"].notna().all()
    report = (REPO / "reports" / "reviewer_phase3_references.md").read_text()
    assert "Unresolved manuscript citation" in report
    assert "not verified" in report


def test_phase3_does_not_modify_manuscript():
    import subprocess
    diff = subprocess.run(["git", "diff", "868ca108d8e14ed89ceca73c77b7086dfd54faf9", "--",
                           "report/manuscript-aging-mechanisms.md"], capture_output=True, text=True, check=True)
    assert not diff.stdout
