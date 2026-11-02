import numpy as np
import pandas as pd

from scripts.extended_atlas_utils import (
    afca_sample_id, clr_table, composition_tests, zebrafish_sample_metadata,
)


def test_afca_sample_id_preserves_full_library_suffix():
    names = pd.Index([
        "AAACCCACAGTGAGCA-1_AFCA_female_head_30_S1",
        "AAACCCACAGTGAGCA-1_FCA10_Male_head_adult",
    ])
    assert afca_sample_id(names).tolist() == [
        "AFCA_female_head_30_S1", "FCA10_Male_head_adult",
    ]


def test_zebrafish_metadata_distinguishes_fish_from_technical_replicate():
    parsed = zebrafish_sample_metadata("GSM8475554_mo5_fish1_rep1_barcodes.tsv.gz")
    assert parsed == {"sample": "mo5_fish1_rep1", "donor": "mo5_fish1", "age": 5, "replicate": 1}


def test_clr_rows_sum_to_zero():
    counts = pd.DataFrame({"a": [10, 20], "b": [5, 0], "c": [0, 10]})
    assert np.allclose(clr_table(counts).sum(axis=1), 0)


def test_composition_uses_samples_not_cells_as_replication():
    counts = pd.DataFrame(
        {"a": [90, 80, 70, 60, 50, 40], "b": [10, 20, 30, 40, 50, 60]},
        index=[f"s{i}" for i in range(6)],
    )
    metadata = pd.DataFrame({"age": [1, 2, 3, 4, 5, 6], "sex_male": 0}, index=counts.index)
    rows = composition_tests(counts, metadata, "Test", "Organ", "Cohort", "day")
    assert len(rows) == 2
    assert all(row["n_samples"] == 6 for row in rows)


def test_bayesian_composition_input_preserves_cell_type_column():
    from scripts.run_extended_bayesian_composition import model_input

    counts = pd.DataFrame({"a": [90, 80], "b": [10, 60]}, index=pd.Index(["s1", "s2"], name="sample"))
    metadata = pd.DataFrame({"age": [1, 2], "sex_male": 0}, index=counts.index)
    result = model_input(counts, metadata, "Test", "Organ", "Cohort")
    assert list(result.columns[:2]) == ["donor", "cell_type"]
    assert set(result["cell_type"]) == {"a", "b"}
    assert result["donor"].tolist() == ["s1", "s1", "s2", "s2"]
