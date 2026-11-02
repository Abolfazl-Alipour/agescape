#!/usr/bin/env python3

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.extended_atlas_utils import clr_table, composition_counts

RESULTS = REPO / "results"


def model_input(counts, metadata, species, organ, cohort):
    counts = counts.loc[:, counts.sum(axis=0) >= 50]
    metadata = metadata.copy()
    metadata["mean_log_umi"] = np.log(counts.sum(axis=1))
    long = clr_table(counts).stack().rename("clr").reset_index()
    long = long.rename(columns={long.columns[0]: "donor", long.columns[1]: "cell_type"})
    long = long.merge(metadata, left_on="donor", right_index=True, validate="many_to_one")
    long["species"] = species
    long["tissue"] = organ
    long["cohort"] = cohort
    long["profile_hellinger"] = long["clr"]
    return long


def afca_inputs():
    import anndata as ad
    from scripts.run_afca_aging import afca_obs

    outputs = []
    for organ in ("head", "body"):
        atlas = ad.read_h5ad(REPO / "data" / "raw" / "afca" / f"afca-{organ}.h5ad", backed="r")
        obs = afca_obs(atlas)
        obs = obs[(obs["dataset"] == "AFCA") & (obs["cell_type"] != "unannotated")]
        counts = composition_counts(obs, "sample", "cell_type")
        metadata = obs.groupby("sample").agg(age=("age", "first"), sex_male=("sex_male", "first"))
        current = model_input(counts, metadata, "Fruit fly", organ, "AFCA")
        current["donor"] = (current["donor"].str.replace("_head_", "_", regex=False)
                            .str.replace("_body_", "_", regex=False))
        outputs.append(current)
    return pd.concat(outputs, ignore_index=True)


def zebrafish_inputs():
    import anndata as ad
    from scripts.run_zebrafish_testis_aging import standard_obs

    atlases = {
        "5-22mo": ad.read_h5ad(REPO / "data" / "interim" / "zebrafish-testis-5-22mo.h5ad", backed="r"),
        "27mo": ad.read_h5ad(REPO / "data" / "interim" / "zebrafish-testis-27mo.h5ad", backed="r"),
    }
    obs = pd.concat([standard_obs(a, name) for name, a in atlases.items()])
    obs = obs[obs["broad_cell_type"] != "unknown"]
    counts = composition_counts(obs, "donor", "broad_cell_type")
    metadata = obs.groupby("donor").agg(age=("age", "first"), sex_male=("sex_male", "first"))
    return model_input(counts, metadata, "Zebrafish", "Testis broad compartments", "GSE275361")


def worm_inputs():
    import anndata as ad
    from scripts.run_gao_worm_aging import normal_obs

    atlas = ad.read_h5ad(REPO / "data" / "interim" / "worm-gao-integrated.h5ad", backed="r")
    obs = normal_obs(atlas)
    counts = composition_counts(obs, "donor", "cell_type")
    metadata = obs.groupby("donor").agg(
        age=("age", "first"), sex_male=("sex_male", "first"),
        platform_nextseq=("platform_nextseq", "first"),
    )
    return model_input(counts, metadata, "C. elegans", "Whole animal", "Gao_GSE229022")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("species", choices=["Fruit fly", "Zebrafish", "C. elegans"])
    parser.add_argument("--skip-diagnostics", action="store_true")
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--use-existing", action="store_true")
    args = parser.parse_args()
    prefixes = {"Fruit fly": "afca", "Zebrafish": "zebrafish", "C. elegans": "worm"}
    input_functions = {"Fruit fly": afca_inputs, "Zebrafish": zebrafish_inputs, "C. elegans": worm_inputs}
    prefix = prefixes[args.species]
    input_path = RESULTS / f"extended-{prefix}-composition-model-input.csv"
    if args.use_existing:
        data = pd.read_csv(input_path)
    else:
        data = input_functions[args.species]()
        data.to_csv(input_path, index=False)
    if args.prepare_only:
        return
    import scripts.run_cross_species_bayesian_vector as bayes

    bayes.ARTIFACTS = RESULTS / "extended-composition-bayesian"
    bayes.ARTIFACTS.mkdir(parents=True, exist_ok=True)
    rows = []
    diagnostics = []
    for depth_adjusted in (False, True):
        posterior, _, diagnostic = bayes.fit_model(
            data, args.species, depth_adjusted, 0.5,
            diagnostics=not args.skip_diagnostics, draws=2000,
            fixed_group_effects=args.species == "Fruit fly",
            covariate_column="platform_nextseq" if args.species == "C. elegans" else "sex_male",
        )
        rows.extend(posterior)
        diagnostics.append(diagnostic)
    pd.DataFrame(rows).to_csv(RESULTS / f"extended-{prefix}-composition-bayesian-posterior.csv", index=False)
    pd.DataFrame(diagnostics).to_csv(
        RESULTS / f"extended-{prefix}-composition-bayesian-diagnostics.csv", index=False
    )


if __name__ == "__main__":
    main()
