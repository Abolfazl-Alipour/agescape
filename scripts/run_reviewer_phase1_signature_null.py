#!/usr/bin/env python3

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import rdata
import scipy.sparse as sp
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.run_m10_pipeline import RDS, ols
from src.aging_programs import load_program_tiers

RESULTS = REPO / "results"


def random_signature_scores(expression, p1_up, p1_down, iterations, rng, available):
    genes = np.asarray(expression["genes"])
    matrix = expression["matrix"]
    means = np.asarray(matrix.mean(axis=1)).ravel()
    deciles = pd.qcut(means, 10, labels=False, duplicates="drop")
    deciles = pd.Series(deciles).fillna(-1).to_numpy()
    rows = []
    for iteration in range(iterations):
        selected_by_direction = []
        for program_genes in (p1_up, p1_down):
            selected = []
            for decile, count in pd.Series(deciles[program_genes]).value_counts().items():
                pool = np.flatnonzero((deciles == decile) & available)
                selected.extend(rng.choice(pool, size=int(count), replace=False))
            selected_by_direction.append(selected)
        up_score = np.asarray(matrix[selected_by_direction[0], :].mean(axis=0)).ravel()
        down_score = np.asarray(matrix[selected_by_direction[1], :].mean(axis=0)).ravel()
        rows.append(up_score - down_score)
    return np.vstack(rows)


def main():
    iterations = 200
    rng = np.random.default_rng(20261105)
    obj = rdata.read_rds(str(RDS))
    cd = obj.colData
    donor_ids = np.array([str(donor) for donor in cd.rownames])
    metadata = {name: np.asarray(cd.listData[name]).ravel()
                for name in ("Age", "Sex", "Brain_bank", "PMI", "TechPC1", "TechPC2", "TechPC3")}
    program = load_program_tiers(REPO, "human")["P1_Clean"]
    actual = pd.read_csv(RESULTS / "m10-celltype-regressions.csv").set_index("cell_type")
    rows = []
    for cell_type in sorted(obj.assays.data.listData.keys()):
        source = obj.assays.data.listData[cell_type]
        genes = np.array([str(gene) for gene in source.Dimnames[0]])
        mat = sp.csc_matrix((np.asarray(source.x), np.asarray(source.i), np.asarray(source.p)),
                            shape=tuple(source.Dim))
        library = np.asarray(mat.sum(axis=0)).ravel()
        keep = library > 0
        mat = mat[:, keep]
        normalized = mat @ sp.diags(1e6 / library[keep])
        normalized = normalized.tocsc()
        normalized.data = np.log1p(normalized.data)
        lookup = {gene: index for index, gene in enumerate(genes)}
        p1_up = np.array([lookup[gene] for gene in program["up"] if gene in lookup], dtype=int)
        p1_down = np.array([lookup[gene] for gene in program["down"] if gene in lookup], dtype=int)
        available = np.ones(len(genes), dtype=bool)
        available[np.concatenate([p1_up, p1_down])] = False
        original_matrix = normalized
        random_scores = random_signature_scores(
            {"matrix": original_matrix, "genes": genes}, p1_up, p1_down, iterations, rng, available
        )
        age = metadata["Age"][keep].astype(float)
        male = (metadata["Sex"][keep] == "Male").astype(float)
        hbcc = (metadata["Brain_bank"][keep] == "HBCC").astype(float)
        covariates = np.column_stack([
            np.ones(len(age)), age, male, hbcc, metadata["PMI"][keep].astype(float),
            metadata["TechPC1"][keep].astype(float), metadata["TechPC2"][keep].astype(float),
            metadata["TechPC3"][keep].astype(float),
        ])
        random_fits = [ols(score, covariates) for score in random_scores]
        random_slopes = np.array([fit["slope"] for fit in random_fits])
        observed = actual.loc[cell_type, "slope_covar_adj"]
        depth_correlations = [pearsonr(score, np.log(library[keep])).statistic for score in random_scores]
        rows.append({
            "cell_type": cell_type, "n_donors": len(age), "n_random_signatures": iterations,
            "actual_slope": observed, "random_mean_slope": random_slopes.mean(),
            "random_sd": random_slopes.std(ddof=1),
            "random_ci95_low": np.quantile(random_slopes, 0.025),
            "random_ci95_high": np.quantile(random_slopes, 0.975),
            "empirical_p_two_sided": (1 + np.sum(np.abs(random_slopes - random_slopes.mean())
                                                   >= abs(observed - random_slopes.mean()))) / (iterations + 1),
            "p1_up_detected": len(p1_up), "p1_down_detected": len(p1_down),
            "p1_coverage": (len(p1_up) + len(p1_down)) / (len(program["up"]) + len(program["down"])),
            "mean_log_library": np.log(library[keep]).mean(),
            "actual_score_depth_r": pearsonr(
                np.asarray(original_matrix[p1_up, :].mean(axis=0)).ravel()
                - np.asarray(original_matrix[p1_down, :].mean(axis=0)).ravel(),
                np.log(library[keep])
            ).statistic,
            "random_score_depth_r_mean": np.mean(depth_correlations),
        })
        print(f"{cell_type}: actual {observed:.6g}, random mean {random_slopes.mean():.6g}", flush=True)
    result = pd.DataFrame(rows)
    result["q_bh"] = multipletests(result["empirical_p_two_sided"], method="fdr_bh")[1]
    result.to_csv(RESULTS / "reviewer-phase1-signature-random-null.csv", index=False)


if __name__ == "__main__":
    main()
