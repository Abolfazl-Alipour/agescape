from __future__ import annotations

import re

import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests


UNKNOWN_LABELS = {"", "unknown", "unassigned", "nan", "none", "na"}


def known_cell_type(values: pd.Series) -> np.ndarray:
    return ~values.astype(str).str.strip().str.lower().isin(UNKNOWN_LABELS).to_numpy()


def common_cell_type(label: str) -> str:
    key = re.sub(r"[^a-z0-9]+", "", str(label).lower())
    targeted = {
        "cd8naive": "CD8 naive T", "cd8naivet": "CD8 naive T",
        "cd8naivetcell": "CD8 naive T", "cd8naivetcells": "CD8 naive T",
        "mait": "MAIT", "maitcell": "MAIT", "maitcells": "MAIT",
        "nonclassicalmonocyte": "non-classical monocytes",
        "nonclassicalmonocytes": "non-classical monocytes",
        "cd16monocyte": "non-classical monocytes", "cd16monocytes": "non-classical monocytes",
        "cd16mc": "non-classical monocytes", "cd16mcs": "non-classical monocytes",
    }
    broad = {
        "cd4t": "CD4 T", "cd4tcell": "CD4 T", "cd4tcells": "CD4 T",
        "cd8t": "CD8 T", "cd8tcell": "CD8 T", "cd8tcells": "CD8 T",
        "b": "B", "bcell": "B", "bcells": "B",
        "nk": "NK", "nkcell": "NK", "nkcells": "NK",
        "mono": "monocytes", "monocyte": "monocytes", "monocytes": "monocytes",
        "dc": "dendritic cells", "dendriticcell": "dendritic cells",
        "dendriticcells": "dendritic cells",
    }
    return targeted.get(key, broad.get(key, ""))


def clr_transform(counts: pd.DataFrame, pseudocount: float = 0.5) -> pd.DataFrame:
    if pseudocount <= 0:
        raise ValueError("CLR pseudocount must be positive")
    adjusted = counts.astype(float) + pseudocount
    log_props = np.log(adjusted.div(adjusted.sum(axis=1), axis=0))
    result = log_props.sub(log_props.mean(axis=1), axis=0)
    if not np.isfinite(result.to_numpy()).all() or not np.allclose(result.sum(axis=1), 0, atol=1e-10):
        raise ValueError("Invalid CLR transform")
    return result


def technical_gene_mask(symbols: np.ndarray) -> np.ndarray:
    genes = pd.Index(symbols).astype(str).str.upper()
    return np.array([
        bool(re.match(r"^(MT-|RP[LS]|MRP[LS])", gene))
        or gene in {"XIST", "TSIX", "DDX3Y", "EIF2S3Y", "UTY", "KDM5D"}
        for gene in genes
    ])


def design_matrix(frame: pd.DataFrame, design: str, variant: str, endpoint: str) -> tuple[np.ndarray, list[str]]:
    n = len(frame)
    term = "age" if design == "exact_age" else "old"
    columns = [np.ones(n), frame[term].to_numpy(float)]
    names = ["intercept", term]
    male = (frame["sex"].astype(str).str.lower() == "male").to_numpy(float)
    if len(np.unique(male)) == 2:
        columns.append(male)
        names.append("male")
    if variant == "primary" and endpoint == "Hellinger":
        variant = "depth"
    if variant == "depth":
        depth = "qc_mean_log_umi" if endpoint == "Hellinger" else "mean_log_umi"
        columns.append(frame[depth].to_numpy(float))
        names.append(depth)
    elif variant == "log_total_cells":
        columns.append(np.log(frame["total_cells"].to_numpy(float)))
        names.append("log_total_cells")
    elif variant in {"cmv", "cmv_depth"}:
        columns.append((frame["cmv"].astype(str).str.lower() == "pos").to_numpy(float))
        names.append("cmv_pos")
        if variant == "cmv_depth":
            columns.append(frame["qc_mean_log_umi"].to_numpy(float))
            names.append("qc_mean_log_umi")
    elif variant in {"batch", "batch_depth", "pool", "pool_depth"}:
        kind = "batch" if variant.startswith("batch") else "pool"
        dummies = pd.get_dummies(frame[kind].astype(str), prefix=kind, drop_first=True, dtype=float)
        if dummies.shape[1]:
            columns.extend(dummies[column].to_numpy() for column in dummies)
            names.extend(dummies.columns)
        if variant.endswith("_depth"):
            columns.append(frame["qc_mean_log_umi"].to_numpy(float))
            names.append("qc_mean_log_umi")
    elif variant not in {"primary", "unadjusted", "log_total_cells"}:
        raise ValueError(f"Unknown model variant: {variant}")
    X = np.column_stack(columns)
    return X, names


def fit_donor_model(frame: pd.DataFrame, outcome: str, design: str, variant: str,
                    endpoint: str) -> dict:
    frame = frame.dropna(subset=[outcome]).copy()
    term = "age" if design == "exact_age" else "old"
    X, names = design_matrix(frame, design, variant, endpoint)
    n, k = X.shape
    base = {
        "n_donors": n, "n_cells": int(frame["n_cells"].sum()),
        "n_young": int((frame[term] == 0).sum()) if term == "old" else np.nan,
        "n_old": int((frame[term] == 1).sum()) if term == "old" else np.nan,
        "covariates": ";".join(names), "effect_term": "old_vs_young" if term == "old" else "age_per_year",
    }
    if n < 6 or n <= k + 1 or np.linalg.matrix_rank(X) < k or not np.isfinite(X).all():
        return {**base, "status": "untestable", "reason": "rank_or_residual_degrees_of_freedom"}
    if term == "old":
        groups = frame[term].to_numpy()
        if min((groups == 0).sum(), (groups == 1).sum()) < 3:
            return {**base, "status": "untestable", "reason": "fewer_than_three_donors_per_group"}
    elif frame["age"].nunique() < 2:
        return {**base, "status": "untestable", "reason": "age_does_not_vary"}
    fit = sm.OLS(frame[outcome].to_numpy(float), X).fit()
    ci = fit.conf_int()[1]
    return {
        **base, "status": "fit", "reason": "", "effect": float(fit.params[1]),
        "se": float(fit.bse[1]), "ci95_low": float(ci[0]), "ci95_high": float(ci[1]),
        "pvalue": float(fit.pvalues[1]), "r2": float(fit.rsquared),
    }


def fit_vectorized(outcomes: np.ndarray, frame: pd.DataFrame, design: str, variant: str,
                   endpoint: str = "pseudobulk") -> dict[str, np.ndarray] | None:
    """Fit gene x donor outcomes with the age/group effect in design column two."""
    X, _ = design_matrix(frame, design, variant, endpoint)
    n, k = X.shape
    if n < 6 or n <= k + 1 or np.linalg.matrix_rank(X) < k or not np.isfinite(X).all():
        return None
    term = "age" if design == "exact_age" else "old"
    if term == "old":
        groups = frame[term].to_numpy()
        if min((groups == 0).sum(), (groups == 1).sum()) < 3:
            return None
    elif frame["age"].nunique() < 2:
        return None
    y = np.asarray(outcomes, dtype=float)
    coefficients, _, _, _ = np.linalg.lstsq(X, y.T, rcond=None)
    fitted = X @ coefficients
    residuals = y.T - fitted
    dof = n - k
    sigma2 = np.sum(residuals ** 2, axis=0) / dof
    inverse = np.linalg.inv(X.T @ X)
    se = np.sqrt(np.maximum(sigma2 * inverse[1, 1], 0))
    effect = coefficients[1]
    with np.errstate(divide="ignore", invalid="ignore"):
        tstat = effect / se
    pvalue = 2 * stats.t.sf(np.abs(tstat), dof)
    critical = stats.t.ppf(0.975, dof)
    return {
        "effect": effect, "se": se, "ci95_low": effect - critical * se,
        "ci95_high": effect + critical * se, "pvalue": pvalue,
    }


def add_bh(frame: pd.DataFrame, groups: list[str], p_column: str = "pvalue",
           q_column: str = "q_bh") -> pd.DataFrame:
    frame = frame.copy()
    frame[q_column] = np.nan
    if frame.empty:
        return frame
    for _, indices in frame.groupby(groups, dropna=False).groups.items():
        valid = frame.loc[indices, p_column].dropna().index
        if len(valid):
            frame.loc[valid, q_column] = multipletests(frame.loc[valid, p_column], method="fdr_bh")[1]
    return frame


def random_signature_scores(matrix: np.ndarray, up: np.ndarray, down: np.ndarray,
                            iterations: int, rng: np.random.Generator,
                            exclude: np.ndarray | None = None) -> np.ndarray:
    means = np.asarray(matrix.mean(axis=1)).ravel()
    deciles = pd.Series(pd.qcut(means, 10, labels=False, duplicates="drop")).fillna(-1).to_numpy()
    available = np.ones(matrix.shape[0], dtype=bool)
    available[np.concatenate([up, down])] = False
    if exclude is not None:
        available[np.asarray(exclude, dtype=int)] = False
    scores = []
    for _ in range(iterations):
        selected_directions = []
        for program in (up, down):
            selected = []
            for decile, count in pd.Series(deciles[program]).value_counts().items():
                pool = np.flatnonzero((deciles == decile) & available)
                if len(pool) < int(count):
                    raise ValueError("Insufficient expression-matched random-signature pool")
                selected.extend(rng.choice(pool, size=int(count), replace=False))
            selected_directions.append(selected)
        scores.append(matrix[selected_directions[0]].mean(axis=0) - matrix[selected_directions[1]].mean(axis=0))
    return np.vstack(scores)


def fast_effect(score: np.ndarray, frame: pd.DataFrame, design: str = "strict_group") -> float:
    X, _ = design_matrix(frame, design, "primary", "P1_transport")
    return float(np.linalg.lstsq(X, np.asarray(score, float), rcond=None)[0][1])


def hellinger_distances(matrix, reference_props: np.ndarray, nontechnical_totals: np.ndarray) -> np.ndarray:
    root = matrix.copy().astype(np.float64)
    root.data = np.sqrt(root.data)
    affinity = np.asarray(root @ np.sqrt(reference_props)).ravel() / np.sqrt(nontechnical_totals)
    return np.sqrt(np.maximum(0.0, 1.0 - affinity))
