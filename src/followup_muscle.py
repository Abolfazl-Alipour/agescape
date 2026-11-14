from __future__ import annotations

import re

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
from statsmodels.stats.multitest import multipletests

MIN_DONORS = 6
MIN_PER_GROUP = 2
MIN_CELLS = 20
MIN_REFERENCE_CELLS = 50
MIN_UMI_HELLINGER = 1500
MIN_RESIDUAL_DF = 2
NULL_ITERATIONS = 200
NULL_SEED = 20261105
UNKNOWN_LABELS = {"", "unknown", "unassigned", "nan", "none"}


def known_cell_type(values: pd.Series) -> np.ndarray:
    return ~values.astype(str).str.strip().str.lower().isin(UNKNOWN_LABELS).to_numpy()


def clr_transform(counts: pd.DataFrame, pseudocount: float = 0.5) -> pd.DataFrame:
    if pseudocount <= 0:
        raise ValueError("CLR pseudocount must be positive")
    adjusted = counts.astype(float) + pseudocount
    log_props = np.log(adjusted.div(adjusted.sum(axis=1), axis=0))
    result = log_props.sub(log_props.mean(axis=1), axis=0)
    if not np.isfinite(result.to_numpy()).all():
        raise ValueError("Non-finite CLR value")
    if not np.allclose(result.sum(axis=1), 0.0, atol=1e-10):
        raise AssertionError("CLR rows must sum to zero")
    return result


def technical_gene_mask(gene_names: np.ndarray) -> np.ndarray:
    genes = pd.Index(gene_names).astype(str).str.upper()
    return np.array([
        bool(re.match(r"^(MT-|RP[LS]|MRP[LS])", gene))
        or gene in {"XIST", "TSIX", "DDX3Y", "EIF2S3Y", "UTY", "KDM5D"}
        for gene in genes
    ])


def count_gate(values: np.ndarray) -> dict[str, int]:
    values = np.asarray(values)
    return {
        "nonfinite": int((~np.isfinite(values)).sum()),
        "negative": int((values < 0).sum()),
        "noninteger": int((values != np.floor(values)).sum()),
    }


def assert_count_gate(values: np.ndarray) -> None:
    failures = count_gate(values)
    if any(failures.values()):
        raise ValueError(f"Raw-count gate failed: {failures}")


def build_design(frame: pd.DataFrame, age: str = "old", depth: bool = True,
                 covariates: tuple[str, ...] = ("male", "dcd")) -> tuple[np.ndarray, list[str], list[str]]:
    if age not in frame or frame[age].nunique() != 2:
        raise ValueError("The age contrast must have exactly two observed groups")
    if (frame[age].value_counts() < MIN_PER_GROUP).any():
        raise ValueError("Each age group needs at least two donors")
    columns = [np.ones(len(frame)), frame[age].to_numpy(float)]
    names = ["intercept", age]
    omitted = []
    for name in covariates:
        if name not in frame or frame[name].nunique() < 2:
            omitted.append(name)
            continue
        candidate = np.column_stack([*columns, frame[name].to_numpy(float)])
        if np.linalg.matrix_rank(candidate) == candidate.shape[1]:
            columns = [candidate[:, i] for i in range(candidate.shape[1])]
            names.append(name)
        else:
            omitted.append(name)
    if depth:
        if "mean_log_umi" not in frame or not np.isfinite(frame["mean_log_umi"]).all():
            raise ValueError("Depth adjustment requires finite mean_log_umi")
        candidate = np.column_stack([*columns, frame["mean_log_umi"].to_numpy(float)])
        if np.linalg.matrix_rank(candidate) != candidate.shape[1]:
            raise ValueError("Depth adjustment is aliased with another design column")
        columns = [candidate[:, i] for i in range(candidate.shape[1])]
        names.append("mean_log_umi")
    design = np.column_stack(columns)
    if not np.isfinite(design).all() or np.linalg.matrix_rank(design) < design.shape[1]:
        raise ValueError("Model design is not finite and full rank")
    if len(frame) < MIN_DONORS or len(frame) - design.shape[1] < MIN_RESIDUAL_DF:
        raise ValueError("Insufficient donors or residual degrees of freedom")
    return design, names, omitted


def fit_ols(outcomes: np.ndarray, design: np.ndarray, effect_column: int = 1) -> dict[str, np.ndarray]:
    outcomes = np.asarray(outcomes, dtype=np.float64)
    if outcomes.ndim == 1:
        outcomes = outcomes[:, None]
    inverse = np.linalg.inv(design.T @ design)
    coefficients = inverse @ design.T @ outcomes
    residual = outcomes - design @ coefficients
    residual_df = len(design) - design.shape[1]
    residual_ss = np.sum(residual ** 2, axis=0)
    variance = residual_ss / residual_df
    effect = coefficients[effect_column]
    se = np.sqrt(variance * inverse[effect_column, effect_column])
    with np.errstate(divide="ignore", invalid="ignore"):
        statistic = effect / se
    pvalue = 2 * stats.t.sf(np.abs(statistic), residual_df)
    critical = stats.t.ppf(0.975, residual_df)
    total_ss = np.sum((outcomes - outcomes.mean(axis=0)) ** 2, axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        r2 = 1 - residual_ss / total_ss
    return {
        "effect": effect,
        "se": se,
        "ci95_low": effect - critical * se,
        "ci95_high": effect + critical * se,
        "pvalue": pvalue,
        "residual_df": np.full(outcomes.shape[1], residual_df),
        "r2": r2,
    }


def add_bh(frame: pd.DataFrame, groups: list[str], p_column: str = "pvalue",
           q_column: str = "q_bh") -> pd.DataFrame:
    result = frame.copy()
    result[q_column] = np.nan
    if result.empty:
        return result
    for _, indices in result.groupby(groups, dropna=False).groups.items():
        valid = result.loc[indices, p_column].dropna().index
        if len(valid):
            result.loc[valid, q_column] = multipletests(
                result.loc[valid, p_column], method="fdr_bh"
            )[1]
    return result


def expression_deciles(matrix: np.ndarray) -> np.ndarray:
    means = np.asarray(matrix.mean(axis=1)).ravel()
    return pd.Series(pd.qcut(means, 10, labels=False, duplicates="drop")).fillna(-1).to_numpy()


def random_signature_scores(matrix: np.ndarray, up: np.ndarray, down: np.ndarray,
                            iterations: int = NULL_ITERATIONS,
                            seed: int = NULL_SEED) -> np.ndarray:
    rng = np.random.default_rng(seed)
    deciles = expression_deciles(matrix)
    available = np.ones(matrix.shape[0], dtype=bool)
    available[np.concatenate([up, down])] = False
    scores = []
    for _ in range(iterations):
        directions = []
        for program_genes in (up, down):
            selected = []
            for decile, count in pd.Series(deciles[program_genes]).value_counts().items():
                pool = np.flatnonzero((deciles == decile) & available)
                if len(pool) < int(count):
                    raise ValueError("Insufficient expression-matched random-signature pool")
                selected.extend(rng.choice(pool, size=int(count), replace=False))
            directions.append(np.asarray(selected, dtype=int))
        scores.append(matrix[directions[0], :].mean(axis=0) - matrix[directions[1], :].mean(axis=0))
    return np.vstack(scores)


def matched_null_effect_pvalue(observed: float, null_effects: np.ndarray) -> float:
    centered = null_effects - null_effects.mean()
    return float((1 + np.sum(np.abs(centered) >= abs(observed - null_effects.mean()))) /
                 (len(null_effects) + 1))


def hellinger_distances(matrix: sp.csr_matrix, reference_props: np.ndarray,
                        nontechnical_totals: np.ndarray) -> np.ndarray:
    root = matrix.copy().astype(np.float64)
    root.data = np.sqrt(root.data)
    affinity = np.asarray(root @ np.sqrt(reference_props)).ravel() / np.sqrt(nontechnical_totals)
    return np.sqrt(np.maximum(0.0, 1.0 - affinity))


def aggregate_block(matrix: sp.csr_matrix, codes: np.ndarray, output: np.ndarray) -> None:
    unique = np.unique(codes)
    remap = np.searchsorted(unique, codes)
    indicator = sp.csr_matrix((
        np.ones(len(codes), dtype=np.float64), (np.arange(len(codes)), remap)
    ), shape=(len(codes), len(unique)))
    block = (indicator.T @ matrix).toarray()
    if np.issubdtype(output.dtype, np.integer):
        if not np.equal(block, np.floor(block)).all():
            raise ValueError("Non-integer aggregate assigned to integer count cache")
        block = block.astype(output.dtype)
    output[unique] += block
