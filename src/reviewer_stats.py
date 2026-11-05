import numpy as np
from scipy import stats


def vectorized_ols(expression, design, age_column=1):
    expression = np.asarray(expression, dtype=float)
    design = np.asarray(design, dtype=float)
    coefficients = np.linalg.pinv(design) @ expression
    residuals = expression - design @ coefficients
    df = design.shape[0] - np.linalg.matrix_rank(design)
    if df <= 0:
        raise ValueError("Insufficient residual degrees of freedom")
    residual_variance = np.sum(residuals ** 2, axis=0) / df
    age_variance = np.linalg.pinv(design.T @ design)[age_column, age_column]
    se = np.sqrt(residual_variance * age_variance)
    slopes = coefficients[age_column]
    statistics = np.divide(slopes, se, out=np.zeros_like(slopes), where=se > 0)
    pvalues = 2 * stats.t.sf(np.abs(statistics), df)
    return slopes, se, pvalues


def pairwise_pearson(effects):
    effects = np.asarray(effects, dtype=float)
    centered = effects - effects.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(centered, axis=1)
    valid = norms > 0
    centered = centered[valid]
    norms = norms[valid]
    correlation = (centered @ centered.T) / norms[:, None] / norms[None, :]
    return correlation[np.triu_indices(correlation.shape[0], 1)]


def random_effects_meta(effects, standard_errors):
    effects = np.asarray(effects, dtype=float)
    standard_errors = np.asarray(standard_errors, dtype=float)
    variances = standard_errors ** 2
    weights = 1 / variances
    fixed = np.sum(weights * effects) / np.sum(weights)
    q = np.sum(weights * (effects - fixed) ** 2)
    df = len(effects) - 1
    denominator = np.sum(weights) - np.sum(weights ** 2) / np.sum(weights)
    tau2 = max(0.0, (q - df) / denominator)
    random_weights = 1 / (variances + tau2)
    pooled = np.sum(random_weights * effects) / np.sum(random_weights)
    se = np.sqrt(1 / np.sum(random_weights))
    z = pooled / se
    return {
        "n_species": len(effects), "pooled_effect": pooled, "se": se,
        "ci95_low": pooled - stats.norm.ppf(0.975) * se,
        "ci95_high": pooled + stats.norm.ppf(0.975) * se,
        "pvalue": 2 * stats.norm.sf(abs(z)), "tau2": tau2, "q": q, "q_df": df,
        "q_pvalue": stats.chi2.sf(q, df),
    }
