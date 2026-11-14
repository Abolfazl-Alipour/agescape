from __future__ import annotations

import numpy as np
import scipy.sparse as sp


def decade(age: int) -> str:
    if not 0 <= age <= 120:
        raise ValueError(f"Invalid age: {age}")
    return f"{age // 10 * 10}s"


def validate_raw_counts(matrix: sp.spmatrix) -> dict:
    values = matrix.data
    if not np.isfinite(values).all():
        raise ValueError("Non-finite raw count")
    if (values < 0).any():
        raise ValueError("Negative raw count")
    if not np.equal(values, np.floor(values)).all():
        raise ValueError("Non-integer raw count")
    return {
        "all_stored_counts_finite_nonnegative_integer": True,
        "stored_nonzero_entries": int(matrix.nnz),
        "raw_umi_sum": int(values.sum(dtype=np.int64)),
    }


def summarize_metric(values: np.ndarray, prefix: str) -> dict:
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or len(values) == 0 or not np.isfinite(values).all():
        raise ValueError(f"Invalid values for {prefix}")
    return {
        f"{prefix}_min": float(values.min()),
        f"{prefix}_p25": float(np.quantile(values, 0.25)),
        f"{prefix}_median": float(np.median(values)),
        f"{prefix}_mean": float(values.mean()),
        f"{prefix}_p75": float(np.quantile(values, 0.75)),
        f"{prefix}_max": float(values.max()),
    }
