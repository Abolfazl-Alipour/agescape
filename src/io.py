"""
I/O and dataset utilities for Conserved Aging Explorer.
"""

import hashlib
import json
from pathlib import Path
from typing import Any, Dict, Optional, Tuple
import anndata as ad
import pandas as pd


def compute_sha256(file_path: Path, chunk_size: int = 65536) -> str:
    """Compute SHA-256 digest of a file."""
    hasher = hashlib.sha256()
    with open(file_path, "rb") as f:
        while chunk := f.read(chunk_size):
            hasher.update(chunk)
    return hasher.hexdigest()


def read_manifest(manifest_path: Path) -> Dict[str, Any]:
    """Read a JSON manifest file."""
    with open(manifest_path, "r") as f:
        return json.load(f)


def verify_dataset(file_path: Path, expected_sha256: str) -> bool:
    """Check if file matches expected SHA-256."""
    if not file_path.exists():
        return False
    return compute_sha256(file_path).lower() == expected_sha256.lower()


def strip_obs_strings(adata: ad.AnnData) -> Tuple[ad.AnnData, Dict[str, Dict[str, int]]]:
    """
    Strip leading and trailing whitespace from all string/object/categorical columns in adata.obs.
    Logs before and after value counts for modified columns.
    """
    log_changes = {}
    for col in adata.obs.columns:
        if isinstance(adata.obs[col].dtype, pd.CategoricalDtype):
            orig_categories = list(adata.obs[col].cat.categories)
            stripped_categories = [str(c).strip() if isinstance(c, str) else c for c in orig_categories]
            if orig_categories != stripped_categories:
                before_counts = adata.obs[col].value_counts().to_dict()
                adata.obs[col] = adata.obs[col].astype(str).str.strip().astype("category")
                after_counts = adata.obs[col].value_counts().to_dict()
                log_changes[col] = {"before": before_counts, "after": after_counts}
        elif pd.api.types.is_string_dtype(adata.obs[col]) or pd.api.types.is_object_dtype(adata.obs[col]):
            before_vals = adata.obs[col].copy()
            adata.obs[col] = adata.obs[col].astype(str).str.strip()
            if not before_vals.equals(adata.obs[col]):
                log_changes[col] = {
                    "before": before_vals.value_counts().to_dict(),
                    "after": adata.obs[col].value_counts().to_dict(),
                }
    return adata, log_changes


def load_anndata(file_path: Path, strip_strings: bool = True, backed: Optional[str] = None) -> ad.AnnData:
    """Load AnnData object, optionally applying string stripping."""
    adata = ad.read_h5ad(file_path, backed=backed)
    if strip_strings and backed is None:
        adata, _ = strip_obs_strings(adata)
    return adata
