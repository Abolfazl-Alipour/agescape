"""
Determinism utilities for Conserved Aging Explorer.
"""

import os
import sys

# CRITICAL: Thread and hash seeds must be set at module load time before C-libraries initialize
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
if "numba" not in sys.modules:
    os.environ.setdefault("NUMBA_NUM_THREADS", "1")
os.environ.setdefault("PYTHONHASHSEED", "42")

import random
import numpy as np


def set_seeds(seed: int = 42) -> None:
    """Set random seeds across Python, NumPy, and environment for reproducible execution."""
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)
    try:
        import torch
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
    except ImportError:
        pass
