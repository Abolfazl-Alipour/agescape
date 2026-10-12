"""
Determinism utilities for Conserved Aging Explorer.
"""

import os
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
