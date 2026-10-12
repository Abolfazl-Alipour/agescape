"""
Consistent publication-grade plotting style for Conserved Aging Explorer.
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns


def set_pub_style() -> None:
    """Set global matplotlib and seaborn publication styling."""
    plt.rcParams.update({
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "figure.titlesize": 12,
        "axes.edgecolor": "#333333",
        "axes.linewidth": 0.8,
        "grid.color": "#E5E5E5",
        "grid.linestyle": "--",
        "grid.linewidth": 0.5,
        "axes.grid": False,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
    })
    sns.set_palette("colorblind")
