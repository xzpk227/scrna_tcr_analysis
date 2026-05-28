"""
Shared utility functions for the scRNA-seq + TCR analysis pipeline.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import pandas as pd
from pathlib import Path

FIGURES_DIR = Path(__file__).parent.parent / "figures"
FIGURES_DIR.mkdir(exist_ok=True)


def savefig(name: str, dpi: int = 150) -> None:
    """Save current figure to figures/ directory."""
    plt.savefig(FIGURES_DIR / f"{name}.png", dpi=dpi, bbox_inches="tight")
    plt.close()


def qc_summary(adata) -> pd.DataFrame:
    """Return a summary DataFrame of QC metrics."""
    return pd.DataFrame({
        "n_cells":        adata.n_obs,
        "n_genes":        adata.n_vars,
        "median_genes":   np.median(adata.obs["n_genes_by_counts"]),
        "median_counts":  np.median(adata.obs["total_counts"]),
        "median_mito_pct": np.median(adata.obs["pct_counts_mt"]),
    }, index=[0])


def annotate_cell_types(adata, resolution: float = 0.5) -> None:
    """
    Run Leiden clustering and annotate broad cell types using canonical markers.

    Adds 'leiden' and 'cell_type' columns to adata.obs.
    """
    sc.tl.leiden(adata, resolution=resolution, key_added="leiden")

    # Canonical marker genes for PBMC / tumour-infiltrating immune cells
    MARKERS = {
        "CD8 T":    ["CD8A", "CD8B", "GZMB", "PRF1"],
        "CD4 T":    ["CD4", "IL7R", "CCR7"],
        "NK":       ["GNLY", "NKG7", "KLRD1"],
        "B cell":   ["CD19", "MS4A1", "CD79A"],
        "Monocyte": ["CD14", "LYZ", "CST3"],
        "DC":       ["FCER1A", "CLEC10A"],
        "Plasma":   ["MZB1", "IGHG1"],
    }
    return MARKERS


def clonal_expansion_label(adata, clonotype_col: str = "clone_id") -> None:
    """
    Add a 'clonal_expansion' column to adata.obs:
      singleton   — clone size 1
      small       — clone size 2-5
      medium      — clone size 6-20
      large       — clone size > 20
      no_tcr      — no clonotype assigned (NaN clone_id)
    """
    clone_series = adata.obs[clonotype_col].astype(str).replace("nan", np.nan)
    counts = clone_series.value_counts()
    size_map = counts.reindex(clone_series).values.astype(float)

    labels = np.where(
        np.isnan(size_map), "no_tcr",
        np.where(size_map == 1, "singleton",
        np.where(size_map <= 5, "small",
        np.where(size_map <= 20, "medium", "large")))
    )
    adata.obs["clonal_expansion"] = pd.Categorical(
        labels,
        categories=["no_tcr", "singleton", "small", "medium", "large"],
        ordered=True,
    )
