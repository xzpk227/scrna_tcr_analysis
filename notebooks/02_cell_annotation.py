# %% [markdown]
# # Notebook 2: Cell Type Annotation
#
# Using canonical marker genes to annotate major immune cell populations,
# with a focus on T cell subsets relevant to tumour immunology.

# %%
import sys
sys.path.insert(0, "..")

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from src.utils import savefig

sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=100, facecolor="white")

adata = sc.read("../data/01_preprocessed.h5ad")
print(f"Loaded: {adata.n_obs:,} cells, {adata.n_vars:,} genes")

# %% [markdown]
# ## 1. Canonical Marker Genes
#
# Marker genes for major immune cell populations in tumour-infiltrating
# lymphocyte datasets.

# %%
MARKERS = {
    "CD8 T cell":   ["CD8A", "CD8B", "GZMB", "PRF1", "IFNG"],
    "CD4 T cell":   ["CD4", "IL7R", "CCR7", "FOXP3"],
    "NK cell":      ["GNLY", "NKG7", "KLRD1", "NCAM1"],
    "B cell":       ["CD19", "MS4A1", "CD79A", "CD79B"],
    "Monocyte":     ["CD14", "LYZ", "CST3", "FCGR3A"],
    "DC":           ["FCER1A", "CLEC10A", "ITGAX"],
    "Plasma cell":  ["MZB1", "IGHG1", "JCHAIN"],
}

# %%
# Filter markers to genes present in the dataset
MARKERS_PRESENT = {
    ct: [g for g in genes if g in adata.var_names]
    for ct, genes in MARKERS.items()
    if any(g in adata.var_names for g in genes)
}

# Dotplot of marker genes across Leiden clusters
sc.pl.dotplot(
    adata,
    var_names=MARKERS_PRESENT,
    groupby="leiden",
    standard_scale="var",
    show=False,
)
savefig("02_marker_dotplot")
plt.show()

# %%
# Violin plots for key discriminating markers
key_markers = ["CD8A", "CD4", "NKG7", "CD19", "CD14", "MZB1"]
sc.pl.stacked_violin(
    adata,
    var_names=key_markers,
    groupby="leiden",
    show=False,
)
savefig("02_marker_violin")
plt.show()

# %% [markdown]
# ## 2. Assign Cell Type Labels
#
# Based on marker expression, assign broad cell type labels to each
# Leiden cluster. Clusters with ambiguous identity are labelled as
# "Unknown" for downstream exclusion or re-clustering.

# %%
# Annotation based on per-cluster mean marker expression (see below), not just
# dotplot inspection. Two corrections to the original dotplot-based reading:
#
# 1. This dataset is CD3-sorted and T-cell-dominated: CD3E is moderately-to-
#    highly expressed in nearly every cluster, so a high score on the
#    "NK cell" marker panel (GNLY/NKG7/KLRD1) does not distinguish true NK
#    cells from cytotoxic CD8 effector T cells, which also express these
#    granzyme/cytotoxicity genes highly.
# 2. Clusters originally labelled Monocyte (7), B cell (8), DC (12), and
#    Plasma cell (13) have essentially zero expression of the
#    lineage-defining markers for those types (CD14/LYZ, CD19/MS4A1/CD79A,
#    FCER1A/CLEC10A, MZB1 respectively) but very high CD3E + CD8A/CD8B —
#    they are CD8 T cell clusters. The true B cell and Plasma cell clusters
#    (high MS4A1/CD79A/CD19 and MZB1 respectively) were instead sitting in
#    the originally-unlabelled clusters 16 and 17.
CLUSTER_ANNOTATION = {
    "0":  "CD8 T cell",
    "1":  "CD4 T cell",
    "2":  "CD4 T cell",
    "3":  "CD4 T cell",
    "4":  "CD8 T cell",
    "5":  "CD4 T cell",
    "6":  "CD8 T cell",
    "7":  "CD8 T cell",
    "8":  "CD8 T cell",
    "9":  "CD8 T cell",
    "10": "CD8 T cell",
    "11": "NK cell",
    "12": "CD8 T cell",
    "13": "CD8 T cell",
    "14": "Monocyte/DC",
    "15": "CD8 T cell",
    "16": "B cell",
    "17": "Plasma cell",
}

adata.obs["cell_type"] = (
    adata.obs["leiden"]
    .map(CLUSTER_ANNOTATION)
    .fillna("Unknown")
    .astype("category")
)

print(adata.obs["cell_type"].value_counts())

# %%
PALETTE = {
    "CD8 T cell":  "#E64B35",
    "CD4 T cell":  "#4DBBD5",
    "NK cell":     "#00A087",
    "B cell":      "#3C5488",
    "Monocyte/DC": "#F39B7F",
    "Plasma cell": "#91D1C2",
    "Unknown":     "#B0B0B0",
}

sc.pl.umap(
    adata,
    color="cell_type",
    palette=PALETTE,
    legend_loc="right margin",
    title="Cell type annotation",
    show=False,
)
savefig("02_umap_celltypes")
plt.show()

# %% [markdown]
# ## 3. T Cell Subset Analysis
#
# Focus on T cells for downstream TCR integration. Identify naive,
# memory, effector, and regulatory T cell subsets.

# %%
T_SUBSET_MARKERS = {
    "Naive":       ["CCR7", "TCF7", "SELL", "LEF1"],
    "Memory":      ["IL7R", "GPR183", "CXCR3"],
    "Effector":    ["GZMB", "GZMK", "PRF1", "IFNG", "TNF"],
    "Exhausted":   ["PDCD1", "HAVCR2", "TIGIT", "LAG3", "CTLA4"],
    "Regulatory":  ["FOXP3", "IL2RA", "IKZF2"],
}

t_cells = adata[adata.obs["cell_type"].isin(["CD8 T cell", "CD4 T cell"])].copy()
print(f"T cells: {t_cells.n_obs:,}")

T_SUBSET_MARKERS_PRESENT = {
    state: [g for g in genes if g in adata.var_names]
    for state, genes in T_SUBSET_MARKERS.items()
    if any(g in adata.var_names for g in genes)
}

sc.pl.dotplot(
    t_cells,
    var_names=T_SUBSET_MARKERS_PRESENT,
    groupby="cell_type",
    standard_scale="var",
    title="T cell subset markers",
    show=False,
)
savefig("02_tcell_subsets_dotplot")
plt.show()

# %%
# Score each cell for T cell functional states
for state, markers in T_SUBSET_MARKERS.items():
    valid = [g for g in markers if g in adata.var_names]
    if valid:
        sc.tl.score_genes(adata, gene_list=valid, score_name=f"{state}_score")

sc.pl.umap(
    adata,
    color=[f"{s}_score" for s in T_SUBSET_MARKERS],
    ncols=3,
    cmap="RdYlBu_r",
    show=False,
)
savefig("02_tcell_state_scores")
plt.show()

# %% [markdown]
# ## 4. Save Annotated Data

# %%
adata.write("../data/02_annotated.h5ad")
print("Saved: data/02_annotated.h5ad")
