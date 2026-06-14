# %% [markdown]
# # Notebook 1: QC and Preprocessing
#
# This notebook performs quality control, normalization, dimensionality
# reduction, and clustering on the Wu et al. (2020) tumour-infiltrating
# lymphocyte dataset, which contains paired scRNA-seq and TCR sequencing
# data from cancer patients.
#
# **Dataset:** Wu TD et al. Nature 2020 — peripheral and tumour-infiltrating
# T cells from patients with basal cell carcinoma.

# %%
import sys
sys.path.insert(0, "..")

import scanpy as sc
import scirpy as ir
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from src.utils import savefig, qc_summary

sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=100, facecolor="white")

# %% [markdown]
# ## 1. Load Dataset
#
# scirpy provides the Wu et al. 2020 dataset as a built-in download.
# It returns a MuData object containing both GEX and TCR modalities.

# %%
mdata = ir.datasets.wu2020()
mdata

# %%
# Extract the gene expression AnnData
adata = mdata["gex"].copy()
print(f"Cells: {adata.n_obs:,}  |  Genes: {adata.n_vars:,}")

# %% [markdown]
# ## 2. Quality Control

# %%
# Identify mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True,
)

# %%
def plot_qc_distributions(adata, suptitle):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    axes[0].hist(adata.obs["n_genes_by_counts"], bins=50, color="steelblue", edgecolor="white")
    axes[0].set_xlabel("Genes per cell")
    axes[0].set_ylabel("Cells")
    axes[0].set_title("Genes per cell")

    axes[1].hist(adata.obs["total_counts"], bins=50, color="salmon", edgecolor="white")
    axes[1].set_xlabel("UMI counts per cell")
    axes[1].set_title("UMI counts per cell")

    axes[2].hist(adata.obs["pct_counts_mt"], bins=50, color="mediumpurple", edgecolor="white")
    axes[2].set_xlabel("% mitochondrial counts")
    axes[2].set_title("Mitochondrial fraction")

    fig.suptitle(suptitle)
    plt.tight_layout()
    return fig

plot_qc_distributions(adata, f"QC distributions — before filtering (n={adata.n_obs:,} cells)")
savefig("01_qc_distributions_before")
plt.show()

# %% [markdown]
# ### QC Thresholds
#
# Standard thresholds for PBMC / tumour-infiltrating lymphocytes:
# - Min genes: 200 (removes empty droplets)
# - Max genes: 6,000 (removes doublets)
# - Max mt%: 20% (removes dying cells)

# %%
print(f"Cells before filtering: {adata.n_obs:,}")

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_cells(adata, max_genes=6000)
sc.pp.filter_genes(adata, min_cells=3)

adata = adata[adata.obs["pct_counts_mt"] < 20].copy()

print(f"Cells after filtering:  {adata.n_obs:,}")
print(qc_summary(adata).to_string(index=False))

# %%
plot_qc_distributions(adata, f"QC distributions — after filtering (n={adata.n_obs:,} cells)")
savefig("01_qc_distributions_after")
plt.show()

# %% [markdown]
# ## 3. Normalization and Feature Selection

# %%
# Save raw counts before normalization
adata.layers["counts"] = adata.X.copy()

# Normalize to 10,000 counts per cell then log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata  # store log-normalized data

# Highly variable gene selection
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=3000,
    flavor="seurat_v3",
    layer="counts",
)

print(f"Highly variable genes selected: {adata.var['highly_variable'].sum():,}")

# %%
sc.pl.highly_variable_genes(adata, show=False)
savefig("01_hvg")
plt.show()

# %% [markdown]
# ## 4. Dimensionality Reduction and Clustering

# %%
adata_hvg = adata[:, adata.var["highly_variable"]].copy()

sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, svd_solver="arpack", n_comps=50)

# %%
sc.pl.pca_variance_ratio(adata_hvg, n_pcs=50, show=False)
savefig("01_pca_variance")
plt.show()

# %%
# Build neighbourhood graph and compute UMAP
sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata_hvg)
sc.tl.leiden(adata_hvg, resolution=0.5, key_added="leiden")

# Transfer embeddings back to full adata
adata.obsm["X_pca"]  = adata_hvg.obsm["X_pca"]
adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]
adata.obsp            = adata_hvg.obsp
adata.uns             = adata_hvg.uns
adata.obs["leiden"]   = adata_hvg.obs["leiden"]

# %%
sc.pl.umap(adata, color="leiden", legend_loc="on data",
           title="Leiden clusters (res=0.5)", show=False)
savefig("01_umap_leiden")
plt.show()

# %% [markdown]
# ## 5. Save Processed Data

# %%
adata.write("../data/01_preprocessed.h5ad")
print("Saved: data/01_preprocessed.h5ad")
