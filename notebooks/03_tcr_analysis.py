# %% [markdown]
# # Notebook 3: TCR Analysis and Clonotype-State Linking
#
# Integrates T cell receptor (TCR) sequencing with scRNA-seq annotations
# to link clonal identity with transcriptional state — a central question
# in tumour immunology and leukemia biology.
#
# Key analyses:
# 1. TCR quality control and chain pairing
# 2. Clonotype definition (identical alpha + beta CDR3)
# 3. Clonal expansion across cell types
# 4. **Clonotype → transcriptional state linkage** (core requirement)

# %%
import sys
sys.path.insert(0, "..")

import scanpy as sc
import scirpy as ir
import muon as mu
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from src.utils import savefig, clonal_expansion_label

sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=100, facecolor="white")
sc.settings.autoshow = False   # prevent auto-display; we save manually

# %% [markdown]
# ## 1. Load Data
#
# In scirpy >= 0.16, AIRR data is stored in `adata.obsm["airr"]`.
# We work with the full MuData object and update the GEX modality
# with our preprocessed, annotated data from Notebook 2.

# %%
# Load our preprocessed + annotated GEX data
adata = sc.read("../data/02_annotated.h5ad")
print(f"Preprocessed GEX: {adata.n_obs:,} cells")

# Load the original wu2020 MuData (cached from Notebook 1)
mdata = ir.datasets.wu2020()
print(f"MuData modalities: {list(mdata.mod.keys())}")

# Run index_chains FIRST — required before any scirpy tl/pl functions
ir.pp.index_chains(mdata)
print("Chain indices built.")

# %%
# Align barcodes
shared = adata.obs_names.intersection(mdata["gex"].obs_names)
print(f"Shared barcodes: {len(shared):,}")

# Subset MuData to shared cells (index_chains already built on full mdata)
mdata_sub = mdata[shared].copy()

# Transfer annotations — write to both gex obs AND global mdata obs
# so scirpy functions that look in either place can find them
ann_cols = ["leiden", "cell_type", "Exhausted_score", "Effector_score", "Naive_score"]
for col in ann_cols:
    if col in adata.obs.columns:
        vals = adata.obs.loc[shared, col].values
        mdata_sub["gex"].obs[col] = vals
        # Also add to global obs (reindexed to full MuData obs)
        mdata_sub.obs[col] = (
            adata.obs[col]
            .reindex(mdata_sub.obs_names)
            .values
        )

# Transfer UMAP embedding
mdata_sub["gex"].obsm["X_umap"] = adata[shared].obsm["X_umap"]

print(f"Working dataset: {mdata_sub.n_obs:,} cells")
print(f"Annotations in gex: {[c for c in ann_cols if c in mdata_sub['gex'].obs.columns]}")
print(f"Annotations in mdata.obs: {[c for c in ann_cols if c in mdata_sub.obs.columns]}")

# %% [markdown]
# ## 2. TCR Chain Quality Control
#
# Not all cells have complete alpha+beta chain pairs. We assess pairing
# quality before defining clonotypes.

# %%
ir.tl.chain_qc(mdata_sub)

# chain_qc writes to mdata["airr"].obs — align by barcode before copying
airr_obs = mdata_sub["airr"].obs
if "receptor_type" in airr_obs.columns:
    # Use reindex to safely align — cells without TCR get NaN
    mdata_sub["gex"].obs["receptor_type"] = (
        airr_obs["receptor_type"]
        .reindex(mdata_sub["gex"].obs_names)
        .values
    )
elif "receptor_type" in mdata_sub.obs.columns:
    mdata_sub["gex"].obs["receptor_type"] = (
        mdata_sub.obs["receptor_type"]
        .reindex(mdata_sub["gex"].obs_names)
        .values
    )

print("Chain QC complete. receptor_type counts:")
print(mdata_sub["gex"].obs["receptor_type"].value_counts(dropna=False))

# %%
ir.pl.group_abundance(
    mdata_sub,
    groupby="receptor_type",
    target_col="cell_type",
)
savefig("03_chain_qc_receptor_type")
plt.show()

# %%
# Keep only cells with paired alpha+beta TCR for clonotype analysis
mdata_paired = mdata_sub[
    mdata_sub["gex"].obs["receptor_type"].astype(str) == "TCR"
].copy()
ir.pp.index_chains(mdata_paired)
print(f"Cells with paired TCR: {mdata_paired.n_obs:,}")

# Force re-transfer all annotations after MuData subset
paired_barcodes = mdata_paired["gex"].obs_names
for col in ann_cols:
    if col in adata.obs.columns:
        mdata_paired["gex"].obs[col] = (
            adata.obs[col].reindex(paired_barcodes).values
        )
print(f"Annotations available: {[c for c in ann_cols if c in mdata_paired['gex'].obs.columns]}")

# %% [markdown]
# ## 3. Clonotype Definition
#
# Cells sharing identical CDR3 amino acid sequences for both TCR alpha
# and beta chains are assigned to the same clonotype.

# %%
ir.tl.define_clonotypes(
    mdata_paired,
    receptor_arms="all",
    dual_ir="primary_only",
)

# In scirpy 0.22, clone_id is stored in mdata["airr"].obs — transfer to gex
for col in ["clone_id", "clone_id_size"]:
    if col in mdata_paired["airr"].obs.columns:
        mdata_paired["gex"].obs[col] = (
            mdata_paired["airr"].obs[col]
            .reindex(mdata_paired["gex"].obs_names)
            .values
        )

n_clonotypes = mdata_paired["gex"].obs["clone_id"].nunique()
print(f"Unique clonotypes defined: {n_clonotypes:,}")
print(f"Cells assigned to a clonotype: "
      f"{mdata_paired['gex'].obs['clone_id'].notna().sum():,}")

# %% [markdown]
# ## 4. Clonal Expansion

# %%
clonal_expansion_label(mdata_paired["gex"], clonotype_col="clone_id")
print(mdata_paired["gex"].obs["clonal_expansion"].value_counts())

# %%
sc.pl.umap(
    mdata_paired["gex"],
    color="clonal_expansion",
    palette={
        "no_tcr":    "#EEEEEE",
        "singleton": "#CCCCCC",
        "small":     "#FDB462",
        "medium":    "#E41A1C",
        "large":     "#8B0000",
    },
    title="Clonal expansion on UMAP",
    show=False,
)
savefig("03_umap_clonal_expansion")
plt.show()

# %%
# Expansion proportion by cell type
exp_ct = (
    mdata_paired["gex"].obs
    .groupby(["cell_type", "clonal_expansion"], observed=True)
    .size()
    .unstack(fill_value=0)
    .apply(lambda x: x / x.sum(), axis=1)
)

exp_ct.plot(
    kind="bar", stacked=True,
    color=["#EEEEEE", "#CCCCCC", "#FDB462", "#E41A1C", "#8B0000"],
    figsize=(8, 5), edgecolor="white",
)
plt.xlabel("Cell type")
plt.ylabel("Proportion")
plt.title("Clonal expansion by cell type")
plt.legend(title="Expansion", bbox_to_anchor=(1.05, 1))
plt.tight_layout()
savefig("03_expansion_by_celltype")
plt.show()

# %% [markdown]
# ## 5. Clonotype → Transcriptional State Linkage
#
# The key analysis: do expanded clones preferentially occupy exhausted,
# effector, or naive transcriptional states?

# %%
gex = mdata_paired["gex"]
T_only = gex[gex.obs["cell_type"].isin(["CD8 T cell", "CD4 T cell"])].copy()

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
expansion_order = ["singleton", "small", "medium", "large"]

for ax, score in zip(axes, ["Exhausted_score", "Effector_score", "Naive_score"]):
    data_plot = [
        T_only.obs.loc[T_only.obs["clonal_expansion"] == exp, score].dropna().values
        for exp in expansion_order
        if exp in T_only.obs["clonal_expansion"].values
    ]
    present_labels = [
        exp for exp in expansion_order
        if exp in T_only.obs["clonal_expansion"].values
    ]
    ax.boxplot(data_plot, labels=present_labels, patch_artist=True,
               boxprops=dict(facecolor="lightblue"),
               medianprops=dict(color="red", linewidth=2))
    ax.set_title(score.replace("_score", "").capitalize() + " score")
    ax.set_xlabel("Clonal expansion")
    ax.set_ylabel("Gene score")

plt.suptitle("Transcriptional state by clonal expansion", fontsize=13)
plt.tight_layout()
savefig("03_state_by_expansion")
plt.show()

# %% [markdown]
# ## 6. Save

# %%
gex.write("../data/03_tcr_integrated.h5ad")
print("Saved: data/03_tcr_integrated.h5ad")
