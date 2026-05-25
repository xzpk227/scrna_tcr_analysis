# %% [markdown]
# # Notebook 3: TCR Analysis and Clonotype-State Linking
#
# This notebook integrates T cell receptor (TCR) sequencing data with the
# scRNA-seq annotations from Notebook 2. The central analysis links T cell
# clonal identity (clonotype) with transcriptional state — a key approach
# for understanding tumour-infiltrating T cell dynamics in leukemia and
# other cancers.
#
# Key analyses:
# 1. TCR data integration with scRNA-seq AnnData
# 2. Clonotype definition and abundance
# 3. Clonal expansion across cell types and functional states
# 4. **Clonotype → transcriptional state linkage** (core UHN requirement)
# 5. VDJ gene usage and CDR3 properties

# %%
import sys
sys.path.insert(0, "..")

import scanpy as sc
import scirpy as ir
import muon as mu
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency

from src.utils import savefig, clonal_expansion_label

sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=100, facecolor="white")

# %%
# Load annotated GEX data and full MuData (contains TCR modality)
adata  = sc.read("../data/02_annotated.h5ad")
mdata  = ir.datasets.wu2020()
print(f"GEX cells: {adata.n_obs:,}")
print(f"TCR cells in MuData: {mdata['airr'].n_obs:,}")

# %% [markdown]
# ## 1. Integrate TCR Data
#
# scirpy merges the AIRR (Adaptive Immune Receptor Repertoire) modality
# with the GEX AnnData by shared cell barcodes. Only cells with both
# GEX and TCR data are retained for clonotype analyses.

# %%
# Transfer AIRR data into adata using shared barcodes
airr_adata = mdata["airr"]
shared_barcodes = adata.obs_names.intersection(airr_adata.obs_names)
print(f"Cells with both GEX and TCR: {len(shared_barcodes):,}")

# Subset to shared cells
adata_tcr = adata[shared_barcodes].copy()
airr_sub  = airr_adata[shared_barcodes].copy()

# Copy AIRR obs columns into GEX adata
for col in airr_sub.obs.columns:
    adata_tcr.obs[col] = airr_sub.obs[col].values

# %% [markdown]
# ## 2. Clonotype Definition
#
# A clonotype is defined by identical CDR3 amino acid sequences for both
# TCR alpha and beta chains. Cells sharing the same clonotype are clonally
# related — they descended from the same ancestral T cell.

# %%
ir.tl.chain_qc(adata_tcr)

# Define clonotypes: cells with identical TRA+TRB CDR3 sequences
ir.tl.define_clonotypes(
    adata_tcr,
    receptor_arms="all",          # require both alpha and beta chains
    dual_ir="primary_only",       # use primary chain only
    sequence="aa",                # amino acid CDR3 matching
)

n_clonotypes = adata_tcr.obs["clone_id"].nunique()
print(f"Unique clonotypes: {n_clonotypes:,}")
print(f"Cells with clonotype: {adata_tcr.obs['clone_id'].notna().sum():,}")

# %%
ir.tl.clonotype_network(adata_tcr, min_cells=2)

# %% [markdown]
# ## 3. Clonal Expansion
#
# Clonal expansion — the proliferation of a single T cell clone — is a
# hallmark of antigen-driven immune responses. Highly expanded clones in
# tumours indicate tumour-reactive T cells.

# %%
# Add clonal expansion category (singleton / small / medium / large)
clonal_expansion_label(adata_tcr, clonotype_col="clone_id")

print(adata_tcr.obs["clonal_expansion"].value_counts())

# %%
# Expansion on UMAP
ir.pl.umap(
    adata_tcr,
    color="clonal_expansion",
    palette={
        "singleton": "#CCCCCC",
        "small":     "#FDB462",
        "medium":    "#E41A1C",
        "large":     "#8B0000",
    },
    title="Clonal expansion",
    show=False,
)
savefig("03_umap_clonal_expansion")
plt.show()

# %%
# Expansion by cell type — stacked bar chart
exp_ct = (
    adata_tcr.obs
    .groupby(["cell_type", "clonal_expansion"], observed=True)
    .size()
    .unstack(fill_value=0)
    .apply(lambda x: x / x.sum(), axis=1)
)

exp_ct.plot(
    kind="bar",
    stacked=True,
    color=["#CCCCCC", "#FDB462", "#E41A1C", "#8B0000"],
    figsize=(8, 5),
    edgecolor="white",
)
plt.xlabel("Cell type")
plt.ylabel("Proportion")
plt.title("Clonal expansion by cell type")
plt.legend(title="Expansion", bbox_to_anchor=(1.05, 1))
plt.tight_layout()
savefig("03_expansion_by_celltype")
plt.show()

# %% [markdown]
# ## 4. Clonotype → Transcriptional State Linkage
#
# **Key analysis:** This links clonal identity with transcriptional
# programme, revealing whether expanded clones preferentially occupy
# specific functional states (e.g., exhausted, effector, naive).
#
# This is the central question in tumour immunology: do expanded,
# tumour-reactive clones become exhausted or remain functional?

# %%
# Score cells for exhaustion and effector states
EXHAUSTION_GENES = ["PDCD1", "HAVCR2", "TIGIT", "LAG3", "CTLA4", "ENTPD1"]
EFFECTOR_GENES   = ["GZMB", "PRF1", "IFNG", "TNF", "GZMK"]
NAIVE_GENES      = ["CCR7", "TCF7", "SELL", "LEF1"]

for name, genes in [
    ("exhaustion", EXHAUSTION_GENES),
    ("effector",   EFFECTOR_GENES),
    ("naive",      NAIVE_GENES),
]:
    valid = [g for g in genes if g in adata_tcr.var_names]
    sc.tl.score_genes(adata_tcr, gene_list=valid, score_name=f"{name}_score")

# %%
# Compare transcriptional scores across expansion groups
T_only = adata_tcr[adata_tcr.obs["cell_type"].isin(["CD8 T cell", "CD4 T cell"])].copy()

fig, axes = plt.subplots(1, 3, figsize=(15, 5))
expansion_order = ["singleton", "small", "medium", "large"]
palette = ["#CCCCCC", "#FDB462", "#E41A1C", "#8B0000"]

for ax, score in zip(axes, ["exhaustion_score", "effector_score", "naive_score"]):
    data_plot = [
        T_only.obs.loc[T_only.obs["clonal_expansion"] == exp, score].dropna().values
        for exp in expansion_order
    ]
    ax.boxplot(data_plot, labels=expansion_order, patch_artist=True,
               boxprops=dict(facecolor="lightblue"),
               medianprops=dict(color="red", linewidth=2))
    ax.set_title(score.replace("_score", "").capitalize() + " score")
    ax.set_xlabel("Clonal expansion")
    ax.set_ylabel("Gene score")

plt.suptitle("Transcriptional state by clonal expansion", fontsize=13)
plt.tight_layout()
savefig("03_state_by_expansion")
plt.show()

# %%
# Heatmap: mean transcriptional scores per clonotype (top 20 expanded)
top_clones = (
    adata_tcr.obs["clone_id"]
    .value_counts()
    .head(20)
    .index
)

clone_scores = (
    adata_tcr.obs[adata_tcr.obs["clone_id"].isin(top_clones)]
    .groupby("clone_id")[["exhaustion_score", "effector_score", "naive_score"]]
    .mean()
)
clone_scores.index = [f"Clone {i+1}" for i in range(len(clone_scores))]

fig, ax = plt.subplots(figsize=(6, 8))
sns.heatmap(
    clone_scores,
    cmap="RdBu_r",
    center=0,
    linewidths=0.5,
    ax=ax,
    cbar_kws={"label": "Mean score"},
)
ax.set_title("Transcriptional state of top 20 expanded clonotypes")
ax.set_xlabel("Functional score")
plt.tight_layout()
savefig("03_clonotype_state_heatmap")
plt.show()

# %% [markdown]
# ## 5. VDJ Gene Usage and CDR3 Properties

# %%
# V gene usage for TRB (beta chain)
ir.pl.vdj_usage(
    adata_tcr,
    full_combination=False,
    show=False,
)
savefig("03_vdj_usage")
plt.show()

# %%
# CDR3 length distribution
ir.pl.spectratype(
    adata_tcr,
    color="cell_type",
    viztype="bar",
    show=False,
)
savefig("03_cdr3_spectratype")
plt.show()

# %% [markdown]
# ## 6. Save TCR-Integrated Data

# %%
adata_tcr.write("../data/03_tcr_integrated.h5ad")
print("Saved: data/03_tcr_integrated.h5ad")
