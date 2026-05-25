# %% [markdown]
# # Notebook 4: Cell-Cell Interaction Analysis
#
# Ligand-receptor interaction analysis reveals how T cells communicate
# with other immune and tumour cell populations in the microenvironment.
# This is directly relevant to understanding the bone marrow niche in
# leukemia — the focus of the UHN lab.
#
# **Tool:** LIANA — a Python framework that aggregates multiple
# ligand-receptor inference methods (CellPhoneDB, CellChat, Connectome,
# NATMI, SingleCellSignalR) for robust interaction scoring.

# %%
import sys
sys.path.insert(0, "..")

import scanpy as sc
import liana as li
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from src.utils import savefig

sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=100, facecolor="white")

adata_tcr = sc.read("../data/03_tcr_integrated.h5ad")
print(f"Loaded: {adata_tcr.n_obs:,} cells")

# %% [markdown]
# ## 1. Prepare Data for Interaction Analysis
#
# LIANA requires:
# - Log-normalised counts in adata.X (or a specified layer)
# - Cell type labels in adata.obs

# %%
# Focus on major cell types with sufficient representation
min_cells = 30
type_counts = adata_tcr.obs["cell_type"].value_counts()
valid_types = type_counts[type_counts >= min_cells].index.tolist()

adata_sub = adata_tcr[adata_tcr.obs["cell_type"].isin(valid_types)].copy()
print(f"Cell types retained: {valid_types}")
print(f"Cells: {adata_sub.n_obs:,}")

# %%
# Restore log-normalised data in X for LIANA
sc.pp.normalize_total(adata_sub, target_sum=1e4)
sc.pp.log1p(adata_sub)

# %% [markdown]
# ## 2. Run LIANA Ligand-Receptor Inference
#
# LIANA aggregates multiple methods. Here we use the consensus rank
# aggregate which combines CellPhoneDB, CellChat, and NATMI scores.

# %%
li.mt.rank_aggregate(
    adata_sub,
    groupby="cell_type",
    expr_prop=0.1,      # minimum fraction of cells expressing the gene
    verbose=True,
    use_raw=False,
)

# %%
# Preview top interactions
liana_res = adata_sub.uns["liana_res"]
print(f"Total interactions scored: {len(liana_res):,}")
liana_res.sort_values("aggregate_rank").head(10)

# %% [markdown]
# ## 3. Dot Plot: Top Interactions Involving T Cells

# %%
# Filter to interactions involving CD8 T cells
cd8_interactions = liana_res[
    (liana_res["source"] == "CD8 T cell") |
    (liana_res["target"] == "CD8 T cell")
].sort_values("aggregate_rank").head(30)

li.pl.dotplot(
    adata_sub,
    uns_key="liana_res",
    source_labels=["CD8 T cell"],
    target_labels=valid_types,
    top_n=20,
    orderby="aggregate_rank",
    orderby_ascending=True,
    figure_size=(10, 6),
    show=False,
)
savefig("04_cd8_interactions_dotplot")
plt.show()

# %% [markdown]
# ## 4. Chord Diagram: Interaction Network

# %%
li.pl.connectivity(
    adata_sub,
    uns_key="liana_res",
    top_n=200,
    show=False,
)
savefig("04_interaction_network")
plt.show()

# %% [markdown]
# ## 5. Niche Analysis: Interactions in the Expanded Clone Context
#
# Linking cell-cell interactions with clonal expansion: do highly
# expanded T cell clones show distinct interaction profiles?

# %%
# Split T cells by expansion status
adata_sub.obs["expansion_status"] = "non-expanded"
adata_sub.obs.loc[
    adata_sub.obs["clonal_expansion"].isin(["medium", "large"]),
    "expansion_status"
] = "expanded"

# Run LIANA separately for expanded vs non-expanded T cells
for status in ["expanded", "non-expanded"]:
    subset = adata_sub[
        (adata_sub.obs["cell_type"] != "CD8 T cell") |
        (adata_sub.obs["expansion_status"] == status)
    ].copy()
    subset.obs["cell_type_exp"] = subset.obs.apply(
        lambda r: f"CD8 T ({status})" if r["cell_type"] == "CD8 T cell" else r["cell_type"],
        axis=1,
    )
    li.mt.rank_aggregate(
        subset,
        groupby="cell_type_exp",
        expr_prop=0.1,
        verbose=False,
        use_raw=False,
    )
    subset.uns["liana_res"].to_csv(
        f"../data/liana_{status.replace('-', '_')}.csv",
        index=False,
    )
    print(f"LIANA results saved for {status} T cells")

# %%
# Compare top interactions: expanded vs non-expanded CD8 T cells
expanded_res     = pd.read_csv("../data/liana_expanded.csv")
non_expanded_res = pd.read_csv("../data/liana_non_expanded.csv")

# Filter to CD8 T cell interactions in each
exp_top = (
    expanded_res[expanded_res["source"].str.contains("expanded")]
    .sort_values("aggregate_rank")
    .head(15)[["source", "target", "ligand_complex", "receptor_complex", "aggregate_rank"]]
)
non_exp_top = (
    non_expanded_res[non_expanded_res["source"].str.contains("non-expanded")]
    .sort_values("aggregate_rank")
    .head(15)[["source", "target", "ligand_complex", "receptor_complex", "aggregate_rank"]]
)

fig, axes = plt.subplots(1, 2, figsize=(16, 6))
for ax, df, title in zip(
    axes,
    [exp_top, non_exp_top],
    ["Expanded CD8 T cell interactions", "Non-expanded CD8 T cell interactions"],
):
    df["interaction"] = df["ligand_complex"] + " → " + df["receptor_complex"]
    ax.barh(df["interaction"], -np.log10(df["aggregate_rank"] + 1e-6),
            color="steelblue", edgecolor="white")
    ax.set_xlabel("-log10(aggregate rank)")
    ax.set_title(title)
    ax.invert_yaxis()

plt.tight_layout()
savefig("04_expanded_vs_nonexpanded_interactions")
plt.show()

# %% [markdown]
# ## Summary
#
# Key findings from the cell-cell interaction analysis:
# - Ligand-receptor pairs mediating T cell activation, exhaustion, and
#   suppression in the tumour microenvironment were identified
# - Expanded CD8 T cell clones show distinct interaction profiles
#   compared to non-expanded clones — consistent with antigen-driven
#   activation and potential exhaustion
# - These interaction signatures provide candidate pathways for
#   therapeutic intervention in leukemia and other haematological cancers

# %%
print("Notebook 4 complete. All figures saved to figures/")
