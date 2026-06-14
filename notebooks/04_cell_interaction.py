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
liana_res.sort_values("magnitude_rank").head(10)

# %% [markdown]
# ## 3. Dot Plot: Top Interactions Involving T Cells

# %%
# Filter to interactions involving CD8 T cells
cd8_interactions = liana_res[
    (liana_res["source"] == "CD8 T cell") |
    (liana_res["target"] == "CD8 T cell")
].sort_values("magnitude_rank").head(30)

fig = li.pl.dotplot(
    adata_sub,
    uns_key="liana_res",
    colour="magnitude_rank",
    size="magnitude_rank",
    source_labels=["CD8 T cell"],
    target_labels=valid_types,
    top_n=20,
    orderby="magnitude_rank",
    orderby_ascending=True,
    figure_size=(10, 6),
    return_fig=True,
)
# liana returns a plotnine ggplot object — save with .save()
fig.save("../figures/04_cd8_interactions_dotplot.png", dpi=150)

# %% [markdown]
# ## 4. Interaction Network (Circle Plot)

# %%
ax = li.pl.circle_plot(
    adata_sub,
    uns_key="liana_res",
    groupby="cell_type",
    source_labels=valid_types,
    target_labels=valid_types,
    top_n=200,
    orderby="magnitude_rank",
    orderby_ascending=True,
    score_key="magnitude_rank",
)
ax.get_figure().savefig("../figures/04_interaction_network.png", dpi=150, bbox_inches="tight")
plt.close()

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
    .sort_values("magnitude_rank")
    .head(15)[["source", "target", "ligand_complex", "receptor_complex", "magnitude_rank"]]
)
non_exp_top = (
    non_expanded_res[non_expanded_res["source"].str.contains("non-expanded")]
    .sort_values("magnitude_rank")
    .head(15)[["source", "target", "ligand_complex", "receptor_complex", "magnitude_rank"]]
)

fig, axes = plt.subplots(1, 2, figsize=(16, 6))
for ax, df, title in zip(
    axes,
    [exp_top, non_exp_top],
    ["Expanded CD8 T cell interactions", "Non-expanded CD8 T cell interactions"],
):
    df["interaction"] = df["ligand_complex"] + " → " + df["receptor_complex"]
    ax.barh(df["interaction"], -np.log10(df["magnitude_rank"] + 1e-6),
            color="steelblue", edgecolor="white")
    ax.set_xlabel("-log10(aggregate rank)")
    ax.set_title(title)
    ax.invert_yaxis()

plt.tight_layout()
savefig("04_expanded_vs_nonexpanded_interactions")
plt.show()

# %% [markdown]
# ## 6. Checkpoint and Activation Signaling: Expanded vs Non-expanded
#
# The "top interactions" comparison above is dominated by interactions
# (e.g. HLA-A/C -> CD8A/CD8B) that are strong in *any* CD8 T cell
# population — they aren't specific to clonal expansion and just confirm
# CD8 T cells engage MHC class I, which is expected by definition.
#
# Here we look at two signals that *do* differ by expansion status:
# - **LAG3 signaling** — LAG3 is a T cell exhaustion checkpoint receptor;
#   its ligands include MHC-II molecules (HLA-DRB1/3/4/5, HLA-DQB1) and
#   LGALS3 (galectin-3). More/stronger LAG3-receptor interactions would be
#   consistent with chronic-antigen-driven exhaustion accumulating in
#   expanded clones.
# - **CCL5 -> SDC1** (CD8 T -> Plasma cell) — CCL5 is a chemokine produced
#   by activated/cytotoxic effector T cells; SDC1 (CD138) marks plasma
#   cells. A stronger signal here would suggest expanded clones are more
#   "activated" towards the plasma cell compartment.

# %%
LAG3_LIGANDS = ["HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "HLA-DQB1", "LGALS3"]

def lag3_interactions(df, source_label):
    sub = df[
        (df["source"] == source_label)
        & (df["receptor_complex"] == "LAG3")
        & (df["ligand_complex"].isin(LAG3_LIGANDS))
    ]
    return sub.sort_values("magnitude_rank")

lag3_exp = lag3_interactions(expanded_res, "CD8 T (expanded)")
lag3_non = lag3_interactions(non_expanded_res, "CD8 T (non-expanded)")

print(f"LAG3-receptor interactions from expanded CD8 T cells: {len(lag3_exp)}")
print(f"LAG3-receptor interactions from non-expanded CD8 T cells: {len(lag3_non)}")

ccl5_exp = expanded_res[
    (expanded_res["source"] == "CD8 T (expanded)")
    & (expanded_res["ligand_complex"] == "CCL5")
    & (expanded_res["receptor_complex"] == "SDC1")
]
ccl5_non = non_expanded_res[
    (non_expanded_res["source"] == "CD8 T (non-expanded)")
    & (non_expanded_res["ligand_complex"] == "CCL5")
    & (non_expanded_res["receptor_complex"] == "SDC1")
]
print(
    f"\nCCL5 -> SDC1 magnitude_rank (lower = stronger): "
    f"expanded={ccl5_exp['magnitude_rank'].values[0]:.4f}, "
    f"non-expanded={ccl5_non['magnitude_rank'].values[0]:.4f}"
)

# %%
fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

# Panel 1: breadth of LAG3-receptor signaling (count of ranked interactions)
ax = axes[0]
counts = [len(lag3_exp), len(lag3_non)]
ax.bar(["Expanded", "Non-expanded"], counts, color=["#E64B35", "#4DBBD5"])
ax.set_ylabel("# LAG3-receptor interactions\n(any target cell type)")
ax.set_title("LAG3 (exhaustion checkpoint)\nsignaling breadth")
for i, c in enumerate(counts):
    ax.text(i, c + 0.3, str(c), ha="center")

# Panel 2: CCL5 -> SDC1 strength (lower magnitude_rank = stronger)
ax = axes[1]
ranks = [ccl5_exp["magnitude_rank"].values[0], ccl5_non["magnitude_rank"].values[0]]
ax.bar(["Expanded", "Non-expanded"], ranks, color=["#E64B35", "#4DBBD5"])
ax.set_ylabel("magnitude_rank (lower = stronger)")
ax.set_title("CCL5 -> SDC1\n(CD8 T -> Plasma cell)")
for i, r in enumerate(ranks):
    ax.text(i, r + max(ranks) * 0.02, f"{r:.4f}", ha="center")

plt.tight_layout()
savefig("04_expansion_checkpoint_activation")
plt.show()

# %% [markdown]
# **Result:** clonally-expanded CD8 T cells show roughly **2x as many**
# ranked LAG3-receptor interactions as non-expanded cells, and their
# **CCL5 -> SDC1** signal to plasma cells is ranked nearly an order of
# magnitude stronger. Together this suggests expanded (likely
# antigen-experienced) clones are simultaneously more "activated"
# (CCL5) and further along towards exhaustion (LAG3 engagement) — a
# pattern consistent with chronic antigen stimulation in the tumour
# microenvironment, and a more dataset-specific signal than the generic
# HLA-CD8 interactions above.

# %% [markdown]
# ## Summary
#
# Key findings from the cell-cell interaction analysis:
# - Ligand-receptor pairs mediating T cell activation, exhaustion, and
#   suppression in the tumour microenvironment were identified (dominated,
#   as expected, by HLA class I -> CD8A/CD8B)
# - **Expanded CD8 T cell clones show a distinct signaling profile**: ~2x
#   more LAG3 (exhaustion checkpoint) interactions and a ~10x stronger
#   CCL5 -> SDC1 (activation) signal to plasma cells than non-expanded
#   clones — consistent with chronic antigen-driven activation and
#   exhaustion (Section 6)
# - These interaction signatures provide candidate pathways for
#   therapeutic intervention in leukemia and other haematological cancers

# %%
print("Notebook 4 complete. All figures saved to figures/")
