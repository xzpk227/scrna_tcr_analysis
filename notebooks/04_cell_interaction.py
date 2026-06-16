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
# Figure A: LAG3 signaling breadth
fig, ax = plt.subplots(figsize=(5, 4.5))
counts = [len(lag3_exp), len(lag3_non)]
ax.bar(["Expanded", "Non-expanded"], counts, color=["#E64B35", "#4DBBD5"])
ax.set_ylabel("# LAG3-receptor interactions\n(any target cell type)")
ax.set_title("LAG3 (exhaustion checkpoint)\nsignaling breadth")
for i, c in enumerate(counts):
    ax.text(i, c + 0.3, str(c), ha="center")
plt.tight_layout()
savefig("04_lag3_signaling_breadth")
plt.show()

# Figure B: CCL5 -> SDC1 strength
fig, ax = plt.subplots(figsize=(5, 4.5))
ranks = [ccl5_exp["magnitude_rank"].values[0], ccl5_non["magnitude_rank"].values[0]]
ax.bar(["Expanded", "Non-expanded"], ranks, color=["#E64B35", "#4DBBD5"])
ax.set_ylabel("magnitude_rank (lower = stronger)")
ax.set_title("CCL5 → SDC1\n(CD8 T → Plasma cell, expanded vs non-expanded)")
for i, r in enumerate(ranks):
    ax.text(i, r + max(ranks) * 0.02, f"{r:.4f}", ha="center")
plt.tight_layout()
savefig("04_ccl5_sdc1_rank")
plt.show()

# %% [markdown]
# **Result:** the CCL5 -> SDC1 signal (CD8 T -> Plasma cell) is ranked
# ~8x stronger in expanded clones (0.010 vs 0.080), with equal group
# sizes (~25k cells each) ruling out a cell-count artefact. Expanded
# clones communicate more strongly with the plasma cell compartment —
# the functional consequence is unclear and requires experimental
# follow-up, but this is a dataset-specific signal beyond the generic
# HLA->CD8 interactions.

# %% [markdown]
# ## 7. CD8 vs CD4 T Cell Sender Profiles
#
# CD8 T cells are primarily cytotoxic effectors; CD4 T cells primarily
# provide "help" via cytokines and co-stimulatory signals to other immune
# cells. Comparing the top ligand-receptor pairs *sent* by each reveals
# whether LIANA recovers this functional division from gene expression alone.

# %%
fig, axes = plt.subplots(1, 2, figsize=(16, 10))

for ax, source in zip(axes, ["CD8 T cell", "CD4 T cell"]):
    if source not in liana_res["source"].values:
        ax.set_title(f"{source} — not in dataset")
        continue
    top = (
        liana_res[liana_res["source"] == source]
        .sort_values("magnitude_rank")
        .head(15)
        .copy()
    )
    top["interaction"] = (
        top["ligand_complex"] + " → " + top["receptor_complex"]
        + " (" + top["target"] + ")"
    )
    ax.barh(top["interaction"], -np.log10(top["magnitude_rank"] + 1e-6),
            color="steelblue", edgecolor="white")
    ax.set_xlabel("-log10(aggregate rank)")
    ax.set_title(f"Top interactions sent by {source}")
    ax.tick_params(axis="y", labelsize=9)
    ax.invert_yaxis()

plt.tight_layout()
savefig("04_cd8_vs_cd4_sender_profiles")
plt.show()

# %% [markdown]
# ## 8. Myeloid → T Cell Signalling
#
# Monocytes and DCs are key microenvironment regulators of T cell function:
# they present antigen (MHC-II → LAG3/CD4), deliver co-stimulatory signals
# (CD80/CD86 → CD28), or suppress T cells via checkpoint ligands
# (PD-L1 → PD-1, LGALS9 → HAVCR2/TIM3). Here we extract the top
# Monocyte/DC → T cell interactions to characterise myeloid input into
# T cell regulation — the "niche" perspective missing from the CD8-centric
# analysis above.

# %%
myeloid_sources = [c for c in ["Monocyte", "DC", "Monocyte/DC"] if c in liana_res["source"].values]
t_targets = [c for c in ["CD8 T cell", "CD4 T cell"] if c in liana_res["target"].values]

print(f"Myeloid sources available: {myeloid_sources}")
print(f"T cell targets available: {t_targets}")

if myeloid_sources and t_targets:
    myeloid_to_t = (
        liana_res[
            liana_res["source"].isin(myeloid_sources) &
            liana_res["target"].isin(t_targets)
        ]
        .sort_values("magnitude_rank")
        .head(20)
        .copy()
    )
    myeloid_to_t["interaction"] = (
        myeloid_to_t["source"] + ": "
        + myeloid_to_t["ligand_complex"] + " → "
        + myeloid_to_t["receptor_complex"]
        + " (" + myeloid_to_t["target"] + ")"
    )

    fig, ax = plt.subplots(figsize=(10, 7))
    color_map = {"Monocyte": "#E64B35", "DC": "#4DBBD5", "Monocyte/DC": "#E64B35"}
    colors = myeloid_to_t["source"].map(color_map)
    ax.barh(myeloid_to_t["interaction"], -np.log10(myeloid_to_t["magnitude_rank"] + 1e-6),
            color=colors, edgecolor="white")
    ax.set_xlabel("-log10(aggregate rank)")
    ax.set_title("Myeloid → T cell interactions\n(Monocyte/DC as sender, CD8/CD4 T cell as receiver)")
    ax.invert_yaxis()

    from matplotlib.patches import Patch
    ax.legend(
        handles=[Patch(facecolor=c, label=l) for l, c in color_map.items()
                 if l in myeloid_sources],
        loc="lower right",
    )
    plt.tight_layout()
    savefig("04_myeloid_to_tcell_interactions")
    plt.show()

    print("\nTop myeloid → T cell interactions:")
    print(myeloid_to_t[
        ["source", "target", "ligand_complex", "receptor_complex", "magnitude_rank"]
    ].to_string(index=False))
else:
    print("Insufficient myeloid or T cell types in dataset for this comparison.")

# %% [markdown]
# ## Summary
#
# Key findings from the cell-cell interaction analysis:
#
# **Expanded vs non-expanded CD8 T cells (Section 6):**
# - CCL5 -> SDC1 (CD8 T -> Plasma cell) is ranked ~8x stronger in expanded
#   clones (0.010 vs 0.080); groups are equal in size (~25k cells each) so
#   this is not a cell-count artefact — expanded clones communicate more
#   strongly with plasma cells, though the functional consequence requires
#   experimental follow-up (hypothesis-generating finding)
#
# **CD8 vs CD4 sender profiles (Section 7):**
# - LIANA recovers the known functional division from expression alone:
#   CD8 T cells dominate MHC-I antigen presentation interactions (HLA ->
#   CD8A/CD8B) while CD4 T cells send broader co-stimulatory and
#   cytokine-type signals — a sanity check confirming the interaction
#   scores are biologically meaningful
#
# **Monocyte/DC -> T cell signalling (Section 8):**
# - Monocyte/DCs play a "double agent" role in the tumour microenvironment:
#   they simultaneously (1) present antigen to CD8 T cells (HLA-A/C ->
#   CD8A/CD8B, the strongest signal), (2) recruit T cells into the tumour
#   via HMGB1 and MIF -> CXCR4 chemotaxis, and (3) suppress T cells via
#   LGALS1 (galectin-1) -> CD69/PTPRC — the same myeloid population both
#   activates and dampens the T cell response
# - This activation/suppression tension from myeloid cells mirrors the
#   LAG3/CCL5 finding from the T cell side: both point to chronic antigen
#   stimulation driving a response that is simultaneously active and
#   increasingly restrained — candidate pathways for microenvironment-
#   targeted therapeutic intervention

# %%
print("Notebook 4 complete. All figures saved to figures/")
