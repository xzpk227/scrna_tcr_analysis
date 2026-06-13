# %% [markdown]
# # Notebook 6: Reproducing Figure 2 (Wu et al. 2020)
#
# Figure 2 of the original paper links **T cell transcriptional clusters**
# to the tissue/blood expansion patterns defined in Figure 1 (Notebook 5).
# The key finding is that clonotypes which are **dual-expanded** (tumour +
# NAT, pattern `D`) AND **peripherally expanded** (detected as a blood
# multiplet, pattern `B`) are enriched for an **effector-like CD8 T cell
# phenotype**, in contrast to dual-expanded clones that are *not* detected
# as blood-expanded, which skew towards more exhausted/naive-like states.
#
# We don't have the paper's full 33-cluster / 16-subtype annotation, but our
# existing Leiden clustering (Notebook 1/2, resolution 0.5) already splits
# CD8 and CD4 T cells into multiple subclusters, and Notebook 2 computed
# per-cell `Exhausted_score`, `Effector_score`, and `Naive_score` gene-set
# scores. We use these to characterise each cluster and test the same
# enrichment claim.
#
# Approach:
# 1. Load `data/03_tcr_integrated.h5ad` and the clone classification table
#    from Notebook 5 (`data/05_figure1_clone_classification.csv`).
# 2. Restrict to **dual-expanded (`D`)** clones in the 4 blood-sample
#    patients, split into:
#    - **D & B** — dual-expanded AND peripherally (blood-)expanded
#    - **D & not-B** — dual-expanded but blood singleton/undetected
# 3. Characterise CD8/CD4 Leiden clusters by mean Exhausted/Effector/Naive
#    scores.
# 4. Compare cluster composition and gene-set scores between D&B and D&not-B
#    cells.

# %%
import sys
sys.path.insert(0, "..")

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency, mannwhitneyu

from src.utils import savefig

sc.settings.verbosity = 1

# %% [markdown]
# ## 1. Load Data and Merge Clone Classification

# %%
adata = sc.read("../data/03_tcr_integrated.h5ad")
clone_class = pd.read_csv("../data/05_figure1_clone_classification.csv")
print(f"Loaded: {adata.n_obs:,} cells")
print(f"Clone classification table: {len(clone_class):,} clonotypes "
      f"(patients: {clone_class['patient'].unique().tolist()})")

# %%
obs = adata.obs.reset_index().rename(columns={"index": "barcode"})
obs["clone_id_int"] = obs["clone_id"].astype("Int64")

merged = obs.merge(
    clone_class[["patient", "clone_id", "tissue_pattern", "blood_pattern"]],
    left_on=["patient", "clone_id_int"], right_on=["patient", "clone_id"],
    how="inner",
).set_index("barcode")

print(f"Cells matched to a classified clonotype: {len(merged):,}")

# %% [markdown]
# ## 2. Define Comparison Groups
#
# Focus on **dual-expanded (`D`)** clones — those detected in both tumour
# and NAT — and split by whether the clone is also detected as
# blood-expanded (`B`, peripherally expanded) or not (`b`/`none`).

# %%
D_cells = merged[merged["tissue_pattern"] == "D"].copy()
D_cells["periph_group"] = np.where(
    D_cells["blood_pattern"] == "B", "D & B (periph. expanded)", "D & not-B"
)

print(D_cells.groupby(["patient", "periph_group"], observed=True).size().unstack(fill_value=0))
print()
print(D_cells["periph_group"].value_counts())

# %% [markdown]
# ## 3. Subcluster CD8 T Cells at Higher Resolution
#
# The dataset-wide Leiden clustering (resolution 0.5, used in Notebook 1)
# only yields **4 CD8 T cell clusters** — too coarse to separate an
# "effector" subtype from "exhausted"/"naive" subtypes the way the paper's
# 33-cluster analysis does (where CD8 alone splits into 6 subtypes:
# 8.1-Teff, 8.2-Tem, 8.3a/b/c-Trm, 8.4-Chrom, 8.5-Mitosis, 8.6-KLRB1).
#
# `adata.X` here holds **raw counts** (transferred from the original
# `wu2020()` MuData in Notebook 3), so we re-normalize, re-select HVGs, and
# re-cluster the CD8 T cell subset alone at a higher resolution to better
# resolve these subtypes.
#
# Pooling all 14 patients' CD8 cells without any batch correction risks
# Leiden clusters that reflect patient-specific variation rather than shared
# cell states (as the original paper avoids via Seurat's
# `FindIntegrationAnchors`/`IntegrateData`). We use Harmony
# (`sc.external.pp.harmony_integrate`, batch key = `patient`) to integrate
# the PCA embedding across patients before computing neighbors/Leiden.

# %%
cd8_adata = adata[adata.obs["cell_type"] == "CD8 T cell"].copy()
print(f"CD8 T cells: {cd8_adata.n_obs:,}")

sc.pp.normalize_total(cd8_adata, target_sum=1e4)
sc.pp.log1p(cd8_adata)
sc.pp.highly_variable_genes(cd8_adata, n_top_genes=2000)

cd8_hvg = cd8_adata[:, cd8_adata.var["highly_variable"]].copy()
sc.pp.scale(cd8_hvg, max_value=10)
sc.tl.pca(cd8_hvg, svd_solver="arpack", n_comps=30)
sc.external.pp.harmony_integrate(cd8_hvg, key="patient")
sc.pp.neighbors(cd8_hvg, n_neighbors=15, n_pcs=30, use_rep="X_pca_harmony")
sc.tl.leiden(cd8_hvg, resolution=1.0, key_added="cd8_subcluster")

cd8_adata.obs["cd8_subcluster"] = cd8_hvg.obs["cd8_subcluster"]
print(f"CD8 subclusters: {cd8_adata.obs['cd8_subcluster'].nunique()}")
print(cd8_adata.obs["cd8_subcluster"].value_counts().sort_index())

# %% [markdown]
# ## 4. Characterise CD8 Subclusters
#
# Mean Exhausted/Effector/Naive gene-set scores per CD8 subcluster. The
# subcluster with the **highest Effector score combined with the lowest
# Exhausted score** is labelled "Effector-like" (analogous to the paper's
# 8.1-Teff); the rest "Other".

# %%
cluster_scores = cd8_adata.obs.groupby("cd8_subcluster", observed=True)[
    ["Exhausted_score", "Effector_score", "Naive_score"]
].mean()
print("CD8 subcluster mean scores:")
print(cluster_scores.round(3))

# Effector-like: high Effector_score, low Exhausted_score -> rank by
# (Effector_score - Exhausted_score), descending
effector_rank = (cluster_scores["Effector_score"] - cluster_scores["Exhausted_score"])
effector_cluster = effector_rank.idxmax()
cd8_cluster_labels = {
    cl: ("Effector-like" if cl == effector_cluster else "Other")
    for cl in cluster_scores.index
}
print(f"\nEffector-like subcluster: {effector_cluster}")
print(f"CD8 subcluster labels: {cd8_cluster_labels}")

# Merge subcluster assignment back onto the cell-level table
D_cells = D_cells.join(cd8_adata.obs[["cd8_subcluster"]], how="left")

# %% [markdown]
# ## 5. Cluster Composition: D&B vs D&not-B (CD8 T cells)
#
# For CD8 T cells in dual-expanded clones, compare the CD8 subcluster
# distribution between peripherally-expanded (`D & B`) and
# non-peripherally-expanded (`D & not-B`) groups. A chi-square test assesses
# whether the distributions differ, with the paper's prediction being
# enrichment of the "Effector-like" subcluster in `D & B`.

# %%
D_cd8 = D_cells[D_cells["cell_type"] == "CD8 T cell"].copy()
D_cd8["cluster_label"] = D_cd8["cd8_subcluster"].map(cd8_cluster_labels)

comp = (
    D_cd8.groupby(["periph_group", "cluster_label"], observed=True)
    .size()
    .unstack(fill_value=0)
)
comp_frac = comp.div(comp.sum(axis=1), axis=0)
print("Cell counts:")
print(comp)
print("\nFractions:")
print(comp_frac.round(3))

chi2, p, dof, _ = chi2_contingency(comp)
print(f"\nChi-square test: chi2={chi2:.2f}, dof={dof}, p={p:.2e}")

# %%
fig, ax = plt.subplots(figsize=(6, 5))
comp_frac.plot(kind="bar", stacked=True, ax=ax, edgecolor="white",
                color=["#E41A1C", "#377EB8"])
ax.set_ylabel("Fraction of CD8 T cells")
ax.set_xlabel("")
ax.set_title(
    "CD8 cluster composition of dual-expanded (D) clones,\n"
    f"by peripheral expansion status (chi2 p={p:.1e})"
)
ax.tick_params(axis="x", rotation=15)
ax.legend(title="Cluster type", bbox_to_anchor=(1.02, 1), loc="upper left")
plt.tight_layout()
savefig("06_cd8_cluster_composition")
plt.show()

# %% [markdown]
# ## 6. Gene-Set Scores: D&B vs D&not-B
#
# Direct comparison of Exhausted/Effector/Naive scores between the two
# groups, for CD8 T cells in dual-expanded clones — a score-based version of
# the same claim that doesn't depend on cluster labels.

# %%
fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
group_order = ["D & B (periph. expanded)", "D & not-B"]

for ax, score in zip(axes, ["Exhausted_score", "Effector_score", "Naive_score"]):
    data_plot = [
        D_cd8.loc[D_cd8["periph_group"] == g, score].dropna().values
        for g in group_order
    ]
    ax.boxplot(data_plot, labels=group_order, patch_artist=True,
               boxprops=dict(facecolor="lightblue"),
               medianprops=dict(color="red", linewidth=2))
    ax.set_title(score.replace("_score", ""))
    ax.tick_params(axis="x", rotation=15)

    u, p = mannwhitneyu(*data_plot, alternative="two-sided")
    ax.text(0.5, 0.98, f"Mann-Whitney p={p:.1e}", transform=ax.transAxes,
            ha="center", va="top", fontsize=8)

plt.suptitle("CD8 T cell transcriptional state in dual-expanded (D) clones,\n"
              "by peripheral expansion status", fontsize=12)
plt.tight_layout()
savefig("06_cd8_scores_by_periph_expansion")
plt.show()

# %% [markdown]
# ## 7. Per-Clone Primary Cluster
#
# Following the paper's approach, assign each clone a **primary cluster**
# (the CD8 subcluster with the most cells in that clone, among its CD8 T
# cells), then compare primary-cluster assignment between D&B and D&not-B
# clones (one row per clone, rather than per cell).

# %%
def primary_cluster(group):
    counts = group["cd8_subcluster"].value_counts()
    if len(counts) >= 2 and counts.iloc[0] == counts.iloc[1]:
        return np.nan  # tie -> no primary cluster
    return counts.idxmax()


clone_primary = (
    D_cd8.groupby(["patient", "clone_id_int", "periph_group"], observed=True)
    .apply(primary_cluster, include_groups=False)
    .reset_index(name="primary_cluster")
    .dropna(subset=["primary_cluster"])
)
clone_primary["primary_label"] = clone_primary["primary_cluster"].map(cd8_cluster_labels)

clone_comp = (
    clone_primary.groupby(["periph_group", "primary_label"], observed=True)
    .size()
    .unstack(fill_value=0)
)
clone_comp_frac = clone_comp.div(clone_comp.sum(axis=1), axis=0)
print("Clone counts by primary cluster label:")
print(clone_comp)
print("\nFractions:")
print(clone_comp_frac.round(3))

chi2_c, p_c, dof_c, _ = chi2_contingency(clone_comp)
print(f"\nChi-square test (clone-level): chi2={chi2_c:.2f}, dof={dof_c}, p={p_c:.2e}")

# %% [markdown]
# ## 8. Panel d Reproduction: D-Fraction by Blood Status, within Effector-like Clones
#
# Sections 5-7 test the paper's claim with the conditioning reversed (start
# from dual-expanded (`D`) clones, ask if blood-expansion enriches for the
# Effector-like subcluster). The paper's actual Fig 2d conditions the other
# way: starting from clones whose **primary cluster is 8.1-Teff**, it asks
# what fraction have tissue_pattern `D`, comparing:
#
# - **Ind** — blood-independent (`blood_pattern == "none"`, no blood cells)
# - **Non** — blood non-expanded (`blood_pattern == "b"`, blood singleton)
# - **Exp** — blood-expanded (`blood_pattern == "B"`, blood multiplet)
#
# with the paper finding D much more common among Exp clones. We reproduce
# this directly using our "Effector-like" CD8 subcluster as the analogue of
# 8.1-Teff, over **all** clones from the 4 blood-sample patients (not just
# the `D` ones).

# %%
cd8_merged = merged[merged["cell_type"] == "CD8 T cell"].join(
    cd8_adata.obs[["cd8_subcluster"]], how="left"
)

all_primary = (
    cd8_merged.groupby(["patient", "clone_id_int"], observed=True)
    .apply(primary_cluster, include_groups=False)
    .reset_index(name="primary_cluster")
    .dropna(subset=["primary_cluster"])
)
all_primary["primary_label"] = all_primary["primary_cluster"].map(cd8_cluster_labels)

clone_info = clone_class[["patient", "clone_id", "tissue_pattern", "blood_pattern"]].rename(
    columns={"clone_id": "clone_id_int"}
)
all_primary = all_primary.merge(clone_info, on=["patient", "clone_id_int"], how="left")

BLOOD_LABELS = {
    "none": "Ind (blood-independent)",
    "b": "Non (blood non-expanded)",
    "B": "Exp (blood-expanded)",
}
all_primary["blood_group"] = all_primary["blood_pattern"].map(BLOOD_LABELS)
all_primary["is_D"] = all_primary["tissue_pattern"] == "D"

effector_clones = all_primary[all_primary["primary_label"] == "Effector-like"].copy()
print(f"Effector-like (primary cluster) clones: {len(effector_clones):,} / {len(all_primary):,}")

BLOOD_ORDER = ["Ind (blood-independent)", "Non (blood non-expanded)", "Exp (blood-expanded)"]
panel_d = (
    effector_clones.groupby("blood_group", observed=True)["is_D"]
    .agg(n_D="sum", n_total="count")
    .reindex(BLOOD_ORDER)
)
panel_d["frac_D"] = panel_d["n_D"] / panel_d["n_total"]
print(panel_d)

# 3-group chi-square
ct = pd.crosstab(effector_clones["blood_group"], effector_clones["is_D"]).reindex(BLOOD_ORDER)
chi2_d, p_d, dof_d, _ = chi2_contingency(ct)
print(f"\nChi-square (Ind vs Non vs Exp): chi2={chi2_d:.2f}, dof={dof_d}, p={p_d:.2e}")

# Paper-style pairwise comparison: (Ind + Non) vs Exp
ind_non = ct.loc[["Ind (blood-independent)", "Non (blood non-expanded)"]].sum()
exp = ct.loc["Exp (blood-expanded)"]
ct_pair = pd.DataFrame([ind_non, exp], index=["Ind+Non", "Exp"])
chi2_p2, p_p2, dof_p2, _ = chi2_contingency(ct_pair)
print(f"Chi-square (Ind+Non vs Exp): chi2={chi2_p2:.2f}, dof={dof_p2}, p={p_p2:.2e}")

# %%
fig, ax = plt.subplots(figsize=(5, 4.5))
panel_d["frac_D"].plot(kind="bar", ax=ax, color="#984EA3", edgecolor="black")
for i, (n_d, n_t) in enumerate(zip(panel_d["n_D"], panel_d["n_total"])):
    ax.text(i, panel_d["frac_D"].iloc[i] + 0.01, f"n={int(n_d)}/{int(n_t)}", ha="center")
ax.set_ylabel("Fraction of Effector-like clones with tissue_pattern == D")
ax.set_title(
    "Panel-d-style reproduction: D-fraction of Effector-like clones,\n"
    f"by blood status (3-group chi2 p={p_d:.1e})"
)
ax.tick_params(axis="x", rotation=15)
plt.tight_layout()
savefig("06_paneld_dfraction_by_blood")
plt.show()

# %% [markdown]
# ## 8b. Panel d Reproduction (v2): Clone Size of D-Pattern 8.1-Teff Clones, by Blood Status
#
# Looking at the actual Fig 2d more closely, it is **not** a bar chart of
# D-fraction — it's a swarm plot. For clones whose primary cluster is
# **8.1-Teff** AND `tissue_pattern == D` (dual-expanded), each dot is one
# clone, x = clone size (tumour + NAT cells), split into Ind / Non / Exp
# blood-status columns. `n` = number of such D-pattern 8.1-Teff clones per
# group (paper: 16 / 22 / 68), and the P-values compare the **clone-size
# distributions** (Mann-Whitney) between Exp and Ind/Non.
#
# The dataset includes the paper's own per-cell subtype labels in
# `cluster_orig` (16 categories, including `8.1-Teff` directly) — so we use
# those instead of our Harmony-derived "Effector-like" proxy from Section 8,
# no re-clustering needed.

# %%
cd8_cells_orig = merged[merged["cell_type"] == "CD8 T cell"]


def primary_cluster_orig(group):
    counts = group["cluster_orig"].value_counts()
    if len(counts) >= 2 and counts.iloc[0] == counts.iloc[1]:
        return np.nan  # tie -> no primary cluster
    return counts.idxmax()


primary_orig = (
    cd8_cells_orig.groupby(["patient", "clone_id_int"], observed=True)
    .apply(primary_cluster_orig, include_groups=False)
    .reset_index(name="primary_cluster_orig")
    .dropna(subset=["primary_cluster_orig"])
)

clone_sizes = clone_class[["patient", "clone_id", "Tumor", "NAT", "tissue_pattern", "blood_pattern"]].rename(
    columns={"clone_id": "clone_id_int"}
)
primary_orig = primary_orig.merge(clone_sizes, on=["patient", "clone_id_int"], how="left")
primary_orig["clone_size_tumor_nat"] = primary_orig["Tumor"] + primary_orig["NAT"]
primary_orig["blood_group"] = primary_orig["blood_pattern"].map(BLOOD_LABELS)

print(f"Clones with a primary cluster_orig label: {len(primary_orig):,}")
print("Primary cluster_orig value counts (top 5):")
print(primary_orig["primary_cluster_orig"].value_counts().head())

teff_D = primary_orig[
    (primary_orig["primary_cluster_orig"] == "8.1-Teff")
    & (primary_orig["tissue_pattern"] == "D")
].copy()

print(f"\n8.1-Teff (primary cluster_orig) clones with tissue_pattern == D: {len(teff_D)}")
size_summary = (
    teff_D.groupby("blood_group", observed=True)["clone_size_tumor_nat"]
    .agg(n="count", median="median", mean="mean")
    .reindex(BLOOD_ORDER)
)
print(size_summary)

exp_sizes = teff_D.loc[teff_D["blood_group"] == "Exp (blood-expanded)", "clone_size_tumor_nat"]
ind_sizes = teff_D.loc[teff_D["blood_group"] == "Ind (blood-independent)", "clone_size_tumor_nat"]
non_sizes = teff_D.loc[teff_D["blood_group"] == "Non (blood non-expanded)", "clone_size_tumor_nat"]

if len(ind_sizes) > 0 and len(exp_sizes) > 0:
    _, p_ie = mannwhitneyu(ind_sizes, exp_sizes, alternative="two-sided")
    print(f"\nMann-Whitney Ind vs Exp: p={p_ie:.2e}")
if len(non_sizes) > 0 and len(exp_sizes) > 0:
    _, p_ne = mannwhitneyu(non_sizes, exp_sizes, alternative="two-sided")
    print(f"Mann-Whitney Non vs Exp: p={p_ne:.2e}")

# %%
rng = np.random.default_rng(0)
fig, axes = plt.subplots(1, 3, figsize=(9, 4), sharey=True, sharex=True)
for ax, grp in zip(axes, BLOOD_ORDER):
    sub = teff_D[teff_D["blood_group"] == grp]
    y = rng.uniform(-0.3, 0.3, size=len(sub))
    ax.scatter(sub["clone_size_tumor_nat"].clip(lower=0.5), y, alpha=0.6,
               color="#FF7F00", edgecolor="black")
    ax.set_title(f"{grp.split(' ')[0]} (n={len(sub)})")
    ax.set_xscale("symlog")
    ax.set_yticks([])
fig.supxlabel("Clone size (tumour + NAT)")
fig.suptitle("Panel-d reproduction (v2): 8.1-Teff (cluster_orig), tissue_pattern == D,\nby blood status")
plt.tight_layout()
savefig("06_paneld_v2_clonesize_swarm")
plt.show()

# %% [markdown]
# ## 9. Save

# %%
D_cells.to_csv("../data/06_dual_expanded_cells.csv", index=True)
clone_primary.to_csv("../data/06_clone_primary_clusters.csv", index=False)
all_primary.to_csv("../data/06_all_clones_primary_clusters.csv", index=False)
primary_orig.to_csv("../data/06_clone_primary_cluster_orig.csv", index=False)
print("Saved: data/06_dual_expanded_cells.csv, data/06_clone_primary_clusters.csv, "
      "data/06_all_clones_primary_clusters.csv, data/06_clone_primary_cluster_orig.csv")
