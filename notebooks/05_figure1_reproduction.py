# %% [markdown]
# # Notebook 5: Reproducing Figure 1 (Wu et al. 2020)
#
# Figure 1 of the original paper ("Parallel dual expansion and peripheral
# clonal expansion") characterises clonotypes using a **per-tissue-compartment**
# expansion scheme, separate from the overall singleton/small/medium/large
# categorisation used in Notebook 3.
#
# Only 4 of the 14 patients have matched **Blood** samples: Renal1, Renal2,
# Renal3, and Lung6. For these patients, each clonotype is classified by:
#
# - **Tissue expansion pattern** (based on clone size in NAT vs Tumor, per
#   the "Clone size Tu/NAT/Bl" legend in the paper's Fig 1b):
#   - `n` — NAT singleton (NAT == 1, Tumor == 0)
#   - `N` — NAT multiplet (NAT > 1, Tumor == 0)
#   - `t` — Tumour singleton (Tumor == 1, NAT == 0)
#   - `T` — Tumour multiplet (Tumor > 1, NAT == 0)
#   - `D` — Dual-expanded (detected, i.e. >=1 cell, in BOTH NAT and tumour —
#     this matches the paper's `n_D` counts to within ~4-10% per patient,
#     vs. ~3x undercount when D required *multiplet* in both compartments)
#   - `x` — not detected in NAT or tumour (blood-only clone)
# - **Blood expansion pattern** (based on clone size in Blood):
#   - `b` — blood singleton (1 cell)
#   - `B` — blood multiplet (>1 cell, "peripherally expanded")
#   - `none` — not detected in blood
#
# This notebook reproduces the key panels of Figure 1:
# - **(a)** Per-patient scatter of clonotype frequency in tumour vs NAT,
#   coloured by tissue expansion pattern
# - **(c)** Peripheral (blood) clonal expansion per patient
# - **(d)** Tissue infiltration pattern, split by blood expansion
# - **(e)** Detection of tissue-resident TCRs in blood, by tissue pattern
# - **(f)** Correlation of blood vs tumour clone size

# %%
import sys
sys.path.insert(0, "..")

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from src.utils import savefig

sc.settings.verbosity = 1

# %% [markdown]
# ## 1. Load Data and Subset to Blood-Sample Patients

# %%
adata = sc.read("../data/03_tcr_integrated.h5ad")
print(f"Loaded: {adata.n_obs:,} cells")

BLOOD_PATIENTS = ["Renal1", "Renal2", "Renal3", "Lung6"]
TISSUE_PATTERN_ORDER = ["D", "N", "T", "n", "t", "x"]
BLOOD_PATTERN_ORDER = ["B", "b", "none"]

TISSUE_COLORS = {
    "D": "#984EA3",   # dual-expanded — purple
    "N": "#377EB8",   # NAT multiplet — blue
    "T": "#E41A1C",   # tumour multiplet — red
    "n": "#A6CEE3",   # NAT singleton — light blue
    "t": "#FB9A99",   # tumour singleton — light red
    "x": "#CCCCCC",   # not in NAT/tumour — grey
}
BLOOD_COLORS = {"B": "#4DAF4A", "b": "#B3DE69", "none": "#CCCCCC"}

obs = adata.obs[adata.obs["patient"].isin(BLOOD_PATIENTS)].copy()
obs = obs[obs["clone_id"].notna()].copy()
print(f"Cells with TCR in blood-sample patients: {obs.shape[0]:,}")
print(obs.groupby(["patient", "source"], observed=True).size().unstack(fill_value=0))

# %% [markdown]
# ## 2. Classify Clonotypes by Tissue and Blood Expansion Pattern
#
# For each patient, count how many cells of each clonotype fall in each
# tissue source (Tumor / NAT / Blood), then assign the expansion patterns
# defined above.

# %%
# Per-patient, per-clone cell counts by source
clone_counts = (
    obs.groupby(["patient", "clone_id", "source"], observed=True)
    .size()
    .unstack(fill_value=0)
    .reindex(columns=["Tumor", "NAT", "Blood"], fill_value=0)
    .reset_index()
)


def tissue_pattern(row):
    # D: detected (>=1 cell) in BOTH NAT and tumour, per Fig 1b's "Tu>=1, NAT>=1" rule
    if row["NAT"] >= 1 and row["Tumor"] >= 1:
        return "D"
    if row["NAT"] > 1:
        return "N"
    if row["Tumor"] > 1:
        return "T"
    if row["NAT"] == 1:
        return "n"
    if row["Tumor"] == 1:
        return "t"
    return "x"


def blood_pattern(row):
    if row["Blood"] > 1:
        return "B"
    if row["Blood"] == 1:
        return "b"
    return "none"


clone_counts["tissue_pattern"] = clone_counts.apply(tissue_pattern, axis=1)
clone_counts["blood_pattern"] = clone_counts.apply(blood_pattern, axis=1)
clone_counts["tissue_pattern"] = pd.Categorical(
    clone_counts["tissue_pattern"], categories=TISSUE_PATTERN_ORDER, ordered=True
)
clone_counts["blood_pattern"] = pd.Categorical(
    clone_counts["blood_pattern"], categories=BLOOD_PATTERN_ORDER, ordered=True
)

print(f"Total clonotypes classified: {len(clone_counts):,}")
print("\nTissue expansion pattern counts:")
print(clone_counts["tissue_pattern"].value_counts())
print("\nBlood expansion pattern counts:")
print(clone_counts["blood_pattern"].value_counts())

# %% [markdown]
# ## 3. Panel (a) — Clonotype Frequency: NAT vs Tumour
#
# Each point is a clonotype. Axes show its frequency (fraction of all TCR+
# cells from that compartment, for that patient) in NAT vs tumour — matching
# the axis orientation of the paper's Fig 1a. Colour = tissue expansion
# pattern. Point size scales with blood clone size.
#
# `n_D`, `r` (log-clone-size correlation for clones shared between blood and
# tumour) and the dual-expanded blood-detection fraction are annotated per
# patient for direct comparison with the paper's Fig 1a/e values.

# %%
PAPER_FIG1 = {
    "Renal2": {"n_D": 139, "r": 0.93, "D_detect": 0.87},
    "Renal1": {"n_D": 108, "r": 0.81, "D_detect": 0.83},
    "Renal3": {"n_D": 119, "r": 0.58, "D_detect": 0.63},
    "Lung6":  {"n_D": 213, "r": 0.31, "D_detect": 0.05},
}
PANEL_A_ORDER = ["Renal2", "Renal1", "Renal3", "Lung6"]

fig, axes = plt.subplots(2, 2, figsize=(11, 10))

for ax, patient in zip(axes.flat, PANEL_A_ORDER):
    sub = clone_counts[clone_counts["patient"] == patient].copy()
    tumor_total = sub["Tumor"].sum()
    nat_total = sub["NAT"].sum()
    sub["tumor_frac"] = sub["Tumor"] / tumor_total if tumor_total else 0
    sub["nat_frac"] = sub["NAT"] / nat_total if nat_total else 0

    for pattern in TISSUE_PATTERN_ORDER:
        pts = sub[sub["tissue_pattern"] == pattern]
        if pts.empty:
            continue
        sizes = 15 + 25 * np.clip(pts["Blood"], 0, 10)
        ax.scatter(
            pts["nat_frac"], pts["tumor_frac"],
            s=sizes, c=TISSUE_COLORS[pattern], label=pattern,
            alpha=0.6, edgecolors="none",
        )

    n_D = (sub["tissue_pattern"] == "D").sum()
    D_sub = sub[sub["tissue_pattern"] == "D"].copy()
    D_sub["tissue_cells"] = D_sub["Tumor"] + D_sub["NAT"]
    D_total_cells = D_sub["tissue_cells"].sum()
    D_detect = (
        D_sub.loc[D_sub["Blood"] > 1, "tissue_cells"].sum() / D_total_cells
        if D_total_cells else float("nan")
    )
    both = sub[(sub["Blood"] > 0) & (sub["Tumor"] > 0)]
    r = np.corrcoef(np.log1p(both["Tumor"]), np.log1p(both["Blood"]))[0, 1] if len(both) > 1 else float("nan")
    paper = PAPER_FIG1[patient]

    ax.set_xlabel("NAT clone frequency")
    ax.set_ylabel("Tumour clone frequency")
    ax.set_title(
        f"{patient}\n"
        f"n_D={n_D} (paper {paper['n_D']}), r={r:.2f} (paper {paper['r']:.2f}),\n"
        f"D-detect={D_detect:.2f} (paper {paper['D_detect']:.2f})"
    )
    ax.set_xscale("symlog", linthresh=1e-3)
    ax.set_yscale("symlog", linthresh=1e-3)

axes[0, 0].legend(
    title="Tissue pattern\n(point size ~ blood count)",
    bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8,
)
plt.suptitle(
    "Clonotype frequency in NAT vs tumour (blood-sample patients only)\n"
    "annotations compare against Wu et al. 2020 Fig 1a/e",
    fontsize=13,
)
plt.tight_layout()
savefig("05_clonotype_tissue_scatter")
plt.show()

# %% [markdown]
# ## 4. Panel (c) — Peripheral (Blood) Clonal Expansion
#
# The paper computes this as the **fraction of blood cells** (not
# clonotypes) belonging to blood-expanded (`B`, multiplet) vs
# blood-non-expanded (`b`, singleton) clonotypes — i.e. cell-weighted,
# so large clones dominate.

# %%
PAPER_EXP_FRAC = {"Renal2": 0.82, "Renal1": 0.74, "Renal3": 0.50, "Lung6": 0.04}

peripheral_cells = []
for patient in BLOOD_PATIENTS:
    sub = clone_counts[clone_counts["patient"] == patient]
    blood_total = sub["Blood"].sum()
    b_cells = sub.loc[sub["blood_pattern"] == "b", "Blood"].sum()
    B_cells = sub.loc[sub["blood_pattern"] == "B", "Blood"].sum()
    peripheral_cells.append({"patient": patient, "b": b_cells, "B": B_cells})
peripheral_cells = pd.DataFrame(peripheral_cells).set_index("patient")
peripheral_frac = peripheral_cells.div(peripheral_cells.sum(axis=1), axis=0)

fig, ax = plt.subplots(figsize=(6, 5))
peripheral_frac.reindex(BLOOD_PATIENTS).plot(
    kind="bar", stacked=True, ax=ax,
    color=[BLOOD_COLORS["b"], BLOOD_COLORS["B"]], edgecolor="white",
)
for i, patient in enumerate(BLOOD_PATIENTS):
    exp = peripheral_frac.loc[patient, "B"]
    ax.text(i, 1.02, f"{exp:.2f}\n(paper {PAPER_EXP_FRAC[patient]:.2f})",
            ha="center", va="bottom", fontsize=8)
ax.set_ylabel("Fraction of blood cells")
ax.set_xlabel("Patient")
ax.set_title("Peripheral clonal expansion\n(blood-expanded vs non-expanded, cell-weighted)")
ax.set_ylim(0, 1.25)
ax.legend(title="Blood pattern", labels=["non-expanded (b)", "expanded (B)"], loc="lower right")
plt.tight_layout()
savefig("05_peripheral_clonal_expansion")
plt.show()

print(peripheral_cells)

# %% [markdown]
# ## 5. Panel (d) — Tissue Infiltration by Blood Expansion Pattern
#
# For each patient, distribution of **tissue cells** (Tumor + NAT) by tissue
# expansion pattern, split into blood-independent (`Ind`, Bl=0),
# non-expanded (`Non`, Bl=1) and expanded (`Exp`, Bl>1) clonotypes.
# The dual-expanded (`D`) fraction of the `Exp` row is annotated, since the
# paper highlights this value.

# %%
PAPER_D_EXP_FRAC = {"Renal2": 0.96, "Renal1": 0.94, "Renal3": 0.80, "Lung6": 0.52}
ROW_LABELS = {"none": "Ind", "b": "Non", "B": "Exp"}

fig, axes = plt.subplots(1, len(BLOOD_PATIENTS), figsize=(16, 4), sharey=True)
for ax, patient in zip(axes, BLOOD_PATIENTS):
    sub = clone_counts[clone_counts["patient"] == patient].copy()
    sub["tissue_cells"] = sub["Tumor"] + sub["NAT"]
    sub = sub[sub["tissue_cells"] > 0]

    infiltration = (
        sub.groupby(["blood_pattern", "tissue_pattern"], observed=True)["tissue_cells"]
        .sum()
        .unstack(fill_value=0)
        .reindex(index=BLOOD_PATTERN_ORDER, columns=TISSUE_PATTERN_ORDER, fill_value=0)
    )
    infiltration_frac = infiltration.div(infiltration.sum(axis=1), axis=0)
    infiltration_frac.index = [ROW_LABELS[i] for i in infiltration_frac.index]

    infiltration_frac.plot(
        kind="bar", stacked=True, ax=ax, legend=False,
        color=[TISSUE_COLORS[p] for p in TISSUE_PATTERN_ORDER], edgecolor="white",
    )
    d_frac = infiltration_frac.loc["Exp", "D"]
    ax.set_title(f"{patient}\nD (Exp) = {d_frac:.2f} (paper {PAPER_D_EXP_FRAC[patient]:.2f})")
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=0)

axes[0].set_ylabel("Fraction of tissue cells")
axes[-1].legend(title="Tissue pattern", bbox_to_anchor=(1.02, 1), loc="upper left")
plt.suptitle("Tissue infiltration pattern (cell-weighted),\nby blood expansion status")
plt.tight_layout()
savefig("05_tissue_infiltration_by_blood_pattern")
plt.show()

# %% [markdown]
# ## 6. Panel (e) — Detection of Tissue TCRs in Blood
#
# For each tissue expansion category, what fraction of **tissue cells**
# belong to clonotypes also observed as blood-expanded (`B`)?
# The dual-expanded (`D`) value is highlighted for comparison with the paper.

# %%
PAPER_D_DETECT = {"Renal2": 0.87, "Renal1": 0.83, "Renal3": 0.63, "Lung6": 0.05}

fig, axes = plt.subplots(1, len(BLOOD_PATIENTS), figsize=(16, 4), sharey=True)
for ax, patient in zip(axes, BLOOD_PATIENTS):
    sub = clone_counts[clone_counts["patient"] == patient].copy()
    sub["tissue_cells"] = sub["Tumor"] + sub["NAT"]
    sub = sub[sub["tissue_cells"] > 0]

    detection = (
        sub.groupby(["tissue_pattern", "blood_pattern"], observed=True)["tissue_cells"]
        .sum()
        .unstack(fill_value=0)
        .reindex(index=TISSUE_PATTERN_ORDER, columns=BLOOD_PATTERN_ORDER, fill_value=0)
    )
    detection_frac = detection.div(detection.sum(axis=1), axis=0)

    detection_frac[["B", "b", "none"]].plot(
        kind="bar", stacked=True, ax=ax, legend=False,
        color=[BLOOD_COLORS["B"], BLOOD_COLORS["b"], BLOOD_COLORS["none"]],
        edgecolor="white",
    )
    d_detect = detection_frac.loc["D", "B"]
    ax.set_title(f"{patient}\nD detect = {d_detect:.2f} (paper {PAPER_D_DETECT[patient]:.2f})")
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=0)

axes[0].set_ylabel("Fraction of tissue cells")
axes[-1].legend(title="Blood pattern", labels=["expanded (B)", "non-expanded (b)", "not detected"],
                bbox_to_anchor=(1.02, 1), loc="upper left")
plt.suptitle("Detection of tumour/NAT TCRs in blood (cell-weighted)")
plt.tight_layout()
savefig("05_tissue_tcr_detection_in_blood")
plt.show()

# %% [markdown]
# ## 7. Panel (f) — Blood vs Tumour Clone Size Correlation
#
# The paper's panel (f) pools clonotypes shared between blood and tumour
# across all 4 patients (Renal1-3, Lung6) into a single scatter
# (n=328, r=0.72, P=6.3e-53), rather than plotting each patient separately.

# %%
both = clone_counts[(clone_counts["Blood"] > 0) & (clone_counts["Tumor"] > 0)].copy()
print(f"Clonotypes detected in both blood and tumour: {len(both):,} (paper n=328)")

corr = np.corrcoef(np.log1p(both["Tumor"]), np.log1p(both["Blood"]))[0, 1]
print(f"r={corr:.2f} (paper r=0.72)")

fig, ax = plt.subplots(figsize=(6, 5))
for patient in BLOOD_PATIENTS:
    sub = both[both["patient"] == patient]
    ax.scatter(sub["Blood"], sub["Tumor"], alpha=0.5, label=patient, edgecolors="none")

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Blood clone size")
ax.set_ylabel("Tumour clone size")
ax.set_title(
    f"Blood vs tumour clone size (Renal1-3, Lung6)\n"
    f"n={len(both)} (paper 328), r={corr:.2f} (paper 0.72)"
)
ax.legend(title="Patient")
plt.tight_layout()
savefig("05_blood_tumor_clone_correlation")
plt.show()

# %% [markdown]
# ## 8. Save Classified Clonotypes

# %%
clone_counts.to_csv("../data/05_figure1_clone_classification.csv", index=False)
print("Saved: data/05_figure1_clone_classification.csv")
print(clone_counts.head())
