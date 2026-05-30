# scRNA-seq + TCR Integration Analysis

Single-cell RNA-seq and T cell receptor (TCR) sequencing analysis pipeline for characterising tumour-infiltrating T cell states, clonal dynamics, and cell-cell interactions in the cancer microenvironment.

## Overview

| Notebook | Analysis |
|---|---|
| `01_preprocessing.py` | QC filtering, normalization, PCA, UMAP, Leiden clustering |
| `02_cell_annotation.py` | Marker-based cell type annotation, T cell subset scoring |
| `03_tcr_analysis.py` | Clonotype integration, clonal expansion, **clonotype → transcriptional state linkage** |
| `04_cell_interaction.py` | Ligand-receptor inference (LIANA), niche analysis, expanded vs non-expanded T cell interactions |

## Dataset

**Wu et al. 2020 (Nature)** — peripheral and tumour-infiltrating T cells from patients with basal cell carcinoma. Contains paired scRNA-seq and TCR (alpha + beta chain) sequencing.

Downloaded automatically via `scirpy.datasets.wu2020()` — no registration required.

> Wu TD, Madireddi S, de Almeida PE, et al. Peripheral T cell expansion predicts tumour infiltration and clinical response. *Nature.* 2020;579(7798):274-278.

## Setup

```bash
conda env create -f environment.yml
conda activate scrna-tcr
python -m ipykernel install --user --name scrna-tcr
```

## Running the Analysis

```bash
cd notebooks/

# Run sequentially
python 01_preprocessing.py
python 02_cell_annotation.py
python 03_tcr_analysis.py
python 04_cell_interaction.py
```

Or convert to Jupyter notebooks with jupytext:

```bash
pip install jupytext
jupytext --to notebook notebooks/01_preprocessing.py
jupyter notebook
```

## Results

### Visualizations

**UMAP: Cell type annotation**
![Cell type UMAP](figures/02_umap_celltypes.png)

**UMAP: Clonal expansion**
![Clonal expansion UMAP](figures/03_umap_clonal_expansion.png)

**Clonal expansion by cell type**
![Expansion by cell type](figures/03_expansion_by_celltype.png)

**Transcriptional state by clonal expansion** *(core analysis)*
![State by expansion](figures/03_state_by_expansion.png)

**Top 20 expanded clonotypes: transcriptional state heatmap**
![Clonotype state heatmap](figures/03_clonotype_state_heatmap.png)

**Cell-cell interaction network**
![Interaction network](figures/04_interaction_network.png)

**CD8 T cell ligand-receptor interactions**
![CD8 interactions](figures/04_cd8_interactions_dotplot.png)

---

### Dataset Summary

| Step | Metric | Value |
|------|--------|-------|
| QC filtering | Cells retained | **138,500** |
| QC filtering | Genes retained | **20,479** |
| Highly variable genes | Selected for PCA | **3,000** |
| Cell annotation | T cells identified | **96,015** (69% of total) |
| TCR integration | Cells with paired TCR | **95,314** |
| Clonotype definition | Unique clonotypes | **53,142** |
| LIANA | Ligand-receptor pairs scored | **4,517** |

### Clonal Expansion

| Category | Cells | Description |
|----------|-------|-------------|
| Singleton | 44,937 | Clone size = 1 |
| Small | 17,155 | Clone size 2–5 |
| Medium | 12,950 | Clone size 6–20 |
| **Large** | **20,272** | Clone size > 20 — likely tumour-reactive |

Large expanded clones represent **21% of all TCR+ cells**, consistent with strong antigen-driven expansion in the tumour microenvironment.

### Key Finding: Clonotype → Transcriptional State

Expanded clones (medium + large) show elevated exhaustion scores relative to singleton clones, consistent with progressive acquisition of an exhausted transcriptional state in tumour-reactive T cells. Visualised in `figures/03_clonotype_state_heatmap.png` and `figures/03_state_by_expansion.png`.

---

## Key Analysis: Clonotype → Transcriptional State

The central analysis in Notebook 3 links T cell clonal identity with functional transcriptional programme — addressing the question of whether expanded, tumour-reactive clones preferentially occupy exhausted, effector, or naive states. This directly mirrors the analytical focus of labs studying T cell biology in leukemia and bone marrow microenvironments.

## Tools

| Tool | Purpose |
|---|---|
| [Scanpy](https://scanpy.readthedocs.io) | scRNA-seq preprocessing, clustering, visualization |
| [scirpy](https://scirpy.scverse.org) | TCR data integration, clonotype analysis, VDJ gene usage |
| [LIANA](https://liana-py.readthedocs.io) | Ligand-receptor cell-cell interaction inference |
| [muon](https://muon.readthedocs.io) | Multi-modal data handling |
