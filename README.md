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

## Key Analysis: Clonotype → Transcriptional State

The central analysis in Notebook 3 links T cell clonal identity with functional transcriptional programme — addressing the question of whether expanded, tumour-reactive clones preferentially occupy exhausted, effector, or naive states. This directly mirrors the analytical focus of labs studying T cell biology in leukemia and bone marrow microenvironments.

## Tools

| Tool | Purpose |
|---|---|
| [Scanpy](https://scanpy.readthedocs.io) | scRNA-seq preprocessing, clustering, visualization |
| [scirpy](https://scirpy.scverse.org) | TCR data integration, clonotype analysis, VDJ gene usage |
| [LIANA](https://liana-py.readthedocs.io) | Ligand-receptor cell-cell interaction inference |
| [muon](https://muon.readthedocs.io) | Multi-modal data handling |
