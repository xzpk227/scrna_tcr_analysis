# scRNA-seq + TCR Integration Analysis

Single-cell RNA-seq and T cell receptor (TCR) sequencing analysis pipeline for characterising tumour-infiltrating T cell states, clonal dynamics, and cell-cell interactions in the cancer microenvironment, using the Wu et al. 2020 (Nature) dataset.

## What this paper is about

Your immune system has a "search and destroy" branch called **T cells**. Each T cell carries a unique receptor (like a barcode/key) that lets it recognize one specific target — a piece of a virus, a tumor protein, etc. When a T cell finds its target, it multiplies — making many copies of itself, all carrying the same receptor. This is called **clonal expansion**, and a group of identical-receptor cells descended from one original cell is a **clone**.

Cancer immunotherapy drugs (like the "checkpoint inhibitors" mentioned in this paper) work by releasing the brakes on T cells so they can attack tumors more effectively. But doctors don't have a good way to predict, in advance, which patients will respond.

**The core question this paper asks:**

> If we look at a patient's **blood**, can we tell which T cell clones are also active *inside the tumor* — without having to biopsy the tumor itself?

This matters because drawing blood is easy and can be done repeatedly; biopsying a tumor is invasive and usually only done once.

**What they did:** they took samples from 14 cancer patients (lung, colorectal, endometrial, and kidney cancers) from three places — **tumor** tissue, **normal tissue right next to the tumor** (NAT), and **blood** (only available for 4 of the patients). For every T cell in every sample, they measured two things:
1. **What "job" is this cell doing right now?** (via RNA sequencing — reading which genes are switched on, which tells you if a cell is "exhausted," "actively killing," "resting," etc.)
2. **What's this cell's receptor barcode?** (via TCR sequencing — reading the unique receptor sequence, which tells you which clone it belongs to)

Combining these lets them ask: *do clones that are expanding in the tumor also show up, expanded, in the blood? And what "job" are those cells doing?*

**What they found:**

1. **"Dual-expanded" clones are the key signal.** Clones found expanding in *both* the tumor and the adjacent normal tissue are disproportionately the same clones found expanded in blood. So blood does carry an echo of what's happening in the tumor — but only for this specific subset of clones.

2. **Those blood-detectable clones tend to be "effector" cells** — the actively-fighting type, not the "exhausted"/burned-out type that dominates deep inside tumors. So blood preferentially shows you the *healthier, more functional* arm of the anti-tumor response.

3. **They built a "gene signature"** (a 30-gene fingerprint) characteristic of these dual-expanded clones, and tested it against data from actual clinical trials of checkpoint-inhibitor drugs. Patients whose tumors showed this signature tended to **respond better to treatment and live longer without disease progression**.

**Why it matters:** the takeaway is a proof-of-concept that a relatively simple blood test — looking for this "dual-expansion" signature among T cell clones — could one day help predict whether a patient will benefit from immunotherapy, without needing a tumor biopsy.

## What the data actually looks like

Every cell in this dataset has two completely different types of measurement attached to it. Both are indexed by the same **cell barcode** (its unique ID, e.g. `RT2_AAACCTGAGGTAGCCA-1`) — the scRNA-seq data has one row per barcode, and the scTCR-seq data has one row per barcode, so matching barcodes across the two tables is how Notebook 3 joins them into a single merged dataset:

**1. scRNA-seq ("gene expression") — what is this cell *doing*?**

For each cell, this is a row of numbers — one per gene (~20,000 genes) — saying how many copies of that gene's mRNA were detected. For example, two real cells' rows might look like:

| Barcode | LOC283788 | RBCK1 | PSMF1 | SNRPB | PCED1A | ... |
|---|---|---|---|---|---|---|
| `RT2_AAACCTGAGGTAGCCA-1` | 1 | 1 | 1 | 2 | 1 | ... |
| `RT2_AAACCTGCAGGTCCAC-1` | 0 | 0 | 2 | 1 | 0 | ... |

Most of the ~20,000 per-gene values are 0 (a gene simply wasn't detected in that cell). This raw gene-count matrix is the *only* input — everything else (the `leiden` cluster from Notebook 1, then `cell_type` and the continuous **state scores** `Exhausted_score`/`Effector_score`/`Naive_score` from Notebook 2) is *derived* from these ~20,000 numbers per cell, not part of the raw data itself.

**2. scTCR-seq ("receptor sequencing") — what is this cell *programmed to recognize*?**

For each cell, this captures the sequenced T cell receptor chain(s) recovered for that cell — one record per chain (alpha and/or beta), with fields like:

| Barcode | locus | v_call | j_call | junction_aa (CDR3) | duplicate_count |
|---|---|---|---|---|---|
| `RT2_AAACCTGAGGTAGCCA-1` | TRB | TRBV2 | TRBJ2-1 | CASSGGYYNEQFF | 6 |
| `RT2_AAACCTGAGGTAGCCA-1` | TRA | TRAV12-2 | TRAJ42 | ... | ... |

`v_call`/`j_call` are which V/J gene segments were used to build the receptor, and `junction_aa` (the CDR3 amino acid sequence) is the actual "barcode/key" sequence that determines what the receptor recognizes.

Cells with the *identical* `junction_aa` sequence(s) are grouped into the same **clonotype** — these are the "clones" described above. Everything below is *derived* from the raw chain records above, computed by Notebook 3:
- `clone_id` — which clonal family this cell belongs to (e.g. `clone_id = 0`), assigned by grouping cells with identical CDR3 sequences
- `clone_id_size` — how many cells share that clone_id across the whole dataset (e.g. `310`)
- `clonal_expansion` — a bucketed version of clone size: `singleton` (44,937 cells, size=1), `small`, `medium`, or `large` (20,272 cells, the most expanded clones)
- `receptor_type` — QC flag derived from how many/which chains were recovered (`"TCR"` = clean paired alpha+beta, vs. `"no IR"`, `"multichain"`, `"ambiguous"`)

**Putting it together:** after joining on the barcode, a single row in the merged dataset (`data/03_tcr_integrated.h5ad`) looks like (the `cell_type`/`Exhausted_score`/`Effector_score` columns come from the scRNA-seq side above; `clone_id`/`clone_id_size`/`clonal_expansion` come from the scTCR-seq side):

| barcode | patient | source | cell_type | Exhausted_score | Effector_score | clone_id | clone_id_size | clonal_expansion |
|---|---|---|---|---|---|---|---|---|
| `RT2_AAACCTGAGGTAGCCA-1` | Renal2 | Tumor | CD8 T cell | -0.39 | 0.81 | 0 | 310 | large |

This is what makes the whole analysis possible: **for the same cell**, you know both *what job it's doing* (from gene expression / state scores) *and* *which clonal family it belongs to and how expanded that family is* (from TCR sequencing) — so you can directly ask "are highly-expanded clones (large `clone_id_size`) more exhausted (`Exhausted_score`) than singletons?"

## Overview

| Notebook | Analysis |
|---|---|
| `01_preprocessing.py` | QC filtering, normalization, PCA, UMAP, Leiden clustering |
| `02_cell_annotation.py` | Marker-based cell type annotation, T cell subset scoring |
| `03_tcr_analysis.py` | Clonotype integration, clonal expansion, **clonotype → transcriptional state linkage** |
| `04_cell_interaction.py` | Ligand-receptor inference (LIANA), niche analysis, expanded vs non-expanded T cell interactions |
| `05_figure1_reproduction.py` | Reproduction of Wu et al. Figure 1 — per-tissue-compartment clonal dynamics (tumour/NAT/blood) for the 4 patients with matched blood samples |
| `06_figure2_cluster_expansion.py` | Reproduction of Wu et al. Figure 2 — CD8 T cell subcluster composition and exhaustion/effector state of dual-expanded clones, split by peripheral (blood) expansion status |

For a checklist of exactly which paper numbers each notebook reproduces (and how close we get), see [`PAPER_SUMMARY.md`](PAPER_SUMMARY.md).

## Dataset

**Wu et al. 2020 (Nature)** — peripheral and tumour-infiltrating T cells from 14 treatment-naive patients across 4 cancer types: non-small-cell lung cancer, endometrial cancer, colorectal cancer, and renal cell carcinoma. Contains paired scRNA-seq and TCR (alpha + beta chain) sequencing across tumour, normal adjacent tissue (NAT), and (for 4 patients) blood.

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
python 05_figure1_reproduction.py
python 06_figure2_cluster_expansion.py
```

Or convert to Jupyter notebooks with jupytext:

```bash
pip install jupytext
jupytext --to notebook notebooks/01_preprocessing.py
jupyter notebook
```

## Notebook walkthroughs

Each notebook reads the `.h5ad` (or `.csv`) file written by the previous step and writes a new one. Notebooks 1-4 build the annotated, clonotype-linked dataset; Notebooks 5-6 use that dataset to reproduce specific figures from the paper.

### Notebook 1 — QC and Preprocessing (`01_preprocessing.py`)

**Goal:** start from raw single-cell gene expression counts and end up with a clean, normalized dataset that's been clustered into groups of similar cells.

1. **Load the data.** scirpy downloads the Wu et al. 2020 dataset directly — a "MuData" object bundling gene expression (GEX) and TCR sequencing data together. We split off just the GEX part for now (138,500 cells x ~33k genes before filtering).
2. **Quality control (QC).** Every cell is really a tiny droplet that may or may not have worked correctly. We compute three QC metrics per cell — number of genes detected, total UMI (transcript) counts, and % of reads from mitochondrial genes — and remove cells that look like:
   - **Empty droplets** (too few genes, < 200) — no real cell was captured.
   - **Doublets** (too many genes, > 6,000) — two cells captured together, counted as one.
   - **Dying/stressed cells** (> 20% mitochondrial reads) — cells that were dying when captured leak cytoplasmic RNA and over-represent mitochondrial transcripts.

   We also drop genes detected in fewer than 3 cells (too rare to be informative). This leaves 138,500 cells and 20,479 genes.

   ![QC distributions before filtering](figures/01_qc_distributions_before.png)
   ![QC distributions after filtering](figures/01_qc_distributions_after.png)

3. **Normalization.** Each gene's count is the number of UMIs (Unique Molecular Identifiers) detected for it — a UMI is a random tag attached to each mRNA molecule before PCR, so identical-UMI reads can be collapsed back into one original molecule, giving an accurate count rather than an inflated PCR-copy count. The raw counts are kept in a separate `counts` layer for later steps (like differential expression) that need them.
   - **Library-size normalization.** Total UMI counts per cell ("library size") vary just from technical noise in capture and sequencing depth — a cell sequenced 2x as deeply will show ~2x bigger numbers for every gene, even if it's biologically identical to another cell. We remove this effect by normalizing each cell's counts to sum to 10,000 ("counts per 10k").
   - **Log transform (`log1p`).** Expression values are heavily skewed, with a few highly-expressed genes (e.g. ribosomal genes) reaching counts in the thousands while most genes sit near 0. On the raw scale, a jump from 5,000→6,000 counts would dwarf a biologically meaningful jump from 2→3 counts. Taking the log (`log(x + 1)`, where the "+1" avoids `log(0)`) makes equal *fold-changes* (e.g. doubling) count equally regardless of starting value, preventing a handful of highly-expressed genes from dominating PCA/clustering and making the data behave more like the roughly-normal distributions most statistical methods assume.
4. **Highly variable gene (HVG) selection.** Most of the ~20k genes are either always-on "housekeeping" genes or always-off, and don't help distinguish cell types. We pick the top 3,000 genes that vary the most across cells — these carry most of the biological signal.
5. **Dimensionality reduction (PCA).** Even 3,000 genes is too many dimensions to cluster directly, and many of those genes are correlated with each other (e.g. all the genes that turn on together when a cell is "activated"). PCA compresses the 3,000-gene expression profile of each cell down to 50 "principal component" numbers that capture most of the variation, with much less redundancy.
6. **Neighbours, UMAP, and clustering.** Using those 50 PCs:
   - We build a "nearest neighbours" graph connecting each cell to the cells most similar to it: for each cell, we measure the (Euclidean) distance to every other cell in the 30-PC space, and connect it to its 15 closest neighbours (`sc.pp.neighbors(n_neighbors=15, n_pcs=30)`). Scanpy uses an approximate nearest-neighbour search (UMAP's algorithm, based on nearest-neighbour-descent) so this is fast even with ~138,500 cells, and converts the raw distances into edge weights ("connectivities") so cells that are very close together are more strongly connected.
   - **UMAP** projects this high-dimensional graph down to 2D so it can be plotted — nearby points are similar cells.
   - **Leiden clustering** uses the same graph to group cells into discrete clusters (here, at `resolution=0.5`) — this is the unsupervised step that finds "groups of similar cells" without knowing in advance what cell types exist.

   ![UMAP Leiden clusters](figures/01_umap_leiden.png)

7. **Save.** The result — UMAP coordinates, PCA, and a `leiden` cluster label per cell — is written to `data/01_preprocessed.h5ad` for Notebook 2.

> **Note:** at this stage the clusters are just numbers (0, 1, 2, ...) grouped by transcriptional similarity. We don't yet know what *kind* of cell each cluster represents — that's Notebook 2's job.

### Notebook 2 — Cell Type Annotation (`02_cell_annotation.py`)

**Goal:** figure out what cell type each Leiden cluster from Notebook 1 actually is, by checking which genes it expresses.

1. **Marker genes.** Decades of immunology have established that specific cell types reliably turn on specific genes — e.g. CD8 T cells express `CD8A`/`CD8B`, B cells express `CD19`/`MS4A1`, monocytes express `CD14`/`LYZ`. We define a dictionary of these "marker genes" for the major immune cell types we expect to see (CD8 T, CD4 T, NK, B cell, monocyte, DC, plasma cell).

   **Why these 7 categories?** Two reasons:
   - **They're the standard major immune lineages.** T cells, B cells, NK cells, myeloid cells (monocyte/DC), and plasma cells (antibody-secreting B cells) are the cell types any blood or tumour-infiltrating immune sample is expected to contain — this is the same set of categories used in general-purpose immune reference atlases (e.g. 10x's PBMC references, Azimuth), not something specific to this dataset.
   - **They empirically cover almost all the data.** 17 of 18 Leiden clusters (132,532 / 138,500 cells, ~95.7%) get a confident match to one of these 7 categories via their marker genes — only cluster 11 (4.3% of cells) didn't fit any of them (see step 4 below). That high coverage is itself evidence the panel is well-matched to what's actually present in this dataset.

   **Reference:** [Azimuth Human PBMC reference](https://azimuth.hubmapconsortium.org/references/#Human%20-%20PBMC) defines marker genes for these same lineages (T/B/NK/Monocyte/DC/Plasma). PBMC is the relevant comparison rather than a tissue-specific reference like kidney: our markers are *immune lineage* markers (e.g. `CD8A`/`CD8B` for CD8 T cells, `CD14` for monocytes), which define a cell's identity regardless of which organ it's sampled from — unlike tissue-specific references (kidney, lung, etc.), which mostly distinguish *non-immune* parenchymal cell types (tubule cells, podocytes, etc.) that aren't relevant here. PBMC also needs far fewer markers per type than a tissue reference because we're only separating ~7 broad, transcriptionally distinct lineages, not dozens of similar subtypes.
2. **Dot plot / violin plot.** We plot average marker gene expression for each Leiden cluster. A cluster that lights up strongly for the CD8 marker panel and not the others is called a CD8 T cell cluster, and so on.

   ![Marker dotplot](figures/02_marker_dotplot.png)

3. **Assign labels.** Each of the 18 Leiden clusters gets a `cell_type` label based on its marker profile (`CLUSTER_ANNOTATION` dict). One subtlety found during this project: this dataset is **CD3-sorted** (i.e. essentially all cells are T cells to begin with), so some "NK marker" genes (`GNLY`, `NKG7`, `KLRD1` — cytotoxic granule genes) are *also* highly expressed by cytotoxic CD8 effector T cells, and several clusters that superficially looked like NK/monocyte/B/DC/plasma cell were actually CD8 T cells. Getting this annotation right mattered a lot — it changed the CD8 T cell count from ~34k to ~80k cells in later notebooks (see Notebook 6 below).

   ![Cell type UMAP](figures/02_umap_celltypes.png)

4. **Cross-check with differential expression.** The marker panel in step 1 comes from prior knowledge, not from this data. `sc.tl.rank_genes_groups` (Wilcoxon, one cluster vs. rest) instead computes each cluster's *actual* top distinguishing genes with no priors, as a sanity check on the marker-based labels. Cluster 11's top DE genes were ribosomal/mitochondrial genes plus `RPS4Y1` (not a cell-identity marker), with no NK markers — its label is `"Unknown"`. For the cell-type panel (step 1, `MARKERS`), this cross-check didn't suggest adding any new major cell types — every other cluster's top DE genes land cleanly into CD8 T / CD4 T / B / Plasma / Monocyte-DC.

   ![Rank genes groups dotplot](figures/02_rank_genes_groups_dotplot.png)

5. **T cell subset scoring.** Within the T cells, we further score each cell for five *functional states* using small marker gene sets — Naive, Memory, Effector, Exhausted, Regulatory (Treg) — via `sc.tl.score_genes`, which essentially averages expression of each gene set per cell (relative to a random background). These scores (`Exhausted_score`, `Effector_score`, `Naive_score`, etc.) are continuous values per cell, not hard labels, and are what later notebooks use to ask "are expanded clones more exhausted or more effector-like?"

   **Worked example for `Exhausted_score`** (gene set = `PDCD1, HAVCR2, TIGIT, LAG3, CTLA4`), for one real cell (`RT2_AAACCTGCAGGTCCAC-1`):

   | Barcode | PDCD1 | HAVCR2 | TIGIT | LAG3 | CTLA4 | mean of gene set | background mean (random genes, similar expression level) | `Exhausted_score` |
   |---|---|---|---|---|---|---|---|---|
   | `RT2_AAACCTGCAGGTCCAC-1` | 1.61 | 1.61 | 2.83 | 2.83 | 1.61 | 2.10 | 0.50 | **1.60** |

   `Exhausted_score = mean(gene set expression) - mean(background gene set expression)`. This is repeated for each of the 5 gene sets, so the **output format** is 5 new numeric columns added to `adata.obs` — one row per cell:

   | Barcode | Exhausted_score | Effector_score | Naive_score | Memory_score | Regulatory_score |
   |---|---|---|---|---|---|
   | `RT2_AAACCTGCAGGTCCAC-1` | 1.60 | 0.38 | ... | ... | ... |

   ![T cell subset markers](figures/02_tcell_subsets_dotplot.png)

6. **Save.** Cell type labels + functional scores are added as new columns and written to `data/02_annotated.h5ad`.

### Notebook 3 — TCR Analysis and Clonotype-State Linking (`03_tcr_analysis.py`)

**Goal:** bring in the T cell receptor (TCR) sequencing data, group cells into clonal families ("clonotypes"), and ask whether clonally-expanded cells look transcriptionally different from non-expanded ones. This is the analytical core of the whole project.

All of the TCR-specific steps below use **[scirpy](https://scirpy.scverse.org)** (conventionally imported as `ir`, for "immune receptor") — a scanpy-compatible package purpose-built for TCR/BCR repertoire analysis.

1. **Load and merge.** We reload the original MuData (which has the TCR/"airr" modality scirpy needs) and transfer over the cell type + score annotations from Notebook 2, matching cells by their barcode (each cell's unique ID).
2. **TCR chain QC (`receptor_type`).** Each T cell's receptor is built from two protein chains — alpha and beta — each encoded by its own piece of rearranged DNA (its CDR3 sequence). Sequencing doesn't always recover both chains for every cell (dropout, multiple chains detected, etc.), so scirpy's `ir.tl.chain_qc` classifies each cell's TCR reads as e.g. "TCR" (clean paired alpha+beta), "ambiguous", "multichain", or "no IR" (no receptor recovered). We restrict everything downstream to the **95,314 cells with a clean paired alpha+beta TCR**, since those are the only ones we can confidently assign a clonotype to.

   ![Chain QC receptor type](figures/03_chain_qc_receptor_type.png)

3. **Define clonotypes.** A "clonotype" = all cells whose alpha+beta CDR3 sequences are identical, meaning they're descendants of the same original T cell (same antigen specificity). Scirpy's `ir.tl.define_clonotypes` groups the 95,314 paired cells into **53,142 unique clonotypes**.
4. **Clonal expansion.** Most clonotypes are seen only once ("singleton" — clone size 1); a few are seen many times ("expanded"). We bucket clonotypes by size: singleton (1), small (2-5), medium (6-20), large (>20). **20,272 cells (21%) belong to "large" clones** — these are the clones that proliferated heavily, presumably because they recognized something (e.g. a tumour antigen) and got triggered to divide.

   | Category | Cells | Description |
   |----------|-------|-------------|
   | Singleton | 44,937 | Clone size = 1 |
   | Small | 17,155 | Clone size 2–5 |
   | Medium | 12,950 | Clone size 6–20 |
   | **Large** | **20,272** | Clone size > 20 — likely tumour-reactive |

   ![Clonal expansion UMAP](figures/03_umap_clonal_expansion.png)
   ![Expansion by cell type](figures/03_expansion_by_celltype.png)

5. **Clonotype → transcriptional state linkage (the core question).** For each expansion bucket (singleton/small/medium/large), we compare the distribution of `Exhausted_score`, `Effector_score`, and `Naive_score` from Notebook 2.

   ![State by expansion](figures/03_state_by_expansion.png)

   Each panel is a box plot: the **red line is the median** score across all T cells in that bucket, the **box is the interquartile range**, and the dots above are individual cells with unusually high scores. Quick refresher on what the scores mean biologically:
   - **Naive** = genes (`CCR7`, `TCF7`, `SELL`, `LEF1`) marking T cells that have never encountered their target antigen yet — "fresh," uncommitted cells.
   - **Effector** = genes (`GZMB`, `GZMK`, `PRF1`, `IFNG`, `TNF`) marking T cells actively "fighting" — producing cytotoxic/killing molecules and inflammatory signals, the hallmark of a T cell that has recently recognized its target.
   - **Exhausted** = genes (`PDCD1`, `HAVCR2`, `TIGIT`, `LAG3`, `CTLA4` — classic immune checkpoint genes) marking T cells that have been stimulated by antigen for so long they've become functionally "burned out" and less effective.

   The **cleanest, most monotonic trend is in Effector and Naive**: median Effector score rises steadily from singleton (≈ -0.82) → small (≈ -0.52) → medium (≈ -0.33) → large (≈ -0.20), while median Naive score falls steadily from singleton (≈ +0.02) down to large (≈ -0.28). In other words, **the bigger a clone is, the more its cells look like actively-fighting effector cells and the less they look like naive, never-triggered cells** — which makes biological sense, since a clone only grows large by being repeatedly triggered to divide by its target antigen.

   The Exhausted score trend is weaker and not strictly monotonic: it rises from singleton (≈ -0.39) to small (≈ -0.27) to medium (≈ -0.21), but then **dips back down for large clones** (≈ -0.30) — close to the "small" value rather than continuing to climb. So expanded clones (small/medium/large) do sit above singletons in exhaustion, but exhaustion doesn't keep increasing with clone size the way the Effector score does; the very largest clones are slightly *less* exhausted than medium ones, possibly because they're dominated by the still-highly-functional effector cells described above.

6. **Save.** The merged, clonotype-annotated dataset (95,314 cells, with `clone_id`, `clone_id_size`, `clonal_expansion`, and all the Notebook 2 scores) is written to `data/03_tcr_integrated.h5ad` — this is the file Notebooks 5 and 6 build on.

### Notebook 4 — Cell-Cell Interaction Analysis (`04_cell_interaction.py`)

**Goal:** go beyond individual cells and ask how different cell types are *communicating* with each other via ligand-receptor signalling — and whether expanded T cell clones have a distinct "conversation" with their neighbours compared to non-expanded clones.

1. **Ligand-receptor inference with LIANA.** Cells signal to each other when one cell expresses a "ligand" (a secreted or surface protein) and a neighbouring cell expresses the matching "receptor". `li.mt.rank_aggregate` scores thousands of known ligand-receptor pairs across every pair of cell types, combining several established methods (CellPhoneDB, CellChat, NATMI) into one consensus ranking — **4,517 ligand-receptor pairs scored** across cell type pairs.
2. **Dot plot of top CD8 T cell interactions.** Shows the strongest-scoring ligand-receptor pairs involving CD8 T cells as either sender or receiver — e.g. checkpoint molecules, cytokines, adhesion molecules exchanged with monocytes/DCs/other T cells.

   ![CD8 interactions](figures/04_cd8_interactions_dotplot.png)

3. **Interaction network.** A circle plot summarizing the top 200 interactions as a network between cell types — a birds-eye view of which cell type pairs talk to each other the most.

   ![Interaction network](figures/04_interaction_network.png)

4. **Niche analysis: expanded vs non-expanded CD8 T cells.** Using the `clonal_expansion` labels from Notebook 3, we split CD8 T cells into "expanded" (medium/large clones) and "non-expanded" (singleton/small), then re-run LIANA separately for each group and compare their top interaction partners. The idea: if expanded clones are undergoing a different kind of immune response (e.g. more activation/exhaustion signalling), their communication profile with surrounding cells should look different too.

   ![Expanded vs non-expanded interactions](figures/04_expanded_vs_nonexpanded_interactions.png)

   The top interactions in both groups are dominated by **HLA class I → CD8A/CD8B** — this is just "CD8 T cells engage MHC-I," expected by definition and not informative about expansion status.

5. **Expansion-specific signals: FASLG → TNFRSF1A and CCL5 → SDC1.** Beyond the generic HLA→CD8 interactions shared by both groups, two signals differ specifically by expansion status.

   **FASLG → TNFRSF1A (FasL/Fas cytotoxic killing signal):**

   ![FASLG TNFRSF1A rank](figures/04_faslg_tnfrsf1a_rank.png)

   **Key finding:** FASLG → TNFRSF1A is detected in the expanded CD8 T cell group (rank 0.040) but completely absent from the non-expanded top interactions. FasL (FASLG/CD95L) is one of the primary cytotoxic effector mechanisms of CD8 T cells — it binds Fas receptor (TNFRSF1A) on target cells and triggers caspase-dependent apoptosis. Its presence only in expanded clones directly corroborates the Notebook 6 finding that expanded+blood-detected (D-class) clones are enriched for an effector/cytotoxic phenotype.

   **CCL5 → SDC1 (CD8 T cell → Plasma cell):**

   ![CCL5 SDC1 rank](figures/04_ccl5_sdc1_rank.png)

   **Key finding:** the CCL5 → SDC1 magnitude rank is **~8× stronger in expanded clones** (0.010 vs 0.080). The expanded and non-expanded groups are nearly equal in size (~25k cells each), ruling out a cell-count artefact. CCL5 is a known tertiary lymphoid structure (TLS) determinant chemokine — TLS are ectopic immune aggregates in tumours containing CD8 T cells and plasma cells, and TLS presence correlates with better immunotherapy response. Expanded clones signalling more strongly to plasma cells via CCL5 is consistent with a TLS coordination role, though experimental validation would be needed.

   Together, these two signals paint a coherent picture: expanded clones are not just "more cells" but are transcriptionally and functionally distinct — engaging both direct cytotoxic killing (FasL/Fas) and chemokine-mediated immune coordination (CCL5→SDC1).

6. **CD8 vs CD4 T cell sender profiles.** Using the same LIANA result, compare the top ligand-receptor pairs *sent* by CD8 T cells vs CD4 T cells. CD8 T cells are primarily cytotoxic effectors; CD4 T cells provide "help" via cytokines and co-stimulation. This comparison tests whether LIANA recovers this known functional division from gene expression alone — a sanity check on the interaction scores and a broader view of T cell communication than the CD8-only analysis above.

   ![CD8 vs CD4 sender profiles](figures/04_cd8_vs_cd4_sender_profiles.png)

7. **Myeloid → T cell signalling.** Extract the top Monocyte/DC → CD8/CD4 T cell interactions — the "other direction" of the microenvironment. Monocytes and DCs are key regulators of T cell fate: they can activate T cells (antigen presentation via MHC-II, co-stimulation via CD80/CD86 → CD28) or suppress them (checkpoint ligands like PD-L1 → PD-1, LGALS9 → HAVCR2/TIM3). This gives the full picture of the niche — not just what T cells are doing, but what the surrounding myeloid cells are sending to shape T cell function.

   ![Myeloid to T cell interactions](figures/04_myeloid_to_tcell_interactions.png)

   **Key findings:**
   - **Antigen presentation dominates** (HLA-A/C → CD8A/CD8B, strongest signal) — Monocyte/DCs are actively presenting antigen to CD8 T cells, their classical role as APCs.
   - **Simultaneous suppression via LGALS1 → CD69/PTPRC** — galectin-1 is an immunosuppressive signal that inhibits T cell activation and promotes apoptosis. The same myeloid population presenting antigen is also dampening the T cell response — the "double agent" role of tumour-associated myeloid cells.
   - **T cell recruitment via HMGB1 + MIF → CXCR4** — both a danger signal (HMGB1, from stressed/dying cells) and a pro-inflammatory myeloid cytokine (MIF) signal through CXCR4 to recruit T cells into the tumour.

   The combined picture: Monocyte/DCs simultaneously recruit, activate, and suppress T cells — the same myeloid population both drives and restrains the anti-tumour T cell response.

This notebook is more exploratory than 1-3 — it doesn't feed numerically into Notebooks 5/6, but it demonstrates the cell-cell communication side of single-cell analysis, which is commonly used to generate hypotheses about the tumour microenvironment.

### Notebook 5 — Reproducing Figure 1 (`05_figure1_reproduction.py`)

**Goal:** reproduce the paper's Figure 1 ("Parallel dual expansion and peripheral clonal expansion"), which classifies clonotypes by a **per-tissue-compartment expansion scheme** — a different, finer-grained categorisation than the singleton/small/medium/large buckets used in Notebook 3.

**The key finding this notebook is after:** "dual-expanded" clones — clones found expanding in *both* the tumour and the adjacent normal tissue — are disproportionately the same clones found expanded in blood. So a blood draw does carry an echo of what's happening inside the tumour, but only for this specific subset of clones (not for clones that expanded in only one tissue compartment, or not at all). Panel (e) below is the direct evidence for this.

Of the 14 patients, only **Renal1, Renal2, Renal3, and Lung6** have matched **blood** samples, so this notebook restricts to those 4.

1. **Load and subset.** Load `data/03_tcr_integrated.h5ad` and keep only cells from the 4 blood-sample patients that have a defined `clone_id`.
2. **Classify each clonotype by tissue and blood expansion pattern**, following the "Clone size Tu/NAT/Bl" scheme in the paper's Fig 1b legend:
   - `n`/`N` — NAT singleton/multiplet (NAT only)
   - `t`/`T` — tumour singleton/multiplet (tumour only)
   - `D` — **dual-expanded**: detected (≥1 cell) in **both** NAT and tumour
   - `x` — not detected in NAT or tumour (blood-only clone)
   - `b`/`B` — blood non-expanded (1 cell) / expanded (>1 cell)

   Two corrections to the methodology brought this reproduction much closer to the paper's reported values:
   1. **`D` requires only ≥1 cell in both NAT and tumour** (not *multiplet* in both, as first implemented) — this alone brought `n_D` counts within 4-10% of the paper.
   2. **Panel (e) is cell-weighted, not clonotype-weighted** — per the Fig 1 legend ("fractions of tissue-resident cells with clonotypes observed in a blood-expanded clone"). Re-weighting by cell count (so large clones dominate, as in the paper) closed most of the remaining gap.

3. **Panel (a) — clonotype frequency in tumour vs NAT**, coloured by tissue expansion pattern, with `n_D`, `r` (correlation), and D-detection annotated per patient:

   ![Clonotype tissue scatter](figures/05_clonotype_tissue_scatter.png)

   Each **dot is one clonotype** from one patient: its x-position is the clone's normalized size (cell fraction) in NAT, and its y-position is its normalized size in tumour, with a bit of random jitter so overlapping points (especially the singletons piled up near the origin) are visible. Colour marks the clone's tissue expansion category (`n`/`N`/`t`/`T`/`D`, as defined above). Clones that fall on or near the diagonal — present and expanded in *both* NAT and tumour, the **dual-expanded (D)** clones — are the ones counted in `n_D`, and `r` is the correlation between their NAT and tumour clone sizes.

   | Patient | n_D (ours/paper) | r (ours/paper) | D-detection (ours/paper) |
   |---|---|---|---|
   | Renal1 | 107 / 108 | 0.80 / 0.81 | 0.79 / 0.83 |
   | Renal2 | 130 / 139 | 0.85 / 0.93 | 0.83 / 0.87 |
   | Renal3 | 110 / 119 | 0.65 / 0.58 | 0.51 / 0.63 |
   | Lung6 | 216 / 213 | -0.05 / 0.31 | 0.03 / 0.05 |

4. **Panel (e) — detection of tissue-resident TCRs in blood**: for each tissue expansion category, what fraction of tissue cells belong to clonotypes also observed as blood-expanded:

   ![Tissue TCR detection in blood](figures/05_tissue_tcr_detection_in_blood.png)

   Within ~5-12% of the paper for every patient, including the directionally-correct Lung6 outlier (D-detection 0.03 vs 0.05, see table in step 3).

5. **Panel (f) — blood vs tumour clone size correlation.** The paper pools clonotypes shared between blood and tumour across all 4 patients into a single scatter (n=328, r=0.72, P=6.3e-53), rather than per-patient plots:

   ![Blood tumour clone correlation](figures/05_blood_tumor_clone_correlation.png)

   Each dot is a clonotype detected in **both** blood and tumour; x = its clone size in blood, y = its clone size in the tumour. The strong correlation means a clone's size in blood is a good predictor of its size in the tumour — the most direct evidence for the "blood test" idea.

   Reproducing this directly: **n=309 (paper 328), r=0.74 (paper 0.72)** — an excellent match.

6. **Save.** The per-clonotype tissue/blood pattern classification is written to `data/05_figure1_clone_classification.csv` for Notebook 6.

**Overall**, the core biological finding — dual-expanded tissue clones dominate the peripherally-expanded blood repertoire, and this relationship is strong in renal patients but essentially absent in Lung6 — reproduces with good quantitative fidelity.

**Why this matters:** the takeaway is a proof-of-concept that a relatively simple blood test — looking for this "dual-expansion" signature among T cell clones — could one day help predict whether a patient will benefit from immunotherapy, without needing a tumour biopsy.

### Notebook 6 — Reproducing Figure 2 (`06_figure2_cluster_expansion.py`)

**Goal:** reproduce the paper's Figure 2 finding — that among **dual-expanded (D)** clones (from Notebook 5), the ones that are also **peripherally (blood-)expanded (B)** skew towards an **effector-like CD8 phenotype (8.1-Teff)**, while D clones without blood detection skew towards more exhausted/Trm-like states.

**The key finding this notebook is after:** those blood-detectable clones tend to be "effector" cells — the actively-fighting type, not the "exhausted"/burned-out type that dominates deep inside tumours. So blood preferentially shows you the *healthier, more functional* arm of the anti-tumour response. Steps 5 and 6 below are the direct evidence for this.

1. **Load and merge.** Load `data/03_tcr_integrated.h5ad` and join in the tissue/blood pattern classification from Notebook 5 (`data/05_figure1_clone_classification.csv`).
2. **Define comparison groups.** Restrict to **dual-expanded (`D`)** clones in the 4 blood-sample patients, split into:
   - **D & B** — dual-expanded AND peripherally (blood-)expanded (3,083 cells)
   - **D & not-B** — dual-expanded but blood singleton/undetected (1,906 cells)
3. **Subcluster CD8 T cells at higher resolution.** The dataset-wide Leiden clustering (resolution 0.5, Notebook 1) only yields 4 CD8 T cell clusters — too coarse to separate an "effector" subtype from "exhausted"/"naive" subtypes the way the paper's 33-cluster analysis does. So we re-normalize, re-select HVGs, and re-cluster the CD8 T cell subset alone (PCA + **Harmony integration across patients** + Leiden at resolution 1.0) into **14 subclusters**.

   > **Two methodology fixes were needed to get a meaningful result here:**
   > - **Cell-type annotation fix** (Notebook 2, step 3): correcting the CD3-cross-reactivity mislabeling raised the CD8 T cell count from 33,949 to 79,866 (51,190 of which belong to the 4 blood-sample patients).
   > - **Patient integration**: an earlier version subclustered CD8 cells with plain PCA/Leiden (no batch correction), pooling all 14 patients. This produced an "Effector-like" subcluster that happened to contain **zero** cells from the 4 blood-sample patients' dual-expanded clones — a patient-driven clustering artifact that made the cluster-composition test degenerate (chi2 dof=0, p=1.0). Adding **Harmony integration** (`sc.external.pp.harmony_integrate`, batch key = `patient`) fixed this.

4. **Characterise CD8 subclusters.** Compute mean Exhausted/Effector/Naive scores per subcluster; the subcluster with the highest (Effector − Exhausted) score difference is labelled "Effector-like" (subcluster 3: Effector_score=-0.11, Exhausted_score=-0.37 — the highest gap of the 14 subclusters).
5. **Cluster composition: D&B vs D&not-B.** For CD8 T cells in dual-expanded clones, compare the CD8 subcluster distribution between the two peripheral-expansion groups:

   ![CD8 cluster composition](figures/06_cd8_cluster_composition.png)

   | Comparison | Result |
   |---|---|
   | Effector-like fraction (cell-level), D&B vs D&not-B | **50.2% vs 6.7%** (chi2 p=4.7e-219) |
   | Effector-like fraction (per-clone primary cluster), D&B vs D&not-B | **55.8% vs 11.4%** (chi2 p=8.2e-16) |

   The "per-clone" row above follows the paper's approach of assigning each clone a **primary cluster** (the CD8 subcluster with the most cells in that clone), then repeating the composition comparison at the clone level instead of the cell level.

   This is the core claim of the paper's Figure 2 and it reproduces convincingly. **Per-patient sanity check** (D&B vs D&not-B Effector-like %):

   | Patient | D&B Effector-like % (n) | D&not-B Effector-like % (n) |
   |---|---|---|
   | Lung6 | 26.1% (69) | 2.2% (1,469) |
   | Renal1 | 62.5% (40) | 7.1% (14) |
   | Renal2 | 52.3% (2,599) | 27.8% (291) |
   | Renal3 | 38.9% (375) | 9.8% (132) |

   Every patient shows D&B > D&not-B by a wide margin (9-55 percentage points), so the core directional claim holds patient-by-patient — but the pooled chi2 p-value reflects this consistent within-patient effect plus some amplification from patient mix (D&B is 84% Renal2 cells; D&not-B is 77% Lung6 cells), not a single homogeneous effect of that exact size.

6. **Panel d reproduction — clone size of D-pattern 8.1-Teff clones, by blood status.**

   **What question is this answering?** Step 5 showed that *cells* in blood-detectable (D&B) clones skew towards an effector-like phenotype. Panel d asks a related but different question, starting from the *other* direction: if you take only the clones that are **already effector-like (8.1-Teff) AND dual-expanded (D)** — i.e., the population step 5 just showed is enriched in D&B — does their **clone size** (how many cells that clone has, in tumour+NAT) depend on whether the clone is also detected in blood?

   In other words: among effector-like dual-expanded clones, are the blood-detectable ones the *biggest* ones? This is the paper's evidence that blood-expanded 8.1-Teff/D clones aren't just "present" in blood — they're the dominant, most-expanded clones driving tumour infiltration.

   **How it's plotted.** The paper's actual Fig 2d is a **swarm plot**: take all clones whose primary cluster is **8.1-Teff** and whose tissue pattern is **D** (dual-expanded). Each dot is one such clone, positioned by its clone size (total cell count across tumour + NAT). Clones are split into three columns by their blood status:
   - **Ind** (blood-independent) — clone wasn't detected in blood at all
   - **Non** (blood non-expanded) — detected in blood, but only as a singleton
   - **Exp** (blood-expanded) — detected in blood as an expanded (multi-cell) clone

   A Mann-Whitney test compares the clone-size distributions between groups — the question is whether **Exp** clones are significantly larger than **Ind**/**Non** clones.

   We discovered that `data/03_tcr_integrated.h5ad` already contains the paper's own per-cell subtype labels in `obs["cluster_orig"]` — so this reproduction uses `cluster_orig == "8.1-Teff"` directly (majority-vote per clone) instead of our Harmony-derived "Effector-like" proxy, avoiding re-clustering entirely:

   ![Panel d swarm](figures/06_paneld_v2_clonesize_swarm.png)

   | Blood group | n (ours / paper) | Median clone size (tumour+NAT) | Mean clone size |
   |---|---|---|---|
   | Ind (blood-independent) | **16 / 16** | 4.0 | 5.2 |
   | Non (blood non-expanded) | 13 / 22 | 6.0 | 6.1 |
   | Exp (blood-expanded) | 37 / 68 | 20.0 | 35.2 |

   - Mann-Whitney Ind vs Exp: **p=4.8e-5** (paper: P=3.0e-5) — excellent match
   - Mann-Whitney Non vs Exp: **p=4.1e-4** (paper: P=1.4e-5) — same direction, ~1 order of magnitude less significant, consistent with our smaller Non group (n=13 vs paper's 22)

   This is the closest reproduction of any Figure 2 panel in this repo: the **Ind group count matches the paper exactly (n=16)**, and the swarm plot visually reproduces the paper's qualitative pattern — blood-expanded (Exp) 8.1-Teff/D clones have substantially larger tumour+NAT clone sizes than blood-independent or non-expanded ones.

7. **Save.** Various clone-level classification tables are written to `data/06_*.csv`.

## Dataset Summary

| Step | Metric | Value |
|------|--------|-------|
| QC filtering | Cells retained | **138,500** |
| QC filtering | Genes retained | **20,479** |
| Highly variable genes | Selected for PCA | **3,000** |
| Cell annotation | T cells identified | **128,654** (93% of total) |
| TCR integration | Cells with paired TCR | **95,314** |
| Clonotype definition | Unique clonotypes | **53,142** |
| LIANA | Ligand-receptor pairs scored | **4,517** |

## Tools

| Tool | Purpose |
|---|---|
| [Scanpy](https://scanpy.readthedocs.io) | scRNA-seq preprocessing, clustering, visualization |
| [scirpy](https://scirpy.scverse.org) | TCR data integration, clonotype analysis, VDJ gene usage |
| [LIANA](https://liana-py.readthedocs.io) | Ligand-receptor cell-cell interaction inference |
| [muon](https://muon.readthedocs.io) | Multi-modal data handling |
