# ScanpyAutoAnalyzer

**ScanpyAutoAnalyzer** is a lightweight Python module providing a curated collection of ready-to-use functions for single-cell RNA-seq data analysis using [Scanpy](https://scanpy.readthedocs.io/), [Anndata](https://anndata.readthedocs.io/), and [Scvelo](https://scvelo.readthedocs.io/).  
It is designed to simplify and standardize preprocessing, visualization, pseudo-bulk aggregation, and cell-type annotation tasks in a reproducible, lab-friendly way.

---

## 📦 Installation

You can install directly from GitHub:

```bash
pip install git+https://github.com/stela2502/ScanpyAutoAnalyzer.git
```

This automatically installs all required dependencies (`scanpy`, `anndata`, `scvelo`, `pandas`, `matplotlib`, `seaborn`, `numpy`, `scipy`, etc.).

---

## 🧬 Quick Start

Once installed, simply import the functions module:

```python
import functions as faa
import scanpy as sc

adata = sc.read_h5ad("my_dataset.h5ad")

# Run the built-in preprocessing workflow
adata = faa.preprocess(adata)
```

This will perform:
- QC metric computation
- Filtering of mitochondrial and ribosomal genes
- Normalization and log-transformation
- Identification of highly variable genes
- UMAP computation and Leiden clustering

---

## 🔬 Example Workflows

### 1. Annotating Cell Types
```python
celltype_dict = {
    "0": "T cells",
    "1": "B cells",
    "2": "Myeloid"
}

adata = faa.add_cell_types_from_dict(adata, cluster_key="leiden", cell_types_dict=celltype_dict)
```

This adds `adata.obs["cell_type"]` and automatically maps colors based on cluster assignments.

---

### 2. Visualizing Data

#### UMAP with Cluster Contours
```python
faa.plot_umap_with_cluster_contours(adata, contour_group="leiden", subset_group="cell_type", group_key_selected="T cells")
```

#### Summary Heatmap
```python
pseudo = faa.compute_cluster_mean_profiles_adata(adata, cluster_key="leiden")
faa.summary_heatmap(pseudo.to_df(), cluster_key="leiden", experiment="myexp")
```

#### Stacked Bar Plot
```python
faa.plot_stacked_bar(adata, value_col="geneSum", index_col="sample", label_col="cell_type", color_map=my_colors)
```

---

### 3. Pseudo-bulk and Differential Expression

Generate pseudo-bulk samples and perform DE analysis per cluster:
```python
faa.write_pseudobulk_stats_tables(adata, key_added="pseudo_bulk_DE")
```

The results are stored in `adata.uns["pseudo_bulk_DE"]` and written to CSV files.

---

### 4. Predicting Cell Types from Marker Genes
```python
marker_groups = {
    "T cells": ["CD3D", "CD3E"],
    "B cells": ["MS4A1"],
    "NK cells": ["NKG7"]
}
faa.predict_celltype(adata, marker_groups, cluster_key="leiden")
```

This computes per-cluster marker expression and generates predicted cell-type labels.

---

## 🧠 Additional Utilities

| Function | Description |
|-----------|-------------|
| `preprocess` | End-to-end preprocessing workflow with QC, normalization, clustering, and UMAP. |
| `add_cell_types_from_dict` | Map cluster IDs to cell types and set colors. |
| `update_obs_colors` | Assign consistent colors to obs categories. |
| `compute_cluster_mean_profiles_adata` | Generate pseudo-cells by averaging expression within clusters. |
| `write_pseudobulk_stats_tables` | Perform DE on pseudo-bulk samples and export results. |
| `plot_stacked_bar`, `plot_dotplot_matrix_legend` | Visual summaries of group-level relationships. |
| `gene_coexpression_store` | Compute and store co-expressed genes for a given marker. |
| `predict_celltype` | Infer likely cell-type labels based on marker genes. |
| `write_data_to_directory` | Export AnnData as 10x-style matrix (MTX + TSVs). |

---

## 🧪 Example Session

```python
import functions as faa
import scanpy as sc

# Load data
adata = sc.read_10x_h5("filtered_feature_bc_matrix.h5")

# Preprocess and cluster
adata = faa.preprocess(adata)

# Assign cell types and visualize
cell_types = {"0": "T cells", "1": "B cells", "2": "Monocytes"}
adata = faa.add_cell_types_from_dict(adata, "leiden", cell_types)

faa.plot_umap_with_cluster_contours(adata, "leiden", "cell_type", "T cells")
```

---

## 🧭 Philosophy

> “Make high-quality single-cell analysis as straightforward as possible — while keeping every step transparent.”

ScanpyAutoAnalyzer is not a framework, but a **curated toolkit**: every function is meant to be inspected, reused, and adapted for your own analysis pipelines.

---

## 📚 Citation

If you use this module in your research, please cite the repository:

```
Lang S. (2025) ScanpyAutoAnalyzer: lightweight Scanpy extensions for streamlined single-cell workflows.
https://github.com/stela2502/ScanpyAutoAnalyzer
```

---

## 🧩 License

This project is licensed under the **BSD 3-Clause License** – you are free to use, modify, and distribute it with attribution.