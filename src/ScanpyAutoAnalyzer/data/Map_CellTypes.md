---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.12.0
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---


```python
import scanpy as sc
from scmappy import common_genes, scmap_annotate

# ----------------------------
# 1. Load reference and query datasets
# ----------------------------
# reference: should already have annotations
adata_ref = sc.read_h5ad("reference.h5ad")
adata_query = sc.read_h5ad("query.h5ad")

# Assume annotations are stored in reference.obs["celltype"]

# ----------------------------
# 2. Align genes between reference and query
# ----------------------------
# "Gene_names" is the key in .var that holds gene identifiers
# replace with the actual column name in your data (often "gene_symbols" or similar)
adata_ref, adata_query = common_genes(
    adata_ref,
    adata_query,
    "Gene_names",
)

# ----------------------------
# 3. Annotate query dataset with scmap
# ----------------------------
scmap_annotate(
    adata_ref,
    adata_query,
    gene_id_col="Gene_names",
    labels_col="celltype",           # must exist in adata_ref.obs
    algorithm_flavor="centroid",     # "centroid" (scmap-cluster) or "nearest_neighbors" (scmap-cell)
    gene_selection_flavor="HVGs",    # highly variable genes
    similarity_threshold=0.7,        # adjust threshold for "unassigned"
    key_added="scmap_annotation"     # results stored here
)

# ----------------------------
# 4. Inspect results
# ----------------------------
print(adata_query.obs["scmap_annotation"].value_counts())

# Each query cell now has a predicted celltype
sc.pl.umap(adata_query, color=["scmap_annotation"])
```