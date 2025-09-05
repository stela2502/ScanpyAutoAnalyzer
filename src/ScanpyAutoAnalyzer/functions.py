import re
import scanpy
import pandas as pd
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import colorsys
import anndata
import scvelo as scv
from scipy import sparse


def add_cell_types_from_dict(adata, cluster_key, cell_types_dict, new_key="cell_type", color_map=None):
    """
    Adds cell type annotations to an AnnData object based on a cluster→cell type dictionary.
    Also assigns colors to the new cell_type column based on the cluster colors.

    Parameters
    ----------
    adata : AnnData
        The AnnData object containing the clustering results.
    cluster_key : str
        The key in adata.obs that contains cluster assignments (e.g., 'leiden').
    cell_types_dict : dict
        Dictionary mapping cluster IDs (as strings) to cell types.
    new_key : str, optional
        The name of the new column in adata.obs where cell types will be stored. Default is 'cell_type'.
    color_map : dict, optional
        Optional dictionary mapping cluster IDs to colors. Overrides automatic cluster color mapping.

    Returns
    -------
    AnnData
        The same AnnData object with new annotations and colors added.
    """
    
    # Ensure cluster_key exists
    if cluster_key not in adata.obs:
        raise KeyError(f"{cluster_key} not found in adata.obs")

    # Map clusters to cell types
    adata.obs[new_key] = adata.obs[cluster_key].map(cell_types_dict).astype("category")

    # Determine colors
    if color_map is not None:
        # Use the provided color_map
        cell_type_to_color = {}
        for cl, color in color_map.items():
            if cl in cell_types_dict:
                cell_type_to_color[cell_types_dict[cl]] = color
        adata.uns[f"{new_key}_colors"] = [cell_type_to_color[ct] for ct in adata.obs[new_key].cat.categories]

    elif f"{cluster_key}_colors" in adata.uns:
        # Copy colors from cluster_key to new_key
        cluster_colors = adata.uns[f"{cluster_key}_colors"]
        clusters = adata.obs[cluster_key].cat.categories

        # Map cell types to the corresponding cluster colors
        cell_type_to_color = {
            cell_types_dict[cl]: color
            for cl, color in zip(clusters, cluster_colors)
            if cl in cell_types_dict
        }

        # Assign colors to new cell_type column
        adata.uns[f"{new_key}_colors"] = [
            cell_type_to_color[ct] for ct in adata.obs[new_key].cat.categories
        ]

    return adata

def update_obs_colors(adata, group_key, color_dict, new_key=None):
    """
    Assign colors to an AnnData obs column based on a mapping dictionary.

    Parameters
    ----------
    adata : AnnData
        The AnnData object.
    group_key : str
        Name of the obs column containing group labels (clusters or cell types).
    color_dict : dict
        Dictionary mapping group labels to color codes (hex or named colors).
    new_key : str, optional
        Name of the color entry in adata.uns. Defaults to f"{group_key}_colors".

    Returns
    -------
    AnnData
        Same AnnData object with colors added to adata.uns.

    Raises
    ------
    KeyError
        If `group_key` is not in adata.obs.
    ValueError
        If any group in adata.obs[group_key] does not exist in `color_dict`.
    """
    import pandas as pd

    if group_key not in adata.obs:
        raise KeyError(f"{group_key} not found in adata.obs")

    if new_key is None:
        new_key = f"{group_key}_colors"

    groups = adata.obs[group_key]

    # Check for missing groups
    missing = set(groups.unique()) - set(color_dict.keys())
    if missing:
        raise ValueError(f"The following groups have no color assigned: {missing}")

    

    # Save colors to adata.uns for plotting
    # Ensure the order matches categories if it's a categorical column
    if pd.api.types.is_categorical_dtype(groups):
        adata.uns[new_key] = [color_dict[cat] for cat in groups.cat.categories]
    else:
        # Map groups to colors
        color_series = groups.map(color_dict)
        adata.uns[new_key] = color_series.unique().tolist()

    return adata 

def generate_distinct_shades(base_hex, n):
    """
    Generate n visually distinct shades based on a base color.
    Adjusts hue and saturation for better spread.
    """
    base_rgb = np.array(mcolors.to_rgb(base_hex))
    h, l, s = colorsys.rgb_to_hls(*base_rgb)
    
    # Spread hues +/- 15 degrees from base hue for multiple colors
    hue_spread = 0.15  # fraction of 1
    if n > 1:
        hues = np.linspace(max(0, h - hue_spread), min(1, h + hue_spread), n)
    else:
        hues = [h]
    
    # Vary lightness slightly for more distinction
    lightness_vals = np.linspace(max(0.4, l - 0.2), min(0.7, l + 0.2), n)
    
    # Keep saturation slightly varied
    saturation_vals = np.linspace(max(0.6, s - 0.2), min(0.9, s + 0.1), n)
    
    shades = [mcolors.to_hex(colorsys.hls_to_rgb(hue, li, sat)) 
              for hue, li, sat in zip(hues, lightness_vals, saturation_vals)]
    return shades


def scale_values(row):
    total_sum = sum(row)
    scaled_values = {key: value * 100 / total_sum for key, value in row.items()}
    return scaled_values


def compute_cluster_mean_profiles_adata(adata, cluster_key="leiden", n_cells_per_cluster=10,
                                        random_state=None):
    """
    Compute mean expression profiles per cluster as pseudo-cells,
    returning an AnnData object with cluster info in obs.

    Parameters:
        adata: AnnData object
        cluster_key: str, column in adata.obs that defines clusters
        n_cells_per_cluster: int, number of pseudo-cells per cluster
        random_state: int, optional random seed for reproducibility

    Returns:
        pseudo_adata: AnnData object
            Rows = pseudo-cells, Columns = genes
            obs includes:
                - cluster_key: cluster assignment
                - subset_id: index of the subset within the cluster
    """
    if random_state is not None:
        np.random.seed(random_state)

    mean_profiles = []
    obs_records = []

    for cluster in adata.obs[cluster_key].unique():
        # Subset cells in the current cluster
        adata_cluster = adata[adata.obs[cluster_key] == cluster]
        n_subsets = min(n_cells_per_cluster, adata_cluster.n_obs)

        if n_subsets == 0:
            print(f"Cluster {cluster} has no cells!")
            continue

        # Shuffle and split indices
        idx = np.random.permutation(adata_cluster.n_obs)
        split_indices = np.array_split(idx, n_subsets)

        for rid, sub_idx in enumerate(split_indices):
            if len(sub_idx) == 0:
                continue
            # Mean expression for this subset
            mean_expr = adata_cluster[sub_idx].X.mean(axis=0)
            mean_profiles.append(np.asarray(mean_expr).flatten())

            # Record metadata
            obs_records.append({
                cluster_key: cluster,
                "subset_id": rid,
                "n_cells": len(sub_idx)
            })

    # Create AnnData
    pseudo_adata = anndata.AnnData(
        X=np.array(mean_profiles),
        obs=pd.DataFrame(obs_records),
        var=pd.DataFrame(index=adata.var_names)
    )

    return pseudo_adata


def summary_heatmap(df, cluster_key='leiden', experiment=None, save_path=None):
    """
    Create a heatmap of mean expression profiles, grouped by clusters.
    You can get this data from compute_cluster_mean_profiles.

    Parameters:
        df: pd.DataFrame
            Rows = pseudo-cells (e.g., 'leiden0#0'), Columns = genes
        cluster_key: str
            Cluster prefix in row names (default='leiden')
        experiment: str, optional
            Used for default save filename if save_path is None
        save_path: str, optional
            File path to save the heatmap

    Returns:
        None (plots the heatmap and optionally saves to file)
    """

    # Extract cluster IDs from row names
    leiden_ids = df.index.str.extract(f"{cluster_key}(\\d+)").iloc[:, 0].astype(int)
    
    # Sort pseudo-cells by cluster
    sorted_idx = leiden_ids.argsort()
    df_sorted = df.iloc[sorted_idx]
    leiden_sorted = leiden_ids.iloc[sorted_idx]

    # Create color palette
    unique_clusters = sorted(leiden_sorted.unique())
    palette = sns.color_palette("hls", len(unique_clusters))
    cluster_colors = dict(zip(unique_clusters, palette))
    row_colors = [cluster_colors[i] for i in leiden_sorted]

    # Convert to AnnData
    adata_pseudo = anndata.AnnData(df_sorted)
    adata_pseudo.obs[cluster_key] = pd.Categorical(leiden_sorted.values, ordered=True)
    adata_pseudo.uns[f"{cluster_key}_colors"] = [cluster_colors[i] for i in unique_clusters]

    # Determine save filename
    if save_path is None and experiment is not None:
        save_path = f"_{experiment}_{cluster_key}_mean_cells.svg"

    # Plot heatmap
    scanpy.pl.heatmap(
        adata_pseudo,
        var_names=adata_pseudo.var_names,
        swap_axes=True,
        show_gene_labels=True,
        groupby=cluster_key,
        cmap='viridis',
        standard_scale='var',
        show=True,
        save=save_path
    )


def plot_dotplot_matrix_legend(adata, value_col, x_col, y_col, color_map,
                               scale_func=None, log_transform=False,
                               s_min=30, s_max=300, figsize=(14,8),
                               title=None, save_path=None, legend_steps=4):
    """
    Dot plot with dot size legend.
    
    Parameters:
        adata: AnnData object
        value_col: column to aggregate
        x_col: column for x-axis (cell types)
        y_col: column for y-axis (clusters)
        color_map: dict, mapping x_col values to hex colors
        scale_func: callable, optional row-wise scaling
        log_transform: bool, apply log1p to counts
        s_min, s_max: min/max dot size
        figsize: figure size
        title: optional plot title
        save_path: optional path to save figure
        legend_steps: number of dots to show in size legend
    """
    # Pivot table
    pivot_tab = adata.obs.pivot_table(values=value_col, index=y_col, columns=x_col,
                                     aggfunc='count', fill_value=0)

    # Row-wise scaling
    if scale_func:
        pivot_tab = pivot_tab.apply(scale_func, axis=1, result_type='expand')

    # Log transform
    if log_transform:
        pivot_tab = np.log1p(pivot_tab)

    # Normalize dot sizes
    val_min = pivot_tab.values.min()
    val_max = pivot_tab.values.max()
    if val_max == val_min:
        dot_sizes = np.full(pivot_tab.shape, (s_min + s_max) / 2)
    else:
        dot_sizes = s_min + (pivot_tab.values - val_min) / (val_max - val_min) * (s_max - s_min)

    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    for i, cluster in enumerate(pivot_tab.index):
        for j, cell_type in enumerate(pivot_tab.columns):
            ax.scatter(j, i, s=dot_sizes[i, j],
                       color=color_map.get(cell_type, "#7f7f7f"),
                       edgecolor="k", alpha=0.8, linewidth=0.5)

    # Axes
    ax.set_xticks(range(len(pivot_tab.columns)))
    ax.set_xticklabels(pivot_tab.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(pivot_tab.index)))
    ax.set_yticklabels(pivot_tab.index)
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    if title:
        ax.set_title(title)

    # Create dot size legend
    if val_max != val_min:
        # select representative values for legend
        legend_values = np.linspace(val_min, val_max, legend_steps)
        legend_sizes = s_min + (legend_values - val_min) / (val_max - val_min) * (s_max - s_min)
        legend_dots = [Line2D([0], [0], marker='o', color='w',
                              markerfacecolor='gray', markersize=np.sqrt(s), label=f"{v:.1f}")
                       for s, v in zip(legend_sizes, legend_values)]
        ax.legend(handles=legend_dots, title=value_col, loc='upper right', bbox_to_anchor=(1.2, 1))
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, bbox_inches='tight', format='svg')
    plt.show()


def plot_stacked_bar(adata, value_col, index_col, label_col, color_map, scale_func=None,
                     figsize=(10,6), title=None, save_path=None, reverse_legend=True):
    """
    Create a stacked bar plot with correctly aligned legend.
    
    Parameters:
        adata: AnnData object
        value_col: str, column to aggregate
        index_col: str, column for x-axis
        label_col: str, column for stacked bars
        color_map: dict, mapping of label -> hex color
        scale_func: callable, optional function to scale values (applied per row)
        figsize: tuple, figure size
        title: str, optional title
        save_path: str, optional path to save figure
        reverse_legend: bool, if True, reverse legend to match stack order
    """
    # Pivot table
    pivot_tab = adata.obs.pivot_table(values=value_col, index=index_col, columns=label_col, aggfunc='count')
    
    # Apply scaling if needed
    if scale_func is not None:
        pivot_tab = pivot_tab.apply(scale_func, axis=1, result_type='expand')
    
    # Determine plotting order
    plot_columns = pivot_tab.columns
    plot_colors = [color_map.get(col, "#7f7f7f") for col in plot_columns]
    
    if reverse_legend:
        # Reverse columns and colors to match stacking order
        plot_columns = plot_columns[::-1]
        plot_colors = plot_colors[::-1]
        pivot_tab = pivot_tab[plot_columns]
    
    # Plot
    ax = pivot_tab.plot.bar(stacked=True, color=plot_colors, figsize=figsize)
    
    # Legend outside
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    # Titles & labels
    if title is not None:
        ax.set_title(title)
    plt.ylabel("Counts" if scale_func is None else "Scaled values")
    plt.xlabel(index_col)
    plt.tight_layout()
    
    # Save figure if needed
    if save_path is not None:
        plt.savefig(save_path, format='svg', bbox_inches='tight')
    
    plt.show()

def plot_umap_with_cluster_contours(
    adata,
    contour_group='leiden',
    subset_group="CellType",
    group_key_selected="T_cells",
    basis='X_umap',
    cells_unselected='lightgray',
    fill_alpha=0.2,
    contour_alpha=0.8,
    contour_linewidth=1.5
):
    """
    Plots UMAP with cluster density contours and highlights a selected cell type.
    
    Parameters:
    - adata: AnnData object
    - contour_group: str, column in adata.obs used for contour grouping (e.g., 'louvain')
    - subset_group: str, column in adata.obs used for subset highlighting (e.g., 'CellType')
    - cells_selected: str, value in subset_group to highlight (e.g., 'T_cells')
    - group_key_selected: str, color for background cells
    - fill_alpha: float, alpha transparency for contour fill
    - contour_alpha: float, alpha transparency for contour edges
    - contour_linewidth: float, line width of contour edges
    """
    if basis not in adata.obsm:
        raise KeyError(
            f"UMAP coordinates '{basis}' not found in adata.obsm. "
            f"Available keys: {list(adata.obsm.keys())}. "
            "Ensure that you have run UMAP and stored the result in adata.obsm."
        )

    x = adata.obsm[basis][:, 0]
    y = adata.obsm[basis][:, 1]

    if subset_group not in adata.obs:
        raise KeyError(
            f"obs key '{subset_group}' not found in adata.obs. "
            f"Available keys: {list(adata.obs.keys())}. "
        )

    if group_key_selected not in adata.obs.keys():
        raise KeyError(
            f"obs group '{subset_group}' key does not exists '{group_key_selected}'"
            f"Available keys: {list(adata.obs[subset_group].unique())}. "
        )


    fig, ax = plt.subplots(figsize=(8,8))

    # Plot all cells in background color
    ax.scatter(x, y, s=5, color=cells_unselected, alpha=0.5)

    # Highlight selected cells
    mask_selected = adata.obs[subset_group] == group_key_selected
    ax.scatter(x[mask_selected], y[mask_selected], s=10, color='red', label=group_key_selected)

    # Get cluster information
    clusters = adata.obs[contour_group].cat.categories
    cluster_colors = adata.uns[f'{contour_group}_colors']

    # Draw filled contours for each cluster
    for cluster, color in zip(clusters, cluster_colors):
        idx = adata.obs[contour_group] == cluster
        if np.sum(idx) < 5:
            continue  # skip tiny clusters
        
        sns.kdeplot(
            data={'x': x[idx], 'y': y[idx]},
            x='x',
            y='y',
            bw_adjust=0.5,         # thighter line around the cluster
            levels=[0.01],         # only outer contour
            fill=False,
            thresh=0.01,
            #alpha=fill_alpha,
            color=color,
            linewidths=contour_linewidth,
            ax=ax
        )

    plt.legend()
    plt.title(f"Highlight: {group_key_selected} with {contour_group} contours")
    plt.xlabel('UMAP1')
    plt.ylabel('UMAP2')
    plt.axis('off')
    plt.tight_layout()
    plt.show()
    return fig, ax

def preprocess( adata, drop_ribosomal=True, drop_mitochondrial=True, minUMI=1000, minFeatures=100, n_highly_significant=300, dimensions=2, resolution=1,
    mt_regex=r'^[Mm][Tt]-',     # Default regex for mitochondrial genes (e.g., MT-CO1)
    rp_regex=r'^R[pP][SsLl]',   # Default regex for ribosomal genes (e.g., RPS12, RPL34)
):
    """
    This function provides a complete, ready-to-use preprocessing workflow for 
    single-cell RNA-seq analysis using Scanpy and Scvelo. It is designed to enforce 
    reasonable defaults that reflect personal best practices, while still allowing 
    customization where appropriate.

    This function:
    1. Calculates quality-control (QC) metrics.
    2. Filters out ribosomal and/or mitochondrial genes (if requested).
    3. Removes low-quality cells (based on minUMI and minFeatures thresholds).
    4. Normalizes counts and removes lowly expressed genes.
    5. Identifies highly variable genes.
    6. Computes neighbors, UMAP embeddings, and Leiden clusters.

    Parameters
    ----------
    adata : AnnData
        The annotated data matrix of shape n_obs × n_vars. 
        Should contain raw counts in `adata.X`.

    drop_ribosomal : bool, optional (default: True)
        If True, filters out ribosomal genes matching `rp_regex`.

    drop_mitochondrial : bool, optional (default: True)
        If True, filters out mitochondrial genes matching `mt_regex`.

    minUMI : int, optional (default: 1000)
        Minimum number of UMIs per cell. Cells below this threshold will be discarded.
        (Think of it as: "No cheap, low-quality cells allowed!")

    minFeatures : int, optional (default: 100)
        Minimum number of detected genes per cell. Cells with fewer genes will be filtered out.
        (Because a cell with 5 genes? Seriously?)

    n_highly_significant : int, optional (default: 300)
        Number of highly variable genes to keep. The more, the merrier—but also slower.
        If the umap looks speckled - use more!

    dimensions : int, optional (default: 2)
        Dimensionality of UMAP embedding. Usually 2 (for plotting) or 3 (for 3D adventures).

    resolution : float, optional (default: 1)
        Resolution parameter for Leiden clustering. Higher values = more granular clusters.
        (Like zooming in on subpopulations.)

    mt_regex : str, optional
        Regular expression to match mitochondrial genes. Default matches "MT-*" style names.

    rp_regex : str, optional
        Regular expression to match ribosomal genes. Default matches "RPS*" or "RPL*" style names.

    Returns
    -------
    adata : AnnData
        The processed AnnData object with:
        - QC metrics (`adata.obs`)
        - Filtered counts
        - UMAP embedding (`adata.obsm['X_umap']`)
        - Leiden clusters (`adata.obs['leiden']`)

    Notes
    -----
    - This function is **opinionated**. It's designed for "reasonably nice" single-cell datasets.
    - If your data is garbage, this will simply throw it away without regret.
    - Regex filtering is optional, because sometimes ribosomal and mitochondrial genes are important!
    """

    # 1. Check AnnData integrity
    if adata.X is None:
        raise ValueError("`adata.X` is empty! Did you load counts into your AnnData object?")

    # 2. Compute QC metrics (n_genes_by_counts, total_counts, pct_counts_mt, etc.)
    scanpy.pp.calculate_qc_metrics(adata, inplace=True)
    print("[INFO] QC metrics calculated.")

    # 3. Compile regex for gene filtering
    MT = re.compile(mt_regex)
    RP = re.compile(rp_regex)

    def keep_gene(gene):
        """Return True if gene passes filtering rules (not MT/RP if drop flags are set)."""
        if drop_mitochondrial and MT.match(gene):
            return False
        if drop_ribosomal and RP.match(gene):
            return False
        return True

    # 4. Apply gene-level filtering
    genes_to_keep = [g for g in adata.var_names if keep_gene(g)]
    adata = adata[:, genes_to_keep]
    print(f"[INFO] Retained {len(genes_to_keep)} genes after MT/RP filtering.")

    # 5. Add geneSum metric for UMI counts
    adata.obs['geneSum'] = adata.X.sum(axis=1).A1 if hasattr(adata.X, 'A') else adata.X.sum(axis=1)

    # 6. Filter cells based on UMI and gene count thresholds
    adata = adata[adata.obs['geneSum'] > minUMI]
    scanpy.pp.filter_cells(adata, min_genes=minFeatures)
    print(f"[INFO] Remaining cells after filtering: {adata.n_obs}")

    # 7. Normalize and log-transform
    scanpy.pp.downsample_counts(adata, counts_per_cell=minUMI)
    scv.pp.log1p(adata)
    print("[INFO] Normalization & log-transformation done.")

    # 8. Drop genes with zero counts post-filtering
    scanpy.pp.filter_genes(adata, min_counts=1)

    # 9. Identify highly variable genes
    scanpy.pp.highly_variable_genes(adata, n_top_genes=n_highly_significant)
    print(f"[INFO] Selected {n_highly_significant} highly variable genes.")

    # 10. Compute neighbors & UMAP embedding
    scanpy.pp.neighbors(adata)
    scanpy.tl.umap(adata, n_components=dimensions)
    print(f"[INFO] UMAP computed with {dimensions} dimensions.")

    # 11. Cluster using Leiden
    scanpy.tl.leiden(adata, resolution=resolution)
    print(f"[INFO] Leiden clustering done (resolution={resolution}).")

    return adata


def gene_coexpression_store(ad, gene, subset_cells=None, top_n=10, store_in_var=True):
    """
    Compute co-expressed genes with a given gene in a Scanpy AnnData object,
    and optionally store the results in ad.uns and ad.var.

    Parameters
    ----------
    ad : AnnData
        The single-cell dataset (log-normalized recommended).
    gene : str
        Gene of interest (e.g., "NCAM1").
    subset_cells : str or list, optional
        Column name in ad.obs to filter cells (e.g., 'cell_type') 
        or boolean array/list to subset cells.
    top_n : int
        Number of top positively co-expressed genes to return.
    store_in_var : bool
        Whether to store correlations in ad.var['coexpression_with_<gene>'].

    Returns
    -------
    corr_df : pd.DataFrame
        DataFrame with genes and their Pearson correlation to the target gene, sorted descending.
    """
    # Optionally subset cells
    if subset_cells is not None:
        if isinstance(subset_cells, str):
            if subset_cells not in ad.obs.columns:
                raise ValueError(f"{subset_cells} not found in ad.obs")
            ad_sub = ad[ad.obs[subset_cells].astype(bool)].copy()
        else:  # assume boolean list/array
            ad_sub = ad[subset_cells].copy()
    else:
        ad_sub = ad

    # Ensure the gene exists
    if gene not in ad_sub.var_names:
        raise ValueError(f"{gene} not found in dataset")

    # Expression matrix as dense numpy array
    X = ad_sub.X.toarray() if hasattr(ad_sub.X, "toarray") else ad_sub.X

    # Gene index
    idx = ad_sub.var_names.get_loc(gene)

    # Vectorized Pearson correlation
    x = X[:, idx]
    x_centered = x - x.mean()
    X_centered = X - X.mean(axis=0)
    corr = (X_centered.T @ x_centered) / (np.sqrt((X_centered**2).sum(axis=0)) * np.sqrt((x_centered**2).sum()))
    
    corr_df = pd.DataFrame({'gene': ad_sub.var_names, 'pearson_corr': corr})
    corr_df = corr_df.sort_values('pearson_corr', ascending=False).reset_index(drop=True)

    # Store results in AnnData
    ad.uns.setdefault('coexpression', {})
    ad.uns['coexpression'][gene] = corr_df

    if store_in_var:
        col_name = f'coexpression_with_{gene}'
        ad.var[col_name] = corr_df.set_index('gene').reindex(ad.var_names)['pearson_corr']

    return corr_df.head(top_n)


def write_pseudobulk_stats_tables( adata, key_added, cluster_key="leiden", min_pseudosamples=10, cells_per_pseudo=50, agg="mean", random_state=42 ):
    """
    Creates pseudo-bulk samples per cluster, performs differential expression,
    writes tables for each cluster, and stores the results in the original adata.uns.
    """
    # --- Step 1: create pseudobulk ---
    pseudo_adata = create_pseudobulk(
        adata,
        cluster_key=cluster_key,
        min_pseudosamples=min_pseudosamples,
        cells_per_pseudo=cells_per_pseudo,
        agg=agg,
        random_state=random_state
    )

    # --- Step 2: run DE on pseudo-bulk ---
    scanpy.tl.rank_genes_groups(
        pseudo_adata,
        groupby=cluster_key,
        key_added=key_added,
        method='wilcoxon',
    )

    # --- Step 3: save overview plot ---
    scanpy.pl.rank_genes_groups(pseudo_adata, key=key_added, save=f"_{key_added}_overview.svg")

    # --- Step 4: write tables per cluster ---
    diff_results = pseudo_adata.uns[key_added]
    columns = ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']
    
    try:
        os.makedirs(f"{key_added}", exist_ok=True)
    except OSError as e:
        print(e)
    
    safe_key = re.sub(r'[<>:"/\\|?*\s]', '_', key_added)
    for cluster in pseudo_adata.obs[cluster_key].unique():
        table = {}
        for col in columns:
            table[col] = pd.DataFrame(diff_results[col])[str(cluster)]
        safe_cluster = re.sub(r'[<>:"/\\|?*\s]', '_', cluster)
        safe_cn = re.sub(r'[<>:"/\\|?*\s]', '_', cluster_key)
        table_df = pd.DataFrame(table)
        table_df.to_csv(f"{safe_key}/{safe_cn}_cluster_{safe_cluster}.csv")

    # --- Step 5: store DE results in original adata ---
    if "pseudo_bulk_DE" not in adata.uns:
        adata.uns["pseudo_bulk_DE"] = {}
    adata.uns["pseudo_bulk_DE"][key_added] = pseudo_adata.uns[key_added]

    print(f"Pseudo-bulk DE analysis complete. Results stored under adata.uns['pseudo_bulk_DE']['{key_added}'].")

def create_pseudobulk(adata, cluster_key='leiden', min_pseudosamples=10, cells_per_pseudo=50, agg='mean', random_state=None):
    """
    Create pseudo-bulk samples from single-cell data.
    
    Parameters
    ----------
    adata : AnnData
        Single-cell AnnData object.
    cluster_key : str
        Column in adata.obs with cluster labels (e.g., 'leiden').
    min_pseudosamples : int
        Minimum number of pseudo-samples per cluster.
    cells_per_pseudo : int
        Target number of cells per pseudo-sample.
    agg : str
        Aggregation method: 'mean' or 'sum'.
    random_state : int or None
        Random seed for reproducibility.
        
    Returns
    -------
    pseudo_adata : AnnData
        AnnData object with pseudo-bulk samples as observations.
    """
    if random_state is not None:
        np.random.seed(random_state)
    
    clusters = adata.obs[cluster_key].unique()
    pseudo_data_list = []
    pseudo_obs = []
    
    for cl in clusters:
        cells = adata.obs_names[adata.obs[cluster_key] == cl]
        n_cells = len(cells)
        
        # Determine number of pseudo-samples
        n_pseudo = max(min_pseudosamples, n_cells // cells_per_pseudo)
        n_pseudo = min(n_pseudo, n_cells)  # cannot exceed number of cells

        # Skip cluster if it cannot generate at least 3 pseudo-samples
        if n_pseudo < 3:
            print(f"Skipping cluster {cl} (only {n_pseudo} pseudo-samples possible)")
            continue
        
        # Shuffle cells
        shuffled = np.random.permutation(cells)
        
        # Split cells into pseudo-samples
        pseudo_samples = np.array_split(shuffled, n_pseudo)
        
        for i, pseudo in enumerate(pseudo_samples):
            # Aggregate expression
            X_subset = adata[pseudo].X
            if agg == 'mean':
                mean_counts = X_subset.mean(axis=0)
            elif agg == 'sum':
                mean_counts = X_subset.sum(axis=0)
            else:
                raise ValueError("agg must be 'mean' or 'sum'")
            
            # Handle sparse matrices
            if hasattr(mean_counts, "A1"):
                mean_counts = mean_counts.A1
            
            pseudo_data_list.append(mean_counts)
            pseudo_obs.append(f"{cl}")
    
    pseudo_matrix = np.vstack(pseudo_data_list)
    pseudo_adata = scanpy.AnnData(X=pseudo_matrix, obs=pd.DataFrame({cluster_key: pseudo_obs}), var=adata.var.copy())
    
    return pseudo_adata

def addSampleDataRust( adata, sname, sampleInfo ):
    """Add the sample information obtained by either running demux10x or quantifyRhapsody to an anndata object"""
    adata.obs['sname'] = sname
    adata.obs_names = [ sname +"."+ n for n in adata.obs_names]
    sampleInfo = pd.read_csv( sampleInfo, sep="\t",index_col= 0)
    sampleInfo['CellIDs'] = [ sname +"."+n+"-1" for n in sampleInfo.index.values]
    sampleInfo.set_index('CellIDs', inplace=True)
    adata.obs['CellIDs'] = adata.obs_names
    obs = adata.obs.merge( sampleInfo, left_index=True, right_index=True)
    print(obs.shape)
    obs.drop_duplicates()
    obs.drop( ['CellIDs'], axis=1, inplace=True)
    print(obs.shape)
    adata = adata[ obs.index.values]
    adata.obs = obs
    ADT = pd.read_csv( 'Cell2Sample.ADT'+sname+'.R1_001.fastq.gz.tsv', sep="\t",index_col= 0)
    ADT['CellIDs'] = [ sname+"."+n+"-1" for n in ADT.index.values]
    ADT.set_index( 'CellIDs' , inplace=True)
    
    obs = adata.obs.merge( ADT, left_index=True, right_index=True )
    #print(obs)
    #print (adata)
    adata = adata[ obs.index.values]
    #print(adata)
    adata.obs = obs
    return (adata)

def testRP(x, RP):
    r= True
    if RP.match(x):
        r = False
    return (r)

def dropRP ( adata, drop = True ):
    """ remove all genes matching to /^R[pP][SsLl]/ """
    RP = re.compile('^R[pP][SsLl]')
    RPgenes = [ x  for x in adata.var.index if not testRP(x, RP)]
    print ( Counter(RPgenes) )
    if  len(RPgenes) > 0:
        adata.obs['RP[LS]sum'] = adata[:,RPgenes].X.sum(axis=1)
    else:
        adata.obs['RP[LS]sum'] = 0
    if drop and len(RPgenes) > 0:
        OK = [ testRP(x, RP) for x in adata.var.index ]
        print(Counter(OK))
        adata._inplace_subset_var( np.array(OK) )
    adata

def testMT(x, MT):
    r= True
    if MT.match(x):
        r = False
    return (r)

def dropMT (adata, drop = True ) :
    """ remove all genes matching to /^[Mm][Tt]-/ """
    MT = re.compile('^[Mm][Tt]-')
    MTgenes = [ x  for x in adata.var.index if not testMT(x, MT)]
    print ( Counter(MTgenes) )
    if len(MTgenes) > 0:
        adata.obs['MTsum'] = adata[:,MTgenes].X.sum(axis=1)
    else:
        adata.obs['MTsum'] = 0
    if drop and len(MTgenes) > 0:
        OK = [ testMT(x, MT) for x in adata.var.index ]
        adata._inplace_subset_var( np.array(OK) )


def testGene(x, MT, RP):
    r= True
    if MT.match(x) or RP.match(x):
        r = False
    return (r)

    
# Example usage:
# Assuming you have an AnnData object `adata` and have run `scanpy.tl.rank_genes_groups` already:
# write_top_genes_per_cluster(adata, n_top_genes=20, output_file="top_genes.csv")
def write_top_genes_per_cluster(adata, stats_name, n_top_genes=20, output_file="top_genes_per_cluster.csv"):
    """
    Extracts the top `n_top_genes` highly expressed genes per cluster from a Scanpy `AnnData` object
    and writes them into a CSV file. Each row in the file contains the top genes for a cluster,
    comma-separated.

    Parameters:
        adata (AnnData): The annotated data matrix from Scanpy containing gene expression data.
        n_top_genes (int): The number of top genes to extract per cluster. Default is 20.
        output_file (str): The output CSV file where top genes per cluster will be saved.
    """


    # Check if the `rank_genes_groups` function has been run
    if stats_name in adata.uns:
        ranked_genes = pd.DataFrame(adata.uns[ stats_name ]['names']).head(n_top_genes)
    elif 'pseudo_bulk_DE' in adata.uns and stats_name in adata.uns["pseudo_bulk_DE"]:
        ranked_genes = pd.DataFrame(adata.uns["pseudo_bulk_DE"][ stats_name ]['names']).head(n_top_genes)
    else:
        raise ValueError(
           f"No {stats_name} results found in the AnnData object. "
            "Please run `scanpy.tl.rank_genes_groups()` or pseudobulk DE first."
        )
    ret = {}
    with open(output_file, 'w') as f:
        # Write top genes for each cluster
        for cluster in ranked_genes.columns:
            top_genes = ranked_genes[cluster].tolist()
            ret[cluster] = top_genes
            f.write(f"{cluster}: {','.join(top_genes)}\n")

    print(f"Top {n_top_genes} genes per cluster have been written to {output_file}")
    return(ret)


def write_stats_tables( adata, key_added,cn ="louvain" ):
    """
    Creates the stats using the 'wilcoxon' method, plots the top genes using scanpy.pl.rank_genes_groups
    and finally writes tab separated tables to a <key_added> foilder.
    """ 

    scanpy.tl.rank_genes_groups(
        adata, 
        groupby   = cn,
        key_added = key_added,
        method    = 'wilcoxon',
    )

    scanpy.pl.rank_genes_groups(adata, key = key_added, save= key_added+"_overview.svg")

    diff_results = adata.uns[key_added]
    columns = ['names', 'scores', 'pvals', 'pvals_adj', 'logfoldchanges']
    try: 
        os.mkdir(f"{key_added}") 
    except OSError as error: 
        print(error)  
    safe_key = re.sub(r'[<>:"/\\|?*\s]', '_', key_added)
    for i in adata.obs[cn].unique():
        table = {}
        for j in columns:
            #print(f"Analyszing column {diff_results[j]} {i}")
            table[j] = pd.DataFrame(diff_results[j])[str(i)]

        safe_cn = re.sub(r'[<>:"/\\|?*\s]', '_', cn)
        safe_i = re.sub(r'[<>:"/\\|?*\s]', '_', i) or "__"
        table = pd.DataFrame(table).to_csv(f"{safe_key}/{safe_cn}_cluster_{safe_i}.csv")



def remove_gzipped_files(path):
    # Ensure the directory exists
    if os.path.exists(path):
        for file_name in os.listdir(path):
            # Check if the file is a gzipped file
            if file_name.endswith(".gz"):
                gzipped_file_path = os.path.join(path, file_name)
                os.remove(gzipped_file_path)
                print(f"Removed gzipped file: {gzipped_file_path}")
    else:
        print(f"The directory {path} does not exist.")

def gzip_files_in_directory(path):
    # Ensure the directory exists
    if os.path.exists(path):
        # Run the gzip command on all files in the directory
        subprocess.run(f"gzip {path}/*", shell=True, check=True)
        print(f"All files in {path} have been gzipped.")
    else:
        print(f"The directory {path} does not exist.")

def write_data_to_directory(adata, output_dir):
    """ 
    exports the anndata.X sparse matrix as 10X formated matrix.mtx barcodes.tsv and features.tsv set 
    """
    # Ensure the directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Remove gzipped files before writing new ones
    remove_gzipped_files(output_dir)

    # Write the matrix data to a .mtx file
    from scipy.io import mmwrite
    mmwrite(os.path.join(output_dir, 'matrix.mtx'), adata.X.T )

    # Write the barcodes (cell identifiers) to a .tsv file
    with open(os.path.join(output_dir, 'barcodes.tsv'), 'w') as f:
        for barcode in adata.obs_names.values:
            f.write(f"{barcode}\n")

    # Write the features (gene identifiers) to a .tsv file
    # Create the dataframe for the features.tsv.gz file
    features_df = pd.DataFrame({
        'gene_id': adata.var.index,  # Gene identifiers (e.g., Ensembl IDs)
        'gene_name': adata.var['gene_name'] if 'gene_name' in adata.var.columns else adata.var.index,  # Gene names (or use gene_id if unavailable)
        'gene_type': adata.var['feature_types'] if 'feature_types' in adata.var.columns else  ['Gene Expression'] * len(adata.var)  # Gene type, usually 'Gene Expression'
    })
    features_df.to_csv(os.path.join(output_dir, 'features.tsv'), sep='\t', index=False, header=False )

    # Optionally, gzip the files in the directory
    gzip_files_in_directory(output_dir)


def impute_expression_data( adata, output_dir ):
    """
    Uses scvi.model.SCVI to impute the expression data and adds that (model.posterior_predictive_sample) as the adata.X sparse matrix.
    In addition the imputed data is saved as 10X matrix.mtx, features.tsv and barcodes.tsv file triplet.
    """

    # Step 1: Get the number of CPU cores
    num_cores = os.cpu_count()
    
    # Step 2: Set num_workers to the number of cores + 1
    num_workers = num_cores + 1
    
    # Convert to CSR format
    adata.X = scipy.sparse.csr_matrix(adata.X)
    # Set up the model
    scvi.model.SCVI.setup_anndata(adata)
    model = scvi.model.SCVI(adata)
    # Train the model
    model.train()
    # Get the imputed expression values for cells
    adata.X = model.posterior_predictive_sample(adata, n_samples=1).to_scipy_sparse()
    write_data_to_directory( adata, output_dir )


def scatter_safe(adata, color, **kwargs):
    """
    Wrapper around scvelo.pl.scatter that:
    1. Skips genes that are not in the dataset.
    2. Skips genes with zero expression in all cells.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    color : list or str
        Gene(s) or feature(s) to plot.
    **kwargs : dict
        Other keyword arguments passed to scv.pl.scatter.
    """
    # Convert single string to list for uniform processing
    if isinstance(color, str):
        color = [color]
    
    color_filtered = []
    for gene in color:
        if gene in adata.obs:
            # Always keep obs features
            color_filtered.append(gene)
        elif gene in adata.var_names:
            # Check if gene has any nonzero expression
            Xcol = adata[:, gene].X
            if hasattr(Xcol, "toarray"):  # sparse
                nonzero = np.any(Xcol.toarray())
            else:  # dense
                nonzero = np.any(Xcol)
            if nonzero:
                color_filtered.append(gene)
        # silently skip missing
    
    if not color_filtered:
        print("No valid genes/features to plot.")
        return
    
    scv.pl.scatter(adata, color=color_filtered, **kwargs)

def converge_cluster_merging( adata, group_key, cutoff=0.9, layer=None, gene_symbols=None, key_added=None, verbose=True ):
    """
    Iteratively merge clusters until no two cluster mean expression profiles
    exceed the given correlation cutoff.

    Parameters
    ----------
    adata : AnnData
        AnnData object with clustering in .obs[group_key].
    group_key : str
        Column in adata.obs with initial cluster labels.
    cutoff : float
        Correlation threshold. Clusters with corr >= cutoff are merged.
    layer : str, optional
        Which layer to use instead of adata.X.
    gene_symbols : str, optional
        Column in adata.var to use as feature index instead of var_names.
    key_added : str, optional
        Where to store the merged cluster labels. Defaults to group_key + "_merged".
    verbose : bool
        Print progress if True.

    Returns
    -------
    dict
        Mapping of old cluster labels → new sequential IDs.
    """

    if key_added is None:
        key_added = group_key + f"_merged{cutoff}"
    adata.obs[key_added] = adata.obs[group_key].astype(str).copy()


    # iterative merging
    while True:
        mean_exp = grouped_obs_mean(adata, key_added, layer, gene_symbols)
        if not mergeClosest( adata, key_added, cutoff ):
            break



def grouped_obs_mean(adata, group_key, layer=None, gene_symbols=None):
    """
    calculate the mean expression for all genes in all groups
    """
    if layer is not None:
        getX = lambda x: x.layers[layer]
    else:
        getX = lambda x: x.X
    if gene_symbols is not None:
        new_idx = adata.var[idx]
    else:
        new_idx = adata.var_names

    grouped = adata.obs.groupby(group_key)
    out = pd.DataFrame(
        np.zeros((adata.shape[1], len(grouped)), dtype=np.float64),
        columns=list(grouped.groups.keys()),
        index=adata.var_names
    )

    for group, idx in grouped.indices.items():
        X = getX(adata[idx])
        out[group] = np.ravel(X.mean(axis=0, dtype=np.float64))
    return out
    reIDgroup( adata, group_key )
    adata.obs[group_key] = adata.obs[group_key+'_renamed']
    del adata.obs[group_key+'_renamed']

def mergeClosest( adata, group = 'louvain11_merged.9', cut = 0.9 ):
    """ 
    create mean expression per group, correlate all mean groups and scan all groups for corr > cut.
    The first group with cor > cut will merge with the best correlating groups available.
    Only one group will be merged at a time.
    """
    #print (str(len(Counter(adata.obs[group])))+ " -", end=" ")
    df = grouped_obs_mean( adata, group )
    df = df.corr()
    for i in range(0,df.shape[0]):
        df.at[df.columns[i],df.columns[i]] = 0
    #print("processing groups")
    #print (df)
    for i in range( 0, len(df.columns)):
        i = df.columns[i]
        col = df[i]
        #print ( "max value="+ str( max( col )) )
        if max( col ) >= cut:
            gr = np.array(adata.obs[group])
            m = [df.columns[a] for a in range(0, len(col)) if col[a] == max(col) ]
            m = m[0]
            print( f"m{i}->{m}", end=",", flush=True)
            #print ("max value="+str(max[col])+ " merge " + str(m) + " to "+ str(i)+".")
            for x in range(0,len(adata.obs[group])):
                    if str(gr[x]) == i:
                        #print ( "changed" + str(x)+ " change " + gr[x]  + " to "+ str(m))
                        gr[x] = m
            #type(gr)
            #print("finished")
            adata.obs[group] = gr
            return True
    print ( "no cluster did match to any other cluster with cor > "+ str(cut) )
    return False


def reIDgroup( adata, group ):
    """
    Re-set the cluster ids to be in the tange of 0-max cluster id
    """
    gr = np.array(adata.obs[group])
    n = 0
    for i in set(adata.obs[group]):
        for x in range(0,len(adata.obs[group])):
            if adata.obs[group][x] == i:
                gr[x] = n
        n = n +1
    adata.obs[group+'_renamed'] = gr



def predict_celltype(
    adata,
    marker_groups,
    cluster_key="leiden",
    out_csv="cluster_marker_summary.csv",
    excluded=None,
):
    """
    Predicts cluster-level cell types in an AnnData object by comparing
    expression of known marker genes against each cluster, with optional
    exclusion markers for sanity checks.

    The function computes per-cluster statistics (mean expression and
    percent-positive cells) for both reference marker groups and
    exclusion/control genes, then proposes a cell type label for each cluster
    based on the strongest evidence.

    Results are written to two CSV files:
      1. `out_csv` (cluster-level marker summary)
      2. `cluster_label_proposals.csv` (predicted cluster labels)

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix, typically containing single-cell RNA-seq counts.
        Must include cluster labels under `adata.obs[cluster_key]`.

    marker_groups : dict[str, list[str]]
        Dictionary mapping cell type names (keys) to lists of known marker genes
        (values). Example:
        ``{"T cells": ["CD3D", "CD3E"], "B cells": ["MS4A1"]}``.

    cluster_key : str, optional (default: "leiden")
        Column in `adata.obs` that contains cluster assignments.

    out_csv : str, optional (default: "cluster_marker_summary.csv")
        Output path for the detailed per-cluster marker summary.

    excluded : list[str] or None, optional
        Genes to treat as “exclusion markers” (e.g., B-cell markers when
        analyzing T-cell datasets). If these are unexpectedly expressed in a
        cluster, a warning suffix is added to the predicted label. If None, no
        exclusion markers are checked.

    Returns
    -------
    mapping_df : pandas.DataFrame
        Table of clusters and their suggested labels.
    proposals : dict[str, str]
        Dictionary mapping cluster IDs to proposed cell type labels.
    """

    # ---- Helper to access expression matrix safely (prefers adata.raw if present) ----
    use_adata = adata.raw if adata.raw is not None else adata
    var_names = list(use_adata.var_names)

    def get_gene_array(gene):
        if gene not in var_names:
            return None
        gi = var_names.index(gene)
        arr = use_adata.X[:, gi]
        if not isinstance(arr, np.ndarray):
            arr = arr.toarray().flatten()
        else:
            arr = arr.flatten()
        return arr

    # ---- Clusters ----
    clusters = list(np.unique(adata.obs[cluster_key].values))
    clusters.sort(key=lambda x: (float(x) if str(x).isdigit() else str(x)))  # sensible sort

    if os.path.exists(out_csv):
        print(f"Using existing file '{out_csv}' - remove the file if this is not wanted")
        summary_df = pd.read_csv(out_csv, dtype={"cluster": str})
    else:
        rows = []
        for cl in clusters:
            adata_sub = adata[adata.obs[cluster_key] == cl]  # subset by cluster
            row = {"cluster": cl, "n_cells": adata_sub.n_obs }
            for grp, genes in marker_groups.items():
                # make sure genes exist
                genes = [g for g in genes if g in adata_sub.var_names]
                if len(genes) == 0:
                    grp_mean = np.nan
                    grp_pct = 0.0
                else:
                    expr = adata_sub[:, genes].X  # matrix of cells x genes
                    if not isinstance(expr, np.ndarray):
                        expr = expr.toarray()
                    grp_mean = float(expr.mean())
                    grp_pct = float((expr > 0).mean() * 100)
                row[f"{grp}_mean"] = grp_mean
                row[f"{grp}_pct_pos"] = grp_pct
            rows.append(row)

        summary_df = pd.DataFrame(rows).sort_values("cluster")
        summary_df.to_csv(out_csv, index=False)
        print(f"Wrote cluster summary to: {out_csv}")

    # ---- Propose label per cluster ----
    proposals = {}
    for _, r in summary_df.iterrows():
        cl = r["cluster"]

        # Flag unexpected expression of exclusion genes
        b_flag = False
        if excluded:
            b_flag = any(
                r[f"{g}_pct_pos"] > 1.0 for g in excluded if r.get(f"{g}_exists", False)
            )

        # choose group with highest group_mean
        grp_means = {grp: r[f"{grp}_mean"] for grp in marker_groups.keys()}
        grp_means_nonan = {k: (v if not pd.isna(v) else -1e9) for k, v in grp_means.items()}
        best_grp = max(grp_means_nonan, key=grp_means_nonan.get)
        best_mean = grp_means_nonan[best_grp]

        # require some minimal evidence
        cond_pct = r[f"{best_grp}_pct_pos"] >= 5.0
        cond_mean = best_mean > 0.01
        if not (cond_pct or cond_mean):
            proposal = "Unknown_or_ambiguous"
        else:
            proposal = best_grp

        if b_flag:
            proposal += " (CHECK_B_MARKERS_PRESENT)"

        proposals[cl] = proposal

    # ---- Print suggestions ----
    print("\nSuggested labels per cluster (best guess):")
    for cl in clusters:
        print(f"Cluster {cl}: {proposals[cl]}")

    # ---- Write mapping ----
    mapping_df = pd.DataFrame(
        [{"cluster": k, "suggested_label": v} for k, v in proposals.items()]
    )
    mapping_df.to_csv("cluster_label_proposals.csv", index=False)
    print("Wrote label proposals to: cluster_label_proposals.csv")

    # ---- Add to AnnData ----
    cluster_to_label = dict(zip(mapping_df["cluster"], mapping_df["suggested_label"]))
    adata.obs["predicted_celltype"] = adata.obs[cluster_key].map(cluster_to_label)
    print("Added 'predicted_celltype' column to adata.obs")

    return mapping_df, proposals


