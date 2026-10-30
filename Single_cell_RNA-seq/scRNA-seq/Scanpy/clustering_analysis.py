#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 2: Dimensionality reduction and clustering analysis
Input: post-QC data
Function: PCA / PCA+Harmony + clustering + UMAP/tSNE
"""

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import os
import argparse
import logging
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
LOGGER = logging.getLogger(__name__)


def run_pca(
    adata,
    n_comps=50,
    n_pcs=30,
    max_value=10
):
    """
    Run PCA dimensionality reduction

    Parameters:
    -----------
    adata : AnnData
        Input AnnData object (already QC'd)
    n_comps : int
        Number of PCA components
    n_pcs : int
        Number of PCs for downstream analysis
    max_value : float
        Maximum value for scaling
    """
    LOGGER.info("Starting PCA dimensionality reduction...")

    # If HVG not yet selected, select them first
    if 'highly_variable' not in adata.var.columns:
        LOGGER.info("Selecting highly variable genes...")
        sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3')

    # Use only highly variable genes
    adata_hvg = adata[:, adata.var['highly_variable']].copy()

    # Scaling
    sc.pp.scale(adata_hvg, max_value=max_value)
    LOGGER.info(f"Scaling complete, max_value={max_value}")

    # PCA
    sc.tl.pca(adata_hvg, n_comps=n_comps)
    LOGGER.info(f"PCA complete, {n_comps} components")

    # Copy the PCA results back to the original adata
    # Note: varm and uns need not be copied; only the PCA coordinates in obsm are needed for downstream analysis
    adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']

    return adata, n_pcs


def run_harmony(
    adata,
    batch_key="batchid",
    n_pcs=30
):
    """
    Run Harmony batch correction

    Parameters:
    -----------
    adata : AnnData
        Input AnnData object (PCA already computed)
    batch_key : str
        Batch key name
    n_pcs : int
        Number of PCs for downstream analysis
    """
    LOGGER.info(f"Starting Harmony batch correction, batch_key={batch_key}...")

    try:
        import harmony_pytorch as harmony

        # Run Harmony
        sc.external.pp.harmony_integrate(
            adata,
            key=batch_key,
            basis='X_pca',
            adjusted_basis='X_pca_harmony'
        )

        LOGGER.info("Harmony batch correction complete")
        use_rep = 'X_pca_harmony'

    except ImportError:
        LOGGER.warning("harmony_pytorch not installed, using plain PCA")
        use_rep = 'X_pca'

    return adata, use_rep


def run_clustering(
    adata,
    n_pcs=30,
    n_neighbors=15,
    resolution=0.8,
    method="leiden",
    use_rep="X_pca"
):
    """
    Run clustering

    Parameters:
    -----------
    adata : AnnData
        Input AnnData object
    n_pcs : int
        Number of PCs for clustering
    n_neighbors : int
        Number of neighbors
    resolution : float
        Clustering resolution
    method : str
        Clustering method (leiden or louvain)
    use_rep : str
        Representation used to build the neighbor graph
    """
    LOGGER.info(f"Starting clustering, method={method}, resolution={resolution}...")

    # Build neighborhood graph
    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        use_rep=use_rep
    )
    LOGGER.info(f"Neighbor graph construction complete, n_neighbors={n_neighbors}, n_pcs={n_pcs}")

    # Clustering
    if method == "leiden":
        sc.tl.leiden(adata, resolution=resolution)
        # If clusters already exists, delete it first
        if 'clusters' in adata.obs.columns:
            adata.obs.drop(columns=['clusters'], inplace=True)
        adata.obs.rename(columns={"leiden": "clusters"}, inplace=True)
    elif method == "louvain":
        sc.tl.louvain(adata, resolution=resolution)
        # If clusters already exists, delete it first
        if 'clusters' in adata.obs.columns:
            adata.obs.drop(columns=['clusters'], inplace=True)
        adata.obs.rename(columns={"louvain": "clusters"}, inplace=True)

    # Ensure clusters is categorical type
    adata.obs['clusters'] = pd.Categorical(adata.obs['clusters'])

    n_clusters = adata.obs['clusters'].nunique()
    LOGGER.info(f"Clustering complete, {n_clusters} clusters in total")

    return adata


def run_embedding(
    adata,
    n_pcs=30,
    min_dist=0.5,
    use_rep="X_pca",
    embedding_method="umap"
):
    """
    Run dimensionality reduction visualization

    Parameters:
    -----------
    adata : AnnData
        Input AnnData object
    n_pcs : int
        Number of PCs for tSNE
    min_dist : float
        min_dist parameter for UMAP
    use_rep : str
        Representation used for tSNE
    embedding_method : str
        Dimensionality reduction method (umap, tsne, both)
    """
    LOGGER.info(f"Starting dimensionality reduction visualization, method={embedding_method}...")

    if embedding_method in ["umap", "both"]:
        sc.tl.umap(adata, min_dist=min_dist)
        LOGGER.info(f"UMAP complete, min_dist={min_dist}")

    if embedding_method in ["tsne", "both"]:
        sc.tl.tsne(adata, n_pcs=n_pcs, use_rep=use_rep)
        LOGGER.info("tSNE complete")

    return adata


def find_marker_genes(
    adata,
    method="wilcoxon",
    n_genes=100
):
    """
    Find marker genes for each cluster

    Parameters:
    -----------
    adata : AnnData
        Input AnnData object
    method : str
        Test method (wilcoxon, t-test, logreg)
    n_genes : int
        Number of genes returned per cluster

    Returns:
    --------
    DataFrame
        Marker gene results
    """
    LOGGER.info(f"Finding marker genes, method={method}...")

    # Run differential analysis
    sc.tl.rank_genes_groups(
        adata,
        groupby='clusters',
        method=method,
        n_genes=n_genes,
        use_raw=False
    )

    # Extract results
    result = adata.uns['rank_genes_groups']

    # Convert to DataFrame
    all_results = []
    for cluster in result['names'].dtype.names:
        cluster_data = pd.DataFrame({
            'cluster': cluster,
            'gene': result['names'][cluster],
            'scores': result['scores'][cluster],
            'logfoldchanges': result['logfoldchanges'][cluster],
            'pvals': result['pvals'][cluster],
            'pvals_adj': result['pvals_adj'][cluster],
        })
        all_results.append(cluster_data)

    marker_df = pd.concat(all_results, ignore_index=True)

    LOGGER.info(f"Found {len(marker_df)} marker genes")

    return marker_df


def add_gene_annotation(df, anno_df, gene_col="gene"):
    """
    Add annotations to genes

    Parameters:
    -----------
    df : DataFrame
        Gene data table
    anno_df : DataFrame
        Gene annotation table
    gene_col : str
        Gene name column name

    Returns:
    --------
    DataFrame
        Data table with annotations added
    """
    # Merge annotations
    df_merged = df.merge(
        anno_df,
        left_on=gene_col,
        right_on="id",
        how="left"
    )

    # Drop the duplicate id column
    if "id" in df_merged.columns:
        df_merged = df_merged.drop(columns=["id"])

    # Fill missing annotations with '--'
    df_merged = df_merged.fillna('--')

    # Reorder columns
    base_cols = ["cluster", "gene", "scores", "logfoldchanges", "pvals", "pvals_adj"]
    anno_cols = ["ensembl_id", "gene_type", "gene_description", "TFs_Family", "GO_id", "GO_term", "KEGG_id", "KEGG_description"]

    # Keep only existing columns
    ordered_cols = [col for col in base_cols + anno_cols if col in df_merged.columns]
    other_cols = [col for col in df_merged.columns if col not in ordered_cols]
    ordered_cols.extend(other_cols)

    df_merged = df_merged[ordered_cols]

    return df_merged


def plot_embeddings(
    adata,
    output_dir,
    color_by=["clusters", "sampleid"],
    embedding_method="umap"
):
    """
    Plot dimensionality reduction visualizations

    Parameters:
    -----------
    adata : AnnData
        Input AnnData object
    output_dir : str
        Output directory
    color_by : list
        Column names to color by
    embedding_method : str
        Dimensionality reduction method (umap, tsne, both)
    """
    LOGGER.info("Plotting dimensionality reduction visualizations...")

    os.makedirs(output_dir, exist_ok=True)

    # Set the scanpy figure output directory
    sc.settings.figdir = output_dir

    for color in color_by:
        if color not in adata.obs.columns:
            LOGGER.warning(f"Column {color} does not exist, skipping")
            continue

        try:
            if embedding_method in ["umap", "both"]:
                # UMAP plot
                sc.pl.umap(
                    adata,
                    color=color,
                    title=f'UMAP - {color}',
                    save=f'_umap_{color}.png',
                    show=False
                )
                LOGGER.info(f"Saved UMAP plot: umap_{color}.png")

            if embedding_method in ["tsne", "both"]:
                # tSNE plot
                sc.pl.tsne(
                    adata,
                    color=color,
                    title=f'tSNE - {color}',
                    save=f'_tsne_{color}.png',
                    show=False
                )
                LOGGER.info(f"Saved tSNE plot: tsne_{color}.png")
        except Exception as e:
            LOGGER.warning(f"Error plotting {color}: {e}, skipping")


def save_results(adata, marker_df, outdir, refgenome=None):
    """
    Save results

    Parameters:
    -----------
    adata : AnnData
        AnnData object after analysis
    marker_df : DataFrame
        Marker gene results
    outdir : str
        Output directory
    refgenome : str
        Reference genome directory path (for gene annotation)
    """
    LOGGER.info(f"Saving results to {outdir}...")

    os.makedirs(outdir, exist_ok=True)

    # Save h5ad
    adata.write_h5ad(os.path.join(outdir, "clustering_result.h5ad"))
    LOGGER.info("Saved clustering_result.h5ad")

    # Save UMAP coordinates (format: Barcode, UMAP_1, UMAP_2)
    if 'X_umap' in adata.obsm:
        umap_df = pd.DataFrame(
            adata.obsm['X_umap'],
            columns=['UMAP_1', 'UMAP_2'],
            index=adata.obs_names
        )
        umap_df.index.name = 'Barcode'
        umap_df.to_csv(os.path.join(outdir, "umap_coordinates.csv"))
        LOGGER.info("Saved umap_coordinates.csv")

    # Save clustering results (format: Barcode, sampleid, group, batchid, clusters)
    cluster_result_df = pd.DataFrame(index=adata.obs_names)
    cluster_result_df.index.name = 'Barcode'

    # Add available columns
    for col in ['sampleid', 'batchid', 'clusters']:
        if col in adata.obs.columns:
            cluster_result_df[col] = adata.obs[col]

    # If a group column exists, add it; otherwise use batchid as group
    if 'group' in adata.obs.columns:
        cluster_result_df['group'] = adata.obs['group']
    elif 'batchid' in adata.obs.columns:
        cluster_result_df['group'] = adata.obs['batchid']

    # Reorder columns
    desired_order = ['sampleid', 'group', 'batchid', 'clusters']
    cols = [col for col in desired_order if col in cluster_result_df.columns]
    cluster_result_df = cluster_result_df[cols]

    cluster_result_df.to_csv(os.path.join(outdir, "clusters_result.csv"))
    LOGGER.info("Saved clusters_result.csv")

    # Add gene annotation and save marker genes
    if refgenome and os.path.exists(os.path.join(refgenome, "gene_annotation.xls")):
        LOGGER.info(f"Loading gene annotation: {refgenome}/gene_annotation.xls")
        anno_df = pd.read_csv(os.path.join(refgenome, "gene_annotation.xls"), sep="\t")
        anno_df["id"] = anno_df["id"].astype(str)

        # Add annotation
        marker_df_annotated = add_gene_annotation(marker_df, anno_df)

        # Save all marker genes
        marker_df_annotated.to_csv(os.path.join(outdir, "marker_genes.csv"), index=False)
        LOGGER.info("Saved marker_genes.csv (with gene annotation)")

        # Save top 10 marker genes per cluster
        top_markers = marker_df_annotated.groupby("cluster").head(10)
        top_markers.to_csv(os.path.join(outdir, "top_marker_genes.csv"), index=False)
        LOGGER.info("Saved top_marker_genes.csv (with gene annotation)")
    else:
        # No annotation file; save raw data
        marker_df.to_csv(os.path.join(outdir, "marker_genes.csv"), index=False)
        LOGGER.info("Saved marker_genes.csv")

        # Save top 10 marker genes per cluster
        top_markers = marker_df.groupby("cluster").head(10)
        top_markers.to_csv(os.path.join(outdir, "top_marker_genes.csv"), index=False)
        LOGGER.info("Saved top_marker_genes.csv")

    LOGGER.info("All results saved")


def main():
    """
    Main function
    """
    parser = argparse.ArgumentParser(
        description="Step 2: Dimensionality reduction and clustering analysis (PCA / PCA+Harmony + clustering + UMAP/tSNE + marker genes)"
    )

    # Input/Output
    parser.add_argument("-i", "--input", required=True, help="Input file path (post-QC .h5ad)")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory")

    # PCA parameters
    parser.add_argument("--n_comps", type=int, default=50, help="Number of PCA components (default: 50)")
    parser.add_argument("--n_pcs", type=int, default=30, help="Number of PCs for downstream analysis (default: 30)")
    parser.add_argument("--max_value", type=float, default=10.0, help="Maximum scaling value (default: 10)")

    # Harmony parameters
    parser.add_argument("--use_harmony", action="store_true", help="Whether to apply Harmony batch correction")
    parser.add_argument("--batch_key", type=str, default="batchid", help="Batch key name (default: batchid)")

    # Clustering parameters
    parser.add_argument("--cluster_method", type=str, default="leiden", choices=["leiden", "louvain"],
                        help="Clustering method (default: leiden)")
    parser.add_argument("--resolution", type=float, default=0.8, help="Clustering resolution (default: 0.8)")
    parser.add_argument("--n_neighbors", type=int, default=15, help="Number of neighbors (default: 15)")

    # Visualization parameters
    parser.add_argument("--embedding", type=str, default="umap", choices=["umap", "tsne", "both"],
                        help="Dimensionality reduction visualization method (default: umap)")
    parser.add_argument("--min_dist", type=float, default=0.5, help="UMAP min_dist (default: 0.5)")

    # Marker gene parameters
    parser.add_argument("--marker_method", type=str, default="wilcoxon",
                        choices=["wilcoxon", "t-test", "logreg"],
                        help="Marker gene test method (default: wilcoxon)")
    parser.add_argument("--n_genes", type=int, default=100, help="Number of genes returned per cluster (default: 100)")
    parser.add_argument("--refgenome", type=str, default=None, help="Reference genome directory path (for gene annotation)")

    args = parser.parse_args()

    # Read data
    LOGGER.info(f"Reading data: {args.input}")
    adata = sc.read_h5ad(args.input)

    # 1. PCA dimensionality reduction
    adata, n_pcs = run_pca(
        adata,
        n_comps=args.n_comps,
        n_pcs=args.n_pcs,
        max_value=args.max_value
    )

    # 2. Harmony batch correction (optional)
    if args.use_harmony:
        adata, use_rep = run_harmony(
            adata,
            batch_key=args.batch_key,
            n_pcs=args.n_pcs
        )
    else:
        use_rep = "X_pca"
        LOGGER.info("Not using Harmony, using plain PCA")

    # 3. Clustering
    adata = run_clustering(
        adata,
        n_pcs=args.n_pcs,
        n_neighbors=args.n_neighbors,
        resolution=args.resolution,
        method=args.cluster_method,
        use_rep=use_rep
    )

    # 4. Dimensionality reduction visualization
    adata = run_embedding(
        adata,
        n_pcs=args.n_pcs,
        min_dist=args.min_dist,
        use_rep=use_rep,
        embedding_method=args.embedding
    )

    # 5. Plot visualizations
    color_by = ["clusters"]
    if "sampleid" in adata.obs.columns:
        color_by.append("sampleid")
    if "batchid" in adata.obs.columns and "batchid" != "sampleid":
        color_by.append("batchid")

    plot_embeddings(
        adata,
        os.path.join(args.outdir, "embeddings"),
        color_by=color_by,
        embedding_method=args.embedding
    )

    # 6. Find marker genes
    marker_df = find_marker_genes(
        adata,
        method=args.marker_method,
        n_genes=args.n_genes
    )

    # 7. Save results
    save_results(adata, marker_df, args.outdir, refgenome=args.refgenome)

    LOGGER.info("Analysis complete!")


if __name__ == "__main__":
    main()
