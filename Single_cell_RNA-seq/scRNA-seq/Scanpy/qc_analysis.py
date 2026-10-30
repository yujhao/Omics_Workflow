#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 1: Quality control (QC)
Input: raw data (10x output directory or multiple samples)
Function: QC filtering + doublet removal
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


def load_gene_list(path):
    """
    Load gene list from a GMT file
    """
    genes = []
    with open(path, "r") as f:
        for line in f:
            genes.extend(line.strip().split())

    genes = [g for g in genes if g and g != "NA"]
    return np.unique(genes)


def load_samples(root_dir, sample_dirs):
    """
    Read 10x data for multiple samples and merge them

    Parameters:
    -----------
    root_dir : str
        Root directory of the sample data
    sample_dirs : list
        List of sample directories

    Returns:
    --------
    AnnData
        Merged AnnData object
    """
    LOGGER.info(f"Reading {len(sample_dirs)} samples from {root_dir}...")

    adata_list = []
    for sample in sample_dirs:
        sample_path = os.path.join(root_dir, sample)

        # Try to read the h5 file
        h5_path = os.path.join(sample_path, "filtered_feature_bc_matrix.h5")
        mtx_path = os.path.join(sample_path, "filtered_feature_bc_matrix")

        if os.path.exists(h5_path):
            LOGGER.info(f"Reading sample {sample}: {h5_path}")
            adata = sc.read_10x_h5(h5_path)
        elif os.path.exists(mtx_path):
            LOGGER.info(f"Reading sample {sample}: {mtx_path}")
            adata = sc.read_10x_mtx(mtx_path, var_names='gene_symbols', cache=False)
        else:
            raise FileNotFoundError(f"Data files not found for sample {sample}")

        # Add sample info
        adata.obs['sampleid'] = sample
        adata.obs['batchid'] = sample  # Batch info is the same as the sample

        # Ensure var_names are unique
        adata.var_names_make_unique()

        adata_list.append(adata)

    # Merge all samples
    LOGGER.info(f"Merging {len(adata_list)} samples...")
    adata_merged = ad.concat(adata_list, join='outer', merge='same')

    # Ensure obs_names are unique
    adata_merged.obs_names_make_unique()

    LOGGER.info(f"Merging complete: {adata_merged.n_obs} cells, {adata_merged.n_vars} genes")

    return adata_merged


def calculate_qc_metrics(
    adata,
    refgenome,
    mt_genelist_file="MT_genelist.gmt",
    hb_genelist_file="HB_genelist.gmt"
):
    """
    Calculate QC metrics - implementation based on qc.py

    Parameters:
    -----------
    adata : AnnData
        Input AnnData object
    refgenome : str
        Path to the reference genome directory
    mt_genelist_file : str
        Mitochondrial gene list file name
    hb_genelist_file : str
        Hemoglobin gene list file name
    """
    LOGGER.info("Computing QC metrics...")

    # Load gene lists
    mt_genes = load_gene_list(os.path.join(refgenome, mt_genelist_file))
    hb_genes = load_gene_list(os.path.join(refgenome, hb_genelist_file))

    LOGGER.info(f"Number of mitochondrial genes: {len(mt_genes)}")
    LOGGER.info(f"Number of hemoglobin genes: {len(hb_genes)}")

    # Flag genes (use isin instead of startswith)
    adata.var["mt"] = adata.var_names.isin(mt_genes)
    adata.var["hb"] = adata.var_names.isin(hb_genes)

    # Compute QC metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "hb"],
        inplace=True
    )

    # Standardize naming (Seurat style)
    adata.obs["nFeature_RNA"] = adata.obs["n_genes_by_counts"]
    adata.obs["nCount_RNA"] = adata.obs["total_counts"]
    adata.obs["percent.mito"] = adata.obs["pct_counts_mt"] / 100
    adata.obs["percent.HB"] = adata.obs["pct_counts_hb"] / 100

    LOGGER.info(f"QC metrics computed for {adata.n_obs} cells, {adata.n_vars} genes")

    return adata


def filter_cells(
    adata,
    mingene=200,
    maxgene=7500,
    minumi=1000,
    maxumi=50000,
    mtfilter=0.2,
    hbfilter=0.05
):
    """
    Filter cells - correct implementation based on qc.py

    Parameters:
    -----------
    adata : AnnData
        Input AnnData object
    mingene : int
        Minimum number of genes per cell
    maxgene : int
        Maximum number of genes per cell
    minumi : int
        Minimum UMI count per cell
    maxumi : int
        Maximum UMI count per cell
    mtfilter : float
        Upper limit of mitochondrial gene fraction
    hbfilter : float
        Upper limit of hemoglobin gene fraction
    """
    LOGGER.info("Filtering cells...")

    n_before = adata.n_obs

    # Filter conditions
    filter_mask = (
        (adata.obs["nFeature_RNA"] >= mingene) &
        (adata.obs["nFeature_RNA"] <= maxgene) &
        (adata.obs["nCount_RNA"] >= minumi) &
        (adata.obs["nCount_RNA"] <= maxumi) &
        (adata.obs["percent.mito"] <= mtfilter) &
        (adata.obs["percent.HB"] <= hbfilter)
    )

    # Apply filtering
    adata = adata[filter_mask].copy()

    n_after = adata.n_obs
    LOGGER.info(f"Filtering complete: {n_before} -> {n_after} cells ({n_before - n_after} removed)")

    return adata


def remove_doublets(adata, random_state=2025):
    """
    Remove doublets

    Parameters:
    -----------
    adata : AnnData
        Input AnnData object
    random_state : int
        Random seed
    """
    LOGGER.info("Starting doublet detection...")

    try:
        import doubletdetection

        n_before = adata.n_obs

        # Save original counts
        adata.layers["counts"] = adata.X.copy()

        # Run doublet detection
        clf = doubletdetection.BoostClassifier(
            n_iters=10,
            random_state=random_state
        )
        doublets = clf.fit(adata.X).predict()

        # Remove doublets
        adata = adata[doublets == 0].copy()

        n_after = adata.n_obs
        LOGGER.info(f"Doublet removal complete: {n_before} -> {n_after} cells ({n_before - n_after} removed)")

    except ImportError:
        LOGGER.warning("doubletdetection is not installed, skipping doublet removal")

    return adata


def plot_qc_violin(adata, output_dir, title_suffix, group_col="sampleid"):
    """
    Plot QC violin plots

    Parameters:
    -----------
    adata : AnnData
        Input AnnData object
    output_dir : str
        Output directory
    title_suffix : str
        Title suffix (e.g. "beforeQC" or "afterQC")
    group_col : str
        Column name for grouping
    """
    LOGGER.info(f"Plot QC violin plots ({title_suffix})...")

    os.makedirs(output_dir, exist_ok=True)

    metrics = ["nFeature_RNA", "nCount_RNA", "percent.mito", "percent.HB"]

    # Check whether the grouping column exists
    if group_col not in adata.obs.columns:
        LOGGER.warning(f"Column {group_col} does not exist, plotting all data")
        adata.obs["_group"] = "All"
        group_col = "_group"

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    for idx, metric in enumerate(metrics):
        if metric not in adata.obs.columns:
            continue

        sns.violinplot(
            data=adata.obs,
            x=group_col,
            y=metric,
            ax=axes[idx],
            inner="quartile",
            palette="Set2"
        )

        axes[idx].set_title(f"{metric} - {title_suffix}", fontsize=12, fontweight='bold')
        axes[idx].set_xlabel(group_col, fontsize=10)
        axes[idx].set_ylabel(metric, fontsize=10)
        axes[idx].tick_params(axis='x', rotation=45)

    plt.tight_layout()

    # Save figure
    plt.savefig(os.path.join(output_dir, f"QC_metrics_{title_suffix}.png"), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, f"QC_metrics_{title_suffix}.pdf"), bbox_inches='tight')
    plt.close()

    LOGGER.info(f"Saved: {output_dir}/QC_metrics_{title_suffix}.png")
    LOGGER.info(f"Saved: {output_dir}/QC_metrics_{title_suffix}.pdf")


def generate_qc_statistics(adata_before, adata_after_filter, adata_after_doublet, outdir, group_col="sampleid"):
    """
    Generate QC statistics files

    Parameters:
    -----------
    adata_before : AnnData
        Data before QC
    adata_after_filter : AnnData
        Data after filtering
    adata_after_doublet : AnnData
        Data after doublet removal
    outdir : str
        Output directory
    group_col : str
        Column name for grouping
    """
    LOGGER.info("Generating QC statistics files...")

    os.makedirs(outdir, exist_ok=True)

    # Get all sample IDs
    if group_col in adata_before.obs.columns:
        sample_ids = adata_before.obs[group_col].unique()
    else:
        sample_ids = ["All"]

    # 1. Generate cell_count_before_after_QC.xls
    cell_count_data = []
    for sample in sample_ids:
        if group_col in adata_before.obs.columns:
            before_data = adata_before.obs[adata_before.obs[group_col] == sample]
            filter_data = adata_after_filter.obs[adata_after_filter.obs[group_col] == sample]
            doublet_data = adata_after_doublet.obs[adata_after_doublet.obs[group_col] == sample]
        else:
            before_data = adata_before.obs
            filter_data = adata_after_filter.obs
            doublet_data = adata_after_doublet.obs

        row = {
            "sampleid": sample,
            "cell_count_beforeQC": len(before_data),
            "cell_count_afterQC": len(filter_data),
            "cell_count_after_rmdoublets": len(doublet_data),
            "nFeature_RNA_max_after_rmdoublets": doublet_data["nFeature_RNA"].max(),
            "nFeature_RNA_min_after_rmdoublets": doublet_data["nFeature_RNA"].min(),
            "nFeature_RNA_mean_after_rmdoublets": doublet_data["nFeature_RNA"].mean(),
            "nFeature_RNA_median_after_rmdoublets": doublet_data["nFeature_RNA"].median(),
            "nCount_RNA_max_after_rmdoublets": doublet_data["nCount_RNA"].max(),
            "nCount_RNA_min_after_rmdoublets": doublet_data["nCount_RNA"].min(),
            "nCount_RNA_mean_after_rmdoublets": doublet_data["nCount_RNA"].mean(),
            "nCount_RNA_median_after_rmdoublets": doublet_data["nCount_RNA"].median(),
            "percent.mito_max_after_rmdoublets": doublet_data["percent.mito"].max(),
            "percent.mito_min_after_rmdoublets": doublet_data["percent.mito"].min(),
            "percent.mito_mean_after_rmdoublets": doublet_data["percent.mito"].mean(),
            "percent.mito_median_after_rmdoublets": doublet_data["percent.mito"].median(),
            "percent.HB_max_after_rmdoublets": doublet_data["percent.HB"].max(),
            "percent.HB_min_after_rmdoublets": doublet_data["percent.HB"].min(),
            "percent.HB_mean_after_rmdoublets": doublet_data["percent.HB"].mean(),
            "percent.HB_median_after_rmdoublets": doublet_data["percent.HB"].median(),
        }
        cell_count_data.append(row)

    cell_count_df = pd.DataFrame(cell_count_data)
    cell_count_df.to_csv(os.path.join(outdir, "cell_count_before_after_QC.xls"), sep="\t", index=False)
    LOGGER.info("Saved cell_count_before_after_QC.xls")

    # 2. Generate cell_statitics_before_after_QC.xls
    cell_stat_data = []
    for sample in sample_ids:
        if group_col in adata_before.obs.columns:
            before_data = adata_before.obs[adata_before.obs[group_col] == sample]
            filter_data = adata_after_filter.obs[adata_after_filter.obs[group_col] == sample]
            doublet_data = adata_after_doublet.obs[adata_after_doublet.obs[group_col] == sample]
        else:
            before_data = adata_before.obs
            filter_data = adata_after_filter.obs
            doublet_data = adata_after_doublet.obs

        row = {
            "sampleid": sample,
            "cell_count_beforeQC": len(before_data),
            "nFeature_RNA_max_beforeQC": before_data["nFeature_RNA"].max() if "nFeature_RNA" in before_data.columns else 0,
            "nFeature_RNA_min_beforeQC": before_data["nFeature_RNA"].min() if "nFeature_RNA" in before_data.columns else 0,
            "nFeature_RNA_mean_beforeQC": before_data["nFeature_RNA"].mean() if "nFeature_RNA" in before_data.columns else 0,
            "nFeature_RNA_median_beforeQC": before_data["nFeature_RNA"].median() if "nFeature_RNA" in before_data.columns else 0,
            "nCount_RNA_max_beforeQC": before_data["nCount_RNA"].max() if "nCount_RNA" in before_data.columns else 0,
            "nCount_RNA_min_beforeQC": before_data["nCount_RNA"].min() if "nCount_RNA" in before_data.columns else 0,
            "nCount_RNA_mean_beforeQC": before_data["nCount_RNA"].mean() if "nCount_RNA" in before_data.columns else 0,
            "nCount_RNA_median_beforeQC": before_data["nCount_RNA"].median() if "nCount_RNA" in before_data.columns else 0,
            "percent.mito_max_beforeQC": before_data["percent.mito"].max() if "percent.mito" in before_data.columns else 0,
            "percent.mito_min_beforeQC": before_data["percent.mito"].min() if "percent.mito" in before_data.columns else 0,
            "percent.mito_mean_beforeQC": before_data["percent.mito"].mean() if "percent.mito" in before_data.columns else 0,
            "percent.mito_median_beforeQC": before_data["percent.mito"].median() if "percent.mito" in before_data.columns else 0,
            "percent.HB_max_beforeQC": before_data["percent.HB"].max() if "percent.HB" in before_data.columns else 0,
            "percent.HB_min_beforeQC": before_data["percent.HB"].min() if "percent.HB" in before_data.columns else 0,
            "percent.HB_mean_beforeQC": before_data["percent.HB"].mean() if "percent.HB" in before_data.columns else 0,
            "percent.HB_median_beforeQC": before_data["percent.HB"].median() if "percent.HB" in before_data.columns else 0,
            "cell_count_afterQC": len(filter_data),
            "nFeature_RNA_max_afterQC": filter_data["nFeature_RNA"].max() if "nFeature_RNA" in filter_data.columns else 0,
            "nFeature_RNA_min_afterQC": filter_data["nFeature_RNA"].min() if "nFeature_RNA" in filter_data.columns else 0,
            "nFeature_RNA_mean_afterQC": filter_data["nFeature_RNA"].mean() if "nFeature_RNA" in filter_data.columns else 0,
            "nFeature_RNA_median_afterQC": filter_data["nFeature_RNA"].median() if "nFeature_RNA" in filter_data.columns else 0,
            "nCount_RNA_max_afterQC": filter_data["nCount_RNA"].max() if "nCount_RNA" in filter_data.columns else 0,
            "nCount_RNA_min_afterQC": filter_data["nCount_RNA"].min() if "nCount_RNA" in filter_data.columns else 0,
            "nCount_RNA_mean_afterQC": filter_data["nCount_RNA"].mean() if "nCount_RNA" in filter_data.columns else 0,
            "nCount_RNA_median_afterQC": filter_data["nCount_RNA"].median() if "nCount_RNA" in filter_data.columns else 0,
            "percent.mito_max_afterQC": filter_data["percent.mito"].max() if "percent.mito" in filter_data.columns else 0,
            "percent.mito_min_afterQC": filter_data["percent.mito"].min() if "percent.mito" in filter_data.columns else 0,
            "percent.mito_mean_afterQC": filter_data["percent.mito"].mean() if "percent.mito" in filter_data.columns else 0,
            "percent.mito_median_afterQC": filter_data["percent.mito"].median() if "percent.mito" in filter_data.columns else 0,
            "percent.HB_max_afterQC": filter_data["percent.HB"].max() if "percent.HB" in filter_data.columns else 0,
            "percent.HB_min_afterQC": filter_data["percent.HB"].min() if "percent.HB" in filter_data.columns else 0,
            "percent.HB_mean_afterQC": filter_data["percent.HB"].mean() if "percent.HB" in filter_data.columns else 0,
            "percent.HB_median_afterQC": filter_data["percent.HB"].median() if "percent.HB" in filter_data.columns else 0,
            "cell_count_after_rmdoublets": len(doublet_data),
            "nFeature_RNA_max_after_rmdoublets": doublet_data["nFeature_RNA"].max() if "nFeature_RNA" in doublet_data.columns else 0,
            "nFeature_RNA_min_after_rmdoublets": doublet_data["nFeature_RNA"].min() if "nFeature_RNA" in doublet_data.columns else 0,
            "nFeature_RNA_mean_after_rmdoublets": doublet_data["nFeature_RNA"].mean() if "nFeature_RNA" in doublet_data.columns else 0,
            "nFeature_RNA_median_after_rmdoublets": doublet_data["nFeature_RNA"].median() if "nFeature_RNA" in doublet_data.columns else 0,
            "nCount_RNA_max_after_rmdoublets": doublet_data["nCount_RNA"].max() if "nCount_RNA" in doublet_data.columns else 0,
            "nCount_RNA_min_after_rmdoublets": doublet_data["nCount_RNA"].min() if "nCount_RNA" in doublet_data.columns else 0,
            "nCount_RNA_mean_after_rmdoublets": doublet_data["nCount_RNA"].mean() if "nCount_RNA" in doublet_data.columns else 0,
            "nCount_RNA_median_after_rmdoublets": doublet_data["nCount_RNA"].median() if "nCount_RNA" in doublet_data.columns else 0,
            "percent.mito_max_after_rmdoublets": doublet_data["percent.mito"].max() if "percent.mito" in doublet_data.columns else 0,
            "percent.mito_min_after_rmdoublets": doublet_data["percent.mito"].min() if "percent.mito" in doublet_data.columns else 0,
            "percent.mito_mean_after_rmdoublets": doublet_data["percent.mito"].mean() if "percent.mito" in doublet_data.columns else 0,
            "percent.mito_median_after_rmdoublets": doublet_data["percent.mito"].median() if "percent.mito" in doublet_data.columns else 0,
            "percent.HB_max_after_rmdoublets": doublet_data["percent.HB"].max() if "percent.HB" in doublet_data.columns else 0,
            "percent.HB_min_after_rmdoublets": doublet_data["percent.HB"].min() if "percent.HB" in doublet_data.columns else 0,
            "percent.HB_mean_after_rmdoublets": doublet_data["percent.HB"].mean() if "percent.HB" in doublet_data.columns else 0,
            "percent.HB_median_after_rmdoublets": doublet_data["percent.HB"].median() if "percent.HB" in doublet_data.columns else 0,
        }
        cell_stat_data.append(row)

    cell_stat_df = pd.DataFrame(cell_stat_data)
    cell_stat_df.to_csv(os.path.join(outdir, "cell_statitics_before_after_QC.xls"), sep="\t", index=False)
    LOGGER.info("Saved cell_statitics_before_after_QC.xls")


def save_results(adata, outdir, adata_before=None, adata_after_filter=None, group_col="sampleid"):
    """
    Save results
    """
    LOGGER.info(f"Saving results to {outdir}...")

    # Create output directory
    os.makedirs(outdir, exist_ok=True)

    # Save h5ad
    adata.write_h5ad(os.path.join(outdir, "qc_result.h5ad"))
    LOGGER.info("Saved qc_result.h5ad")

    # Save QC metrics
    adata.obs.to_csv(os.path.join(outdir, "qc_metrics.csv"), index=False)
    LOGGER.info("Saved qc_metrics.csv")

    # Generate QC statistics files
    if adata_before is not None and adata_after_filter is not None:
        generate_qc_statistics(adata_before, adata_after_filter, adata, outdir, group_col=group_col)

        # Plot QC violin plots (before filtering)
        plot_qc_violin(adata_before, outdir, "beforeQC", group_col=group_col)

        # Plot QC violin plots (after filtering)
        plot_qc_violin(adata_after_filter, outdir, "afterQC", group_col=group_col)

    LOGGER.info("All results saved")


def main():
    """
    Main function
    """
    parser = argparse.ArgumentParser(
        description="Step 1: Quality control (QC filtering + doublet removal)"
    )

    # Input/Output
    parser.add_argument("-i", "--input", required=True, help="Input file path (.h5ad) or root directory containing sample directories")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory")
    parser.add_argument("--samples", nargs="+", default=None, help="List of sample directories (e.g. A1 A2 P1 P2)")

    # Reference genome
    parser.add_argument("--refgenome", required=True, help="Path to the reference genome directory")
    parser.add_argument("--mt_genelist", default="MT_genelist.gmt", help="Mitochondrial gene list file name")
    parser.add_argument("--hb_genelist", default="HB_genelist.gmt", help="Hemoglobin gene list file name")

    # QC filtering parameters
    parser.add_argument("--mingene", type=int, default=200, help="Minimum number of genes per cell (default: 200)")
    parser.add_argument("--maxgene", type=int, default=7500, help="Maximum number of genes per cell (default: 7500)")
    parser.add_argument("--minumi", type=int, default=1000, help="Minimum UMI count per cell (default: 1000)")
    parser.add_argument("--maxumi", type=int, default=50000, help="Maximum UMI count per cell (default: 50000)")
    parser.add_argument("--mtfilter", type=float, default=0.2, help="Upper limit of mitochondrial gene fraction (default: 0.2)")
    parser.add_argument("--hbfilter", type=float, default=0.05, help="Upper limit of hemoglobin gene fraction (default: 0.05)")

    # Doublet detection
    parser.add_argument("--doublet_random_state", type=int, default=2025, help="Random seed for doublet detection (default: 2025)")
    parser.add_argument("--skip_doublet_removal", action="store_true", help="Skip the doublet removal step")

    args = parser.parse_args()

    # Read data
    if args.samples:
        # Read data from multiple sample directories
        LOGGER.info(f"Reading multiple samples from directory: {args.input}")
        adata = load_samples(args.input, args.samples)
    else:
        # Read a single h5ad file
        LOGGER.info(f"Reading data: {args.input}")
        adata = sc.read_h5ad(args.input)

    # 1. Compute QC metrics
    adata = calculate_qc_metrics(
        adata,
        refgenome=args.refgenome,
        mt_genelist_file=args.mt_genelist,
        hb_genelist_file=args.hb_genelist
    )

    # Save pre-filtering data (for statistics)
    adata_before = adata.copy()

    # 2. Filter cells
    adata_after_filter = filter_cells(
        adata,
        mingene=args.mingene,
        maxgene=args.maxgene,
        minumi=args.minumi,
        maxumi=args.maxumi,
        mtfilter=args.mtfilter,
        hbfilter=args.hbfilter
    )
    adata = adata_after_filter

    # 3. Remove doublets (optional)
    if args.skip_doublet_removal:
        LOGGER.info("Skip the doublet removal step")
    else:
        adata = remove_doublets(adata, random_state=args.doublet_random_state)

    # 4. Save results (including QC statistics)
    save_results(
        adata,
        args.outdir,
        adata_before=adata_before,
        adata_after_filter=adata_after_filter,
        group_col="sampleid" if "sampleid" in adata.obs.columns else None
    )

    LOGGER.info("QC analysis complete!")


if __name__ == "__main__":
    main()
