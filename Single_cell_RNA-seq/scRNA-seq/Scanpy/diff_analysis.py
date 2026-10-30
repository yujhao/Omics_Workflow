#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 3: Differential analysis
Input: Clustered data
Function: Run differential analysis for specified groups, output all DEGs + top 20 up- and down-regulated each
"""

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import os
import argparse
import logging
from typing import Optional, List, Union

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
LOGGER = logging.getLogger(__name__)


def load_gene_annotation(annotation_file):
    """
    Load the gene annotation file

    Parameters:
    -----------
    annotation_file : str
        Path to the gene annotation file

    Returns:
    --------
    DataFrame
        Gene annotation data
    """
    LOGGER.info(f"Loading gene annotation: {annotation_file}")

    anno_df = pd.read_csv(annotation_file, sep="\t")

    # Ensure the id column is a string
    anno_df["id"] = anno_df["id"].astype(str)

    LOGGER.info(f"Loaded annotation for {len(anno_df)} genes")

    return anno_df


def add_gene_annotation(df, anno_df, gene_col="gene"):
    """
    Add annotations to the differentially expressed genes

    Parameters:
    -----------
    df : DataFrame
        Differential gene table
    anno_df : DataFrame
        Gene annotation table
    gene_col : str
        Name of the gene column

    Returns:
    --------
    DataFrame
        Differential gene table with annotations added
    """
    # Merge annotations
    df_merged = df.merge(
        anno_df,
        left_on=gene_col,
        right_on="id",
        how="left"
    )

    # Drop the duplicate id column (the id from anno_df)
    if "id_x" in df_merged.columns:
        df_merged = df_merged.drop(columns=["id_x"])
    if "id_y" in df_merged.columns:
        df_merged = df_merged.rename(columns={"id_y": "id"})

    # Fill missing annotations with '--'
    df_merged = df_merged.fillna('--')

    return df_merged


def run_diff_analysis(
    adata,
    groupby="clusters",
    group1=None,
    group2=None,
    method="wilcoxon",
    n_genes=None,
    corr_method="benjamini-hochberg",
    use_raw=True
):
    """
    Run differential expression analysis

    Parameters:
    -----------
    adata : AnnData
        Input AnnData object
    groupby : str
        Grouping key
    group1 : str or list
        First group
    group2 : str or list
        Second group (if None, compare against rest)
    method : str
        Differential analysis method
    n_genes : int
        Number of genes to return (None means all genes)
    corr_method : str
        p-value correction method
    use_raw : bool
        Whether to use raw data

    Returns:
    --------
    DataFrame
        Differential analysis results
    """
    LOGGER.info(f"Starting differential analysis, method={method}...")

    # If raw data exists, use raw
    if use_raw and adata.raw is not None:
        adata_for_test = adata.raw.to_adata()
    else:
        adata_for_test = adata.copy()

    # Normalize
    sc.pp.normalize_total(adata_for_test, target_sum=10000)
    sc.pp.log1p(adata_for_test)

    # Determine the groups to compare
    if group1 is None:
        groups = None
        reference = "rest"
        LOGGER.info(f"Running differential analysis for all groups (vs rest), groupby={groupby}")
    else:
        if isinstance(group1, str):
            group1 = [group1]

        if group2 is None:
            groups = group1
            reference = "rest"
            LOGGER.info(f"Comparing {group1} vs rest, groupby={groupby}")
        else:
            if isinstance(group2, str):
                group2 = [group2]

            # Create a temporary comparison group
            adata_for_test.obs["_comparison_group"] = "rest"
            for g in group1:
                adata_for_test.obs.loc[
                    adata_for_test.obs[groupby] == g,
                    "_comparison_group"
                ] = "_".join(group1)
            for g in group2:
                adata_for_test.obs.loc[
                    adata_for_test.obs[groupby] == g,
                    "_comparison_group"
                ] = "_".join(group2)

            groups = ["_".join(group1)]
            reference = "_".join(group2)
            groupby = "_comparison_group"

            LOGGER.info(f"Comparing {group1} vs {group2}, groupby={groupby}")

    # Run differential analysis
    sc.tl.rank_genes_groups(
        adata_for_test,
        groupby=groupby,
        groups=groups,
        reference=reference,
        method=method,
        n_genes=n_genes,
        corr_method=corr_method,
        use_raw=False,
        pts=True
    )

    # Use scanpy's function to extract results
    all_results = []
    rank_result = adata_for_test.uns["rank_genes_groups"]

    for cluster in rank_result["names"].dtype.names:
        # Use scanpy's function to extract results for each cluster
        de_df = sc.get.rank_genes_groups_df(adata_for_test, group=cluster)

        # If pairwise comparison, add reference group expression percentage
        if "pts" in rank_result and rank_result["pts"].shape[1] == 2:
            pts_df = pd.DataFrame({
                "names": rank_result["pts"].index,
                "pct_nz_reference": rank_result["pts"].iloc[:, 1]
            })
            de_df = de_df.merge(pts_df, on="names", how="left")

        # Add cluster info
        de_df["cluster"] = cluster

        all_results.append(de_df)

    # Combine all results
    df_all = pd.concat(all_results, ignore_index=True)

    # Rename columns
    df_all = df_all.rename(columns={
        "names": "gene",
        "logfoldchanges": "log2FoldChange",
        "pvals": "p-value",
        "pvals_adj": "q-value",
        "pct_nz_group": "pct.1",
        "pct_nz_reference": "pct.2"
    })

    # Compute FoldChange (preserve sign: up-regulated positive, down-regulated negative)
    df_all["FoldChange"] = np.sign(df_all["log2FoldChange"]) * (2 ** df_all["log2FoldChange"].abs())

    # Add id column (same as gene)
    df_all["id"] = df_all["gene"]

    LOGGER.info(f"Differential analysis complete, {len(df_all)} genes in total")

    return df_all


def filter_and_sort_results(
    df,
    qval_cutoff=0.05,
    fc_cutoff=1.5
):
    """
    Filter and sort the results

    Parameters:
    -----------
    df : DataFrame
        Differential analysis results
    qval_cutoff : float
        q-value threshold
    fc_cutoff : float
        FoldChange threshold

    Returns:
    --------
    tuple
        (significant DEG table, top 20 up-regulated, top 20 down-regulated, statistics)
    """
    # Filter significantly differentially expressed genes (compare absolute log2FoldChange with log2(fc_cutoff))
    df_sig = df[
        (df["q-value"] < qval_cutoff) &
        (df["log2FoldChange"].abs() >= np.log2(fc_cutoff))
    ].copy()

    LOGGER.info(f"Number of significant DEGs (q<{qval_cutoff}, FC>={fc_cutoff}): {len(df_sig)}")

    # Add a Regulation column for significant DEGs
    df_sig["Regulation"] = np.where(
        df_sig["log2FoldChange"] > 0,
        "Up",
        np.where(
            df_sig["log2FoldChange"] < 0,
            "Down",
            "NoChange"
        )
    )

    # Count up- and down-regulated genes
    up_count = len(df_sig[df_sig["Regulation"] == "Up"])
    down_count = len(df_sig[df_sig["Regulation"] == "Down"])

    # Top 20 up-regulated genes
    df_up = df_sig[df_sig["Regulation"] == "Up"].copy()
    df_up = df_up.sort_values("log2FoldChange", ascending=False).head(20)

    # Top 20 down-regulated genes
    df_down = df_sig[df_sig["Regulation"] == "Down"].copy()
    df_down = df_down.sort_values("log2FoldChange", ascending=True).head(20)

    LOGGER.info(f"Top 20 up-regulated genes: {len(df_up)}")
    LOGGER.info(f"Top 20 down-regulated genes: {len(df_down)}")

    # Return statistics
    stats = {
        "Up_diff": up_count,
        "Down_diff": down_count,
        "Total_diff": up_count + down_count
    }

    return df_sig, df_up, df_down, stats


def save_results(df_all, df_sig, df_up, df_down, stats, outdir, group1, group2, qval_cutoff, fc_cutoff):
    """
    Save the results

    Parameters:
    -----------
    df_all : DataFrame
        All differentially expressed genes
    df_sig : DataFrame
        Significant differentially expressed genes
    df_up : DataFrame
        Top 20 up-regulated genes
    df_down : DataFrame
        Top 20 down-regulated genes
    stats : dict
        Statistics
    outdir : str
        Output directory
    group1 : list
        First group
    group2 : list
        Second group
    qval_cutoff : float
        q-value threshold
    fc_cutoff : float
        FoldChange threshold
    """
    os.makedirs(outdir, exist_ok=True)

    # Build filename suffix
    group1_name = "_".join(group1) if group1 else "all"
    group2_name = "_".join(group2) if group2 else "rest"
    suffix = f"group_{group1_name}-vs-{group2_name}"

    # Reorder columns
    def reorder_columns(df, include_regulation=False):
        """Reorder the columns"""
        # Columns to drop
        drop_cols = ["cluster", "scores"]

        # Base columns
        base_cols = ["gene", "log2FoldChange", "p-value", "q-value", "pct.1", "pct.2", "FoldChange"]
        if include_regulation and "Regulation" in df.columns:
            base_cols.append("Regulation")

        # Annotation columns
        anno_cols = ["id", "ensembl_id", "gene_type", "gene_description", "TFs_Family", "GO_id", "GO_term", "KEGG_id", "KEGG_description"]

        # Keep only existing columns
        ordered_cols = [col for col in base_cols + anno_cols if col in df.columns and col not in drop_cols]

        # Add any other existing columns (excluding those to be dropped)
        other_cols = [col for col in df.columns if col not in ordered_cols and col not in drop_cols]
        ordered_cols.extend(other_cols)

        return df[ordered_cols]

    # Save all DEGs
    df_all_ordered = reorder_columns(df_all, include_regulation=False)
    df_all_ordered.to_csv(
        os.path.join(outdir, f"{suffix}-all_diffexp_genes_anno.xls"),
        sep="\t",
        index=False
    )
    LOGGER.info(f"Saved all differentially expressed genes: {suffix}-all_diffexp_genes_anno.xls")

    # Save significant DEGs
    df_sig_ordered = reorder_columns(df_sig, include_regulation=True)
    df_sig_ordered.to_csv(
        os.path.join(outdir, f"{suffix}-diff-qval-{qval_cutoff}-FC-{fc_cutoff}_anno.xls"),
        sep="\t",
        index=False
    )
    LOGGER.info(f"Saved significant differentially expressed genes: {suffix}-diff-qval-{qval_cutoff}-FC-{fc_cutoff}_anno.xls")

    # Save top 20 up-regulated genes
    df_up_ordered = reorder_columns(df_up, include_regulation=True)
    df_up_ordered.to_csv(
        os.path.join(outdir, f"{suffix}-up_top20_anno.xls"),
        sep="\t",
        index=False
    )
    LOGGER.info(f"Saved top 20 up-regulated genes: {suffix}-up_top20_anno.xls")

    # Save top 20 down-regulated genes
    df_down_ordered = reorder_columns(df_down, include_regulation=True)
    df_down_ordered.to_csv(
        os.path.join(outdir, f"{suffix}-down_top20_anno.xls"),
        sep="\t",
        index=False
    )
    LOGGER.info(f"Saved top 20 down-regulated genes: {suffix}-down_top20_anno.xls")

    # Save statistics
    stats_df = pd.DataFrame([{
        "Case": group1_name,
        "Control": group2_name,
        "Up_diff": stats["Up_diff"],
        "Down_diff": stats["Down_diff"],
        f"Total_diff(q-value<{qval_cutoff}&FoldChange>{fc_cutoff})": stats["Total_diff"]
    }])
    stats_df.to_csv(
        os.path.join(outdir, "diffexp_results_stat.xls"),
        sep="\t",
        index=False
    )
    LOGGER.info(f"Saved statistics: diffexp_results_stat.xls")


def main():
    """
    Main function
    """
    parser = argparse.ArgumentParser(
        description="Step 3: Differential analysis (specify groups, output all DEGs + top 20 up- and down-regulated each)"
    )

    # Input/Output
    parser.add_argument("-i", "--input", required=True, help="Input file path (clustered .h5ad)")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory")

    # Gene annotation
    parser.add_argument(
        "--annotation",
        type=str,
        default="/hwstorage/oe-scrna/jhyu/Git_lab/sh/yujunhao/refdata-gex-GRCm39-2024-A/gene_annotation.xls",
        help="Path to the gene annotation file"
    )

    # Grouping parameters
    parser.add_argument(
        "--groupby",
        type=str,
        default="clusters",
        help="Grouping key (default: clusters)"
    )
    parser.add_argument(
        "--group1",
        type=str,
        nargs="+",
        default=None,
        help="First group (may be multiple clusters, e.g. 0 1 2)"
    )
    parser.add_argument(
        "--group2",
        type=str,
        nargs="+",
        default=None,
        help="Second group (may be multiple clusters; if None, compare against rest)"
    )

    # Differential analysis parameters
    parser.add_argument(
        "--method",
        type=str,
        default="wilcoxon",
        choices=["wilcoxon", "t-test", "logreg"],
        help="Differential analysis method (default: wilcoxon)"
    )
    parser.add_argument(
        "--n_genes",
        type=int,
        default=None,
        help="Number of genes to return (None means all genes)"
    )
    parser.add_argument(
        "--qval_cutoff",
        type=float,
        default=0.05,
        help="q-value threshold (default: 0.05)"
    )
    parser.add_argument(
        "--fc_cutoff",
        type=float,
        default=1.5,
        help="FoldChange threshold (default: 1.5)"
    )

    args = parser.parse_args()

    # Read data
    LOGGER.info(f"Reading data: {args.input}")
    adata = sc.read_h5ad(args.input)

    # Load gene annotation
    anno_df = load_gene_annotation(args.annotation)

    # Run differential analysis
    df_all = run_diff_analysis(
        adata,
        groupby=args.groupby,
        group1=args.group1,
        group2=args.group2,
        method=args.method,
        n_genes=args.n_genes
    )

    # Add gene annotation
    df_all = add_gene_annotation(df_all, anno_df)

    # Filter and sort
    df_sig, df_up, df_down, stats = filter_and_sort_results(
        df_all,
        qval_cutoff=args.qval_cutoff,
        fc_cutoff=args.fc_cutoff
    )

    # Annotate significant DEGs too (already done)
    # Annotate top 20 too (already done)

    # Save results
    group1_list = args.group1 if args.group1 else ["all"]
    group2_list = args.group2 if args.group2 else ["rest"]

    save_results(
        df_all,
        df_sig,
        df_up,
        df_down,
        stats,
        args.outdir,
        group1_list,
        group2_list,
        args.qval_cutoff,
        args.fc_cutoff
    )

    LOGGER.info("Differential analysis complete!")


if __name__ == "__main__":
    main()
