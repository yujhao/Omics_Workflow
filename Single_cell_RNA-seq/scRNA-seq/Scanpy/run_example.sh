#!/bin/bash
# Single-cell analysis run script


# ============================================
# Step 1: QC Analysis
# ============================================

# Method 1: Read multiple samples from raw 10x data
ROOT_DIR="data/"
SAMPLES="A1 A2 P1 P2"
OUTPUT_QC="vi/1.QC"
REFGENOME="refdata-gex-GRCm39-2024-A/"

python qc_analysis.py \
    -i $ROOT_DIR \
    --samples $SAMPLES \
    -o $OUTPUT_QC \
    --refgenome $REFGENOME \
    --mt_genelist MT_genelist.gmt \
    --hb_genelist HB_genelist.gmt \
    --mingene 200 \
    --maxgene 7500 \
    --minumi 1000 \
    --maxumi 50000 \
    --mtfilter 0.2 \
    --hbfilter 0.05
    # --skip_doublet_removal  # Add this argument to skip doublet removal


# # ============================================
# # Step 2: Dimensionality reduction and clustering analysis
# # ============================================
INPUT_CLUSTER="vi/1.QC/qc_result.h5ad"  # Data after QC
OUTPUT_CLUSTER="vi/2.Clustering"

基本用法（默认 UMAP，PCA，resolution=0.8）
REFGENOME="refdata-gex-GRCm39-2024-A/"

python clustering_analysis.py \
    -i $INPUT_CLUSTER \
    -o $OUTPUT_CLUSTER \
    --cluster_method leiden \
    --resolution 0.8 \
    --embedding umap \
    --use_harmony \
    --batch_key batchid \
    --refgenome $REFGENOME


# # ============================================
# # Step 3: Differential analysis
# # ============================================
INPUT_DIFF="vi/2.Clustering/clustering_result.h5ad"  # Data after clustering
OUTPUT_DIFF="vi/3.diff_output"

# Basic usage: compare sample A1 vs A2
python diff_analysis.py \
    -i $INPUT_DIFF \
    -o $OUTPUT_DIFF \
    --groupby sampleid \
    --group1 A1 \
    --group2 A2 \
    --qval_cutoff 0.05 \
    --fc_cutoff 1.5

# Note: Python only performs the analysis; plotting is still done with R scripts. Convert h5ad to rds; see h5ad_rds.R

