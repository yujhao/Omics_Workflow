#!/bin/bash
# 单细胞分析运行脚本


# ============================================
# 第一步：QC 分析
# ============================================

# 方式 1: 从原始 10x 数据读取多个样本
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
    # --skip_doublet_removal  # 添加此参数可跳过双细胞去除


# # ============================================
# # 第二步：降维聚类分析
# # ============================================
INPUT_CLUSTER="vi/1.QC/qc_result.h5ad"  # QC 后的数据
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
# # 第三步：差异分析
# # ============================================
INPUT_DIFF="vi/2.Clustering/clustering_result.h5ad"  # 聚类后的数据
OUTPUT_DIFF="vi/3.diff_output"

# 基本用法：比较 sample A1 vs A2
python diff_analysis.py \
    -i $INPUT_DIFF \
    -o $OUTPUT_DIFF \
    --groupby sampleid \
    --group1 A1 \
    --group2 A2 \
    --qval_cutoff 0.05 \
    --fc_cutoff 1.5

# 注意，python仅仅分析，绘图还是使用R脚本绘图,h5ad转成rds，参考h5ad_rds.R

