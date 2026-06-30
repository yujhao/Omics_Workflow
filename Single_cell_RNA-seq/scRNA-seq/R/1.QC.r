#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ==============================================================================
# 单细胞数据质控脚本
# 包括数据过滤、线粒体/血红蛋白比例计算、双细胞去除等
# ==============================================================================

# ========================== 加载包 ==========================

suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratWrappers)
    library(DoubletFinder)
    library(batchelor)
    library(GSEABase)
    library(org.Hs.eg.db)
    library(clusterProfiler)
    library(SingleR)
    library(ggplot2)
    library(patchwork)
    library(pheatmap)
    library(tidyverse)
    library(dplyr)
    library(ggstatsplot)
    library(tidyr)
    library(Matrix)
    library(infercnv)
    library(tibble)
    library(RColorBrewer)
    library(ComplexHeatmap)
    library(dittoSeq)
    library(plyranges)
    library(monocle)
    library(ggsci)
    library(igraph)
    library(CellChat)
    library(ggalluvial)
    library(SingleCellExperiment)
    library(magrittr)
})

source("plotting/0.Outliers.r")

# ========================== 配置参数 ==========================

OUTPUT_DIR <- "1.SingleCell_QC"
DATA_FILE <- "data_ob_v3.rds"
MT_GMT_FILE <- "utils/reference/refdata-gex-GRCm39-2024-A/MT_genelist.gmt"

# 过滤参数
FILTER_PARAMS <- c("nFeature_RNA", "nCount_RNA", "percent.mito")
LOWER_THRESHOLD <- c(NA, NA, -Inf)
UPPER_THRESHOLD <- c(NA, NA, 1)

# 双细胞去除参数
DOUBLET_RATE <- 0.008
PC_NUM <- 1:30
RESOLUTION <- 0.3

# 血红蛋白基因列表
HB_GENES <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")

# ========================== 初始化 ==========================

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# ========================== 数据加载 ==========================

scRNA1 <- readRDS(DATA_FILE)

# ========================== 1.1 数据整理 ==========================

# 设置默认 ID
Idents(scRNA1) <- "orig.ident"

# 计算线粒体基因比例
scRNA1[["percent.mito"]] <- PercentageFeatureSet(scRNA1, pattern = "^MT-")

# 计算血红蛋白基因比例
HB_m <- match(HB_GENES, rownames(scRNA1@assays$RNA))
HB_matched <- rownames(scRNA1@assays$RNA)[HB_m]
HB_matched <- HB_matched[!is.na(HB_matched)]
scRNA1[["percent.HB"]] <- PercentageFeatureSet(scRNA1, features = HB_matched)

# 绘制 QC 前小提琴图
beforeQC_vlnplot <- VlnPlot(
    scRNA1,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.HB"),
    ncol = 4,
    pt.size = 0
)
ggsave(file.path(OUTPUT_DIR, "BeforeQC_nFeature_nCount_percent.mito_percent.HB_vlnplot.pdf"),
       plot = beforeQC_vlnplot)
ggsave(file.path(OUTPUT_DIR, "BeforeQC_nFeature_nCount_percent.mito_percent.HB_vlnplot.png"),
       plot = beforeQC_vlnplot)

# ========================== 1.2 数据过滤 ==========================

# 设置过滤阈值
names(LOWER_THRESHOLD) <- FILTER_PARAMS
names(UPPER_THRESHOLD) <- FILTER_PARAMS

# 构建阈值列表
bounds_list <- list()
for (x in FILTER_PARAMS) {
    bounds_list[[x]] <- c(min = LOWER_THRESHOLD[x], max = UPPER_THRESHOLD[x])
}

# 检测异常值
outliers <- FindOutliers(
    scRNA1,
    vars = FILTER_PARAMS,
    var.limit = bounds_list,
    batch = "sampleid",
    type = "both",
    cut.1 = "mean",
    cut.2 = "sd",
    n = 2,
    log = FALSE
)

# 标记异常细胞
outliercells <- do.call(cbind, outliers)
metric_outlier <- apply(outliercells, 1, function(x) any(x == TRUE))
scRNA1 <- AddMetaData(scRNA1, metadata = metric_outlier, col.name = "is_metric_outlier")

# 过滤异常细胞
outlier_variables <- "is_metric_outlier"
is_valid_cell <- !apply(FetchData(scRNA1, vars = outlier_variables), 1, function(x) any(x == TRUE))
scRNA1 <- AddMetaData(scRNA1, metadata = is_valid_cell, col.name = "is_valid")
scRNA1 <- subset(scRNA1, subset = is_valid == TRUE)

# ========================== 线粒体基因过滤 ==========================

# 读取线粒体基因列表
counts <- GetAssayData(scRNA1, "counts")
gmt_list <- unlist(strsplit(MT_GMT_FILE, ",", perl = TRUE))
gset_list <- lapply(gmt_list, function(gmtfile) {
    GSEABase::geneIds(GSEABase::getGmt(con = gmtfile))
})

# 根据线粒体 UMI 过滤细胞
char_vector <- unlist(gset_list[[1]])
mt <- counts[char_vector, ]
mt_umi <- colSums(mt)
median_value <- summary(mt_umi)["Median"]
fust <- as.data.frame(mt_umi)
fust$bak <- fust$mt_umi
filter_cell <- rownames(fust[(fust$mt_umi < median_value * 4), ])
scRNA1 <- scRNA1[, filter_cell]

# 计算 log10 Genes per UMI
scRNA1@meta.data$log10GenesPerUMI <- log10(scRNA1@meta.data$nFeature_RNA) / log10(scRNA1@meta.data$nCount_RNA)

# ========================== 双细胞去除 ==========================

obj <- SplitObject(scRNA1, split.by = "orig.ident")
obj_rm <- list()
doublets_plot <- list()

for (i in names(obj)) {
    print(i)
    obj[[i]] <- NormalizeData(obj[[i]])
    obj[[i]] <- FindVariableFeatures(obj[[i]], selection.method = "vst", nfeatures = 2000)
    obj[[i]] <- ScaleData(obj[[i]])
    obj[[i]] <- RunPCA(obj[[i]])
    obj[[i]] <- RunUMAP(obj[[i]], dims = 1:30)
    obj[[i]] <- FindNeighbors(obj[[i]], dims = 1:30) %>%
        FindClusters(resolution = RESOLUTION)

    # 去除双细胞
    tmp <- RemoveDoublets(obj[[i]], doublet.rate = DOUBLET_RATE, pc.num = PC_NUM)
    obj_rm[[i]] <- tmp$obj
    doublets_plot[[i]] <- tmp$plot
}

# 合并处理后的对象
scRNA1 <- obj_rm[[1]]
if (length(obj_rm) > 1) {
    for (i in 2:length(obj_rm)) {
        scRNA1 <- merge(scRNA1, y = obj_rm[[i]])
    }
}

Idents(scRNA1) <- "orig.ident"

# 绘制 QC 后小提琴图
afterQC_vlnplot <- VlnPlot(
    scRNA1,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.HB"),
    ncol = 4,
    pt.size = 0
)
ggsave(file.path(OUTPUT_DIR, "afterQC_nFeature_nCount_percent.mito_percent.HB_vlnplot.pdf"),
       plot = afterQC_vlnplot, width = 8, height = 6)
ggsave(file.path(OUTPUT_DIR, "afterQC_nFeature_nCount_percent.mito_percent.HB_vlnplot.png"),
       plot = afterQC_vlnplot, width = 8, height = 6)

# ========================== 1.3 数据归一化与标准化 ==========================

scRNA1 <- NormalizeData(scRNA1)
scRNA1 <- FindVariableFeatures(scRNA1, selection.method = "vst")
scRNA1 <- ScaleData(scRNA1, features = VariableFeatures(scRNA1))
