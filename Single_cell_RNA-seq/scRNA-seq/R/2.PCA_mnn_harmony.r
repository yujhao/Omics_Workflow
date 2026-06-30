#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ==============================================================================
# 批次矫正和降维聚类脚本
# 包括 PCA、MNN、Harmony 三种降维方法
# ==============================================================================

library(Seurat)
library(SeuratWrappers)
library(batchelor)
library(ggplot2)

# ========================== 配置参数 ==========================

OUTPUT_DIR <- "2.Clustering"
DATA_FILE <- "data_ob_v3.rds"
RESOLUTION <- 0.4
DIMS <- 1:30
BATCH_VAR <- "sampleid"

# ========================== 初始化 ==========================

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# ========================== 数据加载 ==========================

scRNA1 <- readRDS(DATA_FILE)

# ========================== PCA 降维 ==========================

scRNA_pca <- RunPCA(scRNA1)
scRNA_pca <- FindNeighbors(scRNA_pca, reduction = "pca", dims = DIMS)
scRNA_pca <- FindClusters(scRNA_pca, resolution = RESOLUTION)
scRNA_pca <- RunUMAP(scRNA_pca, reduction = "pca", dims = DIMS)

# 设置聚类标签
scRNA_pca@meta.data$clusters <- as.numeric(scRNA_pca@meta.data$seurat_clusters)
scRNA_pca@meta.data$clusters <- as.factor(scRNA_pca@meta.data$clusters)
Idents(scRNA_pca) <- "clusters"

# 可视化
p1 <- DimPlot(scRNA_pca, pt.size = 0.1, label = TRUE)
ggsave(file.path(OUTPUT_DIR, "pca.pdf"), p1, width = 8, height = 8)
ggsave(file.path(OUTPUT_DIR, "pca.png"), p1, width = 8, height = 8)

# 保存结果
saveRDS(scRNA_pca, file.path(OUTPUT_DIR, "pca.rds"))

# ========================== MNN 降维 ==========================

# 按样本拆分
scRNAlist <- SplitObject(scRNA1, split.by = "orig.ident")
scRNAlist <- lapply(scRNAlist, FUN = function(x) NormalizeData(x))
scRNAlist <- lapply(scRNAlist, FUN = function(x) FindVariableFeatures(x))

# 运行 MNN
scRNA_mnn <- RunFastMNN(object.list = scRNAlist)
scRNA_mnn <- FindVariableFeatures(scRNA_mnn)
scRNA_mnn <- RunUMAP(scRNA_mnn, reduction = "mnn", dims = DIMS)
scRNA_mnn <- FindNeighbors(scRNA_mnn, reduction = "mnn", dims = DIMS)
scRNA_mnn <- FindClusters(scRNA_mnn, resolution = RESOLUTION)

# 设置聚类标签
scRNA_mnn@meta.data$clusters <- as.numeric(scRNA_mnn@meta.data$seurat_clusters)
scRNA_mnn@meta.data$clusters <- as.factor(scRNA_mnn@meta.data$clusters)
Idents(scRNA_mnn) <- "clusters"

# 可视化
p2 <- DimPlot(scRNA_mnn, pt.size = 0.1, label = TRUE)
ggsave(file.path(OUTPUT_DIR, "MNN.pdf"), p2, width = 8, height = 8)
ggsave(file.path(OUTPUT_DIR, "MNN.png"), p2, width = 8, height = 8)

# 保存结果
saveRDS(scRNA_mnn, file.path(OUTPUT_DIR, "MNN.rds"))

# ========================== Harmony 降维 ==========================

scRNA_Harmony <- RunPCA(object = scRNA1)
scRNA_Harmony <- FindNeighbors(scRNA_Harmony, reduction = "pca", dims = DIMS)
scRNA_Harmony <- FindClusters(scRNA_Harmony, resolution = RESOLUTION)
scRNA_Harmony <- RunHarmony(scRNA_Harmony, BATCH_VAR)
scRNA_Harmony <- RunUMAP(scRNA_Harmony, dims = DIMS, reduction = "harmony")

# 设置聚类标签
scRNA_Harmony@meta.data$clusters <- as.numeric(scRNA_Harmony@meta.data$seurat_clusters)
scRNA_Harmony@meta.data$clusters <- as.factor(scRNA_Harmony@meta.data$clusters)
Idents(scRNA_Harmony) <- "clusters"

# 可视化
p3 <- DimPlot(scRNA_Harmony, pt.size = 0.1, label = TRUE)
ggsave(file.path(OUTPUT_DIR, "Harmony.pdf"), p3, width = 8, height = 8)
ggsave(file.path(OUTPUT_DIR, "Harmony.png"), p3, width = 8, height = 8)

# 保存结果
saveRDS(scRNA_Harmony, file.path(OUTPUT_DIR, "Harmony.rds"))
