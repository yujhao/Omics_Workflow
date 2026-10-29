#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ==============================================================================
# Batch correction and dimensionality reduction/clustering script
# Includes three dimensionality reduction methods: PCA, MNN, and Harmony
# ==============================================================================

library(Seurat)
library(SeuratWrappers)
library(batchelor)
library(ggplot2)

# ========================== Configuration parameters ==========================

OUTPUT_DIR <- "2.Clustering"
DATA_FILE <- "data_ob_v3.rds"
RESOLUTION <- 0.4
DIMS <- 1:30
BATCH_VAR <- "sampleid"

# ========================== Initialization ==========================

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# ========================== Data loading ==========================

scRNA1 <- readRDS(DATA_FILE)

# ========================== PCA dimensionality reduction ==========================

scRNA_pca <- RunPCA(scRNA1)
scRNA_pca <- FindNeighbors(scRNA_pca, reduction = "pca", dims = DIMS)
scRNA_pca <- FindClusters(scRNA_pca, resolution = RESOLUTION)
scRNA_pca <- RunUMAP(scRNA_pca, reduction = "pca", dims = DIMS)

# Set cluster labels
scRNA_pca@meta.data$clusters <- as.numeric(scRNA_pca@meta.data$seurat_clusters)
scRNA_pca@meta.data$clusters <- as.factor(scRNA_pca@meta.data$clusters)
Idents(scRNA_pca) <- "clusters"

# Visualization
p1 <- DimPlot(scRNA_pca, pt.size = 0.1, label = TRUE)
ggsave(file.path(OUTPUT_DIR, "pca.pdf"), p1, width = 8, height = 8)
ggsave(file.path(OUTPUT_DIR, "pca.png"), p1, width = 8, height = 8)

# Save results
saveRDS(scRNA_pca, file.path(OUTPUT_DIR, "pca.rds"))

# ========================== MNN dimensionality reduction ==========================

# Split by sample
scRNAlist <- SplitObject(scRNA1, split.by = "orig.ident")
scRNAlist <- lapply(scRNAlist, FUN = function(x) NormalizeData(x))
scRNAlist <- lapply(scRNAlist, FUN = function(x) FindVariableFeatures(x))

# Run MNN
scRNA_mnn <- RunFastMNN(object.list = scRNAlist)
scRNA_mnn <- FindVariableFeatures(scRNA_mnn)
scRNA_mnn <- RunUMAP(scRNA_mnn, reduction = "mnn", dims = DIMS)
scRNA_mnn <- FindNeighbors(scRNA_mnn, reduction = "mnn", dims = DIMS)
scRNA_mnn <- FindClusters(scRNA_mnn, resolution = RESOLUTION)

# Set cluster labels
scRNA_mnn@meta.data$clusters <- as.numeric(scRNA_mnn@meta.data$seurat_clusters)
scRNA_mnn@meta.data$clusters <- as.factor(scRNA_mnn@meta.data$clusters)
Idents(scRNA_mnn) <- "clusters"

# Visualization
p2 <- DimPlot(scRNA_mnn, pt.size = 0.1, label = TRUE)
ggsave(file.path(OUTPUT_DIR, "MNN.pdf"), p2, width = 8, height = 8)
ggsave(file.path(OUTPUT_DIR, "MNN.png"), p2, width = 8, height = 8)

# Save results
saveRDS(scRNA_mnn, file.path(OUTPUT_DIR, "MNN.rds"))

# ========================== Harmony dimensionality reduction ==========================

scRNA_Harmony <- RunPCA(object = scRNA1)
scRNA_Harmony <- FindNeighbors(scRNA_Harmony, reduction = "pca", dims = DIMS)
scRNA_Harmony <- FindClusters(scRNA_Harmony, resolution = RESOLUTION)
scRNA_Harmony <- RunHarmony(scRNA_Harmony, BATCH_VAR)
scRNA_Harmony <- RunUMAP(scRNA_Harmony, dims = DIMS, reduction = "harmony")

# Set cluster labels
scRNA_Harmony@meta.data$clusters <- as.numeric(scRNA_Harmony@meta.data$seurat_clusters)
scRNA_Harmony@meta.data$clusters <- as.factor(scRNA_Harmony@meta.data$clusters)
Idents(scRNA_Harmony) <- "clusters"

# Visualization
p3 <- DimPlot(scRNA_Harmony, pt.size = 0.1, label = TRUE)
ggsave(file.path(OUTPUT_DIR, "Harmony.pdf"), p3, width = 8, height = 8)
ggsave(file.path(OUTPUT_DIR, "Harmony.png"), p3, width = 8, height = 8)

# Save results
saveRDS(scRNA_Harmony, file.path(OUTPUT_DIR, "Harmony.rds"))
