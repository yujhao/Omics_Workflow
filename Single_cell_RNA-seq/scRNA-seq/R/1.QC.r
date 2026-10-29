#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ==============================================================================
# Single-cell data quality control script
# Includes data filtering, mitochondrial/hemoglobin ratio calculation, doublet removal, etc.
# ==============================================================================

# ========================== Load packages ==========================

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

# ========================== Configuration parameters ==========================

OUTPUT_DIR <- "1.SingleCell_QC"
DATA_FILE <- "data_ob_v3.rds"
MT_GMT_FILE <- "utils/reference/refdata-gex-GRCm39-2024-A/MT_genelist.gmt"

# Filtering parameters
FILTER_PARAMS <- c("nFeature_RNA", "nCount_RNA", "percent.mito")
LOWER_THRESHOLD <- c(NA, NA, -Inf)
UPPER_THRESHOLD <- c(NA, NA, 1)

# Doublet removal parameters
DOUBLET_RATE <- 0.008
PC_NUM <- 1:30
RESOLUTION <- 0.3

# Hemoglobin gene list
HB_GENES <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")

# ========================== Initialization ==========================

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# ========================== Data loading ==========================

scRNA1 <- readRDS(DATA_FILE)

# ========================== 1.1 Data processing ==========================

# Set default identities
Idents(scRNA1) <- "orig.ident"

# Calculate mitochondrial gene proportion
scRNA1[["percent.mito"]] <- PercentageFeatureSet(scRNA1, pattern = "^MT-")

# Calculate hemoglobin gene proportion
HB_m <- match(HB_GENES, rownames(scRNA1@assays$RNA))
HB_matched <- rownames(scRNA1@assays$RNA)[HB_m]
HB_matched <- HB_matched[!is.na(HB_matched)]
scRNA1[["percent.HB"]] <- PercentageFeatureSet(scRNA1, features = HB_matched)

# Plot violin plots before QC
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

# ========================== 1.2 Data filtering ==========================

# Set filtering thresholds
names(LOWER_THRESHOLD) <- FILTER_PARAMS
names(UPPER_THRESHOLD) <- FILTER_PARAMS

# Build threshold list
bounds_list <- list()
for (x in FILTER_PARAMS) {
    bounds_list[[x]] <- c(min = LOWER_THRESHOLD[x], max = UPPER_THRESHOLD[x])
}

# Detect outliers
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

# Flag outlier cells
outliercells <- do.call(cbind, outliers)
metric_outlier <- apply(outliercells, 1, function(x) any(x == TRUE))
scRNA1 <- AddMetaData(scRNA1, metadata = metric_outlier, col.name = "is_metric_outlier")

# Filter outlier cells
outlier_variables <- "is_metric_outlier"
is_valid_cell <- !apply(FetchData(scRNA1, vars = outlier_variables), 1, function(x) any(x == TRUE))
scRNA1 <- AddMetaData(scRNA1, metadata = is_valid_cell, col.name = "is_valid")
scRNA1 <- subset(scRNA1, subset = is_valid == TRUE)

# ========================== Mitochondrial gene filtering ==========================

# Read mitochondrial gene list
counts <- GetAssayData(scRNA1, "counts")
gmt_list <- unlist(strsplit(MT_GMT_FILE, ",", perl = TRUE))
gset_list <- lapply(gmt_list, function(gmtfile) {
    GSEABase::geneIds(GSEABase::getGmt(con = gmtfile))
})

# Filter cells by mitochondrial UMI
char_vector <- unlist(gset_list[[1]])
mt <- counts[char_vector, ]
mt_umi <- colSums(mt)
median_value <- summary(mt_umi)["Median"]
fust <- as.data.frame(mt_umi)
fust$bak <- fust$mt_umi
filter_cell <- rownames(fust[(fust$mt_umi < median_value * 4), ])
scRNA1 <- scRNA1[, filter_cell]

# Calculate log10 Genes per UMI
scRNA1@meta.data$log10GenesPerUMI <- log10(scRNA1@meta.data$nFeature_RNA) / log10(scRNA1@meta.data$nCount_RNA)

# ========================== Doublet removal ==========================

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

    # Remove doublets
    tmp <- RemoveDoublets(obj[[i]], doublet.rate = DOUBLET_RATE, pc.num = PC_NUM)
    obj_rm[[i]] <- tmp$obj
    doublets_plot[[i]] <- tmp$plot
}

# Merge processed objects
scRNA1 <- obj_rm[[1]]
if (length(obj_rm) > 1) {
    for (i in 2:length(obj_rm)) {
        scRNA1 <- merge(scRNA1, y = obj_rm[[i]])
    }
}

Idents(scRNA1) <- "orig.ident"

# Plot violin plots after QC
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

# ========================== 1.3 Data normalization and standardization ==========================

scRNA1 <- NormalizeData(scRNA1)
scRNA1 <- FindVariableFeatures(scRNA1, selection.method = "vst")
scRNA1 <- ScaleData(scRNA1, features = VariableFeatures(scRNA1))
