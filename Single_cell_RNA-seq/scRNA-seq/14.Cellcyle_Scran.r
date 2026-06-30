#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ==============================================================================
# 细胞周期分析脚本
# 使用 scran 包进行细胞周期阶段分类（G1/S/G2M）
# ==============================================================================

# ========================== 加载包 ==========================

suppressPackageStartupMessages({
    library(Seurat)
    library(SingleCellExperiment)
    library(scran)
    library(optparse)
    library(pheatmap)
    library(dplyr)
    library(ggplot2)
    library(tibble)
})

rm(list = ls())

# ========================== 配置参数 ==========================

# 路径配置
OUTPUT_DIR <- "14.Cellcyle_Scran.R"
DATA_FILE <- "data_ob_v3.rds"
GENE_NAME_DIR <- "utils/reference/genename2id"

# 分析参数
CELLCYCLE_NAME <- "Scrna_CellCycle"
IDENT_COLUMN <- "clusters"
POINT_SIZE <- 1.5
SPECIES <- "mouse"  # 可选: "human" 或 "mouse"
REDUCTION <- "umap"

# 颜色配置
COLORS <- c("#7fc97f", "#beaed4", "#fdc086")

# ========================== 初始化 ==========================

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ========================== 数据加载 ==========================

seurat_ob <- readRDS(DATA_FILE)
seurat_ob <- SetIdent(seurat_ob, value = IDENT_COLUMN)

# ========================== 加载细胞周期标记基因 ==========================

if (SPECIES %in% c("human", "mouse")) {
    # 从 scran 包读取参考基因对
    ref.pairs <- readRDS(system.file("exdata", paste0(SPECIES, "_cycle_markers.rds"), package = "scran"))

    # 读取基因名到 ID 的映射
    gene_map <- read.table(file.path(GENE_NAME_DIR, paste0(SPECIES, ".txt")), header = FALSE)
    gene_map <- tibble::column_to_rownames(as.data.frame(gene_map), "V1")
    gene_map[, 1] <- as.character(gene_map[, 1])
} else {
    stop("物种不支持，请选择 'human' 或 'mouse'")
}

# 构建参考列表
ref <- list(
    G1 = data.frame(first = gene_map[ref.pairs$G1[, 1], ], second = gene_map[ref.pairs$G1[, 2], ], stringsAsFactors = FALSE),
    S = data.frame(first = gene_map[ref.pairs$S[, 1], ], second = gene_map[ref.pairs$S[, 2], ], stringsAsFactors = FALSE),
    G2M = data.frame(first = gene_map[ref.pairs$G2M[, 1], ], second = gene_map[ref.pairs$G2M[, 2], ], stringsAsFactors = FALSE)
)

# ========================== 细胞周期分类 ==========================

# 检查基因名是否匹配
if (length(intersect(rownames(seurat_ob), ref[[1]]$first)) == 0) {
    stop("基因名不匹配，请检查物种设置")
}

# 转换为 SingleCellExperiment 对象
sce_ob <- as.SingleCellExperiment(seurat_ob)

# 运行 cyclone 分类
assignments <- scran::cyclone(sce_ob, ref)

# 整理结果
scores <- assignments$normalized.scores %>%
    mutate(barcodes = rownames(seurat_ob@meta.data), phase = assignments$phases) %>%
    select(barcodes, everything()) %>%
    tibble::column_to_rownames("barcodes")

colnames(scores) <- paste0("cyclone_", colnames(scores))
colnames(scores)[4] <- CELLCYCLE_NAME

# 添加元数据
seurat_ob <- AddMetaData(seurat_ob, metadata = scores, col.name = colnames(scores))

# 查看结果
head(seurat_ob@meta.data)

# ========================== 可视化 ==========================

# 设置颜色
names(COLORS) <- unique(seurat_ob[[CELLCYCLE_NAME]])

# 绘制 UMAP 图
ggdim <- DimPlot(
    object = seurat_ob,
    dims = c(1, 2),
    reduction = REDUCTION,
    pt.size = POINT_SIZE,
    group.by = CELLCYCLE_NAME
) +
    scale_colour_manual(values = COLORS)

ggsave(
    file.path(OUTPUT_DIR, paste0(REDUCTION, "_visualize_CellCycle_plot.pdf")),
    plot = ggdim, bg = "white"
)

# ========================== 保存结果 ==========================

# 保存元数据
meta.data <- seurat_ob@meta.data %>%
    mutate(cell_barcode = rownames(seurat_ob@meta.data)) %>%
    dplyr::select(cell_barcode, everything())

write.table(
    meta.data,
    file.path(OUTPUT_DIR, "cell_cycle_annotation_result.xls"),
    col.names = TRUE,
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
)
