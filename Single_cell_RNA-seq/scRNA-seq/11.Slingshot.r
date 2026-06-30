#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ==============================================================================
# Slingshot 拟时间分析脚本
# 进行单细胞轨迹推断和基因动态表达分析
# ==============================================================================

# ========================== 加载包 ==========================

library(Seurat)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(Biobase)
library(slingshot)
library(tradeSeq)
library(RColorBrewer)
library(DelayedMatrixStats)
library(scales)
library(paletteer)
library(viridis)

# ========================== 配置参数 ==========================

OUTPUT_DIR <- "11.Slingshot"
DATA_FILE <- "data_ob_v3.rds"
MARKERS_FILE <- "top10_markers_for_each_cluster_anno.xls"
START_CLUS <- "C1"
N_KNOTS <- 5

# ========================== 初始化 ==========================

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# ========================== 数据加载 ==========================

# 读取 Seurat 对象
scRNAsub <- readRDS(DATA_FILE)

# 转换为 SingleCellExperiment 对象
scRNAsub <- as.SingleCellExperiment(scRNAsub, assay = "RNA")

# 设置聚类颜色
colx <- c("#7fc97f", "#beaed4", "#fdc086", "#386cb0")
names(colx) <- unique(scRNAsub$clusters)

# ========================== Slingshot 拟时间推断 ==========================

sce <- slingshot(
    scRNAsub,
    reducedDim = "UMAP",
    start.clus = START_CLUS,
    clusterLabels = scRNAsub$clusters
)

# ========================== 可视化 ==========================

# 绘制聚类与拟时间曲线
pdf(file.path(OUTPUT_DIR, "Clustering_slingshot.pdf"), width = 8)
par(mai = c(1, 1, 1, 2))
plot(reducedDims(sce)[["UMAP"]], col = colx, pch = 16, asp = 1, cex = 1)
xy <- par("usr")
lines(SlingshotDataSet(sce), lwd = 1, col = "black")
legend(
    x = xy[2] + xinch(0.2), y = xy[4], xpd = TRUE,
    legend = names(colx), col = colx, pch = 19
)
dev.off()

# 绘制拟时间梯度
colors <- colorRampPalette(brewer.pal(11, "Spectral")[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks = 100)]

pdf(file.path(OUTPUT_DIR, "slingPseudotime.pdf"), width = 8)
par(mai = c(1, 1, 1, 2))
plot(reducedDims(sce)$UMAP, col = plotcol, pch = 16, asp = 0.5)
lines(SlingshotDataSet(sce), lwd = 2, col = "black")
dev.off()

# ========================== 差异基因与 GAM 拟合 ==========================

# 读取差异基因
diff.genes <- read.delim(MARKERS_FILE)
pseudotimeED <- slingPseudotime(sce, na = FALSE)
cellWeightsED <- slingCurveWeights(sce)

# 准备表达矩阵
counts <- scRNAsub@assays@data$counts
counts <- as.data.frame(counts)
counts <- counts[diff.genes$gene, ]

# 拟合 GAM
sce_slingshot <- fitGAM(
    counts = as.matrix(counts),
    pseudotime = pseudotimeED,
    cellWeights = cellWeightsED,
    nknots = N_KNOTS,
    verbose = TRUE
)
saveRDS(sce_slingshot, file.path(OUTPUT_DIR, "gamlist_sce.rds"))

# 关联检验
ATres <- associationTest(sce_slingshot)
topgenes <- rownames(ATres[order(ATres$pvalue), ])

# ========================== 热图 ==========================

pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
heatdata <- assays(sce)$counts[topgenes, pst.ord]
heatclus <- scRNAsub$clusters[pst.ord]
heatdata <- as.matrix(heatdata)

pdf(file.path(OUTPUT_DIR, "heatmap.pdf"), width = 8)
heatmap(
    log1p(heatdata),
    Colv = NA,
    ColSideColors = brewer.pal(9, "Set1")[heatclus],
    labCol = ""
)
dev.off()

# ========================== 基因动态表达模式 ==========================

# 关联检验结果
mean(rowData(sce_slingshot)$tradeSeq$converged)
rowData(sce_slingshot)$assocRes <- associationTest(sce_slingshot, lineages = TRUE, l2fc = log2(2))
assocRes <- rowData(sce_slingshot)$assocRes

# 绘制基因平滑曲线
gene_dynamic <- list()
genes_plot <- c("Foxp1", "Erbb4")

for (i in seq_along(genes_plot)) {
    p <- plotSmoothers(
        sce_slingshot,
        assays(sce_slingshot)$counts,
        gene = genes_plot[i],
        alpha = 0.6,
        border = TRUE,
        lwd = 2
    ) +
        ggtitle(genes_plot[i])
    gene_dynamic[[i]] <- p
}

ggsave(
    file.path(OUTPUT_DIR, "genes_plot.pdf"),
    gene_dynamic[[1]] + gene_dynamic[[2]],
    height = 2, width = 5
)
