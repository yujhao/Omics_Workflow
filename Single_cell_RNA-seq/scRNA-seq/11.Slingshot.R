# ==============================================================================
# Slingshot 拟时间分析
# ==============================================================================

# 加载包
library(Seurat)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(Biobase)
library(dplyr)
library(slingshot)
library(tradeSeq)
library(RColorBrewer)
library(DelayedMatrixStats)
library(scales)
library(paletteer)
library(viridis)
library(grDevices)

# 读取数据
scRNAsub <- readRDS("data_ob_v3.rds")
output_dir <- "slingshot_analysis"

# 转换为 SingleCellExperiment 对象
scRNAsub <- as.SingleCellExperiment(scRNAsub, assay = "RNA")

# 设置聚类颜色
colx <- c("#7fc97f", "#beaed4", "#fdc086", "#386cb0")
names(colx) <- unique(scRNAsub$clusters)

# ------------------------------------------------------------------------------
# 1. Slingshot 拟时间推断
# ------------------------------------------------------------------------------
sce <- slingshot(scRNAsub,
                 reducedDim = "UMAP",
                 start.clus = "C1",
                 clusterLabels = scRNAsub$clusters)

# 绘制聚类与拟时间曲线
pdf(file.path(output_dir, "Clustering_slingshot.pdf"), width = 8)
par(mai = c(1, 1, 1, 2))
plot(reducedDims(sce)[["UMAP"]], col = colx, pch = 16, asp = 1, cex = 1)
xy <- par("usr")
lines(SlingshotDataSet(sce), lwd = 1, col = "black")
legend(x = xy[2] + xinch(0.2), y = xy[4], xpd = TRUE,
       legend = names(colx), col = colx, pch = 19)
dev.off()

# 绘制拟时间梯度
colors <- colorRampPalette(brewer.pal(11, "Spectral")[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks = 100)]

pdf(file.path(output_dir, "slingPseudotime.pdf"), width = 8)
par(mai = c(1, 1, 1, 2))
plot(reducedDims(sce)$UMAP, col = plotcol, pch = 16, asp = 0.5)
lines(SlingshotDataSet(sce), lwd = 2, col = "black")
dev.off()

# ------------------------------------------------------------------------------
# 2. 差异基因与 GAM 拟合
# ------------------------------------------------------------------------------
diff.genes <- read.delim("top10_markers_for_each_cluster_anno.xls")
pseudotimeED <- slingPseudotime(sce, na = FALSE)
cellWeightsED <- slingCurveWeights(sce)

counts <- scRNAsub@assays@data$counts
counts <- as.data.frame(counts)
counts <- counts[diff.genes$gene, ]

sce_slingshot <- fitGAM(counts = as.matrix(counts),
                        pseudotime = pseudotimeED,
                        cellWeights = cellWeightsED,
                        nknots = 5, verbose = TRUE)
saveRDS(sce_slingshot, file.path(output_dir, "gamlist_sce.rds"))

ATres <- associationTest(sce_slingshot)
topgenes <- rownames(ATres[order(ATres$pvalue), ])

# ------------------------------------------------------------------------------
# 3. 热图
# ------------------------------------------------------------------------------
pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
heatdata <- assays(sce)$counts[topgenes, pst.ord]
heatclus <- scRNAsub$clusters[pst.ord]
heatdata <- as.matrix(heatdata)

pdf(file.path(output_dir, "heatmap.pdf"), width = 8)
heatmap(log1p(heatdata),
        Colv = NA,
        ColSideColors = brewer.pal(9, "Set1")[heatclus],
        labCol = "")
dev.off()

# ------------------------------------------------------------------------------
# 4. 基因动态表达模式
# ------------------------------------------------------------------------------
mean(rowData(sce_slingshot)$tradeSeq$converged)
rowData(sce_slingshot)$assocRes <- associationTest(sce_slingshot, lineages = TRUE, l2fc = log2(2))
assocRes <- rowData(sce_slingshot)$assocRes

gene_dynamic <- list()
genes_plot <- c("Foxp1", "Erbb4")

for (i in seq_along(genes_plot)) {
  p <- plotSmoothers(sce_slingshot, assays(sce_slingshot)$counts,
                     gene = genes_plot[i], alpha = 0.6, border = TRUE, lwd = 2) +
    ggtitle(genes_plot[i])
  gene_dynamic[[i]] <- p
}

ggsave(file.path(output_dir, "genes_plot.pdf"),
       gene_dynamic[[1]] + gene_dynamic[[2]],
       height = 2, width = 5)
