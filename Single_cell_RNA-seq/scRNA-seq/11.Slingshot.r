#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ==============================================================================
# Slingshot pseudotime analysis script
# Performs single-cell trajectory inference and gene dynamic expression analysis
# ==============================================================================

# ========================== Load packages ==========================

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

# ========================== Configuration parameters ==========================

OUTPUT_DIR <- "11.Slingshot"
DATA_FILE <- "data_ob_v3.rds"
MARKERS_FILE <- "top10_markers_for_each_cluster_anno.xls"
START_CLUS <- "C1"
N_KNOTS <- 5

# ========================== Initialization ==========================

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# ========================== Data loading ==========================

# Read Seurat object
scRNAsub <- readRDS(DATA_FILE)

# Convert to SingleCellExperiment object
scRNAsub <- as.SingleCellExperiment(scRNAsub, assay = "RNA")

# Set cluster colors
colx <- c("#7fc97f", "#beaed4", "#fdc086", "#386cb0")
names(colx) <- unique(scRNAsub$clusters)

# ========================== Slingshot pseudotime inference ==========================

sce <- slingshot(
    scRNAsub,
    reducedDim = "UMAP",
    start.clus = START_CLUS,
    clusterLabels = scRNAsub$clusters
)

# ========================== Visualization ==========================

# Plot clusters and pseudotime curves
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

# Plot pseudotime gradient
colors <- colorRampPalette(brewer.pal(11, "Spectral")[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks = 100)]

pdf(file.path(OUTPUT_DIR, "slingPseudotime.pdf"), width = 8)
par(mai = c(1, 1, 1, 2))
plot(reducedDims(sce)$UMAP, col = plotcol, pch = 16, asp = 0.5)
lines(SlingshotDataSet(sce), lwd = 2, col = "black")
dev.off()

# ========================== Differential genes and GAM fitting ==========================

# Read differential genes
diff.genes <- read.delim(MARKERS_FILE)
pseudotimeED <- slingPseudotime(sce, na = FALSE)
cellWeightsED <- slingCurveWeights(sce)

# Prepare expression matrix
counts <- scRNAsub@assays@data$counts
counts <- as.data.frame(counts)
counts <- counts[diff.genes$gene, ]

# Fit GAM
sce_slingshot <- fitGAM(
    counts = as.matrix(counts),
    pseudotime = pseudotimeED,
    cellWeights = cellWeightsED,
    nknots = N_KNOTS,
    verbose = TRUE
)
saveRDS(sce_slingshot, file.path(OUTPUT_DIR, "gamlist_sce.rds"))

# Association test
ATres <- associationTest(sce_slingshot)
topgenes <- rownames(ATres[order(ATres$pvalue), ])

# ========================== Heatmap ==========================

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

# ========================== Gene dynamic expression patterns ==========================

# Association test results
mean(rowData(sce_slingshot)$tradeSeq$converged)
rowData(sce_slingshot)$assocRes <- associationTest(sce_slingshot, lineages = TRUE, l2fc = log2(2))
assocRes <- rowData(sce_slingshot)$assocRes

# Plot gene smoothing curves
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
