#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ==============================================================================
# hdWGCNA (high-dimensional weighted gene co-expression network analysis)
# Builds metacells, constructs a co-expression network, detects gene modules and
# computes module eigengenes / hub genes. Follows the official hdWGCNA tutorial.
# ==============================================================================

# ========================== Load packages ==========================

suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratDisk)
    library(WGCNA)
    library(hdWGCNA)
    library(tidyverse)
    library(cowplot)
    library(patchwork)
    library(UCell)
    library(qlcMatrix)
})

# ========================== Configuration parameters ==========================

INPUT_FILE <- "data/data_ob.h5seurat"
INFORMAT <- "h5seurat"          # h5seurat or rds
OUTPUT_DIR <- "3.hdWGCNA"

GROUP_BY <- "clusters"          # metadata column used to group metacells
REDUCT1 <- "pca"                # reduction used to build metacells (KNN)
REDUCT2 <- "umap"               # reduction used for downstream plots
ASSAY <- "RNA"
TOP_N_HUB <- 10                 # top hub genes by kME

# ========================== Initialization ==========================

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
theme_set(theme_cowplot())
set.seed(12345)
enableWGCNAThreads(nThreads = 8)

# ========================== Data loading ==========================

if (INFORMAT == "h5seurat") {
    seurat_obj <- SeuratDisk::LoadH5Seurat(INPUT_FILE, assays = ASSAY, verbose = FALSE)
} else {
    seurat_obj <- readRDS(INPUT_FILE)
}

# ========================== Metacells ==========================

seurat_obj <- SetupForWGCNA(
    seurat_obj,
    gene_select = "fraction",
    fraction = 0.05,
    wgcna_name = "hdWGCNA"
)

# k should stay in the 25-75 range, scaled with the cell number
n_cells <- dim(seurat_obj)[2]
k <- if (n_cells <= 50000) 25 else if (n_cells <= 80000) 55 else 75

seurat_obj <- MetacellsByGroups(
    seurat_obj = seurat_obj,
    group.by = c(GROUP_BY),
    reduction = REDUCT1,
    k = k,
    max_shared = 30,
    ident.group = GROUP_BY
)
seurat_obj <- NormalizeMetacells(seurat_obj)

# ========================== Co-expression network ==========================

group_name <- unique(seurat_obj@misc$hdWGCNA$wgcna_metacell_obj[[GROUP_BY]][[GROUP_BY]])

seurat_obj <- SetDatExpr(
    seurat_obj,
    group_name = group_name,
    group.by = GROUP_BY,
    assay = ASSAY,
    slot = "data"
)

seurat_obj <- TestSoftPowers(seurat_obj, networkType = "signed")
pdf(file.path(OUTPUT_DIR, "PlotSoftPowers.pdf"))
print(wrap_plots(PlotSoftPowers(seurat_obj), ncol = 2))
dev.off()

seurat_obj <- ConstructNetwork(
    seurat_obj,
    soft_power = NULL,          # selected automatically from TestSoftPowers
    setDatExpr = FALSE,
    tom_outdir = OUTPUT_DIR,
    wgcna_name = "hdWGCNA",
    networkType = "signed",
    minModuleSize = 20,
    deepSplit = 4
)

pdf(file.path(OUTPUT_DIR, "PlotDendrogram.pdf"))
PlotDendrogram(seurat_obj, main = "hdWGCNA Dendrogram")
dev.off()

# ========================== Module eigengenes and connectivity ==========================

# ScaleData must be run before ModuleEigengenes
seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
seurat_obj <- ModuleEigengenes(seurat_obj)

seurat_obj <- ModuleConnectivity(seurat_obj, group.by = GROUP_BY, group_name = group_name)

# ========================== Hub genes ==========================

modules <- GetModules(seurat_obj, wgcna_name = "hdWGCNA")
write.table(modules, file = file.path(OUTPUT_DIR, "gene_modules_allocation.xls"),
            sep = "\t", quote = FALSE, row.names = FALSE)

hub_genes <- GetHubGenes(seurat_obj, n_hubs = TOP_N_HUB, wgcna_name = "hdWGCNA")
write.table(hub_genes, file = file.path(OUTPUT_DIR, "modules_topn_hubgenes.xls"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Per-cell module scores based on the top hub genes
seurat_obj <- ModuleExprScore(seurat_obj, n_genes = TOP_N_HUB, method = "UCell")

# ========================== Plots ==========================

# hMEs on the embedding
plot_list <- ModuleFeaturePlot(seurat_obj, wgcna_name = "hdWGCNA",
                               features = "hMEs", order = TRUE, reduction = REDUCT2)
pdf(file.path(OUTPUT_DIR, "1.ModuleFeaturePlot_hMEs.pdf"))
print(wrap_plots(plot_list, ncol = 4))
dev.off()

# Hub gene module scores on the embedding
plot_list <- ModuleFeaturePlot(seurat_obj, wgcna_name = "hdWGCNA",
                               features = "scores", order = "shuffle",
                               ucell = TRUE, reduction = REDUCT2)
pdf(file.path(OUTPUT_DIR, "2.ModuleFeaturePlot_scores.pdf"))
print(wrap_plots(plot_list, ncol = 4))
dev.off()

# Module-module correlation
pdf(file.path(OUTPUT_DIR, "3.ModuleCorrelogram.pdf"))
ModuleCorrelogram(seurat_obj, sig.level = 0.05)
dev.off()

# hME dot plot across groups
MEs <- GetMEs(seurat_obj, harmonized = TRUE)
mods <- colnames(MEs)
mods <- mods[mods != "grey"]
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

p <- DotPlot(seurat_obj, features = mods, group.by = GROUP_BY) +
    coord_flip() + RotatedAxis() +
    scale_color_gradient2(high = "#005E00", mid = "white", low = "#004A97")
pdf(file.path(OUTPUT_DIR, "4.ModuleDotPlot.pdf"))
print(p)
dev.off()

# ========================== Save ==========================

saveRDS(seurat_obj, file = file.path(OUTPUT_DIR, "hdWGCNA_object.rds"))

message("hdWGCNA finished, results written to ", OUTPUT_DIR)
