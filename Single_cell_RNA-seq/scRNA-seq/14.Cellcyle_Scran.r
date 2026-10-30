#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ==============================================================================
# Cell cycle analysis script
# Classifies cell cycle stages (G1/S/G2M) using the scran package
# ==============================================================================

# ========================== Load packages ==========================

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

# ========================== Configuration parameters ==========================

# Path configuration
OUTPUT_DIR <- "14.Cellcyle_Scran.R"
DATA_FILE <- "data_ob_v3.rds"
GENE_NAME_DIR <- "utils/reference/genename2id"

# Analysis parameters
CELLCYCLE_NAME <- "Scrna_CellCycle"
IDENT_COLUMN <- "clusters"
POINT_SIZE <- 1.5
SPECIES <- "mouse"  # Options: "human" or "mouse"
REDUCTION <- "umap"

# Color configuration
COLORS <- c("#7fc97f", "#beaed4", "#fdc086")

# ========================== Initialization ==========================

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ========================== Data loading ==========================

seurat_ob <- readRDS(DATA_FILE)
seurat_ob <- SetIdent(seurat_ob, value = IDENT_COLUMN)

# ========================== Load cell cycle marker genes ==========================

if (SPECIES %in% c("human", "mouse")) {
    # Read reference gene pairs from the scran package
    ref.pairs <- readRDS(system.file("exdata", paste0(SPECIES, "_cycle_markers.rds"), package = "scran"))

    # Read gene name to ID mapping
    gene_map <- read.table(file.path(GENE_NAME_DIR, paste0(SPECIES, ".txt")), header = FALSE)
    gene_map <- tibble::column_to_rownames(as.data.frame(gene_map), "V1")
    gene_map[, 1] <- as.character(gene_map[, 1])
} else {
    stop("Species not supported, please choose 'human' or 'mouse'")
}

# Build reference list
ref <- list(
    G1 = data.frame(first = gene_map[ref.pairs$G1[, 1], ], second = gene_map[ref.pairs$G1[, 2], ], stringsAsFactors = FALSE),
    S = data.frame(first = gene_map[ref.pairs$S[, 1], ], second = gene_map[ref.pairs$S[, 2], ], stringsAsFactors = FALSE),
    G2M = data.frame(first = gene_map[ref.pairs$G2M[, 1], ], second = gene_map[ref.pairs$G2M[, 2], ], stringsAsFactors = FALSE)
)

# ========================== Cell cycle classification ==========================

# Check whether gene names match
if (length(intersect(rownames(seurat_ob), ref[[1]]$first)) == 0) {
    stop("Gene names do not match, please check the species setting")
}

# Convert to SingleCellExperiment object
sce_ob <- as.SingleCellExperiment(seurat_ob)

# Run cyclone classification
assignments <- scran::cyclone(sce_ob, ref)

# Organize results
scores <- assignments$normalized.scores %>%
    mutate(barcodes = rownames(seurat_ob@meta.data), phase = assignments$phases) %>%
    select(barcodes, everything()) %>%
    tibble::column_to_rownames("barcodes")

colnames(scores) <- paste0("cyclone_", colnames(scores))
colnames(scores)[4] <- CELLCYCLE_NAME

# Add metadata
seurat_ob <- AddMetaData(seurat_ob, metadata = scores, col.name = colnames(scores))

# Inspect results
head(seurat_ob@meta.data)

# ========================== Visualization ==========================

# Set colors
names(COLORS) <- unique(seurat_ob[[CELLCYCLE_NAME]])

# Plot UMAP
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

# ========================== Save results ==========================

# Save metadata
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
