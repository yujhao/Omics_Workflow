#!/usr/bin/env Rscript
# ============================================================================
# GSEA (Gene Set Enrichment Analysis) analysis script
# Function: Perform GSEA enrichment analysis based on differential expression analysis
# Input: Seurat object + gene set GMT file
# Output: GSEA enrichment results, visualization plots
# ============================================================================

# ========================== Load packages ==========================

suppressPackageStartupMessages({
    library(Seurat)
    library(clusterProfiler)
    library(GSEABase)
    library(ggplot2)
    library(enrichplot)
    library(DOSE)
    library(pheatmap)
    library(RColorBrewer)
})

# ========================== Configuration parameters ==========================

OUTPUT_DIR <- "15.GSEA"
DATA_FILE <- "data_ob_v3.rds"
GMT_FILE <- "utils/reference/refdata-gex-GRCm39-2024-A/gmt/mouse_2024_kegg.gmt"

# Analysis parameters
IDENT_COLUMN <- "clusters"
GROUP1 <- "1"
GROUP2 <- "2"
PVAL_CUTOFF <- 0.05
LOGFC_CUTOFF <- 0.5
MIN_GSSIZE <- 1
PVALUE_CUTOFF <- 0.5

# Visualization parameters
GSEA_TERM <- "Influenza A(mmu05164)"
PLOT_WIDTH <- 10
PLOT_HEIGHT <- 8
PLOT_DPI <- 300

# ========================== Initialization ==========================

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# ========================== Data loading ==========================

# Read Seurat object
pbmc <- readRDS(DATA_FILE)
Idents(pbmc) <- IDENT_COLUMN

# Read gene set file
gmt <- read.gmt(GMT_FILE)
gmt$term <- gsub("KEGG_", "", gmt$term)

# ========================== Cell subpopulation extraction ==========================

scedata <- pbmc
cells_sub <- subset(scedata@meta.data, clusters %in% c(GROUP1, GROUP2))
scRNA_sub <- subset(scedata, cells = row.names(cells_sub))

# ========================== Differential expression analysis ==========================

# Compute per-gene logFC between groups
sub.markers <- FindMarkers(
    scRNA_sub,
    group.by = IDENT_COLUMN,
    ident.1 = GROUP1,
    ident.2 = GROUP2
)

# Filter significantly differential genes
deg <- subset(sub.markers, p_val_adj < PVAL_CUTOFF & avg_log2FC > LOGFC_CUTOFF)

# ========================== GSEA enrichment analysis ==========================

# Extract log2FC values and rank them
genelist <- deg$avg_log2FC
names(genelist) <- rownames(deg)
genelist <- sort(genelist, decreasing = TRUE)

# Run GSEA
gsea <- GSEA(
    genelist,
    TERM2GENE = gmt,
    minGSSize = MIN_GSSIZE,
    pvalueCutoff = PVALUE_CUTOFF
)

# ========================== Visualization ==========================

# Plot GSEA figure
p <- gseaplot2(
    gsea,
    geneSetID = GSEA_TERM,
    subplots = c(1, 2, 3),
    color = c("#dd2222"),
    ES_geom = "line",
    pvalue_table = TRUE,
    title = FALSE
)

# Save plot
ggsave(
    file.path(OUTPUT_DIR, "GSEA.pdf"),
    p,
    width = PLOT_WIDTH,
    height = PLOT_HEIGHT,
    dpi = PLOT_DPI
)
