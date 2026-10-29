#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ==============================================================================
# Monocle2 pseudotime analysis script
# Perform single-cell trajectory inference and pseudotime analysis
# ==============================================================================

library(monocle)

# ========================== Configuration parameters ==========================

OUTPUT_DIR <- "8.Monocle2"
MIN_EXPR <- 1
MIN_CELLS_EXPRESSED <- 10
QVAL_CUTOFF <- 0.01
NUM_CLUSTERS <- 4
BRANCH_POINT <- 2
CORES <- 4

# ========================== Initialization ==========================

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# ========================== CDS object creation ==========================

# Convert to CellDataSet
cds <- as.CellDataSet(scRNA_mnn)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = MIN_EXPR)

# Filter expressed genes
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed > MIN_CELLS_EXPRESSED))

# ========================== Differential expression analysis ==========================

# Detect differential genes between clusters
clustering_DEGs <- differentialGeneTest(
    cds[expressed_genes, ],
    fullModelFormulaStr = "~clusters",
    cores = CORES
)

# Update qval
featureData(cds)@data[rownames(clustering_DEGs), "qval"] <- clustering_DEGs$qval

# Select ordering genes
ordering_genes <- row.names(subset(clustering_DEGs, qval < QVAL_CUTOFF))

# ========================== Trajectory construction ==========================

# Set ordering filter
gbm_cds <- setOrderingFilter(cds, ordering_genes = ordering_genes)

# Save ordering gene plot
p <- plot_ordering_genes(gbm_cds)
ggsave(file.path(OUTPUT_DIR, "ordering_genes.pdf"), p, width = 8, height = 4)
ggsave(file.path(OUTPUT_DIR, "ordering_genes.png"), p, width = 8, height = 4)

# Dimensionality reduction and ordering
gbm_cds <- reduceDimension(
    gbm_cds,
    max_components = 2,
    verbose = TRUE,
    check_duplicates = FALSE,
    num_dim = 10
)
gbm_cds <- orderCells(gbm_cds, reverse = FALSE)

# ========================== Trajectory visualization ==========================

# Color by State
p <- plot_cell_trajectory(
    gbm_cds,
    color_by = "State",
    cell_size = 1.5,
    show_branch_points = TRUE
) + scale_color_simpsons()
ggsave(file.path(OUTPUT_DIR, "cell_trajectory_state.pdf"), p, width = 8, height = 4)
ggsave(file.path(OUTPUT_DIR, "cell_trajectory_state.png"), p, width = 8, height = 4)

# Color by clusters
p <- plot_cell_trajectory(
    gbm_cds,
    color_by = "clusters",
    cell_size = 1.5,
    show_branch_points = TRUE
)
ggsave(file.path(OUTPUT_DIR, "cell_trajectory_clusters.pdf"), p, width = 8, height = 4)
ggsave(file.path(OUTPUT_DIR, "cell_trajectory_clusters.png"), p, width = 8, height = 4)

# Facet by clusters
p <- plot_cell_trajectory(
    gbm_cds,
    color_by = "clusters",
    cell_size = 1.5,
    show_branch_points = TRUE
) + facet_wrap(~clusters)
ggsave(file.path(OUTPUT_DIR, "cell_trajectory_facet.pdf"), p, width = 8, height = 4)
ggsave(file.path(OUTPUT_DIR, "cell_trajectory_facet.png"), p, width = 8, height = 4)

# Color by Pseudotime
p <- plot_cell_trajectory(
    gbm_cds,
    color_by = "Pseudotime",
    show_branch_points = FALSE
) + scale_colour_viridis_c(option = "inferno")
ggsave(file.path(OUTPUT_DIR, "cell_trajectory_pseudotime.pdf"), p, width = 8, height = 4)
ggsave(file.path(OUTPUT_DIR, "cell_trajectory_pseudotime.png"), p, width = 8, height = 4)

# ========================== Heatmap analysis ==========================

# Prepare gene data
genes <- as.factor(subset(gbm_cds@featureData@data, use_for_ordering == TRUE)$gene_short_name)
to_be_tested <- row.names(subset(fData(gbm_cds), gene_short_name %in% levels(genes)))
gbm_cds <- gbm_cds[to_be_tested, ]
varMetadata(gbm_cds)[, 1] <- rownames(varMetadata(gbm_cds))
gbm_cds@featureData@varMetadata[, 1] <- rownames(gbm_cds@featureData@varMetadata)

# Plot pseudotime heatmap
p <- plot_pseudotime_heatmap(
    gbm_cds,
    cores = 1,
    cluster_rows = TRUE,
    num_clusters = NUM_CLUSTERS,
    show_rownames = FALSE,
    return_heatmap = TRUE
)
ggsave(file.path(OUTPUT_DIR, "pseudotime_heatmap.pdf"), p$ph_res, width = 8, height = 4)
ggsave(file.path(OUTPUT_DIR, "pseudotime_heatmap.png"), p$ph_res, width = 8, height = 4)

# ========================== Branch point analysis ==========================

# Gene module clustering
gene_clusters <- cutree(p$tree_row, k = NUM_CLUSTERS)
gene_clustering <- data.frame(gene_clusters)
gene_clustering[, 1] <- as.character(gene_clustering[, 1])
colnames(gene_clustering) <- "gene_module"

# Branch expression analysis
BEAM_res <- BEAM(gbm_cds, branch_point = BRANCH_POINT, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval), ]
BEAM_res <- BEAM_res[, c("gene_short_name", "pval", "qval")]

# Branch heatmap
p <- plot_genes_branched_heatmap(
    gbm_cds[row.names(subset(BEAM_res, qval < 1e-4)), ],
    branch_point = BRANCH_POINT,
    num_clusters = NUM_CLUSTERS,
    cores = 1,
    use_gene_short_name = TRUE,
    show_rownames = FALSE,
    return_heatmap = TRUE
)
ggsave(
    file.path(OUTPUT_DIR, "pseudotime_heatmap_branchtime.pdf"),
    plot = p$ph_res, height = 7, width = 7, bg = "white"
)
ggsave(
    file.path(OUTPUT_DIR, "pseudotime_heatmap_branchtime.png"),
    plot = p$ph_res, height = 7, width = 7, dpi = 1000, bg = "white"
)

# ========================== Specified gene visualization ==========================

# Select genes of interest
genes_of_interest <- c("ISG15", "PARK7", "GLUL")
to_be_tested_sub <- row.names(subset(fData(gbm_cds), gene_short_name %in% genes_of_interest))

# Gene jitter plot
p <- plot_genes_jitter(
    gbm_cds[to_be_tested_sub, ],
    grouping = "State",
    min_expr = 0.1,
    color_by = "State",
    cell_size = 1
) + scale_color_simpsons()
ggsave(file.path(OUTPUT_DIR, "genes_jitter.pdf"), p, width = 8, height = 4)
ggsave(file.path(OUTPUT_DIR, "genes_jitter.png"), p, width = 8, height = 4)

# Gene expression along pseudotime
p <- plot_genes_in_pseudotime(
    gbm_cds[to_be_tested_sub, ],
    color_by = "clusters",
    cell_size = 1,
    ncol = 1
) + scale_color_simpsons()
ggsave(file.path(OUTPUT_DIR, "genes_in_pseudotime.pdf"), p, width = 8, height = 4)
ggsave(file.path(OUTPUT_DIR, "genes_in_pseudotime.png"), p, width = 8, height = 4)

# Branched pseudotime
new_cds <- buildBranchCellDataSet(
    gbm_cds[to_be_tested_sub, ],
    branch_point = BRANCH_POINT,
    progenitor_method = "duplicate"
)

cell_fate1 <- unique(pData(new_cds)[which(pData(new_cds)$Branch == unique(pData(new_cds)$Branch)[1]), ]$State)
cell_fate2 <- unique(pData(new_cds)[which(pData(new_cds)$Branch == unique(pData(new_cds)$Branch)[2]), ]$State)
branch_labels <- c(
    paste("State", paste(sort(setdiff(cell_fate1, cell_fate2)), collapse = "-")),
    paste("State", paste(sort(setdiff(cell_fate2, cell_fate1)), collapse = "-"))
)

p <- plot_genes_branched_pseudotime(
    gbm_cds[to_be_tested_sub, ],
    color_by = "clusters",
    branch_point = BRANCH_POINT,
    cell_size = 1,
    ncol = 1,
    branch_labels = branch_labels
) + scale_color_simpsons()
ggsave(file.path(OUTPUT_DIR, "genes_branched_pseudotime.pdf"), p, width = 8, height = 4)
ggsave(file.path(OUTPUT_DIR, "genes_branched_pseudotime.png"), p, width = 8, height = 4)

# Single-gene expression gradient
p <- plot_cell_trajectory(
    gbm_cds,
    markers = "GLUL",
    use_color_gradient = TRUE,
    show_branch_points = FALSE,
    show_tree = FALSE,
    cell_size = 1.5
) +
    theme(legend.text = element_text(size = 10)) +
    scale_color_gradientn(colours = c("grey", "yellow", "red"))
ggsave(file.path(OUTPUT_DIR, "genes_gradient_GLUL.pdf"), p, width = 8, height = 4)
ggsave(file.path(OUTPUT_DIR, "genes_gradient_GLUL.png"), p, width = 8, height = 4)
