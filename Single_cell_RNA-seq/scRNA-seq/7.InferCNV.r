#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ==============================================================================
# InferCNV copy number variation analysis script
# Detect copy number variation in single-cell data
# ==============================================================================

library(infercnv)
library(plyranges)
library(ComplexHeatmap)

# ========================== Configuration parameters ==========================

OUTPUT_DIR <- "7.InferCNV"
GTF_FILE <- "genes.gtf"
REF_GROUP <- c("T_cells")  # Reference cell type
CUTOFF <- 0.1
NUM_THREADS <- 2

# ========================== Initialization ==========================

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# ========================== Data preparation ==========================

# Read GTF file to obtain gene position information
gtf <- plyranges::read_gff(GTF_FILE)
gene.chr <- gtf %>%
    plyranges::filter(type == "gene" & gene_name %in% rownames(scRNA_mnn)) %>%
    as.data.frame() %>%
    dplyr::select(gene_name, seqnames, start, end) %>%
    dplyr::distinct(gene_name, .keep_all = TRUE) %>%
    dplyr::mutate(seqnames = seqnames)

# Obtain expression matrix and cell annotations
count_mat <- GetAssayData(scRNA_mnn, "counts")
cellanno <- FetchData(scRNA_mnn, vars = "new_celltype") %>%
    tibble::rownames_to_column(var = "cellbarcode")

# Save temporary files
tempdir <- tempdir()
cnv_celltyping <- file.path(tempdir, "cnv_celltype_group.xls")
write.table(cellanno, cnv_celltyping, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

gene_order_f <- file.path(tempdir, "gene_order_file.xls")
write.table(gene.chr, gene_order_f, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

# ========================== Run InferCNV ==========================

# Create InferCNV object
infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = count_mat,
    annotations_file = cnv_celltyping,
    delim = "\t",
    gene_order_file = gene_order_f,
    ref_group_names = REF_GROUP
)

# Run analysis
infercnv_obj <- infercnv::run(
    infercnv_obj,
    cutoff = CUTOFF,
    analysis_mode = "subclusters",
    tumor_subcluster_pval = 0.05,
    hclust_method = "ward.D2",
    out_dir = OUTPUT_DIR,
    num_threads = NUM_THREADS,
    cluster_by_groups = TRUE,
    denoise = TRUE,
    HMM = TRUE
)

# ========================== Visualization ==========================

# Plot heatmap
pdf(file.path(OUTPUT_DIR, "heatmap.pdf"), width = 18, height = 12)
ComplexHeatmap::Heatmap(
    t(as.matrix(infercnv_obj@expr.data)),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    name = "CNV level",
    use_raster = TRUE,
    raster_quality = 4
)
dev.off()

# ========================== CNV score calculation ==========================

# Add CNV results to the Seurat object
scRNA_mnn[["CNV"]] <- CreateAssayObject(data = infercnv_obj@expr.data)

# Compute CNV level for each cell
infercnv_level <- apply(as.data.frame(t(infercnv_obj@expr.data)), 1, function(x) {
    x[is.na(x)] <- 0
    return(sum(x))
})
infercnv_level <- round(
    scales::rescale(infercnv_level / nrow(infercnv_obj@expr.data), c(1, 100)),
    0
)
infercnv_level <- infercnv_level[Cells(scRNA_mnn)]
scRNA_mnn@meta.data$cnv_level <- infercnv_level

# Plot CNV level figure
p1 <- FeaturePlot(scRNA_mnn, features = "cnv_level")
p2 <- VlnPlot(scRNA_mnn, features = "cnv_level", group.by = "orig.ident", pt.size = 0)

ggsave(file.path(OUTPUT_DIR, "featureplot.png"), p1, height = 5, width = 5)
ggsave(file.path(OUTPUT_DIR, "vlnplot.png"), p2, height = 4, width = 7)
