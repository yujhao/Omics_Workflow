#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ==============================================================================
# MetaNeighbor analysis
# Quantifies transcriptomic similarity of cell types across datasets / species
# and tests whether cell types are replicable (reproducibility AUROC).
#
# Procedure 1 : variable genes -> AUROC heatmap and best hits
# Procedure 2 : train a signature on the data and re-score with the frozen model
# ==============================================================================

# ========================== Load packages ==========================

suppressPackageStartupMessages({
    library(MetaNeighbor)
    library(Seurat)
    library(SeuratDisk)
    library(SingleCellExperiment)
    library(tibble)
})

# ========================== Configuration parameters ==========================

OUTPUT_DIR <- "5.MetaNeighbor"

REF_FILE <- "data/ref_data_ob.rds"        # reference dataset
QUERY_FILE <- "data/query_data_ob.rds"    # query dataset
INFORMAT <- "rds"                         # rds or h5seurat

REF_LABEL <- "SpeciesA"                   # dataset label added to the cluster names
QUERY_LABEL <- "SpeciesB"
CLUSTER_COL <- "clusters"                 # metadata column holding the cluster id
SPECIES_COL <- "species"                  # metadata column holding the dataset label

TOPHITS_THRESHOLD <- 0.6

# ========================== Initialization ==========================

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ========================== Data loading ==========================

LoadData <- function(file, informat = "rds") {
    if (informat == "h5seurat") {
        return(SeuratDisk::LoadH5Seurat(file, verbose = FALSE))
    }
    readRDS(file)
}

ref_obj <- UpdateSeuratObject(LoadData(REF_FILE, INFORMAT))
query_obj <- UpdateSeuratObject(LoadData(QUERY_FILE, INFORMAT))

# Both datasets must share the same feature space
query_obj <- query_obj[rownames(ref_obj), ]

# Prefix cluster ids so the two datasets stay distinguishable after merging
ref_obj@meta.data$all_clusters <- paste0(REF_LABEL, "_", ref_obj@meta.data[[CLUSTER_COL]])
query_obj@meta.data$all_clusters <- paste0(QUERY_LABEL, "_", query_obj@meta.data[[CLUSTER_COL]])
ref_obj@meta.data[[SPECIES_COL]] <- REF_LABEL
query_obj@meta.data[[SPECIES_COL]] <- QUERY_LABEL

data_obj <- merge(ref_obj, query_obj)
data_obj@meta.data$all_clusters <- factor(data_obj@meta.data$all_clusters)

sce_obj <- as.SingleCellExperiment(data_obj)

# ========================== Procedure 1: variable genes ==========================

proc_dir <- file.path(OUTPUT_DIR, "Procedure_1")
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)

global_hvgs <- variableGenes(dat = sce_obj, exp_labels = sce_obj[[SPECIES_COL]])

aurocs <- MetaNeighborUS(var_genes = global_hvgs,
                         dat = sce_obj,
                         study_id = sce_obj[[SPECIES_COL]],
                         cell_type = sce_obj$all_clusters,
                         fast_version = TRUE)

write.table(tibble::rownames_to_column(as.data.frame(aurocs), "all_clusters"),
            file.path(proc_dir, "aurocs.xls"), sep = "\t", quote = FALSE, row.names = FALSE)

png(file.path(proc_dir, "heatmap.png"), width = 6, height = 6, res = 700, units = "in")
plotHeatmap(aurocs, margins = c(12, 12), cex = 0.8)
dev.off()

top_hits <- topHits(aurocs, dat = sce_obj, study_id = sce_obj[[SPECIES_COL]],
                    cell_type = sce_obj$all_clusters, threshold = TOPHITS_THRESHOLD)
write.table(top_hits, file.path(proc_dir, "topHits.xls"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ========================== Procedure 2: pretrained model ==========================

proc2_dir <- file.path(OUTPUT_DIR, "Procedure_2")
dir.create(proc2_dir, recursive = TRUE, showWarnings = FALSE)

pretrained_model <- MetaNeighbor::trainModel(var_genes = rownames(sce_obj),
                                             dat = sce_obj,
                                             study_id = sce_obj[[SPECIES_COL]],
                                             cell_type = sce_obj$all_clusters)

write.table(pretrained_model, file.path(proc2_dir, "pretrained_model.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE)

aurocs_pre <- MetaNeighborUS(trained_model = pretrained_model,
                             dat = sce_obj,
                             study_id = sce_obj[[SPECIES_COL]],
                             cell_type = sce_obj$all_clusters,
                             fast_version = TRUE)

png(file.path(proc2_dir, "heatmap_pretrained.png"), width = 8, height = 8, res = 600, units = "in")
plotHeatmapPretrained(aurocs_pre, margins = c(12, 12))
dev.off()

message("MetaNeighbor analysis finished, results written to ", OUTPUT_DIR)
