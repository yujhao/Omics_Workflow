#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ==============================================================================
# hdWGCNA module preservation analysis
# Projects the co-expression modules of a reference dataset onto a query dataset
# and evaluates module preservation with permutation based statistics.
# Follows the official hdWGCNA tutorial.
# ==============================================================================

# ========================== Load packages ==========================

suppressPackageStartupMessages({
    library(Seurat)
    library(hdWGCNA)
    library(tidyverse)
    library(patchwork)
})

# ========================== Configuration parameters ==========================

OUTPUT_DIR <- "4.ModulePreservation"

REF_FILE <- "data/ref_hdWGCNA_object.rds"      # reference hdWGCNA object
QUERY_FILE <- "data/query_hdWGCNA_object.rds"  # query hdWGCNA object
GENE_MAPPING_FILE <- NULL                      # optional TSV with query/ref gene pairs,
                                               # required only for cross-species analysis

GROUP_BY <- "clusters"
ASSAY <- "RNA"
REDUCTION <- "umap"
N_PERMUTATIONS <- 200

# ========================== Initialization ==========================

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(12345)

# ========================== Data loading ==========================

seurat_ref <- readRDS(REF_FILE)
seurat_query <- readRDS(QUERY_FILE)

# Cross-species analysis: map query genes onto the reference gene space
gene_mapping <- NULL
if (!is.null(GENE_MAPPING_FILE)) {
    gene_mapping <- read.delim(GENE_MAPPING_FILE, sep = "\t", header = TRUE) %>%
        dplyr::filter(gene_ref != "--") %>%
        dplyr::select(gene_query, gene_ref)
}

# ========================== Module projection ==========================

seurat_pro <- ProjectModules(
    seurat_query,
    seurat_ref = seurat_ref,
    group.by.vars = GROUP_BY,
    gene_mapping = gene_mapping,
    overlap_proportion = 0.5,
    vars.to.regress = c("nCount_RNA", "nFeature_RNA"),
    scale.model.use = "linear",
    wgcna_name = "hdWGCNA",
    wgcna_name_proj = "hdWGCNA",
    assay = ASSAY
)

seurat_pro <- ModuleConnectivity(
    seurat_pro,
    group.by = GROUP_BY,
    group_name = levels(seurat_pro@meta.data[[GROUP_BY]])
)

seurat_pro <- ModuleExprScore(seurat_pro, method = "UCell")

# Projected module eigengenes on the query embedding
plot_list <- ModuleFeaturePlot(seurat_pro, features = "hMEs",
                               order = TRUE, restrict_range = FALSE,
                               reduction = REDUCTION)
ggsave(wrap_plots(plot_list, ncol = 4), filename = file.path(OUTPUT_DIR, "hMEs_project.png"),
       width = 10, height = 5, bg = "white")

saveRDS(seurat_pro, file = file.path(OUTPUT_DIR, "ProjectModules.rds"))

# ========================== Module preservation ==========================

seurat_ref_expr <- SetDatExpr(
    seurat_ref,
    group_name = levels(seurat_ref@meta.data[[GROUP_BY]]),
    group.by = GROUP_BY,
    assay = ASSAY
)

seurat_pro_expr <- SetDatExpr(
    seurat_pro,
    group_name = levels(seurat_pro@meta.data[[GROUP_BY]]),
    group.by = GROUP_BY,
    assay = ASSAY
)

seurat_pre <- ModulePreservation(
    seurat_pro_expr,
    seurat_ref = seurat_ref_expr,
    name = "preserving",
    gene_mapping = gene_mapping,
    verbose = 3,
    n_permutations = N_PERMUTATIONS
)

# ========================== Results ==========================

mod_pres <- GetModulePreservation(seurat_pre, "preserving")
write.table(mod_pres$Z, file = file.path(OUTPUT_DIR, "mod_preservation_Z_summary.xls"),
            sep = "\t", quote = FALSE, row.names = TRUE)
write.table(mod_pres$obs, file = file.path(OUTPUT_DIR, "mod_preservation_obs_summary.xls"),
            sep = "\t", quote = FALSE, row.names = TRUE)

p <- wrap_plots(PlotModulePreservation(seurat_pre, name = "preserving",
                                       statistics = "summary"), ncol = 2)
ggsave(p, filename = file.path(OUTPUT_DIR, "summary.png"), width = 8, height = 4, bg = "white")

p <- wrap_plots(PlotModulePreservation(seurat_pre, name = "preserving",
                                       statistics = "rank"), ncol = 2)
ggsave(p, filename = file.path(OUTPUT_DIR, "rank.png"), width = 8, height = 8, bg = "white")

saveRDS(seurat_pre, file = file.path(OUTPUT_DIR, "ModulePreservation.rds"))

message("Module preservation analysis finished, results written to ", OUTPUT_DIR)
