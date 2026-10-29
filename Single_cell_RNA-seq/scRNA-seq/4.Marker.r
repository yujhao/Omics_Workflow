#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ==============================================================================
# Marker gene identification script
# Identify marker genes for each cluster and generate a heatmap
# ==============================================================================

# ========================== Configuration parameters ==========================

OUTPUT_DIR <- "4.Marker"
GROUPBY <- "clusters"
TOP_N <- 10

# ========================== Initialization ==========================

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# ========================== Marker gene identification ==========================

# Find all marker genes (positive markers only)
all.markers <- FindAllMarkers(scRNA_mnn, only.pos = TRUE)

# Get top N marker genes
top10 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = TOP_N, wt = avg_log2FC)

# ========================== Save results ==========================

# Save all marker genes
write.table(
    all.markers,
    file.path(OUTPUT_DIR, "all_Markers_of_each_clusters.xls"),
    col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE
)

# Save top N marker genes
write.table(
    top10,
    file.path(OUTPUT_DIR, "top10_Markers_of_each_clusters.xls"),
    col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE
)

# ========================== Heatmap visualization ==========================

# Data normalization
scRNA_mnn <- ScaleData(scRNA_mnn, features = row.names(scRNA_mnn))

# Plot marker gene heatmap
heatmap_plot <- DoHeatmap(
    object = scRNA_mnn,
    features = as.vector(top10$gene),
    group.by = GROUPBY,
    group.bar = TRUE,
    size = 3
) +
    theme(axis.text.y = element_text(size = 4))

# Save figures
ggsave(
    file.path(OUTPUT_DIR, "top10_marker_of_each_cluster_heatmap.pdf"),
    width = 12, height = 12, plot = heatmap_plot
)
ggsave(
    file.path(OUTPUT_DIR, "top10_marker_of_each_cluster_heatmap.png"),
    width = 12, height = 12, plot = heatmap_plot
)
