#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ==============================================================================
# Correlation analysis script
# Compute the correlation matrix of inter-cluster gene expression and generate a heatmap
# ==============================================================================

# ========================== Configuration parameters ==========================

OUTPUT_DIR <- "3.Correlation"
GROUPBY <- "clusters"

# ========================== Initialization ==========================

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# ========================== Data computation ==========================

# Compute average expression
groupby_data <- AverageExpression(scRNA_mnn, group.by = GROUPBY)[["RNA"]]

# Save normalized data
data <- tibble::rownames_to_column(as.data.frame(groupby_data), var = "GeneID")
write.table(
    data,
    file.path(OUTPUT_DIR, paste0("normalized_data_groupby_", GROUPBY, ".xls")),
    quote = FALSE, row.names = FALSE, sep = "\t"
)

# ========================== Correlation computation ==========================

# Add cluster prefix
colnames(groupby_data) <- gsub('^', paste0("clusters", "_"), colnames(groupby_data))

# Compute Pearson correlation coefficient matrix
matrix <- cor(groupby_data, method = "pearson")

# ========================== Visualization ==========================

# Dynamically compute figure dimensions
wid <- 5 + 1.5 * log2(length(colnames(data)))
hig <- 5 + 1.5 * log2(length(colnames(data)))

# Dynamically compute font size
fontsize_number <- 10.0 + 0.0001 * log2(length(colnames(data)))
fontsize_row <- 10.0 + 0.0001 * log2(length(colnames(data)))
fontsize_col <- 10.0 + 0.0001 * log2(length(colnames(data)))

# Plot correlation heatmap
coefficient <- pheatmap::pheatmap(
    matrix,
    display_numbers = FALSE,
    border_color = "white",
    scale = "none",
    fontsize_number = fontsize_number,
    number_format = "%.1f",
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    number_color = "black",
    angle_col = 45
)

# Save figures
ggsave(file.path(OUTPUT_DIR, "coefficient_heatmap.pdf"), coefficient, width = 8, height = 8)
ggsave(file.path(OUTPUT_DIR, "coefficient_heatmap.png"), coefficient, width = 8, height = 8)
