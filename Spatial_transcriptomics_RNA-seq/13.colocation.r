# Load required packages
library(Seurat)
library(ggplot2)
library(patchwork)

#' Plot a spatial colocalization blend figure (Feature Blend) of two genes (or metadata columns)
#'
#' @param object Seurat object (must contain spatial image)
#' @param features character vector of length 2: gene names or metadata column names
#' @param cols color vector of length 3 (background, gene1, gene2), default c("lightgrey","red","darkblue")
#' @param blend.threshold blend threshold, 0~1, default 0.2
#' @param image.alpha image transparency, default 1
#' @param pt.size point size (for blend and spatial plots), default 1.2
#' @param crop whether to crop the image, default TRUE
#' @param output_dir output directory (optional; if not provided, print the combined plot directly)
#'
#' @return returns a patchwork combined plot object, and can be saved as PDF/PNG
#' @export
#'
feature_blend_func <- function(object,
                               features,
                               cols = c("lightgrey", "red", "darkblue"),
                               blend.threshold = 0.2,
                               image.alpha = 1,
                               pt.size = 1.2,
                               crop = TRUE,
                               output_dir = NULL) {
  
  # 1. Check input
  stopifnot(length(features) == 2)
  
  # 2. Match features (genes or metadata)
  genes <- Seurat::CaseMatch(features, rownames(object))
  if (length(genes) == 0) {
    genes <- Seurat::CaseMatch(features, colnames(object[[]]))
  }
  if (length(genes) != 2) stop("Could not find both features in object.")
  
  # 3. Generate a standard FeaturePlot (blend = TRUE) to display the blended scatter plot
  p_blend <- Seurat::FeaturePlot(object,
                                 features = genes,
                                 blend = TRUE,
                                 cols = cols,
                                 blend.threshold = blend.threshold,
                                 pt.size = pt.size,
                                 combine = TRUE)
  
  # 4. Generate three spatial subplots (showing gene1, gene2, and blend separately)
  # First obtain the three separate layers with blend = FALSE (to extract the color mapping)
  p_list <- Seurat::FeaturePlot(object,
                                features = genes,
                                blend = TRUE,
                                cols = cols,
                                blend.threshold = blend.threshold,
                                pt.size = pt.size,
                                combine = FALSE)  # returns a list with 4 elements: gene1, gene2, blend, legend
  
  # Extract the color mapping of each subplot (ensure SpatialDimPlot uses the same colors)
  get_cols <- function(plot_obj, group) {
    g <- ggplot_build(plot_obj)
    df <- unique(plot_obj$data[, c(group, "colour")])
    cols <- df$colour
    names(cols) <- as.character(df[[group]])
    return(cols)
  }
  
  # color of gene1
  col1 <- get_cols(p_list[[1]], genes[1])
  # color of gene2
  col2 <- get_cols(p_list[[2]], genes[2])
  # color of the blend channel (group name is "gene1_gene2")
  blend_group <- paste(genes[1], genes[2], sep = "_")
  col_blend <- get_cols(p_list[[3]], blend_group)
  
  # 5. Add gene expression (or metadata) as grouping variables to meta.data for use by SpatialDimPlot
  object <- AddMetaData(object,
                        metadata = cbind(p_list[[1]]$data[, genes[1], drop = FALSE],
                                         p_list[[2]]$data[, genes[2], drop = FALSE],
                                         p_list[[3]]$data[, blend_group, drop = FALSE]))
  
  # 6. Plot the three spatial plots separately (no legend)
  p_sp1 <- SpatialDimPlot(object,
                          group.by = genes[1],
                          cols = col1,
                          image.alpha = image.alpha,
                          pt.size = pt.size,
                          crop = crop) +
    ggtitle(genes[1]) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16)) &
    NoLegend()
  
  p_sp2 <- SpatialDimPlot(object,
                          group.by = genes[2],
                          cols = col2,
                          image.alpha = image.alpha,
                          pt.size = pt.size,
                          crop = crop) +
    ggtitle(genes[2]) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16)) &
    NoLegend()
  
  p_sp3 <- SpatialDimPlot(object,
                          group.by = blend_group,
                          cols = col_blend,
                          image.alpha = image.alpha,
                          pt.size = pt.size,
                          crop = crop) +
    ggtitle(blend_group) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16)) &
    NoLegend()
  
  # 7. Combine: first row is the standard blend scatter plot, second row is the three spatial plots
  p_top <- p_blend
  p_bottom <- wrap_plots(p_sp1, p_sp2, p_sp3, p_list[[4]], nrow = 1, widths = c(1,1,1,0.9))
  p_combined <- wrap_plots(p_top, p_bottom, nrow = 2)
  
  # 8. Save or return
  if (!is.null(output_dir)) {
    ggsave(filename = file.path(output_dir, paste0("FeatureBlend_", genes[1], "_", genes[2], ".png")),
           plot = p_combined, width = 16, height = 9, dpi = 300)
  }
  
  return(p_combined)
}