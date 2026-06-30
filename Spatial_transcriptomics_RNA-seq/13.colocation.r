# 加载必要的包
library(Seurat)
library(ggplot2)
library(patchwork)

#' 绘制两个基因（或元数据列）的空间共定位混合图（Feature Blend）
#'
#' @param object Seurat 对象（需包含空间图像）
#' @param features 长度为2的字符向量，基因名或 metadata 列名
#' @param cols 颜色向量，长度3（背景、基因1、基因2），默认 c("lightgrey","red","darkblue")
#' @param blend.threshold 混合阈值，0~1，默认 0.2
#' @param image.alpha 图像透明度，默认 1
#' @param pt.size 点大小（用于混合图及空间图），默认 1.2
#' @param crop 是否裁剪图像，默认 TRUE
#' @param output_dir 输出目录（可选，若不提供则直接打印组合图）
#'
#' @return 返回一个 patchwork 组合图对象，并可保存为 PDF/PNG
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
  
  # 1. 检查输入
  stopifnot(length(features) == 2)
  
  # 2. 匹配特征（基因或 metadata）
  genes <- Seurat::CaseMatch(features, rownames(object))
  if (length(genes) == 0) {
    genes <- Seurat::CaseMatch(features, colnames(object[[]]))
  }
  if (length(genes) != 2) stop("Could not find both features in object.")
  
  # 3. 生成普通 FeaturePlot (blend = TRUE) 用于展示混合散点图
  p_blend <- Seurat::FeaturePlot(object,
                                 features = genes,
                                 blend = TRUE,
                                 cols = cols,
                                 blend.threshold = blend.threshold,
                                 pt.size = pt.size,
                                 combine = TRUE)
  
  # 4. 生成三个空间子图（分别显示基因1、基因2、混合）
  # 先获取 blend = FALSE 的三个单独图层（用于提取颜色映射）
  p_list <- Seurat::FeaturePlot(object,
                                features = genes,
                                blend = TRUE,
                                cols = cols,
                                blend.threshold = blend.threshold,
                                pt.size = pt.size,
                                combine = FALSE)  # 返回列表，含4个元素：基因1、基因2、混合、图例
  
  # 提取每个子图的颜色映射（确保 SpatialDimPlot 使用相同颜色）
  get_cols <- function(plot_obj, group) {
    g <- ggplot_build(plot_obj)
    df <- unique(plot_obj$data[, c(group, "colour")])
    cols <- df$colour
    names(cols) <- as.character(df[[group]])
    return(cols)
  }
  
  # 基因1的颜色
  col1 <- get_cols(p_list[[1]], genes[1])
  # 基因2的颜色
  col2 <- get_cols(p_list[[2]], genes[2])
  # 混合通道的颜色（组名为 "gene1_gene2"）
  blend_group <- paste(genes[1], genes[2], sep = "_")
  col_blend <- get_cols(p_list[[3]], blend_group)
  
  # 5. 将基因表达量（或 metadata）作为分组变量添加到 meta.data，以便 SpatialDimPlot 使用
  object <- AddMetaData(object,
                        metadata = cbind(p_list[[1]]$data[, genes[1], drop = FALSE],
                                         p_list[[2]]$data[, genes[2], drop = FALSE],
                                         p_list[[3]]$data[, blend_group, drop = FALSE]))
  
  # 6. 分别绘制三个空间图（无图例）
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
  
  # 7. 合并：第一行是普通混合散点图，第二行是三个空间图
  p_top <- p_blend
  p_bottom <- wrap_plots(p_sp1, p_sp2, p_sp3, p_list[[4]], nrow = 1, widths = c(1,1,1,0.9))
  p_combined <- wrap_plots(p_top, p_bottom, nrow = 2)
  
  # 8. 保存或返回
  if (!is.null(output_dir)) {
    ggsave(filename = file.path(output_dir, paste0("FeatureBlend_", genes[1], "_", genes[2], ".png")),
           plot = p_combined, width = 16, height = 9, dpi = 300)
  }
  
  return(p_combined)
}