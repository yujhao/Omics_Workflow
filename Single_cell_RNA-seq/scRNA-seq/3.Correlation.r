#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ==============================================================================
# 相关性分析脚本
# 计算聚类间基因表达的相关性矩阵并生成热图
# ==============================================================================

# ========================== 配置参数 ==========================

OUTPUT_DIR <- "3.Correlation"
GROUPBY <- "clusters"

# ========================== 初始化 ==========================

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# ========================== 数据计算 ==========================

# 计算平均表达量
groupby_data <- AverageExpression(scRNA_mnn, group.by = GROUPBY)[["RNA"]]

# 保存标准化数据
data <- tibble::rownames_to_column(as.data.frame(groupby_data), var = "GeneID")
write.table(
    data,
    file.path(OUTPUT_DIR, paste0("normalized_data_groupby_", GROUPBY, ".xls")),
    quote = FALSE, row.names = FALSE, sep = "\t"
)

# ========================== 相关性计算 ==========================

# 添加聚类前缀
colnames(groupby_data) <- gsub('^', paste0("clusters", "_"), colnames(groupby_data))

# 计算 Pearson 相关系数矩阵
matrix <- cor(groupby_data, method = "pearson")

# ========================== 可视化 ==========================

# 动态计算图片尺寸
wid <- 5 + 1.5 * log2(length(colnames(data)))
hig <- 5 + 1.5 * log2(length(colnames(data)))

# 动态计算字体大小
fontsize_number <- 10.0 + 0.0001 * log2(length(colnames(data)))
fontsize_row <- 10.0 + 0.0001 * log2(length(colnames(data)))
fontsize_col <- 10.0 + 0.0001 * log2(length(colnames(data)))

# 绘制相关性热图
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

# 保存图片
ggsave(file.path(OUTPUT_DIR, "coefficient_heatmap.pdf"), coefficient, width = 8, height = 8)
ggsave(file.path(OUTPUT_DIR, "coefficient_heatmap.png"), coefficient, width = 8, height = 8)
