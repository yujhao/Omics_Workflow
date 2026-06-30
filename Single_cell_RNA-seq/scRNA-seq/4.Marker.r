#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ==============================================================================
# Marker 基因鉴定脚本
# 鉴定各聚类的 Marker 基因并生成热图
# ==============================================================================

# ========================== 配置参数 ==========================

OUTPUT_DIR <- "4.Marker"
GROUPBY <- "clusters"
TOP_N <- 10

# ========================== 初始化 ==========================

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# ========================== Marker 基因鉴定 ==========================

# 查找所有 Marker 基因（仅正标记）
all.markers <- FindAllMarkers(scRNA_mnn, only.pos = TRUE)

# 获取 Top N Marker 基因
top10 <- all.markers %>%
    group_by(cluster) %>%
    top_n(n = TOP_N, wt = avg_log2FC)

# ========================== 保存结果 ==========================

# 保存所有 Marker 基因
write.table(
    all.markers,
    file.path(OUTPUT_DIR, "all_Markers_of_each_clusters.xls"),
    col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE
)

# 保存 Top N Marker 基因
write.table(
    top10,
    file.path(OUTPUT_DIR, "top10_Markers_of_each_clusters.xls"),
    col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE
)

# ========================== 热图可视化 ==========================

# 数据标准化
scRNA_mnn <- ScaleData(scRNA_mnn, features = row.names(scRNA_mnn))

# 绘制 Marker 基因热图
heatmap_plot <- DoHeatmap(
    object = scRNA_mnn,
    features = as.vector(top10$gene),
    group.by = GROUPBY,
    group.bar = TRUE,
    size = 3
) +
    theme(axis.text.y = element_text(size = 4))

# 保存图片
ggsave(
    file.path(OUTPUT_DIR, "top10_marker_of_each_cluster_heatmap.pdf"),
    width = 12, height = 12, plot = heatmap_plot
)
ggsave(
    file.path(OUTPUT_DIR, "top10_marker_of_each_cluster_heatmap.png"),
    width = 12, height = 12, plot = heatmap_plot
)
