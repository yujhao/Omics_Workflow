#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ==============================================================================
# AddModuleScore 基因集评分脚本
# 为细胞添加自定义基因集的模块评分
# ==============================================================================

# ========================== 配置参数 ==========================

OUTPUT_DIR <- "6.AddModuleScore"
GENESET_FILE <- "geneset.addmodeulescore.xls"
SEED <- 1

# ========================== 初始化 ==========================

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# ========================== 基因集加载与匹配 ==========================

# 读取基因集
extra_gene <- read.delim(GENESET_FILE, check.names = FALSE)
if (dim(extra_gene)[2] == 1 && colnames(extra_gene)[1] == "gene") {
    colnames(extra_gene)[1] <- "extra"
}

# 转换为长格式
formated_extra_gene <- as.data.frame(tidyr::gather(extra_gene, key = "cluster", value = "GENE"))

# 匹配基因
match_results <- CaseMatch(
    search = as.vector(formated_extra_gene$GENE),
    match = rownames(scRNA_mnn)
)
match_results <- Filter(function(x) length(x) > 0, match_results)

# 记录未匹配的基因
filtered_gene <- formated_extra_gene$GENE[
    !formated_extra_gene$GENE %in% names(match_results) &
        formated_extra_gene$GENE != ""
]

if (length(filtered_gene) != 0) {
    filtered_gene <- as.data.frame(filtered_gene)
    colnames(filtered_gene) <- "Gene"
    write.table(
        filtered_gene,
        file.path(OUTPUT_DIR, "genes_not_matched.xls"),
        quote = FALSE, row.names = FALSE
    )
    print("有部分基因未匹配到，见：genes_not_matched.xls.")
}

# 使用匹配成功的基因
formated_extra_gene <- formated_extra_gene %>%
    dplyr::filter(GENE %in% names(match_results)) %>%
    rowwise() %>%
    dplyr::mutate(MATCHED = list(match_results[[GENE]])) %>%
    tidyr::unnest(cols = MATCHED) %>%
    dplyr::rename(folder_suffix = cluster, gene = MATCHED)

# ========================== 计算模块评分 ==========================

# 准备基因集列表
topn_markers <- formated_extra_gene
topn_markers2vis <- list()
for (clusterx in unique(topn_markers$folder_suffix)) {
    topn_markers2vis[[clusterx]] <- subset(topn_markers, folder_suffix == clusterx)$gene
}

# 计算模块评分
scRNA_mnn <- AddModuleScore(
    scRNA_mnn,
    features = topn_markers2vis,
    name = names(topn_markers2vis),
    seed = SEED
)

# 修正列名
colnames(scRNA_mnn@meta.data)[
    (dim(scRNA_mnn[[]])[2] - length(topn_markers2vis) + 1):dim(scRNA_mnn[[]])[2]
] <- names(topn_markers2vis)

# ========================== 保存评分结果 ==========================

matrix <- scRNA_mnn@meta.data %>%
    dplyr::rename("Barcode" = "rawbc") %>%
    dplyr::select(Barcode, "clusters", !!names(topn_markers2vis)) %>%
    dplyr::arrange(!!sym("clusters"))

write.table(
    matrix,
    quote = FALSE, sep = "\t", row.names = FALSE,
    file.path(OUTPUT_DIR, "addmodeulescore_plot.xls")
)

# ========================== 可视化 ==========================

# 箱线图/小提琴图
for (i in colnames(matrix)[-(1:2)]) {
    p <- ggstatsplot::ggbetweenstats(
        matrix,
        x = !!sym("clusters"),
        y = !!sym(i),
        plot.type = "boxviolin",
        results.subtitle = FALSE,
        messages = FALSE,
        pairwise.comparisons = FALSE,
        mean.label.size = 0,
        centrality.plotting = FALSE,
        ylab = paste(i)
    ) +
        theme(
            axis.text.x = element_text(size = 8, colour = "black", angle = 30, vjust = 0.85, hjust = 0.75),
            axis.text.y = element_text(size = 8, colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black")
        )

    savenames <- paste0(i, "_ggstatsplot_violin_boxplot")
    ggsave(file.path(OUTPUT_DIR, paste0(savenames, ".pdf")), plot = p, bg = "white")
    ggsave(file.path(OUTPUT_DIR, paste0(savenames, ".png")), plot = p, bg = "white")
}

# FeaturePlot
lapply(names(topn_markers2vis), function(x) {
    plot <- Seurat::FeaturePlot(
        scRNA_mnn,
        features = x,
        cols = c("grey", "red"),
        split.by = NULL,
        reduction = "umap",
        ncol = 1,
        pt.size = 0.5,
        order = TRUE
    ) +
        theme(
            plot.title = element_text(hjust = 0.5),
            legend.position = "right"
        ) +
        theme(aspect.ratio = 1/1)

    ggsave(
        file.path(OUTPUT_DIR, paste0(x, "_score_featureplot.pdf")),
        plot = plot, width = 6, height = 5, limitsize = FALSE
    )
    ggsave(
        file.path(OUTPUT_DIR, paste0(x, "_score_featureplot.png")),
        plot = plot, width = 6, height = 5, bg = "white"
    )
})
