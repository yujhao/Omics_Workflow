#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ==============================================================================
# GSVA (Gene Set Variation Analysis) 分析脚本
# 进行通路活性评分及差异通路分析
# ==============================================================================

# ========================== 加载包 ==========================

suppressPackageStartupMessages({
    library(Seurat)
    library(GSVA)
    library(GSEABase)
    library(optparse)
    library(methods)
    library(limma)
    library(genefilter)
    library(dplyr)
    library(data.table)
    library(doParallel)
    library(foreach)
    library(pheatmap)
    library(ggplot2)
    library(tibble)
    library(tidyr)
    library(stringr)
})

# ========================== 配置参数 ==========================

OUTPUT_DIR <- "16.GSVA"
DATA_FILE <- "data_ob_v3.rds"
GMT_FILE <- "utils/reference/refdata-gex-GRCm39-2024-A/gmt/mouse_2024_kegg.gmt"

# GSVA 参数
METHOD <- "gsva"
KCDF <- "Poisson"
MIN_SZ <- 2
MAX_SZ <- 100000
PARALLEL_SZ <- 6
MX_DIFF <- TRUE
TAU <- 1

# 差异分析参数
CONTRAST <- c("clusters", "all", "all")
PVAL_CUTOFF <- 0.05
TOP_N <- 10

# ========================== 初始化 ==========================

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ========================== 数据加载 ==========================

seurat_ob <- readRDS(DATA_FILE)
seurat_ob <- StashIdent(seurat_ob, save.name = "clusters")

# 读取基因集
gene_sets <- getGmt(GMT_FILE)

# ========================== GSVA 评分 ==========================

# 计算 GSVA 评分
gsva_scores <- gsva(
    as.matrix(GetAssayData(seurat_ob, slot = "counts")),
    gene_sets,
    method = METHOD,
    kcdf = KCDF,
    min.sz = MIN_SZ,
    max.sz = MAX_SZ,
    parallel.sz = PARALLEL_SZ,
    mx.diff = MX_DIFF,
    tau = TAU
)

# 保存 GSVA 评分结果
gsva_scores1 <- as.data.frame(gsva_scores)
gsva_scores1$geneset <- rownames(gsva_scores1)
gsva_scores1 <- gsva_scores1 %>% dplyr::select(geneset, everything())

fwrite(
    gsva_scores1,
    file = file.path(OUTPUT_DIR, "GSVA_enrichment_results.xls"),
    col.names = TRUE,
    row.names = FALSE,
    sep = "\t"
)

# ========================== 差异通路分析 ==========================

# 创建 GSVA Assay
seurat_ob[["GSVA"]] <- CreateAssayObject(counts = gsva_scores)
seurat_ob <- ScaleData(seurat_ob, assay = "GSVA")

# 获取分组信息
assay_metadata <- seurat_ob@meta.data
all_levels <- as.vector(unique(assay_metadata[, CONTRAST[1]]))
groupby <- CONTRAST[1]

# 生成所有对比组合
all_comparisions <- lapply(all_levels, function(x) {
    paste(CONTRAST[1], x, paste0(all_levels[-which(all_levels == x)], collapse = ","), sep = ":")
})
all_comparisions <- unlist(all_comparisions)

# 执行差异分析
gsva_results <- c()
for (contrastx in all_comparisions) {
    contrastsx <- unlist(strsplit(contrastx, ":", perl = TRUE))
    DEG_gsva_tmp <- FindMarkers(
        seurat_ob,
        ident.1 = as.character(unlist(strsplit(contrastsx[2], ",", perl = TRUE))),
        ident.2 = as.character(unlist(strsplit(contrastsx[3], ",", perl = TRUE))),
        test.use = "limma",
        min.pct = 0,
        group.by = contrastsx[1],
        logfc.threshold = -Inf,
        assay = "GSVA",
        slot = "data"
    )
    DEG_gsva_tmp$cluster <- as.character(contrastsx[2])
    DEG_gsva <- DEG_gsva_tmp %>% tibble::rownames_to_column(var = "geneset")
    gsva_results <- rbind(gsva_results, DEG_gsva)
}

# 整理结果
gsva_results <- gsva_results %>%
    dplyr::select(geneset, logFC, AveExpr, t, p_val, adj.P.Val, cluster)
colnames(gsva_results) <- c("geneset", "logFC", "avgExp", "t", "pval", "FDR", "Interested_group")

write.table(
    gsva_results,
    file.path(OUTPUT_DIR, "GSVA_results.xls"),
    quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE
)

# ========================== 筛选 Top 通路 ==========================

# 筛除 p 不显著通路，筛选各细胞类型 top N 通路
plot_term <- gsva_results %>%
    subset(pval < PVAL_CUTOFF) %>%
    group_by(Interested_group) %>%
    top_n(TOP_N, t) %>%
    arrange(Interested_group)

write.table(
    plot_term,
    file.path(OUTPUT_DIR, paste0("GSVA_top", TOP_N, "_results.xls")),
    quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE
)

# ========================== 热图可视化 ==========================

# 制作用于标注显著性星号*的表格
plot_data_p <- tidyr::spread(
    as.data.frame(plot_term[, c("geneset", "Interested_group", "t")]),
    Interested_group, t
) %>%
    tibble::column_to_rownames(var = "geneset")
plot_data_p <- plot_data_p[as.vector(unique(plot_term$geneset)), ]
display_numbers <- ifelse(is.na(plot_data_p), "", "**")

# 整体绘图数据表格
plot_data <- tidyr::spread(
    as.data.frame(gsva_results[, c("geneset", "Interested_group", "t")]),
    Interested_group, t
) %>%
    tibble::column_to_rownames(var = "geneset")
plot_data <- plot_data[as.vector(unique(plot_term$geneset)), ]

# 绘制热图
p <- pheatmap(
    plot_data,
    scale = "row",
    color = colorRampPalette(c("#406AA8", "white", "#D91216"))(100),
    show_colnames = TRUE,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    border_color = "white",
    fontsize_row = 8,
    cellheight = 12,
    cellwidth = 12,
    display_numbers = display_numbers
)

pdf(file.path(OUTPUT_DIR, "heatmap.pdf"))
print(p)
dev.off()

# ========================== 柱状图可视化 ==========================

# 两组间差异通路柱状图
DEG_gsva_tmp <- FindMarkers(
    seurat_ob,
    ident.1 = 1,
    ident.2 = 2,
    test.use = "limma",
    min.pct = 0,
    group.by = CONTRAST[1],
    logfc.threshold = -Inf,
    assay = "GSVA",
    slot = "data"
)

DEG_gsva_tmp <- DEG_gsva_tmp %>%
    tibble::rownames_to_column(var = "geneset") %>%
    dplyr::select(geneset, logFC, AveExpr, t, p_val, adj.P.Val)
colnames(DEG_gsva_tmp) <- c("geneset", "logFC", "avgExp", "t", "pval", "FDR")

write.table(
    DEG_gsva_tmp,
    file.path(OUTPUT_DIR, "diffexp_genesets_GSVA_score_c1-vs-c2.xls"),
    quote = FALSE, sep = "\t", row.names = FALSE
)

# 添加调控方向标记
DEG_gsva_tmp$just <- ifelse(DEG_gsva_tmp$t < 0, 0, 1)
DEG_gsva_tmp$Regulation <- ifelse(DEG_gsva_tmp$just == 1, "Up", "Down")
DEG_gsva_tmp$Regulation <- factor(DEG_gsva_tmp$Regulation, levels = c("Up", "Down"))

# 挑选 p < 0.05 的 top10 上调/下调通路
pathway2vis <- DEG_gsva_tmp %>%
    filter(pval < PVAL_CUTOFF) %>%
    group_by(Regulation) %>%
    arrange(abs(t)) %>%
    top_n(10, abs(t))

# 绘制柱状图
pp <- ggplot(pathway2vis, aes(x = reorder(geneset, t), y = t)) +
    geom_col(aes(fill = Regulation)) +
    scale_fill_manual(values = c("#6CC570", "#2A5078")) +
    coord_flip() +
    labs(x = "Pathway", y = "t value of GSVA") +
    theme_minimal() +
    geom_text(
        aes(x = geneset, y = 0, label = geneset),
        hjust = pathway2vis$just, size = 3.5
    ) +
    theme(
        axis.text.y = element_blank(),
        panel.grid = element_blank()
    ) +
    ylim(c(-max(abs(pathway2vis$t)) - 5, max(abs(pathway2vis$t)) + 5)) +
    labs(fill = "Group")

# 保存图片
ggsave(
    file.path(OUTPUT_DIR, "barplot.pdf"),
    plot = pp, width = 18, height = 8
)
