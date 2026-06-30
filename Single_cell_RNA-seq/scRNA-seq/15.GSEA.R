#!/usr/bin/env Rscript
# ============================================================================
# GSEA (Gene Set Enrichment Analysis) 分析脚本
# 功能: 基于差异表达分析进行 GSEA 富集分析
# 输入: Seurat 对象 + 基因集 GMT 文件
# 输出: GSEA 富集结果、可视化图表
# ============================================================================

# ========================== 加载包 ==========================

suppressPackageStartupMessages({
    library(Seurat)
    library(clusterProfiler)
    library(GSEABase)
    library(ggplot2)
    library(enrichplot)
    library(DOSE)
    library(pheatmap)
    library(RColorBrewer)
})

# ========================== 配置参数 ==========================

OUTPUT_DIR <- "15.GSEA"
DATA_FILE <- "data_ob_v3.rds"
GMT_FILE <- "utils/reference/refdata-gex-GRCm39-2024-A/gmt/mouse_2024_kegg.gmt"

# 分析参数
IDENT_COLUMN <- "clusters"
GROUP1 <- "1"
GROUP2 <- "2"
PVAL_CUTOFF <- 0.05
LOGFC_CUTOFF <- 0.5
MIN_GSSIZE <- 1
PVALUE_CUTOFF <- 0.5

# 可视化参数
GSEA_TERM <- "Influenza A(mmu05164)"
PLOT_WIDTH <- 10
PLOT_HEIGHT <- 8
PLOT_DPI <- 300

# ========================== 初始化 ==========================

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# ========================== 数据加载 ==========================

# 读取 Seurat 对象
pbmc <- readRDS(DATA_FILE)
Idents(pbmc) <- IDENT_COLUMN

# 读取基因集文件
gmt <- read.gmt(GMT_FILE)
gmt$term <- gsub("KEGG_", "", gmt$term)

# ========================== 细胞亚群提取 ==========================

scedata <- pbmc
cells_sub <- subset(scedata@meta.data, clusters %in% c(GROUP1, GROUP2))
scRNA_sub <- subset(scedata, cells = row.names(cells_sub))

# ========================== 差异表达分析 ==========================

# 组间对比计算各基因 logFC
sub.markers <- FindMarkers(
    scRNA_sub,
    group.by = IDENT_COLUMN,
    ident.1 = GROUP1,
    ident.2 = GROUP2
)

# 筛选显著差异基因
deg <- subset(sub.markers, p_val_adj < PVAL_CUTOFF & avg_log2FC > LOGFC_CUTOFF)

# ========================== GSEA 富集分析 ==========================

# 提取 log2FC 值并排序
genelist <- deg$avg_log2FC
names(genelist) <- rownames(deg)
genelist <- sort(genelist, decreasing = TRUE)

# 运行 GSEA
gsea <- GSEA(
    genelist,
    TERM2GENE = gmt,
    minGSSize = MIN_GSSIZE,
    pvalueCutoff = PVALUE_CUTOFF
)

# ========================== 可视化 ==========================

# 绘制 GSEA 图
p <- gseaplot2(
    gsea,
    geneSetID = GSEA_TERM,
    subplots = c(1, 2, 3),
    color = c("#dd2222"),
    ES_geom = "line",
    pvalue_table = TRUE,
    title = FALSE
)

# 保存图片
ggsave(
    file.path(OUTPUT_DIR, "GSEA.pdf"),
    p,
    width = PLOT_WIDTH,
    height = PLOT_HEIGHT,
    dpi = PLOT_DPI
)
