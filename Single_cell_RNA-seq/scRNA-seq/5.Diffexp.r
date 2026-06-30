#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ==============================================================================
# 差异基因鉴定与富集分析脚本
# 进行两组间的差异表达分析及 GO/KEGG 富集分析
# ==============================================================================

library(clusterProfiler)
library(org.Hs.eg.db)

# ========================== 配置参数 ==========================

OUTPUT_DIR <- "5.Diffexp"
CONTRASTS <- c("group", "case", "control")  # 分组变量, 实验组, 对照组
PVALUE_CUTOFF <- 0.05
LOGFC_CUTOFF <- log2(1.5)

# ========================== 初始化 ==========================

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# ========================== 差异表达分析 ==========================

# 执行差异表达分析（Wilcoxon 检验）
Diff_exp <- FindMarkers(
    scRNA_mnn,
    logfc.threshold = 0.25,
    only.pos = FALSE,
    ident.1 = CONTRASTS[2],
    ident.2 = CONTRASTS[3],
    group.by = CONTRASTS[1],
    test.use = "wilcox"
)

# 整理结果格式
Diff_exp <- Diff_exp %>%
    tibble::rownames_to_column(var = "gene") %>%
    dplyr::rename(pvalue = p_val, padj = p_val_adj)

Diff_exp1 <- Diff_exp %>%
    dplyr::mutate(FoldChange = 2^avg_log2FC) %>%
    dplyr::rename(log2FoldChange = avg_log2FC) %>%
    dplyr::select(gene, everything())

colnames(Diff_exp1) <- c("gene", "p-value", "log2FoldChange", "pct.1", "pct.2", "q-value", "FoldChange")

# 筛选显著差异基因
res_Significant <- dplyr::filter(Diff_exp1, `p-value` < PVALUE_CUTOFF, abs(log2FoldChange) > LOGFC_CUTOFF)
res_Significant[which(res_Significant$log2FoldChange > 0), "Regulation"] <- "Up"
res_Significant[which(res_Significant$log2FoldChange < 0), "Regulation"] <- "Down"
colnames(res_Significant) <- c("gene", "p-value", "log2FoldChange", "pct.1", "pct.2", "q-value", "FoldChange", "Regulation")

# ========================== 保存结果 ==========================

write.table(
    Diff_exp1,
    file.path(OUTPUT_DIR, "all_diff.xls"),
    col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE
)

write.table(
    res_Significant,
    file.path(OUTPUT_DIR, "diff_p<0.05_FC>1.5.xls"),
    col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE
)

# ========================== GO 富集分析 ==========================

# 基因符号转换为 Entrez ID
genes_symbol <- as.character(res_Significant$gene)
eg <- bitr(genes_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
id <- as.character(eg[, 2])

# GO 富集分析
ego <- enrichGO(
    gene = id,
    OrgDb = org.Hs.eg.db,
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
)

# 绘制 GO 结果图
GO_dot <- dotplot(ego, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scales = "free")
GO_bar <- barplot(ego, split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scales = "free")
res_plot <- CombinePlots(list(GO_dot, GO_bar), nrow = 1)

ggsave(file.path(OUTPUT_DIR, "GO_results_all.pdf"), plot = res_plot, width = 12, height = 10)
ggsave(file.path(OUTPUT_DIR, "GO_results_all.png"), plot = res_plot, width = 12, height = 10)

# ========================== KEGG 富集分析 ==========================

# 注：如需进行 KEGG 分析，请取消以下注释
# ekegg <- enrichKEGG(
#     gene = id,
#     organism = "hsa",
#     pvalueCutoff = 0.05,
#     qvalueCutoff = 0.05
# )
# ggsave(file.path(OUTPUT_DIR, "KEGG_results.pdf"), plot = dotplot(ekegg), width = 10, height = 8)
