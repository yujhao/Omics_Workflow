#4 差异基因鉴定
dir.create("4.Diffexp")
contrasts = unlist(strsplit(c("group:H_cj1:H_cj2"), ":", perl = T))
numerator_subset = unlist(strsplit(contrasts[2], ",", perl = T))
denominator_subset = unlist(strsplit(contrasts[3], ",", perl = T))
cellmeta = Seurat::FetchData(data_ob, vars = contrasts[1]) %>% tibble::rownames_to_column(var = "barcode")
numerator = cellmeta %>% dplyr::filter(!!rlang::sym(contrasts[1]) %in%
    numerator_subset) %>% dplyr::pull(barcode)
denominator = cellmeta %>% dplyr::filter(!!rlang::sym(contrasts[1]) %in%
    denominator_subset) %>% dplyr::pull(barcode)
Diff_exp = FindMarkers(data_ob, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25, 
                               only.pos = TRUE,
                               ident.1 = numerator_subset,
                               ident.2 = denominator_subset, 
                               group.by = contrasts[1],
                               test.use = "wilcox")
Diff_exp = Diff_exp[abs(Diff_exp$avg_log2FC) > 0.25,]
Diff_exp = Diff_exp %>% tibble::rownames_to_column(var = "gene") %>%
        dplyr::rename(pvalue = p_val, padj = p_val_adj)
numerator_means = Matrix::rowMeans(SeuratObject::GetAssayData(data_ob,
        slot = "data")[Diff_exp$gene, numerator])
denominator_means = Matrix::rowMeans(SeuratObject::GetAssayData(data_ob,
        slot = "data")[Diff_exp$gene, denominator])
Diff_exp1 = Diff_exp %>% dplyr::mutate(FoldChange = 2^avg_log2FC, baseMean = 1/2 *
        (log2(numerator_means) + log2(denominator_means))) %>%
        dplyr::rename(log2FoldChange = avg_log2FC) %>% dplyr::select(gene,
        everything()) %>% dplyr::select(-baseMean)
colnames(Diff_exp1) =c("gene","p-value","log2FoldChange","pct.1","pct.2","q-value","FoldChange")
res_Significant = dplyr::filter(Diff_exp1, `p-value` < 0.05,
            abs(log2FoldChange) > log2(1.5))
res_Significant[which(res_Significant$log2FoldChange > 0),
        "Regulation"] <- "Up"
res_Significant[which(res_Significant$log2FoldChange < 0),
        "Regulation"] <- "Down"
colnames(res_Significant) =c("gene","p-value","log2FoldChange","pct.1","pct.2","q-value","FoldChange","Regulation")
write.table(Diff_exp, 
            "4.Diffexp/all_diff.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")
write.table(res_Significant, 
            "4.Diffexp/diff_p<0.05_FC>1.5.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t")


#4.1 GO富集分析
genes_symbol <- as.character(res_Significant$gene)
eg = bitr(genes_symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
id = as.character(eg[,2])
ego <- enrichGO(gene = id,
                OrgDb = org.Hs.eg.db,
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
GO_dot = dotplot(ego,split = "ONTOLOGY") + facet_grid(ONTOLOGY~.,scales = "free") 
GO_bar = barplot(ego,split = "ONTOLOGY")+ facet_grid(ONTOLOGY~.,scales = "free")
res_plot <- CombinePlots(list(GO_dot,GO_bar), nrow=1)
ggsave("4.Diffexp/GO_results_all.pdf", plot=res_plot, width = 12,height = 10)
ggsave("4.Diffexp/GO_results_all.png", plot=res_plot, width = 12,height = 10)

#4.2 KEGG富集分析
# https://davidbioinformatics.nih.gov/conversion.jsp
# http://bioinfo.org/kobas/genelist/

# 4.3 火山图
library(ggrepel)
df <- read.table("group_Dilated_Zone-vs-Proximal_Zone-volcano-p-val-0.05-FC-1.2.gene_symbol_gene.xls", sep = "\t", header = TRUE)
rownames(df) <- make.unique(as.character(df[[1]]))

df$group <- "NS"
df$group[df$p.value < 0.05 & df$log2FoldChange > 0] <- "Up"
df$group[df$p.value < 0.05 & df$log2FoldChange < 0] <- "Down"
df$group[df$p.value < 0.05 & df$log2FoldChange > log2(1.2)] <- "SigUp"
df$group[df$p.value < 0.05 & df$log2FoldChange < -log2(1.2)] <- "SigDown"
df$group <- factor(df$group, levels = c("SigUp", "SigDown", "Up", "Down", "NS"))

df_up <- df[df$p.value < 0.05 & df$log2FoldChange > log2(1.2), ]
df_down <- df[df$p.value < 0.05 & df$log2FoldChange < -log2(1.2), ]

df_up <- head(df_up[order(df_up$log2FoldChange, decreasing = TRUE), ], 20)
df_down <- head(df_down[order(df_down$log2FoldChange, decreasing = FALSE), ], 20)

df_label <- rbind(df_up, df_down)
df_label$label <- rownames(df_label)

p <- ggplot(df, aes(log2FoldChange, -log10(p.value), color = group)) +
        geom_point() +
        geom_vline(xintercept = c(-log2(1.2), log2(1.2))) +
        geom_hline(yintercept = -log10(0.05), linetype = 2) +
        geom_text_repel(
        data = df_label,
        aes(label = label)
        ) +
        theme_bw()

ggsave("volcano.png", p, width = 10, height = 8, dpi = 300)
ggsave("volcano.pdf", p, width = 10, height = 8)

write.csv(df,"volcano.csv")
  
