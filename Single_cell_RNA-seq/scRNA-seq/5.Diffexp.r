#5 差异基因鉴定
dir.create("5.Diffexp")
contrasts = unlist(strsplit(c("group:case:control"), ":", perl = T))
Diff_exp = FindMarkers(scRNA_mnn, 
                               logfc.threshold = 0.25, 
                               only.pos = TRUE,
                               ident.1 = contrasts[2],
                               ident.2 = contrasts[3], 
                               group.by = contrasts[1],
                               test.use = "wilcox")
Diff_exp = Diff_exp %>% tibble::rownames_to_column(var = "gene") %>%
        dplyr::rename(pvalue = p_val, padj = p_val_adj)
Diff_exp1 = Diff_exp %>% dplyr::mutate(FoldChange = 2^avg_log2FC) %>%
        dplyr::rename(log2FoldChange = avg_log2FC) %>% dplyr::select(gene,
        everything())
colnames(Diff_exp1) =c("gene","p-value","log2FoldChange","pct.1","pct.2","q-value","FoldChange")
res_Significant = dplyr::filter(Diff_exp1, `p-value` < 0.05,
            abs(log2FoldChange) > log2(1.5))
res_Significant[which(res_Significant$log2FoldChange > 0),
        "Regulation"] <- "Up"
res_Significant[which(res_Significant$log2FoldChange < 0),
        "Regulation"] <- "Down"
colnames(res_Significant) =c("gene","p-value","log2FoldChange","pct.1","pct.2","q-value","FoldChange","Regulation")
write.table(Diff_exp1, 
            "5.Diffexp/all_diff.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t",quote = F)
write.table(res_Significant, 
            "5.Diffexp/diff_p<0.05_FC>1.5.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t",quote = F)

# 富集分析
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
GO_dot = dotplot(ego_all,split = "ONTOLOGY") + facet_grid(ONTOLOGY~.,scales = "free") 
GO_bar = barplot(ego_all,split = "ONTOLOGY")+ facet_grid(ONTOLOGY~.,scales = "free")
res_plot <- CombinePlots(list(GO_dot,GO_bar), nrow=1)
ggsave("5.Diffexp/GO_results_all.pdf", plot=res_plot, width = 12,height = 10)
ggsave("5.Diffexp/GO_results_all.png", plot=res_plot, width = 12,height = 10)

#6.2 KEGG富集分析
# https://davidbioinformatics.nih.gov/conversion.jsp
# http://bioinfo.org/kobas/genelist/
