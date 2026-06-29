#3 相关性分析
dir.create("3.Correlation")
groupby = "clusters"
groupby_data=AverageExpression(scRNA_mnn, group.by = groupby)[["RNA"]]
data = tibble::rownames_to_column(as.data.frame(groupby_data),var="GeneID")
write.table(data, file.path("3.Correlation",paste0("normalized_data_groupby_",groupby,".xls")),quote = F, row.names = F, sep = "\t")
colnames(groupby_data) = gsub('^',paste0("clusters","_"),colnames(groupby_data))
matrix<-cor(groupby_data,method="pearson")
wid<-5+1.5*log2(length(colnames(data)))
hig<-5+1.5*log2(length(colnames(data)))
coefficient = pheatmap::pheatmap(matrix,
                      display_numbers = F,
                      border_color = "white",
                      scale = "none",
                      fontsize_number=(10.0+0.0001*log2(length(colnames(data)))),
                      number_format = "%.1f",
                      fontsize_row = (10.0+0.0001*log2(length(colnames(data)))),
                      fontsize_col = (10.0+0.0001*log2(length(colnames(data)))),
                      number_color="black",
                      angle_col=45)
ggsave('3.Correlation/coefficient_heatmap.pdf', coefficient, width=8, height=8)
ggsave('3.Correlation/coefficient_heatmap.png', coefficient, width=8, height=8)
