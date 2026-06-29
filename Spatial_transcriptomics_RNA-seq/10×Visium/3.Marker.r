#3 marker基因鉴定
dir.create("3.Marker")
All_markers = FindAllMarkers(data_ob, 
                             only.pos = TRUE,
                             test.use = "wilcox")
top10 = All_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(All_markers, 
            "3.Marker/all_Markers_of_each_clusters.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t",quote = F)
write.table(top10, 
            "3.Marker/top10_Markers_of_each_clusters.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t",quote = F)
p1 = DoHeatmap(data_ob, 
               features = top10$gene,
               group.by="clusters",
               size = 2)
ggsave("3.Marker/Heatmap.pdf", 
       plot = p1,
       device = "pdf", 
       width = 10, 
       height = 15, 
       dpi = 300)