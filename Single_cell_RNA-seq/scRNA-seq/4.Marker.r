#4 marker基因鉴定
dir.create("4.Marker")
all.markers = FindAllMarkers(scRNA_mnn, only.pos = TRUE)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(all.markers, 
            "4.Marker/all_Markers_of_each_clusters.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t",quote = F)
write.table(top10, 
            "4.Marker/top10_Markers_of_each_clusters.xls", 
            col.names = T, 
            row.names = F, 
            sep = "\t",quote = F)
scRNA_mnn <- ScaleData(scRNA_mnn, features = row.names(scRNA_mnn))
heatmap_plot = DoHeatmap(object = scRNA_mnn, 
                         features = as.vector(top10$gene), 
                         group.by = "clusters", 
                         group.bar = T, 
                         size = 3) +
  theme(axis.text.y = element_text(size = 4))
ggsave("4.Marker/top10_marker_of_each_cluster_heatmap.pdf", width = 12, height = 12,
       plot = heatmap_plot)
ggsave("4.Marker/top10_marker_of_each_cluster_heatmap.png", width = 12, height = 12,
       plot = heatmap_plot)
