

library(schard)
library(ggplot2)
seu_obj = schard::h5ad2seurat('vi/2.Clustering/clustering_result.h5ad')

head(seu_obj)

p = DimPlot(seu_obj,group.by = "clusters")
ggsave('vi/2.Clustering/clustering_result.png',p)
saveRDS(seu_obj,"vi/2.Clustering/clustering_result.rds")
