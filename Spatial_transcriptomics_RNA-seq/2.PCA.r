#2 标准化与降维聚类
dir.create("2.Clustering",recursive = T)
data_ob = Seurat::SCTransform(
  data_ob,
  assay = "Spatial",
  method = "glmGamPoi", 
  vars.to.regress =  "percent.mito",
  variable.features.n = 3000,
  verbose = TRUE,
  return.only.var.genes = FALSE)
data_ob = RunPCA(data_ob, 
               assay = "SCT", 
               features = VariableFeatures(data_ob))
p1 = DimPlot(data_ob, 
             reduction = "pca",
             group.by = "orig.ident")
p2 = ElbowPlot(data_ob, ndims = 50, 
               reduction = "pca")
p3 = wrap_plots(p1,p2)
p3
ggsave("2.Clustering/RunPCA.pdf", 
       plot = p3,
       device = "pdf", 
       width = 10, 
       height = 4, 
       dpi = 300)
data_ob = FindNeighbors(data_ob, 
                      reduction = "pca", 
                      dims = 1:30, 
                      features = VariableFeatures(data_ob))
data_ob = FindClusters(data_ob, 
                     verbose = TRUE,resolution = 0.4)
data_ob = RunUMAP(data_ob, 
                reduction = "pca", 
                dims = 1:30)
data_ob@meta.data$clusters = as.numeric(data_ob@meta.data$seurat_clusters)
data_ob@meta.data$clusters = as.factor(data_ob@meta.data$clusters) 
Idents(data_ob) = "clusters"
p1 = DimPlot(data_ob, 
             reduction = "umap", 
             label = TRUE)
p2 = SpatialPlot(data_ob, 
                 ncol = 2,
                 group.by = "clusters",
                 pt.size.factor = 1.6 )
p3 = wrap_plots(p1,p2)
p3
ggsave("2.Clustering/UMAP.pdf", 
       plot = p3, 
       device = "pdf",
       width = 16, 
       height = 4, 
       dpi = 300)