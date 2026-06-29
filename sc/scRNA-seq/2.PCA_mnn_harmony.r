#2 批次矫正和降维聚类
dir.create("2.Clustering")
scRNA1 = readRDS("data/Harmony.rds")

#pca
scRNA_pca <- RunPCA(scRNA1)
scRNA_pca <- FindNeighbors(scRNA_pca, reduction = "pca", dims = 1:30)
scRNA_pca <- FindClusters(scRNA_pca,resolution = 0.4)
scRNA_pca <- RunUMAP(scRNA_pca, reduction = "pca", dims = 1:30)
scRNA_pca@meta.data$clusters = as.numeric(scRNA_pca@meta.data$seurat_clusters)
scRNA_pca@meta.data$clusters = as.factor(scRNA_pca@meta.data$clusters)
Idents(scRNA_pca) = "clusters"
p1 <- DimPlot(scRNA_pca,  pt.size=0.1,label = T) 
ggsave('2.Clustering/pca.pdf', p1, width=8, height=8)
ggsave('2.Clustering/pca.png', p1, width=8, height=8)
saveRDS(scRNA_pca,"pca.rds")


#mnn
scRNAlist <- SplitObject(scRNA1, split.by = "orig.ident")
scRNAlist <- lapply(scRNAlist, FUN = function(x) NormalizeData(x))
scRNAlist <- lapply(scRNAlist, FUN = function(x) FindVariableFeatures(x))
scRNA_mnn <- RunFastMNN(object.list = scRNAlist) 
scRNA_mnn <- FindVariableFeatures(scRNA_mnn)
scRNA_mnn <- RunUMAP(scRNA_mnn, reduction = "mnn", dims = 1:30)
scRNA_mnn <- FindNeighbors(scRNA_mnn, reduction = "mnn", dims = 1:30)
scRNA_mnn <- FindClusters(scRNA_mnn,resolution = 0.4)
scRNA_mnn@meta.data$clusters = as.numeric(scRNA_mnn@meta.data$seurat_clusters)
scRNA_mnn@meta.data$clusters = as.factor(scRNA_mnn@meta.data$clusters)
Idents(scRNA_mnn) = "clusters"
p2 <- DimPlot(scRNA_mnn,  pt.size=0.1,label = T) 
ggsave('2.Clustering/MNN.pdf', p2, width=8, height=8)
ggsave('2.Clustering/MNN.png', p2, width=8, height=8)
saveRDS(scRNA_mnn,"MNN.rds")

#Harmony
scRNA_Harmony <- RunPCA(object = scRNA1)
scRNA_Harmony <- FindNeighbors(scRNA_Harmony, reduction = "pca", dims = 1:30)
scRNA_Harmony <- FindClusters(scRNA_Harmony,resolution = 0.4)
scRNA_Harmony <- RunHarmony(scRNA_Harmony, "sampleid")
scRNA_Harmony <- RunUMAP(scRNA_Harmony,  dims = 1:30, reduction = "harmony")
scRNA_Harmony@meta.data$clusters = as.numeric(scRNA_Harmony@meta.data$seurat_clusters)
scRNA_Harmony@meta.data$clusters = as.factor(scRNA_Harmony@meta.data$clusters)
Idents(scRNA_Harmony) = "clusters"
p1 <- DimPlot(scRNA_Harmony,  pt.size=0.1,label = T) 
ggsave('2.Clustering/Harmony.pdf', p1, width=8, height=8)
ggsave('2.Clustering/Harmony.png', p1, width=8, height=8)
saveRDS(scRNA_Harmony,"Harmony.rds")

