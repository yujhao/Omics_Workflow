suppressWarnings({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(data.table)
  library(clusterProfiler)
  library(fgsea)
  library(enrichplot)
  library(pheatmap)
  library(igraph)
  library(ggraph)
  library(SPOTlight)
  library(igraph)
  library(RColorBrewer)
  library(spacexr)
  library(doParallel)
  library(glmGamPoi)
  library(tidyverse)
})

#1 10×空间转录组数据读取 
dir.create("1.SpaceRanger_QC")
data_ob = readRDS("data/data_ob.rds")
data_ob = PercentageFeatureSet(data_ob, pattern = "^MT-",
                             col.name = "percent.mito")
Idents(data_ob) = "sampleid"
for (i in c("nCount_Spatial","nFeature_Spatial","percent.mito")){
p1 = VlnPlot(data_ob,
             features = i,
             pt.size = 0.1) + NoLegend()
p2 = SpatialFeaturePlot(data_ob,
                 ncol = 2,
                 features=i,
                 pt.size.factor = 1.4,
                 alpha = 1)
plot = wrap_plots(p1,p2)
plot = plot
ggsave(paste0("1.SpaceRanger_QC/",i,".pdf"),
       plot = plot,
       device = "pdf",
       width = 16, 
       dpi = 300)
ggsave(paste0("1.SpaceRanger_QC/",i,".png"),
       plot = plot,
       device = "png",
       width = 16, 
       dpi = 300)
}