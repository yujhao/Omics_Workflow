suppressWarnings({
  library(SeuratWrappers)
  library(DoubletFinder)
  library(enrichplot)
  library(batchelor)
  library(GO.db)
  library(GSEABase)
  library(org.Hs.eg.db)
  library(DOSE) 
  library(clusterProfiler)
  library(SingleR)
  library(ggplot2)
  library(patchwork) 
  library(pheatmap)
  library(magrittr) 
  library(tidyverse)
  library(Seurat)
  library(dplyr)
  library(ggstatsplot)
  library(tidyr)
  library(Matrix)
  library(infercnv)
  library(tibble)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(dittoSeq)
  library(plyranges)
  library(monocle)
  library(ggsci)
  library(igraph)
  library(CellChat)
  library(ggalluvial)
  library(SingleCellExperiment)
})
source("./0.Outliers.r")

#1 单细胞数据质控
scRNA1 = readRDS('data_ob_v3.rds')
dir.create("1.SingleCell_QC")
#1.1 整理数据
Idents(scRNA1) <- 'orig.ident' 
scRNA1[["percent.mito"]] = PercentageFeatureSet(scRNA1,pattern = "^MT-") 
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") 
HB_m <- match(HB.genes, rownames(scRNA1@assays$RNA)) 
HB.genes <- rownames(scRNA1@assays$RNA)[HB_m]
HB.genes <- HB.genes[!is.na(HB.genes)]
scRNA1[["percent.HB"]]<-PercentageFeatureSet(scRNA1, features=HB.genes)
beforeQC_vlnplot = VlnPlot(scRNA1, 
                           features = c("nFeature_RNA", 
                                        "nCount_RNA", 
                                        "percent.mito",
                                        "percent.HB"), 
                           ncol = 4, 
                           pt.size = 0)
ggsave("1.SingleCell_QC/BeforeQC_nFeature_nCount_percent.mito_percent.HB_vlnplot.pdf", plot = beforeQC_vlnplot)
ggsave("1.SingleCell_QC/BeforeQC_nFeature_nCount_percent.mito_percent.HB_vlnplot.png", plot = beforeQC_vlnplot)

#1.2 数据过滤
filter_params = c("nFeature_RNA","nCount_RNA","percent.mito")
lower_threshold = sapply( unlist(strsplit(c("NULL","NULL"," -Inf"), ",") ),
  function(x) ifelse(x == "NULL", NA, as.numeric(x)) )
if ( length(lower_threshold) != length(filter_params) ){
  stop("The lower threshold setting is not consistant with the parameters in --filters!")}
names(lower_threshold) = filter_params
upper_threshold = sapply( unlist(strsplit(c("NULL","NULL",1), ",") ),
  function(x) ifelse(x == "NULL", NA, as.numeric(x)) )
if ( length(upper_threshold) != length(filter_params) ){
    stop("The upper threshold setting is not consistant with the parameters in --filters!")}
names(upper_threshold) = filter_params
bounds_list = list()
for ( x in filter_params)  bounds_list[[x]] = c(min = unname(lower_threshold[x]), max = unname(upper_threshold[x]) )
outliers = FindOutliers(scRNA1, vars = filter_params,
                                            var.limit = bounds_list, batch = "sampleid",
                                            type = "both", cut.1 = "mean", cut.2 = "sd",
                                            n = 2, log = FALSE )
outliercells = do.call(cbind, outliers)
metric_outlier = apply(outliercells, 1, function(x) any(x == T))
scRNA1 = AddMetaData(scRNA1, metadata = metric_outlier,
                                            col.name = "is_metric_outlier")
outlier_variables = "is_metric_outlier"
is_valid_cell = !apply(FetchData(scRNA1, vars = outlier_variables), 1, function(x) any(x == T))
scRNA1 = AddMetaData(scRNA1, metadata = is_valid_cell, col.name = "is_valid")
scRNA1 = subset(scRNA1, subset = is_valid == TRUE )
counts<-GetAssayData(scRNA1,'counts')
gmt_list = unlist(strsplit("/data/database/cellranger-refdata/refdata-gex-GRCh38-2024-A/MT_genelist.gmt",",",perl = T))
gset_list <- lapply(gmt_list, function(gmtfile){
                            gset_list <- GSEABase::geneIds(GSEABase::getGmt(con=gmtfile))  
                            return(gset_list)
        })
char_vector <- unlist(gset_list[[1]])
mt <- counts[char_vector,]
mt_umi = colSums(mt)
median_value <- summary(mt_umi)["Median"]
fust = as.data.frame(mt_umi)
fust$bak = fust$mt_umi
filter_cell = rownames(fust[(fust$mt_umi < median_value*4),])
scRNA1 = scRNA1[,filter_cell]
scRNA1@meta.data$log10GenesPerUMI <- log10(scRNA1@meta.data$nFeature_RNA)/log10(scRNA1@meta.data$nCount_RNA)
obj = SplitObject(scRNA1, split.by = "orig.ident")
obj_rm=list() 
doublets_plot = list() 
for( i in names(obj)){
    print(i)
    obj[[i]] <- NormalizeData(obj[[i]])
    obj[[i]] <- FindVariableFeatures(obj[[i]], selection.method = "vst", nfeatures = 2000)
    obj[[i]] <- ScaleData(obj[[i]])
    obj[[i]] <- RunPCA(obj[[i]])
    obj[[i]] <- RunUMAP(obj[[i]], dims = 1:30)
    obj[[i]] <- FindNeighbors(obj[[i]], dims = 1:30) %>% FindClusters(resolution = 0.3)
    tmp <- RemoveDoublets(obj[[i]], doublet.rate=0.008,pc.num=1:30)
    obj_rm[[i]] <- tmp$obj 
    doublets_plot[[i]] <- tmp$plot
  }
scRNA1 <- obj_rm[[1]] 
if(length(obj_rm) > 1) { 
  for(i in 2:length(obj_rm)) {
    scRNA1 <- merge(scRNA1, y = obj_rm[[i]]) 
  }
}
Idents(scRNA1) <- 'orig.ident' 
afterQC_vlnplot = VlnPlot(scRNA1, 
                           features = c("nFeature_RNA", 
                                        "nCount_RNA", 
                                        "percent.mito",
                                        "percent.HB"), 
                           ncol = 4, 
                           pt.size = 0) 

ggsave("1.SingleCell_QC/afterQC_nFeature_nCount_percent.mito_percent.HB_vlnplot.pdf",plot = afterQC_vlnplot,width = 8,height = 6)
ggsave("1.SingleCell_QC/afterQC_nFeature_nCount_percent.mito_percent.HB_vlnplot.png", plot = afterQC_vlnplot,width = 8,height = 6)

#1.3 数据归一化与标准化
scRNA1 <- NormalizeData(scRNA1)
scRNA1 <- FindVariableFeatures(scRNA1,  selection.method = "vst") 
scRNA1 <- ScaleData(scRNA1, features = VariableFeatures(scRNA1))
