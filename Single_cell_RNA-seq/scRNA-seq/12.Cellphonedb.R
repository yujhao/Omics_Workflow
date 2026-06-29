
suppressWarnings({
    suppressPackageStartupMessages( library("annotables") )
    suppressPackageStartupMessages( library("future") )
    suppressPackageStartupMessages( library("future.apply") )
    suppressPackageStartupMessages( library("homologene") )
    suppressPackageStartupMessages( library("dplyr") )
    suppressPackageStartupMessages( library("optparse") )
    suppressPackageStartupMessages( library("Seurat") )
    suppressPackageStartupMessages( library("dplyr") )
    suppressPackageStartupMessages( library("optparse") )
    suppressPackageStartupMessages( library("ggplot2") )
    suppressPackageStartupMessages( library("tidyr") )
    suppressPackageStartupMessages( library("network") )
    suppressPackageStartupMessages( library("igraph") )
    suppressPackageStartupMessages( library("circlize") )
    suppressPackageStartupMessages( library("ComplexHeatmap") )
    suppressPackageStartupMessages( library("circlize") ) 
    suppressPackageStartupMessages( library("ComplexHeatmap") )
    suppressPackageStartupMessages( library("RColorBrewer") )
    source("/hwstorage/oe-scrna/jhyu/Git_lab/sh/yujunhao/cellphonedb_analysis/test/vis_cellcomm_source.r")
})


output_dir = 'result'
dir.create(output_dir,recursive=T)
groupby = NULL
column4cell = 'clusters'
threads = 10
iterations = 1000
threshold = 0.1
pvalue = 0.05
genetype = "gene_name"
database = 'cellphonedb.zip'
is_onlySig = TRUE
species = "mouse"
microenvs = ''

# setting the cores for parallization
options(future.globals.maxSize= Inf ) # setting the maxumium mermory usage much bigger in case of big data
plan("multicore", workers = min(availableCores(), threads )) # parallization using specified CPUs start from here

seurat_ob = readRDS('data_ob_v3.rds')


sp <- list()
sp[['all']] <- seurat_ob
if(class(seurat_ob@meta.data[,column4cell]) == "factor"){sp = lapply(sp, function(x){ x@meta.data[,column4cell] = droplevels(x@meta.data[,column4cell] )
                                                                                        return(x)}) }

#get the count matrix
counts = lapply(sp, function(x){
            counts <- GetAssayData(x, slot = 'data') 
            return(counts)
            })
metadata = lapply(sp, function(x){
            if ( is.null(column4cell) ){
                metadata <- Seurat::FetchData(x, vars = "ident") %>%
                            tibble::rownames_to_column(var = "Cell") %>%
                            dplyr::rename("cell_type" = ident )
            }else{
                metadata <- Seurat::FetchData(x, vars = column4cell ) %>%
                            tibble::rownames_to_column(var = "Cell")
                    metadata[,2] = gsub("[^[:alnum:]]", "_", metadata[,2])
                    metadata[,2] = gsub("__+", "_", metadata[,2])
                colnames(metadata) = c("Cell","cell_type")
                if( column4cell == "clusters" ){ metadata$cell_type = paste0("C",metadata$cell_type) }
            }
            return(metadata)
            })
rm(seurat_ob, sp)
gc()

if ( species == "mouse" ){
  mouse2human <- homologene::human2mouse(annotables::grch38$symbol) %>% dplyr::select( mouseGene, humanGene )
  print("物种为小鼠，counts基因名将替换为人的同源基因.")
  write.table(mouse2human, file = file.path(output_dir,'mouse2human.homologene.xls'), 
                                    quote = F, col.names = T, row.names = F, sep = '\t')

  counts <- lapply(counts, function(x){
                    genes2use <- intersect(names(Seurat::CaseMatch(rownames(x), mouse2human$mouseGene)), rownames(x) )
                    counts <- x[genes2use,]
                    rownames(counts) = plyr::mapvalues(rownames(counts),
                                                        from = mouse2human$mouseGene,
                                                        to = mouse2human$humanGene, warn_missing = F)
                    return(counts)
            })
}
cmds = list()
i = "all"
    # tempdir <- tempdir()
if ( is.null(groupby) ) {
    print("开始做不分组cellphoneDB分析.")
    tempdir <- output_dir
} else {
    print(paste0("开始做分组" ,i, "的cellphoneDB分析."))
    tempdir <- file.path(output_dir, i)
    dir.create(tempdir)
}
print("以下细胞群将进入细胞通讯分析：")
print(table(metadata[[i]]$cell_type))
print(paste0(c("基因数是：","细胞数是："),dim(counts[[i]])))
# using vroom to speed up the count matrix writing
counts_out <- tibble::rownames_to_column( as.data.frame(counts[[i]]), var = 'Gene')
vroom::vroom_write(x = counts_out, path = file.path(tempdir,'counts.tsv'), col_names = T, quote = "none", delim = "\t")
write.table(metadata[[i]], file = file.path(tempdir,'metadata.tsv'), 
            quote = F, col.names = T, row.names = F, sep = '\t')
outdir <- file.path(tempdir, "out")

docellphone <- "docellphone.py"
cmds[[i]] <- glue::glue('echo "group {i} begin";cd {tempdir} && module purge &&  cpdb && python {docellphone} --counts={tempdir}/counts.tsv ',
                            '--metadata={tempdir}/metadata.tsv --genetype={genetype} --iterations={iterations}',
                            ' --threshold={threshold} --threads={threads} --pvalue={pvalue} ',
                            '--database {database} --microenvs={microenvs}'  )             

print(glue::glue({cmds[[i]]}))
system(cmds[[i]])
print(paste0("group ",i," cellphoneDB execution finished."))
print(paste0("结果路径：",tempdir))
outdir <- file.path(tempdir, "out")
fileNames <- list.files(path = outdir, full.names = FALSE)
pattern <- "^statistical_analysis_(deconvoluted|means|pvalues|significant_means)_.*\\.txt$"
for (file in fileNames) {
    if (grepl(pattern, file)) {
    new_file <- sub(pattern, "\\1.txt", file)
    new_file <- file.path(outdir, new_file)
    old_file <- file.path(outdir, file)
    file.rename(old_file, new_file)
    cat("Renamed file:", old_file, "to", new_file, "\n")
    }
}
for (file in c('means.txt','pvalues.txt','significant_means.txt')){
    temp_data = read.delim(file.path(outdir, file), header = T, sep = "\t")
    temp_data <- temp_data %>% dplyr::select(-c('directionality', 'classification'))
    colnames(temp_data) = gsub('\\.','|',colnames(temp_data))
    write.table(temp_data, file.path(outdir, file), sep = "\t", col.names = T, row.names = F
        , quote = F)
}
cellphonedb <- ParseCpdb( outdir , pvalue = pvalue )
write.table(cellphonedb$ligrec, file.path(tempdir, "cell_comm_annotation.xls"), sep = "\t", col.names = T, row.names = F, quote = F)
saveRDS(cellphonedb, file.path(tempdir, "cellphonedb_results.rds"))

# plot

    
data = read.table('cell_comm_annotation.xls', header = T, sep = "\t", stringsAsFactors = F )
if(is.numeric(sort(unique(data$ligand_cell))) && is.numeric(sort(unique(data$receptor_cell)))){
        data$ligand_cell = factor(data$ligand_cell,levels = sort(unique(data$ligand_cell)))
        data$receptor_cell = factor(data$receptor_cell,levels = sort(unique(data$receptor_cell)))
}

data$group = "all"

net = data %>% filter(pval < 0.05) %>% group_by(ligand_cell, receptor_cell) %>%
        dplyr::summarize(n=n()) %>% rename(significant_pairs =  n )
write.table(net,file.path(output_dir, "cell_comm_summary.xls"),quote=F,sep="\t",row.names=F)


# "dotplot" 
out = LRDotplot(data = data, is_onlySig = is_onlySig, topn = 5, xsize = 15, xangle = 45,remove.isolate = TRUE)
ggdot <- out[[1]]
plot_data <- out[[2]]
png_width = length(unique(plot_data$clusters))
png_height = length(unique(plot_data$pair))
max_length_width = max(nchar(as.character(plot_data$clusters)))
max_length_height = max(nchar(as.character(plot_data$pair)))
if (max_length_width>40){
    width=png_width*0.45+max_length_width*1.5+8
}else{
    width=png_width*0.45+6
}
height=png_height*0.4+3
ggsave(file.path(output_dir, "cell_comm_dotplot.pdf"), limitsize = F, 
            plot = ggdot, width = width, height = height,bg="white")
ggsave(file.path(output_dir, "cell_comm_dotplot.png"), limitsize = F, 
        plot = ggdot, width = png_width*0.4+6, height = png_height*0.3)


#"network" 
datax = data
celltypes = unique(c(as.character(datax$receptor_cell), as.character(datax$ligand_cell)))
colx = c("#7fc97f","#beaed4","#fdc086","#386cb0")
names(colx) = celltypes

pdf(file.path(output_dir, "cell_comm_network.pdf"),width = 6.5+max(nchar(unique(as.character(c(data$ligand_cell,data$receptor_cell )))))*0.4, height = 7)
LRNetwork( data %>% filter(pval < 0.05),
                    col = colx, edge.label.cex = 0.3, 
                    edge.max.width = 5, arrow.width = 0.5, 
                    vertex.label.cex = 0.5, vertex.size = 10)
dev.off()
png(file.path(output_dir, "cell_comm_network.png"), width = 6.5+max(nchar(unique(as.character(c(data$ligand_cell,data$receptor_cell )))))*0.4, height = 7, res = 500, units = "in" )
LRNetwork( data %>% filter(pval < 0.05),
                    col = colx, edge.label.cex = 0.3, 
                    edge.max.width = 5, arrow.width = 0.5, 
                    vertex.label.cex = 0.5, vertex.size = 10)
dev.off()


#plotx == "circos"

pdf( file.path(output_dir, "cell_comm_circos_plot.pdf"), width = 20, height = 10)
LRCircos( data, gap.degree = 0.05, cell_col = colx, screenvar = "expr",
                topn = 5, labels.cex = 0.35,link.lwd=1, arr.length = 0.1)
dev.off()
png( file.path(output_dir, "cell_comm_circos_plot.png"), width = 20, height = 10, res = 500, units = "in" )
LRCircos( data, gap.degree = 0.05, cell_col = colx, screenvar = "expr",
                topn = 5, labels.cex = 0.35,link.lwd=1, arr.length = 0.1)
dev.off()


# "chorddiagram" 

pdf( file.path(output_dir, "cell_comm_chorddiagram_plot.pdf"), width = 9, height = 7)
par(cex = 0.6, font = 2)
LRChorddiagram(data, grid_orbit_h=0.02,label_h=0.04,diffHeight=4,output_dir=output_dir, colx = colx, extra_legend=TRUE)
dev.off()
png( file.path(output_dir, "cell_comm_chorddiagram_plot.png"), width = 9, height = 7, res = 96, units = "in" )
par(cex = 0.6, font = 2)
LRChorddiagram(data, grid_orbit_h=0.02,label_h=0.04,diffHeight=4,output_dir=output_dir, colx = colx, extra_legend=TRUE)
dev.off()
circos.clear()


#"heatmap" 

pdf( file.path(output_dir, "cell_comm_heatmap_plot.pdf") )
print(LRHeatmap(data))
dev.off()
png( file.path(output_dir, "cell_comm_heatmap_plot.png"), width = 7, height = 7, res = 96, units = "in" )
print(LRHeatmap(data))
dev.off()


#"barplot"
out = LRBarplot(data = data, bar_width = 0.6,colx=colx) + scale_y_continuous(expand=c(0,0))
ggsave(file.path(output_dir, "cell_comm_histogram_plot.pdf"), plot = out, width = ifelse( length(celltypes)<7,7,length(celltypes)),bg="white")
ggsave(file.path(output_dir, "cell_comm_histogram_plot.png"), plot = out, width = ifelse( length(celltypes)<7,7,length(celltypes)), dpi=1000,bg="white")


