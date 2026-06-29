dir.create("7.InferCNV")
gtf= plyranges::read_gff('genes.gtf')
gene.chr = gtf %>% plyranges::filter(type == "gene" & gene_name %in% rownames(scRNA_mnn)) %>%
    as.data.frame() %>%
    dplyr::select(gene_name, seqnames, start, end) %>%
    dplyr::distinct(gene_name, .keep_all=T) %>% 
    dplyr::mutate( seqnames =seqnames)
count_mat = GetAssayData(scRNA_mnn, "counts")
cellanno = FetchData(scRNA_mnn, vars = 'new_celltype' ) %>% tibble::rownames_to_column(var = "cellbarcode")
tempdir = tempdir()
cnv_celltyping = file.path(tempdir, "cnv_celltype_group.xls")
write.table(cellanno, cnv_celltyping, sep = "\t", col.names = F,row.names = F, quote = F)
gene_order_f = file.path(tempdir, "gene_order_file.xls" )
write.table(gene.chr, gene_order_f, col.names =F, row.names =F, sep = "\t", quote =F )
infercnv_obj = CreateInfercnvObject(raw_counts_matrix= count_mat,
                                    annotations_file= cnv_celltyping,
                                    delim="\t",
                                    gene_order_file= gene_order_f,
                                    ref_group_names=c('T_cells'))
output_dir='7.InferCNV'
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff= 0.1, 
                             analysis_mode= 'subclusters',
                             tumor_subcluster_pval=0.05,
                             hclust_method = 'ward.D2', 
                             out_dir= output_dir,
                             num_threads=2,
                             cluster_by_groups=T,
                             denoise=T,
                             HMM=T)
pdf( file.path("7.InferCNV/heatmap.pdf"), width = 18, height = 12 )
ComplexHeatmap::Heatmap(
                     t(as.matrix(infercnv_obj@expr.data)),
                     cluster_rows = FALSE,
                     cluster_columns = FALSE,
                     show_row_names =F,
                     show_column_names = F,
                     name="CNV level",
                     use_raster=TRUE,
                     raster_quality=4 )
dev.off()
scRNA_mnn[['CNV']] = CreateAssayObject(data = infercnv_obj@expr.data)
infercnv_level <- apply(as.data.frame(t(infercnv_obj@expr.data)), 1, function(x) {
          x[is.na(x)] <- 0
          return(sum(x))
        })
infercnv_level <- round(scales::rescale(infercnv_level / nrow(infercnv_obj@expr.data), c(1, 100)), 0)
infercnv_level <- infercnv_level[Cells(scRNA_mnn)]
scRNA_mnn@meta.data$cnv_level = infercnv_level
p1 = FeaturePlot(scRNA_mnn , features = "cnv_level" )
p2 = VlnPlot(scRNA_mnn , features = "cnv_level" , group.by = "orig.ident",pt.size = 0)
ggsave("7.InferCNV/featureplot.png",p1,height = 5,width = 5)
ggsave("7.InferCNV/vlnplot.png",p2,height = 4,width = 7)