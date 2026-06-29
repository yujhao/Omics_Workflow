dir.create("8.Monocle2")
cds <- as.CellDataSet(scRNA_mnn)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds = detectGenes(cds,min_expr = 1)
expressed_genes = row.names(subset(fData(cds), num_cells_expressed > 10))
clustering_DEGs = differentialGeneTest(cds[expressed_genes,],
                                       fullModelFormulaStr ="~clusters",cores = 4)
featureData(cds)@data[rownames(clustering_DEGs),"qval"]=clustering_DEGs$qval
ordering_genes <- row.names (subset(clustering_DEGs, qval < 0.01))
gbm_cds = setOrderingFilter(cds,ordering_genes = ordering_genes)
p = plot_ordering_genes(gbm_cds)
ggsave('8.Monocle2/ordering_genes.pdf', p, width=8, height=4)
ggsave('8.Monocle2/ordering_genes.png', p, width=8, height=4)
gbm_cds = reduceDimension(gbm_cds, max_components = 2, verbose = T,
                          check_duplicates = F, num_dim = 10)
gbm_cds = orderCells(gbm_cds, reverse = F)
p = plot_cell_trajectory(gbm_cds, color_by = "State", cell_size = 1.5, show_branch_points = T)+ scale_color_simpsons()
ggsave('8.Monocle2/cell_trajectory.pdf', p, width=8, height=4)
ggsave('8.Monocle2/cell_trajectory.png', p, width=8, height=4)
p = plot_cell_trajectory(gbm_cds, color_by = "clusters", 
                     cell_size = 1.5, show_branch_points = T)
ggsave('8.Monocle2/cell_trajectory.pdf', p, width=8, height=4)
ggsave('8.Monocle2/cell_trajectory.png', p, width=8, height=4)
p = plot_cell_trajectory(gbm_cds, color_by = "clusters", 
                     cell_size = 1.5, show_branch_points = T)+
     facet_wrap(~clusters)
ggsave('8.Monocle2/cell_trajectory_facet.pdf', p, width=8, height=4)
ggsave('8.Monocle2/cell_trajectory_facet.png', p, width=8, height=4)

p=plot_cell_trajectory(gbm_cds, color_by = "Pseudotime",
                     show_branch_points = F) +
         scale_colour_viridis_c(option = "inferno")
ggsave('8.Monocle2/cell_trajectory_pseudotime.pdf', p, width=8, height=4)
ggsave('8.Monocle2/cell_trajectory_pseudotime.png', p, width=8, height=4)
p = plot_complex_cell_trajectory(gbm_cds, color_by = "clusters", 
                                  show_branch_points = T, cell_size = 1, 
                                  cell_link_size = 0.3) 
ggsave('8.Monocle2/cell_trajectory_facet.pdf', p, width=8, height=4)
ggsave('8.Monocle2/cell_trajectory_facet.png', p, width=8, height=4)
p = plot_complex_cell_trajectory(gbm_cds, color_by = "clusters", 
                                  show_branch_points = T, 
                                  cell_size = 1, 
                                  cell_link_size = 0.3) + facet_wrap(~clusters)
ggsave('8.Monocle2/cell_trajectory_facet.pdf', p, width=8, height=4)
ggsave('8.Monocle2/cell_trajectory_facet.png', p, width=8, height=4)
genes <- as.factor(subset(gbm_cds@featureData@data, use_for_ordering == TRUE)$gene_short_name)
to_be_tested <- row.names(subset(fData(gbm_cds), gene_short_name %in% levels(genes)))
gbm_cds <- gbm_cds[to_be_tested, ]
varMetadata(gbm_cds)[,1] = rownames(varMetadata(gbm_cds))
gbm_cds@featureData@varMetadata[,1] = rownames(gbm_cds@featureData@varMetadata)
p <- plot_pseudotime_heatmap(gbm_cds, cores = 1, 
                             cluster_rows = T, num_clusters = 4, 
                             show_rownames = F, return_heatmap = T)
ggsave('8.Monocle2/pseudotime_heatmap.pdf', p$ph_res, width=8, height=4)
ggsave('8.Monocle2/pseudotime_heatmap.png', p$ph_res, width=8, height=4)
gene_clusters <- cutree(p$tree_row, k = 4)
gene_clustering <- data.frame(gene_clusters)
gene_clustering[, 1] <- as.character(gene_clustering[, 1])
colnames(gene_clustering) <- "gene_module"
BEAM_res <- BEAM(gbm_cds, branch_point = 2, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
p = plot_genes_branched_heatmap(gbm_cds[row.names(subset(BEAM_res,
                                                  qval < 1e-4)),],
                            branch_point = 2,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = F,
                            return_heatmap = T)
ggsave("8.Monocle2/pseudotime_heatmap_branchtime.pdf", plot = p$ph_res,
      height =7,width = 7, bg="white")
ggsave("pseudotime_heatmap_branchtime.png", plot = p$ph_res,
      height =7,width = 7, dpi = 1000,bg="white")
p$heatmap_matrix_ori[1:5,1:10]
gene_clusters <- cutree(p$ph_res$tree_row, k = 4)
gene_clustering <- data.frame(gene_clusters)
gene_clustering[, 1] <- as.character(gene_clustering[, 1])
colnames(gene_clustering) <- "gene_module"
to_be_tested_sub <- row.names(subset(fData(gbm_cds), 
                                     gene_short_name %in% c("ISG15","PARK7", "GLUL")))
p = plot_genes_jitter(gbm_cds[to_be_tested_sub,], grouping = "State", 
                  min_expr = 0.1,color_by = "State",cell_size = 1)+ 
           scale_color_simpsons()
ggsave('8.Monocle2/genes_jitter.pdf', p, width=8, height=4)
ggsave('8.Monocle2/genes_jitter.png', p, width=8, height=4)
p = plot_genes_in_pseudotime(gbm_cds[to_be_tested_sub,],
                         color_by = "clusters", cell_size = 1, ncol = 1) +scale_color_simpsons()
ggsave('8.Monocle2/genes_in_pseudotime.pdf', p, width=8, height=4)
ggsave('8.Monocle2/genes_in_pseudotime.png', p, width=8, height=4)
branchpoint = 2
new_cds <- buildBranchCellDataSet(gbm_cds[to_be_tested_sub,],
                                  branch_point = branchpoint, 
                                  progenitor_method = "duplicate")
cell_fate1 <- unique(pData(new_cds)[which(pData(new_cds)$Branch == unique(pData(new_cds)$Branch)[1]), ]$State)
cell_fate2 <- unique(pData(new_cds)[which(pData(new_cds)$Branch == unique(pData(new_cds)$Branch)[2]), ]$State)
branch_labels <- c(paste("State", paste(sort(setdiff(cell_fate1, cell_fate2)), collapse = "-")), 
                   paste("State", paste(sort(setdiff(cell_fate2, cell_fate1)), collapse = "-")))
p = plot_genes_branched_pseudotime(gbm_cds[to_be_tested_sub,], 
                               color_by = "clusters", 
                               branch_point = branchpoint, 
                               cell_size = 1, ncol = 1, 
                               branch_labels = branch_labels) + scale_color_simpsons()
ggsave('8.Monocle2/genes_branched_pseudotime.pdf', p, width=8, height=4)
ggsave('8.Monocle2/genes_branched_pseudotime.png', p, width=8, height=4)
p <- plot_cell_trajectory(gbm_cds, markers = "GLUL",
                          use_color_gradient = T, show_branch_points = F, 
                          show_tree = F, cell_size = 1.5) + 
  theme(legend.text = element_text(size = 10)) + 
  scale_color_gradientn(colours = c("grey", "yellow", "red"))
ggsave('8.Monocle2/genes_branched_pseudotime.pdf', p, width=8, height=4)
ggsave('8.Monocle2/genes_branched_pseudotime.png', p, width=8, height=4)