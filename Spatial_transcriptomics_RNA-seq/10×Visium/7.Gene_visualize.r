library(RColorBrewer)
library(Seurat)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(ggsci)
library(RColorBrewer)

# 1.dotplot
genes = c("DES","FLNA","PRPH","ERN2")
p = DotPlot(data_ob,features = genes,
    group.by = "clusters",
    cols = colorRampPalette(RColorBrewer::brewer.pal(11,"Spectral"))(100))
ggsave("dotplot.pdf",p,device = "pdf",width = 10,height = 6)
ggsave("dotplot.png",p,device = "png",width = 10,height = 6,dpi = 300)

# 2 .featureplot
features = c("DES","FLNA","PRPH","ERN2")
p2 = SpatialFeaturePlot(data_ob,
                 ncol = 2,
                 features=features,
                 pt.size.factor = 1.4,
                 alpha = 1)

ggsave(paste0("gene_spatial.pdf"),
       plot = p2,
       device = "pdf",
       width = 16, 
       dpi = 300)
ggsave(paste0("gene_spatial.png"),
       plot = p2,
       device = "png",
       width = 16, 
       dpi = 300)

# 3 vlnplot
features = c("FN1","COL4A1","COL1A2","COL6A2","THBS1")
vln_df = data_ob@assays$SCT@data[features, ]
vln_df = as.data.frame(t(vln_df))
vln_df$clusters = data_ob$clusters
vln_df = reshape2::melt(vln_df,
                        id.vars = "clusters",
                        variable.name = "gene",
                        value.name = "expression")

vln_df$gene_padded = vln_df$gene
p3 = ggplot(vln_df, aes_string(x = "clusters", y = "expression")) +
        geom_violin(aes_string(x = "clusters", y = "expression"), 
                    scale = "width", color = "transparent",) +
        geom_boxplot(aes_string(x = "clusters", y = "expression"),
                    width = 0.1, outlier.shape = NA, fill = "white") +
        theme_classic() +
        theme(
          panel.border = element_rect(colour = "black", fill = NA),
          panel.spacing = unit(0, "lines")
        ) + facet_grid(gene_padded ~ ., scales = "free_y") +
         theme( 
            axis.ticks = element_blank(),
            strip.text.y = element_text(angle = 0, hjust = 0),
            strip.background = element_blank()
         )

ggsave("vlnplot.pdf",plot = p3,width = 12,height = 7)
ggsave("vlnplot.png",plot = p3,width = 12,height = 7,dpi = 300)

