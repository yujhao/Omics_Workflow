#5 RCTD
dir.create("5.RCTD")
data_ob <- readRDS("MNN.rds")
counts_ref = GetAssayData(data_ob, 
                          assay = "RNA",
                          slot = "counts")
cell_types_ref <- data_ob@meta.data %>%  dplyr::select("subclass") %>%
  dplyr::mutate(cell_type =str_replace_all(!!!rlang::syms("subclass"), "-", "_") %>%   
  as.factor()) %>% .$cell_type
cell_types_ref <- setNames(cell_types_ref,
                           data_ob@meta.data  %>% rownames)
nUMI_ref <- data_ob@meta.data %>% dplyr::select("nCount_RNA") %>% .$nCount_RNA
nUMI_ref <- setNames(nUMI_ref,
                     data_ob@meta.data %>% rownames)
reference = spacexr::Reference(counts_ref,
                               cell_types_ref,
                               nUMI_ref)
seurat_to_SpatialRNA = function(st_ob, st_assay) {
  images = Seurat::Images(st_ob)  
  coords = do.call(rbind, lapply(images, function(slice) {
    spatial_coord = data.frame(st_ob[[slice]]@coordinates)
    colnames(spatial_coord) = c("in_tissue", "y", "x", "pxl_col_in_fullres", "pxl_row_in_fullres")
    spatial_coord$y = -1 * (spatial_coord$y) + 78
    spatial_coord = spatial_coord %>% dplyr::select(c(x, y))
    return(spatial_coord)
  }))
  counts = Seurat::GetAssayData(st_ob, assay = st_assay, slot = "counts")
  puck = spacexr::SpatialRNA(coords, counts) 
  return(puck)
}
spatialRNA = seurat_to_SpatialRNA(data_ob, "Spatial")
RCTD_ob = spacexr::create.RCTD(
  spatialRNA=spatialRNA,      
  reference=reference,        
  test_mode = FALSE,          
  gene_cutoff = 1.25e-4,     
  fc_cutoff = 0.5,           
  gene_cutoff_reg = 2e-04,  
  fc_cutoff_reg = 0.75, 
  UMI_max = 2e+05,
  counts_MIN = 10,    
  UMI_min_sigma = 0,       
  class_df = NULL,           
  CELL_MIN_INSTANCE = 1 ,     
  cell_type_names = NULL,     
  cell_type_profiles = NULL, 
  MAX_MULTI_TYPES = 4,        
  keep_reference = F,        
  CONFIDENCE_THRESHOLD = 10, 
  DOUBLET_THRESHOLD = 25,    
  max_cores = 20 )
myRCTD = spacexr::run.RCTD(RCTD_ob,
                           doublet_mode = "full")
norm_weights = normalize_weights(myRCTD@results$weights)%>% as.data.frame
data_ob@misc$RCTD_result = norm_weights
maxn = function(n) function(x) order(x, decreasing = TRUE)[!is.na(x)][n]
metadata = norm_weights %>% 
  dplyr::mutate(
    top1_celltype = apply(., 1, function(x) names(x)[maxn(1)(x)]),
    top2_celltype = apply(., 1, function(x) names(x)[maxn(2)(x)])
  )
data_ob = Seurat::AddMetaData(data_ob, metadata = metadata)
mat=normalize_weights(myRCTD@results$weights)%>%as.data.frame()
p = SPOTlight::spatial_scatterpie(se_obj = data_ob,
                                  cell_types_all = colnames(mat),
                                  img_path = "tissue_lowres_image.png",
                                  pie_scale = 0.4)
ggsave("5.RCTD/Celltype_RCTD.pdf", 
       plot = p,
       device = "pdf", 
       width = 14, 
       height = 12,
       dpi = 300)
p = SpatialDimPlot(data_ob,group.by = "top1_celltype")+
  theme( legend.key = element_rect(fill = NA, color = NA),
         legend.text = element_text(size = 10),
         plot.title = element_text(face = "bold", hjust = 0.5))
ggsave("5.RCTD/top1_celltype_RCTD.pdf",
       plot = p, 
       device = "pdf", 
       width = 14,
       height = 12,
       dpi = 300)
