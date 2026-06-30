#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ==============================================================================
# InferCNV 拷贝数变异分析脚本
# 检测单细胞数据中的拷贝数变异
# ==============================================================================

library(infercnv)
library(plyranges)
library(ComplexHeatmap)

# ========================== 配置参数 ==========================

OUTPUT_DIR <- "7.InferCNV"
GTF_FILE <- "genes.gtf"
REF_GROUP <- c("T_cells")  # 参考细胞类型
CUTOFF <- 0.1
NUM_THREADS <- 2

# ========================== 初始化 ==========================

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# ========================== 数据准备 ==========================

# 读取 GTF 文件获取基因位置信息
gtf <- plyranges::read_gff(GTF_FILE)
gene.chr <- gtf %>%
    plyranges::filter(type == "gene" & gene_name %in% rownames(scRNA_mnn)) %>%
    as.data.frame() %>%
    dplyr::select(gene_name, seqnames, start, end) %>%
    dplyr::distinct(gene_name, .keep_all = TRUE) %>%
    dplyr::mutate(seqnames = seqnames)

# 获取表达矩阵和细胞注释
count_mat <- GetAssayData(scRNA_mnn, "counts")
cellanno <- FetchData(scRNA_mnn, vars = "new_celltype") %>%
    tibble::rownames_to_column(var = "cellbarcode")

# 保存临时文件
tempdir <- tempdir()
cnv_celltyping <- file.path(tempdir, "cnv_celltype_group.xls")
write.table(cellanno, cnv_celltyping, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

gene_order_f <- file.path(tempdir, "gene_order_file.xls")
write.table(gene.chr, gene_order_f, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

# ========================== 运行 InferCNV ==========================

# 创建 InferCNV 对象
infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = count_mat,
    annotations_file = cnv_celltyping,
    delim = "\t",
    gene_order_file = gene_order_f,
    ref_group_names = REF_GROUP
)

# 运行分析
infercnv_obj <- infercnv::run(
    infercnv_obj,
    cutoff = CUTOFF,
    analysis_mode = "subclusters",
    tumor_subcluster_pval = 0.05,
    hclust_method = "ward.D2",
    out_dir = OUTPUT_DIR,
    num_threads = NUM_THREADS,
    cluster_by_groups = TRUE,
    denoise = TRUE,
    HMM = TRUE
)

# ========================== 可视化 ==========================

# 绘制热图
pdf(file.path(OUTPUT_DIR, "heatmap.pdf"), width = 18, height = 12)
ComplexHeatmap::Heatmap(
    t(as.matrix(infercnv_obj@expr.data)),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    name = "CNV level",
    use_raster = TRUE,
    raster_quality = 4
)
dev.off()

# ========================== CNV 评分计算 ==========================

# 将 CNV 结果添加到 Seurat 对象
scRNA_mnn[["CNV"]] <- CreateAssayObject(data = infercnv_obj@expr.data)

# 计算每个细胞的 CNV 水平
infercnv_level <- apply(as.data.frame(t(infercnv_obj@expr.data)), 1, function(x) {
    x[is.na(x)] <- 0
    return(sum(x))
})
infercnv_level <- round(
    scales::rescale(infercnv_level / nrow(infercnv_obj@expr.data), c(1, 100)),
    0
)
infercnv_level <- infercnv_level[Cells(scRNA_mnn)]
scRNA_mnn@meta.data$cnv_level <- infercnv_level

# 绘制 CNV 水平图
p1 <- FeaturePlot(scRNA_mnn, features = "cnv_level")
p2 <- VlnPlot(scRNA_mnn, features = "cnv_level", group.by = "orig.ident", pt.size = 0)

ggsave(file.path(OUTPUT_DIR, "featureplot.png"), p1, height = 5, width = 5)
ggsave(file.path(OUTPUT_DIR, "vlnplot.png"), p2, height = 4, width = 7)
