#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ==============================================================================
# CellPhoneDB 细胞通讯分析脚本
# 分析细胞间的配体-受体相互作用
# ==============================================================================

# ========================== 加载包 ==========================

suppressPackageStartupMessages({
    library(annotables)
    library(future)
    library(future.apply)
    library(homologene)
    library(dplyr)
    library(optparse)
    library(Seurat)
    library(ggplot2)
    library(tidyr)
    library(network)
    library(igraph)
    library(circlize)
    library(ComplexHeatmap)
    library(RColorBrewer)
    library(vroom)
    library(glue)
    library(plyr)
})

# 加载可视化源文件
source("utils/plotting/vis_cellcomm_source.r")

# ========================== 配置参数 ==========================

OUTPUT_DIR <- "12.Cellphonedb"
DATA_FILE <- "data_ob_v3.rds"
COLUMN_4_CELL <- "clusters"
THREADS <- 10
ITERATIONS <- 1000
THRESHOLD <- 0.1
PVALUE <- 0.05
GENETYPE <- "gene_name"
DATABASE <- "utils/data/cellphonedb.zip"
IS_ONLYSIG <- TRUE
SPECIES <- "mouse"
GROUPBY <- NULL

# ========================== 初始化 ==========================

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# 设置并行
options(future.globals.maxSize = Inf)
plan("multicore", workers = min(availableCores(), THREADS))

# ========================== 数据加载 ==========================

seurat_ob <- readRDS(DATA_FILE)

# 准备样本列表
sp <- list()
sp[["all"]] <- seurat_ob

# 去除空因子水平
if (class(seurat_ob@meta.data[, COLUMN_4_CELL]) == "factor") {
    sp <- lapply(sp, function(x) {
        x@meta.data[, COLUMN_4_CELL] <- droplevels(x@meta.data[, COLUMN_4_CELL])
        return(x)
    })
}

# 获取表达矩阵
counts <- lapply(sp, function(x) {
    counts <- GetAssayData(x, slot = "data")
    return(counts)
})

# 获取细胞注释
metadata <- lapply(sp, function(x) {
    if (is.null(COLUMN_4_CELL)) {
        metadata <- Seurat::FetchData(x, vars = "ident") %>%
            tibble::rownames_to_column(var = "Cell") %>%
            dplyr::rename("cell_type" = ident)
    } else {
        metadata <- Seurat::FetchData(x, vars = COLUMN_4_CELL) %>%
            tibble::rownames_to_column(var = "Cell")
        metadata[, 2] <- gsub("[^[:alnum:]]", "_", metadata[, 2])
        metadata[, 2] <- gsub("__+", "_", metadata[, 2])
        colnames(metadata) <- c("Cell", "cell_type")
        if (COLUMN_4_CELL == "clusters") {
            metadata$cell_type <- paste0("C", metadata$cell_type)
        }
    }
    return(metadata)
})

rm(seurat_ob, sp)
gc()

# ========================== 物种转换 ==========================

if (SPECIES == "mouse") {
    mouse2human <- homologene::human2mouse(annotables::grch38$symbol) %>%
        dplyr::select(mouseGene, humanGene)

    print("物种为小鼠，counts 基因名将替换为人的同源基因.")
    write.table(
        mouse2human,
        file = file.path(OUTPUT_DIR, "mouse2human.homologene.xls"),
        quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t"
    )

    counts <- lapply(counts, function(x) {
        genes2use <- intersect(
            names(Seurat::CaseMatch(rownames(x), mouse2human$mouseGene)),
            rownames(x)
        )
        counts <- x[genes2use, ]
        rownames(counts) <- plyr::mapvalues(
            rownames(counts),
            from = mouse2human$mouseGene,
            to = mouse2human$humanGene,
            warn_missing = FALSE
        )
        return(counts)
    })
}

# ========================== 运行 CellPhoneDB ==========================

i <- "all"
if (is.null(GROUPBY)) {
    print("开始做不分组 cellphonedb 分析.")
    tempdir <- OUTPUT_DIR
} else {
    print(paste0("开始做分组", i, "的 cellphonedb 分析."))
    tempdir <- file.path(OUTPUT_DIR, i)
    dir.create(tempdir)
}

print("以下细胞群将进入细胞通讯分析：")
print(table(metadata[[i]]$cell_type))
print(paste0(c("基因数是：", "细胞数是："), dim(counts[[i]])))

# 输出计数矩阵
counts_out <- tibble::rownames_to_column(as.data.frame(counts[[i]]), var = "Gene")
vroom::vroom_write(x = counts_out, path = file.path(tempdir, "counts.tsv"),
                   col_names = TRUE, quote = "none", delim = "\t")
write.table(
    metadata[[i]],
    file = file.path(tempdir, "metadata.tsv"),
    quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t"
)

# 构建命令
outdir <- file.path(tempdir, "out")
docellphone <- "docellphone.py"
cmd <- glue::glue(
    'echo "group {i} begin"; cd {tempdir} && module purge && cpdb && python {docellphone}',
    ' --counts={tempdir}/counts.tsv',
    ' --metadata={tempdir}/metadata.tsv',
    ' --genetype={GENETYPE}',
    ' --iterations={ITERATIONS}',
    ' --threshold={THRESHOLD}',
    ' --threads={THREADS}',
    ' --pvalue={PVALUE}',
    ' --database {DATABASE}',
    ' --microenvs={""}'
)

print(glue::glue("{cmd}"))
system(cmd)
print(paste0("group ", i, " cellphonedb execution finished."))
print(paste0("结果路径：", tempdir))

# ========================== 结果处理 ==========================

# 重命名输出文件
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

# 整理结果文件
for (file in c("means.txt", "pvalues.txt", "significant_means.txt")) {
    temp_data <- read.delim(file.path(outdir, file), header = TRUE, sep = "\t")
    temp_data <- temp_data %>% dplyr::select(-c("directionality", "classification"))
    colnames(temp_data) <- gsub("\\.", "|", colnames(temp_data))
    write.table(
        temp_data,
        file.path(outdir, file),
        sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
    )
}

# 解析结果
cellphonedb <- ParseCpdb(outdir, pvalue = PVALUE)
write.table(
    cellphonedb$ligrec,
    file.path(tempdir, "cell_comm_annotation.xls"),
    sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)
saveRDS(cellphonedb, file.path(tempdir, "cellphonedb_results.rds"))

# ========================== 可视化 ==========================

data <- read.table("cell_comm_annotation.xls", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

if (is.numeric(sort(unique(data$ligand_cell))) && is.numeric(sort(unique(data$receptor_cell)))) {
    data$ligand_cell <- factor(data$ligand_cell, levels = sort(unique(data$ligand_cell)))
    data$receptor_cell <- factor(data$receptor_cell, levels = sort(unique(data$receptor_cell)))
}

data$group <- "all"

# 统计通讯对数
net <- data %>%
    filter(pval < 0.05) %>%
    group_by(ligand_cell, receptor_cell) %>%
    dplyr::summarize(n = n()) %>%
    dplyr::rename(significant_pairs = n)

write.table(
    net,
    file.path(OUTPUT_DIR, "cell_comm_summary.xls"),
    quote = FALSE, sep = "\t", row.names = FALSE
)

# 设置颜色
celltypes <- unique(c(as.character(data$receptor_cell), as.character(data$ligand_cell)))
colx <- c("#7fc97f", "#beaed4", "#fdc086", "#386cb0")
names(colx) <- celltypes

# Dotplot
out <- LRDotplot(data = data, is_onlySig = IS_ONLYSIG, topn = 5, xsize = 15, xangle = 45, remove.isolate = TRUE)
ggdot <- out[[1]]
plot_data <- out[[2]]

png_width <- length(unique(plot_data$clusters))
png_height <- length(unique(plot_data$pair))
max_length_width <- max(nchar(as.character(plot_data$clusters)))
max_length_height <- max(nchar(as.character(plot_data$pair)))

if (max_length_width > 40) {
    width <- png_width * 0.45 + max_length_width * 1.5 + 8
} else {
    width <- png_width * 0.45 + 6
}
height <- png_height * 0.4 + 3

ggsave(
    file.path(OUTPUT_DIR, "cell_comm_dotplot.pdf"),
    limitsize = FALSE, plot = ggdot, width = width, height = height, bg = "white"
)
ggsave(
    file.path(OUTPUT_DIR, "cell_comm_dotplot.png"),
    limitsize = FALSE, plot = ggdot, width = png_width * 0.4 + 6, height = png_height * 0.3
)

# Network
datax <- data
pdf(
    file.path(OUTPUT_DIR, "cell_comm_network.pdf"),
    width = 6.5 + max(nchar(unique(as.character(c(data$ligand_cell, data$receptor_cell))))) * 0.4,
    height = 7
)
LRNetwork(
    data %>% filter(pval < 0.05),
    col = colx,
    edge.label.cex = 0.3,
    edge.max.width = 5,
    arrow.width = 0.5,
    vertex.label.cex = 0.5,
    vertex.size = 10
)
dev.off()

png(
    file.path(OUTPUT_DIR, "cell_comm_network.png"),
    width = 6.5 + max(nchar(unique(as.character(c(data$ligand_cell, data$receptor_cell))))) * 0.4,
    height = 7, res = 500, units = "in"
)
LRNetwork(
    data %>% filter(pval < 0.05),
    col = colx,
    edge.label.cex = 0.3,
    edge.max.width = 5,
    arrow.width = 0.5,
    vertex.label.cex = 0.5,
    vertex.size = 10
)
dev.off()

# Circos plot
pdf(file.path(OUTPUT_DIR, "cell_comm_circos_plot.pdf"), width = 20, height = 10)
LRCircos(
    data,
    gap.degree = 0.05,
    cell_col = colx,
    screenvar = "expr",
    topn = 5,
    labels.cex = 0.35,
    link.lwd = 1,
    arr.length = 0.1
)
dev.off()

png(
    file.path(OUTPUT_DIR, "cell_comm_circos_plot.png"),
    width = 20, height = 10, res = 500, units = "in"
)
LRCircos(
    data,
    gap.degree = 0.05,
    cell_col = colx,
    screenvar = "expr",
    topn = 5,
    labels.cex = 0.35,
    link.lwd = 1,
    arr.length = 0.1
)
dev.off()

# Chord diagram
pdf(file.path(OUTPUT_DIR, "cell_comm_chorddiagram_plot.pdf"), width = 9, height = 7)
par(cex = 0.6, font = 2)
LRChorddiagram(
    data,
    grid_orbit_h = 0.02,
    label_h = 0.04,
    diffHeight = 4,
    output_dir = OUTPUT_DIR,
    colx = colx,
    extra_legend = TRUE
)
dev.off()

png(
    file.path(OUTPUT_DIR, "cell_comm_chorddiagram_plot.png"),
    width = 9, height = 7, res = 96, units = "in"
)
par(cex = 0.6, font = 2)
LRChorddiagram(
    data,
    grid_orbit_h = 0.02,
    label_h = 0.04,
    diffHeight = 4,
    output_dir = OUTPUT_DIR,
    colx = colx,
    extra_legend = TRUE
)
dev.off()
circos.clear()

# Heatmap
pdf(file.path(OUTPUT_DIR, "cell_comm_heatmap_plot.pdf"))
print(LRHeatmap(data))
dev.off()

png(
    file.path(OUTPUT_DIR, "cell_comm_heatmap_plot.png"),
    width = 7, height = 7, res = 96, units = "in"
)
print(LRHeatmap(data))
dev.off()

# Barplot
out <- LRBarplot(data = data, bar_width = 0.6, colx = colx) + scale_y_continuous(expand = c(0, 0))
ggsave(
    file.path(OUTPUT_DIR, "cell_comm_histogram_plot.pdf"),
    plot = out,
    width = ifelse(length(celltypes) < 7, 7, length(celltypes)),
    bg = "white"
)
ggsave(
    file.path(OUTPUT_DIR, "cell_comm_histogram_plot.png"),
    plot = out,
    width = ifelse(length(celltypes) < 7, 7, length(celltypes)),
    dpi = 1000,
    bg = "white"
)
