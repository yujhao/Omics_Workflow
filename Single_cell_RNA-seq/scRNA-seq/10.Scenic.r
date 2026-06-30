#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ==============================================================================
# SCENIC 转录因子调控网络分析脚本
# 包含 RSS 计算、CSI 分析、 Regulon 活性评估等功能
# ==============================================================================

# ========================== 加载包 ==========================

suppressPackageStartupMessages({
    library(RcisTarget)
    library(GENIE3)
    library(AUCell)
    library(SCENIC)
    library(tidyverse)
    library(Seurat)
    library(Matrix)
    library(doParallel)
    library(ggplot2)
    library(feather)
    library(DT)
    library(viridis)
    library(pheatmap)
    library(ggrepel)
    library(BiocParallel)
    library(reticulate)
    library(ComplexHeatmap)
    library(reshape2)
    library(data.table)
    library(tictoc)
    library(ff)
    library(doMC)
})

# ========================== 配置参数 ==========================

OUTPUT_DIR <- "10.SCENIC"
DATA_FILE <- "data_ob_v3.rds"
N_CORES <- 10

# 物种选择: "human" 或 "mouse"
SPECIES <- "mouse"

# 数据库路径（根据物种自动选择）
SCENIC_DB_BASE <- "/hwstorage/oe-scrna/jhyu/github_jhyu/Omics_Workflow/utils/data/SCENIC"

if (SPECIES == "human") {
    DB_DIR <- file.path(SCENIC_DB_BASE, "human")
    DBS <- c(
        "utils/data/SCENIC/human/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather",
        "utils/data/SCENIC/human/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
    )
    names(DBS) <- c("500bp", "10kb")
    ORG <- "hgnc"
} else if (SPECIES == "mouse") {
    DB_DIR <- file.path(SCENIC_DB_BASE, "mouse")
    DBS <- c(
        "mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather",
        "mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
    )
    names(DBS) <- c("500bp", "10kb")
    ORG <- "mgi"
} else {
    stop("物种不支持，请选择 'human' 或 'mouse'")
}

# ========================== 初始化 ==========================

dir.create(OUTPUT_DIR, showWarnings = FALSE)
setwd(OUTPUT_DIR)

# ========================== 辅助函数定义 ==========================

# 内部函数：熵计算
.H <- function(pVect) {
    pVect <- pVect[pVect > 0]
    -sum(pVect * log2(pVect))
}

# Jensen-Shannon 散度
calcJSD <- function(pRegulon, pCellType) {
    (.H((pRegulon + pCellType) / 2)) - ((.H(pRegulon) + .H(pCellType)) / 2)
}

# 单个 Regulon 的 RSS 计算
.calcRSS.oneRegulon <- function(pRegulon, pCellType) {
    jsd <- calcJSD(pRegulon, pCellType)
    1 - sqrt(jsd)
}

# 计算 RSS（Regulon Specificity Score）
calcRSS <- function(AUC, cellAnnotation, cellTypes = NULL) {
    if (any(is.na(cellAnnotation))) stop("NAs in annotation")
    if (any(class(AUC) == "aucellResults")) AUC <- getAUC(AUC)
    normAUC <- AUC / rowSums(AUC)
    if (is.null(cellTypes)) cellTypes <- unique(cellAnnotation)

    ctapply <- lapply
    rss <- ctapply(cellTypes, function(thisType) {
        sapply(rownames(normAUC), function(thisRegulon) {
            pRegulon <- normAUC[thisRegulon, ]
            pCellType <- as.numeric(cellAnnotation == thisType)
            pCellType <- pCellType / sum(pCellType)
            .calcRSS.oneRegulon(pRegulon, pCellType)
        })
    })

    rss <- do.call(cbind, rss)
    colnames(rss) <- cellTypes
    return(rss)
}

# RSS 热图内部函数
.plotRSS_heatmap <- function(rss, thr = NULL, row_names_gp = gpar(fontsize = 5),
                             order_rows = TRUE, cluster_rows = FALSE,
                             name = "RSS", verbose = TRUE, ...) {
    if (is.null(thr)) thr <- signif(quantile(rss, p = 0.97), 2)

    rssSubset <- rss[rowSums(rss > thr) > 0, ]
    rssSubset <- rssSubset[, colSums(rssSubset > thr) > 0]

    if (verbose) message("Showing regulons and cell types with any RSS > ", thr,
                        " (dim: ", nrow(rssSubset), "x", ncol(rssSubset), ")")

    if (order_rows) {
        maxVal <- apply(rssSubset, 1, which.max)
        rss_ordered <- rssSubset[0, ]
        for (i in 1:ncol(rssSubset)) {
            tmp <- rssSubset[which(maxVal == i), , drop = FALSE]
            tmp <- tmp[order(tmp[, i], decreasing = FALSE), , drop = FALSE]
            rss_ordered <- rbind(rss_ordered, tmp)
        }
        rssSubset <- rss_ordered
        cluster_rows <- FALSE
    }

    Heatmap(rssSubset, name = name, row_names_gp = row_names_gp,
            cluster_rows = cluster_rows, ...)
}

# 绘制 RSS 图
plotRSS <- function(rss, labelsToDiscard = NULL, zThreshold = 1,
                    cluster_columns = FALSE, order_rows = TRUE, thr = 0.01,
                    varName = "cellType",
                    col.low = "grey90", col.mid = "darkolivegreen3", col.high = "darkgreen",
                    revCol = FALSE, verbose = TRUE) {
    varSize <- "RSS"
    varCol <- "Z"
    if (revCol) {
        varSize <- "Z"
        varCol <- "RSS"
    }

    rssNorm <- scale(rss)
    rssNorm <- rssNorm[, which(!colnames(rssNorm) %in% labelsToDiscard)]
    rssNorm[rssNorm < 0] <- 0

    rssSubset <- rssNorm
    if (!is.null(zThreshold)) rssSubset[rssSubset < zThreshold] <- 0
    tmp <- .plotRSS_heatmap(rssSubset, thr = thr, cluster_columns = cluster_columns,
                           order_rows = order_rows, verbose = verbose)
    rowOrder <- rev(tmp@row_names_param$labels)
    rm(tmp)

    rss.df <- reshape2::melt(rss)
    colnames(rss.df) <- c("Topic", varName, "RSS")
    rssNorm.df <- reshape2::melt(rssNorm)
    colnames(rssNorm.df) <- c("Topic", varName, "Z")
    rss.df <- base::merge(rss.df, rssNorm.df)

    rss.df <- rss.df[which(!rss.df[, varName] %in% labelsToDiscard), ]
    if (nrow(rss.df) < 2) stop("Insufficient rows left to plot RSS.")

    rss.df <- rss.df[which(rss.df$Topic %in% rowOrder), ]
    rss.df[, "Topic"] <- factor(rss.df[, "Topic"], levels = rowOrder)

    p <- dotHeatmap(rss.df,
                   var.x = varName, var.y = "Topic",
                   var.size = varSize, min.size = 0.5, max.size = 5,
                   var.col = varCol, col.low = col.low, col.mid = col.mid, col.high = col.high)

    invisible(list(plot = p, df = rss.df, rowOrder = rowOrder))
}

# 绘制单个细胞类型的 RSS 排名图
plotRSS_oneSet <- function(rss, setName, n = 5) {
    rssThisType <- sort(rss[, setName], decreasing = TRUE)
    thisRss <- data.frame(regulon = names(rssThisType), rank = seq_along(rssThisType), rss = rssThisType)
    thisRss$regulon[(n + 1):nrow(thisRss)] <- NA

    ggplot(thisRss, aes(x = rank, y = rss)) +
        geom_point(color = "blue", size = 1) +
        ggtitle(setName) +
        geom_label_repel(aes(label = regulon),
                        box.padding = 0.35,
                        point.padding = 0.5,
                        segment.color = "grey50",
                        na.rm = TRUE) +
        theme_classic()
}

# Dot Heatmap 函数
dotHeatmap <- function(enrichmentDf,
                      var.x = "Topic", var.y = "ID",
                      var.col = "FC", col.low = "dodgerblue", col.mid = "floralwhite", col.high = "brown1",
                      var.size = "p.adjust", min.size = 1, max.size = 8, ...) {
    colorPal <- grDevices::colorRampPalette(c(col.low, col.mid, col.high))
    p <- ggplot(data = enrichmentDf, mapping = aes_string(x = var.x, y = var.y)) +
        geom_point(mapping = aes_string(size = var.size, color = var.col)) +
        scale_radius(range = c(min.size, max.size)) +
        scale_colour_gradientn(colors = colorPal(10)) +
        theme_bw() +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
              axis.text.x = element_text(angle = 90, hjust = 1), ...)
    return(p)
}

# CSI 模块活性计算
calc_csi_module_activity <- function(clusters_df, regulonAUC, metadata, cell_type_column) {
    metadata$cell_type <- metadata[, cell_type_column]
    cell_types <- unique(metadata$cell_type)
    regulons <- unique(clusters_df$regulon)
    regulonAUC_sub <- regulonAUC@assays@data@listData$AUC
    regulonAUC_sub <- regulonAUC_sub[regulons, ]

    csi_activity_matrix_list <- list()
    csi_cluster_activity <- data.frame(
        "csi_cluster" = c(),
        "mean_activity" = c(),
        "cell_type" = c()
    )

    for (ct in cell_types) {
        cell_type_aucs <- rowMeans(regulonAUC_sub[, rownames(subset(metadata, cell_type == ct))])
        cell_type_aucs_df <- data.frame(
            "regulon" = names(cell_type_aucs),
            "activity" = cell_type_aucs,
            "cell_type" = ct
        )
        csi_activity_matrix_list[[ct]] <- cell_type_aucs_df
    }

    for (ct in names(csi_activity_matrix_list)) {
        for (cluster in unique(clusters_df$csi_module)) {
            csi_regulon <- subset(clusters_df, csi_module == cluster)
            csi_regulon_activity <- subset(csi_activity_matrix_list[[ct]], regulon %in% csi_regulon$regulon)
            csi_activity_mean <- mean(csi_regulon_activity$activity)
            this_cluster_ct_activity <- data.frame(
                "csi_module" = cluster,
                "mean_activity" = csi_activity_mean,
                "cell_type" = ct
            )
            csi_cluster_activity <- rbind(csi_cluster_activity, this_cluster_ct_activity)
        }
    }

    csi_cluster_activity[is.na(csi_cluster_activity)] <- 0
    csi_cluster_activity_wide <- csi_cluster_activity %>%
        tidyr::spread(cell_type, mean_activity)

    rownames(csi_cluster_activity_wide) <- csi_cluster_activity_wide$csi_cluster
    csi_cluster_activity_wide <- as.matrix(csi_cluster_activity_wide[2:ncol(csi_cluster_activity_wide)])

    return(csi_cluster_activity_wide)
}

# CSI 计算
calculate_csi <- function(regulonAUC, calc_extended = FALSE, verbose = FALSE) {
    compare_pcc <- function(vector_of_pcc, pcc) {
        pcc_larger <- length(vector_of_pcc[vector_of_pcc > pcc])
        if (pcc_larger == length(vector_of_pcc)) {
            return(0)
        } else {
            return(length(vector_of_pcc))
        }
    }

    calc_csi <- function(reg, reg2, pearson_cor) {
        test_cor <- pearson_cor[reg, reg2]
        total_n <- ncol(pearson_cor)
        pearson_cor_sub <- subset(pearson_cor, rownames(pearson_cor) == reg | rownames(pearson_cor) == reg2)
        sums <- apply(pearson_cor_sub, MARGIN = 2, FUN = compare_pcc, pcc = test_cor)
        fraction_lower <- length(sums[sums == nrow(pearson_cor_sub)]) / total_n
        return(fraction_lower)
    }

    regulonAUC_sub <- regulonAUC@assays@data@listData$AUC

    if (calc_extended) {
        regulonAUC_sub <- subset(regulonAUC_sub, grepl("extended", rownames(regulonAUC_sub)))
    } else {
        regulonAUC_sub <- subset(regulonAUC_sub, !grepl("extended", rownames(regulonAUC_sub)))
    }

    regulonAUC_sub <- t(regulonAUC_sub)
    pearson_cor <- cor(regulonAUC_sub)
    pearson_cor_df <- as.data.frame(pearson_cor)
    pearson_cor_df$regulon_1 <- rownames(pearson_cor_df)
    pearson_cor_long <- pearson_cor_df %>%
        tidyr::gather(regulon_2, pcc, -regulon_1) %>%
        dplyr::mutate("regulon_pair" = paste(regulon_1, regulon_2, sep = "_"))

    regulon_names <- unique(colnames(pearson_cor))
    num_of_calculations <- length(regulon_names) * length(regulon_names)

    csi_regulons <- data.frame(matrix(nrow = num_of_calculations, ncol = 3))
    colnames(csi_regulons) <- c("regulon_1", "regulon_2", "CSI")

    f <- 0
    for (reg in regulon_names) {
        if (verbose) print(reg)
        for (reg2 in regulon_names) {
            f <- f + 1
            fraction_lower <- calc_csi(reg, reg2, pearson_cor)
            csi_regulons[f, ] <- c(reg, reg2, fraction_lower)
        }
    }
    csi_regulons$CSI <- as.numeric(csi_regulons$CSI)
    return(csi_regulons)
}

# CSI 模块可视化
plot_csi_modules <- function(csi_df, nclust = 10, font_size_regulons = 6) {
    csi_test_mat <- csi_df %>%
        tidyr::spread(regulon_2, CSI)

    future_rownames <- csi_test_mat$regulon_1
    csi_test_mat <- as.matrix(csi_test_mat[, 2:ncol(csi_test_mat)])
    rownames(csi_test_mat) <- future_rownames

    pheatmap(csi_test_mat,
            show_colnames = FALSE,
            color = viridis(n = 10),
            cutree_cols = nclust,
            cutree_rows = nclust,
            fontsize_row = font_size_regulons,
            cluster_cols = TRUE,
            cluster_rows = TRUE,
            treeheight_row = 20,
            treeheight_col = 20,
            clustering_distance_rows = "euclidean",
            clustering_distance_cols = "euclidean")
}

# 大矩阵相关计算
bigcor <- function(x, y = NULL, size = 2000, cores = 8, verbose = TRUE, ...) {
    tictoc::tic()
    if (!is.null(y) & NROW(x) != NROW(y)) stop("'x' and 'y' must have compatible dimensions!")

    NCOL <- ncol(x)
    if (!is.null(y)) YCOL <- NCOL(y)
    REST <- NCOL %% size
    LARGE <- NCOL - REST
    NBLOCKS <- NCOL %/% size

    if (is.null(y)) {
        resMAT <- ff::ff(vmode = "double", dim = c(NCOL, NCOL))
    } else {
        resMAT <- ff::ff(vmode = "double", dim = c(NCOL, YCOL))
    }

    GROUP <- rep(1:NBLOCKS, each = size)
    if (REST > 0) GROUP <- c(GROUP, rep(NBLOCKS + 1, REST))
    SPLIT <- split(1:NCOL, GROUP)
    COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
    COMBS <- t(apply(COMBS, 1, sort))
    COMBS <- unique(COMBS)
    if (!is.null(y)) COMBS <- cbind(1:length(SPLIT), rep(1, length(SPLIT)))

    ncore <- min(future::availableCores(), cores)
    doMC::registerDoMC(cores = ncore)

    results <- foreach(i = 1:nrow(COMBS)) %dopar% {
        COMB <- COMBS[i, ]
        G1 <- SPLIT[[COMB[1]]]
        G2 <- SPLIT[[COMB[2]]]

        if (is.null(y)) {
            if (verbose) message("bigcor: ", sprintf("#%d:Block %s and Block %s (%s x %s) ... ",
                                                    i, COMB[1], COMB[2], length(G1), length(G2)))
            flush.console()
            RES <- do.call("cor", list(x = x[, G1], y = x[, G2], ...))
            resMAT[G1, G2] <- RES
            resMAT[G2, G1] <- t(RES)
        } else {
            if (verbose) message("bigcor: ", sprintf("#%d:Block %s and 'y' (%s x %s) ... ",
                                                    i, COMB[1], length(G1), YCOL))
            flush.console()
            RES <- do.call("cor", list(x = x[, G1], y = y, ...))
            resMAT[G1, ] <- RES
        }
    }

    if (is.null(y)) {
        resMAT <- resMAT[1:ncol(x), 1:ncol(x)]
        colnames(resMAT) <- colnames(x)
        rownames(resMAT) <- colnames(x)
    } else {
        resMAT <- resMAT[1:ncol(x), 1:ncol(y)]
        colnames(resMAT) <- colnames(x)
        rownames(resMAT) <- colnames(y)
    }

    tictoc::toc()
    return(resMAT)
}

# 写入 GMT 文件
write.gmt <- function(geneSet, gmt_file = "kegg2symbol.gmt") {
    sink(gmt_file)
    for (i in 1:length(geneSet)) {
        cat(names(geneSet)[i])
        cat('\t')
        cat(paste(geneSet[[i]], collapse = '\t'))
        cat('\n')
    }
    sink()
}

# ========================== 主分析流程 ==========================

# 初始化 SCENIC
scenicOptions <- initializeScenic(org = ORG, dbDir = DB_DIR, dbs = DBS, nCores = N_CORES)

# 读取数据
scRNA_mnn <- readRDS(file.path("..", DATA_FILE))

# 基因过滤
minCell4gene <- round(0.01 * ncol(scRNA_mnn))
exprMat <- scRNA_mnn@assays$RNA@counts
genesKept <- geneFiltering(
    as.matrix(exprMat),
    scenicOptions = scenicOptions,
    minCountsPerGene = 1,
    minSamples = minCell4gene
)
exprMat_filtered <- exprMat[genesKept, ]

# 获取 TF 名称
tf_names <- getDbTfs(scenicOptions)
tf_names <- CaseMatch(search = tf_names, match = rownames(scRNA_mnn))

# 运行 GENIE3
arb.algo <- import("arboreto.algo")
adjacencies <- arb.algo$grnboost2(
    as.data.frame(t(as.matrix(exprMat_filtered))),
    tf_names = tf_names,
    verbose = TRUE,
    seed = 123L
)
colnames(adjacencies) <- c("TF", "Target", "weight")
saveRDS(adjacencies, file = getIntName(scenicOptions, "genie3ll"))

# 计算相关矩阵
corrMat <- bigcor(t(as.matrix(exprMat_filtered)), size = 2000, cores = 20, method = "spearman")
saveRDS(corrMat, file = getIntName(scenicOptions, "corrMat"))

# 运行 SCENIC 流程
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod = "top10perTarget")

# 重新初始化并评分
scenicOptions <- initializeScenic(org = ORG, dbDir = DB_DIR, dbs = DBS, nCores = 1)
runSCENIC_3_scoreCells(scenicOptions, log2(as.matrix(exprMat_filtered) + 1))

# ========================== 结果提取与保存 ==========================

# 提取 Regulon AUC
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC_mat <- AUCell::getAUC(regulonAUC)
rownames(regulonAUC_mat) <- gsub("_", "-", rownames(regulonAUC_mat))
regulonAUC_mat_out <- regulonAUC_mat[-grep(pattern = "-extended", rownames(regulonAUC_mat)), ]

# 保存 Regulon 活性
write.table(
    as.data.frame(regulonAUC_mat_out) %>% tibble::rownames_to_column(var = "regulon"),
    file.path("regulon_activity.xls"),
    sep = "\t", col.names = TRUE, row.names = FALSE
)

# 添加到 Seurat 对象
scRNA_mnn[["SCENIC"]] <- CreateAssayObject(counts = regulonAUC_mat)
scRNA_mnn <- ScaleData(scRNA_mnn, assay = "SCENIC")
scRNA_mnn@tools$RunAUCell <- regulonAUC

# 保存 TF 靶标信息
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
write.table(
    regulonTargetsInfo,
    file.path("0.1.TF_target_enrichment_annotation.xls"),
    sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

# 保存 Regulon GMT
regulons <- loadInt(scenicOptions, "regulons")
sub_regulons <- gsub(" .*", "", rownames(regulonAUC_mat_out))
regulons <- regulons[sub_regulons]
write.gmt(regulons, gmt_file = "0.2.regulon_annotation.xls")

# ========================== 热图可视化 ==========================

# 准备注释信息
Idents(object = scRNA_mnn) <- scRNA_mnn@meta.data$clusters
cellInfo <- scRNA_mnn@meta.data
col_anno <- as.data.frame(scRNA_mnn@meta.data) %>% tibble::rownames_to_column(var = "barcodes")
col_anno <- col_anno[, c("barcodes", "clusters")]
col_anno <- col_anno %>% dplyr::arrange(clusters) %>% tibble::column_to_rownames(var = "barcodes")
regulonAUC_plotdata <- regulonAUC_mat_out[, rownames(col_anno)]

# 颜色设置
bks <- unique(c(seq(-2.5, 0, length = 100), seq(0, 2.5, length = 100)))
color_use <- c("#7fc97f", "#beaed4", "#fdc086", "#386cb0", "#f0027f", "#a34e3b",
              "#666666", "#1b9e77", "#d95f02", "#7570b3", "#d01b2a", "#43acde")
cluster_levels <- as.character(unique(col_anno$clusters))
names(color_use) <- cluster_levels
annotation_colors <- list(clusters = color_use)

# Regulon 活性热图
pdf("1.1.regulon_activity_heatmap_groupby_cells.pdf", height = 8)
pheatmap::pheatmap(
    regulonAUC_plotdata,
    scale = "row",
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    show_colnames = FALSE,
    color = colorRampPalette(c("#406AA8", "white", "#D91216"))(200),
    annotation_col = col_anno,
    annotation_colors = annotation_colors,
    treeheight_col = 10,
    border_color = NA,
    breaks = bks,
    fontsize_row = 6
)
dev.off()

# 按聚类平均活性
regulonActivity_byclusters <- sapply(split(rownames(cellInfo), cellInfo$clusters), function(cells)
    rowMeans(getAUC(regulonAUC)[, cells]))
regulonActivity_byclusters_Scaled <- t(scale(t(regulonActivity_byclusters), center = TRUE, scale = TRUE))

df <- as.data.frame(regulonActivity_byclusters_Scaled) %>% tibble::rownames_to_column(var = "regulon")
write.table(
    df,
    file.path("1.2.centered_regulon_activity_groupby_design.xls"),
    sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

# Top10 Regulon 热图
regulonAUC_plotdata <- regulonActivity_byclusters_Scaled[1:10, ]
pdf("1.3.regulon_activity_heatmap.pdf")
pheatmap::pheatmap(
    regulonAUC_plotdata,
    cellwidth = 18,
    cellheight = 18,
    color = colorRampPalette(c("#406AA8", "white", "#D91216"))(299),
    angle_col = 45,
    treeheight_col = 20,
    treeheight_row = 20,
    border_color = NA
)
dev.off()

# ========================== RSS 分析 ==========================

# 计算 RSS
rss <- calcRSS(AUC = getAUC(regulonAUC), cellAnnotation = cellInfo[colnames(regulonAUC), "clusters"])

# 绘制 RSS 排名图
pdf("2.2.RSS_ranking_plot.pdf")
setName <- "1"
n <- 5
rssThisType <- sort(rss[, setName], decreasing = TRUE)
thisRss <- data.frame(regulon = names(rssThisType), rank = seq_along(rssThisType), rss = rssThisType)
thisRss$regulon[(n + 1):nrow(thisRss)] <- NA

p4 <- ggplot(thisRss, aes(x = rank, y = rss)) +
    geom_point(color = "grey50", size = 1) +
    ggtitle(setName) +
    geom_point(data = subset(thisRss, rank < n + 1), color = "red", size = 2) +
    geom_label_repel(aes(label = regulon),
                    box.padding = 0.35,
                    point.padding = 0.5,
                    segment.color = "grey50",
                    na.rm = TRUE) +
    theme_classic()
print(p4)
dev.off()

# 保存 RSS 结果
write.table(
    thisRss[1:5, ],
    file.path("2.1.regulon_RSS_annotation.xls"),
    sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)

# ========================== CSI 分析 ==========================

# 计算 CSI
regulons_csi <- calculate_csi(regulonAUC, calc_extended = FALSE)

# CSI 矩阵
csi_csi_wide <- regulons_csi %>% tidyr::spread(regulon_2, CSI)
future_rownames <- csi_csi_wide$regulon_1
csi_csi_wide <- as.matrix(csi_csi_wide[, 2:ncol(csi_csi_wide)])
rownames(csi_csi_wide) <- future_rownames

# 层次聚类
regulons_hclust <- hclust(dist(csi_csi_wide, method = "euclidean"))
nclust <- 4
clusters <- cutree(regulons_hclust, k = nclust)
clusters_df <- data.frame("regulon" = names(clusters), "csi_module" = clusters)

# CSI 模块活性
cellinfo <- scRNA_mnn@meta.data
csi_cluster_activity_wide <- calc_csi_module_activity(
    clusters_df, regulonAUC, cellinfo, cell_type_column = "seurat_clusters"
)
rownames(csi_cluster_activity_wide) <- paste0("module", c(1:nclust))

# CSI 活性热图
plot <- pheatmap::pheatmap(
    csi_cluster_activity_wide,
    show_colnames = TRUE,
    show_rownames = TRUE,
    scale = "row",
    color = viridis::viridis(n = 10),
    cellwidth = 24,
    cellheight = 24,
    cluster_cols = TRUE,
    cluster_rows = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean"
)
ggsave(
    "3.3.csi_module_activity_heatmap.pdf",
    plot = plot, width = 8, height = 8, dpi = 1000, limitsize = FALSE
)

# CSI 相关热图
plot <- plot_csi_modules(regulons_csi, nclust = nclust)
pdf("3.2.regulons_csi_correlation_heatmap.pdf")
print(plot)
dev.off()

# 保存 CSI 模块注释
write.table(
    clusters_df,
    "3.1.csi_module_annotation.xls",
    sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)
