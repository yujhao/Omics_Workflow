# 在空间转录组分析中，欧式距离的意义是把 “二维坐标”→“生物学可解释的空间尺度”，用来：
# 定量描述细胞相对位置；
# 揭示微环境影响；
# 分析空间梯度效应；
# 在多样本之间建立比较基准。

# 计算所有细胞/spot 到某个指定中心细胞类型（例如 L6b）最近的距离，并且在密度图和空间切片上可视化

# 读取已经保存的 Seurat 对象
object = readRDS("./seurat_RCTD.rds")
# 查看对象内容（前6行）
head(object)

# 设置关注的中心细胞类型名称（这里是 "L6"）
center = "L6b"

# 设置细胞类型的列名（meta.data 中存储的分类结果）
celltype_col = "top1_celltype"

# -----------------------------
# 获取所有图像的空间坐标 (row, col)，合并成一个 dataframe
# object@images 里存的是 Seurat Spatial 相关数据
# -----------------------------
locat <- lapply(object@images, function(x) {
        {
            data <- x
        }@coordinates[, c("row", "col")]   # 提取坐标
    }) %>% dplyr::bind_rows()

# 计算所有 spots 之间的欧几里得距离矩阵
distest <- dist(locat, p = 2)

# ---------------------------
# 找到 meta.data 中，细胞类型包含 "center" (L6b) 的 barcode
# -----------------------------
center_barcode <- rownames(object@meta.data)[grep(center, unlist(object[[celltype_col]]))]

# 把 dist 对象转成矩阵，再转成 data.frame
distest2 <- data.frame(as.matrix(distest))

# 将列名中的 "." 替换成 "-"（保持和 barcode 格式一致）
names(distest2) <- gsub("\\.", "-", names(distest2))

# 保留 meta.data 中所有细胞到 center 细胞的距离
distest2 <- distest2[rownames(object@meta.data), center_barcode]


# -----------------------------
# 获取每个细胞的 sampleid（样本信息）
# -----------------------------
barcodelist <- Seurat::FetchData(object, "sampleid")
barcodelist$barcode <- colnames(object)

# 只保留中心细胞的 barcode + sampleid
center_barcode_list <- barcodelist[barcodelist$barcode %in% center_barcode, ]

# 按 sampleid 拆分成多个子表
barcodelist <- split(barcodelist, barcodelist$sampleid)
center_barcode_list <- split(center_barcode_list, center_barcode_list$sampleid)

# -----------------------------
# 计算每个细胞到最近的中心细胞的最小距离
# -----------------------------
dist_list <- list()
for (i in names(barcodelist)) {
        # 取样本 i 中的所有细胞到样本 i 中心细胞的距离矩阵
        data <- distest2[barcodelist[[i]]$barcode, center_barcode_list[[i]]$barcode]
        # 对每个细胞取最小距离
        dist_list[[i]] <- data.frame(mindist = apply(data, 1, min))
    }


# 合并所有样本的最小距离
dist <- dplyr::bind_rows(dist_list)

# 把这个最小距离保存到 Seurat meta.data 中
object[[paste0(center, "_distance")]] <- dist[, 1]


# -----------------------------
# 准备绘图数据：包含 celltype、距离、sample 信息
# -----------------------------
plot_data <- data.frame(
        celltype = factor(
            object@meta.data[, celltype_col],
            levels = unique(object@meta.data[, celltype_col])
        ),
        dist = dist[, 1],
        sample = object$sampleid
    )

# -----------------------------
# 图1：按 celltype 绘制到中心细胞的距离密度曲线
# facet_grid 按 sample 分面
# -----------------------------
p1_density <- ggplot2::ggplot(plot_data) +
        ggplot2::stat_density(ggplot2::aes(x = dist, colour = celltype), 
                              geom = "line", position = "identity", size = 0.3) +
        #ggplot2::scale_color_manual(values = colors) +
        ggplot2::theme_set(ggplot2::theme_bw()) +
        ggplot2::theme(panel.grid = ggplot2::element_blank()) +
        ggplot2::xlab(paste0("distance to ", center, " of ", Seurat::Images(object)[1])) +
        ggplot2::facet_grid(sample ~ .) +
        ggplot2::scale_y_continuous(limits = c(0, 1)) +
        ggplot2::scale_x_continuous(limits = c(0, 15))

# ----------------------------
# 图2：在空间切片上展示每个细胞的中心距离
# 每个 Seurat::Images(object) 代表一个样本切片
# -----------------------------
p1_slice <- list()
for (i in Seurat::Images(object)) {
        p1_slice[[i]] <- SpatialPlot(
            object,
            features = paste0(center, "_distance")
        )[[1]]
    }

# 合并所有切片图，排列在一行，共享图例
p1_slice <- do.call(
        ggpubr::ggarrange,
        c(
            p1_slice,
            list(
                nrow = 1,
                common.legend = TRUE,
                legend = "right",
                align = "hv"
            )
        )
    )