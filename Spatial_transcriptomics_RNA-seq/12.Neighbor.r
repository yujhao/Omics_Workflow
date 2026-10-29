# In spatial transcriptomics analysis, the Euclidean distance converts '2D coordinates'
# into a biologically interpretable spatial scale, used to:
# quantitatively describe relative cell positions;
# reveal microenvironment effects;
# analyze spatial gradient effects;
# establish a comparison baseline across multiple samples.

# Compute the nearest distance from all cells/spots to a specified center cell type (e.g., L6b), and visualize on density plots and spatial slices.

# Read the previously saved Seurat object
object = readRDS("./seurat_RCTD.rds")
# Inspect object contents (first 6 rows)
head(object)

# Set the center cell type of interest (here "L6")
center = "L6b"

# Set the cell-type column name (classification results stored in meta.data)
celltype_col = "top1_celltype"

# -----------------------------
# Retrieve the spatial coordinates (row, col) of all images and combine into a dataframe
# object@images stores Seurat Spatial-related data
# -----------------------------
locat <- lapply(object@images, function(x) {
        {
            data <- x
        }@coordinates[, c("row", "col")]   # extract coordinates
    }) %>% dplyr::bind_rows()

# Compute the Euclidean distance matrix between all spots
distest <- dist(locat, p = 2)

# ---------------------------
# Find barcodes in meta.data whose cell type contains "center" (L6b)
# -----------------------------
center_barcode <- rownames(object@meta.data)[grep(center, unlist(object[[celltype_col]]))]

# Convert the dist object to a matrix, then to a data.frame
distest2 <- data.frame(as.matrix(distest))

# Replace "." in column names with "-" (to match the barcode format)
names(distest2) <- gsub("\\.", "-", names(distest2))

# Keep distances from all cells to center cells in meta.data
distest2 <- distest2[rownames(object@meta.data), center_barcode]


# -----------------------------
# Get the sampleid (sample info) for each cell
# -----------------------------
barcodelist <- Seurat::FetchData(object, "sampleid")
barcodelist$barcode <- colnames(object)

# Keep only the center cells' barcode + sampleid
center_barcode_list <- barcodelist[barcodelist$barcode %in% center_barcode, ]

# Split into multiple subtables by sampleid
barcodelist <- split(barcodelist, barcodelist$sampleid)
center_barcode_list <- split(center_barcode_list, center_barcode_list$sampleid)

# -----------------------------
# Compute the minimum distance from each cell to the nearest center cell
# -----------------------------
dist_list <- list()
for (i in names(barcodelist)) {
        # Get the distance matrix from all cells in sample i to sample i's center cells
        data <- distest2[barcodelist[[i]]$barcode, center_barcode_list[[i]]$barcode]
        # Take the minimum distance for each cell
        dist_list[[i]] <- data.frame(mindist = apply(data, 1, min))
    }


# Combine the minimum distances across all samples
dist <- dplyr::bind_rows(dist_list)

# Save this minimum distance into the Seurat meta.data
object[[paste0(center, "_distance")]] <- dist[, 1]


# -----------------------------
# Prepare plotting data: include celltype, distance, and sample information
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
# Figure 1: density curve of distance to center cells by celltype
# facet_grid facets by sample
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
# Figure 2: display each cell's center distance on the spatial slice
# Each Seurat::Images(object) represents a sample slice
# -----------------------------
p1_slice <- list()
for (i in Seurat::Images(object)) {
        p1_slice[[i]] <- SpatialPlot(
            object,
            features = paste0(center, "_distance")
        )[[1]]
    }

# Combine all slice plots arranged in one row, sharing the legend
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