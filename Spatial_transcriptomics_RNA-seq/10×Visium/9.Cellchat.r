options(stringsAsFactors = FALSE)
library(CellChat)
library(Seurat)
library(tidyverse)
library(viridis)
library(RColorBrewer)

library(CellChat)
library(jsonlite)
library(future)

scale.factors0 = jsonlite::fromJSON(
  file.path("E:/bioinformation/sc_ST/data/Brain_35707680/spatial", "scalefactors_json.json")
)

scale.factors = list(
  spot.diameter = 65,
  spot = scale.factors0$spot_diameter_fullres,
  fiducial = scale.factors0$fiducial_diameter_fullres,
  hires = scale.factors0$tissue_hires_scalef,
  lowres = scale.factors0$tissue_lowres_scalef
)

CellChatDB.use = CellChatDB.human

cellchat_list = list()
data_ob_list = SplitObject(data_ob, split.by = "sampleid")
for(i in seq(1,length(data_ob_list))){

  obj = data_ob_list[[i]]

  data.input = Seurat::GetAssayData(obj, slot = "data", assay = "SCT")

  meta = data.frame(
    labels = Idents(obj),
    row.names = names(Idents(obj))
  )

  spatial.locs = Seurat::GetTissueCoordinates(
    obj,
    scale = NULL,
    cols = c("imagerow", "imagecol")
  )

  cellchat = createCellChat(
    object = data.input,
    meta = meta,
    group.by = "labels",
    datatype = "spatial",
    coordinates = spatial.locs,
    scale.factors = scale.factors
  )

  cellchat@DB = CellChatDB.use
  cellchat = subsetData(cellchat)

  future::plan("multisession", workers = 1)

  cellchat = identifyOverExpressedGenes(cellchat)
  cellchat = identifyOverExpressedInteractions(cellchat)
  cellchat = computeCommunProb(
    cellchat,
    type = "truncatedMean",
    trim = 0.1,
    distance.use = TRUE,
    scale.distance = 0.01
  )
  cellchat = filterCommunication(cellchat, min.cells = 10)
  cellchat = computeCommunProbPathway(cellchat)
  cellchat = aggregateNet(cellchat)

  cellchat_list[[i]] = cellchat
}

names(cellchat_list) = names(data_ob_list)
saveRDS(cellchat_list, file = "cellchat_list.rds")

# compare
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

count.max = getMaxWeight(cellchat_list, attribute = c("idents", "count"))
weight.max = getMaxWeight(cellchat_list, attribute = c("idents", "weight"))

for(i in 1:length(cellchat_list)){

  x = cellchat_list[[i]]
  nm = names(cellchat_list)[i]
  groupSize = as.numeric(table(x@idents))
  sub_col = color_use[names(groupSize)]
  pdf(file.path(output_dir, paste0(nm, "_interaction_number_network.pdf")))
  par(mfrow = c(1,1), xpd = TRUE)
  netVisual_circle(
    x@net$count,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    label.edge = FALSE,
    arrow.size = 0.5,
    edge.weight.max = count.max[2],
    edge.width.max = 12,
    title.name = paste0("Number of interactions - ", nm),
    color.use = sub_col
  )
  dev.off()

  pdf(file.path(output_dir, paste0(nm, "_interaction_number_heatmap.pdf")))
  par(mfrow = c(1,1), xpd = TRUE)
  ht = netVisual_heatmap(
    x,
    measure = "count",
    color.heatmap = "Reds",
    color.use = sub_col,
    title.name = paste0("Number of interactions - ", nm)
  )
  ComplexHeatmap::draw(ht, ht_gap = unit(0.5, "cm"))
  dev.off()

  num_df = cbind(cell = rownames(x@net$count), x@net$count)
  write.table(
    num_df,
    file.path(output_dir, paste0(nm, "_interaction_number_network.xls")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  pdf(file.path(output_dir, paste0(nm, "_interaction_strength_network.pdf")))
  par(mfrow = c(1,1), xpd = TRUE)
  netVisual_circle(
    x@net$weight,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    label.edge = FALSE,
    arrow.size = 0.5,
    edge.weight.max = weight.max[2],
    edge.width.max = 12,
    title.name = paste0("Interaction strength - ", nm),
    color.use = sub_col
  )
  dev.off()

  pdf(file.path(output_dir, paste0(nm, "_interaction_strength_heatmap.pdf")))
  par(mfrow = c(1,1), xpd = TRUE)
  ht = netVisual_heatmap(
    x,
    measure = "weight",
    color.heatmap = "Reds",
    color.use = sub_col,
    title.name = paste0("Interaction strength - ", nm)
  )
  ComplexHeatmap::draw(ht, ht_gap = unit(0.5, "cm"))
  dev.off()

  wt_df = cbind(cell = rownames(x@net$weight), x@net$weight)
  colnames(wt_df)[1] = "ligand_cell/receptor_cell"
  write.table(
    wt_df,
    file.path(output_dir, paste0(nm, "_interaction_strength_network.xls")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

cellchat_merge = mergeCellChat(cellchat_list, add.names = names(cellchat_list))
gg = rankNet(
  cellchat_merge,
  mode = "comparison",
  stacked = FALSE,
  comparison = c(1:length(cellchat_list)),
  color.use = color_use_group2
)

ggsave(
  file.path(output_dir, "information_flow_of_each_signaling_pathway.pdf"),
  plot = gg,
  height = 8,
  width = 5,
  bg = "white"
)

write.table(
  gg$data,
  file.path(output_dir, "information_flow_of_each_signaling_pathway.xls"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)