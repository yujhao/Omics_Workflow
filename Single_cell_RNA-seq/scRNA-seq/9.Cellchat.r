#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# ==============================================================================
# CellChat 细胞通讯分析脚本
# 分析细胞间的通讯相互作用
# ==============================================================================

library(CellChat)

# ========================== 配置参数 ==========================

OUTPUT_DIR <- "9.cellchat"
IDENT_COLUMN <- "new_celltype"
MIN_CELLS <- 10
PATHWAYS_SHOW <- c("MIF")

# ========================== 初始化 ==========================

dir.create(OUTPUT_DIR, showWarnings = FALSE)

# ========================== 全局 CellChat 分析 ==========================

# 设置细胞类型
scRNA_mnn <- SetIdent(scRNA_mnn, value = IDENT_COLUMN)

# 创建 CellChat 对象
cellchat <- createCellChat(object = scRNA_mnn, group.by = IDENT_COLUMN)

# 设置数据库
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use

# 数据处理
cellchat <- CellChat::subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = MIN_CELLS)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# ========================== 全局可视化 ==========================

groupSize <- as.numeric(table(cellchat@idents))

# 互作数量图
pdf(file.path(OUTPUT_DIR, "cellchat_circle_count.pdf"), width = 8, height = 4)
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(
    cellchat@net$count,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = "Number of interactions"
)
dev.off()

# 互作强度图
pdf(file.path(OUTPUT_DIR, "cellchat_circle_weight.pdf"), width = 8, height = 4)
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(
    cellchat@net$weight,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = "Interaction weights/strength"
)
dev.off()

# 聚合图 - Circle 布局
pdf(file.path(OUTPUT_DIR, "cellchat_circle_aggregate.pdf"), width = 8, height = 4)
par(mfrow = c(1, 1))
netVisual_aggregate(cellchat, signaling = PATHWAYS_SHOW, layout = "circle")
dev.off()

# 聚合图 - Hierarchy 布局
pdf(file.path(OUTPUT_DIR, "cellchat_circle_aggregate_hierarchy.pdf"), width = 8, height = 4)
vertex.receiver <- seq(1, 3)
netVisual_aggregate(cellchat, signaling = PATHWAYS_SHOW, vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

# 聚合图 - Chord 布局
pdf(file.path(OUTPUT_DIR, "cellchat_circle_aggregate_chord.pdf"), width = 8, height = 4)
netVisual_aggregate(cellchat, signaling = PATHWAYS_SHOW, layout = "chord")
dev.off()

# 热图 - count
pdf(file.path(OUTPUT_DIR, "cellchat_heatmap_count.pdf"), width = 8, height = 4)
par(mfrow = c(1, 1))
netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Reds")
dev.off()

# 热图 - weight
pdf(file.path(OUTPUT_DIR, "cellchat_heatmap_weight.pdf"), width = 8, height = 4)
par(mfrow = c(1, 1))
netVisual_heatmap(cellchat, measure = "weight", color.heatmap = "Reds")
dev.off()

# 热图 - pathway
pdf(file.path(OUTPUT_DIR, "cellchat_heatmap_pathways.pdf"), width = 8, height = 4)
par(mfrow = c(1, 1))
netVisual_heatmap(cellchat, signaling = PATHWAYS_SHOW, color.heatmap = "Reds")
dev.off()

# 信号角色分析
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

pdf(file.path(OUTPUT_DIR, "cellchat_signalingRole_network.pdf"), width = 8, height = 4)
netAnalysis_signalingRole_network(cellchat, signaling = PATHWAYS_SHOW)
dev.off()

pdf(file.path(OUTPUT_DIR, "cellchat_signalingRole_scatter.pdf"), width = 8, height = 4)
netAnalysis_signalingRole_scatter(cellchat)
dev.off()

pdf(file.path(OUTPUT_DIR, "cellchat_signalingRole_heatmap.pdf"), width = 8, height = 4)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2
dev.off()

# Bubble plot
pdf(file.path(OUTPUT_DIR, "cellchat_bubble.pdf"), width = 8, height = 4)
netVisual_bubble(cellchat, remove.isolate = FALSE, return.data = FALSE)
dev.off()

# Chord gene plot
pdf(file.path(OUTPUT_DIR, "cellchat_chord_gene.pdf"), width = 8, height = 4)
netVisual_chord_gene(cellchat, lab.cex = 0.8, legend.pos.y = 50, legend.pos.x = 50)
dev.off()

# ========================== 分组比较分析 ==========================

# 筛选细胞类型
scRNA_mnn <- subset(scRNA_mnn, new_celltype %in% c("EC1", "EC3", "EC4", "FC1", "FC2", "FC3"))
scRNA_mnn@meta.data$celltype <- droplevels(scRNA_mnn@meta.data$celltype)

# 按样本拆分
sp <- SplitObject(scRNA_mnn, split.by = "orig.ident")

# 创建 CellChat 对象列表
cellchat_list <- list()
for (i in names(sp)) {
    cellchat <- createCellChat(object = sp[[i]], group.by = IDENT_COLUMN)
    cellchat@idents <- droplevels(cellchat@idents)
    cellchat@DB <- CellChatDB.human
    cellchat <- CellChat::subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- computeCommunProb(cellchat)
    cellchat <- filterCommunication(cellchat, min.cells = MIN_CELLS)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    cellchat_list[[i]] <- cellchat
}

# 合并 CellChat 对象
cellchat <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))

# 比较互作
pdf(file.path(OUTPUT_DIR, "cellchat_compare_interaction.pdf"), width = 8, height = 4)
gg1 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2), measure = "weight")
gg1 + gg2
dev.off()

# 差异互作
pdf(file.path(OUTPUT_DIR, "cellchat_compare_interaction_diff.pdf"), width = 8, height = 4)
netVisual_diffInteraction(cellchat, weight.scale = TRUE, comparison = c(2, 1))
netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")
dev.off()

# 比较热图
pdf(file.path(OUTPUT_DIR, "cellchat_compare_interaction_heatmap.pdf"), width = 8, height = 4)
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2
dev.off()
