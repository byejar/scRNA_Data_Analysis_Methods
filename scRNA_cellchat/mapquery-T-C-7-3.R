####all-v2
rm(list=ls())
library(Seurat)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(patchwork)
library(CellChat)
options(stringsAsFactors = FALSE)
data.dir <- '/home/wus/scRNA_for_fig/6-29_CellChat'
#dir.create(data.dir)
setwd(data.dir)


# load(file='/home/wus/c_t_all_12_7_v2.RData')
# Idents(immune.combined.sct)
# immune_combine=immune.combined.sct
# Idents(immune_combine) <- "seurat_clusters"
# immune_combine <- RenameIdents(immune_combine, `0` = "cancer cell", `1` = "cancer cell", `2` = "macrophage", 
#     `3` = "cancer cell", `4` = "neurtophil", `5` = "fibroblast", `6` = "T cell", `7` = "7", `8` = "cancer cell", `9` = "DC", `10` = "B cell", `11` = "11")

# # DimPlot(immune_combine, reduction = "umap", #label = TRUE,
# #     repel = TRUE,pt.size = 1)

# #去掉7 11 
# mac.cells <- subset(immune_combine , idents = c("cancer cell","macrophage","fibroblast","neurtophil", "T cell","DC","B cell"))
# # DimPlot(mac.cells, reduction = "umap", #label = TRUE,
# #     repel = TRUE,pt.size = 1)
# saveRDS(mac.cells,file='scRNA.RData')



#####增加分类
load(file = '/home/wus/for_seurat/8_1/immune_rename_2023-7-4.Rdata')
sce=mac_sce_rename
head(mac_sce_rename)
head(sce)
#mac_sce_re<-subset(x=sce,subset=cell_type==c('M1',"M2","T cell","B cell","cancer_cell","DC",'neutrophil','fib'))


#mac_sce_re<-subset(x=sce,subset=cell_type==c('M1',"M2","T cell","B cell","cancer_cell","DC",'neutrophil','fib'))

#mac_sce_rename$cell_type
#mac_sce_re = sce[,sce@meta.data$seurat_clusters %in% c(0,1,2,3,6,13,14,15)]

DimPlot(mac_sce_rename, reduction = 'umap', group.by = 'cell_type',
        label = TRUE, pt.size = 0.5) + NoLegend()

#rownames(mac_sce_re)
mac_sce_rename$ori.group<-substr(colnames(mac_sce_rename),1,2)
###########抽样
# mac.cells.C=subset(mac.cells,orig.ident %in% c("C1","C2","C3"))
# mac.cells.T=subset(mac.cells,orig.ident %in% c("T1","T2","T3"))

# DimPlot(subset(mac.cells.C, downsample = 10000), reduction = "umap", label = TRUE,
#     repel = TRUE,pt.size = 0.01)

# DimPlot(subset(mac.cells.T, downsample = 10000), reduction = "umap", label = TRUE,
#     repel = TRUE,pt.size = 0.01)


############### cellchat
################ scRNA


#创建目录


# 加载每个数据集的cellchat对象，然后合并在一起
# 用户需要在每个数据集上单独运行 CellChat，然后将不同的 CellChat 对象合并在一起。
#object.list <- SplitObject(mac.cells, split.by = "orig.ident")
#object.list
#cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#cellchat

mac.cells.T1=subset(mac_sce_rename,ori.group %in% c("T1"))
# DimPlot(mac.cells.T1, reduction = "umap", #label = TRUE,
#     repel = TRUE,pt.size = 1)


#cellchat <- createCellChat(data = mac.cells.T1,group.by = "ident")
#将细胞的identity加到matadata上

# head(mac.cells.T1[[]])

# mac.cells.T1 <- StashIdent(mac.cells.T1, save.name = 'cell_type')

# head(mac.cells.T1[[]])

# #####
# mac.cells.T1_nor <- NormalizeData(mac.cells.T1@assays$RNA@counts)
# cellchat <- createCellChat(data = mac.cells.T1_nor,group.by = "ident")

# #normalize的结果在data,但是这里data和counts数值一样，估计是之前没normalize过
# data.input  <- mac.cells.T1@assays$RNA@data
# cellchat <- createCellChat(data = data.input,group.by = "ident")
# identity = data.frame(group =mac.cells.T1$cell_type   , row.names = names(mac.cells.T1$cell_type)) # create a dataframe consisting of the cell labels
# unique(identity$group) # check the cell labels
# cellchat <- createCellChat(data = mac.cells.T1_nor,meta =identity,group.by = "ident")
# cellchat

## sct的归一化方式也支持
Idents(mac.cells.T1)="cell_type"
cellchat <- createCellChat(object = mac.cells.T1@assays$SCT@data,group.by = "ident")
cellchat

#####
#mac.cells.T1[["cell_type"]] <- Idents(object = mac.cells.T1)
identity = data.frame(group =mac.cells.T1$cell_type   , row.names = names(mac.cells.T1$cell_type))

cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
####
# mac.cells.T1[["cell_type"]] <- Idents(object = mac.cells.T1)
# data.input  <- mac.cells.T1@assays$SCT@data
# identity = data.frame(group =mac.cells.T1$cell_type   , row.names = names(mac.cells.T1$cell_type)) # create a dataframe consisting of the cell labels
# unique(identity$group) # check the cell labels
# cellchat <- createCellChat(data = data.input,meta =identity,group.by = "cell_type")
# cellchat

###
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
#showDatabaseCategory(CellChatDB)

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
#showDatabaseCategory(CellChatDB.use)
#用某个
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
#showDatabaseCategory(CellChatDB.use)
# set the used database in the object
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProb(cellchat)  
mycomputeCommunProb <-edit(computeCommunProb)  
environment(mycomputeCommunProb) <- environment(computeCommunProb)
cellchat <- mycomputeCommunProb(cellchat)  

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat@netP$pathways


cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

pathways.show <- c("SPP1") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")


# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object

netAnalysis_contribution(cellchat, signaling = pathways.show)

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
#> Comparing communications on a single object

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)


# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2


library(NMF)
library(ggalluvial)
selectK(cellchat, pattern = "outgoing")
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)

# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
#> Please make sure you have load `library(ggalluvial)` when running this function
# dot plot
netAnalysis_dot(cellchat, pattern = "outgoing")


nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "incoming")
#> Please make sure you have load `library(ggalluvial)` when running this function



cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)


################C1
mac.cells.C1=subset(mac_sce_rename,ori.group %in% c("C1"))