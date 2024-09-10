####all-v2
rm(list=ls())
library(Seurat)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(patchwork)

C1 <- Read10X(data.dir = "/home/wus/scRUN/sampleC1/outs/filtered_feature_bc_matrix")
C1_obj <- CreateSeuratObject(counts = C1, project = "C1", min.cells = 3, min.features = 200)
C1_obj[["percent.mt"]] <- PercentageFeatureSet(C1_obj, pattern = "^mt-")
C1_obj <- subset(C1_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 8200 & percent.mt < 5)
C1_obj


C2 <- Read10X(data.dir = "/home/wus/scRUN/sampleC2/outs/filtered_feature_bc_matrix")
C2_obj <- CreateSeuratObject(counts = C2, project = "C2", min.cells = 3, min.features = 200)
C2_obj[["percent.mt"]] <- PercentageFeatureSet(C2_obj, pattern = "^mt-")
C2_obj <- subset(C2_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 8200 & percent.mt < 5)
C2_obj

C3 <- Read10X(data.dir = "/home/wus/scRUN/sampleC3/outs/filtered_feature_bc_matrix")
C3_obj <- CreateSeuratObject(counts = C3, project = "C3", min.cells = 3, min.features = 200)
C3_obj[["percent.mt"]] <- PercentageFeatureSet(C3_obj, pattern = "^mt-")
C3_obj <- subset(C3_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 8200 & percent.mt < 5)
C3_obj


T1 <- Read10X(data.dir = "/home/wus/scRUN/sampleT1/outs/filtered_feature_bc_matrix")
T1_obj <- CreateSeuratObject(counts = T1, project = "T1", min.cells = 3, min.features = 200)
T1_obj[["percent.mt"]] <- PercentageFeatureSet(T1_obj, pattern = "^mt-")
T1_obj <- subset(T1_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 8200 & percent.mt < 5)
T1_obj

T2 <- Read10X(data.dir = "/home/wus/scRUN/sampleT2/outs/filtered_feature_bc_matrix")
T2_obj <- CreateSeuratObject(counts = T2, project = "T2", min.cells = 3, min.features = 200)
T2_obj[["percent.mt"]] <- PercentageFeatureSet(T2_obj, pattern = "^mt-")
T2_obj <- subset(T2_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 8200 & percent.mt < 5)
T2_obj

T3 <- Read10X(data.dir = "/home/wus/scRUN/sampleT3/outs/filtered_feature_bc_matrix")
T3_obj <- CreateSeuratObject(counts = T3, project = "T3", min.cells = 3, min.features = 200)
T3_obj[["percent.mt"]] <- PercentageFeatureSet(T3_obj, pattern = "^mt-")
T3_obj <- subset(T3_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 8200 & percent.mt < 5)
T3_obj


sc.combined <- merge(C1_obj, y = c(C2_obj,C3_obj,T1_obj,T2_obj,T3_obj), add.cell.ids = c("C1","C2","C3",'T1','T2','T3'), project = "scRNA")
head(sc.combined@meta.data)
pancreas.list <- SplitObject(sc.combined, split.by = "orig.ident")

C_group <- subset(sc.combined,orig.ident %in% c("C1","C2","C3"))

T_group <- subset(sc.combined,orig.ident %in% c('T1','T2','T3'))


ctrl <- SCTransform(C_group, vst.flavor = "v2", verbose = FALSE) %>%
    RunPCA(npcs = 30, verbose = FALSE)# %>%

stim <- SCTransform(T_group, vst.flavor = "v2", verbose = FALSE) %>%
    RunPCA(npcs = 30, verbose = FALSE)


ifnb.list <- list(ctrl = ctrl, stim = stim)
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)

#为了集成这两个数据集，我们使用函数，该函数将 Seurat 对象列表作为输入，并使用这些锚点将两个数据集与 集成在一起。
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",
    anchor.features = features)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
#执行集成分析
immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30, verbose = FALSE)
immune.combined.sct <- FindNeighbors(immune.combined.sct, reduction = "pca", dims = 1:30)
immune.combined.sct <- FindClusters(immune.combined.sct, resolution = 0.3)
#并排可视化这两个条件，使用参数来显示按聚类着色的每个条件。split.by
DimPlot(immune.combined.sct, reduction = "umap", split.by = "orig.ident")

save(immune.anchors,immune.combined.sct,file='/home/wus/c_t_all_12_7_v2.RData')

##以上，比对完成

load(file='/home/wus/c_t_all_12_7_v2.RData')
Idents(immune.combined.sct)
immune_combine=immune.combined.sct
Idents(immune_combine) <- "seurat_clusters"
immune_combine <- RenameIdents(immune_combine, `0` = "cancer cell", `1` = "cancer cell", `2` = "macrophage", 
    `3` = "cancer cell", `4` = "neurtophil", `5` = "fibroblast", `6` = "T cell", `7` = "7", `8` = "cancer cell", `9` = "DC", `10` = "B cell", `11` = "11")

DimPlot(immune_combine, reduction = "umap", #label = TRUE,
    repel = TRUE,pt.size = 1)

#去掉7 11 
mac.cells <- subset(immune_combine , idents = c("cancer cell","macrophage","fibroblast","neurtophil", "T cell","DC","B cell"))
# DimPlot(mac.cells, reduction = "umap", #label = TRUE,
#     repel = TRUE,pt.size = 1)

###########抽样
# mac.cells.C=subset(mac.cells,orig.ident %in% c("C1","C2","C3"))
# mac.cells.T=subset(mac.cells,orig.ident %in% c("T1","T2","T3"))

# DimPlot(subset(mac.cells.C, downsample = 10000), reduction = "umap", label = TRUE,
#     repel = TRUE,pt.size = 0.01)

# DimPlot(subset(mac.cells.T, downsample = 10000), reduction = "umap", label = TRUE,
#     repel = TRUE,pt.size = 0.01)


############### cellchat
################ scRNA

library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
#创建目录
data.dir <- '/home/wus/scRNA_for_fig/6-29_CellChat'
dir.create(data.dir)
setwd(data.dir)

# 加载每个数据集的cellchat对象，然后合并在一起
# 用户需要在每个数据集上单独运行 CellChat，然后将不同的 CellChat 对象合并在一起。
#object.list <- SplitObject(mac.cells, split.by = "orig.ident")
#object.list
#cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#cellchat

mac.cells.T1=subset(mac.cells,orig.ident %in% c("T1"))
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
cellchat <- createCellChat(object = mac.cells.T1@assays$SCT@data,group.by = "ident")
cellchat

#####
mac.cells.T1[["cell_type"]] <- Idents(object = mac.cells.T1)
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
showDatabaseCategory(CellChatDB)

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
showDatabaseCategory(CellChatDB.use)
#用某个
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
showDatabaseCategory(CellChatDB.use)
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

pathways.show <- c("CXCL") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
pathways.show <- c("SPP1")


