################C1
mac.cells.C1=subset(mac_sce_rename,ori.group %in% c("C1"))
Idents(mac.cells.C1)="cell_type"
cellchat <- createCellChat(object = mac.cells.C1@assays$SCT@data,group.by = "ident")
cellchat

#####
#mac.cells.T1[["cell_type"]] <- Idents(object = mac.cells.T1)
identity = data.frame(group =mac.cells.C1$cell_type   , row.names = names(mac.cells.C1$cell_type))

cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents))
groupSize # number of cells in each cell group
table(Idents(mac.cells.C1))
table(cellchat@idents)
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



##############################################################################
mac.cells.C1=subset(mac_sce_rename,ori.group %in% c("C1"))
Idents(mac.cells.C1)="cell_type"
cellchat.C1 <- createCellChat(object = mac.cells.C1@assays$SCT@data,group.by = "ident")

mac.cells.T1=subset(mac_sce_rename,ori.group %in% c("T1"))
Idents(mac.cells.T1)="cell_type"
cellchat.T1 <- createCellChat(object = mac.cells.T1@assays$SCT@data,group.by = "ident")


object.list <- list(C1 = cellchat.C1, T1 = cellchat.T1)
cellchat <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = TRUE)
#> Merge the following slots: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
cellchat

#############
Idents(mac_sce_rename)="cell_type"
#data.input <- GetAssayData(mac_sce_rename, assay = "RNA", slot = "data") # normalized data matrix
mac.cells.C1=subset(mac_sce_rename,ori.group %in% c("C1"))
data.input <-normalizeData(GetAssayData(mac.cells.C1, assay = "RNA", slot = "data"), scale.factor = 10000, do.log = TRUE)
labels <- Idents(mac.cells.C1)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
cellchat.C1 <- createCellChat(object = data.input, meta = meta, group.by = "group")
cellchat.C1 <- addMeta(cellchat, meta = meta, meta.name = "labels")
cellchat.C1 <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat.C1@idents) # show factor levels of the cell labels


cellchat.C1@DB <- CellChatDB.mouse 

# subset the expression data of signaling genes for saving computation cost
cellchat.C1 <- subsetData(cellchat.C1) # This step is necessary even if using the whole database

future::availableCores() #12 #查看几个核可用
future::nbrOfWorkers() #4 当前可用的核有多少个

future::plan("multisession", workers = 4) # do parallel
cellchat.C1 <- identifyOverExpressedGenes(cellchat.C1)
cellchat.C1 <- identifyOverExpressedInteractions(cellchat.C1)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat.C1 <- projectData(cellchat.C1, PPI.human)
cellchat.C1 <- computeCommunProb(cellchat.C1,raw.use = FALSE)
cellchat.C1 <- filterCommunication(cellchat.C1, min.cells = 10)
df.net <- subsetCommunication(cellchat.C1)
write.csv(df.net,"net_c1.csv")

cellchat.C1 <- computeCommunProbPathway(cellchat.C1)
df.netp <- subsetCommunication(cellchat.C1,slot.name="netP")
write.csv(df.netp,"netp_c1.csv")

cellchat.C1 <- aggregateNet(cellchat.C1)
groupSize <- as.numeric(table(cellchat.C1@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat.C1@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.C1@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


cellchat.C1.mac_t <- aggregateNet(cellchat.C1,sources.use=c('M1','M2'), targets.use=c('T cell'),remove.isolate = FALSE)
groupSize <- as.numeric(table(cellchat.C1.mac_t@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat.C1.mac_t@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.C1.mac_t@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


pathways.show <- c("CXCL") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
#vertex.receiver = c('T cell')
netVisual_aggregate(cellchat.C1, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout='hierarchy')
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat.C1, signaling = pathways.show, layout = "circle")
# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat.C1, signaling = pathways.show, layout = "chord")
# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat.C1, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object
netAnalysis_contribution(cellchat.C1, signaling = pathways.show)

######### mac t
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat.C1,sources.use=c('M1','M2'), targets.use='T cell',remove.isolate = FALSE)
#> Comparing communications on a single object
##巨噬细胞收到的
netVisual_bubble(cellchat.C1, targets.use=c('M1','M2'),remove.isolate = FALSE)
##巨噬细胞发出的
netVisual_bubble(cellchat.C1, sources.use=c('M1','M2'),remove.isolate = FALSE)

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
#
# show all the interactions received by T cell
netVisual_chord_gene(cellchat.C1, sources.use=c('M1','M2'), targets.use='T cell', legend.pos.x = 15)

# show all the interactions received by 'M1','M2'
netVisual_chord_gene(cellchat.C1,  targets.use=c('M1','M2'), legend.pos.x = 15)

plotGeneExpression(cellchat.C1, signaling = "CXCL")



# Compute the network centrality scores
cellchat.C1 <- netAnalysis_computeCentrality(cellchat.C1, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat.C1, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat.C1)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat.C1, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.C1, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.C1, pattern = "incoming")
ht1 + ht2


library(NMF)
library(ggalluvial)
selectK(cellchat.C1, pattern = "outgoing")
nPatterns = 4
cellchat.C1 <- identifyCommunicationPatterns(cellchat.C1, pattern = "outgoing", k = nPatterns)

# river plot
netAnalysis_river(cellchat.C1, pattern = "outgoing")
#> Please make sure you have load `library(ggalluvial)` when running this function

cellchat.C1 <- computeNetSimilarity(cellchat.C1, type = "functional")
cellchat.C1 <- netEmbedding(cellchat.C1, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat.C1 <- netClustering(cellchat.C1, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat.C1, type = "functional", label.size = 3.5)

saveRDS(cellchat.C1,'cellchat_C1.Rds')