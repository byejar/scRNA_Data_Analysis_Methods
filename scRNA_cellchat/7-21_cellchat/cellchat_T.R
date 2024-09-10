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


load(file = '/home/wus/for_seurat/8_1/immune_rename_2023-7-20.Rdata')

head(sce)

# DimPlot(sce, reduction = 'umap', group.by = 'cell_type',
#         label = TRUE, pt.size = 0.5) + NoLegend()


Idents(sce)="cell_type"
mac.cells.T=subset(sce,orig.ident %in% c("T"))
data.input <-normalizeData(GetAssayData(mac.cells.T, assay = "RNA", slot = "data"), scale.factor = 10000, do.log = TRUE)
labels <- Idents(mac.cells.T)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
cellchat.T <- createCellChat(object = data.input, meta = meta, group.by = "group")
cellchat.T <- addMeta(cellchat.T, meta = meta, meta.name = "labels")
cellchat.T <- setIdent(cellchat.T, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat.T@idents) # show factor levels of the cell labels


cellchat.T@DB <- CellChatDB.mouse 

# subset the expression data of signaling genes for saving computation cost
cellchat.T <- subsetData(cellchat.T) # This step is necessary even if using the whole database

future::availableCores() #12 #查看几个核可用
future::nbrOfWorkers() #4 当前可用的核有多少个

future::plan("multisession", workers = 16) # do parallel
cellchat.T <- identifyOverExpressedGenes(cellchat.T)
cellchat.T <- identifyOverExpressedInteractions(cellchat.T)


data.dir <- '/home/wus/scRNA_for_fig/7-20_CellChat'
dir.create(data.dir)
setwd(data.dir)


# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat.T <- projectData(cellchat.T, PPI.human)
cellchat.T <- computeCommunProb(cellchat.T,raw.use = FALSE)
cellchat.T <- filterCommunication(cellchat.T, min.cells = 10)
df.net <- subsetCommunication(cellchat.T)
#write.csv(df.net,"net_t.csv")

cellchat.T <- computeCommunProbPathway(cellchat.T)
df.netp <- subsetCommunication(cellchat.T,slot.name="netP")
#write.csv(df.netp,"netp_t.csv")

cellchat.T <- aggregateNet(cellchat.T)
groupSize <- as.numeric(table(cellchat.T@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat.T@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.T@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


# cellchat.T.mac_t <- aggregateNet(cellchat.T,sources.use=c('M1','M2'), targets.use=c('T cell'),remove.isolate = FALSE)
# groupSize <- as.numeric(table(cellchat.T.mac_t@idents))
# par(mfrow = c(1,2), xpd=TRUE)
# netVisual_circle(cellchat.T.mac_t@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
# netVisual_circle(cellchat.T.mac_t@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


pathways.show <- c("SPP1") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
#vertex.receiver = c('T cell')
netVisual_aggregate(cellchat.T, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout='hierarchy')
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat.T, signaling = pathways.show, layout = "circle")
# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat.T, signaling = pathways.show, layout = "chord")
# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat.T, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object
netAnalysis_contribution(cellchat.T, signaling = pathways.show)

######### mac t
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat.T,sources.use=c('M1','M2'), targets.use='T cell',remove.isolate = FALSE)
#> Comparing communications on a single object
##巨噬细胞收到的
netVisual_bubble(cellchat.T, targets.use=c('M1','M2'),remove.isolate = FALSE)
##巨噬细胞发出的
netVisual_bubble(cellchat.T, sources.use=c('M1','M2'),remove.isolate = FALSE)

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
#
# show all the interactions received by T cell
netVisual_chord_gene(cellchat.T, sources.use=c('M1','M2'), targets.use='T cell', legend.pos.x = 15)

# show all the interactions received by 'M1','M2'
netVisual_chord_gene(cellchat.T,  targets.use=c('M1','M2'), legend.pos.x = 15)

plotGeneExpression(cellchat.T, signaling = "CXCL")

plotGeneExpression(cellchat.T, signaling = "SPP1")

# Compute the network centrality scores
cellchat.T <- netAnalysis_computeCentrality(cellchat.T, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat.T, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat.T)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat.T, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.T, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.T, pattern = "incoming")
ht1 + ht2


library(NMF)
library(ggalluvial)
selectK(cellchat.T, pattern = "outgoing")
nPatterns = 4
cellchat.T <- identifyCommunicationPatterns(cellchat.T, pattern = "outgoing", k = nPatterns)

# river plot
netAnalysis_river(cellchat.T, pattern = "outgoing")

selectK(cellchat.T, pattern = "incoming")
nPatterns = 4
cellchat.T <- identifyCommunicationPatterns(cellchat.T, pattern = "incoming", k = nPatterns)
# river plot
netAnalysis_river(cellchat.T, pattern = "incoming")
#> Please make sure you have load `library(ggalluvial)` when running this function



cellchat.T <- computeNetSimilarity(cellchat.T, type = "functional")
cellchat.T <- netEmbedding(cellchat.T, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat.T <- netClustering(cellchat.T, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat.T, type = "functional", label.size = 3.5)

saveRDS(cellchat.T,'cellchat.T.Rds')
