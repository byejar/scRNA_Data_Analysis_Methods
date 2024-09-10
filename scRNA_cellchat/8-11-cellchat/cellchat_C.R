#### scRNA
rm(list=ls())
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(CellChat)
options(stringsAsFactors = FALSE)


load(file = '/home/wus/for_seurat/8_1/immune_rename_2023-7-20.Rdata')

#### 合并
Idents(sce)<-'cell_type'
levels(Idents(sce))
p1 = DimPlot(sce, reduction = "umap", pt.size = 0.5, label = TRUE) 

new.cluster.ids <- c("DC","B cell","cancer_cell","Mac","T_cell",
                     "Mac","undifine","fib","neutrophil","NK")
names(new.cluster.ids) <- levels(sce)
sce <- RenameIdents(sce, new.cluster.ids)
p2 = DimPlot(sce, reduction = "umap", label = TRUE,pt.size = 0.5) 

p1+p2

head(sce)

# Idents(sce)="cell_type"
mac.cells.C=subset(sce,orig.ident %in% c("C"))
data.input <-normalizeData(GetAssayData(mac.cells.C, assay = "RNA", slot = "data"), scale.factor = 10000, do.log = TRUE)
labels <- Idents(mac.cells.C)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
cellchat.C <- createCellChat(object = data.input, meta = meta, group.by = "group")
cellchat.C <- addMeta(cellchat.C, meta = meta, meta.name = "labels")
cellchat.C <- setIdent(cellchat.C, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat.C@idents) # show factor levels of the cell labels


cellchat.C@DB <- CellChatDB.mouse 

# subset the expression data of signaling genes for saving computation cost
cellchat.C <- subsetData(cellchat.C) # This step is necessary even if using the whole database

future::availableCores() #12 #查看几个核可用
future::nbrOfWorkers() #4 当前可用的核有多少个

future::plan("multisession", workers = 16) # do parallel
cellchat.C <- identifyOverExpressedGenes(cellchat.C)
cellchat.C <- identifyOverExpressedInteractions(cellchat.C)


data.dir <- '/home/wus/scRNA_for_fig/8-11_CellChat'
dir.create(data.dir)
setwd(data.dir)


# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat.C <- projectData(cellchat.C, PPI.human)
cellchat.C <- computeCommunProb(cellchat.C,raw.use = FALSE)
cellchat.C <- filterCommunication(cellchat.C, min.cells = 10)
df.net <- subsetCommunication(cellchat.C)
#write.csv(df.net,"net_t.csv")

cellchat.C <- computeCommunProbPathway(cellchat.C)
df.netp <- subsetCommunication(cellchat.C,slot.name="netP")
#write.csv(df.netp,"netp_t.csv")

cellchat.C <- aggregateNet(cellchat.C)
groupSize <- as.numeric(table(cellchat.C@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat.C@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.C@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


# cellchat.C.mac_t <- aggregateNet(cellchat.C,sources.use=c('M1','M2'), targets.use=c('T cell'),remove.isolate = FALSE)
# groupSize <- as.numeric(table(cellchat.C.mac_t@idents))
# par(mfrow = c(1,2), xpd=TRUE)
# netVisual_circle(cellchat.C.mac_t@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
# netVisual_circle(cellchat.C.mac_t@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


# pathways.show <- c("SPP1") 
# # Hierarchy plot
# # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
# vertex.receiver = seq(1,4) # a numeric vector. 
# #vertex.receiver = c('T cell')
# netVisual_aggregate(cellchat.C, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout='hierarchy')
# # Circle plot
# par(mfrow=c(1,1))
# netVisual_aggregate(cellchat.C, signaling = pathways.show, layout = "circle")
# # Chord diagram
# par(mfrow=c(1,1))
# netVisual_aggregate(cellchat.C, signaling = pathways.show, layout = "chord")
# # Heatmap
# par(mfrow=c(1,1))
# netVisual_heatmap(cellchat.C, signaling = pathways.show, color.heatmap = "Reds")
# #> Do heatmap based on a single object
# netAnalysis_contribution(cellchat.C, signaling = pathways.show)

# ######### mac t
# # show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# netVisual_bubble(cellchat.C,sources.use=c('M1','M2'), targets.use='T cell',remove.isolate = FALSE)
# #> Comparing communications on a single object
# ##巨噬细胞收到的
# netVisual_bubble(cellchat.C, targets.use=c('M1','M2'),remove.isolate = FALSE)
# ##巨噬细胞发出的
# netVisual_bubble(cellchat.C, sources.use=c('M1','M2'),remove.isolate = FALSE)

# # show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# # show all the interactions sending from Inflam.FIB
# #
# # show all the interactions received by T cell
# netVisual_chord_gene(cellchat.C, sources.use=c('M1','M2'), targets.use='T cell', legend.pos.x = 15)

# # show all the interactions received by 'M1','M2'
# netVisual_chord_gene(cellchat.C,  targets.use=c('M1','M2'), legend.pos.x = 15)

# plotGeneExpression(cellchat.C, signaling = "CXCL")

# plotGeneExpression(cellchat.C, signaling = "SPP1")

# Compute the network centrality scores
cellchat.C <- netAnalysis_computeCentrality(cellchat.C, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
# netAnalysis_signalingRole_network(cellchat.C, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# gg1 <- netAnalysis_signalingRole_scatter(cellchat.C)
# #> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# # Signaling role analysis on the cell-cell communication networks of interest
# gg2 <- netAnalysis_signalingRole_scatter(cellchat.C, signaling = c("CXCL", "CCL"))
# #> Signaling role analysis on the cell-cell communication network from user's input
# gg1 + gg2

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# ht1 <- netAnalysis_signalingRole_heatmap(cellchat.C, pattern = "outgoing")
# ht2 <- netAnalysis_signalingRole_heatmap(cellchat.C, pattern = "incoming")
# ht1 + ht2


library(NMF)
library(ggalluvial)
selectK(cellchat.C, pattern = "outgoing")
nPatterns = 4
cellchat.C <- identifyCommunicationPatterns(cellchat.C, pattern = "outgoing", k = nPatterns)

# river plot
netAnalysis_river(cellchat.C, pattern = "outgoing")

selectK(cellchat.C, pattern = "incoming")
nPatterns = 4
cellchat.C <- identifyCommunicationPatterns(cellchat.C, pattern = "incoming", k = nPatterns)
# river plot
netAnalysis_river(cellchat.C, pattern = "incoming")
#> Please make sure you have load `library(ggalluvial)` when running this function



cellchat.C <- computeNetSimilarity(cellchat.C, type = "functional")
cellchat.C <- netEmbedding(cellchat.C, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
### K=1，不用多核，接解决了这里报错的问题
cellchat.C <- netClustering(cellchat.C, type = "functional",k=1)
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat.C, type = "functional", label.size = 3.5)

saveRDS(cellchat.C,'cellchat.C.Rds')
