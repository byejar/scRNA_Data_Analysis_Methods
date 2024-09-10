################T1
load(file = '/home/wus/for_seurat/8_1/immune_rename_2023-7-4.Rdata')
sce=mac_sce_rename
head(mac_sce_rename)
head(sce)

mac_sce_rename$ori.group<-substr(colnames(mac_sce_rename),1,2)



##############################################################################


#############
Idents(mac_sce_rename)="cell_type"
#data.input <- GetAssayData(mac_sce_rename, assay = "RNA", slot = "data") # normalized data matrix
mac.cells.T1=subset(mac_sce_rename,ori.group %in% c("T1"))
data.input <-normalizeData(GetAssayData(mac.cells.T1, assay = "RNA", slot = "data"), scale.factor = 10000, do.log = TRUE)
labels <- Idents(mac.cells.T1)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
cellchat.T1 <- createCellChat(object = data.input, meta = meta, group.by = "group")
cellchat.T1 <- addMeta(cellchat.T1, meta = meta, meta.name = "labels")
cellchat.T1 <- setIdent(cellchat.T1, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat.T1@idents) # show factor levels of the cell labels


cellchat.T1@DB <- CellChatDB.mouse 

# subset the expression data of signaling genes for saving computation cost
cellchat.T1 <- subsetData(cellchat.T1) # This step is necessary even if using the whole database

future::availableCores() #12 #查看几个核可用
future::nbrOfWorkers() #4 当前可用的核有多少个

future::plan("multisession", workers = 16) # do parallel
cellchat.T1 <- identifyOverExpressedGenes(cellchat.T1)
cellchat.T1 <- identifyOverExpressedInteractions(cellchat.T1)


data.dir <- '/home/wus/scRNA_for_fig/6-29_CellChat'
#dir.create(data.dir)
setwd(data.dir)


# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
cellchat.T1 <- projectData(cellchat.T1, PPI.human)
cellchat.T1 <- computeCommunProb(cellchat.T1,raw.use = FALSE)
cellchat.T1 <- filterCommunication(cellchat.T1, min.cells = 10)
df.net <- subsetCommunication(cellchat.T1)
write.csv(df.net,"net_t1.csv")

cellchat.T1 <- computeCommunProbPathway(cellchat.T1)
df.netp <- subsetCommunication(cellchat.T1,slot.name="netP")
write.csv(df.netp,"netp_t1.csv")

cellchat.T1 <- aggregateNet(cellchat.T1)
groupSize <- as.numeric(table(cellchat.T1@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat.T1@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.T1@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


cellchat.T1.mac_t <- aggregateNet(cellchat.T1,sources.use=c('M1','M2'), targets.use=c('T cell'),remove.isolate = FALSE)
groupSize <- as.numeric(table(cellchat.T1.mac_t@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat.T1.mac_t@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.T1.mac_t@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


pathways.show <- c("CXCL") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
#vertex.receiver = c('T cell')
netVisual_aggregate(cellchat.T1, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout='hierarchy')
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat.T1, signaling = pathways.show, layout = "circle")
# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat.T1, signaling = pathways.show, layout = "chord")
# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat.T1, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object
netAnalysis_contribution(cellchat.T1, signaling = pathways.show)

######### mac t
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat.T1,sources.use=c('M1','M2'), targets.use='T cell',remove.isolate = FALSE)
#> Comparing communications on a single object
##巨噬细胞收到的
netVisual_bubble(cellchat.T1, targets.use=c('M1','M2'),remove.isolate = FALSE)
##巨噬细胞发出的
netVisual_bubble(cellchat.T1, sources.use=c('M1','M2'),remove.isolate = FALSE)

# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
#
# show all the interactions received by T cell
netVisual_chord_gene(cellchat.T1, sources.use=c('M1','M2'), targets.use='T cell', legend.pos.x = 15)

# show all the interactions received by 'M1','M2'
netVisual_chord_gene(cellchat.T1,  targets.use=c('M1','M2'), legend.pos.x = 15)

plotGeneExpression(cellchat.T1, signaling = "CXCL")



# Compute the network centrality scores
cellchat.T1 <- netAnalysis_computeCentrality(cellchat.T1, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat.T1, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat.T1)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat.T1, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.T1, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.T1, pattern = "incoming")
ht1 + ht2


library(NMF)
library(ggalluvial)
selectK(cellchat.T1, pattern = "outgoing")
nPatterns = 4
cellchat.T1 <- identifyCommunicationPatterns(cellchat.T1, pattern = "outgoing", k = nPatterns)

# river plot
netAnalysis_river(cellchat.T1, pattern = "outgoing")
#> Please make sure you have load `library(ggalluvial)` when running this function

cellchat.T1 <- computeNetSimilarity(cellchat.T1, type = "functional")
cellchat.T1 <- netEmbedding(cellchat.T1, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat.T1 <- netClustering(cellchat.T1, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat.T1, type = "functional", label.size = 3.5)

saveRDS(cellchat.T1,'cellchat.T1.Rds')


########## 合并测试
cellchat.C1 <- readRDS("cellchat_C1.Rds")
object.list <- list(C1 = cellchat.C1, T1 = cellchat.T1)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
cellchat
#> An object of class CellChat created from a merged object with multiple datasets 
#>  555 signaling genes.
#>  7563 cells. 
#> CellChat analysis of single cell RNA-seq data!
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
# visualize the enriched ligands in the first condition
computeEnrichmentScore(net.down, species = 'mouse')

pathways.show <- c("CD45") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))