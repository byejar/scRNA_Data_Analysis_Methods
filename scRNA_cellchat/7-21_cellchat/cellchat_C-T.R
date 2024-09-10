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
data.dir <- '/home/wus/scRNA_for_fig/7-20_CellChat'
#dir.create(data.dir)
setwd(data.dir)

########## 合并测试
cellchat.C <- readRDS("cellchat.C.Rds")
cellchat.T <- readRDS("cellchat.T.Rds")
object.list <- list(C = cellchat.C, T = cellchat.T)
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


par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")


gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2


gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

# visualize the enriched ligands in the first condition
computeEnrichmentScore(net.down, species = 'mouse')

pathways.show <- c("SPP1") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))



vertex.receiver = seq(1,4) # a numeric vector. 
netAnalysis_contribution(cellchat.T, signaling = pathways.show)
netVisual_aggregate(cellchat.T, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout='hierarchy')

netAnalysis_contribution(cellchat.C, signaling = pathways.show)
netVisual_aggregate(cellchat.C, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout='hierarchy')



pathways.show <- c("CD52") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
#vertex.receiver = c('T cell')
netVisual_aggregate(cellchat.C, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout='hierarchy')

pathways.show <- c("PTN") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netAnalysis_contribution(cellchat.T, signaling = pathways.show)
netVisual_aggregate(cellchat.T, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout='hierarchy')

pairLR.CXCL <- extractEnrichedLR(cellchat.T, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
LR.show <- 'PTN_NCL'
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat.T, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver,layout='hierarchy')
#> [[1]]
# Circle plot
netVisual_individual(cellchat.T, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

plotGeneExpression(cellchat.C, signaling = "SPP1")
plotGeneExpression(cellchat.T, signaling = "MK")

########## 验证 spp1表达
load(file = '/home/wus/for_seurat/8_1/immune_rename_2023-7-4.Rdata')

head(mac_sce_rename)

Idents(mac_sce_rename)="cell_type"

VlnPlot(mac_sce_rename, 
        features = c("Spp1","Cd44","Itgav","Itga4","Itga9","Itga5","Itga3","Itga1"),
        pt.size = 0,
        ncol = 3,
         split.by = "orig.ident") 


netVisual_bubble(cellchat, sources.use = 'M1', targets.use = c(5:11),  comparison = c(1, 2), angle.x = 45)
#> Comparing communications on a merged object

netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  comparison = c(1, 2), angle.x = 45)
#> Comparing communications on a merged object


num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)



gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "M2", signaling.exclude = "MIF")
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0
#> The following `from` values were not present in `x`: 0, -1
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "T cell", signaling.exclude = c("MIF"))
#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
patchwork::wrap_plots(plots = list(gg1,gg2))



# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "T"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "T",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "C",ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
# gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
# #> Comparing communications on a merged object
# pairLR.use.down = net.down[, "interaction_name", drop = F]
# gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
# #> Comparing communications on a merged object
# gg1 + gg2


gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 'M1',  comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 'M1',  comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2

gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 'M2',  comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 'M2',  comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2

computeEnrichmentScore(net.down, species = 'mouse')


pathways.show <- c("SPP1") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}



########## CSDN
#如果不需要多线程，可以试试添加do.parallel = FALSE参数
#cellchat_clustering <- netClustering(cellchat, type = "functional",do.parallel = FALSE)


cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2

cellchat <- netClustering(cellchat, type = "functional",do.parallel = FALSE)
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2


cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "structural",do.parallel = FALSE)
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2



netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)
#> 2D visualization of signaling networks from datasets 1 2


rankSimilarity(cellchat, type = "functional")
#> Compute the distance of signaling networks between datasets 1 2

#### 查看两组FGF功能变化
pathways.show <- c("FGF") 
vertex.receiver = seq(1,4) # a numeric vector. 
netAnalysis_contribution(cellchat.T, signaling = pathways.show)
netVisual_aggregate(cellchat.T, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout='hierarchy')


pathways.show <- c("FGF") 
vertex.receiver = seq(1,4) # a numeric vector. 
netAnalysis_contribution(cellchat.C, signaling = pathways.show)
netVisual_aggregate(cellchat.C, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout='hierarchy')


pathways.show <- c("CXCL") 
vertex.receiver = seq(1,4) # a numeric vector. 
netAnalysis_contribution(cellchat.T, signaling = pathways.show)
netVisual_aggregate(cellchat.T, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout='hierarchy')


netAnalysis_contribution(cellchat.C, signaling = pathways.show)
netVisual_aggregate(cellchat.C, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout='hierarchy')



########### 检验
gc<-plotGeneExpression(cellchat.C, signaling = "CD52")
gt<-plotGeneExpression(cellchat.T, signaling = "CD52")
gc+gt

load(file = '/home/wus/for_seurat/8_1/immune_rename_2023-7-4.Rdata')
Idents(mac_sce_rename)="cell_type"
VlnPlot(mac_sce_rename, 
        features = c("Siglecg","Cd52"),
        pt.size = 0,
        ncol = 3,
         split.by = "orig.ident") 

gc<-plotGeneExpression(cellchat.C, signaling = "CLEC")
gt<-plotGeneExpression(cellchat.T, signaling = "CLEC")
gc+gt

### NK
FeaturePlot(mac_sce_rename, features = c('Klrb1c','Gzma'),  max.cutoff = 3,
    cols = c("grey", 'orange',"red"))