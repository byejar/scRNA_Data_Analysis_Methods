
library(Seurat)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(patchwork)


######
rm(list=ls())
load(file='/home/wus/c_t_all_9_3.RData')
DimPlot(immune.combined.sct, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
    repel = TRUE)

immune.combined.sct$celltype.stim <- paste(immune.combined.sct$seurat_clusters, immune.combined.sct$orig.ident,
    sep = "_")
Idents(immune.combined.sct) <-"orig.ident"

immune.combined.sct <- PrepSCTFindMarkers(immune.combined.sct)

save(immune.combined.sct,file='/home/wus/c_t_allchange_step_1_9_13.RData')

head(immune.combined.sct)
b.interferon.response <- FindMarkers(immune.combined.sct, assay = "SCT", ident.1 = "C", ident.2 = "T",
    verbose = FALSE)
head(b.interferon.response, n = 15)
save(immune.combined.sct,b.interferon.response,file='/home/wus/c_t_allchange_step_2_9_13.RData')


####
rm(list=ls())
load(file='/home/wus/c_t_allchange_step_2_9_13.RData')
mac.cells =immune.combined.sct
Idents(mac.cells) <- "orig.ident"
head(mac.cells)
#mac.cells <- NormalizeData(object = mac.cells, assay = "RNA")
#avg.mac.cells <- log1p(AverageExpression(mac.cells, verbose = FALSE)$RNA)
avg.mac.cells <- log1p(AverageExpression(mac.cells, verbose = FALSE)$SCT)
head(avg.mac.cells)

avg.mac.cells1=as.data.frame(avg.mac.cells)
class(avg.mac.cells1)
avg.mac.cells1$gene<-rownames(avg.mac.cells1)


genes.to.label = c("Fn1", "Saa3", "Spp1",'Mrc1','Arg1','Il1b','Nos2')
p1 <- ggplot(avg.mac.cells1, aes(C, T)) + geom_point() + ggtitle("macrophage")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)

plot_grid(p1)


#散点图，带自动注释
#mac.cells <- subset(immune.combined.sct, idents = "2")
#Idents(mac.cells) <- "seurat_clusters"
mac_sub = mac.cells
Idents(mac_sub) <- "orig.ident"
head(mac_sub)
#mac_sub <- NormalizeData(object = mac_sub, assay = "RNA")
avg.mac_sub <- log1p(AverageExpression(mac_sub, verbose = FALSE)$SCT)
head(avg.mac_sub)

avg.mac_sub1=as.data.frame(avg.mac_sub)
class(avg.mac_sub1)
avg.mac_sub1$gene<-rownames(avg.mac_sub1)


#genes.to.label = c("Fn1", "Saa3", "Spp1",'Arg1','Ccl5')
p1 <- ggplot(avg.mac_sub1, aes(C, T)) + geom_point() + ggtitle("macrophage")
#p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)

plot_grid(p1)

avg.mac_sub1.marker=b.interferon.response
avg.mac_sub1.marker %>%
    slice_max(n = 10, order_by = avg_log2FC)-> p_top10

avg.mac_sub1.marker %>%
    slice_min(n = 10, order_by = avg_log2FC)-> n_top10

p_list=rownames(p_top10)
n_list=rownames(n_top10)

p_top10
n_top10

p1 <- LabelPoints(plot = p1, points = p_list, repel = TRUE)
p1 <- LabelPoints(plot = p1, points = n_list, repel = TRUE)
plot_grid(p1)

#保存前50
avg.mac_sub1.marker %>%
    slice_max(n = 50, order_by = avg_log2FC)-> p_top50

avg.mac_sub1.marker %>%
    slice_min(n = 50, order_by = avg_log2FC)-> n_top50
write.csv(p_top50,'/home/wus/for_seurat/9_6/C_T_all_p_top50_9_6.csv')
write.csv(n_top50,'/home/wus/for_seurat/9_6/C_T_all_n_top50_9_6.csv')
