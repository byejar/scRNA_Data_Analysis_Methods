rm(list=ls())
library(Seurat)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(patchwork)
load(file='/home/wus/c_t_all_9_3.RData')
DimPlot(immune.combined.sct, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
    repel = TRUE)

immune.combined.sct$celltype.stim <- paste(immune.combined.sct$seurat_clusters, immune.combined.sct$orig.ident,
    sep = "_")
Idents(immune.combined.sct) <- "celltype.stim"

immune.combined.sct <- PrepSCTFindMarkers(immune.combined.sct)
head(immune.combined.sct)
b.interferon.response <- FindMarkers(immune.combined.sct, assay = "SCT", ident.1 = "4_C", ident.2 = "4_T",
    verbose = FALSE)
head(b.interferon.response, n = 15)

b.interferon.response <- FindMarkers(immune.combined.sct, assay = "SCT", ident.1 = "4_C", ident.2 = "4_T",
    verbose = FALSE)
head(b.interferon.response, n = 15)
write.csv(b.interferon.response,'/home/wus/for_seurat/9_6/c_t_neu_all_markers_gene_9_25.csv')


save(b.interferon.response,immune.combined.sct,file='/home/wus/c_t_neu_9_6.RData')

##
rm(list=ls())
load(file='/home/wus/c_t_neu_9_6.RData')

Idents(immune.combined.sct) <- "seurat_clusters"
mac.cells <- subset(immune.combined.sct, idents = "4")
Idents(mac.cells) <- "orig.ident"
head(mac.cells)
avg.mac.cells <- log2((AverageExpression(mac.cells, verbose = FALSE)$SCT)+1)
head(avg.mac.cells)
#avg.mac.cells <- AverageExpression(mac.cells, verbose = FALSE)$SCT
#head(avg.mac.cells)

avg.mac.cells1=as.data.frame(avg.mac.cells)
class(avg.mac.cells1)
avg.mac.cells1$gene<-rownames(avg.mac.cells1)


p1 <- ggplot(avg.mac.cells1, aes(C, T)) + geom_point() + ggtitle("neu")


p1


avg.mac_sub1=avg.mac.cells1

avg.mac_sub1.marker=b.interferon.response
avg.mac_sub1.marker %>%
    slice_max(n = 10, order_by = avg_log2FC)-> p_top10

avg.mac_sub1.marker %>%
    slice_min(n = 10, order_by = avg_log2FC)-> n_top10

p_list=rownames(p_top10)
n_list=rownames(n_top10)

p1 <- LabelPoints(plot = p1, points = p_list, repel = TRUE)
p1 <- LabelPoints(plot = p1, points = n_list, repel = TRUE)
plot_grid(p1)

p_top10
n_top10

#保存前50
avg.mac_sub1.marker %>%
    slice_max(n = 50, order_by = avg_log2FC)-> p_top50

avg.mac_sub1.marker %>%
    slice_min(n = 50, order_by = avg_log2FC)-> n_top50
write.csv(p_top50,'/home/wus/for_seurat/9_6/C_T_neu_p_top50_9_6-v2.csv')
write.csv(n_top50,'/home/wus/for_seurat/9_6/C_T_neu_n_top50_9_6-v2.csv')

####更改画图
head(avg.mac_sub1)

col_dataset=avg.mac_sub1
col_dataset$log2fc=col_dataset$C-col_dataset$T
col_dataset$up_down=as.factor(ifelse(col_dataset$log2fc<0.6,ifelse(col_dataset$log2fc<(-0.6),'down','nosig'),'up'))

pc<-ggplot(col_dataset, aes(x=C, y=T,colour=up_down)) + 
  geom_point(shape=20,size=3,aes(colour=up_down))+
  scale_colour_manual(values=c('up'='#b4292c', 'nosig'="#c8c8ca",'down'='#5177b9'),name='up_or_down') + 
theme_bw()+
theme(axis.title=element_text(face="bold", size=28,colour = 'black'), #坐标轴标题
                 axis.text=element_text(face="bold", size=24,colour = 'black'))+ #坐标轴标签
scale_x_continuous(expand = c(0,0.005),limits=c(0,7))+scale_y_continuous(expand = c(0,0.01),limits=c(0,7))+
geom_abline(slope = 1,intercept = 0,color='black',size=1)+
geom_abline(slope = 1,intercept = 0.6,color='black',linetype = "dashed",size=1)+
geom_abline(slope = 1,intercept = -0.6,color='black',linetype = "dashed",size=1)+
  ylab("log2 T_mean")+
  xlab("log2 C_mean")
pc+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

pc <- LabelPoints(plot = pc, points = p_list, repel = TRUE)
pc <- LabelPoints(plot = pc, points = n_list, repel = TRUE)
plot_grid(pc)








####T
rm(list=ls())
library(Seurat)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(patchwork)
load(file='/home/wus/c_t_all_9_3.RData')
DimPlot(immune.combined.sct, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
    repel = TRUE)

immune.combined.sct$celltype.stim <- paste(immune.combined.sct$seurat_clusters, immune.combined.sct$orig.ident,
    sep = "_")
Idents(immune.combined.sct) <- "celltype.stim"

immune.combined.sct <- PrepSCTFindMarkers(immune.combined.sct)
head(immune.combined.sct)
b.interferon.response <- FindMarkers(immune.combined.sct, assay = "SCT", ident.1 = "6_C", ident.2 = "6_T",
    verbose = FALSE)
head(b.interferon.response, n = 15)


#总表
b.interferon.response <- FindMarkers(immune.combined.sct, assay = "SCT", ident.1 = "6_C", ident.2 = "6_T",
    verbose = FALSE)
head(b.interferon.response, n = 15)
write.csv(b.interferon.response,'/home/wus/for_seurat/9_6/c_t_tcell_all_markers_gene_9_12.csv')

save(b.interferon.response,immune.combined.sct,file='/home/wus/c_t_T_cell_9_6.RData')


##
rm(list=ls())
load(file='/home/wus/c_t_T_cell_9_6.RData')

#write.csv(b.interferon.response,'/home/wus/for_seurat/9_6/C_T_tcell_findmarker_9_6.csv')

###识别不同条件下的差异表达基因（带自己的注释，散点图）
Idents(immune.combined.sct) <- "seurat_clusters"
mac.cells <- subset(immune.combined.sct, idents = "6")
Idents(mac.cells) <- "orig.ident"
head(mac.cells)
#mac.cells <- NormalizeData(object = mac.cells, assay = "RNA")
avg.mac.cells <- log1p(AverageExpression(mac.cells, verbose = FALSE)$SCT)
head(avg.mac.cells)

avg.mac.cells1=as.data.frame(avg.mac.cells)
class(avg.mac.cells1)
avg.mac.cells1$gene<-rownames(avg.mac.cells1)


#genes.to.label = c("Fn1", "Saa3", "Spp1",'Mrc1','Arg1','Il1b','Nos2')
#p1 <- ggplot(avg.mac.cells1, aes(C, T)) + geom_point() + ggtitle("macrophage")
#p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)

#plot_grid(p1)


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

p1 <- LabelPoints(plot = p1, points = p_list, repel = TRUE)
p1 <- LabelPoints(plot = p1, points = n_list, repel = TRUE)
plot_grid(p1)

p_top10
n_top10
#保存前50
avg.mac_sub1.marker %>%
    slice_max(n = 50, order_by = avg_log2FC)-> p_top50

avg.mac_sub1.marker %>%
    slice_min(n = 50, order_by = avg_log2FC)-> n_top50
write.csv(p_top50,'/home/wus/for_seurat/9_6/C_T_tcell_p_top50_9_6-v2.csv')
write.csv(n_top50,'/home/wus/for_seurat/9_6/C_T_tcell_n_top50_9_6-v2.csv')


#####对t cell进行分群
rm(list=ls())
library(Seurat)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(patchwork)
load(file='/home/wus/c_t_mac_9_6.RData')

mac.cells <- subset(immune.combined.sct, idents = "6")

pancreas.list <- SplitObject(mac.cells, split.by = "orig.ident")
#C_group <- pancreas.list[c("C1","C2","C3")]
C_group <- pancreas.list[["C"]]
#T_group <- pancreas.list[c('T1','T2','T3')]
T_group <- pancreas.list[["T"]]

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
mac.cells <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")




mac.cells <- RunPCA(mac.cells, verbose = FALSE)
mac.cells <- RunUMAP(mac.cells, reduction = "pca", dims = 1:30, verbose = FALSE)
mac.cells <- FindNeighbors(mac.cells, reduction = "pca", dims = 1:30)
mac.cells <- FindClusters(mac.cells, resolution = 0.3)
#并排可视化这两个条件，使用参数来显示按聚类着色的每个条件。split.by
DimPlot(mac.cells, reduction = "umap", split.by = "orig.ident")
DimPlot(mac.cells, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
    repel = TRUE)


##差异表达基因
pancreas.integrated.markers<- FindAllMarkers(mac.cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pancreas.integrated.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)


pancreas.integrated.markers %>%
    group_by(cluster) %>%
    top_n(n = 20, wt = avg_log2FC) -> top20

#DoHeatmap(pancreas.integrated, features = top20$gene) + NoLegend()
DimPlot(mac.cells, reduction = 'umap', group.by = 'seurat_clusters',
           label = TRUE, pt.size = 0.5) + NoLegend()

write.csv(top20,'/home/wus/for_seurat/8_2/CT_t_cluster_merge_top20_gene_for_marker_9_3.csv')
write.csv(pancreas.integrated.markers,'/home/wus/for_seurat/8_2/CT_t_cluster_merge_markers_gene_9_3.csv')

save(mac.cells,pancreas.integrated.markers,file='/home/wus/c_t_t_cluster_step_1_9_3.RData')

#NK
load(file='/home/wus/c_t_t_cluster_step_1_9_3.RData')
Idents(mac.cells)
immune_combine=mac.cells

Idents(immune_combine) <- "seurat_clusters"
immune_combine <- RenameIdents(immune_combine, `0` = "T cell", `1` = "T cell", `2` = "T cell", 
    `3` = "T cell", `4` = "NK", `5` = "T cell")
DimPlot(immune_combine, reduction = "umap", label = TRUE,
    repel = TRUE)
######

DefaultAssay(mac.cells) <- "SCT"
FeaturePlot(mac.cells, features = c('Cd3e','Cd4','Cd8a','Cd8b1','Klrb1c','Gzma'),  max.cutoff = 3,
    cols = c("grey", 'orange',"red"))

##NK的差异基因
mac.cells$celltype.stim <- paste(mac.cells$seurat_clusters, mac.cells$orig.ident,
    sep = "_")
Idents(mac.cells) <- "celltype.stim"
mac.cells <- PrepSCTFindMarkers(mac.cells)
head(mac.cells)
b.interferon.response <- FindMarkers(mac.cells, assay = "SCT", ident.1 = "4_C", ident.2 = "4_T",
    verbose = FALSE)
head(b.interferon.response, n = 15)
#####散点图

mac.cells <- subset(immune_combine, idents = "NK")
Idents(mac.cells) <- "orig.ident"
head(mac.cells)

avg.mac.cells <- log2((AverageExpression(mac.cells, verbose = FALSE)$SCT)+1)

head(avg.mac.cells)


avg.mac.cells1=as.data.frame(avg.mac.cells)
class(avg.mac.cells1)
avg.mac.cells1$gene<-rownames(avg.mac.cells1)


p1 <- ggplot(avg.mac.cells1, aes(C, T)) + geom_point() + ggtitle("mac")


p1


avg.mac_sub1=avg.mac.cells1

avg.mac_sub1.marker=b.interferon.response
avg.mac_sub1.marker %>%
    slice_max(n = 15, order_by = avg_log2FC)-> p_top10

avg.mac_sub1.marker %>%
    slice_min(n = 15, order_by = avg_log2FC)-> n_top10

p_list=rownames(p_top10)
n_list=rownames(n_top10)

p1 <- LabelPoints(plot = p1, points = p_list, repel = TRUE)
p1 <- LabelPoints(plot = p1, points = n_list, repel = TRUE)
plot_grid(p1)

p_top10
n_top10

#保存前50
#avg.mac_sub1.marker %>%
#    slice_max(n = 50, order_by = avg_log2FC)-> p_top50

#avg.mac_sub1.marker %>%
#    slice_min(n = 50, order_by = avg_log2FC)-> n_top50
#write.csv(p_top50,'/home/wus/for_seurat/9_6/C_T_p_cancer_top50_9_6.csv')
#write.csv(n_top50,'/home/wus/for_seurat/9_6/C_T_n_cancer_top50_9_6.csv')

####更改画图
head(avg.mac_sub1)

col_dataset=avg.mac_sub1
col_dataset$log2fc=col_dataset$C-col_dataset$T
col_dataset$up_down=as.factor(ifelse(col_dataset$log2fc<0.6,ifelse(col_dataset$log2fc<(-0.6),'down','nosig'),'up'))
sort<-col_dataset[order(col_dataset$log2fc),]
head(sort)
pc<-ggplot(col_dataset, aes(x=C, y=T,colour=up_down)) + 
  geom_point(shape=20,size=3,aes(colour=up_down))+
  scale_colour_manual(values=c('up'='#b4292c', 'nosig'="#c8c8ca",'down'='#5177b9'),name='up_or_down') + 
theme_bw()+
theme(axis.title=element_text(face="bold", size=28,colour = 'black'), #坐标轴标题
                 axis.text=element_text(face="bold", size=24,colour = 'black'))+ #坐标轴标签
scale_x_continuous(expand = c(0,0.005),limits=c(0,7))+scale_y_continuous(expand = c(0,0.01),limits=c(0,7))+
geom_abline(slope = 1,intercept = 0,color='black',size=1)+
geom_abline(slope = 1,intercept = 0.6,color='black',linetype = "dashed",size=1)+
geom_abline(slope = 1,intercept = -0.6,color='black',linetype = "dashed",size=1)+
  ylab("log2 T_mean")+
  xlab("log2 C_mean")
pc+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

pc <- LabelPoints(plot = pc, points = p_list, repel = TRUE)
pc <- LabelPoints(plot = pc, points = n_list, repel = TRUE)
plot_grid(pc)



#总表
Idents(mac.cells) <- "celltype.stim"
b.interferon.response <- FindMarkers(mac.cells, assay = "SCT", ident.1 = "4_C", ident.2 = "4_T",
    verbose = FALSE)
head(b.interferon.response, n = 15)
write.csv(b.interferon.response,'/home/wus/for_seurat/9_6/c_t_NK_all_markers_gene_9_12.csv')





###T除去NK的t cell
load(file='/home/wus/c_t_t_cluster_step_1_9_3.RData')
Idents(mac.cells)
immune_combine=mac.cells

Idents(immune_combine) <- "seurat_clusters"
immune_combine <- RenameIdents(immune_combine, `0` = "T cell", `1` = "T cell", `2` = "T cell", 
    `3` = "T cell", `4` = "NK", `5` = "T cell")
DimPlot(immune_combine, reduction = "umap", label = TRUE,
    repel = TRUE)

#总表
mac.cells$celltype.stim <- paste(mac.cells$seurat_clusters, mac.cells$orig.ident,
    sep = "_")
Idents(mac.cells) <- "celltype.stim"
mac.cells <- PrepSCTFindMarkers(mac.cells)
head(mac.cells)
b.interferon.response <- FindMarkers(mac.cells, assay = "SCT", ident.1 = c("1_C",'2_C','3_C','5_C'), ident.2 = c("1_T",'2_T','3_T','5_T'),
    verbose = FALSE)
head(b.interferon.response, n = 15)

b.interferon.response <- FindMarkers(mac.cells, assay = "SCT", ident.1 = c("1_C",'2_C','3_C','5_C') ,ident.2 = c("1_T",'2_T','3_T','5_T'),
    verbose = FALSE)
head(b.interferon.response, n = 15)
write.csv(b.interferon.response,'/home/wus/for_seurat/9_6/c_t_tcell_all_markers_gene_9_25.csv')

######

#DefaultAssay(mac.cells) <- "SCT"
#FeaturePlot(mac.cells, features = c('Cd3e','Cd4','Cd8a','Cd8b1','Klrb1c','Gzma'),  max.cutoff = 3,
#    cols = c("grey", 'orange',"red"))

##差异基因


#####散点图
Idents(immune_combine) <- "seurat_clusters"
immune_combine <- RenameIdents(immune_combine, `0` = "T cell", `1` = "T cell", `2` = "T cell", 
    `3` = "T cell", `4` = "NK", `5` = "T cell")
mac.cells <- subset(immune_combine, idents = "T cell")
Idents(mac.cells) <- "orig.ident"
head(mac.cells)

avg.mac.cells <- log2((AverageExpression(mac.cells, verbose = FALSE)$SCT)+1)

head(avg.mac.cells)


avg.mac.cells1=as.data.frame(avg.mac.cells)
class(avg.mac.cells1)
avg.mac.cells1$gene<-rownames(avg.mac.cells1)


p1 <- ggplot(avg.mac.cells1, aes(C, T)) + geom_point() + ggtitle("mac")


p1


avg.mac_sub1=avg.mac.cells1

avg.mac_sub1.marker=b.interferon.response
avg.mac_sub1.marker %>%
    slice_max(n = 15, order_by = avg_log2FC)-> p_top10

avg.mac_sub1.marker %>%
    slice_min(n = 15, order_by = avg_log2FC)-> n_top10

p_list=rownames(p_top10)
n_list=rownames(n_top10)

p1 <- LabelPoints(plot = p1, points = p_list, repel = TRUE)
p1 <- LabelPoints(plot = p1, points = n_list, repel = TRUE)
plot_grid(p1)

p_top10
n_top10


####更改画图
head(avg.mac_sub1)

col_dataset=avg.mac_sub1
col_dataset$log2fc=col_dataset$C-col_dataset$T
col_dataset$up_down=as.factor(ifelse(col_dataset$log2fc<0.6,ifelse(col_dataset$log2fc<(-0.6),'down','nosig'),'up'))
sort<-col_dataset[order(col_dataset$log2fc),]
head(sort)
pc<-ggplot(col_dataset, aes(x=C, y=T,colour=up_down)) + 
  geom_point(shape=20,size=3,aes(colour=up_down))+
  scale_colour_manual(values=c('up'='#b4292c', 'nosig'="#c8c8ca",'down'='#5177b9'),name='up_or_down') + 
theme_bw()+
theme(axis.title=element_text(face="bold", size=28,colour = 'black'), #坐标轴标题
                 axis.text=element_text(face="bold", size=24,colour = 'black'))+ #坐标轴标签
scale_x_continuous(expand = c(0,0.005),limits=c(0,7))+scale_y_continuous(expand = c(0,0.01),limits=c(0,7))+
geom_abline(slope = 1,intercept = 0,color='black',size=1)+
geom_abline(slope = 1,intercept = 0.6,color='black',linetype = "dashed",size=1)+
geom_abline(slope = 1,intercept = -0.6,color='black',linetype = "dashed",size=1)+
  ylab("log2 T_mean")+
  xlab("log2 C_mean")
pc+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

pc <- LabelPoints(plot = pc, points = p_list, repel = TRUE)
pc <- LabelPoints(plot = pc, points = n_list, repel = TRUE)
plot_grid(pc)


