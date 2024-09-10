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



###########抽样
mac.cells.C=subset(mac.cells,orig.ident %in% c("C1","C3"))
mac.cells.T=subset(mac.cells,orig.ident %in% c("T1","T2","T3"))

DimPlot(subset(mac.cells.C, downsample = 10000), reduction = "umap", label = TRUE,
    repel = TRUE,pt.size = 0.01)

DimPlot(subset(mac.cells.T, downsample = 10000), reduction = "umap", label = TRUE,
    repel = TRUE,pt.size = 0.01)


######################################################################################### tumour fig
mac.cells <- subset(immune_combine, idents = "cancer cell")


mac.cells$celltype.stim <- substr(mac.cells$orig.ident, 1, 1) 

Idents(mac.cells) <- "celltype.stim"


# Idents(mac.cells) <- "orig.ident"
# head(mac.cells)
# head(mac.cells)

avg.mac.cells <- log2((AverageExpression(mac.cells, verbose = FALSE)$SCT)+1)
#avg.mac.cells <- mac.cells[['SCT']]@data
head(avg.mac.cells)
#avg.mac.cells <- AverageExpression(mac.cells, verbose = FALSE)$SCT
#head(avg.mac.cells)

avg.mac.cells1=as.data.frame(avg.mac.cells)
class(avg.mac.cells1)
avg.mac.cells1$gene<-rownames(avg.mac.cells1)


p1 <- ggplot(avg.mac.cells1, aes(C, T)) + geom_point() + ggtitle("mac")


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
write.csv(p_top50,'/home/wus/for_seurat/9_6/C_T_p_cancer_top50_9_6.csv')
write.csv(n_top50,'/home/wus/for_seurat/9_6/C_T_n_cancer_top50_9_6.csv')

####更改画图
head(avg.mac_sub1)

col_dataset=avg.mac_sub1
col_dataset$log2fc=col_dataset$C-col_dataset$T
col_dataset$up_down=as.factor(ifelse(col_dataset$log2fc<0.6,ifelse(col_dataset$log2fc<(-0.6),'down','nosig'),'up'))
sort<-col_dataset[order(col_dataset$log2fc),]
head(sort)
pc<-ggplot(col_dataset, aes(x=C, y=T,colour=up_down)) + 
  geom_point(shape=20,size=3,aes(colour=up_down))+
  scale_colour_manual(values=c('up'='#5177b9', 'nosig'="#c8c8ca",'down'='#b4292c'),name='up_or_down') + 
theme_bw()+
theme(axis.title=element_text(face="bold", size=28,colour = 'black'), #坐标轴标题
                 axis.text=element_text(face="bold", size=24,colour = 'black'))+ #坐标轴标签
scale_x_continuous(expand = c(0,0.005),limits=c(0,7))+scale_y_continuous(expand = c(0,0.01),limits=c(0,7))+
#geom_abline(slope = 1,intercept = 0,color='black',size=1)+
geom_abline(slope = 1,intercept = 0.6,color='black',linetype = "dashed",size=1)+
geom_abline(slope = 1,intercept = -0.6,color='black',linetype = "dashed",size=1)+
  ylab("log2 T_mean")+
  xlab("log2 C_mean")
pc+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

#pc <- LabelPoints(plot = pc, points = p_list, repel = TRUE)
#pc <- LabelPoints(plot = pc, points = n_list, repel = TRUE)
plot_grid(pc)

ggsave(pc,filename = "/home/wus/scRNA_for_fig/tumor-v3.png",width = 6,height = 5,dpi=300)
