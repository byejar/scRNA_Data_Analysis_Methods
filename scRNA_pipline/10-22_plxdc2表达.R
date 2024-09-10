library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(cowplot)
library(paletteer)  
library(gplots)
library(ggpubr)    
library(ggsci) 
library(stringr)

rm(list=ls())
load(file='/home/wus/c_t_mac_9_12.RData')
Idents(immune.combined.sct)
immune_combine=immune.combined.sct
Idents(immune_combine) <- "seurat_clusters"
immune_combine <- RenameIdents(immune_combine, `0` = "cancer cell", `1` = "cancer cell", `2` = "macrophage", 
    `3` = "cancer cell", `4` = "neurtophil", `5` = "fibroblast", `6` = "T cell", `7` = "7", `8` = "cancer cell", `9` = "DC", `10` = "B cell", `11` = "11")

DimPlot(immune_combine, reduction = "umap", label = TRUE,
    repel = TRUE)

mac.cells <- subset(immune_combine , idents = c("cancer cell","macrophage","fibroblast","neurtophil", "T cell","DC","B cell"))
DimPlot(mac.cells, reduction = "umap", label = TRUE,
    repel = TRUE)

##voline
plots <- VlnPlot(mac.cells, features = c("Plxdc2"), 
     pt.size = 1, combine = FALSE,assay = "RNA")
wrap_plots(plots = plots, ncol = 1)

AverageExpression(mac.cells, features = c("Plxdc2"))

###提取plxdc2 和 mrc1基因表达
sce<- subset(immune_combine , idents = c("macrophage"))
expr <- sce@assays$SCT@scale.data
gene_name <- c("Plxdc2", "Mrc1")
gene_expression <- expr[gene_name,] %>% 
  as.data.frame()
gene_expression <- as.data.frame(expr[gene_name,])
gene_expression2 <- as.data.frame(t(gene_expression))  

normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
gene_expression2_normal=as.data.frame(lapply(gene_expression2,normalize))
head(gene_expression2_normal)

write.csv(gene_expression2_normal,'/home/wus/mahm_suv/gene_expression2_normal.csv')

source('/home/wus/mahm_suv/cor_addtext.R')

path='/home/wus/mahm_suv/'
setwd=path
scRNA<-read.csv('/home/wus/mahm_suv/gene_expression2_normal.csv')
cor_addtext(scRNA,'scRNA','Plxdc2','Mrc1',path)

####或者使用自带的方法
FeatureScatter(sce,feature1='Plxdc2',feature2='Mrc1')
FeatureScatter(sce,feature1='Plxdc2',feature2='Arg1')
FeatureScatter(sce,feature1='Plxdc2',feature2='Il1b')
FeatureScatter(sce,feature1='Plxdc2',feature2='Ccl5')
#colnames(gene_expression) <- paste0(gene_name)
#identical(colnames(sce),row.names(gene_expression))
#sce$CXCL10 <- gene_expression[,paste0(gene_name)]  # CXCL10, 这里需要修改为目标基因
#identical(sce@meta.data[,paste0(gene_name)],gene_expression[,paste0(gene_name)])


###M1 M2区分并加入

DefaultAssay(immune.combined.sct) <- "SCT"
Idents(immune_combine) <- "seurat_clusters"


mac.cells <- subset(immune.combined.sct, idents = "2")

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
mac.cells <- RunUMAP(mac.cells, reduction = "pca", dims = 1:20, verbose = FALSE)
mac.cells <- FindNeighbors(mac.cells, reduction = "pca", dims = 1:20)
mac.cells <- FindClusters(mac.cells, resolution = 0.2)
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

#write.csv(top20,'/home/wus/for_seurat/8_2/CT_macs_merge_top20_gene_for_marker_10_18.csv')
#write.csv(pancreas.integrated.markers,'/home/wus/for_seurat/8_2/CT_macs_merge_markers_gene_10_18.csv')

#save(mac.cells,pancreas.integrated.markers,file='/home/wus/c_t_macs_step_1_10_18.RData')


#M2 M2 数量
table(mac.cells$seurat_clusters,mac.cells$orig.ident )


##气泡图并排联合展示
Idents(mac.cells) <- "seurat_clusters"
mac.cells <- RenameIdents(mac.cells, `0` = "M2", `1` = "M2", `2` = "M1", 
    `3` = "other")
mac.cells <- subset(mac.cells, idents = c("M2",'M1'))


sub_cells = mac.cells
sce=immune_combine
mac_sce1=immune_combine
head(sub_cells)

## Rename sub cells
metadata <- sub_cells@meta.data
metadata$cluster <- Idents(sub_cells)
metadata$sub_cluster <- paste0('sub', metadata$cluster)
sub_cells@meta.data <- metadata
table(sub_cells$sub_cluster)

## Mapping to old clusters
metadata <- sce@meta.data
metadata$cell_type <- NA
metadata$cell_type <- sub_cells@meta.data[match(rownames(sce@meta.data), rownames(sub_cells@meta.data)), 'sub_cluster']
print(table(metadata$cell_type))
head(metadata)


celltype_names <- NULL
for(i in 1:dim(metadata)[1]){
  # print(i)
  sub_data <- metadata[i,]
  if(is.na(sub_data$cell_type)){
    # print('Change value')
    sub_data$cell_type <- sub_data$seurat_clusters
    celltype_names <- c(celltype_names, sub_data$cell_type)
  }else{
    celltype_names <- c(celltype_names, sub_data$cell_type)
  }
}

print(table(celltype_names))

metadata$cell_type <- celltype_names
metadata$cell_type <- factor(metadata$cell_type)
metadata$cell_type <- factor(levels = c( 1:nlevels(mac_sce1$seurat_clusters),names(table(sub_cells$sub_cluster))),metadata$cell_type)
print(table(metadata$cell_type))

mac_sce1@meta.data <- metadata 
head(mac_sce1)
DimPlot(mac_sce1, reduction='umap', label=T, label.size=5, pt.size=1, group.by='cell_type') 

#
cluster2celltype <- c(
                      "1"="cancer_cell", 
                      "2"="cancer_cell", 
                      "3"= "undifine", 
                      "4"= "cancer cell",
                      "5"= "neutrophil",
                       #"sub0"= "naive",
                      "6"= "fib", 
                      "7"= "T cell",
                      
                      "8"= "undifine",
                      "9"= "cancer cell",
                      '10'='DC',
                      "11"= "B cell",
                      "12"= "undifine",
                      "subM1"= "M1",
                      "subM2"= "M2"
                      
                      )
mac_sce_rename=mac_sce1
mac_sce_rename[['cell_type']] = unname(cluster2celltype[mac_sce_rename@meta.data$cell_type])

DimPlot(mac_sce_rename, reduction = 'umap', group.by = 'cell_type',
        label = TRUE, pt.size = 0.5) + NoLegend()


save(mac_sce_rename,file = '/home/wus/for_seurat/8_1/T_immune_rename_10_18.Rdata')

#取子集
load(file = '/home/wus/for_seurat/8_1/T_immune_rename_10_18.Rdata')
sce=mac_sce_rename
head(mac_sce_rename)
head(sce)
mac_sce_re<-subset(x=sce,subset=cell_type==c('M1',"M2","T cell","B cell","cancer cell","DC",'neutrophil','fib'))

#mac_sce_re = sce[,sce@meta.data$seurat_clusters %in% c(0,1,2,3,6,13,14,15)]

DimPlot(mac_sce_re, reduction = 'umap', group.by = 'cell_type',
        label = TRUE, pt.size = 0.5) + NoLegend()


plots <- VlnPlot(mac_sce_re, features = c("Plxdc2"), group.by = 'cell_type',
     pt.size = 1, combine = FALSE,assay = "SCT")
wrap_plots(plots = plots, ncol = 1)

AverageExpression(mac_sce_re, features = c("Plxdc2"),group.by = 'cell_type')


Idents(mac_sce_re) <- "cell_type"
df <- as.data.frame(GetAssayData(mac_sce_re, assay = "SCT", slot = "counts"))#获取标准化后表达矩阵数据
pbmc_idents <- as.data.frame(t(as.data.frame(Idents(object = mac_sce_re))))#获取细胞ID和亚群的对应关系
identical(colnames(df),colnames(pbmc_idents))#判断顺序
colnames(df) <- pbmc_idents[1,]#把细胞ID替换成亚群名称（有重复所以使用不了everything函数）
x <- data.frame(GeneSymbol = rownames(df))
rownames(df) <- NULL
final <- cbind(x,df)
#得到没有观测名，第一列是Gene symbol后面是细胞亚群对应的表达矩阵
#这个就是进行CIBERSORTX的单细胞参考表达数据
write.table(final,'/home/wus/for_seurat/8_1/final.txt',row.names=FALSE)