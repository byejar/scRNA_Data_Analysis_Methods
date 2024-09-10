####all-v2
rm(list=ls())
library(Seurat)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(patchwork)


C2 <- Read10X(data.dir = "/home/wus/scRUN/sampleC2/outs/filtered_feature_bc_matrix")
C2_obj <- CreateSeuratObject(counts = C2, project = "C", min.cells = 3, min.features = 200)
C2_obj[["percent.mt"]] <- PercentageFeatureSet(C2_obj, pattern = "^mt-")
C2_obj <- subset(C2_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 8200 & percent.mt < 5)
C2_obj

C3 <- Read10X(data.dir = "/home/wus/scRUN/sampleC3/outs/filtered_feature_bc_matrix")
C3_obj <- CreateSeuratObject(counts = C3, project = "C", min.cells = 3, min.features = 200)
C3_obj[["percent.mt"]] <- PercentageFeatureSet(C3_obj, pattern = "^mt-")
C3_obj <- subset(C3_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 8200 & percent.mt < 5)
C3_obj

C1 <- Read10X(data.dir = "/home/wus/scRUN/sampleC1/outs/filtered_feature_bc_matrix")
C1_obj <- CreateSeuratObject(counts = C1, project = "C", min.cells = 3, min.features = 200)
C1_obj[["percent.mt"]] <- PercentageFeatureSet(C1_obj, pattern = "^mt-")
C1_obj <- subset(C1_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 8200 & percent.mt < 5)
C1_obj


T1 <- Read10X(data.dir = "/home/wus/scRUN/sampleT1/outs/filtered_feature_bc_matrix")
T1_obj <- CreateSeuratObject(counts = T1, project = "T", min.cells = 3, min.features = 200)
T1_obj[["percent.mt"]] <- PercentageFeatureSet(T1_obj, pattern = "^mt-")
T1_obj <- subset(T1_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 8200 & percent.mt < 5)
T1_obj

T2 <- Read10X(data.dir = "/home/wus/scRUN/sampleT2/outs/filtered_feature_bc_matrix")
T2_obj <- CreateSeuratObject(counts = T2, project = "T", min.cells = 3, min.features = 200)
T2_obj[["percent.mt"]] <- PercentageFeatureSet(T2_obj, pattern = "^mt-")
T2_obj <- subset(T2_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 8200 & percent.mt < 5)
T2_obj

T3 <- Read10X(data.dir = "/home/wus/scRUN/sampleT3/outs/filtered_feature_bc_matrix")
T3_obj <- CreateSeuratObject(counts = T3, project = "T", min.cells = 3, min.features = 200)
T3_obj[["percent.mt"]] <- PercentageFeatureSet(T3_obj, pattern = "^mt-")
T3_obj <- subset(T3_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 8200 & percent.mt < 5)
T3_obj


sc.combined <- merge(C1_obj, y = c(C2_obj,C3_obj,T1_obj,T2_obj,T3_obj), add.cell.ids = c("C","C","C",'T','T','T'), project = "scRNA")
#head(sc.combined@meta.data)
pancreas.list <- SplitObject(sc.combined, split.by = "orig.ident")
pancreas.list <- pancreas.list[c("C",'T')]

for (i in 1:length(pancreas.list)) {
    pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
    pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000, 
        verbose = FALSE)
}

#SCTransform识别不同条件下的共有可变基因应该更合理
reference.list <- pancreas.list[c("C",'T')]
#pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
#pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors)
DefaultAssay(pancreas.integrated) <- "integrated"



pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
#pa1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "orig.ident")
#pa2 <- DimPlot(pancreas.integrated, reduction = "umap", label = TRUE, repel = TRUE) +  NoLegend()
#pa1 + pa2
DimPlot(pancreas.integrated, reduction = "umap", group.by = "orig.ident")

save(C1_obj, C2_obj,C3_obj,T1_obj,T2_obj,T3_obj,sc.combined,pancreas.list,pancreas.integrated,pancreas.anchors,file='/home/wus/c_t_all_8_31.RData')

#分类,findmarker
rm(list=ls())
load(file='/home/wus/c_t_all_8_23.RData')
# normalize and run dimensionality reduction on control dataset
ctrl <- SCTransform(pancreas.integrated, vst.flavor = "v2", verbose = FALSE) %>%
    RunPCA(npcs = 30, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
    FindClusters(resolution = 0.7, verbose = FALSE)

p1 <- DimPlot(ctrl, label = T, repel = T) + ggtitle("Unsupervised clustering")
p2 <- DimPlot(ctrl, label = T, repel = T, group.by =  "orig.ident") + ggtitle("Annotated celltypes")

p1 | p2


DimPlot(pancreas.integrated, split.by = "orig.ident")  

substitute=pancreas.integrated
substitute <- JackStraw(substitute, num.replicate = 100)
substitute <- ScoreJackStraw(substitute, dims = 1:40)

JackStrawPlot(substitute, dims = 1:40)

ElbowPlot(substitute)

pancreas.integrated <- FindNeighbors(pancreas.integrated, dims = 1:30)
pancreas.integrated <- FindClusters(pancreas.integrated, resolution = 0.5 )
DimPlot(pancreas.integrated, reduction = 'umap', group.by = 'seurat_clusters',
           label = TRUE, pt.size = 0.5) + NoLegend()
pancreas.integrated.markers<- FindAllMarkers(pancreas.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pancreas.integrated.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)


pancreas.integrated.markers %>%
    group_by(cluster) %>%
    top_n(n = 20, wt = avg_log2FC) -> top20

DoHeatmap(pancreas.integrated, features = top20$gene) + NoLegend()
DimPlot(pancreas.integrated, reduction = "umap", group.by = "orig.ident")
#DimPlot(pancreas.integrated, reduction = "ref.umap",label = TRUE, 
#    label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")

VlnPlot(pancreas.integrated, features = c('Adgre1','Cd68','Csfr1','Itgam','Cd14','Mafb'))


write.csv(top20,'/home/wus/for_seurat/8_2/c_t_all_top20_gene_for_marker_8_25.csv')
write.csv(pancreas.integrated.markers,'/home/wus/for_seurat/8_2/c_t_all_markers_gene_8_25.csv')
save(top20,pancreas.integrated.markers,pancreas.integrated,file='/home/wus/c_t_all_8_23-all-feature.RData')

#画图（截取到annotation）
rm(list=ls())
load(file='/home/wus/c_t_all_8_23-all-feature.RData')
table(pancreas.integrated$orig.ident)
table(pancreas.integrated$seurat_clusters)
table(pancreas.integrated$sampleID)
FeaturePlot(pancreas.integrated, features = c('Adgre1','Cd68','Csfr1','Itgam','Cd14','Mafb'),col=c('grey','orange','red'),ncol=3)
###表达高低分
sce=pancreas.integrated
#mac_s = sce[,sce@meta.data$seurat_clusters %in% c(2)]
#DimPlot(mac_s, reduction = "umap")
#highCells=colnames(subset(x = mac_s, subset = S100a8 > 1, slot = 'counts'))
#highCells=colnames(subset(x = mac_s, subset = S100a9 > 0,slot = 'counts'))
highCells=colnames(subset(x = sce, subset = Adgre1 > 0, slot = 'counts'))
highORlow=ifelse(colnames(mac_s) %in% highCells,'high','low')
table(highORlow)
mac_s@meta.data$highORlow=highORlow
markers <- FindMarkers(mac_s, ident.1 = "high", 
                       group.by = 'highORlow', 
                       subset.ident = "2")
head(x = markers)




##### t 一起去分
rm(list=ls())
load(file='/home/wus/c_t_all_8_23-all-feature.RData')

sce=pancreas.integrated
cd4_sce1 = sce[,sce@meta.data$seurat_clusters %in% c(10)]
pancreas.list <- SplitObject(cd4_sce1, split.by = "orig.ident")
pancreas.list <- pancreas.list[c("C",'T')]


reference.list <- pancreas.list[c("C",'T')]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)


DefaultAssay(pancreas.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
pa1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "orig.ident")
pa2 <- DimPlot(pancreas.integrated, reduction = "umap", label = TRUE, repel = TRUE) + 
    NoLegend()
pa1 + pa2


###
pancreas.integrated <- FindNeighbors(pancreas.integrated, dims = 1:30)
pancreas.integrated <- FindClusters(pancreas.integrated, resolution = 0.5 )
DimPlot(pancreas.integrated, reduction = 'umap', group.by = 'seurat_clusters',
           label = TRUE, pt.size = 0.5) + NoLegend()
pancreas.integrated.markers<- FindAllMarkers(pancreas.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pancreas.integrated.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)


pancreas.integrated.markers %>%
    group_by(cluster) %>%
    top_n(n = 20, wt = avg_log2FC) -> top20

DoHeatmap(pancreas.integrated, features = top20$gene) + NoLegend()
DimPlot(pancreas.integrated, reduction = 'umap', group.by = 'seurat_clusters',
           label = TRUE, pt.size = 0.5) + NoLegend()

#write.csv(top20,'/home/wus/for_seurat/8_2/T_t-cell_merge_top20_gene_for_marker_8_25.csv')
#write.csv(pancreas.integrated.markers,'/home/wus/for_seurat/8_2/T_t-cell_merge_markers_gene_8_25.csv')

#VlnPlot(pancreas.integrated, features = c('Cd3e','Cd4','Cd8a','Foxp3','Mki67','Mzb1','Ncr1','Ncam1'))
#VlnPlot(pancreas.integrated, features = c('Sell','Ccr7','Sosc3','Myc','Lef1','Tshz2','Mal','Nosip','Tcf7'))
#VlnPlot(pancreas.integrated, features = c('Cd3e','Cd4','Cd8a','Cd20','Cd19','Nkg7'))
VlnPlot(pancreas.integrated,features = c('Mafb','Adgre1','Cd68','Itgam'))
#选C画图
 C_t_cell=subset(x=pancreas.integrated,subset=orig.ident=='C')
 DimPlot(C_t_cell, reduction = 'umap', group.by = 'seurat_clusters',
           label = TRUE, pt.size = 0.5) + NoLegend()
dim(C_t_cell)
table(C_t_cell@meta.data$seurat_clusters)
Cell_Num=data.frame(table(C_t_cell@meta.data$seurat_clusters))
save(C_t_cell,file = '/home/wus/for_seurat/8_2/C_t_cell_8_26.Rdata')

 T_t_cell=subset(x=pancreas.integrated,subset=orig.ident=='T')
 DimPlot(T_t_cell, reduction = 'umap', group.by = 'seurat_clusters',
           label = TRUE, pt.size = 0.5) + NoLegend()
 dim(T_t_cell)
table(T_t_cell@meta.data$seurat_clusters)
Cell_Num=data.frame(table(T_t_cell@meta.data$seurat_clusters))
save(T_t_cell,file = '/home/wus/for_seurat/8_2/T_t_cell_8_26.Rdata')




##### macs 统一

rm(list=ls())
load(file='/home/wus/c_t_all_8_23-all-feature.RData')

sce=pancreas.integrated
cd4_sce1 = sce[,sce@meta.data$seurat_clusters %in% c(2,5,12,16)]
pancreas.list <- SplitObject(cd4_sce1, split.by = "orig.ident")
pancreas.list <- pancreas.list[c("C",'T')]


reference.list <- pancreas.list[c("C",'T')]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)


DefaultAssay(pancreas.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
pa1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "orig.ident")
pa2 <- DimPlot(pancreas.integrated, reduction = "umap", label = TRUE, repel = TRUE) + 
    NoLegend()
pa1 + pa2
save(pancreas.integrated,pancreas.anchors,file = '/home/wus/for_seurat/8_2/CT_macs__step1_8_23.Rdata')

###
rm(list=ls())
load(file = '/home/wus/for_seurat/8_2/CT_macs__step1_8_23.Rdata')
#pancreas.integrated <- FindNeighbors(pancreas.integrated, dims = 1:30)
#pancreas.integrated <- FindClusters(pancreas.integrated, resolution = 0.5 )
pancreas.integrated <- FindNeighbors(pancreas.integrated, dims = 1:15)
pancreas.integrated <- FindClusters(pancreas.integrated, resolution = 0.3 )
DimPlot(pancreas.integrated, reduction = 'umap', group.by = 'seurat_clusters',
           label = TRUE, pt.size = 0.5) + NoLegend()
pancreas.integrated.markers<- FindAllMarkers(pancreas.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pancreas.integrated.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)


pancreas.integrated.markers %>%
    group_by(cluster) %>%
    top_n(n = 20, wt = avg_log2FC) -> top20

#DoHeatmap(pancreas.integrated, features = top20$gene) + NoLegend()
DimPlot(pancreas.integrated, reduction = 'umap', group.by = 'seurat_clusters',
           label = TRUE, pt.size = 0.5) + NoLegend()

write.csv(top20,'/home/wus/for_seurat/8_2/CT_macs_merge_top20_gene_for_marker_8_25.csv')
write.csv(pancreas.integrated.markers,'/home/wus/for_seurat/8_2/CT_macs_merge_markers_gene_8_25.csv')

save(pancreas.integrated,pancreas.integrated.markers,file = '/home/wus/for_seurat/8_2/CT_macs_step2_8_28.Rdata')

VlnPlot(pancreas.integrated, features = c('H2-Aa','Itgam','H2-Ab1','Clec9a','Xcr1','Itgax','Sirpa','Fcgr1','Irf8','Flt3','Zbtb46','Batf3','Itgae'))

#C画图
 C_t_cell=subset(x=pancreas.integrated,subset=orig.ident=='C')
 DimPlot(C_t_cell, reduction = 'umap', group.by = 'seurat_clusters',
           label = TRUE, pt.size = 0.5) + NoLegend()
dim(C_t_cell)
table(C_t_cell@meta.data$seurat_clusters)
Cell_Num=data.frame(table(C_t_cell@meta.data$seurat_clusters))
VlnPlot(C_t_cell, features = c('Itgam','Ly6c2','Ly6g'))

save(C_t_cell,file = '/home/wus/for_seurat/8_2/C_macs_8_26.Rdata')

#M1 M2
load(file = '/home/wus/for_seurat/8_2/C_macs_8_26.Rdata')
sub.markers<-FindMarkers(C_t_cell,ident.1=c(0,1),ident.2='4')
sub.markers
# M neu
load(file = '/home/wus/for_seurat/8_2/C_macs_8_26.Rdata')
sub.markers<-FindMarkers(C_t_cell,ident.1=c(0,1,4),ident.2=c(2,3,7))
sub.markers




#T画图
 T_t_cell=subset(x=pancreas.integrated,subset=orig.ident=='T')
 DimPlot(T_t_cell, reduction = 'umap', group.by = 'seurat_clusters',
           label = TRUE, pt.size = 0.5) + NoLegend()
dim(T_t_cell)
table(T_t_cell@meta.data$seurat_clusters)
Cell_Num=data.frame(table(T_t_cell@meta.data$seurat_clusters))

save(T_t_cell,file = '/home/wus/for_seurat/8_2/T_macs_8_26.Rdata')









###只对2（mac）重分
rm(list=ls())
load(file='/home/wus/c_t_all_8_23-all-feature.RData')

sce=pancreas.integrated
cd4_sce1 = sce[,sce@meta.data$seurat_clusters %in% c(2)]
pancreas.list <- SplitObject(cd4_sce1, split.by = "orig.ident")
pancreas.list <- pancreas.list[c("C",'T')]


reference.list <- pancreas.list[c("C",'T')]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
DefaultAssay(pancreas.integrated) <- "integrated"

#是否应该重新找高表达基因FindVariableFeatures
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 15, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:15, verbose = FALSE)
pa1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "orig.ident")
pa2 <- DimPlot(pancreas.integrated, reduction = "umap", label = TRUE, repel = TRUE) + 
    NoLegend()
pa1 + pa2
save(pancreas.integrated,pancreas.anchors,file = '/home/wus/for_seurat/8_2/CT_macs__step1_9_1.Rdata')

###新算法
rm(list=ls())
load(file = '/home/wus/for_seurat/8_2/CT_macs__step1_9_1.Rdata')
#pancreas.integrated <- FindNeighbors(pancreas.integrated, dims = 1:30)
#pancreas.integrated <- FindClusters(pancreas.integrated, resolution = 0.5 )
pancreas.integrated <- FindNeighbors(pancreas.integrated, dims = 1:15)
reticulate::use_python('/home/wus/.local/share/r-miniconda/envs/r-reticulate/bin/python')
pancreas.integrated <- FindClusters(pancreas.integrated, resolution = 0.3 ,algorithm=4)
DimPlot(pancreas.integrated, reduction = 'umap', group.by = 'seurat_clusters',
           label = TRUE, pt.size = 0.5,split.by = "orig.ident",ncol = 3) + NoLegend()
pancreas.integrated.markers<- FindAllMarkers(pancreas.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pancreas.integrated.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)


pancreas.integrated.markers %>%
    group_by(cluster) %>%
    top_n(n = 20, wt = avg_log2FC) -> top20

#DoHeatmap(pancreas.integrated, features = top20$gene) + NoLegend()
DimPlot(pancreas.integrated, reduction = 'umap', group.by = 'seurat_clusters',
           label = TRUE, pt.size = 0.5) + NoLegend()

write.csv(top20,'/home/wus/for_seurat/8_2/CT_macs_merge_top20_gene_for_marker_9_2.csv')
write.csv(pancreas.integrated.markers,'/home/wus/for_seurat/8_2/CT_macs_merge_markers_gene_9_2.csv')
save(pancreas.integrated,pancreas.integrated.markers,file = '/home/wus/for_seurat/8_2/CT_macs_step2_9_2.Rdata')

###保守marker
DefaultAssay(pancreas.integrated) <- "RNA"
pancreas.integrated.conservedmarkers <- FindConservedMarkers(pancreas.integrated,
                              ident.1 = 4,
                              grouping.var = "orig.ident",
                              only.pos = TRUE,
                      logfc.threshold = 0.25)
pancreas.integrated.conservedmarkers
pancreas.integrated.conservedmarkers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)


pancreas.integrated.conservedmarkers %>%
    group_by(cluster) %>%
    top_n(n = 20, wt = avg_log2FC) -> ctop20

#DoHeatmap(pancreas.integrated, features = top20$gene) + NoLegend()
DimPlot(pancreas.integrated, reduction = 'umap', group.by = 'seurat_clusters',
           label = TRUE, pt.size = 0.5) + NoLegend()

#write.csv(ctop20,'/home/wus/for_seurat/8_2/CT_macs_merge_ctop20_gene_for_marker_9_2.csv')
#write.csv(pancreas.integrated.conservedmarkers,'/home/wus/for_seurat/8_2/CT_macs_merge_conservedmarkers_gene_9_2.csv')
save(pancreas.integrated.conservedmarkers,file = '/home/wus/for_seurat/8_2/CT_macs_step2_9_2.Rdata')


###计数每组
VlnPlot(pancreas.integrated, features = c('H2-Aa','Itgam','H2-Ab1','Clec9a','Xcr1','Itgax','Sirpa','Fcgr1','Irf8','Flt3','Zbtb46','Batf3','Itgae'))
sce=pancreas.integrated
table(sce$seurat_clusters,sce$orig.ident )

#C T 两组比较
Idents(sce) = paste0(sce$orig.ident )
table(Idents(sce))
FindMarkers(sce[,sce$seurat_clusters==1],ident.1 = 'C',ident.2 = 'T')

FindMarkers(sce[,sce$seurat_clusters==1],ident.1 = 'T',ident.2 = 'C')

#C画图
 C_t_cell=subset(x=pancreas.integrated,subset=orig.ident=='C')
 DimPlot(C_t_cell, reduction = 'umap', group.by = 'seurat_clusters',
           label = TRUE, pt.size = 0.5) + NoLegend()
dim(C_t_cell)
table(C_t_cell@meta.data$seurat_clusters)
Cell_Num=data.frame(table(C_t_cell@meta.data$seurat_clusters))
VlnPlot(C_t_cell, features = c('Itgam','Ly6c2','Ly6g'))

save(C_t_cell,file = '/home/wus/for_seurat/8_2/C_macs_8_26.Rdata')

#M1 M2
load(file = '/home/wus/for_seurat/8_2/C_macs_8_26.Rdata')
sub.markers<-FindMarkers(C_t_cell,ident.1=c(0,1),ident.2='4')
sub.markers
# M neu
load(file = '/home/wus/for_seurat/8_2/C_macs_8_26.Rdata')
sub.markers<-FindMarkers(C_t_cell,ident.1=c(0,1,4),ident.2=c(2,3,7))
sub.markers




#T画图
 T_t_cell=subset(x=pancreas.integrated,subset=orig.ident=='T')
 DimPlot(T_t_cell, reduction = 'umap', group.by = 'seurat_clusters',
           label = TRUE, pt.size = 0.5) + NoLegend()
dim(T_t_cell)
table(T_t_cell@meta.data$seurat_clusters)
Cell_Num=data.frame(table(T_t_cell@meta.data$seurat_clusters))

save(T_t_cell,file = '/home/wus/for_seurat/8_2/T_macs_8_26.Rdata')





