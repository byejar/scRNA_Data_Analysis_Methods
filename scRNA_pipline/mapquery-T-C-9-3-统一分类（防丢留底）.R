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
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
#执行集成分析
immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30, verbose = FALSE)
immune.combined.sct <- FindNeighbors(immune.combined.sct, reduction = "pca", dims = 1:30)
immune.combined.sct <- FindClusters(immune.combined.sct, resolution = 0.3)
#并排可视化这两个条件，使用参数来显示按聚类着色的每个条件。split.by
DimPlot(immune.combined.sct, reduction = "umap", split.by = "orig.ident")

save(immune.anchors,immune.combined.sct,file='/home/wus/c_t_all_9_3.RData')

#识别相同基因，相同cluster，不同条件下的差异表达基因（基础）
load(file='/home/wus/c_t_all_9_3.RData')
DimPlot(immune.combined.sct, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
    repel = TRUE)

immune.combined.sct$celltype.stim <- paste(immune.combined.sct$seurat_clusters, immune.combined.sct$orig.ident,
    sep = "_")
Idents(immune.combined.sct) <- "celltype.stim"

immune.combined.sct <- PrepSCTFindMarkers(immune.combined.sct)
head(immune.combined.sct)
b.interferon.response <- FindMarkers(immune.combined.sct, assay = "SCT", ident.1 = "2_C", ident.2 = "2_T",
    verbose = FALSE)
head(b.interferon.response, n = 15)


###识别不同条件下的差异表达基因（散点图）
mac.cells <- subset(immune.combined.sct, idents = "2")
Idents(mac.cells) <- "orig.ident"
head(mac.cells)
mac.cells <- NormalizeData(object = mac.cells, assay = "RNA")
avg.mac.cells <- log1p(AverageExpression(mac.cells, verbose = FALSE)$RNA)
head(avg.mac.cells)

avg.mac.cells1=as.data.frame(avg.mac.cells)
class(avg.mac.cells1)
avg.mac.cells1$gene<-rownames(avg.mac.cells1)


genes.to.label = c("Fn1", "Saa3", "Spp1")
p1 <- ggplot(avg.mac.cells1, aes(C, T)) + geom_point() + ggtitle("macrophage")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)

plot_grid(p1)

#画图相同基因，相同cluster，在不同条件下的表达（feature和气泡图）
DefaultAssay(immune.combined.sct) <- "SCT"
FeaturePlot(immune.combined.sct, features = c("Mertk", "Csf1r", "Cd86"), split.by = "orig.ident", max.cutoff = 3,
    cols = c("grey", 'orange',"red"))


plots <- VlnPlot(immune.combined.sct, features = c("Fn1", "Saa3", "Spp1"), split.by = "orig.ident",
    group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

#流式
plots <- VlnPlot(immune.combined.sct, features = c("Mrc1", "Arg1", "H2-Ab1",'H2-Aa'), split.by = "orig.ident",
    group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

#不同展现形式（失败）
plots <- VlnPlot(immune.combined.sct, features = c("Mrc1", "Arg1"), split.by = "orig.ident",
    group.by = "seurat_clusters", pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)

#识别在各种条件下保守的规范细胞类型标记基因，鉴定出与刺激条件无关的保守标记的基因
Idents(immune.combined.sct) <- "seurat_clusters"
nk.markers <- FindConservedMarkers(immune.combined.sct, assay = "SCT", ident.1 = "2", grouping.var = "orig.ident",
    verbose = FALSE)
head(nk.markers)

save(immune.combined.sct,file='/home/wus/c_t_all_step_2_9_3.RData')


###对macrophage进行分群
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

write.csv(top20,'/home/wus/for_seurat/8_2/CT_macs_merge_top20_gene_for_marker_9_3.csv')
write.csv(pancreas.integrated.markers,'/home/wus/for_seurat/8_2/CT_macs_merge_markers_gene_9_3.csv')




#识别在各种条件下保守的规范细胞类型标记基因，鉴定出与刺激条件无关的保守标记的基因
Idents(mac.cells) <- "seurat_clusters"
nk.markers <- FindConservedMarkers(mac.cells, assay = "SCT", ident.1 = "2", grouping.var = "orig.ident",
    verbose = FALSE)
head(nk.markers)



###
#识别相同基因，相同cluster，不同条件下的差异表达基因（基础）

mac.cells$celltype.stim <- paste(mac.cells$seurat_clusters, mac.cells$orig.ident,
    sep = "_")
Idents(mac.cells) <- "celltype.stim"

mac.cells <- PrepSCTFindMarkers(mac.cells)
head(mac.cells)
b.interferon.response <- FindMarkers(mac.cells, assay = "SCT", ident.1 = "2_C", ident.2 = "2_T",
    verbose = FALSE)
head(b.interferon.response, n = 15)