
rm(list=ls())
load(file='/home/wus/c_t_all_9_3.RData')
library(Seurat)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(patchwork)

#mac
DefaultAssay(immune.combined.sct) <- "SCT"
FeaturePlot(immune.combined.sct, features = c("Mertk", "Csf1r", "Cd86"),  max.cutoff = 3,
    cols = c("grey", 'orange',"red"))
#neu
FeaturePlot(immune.combined.sct, features = c("S100a9", "S100a8", "Csf3r"), max.cutoff = 3,
    cols = c("grey", 'orange',"red"))

#DC  'H2-Aa','Itgam','H2-Ab1'=dc  'Clec9a','Xcr1'=dc1  'Itgax','Sirpa'=cd2
# B
FeaturePlot(immune.combined.sct, features = c('H2-Aa','H2-Ab1','Cd19','Cd79a'),  max.cutoff = 3,
    cols = c("grey", 'orange',"red"))

#fib

FeaturePlot(immune.combined.sct, features =  c('Dcn'),  max.cutoff = 3,
    cols = c("grey", 'orange',"red"))

#NK（没有）
#DefaultAssay(immune.combined.sct) <- "SCT"
#FeaturePlot(immune.combined.sct, features = c('Cd3e','Cd4','Cd8a','Cd8b1','Klrb1c','Gzma'),  max.cutoff = 3,
#    cols = c("grey", 'orange',"red"))


###看allmarker
rm(list=ls())
library(Seurat)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(patchwork)
load(file='/home/wus/c_t_mac_9_6.RData')

pancreas.integrated.markers<- FindAllMarkers(immune.combined.sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pancreas.integrated.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)


pancreas.integrated.markers %>%
    group_by(cluster) %>%
    top_n(n = 20, wt = avg_log2FC) -> top20


write.csv(top20,'/home/wus/for_seurat/9_6/c_t_all_top20_gene_for_marker_9_12.csv')
write.csv(pancreas.integrated.markers,'/home/wus/for_seurat/9_6/c_t_all_markers_gene_9_12.csv')

save(pancreas.integrated.markers,immune.combined.sct,file='/home/wus/c_t_mac_9_12.RData')

##取出 6，10，9，4，2
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
     pt.size = 0, combine = FALSE,assay = "SCT")
wrap_plots(plots = plots, ncol = 1)


markers.to.plot <-  c("Mertk", "Csf1r", "S100a9", "S100a8",'Dcn','Cd3d','Cd3e','H2-Aa','H2-Ab1','Cd19','Cd79a' )

DotPlot(mac.cells, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, 
    split.by = "orig.ident") + RotatedAxis()


###


#识别macrophage总群不同条件下的差异表达基因（基础）
rm(list=ls())
load(file='/home/wus/c_t_all_9_3.RData')
DimPlot(immune.combined.sct, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
    repel = TRUE)

immune.combined.sct$celltype.stim <- paste(immune.combined.sct$seurat_clusters, immune.combined.sct$orig.ident,
    sep = "_")
Idents(immune.combined.sct) <- "celltype.stim"

immune.combined.sct <- PrepSCTFindMarkers(immune.combined.sct)
head(immune.combined.sct)
b.interferon.response <- FindMarkers(immune.combined.sct, assay = "SCT", ident.1 = "2_C", ident.2 = "2_T",
    verbose = FALSE,slot = 'counts',logfc.threshold =0,min.pct = 0)

head(b.interferon.response, n = 15)

write.csv(b.interferon.response,'/home/wus/for_seurat/9_6/c_t_mac_all_markers_gene_9_12.csv')

save(b.interferon.response,immune.combined.sct,file='/home/wus/c_t_mac_9_13.RData')
