
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
rm(list=ls())
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


load(file='/home/wus/c_t_mac_9_12.RData')
Idents(immune.combined.sct)
immune_combine=immune.combined.sct
Idents(immune_combine) <- "seurat_clusters"
immune_combine <- RenameIdents(immune_combine, `0` = "cancer cell", `1` = "cancer cell", `2` = "macrophage", 
    `3` = "cancer cell", `4` = "neurtophil", `5` = "fibroblast", `6` = "T cell", `7` = "7", `8` = "cancer cell", `9` = "DC", `10` = "B cell", `11` = "11")

DimPlot(immune_combine, reduction = "umap", #label = TRUE,
    repel = TRUE,pt.size = 1)

#去掉7 11 
mac.cells <- subset(immune_combine , idents = c("cancer cell","macrophage","fibroblast","neurtophil", "T cell","DC","B cell"))

DimPlot(mac.cells, reduction = "umap", label = TRUE,
    repel = TRUE)
DimPlot(mac.cells, reduction = "umap", #label = TRUE,
    repel = TRUE,pt.size = 1)

png(file = "/home/wus/scRNA_12_7.png",width = 1200,height = 600)
DimPlot(mac.cells, reduction = "umap", #label = TRUE,
    repel = TRUE,pt.size = 0.01,split.by = "orig.ident")
dev.off()

png(file = "/home/wus/scRNA_12_7-label.png",width = 1200,height = 600)
DimPlot(mac.cells, reduction = "umap", label = TRUE,
    repel = TRUE,pt.size = 0.01,split.by = "orig.ident")
dev.off()

###########抽样
mac.cells.C=subset(mac.cells,orig.ident=='C')
mac.cells.T=subset(mac.cells,orig.ident=='T')

DimPlot(subset(mac.cells.C, downsample = 10000), reduction = "umap",# label = TRUE,
    repel = TRUE,pt.size = 0.01,split.by = "orig.ident")

DimPlot(subset(mac.cells.T, downsample = 10000), reduction = "umap",# label = TRUE,
    repel = TRUE,pt.size = 0.01,split.by = "orig.ident")

#分C T 组的marker气泡图
markers.to.plot <-  c("Mertk", "Csf1r", "S100a9", "S100a8",'Dcn','Cd3d','Cd3e','H2-Aa','H2-Ab1','Cd19','Cd79a' )
DotPlot(mac.cells, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, 
    split.by = "orig.ident") + RotatedAxis()


###计数
sce=mac.cells

table(sce@meta.data$orig.ident)

head(sce[[]])
sce[["idents"]] <- Idents(object = sce)
head(sce[[]])
table(sce@meta.data$idents)

sce$group=ifelse(grepl("C",sce@meta.data$orig.ident),"C","T")
table(sce@meta.data$group)


source("home/wus/custom_seurat_functions.R")
p1 = plot.clusters.group(data = sce,clusters =  "idents", 
                         xlab = "idents", log = FALSE,
                         group = "orig.ident",legend.title = "Sample",
                         widths = c(3,1),color = 2)
p1

##两组各自


data=sce
count_table <- table(data@meta.data[,'idents'], data@meta.data[,'orig.ident'])

T_group<-as.data.frame(count_table[,2])
T_group <- cbind(rownames(T_group),T_group) 
names(T_group)<-c('type',  'Num')
T_group$group='T'

C_group<-as.data.frame(count_table[,1])
C_group <- cbind(rownames(C_group),C_group) 
names(C_group)<-c('type',  'Num')
C_group$group='C'
#C_group

df<-rbind(C_group,T_group)
df %>%
  group_by(group) %>%
  mutate(percent = round(Num/sum(Num),4)) ->df2
df2
ggplot(df2,aes(x=group,y=percent,fill=type))+
  geom_bar(stat="identity")+
  geom_text(aes(label=percent),position=position_stack(vjust=0.5))

df2$output<-paste(as.character(df2$percent),'(',as.character(df2$Num),')')
df2
ggplot(df2,aes(x=group,y=percent,fill=type))+
  geom_bar(stat="identity")+
  geom_text(aes(label=output),position=position_stack(vjust=0.5))

##也可以这么做
source("home/wus/custom_seurat_functions.R")
p1 = plot.clusters.group(data = sce,clusters =  "orig.ident", 
                         xlab = "idents", log = FALSE,
                         group = "idents",legend.title = "Sample",
                         widths = c(3,1),color = 2)
p1
ggsave("origident_celltype_percent.pdf",units = "cm",width = 27,height = 13)