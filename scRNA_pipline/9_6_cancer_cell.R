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
b.interferon.response <- FindMarkers(immune.combined.sct, assay = "SCT", ident.1 = c("0_C","1_C","3_C","8_C"), ident.2 = c("0_T","1_T","3_T","8_T"),
    verbose = FALSE)
head(b.interferon.response, n = 15)

save(b.interferon.response,immune.combined.sct,file='/home/wus/c_t_cancer_cell_9_6.RData')

###zhushi 
rm(list=ls())
load(file='/home/wus/c_t_cancer_cell_9_6.RData')
Idents(immune.combined.sct)
immune_combine=immune.combined.sct
Idents(immune_combine) <- "seurat_clusters"
immune_combine <- RenameIdents(immune_combine, `0` = "cancer cell", `1` = "cancer cell", `2` = "macrophage", 
    `3` = "cancer cell", `4` = "4", `5` = "5", `6` = "6", `7` = "7", `8` = "cancer cell", `9` = "9", `10` = "10", `11` = "11")
DimPlot(immune_combine, reduction = "umap", label = TRUE,
    repel = TRUE)
#####

###识别不同条件下的差异表达基因（带自己的注释，散点图）

mac.cells <- subset(immune_combine, idents = "cancer cell")
Idents(mac.cells) <- "orig.ident"
head(mac.cells)
head(mac.cells)

#log1p=ln(1+x)
#mac.cells <- NormalizeData(object = mac.cells, assay = "RNA")
#AverageExpression(mac.cells, features = 'Saa3')
#FindMarkers(mac.cells, features = 'Saa3', ident.1 = "C", ident.2 = "T")
#log((0.1118246+1)/(1.225565+1))

#log((0.03254314+1)/(0.464605+1))
#avg.mac.cells <- log1p(AverageExpression(mac.cells, verbose = FALSE)$RNA)
#avg.mac.cells <- log2((AverageExpression(mac.cells, verbose = FALSE)$RNA)+1)
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
  scale_colour_manual(values=c('up'='#b4292c', 'nosig'="#c8c8ca",'down'='#5177b9'),name='up_or_down') + 
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

ggsave(pc,filename = "/home/wus/scRNA_for_fig/tumor-v2.png",width = 6,height = 5,dpi=300)
