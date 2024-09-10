###add score
rm(list=ls())
load(file='/home/wus/c_t_all_step_2_9_3.RData')
library(Seurat)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(patchwork)

DimPlot(immune.combined.sct,
        reduction = "umap",
        label = TRUE,
        label.size = 5)


DefaultAssay(immune.combined.sct) <- "SCT"
cd4_sce1<- subset(immune.combined.sct, idents = "6")

DimPlot(cd4_sce1,
        reduction = "umap",
        label = TRUE,
        label.size = 5)

FeaturePlot(cd4_sce1, 
            reduction = "umap", 
            features = c('Cd3e', 'Cd4', 'Cd8a','Nkg7','Mki67','Gzma'), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE
            )        

setwd('/home/wus/2023-5-9_012_PTEN_T-activate/')
# T cell activation 652
f=read.table('GO0042110.tsv',header=T)
typeof(f['name'])
#cd_features <- list(c('Adrm1','Ahsa1','C1galt1c1','Cct6b','Cd37','Cd3d','Cd3e','Cd3g','Cd69','Cd8a','Cetn3','Cse1l','Gemin6','Gnly','Gpt2','Gzma','Gzmh','Gzmk','Il2rb','Lck','Mpzl1','Nkg7','Pik3ip1','Ptrh2','Timm13','Zap70'))
cd_features <-f['name']


# T cell activation 500+
f=read.table('GO0042110-v2.tsv',header=T)
typeof(f['name'])
#cd_features <- list(c('Adrm1','Ahsa1','C1galt1c1','Cct6b','Cd37','Cd3d','Cd3e','Cd3g','Cd69','Cd8a','Cetn3','Cse1l','Gemin6','Gnly','Gpt2','Gzma','Gzmh','Gzmk','Il2rb','Lck','Mpzl1','Nkg7','Pik3ip1','Ptrh2','Timm13','Zap70'))
cd_features <-f['name']


f=read.table('GO0002286.tsv',header=T)
typeof(f['name'])
cd_features <-f['name']



f=read.table('GO0002418.tsv')
typeof(f['V2'])
cd_features <-f['V2']



cd_features <- list(c('Ahr','Il4i1','Ceacam1','Muc4','Havcr2','Klre1','Il4i1b','Tgfb1'))

cd_features <- list( c('Xcl1','Pvr','Raet1b','Mr1','Klrk1','H60a','Crtam',
'Cd160','Gsdme','Ulbp1','Cd274','Slc22a13','Hrg','Fbxo38',
'Cd226','Cd24a','Cd40lg','Hspd1','Il12a','Il12b','Nectin2'))


cd_features <- list( c('Pvr','Klrk1','Crtam','Cd226','Il12a','Il12b','Nectin2'))

cd_features <- list( c('Mr1','Slc22a13','Fbxo38','Cd24a','Hspd1'))


cd4_sce1 <- AddModuleScore(
  object = cd4_sce1,
  features = cd_features,
  ctrl = 5,
  name = 'Activated_CD8_T_cell'
)
head(x = cd4_sce1[])

colnames(cd4_sce1@meta.data)
VlnPlot(cd4_sce1,feature='Activated_CD8_T_cell1',split.by = "orig.ident")
x = cd4_sce1@meta.data
write.csv(x,'/home/wus/for_seurat/9_6/tcell_C_T_addmoudlescore.csv')
#加p值

library(ggpubr)
x=read.csv('/home/wus/for_seurat/9_6/tcell_C_T_addmoudlescore.csv')
#ggviolin(x,x='orig.ident',fill='orig.ident',y='M2_Features1')

p=ggviolin(x,x='orig.ident',fill='orig.ident',y='Activated_CD8_T_cell1',add = "boxplot", add.params = list(fill="white"))

compare_means(Activated_CD8_T_cell1~orig.ident, data = x, method = "wilcox.test")

my_comparisons <- list(c("C", "T"))
p+stat_compare_means(comparisons = my_comparisons,label = "p")+#不同组间的比较, 
stat_compare_means(label.y = 0.5)



#M1
cd_features <- list( c("Lcn2","Tnfaip2",'Lyz2','Fth1','Il1r1','Nos2','Tlr2','Tlr4',
 	'Il1b','Il18','Il12b','Il6','Cxcl2','Ifng',
 	'Cd38','Tnf','Socs3','Ptgs2','Nfkbiz','Lrf5','Gpr18',
 	'Fpr2','Cxcl10','Azin1','Cd68',
 	'Ccl5','Irf1','Irf5','Il1a','Ido1','Il12a','Il23a'))



cd4_sce1 <- AddModuleScore(
  object = cd4_sce1,
  features = cd_features,
  ctrl = 5,
  name = 'M1_Features'
)
head(x = cd4_sce1[])

colnames(cd4_sce1@meta.data)
VlnPlot(cd4_sce1,feature='M1_Features1',split.by = "orig.ident")
x = cd4_sce1@meta.data
write.csv(x,'/home/wus/for_seurat/9_6/C_T_M1_addmoudlescore.csv')

#加p值
library(ggpubr)
x=read.csv('/home/wus/for_seurat/9_6/C_T_M1_addmoudlescore.csv')
ggviolin(x,x='orig.ident',color='black',y='M1_Features1',width=0.5,add = c("boxplot",'jitter') , add.params = list(color="grey",size=1)  )#)#color=边框颜色 不写fill就不填充

ggviolin(x,x='orig.ident',color='black',y='M1_Features1',width=0.5,add = c('jitter'), add.params = list(color="orig.ident",size=4,alpha=0.5)  )

p=ggviolin(x,x='orig.ident',fill='orig.ident',y='M1_Features1',add = "boxplot", add.params = list(fill="white"))

compare_means(M1_Features1~orig.ident, data = x, method = "wilcox.test")


my_comparisons <- list(c("C", "T"))
p+stat_compare_means(comparisons = my_comparisons, label = "p.format",paired=FALSE)+#不同组间的比较
stat_compare_means(label.y = 1.5)



###ggstatsplot画图

library(ggstatsplot)
rt=read.csv('/home/wus/for_seurat/9_6/C_T_M1_addmoudlescore.csv')

pdf("/home/wus/scRNA_for_fig/M1.pdf",width = 6,height = 5) 
#png("/home/wus/scRNA_for_fig/M1.png",width = 1200,height = 1000,units = "px")
#tiff("/home/wus/scRNA_for_fig/M1.tiff",width = 6,height = 5,res=300) 
ggbetweenstats(data = rt,x='orig.ident',y='M1_Features1',title = "M1 expression",messages = FALSE)
dev.off() 


pdf("/home/wus/scRNA_for_fig/M2.pdf",width = 6,height = 5) 
ggbetweenstats(data = rt,x='orig.ident',y='M2_Features1',title = "M2 expression",messages = FALSE)
dev.off() 