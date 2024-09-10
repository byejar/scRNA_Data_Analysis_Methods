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
DefaultAssay(immune.combined.sct) <- "SCT"
cd4_sce1<- subset(immune.combined.sct, idents = "2")

cd_features <- list(c('Myc','Cd83','Mrc1','Arg1','Egr2','Ccr7','Chil4',
	'Pparg','Cd163','Clec7a','Il10','Il4','Irf4','Pdgfb','Stat6','Chil3',
	'Cd40','Csf1','Tlr7','Ccl13','Ccl17','Ccl18','Ccl22','Ccl24','Cd86',
	'Vegfa','Vegfb','Vegfc','Vegfd','Mmp9','Fn1','Egf','Lyve1','Mmp14',
	'Mmp19','Cd276','Fasl','Ctsa','Ctsb','Ctsc','Ctsd','Mgl2','Ear11',
	'Clec10a','Retnla','Ccl10'))



cd4_sce1 <- AddModuleScore(
  object = cd4_sce1,
  features = cd_features,
  ctrl = 5,
  name = 'M2_Features'
)
head(x = cd4_sce1[])

colnames(cd4_sce1@meta.data)
VlnPlot(cd4_sce1,feature='M2_Features1',split.by = "orig.ident")
x = cd4_sce1@meta.data
write.csv(x,'/home/wus/for_seurat/9_6/C_T_addmoudlescore.csv')
#加p值

library(ggpubr)
x=read.csv('/home/wus/for_seurat/9_6/C_T_addmoudlescore.csv')
#ggviolin(x,x='orig.ident',fill='orig.ident',y='M2_Features1')

p=ggviolin(x,x='orig.ident',fill='orig.ident',y='M2_Features1',add = "boxplot", add.params = list(fill="white"))

compare_means(M2_Features1~orig.ident, data = x, method = "wilcox.test")

#my_comparisons <- list(c("C", "T"))
p+stat_compare_means(comparisons = my_comparisons,label = "p")+#不同组间的比较, 
stat_compare_means(label.y = 1.5)



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