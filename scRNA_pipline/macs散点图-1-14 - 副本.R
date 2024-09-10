rm(list=ls())
load(file='/home/wus/c_t_mac_9_12.RData')
library(Seurat)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(patchwork)


###识别不同条件下的差异表达基因（带自己的注释，散点图）

mac.cells <- subset(immune.combined.sct, idents = "2")
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
#write.csv(p_top50,'/home/wus/for_seurat/9_6/C_T_mac_p_top50_9_6.csv')
#write.csv(n_top50,'/home/wus/for_seurat/9_6/C_T_mac_n_top50_9_6.csv')

write.csv(avg.mac_sub1.marker,'/home/wus/for_seurat/9_6/C_T_mac_1-14.csv')

####更改画图  无字版本
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

#pc <- LabelPoints(plot = pc, points = p_list, repel = TRUE)
#pc <- LabelPoints(plot = pc, points = n_list, repel = TRUE)
plot_grid(pc)
ggsave(pc,filename = "/home/wus/scRNA_for_fig/mac_org-v2.png",width = 6,height = 5,dpi=300)







####带标注版本

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

#pc <- LabelPoints(plot = pc, points = p_list, repel = TRUE)
#pc <- LabelPoints(plot = pc, points = n_list, repel = TRUE)
add_list=c('Il1b','Mrc1','Arg1')
pc <- LabelPoints(plot = pc, points = add_list, repel = TRUE,xnudge = 1,ynudge = 1)
plot_grid(pc)

ggsave(pc,filename = "/home/wus/scRNA_for_fig/mac_addtext-v2.png",width = 6,height = 5,dpi=300)