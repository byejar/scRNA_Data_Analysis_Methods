library(ggpubr)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(dplyr)




########### 循环
csv_file='/home/wus/2023-2-23_008_TCGAbiolinks_download/box/merge-v2.txt'

projects <- getGDCprojects()

library(dplyr)
projects <- projects %>% 
as.data.frame() %>% 
select(project_id,tumor) %>% 
  filter(grepl(pattern="TCGA",project_id))
###加载变量
projects

for (i in 1:nrow(projects)) {
	## 0.运行信息
  print(paste0(" number ",i,",project name: ",projects$project_id[i]))

###
#fig_file=paste('/home/wus/2023-2-23_008_TCGAbiolinks_download/pairline/',projects$project_id[i],"_paired_points.png")
####### 配对数据构建
TCGA_rda<-paste0("/home/wus/2023-2-23_008_TCGAbiolinks_download/TCGA_RNA_data/",projects$project_id[i],".rda")

load(TCGA_rda)


geneexp <- assay(data,i = "tpm_unstrand")
#TCGAquery_SampleTypes(barcode, typesample)
##  colnames()=  vector
#### 选择是A的冰冻样本 
barcode <- colnames(geneexp) %>% 
	as.data.frame() %>%
  filter(grepl(pattern="-[0-9]+A-",colnames(geneexp)))
#取出一列就是向量类型
b<-barcode[,1]
#### 选取样本
 samples=TCGAquery_SampleTypes(b, c("NT","TP"))

####
samplesNT=TCGAquery_SampleTypes(b, c("NT"))

print(paste0('samplesNT:  ',samplesNT,'  length: ',length(samplesNT)))

#R=TCGAquery_MatchedCoupledSampleTypes(b,c("NT","TP"))

if(length(samplesNT)==0){
  print(paste0(projects$project_id[i],'have no matched sample'))
  next;}
else{
  

#### 基因名转换
gene_id<-rowData(data)[rowData(data)$gene_name=="USP39",'gene_id']

####
gd<-as.data.frame(geneexp)
samples_express <- gd[gene_id,samples,drop=FALSE]
##转置并加一列type
samples_express<-t(samples_express)
##样本信息
sample_type_information<-colData(data)[ samples  ,c('sample_type','patient','project_id'),drop=FALSE]

###合并
df3<-merge(samples_express, sample_type_information, by = "row.names", all = F)

df3

NT_mean=mean(df3[df3$sample_type == 'Primary Tumor','ENSG00000168883.20'])
TP_mean=mean(df3[df3$sample_type == 'Solid Tissue Normal','ENSG00000168883.20'])

df3$fc=NT_mean/TP_mean

write.table( df3, csv_file,row.names=F,append = TRUE)
}}


#######

library(ggplot2)
data <- read.table("/home/wus/2023-2-23_008_TCGAbiolinks_download/box/merge-v2.txt",header = T)

data2 <- subset(data, patient != 'patient')

data2 <- subset(data2, project_id != 'TCGA-SKCM')

write.csv(data2,'/home/wus/2023-2-23_008_TCGAbiolinks_download/box/data2-v2.csv')
#ggplot(data2,aes(SAMPLE,USP39,fill=COLOR)) + 
#        geom_boxplot(outlier.shape = NA)

#ggboxplot(data2,x="SAMPLE",y="USP39",
#          width = 0.6,fill="COLOR")


#################
#data3 <- read.csv("data3.csv",header = T)

#ggplot(data3,aes(SAMPLE,USP39,fill=COLOR)) + 
#        geom_boxplot(outlier.shape = NA)

#ggboxplot(data3,x="SAMPLE",y="USP39",
#          width = 0.6,fill="COLOR")


#################
data3 <- read.csv("/home/wus/2023-2-23_008_TCGAbiolinks_download/box/data2-v2.csv",header = T)



#p <-  ggboxplot(data3, x="project_id", y="ENSG00000168883.20", color = "sample_type",bxp.errorbar =TRUE,outlier.shape = NA)#,order=)
# p+stat_compare_means(aes(group=sample_type),label = "p.format",method = "t.test")+ theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))


data4<-data3[order(data3$fc, decreasing = TRUE), ]

p <-  ggboxplot(data4, x="project_id", y="ENSG00000168883.20", color = "sample_type",bxp.errorbar =TRUE,outlier.shape = NA)#,order=)
 p+stat_compare_means(aes(group=sample_type),label = "p.format",method = "t.test")+ 
 theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))



 #####
data4$project_id <- as.factor(data4$project_id)
head(data4,n=3)

data4<-data3[order(data3$fc, decreasing = TRUE), ]

stat.test <- data4 %>% 
  group_by(fc) %>%
  t_test(ENSG00000168883.20 ~ sample_type) %>%
  adjust_pvalue(method = "bonferroni") %>% 
  add_significance("p.adj")
stat.test





 bxp <- ggboxplot(
    data4, x="project_id", y="ENSG00000168883.20", color = "sample_type",
    palette=c("blue","pink"))

# add p value 并且确定xy的位置
stat.test <- stat.test %>% 
  add_xy_position(x='project_id',dodge = 1)
  #add_y_position()                                                                                                                                                                                                                                                                  


bxp + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, .5)))


bxp + stat_pvalue_manual(
  stat.test,  label = "{p.adj}{p.adj.signif}", 
  tip.length = 0, hide.ns = TRUE
  )




###########3-5 -测试给NT加factor 并
data4<-data3[order(data3$fc, decreasing = TRUE), ]
data4$fc <- as.factor(data4$fc)
data4$fc <- factor(data4$fc, levels=rev(levels(data4$fc)))

data4$sample_type <-as.factor(data4$sample_type)
data4$sample_type <- factor(data4$sample_type, levels=rev(levels(data4$sample_type)))

head(data4,n=3)

stat.test <- data4 %>% 
  group_by(fc) %>%
  t_test(ENSG00000168883.20 ~ sample_type) %>%
  adjust_pvalue(method = "bonferroni") %>% 
  add_significance("p.adj")

stat.test
stat.test<-stat.test[order(stat.test$fc, decreasing = FALSE), ]
stat.test


#data4<-data3[order(data3$fc, decreasing = TRUE), ]
 bxp <- ggboxplot(
    data4, x="project_id", y="ENSG00000168883.20", color = "sample_type",bxp.errorbar =TRUE,outlier.shape = NA)
#bxp
# add p value 并且确定xy的位置
stat.test <- stat.test %>% 
  add_xy_position(x='fc',dodge = 1)
  #add_y_position() 
as.data.frame(stat.test)


bxp + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
  ) + 
scale_y_continuous(expand = expansion(mult = c(0, .5)))+ 
 theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))

########################  使用自带的方式计算p（阳性对照）

data4<-data3[order(data3$fc, decreasing = TRUE), ]

p <-  ggboxplot(data4, x="project_id", y="ENSG00000168883.20", color = "sample_type",bxp.errorbar =TRUE,outlier.shape = NA)#,order=)
 p+stat_compare_means(aes(group=sample_type),label = "p.format",method = "t.test",angle=45,label.y.npc='center')+ 
 theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))


###########3-5 保存图片（ T N 版本）
data4<-data3[order(data3$fc, decreasing = TRUE), ]
data4$fc <- as.factor(data4$fc)
data4$fc <- factor(data4$fc, levels=rev(levels(data4$fc)))
head(data4,n=3)

stat.test <- data4 %>% 
  group_by(fc) %>%
  t_test(ENSG00000168883.20 ~ sample_type) %>%
  adjust_pvalue(method = "bonferroni") %>% 
  add_significance("p.adj")

stat.test
stat.test<-stat.test[order(stat.test$fc, decreasing = FALSE), ]
stat.test


#data4<-data3[order(data3$fc, decreasing = TRUE), ]
 bxp <- ggboxplot(
    data4, x="project_id", y="ENSG00000168883.20", color = "sample_type",bxp.errorbar =TRUE,outlier.shape = NA)
#bxp
# add p value 并且确定xy的位置
stat.test <- stat.test %>% 
  add_xy_position(x='fc',dodge = 1)
  #add_y_position() 
as.data.frame(stat.test)


fig_file_txt='/home/wus/2023-2-23_008_TCGAbiolinks_download/box/box_txt.png'


bxp + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
  ) + 
scale_y_continuous(expand = expansion(mult = c(0, .5)))+ 
 theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+ylab("USP39_TPM")

  ggsave(fig_file_txt,dpi=300,width=15,height=7)



fig_file='/home/wus/2023-2-23_008_TCGAbiolinks_download/box/box_no_txt.png'
bxp+ 
 theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+ylab("USP39_TPM")  + ylim(0,200)
ggsave(fig_file,dpi=300,width=15,height=7)

####改颜色
fig_file_txt_c='/home/wus/2023-2-23_008_TCGAbiolinks_download/box/box_txt-color.png'


bxp +
    scale_color_manual(values = c("#9a0000", "#0b5d0a"))+ stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
  ) + 
scale_y_continuous(expand = expansion(mult = c(0, .5)))+ 
 theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+ylab("USP39_TPM")

  ggsave(fig_file_txt_c,dpi=300,width=15,height=7)



fig_file_c='/home/wus/2023-2-23_008_TCGAbiolinks_download/box/box_no_txt-color.png'
bxp+
    scale_color_manual(values = c("#9a0000", "#0b5d0a"))+ 
 theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+ylab("USP39_TPM")  + ylim(0,200)
ggsave(fig_file_c,dpi=300,width=15,height=7)



############## N T版本
 
data4<-data3[order(data3$fc, decreasing = TRUE), ]
data4$fc <- as.factor(data4$fc)
data4$fc <- factor(data4$fc, levels=rev(levels(data4$fc)))


data4$sample_type <-as.factor(data4$sample_type)
data4$sample_type <- factor(data4$sample_type, levels=rev(levels(data4$sample_type)))


head(data4,n=3)

stat.test <- data4 %>% 
  group_by(fc) %>%
  t_test(ENSG00000168883.20 ~ sample_type) %>%
  adjust_pvalue(method = "bonferroni") %>% 
  add_significance("p.adj")

stat.test
stat.test<-stat.test[order(stat.test$fc, decreasing = FALSE), ]
stat.test


#data4<-data3[order(data3$fc, decreasing = TRUE), ]
 bxp <- ggboxplot(
    data4, x="project_id", y="ENSG00000168883.20", color = "sample_type",bxp.errorbar =TRUE,outlier.shape = NA)
#bxp
# add p value 并且确定xy的位置
stat.test <- stat.test %>% 
  add_xy_position(x='fc',dodge = 1)
  #add_y_position() 
as.data.frame(stat.test)


fig_file_txt='/home/wus/2023-2-23_008_TCGAbiolinks_download/box/box_txt-NT.png'


bxp + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
  ) + 
scale_y_continuous(expand = expansion(mult = c(0, .5)))+ 
 theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+ylab("USP39_TPM")

  ggsave(fig_file_txt,dpi=300,width=15,height=7)



fig_file='/home/wus/2023-2-23_008_TCGAbiolinks_download/box/box_no_txt-NT.png'
bxp+ 
 theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+ylab("USP39_TPM")  + ylim(0,200)
ggsave(fig_file,dpi=300,width=15,height=7)

####改颜色
fig_file_txt_c='/home/wus/2023-2-23_008_TCGAbiolinks_download/box/box_txt-color-NT.png'


bxp +
    scale_color_manual(values = c( "#0b5d0a","#9a0000"))+ stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
  ) + 
scale_y_continuous(expand = expansion(mult = c(0, .5)))+ 
 theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+ylab("USP39_TPM")

  ggsave(fig_file_txt_c,dpi=300,width=15,height=7)



fig_file_c='/home/wus/2023-2-23_008_TCGAbiolinks_download/box/box_no_txt-color-NT.png'
bxp+
    scale_color_manual(values = c( "#0b5d0a","#9a0000"))+ 
 theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+ylab("USP39_TPM")  + ylim(0,200)
ggsave(fig_file_c,dpi=300,width=15,height=7)


####### 改变长宽

fig_file_c9='/home/wus/2023-2-23_008_TCGAbiolinks_download/box/box_no_txt-color-NT-9.png'
bxp+
    scale_color_manual(values = c( "#0b5d0a","#9a0000"))+ 
 theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+ylab("USP39_TPM")  + ylim(0,200)
ggsave(fig_file_c9,dpi=300,width=13,height=7)


#######统计NT数量（每组的）,保存文本

 
data4<-data3[order(data3$fc, decreasing = TRUE), ]
data4$fc <- as.factor(data4$fc)
data4$fc <- factor(data4$fc, levels=rev(levels(data4$fc)))


data4$sample_type <-as.factor(data4$sample_type)
data4$sample_type <- factor(data4$sample_type, levels=rev(levels(data4$sample_type)))


head(data4,n=3)

stat.test <- data4 %>% 
  group_by(fc,project_id) %>%
  t_test(ENSG00000168883.20 ~ sample_type) %>%
  adjust_pvalue(method = "bonferroni") %>% 
  add_significance("p.adj")

stat.test
stat.test<-stat.test[order(stat.test$fc, decreasing = FALSE), ]
stat.test


#bxp
# add p value 并且确定xy的位置
stat.test <- stat.test %>% 
  add_xy_position(x='fc',dodge = 1)
  #add_y_position() 
s=as.data.frame(stat.test)
d4=as.data.frame(data4)

s

write.csv(s,'/home/wus/2023-2-23_008_TCGAbiolinks_download/box/stat_test.csv')
#write.csv(d4,'/home/wus/2023-2-23_008_TCGAbiolinks_download/box/data4.csv')




