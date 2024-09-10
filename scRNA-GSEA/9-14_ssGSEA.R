############## 使用下载好的数据
#两种方式，使用rda，或者query但是不down，第二种虽然可以但是还是麻烦，就选择使用rda的
library(SummarizedExperiment)
library(TCGAbiolinks)
library(dplyr)
# load("/home/wus/2023-2-23_008_TCGAbiolinks_download/TCGA_RNA_data/TCGA-COAD.rda")
# ls()
# data
# assayNames(data)
# geneexp <- assay(data,i = "fpkm_unstrand")
# samplesTP <- TCGAquery_SampleTypes(colnames(geneexp), typesample = c("TP"))


# USP39 <- geneexp['ENSG00000168883.20',samplesTP]

# View(USP39)

# ###转换出ID
# gene_id<-rowData(data)[rowData(data)$gene_name=="USP39",'gene_id']


# samplesNT <- TCGAquery_SampleTypes(colnames(geneexp), typesample = c("NT"))
# USP39 <- geneexp['ENSG00000168883.20',samplesNT]
# View(USP39)

# Pair_sample = TCGAquery_MatchedCoupledSampleTypes(colnames(geneexp),c("NT","TP"))
# USP39 <- geneexp['ENSG00000168883.20',Pair_sample]
# View(USP39)
# typeof(USP39)
# class(USP39)
# class(geneexp)

# pd<-as.data.frame(USP39)
# class(pd)

# gd<-as.data.frame(geneexp)
# class(gd)
# rownames(gd)
# colnames(gd)
# USP39 <- gd['ENSG00000168883.20',Pair_sample]



# ### 只有一个样本
# gd<-as.data.frame(geneexp)
# USP39 <- gd['ENSG00000168883.20','TCGA-GC-A3WC-11A-11R-A22U-07',drop=FALSE]
# View(USP39)
# typeof(USP39)
# class(USP39)
# class(geneexp)
# tUSP39<-t(USP39)

# pd<-as.data.frame(USP39)
# class(pd)

# colnames(USP39)
# #####

####################患者信息

# class(colData(data))

# patients<-colData(data)[Pair_sample,c('days_to_death','days_to_last_follow_up')]
# typeof(patients)

# patients<-colData(data)[Pair_sample,c('barcode','days_to_death')]
# typeof(patients)
# pad<-as.data.frame(patients)
# class(pad)



# dd<-as.data.frame(colData(data))
# class(dd)
# colnames(dd)

# patients<-dd[Pair_sample,'days_to_death']

# class(patients)

# patients<-dd[Pair_sample,c('days_to_last_follow_up'),drop=FALSE]
# patients
# class(patients)


# df3<-merge(tUSP39, patients, by = "row.names", all = F)

#####1.去掉B，合并T N ，调用pairline函数画图；
###文件格式
###
###______sample_______ | ____USP39_____ | _____type_____
###                                           N or T



###单独测试
load("/home/wus/2023-2-23_008_TCGAbiolinks_download/TCGA_RNA_data/TCGA-COAD.rda")
ls()
geneexp <- assay(data,i = "fpkm_unstrand")
#TCGAquery_SampleTypes(barcode, typesample)
##  colnames()=  vector
#### 选择是A的冰冻样本 
barcode <- colnames(geneexp) %>% 
	as.data.frame() %>%
  filter(grepl(pattern="-[0-9]+A-",colnames(geneexp)))
#取出一列就是向量类型
b<-barcode[,1]
#### 选取样本
# samplesNT=TCGAquery_SampleTypes(b, 'NT')
samplesTP=TCGAquery_SampleTypes(b, 'TP')
# Pair_sample=TCGAquery_MatchedCoupledSampleTypes(b,c("NT","TP"))

#### 基因名转换
gene_id<-rowData(data)[rowData(data)$gene_name=="USP39",'gene_id']

####
gd<-as.data.frame(geneexp)
Pair_sample_express <- gd[gene_id,samplesTP,drop=FALSE]
##转置并加一列type
Pair_sample_express<-t(Pair_sample_express)
##样本信息
sample_type_information<-colData(data)[ samplesTP  ,c('sample_type','patient','days_to_death','days_to_last_follow_up','vital_status'),drop=FALSE]

###合并
df3<-merge(Pair_sample_express, sample_type_information, by = "row.names", all = F)

write.csv(df3,'/home/wus/021-wus-USP39-mTOR/COAD_expr-patient.csv')

# #######画图
# library(ggplot2)
#       ggplot(df3,aes(sample_type,ENSG00000168883.20)) +
#     geom_point(aes(color=sample_type),size=1.5) +
#       geom_line(aes(group = patient),
#               color="grey"
#               )+ theme_bw()+
#   theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.y=element_blank())+
#   labs(title = i)
#     write.csv(dfu,csv_file)
#   ggsave(fig_file)


###加载变量
# projects

# for (i in 1:nrow(projects)) {
# ## 0.运行信息
#   print(paste0("Downloading number ",i,",project name: ",projects$project_id[i]))
#####2.
typeof(data)
rowData(data)
typeof(rowData(data))
rowData(data)$gene_name

duplicated(rowData(data)$gene_name)
uniq_express <- data[!duplicated(rowData(data)$gene_name),samplesTP]
uniq_gd<-as.data.frame(uniq_express)
rownames(uniq_gd)<-rowData(uniq_express)$gene_name

write.csv(uniq_gd,'/home/wus/021-wus-USP39-mTOR/COAD_gene_name.csv')

# conda activate irGSEA
library(GSVA)


geneSets <- GSEABase::getGmt('/home/wus/019_2023-07-20_zc_BMDM/202307BMDM-si/cell_line/10_analysis/HALLMARK_MTORC1_SIGNALING.v2023.1.Hs.gmt')

    mydata <- read.csv(file = '/home/wus/021-wus-USP39-mTOR/COAD_gene_name.csv',header=T,row.names=1,check.names=FALSE)
    # mdt=data.frame(t(mydata))
    # expr=as.matrix(mdt)
    # expr
    expr=as.matrix(mydata)
    gsva_matrix<- gsva(expr, 
                       geneSets,             
                       method='ssgsea',  
                       kcdf='Poisson',  #By default, kcdf="Gaussian" which is suitable when input expression values are continuous, such as microarray fluorescent units in logarithmic scale, RNA-seq log-CPMs, log-RPKMs or log-TPMs. When input expression values are integer counts, such as those derived from RNA-seq experiments, then this argument should be set to kcdf="Poisson".
                       abs.ranking=TRUE)

write.csv(gsva_matrix,'/home/wus/021-wus-USP39-mTOR/COAD_gsva_matrix.csv')

## excel合并


z <- read.csv(file ='/home/wus/021-wus-USP39-mTOR/COAD-MTOR-USP33-8-28.csv',header=T,row.names=1,check.names=FALSE)
path='/home/wus/021-wus-USP39-mTOR/'
#source('/home/wus/ge/pipline/cor_addtext.R')
source('/home/wus/021-wus-USP39-mTOR/cor_addtext.R')
cor_addtext(z,'z','HALLMARK_MTORC1_SIGNALING','USP39',path,'HALLMARK_MTORC1_SIGNALING','USP39')


#### 相关性图
# library(ggstatsplot)
# ggscatterstats(
#   data = z,
#   x = HALLMARK_MTORC1_SIGNALING,
#   y = USP39,
#   type = "p",
#   conf.level = 0.99,
#   messages = TRUE
# )

############### cox
#conda activate R4.0
library("survival") #生存分析
library("survminer") #结果可视化
z <- read.csv(file ='/home/wus/021-wus-USP39-mTOR/COAD-MTOR-USP33-8-28.csv',header=T,row.names=1,check.names=FALSE)
z$type<-as.character(z$type)
# res.cox <- coxph(Surv(os, state) ~ MTORC1_SIGNALING + USP39 + type, data = z)
# summary(res.cox)
# ggforest(res.cox, data = z, 
#          main = "Hazard ratio",
#          cpositions = c(0.10, 0.22, 0.4),
#          fontsize = 1.0) 

res.cox <- coxph(Surv(os, state) ~  type, data = z)
summary(res.cox)
ggforest(res.cox, data = z, 
         main = "Hazard ratio",
         cpositions = c(0.10, 0.22, 0.4),
         fontsize = 1.0) 


res.cox <- coxph(Surv(os, state) ~ MTORC1_SIGNALING + USP39, data = z)
summary(res.cox)
ggforest(res.cox, data = z, 
         main = "Hazard ratio",
         cpositions = c(0.10, 0.22, 0.4),
         fontsize = 1.0) 

# res.cox <- coxph(Surv(os, state) ~ HALLMARK_MTORC1_SIGNALING , data = z)
# summary(res.cox)
# ggforest(res.cox, data = z, 
#          main = "Hazard ratio",
#          cpositions = c(0.10, 0.22, 0.4),
#          fontsize = 1.0) 

# res.cox <- coxph(Surv(os, state) ~ USP39, data = z)
# summary(res.cox)
# ggforest(res.cox, data = z, 
#          main = "Hazard ratio",
#          cpositions = c(0.10, 0.22, 0.4),
#          fontsize = 1.0) 


fit <- survfit(Surv(os,state) ~  type, data = z)
ggsurvplot(fit, data = z)

fit <- survfit(Surv(os,state) ~  U, data = z)
ggsurvplot(fit, data = z)

fit <- survfit(Surv(os,state) ~  t3, data = z)
ggsurvplot(fit, data = z)



######################################### 
###                   1.ssGSEA  2.分为HH HL L  3.画图
#########################################

load("/home/wus/2023-2-23_008_TCGAbiolinks_download/TCGA_RNA_data/TCGA-BRCA.rda")
ls()
geneexp <- assay(data,i = "tpm_unstrand")
barcode <- colnames(geneexp) %>% 
  as.data.frame() %>%
  filter(grepl(pattern="-[0-9]+A-",colnames(geneexp)))
b<-barcode[,1]
#### 选取样本
samplesTP=TCGAquery_SampleTypes(b, 'TP')
#### 基因名转换
gene_id<-rowData(data)[rowData(data)$gene_name=="USP39",'gene_id']
####
gd<-as.data.frame(geneexp)
Pair_sample_express <- gd[gene_id,samplesTP,drop=FALSE]
##转置并加一列type
Pair_sample_express<-t(Pair_sample_express)
##样本信息
sample_type_information<-colData(data)[ samplesTP  ,c('sample_type','patient','days_to_death','days_to_last_follow_up','vital_status'),drop=FALSE]
###合并
df3<-merge(Pair_sample_express, sample_type_information, by = "row.names", all = F)


df2 <- df3
df2 <- df2[!is.na(df2$vital_status),]
# 10.2 将status表示患者结局，1表示删失，2表示死亡
df2[df2$vital_status=='Dead',]$vital_status <- 1
df2[df2$vital_status=='Alive',]$vital_status <- 0
df2$vital_status <- as.numeric(df2$vital_status)

df2$time <- df2$days_to_death
df2$time[which(is.na(df2$time))] <- df2$days_to_last_follow_up[which(is.na(df2$time))]
        
write.csv(df2,'/home/wus/021-wus-USP39-mTOR/BRCA_expr-patient.csv')

# #######画图
# library(ggplot2)
#       ggplot(df3,aes(sample_type,ENSG00000168883.20)) +
#     geom_point(aes(color=sample_type),size=1.5) +
#       geom_line(aes(group = patient),
#               color="grey"
#               )+ theme_bw()+
#   theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.y=element_blank())+
#   labs(title = i)
#     write.csv(dfu,csv_file)
#   ggsave(fig_file)


###加载变量
# projects

# for (i in 1:nrow(projects)) {
# ## 0.运行信息
#   print(paste0("Downloading number ",i,",project name: ",projects$project_id[i]))
#####2.
# typeof(data)
# rowData(data)
# typeof(rowData(data))
# rowData(data)$gene_name

# duplicated(rowData(data)$gene_name)
uniq_express <- data[!duplicated(rowData(data)$gene_name),samplesTP]

uniq_gd<-assay(uniq_express,i = "tpm_unstrand")
rownames(uniq_gd)<-rowData(uniq_express)$gene_name
write.csv(uniq_gd,'/home/wus/021-wus-USP39-mTOR/BRCA_gene_name.csv')

# conda activate irGSEA
library(GSVA)


geneSets <- GSEABase::getGmt('/home/wus/019_2023-07-20_zc_BMDM/202307BMDM-si/cell_line/10_analysis/HALLMARK_MTORC1_SIGNALING.v2023.1.Hs.gmt')
    mydata <- read.csv(file = '/home/wus/021-wus-USP39-mTOR/BRCA_gene_name.csv',header=T,row.names=1,check.names=FALSE)
    expr=as.matrix(mydata)
    gsva_matrix<- gsva(expr, 
                       geneSets,             
                       method='ssgsea',  
                       kcdf='Poisson',  #By default, kcdf="Gaussian" which is suitable when input expression values are continuous, such as microarray fluorescent units in logarithmic scale, RNA-seq log-CPMs, log-RPKMs or log-TPMs. When input expression values are integer counts, such as those derived from RNA-seq experiments, then this argument should be set to kcdf="Poisson".
                       abs.ranking=TRUE)

write.csv(gsva_matrix,'/home/wus/021-wus-USP39-mTOR/BRCA_gsva_matrix.csv')

## excel合并


z <- read.csv(file ='/home/wus/021-wus-USP39-mTOR/BRCA-MTOR-USP33-9-1.csv',header=T,row.names=1,check.names=FALSE)
path='/home/wus/021-wus-USP39-mTOR/'
#source('/home/wus/ge/pipline/cor_addtext.R')
source('/home/wus/021-wus-USP39-mTOR/cor_addtext.R')
cor_addtext(z,'z','HALLMARK_MTORC1_SIGNALING','ENSG00000168883.20',path,'HALLMARK_MTORC1_SIGNALING','ENSG00000168883.20')


#### 相关性图
# library(ggstatsplot)
# ggscatterstats(
#   data = z,
#   x = HALLMARK_MTORC1_SIGNALING,
#   y = USP39,
#   type = "p",
#   conf.level = 0.99,
#   messages = TRUE
# )

############### cox
#conda activate R4.0
library("survival") #生存分析
library("survminer") #结果可视化
z <- read.csv(file ='/home/wus/021-wus-USP39-mTOR/BRCA-MTOR-USP33-9-15.csv',header=T,row.names=1,check.names=FALSE)




fit <- survfit(Surv(time,vital_status) ~  type, data = z)
ggsurvplot(fit, data = z)

res.cox <- coxph(Surv(time,vital_status) ~  type, data = z)
summary(res.cox)
ggforest(res.cox, data = z, 
         main = "Hazard ratio",
         cpositions = c(0.10, 0.22, 0.4),
         fontsize = 1.0) 


z <- read.csv(file ='/home/wus/021-wus-USP39-mTOR/LUAD-MTOR-USP33-9-15.csv',header=T,row.names=1,check.names=FALSE)

fit <- survfit(Surv(time,vital_status) ~  type, data = z)
ggsurvplot(fit, data = z,risk.table = TRUE)

res.cox <- coxph(Surv(time,vital_status) ~  type, data = z)
summary(res.cox)
ggforest(res.cox, data = z, 
         main = "Hazard ratio",
         cpositions = c(0.10, 0.22, 0.4),
         fontsize = 1.0) 



#####################
##              9-3
#####################
# load("/home/wus/2023-2-23_008_TCGAbiolinks_download/TCGA_RNA_data/TCGA-LUAD.rda")
# ls()
# geneexp <- assay(data,i = "tpm_unstrand")
# barcode <- colnames(geneexp) %>% 
#   as.data.frame() %>%
#   filter(grepl(pattern="-[0-9]+A-",colnames(geneexp)))
# b<-barcode[,1]
# #### 选取样本
# samplesTP=TCGAquery_SampleTypes(b, 'TP')

# uniq_express <- data[!duplicated(rowData(data)$gene_name),samplesTP]

# uniq_gd<-assay(uniq_express,i = "unstranded")
# rownames(uniq_gd)<-rowData(uniq_express)$gene_name
# write.table(uniq_gd,'/home/wus/021-wus-USP39-mTOR/LUAD_gene_name-unstranded.gct',sep='\t')


### logstic 
# z <- read.csv(file ='/home/wus/021-wus-USP39-mTOR/LUAD-MTOR-USP33-9-1-v2.csv',header=T,row.names=1,check.names=FALSE)

# training_dataset <- z[1:400, ]
# wdbc_test <- z[401:521,]
# logmodel <- glm(vital_status~USP39 + MTORC1_SIGNALING, training_dataset, family = binomial(link = "logit"))
# train_pre<-predict.glm(logmodel, wdbc_test)# 预测模型

# library(pROC)
# roclong <- plot.roc(wdbc_test$vital_status, lty = 1, train_pre, grid = T,
#                     percent=TRUE,ci=TRUE,col="red",
#                     print.auc = T,print.auc.cex = 0.8,print.auc.x = 60,print.auc.y = 50,
#                     main = "逻辑回归的ROC曲线和AUC值")

#############################
# 9-10
#############################
library("survival") #生存分析
library("survminer") #结果可视化
z <- read.csv(file ='/home/wus/021-wus-USP39-mTOR/LUAD-MTOR-USP33-9-10.csv',header=T,row.names=1,check.names=FALSE)

fit <- survfit(Surv(time,vital_status) ~  S, data = z)
ggsurvplot(fit, data = z,risk.table = TRUE,pval =T)

fit <- survfit(Surv(time,vital_status) ~  type2, data = z)
ggsurvplot(fit, data = z,risk.table = TRUE)

###
z <- read.csv(file ='/home/wus/021-wus-USP39-mTOR/COAD-MTOR-USP33-9-10.csv',header=T,row.names=1,check.names=FALSE)

fit <- survfit(Surv(os,state) ~  S, data = z)
ggsurvplot(fit, data = z,risk.table = TRUE,pval =T)


############# 9-21
#conda activate R4.0
library("survival") #生存分析
library("survminer") #结果可视化
z <- read.csv(file ='/home/wus/021-wus-USP39-mTOR/BRCA-MTOR-USP33-9-21.csv',header=T,row.names=1,check.names=FALSE)




fit <- survfit(Surv(time_5,vital_status_5) ~  type, data = z)
ggsurvplot(fit, data = z,
         pval = TRUE,
         risk.table = TRUE)

res.cox <- coxph(Surv(time_5,vital_status_5) ~  type, data = z)
summary(res.cox)
ggforest(res.cox, data = z, 
         main = "Hazard ratio",
         cpositions = c(0.10, 0.22, 0.4),
         fontsize = 1.0) 
