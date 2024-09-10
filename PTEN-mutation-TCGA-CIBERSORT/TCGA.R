----------------
#超算
----------------
Y3W#uSrF

----------------
#使用下载好的数据
----------------


load("/dssg/home/acct-medwshuai/medwshuai/2024-3-25-TCGA/TCGA-BRCA.rda")
ls()
data
assayNames(data)
#geneexp <- assay(data,i = "unstranded")
geneexp <- assay(data,i = "tpm_unstrand")

BRCA<-geneexp['ENSG00000168883.20',]

write.csv(BRCA,'/dssg/home/acct-medwshuai/medwshuai/2024-3-25-TCGA/TCGA-BRCA.csv')

###############

load("/dssg/home/acct-medwshuai/medwshuai/2024-3-25-TCGA/TCGA-COAD.rda")
ls()
data
assayNames(data)
#geneexp <- assay(data,i = "unstranded")
geneexp <- assay(data,i = "tpm_unstrand")

COAD<-geneexp['ENSG00000168883.20',]

write.csv(COAD,'/dssg/home/acct-medwshuai/medwshuai/2024-3-25-TCGA/TCGA-COAD.csv')

###############

load("/dssg/home/acct-medwshuai/medwshuai/2024-3-25-TCGA/TCGA-LUAD.rda")
ls()
data
assayNames(data)
#geneexp <- assay(data,i = "unstranded")
geneexp <- assay(data,i = "tpm_unstrand")

LUAD<-geneexp['ENSG00000168883.20',]

write.csv(LUAD,'/dssg/home/acct-medwshuai/medwshuai/2024-3-25-TCGA/TCGA-LUAD.csv')

###############
library(dplyr)
library(SummarizedExperiment)
library(TCGAbiolinks)

load("/dssg/home/acct-medwshuai/medwshuai/2024-3-25-TCGA/TCGA-PAAD.rda")
ls()
data
assayNames(data)
#geneexp <- assay(data,i = "unstranded")
geneexp <- assay(data,i = "tpm_unstrand")


barcode <- colnames(geneexp) %>% 
  as.data.frame() %>%
  filter(grepl(pattern="-[0-9]+A-",colnames(geneexp)))

#取出一列就是向量类型
b<-barcode[,1]
#### 选取样本
samplesNT=TCGAquery_SampleTypes(b, 'NT')
samplesTP=TCGAquery_SampleTypes(b, 'TP')

PAAD_NT<-geneexp['ENSG00000168883.20',samplesNT]
PAAD_TP<-geneexp['ENSG00000168883.20',samplesTP]

write.csv(PAAD,'/dssg/home/acct-medwshuai/medwshuai/2024-3-25-TCGA/TCGA-PAAD.csv')


-------------------
# 循环
-------------------
library(dplyr)
library(SummarizedExperiment)
library(TCGAbiolinks)

project_id=c('BRCA','LUAD','PAAD','COAD')


for (i in 1:4) {
## 0.运行信息
  print(paste0(" number ",i,",project name: ",project_id[i]))

### 
#fig_file_T=paste('/home/wus/2023-2-23_008_TCGAbiolinks_download/survival/',projects$project_id[i],"_TCGAbiolinks.png")
#fig_file_S=paste('/home/wus/2023-2-23_008_TCGAbiolinks_download/survival/',projects$project_id[i],"_survival.png")
NT=paste('/dssg/home/acct-medwshuai/medwshuai/2024-3-25-TCGA/',project_id[i],"_NT.csv")
TP=paste('/dssg/home/acct-medwshuai/medwshuai/2024-3-25-TCGA/',project_id[i],"_TP.csv")
####### 配对数据构建
TCGA_rda<-paste0("/dssg/home/acct-medwshuai/medwshuai/2024-3-25-TCGA/TCGA-",project_id[i],".rda")

load(TCGA_rda)
ls()
data
# assayNames(data)
#geneexp <- assay(data,i = "unstranded")
geneexp <- assay(data,i = "tpm_unstrand")


barcode <- colnames(geneexp) %>% 
  as.data.frame() %>%
  filter(grepl(pattern="-[0-9]+A-",colnames(geneexp)))

#取出一列就是向量类型
b<-barcode[,1]
#### 选取样本
samplesNT=TCGAquery_SampleTypes(b, 'NT')
samplesTP=TCGAquery_SampleTypes(b, 'TP')

PAAD_NT<-geneexp['ENSG00000168883.20',samplesNT]
PAAD_TP<-geneexp['ENSG00000168883.20',samplesTP]

write.csv(PAAD_NT,NT)
write.csv(PAAD_TP,TP)
}

------------------------------