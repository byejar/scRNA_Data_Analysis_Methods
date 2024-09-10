library(SummarizedExperiment)
library(TCGAbiolinks)
#存储控制台输出记录
sink("/home/wus/GDCdata/TCGAbiolinks_2-28_output.txt",append=FALSE,split = FALSE)
#版本信息
TCGAbiolinks::getGDCInfo()

###在默认位置存放文件：/home/wus/GDCdata/


query <- GDCquery(project = "TCGA-OV",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts")#GDC在2022.3.28更新了，workflow.type只能选择STAR - Counts，没有其他选项，所有FPKM等都一并被下载，在后续assay（ i=）中选择。
####查看数据内容
#dim(getResults(query))
#getResults(query)[1:4,]
#getResults(query)$sample_type
##样本类型
table(getResults(query)$sample_type)

##具体编号
sample_id<-getResults(query)$sample.submitter_id
table(substr(sample_id,14,16))



GDCdownload(query)
##需要给出保存文件的文件名，默认生成的文件名太长
data <- GDCprepare(query, save = T)

assay(data)

query_info = getResults(query)
TP = TCGAquery_SampleTypes(query_info$sample.submitter_id,"TP")
NT = TCGAquery_SampleTypes(query_info$sample.submitter_id,"NT")
query <- GDCquery(project = c("TCGA-OV"),
                  legacy = FALSE, #default(GDC harmonized database)
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts",
                  barcode = c(TP, NT))
dim(getResults(query))

data <- GDCprepare(query, save = T)

geneexp <- assay(data,i = "unstranded")#tpm_unstrand fpkm_unstrand

geneexp <- assay(data,i = "fpkm_unstrand")#tpm_unstrand fpkm_unstrand
geneexp[1:4,1:4]


samplesTP <- TCGAquery_SampleTypes(colnames(geneexp), typesample = c("NT"))
PGA5 <- geneexp[c("ENSG00000000003.15 "),samplesTP]
View(PGA5)



grep("clinical_", names(data), value = T)

Pair_sample = TCGAquery_MatchedCoupledSampleTypes(query_info$sample.submitter_id,c("NT","TP"))
query <- GDCquery(project = c("TCGA-OV"),
                  legacy = FALSE, #default(GDC harmonized database)
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts",
                  barcode = Pair_sample)
dim(getResults(query))


query <- GDCquery(project = c("TCGA-OV"),
                  legacy = FALSE, #default(GDC harmonized database)
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "STAR - Counts",
                  barcode = c( NT))
dim(getResults(query))


###############
####随机选取PGA5基因为例探索单个基因表达的情况对患者生存的影响
# 1、取特定基因（这里以PGA5基因为例）在癌症样本中的表达
#  2、先读取基因表达矩阵
#  3、提取特定基因在癌症样本中的表达，选出的样本都是肿瘤样本
samplesTP <- TCGAquery_SampleTypes(colnames(dataFilt_STAD_final), typesample = c("TP"))
PGA5 <- dataFilt_STAD_final[c("PGA5"),samplesTP]
View(PGA5)


################# 使用保存好的数据
load("TCGA-mRNA/TCGA-READ_mRNA.Rdata")

se <- data

这个se就是你的对象，含有coldata, rowdata， meta-data，以及最重要的assay，共有6个assay

类型是：RangedSummarizedExperiment
维度：60660行，177列
6个assay以及它们的名字
表达矩阵的行名
行信息（也就是基因信息），比如gene id,gene name,gene type
表达矩阵的列名（也就是样本名）
列信息，也就是样本信息，比如生存时间、生存状态这些

每个assay你可以理解为一个表达矩阵，
————————————————
版权声明：本文为CSDN博主「医学和生信笔记」的原创文章，遵循CC 4.0 BY-SA版权协议，转载请附上原文出处链接及本声明。
原文链接：https://blog.csdn.net/Ayue0616/article/details/126452387


############


#colData(geneexp)[1:4,]

#rowData(geneexp)[1:6,]



rowData(data)$gene_name
colData(data)$barcode



clinical <- GDCquery_clinic(project= 'TCGA-OV',type = "clinical")
clinical[1:4,]



d<-colnames(colData(data))
c<-colnames(clinical)
intersect(d,c)





############  TCGA 实际使用的下载代码
library(SummarizedExperiment)
library(TCGAbiolinks)
#存储控制台输出记录
#sink("/home/wus/GDCdata/TCGAbiolinks_2-28_output.txt",append=FALSE,split = FALSE)
#版本信息
TCGAbiolinks::getGDCInfo()

####下载数据类型
projects <- getGDCprojects()

library(dplyr)
projects <- projects %>% 
as.data.frame() %>% 
select(project_id,tumor) %>% 
  filter(grepl(pattern="TCGA",project_id))

#### 循环  
#dir.create("/home/wus/2023-2-23_008_TCGAbiolinks_download/TCGA_RNA_data/")
for (i in 1:nrow(projects)) {
## 0.运行信息
  print(paste0("Downloading number ",i,",project name: ",projects$project_id[i]))
## 1.查询信息
  query.exp = GDCquery(project = projects$project_id[i], 
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "STAR - Counts")

##输出信息
print(table(getResults(query.exp)$sample_type))
sample_id<-getResults(query.exp)$sample.submitter_id
print(table(substr(sample_id,14,16)))
## 2.正式下载
  GDCdownload(query.exp)
## 3.同一癌种数据合并
  pre.exp = GDCprepare(query = query.exp,
  						save = T,
                   		save.filename = paste0("/home/wus/2023-2-23_008_TCGAbiolinks_download/TCGA_RNA_data/",projects$project_id[i],".rda"))
## 4.提取表达量数据
  #countsdata = SummarizedExperiment::assay(pre.exp)
## 5.保存数据
  #save(countsdata,file = paste0("/home/wus/2023-2-23_008_TCGAbiolinks_download/TCGA_RNA_data/",projects$project_id[i],".Rdata"))
}

######最后两个文件下载失败重新下载
############  TCGA 实际使用的下载代码
library(SummarizedExperiment)
library(TCGAbiolinks)
#存储控制台输出记录
#sink("/home/wus/GDCdata/TCGAbiolinks_2-28_output.txt",append=FALSE,split = FALSE)
#版本信息
TCGAbiolinks::getGDCInfo()

####下载数据类型
projects <- getGDCprojects()

library(dplyr)
projects <- projects %>% 
as.data.frame() %>% 
select(project_id,tumor) %>% 
  filter(grepl(pattern="TCGA",project_id))

projects <-projects[32:33,]


#### 循环  
#dir.create("/home/wus/2023-2-23_008_TCGAbiolinks_download/TCGA_RNA_data/")
for (i in 1:nrow(projects)) {
## 0.运行信息
  print(paste0("Downloading number ",i,",project name: ",projects$project_id[i]))
## 1.查询信息
  query.exp = GDCquery(project = projects$project_id[i], 
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "STAR - Counts")

##输出信息
print(table(getResults(query.exp)$sample_type))
sample_id<-getResults(query.exp)$sample.submitter_id
print(table(substr(sample_id,14,16)))
## 2.正式下载
  GDCdownload(query.exp,
  			method='api',
  			files.per.chunk=6)
## 3.同一癌种数据合并
  pre.exp = GDCprepare(query = query.exp,
  						save = T,
                   		save.filename = paste0("/home/wus/2023-2-23_008_TCGAbiolinks_download/TCGA_RNA_data/",projects$project_id[i],".rda"))
## 4.提取表达量数据
  #countsdata = SummarizedExperiment::assay(pre.exp)
## 5.保存数据
  #save(countsdata,file = paste0("/home/wus/2023-2-23_008_TCGAbiolinks_download/TCGA_RNA_data/",projects$project_id[i],".Rdata"))
}

################ 统计信息
for (i in 1:nrow(projects)) {
## 0.运行信息
  print(paste0(",project name: ",projects$project_id[i]))
## 1.查询信息
  query.exp = GDCquery(project = projects$project_id[i], 
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "STAR - Counts")

##输出信息
print(table(getResults(query.exp)$sample_type))
sample_id<-getResults(query.exp)$sample.submitter_id
print(table(substr(sample_id,14,16)))
}

############## 使用下载好的数据
#两种方式，或者使用rda，或者query但是不down，第二种虽然可以但是还是麻烦，就选择使用rda的

load("/home/wus/2023-2-23_008_TCGAbiolinks_download/TCGA_RNA_data/TCGA-BLCA.rda")
ls()
data
assayNames(data)
geneexp <- assay(data,i = "unstranded")
samplesTP <- TCGAquery_SampleTypes(colnames(geneexp), typesample = c("TP"))
PGA5 <- geneexp['ENSG00000168883.20',samplesTP]
View(PGA5)

###转换出ID
gene_id<-rowData(data)[rowData(data)$gene_name=="USP39",'gene_id']


samplesNT <- TCGAquery_SampleTypes(colnames(geneexp), typesample = c("NT"))
PGA5 <- geneexp['ENSG00000168883.20',samplesNT]
View(PGA5)

Pair_sample = TCGAquery_MatchedCoupledSampleTypes(colnames(geneexp),c("NT","TP"))
PGA5 <- geneexp['ENSG00000168883.20',Pair_sample]
View(PGA5)
typeof(PGA5)
class(PGA5)
class(geneexp)

pd<-as.data.frame(PGA5)
class(pd)

gd<-as.data.frame(geneexp)
class(gd)
rownames(gd)
colnames(gd)
PGA5 <- gd['ENSG00000168883.20',Pair_sample]



### 只有一个样本
gd<-as.data.frame(geneexp)
PGA5 <- gd['ENSG00000168883.20','TCGA-GC-A3WC-11A-11R-A22U-07',drop=FALSE]
View(PGA5)
typeof(PGA5)
class(PGA5)
class(geneexp)
tPGA5<-t(PGA5)

pd<-as.data.frame(PGA5)
class(pd)

colnames(PGA5)
#####

####################患者信息

class(colData(data))

patients<-colData(data)[Pair_sample,c('days_to_death','days_to_last_follow_up')]
typeof(patients)

patients<-colData(data)[Pair_sample,c('barcode','days_to_death')]
typeof(patients)
pad<-as.data.frame(patients)
class(pad)



dd<-as.data.frame(colData(data))
class(dd)
colnames(dd)

patients<-dd[Pair_sample,'days_to_death']

class(patients)

patients<-dd[Pair_sample,c('days_to_last_follow_up'),drop=FALSE]
patients
class(patients)


df3<-merge(tPGA5, patients, by = "row.names", all = F)