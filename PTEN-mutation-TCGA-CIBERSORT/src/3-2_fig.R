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

#####1.去掉B，合并T N ，调用pairline函数画图；
###文件格式
###
###______sample_______ | ____USP39_____ | _____type_____
###                                           N or T



###单独测试
load("/home/wus/2023-2-23_008_TCGAbiolinks_download/TCGA_RNA_data/TCGA-CHOL.rda")
ls()
geneexp <- assay(data,i = "unstranded")
#TCGAquery_SampleTypes(barcode, typesample)
##  colnames()=  vector
#### 选择是A的冰冻样本 
barcode <- colnames(geneexp) %>% 
	as.data.frame() %>%
  filter(grepl(pattern="-[0-9]+A-",colnames(geneexp)))
#取出一列就是向量类型
b<-barcode[,1]
#### 选取样本
 samplesNT=TCGAquery_SampleTypes(b, 'NT')

Pair_sample=TCGAquery_MatchedCoupledSampleTypes(b,c("NT","TP"))

#### 基因名转换
gene_id<-rowData(data)[rowData(data)$gene_name=="USP39",'gene_id']

####
gd<-as.data.frame(geneexp)
Pair_sample_express <- gd[gene_id,Pair_sample,drop=FALSE]
##转置并加一列type
Pair_sample_express<-t(Pair_sample_express)
##样本信息
sample_type_information<-colData(data)[ Pair_sample  ,c('sample_type','patient'),drop=FALSE]

###合并
df3<-merge(Pair_sample_express, sample_type_information, by = "row.names", all = F)

#######画图
library(ggplot2)
      ggplot(df3,aes(sample_type,ENSG00000168883.20)) +
    geom_point(aes(color=sample_type),size=1.5) +
      geom_line(aes(group = patient),
              color="grey"
              )+ theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.y=element_blank())+
  labs(title = i)
    write.csv(dfu,csv_file)
  ggsave(fig_file)


###加载变量
projects

for (i in 1:nrow(projects)) {
## 0.运行信息
  print(paste0("Downloading number ",i,",project name: ",projects$project_id[i]))


















#####2.