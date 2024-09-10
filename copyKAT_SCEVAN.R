#scRNA







################################## 4-16 循环每个样本，将上皮细胞恶性和非恶性写入整合文件
#scRNA
library(Seurat)
library(ggplot2)
library(copykat)
library(tidyverse)
load('/home/wus/2023-3-21_011_AA-signature/GSE199252/scRNA.RData')
setwd('/home/wus/inferCNV/4-16_copyKAT/')
##取出T
scRNA.T=subset(scRNA,group == "T")
pancreas.list <- SplitObject(scRNA.T, split.by = 'orig.ident')
library(future)
#改变内存使用大小限制
options(future.globals.maxSize = 40000 * 1024^2)


for (i in 1:length(pancreas.list)) {
    pancreas.list[[i]] 
    T1<-pancreas.list[[i]]

epi.cells <- subset(T1, idents = c('epi_or_tumor'))


counts <- as.matrix(epi.cells@assays$RNA@counts)



copykat.test <- copykat(rawmat=counts, 
  id.type="S", ngene.chr=5, win.size=25, 
  KS.cut=0.1, sam.name=paste('T',i,sep=''), distance="euclidean",
   norm.cell.names="",output.seg="FLASE", 
   plot.genes="TRUE", genome="hg20",n.cores=40)

pred.test <- data.frame(copykat.test$prediction)
pred.test <- pred.test[-which(pred.test$copykat.pred=="not.defined"),]  ##remove undefined cells
CNA.test <- data.frame(copykat.test$CNAmat)

    }

############################# 单独测试
#scRNA@meta.data

#levels(scRNA)



# mallignant <- read.table("/home/wus/inferCNV/4-16_copyKAT/T2_copykat_prediction.txt",header=1)
# cell.names=mallignant$cell.names

# Idents(scRNA, cells = cell.names) <- mallignant$copykat.pred

########################################### 批量修改idents
for (i in 1:length(pancreas.list)) {
    file_name=paste("/home/wus/inferCNV/4-16_copyKAT/T",i,'_copykat_prediction.txt',sep='')
    mallignant <- read.table(file_name,header=1)
    cell.names=mallignant$cell.names
    Idents(scRNA, cells = cell.names) <- mallignant$copykat.pred
    }
DimPlot(scRNA,label = T)

scRNA_T<-scRNA

save(scRNA_T,file='/home/wus/2023-3-21_011_AA-signature/GSE199252/4_16_scRNA.Rdata')


###############纯肿瘤样本
########################################### 批量修改idents
for (i in 1:length(pancreas.list)) {
    file_name=paste("/home/wus/inferCNV/4-16_copyKAT/T",i,'_copykat_prediction.txt',sep='')
    mallignant <- read.table(file_name,header=1)
    cell.names=mallignant$cell.names
    Idents(scRNA.T, cells = cell.names) <- mallignant$copykat.pred
    }
DimPlot(scRNA.T,label = T)



save(scRNA.T,file='/home/wus/2023-3-21_011_AA-signature/GSE199252/4_16_scRNA.T.Rdata')


################################################### SCEVAN
library(Seurat)
library(SCEVAN)

counts <- Read10X(data.dir = '/home/wus/sc/clear_cellranger/M1015/M1015_TF2103ZJCR4_T11/filtered_feature_bc_matrix')
T1_org <- CreateSeuratObject(counts,project='t1', min.cells = 3, min.features = 200)

pbmc<-T1_org
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap", label=T)
T1_org<-pbmc



counts <- as.matrix(T1_org@assays$RNA@counts)


results <- SCEVAN::pipelineCNA(counts, sample = "T1_org", par_cores = 20, SUBCLONES = TRUE, plotTree = TRUE)

T1_org <- AddMetaData(T1_org, metadata = results)
head(T1_org)

p1 <- DimPlot(T1_org,label = T) #group.by = 'seurat_clusters'
p2 <- DimPlot(T1_org, group.by = "class") + scale_color_manual(values = c("blue", "gray",'red'))
pc <- p1 + p2
pc



###########
# 10-10
###########
#scRNA
library(Seurat)
library(ggplot2)
library(copykat)
library(tidyverse)

l<-load('/home/wus/021-wus-USP39-mTOR/10-9_SMC_Cohort_scRNA_data/SMC_wholetumor.RData')
sce<-mepub
Idents(sce)<-'Cell type'
SMC_Epi<-subset(sce,idents='Epithelial')
DimPlot(SMC_Epi, reduction = "umap", label=T)
# library(future)
# #改变内存使用大小限制
# options(future.globals.maxSize = 40000 * 1024^2)
epi.cells <- SMC_Epi
counts <- as.matrix(epi.cells@assays$RNA@counts)
copykat.test <- copykat(rawmat=counts, 
  id.type="S", ngene.chr=5, win.size=25, 
  KS.cut=0.1, sam.name='scm', distance="euclidean",
   norm.cell.names="",output.seg="FLASE", 
   plot.genes="TRUE", genome="hg20",n.cores=40)

pred.test <- data.frame(copykat.test$prediction)
pred.test <- pred.test[-which(pred.test$copykat.pred=="not.defined"),]  ##remove undefined cells
CNA.test <- data.frame(copykat.test$CNAmat)

# scm_copykat_prediction.txt是上一步生成的结果
setwd('/home/wus/021-wus-USP39-mTOR/10-9_SMC_Cohort_scRNA_data/')
mallignant <- read.delim("scm_copykat_prediction.txt", row.names = 1)
# 把细胞的良恶性信息加入metadata
scRNA <- AddMetaData(epi.cells, metadata = mallignant)
p1 <- DimPlot(scRNA, label = T)
p2 <- DimPlot(scRNA, group.by = "copykat.pred") #+ scale_color_manual(values = c("red", "gray"))
p3 <- DimPlot(scRNA, split.by = "copykat.pred")
pc <- p1 + p2
# ggsave("pred_mallignant.pdf", pc, width = 12, height = 5)
scRNA$USP39_express = ifelse(scRNA@assays$RNA@counts['USP39',]>mean(scRNA@assays$RNA@counts['USP39',]),'USP39_H','USP39_L')
table(scRNA$'copykat.pred',scRNA$USP39_express)

#########全部细胞
setwd('/home/wus/021-wus-USP39-mTOR/10-9_SMC_Cohort_scRNA_data/')
sce<-mepub
counts <- as.matrix(sce@assays$RNA@counts)
copykat.test <- copykat(rawmat=counts, 
  id.type="S", ngene.chr=5, win.size=25, 
  KS.cut=0.1, sam.name='scm_all', distance="euclidean",
   norm.cell.names="",output.seg="FLASE", 
   plot.genes="TRUE", genome="hg20",n.cores=40)

# pred.test <- data.frame(copykat.test$prediction)
# pred.test <- pred.test[-which(pred.test$copykat.pred=="not.defined"),]  ##remove undefined cells
# CNA.test <- data.frame(copykat.test$CNAmat)

mallignant <- read.delim("scm_all_copykat_prediction.txt", row.names = 1)
# 把细胞的良恶性信息加入metadata
scRNA <- AddMetaData(sce, metadata = mallignant)
p1 <- DimPlot(scRNA, label = T)
p2 <- DimPlot(scRNA, group.by = "copykat.pred") #+ scale_color_manual(values = c("red", "gray"))
p3 <- DimPlot(scRNA, split.by = "copykat.pred")
pc <- p1 + p2

scRNA$USP39_express = ifelse(scRNA@assays$RNA@counts['USP39',]>mean(scRNA@assays$RNA@counts['USP39',]),'USP39_H','USP39_L')
table(scRNA$'copykat.pred',scRNA$USP39_express)
Epi<-subset(scRNA,idents='Epithelial')
table(Epi$'copykat.pred',Epi$USP39_express)

##
load("/home/wus/12-17_005_GEO-scRNA/GSE178341_/12_27/GSE178341.Rdata")
# DimPlot(sce, reduction = 'umap',raster=FALSE,label =TRUE)
EpiT<-subset(sce,subset=ClusterMidway=='EpiT')
# DimPlot(Epi, reduction = "umap", label=T)
library(future)
#改变内存使用大小限制
options(future.globals.maxSize = 80000 * 1024^2)
epi.cells <- EpiT
counts <- as.matrix(epi.cells@assays$RNA@counts)
copykat.test <- copykat(rawmat=counts, 
  id.type="S", ngene.chr=5, win.size=25, 
  KS.cut=0.1, sam.name='crc', distance="euclidean",
   norm.cell.names="",output.seg="FLASE", 
   plot.genes="TRUE", genome="hg20",n.cores=40)

pred.test <- data.frame(copykat.test$prediction)
pred.test <- pred.test[-which(pred.test$copykat.pred=="not.defined"),]  ##remove undefined cells
CNA.test <- data.frame(copykat.test$CNAmat)


#### 换个方法
################################################### SCEVAN
library(Seurat)
library(SCEVAN)

load("/home/wus/12-17_005_GEO-scRNA/GSE178341_/12_27/GSE178341.Rdata")
EpiT<-subset(sce,subset=ClusterMidway=='EpiT')
epi.cells <- EpiT
counts <- as.matrix(epi.cells@assays$RNA@counts)

setwd('/home/wus/021-wus-USP39-mTOR/')
results <- SCEVAN::pipelineCNA(counts, sample = "epi.cells", par_cores = 10, SUBCLONES = TRUE, plotTree = TRUE)

epi.cells <- AddMetaData(epi.cells, metadata = results)
head(epi.cells)

p1 <- DimPlot(epi.cells,label = T) #group.by = 'seurat_clusters'
p2 <- DimPlot(epi.cells, group.by = "class") + scale_color_manual(values = c("blue", "gray",'red'))
pc <- p1 + p2
pc


##########
### 10-13 Harvard_Cohort_scRNA_data
##########
library(Seurat)
library(ggplot2)
library(copykat)
library(tidyverse)
setwd('/home/wus/021-wus-USP39-mTOR/10-13_Harvard_Cohort_scRNA_data/')
load("/home/wus/12-17_005_GEO-scRNA/GSE178341_/12_27/GSE178341.Rdata")



pancreas.list <- SplitObject(sce, split.by = 'orig.ident')


for (i in 1:length(pancreas.list)) {
    pancreas.list[[i]] 
    T1<-pancreas.list[[i]]

    # epi.cells <- subset(T1, idents = c('epi_or_tumor'))


    # counts <- as.matrix(epi.cells@assays$RNA@counts)
    counts <- as.matrix(T1@assays$RNA@counts)



copykat.test <- copykat(rawmat=counts, 
  id.type="S", ngene.chr=5, win.size=25, 
  KS.cut=0.1, sam.name=paste('Harvard_',i,sep=''), distance="euclidean",
   norm.cell.names="",output.seg="FLASE", 
   plot.genes="TRUE", genome="hg20",n.cores=40)

# pred.test <- data.frame(copykat.test$prediction)
# pred.test <- pred.test[-which(pred.test$copykat.pred=="not.defined"),]  ##remove undefined cells
# CNA.test <- data.frame(copykat.test$CNAmat)

    }

for (i in 1:length(pancreas.list)) {
    file_name=paste("/home/wus/021-wus-USP39-mTOR/10-13_Harvard_Cohort_scRNA_data/Harvard_",i,'_copykat_prediction.txt',sep='')
    mallignant <- read.table(file_name,header=1)
    cell.names=mallignant$cell.names
    Idents(scRNA, cells = cell.names) <- mallignant$copykat.pred
    }