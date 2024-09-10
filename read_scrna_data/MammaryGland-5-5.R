##scRNA
################ 全部时期
MammaryGland.Involution1_dge.txt
MammaryGland.Involution2_dge.txt.gz
MammaryGland.Lactation1_dge.txt.gz
MammaryGland.Lactation2_dge.txt.gz
MammaryGland.Pregnancy_dge.txt.gz
MammaryGland.Virgin1_dge.txt.gz
MammaryGland.Virgin2_dge.txt.gz
MammaryGland.Virgin3_dge.txt.gz
MammaryGland.Virgin4_dge.txt.gz


library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(cowplot)

set.seed(123)  #设置随机数种子，使结果可重复
##==合并数据集==##
##使用目录向量合并
setwd='/home/wus/2023-4-25_012_gmk_Mammary-Gland/500more_dge/'

dir = c('MammaryGland.Involution1_dge.txt',
        'MammaryGland.Involution2_dge.txt',
        'MammaryGland.Pregnancy_dge.txt',
        'MammaryGland.Virgin1_dge.txt',
        'MammaryGland.Virgin2_dge.txt',
        'MammaryGland.Virgin3_dge.txt',
        'MammaryGland.Virgin4_dge.txt')

samples_name = c('Involution1_dge','Involution2_dge',
	'Pregnancy_dge','Virgin1_dge','Virgin2_dge','Virgin3_dge','Virgin4_dge')

group_name=c('Involution','Involution',
	'Pregnancy','Virgin','Virgin','Virgin','Virgin')

scRNAlist <- list()

for(i in 1:length(dir)){
counts <- read.table(dir[i], sep = " ", header = TRUE)
print(dir[i])
a=length(colnames(counts))
print(a)


scRNAlist[[i]] <- CreateSeuratObject(counts,project=samples_name[i], min.cells = 3, min.features = 200)
print(scRNAlist[[i]])


scRNAlist[[i]]<-RenameCells(scRNAlist[[i]],add.cell.id=samples_name[i])

scRNAlist[[i]][['percent.MT']]<-PercentageFeatureSet(scRNAlist[[i]], pattern = '^mt-')
scRNAlist[[i]][['percent.RP']]<-PercentageFeatureSet(scRNAlist[[i]], pattern = '^Rp[sl]')
scRNAlist[[i]][['group']]=group_name[i]
scRNAlist[[i]][['sample']]=samples_name[i]
#HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
#HB.genes <- CaseMatch(HB.genes, rownames(scRNAlist[[i]])) 
#scRNAlist[[i]][["percent.HB"]]<-PercentageFeatureSet(scRNAlist[[i]], features=HB.genes) 

}
#使用merge函数讲10个seurat对象合并成一个seurat对象
scRNA2 <- merge(scRNAlist[[1]], y=scRNAlist[2:length(scRNAlist)])


#################################
VlnPlot(scRNA2, features = c("nFeature_RNA", "nCount_RNA", "percent.MT","percent.RP"), ncol = 4,group='sample',pt.size = 0)


#save(scRNA2,file='5_2_MammaryGland.Rdata')

#scRNA2[['group']]=c('Involution','Involution','Lactation','Lactation',
#	'Pregnancy','Virgin','Virgin','Virgin','Virgin')
head(scRNA2@meta.data)

table(scRNA2@meta.data$orig.ident)

#####检查 Lactation
#subset(x = scRNA2, subset = group == "Lactation")


#设置质控标准
maxGene=2000
minUMI=200
maxUMI=5000
pctMT=10
#数据质控并绘制小提琴图
scRNA <-subset(scRNA2,subset =nCount_RNA >= minUMI & nCount_RNA <= maxUMI &
                nFeature_RNA < maxGene & percent.MT < pctMT  )
#################################
VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.MT","percent.RP"), ncol = 4,group='sample',pt.size = 0)

scRNA2<-scRNA

pancreas.list <- SplitObject(scRNA2, split.by = 'group')
split_seurat <- pancreas.list [c('Involution','Pregnancy','Virgin')]

library(future)
plan("multiprocess", workers = 24)
#改变内存使用大小限制
options(future.globals.maxSize = 40000 * 1024^2)
for (i in 1:length(split_seurat)) {
    split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
    #split_seurat[[i]] <- CellCycleScoring(split_seurat[[i]], g2m.features=g2m_genes, s.features=s_genes)
    split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("percent.MT"))
    }
# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 
# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)
# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
# Integrate across conditions
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")

seurat_integrated <- RunPCA(object = seurat_integrated)

# Plot PCA
PCAPlot(seurat_integrated,
        split.by = "group")  

# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
           reduction = "pca")

# Plot UMAP                             
#DimPlot(seurat_integrated)  

DimPlot(seurat_integrated,
        split.by = "group")



seurat_integrated <- FindNeighbors(object = seurat_integrated, 
                                dims = 1:40)
                                
# Determine the clusters for various resolutions                                
seurat_integrated <- FindClusters(object = seurat_integrated,
                               resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))
# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.8"
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 5)
# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.6"
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 5)
# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.0.4"
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 5)
# Assign identity of clusters
Idents(object = seurat_integrated) <- "integrated_snn_res.1.4"
# Plot the UMAP
DimPlot(seurat_integrated,
        reduction = "umap",
        label = TRUE,
        label.size = 5)

################################## 按照group画图
DimPlot(seurat_integrated,
        split.by = "group",
        reduction = "umap",
        label = TRUE,
        label.size = 5)
################################## findmarker
pancreas.integrated.markers<- FindAllMarkers(seurat_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pancreas.integrated.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)


pancreas.integrated.markers %>%
    group_by(cluster) %>%
    top_n(n = 20, wt = avg_log2FC) -> top20


write.csv(top20,'top20_gene_for_marker_5_6.csv')
#write.csv(pancreas.integrated.markers,'/home/wus/for_seurat/9_6/c_t_all_markers_gene_9_12.csv')





DimPlot(seurat_integrated,
        split.by = "orig.ident")

seurat_integrated@SCTransform

FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c('Acta2','Etv5','Krt14','Krt4','Pdpn'), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE
            )


FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c('Aldh1a3', 'Cd14', 'c-Kit'), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE
            )


FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c('Fbxo22'), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE
            )


FeaturePlot(seurat_integrated, 
            reduction = "umap", 
            features = c('Epcam','Cd45'), 
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE
            )



################## 看列数


for(i in 1:length(dir)){
TPM <- read.table(dir[i], sep = " ", header = TRUE)

a=length(colnames(TPM))
print(a)
}


########### 





# Compute the relative expression of each gene per cell Use sparse matrix
# operations, if your dataset is large, doing matrix devisions the regular way
# will take a very long time.
par(mar = c(4, 8, 2, 1))
C <- Seurat_Day0@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]

boxplot(t(as.matrix(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell",
col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)

######################## 循环画质控图


#质控前小提琴图
dir.create("QC")
for(i in seq_along(dir)){
	TPM <- read.table(dir[i], sep = " ", header = TRUE)
	Seurat_Day0 <- CreateSeuratObject(counts = TPM, project = "MammaryGland.Involution1_dge")
	Seurat_Day0[["percent.mt"]] <- PercentageFeatureSet(Seurat_Day0, pattern = "^mt-")
	Seurat_Day0[["percent.Rp"]] <-PercentageFeatureSet(Seurat_Day0, pattern = '^Rp[sl]')
    plots=VlnPlot(Seurat_Day0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Rp"), ncol = 4)
	fig_name=paste("QC/vlnplot_before_qc",dir[i],".png",sep='')
	ggsave(fig_name,plot=plots,width=20,height=13)

}


