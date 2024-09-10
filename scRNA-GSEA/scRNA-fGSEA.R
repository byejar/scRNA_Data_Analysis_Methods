### scRNA_v2

rm(list=ls())


load(file='/home/wus/c_t_all_step_2_9_3.RData')
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(cowplot)
library(patchwork)
library(presto)
library(msigdbr)
library(fgsea)
library(tidyverse)



DimPlot(immune.combined.sct) 

#immune.combined.sct[['group']]<-apply(immune.combined.sct[['orig.ident']],1,substr,1,1)

DefaultAssay(immune.combined.sct) <- "SCT"
cd4_sce1<- subset(immune.combined.sct, idents = "2")

DimPlot(cd4_sce1,split.by = "orig.ident") 

#scRNA.T=subset(cd4_sce1,group == "T")
#scRNA.C=subset(cd4_sce1,group == "C")


pbmc.genes <- presto::wilcoxauc(cd4_sce1, "orig.ident")

head(pbmc.genes)

dplyr::count(pbmc.genes, orig.ident)

## 差异基因换成rank
T.genes<-pbmc.genes %>% dplyr::filter(group =="T")%>%arrange(desc(auc))%>%dplyr::select(feature,auc)
ranks<-deframe(T.genes)
head(ranks)




#m_db=msigdbr(species='Homo sapiens')
#head(m_db)
#a= m_db %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)

#m_db=msigdbr(species='Mus musculus',category='C2')

#m_df = msigdbr(species = "Mus musculus", category = "C2")
#head(m_df)
#a= m_df %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)


#将m_d的基因与通路取出并改成一个通路对应相应基因的格式
#fgsea_sets<-m_df %>% split(x =.$gene_symbol,f =.$gs_name)
#以gs_name为Factory对gene_symboli进行分类，统计落在每个gs_name中的gene_symbol的个数，并生成1ist。
#summary(fgsea_sets)


#fgseaRes<-fgsea(fgsea_sets,stats =ranks,nperm = 1000)
#npermi设置的是permutation次数


#fgseaResTidy <-fgseaRes %>% as_tibble()%>%arrange(desc(NES))
#fgseaResTidy %>% dplyr::select(-leadingEdge,-ES,-nMoreExtreme)%>%arrange(padj)%>%head()

###用下载的文件
gmt.file='/home/wus/MSigDB/msigdb_v2023.1.Mm_files_to_download_locally/msigdb_v2023.1.Mm_GMTs/m2.all.v2023.1.Mm.symbols.gmt'

pathways <- gmtPathways(gmt.file)
fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize=500)
head(fgseaRes)


topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways[topPathways], ranks, fgseaRes, 
              gseaParam=0.5)


head(fgseaRes[order(pval), ])
head(fgseaRes[order(NES), ])

#a=as.data.frame(fgseaRes)
#write.csv( a,'fgseaRes.csv')
grep('LPS',fgseaRes)
fwrite(fgseaRes, file="fgseaRes.txt", sep="\t", sep2=c("", " ", ""))



################ findmarker
#每一个细胞类型的GSEA按显著性进行降序排序
Idents(object = cd4_sce1) <- "orig.ident"
cd4_sce1[["percent.Rp"]] <-PercentageFeatureSet(cd4_sce1, pattern = '^Rp[sl]')
VlnPlot(cd4_sce1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Rp"), ncol = 4)
VlnPlot(cd4_sce1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
DefaultAssay(cd4_sce1) <- "RNA"
gesa_TvsC_allgenes<-FindMarkers(cd4_sce1, ident.1 = "T", ident.2 = "C", verbose = FALSE,test.use ="roc",logfc.threshold = 0.01,only.pos =F)
gesa_TvsC_allgenes$gene<-rownames(gesa_TvsC_allgenes)

gsea_genes<-gesa_TvsC_allgenes %>%
  arrange( desc(avg_diff),desc(myAUC)) %>%
dplyr::select(gene,myAUC)

ranks <- deframe(gsea_genes)

write.csv(ranks,'ranks.csv')

gesa_TvsC<-FindMarkers(cd4_sce1, ident.1 = "T", ident.2 = "C", logfc.threshold = 0.01,only.pos =F)
write.csv(gesa_TvsC,'gesa_TvsC.csv')


fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize=500)
head(fgseaRes)

fwrite(fgseaRes, file="fgseaRes.txt", sep="\t", sep2=c("", " ", ""))


topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways[topPathways], ranks, fgseaRes, 
              gseaParam=0.5)

plotEnrichment(pathways[["FOSTER_TOLERANT_MACROPHAGE_DN"]],
               ranks) + labs(title="FOSTER_TOLERANT_MACROPHAGE_DN")


fgseaRes[["FOSTER_TOLERANT_MACROPHAGE_DN"]]
STARK_PREFRONTAL_CORTEX_22Q11_DELETION_DN


################# 第一次删除 （mt rp）
sce<-cd4_sce1

mt.index <- grep(pattern = "^mt-", x = rownames(x = sce@assays$RNA@counts), value = FALSE)
sce <- sce[-mt.index,]


mt.genes<-grep(pattern="^mt-",x=rownames(sce),value=T)
Rp.genes<-grep(pattern="^Rp[sl]",x=rownames(sce),value=T)


Rp.index <- grep(pattern = "^Rp[sl]", x = rownames(x = sce@assays$RNA@counts), value = FALSE)
sce <- sce[-Rp.index,]

### 删除完成
################ findmarker
#每一个细胞类型的GSEA按显著性进行降序排序
cd4_sce1<-sce
Idents(object = cd4_sce1) <- "orig.ident"
cd4_sce1[["percent.mt"]] <-PercentageFeatureSet(cd4_sce1, pattern = '^mt-')
cd4_sce1[["percent.Rp"]] <-PercentageFeatureSet(cd4_sce1, pattern = '^Rp[sl]')
VlnPlot(cd4_sce1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.Rp"), ncol = 4)

DefaultAssay(cd4_sce1) <- "RNA"
gesa_TvsC_allgenes<-FindMarkers(cd4_sce1, ident.1 = "T", ident.2 = "C", verbose = FALSE,test.use ="roc",logfc.threshold = 0.01,only.pos =F)
gesa_TvsC_allgenes$gene<-rownames(gesa_TvsC_allgenes)

gsea_genes<-gesa_TvsC_allgenes %>%
  arrange( desc(avg_diff),desc(myAUC)) %>%
dplyr::select(gene,myAUC)

ranks <- deframe(gsea_genes)

write.csv(ranks,'ranks.csv')

#gesa_TvsC<-FindMarkers(cd4_sce1, ident.1 = "T", ident.2 = "C", logfc.threshold = 0.01,only.pos =F)
#write.csv(gesa_TvsC,'gesa_TvsC.csv')


fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize=500)
head(fgseaRes)

fwrite(fgseaRes, file="fgseaRes.txt", sep="\t", sep2=c("", " ", ""))


topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways[topPathways], ranks, fgseaRes, 
              gseaParam=0.5)


################# 第二次删除 （线粒体）






f=read.csv('Mouse.MitoCarta3.0-v2.csv')
f2=f[,2,drop=F]
a=apply(f2,1,grep,x=rownames(sce), value = FALSE)
MT.index <-unique(unlist(a))


sce <- sce[-MT.index,]


Atp6.index <- grep(pattern = "^Atp6", x = rownames(x = sce@assays$RNA@counts), value = FALSE)
sce <- sce[-Atp6.index,]




cd4_sce1<-sce
Idents(object = cd4_sce1) <- "orig.ident"

DefaultAssay(cd4_sce1) <- "RNA"
gesa_TvsC_allgenes<-FindMarkers(cd4_sce1, ident.1 = "T", ident.2 = "C", verbose = FALSE,test.use ="roc",logfc.threshold = 0.01,only.pos =F)
gesa_TvsC_allgenes$gene<-rownames(gesa_TvsC_allgenes)

gsea_genes<-gesa_TvsC_allgenes %>%
  arrange( desc(avg_diff),desc(myAUC)) %>%
dplyr::select(gene,myAUC)

ranks <- deframe(gsea_genes)

#write.csv(ranks,'ranks.csv')

#gesa_TvsC<-FindMarkers(cd4_sce1, ident.1 = "T", ident.2 = "C", logfc.threshold = 0.01,only.pos =F)
#write.csv(gesa_TvsC,'gesa_TvsC.csv')


fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize=500)
head(fgseaRes)

fwrite(fgseaRes, file="fgseaRes.txt", sep="\t", sep2=c("", " ", ""))


topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways[topPathways], ranks, fgseaRes, 
              gseaParam=0.5)



######################## avg_log2FC




gsea_genes<-gesa_TvsC_allgenes %>%
  arrange( desc(avg_log2FC),desc(myAUC)) %>%
dplyr::select(gene,myAUC)

ranks <- deframe(gsea_genes)

#write.csv(ranks,'ranks.csv')

#gesa_TvsC<-FindMarkers(cd4_sce1, ident.1 = "T", ident.2 = "C", logfc.threshold = 0.01,only.pos =F)
#write.csv(gesa_TvsC,'gesa_TvsC.csv')


fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize=500)
head(fgseaRes)

fwrite(fgseaRes, file="fgseaRes-5-7.txt", sep="\t", sep2=c("", " ", ""))


topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways[topPathways], ranks, fgseaRes, 
              gseaParam=0.5)


plotEnrichment(pathways[["STARK_PREFRONTAL_CORTEX_22Q11_DELETION_DN"]],
               ranks) + labs(title="STARK_PREFRONTAL_CORTEX_22Q11_DELETION_DN")


fgseaRes[["STARK_PREFRONTAL_CORTEX_22Q11_DELETION_DN"]]


#### genesets.v2023.1.Mm_macrophage.gmt
gmt.file='/home/wus/2023-4-25_012_gmk_Mammary-Gland/GSE108097/SRR6428483/genesets.v2023.1.Mm_macrophage.gmt'

pathways <- gmtPathways(gmt.file)

fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize=500)
head(fgseaRes)

fwrite(fgseaRes, file="fgseaRes-5-7-mac.txt", sep="\t", sep2=c("", " ", ""))


topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways[topPathways], ranks, fgseaRes, 
              gseaParam=0.5)


####genesets.v2023.1.Mm_macrophage_M2.gmt

gmt.file='/home/wus/2023-4-25_012_gmk_Mammary-Gland/GSE108097/SRR6428483/genesets.v2023.1.Mm_macrophage_M2.gmt'

pathways <- gmtPathways(gmt.file)

fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize=500)
head(fgseaRes)

fwrite(fgseaRes, file="fgseaRes-5-7-mac-m2.txt", sep="\t", sep2=c("", " ", ""))


topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways[topPathways], ranks, fgseaRes, 
              gseaParam=0.5)