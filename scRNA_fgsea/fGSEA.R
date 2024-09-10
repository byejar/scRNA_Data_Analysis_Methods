## scRNA_v2
### fGSEA
#/home/wus/MSigDB/msigdb_v2023.1.Mm_files_to_download_locally/msigdb_v2023.1.Mm_GMTs/
library(devtools)
library(presto)
library(msigdbr)
library(ggplot2)
library(tidyverse)
library(Seurat)
library(fgsea)

load(file = '/home/wus/for_seurat/8_1/immune_rename_2023-7-4.Rdata')

head(mac_sce_rename)

Idents(mac_sce_rename) <- "cell_type"
mac.cells <- subset(mac_sce_rename , idents = c("M1","M2","fib","neutrophil", "T cell","DC","B cell"))
table(mac.cells$cell_type)

pbmc.genes <- wilcoxauc(mac.cells, 'orig.ident')

head(pbmc.genes)


# msigdbr_species()


# m_df<- msigdbr(species = "Mus musculus")

# head(m_df)

# table(m_df$gs_cat)

# fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

pbmc.genes %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)      #进行降序排序
  
  
#   cluster0.genes<- pbmc.genes %>%
#   dplyr::filter(group == "T") %>%
#   arrange(desc(auc)) %>%
#   dplyr::select(feature, auc)

# ranks<- deframe(cluster0.genes)

# head(ranks)


# fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)



# fgseaResTidy <- fgseaRes %>%
#   as_tibble() %>%
#   arrange(desc(NES))

# fgseaResTidy %>%
#   dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
#   arrange(padj) %>%
#   head()
  
  
#   # 显示top20信号通路
# ggplot(fgseaResTidy %>% filter(padj < 0.05) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
#   geom_col(aes(fill= NES < 7.5)) +
#   coord_flip() +
#   labs(x="Pathway", y="Normalized Enrichment Score",
#        title="m_all pathways NES from GSEA") +
#   theme_minimal() ####以7.5进行绘图填色
  
# #单个画图
#   plotEnrichment(fgsea_sets[["GSE10325_CD4_TCELL_VS_MYELOID_UP"]],
#                ranks) + labs(title="GSE10325 CD4 TCELL VS MYELOID UP")
               
               


#######################
gmt<-gmtPathways("/home/wus/MSigDB/msigdb_v2023.1.Mm_files_to_download_locally/msigdb_v2023.1.Mm_GMTs/m_all-uniq.gmt")
# fgseaRes<- fgsea(gmt, stats = ranks, nperm = 1000)

# fgseaResTidy <- fgseaRes %>%
#   as_tibble() %>%
#   arrange(desc(NES))
  
#   fgseaResTidy %>%
#   dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
#   arrange(padj) %>%
#   head(n=30)

# fgseadataframe<-as.data.frame(fgseaResTidy)
#write.table(fgseadataframe,'FGSEA.txt')

# fwrite(fgseadataframe,'FGSEA.txt')


#  plotEnrichment(gmt[["WP_ELECTRON_TRANSPORT_CHAIN"]],ranks) + labs(title="WORSCHECH_TUMORREJECTION_UP")
  
  ##### 不选
    cluster0.genes<- pbmc.genes %>%
  #dplyr::filter(group == "T") %>%
  arrange(desc(auc)) %>%
  dplyr::select(feature, auc)

ranks<- deframe(cluster0.genes)

head(ranks)


fgseaRes<- fgsea(gmt, stats = ranks, nperm = 1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
  
  # fgseaResTidy %>%
  # dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
  # arrange(padj) %>%
  # head(n=30)

fgseadataframe<-as.data.frame(fgseaResTidy)
#write.table(fgseadataframe,'FGSEA.txt')

fwrite(fgseadataframe,'FGSEA-v2-v2.txt')


# ggplot(fgseaResTidy %>% filter(padj < 0.05) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
#   geom_col(aes(fill= NES < 7.5)) +
#   coord_flip() +
#   labs(x="Pathway", y="Normalized Enrichment Score",
#        title="m_all pathways NES from GSEA") +
#   theme_minimal() ####以7.5进行绘图填色
  

ggplot(fgseaRes %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pval < 0.05) %>% 
         head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES)) +
  coord_flip() +
  labs(x="KEGG", y="Normalized Enrichment Score",title="KEGG gene sets NES from GSEA")  

plotEnrichment(fgseaRes[["ABBUD_LIF_SIGNALING_1_DN"]],ranks) + labs(title="ABBUD_LIF_SIGNALING_1_DN")  

plotEnrichment(gmt[["ABBUD_LIF_SIGNALING_1_DN"]],
               ranks,ticksSize = 20) + 
  labs(title="Programmed Cell Death")

topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(examplePathways[topPathways], exampleRanks, fgseaRes, 
              gseaParam = 0.5)

############### irGSEA
pbmc3k.final <- irGSEA.score(object = pbmc3k.final, assay = "RNA", 
                             slot = "data", seeds = 123, ncores = 1,
                             min.cells = 3, min.feature = 0,
                             custom = F, geneset = NULL, msigdb = T, 
                             species = "Homo sapiens", category = "H",  
                             subcategory = NULL, geneid = "symbol",
                             method = c("AUCell", "UCell", "singscore", 
                                        "ssgsea"),
                             aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                             kcdf = 'Gaussian')

############## 
## scRNA_v2 
### fGSEA findmarker
#/home/wus/MSigDB/msigdb_v2023.1.Mm_files_to_download_locally/msigdb_v2023.1.Mm_GMTs/
# library(devtools)
# library(presto)
# library(msigdbr)
# library(ggplot2)
# library(tidyverse)
# library(Seurat)
# library(fgsea)

# load(file = '/home/wus/for_seurat/8_1/immune_rename_2023-7-4.Rdata')

# head(mac_sce_rename)

# Idents(mac_sce_rename) <- "cell_type"
# mac.cells <- subset(mac_sce_rename , idents = c("M1","M2","fib","neutrophil", "T cell","DC","B cell"))

# head(mac.cells)
# Idents(mac.cells) <- "orig.ident"


# pbmc_merged <- PrepSCTFindMarkers(object = mac.cells)


# cmarkers <- FindMarkers(pbmc_merged ,ident.1 = 'T',ident.2 = "C", only.pos = FALSE, 
#                         min.pct = 0.1, logfc.threshold = 0)

# #基因按logFC排序
# cmarkers$genes = rownames(cmarkers)
# cluster0.genes<- cmarkers %>% arrange(desc(avg_log2FC)) %>% dplyr::select(genes,avg_log2FC)
# ranks<- deframe(cluster0.genes)


# ###
# gmt<-gmtPathways("/home/wus/MSigDB/msigdb_v2023.1.Mm_files_to_download_locally/msigdb_v2023.1.Mm_GMTs/m_all-uniq.gmt")
# #fgsea
# fgseaRes<- fgsea(gmt, stats = ranks, nperm = 1000)


# fgseaResTidy <- fgseaRes %>%
#   as_tibble() %>%
#   arrange(desc(NES))
  
#   fgseaResTidy %>%
#   dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>%
#   arrange(padj) %>%
#   head(n=30)

# fgseadataframe<-as.data.frame(fgseaResTidy)
# #write.table(fgseadataframe,'FGSEA.txt')

# fwrite(fgseadataframe,'FGSEA-v3.txt')

# fgseadataframe[1:30,1:4]


# ggplot(fgseaResTidy %>% filter(padj < 0.05) %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
#   geom_col(aes(fill= NES < 7.5)) +
#   coord_flip() +
#   labs(x="Pathway", y="Normalized Enrichment Score",
#        title="m_all pathways NES from GSEA") +
#   theme_minimal() ####以7.5进行绘图填色


# ggplot(fgseaRes %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pval < 0.05) %>% 
#          head(n= 20), aes(reorder(pathway, NES), NES)) +
#   geom_col(aes(fill= NES)) +
#   coord_flip() +
#   labs(x="KEGG", y="Normalized Enrichment Score",title="KEGG gene sets NES from GSEA")


#### irGSEA
library(GSVA)
library(devtools)
#library(presto)
#library(ggplot2)
#library(tidyverse)
library(Seurat)

load(file = '/home/wus/for_seurat/8_1/immune_rename_2023-7-4.Rdata')

head(mac_sce_rename)

Idents(mac_sce_rename) <- "cell_type"
mac.cells <- subset(mac_sce_rename , idents = c("M1","M2","fib","neutrophil", "T cell","DC","B cell"))
table(mac.cells$cell_type)

#pbmc.genes <- wilcoxauc(mac.cells, 'orig.ident')

head(mac.cells)
Idents(mac.cells) <- "orig.ident"
gene.expr <-  as.matrix(mac.cells[["RNA"]]@data)

dim(gene.expr)

gmt <- GSEABase::getGmt("/home/wus/MSigDB/msigdb_v2023.1.Mm_files_to_download_locally/msigdb.v2023.1.Mm.symbols.txt")
# gmt <- GSEABase::getGmt("/home/wus/MSigDB/msigdb_v2023.1.Mm_files_to_download_locally/msigdb_v2023.1.Mm_GMTs/m1.all.v2023.1.Mm.entrez.gmt")

# library(msigdbr)
# m_df<- msigdbr(species = "Mus musculus")
# fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

gsva.result <- GSVA::gsva(gene.expr, fgsea_sets,kcdf='Gaussian')

mymatrix <- t(gsva.result)
mymatrix <- cbind(mymatrix,as.data.frame(mac.cells$orig.ident))
colnames(mymatrix)[ncol(mymatrix)] <- 'orig.ident'
head(mymatrix)


####
library(UCell)
library(irGSEA)
Idents(mac.cells) <- "cell_type"

pbmc3k.final <- irGSEA.score(object = mac.cells, assay = "RNA", 
                             slot = "data", seeds = 123, ncores = 1,
                             min.cells = 3, min.feature = 0,
                             custom = F, geneset = NULL, msigdb = T, 
                             species = "Mus musculus", #category = "H",  
                             subcategory = NULL, geneid = "symbol",
                             method = c("AUCell", "UCell", "singscore", 
                                        "ssgsea"),
                             aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                             kcdf = 'Gaussian')

result.dge <- irGSEA.integrate(object = pbmc3k.final, 
                               group.by = "orig.ident",
                               metadata = NULL, col.name = NULL,
                               method = c("AUCell","UCell","singscore",
                                          "ssgsea"))


#ggsave(irGSEA.heatmap.plot,'irGSEA.heatmap.plot.png')

setwd('/home/wus/12-8_001_PTEN-scRNA/')

pdf("irGSEA.heatmap.plot.pdf")
irGSEA.heatmap.plot <- irGSEA.heatmap(object = result.dge, 
                                      method = "RRA",
                                      top = 50, 
                                      show.geneset = NULL)
irGSEA.heatmap.plot
#print(irGSEA.heatmap.plot)
dev.off()
​

######################################################################### fgsea 8-2


load(file = '/home/wus/for_seurat/8_1/immune_rename_2023-7-4.Rdata')

head(mac_sce_rename)

Idents(mac_sce_rename) <- "cell_type"
mac.cells <- subset(mac_sce_rename , idents = c("M1","M2","neutrophil", "T cell","DC","B cell"))
table(mac.cells$cell_type)

pbmc.genes <- wilcoxauc(mac.cells, 'orig.ident')

head(pbmc.genes)



pbmc.genes %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)      #进行降序排序
  


#######################
gmt<-gmtPathways("/home/wus/MSigDB/msigdb_v2023.1.Mm_files_to_download_locally/msigdb.v2023.1.Mm.symbols.txt")

  
  ##### 不选
    cluster0.genes<- pbmc.genes %>%
  #dplyr::filter(group == "T") %>%
  arrange(desc(auc)) %>%
  dplyr::select(feature, auc)

ranks<- deframe(cluster0.genes)

head(ranks)


#fgseaRes<- fgsea(gmt, stats = ranks, nperm = 1000)


fgseaRes<-fgseaMultilevel(gmt, stats = ranks, nPermSimple = 1000)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
  


fgseadataframe<-as.data.frame(fgseaResTidy)
#write.table(fgseadataframe,'FGSEA.txt')

fwrite(fgseadataframe,'FGSEA-v2-v2.txt')




ggplot(fgseaRes %>% as_tibble() %>% arrange(desc(NES)) %>% filter(pval < 0.05) %>% 
         head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES)) +
  coord_flip() +
  labs(x="KEGG", y="Normalized Enrichment Score",title="KEGG gene sets NES from GSEA")  

plotEnrichment(gmt[["ABBUD_LIF_SIGNALING_1_DN"]],ranks) + labs(title="ABBUD_LIF_SIGNALING_1_DN")  

plotEnrichment(gmt[["WORSCHECH_TUMOR_REJECTION_UP"]],ranks,ticksSize = 20) + labs(title="WORSCHECH_TUMOR_REJECTION_UP")
plotEnrichment(gmt[["GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MEDIATED"]],ranks,ticksSize = 20) + labs(title="GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MEDIATED")
plotEnrichment(gmt[["GOBP_POSITIVE_REGULATIONOF_CELL_KILLING"]],ranks,ticksSize = 20) + labs(title="GOBP_POSITIVE_REGULATIONOF_CELL_KILLING")
plotEnrichment(gmt[["WORSCHECH_TUMOR_REJECTION_UP"]],ranks,ticksSize = 20) + labs(title="WORSCHECH_TUMOR_REJECTION_UP")

topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(examplePathways[topPathways], exampleRanks, fgseaRes, 
              gseaParam = 0.5)

#################################################################

#################################################################
Idents(object = mac.cells) <- "orig.ident"

DefaultAssay(mac.cells) <- "RNA"
gesa_TvsC_allgenes<-FindMarkers(mac.cells, ident.1 = "T", ident.2 = "C", verbose = FALSE,test.use ="roc",logfc.threshold = 0.01,only.pos =F)
gesa_TvsC_allgenes$gene<-rownames(gesa_TvsC_allgenes)

gsea_genes<-gesa_TvsC_allgenes %>%
  arrange( desc(avg_diff),desc(myAUC)) %>%
dplyr::select(gene,myAUC)

ranks <- deframe(gsea_genes)

#write.csv(ranks,'ranks.csv')

#gesa_TvsC<-FindMarkers(cd4_sce1, ident.1 = "T", ident.2 = "C", logfc.threshold = 0.01,only.pos =F)
#write.csv(gesa_TvsC,'gesa_TvsC.csv')


fgseaRes <- fgsea(gmt, ranks, minSize=15, maxSize=500)
head(fgseaRes)

fwrite(fgseaRes, file="fgseaRes.txt", sep="\t", sep2=c("", " ", ""))


topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(gmt[topPathways], ranks, fgseaRes, 
              gseaParam=0.5)

plotEnrichment(gmt[["WORSCHECH_TUMOR_REJECTION_UP"]],ranks,ticksSize = 20) + labs(title="WORSCHECH_TUMOR_REJECTION_UP")
plotEnrichment(gmt[["GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MEDIATED"]],ranks,ticksSize = 20) + labs(title="GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MEDIATED")
plotEnrichment(gmt[["GOBP_POSITIVE_REGULATIONOF_CELL_KILLING"]],ranks,ticksSize = 20) + labs(title="GOBP_POSITIVE_REGULATIONOF_CELL_KILLING")

fgseaRes[ES > 0][head(order(pval), n=100), pathway]
fgseaRes[ES < 0][head(order(pval), n=100), pathway]

write.table(gesa_TvsC_allgenes,'gesa_TvsC_allgenes.txt',sep='\t')

gsea_genes<-gesa_TvsC_allgenes %>%
  arrange( desc(avg_diff),desc(myAUC)) 
write.table(gsea_genes,'gsea_genes.txt',sep='\t')
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