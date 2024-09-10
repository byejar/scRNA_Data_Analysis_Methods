# library(Seurat)
# dir="/home/wus/027-2023-10-25_tl_HCC-DDIT4/DATA/"
# list.files(dir)
# #[1] "barcodes.tsv.gz" "features.tsv.gz" "matrix.mtx.gz" 

# counts <- Read10X(data.dir = dir)
# class(counts)
# #[1] "dgCMatrix"
# #attr(,"package")
# #[1] "Matrix"

# scRNA <- CreateSeuratObject(counts = counts)
# scRNA
# #An object of class Seurat 
# #33694 features across 9112 samples within 1 assay 
# #Active assay: RNA (33694 features, 0 variable features)

# #info<-read.table('/home/wus/027-2023-10-25_tl_HCC-DDIT4/DATA/GSE151530_Info.txt',header=1,sep='\t',col.names='Sample')
# info<-read.table('/home/wus/027-2023-10-25_tl_HCC-DDIT4/DATA/GSE151530_Info.txt',header=1,sep='\t',row.names=3)

# scRNA <- AddMetaData(object = scRNA,                #seurat对象
#                     metadata = info              #需要添加的metadata
#                     ) 

# Malignant_cells<-subset(scRNA,subset=Type=='Malignant cells')

# ####
# pbmc <- NormalizeData(Malignant_cells, normalization.method = "LogNormalize", scale.factor = 10000)
# #鉴定细胞间表达量高变的基因（feature selection）
# #这一步的目的是鉴定出细胞与细胞之间表达量相差很大的基因，用于后续鉴定细胞类型，
# #我们使用默认参数，即“vst”方法选取2000个高变基因。
# pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# # Identify the 10 most highly variable genes
# all.genes <- rownames(pbmc)
# pbmc <- ScaleData(pbmc, features = all.genes)
# #Perform linear dimensional reduction
# pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# #Examine and visualize PCA results a few different ways
# pbmc <- JackStraw(pbmc, num.replicate = 100)
# pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
# JackStrawPlot(pbmc, dims = 1:15)
# ElbowPlot(pbmc)
# #UMAP
# pbmc <- RunUMAP(pbmc, dims = 1:10, label = T)
# head(pbmc@reductions$umap@cell.embeddings) # 提取UMAP坐标值。
# #note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
# p1 <- DimPlot(pbmc, reduction = "umap")
# #T-SNE
# pbmc <- RunTSNE(pbmc, dims = 1:10)
# head(pbmc@reductions$tsne@cell.embeddings)
# p2 <- DimPlot(pbmc, reduction = "tsne",group.by='Sample')
# p1 + p2
# Malignant_cells<-pbmc
# saveRDS(Malignant_cells, file = "/home/wus/027-2023-10-25_tl_HCC-DDIT4/Malignant_cells.rds")  #保存数据，用于后续个性化分析

# FeaturePlot(Malignant_cells, features =  'DDIT4',
#     cols = c("grey", 'orange',"red"), reduction = "tsne")


############ conda activate scRNA
library(Seurat)
Malignant_cells<-readRDS(file = "/home/wus/027-2023-10-25_tl_HCC-DDIT4/Malignant_cells.rds")
# library("Nebulosa")
# plot_density(Malignant_cells, 'DDIT4', reduction = "tsne")
DimPlot(Malignant_cells, reduction = "umap",group.by='Sample', label = T)

# Malignant_cells$DDIT4_express = ifelse(Malignant_cells@assays$RNA@counts['DDIT4',]>mean(Malignant_cells@assays$RNA@counts['DDIT4',]),'DDIT4_H','DDIT4_L')
# mean(Malignant_cells@assays$RNA@counts['DDIT4',])
# table(Malignant_cells$DDIT4_express)


Malignant_cells$DDIT4_express = ifelse(Malignant_cells@assays$RNA@counts['DDIT4',]>0,'DDIT4_H','DDIT4_L')

# write.table(Malignant_cells[["RNA"]]@data,'/home/wus/027-2023-10-25_tl_HCC-DDIT4/Malignant_cells.txv',sep='\t')

# DDIT4<-as.data.frame(Malignant_cells$DDIT4_express)
# write.table(x = DDIT4,file = "cell_metadata.csv",sep = ",",quote = F)


#### 保存文件

library(Matrix)
######### DDIT4_H
setwd('/home/wus/027-2023-10-25_tl_HCC-DDIT4/Malignant_cells/DDIT4_H')
Malignant_cells_D4H<-subset(Malignant_cells,subset=DDIT4_express=='DDIT4_H')
##构建稀疏矩阵
sparse <- Matrix(Malignant_cells_D4H@assays$RNA@counts,sparse = T)
##输出基因feature
features <- data.frame(ID = sparse@Dimnames[[1]],Name = sparse@Dimnames[[1]],EX = "Gene Expression")
write.table(x = features,file = "features.tsv",sep = "\t",quote = F,col.names = F,row.names = F)
##输出barcode
write(x = sparse@Dimnames[[2]],file = "barcodes.tsv")
##输出表达矩阵
writeMM(sparse,file = "matrix.mtx")

write.table(file="expression.tsv",as.data.frame(Malignant_cells_D4H[["RNA"]]@counts),quote = F,sep='\t')

######### DDIT4_L
setwd('/home/wus/027-2023-10-25_tl_HCC-DDIT4/Malignant_cells/DDIT4_L')
Malignant_cells_D4L<-subset(Malignant_cells,subset=DDIT4_express=='DDIT4_L')
##构建稀疏矩阵
sparse <- Matrix(Malignant_cells_D4L@assays$RNA@counts,sparse = T)
##输出基因feature
features <- data.frame(ID = sparse@Dimnames[[1]],Name = sparse@Dimnames[[1]],EX = "Gene Expression")
write.table(x = features,file = "features.tsv",sep = "\t",quote = F,col.names = F,row.names = F)
##输出barcode
write(x = sparse@Dimnames[[2]],file = "barcodes.tsv")
##输出表达矩阵
writeMM(sparse,file = "matrix.mtx")

write.table(file="expression.tsv",as.data.frame(Malignant_cells_D4L[["RNA"]]@counts),quote = F,sep='\t')

#### conda activate python3
cd /home/wus/027-2023-10-25_tl_HCC-DDIT4/Malignant_cells/DDIT4_H
# compass --data-mtx matrix.mtx features.tsv barcodes.tsv --num-processes 10 --species homo_sapiens --output-dir /home/wus/027-2023-10-25_tl_HCC-DDIT4/compass/

# compass --data-mtx matrix.mtx features.tsv barcodes.tsv --microcluster-size 10 --num-processes 10 --species homo_sapiens --output-dir /home/wus/027-2023-10-25_tl_HCC-DDIT4/compass/

# compass --data expression.tsv --microcluster-size 10 --num-processes 10 --species homo_sapiens --output-dir /home/wus/027-2023-10-25_tl_HCC-DDIT4/compass/

compass --data expression.tsv --microcluster-size 10 --num-processes 10 --species homo_sapiens --output-dir /scratch/wus/compass/DDIT4_H


cd /home/wus/027-2023-10-25_tl_HCC-DDIT4/Malignant_cells/DDIT4_L
compass --data expression.tsv --microcluster-size 10 --num-processes 10 --species homo_sapiens --output-dir /scratch/wus/compass/DDIT4_L

############### conda activate irGSEA
##################### irgsea使用自己的geneset注意：1.必须是list形式 2.命名必须是h.gsets.list
library(Seurat)
library(irGSEA)
library(UCell)
Malignant_cells<-readRDS(file = "/home/wus/027-2023-10-25_tl_HCC-DDIT4/Malignant_cells.rds")
Malignant_cells$DDIT4_express = ifelse(Malignant_cells@assays$RNA@counts['DDIT4',]>0,'DDIT4_H','DDIT4_L')

gmt.file='/home/wus/027-2023-10-25_tl_HCC-DDIT4/genesets.Glucose-Deprived.PMID17409444.gmt'

pathways <- GSEABase::getGmt(gmt.file)

gmt.geneset<- list(c('EIF2AK3','PPEF2','CPR8','DNAJC3','PRKCZ','CDC2','SGKL','ZAP70','TOPK','CDT1','CNK','TRPS1','ARNTL','RELB','PREB','ELF3','ART3','FUT1','NEU1','GOT1','ASNS','GFPT1','PCK2','CBR3','PYCR1','SARS','NDUFB3','RRM2','B3GNT6','PIGA','TXNRD1','PPIB','PRAME','FKBP2','DFNA5','RTP801','SEMG1','H4FE','TFRC','GRO2','GADD45','TGIF','ACATN','CLGN','DNAJB11','CALR','CANX','SLC3A2','UGTREL1','KDELR2','KCNF1','PDIA4','TXNIP','SPS2','PDIA6','GRP58','GFAP','GABARAPL3','OCLN','KIF5C','ANLN'))

names(gmt.geneset) <- c("glucose-deprived")
h.gsets.list<-gmt.geneset

Malignant_cells <- irGSEA.score(object = Malignant_cells, assay = "RNA", slot = "data",geneset =h.gsets.list,msigdb = F,
method = c( "UCell"), kcdf = 'Gaussian')

scatterplot <- irGSEA.density.scatterplot(object = Malignant_cells,
                             method = "UCell",
                             show.geneset = 'glucose-deprived',
                             reduction = "tsne")
scatterplot

library("Nebulosa")
plot_density(Malignant_cells, "DDIT4",reduction = "umap")

plot_density(Malignant_cells, "APOH",reduction = "umap")


scatterplot <- irGSEA.density.scatterplot(object = Malignant_cells,
                             method = "UCell",
                             show.geneset = 'glucose-deprived',
                             reduction = "umap")
scatterplot







########### test

gmt.geneset<- list(c('SOX4','PKM','HILPDA','CD24','MIF','ALDOA','ENO1','INHBE','UPP1','SLC2A1','PPP1R3C','OAZ1','DDIT3','PGAM1','S100A11','SH3BGRL3'))

names(gmt.geneset) <- c("glucose-test")
h.gsets.list<-gmt.geneset

Malignant_cells <- irGSEA.score(object = Malignant_cells, assay = "RNA", slot = "data",geneset =h.gsets.list,msigdb = F,
method = c( "UCell"), kcdf = 'Gaussian')

scatterplot <- irGSEA.density.scatterplot(object = Malignant_cells,
                             method = "UCell",
                             show.geneset = 'glucose-test',
                             reduction = "umap")
scatterplot


####
gmt.geneset<- list(c('TRIM11','PHGDH','RHOF','NFYA','SESN2'))

names(gmt.geneset) <- c("glucose-test")
h.gsets.list<-gmt.geneset

Malignant_cells <- irGSEA.score(object = Malignant_cells, assay = "RNA", slot = "data",geneset =h.gsets.list,msigdb = F,
method = c( "UCell"), kcdf = 'Gaussian')

scatterplot <- irGSEA.density.scatterplot(object = Malignant_cells,
                             method = "UCell",
                             show.geneset = 'glucose-test',
                             reduction = "umap")
scatterplot

################################
gmt.geneset<- list(c('PKM','ALDOA','SLC2A1','DDIT3','PGAM1','STC2','LDHA','NDRG1','TPI1','CTSA'))

names(gmt.geneset) <- c("glucose-test")
h.gsets.list<-gmt.geneset

Malignant_cells <- irGSEA.score(object = Malignant_cells, assay = "RNA", slot = "data",geneset =h.gsets.list,msigdb = F,
method = c( "UCell"), kcdf = 'Gaussian')

scatterplot <- irGSEA.density.scatterplot(object = Malignant_cells,
                             method = "UCell",
                             show.geneset = 'glucose-test',
                             reduction = "umap")
scatterplot


################################ 11-27
# gmt.geneset<- list(c('FUT1+','CHAC1+','DDIT4+','FMOD+','F2RL2+','ASNS+','ATF3+',
#     'SH3BP2+','ADM2+','INHBE+','STC2+','ASNSP1+','ULBP1+','SLC7A11+','DDIT3+','SESN2+',
#     'TUBE1+','GDF15+','SLC22A15+','MTHFD2+','PMAIP1+','TSC22D3+','JDP2+',
#     'TXNIP-','TNFSF10-','PRAME-','CCNE1-','TMEM200B-','SIX4-','KANK2-',
#     'RASL11A-','SAA2-','IFFO2-','GPER1-','FAM111B-','LRRC45-','LPCAT1-','TMEM51-'))
gmt.geneset<- list(c(
    'TXNIP-','TNFSF10-','PRAME-','CCNE1-','TMEM200B-','SIX4-','KANK2-',
    'RASL11A-','SAA2-','IFFO2-','GPER1-','FAM111B-','LRRC45-','LPCAT1-','TMEM51-'))

names(gmt.geneset) <- c("glucose-test")
h.gsets.list<-gmt.geneset

Malignant_cells <- irGSEA.score(object = Malignant_cells, assay = "RNA", slot = "data",geneset =h.gsets.list,msigdb = F,
method = c( "UCell"), kcdf = 'Gaussian')

scatterplot <- irGSEA.density.scatterplot(object = Malignant_cells,
                             method = "UCell",
                             show.geneset = 'glucose-test',
                             reduction = "umap")
scatterplot




library("Nebulosa")
plot_density(Malignant_cells, c('FUT1','CHAC1','DDIT4','FMOD','ASNS','ATF3','SH3BP2','ADM2','INHBE','STC2','ULBP1','SLC7A11','DDIT3','SESN2','TUBE1','GDF15','SLC22A15','MTHFD2','PMAIP1','TSC22D3','JDP2'),reduction = "umap")


# library(Seurat)
# counts <- Read10X(data.dir = '/home/wus/027-2023-10-25_tl_HCC-DDIT4/Malignant_cells/')
# class(counts)

# scRNA <- CreateSeuratObject(counts = counts)
# scRNA

########################################################
############################### conda activate python3
########################################################
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from compass_analysis import cohens_d, wilcoxon_test, get_reaction_consistencies, get_metareactions, labeled_reactions, amino_acid_metab


#####
cell_md = pd.read_csv("/home/wus/027-2023-10-25_tl_HCC-DDIT4/Malignant_cells/cell_metadata.csv", index_col=0)
reaction_penalties = pd.read_csv("/scratch/wus/compass/reactions.tsv", sep="\t", index_col = 0)
micropools = pd.read_csv("/scratch/wus/compass/micropools.tsv", sep="\t", index_col=0)


clusters = {}
for cell in micropools.index:
    mc = micropools.loc[cell, 'microcluster']
    if mc in clusters:
        clusters[mc] += [cell]
    else:
        clusters[mc] = [cell]

Malignant_cells$DDIT4_express'] == 'DDIT4_H

DDIT4_H, DDIT4_L = {cl:0 for cl in clusters}, {cl:0 for cl in clusters}
for cl in clusters:
    for cell in clusters[cl]:
        cell_type = cell_md.loc[cell, 'Malignant_cells$DDIT4_express']
        if cell_type == 'DDIT4_H':
            DDIT4_H[cl] += 1
        elif cell_type == 'DDIT4_L':
            DDIT4_L[cl] += 1
        else:
            print("Should not happen")

pctTh17p = {cl:DDIT4_H[cl] / (DDIT4_H[cl] + DDIT4_L[cl]) for cl in clusters}

def mc_type(pct):
    if pct > 0.9:
        return 'DDIT4_H'
    elif pct < 0.1:
        return 'DDIT4_L'
    else:
        return 'Uncertain'
micropool_md = {'cluster_'+str(cl):mc_type(pctTh17p[cl]) for cl in pctTh17p}
micropool_md = pd.DataFrame.from_dict(micropool_md, orient='index', columns=['cell_type'])
micropool_md.to_csv("/scratch/wus/compass/cluster_metadata.csv")

##############
cell_metadata = pd.read_csv("/scratch/wus/compass/cluster_metadata.csv", index_col = 0)
reaction_penalties = pd.read_csv("/scratch/wus/compass/reactions.tsv", sep="\t", index_col = 0)

Th17p_cells = cell_metadata.index[cell_metadata['cell_type'] == 'DDIT4_H']
Th17n_cells = cell_metadata.index[cell_metadata['cell_type'] == 'DDIT4_L']

def get_reaction_consistencies(compass_reaction_penalties, min_range=1e-3):
	df = -np.log(compass_reaction_penalties + 1)
	df = df[df.max(axis=1) - df.min(axis=1) >= min_range]
	df = df - df.min().min()
	return df


reaction_consistencies = get_reaction_consistencies(reaction_penalties)

wilcox_results = wilcoxon_test(reaction_consistencies, Th17p_cells, Th17n_cells)

reaction_metadata = pd.read_csv("/scratch/wus/compass/Compass/notebooks/extdata/RECON2//reaction_metadata.csv", index_col = 0)

wilcox_results['metadata_r_id'] = ""
for r in wilcox_results.index:
    if r in reaction_metadata.index:
        wilcox_results.loc[r, 'metadata_r_id'] = r
    elif r[:-4] in reaction_metadata.index:
        wilcox_results.loc[r, 'metadata_r_id'] = r[:-4]
    else:
        print("Should not occur")

W = wilcox_results.merge(reaction_metadata, how='left',
                         left_on='metadata_r_id', right_index=True, validate='m:1')
W = W[W['confidence'].isin([0,4])]
W = W[~W['EC_number'].isna()]
W.loc[(W['formula'].map(lambda x: '[m]' not in x)) & (W['subsystem'] == "Citric acid cycle"), 'subsystem'] = 'Other'


def plot_differential_scores(data, title, c):
    plt.figure(figsize=(10,10))
    axs = plt.gca()
    axs.scatter(data['cohens_d'], -np.log10(data['adjusted_pval']), c=c)
    axs.set_xlabel("Cohen's d", fontsize=16)
    axs.set_ylabel("-log10 (Wilcoxon-adjusted p)", fontsize=16)
    #Everything after this should be tweaked depending on your application
    axs.set_xlim(-2.2, 2.2)
    axs.axvline(0, dashes=(3,3), c='black')
    axs.axhline(1, dashes=(3,3), c='black')
    axs.set_title(title, fontdict={'fontsize':20})
    axs.annotate('', xy=(0.5, -0.08), xycoords='axes fraction', xytext=(0, -0.08),
            arrowprops=dict(arrowstyle="<-", color='#348C73', linewidth=4))
    axs.annotate('DDIT4_H', xy=(0.75, -0.12), xycoords='axes fraction', fontsize=16)
    axs.annotate('', xy=(0.5, -0.08), xycoords='axes fraction', xytext=(1, -0.08),
            arrowprops=dict(arrowstyle="<-", color='#E92E87', linewidth=4))
    axs.annotate('DDIT4_L', xy=(0.25, -0.12), xycoords='axes fraction', fontsize=16)
    for r in data.index:
        if r in labeled_reactions:
            x = data.loc[r, 'cohens_d']
            y = -np.log10(data.loc[r, 'adjusted_pval'])
            offset = (20, 0)
            if x < 0:
                offset = (-100, -40)
            axs.annotate(labeled_reactions[r], (x,y), xytext = offset,
                         textcoords='offset pixels', arrowprops={'arrowstyle':"-"})


filtered_data = pd.concat([W[W['subsystem'] == "Glycolysis/gluconeogenesis"],
             W[W['subsystem'] == "Citric acid cycle"],
            W[W['subsystem'].isin(amino_acid_metab)],
           W[W['subsystem'] == "Fatty acid oxidation"]])


data = W[W['subsystem'] == "Glycolysis/gluconeogenesis"]
plot_differential_scores(data, title='Glycolysis', c="#695D73")
plt.show()

data = W[W['subsystem'].isin(amino_acid_metab)].copy()
data['adjusted_pval'] = data['adjusted_pval'].clip(1e-12)
plot_differential_scores(data, "Amino Acid Metabolism", c="#BF1E2E")
plt.show()


############

data = W[~W['subsystem'].isin(["Miscellaneous", "Unassigned"])]
data = data[~data['subsystem'].map(lambda x: "Transport" in x or "Exchange" in x or x == "Other")]
items, counts = np.unique(data['subsystem'], return_counts=True)
items = [items[i] for i in range(len(items)) if counts[i] > 5] #filter(n() > 5) %>%
data = data[data['subsystem'].isin(items)]


plt.figure(figsize=(12,12))
axs = plt.gca()
#Sorts the reactions for plotting
d = data[data['adjusted_pval'] < 0.1].groupby('subsystem')['cohens_d'].median().abs()
axs.scatter(d[d.argsort], d[d.argsort].index, alpha=0)
color = data['cohens_d'].map(lambda x: 'r' if x >= 0 else 'b')
alpha = data['adjusted_pval'].map(lambda x: 1.0 if x < 0.1 else 0.25)
axs.scatter(data['cohens_d'], data['subsystem'], c=color, alpha=alpha)
axs.set_xlabel("Cohen's d")