#https://blog.csdn.net/Jagger_Lin/article/details/104350750

#-----------------------------------------------------------------------------
#                  安装TCGAbiolinks
#-----------------------------------------------------------------------------
#安装稳定版biolinks
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.10")
#BiocManager::install("TCGAbiolinks")

#安装开发版biolinks
#devtools::install_github('BioinformaticsFMRP/TCGAbiolinks')

#-----------------------------------------------------------------------------
#                   手动安装示例区域（自己看着办）
#-----------------------------------------------------------------------------
#常规联网查找安装，
#install.packages("example_name")
#通过bioconductor资源安装
#BiocManager::install("example_name")
#to search through available packages programmatically
#BiocManager::available("example_name")
#to flag packages that are either out-of-date or too new for your version of Bioconductor. 
#BiocManager::valid()
#-----------------------------------------------------------------------------
#                   加载库    环境名：TCGA
#-----------------------------------------------------------------------------
#biolinks使用
library(SummarizedExperiment)
library(TCGAbiolinks)

#分析绘图使用
#BiocManager::install("edgeR")
#install.packages("gplots")
#BiocManager::install("ggplot2")
#library(edgeR)
#library(gplots)
#library(ggplot2)
#require(gplot2)
#-----------------------------------------------------------------------------
#                   创建工作目录
#-----------------------------------------------------------------------------
#创建新目录
#dir.create("c:/TCGA");
#dir.create("c:/TCGA/CESC");
#目录
#work_dir <- "c:/TCGA/CESC"
#将当前目录指向工作目录,并返回显示
#setwd(work_dir)
#getwd()

#-----------------------------------------------------------------------------
#                   创建全局变量（根据实际填写，非全部必须）
#-----------------------------------------------------------------------------
##### GDCquery #####
#project <- "TCGA-CESC"
projects <- TCGAbiolinks::getGDCprojects()$project_id
projects <- projects[grepl('^TCGA', projects, perl=TRUE)]
projects

data_category <- "Transcriptome Profiling" #RNA
data_type <- "Gene Expression Quantification"
workflow_type <- "STAR - Counts" #GDC在2022.3.28更新了，workflow.type只能选择STAR - Counts
legacy <- FALSE #这个参数主要是因为TCGA数据有两个入口可以下载，GDC Legacy Archive 和 GDC Data Portal，区别主要是注释参考基因组版本不同分别是：GDC Legacy Archive（hg19和GDC Data Portal（hg38）。参数默认为FALSE，下载GDC Data Portal（hg38）。这里小编的建议是，下载转录组层面的数据使用hg38，下载DNA层面的数据使用hg19，因为比如做SNP分析的时候很多数据库没有hg38版本的数据，都是hg19的。

##### GDCdownload #####
# TokenFile <- 
method <- 'api'
#download_dir <- paste0(work_dir,"/GDC/",gsub("-","_",project)) #默认则在默认位置
FilesChunk <- 6 #使用API下载大文件的时候，可以把文件分成几个小文件来下载，可以解决下载容易中断的问题。
  
##### GDCprepare #####  
#GDCprepare()将前面GDCquery（）的结果准备成R语言可处理的SE（SummarizedExperiment）文件
save <- TRUE
#SaveName <- paste0(download_dir,"_","RNAsqe_HTSeq_Counts",".rda") #文件名，如果没有设置，系统将默认设置
# summarizedExperiment <- FASLSE
  
#-----------------------------------------------------------------------------
#                   查询(根据实际情况调整参数)
#-----------------------------------------------------------------------------
query <- GDCquery(project = project,
                  data.category = data_category,
                  data.type = data_type,
                  workflow.type = workflow_type,
                  legacy = legacy 
                  )

#-----------------------------------------------------------------------------
#                   下载(根据实际情况调整参数)
#-----------------------------------------------------------------------------
GDCdownload(query = query,
            directory = download_dir,
            files.per.chunk = FilesChunk,
            method = method
            )

#-----------------------------------------------------------------------------
#                   预处理(根据实际情况调整参数)
#-----------------------------------------------------------------------------
data <- GDCprepare(query,
                   directory = download_dir,
                   save = save,
                   save.filename = SaveName
                   )

#-----------------------------------------------------------------------------
#                   提取
#-----------------------------------------------------------------------------
#获取表达矩阵
data_expr <- assay(data)

#查看表达矩阵的维度
dim(data_expr)

#创建表达矩阵文件路径expr_file
expr_file <- paste0(download_dir,"_","All_HTSeq_Counts",".txt")

#输出到表达矩阵文件中,file参数为""可以输出到console中
write.table(data_expr,file = expr_file,sep = "\t",row.names = T,quote = F)

#提取Gene表达量矩阵
gene_info = read.table("C:/Jagger/fannie/test/gene_info.txt",header = 1)
cat("Total Gene annotation:",dim(gene_info)[1])
gene_selected <- row.names(data_expr)%in%gene_info$Gene_id
gene_expr <- data_expr[gene_selected,]
cat("Total Gene Selected:",dim(gene_expr)[1])
gene_expr_file <- paste0(download_dir,"_","Gene_HTSeq_Counts",".txt")
write.table(gene_expr,file = gene_expr_file,sep="\t",row.names = T,quote = F)

#-----------------------------------------------------------------------------
#                   分析
#-----------------------------------------------------------------------------

##### Terminal模式中运行时使用此交互 #####
args <- commandArgs(TRUE)
working_dir <- args[1]
rawdata_file <- args[2]
sample_info_file <- args[3]

##### 编辑器中使用 #####
#工作目录
working_dir <- "C:/TCGA/Analysis"
#基因表达矩阵文件路径
rawdata_file <- "C:/TCGA/CESC/GDC/TCGA_CESC_Gene_HTSeq_Counts.txt"
#
#sample_info_file <- 
  
##### 设置工作目录并返回 #####
setwd(working_dir)
getwd()
  
#设置FDR和 fold change 阈值
fdr = 0.01
fold_change = 2

#读取基因表达矩阵
rawdata=read.table(rawdata_file,header = TRUE,stringsAsFactors = FALSE,row.names = 1)

#原始数据的基因数和样本量
dim(rawdata)

#针对基因的表达量进行过滤，过滤标准设置为：至少有25%的样本，基因的表达量大于2
quant <- apply(rawdata,1,quantile,0.75)
keep <- which((quant >= 2) == 1)
rawdata <- rawdata[keep,]
dim(rawdata)

#样分组信息，可以通过文件传入，也可以采用下面的简单方法，基于barcode判断是癌症还是正常组织
#sample_info =read.table(sample_info_file,header=TRUE,stringAsFactors=FALSE,row.names=1)

#基于样本的编号，根据第14位判断是癌症还是正常样本
#0为癌症样本，1为正常样本

#Barcode 实例： TCGA-W5-AA34-11A-11R-A41I-07

group <- factor(substr(colnames(rawdata),14,14))
summary(group)

#获得基因名称列表
genes=rownames(rawdata)

#制作DEGlist对象
y <- DGEList(counts=rawdata,genes=genes,group=group)
#计算有多少行
nrow(y)

#TMM标准化（基于count数据的标准化）
y <- calcNormFactors(y)

#分布估计（这步较慢）
y <- estimateCommonDisp(y,verbose=TRUE)
y <- estimateTagwiseDisp(y)

#差异分析检验（T检验）
et <- exactTest(y)

#针对FDR，fold change 对基因的显著性表达进行‘多重检验’
#对p value 进行BH检验
et$table$FDR <- p.adjust(et$table$PValue,method = "BH")
#summary显示共上调多少基因，下调多少基因（数据较大则重新调整FDR）
summary(de <- decideTestsDGE(et,adjust.method = "BH",p.value = fdr, lfc = log2(fold_change)))

#标识基因的上下调情况
et$table$regulate = as.factor(ifelse(et$table$FDR < fdr & abs(et$table$logFC) >=log2(fold_change),ifelse(et$table$logFC>log2(fold_change),'Up','Down'),'Normal'))

summary(et$table$regulate)

#绘制前调整右下角大小，过小会报错
plotMD(et)
print('et:')
et

#筛选出差异表达的基因
diff_expr_out <- et$table[et$table$regulate !='Normal',]

#输出结果

#保存结果
write.table(diff_expr_out,file = "DE_genes.txt",quote = FALSE, row.names = T,sep = "\t")

#绘制火山图
#设置FDR，Fold_Change 的辅助线
fdr_line = -log10(fdr)
fc_line = log2(fold_change)

#采用ggplot进行绘图
gp1 <- ggplot(et$table,aes(x=logFC,y=-log10(FDR),colour=regulate))+xlab("log2Fold Change")+ylab("-log10 FDR")+ylim(0,30)
#基于上下调关系，制定不同的颜色
gp2 <- gp1 + geom_point() +scale_color_manual(values = c("green","black","red"))
#增加两条辅助线
gp3 <- gp2 + geom_hline(yintercept = fdr_line,linetype=4)+geom_vline(xintercept = c(-fc_line,fc_line),linetype=4)

#run gp3可直接得到图
gp3

#保存图片
ggsave("DE_gene_Vocano.pdf",width = 16,height = 9)

##################
#绘制聚类热图
normData=y$pseudo.counts
heatmapData <- normData[rownames(diff_expr_out),]

#有些表达量为0，无法进行log转换，增加一个小背景：0.001
heatmapExp=log2(heatmapData+0.001)
heatmapMat=as.matrix(heatmapExp)

#输出的热图的储存途径
pdf(file = "DE_gene_heatmap.pdf",width = 60,height = 90)
par(oma=c(10,3,3,7))
heatmap.2(heatmapMat,col='greenred',trace='none')
dev.off()

#样品相关性绘图
pdf(file="gene_cor_cluster2.pdf",width = 60,height = 90)
par(oma=c(10,3,3,7))

#基于pearson相关性对样品进行聚类
heatmap(cor(normData))
dev.off()

#-----------------------------------------------------------------------------
#                   PPI
#-----------------------------------------------------------------------------
# 安装STRINGdb 软件包
BiocManager::install("STRINGdb")
BiocManager::install("rlist")
############################################################
library(STRINGdb)
library("rlist")
# 设置程序参数
work_dir <- "C:/TCGA/CESC/Analysis" 
deg_file <- "C:/TCGA/CESC/Analysis/DE_genes.txt"
setwd(work_dir)


# 获取物种的分类编号
# get_STRING_species(version="10", species_name=NULL) 
# NCBI数据库中9606 代表人类，700为可信度分值
#version已更新到11，但只有10不报错
string_db <- STRINGdb$new(version="10", species=9606,score_threshold=700, input_directory= work_dir)

# 读取差异表达的文件，获得差异表达基因列表
degs = read.table(deg_file,header=T,comment.char = "",check.names=F)
degs$gene <- rownames(degs)
head(degs)

# 查看有多少差异表达的基因需要分析 
cat("Total deg genes:", dim(degs)[1])

# 将基因的ID map 到string 数据库中， 不一定每个基因都能map上，很慢很慢
deg_mapped <- string_db$map( degs, "gene", removeUnmappedRows = TRUE )

# 查看有多少ID map 上了 
cat("Total String id mapped :", dim(deg_mapped)[1])

#-----------------------------------------------------------------------------
#                   GO KEGG
#-----------------------------------------------------------------------------
BiocManager::install("org.Hs.eg.db")
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("topGO")
install.packages("SparseM")
BiocManager::install("enrichMap")

library(org.Hs.eg.db)
library(clusterProfiler)
library(pathview)
require(pathview)
library(topGO)

#args <- commandArgs(TRUE)
#work_dir <-  args[1]
#deg_file <- args[2]


work_dir <- "C:/TCGA/Annotation"
deg_file <- "C:/TCGA/Analysis/DE_genes.txt"
setwd(work_dir)

# 读取差异表达的文件，获得差异表达基因列表
degs = read.table(deg_file,header=T,comment.char = "",check.names=F)
DEG_list <- rownames(degs)
head(DEG_list)
#GO 富集
# ont:  MF,BP,CC
ego <- enrichGO(gene          = DEG_list,
                keyType       = "ENSEMBL",
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

write.table(ego, file = "GO_enrichment_stat.txt",sep="\t", row.names =F, quote = F)


ego_results<-summary(ego)
ego_results
# 绘制相关图
pdf(file = "ego_barplot.pdf")
barplot(ego, showCategory=20, x = "GeneRatio")
dev.off()

# 点图
pdf(file = "ego_dotplot.pdf")
dotplot(ego,showCategory=20)
dev.off()

pdf(file = "enrich_map.pdf")
enrichMap(ego)
dev.off()

# 绘制topGO图
pdf(file = "topgo.pdf")
plotGOgraph(ego)
dev.off()


# KEGG pathway 富集
gene_ids <- bitr(DEG_list, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
kegg <- enrichKEGG(gene         = gene_ids$ENTREZID,
                   organism     = 'hsa',
                   pvalueCutoff = 0.05)
head(kegg)
write.table(kegg, file = "KEGG_enrichment_stat.txt",sep="\t", row.names =F, quote = F)

# 点图
png(file = "kegg_dotplot.png")
dotplot(kegg,title="Enrichment KEGG_dot")
dev.off()

# 合并表达信息和基因信息
deg_info <- degs['logFC']
deg_info['ENSEMBL'] <- rownames(deg_info)
gene_ids_merge <- merge(gene_ids, deg_info,by='ENSEMBL')
# 去掉ENTREZID 重复的基因
index<-duplicated(gene_ids_merge$ENTREZID)
gene_ids_merge<-gene_ids_merge[!index,]

# 提取FC的值，在map上进行颜色标注
map_ids <- as.matrix(gene_ids_merge['logFC'])
rownames(map_ids) <- gene_ids_merge$ENTREZID

                
#  绘制富集的pathway 以前三个为例
pathway_id <- "hsa04110"
map <- pathview(gene.data  = map_ids[,1],
                pathway.id = pathway_id,
                species    = "hsa", kegg.native = TRUE)
pathway_id <- "hsa03030"
map <- pathview(gene.data  = map_ids[,1],
                pathway.id = pathway_id,
                species    = "hsa", kegg.native = TRUE)
pathway_id <- "hsa04022"
map <- pathview(gene.data  = map_ids[,1],
                pathway.id = pathway_id,
                species    = "hsa", kegg.native = TRUE)
pathway_id <- "hsa04115"
map <- pathview(gene.data  = map_ids[,1],
                pathway.id = pathway_id,
                species    = "hsa", kegg.native = TRUE)
pathway_id <- "hsa04151"
map <- pathview(gene.data  = map_ids[,1],
                pathway.id = pathway_id,
                species    = "hsa", kegg.native = TRUE)
pathway_id <- "hsa04010"
map <- pathview(gene.data  = map_ids[,1],
                pathway.id = pathway_id,
                species    = "hsa", kegg.native = TRUE)

###第二种：MF
ego <- enrichGO(gene          = DEG_list,
                keyType       = "ENSEMBL",
                OrgDb         = org.Hs.eg.db,
                ont           = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

write.table(ego, file = "GO_enrichment_stat.txt",sep="\t", row.names =F, quote = F)
ego_results<-summary(ego)
ego_results
# 绘制相关图
pdf(file = "ego_barplot_MF.pdf")
barplot(ego, showCategory=20, x = "GeneRatio")
dev.off()
# 点图
pdf(file = "ego_dotplot_MF.pdf")
dotplot(ego,showCategory=20)
dev.off()
#效果不好
#pdf(file = "enrich_map.pdf")
#enrichMap(ego)
#dev.off()
# 绘制topGO图
pdf(file = "topgo_MF.pdf")
plotGOgraph(ego)
dev.off()

###第三种：CC
ego <- enrichGO(gene          = DEG_list,
                keyType       = "ENSEMBL",
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

write.table(ego, file = "GO_enrichment_stat.txt",sep="\t", row.names =F, quote = F)
ego_results<-summary(ego)
ego_results
# 绘制相关图
pdf(file = "ego_barplot_CC.pdf")
barplot(ego, showCategory=20, x = "GeneRatio")
dev.off()
# 点图
pdf(file = "ego_dotplot_CC.pdf")
dotplot(ego,showCategory=20)
dev.off()
#效果不好
#pdf(file = "enrich_map.pdf")
#enrichMap(ego)
#dev.off()
# 绘制topGO图
pdf(file = "topgo_CC.pdf")
plotGOgraph(ego)
dev.off()

