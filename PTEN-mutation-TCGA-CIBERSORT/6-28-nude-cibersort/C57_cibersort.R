
### /home/wus/scRNA_program_blackmouseRNA/Nude/6_featureCounts/featurecount_gene_name/6-28_FPKM/

##/home/wus/scRNA_program_blackmouseRNA/blackmouse/6_featureCounts/featurecount_gene_name/6-29_TPM/
rm(list=ls())
options(stringsAsFactors = F) 
library(tidyverse) # ggplot2 stringer dplyr tidyr readr purrr  tibble forcats
###读入文件
a1 <- read.table('m.count', header = T)
a1 <- a1[!duplicated(a1$Geneid),]
### counts矩阵的构建
counts <- a1[,7:ncol(a1)] #截取样本基因表达量的counts部分作为counts 
rownames(counts) <- a1$Geneid #将基因名作为行名
### 从featurecounts 原始输出文件counts.txt中提取Geneid、Length(转录本长度)，
geneid_efflen <- subset(a1,select = c("Geneid","Length"))
       colnames(geneid_efflen) <- c("geneid","efflen")  
geneid_efflen_fc <- geneid_efflen #用于之后比较

### 取出counts中geneid的对应的efflen
dim(geneid_efflen)
efflen <- geneid_efflen[match(rownames(counts),
                              geneid_efflen$geneid),
                        "efflen"]

### 计算 TPM
#TPM (Transcripts Per Kilobase Million)  每千个碱基的转录每百万映射读取的Transcripts
counts2TPM <- function(count=count, efflength=efflen){
  RPK <- count/(efflength/1000)       #每千碱基reads (“per million” scaling factor) 长度标准化
  PMSC_rpk <- sum(RPK)/1e6        #RPK的每百万缩放因子 (“per million” scaling factor ) 深度标准化
  RPK/PMSC_rpk                    
  }  
tpm <- as.data.frame(apply(counts,2,counts2TPM))
colSums(tpm)

colnames(tpm) <- c('c1', 'c2','c3','p1','p2','p3')
#write.csv(tpm,'M_Nude.csv')
#write.table(tpm,'M_Nude.txt',row.names = T,quote = F,sep = "\t")
df <- tibble::rownames_to_column(tpm,"gene")
write.table(df,'C57.txt',row.names = F,quote = F,sep = "\t")

############# cibersort
##/home/wus/016-2023-6-15_PTEN-mutation-TCGA-CIBERSORT/src



library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)
source('/home/wus/scRNA_program_blackmouseRNA/Nude/6_featureCounts/featurecount_gene_name/6-28_FPKM/Cibersort.R')
m_Nude.file <- "./C57.txt"
mice.file<-'/home/wus/scRNA_program_blackmouseRNA/Nude/6_featureCounts/featurecount_gene_name/6-28_FPKM/mics.txt'
CIBERSORT(mice.file ,m_Nude.file, perm = 50, QN = F)
#m_Nude.results<-CIBERSORT(mice.file ,m_Nude.file, perm = 50, QN = F)  
# perm置换次数=1000
# QN如果是芯片设置为T，如果是测序就设置为F

#write.csv(m_Nude.results, "./Output/C57_CIBERSORT_Results.csv")


############# box
library(tidyverse)
library(ggplot2)
library(ggpubr)

#### WIN上操作，删除不需要的行并添加列 列名

a<-read.csv('/home/wus/scRNA_program_blackmouseRNA/blackmouse/6_featureCounts/featurecount_gene_name/6-29_TPM/CIBERSORT-Results.csv')
head(a)
#a[, ] <- sapply(a[, ], as.character)


library(tidyr)
chji_long <- gather(a, key = "CELL", value = "SCORE",
                    -'SAMPLE', -'TYPE')
 bxp <- ggboxplot(
    chji_long, x="CELL", y="SCORE", color = "TYPE",bxp.errorbar =TRUE,outlier.shape = NA)

bxp<-bxp+ theme(axis.text.x = element_text(angle = 90, hjust = 1))

bxp
ggsave(bxp,filename = "/home/wus/scRNA_program_blackmouseRNA/blackmouse/6_featureCounts/featurecount_gene_name/6-29_TPM/c57-cibersort.png",width = 6,height = 5,dpi=300)


#### test1 
# bxp<-bxp+stat_compare_means(aes(group=CELL),label = "p.format",method = "t.test",angle=45,label.y.npc='center')+ 
#  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))

#### test2
 bxp <- ggboxplot(
    chji_long, x="CELL", y="SCORE", color = "TYPE",bxp.errorbar =TRUE,outlier.shape = NA,
    palette = "jco")
#bxp + stat_compare_means(aes(group = TYPE),hide.ns=T,position='jitter',label = "p.format")+ 
# theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))

bxp
ggsave(bxp,filename = "/home/wus/scRNA_program_blackmouseRNA/blackmouse/6_featureCounts/featurecount_gene_name/6-29_TPM/c57-cibersort.png",width = 6,height = 5,dpi=300)
