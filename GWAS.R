----------------
#超算
----------------
Y3W#uSrF

----------------
scRNA
----------------
 library (stringr)
/dssg/home/acct-medwshuai/medwshuai/2024-6-5_GWAS/snp_result.txt
d<-read.table('/dssg/home/acct-medwshuai/medwshuai/2024-6-5_GWAS/THIN_NORMAL_FILTER_MAF_1prc_INFO_0.8.txt',header=T,sep='\t')

d[c(' First ', ' Last ')] <- str_split_fixed(d$MARKER, ':', 2)

s<-read.table('/dssg/home/acct-medwshuai/medwshuai/2024-6-5_GWAS/rs.list',header=F,sep='\t')

m1<-merge(d,s,by.x=' First ',by.y='V1')

write.table(m1,'THIN.txt',sep='\t')
write.table(d,'d.txt')


d<-read.table('/dssg/home/acct-medwshuai/medwshuai/2024-6-5_GWAS/OBESE_NORMAL_FILTER_MAF_1prc_INFO_0.8.txt',header=T,sep='\t')

d[c(' First ', ' Last ')] <- str_split_fixed(d$MARKER, ':', 2)



m2<-merge(d,s,by.x=' First ',by.y='V1')
write.table(m2,'OBESE.txt',sep='\t')
write.table(d,'d_OBESE.txt',sep='\t')


library(qqman)
d<-read.table('/dssg/home/acct-medwshuai/medwshuai/2024-6-5_GWAS/OBESE_NORMAL_FILTER_MAF_1prc_INFO_0.8.txt',header=T,sep='\t')


d[c('SNP', 'BP','OTHERS')] <- str_split_fixed(d$MARKER, ':', 3)

manhattan(d)


-----------
ALK
-----------
