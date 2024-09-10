#########最终出图 8-27
###整体看（小提琴图，去点+feature）
###全面注释
rm(list=ls())
load(file='/home/wus/c_t_all_8_23-all-feature.RData')
mac_sce1=pancreas.integrated





######实际使用
VlnPlot(mac_sce1, features = c('Ptprc'),pt.size=0)
#大致了解
VlnPlot(mac_sce1, features = c('Adgre1','Cd68','Itgam','Cd14','S100a8','S100a9','Cd3e','Cd3d','Cd3g','Cd79a'),pt.size=0)
FeaturePlot(mac_sce1, features = c('Adgre1','Cd68','Itgam','Cd14','S100a8','S100a9','Cd3e','Cd3d','Cd3g','Cd79a'))
#mac neu
VlnPlot(mac_sce1, features = c('Adgre1','Cd68','Itgam','Cd14','S100a8','S100a9'),pt.size=0)
FeaturePlot(mac_sce1, features = c('Adgre1','Cd68','Itgam','Cd14','S100a8','S100a9'))

# 经典的bug+流式+realtime
VlnPlot(mac_sce1, features = c('Il1b','Nos2','Cd86','Tnf','Mrc1','Arg1','Il10','Il6','H2-Ab1'),pt.size=0)
FeaturePlot(mac_sce1, features =c('Il1b','Nos2','Cd86','Tnf','Mrc1','Arg1','Il10'))

#neu
FeaturePlot(mac_sce1, features =c('Fcgr3a','Csf3r','Cclr2','Ly6g','Trem1','Bst1','Il1r2','Ptgs2','Ctsg','Itgax','Mmp9','Pilra'))
VlnPlot(mac_sce1, features = c('Fcgr3a','Csf3r','Cclr2','Ly6g','Trem1','Bst1','Il1r2','Ptgs2','Ctsg','Itgax','Mmp9','Pilra'),pt.size=0)
#neu vs macs
VlnPlot(mac_sce1, features = c('Csf3r','Csf1r','Itgam','Ly6g'),pt.size=0)
FeaturePlot(mac_sce1, features = c('Csf3r','Csf1r','Itgam','Ly6g'))

# B T
VlnPlot(mac_sce1, features = c('Cd3e','Cd3d','Cd3g','Cd79a'),pt.size=0)

#fib #内皮 #erythroid 污染
VlnPlot(mac_sce1, features = c('Dcn','Adgrf5','Clec4g','Ehd3','Kdr','Gypa','Ly76'),pt.size=0)


VlnPlot(mac_sce1, features = c('H2-Aa','Itgam','H2-Ab1','Clec9a','Xcr1','Itgax','Sirpa','Fcgr1','Irf8','Flt3','Zbtb46','Batf3','Itgae'),pt.size=0)






####统一
rm(list=ls())
load(file = '/home/wus/for_seurat/8_2/CT_macs_step2_8_28.Rdata')
mac_sce1=pancreas.integrated

#mac neu
VlnPlot(mac_sce1, features = c('Adgre1','Cd68','Itgam','Cd14','S100a8','S100a9'),pt.size=0)
FeaturePlot(mac_sce1, features = c('Adgre1','Cd68','Itgam','Cd14','S100a8','S100a9'))
#neu vs macs
VlnPlot(mac_sce1, features = c('Csf3r','Csf1r','Itgam','Ly6g'),pt.size=0)
FeaturePlot(mac_sce1, features = c('Csf3r','Csf1r','Itgam','Ly6g'))

# 经典的bug+流式+realtime
VlnPlot(mac_sce1, features = c('Il1b','Nos2','Cd86','Tnf','Mrc1','Arg1','Il10','Il6','H2-Ab1'),pt.size=0)
FeaturePlot(mac_sce1, features =c('Il1b','Nos2','Cd86','Tnf','Mrc1','Arg1','Il10'))

#neu
FeaturePlot(mac_sce1, features =c('Fcgr3a','Csf3r','Cclr2','Ly6g','Trem1','Bst1','Il1r2','Ptgs2','Ctsg','Itgax','Mmp9','Pilra'))
VlnPlot(mac_sce1, features = c('Fcgr3a','Csf3r','Cclr2','Ly6g','Trem1','Bst1','Il1r2','Ptgs2','Ctsg','Itgax','Mmp9','Pilra'),pt.size=0)

#两版整合去重
M1 = c('Il1r1','Nos2','Tlr2','Tlr4','Cd80','Cd86','Csf2','Il1b','Il18','Il12b','Ccr7','Il6','Cxcl2','Ifng','Cd38','Tnf','Socs3','Ptgs2','Nfkbiz','Lrf5','Gpr18','Fpr2','Cxcl9','Cxcl10','Azin1','Cd68','Ccl20','Ccl4','Ccl5','Cxcl11','Irf1','Irf5','Il1a','Ido1','Il12a','Il23a')
M2 = c('Myc','Cd83','Mrc1','Arg1','Egr2','Ccr7','Chil4','Pparg','Cd163','Clec7a','Il10','Il4','Irf4','Pdgfb','Stat6','Cd40','Csf1','Tlr7','Ccl13','Ccl17','Ccl18','Ccl22','Ccl24','Cd86','Vegfa','Vegfb','Vegfc','Vegfd','Mmp9','Fn1','Egf','Lyve1','Mmp14','Mmp19','Cd276','Fasl','Ctsa','Ctsb','Ctsc','Ctsd','Mgl2','Ear11','Clec10a','Retnla','Ccl10')
all_markers<-c(M1,M2)
DotPlot(mac_sce1, features = unique(all_markers),group.by = "seurat_clusters",col.max=1,col.min=-1,cols = c( "green",'black'))+RotatedAxis()

#DC
VlnPlot(mac_sce1, features = c('H2-Aa','H2-Ab1','Itgam','Clec9a','Xcr1','Itgax','Sirpa','Fcgr1','Irf8','Flt3','Zbtb46','Batf3','Itgae'),pt.size=0)

##macs T
rm(list=ls())
load(file = '/home/wus/for_seurat/8_2/T_macs_8_26.Rdata')
mac_sce1=T_t_cell

DimPlot(mac_sce1, reduction = "umap",label = TRUE, 
    label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")

#
FeaturePlot(mac_sce1, features = c('Cd3e','Nkg7','Cd20','Cd19','Ly6g','S100a9','S100a8','Adgre1','Cd68','Csfr1'))
VlnPlot(mac_sce1,features = c('Cd3e','Nkg7','Cd20','Cd19','Ly6g','S100a9','S100a8','Adgre1','Cd68','Csfr1'))
#macro
FeaturePlot(mac_sce1, features = c('Adgre1','Fcgr1','Itgam','Itgax'))
#neu
FeaturePlot(mac_sce1, features =c('Fcgr3a','Csf3r','Cclr2','Ly6g','Trem1','Bst1','Il1r2','Ptgs2','Ctsg','Itgax','Mmp9','Pilra'))
VlnPlot(mac_sce1, features = c('Fcgr3a','Csf3r','Cclr2','Ly6g','Trem1','Bst1','Il1r2','Ptgs2','Ctsg','Itgax','Mmp9','Pilra'))
#neu vs macs
VlnPlot(mac_sce1, features = c('Csf3r','Csf1r'))
FeaturePlot(mac_sce1, features = c('Csf3r','Csf1r'))


#命名与画图

cluster2celltype <- c(
                      "0"='M2-like',
                      "1"="M2-like", 
                      "2"="neu", 
                      "3"= "neu", 
                      "4"= "M1-like",
                      '5'='M0-like',
                      '6'='DC',
                      '7'='neu',
                      '8'='DC',
                      '9'='proliferation'

                      )
mac_sce1_rename=mac_sce1
new_id<-c('M2-like',"M2-like","neu",'neu',"M1-like",'M0-like','dc','neu','dc','proliferation')
names(new_id)<-levels(mac_sce1_rename)
mac_sce1_rename<-RenameIdents(mac_sce1_rename,new_id)
#DimPlot(mac_sce1_rename, reduction = 'umap', 
#        label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(mac_sce1_rename, reduction = "umap",label = TRUE, 
    label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")

save(mac_sce1_rename,file = '/home/wus/for_seurat/8_2/CT_real_macs_rename_8_26.Rdata')


#取子集

#sce=mac_sce1_rename

#mac_sce1_re = sce[,sce@meta.data$seurat_clusters %in% c(0,1,2,4)]

#DimPlot(mac_sce1_re, reduction = 'umap', 
#        label = TRUE, pt.size = 0.5) + NoLegend()

#save(mac_sce1_re,file = '/home/wus/for_seurat/8_2/C_macs_rename_8_22.Rdata')

rm(list=ls())
load(file = '/home/wus/for_seurat/8_2//home/wus/for_seurat/8_2/CT_real_macs_rename_8_26.Rdata')
C_ref=mac_sce1_rename
dim(C_ref)
print(table(Idents(C_ref)))
C_macs_Cell_Num=table(Idents(C_ref))

#macs vs neu
genes_to_check = c('Adgre1','Cd68','Csfr1','Mertk','Csfr3','S100a8','S100a9','Lcn2','Mmp9','Mmp8')
DotPlot(C_ref, features = genes_to_check,
             assay='RNA')
#M1 vs M2
genes_to_check = c('Pf4','Apoe','Mrc1','C1qb','C1qa','C1qc','Hp','Chil3','S100a4','Traf1','Lsp1','Il1b','F10')
DotPlot(C_ref, features = genes_to_check,
             assay='RNA')

#合并
genes_to_check = c('Adgre1','Cd68','Csfr1','Mertk','Pf4','Apoe','Mrc1','C1qb','C1qa','C1qc','Hp','Chil3','S100a4','H2-Aa','H2-Ab1','Ly6c','Lyz2','Lsp1','Il1b','F10','Csfr3','S100a8','S100a9','Lcn2','Mmp9','Mmp8')
DotPlot(C_ref, features = genes_to_check,
             assay='RNA')

#M1-流式
VlnPlot(mac_sce1, features = c('Il1b','Il6','Ccl2','Ccl5'))
VlnPlot(mac_sce1, features = c('Arg1','Mrc1','H2-Ab1','Cd163'))


##### C
##macs C
rm(list=ls())
load(file = '/home/wus/for_seurat/8_2/C_myeloid_8_23.Rdata')
mac_sce1=pancreas.query

DimPlot(mac_sce1, reduction = "ref.umap",label = TRUE, 
    label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")

#
FeaturePlot(mac_sce1, features = c('Cd3e','Nkg7','Cd20','Cd19','Ly6g','S100a9','S100a8','Adgre1','Cd68','Csfr1'))
VlnPlot(mac_sce1,features = c('Cd3e','Nkg7','Cd20','Cd19','Ly6g','S100a9','S100a8','Adgre1','Cd68','Csfr1'))
#macro
FeaturePlot(mac_sce1, features = c('Adgre1','Fcgr1','Itgam','Itgax'))
#neu
FeaturePlot(mac_sce1, features =c('Fcgr3a','Csf3r','Cclr2','Ly6g','Trem1','Bst1','Il1r2','Ptgs2','Ctsg','Itgax','Mmp9','Pilra'))
VlnPlot(mac_sce1, features = c('Fcgr3a','Csf3r','Cclr2','Ly6g','Trem1','Bst1','Il1r2','Ptgs2','Ctsg','Itgax','Mmp9','Pilra'))
#neu vs macs
VlnPlot(mac_sce1, features = c('Csf3r','Csf1r'))
FeaturePlot(mac_sce1, features = c('Csf3r','Csf1r'))


#命名与画图

cluster2celltype <- c(
                      "0"='M2-like',
                      "1"="M2-like", 
                      "2"="M1-like", 
                      "3"= "neu", 
                      "4"= "neu",
                      '5'='M2-like',
                      '6'='dc',
                      '7'='dc'

                      )
mac_sce1_rename=mac_sce1
new_id<-c('M2-like',"M2-like","M1-like","neu",'neu','M2-like','dc','dc')
names(new_id)<-levels(mac_sce1_rename)
mac_sce1_rename<-RenameIdents(mac_sce1_rename,new_id)
#DimPlot(mac_sce1_rename, reduction = 'umap', 
#        label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(mac_sce1_rename, reduction = "ref.umap",label = TRUE, 
    label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")

save(mac_sce1_rename,file = '/home/wus/for_seurat/8_2/C_real_macs_rename_8_22.Rdata')


#取子集

#sce=mac_sce1_rename

#mac_sce1_re = sce[,sce@meta.data$seurat_clusters %in% c(0,1,2,4)]

#DimPlot(mac_sce1_re, reduction = 'umap', 
#        label = TRUE, pt.size = 0.5) + NoLegend()

#save(mac_sce1_re,file = '/home/wus/for_seurat/8_2/C_macs_rename_8_22.Rdata')

rm(list=ls())
load(file = '/home/wus/for_seurat/8_2/C_real_macs_rename_8_22.Rdata')
C_ref=mac_sce1_rename
dim(C_ref)
print(table(Idents(C_ref)))
C_macs_Cell_Num=table(Idents(C_ref))

#macs vs neu
genes_to_check = c('Adgre1','Cd68','Csfr1','Mertk','Csfr3','S100a8','S100a9','Lcn2','Mmp9','Mmp8')
DotPlot(C_ref, features = genes_to_check,
             assay='RNA')
#M1 vs M2
genes_to_check = c('Pf4','Apoe','Mrc1','C1qb','C1qa','C1qc','Hp','Chil3','S100a4','Traf1','Lsp1','Il1b','F10')
DotPlot(C_ref, features = genes_to_check,
             assay='RNA')

#合并
genes_to_check = c('Adgre1','Cd68','Csfr1','Mertk','Pf4','Apoe','Mrc1','C1qb','C1qa','C1qc','Hp','Chil3','S100a4','H2-Aa','H2-Ab1','Ly6c','Lyz2','Lsp1','Il1b','F10','Csfr3','S100a8','S100a9','Lcn2','Mmp9','Mmp8')
DotPlot(C_ref, features = genes_to_check,
             assay='RNA')

#M1-流式
VlnPlot(mac_sce1, features = c('Il1b','Il6','Ccl2','Ccl5'))
VlnPlot(mac_sce1, features = c('Arg1','Mrc1','H2-Ab1','Cd163'))



#####8-26

####t cell
##T
### 判断+ 作图
rm(list=ls())
load(file = '/home/wus/for_seurat/8_2/T_t_cell_8_26.Rdata')
T_ref=T_t_cell

#DoHeatmap(T_ref, features = top20$gene) + NoLegend()



genes_to_check = c('Cd3e','Cd4','Cd8a','Cd8b1','Top2a','Mki67','Klrb1c','Gzma')
DotPlot(T_ref, features = genes_to_check,assay='RNA')



VlnPlot(T_ref, features = c('Cd3e','Cd4','Cd8a','Foxp3','Mki67','Mzb1','Ncr1','Ncam1','Klrb1c'),pt.size=0)
FeaturePlot(T_ref, features = c('Cd3e','Cd4','Cd8a','Foxp3','Mki67','Mzb1','Ncr1','Ncam1','Klrb1c'))
#
cluster2celltype <- c(
                      "0"='CD4',
                      "1"="γσT", 
                      "2"="CD8", 
                      "3"= "proliferation",
                      '4'='NK',
                      '5'='CD8',
                      '6'='other'
                      )
T_ref_rename=T_ref
new_id<-c('CD4',"γσT", "CD8","proliferation","NK",'CD8','other')
names(new_id)<-levels(T_ref_rename)
T_ref_rename<-RenameIdents(T_ref_rename,new_id)

DimPlot(T_ref_rename, reduction = "ref.umap",label = TRUE, 
    label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")

genes_to_check = c('Cd3e','Cd4','Cd8a','Cd8b1','Top2a','Mki67','Klrb1c','Gzma')
DotPlot(T_ref_rename, features = genes_to_check,assay='RNA')

save(T_ref_rename,file = '/home/wus/for_seurat/8_2/T_t-cell_rename_8_24.Rdata')


#取子集
rm(list=ls())
load(file = '/home/wus/for_seurat/8_2/T_t-cell_rename_8_24.Rdata')
sce=T_ref_rename

T_ref_re = sce[,sce@meta.data$seurat_clusters %in% c(0,1,2,3,4)]

DimPlot(T_ref_re, reduction = "umap",label = TRUE, 
    label.size = 3, repel = TRUE) + NoLegend() + ggtitle("T")

genes_to_check = c('Cd3e','Cd4','Cd8a','Cd8b1','Top2a','Mki67','Klrb1c','Gzma')
DotPlot(T_ref_re, features = genes_to_check,assay='RNA')

save(T_ref_re,file = '/home/wus/for_seurat/8_2/T_macs_rename_8_24.Rdata')

rm(list=ls())
load(file = '/home/wus/for_seurat/8_2/T_macs_rename_8_24.Rdata')
genes_to_check = c('Cd3e','Cd4','Cd8a','Cd8b1','Top2a','Mki67','Klrb1c','Gzma')
DotPlot(T_ref_re, features = genes_to_check,assay='RNA')+ coord_flip()



##C
### 判断+ 作图
rm(list=ls())
load(file = '/home/wus/for_seurat/8_2/C_t_cell_8_26.Rdata')
T_ref=C_t_cell

#DoHeatmap(T_ref, features = top20$gene) + NoLegend()

genes_to_check = c('Cd3e','Cd4','Cd8a','Cd8b1','Top2a','Mki67','Klrb1c','Gzma')
DotPlot(T_ref, features = genes_to_check,assay='RNA')

#
cluster2celltype <- c(
                      "0"='CD4',
                      "1"="γσT", 
                      "2"="CD8", 
                      "3"= "proliferation",
                      '4'='NK',
                      '5'='CD8',
                      '6'='other'
                      )
T_ref_rename=T_ref
new_id<-c('CD4',"γσT", "CD8","proliferation","NK",'CD8','other')
names(new_id)<-levels(T_ref_rename)
T_ref_rename<-RenameIdents(T_ref_rename,new_id)

DimPlot(T_ref_rename, reduction = "ref.umap",label = TRUE, 
    label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")

genes_to_check = c('Cd3e','Cd4','Cd8a','Cd8b1','Top2a','Mki67','Klrb1c','Gzma')
DotPlot(T_ref_rename, features = genes_to_check,assay='RNA')

save(T_ref_rename,file = '/home/wus/for_seurat/8_2/T_t-cell_rename_8_24.Rdata')


#取子集
rm(list=ls())
load(file = '/home/wus/for_seurat/8_2/T_t-cell_rename_8_24.Rdata')
sce=T_ref_rename

T_ref_re = sce[,sce@meta.data$seurat_clusters %in% c(0,1,2,3,4)]

DimPlot(T_ref_re, reduction = "umap",label = TRUE, 
    label.size = 3, repel = TRUE) + NoLegend() + ggtitle("C")

genes_to_check = c('Cd3e','Cd4','Cd8a','Cd8b1','Top2a','Mki67','Klrb1c','Gzma')
DotPlot(T_ref_re, features = genes_to_check,assay='RNA')

save(T_ref_re,file = '/home/wus/for_seurat/8_2/T_macs_rename_8_24.Rdata')

#rm(list=ls())
#load(file = '/home/wus/for_seurat/8_2/T_macs_rename_8_24.Rdata')
#genes_to_check = c('Cd3e','Cd4','Cd8a','Cd8b1','Top2a','Mki67','Klrb1c','Gzma')
#DotPlot(T_ref_re, features = genes_to_check,assay='RNA')