测试过的marker




#######marker
FeaturePlot(cd4_sce1, features = c('Cd3e','Cd4','Cd8a','Foxp3','Mki67','Mzb1','Ncr1','Ncam1','Klrb1c'))
VlnPlot(cd4_sce1, features = c('Cd3e','Cd4','Cd8a','Foxp3','Mki67','Mzb1','Ncr1','Ncam1'))
VlnPlot(cd4_sce1, features = c('Cd4','Cd8a','Klrb1c','Gzmk','Il7r'))
VlnPlot(cd4_sce1, features = c('Lef1','Ccr7','Tcf7','Il7r'))
#5.细分T 
#Naive SELL, CCR7, SOCS3, MYC, LEF1, TSHZ2, MAL, NOSIP, TCF7
VlnPlot(cd4_sce1, features = c('Sell','Ccr7','Sosc3','Myc','Lef1','Tshz2','Mal','Nosip','Tcf7'))
#Treg    FOXP3, RTKN2, CTLA4, FANK1, IL2RA, ICA1, ARID5B, CORO1B, GBP2
VlnPlot(cd4_sce1, features = c('Foxp3','Rtkn2','Ctla4','Fank1','Il2ra','Ica1','Arid5b','Coro1b','Gbp2'))
VlnPlot(cd4_sce1, features = c('Ikzf2','Foxp3','Ctla4','Entbd1','Il2ra','Itgae','Il1r1'))
#nk
VlnPlot(cd4_sce1, features = c('Klrd1','Gnly','Nkg7','Gzma','Hopx','Klrb1c'))
#nk t cell
VlnPlot(cd4_sce1, features = c('Il2rb','Ncam1','Il12rb2','Cxcr4','Il17ra'))
#tem
VlnPlot(cd4_sce1, features = c('Il7r','Cxcr3'))
#tcm
VlnPlot(cd4_sce1, features = c('Sell','Ccr7','Lef1','Tcf7'))





#cd45 Ptprc
VlnPlot(mac_sce1, features = c('Ptprc'),pt.size=0)
#macs
VlnPlot(T_obj, features = c('Adgre1','Cd68','Itgam','Cd14'))
FeaturePlot(T_obj, features = c('Adgre1','Cd68','Itgam','Cd14'))
#Mdsc
##Itgam/CD11b, Ly6c2/Ly6C, 
#M-MDSC  cd11b ly6c
VlnPlot(mac_sce1, features = c('Itgam','Ly6c2'))
#G-MDSC cd11b ly6g
VlnPlot(mac_sce1, features = c('Itgam','Ly6g'))
#cd11b ly6c ly6g
VlnPlot(mac_sce1, features = c('Itgam','Ly6c2','Ly6g'))
FeaturePlot(mac_sce1, features = c('Itgam','Ly6c2','Ly6g'))


#DC  'H2-Aa','Itgam','H2-Ab1'=dc  'Clec9a','Xcr1'=dc1  'Itgax','Sirpa'=cd2
VlnPlot(mac_sce1, features = c('H2-Aa','H2-Ab1','Itgam','Clec9a','Xcr1','Itgax','Sirpa','Fcgr1','Irf8','Flt3','Zbtb46','Batf3','Itgae'))

#cDC1 Flt3 Irf8 Xcr1
VlnPlot(mac_sce1, features = c('Flt3','Irf8','Xcr1'))
#cDC1  Irf8 Batf3,Nfil3,Bcl6
VlnPlot(mac_sce1, features = c('Flt3','Nfil3','Batf3','Bcl6','Irf8'))
#cDC2 Flt3 Cd209
VlnPlot(mac_sce1, features = c('Flt3','Cd209a'))

#4.其他免疫细胞和污染

VlnPlot(mac_sce1, features = c('Dcn'))
#T
VlnPlot(mac_sce1, features = c('Cd3c','Cd3d','Cd3e','Cd3g'))
#B
VlnPlot(mac_sce1, features = c('Cd19','Cd79a'))
#内皮
VlnPlot(mac_sce1, features = c('Adgrf5','Clec4g','Ehd3','Kdr'))
#erythroid
FeaturePlot(mac_sce1, features = c('Gypa','Ly76'))



### t cell


VlnPlot(cd4_sce1, features = c('Cd4','Cd8a','Klrb1c','Gzmk','Il7r'))
VlnPlot(cd4_sce1, features = c('Lef1','Ccr7','Tcf7','Il7r'))

FeaturePlot(cd4_sce1, features = c('Cd3e','Cd4','Cd8a','Foxp3','Mki67','Mzb1','Ncr1','Ncam1'))
VlnPlot(cd4_sce1, features = c('Cd3e','Cd4','Cd8a','Foxp3','Mki67','Mzb1','Ncr1','Ncam1'))

VlnPlot(cd4_sce1, features = c('Foxp3','Il2ra'))



VlnPlot(cd4_sce1, features = c('Cd3e','Cd4','Cd8a','Cd20','Cd19','Nkg7'))



FeaturePlot(pancreas.query, features = c('Cd3e','Cd4','Cd8a','Foxp3','Mki67','Mzb1','Ncr1','Ncam1'))
VlnPlot(cd4_sce1, features = c('Cd3e','Cd4','Cd8a','Foxp3','Mki67','Mzb1','Ncr1','Ncam1'))

VlnPlot(cd4_sce1, features = c('Cd4','Cd8a','Klrb1c','Gzmk','Il7r'))
VlnPlot(cd4_sce1, features = c('Lef1','Ccr7','Tcf7','Il7r'))










#2，再确定M1 M2，用最经典的
VlnPlot(mac_sce1, features = c('Il1b','Nos2','Cd68','Cd80','Tnf','Cd86','Mrc1','Arg1','Il10'))
FeaturePlot(mac_sce1, features = c('Il1b','Nos2','Cd68','Cd80','Tnf','Cd86','Mrc1','Arg1','Il10'))
#2.5，如果不行，换最全面的
#M1
M1 = c('Il1r1','Nos2','Tlr2','Tlr4','Cd80','Cd86','Csf2','Il1b','Il18','Il12b','Ccr7','Il6','Cxcl2','Ifng','Cd38','Tnf','Socs3','Ptgs2','Nfkbiz','Lrf5','Gpr18','Fpr2','Cxcl9','Cxcl10','Azin1')
#M2
M2 = c('Myc','Cd83','Mrc1','Arg1','Egr2','Ccr7','Chil4','Pparg','Cd163','Clec7a','Il10','Il4','Irf4','Pdgfb','Stat6','Chil3','Cd40','Csf1','Tlr7')
all_markers<-c(M1,M2)
DotPlot(mac_sce1, features = unique(all_markers),group.by = "seurat_clusters",col.max=1,col.min=-1,cols = c( "green",'black'))+RotatedAxis()

##另一版
M1 = c('Nos2','Tlr4','Cd80','Cd68','Csf2','Tnf','Il1b','Ccr7','Il6','Cxcl2','Ifng','Cd38','Ccl20','Ccl4','Ccl5','Cd80','Cxcl10','Cxcl11','Cxcl9','Il1b','Nos2','Tnf','Il6','Irf1','Irf5','Il1a','Ido1','Il12a','Il12b','Il23a')
#M2
M2 = c('Arg1','Ccl13','Ccl17','Ccl18','Ccl22','Ccl24','Cd163','Cd86','Il10','Mrc1','Vegfa','Vegfb','Vegfc','Vegfd','Mmp9','Irf4','Fn1','Egf','Lyve1','Mmp14','Mmp19','Cd276','Fasl','Ctsa','Ctsb','Ctsc','Ctsd','Cd83','Mrc1','Ccr7','Arg1','Mgl2','Ccl17','Ear11','Pparg','Cd163','Clec10a','Clec7a','Retnla','Ccl10','Il10','Il4','Irf4','Pdgfb','Stat6','Chil3')

#monocytes=c('Cfp','Apobec3a','Lyz','Hla-dra','Cd7','Tet2','Cd40','Dysf','Mefv')
all_markers<-c(M1,M2)
DotPlot(mac_sce1, features = unique(all_markers),group.by = "seurat_clusters",col.max=1,col.min=-1,cols = c( "green",'black'))+RotatedAxis()

#两版整合去重
M1 = c('Il1r1','Nos2','Tlr2','Tlr4','Cd80','Cd86','Csf2','Il1b','Il18','Il12b','Ccr7','Il6','Cxcl2','Ifng','Cd38','Tnf','Socs3','Ptgs2','Nfkbiz','Lrf5','Gpr18','Fpr2','Cxcl9','Cxcl10','Azin1','Cd68','Ccl20','Ccl4','Ccl5','Cxcl11','Irf1','Irf5','Il1a','Ido1','Il12a','Il23a')
M2 = c('Myc','Cd83','Mrc1','Arg1','Egr2','Ccr7','Chil4','Pparg','Cd163','Clec7a','Il10','Il4','Irf4','Pdgfb','Stat6','Chil3','Cd40','Csf1','Tlr7','Ccl13','Ccl17','Ccl18','Ccl22','Ccl24','Cd86','Vegfa','Vegfb','Vegfc','Vegfd','Mmp9','Fn1','Egf','Lyve1','Mmp14','Mmp19','Cd276','Fasl','Ctsa','Ctsb','Ctsc','Ctsd','Mgl2','Ear11','Clec10a','Retnla','Ccl10')
all_markers<-c(M1,M2)
DotPlot(mac_sce1, features = unique(all_markers),group.by = "seurat_clusters",col.max=1,col.min=-1,cols = c( "green",'black'))+RotatedAxis()

#3，排除混杂有neu，mono的干扰
#neu
FeaturePlot(mac_sce1, features = c('Itgam','Ly6g','S100a9','S100a8'))
VlnPlot(mac_sce1,features = c('Itgam','Ly6g','S100a9','S100a8'))

FeaturePlot(mac_sce1, features = c('Itgam','Ly6g','Cd68','Cd80','Cd86','Mertk'))
VlnPlot(mac_sce1,features = c('Itgam','Ly6g','Cd68','Cd80','Cd86','Mertk'))
#monocyte
FeaturePlot(mac_sce1, features = c('Itgam','Cx3cr1','S100a4','Csf1r','Ly6c2'))
VlnPlot(mac_sce1, features = c('Itgam','Cx3cr1','S100a4','Csf1r','Ly6c2'),pt.size=0)