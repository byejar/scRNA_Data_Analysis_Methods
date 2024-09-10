<!--
# scRNA_pipline

scRNA_pipline+fig 
scRNA+public data（读取不同类型的公开数据：h5 excel radta）
scrna+打分函数（ucell aucell fgsea/javagsea/ssgsea）
单细胞研究结果结合生存期：scRNA+TCGA survival（）
通讯、发育轨迹：scRNA+ cellchat、Monocle 3
肿瘤细胞：scRNA+CopyKAT/infercnv
空间通讯：scRNA+spRNA+cellchat——v2
代谢预测：scRNA+compass
-->
## Single-Cell RNA (scRNA) Data Analysis Methods

1. **scRNA_pipline+fig**: 
   - Standard scRNA-seq analysis pipeline with figure generation for visualization of results.

2. **scRNA+public data**: 
   - Incorporates publicly available datasets in various formats such as **h5**, **Excel**, and **rdata** for comparative analysis.

3. **scRNA+Scoring Functions**: 
   - Uses scoring functions like **UCell**, **AUCell**, and various GSEA methods (**fgsea**, **javaGSEA**, **ssGSEA**) to quantify pathway activity in single cells.

4. **scRNA+TCGA survival**: 
   - Integrates scRNA-seq results with **TCGA survival analysis** to associate cellular heterogeneity with patient outcomes.

5. **scRNA+cellchat, Monocle 3**: 
   - For analyzing **cell-cell communication** using **CellChat** and **developmental trajectory analysis** with **Monocle 3**.

6. **scRNA+CopyKAT/inferCNV**: 
   - For tumor cell analysis, applies **CopyKAT** or **inferCNV** to detect copy number variations (CNVs) in scRNA-seq data.

7. **scRNA+spRNA+cellchat v2**: 
   - Integrates **spatial transcriptomics** (ST) and **CellChat v2** to study spatial cell communication networks.

8. **scRNA+COMPASS**: 
   - Predicts metabolic activity and diversity using **COMPASS** in scRNA-seq data.


