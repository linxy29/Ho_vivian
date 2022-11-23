## This file is used to check the similarity of different cell types 

library(Seurat)
library(readxl)
library(tidyverse)

## Read file
integrateObj = readRDS("/usersdata/share/linxy/2022-11-15_NCC.MES.EU.fastMNN.annot.RDS")
matrisome_list <- read_excel("/usersdata/share/linxy/matrisome_list.xlsx")
TF_names_v_3 <- read.table("/usersdata/share/linxy/TF_names_v_3.txt", quote="\"", comment.char="")

## Gene list
matrisomeL = matrisome_list$`Gene Symbol`
tfL = TF_names_v_3$V1
celltypeL = Idents(integrateObj) %>% unique()

## calculate similarity
integrateObj <- NormalizeData(object = integrateObj, assay = 'SCT.CC.reg.subset')
matrisome_comp = intersect(rownames(fastMNN_data), matrisomeL)
tf_comp = intersect(rownames(fastMNN_data), tfL)

### NCC vs iPSC
cell_numC = table(integrateObj@meta.data$orig.ident, Idents(integrateObj))
# function 8: calculate p-values of fisher's exact test
getFisher = function(counts1, counts2){
  ### Check counts1 and counts2 shape
  if (length(counts1) == 1 || is.null(dim(counts1)) || length(dim(counts1)) < 2) {
    counts1 = matrix(counts1, nrow=1)
  }
  if (length(counts2) == 1 || is.null(dim(counts2)) ||
      length(dim(counts2)) < 2) {
    counts2 = matrix(counts2, nrow=1)
  }
  ### set up
  totalCounts1 = colSums(counts1)
  totalCounts2 = colSums(counts2)
  count_mtrx = rbind(totalCounts1, totalCounts2)
  ### test 
  fisher_pvals = rep(NA,ncol(count_mtrx))
  for (i in 1:ncol(count_mtrx)) {
    tested_table = cbind(count_mtrx[, i], rowSums(count_mtrx)-count_mtrx[, i])
    fisher_pvals[i] = fisher.test(tested_table)$p.value
  }
  return(fisher_pvals)
}

ipsc_mes_fisher = getFisher(cell_numC["ipsc_nc_diff",], cell_numC["MES",])
ipsc_ncc_fisher = getFisher(cell_numC["ipsc_nc_diff",], cell_numC["NCC",])
fisher_res = rbind(ipsc_mes_fisher, ipsc_ncc_fisher)
colnames(fisher_res) = colnames(cell_numC)

#### one cell type
cell_type = celltypeL[1]
cell_numM = matrix(nrow = length(matrisome_comp), ncol = 2)
celltype_subset = subset(integrateObj, idents = cell_type)
ncc_celltype_subset = subset(integrateObj, ) 
fastMNN_data = integrateObj@assays$SCT.CC.reg.subset@data

matrisome_data = fastMNN_data[intersect(rownames(fastMNN_data), matrisomeL),]
tf_data = fastMNN_data[intersect(rownames(fastMNN_data), tfL),]

