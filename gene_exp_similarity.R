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

## DE genes
code_ident = Idents(integrateObj) %>% as.character()
code_ident = ifelse(code_ident == "Cell cycle", "cell_cycle", code_ident)
code_ident = ifelse(code_ident == "IPSC meso, MES", "ipsc_meso_mes", code_ident)
code_ident = ifelse(code_ident == "Dev", "dev", code_ident)
code_ident = ifelse(code_ident == "NLC", "nlc", code_ident)
code_ident = ifelse(code_ident == "Cardiac mesoderm", "cardiac_mesoderm", code_ident)
code_ident = ifelse(code_ident == "Endoderm", "endoderm", code_ident)
code_ident = ifelse(code_ident == "MES, dev", "mes_dev", code_ident)
code_ident %>% unique()

code_ident = str_c(code_ident, "_", integrateObj@meta.data$orig.ident)
integrateObj$ident_sample = code_ident
Idents(integrateObj) = 'ident_sample'

## ncc vs iPSC
ncc_markerL = list()
ncc_markerV = rep(NA, 3)
ncc_celltype = c("cell_cycle", "nlc", "mes_dev")
i = 1
for (celltype in ncc_celltype) {
  temp = FindMarkers(integrateObj, ident.1 = str_c(celltype, "_ipsc_nc_diff"), ident.2 = str_c(celltype, "_NCC"), min.pct = 0.25)
  ncc_markerL[[celltype]] = temp
  ncc_markerV[i] = temp %>% 
    filter(p_val_adj < 0.05) %>% 
    rownames_to_column("gene") %>% 
    .$gene %>% 
    paste(collapse=" ")
  i = i+1
}

## mes vs iPSC
mes_markerL = list()
mes_markerV = rep(NA, 3)
mes_celltype = c("ipsc_meso_mes", "dev", "nlc", "cardiac_mesoderm", "endoderm")
i = 1
for (celltype in mes_celltype) {
  temp = FindMarkers(integrateObj, ident.1 = str_c(celltype, "_ipsc_nc_diff"), ident.2 = str_c(celltype, "_MES"), min.pct = 0.25)
  mes_markerL[[celltype]] = temp
  mes_markerV[i] = temp %>% 
    filter(p_val_adj < 0.05) %>% 
    rownames_to_column("gene") %>% 
    .$gene %>% 
    paste(collapse=" ")
  i = i+1
}

save(ncc_markerL, ncc_markerV, mes_markerL, mes_markerV, file = "similarity_analysis.RData")