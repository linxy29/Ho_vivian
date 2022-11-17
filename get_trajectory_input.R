library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(tidyverse)


pathSPC = "/storage/holab/linxy/vivian/veloAE/"  ## for server
#pathSPC = "/mnt/d/Data"   ## for WSL
#pathSPC = "D:/Data"    ## for windows
options(future.globals.maxSize = 40000 * 1024^2)


nccmeseu = readRDS("/usersdata/share/linxy/2022-11-15_NCC.MES.EU.fastMNN.annot.RDS")
DimPlot(nccmeseu, label = T)


head(Cells(nccmeseu))
write.csv(Cells(nccmeseu), file = str_c(pathSPC, "nccmeseu_cellID_obs.csv"), row.names = FALSE)
head(Embeddings(nccmeseu, reduction = "umap"))
write.csv(Embeddings(nccmeseu, reduction = "umap"), file = str_c(pathSPC, "nccmeseu_cell_embeddings.csv"))
head(Idents(nccmeseu))
write.csv(Idents(nccmeseu), file = str_c(pathSPC, "nccmeseu_clusters.csv"))