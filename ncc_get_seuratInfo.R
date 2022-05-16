library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(tidyverse)


pathSPC = "/storage/holab/linxy/vivian/veloAE/"  ## for server
#pathSPC = "/mnt/d/Data"   ## for WSL
#pathSPC = "D:/Data"    ## for windows
options(future.globals.maxSize = 40000 * 1024^2)


ncc_seuratObj = readRDS("/storage/holab/linxy/vivian/seuratObj/2021-11-09_NCC.rds")
DimPlot(ncc_seuratObj, reduction = "umap", label = T)


head(Cells(ncc_seuratObj))
write.csv(Cells(ncc_seuratObj), file = str_c(pathSPC, "ncc_cellID_obs.csv"), row.names = FALSE)
head(Embeddings(ncc_seuratObj, reduction = "umap"))
write.csv(Embeddings(ncc_seuratObj, reduction = "umap"), file = str_c(pathSPC, "ncc_cell_embeddings.csv"))
head(Idents(ncc_seuratObj))
write.csv(Idents(ncc_seuratObj), file = str_c(pathSPC, "ncc_clusters.csv"))