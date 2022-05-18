library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(tidyverse)


pathSPC = "/storage/holab/linxy/vivian/veloAE/"  ## for server
#pathSPC = "/mnt/d/Data"   ## for WSL
#pathSPC = "D:/Data"    ## for windows
options(future.globals.maxSize = 40000 * 1024^2)


load("/storage/holab/linxy/vivian/seuratObj/clustered_ncc.rdata")
DimPlot(ncc, reduction = "umap", label = T)


head(Cells(ncc))
write.csv(Cells(ncc), file = str_c(pathSPC, "ncc_cellID_obs.csv"), row.names = FALSE)
head(Embeddings(ncc, reduction = "umap"))
write.csv(Embeddings(ncc, reduction = "umap"), file = str_c(pathSPC, "ncc_cell_embeddings.csv"))
head(Idents(ncc))
write.csv(Idents(ncc), file = str_c(pathSPC, "ncc_clusters.csv"))