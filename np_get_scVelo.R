library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(tidyverse)



pathSPC = "/storage/holab/linxy"  ## for server
#pathSPC = "/mnt/d/Data"   ## for WSL
#pathSPC = "D:/Data"    ## for windows
options(future.globals.maxSize = 40000 * 1024^2)


nccSeurat = readRDS("/usersdata/share/Vivian/2021-11-09_NCC.rds")
DimPlot(nccSeurat, reduction = "umap", label = T)



head(Cells(nccSeurat))
write.csv(Cells(nccSeurat), file = str_c(pathSPC, "/vivian/ncc_cellID_obs.csv"), row.names = FALSE)
head(Embeddings(nccSeurat, reduction = "umap"))
write.csv(Embeddings(nccSeurat, reduction = "umap"), file = str_c(pathSPC, "/GCTB/trajectory/GCTB6_cell_embeddings.csv"))
head(Idents(GCTB6.combined.V2))
write.csv(Idents(GCTB6.combined.V2), file = str_c(pathSPC, "/GCTB/trajectory/GCTB6_clusters.csv"))