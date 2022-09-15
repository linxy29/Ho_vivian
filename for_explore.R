library(Seurat)
library(tidyverse)

seuratObj = readRDS("/storage/holab/linxy/iPSC/seuratObj/2022-08-25_combined.NCC.MES.EU.clustered.RDS")

DefaultAssay(seuratObj) <- 'RNA'
markers <- FindAllMarkers(seuratObj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10genes = markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

DefaultAssay(seuratObj) <- 'integrated'

DoHeatmap(seuratObj, features = top10genes$gene, size = 3)

