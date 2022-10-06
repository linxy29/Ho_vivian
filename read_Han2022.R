library(Seurat)

expression_mtx = read.table(gzfile("/usersdata/share/linxy/Han_2022/ctrl_matrix.tsv.gz"))
dim(expression_mtx) ## [1] 18108   251
seurat_object <- CreateSeuratObject(counts = expression_mtx)

expression_mtx = read.table(gzfile("/usersdata/share/linxy/Han_2022/NP2_matrix.tsv.gz"))
dim(expression_mtx) ## [1] 26546  4037
seurat_object <- CreateSeuratObject(counts = expression_mtx)
