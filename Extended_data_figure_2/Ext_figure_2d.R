library(Seurat)
library(ggplot2)

hfile <- readRDS("F:/GIGA - PHD/h5_test_OBfree.rds")
seurat_cds<- readRDS("F:/GIGA - PHD/250514_seurat_cds.rds")

seurat_pseudotime <- subset(hfile, ident=c("Endothelial cells", "Microglia", "Pericytes", "CGE Interneurons", "MGE Interneurons", "Cajal Retzius Cells"), invert=T)
seurat_pseudotime <- AddMetaData(object = seurat_pseudotime,
                                 metadata = seurat_cds@principal_graph_aux@listData$UMAP$pseudotime, 
                                 col.name = "pseudotime")

##### Pseudotime Ridge Plot #####
test1<-seurat_pseudotime
test1@meta.data[["pseudotime"]]<- replace(test1@meta.data[["pseudotime"]], test1@meta.data[["pseudotime"]]=='Inf', 60)

RidgePlot(test1, features = 'pseudotime') + scale_fill_manual(values = alpha(c("#7bafde", "#4eb265", "#1965b0", "#ba8bd4",  "#994f88",  "#d9cce3"), c(0.5)))
