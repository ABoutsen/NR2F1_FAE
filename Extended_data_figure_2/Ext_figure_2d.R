library(Seurat)
library(ggplot2)

seurat_OBfree <- readRDS("../../seurat_OBfree.rds")

# Pseudotime Analysis
seurat_OBfree <- subset(seurat_OBfree, ident=c("Endothelial cells", "Microglia", "Pericytes", "MGE Interneurons", "CGE Interneurons", "Cajal Retzius Cells"), invert=T)
DimPlot(seurat_OBfree, label = T)
seurat_cds <- SeuratWrappers::as.cell_data_set(seurat_OBfree)
seurat_cds <- cluster_cells(cds = seurat_cds, reduction_method = "UMAP")
seurat_cds <- learn_graph(seurat_cds, use_partition = T, close_loop = T)

seurat_pseudotime <- subset(seurat_OBfree, ident=c("Endothelial cells", "Microglia", "Pericytes", "CGE Interneurons", "MGE Interneurons", "Cajal Retzius Cells"), invert=T)
seurat_pseudotime <- AddMetaData(object = seurat_pseudotime,
                                 metadata = seurat_cds@principal_graph_aux@listData$UMAP$pseudotime, 
                                 col.name = "pseudotime")

##### Pseudotime Ridge Plot #####
test1<-seurat_pseudotime
test1@meta.data[["pseudotime"]]<- replace(test1@meta.data[["pseudotime"]], test1@meta.data[["pseudotime"]]=='Inf', 60)

RidgePlot(test1, features = 'pseudotime') + scale_fill_manual(values = alpha(c("#7bafde", "#4eb265", "#1965b0", "#ba8bd4",  "#994f88",  "#d9cce3"), c(0.5)))
