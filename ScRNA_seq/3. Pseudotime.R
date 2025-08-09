library(Seurat)
library(monocle3)

seurat_OBfree <- readRDS("../../seurat_OBfree.rds")
seurat_OBfree <- subset(seurat_OBfree, ident=c("Endothelial cells", "Microglia", "Pericytes", "MGE Interneurons", "CGE Interneurons", "Cajal Retzius Cells"), invert=T)
DimPlot(seurat_OBfree, label = T)

seurat_cds <- SeuratWrappers::as.cell_data_set(seurat_OBfree)
seurat_cds <- cluster_cells(cds = seurat_cds, reduction_method = "UMAP")
seurat_cds <- learn_graph(seurat_cds, use_partition = T, close_loop = T)

seurat_cds <- order_cells(seurat_cds, reduction_method = "UMAP")

plot_cells(cds = seurat_cds,
           color_cells_by = "pseudotime",  
           show_trajectory_graph = TRUE,
           label_leaves =  F, 
           label_branch_points = F, 
           label_cell_groups = F, 
) 
