library(Seurat)
library(monocle3)

seurat_OBfree <- readRDS("../../seurat_OBfree.rds")

# Pseudotime Analysis
seurat_OBfree <- subset(seurat_OBfree, ident=c("Endothelial cells", "Microglia", "Pericytes", "MGE Interneurons", "CGE Interneurons", "Cajal Retzius Cells"), invert=T)
DimPlot(seurat_OBfree, label = T)

seurat_cds <- SeuratWrappers::as.cell_data_set(seurat_OBfree)
seurat_cds <- cluster_cells(cds = seurat_cds, reduction_method = "UMAP")
seurat_cds <- learn_graph(seurat_cds, use_partition = T, close_loop = T)

seurat_cds <- order_cells(seurat_cds, reduction_method = "UMAP")

# Figure generation
plot_cells(cds = seurat_cds,
           color_cells_by = "SubCluster_Apical", 
           show_trajectory_graph = TRUE,
           label_leaves =  F, 
           label_branch_points = F, 
           label_cell_groups = F, 
) + scale_color_manual(breaks = c("Migrating Neurons", 
                                  "Projection Neurons",
                                  "Apical Progenitors", 
                                  "Astrocytes", 
                                  "Intermediate Progenitors", 
                                  "OPCs"), 
                       values=c("#7bafde", 
                                "#4eb265",
                                "#ba8bd4", 
                                "#d9cce3",
                                "#1965b0",
                                "#994f88"))

hfile <- readRDS("F:/GIGA - PHD/h5_test_OBfree.rds")

DimPlot(hfile, reduction = "umap", group.by = "SubCluster_Apical", order = T, cols = c("Migrating PNs" = "lightgrey",
                                                                                       "Differenciating PNs" = "lightgrey",
                                                                                       "Interneurons" = "lightgrey",
                                                                                       "Apical Progenitors" = "lightgrey",
                                                                                       "Astrocytes" ="lightgrey",
                                                                                       "Intermediate Progenitors" = "lightgrey",
                                                                                       "OPCs" = "lightgrey",
                                                                                       "Pericytes" = "lightgrey",
                                                                                       "Endothelial cells" = "lightgrey",
                                                                                       "Microglia" = "lightgrey",
                                                                                       "Cajal Retzius Cells" = "lightgrey")) + NoLegend()
