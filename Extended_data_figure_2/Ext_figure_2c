library(Seurat)
library(monocle3)

seurat_cds<- readRDS("F:/GIGA - PHD/250514_seurat_cds.rds")

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


DimPlot(hfile, reduction = "umap", group.by = "SubCluster_Apical", order = T, cols = c("Migrating Neurons" = "lightgrey",
                                                                                             "Projection Neurons" = "lightgrey",
                                                                                             "CGE Interneurons" = "lightgrey",
                                                                                             "MGE Interneurons" = "lightgrey",
                                                                                             "OB Interneurons" = "lightgrey",
                                                                                             "Apical Progenitors" = "lightgrey",
                                                                                             "Astrocytes" ="lightgrey",
                                                                                             "Intermediate Progenitors" = "lightgrey",
                                                                                             "OPCs" = "lightgrey",
                                                                                             "Pericytes" = "lightgrey",
                                                                                             "Endothelial cells" = "lightgrey",
                                                                                             "Microglia" = "lightgrey",
                                                                                             "Cajal Retzius Cells" = "lightgrey")) + NoLegend()
