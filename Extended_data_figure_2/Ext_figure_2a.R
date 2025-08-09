library(Seurat)

seurat_OBfree <- readRDS("../../seurat_OBfree.rds")

#Order of cell on top of one another 
DimPlot(seurat_OBfree, reduction = "umap", 
        order = c( 
          "Astrocytes",
          "Apical Progenitors",
          "OPCs",
          "Intermediate Progenitors",
          "Migrating PNs",
          "Differenciating PNs",
          "Interneurons",
          "Cajal Retzius Cells",
          "Endothelial cells",
          "Microglia",        
          "Pericytes" ),
        cols = c("Apical Progenitors" = "#ba8bd4",
                 "Astrocytes" = "#e8ecfb",
                 "Cajal Retzius Cells" = "deeppink",
                 "Interneurons" = "#cae0ab", 
                 "Differenciating PNs" = "#4eb265",
                 "Intermediate Progenitors" = "#1965b0",
                 "Endothelial cells" = "#dc050c",
                 "Migrating PNs" = "#7bafde",
                 "Microglia" = "#ee8026",
                 "OPCs" = "#994f88",
                 "Pericytes" = "#f7cb45")) 

#Order of cell on top of one another 
DimPlot(seurat_OBfree, reduction = "umap", 
        order = c( 
          "Astrocytes",
          "Apical Progenitors",
          "OPCs",
          "Intermediate Progenitors",
          "Migrating PNs",
          "Differenciating PNs",
          "Interneurons",
          "Cajal Retzius Cells",
          "Endothelial cells",
          "Microglia",        
          "Pericytes" ),
        cols = c("Apical Progenitors" = "#ba8bd4",
                 "Astrocytes" = "#e8ecfb",
                 "Cajal Retzius Cells" = "deeppink",
                 "Interneurons" = "#cae0ab", 
                 "Differenciating PNs" = "#4eb265",
                 "Intermediate Progenitors" = "#1965b0",
                 "Endothelial cells" = "#dc050c",
                 "Migrating PNs" = "#7bafde",
                 "Microglia" = "#ee8026",
                 "OPCs" = "#994f88",
                 "Pericytes" = "#f7cb45"), split.by = "type") 
