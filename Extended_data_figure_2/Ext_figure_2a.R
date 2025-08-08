library(Seurat)

hfile <- readRDS("F:/GIGA - PHD/h5_test_OBfree.rds")
hfile <- SetIdent(hfile, value = "SubCluster_Apical")

# Ajout type All SubCluster
hfile$Population_1 = 'NA'
hfile[["Population_1"]][WhichCells(hfile, idents = c("Astrocytes")),] <- "Astrocytes"
hfile[["Population_1"]][WhichCells(hfile, idents = c("Apical Progenitors")),] <- "Apical progenitors"
hfile[["Population_1"]][WhichCells(hfile, idents = c("OPCs")),] <- "OPCs"
hfile[["Population_1"]][WhichCells(hfile, idents = c("Migrating Neurons")),] <- "Migrating PNs"
hfile[["Population_1"]][WhichCells(hfile, idents = c("Projection Neurons")),] <- "Differentiating PNs"
hfile[["Population_1"]][WhichCells(hfile, idents = c("MGE Interneurons", "CGE Interneurons")),] <- "Interneurons"
hfile[["Population_1"]][WhichCells(hfile, idents = c("Intermediate Progenitors")),] <- "Intermediate Progenitors"
hfile[["Population_1"]][WhichCells(hfile, idents = c("Pericytes")),] <- "Pericytes"
hfile[["Population_1"]][WhichCells(hfile, idents = c("Microglia")),] <- "Microglia"
hfile[["Population_1"]][WhichCells(hfile, idents = c("Cajal Retzius Cells")),] <- "Cajal Retzius Cells"
hfile[["Population_1"]][WhichCells(hfile, idents = c("Endothelial cells")),] <- "Endothelial cells"

# DimPlot(hfile, group.by = "Population_1")


##test colo3
#Order of cell on top of one another 
DimPlot(hfile, group.by = 'Population_1', reduction = "umap", 
        order = c( 
          "Astrocytes",
          "Apical progenitors",
          "OPCs",
          "Intermediate Progenitors",
          "Migrating PNs",
          "Differentiating PNs",
          "Interneurons",
          "Cajal Retzius Cells",
          "Endothelial cells",
          "Microglia",        
          "Pericytes" ),
        cols = c("Apical progenitors" = "#ba8bd4",
                 "Astrocytes" = "#e8ecfb",
                 "Cajal Retzius Cells" = "deeppink",
                 "Interneurons" = "#cae0ab", 
                 "Differentiating PNs" = "#4eb265",
                 "Intermediate Progenitors" = "#1965b0",
                 "Endothelial cells" = "#dc050c",
                 "Migrating PNs" = "#7bafde",
                 "Microglia" = "#ee8026",
                 "OPCs" = "#994f88",
                 "Pericytes" = "#f7cb45")) 
##test colo3
#Order of cell on top of one another 
DimPlot(hfile, group.by = 'Population_1', reduction = "umap", 
        order = c( 
          "Astrocytes",
          "Apical progenitors",
          "OPCs",
          "Intermediate Progenitors",
          "Migrating PNs",
          "Differentiating PNs",
          "Interneurons",
          "Cajal Retzius Cells",
          "Endothelial cells",
          "Microglia",        
          "Pericytes" ),
        cols = c("Apical progenitors" = "#ba8bd4",
                 "Astrocytes" = "#e8ecfb",
                 "Cajal Retzius Cells" = "deeppink",
                 "Interneurons" = "#cae0ab", 
                 "Differentiating PNs" = "#4eb265",
                 "Intermediate Progenitors" = "#1965b0",
                 "Endothelial cells" = "#dc050c",
                 "Migrating PNs" = "#7bafde",
                 "Microglia" = "#ee8026",
                 "OPCs" = "#994f88",
                 "Pericytes" = "#f7cb45"), split.by = "type") 
