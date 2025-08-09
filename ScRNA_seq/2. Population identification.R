library(Seurat)

#### CLUSTERS IDENTIFICATION ####

### Overclustering and UMAP ###
# Perform clustering with high resolution (overclustering)
combined_seurat <- FindClusters(combined_seurat, algorithm = 1, verbose = FALSE, resolution = 5)

# Run UMAP for visualization using the first 25 PCs
combined_seurat <- RunUMAP(combined_seurat, dims = 1:25)

# Visualize clusters without labels
DimPlot(combined_seurat, label = FALSE)


## 1. Apical Progenitors (VZ) ## 
genes <- c("Hes5", "Pax6", "Sox2", "Ednrb", "Hes1", "Aldoc",
           "Vim", "Ddah1", "Tspan12", "Mfge8", "Mt1", "Ndrg2",
           "Notch2", "Psat1", "Psph", "Wwtr1", "Pdpn")
FeaturePlot(combined_seurat , features = genes)

combined_seurat  <- AddModuleScore(object = combined_seurat , features = list(genes), name = c("Apcs"))

FeaturePlot(object = combined_seurat , features = c("Apcs1"), order = T)
DotPlot(object = combined_seurat , features = c("Apcs1")) + RotatedAxis()


## 2. Intermediate Progenitors (SVZ-IZ) ## 
genes <- c("Eomes", "Neurog2", "Btg2", "Celsr1", "Chd7",
                                         "Ezr", "Gadd45g", "Heg1", "Kif26b", "Mfap4",
                                         "Myo10", "Rhbdl3", "Slc16a2", "Mfng", "Coro1c")
FeaturePlot(combined_seurat , features = genes)

combined_seurat  <- AddModuleScore(object = combined_seurat , features = list(genes), name = c("Ipcs"))

FeaturePlot(object = combined_seurat , features = c("Ipcs1"), order = T)
DotPlot(object = combined_seurat , features = c("Ipcs1")) 


## 3. Cajal retzius cells ## 
genes <- c("Lhx1", "Lhx5", "Reln", "Trp73")
FeaturePlot(combined_seurat , features = genes)

combined_seurat  <- AddModuleScore(object = combined_seurat , features = list(genes), name = "CRC")

FeaturePlot(object = combined_seurat , features = c("CRC1"), order = T)
DotPlot(object = combined_seurat , features = c("CRC1")) 


## 4. Interneurons ## 
genes <- c("Dlx1", "Dlx2", "Dlx5", "Dlx6", "Gad1", "Gad2")
FeaturePlot(combined_seurat , features = genes)

combined_seurat  <- AddModuleScore(object = combined_seurat , features = list(genes), name = c("IN"))

FeaturePlot(object = combined_seurat , features = c("IN1"), order = T)
DotPlot(object = combined_seurat , features = c("IN1")) 


## 5. OPCS ## 
genes <- c("Olig1", "Olig2", "Pdgfra", "Sox10")
FeaturePlot(combined_seurat , features = genes)

combined_seurat  <- AddModuleScore(object = combined_seurat , features = list(genes), name = c("OPCs"))

FeaturePlot(object = combined_seurat , features = c("OPCs1"), order = T)
DotPlot(object = combined_seurat , features = c("OPCs1")) 


## 6. Astrocytes ## 
genes <- c("Aldh1l1", "Apoe", "Aqp4", "Gfap", "Slc1a3", "Sparcl1")
FeaturePlot(combined_seurat , features = genes)

combined_seurat  <- AddModuleScore(object = combined_seurat , features = list(genes), name = c("Astrocytes"))

FeaturePlot(object = combined_seurat , features = c("Astrocytes1"), order = T)
DotPlot(object = combined_seurat , features = c("Astrocytes1")) 


## 7. Endothelial cells ## 
genes <- c("Adgrf5", "Adgrl4", "Cldn5", "Igfbp7", "Mcam")
FeaturePlot(combined_seurat , features = genes)

combined_seurat  <- AddModuleScore(object = combined_seurat , features = list(genes), name = c("EC"))

FeaturePlot(object = combined_seurat , features = c("EC1"), order = T)
DotPlot(object = combined_seurat , features = c("EC1")) 


## 8. Pericytes ## 
genes <- c("Cspg4", "Pdgfrb", "Kcnj8", "Rgs5", "Vtn")
FeaturePlot(combined_seurat , features = genes)

combined_seurat  <- AddModuleScore(object = combined_seurat , features = list(genes), name = c("Pericytes"))

FeaturePlot(object = combined_seurat , features = c("Pericytes1"), order = T)
DotPlot(object = combined_seurat , features = c("Pericytes1")) 


## 9. Microglia ## 
genes <- c("Aif1", "C1qb", "C1qc", "Csf1r", "Fcer1g", 
                                         "Fcrls", "Hexb", "Itgam", "Lgals9", "Mrc1",
                                         "P2ry12", "Ptprc", "Siglech", "Tmem119", "Trem2")  
FeaturePlot(combined_seurat , features = genes)


combined_seurat  <- AddModuleScore(object = combined_seurat , features = list(genes), name = c("Microglia"))

FeaturePlot(object = combined_seurat , features = c("Microglia1"), order = T)
DotPlot(object = combined_seurat , features = c("Microglia1")) 


## 10. Migrating newborns projection neurons ## 
genes <- c("Neurod1", "Pcp4", "Sema3c", "Unc5d") 
FeaturePlot(combined_seurat , features = genes)

combined_seurat  <- AddModuleScore(object = combined_seurat , features = list(genes), name = "Migrating_newborns_projection_neurons")

# Plot multiple module scores together
FeaturePlot(object = combined_seurat , features = c("Migrating_newborns_projection_neurons1"), order = T)
DotPlot(object = combined_seurat , features = c("Migrating_newborns_projection_neurons1")) 


## 11. Cortical projection neurons - Excitatory neurons  ##
genes <- c("Neurod2", "Neurod6", "Tubb3") 
FeaturePlot(combined_seurat , features = genes)

combined_seurat  <- AddModuleScore(object = combined_seurat , features = list(genes), name = c("Projection"))

FeaturePlot(object = combined_seurat , features = c("Projection1"), order = T)
DotPlot(object = combined_seurat , features = c("Projection1"))  


# Plot multiple module scores together
FeaturePlot(object = combined_seurat , features = c( "Astrocytes1", "Apcs1", "Ipcs1", "Migrating_newborns_projection_neurons1", "Projection1", "IN1", "CRC1", "OPCs1","EC1", "Pericytes1", "Microglia1")) 
DotPlot(object = combined_seurat , features = c( "Astrocytes1", "Apcs1", "Ipcs1", "Migrating_newborns_projection_neurons1", "Projection1", "IN1", "CRC1", "OPCs1","EC1", "Pericytes1", "Microglia1"))  + RotatedAxis()



##### Cluster Identification #####
# Select cells by cluster IDs for each cell type and Plot highlighted cells for each population

## 1. Apical Progenitors ##
cells_of_interest <- WhichCells(combined_seurat , idents = c(21, 32, 42, 45, 48, 51, 53, 56)) 
DimPlot(combined_seurat , label=F, cells.highlight=list(Aps), cols.highlight = c("darkblue"), cols="grey") + NoLegend()

## 2. Intermediate Progenitors ## 
cells_of_interest <- WhichCells(combined_seurat , idents = c(13, 27, 43, 46, 67, 41, 39))
DimPlot(combined_seurat , label=F, cells.highlight=list(cells_of_interest), cols.highlight = c("darkblue"), cols="grey") + NoLegend()

## 3. Cajal Retzius Cells ## 
cells_of_interest <- WhichCells(combined_seurat , idents = c(74)) 
DimPlot(combined_seurat , label=F, cells.highlight=list(cells_of_interest), cols.highlight = c("darkblue"), cols="grey") + NoLegend()

## 4. Interneurons ## 
cells_of_interest <- WhichCells(combined_seurat , idents = c(8, 11, 16, 26, 35, 66, 68, 70)) 
DimPlot(combined_seurat , label=F, cells.highlight=list(cells_of_interest), cols.highlight = c("darkblue"), cols="grey") + NoLegend()

## 5. OPCs ## 
cells_of_interest <- WhichCells(combined_seurat , idents = c(64,55, 65)) 
DimPlot(combined_seurat , label=F, cells.highlight=list(cells_of_interest), cols.highlight = c("darkblue"), cols="grey") + NoLegend()

## 6. Astrocytes ##
cells_of_interest <- WhichCells(combined_seurat , idents = c(24)) 
DimPlot(combined_seurat , label=F, cells.highlight=list(cells_of_interest), cols.highlight = c("darkblue"), cols="grey") + NoLegend()

## 7. Endothelial cells ## 
cells_of_interest <- WhichCells(combined_seurat , idents = c(69)) 
DimPlot(combined_seurat , label=F, cells.highlight=list(cells_of_interest), cols.highlight = c("darkblue"), cols="grey") + NoLegend()

## 8. Pericytes ## 
cells_of_interest <- WhichCells(combined_seurat , idents = c(62)) 
DimPlot(combined_seurat , label=F, cells.highlight=list(cells_of_interest), cols.highlight = c("darkblue"), cols="grey") + NoLegend()

## 9. Microglia ## 
cells_of_interest <- WhichCells(combined_seurat , idents = c(71, 73)) 
DimPlot(combined_seurat , label=F, cells.highlight=list(cells_of_interest), cols.highlight = c("darkblue"), cols="grey") + NoLegend()

## 10. Migrating Neurons ## 
cells_of_interest <- WhichCells(combined_seurat , idents = c(0, 1, 17, 20, 22, 30, 34, 47, 50, 57, 61, 13))
DimPlot(combined_seurat , label=F, cells.highlight=list(cells_of_interest), cols.highlight = c("darkblue"), cols="grey") + NoLegend()

## 11. Cortical Projection Neurons ## 
cells_of_interest <- WhichCells(combined_seurat , idents = c(2, 3, 4, 5, 6, 7, 9, 10, 12, 18, 19, 23, 25, 28, 29, 31, 33, 36, 37, 38, 40, 44, 49, 52, 54, 58, 59, 60, 63, 72, 75))
DimPlot(combined_seurat , label=F, cells.highlight=list(cells_of_interest), cols.highlight = c("darkblue"), cols="grey") + NoLegend()


##### Rename clusters with meaningful cell type labels #####
New.cluster.ids <- c("Migrating Neurons", #0
                     "Migrating Neurons", #1
                     "Projection Neurons", #2 
                     "Projection Neurons", #3 
                     "Projection Neurons", #4 
                     "Projection Neurons", #5  
                     "Projection Neurons", #6 
                     "Projection Neurons", #7 
                     "INs", #8
                     "Projection Neurons", #9 
                     "Projection Neurons", #10 
                     "INs", #11
                     "Projection Neurons", #12 
                     "Migrating Neurons", #13 
                     "Migrating Neurons", #14  
                     "Migrating Neurons", #15  
                     "INs", #16
                     "Migrating Neurons", #17
                     "Projection Neurons", #18 
                     "Projection Neurons", #19 
                     "Migrating Neurons", #20  
                     "Apical Progenitors", #21
                     "Migrating Neurons", #22
                     "Projection Neurons", #23 
                     "Astrocytes", #24
                     "Projection Neurons", #25 
                     "INs", #26
                     "Intermediate Progenitors", #27
                     "Projection Neurons", #28 
                     "Projection Neurons", #29 
                     "Migrating Neurons", #30
                     "Projection Neurons", #31 
                     "Apical Progenitors", #32
                     "Projection Neurons", #33 
                     "Migrating Neurons", #34
                     "INs", #35
                     "Projection Neurons", #36 
                     "Projection Neurons", #37 
                     "Projection Neurons", #38 
                     "Intermediate Progenitors", #39
                     "Projection Neurons", #40 
                     "Intermediate Progenitors", #41
                     "Apical Progenitors", #42
                     "Intermediate Progenitors", #43
                     "Projection Neurons", #44 
                     "Apical Progenitors", #45
                     "Intermediate Progenitors", #46 
                     "Migrating Neurons", #47
                     "Apical Progenitors", #48
                     "Projection Neurons", #49 
                     "Migrating Neurons", #50
                     "Apical Progenitors", #51
                     "Projection Neurons", #52 
                     "Apical Progenitors", #53
                     "Projection Neurons", #54 
                     "OPCs", #55
                     "Apical Progenitors", #56
                     "Migrating Neurons", #57
                     "Projection Neurons", #58 
                     "Projection Neurons", #59 
                     "Projection Neurons", #60 
                     "Migrating Neurons", #61
                     "Pericytes", #62
                     "Projection Neurons", #63 
                     "OPCs", #64
                     "OPCs", #Apical Progenitors", #65
                     "INs", #66
                     "Intermediate Progenitors", #67 
                     "INs", #68
                     "Endothelial cells", #69
                     "INs", #70
                     "Microglia", #71
                     "Projection Neurons", #72 
                     "Microglia", #73
                     "Cajal Retzius Cells", #74
                     "Projection Neurons") #75 

names(New.cluster.ids) <- levels(combined_seurat )
combined_seurat  <- RenameIdents(combined_seurat , New.cluster.ids)

# Plot clusters with new labels
DimPlot(combined_seurat, label = TRUE)

# Final DotPlot for all key cell types
DotPlot(combined_seurat, features = c("Migrating_newborns_projection_neurons1", "Projection1", "IN1", "Apcs1", "Astrocytes1", "Ipcs1", "OPCs1", "Pericytes1", "EC1", "Microglia1", "CRC1")) + RotatedAxis()



#### SUBSETTING INTERNEURONS ####
## filter them out:
sub_population <- subset(combined_seurat, ident="INs")
DimPlot(sub_population)

## Processing the interneurons populations ## 
sub_population <- SCTransform(sub_population, vars.to.regress = c("percent.mt", 'pct_chrX', 'pct_chrY', 'percent.hbb'), verbose = TRUE)
sub_population <- RunPCA(sub_population, features = VariableFeatures(object = sub_population), verbose = F)
sub_population <- FindNeighbors(sub_population, dims = 1:25)
sub_population <- FindClusters(sub_population, algorithm = 1, verbose = FALSE, resolution = 0.8)
sub_population <- RunUMAP(sub_population, dims = 1:25)

DimPlot(sub_population, label=T) 

## Interneurons Markers ##
IN_signature_gene_list <- list(c("Dlx1", "Dlx2", "Dlx5", "Dlx6", "Gad1", "Gad2"))

CGE_PCOA_Interneurons_signature_gene_list <- list(c("Cxcl14", "Htr3a", "Prox1", "Sp8"))
FeaturePlot(combined_seurat, features =  c("Cxcl14", "Htr3a", "Prox1")) 

MGE_Interneurons_signature_gene_list <- list(c("Lhx6", "Npy", "Nxph1", "Sst", "Nxph2"))
FeaturePlot(combined_seurat, features =  c("Lhx6", "Npy", "Nxph1", "Sst", "Sox6")) 

OB_Interneurons_signature_gene_list <- list(c("Meis2", "Sp8", "Etv1", "Sp9", "Pbx1", "Pbx3", "Rbfox3", "Tshz1", "Prokr2", "Vax1", "Gsx2", "Nr2e1", "Ascl1")) 
FeaturePlot(combined_seurat, features =  c("Meis2", "Sp8", "Etv1", "Sp9", "Pbx1", "Pbx3", "Rbfox3", "Tshz1", "Prokr2", "Vax1", "Gsx2", "Nr2e1", "Ascl1")) 

sub_population <- AddModuleScore(object = sub_population, features = c(IN_signature_gene_list), name = c("Interneurons"))
sub_population <- AddModuleScore(object = sub_population, features = c(CGE_PCOA_Interneurons_signature_gene_list), name = c("CGE_POA_Interneurons"))
sub_population <- AddModuleScore(object = sub_population, features = c(MGE_Interneurons_signature_gene_list), name = c("MGE_Interneurons"))
sub_population <- AddModuleScore(object = sub_population, features = c(OB_Interneurons_signature_gene_list), name = c("OB_Interneurons"))

FeaturePlot(object = sub_population, features = c("CGE_POA_Interneurons1", "MGE_Interneurons1", "OB_Interneurons1", "Interneurons1"), order = T)
DotPlot(object = sub_population, features = c("CGE_POA_Interneurons1", "MGE_Interneurons1", "OB_Interneurons1")) + RotatedAxis() 

## Reassignation of the clusters to a cell types ##
new.cluster.ids <- c("CGE interneurons", 
                        "MGE interneurons", 
                        "MGE interneurons", 
                        "MGE interneurons",  
                        "OB interneurons",  
                        "OB interneurons", 
                        "MGE interneurons",  
                        "OB interneurons", 
                        "MGE interneurons", 
                        "OB interneurons", 
                        "MGE interneurons")

names(new.cluster.ids) <- levels(sub_population)
sub_population <- RenameIdents(sub_population, new.cluster.ids)

DimPlot(sub_population, label=F)

MGE <- WhichCells(sub_population, idents = c("MGE interneurons")) 
CGE <- WhichCells(sub_population, idents = c("CGE interneurons")) 
OB <- WhichCells(sub_population, idents = c("OB interneurons")) 




#### SUBSETTING OPCS ####
## filter them out:
sub_population <- subset(combined_seurat, ident="OPCs")
DimPlot(sub_population)

## Processing the OPCs populations ## 
sub_population <- SCTransform(sub_population, vars.to.regress = c("percent.mt", 'pct_chrX', 'pct_chrY', 'percent.hbb'), verbose = TRUE)
sub_population <- RunPCA(sub_population, features = VariableFeatures(object = sub_population), verbose = F)
sub_population <- FindNeighbors(sub_population, dims = 1:25)
sub_population <- FindClusters(sub_population, algorithm = 1, verbose = FALSE, resolution = 0.8)
sub_population <- RunUMAP(sub_population, dims = 1:25)

DimPlot(sub_population, label=T) 

## OPCs Markers ##
OPCS_signature_gene_list <- list(c("Olig1", "Olig2", "Pdgfra", "Sox10"))
FeaturePlot(sub_population, features = c("Olig1", "Olig2", "Pdgfra", "Sox10")) 

sub_population <- AddModuleScore(object = sub_population, features = c(OPCS_signature_gene_list), name = c("OPCs"))

## Apical Progenitors MArkers ##
Apcs_signature_gene_list <- list(c("Hes5", "Pax6", "Sox2", "Ednrb", "Hes1", "Aldoc", "Pdpn", "Vim", "Ddah1", "Tspan12", "Mfge8", "Mt1", "Ndrg2", "Notch2", "Psat1", "Psph", "Wwtr1"))
FeaturePlot(sub_population, features=c("Hes5", "Pax6", "Sox2", "Ednrb", "Hes1", "Aldoc", "Pdpn", "Vim", "Ddah1", "Tspan12", "Mfge8", "Mt1", "Ndrg2", "Notch2", "Psat1", "Psph", "Wwtr1")) 

sub_population <- AddModuleScore(object = sub_population, features = c(Apcs_signature_gene_list), name = c("Apcs"))


FeaturePlot(object = sub_population, features = c("OPCs1", "Apcs1"), order = T)
DotPlot(object = sub_population, features = c("OPCs1", "Apcs1")) + RotatedAxis() 

new.cluster.ids <- c("Apical Prog", #OPCs", #0
                          "OPCs", #1
                          "OPCs", #2
                          "Apical Prog", #3
                          "OPCs", #Apical Prog", #4  
                          "OPCs", #5
                          "OPCs", #6
                          "Apical Prog", #7  
                          "OPCs", #8
                          "OPCs") #13

names(new.cluster.ids) <- levels(sub_population)
sub_population_1 <- RenameIdents(sub_population, new.cluster.ids)

DimPlot(sub_population_1, label=F)# + NoLegend()  

OP <- WhichCells(sub_population_1, idents = c("OPCs")) 
AP <- WhichCells(sub_population_1, idents = c("Apical Prog")) 




#### SUBSETTING PROJECTION NEURONS ####
## filter them out:
sub_population <- subset(combined_seurat, ident="Projection Neurons")
DimPlot(sub_population)

## Processing the Projection Neurons populations ## 
sub_population <- SCTransform(sub_population, vars.to.regress = c("percent.mt", 'pct_chrX', 'pct_chrY', 'percent.hbb'), verbose = TRUE)
sub_population <- RunPCA(sub_population, features = VariableFeatures(object = sub_population), verbose = F)
sub_population <- FindNeighbors(sub_population, dims = 1:25)
sub_population <- FindClusters(sub_population, algorithm = 1, verbose = FALSE, resolution = 0.8)
sub_population <- RunUMAP(sub_population, dims = 1:25)

DimPlot(sub_population, label=T) 

## Projection Neurons Markers ##
L23_signature_gene_list <- list(c("Cux1", "Cux2", "Satb2"))
L4_signature_gene_list <- list(c("Satb2", "Tbr1"))
L5_signature_gene_list <- list(c("Bcl11b", "Sox5", "Fezf2", "Satb2"))
L6_signature_gene_list <- list(c("Bcl11b", "Sox5", "Tbr1"))

sub_population <- AddModuleScore(object = sub_population, features = c(L23_signature_gene_list), name = c("L23"))
sub_population <- AddModuleScore(object = sub_population, features = c(L4_signature_gene_list), name = c("L4"))
sub_population <- AddModuleScore(object = sub_population, features = c(L5_signature_gene_list), name = c("L5"))
sub_population <- AddModuleScore(object = sub_population, features = c(L6_signature_gene_list), name = c("L6"))
 
FeaturePlot(object = sub_population, features = c("L231", "L41", "L51", "L61"), order = T)
DotPlot(object = sub_population, features = c("L231", "L41", "L51", "L61")) 

new.cluster.ids <- c("L5&6", #0
                       "L5&6", #1
                       "L2&3", #2
                       "L2&3", #3  
                       "L2&3", #4
                       "L2&3", #5
                       "L5&6", #6
                       "L4", #7
                       "L5&6", #8
                       "L5&6", #9
                       "L5&6", #10
                       "L2&3", #11
                       "L5&6", #11
                       "L5&6") #12

names(new.cluster.ids) <- levels(sub_population)
sub_population_2 <- RenameIdents(sub_population, new.cluster.ids)

DimPlot(sub_population_2, label=F)

L23 <- WhichCells(sub_population_2, idents = c("L2&3")) 
L4  <- WhichCells(sub_population_2, idents = c("L4")) 
L56 <- WhichCells(sub_population_2, idents = c("L5&6"))




#### SUBSETTING Astrocytes VS APICAL PROGENITORS ####
sub_population <- subset(combined_seurat, ident= c('Apical Progenitors', 'Astrocytes'))
DimPlot(sub_population)

## Processing the Projection Neurons populations ## 
sub_population <- SCTransform(sub_population, vars.to.regress = c("percent.mt", 'pct_chrX', 'pct_chrY', 'percent.hbb'), verbose = TRUE)
sub_population <- RunPCA(sub_population, features = VariableFeatures(object = sub_population), verbose = F)
sub_population <- FindNeighbors(sub_population, dims = 1:25)
sub_population <- FindClusters(sub_population, algorithm = 1, verbose = FALSE, resolution = 0.8)
sub_population <- RunUMAP(sub_population, dims = 1:25)

DimPlot(sub_population, label=T) 

## Apical progenitors and astrocytes markers Markers ##
Astrocytes_signature_gene_list <- list(c("Aldh1l1", "Apoe", "Aqp4", "Gfap", "Slc1a3", "Sparcl1"))
Apcs_signature_gene_list <- list(c("Hes5", "Pax6", "Sox2", "Ednrb", "Hes1", "Aldoc", "Vim", "Ddah1", "Tspan12", "Mfge8", "Mt1", "Ndrg2", "Notch2", "Psat1", "Psph", "Wwtr1", "Pdpn"))

sub_population <- AddModuleScore(object = sub_population, features = c(Astrocytes_signature_gene_list), name = c("Astro"))
sub_population <- AddModuleScore(object = sub_population, features = c(Apcs_signature_gene_list), name = c("Apcs"))

FeaturePlot(object = sub_population, features = c("Astro1", "Apcs1"), order = T)
DotPlot(object = sub_population, features = c("Astro1", "Apcs1"))  

new.cluster.ids <- c("Apical Progenitors", #0
                        "Apical Progenitors", #1
                        "Apical Progenitors", #2
                        "Apical Progenitors", #3
                        "Apical Progenitors", #4  
                        "Apical Progenitors", #5
                        "Apical Progenitors", #6
                        "Apical Progenitors", #7
                        "Apical Progenitors", #8
                        "Apical Progenitors", #9
                        "Astrocytes", #10 
                        "Apical Progenitors", #11
                        "Astrocytes", #12
                        "Astrocytes", #13 , 
                        "Apical Progenitors") #14

names(new.cluster.ids) <- levels(sub_population)
sub_population_1 <- RenameIdents(sub_population, new.cluster.ids)

DimPlot(sub_population_1, label=F) 

Atsro_cells <- WhichCells(sub_population_1, idents = c("Astrocytes")) 
Apical_cells <- WhichCells(sub_population_1, idents = c("Apical Progenitors")) 



### Generate Population Meta_Data ###

## SubCluster and removal of OB interneurons ##

Migrating.Neurons <- WhichCells(combined_seurat, idents = c("Migrating Neurons")) 
Projection.Neurons <- WhichCells(combined_seurat, idents = c("Projection Neurons")) 
Interneurons <- WhichCells(combined_seurat, idents = c("INs")) 
Apical.Progenitors <- WhichCells(combined_seurat, idents = c("Apical Progenitors")) 
Astrocytes.cells <- WhichCells(combined_seurat, idents = c("Astrocytes")) 
Intermediate.Progenitors <- WhichCells(combined_seurat, idents = c("Intermediate Progenitors")) 
Oligo.Progenitors <- WhichCells(combined_seurat, idents = c("OPCs")) 
Pericytes.cells <- WhichCells(combined_seurat, idents = c("Pericytes")) 
Endothelial.cells <- WhichCells(combined_seurat, idents = c("Endothelial cells")) 
Microglia.cells <- WhichCells(combined_seurat, idents = c("Microglia")) 
Cajal.Retzius.cells <- WhichCells(combined_seurat, idents = c("Cajal Retzius Cells")) 
Apical.Progenitors <- WhichCells(seurat_OBfree, idents = c("Apical Progenitors"))
Astrocytes.cells <- WhichCells(seurat_OBfree, idents = c("Astrocytes"))

# Generation of metadata containing all subclusters
combined_seurat$Subcluster = 'NA'
combined_seurat[["Subcluster"]][Atsro_cells,] <- "Astrocytes"
combined_seurat[["Subcluster"]][Apical_cells,] <- "Apical Progenitors"
combined_seurat[["Subcluster"]][Oligo.Progenitors,] <- "OPCs"
combined_seurat[["Subcluster"]][MGE,] <- "MGE Interneurons"
combined_seurat[["Subcluster"]][CGE,] <- "CGE Interneurons"
combined_seurat[["Subcluster"]][OB,] <- "OB Interneurons"
combined_seurat[["Subcluster"]][L23,] <- "Projection Neurons L2&3"
combined_seurat[["Subcluster"]][L56,] <- "Projection Neurons L5&6"
combined_seurat[["Subcluster"]][L4,] <- "Projection Neurons L4"
combined_seurat[["Subcluster"]][Migrating.Neurons,] <- "Migrating Neurons"
combined_seurat[["Subcluster"]][Intermediate.Progenitors,] <- "Intermediate Progenitors"
combined_seurat[["Subcluster"]][Pericytes.cells,] <- "Pericytes"
combined_seurat[["Subcluster"]][Endothelial.cells,] <- "Endothelial cells"
combined_seurat[["Subcluster"]][Microglia.cells,] <- "Microglia"
combined_seurat[["Subcluster"]][Cajal.Retzius.cells,] <- "Cajal Retzius Cells"

DimPlot(combined_seurat, group.by = "Subcluster")

# Generation of metadata containing all subclusters without the PN layers and the Gangionic eminence interneurons combined 
combined_seurat$SubCluster_Apical = 'NA'
combined_seurat[["SubCluster_Apical"]][Atsro_cells,] <- "Astrocytes"
combined_seurat[["SubCluster_Apical"]][c(AP, Apical_cells),] <- "Apical Progenitors"
combined_seurat[["SubCluster_Apical"]][OP,] <- "OPCs"
combined_seurat[["SubCluster_Apical"]][Projection.Neurons,] <- "Differenciating PNs"
combined_seurat[["SubCluster_Apical"]][Migrating.Neurons,] <- "Migrating PNs"
combined_seurat[["SubCluster_Apical"]][Intermediate.Progenitors,] <- "Intermediate Progenitors"
combined_seurat[["SubCluster_Apical"]][Pericytes.cells,] <- "Pericytes"
combined_seurat[["SubCluster_Apical"]][Endothelial.cells,] <- "Endothelial cells"
combined_seurat[["SubCluster_Apical"]][Microglia.cells,] <- "Microglia"
combined_seurat[["SubCluster_Apical"]][Cajal.Retzius.cells,] <- "Cajal Retzius Cells"
combined_seurat[["SubCluster_Apical"]][MGE,] <- "Interneurons"
combined_seurat[["SubCluster_Apical"]][CGE,] <- "Interneurons"
combined_seurat[["SubCluster_Apical"]][OB, ] <- "OB Interneurons"

DimPlot(combined_seurat, group.by = "SubCluster_Apical")

# Generation of the dataset, removing th OB interneurons population
seurat_OBfree <- subset(combined_seurat, subset = SubCluster_Apical == "OB Interneurons", invert=T)
# saveRDS("../../seurat_OBfree.rds")
