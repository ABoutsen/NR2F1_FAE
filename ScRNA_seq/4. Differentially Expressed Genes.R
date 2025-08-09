library(Seurat)

seurat_OBfree <- readRDS("../../seurat_OBfree.rds")

#### DEG ANALYSIS ####
GAD <- GetAssayData(object = seurat_OBfree, slot="counts") #GAD@Dimnames[[1]]
row.names.remove.mito <- stringr::str_subset(GAD@Dimnames[[1]], "^mt-")  
row.names.remove.hbb  <- stringr::str_subset(GAD@Dimnames[[1]], "^Hbb-") 
row.names.remove.ribo <- stringr::str_subset(GAD@Dimnames[[1]], "^Rps|^Rpl")
row.names.remove.sex  <- c("Xist", "Tsix", "Uba1y", "Gm28587", "Kdm5d", "Eif2s3y","Gm29650", "Uty", "Ddx3y", "Gm21860", "Gm47283")

idents.DGE <- levels(Idents(seurat_OBfree))

for (f in 1:length(idents.DGE)) {
  
  group <- WhichCells(seurat_OBfree, idents = idents.DGE[f])
  
    cluster.markers <- FindMarkers(seurat_OBfree, ident.1="EtOH", ident.2="Ctrl", group.by = 'type', subset.ident = idents.DGE[f], verbose=F, test.use = "MAST", min.pct = 0.1, logfc.threshold = "0.1")
    cluster.markers <- cluster.markers[!(row.names(cluster.markers) %in% c(row.names.remove.sex, row.names.remove.mito, row.names.remove.ribo)), ]
    
    ### Volcano plot of the DEG between control and mutant ffr migrating projection neurons population ###
    cluster.markers$threshold = factor(ifelse(cluster.markers$avg_log2FC > 0.1 & cluster.markers$p_val_adj < 0.001, 1, 
                                              ifelse(cluster.markers$avg_log2FC < -0.1 & cluster.markers$p_val_adj < 0.001, -1, 0)))
    cluster.markers$p_val_adj <- as.numeric(as.character(cluster.markers$p_val_adj))										 
    
    print(table(cluster.markers$threshold))
    
}

### Finding DEGs for migrating projection neurons subpopulations ###

cluster.markers.migr <- FindMarkers(seurat_OBfree, ident.2="EtOH", ident.1="Ctrl", group.by = 'type', subset.ident = "Migrating Neurons", min.pct = 0.1, logfc.threshold = "0.1", test.use = "MAST", )
cluster.markers.migr <- cluster.markers.migr[!(row.names(cluster.markers.migr) %in% c(row.names.remove.sex, row.names.remove.mito, row.names.remove.ribo)), ]

cluster.markers.migr$threshold = factor(ifelse(cluster.markers.migr$avg_log2FC > 0.1 & cluster.markers.migr$p_val_adj < 0.001, 1, 
                                               ifelse(cluster.markers.migr$avg_log2FC < -0.1 & cluster.markers.migr$p_val_adj < 0.001, -1, 0)))
cluster.markers.migr$p_val_adj <- as.numeric(as.character(cluster.markers.migr$p_val_adj))										 

table(cluster.markers.migr$threshold)
