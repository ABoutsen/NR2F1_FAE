library(Seurat)
library(ggplot2)

seurat_OBfree <- readRDS("../../seurat_OBfree.rds")

cluster.markers <- FindMarkers(seurat_OBfree, ident.2="Ctrl", ident.1="EtOH", group.by = 'type', subset.ident = "Migrating Neurons", min.pct = 0.1, logfc.threshold = "0", test.use = "MAST", verbose = T)

GAD <- GetAssayData(object = seurat_OBfree, slot="counts") 
row.names.remove.mito <- stringr::str_subset(GAD@Dimnames[[1]], "^mt-")  
row.names.remove.hbb  <- stringr::str_subset(GAD@Dimnames[[1]], "^Hbb-") 
row.names.remove.ribo  <- stringr::str_subset(GAD@Dimnames[[1]], "^Rps|^Rpl")
row.names.remove.sex  <- c("Xist", "Tsix", "Uba1y", "Gm28587", "Kdm5d", "Eif2s3y","Gm29650", "Uty", "Ddx3y", "Gm21860", "Gm47283")
cluster.markers <- cluster.markers[!(row.names(cluster.markers) %in% c(row.names.remove.sex, row.names.remove.mito, row.names.remove.ribo)), ]

cluster.markers$threshold = factor(ifelse(cluster.markers$avg_log2FC > 0.1 & cluster.markers$p_val_adj < 0.01, 1, 
                                          ifelse(cluster.markers$avg_log2FC < -0.1 & cluster.markers$p_val_adj < 0.01, -1, 0)))
cluster.markers$p_val_adj <- as.numeric(as.character(cluster.markers$p_val_adj))										 

table(cluster.markers$threshold)

ggplot(data=cluster.markers, aes(x=avg_log2FC, y=-log10(p_val_adj))) + 
  geom_point(aes(color=threshold), alpha=0.99, size=0.5) +
  scale_colour_manual(name = "Threshold", values = c("blue", "grey", "red")) +
  geom_hline(yintercept=2, color="grey", alpha=1.0) +
  geom_vline(xintercept=c(-0.1,0.1), color="grey", alpha=1.0) +
  xlab("Average log 2 fold change") +
  ylab("-log10 (pval_adj)") +
  theme_minimal()   
