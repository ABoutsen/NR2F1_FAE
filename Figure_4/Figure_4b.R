library(ggplot2)

cluster.markers <- readr::read_table("../Data For Analysis/01_MAST_DEG_EtOH_vs_Ctrl_Migrating Neurons_filtered.tsv")
colnames(cluster.markers) <- c("gene_name", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "threshold")

cluster.markers_10 <- cluster.markers[1:10,]

cluster.markers_10$X8 = factor(ifelse(cluster.markers_10$threshold > 0, "Up", "Down"))
cluster.markers_10$X9 = "DEG"

mid <- 0.00

ggplot(cluster.markers_10, aes(x=X9, y=gene_name, size=-log10(p_val_adj), color=avg_log2FC)) + 
  geom_point() + 
  theme_classic() + 
  scale_color_gradient2(midpoint = 0.00, low = "blue", mid = "white", high = "red")+
  scale_size(range = c(1, 10))
