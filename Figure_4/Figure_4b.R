library(ggplot2)

cluster.markers <- readr::read_table("F:/GIGA - PHD/01_MAST_DEG_EtOH_vs_Ctrl_Migrating Neurons_filtered.tsv")
cluster.markers_10 <- cluster.markers[1:10,]

cluster.markers_10$X8 = factor(ifelse(cluster.markers_10$threshold > 0, "Up", "Down"))
cluster.markers_10$X9 = "DEG"

mid <- 0.00

ggplot(cluster.markers_10, aes(x=X9, y=gene_name, size=-log10(p_val_adj), color=avg_log2FC)) + 
  geom_point() + 
  theme_classic() + 
  scale_color_gradient2(midpoint = 0.00, low = "red", mid = "white", high = "blue")+
  scale_size(range = c(1, 10))
