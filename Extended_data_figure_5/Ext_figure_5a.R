library(Seurat)
library(ggplot2)

seurat_OBfree <- readRDS("../../seurat_OBfree.rds")

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

colors_PN <- ggplotColours(n=3)

num_control <- table(seurat_OBfree@meta.data[which(seurat_OBfree@meta.data$type=="Ctrl"), "Subcluster"])
num_etoh <- table(seurat_OBfree@meta.data[which(seurat_OBfree@meta.data$type=="EtOH"), "Subcluster"])
num_total <- table(seurat_OBfree@meta.data["Subcluster"] )

numcells <- as.data.frame(cbind(Control = num_control,EtOH = num_etoh, Prop_Control = num_control/sum(num_control)*100, Prop_EtOH = num_etoh/sum(num_etoh)*100))
numcells <- numcells[order(-numcells$Control),]

numcells$Cluster <- rownames(numcells)
numcells2 <- numcells
numcells <- numcells[, c(5,3,4)]
colnames(numcells) <- c("Cluster", "Control", "EtOH")

data_long <- tidyr::gather(numcells, condition, proportion, Control:EtOH, factor_key=TRUE)
data_long$Cluster <- factor(data_long$Cluster, levels = numcells$Cluster)

max_y <- ceiling(max(data_long$proportion)+10)

DimPlot(seurat_OBfree, cells.highlight = (WhichCells(seurat_OBfree, idents = "Migrating Neurons")), cols.highlight = "#7bafde", pt.size = 0.1, sizes.highlight = 0.1)

print(ggplot(data_long[c(1, 15),], aes(fill=Cluster, y=proportion, x=Cluster, alpha = condition)) +  #(fill=condition, y=proportion, x=Cluster))
        geom_bar(position="dodge", stat="identity", color="black", width = 0.70) +
        scale_alpha_discrete(range= c(1, 0.65)) +
        theme_classic() +
        ylim(0,max_y) +
        ylab("Proportion of cells") +
        scale_fill_manual(values=c("#7bafde"))+
        ggtitle("scRNA-Seq - Proportion of cells per cell type") +
        geom_text(aes(label = paste0(round(proportion, 2), " %")),
                  position = position_dodge(width = 1), angle=90, hjust = -0.25, size = 3.2)  + RotatedAxis() + theme(legend.title=element_blank()))
