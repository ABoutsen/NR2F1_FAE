library(ggplot2)
library(Seurat)

seurat_OBfree <- readRDS("../../seurat_OBfree.rds")
seurat_OBfree <- SetIdent(seurat_OBfree, value = "Subcluster")

MGE <- WhichCells(seurat_OBfree, idents = c("MGE Interneurons")) 
CGE  <- WhichCells(seurat_OBfree, idents = c("CGE Interneurons")) 

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

ggplot(data_long[c(5,8,19,22),], aes(fill=Cluster, y=proportion, x=Cluster, alpha = condition)) +
        geom_bar(position="dodge", stat="identity", color="black", width = 0.70) +
        scale_alpha_discrete(range= c(1, 0.65)) +
        theme_classic() +
        ylim(0,10) +
        ylab("Proportion of cells") +
        scale_fill_manual(values=c(colors_PN[1], colors_PN[2]))+
        ggtitle("scRNA-Seq - Proportion of cells per cell type") +
        geom_text(aes(label = paste0(round(proportion, 2), " %")),
                  position = position_dodge(width = 1), angle=90, hjust = -0.25, size = 3.2)  + RotatedAxis() + theme(legend.title=element_blank())

DimPlot(seurat_OBfree, label=F, sizes.highlight = 0.3, cells.highlight=list("MGE"=MGE, "CGE"=CGE), cols.highlight = c(colors_PN[1], colors_PN[2]), cols="grey")
