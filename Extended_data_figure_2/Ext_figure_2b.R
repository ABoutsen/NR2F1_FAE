
library(ggplot2)
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

num_control <- table(hfile@meta.data[which(hfile@meta.data$type=="Ctrl"), "Population_1"])
num_etoh <- table(hfile@meta.data[which(hfile@meta.data$type=="EtOH"), "Population_1"])
num_total <- table(hfile@meta.data["Population_1"] )

numcells <- as.data.frame(cbind(Control = num_control,EtOH = num_etoh, Prop_Control = num_control/sum(num_control)*100, Prop_EtOH = num_etoh/sum(num_etoh)*100))
numcells <- numcells[order(-numcells$Control),]

numcells$Cluster <- rownames(numcells)
numcells2 <- numcells
numcells <- numcells[, c(5,3,4)]
colnames(numcells) <- c("Cluster", "Control", "EtOH")

data_long <- tidyr::gather(numcells, condition, proportion, Control:EtOH, factor_key=TRUE)
data_long$Cluster <- factor(data_long$Cluster, levels = numcells$Cluster)

max_y <- ceiling(max(data_long$proportion)+10)

print(ggplot(data_long, aes(fill=Cluster, y=proportion, x=Cluster, alpha = condition)) +  #(fill=condition, y=proportion, x=Cluster))
        geom_bar(position="dodge", stat="identity", color="black", width = 0.70) +
        scale_alpha_discrete(range= c(1, 0.65)) +
        theme_classic() +
        ylim(0,max_y) +
        ylab("Proportion of cells") +
        scale_fill_manual(values=c("#4eb265", "#7bafde", "#ba8bd4", "#cae0ab", "#1965b0", "#994f88", "#ee8026", "#f7cb45", "#e8ecfb", "#dc050c", "deeppink"))+
        ggtitle("scRNA-Seq - Proportion of cells per cell type") +
        geom_text(aes(label = paste0(round(proportion, 2), " %")),
                  position = position_dodge(width = 1), angle=90, hjust = -0.25, size = 3.2)  + RotatedAxis() + theme(legend.title=element_blank()))
