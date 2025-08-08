# NR2F1_FAE
Code and annotations for the "NR2F1 is a Cross-Species Mediator of Cerebral Cortex Defects Induced by Fetal Alcohol Intoxication" paper single-cell transcriptomic dataset. 

# Single Cell RNA sequencing

## Installation - R

To re-do the analyses, R 4.1.1 and Rtools40 should be used. And the following packages should be installed 

```r
remotes::install_version("spatstat.utils", version = "2.2-0")
remotes::install_version("spatstat.data", version = "2.1-0")
remotes::install_version("spatstat.geom", version = "2.3-0")
remotes::install_version("spatstat.core", version = "2.3-0")
remotes::install_version("Seurat", version = "4.0.5")
remotes::install_github("satijalab/seurat-wrappers")

install.packages(c("ggplot2", "tidyverse", "ggrepel", "dplyr", "RColorBrewer", "monocle3"))
BiocManager::install("biomaRt")
devtools::install_github('cole-trapnell-lab/monocle3')
```

## Getting started

### From raw fastQ

If you want to start from the raw fastQ files, then first download the data from .... 

### From "raw" filtered feature matrices

If you want to start from the preprocessed filtered features matrixes you can download data from ... . and follow the codes provided in this deposetory to perform the analysis. 

# Chromatin Immunoprecipitation sequencing (ChIP-seq)
