#Installs the packages required to perform the analysis

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("limma",
                       "edgeR",
                       "piano", 
                       "EnhancedVolcano"))

install.packages(c("tidyverse",
                   "RColorBrewer",
                   "ggpubr",
                   "gplots",
                   "gridExtra",
                   "eulerr",
                   "snowfall"))
