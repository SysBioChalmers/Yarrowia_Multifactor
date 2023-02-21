#Installs the packages required to perform the analysis

# Install and load packages from bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("edgeR", quietly = TRUE)) BiocManager::install("edgeR"); library(edgeR)
if (!require("limma", quietly = TRUE)) BiocManager::install("limma"); library(limma)
if (!require("piano", quietly = TRUE)) BiocManager::install("piano"); library(piano)
if (!require("EnhancedVolcano", quietly = TRUE)) BiocManager::install("EnhancedVolcano"); library(EnhancedVolcano)

#Install and load packages from CRAN
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('RColorBrewer')) install.packages('RColorBrewer'); library('RColorBrewer')
if (!require('ggpubr')) install.packages('ggpubr'); library('ggpubr')
if (!require('gplots')) install.packages('gplots'); library('gplots')
if (!require('gridExtra')) install.packages('gridExtra'); library('gridExtra')
if (!require('eulerr')) install.packages('eulerr'); library('eulerr')
if (!require('snowfall')) install.packages('snowfall'); library('snowfall')