library(readxl)
library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set wd

uniprot_table <- read.delim("../data/Uniprot_Yali.txt")
excel_input <- read_excel("../../00_Paper/Urea_genes.xlsx")

genes_to_extract <- as.character(as.vector(excel_input$YALI1))

urea_genes <- 
  uniprot_table %>%
  filter(Gene.names %in% unique(grep(paste(genes_to_extract,collapse="|"), 
                                     uniprot_table$Gene.names, value=TRUE)))
urea_genes$Gene.names <- str_extract(urea_genes$Gene.names, ".{0,0}YALI1_.{0,7}") #only keep YALI1_ gene names

#merge with FC and p value
comparison.list <- c("OKYL029_AS_116_0.1vsOKYL029_U_116_0.1",
                     "OKYL049_AS_116_0.1vsOKYL049_U_116_0.1",
                     "JFYL007_AS_116_0.1vsJFYL007_U_116_0.1",
                     "OKYL029_AS_3_0.1vsOKYL029_U_3_0.1",
                     "OKYL049_AS_3_0.1vsOKYL049_U_3_0.1",
                     "JFYL007_AS_3_0.1vsJFYL007_U_3_0.1")
comparison.list <- paste0(comparison.list, ".csv")

setwd("../results/Comparison_Output/")
df.list <- sapply(comparison.list, read.csv2, simplify = FALSE) #import comparison outputs
names(df.list) <- gsub(".csv", "", names(df.list)) #remove .csv from the comparison list

df.list <- 
  df.list %>%
  map(~filter(.x, genes %in% genes_to_extract)) %>%
  map(~select(.x, c("genes", "P.Value", "logFC")))
df.list <- mapply(cbind, df.list, "ComparisonID"= names(df.list), SIMPLIFY=F)

df <- bind_rows(df.list)

df <- merge(x = df, y = urea_genes, by.x = "genes", by.y = "Gene.names")
