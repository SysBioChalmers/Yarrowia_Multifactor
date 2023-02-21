#This script runs all the R markdown files in the folder, generating 
#the folder structure, the files, and the plots used for analysis

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set working directory
files <- list.files(pattern = "[.]Rmd$") #get Rmd file names

if (!dir.exists("../results/")){dir.create("../results/")}
if (!dir.exists("../results/Reports")){dir.create("../results/Reports")}

#render the Rmd files and perform the analysis
for (f in files) 
  rmarkdown::render(f, 
                    encoding = encoding, 
                    output_dir = "../results/Reports") 