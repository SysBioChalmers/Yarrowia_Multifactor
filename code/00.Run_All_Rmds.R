setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set wd
files <- list.files(pattern = "[.]Rmd$") #get Rmd file names

#render the Rmd file names
for (f in files) 
  rmarkdown::render(f, 
                    encoding = encoding, 
                    output_dir = "../results/Reports") 
