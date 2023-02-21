# Description:
#     Performs consensus GSA using seven common (and relatively fast) methods that are included in Piano.
#     Runs on maximum number of CPU nodes. 
#     Gives resList as output, which can be used to make a plot using consGSAplot.
#
# Input:      TopTable  output from limma's topTable
#             gsc       gene-set, as loaded using Piano
#             nPerm     the number of permutations to use for gene sampling, default = 1000
#             gsSizeDn  cutoff-value for minimum number of genes in gene-sets, default = 5
#             gsSizeUp  cutoff-value for maximum number of genes in gene-sets, default = 500
#             title     short title used in plot and filename
#
# Output:     resList   results list
# 
# 2021-03-24 Eduard Kerkhoven
# 2022-06-29 Simone Zaghen

consGSA <-
  function(TopTable,
           gsc,
           nPerm = 1000,
           gsSizeDn = 5,
           gsSizeUp = 500) 
{

  suppressPackageStartupMessages({
    require(piano)
    require(parallel)
    require(snowfall)
    require(tidyr)
  })

  # Find out number of processor cores, run GSA on all cores.
    cores <- as.numeric(detectCores())
    nPerm = round(nPerm / cores) * cores
    
  # Parameters for GSA
    Pval <- TopTable$adj.P.Val
    names(Pval) <- rownames(TopTable)
    FC <- TopTable$logFC
    names(FC) <- rownames(TopTable)
    gsSizeLim = c(gsSizeDn, gsSizeUp) #Set gene set minimum and maximum size
    
  #Which GSA methods to run
    GSA_methods <- c("mean", "median", "sum", "stouffer", "tailStrength", "fisher", "maxmean")

  #Run consensous GSA  
  resList <- list() #create list in which to store data

  for(i in 1:length(GSA_methods)){
    
    method <- GSA_methods[i] #select method to run
    
    #maxmean does not accept Pval as input; if-else requirement needed for this method
      if(GSA_methods[i] != "maxmean"){
        gsaRes <- runGSA(Pval, 
                         FC, 
                         geneSetStat = method, 
                         gsc = gsc, 
                         gsSizeLim = gsSizeLim,
                         nPerm = nPerm, 
                         ncpus = cores, 
                         verbose = FALSE)
        }
      else { #else statement to run the maxmean
        gsaRes <- runGSA(FC, 
                         geneSetStat = method, 
                         gsc = gsc, 
                         gsSizeLim = gsSizeLim,
                         nPerm = nPerm, 
                         ncpus = cores, 
                         verbose = FALSE)}
          resList[[i]]<- gsaRes
        }

    resList <-setNames(resList, GSA_methods)
    
  return(resList)
  }