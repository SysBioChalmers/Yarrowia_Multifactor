# Description:
#     Performs consensus GSA using eight common (and relatively fast) methods that are
#     included in Piano.
#     Runs on maximum number of CPU nodes. Gives resList as output, which can be
#     used to make a plot using consGSAplot.
#
# Input:      TopTable  output from limma's topTable
#             gsc       gene-set, as loaded using Piano
#             nPerm     the number of permutations to use for gene sampling, default = 1000
#             gsSizeDn  cutoff-value for minimum number of genes in gene-sets, default = 5
#             gsSizeUp  cutoff-value for maximum number of genes in gene-sets, default = 500
#
# Output:     resList   results list
#
# 2021-03-24  Eduard Kerkhoven
# 2022-02-23  Modified by Simone Zaghen



consGSA <-
  function(TopTable,
           gsc,
           nPerm = 1000,
           gsSizeDn = 5,
           gsSizeUp = 500) {
    require(piano)
    require(parallel)
    require(snowfall)
    require(tidyr)
    
    Pval<-TopTable$adj.P.Val
    names(Pval)<-rownames(TopTable)
    FC<-TopTable$logFC
    names(FC)<-rownames(TopTable)
    
    # Find out number of processor cores, run GSA on all cores.
    cores <- as.numeric(detectCores())
    #cat("Run GSA on", cores, "CPU cores.\n")
    nPerm = round(nPerm / cores) * cores
    gsSizeLim = c(gsSizeDn, gsSizeUp)
    
    #cat("Running GSA 1/7 (mean)\n")
    gsaRes1 <- runGSA(Pval, FC, geneSetStat = "mean", gsc = gsc, gsSizeLim = gsSizeLim,
                      nPerm = nPerm, ncpus = cores, verbose = FALSE)
    #cat("Running GSA 2/7 (median)\n")
    gsaRes2 <- runGSA(Pval, FC, geneSetStat = "median", gsc = gsc, gsSizeLim = gsSizeLim,
                      nPerm = nPerm, ncpus = cores, verbose = FALSE)
    #cat("Running GSA 3/7 (sum)\n")
    gsaRes3 <- runGSA(Pval, FC, geneSetStat = "sum", gsc = gsc, gsSizeLim = gsSizeLim,
                      nPerm = nPerm, ncpus = cores, verbose = FALSE)
    #cat("Running GSA 4/7 (stouffer)\n")
    gsaRes4 <- runGSA(Pval, FC, geneSetStat = "stouffer", gsc = gsc, gsSizeLim = gsSizeLim,
                      nPerm = nPerm, ncpus = cores, verbose = FALSE)
    #cat("Running GSA 5/7 (tailStrength)\n")
    gsaRes5 <- runGSA(Pval, FC, geneSetStat = "tailStrength", gsc = gsc, gsSizeLim = gsSizeLim,
                      nPerm = nPerm, ncpus = cores, verbose = FALSE)
    # cat("Running GSA 6/9 (gsea)\n")
    # gsaRes6 <- runGSA(FC, geneSetStat = "gsea", gsc = gsc, gsSizeLim = gsSizeLim,
    #                   nPerm = nPerm, ncpus = cores)
    #cat("Running GSA 6/7 (fisher)\n")
    gsaRes6 <- runGSA(Pval, FC, geneSetStat = "fisher", gsc = gsc, gsSizeLim = gsSizeLim,
                      nPerm = nPerm, ncpus = cores, verbose = FALSE)
    #cat("Running GSA 7/7 (maxmean)\n")
    gsaRes7 <- runGSA(FC, geneSetStat = "maxmean", gsc = gsc, gsSizeLim = gsSizeLim,
                      nPerm = nPerm, ncpus = cores, verbose = FALSE)
    # cat("Running GSA 9/10 (fgsea)\n")
    # gsaRes9 <- runGSA(FC, geneSetStat = "fgsea", gsc = gsc, gsSizeLim = gsSizeLim,
    #                   ncpus = cores)
    # cat("Running GSA 9/9 (page)\n")
    # gsaRes9 <- runGSA(FC, geneSetStat = "page", gsc = gsc, gsSizeLim = gsSizeLim,
    #                   nPerm = nPerm, ncpus = cores)
    # No wilcoxon, and fgsea too slow. No gsea and page, no directional P-value as output for plot
    #cat("Reorganizing data and prepare for plotting")
    # Combine results in list
    resList <- list(gsaRes1, gsaRes2, gsaRes3, gsaRes4, gsaRes5, gsaRes6, gsaRes7)
    resList <-
      setNames(resList,
               c("mean", "median", "sum", "stouffer",
                 "tailStrength", "fisher", "maxmean"))
    
    return(resList)
  }