# Description:
#     Makes a plot showing the most significant change genesets, as identified
#     from consensus geneset analysis using Piano. For the highest ranking
#     genesets, a barchart is plotted showing the number of up- and down
#     regulated genes. A p-value cutoff is used to signify what genes have
#     significantly changed expression.
#
# Input:      resList   result list from consensus GSA, ideally from consGSA.R
#             rankScore cutoff for non-directional consensus score, whether
#                       genesets should be included in barchart.
#             Pcutoff   cutoff for adjusted P-values of differentially expressed
#                       genes. default is 0.05.
#             distinct  'distinct' if distinct-up and -down should be used for rank
#                       cutoffs. 'distinct' will give less results, but more clearly
#                       up or down regulated. 'mixed' will give more results, but
#                       this will also include GO-terms with mixed regulation.
#             savePlot  true if plot should be saved as file, false if not
#             title     short title used in plot and filename
# Output:     plot      the plot
#
# 2018-03-13  Eduard Kerkhoven
# Update: gives the table as an output instead of the plot - Oliver Konzock


consGSAtable <-
  function(resList,
           rankScore,
           Pcutoff = 0.05,
           distinct = 'mixed') {
    library(ggplot2)
    library(piano)
    library(parallel)
    library(snowfall)
    library(tidyr)
    library(scales)
    
    # Extract non-directional, distinct up and down genesets.
    non <- consensusScores(resList, class = "non", plot = F)
    df <- data.frame(Name=rownames(non$rankMat[non$rankMat[,1] < rankScore + 1,])) # Select any non-directional with rank =< rankScore.
    
    # including the RankScore in the XXX to be able to plot the RankScore in the final figure as well
    df4 <- data.frame(Name=rownames(non$rankMat[non$rankMat[,2] < rankScore + 1,])) #I had to change the 1 to a 2 here to fix my bug. 
    rank <- as.data.frame(non$rankMat[non$rankMat[, "ConsScore"] < rankScore + 1 , , drop = F])
    df4$rank <- rank$ConsScore
    
    colnames(df4) <- c("Name", "RankScore")
    
    Pval <- resList[[1]]$geneLevelStats
    FC <- resList[[1]]$directions
    GS <- resList[[1]]$gsc

    
if (distinct=='distinct') {
  up <-
    consensusScores(resList,
                    class = "distinct",
                    direction = "up",
                    plot = F)
  dn <-
    consensusScores(resList,
                    class = "distinct",
                    direction = "down",
                    plot = F)
  non <-
    non$rankMat[non$rankMat[, "ConsScore"] < rankScore + 1 , , drop = F]
  dn <-
    dn$rankMat[dn$rankMat[, "ConsScore"] < rankScore + 1, , drop = F]
  up <-
    up$rankMat[up$rankMat[, "ConsScore"] < rankScore + 1, , drop = F]
  dn <- dn[rownames(dn) %in% rownames(non), , drop = F]
  up <- up[rownames(up) %in% rownames(non), , drop = F]

  df <- data.frame(Name = unique(c(rownames(up), rownames(dn))))
}

    sumTable<-data.frame(Name=names(resList[[1]]$gsc), up=resList[[1]]$nGenesUp, dn=resList[[1]]$nGenesDn, tot=resList[[1]]$nGenesTot)
    df <- merge(df, sumTable)
    df <- merge(df,df4)
    colnames(df) <- c("geneset", "up", "dn", "tot", "rank")
    
    df$highup <- 0
    df$highdown <- 0
    df$lowup <- 0
    df$lowdown <- 0
    
    tmp <- GS[names(GS) %in% df$geneset]
    tmp <- tmp[order(names(tmp))]
    
    for (gset in 1:length(tmp)) {
      df[gset,]$highup <- sum(tmp[[gset]] %in% names(Pval[Pval < Pcutoff,]) &
                                tmp[[gset]] %in% names(FC[FC > 0,]))  # Highly significant up
      df[gset,]$highdown <- sum(tmp[[gset]] %in% names(Pval[Pval < Pcutoff,]) &
                                  tmp[[gset]] %in% names(FC[FC < 0,]))  # Highly significant down
      df[gset,]$lowup <- sum(tmp[[gset]] %in% names(Pval[Pval > Pcutoff,]) &
                               tmp[[gset]] %in% names(FC[FC > 0,]))  # Low significant up
      df[gset,]$lowdown <- sum(tmp[[gset]] %in% names(Pval[Pval > Pcutoff,]) &
                                 tmp[[gset]] %in% names(FC[FC < 0,]))  # Low significant down
    }
    
    df$geneset <- factor(df$geneset, levels = df[order(df$up / df$dn),1])  # Order from mostly up to mostly down.
    #df$geneset<-gsub('(.{100})(.+)','\\1...',df$geneset)
    
    # Add ALL genes
    df3 <- data.frame(geneset = 'All')
    df3$geneset <- 'All'
    df3$up <- sum(FC > 0)
    df3$dn <- sum(FC < 0)
    df3$tot <- length(FC)
    df3$highup <- sum(Pval < Pcutoff & FC > 0) 
    df3$highdown <- sum(Pval < Pcutoff & FC < 0) 
    df3$lowup <- sum(Pval > Pcutoff & FC > 0) 
    df3$lowdown <- sum(Pval > Pcutoff & FC < 0)
    df3$rank <- ""
    df<-rbind(df,df3)

    return(df)
  }

