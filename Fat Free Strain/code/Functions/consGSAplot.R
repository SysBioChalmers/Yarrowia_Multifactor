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
# Output:     plot      the plot
#
# 2018-03-13  Eduard Kerkhoven
# Update: included RankScores on the side of the plot - Oliver Konzock


consGSAplot <-
  function(resList,
           rankScore,
           Pcutoff = 0.05,
           distinct = 'mixed',
           title = 'Consensus GSA barchart'){
    
    suppressPackageStartupMessages({
      require(ggplot2)
      require(piano)
      require(parallel)
      require(snowfall)
      require(tidyr)
      require(scales)
    })
    
  non <- consensusScores(resList, class = "non", plot = F) # Extract non-directional

    if (distinct=='mixed') { 
      
      #if mixed, select any non-directional with ConsScore =< rankScore.
      df <- data.frame(Name=rownames(non$rankMat[non$rankMat[, "ConsScore"] <= rankScore,])) 

      #if distinct, select up and down with rank score =< rankScore
    } else { 
      up <- consensusScores(resList, class = "distinct", direction = "up", plot = F)
      dn <- consensusScores(resList, class = "distinct", direction = "down", plot = F)
      
      #fitler for ConsScore <= rankScore
      dn <- dn$rankMat[dn$rankMat[, "ConsScore"] <= rankScore, , drop = F]
      up <- up$rankMat[up$rankMat[, "ConsScore"] <= rankScore, , drop = F] 
      
      #remove non-directional
      to_keep <- rownames(non$rankMat[non$rankMat[, "ConsScore"] <= rankScore, , drop = F])
      dn <- dn[rownames(dn) %in% to_keep, , drop = F]
      up <- up[rownames(up) %in% to_keep, , drop = F]
      
      df <- data.frame(Name = unique(c(rownames(up), rownames(dn))))
    }
    
  #extract number of up, down, and total genes in each geneset
  sumTable <- data.frame(Name = names(resList$mean$gsc),
                         up = resList$mean$nGenesUp, 
                         dn = resList$mean$nGenesDn, 
                         tot = resList$mean$nGenesTot)
    
  df <- merge(df, sumTable) #merge the two dataframes, only keeping common elements
    
  #include the RankScore value
  df_rankScore <- data.frame(Name=rownames(non$rankMat[non$rankMat[, "ConsScore"] <= rankScore,]),
                             RankScore = non$rankMat[non$rankMat[, "ConsScore"] <= rankScore, , drop = F][, "ConsScore"])
  rownames(df_rankScore) <- NULL
    
  #merge the dataframes to include the rankScore information
  df <- merge(df, df_rankScore)
  colnames(df) <- c("geneset", "up", "dn", "tot", "rank")
  
  #Now count how many significative genes are up or down regulated in each geneset
    #add additional columns in which to store number of significative up/down genes
    df <- cbind(df, highup = NA, highdown = NA, lowup = NA, lowdown = NA) 

    Pval <- resList$mean$geneLevelStats #extract the pvalues
    FC <- resList$mean$directions #extract the directionality of fold change
    GS <- resList$mean$gsc #extract the gene sets
    tmp <- GS[names(GS) %in% df$geneset] #only keep gene sets present in our analysis
    tmp <- tmp[order(names(tmp))] #sort them by name
    
    #populate the NAs with the actual values of how many genes are significantly up/dn
    
    for (i in 1:length(tmp)) {
      #Highly significant up
      df[i,]$highup <- sum(tmp[[i]] %in% names(Pval[Pval < Pcutoff,]) & tmp[[i]] %in% names(FC[FC > 0,]))
      # Highly significant down
      df[i,]$highdown <- sum(tmp[[i]] %in% names(Pval[Pval < Pcutoff,]) & tmp[[i]] %in% names(FC[FC < 0,]))  
      # Low significant up
      df[i,]$lowup <- sum(tmp[[i]] %in% names(Pval[Pval > Pcutoff,]) & tmp[[i]] %in% names(FC[FC > 0,]))
      # Low significant down
      df[i,]$lowdown <- sum(tmp[[i]] %in% names(Pval[Pval > Pcutoff,]) & tmp[[i]] %in% names(FC[FC < 0,]))  
    }
    
  #Add statistics for All genes
  df_all <- data.frame(geneset = 'All',
                       up = sum(FC > 0), 
                       dn = sum(FC < 0),
                       tot = length(FC),
                       highup = sum(Pval < Pcutoff & FC > 0),
                       highdown = sum(Pval < Pcutoff & FC < 0),
                       lowup = sum(Pval > Pcutoff & FC > 0),
                       lowdown = sum(Pval > Pcutoff & FC < 0),
                       rank = "")
  df <- rbind(df,df_all)
  
  
  #Plot
    df_plot <- gather(df[, c("geneset", "highup", "highdown", "lowup","lowdown")], "directAndSignif", "genes", 2:5)

    #Order things from most up to most down
    df_plot$geneset <- factor(df_plot$geneset, levels = df_plot[order(df$up / df$dn),1])  #Order from mostly up to mostly down.
    df_plot$directAndSignif <- factor(df_plot$directAndSignif, levels = c("highup", "lowup", "lowdown", "highdown"))

    #Prepare the plot
    plot <-
      ggplot(df_plot, aes(x = geneset, y = genes, fill = directAndSignif)) +
      geom_bar(position = position_fill(reverse = TRUE), stat = "identity") +
      coord_flip() +
      scale_x_discrete(labels = label_wrap(50)) +
      scale_y_continuous(labels = percent_format(), breaks = c(0,0.25, 0.5, 0.75, 1)) +
      scale_fill_manual(
        values = c("#963836", "#b57372", "#abbdd2", "#4682b4"),
        name = "Direction",
        labels = c(
          paste0("Up (p<", Pcutoff, ")"),
          paste0("Up (p>", Pcutoff, ")"),
          paste0("Down (p>", Pcutoff, ")"),
          paste0("Down (p<", Pcutoff, ")")))  +
      labs(x = "GO terms", y = "Genes") +
      theme_set(theme_bw()) +
      theme(
        text = element_text(size = 11),
        legend.position = "bottom",
        axis.text = element_text(colour = "black"),
        panel.grid = element_blank(),
        line = element_line(size = 0.25),
        plot.margin = unit(c(1,1, 1, 1), "points"),
        legend.key=element_blank(),
        panel.border=element_blank(),
        axis.ticks.y=element_blank()) +
      geom_text( #adds the total number of genes in the gene sets on the right
        data = df,
        aes(x = geneset,
            y = 1.04, label = tot),
        inherit.aes = FALSE,
        size = 3) +
      geom_text( # adds the rankScore of the geneset on the left
        data = df,
        aes(x = geneset,
            y = -0.04, label = rank),
        inherit.aes = FALSE,
        size = 3) +
      labs(title = title,
           subtitle="            <-- RankScore        /         total number of genes --> ")
  

  output = list(df = df,
                plot = plot)
    
  return(output)
    
  }

