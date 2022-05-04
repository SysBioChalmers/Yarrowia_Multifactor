# Description:
#     Generates a list containing genes up or down-regulated in the conditions specified
#     by the input Group. It also generates a vector with the up and down regulated gene names,
#     ready to be used for subletting a whole gene count dataframe.
#
# Input:    df.list         List of dataframes. Each dataframe is the result of a pairwise comparison
#                           and should contain at least the names of the genes, the p value and the 
#                           log fold change. 
#           Group           Optional. Character vector containing the names of the dataframes that need to be 
#                           selected from df.list.
#           pval_cutoff     Cutoff for  P-values of differentially expressed genes.
#           neg_FC          Threshold for the low bound of fold change. 
#           pos_FC          Threshold for the high bound of fold change. 
#
# Output:   output          List containing a list with the up-regulated genes and a list for the 
#                           downregulated genes. Both these lists will contain dataframes with the
#                           name of the comparison used, the genes, the p value, and the fold change
#                           The function will also generate a vector containing the genes up or down
#                           regulated in all the conditions.
#                           The function will generate two plot_ objects which are eulerr diagrams,
#                           one for the up, and one for the down - regulated genes. 
#
# 2021/12/27  Simone Zaghen

#Function that receive as input a list containing multiple data frames, 
#a Group vector (vector of characters) a p-value cutoff, 
#and a lower (neg_FC) and a higher (pos_FC) fold change. 
#Based on these values, it  
mining <- function (df.list, 
                     Group, 
                     pval_cutoff, 
                     neg_FC, 
                     pos_FC)
{
  if (!require('eulerr')) install.packages('eulerr'); library('eulerr')
  if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')

  #Function takes df.list as an input, and it needs to be a list of dataframes
  if (is.list(data.frame(df.list)) != T) {
    stop("input must be a list of dataframes")}
  
  #Check if a Group is defined to filter the original dataframe list
  if (exists("Group")) {
    tmp = df.list[Group] #select only conditions specified in the Group vector
  }
  
  tmp = map(tmp, ~ select(.x, c("genes", "logFC", "P.Value", "adj.P.Val"))) #select only relevant columns
  
  #filtering step: filter up regulated in all conditions
  up = map(tmp, ~ filter(.x, 
                         adj.P.Val < pval_cutoff & 
                           logFC > pos_FC)) 
  
  #filtering step: filter down regulated in all conditions
  dn = map(tmp, ~ filter(.x, 
                         adj.P.Val < pval_cutoff & 
                           logFC < neg_FC))
  
  #extract genes that are up in all conditions
  list_up = list()
  for (i in 1:length(up)) {
    tmp_up <- up[[i]][[1]]
    name <- names(up[i])
    list_up[[name]] <- tmp_up
  }
  common_up = Reduce(intersect, list_up)
  all_up = Reduce(union, list_up)
  
  #extract genes that are down in all conditions
  list_dn = list()
  for (i in 1:length(dn)) {
    tmp_dn <- dn[[i]][[1]]
    name <- names(dn[i])
    list_dn[[name]] <- tmp_dn
  }
  common_down = Reduce(intersect, list_dn)
  all_down = Reduce(union, list_dn)
  
  all_DE = union(all_up, all_down)
  
  plot_up = plot(euler(list_up), 
                 quantities = TRUE,
                 main = "Upregulated",
                 legend = list(side="top"))
  
  plot_dn = plot(euler(list_dn),
                 quantities = TRUE,
                 main = "Downregulated",
                 legend = list(side="top"))
  
  #Generate the output list 
  output = list(up=up, 
                dn=dn,
                common_up = common_up,
                all_up = all_up,
                all_down = all_down,
                all_DE = all_DE,
                common_down = common_down,
                plot_up = plot_up,
                plot_dn = plot_dn)
  
  return(output) #only return the "output" object
  
}