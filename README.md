# Yarrowia Multifactor: study on urea and ammonium sulphate utilization

## Abstract
Media components, including the nitrogen source, are significant cost factors in cultivation processes. The nitrogen source also influences the cell behaviour and its production performance. Ammonium sulphate is a widely used nitrogen source for the cultivation of microorganisms. Urea is a sustainable and cheap alternative nitrogen source. We investigated the influence of urea as a nitrogen source compared to ammonium sulphate by cultivating three phenotypically very different <i>Yarrowia lipolytica</i> strains in chemostats under carbon- or nitrogen limiting conditions. We found no significant coherent changes in growth and lipid production between nitrogen sources. RNA sequencing revealed no significant concerted changes in the transcriptome due to the different nitrogen sources. We also showed that the genes involved in urea uptake and degradation are not up-regulated on a transcriptional level. Our findings support urea usage as a nitrogen source, indicating that previous metabolic engineering efforts where ammonium sulfate was used are likely translatable to the usage of urea and can ease the way for urea as a cheap and sustainable nitrogen source in more applications. 

## Repository structure
This repository contains the transcriptome data and the scripts to reproduce the analysis published in the paper "Coherent transcriptional response of <i>Yarrowia lipolytica</i> to carbon and nitrogen limitation in urea and ammonium sulphate". [Insert DOI:]

Folder structure:
1. <b>code</b>: the code folder contains all the scripts that were used to generate the analysis. This folder contains a "Functions" subfolder in which functions used in the scripts are saved. The script "00.Run_All_Rmds.R" will run all the Rmd files present in the folder which will generate the folder structure of /results and will save the output of the analysis. 
2. <b>data</b>: stored in this folder are all the necessary data used to run the scripts and to reproduce the output of the analysis. 
3. <b>docs</b>: in this folder documentation on the packages, functions, and pipelines used can be found.
4. <b>results</b>: stored in this folder are the results from the scripts
