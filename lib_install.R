##Bioconductor: https://bioconductor.org/

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("impute")
BiocManager::install("biomaRt")
BiocManager::install("topGO")
BiocManager::install("KEGGREST")


##Other

install.packages("ggplot2")
install.packages("reshape2")
install.packages("doParallel")
install.packages("WGCNA")

                     