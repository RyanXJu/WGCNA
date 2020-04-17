###### Analyse the module of intrest #######
# ****** need output data from WGCNA_pipeline_moduleDetection.R *****


library(ggplot2)
library(gplots)
library(reshape2)
library(doParallel)
library(WGCNA)
options(stringsAsFactors = FALSE)


setwd("/u/juxiao/AML_WGCNA")
getwd()
############### load data ###########################################################
datExpr11 <- read.delim ("WGCNA_expressionData.tsv", header=TRUE, row.names=1)

datTraits <- read.delim ("WGCNA_traitData.tsv", header=TRUE, row.names=1)

netInfo <- load("WGCNA_network.RData")

netInfo



############### Gene Module Membership ##############################################
cat("----------------- Calculate gene ModuleMembership --------------------")
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))  # correlation of the gene and the ME
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),  nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")