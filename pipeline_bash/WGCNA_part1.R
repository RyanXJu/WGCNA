# **  WGCNA pipeline part 1:   **

# remove genes with too many missing values, or low variance
# build sample tree with kept genes to detect outliers

# ###### input: ###########################################################
# datExpr: normalized gene expression data (rownames--genes, colnames--samples)
# datTrait: trait data of all the samples (rownames--samples, colnames--trait)
#################################################################################


## module load R/3.6.1
# .libPaths("/u/juxiao/R/x86_64-pc-linux-gnu-library/3.6" )
# .libPaths()

cat("\n
    *********************************\n
    *********  WGCNA part1  *********\n
    ******* Sample clustering *******\n
    *********************************")

library(ggplot2)
library(gplots)
library(reshape2)
library(doParallel)
library(WGCNA)
options(stringsAsFactors = FALSE)


cat("Enter working directory: ")
directory <- readLines("stdin", 1)

dir.create(file.path(directory), showWarnings = FALSE)
setwd(file.path(directory))
getwd()

################## load input files ##################################


# cat("Enter path of the expression data : ")
# fileExpr <- readLines("stdin", 1)

fileExpr <- ""
while (!file.exists(fileExpr)) {
  cat("Enter path of the expression data : ")
  fileExpr <<- readLines("stdin", 1)
  if (!file.exists(fileExpr)) { cat("File doesn't exist. ")}
}
print(fileExpr)

cat("\n-------------Loading expression data -----------\n")
datExpr <- read.delim(fileExpr,header=TRUE, row.names=1) 
dim(datExpr)

# cat("Enter path of the trait data ")
# fileTraits <- readLines("stdin", 1)

fileTraits <- ""
while (!file.exists(fileTraits)) {
  cat("Enter path of the trait data : ")
  fileTraits <<- readLines("stdin", 1)
  if (!file.exists(fileTraits)) { cat("File doesn't exist. ")}
}
print(fileTraits)
cat("\n-------------Loading trait data -----------------\n")
datTraits <- read.delim(fileTraits, header = TRUE, row.names = 1)
dim(datTraits)

# ask if user want to log tranform expression data
# WGCNA recommends a variance-stabilizing transformation with normalized counts (or RPKM/FPKM/TPM data) using log2(x+1). 
# For highly expressed features, the differences between full variance stabilization and a simple log transformation are small.
useLog = ""
while (!(useLog %in% c("Y", "y", "N", "n"))) {
  cat("Log tranform expression data [Y/N]: ")
  useLog <- readLines("stdin", 1)
  if (!(useLog %in% c("Y", "y", "N", "n"))) { cat("Log tranform expression data [Y/N]:")}
}

if (useLog == "Y"| useLog == "y"){
  cat("logarithm base [2/10]: ")
  log.base = as.numeric(readLines("stdin", 1))
  cat("add constant: ")
  log.add = as.numeric(readLines("stdin", 1))
  datExpr = log(datExpr+log.add, base=log.base)
}

datExpr1 <- t(datExpr)

################### data verification ################################
cat("\n---------------------Remove no good genes --------------------\n")

# verify expression data
gsg = goodSamplesGenes(datExpr1, verbose = 3)
# goodSamplesGenes() default parameters
# minFraction = 1/2 :minimum fraction of non-missing samples for a gene
# minNSamples =4 minimum number of non-missing samples for a gene
# remove genes with zero variance

gsg$allOK

## remove no good genes
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(colnames(datExpr1)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));

  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}

####################### remove outlier samples ######################
cat("\n----------------- Build sample tree --------------------------\n")

sampleTree = hclust(dist(datExpr1, method = "euclidean"), method = "average")

# png("Part1_SampleClustering.png", width=1000, height=500)
# par(mar = c(1,6,2,1))
# plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
#      cex.axis = 1.2, cex.main = 2)

pdf("Part1_SampleClustering.pdf", width=12, height=8)
par(mar = c(1,6,2,1))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.2, cex.main = 2)

save(datExpr1, datTraits, sampleTree, file = "sampleTree.Rdata")
cat("\n------------------- Part1 Done ------------------------------\n")