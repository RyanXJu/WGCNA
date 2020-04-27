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

directory <- "~/AML_WGCNA/test_pipeline"

library(ggplot2)
library(gplots)
library(reshape2)
library(doParallel)
library(WGCNA)
options(stringsAsFactors = FALSE)

dir.create(file.path(directory), showWarnings = FALSE)
setwd(file.path(directory))
getwd()


################## load input files ##################################
cat("\n---------------------Loading expression data and trait data--------------------\n ")

datExpr <- read.delim("datExpr.tsv",header=TRUE, row.names=1) 
dim(datExpr)
datTraits <- read.delim("datTraits.tsv", header = TRUE, row.names = 1)
dim(datTraits)

# ask if user want to log tranform expression data
# WGCNA recommends a variance-stabilizing transformation with normalized counts (or RPKM/FPKM/TPM data) using log2(x+1). 
# For highly expressed features, the differences between full variance stabilization and a simple log transformation are small.
useLog = ""

while (!(useLog %in% c("Y", "y", "N", "n"))) {
  useLog <- readline(prompt="Log tranform expression data [Y/N]: ")
  if (!(useLog %in% c("Y", "y", "N", "n"))) { cat("Log tranform expression data [Y/N]:")}
}

if (useLog == "Y"| useLog == "y"){
  log.base = as.numeric(readline(prompt="logarithm base [2/10]: "))
  log.add = as.numeric(readline(prompt = "add constant: "))
  datExpr = log(datExpr+log.add, base=log.base)
}

datExpr1 <- t(datExpr)
  

################### data verification ################################
cat("\n---------------------Remove no good genes and outlier samples--------------------\n ")

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
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(colnames(datExpr1)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
}
dim(datExpr1)

####################### remove outlier samples ######################
# build sample tree based on chosen genes 
sampleTree = hclust(dist(datExpr1, method = "euclidean"), method = "average")

# Plot the sample tree
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6)
par(mar = c(1,6,2,1))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

png("Sample_GE_clustering.png", width=1000, height=500)
par(mar = c(1,6,2,1))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.2, cex.main = 2)
dev.off()

pdf("Sample_GE_clustering.pdf", width=12, height=8)
par(mar = c(1,6,2,1))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.2, cex.main = 2)
dev.off ()


save(datExpr1, datTraits, sampleTree, file = "sampleTree.Rdata")
