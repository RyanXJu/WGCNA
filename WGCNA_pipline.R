###### WGCNA #######
#----------------------------------------------------
###### input files: #######
# datExpr: normalized gene expression data (rownames--genes, colnames--samples)
# datTrait: trait data of all the samples (rownames--samples, colnames--trait)

###### arguments:  ########
# sft.auto : If soft threshold is selected automaticly
# outlier.rm: If outlier samples are removed automaticly (WGCNA suggest remove outliers)
###########################

library(ggplot2)
library(reshape2)
library(WGCNA)
options(stringsAsFactors = TRUE)

setwd("/u/juxiao/AML_WGCNA")
getwd()

sft.auto <- TRUE
outlier.rm <- TRUE

################## load input files ##################################

datExpr <- datExpr[keep,]
dim(datExpr)

dim(datTraits)

############ data verification ######################
# take log2 +1, transpose datExpr for WGCNA
datExpr0 <- t(log2(datExpr+1))

# verify expression data
gsg = goodSamplesGenes(datExpr0, verbose = 3)
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
    printFlush(paste("Removing genes:", paste(colnames(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
dim(datExpr0)

##################### check outlier samples #################
# build sample tree based on chosen genes 
sampleTree = hclust(dist(datExpr0, method = "euclidean"), method = "average")

# Plot the sample tree
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# decide cut hight
# ------------------------ choose cut hight of samples -------------------
cutHeight <- as.numeric(readline(prompt="Enter the hight to cut outliers: "))

# Plot the cut line
abline(h = cutHeight, col = "red");

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = cutHeight, minSize = 10)
table(clust)

# choose the cluster with the most numbers of samples
keepSamples = (clust== names(which.max(table(clust))))
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
nGenes
nSamples

# clust 1 contains the samples we want to keep.
outlierSamples <- !keepSamples
outlierSamples_names<- rownames(datExpr0[!keepSamples, ])
