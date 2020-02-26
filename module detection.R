#####  3. Automatic network construction and module detection  ########

library(WGCNA)

setwd("/u/juxiao/AML_WGCNA")
getwd()
options(stringsAsFactors = FALSE)


# Load the data saved in the 2.data cleaning
lnames = load(file = "TPM_resm_dataInput_avelog.RData")
#The variable lnames contains the names of loaded variables.
lnames

####################
# 1-step, Automatic network construction and module detection
####################

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

## Allow and disable multi-threading for certain WGCNA calculations
# enableWGCNAThreads(nThreads = NULL)
# disableWGCNAThreads()

# Call the network topology analysis function
# ********** for 54927 genes, takes about 1 hour *********
# ********** for 16000 genes, takes about 10 minutes
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)


# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.75,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

## ******** use auto mudule detect ********
## ******** with defult setting 16k genes takes 35 minutes ******
# cat(paste("net blsize 5000 begin   :", Sys.time(),sep=""),
#     file="net_test_20200219.txt",sep="\n", append=TRUE)
# WGCNA has its own cor function
cor <- WGCNA::cor
net = blockwiseModules(datExpr, power = 9,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "TPM_rsem_auto", 
                       verbose = 3)
# cat(paste("net blsize 5000 end   :", Sys.time(),sep=""),
#     file="net_test_20200219.txt",sep="\n", append=TRUE)
cor<-stats::cor
table(net$colors) ##%%%%%%%%%%%%%%%%%%%%%% found 26 modules

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


## save result
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "TPM_avelog-auto.RData")


###************** Test maxBlockSize = 20000 ****************
## ******** with maxBlockSize = 20000 setting 16k genes takes 3.5hours ******
# cat(paste("net blsize 20k begin   :", Sys.time(),sep=""),
#     file="net_test_20200219.txt",sep="\n", append=TRUE)
# WGCNA has its own cor function
cor <- WGCNA::cor
net20k = blockwiseModules(datExpr, power = 7,
                       maxBlockSize = 20000,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "TPM_rsem_auto_20k", 
                       verbose = 3)
# cat(paste("net blsize 20k end   :", Sys.time(),sep=""),
#     file="net_test_20200219.txt",sep="\n", append=TRUE)
cor<-stats::cor
table(net20k$colors) ##%%%%%%%%%%%%%%%%%%%%%% found 21 modules

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net20k$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net20k$dendrograms[[1]], mergedColors[net20k$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


## save result
moduleLabels = net20k$colors
moduleColors = labels2colors(net20k$colors)
MEs = net20k$MEs;
geneTree = net20k$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "TPM_avelog-auto20k.RData")
