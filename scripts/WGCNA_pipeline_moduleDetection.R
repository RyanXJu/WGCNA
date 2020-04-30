###### WGCNA pipeline part1: network construction and module detection #######


#----------------------------------------------------
###### input files: #######
# datExpr: normalized gene expression data (rownames--genes, colnames--samples)
# datTrait: trait data of all the samples (rownames--samples, colnames--trait)

###### arguments:  ########
# sft.auto : If soft threshold is selected automaticlly
# outlier.rm: If outlier samples are removed automaticlly (WGCNA suggest remove outliers)
###########################

## module load R/3.6.1
# .libPaths("/u/juxiao/R/x86_64-pc-linux-gnu-library/3.6" )
# .libPaths()

library(ggplot2)
library(gplots)
library(reshape2)
library(doParallel)
library(WGCNA)
options(stringsAsFactors = FALSE)

directory <- "~/AML_WGCNA/test_pipeline"
dir.create(file.path(directory), showWarnings = FALSE)
setwd(file.path(directory))

getwd()

sft.auto <- TRUE
outlier.rm <- TRUE

################## load input files ##################################
cat("\n---------------------Loading expression data and trait data--------------------\n ")

datExpr <- read.delim("datExpr.tsv",header=TRUE, row.names=1) 
dim(datExpr)

datTraits <- read.delim("datTraits.tsv", header = TRUE, row.names = 1)
dim(datTraits)

################### data verification ################################
cat("\n---------------------Remove no good genes and outlier samples--------------------\n ")
# ??????????????????????????????????????????????????????????????  Check if data is already logged ??????????????????????????????
# take log2 +1, transpose datExpr for WGCNA
# datExpr0 <- t(log2(datExpr+1))
datExpr0 <- t(datExpr)

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

####################### remove outlier samples ######################
# build sample tree based on chosen genes 
sampleTree = hclust(dist(datExpr0, method = "euclidean"), method = "average")

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


# decide cut hight
# ------------------------ choose cut hight of samples -------------------
cutHeight <- as.numeric(readline(prompt="Enter the hight to cut outliers: "))

# Plot the cut line
abline(h = cutHeight, col = "red");

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = cutHeight, minSize = 10)
table(clust)

# choose the cluster with the biggest numbers of samples
keepSamples = (clust== names(which.max(table(clust))))

datExpr = datExpr0[keepSamples, ]
datTraits = datTraits[keepSamples, ] 

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
nGenes
nSamples

outlierSamples <- !keepSamples
outlierSamples_names<- rownames(datExpr0[!keepSamples, ])
outlierSamples_names
cat("Romved outlier samples: \n")
cat( outlierSamples_names)

# Re-cluster samples
cat("\n ------------------- Recalculate sampleTree ------------------- \n")
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")


############ network construction (soft-threshold selection) #########

cat("\n---------------------Select soft-threshold--------------------\n ")

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))

# parallize if possible
registerDoParallel()

# Call the network topology analysis function
# ********** for 55000 genes, takes about 1 hour *********
# ********** for 16000 genes, takes about 10 minutes ******
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
# abline(h=0.75,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# save plot as png and pdf 
png("ScaleFree_topology.png", width=1000, height=500)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

pdf("ScaleFree_topology.pdf", width=12, height=8)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off ()


# decide soft threshold power
# ------------------------ choose cut hight of samples -------------------
sft.power <- as.numeric(readline(prompt="Enter the soft threshold power to use: "))

# use WGCNA's cor() function
cat("\n---------------------module detection--------------------\n ")
cor <- WGCNA::cor
net = blockwiseModules(datExpr, power = sft.power,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "TPM_rsem_auto", 
                       verbose = 3)
cor<-stats::cor
table(net$colors) 

print("---------------finish detecting modules--------------------- ")

# Plot the dendrogram and the module colors underneath
sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

png("Gene_Cluster_Dendrogram.png", width=1000, height=500)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

pdf("Gene_Cluster_Dendrogram.pdf", width=12, height=8)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off ()

## save result
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
table(moduleColors)

# # recalculate module eigengenes by color
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "WGCNA_network.RData")


#################### test module - traits correlation #################
cat("\n ---------------------- Calculate Modules - traits correlation -------------------\n ")
cor <- stats::cor
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(3, 10, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

png("Module_trait_correlation.png", width=1000, height=500)
par(mar = c(3, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

pdf("Module_trait_correlation.png.pdf", width=12, height=8)
par(mar = c(3, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off ()

cat("Number of genes in each module: \n")
print(as.data.frame(table(moduleColors)))

# find the top hub gene in each module
topHubs <- chooseTopHubInEachModule(datExpr=datExpr, colorh = moduleColors, power = 2)
print(topHubs)

# save WGCNA data
write.table(datExpr, quote=FALSE, file="WGCNA_expressionData.tsv", col.names=NA, sep="\t")
write.table(datTraits, quote=FALSE, file="WGCNA_traitData.tsv", col.names=NA, sep="\t")
write.table(MEs, quote=FALSE, file="WGCNA_moduleEigengenes.tsv", col.names=NA, sep="\t")

