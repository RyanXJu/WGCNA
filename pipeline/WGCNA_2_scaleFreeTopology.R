# **  WGCNA pipeline part 2:   **

# offer topology information for soft threshold choosing

# ###### input: ###########################################################
# samplTree.Rdata by WGCNA_1_detectOutliers.R
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

# Load the data saved in the first part
lnames = load("sampleTree.Rdata")
lnames

# Plot the sample tree
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6)
par(mar = c(1,6,2,1))
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

# choose the cluster with the biggest numbers of samples
keepSamples = (clust== names(which.max(table(clust))))

datExpr2 = datExpr1[keepSamples, ]
datTraits = datTraits[keepSamples, ] 

# nGenes = ncol(datExpr2)
# nSamples = nrow(datExpr2)
# nGenes
# nSamples

outlierSamples <- !keepSamples
outlierSamples_names<- rownames(datExpr1[!keepSamples, ])
outlierSamples_names
cat("Romved outlier samples: \n")
cat( outlierSamples_names)

# Re-cluster samples
cat("\n ------------------- Recalculate sampleTree ------------------- \n")
sampleTree2 = hclust(dist(datExpr2), method = "average")
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
sft = pickSoftThreshold(datExpr2, powerVector = powers, verbose = 5)

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


save(datExpr2, datTraits, powers, sft,  file = "topology.Rdata")