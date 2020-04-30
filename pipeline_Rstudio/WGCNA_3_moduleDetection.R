# **  WGCNA pipeline part 3:   **

# choose soft threshold power to creat scale-free network; module detection: module-traits correlation

# ###### input: ###########################################################
# datExpr2: normalized gene expression data (rownames--genes, colnames--samples)
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


# Load the data saved in the second part
lnames = load("topology.Rdata")
lnames

# Plot the topology results:
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
# 
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


# decide soft threshold power
# ------------------------ choose cut hight of samples -------------------
sft.power <- as.numeric(readline(prompt="Enter the soft threshold power to use: "))

# replot the topology results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=-sign(sft$fitIndices[sft.power,3])*sft$fitIndices[sft.power,2],col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=sft$fitIndices[sft.power,5],col="red")

# use WGCNA's cor() function
cat("\n---------------------module detection--------------------\n ")
cor <- WGCNA::cor
net = blockwiseModules(datExpr2, power = sft.power,
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
MEs0 = moduleEigengenes(datExpr2, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "WGCNA_network.RData")


#################### test module - traits correlation #################
cat("\n ---------------------- Calculate Modules - traits correlation -------------------\n ")
cor <- stats::cor
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(datExpr2))

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
               colors = blueWhiteRed(50),
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
               colors = blueWhiteRed(50),
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
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off ()

cat("Number of genes in each module: \n")
print(as.data.frame(table(moduleColors)))

# find the top hub gene in each module
topHubs <- chooseTopHubInEachModule(datExpr=datExpr2, colorh = moduleColors, power = 2)
print(topHubs)

# save WGCNA data
write.table(datExpr2, quote=FALSE, file="WGCNA_expressionData.tsv", col.names=NA, sep="\t")
write.table(datTraits, quote=FALSE, file="WGCNA_traitData.tsv", col.names=NA, sep="\t")
write.table(MEs, quote=FALSE, file="WGCNA_moduleEigengenes.tsv", col.names=NA, sep="\t")
