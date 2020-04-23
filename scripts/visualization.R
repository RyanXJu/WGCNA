############ 5 module visualization #####

library(WGCNA)

setwd("/u/juxiao/AML_WGCNA")
getwd()
options(stringsAsFactors = FALSE)

lnames = load(file = "TPM_resm_dataInput_avelog.RData");
lnames
lnames = load(file = "TPM_avelog-auto20k.RData");
lnames
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


### plot heatmap of gene correlation ##
# *********** This can visualize the TOM correlation of genes *******
# *********** but doesn't affect the analysis *********************
# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
# ****** TOM will take about 3.5 hours *********
# ****** TOM is about 3G big ************
TOM <- TOMsimilarityFromExpr(datExpr, power = 7)
#save(TOM, file = "TPM_flterByExpr_TOM.RData")

dissTOM = 1-TOM
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function

#### ********************************
#sizeGrWindow(9,9)
#TOMplot(plotTOM, net$dendrograms, moduleColors, main = "Network heatmap plot, all genes")
####*********************************

####### plot only 400 gene
nSelect = 400
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
# https://www.biostars.org/p/394615/
TOMplot(1-plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")



#####*************Eigengene dendrogram **********
# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Isolate LIC from the clinical traits
LIC = as.data.frame(datTraits$LIC.frequency.absolute);
names(LIC) = "LIC"
# Isolate LSC17 from the clinical traits
LSC17 = as.data.frame(datTraits$LSC17);
names(LSC17) = "LSC17"

# Add the LSC17 to existing module eigengenes
MET = orderMEs(cbind(MEs,LIC, LSC17))

# Plot the relationships among the eigengenes and the trait
sizeGrWindow(7,9);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(5,6,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)

