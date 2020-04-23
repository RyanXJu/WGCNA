############ 4 find correlation between module and Traits #####

library(WGCNA)

setwd("/u/juxiao/AML_WGCNA")
getwd()

lnames = load(file = "TPM_resm_dataInput.RData");
lnames
lnames = load(file = "TPM_filterByExpr-auto.RData");
lnames
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
# change module names with their colors
MEs = orderMEs(MEs0)

# module vs traits
cor <- stats::cor
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))



#################### LIC genes-traits ################################
# Define variable LIC containing the LIC column of datTrait
LIC = as.data.frame(datTraits$LIC.frequency.absolute);
names(LIC) = "LIC"
# names (colors) of the modules
modNames = substring(names(MEs), 3)


geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), 
                                          nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

# GS: Gene Significance
geneTraitSignificance = as.data.frame(cor(datExpr, LIC, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(LIC), sep="");
names(GSPvalue) = paste("p.GS.", names(LIC), sep="")

geneTraitSignificance_LIC <- geneTraitSignificance
GSPvalue_LIC <- GSPvalue

#### module vs LIC and pvalue
moduleLIC <- cbind(moduleTraitCor[,4], moduleTraitPvalue[,4])
colnames(moduleLIC) <- c("moduleLICcor", "moduleLICpvalue")
moduleLIC
# the most significant correlated module
# MEgreen          0.402761161    8.625843e-28

module = "green"
column = match(module, modNames)
moduleGenes = moduleColors==module

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership (gene module cor) in", module, "module"),
                   ylab = "Gene significance (gene trait cor) for LIC",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

names(datExpr)
green_gene_id <-names(datExpr)[moduleColors=="green"] #540


##################### LSC17 genes-traits #############################
# Define variable LSC17 containing the LSC17 column of datTrait
LSC17 = as.data.frame(datTraits$LSC17);
names(LSC17) = "LSC17"
# names (colors) of the modules
modNames = substring(names(MEs), 3)


geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), 
                                          nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

# GS: Gene Significance
geneTraitSignificance = as.data.frame(cor(datExpr, LSC17, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(LSC17), sep="");
names(GSPvalue) = paste("p.GS.", names(LSC17), sep="")

#### module vs LSC17 and pvalue
moduleLSC17 <- cbind(moduleTraitCor[,5], moduleTraitPvalue[,5])
colnames(moduleLSC17) <- c("moduleLSC17cor", "moduleLSC17pvalue")
moduleLSC17
# the most significant correlated module
# MEgreen            0.745438041     5.155368e-121

module = "green"
column = match(module, modNames)
moduleGenes = moduleColors==module

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership (gene module cor) in", module, "module"),
                   ylab = "Gene significance (gene trait cor) for LIC",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

names(datExpr)
green_gene_id <-names(datExpr)[moduleColors=="green"] #540


# source("/u/juxiao/AML_datasets/probe2gene.R")
# green_gene <- probe2gene(names(datExpr)[ModuleColors=="green"],
#                          srcType = "ensembl_gene_id_version" ) 

############  Convert Ensembl_gene_id to gene_symbol ################
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
allLLIDs <- getBM(
  filters="ensembl_gene_id_version",
  attributes=c("ensembl_gene_id_version", "entrezgene_id", "hgnc_symbol"),
  values=names(datExpr),
  mart=mart) # 438 have symbol in green module, 91 have no symbol(symbol="")

# Create the starting data frame
geneInfo0 = data.frame(gene_id = names(datExpr),
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue,
                       geneTraitSignificance_LIC,
                       GSPvalue_LIC)

gene_id_green <- geneInfo0$gene_id[geneInfo0$moduleColor == "green"]

# Order modules by their significance for LSC17
modOrder = order(-abs(cor(MEs, datTraits$LSC17, use = "p")))

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance

gene_annot <- geneInfo0$gene_id %in% allLLIDs$ensembl_gene_id_version
geneInfo0 <- geneInfo0[allLLIDs$ensembl_gene_id_version,]
geneInfo0$geneSymbol<-allLLIDs$hgnc_symbol
geneInfo0$LocusLinkID<-allLLIDs$entrezgene_id

geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.LSC17));
geneInfo = geneInfo0[geneOrder, ]

## verify the genes in LSC17 signature
## change "GPR56" to "ADGRG1", "KIAA0125" to "FAM30A"
LSC17_genes <- c("DNMT3B","ZBTB46", "NYNRIN", "ARHGAP22", 
           "LAPTM4B", "MMRN1", "DPYSL3","FAM30A", 
           "CDK6", "CPXM1", "SOCS2", "SMIM24", 
           "EMP1", "NGFRAP1", "CD34","AKR1C3", 
           "ADGRG1")

verified <- geneInfo[geneInfo$geneSymbol %in% LSC17_genes, c(2:6,55)]
verified
# miss "NGFRAP1

geneInfo_green <- geneInfo[geneInfo$moduleColor =="green" & geneInfo$geneSymbol!="",]
write.table(geneInfo_green$geneSymbol, file = "/u/juxiao/AML_WGCNA/green_geneid.txt",
            sep = ",", quote = FALSE, row.names = FALSE)

write.csv(geneInfo, file = "geneInfo.csv")



##### correlation between cytogenetic subgroup and modules
cyto <- load("cyto_group.R")
cyto

cyto_group <- cyto_group[rownames(MEs),]
cor <- stats::cor
moduleCytoCor = cor(MEs, cyto_group, use = "p");
moduleCytoPvalue = corPvalueStudent(moduleCytoCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleCytoCor, 2), "\n(",
                    signif(moduleCytoPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleCytoCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleCytoCor,
               xLabels = colnames(cyto_group),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-cytogenetic relationships"))




# ############################ GO enrichment ################# 
# # GOenrichmentAnalysis is offred by WGCNA
# # new function called: enrichmentAnalysis 
# # in R package anRichment
# GOenr = GOenrichmentAnalysis(moduleColors, 
#                              allLLIDs$entrezgene_id, 
#                              organism = "human", 
#                              nBestP = 10)
# 
# tab = GOenr$bestPTerms[[4]]$enrichment
# names(tab)
# write.table(tab, file = "/u/juxiao/AML_WGCNA/GOEnrichmentTable_WGCNA.csv",
#             sep = ",", quote = TRUE, row.names = FALSE)


# 
# ################## A few queries #############
# geneInfo_green <- geneInfo[geneInfo$moduleColor == "green"]
# # most LIC correlated genes 
# geneInfo[order(geneInfo$GS.LIC, decreasing = TRUE)[1:30],c(2:6,61)]
# # most LIC correlated genes in green module
# geneInfo_green[order(geneInfo_green$GS.LIC, decreasing = TRUE),c(2:6,61)]
# 
# # most LSC17 correlated genes
# geneInfo[order(geneInfo$GS.LSC17, decreasing = TRUE)[1:30],c(2:6,61)]

