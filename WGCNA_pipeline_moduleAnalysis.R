###### WGCNA pipeline part2: Analyse the module of intrest #######
# ****** need output data from WGCNA_pipeline_moduleDetection.R *****

#----------------------------------------------------
###### input files: data used to build network (with outlier samples and no good genes remvoed) #######
## ***** these input files are automaticly created and saved by "WGCNA_pipeline_moduleDetection.R"
# datExpr: expression data  (rownames--genes, colnames--samples)
# datTrait: trait data used to build network (rownames--samples, colnames--trait)
# netInfo : Rdata contains MEs, moduleLabels, moduleColors, geneTree of the network

###### arguments:  ########
# geneID.sel : geneID type used in expression data (needed for GO and KEGG analysis) 
# module.sel : color of the module of intrest
# trait.sel : trait of intrest
############################

library(ggplot2)
library(gplots)
library(reshape2)
library(doParallel)
library(WGCNA)
library(topGO)

options(stringsAsFactors = FALSE)

directory <- "/u/juxiao/AML_WGCNA/WGCNA/pipeline"
dir.create(file.path(directory), showWarnings = FALSE)
setwd(file.path(directory))

geneID.sel <- ""  
module.sel <- ""
trait.sel <- ""

############### load data ###########################################################
datExpr <- read.delim ("WGCNA_expressionData.tsv", header=TRUE, row.names=1)
dim(datExpr)

datTraits <- read.delim ("WGCNA_traitData.tsv", header=TRUE, row.names=1)
dim(datTraits)

netInfo <- load("WGCNA_network.RData")
netInfo

## user enter : geneID type (since KEGGREST can only take symobl, this info is needed to convert geneids later)
cat("--------------- Confirm the geneID type used in the expression data : ---------------- \n ")
cat(" this pipeline can only take \n [1]\"ensembl_gene_id\" 
    \n [2]\"ensembl_gene_id_version\" 
    \n [3]\"hgnc_symbol\"")
geneID.type.list = c("ensembl_gene_id", "ensembl_gene_id_version", "hgnc_symbol")
while (!(geneID.sel %in% c(1:3))) {
  geneID.sel <- readline(prompt="Enter the number of expression data geneID type : ")
  if (!(geneID.sel %in%  c(1:3))) { cat("The geneID type is not supported, please choose 
                                         \n [1] for \"ensembl_gene_id\" 
                                         \n [2] for \"ensembl_gene_id_version\"
                                         \n [3] for \"hgnc_symbol\"")}
  geneID.type <<- geneID.type.list[as.integer(geneID.sel)]
}
geneID.type

## user enter : color of module of intrest
cat("---------------Number of genes in each module: ---------------- \n ")
print(as.data.frame(table(moduleColors)))

while (!(module.sel %in% unique(moduleColors))) {
  module.sel <- readline(prompt="Enter the color of the module: ")
  if (!(module.sel %in% unique(moduleColors))) { cat("The entered color module does'nt exist, please Enter a color")}
}

## user enter : trait of intrest
cat("---------------Traits : ---------------- \n ")
print(colnames(datTraits))

while (!(trait.sel %in% colnames(datTraits))) {
  trait.sel <- readline(prompt="Enter the trait of intrest: ")
  if (!(trait.sel %in% colnames(datTraits))) { cat("The entered trait does'nt exist, please Enter a trait")}
}

trait <- as.data.frame(datTraits[,trait.sel])
names(trait) <- trait.sel

############### Eigengene + trait dendrogram
MET = orderMEs(cbind(MEs, trait))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5)
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2),
                      marHeatmap = c(5,4,1,2), cex.lab = 0.8,
                      xLabelsAngle= 90)

png("Eigengene_trait_dendrogram.png", width=1000, height=800)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2),
                      marHeatmap = c(5,4,1,2), cex.lab = 0.8,
                      xLabelsAngle= 90)
dev.off()

pdf("Eigengene_trait_dendrogram.pdf", width=12, height=8)
par(mar = c(6, 6, 3, 3))
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2),
                      marHeatmap = c(5,4,1,2), cex.lab = 0.8,
                      xLabelsAngle= 90)
dev.off ()



############### Gene Module Membership ##############################################
cat("\n ----------------- Calculate gene ModuleMembership --------------------\n ")
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))  # correlation of the gene and the ME
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership),  nSamples))

modNames = gsub("ME","",colnames(MEs))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")


################ GS: Gene Significance ###############################################
cat("\n ----------------- Calculate gene significance for the chosen trait ----\n ")

geneTraitSignificance = as.data.frame(cor(datExpr, trait, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), dim(datTraits)[1]))

colnames(geneTraitSignificance) = paste("GS.", trait.sel, sep="")
colnames(GSPvalue) = paste("p.GS.", trait.sel, sep="")


sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[ moduleColors==module.sel , paste0("MM", module.sel, sep="")]),
                   abs(geneTraitSignificance[moduleColors==module.sel, 1]),
                   xlab = paste("Module Membership (gene module cor) in", module.sel, "module"),
                   ylab = paste("Gene significance (gene trait cor) for ", trait.sel, sep = ""),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module.sel)


png("geneModuleMemebership_vs_geneTraitSignificance.png", width=1000, height=500)
par(mar = c(6, 6, 3, 3))
verboseScatterplot(abs(geneModuleMembership[ moduleColors==module.sel , paste0("MM", module.sel, sep="")]),
                   abs(geneTraitSignificance[moduleColors==module.sel, 1]),
                   xlab = paste("Module Membership (gene module cor) in", module.sel, "module"),
                   ylab = paste("Gene significance (gene trait cor) for ", trait.sel, sep = ""),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module.sel)
dev.off()

pdf("geneModuleMemebership_vs_geneTraitSignificance.pdf", width=12, height=8)
par(mar = c(6, 6, 3, 3))
verboseScatterplot(abs(geneModuleMembership[ moduleColors==module.sel , paste0("MM", module.sel, sep="")]),
                   abs(geneTraitSignificance[moduleColors==module.sel, 1]),
                   xlab = paste("Module Membership (gene module cor) in", module.sel, "module"),
                   ylab = paste("Gene significance (gene trait cor) for ", trait.sel, sep = ""),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module.sel)
dev.off ()



# genes in the selected module
geneInModule <- colnames(datExpr)[moduleColors == module.sel]
geneInModule

# verify geneID type, build annotation dataframe
library("biomaRt")
geneAnn<- getBM(attributes=c('ensembl_gene_id_version', 
                                  'ensembl_gene_id', 
                                  'hgnc_symbol'), 
                     filters = geneID.type, 
                     values = colnames(datExpr), 
                     mart = useMart("ensembl",dataset="hsapiens_gene_ensembl"))
geneAnn


######################### GO term enrichment analysis ########################
cat("\n ---------------GO enrichment analysis of module : ---------------- \n ")
all_genes <- geneAnn$ensembl_gene_id

geneList <- factor(as.integer (geneAnn[[geneID.type]] %in% geneInModule))
geneList
names (geneList) <- all_genes

geneSelFunc <- function (x) {
  return(x == 1 )
}

# Create topGOData object (topGO can use EnsemblID or Symbol)
GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList,
              geneSelectionFun = geneSelFunc,
              nodeSize = 5,
              annot = annFUN.org, mapping = "org.Hs.eg.db",
              ID = "Ensembl")  # ID = c("symbol", "Ensembl")

# retrieve genes2GO list from the "expanded" annotation in GOdata
allGO = genesInTerm(GOdata)
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "fisher") 
# algorithm can be :  "classic", "elim", "weight01"
# statistic can be : "fisher", "ks"

tab <- GenTable(GOdata, raw.p.value = resultKS, topNodes = length(resultKS@score), numChar = 120)
print(head(tab, 10))
      
      
######################### GSEA TODO ************########################


#separate samples into two groups (High vs. Low) based on the ME value of intrest module,
#in comparison to the mean ME level of the module of all samples. 
# GSEA was then performed between the two groups
# GSEA was used to validate the results of GO and KEGG analysis of the blue module. 
# The cut-off criterion for GSEA was FDR < 0.05.

# https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html

######################### KEGG  #################################################
cat("\n ---------------KEGG analysis of module : ---------------- \n ")

# BiocManager::install("KEGGREST")
library(KEGGREST)

# all databases
listDatabases()
# keggList("organism")
# keggList("hsa")

# Pull all pathways for AT  
pathways.list <- keggList("pathway", "hsa")
head(pathways.list)

# Pull all genes for each pathway
pathway.codes <- sub("path:", "", names(pathways.list)) 
names(pathways.list) <- pathway.codes
genes.by.pathway <- sapply(pathway.codes,
                           function(pwid){
                             pw <- keggGet(pwid)
                             if (is.null(pw[[1]]$GENE)) return(NA)
                             pw2 <- pw[[1]]$GENE[c(FALSE, TRUE)] # may need to modify this to c(FALSE, TRUE) for other organisms
                             pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
                             return(pw2)
                           }
)
head(genes.by.pathway)

geneList_KEGG <- geneModuleMembership[geneInModule, paste0("MM", module.sel, sep="")]
names(geneList_KEGG) <- geneAnn$hgnc_symbol[match(geneInModule, geneAnn[[geneID.type]] )]
head(geneList_KEGG)

# Wilcoxon test for each pathway
pVals.by.pathway <- t(sapply(names(genes.by.pathway),
    function(pathway) {
        # print(pathway)
        pathway.genes <- genes.by.pathway[[pathway]]
        
        if (is.na(pathway.genes[1])){
          p.value <- NA
          list.genes.in.pathway <- 0
        }
        
        else{
          list.genes.in.pathway <- intersect(names(geneList_KEGG), pathway.genes)
          list.genes.not.in.pathway <- setdiff(names(geneList_KEGG), list.genes.in.pathway)
          scores.in.pathway <- geneList_KEGG[list.genes.in.pathway]
          scores.not.in.pathway <- geneList_KEGG[list.genes.not.in.pathway]
          
          if (length(scores.in.pathway) > 0){
            p.value <- wilcox.test(scores.in.pathway, scores.not.in.pathway, alternative = "less")$p.value
          }
          else{
            p.value <- NA
          }
        }
        
        return(c(p.value = p.value , Annotated = length(list.genes.in.pathway)))
    }
))

# Assemble output table
outdat <- data.frame(pathway.code = rownames(pVals.by.pathway))
outdat$pathway.name <- pathways.list[outdat$pathway.code]
outdat$p.value <- pVals.by.pathway[,"p.value"]
outdat$Annotated <- pVals.by.pathway[,"Annotated"]
outdat <- outdat[order(outdat$p.value),]
print(head(outdat,10))


# ################# hub genes based on kIM (intramodular connectivity)  ###########
# adj_inModule = adjacency(datExpr[, moduleColors == module.sel], 
#                          power = 2, 
#                          type = "unsigned")
# # suggest: power = 2 for unsigned, 4 for signed
# 
# # the most conneted gene
# hub = which.max(rowSums(adj_inModule))
# hub
# # define numbers of hub wanted
# numberOfHubs <- 30
# hubs = names(sort(rowSums(adj_inModule),decreasing = TRUE)[1:numberOfHubs])
# hubs_Symbol_kIM <- geneInfo[hubs, "geneSymbol"]
# 
# 
# # heat map of the connectivity between genes in the module
# kIM <- adj_inModule[hubs, hubs]
# 
# heatmap(adj_inModule)
# 
