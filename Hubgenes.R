################ 6. Find hub genes in the module of interest #############

library(WGCNA)
library(ggplot2)

setwd("/u/juxiao/AML_WGCNA")
getwd()
options(stringsAsFactors = FALSE)


lnames = load(file = "TPM_resm_dataInput.RData");
lnames
lnames = load(file = "TPM_filterByExpr-auto.RData");
lnames
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# find the top hub gene in each module
chooseTopHubInEachModule(datExpr=datExpr, colorh = moduleColors, power = 2)


############ kIM and kME #############
# https://support.bioconductor.org/p/69858/

################# hub genes based on kIM ###########
### kIM: intramodular connectivity ##########
### somme of connectivity with all the other genes ####
module <- "green"
LSC17_ensembl = geneInfo[geneInfo$geneSymbol %in% LSC17_genes, 1]
geneInModule <- names(datExpr)[moduleColors == module]

adj_inModule = adjacency(datExpr[, moduleColors == module], 
                         power = 2, 
                         type = "unsigned")

# the most conneted gene
hub = which.max(rowSums(adj_inModule))
hub
# define numbers of hub wanted
numberOfHubs <- 30
hubs = names(sort(rowSums(adj_inModule),decreasing = TRUE)[1:numberOfHubs])
hubs_Symbol_kIM <- geneInfo[hubs, "geneSymbol"]


# heat map of the connectivity between genes in the module
kIM <- adj_inModule[hubs, hubs]

heatmap(kIM)

hubs_ensembl <-geneInfo[(geneInfo$gene_id %in% hubs) & geneInfo$geneSymbol !="" & !is.na(geneInfo$geneSymbol) ,
                       "gene_id"]
hubs_Symbol <- geneInfo[(geneInfo$gene_id %in% hubs) & geneInfo$geneSymbol !="" & !is.na(geneInfo$geneSymbol) ,
                        "geneSymbol"]
hubs_plot <- adj_inModule[hubs_ensembl, hubs_ensembl]
rownames(hubs_plot) <- hubs_Symbol
colnames(hubs_plot) <- hubs_Symbol
dim(hubs_plot)
heatmap(hubs_plot)


###### kME based hub genes ##########
## kME: correlation between gene and eigengene

kME_green <- signedKME(datExpr = datExpr[, gene_id_green ], datME = MEs)[,"kMEgreen"]
names(kME_green) <- gene_id_green

kME_green_hub <-sort(kME_green, decreasing = TRUE)[1:30]

kME_green_hub_ensembl <- names(kME_green_hub)
hubs_Symbol_kME <- geneInfo[kME_green_hub_ensembl, "geneSymbol"]



### difference between kIM and kME
hubs_Symbol_kIM
hubs_Symbol_kME

setdiff(hubs_Symbol_kIM, hubs_Symbol_kME)
setdiff(hubs_Symbol_kME, hubs_Symbol_kIM)

### plot ###########
library(igraph)

thr_plot <- 0.5

adj_plot <- hubs_plot
adj_plot[adj_plot<=thr_plot] <- 0
ge1=graph.adjacency(adjmatrix=adj_plot,mode="undirected",
                    weighted=TRUE,diag=FALSE)

# Compute node degrees (#links) and use that to set node size:
deg <-degree(ge1, mode="all")
V(ge1)$size <- deg*2

plot(ge1,edge.label=round(E(ge1)$weight,3), edge.arrow.size = 0.0)
plot(ge1)

################### exporting to Cytoscape ###############3333
# exporting to Cytoscape
# Select modules
modules = c("green")

# Select module probes
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule]
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);