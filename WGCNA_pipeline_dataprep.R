library(ggplot2)
library(reshape2)
library(WGCNA)
options(stringsAsFactors = FALSE)

setwd("/u/juxiao/AML_WGCNA")
getwd()


# ********************** load gene expression data ********************************
# load in TPM or FPKM or RPKM file
exprFile = "/u/juxiao/Datasets/leucegene_data/genes_TPM.unstranded.annotated.tsv"
datExpr = read.delim(exprFile, header = TRUE, row.names = 1)
datExpr[1,]
dim(datExpr)

datAnn = datExpr[,(dim(datExpr)[2]-3):(dim(datExpr)[2])]
colnames(datAnn)=c("EnsemblID", "Gene", "Category", "Location")
datAnn[1,]

datExpr = datExpr[,2:(dim(datExpr)[2]-4)]
colnames(datExpr)
rownames(datExpr)

#************************* load in trait data **************************
traitData = read.csv("/u/juxiao/Datasets/leucegene_data/IRIC - Leucegene Database.csv", header = TRUE);
dim(traitData)
colnames(traitData)

traitData$sample_id <- paste("X",traitData$sample_id,sep = "")

allTraits <- traitData[match(colnames(datExpr),traitData$sample_id),]

names(allTraits)
dim(allTraits)

datTraits <- allTraits[, c("Sex","tissue")]
dim(datTraits)

##************************** remove low read count genes **************
countFile = "/u/juxiao/Datasets/leucegene_data/star_genes_readcount.unstranded.annotated.tsv"
countData = read.delim(countFile, header = TRUE, row.names=1)
countData = countData[-c(1:4), 1:(dim(countData)[2]-4)]
rownames(countData)

# ******* WGCNA tutorial ********
# We suggest removing features whose counts are consistently low 
# (for example, removing all features that have a count of 
# less than say 10 in more than 90% of the samples) 


RemoveLowCountGene <- function(y , method = NULL, 
                       min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.7,
                       prior.count = 2, avelogCPM.cutoff = 0)
{
  # y: numeric matrix containing counts. Rows for genes and columns for libraries.
  # for the detials of other arguments, please see filterByExpr() and avelogCPM() functions in R package "edgeR" 
  if (method == "filterByCounts")
  {
    
  }
  
  if (method == "filterByExpr")
    {
    cat("Use edgeR function: filterByExpr()")
    library(edgeR)
    # create Digital Gene Expression Data
    dge <- DGEList(counts=y)
    genes_filterByExpr <- filterByExpr(dge, 
                                       min.count = min.count,
                                       min.total.count = min.total.count,
                                       large.n = large.n,
                                       min.prop = min.prop)
    # table(genes_filterByExpr)
    keep <- names(genes_filterByExpr[genes_filterByExpr])
  }
  
  if (method == "avelogCPM")
    {
    cat("Use edgeR function: avelogCPM()")
    library(edgeR)
    avelog <- aveLogCPM(y)
    keep <- rownames(y)[avelog >= avelogCPM.cutoff]
  }
  
  cat(paste0("\nRemoved ", dim(y)[1]-length(keep)," low-count genes, ", length(keep), " remaining.", sep=""))
  
  return(keep)
}


keep <- RemoveLowCountGene(countData, method = "avelogCPM")
#keep <- RemoveLowCountGene(countData, method = "filterByExpr")

# density plot
counts <- countData[keep,]
log_counts <- log2(counts + 1)

x = melt(as.matrix(log_counts))
colnames(x) = c('gene_id', 'sample', 'value')
ggplot(x, aes(x=value, color=sample))+
  geom_density()+
  theme(legend.position = "none")


### Prepare epression data for WGCNA
datExpr0 <- datExpr[keep,]
dim(datExpr0)

dim(datTraits)


############################################################ old code ###################################################

### *********** Remove all rows with less than n counts across all samples, where n=#samples ***********
# min_counts <- ncol(countData)*5
# low_count_mask <- rowSums(countData) < min_counts
# sprintf("Removing %d low-count genes (%d remaining).", sum(low_count_mask), 
#         sum(!low_count_mask))
# # [1] "Removing 24449 low-count genes (36160 remaining)."
# counts_cleaned <- countData[!low_count_mask,]
# 
# log_counts <- log2(counts_cleaned + 1)
# 
# x = melt(as.matrix(log_counts))
# colnames(x) = c('gene_id', 'sample', 'value')
# ggplot(x, aes(x=value, color=sample))+
#   geom_density()+
#   theme(legend.position = "none")


# ### *********** Use edgeR ************* 
# library(edgeR)
# 
# # create Digital Gene Expression Data
# dge <- DGEList(counts=countData)
# 
# ### ****** Use fliterByExpr() ******
# # filterByExpr(y, design = NULL, group = NULL, lib.size = NULL,
# # min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.7, ...)
# genes_filterByExpr <- filterByExpr(dge)
# table(genes_filterByExpr)
# # genes_filterByExpr
# # FALSE  TRUE 
# # 41491 19118
# keep_filterByExpr <- genes_filterByExpr[genes_filterByExpr==TRUE]
# 
# # density plot
# counts_filterByExpr <- countData[genes_filterByExpr,]
# log_counts <- log2(counts_filterByExpr + 1)
# 
# x = melt(as.matrix(log_counts))
# colnames(x) = c('gene_id', 'sample', 'value')
# ggplot(x, aes(x=value, color=sample))+
#   geom_density()+
#   theme(legend.position = "none")
# 
# 
# 
# ###****** Use aveLogCPM() ******
# # For samples that mostly have the same set of genes expressed, 
# # this is probably the best approach.
# # https://support.bioconductor.org/p/86215/
# 
# #https://support.bioconductor.org/p/83110/
# #"Genes were filtered from the analysis if their average log2 count
# #per million (as computed by edgeR's aveLogCPM function) was negative. 
# #This had the effect of keeping genes with an average count 
# #of about 5 or more per sample."
# 
# avelog <- aveLogCPM(countData)
# names(avelog) <- rownames(countData)
# keep_avelog <- names(avelog)[avelog >= 0]
# length(keep_avelog)
# # [1] 16006
# 
# # density plot
# counts_avelog <- countData[keep_avelog,]
# log_counts <- log2(counts_avelog + 1)
# 
# x = melt(as.matrix(log_counts))
# colnames(x) = c('gene_id', 'sample', 'value')
# ggplot(x, aes(x=value, color=sample))+
#   geom_density()+
#   theme(legend.position = "none")
