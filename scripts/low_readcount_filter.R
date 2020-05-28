####### For the analysis, it is ########
### use edgeR package functions : filterByExpr or avelogCPM ##########

library(ggplot2)
library(reshape2)

setwd("/u/juxiao/AML_WGCNA")
getwd()

fn_read_counts = "/u/juxiao/Datasets/leucegene_data/star_genes_readcount.unstranded.annotated.tsv"
read_counts = read.delim(fn_read_counts)
read_counts = read_counts[-c(1:4), 1:(length(read_counts)-4)]
rownames(read_counts) <- read_counts$X
# dim(read_counts)
# [1] 60609   692


### density plot of read_counts
counts <- as.matrix(read_counts[,2:692])

log_counts <- log2(counts + 1)

x = melt(as.matrix(log_counts))
colnames(x) = c('gene_id', 'sample', 'value')
ggplot(x, aes(x=value, color=sample))+
  geom_density()+
  theme(legend.position = "none")

###################################################################################
###### OPTION 1 : 
#### Remove all rows with less than n counts across all samples, where n=#samples
###################################################################################
low_count_mask <- rowSums(counts) < ncol(counts)
sprintf("Removing %d low-count genes (%d remaining).", sum(low_count_mask), 
        sum(!low_count_mask))
# [1] "Removing 24449 low-count genes (36160 remaining)."
counts_cleaned <- counts[!low_count_mask,]
log_counts <- log2(counts_cleaned + 1)
x = melt(as.matrix(log_counts))
colnames(x) = c('gene_id', 'sample', 'value')
ggplot(x, aes(x=value, color=sample))+
  geom_density()+
  theme(legend.position = "none")


###################################################################################
###### OPTION 2 : 
### *********** Use edgeR ************* 
###################################################################################
library(edgeR)
counts <- as.matrix(read_counts[,2:692])
rownames(counts) <- read_counts$X
# create Digital Gene Expression Data
dge <- DGEList(counts=counts)

### ****** Use fliterByExpr() ******
# filterByExpr(y, design = NULL, group = NULL, lib.size = NULL,
# min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.7, ...)
genes_filterByExpr <- filterByExpr(dge)
table(genes_filterByExpr)
# genes_filterByExpr
# FALSE  TRUE 
# 41491 19118
keep_filterByExpr <- names(genes_filterByExpr[genes_filterByExpr==TRUE])

# density plot
counts_filterByExpr <- counts[genes_filterByExpr,]
log_counts <- log2(counts_filterByExpr + 1)

x = melt(as.matrix(log_counts))
colnames(x) = c('gene_id', 'sample', 'value')
ggplot(x, aes(x=value, color=sample))+
  geom_density()+
  theme(legend.position = "none")


###################################################################################
###### OPTION 3 : 
###****** Use aveLogCPM() ******
# For samples that mostly have the same set of genes expressed, 
# this is probably the best approach.
# https://support.bioconductor.org/p/86215/

#https://support.bioconductor.org/p/83110/
#"Genes were filtered from the analysis if their average log2 count
#per million (as computed by edgeR's aveLogCPM function) was negative. 
#This had the effect of keeping genes with an average count 
#of about 5 or more per sample."
####################################################################################
avelog <- aveLogCPM(counts)
names(avelog) <- rownames(counts)
keep_avelog <- names(avelog)[avelog >= 0]
length(keep_avelog)
# [1] 16006

# density plot
counts_avelog <- counts[avelog>0,]
log_counts <- log2(counts_avelog + 1)

x = melt(as.matrix(log_counts))
colnames(x) = c('gene_id', 'sample', 'value')
ggplot(x, aes(x=value, color=sample))+
  geom_density()+
  theme(legend.position = "none")


save(keep_filterByExpr, keep_avelog, file = "/u/juxiao/Datasets/leucegene_data/genes_kept_edgeR.RData")

