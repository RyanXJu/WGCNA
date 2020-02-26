####### 1st step in WGCNA ########
####### remove genes with low count reads #######

setwd("/u/juxiao/AML_WGCNA")
getwd()

fn_read_counts = "/u/juxiao/Datasets/leucegene_data/star_genes_readcount.unstranded.annotated.tsv"
read_counts = read.delim(fn_read_counts)
read_counts = read_counts[-c(1:4), 1:(length(read_counts)-4)]
rownames(read_counts) <- read_counts$X
# dim(read_counts)
# [1] 60609   692


# ******* in WGCNA tutorial ********
# We suggest removing features whose counts are consistently low 
# (for example, removing all features that have a count of 
# less than say 10 in more than 90% of the samples) 
# for (g in read_counts$X) {
#   print(g)
#   read_counts[g, "keep"] <-(sum(read_counts[g,1:691]<5) > 691*0.9 )
# }
# 
# sum(read_counts$keep)
# genes <- read_counts$X[read_counts$keep == TRUE]
# genes


### *********** Use edgeR ************* 
library(edgeR)
counts <- as.matrix(read_counts[,2:692])
rownames(counts) <- read_counts$X
# create Digital Gene Expression Data
dge <- DGEList(counts=counts)

### ****** Use fliterByExpr() ******
# filterByExpr(y, design = NULL, group = NULL, lib.size = NULL,
# min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.7, ...)
keep_filterByExpr <- filterByExpr(dge)
table(keep_filterByExpr)
# keep_filterByExpr
# FALSE  TRUE 
# 41491 19118
keep_filterByExpr <- keep_filterByExpr[keep_filterByExpr==TRUE]


###****** Use aveLogCPM() ******
# For samples that mostly have the same set of genes expressed, 
# this is probably the best approach.
# https://support.bioconductor.org/p/86215/

#https://support.bioconductor.org/p/83110/
#"Genes were filtered from the analysis if their average log2 count
#per million (as computed by edgeR's aveLogCPM function) was negative. 
#This had the effect of keeping genes with an average count 
#of about 5 or more per sample."

avelog <- aveLogCPM(counts)
names(avelog) <- rownames(counts)
keep_avelog <- names(avelog)[avelog >= 0]
length(keep_avelog)
# [1] 16006

save(keep_filterByExpr, keep_avelog, file = "/u/juxiao/AML_WGCNA/genes_kept.RData")
