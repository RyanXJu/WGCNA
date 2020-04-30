library(dplyr)
options(stringsAsFactors = FALSE)

# directory <- "~/AML_WGCNA/WGCNA/test"
directory <- "~/AML_WGCNA/test_pipeline"
dir.create(file.path(directory), showWarnings = FALSE)
setwd(file.path(directory))

getwd()

#########################
## expression data
#########################

LSC_DE_genes = read.delim("~/AML_PCA/LSC_DE_genes.csv", header = TRUE, row.names = 1)
dim(LSC_DE_genes)
colnames(LSC_DE_genes)
rownames(LSC_DE_genes)
datExpr0 <- LSC_DE_genes[,16:48]
dim(datExpr0)

datTraits0 <- read.delim("~/AML_LSC_DE/LSC_condition.csv", header = TRUE, row.names = 1)
datTraits0$condition = replace(datTraits0$condition, datTraits0$condition == "Low", 0 )
datTraits0$condition = replace(datTraits0$condition, datTraits0$condition == "High", 1 )
datTraits0$condition = as.numeric(datTraits0$condition)
datTraits0 = datTraits0[1]
dim(datTraits0)
colnames(datTraits0) <- "LSC.Freq"

#########################
## clinical trait
#########################
traitData = read.csv("/u/juxiao/Datasets/leucegene_data/IRIC - Leucegene Database.csv");
dim(traitData)
names(traitData)

traitData$sample_id <- paste("X",traitData$sample_id,sep = "")
allTraits <- traitData[traitData$sample_id %in% colnames(datExpr0),]

names(allTraits)
dim(allTraits)

datTraits0$LIC.abs = allTraits$LIC.frequency.absolute[match(rownames(datTraits0), allTraits$sample_id)]

print(version)

write.table(datExpr0, quote = FALSE, file = "datExpr.tsv", col.names=NA, sep="\t")
write.table(datTraits0, quote = FALSE, file = "datTraits.tsv", col.names = NA, sep="\t")
