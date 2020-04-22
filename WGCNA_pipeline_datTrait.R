## Prepare trait data (select samples and their trait data)


setwd("/u/juxiao/AML_WGCNA")
getwd()

###########************************* load in trait data **************************
traitData = read.csv("/u/juxiao/Datasets/leucegene_data/IRIC - Leucegene Database.csv", header = TRUE);
dim(traitData)
colnames(traitData)

traitData$sample_id <- paste("X",traitData$sample_id,sep = "")

allTraits <- traitData[traitData$Group == "AML",]  # only AML samples have expression data
colnames(allTraits)
dim(allTraits)

# datTraits0 <- allTraits[, c("Sex","tissue")]
# dim(datTraits0)

# choose original traits
datTraits0 <- allTraits[, c("Sex", "tissue", "WHO.2008", "dx_FAB", "cytogenetic.subgroup")]

# encode traits
datTraits0$M3 <- as.integer(as.logical(datTraits0$dx_FAB == "AML-M3"))
datTraits0$MLL <- as.integer(as.logical(datTraits0$cytogenetic.subgroup == "MLL translocations (+MLL FISH positive) (Irrespective of additional cytogenetic abnormalities)"))
datTraits0$t8_21 <- as.integer(as.logical(datTraits0$WHO.2008 == "AML with t(8;21)(q22;q22); RUNX1-RUNX1T1"))

datTraits0$Sex[datTraits0$Sex == "F"] <- 0
datTraits0$Sex[datTraits0$Sex == "M"] <- 1
datTraits0$Sex <- as.numeric(datTraits0$Sex)

datTraits0$tissue[datTraits0$tissue == "Blood" ] <- 0
datTraits0$tissue[datTraits0$tissue == "Bone marrow" ] <- 1
datTraits0$tissue <- as.numeric(datTraits0$tissue)

# keep transformed traits only
datTraits0 = datTraits0[, c("Sex", "tissue", "M3", "MLL", "t8_21")]
rownames(datTraits0) <- allTraits$sample_id
dim(datTraits0)
