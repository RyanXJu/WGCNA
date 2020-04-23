################ 2. Find genes modules in leucegene dataset #############

# install.packages("WGCNA")
library(WGCNA)
library(ggplot2)

setwd("/u/juxiao/AML_WGCNA")
getwd()
options(stringsAsFactors = FALSE)


#########################
## gene expression data
#########################
# load rna-seq data (normalized)
fn_TPM_resm = "/u/juxiao/Datasets/leucegene_data/genes_TPM.unstranded.annotated.tsv"
TPM_resm = read.delim(fn_TPM_resm)
# dim(TPM_resm)
# [1] 60609   697
# note: there are only 691 samples in the 

# Load the filterd gene info
genes<-load(file="/u/juxiao/AML_WGCNA/genes_kept.RData")
genes

# choose kept gene info to use
gene_keep <- keep_avelog 
length(gene_keep)

sample_id <- colnames(TPM_resm)[startsWith(colnames(TPM_resm), "X")]
datExpr0 <- as.data.frame(t(TPM_resm[,sample_id]))
names(datExpr0) <- TPM_resm$gene_id
datExpr0[1:3,1:3]
# ENSG00000000003.15 ENSG00000000005.6 ENSG00000000419.12
# X01H001               0.02              0.00              35.68
# X01H002               0.05              0.00              20.44
# X02H003               0.08              0.01              32.15

# keep only flitered genes, log2(x+1)
datExpr0 <- log2(datExpr0[,gene_keep]+1)
rownames(datExpr0) <- sample_id
dim(datExpr0)

# verify expression data
gsg = goodSamplesGenes(datExpr0, verbose = 3)
# goodSamplesGenes() default parameters
# minFraction = 1/2 :minimum fraction of non-missing samples for a gene
# minNSamples =4 minimum number of non-missing samples for a gene

gsg$allOK

## remove no good genes
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
dim(datExpr0)

# ********** check outlier samples will be very long *************
sampleTree = hclust(dist(datExpr0), method = "average")
# dist(x, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)

# Plot the sample tree
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
## remove outliers
# Plot a line to show the cut
abline(h = 200, col = "red");


# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 200, minSize = 10)
table(clust)
# clust with keep_filterByExpr, cutHeight =200
# 0   1 
# 14 677 
# clust with keep_avelog, cutHeight =200
# 0   1 
# 15 676 

# clust 2 contains the samples we want to keep.
outlierSamples = (clust == 0)
outlierSamples <- datExpr0[outlierSamples,]
outlierSamples_names<- rownames(outlierSamples)
write.table(outlierSamples_names, file = "/u/juxiao/AML_WGCNA/outlierSamples_id.txt",
            sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)


keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
nGenes
nSamples


#########################
## clinical trait
#########################
traitData = read.csv("/u/juxiao/Datasets/leucegene_data/IRIC - Leucegene Database.csv");
dim(traitData)
names(traitData)

traitData$sample_id <- paste("X",traitData$sample_id,sep = "")
allTraits <- traitData[traitData$sample_id %in% rownames(datExpr0),]

names(allTraits)
dim(allTraits)

###****** ADD LSC17 and mRNAsi data ******
fn_TPM_lsc17  = "/u/juxiao/LSC17/TPM_lsc17.tsv"
TPM_lsc17  = read.delim(fn_TPM_lsc17, as.is=TRUE, check.names=FALSE, header=FALSE) 
colnames(TPM_lsc17) <- c("sample_id", "LSC17")
TPM_lsc17$sample_id <- paste("X", TPM_lsc17$sample_id, sep="")

fn_TPM_mRNAsi  = "/u/juxiao/mRNAsi/TPM_StemScore.tsv"
TPM_mRNAsi = read.delim(fn_TPM_mRNAsi, as.is=TRUE, check.names=FALSE, header=FALSE)
colnames(TPM_mRNAsi) <- c("sample_id", "mRNAsi")
TPM_mRNAsi$sample_id <- paste("X", TPM_mRNAsi$sample_id, sep="")

allTraits <- merge(allTraits, TPM_lsc17, by="sample_id")
allTraits <- merge(allTraits, TPM_mRNAsi, by="sample_id")
dim(allTraits)


# Form a data frame analogous to expression data that 
# will hold the clinical traits.
AML_Samples = rownames(datExpr)
traitRows = match(rownames(datExpr), allTraits$sample_id)


### ********* create OS feature ************
# combine overall survival info to a new feature
OS_days <- allTraits$Overall_Survival_Time_days
OS_mean <- mean((OS_days)[!is.na(OS_days)])
OS_mean
# for 452 patients, the mean survival days is 824.19 days
OS_median <- median((OS_days)[!is.na(OS_days)])
OS_median
# median survival days is 300 days

hist((OS_days)[!is.na(OS_days)], breaks = 20)
abline(v = OS_median, col = "blue", lwd = 2)
abline(v = 365, col = "red", lwd = 2)

# female
OS_mean_f = mean(OS_days[!is.na(OS_days) & allTraits$sex == "F"]) #791 days
OS_median_f = median(OS_days[!is.na(OS_days) & allTraits$sex == "F"]) #304days
hist(OS_days[!is.na(OS_days) & allTraits$sex == "F"], breaks = 20)
abline(v = OS_median, col = "blue", lwd = 2)
abline(v = 365, col = "red", lwd = 2)

# male
OS_mean_m = mean(OS_days[!is.na(OS_days) & allTraits$sex == "M"]) #850days
OS_median_m = median(OS_days[!is.na(OS_days) & allTraits$sex == "M"]) #296days
hist(OS_days[!is.na(OS_days) & allTraits$sex == "M"], breaks = 20)
abline(v = OS_median, col = "blue", lwd = 2)
abline(v = 365, col = "red", lwd = 2)

# OS_Status: 1=decease 0=alive
OS_rate <- sum(allTraits$Overall_Survival_Status[!is.na(allTraits$Overall_Survival_Status)])
# in total 338 are dead

# set threshold (if OS_threshold is big enough, we should have 338 OS= 1)
OS_threshold <- 365*1 # 3year(n=427), 2year(n=443), 1year(n=449)

# allTraits$OS: if patient pass the OS_threshold
allTraits$OS <- NA
# patients dead within the threshold time
allTraits$OS[allTraits$Overall_Survival_Status == 1 & 
               allTraits$Overall_Survival_Time_days <= OS_threshold] <- 0
# patient live longer than threshold
allTraits$OS[allTraits$Overall_Survival_Time_days > OS_threshold] <- 1
sum(!is.na(allTraits$OS))

##################### prepare datTraits #########################
# select traits of intrest
# allTraits col: 18-FAB, 22-cytogenetic 30-OS_days 34-WHO
# 72-OS : if alive after the chosen OS threshold
datTraits = allTraits[traitRows, c(30,72,44,59,70,71) ]
rownames(datTraits) = allTraits[traitRows, 1]

datTraits$Sex[datTraits$Sex == "F"] <- 0
datTraits$Sex[datTraits$Sex == "M"] <- 1
datTraits$Sex <- as.numeric(datTraits$Sex)
dim(datTraits)
datTraits[1:3,1:6]
# Overall_Survival_Time_days OS Sex LIC.frequency.absolute  LSC17     mRNAsi
# X01H001                       2169  1   0                     NA 184.55 0.23628790
# X01H002                       5161  1   1                     NA 223.79 0.39975449
# X02H003                       1396  1   0                     NA 490.80 0.06926346

# sum(is.na(datTraits$Overall_Survival_Time_days))

collectGarbage()


# # Re-cluster kept samples
sampleTree2 = hclust(dist(datExpr), method = "average")

# Convert traits to a color representation: 
# white-red: low to high, grey: missing entry
traitColors = numbers2colors(datTraits, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

save(datExpr, datTraits, file = "TPM_resm_dataInput_avelog.RData")



##################### prepare datTraits for FAB #########################
# select traits of intrest
# allTraits col: 18-FAB, 22-cytogenetic 30-OS_days 34-WHO
datTraits_FAB = allTraits[traitRows, c(30, 72,18) ]
table(datTraits_FAB$dx_FAB)
datTraits_FAB$dx_FAB[datTraits_FAB$dx_FAB == "AML-M0"] <- 0
datTraits_FAB$dx_FAB[datTraits_FAB$dx_FAB == "AML-M1"] <- 1
datTraits_FAB$dx_FAB[datTraits_FAB$dx_FAB == "AML-M2"] <- 2
datTraits_FAB$dx_FAB[datTraits_FAB$dx_FAB %in% c("AML-M3", "AML-M3V")] <- 3
datTraits_FAB$dx_FAB[datTraits_FAB$dx_FAB %in% c("AML-M4","AML-M4Eo")] <- 4
datTraits_FAB$dx_FAB[datTraits_FAB$dx_FAB %in% c("AML-M5","AML-M5A","AML-M5B")] <-5
datTraits_FAB$dx_FAB[datTraits_FAB$dx_FAB %in% c("AML-M6","AML-M6B")] <- 6
datTraits_FAB$dx_FAB[datTraits_FAB$dx_FAB == "AML-M7"] <- 7
datTraits_FAB$dx_FAB[datTraits_FAB$dx_FAB %in% c("AREB-T",
                                                 "LMA-NOS",
                                                 "LMAMD", 
                                                 "MDS-EB2",
                                                 "Not classifiable by FAB criteria")] <- 8
datTraits_FAB$dx_FAB <- as.numeric(datTraits_FAB$dx_FAB)

datTraits_FAB$M0 <- as.numeric(datTraits_FAB$dx_FAB == 0)
datTraits_FAB$M1 <- as.numeric(datTraits_FAB$dx_FAB == 1)
datTraits_FAB$M2 <- as.numeric(datTraits_FAB$dx_FAB == 2)
datTraits_FAB$M3 <- as.numeric(datTraits_FAB$dx_FAB == 3)
datTraits_FAB$M4 <- as.numeric(datTraits_FAB$dx_FAB == 4)
datTraits_FAB$M5 <- as.numeric(datTraits_FAB$dx_FAB == 5)
datTraits_FAB$M6 <- as.numeric(datTraits_FAB$dx_FAB == 6)
datTraits_FAB$M7 <- as.numeric(datTraits_FAB$dx_FAB == 7)
datTraits_FAB$M8 <- as.numeric(datTraits_FAB$dx_FAB == 8)

rownames(datTraits_FAB) = allTraits[traitRows, 1]
dim(datTraits_FAB)
datTraits_FAB[1:3,1:9]
# Overall_Survival_Time_days OS dx_FAB M0 M1 M2 M3 M4 M5
# X01H001                       2169  1      5  0  0  0  0  0  1
# X01H002                       5161  1      3  0  0  0  1  0  0
# X02H003                       1396  1      0  1  0  0  0  0  0

collectGarbage()

sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(datTraits_FAB, signed = FALSE)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits_FAB),
                    main = "Sample dendrogram and FAB heatmap")

save(datExpr, datTraits_FAB, file = "TPM_resm_dataInput_avelog_FAB.RData")


##################### prepare datTraits for WHO #########################
# select traits of intrest
# allTraits col: 18-FAB, 22-cytogenetic 30-OS_days 34-WHO
datTraits_WHO = allTraits[traitRows, c(30, 72,34) ]
table(datTraits_WHO$WHO.2008)
# put all subclass with less than 15 samples into one group
datTraits_WHO$WHO.2008[datTraits_WHO$WHO.2008 == "Acute monoblastic and monocytic leukaemia"] <- 0
datTraits_WHO$WHO.2008[datTraits_WHO$WHO.2008 == "Acute myeloid leukaemia, NOS"] <- 1
datTraits_WHO$WHO.2008[datTraits_WHO$WHO.2008 == "Acute myelomonocytic leukaemia"] <- 2
datTraits_WHO$WHO.2008[datTraits_WHO$WHO.2008 == "Acute promyelocytic leukaemia with t(15;17)(q24;q21); PML-RARA"] <- 3
datTraits_WHO$WHO.2008[datTraits_WHO$WHO.2008 == "AML with inv(16)(p13.1q22) or t(16;16)(p13.1;q22); CBFB-MYH11"] <- 4
datTraits_WHO$WHO.2008[datTraits_WHO$WHO.2008 == "AML with maturation"] <- 5
datTraits_WHO$WHO.2008[datTraits_WHO$WHO.2008 == "AML with minimal differentiation"] <- 6
datTraits_WHO$WHO.2008[datTraits_WHO$WHO.2008 == "AML with myelodysplasia-related changes"] <- 7
datTraits_WHO$WHO.2008[datTraits_WHO$WHO.2008 == "AML with t(8;21)(q22;q22); RUNX1-RUNX1T1"] <- 8
datTraits_WHO$WHO.2008[datTraits_WHO$WHO.2008 == "AML without maturation"] <- 9
datTraits_WHO$WHO.2008[datTraits_WHO$WHO.2008 == "Therapy-related myeloid neoplasms"] <- 10
# put all subclass with less than 15 samples into one group
datTraits_WHO$WHO.2008[datTraits_WHO$WHO.2008 %in% c("Acute erythroid leukaemia",
                                                 "Acute megakaryoblastic leukaemia",
                                                 "AML with inv(3)(q21.3q26.2)  or t(3;3)(q21.3;q26.2); GATA2, MECOM_WHO2016", 
                                                 "AML with inv(3)(q21q26.2) or t(3;3)(q21;q26.2); RPN1-EVI1",
                                                 "AML with myelodysplasia-related changes (WHO2017)",
                                                 "AML with t(6;9)(p23;q34); DEK-NUP214",
                                                 "AML with t(9;11)(p22;q23); MLLT3-MLL",
                                                 "Myelodysplastic syndrome with excess blasts (WHO 2016)")] <- 11
datTraits_WHO$WHO.2008 <- as.numeric(datTraits_WHO$WHO.2008)

datTraits_WHO$monob_monoc_0 <- as.numeric(datTraits_WHO$WHO.2008 == 0)
datTraits_WHO$myeloid_NOS_1 <- as.numeric(datTraits_WHO$WHO.2008 == 1)
datTraits_WHO$myelomono_2 <- as.numeric(datTraits_WHO$WHO.2008 == 2)
datTraits_WHO$PML_RARA_3 <- as.numeric(datTraits_WHO$WHO.2008 == 3)
datTraits_WHO$CBFB_MYH11_4 <- as.numeric(datTraits_WHO$WHO.2008 == 4)
datTraits_WHO$maturation_5 <- as.numeric(datTraits_WHO$WHO.2008 == 5)
datTraits_WHO$minDiff_6 <- as.numeric(datTraits_WHO$WHO.2008 == 6)
datTraits_WHO$myelodysplasia_7 <- as.numeric(datTraits_WHO$WHO.2008 == 7)
datTraits_WHO$RUNX1_RUNX1T1_8 <- as.numeric(datTraits_WHO$WHO.2008 == 8)
datTraits_WHO$noMaturation_9 <- as.numeric(datTraits_WHO$WHO.2008 == 9)
datTraits_WHO$Therapy_10<- as.numeric(datTraits_WHO$WHO.2008 == 10)
datTraits_WHO$Others_11 <- as.numeric(datTraits_WHO$WHO.2008 == 11)

rownames(datTraits_WHO) = allTraits[traitRows, 1]
dim(datTraits_WHO)
datTraits_WHO[1:3,1:9]
# Overall_Survival_Time_days OS WHO.2008 monob_monoc myeloid_NOS myelomono
# X01H001                       2169  1       10           0           0         0
# X01H002                       5161  1        3           0           0         0
# X02H003                       1396  1        6           0           0         0

collectGarbage()

sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(datTraits_WHO, signed = FALSE)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits_WHO),
                    main = "Sample dendrogram and trait heatmap")

save(datExpr, datTraits_WHO, file = "TPM_resm_dataInput_avelog_WHO.RData")
