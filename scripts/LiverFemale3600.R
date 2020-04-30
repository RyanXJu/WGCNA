
# Display the current working directory
getwd()
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "~/WGCNA"
setwd(workingDir); 

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

################# import annotation data ######################
annot = read.csv(file = "~/WGCNA/GeneAnnotation.csv")



################# prepare expression data #####################
#Read in the female liver data set
femData = read.csv("LiverFemale3600.csv");
# Take a quick look at what is in the data set:
dim(femData);
names(femData);

datExpr0 = as.data.frame(t(femData[, -c(1:8)]))
names(datExpr0) = femData$gene_symbol
rownames(datExpr0) = names(femData)[-c(1:8)]

################## prepare Trait data ##########################
traitData = read.csv("ClinicalTraits.csv")
dim(traitData)
names(traitData)

# remove columns that hold information we do not need.
allTraits = traitData[, -c(31, 16)]
allTraits = allTraits[, c(2, 11:36) ]
dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.

femaleSamples = rownames(datExpr0)
traitRows = match(femaleSamples, allTraits$Mice)
datTraits = allTraits[traitRows, -1]
rownames(datTraits) = allTraits[traitRows, 1]

write.table(datExpr0, quote=FALSE, file="datExpr_liverFemale3600.tsv", col.names=NA, sep="\t")
write.table(datTraits, quote=FALSE, file="datTraits_liverFemale3600.tsv", col.names=NA, sep="\t")


