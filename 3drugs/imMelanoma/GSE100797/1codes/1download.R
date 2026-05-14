#=======================================================

#=======================================================

rm(list = ls())

library(stringr)
library(GEOquery)
library(dplyr)
library(tidyr)
library(Biobase)
library(limma)
library(data.table)
library(tibble)
library(ggplot2)
library(biomaRt)
library(gdata)
set.seed(1234)
library(org.Hs.eg.db)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)

setwd("/proj/c.zihao/work1/3drugs/immune/GSE100797/data/")

#=======================================================

#=======================================================


express <- fread("GSE100797_ProcessedData.txt")
express <- as.data.frame(express)
express[1:4,1:4]
names(express)
names(express)[1] <- 'gene'
str(express)

exprSet <- aggregate(x = express[,2:ncol(express)],
                     by = list(express$gene),
                     FUN = max)
head(exprSet)[1:5,1:5]

names(exprSet)[1] <- 'ID'
rownames(exprSet) <- exprSet$ID
exprSet$ID <- NULL
exprSet[1:5,1:5]

min(exprSet)
max(exprSet)

names(exprSet)


#=======================================================

#=======================================================

gsename <- "GSE100797"
gse <- getGEO(gsename, destdir = ".")
pd <- pData(gse$GSE100797_series_matrix.txt.gz)
names(pd)
head(pd)

pd <- subset(pd, select=c('title',
                          'ajcc.stage:ch1',
                          'recist:ch1',
                          'pfs.event:ch1', 
                          'pfs.time:ch1'))


names(pd) <- c("Sample", "Stage", "Response", "PFS", "PFS_Time")
str(pd)


table(pd$PFS)
pd$PFS  <-  as.numeric(as.character(pd$PFS))
table(pd$PFS)
str(pd)

table(pd$PFS_Time)
pd$PFS_Time  <-  as.numeric(as.character(pd$PFS_Time))
pd$PFS_Time <- pd$PFS_Time / 12
str(pd)

max(pd$PFS_Time)

pd$Sample
colnames(exprSet)


pd$Sample <- gsub('_', "", pd$Sample)
colnames(exprSet) <- gsub('_', "", colnames(exprSet))

pd$Sample
colnames(exprSet)

pd <- subset(pd, select=c("Sample", 'Response'))
pd <- na.omit(pd)

print(table(pd$Response))


pd$Response[pd$Response=="CR"] <- 4
pd$Response[pd$Response=="PR"] <- 3
pd$Response[pd$Response=="PD"] <- 1
pd$Response[pd$Response=="SD"] <- 2
print(table(pd$Response))

pd$Response <- as.numeric(as.character(pd$Response))

#=======================================================

#=======================================================

samplesname <- intersect(pd$Sample, colnames(exprSet))
samplesname <- unique(samplesname)
pd <- pd[which(pd$Sample %in% samplesname),]
exprSet <- exprSet[,which(colnames(exprSet) %in% samplesname)]
colnames(exprSet)

print(dim(exprSet))

#=======================================================

#=======================================================

print(exprSet[1:5,1:5])
write.csv(exprSet, file = "exprSet.csv")
write.csv(pd, file = "pd.csv")

