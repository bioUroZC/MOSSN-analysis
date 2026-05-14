#=======================================================

#=======================================================

rm(list = ls())

library(GEOquery)
library(dplyr)
library(tidyr)
library(Biobase)
library(limma)
library(data.table)
library(tibble)
library(biomaRt)
library(gdata)
library(stringr)
library(readxl)
library(org.Hs.eg.db)
library(AnnotationDbi)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)

setwd("/proj/c.zihao/work1/3drugs/Paclitaxel/GSE194040/data")
set.seed(1234)

#=======================================================

#=======================================================

gsename <- "GSE194040"
gse <- getGEO(gsename, destdir = ".")


express <- fread('GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988.txt')
express <- as.data.frame(express)
express[1:5,1:5]
names(express)[1] <- 'gene'


express[which(is.na(express),arr.ind = T)]<-0 
express[1:5,1:5]

exprSet <- aggregate(x = express[,2:ncol(express)],
                     by = list(express$gene),
                     FUN = max)
head(exprSet)[1:5,1:5]

exprSet <- as.data.frame(exprSet)
names(exprSet)[1] <- 'ID'
rownames(exprSet) <- exprSet$ID
exprSet$ID <- NULL
exprSet[1:5,1:5]

min(exprSet)
max(exprSet)

colnames(exprSet) <- paste0("Sam", colnames(exprSet))

#=======================================================

#=======================================================


clinical <- read_excel('NIHMS1829047-supplement-3.xlsx')
clinical <- as.data.frame(clinical)
clinical <- subset(clinical, select=c("Patient Identifier",
                                      "Arm", "pCR" ))
names(clinical) <- c("Sample", "Arm", "Response")
clinical$Sample <- paste0("Sam", clinical$Sample)
head(clinical)

table(clinical$Arm)

pd <- clinical[grepl("Paclitaxel", clinical$Arm), ]

print(table(pd$Response))
pd$Response[pd$Response=='0'] <- 2
pd$Response[pd$Response=='1'] <- 3
print(table(pd$Response))


pd$Response <- as.numeric(as.character(pd$Response))
print(pd$Response)


pd <- subset(pd, select=c("Sample", "Response"))
pd <- na.omit(pd)


pd$Sample[1:10]
colnames(exprSet)[1:10]


#=======================================================

#=======================================================

samplesname <- intersect(pd$Sample, colnames(exprSet))
samplesname <- unique(samplesname)
pd <- pd[which(pd$Sample %in% samplesname),]
exprSet <- exprSet[,which(colnames(exprSet) %in% samplesname)]
colnames(exprSet)

#=======================================================

#=======================================================

exprSet[1:5,1:5]
write.csv(exprSet, file = "exprSet.csv")
write.csv(pd, file = "pd.csv")
