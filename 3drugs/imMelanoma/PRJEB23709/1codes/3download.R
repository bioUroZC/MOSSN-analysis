
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

setwd("/proj/c.zihao/work1/3drugs/immune/PRJEB23709/data/")


exprSet1 <- read.csv("exprSet1.csv", header = T, row.names = 1)
exprSet1[1:5,1:5]
exprSet1$Sample <- rownames(exprSet1)

pd1 <- read.csv('pd1.csv', header = T, row.names = 1)
head(pd1)


exprSet2 <- read.csv("exprSet2.csv", header = T, row.names = 1)
exprSet2[1:5,1:5]
exprSet2$Sample <- rownames(exprSet2)


pd2 <- read.csv('pd2.csv', header = T, row.names = 1)
head(pd2)


exprSet <- merge(exprSet1, exprSet2, by="Sample")
pd <- rbind(pd1, pd2)


rownames(exprSet) <- exprSet$Sample
exprSet$Sample <- NULL

exprSet[1:5,1:5]

min(exprSet)
max(exprSet)


table(pd$Response)
pd$Response[pd$Response=="Res"] <- 3
pd$Response[pd$Response=="Non"] <- 2

pd$Response <- as.numeric(as.character(pd$Response))
table(pd$Response)
pd <- na.omit(pd)

#=======================================================

#=======================================================

samplesname <- intersect(pd$Sample, colnames(exprSet))
samplesname <- unique(samplesname)
pd <- pd[which(pd$Sample %in% samplesname),]
exprSet <- exprSet[,which(colnames(exprSet) %in% samplesname)]

colnames(exprSet)
pd$Sample

#=======================================================

#=======================================================

exprSet[1:5,1:5]
write.csv(exprSet, file = "exprSet.csv")
write.csv(pd, file = "pd.csv")

