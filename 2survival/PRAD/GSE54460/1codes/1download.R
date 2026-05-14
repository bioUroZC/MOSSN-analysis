#=======================================================

#=======================================================


rm(list=ls())
library(dplyr)
library(tibble)
library(tidyr)
library(GEOquery)
library(gdata)
library(data.table)

options(stringsAsFactors = FALSE)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)
setwd("/proj/c.zihao/work1/1survival/PRAD/GSE54460/data")

#=======================================================

#=======================================================

gsename = "GSE54460"
gse<- getGEO(gsename, destdir = ".") 

express <- fread("GSE54460_FPKM-genes-TopHat2-106samples-12-4-13.txt")
express <- express[-c(2:16),-1] 
str(express)


express <- as.data.table(express)
colnames(express) <- as.character(express[1, ])
express <- express[-1, ]
express[1:5, 1:5, with = FALSE]
names(express)[1] <- 'gene'
express <- as.data.frame(express)
express[, 2:ncol(express)] <- lapply(express[, 2:ncol(express)], as.numeric)
head(express[, 1:5])


express[which(is.na(express),arr.ind = T)]<-0 
express[1:5,1:5]
str(express)

exprSet <- aggregate(x = express[,2:ncol(express)],
                     by = list(express$gene),
                     FUN = max)
head(exprSet)[1:5,1:5]

exprSet <- as.data.frame(exprSet)
exprSet <-exprSet[-1,]
names(exprSet)[1] <- 'ID'
rownames(exprSet) <- exprSet$ID
exprSet$ID <- NULL
exprSet[1:5,1:5]

min(exprSet)
max(exprSet)

#=======================================================

#=======================================================

pd <- gse[["GSE54460_series_matrix.txt.gz"]]@phenoData@data
names(pd)
head(pd)

pd <- subset(pd, select=c("characteristics_ch1.1",  'pstage:ch1',
                          "bcr:ch1",   "months to bcr:ch1", 'months total f/u:ch1' ))

names(pd) <- c("Sample", "Stage", "OS",  "Time1", "Time2")

head(pd)
table(pd$Time1)
pd$Sample <- substr(x=pd$Sample, start = 10, stop = 10000)


str(pd$Time1)
table(pd$Time1, useNA = "ifany")
pd$Time1[pd$Time1 == "NA"] <- pd$Time2[pd$Time1 == "NA"]

pd$OS_Time <- pd$Time1
pd$Time1 <- NULL
pd$Time2 <- NULL

str(pd)

table(pd$OS)
pd$OS <- as.numeric(as.character(pd$OS))
table(pd$OS)

table(pd$OS_Time)
pd$OS_Time <- as.numeric(as.character(pd$OS_Time))
pd$OS_Time <- pd$OS_Time / 12
pd$OS_Time <- round(pd$OS_Time, 2)
table(pd$OS_Time)


pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]
pd <- subset(pd, pd$OS_Time > 0)

table(pd$OS)
mean(pd$OS_Time)

pd$Sample
pd$Sample <- gsub("[.-]", "", pd$Sample)
pd$Sample

colnames(exprSet)
colnames(exprSet) <- gsub("[.-]", "", colnames(exprSet))
colnames(exprSet)

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


