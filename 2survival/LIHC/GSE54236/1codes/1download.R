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
library(ggplot2)
library(biomaRt)
library(RColorBrewer)
library(gdata)
set.seed(1234)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)

setwd("/proj/c.zihao/work1/1survival/LIHC/GSE54236/data/")

#=======================================================

#=======================================================

gse <- getGEO("GSE54236", destdir = ".")
gpl <- getGEO("GPL6480", destdir = ".")

colnames(Table(gpl))
Table(gpl)[1:10, 1:6]

gpl <- gpl@dataTable@table
colnames(gpl)
gpl <- gpl %>% dplyr::select(ID,  "GENE_SYMBOL" )
write.csv(gpl, file = 'gpl.csv')

genename <- read.csv('gpl.csv', header = TRUE, row.names = 1)
head(genename)

genename$gene <- gdata::trim(genename$GENE_SYMBOL)
genename$GENE_SYMBOL <- NULL

#=======================================================

#=======================================================


exprSet <- as.data.frame(exprs(gse$GSE54236_series_matrix.txt.gz)) 
exprSet$ID = rownames(exprSet)
express = merge( x=genename, y=exprSet, by="ID")
express$ID = NULL
express[which(is.na(express),arr.ind = T)]<-0 
express[1:4,1:4]

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

exprSet <- exprSet + 3

min(exprSet)
max(exprSet)

#=======================================================

#=======================================================

pd <- pData(gse$GSE54236_series_matrix.txt.gz)
head(pd)

names(pd)
head(pd)

pd <- subset(pd, select=c("geo_accession",
                          "tissue type:ch1",
                          "gender:ch1"  ,
                          "survival time(months):ch1" ))

head(pd)
names(pd) <- c("Sample", 'type', "Gender", "OS_Time")

table(pd$type)


pd <- subset(pd, pd$type== "Biopsy of tumor tissue")
pd$type <- NULL

table(pd$Gender)
pd$Gender[pd$Gender=="male"] <- "Male"
pd$Gender[pd$Gender=="female"] <- "Female"
table(pd$Gender)

str(pd)

pd$OS <- 1
pd$OS_Time <- as.numeric(as.character(pd$OS_Time))
pd$OS_Time <- pd$OS_Time / 12
pd$OS_Time <- round(pd$OS_Time, 2)
pd <- subset(pd, pd$OS_Time > 0)
str(pd)
pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]
str(pd)

table(pd$Gender)
table(pd$OS)
mean(pd$OS_Time)
names(pd)

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


