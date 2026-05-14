#=======================================================

#=======================================================

rm(list=ls())

library(GEOquery)
library(Biobase)
library(limma)
library(dplyr)
library(tidyr)
library(gdata)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)

setwd("/proj/c.zihao/work1/1survival/CRC/GSE28722/data/")

#=======================================================

#=======================================================

gse <- getGEO("GSE28722", destdir = ".")
gpl <- getGEO('GPL13425', destdir = ".")


colnames(Table(gpl))
Table(gpl)[1:10, 1:6]

gpl <- gpl@dataTable@table
colnames(gpl)
gpl <- gpl %>% dplyr::select(ID, "ORF")
write.csv(gpl, file = 'gpl.csv')

gpl <- read.csv('gpl.csv', header = TRUE, row.names = 1)
head(gpl)
gpl$gene <- gdata::trim(gpl$ORF)

gpl$ORF <- NULL

#=======================================================

#=======================================================

exprSet <- as.data.frame(exprs(gse$GSE28722_series_matrix.txt.gz))
str(exprSet)

exprSet$ID = rownames(exprSet)
express = merge( x=gpl, y=exprSet, by="ID")
express$ID = NULL

express[which(is.na(express),arr.ind = T)]<-0 
express[1:5,1:5]

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

exprSet <- exprSet + 2
exprSet <- log2(exprSet+1)

min(exprSet)
max(exprSet)

#=======================================================

#=======================================================


pd <- pData(gse$GSE28722_series_matrix.txt.gz)

head(pd)
names(pd)
pd <- subset(pd, select=c("geo_accession", "age at diagnosis:ch2" ,"Stage:ch2", 
                          "overall survival censor (1-censored,0-non-censored):ch2",
                          "overall survival (years):ch2" ))

head(pd)
names(pd) <- c("Sample", "Age", "Stage", "OS", "OS_Time")
head(pd)

str(pd)


pd$Age  <-  as.numeric(as.character(pd$Age))
pd$OS  <-  as.numeric(as.character(pd$OS))
pd$OS_Time <-  as.numeric(as.character(pd$OS_Time))
pd$OS_Time <- round(pd$OS_Time, 2)

pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]
pd <- subset(pd, pd$OS_Time > 0)
str(pd)


str(pd)
mean(pd$Age)
table(pd$OS)
mean(pd$OS_Time)


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

