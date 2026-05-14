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

setwd("/proj/c.zihao/work1/1survival/OV/GSE140082/data/")
set.seed(1234)

#=======================================================

#=======================================================

gse<- getGEO("GSE140082", destdir = ".") 
gpl<- getGEO('GPL14951', destdir = ".") 

colnames(Table(gpl))
gpl <- gpl@dataTable@table
colnames(gpl)
gpl <- gpl %>% dplyr::select(ID, "Symbol")
write.csv(gpl,"GPL.csv", row.names = F)
genename = read.csv("GPL.csv")

genename$Symbol <- gdata::trim(genename$Symbol)

#=======================================================

#=======================================================


exprSet <- as.data.frame(exprs(gse$GSE140082_series_matrix.txt.gz)) 
exprSet$ID = rownames(exprSet)
express = merge( x=genename, y=exprSet, by="ID")
express$ID = NULL
express[which(is.na(express),arr.ind = T)]<-0 
express[1:4,1:4]

exprSet <- aggregate(x = express[,2:ncol(express)],
                     by = list(express$Symbol),
                     FUN = max)
head(exprSet)
exprSet <- as.data.frame(exprSet)
exprSet[1:4,1:4]

exprSet <-exprSet[-1,]
names(exprSet)[1] <- 'ID'
rownames(exprSet) <- exprSet$ID
exprSet$ID <- NULL
exprSet[1:4,1:4]

min(exprSet)
max(exprSet)

#=======================================================

#=======================================================

pd <- pData(gse$GSE140082_series_matrix.txt.gz)
head(pd)
names(pd)

table(pd$'tissue:ch1')

pd <- subset(pd, select=c( "geo_accession", 
                          "age:ch1"  ,
                          'figo_stage:ch1',
                          "final_osid:ch1"  ,
                          "final_ostm:ch1" 
))

head(pd)
colnames(pd) <- c("Sample", "Age", "Stage",  "OS", "OS_Time")
head(pd)

str(pd)

table(pd$Age)
pd$Age <- as.numeric(as.character(pd$Age))
table(pd$Age)


table(pd$OS)
pd$OS[pd$OS=="1"] <- 1
pd$OS[pd$OS=="0"] <- 0
pd$OS <-  as.numeric(as.character(pd$OS))
table(pd$OS)


str(pd)
table(pd$OS_Time)
pd$OS_Time <-  as.numeric(as.character(pd$OS_Time))
pd$OS_Time <- pd$OS_Time / 365
pd$OS_Time <- round(pd$OS_Time, 2)
str(pd)

pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]
pd <- subset(pd, pd$OS_Time>0)
head(pd)

str(pd)
table(pd$OS)
mean(pd$OS_Time)

head(pd)

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
