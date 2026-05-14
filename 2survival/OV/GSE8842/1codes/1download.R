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

setwd("/proj/c.zihao/work1/1survival/OV/GSE8842/data")
set.seed(1234)

#=======================================================

#=======================================================


gsename <- "GSE8842"
gse <- getGEO(gsename, destdir = ".")
gpl <- getGEO('GPL5689', destdir = ".")
colnames(Table(gpl))
Table(gpl)[1:10, 1:6]

gpl <- gpl@dataTable@table
gpl <- subset(gpl, select=c("ID", "GB_ACC"))
gpl$ID <- paste0('P', gpl$ID)


anno=read.table("refGene.txt", sep="\t")
anno=anno[, c(2, 13)]
names(anno)=c("GB_ACC", "gene")
head(anno)
dim(anno)
anno <- anno %>%
  dplyr::distinct(GB_ACC, .keep_all = T)
dim(anno)

gpl <- merge(gpl, anno, by="GB_ACC")
gpl$GB_ACC <- NULL
head(gpl)


#=======================================================

#=======================================================

exprSet <- as.data.frame(exprs(gse$GSE8842_series_matrix.txt.gz))
str(exprSet)

exprSet$ID = rownames(exprSet)
exprSet$ID <- paste0("P", exprSet$ID)
express = merge( x=gpl, y=exprSet, by="ID")
express$ID = NULL

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

exprSet <- exprSet + 2

min(exprSet)
max(exprSet)


#=======================================================

#=======================================================

pd <- pData(gse$GSE8842_series_matrix.txt.gz)
head(pd)
names(pd)

table(pd$'tissue:ch1')

pd <- subset(pd, select=c( "geo_accession", 
                           'Age(years):ch1',
                           'FIGO Stage:ch1',
                           "Status:ch1"  ,
                           "Overall Survival(days):ch1" 
))

head(pd)
colnames(pd) <- c("Sample", "Age",  "Stage",  "OS", "OS_Time")
head(pd)

str(pd)



table(pd$Age)
table(pd$Age)
pd$Age <- as.numeric(as.character(pd$Age))
table(pd$Age)



table(pd$OS)
pd$OS[pd$OS=="death unrelated to cancer"] <- 1
pd$OS[pd$OS=="death related to cancer"] <- 1
pd$OS[pd$OS=="alive no evidence of disease"] <- 0
pd$OS[pd$OS=="in progression"] <- 0
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
