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
setwd("/proj/c.zihao/work1/1survival/LUAD/GSE42127/data/")

#=======================================================

#=======================================================

gsename <- "GSE42127"
gse <- getGEO(gsename, destdir = ".")
gpl <- getGEO('GPL6884', destdir = ".")
colnames(Table(gpl))
Table(gpl)[1:10, 1:6]

gpl <- gpl@dataTable@table
colnames(gpl)
gpl <- gpl %>% dplyr::select(ID, "Symbol" )
write.csv(gpl, file = 'gpl.csv')

gpl <- read.csv('gpl.csv', header = TRUE, row.names = 1)
head(gpl)

gpl$Symbol <- gdata::trim(gpl$Symbol)

#=======================================================

#=======================================================

exprSet <- as.data.frame(exprs(gse$GSE42127_series_matrix.txt.gz))
str(exprSet)

exprSet$ID = rownames(exprSet)
express = merge( x=gpl, y=exprSet, by="ID")
express$ID = NULL

express[which(is.na(express),arr.ind = T)]<-0 
express[1:5,1:5]

exprSet <- aggregate(x = express[,2:ncol(express)],
                     by = list(express$Symbol),
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

pd <- pData(gse$GSE42127_series_matrix.txt.gz)
table(pd$"histology:ch1")

pd <- subset(pd, pd$"histology:ch1"  == "Adenocarcinoma")
head(pd)

pd <- subset(pd, select=c("geo_accession", "gender:ch1", "age at surgery:ch1",
                          "final.pat.stage:ch1", 
                          "survival status:ch1", "overall survival months:ch1"))

head(pd)
names(pd) <- c("Sample", 'Gender', "Age", "Stage", "OS", "OS_Time")
head(pd)
str(pd)

table(pd$Gender)
pd$Gender[pd$Gender=="F"] <- "Female"
pd$Gender[pd$Gender=="M"] <- "Male"
table(pd$Gender)

table(pd$OS)
pd$OS[pd$OS=="A"] <- 0
pd$OS[pd$OS=="D"] <- 1
table(pd$OS)
str(pd)

pd$OS  <-  as.numeric(as.character(pd$OS))
pd$OS_Time <-  as.numeric(as.character(pd$OS_Time))
pd$Age  <-  as.numeric(as.character(pd$Age))

str(pd)

pd$OS_Time <- pd$OS_Time/12

str(pd)

pd$Age <- round(pd$Age, 0)
pd$OS_Time <- round(pd$OS_Time, 2)

pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]
pd <- subset(pd, pd$OS_Time>0)
head(pd)

str(pd)
table(pd$Gender)
mean(pd$Age)
table(pd$T_stages)
table(pd$N_stages)
table(pd$M_stages)
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
