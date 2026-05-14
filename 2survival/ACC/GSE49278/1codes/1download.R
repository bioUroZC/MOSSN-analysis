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

setwd("/proj/c.zihao/work1/1survival/ACC/GSE49278/data/")

#=======================================================

#=======================================================

gsename <- "GSE49278"
gse <- getGEO(gsename, destdir = ".")
gpl <- getGEO('GPL16686', destdir = ".")
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

exprSet <- as.data.frame(exprs(gse$GSE49278_series_matrix.txt.gz))
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

#=======================================================

#=======================================================


dt <- pData(gse$GSE49278_series_matrix.txt.gz)
dt[1:5,1:5]
dt <- subset(dt, select=c('geo_accession', 'title'))
names(dt) <- c('Sample', 'name')

pd <- read.csv('clinical.csv', header = T)
str(pd)


pd <- merge(dt, pd, by='name')
head(pd)
pd$name <- NULL

table(pd$Gender)
pd$Gender[pd$Gender=="F"] <- "Female"
pd$Gender[pd$Gender=="M"] <- "Male"
table(pd$Gender)

table(pd$Age)
pd$Age <- as.numeric(as.character(pd$Age))
table(pd$Age)

table(pd$OS)
pd$OS[pd$OS=="no"] <- 0
pd$OS[pd$OS=="yes"] <- 1
table(pd$OS)
pd$OS  <-  as.numeric(as.character(pd$OS))
str(pd)

table(pd$OS_Time)
pd$OS_Time <-  as.numeric(as.character(pd$OS_Time))
pd$OS_Time <- pd$OS_Time / 12
pd$OS_Time <- round(pd$OS_Time, 2)
str(pd)

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
