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

setwd("/proj/c.zihao/work1/1survival/OV/GSE53963/data/")
set.seed(1234)

#=======================================================

#=======================================================

gsename <- "GSE53963"
gse <- getGEO(gsename, destdir = ".")
gpl <- getGEO('GPL6480', destdir = ".")
colnames(Table(gpl))
Table(gpl)[1:10, 1:6]
gpl <- gpl@dataTable@table
colnames(gpl)


gpl <- gpl %>% dplyr::select(ID, "GENE_SYMBOL")
write.csv(gpl, file = 'gpl.csv')

gpl <- read.csv('gpl.csv', header = TRUE, row.names = 1)
head(gpl)

gpl$gene <- gdata::trim(gpl$GENE_SYMBOL)
colnames(gpl)
gpl <- subset(gpl, select=c("ID", "gene"))

#=======================================================

#=======================================================


exprSet <- as.data.frame(exprs(gse$GSE53963_series_matrix.txt.gz))
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

min(exprSet)
max(exprSet)

#=======================================================

#=======================================================

pd <- pData(gse$GSE53963_series_matrix.txt.gz)
head(pd)
names(pd)

table(pd$'tcga_sampleid:ch2')
pd <- pd[!grepl("TCGA", pd$'tcga_sampleid:ch2'), ]



pd <- subset(pd, select=c( "geo_accession", 
                           'age_at_dx:ch2',
                           'Stage:ch2',
                           "vital_status:ch2"  ,
                           "time_fu_months:ch2" 
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
pd$OS[pd$OS=="Dead"] <- 1
pd$OS[pd$OS=="Alive"] <- 0
pd$OS <-  as.numeric(as.character(pd$OS))
table(pd$OS)


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
