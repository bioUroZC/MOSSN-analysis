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

setwd("/proj/c.zihao/work1/1survival/GBM/GSE108474/data")


#=======================================================

#=======================================================

gsename <- "GSE108474"
gse <- getGEO(gsename, destdir = ".")
gpl <- getGEO('GPL570', destdir = ".")
colnames(Table(gpl))
Table(gpl)[1:10, 1:6]


gpl <- gpl@dataTable@table
colnames(gpl)
gpl <- gpl %>% dplyr::select(ID, "Gene Symbol")
write.csv(gpl, file = 'gpl.csv')

gpl <- read.csv('gpl.csv', header = TRUE, row.names = 1)
head(gpl)

gpl <- gpl %>%
  tidyr::separate("Gene.Symbol", c('gene', 'symbol'), sep = '\\///') %>%
  dplyr::select("ID", 'gene')

gpl$gene <- gdata::trim(gpl$gene)

#=======================================================

#=======================================================


exprSet <- as.data.frame(exprs(gse$GSE108474_series_matrix.txt.gz))
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

exprSet <- exprSet + 4

min(exprSet)
max(exprSet)

#=======================================================

#=======================================================

pd <- pData(gse$GSE108474_series_matrix.txt.gz)
head(pd)
names(pd)

table(pd$`tumor grade:ch1`)

table(pd$`disease:ch1`)

pd <- subset(pd, pd$`disease:ch1` == "glioblastoma multiforme")

table(pd$`tumor grade:ch1`)

pd <- subset(pd, select=c("title", "geo_accession" ))

clinical <- fread("GSE108474_REMBRANDT_clinical.data.txt")
clinical <- as.data.frame(clinical)
head(clinical)

clinical <- subset(clinical, select=c("SUBJECT_ID", "AGE_RANGE", "GENDER",
                                      'EVENT_OS', 'OVERALL_SURVIVAL_MONTHS', 
                                      "WHO_GRADE"))
intersect(pd$title, clinical$SUBJECT_ID)
names(clinical)[1] <- "title"

pd <- merge(pd, clinical, by="title")

head(pd)

pd$title <- NULL

names(pd) <- c("Sample", "Age", "Gender", "OS", "OS_Time", "Grade")

head(pd)

pd$Age[pd$Age == ""] <- NA

table(pd$Gender)
pd$Gender[pd$Gender=="FEMALE"] <- "Female"
pd$Gender[pd$Gender=="MALE"] <- "Male"
table(pd$Gender)


str(pd)

pd$OS_Time <- pd$OS_Time / 12
pd$OS_Time <- round(pd$OS_Time, 2)
pd <- subset(pd, pd$OS_Time>0)

table(pd$Grade)
table(pd$Age)
sum(is.na(pd$Age))

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


