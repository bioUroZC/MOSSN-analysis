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

setwd("/proj/c.zihao/work1/1survival/OV/GSE9891/data")
set.seed(1234)

#=======================================================

#=======================================================

gsename <- "GSE9891"
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

exprSet <- as.data.frame(exprs(gse$GSE9891_series_matrix.txt.gz))
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

#=======================================================

#=======================================================

pd <- read.csv('GSE9891_eset_clindata.csv', header = T, row.names = 1)
head(pd)
names(pd)


pd$geo <- rownames(pd)

pd <- subset(pd, select=c( 'geo',
                           "sample_type", 
                           "age_at_initial_pathologic_diagnosis",
                           "primarysite" ,
                           "arrayedsite" ,
                           "tumorstage"  ,
                           "vital_status" ,
                           "days_to_death" 
))


table(pd$sample_type)
pd <- subset(pd, pd$sample_type == 'tumor')

table(pd$primarysite)
pd <- subset(pd, pd$primarysite == 'ov')

table(pd$arrayedsite)
pd <- subset(pd, pd$arrayedsite == 'ov')


pd <- subset(pd, select=c( 'geo',
                           "age_at_initial_pathologic_diagnosis",
                           "tumorstage"  ,
                           "vital_status" ,
                           "days_to_death" 
))


head(pd)
colnames(pd) <- c("Sample", "Age",  "Stage",  "OS", "OS_Time")
head(pd)

str(pd)


table(pd$OS)
pd$OS[pd$OS == "living"] <- 0
pd$OS[pd$OS == "deceased"] <- 1
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
