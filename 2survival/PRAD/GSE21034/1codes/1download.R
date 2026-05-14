#=======================================================

#=======================================================

rm(list=ls())

library(dplyr)
library(tibble)
library(tidyr)
library(GEOquery)
library(gdata)
library(data.table)
library(biomaRt)
library(readxl)
options(stringsAsFactors = FALSE)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)
setwd("/proj/c.zihao/work1/1survival/PRAD/GSE21034/data")

#=======================================================

#=======================================================


gpl <- getGEO('GPL10264', destdir = ".") 
colnames(Table(gpl))
gpl <- gpl@dataTable@table
colnames(gpl)
head(gpl)
names(gpl)[2] <- 'refseq_mrna'


mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
nm_ids <- gpl$refseq_mrna
anno <- getBM(attributes = c('refseq_mrna', 'hgnc_symbol'),
             filters = 'refseq_mrna',
             values = nm_ids,
             mart = mart)

head(anno)
head(gpl)

anno_unique <- anno[!duplicated(anno$refseq_mrna), ]
gpl_unique  <- gpl[!duplicated(gpl$refseq_mrna), ]
merged <- merge(gpl_unique, anno_unique, by = "refseq_mrna", all.x = TRUE)

head(merged)
merged$refseq_mrna <- NULL

#=======================================================

#=======================================================


gse <- getGEO(filename = "GSE21034-GPL10264_series_matrix.txt.gz", GSEMatrix=TRUE)
exprSet <- as.data.frame(exprs(gse))
str(exprSet)

exprSet$ID = rownames(exprSet)
express = merge( x=merged, y=exprSet, by="ID")
express$ID = NULL

express[which(is.na(express),arr.ind = T)]<-0 
express[1:5,1:5]
names(express)[1] <- 'gene'


exprSet <- aggregate(x = express[,2:ncol(express)],
                     by = list(express$gene),
                     FUN = max)
head(exprSet)

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

pd <- pData(gse)
head(pd)
names(pd)
table(pd$`disease status:ch1`)
pd <- subset(pd, pd$`disease status:ch1` == "prostate cancer")

pd <- subset(pd, select=c("sample id:ch1","geo_accession", 'clint_stage:ch1' ))
names(pd) <- c("ID", "Sample", "Stage")

df <- read_excel("NIHMS271566-supplement-05.xls")
df <- as.data.frame(df)
names(df)
df <- subset(df, select=c("Sample ID",  "BCR_Event", "BCR_FreeTime" ))
names(df) <- c("ID", "OS", "OS_Time")
head(df)

pd <- merge(pd, df, by="ID")

head(pd)
str(pd)

pd$ID <- NULL

table(pd$OS)
pd$OS[pd$OS=='NA'] <- NA
pd$OS[pd$OS=='NO'] <- 0
pd$OS[pd$OS=='BCR_Algorithm'] <- 1
pd$OS <- as.numeric(as.character(pd$OS))
table(pd$OS)


table(pd$OS_Time)
pd$OS_Time[pd$OS_Time=='NA'] <- NA
pd$OS_Time <- as.numeric(as.character(pd$OS_Time))
pd$OS_Time <- pd$OS_Time / 12
pd$OS_Time <- round(pd$OS_Time, 2)
table(pd$OS_Time)

pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]
pd <- subset(pd, pd$OS_Time > 0)

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

