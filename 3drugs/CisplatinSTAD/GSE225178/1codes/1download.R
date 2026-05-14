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

setwd("/proj/c.zihao/work1/3drugs/CisplatinSTAD/GSE225178/data/")
set.seed(1234)

#=======================================================

#=======================================================

gsename <- "GSE225178"
gse <- getGEO(gsename, destdir = ".")

gpl <- getGEO('GPL21185', destdir = ".")
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


exprSet <- as.data.frame(exprs(gse$GSE225178_series_matrix.txt.gz))
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

exprSet <- exprSet + 16

min(exprSet)
max(exprSet)

#=======================================================

#=======================================================

pd <- pData(gse$GSE225178_series_matrix.txt.gz)
head(pd)
names(pd)

table(pd$'histologic response based on papers (pmid:ch1')

pd$grade <- sub(".*: ", "", pd$`histologic response based on papers (pmid:ch1`)

pd$Response <- ifelse(pd$grade %in% c("2", "3"), "Sensitive",
                      ifelse(pd$grade %in% c("0", "1a", "1b"), "Resistant", NA))



pd <- subset(pd, select=c( "geo_accession", 
                           "Response" 
))

pd <- na.omit(pd)

head(pd)
colnames(pd) <- c("Sample",  "Response")
head(pd)

table(pd$Response)
pd$Response[pd$Response=="Sensitive"] <- 3
pd$Response[pd$Response=="Resistant"] <- 2
table(pd$Response)

str(pd)

pd$Response <- as.numeric(as.character(pd$Response))

pd <- na.omit(pd)
print(table(pd$Response))

#=======================================================

#=======================================================

samplesname <- intersect(pd$Sample, colnames(exprSet))
samplesname <- unique(samplesname)
pd <- pd[which(pd$Sample %in% samplesname),]
exprSet <- exprSet[,which(colnames(exprSet) %in% samplesname)]
colnames(exprSet)

print(dim(exprSet))

#=======================================================

#=======================================================

exprSet[1:5,1:5]
write.csv(exprSet, file = "exprSet.csv")
write.csv(pd, file = "pd.csv")
