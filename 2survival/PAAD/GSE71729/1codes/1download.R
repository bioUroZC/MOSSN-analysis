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
set.seed(1234)
library(stringr)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)

setwd("/proj/c.zihao/work1/1survival/PAAD/GSE71729/data/")

#=======================================================

#=======================================================

gsename <- "GSE71729"
gse <- getGEO(gsename, destdir = ".")
gpl <- getGEO('GPL20769', destdir = ".")
colnames(Table(gpl))

gpl <- gpl@dataTable@table
colnames(gpl)
gpl <- gpl %>% dplyr::select(ID,  "ORF")
write.csv(gpl, file = 'gpl.csv')

gpl <- read.csv('gpl.csv', header = TRUE, row.names = 1)
head(gpl)

gpl$gene <- gdata::trim(gpl$ORF)
colnames(gpl)
gpl <- subset(gpl, select=c("ID", "gene"))


#=======================================================

#=======================================================

exprSet <- as.data.frame(exprs(gse$GSE71729_series_matrix.txt.gz))
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


pd <- pData(gse$GSE71729_series_matrix.txt.gz)
head(pd)
names(pd)

table(pd$"tissue type:ch2")
pd <- subset(pd, pd$"tissue type:ch2" == 'Primary')


pd <- subset(pd, select=c("geo_accession", 
                          "death_event_1death_0censor:ch2"  ,
                          "survival_months:ch2"
))

head(pd)
colnames(pd) <- c("Sample",  "OS", "OS_Time")
head(pd)

str(pd)

table(pd$OS)
pd$OS[pd$OS=="1"] <- 1
pd$OS[pd$OS=="0"] <- 0
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
