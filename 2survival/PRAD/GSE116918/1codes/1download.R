#=======================================================

#=======================================================

rm(list=ls())

library(dplyr)
library(tibble)
library(tidyr)
library(GEOquery)
library(gdata)
library(data.table)

options(stringsAsFactors = FALSE)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)
setwd("/proj/c.zihao/work1/1survival/PRAD/GSE116918/data")


#=======================================================

#=======================================================

gsename = "GSE116918"
gse<- getGEO(gsename, destdir = ".") 

gpl<- getGEO('GPL25318', destdir = ".") 
colnames(Table(gpl))
Table(gpl)[1:10,1:6]
gpl <- gpl@dataTable@table
colnames(gpl)
gpl <- gpl %>% dplyr::select(ID,  "Gene Symbol" )

write.csv(gpl, file = 'gpl.csv')

gpl <- read.csv('gpl.csv', header = T, row.names = 1)
head(gpl)

gpl <- gpl %>%
  tidyr::separate(Gene.Symbol, c('gene', 'symbol'), sep='\\///'   )%>%
  dplyr::select("ID", 'gene')

gpl$gene <- gdata::trim(gpl$gene)


#=======================================================

#=======================================================

exprSet <- as.data.frame(exprs(gse$GSE116918_series_matrix.txt.gz))
str(exprSet)

exprSet$ID = rownames(exprSet)
express = merge( x=gpl, y=exprSet, by="ID")
express$ID = NULL

express[which(is.na(express),arr.ind = T)]<-0 
express[1:5,1:5]

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

pd <- pData(gse$GSE116918_series_matrix.txt.gz)
names(pd)
head(pd)

table(pd$`tissue:ch1`)

pd <- subset(pd, select=c("geo_accession",  
                          'patient age (years):ch1', 
                          "bcr event (1=yes, 0=no):ch1", 
                          'follow-up time (bcr, months):ch1'))

names(pd) <- c("Sample", "Age", "OS",  "OS_Time")

str(pd)


table(pd$Age)
pd$Age <- as.numeric(as.character(pd$Age))
table(pd$Age)
str(pd)


table(pd$OS)
pd$OS <- as.numeric(as.character(pd$OS))
table(pd$OS)

table(pd$OS_Time)
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


