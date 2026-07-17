#=======================================================

#=======================================================

rm(list=ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")

library(dplyr)
library(tibble)
library(tidyr)
library(GEOquery)
library(gdata)
library(data.table)

options(stringsAsFactors = FALSE)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)
setwd(paste0(PROJ_ROOT, "/1survival/PRAD/GSE70769/data"))

#=======================================================

#=======================================================

gsename = "GSE70769"
gse<- getGEO(gsename, destdir = ".") 

gpl<- getGEO('GPL10558', destdir = ".") 
colnames(Table(gpl))
Table(gpl)[1:10,1:6]
gpl <- gpl@dataTable@table
colnames(gpl)
gpl <- gpl %>% dplyr::select(ID,   "ILMN_Gene"  )
head(gpl)

gpl$gene <- gdata::trim(gpl$ILMN_Gene)
gpl$ILMN_Gene <- NULL

head(gpl)

#=======================================================

#=======================================================

exprSet <- as.data.frame(exprs(gse$GSE70769_series_matrix.txt.gz))
str(exprSet)
head(exprSet)[1:5,1:5]

exprSet <- as.data.frame(exprSet)
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
exprSet <-exprSet[-(1:2),]
names(exprSet)[1] <- 'ID'
rownames(exprSet) <- exprSet$ID
exprSet$ID <- NULL
exprSet[1:5,1:5]

min(exprSet)
max(exprSet)

#=======================================================

#=======================================================

pd <- pData(gse$GSE70769_series_matrix.txt.gz)
names(pd)
head(pd)


pd <- subset(pd, select=c("geo_accession",  
                          'clinical stage:ch1', 
                          "biochemical relapse (bcr):ch1", 
                          "time to bcr (months):ch1", 
                          'total follow up (months):ch1'))

names(pd) <- c("Sample", "Stage", "OS",  "Time1", "Time2")



str(pd$Time1)
table(pd$Time1, useNA = "ifany")
pd$Time1[pd$Time1 == "N/A"] <- pd$Time2[pd$Time1 == "N/A"]

pd$OS_Time <- pd$Time1
pd$Time1 <- NULL
pd$Time2 <- NULL

str(pd)


table(pd$OS)
pd$OS[pd$OS=='N/A'] <- NA
pd$OS[pd$OS=='N'] <- 0
pd$OS[pd$OS=='Y'] <- 1
pd$OS <- as.numeric(as.character(pd$OS))
table(pd$OS)

table(pd$OS_Time)
pd$OS_Time[pd$OS_Time=='N/A'] <- NA
pd$OS_Time[pd$OS_Time=='UNKNOWN'] <- NA
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


