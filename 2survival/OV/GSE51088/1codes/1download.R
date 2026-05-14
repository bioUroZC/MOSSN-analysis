#=======================================================

#=======================================================

rm(list=ls())
library(dplyr)
library(tibble)
library(tidyr)
library(GEOquery)
library(gdata)

options(stringsAsFactors = FALSE)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)
setwd("/proj/c.zihao/work1/1survival/OV/GSE51088/data/")

#=======================================================

#=======================================================

gse<- getGEO("GSE51088", destdir = ".") 
gpl<- getGEO('GPL7264', destdir = ".") 

colnames(Table(gpl))
gpl <- gpl@dataTable@table
colnames(gpl)
gpl <- gpl %>% dplyr::select(ID, "GENE_SYMBOL" )
write.csv(gpl,"GPL.csv", row.names = F)
genename = read.csv("GPL.csv")

genename$gene <- gdata::trim(genename$GENE_SYMBOL)
genename$GENE_SYMBOL <- NULL


#=======================================================

#=======================================================

exprSet <- as.data.frame(exprs(gse$GSE51088_series_matrix.txt.gz)) 
exprSet$ID = rownames(exprSet)
express = merge( x=genename, y=exprSet, by="ID")
express$ID = NULL
express[which(is.na(express),arr.ind = T)]<-0 
express[1:4,1:4]

exprSet <- aggregate(x = express[,2:ncol(express)],
                     by = list(express$gene),
                     FUN = max)
head(exprSet)
exprSet <- as.data.frame(exprSet)
exprSet[1:4,1:4]

exprSet <-exprSet[-1,]
names(exprSet)[1] <- 'ID'
rownames(exprSet) <- exprSet$ID
exprSet$ID <- NULL
exprSet[1:4,1:4]

min(exprSet)
max(exprSet)

exprSet <- exprSet + 2

min(exprSet)
max(exprSet)

#=======================================================

#=======================================================

pd <- pData(gse$GSE51088_series_matrix.txt.gz)
head(pd)
names(pd)
table(pd$`tcga id:ch2`)
pd <- pd[!grepl("TCGA", pd$`tcga id:ch2`), ]

table(pd$'tissue status:ch2')
pd <- subset(pd, pd$`tissue status:ch2` == 'Primary')

pd <- subset(pd, select=c( "geo_accession" ,
                           "age at surgery:ch2", 
                           'patient status:ch2',
                           "follow up months:ch2"))

names(pd) <- c("Sample",  "Age", "OS", "OS_Time")

str(pd)

pd$Age <- as.numeric(as.character(pd$Age))



table(pd$OS)
pd$OS[pd$OS=="Dead"] <- 1
pd$OS[pd$OS=="Alive"] <- 0
pd$OS[pd$OS=="Unknown"] <- NA
pd$OS <-  as.numeric(as.character(pd$OS))
table(pd$OS)


str(pd)
table(pd$OS_Time)
pd$OS_Time[pd$OS_Time=="Unknown"] <- NA
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

