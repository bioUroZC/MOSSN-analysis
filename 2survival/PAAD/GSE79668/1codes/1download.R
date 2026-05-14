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
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)

setwd("/proj/c.zihao/work1/2survival/PAAD/GSE79668/data/")

#=======================================================

#=======================================================

gsename <- "GSE79668"
gse <- getGEO(gsename, destdir = ".")

express <- fread('GSE79668_51_tumors_sharedgenecounts.txt')
express$Gene_EnsembleID[1:5]

express<- separate(express, Gene_EnsembleID, into = c("gene", "EnsembleID"), sep = "_")
express$EnsembleID <- NULL

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


pd <- pData(gse$GSE79668_series_matrix.txt.gz)
head(pd)
names(pd)

table(pd$`patient:ch1`)

pd <- subset(pd, select=c("title", 
                          "age:ch1", "patient gender:ch1",
                          "t (tnm score):ch1",
                          "n (tnm score):ch1" ,
                          "m (tnm score):ch1"  , 
                          "patient survival status:ch1",
                          "survival time (days):ch1" 
))

head(pd)
colnames(pd) <- c("Sample",  "Age",  "Gender",
                  "Tstage", "Nstage", "Mstage",
                  "OS", "OS_Time")
head(pd)

str(pd)

str(pd)
table(pd$Age)
pd$Age <- as.numeric(as.character(pd$Age))
table(pd$Age)


table(pd$Gender)
pd$Gender[pd$Gender=="female"] <- "Female"
pd$Gender[pd$Gender=="male"] <- "Male"
table(pd$Gender)


str(pd)
table(pd$OS)
pd$OS[pd$OS=="Dead"] <- 1
pd$OS[pd$OS=="Alive"] <- 0
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


pd$Sample <- sub("_SL.*$", "", pd$Sample)
pd$Sample <- sub("_T$", "", pd$Sample)


pd$Sample
colnames(exprSet)


pd$Sample <- gsub("_", "", pd$Sample)
colnames(exprSet) <- gsub("_", "", colnames(exprSet))

pd$Sample <- paste0("A", pd$Sample)
colnames(exprSet) <- paste0("A", colnames(exprSet))

pd$Sample
colnames(exprSet)

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
