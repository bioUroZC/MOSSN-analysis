# ===============================================================

# ===============================================================

rm(list=ls())
setwd("/proj/c.zihao/work1/1survival/BLCA/GSE32894/data/")
library(dplyr)
library(tibble)
library(tidyr)
library(GEOquery)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)

# ===============================================================
#
# ===============================================================

gse<- getGEO("GSE32894", destdir = ".") 
gpl<- getGEO('GPL6947', destdir = ".") 

colnames(Table(gpl))
Table(gpl)[1:10,1:6] 
gpl <- gpl@dataTable@table
colnames(gpl)
gpl <- gpl %>% dplyr::select(ID, "Symbol")
write.csv(gpl,"GPL.csv", row.names = F)
genename = read.csv("GPL.csv")

# ===============================================================
#
# ===============================================================

exprSet <- as.data.frame(exprs(gse$GSE32894_series_matrix.txt.gz)) 
exprSet$ID = rownames(exprSet)
express = merge( x=genename, y=exprSet, by="ID")
express$ID = NULL
express[which(is.na(express),arr.ind = T)]<-0 
express[1:4,1:4]

exprSet <- aggregate(x = express[,2:ncol(express)],
                     by = list(express$Symbol),
                     FUN = max)
head(exprSet)[1:4,1:4]
exprSet <- as.data.frame(exprSet)
exprSet[1:4,1:4]

exprSet <-exprSet[-1,]
names(exprSet)[1] <- 'ID'
rownames(exprSet) <- exprSet$ID
exprSet$ID <- NULL

min(exprSet)
max(exprSet)

exprSet <- exprSet + 8

min(exprSet)
max(exprSet)

# ===============================================================
#
# ===============================================================

pd <- pData(gse$GSE32894_series_matrix.txt.gz)
head(pd)

pd <- subset(pd, select=c( "geo_accession" ,
                           "age:ch1",  "gender:ch1",   "tumor_stage:ch1" , 
                         "dod_event_(yes/no):ch1" ,  "time_to_dod_(months):ch1"  ))
names(pd) <- c("Sample", "Age", "Gender", "Stage", "OS", "OS_Time")
table(pd$OS)
str(pd)

pd$Age <- as.numeric(as.character(pd$Age))
table(pd$Gender)
pd$Gender <- ifelse(pd$Gender == 'F', "Female", "Male")
table(pd$Gender)

table(pd$OS)
pd$OS <- ifelse(pd$OS == 'no', 0, 1)
pd$OS <- as.numeric(as.character(pd$OS))
table(pd$OS)


pd$OS_Time <- as.numeric(as.character(pd$OS_Time))
pd$OS_Time <- pd$OS_Time / 12
pd$OS_Time <- round(pd$OS_Time, 2)
pd$OS_Time[1:5]
str(pd)

pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]
pd <- subset(pd, pd$OS_Time > 0)
str(pd)

table(pd$Gender)
mean(pd$Age)
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

