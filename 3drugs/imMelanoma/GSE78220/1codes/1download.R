#=======================================================

#=======================================================

rm(list=ls())
library(dplyr)
library(tibble)
library(tidyr)
library(GEOquery)
library(gdata)
library(data.table)
library(readxl)
options(stringsAsFactors = FALSE)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)
setwd("/proj/c.zihao/work1/3drugs/immune/GSE78220/data/")

#=======================================================

#=======================================================


gsename = "GSE78220"
gse<- getGEO(gsename, destdir = ".") 


express <- read_excel('GSE78220_PatientFPKM.xlsx')
express <- as.data.frame(express)
str(express)
express[which(is.na(express),arr.ind = T)]<-0 
express[1:5,1:5]

exprSet <- aggregate(x = express[,2:ncol(express)],
                     by = list(express$Gene),
                     FUN = max)
head(exprSet)[1:5,1:5]

names(exprSet)[1] <- 'ID'
rownames(exprSet) <- exprSet$ID
exprSet$ID <- NULL
exprSet[1:5,1:5]

min(exprSet)
max(exprSet)

#=======================================================

#=======================================================

pd <- pData(gse$GSE78220_series_matrix.txt.gz)
head(pd)
names(pd)

table(pd$'tissue:ch1')

pd <- subset(pd, select=c("title"  , 
                          'age (yrs):ch1',
                          'gender:ch1',
                          'anti-pd-1 response:ch1',
                          'disease status:ch1', 
                          "vital status:ch1",
                          'overall survival (days):ch1' ))

names(pd) <- c("Sample", 'Age', "Gender", "Response",
               "Stage", "OS", "OS_Time")


str(pd)

table(pd$Age)
pd$Age <- as.numeric(as.character(pd$Age))
table(pd$Age)


table(pd$Gender)
pd$Gender[pd$Gender=="F"] <- 'Female'
pd$Gender[pd$Gender=="M"] <- 'Male'
table(pd$Gender)

table(pd$Response)

table(pd$OS)
pd$OS[pd$OS=="Alive"] <- 0
pd$OS[pd$OS=="Dead"] <- 1
pd$OS <- as.numeric(as.character(pd$OS))
table(pd$OS)

table(pd$OS_Time)
pd$OS_Time[pd$OS_Time=="NA"] <- NA
pd$OS_Time <- as.numeric(as.character(pd$OS_Time))
pd$OS_Time <- pd$OS_Time / 365
pd$OS_Time <- round(pd$OS_Time, 2)
table(pd$OS_Time)

# pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]
# pd <- subset(pd, pd$OS_Time > 0)

table(pd$OS)
mean(pd$OS_Time)


pd$Sample
colnames(exprSet)
colnames(exprSet) <- sub("\\..*$", "", colnames(exprSet))
colnames(exprSet)


pd <- subset(pd, select=c("Sample", 'Response'))
pd <- na.omit(pd)


print(table(pd$Response))
pd$Response[pd$Response=="Complete Response"] <- 4
pd$Response[pd$Response=="Partial Response"] <- 3
pd$Response[pd$Response=="Progressive Disease"] <- 1
print(table(pd$Response))

pd$Response <- as.numeric(as.character(pd$Response))

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


