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
set.seed(1234)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)

setwd("/proj/c.zihao/work1/1survival/LUAD/GSE41271/data/")

#=======================================================

#=======================================================

gsename <- "GSE41271"
gse <- getGEO(gsename, destdir = ".")
gplname <- 'GPL6884'
gpl <- getGEO(gplname, destdir = ".")

colnames(Table(gpl))
Table(gpl)[1:10, 1:6]
gpl <- gpl@dataTable@table
colnames(gpl)
gpl <- gpl %>% dplyr::select(ID, "Symbol" )
write.csv(gpl, file = 'gpl.csv')

gpl <- read.csv('gpl.csv', header = TRUE, row.names = 1)
head(gpl)

gpl$Symbol <- gdata::trim(gpl$Symbol)

#=======================================================

#=======================================================

exprSet <- as.data.frame(exprs(gse$GSE41271_series_matrix.txt.gz))
str(exprSet)

exprSet$ID = rownames(exprSet)
express = merge( x=gpl, y=exprSet, by="ID")
express$ID = NULL

express[which(is.na(express),arr.ind = T)]<-0 
express[1:5,1:5]

exprSet <- aggregate(x = express[,2:ncol(express)],
                     by = list(express$Symbol),
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

pd <- pData(gse$GSE41271_series_matrix.txt.gz)
table(pd$`histology:ch1`)
pd <-  subset(pd,pd$`histology:ch1`=="Adenocarcinoma" )
head(pd)

pd <- subset(pd, select=c("geo_accession", "gender:ch1", "final patient stage:ch1",
                          "date of birth:ch1", "date of surgery:ch1", 
                          "vital statistics:ch1", "last follow-up survival:ch1"))

head(pd)

pd$`date of birth:ch1` <- as.Date(pd$`date of birth:ch1`, format="%Y-%m-%d")
pd$`date of surgery:ch1` <- as.Date(pd$`date of surgery:ch1`, format="%Y-%m-%d")
pd$`last follow-up survival:ch1` <- as.Date(pd$`last follow-up survival:ch1`, format="%Y-%m-%d")


pd$Age  <- difftime(pd$`date of surgery:ch1`, pd$`date of birth:ch1`, units="days")
pd$OS.Time <- difftime(pd$`last follow-up survival:ch1`, pd$`date of surgery:ch1`, units="days")
pd$`date of birth:ch1` <- NULL
pd$`date of surgery:ch1` <- NULL
pd$`last follow-up survival:ch1`<- NULL

head(pd)

names(pd) <- c("Sample", 'Gender', "Stage", "OS", "Age", "OS_Time")
head(pd)


table(pd$Gender)
pd$Gender[pd$Gender=="F"] <- "Female"
pd$Gender[pd$Gender=="M"] <- "Male"
table(pd$Gender)


table(pd$OS)
pd$OS[pd$OS=="A"] <- 0
pd$OS[pd$OS=="D"] <- 1
table(pd$OS)
str(pd)

head(pd)

pd$Age <- as.numeric(pd$Age)
pd$OS_Time <- gsub("days", "", pd$OS_Time)
table(pd$OS_Time)
pd$OS_Time <- as.numeric(as.character(pd$OS_Time))

pd$Age <- pd$Age / 365
pd$OS_Time <- pd$OS_Time / 365
head(pd)
str(pd)

pd <- subset(pd, select=c( 'Sample',  'Gender',  'Age', "Stage", 'OS',  'OS_Time' ))
pd$OS  <-  as.numeric(as.character(pd$OS))
pd$Age <- round(pd$Age, 0)
pd$OS_Time <- round(pd$OS_Time, 2)

str(pd)

pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]
pd <- subset(pd, pd$OS_Time>0)
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
