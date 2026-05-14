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
library(openxlsx)

set.seed(1234)


Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)

setwd("/proj/c.zihao/work1/1survival/STAD/GSE13861/data")


#=======================================================

#=======================================================

gsename <- "GSE13861"
gse <- getGEO(gsename, destdir = ".")
gpl <- getGEO('GPL6884', destdir = ".")
colnames(Table(gpl))
Table(gpl)[1:10, 1:6]
gpl <- gpl@dataTable@table
colnames(gpl)


gpl <- gpl %>% dplyr::select(ID, "Symbol")
write.csv(gpl, file = 'gpl.csv')

gpl <- read.csv('gpl.csv', header = TRUE, row.names = 1)
head(gpl)

gpl$gene <- gdata::trim(gpl$Symbol)
colnames(gpl)
gpl <- subset(gpl, select=c("ID", "gene"))

#=======================================================

#=======================================================

exprSet <- as.data.frame(exprs(gse$GSE13861_series_matrix.txt.gz))
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


pd <- pData(gse$GSE13861_series_matrix.txt.gz)
head(pd)
names(pd)

table(pd$`patient:ch1`)

pd <- subset(pd, select=c("geo_accession", 
                          "title"))

head(pd)

yg_parts <- str_extract(pd$title, "YG[0-9A-Za-z_]+")
yg_clean <- str_replace(yg_parts, "_.*", "")
yg_clean <- gsub('V3', '', yg_clean)
pd$Patient <- yg_clean

clinical <- read.xlsx("41467_2018_4179_MOESM5_ESM.xlsx", sheet = 2)

names(clinical)[1] <- "Patient"

intersect(pd$Patient, clinical$Patient)


pd <- merge(pd, clinical, by="Patient")
head(pd)
names(pd)

pd <- subset(pd, select=c( "geo_accession" ,  "Age"  , "Sex",
                           "AJCC6" ,"Death.(1=yes,.0=no)"  ,    "OS.m"        ))

head(pd)
colnames(pd) <- c("Sample",  "Age", "Gender", "Stage", "OS", "OS_Time")
head(pd)

str(pd)

table(pd$Gender)
pd$Gender[pd$Gender=="F"] <- "Female"
pd$Gender[pd$Gender=="M"] <- "Male"
table(pd$Gender)

table(pd$Stage)

str(pd)
table(pd$OS)

str(pd)
table(pd$OS_Time)
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
