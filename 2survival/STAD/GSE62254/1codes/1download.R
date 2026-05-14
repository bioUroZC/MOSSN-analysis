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
library(readxl)
set.seed(1234)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)

setwd("/proj/c.zihao/work1/1survival/STAD/GSE62254/data")


#=======================================================

#=======================================================

gsename <- "GSE62254"
gse <- getGEO(gsename, destdir = ".")
gpl <- getGEO('GPL570', destdir = ".")
colnames(Table(gpl))
Table(gpl)[1:10, 1:6]


gpl <- gpl@dataTable@table
colnames(gpl)
gpl <- gpl %>% dplyr::select(ID, "Gene Symbol")
write.csv(gpl, file = 'gpl.csv')

gpl <- read.csv('gpl.csv', header = TRUE, row.names = 1)
head(gpl)

gpl <- gpl %>%
  tidyr::separate("Gene.Symbol", c('gene', 'symbol'), sep = '\\///') %>%
  dplyr::select("ID", 'gene')

gpl$gene <- gdata::trim(gpl$gene)

#=======================================================

#=======================================================


exprSet <- as.data.frame(exprs(gse$GSE62254_series_matrix.txt.gz))
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


pd <- pData(gse[["GSE62254_series_matrix.txt.gz"]])
head(pd)
names(pd)

pd <- pd %>%
  dplyr::select( title, geo_accession, 
                 "tissue:ch1" )

table(pd$`tissue:ch1`)

str(pd)

library(readxl)
clinical <- read_excel("41591_2015_BFnm3850_MOESM34_ESM.xls", sheet = "FINAL")
names(clinical)

clinical <- subset(clinical, select=c("Tumor ID", "Sample\nName","age", "sex" ,"pStage" ,
                                      "FU status0=alive without ds, 1=alive with recurren ds, 2=dead without ds, 3=dead d/t recurrent ds, 4=dead, unknown, 5= FU loss" ,
                                      "OS\n(months)" ))

table(clinical$`FU status0=alive without ds, 1=alive with recurren ds, 2=dead without ds, 3=dead d/t recurrent ds, 4=dead, unknown, 5= FU loss`)

clinical$OS <- ifelse(clinical$`FU status0=alive without ds, 1=alive with recurren ds, 2=dead without ds, 3=dead d/t recurrent ds, 4=dead, unknown, 5= FU loss` %in% c(0,1), 
                                     "alive", 
                                     "death")

clinical$`FU status0=alive without ds, 1=alive with recurren ds, 2=dead without ds, 3=dead d/t recurrent ds, 4=dead, unknown, 5= FU loss` <- NULL

clinical$"Tumor ID"

clinical$title <- paste0("T", clinical$`Tumor ID`)

intersect(clinical$title, pd$title)

pd <- merge(pd, clinical, by='title')
names(pd)
pd <- subset(pd, select=c( "geo_accession", "age" , "sex"   ,   "pStage"     ,
                           "OS", "OS\n(months)" ))

names(pd) <- c("Sample", "Age", "Gender", "Stage", "OS", "OS_Time")
str(pd)

table(pd$Gender)
pd$Gender[pd$Gender=="F"] <- 'Female'
pd$Gender[pd$Gender=="M"] <- 'Male'
table(pd$Gender)

table(pd$Stage)
pd$Stage_simple <- NA  
pd$Stage_simple[grepl("^I$", pd$Stage)] <- "I"
pd$Stage_simple[grepl("^II$", pd$Stage)] <- "II"
pd$Stage_simple[grepl("^III$", pd$Stage)] <- "III"
pd$Stage_simple[grepl("^IV$", pd$Stage)] <- "IV"
pd$Stage_simple[grepl("II", pd$Stage) & is.na(pd$Stage_simple)] <- "II"
pd$Stage_simple[grepl("III", pd$Stage) & is.na(pd$Stage_simple)] <- "III"
pd$Stage_simple[grepl("IV", pd$Stage) & is.na(pd$Stage_simple)] <- "IV"
pd$Stage_simple[grepl("I", pd$Stage) & is.na(pd$Stage_simple)] <- "I"
table(pd$Stage_simple, useNA = "ifany")

str(pd)     
pd$Stage <- pd$Stage_simple   
pd$Stage_simple    <- NULL       


str(pd)
table(pd$OS)
pd$OS[pd$OS=="death"] <- 1
pd$OS[pd$OS=="alive"] <- 0
table(pd$OS)

str(pd)
table(pd$OS_Time)
pd$OS_Time <- pd$OS_Time / 12
pd$OS_Time <- round(pd$OS_Time, 2)

pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]
pd <- subset(pd, pd$OS_Time>0)
head(pd)

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


