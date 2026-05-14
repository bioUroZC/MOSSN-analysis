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

setwd("/proj/c.zihao/work1/1survival/STAD/GSE29272/data")

#=======================================================

#=======================================================

gsename <- "GSE29272"
gse <- getGEO(gsename, destdir = ".")
gpl <- getGEO('GPL96', destdir = ".")
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


exprSet <- as.data.frame(exprs(gse$GSE29272_series_matrix.txt.gz))
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


pd <- pData(gse[["GSE29272_series_matrix.txt.gz"]])
head(pd)
names(pd)

pd <- pd %>%
  dplyr::select( title, geo_accession,
                 "tissue:ch1")
str(pd)

table(pd$`tissue:ch1`)

pd <- subset(pd, pd$`tissue:ch1` != "adjacent tissue normal gastric glands")

pd$title[1:5]

pd$title <- sub(".*\\s+", "", pd$title)

pd$title[1:5]



clinical1 <- read.xlsx("pone.0063826.s001.xlsx", sheet = "S1a GCA")
names(clinical1)

clinical2 <- read.xlsx("pone.0063826.s001.xlsx", sheet = "S1b GNCA")
names(clinical2)

clinical <- rbind(clinical1, clinical2)
names(clinical)


clinical <- clinical %>%
  dplyr::select("ID", "Sex",  "Age" , "Tumor.Stage", 
                "Surivval.during.follow.up" , "Survival.(days)" ,
                "Cause.of.death" )

names(clinical)[1] <- "title"

clinical$title <- paste0(clinical$title, "T")

pd <- merge(pd, clinical, by="title")

str(pd)


pd <- subset(pd, select=c( "geo_accession" ,   "Age" , "Sex" , "Tumor.Stage" ,
                           "Surivval.during.follow.up"  , 
                           "Survival.(days)"  ))

str(pd)

names(pd) <- c("Sample", "Age", "Gender", "Stage", "OS", "OS_Time")
str(pd)

pd$Gender <- gdata::trim(pd$Gender)
table(pd$Gender)
pd$Gender[pd$Gender=="F"] <- 'Female'
pd$Gender[pd$Gender=="M"] <- 'Male'
table(pd$Gender)

table(pd$Stage)

str(pd)
table(pd$OS)
pd$OS[pd$OS=="N"] <- 1
pd$OS[pd$OS=="Y"] <- 0
pd$OS[pd$OS=="Unknown"] <- NA
table(pd$OS)

str(pd)
table(pd$OS_Time)
pd$OS_Time[pd$OS_Time=="Unknown"] <- NA
pd$OS_Time <- as.numeric(as.character(pd$OS_Time))
pd$OS_Time <- pd$OS_Time / 365
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


