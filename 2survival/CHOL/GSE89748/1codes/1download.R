#=======================================================

#=======================================================

rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


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
library(readxl)
library(org.Hs.eg.db)
library(AnnotationDbi)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)

setwd(paste0(PROJ_ROOT, "/1survival/CHOL/GSE89748/data/"))
set.seed(1234)

#=======================================================

#=======================================================

gsename <- "GSE89748"
gse <- getGEO(gsename, destdir = ".")
gpl <- getGEO('GPL10558', destdir = ".")
colnames(Table(gpl))
Table(gpl)[1:10, 1:6]

gpl <- gpl@dataTable@table
colnames(gpl)
gpl <- gpl %>% dplyr::select(ID, "ILMN_Gene" )
write.csv(gpl, file = 'gpl.csv')

gpl <- read.csv('gpl.csv', header = TRUE, row.names = 1)
head(gpl)

gpl$gene <- gdata::trim(gpl$ILMN_Gene)
gpl$ILMN_Gene <- NULL

#=======================================================

#=======================================================


exprSet <- fread('GSE89748_CCA_batch02_illumina_Gene_expression_noNorm_noBKGD.txt')
exprSet <- as.data.frame(exprSet)
exprSet[1:5,1:5]
names(exprSet)[1] <- 'ID'
exprSet <- exprSet %>% dplyr::select(-matches("Pval", ignore.case = TRUE))

express = merge( x=gpl, y=exprSet, by="ID")
express$ID = NULL

express[which(is.na(express),arr.ind = T)]<-0 
express[1:5,1:5]
str(express)


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

pd <- read_excel("NIHMS890532-supplement-3.xlsx")
pd <- as.data.frame(pd)
head(pd)
names(pd)

pd <- subset(pd, select=c( "Sample ID" , 
                           "Age at surgery"  ,
                           "Sex",
                           "Stage" ,
                           "Vital state (1=Dead)"    ,
                           "Overall survival (days)" 
))

head(pd)
colnames(pd) <- c("Sample", "Age", "Gender",  "Stage",  "OS", "OS_Time")
head(pd)

str(pd)


table(pd$Age)
pd$Age[pd$Age=='N/A'] <- NA
pd$Age <- as.numeric(as.character(pd$Age))
pd$Age <- round(pd$Age, 1)
table(pd$Age)


table(pd$Gender)
pd$Gender[pd$Gender=="F"] <- "Female"
pd$Gender[pd$Gender=="M"] <- "Male"
table(pd$Gender)


table(pd$OS)
pd$OS[pd$OS=="N/A"] <- NA
pd$OS[pd$OS=="1"] <- 1
pd$OS[pd$OS=="0"] <- 0
pd$OS <-  as.numeric(as.character(pd$OS))
table(pd$OS)


str(pd)
table(pd$OS_Time)
pd$OS_Time[pd$OS_Time=="N/A"] <- NA
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

head(pd)

colnames(exprSet) <- gsub("_", "", colnames(exprSet))
pd$Sample <- gsub("_", "", pd$Sample)

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
