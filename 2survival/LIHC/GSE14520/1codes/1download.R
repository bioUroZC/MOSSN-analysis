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
library(ggplot2)
library(biomaRt)
library(RColorBrewer)
library(gdata)
set.seed(1234)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)

setwd(paste0(PROJ_ROOT, "/1survival/LIHC/GSE14520/data/"))

#=======================================================

#=======================================================

gse <- getGEO("GSE14520", destdir = ".")
gpl <- getGEO("GPL3921", destdir = ".")

colnames(Table(gpl))
Table(gpl)[1:10, 1:6]

gpl <- gpl@dataTable@table
colnames(gpl)
gpl <- gpl %>% dplyr::select(ID, "Gene Symbol" )
write.csv(gpl, file = 'gpl.csv')

genename <- read.csv('gpl.csv', header = TRUE, row.names = 1)
head(genename)
genename$gene <- sub("///.*", "", genename$"Gene.Symbol")
genename$gene <- gdata::trim(genename$gene)
genename$"Gene.Symbol" <- NULL

#=======================================================

#=======================================================

exprSet <- as.data.frame(exprs(gse$`GSE14520-GPL3921_series_matrix.txt.gz`)) 
exprSet$ID = rownames(exprSet)
express = merge( x=genename, y=exprSet, by="ID")
express$ID = NULL
express[which(is.na(express),arr.ind = T)]<-0 
express[1:4,1:4]

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

pd <- fread("GSE14520_Extra_Supplement.txt")
head(pd)

table(pd$`Tissue Type`)
pd <- subset(pd, pd$`Tissue Type` == "Tumor")

pd <- subset(pd, select=c("Affy_GSM",  "Gender", "Age", "TNM staging", "Survival status", "Survival months"))
head(pd)

names(pd) <- c("Sample", 'Gender', "Age", "Stage", "OS", "OS_Time")
head(pd)

table(pd$Gender)
pd$Gender[pd$Gender=="M"] <- "Male"
pd$Gender[pd$Gender=="F"] <- "Female"
table(pd$Gender)

str(pd)

table(pd$Stage)
pd$Stage[pd$Stage=="."] <-NA
pd$Stage[pd$Stage==" "] <-NA
pd$Stage[pd$Stage=="I"] <-"I"
pd$Stage[pd$Stage=="II"] <- "II"
pd$Stage[pd$Stage=="III"] <- "III"
pd$Stage[pd$Stage=="IIIA"] <- "III"
pd$Stage[pd$Stage=="IIIB"] <- "III"
pd$Stage[pd$Stage=="IIIC"] <- "III"
table(pd$Stage)

pd$OS  <-  as.numeric(as.character(pd$OS ))

pd$OS_Time <-  as.numeric(as.character(pd$OS_Time))
pd$OS_Time <- pd$OS_Time/ 12
pd$OS_Time <- round(pd$OS_Time, 2)
pd <- subset(pd, pd$OS_Time > 0)
str(pd)

pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]
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


