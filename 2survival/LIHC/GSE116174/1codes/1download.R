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
library(readxl)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)

setwd(paste0(PROJ_ROOT, "/1survival/LIHC/GSE116174/data/"))

# ===================================================

# ===================================================

gse <- getGEO("GSE116174", destdir = ".")
gpl<- getGEO("GPL13158", destdir = ".") 

colnames(Table(gpl))
Table(gpl)[1:10, 1:6]

gpl <- gpl@dataTable@table
colnames(gpl)
gpl <- gpl %>% dplyr::select(ID, "Gene Symbol")
write.csv(gpl, file = 'gpl.csv')

genename <- read.csv('gpl.csv', header = TRUE, row.names = 1)
head(genename)

genename$gene <- gdata::trim(genename$"Gene.Symbol")
genename$'Gene.Symbol' <- NULL
head(genename)

# ===================================================

# ===================================================


exprSet <- as.data.frame(exprs(gse$GSE116174_series_matrix.txt.gz)) 
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

# ===================================================

# ===================================================

pd <- pData(gse$GSE116174_series_matrix.txt.gz)
pd <- subset(pd, select=c("geo_accession", "source_name_ch1"))
names(pd) <- c("Sample", "sampleID")

OS <- read_excel("GSE116174_HCC-64-u133_plus_2_clinical_data.xls")
OS <- as.data.frame(OS)

pd <- merge(pd, OS, by="sampleID")
head(pd)

names(pd)

pd <- pd %>%
  dplyr::select( "Sample", "Age" , "Gender"  , "clinstage" , "Event_death",   "Follow_up time(month)"  )

names(pd) <- c( "Sample", "Age" , "Gender", "Stage", "OS", "OS_Time")
str(pd)

table(pd$Gender)
pd$Gender <- ifelse(pd$Gender == 'F', "Female", "Male")
table(pd$Gender)


str(pd)


table(pd$Stage)
pd$Stage[pd$Stage=="I"] <-"I"
pd$Stage[pd$Stage=="Ⅰ"] <-"I"
pd$Stage[pd$Stage=="Ⅱ"] <- "II"
pd$Stage[pd$Stage=="Ⅲ"] <- "III"
table(pd$Stage)


pd$OS_Time <- pd$OS_Time / 12
pd <- subset(pd, pd$OS_Time > 0)


pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]

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

