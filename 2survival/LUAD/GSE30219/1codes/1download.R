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

setwd(paste0(PROJ_ROOT, "/1survival/LUAD/GSE30219/data/"))


#=======================================================

#=======================================================

gsename <- "GSE30219"
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

exprSet <- as.data.frame(exprs(gse$GSE30219_series_matrix.txt.gz))
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

gene_iqr <- apply(exprSet, 1, IQR)
mean_iqr <- mean(gene_iqr)
median_iqr <- median(gene_iqr)
mean_iqr
median_iqr

#=======================================================

#=======================================================

pd <- pData(gse$GSE30219_series_matrix.txt.gz)
pd <- subset(pd, pd$characteristics_ch1.3 == "histology: ADC")
head(pd)

pd <- subset(pd, select=c("geo_accession", "gender:ch1", "age at surgery:ch1",
                          "status:ch1", "follow-up time (months):ch1",
                          "pt stage:ch1" , "pn stage:ch1",  "pm stage:ch1"))

head(pd)
names(pd) <- c("Sample", 'Gender', "Age","OS", "OS_Time", 
               "T_stages", "N_stages", "M_stages")
head(pd)

table(pd$Gender)
pd$Gender[pd$Gender=="F"] <- "Female"
pd$Gender[pd$Gender=="M"] <- "Male"
table(pd$Gender)

table(pd$OS)
pd$OS[pd$OS=="ALIVE"] <- 0
pd$OS[pd$OS=="DEAD"] <- 1
table(pd$OS)
str(pd)

pd$Age <- as.numeric(as.character(pd$Age))
pd$OS <- as.numeric(as.character(pd$OS))
pd$OS_Time <- as.numeric(as.character(pd$OS_Time))

pd$OS_Time <- pd$OS_Time / 12
pd$OS_Time <- round(pd$OS_Time, 2)

str(pd)

pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]
pd <- subset(pd, pd$OS_Time>0)
str(pd)

table(pd$Gender)
mean(pd$Age)
table(pd$T_stages)
table(pd$N_stages)
table(pd$M_stages)
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
