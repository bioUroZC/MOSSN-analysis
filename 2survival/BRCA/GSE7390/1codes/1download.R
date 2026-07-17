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

setwd(paste0(PROJ_ROOT, "/1survival/BRCA/GSE7390/data"))

#=======================================================

#=======================================================

gsename <- "GSE7390"
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

exprSet <- as.data.frame(exprs(gse$GSE7390_series_matrix.txt.gz))
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

exprSet <- exprSet + 3

#=======================================================

#=======================================================

pd <- pData(gse$GSE7390_series_matrix.txt.gz)
head(pd)
names(pd)

pd <- subset(pd, select=c("geo_accession", 
                          "e.dmfs:ch1",
                          "t.dmfs:ch1"))

head(pd)
names(pd) <- c("Sample", "EFS", "EFS_Time")
head(pd)

str(pd)

table(pd$EFS)
pd$EFS  <-  as.numeric(as.character(pd$EFS))
str(pd)

table(pd$EFS_Time)
table(pd$EFS_Time)
pd$EFS_Time <-  as.numeric(as.character(pd$EFS_Time))
pd$EFS_Time <- pd$EFS_Time / 365
pd$EFS_Time <- round(pd$EFS_Time, 2)
str(pd)

pd <- pd[!is.na(pd$EFS) & !is.na(pd$EFS_Time), ]
pd <- subset(pd, pd$EFS_Time>0)
head(pd)

str(pd)
table(pd$EFS)
mean(pd$EFS_Time)

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
