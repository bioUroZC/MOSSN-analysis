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


setwd("/proj/c.zihao/work1/1survival/BRCA/GSE45255/data")


#=======================================================

#=======================================================

gsename <- "GSE45255"
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

exprSet <- as.data.frame(exprs(gse$GSE45255_series_matrix.txt.gz))
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

min(exprSet)
max(exprSet)

#=======================================================

#=======================================================

pd <- pData(gse$GSE45255_series_matrix.txt.gz)
head(pd)
names(pd)

pd <- subset(pd, select=c("geo_accession", 
                          "patient age:ch1",
                          "dfs event (defined as any type of recurrence or death from breast cancer):ch1",
                         "dmfs event (defined as distant metastasis or death from breast cancer):ch1",
                           "dfs time:ch1",
                         "dmfs time:ch1"
                          ))

head(pd)
colnames(pd) <- c("Sample", "Age", "dfs_event", "dmfs_event", "t_dfs", "t_dmfs")
head(pd)

str(pd)

table(pd$Age)
pd$Age[pd$Age=="NA"] <- NA
pd$Age <- as.numeric(as.character(pd$Age))

table(pd$dfs_event)
pd$dfs_event[pd$dfs_event=="NA"] <- NA
pd$dfs_event <- as.numeric(pd$dfs_event)
table(pd$dfs_event)

table(pd$dmfs_event)
pd$dmfs_event[pd$dmfs_event=="NA"] <- NA
pd$dmfs_event <- as.numeric(pd$dmfs_event)
table(pd$dmfs_event)




table(pd$t_dfs)
pd$t_dfs[pd$t_dfs=="NA"] <- NA
pd$t_dfs <- as.numeric(as.character(pd$t_dfs))
table(pd$t_dfs)

table(pd$t_dmfs)
pd$t_dmfs[pd$t_dmfs=="NA"] <- NA
pd$t_dmfs <- as.numeric(pd$t_dmfs)
table(pd$t_dmfs)



head(pd)
pd$EFS <- with(pd, ifelse(dfs_event == 1 | dmfs_event == 1, 1, 0))


pd <- pd[!(is.na(pd$t_dfs) & is.na(pd$t_dmfs)), ]


pd$EFS_Time <- apply(
  pd[, c("t_dfs", "t_dmfs")], 
  1, 
  function(x) min(x, na.rm = TRUE)
)

head(pd)

pd <- subset(pd, select=c("Sample", "Age", "EFS", "EFS_Time" ))

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
