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

setwd(paste0(PROJ_ROOT, "/1survival/BRCA/GSE20711/data"))

#=======================================================

#=======================================================

gsename <- "GSE20711"
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

exprSet <- as.data.frame(exprs(gse$GSE20711_series_matrix.txt.gz))
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

pd <- pData(gse$GSE20711_series_matrix.txt.gz)
head(pd)
names(pd)

pd <- subset(pd, select=c("geo_accession", 
                          "e.os:ch1",
                          "e.rfs:ch1",
                          "t.os:ch1",
                          "t.rfs:ch1"))

head(pd)
str(pd)

colnames(pd) <- c("geo_accession", "os_event", "rfs_event", "t_os", "t_rfs")


table(pd$t_os)
pd$t_os[pd$t_os=="NA"] <- NA
pd$t_os <- as.numeric(sub(" y", "", pd$t_os))
table(pd$t_os)

table(pd$t_rfs)
pd$t_rfs[pd$t_rfs=="NA"] <- NA
pd$t_rfs <- as.numeric(pd$t_rfs)
table(pd$t_rfs)

table(pd$os_event)
pd$os_event[pd$os_event=="NA"] <- NA
pd$os_event <- as.numeric(pd$os_event)
table(pd$os_event)

table(pd$rfs_event)
pd$rfs_event[pd$rfs_event=="NA"] <- NA
pd$rfs_event <- as.numeric(pd$rfs_event)
table(pd$rfs_event)

head(pd)
pd$EFS <- with(pd, ifelse(os_event == 1 | rfs_event == 1, 1, 0))


pd <- pd[!(is.na(pd$t_os) & is.na(pd$t_rfs)), ]


pd$EFS_Time <- apply(
  pd[, c("t_os", "t_rfs")], 
  1, 
  function(x) min(x, na.rm = TRUE)
)

head(pd)

pd <- subset(pd, select=c("geo_accession", "EFS", "EFS_Time" ))
names(pd) <- c("Sample", "EFS", "EFS_Time" )

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
