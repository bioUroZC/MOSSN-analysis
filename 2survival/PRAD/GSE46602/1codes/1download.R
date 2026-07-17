#=======================================================

#=======================================================


rm(list=ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")

library(dplyr)
library(tibble)
library(tidyr)
library(GEOquery)
library(gdata)

options(stringsAsFactors = FALSE)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)
setwd(paste0(PROJ_ROOT, "/1survival/PRAD/GSE46602/data"))

#=======================================================

#=======================================================

gsename = "GSE46602"
gse<- getGEO(gsename, destdir = ".") 

gpl<- getGEO('GPL570', destdir = ".") 
colnames(Table(gpl))
Table(gpl)[1:10,1:6]
gpl <- gpl@dataTable@table
colnames(gpl)
gpl <- gpl %>% dplyr::select(ID,  "Gene Symbol" )


write.csv(gpl, file = 'gpl.csv')
gpl <- read.csv('gpl.csv', header = T, row.names = 1)
head(gpl)

gpl <- gpl %>%
  tidyr::separate(Gene.Symbol, c('gene', 'symbol'), sep='\\///'   )%>%
  dplyr::select("ID", 'gene')

gpl$gene <- gdata::trim(gpl$gene)

#=======================================================

#=======================================================

exprSet <- as.data.frame(exprs(gse$GSE46602_series_matrix.txt.gz))
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

pd <- pData(gse$GSE46602_series_matrix.txt.gz)
names(pd)
head(pd)

table(pd$`tissue:ch1`)

pd <- subset(pd, pd$`tissue:ch1` == 'prostate tumor')

pd <- subset(pd, select=c("geo_accession" , 'age:ch1', 't-stage:ch1',
                          "bcr:ch1", 'bcr_free_time:ch1'))
names(pd) <- c("Sample", "Age", "Stage", "OS", "OS_Time")
head(pd)
str(pd)

pd$Age <- as.numeric(as.character(pd$Age))


table(pd$OS)
pd$OS[pd$OS=="NO"] <- 0
pd$OS[pd$OS=="YES"] <- 1
pd$OS <- as.numeric(as.character(pd$OS))
table(pd$OS)

table(pd$OS_Time)
pd$OS_Time <- as.numeric(as.character(pd$OS_Time))
pd$OS_Time <- pd$OS_Time / 12
pd$OS_Time <- round(pd$OS_Time, 2)
table(pd$OS_Time)

pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]
pd <- subset(pd, pd$OS_Time > 0)

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


