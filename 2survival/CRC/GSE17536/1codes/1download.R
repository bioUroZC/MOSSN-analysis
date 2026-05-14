# ===============================================================

# ===============================================================

rm(list=ls())
setwd("/proj/c.zihao/work1/1survival/CRC/GSE17536/data/")
library(dplyr)
library(tibble)
library(tidyr)
library(GEOquery)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)

# ===============================================================
#
# ===============================================================

gse<- getGEO("GSE17536", destdir = ".") 
gpl<- getGEO('GPL570', destdir = ".") 

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


# ===============================================================
#
# ===============================================================

exprSet <- as.data.frame(exprs(gse$GSE17536_series_matrix.txt.gz)) 
exprSet$ID = rownames(exprSet)
express = merge(x=gpl, y=exprSet, by="ID")
express$ID = NULL
express[which(is.na(express),arr.ind = T)]<-0 
express[1:4,1:4]

exprSet <- aggregate(x = express[,2:ncol(express)],
                     by = list(express$gene),
                     FUN = max)
head(exprSet)
exprSet <- as.data.frame(exprSet)
exprSet[1:4,1:4]

exprSet <-exprSet[-1,]
names(exprSet)[1] <- 'ID'
rownames(exprSet) <- exprSet$ID
exprSet$ID <- NULL

min(exprSet)
max(exprSet)

# ===============================================================
#
# ===============================================================

pd <- pData(gse$GSE17536_series_matrix.txt.gz)
head(pd)

pd <- subset(pd, select=c( "geo_accession" ,
                           "age:ch1",  "gender:ch1",   "ajcc_stage:ch1" , 
                         "overall_event (death from any cause):ch1" ,  "overall survival follow-up time:ch1"  ))
names(pd) <- c("Sample",  "Age", "Gender", "Stage", "OS", "OS_Time")
table(pd$OS)
str(pd)

pd$Age <- as.numeric(as.character(pd$Age))
table(pd$Gender)
pd$Gender <- ifelse(pd$Gender == 'female', "Female", "Male")
table(pd$Gender)

table(pd$OS)
pd$OS[pd$OS=="no death"] <- 0
pd$OS[pd$OS=="death"] <- 1
table(pd$OS)
pd$OS <- as.numeric(as.character(pd$OS))
table(pd$OS)


pd$OS_Time <- as.numeric(as.character(pd$OS_Time))
pd$OS_Time <- pd$OS_Time / 12
pd$OS_Time <- round(pd$OS_Time, 2)
str(pd)

table(pd$Stage)

pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]
pd <- subset(pd, pd$OS_Time > 0)
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

