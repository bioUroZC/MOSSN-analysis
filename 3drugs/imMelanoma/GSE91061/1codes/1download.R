#=======================================================

#=======================================================

rm(list=ls())
library(dplyr)
library(tibble)
library(tidyr)
library(GEOquery)
library(gdata)
library(data.table)
library(readxl)
options(stringsAsFactors = FALSE)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)
setwd("/proj/c.zihao/work1/3drugs/immune/GSE91061/data/")

#=======================================================

#=======================================================

gsename = "GSE91061"
gse<- getGEO(gsename, destdir = ".") 


express <- fread('GSE91061_raw_counts_GRCh38.p13_NCBI.tsv')
express <- as.data.frame(express)
str(express)



gpl <- fread('Human.GRCh38.p13.annot.tsv')
gpl <- as.data.frame(gpl)
colnames(gpl)
head(gpl)
gpl <- gpl %>% dplyr::select(GeneID, "GeneType",   "Symbol" )
head(gpl)
table(gpl$GeneType)
gpl <- subset(gpl, gpl$GeneType == 'protein-coding' )
gpl$GeneType <- NULL

express <- merge(gpl, express,by="GeneID")
express$GeneID <- NULL
str(express)
express[which(is.na(express),arr.ind = T)]<-0 
express[1:5,1:5]


exprSet <- aggregate(x = express[,2:ncol(express)],
                     by = list(express$Symbol),
                     FUN = max)
head(exprSet)[1:5,1:5]

names(exprSet)[1] <- 'ID'
rownames(exprSet) <- exprSet$ID
exprSet$ID <- NULL
exprSet[1:5,1:5]

min(exprSet)
max(exprSet)


#=======================================================

#=======================================================

pd <- pData(gse$GSE91061_series_matrix.txt.gz)
head(pd)
names(pd)

table(pd$'tissue:ch1')
table(pd$'visit (pre or on treatment):ch1')

pd <- subset(pd, select=c("geo_accession"  , 
                          'visit (pre or on treatment):ch1',
                          'response:ch1' ))

names(pd) <- c("Sample", 'group', "Response")
pd <- subset(pd, pd$group == "Pre")

str(pd)

table(pd$Response)

pd <- subset(pd, pd$Response != "UNK")

pd <- subset(pd, select=c("Sample", 'Response'))
pd <- na.omit(pd)

print(table(pd$Response))


pd$Response[pd$Response=="PRCR"] <- 3
pd$Response[pd$Response=="PD"] <- 1
pd$Response[pd$Response=="SD"] <- 2
print(table(pd$Response))

pd$Response <- as.numeric(as.character(pd$Response))

#=======================================================

#=======================================================

samplesname <- intersect(pd$Sample, colnames(exprSet))
samplesname <- unique(samplesname)
pd <- pd[which(pd$Sample %in% samplesname),]
exprSet <- exprSet[,which(colnames(exprSet) %in% samplesname)]
colnames(exprSet)

print(dim(exprSet))

#=======================================================

#=======================================================

exprSet[1:5,1:5]
write.csv(exprSet, file = "exprSet.csv")
write.csv(pd, file = "pd.csv")

min(exprSet)
max(exprSet)

