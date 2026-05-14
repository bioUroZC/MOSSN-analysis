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
setwd("/proj/c.zihao/work1/3drugs/Paclitaxel/GSE241876/data")

#=======================================================

#=======================================================


gsename = "GSE241876"
gse<- getGEO(gsename, destdir = ".") 


express <- fread('GSE241876_raw_counts_GRCh38.p13_NCBI.tsv')
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

pd <- pData(gse$GSE241876_series_matrix.txt.gz)
head(pd)
names(pd)

table(pd$'tissue:ch1')
pd <- pd[!grepl("post", pd$`tissue:ch1`, ignore.case = TRUE), ]
table(pd$'tissue:ch1')

pd <- subset(pd, select=c("geo_accession"  , 
                          'overall response:ch1',
                          'pfs event:ch1',
                          "pfs:ch1"))

names(pd) <- c("Sample", "Response", "PFS", "PFS_Time")
str(pd)


table(pd$PFS)
pd$PFS <- as.numeric(as.character(pd$PFS))
table(pd$PFS)


table(pd$PFS_Time)
pd$PFS_Time <- as.numeric(as.character(pd$PFS_Time))
pd$PFS_Time <- pd$PFS_Time / 12
pd$PFS_Time <- round(pd$PFS_Time, 2)
table(pd$PFS_Time)


print(table(pd$Response))

pd$Response[pd$Response=='Progressive Disease'] <- 1
pd$Response[pd$Response=='Stable Disease'] <- 2

pd$Response[pd$Response=='Complete Response'] <- 4
pd$Response[pd$Response=='Partial Response'] <- 3
table(pd$Response)

pd <- subset(pd, select=c("Sample", "Response"))


pd$Response <- as.numeric(as.character(pd$Response))
print(table(pd$Response))

pd <- na.omit(pd)

print(pd$Response)

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


min(exprSet)
max(exprSet)


