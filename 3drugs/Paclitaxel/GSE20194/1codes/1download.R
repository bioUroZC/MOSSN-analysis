#=======================================================

#=======================================================


rm(list=ls())
library(dplyr)
library(tibble)
library(tidyr)
library(GEOquery)
library(gdata)

options(stringsAsFactors = FALSE)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)
setwd("/proj/c.zihao/work1/3drugs/Paclitaxel/GSE20194/data")

#=======================================================

#=======================================================

gsename = "GSE20194"
gse<- getGEO(gsename, destdir = ".") 

gpl<- getGEO('GPL96', destdir = ".") 
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

exprSet <- as.data.frame(exprs(gse$GSE20194_series_matrix.txt.gz))
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

exprSet <- exprSet + 4

min(exprSet)
max(exprSet)


#=======================================================

#=======================================================

pd <- pData(gse$GSE20194_series_matrix.txt.gz)
names(pd)
head(pd)
str(pd)

table(pd$'treatment code:ch1')

pd <- subset(pd, `treatment code:ch1` %in% c("TFAC", "TFEC", 
                                             "TH/FAC", "TH/FEC", 
                                             "Tonly", "TXFAC"))
table(pd$`treatment code:ch1`)


pd <- subset(pd, select=c("geo_accession" ,
                          "age:ch1" ,
                          "bmngrd:ch1" ,
                          "pcr_vs_rd:ch1" ))

names(pd) <- c("Sample", "Age", "Grade", "Response")
head(pd)
str(pd)

pd <- na.omit(pd)

table(pd$Response)

print(pd$Response)

pd$Response[pd$Response=="RD"] <- 1
pd$Response[pd$Response=="pCR"] <- 4
print(table(pd$Response))

pd$Response <- as.numeric(as.character(pd$Response))
print(pd$Response)


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


