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
setwd("/proj/c.zihao/work1/3drugs/Paclitaxel/GSE41998/data")

#=======================================================

#=======================================================


gsename = "GSE41998"
gse<- getGEO(gsename, destdir = ".") 

gpl<- getGEO('GPL571', destdir = ".") 
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

exprSet <- as.data.frame(exprs(gse$GSE41998_series_matrix.txt.gz))
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

pd <- pData(gse$GSE41998_series_matrix.txt.gz)
names(pd)
head(pd)


table(pd$`treatment arm:ch1`)
pd <- subset(pd, pd$`treatment arm:ch1` == 'Paclitaxel')

pd <- subset(pd, select=c("geo_accession" , 'ac response:ch1'))
names(pd) <- c("Sample", "Response")
table(pd$Reseponse)

str(pd)



print(table(pd$Response))
pd$Response[pd$Response=='progressive disease'] <- 1
pd$Response[pd$Response=='stable disease'] <- 2

pd$Response[pd$Response=='partial response'] <- 3
pd$Response[pd$Response=='complete response'] <- 4

pd$Response[pd$Response=='unable to determine'] <- NA
pd <- na.omit(pd)

print(table(pd$Response))
pd$Response <- as.numeric(as.character(pd$Response))
print(pd$Response)

pd$Sample[1:10]
colnames(exprSet)[1:10]

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


