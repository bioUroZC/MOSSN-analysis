
#=======================================================

#=======================================================

rm(list = ls())

library(stringr)
library(GEOquery)
library(dplyr)
library(tidyr)
library(Biobase)
library(limma)
library(data.table)
library(tibble)
library(ggplot2)
library(biomaRt)
library(gdata)
set.seed(1234)
library(org.Hs.eg.db)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)

setwd("/proj/c.zihao/work1/3drugs/immune/PRJEB23709/data/")

#=======================================================

#=======================================================


express <- read.delim("ICB.Gide2019_Pembrolizumab-Nivolumab+Ipilimumab_Melanoma.self_subtract",
                      header=TRUE)
express$gene_id <- rownames(express)
express[1:4,1:4]



gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = as.character(express$gene_id), 
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)
head(gene_symbols)


express$gene <- gene_symbols
express$gene_id <- NULL
express <- express %>%
  dplyr::select(gene, everything())

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

exprSet <- exprSet + 10

min(exprSet)
max(exprSet)

#=======================================================

#=======================================================


pd <- read.delim('ICB.Gide2019_Pembrolizumab-Nivolumab+Ipilimumab_Melanoma.clinical', header = T)
head(pd)
table(pd$Response)
names(pd)[1] <- "Sample"
pd$Sample <- paste0('X', pd$Sample)
str(pd)


pd <- subset(pd, select=c("Sample", "Age", "Gender", 
                          "Site", "Response", "PFS.Event", "PFS"))
names(pd) <- c("Sample", "Age", "Gender", 
               "Site", "Response", "PFS", "PFS_Time")
str(pd)

pd$PFS_Time <- pd$PFS_Time / 365

pd$Sample
colnames(exprSet)

pd$Sample <- paste0('PNI', pd$Sample)
colnames(exprSet) <- paste0('PNI', colnames(exprSet))


pd <- subset(pd, select=c("Sample", "Response"))

table(pd$Response)
str(pd)

pd$Response[pd$Response==1] <- "Res"
pd$Response[pd$Response==0] <- "Non"


table(pd$Response)
pd <- na.omit(pd)



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
write.csv(exprSet, file = "exprSet2.csv")
write.csv(pd, file = "pd2.csv")
