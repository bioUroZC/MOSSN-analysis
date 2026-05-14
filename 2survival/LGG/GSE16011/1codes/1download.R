#=======================================================

#=======================================================

rm(list = ls())

library(GEOquery)
library(dplyr)
library(tidyr)
library(Biobase)
library(limma)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tibble)
library(ggplot2)
library(biomaRt)
library(RColorBrewer)
library(gdata)
set.seed(1234)
library(readxl)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)

setwd("/proj/c.zihao/work1/1survival/LGG/GSE16011/data/")

#=======================================================

#=======================================================

gsename <- "GSE16011"
gse <- getGEO(gsename, destdir = ".")
gpl <- getGEO('GPL8542', destdir = ".")
colnames(Table(gpl))
Table(gpl)[1:10, 1:6]

gpl <- gpl@dataTable@table
colnames(gpl)
gpl <- gpl %>% dplyr::select(ID, "ORF")
write.csv(gpl, file = 'gpl.csv')

gpl <- read.csv('gpl.csv', header = TRUE, row.names = 1)
head(gpl)


gpl$ORF <- as.character(gpl$ORF)
e2s <- bitr(gpl$ORF,
            fromType = "ENTREZID",
            toType = "SYMBOL",
            OrgDb = org.Hs.eg.db)

gpl_anno <- merge(gpl, e2s, by.x = "ORF", by.y = "ENTREZID", all.x = TRUE)
gpl_anno <- gpl_anno[, c("ID", "SYMBOL")]
head(gpl_anno)
gpl <- gpl_anno
gpl$gene <- gdata::trim(gpl$SYMBOL)
gpl$SYMBOL <- NULL
head(gpl)

#=======================================================

#=======================================================

exprSet <- as.data.frame(exprs(gse$GSE16011_series_matrix.txt.gz))
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

pd <- pData(gse$GSE16011_series_matrix.txt.gz)
head(pd)
names(pd)

table(pd$`histology:ch1`)




pd <- subset(pd, select = c('title', 'geo_accession', 
                            "age at diagnosis:ch1" ,
                            'histology:ch1'))

head(pd)


LGG_types <- c(
  "A (grade II)", "A (grade III)",
  "OA (grade II)", "OA (grade III)",
  "OD (grade II)", "OD (grade III)"
)

pd <- pd[pd$`histology:ch1` %in% LGG_types, ]

stabs_1_6 <- read_excel("stabs_1-6.xlsx")
stabs_1_6 <- subset(stabs_1_6, select = c("Database number", "Gender", "Alive", "Survival (years)"))
stabs_1_6$title <- paste0('glioma ', stabs_1_6$`Database number`)
pd <- merge(pd, stabs_1_6, by = 'title')

names(pd) <- c("title", "Sample", "Age", "Grade", "number", "Gender", "OS", "OS_Time" )
pd$title <- NULL

pd$Age <- as.numeric(as.character(pd$Age))
pd$Age <- round(pd$Age, 0)

table(pd$Grade)
pd$Grade <- NULL
pd$number <- NULL

table(pd$Gender)
pd$Gender  <- tools::toTitleCase(tolower(gsub(" ", "", pd$Gender)))

str(pd)

table(pd$OS)
pd$OS <- ifelse(pd$OS == 'Dead', 1, 0)
pd$OS <- as.numeric(as.character(pd$OS))
table(pd$OS)


pd$OS_Time <- chartr(old = ',', new = '.', x = pd$OS_Time)
table(pd$OS_Time)
pd$OS_Time <- as.numeric(as.character(pd$OS_Time))
table(pd$OS_Time)

pd <- subset(pd, pd$OS_Time > 0)

head(pd)
str(pd)

pd <- subset(pd, select=c( "Sample", "Gender","Age","OS", "OS_Time"))

head(pd)

str(pd)

str(pd)
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
