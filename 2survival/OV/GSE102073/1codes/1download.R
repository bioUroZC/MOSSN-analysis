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
library(biomaRt)
library(gdata)
library(stringr)
library(readxl)
library(org.Hs.eg.db)
library(AnnotationDbi)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)

setwd("/proj/c.zihao/work1/1survival/OV/GSE102073/data/")
set.seed(1234)

#=======================================================

#=======================================================

express <- read_excel("GSE102073_RNAseq_exp_dat.xlsx")
express[1:5,1:5]

express$Ensembl_gene_id <- sub("\\..*", "", express$ID_REF)

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = express$Ensembl_gene_id,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

express$gene <- gene_symbols

express$Ensembl_gene_id <- NULL
express$ID_REF <- NULL 

express <- express %>%
  dplyr::select(gene, everything())

express[1:5, 1:5]


#=======================================================

#=======================================================


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

exprSet <- exprSet + 8


min(exprSet)
max(exprSet)



#=======================================================

#=======================================================

gsename <- "GSE102073"
gse <- getGEO(gsename, destdir = ".")

pd <- pData(gse$GSE102073_series_matrix.txt.gz)
head(pd)
names(pd)


pd <- subset(pd, select=c("title", 
                          "age:ch1"  ,
                          'path.stage:ch1',
                          "vitalstatus:ch1"  ,
                          "os.mos..:ch1" 
))

head(pd)
colnames(pd) <- c("Sample", "Age", "Stage",  "OS", "OS_Time")
head(pd)

str(pd)

table(pd$Age)
pd$Age <- as.numeric(as.character(pd$Age))
table(pd$Age)


table(pd$OS)
pd$OS[pd$OS=="Not available"] <- NA
pd$OS[pd$OS=="Dead"] <- 1
pd$OS[pd$OS=="Alive"] <- 0
pd$OS <-  as.numeric(as.character(pd$OS))
table(pd$OS)


str(pd)
table(pd$OS_Time)
pd$OS_Time[pd$OS_Time=="Not available"] <- NA
pd$OS_Time <-  as.numeric(as.character(pd$OS_Time))
pd$OS_Time <- pd$OS_Time / 12
pd$OS_Time <- round(pd$OS_Time, 2)
str(pd)

pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]
pd <- subset(pd, pd$OS_Time>0)
head(pd)

str(pd)
table(pd$OS)
mean(pd$OS_Time)

head(pd)

pd$Sample <- gsub("_", "", pd$Sample)
colnames(exprSet) <- paste0('HGSC', colnames(exprSet))

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
