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
library(biomaRt)
library(gdata)
library(stringr)
library(readxl)
library(org.Hs.eg.db)
library(AnnotationDbi)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)

setwd(paste0(PROJ_ROOT, "/1survival/CHOL/GSE244807/data/"))
set.seed(1234)

#=======================================================

#=======================================================

gsename <- "GSE244807"
gse <- getGEO(gsename, destdir = ".")

#=======================================================

#=======================================================


exprSet <- fread("GSE244807_raw_counts.txt")
exprSet <- as.data.frame(exprSet)
exprSet[1:5,1:5]
names(exprSet)[1] <- "ensembl_gene_id"
exprSet[1:5,1:5]

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_ids <- exprSet$ensembl_gene_id
geneAnno <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = gene_ids,
  mart = mart
)
head(geneAnno)

express <- merge(geneAnno, exprSet, by = "ensembl_gene_id")
express[1:5,1:5]

express$ensembl_gene_id <- NULL
names(express)[1] <- 'gene'
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

pd <- pData(gse$GSE244807_series_matrix.txt.gz)
head(pd)
names(pd)

table(pd$'tissue:ch1')

pd <- subset(pd, select=c( "sample id:ch1",
                           'age at_the_time_of_sampling:ch1',
                           'gender:ch1',
                           "death:ch1"  ,
                           "os:ch1" 
))

head(pd)
colnames(pd) <- c("Sample", "Age", "Gender",  "OS", "OS_Time")
head(pd)

str(pd)

table(pd$Age)
pd$Age <- as.numeric(as.character(pd$Age))
table(pd$Age)

table(pd$Gender)

table(pd$OS)
pd$OS[pd$OS=="1"] <- 1
pd$OS[pd$OS=="0"] <- 0
pd$OS <-  as.numeric(as.character(pd$OS))
table(pd$OS)


str(pd)
table(pd$OS_Time)
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
