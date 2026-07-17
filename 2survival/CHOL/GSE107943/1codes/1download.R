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

setwd(paste0(PROJ_ROOT, "/1survival/CHOL/GSE107943/data/"))
set.seed(1234)

#=======================================================

#=======================================================

gsename <- "GSE107943"
gse <- getGEO(gsename, destdir = ".")

#=======================================================

#=======================================================


exprSet <- fread('GSE107943_RPKM.txt')
exprSet <- as.data.frame(exprSet)
exprSet[1:5,1:5]

exprSet$No <- NULL
exprSet$Chr <- NULL
exprSet$Ensenble <- NULL
exprSet$Start <- NULL
exprSet$Stop <- NULL
exprSet$CodingLength  <- NULL


new_colnames <- as.character(exprSet[1, ])
colnames(exprSet) <- new_colnames
exprSet <- exprSet[-1, ]
head(exprSet)


names(exprSet)[1] <- 'gene'
express <- exprSet
express <- express %>%
  mutate(across(-gene, as.numeric))

express[which(is.na(express),arr.ind = T)]<-0 
express[1:5,1:5]

str(express)

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

pd <- pData(gse$GSE107943_series_matrix.txt.gz)
head(pd)
names(pd)

table(pd$'tissue:ch1')

pd <- subset(pd, pd$`tissue:ch1` == 'Tumor (Intrahepatic cholangiocarcinoma)')

pd <- subset(pd, select=c( "title" , 
                           'age:ch1',
                           'Sex:ch1',
                           'stageajcc:ch1',
                           "death:ch1"  ,
                           "survival(mo):ch1" 
))

head(pd)
colnames(pd) <- c("Sample", "Age", "Gender",  "Stage",  "OS", "OS_Time")
head(pd)

str(pd)


table(pd$Age)
pd$Age <- as.numeric(as.character(pd$Age))
table(pd$Age)



table(pd$Gender)
pd$Gender[pd$Gender=="female"] <- "Female"
pd$Gender[pd$Gender=="male"] <- "Male"
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

pd$Sample <- paste0('GSE107943', pd$Sample)
colnames(exprSet) <- paste0('GSE107943', colnames(exprSet))

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
