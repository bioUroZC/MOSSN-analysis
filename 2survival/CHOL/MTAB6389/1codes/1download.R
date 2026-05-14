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

setwd("/proj/c.zihao/work1/1survival/CHOL/MTAB6389/data/")
set.seed(1234)

#=======================================================

#=======================================================

# download.file(
#   url = "ftp://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/389/E-MTAB-6389/Files/data_exp_icc.txt",
#   destfile = "data_exp_icc.txt",
#   mode = "wb",
#   method = "curl"
# )

gpl <- getGEO('GPL17585', destdir = ".")
colnames(Table(gpl))
Table(gpl)[1:10, 1:6]

gpl <- gpl@dataTable@table
colnames(gpl)


gpl <- gpl %>% dplyr::select(ID, "gene_symbols")
write.csv(gpl, file = 'gpl.csv')
gpl <- read.csv('gpl.csv', header = TRUE, row.names = 1)
head(gpl)
names(gpl)[2] <- 'gene'
head(gpl)



#=======================================================

#=======================================================

exprSet <- fread('data_exp_icc.txt')
exprSet <- as.data.frame(exprSet)
exprSet[1:5,1:5]
names(exprSet)[1] <- "ID"
exprSet[1:5,1:5]


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

pd <- fread('https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/389/E-MTAB-6389/Files/E-MTAB-6389.sdrf.txt')
pd <- as.data.frame(pd)
colnames(pd)

table(pd$"Characteristics[sampling site]")

pd <- subset(pd, pd$`Characteristics[sampling site]` == 'tumor tissue')


pd <- subset(pd, select=c( "Source Name" ,
                           "Characteristics[sex]" ,
                           "Characteristics[event death]"  ,
                           "Characteristics[overall survival]" 
))

head(pd)
colnames(pd) <- c("Sample", "Gender",  "OS", "OS_Time")
head(pd)

str(pd)

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



