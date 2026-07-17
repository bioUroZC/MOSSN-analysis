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
library(ggplot2)
library(biomaRt)
library(RColorBrewer)
library(gdata)
library(readxl)
set.seed(1234)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)

setwd(paste0(PROJ_ROOT, "/1survival/STAD/GSE57303/data"))

#=======================================================

#=======================================================

gsename <- "GSE57303"
gse <- getGEO(gsename, destdir = ".")
gpl <- getGEO('GPL570', destdir = ".")
colnames(Table(gpl))
Table(gpl)[1:10, 1:6]


gpl <- gpl@dataTable@table
colnames(gpl)
gpl <- gpl %>% dplyr::select(ID, "Gene Symbol")
write.csv(gpl, file = 'gpl.csv')

gpl <- read.csv('gpl.csv', header = TRUE, row.names = 1)
head(gpl)

gpl <- gpl %>%
  tidyr::separate("Gene.Symbol", c('gene', 'symbol'), sep = '\\///') %>%
  dplyr::select("ID", 'gene')

gpl$gene <- gdata::trim(gpl$gene)

#=======================================================

#=======================================================


exprSet <- as.data.frame(exprs(gse$GSE57303_series_matrix.txt.gz))
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


pd <- pData(gse[["GSE57303_series_matrix.txt.gz"]])
head(pd)
names(pd)

pd <- pd %>%
  dplyr::select( title, geo_accession, 
                 "phenotype:ch1" )

table(pd$`phenotype:ch1`)
pd$title <- gsub("Human_BTH_", '', pd$title)
str(pd)

library(readxl)
clinical <- read_excel("gcc22196-sup-0006-supptable1.xlsx", sheet = "Sheet1")
names(clinical)

clinical$Gene
pd$title

intersect(clinical$Gene, pd$title)



clinical <- subset(clinical, select=c("Gene", "Age" ,  "Gender"  , "Clinical Stage"  ,
                                      "TNM"  ,
                                      "Survival State" ,
                                      "Survival Length" ))

names(clinical)[1] <- 'title'



pd <- merge(pd, clinical, by='title')
names(pd)
pd <- subset(pd, select=c( "geo_accession", "Age" ,  "Gender"  ,   "Clinical Stage"   ,
                           "TNM" ,   "Survival State" , "Survival Length" ))

names(pd) <- c("Sample", "Age", "Gender", "Stage", "TNM", "OS", "OS_Time")
str(pd)

table(pd$Gender)


table(pd$Stage)
pd$Stage[pd$Stage=="IB"] <- "I"
pd$Stage[pd$Stage=="IIIA"] <- "III"
pd$Stage[pd$Stage=="IIIB"] <- "III"
table(pd$Stage)


str(pd)     

str(pd)
table(pd$OS)
pd$OS[pd$OS=="dead"] <- 1
pd$OS[pd$OS=="survival"] <- 0
pd$OS[pd$OS=="Survival"] <- 0
pd$OS <-  as.numeric(as.character(pd$OS))
table(pd$OS)

str(pd)
table(pd$OS_Time)
pd$OS_Time <- pd$OS_Time / 12
pd$OS_Time <- round(pd$OS_Time, 2)

pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]
pd <- subset(pd, pd$OS_Time>0)
head(pd)

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


