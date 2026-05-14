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
library(ggplot2)
library(biomaRt)
library(RColorBrewer)
library(gdata)
set.seed(1234)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)

setwd("/proj/c.zihao/work1/1survival/LUAD/GSE68465/data/")

#=======================================================

#=======================================================

gsename <- "GSE68465"
gse <- getGEO(gsename, destdir = ".")
gpl <- getGEO('GPL96', destdir = ".")
colnames(Table(gpl))
Table(gpl)[1:10, 1:6]

gpl <- gpl@dataTable@table
colnames(gpl)
gpl <- gpl %>% dplyr::select(ID, "Gene Symbol")
write.csv(gpl, file = 'gpl.csv')

gpl <- read.csv('gpl.csv', header = TRUE, row.names = 1)
head(gpl)

gpl <- gpl %>%
  tidyr::separate("Gene.Symbol", c('gene', 'symbol'), sep = '\\///', fill = "right") %>%
  dplyr::select("ID", 'gene')

gpl$gene <- gdata::trim(gpl$gene)

#=======================================================

#=======================================================

exprSet <- as.data.frame(exprs(gse[["GSE68465_series_matrix.txt.gz"]]))
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

gene_iqr <- apply(exprSet, 1, IQR)
mean_iqr <- mean(gene_iqr)
median_iqr <- median(gene_iqr)
mean_iqr
median_iqr

min(exprSet)
max(exprSet)

#=======================================================

#=======================================================

pd <- pData(gse[["GSE68465_series_matrix.txt.gz"]])
table(pd$characteristics_ch1)

pd <- subset(pd, pd$characteristics_ch1 == "disease_state: Lung Adenocarcinoma")

head(pd)

pd <- subset(pd, select=c("geo_accession", "Sex:ch1", "age:ch1",
                          "disease_stage:ch1", 
                          "vital_status:ch1", "months_to_last_contact_or_death:ch1"))

names(pd) <- c("Sample", 'Gender', "Age", "Stage", "OS", "OS_Time")
head(pd)


table(pd$Gender)
table(pd$Stage)


pd$Stage    <- gsub("disease_stage: ", "", pd$Stage)

table(pd$OS)
pd$OS[pd$OS=="Alive"] <- 0
pd$OS[pd$OS=="Dead"] <- 1
table(pd$OS)
str(pd)

pd$OS  <-  as.numeric(as.character(pd$OS))
pd$OS_Time <-  as.numeric(as.character(pd$OS_Time))
pd$Age  <-  as.numeric(as.character(pd$Age))

pd$OS_Time <- pd$OS_Time / 12
pd$OS_Time <- round(pd$OS_Time, 2)

pd$T_stages <- gsub(".*pT(\\d+).*", "\\1", pd$Stage)
pd$T_stages <- as.numeric(as.character(pd$T_stages))

pd$N_stages <- gsub(".*pN(\\d+).*", "\\1", pd$Stage)
pd$N_stages <- as.numeric(as.character(pd$N_stages))

str(pd)

pd$Stage <- NULL

pd$T_stages <- paste0("T", pd$T_stages)
pd$N_stages <- paste0("N", pd$N_stages)

pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]
pd <- subset(pd, pd$OS_Time>0)
str(pd)
table(pd$Gender)
mean(pd$Age)
table(pd$T_stages)
table(pd$N_stages)
table(pd$M_stages)
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
