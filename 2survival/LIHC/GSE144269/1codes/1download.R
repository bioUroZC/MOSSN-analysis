#=======================================================

#=======================================================

rm(list=ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")

library(dplyr)
library(tibble)
library(tidyr)
library(GEOquery)
library(gdata)
library(stringr)
library(openxlsx)

options(stringsAsFactors = FALSE)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)
setwd(paste0(PROJ_ROOT, "/1survival/LIHC/GSE144269/data/"))

#=======================================================

#=======================================================

express <- read.table("GSE144269_RSEM_GeneCounts.txt", header = T)

express[1:5,1:5]

# Assuming 'express' is your data frame
express <- express %>%
  tidyr::separate(entrez_id, into = c("ensembl_id", "gene_symbol"), sep = "\\|") %>%
  dplyr::select(-ensembl_id) %>%
  dplyr::relocate(gene_symbol)

express[1:5,1:5]

express[which(is.na(express),arr.ind = T)]<-0 
express[1:4,1:4]

exprSet <- aggregate(x = express[,2:ncol(express)],
                     by = list(express$gene_symbol),
                     FUN = max)
head(exprSet)
exprSet <- as.data.frame(exprSet)
exprSet[1:4,1:4]

names(exprSet)[1] <- 'ID'
rownames(exprSet) <- exprSet$ID
exprSet$ID <- NULL
exprSet[1:4,1:4]


min(exprSet)
max(exprSet)

#=======================================================

#=======================================================

gse<- getGEO("GSE144269", destdir = ".") 

annodata <- pData(gse$GSE144269_series_matrix.txt.gz)

annodata <- subset(annodata, select=c("title" ,  "tumor/non-tumor:ch1"  ))
names(annodata)[2] <- 'group'
table(annodata$group)
annodata <- subset(annodata, annodata$group == "tumor")

annodata$title[1:5]
annodata <- annodata %>%
  tidyr::separate(title, into = c("ID", "Sample"), sep = " \\[|\\]", remove = TRUE)

annodata$ID[1:5]
annodata$ID <- str_extract(annodata$ID, "\\d+")
annodata$ID <- paste0('Patient', annodata$ID)

#=======================================================

#=======================================================

pd <- read.xlsx("41467_2020_18186_MOESM4_ESM.xlsx", sheet = 1)
colnames(pd)

pd <- subset(pd, select=c("PatientID", "TNMstaging", "Status", "Time"  ))

names(pd) <- c("ID", "Stage", "OS", "OS_Time")

pd$ID <- paste0("Patient", pd$ID)

str(pd)

table(pd$Stage)
pd$Stage[pd$Stage==1] <-"I"
pd$Stage[pd$Stage==2] <-"II"
pd$Stage[pd$Stage==3] <- "III"
pd$Stage[pd$Stage==4] <- "IV"
table(pd$Stage)

str(pd)

table(pd$OS)

pd$OS[pd$OS=="Alive"] <- 0
pd$OS[pd$OS=="Deceased"] <- 1
pd$OS <- as.numeric(as.character(pd$OS))
str(pd)



pd$OS_Time <- pd$OS_Time / 365
pd$OS_Time <- round(pd$OS_Time, 2)
pd <- subset(pd, pd$OS_Time > 0)

pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]

table(pd$OS)
mean(pd$OS_Time)

pd <- merge(pd, annodata, by="ID")
colnames(pd)

pd<- pd %>%
  dplyr::select("Sample",   "Stage" ,  "OS" ,     "OS_Time")

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

