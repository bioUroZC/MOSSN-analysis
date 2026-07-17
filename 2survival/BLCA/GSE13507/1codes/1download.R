# ===============================================================

# ===============================================================

rm(list=ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")

setwd(paste0(PROJ_ROOT, "/1survival/BLCA/GSE13507/data/"))
library(dplyr)
library(tibble)
library(tidyr)
library(GEOquery)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)

# ===============================================================
#
# ===============================================================

gse<- getGEO("GSE13507", destdir = ".") 
gpl<- getGEO('GPL6102', destdir = ".") 

colnames(Table(gpl))
Table(gpl)[1:10,1:6] 
gpl <- gpl@dataTable@table
colnames(gpl)
gpl <- gpl %>% dplyr::select(ID, "Symbol")
write.csv(gpl,"GPL.csv", row.names = F)
genename = read.csv("GPL.csv")

# ===============================================================
#
# ===============================================================

exprSet <- as.data.frame(exprs(gse$GSE13507_series_matrix.txt.gz)) 
exprSet$ID = rownames(exprSet)
express = merge( x=genename, y=exprSet, by="ID")
express$ID = NULL
express[which(is.na(express),arr.ind = T)]<-0 
express[1:4,1:4]

exprSet <- aggregate(x = express[,2:ncol(express)],
                     by = list(express$Symbol),
                     FUN = max)
head(exprSet)
exprSet <- as.data.frame(exprSet)
exprSet[1:4,1:4]

exprSet <-exprSet[-1,]
names(exprSet)[1] <- 'ID'
rownames(exprSet) <- exprSet$ID
exprSet$ID <- NULL

min(exprSet)
max(exprSet)

# ===============================================================
#
# ===============================================================

pd <- pData(gse$GSE13507_series_matrix.txt.gz)
pd <- pd[grepl("Primary bladder cancer", pd$title), ]
head(pd)

pd <- subset(pd, select=c( "geo_accession" ,'characteristics_ch1.3',
                           "AGE:ch1",  "SEX:ch1",   "stage:ch1" , 
                         "overall survival:ch1" ,  "survival month:ch1"  ))
names(pd) <- c("Sample", 'Group', "Age", "Gender", "Stage", "OS", "OS_Time")
table(pd$OS)
str(pd)

pd$Age <- as.numeric(as.character(pd$Age))
table(pd$Gender)
pd$Gender <- ifelse(pd$Gender == 'F', "Female", "Male")
table(pd$Gender)

table(pd$OS)
pd$OS <- ifelse(pd$OS == 'survival', 0, 1)
pd$OS <- as.numeric(as.character(pd$OS))
table(pd$OS)


pd$OS_Time <- as.numeric(as.character(pd$OS_Time))
pd$OS_Time <- pd$OS_Time / 12
pd$OS_Time <- round(pd$OS_Time, 2)
pd$OS_Time[1:5]
str(pd)


pd$T_stages <- gsub(".*T(\\d+).*", "\\1", pd$Stage)
pd$T_stages <- as.numeric(as.character(pd$T_stages))

pd$N_stages <- gsub(".*N(\\d+).*", "\\1", pd$Stage)
pd$N_stages <- as.numeric(as.character(pd$N_stages))

pd$M_stages <- gsub(".*M(\\d+).*", "\\1", pd$Stage)
pd$M_stages <- as.numeric(as.character(pd$M_stages))

pd$T_stages <- paste0("T", pd$T_stages)
pd$N_stages <- paste0("N", pd$N_stages)
pd$M_stages <- paste0("M", pd$M_stages)

str(pd)

pd$Stage <- NULL

pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]
pd <- subset(pd, pd$OS_Time > 0)
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

