#===================================================

#===================================================

rm(list = ls())

library(data.table)
library(GEOquery)
library(dplyr)
library(tidyr)
library(Biobase)
library(limma)
library(tibble)
library(stringr)
library(gdata)
Sys.setenv("VROOM_CONNECTION_SIZE"=131072*16)

#===================================================

#===================================================

setwd("/proj/c.zihao/work1/1survival/KIRC/GSE167573/data/")

gse <- getGEO("GSE167573", destdir = ".") 


exprSet <- fread("GSE167573_Processed_normalized_count_matrix.txt")
exprSet <- as.data.frame(exprSet)
exprSet[1:4,1:4]
rownames(exprSet) <- exprSet$GeneSymbol
exprSet$GeneSymbol <- NULL
exprSet[1:4,1:4]

min(exprSet)
max(exprSet)

#===================================================

#===================================================

pd <- pData(gse$GSE167573_series_matrix.txt.gz)
names(pd)
head(pd)
table(pd$'tissue:ch1')
pd <- subset(pd, pd$'tissue:ch1' == "Renal cell carcinoma")



pd <- subset(pd, select = c("title", 
                            "age(0:ch1",
                            "genotype(1:ch1",
                            "t stage (1:ch1",
                            "n stage (0:ch1",
                            "m stage (0:ch1",
                            "OS status:ch1" , 
                            "os time (months):ch1" ))

names(pd) <- c("Sample", "Age", "Gender", "Tstage", "Nstage", "Mstage",  "OS", "OS_Time")

str(pd)

table(pd$Age)
pd$Age <- as.integer(gsub(".*\\):\\s*", "", pd$Age))
table(pd$Age)

table(pd$Gender)
pd$Gender[pd$Gender=="male; 2: female): 1"] <- "Male"
pd$Gender[pd$Gender=="male; 2: female): 2"] <- "Female"
table(pd$Gender)

table(pd$Tstage)
pd$Tstage[pd$Tstage=="t1; 2: t2; 3: t3; 4:t4): 1"] <- "T1"
pd$Tstage[pd$Tstage=="t1; 2: t2; 3: t3; 4:t4): 2"] <- "T2"
pd$Tstage[pd$Tstage=="t1; 2: t2; 3: t3; 4:t4): 3"] <- "T3"
pd$Tstage[pd$Tstage=="t1; 2: t2; 3: t3; 4:t4): 4"] <- "T4"
table(pd$Tstage)


table(pd$Nstage)
pd$Nstage[pd$Nstage=="n0; 1:n1): 0"] <- "N0"
pd$Nstage[pd$Nstage=="n0; 1:n1): 1"] <- "N1"
table(pd$Nstage)


table(pd$Mstage)
pd$Mstage[pd$Mstage=="m0; 1:m1): 0"] <- "M0"
pd$Mstage[pd$Mstage=="m0; 1:m1): 1"] <- "M1"
table(pd$Mstage)

str(pd)

table(pd$OS)
pd$OS[pd$OS=='-'] <- NA
pd$OS <- as.numeric(as.character(pd$OS))
table(pd$OS)

table(pd$OS_Time)
pd$OS_Time[pd$OS_Time=='-'] <- NA
pd$OS_Time <- as.numeric(as.character(pd$OS_Time))
pd$OS_Time <- pd$OS_Time / 12
pd$OS_Time <- round(pd$OS_Time, 2)


pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]
pd <- subset(pd, pd$OS_Time > 0)

str(pd)

table(pd$Age)
table(pd$Tstage)
table(pd$Nstage)
table(pd$Mstage)
table(pd$Gender)
table(pd$OS)
mean(pd$OS_Time)


pd$Sample <- gsub('-', '', pd$Sample)
colnames(exprSet) <- gsub('-', '', colnames(exprSet))

# ===================================================

# ===================================================

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

