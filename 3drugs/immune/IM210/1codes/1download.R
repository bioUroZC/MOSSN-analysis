# ===================================================

# ===================================================


rm(list = ls())
library(dplyr)
library(tidyr)

# ===================================================

# ===================================================

setwd("/proj/c.zihao/work1/3drugs/immune/IM210/data/")

exprSet <- read.csv("IMvigor210_exprSet.csv", header = T, row.names = 1)

# ===================================================

# ===================================================

pheno <- read.csv('IMvigor210_FollowUp.csv', header = T, row.names = 1)

colnames(pheno)

pd <- subset(pheno, select=c("Sex","Tissue" , "Met.Disease.Status" ,
                             "Best.Confirmed.Overall.Response", "binaryResponse",
                             "os", "censOS"))

names(pd) <- c("Gender", "Tissue", "Met", "Response", "ResponseGroup",  "OS_Time", "OS")
str(pd)
pd$OS_Time <- round(pd$OS_Time, 2)
table(pd$Tissue)
pd$Sample <- rownames(pd)

pd <- pd %>%
  dplyr::select(Sample, Response)

head(pd)



table(pd$Response)
pd$Response[pd$Response=="NE"] <- NA
pd$Response[pd$Response=="CR"] <- 4
pd$Response[pd$Response=="PR"] <- 3
pd$Response[pd$Response=="PD"] <- 1
pd$Response[pd$Response=="SD"] <- 2

pd <- na.omit(pd)
pd$Response <- as.numeric(as.character(pd$Response))

pd$Sample[1:10]
colnames(exprSet)[1:10]

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


min(exprSet)
max(exprSet)
