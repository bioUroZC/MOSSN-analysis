#=======================================================

#=======================================================

rm(list=ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


library(GEOquery)
library(Biobase)
library(limma)
library(dplyr)
library(tidyr)
library(gdata)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)

setwd(paste0(PROJ_ROOT, "/1survival/CRC/GSE41258/data/"))

#=======================================================

#=======================================================

gse <- getGEO("GSE41258", destdir = ".")
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
  tidyr::separate("Gene.Symbol", c('gene', 'symbol'), sep = '\\///') %>%
  dplyr::select("ID", 'gene')

gpl$gene <- gdata::trim(gpl$gene)

#=======================================================

#=======================================================

exprSet <- as.data.frame(exprs(gse$GSE41258_series_matrix.txt.gz))
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

exprSet <- log2(exprSet+1)

min(exprSet)
max(exprSet)

#=======================================================

#=======================================================


pd <- pData(gse$GSE41258_series_matrix.txt.gz)

head(pd)
names(pd)

table(pd$`tissue:ch1`)


pd <- subset(pd, select=c("geo_accession","age:ch1" ,  "gender:ch1"  ,
                          "anatomic location:ch1"  , "group stage:ch1" ,
                          "t:ch1" , "n:ch1" , "m:ch1"  ,
                          "fup status (ned,  no evidence of disease, awd, alive with disease, aun, alive unknown, dod, dead of disease, doc, dead of other cause, dun, dead, cause unknown):ch1",
                          "fup interval:ch1" ))

head(pd)

names(pd) <- c("Sample",  "Age", "Gender", "Location", "Stage", 
               "Tstage", "Nstage", "Mstage",  "OS", "OS_Time")
head(pd)

str(pd)

table(pd$Location)

pd <- pd[!is.na(pd$Location), ]

str(pd)

table(pd$Gender)

pd$Gender[pd$Gender=="F"] <- "Female"
pd$Gender[pd$Gender=="M"] <- "Male"

table(pd$OS)
pd$OS <- ifelse(pd$OS %in% c("DOD", "DOC", "DUN"), 1, 0)
pd$OS  <-  as.numeric(as.character(pd$OS))
table(pd$OS)


pd$OS_Time <-  as.numeric(as.character(pd$OS_Time))
pd$OS_Time <- pd$OS_Time / 12
pd$OS_Time <- round(pd$OS_Time, 2)


pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]
pd <- subset(pd, pd$OS_Time > 0)
str(pd)
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

