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
set.seed(1234)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)

setwd(paste0(PROJ_ROOT, "/1survival/BRCA/GSE20685/data"))


#=======================================================

#=======================================================

gsename <- "GSE20685"
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

exprSet <- as.data.frame(exprs(gse$GSE20685_series_matrix.txt.gz))
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

pd <- pData(gse$GSE20685_series_matrix.txt.gz)
head(pd)
names(pd)

pd <- subset(pd, select=c("geo_accession", 
                          "age at diagnosis:ch1",
                          "event_death:ch1",
                          "event_metastasis:ch1",
                          "regional_relapse:ch1",
                          
                          "follow_up_duration (years):ch1",
                          "time_to_metastasis (years):ch1",
                          "time_to_relapse (years):ch1"))

colnames(pd) <- c(
  "geo_accession", "age", "event_death", "event_metastasis",
  "regional_relapse", "follow_up", "time_to_metastasis", "time_to_relapse"
)

head(pd)
str(pd)

num_cols <- c("age", "event_death", "event_metastasis", 
              "regional_relapse", "follow_up", 
              "time_to_metastasis", "time_to_relapse")

pd[num_cols] <- lapply(pd[num_cols], as.numeric)


pd$EFS <- with(pd, ifelse(event_death == 1 | event_metastasis == 1 | regional_relapse == 1, 1, 0))

pd$EFS_Time <- apply(
  pd[, c("time_to_metastasis", "time_to_relapse", "follow_up")], 
  1, 
  function(x) min(x, na.rm = TRUE)
)


str(pd)

pd <- subset(pd, select=c("geo_accession", "age", "EFS", "EFS_Time" ))
names(pd) <- c("Sample", "Age", "EFS", "EFS_Time" )

pd$EFS_Time <- round(pd$EFS_Time, 2)
str(pd)

pd <- pd[!is.na(pd$EFS) & !is.na(pd$EFS_Time), ]
pd <- subset(pd, pd$EFS_Time>0)
head(pd)

str(pd)
table(pd$EFS)
mean(pd$EFS_Time)

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
