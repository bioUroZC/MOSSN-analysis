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
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(ArrayExpress)
library(RCurl)
library(affy)
library(hgu133plus2.db)
library(AnnotationDbi)
setwd("/proj/c.zihao/work1/1survival/LGG/MTAB3892/data/")

#=======================================================

#=======================================================

ftp_url <- "ftp://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/892/E-MTAB-3892/Files/"
filenames <- getURL(ftp_url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
files <- strsplit(filenames, "\r*\n")[[1]]
files <- files[nchar(files) > 0] 

if(!dir.exists("MTAB3892_all_files")) dir.create("MTAB3892_all_files")

for (f in files) {
  destfile <- file.path("MTAB3892_all_files", f)
  cat("Downloading", f, "\n")
  tryCatch({
    download.file(paste0(ftp_url, f), destfile = destfile, mode = "wb")
  }, error = function(e) {
    cat("Failed to download:", f, "\n")
  })
}

#=======================================================

#=======================================================

setwd("/proj/c.zihao/work1/1survival/LGG/MTAB3892/data/MTAB3892_all_files")
data_raw <- ReadAffy()  

eset <- rma(data_raw)   
exprs_matrix <- exprs(eset)
dim(exprs_matrix)       
head(exprs_matrix)[, 1:5]  


pheno <- read.table("E-MTAB-3892.sdrf.txt", sep = "\t",
                    header = TRUE, stringsAsFactors = FALSE)
head(pheno)



probes <- rownames(exprs_matrix)
gene_symbols <- mapIds(hgu133plus2.db, 
                       keys = probes,
                       column = "SYMBOL",
                       keytype = "PROBEID",
                       multiVals = "first")

exprs_annot <- data.frame(GeneSymbol = gene_symbols, exprs_matrix, check.names = FALSE)
head(exprs_annot)[1:5,1:5]
colnames(exprs_annot)


exprSet <- aggregate(x = exprs_annot[,2:ncol(exprs_annot)],
                     by = list(exprs_annot$GeneSymbol),
                     FUN = max)

exprSet[1:5,1:5]

exprSet <- as.data.frame(exprSet)
names(exprSet)[1] <- 'ID'
rownames(exprSet) <- exprSet$ID
exprSet$ID <- NULL
exprSet[1:5,1:5]

min(exprSet)
max(exprSet)


#=======================================================

#=======================================================

pheno <- read.table("E-MTAB-3892.sdrf.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE, quote="")
head(pheno)
dim(pheno)


pheno <- subset(pheno, select=c("Array.Data.File",
                                "Characteristics.sex.",
                                "Characteristics.age.",
                                "Factor.Value.histology.",
                                "Factor.Value.histology_grade.",
                                "Characteristics.os.event.",
                                "Characteristics.os.delay."))
colnames(pheno) <- c("Sample", "Gender", "Age", "Histology", "Grade", "OS", "OS_Time")

head(pheno)
str(pheno)

table(pheno$Histology)
table(pheno$Grade)


LGG_histology <- c(
  "Low-grade Oligodendroglioma",
  "Low-grade Oligoastrocytoma",
  "High-grade Oligodendroglioma",
  "High-grade Oligoastrocytoma",
  "Diffuse astrocytoma"
)


pd <- pheno[pheno$Grade %in% c(2,3) & pheno$Histology %in% LGG_histology, ]


table(pd$Gender)
pd$Gender[pd$Gender=="female"] <- "Female"
pd$Gender[pd$Gender=="male"] <- "Male"
table(pd$Gender)

table(pd$Age)
table(pd$Grade)

table(pd$OS)

table(pd$OS_Time)

pd$OS_Time <- pd$OS_Time / 12
pd$OS_Time <- round(pd$OS_Time, 2)
str(pd)

pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]
pd <- subset(pd, pd$OS_Time>0)
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
getwd()
setwd("/proj/c.zihao/work1/1survival/LGG/MTAB3892/data/")

write.csv(exprSet, file = "exprSet.csv")
write.csv(pd, file = "pd.csv")

