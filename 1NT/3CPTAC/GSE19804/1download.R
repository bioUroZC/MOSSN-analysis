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

setwd(paste0(PROJ_ROOT, "/1NT/3CPTAC/GSE19804"))

#=======================================================

#=======================================================

gsename <- "GSE19804"
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

head(gpl)

#=======================================================

#=======================================================

exprSet <- as.data.frame(exprs(gse$GSE19804_series_matrix.txt.gz))
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

cat("Min expression value:", min(exprSet), "\n")
cat("Max expression value:", max(exprSet), "\n")

#=======================================================

#=======================================================

pd_raw <- pData(gse$GSE19804_series_matrix.txt.gz)

# Title format: "Lung Cancer 2T" / "Lung Normal 2N"
# source_name_ch1: "frozen tissue of primary tumor" / "frozen tissue of adjacent normal"
pd <- data.frame(
  Sample      = rownames(pd_raw),
  patient_id  = sub("Lung (?:Cancer|Normal) ([0-9]+)[TN]$", "\\1", pd_raw$title, perl = TRUE),
  tissue_type = ifelse(grepl("primary tumor", pd_raw$source_name_ch1), "Cancer", "Normal"),
  stringsAsFactors = FALSE
)

pd$patient_id <- paste0('Sam', pd$patient_id)
pd$sample_id <- paste(pd$patient_id, pd$tissue_type, sep = "_")

head(pd)
cat("Total samples  :", nrow(pd), "\n")
cat("Cancer samples :", sum(pd$tissue_type == "Cancer"), "\n")
cat("Normal samples :", sum(pd$tissue_type == "Normal"), "\n")


# Check 1:1 pairing
pair_check <- table(pd$patient_id, pd$tissue_type)
unpaired <- rownames(pair_check)[!(pair_check[, "Cancer"] == 1 & pair_check[, "Normal"] == 1)]
if (length(unpaired) == 0) {
  cat("Pairing check: all", nrow(pair_check), "patients have exactly 1 Normal + 1 Cancer sample.\n")
} else {
  cat("Unpaired patients:", length(unpaired), "\n")
  print(pair_check[unpaired, ])
}



#=======================================================

#=======================================================

samplesname <- intersect(pd$Sample, colnames(exprSet))
samplesname <- unique(samplesname)
pd <- pd[which(pd$Sample %in% samplesname),]
exprSet <- exprSet[,which(colnames(exprSet) %in% samplesname)]
colnames(exprSet)

#=======================================================

#=======================================================

colnames(exprSet) <- pd$sample_id[match(colnames(exprSet), pd$Sample)]

exprSet[1:5,1:5]
exprSet <- round(exprSet, 5)
write.csv(exprSet, file = "exprSet.csv")
write.csv(pd, file = "pd.csv")


cat("Min expression value:", min(exprSet), "\n")
cat("Max expression value:", max(exprSet), "\n")



link_df <- fread(paste0(PROJ_ROOT, "/1NT/1data/string/links.csv"))
link_genes <- unique(c(as.character(link_df$protein1), as.character(link_df$protein2)))
keep_genes <- intersect(rownames(mergedData), link_genes)
mergedData <- mergedData[keep_genes, , drop = FALSE]

cat("Genes retained after STRING link filtering:", nrow(mergedData), "\n")

print(min(mergedData))
print(max(mergedData))

mergedData <- log2(mergedData + 1)


print(min(mergedData))
print(max(mergedData))

write.csv(mergedData, file = "eLUAD_exprSet_filtered.csv", row.names = TRUE)
print(dim(mergedData))

print(mergedData[1:5,1:5])
print(colnames(mergedData))
