#===================================================

#===================================================

rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")

library(dplyr)
library(TCGAbiolinks)
library(data.table)
library(SummarizedExperiment)
library(openxlsx)
library(survival)
library(survminer)
library(dplyr)

#===================================================

#===================================================

# Create a vector of TCGA project IDs
x <- c(
  "TCGA-LUAD", 
  "TCGA-LUSC", 
  "TCGA-COAD",  
  "TCGA-READ",  
  "TCGA-LIHC", 
  
  "TCGA-BRCA", 
  "TCGA-STAD", 
  "TCGA-ESCA", 
  "TCGA-PRAD",
  
  "TCGA-BLCA", 
  "TCGA-HNSC", 
  "TCGA-KIRC", 
  "TCGA-CHOL"
)



# Set working directory to save output files
setwd(paste0(PROJ_ROOT, '/1NT/1data/TCGA/'))


# Loop through each cancer project
for (pr in 1:length(x)) {
  
  project_id  <- x[pr]
  cat("Processing project:", project_id, "\n")
  
  query <- GDCquery(
    project = project_id,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  
  GDCdownload(query, method = "api", files.per.chunk = 10)
  
  #=======================================================
  
  expr <- GDCprepare(query = query)
  TPM <- as.data.frame(assay(expr, i = "tpm_unstrand"))
  
  #Gene annotation
  anno <- as.data.frame(expr@rowRanges@elementMetadata)
  names(anno)
  anno <- subset(anno, select = c('gene_name', 'gene_id',  "gene_type"))
  table(anno$gene_type)
  anno <- subset(anno, anno$gene_type == "protein_coding")
  anno$gene_type <- NULL
  
  TPM$gene_id <- rownames(TPM)
  TPM <- merge(anno, TPM, by = 'gene_id')
  TPM$gene_id <- NULL
  
  #Missing values
  TPM[which(is.na(TPM), arr.ind = TRUE)] <- 0
  TPM[1:5,1:5]
  
  #=======================================================
  
  #genes with multiple probes
  exprSet <- aggregate(x = TPM[, 2:ncol(TPM)], 
                       by = list(TPM$gene_name), 
                       FUN = max)
  exprSet[1:5,1:5]
  
  names(exprSet)[1] <- "ID"
  rownames(exprSet) <- exprSet$ID
  exprSet$ID <- NULL
  exprSet[1:5,1:5]
  max(exprSet)
  min(exprSet)
  
  dim(exprSet)
  
  exprSet <- na.omit(exprSet)
  colnames(exprSet) <- chartr(old = '-', new = "_", colnames(exprSet))
  exprSet[1:4,1:4]
  
  dim(exprSet)
  
  #=======================================================
  
  #remove normal samples
  metad <- data.frame(names(exprSet))
  colnames(metad)[1] <- 'sample'
  metad$id <- substr(metad$sample, start = 1, stop = 12)
  metad$tape <- substr(metad$sample, start = 14, stop = 16)
  metad <- metad[order(metad$id), ]
  table(metad$tape)
  print(table(metad$tape))
  
  mNormal <- subset(metad, metad$tape == "11A")
  mNormal <- mNormal %>%
    dplyr::distinct(id, .keep_all = T)
  
  idNormal <- mNormal$id
  
  mTumor <- subset(metad, metad$tape == "01A")
  mTumor <- mTumor %>%
    dplyr::distinct(id, .keep_all = T)
  
  idTumor <- mTumor$id
  
  idPairs <- intersect(idNormal, idTumor)
  
  cat("length of pairs: ", length(idPairs), "\n")
  
  mNormal <- mNormal[which(mNormal$id %in% idPairs),]
  mTumor <-  mTumor[which(mTumor$id %in% idPairs),]
  
  
  metad <- rbind(mNormal, mTumor)
  exprSet <- exprSet[, which(colnames(exprSet) %in% metad$sample)]
  colnames(exprSet) <- substr(colnames(exprSet), start = 1, stop = 16)
  exprSet[1:5,1:5]
  
  
  # Save the processed data to a CSV file
  output_filename <- paste0(project_id, ".csv")
  write.csv(exprSet, file = output_filename)
  
  cat("Data for", project_id, "saved to", output_filename, "\n")
}

cat("All projects processed successfully!\n")