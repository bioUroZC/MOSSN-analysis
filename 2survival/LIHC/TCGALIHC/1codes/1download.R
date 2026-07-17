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
  "TCGA-LIHC"
)


exprdir <- paste0(PROJ_ROOT, '/1survival/LIHC/TCGALIHC/data/')
surdir <- paste0(PROJ_ROOT, "/1survival/")
cancername <- "LIHC"

# Set working directory to save output files
setwd(exprdir)


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
  metad <- subset(metad, metad$tape == "01A")
  metad <- metad %>% distinct(id, .keep_all = TRUE)
  exprSet <- exprSet[, which(colnames(exprSet) %in% metad$sample)]
  colnames(exprSet) <- substr(colnames(exprSet), start = 1, stop = 12)
  exprSet[1:5,1:5]
  
}

exprSet[1:5,1:5]


gene_iqr <- apply(exprSet, 1, IQR)
mean_iqr <- mean(gene_iqr)
median_iqr <- median(gene_iqr)
mean_iqr
median_iqr

min(exprSet)
max(exprSet)


# ===================================================

# ===================================================

setwd(surdir)
survivaldata <- read.xlsx("mmc1.xlsx", rowNames = TRUE)
survivaldata <- subset(survivaldata, select = c("bcr_patient_barcode", "type",
                                                "age_at_initial_pathologic_diagnosis",
                                                "gender",
                                                "ajcc_pathologic_tumor_stage",
                                                "clinical_stage",
                                                "OS", "OS.time"))

survivaldata <- subset(survivaldata, select = c("bcr_patient_barcode", "type",
                                                "ajcc_pathologic_tumor_stage",
                                                "OS", "OS.time"))
names(survivaldata) <- c("Sample", "cancertype", "Stage", "OS", "OS_Time")

survivaldata <- subset(survivaldata, survivaldata$cancertype == cancername)
head(survivaldata)
table(survivaldata$cancertype)
survivaldata$cancertype <- NULL
dim(survivaldata)
survivaldata <- survivaldata %>% distinct(Sample, .keep_all = TRUE)
dim(survivaldata)
survivaldata$Sample <- chartr(old = '-', new = '_', x=survivaldata$Sample)
str(survivaldata)

table(survivaldata$Stage)
survivaldata$Stage <- ifelse(grepl("^Stage I($|[^A-Z])", survivaldata$Stage), "I",
                             ifelse(grepl("^Stage II($|[^I])", survivaldata$Stage), "II",
                                    ifelse(grepl("^Stage III", survivaldata$Stage), "III",
                                           ifelse(grepl("^Stage IV", survivaldata$Stage), "IV", NA))))
table(survivaldata$Stage)

survivaldata <- as.data.frame(survivaldata)
survivaldata$OS_Time <- survivaldata$OS_Time / 365
pd <- survivaldata

pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]
pd$OS_Time <- round(pd$OS_Time, 2)
pd <- subset(pd, pd$OS_Time>0)

# ===================================================

# ===================================================


samplesname <- intersect(pd$Sample, colnames(exprSet))
samplesname <- unique(samplesname)
pd <- pd[which(pd$Sample %in% samplesname),]
exprSet <- exprSet[,which(colnames(exprSet) %in% samplesname)]
colnames(exprSet)

# ===================================================

# ===================================================

setwd(exprdir)
getwd()

exprSet[1:5,1:5]
write.csv(exprSet, file = "exprSet.csv")
write.csv(pd, file = "pd.csv")
cat("Data for", project_id, "saved to", "exprSet.csv", "\n")
