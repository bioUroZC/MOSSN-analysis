# ====================================================================================

# ====================================================================================

rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


library(dplyr)
library(igraph)
library(e1071)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)

# ====================================================================================

# ====================================================================================

setwd(paste0(PROJ_ROOT, "/1survival/KIRC/PCWAG/data/"))

links <- 'https://pcawg-hub.s3.us-east-1.amazonaws.com/download/tophat_star_fpkm_uq.v2_aliquot_gl.sp.log'

exprsess <- fread(links)
exprsess <- as.data.frame(exprsess)
exprsess[1:5,1:5]


names(exprsess)[1] <- 'gene_id'

exprsess$gene_id <- sub("\\..*", "", exprsess$gene_id)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_ids <- exprsess$gene_id
genes_info <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'description'),
                    filters = 'ensembl_gene_id', 
                    values = gene_ids, 
                    mart = mart)
head(genes_info)

genes_info <- subset(genes_info, select=c('ensembl_gene_id', 'external_gene_name'))
names(genes_info) <- c('gene_id', 'gene_name')
data <- merge(genes_info, exprsess, by="gene_id")
data$gene_id <- NULL
data[1:4,1:6]

# ====================================================================================

# ====================================================================================

pd <- fread("survival_sp")
pd <- as.data.frame(pd)
head(pd)

project  <- fread("project_code_sp")

table(project$dcc_project_code)

grep("KIRC", project$dcc_project_code, value = TRUE)

project <- subset(project, project$dcc_project_code == "KIRC-US")

pd <- pd[which(pd$xena_sample %in% project$icgc_specimen_id),]

pd <- pd %>%
  dplyr::select('xena_sample', '_EVENT', '_TIME_TO_EVENT')

head(pd)

names(pd) <- c("Sample", "OS", "OS_Time")

str(pd)

table(pd$OS)
pd$OS[pd$OS=="alive"] <- 0
pd$OS[pd$OS=="deceased"] <- 1
table(pd$OS)

pd$OS_Time <- pd$OS_Time / 365
pd$OS_Time <- round(pd$OS_Time, 2)


pd <- pd[!is.na(pd$OS) & !is.na(pd$OS_Time), ]
pd <- subset(pd, pd$OS_Time > 0)



#===================================================

#===================================================


keep_samples <- pd$Sample
data_filtered <- data[, c("gene_name", intersect(colnames(data), keep_samples))]
data_filtered[1:5,1:5]

exprSet <- aggregate(x = data_filtered[, 2:ncol(data_filtered)],
                     by = list(data_filtered$gene_name),
                     FUN = max)

exprSet[1:5,1:5]
exprSet <- as.data.frame(exprSet)
exprSet <- exprSet[-1, ]
names(exprSet)[1] <- 'ID'
rownames(exprSet) <- exprSet$ID
exprSet$ID <- NULL
exprSet[1:5,1:5]

min(exprSet)
max(exprSet)

exprSet <- exprSet + 10

min(exprSet)
max(exprSet)


#===================================================

#===================================================


samplesname <- intersect(pd$Sample, colnames(exprSet))
samplesname <- unique(samplesname)
pd <- pd[which(pd$Sample %in% samplesname),]
exprSet <- exprSet[,which(colnames(exprSet) %in% samplesname)]
colnames(exprSet)

pd$Sample <- gsub('-', '', pd$Sample)
colnames(exprSet) <- gsub('-', '', colnames(exprSet))

pd$Sample
colnames(exprSet)

#=======================================================

#=======================================================

exprSet[1:5,1:5]
write.csv(exprSet, file = "exprSet.csv")
write.csv(pd, file = "pd.csv")

