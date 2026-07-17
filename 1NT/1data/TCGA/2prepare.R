# ===================================================

# ===================================================

rm(list = ls())

PROJ_ROOT <- Sys.getenv("MOSSN_ROOT", "/proj/c.zihao/work1")


library(dplyr)
library(TCGAbiolinks)
library(data.table)
library(SummarizedExperiment)
library(rvest)
library(tidyr)

# =======================================================

# =======================================================

COAD <- read.csv(paste0(PROJ_ROOT, "/1NT/1data/TCGA/TCGA-COAD.csv"), header = T, row.names = 1)
READ <- read.csv(paste0(PROJ_ROOT, "/1NT/1data/TCGA/TCGA-READ.csv"), header = T, row.names = 1)

COAD$genes <- rownames(COAD)
READ$genes <- rownames(READ)

CRC <- merge(COAD, READ, by='genes')
rownames(CRC) <- CRC$genes
CRC$genes <- NULL

setwd(paste0(PROJ_ROOT, '/1NT/1data/TCGA/'))
write.csv(CRC, file = 'TCGA-CRC.csv')
