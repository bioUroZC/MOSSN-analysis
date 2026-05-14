#=======================================================

#=======================================================

rm(list=ls())
library(dplyr)
library(tibble)
library(tidyr)
library(GEOquery)
library(gdata)
library(data.table)

options(stringsAsFactors = FALSE)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 12)

setwd('/proj/c.zihao/work1/1NT/2CPTAC/KIRC')

#=======================================================

#=======================================================

gsename = "GSE40435"
gse <- getGEO(gsename, destdir = ".")

gpl <- getGEO('GPL10558', destdir = ".")
colnames(Table(gpl))
Table(gpl)[1:10, 1:6]
gpl <- gpl@dataTable@table
colnames(gpl)

gpl <- gpl %>% dplyr::select(ID, "ILMN_Gene")
gpl$gene <- gdata::trim(gpl$ILMN_Gene)
gpl$ILMN_Gene <- NULL
gpl$gene <- as.character(gpl$gene)
gpl$gene[gpl$gene == "" | is.na(gpl$gene)] <- NA
gpl$gene <- sub(" ///.*$", "", gpl$gene)
gpl$gene <- sub(" //.*$", "", gpl$gene)
gpl <- gpl[!is.na(gpl$gene), ]

head(gpl)

#=======================================================

#=======================================================

exprSet <- as.data.frame(exprs(gse[[1]]))
str(exprSet)

exprSet$ID <- rownames(exprSet)
express <- merge(x = gpl, y = exprSet, by = "ID")
express$ID <- NULL

express[which(is.na(express), arr.ind = TRUE)] <- 0

exprSet <- aggregate(x = express[, 2:ncol(express)],
                     by = list(express$gene),
                     FUN = max)

exprSet <- as.data.frame(exprSet)
names(exprSet)[1] <- 'ID'

exprSet$ID <- as.character(exprSet$ID)
exprSet$ID <- sub(" ///.*$", "", exprSet$ID)
exprSet$ID <- sub(" //.*$", "", exprSet$ID)
exprSet <- exprSet[!is.na(exprSet$ID) & exprSet$ID != "", ]
exprSet <- exprSet[!duplicated(exprSet$ID), ]
rownames(exprSet) <- exprSet$ID
exprSet$ID <- NULL
exprSet[1:5, 1:5]

min(exprSet)
max(exprSet)

#=======================================================

#=======================================================

pd <- pData(gse[[1]])
names(pd)
head(pd)
write.csv(pd, file = "pd_info.csv", row.names = TRUE)

tissue  <- pd$`characteristics_ch1.4`
patient <- sub("patient (\\d+),.*", "\\1", pd$title)

tumor_idx  <- which(tissue == "tissue type: clear cell renal cell carcinoma")
normal_idx <- which(tissue == "tissue type: adjacent non-tumour renal tissue")

paired_patients <- intersect(patient[tumor_idx], patient[normal_idx])
cat("Tumor samples:", length(tumor_idx),
    "\nNormal samples:", length(normal_idx),
    "\nPaired patients:", length(paired_patients), "\n")

tumor_keep  <- tumor_idx[patient[tumor_idx]  %in% paired_patients]
normal_keep <- normal_idx[patient[normal_idx] %in% paired_patients]

tumor_samples  <- rownames(pd[tumor_keep, ])
normal_samples <- rownames(pd[normal_keep, ])

tumor_patient_order  <- patient[tumor_keep]
normal_patient_order <- patient[normal_keep]
tumor_samples  <- tumor_samples[order(tumor_patient_order)]
normal_samples <- normal_samples[order(normal_patient_order)]
paired_patients_ordered <- sort(paired_patients)

tumor_expr  <- exprSet[, tumor_samples,  drop = FALSE]
normal_expr <- exprSet[, normal_samples, drop = FALSE]

colnames(tumor_expr)  <- paste0("patient", paired_patients_ordered, "_Tumor")
colnames(normal_expr) <- paste0("patient", paired_patients_ordered, "_Normal")

mergedData <- cbind(tumor_expr, normal_expr)
colnames(mergedData) <- gsub("\\.", "_", colnames(mergedData))
mergedData[is.na(mergedData)] <- 0

link_df <- fread("/proj/c.zihao/work1/1NT/1data/string/links.csv")
link_genes <- unique(c(as.character(link_df$protein1), as.character(link_df$protein2)))
keep_genes <- intersect(rownames(mergedData), link_genes)
mergedData <- mergedData[keep_genes, , drop = FALSE]

cat("Genes retained after STRING link filtering:", nrow(mergedData), "\n")

dim(mergedData)
print(mergedData[1:5, 1:4])

write.csv(mergedData, file = "KIRC_exprSet_filtered.csv", row.names = TRUE)
print(dim(mergedData))

print(min(mergedData))
print(max(mergedData))
