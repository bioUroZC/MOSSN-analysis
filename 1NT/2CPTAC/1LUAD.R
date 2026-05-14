rm(list = ls())
library(data.table)

tumor_url  <- "https://linkedomics.org/data_download/CPTAC-LUAD/HS_CPTAC_LUAD_rnaseq_uq_rpkm_log2_NArm_TUMOR.cct"
normal_url <- "https://linkedomics.org/data_download/CPTAC-LUAD/HS_CPTAC_LUAD_rnaseq_uq_rpkm_log2_NArm_NORMAL.cct"
output_dir <- "/proj/c.zihao/work1/1NT/2CPTAC/LUAD"

tumorData <- fread(tumor_url, data.table = FALSE, check.names = FALSE)
normalData <- fread(normal_url, data.table = FALSE, check.names = FALSE)

print(tumorData[1:5, 1:5])
print(normalData[1:5, 1:5])

tumor_gene_col <- intersect(c("GeneSymbol", "gene", "NAME", "Name", "id", "ID"), colnames(tumorData))[1]
normal_gene_col <- intersect(c("GeneSymbol", "gene", "NAME", "Name", "id", "ID"), colnames(normalData))[1]

if (is.na(tumor_gene_col) || is.na(normal_gene_col)) {
    stop("Cannot find gene identifier column in tumor or normal matrix.")
}

tumor_sample_cols <- colnames(tumorData)[sapply(tumorData, is.numeric)]
normal_sample_cols <- colnames(normalData)[sapply(normalData, is.numeric)]

tumor_sample_cols <- setdiff(tumor_sample_cols, tumor_gene_col)
normal_sample_cols <- setdiff(normal_sample_cols, normal_gene_col)

if (length(tumor_sample_cols) == 0 || length(normal_sample_cols) == 0) {
    stop("Cannot find numeric expression sample columns in tumor or normal matrix.")
}

tumorData <- tumorData[, c(tumor_gene_col, tumor_sample_cols), drop = FALSE]
normalData <- normalData[, c(normal_gene_col, normal_sample_cols), drop = FALSE]

colnames(tumorData)[1] <- "GeneSymbol"
colnames(normalData)[1] <- "GeneSymbol"

tumorData$GeneSymbol <- trimws(as.character(tumorData$GeneSymbol))
normalData$GeneSymbol <- trimws(as.character(normalData$GeneSymbol))

tumorData <- tumorData[!is.na(tumorData$GeneSymbol) & tumorData$GeneSymbol != "", , drop = FALSE]
normalData <- normalData[!is.na(normalData$GeneSymbol) & normalData$GeneSymbol != "", , drop = FALSE]

tumorData <- tumorData[!duplicated(tumorData$GeneSymbol), , drop = FALSE]
normalData <- normalData[!duplicated(normalData$GeneSymbol), , drop = FALSE]

colnames(tumorData)[-1] <- paste0(colnames(tumorData)[-1], "_Tumor")
colnames(normalData)[-1] <- paste0(colnames(normalData)[-1], "_Normal")

tumor_ids <- sub("_Tumor$", "", colnames(tumorData)[-1])
normal_ids <- sub("_Normal$", "", colnames(normalData)[-1])
paired_ids <- intersect(tumor_ids, normal_ids)

cat(
    "Tumor samples:", length(tumor_ids),
    "\nNormal samples:", length(normal_ids),
    "\nPaired samples:", length(paired_ids), "\n"
)

tumor_keep <- c("GeneSymbol", paste0(paired_ids, "_Tumor"))
normal_keep <- c("GeneSymbol", paste0(paired_ids, "_Normal"))

tumor_paired <- tumorData[, tumor_keep, drop = FALSE]
normal_paired <- normalData[, normal_keep, drop = FALSE]

mergedData <- merge(tumor_paired, normal_paired, by = "GeneSymbol", sort = FALSE)
rownames(mergedData) <- mergedData$GeneSymbol
mergedData$GeneSymbol <- NULL

colnames(mergedData) <- gsub("\\.", "_", colnames(mergedData))
mergedData[] <- lapply(mergedData, as.numeric)
mergedData[is.na(mergedData)] <- 0

link_df <- fread("/proj/c.zihao/work1/1NT/1data/string/links.csv")
link_genes <- unique(c(as.character(link_df$protein1), as.character(link_df$protein2)))
keep_genes <- intersect(rownames(mergedData), link_genes)
mergedData <- mergedData[keep_genes, , drop = FALSE]

cat("Genes retained after STRING link filtering:", nrow(mergedData), "\n")

min_value <- min(as.matrix(mergedData), na.rm = TRUE)
max_value <- max(as.matrix(mergedData), na.rm = TRUE)
cat("Before shift - min:", min_value, "max:", max_value, "\n")

if (min_value < 0) {
    mergedData <- mergedData - min_value
}

final_min <- min(as.matrix(mergedData), na.rm = TRUE)
final_max <- max(as.matrix(mergedData), na.rm = TRUE)
cat("Final matrix - min:", final_min, "max:", final_max, "\n")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
    mergedData,
    file = file.path(output_dir, "LUAD_exprSet_filtered.csv"),
    row.names = TRUE
)

print(dim(mergedData))
print(mergedData[1:5, 1:5])
